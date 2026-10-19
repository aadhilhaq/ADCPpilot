"""The web front end: a local page that sets a run up and watches it.

By Aadhil Haq

Standard library only, and bound to localhost by default. A tool people run
on a login node should not open a port to the rest of the cluster, and the
token in the URL means a second person on the same node cannot drive a user
else's session by guessing the port.

Nothing heavy happens in here. Every request either reads a file, validates
something, or hands work to SLURM.
"""
import csv
import glob
import json
import os
import secrets
import re
import shutil
import stat
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from . import agfr_target, fetch, job as jobmod, page, peptide, pipeline, report
from . import sizing
from . import slurm, structure
from . import tools as toolmod
from . import adcp11 as adcp11_mod
from . import mac10 as mac10_mod

MAX_BODY = 4 << 20
# Uploads have their own ceiling. Every other request is a small JSON object,
# for which 4 MB is generous; a receptor with waters and several models can
# run to tens of megabytes.
MAX_UPLOAD = 64 << 20

# A PDBQT is a receptor already prepared for AutoDock. This workflow prepares
# its receptor itself, from a PDB or mmCIF, and has no step that builds maps
# from a PDBQT as it stands, so one is refused with the way that does work.
PDBQT_REFUSED = ("%s is a PDBQT receptor, which this tool cannot use as it "
                 "stands: it prepares the receptor itself, from a PDB or "
                 "mmCIF file. To dock into a receptor you prepared yourself, "
                 "make a target file from it in AGFR or AGFRGUI and give that "
                 ".trg in the target file box instead.")


# Pending reasons that will not clear on their own. A job held for one of
# these is broken rather than queued, and showing it the same amber dot as a
# waiting job is misleading.
STUCK_REASONS = {
    "PartitionTimeLimit",      # asks for longer than the partition allows
    "PartitionNodeLimit",
    "PartitionConfig",
    "QOSMaxWallDurationPerJobLimit",
    "QOSMaxNodePerJobLimit",
    "InvalidAccount",
    "InvalidQOS",
    "AccountNotAllowed",
    "PartitionDown",
    "PartitionInactive",
    "ReqNodeNotAvail",
    "BadConstraints",
}


# Where a problem report is sent: the addresses the credit line names.
REPORT_TO = "aadhil.haq@tamu.edu"
REPORT_CC = "sfernando@tamu.edu"


def _prep_detail(value):
    """The detail of a preparation log step, after "done" or "skipped" and
    its separator ("; " in current logs, " - " in older ones)."""
    import re
    return re.sub(r"^(done|skipped)\s*(;|-)\s*", "", value)[:160]


def prep_warnings(rundir, j):
    """What the receptor preparation did not do that a result depends on.

    Preparation degrades rather than stops, a pdb2pqr that cannot run is
    logged and the docking proceeds with hydrogens added geometrically, and
    the log is not read before a results table. A run docked without
    protonation at its pH is still a result, but it must not look like a
    clean one.
    """
    log = os.path.join(rundir, "receptor", "preparation.log")
    if not os.path.isfile(log):
        return []
    steps = {}
    for line in open(log, encoding="utf-8", errors="replace"):
        s = line.strip()
        for name in ("protonation", "protonation retry", "minimisation",
                     "receptor pdbqt", "reduce"):
            if s.startswith(name + ":"):
                steps[name] = s[len(name) + 1:].strip()
    out = []
    red = steps.get("reduce", "")
    if red and not red.startswith("done"):
        out.append("The receptor was not protonated with reduce, as the run "
                   "asked (%s). PDB2PQR with PROPKA at pH %g protonated it "
                   "instead. See receptor/preparation.log."
                   % (_prep_detail(red), float(j.get("ph", 7.4))))
    prot = steps.get("protonation", "")
    retry = steps.get("protonation retry", "")
    if prot and not prot.startswith("done") and not retry.startswith("done"):
        out.append("The receptor was NOT protonated by pdb2pqr (%s). Hydrogens "
                   "were added geometrically, so the run's pH is not reflected "
                   "in the receptor. See receptor/preparation.log."
                   % _prep_detail(prot))
    elif prot and not prot.startswith("done") and retry.startswith("done"):
        out.append("The receptor was protonated without PROPKA (%s), so its "
                   "titratable residues have the force field's standard states "
                   "and the run's pH is not reflected in them. See "
                   "receptor/preparation.log." % _prep_detail(prot))
    pdbqt = steps.get("receptor pdbqt", "")
    if pdbqt and not pdbqt.startswith("done"):
        out.append("The docking receptor was written from the structure as "
                   "read, not from the prepared one: prepare_receptor could "
                   "not use the prepared structure. Rebuilt side chains, caps, "
                   "protonation and minimisation are not in what was docked "
                   "against. See the prepare log.")
    mini = steps.get("minimisation", "")
    if j.get("minimise", True) and mini and not mini.startswith("done"):
        out.append("The restrained minimisation did not run (%s)."
                   % _prep_detail(mini))
    return out


def fill_names(results, j):
    """Give every results row its peptide's name, from the job file.

    Runs docked before the dock stage wrote names have a results.json without
    them, and an old run's report is never rewritten. The job file is the
    authority for names, and is read whenever the page asks for results.

    A row that already has a name keeps it. Matching is by sequence.
    """
    if not results:
        return results
    rows = results if isinstance(results, list) else results.get("rows")
    if not isinstance(rows, list):
        return results
    by_seq = {}
    for n, s in jobmod.peptide_list(j):
        by_seq.setdefault(s, n)
    for r in rows:
        if isinstance(r, dict) and not r.get("name"):
            r["name"] = by_seq.get(r.get("peptide") or r.get("sequence"), "")
    return results


def fill_rings(results, j):
    """Give every results row the ring its peptide was docked with.

    The same idea as fill_names: a report written before rows carried the
    ring is read with it, from the job file, which for such a run holds one
    ring for every peptide. A row that already says keeps what it says.
    """
    if not results:
        return results
    rows = results if isinstance(results, list) else results.get("rows")
    if not isinstance(rows, list):
        return results
    for r in rows:
        if isinstance(r, dict) and not r.get("ring"):
            mode = jobmod.ring_for(j, r.get("peptide") or r.get("sequence")
                                   or "")
            r["ring"] = peptide.RING_NAMES.get(mode, mode)
    return results


# Energy terms read for a run, keyed by the run and the modification times of
# what they are read from, so a page polling an open run reads the files once.
_TERMS_CACHE = {}


def fill_terms(results, rundir):
    """Give an older run's rows the four MM-GBSA terms, read from its files.

    The same idea as fill_names: a run docked before the panel existed has a
    results.json without them, and its report is not rewritten just because
    someone looked at it. The terms are in the MM-GBSA output the run already
    holds, so they are read here, in memory, whenever the page asks.

    Also says how many peptides have no interface columns, and how many of
    those could still be measured from the files on disk, so the page can
    offer to measure them.
    """
    if not results or not rundir:
        return results
    rows = results.get("rows") if isinstance(results, dict) else None
    if not isinstance(rows, list):
        return results
    try:
        stamp = (os.path.getmtime(os.path.join(rundir, "results.json")),
                 os.path.getmtime(os.path.join(rundir, "mmgbsa"))
                 if os.path.isdir(os.path.join(rundir, "mmgbsa")) else 0)
    except OSError:
        stamp = None
    got = _TERMS_CACHE.get(rundir)
    if got is None or got[0] != stamp:
        extra, clusters, missing, can = {}, {}, 0, 0
        for f in sorted(glob.glob(os.path.join(rundir, "results", "*.json"))):
            try:
                r = json.load(open(f))
            except (OSError, ValueError):
                continue
            seq = r.get("peptide")
            if seq and not r.get("cluster_total"):
                # A run docked before the dock stage kept the cluster sizes
                # still has them in ADCP's log.
                cc = pipeline.cluster_counts(rundir, r)
                if cc:
                    clusters[seq] = cc
            if not seq or not r.get("mmgbsa"):
                continue
            panel = r.get("panel") or {}
            # From the panel where one was measured, since its report may not
            # have been rebuilt yet; from the MM-GBSA output where not.
            t = ({k: panel[k] for k in ("pose", "mm_vdw", "mm_eel", "mm_egb",
                                        "mm_esurf") if k in panel}
                 if panel else pipeline.energy_terms(rundir, r))
            if t:
                extra[seq] = t
            if panel.get("bsa") is None:
                missing += 1
                if pipeline.measurable(rundir, r):
                    can += 1
        got = (stamp, extra, clusters, missing, can)
        _TERMS_CACHE[rundir] = got
    _stamp, extra, clusters, missing, can = got
    for r in rows:
        if not isinstance(r, dict):
            continue
        cc = clusters.get(r.get("peptide"))
        if cc and r.get("cluster_top") is None:
            r["cluster_top"] = cc["cluster_top"]
            r["cluster_modes"] = cc["cluster_modes"]
            r["cluster_frac"] = round(cc["cluster_top"]
                                      / float(cc["cluster_total"]), 3)
        t = extra.get(r.get("peptide"))
        if not t:
            continue
        if r.get("panel_pose") is None:
            r["panel_pose"] = t.get("pose")
        for k in ("mm_vdw", "mm_eel", "mm_egb", "mm_esurf"):
            if r.get(k) is None and t.get(k) is not None:
                r[k] = t[k]
    results["panel_missing"] = missing
    results["panel_measurable"] = can
    return results


def _combined():
    from ._vendor import combined
    return combined


def fill_combined(results, rundir):
    """Give a run's rows the combined score and its components, in memory.

    The same idea as fill_terms, and the reason the page can offer to reorder
    a run that finished long before this metric existed: the score is
    arithmetic over columns the run already carries, so it costs nothing to
    work out when somebody looks. Nothing is written here - a run is not
    rewritten because a page was opened - and results.tsv keeps the order it
    was reported in until somebody asks for the other one.

    Also says which terms the run could be scored on, so the page can account
    for the order instead of showing a number that explains nothing.
    """
    rows = results.get("rows") if isinstance(results, dict) else None
    if not isinstance(rows, list) or not rows:
        return results
    info = report.combined_scores([r for r in rows if isinstance(r, dict)])
    if info is None:
        return results
    mod = _combined()
    terms = []
    for name in ("energy",) + tuple(t[0] for t in mod.TERMS) + ("clash",):
        if name not in (info.get("used") or {}):
            continue
        what = next((t[4] for t in mod.TERMS if t[0] == name), None)
        if name == "energy":
            what = "binding energy (%s)" % " and ".join(info["used"]["energy"])
        elif name == "clash":
            what = "steric clashes, as a penalty"
        terms.append({"name": name, "weight": info["weights"].get(name, 0.0),
                      "what": what})
    results["combined"] = {
        "terms": terms,
        "constant": info.get("constant") or [],
        "off": info.get("off") or [],
        "on_disk": os.path.isfile(os.path.join(rundir,
                                               "results_combined.tsv")),
    }
    return results


def _open_detached(cmd):
    """Start a file manager and return at once; its window is the user's."""
    import subprocess
    subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                     stdin=subprocess.DEVNULL, start_new_session=True)


# The function that starts the file manager: a name, so the tests can put
# a recorder in its place instead of opening windows.
OPENER = _open_detached


def _mtime(path):
    """When a path was last modified, in seconds, or None if it cannot be
    read."""
    try:
        return os.path.getmtime(path)
    except OSError:
        return None


#: Where counting a folder's entries stops. A folder of unpacked trajectories
#: holds thousands of files, and on a shared filesystem each listing of a
#: subfolder costs time the page spends waiting.
COUNT_MOST = 5000


def _count(folder, most=COUNT_MOST):
    """How many entries a folder holds, hidden ones left out, as the
    listing leaves them out, counted up to `most`. None when it cannot be
    read."""
    n = 0
    try:
        with os.scandir(folder) as it:
            for e in it:
                if not e.name.startswith("."):
                    n += 1
                    if n >= most:
                        break
    except OSError:
        return None
    return n


def _attachment(name):
    """A Content-Disposition for a download, safe for any file name: the
    plain form with only printable ASCII and no quotes, and the UTF-8 form
    beside it for a browser that reads it."""
    plain = "".join(c if (32 <= ord(c) < 127 and c not in '"\\') else "_"
                    for c in name) or "download"
    return ('attachment; filename="%s"; filename*=UTF-8\'\'%s'
            % (plain, urllib.parse.quote(name, safe="")))


class App:
    def __init__(self, runs, tools_roots=()):
        self.runs = os.path.abspath(os.path.expanduser(runs))
        os.makedirs(self.runs, exist_ok=True)
        self.tools_roots = [os.path.abspath(r) for r in tools_roots]
        self.token = secrets.token_urlsafe(16)
        self.lock = threading.Lock()

    # ---- the pieces the page asks for -----------------------------------
    def tools(self):
        T = toolmod.discover(roots=self.tools_roots)
        return {"tools": [{"name": n, "path": p} for n, p in T.report()],
                "sbatch": slurm.available(),
                # Local execution is reported in the header
                # with the core count, because a search sized
                # for a compute node takes far longer on a
                # laptop and the person choosing the settings
                # should know which they are on.
                "local": slurm.local_mode(),
                "cores": slurm.usable_cores(),
                # The cores on a node of the default partition, so the page's
                # Cores per task starts at something a node here can give.
                "node_cores": None if slurm.local_mode() else slurm.node_cores(),
                "roots": self.tools_roots,
                # The Slurm account and options every job on this cluster
                # needs, when activate.sh or the desktop launcher says: the
                # page fills its own boxes with them when they are empty.
                "slurm_defaults": {
                    "account": os.environ.get("ADCP_SLURM_ACCOUNT", "").strip(),
                    "extra": os.environ.get("ADCP_SLURM_EXTRA", "").strip()},
                # Required tools only. Optional ones are listed apart, so a
                # machine without AmberTools reads as able to dock.
                "missing": T.missing_required(),
                "optional_missing": [n for n, p in T.report()
                                     if not p and n not in toolmod.REQUIRED],
                "mmgbsa": T.can_rescore(),
                # Which MM-GBSA engines this machine can actually run, with
                # what each one computes. The page offers only these: an
                # engine offered but not installed produces a run that
                # rescores nothing and a table that cannot be ranked on.
                "gb_engines": self._gb_engines(),
                # The ADCP 1.1 suite, which every run docks with: its
                # folder, its version, and whether it carries its own AGFR,
                # which builds the maps. Without it nothing can dock.
                "adcp11": self._adcp11()[0],
                "adcp11_version": self._adcp11()[1],
                "suite": (adcp11_mod.home_of(self._adcp11()[0])
                          if self._adcp11()[0] else ""),
                "agfr11": self._agfr11(),
                # On a Mac, where Scripps builds no ADCP 1.1, ADFRsuite 1.0
                # docks instead (mac10); the page says so where the engine
                # is shown. Never engaged on any other system.
                "mac10": ({"engaged": True, "label": mac10_mod.LABEL,
                           "note": mac10_mod.NOTE,
                           "max_poses": mac10_mod.MAX_POSES}
                          if mac10_mod.engaged() else {"engaged": False}),
                # ADapt, a separate tool for preparing the receptor, when it
                # is installed where this can find it.
                "adapt": self._adapt(),
                # The two programs taken from the suite as well: reduce, the
                # default protonation, and MSMS, for the viewers' surfaces.
                "reduce": T.suite_program("reduce") or "",
                "msms": self._msms_said(T)}

    @staticmethod
    def _msms_said(T):
        from . import surface as surf
        return surf.describe(toolmod.msms_command(T))

    def adapt_warm(self):
        """Start ADapt in the background, once, if it can run here.

        A first start imports AGFR's libraries from a shared filesystem, and
        on a busy login node that took over a minute: long enough that a
        person pressing Open ADapt saw a blank tab and gave up. The page asks
        for this when it shows the button, so the button finds ADapt running.
        Asked for, rather than done whenever the tools are listed, so nothing
        that merely lists them starts a server.
        """
        if not self._adapt().get("available"):
            return {"warming": False}
        with _ADAPT_LOCK:
            if getattr(self, "_adapt_warming", False):
                return {"warming": True}
            self._adapt_warming = True
        threading.Thread(target=self.adapt_open, daemon=True).start()
        return {"warming": True}

    def adapt_status(self):
        """How far ADapt's start has got, for the page's loading screen: its
        process (none, starting, running or failed, with why), and whether
        it answers a request. Read without waiting on a start in progress,
        and never starts anything itself."""
        from . import adapt as adaptmod
        runner = getattr(self, "_adapt_runner", None)
        out = runner.phase() if runner is not None else {"process": "none"}
        out["answering"] = bool(out["process"] == "running"
                                and adaptmod.answering(runner.url))
        return out

    def _adapt(self):
        """Whether ADapt can be opened from here, and if not, why."""
        from . import adapt as adaptmod
        home = adaptmod.find(self.tools_roots)
        if not home:
            return {"available": False, "repo": adaptmod.REPO,
                    "why": "ADapt is not installed here."}
        if not self._adcp11()[0]:
            return {"available": False, "repo": adaptmod.REPO,
                    "why": "ADapt runs on ADCP 1.1, which is not installed here."}
        if not adaptmod.embeddable(home):
            return {"available": False, "repo": adaptmod.REPO,
                    "why": "The ADapt in %s (%s) is older than the version this "
                           "interface works with (%s or later). Update it: ./update.sh "
                           "in the ADCP_docking checkout updates ADapt too, or git -C "
                           "%s pull" % (home, adaptmod.commit(home) or "unknown",
                                        adaptmod.EMBED_COMMIT, home)}
        return {"available": True, "repo": adaptmod.REPO, "why": "",
                "commit": adaptmod.commit(home)}

    def _adapt_running(self):
        """The running ADapt's address with its token, starting it if it is
        not running. Raises RuntimeError saying why it cannot."""
        from . import adapt as adaptmod
        from . import adcp11
        home = adaptmod.find(self.tools_roots)
        exe = self._adcp11()[0]
        st = self._adapt()
        if not st.get("available"):
            raise RuntimeError(st["why"])
        with _ADAPT_LOCK:
            if getattr(self, "_adapt_runner", None) is None:
                import atexit
                self._adapt_runner = adaptmod.Runner()
                atexit.register(self._adapt_runner.stop)
        try:
            return self._adapt_runner.start(
                home, adcp11.home_of(exe), self.adapt_outdir())
        except OSError as exc:
            raise RuntimeError(str(exc))

    def adapt_outdir(self):
        """Where ADapt keeps its sessions: in the runs folder, not in the
        home directory, which is where ADapt keeps them by default."""
        return os.path.join(self.runs, "adapt_sessions")

    def adapt_open(self, body=None):
        """Start ADapt if it is not running, and return its address.

        With a receptor (body["path"]), ADapt is handed the file through its
        own api/upload and the page is opened on that session (session= in
        the URL, as its README's "Embedding ADapt" describes). ADapt reads
        PDB files only; for an mmCIF fetched by its id, ADapt is asked to
        fetch the entry itself (its api/fetch), which opens a session the
        same way, so this server knows the session whichever way it came.
        Only if that fails is the page sent to ADapt with pdb= instead. The
        answer names the session and its folder, which the page then watches
        for changes and saves.
        """
        from . import adapt as adaptmod
        body = body or {}
        try:
            url = self._adapt_running()
        except RuntimeError as exc:
            return {"error": str(exc)}
        out = {"direct": url}
        path = str(body.get("path") or "").strip()
        code = str(body.get("code") or "").strip().upper()
        session = None
        # A session this page opened before, while ADapt still holds it, is
        # opened again with what was chosen in it rather than started anew.
        again = str(body.get("session") or "").strip()
        if again:
            try:
                got = adaptmod.call(url, "api/session", {"s": again})
                session = again
                out.update(session=again, folder=got.get("folder") or "",
                           tab=got.get("current") or "",
                           pdbid=got.get("pdbid") or "")
            except RuntimeError:
                session = None
        if session:
            pass
        elif path and os.path.isfile(path):
            low = path.lower()
            if low.endswith(adaptmod.READS):
                # Named by the entry's id where the file is one, in lower
                # case as the RCSB names its files: ADFR takes the PDB id
                # from the file name and looks the entry's ligands up by it.
                name = os.path.basename(path)
                if fetch.looks_like_code(name.split(".")[0]):
                    name = name.lower()
                try:
                    got = adaptmod.upload(url, path, name)
                except (OSError, RuntimeError) as exc:
                    return {"error": "ADapt could not open %s: %s"
                                     % (os.path.basename(path), exc)}
                session = got.get("session")
                out.update(session=session, folder=got.get("folder") or "",
                           tab=got.get("current") or "",
                           pdbid=got.get("pdbid") or "")
            elif fetch.looks_like_code(code):
                session = self._adapt_fetch(url, code, out)
            else:
                out["note"] = ("ADapt reads PDB files only, and this receptor "
                               "is %s; load a PDB file in ADapt with Open PDB "
                               "or Fetch." % os.path.basename(path))
        elif fetch.looks_like_code(code):
            session = self._adapt_fetch(url, code, out)
        if session:
            self._adapt_last = {"session": session, "folder": out["folder"]}
            self._adapt_folders()[session] = out["folder"]
        # ADapt listens on a port of its own on this machine, which a browser
        # elsewhere cannot reach through the one port the launcher or an ssh
        # tunnel forwards. So the page is sent to /adapt/ on THIS server,
        # which passes each request on (Handler._adapt_proxy). ADapt's own
        # token travels with it and ADapt checks it, as it always did.
        out["url"] = adaptmod.page_url(url, session=session, pdb=out.get("pdb"))
        return out

    def _adapt_fetch(self, url, code, out):
        """An entry opened in ADapt by its id, through ADapt's own api/fetch:
        the session it made, recorded in `out`, or None, with the page sent
        to ADapt's pdb= instead so ADapt says itself what went wrong."""
        from . import adapt as adaptmod
        try:
            got = adaptmod.call(url, "api/fetch", {"code": code.upper()}, data=b"")
        except RuntimeError as exc:
            out["pdb"] = code
            out["note"] = "ADapt could not fetch %s here (%s)" % (code, exc)
            return None
        out.update(session=got.get("session"), folder=got.get("folder") or "",
                   tab=got.get("current") or "", pdbid=got.get("pdbid") or "")
        return got.get("session")

    def _adapt_folders(self):
        """Each session opened from here, with its folder."""
        if not hasattr(self, "_adapt_folder_map"):
            self._adapt_folder_map = {}
        return self._adapt_folder_map

    def _adapt_folder(self, session, folder=""):
        """A session's folder: as given, as opened from here, or as ADapt
        lists it (a structure opened inside ADapt itself)."""
        from . import adapt as adaptmod
        folder = str(folder or "").strip()
        session = str(session or "").strip()
        if folder or not session:
            return folder
        folder = self._adapt_folders().get(session, "")
        if folder:
            return folder
        runner = getattr(self, "_adapt_runner", None)
        if runner is not None and runner.running():
            folder = adaptmod.session_folder(runner.url, session)
            if folder:
                self._adapt_folders()[session] = folder
        return folder

    def adapt_state(self, session="", folder=""):
        """What has been changed in an ADapt session so far, read from the
        session.json ADapt rewrites after every edit (adapt.read_state).

        The page asks every few seconds while ADapt is shown, and once the
        number of edits grows it says that docking will use the receptor as
        set in ADapt. The saves this server asked for are marked, with the
        exact number of edits made before each, so the page can tell whether
        the receptor written last is still the receptor as set."""
        from . import adapt as adaptmod
        folder = self._adapt_folder(session, folder)
        if not folder:
            return {"found": False,
                    "why": "ADapt does not hold that session any more"}
        try:
            got = adaptmod.read_state(folder)
        except adaptmod.NoSave as exc:
            return {"found": False, "why": str(exc), "folder": folder}
        mine = getattr(self, "_adapt_written", {})
        for sv in got["saves"]:
            k = (got["folder"], sv["n"])
            if k in mine:
                sv["edits_before"] = mine[k]
                sv["written_here"] = True
        got["found"] = True
        got["session"] = str(session or "")
        return got

    def adapt_write(self, body=None):
        """Ask ADapt to write the receptor of a session as it is set there,
        the way its own receptor PDBQT button does (adapt.write_receptor),
        and read the save back as the page reads one it watched for
        (adapt.read_save): the same PDBQT, _rec_H.pdb and checksums that the
        run then records under `adapt` in job.json.

        A receptor this server wrote is used again rather than written anew
        when nothing has been changed in the session since."""
        from . import adapt as adaptmod
        body = body or {}
        session = str(body.get("session") or "").strip()
        if not session:
            return {"found": False, "why": "no ADapt session is open on this page"}
        runner = getattr(self, "_adapt_runner", None)
        if runner is None or not runner.running():
            return {"found": False,
                    "why": "ADapt is not running; open the receptor in ADapt again"}
        folder = self._adapt_folder(session, body.get("folder"))
        if not folder:
            return {"found": False,
                    "why": "ADapt does not hold that session any more; open the "
                           "receptor in ADapt again"}
        with _ADAPT_WRITE_LOCK:
            try:
                state = adaptmod.read_state(folder)
            except adaptmod.NoSave as exc:
                return {"found": False, "why": str(exc), "folder": folder}
            if not hasattr(self, "_adapt_written"):
                self._adapt_written = {}
            reuse = None
            for sv in reversed(state["saves"]):
                if self._adapt_written.get((state["folder"], sv["n"])) \
                        == state["edits"]:
                    reuse = sv["n"]
                    break
            if reuse is None:
                t0 = time.time()
                try:
                    done = adaptmod.write_receptor(runner.url, session)
                except RuntimeError as exc:
                    return {"found": False, "folder": folder,
                            "why": "ADapt could not write the receptor: %s" % exc}
                reuse = int(done["save"])
                self._adapt_written[(state["folder"], reuse)] = state["edits"]
                wrote = round(time.time() - t0, 1)
            else:
                wrote = None
            try:
                got = adaptmod.read_save(folder, reuse)
            except adaptmod.NoSave as exc:
                return {"found": False, "why": str(exc), "folder": folder}
        got.update(found=True, session=session, edits=state["edits"],
                   written=wrote is not None, write_seconds=wrote)
        self._adapt_last = {"session": session, "folder": got["folder"]}
        return got

    def adapt_saved(self, session="", folder="", save=None):
        """The receptor ADapt saved, read back from the session folder it
        wrote: its PDBQT, its _rec_H.pdb, and what session.json says was
        chosen (adapt.read_save). The folder is asked of ADapt when only the
        session is known, and taken from the last session opened from this
        page when neither is."""
        from . import adapt as adaptmod
        folder = str(folder or "").strip()
        session = str(session or "").strip()
        if not folder and session:
            runner = getattr(self, "_adapt_runner", None)
            if runner is not None and runner.running():
                try:
                    folder = adaptmod.call(runner.url, "api/session",
                                           {"s": session}).get("folder") or ""
                except RuntimeError:
                    folder = ""
            last = getattr(self, "_adapt_last", None) or {}
            if not folder and last.get("session") == session:
                folder = last.get("folder") or ""
        if not folder:
            folder = (getattr(self, "_adapt_last", None) or {}).get("folder") or ""
        if not folder:
            return {"found": False,
                    "why": "no ADapt session has been opened from this page"}
        try:
            got = adaptmod.read_save(folder, save)
        except adaptmod.NoSave as exc:
            return {"found": False, "why": str(exc), "folder": folder}
        got["found"] = True
        got["session"] = session
        return got

    def surface(self, q):
        """A molecular surface for a viewer, computed by MSMS (surface.py).

        For the site view, the receptor file as the run will keep it: the
        chains asked for, without waters, and without hetero groups unless
        they are kept. A receptor ADapt wrote is kept as ADapt wrote it, so
        the page asks for its hetero groups and waters too (het=1, water=1).
        For a run, its prepared receptor as it was docked, or one pose's
        peptide. The answer is cached in the runs folder by what it was
        computed from. Where MSMS cannot give one, the answer says why and
        the viewer draws its own.

        Every separate piece of the receptor has its outer surface: two
        chains that do not touch, or a domain whose linker is missing, are
        each surfaced (surface.compute).
        """
        from . import surface as surf
        exe = toolmod.msms_command(toolmod.discover(roots=self.tools_roots))
        one = lambda k: (q.get(k) or [""])[0]                # noqa: E731
        rundir, what = one("dir"), one("what") or "receptor"
        try:
            if rundir:
                if not os.path.isfile(jobmod.path_of(rundir)):
                    return {"source": "fallback", "why": "no such run"}
                if what == "pose":
                    text, err = self.structure_text(
                        rundir, one("peptide") or None, one("rank") or "1",
                        "pose")
                else:
                    text, err = self.structure_text(rundir, None, 1, "receptor")
                if err:
                    return {"source": "fallback", "why": err}
                import tempfile
                with tempfile.NamedTemporaryFile("w", suffix=".pdb",
                                                 delete=False) as fh:
                    fh.write(text)
                    tmp = fh.name
                try:
                    # The receptor as it was docked: a water ADapt was told
                    # to keep is part of it.
                    atoms = structure.read(tmp, keep_water=True)
                finally:
                    os.unlink(tmp)
            else:
                path = one("path")
                if not path or not os.path.isfile(path) \
                        or not structure.is_structure(path):
                    return {"source": "fallback", "why": "no such structure"}
                chains = [c for c in one("chains").split(",") if c]
                het = one("het") in ("1", "true", "yes")
                water = one("water") in ("1", "true", "yes")
                # Modified residues the preparation reads as the residue they
                # stand in for stay, as they stay in the receptor.
                from ._vendor import prepare as _prep
                part = set(structure.AA3) | set(_prep.STANDIN)
                atoms = [a for a in structure.read(path, keep_water=water)
                         if (not chains or a.chain in chains)
                         and (het or not a.hetatm or a.resname in part)]
            got = surf.cached(atoms, exe, os.path.join(self.runs, "_surfaces"))
            got["msms"] = exe
            return got
        except surf.SurfaceError as exc:
            return {"source": "fallback", "why": str(exc)}
        except Exception as exc:                                # noqa: BLE001
            return {"source": "fallback",
                    "why": "the surface could not be computed (%s)" % exc}

    def adapt_stop(self):
        runner = getattr(self, "_adapt_runner", None)
        if runner is not None:
            runner.stop()

    def _gb_engines(self):
        """Every MM-GBSA engine, described, with whether it is installed.

        Described rather than just named, because the three do not compute
        the same quantity and the person choosing has to be told so before
        they choose rather than after they compare two columns.
        """
        from ._vendor import gbengines
        out = []
        for row in gbengines.catalogue():
            if row["name"] == "schrodinger":
                # The licence, not just the program: Prime needs 8 rescoring
                # tokens a job and the modern force field has seats of its
                # own, so what is free decides what a run would actually get.
                try:
                    row["licence"] = gbengines.SCHRODINGER.plan(timeout=3)
                except Exception:                              # noqa: BLE001
                    row["licence"] = {"readable": False}
                # Where it looked, when it did not find it. "Not installed"
                # on a cluster that has the suite means a path this code
                # does not know, and saying what was tried is the only way
                # to learn which one.
                if not row.get("available"):
                    try:
                        row["search"] = gbengines.schrodinger_search_report()
                    except Exception:                          # noqa: BLE001
                        row["search"] = None
            out.append(row)
        return out

    def _adcp11(self):
        """(executable or None, version), looked up once per server."""
        if not hasattr(self, "_a11"):
            from . import adcp11
            exe = adcp11.find(self.tools_roots)
            self._a11 = (exe, adcp11.version(exe) if exe else "")
            self._a11lib = adcp11.library(exe) if exe else set()
        return self._a11

    def _agfr11(self):
        """Whether the ADCP 1.1 suite here carries its own AGFR."""
        from . import adcp11
        exe = self._adcp11()[0]
        return bool(exe and adcp11.agfr_script(exe))

    def modifications(self):
        """What the Modifications dialog may offer: ADCP 1.1's non-standard
        side chains that this installation has, grouped by parent residue in
        the page. Empty, and not available, without ADCP 1.1."""
        exe = self._adcp11()[0]
        if not exe:
            return {"available": False, "modifications": []}
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "static", "modifications.json")
        try:
            rows = json.load(open(path))["modifications"]
        except (OSError, ValueError, KeyError):
            rows = []
        lib = self._a11lib
        return {"available": True,
                "modifications": [r for r in rows if r["code"] in lib]}

    def _pepcheck(self, body):
        """The keyword arguments a peptide check needs: ADCP 1.1's syntax,
        and the side chains this installation's library has."""
        mx = int(body.get("max_length") or peptide.MAX_LENGTH)
        exe = self._adcp11()[0]
        return {"max_length": mx, "library": self._a11lib if exe else None,
                "extended": not mac10_mod.engaged()}

    def roots(self):
        """Where the browser should offer to start.

        Shared group space first, then the user's own: on a cluster the
        structures and runs people care about are usually in a group's or a
        project's space. Clusters lay these out in a few common ways, and only
        folders that exist and can be read are offered.
        """
        out, seen = [], set()

        def add(path, label):
            if not path:
                return
            real = os.path.realpath(os.path.expanduser(path))
            if real in seen or not os.path.isdir(real):
                return
            if not os.access(real, os.R_OK | os.X_OK):
                return
            seen.add(real)
            out.append({"path": real, "label": label})

        import grp
        import pwd
        try:
            names = [grp.getgrgid(g).gr_name for g in os.getgroups()]
        except (OSError, KeyError):
            names = []
        try:
            user = pwd.getpwuid(os.getuid()).pw_name
        except (OSError, KeyError):
            user = os.environ.get("USER", "")
        for g in names:
            # A user's own primary group is not a shared
            # allocation.
            if g == user:
                continue
            for top in ("/", "/scratch", "/work"):
                for kind in ("group", "groups", "project", "projects"):
                    add(os.path.join(top, kind, g), "%s (group)" % g)
        add(os.environ.get("SCRATCH"), "your scratch")
        if user:
            for kind in ("", "user", "users"):
                add(os.path.join("/scratch", kind, user), "your scratch")
        add(os.path.expanduser("~"), "home")
        add(self.runs, "runs")
        return {"roots": out}

    def make_dir(self, parent, name):
        """Create one folder inside `parent`. Returns its path.

        One level only, and the name is taken apart rather than trusted: "a/b" or
        "../elsewhere" typed into a New folder box means the first and would
        otherwise produce the second.
        """
        parent = os.path.abspath(os.path.expanduser((parent or "").strip()))
        name = (name or "").strip()
        if not os.path.isdir(parent):
            return {"error": "no such directory: %s" % parent}
        if not name:
            return {"error": "give the folder a name"}
        if name in (".", "..") or "/" in name or "\\" in name:
            return {"error": "a name, not a path: no slashes, and not .. "
                             "(you can browse to where you want it first)"}
        if name.startswith("."):
            return {"error": "a name starting with a dot is hidden from "
                             "this browser, so you would not see it again"}
        full = os.path.join(parent, name)
        if os.path.exists(full):
            return {"error": "%s already exists here" % name} \
                if not os.path.isdir(full) else {"path": full, "existed": True}
        if not os.access(parent, os.W_OK | os.X_OK):
            return {"error": "%s is not writable by you" % parent}
        try:
            os.mkdir(full)
        except OSError as exc:
            return {"error": "could not create it: %s" % exc}
        return {"path": full, "existed": False}

    def set_runs(self, path):
        """Move where new runs are written.

        The launcher does not ask for this: a user signing in from a laptop has
        no way to know which of a cluster's scratch paths to use.
        """
        raw = (path or "").strip()
        # Checked before abspath. os.path.abspath("") is the current
        # working directory (a real, writable, wrong answer) so a
        # guard that runs afterwards accepts an empty box and writes
        # runs wherever the server was started.
        if not raw:
            return {"error": "give a directory"}
        path = os.path.abspath(os.path.expanduser(raw))
        if path == "/":
            return {"error": "that is the filesystem root"}
        if os.path.exists(path) and not os.path.isdir(path):
            return {"error": "%s is a file, not a directory" % path}
        try:
            os.makedirs(path, exist_ok=True)
        except OSError as exc:
            return {"error": "could not create %s: %s" % (path, exc)}
        if not os.access(path, os.W_OK | os.X_OK):
            return {"error": "%s is not writable by you" % path}
        self.runs = path
        return {"root": self.runs, "runs": self.list_runs()["runs"]}

    #: What the file browser can be asked to list. A receptor is the common
    #: case; the target file box asks for AGFR's .trg files instead. Each is
    #: its own filter, so asking for one never loosens the other.
    BROWSE_KINDS = ("structure", "trg")

    @staticmethod
    def _listed(name, kind):
        if kind == "trg":
            return agfr_target.is_target(name)
        return structure.is_structure(name)

    def browse(self, path, dirs_only=False, kind="structure",
               every_file=False):
        """A directory listing. Structures first, then directories: the
        person is looking for a receptor, not for a filesystem. With
        kind="trg", target files instead of structures.

        every_file lists every regular file, whatever its name. Only the
        Folders view asks for that, through runs_folder, which has already
        held the path to the runs folder; /api/browse never passes it, so
        the picker still lists only what it was asked for."""
        kind = kind or "structure"
        if kind not in self.BROWSE_KINDS:
            return {"error": "the browser cannot list %r files" % kind}
        path = os.path.abspath(os.path.expanduser(path or self.runs))
        if os.path.isfile(path):
            path = os.path.dirname(path)
        if not os.path.isdir(path):
            return {"error": "not a directory: %s" % path}
        dirs, files = [], []
        try:
            for n in sorted(os.listdir(path)):
                if n.startswith("."):
                    continue
                full = os.path.join(path, n)
                if os.path.isdir(full):
                    # The time costs a second look at each folder, so only
                    # the Folders view, which shows it, pays for it.
                    dirs.append({"name": n, "path": full,
                                 "mtime": _mtime(full) if every_file
                                 else None})
                elif every_file or self._listed(n, kind):
                    # A broken symlink must not empty the
                    # listing: stat raises
                    # FileNotFoundError for a link with
                    # no target, so that error is caught
                    # for each entry.
                    try:
                        st = os.stat(full)
                    except OSError:
                        continue
                    # A pipe or a socket is not a file anyone can open
                    # or download.
                    if not stat.S_ISREG(st.st_mode):
                        continue
                    files.append({"name": n, "path": full,
                                  "size": st.st_size, "mtime": st.st_mtime})
        except PermissionError:
            return {"error": "no permission to read %s" % path}
        return {"path": path, "parent": os.path.dirname(path),
                "dirs": dirs, "files": [] if dirs_only else files,
                "writable": os.access(path, os.W_OK | os.X_OK)}

    # ---- the runs folder, walked ----------------------------------------
    #: The most entries the Folders view lists from one folder. A docking
    #: folder whose trajectories were not packed holds thousands of files,
    #: and nobody reads a listing that long in a page.
    FOLDER_MAX = 2000
    #: The largest text or structure the page shows itself. Anything larger
    #: is downloaded instead, so a log of hundreds of megabytes is never
    #: pulled into a browser tab to be drawn.
    VIEW_MAX = 8 << 20
    #: Files the page shows as text.
    TEXT_VIEW = (".txt", ".log", ".out", ".err", ".tsv", ".csv", ".json",
                 ".sbatch", ".sh", ".bib", ".ris", ".fasta", ".fa", ".seq",
                 ".md", ".dlg", ".gpf", ".dat", ".in", ".py")
    #: Files the page draws in its 3D viewer besides the structures
    #: structure.is_structure names: what docking and preparation write.
    #: The structure suffixes themselves are structure's, never a second
    #: list here.
    DOCKED_VIEW = (".pdbqt", ".pqr", ".sdf", ".mol2")
    #: A run's report: opened as the run's own results table.
    REPORT_FILES = ("results.txt", "results.tsv", "results.json",
                    "results_combined.tsv")

    def _in_runs(self, raw):
        """(real path, path relative to the runs folder), or (None, why).

        The path comes from a query string, so it is resolved before it is
        trusted. Symbolic links and ".." are followed to where they really
        lead, and only the runs folder itself or something inside it is
        accepted. A link inside the runs folder that points out of it is
        refused like any other path outside it. A relative path is taken
        from the runs folder; an absolute one must still end up inside it.
        """
        root = os.path.realpath(self.runs)
        raw = "" if raw is None else str(raw).strip()
        if "\x00" in raw:
            return None, "that is not a path"
        cand = raw or root
        if not os.path.isabs(cand):
            cand = os.path.join(root, cand)
        try:
            real = os.path.realpath(cand)
            inside = os.path.commonpath([root, real]) == root
        except (OSError, ValueError):
            return None, "that is not a path"
        if not inside:
            return None, ("only the runs folder, %s, and what is inside it "
                          "can be shown here" % self.runs)
        rel = os.path.relpath(real, root)
        return real, ("" if rel == os.curdir else rel)

    def _run_label(self, folder):
        """A run folder's name for itself, or None when it is not a run."""
        jf = jobmod.path_of(folder)
        if not os.path.isfile(jf):
            return None
        try:
            label = jobmod.load(folder).get("label")
        except (OSError, ValueError):
            label = None
        return str(label or os.path.basename(folder))

    def _view_kind(self, name, size, in_run):
        """How the page shows a file: as the run's results ("report"), in
        its 3D viewer ("structure"), as text ("text"), or not at all ("")."""
        low = name.lower()
        if in_run and low in self.REPORT_FILES:
            return "report"
        if size > self.VIEW_MAX:
            return ""
        # A gzipped structure is read by the receptor box, not by the
        # viewer, so it is downloaded.
        if ((structure.is_structure(low) and not low.endswith(".gz"))
                or low.endswith(self.DOCKED_VIEW)):
            return "structure"
        if low.endswith(self.TEXT_VIEW):
            return "text"
        return ""

    def runs_folder(self, path=""):
        """One folder of the runs folder, for the Folders view: the folders
        in it, then every file in it, with sizes and modification times, and
        whether it, or any folder in it, is a run.

        The listing is the file picker's own (browse), asked for every file
        instead of structures. Nothing outside the runs folder is listed:
        the folder asked for is held to it first, and an entry that is a
        link leading out of it is left out.
        """
        real, rel = self._in_runs(path)
        if real is None:
            return {"error": rel}
        if not os.path.isdir(real):
            return {"error": "there is no folder %s in the runs folder"
                             % (rel or os.curdir)}
        got = self.browse(real, every_file=True)
        if got.get("error"):
            return got
        root = os.path.realpath(self.runs)

        def inside(p):
            try:
                r = os.path.realpath(p)
                return r if os.path.commonpath([root, r]) == root else None
            except (OSError, ValueError):
                return None

        here_run = self._run_label(real)
        dirs, files = [], []
        for d in got["dirs"]:
            r = inside(d["path"])
            if r is None:
                continue
            drel = os.path.relpath(r, root)
            dirs.append({"name": d["name"], "rel": drel,
                         "path": os.path.join(self.runs, drel),
                         "mtime": d.get("mtime"), "items": _count(r),
                         "run": self._run_label(r)})
        for f in got["files"]:
            r = inside(f["path"])
            if r is None:
                continue
            files.append({"name": f["name"],
                          "rel": os.path.relpath(r, root),
                          "size": f["size"], "mtime": f.get("mtime"),
                          "view": self._view_kind(f["name"], f["size"],
                                                  bool(here_run))})
        more = max(0, len(dirs) + len(files) - self.FOLDER_MAX)
        if more:
            dirs = dirs[:self.FOLDER_MAX]
            files = files[:max(0, self.FOLDER_MAX - len(dirs))]
        crumbs = [{"name": os.path.basename(self.runs.rstrip(os.sep))
                   or self.runs, "rel": ""}]
        at = ""
        for part in (rel.split(os.sep) if rel else []):
            at = os.path.join(at, part) if at else part
            crumbs.append({"name": part, "rel": at})
        return {"root": self.runs, "rel": rel,
                "path": os.path.join(self.runs, rel) if rel else self.runs,
                "parent": (os.path.dirname(rel) if rel else None),
                "crumbs": crumbs, "run": here_run,
                "dirs": dirs, "files": files, "more": more,
                "count_most": COUNT_MOST}

    def run_file(self, path):
        """(real path, path inside the runs folder, status) for a file the
        Folders view may send, or (None, why, status) when it may not."""
        real, rel = self._in_runs(path)
        if real is None:
            return None, rel, 403
        if not os.path.isfile(real):
            return None, ("there is no file %s in the runs folder"
                          % (rel or os.curdir)), 404
        return real, rel, 200

    def receptor(self, path):
        """Chains and residues, so the page can offer them."""
        if not path or not os.path.isfile(path):
            return {"error": "no such file: %s" % path}
        if path.lower().endswith(".pdbqt"):
            return {"error": PDBQT_REFUSED % os.path.basename(path)}
        try:
            atoms = structure.read(path)
        except Exception as exc:                            # noqa: BLE001
            return {"error": "could not read the structure: %s" % exc}
        if not atoms:
            return {"error": "no atoms in %s; is it a structure file?"
                            % os.path.basename(path)}
        ch = structure.chains(atoms)
        return {"path": path, "atoms": len(atoms), "chains": ch,
                "residues": {c["id"]: structure.residues(atoms, c["id"])
                             for c in ch if c["kind"] == "protein"},
                # What a box could be built around instead of residues: the
                # peptides and hetero groups the file holds, waters left out.
                "ligands": structure.ligands(atoms)}

    def target_file(self, path):
        """A target file's receptor and box, for the site view.

        Held to the check the run itself applies (agfr_target.preview), so a
        file shown here is one the run will take, and one it would refuse is
        refused now, with the reason, rather than at the start of the run.
        """
        raw = (path or "").strip()
        if not raw:
            return {"error": "give the path to a target file"}
        path = os.path.abspath(os.path.expanduser(raw))
        name = os.path.basename(path)
        if not agfr_target.is_target(path):
            return {"error": "%s is not a target file. AGFR and AGFRGUI write "
                             "targets as .trg files; a structure goes in the "
                             "box above instead." % name}
        if not os.path.isfile(path):
            return {"error": "There is no file at %s." % path}
        try:
            pdb, centre, size = agfr_target.preview(path)
        except Exception as exc:                            # noqa: BLE001
            return {"error": "%s cannot be used as a target file: %s. Make it "
                             "again in AGFR or AGFRGUI, or leave this box "
                             "empty and prepare the receptor here."
                             % (name, str(exc).rstrip("."))}
        from . import target as _t
        atoms = [L for L in pdb.splitlines() if L.startswith(("ATOM", "HETATM"))]
        pts = _t.translation_points(path)
        return {"path": path, "name": name, "pdb": pdb,
                "centre": [round(float(c), 3) for c in centre],
                "size": [round(float(v), 3) for v in size],
                "atoms": len(atoms),
                "chains": sorted({L[21] for L in atoms if L[21].strip()}),
                "maps": _t.maps_in(path),
                # The AutoSite points the docking starts from, as agfr put
                # them in the file: drawn in the view so the pocket can be
                # seen. Thinned past a few thousand, which a view cannot
                # tell apart anyway.
                "points": pts[::max(1, len(pts) // 3000)],
                "n_points": len(pts)}

    def fetch_code(self, code):
        """Pull an entry from the RCSB into the cache and read it back."""
        cache = os.path.join(self.runs, "_structures")
        try:
            path = fetch.fetch(code, cache, log=lambda *a: None)
        except fetch.FetchError as exc:
            return {"error": str(exc)}
        out = self.receptor(path)
        out["fetched"] = True
        return out

    @staticmethod
    def receptor_cost(path, chains=()):
        """(heavy atoms, distinct non-water heteroatom kinds), counted over
        the chains that will be kept (sizing.receptor_cost)."""
        return sizing.receptor_cost(path, chains)

    def estimate_run(self, body):
        """What to ask the queue for, from what the form currently holds."""
        seqs = [peptide.clean(p.get("sequence", ""))
                for p in body.get("peptides") or []]
        # Residues, not characters: N<BHD>PIS&DVD is seven residues.
        seqs = [len(peptide.residues(s)) for s in seqs if s]
        if not seqs:
            return {"error": "add a peptide first"}
        # What the receptor itself costs, read from the file rather
        # than assumed. With no receptor chosen the floor stands.
        atoms, kinds = 0, 0
        rec = body.get("receptor")
        if rec and os.path.isfile(rec):
            atoms, kinds = self.receptor_cost(rec, body.get("chains") or [])
        box = body.get("box") or None
        return sizing.suggest(
            seqs, int(body.get("replicas") or 50),
            int(body.get("steps") or 2_500_000),
            cores=int(body.get("cores") or sizing.CORES_DEFAULT),
            mmgbsa=bool(body.get("mmgbsa")),
            mm_poses=int(body.get("mmgbsa_poses") or 5),
            throttle=int(body.get("throttle") or 0),
            atoms=atoms, hetatm_kinds=kinds,
            fill_loops=bool(body.get("fill_loops")), box=box,
            options=self.cost_options(body),
            prep_booked=self._prep_booked(body))

    def _prep_booked(self, body):
        """What the prepare job will book, in minutes, worked out the way
        the submission works it out (pipeline.prep_time), so the total the
        page states is what the queue is actually asked for. None when it
        cannot be worked out, and the estimate's own booking stands."""
        try:
            j = jobmod.default(**{k: v for k, v in body.items()
                                  if k in jobmod.default()})
            j["tools_roots"] = self.tools_roots
            # Capped at the partition's limit, as the submission caps it.
            booked = slurm.within_limit(pipeline.prep_time(j), j.get("partition"),
                                        j.get("sbatch_extra"))
            h, m, s = (int(x) for x in booked.split(":"))
            return h * 60 + m + (1 if s else 0)
        except Exception:                                      # noqa: BLE001
            return None

    def cost_options(self, body):
        """The settings that change what a run costs, as the page holds
        them."""
        j = jobmod.default(**{k: v for k, v in body.items()
                              if k in jobmod.default()})
        j["tools_roots"] = self.tools_roots
        return {"autosite": str(j.get("autosite") or "1.1"),
                "receptor_gradient": bool(j.get("receptor_gradient")),
                "site_mode": str(j.get("site_mode") or "residues"),
                "cyclic": str(j.get("cyclic") or "none"),
                # Each peptide's ring, its own or the run's, for the
                # peptides estimate_run prices: those with a sequence, in
                # the same order.
                "rings": [jobmod.ring_of(j, p) for p in j.get("peptides") or []
                          if isinstance(p, dict)
                          and peptide.clean(p.get("sequence", ""))],
                "helical": j.get("helical", 0),
                "poses": int(j.get("poses") or 10),
                "keep_trajectories": bool(j.get("keep_trajectories")),
                "box_max": float(j.get("box_max") or 30.0),
                "protonation": jobmod.protonation_of(j),
                "from_adapt": bool((j.get("adapt") or {}).get("folder"))}

    def sequence(self, path, chain):
        """The chain as residues, for the strip people click residues in."""
        if not path or not os.path.isfile(path):
            return {"error": "no such file"}
        atoms = structure.read(path)
        return {"chain": chain,
                "residues": structure.residues(atoms, chain or None)}

    def check_run(self, rundir):
        """None when this is a run directory that may be read, else why.

        The path comes from a query string, so it is checked rather than trusted:
        it must be a directory and must contain a job.json.
        """
        if not rundir:
            return "no run given"
        if not os.path.isdir(rundir):
            return "no such directory: %s" % rundir
        if not os.path.isfile(jobmod.path_of(rundir)):
            return ("%s is not a run directory: it holds no job.json"
                    % rundir)
        return None

    def _inside_runs(self, rundir):
        """Whether a path is a run folder inside the runs folder itself."""
        root = os.path.realpath(self.runs)
        real = os.path.realpath(rundir)
        return real != root and real.startswith(root + os.sep)

    def delete_run(self, rundir):
        """Remove a finished run's folder.

        Only a run folder inside the runs folder, and only when none of its
        jobs is queued or running: deleting under a running job leaves it
        writing into a folder that is gone.
        """
        bad = self.check_run(rundir)
        if bad:
            return {"error": bad}
        if not self._inside_runs(rundir):
            return {"error": "only runs inside %s can be deleted here"
                             % self.runs}
        j = jobmod.load(rundir)
        live = []
        for k, v in (j.get("jobs") or {}).items():
            word = self.summarise_state(slurm.state(v, rundir=rundir))
            if word in ("running", "queued"):
                live.append("%s %s" % (k, v))
        if live:
            return {"error": "this run still has jobs %s (%s); cancel them "
                             "or let them finish first"
                             % ("running" if len(live) == 1 else "queued or "
                                "running", ", ".join(live))}
        import shutil
        shutil.rmtree(rundir)
        return {"deleted": rundir}

    def open_folder(self, rundir):
        """Open a run's folder in the file manager, where that is possible.

        Only when the interface runs on the computer in front of the user
        (adcp-dock local, or a Mac or Windows machine): on a cluster or in
        Colab the folder is on another machine, and the page shows its path
        and how to reach it instead.
        """
        import platform
        import subprocess
        bad = self.check_run(rundir)
        if bad:
            return {"error": bad}
        path = os.path.realpath(rundir)
        host = platform.node()
        if os.environ.get("COLAB_RELEASE_TAG"):
            return {"opened": False, "path": path, "where": "colab",
                    "hint": "In Colab, open the Files panel (the folder icon "
                            "on the left) and go to this folder. The "
                            "notebook's Download a run step saves a run as a "
                            "zip."}
        system = platform.system()
        desktop = (os.environ.get("ADCP_DESKTOP") == "1"
                   or system in ("Darwin", "Windows"))
        if desktop:
            cmd = None
            if system == "Darwin":
                cmd = ["open", path]
            elif system == "Windows":
                cmd = ["explorer", path]
            elif "microsoft" in platform.release().lower():
                try:
                    win = subprocess.run(["wslpath", "-w", path],
                                         capture_output=True, text=True,
                                         timeout=10).stdout.strip()
                    cmd = ["explorer.exe", win] if win else None
                except (OSError, subprocess.SubprocessError):
                    cmd = None
            elif os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"):
                cmd = ["xdg-open", path]
            if cmd:
                try:
                    OPENER(cmd)
                    return {"opened": True, "path": path}
                except OSError as exc:
                    return {"opened": False, "path": path, "where": "local",
                            "hint": "The file manager could not be opened "
                                    "(%s); the folder is below." % exc}
        return {"opened": False, "path": path, "where": "server",
                "hint": "This run is on %s, the machine running the "
                        "interface. Open it there in a terminal (cd to the "
                        "folder), or with a file transfer program such as "
                        "WinSCP, FileZilla or Cyberduck connected to that "
                        "machine." % (host or "the server")}

    def run_sizes(self, rundir):
        """What each download would cost, so it can be said on the button."""
        from . import bundle as _b
        bad = self.check_run(rundir)
        if bad:
            return {"error": bad}
        try:
            s = _b.sizes(rundir)
        except OSError as exc:
            return {"error": str(exc)}
        return {"results": {"files": s["results"][0], "bytes": s["results"][1],
                            "human": _b.human_size(s["results"][1])},
                "all": {"files": s["all"][0], "bytes": s["all"][1],
                        "human": _b.human_size(s["all"][1])}}

    def report_problem(self, body):
        """Mail the logs to the maintainers. Returns what happened.

        The logs rather than the results: a problem report does not need the poses,
        and a large attachment is one a mail server refuses. Where the machine
        cannot send mail this does not fail silently.
        """
        from . import bundle as _b
        rundir = (body.get("rundir") or "").strip()
        bad = self.check_run(rundir)
        if bad:
            return {"error": bad}
        note = (body.get("note") or "").strip()
        atts = _b.logs_for(rundir)
        extra = dict(body.get("extra") or {})
        if (body.get("email") or "").strip():
            extra["email"] = body["email"].strip()
        text = _b.report_body(rundir, note, extra=extra)
        if not _b.can_send():
            return {"sent": False,
                    "error": "this machine has no way to send mail.",
                    "fallback": ("Download the logs with the button beside "
                                 "this one and send them to %s. The note "
                                 "you typed is in the message below."
                                 % REPORT_TO),
                    "body": text}
        who = (body.get("email") or "").strip()
        ok, why = _b.send(REPORT_TO, "GUI to ADCP Docking: a run went wrong",
                          text, attachments=atts, cc=REPORT_CC,
                          reply_to=who or None)
        if not ok:
            return {"sent": False, "error": why,
                    "fallback": ("Download the logs with the button beside "
                                 "this one and send them to %s." % REPORT_TO),
                    "body": text}
        return {"sent": True,
                "note": "%s, with %d log file(s)." % (why, len(atts))}

    def residues(self, rundir, peptide):
        """The per residue map of one peptide, for the viewer.

        Read from metrics/, which the interface panel writes: every receptor
        residue the pose touches and every residue of the peptide, with the
        kind of contact each makes and its share of the MM-GBSA energy. A run
        whose panel was never measured has no map, and says so.
        """
        if not rundir or not os.path.isfile(jobmod.path_of(rundir)):
            return {"error": "no such run"}
        key = jobmod.key(peptide or "")
        tsv = os.path.join(rundir, "metrics", "%s_residues.tsv" % key)
        if not os.path.isfile(tsv):
            return {"rows": [], "measured": False}
        rows = []
        with open(tsv) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                # Maps written before the peptide side was added have no side
                # column; every row in them is the receptor's.
                r.setdefault("side", "receptor")
                for k in ("dg", "d_min"):
                    try:
                        r[k] = float(r[k]) if r.get(k) not in ("", None) \
                            else None
                    except ValueError:
                        r[k] = None
                rows.append(r)
        pose = None
        js = os.path.join(rundir, "metrics", "%s.json" % key)
        if os.path.isfile(js):
            try:
                pose = json.load(open(js)).get("best_pose")
            except (OSError, ValueError):
                pose = None
        return {"rows": rows, "measured": True, "pose": pose}

    def panel_run(self, rundir):
        """Measure the interface panel for a run docked before it existed.

        Submitted to the scheduler like every other stage: it runs MMPBSA.py
        once per pose and rolls a probe over each complex, which is not work
        for the machine serving this page. Nothing is docked again.
        """
        if not rundir or not os.path.isfile(jobmod.path_of(rundir)):
            return {"error": "no such run"}
        busy = self._tidying(rundir)
        if busy:
            return {"error": busy}
        try:
            return {"jobs": pipeline.submit_panel(rundir, log=lambda *a: None)}
        except SystemExit as exc:
            return {"error": str(exc)}
        except Exception as exc:                                # noqa: BLE001
            return {"error": "could not submit the measurement: %s" % exc}

    def rerank_run(self, rundir):
        """Write a finished run's combined order into the run folder.

        Unlike the panel and the rescoring this needs no scheduler and no
        compute: the measurements are already in the run's own files and the
        score is arithmetic over them, so it is done here and finishes before
        the page has finished asking.

        results.tsv is NOT touched. Numbers from these runs have been quoted
        in reports and in conversation, and a table that says something
        different the next time somebody opens it is worse than no table. The
        combined order goes in results_combined.tsv beside it, every term in
        its own column, and the page can show either.
        """
        if not rundir or not os.path.isfile(jobmod.path_of(rundir)):
            return {"error": "no such run"}
        rj = os.path.join(rundir, "results.json")
        if not os.path.isfile(rj):
            return {"error": "this run has no results yet"}
        try:
            rows = json.load(open(rj)).get("rows") or []
        except (OSError, ValueError) as exc:
            return {"error": "could not read the run's results (%s)" % exc}
        rows = [r for r in rows if isinstance(r, dict)]
        if not rows:
            return {"error": "this run has no results yet"}
        mod = _combined()
        was = [r.get("peptide") for r in rows]
        info = report.combined_scores(rows)
        if info is None:
            return {"error": "only one measurement is available for every "
                             "peptide, so there is nothing to combine - "
                             "measure the interface first."}
        done, _i = mod.score(rows)
        spans = mod.spans(done, mod.poses_from(rundir), info=info)
        pairs = mod.covers(done, spans) if spans else {}
        for r in done:
            r["not_separated"] = ", ".join(
                pairs.get(r.get("sequence") or r.get("peptide")) or [])
        out = os.path.join(rundir, "results_combined.tsv")
        try:
            mod.write_table(done, out, info, was, spans)
        except Exception as exc:                                # noqa: BLE001
            return {"error": "could not write the table: %s" % exc}
        groups = mod.unseparated(done, spans) if spans else []
        now = [r.get("peptide") for r in done]
        # The caution travels with the answer rather than waiting in a file
        # somebody may never open: on every run measured so far the pose
        # spread covers most of the order, and a reordered table that does
        # not say so reads far more settled than it is.
        return {"wrote": out, "order": now,
                "moved": sum(1 for a, b in zip(was, now) if a != b),
                "chains": [[x.get("peptide") for x in g] for g in groups]}

    def rescore_run(self, rundir, engine):
        """Score a finished run's poses again with another MM-GBSA engine.

        Submitted to the scheduler, like the panel: it minimises every pose,
        which is not work for the machine serving this page. The run's own
        values, the ones it is ranked on, are left alone; the new engine's
        numbers appear beside them as their own columns once the report
        chained behind it has run.
        """
        from ._vendor import gbengines
        if not rundir or not os.path.isfile(jobmod.path_of(rundir)):
            return {"error": "no such run"}
        engine = str(engine or "").strip().lower()
        if engine not in gbengines.ORDER:
            return {"error": "choose an engine: %s"
                             % ", ".join(gbengines.ORDER)}
        eng = gbengines.get(engine)
        if not eng.available():
            return {"error": "%s is not installed here" % eng.label}
        busy = self._tidying(rundir)
        if busy:
            return {"error": busy}
        try:
            return {"jobs": pipeline.submit_rescore(rundir, engine,
                                                    log=lambda *a: None)}
        except SystemExit as exc:
            return {"error": str(exc)}
        except Exception as exc:                                # noqa: BLE001
            return {"error": "could not submit the rescoring: %s" % exc}

    def resume_run(self, rundir):
        if not rundir or not os.path.isfile(jobmod.path_of(rundir)):
            return {"error": "no such run"}
        busy = self._tidying(rundir)
        if busy:
            return {"error": busy}
        try:
            return pipeline.resume(rundir, log=lambda *a: None)
        except SystemExit as exc:
            return {"error": str(exc)}
        except Exception as exc:                                # noqa: BLE001
            return {"error": "could not resubmit: %s" % exc}

    def pockets(self, body):
        """AutoSite's pockets for the site as the page has it: inside the box
        the site check builds, or over the whole receptor for a box built
        around a pocket (pockets.py). Seconds on the server, cached."""
        from . import pockets as pk
        exe = self._adcp11()[0]
        if not exe:
            return {"error": "finding pockets needs AutoSite 1.1, which comes "
                             "with the ADCP 1.1 suite; it is not installed here"}
        path = str(body.get("receptor") or "").strip()
        mode = str(body.get("site_mode") or "residues")
        drop, centre, size = [], None, None
        if mode != "pocket":
            box = self.site(body)
            if box.get("error"):
                return box
            centre, size = box["centre"], box["size"]
            if mode == "ligand":
                p = structure.parse_ligand(body.get("site_ligand"))
                if p and p[0] == "chain":
                    drop = [p[1]]
        try:
            found = pk.find(path, os.path.join(self.runs, "_pockets"), exe,
                            chains=body.get("chains") or [], drop=drop,
                            keep_het=bool(body.get("keep_hetatms")),
                            centre=centre, size=size,
                            peptide=str(body.get("pocket_score") or "peptide") == "peptide")
        except pk.PocketError as exc:
            return {"error": str(exc)}
        out = pk.for_page(found)
        out["mode"] = mode
        return out

    def site(self, body):
        """Where the box would go, before anything is submitted."""
        path = body.get("receptor")
        if not path or not os.path.isfile(path):
            return {"error": "choose a receptor first"}
        all_atoms = structure.read(path)
        keep = body.get("chains") or []
        mode = str(body.get("site_mode") or "residues")
        if mode == "pocket":
            # AutoSite runs on the prepared receptor, on a compute node, so
            # there is no box to draw before the run has prepared it.
            return {"mode": "pocket", "pending": True,
                    "note": "Without a choice below, AutoSite places this box "
                            "when the receptor has been prepared, on a compute "
                            "node, around the %s. Find pockets to see them now, "
                            "ranked, and choose."
                    % ("best scoring pocket it finds"
                       if (body.get("site_pocket") or "best") == "best"
                       else "pockets it finds, all of them")}
        exclude = set()
        contacts = []
        in_receptor = False
        if mode == "ligand":
            spec = str(body.get("site_ligand") or "").strip()
            if not spec:
                return {"error": "choose the ligand or peptide the box should "
                                 "be built around"}
            # With ADapt's receptor the ligand is read from the structure
            # ADapt was opened on: ADapt may have left it out of its receptor.
            src = str(body.get("site_ligand_from") or "").strip()
            from_adapt = bool(src) and src != path and os.path.isfile(src)
            if from_adapt:
                rec_atoms = all_atoms
                all_atoms = structure.read(src)
                in_receptor = bool(structure.ligand_atoms(rec_atoms, spec))
            probs = structure.ligand_problems(
                all_atoms, spec, [] if from_adapt else keep,
                False if from_adapt else bool(body.get("keep_hetatms")))
            if probs:
                return {"error": probs[0]}
            lig = structure.ligand_atoms(all_atoms, spec)
            p = structure.parse_ligand(spec)
            if p and p[0] == "chain":
                exclude = {p[1]}
            chains_rec = keep or sorted({a.chain for a in all_atoms} - exclude)
            centre, size, note = structure.ligand_box(
                lig, padding=float(body.get("box_padding", 8.0)),
                minimum=float(body.get("box_min", 20.0)),
                maximum=float(body.get("box_max", 30.0)))
            contacts = structure.contacts(all_atoms, lig,
                                          chains_keep=chains_rec)
            hs = contacts
            keep = chains_rec
            site_pts = [(a.x, a.y, a.z) for a in lig if not a.is_hydrogen()]
        else:
            atoms = ([a for a in all_atoms if a.chain in keep] if keep
                     else all_atoms)
            hs = structure.parse_hotspots(body.get("hotspots"),
                                          (keep or [None])[0])
            if not hs:
                return {"error": "no residues chosen"}
            centre, size, note = structure.box_for(
                atoms, hs, padding=float(body.get("box_padding", 8.0)),
                minimum=float(body.get("box_min", 20.0)),
                maximum=float(body.get("box_max", 30.0)))
            want = {(c, int(n)) for c, n in hs}
            site_pts = [(a.x, a.y, a.z) for a in atoms
                        if not a.is_hydrogen()
                        and ((a.chain, a.resnum) in want
                             or (None, a.resnum) in want)]
        if centre is None:
            return {"error": note}
        # Where the site and its padding are wider than the largest box,
        # the box is held to that size and part of the site falls outside
        # it, which nothing else on the page would say.
        clip = structure.clipped(site_pts, float(body.get("box_padding", 8.0)),
                                 float(body.get("box_max", 30.0)))
        # Which other chains the box reaches, so that the decision
        # about them is an informed one. A box around an interface
        # groove routinely covers the chain on the far side.
        others = []
        try:
            from ._vendor import receptor as _rec
            # The chains the run will load: the ticked ones, or every
            # chain when none is ticked. Only chains outside that set are
            # "other"; with none ticked there are none.
            loaded = set(keep) if keep else {a.chain for a in all_atoms}
            found = _rec.chains_in_box(path, centre, size,
                                       exclude=loaded | exclude,
                                       min_atoms=20)
            # Both numbers are reported: the box may reach a few
            # hundred atoms of another chain, but including it
            # loads the whole chain, since half a chain is not a
            # molecule.
            total = {}
            for a in all_atoms:
                total[a.chain] = total.get(a.chain, 0) + 1
            others = [{"chain": c, "atoms": n, "chain_atoms": total.get(c, n)}
                      for c, n in sorted(found.items(), key=lambda x: -x[1])]
        except Exception:                                       # noqa: BLE001
            others = []
        out = {"centre": [round(c, 3) for c in centre],
               "size": list(size), "note": note, "n": len(hs),
               "other_chains": others, "mode": mode,
               "clipped": [{"axis": a, "span": sp, "padding": pd}
                           for a, sp, pd in clip],
               "box_padding": float(body.get("box_padding", 8.0)),
               "ligand_in_receptor": in_receptor,
               "box_max": float(body.get("box_max", 30.0))}
        if mode == "ligand":
            # The residues the ligand touches: the site every pose is then
            # measured against.
            out["contacts"] = ["%s:%d" % (c, n) for c, n in contacts]
        return out

    def peptides(self, body):
        """Live validation. Positions come back so the box can underline.

        Each row is checked with its own ring, or the run's when the row
        has none, and says which ring that is ("cyclic")."""
        kw = self._pepcheck(body)
        out = []
        run = {"cyclic": str(body.get("cyclic") or "none")}
        for row in body.get("peptides") or []:
            seq, probs = peptide.check(row.get("sequence", ""), **kw)
            # A ring this peptide cannot close, said at the box: ADCP would
            # dock it as a linear peptide and say nothing. And a ring the
            # row names that is not one of the three, said the same way.
            mode = jobmod.ring_of(run, row)
            if mode not in peptide.CYCLIC:
                probs = list(probs) + [peptide.ring_problem(row.get("cyclic"))]
            else:
                ring = peptide.cyclic_problem(seq, mode) if seq else None
                if ring:
                    probs = list(probs) + [ring]
            out.append({"name": row.get("name", ""), "sequence": seq,
                        "length": len(peptide.residues(seq)),
                        "cyclic": mode,
                        "ok": not probs,
                        "problems": [p.as_dict() for p in probs]})
        return {"peptides": out,
                "ready": bool(out) and all(p["ok"] for p in out)}

    def paste(self, body):
        """A pasted block or FASTA, split into rows the page can show.

        A row that names its ring carries it ("cyclic"); one that does not
        carries "", and is docked as the run is set."""
        rows = peptide.parse_rows(body.get("text", ""), **self._pepcheck(body))
        return {"peptides": [{"name": r["name"], "sequence": r["sequence"],
                              "cyclic": r["cyclic"],
                              "ok": not r["problems"],
                              "problems": [x.as_dict() for x in r["problems"]]}
                             for r in rows]}

    def upload(self, body):
        """A file chosen in the browser, written into the runs directory.

        The page is served from a login node through the launcher's tunnel, so the
        browser is on the user's own machine and the cluster cannot see its disk.
        The browser reads the file and sends the contents here.

        Structures and peptide lists only, under the runs directory, and by
        basename: a name carrying a path would otherwise choose where it landed.
        """
        name = os.path.basename((body.get("name") or "").strip())
        text = body.get("text")
        if not name or text is None:
            return {"error": "no file was sent"}
        # Dotfiles and names that are only punctuation are not what
        # anyone meant to upload.
        safe = "".join(c if (c.isalnum() or c in "-_.") else "_"
                       for c in name).lstrip(".")[:120]
        if not safe:
            return {"error": "that filename cannot be used here"}
        lists = (".txt", ".csv", ".fasta", ".fa", ".seq")
        if safe.lower().endswith(".pdbqt"):
            return {"error": PDBQT_REFUSED % safe}
        if not (structure.is_structure(safe) or safe.lower().endswith(lists)):
            return {"error": "%s is not a structure or a peptide list. "
                             "Structures: .pdb, .ent, .cif or .mmcif; lists: "
                             "%s" % (safe, ", ".join(lists))}
        data = text.encode("utf-8", "replace")
        if len(data) > MAX_UPLOAD:
            return {"error": "that file is %.0f MB; the limit is %d MB"
                             % (len(data) / 1e6, MAX_UPLOAD >> 20)}
        with self.lock:
            d = os.path.join(self.runs, "uploads")
            try:
                os.makedirs(d, exist_ok=True)
            except OSError as exc:
                return {"error": "could not make %s: %s" % (d, exc)}
            # A timestamp rather than an overwrite: two people
            # uploading receptor.pdb, or one person uploading a
            # correction, must not replace a file a queued job is
            # about to read.
            import time
            stem, ext = os.path.splitext(safe)
            stamp = time.strftime("%Y%m%d_%H%M%S")
            # The stamp resolves to the second, and two uploads
            # within one second are ordinary, a corrected file
            # sent straight after the first, so the name is made
            # unique rather than assumed to be.
            full = os.path.join(d, "%s_%s%s" % (stem, stamp, ext))
            n = 1
            while os.path.exists(full):
                full = os.path.join(d, "%s_%s_%d%s" % (stem, stamp, n, ext))
                n += 1
            try:
                with open(full, "wb") as fh:
                    fh.write(data)
            except OSError as exc:
                return {"error": "could not write it: %s" % exc}
        return {"path": full, "name": safe, "bytes": len(data)}

    # ---- runs ------------------------------------------------------------
    def create(self, body):
        with self.lock:
            label = (body.get("label") or "run").strip() or "run"
            safe = "".join(c if c.isalnum() or c in "-_" else "_"
                           for c in label)[:40]
            import time
            rundir = os.path.join(self.runs, "%s_%s"
                                  % (safe, time.strftime("%Y%m%d_%H%M%S")))
            j = jobmod.default(**{k: v for k, v in body.items()
                                  if k in jobmod.default()})
            j["label"] = label
            j["tools_roots"] = self.tools_roots
            # On a Mac the run docks with ADCP 1.0, and says so.
            jobmod.for_this_machine(j)
            # A row left "as the run" keeps no ring of its own, so the run's
            # choice applies to it; a row's own is written the one way.
            j["peptides"] = jobmod.tidy_peptides(j.get("peptides"))
            problems = jobmod.take_adapt(j)
            problems += jobmod.validate(j) if not problems else []
            if j.get("mmgbsa") and not toolmod.discover(
                    roots=self.tools_roots).can_rescore():
                problems.append("MM-GBSA needs AmberTools (tleap, sander, "
                                "MMPBSA.py), which was not found. Untick "
                                "MM-GBSA, or install AmberTools.")
            if problems:
                return {"error": "this run is not ready", "problems": problems}
            # Pockets chosen on the page: their points go into the run, and
            # for a box built around a pocket the box is set from them.
            sp = body.get("site_pockets") or {}
            j["site_pockets"], j["site_points"], j["site_box"] = {}, "", None
            if sp.get("key") and sp.get("ranks"):
                from . import pockets as pk
                try:
                    sel, pts, found = pk.chosen(os.path.join(self.runs, "_pockets"),
                                                sp["key"], sp["ranks"])
                except pk.PocketError as exc:
                    return {"error": str(exc)}
                pk.save_points(pts, os.path.join(rundir, "target", "site_points.npy"))
                j["site_points"] = os.path.join("target", "site_points.npy")
                j["site_pockets"] = {
                    "ranks": [p["rank"] for p in sel], "points": len(pts),
                    "scored": found.get("scored"),
                    "searched": "box" if found.get("box") else "receptor",
                    "table": [{k: p[k] for k in ("rank", "score", "energy", "points",
                                                 "radius", "buriedness")}
                              for p in found["pockets"][:20]]}
                if j.get("site_mode") == "pocket":
                    c, s, _clip = pk.box_around(
                        pts, float(j.get("box_padding", 8.0)),
                        float(j.get("box_min", 20.0)), float(j.get("box_max", 30.0)))
                    j["site_box"] = {"centre": c, "size": s}
            jobmod.save(rundir, j)
            try:
                ids = pipeline.submit(rundir, j, log=lambda *a: None)
            except SystemExit as exc:
                return {"error": str(exc), "rundir": rundir}
            except Exception as exc:                        # noqa: BLE001
                return {"error": "could not submit: %s" % exc, "rundir": rundir}
            return {"rundir": rundir, "jobs": ids}

    def list_runs(self):
        out = []
        for n in sorted(os.listdir(self.runs), reverse=True):
            d = os.path.join(self.runs, n)
            if not os.path.isfile(jobmod.path_of(d)):
                continue
            try:
                j = jobmod.load(d)
            except (OSError, ValueError):
                continue
            out.append({"rundir": d, "label": j.get("label") or n,
                        "created": j.get("created"), "state": j.get("state"),
                        "peptides": len(j.get("peptides") or []),
                        "mmgbsa": bool(j.get("mmgbsa"))})
        return {"runs": out[:100], "root": self.runs}

    def structure_text(self, rundir, peptide=None, rank=1, what="complex"):
        """The receptor, a pose, or both as one PDB.

        The viewer asks for a complex by default, because a pose on its own is a
        peptide floating in space.
        """
        if not rundir or not os.path.isfile(jobmod.path_of(rundir)):
            return None, "no such run"
        j = jobmod.load(rundir)
        rec = (j.get("receptor_info") or {}).get("prepared_pdb")
        if not rec or not os.path.isfile(rec):
            return None, "the receptor has not been prepared yet"
        if what == "receptor" or not peptide:
            return open(rec).read(), None
        pose = os.path.join(rundir, "docking", jobmod.key(peptide),
                            "out_ranked_%d.pdb" % int(rank))
        if not os.path.isfile(pose):
            return None, "no pose %s for %s" % (rank, peptide)
        if what == "pose":
            return open(pose).read(), None
        # Receptor first, peptide as chain Z: the viewer
        # distinguishes them by chain, so the letter must be one the
        # receptor does not use.
        used = {L[21] for L in open(rec) if L[:4] in ("ATOM", "HETA")}
        lig = next((c for c in "ZYXWV" if c not in used), "Z")
        # TER records are kept, and are why the cartoon renders as a
        # protein: prepared.pdb carries one at every unmodelled loop,
        # and filtering to ATOM/HETATM alone discards them.
        out = [L for L in open(rec)
               if L[:4] in ("ATOM", "HETA") or L[:3] == "TER"]
        out.append("TER\n")
        for L in open(pose):
            if L[:4] in ("ATOM", "HETA"):
                out.append(L[:21] + lig + L[22:])
        out.append("TER\nEND\n")
        return "".join(out), None

    def pose_list(self, rundir, peptide):
        """Which ranked poses exist, with their scores, for the pose picker."""
        import glob
        d = os.path.join(rundir, "docking", jobmod.key(peptide))
        if not os.path.isdir(d):
            return {"poses": []}
        from . import adcp11
        modes = {m["rank"]: m for m in adcp11.read_modes(d)}
        out = []
        for f in sorted(glob.glob(os.path.join(d, "out_ranked_*.pdb")),
                        key=lambda p: int(p.rsplit("_", 1)[1][:-4])):
            r = int(f.rsplit("_", 1)[1][:-4])
            # read_modes names the column "affinity"; the table
            # and the pose picker both call it the docking score.
            out.append({"rank": r,
                        "score": (modes.get(r) or {}).get("affinity")})
        mm = None
        rj = os.path.join(rundir, "results", "%s.json" % jobmod.key(peptide))
        if os.path.isfile(rj):
            try:
                mm = json.load(open(rj)).get("mmgbsa")
            except ValueError:
                mm = None
        if mm:
            for i, row in enumerate(out):
                row["mmgbsa"] = mm[i] if i < len(mm) else None
        return {"poses": out}

    # How a SLURM state maps to the one word and colour a person needs.
    # Anything unrecognised is "queued" rather than "done": calling an
    # unknown state finished stops the page watching something that is
    # still running.
    STATE_DOT = {
        "RUNNING": "running", "COMPLETING": "running",
        "PENDING": "queued", "CONFIGURING": "queued", "REQUEUED": "queued",
        "RESIZING": "queued", "SUSPENDED": "queued",
        "COMPLETED": "done",
        "FAILED": "failed", "TIMEOUT": "failed", "NODE_FAIL": "failed",
        "OUT_OF_MEMORY": "failed", "BOOT_FAIL": "failed",
        "CANCELLED": "gone", "DEADLINE": "gone", "PREEMPTED": "gone",
        "REVOKED": "gone", "SPECIAL_EXIT": "gone",
    }

    @classmethod
    def summarise_state(cls, counts):
        """{SLURM state: n} -> one word for a dot.

        Worst first, with RUNNING ahead of everything except a failure: an array
        where one task failed and thirty continue is a run that needs attention.
        """
        if not counts:
            return "gone"
        kinds = set()
        for raw, n in counts.items():
            if not n:
                continue
            kinds.add(cls.STATE_DOT.get(str(raw).split()[0].upper(), "queued"))
        for want in ("failed", "running", "queued", "done", "gone"):
            if want in kinds:
                return want
        return "queued"

    def monitor(self):
        """Every run, with one word for its state.

        One call rather than one per run: a page polling twenty runs separately
        makes twenty squeue calls every few seconds.
        """
        out = []
        for row in self.list_runs()["runs"]:
            j = {}
            try:
                j = jobmod.load(row["rundir"])
            except Exception:                                   # noqa: BLE001
                pass
            stages, words = [], []
            for k, v in (j.get("jobs") or {}).items():
                counts = slurm.state(v, rundir=row["rundir"])
                word = self.summarise_state(counts)
                words.append(word)
                # Why it is waiting, not only that it is:
                # "PartitionTimeLimit" means the job asks for
                # longer than the partition allows and will
                # never start, while "Priority" means it is
                # in a queue. The two are otherwise
                # indistinguishable.
                why = (slurm.why_pending(v, rundir=row["rundir"])
                       if word == "queued" else "")
                stages.append({"stage": k, "id": v, "dot": word,
                               "counts": counts, "why": why,
                               "stuck": why in STUCK_REASONS})
            dot = "gone"
            for want in ("failed", "running", "queued", "done"):
                if want in words:
                    dot = want
                    break
            # A run held on something that will never clear is
            # not "queued".
            if any(st.get("stuck") for st in stages):
                dot = "failed"
            if not stages:
                dot = "new"
            row = dict(row, dot=dot, stages=stages)
            out.append(row)
        live = any(r["dot"] in ("running", "queued") for r in out)
        return {"runs": out, "root": self.runs, "live": live}

    def run_status(self, rundir):
        if not rundir or not os.path.isfile(jobmod.path_of(rundir)):
            return {"error": "no such run"}
        j = jobmod.load(rundir)
        jobs = {}
        for k, v in (j.get("jobs") or {}).items():
            jobs[k] = {"id": v, "state": slurm.state(v, rundir=rundir)}
        results = None
        rj = os.path.join(rundir, "results.json")
        if os.path.isfile(rj):
            try:
                results = json.load(open(rj))
            except ValueError:
                results = None
        results = fill_names(results, j)
        results = fill_rings(results, j)
        results = fill_terms(results, rundir)
        results = fill_combined(results, rundir)
        done = len([f for f in os.listdir(os.path.join(rundir, "results"))
                    if f.endswith(".json")]) if os.path.isdir(
                        os.path.join(rundir, "results")) else 0
        try:
            todo = pipeline.outstanding(rundir, j)
        except Exception:                                       # noqa: BLE001
            todo = []
        from . import tidy as tidymod
        # A run docked on a Mac with ADCP 1.0 says so above its table.
        mac = (j.get("adcp_engine") == mac10_mod.ENGINE
               or any(r.get("engine") == mac10_mod.LABEL
                      for r in ((results or {}).get("rows") or [])))
        return {"rundir": rundir, "job": j, "jobs": jobs,
                "prep_warnings": prep_warnings(rundir, j),
                "engine_note": mac10_mod.NOTE if mac else "",
                "outstanding": [{"index": i, "name": n, "sequence": s2}
                                for i, n, s2 in todo],
                "done": done, "total": len(j.get("peptides") or []),
                "results": results,
                # Whether the run is tidied, and its tidy job: the page's
                # Tidy button reads these. Neither walks the run.
                "tidy": dict(tidymod.status(rundir),
                             job=tidymod.job_state(rundir)),
                "log_tail": _tail(os.path.join(rundir, "logs"))}

    def tidy_run(self, rundir, undo=False):
        """Tidy a finished run, or undo it, in a short job of its own."""
        bad = self.check_run(rundir)
        if bad:
            return {"error": bad}
        from . import tidy as tidymod
        try:
            return tidymod.submit(rundir, undo_=bool(undo))
        except (RuntimeError, OSError, ValueError) as exc:
            return {"error": "could not submit the %s: %s"
                             % ("undo" if undo else "tidy", exc)}

    @staticmethod
    def _tidying(rundir):
        """Why a job that reads the run's files must wait, or None: a tidy
        moves them while it runs."""
        from . import tidy as tidymod
        job = tidymod.job_state(rundir)
        if job.get("live"):
            return ("this run is being %s (job %s); try again once that has "
                    "finished" % ("put back" if job.get("undo") else "tidied",
                                  job.get("job_id")))
        return None


def _tail(logdir, n=40):
    if not os.path.isdir(logdir):
        return ""
    files = [os.path.join(logdir, f) for f in os.listdir(logdir)]
    if not files:
        return ""
    newest = max(files, key=os.path.getmtime)
    try:
        lines = open(newest, errors="replace").read().splitlines()
    except OSError:
        return ""
    return "%s\n%s" % (os.path.basename(newest), "\n".join(lines[-n:]))


_ADAPT_LOCK = threading.Lock()
#: One receptor written by ADapt at a time: a second press, or a submission
#: while the first is still being written, waits for it and reuses it.
_ADAPT_WRITE_LOCK = threading.Lock()

#: How long one request to ADapt may take. Reading a large structure runs
#: AGFR's own processing and is not quick.
ADAPT_PROXY_SECONDS = 600


class _NoRedirect(urllib.request.HTTPRedirectHandler):
    """A redirect from ADapt goes back to the browser as it is. ADapt sends
    /adapt?t=... to /adapt/?t=..., because its page's paths are relative and
    need the slash; followed here instead, the browser would keep the
    address without it."""

    def redirect_request(self, *a, **kw):
        return None


#: No proxy for a connection to this same machine. A cluster that sets
#: http_proxy for the outside world would otherwise send 127.0.0.1 there.
_ADAPT_OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}),
                                            _NoRedirect)


class Handler(BaseHTTPRequestHandler):
    server_version = "adcp_dock"

    def log_message(self, *a):                              # quiet by default
        pass

    # ---- plumbing --------------------------------------------------------
    def _send(self, obj, code=200, ctype="application/json"):
        body = (obj if isinstance(obj, bytes)
                else json.dumps(obj).encode() if ctype == "application/json"
                else obj.encode())
        self.send_response(code)
        self.send_header("Content-Type", ctype + "; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(body)

    def _auth(self, q):
        return (q.get("t", [""])[0] or
                self.headers.get("X-Token", "")) == self.server.app.token

    def _send_file(self, real, download):
        """A file from the runs folder, already held to it by run_file.

        Opened without following a link at its last step, and sent only if
        what was opened is a regular file, so a file swapped for a link
        after it was checked is not followed. To be shown, it goes as plain
        text whatever it holds, so a page in the runs folder can never run
        as part of this one; to be saved, as an attachment. Neither is
        allowed to load anything.
        """
        flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
        try:
            fd = os.open(real, flags)
        except OSError as exc:
            return self._send("could not open it: %s" % exc.strerror, 404,
                              "text/plain")
        with os.fdopen(fd, "rb") as fh:
            st = os.fstat(fh.fileno())
            if not stat.S_ISREG(st.st_mode):
                return self._send("that is not a file", 404, "text/plain")
            if not download and st.st_size > App.VIEW_MAX:
                return self._send("that file is too large to show here; "
                                  "download it instead", 413, "text/plain")
            self.send_response(200)
            if download:
                self.send_header("Content-Type", "application/octet-stream")
                self.send_header("Content-Disposition",
                                 _attachment(os.path.basename(real)))
            else:
                self.send_header("Content-Type", "text/plain; charset=utf-8")
            self.send_header("Content-Length", str(st.st_size))
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Content-Security-Policy",
                             "default-src 'none'; sandbox")
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            shutil.copyfileobj(fh, self.wfile, 1 << 20)

    def _adapt_proxy(self, method, u):
        """Pass a request under /adapt/ to ADapt, and its answer back.

        ADapt runs with --base-path /adapt and its page uses relative paths
        only, so the path is forwarded as it is, prefix and all, and nothing
        in the answer is rewritten (its README, "Embedding ADapt").

        This tool's own token is not asked for here: ADapt's pages carry
        ADapt's token, and ADapt checks it on everything but its static
        files, exactly as when it is reached directly.
        """
        runner = getattr(self.server.app, "_adapt_runner", None)
        if runner is None or not runner.running():
            return self._send("ADapt is not running. Go back to the docking page "
                              "and press Open in ADapt.", 503, "text/plain")
        base = urllib.parse.urlparse(runner.url)
        target = "%s://%s%s%s" % (base.scheme, base.netloc, u.path or "/",
                                  ("?" + u.query) if u.query else "")
        data = None
        if method == "POST":
            n = int(self.headers.get("Content-Length") or 0)
            if n > MAX_UPLOAD:
                return self._send({"error": "that file is too large"}, 413)
            data = self.rfile.read(n) if n > 0 else b""
        hdrs = {k: self.headers[k] for k in ("Content-Type", "Accept", "X-Token")
                if self.headers.get(k)}
        req = urllib.request.Request(target, data=data, method=method, headers=hdrs)
        try:
            r = _ADAPT_OPENER.open(req, timeout=ADAPT_PROXY_SECONDS)
            code, body, rh = r.status, r.read(), r.headers
        except urllib.error.HTTPError as exc:
            code, body, rh = exc.code, exc.read(), exc.headers
        except (urllib.error.URLError, OSError) as exc:
            return self._send("ADapt did not answer (%s). Press Open in ADapt again."
                              % getattr(exc, "reason", exc), 502, "text/plain")
        ctype = rh.get("Content-Type", "application/octet-stream")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        for k in ("Content-Disposition", "Cache-Control", "Last-Modified",
                  "Location"):
            if rh.get(k):
                self.send_header(k, rh.get(k))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(body)

    def _body(self, limit=MAX_BODY):
        n = int(self.headers.get("Content-Length") or 0)
        if n <= 0 or n > limit:
            return {}
        try:
            return json.loads(self.rfile.read(n).decode())
        except ValueError:
            return {}

    def do_GET(self):
        u = urllib.parse.urlparse(self.path)
        q = urllib.parse.parse_qs(u.query)
        app = self.server.app
        if u.path == "/adapt" or u.path.startswith("/adapt/"):
            return self._adapt_proxy("GET", u)
        if u.path == "/docs" or u.path.startswith("/docs/"):
            return self._docs(u.path)
        if u.path.startswith("/static/"):
            return self._static(u.path)
        if u.path == "/":
            if not self._auth(q):
                return self._send("open the URL this tool printed; it "
                                  "carries the token", 403, "text/plain")
            return self._send(page.render(app.token), ctype="text/html")
        if not self._auth(q):
            return self._send({"error": "bad or missing token"}, 403)
        if u.path == "/cite":
            # The citation as a file, in three formats: .txt to
            # read, .bib for LaTeX, .ris for EndNote, Zotero and
            # Mendeley. Generated from cite.REFERENCES rather
            # than kept as files, which would drift. With ?run=
            # it is scoped to that run's job.json.
            fmt = (q.get("fmt") or ["txt"])[0]
            run = (q.get("run") or [""])[0]
            from . import cite as _c
            if fmt not in _c.FILENAMES:
                return self._send("unknown citation format %r" % fmt,
                                  400, "text/plain")
            j = None
            if run:
                bad = app.check_run(run)
                if bad:
                    return self._send(bad, 400, "text/plain")
                try:
                    from . import job as _j
                    j = _j.load(run)
                except Exception:                         # noqa: BLE001
                    j = None                              # fall back to all
            try:
                body = _c.formatted(fmt, j)
            except Exception as exc:                      # noqa: BLE001
                return self._send("could not build the citation: %s" % exc,
                                  500, "text/plain")
            if (q.get("dl") or [""])[0]:
                data = body.encode()
                self.send_response(200)
                self.send_header("Content-Type",
                                 "text/plain; charset=utf-8")
                self.send_header("Content-Disposition",
                                 'attachment; filename="%s"'
                                 % _c.FILENAMES[fmt])
                self.send_header("Content-Length", str(len(data)))
                self.end_headers()
                return self.wfile.write(data)
            return self._send(body, ctype="text/plain")
        if u.path == "/download":
            # A ZIP streamed back rather than written into the
            # run directory: a run that grows a file on every
            # download is its own problem, and deciding when to
            # delete it is another.
            run = (q.get("run") or [""])[0]
            what = (q.get("what") or ["results"])[0]
            bad = app.check_run(run)
            if bad:
                return self._send(bad, 400, "text/plain")
            try:
                from . import bundle as _b
                name, data = _b.zip_bytes(run, "all" if what == "all"
                                          else "results")
            except Exception as exc:                      # noqa: BLE001
                return self._send("could not build the bundle: %s" % exc,
                                  500, "text/plain")
            self.send_response(200)
            self.send_header("Content-Type", "application/zip")
            self.send_header("Content-Disposition",
                             'attachment; filename="%s"' % name)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            return self.wfile.write(data)
        if u.path == "/api/sizes":
            return self._send(app.run_sizes((q.get("run") or [""])[0]))
        if u.path == "/api/monitor":
            return self._send(app.monitor())
        if u.path == "/api/roots":
            return self._send(app.roots())
        if u.path == "/api/tools":
            return self._send(app.tools())
        if u.path == "/api/browse":
            if (q.get("dirs") or [""])[0] in ("1", "true", "yes"):
                return self._send(app.browse(
                    (q.get("path") or [""])[0], dirs_only=True))
            return self._send(app.browse(q.get("path", [""])[0],
                                         kind=(q.get("kind") or [""])[0]))
        if u.path == "/api/target":
            return self._send(app.target_file((q.get("path") or [""])[0]))
        if u.path == "/api/receptor":
            return self._send(app.receptor(q.get("path", [""])[0]))
        if u.path == "/api/fetch":
            return self._send(app.fetch_code(q.get("code", [""])[0]))
        if u.path == "/api/rawstructure":
            # The site picker needs the receptor as PDB before
            # any run exists.
            path = q.get("path", [""])[0]
            if (not path or not os.path.isfile(path)
                    or not structure.is_structure(path)):
                return self._send("not a structure file", 404, "text/plain")
            try:
                # One conformer per residue, the one the surface is
                # computed from (structure.one_conformer).
                atoms = list(structure.one_conformer(structure.read(path)))
            except Exception as exc:                            # noqa: BLE001
                return self._send("could not read it: %s" % exc, 400,
                                  "text/plain")
            # Written back out as PDB whatever came in, because
            # the viewer reads PDB and an mmCIF would render
            # nothing.
            import tempfile
            with tempfile.NamedTemporaryFile("w+", suffix=".pdb",
                                             delete=False) as fh:
                tmp = fh.name
            structure.write_pdb(atoms, tmp)
            with open(tmp) as fh:
                text = fh.read()
            os.unlink(tmp)
            return self._send(text, ctype="text/plain")
        if u.path == "/api/sequence":
            return self._send(app.sequence(q.get("path", [""])[0],
                                           q.get("chain", [""])[0]))
        if u.path == "/api/modifications":
            return self._send(app.modifications())
        if u.path == "/api/depict":
            # One residue drawn for the Modifications dialog: res=<letter>
            # or code=<ADCP code>, d=1 for the D form, dark=1 for a dark page.
            one = lambda k: (q.get(k) or [""])[0]            # noqa: E731
            try:
                from . import depict
                svg = depict.svg(one("res").upper(), one("code").upper(),
                                 one("d") == "1", one("dark") == "1")
            except (KeyError, ValueError) as e:
                return self._send({"error": "no drawing for that residue: %s" % e}, 404)
            except ImportError:
                return self._send({"error": "drawing needs RDKit"}, 501)
            return self._send(svg, ctype="image/svg+xml")
        if u.path == "/api/runs":
            return self._send(app.list_runs())
        if u.path == "/api/folder":
            # The Folders view: one folder of the runs folder, and nothing
            # outside it.
            return self._send(app.runs_folder((q.get("path") or [""])[0]))
        if u.path == "/api/runfile":
            # A file from the runs folder, to show in the page or, with
            # dl=1, to save. Held to the runs folder the same way.
            real, why, code = app.run_file((q.get("path") or [""])[0])
            if real is None:
                return self._send(why, code, "text/plain")
            return self._send_file(real, (q.get("dl") or [""])[0]
                                   in ("1", "true", "yes"))
        if u.path == "/api/run":
            return self._send(app.run_status(q.get("dir", [""])[0]))
        if u.path == "/api/residues":
            return self._send(app.residues(q.get("dir", [""])[0],
                                           q.get("peptide", [""])[0]))
        if u.path == "/api/poses":
            return self._send(app.pose_list(q.get("dir", [""])[0],
                                            q.get("peptide", [""])[0]))
        if u.path == "/api/surface":
            return self._send(app.surface(q))
        if u.path == "/api/adapt/status":
            return self._send(app.adapt_status())
        if u.path == "/api/adapt/state":
            return self._send(app.adapt_state((q.get("session") or [""])[0],
                                              (q.get("folder") or [""])[0]))
        if u.path == "/api/adapt/saved":
            one = lambda k: (q.get(k) or [""])[0]            # noqa: E731
            return self._send(app.adapt_saved(one("session"), one("folder"),
                                              one("save") or None))
        if u.path == "/api/structure":
            text, err = app.structure_text(
                q.get("dir", [""])[0], q.get("peptide", [""])[0] or None,
                q.get("rank", ["1"])[0], q.get("what", ["complex"])[0])
            if err:
                return self._send({"error": err}, 404)
            return self._send(text, ctype="text/plain")
        return self._send({"error": "no such endpoint"}, 404)

    def _docs(self, path):
        """The documentation: its pages, their stylesheet and their pictures.

        No token: it is the same help for everyone, it holds nothing about
        anybody's runs, and a link to it should just open. What to send is
        decided by the shared documentation module, which the three
        interfaces answer the same way; this writes it out.
        """
        from . import docs
        code, ctype, body, location = docs.respond(path)
        self.send_response(code)
        if location:
            self.send_header("Location", location)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(body)

    def _static(self, path):
        """The bundled JavaScript assets, served from this package so the viewer works without network access.
        """
        name = os.path.basename(path)
        here = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static")
        full = os.path.join(here, name)
        if os.path.dirname(os.path.abspath(full)) != os.path.abspath(here) \
                or not os.path.isfile(full):
            return self._send({"error": "no such file"}, 404)
        ctype = ("application/javascript" if name.endswith(".js")
                 else "text/css" if name.endswith(".css") else "text/plain")
        with open(full, "rb") as fh:
            body = fh.read()
        self.send_response(200)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "max-age=86400")
        self.end_headers()
        self.wfile.write(body)

    def do_POST(self):
        u = urllib.parse.urlparse(self.path)
        q = urllib.parse.parse_qs(u.query)
        if u.path == "/adapt" or u.path.startswith("/adapt/"):
            return self._adapt_proxy("POST", u)
        if not self._auth(q):
            return self._send({"error": "bad or missing token"}, 403)
        app = self.server.app
        body = self._body(MAX_UPLOAD if u.path == "/api/upload"
                          else MAX_BODY)
        if u.path == "/api/adapt/warm":
            return self._send(app.adapt_warm())
        if u.path == "/api/site":
            return self._send(app.site(body))
        if u.path == "/api/pockets":
            return self._send(app.pockets(body))
        if u.path == "/api/peptides":
            return self._send(app.peptides(body))
        if u.path == "/api/upload":
            return self._send(app.upload(body))
        if u.path == "/api/paste":
            return self._send(app.paste(body))
        if u.path == "/api/create":
            return self._send(app.create(body))
        if u.path == "/api/estimate":
            return self._send(app.estimate_run(body))
        if u.path == "/api/delete":
            return self._send(app.delete_run(body.get("rundir")))
        if u.path == "/api/open_folder":
            return self._send(app.open_folder(body.get("rundir")))
        if u.path == "/api/mkdir":
            return self._send(app.make_dir(body.get("parent"),
                                           body.get("name")))
        if u.path == "/api/runsroot":
            return self._send(app.set_runs(body.get("path")))
        if u.path == "/api/report":
            return self._send(app.report_problem(body))
        if u.path == "/api/resume":
            return self._send(app.resume_run(body.get("rundir")))
        if u.path == "/api/panel":
            return self._send(app.panel_run(body.get("rundir")))
        if u.path == "/api/rerank":
            return self._send(app.rerank_run(body.get("rundir")))
        if u.path == "/api/rescore":
            return self._send(app.rescore_run(body.get("rundir"),
                                              body.get("engine")))
        if u.path == "/api/tidy":
            return self._send(app.tidy_run(body.get("rundir"),
                                           body.get("undo")))
        if u.path == "/api/adapt":
            return self._send(app.adapt_open(body))
        if u.path == "/api/adapt/write":
            return self._send(app.adapt_write(body))
        return self._send({"error": "no such endpoint"}, 404)


def build(host="127.0.0.1", port=8765, runs="~/adcp_runs", tools_roots=()):
    httpd = ThreadingHTTPServer((host, port), Handler)
    httpd.app = App(runs, tools_roots)
    httpd.token = httpd.app.token
    return httpd


def serve(host="127.0.0.1", port=8765, runs="~/adcp_runs", tools_roots=()):
    httpd = build(host, port, runs, tools_roots)
    # The port the socket actually received, which is not `port` when 0
    # was requested.
    bound = httpd.server_address[1]
    url = "http://%s:%d/?t=%s" % (host, bound, httpd.token)
    from . import credit
    print(credit())
    from . import APP as _app
    print("\n%s is at\n\n    %s\n" % (_app, url))
    print("runs are kept in %s" % httpd.app.runs)
    miss = httpd.app.tools()["missing"]
    if miss:
        print("\nNOT FOUND: %s\nSet ADCPSUITE_HOME and AMBERHOME, or pass "
              "--tools-root; ./install.sh installs the ADCP 1.1 suite."
              % ", ".join(miss))
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    finally:
        httpd.app.adapt_stop()
