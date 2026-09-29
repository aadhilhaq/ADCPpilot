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
import threading
import urllib.parse
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

from . import fetch, job as jobmod, page, peptide, pipeline
from . import sizing
from . import slurm, structure
from . import tools as toolmod

MAX_BODY = 4 << 20
# Uploads have their own ceiling. Every other request is a small JSON object,
# for which 4 MB is generous; a receptor with waters and several models can
# run to tens of megabytes.
MAX_UPLOAD = 64 << 20


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
                     "receptor pdbqt"):
            if s.startswith(name + ":"):
                steps[name] = s[len(name) + 1:].strip()
    out = []
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


def _open_detached(cmd):
    """Start a file manager and return at once; its window is the user's."""
    import subprocess
    subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                     stdin=subprocess.DEVNULL, start_new_session=True)


# The function that starts the file manager: a name, so the tests can put
# a recorder in its place instead of opening windows.
OPENER = _open_detached


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
                "roots": self.tools_roots,
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
                # ADCP 1.1, when it is installed: the page offers it, its
                # non-standard residue syntax and OpenMM minimization.
                "adcp11": self._adcp11()[0],
                "adcp11_version": self._adcp11()[1]}

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
                # The licence, not just the program: Prime needs 8 of the 25
                # rescoring tokens and the modern force field has one seat,
                # so what is free decides what a run would actually get.
                try:
                    row["licence"] = gbengines.SCHRODINGER.plan(timeout=3)
                except Exception:                              # noqa: BLE001
                    row["licence"] = {"readable": False}
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
        """The keyword arguments a peptide check needs for the chosen engine."""
        mx = int(body.get("max_length") or peptide.MAX_LENGTH)
        engine = str(body.get("adcp_engine") or "auto")
        exe = self._adcp11()[0] if engine != "1.0" else None
        return {"max_length": mx, "extended": bool(exe),
                "library": self._a11lib if exe else None}

    def roots(self):
        """Where the browser should offer to start.

        Shared group space first, then the user's own: a shared installation means
        the structures and runs people care about are usually in the group
        allocation.
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
            for base in ("/scratch/group", "/group", "/work/group"):
                add(os.path.join(base, g), "%s (group)" % g)
        add(os.environ.get("SCRATCH"), "your scratch")
        add(os.path.join("/scratch/user", user) if user else "", "your scratch")
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

    def browse(self, path, dirs_only=False):
        """A directory listing. Structures first, then directories: the
        person is looking for a receptor, not for a filesystem."""
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
                    dirs.append({"name": n, "path": full})
                elif structure.is_structure(n):
                    # A broken symlink must not empty the
                    # listing: getsize raises
                    # FileNotFoundError for a link with
                    # no target, so that error is caught
                    # for each entry.
                    try:
                        size = os.path.getsize(full)
                    except OSError:
                        continue
                    files.append({"name": n, "path": full, "size": size})
        except PermissionError:
            return {"error": "no permission to read %s" % path}
        return {"path": path, "parent": os.path.dirname(path),
                "dirs": dirs, "files": [] if dirs_only else files,
                "writable": os.access(path, os.W_OK | os.X_OK)}

    def receptor(self, path):
        """Chains and residues, so the page can offer them."""
        if not path or not os.path.isfile(path):
            return {"error": "no such file: %s" % path}
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
                             for c in ch if c["kind"] == "protein"}}

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
        """(heavy atoms, distinct non-water heteroatom kinds).

        Counted over the chains that will be kept rather than the whole file: a
        structure is routinely several copies of a complex, and costing a run by
        the file means costing work it will never do.
        """
        want = {c.strip() for c in (chains or []) if c and c.strip()}
        atoms, kinds = 0, set()
        try:
            with open(path, errors="replace") as fh:
                for line in fh:
                    if line.startswith("ENDMDL"):
                        break
                    if line[:4] not in ("ATOM", "HETA"):
                        continue
                    if line[16] not in (" ", "A"):     # one altloc only
                        continue
                    chain = line[21].strip() or "_"
                    if want and chain not in want:
                        continue
                    res = line[17:20].strip()
                    if line.startswith("HETATM"):
                        if res in ("HOH", "WAT", "DOD", "D2O"):
                            continue
                        kinds.add(res)
                    if (line[76:78].strip() or res[:1]).upper() == "H":
                        continue
                    atoms += 1
        except OSError:
            return 0, 0
        return atoms, len(kinds)

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
            fill_loops=bool(body.get("fill_loops")), box=box)

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
                            "on the left) and go to this folder. Step 3 of "
                            "the notebook downloads a run as a zip."}
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
        try:
            return {"jobs": pipeline.submit_panel(rundir, log=lambda *a: None)}
        except SystemExit as exc:
            return {"error": str(exc)}
        except Exception as exc:                                # noqa: BLE001
            return {"error": "could not submit the measurement: %s" % exc}

    def resume_run(self, rundir):
        if not rundir or not os.path.isfile(jobmod.path_of(rundir)):
            return {"error": "no such run"}
        try:
            return pipeline.resume(rundir, log=lambda *a: None)
        except SystemExit as exc:
            return {"error": str(exc)}
        except Exception as exc:                                # noqa: BLE001
            return {"error": "could not resubmit: %s" % exc}

    def site(self, body):
        """Where the box would go, before anything is submitted."""
        path = body.get("receptor")
        if not path or not os.path.isfile(path):
            return {"error": "choose a receptor first"}
        all_atoms = structure.read(path)
        keep = body.get("chains") or []
        atoms = [a for a in all_atoms if a.chain in keep] if keep else all_atoms
        hs = structure.parse_hotspots(body.get("hotspots"),
                                      (keep or [None])[0])
        if not hs:
            return {"error": "no residues chosen"}
        centre, size, note = structure.box_for(
            atoms, hs, padding=float(body.get("box_padding", 8.0)),
            minimum=float(body.get("box_min", 20.0)),
            maximum=float(body.get("box_max", 30.0)))
        if centre is None:
            return {"error": note}
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
            found = _rec.chains_in_box(path, centre, size, exclude=loaded,
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
        return {"centre": [round(c, 3) for c in centre],
                "size": list(size), "note": note, "n": len(hs),
                "other_chains": others}

    def peptides(self, body):
        """Live validation. Positions come back so the box can underline."""
        kw = self._pepcheck(body)
        out = []
        for row in body.get("peptides") or []:
            seq, probs = peptide.check(row.get("sequence", ""), **kw)
            out.append({"name": row.get("name", ""), "sequence": seq,
                        "length": len(peptide.residues(seq)),
                        "ok": not probs,
                        "problems": [p.as_dict() for p in probs]})
        return {"peptides": out,
                "ready": bool(out) and all(p["ok"] for p in out)}

    def paste(self, body):
        """A pasted block or FASTA, split into rows the page can show."""
        rows = peptide.parse_many(body.get("text", ""), **self._pepcheck(body))
        return {"peptides": [{"name": n, "sequence": s,
                              "ok": not p,
                              "problems": [x.as_dict() for x in p]}
                             for n, s, p in rows]}

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
        if not (structure.is_structure(safe) or safe.lower().endswith(lists)):
            return {"error": "%s is not a structure or a peptide list. "
                             "Structures: .pdb, .cif, .ent, .pdbqt; lists: "
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
            problems = jobmod.validate(j)
            if j.get("mmgbsa") and not toolmod.discover(
                    roots=self.tools_roots).can_rescore():
                problems.append("MM-GBSA needs AmberTools (tleap, sander, "
                                "MMPBSA.py), which was not found. Untick "
                                "MM-GBSA, or install AmberTools.")
            if problems:
                return {"error": "this run is not ready", "problems": problems}
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
        results = fill_terms(results, rundir)
        done = len([f for f in os.listdir(os.path.join(rundir, "results"))
                    if f.endswith(".json")]) if os.path.isdir(
                        os.path.join(rundir, "results")) else 0
        try:
            todo = pipeline.outstanding(rundir, j)
        except Exception:                                       # noqa: BLE001
            todo = []
        return {"rundir": rundir, "job": j, "jobs": jobs,
                "prep_warnings": prep_warnings(rundir, j),
                "outstanding": [{"index": i, "name": n, "sequence": s2}
                                for i, n, s2 in todo],
                "done": done, "total": len(j.get("peptides") or []),
                "results": results,
                "log_tail": _tail(os.path.join(rundir, "logs"))}


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
            return self._send(app.browse(q.get("path", [""])[0]))
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
                atoms = structure.read(path)
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
        if u.path == "/api/run":
            return self._send(app.run_status(q.get("dir", [""])[0]))
        if u.path == "/api/residues":
            return self._send(app.residues(q.get("dir", [""])[0],
                                           q.get("peptide", [""])[0]))
        if u.path == "/api/poses":
            return self._send(app.pose_list(q.get("dir", [""])[0],
                                            q.get("peptide", [""])[0]))
        if u.path == "/api/structure":
            text, err = app.structure_text(
                q.get("dir", [""])[0], q.get("peptide", [""])[0] or None,
                q.get("rank", ["1"])[0], q.get("what", ["complex"])[0])
            if err:
                return self._send({"error": err}, 404)
            return self._send(text, ctype="text/plain")
        return self._send({"error": "no such endpoint"}, 404)

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
        if not self._auth(q):
            return self._send({"error": "bad or missing token"}, 403)
        app = self.server.app
        body = self._body(MAX_UPLOAD if u.path == "/api/upload"
                          else MAX_BODY)
        if u.path == "/api/site":
            return self._send(app.site(body))
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
        print("\nNOT FOUND: %s\nSet ADFRSUITE_HOME and AMBERHOME, or pass "
              "--tools-root." % ", ".join(miss))
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
