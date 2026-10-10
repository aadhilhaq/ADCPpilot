"""ADapt, a separate tool for preparing a receptor before docking.

By Aadhil Haq

ADapt, by Michel Sanner, shows what a PDB file holds: each biological
assembly, alternate locations, side chains with missing atoms, missing
residues, ligands, cofactors and metals. It writes receptor and ligand PDBQT
files the way AGFR does, and beside each receptor PDBQT a PDB of the same
receptor with every hydrogen reduce placed (<id>_rec_H.pdb).

It drives AGFR's own code, so it runs on the ADCP 1.1 environment's Python and
is started here as its own program with that environment named. Nothing of it
is imported or shipped with this package; it is found where it was installed.
install.sh and update.sh fetch it from

    https://github.com/michelsanner/adapt

into a folder named adapt beside ADFRsuite, anonymously, once that repository
is public (scripts/adapt_fetch.sh). It is also found at ADAPT_HOME, or beside
this checkout.

It is embedded the way its README's "Embedding ADapt" section describes: it
is started with --base-path /adapt, so this interface's server forwards
/adapt/... to it unchanged; a receptor is handed to it with its own
api/upload, and the page is opened on that session with session= in the URL.
What it saves is read back from the session folder it wrote: the receptor
PDBQT, the _rec_H.pdb and session.json.
"""
import hashlib
import json
import os
import re
import signal
import subprocess
import threading
import time
import urllib.error
import urllib.parse
import urllib.request

REPO = "https://github.com/michelsanner/adapt"
# The port ADapt listens on. Zero means the operating system picks a free one,
# which is right on a shared login node. Anyone reaching this interface through
# an ssh tunnel wants the same port every time instead, since the tunnel names
# it: ADCP_ADAPT_PORT=40805 pins it.
PORT_ENV = "ADCP_ADAPT_PORT"

# The line ADapt prints once it is serving: its address with its token.
_URL = re.compile(r"(https?://\S+/\?t=\S+)")

#: The path ADapt is served under on this interface's own server. ADapt is
#: started with --base-path set to it, so a request for /adapt/api/... is
#: forwarded as it is and ADapt takes the prefix off itself.
BASE_PATH = "/adapt"

#: The commit of ADapt that added what embedding needs: --base-path, the
#: session= and pdb= URL parameters, api/session and _rec_H.pdb.
EMBED_COMMIT = "6dae90b"

#: What ADapt reads: PDB files, plain or gzipped.
READS = (".pdb", ".ent", ".pdb.gz", ".ent.gz")

#: No proxy for a connection to this same machine.
_OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))


def candidates(roots=()):
    """Where an ADapt checkout is looked for, in order."""
    out = []
    env = os.environ.get("ADAPT_HOME")
    if env:
        out.append(env)
    adfr = os.environ.get("ADFRSUITE_HOME")
    inst = os.environ.get("ADCP_HOME")
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for root in (list(roots)
                 + ([os.path.dirname(os.path.abspath(adfr))] if adfr else [])
                 + ([os.path.join(inst, "tools", "docking")] if inst else [])
                 + [os.path.dirname(here)]):
        out += [os.path.join(root, "adapt"), os.path.join(root, "ADapt")]
    out.append(os.path.expanduser("~/.adcp-dock/adapt"))
    return out


def is_adapt(folder):
    return (os.path.isfile(os.path.join(folder, "bin", "adapt"))
            and os.path.isfile(os.path.join(folder, "adapt", "__main__.py")))


def find(roots=()):
    """The ADapt checkout, or None."""
    for d in candidates(roots):
        if d and is_adapt(d):
            return os.path.abspath(d)
    return None


def embeddable(home):
    """Whether this checkout can be embedded: it takes --base-path, which
    came with the session= URL parameter, api/session and _rec_H.pdb."""
    try:
        text = open(os.path.join(home, "adapt", "__main__.py")).read()
    except OSError:
        return False
    return "--base-path" in text


def commit(home):
    """The checkout's commit, short, read from .git without running git;
    '' when it cannot be read."""
    git = os.path.join(home or "", ".git")
    try:
        ref = open(os.path.join(git, "HEAD")).read().strip()
        if ref.startswith("ref: "):
            name = ref[5:]
            path = os.path.join(git, name)
            if os.path.isfile(path):
                return open(path).read().strip()[:7]
            for line in open(os.path.join(git, "packed-refs")):
                if line.strip().endswith(" " + name):
                    return line.split()[0][:7]
            return ""
        return ref[:7]
    except OSError:
        return ""


# --- talking to a running ADapt directly -----------------------------------

def _base(url):
    """scheme://host:port/adapt from the address ADapt printed."""
    u = urllib.parse.urlparse(url)
    path = u.path.rstrip("/")
    if path.endswith(BASE_PATH):
        path = path[:-len(BASE_PATH)]
    return "%s://%s%s%s" % (u.scheme, u.netloc, path, BASE_PATH)


def token_of(url):
    q = urllib.parse.parse_qs(urllib.parse.urlparse(url or "").query)
    return (q.get("t") or [""])[0]


def call(url, path, params=None, data=None, timeout=600):
    """One request to ADapt's own API, made by this server, as a dict.

    `url` is the address ADapt printed; `path` is "api/upload" and so on.
    Raises RuntimeError with ADapt's own message when it refuses.
    """
    q = dict(params or {})
    q["t"] = token_of(url)
    full = "%s/%s?%s" % (_base(url), path.lstrip("/"), urllib.parse.urlencode(q))
    req = urllib.request.Request(
        full, data=data, method="POST" if data is not None else "GET",
        headers={"Content-Type": "application/octet-stream"} if data is not None else {})
    try:
        with _OPENER.open(req, timeout=timeout) as r:
            body = r.read()
    except urllib.error.HTTPError as exc:
        body = exc.read()
        try:
            why = json.loads(body.decode()).get("error")
        except (ValueError, AttributeError):
            why = None
        raise RuntimeError(why or "ADapt answered %s" % exc.code)
    except (urllib.error.URLError, OSError) as exc:
        raise RuntimeError("ADapt did not answer (%s)" % getattr(exc, "reason", exc))
    try:
        return json.loads(body.decode())
    except ValueError:
        raise RuntimeError("ADapt's answer was not JSON")


def answering(url, timeout=3):
    """Whether a running ADapt answers a request: its list of sessions, the
    cheapest thing it serves. Any answer counts, a refusal included; only no
    answer at all does not."""
    if not url:
        return False
    try:
        call(url, "api/sessions", timeout=timeout)
        return True
    except RuntimeError as exc:
        return not str(exc).startswith("ADapt did not answer")


def upload(url, path, name=None, timeout=600):
    """Hand a structure file to a running ADapt (its api/upload): a new
    session, as a dict with its id ("session") and folder ("folder")."""
    name = name or os.path.basename(path)
    with open(path, "rb") as fh:
        data = fh.read()
    return call(url, "api/upload", {"name": name}, data=data, timeout=timeout)


def page_url(url, session=None, pdb=None, tab=None):
    """The address of ADapt's page under this server, with its token and
    what it should open."""
    q = [("t", token_of(url))]
    if session:
        q.append(("session", session))
    if pdb:
        q.append(("pdb", pdb))
    if tab:
        q.append(("tab", tab))
    return BASE_PATH + "/?" + urllib.parse.urlencode(q)


# --- what ADapt saved --------------------------------------------------------

class NoSave(ValueError):
    """A session folder holds no receptor ADapt saved, or not the one asked
    for, or its files are not what session.json says they are."""


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _header(path):
    """The REMARK 850 lines ADapt and AGFR write at the top of a receptor."""
    out = []
    try:
        for line in open(path, errors="replace"):
            if line.startswith("REMARK 850"):
                out.append(line.rstrip())
            elif line.startswith(("ATOM", "HETATM")):
                break
    except OSError:
        pass
    return out


WATERS = {"HOH", "WAT", "DOD", "H2O"}


def _hetero(path):
    """The hetero groups in a receptor file, as 'ZN A:401', waters left out."""
    seen, out = set(), []
    for line in open(path, errors="replace"):
        if not line.startswith("HETATM"):
            continue
        name = line[17:20].strip()
        if name in WATERS:
            continue
        key = "%s %s:%s%s" % (name, line[21].strip() or "_", line[22:26].strip(),
                               line[26].strip())
        if key not in seen:
            seen.add(key)
            out.append(key)
    return out


def read_save(folder, save=None):
    """The receptor ADapt saved in a session folder, from its session.json.

    `save` is a save number; without one, the newest receptor save. Returns a
    dict: the folder, the save number and its tab, the PDBQT and _rec_H.pdb
    with their md5, the biomolecule line, the hetero groups the receptor
    keeps, the ion charges, alternate locations and mutations chosen, and
    reduce's version and histidine decisions. Raises NoSave saying why.
    """
    folder = os.path.abspath(os.path.expanduser(str(folder or "").strip()))
    sj = os.path.join(folder, "session.json")
    try:
        rec = json.load(open(sj))
    except OSError:
        raise NoSave("%s holds no session.json, so it is not an ADapt session "
                     "folder" % folder)
    except ValueError:
        raise NoSave("%s could not be read" % sj)
    saves = [s for s in rec.get("saves") or []
             if isinstance(s, dict) and s.get("kind") == "receptor"
             and s.get("files")]
    if not saves:
        raise NoSave("ADapt has not saved a receptor in this session yet: "
                     "press receptor PDBQT in ADapt")
    if save not in (None, "", 0, "0"):
        try:
            want = int(save)
        except (TypeError, ValueError):
            raise NoSave("%r is not a save number" % (save,))
        pick = next((s for s in saves if int(s.get("n", 0)) == want), None)
        if pick is None:
            raise NoSave("this session has no receptor save %d; it has %s"
                         % (want, ", ".join(str(s.get("n")) for s in saves)))
    else:
        pick = max(saves, key=lambda s: int(s.get("n", 0)))
    sub = os.path.join(folder, pick.get("folder") or "")
    files = {f.get("name"): f.get("md5") for f in pick.get("files") or []}
    pdbqt = next((n for n in files if n.endswith("_rec.pdbqt")), None)
    pdbh = next((n for n in files if n.endswith("_rec_H.pdb")), None)
    if not pdbqt:
        raise NoSave("save %s has no receptor PDBQT" % pick.get("n"))
    if not pdbh:
        raise NoSave("save %s has no _rec_H.pdb: this ADapt is older than "
                     "the version that writes it (%s); update ADapt and save "
                     "the receptor again" % (pick.get("n"), EMBED_COMMIT))
    paths = {}
    for n in (pdbqt, pdbh):
        p = os.path.join(sub, n)
        if not os.path.isfile(p):
            raise NoSave("%s is missing from the session folder" % p)
        got = _md5(p)
        if files.get(n) and got != files[n]:
            raise NoSave("%s has changed since ADapt wrote it" % p)
        paths[n] = (p, got)
    red = pick.get("reduce") or {}
    his = {}
    for d in red.get("decisions") or []:
        if str(d.get("residue", "")).split(":")[1:2] == ["HIS"]:
            his[d.get("residue")] = d.get("decision")
    biomol = next((l.split("biomolecule", 1)[1].strip()
                   for l in _header(paths[pdbqt][0])
                   if "AGFR biomolecule" in l), "")
    env = rec.get("environment") or {}
    tool = rec.get("tool") or {}
    return {
        "folder": folder, "save": int(pick.get("n", 0)),
        "saves": [int(s.get("n", 0)) for s in saves],
        "tab": pick.get("tab") or "", "at": pick.get("at") or "",
        "pdbid": rec.get("pdbid") or "",
        "input": (rec.get("input") or {}).get("name") or "",
        "pdbqt": paths[pdbqt][0], "pdbqt_md5": paths[pdbqt][1],
        "pdb": paths[pdbh][0], "pdb_md5": paths[pdbh][1],
        "session_json": sj,
        "biomolecule": biomol,
        "hetero": _hetero(paths[pdbh][0]),
        "charges": pick.get("charges") or {},
        "altlocs": pick.get("altlocs") or {},
        "mutations": pick.get("mutations") or {},
        "reduce": red.get("version") or "",
        "his": his,
        "adapt_version": "%s %s" % (tool.get("name") or "ADapt",
                                    tool.get("version") or ""),
        "adapt_commit": (env.get("adapt_commit") or "")[:7],
        "adfr": env.get("adfr_dist") or "",
    }


# --- what has been changed in a session, and writing its receptor -------------

def _describe(ev):
    """One edit from session.json's events, in a few plain words."""
    what = ev.get("what")
    node = str(ev.get("node") or ev.get("key") or "an item")
    if len(node) <= 2 and node.isalnum():
        # ADapt names a chain's row by its letter alone.
        node = "chain " + node
    if what == "check":
        col = "receptor" if ev.get("column") == "rec" else "ligand"
        state = ev.get("state")
        verb = ("in the " if state == "checked" else
                "partly in the " if state == "partial" else "out of the ")
        return "%s %s%s" % (node, verb, col)
    if what == "tab":
        return "assembly %s" % (ev.get("name") or ev.get("tab") or "")
    if what == "altloc":
        return "%s at alternate location %s" % (node, ev.get("altloc"))
    if what == "charge":
        return "%s charge %s" % (node, ev.get("charge"))
    if what == "sidechain":
        return "%s side chain %s" % (
            node, "completed" if ev.get("action") == "complete" else "as it is")
    if what == "mutate":
        return "%s mutated to %s" % (node, str(ev.get("choice") or "").replace(
            "Base Res:", "").strip())
    return "%s %s" % (what or "edit", node)


def read_state(folder):
    """What has been done in an ADapt session, from the session.json ADapt
    rewrites after every edit and every save.

    ADapt records each edit in order (a box ticked or unticked, another
    assembly, an alternate location, a charge, a side chain, a mutation), so
    the number of edits says whether anything has been changed since the
    session was opened or since a save. Looking only at the file means this
    never touches the running session: ADapt's own api/session hands out the
    log lines its page has not shown yet, and asking it would take them from
    the page.

    Returns {folder, pdbid, input, tab, edits, changes, saves}; `changes` is
    the last few edits in words, `saves` each receptor save with the number
    of edits made before it. Raises NoSave when there is no session.json.
    """
    folder = os.path.abspath(os.path.expanduser(str(folder or "").strip()))
    sj = os.path.join(folder, "session.json")
    try:
        rec = json.load(open(sj))
    except OSError:
        raise NoSave("%s holds no session.json, so it is not an ADapt session "
                     "folder" % folder)
    except ValueError:
        # Caught between ADapt's write and its rename: the next look reads it.
        raise NoSave("%s could not be read" % sj)
    events = [e for e in rec.get("events") or [] if isinstance(e, dict)]
    tab = rec.get("opened_tab") or ""
    for e in events:
        if e.get("what") == "tab" and e.get("name"):
            tab = e["name"]
    saves = []
    for s in rec.get("saves") or []:
        if not isinstance(s, dict) or s.get("kind") != "receptor" \
                or not s.get("files"):
            continue
        at = str(s.get("at") or "")
        # Seconds are all ADapt records, so an edit in the same second as a
        # save counts as made before it.
        before = sum(1 for e in events if str(e.get("at") or "") <= at)
        saves.append({"n": int(s.get("n", 0)), "at": at, "tab": s.get("tab") or "",
                      "edits_before": before})
    return {"folder": folder, "pdbid": rec.get("pdbid") or "",
            "input": (rec.get("input") or {}).get("name") or "",
            "tab": tab, "edits": len(events),
            "changes": [_describe(e) for e in events[-6:]],
            "saves": saves}


def session_folder(url, session):
    """The folder of a session ADapt holds, from its list of sessions (its
    api/sessions, which changes nothing). '' when it does not hold it."""
    try:
        got = call(url, "api/sessions", timeout=30)
    except RuntimeError:
        return ""
    for s in got.get("sessions") or []:
        if s.get("session") == session:
            return os.path.expanduser(s.get("folder_display") or "")
    return ""


def write_receptor(url, session, timeout=600):
    """Ask a running ADapt to write the receptor of a session, as its own
    "receptor PDBQT" button does: POST api/save with the session and kind
    "receptor". ADapt writes the PDBQT and the _rec_H.pdb into a new save
    folder and records the save in session.json. Returns its answer, whose
    "save" is the save's number; raises RuntimeError with ADapt's reason."""
    body = json.dumps({"s": session, "kind": "receptor"}).encode()
    got = call(url, "api/save", data=body, timeout=timeout)
    if not got.get("save"):
        raise RuntimeError("ADapt wrote no receptor")
    return got


def wanted_port():
    """The port asked for in the environment, or 0 to let the system choose.

    A value that is not a port number is ignored rather than refused: a bad
    setting should not stop ADapt from opening at all.
    """
    raw = (os.environ.get(PORT_ENV) or "").strip()
    try:
        n = int(raw)
    except ValueError:
        return 0
    return n if 0 <= n <= 65535 else 0


def address(url):
    """The port in an ADapt address, as an integer, or None."""
    m = re.search(r"://[^/]*:(\d+)", url or "")
    return int(m.group(1)) if m else None


class Runner:
    """One ADapt server for this interface: started on request, then reused,
    and stopped when the interface stops."""

    # A first start imports AGFR and its libraries: 15 seconds on a quiet
    # cluster node with a cold file cache, and over a minute on a busy login
    # node, where a two minute limit was too tight to be safe.
    START_SECONDS = 300

    def __init__(self):
        self.proc = None
        self.url = None
        # A process still starting, so stop() can end it without waiting for
        # the start to finish: that can take minutes, and the interface may
        # be shutting down meanwhile.
        self.starting = None
        # When the start now under way began, and why the last one failed:
        # what the page shows while it waits, read without the lock, which a
        # start holds for as long as it waits.
        self.started_at = None
        self.failed = ""
        self.lock = threading.Lock()

    def running(self):
        return bool(self.url) and self.proc is not None and self.proc.poll() is None

    def phase(self):
        """Where a start has got to: "running" once ADapt has printed its
        address, "starting" while its process is up and has not yet,
        "failed" when the last start gave up (with why), else "none"."""
        if self.running():
            return {"process": "running"}
        proc, began = self.starting, self.started_at
        if proc is not None and proc.poll() is None:
            return {"process": "starting",
                    "seconds": round(time.time() - began, 1) if began else None}
        if self.failed:
            return {"process": "failed", "why": self.failed}
        return {"process": "none"}

    def start(self, home, env_root, outdir, timeout=None, port=None,
              base_path=BASE_PATH):
        """ADapt's address with its token, starting it if it is not running.

        ADapt is started with --base-path, so it answers under /adapt/ as
        well as at its root, and prints its address with the prefix: this
        interface's server forwards /adapt/... to it unchanged.

        Raises RuntimeError with ADapt's own last lines when it does not come
        up, since that is where the reason is.
        """
        with self.lock:
            if self.running():
                return self.url
            os.makedirs(outdir, exist_ok=True)
            want = wanted_port() if port is None else port
            log_path = os.path.join(outdir, "adapt_server.log")
            env = dict(os.environ, ADCPSUITE=env_root, ADAPT_OUTDIR=outdir,
                       # Its address is read from its output, which Python
                       # holds back when it is not writing to a terminal.
                       PYTHONUNBUFFERED="1")
            proc = subprocess.Popen(
                ["sh", os.path.join(home, "bin", "adapt"), "serve",
                 "--host", "127.0.0.1", "--port", str(want), "--outdir", outdir]
                + (["--base-path", base_path] if base_path else []),
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL, env=env, text=True,
                start_new_session=True)
            self.starting = proc
            found, ready, tail = {}, threading.Event(), []

            def pump():
                # Keeps reading for as long as ADapt runs, so its output never
                # fills the pipe and stalls it.
                with open(log_path, "a") as log:
                    for line in proc.stdout:
                        log.write(line)
                        log.flush()
                        tail.append(line.rstrip())
                        del tail[:-8]
                        m = _URL.search(line)
                        if m and "url" not in found:
                            found["url"] = m.group(1)
                            ready.set()
                ready.set()
            threading.Thread(target=pump, daemon=True).start()
            ready.wait(timeout or self.START_SECONDS)
            if "url" not in found:
                self.starting = None
                _stop(proc)
                why = [l for l in tail if l.strip()][-3:]
                self.failed = ("ADapt did not start. Its last words: %s"
                               % (" / ".join(why) or "none; see %s" % log_path))
                raise RuntimeError(self.failed)
            # Running before no longer starting, so the page never sees a
            # moment in which neither is true.
            self.proc, self.url = proc, found["url"]
            self.starting = None
            return self.url

    def stop(self):
        # A start in progress is ended first and without the lock, which the
        # start holds for as long as it waits.
        starting = self.starting
        if starting is not None:
            _stop(starting)
        with self.lock:
            proc, self.proc, self.url = self.proc, None, None
        if proc is not None:
            _stop(proc)


def _stop(proc):
    if proc.poll() is not None:
        return
    # The whole process group: bin/adapt is a shell script that execs Python.
    # Windows has no process groups, and there the process itself is stopped.
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except (OSError, AttributeError):
        proc.terminate()
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (OSError, AttributeError):
            proc.kill()
