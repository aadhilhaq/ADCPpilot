"""ADapt, a separate tool for looking inside a structure before docking.

By Aadhil Haq

ADapt, by Michel Sanner, shows what a PDB file holds: each biological
assembly, alternate locations, side chains with missing atoms, missing
residues, ligands, cofactors and metals. It writes receptor and ligand PDBQT
files the way AGFR does.

It drives AGFR's own code, so it runs on the ADCP 1.1 environment's Python and
is started here as its own program with that environment named. Nothing of it
is imported or shipped with this package; it is found where it was installed:

    git clone https://github.com/michelsanner/adapt

and then ADAPT_HOME set to that folder, or the folder placed beside this
checkout or beside ADFRsuite.
"""
import os
import re
import signal
import subprocess
import threading

REPO = "https://github.com/michelsanner/adapt"
# The port ADapt listens on. Zero means the operating system picks a free one,
# which is right on a shared login node. Anyone reaching this interface through
# an ssh tunnel wants the same port every time instead, since the tunnel names
# it: ADCP_ADAPT_PORT=40805 pins it.
PORT_ENV = "ADCP_ADAPT_PORT"

# The line ADapt prints once it is serving: its address with its token.
_URL = re.compile(r"(https?://\S+/\?t=\S+)")


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
    # cluster node with a cold file cache, and over a minute on a busy Grace
    # login node, where a two minute limit was too tight to be safe.
    START_SECONDS = 300

    def __init__(self):
        self.proc = None
        self.url = None
        # A process still starting, so stop() can end it without waiting for
        # the start to finish: that can take minutes, and the interface may
        # be shutting down meanwhile.
        self.starting = None
        self.lock = threading.Lock()

    def running(self):
        return bool(self.url) and self.proc is not None and self.proc.poll() is None

    def start(self, home, env_root, outdir, timeout=None, port=None):
        """ADapt's address with its token, starting it if it is not running.

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
                 "--host", "127.0.0.1", "--port", str(want), "--outdir", outdir],
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
            self.starting = None
            if "url" not in found:
                _stop(proc)
                why = [l for l in tail if l.strip()][-3:]
                raise RuntimeError("ADapt did not start. Its last words: %s"
                                   % (" / ".join(why) or "none; see %s" % log_path))
            self.proc, self.url = proc, found["url"]
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
