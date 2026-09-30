"""Run the interface on THIS computer instead of on a cluster.

By Aadhil Haq

The launcher is a small SSH client; it does not carry the docking stack
itself (RDKit, OpenMM, ADFRsuite are far too large to bundle). So "This
computer" starts `adcp-dock local` (installed separately, once) reads the
address it prints, and opens it, exactly as the cluster path reads the
address `adcp-dock serve` prints on a login node.

Where it looks, in order: `adcp-dock` on PATH, then the virtual environment
INSTALL.md tells people to make (~/.adcp-dock/venv). A venv, not a plain
`pip install`, because current Ubuntu, WSL's default, and Homebrew's
Python refuse to install into the system interpreter at all.

On Windows it runs inside WSL: there is no native Windows build of what the
docking needs, and WSL shares localhost with Windows, so the address it
prints opens in the Windows browser unchanged.
"""
import os
import platform
import re
import shutil
import subprocess
import time

URL = re.compile(r"http://127\.0\.0\.1:(\d+)/\?t=([A-Za-z0-9_\-]+)")
VENV_BIN = "~/.adcp-dock/venv/bin/adcp-dock"
ARGS = ["local", "--no-browser", "--port", "0"]

SETUP = {
    "windows": (
        "The docking runs inside WSL. Once:\n\n"
        "  1. Install Ubuntu from the Microsoft Store and open it.\n"
        "  2. In Ubuntu:\n"
        "       sudo apt install -y python3-venv\n"
        "       python3 -m venv ~/.adcp-dock/venv\n"
        "       ~/.adcp-dock/venv/bin/pip install \"adcp-dock[local] @ "
        "git+https://github.com/aadhilhaq/ADCP_docking\"\n\n"
        "Then choose This computer again."),
    "darwin": (
        "Once, in Terminal:\n\n"
        "  softwareupdate --install-rosetta      (Apple silicon only)\n"
        "  python3 -m venv ~/.adcp-dock/venv\n"
        "  ~/.adcp-dock/venv/bin/pip install \"adcp-dock[local] @ "
        "git+https://github.com/aadhilhaq/ADCP_docking\"\n\n"
        "Then choose This computer again."),
    "linux": (
        "Once, in a terminal:\n\n"
        "  python3 -m venv ~/.adcp-dock/venv\n"
        "  ~/.adcp-dock/venv/bin/pip install \"adcp-dock[local] @ "
        "git+https://github.com/aadhilhaq/ADCP_docking\"\n\n"
        "Then choose This computer again."),
}


def setup_text(system=None):
    return SETUP.get((system or platform.system()).lower(), SETUP["linux"])


def command(system=None, which=shutil.which, exists=os.path.isfile):
    """The argv that starts `adcp-dock local` here, or None if it is not
    installed. `which` and `exists` are parameters so this can be tested on
    any machine for any other."""
    system = (system or platform.system()).lower()
    if system == "windows":
        wsl = which("wsl.exe") or which("wsl")
        if not wsl:
            return None
        # A login shell, so ~/.local/bin and a profile's PATH apply; then the
        # documented venv if adcp-dock is not on PATH.
        inner = ("exec \"$(command -v adcp-dock || echo %s)\" %s"
                 % (VENV_BIN, " ".join(ARGS)))
        return [wsl, "-e", "bash", "-lc", inner]
    exe = which("adcp-dock")
    if not exe:
        venv = os.path.expanduser(VENV_BIN)
        exe = venv if exists(venv) else None
    return [exe] + ARGS if exe else None


class Local:
    """What the launcher holds while the interface runs here. Has the same
    two things the cluster session has: a url, and stop()."""

    def __init__(self, proc, url):
        self.proc, self.url = proc, url

    def stop(self):
        if self.proc.poll() is None:
            self.proc.terminate()
            try:
                self.proc.wait(10)
            except subprocess.TimeoutExpired:
                self.proc.kill()


def start(on_output=None, timeout=1800, cmd=None, popen=subprocess.Popen):
    """Start it, and return a Local once it prints its address.

    The timeout is long on purpose: the FIRST start downloads and installs
    ADFRsuite, which takes minutes, and everything it says meanwhile is
    passed to on_output so the window is not silent while it happens.
    """
    say = on_output or (lambda *a: None)
    cmd = cmd or command()
    if not cmd:
        raise RuntimeError("adcp-dock is not installed on this computer.\n\n"
                           + setup_text())
    try:
        proc = popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                     text=True, bufsize=1)
    except OSError as exc:
        # NOT A TRACEBACK. The command could not be started at all: on
        # Windows, typically wsl.exe present with no Linux installed in it.
        # That is a setup step, and the window should say which one.
        raise RuntimeError("Could not start adcp-dock on this computer (%s)."
                           "\n\n%s" % (exc, setup_text()))
    # A reader thread, so the timeout can fire. Reading the pipe directly
    # blocks until a line arrives, and a child that prints nothing,
    # which is what a buffered stdout looks like, would block
    # indefinitely with the deadline never checked.
    import queue
    import threading
    lines = queue.Queue()

    def pump():
        for ln in proc.stdout:
            lines.put(ln)
        lines.put(None)
    threading.Thread(target=pump, daemon=True).start()

    tail, deadline = [], time.time() + timeout
    while True:
        left = deadline - time.time()
        if left <= 0:
            proc.terminate()
            raise RuntimeError(
                "adcp-dock local printed no address in %d minute(s). What it "
                "said last:\n\n%s" % (max(1, timeout // 60),
                                         "\n".join(tail) or "(nothing)"))
        try:
            line = lines.get(timeout=min(left, 1.0))
        except queue.Empty:
            continue
        if line is None:
            break
        line = line.rstrip()
        if line:
            say(line)
            tail = (tail + [line])[-15:]
        m = URL.search(line)
        if m:
            return Local(proc, m.group(0))
    proc.wait()
    raise RuntimeError("adcp-dock local stopped before it was ready (exit "
                       "%s). What it said last:\n\n%s"
                       % (proc.returncode, "\n".join(tail)))
