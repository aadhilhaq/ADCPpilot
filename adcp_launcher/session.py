"""Starting the tool on the cluster and getting a URL back.

By Aadhil Haq

The server prints one line holding both the port it got and the token that
gates it. That line is the whole handshake: nothing is guessed, nothing is
assumed about which port was free, and a server that failed to start says so
in its own words instead of leaving the launcher waiting on a port that will
never open.
"""
import re
import shlex
import time

from . import clusters

URL = re.compile(r"http://127\.0\.0\.1:(\d+)/\?t=([A-Za-z0-9_\-]+)")

#: Printed by the remote command when the Python it found is too old. A
#: marker rather than a sentence, so the launcher can recognize it without
#: matching on wording that may be translated or reworded.
OLD_PYTHON_MARKER = "ADCP_PYTHON_TOO_OLD"

#: Where the tool might be, when the person has not said. Tried in order.
GUESSES = (
    "~/ADCP_docking",
    "$SCRATCH/ADCP_docking",
)


#: The oldest Python this tool's code can even be PARSED by. The engine
#: modules open with `from __future__ import annotations`, which is a
#: SyntaxError before 3.7, and the package declares 3.9.
MIN_PYTHON = (3, 9)


def command(install, runs="", python="", port=0, tools_roots=(), env=None,
            defaults=None):
    """The line run on the far end.

    The server is started on the login node: it reads files, validates a form
    and submits batch jobs, and everything heavy it starts goes to the queue.

    It sources activate.sh rather than naming a Python. install.sh writes that
    file, so it is the one place the installed layout is recorded.

    The Python version is checked first, by a one-line program that any
    Python can parse. A cluster whose default python3 is too old would
    otherwise fail on a SyntaxError deep in an import, which reaches the user
    as a traceback instead of a sentence.

    `defaults` are set only where activate.sh has not set them already: what
    a named cluster's table suggests (where runs go, the Slurm account),
    which the installation's own settings override.
    """
    where = clusters.remote_path(install)
    bits = ["-m", "adcp_pilot.cli", "serve",
            "--host", "127.0.0.1", "--port", str(port)]
    quoted = [shlex.quote(b) for b in bits]
    if runs:
        quoted += ["--runs", clusters.remote_path(runs)]
    for r in tools_roots:
        quoted += ["--tools-root", clusters.remote_path(r)]
    # A named interpreter takes precedence, for environments this does
    # not describe.
    pick = ('PY=%s' % shlex.quote(python) if python
            else 'PY="${ADCP_PYTHON:-python3}"')
    return (
        'cd %(where)s || exit 1; '
        '[ -f ./activate.sh ] && . ./activate.sh; '
        '%(pick)s; '
        '"$PY" -c %(check)s 2>/dev/null || { '
        'echo "%(marker)s $("$PY" -V 2>&1)"; exit 78; }; '
        '%(defaults)s%(env)sPYTHONPATH=%(where)s exec "$PY" %(args)s'
        % {"where": where, "pick": pick, "marker": OLD_PYTHON_MARKER,
           "defaults": "".join(
               '[ -n "${%s:-}" ] || export %s=%s; '
               % (k, k, clusters.remote_path(v) if k == "ADCP_RUNS" else shlex.quote(v))
               for k, v in sorted((defaults or {}).items()) if v),
           "env": "".join("export %s=%s; " % (k, shlex.quote(v))
                          for k, v in sorted((env or {}).items()) if v),
           "check": shlex.quote(
               "import sys; sys.exit(0 if sys.version_info[:2] >= %r else 1)"
               % (MIN_PYTHON,)),
           "args": " ".join(quoted)})


class Session:
    """A running server on the far end, and the local port that reaches it."""

    def __init__(self, channel, forwarder, local_port, remote_port, token,
                 host):
        self.channel = channel
        self.forwarder = forwarder
        self.local_port = local_port
        self.remote_port = remote_port
        self.token = token
        self.host = host

    @property
    def url(self):
        return "http://127.0.0.1:%d/?t=%s" % (self.local_port, self.token)

    def stop(self):
        for close in (getattr(self.channel, "close", None),
                      getattr(self.forwarder, "shutdown", None),
                      getattr(self.forwarder, "server_close", None)):
            try:
                if close:
                    close()
            except Exception:                                   # noqa: BLE001
                pass


def start(t, install, runs="", python="python3", tools_roots=(),
          timeout=90, on_output=None, forward=None, free_port=None,
          port=0, env=None, defaults=None):
    """Run the server, read its URL, forward a local port to it."""
    say = on_output or (lambda *a: None)
    # Port 0 means the server picks a free one and prints it, which is
    # the only choice that cannot collide on a shared login node. A
    # number here is a deliberate request for a stable one.
    line = command(install, runs, python, port, tools_roots, env=env, defaults=defaults)
    say("$ " + line)
    ch = t.open_session()
    ch.get_pty()          # so the server dies with the connection
    ch.exec_command(line)
    ch.settimeout(1.0)

    buf, deadline = "", time.time() + timeout
    while time.time() < deadline:
        if ch.exit_status_ready() and not ch.recv_ready():
            raise RuntimeError(_why(buf, install))
        try:
            chunk = ch.recv(65536).decode("utf8", "replace")
        except Exception:                                       # noqa: BLE001
            chunk = ""
        if chunk:
            buf += chunk
            for ln in chunk.splitlines():
                if ln.strip():
                    say(ln.rstrip())
        m = URL.search(buf)
        if m:
            remote_port, token = int(m.group(1)), m.group(2)
            local = (free_port or (lambda: 0))()
            srv = forward(t, local, remote_port)
            local = srv.server_address[1]
            say("forwarding 127.0.0.1:%d to the cluster's %d"
                % (local, remote_port))
            return Session(ch, srv, local, remote_port, token, install)
        time.sleep(0.2)
    raise RuntimeError(_why(buf, install) or
                       "the server did not report a URL within %ds. What it "
                       "said:\n%s" % (timeout, _tail(buf)))


def _why(buf, install):
    """What went wrong, in terms that can be acted on."""
    low = buf.lower()
    if OLD_PYTHON_MARKER in buf:
        found = buf.split(OLD_PYTHON_MARKER, 1)[1].strip().splitlines()[:1]
        return ("the Python on the cluster is too old: %s, and this needs "
                "%d.%d or newer.\n\n"
                "The fix is to finish the install, which writes an "
                "activate.sh naming a new enough Python:\n"
                "    cd %s && ./install.sh <shared-prefix>\n\n"
                "Then connect again; the launcher sources that file."
                % (found[0] if found else "unknown", MIN_PYTHON[0],
                   MIN_PYTHON[1], install))
    if "future feature annotations is not defined" in low:
        return ("the cluster ran this with a Python older than %d.%d, which "
                "cannot parse it.\n\nFinish the install so activate.sh "
                "exists and names the right interpreter:\n"
                "    cd %s && ./install.sh <shared-prefix>"
                % (MIN_PYTHON[0], MIN_PYTHON[1], install))
    if "no such file or directory" in low and "activate.sh" in low:
        return ("the tool is there but has not been installed: no "
                "activate.sh in %s.\n    cd %s && ./install.sh "
                "<shared-prefix>" % (install, install))
    if "no module named 'adcp_pilot'" in low or "no module named adcp_pilot" in low:
        return ("the tool was not found at %s. Check the install path, or "
                "clone it there:\n    git clone "
                "https://github.com/aadhilhaq/ADCP_docking %s" % (install, install))
    if "no such file or directory" in low and "cd " in low:
        return "there is no directory %s on the cluster" % install
    if "address already in use" in low:
        return "the port the server chose was taken; try again"
    if "permission denied" in low:
        return "permission denied on the far end:\n%s" % _tail(buf)
    if buf.strip():
        return "the server stopped before it printed a URL:\n%s" % _tail(buf)
    return ""


def _tail(buf, n=14):
    lines = [l for l in buf.splitlines() if l.strip()]
    return "\n".join(lines[-n:])
