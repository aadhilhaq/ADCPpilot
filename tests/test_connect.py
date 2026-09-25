"""The laptop-to-cluster connection, end to end, against a fake SSH server.

By Aadhil Haq

A real login needs a password and a second-factor approval from the
account owner, so it cannot be tested here. Everything around it can:
answering a two-stage challenge, passing the server's wording through
unchanged, not opening a dialog for a message that asks no question,
reading the server's URL line from its output, and carrying HTTP over a
forwarded port.

This stands up an SSH server in-process that uses keyboard-interactive
authentication with a password and a second factor, and a real HTTP server
behind it. The client under test is the real launcher, and the last check
is an actual GET through the forwarded port.

Needs paramiko, the launcher's only dependency. Skipped rather than
failed where it is absent, because the cluster side does not need it.

    python3 tests/test_connect.py
"""
import os
import pathlib
import socket
import sys
import tempfile
import threading
import time
import http.server

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

try:
    import paramiko
except ImportError:
    print("paramiko is not installed; the launcher's connection tests are "
          "skipped. pip install 'adcp-dock[launcher]'")
    raise SystemExit(77)

from adcp_launcher import session, tunnel                      # noqa: E402

# The host key is generated fresh each run and the port is whatever was
# free, so trusting this "host" must not append a dead entry to the user's
# real known_hosts. It goes to a scratch file that is thrown away.
_TMP = tempfile.mkdtemp(prefix="adcp-connect-test-")
tunnel.KNOWN_HOSTS = os.path.join(_TMP, "known_hosts")

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


DUO = ("Duo two-factor login for someone\n\n"
       "Enter a passcode or select one of the following options:\n\n"
       " 1. Duo Push to XXX-XXX-0000\n 2. Phone call to XXX-XXX-0000\n\n"
       "Passcode or option (1-2): ")

#: What the real server sends after the push is approved: a message, and NO
#: prompts. RFC 4256 3.3 says display it and reply empty: the client "MUST
#: NOT prompt the user" for it.
DONE = "Success. Logging you in..."

HOST_KEY = paramiko.RSAKey.generate(2048)


class Server(paramiko.ServerInterface):
    """Keyboard-interactive only: a password, then a second factor."""

    def __init__(self, url_line):
        self.stage = 0
        self.seen = []
        self.dests = {}
        self.commands = []
        self.url_line = url_line

    def get_allowed_auths(self, username):
        return "keyboard-interactive"

    def check_auth_none(self, username):
        return paramiko.AUTH_FAILED

    def check_auth_password(self, username, password):
        return paramiko.AUTH_FAILED          # the site does not offer this

    def check_auth_interactive(self, username, submethods):
        self.stage = 1
        return paramiko.InteractiveQuery("", "", ("Password: ", False))

    def check_auth_interactive_response(self, responses):
        self.seen.append(list(responses))
        if self.stage == 1:
            if responses[0] != "correct horse":
                return paramiko.AUTH_FAILED
            self.stage = 2
            return paramiko.InteractiveQuery("", DUO, ("", True))
        if self.stage == 2:
            if responses[0] != "1":
                return paramiko.AUTH_FAILED
            self.stage = 3
            return paramiko.InteractiveQuery("", DONE)
        return paramiko.AUTH_SUCCESSFUL

    def check_channel_request(self, kind, chanid):
        return paramiko.OPEN_SUCCEEDED

    def check_channel_pty_request(self, *a, **k):
        return True

    def check_channel_direct_tcpip_request(self, chanid, origin, destination):
        # sshd does not merely accept a forwarding channel, it connects to
        # the destination and pumps bytes between the two. Accepting without
        # bridging makes the client look broken when the fake server is.
        self.dests[chanid] = destination
        return paramiko.OPEN_SUCCEEDED

    def check_channel_exec_request(self, channel, command):
        self.commands.append(command.decode())
        threading.Thread(target=self._exec, args=(channel, command),
                         daemon=True).start()
        return True

    def _exec(self, channel, command):
        # Output, then the status, then HALF-close. Not channel.close():
        # paramiko replies to the exec request after this thread has been
        # started, so closing from here can land before that reply reaches
        # the client, and exec_command raises "Channel closed"; a race that
        # shows up on one platform and not another.
        try:
            channel.sendall(b"starting the ADCP Docking server\n")
            time.sleep(0.05)
            channel.sendall((self.url_line + "\n").encode())
            # A real server keeps running; it does not exit after printing.
            for _ in range(200):
                if not channel.active:
                    return
                time.sleep(0.05)
        except Exception:                                       # noqa: BLE001
            pass


def bridge(ch, dest):
    """Carry bytes between a forwarding channel and where it asked for."""
    import select
    try:
        s = socket.create_connection(dest, 10)
    except OSError:
        ch.close()
        return
    try:
        while True:
            r, _, _ = select.select([ch, s], [], [], 1.0)
            if ch in r:
                d = ch.recv(16384)
                if not d:
                    break
                s.sendall(d)
            if s in r:
                d = s.recv(16384)
                if not d:
                    break
                ch.sendall(d)
    except OSError:
        pass
    finally:
        s.close()
        ch.close()


def serve_one(port, server):
    ls = socket.socket()
    ls.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    ls.bind(("127.0.0.1", port))
    ls.listen(1)

    def go():
        sock, _ = ls.accept()
        t = paramiko.Transport(sock)
        t.add_server_key(HOST_KEY)
        try:
            t.start_server(server=server)
        except Exception:                                       # noqa: BLE001
            return
        while t.is_active():
            ch = t.accept(0.2)
            if ch is None:
                continue
            dest = server.dests.pop(ch.get_id(), None)
            if dest:
                threading.Thread(target=bridge, args=(ch, dest),
                                 daemon=True).start()

    threading.Thread(target=go, daemon=True).start()
    return ls


# --- a real HTTP server, standing in for the one on the login node --------
BODY = b"<!doctype html><title>ADCP Docking</title>the page"


class Page(http.server.BaseHTTPRequestHandler):
    def do_GET(self):                                           # noqa: N802
        self.send_response(200)
        self.send_header("Content-Type", "text/html")
        self.send_header("Content-Length", str(len(BODY)))
        self.end_headers()
        self.wfile.write(BODY)

    def log_message(self, *a):
        pass


httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Page)
REMOTE_PORT = httpd.server_address[1]
threading.Thread(target=httpd.serve_forever, daemon=True).start()

TOKEN = "Xq4vT_ab-9"
URL_LINE = "http://127.0.0.1:%d/?t=%s" % (REMOTE_PORT, TOKEN)

# --- the two-stage challenge, answered ------------------------------------
port = tunnel.free_port()
srv = Server(URL_LINE)
listener = serve_one(port, srv)
time.sleep(0.3)

shown = []


def ask(title, instructions, prompts):
    """Stand in for the person: record what was shown, answer in order."""
    shown.append((title, instructions, [p for p, _ in prompts]))
    if instructions and "Duo" in instructions:
        return ["1"]
    return ["correct horse"]


noted = []
t = tunnel.connect("127.0.0.1", "someone", ask, port=port, try_keys=False,
                   log=noted.append)
check(t.is_authenticated(),
      "the launcher did not get through a password-then-Duo challenge")
check(srv.seen == [["correct horse"], ["1"], []],
      "the server saw %r; it should see the password, the Duo choice, then "
      "an EMPTY reply to the message that carries no question" % (srv.seen,))
check(len(shown) == 2,
      "%d dialog(s) were opened for 3 requests. The third carries a message "
      "and no prompts, and RFC 4256 3.3 says show it without prompting: a "
      "box here appears AFTER the person has already approved on their "
      "phone, and the seconds spent dismissing it come out of the same "
      "auth_timeout budget as the approval did" % len(shown))
check(any("Duo" in (i or "") for _t, i, _p in shown),
      "Duo's own wording did not reach the dialog. It lists the methods and "
      "their numbers; a launcher that paraphrases it is asking a user to "
      "choose between options it invented")
check(tunnel.AUTH_PATIENCE >= 120,
      "the two-factor timeout is %ds. It has to cover reading the menu, choosing, "
      "finding a phone and tapping Approve, and 120 is what OpenSSH's own "
      "LoginGraceTime allows" % tunnel.AUTH_PATIENCE)

# --- the server is started, its URL read, and a port forwarded to it ------
INSTALL = "/home/someone/ADCP_docking"
out = []
s = session.start(t, INSTALL, runs="~/adcp_runs",
                  tools_roots=("/opt/adfr", "/opt/amber"),
                  on_output=out.append, forward=tunnel.forward,
                  free_port=tunnel.free_port)
check(srv.commands, "no command was run on the far end")
cmd = srv.commands[0] if srv.commands else ""
check("adcp_dock.cli serve" in cmd,
      "the launcher ran %r, which does not start the server" % cmd[:90])
check("--host 127.0.0.1" in cmd,
      "the server would be started on a public interface; anyone on the "
      "login node could reach the session")
check(cmd.count("--tools-root") == 2,
      "the ADFRsuite and Amber roots did not reach the far end, so the run "
      "would be told those tools are missing")
check(INSTALL in cmd and "PYTHONPATH=" in cmd,
      "the command does not put the checkout on PYTHONPATH")
check(s.token == TOKEN,
      "the token was read as %r, not %r: every request would be refused"
      % (s.token, TOKEN))
check(s.remote_port == REMOTE_PORT,
      "the launcher forwarded to port %r, the server said %r"
      % (s.remote_port, REMOTE_PORT))
check(s.local_port and s.local_port != REMOTE_PORT,
      "no distinct local port was chosen")
check(s.url.startswith("http://127.0.0.1:%d/?t=" % s.local_port),
      "the URL handed to the browser is %r" % s.url)

# --- and HTTP actually crosses it -----------------------------------------
import urllib.request                                          # noqa: E402
try:
    with urllib.request.urlopen(s.url, timeout=20) as r:
        body = r.read()
    check(body == BODY,
          "the forwarded port returned %r; the tunnel carries the "
          "connection but not the bytes" % body[:60])
except Exception as exc:                                        # noqa: BLE001
    check(False, "nothing could be fetched through the forwarded port (%s: "
                 "%s). This is what the launcher provides: the page has "
                 "to be reachable on the laptop." % (type(exc).__name__, exc))

# --- and closing it takes the far end with it -----------------------------
s.stop()
time.sleep(0.4)
try:
    urllib.request.urlopen(s.url, timeout=3)
    check(False, "the forwarded port is still open after stop(); a launcher "
                 "that leaves a listener behind leaves a server running on a "
                 "shared login node too")
except Exception:                                               # noqa: BLE001
    check(True, "")

t.close()
httpd.shutdown()
listener.close()

if fails:
    print("test_connect FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the launcher signs in through a two-factor challenge, starts the server, forwards a "
      "port and carries HTTP over it, %d checks" % ran)
