"""SSH to the cluster, through two-factor authentication, with a forwarded port.

By Aadhil Haq

Three things here differ from paramiko's defaults:

Two-factor approval needs more than thirty seconds. paramiko's
auth_timeout must cover reading the prompt, choosing a method, and
approving on a phone, so it is raised to two minutes.

A keyboard-interactive challenge carrying no prompts is a message to
display, not something to prompt for: RFC 4256 section 3.3 says the client
must not prompt in that case, so it is written to the log.

Host keys are checked against the user's known_hosts and a change is
reported rather than accepted silently.
"""
import os
import socket
import socketserver
import time
import threading

try:
    import paramiko
except ImportError:                                            # pragma: no cover
    paramiko = None

#: OpenSSH's LoginGraceTime, which is what the far end allows.
AUTH_PATIENCE = 120
KNOWN_HOSTS = os.path.expanduser("~/.ssh/known_hosts")


class TunnelError(RuntimeError):
    pass


def need_paramiko():
    if paramiko is None:
        raise TunnelError(
            "paramiko is not installed. The launcher needs it to make the SSH "
            "connection:\n\n    pip install paramiko")


def _host_id(host, port):
    return host if int(port) == 22 else "[%s]:%d" % (host, int(port))


def check_host_key(t, host, port=22, log=None):
    """Compare the key offered against ~/.ssh/known_hosts.

    Reported rather than enforced: a cluster that has been rebuilt offers a
    new key legitimately, and a launcher that refuses to connect and does not
    say why is worse than one that says what changed.
    """
    key = t.get_remote_server_key()
    try:
        hosts = paramiko.hostkeys.HostKeys(KNOWN_HOSTS)
    except (OSError, paramiko.SSHException):
        return None
    entry = hosts.lookup(_host_id(host, port)) or hosts.lookup(host)
    if not entry:
        return ("%s is not in known_hosts yet; its key fingerprint is %s"
                % (host, key.get_fingerprint().hex()))
    known = entry.get(key.get_name())
    if known is None:
        return None
    if known.asbytes() != key.asbytes():
        return ("THE HOST KEY FOR %s HAS CHANGED. This is either a rebuilt "
                "cluster or something between you and it. Check with the "
                "administrators before continuing." % host)
    if log:
        log("host key matches known_hosts")
    return None


def _handler(ask, log=None):
    """The callback paramiko hands every keyboard-interactive request."""
    def handler(title, instructions, prompts):
        if not prompts:
            if instructions and instructions.strip() and log:
                log(instructions.strip())
            return []
        return ask(title, instructions, [(p, echo) for p, echo in prompts])
    return handler


def _patience(t, timeout):
    t.banner_timeout = timeout
    t.handshake_timeout = timeout
    t.auth_timeout = AUTH_PATIENCE


def _keys():
    out = []
    d = os.path.expanduser("~/.ssh")
    for name in ("id_ed25519", "id_rsa", "id_ecdsa"):
        p = os.path.join(d, name)
        if os.path.isfile(p):
            out.append(p)
    return out


def connect(host, user, ask, port=22, timeout=40, try_keys=True, log=None):
    """An authenticated transport, or a TunnelError saying why not.

    `ask(title, instructions, prompts)` is called for every challenge. Prompts
    are (text, echo) pairs, echo False for a password; it returns one answer
    per prompt. Keys are tried first, because a site that accepts a key usually
    skips the second factor.
    """
    need_paramiko()
    say = log or (lambda *a: None)
    try:
        sock = socket.create_connection((host, int(port)), timeout=timeout)
    except OSError as exc:
        raise TunnelError("could not reach %s port %s: %s" % (host, port, exc))
    t = paramiko.Transport(sock)
    _patience(t, timeout)
    try:
        t.start_client(timeout=timeout)
    except Exception as exc:                                    # noqa: BLE001
        raise TunnelError("the SSH handshake with %s failed: %s" % (host, exc))
    warn = check_host_key(t, host, port, log=say)
    if warn:
        say(warn)
    if try_keys:
        for path in _keys():
            for cls in (paramiko.Ed25519Key, paramiko.RSAKey,
                        paramiko.ECDSAKey):
                try:
                    k = cls.from_private_key_file(path)
                except Exception:                               # noqa: BLE001
                    continue
                try:
                    t.auth_publickey(user, k)
                    say("signed in with %s" % os.path.basename(path))
                    return t
                except paramiko.SSHException:
                    break
    try:
        say("asking %s for a challenge" % host)
        t.auth_interactive(user, _handler(ask, say))
        if t.is_authenticated():
            return t
    except paramiko.BadAuthenticationType as exc:
        if "password" not in (exc.allowed_types or []):
            raise TunnelError("%s will not accept a password or a challenge "
                              "from here (it allows: %s)"
                              % (host, ", ".join(exc.allowed_types or [])))
    except paramiko.AuthenticationException as exc:
        raise TunnelError("authentication failed: %s" % exc)
    except Exception as exc:                                    # noqa: BLE001
        raise TunnelError("authentication failed: %s" % exc)
    try:
        answers = ask("Password", "", [("Password: ", False)])
        t.auth_password(user, answers[0])
    except Exception as exc:                                    # noqa: BLE001
        raise TunnelError("authentication failed: %s" % exc)
    if not t.is_authenticated():
        raise TunnelError("authentication failed")
    return t


def run(t, command, timeout=60):
    """One command, its output, and its exit status."""
    ch = t.open_session(timeout=timeout)
    ch.settimeout(timeout)
    ch.exec_command(command)
    out = []
    while True:
        if ch.recv_ready():
            out.append(ch.recv(65536).decode("utf8", "replace"))
        elif ch.exit_status_ready():
            while ch.recv_ready():
                out.append(ch.recv(65536).decode("utf8", "replace"))
            break
        else:
            # a short sleep, not select() on the channel: that makes paramiko keep a pipe
            # for it, which its own thread can write into after the channel has closed
            time.sleep(0.05)
    status = ch.recv_exit_status()
    ch.close()                  # its session is given back now, not when collected
    return "".join(out), status


class _Handler(socketserver.BaseRequestHandler):
    """One browser connection, carried over a channel of the shared connection.

    A thread each way, each waiting on its own side. It used to wait on both with select(),
    which asks paramiko for the channel's file number, and for that paramiko keeps a pipe beside
    the channel that its own thread writes to whenever data arrives. A browser that gives up on
    a request - it does all the time, whenever a page loads something else instead - closed the
    channel, and that pipe, while the rest of the answer was still arriving; paramiko's thread
    then wrote into a pipe no longer there, took the error for its connection failing ("Broken
    pipe", "Bad file descriptor") and ended it, and the page stopped answering as if signed out.
    Nothing here asks for a file number, so paramiko makes no pipe.
    """

    def handle(self):
        try:
            ch = self.server.transport.open_channel(
                "direct-tcpip", (self.server.remote_host, self.server.remote_port),
                self.request.getpeername())
        except Exception:                                       # noqa: BLE001
            return                      # the far end is gone; the browser retries
        if ch is None:
            return
        browser = self.request

        def down():
            # The answer, to the browser, until the far end or the browser is done.
            try:
                while True:
                    data = ch.recv(65536)
                    if not data:
                        break
                    browser.sendall(data)
            except Exception:                                   # noqa: BLE001
                pass                    # a closed browser tab is not an error
            finally:
                try:
                    browser.shutdown(socket.SHUT_WR)
                except OSError:
                    pass

        answer = threading.Thread(target=down, daemon=True)
        answer.start()
        try:
            while True:
                data = browser.recv(65536)
                if not data:
                    break
                ch.sendall(data)
        except Exception:                                       # noqa: BLE001
            pass
        finally:
            try:
                ch.shutdown_write()     # the browser has said all it will: the far end may finish
            except Exception:                                   # noqa: BLE001
                pass
        # The rest of the answer is waited for, and a far end that keeps a connection the browser
        # has left open for no reason is let go of after a minute: closing the channel ends the wait.
        answer.join(60)
        try:
            ch.close()
        except Exception:                                       # noqa: BLE001
            pass


class _Forwarder(socketserver.ThreadingTCPServer):
    daemon_threads = True
    allow_reuse_address = True


def forward(t, local_port, remote_port, remote_host="127.0.0.1"):
    """A local port that reaches a port on the far end. Returns the server."""
    srv = _Forwarder(("127.0.0.1", int(local_port)), _Handler)
    srv.transport = t
    srv.remote_host = remote_host
    srv.remote_port = int(remote_port)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv


def free_port():
    s = socket.socket()
    s.bind(("127.0.0.1", 0))
    p = s.getsockname()[1]
    s.close()
    return p
