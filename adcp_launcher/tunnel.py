"""SSH to the cluster, through two-factor authentication, with a forwarded port.

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
import select
import socket
import socketserver
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
            select.select([ch], [], [], 0.2)
    return "".join(out), ch.recv_exit_status()


class _Handler(socketserver.BaseRequestHandler):
    def handle(self):
        try:
            ch = self.server.transport.open_channel(
                "direct-tcpip", (self.server.remote_host,
                                 self.server.remote_port),
                self.request.getpeername())
        except Exception:                                       # noqa: BLE001
            return
        if ch is None:
            return
        try:
            while True:
                r, _w, _x = select.select([self.request, ch], [], [], 1.0)
                if self.request in r:
                    data = self.request.recv(16384)
                    if not data:
                        break
                    ch.sendall(data)
                if ch in r:
                    data = ch.recv(16384)
                    if not data:
                        break
                    self.request.sendall(data)
        except OSError:
            pass
        finally:
            try:
                ch.close()
            except Exception:                                   # noqa: BLE001
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
