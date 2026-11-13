"""ADapt is found, started with the right environment, reused and stopped,
reached through this server under its own base path, handed the receptor, and
what it saves is read back.

By Aadhil Haq

ADapt is a separate program that the interface starts on request, so what is
tested is the handover: which folder counts as ADapt, what it is started with
(the ADCP 1.1 environment, a sessions folder that is not the home directory,
output that is not held back), that its address is read from its output, that
a second request reuses it, and that it stops. A stand-in script plays ADapt.

    python3 tests/test_adapt.py
"""
import os
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_pilot import adapt  # noqa: E402

import threading  # noqa: E402

fails, ran = [], 0


def _try(fn):
    """Run fn and return its result or the exception it raised."""
    try:
        return fn()
    except Exception as exc:                                    # noqa: BLE001
        return exc


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


def fake_adapt(folder, speaks=True, embeddable=True):
    """A folder that looks like an ADapt checkout, whose bin/adapt records
    how it was started and then prints an address, or dies without one. An
    embeddable one takes --base-path, as ADapt has since 6dae90b."""
    os.makedirs(os.path.join(folder, "bin"), exist_ok=True)
    os.makedirs(os.path.join(folder, "adapt"), exist_ok=True)
    open(os.path.join(folder, "adapt", "__main__.py"), "w").write(
        "s.add_argument('--base-path')\n" if embeddable else "")
    script = ["#!/bin/sh",
              'mkdir -p "$ADAPT_OUTDIR"',
              'echo "$ADCPSUITE|$ADAPT_OUTDIR|$PYTHONUNBUFFERED|$*" > "$ADAPT_OUTDIR/started_with.txt"']
    if speaks:
        script += ['echo "ADapt 0.1.0 is at"', 'echo ""',
                   'echo "    http://127.0.0.1:45678/adapt/?t=abc123"', "exec sleep 60"]
    else:
        script += ['echo "no ADCP environment found" >&2', "exit 2"]
    path = os.path.join(folder, "bin", "adapt")
    open(path, "w").write("\n".join(script) + "\n")
    os.chmod(path, 0o755)
    return folder


tmp = tempfile.mkdtemp(prefix="adapt_test_")
os.environ["HOME"] = tmp                  # nothing may fall back to a real home
for k in ("ADAPT_HOME", "ADCP_HOME", "ADCPSUITE_HOME"):
    os.environ.pop(k, None)
# A checkout on the lab's cluster sits beside tools/docking, whose suite and
# ADapt are found where a test says there are none; a copy in a folder with
# nothing beside it is what "not installed" means here.
from adcp_pilot import adcp11  # noqa: E402

adcp11._HERE = adapt._HERE = os.path.join(tmp, "alone", "ADCP_docking")

# --- which folder counts as ADapt -----------------------------------------
home = fake_adapt(os.path.join(tmp, "adapt"))
check(adapt.is_adapt(home), "a folder with bin/adapt and adapt/__main__.py is ADapt")
check(not adapt.is_adapt(tmp), "a folder without them is not")
os.environ["ADAPT_HOME"] = home
check(adapt.find() == home, "ADAPT_HOME is where it is looked for first")
os.environ.pop("ADAPT_HOME")
check(adapt.find(roots=[tmp]) == home, "an adapt folder in a tools folder is found")

# --- started with the right environment, read, reused, stopped ------------
out = os.path.join(tmp, "runs", "adapt_sessions")
r = adapt.Runner()
t0 = time.time()
url = r.start(home, "/envs/adcpsuite", out, timeout=20)
check(url == "http://127.0.0.1:45678/adapt/?t=abc123",
      "the address is read from ADapt's output, base path and all (%r)" % url)
check(time.time() - t0 < 15, "and read as soon as it is printed")
check(r.running(), "it counts as running once it has printed its address")
for _ in range(50):
    if os.path.isfile(os.path.join(out, "started_with.txt")):
        break
    time.sleep(0.1)
env_root, outdir, unbuffered, args = open(os.path.join(out, "started_with.txt")).read().strip().split("|")
check(env_root == "/envs/adcpsuite", "ADapt is given the ADCP 1.1 environment as ADCPSUITE")
check(outdir == out and "--outdir %s" % out in args,
      "its sessions go where the interface says, not to the home directory")
check(unbuffered == "1", "its output is not held back, or its address never arrives")
check("serve" in args and "--host 127.0.0.1" in args and "--port 0" in args,
      "it serves on this machine only, on a free port (%s)" % args)
check("--base-path /adapt" in args,
      "it is started under this server's /adapt/, as its README says to embed it (%s)" % args)
check(os.path.isfile(os.path.join(out, "adapt_server.log")), "its output is kept in a log")
pid = r.proc.pid
check(r.start(home, "/envs/adcpsuite", out) == url and r.proc.pid == pid,
      "a second request reuses the running ADapt instead of starting another")
proc = r.proc
r.stop()
check(proc.poll() is not None, "stopping it ends the process")
check(not r.running(), "and it no longer counts as running")

# --- the port a tunnel can name -------------------------------------------
# Anyone reaching the interface through ssh forwards one port by name, so a
# port that changes at every start cannot be forwarded in advance.
check(adapt.wanted_port() == 0, "with nothing set, the system picks the port")
check(adapt.address("http://127.0.0.1:40805/?t=x") == 40805, "the port is read back from the address")
check(adapt.address("not an address") is None, "and nothing is read from a line that has none")
os.environ[adapt.PORT_ENV] = "40805"
check(adapt.wanted_port() == 40805, "%s pins it" % adapt.PORT_ENV)
os.environ[adapt.PORT_ENV] = "not a port"
check(adapt.wanted_port() == 0, "a setting that is not a port is ignored, not refused")
os.environ[adapt.PORT_ENV] = "40806"
pinned_out = os.path.join(tmp, "pinned_out")
pr = adapt.Runner()
pr.start(home, "/envs/adcpsuite", pinned_out, timeout=20)
check("--port 40806" in open(os.path.join(pinned_out, "started_with.txt")).read(),
      "and the pinned port is the one ADapt is asked for")
pr.stop()
del os.environ[adapt.PORT_ENV]

# --- a start that fails says why ------------------------------------------
bad = fake_adapt(os.path.join(tmp, "bad", "adapt"), speaks=False)
try:
    adapt.Runner().start(bad, "/envs/adcpsuite", os.path.join(tmp, "bad_out"), timeout=20)
    check(False, "a start that prints no address is reported as a failure")
except RuntimeError as exc:
    check("no ADCP environment found" in str(exc),
          "the failure carries ADapt's own last words (%s)" % exc)

# --- the interface: offered only when it can run ---------------------------
from adcp_pilot import server  # noqa: E402

app = server.App(os.path.join(tmp, "runs"))
st = app._adapt()
check(st["available"] is False and "not installed" in st["why"],
      "without ADapt the interface says it is not installed (%s)" % st)
check("error" in app.adapt_open(), "and opening it returns an error, not a crash")

os.environ["ADAPT_HOME"] = home
app = server.App(os.path.join(tmp, "runs"))
st = app._adapt()
check(st["available"] is False and "ADCP 1.1" in st["why"],
      "with ADapt but without ADCP 1.1 it says what is missing (%s)" % st)

# A stand-in ADCP 1.1 environment, enough for the interface to find one.
suite = os.path.join(tmp, "adcpsuite")
os.makedirs(os.path.join(suite, "bin"))
site = os.path.join(suite, "lib", "python3.7", "site-packages", "ADCP")
os.makedirs(site)
open(os.path.join(suite, "bin", "adcp"), "w").write("#!/bin/sh\n")
open(os.path.join(site, "runADCP.py"), "w").write('version="1.1.%d"\n')
open(os.path.join(site, "utils.py"), "w").write("BUILD_VERSION = 21\n")
os.environ["ADCPSUITE_HOME"] = suite
app = server.App(os.path.join(tmp, "runs"))
check(app._adapt()["available"] is True, "with both it is offered")
res = app.adapt_open()
check(res.get("url") == "/adapt/?t=abc123",
      "and opening it returns a path on THIS server, with ADapt's token (%s)" % res)
check(res.get("direct") == "http://127.0.0.1:45678/adapt/?t=abc123",
      "and ADapt's own address beside it")
env_root = open(os.path.join(tmp, "runs", "adapt_sessions", "started_with.txt")).read().split("|")[0]
check(os.path.realpath(env_root) == os.path.realpath(suite),
      "the interface hands ADapt the ADCP 1.1 environment it found (%s)" % env_root)
check(not any(os.path.basename(x["rundir"]) == "adapt_sessions" for x in app.list_runs()["runs"]),
      "ADapt's sessions folder is not listed as a run")
app.adapt_stop()
check(not app._adapt_runner.running(), "the interface stops ADapt when it stops")

old_home = fake_adapt(os.path.join(tmp, "old", "adapt"), embeddable=False)
check(not adapt.embeddable(old_home) and adapt.embeddable(home),
      "a checkout without --base-path is told apart from one with it")
os.environ["ADAPT_HOME"] = old_home
st_old = server.App(os.path.join(tmp, "runs_old"))._adapt()
check(st_old["available"] is False and "older" in st_old["why"] and "git -C" in st_old["why"],
      "an ADapt older than the embedding interface is not offered, with how to update it (%s)"
      % st_old)
os.environ["ADAPT_HOME"] = home

# --- started ahead of the click, but only when the page asks ----------------
warm_app = server.App(os.path.join(tmp, "runs_warm"))
warm_app.tools()
check(getattr(warm_app, "_adapt_runner", None) is None,
      "listing the tools does not start ADapt: only the page's request does")
check(warm_app.adapt_warm() == {"warming": True}, "the page's request starts it")
for _ in range(100):
    r_ = getattr(warm_app, "_adapt_runner", None)
    if r_ is not None and r_.running():
        break
    time.sleep(0.1)
check(warm_app._adapt_runner.running(), "and it is running by the time anyone clicks")
pid_ = warm_app._adapt_runner.proc.pid
check(warm_app.adapt_warm() == {"warming": True} and warm_app._adapt_runner.proc.pid == pid_,
      "asking again starts nothing new")
check(warm_app.adapt_open().get("url") == "/adapt/?t=abc123",
      "and the button is answered at once by the one already running")
warm_app.adapt_stop()

# --- a start still in progress can be stopped at once ------------------------
# Stopping used to wait for the start's lock, which a slow start holds for up
# to START_SECONDS, and a process still starting was not tracked, so closing
# the interface could leave it behind on a login node.
slow = os.path.join(tmp, "slow_adapt")
fake_adapt(slow)
open(os.path.join(slow, "bin", "adapt"), "w").write(
    '#!/bin/sh\nmkdir -p "$ADAPT_OUTDIR"\necho "loading AGFR"\nexec sleep 120\n')
sr = adapt.Runner()
err_ = []
th_ = threading.Thread(target=lambda: err_.append(
    _try(lambda: sr.start(slow, "/envs/adcpsuite", os.path.join(tmp, "slow_out"), timeout=100))))
th_.start()
for _ in range(50):
    if sr.starting is not None:
        break
    time.sleep(0.1)
starting_ = sr.starting
t0_ = time.time()
sr.stop()
th_.join(20)
check(starting_ is not None and starting_.poll() is not None,
      "stopping ends a process that is still starting")
check(time.time() - t0_ < 15 and not th_.is_alive(),
      "without waiting for the start to give up (%.1f s)" % (time.time() - t0_))
check(err_ and isinstance(err_[0], RuntimeError), "and the start reports that it did not come up")

# --- the page offers it only when the server does --------------------------
from adcp_pilot import page  # noqa: E402

html = page.PAGE if hasattr(page, "PAGE") else open(page.__file__).read()
check("id=adaptrow" in html and "id=adaptopen" in html, "the Receptor card has an ADapt row and button")
appjs = open(os.path.join(ROOT, "adcp_pilot", "static", "app.js")).read()
check("d.adapt && d.adapt.available" in appjs, "the row is shown only when the server offers ADapt")
check("api('/api/adapt'" in appjs, "the button asks the server to start it")
check("Open in ADapt" in html and "id=adaptuse" in html and "id=adaptdrop" in html,
      "the Receptor card opens the receptor in ADapt and shows the receptor read back")
check("/api/adapt/saved" in appjs and "/api/adapt/state" in appjs
      and "async function adaptPoll" in appjs,
      "the page watches the session for changes and for a receptor save")
check("window.open('', '_blank')" in appjs,
      "on a phone the tab is opened inside the click, where a browser allows it")

# --- ADapt in the card, opened by itself -------------------------------------
import re as _re  # noqa: E402
frame = _re.search(r"<iframe id=adaptframe[^>]*>", html)
check(frame is not None, "the Receptor card holds a frame for ADapt")
sandbox = _re.search(r'sandbox="([^"]*)"', frame.group(0)).group(1).split() if frame else []
check("allow-scripts" in sandbox and "allow-same-origin" in sandbox,
      "ADapt runs in the frame, on this page's origin through /adapt/ (%s)" % sandbox)
check(not any(t.startswith("allow-top-navigation") for t in sandbox),
      "and the frame cannot navigate this page (%s)" % sandbox)
check("allow-popups-to-escape-sandbox" in sandbox and "allow-downloads" in sandbox,
      "its Docs link and its downloads still work")
check("adaptAuto(d);" in appjs and "function adaptAuto" in appjs,
      "a structure loaded opens in ADapt by itself")
check("matchMedia('(max-width: 700px)')" in appjs and "id=adaptopen" in html,
      "on a phone the button stands in for the frame")
check("id=adaptfull" in html and "Open full screen" in html and "id=adaptfold" in html,
      "the frame can be opened full screen, and folded away")
check("Docking will use <b>this page\\u2019s preparation</b>" in appjs
      and "Docking will use <b>the receptor as set in ADapt</b>" in appjs,
      "the card always says which receptor will be docked")
check("Use this page&rsquo;s preparation instead" in html
      and "Use ADapt&rsquo;s receptor now" in html,
      "with a way back, and a way to have ADapt write its receptor at once")
check("/api/adapt/write" in appjs and "await adaptWriteNow()" in appjs,
      "submitting with ADapt's receptor has ADapt write it first")
check("localStorage.setItem('adapt.theme'" in appjs and "adaptFrameTheme();" in appjs,
      "the frame follows this page's light or dark theme")
check("ssh -L" not in appjs,
      "nobody is asked to forward a second port: ADapt comes through this page's own")
check("new URL(r.url, location.href)" in appjs,
      "the tab is sent to an absolute address, since it began as about:blank")
check("Starting ADapt" in appjs and "setInterval(tick" in appjs,
      "the wait is shown, in the tab and under the button, not a blank page")
check("api('/api/adapt/warm'" in appjs,
      "the page asks for ADapt to start as soon as it shows the button")

# --- ADapt through this server's own port -----------------------------------
# On a cluster the page is reached through one forwarded port. ADapt listened
# on a second one, and the tab said "refused to connect". Now /adapt/ on this server
# passes requests on. ADapt runs with --base-path /adapt and its page uses
# relative paths only, so the path is forwarded as it is and nothing in the
# answer is rewritten. A small HTTP server plays ADapt here the way its README
# describes: it takes the prefix off, redirects /adapt to /adapt/, checks its
# own token, and answers the requests the embedding uses.
import hashlib  # noqa: E402
import http.server  # noqa: E402
import json  # noqa: E402
import urllib.parse  # noqa: E402
import urllib.request  # noqa: E402
import urllib.error  # noqa: E402

FAKE_HTML = (b'<!doctype html><link rel="stylesheet" href="static/app.css">'
             b'<a href="docs/index.html">Docs</a><a href="//elsewhere.example/x">x</a>'
             b'<script src="static/3Dmol-min.js"></script><script src="static/app.js"></script>')
FAKE_JS = b"const u = 'api/status'; fetch(\"api/upload?name=a\"); const d = `docs/`;"
FAKE_CSS = b"body{background:url(bg.png)}"
seen_posts, seen_paths = [], []

# A session folder as ADapt writes one: the input, a receptor save with its
# PDBQT and _rec_H.pdb, and session.json recording both with their md5.
SESSION = os.path.join(tmp, "runs_proxy", "adapt_sessions", "1ycr_20261009-010203_abcdef")
SAVE = os.path.join(SESSION, "save_001_rec")
os.makedirs(SAVE)
PDBQT = ("REMARK 850 AGFR source 1ycr.pdb\nREMARK 850 AGFR biomolecule Biomol 1: chains A\n"
         "ATOM      1  N   GLU A  25      13.000  14.000  15.000  1.00  0.00    -0.3 N \n"
         "HETATM    2 ZN    ZN A 401      10.000  10.000  10.000  1.00  0.00     2.0 Zn\n")
PDBH = ("REMARK 850 AGFR biomolecule Biomol 1: chains A\n"
        "USER  MOD reduce.3.23.130521 H: found=0, std=0, add=1, rem=0, adj=1\n"
        "ATOM      1  N   GLU A  25      13.000  14.000  15.000  1.00  0.00           N\n"
        "ATOM      2  H1  GLU A  25      13.500  14.500  15.500  1.00  0.00           H\n"
        "HETATM    3 ZN    ZN A 401      10.000  10.000  10.000  1.00  0.00          ZN\n"
        "HETATM    4  O   HOH A 501      19.000  10.000  10.000  1.00  0.00           O\n")
open(os.path.join(SAVE, "1ycr_rec.pdbqt"), "w").write(PDBQT)
open(os.path.join(SAVE, "1ycr_rec_H.pdb"), "w").write(PDBH)
open(os.path.join(SESSION, "1ycr.pdb"), "w").write(PDBH)
md5 = lambda t: hashlib.md5(t.encode()).hexdigest()               # noqa: E731
SESSION_JSON = {
    "tool": {"name": "ADapt", "version": "0.1.0"}, "pdbid": "1ycr",
    "input": {"name": "1ycr.pdb"},
    "environment": {"adapt_commit": "1005e87a9b16", "adfr_dist": "adfr-0.0.22"},
    "events": [], "saves": [
        {"n": 1, "kind": "receptor", "tab": "Biomol 1", "at": "2026-10-09T01:02:03",
         "folder": "save_001_rec",
         "files": [{"name": "1ycr_rec.pdbqt", "md5": md5(PDBQT)},
                   {"name": "1ycr_rec_H.pdb", "md5": md5(PDBH)}],
         "mutations": {}, "altlocs": {}, "charges": {"A:ZN:401:_": 2.0},
         "reduce": {"version": "reduce.3.23.130521",
                    "decisions": [{"residue": "A:HIS:96:_", "decision": "no HD1"}]}},
        {"n": 2, "kind": "ligand", "tab": "Biomol 1", "folder": "save_002_lig",
         "files": [{"name": "1ycr_lig.pdbqt", "md5": "x"}]}]}
json.dump(SESSION_JSON, open(os.path.join(SESSION, "session.json"), "w"))


class FakeADapt(http.server.BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def _out(self, body, ctype, code=200, extra=()):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        for k, v in extra:
            self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def _path(self):
        seen_paths.append(self.path)
        path, _, query = self.path.partition("?")
        if path == "/adapt":
            # ADapt's own redirect: its page's paths need the slash.
            self._out(b"", "text/plain", 301,
                      extra=[("Location", "/adapt/" + ("?" + query if query else ""))])
            return None, None
        if path.startswith("/adapt/"):
            path = path[len("/adapt"):]
        return path, query

    def do_GET(self):
        path, query = self._path()
        if path is None:
            return
        if path == "/static/app.js":
            return self._out(FAKE_JS, "application/javascript")
        if path == "/static/app.css":
            return self._out(FAKE_CSS, "text/css")
        if "t=TOK" not in query:
            return self._out(b'{"error": "bad or missing token"}', "application/json", 403)
        if path == "/":
            return self._out(FAKE_HTML, "text/html; charset=utf-8")
        if path == "/api/session":
            if "s=s1" not in query:
                return self._out(b'{"error": "no such session; load the PDB file again"}',
                                 "application/json", 400)
            return self._out(json.dumps({"session": "s1", "folder": SESSION,
                                         "pdbid": "1ycr", "current": "Biomol 1"}).encode(),
                             "application/json")
        if path == "/api/sessions":
            return self._out(json.dumps({"sessions": [
                {"session": "s1", "pdbid": "1ycr", "folder_display": SESSION},
                {"session": "s9", "pdbid": "4hhb",
                 "folder_display": "~/elsewhere/4hhb_session"}]}).encode(),
                "application/json")
        if path == "/api/download":
            return self._out(b"ATOM", "chemical/x-pdbqt",
                             extra=[("Content-Disposition", 'attachment; filename="r.pdbqt"')])
        return self._out(b'{"ok": true}', "application/json")

    def do_POST(self):
        path, query = self._path()
        n = int(self.headers.get("Content-Length") or 0)
        seen_posts.append((self.path, self.rfile.read(n), self.headers.get("Content-Type")))
        if "t=TOK" not in query:
            return self._out(b'{"error": "bad or missing token"}', "application/json", 403)
        if path == "/api/fetch" and "code=9ZZZ" in query:
            return self._out(b'{"error": "9ZZZ has no PDB format file"}',
                             "application/json", 400)
        if path == "/api/save":
            # As ADapt's Session.save does: a new save folder with the
            # receptor PDBQT and _rec_H.pdb, recorded in session.json.
            body = json.loads(seen_posts[-1][1].decode())
            if body.get("s") != "s1":
                return self._out(b'{"error": "no such session; load the PDB file again"}',
                                 "application/json", 400)
            rec = json.load(open(os.path.join(SESSION, "session.json")))
            k = len(rec["saves"]) + 1
            sub = "save_%03d_rec" % k
            os.makedirs(os.path.join(SESSION, sub))
            files = []
            for name, text in (("1ycr_rec.pdbqt", PDBQT), ("1ycr_rec_H.pdb", PDBH)):
                open(os.path.join(SESSION, sub, name), "w").write(text)
                files.append({"name": name, "md5": md5(text)})
            rec["saves"].append({"n": k, "kind": "receptor", "tab": "Biomol 1",
                                 "at": "2026-10-09T02:00:%02d" % k, "folder": sub,
                                 "files": files, "mutations": {}, "altlocs": {},
                                 "charges": {}, "reduce": {"version": "reduce.3.23.130521",
                                                           "decisions": []}})
            json.dump(rec, open(os.path.join(SESSION, "session.json"), "w"))
            return self._out(json.dumps({"save": k, "files": files, "log": []}).encode(),
                             "application/json")
        return self._out(json.dumps({"session": "s1", "folder": SESSION, "pdbid": "1ycr",
                                     "current": "Biomol 1"}).encode(), "application/json")


fake = http.server.ThreadingHTTPServer(("127.0.0.1", 0), FakeADapt)
threading.Thread(target=fake.serve_forever, daemon=True).start()


class StubRunner:
    url = "http://127.0.0.1:%d/adapt/?t=TOK" % fake.server_address[1]

    def running(self):
        return True

    def stop(self):
        pass


proxy_app = server.App(os.path.join(tmp, "runs_proxy"))
httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
httpd.app = proxy_app
httpd.token = proxy_app.token
threading.Thread(target=httpd.serve_forever, daemon=True).start()
BASE = "http://127.0.0.1:%d" % httpd.server_address[1]
# A cluster that sets a web proxy must not have 127.0.0.1 sent to it.
os.environ["http_proxy"] = os.environ["HTTP_PROXY"] = "http://127.0.0.1:9"


class _NoFollow(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, *a, **kw):
        return None


def get(path, data=None, headers=None):
    req = urllib.request.Request(BASE + path, data=data, headers=headers or {},
                                 method="POST" if data is not None else "GET")
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}), _NoFollow)
    try:
        with opener.open(req, timeout=20) as r:
            return r.status, r.read(), r.headers
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read(), exc.headers


code, body, _h = get("/adapt/?t=TOK")
check(code == 503 and b"Open in ADapt" in body,
      "before ADapt runs, /adapt/ says to press Open in ADapt (%s %r)" % (code, body[:80]))
proxy_app._adapt_runner = StubRunner()
code, body, h = get("/adapt/?t=TOK")
check(code == 200 and body == FAKE_HTML,
      "ADapt's page arrives exactly as ADapt sent it, nothing rewritten (%r)" % body[:160])
check(seen_paths and seen_paths[-1] == "/adapt/?t=TOK",
      "the request reaches ADapt with its /adapt prefix, for --base-path to take off (%s)"
      % seen_paths[-1:])
code, body, h = get("/adapt?t=TOK")
check(code == 301 and h.get("Location") == "/adapt/?t=TOK",
      "ADapt's redirect to /adapt/ goes back to the browser, not followed here (%s %s)"
      % (code, h.get("Location")))
code, body, _h = get("/adapt/static/app.js")
check(body == FAKE_JS, "its script passes through untouched (%r)" % body)
code, body, _h = get("/adapt/static/app.css")
check(body == FAKE_CSS, "and its stylesheet")
code, body, _h = get("/adapt/api/status?t=TOK")
check(code == 200 and b'"ok"' in body, "ADapt's API answers through this server")
code, body, _h = get("/adapt/api/status")
check(code == 403, "and ADapt still refuses a request without its token (%s)" % code)
code, body, _h = get("/adapt/api/upload?t=TOK&name=a.pdb", data=b"ATOM      1  N",
                     headers={"Content-Type": "application/octet-stream"})
check(code == 200 and seen_posts and seen_posts[-1][1] == b"ATOM      1  N"
      and seen_posts[-1][0].startswith("/adapt/api/upload?"),
      "an uploaded file reaches ADapt intact, under the prefix (%r)" % (seen_posts[-1:],))
code, body, h = get("/adapt/api/download?t=TOK&s=1")
check(h.get("Content-Disposition", "").startswith("attachment"),
      "a PDBQT download keeps its file name")
code, body, _h = get("/api/tools")
check(code == 403, "this tool's own API still needs its own token (%s)" % code)
check(not hasattr(server, "_rebase_adapt"),
      "the path rewriting that --base-path made unnecessary is gone")

# --- the receptor goes to ADapt through its own api/upload --------------------
proxy_app._adapt_running = lambda: StubRunner.url
rec = os.path.join(tmp, "1YCR.pdb")
open(rec, "w").write(PDBH)
seen_posts.clear()
res = proxy_app.adapt_open({"path": rec})
check(res.get("session") == "s1" and res.get("folder") == SESSION,
      "the receptor is handed to ADapt and its session named (%s)" % res)
check(res.get("url") == "/adapt/?t=TOK&session=s1",
      "the tab opens ADapt on that session, under this server (%s)" % res.get("url"))
check(seen_posts and seen_posts[-1][0].startswith("/adapt/api/upload?")
      and "name=1ycr.pdb" in seen_posts[-1][0] and seen_posts[-1][1] == PDBH.encode(),
      "the file goes whole, named in lower case by its id, as ADFR looks ligands up by it (%s)"
      % (seen_posts[-1][0] if seen_posts else None))
res = proxy_app.adapt_open({"path": rec, "session": "s1"})
check(res.get("url") == "/adapt/?t=TOK&session=s1" and len(seen_posts) == 1,
      "a session ADapt still holds is opened again rather than uploaded anew")
cif = os.path.join(tmp, "7a3r.cif")
open(cif, "w").write("data_7A3R\n")
seen_posts.clear()
res = proxy_app.adapt_open({"path": cif, "code": "7A3R"})
check(res.get("url") == "/adapt/?t=TOK&session=s1" and seen_posts
      and seen_posts[-1][0].startswith("/adapt/api/fetch?") and "code=7A3R" in seen_posts[-1][0],
      "an mmCIF fetched by its id is fetched by ADapt itself, so its session is "
      "known here (%s, %s)" % (res.get("url"), seen_posts[-1:]))
res = proxy_app.adapt_open({"path": cif, "code": "9ZZZ"})
check(res.get("url") == "/adapt/?t=TOK&pdb=9ZZZ" and "could not fetch" in res.get("note", ""),
      "an id ADapt cannot fetch here sends the page to ADapt by the id, with why (%s)"
      % res)
res = proxy_app.adapt_open({"path": cif})
check("PDB files only" in res.get("note", "") and "session=" not in res.get("url", ""),
      "an mmCIF without an id is not uploaded, and the page is told why (%s)" % res)

# --- what ADapt saved, read back ------------------------------------------------
got = adapt.read_save(SESSION)
check(got["save"] == 1 and got["tab"] == "Biomol 1" and got["pdbid"] == "1ycr",
      "the newest receptor save is read, not the ligand one (%s)" % {k: got[k] for k in ("save", "tab")})
check(got["pdbqt"].endswith("1ycr_rec.pdbqt") and got["pdb"].endswith("1ycr_rec_H.pdb"),
      "the receptor PDBQT and the _rec_H.pdb are found")
check(got["hetero"] == ["ZN A:401"], "the hetero groups kept are listed, waters left out (%s)"
      % got["hetero"])
check(got["biomolecule"] == "Biomol 1: chains A", "the biomolecule line is read (%r)" % got["biomolecule"])
check(got["charges"] == {"A:ZN:401:_": 2.0} and got["reduce"] == "reduce.3.23.130521"
      and got["his"] == {"A:HIS:96:_": "no HD1"},
      "the charges, reduce's version and its histidine decisions are read from session.json")
check(got["adapt_commit"] == "1005e87", "and the ADapt that wrote it (%r)" % got["adapt_commit"])
try:
    adapt.read_save(SESSION, 2)
    check(False, "a ligand save is taken as a receptor")
except adapt.NoSave as exc:
    check("no receptor save 2" in str(exc), "a save that is not a receptor is refused (%s)" % exc)
try:
    adapt.read_save(tmp)
    check(False, "a folder without session.json is read")
except adapt.NoSave as exc:
    check("not an ADapt session folder" in str(exc), "a folder that is not a session says so")
r = proxy_app.adapt_saved(session="s1")
check(r.get("found") and r.get("folder") == SESSION and r.get("session") == "s1",
      "the page reads the save back through the session ADapt names (%s)" % r.get("why"))
open(os.path.join(SAVE, "1ycr_rec.pdbqt"), "a").write("REMARK changed\n")
r = proxy_app.adapt_saved(folder=SESSION)
check(not r.get("found") and "changed since ADapt wrote it" in r.get("why", ""),
      "a file changed since ADapt wrote it is refused (%s)" % r.get("why"))
open(os.path.join(SAVE, "1ycr_rec.pdbqt"), "w").write(PDBQT)

# A run set up from the save records it, and its receptor is the _rec_H.pdb.
from adcp_pilot import job as jobmod  # noqa: E402
j = jobmod.default(adapt={"folder": SESSION, "session": "s1"}, protonation="propka")
check(jobmod.take_adapt(j) == [] and j["receptor"] == got["pdb"]
      and j["adapt"]["save"] == 1 and j["adapt"]["pdbqt_md5"] == md5(PDBQT)
      and j["protonation"] == "reduce" and j["adapt"]["session"] == "s1",
      "a run from ADapt records the save and its md5, and is protonated as ADapt did (%s)"
      % j.get("adapt"))
check(jobmod.adapt_problems(j) == [], "an unchanged save is fine")
open(os.path.join(SAVE, "1ycr_rec_H.pdb"), "a").write("REMARK changed\n")
check(jobmod.adapt_problems(j), "a save changed since the run was set up is refused")
open(os.path.join(SAVE, "1ycr_rec_H.pdb"), "w").write(PDBH)
# --- what has been changed in a session, read from its session.json ------------
# ADapt rewrites session.json after every edit and every save. Reading the file
# leaves ADapt's running session alone: its api/session hands out the log lines
# its own page has not shown yet.
st = adapt.read_state(SESSION)
check(st["edits"] == 0 and st["tab"] == "" and [x["n"] for x in st["saves"]] == [1],
      "a session nothing was changed in has no edits, and its receptor save is listed "
      "(%s)" % {k: st[k] for k in ("edits", "tab", "saves")})
rec_ = json.load(open(os.path.join(SESSION, "session.json")))
rec_["opened_tab"] = "Biomol 1"
rec_["events"] = [
    {"what": "check", "node": "B", "key": "B", "column": "rec", "state": "unchecked",
     "tab": "Biomol 1", "at": "2026-10-09T01:00:00"},
    {"what": "check", "node": "water (55)", "column": "rec", "state": "checked",
     "tab": "Biomol 1", "at": "2026-10-09T01:01:00"},
    {"what": "charge", "node": "ZN 401", "charge": 2.0, "tab": "Biomol 1",
     "at": "2026-10-09T03:00:00"}]
json.dump(rec_, open(os.path.join(SESSION, "session.json"), "w"))
st = adapt.read_state(SESSION)
check(st["edits"] == 3 and st["tab"] == "Biomol 1",
      "every edit ADapt recorded is counted (%s)" % st["edits"])
check(st["changes"] == ["chain B out of the receptor", "water (55) in the receptor",
                        "ZN 401 charge 2.0"],
      "and said in words (%s)" % st["changes"])
check(st["saves"][0]["edits_before"] == 2,
      "a save knows how many edits came before it, by ADapt's own times (%s)" % st["saves"])
try:
    adapt.read_state(tmp)
    check(False, "a folder without session.json has a state")
except adapt.NoSave:
    check(True, "a folder without session.json has none")

# Through the server, by the session alone: its folder is known from the open.
proxy_app._adapt_folders()["s1"] = SESSION
got_st = proxy_app.adapt_state(session="s1")
check(got_st.get("found") and got_st.get("edits") == 3 and got_st.get("folder") == SESSION,
      "the page is told what has been changed in its session (%s)" % got_st.get("why"))
proxy_app._adapt_folders().clear()
got_st = proxy_app.adapt_state(session="s1")
check(got_st.get("found") and got_st.get("folder") == SESSION,
      "a session opened inside ADapt is found in ADapt's own list (%s)" % got_st.get("why"))
check(adapt.session_folder(StubRunner.url, "s9")
      == os.path.expanduser("~/elsewhere/4hhb_session"),
      "a folder ADapt shows from the home folder is read as that folder")
check(not proxy_app.adapt_state(session="gone").get("found"),
      "a session ADapt does not hold is said to be gone")
code, body, _h = get("/api/adapt/state?session=s1&t=" + proxy_app.token)
check(code == 200 and json.loads(body).get("edits") == 3,
      "and asked for over HTTP with this page's token (%s)" % code)
code, body, _h = get("/api/adapt/state?session=s1")
check(code == 403, "not without it (%s)" % code)

# --- ADapt writes the receptor when asked, the way its own button does ------------
seen_posts.clear()
w1 = proxy_app.adapt_write({"session": "s1", "folder": SESSION})
check(w1.get("found") and w1.get("written") and w1.get("save") == 3,
      "ADapt is asked to write the receptor and its save is read back (%s)"
      % {k: w1.get(k) for k in ("found", "written", "save", "why")})
saves_posted = [x for x in seen_posts if x[0].startswith("/adapt/api/save?")]
check(len(saves_posted) == 1
      and json.loads(saves_posted[0][1].decode()) == {"s": "s1", "kind": "receptor"},
      "with the call ADapt's receptor PDBQT button makes: api/save, the session, "
      "kind receptor (%s)" % saves_posted)
check(w1.get("pdbqt", "").endswith("save_003_rec/1ycr_rec.pdbqt")
      and w1.get("pdbqt_md5") == md5(PDBQT) and w1.get("pdb_md5") == md5(PDBH),
      "read back as any save is, with the checksums the run records")
check(w1.get("edits") == 3, "and with the number of edits it holds")
st3 = proxy_app.adapt_state(session="s1", folder=SESSION)
mine3 = [x for x in st3["saves"] if x["n"] == 3]
check(mine3 and mine3[0].get("written_here") and mine3[0]["edits_before"] == 3,
      "the save this server asked for is marked, with the exact edits before it (%s)"
      % mine3)
w2 = proxy_app.adapt_write({"session": "s1", "folder": SESSION})
check(w2.get("found") and not w2.get("written") and w2.get("save") == 3
      and len([x for x in seen_posts if x[0].startswith("/adapt/api/save?")]) == 1,
      "nothing changed since, so the same save is used rather than another written")
rec_ = json.load(open(os.path.join(SESSION, "session.json")))
rec_["events"].append({"what": "altloc", "node": "SER 5", "altloc": "B",
                       "tab": "Biomol 1", "at": "2026-10-09T04:00:00"})
json.dump(rec_, open(os.path.join(SESSION, "session.json"), "w"))
w3 = proxy_app.adapt_write({"session": "s1", "folder": SESSION})
check(w3.get("found") and w3.get("written") and w3.get("save") == 4,
      "a change since the last save is written anew (%s)" % w3.get("save"))
j3 = jobmod.default(adapt={"folder": w3["folder"], "save": w3["save"], "session": "s1"})
check(jobmod.take_adapt(j3) == [] and j3["adapt"]["save"] == 4
      and j3["adapt"]["pdbqt_md5"] == md5(PDBQT) and j3["receptor"] == w3["pdb"],
      "and a run set up from it records that save, as job.json's adapt key")
check(not proxy_app.adapt_write({"session": "nope"}).get("found"),
      "a session whose folder is not known is not written")
rec_ = json.load(open(os.path.join(SESSION, "session.json")))
rec_["events"].append({"what": "tab", "name": "Biomol 2", "tab": "Biomol 2",
                       "at": "2026-10-09T05:00:00"})
json.dump(rec_, open(os.path.join(SESSION, "session.json"), "w"))
bad_w = proxy_app.adapt_write({"session": "nope", "folder": SESSION})
check(not bad_w.get("found") and "no such session" in bad_w.get("why", ""),
      "a session ADapt no longer holds is refused with ADapt's reason (%s)" % bad_w)
check(not proxy_app.adapt_write({}).get("found"), "no session, nothing written")
code, body, _h = get("/api/adapt/write?t=" + proxy_app.token,
                     data=json.dumps({"session": "s1", "folder": SESSION}).encode(),
                     headers={"Content-Type": "application/json"})
check(code == 200 and json.loads(body).get("save") == 5,
      "and asked for over HTTP (%s %r)" % (code, body[:120]))

# --- ADapt may be shown in a frame: nothing it sends forbids it ---------------------
code, body, h = get("/adapt/?t=TOK")
check(not h.get("X-Frame-Options") and "frame-ancestors" not in (h.get(
      "Content-Security-Policy") or ""),
      "ADapt's page through this server carries nothing that forbids a frame")

httpd.shutdown()
fake.shutdown()
for k in ("http_proxy", "HTTP_PROXY"):
    os.environ.pop(k, None)

# --- the interface's terminal hangs up: ADapt goes with it --------------------------
# On a cluster the interface runs on the terminal of an SSH session, and the session ending hangs
# that terminal up. SIGHUP ended the interface where it stood, before the clean-up that stops
# ADapt, and ADapt, in a session of its own, went on running on the login node: one was found
# there two days later. The interface is run here as `serve` runs it, ADapt (the stand-in) is
# started as the page starts it, and the interface is hung up on.
import signal  # noqa: E402
import subprocess  # noqa: E402

hang = os.path.join(tmp, "hangup.py")
open(hang, "w").write("""
import os, sys, threading
sys.path.insert(0, %(root)r)
from adcp_pilot import adapt, adcp11, server
adcp11._HERE = adapt._HERE = %(alone)r
built = server.build
def build(*a, **k):
    h = built(*a, **k)
    def warm():
        h.app.adapt_open()
        print("ADAPT-PID %%d" %% h.app._adapt_runner.proc.pid, flush=True)
    threading.Timer(0.3, warm).start()
    return h
server.build = build
server.serve(host="127.0.0.1", port=0, runs=%(runs)r)
""" % {"root": ROOT, "alone": adapt._HERE, "runs": os.path.join(tmp, "runs_hang")})
env_ = dict(os.environ, ADAPT_HOME=home, ADCPSUITE_HOME=suite, HOME=tmp)
iface = subprocess.Popen([sys.executable, "-u", hang], stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, env=env_, universal_newlines=True)
adapt_pid, said_ = None, []
end_ = time.time() + 60
while time.time() < end_ and adapt_pid is None:
    line = iface.stdout.readline()
    if not line:
        break
    said_.append(line.rstrip())
    if line.startswith("ADAPT-PID "):
        adapt_pid = int(line.split()[1])
check(adapt_pid is not None, "the interface started ADapt (%s)" % " / ".join(said_[-4:]))


def _running(pid):
    try:
        with open("/proc/%d/stat" % pid) as fh:
            return fh.read().rsplit(")", 1)[1].split()[0] != "Z"
    except OSError:
        return False


if adapt_pid:
    check(_running(adapt_pid), "and ADapt is running")
    iface.send_signal(signal.SIGHUP)
    try:
        iface.wait(timeout=20)
    except subprocess.TimeoutExpired:
        iface.kill()
    for _ in range(100):
        if not _running(adapt_pid):
            break
        time.sleep(0.1)
    check(iface.returncode is not None, "a hang-up ends the interface")
    check(not _running(adapt_pid), "and ADapt with it, rather than leaving it on the login node")
    if _running(adapt_pid):
        os.killpg(adapt_pid, signal.SIGKILL)
else:
    iface.kill()

if fails:
    print("test_adapt FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  ADapt is found, started with the ADCP 1.1 environment, reused and stopped, %d checks" % ran)
