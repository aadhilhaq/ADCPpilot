"""ADapt is found, started with the right environment, reused and stopped.

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

from adcp_dock import adapt  # noqa: E402

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


def fake_adapt(folder, speaks=True):
    """A folder that looks like an ADapt checkout, whose bin/adapt records
    how it was started and then prints an address, or dies without one."""
    os.makedirs(os.path.join(folder, "bin"), exist_ok=True)
    os.makedirs(os.path.join(folder, "adapt"), exist_ok=True)
    open(os.path.join(folder, "adapt", "__main__.py"), "w").write("")
    script = ["#!/bin/sh",
              'mkdir -p "$ADAPT_OUTDIR"',
              'echo "$ADCPSUITE|$ADAPT_OUTDIR|$PYTHONUNBUFFERED|$*" > "$ADAPT_OUTDIR/started_with.txt"']
    if speaks:
        script += ['echo "ADapt 0.1.0 is at"', 'echo ""',
                   'echo "    http://127.0.0.1:45678/?t=abc123"', "exec sleep 60"]
    else:
        script += ['echo "no ADCP environment found" >&2', "exit 2"]
    path = os.path.join(folder, "bin", "adapt")
    open(path, "w").write("\n".join(script) + "\n")
    os.chmod(path, 0o755)
    return folder


tmp = tempfile.mkdtemp(prefix="adapt_test_")
os.environ["HOME"] = tmp                  # nothing may fall back to a real home
for k in ("ADAPT_HOME", "ADFRSUITE_HOME", "ADCP_HOME", "ADCPSUITE_HOME"):
    os.environ.pop(k, None)

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
check(url == "http://127.0.0.1:45678/?t=abc123",
      "the address is read from ADapt's output (%r)" % url)
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
from adcp_dock import server  # noqa: E402

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
check(res.get("direct") == "http://127.0.0.1:45678/?t=abc123",
      "and ADapt's own address beside it")
env_root = open(os.path.join(tmp, "runs", "adapt_sessions", "started_with.txt")).read().split("|")[0]
check(os.path.realpath(env_root) == os.path.realpath(suite),
      "the interface hands ADapt the ADCP 1.1 environment it found (%s)" % env_root)
check(not any(os.path.basename(x["rundir"]) == "adapt_sessions" for x in app.list_runs()["runs"]),
      "ADapt's sessions folder is not listed as a run")
app.adapt_stop()
check(not app._adapt_runner.running(), "the interface stops ADapt when it stops")

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
from adcp_dock import page  # noqa: E402

html = page.PAGE if hasattr(page, "PAGE") else open(page.__file__).read()
check("id=adaptrow" in html and "id=adaptopen" in html, "the Receptor card has an ADapt row and button")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("d.adapt && d.adapt.available" in appjs, "the row is shown only when the server offers ADapt")
check("api('/api/adapt'" in appjs, "the button asks the server to start it")
check("window.open('', '_blank')" in appjs,
      "the tab is opened inside the click, where a browser allows it")
check("ssh -L" not in appjs,
      "nobody is asked to forward a second port: ADapt comes through this page's own")
check("new URL(r.url, location.href)" in appjs,
      "the tab is sent to an absolute address, since it began as about:blank")
check("Starting ADapt" in appjs and "setInterval(tick" in appjs,
      "the wait is shown, in the tab and under the button, not a blank page")
check("api('/api/adapt/warm'" in appjs,
      "the page asks for ADapt to start as soon as it shows the button")

# --- ADapt through this server's own port -----------------------------------
# On Grace the page is reached through one forwarded port. ADapt listened on a
# second one, and the tab said "refused to connect". Now /adapt/ on this server
# passes requests on. A small HTTP server plays ADapt here: its page and script
# use absolute /static/, /docs/ and /api/ paths, as ADapt's do, and its API
# checks its own token.
import http.server  # noqa: E402
import threading  # noqa: E402
import urllib.request  # noqa: E402
import urllib.error  # noqa: E402

FAKE_HTML = (b'<!doctype html><link rel="stylesheet" href="/static/app.css">'
             b'<a href="/docs/index.html">Docs</a><a href="//elsewhere.example/x">x</a>'
             b'<script src="/static/3Dmol-min.js"></script><script src="/static/app.js"></script>')
FAKE_JS = b"const u = '/api/status'; fetch(\"/api/upload?name=a\"); const d = `/docs/`;"
FAKE_3DMOL = b'var x = "/api/keep-me-as-i-am";'
FAKE_CSS = b"body{background:url(/static/bg.png)} a{background:url(//cdn.example/y.png)}"
seen_posts = []


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

    def do_GET(self):
        path, _, query = self.path.partition("?")
        if path == "/static/app.js":
            return self._out(FAKE_JS, "application/javascript")
        if path == "/static/3Dmol-min.js":
            return self._out(FAKE_3DMOL, "application/javascript")
        if path == "/static/app.css":
            return self._out(FAKE_CSS, "text/css")
        if "t=TOK" not in query:
            return self._out(b'{"error": "bad or missing token"}', "application/json", 403)
        if path == "/":
            return self._out(FAKE_HTML, "text/html; charset=utf-8")
        if path == "/api/download":
            return self._out(b"ATOM", "chemical/x-pdbqt",
                             extra=[("Content-Disposition", 'attachment; filename="r.pdbqt"')])
        return self._out(b'{"ok": true}', "application/json")

    def do_POST(self):
        n = int(self.headers.get("Content-Length") or 0)
        seen_posts.append((self.path, self.rfile.read(n), self.headers.get("Content-Type")))
        if "t=TOK" not in self.path:
            return self._out(b'{"error": "bad or missing token"}', "application/json", 403)
        return self._out(b'{"session": "s1"}', "application/json")


fake = http.server.ThreadingHTTPServer(("127.0.0.1", 0), FakeADapt)
threading.Thread(target=fake.serve_forever, daemon=True).start()


class StubRunner:
    url = "http://127.0.0.1:%d/?t=TOK" % fake.server_address[1]

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


def get(path, data=None, headers=None):
    req = urllib.request.Request(BASE + path, data=data, headers=headers or {},
                                 method="POST" if data is not None else "GET")
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    try:
        with opener.open(req, timeout=20) as r:
            return r.status, r.read(), r.headers
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read(), exc.headers


code, body, _h = get("/adapt/?t=TOK")
check(code == 503 and b"Open ADapt" in body,
      "before ADapt runs, /adapt/ says to press Open ADapt (%s %r)" % (code, body[:80]))
proxy_app._adapt_runner = StubRunner()
code, body, h = get("/adapt/?t=TOK")
check(code == 200 and b'href="/adapt/static/app.css"' in body and b'src="/adapt/static/app.js"' in body,
      "ADapt's page arrives with its paths moved under /adapt/ (%r)" % body[:160])
check(b'href="/adapt/docs/index.html"' in body, "its Docs link too")
check(b'href="//elsewhere.example/x"' in body, "and a link to another host is left alone")
code, body, _h = get("/adapt/static/app.js")
check(b"'/adapt/api/status'" in body and b'"/adapt/api/upload' in body and b"`/adapt/docs/`" in body,
      "its script's API paths are moved too, whatever the quote (%r)" % body)
code, body, _h = get("/adapt/static/3Dmol-min.js")
check(body == FAKE_3DMOL, "the 3Dmol library is passed through untouched")
code, body, _h = get("/adapt/static/app.css")
check(b"url(/adapt/static/bg.png)" in body and b"url(//cdn.example/y.png)" in body,
      "stylesheet paths are moved, other hosts are not (%r)" % body)
code, body, _h = get("/adapt/api/status?t=TOK")
check(code == 200 and b'"ok"' in body, "ADapt's API answers through this server")
code, body, _h = get("/adapt/api/status")
check(code == 403, "and ADapt still refuses a request without its token (%s)" % code)
code, body, _h = get("/adapt/api/upload?t=TOK&name=a.pdb", data=b"ATOM      1  N",
                     headers={"Content-Type": "application/octet-stream"})
check(code == 200 and seen_posts and seen_posts[-1][1] == b"ATOM      1  N"
      and seen_posts[-1][0].startswith("/api/upload?"),
      "an uploaded file reaches ADapt intact, at ADapt's own path (%r)" % (seen_posts[-1:],))
code, body, h = get("/adapt/api/download?t=TOK&s=1")
check(h.get("Content-Disposition", "").startswith("attachment"),
      "a PDBQT download keeps its file name")
code, body, _h = get("/api/tools")
check(code == 403, "this tool's own API still needs its own token (%s)" % code)
httpd.shutdown()
fake.shutdown()
for k in ("http_proxy", "HTTP_PROXY"):
    os.environ.pop(k, None)

if fails:
    print("test_adapt FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  ADapt is found, started with the ADCP 1.1 environment, reused and stopped, %d checks" % ran)
