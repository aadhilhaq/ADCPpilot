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

fails, ran = [], 0


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
check(res.get("url") == "http://127.0.0.1:45678/?t=abc123", "and opening it returns its address (%s)" % res)
env_root = open(os.path.join(tmp, "runs", "adapt_sessions", "started_with.txt")).read().split("|")[0]
check(os.path.realpath(env_root) == os.path.realpath(suite),
      "the interface hands ADapt the ADCP 1.1 environment it found (%s)" % env_root)
check(res.get("port") == 45678 and res.get("host"),
      "with the machine and the port, so the page can say what to forward (%s)" % res)
check(res.get("pinned") is False and res.get("port_env") == adapt.PORT_ENV,
      "and whether that port will be the same next time")
check(not any(os.path.basename(x["rundir"]) == "adapt_sessions" for x in app.list_runs()["runs"]),
      "ADapt's sessions folder is not listed as a run")
app.adapt_stop()
check(not app._adapt_runner.running(), "the interface stops ADapt when it stops")

# --- the page offers it only when the server does --------------------------
from adcp_dock import page  # noqa: E402

html = page.PAGE if hasattr(page, "PAGE") else open(page.__file__).read()
check("id=adaptrow" in html and "id=adaptopen" in html, "the Receptor card has an ADapt row and button")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("d.adapt && d.adapt.available" in appjs, "the row is shown only when the server offers ADapt")
check("api('/api/adapt'" in appjs, "the button asks the server to start it")
check("window.open('', '_blank')" in appjs,
      "the tab is opened inside the click, where a browser allows it")
check("ssh -L ' + esc(port) + ':127.0.0.1:' + esc(port)" in appjs,
      "the note gives the forwarding command with the real port")
check("esc(host)" in appjs and "r.host" in appjs,
      "and names the machine, rather than leaving the user to guess it")
check("r.port_env" in appjs, "and says how to keep the same port next time")

if fails:
    print("test_adapt FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  ADapt is found, started with the ADCP 1.1 environment, reused and stopped, %d checks" % ran)
