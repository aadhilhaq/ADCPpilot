"""The launcher, without a cluster, a display or paramiko.

By Aadhil Haq

What can be checked here is the shape: the command it would run, the line
it would recognize coming back, the timeout it allows for two-factor
approval, and that the authorship is written in one place that everything
else reads. Whether a cluster accepts the key needs a real cluster.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import credit as tool_credit, page               # noqa: E402
from adcp_launcher import credit as launcher_credit, entry      # noqa: E402
from adcp_launcher import session, tunnel                       # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


# --- the command it runs on the far end -----------------------------------
line = session.command("/scratch/x/ADCP_docking", runs="~/adcp_runs",
                       tools_roots=["/opt/adcpsuite", "/opt/amber"])
check("adcp_dock.cli serve" in line, "the launcher does not start the server")
check("--host 127.0.0.1" in line,
      "the server would be started on a public interface, not loopback")
check("--port 0" in line,
      "a fixed port is asked for; two people on the same login node would "
      "then collide")
check(line.count("--tools-root") == 2, "the tool roots are not passed through")
check(line.startswith("cd "),
      "the command does not change into the install directory first")
check("PYTHONPATH=" in line,
      "PYTHONPATH is not set, so the module would not be importable")

# --- a tilde path still expands on the far end -----------------------------
#
# shlex.quote puts a path in single quotes, which is correct for a space but
# stops the remote shell expanding a tilde, so cd fails on a directory that
# exists.
tilde = session.command("~/ADCP_docking", runs="~/adcp_runs")
check("'~/ADCP_docking'" not in tilde,
      "the tool directory is still single-quoted with the tilde inside, so "
      "the remote shell cannot expand it and cd fails on a directory that "
      "exists: %s" % tilde[:120])
check("'~/adcp_runs'" not in tilde,
      "the runs directory has the same problem: %s" % tilde[:120])
check('"$HOME"/ADCP_docking' in tilde,
      "the tilde is not being turned into something the far end expands: %s"
      % tilde[:120])
# the spelling of the separator is not the point: that it changes into the
# checkout FIRST, and stops if it cannot, is.
check(tilde.startswith('cd "$HOME"/ADCP_docking'),
      "the command does not start by changing into the checkout: %s"
      % tilde[:80])
check("|| exit" in tilde.split(";")[0],
      "the command carries on after a failed cd, so the tool would be run "
      "from whatever directory the login shell left it in: %s" % tilde[:80])

# and a path with a space is still quoted, which is why the tilde case could
# not simply be left unquoted
spacey = session.command("~/my tools", runs="/scratch/a b/runs")
check("'my tools'" in spacey,
      "a space in a tilde path is no longer quoted, so the command would "
      "split: %s" % spacey[:120])
check("'/scratch/a b/runs'" in spacey,
      "a space in an absolute path is no longer quoted: %s" % spacey[:120])

# an absolute path is untouched: no $HOME, no needless quoting games
absolute = session.command("/shared/apps/ADCP_docking")
check("$HOME" not in absolute,
      "an absolute path was rewritten through $HOME: %s" % absolute[:120])
check(absolute.startswith("cd /shared/apps/ADCP_docking"),
      "an absolute path is not used as given: %s" % absolute[:120])

# --- the cluster's own python, not whatever python3 happens to be ----------
#
# A cluster's system python3 may be too old to parse this package (for
# example 3.6, which fails on `from __future__ import annotations`). The
# command sources activate.sh and checks the version first, so an old Python
# produces a sentence rather than a traceback.
line2 = session.command("/shared/apps/ADCP_docking")
check("activate.sh" in line2,
      "the command does not source activate.sh, so it uses whatever python3 "
      "is on the login node's PATH, which may be too old to parse "
      "this code: %s" % line2[:140])
check("ADCP_PYTHON" in line2,
      "nothing names the interpreter the installer chose")
check("version_info" in line2,
      "the Python version is not checked before the import that would fail, "
      "so the failure arrives as a SyntaxError traceback")
check(session.OLD_PYTHON_MARKER in line2,
      "there is no marker for the launcher to recognize an old Python by")
# Order matters, but only check it when both parts are present: .index()
# raises when a part is missing, and a traceback is not a failure message.
def _before(text, first, second, why):
    if first in text and second in text:
        check(text.index(first) < text.index(second), why)


_before(line2, "activate.sh", "adcp_dock.cli",
        "activate.sh is sourced after the tool is started, too late to "
        "choose the interpreter that runs it")
_before(line2, "version_info", "adcp_dock.cli",
        "the Python version is checked after the tool is started, so the "
        "SyntaxError happens first anyway")
# The check program must be parseable by the old Python it is checking for.
# It is extracted defensively, so a command with no -c at all is reported
# rather than crashing this file.
_prog = ""
if "-c " in line2:
    _prog = line2.split("-c ", 1)[1].split("' ", 1)[0].strip("'")
check(_prog, "there is no inline version-check program in the command")
if _prog:
    import ast as _ast                                          # noqa: E402
    try:
        _ast.parse(_prog)
        _parses = True
    except SyntaxError:
        _parses = False
    check(_parses, "the version check is not valid Python: %r" % _prog)
    check("f'" not in _prog and ':=' not in _prog and "match " not in _prog,
          "the version check uses syntax an old Python cannot parse, so it "
          "would fail with the very error it exists to replace: %r" % _prog)

# and a named python still wins, for anyone this does not describe
check("/my/python" in session.command("/x", python="/my/python"),
      "an explicitly named interpreter is ignored")

# the launcher turns that marker into something actionable
_msg = session._why(session.OLD_PYTHON_MARKER + " Python 3.6.8", "/x/y")
check("3.6.8" in _msg and "install.sh" in _msg,
      "an old Python is not explained in terms of what to do: %r" % _msg)
check("Traceback" not in _msg,
      "the explanation is still a traceback")
_msg2 = session._why("SyntaxError: future feature annotations is not "
                     "defined", "/x/y")
check("install.sh" in _msg2,
      "the raw SyntaxError is not recognized: %r" % _msg2)

# --- the line it reads back --------------------------------------------------
m = session.URL.search("ADCP docking is at\n\n"
                       "    http://127.0.0.1:40515/?t=yDCQJyk0-_abc\n")
check(m is not None, "the server's own URL line is not recognized")
if m:
    check(m.group(1) == "40515" and m.group(2) == "yDCQJyk0-_abc",
          "the port and token were read as %r / %r" % (m.group(1), m.group(2)))
check(session.URL.search("http://10.0.0.5:40515/?t=abc") is None,
      "a URL on some other host would be accepted as the tunnel's own")

# --- what the server prints is what the launcher looks for ---------------
src = open(os.path.join(ROOT, "adcp_dock", "server.py")).read()
check('"http://%s:%d/?t=%s"' in src,
      "the server's URL format changed; the launcher's pattern would miss it")

# --- two-factor challenges ---------------------------------------------------
check(tunnel.AUTH_PATIENCE >= 120,
      "the two-factor timeout is %ds. paramiko's 30 has to cover reading the menu, "
      "choosing, then finding a phone and approving; 120 is OpenSSH's own "
      "LoginGraceTime" % tunnel.AUTH_PATIENCE)
h = tunnel._handler(lambda *a: ["should not be asked"], log=lambda *a: None)
check(h("", "The server says hello", []) == [],
      "a keyboard-interactive request with NO prompts was answered with a "
      "dialog. RFC 4256 3.3 says that is a message to display and the client "
      "must not prompt, and the seconds spent on it come out of the same "
      "budget as the two-factor approval")
asked = []
h2 = tunnel._handler(lambda t, i, p: asked.append(p) or ["x"] * len(p))
check(h2("", "", [("Passcode: ", False)]) == ["x"],
      "a real prompt was not passed to the person")
check(asked and asked[0] == [("Passcode: ", False)],
      "the prompt reached the dialog altered: %r" % (asked or None))

# and it must not be READING the prompt text to decide what to do. Checked
# against the handler's own body rather than the file, because the file also
# contains the comment explaining why it does not, which is what an earlier
# version of this check matched, and failed on.
import inspect                                                  # noqa: E402
body = inspect.getsource(tunnel._handler)
body = body.split('"""', 2)[-1] if body.count('"""') >= 2 else body
for word in ("password", "passcode", "push", "lower()"):
    check(word not in body.lower(),
          "the handler looks at the prompt text (%r in its body). The site's "
          "sign-in service writes that wording and changes it; a handler that "
          "recognized 'password' in a hidden prompt would one day hand the "
          "account password to the second step" % word)

# --- selftest and version -------------------------------------------------
check(entry.main(["x", "--selftest"]) == 0, "--selftest does not pass")
check(entry.main(["x", "--version"]) == 0, "--version does not work")

# --- one place says who made it -------------------------------------------
#
# ADCP and the suite it comes in are the Scripps Research Institute's; this is
# a GUI that drives them. The credit line is what a screenshot carries, so the
# distinction has to be in it.
for name, fn in (("the tool", tool_credit), ("the launcher", launcher_credit)):
    c = fn()
    for want in ("HOPE Lab", "ADCP", "Scripps Research Institute"):
        check(want in c, "%s's credit line is missing %r: %r" % (name, want, c))
    check("a project by HOPE Lab" not in c,
          "%s's credit line still reads as though ADCP were ours: %r"
          % (name, c))
check(tool_credit() == launcher_credit(),
      "the tool and the launcher disagree about who made them:\n  %s\n  %s"
      % (tool_credit(), launcher_credit()))
rendered = page.render("token")
# As a reader sees it. The header sets the sentence with a link on "ADCP" and
# bold on the owner, so the literal string is not in the source; the
# tag-stripped text is compared.
import re as _re3                                               # noqa: E402


def _flat(html):
    """The page as a reader sees it: tags gone, whitespace collapsed, and
    no space left in front of punctuation; a tag between a word and the
    comma after it would otherwise separate them."""
    t = " ".join(_re3.sub(r"<[^>]+>", " ", html).split())
    return _re3.sub(r"\s+([,.;:!?])", r"\1", t)


_text = _flat(rendered)
check(tool_credit() in _text, "the page does not carry the credit")
check("__CREDIT__" not in rendered,
      "the credit placeholder was left in the rendered page")

# The page itself has to name the engine and its owner where a screenshot of
# the header or footer will catch them, and say whose interface this is.
from adcp_dock import engine_credit                             # noqa: E402
for want in ("Scripps Research Institute", "Powered by HOPE Lab",
             "ccsb.scripps.edu"):
    check(want in rendered, "the page does not say %r anywhere" % want)
check(engine_credit() in _text or
      (tool_credit() in _text and "Not affiliated" in _text),
      "the page does not carry the engine credit and its disclaimer")
check("__ENGINE_CREDIT__" not in rendered,
      "the engine-credit placeholder was left in the rendered page")
# And citing it has to be reachable from the page, not only the README.
check("adcp-dock cite" in rendered,
      "the page never tells anyone how to cite this or what it ran")

# --- the credit line reads in the right order -----------------------------
#
# The engine and its owner come first; "Powered by HOPE Lab" ahead of them
# would read as though the lab were behind ADCP.
import adcp_launcher as _L                                      # noqa: E402
from adcp_dock import credit as _tool_credit                    # noqa: E402
_c = _L.credit()
check(_c.index("graphical interface to ADCP") < _c.index("Powered by HOPE Lab"),
      "the credit line does not lead with the interface-to-ADCP sentence")
check("Scripps Research Institute" in _c,
      "the credit line does not name Scripps")
# Nothing appended. Authors and PI belong in the page header and the
# launcher's own header, not at the end of the credit sentence.
check(_c.rstrip().endswith("Powered by HOPE Lab."),
      "something is appended after the credit line's sentence: %r" % _c)
for _who in ("Aadhil Haq", "Sandun Fernando"):
    check(_who not in _c,
          "%s is appended to the credit line instead of sitting in a "
          "header of its own" % _who)
# The people are still named where they belong.
from adcp_dock import page as _pg                               # noqa: E402
_r = _pg.render("tok")
for _who in ("Aadhil Haq", "Sandun Fernando"):
    check(_who in _r, "the page does not name %s anywhere" % _who)
check(_c == _tool_credit(),
      "the tool and the launcher print different credit lines:\n  %s\n  %s"
      % (_tool_credit(), _c))
check(_L.ENGINE_CREDIT in _c,
      "the launcher's credit line is not built from ENGINE_CREDIT, so the two "
      "can drift")
check("HOW_TO_CITE.txt" in _L.HOW_TO_CITE and "CITATION.bib" in _L.HOW_TO_CITE,
      "the launcher never says the citation files come with a download")
# A facility's acknowledgment is the cluster's to set, in activate.sh, which
# the command sources; the launcher's own table names no facility.
from adcp_launcher import clusters as _C                         # noqa: E402
check(not any(s.get("ack") for s in _C.SITES.values()),
      "the launcher's site table carries a facility's acknowledgment")
check(_C.sites()[-2:] == [_C.OTHER, _C.THIS] and not _C.SITES[_C.OTHER]["host"]
      and not _C.SITES[_C.OTHER]["runs"],
      "any other cluster and this computer do not come last, or the other "
      "cluster fills in an address or a runs folder: %r" % _C.sites())
_named = [n for n in _C.sites() if n not in (_C.OTHER, _C.THIS)]
check(_named == ["Grace (TAMU)", "ACES (TAMU)", "Garibaldi (Scripps)"]
      and all(_C.SITES[n]["host"] and _C.SITES[n]["install"] and _C.SITES[n]["runs"]
              for n in _named),
      "each named cluster brings its address, install and runs folder: %r" % _named)
_g = _C.defaults_for("Garibaldi (Scripps)", "ksojitra")
check(_g["host"] == "garibaldi.scripps.edu"
      and _g["install"] == "/gpfs/group/sanner/ksojitra/ADCP_docking_v2"
      and _g["runs"] == "/gpfs/group/sanner/ksojitra/adcp_runs_v2",
      "Garibaldi fills in the Sanner lab's address and folders: %r" % _g)
check(_C.server_env("Garibaldi (Scripps)", "ksojitra") == {
          "ADCP_RUNS": "/gpfs/group/sanner/ksojitra/adcp_runs_v2",
          "ADCP_SLURM_ACCOUNT": "sanner", "ADCP_SLURM_EXTRA": "--partition=highcpu"},
      "Garibaldi suggests the lab's account and partition to the server: %r"
      % _C.server_env("Garibaldi (Scripps)", "ksojitra"))
check(_C.server_env(_C.OTHER, "someone") == {} and _C.server_env("Grace (TAMU)", "") == {},
      "a cluster the launcher does not know, or one with no username yet, is told nothing")
_d = session.command("/x", defaults=_C.server_env("Garibaldi (Scripps)", "ksojitra"))
check('[ -n "${ADCP_SLURM_ACCOUNT:-}" ] || export ADCP_SLURM_ACCOUNT=sanner;' in _d
      and '[ -n "${ADCP_RUNS:-}" ] || export ADCP_RUNS=/gpfs/group/sanner/ksojitra/adcp_runs_v2;' in _d
      and "export ADCP_SLURM_EXTRA=--partition=highcpu;" in _d
      and _d.index(". ./activate.sh") < _d.index("ADCP_SLURM_ACCOUNT"),
      "the suggestions are set after activate.sh, and only where it set none: %s" % _d[:400])
check("'~" not in session.command("/x", defaults={"ADCP_RUNS": "~/adcp_runs"}),
      "a suggested runs folder under ~ is quoted so the tilde never expands")
check(_C.site_of("Some cluster an older version offered") == _C.OTHER
      and _C.site_of(_C.THIS) == _C.THIS,
      "a site saved by an older version is not read as the cluster")
_line = session.command("/x", env={"ADCP_ACKNOWLEDGEMENT": "Thanks, HPC."})
check("export ADCP_ACKNOWLEDGEMENT='Thanks, HPC.';" in _line,
      "the acknowledgment is not passed to the server: %s" % _line[:160])
check("ADCP_ACKNOWLEDGEMENT" not in session.command("/x"),
      "an empty acknowledgment is still exported")

# --- one name, read from each of three files ---------------------------------
#
# The launcher's header reads adcp_dock/brand.py; the dialogs and the page
# read APP; the launcher has its own APP because it installs on computers
# with no copy of adcp_dock. All three must agree.
from adcp_dock import brand as _brand, APP as _tool_app         # noqa: E402
check(_brand.NAME == _tool_app == _L.APP,
      "three names disagree: brand.NAME=%r adcp_dock.APP=%r launcher.APP=%r"
      % (_brand.NAME, _tool_app, _L.APP))
check("ADCP" in _brand.NAME and "GUI" in _brand.NAME,
      "the name no longer says this is a GUI to ADCP: %r" % _brand.NAME)
check("Scripps Research Institute" in _brand.TAGLINE
      and "HOPE Lab" in _brand.TAGLINE,
      "the tagline does not name Scripps and HOPE Lab: %r" % _brand.TAGLINE)
check(_brand.NAME in _r,
      "the page's heading is not the name everything else uses")


# --- the window says what is happening, not what it is running ------------
# The ssh command line, the host-key check and the server's stdout are the
# right information for one person, whoever is debugging a connection, and
# noise for everyone else. They go to the detail channel, which is off by
# default and is turned on automatically when something fails.
_app = open(os.path.join(ROOT, "adcp_launcher", "app.py")).read()
check("log=self.detail" in _app,
      "the sign-in machinery is still written into the visible log")
check("on_output=self.from_cluster" in _app,
      "the server's raw stdout still goes straight to the visible log")
check('self.q.put(("details", ""))' in _app,
      "a failed connection does not reveal the detail it needs reported")
check("Show technical details" in _app,
      "there is no way to see the technical detail at all")
check('self.say("Connecting to %s as %s..." % (host, user))' in _app,
      "the log does not open by saying it is connecting")
check("self.say(HOW_TO_CITE)" in _app,
      "the log never shows HOW_TO_CITE once connected")
# The sentence appears once in a window: credit() and engine_credit() are the
# same sentence, so a window that shows the label and also logs the sentence
# would show it twice.
check("self.say(ENGINE_CREDIT)" not in _app,
      "the credit sentence is logged as well as shown on the bottom label, "
      "so the launcher prints it twice")
check("text=credit()" in _app,
      "the bottom label no longer carries the credit line at all")

# The page footer likewise: the sentence appears once on the whole page.
_sentence = "A graphical interface to ADCP, a product of the Scripps"
_t2 = _flat(_pg.render("tok"))
check(_t2.count(_sentence) == 1,
      "the credit sentence appears %d times on the page; it belongs in the "
      "header and nowhere else" % _t2.count(_sentence))
_foot = _r[_r.index("<footer"):]
check("Not affiliated with or endorsed by Scripps" in _foot,
      "the footer dropped the Scripps disclaimer")
check(_sentence not in _foot,
      "the footer repeats the sentence the header already made")

# --- "This computer": adcp-dock local, found, started, read and stopped -----
from adcp_launcher import localrun, clusters as _cl               # noqa: E402
import tempfile as _tf2, stat as _st2                            # noqa: E402

check(_cl.sites()[-1] == _cl.THIS and _cl.SITES[_cl.THIS].get("local"),
      "This computer is not offered as a site, or not marked local")
W = lambda found: (lambda n: found.get(n))                       # noqa: E731
c = localrun.command("linux", which=W({"adcp-dock": "/x/adcp-dock"}))
check(c == ["/x/adcp-dock", "local", "--no-browser", "--port", "0"],
      "on Linux/macOS it does not run adcp-dock from PATH: %r" % c)
c = localrun.command("darwin", which=W({}), exists=lambda p: True)
check(c and c[0].endswith(".adcp-dock/venv/bin/adcp-dock"),
      "it does not fall back to the documented venv: %r" % c)
check(localrun.command("linux", which=W({}), exists=lambda p: False) is None,
      "with adcp-dock nowhere, it still returned a command")
c = localrun.command("windows", which=W({"wsl.exe": "C:/wsl.exe"}))
check(c and c[:4] == ["C:/wsl.exe", "-e", "bash", "-lc"]
      and "local --no-browser --port 0" in c[4] and "venv" in c[4],
      "on Windows it does not run through WSL with the venv fallback: %r" % c)
check(localrun.command("windows", which=W({})) is None,
      "without WSL, Windows still got a command")
for _os in ("windows", "darwin", "linux"):
    check("venv" in localrun.setup_text(_os) and "pip install" in
          localrun.setup_text(_os), "no setup steps for %s" % _os)

# The real adcp-dock, through a pipe. Python block-buffers stdout on a pipe,
# so a server that prints its address and then blocks must flush that line; a
# shell stand-in would not buffer and would hide the problem. This runs
# `serve`, which prints the address and blocks like `local` does, without
# `local`'s first-run installation of the ADCP 1.1 suite.
_d = _tf2.mkdtemp()
_seen = []
# Only where the whole tool is installed. The launcher build machines carry
# the launcher alone (no numpy, no RDKit), so the CLI cannot import there; it
# runs on the cluster and is skipped with a message elsewhere.
import importlib.util as _iu                                       # noqa: E402
_full = all(_iu.find_spec(m) for m in ("numpy", "rdkit"))
if not _full:
    print("  (the real-CLI pipe check needs numpy and RDKit; not on this "
          "machine, so it was skipped)")
try:
    if not _full:
        raise StopIteration
    _real = localrun.start(on_output=_seen.append, timeout=60,
                           cmd=[sys.executable, "-m", "adcp_dock.cli",
                                "serve", "--port", "0",
                                "--runs", os.path.join(_d, "runs")])
    check(_real.url.startswith("http://127.0.0.1:") and "?t=" in _real.url,
          "the real CLI's address came through the pipe malformed: %r"
          % _real.url)
    _real.stop()
except StopIteration:
    pass
except RuntimeError as _e:
    check(False, "the real CLI's address never came through a pipe: its "
                 "output is buffered, so the launcher would wait for ever: "
                 "%s" % str(_e)[:200])

# The address-reading, against a Python program that prints and then blocks
# without flushing by hand: the exact shape of the real server.
_fake = os.path.join(_d, "fake_local.py")
open(_fake, "w").write(
    "import sys, time\n"
    "sys.stdout.reconfigure(line_buffering=True)\n"
    "print('installing the ADCP 1.1 suite...')\n"
    "print('    http://127.0.0.1:43210/?t=AbC_d-9')\n"
    "time.sleep(60)\n")
_seen = []
_run = localrun.start(on_output=_seen.append, cmd=[sys.executable, _fake])
check(_run.url == "http://127.0.0.1:43210/?t=AbC_d-9",
      "the address was not read from what it printed: %r" % _run.url)
check("installing the ADCP 1.1 suite..." in _seen,
      "what it printed while starting was not passed on to the window")
_run.stop()
check(_run.proc.poll() is not None, "stop() left adcp-dock local running")

# A PYTHON stand-in, run by this interpreter: a #!/bin/sh script cannot be
# executed on Windows at all, by which this check first failed there.
_bad = os.path.join(_d, "dies.py")
open(_bad, "w").write("print('Missing: openmm')\nraise SystemExit(1)\n")
try:
    localrun.start(cmd=[sys.executable, _bad])
    check(False, "a start that died without an address was reported ready")
except RuntimeError as _e:
    check("Missing: openmm" in str(_e),
          "a failed start does not say what it printed: %s" % _e)

# A command that cannot be started at all (on Windows, wsl.exe with no Linux
# installed) is a setup message, not a traceback.
try:
    localrun.start(cmd=[os.path.join(_d, "no-such-program")])
    check(False, "a command that does not exist was reported ready")
except RuntimeError as _e:
    check("pip install" in str(_e),
          "an unstartable command does not say how to set up: %s" % _e)
except OSError as _e:
    check(False, "an unstartable command raised %s instead of explaining "
                 "what to install" % type(_e).__name__)
try:
    localrun.start(cmd=None) if localrun.command() is None else None
    check(localrun.command() is not None,
          "no adcp-dock here, and start() did not refuse")
except RuntimeError as _e:
    check("pip install" in str(_e), "the refusal does not say how to install")

# A child that prints nothing times out rather than holding the launcher
# indefinitely: the pipe is read on a separate thread, so the deadline is
# checked.
import time as _t3                                                # noqa: E402
_silent = os.path.join(_d, "silent.py")
open(_silent, "w").write("import time\ntime.sleep(60)\n")
_t0 = _t3.time()
try:
    localrun.start(cmd=[sys.executable, _silent], timeout=3)
    check(False, "a child that printed nothing was reported ready")
except RuntimeError as _e:
    check(_t3.time() - _t0 < 15,
          "a silent child held start() for %.0f s despite a 3 s timeout"
          % (_t3.time() - _t0))
    check("(nothing)" in str(_e), "the timeout does not say nothing was printed")

_appsrc = open(os.path.join(ROOT, "adcp_launcher", "app.py")).read()
check("clusters.THIS" in _appsrc and "def work_here" in _appsrc,
      "the window never takes the This computer path")

if fails:
    print("test_launcher FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the launcher's command, its URL pattern, its two-factor handling and the "
      "credit line all hold, %d checks" % ran)
