"""The loading screen, and the status card for ADapt's receptor.

By Aadhil Haq

While ADapt starts, the space its frame will fill shows a loading stage: a
quiet outline of ADapt's layout, and over it a card with a small molecule
drawing itself, the real steps ticking as they happen, the time so far, and a
line that the rest of the page can be used meanwhile. If the start fails the
card says why in plain words, with Try again. The same component, in a
compact form, stands in for every other wait that showed only a line of text.

What is checked here is what can be checked without a browser: the rules that
decide each step and each error, and the status card's state, run under Node
where Node can be found (ADCP_NODE, or node on the PATH); the server's word on
how far a start has got, and that a failed start can be tried again; and the
page's markup and style (role=status, reduced motion, the page's own colour
tokens, a phone). The screenshot script drives the page in a browser.

    python3 tests/test_waiting.py
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import adapt, page, server  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


html = page.render("TOKEN")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
css = html[html.index("<style>"):html.index("</style>")]
wcss = css[css.index("---- waiting"):css.index("---- ADapt's receptor, as a status card")]
cardcss = css[css.index("---- ADapt's receptor, as a status card"):
              css.index("---- column widths")]


def between(text, start, end):
    i = text.index(start)
    return text[i:text.index(end, i)]


wait_js = between(appjs, "/* ===================== waiting ====", "/* ===================== column widths")
adapt_wait_js = between(appjs, "/* ===================== waiting for ADapt", "/* The frame onto the session")

# --- one component, in two sizes ----------------------------------------------
check("function waitFor(o)" in appjs and "function waitIn(el, label, opt)" in appjs,
      "there is no shared loading component")
for part in ("ld-skel", "ld-card", "ld-mark", "ld-steps", "ld-time", "ld-slow", "ld-keep",
             "ld-errbox", "ld-retry"):
    check(part in wait_js, "the stage has no %s" % part)
check("'stage'" in wait_js and "'inline'" in wait_js,
      "the component does not come in a stage and an inline size")
# Drawn, not a picture: a small molecule in SVG with CSS animation.
check("<svg class=ld-mark" in wait_js and ".gif" not in wait_js
      and "@keyframes ld-draw" in wcss and "@keyframes ld-shimmer" in wcss,
      "the animated mark is not a drawing of its own")
# A screen reader hears each change of step through one polite line, and the
# clock is hidden from it.
check("role=status aria-live=polite" in wait_js,
      "the stage has no polite status line for a screen reader")
check("'aria-live', 'polite'" in wait_js and "setAttribute('role', 'status')" in wait_js,
      "the inline wait's place is not a polite status region")
check("class=ld-time aria-hidden=true" in wait_js
      and "<p class=ld-meta aria-hidden=true>" in wait_js,
      "the clock would be read out every second")
# Try again shows the keyboard's focus.
check(".ld-retry:focus-visible{outline" in wcss,
      "Try again does not show the keyboard's focus")
# Less motion asked for: the shimmer, the drawing and the spinners stop.
calm = re.search(r"@media \(prefers-reduced-motion:reduce\)\{(.*?)\n\}", wcss, re.S)
check(calm is not None, "the loading screen ignores prefers-reduced-motion")
if calm:
    for part in (".ld-sk-tree i", ".ld-mark .ld-bond", ".ld-mark .ld-spin",
                 ".ld-step.ld-now .ld-ic", ".ld-mini circle"):
        check(part in calm.group(1), "%s still moves with less motion asked for" % part)
    check("animation:none" in calm.group(1) and "stroke-dashoffset:0" in calm.group(1),
          "with less motion the mark is not drawn still and whole")
check("matchMedia('(prefers-reduced-motion: reduce)')" in wait_js,
      "the fade does not stop for prefers-reduced-motion")
# Colours from the page's tokens, so both themes follow.
for name, block in (("loading screen", wcss), ("status card", cardcss)):
    check(not re.search(r"#[0-9a-fA-F]{3,6}\b", block),
          "the %s has a colour of its own instead of a token: %r"
          % (name, re.findall(r"#[0-9a-fA-F]{3,6}\b", block)))
    for tok in set(re.findall(r"var\(--([\w-]+)\)", block)):
        check("--%s:" % tok in page.DARK_VARS and "--%s:" % tok in css,
              "the %s uses --%s, which a theme does not define" % (name, tok))
# A phone: the card fits, the outline of ADapt's tree goes.
check(re.search(r"@media \(max-width:600px\)\{\s*\.ld-skel\{grid-template-columns:"
                r"minmax\(0,1fr\)\}", wcss) is not None
      and ".ld-card{position:relative;width:min(440px,100%)" in wcss,
      "the loading card does not fit a phone")

# --- ADapt's start, as it is shown -------------------------------------------
check("adaptWaitBegin(src, session, gen);" in appjs,
      "opening ADapt does not show the loading stage")
check("w.show(wrap, 'stage');" in adapt_wait_js,
      "the stage is not in the space the frame will fill")
check("'/api/adapt/status'" in adapt_wait_js,
      "the steps are not driven by the server's word on the start")
check("adaptViewDrawn(doc);" in appjs and "addEventListener('load'" in appjs,
      "the last step is not the frame's own load")
check("You can keep going below while it starts." in adapt_wait_js,
      "the card does not say the page can be used meanwhile")
check("A first start on a busy login node takes 1 to 2 minutes" in adapt_wait_js
      and "slow: [20," in adapt_wait_js,
      "there is no reassurance after twenty seconds")
check("GUIDE_BUSY.receptor" in adapt_wait_js and "li.classList.toggle('busy'" in appjs,
      "the step guide does not show ADapt starting")
check("id=recbusy" in adapt_wait_js and "#step-receptor > h2" in adapt_wait_js,
      "the Receptor card's heading has no badge while ADapt starts")
check(".guide li.busy .gnum::after" in css,
      "the strip has no busy state for its first step")
# Nothing is blocked while it starts.
for bad in ("inert", "pointerEvents", ".disabled = true", "showModal"):
    check(bad not in adapt_wait_js, "starting ADapt blocks the page (%s)" % bad)
# The full screen tab waits on the same card.
check("adaptTabOpen()" in appjs and "class=ld-page" in appjs
      and "A.tabs.forEach(t => { try { t.location.href = f.url; }" in appjs,
      "the full screen tab does not wait on the same card")
# The error and retry path: Try again opens the structure again, and a
# newer wait in the same place replaces the one before.
check("e.retry ? () => adaptOpenFor(A.src, A.session) : null" in adapt_wait_js,
      "Try again does not start ADapt again")
check("b.onclick = () => { const again = e.retry; w.drop(); if(again) again(); };"
      in wait_js, "the Try again button does not call the retry")
check("if(was) was.w.drop();" in adapt_wait_js,
      "a retry leaves the failed card behind")
check("ADAPT_LIMIT_MS" in appjs and "{late: true}" in appjs,
      "a start that never answers is waited on for ever")
# The old line of grey text is gone.
check("A first start can take a " not in appjs and "function adaptWait(" not in appjs,
      "the line of text the loading screen replaces is still there")

# --- the other waits use the same component --------------------------------------
for host, words in (("recinfo", "Fetching ' + c + ' from the RCSB"),
                    ("sitenote", "Checking the site"),
                    ("surfnote", "Computing the surface with MSMS"),
                    ("targetnote", "Reading the target file"),
                    ("estnote", "Working out how long this will take"),
                    ("submsg", "Submitting to the queue"),
                    ("viewmsg", "Loading pose ' + rank")):
    check(words in appjs, "the wait for %s does not use the component" % host)
check("busy.join(' and the ')" in appjs and "waitKeep(el," in appjs,
      "the pose view's surfaces do not use the component")
for old in ("'fetching ' + c + ' from the RCSB...'", "'computing the surface with MSMS...'",
            "'reading the target file…'", "'submitting…'", "'reading…'"):
    check(old not in appjs, "a wait still shows only text: %s" % old)

# --- no dash punctuation in what was added -----------------------------------------
added = wait_js + adapt_wait_js + between(appjs, "function adaptCardState", "/* The pH belongs")
strings = re.findall(r"'([^'\n]*)'", added)
for t in strings:
    check("\u2014" not in t and "\u2013" not in t and " - " not in t,
          "dash punctuation in %r" % t)

# --- the rules, under Node -------------------------------------------------------------
NODE = os.environ.get("ADCP_NODE") or shutil.which("node")


def run_js(fns, calls):
    """Run named functions from app.js on each case; their answers, or None."""
    src = []
    for f in fns:
        if f.startswith("const "):
            src.append(re.search(r"%s = [^;]*;" % re.escape(f), appjs).group(0))
            continue
        i = appjs.index("function %s(" % f)
        src.append(appjs[i:appjs.index("\n}\n", i) + 3])
    prog = ("\n".join(src) + "\nconst cases = JSON.parse(process.argv[1]);"
            "\nprocess.stdout.write(JSON.stringify(cases.map(c => " + calls + ")));\n")
    out = subprocess.run([NODE, "-e", prog, json.dumps(CASES_NOW)],
                         capture_output=True, text=True, timeout=60)
    if out.returncode:
        check(False, "app.js could not be run under Node: %s" % out.stderr[-400:])
        return None
    return json.loads(out.stdout or "null")


if NODE and os.path.isfile(NODE if os.sep in NODE else shutil.which(NODE) or ""):
    # adaptStartSteps(facts): each step from its own event, a later event
    # implying the earlier ones, and the first step not done the one under way.
    STEPCASES = [
        ({}, {"start": "now", "answer": "wait", "hand": "wait", "view": "wait"}, False),
        ({"process": "starting"},
         {"start": "ok", "answer": "now", "hand": "wait", "view": "wait"}, False),
        ({"process": "running"},
         {"start": "ok", "answer": "now", "hand": "wait", "view": "wait"}, False),
        ({"process": "running", "answering": True},
         {"start": "ok", "answer": "ok", "hand": "now", "view": "wait"}, False),
        ({"session": True}, {"start": "ok", "answer": "ok", "hand": "ok", "view": "now"}, False),
        ({"session": True, "view": True},
         {"start": "ok", "answer": "ok", "hand": "ok", "view": "ok"}, True),
        # A phone opens ADapt in a tab: three steps, done once handed over.
        ({"session": True, "inline": False}, {"start": "ok", "answer": "ok", "hand": "ok"}, True),
        # Failed: the step under way is the one marked, nothing after it begun.
        ({"error": True}, {"start": "no", "answer": "wait", "hand": "wait", "view": "wait"}, False),
        ({"process": "running", "answering": True, "error": True},
         {"start": "ok", "answer": "ok", "hand": "no", "view": "wait"}, False),
        # A process that died is not a step done.
        ({"process": "failed"}, {"start": "now", "answer": "wait", "hand": "wait", "view": "wait"},
         False),
    ]
    CASES_NOW = [c[0] for c in STEPCASES]
    got = run_js(["const ADAPT_STEPS", "adaptStartSteps"], "adaptStartSteps(c)")
    if got is not None:
        for (facts, steps, ready), g in zip(STEPCASES, got):
            check(g["steps"] == steps and g["ready"] == ready,
                  "adaptStartSteps(%r) gave %r, ready %r; wanted %r, ready %r"
                  % (facts, g["steps"], g["ready"], steps, ready))
        check(got[7]["failed"] == "start" and got[8]["failed"] == "hand",
              "a failed start does not say which step failed: %r" % [got[7], got[8]])

    # adaptStartError(answer, what): plain words, and whether to try again.
    ERRCASES = [
        ({"error": "ADapt did not start. Its last words: no ADCP environment found"},
         "ADapt did not start", True, "no ADCP environment found"),
        ({"error": "ADapt did not answer (Connection refused)"}, "ADapt stopped answering",
         True, "Connection refused"),
        ({"note": "ADapt reads PDB files only, and this receptor is x.cif; load a PDB"},
         "ADapt cannot open this file", False, ""),
        ({"note": "ADapt could not fetch 9ZZZ here (no PDB format file)"},
         "ADapt could not fetch 7A3R", True, "9ZZZ"),
        ({"error": "ADapt could not open x.pdb: bad file"}, "ADapt could not open 7A3R", True,
         "bad file"),
        ({"lost": "Failed to fetch"}, "The page lost touch with its server", True,
         "Failed to fetch"),
        ({"late": True}, "ADapt is taking far too long", True, ""),
        ({"error": "ADapt runs on ADCP 1.1, which is not installed here."},
         "ADapt cannot run here", False, ""),
        ({}, "ADapt did not open 7A3R", True, ""),
    ]
    CASES_NOW = [c[0] for c in ERRCASES]
    got = run_js(["adaptStartError"], "adaptStartError(c, '7A3R')")
    if got is not None:
        for (answer, what, retry, detail), g in zip(ERRCASES, got):
            check(g["what"] == what and g["retry"] is retry and detail in g["detail"],
                  "adaptStartError(%r) gave %r" % (answer, g))
            check(g["why"] and not re.search("[\u2013\u2014]| - ", g["why"] + g["what"]),
                  "the error for %r is not plain words: %r" % (answer, g))

    # adaptCardState(facts): the status card's state, and so its one action.
    CARDCASES = [
        ({"rec": None}, {"state": "none"}),
        ({"rec": {"save": 1}}, {"state": "current"}),
        ({"rec": {"save": 1}, "dirty": True, "edits": 2}, {"state": "pending", "edits": 2}),
        ({"rec": {"save": 1}, "dirty": True}, {"state": "pending", "edits": 0}),
        ({"rec": {"save": 1}, "unread": {"n": 2, "why": "its PDBQT is missing"}},
         {"state": "newer", "n": 2, "why": "its PDBQT is missing"}),
        # a save already in use is not newer than itself
        ({"rec": {"save": 2}, "unread": {"n": 2, "why": "x"}}, {"state": "current"}),
        # a newer save to read comes before changes made since
        ({"rec": {"save": 1}, "dirty": True, "unread": {"n": 3, "why": ""}},
         {"state": "newer", "n": 3, "why": ""}),
        ({"rec": {"save": 1}, "gone": True}, {"state": "earlier"}),
        ({"rec": {"save": 1}, "gone": True, "dirty": True, "edits": 1},
         {"state": "pending", "edits": 1}),
    ]
    CASES_NOW = [c[0] for c in CARDCASES]
    got = run_js(["adaptCardState"], "adaptCardState(c)")
    if got is not None:
        for (facts, want), g in zip(CARDCASES, got):
            check(g == want, "adaptCardState(%r) gave %r, not %r" % (facts, g, want))

    # ldTime: the time so far as a person reads it.
    CASES_NOW = [0, 7.9, 59, 60, 65, 3725]
    got = run_js(["ldTime"], "ldTime(c)")
    if got is not None:
        check(got == ["0 s", "7 s", "59 s", "1 min 00 s", "1 min 05 s", "62 min 05 s"],
              "the clock reads %r" % got)
else:
    print("note  no Node here (set ADCP_NODE), so the rules were checked by reading "
          "them only; the screenshot script drives them in a browser")

# --- the server's word on how far a start has got ------------------------------------
tmp = tempfile.mkdtemp(prefix="waiting_test_")
os.environ["HOME"] = tmp


def fake(folder, script):
    os.makedirs(os.path.join(folder, "bin"), exist_ok=True)
    os.makedirs(os.path.join(folder, "adapt"), exist_ok=True)
    open(os.path.join(folder, "adapt", "__main__.py"), "w").write("s.add_argument('--base-path')\n")
    path = os.path.join(folder, "bin", "adapt")
    open(path, "w").write("#!/bin/sh\n" + script + "\n")
    os.chmod(path, 0o755)
    return folder


r = adapt.Runner()
check(r.phase() == {"process": "none"}, "a runner never started says none: %r" % r.phase())
slow = fake(os.path.join(tmp, "slow"), 'mkdir -p "$ADAPT_OUTDIR"\nsleep 2\n'
            'echo "    http://127.0.0.1:45679/adapt/?t=abc"\nexec sleep 60')
import threading  # noqa: E402
th = threading.Thread(target=lambda: r.start(slow, "/envs/x", os.path.join(tmp, "slow_out"),
                                             timeout=30))
th.start()
seen = set()
for _ in range(60):
    seen.add(r.phase()["process"])
    if r.running():
        break
    time.sleep(0.1)
th.join(30)
check("starting" in seen, "a start under way is not reported as starting: %r" % seen)
check(r.phase() == {"process": "running"}, "a started runner is not running: %r" % r.phase())
r.stop()

bad_flag = os.path.join(tmp, "may_start")
flaky = fake(os.path.join(tmp, "flaky"),
             'if [ ! -e "%s" ]; then echo "ADapt: no ADCP environment found" >&2; exit 2; fi\n'
             'mkdir -p "$ADAPT_OUTDIR"\necho "    http://127.0.0.1:45680/adapt/?t=abc"\n'
             'exec sleep 60' % bad_flag)
r2 = adapt.Runner()
try:
    r2.start(flaky, "/envs/x", os.path.join(tmp, "flaky_out"), timeout=20)
    check(False, "a start that dies is not reported")
except RuntimeError as exc:
    check(str(exc).startswith("ADapt did not start") and "no ADCP environment" in str(exc),
          "the failure does not carry ADapt's words: %s" % exc)
ph = r2.phase()
check(ph.get("process") == "failed" and "no ADCP environment" in ph.get("why", ""),
      "a failed start is not reported as failed, with why: %r" % ph)
# Try again really starts it again, and a start that works clears the failure.
open(bad_flag, "w").write("1")
url = r2.start(flaky, "/envs/x", os.path.join(tmp, "flaky_out"), timeout=20)
check(url.endswith("?t=abc") and r2.phase() == {"process": "running"},
      "trying again after a failed start does not start ADapt: %r" % r2.phase())
r2.stop()

# Through the interface: /api/adapt/status reads the runner and asks a
# running ADapt to answer, and never starts anything itself.
app = server.App(os.path.join(tmp, "runs"))
check(app.adapt_status() == {"process": "none", "answering": False},
      "with nothing started the status is not none: %r" % app.adapt_status())
check(getattr(app, "_adapt_runner", None) is None, "asking the status started ADapt")
import http.server  # noqa: E402


class Answers(http.server.BaseHTTPRequestHandler):
    def log_message(self, *a):
        pass

    def do_GET(self):
        body = b'{"sessions": []}'
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Answers)
threading.Thread(target=srv.serve_forever, daemon=True).start()
os.environ["http_proxy"] = os.environ["HTTP_PROXY"] = "http://127.0.0.1:9"


class Up:
    url = "http://127.0.0.1:%d/adapt/?t=TOK" % srv.server_address[1]

    def phase(self):
        return {"process": "running"}

    def stop(self):
        pass


app._adapt_runner = Up()
check(app.adapt_status() == {"process": "running", "answering": True},
      "a running ADapt that answers is not reported as answering: %r" % app.adapt_status())
check(adapt.answering(Up.url), "ADapt's answer is not counted")
srv.shutdown()
srv.server_close()
check(not adapt.answering(Up.url, timeout=2), "no answer at all is counted as one")
check("/api/adapt/status" in open(os.path.join(ROOT, "adcp_dock", "server.py")).read(),
      "the page cannot ask how far the start has got")

# --- the status card --------------------------------------------------------------------
card = between(html, "<div id=adaptuse", "<!-- The frame:")
check("Docking will use ADapt&rsquo;s receptor" in card,
      "the card's headline is not what docking will use")
check('class="primary hide" id=adaptupdate' in card
      and "Use ADapt&rsquo;s current receptor" in card,
      "edits since the save have no primary button to use them")
check('class="primary hide" id=adaptread' in card,
      "a newer save has no button to read it")
check("class=secondary id=adaptback" in card
      and "Use this page&rsquo;s preparation instead" in card,
      "the way back is not a clear secondary button")
check("<details class=aumore" in card and "<summary>Details</summary>" in card
      and "id=adaptwhat" in card,
      "the long details are not folded under Details")
check("id=auchips" in card and "Chains" in appjs and "Hetero groups kept" in appjs
      and "Changes in this session" in appjs, "the key facts are not chips")
check("Up to date with ADapt" in appjs, "an up to date card does not say so")
check("Read the newest save again" not in html + appjs,
      "the button nobody saw is still there")
check(html.index("id=adaptuse") < html.index("id=adaptpanel"),
      "the status card is below the frame, out of sight")
check("f.unread = {n: newest.n" in appjs,
      "a newer save that could not be read is not remembered")

shutil.rmtree(tmp, ignore_errors=True)
if fails:
    print("test_waiting FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the loading screen, its steps and errors, and ADapt's status card, %d checks" % ran)
