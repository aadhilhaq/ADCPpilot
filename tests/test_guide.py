"""The step guide, and the small things fixed beside it.

By Aadhil Haq

The guide is five numbers across the top of the page and the five cards they
stand for. What is checked here is what can be checked without a browser:
that the strip and the cards agree on the steps and their order, that a done
step is marked by more than its color, that gray cards are only gray (never
disabled), that motion stops for somebody who asked for less of it, and the
rule that decides each step. That rule is a pure function in app.js; it is
run under Node where Node can be found (ADCP_NODE, or node on the PATH), and
the page itself is driven in a browser by the screenshot script.

Then the smaller fixes: the site line in angstroms without repeating the
box, ranges written "23 to 29", the estimate stating what the jobs actually
book, the target file's hint, and the phone layout.

    python3 tests/test_guide.py
"""
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_pilot import page, pipeline, server, sizing, slurm, structure  # noqa: E402
from adcp_pilot import job as jobmod                               # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


html = page.render("TOKEN")
appjs = open(os.path.join(ROOT, "adcp_pilot", "static", "app.js")).read()
css = html[html.index("<style>"):html.index("</style>")]
STEPS = ["receptor", "site", "peptides", "docking", "submit"]
LABELS = ["Receptor", "Site", "Peptides", "Docking", "Submit"]

# --- the strip and the cards agree --------------------------------------------
nav = re.search(r"<nav class=guide id=guide[^>]*>(.*?)</nav>", html, re.S)
check(nav is not None, "the page has no step guide")
nav = nav.group(1) if nav else ""
check(re.findall(r"<li data-step=(\w+)>", nav) == STEPS,
      "the strip's steps are %r" % re.findall(r"<li data-step=(\w+)>", nav))
check(re.findall(r"class=glab>([^<]+)<", nav) == LABELS,
      "the strip reads %r" % re.findall(r"class=glab>([^<]+)<", nav))
check(re.findall(r"class=gn>(\d)<", nav) == ["1", "2", "3", "4", "5"],
      "the strip is not numbered 1 to 5")
check(re.findall(r"data-goto=(\w+)", nav) == STEPS,
      "a number does not lead to its own card")
# Done is shown by a tick drawn in the number, and said in words for a
# screen reader, not by the green alone.
check(nav.count("class=gtick") == 5,
      "not every number has a tick to show when it is done")
check(nav.count('class="vh gstate"') == 5,
      "not every number says its state in words")
cards = re.findall(r'<div class=card id=step-(\w+) data-step=(\w+)>', html)
check([c[0] for c in cards] == STEPS and [c[1] for c in cards] == STEPS,
      "the cards that are steps are %r, not the strip's steps in its order"
      % cards)
check(html.index("id=guide") < html.index("id=step-receptor"),
      "the strip is not above the cards")
for k, label in zip(STEPS, ["Receptor", "Site", "Peptides",
                            "Docking and rescoring", "Submit"]):
    i = html.find("id=step-%s " % k)
    check(i >= 0 and html.find("<h2>%s</h2>" % label, i) == html.find("<h2>", i),
          "the %s card's heading is not %r" % (k, label))
# The submit button is in the Submit card, the search in the Docking card.
_sub = html[html.index("id=step-submit"):html.index("<h2>Runs</h2>")]
for want in ("id=submit", "id=label", "id=outdir", "id=account"):
    check(want in _sub, "%s is not in the Submit card" % want)
_dock = html[html.index("id=step-docking"):html.index("id=step-submit")]
for want in ("id=poses", "id=replicas", "id=mmgbsa", "id=walltime",
             "id=autofill", "id=estnote"):
    check(want in _dock, "%s is not in the Docking card" % want)

# --- gray is a cue, never a lock ------------------------------------------------
for rule in re.findall(r"\.card\.step-later[^{]*\{[^}]*\}", css):
    check("pointer-events" not in rule and "display:none" not in rule
          and "visibility" not in rule,
          "a gray card is made unusable: %s" % rule)
check(".card.step-later:focus-within" in css and ".card.step-later:hover" in css,
      "a gray card does not come back to full color when it is used")
_guide_js = appjs[appjs.index("function guide(){"):
                  appjs.index("function guideSoon(")]
for bad in (".disabled", "inert", "pointerEvents", "pointer-events",
            "readOnly"):
    check(bad not in _guide_js, "the guide sets %s on something" % bad)
# Motion stops for somebody who asked for less of it, in the strip, the
# cards and the scrolling.
_calm = re.search(r"@media \(prefers-reduced-motion:reduce\)\{([^@]*?)\}\s*\}",
                  css[css.index("---- the cards follow the guide"):])
check(_calm and ".gnum" in _calm.group(1) and ".card[data-step]" in _calm.group(1),
      "the guide's transitions do not stop for prefers-reduced-motion")
check("prefers-reduced-motion: reduce" in appjs[appjs.index("data-goto"):
                                                appjs.index("let ADCP11")],
      "a step number scrolls smoothly even when less motion was asked for")
# Colors from the page's own tokens, so both themes follow.
_gcss = css[css.index("---- the step guide"):css.index("---- the runs: two views")]
check(not re.search(r"#[0-9a-fA-F]{3,6}\b", _gcss),
      "the guide's style has a color of its own instead of a token: %r"
      % re.findall(r"#[0-9a-fA-F]{3,6}\b", _gcss))
for tok in set(re.findall(r"var\(--([\w-]+)\)", _gcss)):
    check("--%s:" % tok in page.DARK_VARS and "--%s:" % tok in css,
          "the guide uses --%s, which a theme does not define" % tok)
# A phone: the strip compacts and no longer follows the page down.
check(re.search(r"@media \(max-width:600px\)[^{]*\{\.guide\{position:static\}",
                css) is not None,
      "the strip still follows the page down on a phone")
check(".guide li{flex:1 1 0;min-width:0" in css,
      "the strip's steps can grow wider than a phone")

# The strip leaves with its cards: it is sticky inside a box holding only the
# strip and the five cards, so it never covers the runs or their results.
_flow = html[html.index("<div class=stepsflow>"):html.index("</div><!-- .stepsflow -->")]
check("id=guide" in _flow and all("id=step-%s " % k in _flow for k in STEPS)
      and "<h2>Runs</h2>" not in _flow,
      "the strip and its five cards are not in one box of their own")

# --- the rule that decides each step --------------------------------------------
check("function guideSteps(f)" in appjs, "app.js has no guideSteps")
for fact in ("TARGET.path === trg", "RECEPTOR.path", "SITE_OK === siteMode()",
             "PEPS_READY", "!dockingProblem()", "SUBMITTED"):
    check(fact in appjs[appjs.index("function guideFacts"):
                        appjs.index("function guide()")],
          "a step no longer reads %s" % fact)
check("PEPS_READY = !!d.ready && n > 0;" in appjs,
      "the peptides step is not the server's own verdict on the list")
check("SITE_OK = d.error ? null : mode;" in appjs,
      "the site step is not the server's check of the site")
check("SUBMITTED = true;" in appjs[appjs.index("$('submit').onclick"):],
      "a submission does not mark the last step")

NODE = os.environ.get("ADCP_NODE") or shutil.which("node")
CASES = [
    # (facts, the steps done, the step to do now)
    ({}, [], "receptor"),
    # Docking's defaults are valid, but it is not done on an empty page.
    ({"docking": True}, [], "receptor"),
    ({"receptor": True, "docking": True}, ["receptor"], "site"),
    # A target file: receptor and site at once.
    ({"receptor": True, "site": True, "docking": True},
     ["receptor", "site"], "peptides"),
    # Peptides before a receptor: done, but the receptor is still next.
    ({"peptides": True, "docking": True}, ["peptides"], "receptor"),
    ({"receptor": True, "site": True, "peptides": True, "docking": True},
     ["receptor", "site", "peptides", "docking"], "submit"),
    ({"receptor": True, "site": True, "peptides": True, "docking": False},
     ["receptor", "site", "peptides"], "docking"),
    ({"receptor": True, "site": True, "peptides": True, "docking": True,
      "submitted": True}, STEPS, ""),
]
if NODE and os.path.isfile(NODE if os.sep in NODE else shutil.which(NODE) or ""):
    steps_src = re.search(r"const STEPS = \[[^\]]*\];", appjs).group(0)
    fn = appjs[appjs.index("function guideSteps(f){"):]
    fn = fn[:fn.index("\n}\n") + 3]
    prog = (steps_src + "\n" + fn + "\nconst cases = JSON.parse(process.argv[1]);"
            "\nprocess.stdout.write(JSON.stringify(cases.map(guideSteps)));\n")
    try:
        out = subprocess.run([NODE, "-e", prog,
                              json.dumps([c[0] for c in CASES])],
                             capture_output=True, text=True, timeout=60)
        got = json.loads(out.stdout or "null")
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        got = None
        check(False, "guideSteps could not be run under Node: %s" % exc)
    if got is not None:
        for (facts, done, now), g in zip(CASES, got):
            check(sorted(k for k, v in g["done"].items() if v) == sorted(done)
                  and g["now"] == now,
                  "guideSteps(%r) gave done %r and now %r, not %r and %r"
                  % (facts, sorted(k for k, v in g["done"].items() if v),
                     g["now"], sorted(done), now))
else:
    print("note  no Node here (set ADCP_NODE), so guideSteps was checked by "
          "reading it only; the screenshot script drives it in a browser")

# --- the site line: angstroms, and the box said once --------------------------------
REC = os.environ.get("ADCP_TEST_RECEPTOR") or os.path.join(ROOT, "examples", "receptor.pdb")
if not os.path.isfile(REC):
    # examples/receptor.pdb is the user's own and ignored by git, so a fresh clone has none.
    # The package's own example ships MDM2 (chain A) with p53's peptide (chain B), a site in A
    # and a ligand in B, which is what these checks ask of a receptor.
    REC = os.path.join(ROOT, "adcp_pilot", "example_run", "receptor", "1YCR.pdb")
atoms = structure.read(REC)
_c, size, note = structure.box_for(
    [a for a in atoms if a.chain == "A"],
    structure.parse_hotspots("A:54, A:58, A:61, A:62", "A"))
check("Å" in note and not re.search(r"\d A\b", note),
      "the site line writes A for angstrom: %r" % note)
check("box" not in note and "%.1f x" % size[0] not in note,
      "the site line repeats the box, which the next line states: %r" % note)
lig = structure.ligand_atoms(atoms, "B")
_c, size, note = structure.ligand_box(lig)
check(note.startswith("%d ligand atoms" % len(lig)) and "Å" in note
      and "box" not in note and not re.search(r"\d A\b", note),
      "the ligand site line is %r" % note)
# The page states the box on the line after it, from the same numbers.
check("'; box '" in appjs and "d.size.map(x => x.toFixed(1)).join(' × ')"
      in appjs, "the page no longer states the box after the site line")
# The prepare log still states the ligand box, beside the line.
_pl = open(os.path.join(ROOT, "adcp_pilot", "pipeline.py")).read()
check('site["ligand_note"]) + tuple(size)' in _pl,
      "the prepare log lost the ligand box's size")

# --- ranges read "23 to 29" ----------------------------------------------------------
w = sizing.box_warning(13, [26, 26, 26])
check(w and " to " in w and not re.search(r"\d\s*[-–—]\s*\d", w),
      "the box warning writes a range with a dash: %r" % w)
check("1-3 million" not in appjs and "1 to 3 million" in appjs,
      "the steps guideline writes its range with a hyphen")

# --- the estimate states what the jobs book -------------------------------------------
def mins(t):
    h, m, s = (int(x) for x in t.split(":"))
    return h * 60 + m + (1 if s else 0)


for throttle in (0, 1):
    e = sizing.suggest([13, 11, 9], 50, 2500000, 48, mmgbsa=True, mm_poses=5,
                       throttle=throttle)
    b = e.get("booked") or {}
    check(set(b) >= {"prepare", "dock", "waves", "report", "total"},
          "the estimate does not say what each job books: %r" % b)
    if b:
        check(b["waves"] == (3 if throttle == 1 else 1),
              "three peptides one at a time are %r waves" % b["waves"])
        check(mins(b["total"]) == mins(b["prepare"]) + mins(b["dock"])
              * b["waves"] + mins(b["report"]),
              "the total booked is not the sum of what the jobs book: %r" % b)
        check(e["end_to_end"] == b["total"],
              "end_to_end is %r and the jobs book %r" % (e["end_to_end"],
                                                         b["total"]))
        check(b["dock"] == e["walltime"] and b["prepare"] == e["prep_walltime"],
              "the parts are not the walltimes the form and the prepare job "
              "are given: %r" % b)
        check(b["report"] == sizing.REPORT_WALLTIME,
              "the report is counted at %r and booked at %r"
              % (b["report"], sizing.REPORT_WALLTIME))
check('time=sizing.REPORT_WALLTIME' in _pl and 'time="01:00:00"' not in _pl,
      "the report job is not booked from sizing.REPORT_WALLTIME")
# Through the server, the prepare job is counted as the submission books it.
app = server.App(os.path.join(HERE, "_work", "guide_runs"))
body = {"peptides": [{"sequence": "ETFSDLWKLLPEN"}], "replicas": 50,
        "steps": 2500000, "cores": 48, "receptor": REC, "chains": ["A"],
        "box_max": 30, "mmgbsa": False}
e = app.estimate_run(dict(body))
j = jobmod.default(**{k: v for k, v in body.items() if k in jobmod.default()})
j["tools_roots"] = app.tools_roots
# What the submission books is prep_time capped at the partition's own limit (slurm.within_limit),
# so on a cluster whose default partition allows less than the four hours, both say less.
booked = slurm.within_limit(pipeline.prep_time(j), j.get("partition"), j.get("sbatch_extra"))
check(e.get("prep_walltime") == booked,
      "the estimate counts the prepare job at %r and the submission books %r"
      % (e.get("prep_walltime"), booked))
e = app.estimate_run(dict(body, target_file="/nowhere/target.trg"))
check(e.get("prep_walltime") == slurm.within_limit("04:00:00"),
      "a target file's prepare job is counted at %r" % e.get("prep_walltime"))
check("of queue time booked in total" not in appjs and "b.total" in appjs,
      "the page still states the chain's estimate as if it were one booking")

# --- a target file's site view does not ask for a residue -------------------------------
_ph = appjs[appjs.index("function pickHints()"):
            appjs.index("function siteModeChanged()")]
_tbranch = _ph[_ph.index("TARGET"):_ph.index(": (PICKHINT[m]")]
check("Click" not in _tbranch and "no residues are picked" in _tbranch,
      "with a target file the hint still invites a click: %r" % _tbranch)
_tc = appjs[appjs.index("function targetChanged()"):
            appjs.index("function targetSite()")]
check("pickHints()" in _tc, "loading a target file does not change the hint")

# --- a phone ------------------------------------------------------------------------------
check("#runs{overflow-x:auto}" in css,
      "the table of runs pushes the page sideways on a phone")
check(re.search(r"@media \(max-width:600px\)\{\s*\.seq2 \.ln>span:last-child"
                r"\{white-space:normal", css) is not None,
      "the residue strip runs out of its box on a phone")
check(".ftable{overflow:auto" in css,
      "the Folders table pushes the page sideways on a phone")

if fails:
    print("test_guide FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the step guide follows the page's checks, and the small fixes "
      "beside it hold, %d checks" % ran)
