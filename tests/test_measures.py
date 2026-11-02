"""Add measures: what finished runs lack, read for every run, and added to the
runs chosen in one press.

By Aadhil Haq

GET /api/measures lists every run the page lists, with what each has of the
peptide strain, the interface measures and the MM-GBSA engines installed here,
what it lacks and what it cannot have, and why no measure can be added to it at
all; read from job.json and the result files only. POST /api/measures adds
the chosen measures to the chosen runs: refused, with the reason, while a tidy
moves a run's files or another of its jobs is queued or running; otherwise one
short job for the light measures a run lacks, one job for each engine, the
report rebuilt after all of them, and the run's own tidy queued again behind.

What is checked, with the scheduler stood in for: the states and the reasons;
the jobs each run gets, their cores, walltime and stage, and what they wait
for; that nothing a run already has is queued again; the local executor; the
routes; the run's own page carrying what it lacks; and the page itself: the
dialog's markup, its words and its sizes, the controls it replaces gone, and
its rules (what a run lacks, what Add would queue, what the lines at the top
and over a run's table say) run under Node where Node is found (ADCP_NODE, or
node on the PATH). tests/test_measures_browser.py drives it in Chrome.

    python3 tests/test_measures.py
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import urllib.error
import urllib.parse
import urllib.request
from http.server import ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

try:
    from adcp_dock import job as jobmod, page, pipeline, server, sizing, slurm, tidy  # noqa: E402
    from adcp_dock._vendor import gbengines, strain as strainmod  # noqa: E402
except ImportError as exc:
    print("skip: the package cannot be imported here (%s)" % exc)
    sys.exit(77)

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


WORK = tempfile.mkdtemp(prefix="adcp_measures_")
RUNS = os.path.join(WORK, "runs")
os.makedirs(RUNS)
FILES = ("_com.parm7", "_rec.parm7", "_lig.parm7", "_min.rst7")


def make(name, peptides, engine="amber", ran=None, cpus=24, files=True,
         results=True):
    """A finished run: [(sequence, ring, scored poses, has strain, has panel)]."""
    rd = os.path.join(RUNS, name)
    j = jobmod.default(mmgbsa=True, mmgbsa_engine=engine, label=name, cpus=cpus,
                       peptides=[{"name": s, "sequence": s, "cyclic": ring}
                                 for s, ring, *_ in peptides])
    j["created"] = "2026-09-22 13:02:24"
    if ran:
        j["mmgbsa_engines_run"] = ran
    j["state"] = "done"
    jobmod.save(rd, j)
    if not results:
        return rd
    os.makedirs(os.path.join(rd, "results"))
    for seq, ring, n, has_strain, has_panel in peptides:
        r = {"peptide": seq, "cyclic": ring, "modes": [{"rank": 1, "score": -10.0}],
             "mmgbsa": [-30.0 - k for k in range(n)],
             "mmgbsa_ranks": list(range(1, n + 1))}
        if has_strain:
            r.update(strain=[{"pose": k, "strain": 5.0, "why": None}
                             for k in range(1, n + 1)],
                     strain_engine=engine, strain_version=strainmod.VERSION)
        if has_panel:
            r["panel"] = {"pose": 1, "bsa": 700.0}
        json.dump(r, open(os.path.join(rd, "results", "%s.json" % seq), "w"), indent=1)
        for k in range(1, n + 1):
            pd = os.path.join(rd, "mmgbsa", seq, "pose%d" % k)
            os.makedirs(pd)
            for t in FILES if files else ():
                open(os.path.join(pd, seq + t), "w").write("x\n")
    return rd


lacks = make("a_lacks", [("LLKWHY", "none", 5, False, False),
                         ("WHWHLLPF", "none", 2, False, False)])
hasall = make("b_has", [("GGWWGG", "none", 3, True, True)], ran=["amber", "openmm"])
ring = make("c_ring", [("CAAAAC", "disulfide", 3, True, True)])
unscored = make("d_unscored", [("KKKK", "none", 0, False, False)])
empty = make("e_empty", [("EEEE", "none", 0, False, False)], results=False)
gone = make("f_gone", [("FFFFF", "none", 4, False, False)], files=False)
outside = make("g_outside", [("YYYYY", "none", 1, False, True)])
# A run outside the runs folder, as a flow's run opened from HOPE Labs is.
flow = os.path.join(WORK, "flow_run")
shutil.move(outside, flow)

app = server.App(RUNS)
real_avail = (gbengines.OPENMM.available, gbengines.SCHRODINGER.available)
gbengines.OPENMM.available = lambda: True
gbengines.SCHRODINGER.available = lambda: True

# --- GET: what each run has ---------------------------------------------------
L = app.measures_list()
by = {r["name"]: r for r in L["runs"]}
check([m["k"] for m in L["measures"]] == ["strain", "panel", "openmm", "prime"],
      "the measures are not the four, in the dialog's order: %s" % L["measures"])
check([m["name"] for m in L["measures"]] == [
    "Peptide strain", "Interface measures", "MM-GBSA by OpenMM",
    "MM-GBSA by Schrödinger Prime"], "the measures' names are not the approved words")
a = by.get("a_lacks", {})
check(a.get("poses") == 7 and a.get("peptides") == 2 and not a.get("blocked")
      and a.get("has") == [] and a.get("missing") == {"strain": 2, "panel": 2,
                                                       "openmm": 2, "prime": 2}
      and a.get("created") == "2026-09-22 13:02:24",
      "a run lacking everything is not read as lacking everything: %s" % a)
b = by.get("b_has", {})
check(sorted(b.get("has") or []) == ["openmm", "panel", "strain"]
      and b.get("missing") == {"prime": 1},
      "a run with the strain, the panel and OpenMM is not read so: %s" % b)
c = by.get("c_ring", {})
check(c.get("cannot", {}).get("prime") == "ring peptide" and "prime" not in c.get("missing", {}),
      "Prime's refusal of a ring peptide is not the run's reason: %s" % c)
check(by.get("d_unscored", {}).get("blocked") == "no poses were scored"
      and by.get("d_unscored", {}).get("poses") == 0,
      "a run with no scored poses is not blocked with the reason: %s" % by.get("d_unscored"))
check(by.get("e_empty", {}).get("blocked") == "no peptide has finished docking yet",
      "a run with no results is not blocked with the reason: %s" % by.get("e_empty"))
f = by.get("f_gone", {})
check(f.get("missing", {}).get("strain") == 1 and not f.get("cannot"),
      "the cheap reading looked into the pose folders: %s" % f)
check("g_outside" not in by and not any(r.get("asked") for r in L["runs"]),
      "a run outside the runs folder is listed without being asked for")
L2 = app.measures_list(flow)
check(L2["runs"][0]["rundir"] == flow and L2["runs"][0].get("asked")
      and L2["runs"][0].get("missing") == {"strain": 1, "openmm": 1, "prime": 1},
      "the run the page has open, outside the runs folder, is not listed first: %s"
      % L2["runs"][:1])
L3 = app.measures_list(gone)
f3 = {r["name"]: r for r in L3["runs"]}["f_gone"]
check(f3.get("asked") and f3.get("cannot", {}).get("strain") == server.FILES_GONE
      and f3.get("cannot", {}).get("panel") == server.FILES_GONE,
      "the open run's measures whose files are gone are offered: %s" % f3)
gbengines.SCHRODINGER.available = lambda: False
check([m["k"] for m in app.measures_list()["measures"]] == ["strain", "panel", "openmm"],
      "an engine that is not installed here is offered")
check("cannot be added here" in app.measures_submit(
    {"runs": [lacks], "measures": ["prime"]}).get("error", ""),
      "an engine that is not installed here can be asked for")
gbengines.SCHRODINGER.available = lambda: True
st = app.run_status(lacks)
check((st.get("measures") or {}).get("missing", {}).get("strain") == 2
      and "offered" in st["measures"] and st["measures"].get("local") is False,
      "the run's own page is not told what it lacks: %s" % st.get("measures"))

# --- POST: the jobs each run gets, with the scheduler stood in for -------------
asked = []
ids = iter(range(5001, 6000))
real = (slurm.available, slurm.submit, slurm.state, slurm.local_mode,
        slurm.within_limit, slurm.cancel, slurm.submit_local)


def fake_submit(script, depends_on=None, kind="afterok"):
    jid = str(next(ids))
    asked.append({"script": script, "after": depends_on, "kind": kind, "id": jid})
    return jid


states = {}
slurm.available = lambda: True
slurm.submit = fake_submit
slurm.state = lambda jid, rundir=None: states.get(str(jid), {})
slurm.local_mode = lambda: False
slurm.within_limit = lambda t, partition="", extra="": t
slurm.cancel = lambda jid: ""
try:
    got = app.measures_submit({"runs": [lacks, hasall, ring, unscored, gone],
                               "measures": ["strain", "panel", "openmm", "prime"]})
    res = {os.path.basename(r["rundir"]): r for r in got.get("runs") or []}
    ra = res.get("a_lacks", {})
    jobs = ra.get("jobs") or {}
    check(sorted(ra.get("adds") or []) == ["openmm", "panel", "prime", "strain"]
          and sorted(jobs) == ["measures", "report", "rescore_openmm", "rescore_schrodinger"],
          "a run lacking everything did not get one short job, a job per engine and "
          "its report: %s" % ra)
    sub = {os.path.basename(x["script"]): x for x in asked
           if os.path.dirname(x["script"]) == lacks}
    m = open(os.path.join(lacks, "04_measures.sbatch")).read() \
        if os.path.isfile(os.path.join(lacks, "04_measures.sbatch")) else ""
    check("--cpus-per-task=16" in m and m.rstrip().endswith(" panel+strain")
          and "--time=%s" % sizing.hhmmss(sizing.measures_minutes(7, 16)) in m
          and "peptide strain" in m and "docked or rescored again" in m,
          "the short job is not min(16, the run's cores) for 2 x (poses x 1 min / "
          "cores) + 30 min, running the light stages:\n%s" % m)
    pr = open(os.path.join(lacks, "05_rescore_schrodinger.sbatch")).read()
    om = open(os.path.join(lacks, "05_rescore_openmm.sbatch")).read()
    check("--cpus-per-task=1\n" in pr and pr.rstrip().endswith("rescore:schrodinger")
          and "--cpus-per-task=16" in om and om.rstrip().endswith("rescore:openmm"),
          "each engine is not a job of its own, as a rescoring is")
    rep = sub.get("03_report.sbatch") or {}
    want_after = ":".join(jobs[k] for k in ("measures", "rescore_openmm",
                                            "rescore_schrodinger"))
    check(rep.get("after") == want_after and rep.get("kind") == "afterany",
          "the report does not wait on every one of the run's jobs, however they "
          "end: %s" % rep)
    tid = sub.get(tidy.SCRIPT) or {}
    check(tid.get("kind") == "afterany" and set((tid.get("after") or "").split(":"))
          >= set(jobs.values()) and (ra.get("tidy") or {}).get("armed"),
          "the run's own tidy is not queued again behind every new job: %s %s"
          % (tid, ra.get("tidy")))
    j = jobmod.load(lacks)
    check(all(j["jobs"].get(k) == v for k, v in jobs.items()),
          "the jobs are not recorded in job.json, where the page follows them")
    check([q["adds"] for q in ra.get("queued") or []] == [["strain", "panel"],
                                                          ["openmm"], ["prime"]]
          and ra.get("report") == jobs["report"],
          "the answer does not say what each job adds: %s" % ra.get("queued"))
    rb = res.get("b_has", {})
    check(rb.get("adds") == ["prime"] and sorted(rb.get("jobs") or {}) ==
          ["report", "rescore_schrodinger"]
          and not os.path.isfile(os.path.join(hasall, "04_measures.sbatch")),
          "what a run already has was queued again: %s" % rb)
    rc = res.get("c_ring", {})
    check(rc.get("adds") == ["openmm"],
          "a ring run was queued for Prime, which refuses it: %s" % rc)
    check(res.get("d_unscored", {}).get("error") == "no poses were scored"
          and not os.path.isfile(os.path.join(unscored, "03_report.sbatch")),
          "a run with no scored poses was queued, or not told why: %s" % res.get("d_unscored"))
    rf = res.get("f_gone", {})
    check(rf.get("adds") == ["openmm", "prime"],
          "a run whose files are gone was queued for the strain or the panel: %s" % rf)
    check(got.get("queued") == 4 and got.get("refused") == 1,
          "the totals are not the runs queued and refused: %s" % got)
    n = len(asked)
    again = app.measures_submit({"runs": [hasall], "measures": ["strain", "panel"]})
    check("already has" in (again["runs"][0].get("error") or "") and len(asked) == n,
          "a run that has what was chosen was queued: %s" % again)

    # The guards: a tidy moving the run's files, a job still going on it.
    rd = make("h_guard", [("GGGGG", "none", 2, False, False)])
    states["9999"] = {"RUNNING": 1}
    json.dump({"job_id": "9999", "undo": False}, open(os.path.join(rd, tidy.JOBFILE), "w"))
    said = app.measures_submit({"runs": [rd], "measures": ["strain"]})["runs"][0]
    check("tidied" in (said.get("error") or "") and len(asked) == n,
          "a run being tidied was queued: %s" % said)
    os.remove(os.path.join(rd, tidy.JOBFILE))
    j = jobmod.load(rd)
    j["jobs"] = {"dock": "8888"}
    jobmod.save(rd, j)
    states["8888"] = {"RUNNING": 3}
    said = app.measures_submit({"runs": [rd], "measures": ["strain"]})["runs"][0]
    check("still queued or running" in (said.get("error") or "") and "8888" in said["error"]
          and len(asked) == n,
          "a run with a job still running was queued: %s" % said)
    states["8888"] = {"COMPLETED": 3}
    said = app.measures_submit({"runs": [rd], "measures": ["strain"]})["runs"][0]
    check(sorted(said.get("jobs") or {}) == ["measures", "report"]
          and open(os.path.join(rd, "04_measures.sbatch")).read().rstrip().endswith(" strain"),
          "the strain alone is not one short job and its report: %s" % said)
    check("choose a measure" in app.measures_submit({"runs": [rd]}).get("error", "")
          and "choose a run" in app.measures_submit({"measures": ["strain"]}).get("error", ""),
          "an empty request is not refused with what is missing")

    # The local executor: the same stages, in turn, on this machine.
    ran_local = []
    slurm.local_mode = lambda: True
    slurm.submit_local = lambda rundir, jobfile, stages, cpus=1, python=None, env=None: (
        ran_local.append(stages) or {k: slurm.LOCAL_PREFIX + k for k, _t in stages})
    rl = make("i_local", [("IIIII", "none", 2, False, False)])
    said = app.measures_submit({"runs": [rl], "measures": ["strain", "openmm"]})["runs"][0]
    check(ran_local == [[("measures", ["strain"]), ("rescore_openmm", ["rescore:openmm"]),
                         ("report", ["report"])]]
          and said.get("jobs", {}).get("measures") == "L-measures"
          and not (said.get("tidy") or {}).get("armed"),
          "added measures do not run here, in turn, with the local executor: %s %s"
          % (ran_local, said))
    check(app.measures_list()["local"] is True,
          "the dialog is not told its jobs run on this computer")
finally:
    (slurm.available, slurm.submit, slurm.state, slurm.local_mode,
     slurm.within_limit, slurm.cancel, slurm.submit_local) = real

# --- the routes ------------------------------------------------------------------
httpd = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
httpd.app = app
httpd.token = app.token
threading.Thread(target=httpd.serve_forever, daemon=True).start()
opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
BASE = "http://127.0.0.1:%d" % httpd.server_address[1]
heard = []
real_list, real_submit = app.measures_list, app.measures_submit
app.measures_list = lambda also="": heard.append(("get", also)) or {"runs": []}
app.measures_submit = lambda body: heard.append(("post", body)) or {"runs": []}
try:
    with opener.open("%s/api/measures?dir=%s&t=%s" % (BASE, urllib.parse.quote(flow),
                                                       app.token), timeout=30) as r:
        body = json.loads(r.read())
    req = urllib.request.Request(BASE + "/api/measures",
                                 data=json.dumps({"runs": [lacks], "measures": ["strain"]}).encode(),
                                 headers={"Content-Type": "application/json",
                                          "X-Token": app.token})
    with opener.open(req, timeout=30) as r:
        body2 = json.loads(r.read())
    check(heard == [("get", flow), ("post", {"runs": [lacks], "measures": ["strain"]})]
          and body == {"runs": []} and body2 == {"runs": []},
          "GET and POST /api/measures do not reach the App: %s" % heard)
    try:
        opener.open(BASE + "/api/measures", timeout=30)
        check(False, "GET /api/measures answers without the token")
    except urllib.error.HTTPError as e:
        check(e.code == 403, "GET /api/measures without the token is %s" % e.code)
finally:
    app.measures_list, app.measures_submit = real_list, real_submit
    httpd.shutdown()
gbengines.OPENMM.available, gbengines.SCHRODINGER.available = real_avail

# --- the page ----------------------------------------------------------------------
html = page.render("TOKEN")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js"), encoding="utf-8").read()
css = html[html.index("<style>"):html.index("</style>")]
mcss = css[css.index("---- adding measures"):css.index("---- column widths")]
dlg = html[html.index("<dialog id=measuresdlg"):html.index("</dialog>",
                                                           html.index("<dialog id=measuresdlg"))]
flat = " ".join(dlg.split())
for want, why in (
        ("aria-labelledby=mstitle", "the dialog is not named by its heading"),
        ("<h2 id=mstitle>Add measures</h2>", "the dialog's heading is not Add measures"),
        ("id=mslead role=status aria-live=polite", "the lead line is not announced"),
        ("<label class=msall><input type=checkbox id=msmall> Select all</label>",
         "the measures' Select all is not a real checkbox with its label"),
        ("<label class=msall><input type=checkbox id=msrall> <span id=msralllab>Select all</span></label>",
         "the runs' Select all is not a real checkbox with its label"),
        ('aria-label="Which runs to list"', "the filter is not a named group"),
        ("id=msfmissing aria-pressed=true>Missing something</button>",
         "the filter has no Missing something, pressed"),
        ("id=msfall aria-pressed=false>Every run</button>", "the filter has no Every run"),
        ('type=search placeholder="Find a run" aria-label="Find a run"',
         "the search is not Find a run, labeled"),
        ("<th>Run</th>", "the runs table has no Run column"),
        ("<th class=msnum>Poses</th><th>Adds</th>", "the runs table has no Poses and Adds"),
        ("Nothing already measured is changed. Shift-click to choose a range of runs.",
         "the footer does not say nothing measured is changed"),
        ("class=ghost id=mscancel>Cancel</button>", "there is no Cancel"),
        ("class=primary id=msgo disabled>Add</button>", "Add is not the dialog's one primary"),
        ("aria-label=Close", "the close button has no name")):
    check(want in flat, why)
check("class=secondary id=measuresopen>Add measures&hellip;</button>" in html,
      "the Runs card has no Add measures button")
small = [m.group(0) for m in re.finditer(r"font-size:\s*(\d+(?:\.\d+)?)px", mcss)
         if float(m.group(1)) < 13]
check(not small, "the dialog sets text under 13 px: %s" % small)
check("var(--accent)" in mcss and "var(--warn)" in mcss and "var(--good)" in mcss
      and "#0d5c6b" not in mcss and "--hot" not in mcss,
      "the dialog is not drawn in this page's own tokens")
# The controls it replaces are gone, and nothing reads what they needed.
for gone_ in ("panelgo", "panelmsg", "rescoreControl", "rescorego", "rescoremsg",
              "rescoreeng", "GBAVAIL", "Measure them"):
    check(gone_ not in appjs, "%s is still in the page's script" % gone_)
check("api('/api/panel'" not in appjs and "api('/api/rescore'" not in appjs,
      "the page still asks the old endpoints")
ms = appjs[appjs.index("/* ===================== adding measures"):
           appjs.index("/* ===================== the 3D panel")]
for want, why in (("box.indeterminate = on > 0 && on < of", "Select all has no third state"),
                  ("dlg.showModal()", "the dialog is not modal, so Esc does not close it"),
                  ("first.focus()", "the dialog does not take the focus"),
                  ("ask('/api/measures' + (at", "the dialog does not read the runs"),
                  ("ask('/api/measures', {runs:", "Add does not post the choice"),
                  ("e.shiftKey", "shift-click does not choose a range"),
                  ("aria-pressed", "the filter does not say which is pressed"),
                  ("$('measuresopen').onclick", "the Runs card's button opens nothing")):
    check(want in ms, why)
check("const ma = $('msadd')" in appjs and "openMeasures({run: dir" in appjs
      and "+ measuresLine(d.measures, d.results, " in appjs,
      "the line over a run's table is not drawn, or does not open the dialog on it")
check(not re.search(r"\b(colour|minimis|optimis|analys|behaviour|favour|licence|"
                    r"recognis|organis|centre)", ms + flat, re.I),
      "the dialog's words are not US English")

# --- the rules, under Node -----------------------------------------------------------
NODE = os.environ.get("ADCP_NODE") or shutil.which("node")
if NODE and os.path.isfile(NODE if os.sep in NODE else shutil.which(NODE) or ""):
    names = ("msLacks", "msAdds", "msSelectable", "msShown", "msName", "msShort",
             "msLower", "msWords", "msCount", "msPlan", "msLead", "msTicked", "msChips",
             "msQueued", "measuresLine", "msWhen")
    parts = [re.search(r"const esc = .*?\}\[c\]\)\);", appjs, re.S).group(0),
             re.search(r"const MSLIGHT = [^;]*;", appjs).group(0)]
    for nm in names:
        # A function on one line ends with it; any other at its closing brace.
        i = appjs.index("function %s(" % nm)
        line = appjs[i:appjs.index("\n", i)]
        parts.append(line if line.rstrip().endswith("}")
                     else appjs[i:appjs.index("\n}\n", i) + 3])
    LJ = {"local": False,
          "measures": [{"k": k, "name": n, "what": "", "job": "light" if k in ("strain", "panel")
                        else "engine"}
                       for k, n in (("strain", "Peptide strain"), ("panel", "Interface measures"),
                                    ("openmm", "MM-GBSA by OpenMM"),
                                    ("prime", "MM-GBSA by Schrödinger Prime"))],
          "runs": [{"rundir": "/r/a", "name": "round2", "poses": 35, "has": ["panel"],
                    "cannot": {}},
                   {"rundir": "/r/b", "name": "peptides6B88", "poses": 5, "has": [],
                    "cannot": {}},
                   {"rundir": "/r/c", "name": "drop10", "poses": 3, "has": ["panel"],
                    "cannot": {"prime": "ring peptide"}},
                   {"rundir": "/r/d", "name": "test6", "poses": 2,
                    "has": ["strain", "panel"], "cannot": {}},
                   {"rundir": "/r/e", "name": "empty", "poses": 0, "has": [], "cannot": {},
                    "blocked": "no poses were scored"}]}
    prog = ("\n".join(parts) + r"""
const L = JSON.parse(process.argv[1]);
const S = (m, r, extra) => Object.assign({m: new Set(m), r: new Set(r), filter: 'missing',
                                          find: ''}, extra || {});
const out = {};
out.tick = msLead(S([], []), L);
out.choose = msLead(S(['strain', 'panel'], []), L);
out.ready = msLead(S(['panel', 'strain'], ['/r/a', '/r/b']), L);
out.heavy = msLead(S(['strain', 'openmm'], ['/r/a']), L);
out.local = msLead(S(['strain'], ['/r/a']), Object.assign({}, L, {local: true}));
out.plan = msPlan(S(['strain', 'panel', 'openmm', 'prime'], ['/r/a', '/r/c', '/r/e']), L);
out.plan.runs = out.plan.runs.map(r => r.name);
out.shown = msShown(S(['strain'], []), L.runs).map(r => r.name);
out.every = msShown(S(['strain'], [], {filter: 'all'}), L.runs).map(r => r.name);
out.found = msShown(S(['strain'], [], {filter: 'all', find: 'TEST'}), L.runs).map(r => r.name);
out.chips = L.runs.map(r => msChips(S(['strain', 'prime'], []), L, r));
out.none = msChips(S([], []), L, L.runs[0]);
out.had = msChips(S(['panel'], []), L, L.runs[3]);
out.line2 = measuresLine({peptides: 5, missing: {strain: 5, panel: 5, openmm: 5}, has: [],
                          cannot: {}, offered: ['strain', 'panel', 'openmm']}, {});
out.line1 = measuresLine({peptides: 5, missing: {strain: 5}, has: ['panel'], cannot: {},
                          local: true}, {});
out.linep = measuresLine({peptides: 5, missing: {panel: 2}, has: ['strain'], cannot: {}}, {});
out.lineo = measuresLine({peptides: 5, missing: {openmm: 5}, has: ['strain', 'panel'],
                          cannot: {}}, {});
out.lineb = measuresLine({blocked: 'no poses were scored'}, {});
out.linebusy = measuresLine({peptides: 5, missing: {strain: 5}, has: [], cannot: {}}, {}, true);
out.lineg = measuresLine({peptides: 3, missing: {}, has: ['strain'],
                          cannot: {panel: 'gone'}}, {panel_missing: 3});
out.queued = msQueued(L, {runs: [
  {name: 'round2', queued: [{id: '101', adds: ['strain', 'panel']},
                            {id: '102', adds: ['openmm']}], report: '103'},
  {name: 'test6', queued: [{id: '104', adds: ['prime']}], report: '105'},
  {name: 'empty', error: 'no poses were scored'}]});
out.nothing = msQueued(L, {runs: [{name: 'empty', error: 'no poses were scored'}]});
const y = new Date().getFullYear();
out.when = [msWhen(y + '-09-22 13:02:24'), msWhen('2024-01-05 08:00:00'), msWhen('')];
process.stdout.write(JSON.stringify(out));
""")
    run = subprocess.run([NODE, "-e", prog, json.dumps(LJ)], capture_output=True,
                         text=True, timeout=60)
    check(run.returncode == 0, "the dialog's rules could not be run under Node: %s"
          % run.stderr[-600:])
    o = json.loads(run.stdout or "{}") if run.returncode == 0 else {}
    if o:
        check(not o["tick"]["ready"] and o["tick"]["html"].startswith("<b>Tick what to add.</b>")
              and "Each measure says how many of your runs lack it." in o["tick"]["html"],
              "with nothing ticked the line does not say to tick: %s" % o["tick"])
        check(o["choose"]["html"].startswith("<b>Now choose the runs.</b>")
              and "3 of them lack peptide strain and interface measures. Select all takes "
                  "every one listed." in o["choose"]["html"],
              "with nothing chosen the line does not say to choose the runs: %s" % o["choose"])
        check(o["ready"]["ready"] and "<b>Peptide strain and Interface measures for 2 runs "
              "(40 poses).</b>" in o["ready"]["html"]
              and "2 jobs on CPU nodes, charged to your usual account. Short jobs."
              in o["ready"]["html"],
              "the ready line does not say what Add will queue: %s" % o["ready"])
        check("An MM-GBSA engine takes about as long as the run’s own MM-GBSA did."
              in o["heavy"]["html"] and "2 jobs" in o["heavy"]["html"],
              "an engine is not said to take as long as the run's own MM-GBSA: %s" % o["heavy"])
        check("1 job on this computer, one run after another." in o["local"]["html"]
              and "(35 poses)" in o["local"]["html"],
              "the local executor's line does not say the jobs run here: %s" % o["local"])
        check(o["plan"] == {"runs": ["round2", "drop10"], "jobs": 5, "poses": 38},
              "Add would not queue one short job a run and a job per engine it lacks, "
              "and nothing for a blocked run: %s" % o["plan"])
        check(o["shown"] == ["round2", "peptides6B88", "drop10"]
              and o["every"] == ["round2", "peptides6B88", "drop10", "test6", "empty"]
              and o["found"] == ["test6"],
              "Missing something, Every run or Find a run list the wrong runs: %s %s %s"
              % (o["shown"], o["every"], o["found"]))
        ch = o["chips"]
        check("<span class=mschip>Peptide strain</span>" in ch[0]
              and "<span class=mschip>Schrödinger Prime</span>" in ch[0]
              and "no Schrödinger Prime: ring peptide" in ch[2]
              and "has what you ticked" not in ch[3] and "Schrödinger Prime" in ch[3]
              and "no poses were scored" in ch[4] and "tick a measure" in o["none"]
              and "has what you ticked" in o["had"],
              "the Adds column does not say what each run gets: %s %s %s" % (ch, o["none"],
                                                                            o["had"]))
        check("<b>This run has no peptide strain or interface measures.</b> One short job "
              "on a CPU node adds both." in o["line2"] and ">Add them…</button>" in o["line2"]
              and 'data-ms="strain panel"' in o["line2"] and "OpenMM" not in o["line2"],
              "the line over a run's table does not say what it lacks: %s" % o["line2"])
        check("<b>This run has no peptide strain.</b> One short job on this computer adds it."
              in o["line1"] and ">Add it…</button>" in o["line1"],
              "the line for the strain alone is not adapted to it: %s" % o["line1"])
        check("2 of this run’s 5 scored peptides have no interface measures." in o["linep"]
              and "adds them" in o["linep"],
              "the line does not count the peptides that lack a measure: %s" % o["linep"])
        check(o["lineo"] == "" and o["lineb"] == "" and o["linebusy"] == "",
              "a line stands over a run that lacks only an engine, is blocked, or has "
              "its jobs running: %s %s %s" % (o["lineo"], o["lineb"], o["linebusy"]))
        check("3 peptide(s) have no interface columns" in o["lineg"] and "msadd" not in o["lineg"],
              "a run whose interface files are gone is offered the interface: %s" % o["lineg"])
        q = o["queued"]
        check(q["ready"] and q["lead"].startswith("<b>Queued 3 jobs for 2 runs; 1 run could not be.</b>")
              and "job 101, peptide strain and interface measures; job 102, MM-GBSA by OpenMM; "
                  "then its report, job 103" in q["body"]
              and "<span class=msrun>empty</span>: not queued, no poses were scored" in q["body"],
              "after Add the dialog does not say what was queued and what was not: %s" % q)
        check(not o["nothing"]["ready"] and "Nothing was queued." in o["nothing"]["lead"],
              "a press that queued nothing does not say so: %s" % o["nothing"])
        check(o["when"][0] == "22 Sep" and o["when"][1] == "5 Jan 2024" and o["when"][2] == "",
              "a run's date is not shown as 22 Sep, with its year when not this one: %s"
              % o["when"])
else:
    print("note  no Node here (ADCP_NODE or node on the PATH), so the dialog's rules "
          "were not run")

shutil.rmtree(WORK, ignore_errors=True)
if fails:
    print("test_measures FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  every run's measures are read cheaply, each run gets only what it lacks, "
      "refused while busy, and the dialog says what it will do, %d checks" % ran)
