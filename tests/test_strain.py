"""The peptide strain: the row, the table files, a finished run brought up to
date, a new run measuring it, and what the run says it did.

By Aadhil Haq

_vendor/strain.py measures one pose (its own tests are HOPE-pipelines'; the
copy is held to hope_core's by tests/test_vendored_strain.py). What is checked
here is everything around it, with the measurement stood in for:

  * the row's Peptide strain is the strain of the pose with the lowest MM-GBSA
    energy among those scored, the pose the row's other per pose numbers
    describe; a failure is an empty cell with its reason, never 0; a run Prime
    scored shows Prime's own; the strain is in no ranking, combined score or
    z-score;
  * results.tsv, results.json and results.txt carry it, right after the
    MM-GBSA columns, under the approved heading and help;
  * stage_strain measures only the scored poses of the peptides that have no
    strain yet, finds each pose's folder by its energy, writes the strain into
    that peptide's result file and changes nothing else in it, and tries again
    a pose whose measurement could not start; a Prime run computes nothing;
  * a new run measures it behind mmgbsa_strain (on by default, off with
    --no-strain), and a job file written before it does not;
  * the worker's stages, the command line, the estimate, the job script's
    notes, the citation and the report's methods say so.

    python3 tests/test_strain.py
"""
import contextlib
import csv
import io
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

try:
    from adcp_pilot import cite, cli, job as jobmod, pipeline, report, sizing, worker  # noqa: E402
    from adcp_pilot._vendor import strain as strainmod  # noqa: E402
except ImportError as exc:
    print("skip: the package cannot be imported here (%s)" % exc)
    sys.exit(77)

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


WORK = tempfile.mkdtemp(prefix="adcp_strain_")

# --- the row -------------------------------------------------------------------
MM = [-30.0, -35.5, -28.0]
RANKS = [1, 2, 4]
S = [{"pose": 1, "strain": 21.5, "why": None},
     {"pose": 2, "strain": 18.25, "why": None},
     {"pose": 4, "strain": 40.0, "why": None}]
row = report.summarise("LLKWHY", [{"score": -14.0}], mm=MM, mm_ranks=RANKS,
                       strain=S, strain_engine="amber")
check(row.get("pep_strain") == 18.25 and row.get("pep_strain_pose") == 2,
      "the row's strain is not that of the lowest energy pose (ADCP pose 2): %s"
      % {k: v for k, v in row.items() if k.startswith("pep_")})
check("pep_strain_why" not in row, "a measured strain carries a reason: %s" % row)
bad = [dict(s) for s in S]
bad[1].update(strain=None, why="sander did not finish the relaxation")
row2 = report.summarise("LLKWHY", [{"score": -14.0}], mm=MM, mm_ranks=RANKS,
                        strain=bad)
check(row2.get("pep_strain") is None and "pep_strain" in row2
      and row2.get("pep_strain_why") == "sander did not finish the relaxation",
      "a failed strain is not an empty cell with its reason: %s" % row2)
check(row2.get("pep_strain") != 0, "a failed strain reads as 0")
row3 = report.summarise("LLKWHY", [{"score": -14.0}], mm=MM, mm_ranks=RANKS,
                        strain=S[:1] + S[2:])
check(row3.get("pep_strain") is None and row3.get("pep_strain_why"),
      "the row borrowed another pose's strain when its own pose has none: %s" % row3)
tie = report.summarise("X", [], mm=[-30.0, -30.0], mm_ranks=[3, 5],
                       strain=[{"pose": 3, "strain": 1.0}, {"pose": 5, "strain": 2.0}])
check(tie.get("pep_strain") == 1.0, "of two equal energies the first pose is not taken")
legacy = report.summarise("X", [], mm=[-20.0, -25.0],
                          strain=[{"pose": 1, "strain": 3.0}, {"pose": 2, "strain": 4.0}])
check(legacy.get("pep_strain") == 4.0,
      "a result with no ranks does not take the values' own positions: %s" % legacy)
none = report.summarise("X", [], mm=MM, mm_ranks=RANKS)
check("pep_strain" not in none, "a peptide never measured has a strain cell: %s" % none)
prime = report.summarise("X", [], mm=MM, mm_ranks=RANKS, strain=S,
                         strain_engine="schrodinger")
check(prime.get("pep_strain_engine") == "schrodinger",
      "a Prime run's row does not say its strain is Prime's")
check(report.strain_help([prime]).endswith(
    "This run was scored by Schrödinger Prime, so this is Prime's own ligand strain.")
      and not report.strain_help([row]).endswith("ligand strain."),
      "the help does not add the Prime sentence for a Prime run, and only then")
check(report.STRAIN_HELP.startswith("What the peptide pays to hold its docked shape: "
                                    "its energy in the pose, minus the same peptide")
      and report.STRAIN_HELP.endswith("It is not part of the MM-GBSA ranking.")
      and report.STRAIN_HEADING == "Peptide strain",
      "the heading or the help is not the approved words")

# The page carries the same words (app.js), word for word.
appjs = open(os.path.join(ROOT, "adcp_pilot", "static", "app.js"), encoding="utf-8").read()


def js_string(name):
    """A const string assembled from '...' + '...' pieces in app.js."""
    i = appjs.index("const %s = " % name)
    j = appjs.index(";\n", i)
    parts = []
    import re
    for m in re.finditer(r"'((?:[^'\\]|\\.)*)'", appjs[i:j]):
        parts.append(m.group(1).encode().decode("unicode_escape"))
    return "".join(parts)


check(js_string("STRAIN_HELP") == report.STRAIN_HELP
      and js_string("STRAIN_HELP_PRIME") == report.STRAIN_HELP_PRIME
      and js_string("STRAIN_HEADING") == report.STRAIN_HEADING,
      "the page's help for the strain is not report.STRAIN_HELP: %r"
      % js_string("STRAIN_HELP"))

# Never ranked on: the order, the combined score and its z-scores are the same
# whatever the strain says.
rows = []
for i, (pep, mm, bsa, hb) in enumerate((("AAAA", [-30.0, -31.0], 900.0, 4),
                                        ("CCCC", [-35.0, -20.0], 700.0, 6),
                                        ("DDDD", [-25.0, -27.0], 800.0, 2))):
    r = report.summarise(pep, [{"score": -10.0 - i}], mm=mm,
                         panel={"pose": 1, "bsa": bsa, "bsa_per_res": bsa / 4,
                                "n_hbond": hb, "site_dg": -10.0 - i})
    rows.append(r)
plain = json.loads(json.dumps(rows))
strained = json.loads(json.dumps(rows))
for r, v in zip(strained, (5.0, 99.0, 0.5)):
    r["pep_strain"] = v
report.combined_scores(plain)
report.combined_scores(strained)
check([r.get("combined_z") for r in plain] == [r.get("combined_z") for r in strained]
      and all(not k.startswith("z_pep") for r in strained for k in r),
      "the strain moves the combined score or has a z of its own")
check([r["peptide"] for r in report.rank(plain)[0]]
      == [r["peptide"] for r in report.rank(strained)[0]],
      "the strain changes the order")

# --- the table files -----------------------------------------------------------
out = os.path.join(WORK, "files")
rows = [report.summarise("LLKWHY", [{"score": -14.0}], mm=MM, mm_ranks=RANKS,
                         strain=S, mmgbsa_engine="amber"),
        report.summarise("GGGG", [{"score": -9.0}], mm=[-20.0], mm_ranks=[1],
                         strain=[{"pose": 1, "strain": None, "why": "x"}],
                         mmgbsa_engine="amber")]
report.write(rows, out)
head = open(os.path.join(out, "results.tsv")).read().splitlines()
cols = head[0].split("\t")
check("pep_strain" in cols and cols.index("pep_strain") == cols.index("mmgbsa_engine") + 1,
      "results.tsv does not carry pep_strain right after the MM-GBSA columns: %s" % cols)
cells = {line.split("\t")[cols.index("peptide")]: line.split("\t")[cols.index("pep_strain")]
         for line in head[1:]}
check(cells.get("LLKWHY") == "18.25" and cells.get("GGGG") == "",
      "results.tsv's strain cells are not the value and an empty cell: %s" % cells)
doc = json.load(open(os.path.join(out, "results.json")))
got = {r["peptide"]: r for r in doc["rows"]}
check(got["LLKWHY"].get("pep_strain") == 18.25 and got["GGGG"].get("pep_strain") is None
      and got["GGGG"].get("pep_strain_why") == "x",
      "results.json's rows do not carry the strain and the reason")
txt = open(os.path.join(out, "results.txt")).read()
check("Peptide strain" in txt and "18.25" in txt,
      "results.txt has no Peptide strain column:\n%s" % txt[:600])
report.write([report.summarise("LLKWHY", [{"score": -14.0}], mm=MM)], out)
check("pep_strain" not in open(os.path.join(out, "results.tsv")).readline()
      and "Peptide strain" not in open(os.path.join(out, "results.txt")).read(),
      "a run with no strain measured has an empty strain column")

# --- a finished run brought up to date -----------------------------------------
DAT = """
Differences (Complex - Receptor - Ligand):
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
VDWAALS                    -40.0000                0.0000              0.0000
EEL                         -5.0000                0.0000              0.0000
EGB                         20.0000                0.0000              0.0000
ESURF                       -5.0000                0.0000              0.0000

DELTA TOTAL                %s                0.0000              0.0000
"""


def make_run(name, engine="amber", strain_key=None):
    rd = os.path.join(WORK, name)
    j = jobmod.default(mmgbsa=True, mmgbsa_engine=engine, label=name,
                       peptides=[{"name": "a", "sequence": "AAAAA"}])
    if strain_key is None:
        j.pop("mmgbsa_strain")              # a job file from before strain
    else:
        j["mmgbsa_strain"] = strain_key
    jobmod.save(rd, j)
    os.makedirs(os.path.join(rd, "results"))
    mm = os.path.join(rd, "mmgbsa")
    # AAAAA: four poses given to MM-GBSA, the third discarded (its value is in
    # the cache and not in the result), so ADCP poses 1, 2 and 4 are scored,
    # in folders 1, 2 and 4.
    for n in (1, 2, 3, 4):
        os.makedirs(os.path.join(mm, "AAAAA", "pose%d" % n))
    with open(os.path.join(mm, "AAAAA" + (".dg" if engine == "amber" else ".%s.dg" % engine)),
              "w") as fh:
        fh.write("# igb 5 saltcon 0.150 intdiel 1 radii mbondi2\n"
                 "1 -30.000000\n2 -35.500000\n3 -99.000000\n4 -28.000000\n")
    a = {"name": "a", "peptide": "AAAAA", "engine": "ADCP 1.1.21",
         "modes": [{"rank": k, "score": -10.0 - k} for k in range(1, 6)],
         "mmgbsa": MM, "mmgbsa_ranks": RANKS, "mmgbsa_engine": engine,
         "panel": {"pose": 2, "bsa": 800.0}, "note": "café Å"}
    # CCCCC: measured already; GGGGG: never rescored; DDDDD: written before
    # ranks were recorded, its energies only in each folder's own output.
    c = {"peptide": "CCCCC", "mmgbsa": [-20.0], "mmgbsa_ranks": [1],
         "strain": [{"pose": 1, "strain": 9.0, "why": None}],
         "strain_engine": "amber", "strain_version": strainmod.VERSION}
    g = {"peptide": "GGGGG", "mmgbsa": [], "mmgbsa_skipped": "non-standard residues"}
    dd = {"peptide": "DDDDD", "mmgbsa": [-22.0, -26.5]}
    for n, v in ((1, "-22.0000"), (2, "-26.5000")):
        p = os.path.join(mm, "DDDDD", "pose%d" % n)
        os.makedirs(p)
        open(os.path.join(p, "DDDDD_mmgbsa.dat"), "w").write(DAT % v)
    files = {}
    for r in (a, c, g, dd):
        f = os.path.join(rd, "results", "%s.json" % r["peptide"])
        with open(f, "w") as fh:
            json.dump(r, fh, indent=1)
        files[r["peptide"]] = f
    return rd, files


calls = []
answers = {}
cached = {}
real = (strainmod.measure, strainmod.read, strainmod.measurable)


def fake_measure(posedir, prefix=None, home=None, env=None, log=None, force=False):
    calls.append(os.path.relpath(posedir, WORK))
    got = answers.get(os.path.basename(posedir), {"strain": 12.345, "why": None})
    if got.get("strain") is not None or got.get("cached"):
        cached[posedir] = got
    return dict(got)


def fake_read(posedir, prefix=None):
    return cached.get(posedir)


def fake_measurable(posedir, prefix=None):
    return True, None


strainmod.measure, strainmod.read, strainmod.measurable = (fake_measure, fake_read,
                                                            fake_measurable)
saved_env = os.environ.get("SLURM_CPUS_PER_TASK")
os.environ["SLURM_CPUS_PER_TASK"] = "2"
try:
    rd, files = make_run("old")
    before = {k: open(f).read() for k, f in files.items()}
    said = []
    n = pipeline.stage_strain(rd, log=said.append)
    measured = sorted(calls)
    check(n == 2 and measured == sorted([os.path.join("old", "mmgbsa", "AAAAA", "pose1"),
                                         os.path.join("old", "mmgbsa", "AAAAA", "pose2"),
                                         os.path.join("old", "mmgbsa", "AAAAA", "pose4"),
                                         os.path.join("old", "mmgbsa", "DDDDD", "pose1"),
                                         os.path.join("old", "mmgbsa", "DDDDD", "pose2")]),
          "stage_strain did not measure exactly the scored poses still missing it "
          "(%s): %s" % (n, measured))
    new = json.load(open(files["AAAAA"]))
    check([s["pose"] for s in new["strain"]] == [1, 2, 4]
          and all(s["strain"] == 12.35 for s in new["strain"])
          and new["strain_engine"] == "amber" and new["strain_version"] == strainmod.VERSION,
          "AAAAA's strain is not one entry per scored pose, by its ADCP rank: %s"
          % new.get("strain"))
    rest = {k: v for k, v in new.items()
            if k not in ("strain", "strain_engine", "strain_version")}
    check(json.dumps(rest, indent=1) == before["AAAAA"],
          "stage_strain changed something else in the result file")
    check(open(files["CCCCC"]).read() == before["CCCCC"]
          and open(files["GGGGG"]).read() == before["GGGGG"],
          "a peptide measured already, or never rescored, had its file rewritten")
    dd = json.load(open(files["DDDDD"]))
    check([s["pose"] for s in dd["strain"]] == [1, 2],
          "a result with no ranks, its energies only in each folder, was not "
          "matched to its folders: %s" % dd.get("strain"))
    check(jobmod.load(rd).get("mmgbsa_strain") is True,
          "the run does not record that it measures strain now, so a resume would not")
    calls[:] = []
    check(pipeline.stage_strain(rd, log=said.append) == 0 and not calls,
          "a second stage_strain measured again: %s" % calls)
    # The report shows the column, the lowest energy pose's, and the rest of
    # the table as before.
    pipeline.stage_report(rd, log=lambda *a: None)
    tsv = open(os.path.join(rd, "results.tsv")).read().splitlines()
    c = tsv[0].split("\t")
    vals = {l.split("\t")[c.index("peptide")]: l.split("\t")[c.index("pep_strain")]
            for l in tsv[1:]}
    check(vals.get("AAAAA") == "12.35" and vals.get("GGGGG") == "",
          "the rebuilt report does not show the strain: %s" % vals)

    # A pose whose measurement could not start is tried again; one that
    # could not be matched to its energy says so and is not.
    rd2, files2 = make_run("retry")
    answers.clear()
    answers["pose2"] = {"strain": None, "why": "AmberTools' sander was not found"}
    calls[:] = []
    pipeline.stage_strain(rd2, log=lambda *a: None)
    a2 = json.load(open(files2["AAAAA"]))
    e2 = {s["pose"]: s for s in a2["strain"]}
    check(e2[2].get("retry") and e2[2]["strain"] is None and not pipeline.strain_done(a2),
          "a pose whose measurement could not start is not left to try again: %s" % e2)
    answers.clear()
    calls[:] = []
    pipeline.stage_strain(rd2, log=lambda *a: None)
    check(os.path.join("retry", "mmgbsa", "AAAAA", "pose2") in calls
          and pipeline.strain_done(json.load(open(files2["AAAAA"]))),
          "it was not tried again: %s" % calls)
    r = json.load(open(files2["DDDDD"]))
    r.pop("strain")
    r["mmgbsa"] = [-22.0, -77.0]
    json.dump(r, open(files2["DDDDD"], "w"), indent=1)
    pipeline.stage_strain(rd2, log=lambda *a: None)
    dd2 = {s["pose"]: s for s in json.load(open(files2["DDDDD"]))["strain"]}
    check(dd2[1]["strain"] is not None and dd2[2]["strain"] is None
          and dd2[2]["why"] == pipeline.UNMATCHED and not dd2[2].get("retry"),
          "a pose matching no folder is not said to be unmatched: %s" % dd2)

    # A run Prime scored: Prime's own strain, read, and nothing measured.
    rd3, files3 = make_run("prime", engine="schrodinger", strain_key=True)
    for n, v in ((1, "11.5"), (2, "7.25"), (4, "")):
        with open(os.path.join(rd3, "mmgbsa", "AAAAA", "pose%d" % n, "AAAAA-out.csv"),
                  "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["title", "r_psp_MMGBSA_dG_Bind(NS)", strainmod.PRIME_KEY])
            w.writerow(["x", "-60.0", v])
    calls[:] = []
    pipeline.stage_strain(rd3, log=lambda *a: None)
    p3 = json.load(open(files3["AAAAA"]))
    e3 = {s["pose"]: s for s in p3["strain"]}
    check(not [c_ for c_ in calls if "prime" in c_] and p3["strain_engine"] == "schrodinger"
          and e3[1]["strain"] == 11.5 and e3[2]["strain"] == 7.25
          and e3[4]["strain"] is None and e3[4]["why"],
          "a Prime run's strain is not Prime's own, read and not computed: %s %s"
          % (calls, e3))
    pipeline.stage_report(rd3, log=lambda *a: None)
    pr = {r["peptide"]: r for r in json.load(open(os.path.join(rd3, "results.json")))["rows"]}
    check(pr["AAAAA"].get("pep_strain") == 7.25
          and pr["AAAAA"].get("pep_strain_engine") == "schrodinger",
          "a Prime run's row does not show Prime's strain of its lowest energy pose: %s"
          % pr["AAAAA"])

    # --- a new run ----------------------------------------------------------------
    import adcp_pilot.pipeline as P
    saved = (P._mm_inputs, P._rec_for_mmgbsa, P.mmgbsa.run, P.verify.mmgbsa_split,
             P._panel_scored)
    P._mm_inputs = lambda j, seq, pdir, n, log=print, ring=None: (["m1", "m2", "m3"],
                                                                  [1, 2, 3], [])
    P._rec_for_mmgbsa = lambda j, log=print: None

    def fake_run(constructs, rec, poses, workdir, log=print, details=None, **kw):
        details["AAAAA"] = {"values": [-30.0, -35.0, -28.0], "ranks": [1, 2, 3],
                            "radii": "mbondi2", "model": "obc2"}
        return {"AAAAA": -31.0}
    P.mmgbsa.run = fake_run
    P.verify.mmgbsa_split = lambda dat: (True, None)
    P._panel_scored = lambda *a, **k: None
    try:
        for want in (True, False, None):
            j = jobmod.default(mmgbsa=True, cpus=2,
                               peptides=[{"name": "a", "sequence": "AAAAA"}])
            if want is None:
                j.pop("mmgbsa_strain")
            else:
                j["mmgbsa_strain"] = want
            o = {"modes": [{"rank": 1}]}
            calls[:] = []
            rdn = os.path.join(WORK, "new_%s" % want)
            os.makedirs(rdn, exist_ok=True)
            P._stage_mmgbsa(rdn, j, "AAAAA", o, os.path.join(rdn, "docking"), 10, 2,
                            {"cyclic": "none"}, log=lambda *a: None)
            if want:
                check([s["pose"] for s in o.get("strain") or []] == [1, 2, 3]
                      and o.get("strain_engine") == "amber" and len(calls) == 3,
                      "a new run does not measure the strain of its scored poses: %s %s"
                      % (o.get("strain"), calls))
            else:
                check("strain" not in o and not calls,
                      "a run with mmgbsa_strain %s measured strain anyway" % want)
    finally:
        (P._mm_inputs, P._rec_for_mmgbsa, P.mmgbsa.run, P.verify.mmgbsa_split,
         P._panel_scored) = saved
finally:
    strainmod.measure, strainmod.read, strainmod.measurable = real
    if saved_env is None:
        os.environ.pop("SLURM_CPUS_PER_TASK", None)
    else:
        os.environ["SLURM_CPUS_PER_TASK"] = saved_env

check(jobmod.default()["mmgbsa_strain"] is True,
      "a new run's job file does not measure strain by default")

# --- the command line, the worker, the estimate, the notes ---------------------
def run_cli(*argv):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        try:
            rc = cli.main(list(argv))
        except SystemExit as e:
            rc = e.code
    return rc, buf.getvalue()


REC = os.path.join(ROOT, "adcp_pilot", "example_run", "receptor", "1YCR.pdb")
for flag, want in (([], True), (["--no-strain"], False)):
    rd = os.path.join(WORK, "cli_%s" % want)
    rc, said = run_cli("new", rd, "-r", REC, "--chain", "A", "-s", "A:54",
                       "-p", "p53: ETFSDLWKLLPEN", "--mmgbsa", *flag)
    check(os.path.isfile(jobmod.path_of(rd))
          and jobmod.load(rd).get("mmgbsa_strain") is want,
          "adcp-pilot new %s does not write mmgbsa_strain %s (rc %r): %s"
          % (" ".join(flag), want, rc, said[-300:]))
rc, said = run_cli("strain", "--help")
check(rc == 0 and "--here" in said, "the command line has no strain subcommand: %s" % said)

heard = []
real_stages = (pipeline.stage_panel, pipeline.stage_strain)
pipeline.stage_panel = lambda rd, j=None, log=print: heard.append("panel")
pipeline.stage_strain = lambda rd, j=None, log=print: heard.append("strain")
try:
    jf = jobmod.path_of(os.path.join(WORK, "old"))
    check(worker.main([jf, "strain"]) == 0 and heard == ["strain"],
          "the worker has no strain stage: %s" % heard)
    heard[:] = []
    check(worker.main([jf, "panel+strain"]) == 0 and heard == ["panel", "strain"],
          "the worker does not run the light stages in turn: %s" % heard)
    heard[:] = []

    def stops(rd, j=None, log=print):
        raise RuntimeError("no panel")
    pipeline.stage_panel = stops
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        rc = worker.main([jf, "panel+strain"])
    check(rc == 1 and heard == ["strain"],
          "one light stage that stops keeps the next from running, or the job "
          "does not fail: %s %s" % (rc, heard))
finally:
    pipeline.stage_panel, pipeline.stage_strain = real_stages

base = dict(poses=10)
p0 = sizing.plan([9], 50, 2500000, 16, mmgbsa=True, mm_poses=5, options=base)
p1 = sizing.plan([9], 50, 2500000, 16, mmgbsa=True, mm_poses=5,
                 options=dict(base, strain=True))
check(abs(p1["task_minutes"] - p0["task_minutes"] - sizing.strain_minutes(5, 16)) < 0.11
      and p1["strain_minutes"] == 1.0 and p0["strain_minutes"] == 0,
      "the estimate does not add the strain of the scored poses: %s %s"
      % (p0["task_minutes"], p1["task_minutes"]))
check(sizing.strain_minutes(20, 8) == 3.0 and sizing.strain_minutes(0, 8) == 0,
      "the strain's minutes are not waves of a minute a pose")
check(abs(sizing.measures_minutes(35, 16) - (2 * 35 / 16.0 + 30)) < 1e-9,
      "the short job's booking is not 2 x (poses x 1 min / cores) + 30 min")
check("peptide strain" in sizing.suggest([9], 50, 2500000, 16, mmgbsa=True, mm_poses=5,
                                         options=dict(base, strain=True))["why"],
      "the estimate's words do not name the strain")
jn = jobmod.default(mmgbsa=True)
notes = "\n".join(pipeline.batch_notes(jn, "dock"))
check("peptide strain" in notes and "conjugate gradient" in notes,
      "the dock job's notes do not say the strain follows the MM-GBSA:\n%s" % notes)
light = "\n".join(pipeline.batch_notes(jn, "panel+strain"))
check("interface panel" in light and "peptide strain" in light and "docking" not in light,
      "the short job's notes do not say what it does:\n%s" % light)
prime_notes = "\n".join(pipeline.batch_notes(dict(jn, mmgbsa_engine="schrodinger"), "strain"))
check("Prime's own" in prime_notes, "a Prime run's notes say its strain is computed")

# --- what the run says it did --------------------------------------------------
told = cite.as_text(jobmod.default(mmgbsa=True))
check("relaxation of the peptide on its own for its strain" in told,
      "the citation list does not say AmberTools' sander relaxes the peptide")
check("for its strain" not in cite.as_text(dict(jobmod.default(mmgbsa=True),
                                                mmgbsa_strain=False))
      and "for its strain" not in cite.as_text(dict(jobmod.default(mmgbsa=True),
                                                    mmgbsa_engine="schrodinger")),
      "a run that did not relax the peptide with sander is told it did")
from adcp_pilot import runreport  # noqa: E402
fake = type("R", (), {})()
fake.rows = [{"pep_strain": 3.0, "pep_strain_engine": "amber"}]
said = runreport.Report.strain_method(fake, lambda *k: " [4]")
check("5,000 steps of conjugate gradient in AmberTools' sander [4]" in said
      and "not part of the ranking" in said,
      "the report's methods do not describe the relaxation: %s" % said)
fake.rows = [{"pep_strain": 3.0, "pep_strain_engine": "schrodinger"}]
check("Prime's own ligand strain" in runreport.Report.strain_method(fake, lambda *k: ""),
      "a Prime run's methods say the strain was computed")
fake.rows = [{}]
check(runreport.Report.strain_method(fake, lambda *k: "") == "",
      "a run with no strain is given a methods sentence about it")

# --- tidy keeps it ---------------------------------------------------------------
from adcp_pilot import tidy  # noqa: E402
check(tidy.pose_kind("AAAAA_strain.json", "AAAAA", WORK) is None
      and "_strain.json" in tidy.__doc__,
      "the tidy would pack a pose's strain, or does not say it stays")

shutil.rmtree(WORK, ignore_errors=True)
if fails:
    print("test_strain FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the strain is the lowest energy pose's, never ranked on, in every table "
      "file, added to a finished run without touching the rest, %d checks" % ran)
