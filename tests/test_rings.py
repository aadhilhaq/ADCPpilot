"""A ring for each peptide: one run docking some peptides linear and others
closed head to tail or through a disulfide.

By Aadhil Haq

The run's choice (job["cyclic"]) is the default, and a peptide may carry its
own ("cyclic" on the peptide). What is checked here, without docking: the
ring column of a peptide file and a pasted list, with every name people
write for a ring and a plain refusal of anything else; validation peptide by
peptide, each with its own ring, naming the peptide; the ADCP command each
peptide's task builds, read back from the command the task ran; MM-GBSA following the peptide's own ring; the ring in the result
files, the report, the page's results and the estimate; the command line's
file column; and that a job file from before any of this docks and reports
exactly as it did.

    python3 tests/test_rings.py
"""
import contextlib
import io
import json
import os
import shutil
import stat
import sys
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_pilot import (adcp11, cli, job as jobmod, mmforms, peptide,  # noqa: E402
                       pipeline, report, server, sizing)

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)
        print("  FAIL " + why)


def run_cli(*argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            rc = cli.main(list(argv))
        except SystemExit as e:
            rc = e.code
        except Exception as e:                                  # noqa: BLE001
            rc = "CRASH %s: %s" % (type(e).__name__, e)
    return rc, out.getvalue()


TMP = tempfile.mkdtemp(prefix="adcp_rings_")
REC = os.path.join(ROOT, "adcp_pilot", "example_run", "receptor", "1YCR.pdb")
P53, PMI, ONE_C = "ETFSDLWKLLPEN", "TSFAEYWNLLS", "RFMDYWEGLC"
SS = "CTSFAEYWNLLC"

# --- the words for a ring ----------------------------------------------------
for words, want in ((("linear", "none", "Linear", "NONE"), "none"),
                    (("head to tail", "Head-to-tail", "head_to_tail",
                      "backbone", "cyclic", "Cyclic", "-cyc", "backbone cyclic"),
                     "backbone"),
                    (("disulfide", "disulphide", "Disulfide", "-cys"),
                     "disulfide"),
                    (("", "  ", None, "as the run", "default"), "")):
    for w in words:
        check(peptide.ring_choice(w) == want,
              "%r should read as %r, not %r" % (w, want, peptide.ring_choice(w)))
for w in ("circular", "ring", "lariat", "s-s bond", "yes"):
    check(peptide.ring_choice(w) is None, "%r was read as a ring" % w)
msg = peptide.ring_problem("circular").message
check("'circular'" in msg and "linear" in msg and "head to tail" in msg
      and "disulfide" in msg and "empty" in msg,
      "the refusal of a ring does not say what it was given and what the "
      "choices are: %s" % msg)

# --- a peptide file's ring column -------------------------------------------
rows = peptide.parse_rows("name,sequence,ring\n"
                          "p53,%s,head to tail\n"
                          "PMI,%s,\n"
                          "ss_PMI,%s,disulfide\n"
                          "lin,RWGPGAWGSGSGHY,Linear\n" % (P53, PMI, SS))
check([(r["name"], r["sequence"], r["cyclic"]) for r in rows]
      == [("p53", P53, "backbone"), ("PMI", PMI, ""), ("ss_PMI", SS, "disulfide"),
          ("lin", "RWGPGAWGSGSGHY", "none")]
      and not any(r["problems"] for r in rows),
      "a header row with a ring column was read as %s"
      % [(r["name"], r["cyclic"], [p.message for p in r["problems"]]) for r in rows])
rows = peptide.parse_rows("ring;sequence;name\nbackbone;%s;p53\n;%s;PMI\n"
                          % (P53, PMI))
check([(r["name"], r["cyclic"]) for r in rows] == [("p53", "backbone"),
                                                   ("PMI", "")],
      "a header in another order, split on semicolons, was read as %s"
      % [(r["name"], r["cyclic"]) for r in rows])
rows = peptide.parse_rows("sequence\tring\n%s\tcyclic\n" % P53)
check(len(rows) == 1 and rows[0]["cyclic"] == "backbone"
      and rows[0]["sequence"] == P53,
      "a tab separated sequence and ring with a header was read as %s" % rows)
# No header: a third field, in whichever place the ring is.
rows = peptide.parse_rows("p53,%s,head to tail\n%s,PMI,linear\n"
                          "ss_PMI,%s,Disulphide\n" % (P53, PMI, SS))
check([(r["name"], r["sequence"], r["cyclic"]) for r in rows]
      == [("p53", P53, "backbone"), ("PMI", PMI, "none"),
          ("ss_PMI", SS, "disulfide")],
      "three columns with no header were read as %s"
      % [(r["name"], r["sequence"], r["cyclic"]) for r in rows])
# A sequence and its ring, with no name.
rows = peptide.parse_rows("%s,head to tail\n%s,disulfide\n" % (P53, SS))
check([(r["sequence"], r["cyclic"], r["named"]) for r in rows]
      == [(P53, "backbone", False), (SS, "disulfide", False)],
      "a sequence beside its ring was read as %s"
      % [(r["name"], r["sequence"], r["cyclic"]) for r in rows])
# What is wrong is said, never guessed.
rows = peptide.parse_rows("p53,%s,circular\n" % P53)
check(len(rows) == 1 and rows[0]["cyclic"] == ""
      and [p.kind for p in rows[0]["problems"]] == ["ring"]
      and "'circular'" in rows[0]["problems"][0].message,
      "a ring that is no ring was not refused plainly: %s"
      % [(r["cyclic"], [p.message for p in r["problems"]]) for r in rows])
rows = peptide.parse_rows("name,sequence,ring\nx,%s,lariat\n" % P53)
check(rows and any(p.kind == "ring" for p in rows[0]["problems"]),
      "a ring column under a header that says lariat was accepted")
rows = peptide.parse_rows("one_c,%s,disulfide\n" % ONE_C)
check(rows and any("two cysteines" in p.message for p in rows[0]["problems"]),
      "a file row asking for a disulfide with one cysteine was accepted")
# What it always read, it still reads the same.
for text, want in (("TSFAEYWNLLS,PMI", [("PMI", PMI)]),
                   ("PMI,TSFAEYWNLLS", [("PMI", PMI)]),
                   ("sequence,name\n%s,p53" % P53, [("p53", P53)]),
                   ("p53: %s" % P53, [("p53", P53)]),
                   (">p53\n%s\n" % P53, [("p53", P53)]),
                   (P53, [(P53[:12], P53)])):
    got = [(n, s) for n, s, _p in peptide.parse_many(text)]
    check(got == want, "%r used to read as %s and now reads as %s"
          % (text, want, got))
check(all(len(t) == 3 for t in peptide.parse_many("p53,%s,linear" % P53)),
      "parse_many no longer gives (name, sequence, problems)")
tied = peptide.parse_rows("NATIVE,%s" % P53)
check(tied and any(p.kind == "ambiguous_columns" for p in tied[0]["problems"]),
      "two columns that both read as sequences are no longer said to be tied")

# --- the job: each peptide's ring -------------------------------------------
j = jobmod.default(receptor=REC, hotspots="A:54", cyclic="none",
                   peptides=[{"name": "p53", "sequence": P53,
                              "cyclic": "backbone"},
                             {"name": "PMI", "sequence": PMI},
                             {"name": "ss_PMI", "sequence": SS,
                              "cyclic": "disulfide"}])
check(jobmod.rings(j) == ["backbone", "none", "disulfide"],
      "the peptides' rings are %s" % jobmod.rings(j))
check(jobmod.ring_for(j, PMI) == "none" and jobmod.ring_for(j, SS) == "disulfide"
      and jobmod.ring_for(j, "AAAAAAA") == "none",
      "ring_for does not find a peptide's ring by its sequence")
check(jobmod.ring_of(dict(j, cyclic="backbone"), {"sequence": PMI}) == "backbone"
      and jobmod.ring_of(j, {"sequence": PMI, "cyclic": "Head-to-tail"})
      == "backbone"
      and jobmod.ring_of(j, {"sequence": PMI, "cyclic": "as the run"}) == "none"
      and jobmod.ring_of(j, {"sequence": PMI, "cyclic": "lariat"}) == "lariat",
      "ring_of does not give a peptide its own ring, else the run's")
check(jobmod.tidy_peptides([{"sequence": "A", "cyclic": ""},
                            {"sequence": "B", "cyclic": "Head to tail"},
                            {"sequence": "C", "cyclic": "as the run"},
                            {"sequence": "D", "cyclic": "lariat"},
                            {"sequence": "E"}])
      == [{"sequence": "A"}, {"sequence": "B", "cyclic": "backbone"},
          {"sequence": "C"}, {"sequence": "D", "cyclic": "lariat"},
          {"sequence": "E"}],
      "tidy_peptides does not keep a ring the one way, drop an empty one and "
      "keep a wrong one for validation to name")
probs = jobmod.validate(j)
check(not [p for p in probs if "ring" in p.lower() or "cyste" in p.lower()],
      "a run with a linear, a head to tail and a disulfide peptide was "
      "refused: %s" % probs)
# A disulfide refused for the one peptide that asks for it, by name.
bad = dict(j, peptides=j["peptides"] + [{"name": "12/1", "sequence": ONE_C,
                                         "cyclic": "disulfide"}])
probs = [p for p in jobmod.validate(bad) if "cysteine" in p]
check(len(probs) == 1 and probs[0].startswith("Peptide 4 (12/1)")
      and "has one" in probs[0],
      "a disulfide with one cysteine was not refused for that peptide alone, "
      "by name: %s" % probs)
# The run's disulfide does not reach a row that is its own linear chain...
own = jobmod.default(receptor=REC, hotspots="A:54", cyclic="disulfide",
                     peptides=[{"name": "p53", "sequence": P53, "cyclic": "none"},
                               {"name": "ss_PMI", "sequence": SS}])
check(not [p for p in jobmod.validate(own) if "cysteine" in p],
      "a linear row in a disulfide run was refused for its cysteines: %s"
      % jobmod.validate(own))
# ...and does reach a row that takes the run's ring.
own["peptides"][0].pop("cyclic")
probs = [p for p in jobmod.validate(own) if "cysteine" in p]
check(len(probs) == 1 and probs[0].startswith("Peptide 1 (p53)"),
      "a row taking the run's disulfide was not checked for cysteines: %s"
      % probs)
probs = jobmod.validate(dict(j, peptides=[{"name": "x", "sequence": P53,
                                           "cyclic": "lariat"}]))
check(any(p.startswith("Peptide 1 (x)") and "'lariat'" in p for p in probs),
      "a ring that is no ring was not refused by the peptide's name: %s" % probs)
probs = jobmod.validate(dict(j, peptides=[
    {"name": "lin", "sequence": P53, "cyclic": "none"},
    {"name": "ring", "sequence": P53, "cyclic": "backbone"}]))
check(any("same sequence" in p and "as a linear chain" in p
          and "as a ring head to tail" in p and "run of its own" in p
          for p in probs),
      "one sequence asked as a chain and as a ring is not explained: %s" % probs)

# --- the ADCP 1.1 command each peptide's task builds -------------------------
# A stand in for the suite: its "python" writes down the arguments it was
# given, which is what ADCP would have run on.
home = os.path.join(TMP, "adcpsuite")
pkg = os.path.join(home, "lib", "python3.9", "site-packages", "ADCP")
os.makedirs(pkg)
os.makedirs(os.path.join(home, "bin"))
open(os.path.join(pkg, "runADCP.py"), "w").write('version="1.1.%d"\n')
open(os.path.join(pkg, "utils.py"), "w").write("BUILD_VERSION = 21\n")
for name in ("adcp", "python"):
    f = os.path.join(home, "bin", name)
    open(f, "w").write('#!/bin/sh\nprintf "%s\\n" "$@" > argv.txt\n')
    os.chmod(f, os.stat(f).st_mode | stat.S_IEXEC)
EXE11 = os.path.join(home, "bin", "adcp")


def make_run(name, j):
    rd = os.path.join(TMP, name)
    os.makedirs(os.path.join(rd, "target"))
    trg = os.path.join(rd, "target", "receptor.trg")
    open(trg, "w").write("a target")
    jobmod.save(rd, dict(j, target=trg))
    return rd


_check, _engine = pipeline.target.check, jobmod.resolve_engine
pipeline.target.check = lambda trg, *a, **k: None
try:
    jobmod.resolve_engine = lambda jj: ("1.1", EXE11)
    rd = make_run("dock11", dict(j, replicas=8, steps=1000, poses=5, cpus=2,
                                 helical=25))
    logs = {}
    for i in range(3):
        lines = []
        pipeline.stage_dock(rd, i, log=lines.append)
        logs[i] = lines
    argv = {}
    for i, seq in enumerate((P53, PMI, SS)):
        f = os.path.join(rd, "docking", seq, "argv.txt")
        argv[i] = open(f).read().split() if os.path.isfile(f) else []
        run_log = os.path.join(rd, "docking", seq, "adcp_run.log")
        first = open(run_log).readline() if os.path.isfile(run_log) else ""
        check(all(a in first for a in argv[i][1:]) and argv[i],
              "adcp_run.log does not begin with the command %s ran" % seq)
    check("-cyc" in argv[0] and "-cys" not in argv[0],
          "p53, head to tail, was not given -cyc alone: %s" % argv[0])
    check("-cyc" not in argv[1] and "-cys" not in argv[1],
          "PMI, linear, was given a ring: %s" % argv[1])
    check("-cys" in argv[2] and "-cyc" not in argv[2],
          "ss_PMI, a disulfide, was not given -cys alone: %s" % argv[2])
    check(all(argv[i][argv[i].index("-p") + 1] == "25" for i in range(3)),
          "the run's other settings did not reach every peptide: %s" % argv)
    for i, flag in ((0, "-cyc"), (2, "-cys")):
        said = [l for l in logs[i] if l.strip().startswith("adcp:")]
        check(len(said) == 1 and flag in said[0],
              "peptide %d's task log does not carry its ADCP command with %s: "
              "%s" % (i + 1, flag, said))
    check("backbone cyclic, head to tail" in logs[0][0]
          and "; linear;" in logs[1][0]
          and "cyclic through a disulfide" in logs[2][0],
          "each task's first line does not name its own peptide's ring: %s"
          % [logs[i][0] for i in range(3)])
    got = {s: json.load(open(os.path.join(rd, "results", s + ".json")))
           .get("cyclic") for s in (P53, PMI, SS)}
    check(got == {P53: "backbone", PMI: "none", SS: "disulfide"},
          "the result files do not record how each peptide was docked: %s" % got)
    notes = "\n".join(pipeline.batch_notes(jobmod.load(rd), "dock"))
    check("rings by peptide: PMI linear; p53 head to tail; ss_PMI disulfide"
          in notes,
          "the dock script's notes do not say which peptide has which ring:\n%s"
          % notes)

    # A task given a ring no one can close stops, and says so.
    rdx = make_run("dockx", dict(j, peptides=[{"name": "x", "sequence": P53,
                                               "cyclic": "lariat"}]))
    try:
        pipeline.stage_dock(rdx, 0, log=lambda *a: None)
        stopped = ""
    except SystemExit as exc:
        stopped = str(exc)
    check("lariat" in stopped and not os.path.isdir(
        os.path.join(rdx, "docking", P53)),
        "a task with a ring that is no ring docked anyway (%r)" % stopped)
finally:
    pipeline.target.check, jobmod.resolve_engine = _check, _engine

# --- the report: each peptide's ring ----------------------------------------
pipeline.stage_report(rd, jobmod.load(rd), log=lambda *a: None)
tsv = open(os.path.join(rd, "results.tsv")).read().splitlines()
head = tsv[0].split("\t")
rings = {}
if "ring" in head:
    for line in tsv[1:]:
        f = line.split("\t")
        rings[f[head.index("peptide")]] = f[head.index("ring")]
check(rings == {P53: "head to tail", PMI: "linear", SS: "disulfide"},
      "results.tsv does not carry each peptide's ring: %s" % rings)
txt = open(os.path.join(rd, "results.txt")).read()
check(" ring " in txt.splitlines()[2] and "head to tail" in txt
      and "disulfide" in txt,
      "results.txt does not show each peptide's ring:\n%s" % txt[:600])
res = json.load(open(os.path.join(rd, "results.json")))
check({r["peptide"]: r.get("ring") for r in res["rows"]}
      == {P53: "head to tail", PMI: "linear", SS: "disulfide"},
      "results.json rows do not carry the ring")
# A run of linear peptides has no ring column, as before.
lin_rows = [report.summarise(P53, [{"score": -10.0}], ring="none"),
            report.summarise(PMI, [{"score": -9.0}], ring="none")]
out = os.path.join(TMP, "linear_report")
report.write(lin_rows, out)
check("ring" not in open(os.path.join(out, "results.tsv")).readline().split("\t")
      and " ring " not in open(os.path.join(out, "results.txt")).read(),
      "a run of linear peptides gained a ring column")

# --- the page's results, for a run reported before rows carried a ring -----
old = {"rows": [{"peptide": P53, "name": "p53"}, {"peptide": PMI}]}
server.fill_rings(old, {"cyclic": "backbone",
                        "peptides": [{"name": "p53", "sequence": P53},
                                     {"name": "PMI", "sequence": PMI,
                                      "cyclic": "none"}]})
check([r.get("ring") for r in old["rows"]] == ["head to tail", "linear"],
      "fill_rings does not give an older report's rows their rings: %s"
      % old["rows"])
kept = {"rows": [{"peptide": P53, "ring": "disulfide"}]}
server.fill_rings(kept, {"cyclic": "backbone"})
check(kept["rows"][0]["ring"] == "disulfide",
      "fill_rings overwrote the ring a row already recorded")

# --- MM-GBSA follows the peptide's own ring ---------------------------------


def write_open_pose(pdir, k, n=6, gap=4.3):
    """A pose whose first N and last C are `gap` A apart: an open ring."""
    lines, a = [], 0
    for r in range(1, n + 1):
        for nm, dx in (("N", 0.0), ("CA", 1.0), ("C", 2.0), ("O", 3.0)):
            a += 1
            x, y = 3.8 * r, dx
            if r == n and nm == "C":
                x, y = 3.8 + gap, 0.0
            lines.append("ATOM  %5d  %-3s ALA A%4d    %8.3f%8.3f%8.3f  1.00  "
                         "0.00           %s\n" % (a, nm, r, x, y, 0.0, nm[0]))
    open(os.path.join(pdir, "out_ranked_%d.pdb" % k), "w").writelines(
        lines + ["END\n"])


ring_seq, lin_seq = "AAAAAA", "AAAAAG"
mmrun = os.path.join(TMP, "mm")
for seq in (ring_seq, lin_seq):
    d = os.path.join(mmrun, "docking", seq)
    os.makedirs(d)
    for k in (1, 2, 3):
        write_open_pose(d, k)
mj = jobmod.default(cyclic="none", mmgbsa=True, mmgbsa_poses=3,
                    peptides=[{"name": "ring", "sequence": ring_seq,
                               "cyclic": "backbone"},
                              {"name": "lin", "sequence": lin_seq}])
_m, _r, skipped = pipeline._mm_inputs(mj, ring_seq,
                                      os.path.join(mmrun, "docking", ring_seq),
                                      3, log=lambda *a: None)
check(len(skipped) == 3 and all(s["why"].startswith("ring not closed")
                                for s in skipped),
      "the head to tail peptide's open poses were not held to its ring: %s"
      % skipped)
_m, _r, skipped = pipeline._mm_inputs(mj, lin_seq,
                                      os.path.join(mmrun, "docking", lin_seq),
                                      3, log=lambda *a: None)
check(not any(s["why"].startswith("ring") for s in skipped),
      "the linear peptide in the same run was held to a ring: %s" % skipped)
out = {"modes": [{"rank": 1}, {"rank": 2}, {"rank": 3}]}
pipeline._stage_mmgbsa(mmrun, mj, ring_seq, out,
                       os.path.join(mmrun, "docking", ring_seq), 3, 1,
                       dict(pipeline.docking_options(mj), cyclic="backbone"),
                       log=lambda *a: None)
check(out.get("mmgbsa_skipped") == "ring not closed in any of its top 3 poses"
      and not os.path.isdir(os.path.join(mmrun, "mmgbsa")),
      "a ring of the peptide's own, open in every pose, was not refused with "
      "the reason, or was built anyway: %s" % out.get("mmgbsa_skipped"))
# The run's ring does not reach a peptide that is linear by its own choice.
mj2 = dict(mj, cyclic="backbone",
           peptides=[{"name": "lin", "sequence": lin_seq, "cyclic": "none"}])
_m, _r, skipped = pipeline._mm_inputs(mj2, lin_seq,
                                      os.path.join(mmrun, "docking", lin_seq),
                                      3, log=lambda *a: None)
check(not any(s["why"].startswith("ring") for s in skipped),
      "a linear peptide in a head to tail run was held to the run's ring: %s"
      % skipped)
check(pipeline.docking_options(mj2, seq=lin_seq)["cyclic"] == "none"
      and pipeline.docking_options(mj2)["cyclic"] == "backbone",
      "docking_options does not give a peptide its own ring")
check(mmforms.engine_refusal("schrodinger", ring_seq,
                             jobmod.ring_for(mj, ring_seq))
      and mmforms.engine_refusal("schrodinger", lin_seq,
                                 jobmod.ring_for(mj, lin_seq)) is None,
      "Prime does not refuse the ring alone in a run of both")
check(pipeline.result_ring(mj, {"peptide": lin_seq, "cyclic": "disulfide"})
      == "disulfide"
      and pipeline.result_ring(mj, {"peptide": ring_seq}) == "backbone",
      "a rescoring does not read the ring the result recorded, else the job's")
old_mm = dict(mj)
old_mm.pop(mmforms.JOB_KEY)
check(mmforms.unscoreable(ring_seq, old_mm) == mmforms.SKIP_OLD_RING
      and mmforms.unscoreable(lin_seq, old_mm) is None,
      "a run from before rings were rescored judged a peptide by another's "
      "ring")

# --- the estimate prices each peptide with its own ring --------------------
base = dict(agfr="1.1", autosite="1.0", receptor_gradient=False,
            site_mode="residues", cyclic="none", helical=0, poses=10,
            keep_trajectories=False)
long_lin = sizing.plan([16, 10], 50, 2500000, 48, options=dict(base))
short_ring = sizing.plan([16, 10], 50, 2500000, 48,
                         options=dict(base, rings=["none", "backbone"]))
check(short_ring["task_minutes"] > long_lin["task_minutes"]
      and short_ring["sized_for"] == {"length": 10, "ring": "backbone"}
      and short_ring["rings_differ"],
      "a short peptide closed head to tail was not priced as the slower "
      "task: %s" % short_ring.get("sized_for"))
check(short_ring["core_hours"] > long_lin["core_hours"],
      "the ring was not counted in the core hours")
same = sizing.plan([16, 10], 50, 2500000, 48,
                   options=dict(base, cyclic="backbone"))
again = sizing.plan([16, 10], 50, 2500000, 48,
                    options=dict(base, cyclic="backbone",
                                 rings=["backbone", ""]))
check(same["task_minutes"] == again["task_minutes"] and not same["rings_differ"]
      and same["sized_for"]["length"] == 16,
      "a run whose peptides all take the run's ring is priced differently")
why = sizing.suggest([16, 10], 50, 2500000, 48,
                     options=dict(base, rings=["none", "backbone"]))["why"]
check("takes longest to dock (10 residues, head to tail)" in why,
      "the estimate does not say it was sized for the ring: %s" % why)
check("for the longest peptide (16 residues)" in sizing.suggest(
    [16, 10], 50, 2500000, 48, options=dict(base))["why"],
      "a run of linear peptides no longer says it was sized for the longest")

# --- the page's server: live checks, a paste, the estimate, a new run ------
app = server.App(os.path.join(TMP, "runs"))
got = app.peptides({"cyclic": "disulfide", "peptides": [
    {"sequence": P53, "cyclic": "none"},        # its own linear chain
    {"sequence": SS},                           # the run's disulfide
    {"sequence": ONE_C},                        # the run's, one cysteine
    {"sequence": PMI, "cyclic": "backbone"},
    {"sequence": PMI + "K", "cyclic": "lariat"}]})["peptides"]
check([p["ok"] for p in got] == [True, True, False, True, False]
      and [p["cyclic"] for p in got] == ["none", "disulfide", "disulfide",
                                         "backbone", "lariat"],
      "the page's live check does not judge each row by its own ring: %s"
      % [(p["ok"], p["cyclic"]) for p in got])
check(any(x["kind"] == "cyclic" and "two cysteines" in x["message"]
          for x in got[2]["problems"])
      and any(x["kind"] == "ring" for x in got[4]["problems"]),
      "the rows refused do not say why")
pasted = app.paste({"text": "p53,%s,head to tail\nPMI,%s\nx,%s,lariat"
                    % (P53, PMI, P53)})["peptides"]
check([(p["name"], p["cyclic"], p["ok"]) for p in pasted]
      == [("p53", "backbone", True), ("PMI", "", True), ("x", "", False)],
      "a paste does not carry each row's ring: %s" % pasted)
est = app.estimate_run({"peptides": [{"sequence": "A" * 16},
                                     {"sequence": PMI[:10],
                                      "cyclic": "backbone"}],
                        "replicas": 50, "steps": 2500000, "cores": 48})
check(est.get("sized_for") == {"length": 10, "ring": "backbone"},
      "the page's estimate does not price a row's own ring: %s"
      % est.get("sized_for"))
seen = []
_submit = pipeline.submit
pipeline.submit = lambda rundir, j, **kw: seen.append(j) or {"dock": "1"}
try:
    made = app.create({"label": "rings", "receptor": REC, "hotspots": "A:54",
                       "cyclic": "none",
                       "peptides": [{"name": "p53", "sequence": P53,
                                     "cyclic": "backbone"},
                                    {"name": "PMI", "sequence": PMI,
                                     "cyclic": ""},
                                    {"name": "ss_PMI", "sequence": SS,
                                     "cyclic": "disulfide"}]})
finally:
    pipeline.submit = _submit
saved = (jobmod.load(made["rundir"])["peptides"]
         if made.get("rundir") and os.path.isfile(
             jobmod.path_of(made["rundir"])) else [])
check(not made.get("error") and saved == [
    {"name": "p53", "sequence": P53, "cyclic": "backbone"},
    {"name": "PMI", "sequence": PMI},
    {"name": "ss_PMI", "sequence": SS, "cyclic": "disulfide"}],
      "a run made from the page does not keep each row's ring, or keeps an "
      "empty one: %s %s" % (made.get("error"), saved))
made = app.create({"label": "rings_bad", "receptor": REC, "hotspots": "A:54",
                   "cyclic": "none",
                   "peptides": [{"name": "12/1", "sequence": ONE_C,
                                 "cyclic": "disulfide"}]})
check(made.get("error") and any("12/1" in p and "two cysteines" in p
                                for p in made.get("problems") or []),
      "the page's server accepted a disulfide with one cysteine: %s" % made)

# --- the command line: the file's ring column, with --cyclic the default ----
pf = os.path.join(TMP, "peptides.csv")
open(pf, "w").write("name,sequence,ring\np53,%s,head to tail\nPMI,%s,\n"
                    "ss_PMI,%s,disulfide\n" % (P53, PMI, SS))
rd = os.path.join(TMP, "cli_rings")
rc, said = run_cli("new", rd, "-r", REC, "--chain", "A", "-s", "A:54",
                   "-f", pf)
jf = os.path.join(rd, "job.json")
jc = json.load(open(jf)) if os.path.isfile(jf) else {}
check(jobmod.rings(jc) == ["backbone", "none", "disulfide"]
      and jc.get("cyclic") == "none"
      and "cyclic" not in jc["peptides"][1],
      "adcp-pilot new -f did not keep the file's rings with --cyclic the "
      "default (rc %r): %s" % (rc, jc.get("peptides")))
rd = os.path.join(TMP, "cli_rings_ss")
rc, said = run_cli("new", rd, "-r", REC, "--chain", "A", "-s", "A:54",
                   "-f", pf, "--cyclic", "disulfide")
said = (rc if isinstance(rc, str) else "") + said
check(rc not in (0, None) and "PMI" in said and "two cysteines" in said
      and "p53" not in said.split("cannot be docked", 1)[-1]
      and not os.path.isfile(os.path.join(rd, "job.json")),
      "--cyclic disulfide was not refused for the one row it reaches, by "
      "name (rc %r): %s" % (rc, said[-400:]))
bf = os.path.join(TMP, "bad.csv")
open(bf, "w").write("p53,%s,lariat\n" % P53)
rd = os.path.join(TMP, "cli_rings_bad")
rc, said = run_cli("new", rd, "-r", REC, "--chain", "A", "-s", "A:54",
                   "-f", bf)
said = (rc if isinstance(rc, str) else "") + said
check(rc not in (0, None) and "'lariat'" in said
      and not os.path.isfile(os.path.join(rd, "job.json")),
      "a file whose ring column says lariat was accepted (rc %r): %s"
      % (rc, said[-300:]))

# --- a job file from before peptides had rings of their own ----------------
ex = os.path.join(ROOT, "adcp_pilot", "example_run")
# The shipped example as a job file from before rings: without the run's
# ring, and with no peptide carrying one.
old_job = json.load(open(os.path.join(ex, "job.json")))
old_job.pop("cyclic", None)
check(not any("cyclic" in p for p in old_job["peptides"]),
      "the shipped example's peptides carry rings of their own")
check(jobmod.rings(old_job) == ["none"] * len(old_job["peptides"])
      and pipeline.describe_rings(old_job) == "linear"
      and pipeline.docking_options(old_job)["cyclic"] == "none",
      "an old job file no longer docks every peptide linear")
old_ring = dict(old_job, cyclic="backbone")
check(jobmod.rings(old_ring) == ["backbone"] * len(old_job["peptides"])
      and pipeline.describe_rings(old_ring) == "backbone cyclic, head to tail",
      "an old run with one ring for all no longer docks them all with it")
copy = os.path.join(TMP, "example")
shutil.copytree(ex, copy)
pipeline.stage_report(copy, json.load(open(os.path.join(copy, "job.json"))),
                      log=lambda *a: None)
check("ring" not in open(os.path.join(copy, "results.tsv")).readline()
      .split("\t"),
      "the shipped example, all linear, gained a ring column")
check(all(r.get("ring") == "linear" for r in
          json.load(open(os.path.join(copy, "results.json")))["rows"]),
      "the shipped example's rows do not say they were docked linear")
# An old ring run, reported again: its one ring is every peptide's. Its
# result files, written before a result recorded its own ring, say none.
j_old = json.load(open(os.path.join(copy, "job.json")))
j_old["cyclic"] = "backbone"
jobmod.save(copy, j_old)
for _f in os.listdir(os.path.join(copy, "results")):
    _p = os.path.join(copy, "results", _f)
    _r = json.load(open(_p))
    _r.pop("cyclic", None)
    json.dump(_r, open(_p, "w"), indent=1)
pipeline.stage_report(copy, j_old, log=lambda *a: None)
rows = json.load(open(os.path.join(copy, "results.json")))["rows"]
check([r.get("ring") for r in rows] == ["head to tail"] * len(rows),
      "an old run with one ring for every peptide does not report it on "
      "every row: %s" % [r.get("ring") for r in rows])

# --- the page's script -------------------------------------------------------
appjs = open(os.path.join(ROOT, "adcp_pilot", "static", "app.js")).read()
for want, why in (("class=pr", "a row has no ring of its own to choose"),
                  ("['', 'as the run']", "a row's ring does not default to the run's"),
                  ("cyclic: p.cyclic || ''", "a reused run, a paste or an "
                   "upload drops the rows' rings"),
                  ("The rows differ", "the note does not say what the run "
                   "will do when the rows differ"),
                  ("x.kind === 'ring'", "an upload with a ring that is no "
                   "ring is not refused")):
    check(want in appjs, why)

shutil.rmtree(TMP, ignore_errors=True)
if fails:
    print("test_rings FAILED (%d of %d)" % (len(fails), ran))
    sys.exit(1)
print("ok  each peptide its own ring: files, checks, commands, MM-GBSA, "
      "reports, estimates and old jobs, %d checks" % ran)
