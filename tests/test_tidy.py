"""Tidy: what a finished run packs, what it keeps, and that Undo puts it all
back.

By Aadhil Haq

The fixture is the bundled example run (three peptides docked and rescored),
given what real runs leave beside their results: the archive of an ADCP 1.0
docking converted to today's layout for one peptide, ADCP 1.1's working folder
and packed trajectories for another, the copy of the grid target in each, MM-GBSA's
working files in every pose folder, Prime's output in one, a pose folder with
no prepared complex, one peptide whose panel was never measured, and a peptide
with no result yet.

What is checked: the plan names exactly the intermediates and none of what a
reader opens; a dry run changes nothing; a run with a live job, and a folder
that is not a run, are refused; after the tidy the results page, the pose
picker, the viewer's structures, the residue map, the energy terms, the
re-ranking, the rebuilt report and the downloads give what they gave before;
a second tidy finds nothing; Undo restores every file byte for byte, with its
permissions; the worker's tidy and untidy stages and the local executor work
without Slurm; a Slurm submission writes its script and records its job apart
from job.json; the page's route and payload are wired; and the button shows
each of its states, under Node where it is found.

    python3 tests/test_tidy.py

Exit 0 passes, 77 skips (the package cannot be imported here), anything else
fails.
"""
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import time
import urllib.request
from http.server import ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

try:
    from adcp_pilot import (bundle, example, pipeline, server, slurm,  # noqa: E402
                           tidy, worker)
    from adcp_pilot._vendor import adcp as adcpvendor, mmgbsa     # noqa: E402
except ImportError as exc:
    print("skip: the package cannot be imported here (%s)" % exc)
    sys.exit(77)

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


def write(path, text="x\n", mode=0o664):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, mode)


def snapshot(root):
    """{path in the run: (sha1, permissions)} for every file, hidden ones too."""
    out = {}
    for d, _dirs, files in os.walk(root):
        for f in files:
            p = os.path.join(d, f)
            with open(p, "rb") as fh:
                out[os.path.relpath(p, root)] = (
                    hashlib.sha1(fh.read()).hexdigest(),
                    stat.S_IMODE(os.lstat(p).st_mode))
    return out


MMPBSA_IN = ("Single-structure MM-GBSA rescoring of a docked pose\n&general\n"
             "  startframe=1, endframe=1, interval=1, verbose=2, keep_files=0,\n"
             "  receptor_mask=\":1-85\", ligand_mask=\":86-96\",\n/\n&gb\n"
             "  igb=5, saltcon=0.150,\n/\n")
#: What a pose folder keeps and what it packs, by what follows its name.
KEEP_TAILS = ("_mmpbsa.in", "_decomp.dat", "_min.out", "_com.parm7",
              "_rec.parm7", "_lig.parm7", "_min.rst7", "_com.rst7",
              "_scored.pdb", "_scored_on.txt", "_com_forleap.pdb",
              "_com_forleap.json", "_com_prepped.pdb", "_chains.json",
              "_checks.json", "_com.log", "_mirror_leap.log", "_strain.json")
PACK_TAILS = (".leap", "_rec.pdb", "_lig.pdb", "_com.pdb", "_com.pqr",
              "_com_reduce.pdb", "_min.in", "_decomp.in", "_decomp_total.dat",
              "_mirror.leap", "_mirror.parm7", "_mirror.rst7",
              "_mirror_pep.pdb", "_rec_forleap.pdb", "_pep_forleap.pdb",
              "_com_tleap.rst7")
PACK_NAMES = ("mdinfo", "reference.frc", "_MMPBSA_complex_gb.mdout")

A, B, C = "ETFSDLWKLLPEN", "TSFAEYWNLLS", "RFMDYWEGL"
work = tempfile.mkdtemp(prefix="adcp_tidy_")
runs = os.path.join(work, "runs")
run = example.install(runs)
J = lambda *p: os.path.join(run, *p)                          # noqa: E731
expect = set()

# ---- the fixture ------------------------------------------------------------
# A: an ADCP 1.0 docking converted to today's layout: its replicas packed
# into an archive of their own, which a tidy leaves alone, and a replica file
# of the old kind left loose, which nothing reads and a tidy no longer knows.
write(J("docking", A, "adcp10_replicas.tar.gz"), "replicas\n")
write(J("docking", A, "out_7.pdb"), "replica 7\n")
write(J("docking", A, "receptor.trg"), "receptor.trg\n")
expect.add(os.path.join("docking", A, "receptor.trg"))
write(J("docking", A, ".adcp_finished"), "25 replicas\n")
# B: ADCP 1.1, its trajectories packed by the docking, a working folder left.
for n in ("modes.json", "out_out.pdb", "out_summary.dlg", "adcp_run.log",
          "openmm.json", ".adcp_finished", ".trajectories_packed",
          "%s_trajectories.tar.gz" % B):
    write(J("docking", B, n), "[]\n" if n.endswith(".json") else n + "\n")
for n in ("receptor.trg", os.path.join("out", "out_1.pdb"),
          os.path.join("out", "out_1.log"), os.path.join("out", "sub", "x.dat")):
    write(J("docking", B, n), n + "\n")
    expect.add(os.path.join("docking", B, n))
write(J("docking", "notes.txt"), "a file of the person's own\n")
# Every pose folder: what MM-GBSA and the panel leave.
for pep in (A, B, C):
    for k in range(1, 6):
        pd = J("mmgbsa", pep, "pose%d" % k)
        for t in KEEP_TAILS:
            if pep == C and k == 5 and t.startswith("_com_forleap"):
                continue                # no prepared complex: _com.pdb stays
            write(os.path.join(pd, pep + t),
                  MMPBSA_IN if t == "_mmpbsa.in" else
                  "ATOM      1  N   ALA A   1\n" if t.endswith(".pdb") else
                  "{}\n" if t.endswith(".json") else t + "\n")
        write(os.path.join(pd, "leap.log"), "tleap\n")
        for t in PACK_TAILS:
            write(os.path.join(pd, pep + t), t + "\n")
            if not (pep == C and k == 5 and t == "_com.pdb"):
                expect.add(os.path.relpath(os.path.join(pd, pep + t), run))
        for n in PACK_NAMES:
            write(os.path.join(pd, n), n + "\n")
            expect.add(os.path.relpath(os.path.join(pd, n), run))
for n in ("FINAL_RESULTS_MMGBSA.csv", B + "-out.maegz", B + "_com.mae"):
    write(J("mmgbsa", B, "pose1", n), "Prime\n")
write(J("mmgbsa", B, "summary.csv"), "pose,dg\n")
write(J("receptor", "tidied.pdb"), "END\n")
write(J("receptor", "protonated.pdb"), "END\n")
write(J("target", "receptor.trg"), "maps\n")
# C was docked before the panel existed: its terms come from the files.
rc = J("results", C + ".json")
r = json.load(open(rc))
r.pop("panel", None)
r.pop("panel_poses", None)
json.dump(r, open(rc, "w"), indent=1)
# A peptide still docking: no result file.
for n in ("out_1.pdb", "out_1.out", "ramaprob.data", "receptor.trg"):
    write(J("docking", "GGGGG", n), n + "\n")
for n in ("GGGGG.leap", "GGGGG_rec.pdb", "mdinfo"):
    write(J("mmgbsa", "GGGGG", "pose1", n), n + "\n")
os.chmod(J("docking", A, "receptor.trg"), 0o640)  # an unusual permission, kept

# ---- the plan -----------------------------------------------------------------
files, left = tidy.plan(run)
got = set(files)
check(got == expect, "the plan is not exactly the intermediates:\n      extra "
      "%s\n      missing %s" % (sorted(got - expect)[:6], sorted(expect - got)[:6]))
check(len(files) == len(got), "the plan names a file twice")
check(sorted(f for f, _w in left) == [os.path.join("docking", "GGGGG"),
                                      os.path.join("mmgbsa", "GGGGG")],
      "a peptide with no result is not left whole, with the reason: %s" % left)
check(not any("GGGGG" in f for f in files), "a peptide with no result is packed")
LOAD_BEARING = ["job.json", "results.json", "results.tsv", "results.txt",
                os.path.join("results", A + ".json"),
                os.path.join("metrics", A + ".json"),
                os.path.join("receptor", "prepared.pdb"),
                os.path.join("receptor", "prepared.json"),
                os.path.join("receptor", "preparation.log"),
                os.path.join("receptor", "tidied.pdb"),
                os.path.join("target", "receptor.trg"),
                os.path.join("mmgbsa", "RADII"), os.path.join("mmgbsa", A + ".dg"),
                os.path.join("mmgbsa", B, "summary.csv"),
                os.path.join("docking", A, "adcp.log"),
                os.path.join("docking", A, ".adcp_finished"),
                os.path.join("docking", A, "adcp10_replicas.tar.gz"),
                os.path.join("docking", A, "out_7.pdb"),
                os.path.join("docking", "notes.txt"),
                os.path.join("mmgbsa", C, "pose5", C + "_com.pdb"),
                os.path.join("mmgbsa", A, "pose1", "leap.log")]
LOAD_BEARING += [os.path.join("docking", B, n) for n in
                 ("modes.json", "out_out.pdb", "out_summary.dlg", "adcp_run.log",
                  "openmm.json", ".adcp_finished", ".trajectories_packed",
                  "%s_trajectories.tar.gz" % B)]
LOAD_BEARING += [os.path.join("docking", p, "out_ranked_%d.pdb" % k)
                 for p in (A, B, C) for k in range(1, 11)]
LOAD_BEARING += [os.path.join("mmgbsa", p, "pose%d" % k, p + t)
                 for p in (A, B, C) for k in range(1, 5)
                 for t in KEEP_TAILS + ("_mmgbsa.dat",)]
LOAD_BEARING += [os.path.join("mmgbsa", B, "pose1", n) for n in
                 ("FINAL_RESULTS_MMGBSA.csv", B + "-out.maegz", B + "_com.mae")]
missing = [f for f in LOAD_BEARING if not os.path.isfile(J(f))]
check(not missing, "the fixture lacks %s" % missing[:4])
check(not got & set(LOAD_BEARING),
      "the plan packs what a reader opens: %s" % sorted(got & set(LOAD_BEARING))[:6])
kinds = tidy.breakdown(run)
check(kinds.get("ADCP 1.1's working folder") == 3
      and kinds.get("the copy of the grid target ADCP docked from") == 2,
      "the plan's breakdown is not what the fixture holds: %s" % kinds)

# ---- what the readers give before -------------------------------------------
app = server.App(runs)


def readings(rundir):
    """Everything the page, the viewer, the report and the downloads read."""
    st = app.run_status(rundir)
    out = {"job": st["job"], "done": st["done"], "outstanding": st["outstanding"],
           "prep": st["prep_warnings"], "results": st["results"],
           "log": st["log_tail"]}
    for p in (A, B, C):
        out["poses " + p] = app.pose_list(rundir, p)
        out["residues " + p] = app.residues(rundir, p)
        out["complex " + p] = app.structure_text(rundir, p, 1, "complex")
        out["pose " + p] = app.structure_text(rundir, p, 3, "pose")
        res = json.load(open(os.path.join(rundir, "results", p + ".json")))
        out["terms " + p] = pipeline.energy_terms(rundir, res)
        out["clusters " + p] = pipeline.cluster_counts(rundir, res)
        out["measurable " + p] = pipeline.measurable(rundir, res)
        out["engines " + p] = mmgbsa.engines_on_disk(
            os.path.join(rundir, "mmgbsa"), p)
        out["read_pose " + p] = [None if q is None else q["xyz"].tolist()
                                 for q in (adcpvendor.read_pose(
                                     os.path.join(rundir, "docking", p), k)
                                     for k in (1, 5))]
    out["receptor"] = app.structure_text(rundir, None, 1, "receptor")
    out["rerank"] = app.rerank_run(rundir)
    out["radii"] = mmgbsa.radii_for(os.path.join(rundir, "mmgbsa"))
    out["logs"] = sorted(r_ for r_, _p in bundle.walk(rundir, bundle.LOGS))
    return out


before_reads = readings(run)
results_dl = set(r_ for r_, _p in bundle.walk(run, bundle.RESULTS))
check(before_reads["results"] and before_reads["results"].get("rows"),
      "the fixture's results page has no rows")
check(before_reads["terms " + C] and before_reads["measurable " + C],
      "the fixture does not exercise the energy terms and the panel offer: %s"
      % before_reads["terms " + C])
before = snapshot(run)

# ---- refusals -------------------------------------------------------------------
r = tidy.tidy(run, dry_run=True)
check(r.get("dry_run") and r["packed"] == len(files) and snapshot(run) == before,
      "a dry run changes nothing and says what it would pack: %s"
      % {k: v for k, v in r.items() if k != "files"})
check(r.get("inodes_after") == r.get("inodes_before") - len(files),
      "a dry run's count after is not the count before less what it packs")
nonrun = tempfile.mkdtemp(dir=work)
write(os.path.join(nonrun, "out_1.pdb"))
for what, got_ in (("tidy", tidy.tidy(nonrun)), ("undo", tidy.undo(nonrun)),
                   ("submit", tidy.submit(nonrun))):
    check("error" in got_ and "not a run" in got_["error"],
          "a folder that is not a run is not refused by %s: %s" % (what, got_))
check(os.path.isfile(os.path.join(nonrun, "out_1.pdb")),
      "a folder that is not a run was touched")
real_live = tidy.live_jobs
tidy.live_jobs = lambda d: ["dock 4242"]
r = tidy.tidy(run)
s = tidy.submit(run)
u = tidy.submit(run, undo_=True)
tidy.live_jobs = real_live
check("error" in r and "4242" in r["error"] and snapshot(run) == before,
      "a run with a job still live is left alone: %s" % r)
check("error" in s and "4242" in s["error"] and "error" in u,
      "the page's submission refuses a run with a live job: %s %s" % (s, u))
check(not os.path.exists(J(tidy.JOBFILE)), "a refused submission recorded a job")
check(tidy.live_jobs(run) == [],
      "the example's own local jobs, long finished, read as live")
u = tidy.submit(run, undo_=True)
check("error" in u and "not been tidied" in u["error"],
      "undo is offered on a run that was never tidied: %s" % u)

# ---- a report rebuilt from a tidied copy matches one from the untidied run ----
twin = os.path.join(work, "twin")
shutil.copytree(run, twin, symlinks=True)
quiet = lambda *a: None                                         # noqa: E731
saved_env = {k: os.environ.pop(k, None) for k in ("SLURM_JOB_ID", "ADCP_EXECUTOR")}
try:
    rt = tidy.tidy(twin)
    check(not rt.get("error") and rt.get("packed") == len(files),
          "the copy could not be tidied: %s" % rt)
    pipeline.stage_report(run, log=quiet)
    pipeline.stage_report(twin, log=quiet)
    a_tsv = open(J("results.tsv")).read()
    b_tsv = open(os.path.join(twin, "results.tsv")).read()
    check(a_tsv == b_tsv and a_tsv.count("\n") >= 4,
          "the report rebuilt from a tidied run differs from the untidied one's")
    check(json.load(open(J("results.json")))["rows"]
          == json.load(open(os.path.join(twin, "results.json")))["rows"],
          "the report's rows differ once the run is tidied")
    pv = sorted(os.listdir(J("poseviewer"))) if os.path.isdir(J("poseviewer")) else []
    check(pv == (sorted(os.listdir(os.path.join(twin, "poseviewer")))
                 if os.path.isdir(os.path.join(twin, "poseviewer")) else []),
          "the pose viewer files differ once the run is tidied")
finally:
    for k, v in saved_env.items():
        if v is not None:
            os.environ[k] = v
shutil.rmtree(twin)
# The report rewrote the run's tables and viewer files: they are the baseline
# from here on.
before_reads = readings(run)
before = snapshot(run)

# ---- the tidy -------------------------------------------------------------------
r = tidy.tidy(run)
check(not r.get("error") and r.get("packed") == len(files),
      "tidy did not pack what the plan named: %s" % r)
arc = J("archive", r.get("archive") or "?")
check(os.path.isfile(arc) and os.path.isfile(J("archive", "tidy.json")),
      "no archive and record were written")
check(not os.path.exists(arc + ".part"), "the part-written archive was left")
check(r["inodes_after"] < r["inodes_before"] - len(files) + 3,
      "the run does not hold fewer files: %s" % r)
gone = [f for f in files if os.path.exists(J(f))]
check(not gone, "packed files are still in the run: %s" % gone[:4])
check(not os.path.exists(J("docking", B, "out")),
      "ADCP 1.1's emptied working folder was left behind")
check(all(os.path.isdir(J("mmgbsa", p, "pose%d" % k)) for p in (A, B, C)
          for k in range(1, 6)) and os.path.isdir(J("docking", A)),
      "a peptide's own folder went with its files")
missing = [f for f in LOAD_BEARING if not os.path.isfile(J(f))]
check(not missing, "a load-bearing file went: %s" % missing[:4])
after_reads = readings(run)
for k in sorted(before_reads):
    check(after_reads[k] == before_reads[k],
          "what the page reads for %r changed once the run was tidied" % k)
check(after_reads["results"]["rows"] == before_reads["results"]["rows"]
      and len(after_reads["results"]["rows"]) == 3,
      "the results page does not return the same rows once tidied")
viewer = [J("receptor", "prepared.pdb")] + [
    J("poseviewer", n) for n in (os.listdir(J("poseviewer"))
                                 if os.path.isdir(J("poseviewer")) else [])]
check(all(os.path.isfile(p) for p in viewer),
      "a file the viewer serves is gone")
dl = set(r_ for r_, _p in bundle.walk(run, bundle.RESULTS))
lost = results_dl - dl
check(all(os.path.basename(f).startswith("out_") or os.path.basename(f) == "out.pdb"
          for f in lost) and all("docking" in f for f in lost),
      "the results download lost more than replica trajectories: %s"
      % sorted(lost)[:4])
check(any(f.startswith("archive" + os.sep) for f, _p in bundle.walk(run)),
      "the download of everything does not carry the archive")
st = tidy.status(run)
check(st["tidied"] and st["packed"] == len(files) and st["inodes_after"],
      "status does not say the run is tidied: %s" % st)
check(app.run_status(run).get("tidy", {}).get("tidied"),
      "the run's payload does not carry its tidy status")
again = tidy.tidy(run)
check(again.get("packed") == 0 and again.get("note") == "nothing to pack",
      "a second tidy finds something: %s" % again)
s = tidy.submit(run)
check("error" in s and "nothing to tidy" in s["error"],
      "the page offers to tidy a run with nothing left to pack: %s" % s)

# ---- undo -----------------------------------------------------------------------
u = tidy.undo(run)
check(u.get("restored") == len(files) and u.get("kept") == 0,
      "undo does not restore every file: %s" % u)
now = snapshot(run)
check(now == before, "undo does not put every file back as it was: %s"
      % sorted(k for k in set(now) | set(before) if now.get(k) != before.get(k))[:6])
check(stat.S_IMODE(os.stat(J("docking", A, "receptor.trg")).st_mode) == 0o640
      and stat.S_IMODE(os.stat(J("docking", B, "receptor.trg")).st_mode) == 0o664,
      "undo changed the files' permissions")
check(not os.path.exists(J("archive")), "the archive folder was left")
check(not tidy.status(run)["tidied"], "status still says the run is tidied")
check(readings(run) == before_reads, "what the page reads differs after Undo")

# A file written again after the tidy, as a rescoring writes its working
# files, is newer than the archived copy and stays.
tidy.tidy(run)
redo = J("mmgbsa", A, "pose1", A + "_min.in")
write(redo, "written again by a rescoring\n")
u = tidy.undo(run)
check(u.get("kept") == 1 and open(redo).read() == "written again by a rescoring\n",
      "undo overwrote a file written since the tidy: %s" % u)
write(redo, "_min.in\n")
check(snapshot(run) == before, "the run is not back to its baseline")

# ---- the job's own path: what the Slurm job runs, without Slurm ------------------
job = J("job.json")
rc = worker.main([job, "tidy"])
check(rc == 0 and tidy.status(run)["tidied"],
      "the worker's tidy stage does not tidy the run (rc %s)" % rc)
rc = worker.main([job, "untidy"])
check(rc == 0 and snapshot(run) == before,
      "the worker's untidy stage does not put it all back (rc %s)" % rc)
rc = worker.main([os.path.join(nonrun, "job.json"), "tidy"])
check(rc == 1, "the worker's tidy stage of a folder that is not a run exits %s" % rc)


def settle(rundir, most=60):
    for _i in range(most * 10):
        if not tidy.job_state(rundir).get("live"):
            return True
        time.sleep(0.1)
    return False


# The local executor: no scheduler, so the page's own process does it.
os.environ["ADCP_EXECUTOR"] = "local"
try:
    s = tidy.submit(run)
    check(s.get("local") and s.get("job"), "a local tidy was not started: %s" % s)
    check(settle(run) and tidy.status(run)["tidied"]
          and tidy.job_state(run).get("state") == "COMPLETED",
          "the local tidy did not finish: %s" % tidy.job_state(run))
    s = tidy.submit(run, undo_=True)
    check(settle(run) and not tidy.status(run)["tidied"],
          "the local undo did not finish: %s" % tidy.job_state(run))
finally:
    os.environ.pop("ADCP_EXECUTOR", None)
    if saved_env.get("ADCP_EXECUTOR"):
        os.environ["ADCP_EXECUTOR"] = saved_env["ADCP_EXECUTOR"]
os.remove(J(tidy.JOBFILE))
check(snapshot(run) == before, "the local tidy and undo did not leave the run as "
      "it was")

# ---- a Slurm submission, with sbatch and the queue stood in for -----------------
real = (slurm.available, slurm.submit, slurm.state, slurm.local_mode)
asked = []
slurm.available = lambda: True
slurm.submit = lambda script, **k: asked.append(script) or "12345"
slurm.state = lambda jid, rundir=None: {"RUNNING": 1} if str(jid) == "12345" else {}
slurm.local_mode = lambda: False
try:
    s = tidy.submit(run)
    check(s.get("job") == "12345" and asked == [J(tidy.SCRIPT)],
          "the tidy was not submitted as a job of its own: %s" % s)
    text = open(J(tidy.SCRIPT)).read() if os.path.isfile(J(tidy.SCRIPT)) else ""
    check("adcp_pilot.worker" in text and text.rstrip().endswith(" tidy")
          and "--mem=2G" in text and "--cpus-per-task=1" in text,
          "the tidy job's script is not a short single core job running the "
          "worker's tidy stage:\n%s" % text)
    rec = json.load(open(J(tidy.JOBFILE)))
    check(rec.get("job_id") == "12345" and "12345" not in json.dumps(
        json.load(open(job)).get("jobs")),
          "the tidy job is not recorded apart from job.json: %s" % rec)
    js = tidy.job_state(run)
    check(js.get("live") and app.run_status(run)["tidy"]["job"].get("live"),
          "a queued tidy does not read as live: %s" % js)
    s2 = tidy.submit(run)
    check("error" in s2 and "already" in s2["error"],
          "a second tidy is accepted while one is queued: %s" % s2)
    for what, said in (("rescore", app.rescore_run(run, "amber")),
                       ("panel", app.panel_run(run)),
                       ("resume", app.resume_run(run))):
        check("error" in said and ("tidied" in said["error"]
                                   or "not installed" in said["error"]),
              "a %s is accepted while the run is being tidied: %s" % (what, said))
finally:
    slurm.available, slurm.submit, slurm.state, slurm.local_mode = real
for f in (tidy.SCRIPT, tidy.JOBFILE):
    os.remove(J(f))
check(snapshot(run) == before, "the submission touched the run beyond its "
      "script and record")

# ---- the page's route and payload ------------------------------------------------
src = open(os.path.join(ROOT, "adcp_pilot", "server.py")).read()
check('u.path == "/api/tidy"' in src and "def tidy_run(self, rundir" in src,
      "the server does not route POST /api/tidy to App.tidy_run")
appjs = open(os.path.join(ROOT, "adcp_pilot", "static", "app.js")).read()
check("api('/api/tidy'" in appjs and "'Tidying…'" in appjs
      and ">Undo tidy</button>" in appjs and "Tidy</button>" in appjs,
      "the run's page has no Tidy and Undo tidy button wired to /api/tidy")
check("confirm(" not in appjs[appjs.index("function tidyLine"):
                              appjs.index("async function openRun")],
      "Tidy asks for confirmation, though Undo reverses it")
heard = []
real_submit = tidy.submit
tidy.submit = lambda d, undo_=False: heard.append((d, undo_)) or {"job": "7"}
httpd = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
httpd.app = app
httpd.token = app.token
threading.Thread(target=httpd.serve_forever, daemon=True).start()
opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
try:
    for undo in (False, True):
        req = urllib.request.Request(
            "http://127.0.0.1:%d/api/tidy" % httpd.server_address[1],
            data=json.dumps({"rundir": run, "undo": undo}).encode(),
            headers={"Content-Type": "application/json", "X-Token": app.token})
        with opener.open(req, timeout=30) as resp:
            body = json.loads(resp.read())
        check(body == {"job": "7"}, "POST /api/tidy answered %s" % body)
    check(heard == [(run, False), (run, True)],
          "POST /api/tidy did not reach tidy.submit as asked: %s" % heard)
    with opener.open("http://127.0.0.1:%d/api/run?dir=%s&t=%s"
                     % (httpd.server_address[1], urllib.request.quote(run),
                        app.token), timeout=30) as resp:
        payload = json.loads(resp.read())
    check("tidy" in payload and payload["tidy"].get("tidied") is False,
          "GET /api/run does not carry the run's tidy status: %s"
          % payload.get("tidy"))
finally:
    tidy.submit = real_submit
    httpd.shutdown()
check(tidy.main([run, "--dry-run"]) == 0 and snapshot(run) == before,
      "the command line's dry run fails or changes the run")
check(tidy.main([run, "--nonsense"]) == 2, "an unknown option is not refused")

# ---- the button, in each state, under Node where it is found ------------------
NODE = os.environ.get("ADCP_NODE") or shutil.which("node")
if NODE and os.path.isfile(NODE if os.sep in NODE else shutil.which(NODE) or ""):
    parts = [re.search(r"const esc = .*?\}\[c\]\)\);", appjs, re.S).group(0),
             re.search(r"const TIDY_BUSY = [^;]*;", appjs).group(0),
             "const TIDYMSG = {};"]
    i = appjs.index("function tidyLine(")
    parts.append(appjs[i:appjs.index("\n}\n", i) + 3])
    run_ = {"dock": {"id": "1", "state": {"COMPLETED": 3}}}
    busy_ = {"dock": {"id": "1", "state": {"RUNNING": 1, "COMPLETED": 2}}}
    done_ = {"tidied": True, "packed": 3185, "inodes_before": 3879,
             "inodes_after": 694, "job": {"job_id": "9", "state": "COMPLETED"}}
    cases = [{"done": 3, "jobs": run_, "tidy": {"tidied": False}},
             {"done": 3, "jobs": busy_, "tidy": {"tidied": False}},
             {"done": 0, "jobs": {}, "tidy": {"tidied": False}},
             {"done": 3, "jobs": run_, "tidy": {"job": {"live": True}}},
             {"done": 3, "jobs": run_, "tidy": {"tidied": True,
                                                "job": {"live": True, "undo": True}}},
             {"done": 3, "jobs": run_, "tidy": done_},
             {"done": 3, "jobs": busy_, "tidy": done_},
             {"done": 3, "jobs": run_, "tidy": {"tidied": False, "job": {
                 "job_id": "9", "state": "FAILED", "live": False}}}]
    prog = ("\n".join(parts) + "\nconst cases = JSON.parse(process.argv[1]);"
            "\nprocess.stdout.write(JSON.stringify(cases.map(c => "
            "tidyLine('/r', c))));\n")
    out = subprocess.run([NODE, "-e", prog, json.dumps(cases)],
                         capture_output=True, text=True, timeout=60)
    check(out.returncode == 0, "tidyLine could not be run under Node: %s"
          % out.stderr[-400:])
    h = json.loads(out.stdout or "[]") if out.returncode == 0 else [""] * 8
    check(">Tidy</button>" in h[0] and "disabled" not in h[0]
          and "data-undo" not in h[0] and "quota counts files" in h[0],
          "a finished, idle run is not offered Tidy with its tooltip: %s" % h[0])
    check(h[1] == "" and h[2] == "",
          "Tidy is offered while a job runs, or with no peptide finished: %s" % h[1:3])
    check("disabled>Tidying…</button>" in h[3],
          "a running tidy does not show a disabled Tidying…: %s" % h[3])
    check("disabled>Putting the files back…</button>" in h[4],
          "a running undo does not say so: %s" % h[4])
    check("data-undo=1" in h[5] and ">Undo tidy</button>" in h[5]
          and " disabled" not in h[5] and "3,185 files packed" in h[5],
          "a tidied run is not offered Undo tidy with its count: %s" % h[5])
    check(">Undo tidy</button>" in h[6] and " disabled" in h[6],
          "Undo is offered while the run's jobs run: %s" % h[6])
    check(">Tidy</button>" in h[7] and "did not finish (failed)" in h[7],
          "a failed tidy job is not said: %s" % h[7])
else:
    print("note  no Node here (ADCP_NODE or node on the PATH), so the button's "
          "states were not run")

# ---- the run's own tidy, queued with it ------------------------------------
from adcp_pilot import job as jobmod  # noqa: E402

armrun = tempfile.mkdtemp(prefix="adcp_arm_")
jobmod.save(armrun, jobmod.default())
asked, cancelled = [], []
saved = (slurm.submit, slurm.cancel, slurm.state, slurm.local_mode, slurm.available,
         tidy.live_jobs)
slurm.submit = lambda script, depends_on=None, kind="afterok": asked.append(
    (depends_on, kind)) or "777"
slurm.cancel = lambda jid: cancelled.append(jid)
slurm.state = lambda jid, rundir=None: {"PENDING": 1}
slurm.local_mode = lambda: False
slurm.available = lambda: True
tidy.live_jobs = lambda d: ["400"]
got = tidy.arm(armrun, also=["500", "501"])
check(got.get("armed") and asked[-1] == ("500:501:400", "afterany"),
      "a run queues its own tidy behind every job it has, however they end: %s %s"
      % (got, asked))
rec = json.load(open(os.path.join(armrun, tidy.JOBFILE)))
check(rec.get("auto") and rec.get("after") == "500:501:400",
      "and records it as its own: %s" % rec)
waiting = tidy.job_state(armrun)
check(tidy.waiting_auto(waiting) and server.App._tidying(armrun) is None,
      "while it waits, a rescoring or a panel may still be queued: %s" % waiting)
tidy.arm(armrun, also=["600"])
check(cancelled == ["777"] and asked[-1] == ("600:400", "afterany"),
      "a tidy still waiting is replaced by one behind the new job too: %s" % cancelled)
tidy.live_jobs = lambda d: []
n = len(asked)
got = tidy.submit(armrun)
check("nothing to tidy" in got.get("error", "") and cancelled == ["777"] and len(asked) == n,
      "Tidy pressed with nothing to pack leaves the waiting one alone: %s" % got)
j = jobmod.load(armrun)
j["keep_files"] = True
jobmod.save(armrun, j)
got = tidy.arm(armrun, also=["700"])
check(not got.get("armed") and len(asked) == n, "with every file kept, none is queued: %s" % got)
j["keep_files"] = False
jobmod.save(armrun, j)
slurm.local_mode = lambda: True
got = tidy.arm(armrun, also=["local:dock"])
check(not got.get("armed") and len(asked) == n,
      "nor on a computer with no scheduler, where Tidy is pressed by hand: %s" % got)
(slurm.submit, slurm.cancel, slurm.state, slurm.local_mode, slurm.available,
 tidy.live_jobs) = saved
src = open(os.path.join(ROOT, "adcp_pilot", "server.py")).read()
check(src.count("self._arm(") == 4 and '"keep_files"' in open(
    os.path.join(ROOT, "adcp_pilot", "job.py")).read(),
      "a new run, a rescoring, a panel and added measures all queue the run's tidy, "
      "and the form's Keep every file reaches the job")
check('"--keep-files"' in open(os.path.join(ROOT, "adcp_pilot", "cli.py")).read()
      and jobmod.default(keep_files=True)["keep_files"] is True,
      "the command line takes --keep-files, as a flow's card sends it, and the job keeps it")
shutil.rmtree(armrun, ignore_errors=True)

shutil.rmtree(work, ignore_errors=True)
if fails:
    print("test_tidy FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  a tidied run packs only its intermediates, every reader gives what it "
      "gave before, and Undo restores it byte for byte, %d checks" % ran)
