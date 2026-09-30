"""Running here instead of on Slurm: for Colab and for workstations.

By Aadhil Haq

Two halves. The first needs nothing: local mode is opt-in, a cluster
without sbatch still refuses, a dead runner is reported as dead, and the
peptide parser reads the forms the page says it reads.

The second is a real run, driven through the server as the page drives it:
create, poll, read the results. It needs ADFRsuite and a receptor, and
skips (exit 77) where either is missing. Only a real run exercises the
pipeline with the data types it actually passes between stages.
"""
import json
import os
import sys
import tempfile
import threading
import time
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import job as jobmod, peptide, pipeline, slurm   # noqa: E402

ENGINES = []
fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


# --- opt-in only -------------------------------------------------------------
#
# A missing sbatch on a cluster is almost always a login shell with an
# incomplete PATH, and falling back to local execution would start a docking
# on the login node. So local mode is chosen by name and never inferred.
_saved = os.environ.pop("ADCP_EXECUTOR", None)
check(not slurm.local_mode(), "local mode is on without being asked for")
os.environ["ADCP_EXECUTOR"] = "local"
check(slurm.local_mode(), "ADCP_EXECUTOR=local did not turn local mode on")
os.environ["ADCP_EXECUTOR"] = "LOCAL "
check(slurm.local_mode(), "local mode is fussy about case or whitespace")
os.environ.pop("ADCP_EXECUTOR", None)

# And without it, a machine with no sbatch still refuses to run anything.
_avail = slurm.available
slurm.available = lambda: False
try:
    d = tempfile.mkdtemp()
    j = jobmod.default()
    j.update({"receptor": "/nonexistent.pdb", "hotspots": "A:1",
              "peptides": [{"name": "x", "sequence": "AAAAA"}]})
    jobmod.save(d, j)
    try:
        pipeline.submit(d, j, log=lambda *a: None)
        check(False, "with no sbatch and no ADCP_EXECUTOR, submit() ran "
                     "something: on a cluster that is the login node")
    except SystemExit:
        check(True, "")
finally:
    slurm.available = _avail

# --- a runner that died is not reported as running ---------------------------
d = tempfile.mkdtemp()
slurm._write_local(d, {"pid": 2 ** 22 + 12345,       # not a live process
                       "stages": {"dock": {"0": "COMPLETED", "1": "RUNNING",
                                           "2": "PENDING"}}})
st = slurm.local_state("L-dock", d)
check(st.get("COMPLETED") == 1, "a finished task lost its state: %r" % st)
check(st.get("CANCELLED") == 2 and "RUNNING" not in st,
      "a dead runner's tasks still read as running or pending: %r" % st)
check(slurm.state("L-dock", d) == st,
      "slurm.state() does not route a local id to the run directory")
check(slurm.why_pending("L-dock", d) == "",
      "why_pending() asked squeue about a local job")

# --- cores: what this process may use, never more ----------------------------
n = slurm.usable_cores()
check(1 <= n <= (os.cpu_count() or n),
      "usable_cores() is out of range: %d" % n)
j = {"cpus": 10 ** 6, "mm_cpus": 10 ** 6, "mm_workers": 10 ** 6}
pipeline._fit_to_machine(j)
check(j["cpus"] == n and j["mm_cpus"] == n and j["mm_workers"] == n,
      "a job asking for more cores than this machine has was not cut down")

# --- the parser reads what the page says it reads ---------------------------
rows = peptide.parse_many(open(os.path.join(ROOT, "examples",
                                            "peptides.txt")).read())
check([r[0] for r in rows] == ["control_1", "control_2", "longer"]
      and not any(r[2] for r in rows),
      "examples/peptides.txt, the file the docs point at, does not parse "
      "cleanly: %r" % [(r[0], r[1], bool(r[2])) for r in rows])
r = peptide.parse_many("p53: ETFSDLWKLLP")[0]
check(r[0] == "p53" and r[1] == "ETFSDLWKLLP" and not r[2],
      "'name: SEQUENCE' is not read as a named peptide: %r" % (r,))
check(not peptide.parse_many("# a comment\n"),
      "a # comment line became a peptide")
r = peptide.parse_many("alpha SWQWWWQPQRDCR")[0]
check(r[2] and ":" in r[2][-1].message,
      "a name before a space is not refused with the colon form to use")

# --- a real run, through the server -----------------------------------------
from adcp_dock import tools as toolmod                          # noqa: E402

RECEPTOR = os.environ.get("ADCP_TEST_RECEPTOR",
                          os.path.join(ROOT, "examples", "receptor.pdb"))
T = toolmod.discover()
found = dict(T.report())
skip = None
if not all(found.get(b) for b in ("adcp", "agfr", "prepare_receptor")):
    skip = "ADFRsuite is not on this machine"
elif not os.path.isfile(RECEPTOR):
    skip = "no receptor at %s (set ADCP_TEST_RECEPTOR)" % RECEPTOR
elif os.environ.get("ADCP_SKIP_REAL_RUN"):
    skip = "ADCP_SKIP_REAL_RUN is set"
# Not on a login node. This half runs a real preparation and docking, and
# run_all.sh's own guard does not cover running this file directly. Inside a
# Slurm job it runs; on a workstation or in Colab, set ADCP_ALLOW_REAL_RUN=1;
# anywhere else it is skipped with a reason.
elif not (os.environ.get("SLURM_JOB_ID")
          or os.environ.get("ADCP_ALLOW_REAL_RUN")):
    skip = ("not inside a SLURM job, and a real docking does not run on a "
            "login node: set ADCP_ALLOW_REAL_RUN=1 on a workstation")

if skip is None:
    from adcp_dock import server
    os.environ["ADCP_EXECUTOR"] = "local"
    runs = tempfile.mkdtemp(prefix="adcp-local-")
    httpd = server.build(host="127.0.0.1", port=0, runs=runs)
    port, token = httpd.server_address[1], httpd.app.token
    threading.Thread(target=httpd.serve_forever, daemon=True).start()

    def call(path, body=None, **q):
        q["t"] = token
        url = "http://127.0.0.1:%d%s?%s" % (port, path, urllib.parse.urlencode(q))
        req = urllib.request.Request(
            url, data=json.dumps(body).encode() if body is not None else None,
            headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            return json.loads(resp.read().decode())

    tl = call("/api/tools")
    check(tl.get("local") is True,
          "the page is not told it is running locally: %r" % tl.get("local"))
    check(int(tl.get("cores") or 0) >= 1, "the page is not told the core count")

    # The receptor this is written for is 1YCR: MDM2 on chain A, a bound p53
    # peptide on chain B. Chain A only, and redock p53's core into the cleft
    # it came from. The budget is tiny: this proves the path runs, not that
    # the search is good.
    def run_local(body):
        """Create a run through the page's API, wait for it, check it, and
        return its result rows."""
        made = call("/api/create", body)
        check("error" not in made, "the run was not created: %r" % made)
        rundir = made.get("rundir", "")
        jobs = made.get("jobs") or {}
        check(set(jobs) == {"prepare", "dock", "report"}
              and all(str(v).startswith("L-") for v in jobs.values()),
              "local submission did not return local ids: %r" % jobs)
        t0, done, last = time.time(), False, {}
        while time.time() - t0 < 1500:
            st = call("/api/run", dir=rundir)
            last = {k: v.get("state") for k, v in (st.get("jobs") or {}).items()}
            live = any(s in ("RUNNING", "PENDING")
                       for c in last.values() for s in (c or {}))
            if last and not live:
                done = True
                break
            time.sleep(5)
        check(done, "the local run did not finish in 25 minutes: %r" % last)
        check(all(set(c) <= {"COMPLETED"} for c in last.values() if c),
              "a stage did not complete: %r; see %s/logs/local_*.out"
              % (last, rundir))
        # The preparation was complete, not merely finished. A stage that
        # could not protonate still exits 0 and the run still ranks peptides,
        # so the run's preparation warnings are checked too.
        st = call("/api/run", dir=rundir)
        check(st.get("prep_warnings") == [],
              "the preparation was degraded: %r" % st.get("prep_warnings"))
        for f in ("HOW_TO_CITE.txt", "CITATION.bib", "CITATION.ris"):
            check(os.path.isfile(os.path.join(rundir, f)),
                  "%s was not written into the run" % f)
        res = os.path.join(rundir, "results.json")
        check(os.path.isfile(res), "no results.json was written")
        if not os.path.isfile(res):
            return rundir, []
        rows = json.load(open(res))
        return rundir, (rows if isinstance(rows, list) else rows.get("rows", []))

    # The receptor this is written for is 1YCR: MDM2 on chain A, a bound p53
    # peptide on chain B. Chain A only, and redock p53's core into the cleft
    # it came from. The budget is tiny: this proves the path runs, not that
    # the search is good. ADCP 1.0 first, named, since "auto" would pick 1.1
    # wherever it is installed.
    body = {"label": "localtest", "receptor": RECEPTOR, "chains": ["A"],
            "hotspots": "A:54, A:62, A:67, A:96, A:99",
            "peptides": [{"name": "p53_core", "sequence": "FSDLWKLL"}],
            "replicas": 4, "steps": 400000, "poses": 3, "adcp_engine": "1.0",
            "cpus": min(4, slurm.usable_cores()), "mmgbsa": False}
    rundir, rows = run_local(body)
    ENGINES.append("ADCP 1.0")
    check(len(rows) == 1 and rows[0].get("name") == "p53_core",
          "the result does not carry its peptide's name: %r"
          % [(r.get("name"), r.get("peptide")) for r in rows])
    check(rows and rows[0].get("adcp_best") is not None
          and rows[0].get("engine") == "ADCP 1.0",
          "the ADCP 1.0 result has no score or no engine: %r" % rows[:1])

    # ADCP 1.1 where it is installed: a non-standard and a D residue, and the
    # top pose minimized with OpenMM.
    if tl.get("adcp11"):
        ENGINES.append("ADCP 1.1")
        body = dict(body, label="localtest11", adcp_engine="1.1",
                    omm_minimize=1, steps=200000,
                    peptides=[{"name": "p53_core", "sequence": "FSDLWKLL"},
                              {"name": "nonstd", "sequence": "FSN<BHD>LW&KLL"}])
        rundir, rows = run_local(body)
        by = {r.get("name"): r for r in rows}
        check(set(by) == {"p53_core", "nonstd"}
              and all(str(r.get("engine", "")).startswith("ADCP 1.1")
                      for r in rows),
              "the ADCP 1.1 run did not dock both peptides with 1.1: %r"
              % [(r.get("name"), r.get("engine")) for r in rows])
        check(by.get("nonstd", {}).get("length") == 8,
              "the non-standard peptide is not counted as 8 residues: %r"
              % by.get("nonstd", {}).get("length"))
        check(any(r.get("omm_n") for r in rows),
              "no OpenMM minimization result was reported: %r" % rows)
        poses = call("/api/poses", dir=rundir, peptide="FSN<BHD>LW&KLL")
        check(len(poses.get("poses") or []) >= 1,
              "the non-standard peptide's poses cannot be listed: %r" % poses)
    else:
        print("skip  ADCP 1.1 is not installed here; only ADCP 1.0 was run")
    httpd.shutdown()
    os.environ.pop("ADCP_EXECUTOR", None)

if _saved is not None:
    os.environ["ADCP_EXECUTOR"] = _saved

if fails:
    print("test_local FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
if skip:
    print("ok  local mode is opt-in and reads a dead runner as dead, %d checks; "
          "the real run was skipped: %s" % (ran, skip))
else:
    print("ok  real runs went prepare -> dock -> report on this machine "
          "through the server (%s), %d checks" % (" and ".join(ENGINES), ran))
