"""Submitting the computational stages as batch jobs.

By Aadhil Haq

Nothing here runs a docking directly. The login node of a shared cluster is
for setting a run up and reading it back; a tool that starts a multi-core
search there risks the user's account.
"""
import os
import re
import shutil
import subprocess

ARRAY = """#!/bin/bash
#SBATCH --job-name={name}
{partition}#SBATCH --nodes=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --array={indices}{throttle}
#SBATCH --output={logdir}/{name}_%a_%A.out
{extra}
cd {workdir}
{env}
{python} -m adcp_dock.worker {jobfile} ${{SLURM_ARRAY_TASK_ID}}
"""

SINGLE = """#!/bin/bash
#SBATCH --job-name={name}
{partition}#SBATCH --nodes=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}
#SBATCH --time={time}
#SBATCH --output={logdir}/{name}_%j.out
{extra}
cd {workdir}
{env}
{python} -m adcp_dock.worker {jobfile} {stage}
"""


def partitions():
    """Every partition, its time limit in minutes, and which is default.

    Returns (rows, default) where a row is {"name", "limit"} and `limit` is
    None for UNLIMITED or anything sinfo will not parse.
    """
    try:
        out = subprocess.run(["sinfo", "-h", "-o", "%P|%l"],
                             capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return [], ""
    rows, default, seen = [], "", set()
    for line in (out.stdout or "").splitlines():
        if "|" not in line:
            continue
        n, _, lim = line.partition("|")
        n, lim = n.strip(), lim.strip()
        if not n:
            continue
        if n.endswith("*"):
            n = n[:-1]
            default = n
        if n in seen:
            continue
        seen.add(n)
        rows.append({"name": n, "limit": _limit_minutes(lim)})
    return rows, default or (rows[0]["name"] if rows else "")


def _limit_minutes(text):
    """SLURM's "3-00:00:00" or "2:00:00" as minutes, or None."""
    text = (text or "").strip()
    if not text or text.lower().startswith(("infinite", "unlimited", "n/a")):
        return None
    days = 0
    if "-" in text:
        d, _, text = text.partition("-")
        try:
            days = int(d)
        except ValueError:
            return None
    parts = text.split(":")
    try:
        nums = [int(x) for x in parts]
    except ValueError:
        return None
    if len(nums) == 3:
        h, m, sec = nums
    elif len(nums) == 2:
        h, m, sec = 0, nums[0], nums[1]
    elif len(nums) == 1:
        h, m, sec = 0, nums[0], 0
    else:
        return None
    return days * 1440 + h * 60 + m + (1 if sec else 0)


def target_partition(partition="", extra=""):
    """The partition a run's jobs land in: the one named on the run, or in its
    other Slurm options, or else the cluster's default."""
    name = (partition or "").strip()
    if not name:
        found = re.search(r"(?:--partition=|-p\s*)(\S+)", extra or "")
        name = found.group(1) if found else ""
    return name or partitions()[1]


def within_limit(hhmmss, partition="", extra=""):
    """`hhmmss`, or the partition's own time limit when that is shorter.

    For the stages whose length this package chooses rather than the person:
    the preparation asks for four hours to be safe and takes minutes, and on a
    cluster whose default partition allows two a four-hour request is refused
    before anything has run. A limit that cannot be read leaves the request
    as it is.
    """
    rows, _default = partitions()
    name = target_partition(partition, extra)
    limit = next((r["limit"] for r in rows if r["name"] == name), None)
    asked = _limit_minutes(hhmmss)
    if limit is None or asked is None or asked <= limit:
        return hhmmss
    return "%d:%02d:00" % (limit // 60, limit % 60)


def fits(partition, minutes):
    """Whether `minutes` fits inside that partition's time limit.

    None when the partition is unknown here: absence of evidence, and the
    caller should not refuse on it.
    """
    rows, _d = partitions()
    for r in rows:
        if r["name"] == partition:
            return True if r["limit"] is None else int(minutes) <= r["limit"]
    return None


def node_cores():
    """The cores on one node of the partition a job lands in when it names none, or None.

    The default is the largest node of the default partition, because that is
    where a job with no partition goes; a cluster without one marked default
    gives the largest node it has. None when sinfo is not here or says
    nothing, which leaves the caller's own default alone.
    """
    try:
        out = subprocess.run(["sinfo", "-h", "-o", "%P|%c"], capture_output=True,
                             text=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    default, any_ = [], []
    for line in out.splitlines():
        name, _, cpus = line.strip().partition("|")
        found = re.match(r"(\d+)", cpus.strip())
        if not found:
            continue
        any_.append(int(found.group(1)))
        if name.endswith("*"):
            default.append(int(found.group(1)))
    pick = default or any_
    return max(pick) if pick else None


def cores_for_cluster(wanted):
    """`wanted` cores, or fewer when no node of the default partition has that many.

    A task asking for more cores than a node has is refused by sbatch, so a
    default written for one cluster's 48-core nodes would stop every run on a
    cluster with 32.
    """
    have = node_cores()
    return min(int(wanted), have) if have else int(wanted)


def tool_exports(roots=()):
    """Lines for a job script that name where the programs are, read now on the login node.

    A job usually inherits the environment it was submitted from, which is
    where ADFRSUITE_HOME and AMBERHOME come from. Some clusters start every
    job with a clean environment instead, and there the job would not find
    ADFRsuite at all; so the script says it outright, the same way on both.
    """
    import shlex
    from . import tools as toolmod
    found = toolmod.discover(roots=roots)
    lines = []
    for var in ("ADFRSUITE_HOME", "AMBERHOME", "ADCPSUITE_HOME", "ADAPT_HOME",
                "ADCP_DOCK_TOOLS", "ADCP_ACKNOWLEDGEMENT"):
        value = os.environ.get(var)
        if value:
            lines.append("export %s=%s" % (var, shlex.quote(value)))
    if not os.environ.get("ADFRSUITE_HOME") and found.get("adcp"):
        lines.append("export ADFRSUITE_HOME=%s" % shlex.quote(
            os.path.dirname(os.path.dirname(found.get("adcp")))))
    if not os.environ.get("AMBERHOME"):
        amber = toolmod.amber_home(found)
        if amber:
            lines.append("export AMBERHOME=%s" % shlex.quote(amber))
    dirs = []
    for name in toolmod.NEEDED:
        where = found.get(name)
        if where and os.path.dirname(where) not in dirs:
            dirs.append(os.path.dirname(where))
    if dirs:
        lines.append("export PATH=%s:$PATH" % ":".join(shlex.quote(d) for d in dirs))
    return "\n".join(lines)


def why_pending(job_id, rundir=None):
    """What the scheduler says is holding a job up, in one word."""
    if is_local(job_id):
        return ""
    try:
        out = subprocess.run(["squeue", "-j", str(job_id), "-h", "-o", "%r"],
                             capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return ""
    return (out.stdout or "").strip().splitlines()[0].strip() \
        if (out.stdout or "").strip() else ""


def available():
    return shutil.which("sbatch") is not None


def index_spec(indices):
    """SLURM's array spec, collapsed into ranges. A resume runs a handful of
    scattered tasks (3,9,14) and spelling them out is how it only re-runs
    what is missing rather than the whole list again."""
    idx = sorted(set(int(i) for i in indices))
    if not idx:
        return ""
    runs, a, b = [], idx[0], idx[0]
    for i in idx[1:]:
        if i == b + 1:
            b = i
        else:
            runs.append((a, b))
            a = b = i
    runs.append((a, b))
    return ",".join("%d-%d" % r if r[0] != r[1] else str(r[0]) for r in runs)


def partition_line(name=""):
    """The --partition line, which is normally absent.

    Schedulers route a job to a partition from what it requests, and many
    clusters ask users not to name one. A partition is emitted only when the
    site requires it.
    """
    name = (name or "").strip()
    return "#SBATCH --partition=%s\n" % name if name else ""


import re as _re

# What a site-specific option may contain. These become #SBATCH lines in a
# bash script, so anything that could end the line or start a command is
# refused rather than escaped.
_SAFE = _re.compile(r"^-{1,2}[A-Za-z][A-Za-z0-9_-]*(=[A-Za-z0-9_=:,.@%+/-]+)?$")
# Options the tool sets itself and that must not be overridden from the form:
# two --output lines would put the log where the page does not look.
_OURS = {"--job-name", "-J", "--output", "-o", "--error", "-e", "--array",
         "-a", "--cpus-per-task", "-c", "--time", "-t", "--mem",
         "--nodes", "-N", "--dependency", "-d", "--parsable"}


def site_options(account="", extra=""):
    """(#SBATCH lines, [problems]) for a cluster's own requirements.

    Many clusters refuse a job with no --account, and some need a QOS, a
    constraint or a non-default partition. Those are properties of the site,
    not of the docking, so they are typed once, "--qos=normal
    --partition=short", and carried into every script this run writes.
    """
    lines, bad = [], []
    account = (account or "").strip()
    if account:
        if _SAFE.match("--account=" + account):
            lines.append("#SBATCH --account=%s" % account)
        else:
            bad.append("the account %r has characters that cannot go into a "
                       "job script" % account)
    for tok in (extra or "").split():
        name = tok.split("=", 1)[0]
        if not _SAFE.match(tok):
            bad.append("%r is not an sbatch option this can pass on; write "
                       "each as --name=value, e.g. --qos=normal" % tok)
        elif name in _OURS:
            bad.append("%s is set by the tool itself and cannot be changed "
                       "here" % name)
        else:
            lines.append("#SBATCH %s" % tok)
    return "\n".join(lines), bad


def write_script(path, **kw):
    kw["partition"] = partition_line(kw.get("partition"))
    kw.setdefault("extra", "")
    kw.setdefault("env", "")
    kw.setdefault("throttle", "")
    tmpl = ARRAY if "indices" in kw else SINGLE
    with open(path, "w") as fh:
        fh.write(tmpl.format(**kw))
    os.chmod(path, 0o755)
    return path


def submit(script, depends_on=None, kind="afterok"):
    """Returns the job id. `depends_on` chains a stage behind another."""
    cmd = ["sbatch", "--parsable"]
    if depends_on:
        cmd.append("--dependency=%s:%s" % (kind, depends_on))
    cmd.append(script)
    p = subprocess.run(cmd, capture_output=True, text=True)
    out = (p.stdout or "").strip().splitlines()
    for line in out:
        m = re.match(r"^(\d+)", line.strip())
        if m:
            return m.group(1)
    raise RuntimeError("sbatch did not return a job id: %s %s"
                       % (p.stdout.strip(), p.stderr.strip()))


def state(job_id, rundir=None):
    """{state: count} over the tasks of a job, from squeue then sacct.

    A local job (see LOCAL EXECUTION below) is read from the run directory
    instead, in the same shape, so nothing that summarises it can tell the
    difference.
    """
    if is_local(job_id):
        return local_state(job_id, rundir)
    p = subprocess.run(["squeue", "-j", str(job_id), "-h", "-r", "-o", "%T"],
                       capture_output=True, text=True)
    live = [l.strip() for l in p.stdout.splitlines() if l.strip()]
    if live:
        out = {}
        for s in live:
            out[s] = out.get(s, 0) + 1
        return out
    p = subprocess.run(["sacct", "-j", str(job_id), "-X", "-n", "-o", "State"],
                       capture_output=True, text=True)
    out = {}
    for line in p.stdout.splitlines():
        s = line.strip().split()[0] if line.strip() else ""
        if s:
            out[s] = out.get(s, 0) + 1
    return out


def finished(job_id, rundir=None):
    st = state(job_id, rundir)
    return bool(st) and not any(k in st for k in
                                ("PENDING", "RUNNING", "CONFIGURING",
                                 "COMPLETING", "REQUEUED", "SUSPENDED"))


# LOCAL EXECUTION
#
# The same stages, run on this machine instead of submitted to SLURM, for a
# notebook session or a workstation with no scheduler.
#
# Opt-in only, through ADCP_EXECUTOR=local; it is never selected because
# sbatch happens to be missing. On a cluster a missing sbatch usually means
# a login shell with an incomplete PATH, and falling back to local execution
# would start a docking on the login node.
#
# Nothing new is run: every SLURM task is `python -m adcp_dock.worker
# <job.json> <stage>`, and this runs those same commands in order.
#
# State is kept in the run directory (local_state.json) in the shape state()
# returns for a SLURM job, so the page's monitor reads it unchanged. The
# runner is a thread in the server; if the server exits mid-run, the recorded
# pid is no longer alive and anything left RUNNING or PENDING is reported as
# CANCELLED.

import json as _json
import sys as _sys
import threading as _threading

LOCAL_PREFIX = "L-"
LOCAL_STATE = "local_state.json"
# One run at a time. A second run submitted while one is in progress waits
# its turn, as it would in a queue: two dockings sharing the same cores
# finish later than one after the other.
_LOCAL_LOCK = _threading.Lock()


def stage_pythonpath():
    """What a stage needs on PYTHONPATH to import this package: its source checkout when running from one, and nothing when installed.

    For an installed package the parent directory is site-packages, the whole
    environment, which every program a stage starts would then inherit.
    """
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    if os.path.basename(here.rstrip(os.sep)) in ("site-packages",
                                                 "dist-packages"):
        return ""
    return here


def local_mode():
    return os.environ.get("ADCP_EXECUTOR", "").strip().lower() == "local"


def usable_cores():
    """The cores this process may actually run on.

    Not os.cpu_count(), which counts every core in the machine: inside a SLURM
    allocation or a container the process is pinned to a subset, and asking for
    the whole machine oversubscribes the few it has.
    """
    try:
        return len(os.sched_getaffinity(0)) or 1
    except (AttributeError, OSError):
        return os.cpu_count() or 1


def is_local(job_id):
    return str(job_id or "").startswith(LOCAL_PREFIX)


def _local_path(rundir):
    return os.path.join(rundir, LOCAL_STATE)


def _read_local(rundir):
    try:
        with open(_local_path(rundir)) as fh:
            return _json.load(fh)
    except (OSError, ValueError):
        return {}


def _write_local(rundir, st):
    # Written whole and renamed into place, so a page polling mid-write
    # never reads a partial file.
    tmp = _local_path(rundir) + ".tmp"
    with open(tmp, "w") as fh:
        _json.dump(st, fh, indent=1)
    os.replace(tmp, _local_path(rundir))


def _alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, TypeError, ValueError):
        return False


def local_state(job_id, rundir):
    """{state: count} for a local stage, e.g. {"COMPLETED": 3, "RUNNING": 1}."""
    if not rundir:
        return {}
    st = _read_local(rundir)
    tasks = (st.get("stages") or {}).get(str(job_id)[len(LOCAL_PREFIX):]) or {}
    dead = not _alive(st.get("pid"))
    out = {}
    for s in tasks.values():
        if dead and s in ("RUNNING", "PENDING"):
            s = "CANCELLED"
        out[s] = out.get(s, 0) + 1
    return out


def submit_local(rundir, jobfile, stages, cpus=1, python=None, env=None):
    """Run [(stage, [task, ...]), ...] here, in order, in the background.

    Returns {stage: job id}. Dependencies follow the SLURM submission: the dock
    tasks run only if prepare succeeded, and the report runs whatever the dock
    tasks did.
    """
    python = python or _sys.executable
    logdir = os.path.join(rundir, "logs")
    os.makedirs(logdir, exist_ok=True)
    st = {"pid": os.getpid(), "cores": cpus,
          "stages": {k: {t: "PENDING" for t in ts} for k, ts in stages}}
    _write_local(rundir, st)

    base = dict(os.environ)
    base.update(env or {})
    pp = stage_pythonpath()
    if pp:
        base["PYTHONPATH"] = pp + os.pathsep + base.get("PYTHONPATH", "")
    base["SLURM_CPUS_PER_TASK"] = str(cpus)
    base.setdefault("PYTHONUNBUFFERED", "1")

    def mark(stage, task, word):
        cur = _read_local(rundir) or st
        cur.setdefault("stages", {}).setdefault(stage, {})[task] = word
        _write_local(rundir, cur)

    def one(stage, task):
        mark(stage, task, "RUNNING")
        e = dict(base)
        if stage == "dock":
            e["SLURM_ARRAY_TASK_ID"] = task
        log = os.path.join(logdir, "local_%s%s.out"
                           % (stage, "" if task == stage else "_" + task))
        with open(log, "w") as fh:
            rc = subprocess.run([python, "-m", "adcp_dock.worker", jobfile,
                                 task], cwd=rundir, env=e, stdout=fh,
                                stderr=subprocess.STDOUT).returncode
        mark(stage, task, "COMPLETED" if rc == 0 else "FAILED")
        return rc == 0

    def run():
        with _LOCAL_LOCK:
            ok_so_far = True
            for stage, tasks in stages:
                if stage == "dock" and not ok_so_far:
                    for t in tasks:
                        mark(stage, t, "CANCELLED")
                    continue
                results = [one(stage, t) for t in tasks]
                if stage == "prepare":
                    ok_so_far = all(results)

    _threading.Thread(target=run, daemon=True,
                      name="adcp-local-%s" % os.path.basename(rundir)).start()
    return {k: LOCAL_PREFIX + k for k, _ts in stages}
