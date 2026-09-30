"""ADCP re-docking of stitched constructs.

AutoDock CrankPep folds the peptide inside the receptor field instead of docking a
fixed conformer, so a construct is not required to reproduce the poses of the two
anchors it was built from. This is why a pair whose parent poses overlap can still
be worth testing, and why the constructs are re-ranked here rather than inheriting
the anchors' scores.

IMPORTANT, which output file to read.
  out_<n>.pdb          per-replica trajectory in ADCP's REDUCED representation.
                       Side chains are truncated at CG; there are no aromatic ring
                       atoms. Contact counts computed on these files are wrong.
  out_ranked_<n>.pdb   the all-atom ranked poses. Always use these.
"""
from __future__ import annotations
import concurrent.futures as cf
import glob, os, re, shutil, subprocess, threading
import numpy as np
from .receptor import is_hydrogen

_MODE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+\S+\s+\S+\s+(\d+)")

# --- the engine contract ------------------------------------------------------
# Each docking engine module provides the same three things, so that choosing
# an engine is a config value rather than a different code path.

def engine_available(tools=None):
    """Whether this machine can run ADCP."""
    if tools is not None:
        return bool(getattr(tools, "adcp", None) and getattr(tools, "agfr", None))
    return bool(shutil.which("adcp") and shutil.which("agfr"))


def engine_prepare(rec_pdb, rec_pdbqt, center, size, workdir, tools, log=print):
    """ADCP needs grid maps built over the site before anything can dock."""
    return build_target(rec_pdbqt, center, size, workdir, tools)


def build_target(rec_pdbqt: str, center, size, workdir: str, tools, name="target") -> str:
    """agfr grid maps -> .trg consumed by adcp."""
    os.makedirs(workdir, exist_ok=True)
    trg = os.path.join(workdir, f"{name}.trg")
    if os.path.isfile(trg):
        return trg
    cmd = [tools.agfr, "-r", os.path.abspath(rec_pdbqt), "-o", name,
           "-b", "user", f"{center[0]:.3f}", f"{center[1]:.3f}", f"{center[2]:.3f}",
           f"{size[0]:.1f}", f"{size[1]:.1f}", f"{size[2]:.1f}"]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
    if not os.path.isfile(trg):
        cand = sorted(glob.glob(os.path.join(workdir, f"{name}*.trg")))
        if cand:
            return cand[0]
        raise RuntimeError(f"agfr failed to build a target:\n{p.stdout}\n{p.stderr}")
    return trg

def dock(sequence: str, trg: str, workdir: str, tools,
         replicas: int = 100, steps: int = 5_000_000, cpus: int = 0, log=print):
    """Run ADCP on one construct. Returns (modes, workdir)."""
    wd = os.path.join(workdir, sequence)
    os.makedirs(wd, exist_ok=True)
    logf = os.path.join(wd, "adcp.log")
    # ADCP builds internal paths with os.path.join('./tmp_<job>/', target[:-4], ...).
    # An absolute -t path makes os.path.join discard the './tmp_' prefix, so the
    # unpacked maps are looked for in the wrong directory and every run fails.
    # Copy the target next to the job and pass a bare basename.
    local_trg = os.path.join(wd, os.path.basename(trg))
    # The copy is refreshed whenever it differs from the target it was
    # made from, not only when it is missing. An interrupted run can
    # leave a stale or truncated copy here, and copying only when absent
    # would keep docking into the old maps after the target is rebuilt.
    src = os.path.abspath(trg)
    if (not os.path.isfile(local_trg)
            or os.path.getsize(local_trg) != os.path.getsize(src)):
        shutil.copy(src, local_trg)
    # Completion is marked by a file this code writes rather than by
    # adcp.log, which is created when the run starts and would accept a
    # run interrupted part way through.
    marker = os.path.join(wd, ".adcp_finished")
    done = os.path.isfile(marker) or os.path.isfile(
        os.path.join(wd, "out_ranked_1.pdb"))
    if not done:
        env = dict(os.environ)
        env.pop("OMP_NUM_THREADS", None)      # ADCP deadlocks if this is set
        cmd = [tools.adcp, "-t", os.path.basename(trg), "-s", sequence.lower(),
               "-N", str(replicas), "-n", str(steps), "-O", "-o", "out"]
        if cpus:
            cmd += ["-c", str(cpus)]
        with open(logf, "w") as fh:
            r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                               cwd=wd, env=env)
        # Only when ADCP exited cleanly. A killed run leaves no
        # marker and is run again.
        if r.returncode == 0:
            try:
                with open(marker, "w") as fh:
                    fh.write(f"{replicas} replicas, {steps} steps\n")
            except OSError:
                pass
        if not os.path.isfile(os.path.join(wd, "out_ranked_1.pdb")):
            # Why there is no pose, when ADCP says so itself.
            #
            # "No pose" can mean that the search ran and nothing
            # scored well enough, or that ADCP failed. ADCP names
            # the failure in its own output, for example:
            #
            #     Building protein from sequence: nwr (3 amino
            # acids, 1 chains).
            #     ERROR! Something has gone wrong.  Maybe too
            # short chains?
            #
            # A peptide shorter than five residues fails this way
            # in every replica. Reporting a failure as an empty
            # result would hide it, so the reason is passed on.
            _why = _failure_reason(wd)
            if _why:
                log(f"      {sequence}: ADCP could not dock this at all: "
                    f"{_why}")
            else:
                log(f"      {sequence}: ADCP ranked no pose; the search ran and "
                f"clustered nothing above its threshold")
    return read_modes(wd), wd

# Running several dockings at once.
#
# ADCP runs one Monte Carlo replica per process, as many at a time as it has
# cores, so a single docking stops getting faster once one wave holds every
# replica. Measured at 20 replicas: 127 s on one core, 21 s on eight, 17 s on
# sixteen and 14 s on forty-eight.
#
# The dockings are independent, so a large allocation is used better by
# running several side by side, each with a moderate number of cores.

from .estimate import CORES_PER_DOCKING, plan_concurrency   # noqa: E402


def dock_many(sequences, trg, workdir, tools, replicas=100, steps=5_000_000,
              cpus=0, log=print, per_docking=CORES_PER_DOCKING):
    """Dock every sequence, several at a time. Returns {sequence: (modes, wd)}.

    Order of completion is not the order given, so each line says which
    sequence it belongs to; the results are keyed by sequence rather than
    returned in a list, so the caller cannot accidentally depend on the order.
    """
    seqs = list(dict.fromkeys(sequences))          # keep order, drop repeats
    if not seqs:
        return {}
    workers, cores = plan_concurrency(cpus, len(seqs), per_docking,
                                      replicas=replicas)
    if workers == 1:
        out = {}
        for s in seqs:
            out[s] = dock(s, trg, workdir, tools, replicas, steps, cores, log)
        return out

    log(f"    docking {len(seqs)} sequence(s), {workers} at a time on "
        f"{cores} core(s) each")
    lock = threading.Lock()
    out, done = {}, [0]

    def one(seq):
        r = dock(seq, trg, workdir, tools, replicas, steps, cores, lambda *a: None)
        with lock:
            done[0] += 1
            n = len(r[0]) if r and r[0] else 0
            # One line per docking as it lands, naming the sequence: they do
            # not finish in the order they were started.
            log(f"      [{done[0]}/{len(seqs)}] {seq}: "
                + (f"{n} pose(s)" if n else "no poses"))
        return seq, r

    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        for seq, r in pool.map(one, seqs):
            out[seq] = r
    return out


def _failure_reason(wd):
    """What ADCP said went wrong, if it said anything, in one line.

    Its per-replica .out files carry the reason. Reading one is enough: when
    ADCP cannot build the peptide at all, every replica fails the same way.
    """
    import glob as _glob
    # The replica tally first. ADCP writes "N ENDED WITH ERROR" into its
    # own log for each replica that failed, and that is the only record
    # of it: the per-replica .out files carry no error line, so reading
    # only those would report failed replicas as a search that found
    # nothing above its threshold.
    logf = os.path.join(wd, "adcp.log")
    n_err = 0
    try:
        with open(logf, errors="replace") as fh:
            text = fh.read()
        n_err = text.count("ENDED WITH ERROR")
    except OSError:
        text = ""
    if n_err:
        n_res = None
        for line in text.splitlines():
            if "amino acids" in line and "Building protein" in line:
                n_res = line.split("(")[-1].split()[0]
        detail = ("every replica failed" if n_err >= 20
                  else "%d replica(s) failed" % n_err)
        if n_res:
            return "%s; it was given %s residues" % (detail, n_res)
        return detail

    for path in sorted(_glob.glob(os.path.join(wd, "out_*.out")))[:3]:
        try:
            text = open(path, errors="replace").read()
        except OSError:
            continue
        for line in text.splitlines():
            if line.strip().startswith("ERROR!"):
                short = line.strip().lstrip("ERROR!").strip()
                n = None
                for l2 in text.splitlines():
                    if "amino acids" in l2 and "Building protein" in l2:
                        n = l2.split("(")[-1].split()[0]
                if n and "short" in short.lower():
                    return ("%s (it was given %s residues)" % (short, n))
                return short
    return ""


def read_modes(wd: str) -> list[dict]:
    logf = os.path.join(wd, "adcp.log")
    if not os.path.isfile(logf):
        return []
    out = []
    for l in open(logf):
        m = _MODE.match(l)
        if m:
            out.append(dict(rank=int(m.group(1)), affinity=float(m.group(2)),
                            cluster=int(m.group(3))))
    return out

def read_pose(wd: str, rank: int = 1) -> dict | None:
    """All-atom ranked pose. See the module docstring for why out_ranked_* is used."""
    p = os.path.join(wd, f"out_ranked_{rank}.pdb")
    if not os.path.isfile(p):
        return None
    xyz, resnum, resname, name = [], [], [], []
    for l in open(p):
        if not l.startswith("ATOM"):
            continue
        nm = l[12:16].strip()
        el = l[76:78].strip()
        if is_hydrogen(nm, el):
            continue
        xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        resnum.append(int(l[22:26])); resname.append(l[17:20].strip()); name.append(nm)
    if not xyz:
        return None
    return dict(xyz=np.array(xyz), resnum=np.array(resnum),
                resname=resname, name=name)
