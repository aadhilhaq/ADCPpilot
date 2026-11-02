"""Tidy a finished run: the files its results no longer read are packed into
one archive inside the run and removed, and Undo puts them back.

By Aadhil Haq

The scratch file systems count files as well as bytes (Grace: 250 000 per
user). MM-GBSA leaves some 23 working files in every pose folder, and a run
made before keep_trajectories existed may still hold ADCP's working folder,
hundreds of files per peptide. Each was read once, by the next program, and
never again.

What is packed, and only for a peptide with a result file (results/<key>.json):

  docking/<key>/
    out/                                ADCP 1.1's working folder, where one
                                        was left (a new run packs or removes
                                        it itself, trajectories.finish)
    *.trg                               the copy of target/receptor.trg that
                                        ADCP docked from (adcp11.dock copies
                                        it in to run)
  mmgbsa/<name>/pose<n>/, <p> the folder's own name (mmforms.name)
    <p>.leap, <p>_rec.pdb, <p>_lig.pdb, <p>_rec_forleap.pdb,
    <p>_pep_forleap.pdb, <p>_mirror_pep.pdb, <p>_mirror.leap,
    <p>_mirror.parm7, <p>_mirror.rst7, <p>_com_tleap.rst7
                                        what tleap was given while the
                                        topologies were built
                                        (mmgbsa._topologies, _ring_leap,
                                        mirror_hydrogens)
    <p>_com.pdb                         the complex as first written, only
                                        when <p>_com_forleap.pdb is beside it
                                        (metrics.for_pose reads it only
                                        without one)
    <p>_com.pqr, <p>_com_reduce.pdb     what PDB2PQR and reduce wrote on the
                                        way to <p>_com_forleap.pdb
    <p>_min.in, mdinfo                  sander's input and progress file
    <p>_decomp.in, <p>_decomp_total.dat, reference.frc, _MMPBSA_*
                                        MMPBSA.py's input, the decomposition
                                        run's own total (the energy the table
                                        reports is <p>_mmgbsa.dat) and its
                                        leftovers

What stays, with what reads it:

  job.json                      everything, the run page and reuse settings
                                first (server.run_status, app.js reuseRun)
  results.json, results.tsv, results.txt, results_combined.tsv
                                the results page (server.run_status), the
                                re-ranking (server.rerank_run), the downloads
  results/<key>.json            whether a peptide is done (pipeline.
                                outstanding, so a resume never docks it
                                again), the report (stage_report), the
                                rescoring (stage_rescore), the pose picker
  metrics/                      the residue map (server.residues), the
                                re-ranking's pose spans (combined.poses_from)
  poseviewer/                   the structures the page and the downloads offer
  receptor/                     prepared.pdb, which the viewer, MM-GBSA and
                                the pose viewer files are built on
                                (receptor_info.prepared_pdb); prepared.json,
                                the key that stops a resumed preparation
                                preparing again (prepare.run);
                                preparation.log, which the results page reads
                                (server.prep_warnings). Kept whole: a handful
                                of files.
  target/                       receptor.trg, the grid maps every docking
                                reads and a resume checks for (pipeline.resume)
  logs/, *.sbatch, CITATION.*, HOW_TO_CITE.txt
                                the log tail on the page, the downloads
  docking/<key>/out_ranked_<n>.pdb
                                the poses: the viewer (server.structure_text),
                                the pose picker (server.pose_list), MM-GBSA
                                and its rescoring (adcp.read_pose), the pose
                                viewer files (stage_poseviewer)
  docking/<key>/adcp.log, modes.json, out_summary.dlg, out_out.pdb
                                the mode table and cluster sizes
                                (adcp11.read_modes, read_openmm,
                                pipeline.cluster_counts), the downloads
  docking/<key>/.adcp_finished, .trajectories_packed, <key>_trajectories.tar.gz
                                what says the docking and its packing are done
                                (adcp11.dock, trajectories.finish)
  mmgbsa/RADII, mmgbsa/<name>.dg, mmgbsa/<name>.<engine>.dg
                                the radii the run keeps (mmgbsa.radii_for) and
                                each engine's energies, which a rescoring
                                reuses rather than computes again
                                (mmgbsa.run) and the report reads
                                (mmgbsa.engines_on_disk)
  pose<n>/<p>_mmgbsa.dat, <p>_mmpbsa.in
                                the energy terms (pipeline.energy_terms,
                                metrics.for_pose), the downloads
  pose<n>/<p>_com.parm7, _rec.parm7, _lig.parm7, _min.rst7, _com.rst7,
  <p>_scored_on.txt, <p>_scored.pdb, <p>_decomp.dat, <p>_min.out,
  <p>_chains.json, <p>_checks.json
                                the interface panel, measured from these
                                without docking or minimizing again
                                (metrics.for_pose, pipeline.measurable, which
                                decides whether the page offers to measure it)
  pose<n>/<p>_com_forleap.pdb, <p>_com_forleap.json, <p>_com_prepped.pdb
                                the prepared complex and the key that lets a
                                rescoring reuse it (mmgbsa._prepared_for);
                                the residue numbering of the panel
                                (metrics.tleap_map); what Prime is given
                                (mmgbsa.prime_dg)
  pose<n>/<p>_strain.json       the pose's peptide strain and what it was
                                measured on (_vendor/strain.measure), read
                                back rather than measured again when strain
                                is added once more (pipeline.stage_strain);
                                it is measured from the topologies and the
                                scored frame above, which stay with it
  every log outside out/        the problem report (bundle.LOGS): leap.log,
                                <p>_com.log and the rest stay
  anything not named above      a file this module does not know is never
                                packed: Prime's own output among them

A peptide with no result file is left whole, docking and MM-GBSA folders
both: a resume docks it again from what is there. A run with a job still
queued or running is left alone. Nothing is removed until the archive has
been read back and every file found in it at its own size. Undo puts back
every file that is not there again: one a rescoring has written since is
newer than the copy in the archive, and stays.

    python -m adcp_dock.tidy RUNDIR [--dry-run] [--undo] [--force]
"""
import json
import os
import re
import sys
import tarfile
import threading
import time

from . import trajectories

ARCHIVE_DIR = "archive"
STATUS = "tidy.json"
#: The tidy job of a run, kept apart from job.json so that the check for live
#: jobs never waits on the tidy itself.
JOBFILE = "tidy_job.json"
SCRIPT = "06_tidy.sbatch"
LIVE = ("PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "REQUEUED",
        "SUSPENDED")
#: gzip's own default: most of what is packed is text, and level 9 takes about
#: twice as long for a few percent.
LEVEL = 6

#: What a pose folder packs, by what follows the folder's own name.
POSE_TAILS = (
    (".leap", "tleap's input"), ("_rec.pdb", "tleap's input"),
    ("_lig.pdb", "tleap's input"), ("_rec_forleap.pdb", "tleap's input"),
    ("_pep_forleap.pdb", "tleap's input"), ("_mirror_pep.pdb", "tleap's input"),
    ("_mirror.leap", "tleap's input"), ("_mirror.parm7", "tleap's input"),
    ("_mirror.rst7", "tleap's input"), ("_com_tleap.rst7", "tleap's input"),
    ("_com.pqr", "the protonation's own output"),
    ("_com_reduce.pdb", "the protonation's own output"),
    ("_min.in", "sander's and MMPBSA.py's working files"),
    ("_decomp.in", "sander's and MMPBSA.py's working files"),
    ("_decomp_total.dat", "sander's and MMPBSA.py's working files"),
)
POSE_NAMES = {"mdinfo": "sander's and MMPBSA.py's working files",
              "reference.frc": "sander's and MMPBSA.py's working files"}
_POSE = re.compile(r"^pose\d+$")


def _read(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def _write(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=1)
    os.replace(tmp, path)


def is_run(rundir):
    from . import job as jobmod
    return os.path.isfile(jobmod.path_of(rundir))


def count_inodes(rundir):
    """Files and folders under the run, the run itself counted."""
    n = 1
    for _root, dirs, files in os.walk(rundir):
        n += len(dirs) + len(files)
    return n


def finished(rundir):
    """({docking folder name}, {MM-GBSA folder name}) of every peptide with a
    result file: the peptides a resume will not dock again."""
    from . import job as jobmod, mmforms
    keys, names = set(), set()
    d = os.path.join(rundir, "results")
    for f in sorted(os.listdir(d)) if os.path.isdir(d) else []:
        if not f.endswith(".json"):
            continue
        seq = _read(os.path.join(d, f)).get("peptide")
        if seq:
            keys.add(jobmod.key(seq))
            names.add(mmforms.name(seq))
    return keys, names


def _file(path):
    return os.path.isfile(path) and not os.path.islink(path)


def docking_kind(name):
    """Why a file directly in a peptide's docking folder is packed, or None."""
    if name.endswith(".trg"):
        return "the copy of the grid target ADCP docked from"
    return None


def pose_kind(name, prefix, folder):
    """Why a file in an MM-GBSA pose folder is packed, or None."""
    if name in POSE_NAMES:
        return POSE_NAMES[name]
    if name.startswith("_MMPBSA_"):
        return "sander's and MMPBSA.py's working files"
    if not name.startswith(prefix):
        return None
    tail = name[len(prefix):]
    for t, why in POSE_TAILS:
        if tail == t:
            return why
    if tail == "_com.pdb":
        fl = os.path.join(folder, prefix + "_com_forleap.pdb")
        if os.path.isfile(fl) and os.path.getsize(fl) > 0:
            return "tleap's input"
    return None


def _walk_files(top):
    for root, dirs, names in os.walk(top):
        dirs.sort()
        for n in sorted(names):
            p = os.path.join(root, n)
            if _file(p):
                yield p


def survey(rundir):
    """([(file relative to the run, why it is packed)], [(folder, why it is
    left whole)]). Nothing is touched."""
    rundir = os.path.abspath(rundir)
    keys, names = finished(rundir)
    out, left = [], []
    rel = lambda p: os.path.relpath(p, rundir)               # noqa: E731
    dock = os.path.join(rundir, "docking")
    for key in sorted(os.listdir(dock)) if os.path.isdir(dock) else []:
        wd = os.path.join(dock, key)
        if not os.path.isdir(wd) or os.path.islink(wd):
            continue
        if key not in keys:
            left.append((rel(wd), "no result yet, so a resume docks it again "
                                  "from what is here"))
            continue
        for n in sorted(os.listdir(wd)):
            p = os.path.join(wd, n)
            if n == trajectories.WORKING and os.path.isdir(p) \
                    and not os.path.islink(p):
                out.extend((rel(f), "ADCP 1.1's working folder")
                           for f in _walk_files(p))
                continue
            why = docking_kind(n) if _file(p) else None
            if why:
                out.append((rel(p), why))
    mm = os.path.join(rundir, "mmgbsa")
    for name in sorted(os.listdir(mm)) if os.path.isdir(mm) else []:
        md = os.path.join(mm, name)
        if not os.path.isdir(md) or os.path.islink(md):
            continue                    # the caches and RADII, kept
        if name not in names:
            left.append((rel(md), "no result yet, so a resume scores it again "
                                  "from what is here"))
            continue
        for pose in sorted(os.listdir(md), key=lambda s: (len(s), s)):
            pd = os.path.join(md, pose)
            if not _POSE.match(pose) or not os.path.isdir(pd) \
                    or os.path.islink(pd):
                continue
            for n in sorted(os.listdir(pd)):
                p = os.path.join(pd, n)
                why = pose_kind(n, name, pd) if _file(p) else None
                if why:
                    out.append((rel(p), why))
    return out, left


def plan(rundir):
    """(files to pack, relative to the run, in a stable order; [(folder, why
    it is left whole)]). Nothing is touched."""
    got, left = survey(rundir)
    return [f for f, _w in got], left


def breakdown(rundir):
    """{why: count} of what a tidy would pack, for a dry run to report."""
    out = {}
    for _f, why in survey(rundir)[0]:
        out[why] = out.get(why, 0) + 1
    return out


def job_ids(j):
    """Every job the run recorded in job.json: prepare, dock, report, panel,
    rescore and a resume's, by id."""
    return [str(v) for _k, v in sorted((j.get("jobs") or {}).items()) if v]


def live_jobs(rundir):
    """The run's jobs still queued or running, as "stage id"; empty when the
    scheduler is not there to ask."""
    from . import job as jobmod, slurm
    try:
        j = jobmod.load(rundir)
    except (OSError, ValueError):
        return []
    live = []
    for k, v in sorted((j.get("jobs") or {}).items()):
        if not v:
            continue
        try:
            st = slurm.state(v, rundir=rundir)
        except Exception:                                       # noqa: BLE001
            return []                   # no squeue or sacct here
        if any(str(s).split()[0].upper() in LIVE for s in st):
            live.append("%s %s" % (k, v))
    return live


def _stamp(adir):
    base = time.strftime("%Y%m%d_%H%M%S")
    name, n = "tidy_%s.tar.gz" % base, 1
    while os.path.exists(os.path.join(adir, name)) \
            or os.path.exists(os.path.join(adir, name + ".part")):
        n += 1
        name = "tidy_%s_%d.tar.gz" % (base, n)
    return name


def tidy(rundir, dry_run=False, force=False):
    """Pack and remove. Nothing is removed until the archive is read back and
    found to hold every file it was given, each at its own size."""
    rundir = os.path.abspath(rundir)
    if not is_run(rundir):
        return {"error": "%s is not a run folder (no job.json)" % rundir}
    before = count_inodes(rundir)
    files, left = plan(rundir)
    base = {"inodes_before": before,
            "left": [{"folder": f, "why": w} for f, w in left]}
    if not files:
        return dict(base, packed=0, inodes_after=before, note="nothing to pack")
    if not force:
        live = live_jobs(rundir)
        if live:
            return {"error": "jobs still queued or running (%s); tidy once "
                             "the run is done" % ", ".join(live)}
    if dry_run:
        return dict(base, packed=len(files), inodes_after=before - len(files),
                    dry_run=True, files=files)
    adir = os.path.join(rundir, ARCHIVE_DIR)
    os.makedirs(adir, exist_ok=True)
    name = _stamp(adir)
    archive = os.path.join(adir, name)
    part = archive + ".part"
    sizes = {}
    try:
        with tarfile.open(part, "w:gz", compresslevel=LEVEL) as tar:
            for rel in files:
                path = os.path.join(rundir, rel)
                sizes[rel] = os.path.getsize(path)
                tar.add(path, arcname=rel, recursive=False)
        with tarfile.open(part, "r:gz") as tar:
            held = {m.name: m.size for m in tar.getmembers() if m.isfile()}
    except (OSError, tarfile.TarError) as exc:
        if os.path.exists(part):
            os.remove(part)
        return {"error": "the archive could not be written (%s); nothing "
                         "removed" % exc}
    missing = [f for f in files if held.get(f) != sizes[f]]
    if missing:
        os.remove(part)
        return {"error": "the archive missed %d file(s) (%s ...); nothing "
                         "removed" % (len(missing), missing[0])}
    os.replace(part, archive)
    # Recorded before anything is removed, so a job stopped part way through
    # the removal still leaves a run that says it is tidied and can be undone.
    record = {"archive": name, "when": time.strftime("%Y-%m-%d %H:%M:%S"),
              "packed": len(files), "bytes": os.path.getsize(archive),
              "inodes_before": before}
    sfile = os.path.join(adir, STATUS)
    hist = _read(sfile)
    hist.setdefault("archives", []).append(record)
    _write(sfile, hist)
    for rel in files:
        try:
            os.remove(os.path.join(rundir, rel))
        except FileNotFoundError:
            pass
    _prune(rundir, files)
    record["inodes_after"] = count_inodes(rundir)
    _write(sfile, hist)
    return dict(base, **record)


def _prune(rundir, files):
    """Empty folders the packing left, deepest first: only ADCP 1.1's working
    folder and what is inside it. A peptide's own folders always stay."""
    dirs = set()
    for rel in files:
        parts = rel.split(os.sep)
        if len(parts) > 3 and parts[0] == "docking" \
                and parts[2] == trajectories.WORKING:
            for i in range(3, len(parts)):
                dirs.add(os.sep.join(parts[:i]))
    for d in sorted(dirs, key=lambda x: x.count(os.sep), reverse=True):
        try:
            os.rmdir(os.path.join(rundir, d))
        except OSError:
            pass


def _archives(adir):
    if not os.path.isdir(adir):
        return []
    return sorted(f for f in os.listdir(adir)
                  if f.startswith("tidy_") and f.endswith(".tar.gz"))


def _safe(m, root):
    """A regular file whose place is inside the run, followed through any
    link on the way. Checked here rather than by tarfile's own filters,
    which would take the group's write permission off every file."""
    name = m.name
    if not (m.isfile() and name and not name.startswith(("/", ".."))
            and ".." not in name.split("/")):
        return False
    where = os.path.realpath(os.path.join(root, os.path.dirname(name)))
    return where == root or where.startswith(root + os.sep)


def undo(rundir):
    """Every tidy archive back where it came from, then the archives go.

    The newest archive first, and a file that is there again is left as it
    is: one written since the tidy, by a rescoring, is newer than the copy in
    any archive."""
    rundir = os.path.abspath(rundir)
    if not is_run(rundir):
        return {"error": "%s is not a run folder (no job.json)" % rundir}
    adir = os.path.join(rundir, ARCHIVE_DIR)
    arcs = _archives(adir)
    if not arcs:
        return {"restored": 0, "kept": 0, "note": "no tidy archive here"}
    n = kept = 0
    root = os.path.realpath(rundir)
    # Every member is checked by _safe, so tarfile is told not to filter:
    # its filters change the permissions the files were packed with.
    how = ({"filter": "fully_trusted"}
           if hasattr(tarfile, "fully_trusted_filter") else {})
    for a in reversed(arcs):
        path = os.path.join(adir, a)
        with tarfile.open(path, "r:gz") as tar:
            todo = []
            for m in tar.getmembers():
                if not _safe(m, root):
                    continue
                if os.path.lexists(os.path.join(rundir, m.name)):
                    kept += 1
                    continue
                todo.append(m)
            tar.extractall(rundir, members=todo, **how)
            n += len(todo)
        os.remove(path)
    for leftover in (STATUS,):
        if os.path.exists(os.path.join(adir, leftover)):
            os.remove(os.path.join(adir, leftover))
    try:
        os.rmdir(adir)
    except OSError:
        pass
    return {"restored": n, "kept": kept, "archives": len(arcs),
            "inodes_now": count_inodes(rundir)}


def status(rundir):
    """What the page shows: whether the run is tidied, and the counts."""
    arcs = _read(os.path.join(rundir, ARCHIVE_DIR, STATUS)).get("archives", [])
    return {"tidied": bool(arcs), "archives": len(arcs),
            "packed": sum(a.get("packed", 0) for a in arcs),
            "inodes_before": arcs[0].get("inodes_before") if arcs else None,
            "inodes_after": arcs[-1].get("inodes_after") if arcs else None,
            "last": arcs[-1].get("when") if arcs else None}


def _alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, TypeError, ValueError):
        return False


def job_state(rundir):
    """The last tidy job and whether it is still going, for the page."""
    rec = _read(os.path.join(rundir, JOBFILE))
    if not rec.get("job_id"):
        return {}
    if rec.get("local"):
        word = rec.get("local_state") or "unknown"
        if word == "RUNNING" and not _alive(rec.get("pid")):
            word = "CANCELLED"
        rec["live"] = word == "RUNNING"
        rec["state"] = word
        return rec
    try:
        from . import slurm
        st = slurm.state(rec["job_id"])
    except Exception:                                           # noqa: BLE001
        st = {}
    rec["live"] = any(str(s).split()[0].upper() in LIVE for s in st)
    rec["state"] = ", ".join(sorted(st)) or "unknown"
    return rec


def _here(rundir, undo_, rec):
    """The tidy, or its undo, on this machine: the local executor's way, where
    there is no scheduler and the page runs on the person's own computer."""
    from . import slurm
    path = os.path.join(rundir, JOBFILE)
    rec.update({"job_id": slurm.LOCAL_PREFIX + ("untidy" if undo_ else "tidy"),
                "local": True, "pid": os.getpid(), "local_state": "RUNNING"})
    _write(path, rec)

    def go():
        try:
            r = undo(rundir) if undo_ else tidy(rundir)
            word = "FAILED" if r.get("error") else "COMPLETED"
        except Exception:                                       # noqa: BLE001
            word = "FAILED"
        done = _read(path)
        done["local_state"] = word
        _write(path, done)

    threading.Thread(target=go, daemon=True,
                     name="adcp-tidy-%s" % os.path.basename(rundir)).start()


def submit(rundir, undo_=False, after=None, auto=False):
    """The tidy, or its undo, as a short job: packing thousands of files is
    work for a compute node, not the one serving the page. Recorded in
    tidy_job.json, apart from job.json, so the live job check never waits on
    the tidy itself.

    `auto` is the run's own tidy, queued when the run is: `after` holds it
    until every one of those jobs has ended, however it ended, and what there
    is to pack is worked out when it starts, not now, when nothing has run."""
    from . import job as jobmod, pipeline, slurm
    rundir = os.path.abspath(rundir)
    if not is_run(rundir):
        return {"error": "%s is not a run folder (no job.json)" % rundir}
    if not auto:
        live = live_jobs(rundir)
        if live:
            return {"error": "jobs still queued or running (%s); %s once the run "
                             "is done" % (", ".join(live),
                                          "undo the tidy" if undo_ else "tidy")}
        mine = job_state(rundir)
        if mine.get("live") and not waiting_auto(mine):
            return {"error": "a tidy of this run is already queued or running "
                             "(job %s)" % mine["job_id"]}
    if undo_:
        if not _archives(os.path.join(rundir, ARCHIVE_DIR)):
            return {"error": "this run has not been tidied, so there is "
                             "nothing to put back"}
    elif not auto and not plan(rundir)[0]:
        return {"error": "nothing to tidy: every file here is one the "
                         "results read"}
    if not auto and mine.get("live"):
        # pressed by hand while the run's own tidy waits: this one does it now
        slurm.cancel(mine["job_id"])
    stage = "untidy" if undo_ else "tidy"
    rec = {"undo": bool(undo_), "at": time.strftime("%Y-%m-%d %H:%M:%S")}
    if auto:
        rec.update(auto=True, after=after or "")
    if slurm.local_mode():
        _here(rundir, undo_, rec)
        return {"job": rec["job_id"], "undo": bool(undo_), "local": True}
    if not slurm.available():
        return {"error": "sbatch was not found. This tool submits its work to "
                         "SLURM rather than running it where you are sitting."}
    j = jobmod.load(rundir)
    common = pipeline._scripts(rundir, j)
    common["mem"] = "2G"
    script = slurm.write_script(
        os.path.join(rundir, SCRIPT), name=pipeline.job_name(rundir, j, stage),
        cpus=1, time=slurm.within_limit("00:30:00", j.get("partition"),
                                        j.get("sbatch_extra")),
        stage=stage,
        notes=["ADCP docking: %s" % ("put the tidied files back" if undo_ else
                                     "pack the files the results no longer "
                                     "read")], **common)
    rec["job_id"] = slurm.submit(script, depends_on=after, kind="afterany")
    _write(os.path.join(rundir, JOBFILE), rec)
    return {"job": rec["job_id"], "undo": bool(undo_), "auto": bool(auto),
            "after": after or ""}


def waiting_auto(state):
    """Whether a tidy job is the run's own, queued with it and still waiting:
    one that has not started holds nothing, and is queued again behind
    whatever is added to the run (arm)."""
    return bool(state.get("auto")) and "PENDING" in str(state.get("state") or "") \
        and "RUNNING" not in str(state.get("state") or "")


def arm(rundir, also=()):
    """Queue the run's own tidy behind every job it has queued, the ones just
    submitted (`also`) among them, as Tidy would be pressed by hand once they
    had all ended. An automatic tidy still waiting is canceled first: the run
    has more jobs than it was waiting for. Keep every file (keep_files) leaves
    the run alone, and so does the local executor, whose stages run one after
    another on this machine, where Tidy is pressed by hand."""
    from . import job as jobmod, slurm
    rundir = os.path.abspath(rundir)
    j = jobmod.load(rundir)
    if j.get("keep_files"):
        return {"armed": False, "why": "every file of this run is kept, as the form asked"}
    if slurm.local_mode():
        return {"armed": False, "why": "runs here are tidied by hand, with Tidy"}
    old = job_state(rundir)
    if old.get("live") and waiting_auto(old):
        slurm.cancel(old["job_id"])
    waits = []
    for jid in [str(x) for x in also if x] + live_jobs(rundir):
        if jid not in waits:
            waits.append(jid)
    got = submit(rundir, after=":".join(waits) or None, auto=True)
    got["armed"] = "error" not in got
    return got


def main(argv=None):
    args = list(sys.argv[1:] if argv is None else argv)
    flags = {a for a in args if a.startswith("--")}
    dirs = [a for a in args if not a.startswith("--")]
    if not dirs or flags - {"--dry-run", "--undo", "--force"}:
        print(__doc__)
        return 2
    rc = 0
    for d in dirs:
        if "--undo" in flags:
            r = undo(d)
            if r.get("error"):
                print("tidy undo %s: %s" % (d, r["error"]))
                rc = 1
                continue
            print("tidy undo %s: %s" % (d, r.get("note") or (
                "%d file(s) restored from %d archive(s), %d inodes now%s"
                % (r["restored"], r["archives"], r["inodes_now"],
                   "; %d written again since the tidy were kept as they are"
                   % r["kept"] if r["kept"] else ""))))
            continue
        r = tidy(d, dry_run="--dry-run" in flags, force="--force" in flags)
        if r.get("error"):
            print("tidy %s: %s" % (d, r["error"]))
            rc = 1
            continue
        for p in r.get("left") or []:
            print("tidy %s: %s left whole: %s" % (d, p["folder"], p["why"]))
        if r.get("dry_run"):
            print("tidy %s (dry run): %d file(s) would be packed, %d -> %d "
                  "inodes" % (d, r["packed"], r["inodes_before"],
                              r["inodes_after"]))
            for why, n in sorted(breakdown(d).items(), key=lambda kv: -kv[1]):
                print("    %6d  %s" % (n, why))
        else:
            print("tidy %s: %d file(s) packed into %s; %d -> %d inodes%s"
                  % (d, r["packed"], r.get("archive") or "nothing",
                     r["inodes_before"], r["inodes_after"],
                     " (%s)" % r["note"] if r.get("note") else ""))
    return rc


if __name__ == "__main__":
    sys.exit(main())
