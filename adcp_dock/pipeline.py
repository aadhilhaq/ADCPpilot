"""The computational stages, and how a run is assembled from them.

    prepare   receptor -> prepared PDB -> PDBQT -> grid maps (.trg)
    dock      one peptide per array task, then MM-GBSA on its kept poses
    report    gather the scores into a table

The preparation, docking and rescoring code is vendored in _vendor; this
module supplies the order, the checks between stages and the job scripts.
"""
import os
import re
import types

from . import job as jobmod, report, slurm, target, verify
from . import tools as toolmod
from ._vendor import adcp, mmgbsa, poseviewer, prep, receptor


def tools_for(j):
    return toolmod.discover(roots=j.get("tools_roots") or [])


def _reccfg(j, prepared=None):
    """The fields receptor.load and receptor.box read, taken from the job file.

    A namespace rather than a dict, because the vendored modules read attributes
    while this tool stores its settings as JSON.
    """
    return types.SimpleNamespace(
        receptor=prepared or j["receptor"],
        chain=",".join(j.get("chains") or []),
        # The hotspot string is passed exactly as entered: receptor's
        # own parser reads both "A:30,A:31" and "98-112", and parsing
        # it here as well would let the two spellings diverge.
        hotspots=j.get("hotspots") or "",
        shell_radius=float(j.get("shell_radius", 10.0)),
        box_padding=float(j.get("box_padding", 8.0)),
        box_min=float(j.get("box_min", 20.0)),
        box_max=float(j.get("box_max", 30.0)),
        box_chains=j.get("box_chains", "ignore"),
        box_chains_min_atoms=int(j.get("box_chains_min_atoms", 20)),
        forcefield=j.get("forcefield", "AMBER"))


def stage_prepare(rundir, j=None, log=print):
    """Receptor preparation and the grid maps."""
    j = j or jobmod.load(rundir)
    T = tools_for(j)
    recdir = os.path.join(rundir, "receptor")
    cfg = _reccfg(j)

    # Load the receptor and build the box together: load_for_box reads
    # the chains, refuses a hotspot that is not in the structure, builds
    # the box from the hotspots, and reports any other chain the box
    # reaches.
    log("reading the receptor and sizing the box")
    rec, centre, size = receptor.load_for_box(cfg, log=log)
    log("  site: %d hotspot residue(s) on chain(s) %s, box %.1f x %.1f x %.1f A"
        % (len(rec.hotspot_keys), ", ".join(rec.chains),
           size[0], size[1], size[2]))
    # Reported where the box is decided rather than only in the form's
    # estimate: a command-line run, a resumed run and a reused job.json
    # all arrive here without passing through the page.
    from . import sizing as _sz
    _longest = max((len(p.get("sequence") or "")
                    for p in (j.get("peptides") or [])), default=0)
    _w = _sz.box_warning(_longest, size, j.get("box_max"))
    if _w:
        log("  BOX: %s" % _w)

    # A site larger than box_max is clipped silently, and every pose is
    # then scored against a pocket that was not requested.
    outside = receptor.hotspots_outside_box(rec, centre, size)
    if outside:
        log("  WARNING the box does not contain %d of the residues it was "
            "built from: %s"
            % (len(outside), ", ".join("%s by %.1f A" % (a, b)
                                       for a, b in outside[:6])))
        log("    no pose can reach them. Raise box_max, or split the site "
            "into two runs.")

    # The full preparation: alternate locations, modified residues,
    # waters, cofactors, missing side chains, termini, disulphides,
    # protonation at the chosen pH, then a restrained minimisation.
    log("preparing the receptor (this is the full Amber preparation)")
    prepared = prep.run(
        j["receptor"], recdir, ph=float(j.get("ph", 7.4)),
        forcefield=j.get("forcefield", "AMBER"), log=log,
        chain=",".join(rec.chains),
        keep_hetatms=bool(j.get("keep_hetatms", False)),
        # Ions go with the other heteroatoms unless those are kept, as the
        # page says: a kept ion has no template in the minimisation's force
        # field, and the minimisation is then skipped.
        keep_ions=bool(j.get("keep_hetatms", False)),
        fill_loops=bool(j.get("fill_loops", False)),
        minimise_structure=bool(j.get("minimise", True)),
        restraint=float(j.get("restraint", 10.0)))

    # PDBQT from the prepared structure, then the grid maps over the box.
    # The prepared file is passed in beside the loaded structure rather
    # than reloaded, so the chains and hotspots remain the validated
    # ones.
    # prepare_receptor falls back to the chain as read when it cannot use
    # the prepared structure, and says so only in this stage's log. The
    # fallback is noted in preparation.log too, which the results page
    # reads, so a receptor docked without its preparation is not presented
    # as a prepared one.
    fell_back = []

    def _log(*a):
        if FALLBACK_SIGN in " ".join(str(x) for x in a):
            fell_back.append(True)
        log(*a)

    recq = receptor.prepare_pdbqt(rec, recdir, T, prepared_pdb=prepared,
                                  log=_log)
    if fell_back:
        with open(os.path.join(recdir, "preparation.log"), "a") as fh:
            fh.write("  receptor pdbqt: skipped; prepare_receptor could not "
                     "use the prepared structure, so the docking receptor was "
                     "written from the chain as read, with hydrogens added "
                     "geometrically\n")
    # target.build is the vendored builder plus a check of what it
    # produced: a truncated or mismatched .trg still docks, into empty
    # maps.
    trg = target.build(recq, centre, size, os.path.join(rundir, "target"), T,
                       name="receptor", log=log)

    j["site"] = {"centre": [float(v) for v in centre],
                 "size": [float(v) for v in size],
                 "hotspots": [[c, int(n)] for c, n in rec.hotspot_keys],
                 "chains": list(rec.chains),
                 "outside_box": outside}
    j["receptor_info"] = {"prepared_pdb": prepared, "pdbqt": recq}
    j["target"] = trg
    j["state"] = "prepared"
    jobmod.save(rundir, j)
    return trg


def _rec_for_mmgbsa(j, log=print):
    """The receptor MM-GBSA scores against: the prepared one where there is one, the input structure otherwise.

    The pose was docked into maps built from the prepared structure, so the
    complex must be assembled against that same structure.
    """
    prepared = (j.get("receptor_info") or {}).get("prepared_pdb")
    if prepared and os.path.isfile(prepared):
        try:
            rec = receptor.load(_reccfg(j, prepared))
            log("  MM-GBSA uses the prepared receptor")
            return rec
        except Exception as exc:                      # noqa: BLE001
            log("  could not load the prepared receptor (%s); using the "
                "input structure" % exc)
    return receptor.load(_reccfg(j))


# What receptor.prepare_pdbqt logs when it falls back to the unprepared
# chain (tests/test_api.py checks the vendored code still says it).
FALLBACK_SIGN = "Falling back to the chain as it was read"


def stage_dock(rundir, index, j=None, log=print):
    """One peptide, and its MM-GBSA if that was asked for."""
    j = j or jobmod.load(rundir)
    T = tools_for(j)
    peps = jobmod.peptide_list(j)
    if index >= len(peps):
        log("nothing to do for task %d" % index)
        return None
    name, seq = peps[index]
    trg = j.get("target") or os.path.join(rundir, "target", "receptor.trg")
    if not os.path.isfile(trg):
        raise SystemExit("no grid target at %s; run the prepare stage" % trg)
    why = target.check(trg)
    if why:
        raise SystemExit("the grid target is unusable (%s). Docking into it "
                         "would produce scores against empty maps." % why)

    wd = os.path.join(rundir, "docking")
    cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", j.get("cpus", 0)) or 0)
    log("docking %s (%s), %d replicas x %d steps on %d core(s)"
        % (name, seq, int(j.get("replicas", 50)), int(j.get("steps", 2500000)),
           cpus or 1))
    modes, pdir = adcp.dock(seq, trg, wd, T,
                            replicas=int(j.get("replicas", 50)),
                            steps=int(j.get("steps", 2_500_000)),
                            cpus=cpus, log=log)
    keep = int(j.get("poses", 10))
    out = {"name": name, "peptide": seq,
           "modes": [dict(rank=m.get("rank"), score=m.get("affinity"))
                     for m in (modes or [])][:keep]}
    log("  %d ranked pose(s)" % len(out["modes"]))

    if j.get("mmgbsa") and out["modes"]:
        rec = _rec_for_mmgbsa(j, log)
        n = min(int(j.get("mmgbsa_poses", 5)), keep)
        # Docked coordinates laid onto a clean topology, matched by
        # (residue number, atom name). ADCP writes a residue's atoms
        # in two blocks, so a pose copied in file order has
        # duplicated residues; transplant refuses anything it cannot
        # match atom for atom.
        mols, ranks = [], []
        for k in range(1, n + 1):
            pose = adcp.read_pose(pdir, k)
            if pose is None:
                break
            # At the run's pH, as the docking was. The complex is
            # re-protonated by pdb2pqr before the energy in
            # either case; this is passed so that anything
            # written from this molecule carries the right
            # charges.
            mol = poseviewer.transplant(seq, pose,
                                        ph=float(j.get("ph", 7.4)))
            if mol is None:
                log("    pose %d could not be matched atom for atom, skipped"
                    % k)
                continue
            mols.append(mol)
            ranks.append(k)
        if mols:
            mmdir = os.path.join(rundir, "mmgbsa")
            detail = {}
            # The batch entry point rather than a loop: it caches
            # per pose, so a resumed task recomputes only what is
            # missing, and it runs several poses at once on the
            # cores the task already holds.
            workers = max(1, min(cpus or 1, len(mols),
                                 int(j.get("mm_workers", 8))))
            mmgbsa.run([{"sequence": seq}], rec, {seq: mols}, mmdir, log=log,
                       prep_ph=float(j.get("ph", 7.4)), n_poses=len(mols),
                       workers=workers,
                       aggregate=j.get("mmgbsa_aggregate", "mean"),
                       details=detail)
            d = detail.get(seq) or {}
            vals = list(d.get("values") or [])
            got = list(d.get("ranks") or [])
            # Read MMPBSA.py's output back. A complex split
            # through a molecule rather than at the interface
            # still produces a number that looks like any other.
            # Kept by position, not by value: two poses can share an energy.
            good, kept = [], []
            for i, r in enumerate(got):
                dat = os.path.join(mmdir, seq, "pose%d" % r,
                                   "%s_mmgbsa.dat" % seq)
                ok, bad = verify.mmgbsa_split(dat)
                if ok:
                    good.append(vals[i])
                    kept.append(r)
                else:
                    log("    pose %d discarded: %s" % (r, verify.explain(bad)))
            out["mmgbsa"] = good
            out["mmgbsa_ranks"] = kept
        else:
            log("  no pose could be laid on a clean topology; no MM-GBSA")
            out["mmgbsa"] = []
    _save_result(rundir, seq, out)
    return out


def _save_result(rundir, seq, out):
    import json
    d = os.path.join(rundir, "results")
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "%s.json" % seq), "w") as fh:
        json.dump(out, fh, indent=1)


def stage_report(rundir, j=None, log=print):
    import glob
    import json
    j = j or jobmod.load(rundir)
    rows = []
    # The name is taken from the job file when the result has none: a run
    # docked before the dock stage wrote names has results without them,
    # and the table would otherwise show a column of dashes for peptides
    # that were named.
    by_seq = {}
    for _n, _s in jobmod.peptide_list(j):
        by_seq.setdefault(_s, _n)
    for f in sorted(glob.glob(os.path.join(rundir, "results", "*.json"))):
        r = json.load(open(f))
        rows.append(report.summarise(r["peptide"], r.get("modes"),
                                     r.get("mmgbsa"),
                                     top_n=int(j.get("top_n", 5)),
                                     name=(r.get("name")
                                           or by_seq.get(r["peptide"], ""))))
    if not rows:
        log("no results to report yet")
        return None
    # The citation list is written into the run, since a method section
    # is usually written later from the output directory. It names only
    # what this run used, so a run without MM-GBSA is not given an Amber
    # citation.
    try:
        from . import cite as citemod
        # Three formats: .bib for LaTeX, .ris for EndNote, Zotero and
        # Mendeley, .txt to read. The names come from cite.FILENAMES,
        # which the page's download buttons and bundle.RESULTS also
        # use.
        for _fmt, _name in sorted(citemod.FILENAMES.items()):
            with open(os.path.join(rundir, _name), "w") as fh:
                fh.write(citemod.formatted(_fmt, j))
        log("  wrote HOW_TO_CITE.txt, CITATION.bib and CITATION.ris; "
            "what this run used, with DOIs")
    except Exception as exc:                                    # noqa: BLE001
        log("  could not write the citation list: %s" % exc)
    tsv = report.write(rows, rundir, top_n=int(j.get("top_n", 5)))
    log(report.table(rows, top_n=int(j.get("top_n", 5))))
    j["state"] = "done"
    jobmod.save(rundir, j)
    return tsv


def outstanding(rundir, j=None):
    """[(index, name, sequence)] for the peptides with no result yet."""
    j = j or jobmod.load(rundir)
    out = []
    for i, (name, seq) in enumerate(jobmod.peptide_list(j)):
        if not os.path.isfile(os.path.join(rundir, "results", "%s.json" % seq)):
            out.append((i, name, seq))
    return out


def job_name(rundir, j, stage):
    """A SLURM job name that identifies the run.

    Generic names tell a user with several runs in the queue nothing about
    which is which, so the run's own label is used.
    """
    name = (j.get("label") or "").strip() or os.path.basename(
        os.path.normpath(rundir))
    name = re.sub(r"[^A-Za-z0-9_.-]+", "-", name).strip("-_.") or "adcp"
    # Trim after the cut, or a boundary that lands on a separator leaves
    # a trailing dash.
    return "%s-%s" % (name[:18].rstrip("-_.") or "adcp", stage)


def _scripts(rundir, j, python=None):
    import sys
    python = python or sys.executable
    logdir = os.path.join(rundir, "logs")
    os.makedirs(logdir, exist_ok=True)
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return dict(partition=j.get("partition", "cpu"), workdir=rundir,
                logdir=logdir, python=python,
                jobfile=jobmod.path_of(rundir),
                env=("export PYTHONPATH=%s:${PYTHONPATH}" % slurm.stage_pythonpath()
                     if slurm.stage_pythonpath() else ""),
                extra=slurm.site_options(j.get("account"),
                                         j.get("sbatch_extra"))[0],
                mem=j.get("mem", "32G"))


def _fit_to_machine(j):
    """Cores requested, reduced to the cores this machine has.

    A job file written for a 48-core node would otherwise ask ADCP for 48
    threads on a two-core machine, which is slower rather than faster.
    """
    n = slurm.usable_cores()
    for k in ("cpus", "mm_cpus", "mm_workers"):
        if j.get(k):
            j[k] = max(1, min(int(j[k]), n))
    return n


def _submit_local(rundir, j, stages, python, log):
    """Run the stages here (slurm.local_mode()), recorded in job.json under the same keys a SLURM submission uses.
    """
    import sys as _sys
    _fit_to_machine(j)
    ids = {k: slurm.LOCAL_PREFIX + k for k, _t in stages}
    # Saved before the thread starts: the worker reads job.json as it
    # begins, and must see the core count this machine can actually
    # provide.
    j.setdefault("jobs", {}).update(ids)
    j["state"] = "submitted"
    jobmod.save(rundir, j)
    slurm.submit_local(rundir, jobmod.path_of(rundir), stages,
                       cpus=int(j.get("cpus") or 1),
                       python=python or _sys.executable)
    for k, v in ids.items():
        log("%s: %s (running on this machine)" % (k, v))
    return ids


def submit(rundir, j=None, python=None, log=print):
    """prepare -> dock array -> report, each waiting on the one before."""
    j = j or jobmod.load(rundir)
    problems = jobmod.validate(j)
    if problems:
        raise SystemExit("this run is not ready:\n  " + "\n  ".join(problems))
    if slurm.local_mode():
        n = len(jobmod.peptide_list(j))
        return _submit_local(rundir, j, [("prepare", ["prepare"]),
                                         ("dock", [str(i) for i in range(n)]),
                                         ("report", ["report"])],
                             python, log)
    if not slurm.available():
        raise SystemExit("sbatch was not found. This tool submits its work to "
                         "SLURM rather than running it where you are sitting.")
    common = _scripts(rundir, j, python)
    s1 = slurm.write_script(os.path.join(rundir, "01_prepare.sbatch"),
                            name=job_name(rundir, j, "prep"), cpus=min(16, int(j.get("cpus", 48))),
                            time=j.get("prep_time", "04:00:00"),
                            stage="prepare", **common)
    id1 = slurm.submit(s1)
    log("prepare: %s" % id1)

    n = len(jobmod.peptide_list(j))
    thr = int(j.get("throttle") or 0)
    s2 = slurm.write_script(
        os.path.join(rundir, "02_dock.sbatch"),
        # The job name reflects what the task contains: the docking
        # and the rescoring are one task, so naming it for the
        # docking alone invites comparison with an estimate that
        # covers only half of it.
        name=job_name(rundir, j, "dock+mm" if j.get("mmgbsa") else "dock"),
        cpus=int(j.get("cpus", 48)), time=j.get("time", "08:00:00"),
        indices=slurm.index_spec(range(n)),
        throttle=("%%%d" % thr) if thr else "", **common)
    id2 = slurm.submit(s2, depends_on=id1)
    log("dock: %s (%d task%s)" % (id2, n, "" if n == 1 else "s"))

    s3 = slurm.write_script(os.path.join(rundir, "03_report.sbatch"),
                            name=job_name(rundir, j, "report"), cpus=2, time="01:00:00",
                            stage="report", **common)
    id3 = slurm.submit(s3, depends_on=id2, kind="afterany")
    log("report: %s" % id3)
    j["jobs"] = {"prepare": id1, "dock": id2, "report": id3}
    j["state"] = "submitted"
    jobmod.save(rundir, j)
    return j["jobs"]


def resume(rundir, j=None, python=None, log=print):
    """Submit only the peptides that have no result, then report again."""
    j = j or jobmod.load(rundir)
    todo = outstanding(rundir, j)
    if not todo:
        log("nothing outstanding: every peptide has a result")
        return {"todo": 0}
    trg = j.get("target")
    need_prep = not (trg and os.path.isfile(trg))
    if slurm.local_mode():
        stages = ([("prepare", ["prepare"])] if need_prep else []) + [
            ("dock", [str(i) for i, _n, _s in todo]),
            ("report", ["report"])]
        ids = _submit_local(rundir, j, stages, python, log)
        return {"todo": len(todo), "jobs": ids,
                "indices": ",".join(str(i) for i, _n, _s in todo),
                "prepared_again": need_prep}
    if not slurm.available():
        raise SystemExit("sbatch was not found")
    common = _scripts(rundir, j, python)
    dep = None
    if need_prep:
        log("the grid target is missing, so preparation runs again first")
        s1 = slurm.write_script(os.path.join(rundir, "01_prepare.sbatch"),
                                name=job_name(rundir, j, "prep"),
                                cpus=min(16, int(j.get("cpus", 48))),
                                time=j.get("prep_time", "04:00:00"),
                                stage="prepare", **common)
        dep = slurm.submit(s1)
        log("prepare: %s" % dep)
    thr = int(j.get("throttle") or 0)
    spec = slurm.index_spec(i for i, _n, _s in todo)
    s2 = slurm.write_script(os.path.join(rundir, "02_dock_resume.sbatch"),
                            name=job_name(rundir, j, "dock"), cpus=int(j.get("cpus", 48)),
                            time=j.get("time", "08:00:00"), indices=spec,
                            throttle=("%%%d" % thr) if thr else "", **common)
    id2 = slurm.submit(s2, depends_on=dep) if dep else slurm.submit(s2)
    log("dock: %s for %d peptide(s), array %s" % (id2, len(todo), spec))
    s3 = slurm.write_script(os.path.join(rundir, "03_report.sbatch"),
                            name=job_name(rundir, j, "report"), cpus=2, time="01:00:00",
                            stage="report", **common)
    id3 = slurm.submit(s3, depends_on=id2, kind="afterany")
    log("report: %s" % id3)
    j.setdefault("jobs", {})
    j["jobs"].update({"dock_resume": id2, "report": id3})
    if dep:
        j["jobs"]["prepare"] = dep
    j["state"] = "resubmitted"
    jobmod.save(rundir, j)
    return {"todo": len(todo), "jobs": j["jobs"], "indices": spec,
            "prepared_again": need_prep}
