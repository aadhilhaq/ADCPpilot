"""What one SLURM task runs. Not meant to be called by hand.

By Aadhil Haq

    python -m adcp_pilot.worker /path/to/job.json prepare
    python -m adcp_pilot.worker /path/to/job.json 7        (peptide index)
    python -m adcp_pilot.worker /path/to/job.json report
    python -m adcp_pilot.worker /path/to/job.json panel   (an older run)
    python -m adcp_pilot.worker /path/to/job.json strain  (an older run)
    python -m adcp_pilot.worker /path/to/job.json panel+strain  (both, in turn)
    python -m adcp_pilot.worker /path/to/job.json rescore:openmm
    python -m adcp_pilot.worker /path/to/job.json tidy    (a finished run)
    python -m adcp_pilot.worker /path/to/job.json untidy  (undo its tidy)
"""
import os
import sys

from . import job as jobmod, pipeline


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) < 2:
        print(__doc__)
        return 2
    jobfile, stage = argv[0], argv[1]
    rundir = os.path.dirname(os.path.abspath(jobfile))
    if stage in ("tidy", "untidy"):
        # The run's tidy, or its undo, in a short job of its own. Its exit
        # code is the job's, so a tidy that stopped shows as failed.
        from . import tidy
        return tidy.main([rundir] + (["--undo"] if stage == "untidy" else []))
    j = jobmod.load(jobfile)

    def log(*a):
        print(*a, flush=True)

    # The measures a finished run can have added in one short job, joined
    # with "+" when there are several (pipeline.submit_measures), run in turn.
    light = {"panel": pipeline.stage_panel, "strain": pipeline.stage_strain}
    parts = stage.split("+")
    if stage == "prepare":
        pipeline.stage_prepare(rundir, j, log=log)
    elif stage == "report":
        pipeline.stage_report(rundir, j, log=log)
    elif all(p in light for p in parts):
        # Each on its own: one that stops leaves the others to do theirs,
        # and the job is failed afterwards so the stop is seen.
        bad = 0
        for p in parts:
            try:
                light[p](rundir, j, log=log)
            except (Exception, SystemExit) as exc:              # noqa: BLE001
                import traceback
                traceback.print_exc()
                log("%s stopped: %s" % (p, exc))
                bad += 1
        return 1 if bad else 0
    elif stage.startswith("rescore:"):
        pipeline.stage_rescore(rundir, stage.split(":", 1)[1], j, log=log)
    else:
        pipeline.stage_dock(rundir, int(stage), j, log=log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
