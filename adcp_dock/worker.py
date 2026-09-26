"""What one SLURM task runs. Not meant to be called by hand.

By Aadhil Haq

    python -m adcp_dock.worker /path/to/job.json prepare
    python -m adcp_dock.worker /path/to/job.json 7        (peptide index)
    python -m adcp_dock.worker /path/to/job.json report
    python -m adcp_dock.worker /path/to/job.json panel   (an older run)
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
    j = jobmod.load(jobfile)

    def log(*a):
        print(*a, flush=True)

    if stage == "prepare":
        pipeline.stage_prepare(rundir, j, log=log)
    elif stage == "report":
        pipeline.stage_report(rundir, j, log=log)
    elif stage == "panel":
        pipeline.stage_panel(rundir, j, log=log)
    else:
        pipeline.stage_dock(rundir, int(stage), j, log=log)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
