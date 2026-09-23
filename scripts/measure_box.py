#!/usr/bin/env python3
"""Re-derive sizing.BOX_PER_RESIDUE and BOX_INTERCEPT from the structures.

The constants in sizing.py are a least-squares fit over the bound peptides
in the HOPE peptide pipeline's benchmark set (hope_core.benchmarks). Point
HOPE_PIPELINES at a checkout of that pipeline, run this, and compare.

    HOPE_PIPELINES=/path/to/HOPE-pipelines python scripts/measure_box.py [dir]
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if os.environ.get("HOPE_PIPELINES"):
    sys.path.insert(0, os.environ["HOPE_PIPELINES"])


def main(argv):
    from hope_core import benchmarks as B
    from adcp_dock import sizing
    d = (argv[0] if argv else
         os.path.join(os.environ.get("SCRATCH", os.path.expanduser("~")),
                      "hope_benchmarks"))
    rows = []
    print("%-6s %-10s %4s %8s" % ("PDB", "bracket", "len", "cube A"))
    for pid, c in sorted(B.BENCHMARKS.items()):
        try:
            model = B.model_of(pid, d)
            res = B.residues(model, c.peptide_chain)
            xyz = B.heavy(res)
        except Exception as exc:                               # noqa: BLE001
            print("%-6s  skipped (%s)" % (pid, str(exc)[:50]))
            continue
        if not len(xyz):
            continue
        cube = float((xyz.max(axis=0) - xyz.min(axis=0)).max())
        rows.append((len(res), cube))
        print("%-6s %-10s %4d %8.1f" % (pid, c.bracket, len(res), cube))
    if len(rows) < 3:
        print("\nnot enough structures here to fit anything. Fetch them "
              "first:  hope-bench fetch")
        return 1
    a = np.array(rows, float)
    m, b = np.linalg.lstsq(
        np.vstack([a[:, 0], np.ones(len(a))]).T, a[:, 1], rcond=None)[0]
    resid = a[:, 1] - (m * a[:, 0] + b)
    rms = float(np.sqrt((resid ** 2).mean()))
    print("\n  %d peptides, %d-%d residues" % (len(a), a[:, 0].min(),
                                               a[:, 0].max()))
    print("  measured   edge = %.2f N + %.2f, residual RMS %.1f A" % (m, b, rms))
    print("  sizing.py  edge = %.2f N + %.2f, residual RMS %.1f A"
          % (sizing.BOX_PER_RESIDUE, sizing.BOX_INTERCEPT,
             sizing.BOX_RESIDUAL_RMS))
    off = (abs(m - sizing.BOX_PER_RESIDUE) > 0.05
           or abs(b - sizing.BOX_INTERCEPT) > 0.5)
    print("\n  %s" % ("THESE DISAGREE: sizing.py needs updating, or this "
                      "ran against a different set." if off
                      else "they agree."))
    return 1 if off else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
