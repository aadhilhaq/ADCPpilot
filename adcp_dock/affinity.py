"""`adcp-dock affinity`: what a target affinity corresponds to, in the units
this tool reports.

By Aadhil Haq

This prints a scale rather than a Kd column in the results table, because
an MM-GBSA score cannot be converted into a Kd:

  * An MM-GBSA interaction energy has no configurational entropy term. It
    is an interaction enthalpy plus an implicit-solvent term, so it is not
    the quantity in dG = RT ln(Kd).
  * Its magnitude follows the size of the interface and the receptor more
    than the affinity, so the same peptide set can score very differently
    on two receptors.
  * Converted anyway, typical scores give dissociation constants far below
    anything physically possible.

So the scale is printed and the reader makes the comparison. dg_for_kd is
in _vendor/mmgbsa.py; there is deliberately no kd_for_dg, and a test checks
that it stays absent.
"""
from __future__ import annotations

from ._vendor import mmgbsa


def table(temperature=298.15, spread=None, out=print):
    out("dG = RT ln(Kd)    RT = %.4f kcal/mol at %.2f K\n"
        % (0.0019872041 * temperature, temperature))
    out("  %-9s %12s" % ("Kd", "true dG"))
    for lab, _kd, dg in mmgbsa.affinity_table(temperature):
        out("  %-9s %9.2f kcal/mol" % (lab, dg))
    window = (mmgbsa.dg_for_kd(1e-9, temperature)
              - mmgbsa.dg_for_kd(1e-12, temperature))
    out("\n  The whole picomolar-to-nanomolar range is %.1f kcal/mol."
        % window)
    if spread:
        out("  The pose spread on one peptide in this run is %.1f." % spread)
        if spread > window:
            out("  The noise is wider than the window: an MM-GBSA number "
                "cannot\n  place a peptide inside it.")
        else:
            out("  That is inside the window, which is unusual; check "
                "mmgbsa_n\n  before reading anything into it.")
    else:
        out("  The spread between poses of one peptide is usually several")
        out("  times wider than this window, which is why the results table")
        out("  reports the spread beside the score and never a Kd. Give a")
        out("  run directory to compare against that run's own spread.")
    out("")


def spread_of(rundir):
    """The widest pose-to-pose spread in a finished run, or None.

    Read from the run's own results rather than typed, so `--spread` with no
    number answers the question about the run in front of you.
    """
    import glob
    import json
    import os
    worst = None
    for f in glob.glob(os.path.join(rundir, "results", "*.json")):
        try:
            vals = [v for v in (json.load(open(f)).get("mmgbsa") or [])
                    if v is not None]
        except (OSError, ValueError):
            continue
        if len(vals) > 1:
            s = max(vals) - min(vals)
            worst = s if worst is None else max(worst, s)
    return worst
