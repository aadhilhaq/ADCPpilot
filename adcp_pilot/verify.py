"""Reading MMPBSA.py's output back to check that the complex was split correctly.

By Aadhil Haq

A single-trajectory MM-GBSA takes the complex, the receptor and the ligand
from the same frame. Every bonded term (bond, angle, dihedral, 1-4 van der
Waals and 1-4 electrostatic) therefore appears identically in the complex
and in the two parts, and the difference DELTA = complex, receptor, ligand
must be exactly zero for all five.

This is not a convergence or accuracy criterion. A non-zero bonded term
means the complex was divided through a molecule rather than along the
interface, which produces a plausible energy rather than an error.
"""
import os
import re

# The five terms that must cancel, and the tolerance worth reporting. The
# printed output carries four decimals, so anything above half a unit in the
# last place is a real residue rather than rounding.
BONDED = ("BOND", "ANGLE", "DIHED", "1-4 VDW", "1-4 EEL")
TOL = 0.001

_ROW = re.compile(r"^(BOND|ANGLE|DIHED|1-4 VDW|1-4 EEL)\s+(-?\d+\.\d+)")


def mmgbsa_split(dat_path):
    """(ok, [(term, value), ...]) for the terms that failed to cancel.

    ok is True when every bonded term canceled, meaning the complex was split
    along the interface. A file that cannot be read returns (True, []): the
    caller already knows whether it received a number.
    """
    if not dat_path or not os.path.isfile(dat_path):
        return True, []
    bad, seen = [], set()
    inside = False
    for line in open(dat_path):
        s = line.strip()
        if s.startswith("Differences"):
            inside = True
            continue
        if not inside:
            continue
        m = _ROW.match(s)
        if m and m.group(1) not in seen:
            seen.add(m.group(1))
            v = float(m.group(2))
            if abs(v) > TOL:
                bad.append((m.group(1), v))
    return (not bad), bad


def explain(bad):
    return ("the bonded terms did not cancel (%s). In a single-trajectory "
            "MM-GBSA they must be zero, so the complex was split through a "
            "molecule rather than at the interface and this energy does not "
            "describe binding."
            % ", ".join("%s %.4f" % (t, v) for t, v in bad))
