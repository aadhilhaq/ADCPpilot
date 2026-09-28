"""The MM-GBSA engines, and what each one's number actually means.

Three programs can rescore a pose, and they do not compute the same thing.
Each carries its own force field and its own implicit solvent model, so the
numbers are on different scales and an absolute value from one cannot be put
beside an absolute value from another:

    AmberTools   ff14SB + OBC(II), igb=5, mbondi2 radii
    OpenMM       ff14SB + OBC(II), the same model on the same topology
    Schrodinger  OPLS_2005 + VSGB 2.0, a different force field entirely

AmberTools and OpenMM at least share a force field and a topology, so the
difference between them is the solver and the minimisation. Schrodinger
shares neither. On this project's own poses, AmberTools and OpenMM agreed on
the ORDER of five poses exactly while differing by about 32 kcal/mol on the
value, which is the shape of the problem: these engines are comparable by
rank and not by number.

So each engine writes its own column, none of them is merged into another,
and the table says which force field and which solvent model produced each
one. A run may carry all three.

WHAT SWAPPING THE ENGINE DOES NOT FIX. On 1BE9, five docked poses of one
peptide span 21.7 kcal/mol, against a 37.6 kcal/mol range across all twenty
constructs: the noise within a construct is 58 per cent of the signal
between constructs. That is a property of the poses, not of the rescoring,
and every engine here inherits it. More poses would help; a different solver
does not.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess


class Engine:
    """One way to compute an MM-GBSA energy, described honestly.

    `column` is where its numbers land, so two engines never overwrite each
    other. `comparable` names the engines whose ABSOLUTE values can sit in
    the same sentence as this one's; rank is comparable across all of them.
    """

    def __init__(self, name, label, forcefield, solvation, radii,
                 column, comparable, note=""):
        self.name = name
        self.label = label
        self.forcefield = forcefield
        self.solvation = solvation
        self.radii = radii
        self.column = column
        self.comparable = tuple(comparable)
        self.note = note

    def as_dict(self, available=None):
        return {"name": self.name, "label": self.label,
                "forcefield": self.forcefield, "solvation": self.solvation,
                "radii": self.radii, "column": self.column,
                "comparable": list(self.comparable), "note": self.note,
                "available": self.available() if available is None
                             else bool(available)}

    def available(self):                                    # pragma: no cover
        return False

    def describe(self):
        """One line a results table can show under the column heading."""
        bits = [self.forcefield, self.solvation]
        if self.radii and self.radii != "n/a":
            bits.append(self.radii + " radii")
        return ", ".join(bits)


class _Amber(Engine):
    def available(self):
        return all(shutil.which(x) or _in_amberhome(x)
                   for x in ("tleap", "sander", "MMPBSA.py"))


class _OpenMM(Engine):
    def available(self):
        try:
            import openmm                                    # noqa: F401
            import openmm.app                                # noqa: F401
        except Exception:                                    # noqa: BLE001
            return False
        # OpenMM reads the topology AmberTools built, so tleap has to be
        # there even though sander and MMPBSA.py do not.
        return bool(shutil.which("tleap") or _in_amberhome("tleap"))


class _Schrodinger(Engine):
    # One MM-GBSA job takes 8 PSP_PLOP tokens out of a pool of 25, so three
    # is the most that can run at once and a Glide job elsewhere takes five
    # more. The count is read before a job is dispatched rather than when it
    # is queued, because the pool moves while a run waits.
    TOKENS_PER_JOB = 8
    FEATURE = "PSP_PLOP"

    def available(self):
        root = os.environ.get("SCHRODINGER", "")
        return bool(root and os.path.isfile(os.path.join(root, "prime_mmgbsa")))

    def licences(self, timeout=30):
        """(free, total) PSP_PLOP tokens, or (None, None) if unreadable.

        license_status is tried first because it is what answers on this
        site; licadmin is the fallback for installations that only have it.
        Either way an unreadable pool returns None rather than a guess, and
        the caller treats that as "do not dispatch" rather than "plenty".
        """
        root = os.environ.get("SCHRODINGER", "")
        tries = [["license_status", "-v", "Schrodinger"]]
        if root:
            tries.append([os.path.join(root, "licadmin"), "STAT"])
        for cmd in tries:
            try:
                out = subprocess.run(cmd, capture_output=True, text=True,
                                     timeout=timeout).stdout
            except Exception:                                # noqa: BLE001
                continue
            free, total = _parse_tokens(out, self.FEATURE)
            if free is not None:
                return free, total
        return None, None

    def slots(self):
        """How many MM-GBSA jobs the licence pool can take right now."""
        free, _total = self.licences()
        if free is None:
            return None
        return max(0, free // self.TOKENS_PER_JOB)


def _in_amberhome(exe):
    home = os.environ.get("AMBERHOME", "")
    return bool(home and os.path.isfile(os.path.join(home, "bin", exe)))


def _parse_tokens(text, feature):
    """(free, total) for one feature, read from license_status output.

    The line looks like

        Users of PSP_PLOP: (Total of 25 licenses issued; Total of 0 licenses
        in use; Total of 25 licenses available)

    and it states the free count directly, so that is what is read rather
    than subtracting. The colon after the name matters: without it the
    pattern also matches PSP_PLOP_MEMBRANE, which is a different pool.
    """
    if not text:
        return None, None
    m = re.search(
        r"Users of %s:\s*\(Total of (\d+) licenses? issued;"
        r".*?Total of (\d+) licenses? available" % re.escape(feature),
        text, re.I | re.S)
    if m:
        return int(m.group(2)), int(m.group(1))
    # Older output says how many are in use rather than how many are free.
    m = re.search(
        r"Users of %s:\s*\(Total of (\d+) licenses? issued;"
        r".*?(\d+) licenses? in use" % re.escape(feature),
        text, re.I | re.S)
    if m:
        total, used = int(m.group(1)), int(m.group(2))
        return max(0, total - used), total
    return None, None


AMBER = _Amber(
    name="amber", label="AmberTools",
    forcefield="ff14SB", solvation="OBC(II), igb=5", radii="mbondi2",
    column="mmgbsa", comparable=("amber", "openmm"),
    note="The reference. Every number this project has published so far "
         "came from here, so it is the one an older run can be compared "
         "with.")

OPENMM = _OpenMM(
    name="openmm", label="OpenMM",
    forcefield="ff14SB", solvation="OBC(II)", radii="mbondi2",
    column="mmgbsa_openmm", comparable=("amber", "openmm"),
    note="The same force field and the same topology as AmberTools, "
         "minimised to a gradient rather than to a step count. Roughly 30 "
         "times faster. Its values sit about 30 kcal/mol above AmberTools' "
         "on this project's poses while giving the same order, so read the "
         "rank and not the number.")

SCHRODINGER = _Schrodinger(
    name="schrodinger", label="Schrodinger Prime",
    forcefield="OPLS_2005", solvation="VSGB 2.0", radii="n/a",
    column="mmgbsa_schrodinger", comparable=("schrodinger",),
    note="A different force field and a different solvent model, so its "
         "numbers are not on the same scale as the other two and are never "
         "merged with them. One job takes 8 of the 25 PSP_PLOP tokens, so "
         "at most three run at once, and fewer while Glide is running.")

ENGINES = {e.name: e for e in (AMBER, OPENMM, SCHRODINGER)}
ORDER = ["amber", "openmm", "schrodinger"]
DEFAULT = "amber"


def get(name):
    """The engine by that name, defaulting to AmberTools."""
    return ENGINES.get(str(name or DEFAULT).strip().lower(), ENGINES[DEFAULT])


def installed():
    """The engines this machine can actually run."""
    return [n for n in ORDER if ENGINES[n].available()]


def catalogue():
    """Every engine, described, with whether it is available here."""
    return [ENGINES[n].as_dict() for n in ORDER]


def columns():
    """Every column an engine can write, in a fixed order."""
    return [ENGINES[n].column for n in ORDER]


def comparable(a, b):
    """Whether two engines' ABSOLUTE values may be put side by side."""
    return b in get(a).comparable and a in get(b).comparable
