"""The MM-GBSA engines, and what each one's number actually means.

Three programs can rescore a pose, and they do not compute the same thing.
Each carries its own force field and its own implicit solvent model, so the
numbers are on different scales and an absolute value from one cannot be put
beside an absolute value from another:

    AmberTools   ff14SB + OBC(II), igb=5, mbondi2 radii
    OpenMM       ff14SB + OBC(II), the same model on the same topology
    Schrodinger  S-OPLS + VSGB 2.0, or OPLS_2005 when the single force
                 field licence is taken, which is a different answer again

AmberTools and OpenMM at least share a force field and a topology, so the
difference between them is the solver, the minimisation and the generalised
Born implementation. Schrodinger shares none of it.

What was measured, on topologies that were byte for byte the same. The
variable is how far the pose is relaxed before the energy is taken:

    OpenMM on the geometry sander minimised      1.3 kcal/mol from
                                                 MMPBSA.py, order unchanged
    OpenMM given 1000 more steps from there      8.7 out
    OpenMM minimising the raw pose itself,       29.5 out
    1000 steps
    OpenMM minimising the raw pose to            29.6 out, and on a 19
    convergence                                  residue peptide the pose
                                                 order reversed

So the two energy models agree and the MINIMISATION was the entire
difference. sander's maxcyc=1000 is a shallow relaxation while OpenMM's
L-BFGS converges well inside 1000 iterations: capping it at 1000 steps and
letting it run to convergence land 29.5 and 29.6 out, which is the same
answer, so matching the step count cannot recover the reference. The extra
relaxation loosens the interface and weakens the energy.

The OpenMM engine therefore scores sander's geometry rather than making its
own, and declines to run at all where sander is absent. Schrodinger is a
different matter: a different force field and a different solvent model, on
its own scale, and nothing calibrates it against these two.

Each engine writes its own column, none is merged into another, none is
compared with another by value, and the table says which force field and
which solvent model produced each one. A run may carry all three.

WHAT SWAPPING THE ENGINE DOES NOT FIX. On 1BE9, five docked poses of one
peptide span 21.7 kcal/mol, against a 37.6 kcal/mol range across all twenty
constructs: the noise within a construct is 58 per cent of the signal
between constructs. That is a property of the poses, not of the rescoring,
and every engine here inherits it. More poses would help; a different solver
does not.
"""
from __future__ import annotations

import glob
import os
import re
import shutil
import subprocess
import time


#: Where the Schrodinger suite lives when $SCHRODINGER is not set.
#:
#: Requiring the variable reported "not installed" on machines that have the
#: whole suite: it is site-installed on both clusters this runs on, and the
#: module that exports $SCHRODINGER is not loaded inside a batch job unless
#: the job loads it. So the install is looked for as well as asked for.
SCHRODINGER_ROOTS = (
    "/sw/restricted/Schrodinger",      # TAMU HPRC, ACES and Grace
    "/sw/eb/sw/Schrodinger",
    "/opt/schrodinger",
    "/usr/local/schrodinger",
)


def _version_key(path):
    """Sort key for a release directory, newest last.

    Parsed rather than sorted as text: 2025-10 sorts BEFORE 2025-3 as a
    string, which would quietly pick an older suite once a tenth release
    exists.
    """
    name = os.path.basename(path.rstrip("/"))
    m = re.match(r"(\d{4})-(\d+)", name)
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0)


def find_schrodinger():
    """The newest Schrodinger install this account can actually run, or None.

    Executability is part of the question, not an afterthought: the suite is
    restricted software on these clusters, so a directory can be there and
    readable while this account is not in the group that may run it. Offering
    the engine in that case produces a run that fails at the first pose.
    """
    root = os.environ.get("SCHRODINGER", "")
    if root and os.access(os.path.join(root, "prime_mmgbsa"), os.X_OK):
        return root
    found = []
    for base in SCHRODINGER_ROOTS:
        if os.access(os.path.join(base, "prime_mmgbsa"), os.X_OK):
            found.append(base)
        for d in glob.glob(os.path.join(base, "*")):
            if os.access(os.path.join(d, "prime_mmgbsa"), os.X_OK):
                found.append(d)
    if not found:
        return None
    return sorted(found, key=_version_key)[-1]


class Engine:
    """One way to compute an MM-GBSA energy, described honestly.

    `column` is where its numbers land, so two engines never overwrite each
    other. `comparable` names the engines whose ABSOLUTE values can sit in
    the same sentence as this one's. On the measurements below no engine
    passes that test against any other, so each one's list holds only
    itself; what the engines share is the RANK, and only within one target.
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
        # tleap builds the topology and sander produces the geometry, and
        # this engine needs both. It used to require only tleap, on the
        # assumption that OpenMM would minimise for itself; that was measured
        # to put the energy about 30 kcal/mol out, so the geometry is sander's
        # and sander is a requirement. MMPBSA.py is the part not needed.
        return all(shutil.which(x) or _in_amberhome(x)
                   for x in ("tleap", "sander"))


class _Schrodinger(Engine):
    """Prime MM-GBSA, which needs two different licences to do its best work.

    The rescoring itself takes 8 PSP_PLOP tokens out of 25, so three jobs
    fit and fewer while Glide is running. But the force field is a separate
    licence: prime_mmgbsa's own help says the default is S-OPLS "if the
    proper license is present", and OPLS_VERSION=OPLS_2005 is how you step
    down from it deliberately.

    On this site FFLD_OPLS_MAIN has ONE seat. So the honest capacity is
    three jobs on OPLS_2005, or one job on S-OPLS, and a run that assumed
    three of the modern force field would queue two that quietly fall back
    to a twenty year old one without saying so.

    Both pools are therefore read, and the force field a job will actually
    get is reported rather than assumed.
    """

    TOKENS_PER_JOB = 8
    FEATURE = "PSP_PLOP"
    FF_FEATURE = "FFLD_OPLS_MAIN"
    FF_MODERN = "S-OPLS"
    FF_LEGACY = "OPLS_2005"

    def available(self):
        return find_schrodinger() is not None

    def root(self):
        """Where prime_mmgbsa is, whether or not $SCHRODINGER was set."""
        return find_schrodinger()

    def _status(self, timeout=30):
        """The licence server's report, or "" if it cannot be read.

        `timeout` is the TOTAL budget, shared between the commands tried,
        not a budget each. It used to be per command, so two attempts at 5
        seconds could block a caller for 10: the page asks for this while it
        loads, and on a compute node that cannot reach the licence server
        both attempts ran to the full timeout and the whole request stalled.
        """
        root = find_schrodinger() or ""
        tries = [["license_status", "-v", "Schrodinger"]]
        if root:
            tries.append([os.path.join(root, "licadmin"), "STAT"])
        left = float(timeout)
        for cmd in tries:
            if left <= 0:
                break
            started = time.monotonic()
            try:
                out = subprocess.run(cmd, capture_output=True, text=True,
                                     timeout=left).stdout
            except Exception:                                # noqa: BLE001
                out = ""
            left -= time.monotonic() - started
            if out and "Users of" in out:
                return out
        return ""

    def licences(self, timeout=30):
        """(free, total) PSP_PLOP tokens, or (None, None) if unreadable."""
        return _parse_tokens(self._status(timeout), self.FEATURE)

    def forcefield_licences(self, timeout=30):
        """(free, total) seats for the modern force field."""
        return _parse_tokens(self._status(timeout), self.FF_FEATURE)

    def slots(self, timeout=30):
        """How many MM-GBSA jobs the rescoring pool can take right now."""
        free, _total = self.licences(timeout)
        if free is None:
            return None
        return max(0, free // self.TOKENS_PER_JOB)

    #: The last licence reading, and when it was taken. Reading it means
    #: talking to a licence server, which is fast when that works and a
    #: multi-second stall when it does not, so a page that asks on every
    #: load asks this instead.
    _CACHE = {"at": 0.0, "value": None}
    CACHE_SECONDS = 60.0

    def plan(self, timeout=30, max_age=None):
        """What can actually run, and with which force field.

        Returned rather than decided here, because a caller may reasonably
        prefer three OPLS_2005 jobs to one S-OPLS job or the other way
        round, and that is not this module's choice to make.
        """
        age = self.CACHE_SECONDS if max_age is None else max_age
        now = time.monotonic()
        cached = self._CACHE["value"]
        if cached is not None and (now - self._CACHE["at"]) < age:
            return cached
        out = self._status(timeout)
        free, total = _parse_tokens(out, self.FEATURE)
        ff_free, ff_total = _parse_tokens(out, self.FF_FEATURE)
        def remember(value):
            self._CACHE["at"], self._CACHE["value"] = now, value
            return value
        if free is None:
            return remember({"readable": False, "jobs": None,
                             "modern_jobs": None, "forcefield": None,
                             "why": "the licence pool could not be read"})
        jobs = max(0, free // self.TOKENS_PER_JOB)
        modern = 0 if ff_free is None else min(jobs, ff_free)
        return remember({
                "readable": True,
                "jobs": jobs, "tokens_free": free, "tokens_total": total,
                "modern_jobs": modern,
                "ff_free": ff_free, "ff_total": ff_total,
                "forcefield": self.FF_MODERN if modern else self.FF_LEGACY,
                "why": ("%d job(s) fit; %d of them can use %s, the rest fall "
                        "back to %s"
                        % (jobs, modern, self.FF_MODERN, self.FF_LEGACY))})


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
    column="mmgbsa", comparable=("amber",),
    note="The reference. Every number this project has published so far "
         "came from here, so it is the one an older run can be compared "
         "with. Its dG includes the surface term (ESURF), which is part of "
         "what the SA in MM-GBSA means.")

OPENMM = _OpenMM(
    name="openmm", label="OpenMM",
    forcefield="ff14SB", solvation="OBC(II)", radii="mbondi2",
    column="mmgbsa_openmm", comparable=("openmm",),
    note="A second opinion on the same pose. It replaces MMPBSA.py and not "
         "sander: the geometry is the one sander minimised, and OpenMM "
         "computes that geometry's energy with its own generalised Born "
         "code. On this project's poses it lands within 1.3 kcal/mol of "
         "AmberTools with the pose order unchanged, so it is a check on the "
         "solver rather than a different estimate. Letting it minimise for "
         "itself was tried and put the energy about 30 kcal/mol out, so it "
         "does not.")

SCHRODINGER = _Schrodinger(
    name="schrodinger", label="Schrodinger Prime",
    # S-OPLS when the single FFLD_OPLS_MAIN seat is free, OPLS_2005 when it
    # is not. Which one a job got is recorded with the result rather than
    # assumed from here.
    forcefield="S-OPLS, or OPLS_2005 without the licence",
    solvation="VSGB 2.0", radii="n/a",
    column="mmgbsa_schrodinger", comparable=("schrodinger",),
    note="A different force field and a different solvent model, so its "
         "numbers are not on the same scale as the other two and are never "
         "merged with them: on the same 1BE9 pose it gave -72.88 kcal/mol "
         "against AmberTools' -48.07. Two licences are needed and they do not have "
         "the same depth: the rescoring takes 8 of 25 PSP_PLOP tokens, so "
         "three jobs fit, while the modern force field has a single seat. "
         "Three jobs on OPLS_2005 or one on S-OPLS, and which one a result "
         "came from is recorded with it.")

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
