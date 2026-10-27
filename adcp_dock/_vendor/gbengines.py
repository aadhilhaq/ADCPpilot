"""The MM-GBSA engines, and what each one's number actually means.

Three programs can rescore a pose, and they do not compute the same thing.
Each carries its own force field and its own implicit solvent model, so the
numbers are on different scales and an absolute value from one cannot be put
beside an absolute value from another:

    AmberTools   ff14SB + OBC(II), igb=5, mbondi2 radii
    OpenMM       ff14SB + OBC(II), the same model on the same topology
    Schrodinger  S-OPLS + VSGB 2.0, or OPLS_2005 when the single force
                 field license is taken, which is a different answer again

AmberTools and OpenMM at least share a force field and a topology, so the
difference between them is the solver, the minimization and the generalized
Born implementation. Schrodinger shares none of it.

What was measured, on topologies that were byte for byte the same. The
variable is how far the pose is relaxed before the energy is taken:

    OpenMM on the geometry sander minimized      1.3 kcal/mol from
                                                 MMPBSA.py, order unchanged
    OpenMM given 1000 more steps from there      8.7 out
    OpenMM minimizing the raw pose itself,       29.5 out
    1000 steps
    OpenMM minimizing the raw pose to            29.6 out, and on a 19
    convergence                                  residue peptide the pose
                                                 order reversed

So the two energy models agree and the MINIMIZATION was the entire
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

WHAT SWAPPING THE ENGINE DOES NOT FIX. The docked poses of one peptide can
differ from each other by more than half the range across the peptides being
compared, so the noise within a construct is a large share of the signal
between constructs. That is a property of the poses, not of the rescoring,
and every engine here inherits it. More poses help; a different solver does
not.
"""
from __future__ import annotations

import glob
import importlib.util
import os
import re
import shutil
import subprocess
import sys
import time


#: Where the Schrodinger suite lives when $SCHRODINGER is not set.
#:
#: Requiring the variable reported "not installed" on machines that have the
#: whole suite: clusters commonly install it site wide as a module, and the
#: module that exports $SCHRODINGER is not loaded inside a batch job unless
#: the job loads it. So the install is looked for as well as asked for.
SCHRODINGER_ROOTS = (
    "/sw/restricted/Schrodinger",      # a site install, one folder a release
    "/sw/restricted/schrodinger",
    "/sw/eb/sw/Schrodinger",
    "/sw/software/Schrodinger",
    "/software/Schrodinger",
    "/opt/schrodinger",
    "/opt/Schrodinger",
    "/usr/local/schrodinger",
    "/usr/local/Schrodinger",
)

#: Directories that hold one installation each, searched one level down as
#: well: <root>/<release> is the usual shape of a site install, and a
#: cluster that keeps releases under a vendor directory is the same shape one
#: level further in.
SCHRODINGER_GLOBS = (
    "/sw/restricted/*/Schrodinger*",
    "/sw/restricted/*/schrodinger*",
    # Some sites keep restricted software one level further in, under the
    # group that licenses it: /sw/restricted/<group>/sw/Schrodinger/<release>.
    # The shape above stops short of it, so Prime would show as not installed
    # while "module load Schrodinger" finds it at once.
    "/sw/restricted/*/sw/Schrodinger*",
    "/sw/restricted/*/sw/schrodinger*",
)
# Not "/sw/eb/sw/*/Schrodinger*": EasyBuild keeps a package at
# /sw/eb/sw/<Name>/<version>, which SCHRODINGER_ROOTS already names, and
# listing every package directory under /sw/eb/sw to look one level inside
# each took long enough on a cold compute node to time out the page.


def _version_key(path):
    """Sort key for a release directory, newest last.

    Parsed rather than sorted as text: 2025-10 sorts BEFORE 2025-3 as a
    string, which would quietly pick an older suite once a tenth release
    exists.
    """
    name = os.path.basename(path.rstrip("/"))
    m = re.match(r"(\d{4})-(\d+)", name)
    return (int(m.group(1)), int(m.group(2))) if m else (0, 0)


def _from_module_system():
    """Where an LMOD module would point, without loading it.

    Clusters install this as a module, and the module is what sets
    $SCHRODINGER. Inside a batch job nothing has loaded it, so the modulefile
    is read for the prefix it would export rather than the variable being
    relied on. This is what makes the search work on a cluster whose install
    path nobody wrote down here.
    """
    out = []
    mp = os.environ.get("MODULEPATH", "")
    for base in [d for d in mp.split(":") if d]:
        # One level in, where LMOD keeps a package's module files
        # (<modulepath>/<Name>/<version>.lua). A recursive search here walked
        # whole module trees on every call.
        for f in glob.glob(os.path.join(base, "[Ss]chrodinger*", "*"))[:40]:
            if os.path.isdir(f):
                continue
            try:
                text = open(f, errors="replace").read()
            except OSError:
                continue
            # setenv("SCHRODINGER","/path"), setenv SCHRODINGER /path, and
            # the prefix a Lua modulefile builds its paths from.
            for m in re.finditer(
                    r'(?:SCHRODINGER["\'\s,]+|prefix\s*=\s*")(/[^"\'\s,)]+)',
                    text):
                out.append(m.group(1))
            # Some sites' modulefiles name the install in a local variable
            # and then setenv("SCHRODINGER", root), which the pattern above
            # cannot follow. Any quoted absolute path that names the suite
            # is a candidate; newest_under decides which ones are installs.
            for m in re.finditer(r'["\'](/[^"\'\s]*[Ss]chrodinger[^"\'\s]*)["\']',
                                 text):
                if m.group(1) not in out:
                    out.append(m.group(1))
    return out


_FOUND = {}          # $SCHRODINGER at the time -> what the search found


def find_schrodinger(trail=None):
    """The newest Schrodinger install this account can actually run, or None.

    `trail`, if a list is passed, is filled with what was tried, so a page
    that reports "not installed" can say where it looked instead of leaving
    somebody to guess. That matters because clusters install the suite in
    different places, and no list of paths can know them all.

    Executability is part of the question, not an afterthought: the suite is
    restricted software, so a directory can be there and readable while this
    account is not in the group that may run it. Offering the engine in that
    case produces a run that fails at the first pose.
    """
    def note(what, ok=False):
        if trail is not None:
            trail.append(("found " if ok else "looked at ") + what)

    # Remembered for the life of the process: an install does not move while
    # a server runs, and the page asks on every load. The report function
    # passes a trail to get a fresh search it can show.
    key = os.environ.get("SCHRODINGER", "")
    if trail is None and key in _FOUND:
        return _FOUND[key]

    def usable(d):
        return os.access(os.path.join(d, "prime_mmgbsa"), os.X_OK)

    def newest_under(base):
        hits = [base] if usable(base) else []
        try:
            names = os.listdir(base)
        except OSError:
            names = []
        hits += [os.path.join(base, n) for n in names
                 if usable(os.path.join(base, n))]
        return sorted(hits, key=_version_key)[-1] if hits else None

    def done(value):
        if trail is None:
            _FOUND[key] = value
        return value

    root = key
    if root:
        note("$SCHRODINGER=%s" % root, usable(root))
        if usable(root):
            return done(root)
    else:
        note("$SCHRODINGER (not set)")

    # Already on PATH, which is the case once the module is loaded.
    exe = shutil.which("prime_mmgbsa")
    if exe:
        cand = os.path.dirname(os.path.realpath(exe))
        note("prime_mmgbsa on PATH at %s" % cand, True)
        return done(cand)
    note("prime_mmgbsa on PATH (not there)")

    # Cheapest first, and stop at the first place that has one: the known
    # roots, then what the module system would set, then the wider patterns.
    for base in SCHRODINGER_ROOTS:
        if not os.path.isdir(base):
            continue
        got = newest_under(base)
        note("%s%s" % (base, "" if got else " (no usable install)"), bool(got))
        if got:
            return done(got)
    for cand in _from_module_system():
        got = newest_under(cand)
        note("%s, from the module system" % cand, bool(got))
        if got:
            return done(got)
    for pat in SCHRODINGER_GLOBS:
        for base in sorted(glob.glob(pat)):
            got = newest_under(base)
            note("%s%s" % (base, "" if got else " (no usable install)"),
                 bool(got))
            if got:
                return done(got)
    return done(None)


def schrodinger_search_report():
    """Where find_schrodinger looked, for a page that has to explain itself."""
    trail = []
    got = find_schrodinger(trail)
    return {"found": got, "looked": trail}


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


def _have(tool):
    """Whether an executable can be found the way the engines will look.

    PATH, then the running interpreter's own bin, then $AMBERHOME/bin: the
    same order as prepare.find_tool, which is what mmgbsa uses to run them.
    shutil.which alone said AmberTools and OpenMM were missing whenever the
    environment's python was run without its bin on PATH, which is an
    ordinary way to start a job or a monitor. mmgbsa then asked this module
    whether the engine it had just found was installed, was told no, and
    skipped the rescoring altogether.
    """
    if shutil.which(tool):
        return True
    here = os.path.dirname(os.path.abspath(sys.executable))
    cand = os.path.join(here, tool)
    if os.path.isfile(cand) and os.access(cand, os.X_OK):
        return True
    return _in_amberhome(tool)


class _Amber(Engine):
    def available(self):
        return all(_have(x) for x in ("tleap", "sander", "MMPBSA.py"))


class _OpenMM(Engine):
    def available(self):
        # Found, not imported. Importing OpenMM reads hundreds of files, and
        # on a busy shared filesystem that took longer than the page waits
        # for its first answer, which asks this of every engine.
        try:
            if importlib.util.find_spec("openmm") is None:
                return False
        except (ImportError, ValueError):
            return False
        # tleap builds the topology and sander produces the geometry, and
        # this engine needs both. It used to require only tleap, on the
        # assumption that OpenMM would minimize for itself; that was measured
        # to put the energy about 30 kcal/mol out, so the geometry is sander's
        # and sander is a requirement. MMPBSA.py is the part not needed.
        return all(_have(x) for x in ("tleap", "sander"))


class _Schrodinger(Engine):
    """Prime MM-GBSA, which needs two different licenses to do its best work.

    The rescoring itself takes 8 PSP_PLOP tokens a job, so the size of that
    pool, less what other Schrodinger jobs hold, sets how many run at once.
    But the force field is a separate license: prime_mmgbsa's own help says
    the default is S-OPLS "if the proper license is present", and
    OPLS_VERSION=OPLS_2005 is how you step down from it deliberately.

    FFLD_OPLS_MAIN can have fewer seats than the rescoring pool has room
    for jobs, often just one. A run that assumed every job would get the
    modern force field would then queue jobs that quietly fall back to a
    twenty year old one without saying so.

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
        """The license server's report, or "" if it cannot be read.

        `timeout` is the TOTAL budget, shared between the commands tried,
        not a budget each. It used to be per command, so two attempts at 5
        seconds could block a caller for 10: the page asks for this while it
        loads, and on a compute node that cannot reach the license server
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

    #: The last license reading, and when it was taken. Reading it means
    #: talking to a license server, which is fast when that works and a
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
                             "why": "the license pool could not be read"})
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
    note="The reference, and the default. Runs made before other engines "
         "could be chosen used it, so it is the one an older run can be "
         "compared with. Its dG includes the surface term (ESURF), which "
         "is part of what the SA in MM-GBSA means.")

OPENMM = _OpenMM(
    name="openmm", label="OpenMM",
    forcefield="ff14SB", solvation="OBC(II)", radii="mbondi2",
    column="mmgbsa_openmm", comparable=("openmm",),
    note="A second opinion on the same pose. It replaces MMPBSA.py and not "
         "sander: the geometry is the one sander minimized, and OpenMM "
         "computes that geometry's energy with its own generalized Born "
         "code. In testing it lands within about 1.3 kcal/mol of "
         "AmberTools with the pose order unchanged, so it is a check on the "
         "solver rather than a different estimate. Letting it minimize for "
         "itself was tried and put the energy about 30 kcal/mol out, so it "
         "does not.")

SCHRODINGER = _Schrodinger(
    name="schrodinger", label="Schrodinger Prime",
    # S-OPLS when the single FFLD_OPLS_MAIN seat is free, OPLS_2005 when it
    # is not. Which one a job got is recorded with the result rather than
    # assumed from here.
    forcefield="S-OPLS, or OPLS_2005 without the license",
    solvation="VSGB 2.0", radii="n/a",
    column="mmgbsa_schrodinger", comparable=("schrodinger",),
    note="A different force field and a different solvent model, so its "
         "numbers are not on the same scale as the other two and are never "
         "merged with them: on one test pose it gave -81.79 kcal/mol where "
         "AmberTools gave -48.07. The value is Prime's no-strain binding "
         "energy, taken like the others from the complex's own geometry; its "
         "headline figure also subtracts the strain of a separately relaxed "
         "free peptide, which grew from 9 to 59 kcal/mol between a 5 and a 15 "
         "residue peptide. Two licenses are needed: each job takes 8 PSP_PLOP "
         "tokens, and S-OPLS also needs a seat of its own, which a site often "
         "has fewer of. A job that finds no seat free runs on OPLS_2005, and "
         "which force field a result came from is recorded with it.")

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
