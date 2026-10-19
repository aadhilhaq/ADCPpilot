"""Locate the external executables, without hard-coding any one installation.

By Aadhil Haq

The docking comes from the ADCP 1.1 suite, a conda environment the Scripps
Research Institute's own recipe builds (adcp11_install): adcp, agfr,
prepare_receptor and reduce are taken from its bin folder and from nowhere
else, and MSMS from its MSMS library (surface.find_msms). Where the suite is
looked for is adcp11.candidates: ADCPSUITE_HOME, the roots passed in, the
installation's tools folder and the one beside this checkout.

The rest (pdb2pqr for protonation at a pH, AmberTools for MM-GBSA) is
searched for in this interpreter's own environment, AMBERHOME, PATH, then
conventional subdirectories under any roots passed in.

On a Mac, where the Scripps Research Institute builds no ADCP 1.1, the
docking programs, reduce and msms come from ADFRsuite 1.0 instead (mac10);
that is asked of mac10.engaged(), which is False on every other system.

A missing tool is reported by name, with what it is needed for, rather than
as a FileNotFoundError from a job that is already running.
"""
import os
import shutil
import sys

# what each one is for, used in the message when it cannot be found
NEEDED = {
    "adcp": "the peptide docking search",
    "agfr": "building the grid maps ADCP docks into",
    "prepare_receptor": "writing the receptor as PDBQT",
    "pdb2pqr": "protonating the receptor at a chosen pH",
    "tleap": "building Amber topologies for MM-GBSA",
    "sander": "the minimisation before MM-GBSA",
    "MMPBSA.py": "the MM-GBSA calculation itself",
    "ante-MMPBSA.py": "splitting the complex topology for MM-GBSA",
}

#: The programs taken from the ADCP 1.1 suite's bin folder, never from PATH:
#: a copy found there may be another version, from another installation.
SUITE_BINS = ("adcp", "agfr", "prepare_receptor")

# Two more from the suite, optional: the page says what is missing. reduce is
# its bin/reduce (AmberTools has a reduce of its own, which is not this one),
# the same program ADFR's own preparation and ADapt run; MSMS is its MSMS
# library, run through msms_mslib.py (msms_command).
SUITE_PROGRAMS = {
    "reduce": "protonating the receptor (the default) and, with it, each "
              "complex for MM-GBSA",
    "msms": "the molecular surfaces the viewers draw",
}

# Docking cannot run without these. The rest are optional: the Amber tools
# are needed only for MM-GBSA rescoring.
REQUIRED = ("adcp", "agfr", "prepare_receptor", "pdb2pqr")
MMGBSA = ("tleap", "sander", "MMPBSA.py", "ante-MMPBSA.py")

# subdirectories worth searching under any root passed in
SUBDIRS = ("bin", "third_party/bin", "amber/bin", "amber22/bin")


class Missing(Exception):
    """A tool the requested work cannot proceed without."""


class Tools:
    def __init__(self, found, roots):
        self._found = found
        self.roots = list(roots)

    def __getattr__(self, name):
        try:
            return self._found[name]
        except KeyError:
            raise AttributeError(name)

    def get(self, name):
        return self._found.get(name)

    def require(self, *names):
        """Every one of these, or a message naming the ones that are absent."""
        gone = [n for n in names if not self._found.get(n)]
        if gone:
            raise Missing(
                "cannot run: " + "; ".join(
                    "%s not found (needed for %s)" % (n, NEEDED.get(n, "this step"))
                    for n in gone)
                + ".\nLooked on PATH and under: "
                + (", ".join(self.roots) if self.roots else "(no extra roots given)")
                + "\nSet ADCPSUITE_HOME and AMBERHOME, or pass --tools-root.")
        return [self._found[n] for n in names]

    def report(self):
        return [(n, self._found.get(n) or "") for n in sorted(NEEDED)]

    def suite_program(self, name):
        """reduce, or another program, from the ADCP 1.1 suite's bin folder,
        or None. Never from PATH: see SUITE_PROGRAMS."""
        return suite_program(name, self, self.roots)

    def missing_required(self):
        """The required tools not found."""
        return [n for n in REQUIRED if not self._found.get(n)]

    def can_rescore(self):
        """True when every tool MM-GBSA needs was found."""
        return all(self._found.get(n) for n in MMGBSA)


def suite_home(roots=()):
    """The ADCP 1.1 suite's folder, when it is installed with its own AGFR,
    or None."""
    try:
        from . import adcp11
        exe, script = adcp11.find_agfr(list(roots))
        return adcp11.home_of(exe) if exe and script else None
    except Exception:                                          # noqa: BLE001
        return None


def _executable(p):
    return bool(p) and os.path.isfile(p) and os.access(p, os.X_OK)


def discover(roots=(), env=True):
    """Tools(), with every executable this package knows how to use."""
    from . import mac10
    found = {}
    home = mac10.home() or suite_home(roots)
    if home:
        for n in SUITE_BINS:
            p = os.path.join(home, "bin", n)
            if _executable(p):
                found[n] = p
    if env and os.environ.get("AMBERHOME"):
        for n in MMGBSA:
            p = os.path.join(os.environ["AMBERHOME"], "bin", n)
            if _executable(p):
                found[n] = p
    # This interpreter's own environment is searched before PATH. pdb2pqr
    # is itself a Python program, and the only copy guaranteed to match
    # this environment's numpy is the one installed beside the
    # interpreter; a copy found on PATH may belong to another environment
    # and fail to import. It is also where a batch script that calls
    # <env>/bin/python without putting <env>/bin on PATH finds its tools.
    rest = [n for n in NEEDED if n not in SUITE_BINS]
    _here = os.path.dirname(os.path.abspath(sys.executable))
    for n in rest:
        if n not in found and _executable(os.path.join(_here, n)):
            found[n] = os.path.join(_here, n)
    for n in rest:
        if n not in found:
            p = shutil.which(n)
            if p:
                found[n] = p
    for root in roots:
        if not root:
            continue
        for sub in SUBDIRS:
            d = os.path.join(root, sub)
            if not os.path.isdir(d):
                continue
            for n in rest:
                if n not in found and _executable(os.path.join(d, n)):
                    found[n] = os.path.join(d, n)
    return Tools(found, roots)


def suite_program(name, tools=None, roots=()):
    """A program in the ADCP 1.1 suite's bin folder (SUITE_PROGRAMS), or
    None when the suite is not installed or has no such program. On a Mac,
    ADFRsuite 1.0's (mac10)."""
    from . import mac10
    if mac10.engaged():
        return mac10.program(name)
    roots = list(roots) or list(getattr(tools, "roots", None) or [])
    home = suite_home(roots)
    p = os.path.join(home, "bin", name) if home else None
    return p if _executable(p) else None


def msms_command(tools=None, roots=()):
    """How MSMS is run here, as surface.find_msms finds it: $MSMS_EXE, an
    msms program in a root's bin folder, or the ADCP 1.1 suite's MSMS
    library (mslib) through msms_mslib.py and the suite's own Python, which
    take the program's arguments and write its files (surface.run_msms runs
    either). None when there is none of them."""
    from . import mac10, surface
    if mac10.engaged():
        # ADFRsuite 1.0's msms program, in its bin folder.
        return surface.find_msms([mac10.home()])
    roots = list(roots) or list(getattr(tools, "roots", None) or [])
    home = suite_home(roots)
    return surface.find_msms(([home] if home else []) + roots)


def amber_home(tools=None):
    """AMBERHOME, worked out from tleap's location when it is not set."""
    v = os.environ.get("AMBERHOME")
    if v and os.path.isfile(os.path.join(v, "bin", "tleap")):
        return v
    t = (tools.get("tleap") if tools else None) or shutil.which("tleap")
    if t:
        return os.path.dirname(os.path.dirname(t))
    return None
