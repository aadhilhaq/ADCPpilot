"""Locate the external executables, without hard-coding any one installation.

By Aadhil Haq

ADCP and agfr ship inside ADFRsuite, which is unpacked wherever the
installer chose. Search order: this interpreter's own environment, the
environment variables ADFRSUITE_HOME, AMBERHOME and ADCP_DOCK_TOOLS,
PATH, then conventional subdirectories under any roots passed in.

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

# Programs that ship inside ADFRsuite and are taken from there and nowhere
# else: another copy on PATH may be another version (AmberTools has its own
# reduce), and the point of using these is that they are the ones ADFR's own
# preparation and ADapt use. Optional: the page says what is missing.
ADFR_PROGRAMS = {
    "reduce": "protonating the receptor (the default) and, with it, each "
              "complex for MM-GBSA",
    "msms": "the molecular surfaces the viewers draw",
}

# Docking cannot run without these. The rest are optional: the Amber tools
# are needed only for MM-GBSA rescoring.
REQUIRED = ("adcp", "agfr", "prepare_receptor", "pdb2pqr")
MMGBSA = ("tleap", "sander", "MMPBSA.py", "ante-MMPBSA.py")

# subdirectories worth searching under any root passed in
SUBDIRS = ("bin", "ADFRsuite/bin", "third_party/ADFRsuite/bin",
           "third_party/bin", "adfrsuite/bin", "amber/bin", "amber22/bin")


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
                + "\nSet ADFRSUITE_HOME and AMBERHOME, or pass --tools-root.")
        return [self._found[n] for n in names]

    def report(self):
        return [(n, self._found.get(n) or "") for n in sorted(NEEDED)]

    def adfr_program(self, name):
        """reduce or msms from the ADFRsuite whose prepare_receptor, agfr or
        adcp this found, or from the ADCP 1.1 suite, or None. Never from
        PATH: see ADFR_PROGRAMS."""
        return adfr_program(name, self, self.roots)

    def missing_required(self):
        """The required tools not found. ADCP and AGFR count as found where
        the ADCP 1.1 suite is installed, which runs them with its own Python
        (adcp11.command, agfr_command) rather than through bin/adcp and
        bin/agfr; those names are kept for ADFRsuite 1.0's programs."""
        gone = [n for n in REQUIRED if not self._found.get(n)]
        if ("adcp" in gone or "agfr" in gone) and suite11(self.roots):
            gone = [n for n in gone if n not in ("adcp", "agfr")]
        return gone

    def can_rescore(self):
        """True when every tool MM-GBSA needs was found."""
        return all(self._found.get(n) for n in MMGBSA)


def _from_env():
    out = {}
    for var, names in (("ADFRSUITE_HOME", ("adcp", "agfr", "prepare_receptor")),
                       ("AMBERHOME", ("tleap", "sander", "MMPBSA.py",
                                      "ante-MMPBSA.py"))):
        root = os.environ.get(var)
        if not root:
            continue
        for n in names:
            p = os.path.join(root, "bin", n)
            if os.path.isfile(p) and os.access(p, os.X_OK):
                out[n] = p
    return out


def suite11(roots=()):
    """The ADCP 1.1 suite's folder, when it is installed with its own AGFR,
    or None."""
    try:
        from . import adcp11
        exe, script = adcp11.find_agfr(list(roots))
        return adcp11.home_of(exe) if exe and script else None
    except Exception:                                          # noqa: BLE001
        return None


def discover(roots=(), env=True):
    """Tools(), with every executable this package knows how to use."""
    found = _from_env() if env else {}
    # This interpreter's own environment is searched before PATH. pdb2pqr
    # is itself a Python program, and the only copy guaranteed to match
    # this environment's numpy is the one installed beside the
    # interpreter; a copy found on PATH may belong to another environment
    # and fail to import. It is also where a batch script that calls
    # <env>/bin/python without putting <env>/bin on PATH finds its tools.
    _here = os.path.dirname(os.path.abspath(sys.executable))
    for n in NEEDED:
        if n in found:
            continue
        p = os.path.join(_here, n)
        if os.path.isfile(p) and os.access(p, os.X_OK):
            found[n] = p
    for n in NEEDED:
        if n in found:
            continue
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
            for n in NEEDED:
                if n in found:
                    continue
                p = os.path.join(d, n)
                if os.path.isfile(p) and os.access(p, os.X_OK):
                    found[n] = p
    # Where ADFRsuite 1.0 is not installed, the receptor is written as PDBQT
    # by the ADCP 1.1 suite's prepare_receptor, the same AutoDockTools script
    # ported to Python 3. Its launcher runs whichever python is first on
    # PATH, so receptor.prepare_pdbqt runs it with the suite's bin first.
    if not found.get("prepare_receptor"):
        home = suite11(roots)
        if home:
            p = os.path.join(home, "bin", "prepare_receptor")
            if os.path.isfile(p) and os.access(p, os.X_OK):
                found["prepare_receptor"] = p
    return Tools(found, roots)


def adfr_program(name, tools=None, roots=()):
    """A program that ships in ADFRsuite's bin folder (ADFR_PROGRAMS), found
    beside the ADFRsuite programs this package already uses, or in the ADCP
    1.1 suite's bin folder, which carries the same build. None when neither
    has it."""
    tools = tools or discover(roots=roots)
    dirs = []
    for n in ("prepare_receptor", "agfr", "adcp"):
        p = tools.get(n)
        if p:
            d = os.path.dirname(os.path.realpath(p))
            if d not in dirs:
                dirs.append(d)
    home = os.environ.get("ADFRSUITE_HOME")
    if home:
        dirs.append(os.path.join(home, "bin"))
    try:
        from . import adcp11
        exe = adcp11.find(list(roots) or list(getattr(tools, "roots", [])))
        if exe:
            dirs.append(os.path.join(adcp11.home_of(exe), "bin"))
    except Exception:                                          # noqa: BLE001
        pass
    for d in dirs:
        p = os.path.join(d, name)
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return None


def msms_command(tools=None, roots=()):
    """How MSMS is run here: ADFRsuite's msms program, or, where only the
    ADCP 1.1 suite is installed, its MSMS library (mslib) through
    msms_mslib.py and the suite's own Python, which take the program's
    arguments and write its files (surface.run_msms runs either). None when
    neither is installed."""
    import glob
    from . import adcp11
    tools = tools or discover(roots=roots)
    prog = tools.adfr_program("msms")
    if prog:
        return prog
    exe = adcp11.find(list(roots) or list(getattr(tools, "roots", [])))
    if not exe:
        return None
    home = adcp11.home_of(exe)
    py = os.path.join(home, "bin", "python")
    lib = glob.glob(os.path.join(home, "lib", "python3*", "site-packages",
                                 "mslib", "__init__.py"))
    script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "msms_mslib.py")
    if lib and os.path.isfile(py) and os.path.isfile(script):
        return [py, script]
    return None


def amber_home(tools=None):
    """AMBERHOME, worked out from tleap's location when it is not set."""
    v = os.environ.get("AMBERHOME")
    if v and os.path.isfile(os.path.join(v, "bin", "tleap")):
        return v
    t = (tools.get("tleap") if tools else None) or shutil.which("tleap")
    if t:
        return os.path.dirname(os.path.dirname(t))
    return None
