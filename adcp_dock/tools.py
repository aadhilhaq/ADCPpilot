"""Locate the external executables, without hard-coding any one installation.

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
    "obabel": "converting a peptide to PDB when no template is available",
}

# Docking cannot run without these. The rest are optional: the Amber tools
# are needed only for MM-GBSA rescoring, and obabel only as a fallback.
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

    def missing_required(self):
        return [n for n in REQUIRED if not self._found.get(n)]

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
    return Tools(found, roots)


def amber_home(tools=None):
    """AMBERHOME, worked out from tleap's location when it is not set."""
    v = os.environ.get("AMBERHOME")
    if v and os.path.isfile(os.path.join(v, "bin", "tleap")):
        return v
    t = (tools.get("tleap") if tools else None) or shutil.which("tleap")
    if t:
        return os.path.dirname(os.path.dirname(t))
    return None
