"""ADCP 1.1: finding it, running it, and writing what it produces in the
layout the rest of this package reads.

By Aadhil Haq

ADCP 1.1 is installed by the Scripps Research Institute's own recipe as a
separate conda environment (called adcpsuite), with its own Python. It docks
standard, non-standard and D amino acids, writes every pose with all its
atoms in one multi-model file, and can minimize and re-rank the top poses
with OpenMM.

Everything downstream (the results table, the viewer, MM-GBSA, resume) reads
out_ranked_<n>.pdb and the mode table in modes.json. So after a docking the
multi-model file is split into those files, the mode table is written as
modes.json and ADCP's summary is copied to adcp.log; ADCP's own files are
kept beside them.
"""
from __future__ import annotations

import glob
import json
import os
import re
import shutil
import subprocess

from . import peptide

#: ADCP 1.1's executable, relative to the environment it is installed in.
BIN = os.path.join("bin", "adcp")
# The OpenMM table ADCP 1.1 prints after -nmin: model, OpenMM rank, ADCP rank,
# E(complex) - E(receptor), E(complex) - E(receptor) - E(peptide), affinity.
# ADCP 1.1's mode table: mode, affinity, ref. fnc, cluster size, rmsd stdv,
# energy stdv, best run.
_MODE11 = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?[\d.]+|NA)\s+(\d+)\s+")
_OMM = re.compile(r"^OMM Ranking:\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?[\d.]+)\s+"
                  r"(-?[\d.]+)\s+(-?[\d.]+)")


#: The folder this package is in: a checkout, or site-packages.
_HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def candidates(roots=()):
    """Where an ADCP 1.1 environment is looked for, in order: ADCPSUITE_HOME;
    each root given, as the suite itself, a folder holding it as adcpsuite,
    or a tools folder holding it as tools/docking/adcpsuite; the same inside
    the installation (ADCP_HOME) and beside it, which is where install.sh
    puts it; then the folder `adcp-dock local` installs it in, and Colab's."""
    out = []
    env = os.environ.get("ADCPSUITE_HOME")
    if env:
        out.append(env)
    inst = os.environ.get("ADCP_HOME")
    homes = ([inst] if inst else []) + [_HERE]
    bases = (list(roots) + homes
             + [os.path.dirname(os.path.abspath(h)) for h in homes])
    for root in bases:
        if root:
            out += [root, os.path.join(root, "adcpsuite"),
                    os.path.join(root, "tools", "docking", "adcpsuite")]
    out += [os.path.expanduser("~/.adcp-dock/adcpsuite"), "/content/adcpsuite"]
    seen, uniq = set(), []
    for d in out:
        if d not in seen:
            seen.add(d)
            uniq.append(d)
    return uniq


def at(home):
    """The ADCP 1.1 executable of the suite in exactly this folder, or None."""
    if not home:
        return None
    exe = os.path.join(home, BIN)
    site = glob.glob(os.path.join(home, "lib", "python3*", "site-packages",
                                  "ADCP", "runADCP.py"))
    return exe if os.path.isfile(exe) and site else None


def find(roots=()):
    """The ADCP 1.1 executable, or None."""
    for home in candidates(roots):
        exe = at(home)
        if exe:
            return exe
    return None


def home_of(exe):
    return os.path.dirname(os.path.dirname(os.path.abspath(exe)))


def command(exe):
    """How to start ADCP 1.1: the environment's own Python on runADCP.py.

    Not the bin/adcp launcher, which runs whichever python is first on PATH
    and passes its arguments unquoted through a shell.
    """
    home = home_of(exe)
    script = glob.glob(os.path.join(home, "lib", "python3*", "site-packages",
                                    "ADCP", "runADCP.py"))[0]
    return [os.path.join(home, "bin", "python"), script]


def version(exe):
    """ADCP's version number, e.g. "1.1.21". Read from its source, where
    utils.py sets BUILD_VERSION and runADCP.py prints "1.1.%d": asking ADCP
    imports its whole modeling stack, which takes half a minute on a node
    whose file cache is cold, and the page waits for the answer."""
    for pkg in glob.glob(os.path.join(home_of(exe), "lib", "python3*",
                                      "site-packages", "ADCP")):
        try:
            build = re.search(r"^BUILD_VERSION\s*=\s*(\d+)",
                              open(os.path.join(pkg, "utils.py")).read(), re.M)
            major = re.search(r'version="(\d+\.\d+)\.%d"',
                              open(os.path.join(pkg, "runADCP.py")).read())
        except OSError:
            continue
        if build and major:
            return "%s.%s" % (major.group(1), build.group(1))
    try:
        p = subprocess.run(command(exe) + ["--version"], capture_output=True,
                           text=True, timeout=120, env=_env(exe))
        m = re.search(r"\d+\.\d+\.\d+", p.stdout + p.stderr)
        return m.group(0) if m else ""
    except (OSError, subprocess.SubprocessError):
        return ""


def library(exe):
    """Names of the non-standard side chains ADCP 1.1 ships (L and D forms)."""
    names = set()
    for d in glob.glob(os.path.join(home_of(exe), "lib", "python3*",
                                    "site-packages", "ADCP", "data",
                                    "rotamers")):
        lst = os.path.join(d, "SWISS_ROTAMER_LIST.dat")
        if os.path.isfile(lst):
            for line in open(lst):
                if line.startswith("#") or not line.split():
                    continue
                parts = line.split()
                names.add(parts[0].upper())
                if len(parts) > 3:
                    names.add(parts[3].upper())
        for lib in glob.glob(os.path.join(d, "*.lib")):
            for line in open(lib, errors="replace"):
                m = re.match(r"^\s*([A-Z0-9]{1,4})\s*$", line)
                if m:
                    names.add(m.group(1))
    return names


def openmm_sets(exe):
    """{side chain: {"swiss": (inside, N-terminal, C-terminal), "sannerlab":
    (...)}} from ADCP's AVAILABLE_PARAMETER.dat: which OpenMM parameter set
    covers each non-standard side chain, and in which position."""
    out = {}
    for f in glob.glob(os.path.join(home_of(exe), "lib", "python3*", "site-packages",
                                    "ADCP", "data", "openMMff",
                                    "AVAILABLE_PARAMETER.dat")):
        for line in open(f):
            p = line.split()
            if line.startswith("#") or len(p) < 11 or not p[0].isdigit():
                continue
            yes = [x == "Y" for x in p[2:11]]
            out[p[1].upper()] = {"swiss": tuple(yes[0:3]),
                                 "sannerlab": tuple(yes[3:6])}
    return out


def openmm_set_for(exe, seq):
    """The parameter set (-F) that covers every non-standard residue of `seq`
    where it sits, "sannerlab" (ADCP's default) first; None when neither does,
    and "" when the peptide has no non-standard residue."""
    res = peptide.residues(seq)
    codes = [(i, re.search(r"<([A-Z0-9]{1,4})>", r)) for i, r in enumerate(res)]
    codes = [(i, m.group(1)) for i, m in codes if m]
    if not codes:
        return ""
    table = openmm_sets(exe)
    for name in ("sannerlab", "swiss"):
        ok = True
        for i, code in codes:
            col = 1 if i == 0 else (2 if i == len(res) - 1 else 0)
            cover = (table.get(code) or {}).get(name)
            if not cover or not cover[col]:
                ok = False
                break
        if ok:
            return name
    return None


def adcp_sequence(seq, case="lower"):
    """The sequence with every letter in one case, except the library names
    inside angle brackets, which ADCP matches as written. ADCP reads a lower
    case letter as a residue that starts as coil and an upper case one as a
    residue that starts as helix."""
    out, inside = [], False
    for c in seq:
        if c == "<":
            inside = True
        out.append(c if inside else (c.upper() if case == "upper"
                                     else c.lower()))
        if c == ">":
            inside = False
    return "".join(out)


def start_arguments(seq, helical):
    """(sequence as ADCP is given it, [flags]) for a share of helical starts.

    runADCP.py reads -p as a percentage: it clamps it
    to 0 to 100, works out numHelix = replicas * p / 100, and only when
    0 < p < 100 rewrites the sequence for every replica it starts, upper
    case for replica numbers up to numHelix and lower case for the rest
    (pepSeqStr in 1.1 leaves the names in angle brackets alone). At 0 and at
    100 it does nothing, and the case of the sequence as passed decides, one
    residue at a time. So 0 and 100 are spelled out here rather than passed
    as -p: all lower case starts every replica as coil, which is what every
    run before this setting did, and all upper case starts every replica as
    a helix.
    """
    try:
        h = int(round(float(helical)))
    except (TypeError, ValueError):
        h = 0
    if h <= 0:
        return adcp_sequence(seq, "lower"), []
    if h >= 100:
        return adcp_sequence(seq, "upper"), []
    # The case is ADCP's to set, replica by replica, so the sequence goes
    # as it is written.
    return seq, ["-p", str(h)]


def helical_replicas(replicas, helical):
    """How many of `replicas` start as a helix at `helical` percent: ADCP
    starts replica numbers 1 to replicas * p / 100 as a helix (jobNum <=
    numHelix in runADCP.py), so the count is rounded down."""
    try:
        h = int(round(float(helical)))
    except (TypeError, ValueError):
        h = 0
    if h <= 0:
        return 0
    if h >= 100:
        return int(replicas)
    return int(int(replicas) * h / 100.0)


#: ADCP's flags for a cyclic peptide. runADCP.py adds a head to tail restraint for -cyc (external2=4) and a disulfide
#: potential between cysteines for -cys (SSbond=50,2.2,50,0.35).
CYCLIC_FLAGS = {"none": [], "backbone": ["-cyc"], "disulfide": ["-cys"]}


def cyclic_flags(mode):
    return list(CYCLIC_FLAGS.get(str(mode or "none"), []))


# ---- AGFR from the same environment ------------------------------------
#
# The ADCP 1.1 suite carries its own AGFR: it takes a PDB as well as a
# PDBQT, offers AutoSite 1.1 (-asv 1.1) beside AutoSite 1.0 (-asv 1.0, its
# default), and adds the receptor gradient to the maps only when asked (-g).
# It is found from the same installation as adcp, so wherever ADCP 1.1 is,
# its AGFR is too.

def agfr_script(exe):
    """runAGFR.py inside the environment `exe` belongs to, or None."""
    if not exe:
        return None
    found = glob.glob(os.path.join(home_of(exe), "lib", "python3*",
                                   "site-packages", "ADFR", "bin",
                                   "runAGFR.py"))
    return found[0] if found else None


def find_agfr(roots=()):
    """(the ADCP 1.1 executable, its runAGFR.py), or (None, None)."""
    exe = find(roots)
    script = agfr_script(exe)
    return (exe, script) if script else (None, None)


def agfr_command(exe):
    """How to start this environment's AGFR: its own Python on runAGFR.py,
    for the same reason command() does not use the bin/adcp launcher."""
    return [os.path.join(home_of(exe), "bin", "python"), agfr_script(exe)]


def agfr_env(exe):
    """The environment AGFR 1.1 runs in: the one ADCP 1.1 runs in."""
    return _env(exe)


def agfr_version(exe):
    """The version string AGFR writes into the targets it makes, read from
    Support/version.py rather than by importing it."""
    for f in glob.glob(os.path.join(home_of(exe), "lib", "python3*",
                                    "site-packages", "Support",
                                    "version.py")):
        try:
            m = re.search(r"__version__\s*=\s*['\"]([^'\"]+)", open(f).read())
        except OSError:
            continue
        if m:
            return m.group(1)
    return ""


def private_home(suite):
    """The home folder the suite's own Python is given, which is never $HOME.

    The suite's modules keep files under the home folder and look nowhere
    else: MolKit2 and AGFR make ~/.adfrsuite (with ProDy's settings in it),
    and pybel makes ~/.data/pybel and ~/.config. A cluster's home folder holds
    few files, so the suite's programs alone are given a folder of each
    person's own instead: $HOPE_CACHE/adcpsuite-home when that is set, else
    <lab folder>/.cache/<user>/adcpsuite-home, the lab folder being the one
    holding tools/docking/adcpsuite (or the folder the suite is in, when it is
    installed anywhere else). A folder in the system's temporary space is
    used when that one cannot be written. Nothing else is given this HOME:
    ssh, git and the scheduler keep the person's own.
    """
    import getpass
    import tempfile
    user = getpass.getuser()
    base = (os.environ.get("HOPE_CACHE") or "").strip()
    if not base:
        top = os.path.dirname(os.path.abspath(suite))
        if os.path.basename(top) == "docking" and \
                os.path.basename(os.path.dirname(top)) == "tools":
            top = os.path.dirname(os.path.dirname(top))
        base = os.path.join(top, ".cache", user)
    for where in (os.path.join(base, "adcpsuite-home"),
                  os.path.join(tempfile.gettempdir(), "adcpsuite-home-" + user)):
        try:
            os.makedirs(where, exist_ok=True)
        except OSError:
            continue
        if os.access(where, os.W_OK):
            return where
    return os.environ.get("HOME", "")


def _env(exe):
    # The environment's own bin first: its launcher runs whichever python
    # comes first on PATH, and only this one has ADCP's modules.
    home = home_of(exe)
    env = dict(os.environ)
    env["HOME"] = private_home(home)
    env["PATH"] = os.path.join(home, "bin") + os.pathsep + env.get("PATH", "")
    # What the bin/adcp launcher sets.
    env["ADS_ROOT"] = home
    env["PYTHONNOUSERSITE"] = "disable"
    env["LD_LIBRARY_PATH"] = (env.get("LD_LIBRARY_PATH", "") + os.pathsep
                              + os.path.join(home, "lib")).strip(os.pathsep)
    env.pop("OMP_NUM_THREADS", None)          # ADCP deadlocks if this is set
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    return env


#: The folder ADCP 1.1 searches in, inside the peptide's docking folder: -w
#: is "." and -o is "out", and runADCP.py works in <-w>/<-o>. It is deleted
#: when the docking ends unless -k is given.
WORKING = "out"


def dock_arguments(sequence, trg_name, replicas, steps, poses, cpus=0,
                   helical=0, cyclic="none", keep=False):
    """runADCP.py's arguments for one peptide, before any OpenMM flags.

    -m is the number of poses ADCP writes (its default is 100). -k keeps the
    working folder with every replica's trajectory, which the caller packs
    into one archive afterwards.
    """
    seq_arg, start = start_arguments(sequence, helical)
    args = ["-T", trg_name, "-s", seq_arg, "-N", str(replicas),
            "-n", str(steps), "-m", str(poses), "-o", WORKING, "-w", ".", "-O"]
    args += start
    args += cyclic_flags(cyclic)
    if keep:
        args += ["-k"]
    if cpus:
        args += ["-c", str(cpus)]
    if peptide.nonstandard(sequence):
        args += ["-L", "swiss"]
    return args


def dock(sequence, trg, workdir, exe, key=None, replicas=100, steps=5_000_000,
         poses=10, cpus=0, minimize=0, omm_env="implicit", log=print,
         helical=0, cyclic="none", keep=False):
    """Run ADCP 1.1 for one peptide. Returns (modes, the peptide's docking
    folder)."""
    wd = os.path.join(workdir, key or sequence)
    os.makedirs(wd, exist_ok=True)
    local_trg = os.path.join(wd, os.path.basename(trg))
    src = os.path.abspath(trg)
    if (not os.path.isfile(local_trg)
            or os.path.getsize(local_trg) != os.path.getsize(src)):
        shutil.copy(src, local_trg)
    marker = os.path.join(wd, ".adcp_finished")
    if not (os.path.isfile(marker)
            and os.path.isfile(os.path.join(wd, "out_ranked_1.pdb"))):
        # A working folder left by an attempt that was stopped part way
        # would be searched in again and its old replicas packed with the
        # new ones; and an archive packed before describes another search.
        from . import trajectories
        shutil.rmtree(os.path.join(wd, WORKING), ignore_errors=True)
        trajectories.forget(wd)
        cmd = command(exe) + dock_arguments(
            sequence, os.path.basename(trg), replicas, steps, poses, cpus=cpus,
            helical=helical, cyclic=cyclic, keep=keep)
        ffset = openmm_set_for(exe, sequence) if minimize else ""
        if minimize and ffset is None:
            # ADCP stops the whole docking, not just the minimization, when
            # no parameter set covers a residue, so the docking runs without.
            log("      %s: no OpenMM parameters cover its non-standard residues "
                "where they sit; docking without minimization" % sequence)
            minimize = 0
        base = list(cmd)
        if minimize:
            # -dr keeps ADCP's order, so out_ranked_<n> stays ADCP's rank n
            # and the OpenMM energies are reported beside it.
            cmd += ["-nmin", str(minimize), "-env", omm_env, "-dr"]
            if ffset:
                cmd += ["-F", ffset]
        # The command as it runs, in the task's log as well as at the top of
        # adcp_run.log: which flags a peptide was given, its ring among them,
        # is then on record where the task's other lines are.
        log("  adcp: %s" % " ".join(cmd[len(command(exe)):]))
        _run(cmd, wd, exe)
        n = split_poses(wd)
        if not n and minimize:
            log("      %s: the docking with OpenMM minimization gave no pose; "
                "docking again without it" % sequence)
            _run(base, wd, exe, append=True)
            n = split_poses(wd)
        elif n and minimize and not read_openmm(wd):
            # ADCP writes the poses before it minimizes them, and exits with
            # an error when OpenMM cannot use a residue's template: ADCP
            # 1.1.21 lists parameters for FVAL, FPG2, FPG3, FTR and FT6 that
            # do not match the residues it builds. The docking stands.
            log("      %s: docked, but ADCP's OpenMM minimization failed, so "
                "its poses are not minimized; see %s"
                % (sequence, os.path.join(wd, "adcp_run.log")))
        if n:
            with open(marker, "w") as fh:
                fh.write("ADCP 1.1: %d replicas, %d steps\n" % (replicas, steps))
        else:
            log("      %s: ADCP 1.1 produced no pose; see %s"
                % (sequence, os.path.join(wd, "adcp_run.log")))
    omm = read_openmm(wd)
    if omm:
        with open(os.path.join(wd, "openmm.json"), "w") as fh:
            json.dump(omm, fh, indent=1)
    return read_modes(wd), wd


def read_modes(wd):
    """[{rank, affinity, cluster}] from a peptide's docking folder.

    modes.json, written when its poses are split; where that is missing,
    ADCP's own summary, out_summary.dlg, or its copy in adcp.log. A folder
    with neither gives [].
    """
    js = os.path.join(wd, "modes.json")
    if os.path.isfile(js):
        try:
            return json.load(open(js))
        except ValueError:
            pass
    for name in ("out_summary.dlg", "adcp.log"):
        summary = os.path.join(wd, name)
        if os.path.isfile(summary):
            got = _modes11(summary)
            if got:
                return got
    return []


def _modes11(summary):
    out = []
    for line in open(summary, errors="replace"):
        if line.startswith("OMM"):
            continue
        m = _MODE11.match(line)
        if m:
            out.append({"rank": int(m.group(1)), "affinity": float(m.group(2)),
                        "cluster": int(m.group(4))})
    # The table is printed once; a second copy would restart at rank 1.
    seen, first = set(), []
    for row in out:
        if row["rank"] in seen:
            break
        seen.add(row["rank"])
        first.append(row)
    return first


def _run(cmd, wd, exe, append=False):
    with open(os.path.join(wd, "adcp_run.log"), "a" if append else "w") as fh:
        fh.write(" ".join(cmd) + "\n")
        fh.flush()
        return subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=wd,
                              env=_env(exe))


def split_poses(wd, job="out"):
    """out_out.pdb (one MODEL per pose) -> out_ranked_<n>.pdb, the summary's
    mode table -> modes.json, and the summary itself -> adcp.log, the log
    the downloads and a problem report carry."""
    multi = os.path.join(wd, "%s_out.pdb" % job)
    summary = os.path.join(wd, "%s_summary.dlg" % job)
    if os.path.isfile(summary):
        shutil.copy(summary, os.path.join(wd, "adcp.log"))
        with open(os.path.join(wd, "modes.json"), "w") as fh:
            json.dump(_modes11(summary), fh, indent=1)
    if not os.path.isfile(multi):
        return 0
    n, cur, model = 0, [], None
    for line in open(multi):
        if line.startswith("MODEL"):
            model, cur = int(line.split()[1]), []
        elif line.startswith("ENDMDL"):
            if model is not None:
                with open(os.path.join(wd, "out_ranked_%d.pdb" % model), "w") as fh:
                    fh.writelines(cur + ["END\n"])
                n += 1
            model = None
        elif model is not None:
            cur.append(line)
    return n


def read_openmm(wd):
    """The OpenMM re-ranking table, keyed by ADCP rank, or [] when -nmin was
    not used."""
    out = []
    summary = os.path.join(wd, "adcp.log")
    if not os.path.isfile(summary):
        return out
    for line in open(summary, errors="replace"):
        m = _OMM.match(line)
        if m:
            out.append({"model": int(m.group(1)), "rank_openmm": int(m.group(2)),
                        "rank_adcp": int(m.group(3)),
                        "e_complex_receptor": float(m.group(4)),
                        "e_interaction": float(m.group(5)),
                        "affinity": float(m.group(6))})
    return out
