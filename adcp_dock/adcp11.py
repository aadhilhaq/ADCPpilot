"""ADCP 1.1: finding it, running it, and writing what it produces in the
layout the rest of this package reads.

ADCP 1.1 is installed by the Scripps Research Institute's own recipe as a
separate conda environment (called adcpsuite), with its own Python. It adds
non-standard and D amino acids, all-atom poses in one multi-model file, and
an optional OpenMM minimization and re-ranking of the top poses.

Everything downstream (the results table, the viewer, MM-GBSA, resume) reads
out_ranked_<n>.pdb and the mode table in adcp.log, as ADCP 1.0 writes them.
So after a 1.1 docking the multi-model file is split into those files and
the summary is copied to adcp.log; ADCP's own files are kept beside them.
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
# energy stdv, best run. ADCP 1.0 has other columns, so 1.0's reader does not
# read it.
_MODE11 = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?[\d.]+|NA)\s+(\d+)\s+")
_OMM = re.compile(r"^OMM Ranking:\s+(\d+)\s+(\d+)\s+(\d+)\s+(-?[\d.]+)\s+"
                  r"(-?[\d.]+)\s+(-?[\d.]+)")


def candidates(roots=()):
    """Where an ADCP 1.1 environment is looked for, in order."""
    out = []
    env = os.environ.get("ADCPSUITE_HOME")
    if env:
        out.append(env)
    adfr = os.environ.get("ADFRSUITE_HOME")
    inst = os.environ.get("ADCP_HOME")
    for root in (list(roots) + ([os.path.dirname(adfr)] if adfr else [])
                 + ([os.path.join(inst, "tools", "docking")] if inst else [])):
        out += [os.path.join(root, "adcpsuite"), root]
    out += [os.path.expanduser("~/.adcp-dock/adcpsuite"), "/content/adcpsuite"]
    return out


def find(roots=()):
    """The ADCP 1.1 executable, or None."""
    for home in candidates(roots):
        exe = os.path.join(home, BIN)
        site = glob.glob(os.path.join(home, "lib", "python3*", "site-packages",
                                      "ADCP", "runADCP.py"))
        if os.path.isfile(exe) and site:
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
    imports its whole modelling stack, which takes half a minute on a node
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


def adcp_sequence(seq):
    """Lower case, which starts every replica as coil, except the library
    names, which ADCP matches as written."""
    out, inside = [], False
    for c in seq:
        if c == "<":
            inside = True
        out.append(c if inside else c.lower())
        if c == ">":
            inside = False
    return "".join(out)


def _env(exe):
    # The environment's own bin first: its launcher runs whichever python
    # comes first on PATH, and only this one has ADCP's modules.
    home = home_of(exe)
    env = dict(os.environ)
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


def dock(sequence, trg, workdir, exe, key=None, replicas=100, steps=5_000_000,
         poses=10, cpus=0, minimize=0, omm_env="implicit", log=print):
    """Run ADCP 1.1 for one peptide. Returns (modes, workdir) like adcp.dock."""
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
        cmd = command(exe) + ["-T", os.path.basename(trg), "-s", adcp_sequence(sequence),
               "-N", str(replicas), "-n", str(steps), "-m", str(poses),
               "-o", "out", "-w", ".", "-O"]
        if cpus:
            cmd += ["-c", str(cpus)]
        if peptide.nonstandard(sequence):
            cmd += ["-L", "swiss"]
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
    """[{rank, affinity, cluster}] from a docking folder, 1.0 or 1.1.

    A 1.1 folder carries modes.json, written when its poses are split; any
    other folder is read by ADCP 1.0's own reader.
    """
    js = os.path.join(wd, "modes.json")
    if os.path.isfile(js):
        try:
            return json.load(open(js))
        except ValueError:
            pass
    from ._vendor import adcp as adcp10
    return adcp10.read_modes(wd)


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
    """out_out.pdb (one MODEL per pose) -> out_ranked_<n>.pdb, and the summary
    -> adcp.log, which is where the 1.0 layout keeps the mode table."""
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
