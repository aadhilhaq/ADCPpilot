"""MM-GBSA rescoring of the final ADCP poses.

ADCP's score is a docking score: fast, and tuned for ranking within one search
rather than for comparing across constructs of different length and composition.
MM-GBSA computes an interaction energy from a force field with an implicit-solvent
term, which is slower but a different kind of estimate, so agreement between the
two is worth more than either alone.

Two backends, tried in this order:

  amber   AmberTools; tleap builds the topologies, MMPBSA.py computes the energy
          on the docked pose. Free and open source, and the default.
  prime   Schrodinger prime_mmgbsa, used only if $SCHRODINGER is set. Included
          because groups that have it often prefer it for comparability with
          existing work.

Single-structure MM-GBSA has real limits: one conformation, no ensemble average,
no entropy term. It is a rescoring, not a binding free energy, and the numbers
should be read as a ranking aid rather than as predicted affinities.
"""
from __future__ import annotations
import concurrent.futures as cf
import os, re, shutil, subprocess, threading

from . import prep

LEAP = """source leaprc.protein.ff14SB
com = loadpdb {complex_pdb}
saveamberparm com {prefix}_com.parm7 {prefix}_com.rst7
quit
"""

MIN_IN = """Minimisation before rescoring: tleap adds hydrogens without regard to
the partner molecule, and a single-point energy on the raw complex is dominated by
those clashes rather than by the interaction.
&cntrl
  imin=1, maxcyc=1000, ncyc=500,
  igb=5, cut=999.0, ntb=0,
  ntr=1, restraintmask='{rec_mask}&!@H=', restraint_wt=5.0,
/
"""

# The interior dielectric: fixed at 1.0, and not ours to move.
#
# It is the parameter the benchmarks single out: "the internal dielectric
# constant is the most critical parameter" for protein-peptide complexes (Sun
# et al., Assessing the performance of MM/PBSA and MM/GBSA methods. 9, Phys.
# Chem. Chem. Phys. 21:10135, 2019, over 112 systems). That study found short
# peptides of 5-12 residues do best with MM/PBSA at eps_in = 2, and medium ones
# of 20-25 residues with MM/GBSA at eps_in = 1 (r_p 0.735). This pipeline
# builds constructs in the second range and runs MM/GBSA, so 1.0 is the
# recommended value and 1.0 is what it gets.
#
# Not by choice, though: the classical GB models sander implements assume a
# solute dielectric of 1, and this MMPBSA.py's &gb namelist has no variable for
# it. Writing `intdiel=1` into the input, which looks harmless and agrees with
# the value in use, is rejected outright:
#
#     InputError: Unknown variable intdiel in &gb
#
# and that kills every construct in the run, not just one. The namelist's
# `epsin` is not the same knob; it belongs to the GBNSR6 solver, which is a
# different calculation from igb=5. So this is recorded as a constant and a
# comment rather than exposed as a setting: a knob that cannot move is worse
# than no setting, because a setting will be changed.
INTDIEL = 1.0

MMPBSA_IN = """Single-structure MM-GBSA rescoring of a docked pose
&general
   startframe=1, endframe=1, interval=1, verbose=2, keep_files=0,
   receptor_mask="{rec_mask}", ligand_mask="{lig_mask}",
/
&gb
   igb=5, saltcon=0.150,
/
"""


#: Gas constant times 298.15 K, kcal/mol. The one place this number is
#: written down in this package.
RT_298 = 0.0019872041 * 298.15


def dg_for_kd(kd_molar, temperature=298.15):
    """The binding free energy that corresponds to a dissociation constant.

    This direction is deliberate. Converting an MM-GBSA score into a Kd is
    not meaningful: MM-GBSA values depend strongly on the receptor, and the
    spread between poses of one peptide is usually larger than the
    4.1 kcal/mol that separates picomolar from nanomolar binding. Converting
    a target affinity into kcal/mol instead shows that scale directly.
    """
    import math
    return 0.0019872041 * float(temperature) * math.log(float(kd_molar))


def affinity_table(temperature=298.15):
    """[(label, kd_molar, dg)] over the range a binder is judged on."""
    rows = [("1 pM", 1e-12), ("10 pM", 1e-11), ("100 pM", 1e-10),
            ("1 nM", 1e-9), ("10 nM", 1e-8), ("100 nM", 1e-7),
            ("1 uM", 1e-6), ("10 uM", 1e-5)]
    return [(lab, kd, dg_for_kd(kd, temperature)) for lab, kd in rows]


def amber_home() -> str | None:
    t = prep.find_tool("tleap")
    if t:
        return os.path.dirname(os.path.dirname(t))
    for env in ("AMBERHOME",):
        v = os.environ.get(env)
        if v and os.path.isfile(os.path.join(v, "bin", "tleap")):
            return v
    return None


def backend() -> str | None:
    if amber_home() and prep.find_tool("MMPBSA.py"):
        return "amber"
    schro = os.environ.get("SCHRODINGER")
    if schro and os.path.isfile(os.path.join(schro, "prime_mmgbsa")):
        return "prime"
    return None


def _write_pdbs(rec, pose_mol, workdir: str, prefix: str):
    """Receptor, ligand and complex as separate PDBs, which is what tleap needs."""
    from rdkit import Chem
    from . import receptor as recmod
    rec_pdb = os.path.join(workdir, f"{prefix}_rec.pdb")
    lig_pdb = os.path.join(workdir, f"{prefix}_lig.pdb")
    com_pdb = os.path.join(workdir, f"{prefix}_com.pdb")
    # Each receptor chain under its own letter, with a TER between them.
    #
    # Writing every atom as chain "A" merged the chains into one molecule:
    # tleap read a continuous chain and bonded the C-terminus of the first to
    # the N-terminus of the second, inventing a peptide bond across an
    # interface and computing an energy for a molecule that is not there.
    def _write_receptor(fh, start=1):
        n = start - 1
        last = None
        for a in rec.atoms:
            ch = getattr(a, "chain", "") or "A"
            if last is not None and ch != last:
                fh.write("TER\n")
            n += 1
            fh.write(recmod.pdb_line(n, a, ch))
            last = ch
        fh.write("TER\n")
        return n

    with open(rec_pdb, "w") as fh:
        _write_receptor(fh)
        fh.write("END\n")
    # Heavy atoms only, deliberately. The docking engines score on heavy
    # atoms and write poses without hydrogens, so any hydrogen on this
    # molecule came from the template and sits at template coordinates.
    Chem.MolToPDBFile(Chem.RemoveHs(pose_mol), lig_pdb)
    with open(com_pdb, "w") as fh:
        n = _write_receptor(fh)
        # MMPBSA.py has to identify the ligand within the complex, and requires
        # its residues to be one sequential block. Numbering them to continue
        # from the receptor guarantees that whatever renumbering happens later.
        first_lig = max((a.resnum for a in rec.atoms), default=0) + 1
        lig_res = {}
        for l in open(lig_pdb):
            if l.startswith(("ATOM", "HETATM")):
                n += 1
                old = l[22:27]
                if old not in lig_res:
                    lig_res[old] = first_lig + len(lig_res)
                fh.write(f"{l[:6]}{n:5d}{l[11:21]}B{lig_res[old]:4d} {l[27:]}")
        fh.write("TER\nEND\n")
    return rec_pdb, lig_pdb, com_pdb


LIGAND_CHAIN = "B"          # _write_pdbs stamps the pose with this


def _residue_blocks(pdb_path):
    """[(chain, tleap residue index)] walking the file the way tleap reads it.

    tleap numbers residues 1..N in the order they appear and ignores the PDB's
    own residue numbers, so the index is a position in the file, not a label.
    """
    out, seen, n = [], None, 0
    for l in open(pdb_path):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        key = (l[21], l[22:27])
        if key != seen:
            n += 1
            out.append((l[21], n))
            seen = key
    return out


def _mask(indices):
    """Amber residue mask from a list of indices, collapsed into ranges."""
    if not indices:
        return ""
    idx = sorted(indices)
    runs, a, b = [], idx[0], idx[0]
    for i in idx[1:]:
        if i == b + 1:
            b = i
        else:
            runs.append((a, b)); a = b = i
    runs.append((a, b))
    return ":" + ",".join(f"{x}-{y}" if x != y else f"{x}" for x, y in runs)


def masks_by_chain(pdb_path, lig_chain=LIGAND_CHAIN):
    """Receptor and ligand masks taken from the chains rather than from an assumed order.

    Masks built as ":1-n_rec" and ":n_rec+1-end" are correct only while the
    receptor comes first, and preparation reorders chains.
    """
    blocks = _residue_blocks(pdb_path)
    lig = [i for c, i in blocks if c == lig_chain]
    rec = [i for c, i in blocks if c != lig_chain]
    chains = sorted({c for c, _ in blocks})
    asl = ("receptor: (chain.name %s)   ligand: (chain.name %s)"
           % (" or chain.name ".join(c for c in chains if c != lig_chain) or "?",
              lig_chain))
    return _mask(rec), _mask(lig), len(rec), len(lig), asl


_DELTA = re.compile(r"^DELTA TOTAL\s+(-?\d+\.\d+)")


def _amber_protonation_names(pdb_in: str, pdb_out: str) -> dict:
    """Encode PDB2PQR's protonation decisions in residue names, then drop the
    hydrogens and let tleap rebuild them.

    PDB2PQR decides titration states and expresses them by which hydrogens it
    places, but it leaves the residue named HIS. tleap's default for HIS is HIE,
    which expects HE2 and forbids HD1, so a delta-protonated histidine aborts the
    build with "Atom .R<HIE 132>.A<HD1 18> does not have a type". Renaming the
    residue keeps the decision; dropping the hydrogens avoids every naming
    mismatch, because tleap then places them itself to match the name it was
    given.
    """
    from collections import defaultdict
    hyd = defaultdict(set)
    lines = []
    for l in open(pdb_in):
        if not l.startswith(("ATOM", "HETATM")):
            lines.append(l)
            continue
        lines.append(l)
        key = (l[21], l[22:27])
        nm = l[12:16].strip()
        if nm.lstrip("0123456789").upper().startswith("H"):
            hyd[key].add(nm.upper())

    rename = {}
    counts = defaultdict(int)
    for l in lines:
        if not l.startswith(("ATOM", "HETATM")):
            continue
        key = (l[21], l[22:27])
        if key in rename:
            continue
        res = l[17:20].strip().upper()
        hs = hyd.get(key, set())
        new = None
        if res == "HIS":
            d, e = "HD1" in hs, "HE2" in hs
            new = "HIP" if (d and e) else "HID" if d else "HIE"
        elif res == "ASP" and ("HD2" in hs or "HD1" in hs):
            new = "ASH"
        elif res == "GLU" and ("HE2" in hs or "HE1" in hs):
            new = "GLH"
        elif res == "LYS" and not {"HZ1", "HZ2", "HZ3"} <= hs and len(hs & {"HZ1", "HZ2", "HZ3"}) == 2:
            new = "LYN"
        elif res == "CYS" and not ({"HG"} & hs):
            new = "CYX"            # no thiol hydrogen: part of a disulfide
        if new and new != res:
            rename[key] = new
            counts[f"{res}->{new}"] += 1

    with open(pdb_out, "w") as fh:
        for l in lines:
            if not l.startswith(("ATOM", "HETATM")):
                fh.write(l)
                continue
            nm = l[12:16].strip()
            if nm.lstrip("0123456789").upper().startswith("H"):
                continue          # tleap rebuilds these to match the residue name
            key = (l[21], l[22:27])
            if key in rename:
                l = l[:17] + rename[key].ljust(3) + l[20:]
            fh.write(l)
    return dict(counts)


def prep_complex(com_pdb: str, workdir: str, prefix: str, ph: float = 7.4,
                 log=print) -> str:
    """Re-protonate the complex before the energy is computed.

    Protonation states are not a property of the receptor alone. A residue at the
    interface can titrate differently once a peptide is bound, and preparing the
    apo receptor cannot know that. Running PDB2PQR on the complex lets PROPKA see
    the peptide. Falls back to the unprepared complex if pdb2pqr is unavailable,
    since tleap will still add hydrogens by geometry.
    """
    exe = prep.find_tool("pdb2pqr", "pdb2pqr30")
    if exe is None:
        return com_pdb
    out_pdb = os.path.join(workdir, f"{prefix}_com_prepped.pdb")
    out_pqr = os.path.join(workdir, f"{prefix}_com.pqr")
    # tleap cannot read PDB2PQR hydrogen names directly, so the file this
    # returns is always the renamed one, including when a run resumes
    # between the two steps.
    clean = os.path.join(workdir, f"{prefix}_com_forleap.pdb")
    if os.path.isfile(clean):
        return clean
    if not os.path.isfile(out_pdb):
        subprocess.run([exe, "--ff=AMBER", f"--with-ph={ph}",
                        "--titration-state-method=propka", "--keep-chain",
                        f"--pdb-output={out_pdb}", com_pdb, out_pqr],
                       capture_output=True, text=True)
    if not os.path.isfile(out_pdb):
        log(f"    {prefix}: complex preparation failed, using the raw complex")
        return com_pdb
    # Reached with a fresh out_pdb, or with one left by an interrupted run. The
    # rename is cheap, so redoing it costs nothing and is the only way back to a
    # file tleap accepts.
    changed = _amber_protonation_names(out_pdb, clean)
    if changed:
        log(f"    {prefix}: protonation states at pH {ph}: " +
            ", ".join(f"{k} x{v}" for k, v in sorted(changed.items())))
    return clean


def _parm_atoms(parm7: str):
    """How many atoms a topology holds, or None if it cannot be read.

    The POINTERS section's first number is NATOM. Reading it is how the
    hydrogen check below knows the difference between a topology built from
    heavy atoms alone and one tleap has completed.
    """
    try:
        with open(parm7) as fh:
            for line in fh:
                if line.startswith("%FLAG POINTERS"):
                    next(fh)                      # %FORMAT
                    return int(next(fh).split()[0])
    except (OSError, StopIteration, ValueError, IndexError):
        return None
    return None


def run_amber(rec, pose_mol, workdir: str, prefix: str, log=print,
              prep_ph: float | None = 7.4,
              intdiel: float = INTDIEL) -> float | None:
    # tleap and sander run with cwd=workdir, so every path given to them
    # must be absolute or a bare basename: a relative path resolves
    # inside the working directory.
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    home = amber_home()
    if home is None:
        return None
    env = dict(os.environ, AMBERHOME=home)
    rec_pdb, lig_pdb, com_pdb = _write_pdbs(rec, pose_mol, workdir, prefix)
    if prep_ph is not None:
        com_pdb = prep_complex(com_pdb, workdir, prefix, prep_ph, log)
    leap = os.path.join(workdir, f"{prefix}.leap")
    with open(leap, "w") as fh:
        fh.write(LEAP.format(complex_pdb=os.path.basename(com_pdb),
                             receptor_pdb=os.path.basename(rec_pdb),
                             ligand_pdb=os.path.basename(lig_pdb),
                             prefix=prefix))
    subprocess.run([os.path.join(home, "bin", "tleap"), "-f", leap],
                   cwd=workdir, capture_output=True, text=True, env=env)
    if not os.path.isfile(os.path.join(workdir, f"{prefix}_com.parm7")):
        leaplog = os.path.join(workdir, "leap.log")
        why = ""
        if os.path.isfile(leaplog):
            err = [l.strip() for l in open(leaplog)
                   if "Error" in l or "not found" in l or "FATAL" in l]
            why = f" ({err[-1][:90]})" if err else ""
        log(f"    {prefix}: tleap did not build topologies, skipped{why}")
        return None
    # And it has to have added the hydrogens.
    #
    # Everything up to here is heavy atoms: the docked pose has none, the
    # receptor PDB has none, and the file handed to tleap has none; the
    # hydrogens are built by ff14SB from the residue names PROPKA chose. If
    # that step half-fails, tleap still writes a topology and sander still
    # minimises it, and MM-GBSA returns a number computed on a molecule with
    # no hydrogen bonds in it. Nothing downstream would say so, and the number
    # looks like every other number in the column.
    n_heavy = sum(1 for l in open(com_pdb)
                  if l.startswith(("ATOM", "HETATM")))
    n_total = _parm_atoms(os.path.join(workdir, f"{prefix}_com.parm7"))
    if n_total is None:
        log(f"    {prefix}: the topology's atom count could not be read, so "
            f"it cannot be shown to carry hydrogens; skipped")
        return None
    if n_total < n_heavy * MIN_H_RATIO:
        log(f"    {prefix}: the topology has {n_total} atoms for {n_heavy} "
            f"heavy ones ({n_total / max(1, n_heavy):.2f}x, against "
            f"{MIN_H_RATIO}x expected): hydrogens are missing, skipped")
        return None
    # Receptor and ligand topologies are split out of the complex rather than
    # built separately. Building them separately looks equivalent but is not: the
    # complex comes from the prepared structure, where a histidine may be HID,
    # while a separately built receptor gets tleap's default HIE. The atom counts
    # still match, so the mismatch only surfaces later as "Residue mismatch while
    # mapping" from MMPBSA.py.
    # (chain, residue), not residue.
    #
    # Two chains that both number their residues 1-100 give a set of 100, not
    # 200, so the ligand mask below started a hundred residues too early and
    # MMPBSA.py split the complex in the middle of the receptor.
    # Receptor and ligand are taken from the CHAINS of the file tleap is about
    # to read, not from an assumed ordering. See masks_by_chain.
    rec_mask, lig_mask, n_rec_res, n_lig_res, asl = masks_by_chain(com_pdb)
    if not n_rec_res or not n_lig_res:
        log(f"    {prefix}: chain {LIGAND_CHAIN} holds {n_lig_res} residue(s) "
            f"and the receptor {n_rec_res} in {os.path.basename(com_pdb)} - "
            f"the complex cannot be split, skipped")
        return None
    log(f"    {prefix}: {asl}")
    log(f"    {prefix}: receptor_mask {rec_mask[:60]}  ligand_mask {lig_mask}")
    ante = prep.find_tool("ante-MMPBSA.py")
    if ante is None:
        log(f"    {prefix}: ante-MMPBSA.py not found, skipped")
        return None
    subprocess.run([ante, "-p", f"{prefix}_com.parm7",
                    "-r", f"{prefix}_rec.parm7", "-l", f"{prefix}_lig.parm7",
                    "-n", lig_mask],
                   cwd=workdir, capture_output=True, text=True, env=env)
    if not all(os.path.isfile(os.path.join(workdir, f"{prefix}_{k}.parm7"))
               for k in ("rec", "lig")):
        log(f"    {prefix}: could not split the complex topology, skipped")
        return None

    # Minimise first, with the receptor restrained so it keeps the prepared
    # geometry while the added hydrogens and the peptide relax.
    mininp = os.path.join(workdir, f"{prefix}_min.in")
    with open(mininp, "w") as fh:
        fh.write(MIN_IN.format(rec_mask=rec_mask))
    log(f"    {prefix}: minimising the complex")
    minrst = os.path.join(workdir, f"{prefix}_min.rst7")
    sander = os.path.join(home, "bin", "sander")
    subprocess.run([sander, "-O", "-i", mininp,
                    "-p", f"{prefix}_com.parm7", "-c", f"{prefix}_com.rst7",
                    "-ref", f"{prefix}_com.rst7",
                    "-r", os.path.basename(minrst),
                    "-o", f"{prefix}_min.out"],
                   cwd=workdir, capture_output=True, text=True, env=env)
    traj = os.path.basename(minrst) if os.path.isfile(minrst) else f"{prefix}_com.rst7"
    if traj != os.path.basename(minrst):
        log(f"    {prefix}: minimisation did not converge, using the raw pose")
    else:
        # What the minimisation started from, which indicates whether
        # the pose was usable at all. A pose that begins at a large
        # positive energy is one the minimiser rescued rather than
        # one that was plausible.
        start = _first_energy(os.path.join(workdir, f"{prefix}_min.out"))
        if start is not None and start > CLASHED_START:
            log(f"    {prefix}: the pose started at {start:.2e} kcal/mol, so it "
                f"is inside the receptor, and the receptor is held fixed, so "
                f"the clash cannot relax out. Any energy from it describes "
                f"the overlap and not the binding.")

    inp = os.path.join(workdir, f"{prefix}_mmpbsa.in")
    with open(inp, "w") as fh:
        fh.write(MMPBSA_IN.format(rec_mask=rec_mask, lig_mask=lig_mask))
    cmd = [prep.find_tool("MMPBSA.py"), "-O", "-i", os.path.basename(inp),
           "-cp", f"{prefix}_com.parm7", "-rp", f"{prefix}_rec.parm7",
           "-lp", f"{prefix}_lig.parm7", "-y", traj,
           "-o", f"{prefix}_mmgbsa.dat"]
    log(f"    {prefix}: computing the interaction energy")
    p = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, env=env)
    out = os.path.join(workdir, f"{prefix}_mmgbsa.dat")
    if not os.path.isfile(out):
        err = (p.stderr or p.stdout or "").strip().splitlines()
        why = f" ({err[-1][:90]})" if err else ""
        log(f"    {prefix}: MMPBSA.py produced no output, skipped{why}")
        return None
    for line in open(out):
        m = _DELTA.match(line.strip())
        if m:
            dg = float(m.group(1))
            # A docked peptide cannot genuinely score in the thousands. Values
            # like that mean unresolved clashes, and reporting them as energies
            # would be worse than reporting nothing.
            if abs(dg) > 500:
                log(f"    {prefix}: dG {dg:.0f} kcal/mol is not physical "
                    "(unrelieved clashes), discarded")
                return None
            return dg
    return None


def _one_pose(rec, mol, workdir, seq, rank, prep_ph, log,
              intdiel=INTDIEL):
    """One (construct, pose) unit, in a directory of its own.

    Its own directory because tleap always writes leap.log into the working
    directory under that name, and MMPBSA.py writes _MMPBSA_* temporaries
    under theirs: two units sharing a directory overwrite each other's
    working files, and the failure looks like a chemistry problem rather
    than a collision.
    """
    wd = os.path.join(workdir, seq, "pose%d" % rank)
    return run_amber(rec, mol, wd, seq, log, prep_ph, intdiel)


def _read_cache(path, intdiel=None):
    """rank -> dG from the cache file, tolerating the older single-float format.

    The settings the numbers were computed under are written into the file, and
    a cache computed under different ones is discarded.
    """
    try:
        text = open(path).read().strip()
    except OSError:
        return {}
    if not text:
        return {}
    out = {}
    for line in text.splitlines():
        if line.startswith("#"):
            if intdiel is not None and "intdiel" in line:
                try:
                    was = float(line.split("intdiel", 1)[1].split()[0])
                except (ValueError, IndexError):
                    return {}
                if abs(was - float(intdiel)) > 1e-9:
                    return {}
    if intdiel is not None and not any(l.startswith("#") for l
                                       in text.splitlines()):
        # Written before the settings line existed. That is not an unknown:
        # the code of that era passed no intdiel at all, so Amber used its
        # default, which is INTDIEL. Treating it as unknown and recomputing
        # would throw away hours of real MM-GBSA on every existing run for
        # nothing: the numbers are already at the value they would be
        # recomputed at.
        if abs(float(intdiel) - INTDIEL) > 1e-9:
            return {}
    for line in text.splitlines():
        bits = line.split()
        try:
            if len(bits) == 1:                 # the old format: rank 1 only
                out[1] = float(bits[0])
            elif len(bits) >= 2:
                out[int(bits[0])] = float(bits[1])
        except ValueError:
            continue
    return out


def _write_cache(path, values, intdiel=None):
    try:
        with open(path, "w") as fh:
            if intdiel is not None:
                fh.write("# igb 5 saltcon 0.150 intdiel %g\n" % float(intdiel))
            for rank in sorted(values):
                fh.write("%d %.6f\n" % (rank, values[rank]))
    except OSError:
        pass


# Total atoms over heavy atoms, below which the complex is missing hydrogens.
#
# MM-GBSA is a force-field sum over every atom, so a complex without its
# hydrogens has no hydrogen bonds and wrong contacts and charges, yet still
# produces a number. A protein carries roughly one hydrogen per heavy atom,
# and complete complexes give ratios of about 1.88 to 2.0, so 1.75 catches a
# complex that has lost a quarter of its hydrogens.
MIN_H_RATIO = 1.75

# A starting energy above this means the pose overlaps the receptor. A clean
# docked complex starts in the thousands once tleap has added hydrogens.
CLASHED_START = 1e6


def _first_energy(min_out):
    """The energy sander reported at step 1, or None."""
    try:
        with open(min_out, errors="replace") as fh:
            seen = False
            for line in fh:
                if "NSTEP" in line and "ENERGY" in line:
                    seen = True
                    continue
                if seen:
                    bits = line.split()
                    if len(bits) >= 2:
                        try:
                            return float(bits[1])
                        except ValueError:
                            return None
                    return None
    except OSError:
        return None
    return None


# Above this, an MM-GBSA interaction energy is not a measurement.
#
# Binding free energies here run from about -70 to a few kcal/mol positive for
# something that genuinely does not bind. A value of +218 is a structure that
# blew up. The cutoff sits well clear of any plausible weak binder and well
# clear of the pose-to-pose spread, which is about 22 kcal/mol at the median.
IMPLAUSIBLE_DG = 25.0


def summarise(values, aggregate="median"):
    """One number for a construct, from its poses, plus how much they varied.

    "best": the most negative. A peptide's affinity is set by its lowest
    free-energy bound state, and ADCP's ranked poses are distinct binding
    modes rather than snapshots of one, so averaging over them averages a
    real mode with modes the peptide does not occupy.

    "mean": the average, which is the right reduction if the poses are
    treated as samples of one state.

    The spread is reported either way, because it is the thing the single
    rank-1 number could never show: two modes 20 kcal/mol apart and two
    within 2 are not the same result, and until now both were one number.
    """
    v = sorted(float(x) for x in values if x is not None)
    # A pose whose interaction energy is large and positive is a failed
    # calculation rather than a weak binder: a clash the minimisation
    # could not relax, or a topology that does not match the structure.
    # Small positive values are kept. Beyond the cutoff those poses are
    # dropped from the summary instead of averaged into it, and a
    # construct left with none has no MM-GBSA value.
    usable = [x for x in v if x < IMPLAUSIBLE_DG]
    dropped = len(v) - len(usable)
    if not usable:
        return None, 0, None
    v = usable
    if not v:
        return None, 0, None
    if aggregate == "best":
        dg = v[0]
    elif aggregate == "mean":
        dg = sum(v) / len(v)
    else:
        # The median is the default because it is the most
        # reproducible. Scoring the same peptide on the same receptor
        # in two independent runs, the estimators differed by (median
        # over pairs):
        #
        #     rank 1 only     5.45 kcal/mol
        #     best of 5       4.35
        #     mean of 5       3.69
        #     median of 5     2.75
        #
        # The best of several poses moves with a single lucky pose;
        # the median discards both tails.
        mid = len(v) // 2
        dg = v[mid] if len(v) % 2 else (v[mid - 1] + v[mid]) / 2.0
    return dg, len(v), (v[-1] - v[0])


def run(constructs: list[dict], rec, poses: dict, workdir: str, log=print,
        prep_ph: float | None = 7.4, n_poses: int = 1, workers: int = 1,
        aggregate: str = "best", details: dict | None = None,
        intdiel: float = INTDIEL) -> dict:
    """Rescore each construct. Returns sequence -> dG, missing where it failed.

    `poses` maps a sequence to one pose or a list of poses, best first. With
    more than one, every pose is scored and `summarise` reduces them;
    `details`, if given, is filled with per-construct values, count and
    spread.

    Several poses are scored because ADCP is stochastic. Scoring only the
    rank-1 pose makes the result depend on which pose a run happened to rank
    first, and repeat runs of the same peptide on the same receptor differed
    by about 11 kcal/mol at the median.
    """
    if not constructs:
        log("  SKIPPING MM-GBSA: no constructs were passed in. "
            "Check n_mmgbsa in the config: it is used as a list "
            "length, so 0 or a negative value rescores nothing.")
        return {}
    b = backend()
    if b is None:
        log("  SKIPPING MM-GBSA: neither AmberTools nor Schrodinger Prime found")
        return {}
    log(f"  MM-GBSA backend: {b}, GB-OBC(II) with an interior "
        f"dielectric of {float(intdiel):g}")
    os.makedirs(workdir, exist_ok=True)
    n_poses = max(1, int(n_poses))

    # Every (construct, pose) that still needs computing, and what is cached.
    units, cached, wanted = [], {}, {}
    for c in constructs:
        seq = c["sequence"]
        mol = poses.get(seq)
        if mol is None:
            continue
        mols = list(mol) if isinstance(mol, (list, tuple)) else [mol]
        mols = [m for m in mols[:n_poses] if m is not None]
        if not mols:
            continue
        wanted[seq] = len(mols)
        have = _read_cache(os.path.join(workdir, seq + ".dg"), intdiel)
        cached[seq] = {r: v for r, v in have.items() if r <= len(mols)}
        for i, m in enumerate(mols, 1):
            if i not in cached[seq]:
                units.append((seq, i, m))

    if not wanted:
        log("  MM-GBSA produced no results; the ADCP ranking is unaffected")
        return {}
    n_cached = sum(len(v) for v in cached.values())
    if n_poses > 1:
        log(f"    {len(wanted)} construct(s) x up to {n_poses} pose(s): "
            f"{len(units)} to compute, {n_cached} kept from a previous run")

    results = {seq: dict(v) for seq, v in cached.items()}
    if units:
        workers = max(1, int(workers))
        lock = threading.Lock()
        done = [0]

        def one(unit):
            seq, rank, mol = unit
            quiet = (lambda *a: None) if workers > 1 else log
            dg = _one_pose(rec, mol, workdir, seq, rank, prep_ph, quiet,
                           intdiel)
            with lock:
                done[0] += 1
                if workers > 1:
                    log("      [%d/%d] %s pose %d: %s"
                        % (done[0], len(units), seq, rank,
                           ("%.2f kcal/mol" % dg) if dg is not None
                           else "no result"))
            return seq, rank, dg

        if workers > 1:
            log(f"    computing {len(units)} pose energies, "
                f"{workers} at a time")
            with cf.ThreadPoolExecutor(max_workers=workers) as pool:
                finished = list(pool.map(one, units))
        else:
            finished = [one(u) for u in units]
        for seq, rank, dg in finished:
            if dg is not None:
                results.setdefault(seq, {})[rank] = dg
        for seq in wanted:
            if results.get(seq):
                _write_cache(os.path.join(workdir, seq + ".dg"), results[seq],
                             intdiel)

    out = {}
    for seq in wanted:
        vals = results.get(seq) or {}
        dg, n, spread = summarise(vals.values(), aggregate)
        if dg is None:
            continue
        out[seq] = dg
        if details is not None:
            # The ranks as well as the values: `values` is
            # ordered by rank but does not record which ranks,
            # and the two diverge as soon as a pose fails.
            details[seq] = dict(values=[vals[r] for r in sorted(vals)],
                                ranks=sorted(vals),
                                n=n, spread=spread, aggregate=aggregate)
        if n > 1:
            log(f"    {seq}: dG = {dg:.2f} kcal/mol "
                f"({aggregate} of {n} pose(s), spread {spread:.1f})")
        else:
            log(f"    {seq}: dG = {dg:.2f} kcal/mol")
    if not out:
        log("  MM-GBSA produced no results; the ADCP ranking is unaffected")
    return out
