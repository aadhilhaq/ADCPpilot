"""MM-GBSA rescoring of the final ADCP poses.

By Aadhil Haq

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
import glob, os, re, shutil, subprocess, threading

from . import prep

LEAP = """source leaprc.protein.ff14SB
set default PBRadii {radii}
com = loadpdb {complex_pdb}
saveamberparm com {prefix}_com.parm7 {prefix}_com.rst7
quit
"""

# A ring peptide. The receptor and the peptide are loaded apart and combined,
# so that the peptide's residues are numbered from 1 and the closing bond can
# name them. The combined topology has the same atoms in the same order, with
# the same coordinates, as loading the complex in one piece; only the order of
# tleap's parameter tables differs. A head to tail ring is read against the
# unit names given (loadPdbUsingSeq): read with loadpdb, its first and last
# residues would get tleap's terminal templates, a free amine and a
# carboxylate, which a ring does not have.
LEAP_RING = """source leaprc.protein.ff14SB
set default PBRadii {radii}
rec = loadpdb {rec_pdb}
pep = {load_pep}
{bonds}
com = combine {{ rec pep }}
saveamberparm com {prefix}_com.parm7 {prefix}_com.rst7
quit
"""

# The peptide alone, reflected, for the hydrogens of its D residues
# (mirror_hydrogens).
LEAP_MIRROR = """source leaprc.protein.ff14SB
set default PBRadii {radii}
pep = {load_pep}
{bonds}
saveamberparm pep {prefix}_mirror.parm7 {prefix}_mirror.rst7
quit
"""

#: What a pose's molecule carries to say how its peptide is built: the
#: closing bonds of a ring and the positions of its D residues
#: (adcp_dock.mmforms.mark writes them).
RING_PROP = "mm_ring"
D_PROP = "mm_d"

#: The longest a closing bond may be after the minimisation and still count
#: as the bond: its length in ff14SB plus this. A bond the minimisation did
#: not bring to length means the energy belongs to a strained ring.
CLOSING_SLACK = 0.25

MIN_IN = """Minimisation before rescoring: tleap adds hydrogens without regard to
the partner molecule, and a single-point energy on the raw complex is dominated by
those clashes rather than by the interaction.
&cntrl
  imin=1, maxcyc=1000, ncyc=500,
  igb={igb}, cut=999.0, ntb=0,
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

# The Born radii every topology here is built with. GB-OBC(II), igb=5, was
# fitted with the mbondi2 set (Onufriev, Bashford and Case, Proteins 55:383,
# 2004), and the Amber manual pairs the two. tleap's own default is mbondi,
# and every topology built before this setting existed carries it. A run
# keeps the radii its MM-GBSA started with (radii_for), so a resumed run is
# never scored with one set for some constructs and the other for the rest.
RADII = "mbondi2"
LEGACY_RADII = "mbondi"

# The generalised Born models on offer, each with the radii it was fitted
# with. OBC(II) is the default and what every earlier run used. GBn2 (Nguyen,
# Roe and Simmerling, J Chem Theory Comput 9:2020, 2013) is Amber's igb=8 and
# only means anything with its own mbondi3 radii, which tleap also provides,
# so it runs here in AmberTools with the per residue decomposition intact.
#
# The model is read FROM the radii, never stored beside them: mbondi3 means
# GBn2 and every other set means OBC(II). A run already records its radii and
# keeps them when it is resumed (radii_for), so it keeps its model the same
# way, and one run can never hold energies from two models.
MODELS = {
    "obc2": {"igb": 5, "radii": "mbondi2", "label": "GB-OBC(II)"},
    "gbn2": {"igb": 8, "radii": "mbondi3", "label": "GBn2"},
}
DEFAULT_MODEL = "obc2"


def model_of(radii):
    """(model name, igb) for a set of radii."""
    if (radii or "").strip() == MODELS["gbn2"]["radii"]:
        return "gbn2", MODELS["gbn2"]["igb"]
    return "obc2", MODELS["obc2"]["igb"]

MMPBSA_IN = """Single-structure MM-GBSA rescoring of a docked pose
&general
   startframe=1, endframe=1, interval=1, verbose=2, keep_files=0,
   receptor_mask="{rec_mask}", ligand_mask="{lig_mask}",
/
&gb
   igb={igb}, saltcon=0.150,
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


#: What backend() calls an engine, against what gbengines calls it. Kept
#: because older configs say "prime".
BACKEND_ENGINE = {"amber": "amber", "prime": "schrodinger",
                  "openmm": "openmm"}


def backend() -> str | None:
    """The engine used when a run does not name one. None if none can run.

    Order matters: AmberTools first, because runs made before other engines
    could be chosen used it, and a run that asks for nothing else should stay
    comparable with those. OpenMM second, since it reads the same topologies.
    Prime last, because its numbers are on their own scale.

    This used to return "prime" whenever $SCHRODINGER was set, while the code
    below it ran the AmberTools path regardless. On a machine with Schrodinger
    and no AmberTools it logged Amber's solvent model and then failed inside
    tleap, which reads as a chemistry problem rather than a missing program.
    Every name it returns now has an engine behind it.
    """
    if amber_home() and prep.find_tool("MMPBSA.py"):
        return "amber"
    try:
        import openmm                                          # noqa: F401
        if amber_home():
            return "openmm"
    except Exception:                                          # noqa: BLE001
        pass
    schro = os.environ.get("SCHRODINGER")
    if schro and os.path.isfile(os.path.join(schro, "prime_mmgbsa")):
        return "prime"
    return None


def engine_runner(engine: str):
    """The function that computes one pose's energy for this engine."""
    return {"amber": run_amber, "openmm": run_openmm,
            "schrodinger": run_prime}.get(engine, run_amber)


def cache_name(seq: str, engine: str = "amber") -> str:
    """The cache file for this construct under this engine.

    AmberTools keeps the bare name it has always had, so every cache on disk
    from every previous run is still read. A second engine writing the same
    file would have its numbers read back as AmberTools': different force
    field, different scale, and nothing in the file to tell them apart. The
    others are suffixed.
    """
    return seq + (".dg" if engine == "amber" else ".%s.dg" % engine)


def receptor_chains(rec):
    """The chain letters the receptor is written under: each atom's own, and
    "A" for one with none, as _write_pdbs writes them."""
    return {(getattr(a, "chain", "") or "A") for a in rec.atoms}


def _write_pdbs(rec, pose_mol, workdir: str, prefix: str, lig_chain=None):
    """Receptor, ligand and complex as separate PDBs, which is what tleap needs.

    The peptide is written under `lig_chain`, by default the letter
    peptide_chain gives for this receptor: one the receptor does not use.
    """
    from rdkit import Chem
    from . import receptor as recmod
    rec_pdb = os.path.join(workdir, f"{prefix}_rec.pdb")
    lig_pdb = os.path.join(workdir, f"{prefix}_lig.pdb")
    com_pdb = os.path.join(workdir, f"{prefix}_com.pdb")
    if lig_chain is None:
        lig_chain = peptide_chain(receptor_chains(rec))
    # Each receptor chain under its own letter, with a TER between them.
    #
    # Writing every atom as chain "A" merged the chains into one molecule:
    # tleap read a continuous chain and bonded the C-terminus of the first to
    # the N-terminus of the second, inventing a peptide bond across an
    # interface and computing an energy for a molecule that is not there.
    # separate_chains then puts a TER wherever one chain is broken as well.
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
                fh.write(f"{l[:6]}{n:5d}{l[11:21]}{lig_chain}"
                         f"{lig_res[old]:4d} {l[27:]}")
        fh.write("TER\nEND\n")
    for path in (rec_pdb, com_pdb):
        separate_chains(path, lig_chain)
    return rec_pdb, lig_pdb, com_pdb


#: The chain the peptide is written under whenever the receptor does not use
#: it, which is every receptor every run before this could choose one had.
LIGAND_CHAIN = "B"
# The letters tried after it, in order, for a receptor that uses B.
_SPARE_CHAINS = ("BCDEFGHIJKLMNOPQRSTUVWXYZ" "abcdefghijklmnopqrstuvwxyz"
                 "0123456789")


def peptide_chain(used):
    """The chain letter the peptide is written under in the complex: B, as
    it always was, unless the receptor uses B, and then the first letter the
    receptor does not use.

    The receptor's chain and the peptide's must differ. Everything that
    tells the two apart reads the letter: the split of the complex into
    receptor and ligand for MMPBSA.py, the reduce step that protonates the
    peptide, the ring and D residue builds, and the interface panel. A
    receptor on chain B with the peptide stamped B as well was one molecule
    to all of them: tleap bonded the receptor's last residue to the
    peptide's first across 40 A and stopped.
    """
    used = {str(c) for c in used}
    for c in _SPARE_CHAINS:
        if c not in used:
            return c
    raise ValueError("the receptor uses every chain letter a PDB file can "
                     "hold, so the peptide cannot be given one of its own")


#: Two residues of one chain whose C and N are further apart than this are
#: not bonded: a peptide bond is 1.33 A, and anything past 2 A is a gap, an
#: unmodelled loop or the end of one piece of a chain and the start of
#: another.
BREAK_A = 2.0


def separate_chains(pdb_path, lig_chain=LIGAND_CHAIN, log=None):
    """Write the TER records of a complex so tleap joins what is bonded and
    nothing else, in place, and return the breaks found inside the receptor.

    A TER goes between chains, between receptor and peptide, and wherever a
    receptor chain is broken (its C and the next residue's N more than BREAK_A
    apart): tleap otherwise reads one continuous chain and bonds across the
    gap, a bond of tens of angstroms that the minimisation then pulls on, or
    that stops the build. The peptide's own residues are never split. Every
    TER the file had is replaced, so what reduce or PDB2PQR did with them
    does not matter, and running it again changes nothing.
    """
    import numpy as np
    lines = open(pdb_path).readlines()
    head, atoms = [], []
    for l in lines:
        if _is_atom(l):
            atoms.append(l if l.endswith("\n") else l + "\n")
        elif not atoms and not l.startswith(("TER", "END", "MODEL")):
            head.append(l)
    residues = []
    for l in atoms:
        key = (l[21], l[17:27])
        if not residues or residues[-1]["key"] != key:
            residues.append({"key": key, "chain": l[21], "lines": [],
                             "N": None, "C": None, "label": l[22:27].strip()})
        r = residues[-1]
        r["lines"].append(l)
        nm = l[12:16].strip()
        if nm in ("N", "C"):
            try:
                r[nm] = np.array([float(l[30:38]), float(l[38:46]),
                                  float(l[46:54])])
            except ValueError:
                pass
    out, breaks = list(head), []
    for i, r in enumerate(residues):
        if i:
            prev = residues[i - 1]
            if prev["chain"] != r["chain"]:
                out.append("TER\n")
            elif (r["chain"] != lig_chain and prev["C"] is not None
                  and r["N"] is not None
                  and float(np.linalg.norm(r["N"] - prev["C"])) > BREAK_A):
                out.append("TER\n")
                breaks.append((r["chain"], prev["label"], r["label"]))
        out.extend(r["lines"])
    out += ["TER\n", "END\n"]
    with open(pdb_path, "w") as fh:
        fh.writelines(out)
    if log and breaks:
        log("    the receptor is broken in %d place(s) (%s), and each piece "
            "is built with ends of its own rather than bonded across the gap"
            % (len(breaks), ", ".join("%s:%s to %s" % b for b in breaks[:4])
                + (" and more" if len(breaks) > 4 else "")))
    return breaks


def peptide_chain_in(lines_or_path):
    """The chain of the last atom of a complex: the peptide's, since
    _write_pdbs writes it last. For a caller that does not say."""
    lines = (open(lines_or_path).readlines()
             if isinstance(lines_or_path, str) else lines_or_path)
    for l in reversed(lines):
        if l.startswith(("ATOM", "HETATM")):
            return l[21]
    return LIGAND_CHAIN


def peptide_chain_of(posedir, prefix):
    """The chain the peptide was written under for one scored pose, read
    from what _topologies recorded beside it; B for a pose scored before the
    chain could be anything else."""
    import json
    try:
        got = json.load(open(os.path.join(posedir,
                                          f"{prefix}_chains.json")))
        return str(got.get("peptide") or LIGAND_CHAIN)
    except (OSError, ValueError, AttributeError):
        return LIGAND_CHAIN


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


def _amber_name(res: str, hs: set) -> str | None:
    """The Amber residue name the hydrogens on a residue call for, or None.

    A histidine is named for its ring protons, a carboxylate with a proton on
    it for its neutral form, a lysine with two amine hydrogens LYN, and a
    cysteine without its thiol hydrogen CYX, since that is a disulphide.
    """
    if res == "HIS":
        d, e = "HD1" in hs, "HE2" in hs
        return "HIP" if (d and e) else "HID" if d else "HIE"
    if res == "ASP" and ("HD2" in hs or "HD1" in hs):
        return "ASH"
    if res == "GLU" and ("HE2" in hs or "HE1" in hs):
        return "GLH"
    if res == "LYS" and not {"HZ1", "HZ2", "HZ3"} <= hs \
            and len(hs & {"HZ1", "HZ2", "HZ3"}) == 2:
        return "LYN"
    if res == "CYS" and not ({"HG"} & hs):
        return "CYX"            # no thiol hydrogen: part of a disulfide
    return None


def _is_h(name: str) -> bool:
    return name.lstrip("0123456789").upper().startswith("H")


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
        if _is_h(nm):
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
        new = _amber_name(res, hyd.get(key, set()))
        if new and new != res:
            rename[key] = new
            counts[f"{res}->{new}"] += 1

    with open(pdb_out, "w") as fh:
        for l in lines:
            if not l.startswith(("ATOM", "HETATM")):
                fh.write(l)
                continue
            nm = l[12:16].strip()
            if _is_h(nm):
                continue          # tleap rebuilds these to match the residue name
            key = (l[21], l[22:27])
            if key in rename:
                l = l[:17] + rename[key].ljust(3) + l[20:]
            fh.write(l)
    return dict(counts)


def residue_states(pdb_path: str) -> dict:
    """{(chain, number): Amber name} for the residues of a protonated PDB
    whose hydrogens call for another name (HIS to HID, HIE or HIP, CYS to
    CYX, and so on: _amber_name). Keyed as the complex is written: a blank
    chain is "_", and the number without its insertion code. Empty for a
    file with no hydrogens, which says nothing about any residue's state."""
    from collections import defaultdict
    hyd, resn, n_h = defaultdict(set), {}, 0
    for l in open(pdb_path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        try:
            key = (l[21] if l[21] != " " else "_", int(l[22:26]))
        except ValueError:
            continue
        resn.setdefault(key, l[17:20].strip().upper())
        nm = l[12:16].strip()
        if _is_h(nm):
            hyd[key].add(nm.upper())
            n_h += 1
    if not n_h:
        return {}
    out = {}
    for key, res in resn.items():
        new = _amber_name(res, hyd.get(key, set()))
        if new and new != res:
            out[key] = new
    return out


def prep_complex_reduce(com_pdb: str, workdir: str, prefix: str,
                        protonation: dict, log=print,
                        lig_chain: str | None = None) -> str:
    """The complex for tleap when the receptor was protonated by reduce.

    The receptor's residues are named for the protonation the docked
    receptor has (protonation["states"], read from its hydrogens by
    residue_states): HID, HIE or HIP from where reduce put the ring protons,
    CYX in a disulphide. Asn, Gln and His that reduce turned over are turned
    over in the receptor's coordinates, which are the prepared ones, so they
    stay turned. The peptide is protonated by reduce on the whole complex,
    so its histidines are chosen with the receptor around them, and its
    coordinates are taken from reduce's output, with any amide it turned.
    Every hydrogen is then dropped and tleap builds them to match the
    names, as for PDB2PQR.

    The peptide is the chain `lig_chain`, which the receptor never uses
    (peptide_chain): reduce's output is read back by it, and a receptor on
    the same letter would be read back as part of the peptide.
    """
    clean = os.path.join(workdir, f"{prefix}_com_forleap.pdb")
    if os.path.isfile(clean):
        return clean
    if lig_chain is None:
        lig_chain = peptide_chain_in(com_pdb)
    states = dict(protonation.get("states") or {})
    exe = protonation.get("reduce")
    lig_lines, lig_states = None, {}
    if exe and os.path.isfile(exe):
        firsts = {}
        for l in open(com_pdb):
            if l.startswith(("ATOM", "HETATM")) and l[21] not in firsts:
                try:
                    firsts[l[21]] = int(l[22:26])
                except ValueError:
                    pass
        cmd = [exe, "-BUILD", "-Nterm%d" % max([1] + list(firsts.values()))]
        het = exe + "_wwPDB_het_dict.txt"
        if os.path.isfile(het):
            cmd += ["-DB", het]
        cmd.append(com_pdb)
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            text = p.stdout or ""
        except (OSError, subprocess.SubprocessError):
            text = ""
        if "USER  MOD reduce" in text[:4000]:
            red = os.path.join(workdir, f"{prefix}_com_reduce.pdb")
            with open(red, "w") as fh:
                fh.write(text)
            lig_states = {k: v for k, v in residue_states(red).items()
                          if k[0] == lig_chain}
            lig_lines = [l + "\n" for l in text.splitlines()
                         if l.startswith(("ATOM", "HETATM"))
                         and l[21] == lig_chain
                         and not _is_h(l[12:16].strip())]
        else:
            log(f"    {prefix}: reduce could not protonate the complex; the "
                f"peptide keeps tleap's standard states")
    else:
        log(f"    {prefix}: reduce was not found; the peptide keeps tleap's "
            f"standard states")

    def named(l, table):
        try:
            key = (l[21], int(l[22:26]))
        except ValueError:
            return l
        new = table.get(key)
        if new and new != l[17:20].strip():
            return l[:17] + new.ljust(3) + l[20:]
        return l

    out, done = [], False
    for l in open(com_pdb):
        if not l.startswith(("ATOM", "HETATM")):
            out.append(l)
            continue
        if _is_h(l[12:16].strip()):
            continue
        if l[21] == lig_chain:
            if lig_lines is not None:
                if not done:
                    out.extend(named(x, lig_states) for x in lig_lines)
                    done = True
                continue
            out.append(named(l, lig_states))
            continue
        out.append(named(l, states))
    with open(clean, "w") as fh:
        fh.writelines(out)
    per_res = {}
    for l in out:
        if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() in (
                "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN"):
            per_res[(l[21], l[22:27])] = l[17:20].strip()
    tally = {}
    for v in per_res.values():
        tally[v] = tally.get(v, 0) + 1
    if tally:
        log(f"    {prefix}: protonation as reduce placed it: " +
            ", ".join(f"{k} x{v}" for k, v in sorted(tally.items())))
    return clean


def _prepared_for(com_pdb, workdir, prefix, how):
    """Whether the complex already prepared for tleap in `workdir` was made
    from this complex, prepared this way; (fresh, stamp).

    A prepared complex is kept so that a resumed run does not protonate it
    again. One made from a different complex, or by code that wrote the
    complex differently, must not stand in for this one: a receptor on chain
    B was once written with the peptide on B as well, and its prepared
    complex, kept, would stop every rescoring of it the same way. So what it
    was made from is recorded beside it, and anything prepared from
    something else is removed and made again.
    """
    import hashlib
    import json
    stamp = {"complex": hashlib.md5(open(com_pdb, "rb").read()).hexdigest(),
             "how": how}
    rec = os.path.join(workdir, f"{prefix}_com_forleap.json")
    clean = os.path.join(workdir, f"{prefix}_com_forleap.pdb")
    try:
        old = json.load(open(rec))
    except (OSError, ValueError):
        old = None
    if old == stamp and os.path.isfile(clean):
        return True, stamp
    for tail in ("_com_forleap.pdb", "_com_forleap.json", "_com_prepped.pdb",
                 "_com.pqr", "_com_reduce.pdb"):
        stale = os.path.join(workdir, prefix + tail)
        if os.path.isfile(stale):
            os.remove(stale)
    return False, stamp


def prep_complex(com_pdb: str, workdir: str, prefix: str, ph: float = 7.4,
                 log=print, protonation: dict | None = None,
                 lig_chain: str | None = None) -> str:
    """Re-protonate the complex before the energy is computed.

    Protonation states are not a property of the receptor alone. A residue at the
    interface can titrate differently once a peptide is bound, and preparing the
    apo receptor cannot know that. Running PDB2PQR on the complex lets PROPKA see
    the peptide. Falls back to the unprepared complex if pdb2pqr is unavailable,
    since tleap will still add hydrogens by geometry.

    A receptor protonated by reduce carries its states in `protonation`
    (pipeline._carry_protonation), and the complex is then named from them
    and from reduce instead: prep_complex_reduce.

    Either way, the file returned has its TER records written by
    separate_chains: between chains, between receptor and peptide (chain
    `lig_chain`), and at every break in a receptor chain.
    """
    import json
    if lig_chain is None:
        lig_chain = peptide_chain_in(com_pdb)
    reduce_ = bool(protonation and protonation.get("method") == "reduce")
    how = ({"method": "reduce",
            "states": sorted("%s:%s=%s" % (k[0], k[1], v) for k, v in
                             (protonation.get("states") or {}).items())}
           if reduce_ else {"method": "propka", "ph": float(ph)})
    how["peptide_chain"] = lig_chain
    fresh, stamp = _prepared_for(com_pdb, workdir, prefix, how)
    clean = os.path.join(workdir, f"{prefix}_com_forleap.pdb")
    if fresh:
        return clean
    if reduce_:
        got = prep_complex_reduce(com_pdb, workdir, prefix, protonation, log,
                                  lig_chain=lig_chain)
    else:
        got = _prep_complex_propka(com_pdb, workdir, prefix, ph, log)
    separate_chains(got, lig_chain, log=lambda m: log(f"    {prefix}: "
                                                      + m.strip()))
    if os.path.abspath(got) == os.path.abspath(clean):
        with open(os.path.join(workdir, f"{prefix}_com_forleap.json"),
                  "w") as fh:
            json.dump(stamp, fh)
    return got


def _prep_complex_propka(com_pdb, workdir, prefix, ph, log):
    """The complex protonated by PDB2PQR with PROPKA at `ph`, renamed for
    tleap (_amber_protonation_names); the raw complex when PDB2PQR is not
    here or fails."""
    exe = prep.find_tool("pdb2pqr", "pdb2pqr30")
    if exe is None:
        return com_pdb
    out_pdb = os.path.join(workdir, f"{prefix}_com_prepped.pdb")
    out_pqr = os.path.join(workdir, f"{prefix}_com.pqr")
    # tleap cannot read PDB2PQR hydrogen names directly, so the file this
    # returns is always the renamed one.
    clean = os.path.join(workdir, f"{prefix}_com_forleap.pdb")
    if not os.path.isfile(out_pdb):
        subprocess.run([exe, "--ff=AMBER", f"--with-ph={ph}",
                        "--titration-state-method=propka", "--keep-chain",
                        f"--pdb-output={out_pdb}", com_pdb, out_pqr],
                       capture_output=True, text=True)
    if not os.path.isfile(out_pdb):
        log(f"    {prefix}: complex preparation failed, using the raw complex")
        return com_pdb
    # The rename is cheap, so redoing it costs nothing and is the only way
    # back to a file tleap accepts.
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


# ---- peptides beyond a linear chain of L residues -------------------------

def peptide_spec(mol):
    """How this pose's peptide is to be built: {"ring": None or {"kind",
    "pairs", "gap"}, "d": [1-based positions of D residues]}. A linear
    peptide of L residues carries neither, and is built as it always was."""
    import json
    ring, d = None, []
    if mol is not None and mol.HasProp(RING_PROP):
        ring = json.loads(mol.GetProp(RING_PROP))
    if mol is not None and mol.HasProp(D_PROP):
        d = [int(x) for x in mol.GetProp(D_PROP).split(",") if x.strip()]
    return {"ring": ring, "d": d}


def _is_atom(line):
    return line.startswith(("ATOM", "HETATM"))


def _peptide_lines(lines, renumber=True, mirror=False, lig_chain=None):
    """The peptide's atom lines, numbered from 1 when asked and reflected
    through the yz plane when asked, with a TER and an END."""
    if lig_chain is None:
        lig_chain = peptide_chain_in(lines)
    order, out = {}, []
    for l in lines:
        if not _is_atom(l) or l[21] != lig_chain:
            continue
        key = l[22:27]
        if key not in order:
            order[key] = len(order) + 1
        if renumber:
            l = l[:22] + "%4d " % order[key] + l[27:]
        if mirror:
            l = l[:30] + "%8.3f" % (-float(l[30:38])) + l[38:]
        out.append(l if l.endswith("\n") else l + "\n")
    return out + ["TER\n", "END\n"]


def _units(pep_lines):
    """The residue names of a peptide file, in order."""
    out, last = [], None
    for l in pep_lines:
        if _is_atom(l) and l[22:27] != last:
            out.append(l[17:20].strip())
            last = l[22:27]
    return out


def _ring_ready(com_pdb, ring, lig_chain=None):
    """Make the complex file what a ring is built from, in place, and return
    the peptide residue names in order.

    Head to tail: no OXT, which the protonation step may have added to what it
    took for a C terminus. Disulfide: the joined cysteines named CYX, and any
    other cysteine the protonation named CYX named CYS again, since only the
    joined ones lose their thiol hydrogen. Running it twice changes nothing,
    so a file kept from an interrupted run is fine.
    """
    lines = open(com_pdb).readlines()
    if lig_chain is None:
        lig_chain = peptide_chain_in(lines)
    order = {}
    for l in lines:
        if _is_atom(l) and l[21] == lig_chain and l[22:27] not in order:
            order[l[22:27]] = len(order) + 1
    n = len(order)
    pairs = [(int(a), int(b)) for a, b in ring.get("pairs") or []]
    if not pairs:
        raise ValueError("a ring with no closing bond")
    if ring["kind"] == "backbone" and pairs != [(1, n)]:
        raise ValueError("a head to tail ring joins residue 1 to %d, not %s"
                         % (n, pairs))
    joined = {p for pr in pairs for p in pr}
    out = []
    for l in lines:
        if _is_atom(l) and l[21] == lig_chain:
            pos = order[l[22:27]]
            res = l[17:20].strip()
            if ring["kind"] == "backbone" and l[12:16].strip() == "OXT":
                continue
            if ring["kind"] == "disulfide":
                if pos in joined:
                    if res not in ("CYS", "CYX", "CYM"):
                        raise ValueError("residue %d is %s, not a cysteine"
                                         % (pos, res))
                    l = l[:17] + "CYX" + l[20:]
                elif res == "CYX":
                    l = l[:17] + "CYS" + l[20:]
        out.append(l)
    with open(com_pdb, "w") as fh:
        fh.writelines(out)
    return _units([l for l in out if _is_atom(l) and l[21] == lig_chain])


def _pep_leap(ring, pep_file, units):
    """(how tleap loads the peptide, the bond commands) for a ring or a
    linear peptide whose residues are numbered from 1."""
    if not ring:
        return "loadpdb %s" % pep_file, ""
    pairs = [(int(a), int(b)) for a, b in ring["pairs"]]
    if ring["kind"] == "backbone":
        return ("loadPdbUsingSeq %s { %s }" % (pep_file, " ".join(units)),
                "bond pep.1.N pep.%d.C" % len(units))
    return ("loadpdb %s" % pep_file,
            "\n".join("bond pep.%d.SG pep.%d.SG" % pr for pr in pairs))


def _ring_leap(com_pdb, ring, workdir, prefix, radii, lig_chain=None):
    """The tleap input for a ring peptide, and the two files it reads."""
    if lig_chain is None:
        lig_chain = peptide_chain_in(com_pdb)
    units = _ring_ready(com_pdb, ring, lig_chain)
    lines = open(com_pdb).readlines()
    rec, last_ter = [], False
    for l in lines:
        if _is_atom(l) and l[21] != lig_chain:
            rec.append(l)
            last_ter = False
        elif l.startswith("TER") and rec and not last_ter:
            rec.append("TER\n")
            last_ter = True
    rec_file, pep_file = f"{prefix}_rec_forleap.pdb", f"{prefix}_pep_forleap.pdb"
    with open(os.path.join(workdir, rec_file), "w") as fh:
        fh.writelines(rec + ([] if last_ter else ["TER\n"]) + ["END\n"])
    with open(os.path.join(workdir, pep_file), "w") as fh:
        fh.writelines(_peptide_lines(lines, lig_chain=lig_chain))
    load, bonds = _pep_leap(ring, pep_file, units)
    return LEAP_RING.format(radii=radii, rec_pdb=rec_file, load_pep=load,
                            bonds=bonds, prefix=prefix)


def _is_h_name(name):
    return name.lstrip("0123456789").upper().startswith("H")


def mirror_hydrogens(workdir, prefix, com_pdb, spec, n_rec_res, home, env,
                     radii=RADII, log=print, lig_chain=None):
    """Put the atoms tleap adds to a D residue where a D residue has them.

    tleap adds every hydrogen from its residue template's own geometry, so
    on a D residue built with the L template it puts the alpha hydrogen, and
    each hydrogen whose place depends on handedness, on the wrong side: on
    the 13 residue p53 peptide reflected, 81 of 108 hydrogens landed 0.1 to
    1.8 A from their mirror images, and the angle energy rose from 1065 to
    3113 kcal/mol. A minimisation that started there would push the heavy
    atoms, which are free, towards the L form.

    So the peptide is reflected, which makes each D residue an L one, built
    alone by tleap the same way, and reflected back: each D residue's
    hydrogens are then the mirror image of the hydrogens tleap gives its L
    form, which is where they belong. The same goes for the one heavy atom
    tleap adds, the OXT of a C terminal residue, which ADCP does not write:
    its place is taken from the template too. L residues keep tleap's own.
    Returns a short report, or None when the two builds do not match atom
    for atom.
    """
    import numpy as np
    from . import amberfiles as af
    lines = open(com_pdb).readlines()
    if lig_chain is None:
        lig_chain = peptide_chain_in(lines)
    mfile = f"{prefix}_mirror_pep.pdb"
    pep = _peptide_lines(lines, mirror=True, lig_chain=lig_chain)
    with open(os.path.join(workdir, mfile), "w") as fh:
        fh.writelines(pep)
    load, bonds = _pep_leap(spec.get("ring"), mfile, _units(pep))
    leapf = os.path.join(workdir, f"{prefix}_mirror.leap")
    with open(leapf, "w") as fh:
        fh.write(LEAP_MIRROR.format(radii=radii, load_pep=load, bonds=bonds,
                                    prefix=prefix))
    # tleap always writes leap.log where it runs; the complex's is kept.
    main_log = os.path.join(workdir, "leap.log")
    kept = main_log + ".complex"
    if os.path.isfile(main_log):
        shutil.copy(main_log, kept)
    subprocess.run([os.path.join(home, "bin", "tleap"), "-f", leapf],
                   cwd=workdir, capture_output=True, text=True, env=env)
    if os.path.isfile(main_log):
        shutil.move(main_log, os.path.join(workdir, f"{prefix}_mirror_leap.log"))
    if os.path.isfile(kept):
        shutil.move(kept, main_log)
    mparm = os.path.join(workdir, f"{prefix}_mirror.parm7")
    mrst = os.path.join(workdir, f"{prefix}_mirror.rst7")
    if not (os.path.isfile(mparm) and os.path.isfile(mrst)):
        log(f"    {prefix}: tleap could not build the reflected peptide, so "
            f"its D residues' hydrogens cannot be placed; skipped")
        return None
    cparm = af.read_parm7(os.path.join(workdir, f"{prefix}_com.parm7"))
    crst = os.path.join(workdir, f"{prefix}_com.rst7")
    xc = af.read_rst7(crst)
    xm = af.read_rst7(mrst)
    lig = af.residues(cparm)[n_rec_res:]
    mres = af.residues(af.read_parm7(mparm))
    if len(lig) != len(mres):
        log(f"    {prefix}: the reflected peptide has {len(mres)} residues "
            f"against {len(lig)}; skipped")
        return None
    # The atoms each peptide residue came with, by position: the rest were
    # added by tleap.
    docked, order = {}, {}
    for l in lines:
        if _is_atom(l) and l[21] == lig_chain:
            pos = order.setdefault(l[22:27], len(order) + 1)
            docked.setdefault(pos, set()).add(l[12:16].strip())
    moved = 0
    for p in spec["d"]:
        if not 0 < p <= len(lig):
            log(f"    {prefix}: D residue {p} is not in the peptide; skipped")
            return None
        _, lab_c, at_c = lig[p - 1]
        _, lab_m, at_m = mres[p - 1]
        if lab_c != lab_m or set(at_c) != set(at_m):
            log(f"    {prefix}: residue {p} is {lab_c} in the complex and "
                f"{lab_m} reflected, atom for atom different; skipped")
            return None
        # Every atom tleap added: the hydrogens, and a C terminal OXT, which
        # ADCP does not write and tleap places from the template as well.
        given = docked.get(p, set())
        for nm, ic in at_c.items():
            if _is_h_name(nm) or nm not in given:
                x = xm[at_m[nm]].copy()
                x[0] = -x[0]
                xc[ic] = x
                moved += 1
    # A centre the mirror image does not cover. ADCP 1.1 builds &T and &I
    # with the alpha carbon inverted and the beta carbon as in the L residue
    # (the allo forms), so the reflected residue is not the L template's
    # either, and its beta hydrogen lands beside its neighbours. Any
    # hydrogen on a stereocentre of a D residue that is still on the side of
    # its three heavy neighbours is put where a tetrahedral centre has it:
    # 1.09 A out, opposite the mean of the three bond directions.
    turned = 0
    first = n_rec_res + 1
    centres = af.stereo_centres(cparm, [first + p - 1 for p in spec["d"]])
    for c, h in zip(centres, af.handedness(centres, xc)):
        ic, ia, ib, iq, ih = c[3]
        if ih is None or h["h_opposite"]:
            continue
        u = sum((xc[i] - xc[ic]) / np.linalg.norm(xc[i] - xc[ic])
                for i in (ia, ib, iq))
        xc[ih] = xc[ic] - 1.09 * u / np.linalg.norm(u)
        turned += 1
    tleap_copy = os.path.join(workdir, f"{prefix}_com_tleap.rst7")
    if not os.path.isfile(tleap_copy):
        shutil.copy(crst, tleap_copy)
    af.write_rst7(crst, xc, like=tleap_copy)
    return {"d_residues": list(spec["d"]), "atoms_placed": moved,
            "hydrogens_turned": turned}


def check_peptide(topo, traj, log=print):
    """What the minimisation did to the peptide's stereocentres and to its
    ring, from the topology and the coordinates before and after.

    Returns (why not, report). "why not" is None when the pose may be
    scored, and a reason otherwise:

      a residue whose alpha carbon does not have the handedness the sequence
      gives it in the docked pose, or loses it in the minimisation, or a
      threonine or isoleucine beta carbon that changes, or a hydrogen on a
      stereocentre left on the same side as its three heavy neighbours;

      a closing bond that is not in the topology, or that the minimisation
      left longer than its length in the force field plus CLOSING_SLACK.

    The report is written beside the pose as <prefix>_checks.json, with the
    closing bond's length and its energy term before and after, so the run
    can say what the minimisation did to each ring.
    """
    import json
    from . import amberfiles as af
    spec = getattr(topo, "spec", None) or {"ring": None, "d": []}
    parm = af.read_parm7(topo.path("com", "parm7"))
    before = af.read_rst7(topo.path("com", "rst7"))
    tp = traj if os.path.isabs(traj) else os.path.join(topo.workdir, traj)
    after = af.read_rst7(tp) if os.path.isfile(tp) else before
    first = topo.n_rec_res + 1
    lig = list(range(first, first + topo.n_lig_res))
    centres = af.stereo_centres(parm, lig)
    hb, ha = af.handedness(centres, before), af.handedness(centres, after)
    d = set(spec.get("d") or [])
    problems, stereo = [], []
    for b, a in zip(hb, ha):
        pos = b["residue"] - first + 1
        want = None
        if b["centre"] == "CA":
            want = -1 if pos in d else 1
        row = {"position": pos, "residue": b["name"], "centre": b["centre"],
               "before": b["volume"], "after": a["volume"],
               "h_opposite": a["h_opposite"]}
        stereo.append(row)
        form = "D" if pos in d else "L"
        if want is not None and b["volume"] * want <= 0:
            problems.append("residue %d (%s) is not %s in the docked pose"
                            % (pos, b["name"], form))
        elif b["volume"] * a["volume"] <= 0:
            problems.append("residue %d (%s) changed handedness at %s in the "
                            "minimisation" % (pos, b["name"], b["centre"]))
        elif a["h_opposite"] is False:
            problems.append("residue %d (%s) has its %s hydrogen on the wrong "
                            "side after the minimisation"
                            % (pos, b["name"], b["centre"]))
    closing = []
    ring = spec.get("ring")
    if ring:
        res = af.residues(parm)
        for a_pos, b_pos in ring["pairs"]:
            ra, rb = res[first + int(a_pos) - 2], res[first + int(b_pos) - 2]
            if ring["kind"] == "backbone":
                ia, ib, what = ra[2].get("N"), rb[2].get("C"), "C to N"
            else:
                ia, ib, what = ra[2].get("SG"), rb[2].get("SG"), "S to S"
            if ia is None or ib is None:
                problems.append("the closing bond's atoms are not in the "
                                "topology")
                continue
            t0 = af.bond_term(parm, before, ia, ib)
            t1 = af.bond_term(parm, after, ia, ib)
            closing.append({"bond": what, "residues": [int(a_pos), int(b_pos)],
                            "docked_gap": ring.get("gap"),
                            "r0": t0["r0"], "k": t0["k"],
                            "before": t0["length"], "after": t1["length"],
                            "energy_before": t0["energy"],
                            "energy_after": t1["energy"]})
            if not t0["present"]:
                problems.append("the closing %s bond is not in the topology"
                                % what)
            elif t1["length"] > t0["r0"] + CLOSING_SLACK:
                problems.append("the closing %s bond is %.2f A after the "
                                "minimisation, against %.2f A for the bond"
                                % (what, t1["length"], t0["r0"]))
    why = "; ".join(problems) or None
    report = {"why": why, "scored_on": os.path.basename(tp),
              "d_residues": sorted(d), "stereocentres": stereo,
              "closing": closing}
    try:
        with open(topo.path("checks", "json"), "w") as fh:
            json.dump(report, fh, indent=1)
    except OSError:
        pass
    for c in closing:
        log(f"    {topo.prefix}: closing {c['bond']} bond {c['before']:.2f} A "
            f"({c['energy_before']:.1f} kcal/mol) before the minimisation, "
            f"{c['after']:.2f} A ({c['energy_after']:.2f} kcal/mol) after")
    if d and not why:
        log(f"    {topo.prefix}: every D centre is still D after the "
            f"minimisation ({len(d)} residue(s))")
    return why, report


class _Topo:
    """The tleap built topologies one pose's energy is computed from.

    Built once and handed to whichever engine computes the energy. That is
    the point of it: the engines then provably read the same atoms, the same
    force field and the same protonation, so when their numbers differ the
    topology is not a candidate explanation.
    """

    def __init__(self, workdir, prefix, com_pdb, rec_mask, lig_mask,
                 n_rec_res, n_lig_res, asl, home, env, radii=RADII,
                 lig_chain=LIGAND_CHAIN):
        self.workdir, self.prefix, self.com_pdb = workdir, prefix, com_pdb
        self.rec_mask, self.lig_mask, self.asl = rec_mask, lig_mask, asl
        self.n_rec_res, self.n_lig_res = n_rec_res, n_lig_res
        self.home, self.env, self.radii = home, env, radii
        self.spec = {"ring": None, "d": []}
        # The chain the peptide is written under (peptide_chain).
        self.lig_chain = lig_chain

    def path(self, kind, ext):
        return os.path.join(self.workdir, "%s_%s.%s" % (self.prefix, kind, ext))


def _topologies(rec, pose_mol, workdir, prefix, log=print,
                prep_ph: float | None = 7.4, radii: str = RADII):
    """Build the complex, receptor and ligand topologies. None if it failed.

    Every check here guards one thing: that the molecule the energy is
    computed on is the molecule that was docked. They run once per pose
    rather than once per engine, so two engines cannot disagree about it.
    """
    # tleap and sander run with cwd=workdir, so every path given to them
    # must be absolute or a bare basename: a relative path resolves
    # inside the working directory.
    workdir = os.path.abspath(workdir)
    os.makedirs(workdir, exist_ok=True)
    home = amber_home()
    if home is None:
        return None
    env = dict(os.environ, AMBERHOME=home)
    spec = peptide_spec(pose_mol)
    # The peptide's own chain, one the receptor does not use, and recorded
    # beside the pose for whatever reads the complex later (the panel).
    used = receptor_chains(rec)
    lig_chain = peptide_chain(used)
    if lig_chain != LIGAND_CHAIN:
        log(f"    {prefix}: the receptor uses chain {LIGAND_CHAIN}, so the "
            f"peptide is chain {lig_chain} in the complex")
    rec_pdb, lig_pdb, com_pdb = _write_pdbs(rec, pose_mol, workdir, prefix,
                                            lig_chain=lig_chain)
    import json as _json
    with open(os.path.join(workdir, f"{prefix}_chains.json"), "w") as fh:
        _json.dump({"peptide": lig_chain, "receptor": sorted(used)}, fh)
    if prep_ph is not None:
        com_pdb = prep_complex(com_pdb, workdir, prefix, prep_ph, log,
                               protonation=getattr(rec, "protonation", None),
                               lig_chain=lig_chain)
    leap = os.path.join(workdir, f"{prefix}.leap")
    if spec["ring"]:
        # The closing bond, and for a head to tail ring no termini.
        try:
            text = _ring_leap(com_pdb, spec["ring"], workdir, prefix, radii,
                              lig_chain=lig_chain)
        except (ValueError, KeyError) as exc:
            log(f"    {prefix}: the ring cannot be built ({exc}), skipped")
            return None
    else:
        text = LEAP.format(complex_pdb=os.path.basename(com_pdb),
                           receptor_pdb=os.path.basename(rec_pdb),
                           ligand_pdb=os.path.basename(lig_pdb),
                           prefix=prefix, radii=radii)
    with open(leap, "w") as fh:
        fh.write(text)
    # A topology left by an earlier attempt must not stand in for this one.
    for ext in ("parm7", "rst7"):
        stale = os.path.join(workdir, f"{prefix}_com.{ext}")
        if os.path.isfile(stale):
            os.remove(stale)
    for name in (f"{prefix}_com_tleap.rst7", f"{prefix}_checks.json"):
        stale = os.path.join(workdir, name)
        if os.path.isfile(stale):
            os.remove(stale)
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
    rec_mask, lig_mask, n_rec_res, n_lig_res, asl = masks_by_chain(
        com_pdb, lig_chain)
    if not n_rec_res or not n_lig_res:
        log(f"    {prefix}: chain {lig_chain} holds {n_lig_res} residue(s) "
            f"and the receptor {n_rec_res} in {os.path.basename(com_pdb)} - "
            f"the complex cannot be split, skipped")
        return None
    log(f"    {prefix}: {asl}")
    log(f"    {prefix}: receptor_mask {rec_mask[:60]}  ligand_mask {lig_mask}")
    if spec["d"]:
        placed = mirror_hydrogens(workdir, prefix, com_pdb, spec, n_rec_res,
                                  home, env, radii, log, lig_chain=lig_chain)
        if placed is None:
            return None
        log(f"    {prefix}: D residue(s) {', '.join(map(str, spec['d']))} "
            f"built with the L template; {placed['atoms_placed']} "
            f"added atom(s) placed as their mirror image"
            + (f", {placed['hydrogens_turned']} more turned to the open side "
               f"of their centre" if placed["hydrogens_turned"] else ""))
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

    topo = _Topo(workdir, prefix, com_pdb, rec_mask, lig_mask,
                 n_rec_res, n_lig_res, asl, home, env, radii,
                 lig_chain=lig_chain)
    topo.spec = spec
    return topo


def _minimise(topo, log=print):
    """Relax the pose with sander. Returns the trajectory to score.

    The restraint holds the receptor so the added hydrogens and the peptide
    relax while the prepared geometry stays put, and the generalised Born
    model is the one the radii select. On failure the raw pose is returned,
    so a caller always has something to score.

    Extracted so every engine relaxes the pose the SAME way, which turned out
    to matter more than the energy model. On sander's own minimum OpenMM
    lands within 1.3 kcal/mol of MMPBSA.py; OpenMM minimising the same pose
    itself lands about 30 out. sander's 1000 cycles is a shallow relaxation
    and OpenMM's L-BFGS converges well inside 1000 iterations, so matching
    the step count does not match the amount of relaxation, and the extra
    relaxation loosens the interface. The geometry is produced once, here.
    """
    workdir, prefix = topo.workdir, topo.prefix
    radii = getattr(topo, "radii", RADII)
    mininp = os.path.join(workdir, f"{prefix}_min.in")
    with open(mininp, "w") as fh:
        fh.write(MIN_IN.format(rec_mask=topo.rec_mask,
                               igb=model_of(radii)[1]))
    log(f"    {prefix}: minimising the complex")
    minrst = os.path.join(workdir, f"{prefix}_min.rst7")
    sander = os.path.join(topo.home, "bin", "sander")
    subprocess.run([sander, "-O", "-i", mininp,
                    "-p", f"{prefix}_com.parm7", "-c", f"{prefix}_com.rst7",
                    "-ref", f"{prefix}_com.rst7",
                    "-r", os.path.basename(minrst),
                    "-o", f"{prefix}_min.out"],
                   cwd=workdir, capture_output=True, text=True,
                   env=topo.env)
    return (os.path.basename(minrst) if os.path.isfile(minrst)
            else f"{prefix}_com.rst7")


def run_amber(rec, pose_mol, workdir: str, prefix: str, log=print,
              prep_ph: float | None = 7.4,
              intdiel: float = INTDIEL,
              radii: str = RADII) -> float | None:
    """The reference engine: sander minimises, MMPBSA.py decomposes.

    ff14SB with whichever generalised Born model the radii select. It is
    the default, and the engine every run used before others could be
    chosen.
    """
    topo = _topologies(rec, pose_mol, workdir, prefix, log, prep_ph, radii)
    if topo is None:
        return None
    workdir, home, env = topo.workdir, topo.home, topo.env
    rec_mask, lig_mask = topo.rec_mask, topo.lig_mask
    minrst = os.path.join(workdir, f"{prefix}_min.rst7")
    traj = _minimise(topo, log)
    # Which coordinates MMPBSA.py is about to be given, written down beside the
    # pose. The metric panel measures the geometry the energy belongs to, and
    # it cannot work that out afterwards from the files alone.
    with open(os.path.join(workdir, f"{prefix}_scored_on.txt"), "w") as fh:
        fh.write(traj + "\n")
    # The peptide that is scored must be the one that was docked: its
    # stereocentres as docked, and a ring's closing bond a bond.
    why, _rep = check_peptide(topo, traj, log)
    if why:
        log(f"    {prefix}: not scored: {why}")
        return None
    if traj != os.path.basename(minrst):
        # sander wrote no restart file at all, so the energy below is computed
        # on the complex as tleap built it: hydrogens unrelaxed, and the pose
        # where the docking left it. That is a different quantity from every
        # other pose's energy. It used to be logged as "did not converge",
        # which reads as a detail about a minimisation that never produced
        # anything.
        log(f"    {prefix}: sander wrote no minimised coordinates, so this "
            f"energy is for the UNMINIMISED complex and is not comparable "
            f"with the other poses")
    else:
        # Whether the minimisation reached a gradient or ran out of cycles.
        # The panel reports this, and the energy means less when the
        # structure was still moving.
        _m = _min_state(os.path.join(workdir, f"{prefix}_min.out"))
        if _m:
            log(f"    {prefix}: minimised to RMS gradient {_m['rms']:.3f} "
                f"kcal/mol/A in {_m['cycles']} cycle(s)"
                + (", stopped at maxcyc" if _m["maxed"] else ""))
        # What the minimisation started from, which indicates whether
        # the pose was usable at all. A pose that begins at a large
        # positive energy is one the minimiser rescued rather than
        # one that was plausible.
        start = _first_energy(os.path.join(workdir, f"{prefix}_min.out"))
        junction = _at_ring_junction(topo)
        if start is not None and start > CLASHED_START and junction:
            # A ring's two ends as the docking leaves them: its restraint
            # pulls the end residues' alpha carbons to about 2.2 A, where a
            # peptide bond holds them 3.8 A apart, so the hydrogens tleap adds
            # there overlap. That is inside the peptide, which is free to
            # move, and the closing bond and the stereocentres are checked
            # after the minimisation (check_peptide).
            log(f"    {prefix}: the pose started at {start:.2e} kcal/mol, "
                f"from {junction} where the ring closes, crowded as the "
                f"docking leaves it; the peptide is free to move there")
        elif start is not None and start > CLASHED_START:
            log(f"    {prefix}: the pose started at {start:.2e} kcal/mol, so it "
                f"is inside the receptor, and the receptor is held fixed, so "
                f"the clash cannot relax out. Any energy from it describes "
                f"the overlap and not the binding.")

    inp = os.path.join(workdir, f"{prefix}_mmpbsa.in")
    with open(inp, "w") as fh:
        fh.write(MMPBSA_IN.format(rec_mask=rec_mask, lig_mask=lig_mask,
                                  igb=model_of(radii)[1]))
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


# The restraint weight sander applies through restraint_wt, kept here only so
# the two protocols can be compared in one place. The minimisation is
# sander's.
OMM_RESTRAINT_KCAL_A2 = 5.0


def _omm_model(radii):
    """The OpenMM implicit solvent matching the model these radii select.

    The model follows the radii exactly as it does for AmberTools, so a run
    that chose GBn2 is scored under GBn2 by either engine. Hardcoding OBC(II)
    would score a GBn2 run under a different model and still return a number.
    """
    import openmm.app as app
    name, _igb = model_of(radii)
    return (app.GBn2 if name == "gbn2" else app.OBC2), name


def run_openmm(rec, pose_mol, workdir: str, prefix: str, log=print,
               prep_ph: float | None = 7.4,
               intdiel: float = INTDIEL,
               radii: str = RADII) -> float | None:
    """The same topology and the same geometry as AmberTools, scored again.

    What this replaces is MMPBSA.py, not sander. The pose is relaxed by
    sander, exactly as run_amber relaxes it, and OpenMM then computes that
    geometry's energy with its own generalised Born code.

    WHY IT DOES NOT MINIMISE FOR ITSELF. It was written that way first and
    came out 24 to 32 kcal/mol from AmberTools on a 5 residue peptide and 52
    to 65 on a 19 residue one, with the pose order reversed on the longer one.
    Measured back to its cause: on sander's own minimum OpenMM agrees with
    MMPBSA.py to 1.3 kcal/mol across five poses with the order unchanged, so
    the energy models agree and the minimisation was the whole difference.
    Minimising further made the answer worse, so it is not offered.
    """
    try:
        import openmm                                          # noqa: F401
    except Exception as exc:                                   # noqa: BLE001
        log(f"    {prefix}: openmm is not installed ({exc}), skipped")
        return None
    if amber_home() is None:
        log(f"    {prefix}: the OpenMM engine scores the geometry sander "
            f"minimises, and AmberTools was not found, so there is no "
            f"calibrated geometry to score. Skipped rather than minimised "
            f"here, which would put the energy about 30 kcal/mol out.")
        return None
    topo = _topologies(rec, pose_mol, workdir, prefix, log, prep_ph, radii)
    if topo is None:
        return None
    traj = _minimise(topo, log)
    if traj != f"{prefix}_min.rst7":
        log(f"    {prefix}: minimisation did not converge, using the raw pose")
    why, _rep = check_peptide(topo, traj, log)
    if why:
        log(f"    {prefix}: not scored: {why}")
        return None
    return openmm_dg(topo, log, intdiel, traj)


def openmm_dg(topo, log=print, intdiel: float = INTDIEL,
              rst: str | None = None) -> float | None:
    """OpenMM's MM-GBSA energy for one already built, already relaxed pose.

    `rst` is the coordinate file to score; without one the minimised
    coordinates are used if present and the raw pose otherwise. Nothing is
    minimised here.

    Single trajectory, like MMPBSA.py: the receptor and ligand energies come
    from the SAME frame as the complex by slicing its coordinates, so the
    difference is a binding energy and not a folding one.
    """
    import openmm as mm
    import openmm.app as app
    import openmm.unit as u

    prefix = topo.prefix
    kcal = u.kilocalorie_per_mole
    gb, model = _omm_model(getattr(topo, "radii", RADII))

    def energy(parm7, coords):
        prm = app.AmberPrmtopFile(parm7)
        system = prm.createSystem(
            implicitSolvent=gb, soluteDielectric=float(intdiel),
            solventDielectric=78.5,
            implicitSolventSaltConc=0.15 * u.molar,
            nonbondedMethod=app.NoCutoff, constraints=None)
        integ = mm.VerletIntegrator(1.0 * u.femtosecond)
        try:
            sim = app.Simulation(prm.topology, system, integ,
                                 mm.Platform.getPlatformByName("CPU"))
        except Exception:                                      # noqa: BLE001
            sim = app.Simulation(prm.topology, system, integ)
        sim.context.setPositions(coords)
        return sim.context.getState(getEnergy=True
                                    ).getPotentialEnergy().value_in_unit(kcal)

    com7 = topo.path("com", "parm7")
    rec7, lig7 = topo.path("rec", "parm7"), topo.path("lig", "parm7")
    if rst:
        start = rst if os.path.isabs(rst) else os.path.join(topo.workdir, rst)
    else:
        start = topo.path("min", "rst7")
        if not os.path.isfile(start):
            start = topo.path("com", "rst7")
    if not all(os.path.isfile(f) for f in (com7, rec7, lig7, start)):
        log(f"    {prefix}: the topologies are incomplete, skipped")
        return None

    pos = app.AmberInpcrdFile(start).positions
    n_rec = app.AmberPrmtopFile(rec7).topology.getNumAtoms()
    n_lig = app.AmberPrmtopFile(lig7).topology.getNumAtoms()
    if n_rec + n_lig != len(pos):
        # The split topologies must account for every atom of the complex, or
        # the slice takes the wrong atoms and still returns a number.
        log(f"    {prefix}: the split topologies hold {n_rec} + {n_lig} atoms "
            f"against {len(pos)} in the complex, so the frame cannot be "
            f"sliced, skipped")
        return None

    dg = (energy(com7, pos) - energy(rec7, pos[:n_rec])
          - energy(lig7, pos[n_rec:]))
    log(f"    {prefix}: dG {dg:.2f} kcal/mol "
        f"({os.path.basename(start)}, OpenMM {model})")
    if abs(dg) > 500:
        log(f"    {prefix}: dG {dg:.0f} kcal/mol is not physical "
            "(unrelieved clashes), discarded")
        return None
    return dg


# The binding energy Prime reports that means what the other two engines'
# numbers mean: dG_Bind(NS), "no strain", taken with ligand and receptor in
# their complex geometry, which is how single-trajectory MM-GBSA is defined.
# Prime's headline dG_Bind adds a ligand strain term from minimising the free
# peptide separately, and for a flexible peptide the free peptide collapses in
# implicit solvent into a shape it does not hold: the strain was +8.9 kcal/mol
# on a 5 residue peptide and +59.3 on a 15 residue one, turning -62.2 into
# -2.9 and ranking longer peptides down for reasons unrelated to binding.
_PRIME_DG = "r_psp_MMGBSA_dG_Bind(NS)"


def prime_forcefield(log=print, timeout=30):
    """Which OPLS version a Prime job started now would actually get.

    Returns (flags, name). The force field licence is a separate pool from
    the rescoring one, often with a single seat, so a job started while no
    seat is free falls back to OPLS_2005 on its own. Asking first means
    the result can say which force field produced it, instead of the table
    implying every value came from the same one.
    """
    from . import gbengines
    eng = gbengines.SCHRODINGER
    free, total = eng.forcefield_licences(timeout)
    if free is None:
        # Unreadable is not the same as absent. Let Prime pick its own
        # default and record that the choice was not verified.
        return [], eng.FF_MODERN + " (assumed; the licence was not readable)"
    if free > 0:
        return [], eng.FF_MODERN
    log(f"    the {eng.FF_FEATURE} licence is fully in use "
        f"({total} seat(s)), so this runs on {eng.FF_LEGACY}")
    return ["-prime_opt", "OPLS_VERSION=" + eng.FF_LEGACY], eng.FF_LEGACY


def run_prime(rec, pose_mol, workdir: str, prefix: str, log=print,
              prep_ph: float | None = 7.4,
              intdiel: float = INTDIEL,
              radii: str = RADII) -> float | None:
    """Schrodinger Prime MM-GBSA: OPLS with VSGB 2.0.

    A different force field and a different solvent model from the other two
    engines, so its numbers belong in their own column and are never averaged
    or compared by value with them. It is here because groups holding the
    licence often need it for comparability with their existing work.

    The topologies are still built by tleap, for one reason: the protonation
    and the atoms are then identical to what the other engines read, so a
    difference in the answer belongs to the force field and not to the
    preparation. Prime is given the prepared complex, not the raw pose.
    """
    from . import gbengines
    spec = peptide_spec(pose_mol)
    if spec["ring"] or spec["d"]:
        # The caller refuses these before they get here; this is the guard
        # that keeps a ring from ever being scored as a linear chain.
        log(f"    {prefix}: Prime does not score a peptide with D residues or "
            f"a ring here, skipped")
        return None
    schro = gbengines.find_schrodinger()
    if not schro:
        log(f"    {prefix}: no Schrodinger install with prime_mmgbsa was "
            f"found. Looked at $SCHRODINGER and "
            f"{', '.join(gbengines.SCHRODINGER_ROOTS)}. Skipped.")
        return None
    # The suite's own tools read $SCHRODINGER, so it is exported to them
    # whether or not the module was loaded in this shell. Inside a batch job
    # nothing has loaded it, and the engine used to report itself missing on
    # a machine that has the whole suite.
    env = dict(os.environ, SCHRODINGER=schro)
    plan = gbengines.SCHRODINGER.plan()
    if plan.get("readable") and not plan.get("jobs"):
        log(f"    {prefix}: every {gbengines.SCHRODINGER.FEATURE} licence is "
            f"in use, and one job needs "
            f"{gbengines.SCHRODINGER.TOKENS_PER_JOB}, skipped")
        return None

    topo = _topologies(rec, pose_mol, workdir, prefix, log, prep_ph, radii)
    if topo is None:
        return None
    return prime_dg(topo, log, schro, env)


def prime_dg(topo, log=print, schro=None, env=None) -> float | None:
    """The Prime half of run_prime, on a complex that is already prepared.

    Separate for the same reason openmm_dg is: it can then run against a
    complex an earlier AmberTools run left on disk, which is how this engine
    was exercised without repeating the preparation. On one test pose it
    returned -81.79 kcal/mol (no strain) on S-OPLS where AmberTools gave
    -48.07, which is the scale difference the separate column exists for.
    """
    from . import gbengines
    prefix = topo.prefix
    schro = schro or gbengines.find_schrodinger()
    if not schro:
        log(f"    {prefix}: no Schrodinger install found, skipped")
        return None
    env = env or dict(os.environ, SCHRODINGER=schro)

    # Prime reads Maestro formats, so the prepared complex is converted
    # rather than rebuilt. It is given the PDB2PQR output, which has standard
    # residue names and the run's protonation: the file tleap reads carries
    # Amber's own names (CYX, HID, HIE, LYN), which are not Schrodinger's,
    # and the one Prime run that has been checked end to end used this file.
    src = os.path.join(topo.workdir, prefix + "_com_prepped.pdb")
    if not os.path.isfile(src):
        src = topo.com_pdb
    mae = os.path.join(topo.workdir, prefix + "_com.mae")
    conv = os.path.join(schro, "utilities", "structconvert")
    if not os.path.isfile(conv):
        log(f"    {prefix}: structconvert not found under $SCHRODINGER, "
            f"skipped")
        return None
    c = subprocess.run([conv, src, mae], cwd=topo.workdir,
                       capture_output=True, text=True, env=env)
    if not os.path.isfile(mae):
        err = (c.stderr or c.stdout or "").strip().splitlines()
        log(f"    {prefix}: the complex could not be converted for Prime"
            f"{(' (' + err[-1][:80] + ')') if err else ''}, skipped")
        return None

    flags, ff = prime_forcefield(log)
    log(f"    {prefix}: Prime MM-GBSA on {ff} with VSGB 2.0")
    cmd = [os.path.join(schro, "prime_mmgbsa"), os.path.basename(mae),
           "-ligand", "chain.name " + getattr(topo, "lig_chain",
                                               LIGAND_CHAIN),
           "-jobname", prefix, "-WAIT", "-NOJOBID", "-OVERWRITE"] + flags
    r = subprocess.run(cmd, cwd=topo.workdir, capture_output=True, text=True,
                       env=env)

    csv = os.path.join(topo.workdir, prefix + "-out.csv")
    if not os.path.isfile(csv):
        err = (r.stderr or r.stdout or "").strip().splitlines()
        log(f"    {prefix}: Prime produced no output table"
            f"{(' (' + err[-1][:110] + ')') if err else ''}, skipped")
        return None

    dg = _read_prime_csv(csv)
    if dg is None:
        # Not the headline dG_Bind in its place: that is a different
        # quantity, and a column mixing the two would be read as one.
        log(f"    {prefix}: {_PRIME_DG} was not in Prime's output, skipped")
        return None
    if abs(dg) > 500:
        log(f"    {prefix}: dG {dg:.0f} kcal/mol is not physical "
            "(unrelieved clashes), discarded")
        return None
    log(f"    {prefix}: dG {dg:.2f} kcal/mol ({ff})")
    # Which force field this number came from, beside the number itself. The
    # column holds two different quantities otherwise, with nothing saying so.
    try:
        with open(os.path.join(topo.workdir,
                               prefix + ".forcefield"), "w") as fh:
            fh.write(ff + "\n")
    except OSError:
        pass
    return dg


def _read_prime_csv(path):
    """The binding energy from Prime's output table, or None."""
    import csv as _csv
    try:
        with open(path, newline="") as fh:
            for row in _csv.DictReader(fh):
                v = row.get(_PRIME_DG)
                if v not in (None, ""):
                    try:
                        return float(v)
                    except ValueError:
                        return None
    except OSError:
        return None
    return None


def _one_pose(rec, mol, workdir, seq, rank, prep_ph, log,
              intdiel=INTDIEL, radii=RADII, engine="amber"):
    """One (construct, pose) unit, in a directory of its own.

    Its own directory because tleap always writes leap.log into the working
    directory under that name, and MMPBSA.py writes _MMPBSA_* temporaries
    under theirs: two units sharing a directory overwrite each other's
    working files, and the failure looks like a chemistry problem rather
    than a collision.
    """
    wd = os.path.join(workdir, seq, "pose%d" % rank)
    return engine_runner(engine)(rec, mol, wd, seq, log, prep_ph, intdiel,
                                 radii)


def _read_cache(path, intdiel=None, radii=None):
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
    if radii is not None:
        # A cache written before the radii were recorded holds energies
        # computed with tleap's default, LEGACY_RADII.
        was = LEGACY_RADII
        for line in text.splitlines():
            bits = line.split()
            if line.startswith("#") and "radii" in bits[:-1]:
                was = bits[bits.index("radii") + 1]
        if was != radii:
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


def _write_cache(path, values, intdiel=None, radii=None):
    try:
        with open(path, "w") as fh:
            if intdiel is not None:
                fh.write("# igb %d saltcon 0.150 intdiel %g%s\n"
                         % (model_of(radii)[1], float(intdiel),
                            (" radii %s" % radii) if radii else ""))
            for rank in sorted(values):
                fh.write("%d %.6f\n" % (rank, values[rank]))
    except OSError:
        pass


def engines_on_disk(workdir, seq, intdiel=None, radii=None):
    """{engine: [dG, best first]} for every engine that has scored this
    construct in this folder.

    Read from the caches rather than from a result file, so a rescoring done
    later with another engine shows up without anything having to remember
    to record it. That is what makes a second engine visible in the table
    after the fact: the numbers are already on disk under their own name.
    """
    from . import gbengines
    out = {}
    for name in gbengines.ORDER:
        path = os.path.join(workdir, cache_name(seq, name))
        if not os.path.isfile(path):
            continue
        # Settings are not checked here. A cache written under different
        # radii is still that engine's answer for this construct, and the
        # per engine column is a comparison, not the ranking key.
        got = _read_cache(path)
        vals = [got[r] for r in sorted(got)]
        if vals:
            out[name] = vals
    return out


def radii_for(workdir, default=RADII):
    """The radii a run's MM-GBSA uses: the ones it started with.

    They are recorded in a RADII file in the working directory. Without one,
    energies already in the folder, or in another MM-GBSA folder of the same
    run beside it, decide: the radii their cache names, or tleap's default for
    energies computed before caches named them. Every other run gets RADII.
    """
    marker = os.path.join(workdir, "RADII")
    try:
        got = open(marker).read().strip()
    except OSError:
        got = ""
    if got:
        return got
    radii = default
    here = os.path.abspath(workdir)
    parent = os.path.dirname(here)
    for d in [here] + [x for x in sorted(glob.glob(os.path.join(parent, "mmgbsa*")))
                       if os.path.abspath(x) != here]:
        if not os.path.isdir(d):
            continue
        try:
            other = open(os.path.join(d, "RADII")).read().strip()
        except OSError:
            other = ""
        if other:
            radii = other
            break
        dgs = sorted(glob.glob(os.path.join(d, "*.dg")))
        if dgs:
            # Energies there already: the radii their cache names, or, from
            # before caches named them, tleap's default.
            radii = LEGACY_RADII
            for f in dgs:
                try:
                    bits = open(f).readline().split()
                except OSError:
                    continue
                if bits[:1] == ["#"] and "radii" in bits[:-1]:
                    radii = bits[bits.index("radii") + 1]
                    break
            break
    try:
        os.makedirs(workdir, exist_ok=True)
        with open(marker, "w") as fh:
            fh.write(radii + "\n")
    except OSError:
        pass
    return radii


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


def _at_ring_junction(topo):
    """The atom with the largest force at the first step of the
    minimisation, as "atom HA of residue 1", when the peptide is a head to
    tail ring and that atom is in its first or last residue; else ""."""
    spec = getattr(topo, "spec", None) or {}
    ring = spec.get("ring") or {}
    if ring.get("kind") != "backbone":
        return ""
    try:
        txt = open(topo.path("min", "out"), errors="replace").read()
        i = txt.index("NSTEP       ENERGY")
        bits = txt[i:].splitlines()[1].split()
        atom = int(bits[-1]) - 1
        from . import amberfiles as af
        res = af.residues(af.read_parm7(topo.path("com", "parm7")))
    except (OSError, ValueError, IndexError):
        return ""
    first = topo.n_rec_res + 1
    for r, _lab, atoms in res:
        for nm, k in atoms.items():
            if k == atom:
                pos = r - first + 1
                if pos in (1, topo.n_lig_res):
                    return "atom %s of residue %d" % (nm, pos)
                return ""
    return ""


def _min_state(min_out):
    """Where the minimisation stopped, or None.

    sander stops at maxcyc rather than at a gradient for almost every pose,
    which is expected of a restrained rescoring minimisation and is not the
    same thing as having converged. Nothing read this file before, and the
    energy of a pose that was still moving means less than one that settled.
    """
    last, maxed = None, False
    try:
        with open(min_out, errors="replace") as fh:
            for line in fh:
                bits = line.split()
                if (len(bits) >= 4 and bits[0].isdigit()
                        and all("E" in b for b in bits[1:4])):
                    try:
                        last = (int(bits[0]), float(bits[2]))
                    except ValueError:
                        pass
                elif "Maximum number of minimization cycles reached" in line:
                    maxed = True
    except OSError:
        return None
    if last is None:
        return None
    return {"cycles": last[0], "rms": last[1], "maxed": maxed}


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
        intdiel: float = INTDIEL, radii: str | None = None,
        model: str | None = None, engine: str = "") -> dict:
    """Rescore each construct. Returns sequence -> dG, missing where it failed.

    `model` ("obc2" or "gbn2") picks the generalised Born model for a run whose
    MM-GBSA has not started. A run that has started keeps the model it began
    with, whatever is asked for, because its model is its radii and a run
    keeps those; see MODELS.

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
        log("  SKIPPING MM-GBSA: no MM-GBSA engine is installed. AmberTools "
            "needs tleap and MMPBSA.py on PATH; Prime needs $SCHRODINGER set, "
            "which on a cluster usually means module load Schrodinger.")
        return {}
    # An engine named by the caller is honoured; otherwise the default above.
    from . import gbengines
    engine = (engine or "").strip().lower() or BACKEND_ENGINE.get(b, b)
    eng = gbengines.get(engine)
    engine = eng.name
    if not eng.available():
        usable = gbengines.installed()
        log(f"  SKIPPING MM-GBSA: the {eng.label} engine was asked for and is "
            f"not installed here. Available: "
            f"{', '.join(usable) if usable else 'none'}.")
        return {}
    if radii:
        # Asked for by name: the folder's energies are recomputed under it,
        # and the record says so.
        os.makedirs(workdir, exist_ok=True)
        with open(os.path.join(workdir, "RADII"), "w") as fh:
            fh.write(radii + "\n")
    else:
        want = MODELS.get(model or DEFAULT_MODEL, MODELS[DEFAULT_MODEL])
        radii = radii_for(workdir, default=want["radii"])
        if model and model_of(radii)[0] != model:
            log(f"  this run's MM-GBSA started with "
                f"{MODELS[model_of(radii)[0]]['label']} and keeps it, so its "
                f"energies stay comparable with each other; start a new run "
                f"for {MODELS[model]['label']}")
    mname, igb = model_of(radii)
    if engine == "schrodinger":
        log(f"  MM-GBSA engine: {eng.label}, {eng.solvation}. Its numbers are "
            f"on their own scale and go in the {eng.column} column.")
    elif engine == "openmm":
        log(f"  MM-GBSA engine: {eng.label}, {eng.forcefield} with "
            f"{MODELS[mname]['label']}, {radii} radii and an interior "
            f"dielectric of {float(intdiel):g}, scoring the geometry sander "
            f"minimised")
    else:
        log(f"  MM-GBSA engine: {eng.label}, {MODELS[mname]['label']} "
            f"(igb={igb}) with {radii} radii and an interior dielectric of "
            f"{float(intdiel):g}")
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
        have = _read_cache(os.path.join(workdir, cache_name(seq, engine)),
                           intdiel,
                           radii)
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
                           intdiel, radii, engine)
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
                _write_cache(os.path.join(workdir, cache_name(seq, engine)),
                             results[seq],
                             intdiel, radii)

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
                                n=n, spread=spread, aggregate=aggregate,
                                radii=radii, model=model_of(radii)[0])
        if n > 1:
            log(f"    {seq}: dG = {dg:.2f} kcal/mol "
                f"({aggregate} of {n} pose(s), spread {spread:.1f})")
        else:
            log(f"    {seq}: dG = {dg:.2f} kcal/mol")
    if not out:
        log("  MM-GBSA produced no results; the ADCP ranking is unaffected")
    return out
