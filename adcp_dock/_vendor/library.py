"""Capped peptide library: ACE-X1..Xn-NME, built with RDKit.

By Aadhil Haq

Length is a parameter throughout, and every caller passes it from the config,
so nothing here is specific to one pipeline: 2 gives the 400 dipeptides, 3 the
8,000 tripeptides, 4 the 160,000 tetrapeptides."""
from __future__ import annotations
import itertools, os, random
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")

AA = "ACDEFGHIKLMNPQRSTVWY"

# Backbone unit for residue i, written from N to C. The side chain is
# substituted into {R}. Glycine has no CB; proline's nitrogen is part of its
# ring.
_GENERIC = "N[C@@H]({R})C(=O)"
_GLY     = "NCC(=O)"
_PRO     = "N1CCC[C@H]1C(=O)"

SIDECHAIN = {
    "A": "C",
    "C": "CS",
    "D": "CC(=O)O",
    "E": "CCC(=O)O",
    "F": "Cc1ccccc1",
    "H": "Cc1c[nH]cn1",
    "I": "[C@@H](C)CC",
    "K": "CCCCN",
    "L": "CC(C)C",
    "M": "CCSC",
    "N": "CC(N)=O",
    "Q": "CCC(N)=O",
    "R": "CCCNC(N)=N",
    "S": "CO",
    "T": "[C@@H](C)O",
    "V": "C(C)C",
    "W": "Cc1c[nH]c2ccccc12",
    "Y": "Cc1ccc(O)cc1",
}

def _unit(a: str) -> str:
    if a == "G": return _GLY
    if a == "P": return _PRO
    return _GENERIC.format(R=SIDECHAIN[a])

# Heavy-atom layout of each fragment, in SMILES writing order, which is also
# RDKit's atom order. Needed because Meeko emits generic 'UNL' atoms otherwise,
# which would make backbone/side-chain classification impossible downstream.
_LAYOUT = {
    "ACE": [("CH3", True), ("C", True), ("O", True)],          # CC(=O)
    "NME": [("N", True), ("CH3", True)],                        # NC
    "GLY": [("N", True), ("CA", True), ("C", True), ("O", True)],
    "PRO": [("N", True), ("CD", False), ("CG", False), ("CB", False),
            ("CA", True), ("C", True), ("O", True)],
}

# Standard PDB side-chain atom names, in the order the atoms appear in the
# side-chain SMILES above (RDKit's atom order). These must be real names: a
# placeholder such as "S1" makes viewers read the atom as sulfur and draw
# disulfides.
SIDECHAIN_NAMES = {
    "A": ["CB"],
    "C": ["CB", "SG"],
    "D": ["CB", "CG", "OD1", "OD2"],
    "E": ["CB", "CG", "CD", "OE1", "OE2"],
    "F": ["CB", "CG", "CD1", "CE1", "CZ", "CE2", "CD2"],
    "H": ["CB", "CG", "CD2", "NE2", "CE1", "ND1"],
    "I": ["CB", "CG2", "CG1", "CD1"],
    "K": ["CB", "CG", "CD", "CE", "NZ"],
    "L": ["CB", "CG", "CD1", "CD2"],
    "M": ["CB", "CG", "SD", "CE"],
    "N": ["CB", "CG", "ND2", "OD1"],
    "Q": ["CB", "CG", "CD", "NE2", "OE1"],
    "R": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "S": ["CB", "OG"],
    "T": ["CB", "CG2", "OG1"],
    "V": ["CB", "CG1", "CG2"],
    "W": ["CB", "CG", "CD1", "NE1", "CE2", "CZ2", "CH2", "CZ3", "CE3", "CD2"],
    "Y": ["CB", "CG", "CD1", "CE1", "CZ", "OH", "CE2", "CD2"],
}

def _sidechain_natoms(a: str) -> int:
    from rdkit import Chem as _C
    m = _C.MolFromSmiles(SIDECHAIN[a])
    return m.GetNumHeavyAtoms()

def atom_labels(seq: str, capped: bool = True) -> list[tuple[int, str, str, bool]]:
    """(resnum, resname, atomname, is_backbone) for every heavy atom, in RDKit order.

    Capped:   residue 1 is ACE, residues 2..N+1 the peptide, residue N+2 is NME.
    Uncapped: residues 1..N, with OXT on the last residue. ADCP writes uncapped
              peptides numbered from 1, so its poses need this form.
    """
    three = {"A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE","G":"GLY","H":"HIS",
             "I":"ILE","K":"LYS","L":"LEU","M":"MET","N":"ASN","P":"PRO","Q":"GLN",
             "R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR"}
    out = []
    first = 2 if capped else 1
    if capped:
        for nm, bb in _LAYOUT["ACE"]:
            out.append((1, "ACE", nm, True))
    for i, a in enumerate(seq, start=first):
        if a == "G":
            layout = _LAYOUT["GLY"]
        elif a == "P":
            layout = _LAYOUT["PRO"]
        else:
            # N, CA, <side chain>, C, O
            names = SIDECHAIN_NAMES[a]
            layout = ([("N", True), ("CA", True)] +
                      [(nm, False) for nm in names] +
                      [("C", True), ("O", True)])
        for nm, bb in layout:
            out.append((i, three[a], nm, bb))
    if capped:
        # NME follows the last residue; hardcoding a number only works for tetramers
        for nm, bb in _LAYOUT["NME"]:
            out.append((len(seq) + 2, "NME", nm, True))
    else:
        out.append((len(seq), three[seq[-1]], "OXT", True))   # free C-terminus
    return out


def free_smiles(seq: str) -> str:
    """Uncapped peptide, H2N-...-COOH. Used for validation."""
    return "".join(_unit(a) for a in seq) + "O"

def capped_smiles(seq: str) -> str:
    """ACE-seq-NME."""
    return "CC(=O)" + "".join(_unit(a) for a in seq) + "NC"

def validate_templates(verbose: bool = False) -> list[str]:
    """Every residue template must reproduce RDKit's own sequence parser.

    This catches inverted side-chain stereocentres, which are otherwise silent
    and corrupt every peptide containing that residue.
    """
    bad = []
    for a in AA:
        seq = a * 2
        mine = Chem.MolFromSmiles(free_smiles(seq))
        ref  = Chem.MolFromSequence(seq)
        if mine is None:
            bad.append(f"{a}: template does not parse"); continue
        if ref is None:
            continue                       # parser lacks this residue; skip
        cm = Chem.MolToSmiles(mine); cr = Chem.MolToSmiles(ref)
        if cm != cr:
            bad.append(f"{a}: {cm}  !=  {cr}")
        elif verbose:
            print(f"  {a} ok")
    return bad

# Only ever the signature default. pipeline.py passes
# length=cfg.peptide_length on the one call that matters, which is
# why this module is the same for every pipeline.
DEFAULT_LENGTH = 2

def sequences(n: int | None = None, seed: int = 1234,
              length: int = DEFAULT_LENGTH,
              max_aromatic: int = 0) -> list[str]:
    """The exhaustive 20^length library, or a reproducible random sample of n.

    Sampling must be random, not evenly spaced. The library is generated in
    lexicographic order, so taking every k-th sequence aliases against the 20^4
    structure: with n = 200 the step is exactly 800, a multiple of 400, and every
    sampled sequence ends in 'AA'. A seeded shuffle keeps runs reproducible while
    covering sequence space evenly.
    """
    allseq = ["".join(t) for t in itertools.product(AA, repeat=length)]
    # Vina's score largely tracks buried side-chain surface, so an
    # exhaustive scan returns a Trp/Tyr ladder whatever the site is.
    # Capping aromatics per peptide makes the search reward shape
    # complementarity rather than bulk. Off unless asked for.
    if max_aromatic:
        allseq = [q for q in allseq
                  if sum(ch in "FWY" for ch in q) <= max_aromatic]
    if n is None or n >= len(allseq):
        return allseq
    rng = random.Random(seed)
    return sorted(rng.sample(allseq, n))

# ------------------------------------------------------- pH protonation states

# Standard side-chain pKa values (EMBOSS set). Which atom carries the charge
# is fixed by SIDECHAIN_NAMES above.
IONISABLE = {
    #  resname:  (atom name, pKa, sign)   sign +1 = protonated below the pKa
    "ASP": ("OD2",  3.9, -1),
    "GLU": ("OE2",  4.1, -1),
    "CYS": ("SG",   8.5, -1),
    "TYR": ("OH",  10.1, -1),
    "LYS": ("NZ",  10.8, +1),
    "ARG": ("NH2", 12.5, +1),      # the imine nitrogen of the guanidine
    "HIS": ("ND1",  6.5, +1),      # the ring N that carries no hydrogen
}
NTERM_PKA, CTERM_PKA = 8.6, 3.6


def apply_ph_states(mol: "Chem.Mol", ph: float = 7.4, log=None) -> int:
    """Set side-chain protonation to match the pH the run is being done at.

    The SMILES in SIDECHAIN are written neutral: aspartate as COOH, lysine as
    NH2, arginine as a neutral guanidine. Docked as written, a probe carrying a
    formally neutral aspartate makes the wrong hydrogen bonds and is scored with
    the wrong desolvation, and at pH 5.5 a neutral histidine misses a cation the
    site may well be built around. Heavy atoms are untouched, so the atom
    labelling, the layout tables and the pose transplant all still line up; only
    formal charges change, and RDKit re-derives the hydrogens.

    Returns the number of residues changed.
    """
    n = 0
    for at in mol.GetAtoms():
        pi = at.GetMonomerInfo()
        if pi is None:
            continue
        rule = IONISABLE.get(pi.GetResidueName().strip())
        if rule is None:
            continue
        name, pka, sign = rule
        if pi.GetName().strip() != name:
            continue
        charged = (ph < pka) if sign > 0 else (ph > pka)
        if not charged:
            continue
        at.SetFormalCharge(sign)
        at.SetNoImplicit(False)
        if sign < 0:
            at.SetNumExplicitHs(0)
        n += 1
    if n:
        try:
            Chem.SanitizeMol(mol)
        except Exception:
            return 0
    return n


def build_3d(seq: str, seed: int = 0xF00D, capped: bool = True,
             ph: float | None = None) -> Chem.Mol | None:
    """Embed and MMFF-minimise one peptide of any length, capped or free-terminus.

    `ph` sets the side-chain protonation states before embedding. None keeps the
    neutral SMILES as written, which is what the geometry-only callers want.
    """
    mol = Chem.MolFromSmiles(capped_smiles(seq) if capped else free_smiles(seq))
    if mol is None:
        return None
    labels = atom_labels(seq, capped)
    if len(labels) != mol.GetNumHeavyAtoms():
        raise RuntimeError(
            f"{seq}: atom_labels produced {len(labels)} labels for "
            f"{mol.GetNumHeavyAtoms()} heavy atoms")
    for at, (rn, rname, aname, _bb) in zip(mol.GetAtoms(), labels):
        ri = Chem.AtomPDBResidueInfo()
        ri.SetResidueNumber(rn); ri.SetResidueName(rname)
        ri.SetName(aname[:4].ljust(4)); ri.SetIsHeteroAtom(False)
        at.SetMonomerInfo(ri)
    if ph is not None:
        apply_ph_states(mol, ph)
    mol = Chem.AddHs(mol, addCoords=False)
    hcount: dict[int, int] = {}                   # unique H names within each residue
    for at in mol.GetAtoms():                     # hydrogens inherit their parent residue
        if at.GetAtomicNum() != 1 or at.GetMonomerInfo() is not None:
            continue
        nb = at.GetNeighbors()
        if not nb:
            continue
        pi = nb[0].GetMonomerInfo()
        if pi is None:
            continue
        rn = pi.GetResidueNumber()
        hcount[rn] = hcount.get(rn, 0) + 1
        ri = Chem.AtomPDBResidueInfo()
        ri.SetResidueNumber(rn); ri.SetResidueName(pi.GetResidueName())
        ri.SetName(f"H{hcount[rn]}".ljust(4)); ri.SetIsHeteroAtom(False)
        at.SetMonomerInfo(ri)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = seed
    ps.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mol, ps) != 0:
        ps.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, ps) != 0:
            return None
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass
    mol.SetProp("_Name", seq)
    return mol

def write_pdbqt(mol: Chem.Mol, path: str) -> bool:
    """Ligand PDBQT via Meeko."""
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
    except ImportError as e:
        # meeko >= 0.7 imports gemmi at package load but does not declare it as a
        # dependency, so a bare "install meeko" message points at the wrong thing.
        raise RuntimeError(
            f"ligand preparation needs meeko and its dependencies; import failed with: {e}\n"
            "install with:  conda env create -f environment.yml   "
            "(or: pip install meeko gemmi)") from e
    prep = MoleculePreparation()
    setups = prep.prepare(mol)
    if not setups:
        return False
    txt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return False
    with open(path, "w") as fh:
        fh.write(txt)
    return True
