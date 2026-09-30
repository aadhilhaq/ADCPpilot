"""Drawings for the Modifications dialog: one residue as it sits in a
peptide, as an SVG, with wavy bonds where the chain continues.

By Aadhil Haq

The twenty standard residues are written here; the non-standard ones come
from static/modifications.json, whose structures were taken from ADCP 1.1's
own templates and rotamers (scripts/modifications_table.py). Charged side
chains are drawn charged, as at pH 7 and as ADCP models them, so that
variants that differ only in charge (SEP and SEP2) look different. A D
residue is drawn as the mirror image of the L one, which is how ADCP builds
it: every stereocentre is inverted, double bonds are not.
"""
import functools
import json
import os

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog("rdApp.*")

# L residues, * where the chain continues. Histidine is the N-epsilon-H
# tautomer; Lys, Arg, Asp and Glu are charged.
STANDARD = {
    "A": "*N[C@@H](C)C(*)=O",
    "R": "*N[C@@H](CCCNC(N)=[NH2+])C(*)=O",
    "N": "*N[C@@H](CC(N)=O)C(*)=O",
    "D": "*N[C@@H](CC(=O)[O-])C(*)=O",
    "C": "*N[C@@H](CS)C(*)=O",
    "Q": "*N[C@@H](CCC(N)=O)C(*)=O",
    "E": "*N[C@@H](CCC(=O)[O-])C(*)=O",
    "G": "*NCC(*)=O",
    "H": "*N[C@@H](Cc1c[nH]cn1)C(*)=O",
    "I": "*N[C@H](C(*)=O)[C@@H](C)CC",
    "L": "*N[C@@H](CC(C)C)C(*)=O",
    "K": "*N[C@@H](CCCC[NH3+])C(*)=O",
    "M": "*N[C@@H](CCSC)C(*)=O",
    "F": "*N[C@@H](Cc1ccccc1)C(*)=O",
    "P": "*C(=O)[C@@H]1CCCN1*",
    "S": "*N[C@@H](CO)C(*)=O",
    "T": "*N[C@H](C(*)=O)[C@@H](C)O",
    "W": "*N[C@@H](Cc1c[nH]c2ccccc12)C(*)=O",
    "Y": "*N[C@@H](Cc1ccc(O)cc1)C(*)=O",
    "V": "*N[C@H](C(*)=O)C(C)C",
}
WIDTH, HEIGHT = 320, 240

# Every residue is laid out with its backbone the same way round: N on the
# left, the carbonyl on the right, the side chain up from the alpha carbon.
_TEMPLATE = Chem.MolFromSmiles("*NC(C)C(*)=O")
_XY = [(-2.598, 0.0), (-1.299, -0.75), (0.0, 0.0), (0.0, 1.5), (1.299, -0.75),
       (2.598, 0.0), (1.299, -2.25)]
_conf = Chem.Conformer(_TEMPLATE.GetNumAtoms())
for _i, (_x, _y) in enumerate(_XY):
    _conf.SetAtomPosition(_i, (_x, _y, 0.0))
_TEMPLATE.AddConformer(_conf)
_BACKBONE = Chem.MolFromSmarts("[#0][#7][CX4][CX3](=O)[#0]")


@functools.lru_cache(maxsize=1)
def library():
    """{code: entry} from static/modifications.json."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "static",
                        "modifications.json")
    try:
        rows = json.load(open(path))["modifications"]
    except (OSError, ValueError, KeyError):
        rows = []
    return {r["code"]: r for r in rows if r.get("smiles")}


def mirror(mol):
    """The mirror image: every tetrahedral centre inverted."""
    m = Chem.Mol(mol)
    for a in m.GetAtoms():
        t = a.GetChiralTag()
        if t == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    return m


def residue(letter="", code="", d=False):
    """The residue as an RDKit molecule. KeyError for a residue that is not
    known here; ValueError for a D form that does not exist (glycine)."""
    if code:
        smiles = library()[code]["smiles"]
    else:
        smiles = STANDARD[letter]
        if d and letter == "G":
            raise ValueError("glycine has no D form")
    mol = Chem.MolFromSmiles(smiles)
    return mirror(mol) if d else mol


def _crowded(mol):
    """True when two atoms that are not bonded sit almost on each other, as a
    large side chain can when the backbone is held in place."""
    P = mol.GetConformer().GetPositions()[:, :2]
    d = np.linalg.norm(P[:, None, :] - P[None, :, :], axis=2)
    for b in mol.GetBonds():
        d[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] = d[b.GetEndAtomIdx(), b.GetBeginAtomIdx()] = 9
    np.fill_diagonal(d, 9)
    return bool(d.min() < 0.75)


def _upright(mol):
    """Turn a free layout so the backbone reads left to right, N first, with
    the side chain above it."""
    hit = mol.GetSubstructMatch(_BACKBONE)
    if not hit:
        return
    n, ca, c = hit[1], hit[2], hit[3]
    conf = mol.GetConformer()
    P = conf.GetPositions()
    v = P[c] - P[n]
    t = -np.arctan2(v[1], v[0])
    R = np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])
    P = (P - P[ca]) @ R.T
    side = [x.GetIdx() for x in mol.GetAtomWithIdx(ca).GetNeighbors()
            if x.GetIdx() not in (n, c)]
    if side and P[side[0]][1] < 0:
        P[:, 1] *= -1
    for i, p in enumerate(P):
        conf.SetAtomPosition(i, (float(p[0]), float(p[1]), 0.0))


def _layout(mol):
    try:
        rdDepictor.GenerateDepictionMatching2DStructure(mol, _TEMPLATE,
                                                        acceptFailure=True)
        if not _crowded(mol):
            return
    except (ValueError, RuntimeError):
        pass
    rdDepictor.Compute2DCoords(mol)
    _upright(mol)


@functools.lru_cache(maxsize=512)
def svg(letter="", code="", d=False, dark=False):
    """An SVG drawing of the residue, stereocentres wedged and labelled R or
    S, on a transparent background."""
    mol = residue(letter, code, d)
    _layout(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(WIDTH, HEIGHT)
    opts = drawer.drawOptions()
    if dark:
        try:
            rdMolDraw2D.SetDarkMode(opts)
        except (AttributeError, TypeError):
            pass
    opts.clearBackground = False
    opts.dummiesAreAttachments = True
    opts.addStereoAnnotation = True
    opts.fixedBondLength = 40        # one scale for every residue, unless it must shrink to fit
    opts.padding = 0.05
    opts.bondLineWidth = 2
    opts.minFontSize = 14
    opts.annotationFontScale = 0.8   # the R and S labels, readable at card size
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()
