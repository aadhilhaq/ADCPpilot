"""Build adcp_pilot/static/modifications.json, the list the Modifications
dialog offers, from an ADCP 1.1 installation.

By Aadhil Haq

    OPSIN_JAR=/path/to/opsin.jar python scripts/modifications_table.py /path/to/adcpsuite

For each non-standard side chain in ADCP's Swiss library it records the code,
a name, the parent residue, ADCP's coarse residue, the code of the D form,
the net charge ADCP models it with, and its structure as a residue SMILES
with the two chain attachments written as "*".

The structure comes from ADCP itself: atoms, elements and bonds from its
OpenMM templates (bond orders follow from the explicit hydrogens and the
template's net charge), and side-chain stereochemistry from the first
rotamer in its rotamer library. The alpha carbon is L. The few side chains
ADCP has no template for are taken from the RCSB Chemical Component
Dictionary, after checking their atoms against ADCP's rotamer library.

Names come from the Chemical Component Dictionary where the code is one, and
otherwise from PubChem, looked up by structure. A name is kept only when the
stereochemistry it states agrees with ADCP's residue and, when OPSIN (a
name-to-structure parser; Java) is given, when it describes the same
molecule. Structures that are also in the dictionary are compared with it
and the differences printed. Needs RDKit and network access.
"""
import concurrent.futures as cf
import collections
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdCIPLabeler, rdDetermineBonds
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Geometry import Point3D

RDLogger.DisableLog("rdApp.*")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
         "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
         "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V"}
# AutoDock atom types in the rotamer library, as elements
AD_ELEMENT = {"A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O",
              "OS": "O", "O": "O", "S": "S", "SA": "S", "F": "F", "Cl": "Cl",
              "CL": "Cl", "Br": "Br", "BR": "Br", "I": "I", "P": "P", "Se": "Se",
              "H": "H", "HD": "H", "HS": "H"}
# typical bond lengths to CB, for placing CB in the rotamer frame
TO_CB = {"C": 1.53, "N": 1.47, "O": 1.43, "S": 1.82, "F": 1.35, "Cl": 1.77,
         "Br": 1.94, "Se": 1.95}
UNCHARGER = rdMolStandardize.Uncharger()
ALPHA = Chem.MolFromSmarts("[#7][CX4:1]([!#1])[CX3](=O)")


# --- ADCP's files -------------------------------------------------------------

def load_templates(paths):
    """{residue: (atoms [(name, element)], bonds, external atom names,
    net charge)} from OpenMM force field files, the first file winning."""
    out = {}
    for path in paths:
        root = ET.parse(path).getroot()
        elem = {t.get("name"): t.get("element") for t in root.find("AtomTypes")}
        charge = {a.get("type"): float(a.get("charge"))
                  for a in root.find("NonbondedForce").findall("Atom")
                  if a.get("charge") is not None}
        for r in root.find("Residues"):
            name = r.get("name")
            if name in out:
                continue
            atoms = [(a.get("name"), elem.get(a.get("type"))) for a in r.findall("Atom")]
            q = [charge.get(a.get("type")) for a in r.findall("Atom")]
            if None in [e for _, e in atoms] or None in q:
                continue
            bonds = [(b.get("atomName1"), b.get("atomName2")) for b in r.findall("Bond")]
            ext = [atoms[int(e.get("from"))][0] for e in r.findall("ExternalBond")]
            out[name] = (atoms, bonds, ext, int(round(sum(q))))
    return out


def load_rotamers(path):
    """{residue: (atom names, elements, first rotamer {name: xyz}, every
    rotamer)}. The frame has the alpha carbon at the origin."""
    out = {}
    lines = open(path).read().split("\n")
    for i, line in enumerate(lines):
        if not line.startswith("rotamer "):
            continue
        p = line.split()
        nrot, nat = int(p[3]), int(p[4])
        types = lines[i + 2].split()
        names = lines[i + 4].split()
        xyz = [list(map(float, lines[i + 5 + k].split())) for k in range(nrot * nat)]
        rots = [dict(zip(names, xyz[k * nat:(k + 1) * nat])) for k in range(nrot)]
        out[p[1]] = (names, [AD_ELEMENT.get(t, t) for t in types], rots[0], rots)
    return out


def fit_cb(templates, rotamers):
    """CB is not in the rotamer library; it is one fixed point in its frame.
    Least squares on the distance from CB to every atom bonded to it, in
    every rotamer, and 1.53 A from the alpha carbon."""
    pts, dist = [], []
    for code, (atoms, bonds, _ext, _q) in templates.items():
        if code not in rotamers:
            continue
        el = dict(atoms)
        for a1, a2 in bonds:
            if "CB" not in (a1, a2):
                continue
            x = a2 if a1 == "CB" else a1
            if x == "CA" or el.get(x) in (None, "H"):
                continue
            for rot in rotamers[code][3]:
                if x in rot:
                    pts.append(rot[x])
                    dist.append(TO_CB.get(el[x], 1.5))
    X, D = np.array(pts), np.array(dist)
    cb = np.array([0.5, 0.0, 1.4])
    for _ in range(50):
        v = cb - X
        r = np.linalg.norm(v, axis=1)
        J = np.vstack([v / r[:, None], cb / np.linalg.norm(cb)])
        res = np.concatenate([r - D, [np.linalg.norm(cb) - 1.53]])
        cb = cb - np.linalg.lstsq(J, res, rcond=None)[0]
    rms = float(np.sqrt(((np.linalg.norm(cb - X, axis=1) - D) ** 2).mean()))
    return cb, rms, len(X)


def _rotation(axis, t):
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(t) * K + (1 - np.cos(t)) * K @ K


def backbone(cb):
    """N, C and O of an L residue in the rotamer frame: a real L-alanine,
    turned so that its CB lies on the fitted CB and its N in the frame's
    mirror plane (y = 0), through which ADCP reflects its D rotamers."""
    ala = Chem.AddHs(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))
    AllChem.EmbedMolecule(ala, randomSeed=7)
    P = ala.GetConformer().GetPositions()
    pts = {"N": P[0] - P[1], "CB": P[2] - P[1], "C": P[3] - P[1], "O": P[4] - P[1]}
    a = pts["CB"] / np.linalg.norm(pts["CB"])
    b = cb / np.linalg.norm(cb)
    ax = np.cross(a, b)
    if np.linalg.norm(ax) > 1e-9:
        R = _rotation(ax / np.linalg.norm(ax), np.arccos(np.clip(a @ b, -1, 1)))
        pts = {k: R @ v for k, v in pts.items()}
    best = None
    for deg in range(360):
        n = _rotation(b, np.radians(deg)) @ pts["N"]
        score = abs(n[1]) + (0 if n[0] < 0 else 9)
        if best is None or score < best[0]:
            best = (score, deg)
    R = _rotation(b, np.radians(best[1]))
    out = {k: R @ v for k, v in pts.items()}
    out["CB"], out["CA"] = cb, np.zeros(3)
    return out


# --- structures ------------------------------------------------------------

def alpha_is_l(mol):
    """True for an L alpha carbon, False for D, None when there is no
    alpha hydrogen to judge by. Anchored at the chain ends of a residue (the
    two *), or at the carboxylic acid of a free amino acid, so that a
    second amino acid in a side chain (aminomalonate, diaminopimelate) is
    not read as the backbone."""
    ends = "[#0]" if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()) else ""
    tail = "[#0]" if ends else "[OX2H1]"
    for chir, verdict in (("@@", True), ("@", False)):
        patt = Chem.MolFromSmarts("%s[#7][C%sH]([!#1;!#0])[CX3](=O)%s" % (ends, chir, tail))
        if mol.HasSubstructMatch(patt, useChirality=True):
            return verdict
    return None


def sulfoxides_as_double_bonds(mol):
    """[S+]-[O-] as S=O, the way the structure is usually drawn."""
    rw = Chem.RWMol(mol)
    for b in rw.GetBonds():
        for s, o in ((b.GetBeginAtom(), b.GetEndAtom()), (b.GetEndAtom(), b.GetBeginAtom())):
            if (s.GetSymbol() in ("S", "Se") and s.GetFormalCharge() == 1
                    and o.GetSymbol() == "O" and o.GetFormalCharge() == -1
                    and b.GetBondType() == Chem.BondType.SINGLE):
                b.SetBondType(Chem.BondType.DOUBLE)
                s.SetFormalCharge(0)
                o.SetFormalCharge(0)
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    return m


def from_template(code, templates, rotamers, bb):
    """The residue from ADCP's template, with * where the chain continues."""
    atoms, bonds, ext, q = templates[code]
    rw = Chem.RWMol()
    idx = {}
    for name, el in atoms:
        a = Chem.Atom(el)
        a.SetNoImplicit(True)
        a.SetProp("_n", name)
        idx[name] = rw.AddAtom(a)
    for a1, a2 in bonds:
        rw.AddBond(idx[a1], idx[a2], Chem.BondType.SINGLE)
    ends = {rw.GetAtomWithIdx(idx[e]).GetSymbol(): e for e in ext}
    if set(ends) != {"N", "C"}:
        raise ValueError("external bonds are not on N and C: %s" % ext)
    # methyl caps stand in for the rest of the chain while bond orders are
    # worked out, and become * afterwards
    for end in ("N", "C"):
        c = Chem.Atom(6)
        c.SetNoImplicit(True)
        c.SetProp("_n", "CAP" + end)
        ci = rw.AddAtom(c)
        rw.AddBond(idx[ends[end]], ci, Chem.BondType.SINGLE)
        for _ in range(3):
            h = Chem.Atom(1)
            h.SetNoImplicit(True)
            rw.AddBond(ci, rw.AddAtom(h), Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()), assignId=True)
    rdDetermineBonds.DetermineBondOrders(mol, charge=q, embedChiral=False)
    Chem.SanitizeMol(mol)
    heavy = Chem.RemoveHs(mol)
    # positions: the side chain from the first rotamer, the backbone placed
    pos = {k: np.asarray(v) for k, v in bb.items()}
    pos.update({k: np.asarray(v) for k, v in rotamers[code][2].items()})
    pos["CAPN"] = pos["N"] + np.array([-1.0, 0.2, 0.3])
    pos["CAPC"] = pos["C"] + np.array([1.0, -0.2, 0.4])
    conf = Chem.Conformer(heavy.GetNumAtoms())
    missing = []
    for a in heavy.GetAtoms():
        n = a.GetProp("_n")
        if n in pos:
            conf.SetAtomPosition(a.GetIdx(), Point3D(*map(float, pos[n])))
        else:
            missing.append(n)
    if missing:
        raise ValueError("no coordinates for %s" % missing)
    conf.Set3D(True)
    heavy.RemoveAllConformers()
    heavy.AddConformer(conf)
    Chem.AssignStereochemistryFrom3D(heavy)
    # A phosphate or sulfonate center looks chiral once its oxygens are
    # written as =O, O- and OH, but they are equivalent.
    for a in heavy.GetAtoms():
        if sum(1 for n in a.GetNeighbors()
               if n.GetSymbol() == "O" and n.GetDegree() == 1) >= 2:
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    rw = Chem.RWMol(heavy)
    for a in rw.GetAtoms():
        if a.GetProp("_n").startswith("CAP"):
            a.SetAtomicNum(0)
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(0)
    mol = rw.GetMol()
    mol.RemoveAllConformers()
    Chem.SanitizeMol(mol)
    mol = sulfoxides_as_double_bonds(mol)
    mol = without_cage_tags(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if alpha_is_l(mol) is False:
        raise ValueError("the alpha carbon came out D")
    return mol, q


def without_cage_tags(mol):
    """A bridgehead of a cage (adamantane) has the only configuration the
    cage allows. RDKit still reads one from 3D, and CIP even labels it, but
    drawn it only adds wedges that say nothing. Tags in a bridged ring
    system (two rings sharing three or more atoms) are cleared; ADCP's
    library has no bridged side chain whose stereochemistry is free."""
    m = Chem.Mol(mol)
    rings = [set(r) for r in m.GetRingInfo().AtomRings()]
    for i, r1 in enumerate(rings):
        for r2 in rings[i + 1:]:
            if len(r1 & r2) >= 3:
                for k in r1 | r2:
                    m.GetAtomWithIdx(k).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    return m


def from_dictionary(code, smiles, rotamers):
    """The residue from the dictionary's SMILES, for a side chain ADCP has
    no template for: checked atom for atom against the rotamer library,
    and with the alpha carbon made L."""
    free = UNCHARGER.uncharge(Chem.MolFromSmiles(smiles))
    rxn = AllChem.ReactionFromSmarts(
        "[N;!H0:1]-[C:2]-[C:3](=[O:4])-[OX2H1]>>*-[N:1]-[C:2]-[C:3](=[O:4])-*")
    prods = rxn.RunReactants((free,))
    if not prods:
        raise ValueError("no amino acid backbone in %s" % smiles)
    mol = prods[0][0]
    Chem.SanitizeMol(mol)
    if alpha_is_l(mol) is False:
        mol.GetAtomWithIdx(mol.GetSubstructMatch(ALPHA)[1]).InvertChirality()
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    # the atoms: backbone N, CA, C, O, then CB and the rotamer atoms
    elements = [e for e in rotamers[code][1] if e != "H"]
    want = collections.Counter(["N", "C", "C", "O", "C"] + elements)
    have = collections.Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
    if want != have:
        raise ValueError("atoms differ from ADCP's: %s vs %s" % (dict(have), dict(want)))
    return mol, int(sum(a.GetFormalCharge() for a in mol.GetAtoms()))


def free_acid(res):
    """The residue as a free, neutral amino acid: * on N becomes H, * on
    the carbonyl becomes OH."""
    rw = Chem.RWMol(res)
    for a in rw.GetAtoms():
        if a.GetAtomicNum() == 0:
            if a.GetNeighbors()[0].GetSymbol() == "N":
                a.SetAtomicNum(1)
            else:
                a.SetAtomicNum(8)
                a.SetNoImplicit(False)
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    return UNCHARGER.uncharge(Chem.RemoveHs(m))


# --- names -------------------------------------------------------------------

def tidy(name):
    """The dictionary writes some names in capitals; restore the letters
    that carry meaning once they are lowercased (D-, L-, N-, O-, S-, N5-,
    4H-), and write stereo descriptors in capitals whatever the source."""
    if name.isupper():
        name = name.lower()
    name = re.sub(r"(^|[\s,(\-])([dl])-(?=[a-z0-9(\[])",
                  lambda m: m.group(1) + m.group(2).upper() + "-", name)
    name = re.sub(r"(^|[\s,(\-\d'])([nos])(\d*)-(?=[a-z(\[])",
                  lambda m: m.group(1) + m.group(2).upper() + m.group(3) + "-", name)
    name = re.sub(r"(\d)h-", lambda m: m.group(1) + "H-", name)   # 4H-thieno
    name = re.sub(r"^((?:\d+[rs],)*\d+[rs])-", lambda m: m.group(1).upper() + "-", name)
    return re.sub(r"\(((?:\d+[rsez]|[rsez])(?:,\s*(?:\d+[rsez]|[rsez]))*)\)",
                  lambda m: "(" + m.group(1).upper() + ")", name, flags=re.I)


def soften(title):
    """PubChem capitalizes words in a title ("(S)-2-Amino-3-..."); write them
    the way the dictionary's names read. Element symbols, stereo prefixes and
    locants are left alone."""
    return re.sub(r"(^|[\s\-(\[])([A-Z])(?=[a-z]{2})",
                  lambda m: m.group(1) + m.group(2).lower(), title)


def usable(title):
    """A PubChem title that names this one compound: not a racemate, an
    optical rotation, an index entry ("L-valine, 3-fluoro-"), a catalog
    number (SCHEMBL123456, NPC168081), or written in capitals."""
    return not (re.search(r"\bDL-|\brac\b|\((?:\+|-|\u00b1|\+/-|RS)\)|azaniumyl|,\s*\S+-$",
                          title)
                or re.search(r"[A-Z]{2,}", title))


def plain(name):
    """Remove annotations a source appends to a name."""
    name = re.sub(r"\s*\((?:non-preferred name|preferred name)\)\s*$", "", name)
    return re.sub(r"^M-(?=[a-z])", "m-", name)          # meta, not a locant


def opsin(names, jar):
    """{name: SMILES} for the names OPSIN can read; empty without Java or
    the jar (then names are checked by their stereochemistry only)."""
    names = sorted(set(n for n in names if n))
    if not (jar and os.path.isfile(jar) and shutil.which("java") and names):
        return {}
    p = subprocess.run(["java", "-jar", jar, "-osmi"], input="\n".join(names) + "\n",
                       capture_output=True, text=True, timeout=900)
    return {n: s.strip() for n, s in zip(names, p.stdout.split("\n")) if s.strip()}


def describes(name, free, parsed):
    """False when OPSIN reads the name as a different molecule: another
    skeleton, or every stereocenter given and one of them different."""
    smiles = parsed.get(name)
    if not smiles:
        return True
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return True
    m = UNCHARGER.uncharge(m)
    ours, theirs = Chem.MolToInchiKey(free), Chem.MolToInchiKey(m)
    if ours.split("-")[0] != theirs.split("-")[0]:
        return False
    centres = Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                        useLegacyImplementation=False)
    if centres and all(c != "?" for _, c in centres):
        return ours.split("-")[1] == theirs.split("-")[1]
    return True


def cip(free):
    """{atom index: 'R' or 'S'}, and the alpha carbon's index or None."""
    m = Chem.Mol(free)
    rdCIPLabeler.AssignCIPLabels(m)
    labels = {a.GetIdx(): a.GetProp("_CIPCode") for a in m.GetAtoms()
              if a.HasProp("_CIPCode") and a.GetProp("_CIPCode") in ("R", "S")}
    hit = m.GetSubstructMatch(Chem.MolFromSmarts("[#7][CX4:1][CX3](=O)[OX2H1]"))
    return labels, (hit[1] if hit else None)


def stereo_agrees(name, free):
    """False when the name states a configuration the structure does not
    have: D or L for the alpha carbon, or R or S at locant 2."""
    labels, alpha = cip(free)
    is_l = alpha_is_l(free)
    for m in re.finditer(r"(?:^|[\s,(\-])([DL])-(?=[A-Za-z])", name):
        if is_l is not None and (m.group(1) == "L") != is_l:
            return False
    for m in re.finditer(r"\(([0-9A-Za-z,\s']+)\)", name):
        for item in (x.strip().upper() for x in m.group(1).split(",")):
            want = None
            if re.fullmatch(r"2[RS]", item):
                want = item[1]
            elif item in ("R", "S") and len(labels) == 1:
                want = item
            if want and alpha is not None and labels.get(alpha) not in (None, want):
                return False
    return True


def swap_alpha(name):
    """A dictionary name for the other configuration at the alpha carbon:
    D and L swapped, and R and S at locant 2 or in a bare (R) or (S)."""
    name = re.sub(r"(^|[\s,(\-])([DL])-(?=[A-Za-z])",
                  lambda m: m.group(1) + ("L" if m.group(2) == "D" else "D") + "-", name)
    return re.sub(r"\((2?)([RS])\)|\(2([RS]),",
                  lambda m: ("(%s%s)" % (m.group(1), "S" if m.group(2) == "R" else "R"))
                  if m.group(2) else "(2%s," % ("S" if m.group(3) == "R" else "R"), name)


def mirror_name(name):
    """The name of the mirror image, when the name states configurations:
    every R and S, and D and L, swapped (reflection inverts every center;
    E, Z, cis, trans, allo and threo are unchanged). None when the name
    states none, so the page says "D form of" instead."""
    flip = {"R": "S", "S": "R"}
    swapped = [False]

    def item(m):
        swapped[0] = True
        return m.group(1) + flip[m.group(2)]

    def paren(m):
        return "(" + re.sub(r"(\d*|alpha|beta|gamma|delta)([RS])\b", item, m.group(1)) + ")"

    out = re.sub(r"\(((?:(?:\d*|alpha|beta|gamma|delta)[RSEZ])(?:,\s*(?:\d*|alpha|beta|gamma|delta)[RSEZ])*)\)",
                 paren, name)
    out = re.sub(r"^((?:\d+[RS],)*\d+[RS])-",
                 lambda m: re.sub(r"(\d+)([RS])", item, m.group(1)) + "-", out)

    def dl(m):
        swapped[0] = True
        return m.group(1) + ("L" if m.group(2) == "D" else "D") + "-"
    out = re.sub(r"(^|[\s,(\-])([DL])-(?=[A-Za-z0-9(\[])", dl, out)
    return out if swapped[0] else None


def rcsb(code):
    try:
        with urllib.request.urlopen(
                "https://data.rcsb.org/rest/v1/core/chemcomp/" + code, timeout=30) as r:
            d = json.load(r)
    except Exception:                                           # noqa: BLE001
        return {}
    cc = d.get("chem_comp") or {}
    return {"name": cc.get("name"), "parents": cc.get("mon_nstd_parent_comp_id") or [],
            "smiles": (d.get("rcsb_chem_comp_descriptor") or {}).get("SMILES_stereo")}


def pubchem(smiles):
    """(title, IUPAC name) of the compound with exactly this structure.
    Retried on a busy server; a structure PubChem does not have is (None,
    None) at once."""
    data = urllib.parse.urlencode({"smiles": smiles}).encode()
    for attempt in range(4):
        try:
            with urllib.request.urlopen(urllib.request.Request(
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
                    "property/Title,IUPACName/JSON", data=data), timeout=30) as r:
                p = json.load(r)["PropertyTable"]["Properties"][0]
            return p.get("Title"), p.get("IUPACName")
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None, None
        except Exception:                                       # noqa: BLE001
            pass
        finally:
            time.sleep(0.25 + attempt)   # PubChem asks for no more than 5 a second
    return None, None


# --- the table -------------------------------------------------------------

def main(adcpsuite, opsin_jar=None):
    data = glob.glob(os.path.join(adcpsuite, "lib", "python3*", "site-packages",
                                  "ADCP", "data"))[0]
    rot_dir, ff_dir = os.path.join(data, "rotamers"), os.path.join(data, "openMMff")
    templates = load_templates([os.path.join(ff_dir, "swiss_ff.xml"),
                                os.path.join(ff_dir, "sannerlab_ff.xml")])
    rotamers = load_rotamers(os.path.join(rot_dir, "swiss.lib"))
    cb, rms, n = fit_cb(templates, rotamers)
    print("CB in the rotamer frame %s, rms %.3f A over %d atoms" % (cb.round(3), rms, n))
    bb = backbone(cb)

    coarse = {}
    for line in open(os.path.join(rot_dir, "swiss.lib")):
        if line.startswith("rotamer "):
            p = line.split()
            coarse[p[1]] = p[2].upper()
    rows = []
    for line in open(os.path.join(rot_dir, "SWISS_ROTAMER_LIST.dat")):
        if line.startswith("#") or not line.split():
            continue
        p = line.split()
        # A side chain ADCP gives no coarse residue cannot be written in its
        # sequence syntax, so it is not offered.
        if not coarse.get(p[0]) or coarse[p[0]] == "O":
            continue
        rows.append({"code": p[0], "d": p[3] if p[2] == "Y" else None,
                     "coarse": coarse[p[0]]})
    with cf.ThreadPoolExecutor(8) as ex:
        ccd = dict(zip([r["code"] for r in rows], ex.map(rcsb, [r["code"] for r in rows])))

    report = collections.defaultdict(list)
    flat = lambda m: Chem.MolToSmiles(m, isomericSmiles=False)       # noqa: E731
    for r in rows:
        code, info = r["code"], ccd.get(r["code"]) or {}
        mol, q, source = None, 0, None
        try:
            if code in templates and code in rotamers:
                mol, q = from_template(code, templates, rotamers, bb)
                source = "template"
            elif info.get("smiles"):
                mol, q = from_dictionary(code, info["smiles"], rotamers)
                source = "dictionary"
                report["taken from the dictionary"].append(code)
            else:
                raise ValueError("no template and no dictionary entry")
        except Exception as e:                                  # noqa: BLE001
            report["no structure"].append("%s: %s" % (code, e))
        r["smiles"] = Chem.MolToSmiles(mol) if mol is not None else None
        r["charge"] = q
        r["_free"] = free = free_acid(mol) if mol is not None else None
        ref = Chem.MolFromSmiles(info["smiles"]) if info.get("smiles") else None
        if free is not None and ref is not None and source == "template":
            ref = UNCHARGER.uncharge(ref)
            if flat(ref) != flat(free):
                report["differs from the dictionary in atoms or bonds"].append(
                    "%s: %s vs %s" % (code, flat(ref), flat(free)))
            elif Chem.MolToSmiles(ref) != Chem.MolToSmiles(free):
                report["differs from the dictionary in stereochemistry"].append(
                    "%s: %s vs %s" % (code, Chem.MolToSmiles(ref), Chem.MolToSmiles(free)))
            else:
                report["same as the dictionary"].append(code)
        parent = next((THREE[x] for x in info.get("parents", []) if x in THREE), None)
        r["parent"] = parent or r["coarse"]
        name = soften(tidy(info["name"])) if info.get("name") else None
        if name and free is not None and not stereo_agrees(name, free):
            fixed = swap_alpha(name)
            report["dictionary name for the other alpha configuration"].append(
                "%s: %s, now %s" % (code, name, fixed))
            name = fixed
        r["name"] = name
    # Names, in order: the dictionary's, PubChem's title for exactly this
    # structure, PubChem's IUPAC name. Each must agree with the residue's
    # stereochemistry and, where OPSIN can read it, describe this molecule:
    # the dictionary's names are not always right (BIF is biphenylalanine,
    # named as a cyclohexyl), and PubChem's titles include catalog numbers.
    cache = {}
    for r in rows:
        if r["_free"] is not None:
            key = Chem.MolToSmiles(r["_free"])
            if key not in cache:
                cache[key] = pubchem(key)
            title, iupac = cache[key]
            r["_cands"] = [(r["name"] and plain(r["name"]), "dictionary"),
                           (title and usable(title) and plain(tidy(soften(title))), "PubChem"),
                           (iupac and usable(iupac) and plain(iupac), "PubChem IUPAC")]
            for t, how in ((title, "PubChem"), (iupac, "PubChem IUPAC")):
                if t and not usable(t):
                    report["%s name refused as written" % how].append("%s: %s" % (r["code"], t))
    parsed = opsin([c for r in rows for c, _ in r.get("_cands", [])], opsin_jar)
    if not parsed:
        print("OPSIN not available: names are checked by stereochemistry only")
    for r in rows:
        dictionary_name, r["name"] = r["name"], None
        for cand, how in r.get("_cands", []):
            if not cand:
                continue
            if not stereo_agrees(cand, r["_free"]):
                report["name refused, stereochemistry disagrees"].append(
                    "%s: %s (%s)" % (r["code"], cand, how))
                continue
            if not describes(cand, r["_free"], parsed):
                report["name refused, OPSIN reads another molecule"].append(
                    "%s: %s (%s)" % (r["code"], cand, how))
                continue
            r["name"], r["_how"] = cand, how
            if how != "dictionary":
                report["named by %s" % how + (" instead of the dictionary"
                                              if dictionary_name else "")].append(
                    "%s: %s" % (r["code"], cand))
            break
        if not r["name"]:
            report["no name found"].append(r["code"])
    # Variants of one molecule (SEP and SEP2 differ only in charge, 2HF1
    # and 2HF2 in tautomer) take one name, the dictionary's where there is
    # one; the option label then tells them apart by charge.
    first = {}
    for r in sorted(rows, key=lambda x: x.get("_how") != "dictionary"):
        if r["name"] and r["_free"] is not None:
            first.setdefault(Chem.MolToSmiles(r["_free"], isomericSmiles=True), r["name"])
    for r in rows:
        if r["_free"] is not None:
            same = first.get(Chem.MolToSmiles(r["_free"], isomericSmiles=True))
            if same and same != r["name"]:
                report["named as the same molecule under another code"].append(
                    "%s: %s, was %s" % (r["code"], same, r["name"]))
                r["name"] = same
    # The D form's name: the L name mirrored, kept when OPSIN, where it can
    # read it, finds the mirror image of the residue.
    dnames = {r["code"]: mirror_name(r["name"]) for r in rows if r["name"] and r["d"]}
    parsed_d = opsin([n for n in dnames.values() if n], opsin_jar)
    for r in rows:
        r["dname"] = None
        n = dnames.get(r["code"])
        if n and r["_free"] is not None:
            mirrored = Chem.Mol(r["_free"])
            for a in mirrored.GetAtoms():
                a.InvertChirality()
            Chem.AssignStereochemistry(mirrored, cleanIt=True, force=True)
            if describes(n, mirrored, parsed_d):
                r["dname"] = n
            else:
                report["D name refused, OPSIN reads another molecule"].append(
                    "%s: %s" % (r["code"], n))
    out = [{"code": r["code"], "name": r["name"] or r["code"], "dname": r["dname"],
            "parent": r["parent"], "coarse": r["coarse"], "d": r["d"], "charge": r["charge"],
            "smiles": r["smiles"]} for r in rows]
    out.sort(key=lambda x: (x["parent"], x["name"].lower(), x["code"]))
    dest = os.path.join(ROOT, "adcp_pilot", "static", "modifications.json")
    json.dump({"source": "ADCP 1.1 Swiss side-chain library and OpenMM templates; "
                         "names from the RCSB Chemical Component Dictionary and PubChem",
               "modifications": out}, open(dest, "w"), indent=1)
    for k, v in report.items():
        print("%s: %d" % (k, len(v)))
        if k != "same as the dictionary":
            for x in v:
                print("   " + x)
    print("%d modifications written, %d with a structure" % (
        len(out), sum(1 for x in out if x["smiles"])))


if __name__ == "__main__":
    main(sys.argv[1], os.environ.get("OPSIN_JAR"))
