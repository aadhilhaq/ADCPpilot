"""What a pose is, beside what it scores.

By Aadhil Haq

A docking score orders poses. It does not say why one pose is better than
another, whether it sits on the site that was asked for, or whether the
number came from a large greasy contact or from two salt bridges. This
module computes that panel from what a run already has on disk:

  energy terms   the four MM-GBSA components, read out of the same MMPBSA.py
                 output the reported dG comes from, so they sum to it
  per residue    MMPBSA.py's own decomposition of the same frame, the same
                 topologies and the same GB model, rather than a second
                 estimate from a second force field. It agrees with the
                 reported total to within a few per cent, not exactly: see
                 decomposition() for why, and never present a per residue
                 energy as an exact share of the reported number
  interface      buried area, the share of the peptide that is buried, how
                 much of that area is apolar, typed contacts, clashes
  site           what the nominated residues contribute, and how much of the
                 nominated site the peptide actually touches

ONE GEOMETRY, for every number here: the sander minimised complex where
MM-GBSA ran, and the docked pose against the prepared receptor where it did
not. A metric measured on a different structure from the energy is not a
property of the same pose, and reporting the two side by side invites the
reader to relate them.

Nothing in here re-docks or re-minimises anything, and nothing here changes
the ranking. The energies are read, the geometry is measured.
"""
from __future__ import annotations

import math
import os
import re

import numpy as np

# Bondi van der Waals radii, and the water probe every SASA in this module is
# rolled with. Written down here because the numbers they produce are
# meaningless without them, and reported in the run's output for the same
# reason.
BONDI = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
         "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90,
         "ZN": 1.39, "MG": 1.73, "NA": 2.27, "K": 2.75, "CA": 2.31}
BONDI_DEFAULT = 1.70
PROBE = 1.4
# Points per sphere. Shrake and Rupley used 92.
#
# Measured, not assumed. Against an analytic two sphere area the error is 2.0%
# at 92 points, 0.56% at 256, 0.24% at 1024 and 0.10% at 4096; against
# freesasa on the same Bondi radii and probe, this converges to within 0.01%.
# The number that decides it is not accuracy but ROTATIONAL STABILITY, because
# the point lattice is fixed in space while the molecule is not: rotating a
# real complex moved its buried area by 15 A^2 at 256 points and by 3.6 A^2 at
# 1024. A buried area reported to a tenth of an A^2 from 256 points would be
# quoting four digits of a number whose third is noise, so this is 1024 and
# areas are reported to the nearest A^2.
POINTS = 1024

# The contact rules, identical to the ones the pose viewer paints with. See
# hope_monitor/interactions_kind.py in the pipelines, which holds the same
# constants and the symmetry tests for them: a rule written one way round
# reports half of the interactions it should.
HBOND_A = 3.5
IONIC_A = 4.0
STACK_A = 5.0
VDW_A = 4.0
# A pair closer than this fraction of the sum of its two radii. On a
# minimised complex the count should be zero; anything above it is a pose the
# minimiser could not pull apart, because the receptor is restrained.
CLASH_F = 0.80

POLAR = ("N", "O")
# The CHARGED forms only. Amber's neutral names are deliberately absent: HID
# and HIE are neutral histidines, LYN a neutral lysine, ASH and GLH neutral
# acids, and counting a pair of them as a salt bridge would report an ionic
# interaction between two uncharged side chains. HIS is kept because a
# structure that has not been through tleap names every histidine HIS, and
# then the protonation is unknown; CONSTANTS records that assumption.
BASIC = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",),
         "HIS": ("ND1", "NE2"), "HIP": ("ND1", "NE2")}
ACIDIC = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
AROM = {"PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
        "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
        "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
        "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
        "HID": ("CG", "ND1", "CD2", "CE1", "NE2"),
        "HIE": ("CG", "ND1", "CD2", "CE1", "NE2"),
        "HIP": ("CG", "ND1", "CD2", "CE1", "NE2")}

# Carbon and sulphur count as apolar, and a hydrogen counts as whatever it
# hangs off: see burial(), which folds each hydrogen's area onto its parent
# heavy atom rather than classing it on its own.
APOLAR_ELEM = ("C", "S")

#: The surface tension MMPBSA.py's GB nonpolar term uses, kcal/mol/A^2. It is
#: the conversion between the ESURF difference and a buried area, and it is
#: the default this package never changes.
SURFTEN = 0.0072

CONSTANTS = {
    "radii_set": "Bondi",
    "probe_radius_A": PROBE,
    "sphere_points": POINTS,
    "hbond_A": HBOND_A, "ionic_A": IONIC_A, "stack_A": STACK_A,
    "contact_A": VDW_A, "clash_fraction": CLASH_F,
    "apolar_elements": list(APOLAR_ELEM),
    "surften_kcal_per_A2": SURFTEN,
    "contact_rules": "this tool's own geometric definitions, not PLIP's: "
                     "proximity only, no hydrogen positions and no ring "
                     "orientation are tested",
    "histidine": "a residue named HIS is taken as protonated; HID and HIE "
                 "are read as the neutral forms they are",
}


# ----------------------------------------------------------------- structure

class Structure:
    """The atoms of one PDB, as parallel arrays.

    A class rather than a dict of arrays because every function here needs
    the same six fields and an index into them, and because a missing field
    should fail where it is read rather than three functions later.
    """

    def __init__(self, xyz, name, resname, chain, resnum, element,
                 icode=None):
        # reshape and the explicit dtypes are for the empty case: numpy infers
        # float64 from an empty list, and a float mask cannot be inverted or
        # and-ed, so every function here raised on a structure with no atoms.
        self.xyz = np.asarray(xyz, dtype=float).reshape(-1, 3)
        self.name = list(name)
        self.resname = list(resname)
        self.chain = list(chain)
        self.resnum = list(resnum)
        self.element = list(element)
        # The insertion code, which is part of a residue's identity: residues
        # 54 and 54A are two residues, and the rest of this package keys on
        # columns 23 to 27 for exactly that reason.
        self.icode = list(icode) if icode is not None else [""] * len(self.name)
        self.radii = np.array([BONDI.get(e, BONDI_DEFAULT)
                               for e in self.element], dtype=float)
        self.heavy = np.array([e != "H" for e in self.element], dtype=bool)

    def __len__(self):
        return len(self.name)

    def key(self, i):
        """(chain, residue number, residue name) for atom i.

        The number carries the insertion code where there is one, as "54A",
        so two residues sharing a number stay two residues.
        """
        return (self.chain[i], self.num(i), self.resname[i])

    def num(self, i):
        """The residue number of atom i, with its insertion code."""
        ic = (self.icode[i] or "").strip()
        return ("%d%s" % (self.resnum[i], ic)) if ic else self.resnum[i]

    def residues(self):
        """The residue keys in file order, without repeats."""
        out, seen = [], set()
        for i in range(len(self)):
            k = self.key(i)
            if k not in seen:
                seen.add(k)
                out.append(k)
        return out

    def resid(self, i):
        """(chain, numbered residue) for atom i, without the residue name."""
        return (self.chain[i], self.num(i))

    def select(self, mask):
        m = np.asarray(mask, dtype=bool)
        idx = np.nonzero(m)[0]
        return Structure(self.xyz[m], [self.name[i] for i in idx],
                         [self.resname[i] for i in idx],
                         [self.chain[i] for i in idx],
                         [self.resnum[i] for i in idx],
                         [self.element[i] for i in idx],
                         [self.icode[i] for i in idx])


# Every residue name that means "part of a protein chain", including the forms
# tleap and pdb2pqr rename to. Inside one of these, an atom name is a protein
# atom name and its element is its first letter: CA is the alpha carbon, not
# calcium, and CD, CE, CG, NZ and SD are likewise one letter elements.
AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX", "CYM", "TYM", "ARN",
    "ACE", "NME", "NHE", "NH2", "MSE", "SEC", "PYL", "HYP", "ASX", "GLX",
}
TWO_LETTER = ("CL", "BR", "ZN", "MG", "NA", "SE", "CA", "FE", "MN", "CU",
              "NI", "CO", "CD", "HG", "LI", "RB", "CS", "SR", "BA", "AL")


def _element(name, col, resname=""):
    """The element of an atom, from the PDB's own column where it has one.

    ambpdb leaves columns 77 to 78 empty, and Amber names a hydrogen 1HB as
    readily as HB1, so the name is read the way the PDB standard says: a
    leading digit is not part of the element. Inside an amino acid the first
    letter is the element outright, because guessing two letter elements there
    turns every CA into a calcium ion, with a 2.31 A radius, counted as polar,
    and inventing clashes against everything near it.
    """
    col = (col or "").strip().upper()
    if col:
        return col
    n = (name or "").strip().upper()
    n = n.lstrip("0123456789")
    if not n:
        return "H"
    if (resname or "").strip().upper() in AA3:
        return n[0]
    if n[:2] in TWO_LETTER:
        return n[:2]
    return n[0]


def read_pdb(path, chains=None):
    """Every ATOM and HETATM of a PDB, alternate locations other than the
    first dropped.

    The insertion code is part of a residue's identity, not decoration: two
    residues numbered 100 and 100A are different residues, and a key without
    it merges them and silently drops the second one's atoms.
    """
    xyz, name, resname, chain, resnum, element, icode = ([], [], [], [], [],
                                                         [], [])
    seen = set()
    with open(path) as fh:
        for l in fh:
            if not l.startswith(("ATOM", "HETATM")):
                continue
            nm = l[12:16].strip()
            ch = l[21]
            if chains and ch not in chains:
                continue
            try:
                rn = int(l[22:26])
            except ValueError:
                continue
            # The FIRST alternate location of each atom, whatever letter it
            # carries. Whitelisting " ", "" and "A" dropped an altLoc of "1",
            # which some depositions use, along with the ones it meant to.
            k = (ch, rn, l[26], nm)
            if k in seen:
                continue
            seen.add(k)
            rnm = l[17:20].strip().upper()
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
            name.append(nm)
            resname.append(rnm)
            chain.append(ch)
            resnum.append(rn)
            icode.append(l[26].strip())
            element.append(_element(nm, l[76:78] if len(l) > 77 else "", rnm))
    return Structure(xyz, name, resname, chain, resnum, element, icode)


def write_pdb(st, path):
    """A structure back out as a PDB, for the pose viewer and for a reader who
    wants to see the geometry the numbers were measured on."""
    with open(path, "w") as fh:
        last = None
        for i in range(len(st)):
            if last is not None and st.chain[i] != last:
                fh.write("TER\n")
            last = st.chain[i]
            nm = st.name[i]
            nm = (" %-3s" % nm) if len(nm) < 4 else nm
            fh.write("ATOM  %5d %s %-3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                     "          %2s\n"
                     % (i + 1, nm, st.resname[i], st.chain[i], st.resnum[i],
                        st.xyz[i][0], st.xyz[i][1], st.xyz[i][2],
                        st.element[i][:2]))
        fh.write("TER\nEND\n")
    return path


# ---------------------------------------------------------------------- SASA

def _sphere(n=POINTS):
    """n points spread evenly over the unit sphere, by the golden spiral."""
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = math.pi * (1.0 + 5.0 ** 0.5) * i
    return np.column_stack((np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi), np.cos(phi)))


def sasa(xyz, radii, probe=PROBE, points=POINTS):
    """Per atom solvent accessible surface area, in A^2 (Shrake and Rupley).

    Numerical, not analytical: each atom's expanded sphere carries `points`
    test points and the area is the share of them no other expanded sphere
    covers. Numpy over every neighbour at once, so a 2000 atom complex takes
    well under a second and the three states of a pose are free.
    """
    xyz = np.asarray(xyz, dtype=float)
    r = np.asarray(radii, dtype=float) + probe
    n = len(xyz)
    if n == 0:
        return np.zeros(0)
    unit = _sphere(points)
    out = np.zeros(n)
    # Neighbours found a block of atoms at a time. A full n by n matrix is
    # fine for a peptide complex and is 3 GB for a 20000 atom one, which is a
    # receptor this tool is expected to dock against.
    BLOCK = 256
    for lo in range(0, n, BLOCK):
        hi = min(lo + BLOCK, n)
        d = np.linalg.norm(xyz[lo:hi, None, :] - xyz[None, :, :], axis=-1)
        near_blk = d < (r[lo:hi, None] + r[None, :])
        near_blk[np.arange(hi - lo), np.arange(lo, hi)] = False
        for bi in range(hi - lo):
            i = lo + bi
            j = np.nonzero(near_blk[bi])[0]
            pts = xyz[i] + unit * r[i]
            if len(j):
                dd = np.linalg.norm(pts[:, None, :] - xyz[j][None, :, :],
                                    axis=-1)
                free = ~(dd < r[j][None, :]).any(axis=1)
            else:
                free = np.ones(len(pts), dtype=bool)
            out[i] = 4.0 * math.pi * r[i] ** 2 * free.sum() / len(pts)
    return out


def parents(st):
    """For each hydrogen, the index of the heavy atom it hangs off.

    A hydrogen's surface belongs to its side chain, not to a species of its
    own, so burial() folds each hydrogen's area onto the atom named here.
    """
    out = np.full(len(st), -1, dtype=int)
    heavy = np.nonzero(st.heavy)[0]
    if not len(heavy):
        return out
    # Grouped once. Walking every heavy atom for every hydrogen is quadratic
    # and was the slowest thing in the panel on a large receptor.
    by_res = {}
    for k in heavy:
        by_res.setdefault(st.resid(k), []).append(k)
    for i in np.nonzero(~st.heavy)[0]:
        # Its own residue and nowhere else. A hydrogen is never closer to the
        # next residue's carbon than to its own, and ff14SB always puts it in
        # the same residue, but a residue holding hydrogens ALONE, which a
        # stripped water or a mis-split heteroatom can produce, used to reach
        # across the structure and fold a ligand hydrogen's area onto a
        # receptor atom. Unparented is better: burial leaves the area where it
        # is and the total is still conserved.
        pool = by_res.get(st.resid(i))
        if not pool:
            continue
        pool = np.array(pool)
        d = np.linalg.norm(st.xyz[pool] - st.xyz[i], axis=1)
        out[i] = pool[int(np.argmin(d))]
    return out


def burial(complex_st, lig_chain="B", points=None):
    """Buried area, and what kind of surface it is.

    Three SASA evaluations on the one geometry: the complex, and each partner
    alone in the same coordinates. The partners are never moved, so the only
    difference between the states is the other molecule's presence.

    bsa is the total area both partners lose, which is the convention a
    structure paper reports; half of it is the interface area on one side.
    """
    lig = np.array([c == lig_chain for c in complex_st.chain], dtype=bool)
    rec = ~lig
    if not lig.any() or not rec.any():
        # Zeros and Nones rather than missing keys. A vanished column makes a
        # ragged row in the TSV, and a reader cannot tell "no interface" from
        # "not measured". This is the same shape the answer takes when the
        # peptide is simply far away.
        return dict(EMPTY_BURIAL)
    # Passed through rather than defaulted inside sasa, so a caller asking for
    # a different point count actually gets one.
    n_pts = int(points or POINTS)
    a_com = sasa(complex_st.xyz, complex_st.radii, points=n_pts)
    a_lig = np.zeros(len(complex_st))
    a_rec = np.zeros(len(complex_st))
    a_lig[lig] = sasa(complex_st.xyz[lig], complex_st.radii[lig],
                      points=n_pts)
    a_rec[rec] = sasa(complex_st.xyz[rec], complex_st.radii[rec],
                      points=n_pts)
    free = a_lig + a_rec
    lost = free - a_com
    # Hydrogen area folded onto the parent heavy atom, so "apolar" is a
    # property of a chemical group and not of whichever atom of it the
    # solvent happened to touch.
    par = parents(complex_st)
    folded = lost.copy()
    for i in np.nonzero(~complex_st.heavy)[0]:
        if par[i] >= 0:
            folded[par[i]] += folded[i]
            folded[i] = 0.0
    heavy = complex_st.heavy
    apolar = np.array([e in APOLAR_ELEM for e in complex_st.element])
    pep = lig & heavy
    pep_lost = folded[pep].sum()
    pep_apolar = folded[pep & apolar].sum()
    lig_free = a_lig[lig].sum()
    n_res = len({complex_st.resid(i) for i in np.nonzero(lig)[0]})
    out = {
        # WHICH ATOMS this area was rolled over, carried beside the area and
        # not inferred later. A complex that has been through tleap has its
        # hydrogens and a docked pose does not, and the same pose measured
        # both ways differs by 8 per cent in buried area and 0.05 in the
        # apolar share. The two are not comparable, so the basis travels with
        # the number.
        "bsa_basis": ("heavy atoms and hydrogens" if (~complex_st.heavy).any()
                      else "heavy atoms only"),
        # To the nearest A^2: see POINTS for why not to a tenth.
        "bsa": round(float(lost.sum())),
        "bsa_pep": round(float(pep_lost)),
        "bsa_per_res": round(float(pep_lost / n_res)) if n_res else None,
        "f_pep_buried": (round(float(1.0 - a_com[lig].sum() / lig_free), 3)
                         if lig_free > 0 else None),
        "f_hydrophobic": (round(float(pep_apolar / pep_lost), 3)
                          if pep_lost > 0 else None),
        "sasa_pep_free": round(float(lig_free)),
        "sasa_pep_bound": round(float(a_com[lig].sum())),
    }
    return out


# ------------------------------------------------------------------ contacts

#: What "no interface" looks like, so it looks the same however it arose.
EMPTY_BURIAL = {"bsa": 0.0, "bsa_pep": 0.0, "bsa_per_res": None,
                "bsa_basis": None,
                "f_pep_buried": None, "f_hydrophobic": None,
                "sasa_pep_free": None, "sasa_pep_bound": None}
EMPTY_CONTACTS = {"n_hbond": 0, "n_ionic": 0, "n_stack": 0, "n_contact": 0,
                  "n_clash": 0, "n_rec_res": 0, "n_pep_res": 0}


def _on(table, resn, atom):
    return atom in table.get(resn, ())


def contacts(st, lig_chain="B"):
    """Typed contacts across the interface, and the clashes.

    What is counted, and at which level, because the level is the whole
    meaning of the number:

      n_hbond    ATOM PAIRS of N or O within 3.5 A. A donor and acceptor
                 proximity test: one residue pair can make more than one
                 hydrogen bond, and each counts.
      n_ionic    RESIDUE PAIRS where a basic side chain tip is within 4.0 A
                 of an acidic one. A bidentate arginine is one salt bridge,
                 not four.
      n_stack    RESIDUE PAIRS of aromatic rings within 5.0 A, ring atom to
                 ring atom. Proximity, not geometry: an edge to face pair
                 and a parallel stack are not told apart.
      n_contact  RESIDUE PAIRS with any heavy atom pair within 4.0 A.
      n_clash    heavy atom pairs closer than 0.80 of the sum of their radii.

    Hydrogens are ignored in every one of these counts, because a run without
    MM-GBSA has none, and a count that changed meaning depending on which
    stages ran would be worse than one that is consistently conservative.
    Buried area cannot follow that rule, since a hydrogen genuinely occupies
    surface, so burial() reports the basis it used instead.
    """
    lig = np.array([c == lig_chain for c in st.chain], dtype=bool) & st.heavy
    rec = np.array([c != lig_chain for c in st.chain], dtype=bool) & st.heavy
    li, ri = np.nonzero(lig)[0], np.nonzero(rec)[0]
    if not len(li) or not len(ri):
        return dict(EMPTY_CONTACTS, per_residue={})
    d = np.linalg.norm(st.xyz[li][:, None, :] - st.xyz[ri][None, :, :],
                       axis=-1)
    n_hbond = 0
    pairs_ionic, pairs_stack, pairs_contact = set(), set(), set()
    n_clash = 0
    per_res = {}
    close = np.argwhere(d <= STACK_A)
    for a, b in close:
        i, k = int(li[a]), int(ri[b])
        dist = float(d[a, b])
        an, ar = st.name[i].upper(), st.resname[i]
        bn, br = st.name[k].upper(), st.resname[k]
        ae, be = st.element[i], st.element[k]
        rp = (st.key(i), st.key(k))
        kinds = set()
        if dist <= VDW_A:
            kinds.add("vdw")
            pairs_contact.add(rp)
        if dist <= HBOND_A and ae in POLAR and be in POLAR:
            kinds.add("hbond")
            n_hbond += 1
        if dist <= IONIC_A and ((_on(BASIC, ar, an) and _on(ACIDIC, br, bn))
                                or (_on(ACIDIC, ar, an)
                                    and _on(BASIC, br, bn))):
            kinds.add("ionic")
            pairs_ionic.add(rp)
        if dist <= STACK_A and _on(AROM, ar, an) and _on(AROM, br, bn):
            kinds.add("stack")
            pairs_stack.add(rp)
        if dist < CLASH_F * (st.radii[i] + st.radii[k]):
            kinds.add("clash")
            n_clash += 1
        if kinds:
            e = per_res.setdefault(st.key(k), {"kinds": set(), "d_min": 99.0,
                                               "partners": set()})
            e["kinds"] |= kinds
            e["d_min"] = min(e["d_min"], round(dist, 2))
            e["partners"].add(st.key(i))
    return {
        "n_hbond": n_hbond,
        "n_ionic": len(pairs_ionic),
        "n_stack": len(pairs_stack),
        "n_contact": len(pairs_contact),
        "n_clash": n_clash,
        "n_rec_res": len({b for _a, b in pairs_contact}),
        "n_pep_res": len({a for a, _b in pairs_contact}),
        "per_residue": per_res,
    }


# ------------------------------------------------- MM-GBSA output, read back

_DIFF = "Differences (Complex - Receptor - Ligand):"
# The value must end where the field ends. Without the lookahead, an
# exponent form such as -8.1825E+01 would be read as -8.1825, out by a
# factor of ten and with nothing to show it.
_TERM = re.compile(r"^(BOND|ANGLE|DIHED|VDWAALS|EEL|1-4 VDW|1-4 EEL|EGB|"
                   r"ESURF|DELTA G gas|DELTA G solv|DELTA TOTAL)\s+"
                   r"(-?\d+\.?\d*)(?=\s|$)")
_KEY = {"VDWAALS": "vdw", "EEL": "eel", "EGB": "egb", "ESURF": "esurf",
        "DELTA G gas": "gas", "DELTA G solv": "solv", "DELTA TOTAL": "total"}


_SURF = re.compile(r"\b(surften|surfoff)\s*=\s*(-?\d*\.?\d+)")


def surface_constants(mmpbsa_in):
    """(surften, surfoff) as the pose was scored with, not as we assume.

    They convert the ESURF difference into a buried area, so a run that set
    them would otherwise be reported against ours.
    """
    ften, foff = SURFTEN, 0.0
    if mmpbsa_in and os.path.isfile(mmpbsa_in):
        for m in _SURF.finditer(open(mmpbsa_in).read()):
            if m.group(1) == "surften":
                ften = float(m.group(2))
            else:
                foff = float(m.group(2))
    return ften, foff


def terms(dat, mmpbsa_in=None):
    """The four MM-GBSA components of the reported dG, from the output it was
    reported from.

    They are already in every .dat file a run has written: only DELTA TOTAL
    was ever read out of it. Reading the rest costs nothing and cannot
    disagree with the total, because it is the same calculation.

    Also the buried area the GB nonpolar term implies, for comparison with
    the one burial() rolls a probe over. The two are different surface
    definitions and will not match to the last A^2; a large disagreement
    means the geometry measured is not the geometry scored.
    """
    if not dat or not os.path.isfile(dat):
        return {}
    out, inside = {}, False
    for line in open(dat):
        s = line.strip()
        if s.startswith(_DIFF):
            inside = True
            continue
        if not inside:
            continue
        m = _TERM.match(s)
        if m and m.group(1) in _KEY:
            out[_KEY[m.group(1)]] = float(m.group(2))
        if s.startswith("DELTA TOTAL"):
            break
    out["terms_read"] = sorted(k for k in out if k != "terms_read")
    if "esurf" in out:
        # ESURF = surften * SASA + surfoff per state, and the difference is
        # complex minus receptor minus ligand, so surfoff does NOT cancel: it
        # is subtracted once. Harmless at its default of 0, wrong as a reason.
        ften, foff = surface_constants(mmpbsa_in)
        if ften:
            area = -(out["esurf"] + foff) / ften
            # A positive ESURF difference would mean the complex exposes more
            # surface than its parts, which is not a buried area.
            out["bsa_lcpo"] = round(area) if area > 0 else None
    return out


# tleap numbers the residues 1 to N in the order they appear in the file it
# reads and ignores the PDB's own numbers, and ambpdb writes that numbering
# back out with no chain identifier at all. Both the decomposition and the
# minimised geometry therefore speak in positions, while the site, the log and
# the receptor the user chose speak in chain and residue number. This is the
# dictionary between them, read from the file tleap was actually given.
def tleap_map(forleap_pdb):
    """{tleap index: (chain, residue number, residue name, insertion code)}.

    The insertion code is carried because it is part of the identity: a
    receptor with residues 54 and 54A has two residues there, and a map that
    drops the code makes them one key, silently losing a row of the
    decomposition and mis-answering whether the site was touched.
    """
    out, seen, n = {}, None, 0
    if not forleap_pdb or not os.path.isfile(forleap_pdb):
        return out
    for l in open(forleap_pdb):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        key = (l[21], l[22:27])
        if key != seen:
            n += 1
            seen = key
            try:
                out[n] = (l[21], int(l[22:26]), l[17:20].strip().upper(),
                          l[26].strip())
            except ValueError:
                continue
    return out


class Mismatch(Exception):
    """The minimised file and the file tleap read are not the same complex."""


def read_minimised(min_pdb, tmap):
    """The minimised complex, with its chains and residue numbers put back.

    ambpdb writes tleap's numbering and no chain column, so an unrelabelled
    file cannot be asked which atoms are the peptide, and a hotspot given as
    A:54 cannot be found in it.

    The labels come from position, so they are only right while the two files
    hold the same residues in the same order. One extra residue in either,
    from a cap, an ion or a water, shifts every label after it and every
    number downstream: on a test where one residue was prepended, the clash
    count went from 0 to 8 and the site coverage from 1.0 to 0.67, all of it
    plausible and all of it wrong. So the fit is checked rather than assumed,
    on the count and on every residue name, and a mismatch raises instead of
    returning a labelled structure nobody can tell is mislabelled.
    """
    st = read_pdb(min_pdb)
    if not tmap:
        return st
    idx, seen, n = [], None, 0
    for i in range(len(st)):
        k = (st.chain[i], st.resnum[i], st.icode[i])
        if k != seen:
            n += 1
            seen = k
        idx.append(n)
    if n != len(tmap):
        raise Mismatch("the minimised file holds %d residues and the file "
                       "tleap read holds %d, so the two are not the same "
                       "complex" % (n, len(tmap)))
    for i, k in enumerate(idx):
        got = tmap.get(k)
        if not got:
            raise Mismatch("residue %d of the minimised file is not in the "
                           "map" % k)
        if got[2] != st.resname[i]:
            raise Mismatch("residue %d is %s in the minimised file and %s in "
                           "the file tleap read" % (k, st.resname[i], got[2]))
    for i, k in enumerate(idx):
        got = tmap[k]
        st.chain[i], st.resnum[i], st.resname[i] = got[0], got[1], got[2]
        st.icode[i] = got[3] if len(got) > 3 else ""
    return st


# One row per residue, comma separated, the residue named "<resname> <tleap
# index>" and its location marked R or L. Six terms follow, each with an
# average, a standard deviation and a standard error, so the total's average
# is the third number from the end.
# Amber names an ion Na+ or Cl-, so the name is not letters alone.
_RES = re.compile(r"^\s*([A-Z][A-Za-z0-9+-]{1,3})\s+(\d+)\s*,\s*([RL])\b")


def decomposition(path, tmap=None):
    """MMPBSA.py's per residue decomposition, as {(chain, resnum): terms}.

    Why MMPBSA.py's own decomposition and not a second calculation: the per
    residue energies then come from the same topology, the same frame, the
    same GB model and the same radii as the number the run reports. A per
    residue energy from a different force field can be neither added to nor
    compared with the reported total.

    IT DOES NOT SUM TO THE REPORTED TOTAL, and it is worth knowing why before
    anyone presents a hotspot energy as a share of one. MMPBSA.py forces
    gbsa=2 for any decomposition, so the per residue nonpolar term is an ICOSA
    per atom surface, while the totals table keeps the LCPO surface. Measured
    over the five scored poses of one peptide, the decomposition sums to
    between 1.5 and 5.3 per cent more negative than the reported total, 0.5 to
    1.4 kcal/mol. The van der Waals, electrostatic and polar sums agree with
    the reported terms to within 0.05 kcal/mol in every one of them: the whole
    difference is the two surface models. So expect agreement to a few per
    cent, always in that direction, and read a site energy as a decomposition
    of that pose rather than as an exact slice of its total.

    How a residue's share is arrived at: van der Waals and electrostatics are
    pairwise and SUM exactly, but each pair is split half to each residue by
    convention, so no single residue's share is more exact than any other's.
    With idecomp=2 the 1-4 terms are folded into those two columns, which
    makes them a different quantity from the DELTA VDWAALS and DELTA EEL rows
    they sit beside. The GB polar and surface terms are apportioned per atom.
    A residue does not own a binding energy; this is a bookkeeping convention
    inside an approximate model.
    """
    if not path or not os.path.isfile(path):
        return {}
    out = {}
    for line in open(path):
        m = _RES.match(line)
        if not m:
            continue
        resname, num, side = m.group(1), int(m.group(2)), m.group(3)
        nums = []
        for b in line.split(",")[2:]:
            try:
                nums.append(float(b))
            except ValueError:
                continue
        # Six terms, each an average, a standard deviation and a standard
        # error: 18 numbers, as the file's own header declares. A shorter line
        # is a truncated or half written file, and reading its columns by
        # position would put a total where a surface term belongs.
        if len(nums) < 18:
            continue
        if tmap and num in tmap:
            t = tmap[num]
            ic = (t[3] if len(t) > 3 else "") or ""
            key = (t[0], ("%d%s" % (t[1], ic)) if ic else t[1])
        else:
            key = ("", num)
        out[key] = {
            "resname": resname, "side": side, "tleap": num,
            "total": round(nums[-3], 3),
            "internal": round(nums[0], 3), "vdw": round(nums[3], 3),
            "eel": round(nums[6], 3), "polar": round(nums[9], 3),
            "nonpolar": round(nums[12], 3),
        }
    return out


# ---------------------------------------------------------------------- site

def site_metrics(hotspots, per_residue, decomp=None):
    """What the nominated site contributes, and how much of it is touched.

    `hotspots` is [(chain, resnum)] as the run recorded them. Coverage is
    the share of those residues that the pose actually contacts, and it is
    quantised by the size of the site: on a site of four residues it can only
    ever be 0, 0.25, 0.5, 0.75 or 1.
    """
    if not hotspots:
        return {}
    touched = 0
    energies = []
    for ch, num in hotspots:
        hit = [k for k in (per_residue or {})
               if k[0] == ch and k[1] == num]
        if hit:
            touched += 1
        if decomp:
            d = decomp.get((ch, num)) or decomp.get(("", num))
            if d is not None:
                energies.append(d["total"])
    out = {"site_n": len(hotspots),
           "site_touched": touched,
           "site_coverage": round(touched / len(hotspots), 3)}
    if energies:
        out["site_dg"] = round(sum(energies), 2)
        out["site_dg_mean"] = round(sum(energies) / len(energies), 2)
        out["site_dg_n"] = len(energies)
    return out


# ------------------------------------------------------------------ sequence

# Apolar and aromatic side chains, plus cysteine. Used for the longest apolar
# run, which is a synthesis and aggregation risk rather than a binding one.
# NOT Kyte and Doolittle's positive set, which excludes W at -0.9 and Y at
# -1.3: a run of aromatics stalls solid phase coupling whatever its hydropathy
# index says.
APOLAR_AA = set("AVILMFWYC")

# A sequence holding a non-standard or D residue, in the notation the ADCP 1.1
# library uses: "&" before a D residue and "<NLE>" for a modified one. The
# authority is peptide.EXTENDED; this is the narrow test of whether any
# sequence-level number can be computed at all, kept here so this module does
# not reach back into the package around it.
_NONSTANDARD = re.compile(r"[<&>]")


def apolar_run(seq):
    """The longest run of consecutive apolar residues."""
    best = run = 0
    for a in (seq or "").upper():
        run = run + 1 if a in APOLAR_AA else 0
        best = max(best, run)
    return best


def developability(seq, ph=7.4):
    """The sequence panel, from HOPE-PE.

    Vendored unchanged from the pipelines' hope_core so the two tools report
    the same numbers for the same peptide. Sequence only: nothing here
    depends on the pose, and nothing here is defined for a non-standard or D
    residue, so a peptide holding one gets no developability row rather than
    a wrong one.

    That last sentence is enforced here rather than assumed. hope_pe keeps
    every standard letter it finds and drops the rest, so N<BHD>PIS&DVD, which
    is seven residues with one modified side chain, would otherwise come back
    as the nine letter NHDPISDVD: a hydropathy, a charge and a verdict for a
    peptide that does not exist, disagreeing with the length column beside it.
    """
    try:
        from . import hope_pe
    except ImportError:                                       # noqa: BLE001
        return {}
    if _NONSTANDARD.search(seq or ""):
        return {"pe_note": "no sequence panel: non-standard or D residues"}
    ev = hope_pe.evaluate(seq, ph=ph)
    if ev.get("pe_score") is None:
        return {}
    keep = ("gravy", "charge_ph7", "pi", "aromaticity", "aliphatic_index",
            "boman", "instability", "hydrophobic_moment", "mw",
            "pe_score", "verdict", "pe_flags")
    out = {k: ev[k] for k in keep if k in ev}
    out["apolar_run"] = apolar_run(hope_pe._clean(seq))
    out["liabilities"] = {k: ev[k] for k in hope_pe.liabilities("A")
                          if k in ev and ev[k]}
    return out


# --------------------------------------------------- getting to one geometry

DECOMP_IN = """Per residue decomposition of a pose already scored
&general
   startframe=1, endframe=1, interval=1, verbose=2, keep_files=0,
   receptor_mask="{rec_mask}", ligand_mask="{lig_mask}",
/
&gb
   igb=5, saltcon=0.150,
/
&decomp
   idecomp=2, dec_verbose=0, print_res="all",
/
"""

_MASK = re.compile(r'(receptor|ligand)_mask\s*=\s*"([^"]+)"')


def masks_of(mmpbsa_in):
    """The masks a pose was scored with, from the input file it kept.

    Rebuilding them here would be a second implementation of the rule that
    decides which residues are the receptor, and the two could disagree. The
    file that MMPBSA.py was actually given cannot.
    """
    got = {}
    if mmpbsa_in and os.path.isfile(mmpbsa_in):
        for m in _MASK.finditer(open(mmpbsa_in).read()):
            got[m.group(1)] = m.group(2)
    return got.get("receptor"), got.get("ligand")


def scored_rst(posedir, prefix):
    """(restart file, what it is) for the coordinates MMPBSA.py was given.

    Almost always the minimised complex. When sander writes no restart file,
    run_amber hands MMPBSA.py the complex as tleap built it instead, and the
    energy in the table is then a single point on an unminimised pose. The
    panel has to measure THAT geometry, or the columns describe a structure
    the energy was never computed on. run_amber records its choice beside the
    pose; the file order here is the fallback for a pose scored before it did.
    """
    told = os.path.join(posedir, "%s_scored_on.txt" % prefix)
    if os.path.isfile(told):
        name = open(told).read().strip()
        p = os.path.join(posedir, os.path.basename(name))
        if os.path.isfile(p):
            return p, ("minimised complex" if "_min.rst7" in name
                       else "complex as built, scored without minimisation")
    for name, what in (("%s_min.rst7" % prefix, "minimised complex"),
                       ("%s_com.rst7" % prefix,
                        "complex as built, scored without minimisation")):
        p = os.path.join(posedir, name)
        if os.path.isfile(p):
            return p, what
    return None, None


def minimised_pdb(posedir, prefix, log=print):
    """The coordinates the energy belongs to, as a PDB, written once and kept.

    Measuring the docked pose instead would describe a structure that was
    scored and then moved.
    """
    out = os.path.join(posedir, "%s_scored.pdb" % prefix)
    parm = os.path.join(posedir, "%s_com.parm7" % prefix)
    rst, what = scored_rst(posedir, prefix)
    if os.path.isfile(out) and os.path.getsize(out) > 0:
        return out, what
    if not (os.path.isfile(parm) and rst):
        return None, None
    from . import mmgbsa, prep
    home = mmgbsa.amber_home()
    tool = (os.path.join(home, "bin", "ambpdb") if home
            else prep.find_tool("ambpdb"))
    if not tool or not os.path.isfile(tool):
        return None, None
    import subprocess
    with open(out, "w") as fh:
        p = subprocess.run([tool, "-p", parm, "-c", rst], stdout=fh,
                           stderr=subprocess.PIPE, text=True)
    if p.returncode != 0 or os.path.getsize(out) == 0:
        log("    %s: ambpdb could not write the scored geometry" % prefix)
        os.path.isfile(out) and os.remove(out)
        return None, None
    return out, what


def decomp_file(posedir, prefix, log=print, make=True):
    """The per residue file for a pose, computed on demand if it is missing.

    A run scored before the panel existed has its topologies and its
    minimised frame but no decomposition, and recomputing the decomposition
    is seconds against the minutes the pose already cost. It is written
    beside the pose so it is computed once.
    """
    out = os.path.join(posedir, "%s_decomp.dat" % prefix)
    if os.path.isfile(out) and os.path.getsize(out) > 0:
        return out
    if not make:
        return None
    rec_mask, lig_mask = masks_of(os.path.join(posedir,
                                               "%s_mmpbsa.in" % prefix))
    # The coordinates the energy was computed on, which is not always the
    # minimised ones: a pose scored without minimisation gets a decomposition
    # of what it was actually scored on rather than no decomposition at all.
    rst, _what = scored_rst(posedir, prefix)
    need = [os.path.join(posedir, "%s_%s" % (prefix, s))
            for s in ("com.parm7", "rec.parm7", "lig.parm7")]
    if not (rec_mask and lig_mask and rst
            and all(os.path.isfile(p) for p in need)):
        return None
    from . import mmgbsa, prep
    home = mmgbsa.amber_home()
    tool = prep.find_tool("MMPBSA.py")
    if not (home and tool):
        return None
    import subprocess
    inp = os.path.join(posedir, "%s_decomp.in" % prefix)
    with open(inp, "w") as fh:
        fh.write(DECOMP_IN.format(rec_mask=rec_mask, lig_mask=lig_mask))
    log("    %s: per residue decomposition of the pose already scored"
        % prefix)
    subprocess.run([tool, "-O", "-i", os.path.basename(inp),
                    "-cp", "%s_com.parm7" % prefix,
                    "-rp", "%s_rec.parm7" % prefix,
                    "-lp", "%s_lig.parm7" % prefix,
                    "-y", os.path.basename(rst),
                    "-o", "%s_decomp_total.dat" % prefix,
                    "-do", os.path.basename(out)],
                   cwd=posedir, capture_output=True, text=True,
                   env=dict(os.environ, AMBERHOME=home))
    if os.path.isfile(out) and os.path.getsize(out) > 0:
        return out
    log("    %s: the decomposition produced no output" % prefix)
    return None


_NSTEP = re.compile(r"^\s*(\d+)\s+(-?\d\.\d+E[+-]\d+)\s+(\d\.\d+E[+-]\d+)\s+"
                    r"(\d\.\d+E[+-]\d+)")
MAXCYC_MSG = "Maximum number of minimization cycles reached"


def minimisation(min_out):
    """What the minimisation before the energy actually did.

    Nothing read this before, and the pose's energy depends on it entirely.
    sander stops at maxcyc rather than at a gradient in almost every case,
    which is expected for a restrained rescoring minimisation and is not the
    same thing as having converged. Reporting the final gradient and where it
    stopped lets a reader tell a relaxed pose from one that was still moving.
    """
    if not min_out or not os.path.isfile(min_out):
        return {}
    last, maxed = None, False
    for line in open(min_out):
        m = _NSTEP.match(line)
        if m:
            last = m
        if MAXCYC_MSG in line:
            maxed = True
    if last is None:
        return {"minimised": False}
    return {"minimised": True,
            "min_cycles": int(last.group(1)),
            "min_energy": round(float(last.group(2)), 2),
            "min_rms_grad": round(float(last.group(3)), 4),
            "min_stopped_at_maxcyc": maxed}


def from_atoms(rec_atoms, pose, lig_chain="B"):
    """A complex Structure from the receptor as read and a docked pose.

    The path for a run with no MM-GBSA: no Amber installed, or a peptide with
    non-standard or D residues, which has no protein force field topology and
    so can never have one. Heavy atoms only, and no minimisation, so the
    energy columns stay empty while the interface columns are still measured.
    """
    xyz, name, resname, chain, resnum, element = [], [], [], [], [], []
    for a in rec_atoms:
        xyz.append(tuple(a.xyz))
        name.append(a.name)
        resname.append((a.resname or "").upper())
        chain.append(getattr(a, "chain", "") or "A")
        resnum.append(int(a.resnum))
        element.append(_element(a.name, getattr(a, "element", ""),
                                getattr(a, "resname", "")))
    for i in range(len(pose["name"])):
        xyz.append(tuple(pose["xyz"][i]))
        name.append(pose["name"][i])
        resname.append((pose["resname"][i] or "").upper())
        chain.append(lig_chain)
        resnum.append(int(pose["resnum"][i]))
        # The docked pose carries no element column, so the residue name
        # decides: inside an amino acid, CA is a carbon.
        element.append(_element(pose["name"][i], "", pose["resname"][i]))
    return Structure(xyz, name, resname, chain, resnum, element)


# --------------------------------------------------------------- the panel

def for_structure(st, hotspots=None, lig_chain="B", dat=None, decomp=None,
                  mmpbsa_in=None):
    """Every geometric metric for one complex, plus the energy terms of the
    MM-GBSA output that belongs to it."""
    c = contacts(st, lig_chain=lig_chain)
    out = {}
    out.update({k: v for k, v in c.items() if k != "per_residue"})
    out.update(burial(st, lig_chain=lig_chain))
    t = terms(dat, mmpbsa_in) if dat else {}
    for k, v in t.items():
        # Named so a partial read is visible: a .dat holding NaN or an
        # overflowed field drops a term, and a missing column would otherwise
        # look like a term that is zero.
        out[("bsa_lcpo" if k == "bsa_lcpo" else "mm_" + k)] = v
    out.update(site_metrics(hotspots, c.get("per_residue") or {},
                            decomp or {}))
    out["has_hydrogens"] = bool((~st.heavy).any())
    return out, (c.get("per_residue") or {})


def for_pose(posedir, prefix, hotspots=None, lig_chain="B", log=print,
             make_decomp=True):
    """The panel for one scored pose, from the files that pose kept.

    Returns (metrics, per residue rows), or (None, []) when the pose has too
    little on disk to measure. Nothing is recomputed that is already there.
    """
    tmap = tleap_map(os.path.join(posedir, "%s_com_forleap.pdb" % prefix))
    if not tmap:
        tmap = tleap_map(os.path.join(posedir, "%s_com.pdb" % prefix))
    mp, what = minimised_pdb(posedir, prefix, log=log)
    if not mp or not tmap:
        return None, []
    st = read_minimised(mp, tmap)
    # The ligand is whichever chain the complex was split on, which
    # _write_pdbs stamped as B.
    dat = os.path.join(posedir, "%s_mmgbsa.dat" % prefix)
    dec = decomposition(decomp_file(posedir, prefix, log=log,
                                    make=make_decomp), tmap)
    m, per_res = for_structure(st, hotspots=hotspots, lig_chain=lig_chain,
                               dat=dat if os.path.isfile(dat) else None,
                               decomp=dec,
                               mmpbsa_in=os.path.join(
                                   posedir, "%s_mmpbsa.in" % prefix))
    # Which coordinates these numbers describe, and what the minimisation
    # that produced them did. Recorded rather than assumed: run_amber falls
    # back to the unminimised pose when sander writes no restart file, and
    # that pose's energy and its interface are both a different thing.
    mi = minimisation(os.path.join(posedir, "%s_min.out" % prefix))
    m.update(mi)
    m["geometry"] = what or "unknown"
    m["decomposed"] = bool(dec)
    return m, residue_table(per_res, dec, hotspots)


def residue_table(per_residue, decomp, hotspots=None):
    """The per residue map, as rows ready for a TSV.

    Its own file rather than a column: a receptor has as many rows as it has
    residues at the interface, and a table with one row per peptide has
    nowhere to put them.
    """
    hot = {(c, n) for c, n in (hotspots or [])}
    rows = []
    for (ch, num, resname), e in sorted(per_residue.items(),
                                        key=lambda kv: (kv[0][0], kv[0][1])):
        d = (decomp or {}).get((ch, num)) or (decomp or {}).get(("", num))
        rows.append({
            "chain": ch, "resnum": num, "resname": resname,
            "site": int((ch, num) in hot),
            "kinds": ",".join(sorted(e["kinds"])),
            "d_min": e["d_min"],
            "n_partners": len(e["partners"]),
            "dg": None if d is None else d["total"],
        })
    return rows
