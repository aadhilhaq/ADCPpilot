"""Reading the Amber files MM-GBSA leaves, for the checks made on them.

By Aadhil Haq

Only what the checks need: the atom and residue tables and the bond list of
a topology (parm7), and the coordinates of a restart file (rst7). Text in,
numbers out, no Amber program run, so the checks cost nothing and can be
made on any run whose files are on disk.

What the checks are for. A ring peptide is rescored with a bond tleap was
told to make, and a D residue with the template of its L form, so two
things have to be shown on every such pose rather than assumed: that the
bond is in the topology and is a bond after the minimisation (its length
and its energy term, before and after), and that every stereocentre of the
peptide has the handedness the docking gave it once the minimisation is
done (stereo_centres, handedness).
"""
from __future__ import annotations

import math
import re

import numpy as np


def read_parm7(path, flags=("ATOM_NAME", "RESIDUE_LABEL", "RESIDUE_POINTER",
                            "BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN",
                            "BOND_FORCE_CONSTANT", "BOND_EQUIL_VALUE",
                            "DIHEDRAL_PHASE", "POINTERS", "CHARGE")):
    """{flag: [values]} for the sections asked for, from a parm7 file.

    Every format the sections use is fixed width (20a4, 10I8, 5E16.8), and
    each value is read by its width, never by splitting on spaces: atom names
    such as "HD21" fill the field with no space either side.
    """
    want = set(flags)
    out, cur, width, kind = {}, None, None, None
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith("%FLAG"):
                name = line.split()[1]
                cur = name if name in want else None
                if cur:
                    out[cur] = []
                width = None
                continue
            if line.startswith("%FORMAT"):
                m = re.search(r"\((\d+)([aIE])(\d+)", line)
                if m:
                    kind, width = m.group(2), int(m.group(3))
                continue
            if cur is None or width is None or line.startswith("%"):
                continue
            text = line.rstrip("\n")
            for k in range(0, len(text), width):
                field = text[k:k + width]
                if not field.strip() and kind != "a":
                    continue
                if kind == "a":
                    out[cur].append(field.strip())
                elif kind == "I":
                    out[cur].append(int(field))
                else:
                    out[cur].append(float(field))
    if "ATOM_NAME" in out and "POINTERS" in out:
        out["ATOM_NAME"] = out["ATOM_NAME"][:out["POINTERS"][0]]
    if "RESIDUE_LABEL" in out and "POINTERS" in out:
        out["RESIDUE_LABEL"] = out["RESIDUE_LABEL"][:out["POINTERS"][11]]
    return out


def read_rst7(path):
    """The coordinates of a restart or inpcrd file, as an (n, 3) array.

    Either form sander writes: the text one tleap writes, or the NetCDF one
    sander writes by default for its own restart (it starts "CDF")."""
    with open(path, "rb") as fh:
        magic = fh.read(3)
    if magic == b"CDF":
        from scipy.io import netcdf_file
        with netcdf_file(path, "r", mmap=False) as nc:
            return np.array(nc.variables["coordinates"][:], dtype=float
                            ).reshape(-1, 3)
    with open(path, errors="replace") as fh:
        fh.readline()
        n = int(fh.readline().split()[0])
        vals = []
        for line in fh:
            line = line.rstrip("\n")
            for k in range(0, len(line), 12):
                field = line[k:k + 12]
                if field.strip():
                    vals.append(float(field))
            if len(vals) >= 3 * n:
                break
    return np.array(vals[:3 * n]).reshape(n, 3)


def write_rst7(path, xyz, title="coordinates", like=None):
    """An inpcrd file sander reads, for coordinates computed here. With
    `like`, the title and count lines are copied from that file, so the
    header is exactly the one tleap wrote."""
    xyz = np.asarray(xyz, dtype=float)
    head = None
    if like:
        with open(like, errors="replace") as fh:
            head = [fh.readline(), fh.readline()]
    with open(path, "w") as fh:
        if head:
            fh.writelines(head)
        else:
            fh.write(title[:80] + "\n")
            fh.write("%6d\n" % len(xyz))
        flat = xyz.reshape(-1)
        for k in range(0, len(flat), 6):
            fh.write("".join("%12.7f" % v for v in flat[k:k + 6]) + "\n")


def residues(parm):
    """[(index from 1, label, {atom name: atom index from 0})]."""
    names = parm["ATOM_NAME"]
    ptr = list(parm["RESIDUE_POINTER"]) + [len(names) + 1]
    out = []
    for r, label in enumerate(parm["RESIDUE_LABEL"]):
        atoms = {names[i]: i for i in range(ptr[r] - 1, ptr[r + 1] - 1)}
        out.append((r + 1, label, atoms))
    return out


def bonds(parm):
    """[(atom i, atom j, k, r0)] for every bond, atoms indexed from 0."""
    k = parm.get("BOND_FORCE_CONSTANT") or []
    r0 = parm.get("BOND_EQUIL_VALUE") or []
    out = []
    for flag in ("BONDS_INC_HYDROGEN", "BONDS_WITHOUT_HYDROGEN"):
        v = parm.get(flag) or []
        for n in range(0, len(v) - 2, 3):
            t = v[n + 2] - 1
            out.append((v[n] // 3, v[n + 1] // 3, k[t], r0[t]))
    return out


def bond_term(parm, xyz, i, j):
    """{"present", "length", "k", "r0", "energy"} for the bond i-j.

    "present" is whether the topology holds the bond at all: a pair the
    force field does not join is two atoms pressed together, and its length
    alone would not show the difference. The energy is k (r - r0)^2, Amber's
    harmonic bond term, in kcal/mol.
    """
    r = float(np.linalg.norm(np.asarray(xyz[i]) - np.asarray(xyz[j])))
    for a, b, k, r0 in bonds(parm):
        if {a, b} == {i, j}:
            return {"present": True, "length": round(r, 3), "k": k,
                    "r0": r0, "energy": round(k * (r - r0) ** 2, 3)}
    return {"present": False, "length": round(r, 3), "k": None, "r0": None,
            "energy": None}


def signed_volume(a, centre, b, c):
    """(a - centre) . ((b - centre) x (c - centre)): its sign is the
    handedness of the three substituents seen from the centre."""
    centre = np.asarray(centre)
    return float(np.dot(np.asarray(a) - centre,
                        np.cross(np.asarray(b) - centre,
                                 np.asarray(c) - centre)))


# The stereocentres of the twenty amino acids, as (centre, three heavy atom
# substituents, the hydrogen on it). The alpha carbon of every residue but
# glycine, and the beta carbon of threonine and isoleucine. For an L residue
# the alpha volume (N, C, CB) is positive.
CENTRES = {"CA": ("N", "C", "CB", "HA")}
BETA = {"THR": ("CB", "CA", "OG1", "CG2", "HB"),
        "ILE": ("CB", "CA", "CG1", "CG2", "HB")}
_PROTEIN = {"ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYX", "CYM", "GLN",
            "GLU", "GLH", "HID", "HIE", "HIP", "HIS", "ILE", "LEU", "LYS",
            "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def _base(label):
    """The residue a label names, without tleap's terminal prefix."""
    label = label.strip().upper()
    if len(label) == 4 and label[0] in "NC" and label[1:] in _PROTEIN:
        return label[1:]
    return label


def stereo_centres(parm, residue_indices):
    """[(residue index, residue, centre, (a, b, c, h) atom indices)] for the
    stereocentres of the given residues, from the topology alone."""
    out = []
    table = {r: (lab, atoms) for r, lab, atoms in residues(parm)}
    for r in residue_indices:
        lab, atoms = table.get(r, ("", {}))
        base = _base(lab)
        specs = []
        if all(n in atoms for n in ("N", "CA", "C", "CB")):
            specs.append(("CA", "N", "C", "CB", "HA"))
        if base in BETA:
            specs.append(BETA[base])
        for centre, a, b, c, h in specs:
            if all(n in atoms for n in (centre, a, b, c)):
                out.append((r, base, centre,
                            (atoms[centre], atoms[a], atoms[b], atoms[c],
                             atoms.get(h))))
    return out


def handedness(centres, xyz):
    """[{"residue", "name", "centre", "volume", "h_opposite"}] at these
    coordinates. "volume" is the signed volume of the three heavy
    substituents, "h_opposite" whether the hydrogen on the centre sits on
    the side away from them, which is where a tetrahedral centre has it."""
    out = []
    for r, base, centre, (ic, ia, ib, iq, ih) in centres:
        v = signed_volume(xyz[ia], xyz[ic], xyz[ib], xyz[iq])
        h_ok = None
        if ih is not None:
            mid = (np.asarray(xyz[ia]) + xyz[ib] + xyz[iq]) / 3.0
            h_ok = bool(np.dot(np.asarray(xyz[ih]) - xyz[ic],
                               mid - np.asarray(xyz[ic])) < 0)
        out.append({"residue": r, "name": base, "centre": centre,
                    "volume": round(v, 3), "h_opposite": h_ok})
    return out


def dihedral_phases(parm):
    """Every distinct torsion phase in the topology, in degrees."""
    return sorted({round(math.degrees(p), 3) % 360.0
                   for p in parm.get("DIHEDRAL_PHASE") or []})
