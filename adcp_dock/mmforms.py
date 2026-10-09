"""Which peptides MM-GBSA can build, and which of their poses.

By Aadhil Haq

MM-GBSA builds its Amber topology from the peptide's sequence and lays the
docked coordinates on it. Three kinds of peptide go beyond a linear chain of
standard L residues, and each is handled here:

D RESIDUES ("&K") are built with the template of their L form, and the
docked coordinates carry the handedness. That is exact for this force field
and not an approximation: ff14SB has no term that tells a residue from its
mirror image. Its torsions all have phases of 0 or 180 degrees, so each is
an even function of the angle and reflecting the molecule leaves it
unchanged, and it has no CMAP grid, which is the one protein term that is
not even. Bonds, angles, van der Waals, electrostatics and the generalised
Born term depend on distances and angles only. Both claims are checked: the
torsion phases in the parameter files and in every topology built, and the
energy of a D peptide against its exact mirror image as an L peptide, which
agree to the last printed digit (tests/test_mmgbsa_forms.py). The one place
the L template does matter is where tleap puts the hydrogens it adds, which
it takes from the template's own geometry and so puts on the wrong side of a
D centre; mmgbsa.mirror_hydrogens places them as the mirror image of where
tleap puts them on the mirrored residue.

RINGS are rescored only where the docked pose actually closes them. ADCP's
ring is a restraint its search pulls on, not a bond, and a pose can leave the
two ends apart. Building such a pose with the closing bond would have the
minimisation drag the ends together, and the energy would describe a pose
that was never docked; building it without the bond would score a linear
peptide in place of the ring. Neither is done: a pose whose ends are further
apart than the cutoffs below is not rescored, and says so.

NON-STANDARD RESIDUES ("N<BHD>") are not rescored at all. ff14SB has no
parameters for their side chains.

A run keeps the rules it started with. A job file written before D residues
and rings could be rescored has no JOB_KEY, and its peptides are treated as
they were then: D residues and rings skipped, with a geometry panel only.
Otherwise a run resumed today would rescore some of its ring peptides and
not the ones that finished before, and rank them against each other.
"""
import json
import re

import numpy as np

from . import peptide

#: In the job file. True for every run made since D residues and closed
#: rings are rescored; a job file without it was made before, and keeps
#: skipping them.
JOB_KEY = "mmgbsa_d_rings"

#: How the rescoring describes each kind of peptide it does not score.
SKIP_NONSTANDARD = "non-standard residues"
SKIP_OLD_D = "non-standard or D residues"         # what older runs recorded
SKIP_OLD_RING = "cyclic peptide"                   # what older runs recorded

# ---- the closure cutoffs ----------------------------------------------------
#
# The closing bond is C to N (1.335 A in ff14SB) for a head to tail ring and
# S to S (2.038 A) for a disulfide. A pose is rescored when the two atoms are
# no further apart than these. Measured on 300 ring poses docked into MDM2
# (1YCR) on 2026-10-09: five peptides docked head to tail with ADCP 1.1 (20
# poses each) and 1.0 (10 each), and five with a cysteine at each end
# (mostly the same cores) docked as disulfides with both. Every pose up to 6 A apart was then
# built with its closing bond and minimised as MM-GBSA minimises, whatever
# the cutoff would have said.
#
# HEAD TO TAIL, 2.5 A. ADCP's ring is no help here: its ends sit anywhere
# from 0.39 to 4.65 A (median 2.04, quartiles 1.46 and 2.53), in one broad
# spread with no gap in it, and its restraint holds the end residues' alpha
# carbons about 2.2 A apart, where a peptide bond puts them 3.8 A apart. What
# decides it is what closing the bond does. Below 2.40 A not one pose of 100
# had a stereocentre change in the minimisation; from 2.40 to 3.0 A, 2 of 27
# had an end residue's alpha carbon turned into its mirror image, and above
# 3.0 A, 8 of 23. A carbon and a nitrogen that are not bonded touch at about
# 3.25 A, so at 2.5 A they are held by the restraint, not resting in contact.
# The bond closed to 1.33 to 1.36 A in every pose, its energy term under
# 0.3 kcal/mol, and the peptide moved 1.0 A RMSD at the median.
#
# DISULFIDE, 3.5 A. Here the poses fall in two groups. ADCP 1.1's sulfurs
# are 2.03 to 3.19 A apart or 3.94 A and more (to 22 A), with none between;
# ADCP 1.0's run from 2.13 to 3.54 A, then 4.25 A and more. Two sulfurs that
# are not bonded touch at about 3.6 A, twice sulfur's van der Waals radius,
# so a pair closer than 3.5 A is held by the restraint; and 3.5 A is where
# ADCP 1.1's own OpenMM minimisation adds the bond. Closing it was harmless
# at any distance up to 6 A (no stereocentre changed, 1.0 A RMSD at the
# median, the bond 2.02 to 2.06 A after), so it is the distribution that
# sets this one.
BACKBONE_CLOSED = 2.5
DISULFIDE_CLOSED = 3.5


def enabled(j):
    """Whether this run rescores D residues and closed rings."""
    return bool((j or {}).get(JOB_KEY))


def name(seq):
    """The folder and file prefix a peptide's MM-GBSA output is written
    under.

    The sequence itself for a standard peptide, which is what every run
    before this used. A D residue's "&" becomes "_": no standard sequence
    holds one, so a D peptide can never share a folder with an L one, even
    on a file system that ignores case, and tleap is not handed a file name
    with a shell character in it.
    """
    s = str(seq or "")
    if s.isalpha():
        return s
    return re.sub(r"[^A-Za-z0-9_]", "_", s.replace("&", "_"))


def letters(seq):
    """The standard one-letter sequence with the D marks taken off: what the
    topology is built from."""
    return "".join(r.lstrip("&") for r in peptide.residues(seq))


def d_positions(seq):
    """The 1-based positions of the D residues."""
    return [i + 1 for i, r in enumerate(peptide.residues(seq))
            if r.startswith("&")]


def cysteine_positions(seq):
    """The 1-based positions of the residues a disulfide can join."""
    return [i + 1 for i, r in enumerate(peptide.residues(seq))
            if r in ("C", "&C")]


def unscoreable(seq, j):
    """Why this peptide is not rescored at all, or None when it can be.

    A run made before JOB_KEY gives the reasons it gave then, word for word,
    so its result files read the same whether written before or after.
    """
    mode = str((j or {}).get("cyclic") or "none")
    if not enabled(j):
        if peptide.is_extended(seq):
            return SKIP_OLD_D
        if mode != "none":
            return SKIP_OLD_RING
        return None
    if peptide.nonstandard(seq) or "<" in str(seq):
        return SKIP_NONSTANDARD
    return None


def engine_refusal(engine, seq, mode):
    """Why this MM-GBSA engine cannot score this peptide, or None.

    AmberTools and OpenMM read the topology tleap builds, which carries the
    closing bond and the D residues' coordinates, so both can. Prime is
    given the prepared complex as a PDB and builds its own: nothing here has
    shown that it keeps a ring's closing bond or leaves a D centre as it is,
    so it refuses rather than return a number for some other molecule.
    """
    mode = str(mode or "none")
    if str(engine or "amber").lower() != "schrodinger":
        return None
    d, ring = bool(d_positions(seq)), mode != "none"
    if not (d or ring):
        return None
    what = ("D residues and a ring" if d and ring
            else "D residues" if d else "a ring")
    return ("the Prime engine does not score a peptide with %s: it builds "
            "its own structure from a PDB file, and nothing has shown that "
            "it keeps them as they were docked" % what)


# ---- whether a pose closes its ring ---------------------------------------

def _atoms_by_residue(pose):
    """{residue number: {atom name: xyz}} from a pose read by adcp.read_pose."""
    out = {}
    for i, nm in enumerate(pose["name"]):
        out.setdefault(int(pose["resnum"][i]), {})[str(nm).strip()] = \
            np.asarray(pose["xyz"][i], dtype=float)
    return out


def closure(pose, seq, mode):
    """Whether this pose closes the ring it was docked as.

    Returns {"kind", "closed", "gap", "pairs", "why"}. "gap" is the distance,
    in A, between the atoms the closing bond joins: the first residue's N and
    the last residue's C for a head to tail ring, and the two SG atoms for a
    disulfide. "pairs" are the residues to join, as 1-based positions in the
    peptide: (1, n) for a head to tail ring, one (i, j) per disulfide. "why"
    is the reason a pose that is not closed gives, and None for one that is.

    With more than two cysteines, the closest pair under the cutoff is joined
    first, then the closest of what is left, and so on: ADCP pulls every pair
    of cysteines together, and each sulfur can make one bond.
    """
    mode = str(mode or "none")
    if mode == "none":
        return {"kind": "none", "closed": True, "gap": None, "pairs": [],
                "why": None}
    res = _atoms_by_residue(pose)
    order = sorted(res)
    if mode == "backbone":
        first, last = res[order[0]], res[order[-1]]
        if "N" not in first or "C" not in last:
            return {"kind": mode, "closed": False, "gap": None, "pairs": [],
                    "why": "ring not closed in this pose: its end atoms are "
                           "missing"}
        gap = float(np.linalg.norm(first["N"] - last["C"]))
        ok = gap <= BACKBONE_CLOSED
        return {"kind": mode, "closed": ok, "gap": round(gap, 2),
                "pairs": [[1, len(order)]] if ok else [],
                "why": None if ok else
                "ring not closed in this pose: ends %.2f A apart" % gap}
    # a disulfide
    pos = cysteine_positions(seq)
    where = {p: order[p - 1] for p in pos if 0 < p <= len(order)}
    sg = {p: res[n]["SG"] for p, n in where.items() if "SG" in res[n]}
    pairs = []
    for a in sorted(sg):
        for b in sorted(sg):
            if a < b:
                pairs.append((float(np.linalg.norm(sg[a] - sg[b])), a, b))
    pairs.sort()
    if not pairs:
        return {"kind": mode, "closed": False, "gap": None, "pairs": [],
                "why": "ring not closed in this pose: it has no two cysteine "
                       "sulfurs"}
    used, joined = set(), []
    for d, a, b in pairs:
        if d <= DISULFIDE_CLOSED and a not in used and b not in used:
            joined.append([a, b])
            used.update((a, b))
    gap = pairs[0][0]
    ok = bool(joined)
    return {"kind": mode, "closed": ok, "gap": round(gap, 2),
            "pairs": joined,
            "why": None if ok else
            "ring not closed in this pose: sulfurs %.2f A apart" % gap}


def none_closed(n):
    """The peptide's reason when not one of its top n poses closes the
    ring."""
    if n == 1:
        return "ring not closed in its top pose"
    return "ring not closed in any of its top %d poses" % n


# ---- what the topology builder is told ------------------------------------

# The properties a pose's molecule carries to mmgbsa, which builds from them.
from ._vendor.mmgbsa import D_PROP, RING_PROP  # noqa: E402


def mark(mol, seq, ring):
    """Attach to a transplanted pose what the topology has to be built with:
    the closing bonds and the D positions."""
    if mol is None:
        return None
    if ring and ring.get("kind") in ("backbone", "disulfide") \
            and ring.get("pairs"):
        mol.SetProp(RING_PROP, json.dumps({"kind": ring["kind"],
                                           "pairs": ring["pairs"],
                                           "gap": ring.get("gap")}))
    d = d_positions(seq)
    if d:
        mol.SetProp(D_PROP, ",".join(str(i) for i in d))
    return mol
