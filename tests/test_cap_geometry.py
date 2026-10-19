#!/usr/bin/env python3
"""The caps put on a chain break: whole trans amides, built where an amide goes.

By Aadhil Haq

What is checked, on a fragment of 7A3R chain A (residues 170-200, with 191-194
missing as deposited, and 182-183 taken out here so that there are two breaks,
one of them two residues long):
- _dihedral measures a torsion in the IUPAC sense, the sense _place builds one
  in: 0 cis, 180 trans, positive clockwise looking down the middle bond. It is
  checked against torsions set out by hand, against what _place builds, and
  against the fragment itself, whose peptide bonds are all trans and whose
  L-amino acids have a negative phi almost everywhere;
- cap_breaks puts an NME on residues 181 and 190 and an ACE on 184 and 195,
  each joined by a peptide bond of 1.30-1.36 A with O-C-N at 120-125 degrees,
  omega trans, the amide hydrogen trans to the carbonyl oxygen, and no heavy
  atom of a cap within 2.5 A of a heavy atom that is neither bonded to it nor
  one bond beyond. Built extended, the ACE on 184 stood 1.4 A from that
  residue's own OD1;
- through the whole preparation (PDB2PQR, the caps, then the minimisation,
  which leaves the caps free to move) the same holds before the minimisation,
  each capped residue has one hydrogen on its nitrogen, and both caps come out
  of the minimisation as trans amides.

The torsion used to be measured as 180 degrees less than the one _place
builds, so an NME's nitrogen was placed by the mirror of its true torsion: on
7A3R as close as 0.65 A to the carbonyl oxygen, with O-C-N as small as 29
degrees. And PDB2PQR names the first of an N terminus's three hydrogens H,
which was kept as the amide hydrogen: 0.41 A from the ACE's carbon on 7A3R.

    python3 tests/test_cap_geometry.py
"""
import importlib.util
import math
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock._vendor import prepare  # noqa: E402

FIXTURE = os.path.join(HERE, "data", "7a3r_loop_and_tail.pdb")
CAPS = ("ACE", "NME")
fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


def wrap(x):
    return (x + 180.0) % 360.0 - 180.0


def torsion(a, b, c, d):
    """The a-b-c-d torsion, IUPAC's sign, worked out here rather than taken
    from the module under test."""
    b1, b2, b3 = b - a, c - b, d - c
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return math.degrees(math.atan2(float(np.linalg.norm(b2) * np.dot(b1, n2)),
                                   float(np.dot(n1, n2))))


def angle(a, b, c):
    u, v = a - b, c - b
    cos = float(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def atoms(path):
    out = []
    for line in open(path, errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            el = line[76:78].strip() or name.lstrip("0123456789")[:1]
            out.append({"name": name, "res": line[17:20].strip(),
                        "num": int(line[22:26]), "h": el == "H",
                        "xyz": np.array([float(line[30:38]), float(line[38:46]),
                                         float(line[46:54])])})
    return out


def residues(at):
    """{(number, name): {atom name: index}}, and the numbers in file order."""
    by, order = {}, []
    for i, x in enumerate(at):
        k = (x["num"], x["res"])
        if k not in by:
            order.append(k)
        by.setdefault(k, {})[x["name"]] = i
    return by, order


def caps(path, where, built=True):
    """Check every cap in `path`. Returns {"NME": n, "ACE": n}.

    As built (built=True) its angle and its room are checked too. After the
    minimisation those are the force field's, and what has to have survived
    is the peptide bond, trans."""
    at = atoms(path)
    xyz = np.array([x["xyz"] for x in at])
    heavy = np.array([not x["h"] for x in at])
    by, order = residues(at)
    cap = {i for i, x in enumerate(at) if x["res"] in CAPS}
    # bonds: by distance between the measured atoms, by template for the caps
    nb = {i: set() for i in range(len(at))}
    for i in range(len(at)):
        if i in cap:
            continue
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        for j in np.nonzero(d < (1.95 if heavy[i] else 1.25))[0]:
            j = int(j)
            if j != i and j not in cap and (heavy[i] or heavy[j]):
                nb[i].add(j)
                nb[j].add(i)
    seen = {"NME": 0, "ACE": 0}
    for n, k in enumerate(order):
        rn = k[1]
        if rn not in CAPS:
            continue
        seen[rn] += 1
        R = by[k]
        X = {nm: xyz[i] for nm, i in R.items()}
        methyl = R.get("CH3", R.get("C") if rn == "NME" else None)
        tag = "%s %d in %s" % (rn, k[0], where)
        if rn == "NME":
            P = by[order[n - 1]]
            bond = float(np.linalg.norm(xyz[P["C"]] - X["N"]))
            oca = angle(xyz[P["O"]], xyz[P["C"]], X["N"])
            omega = torsion(xyz[P["CA"]], xyz[P["C"]], X["N"], xyz[methyl])
            hn = (torsion(xyz[P["O"]], xyz[P["C"]], X["N"], X["H"])
                  if "H" in X else None)
            links = [(R["N"], P["C"]), (R["N"], methyl)] + (
                [(R["N"], R["H"])] if "H" in R else [])
        else:
            Q = by[order[n + 1]]
            bond = float(np.linalg.norm(X["C"] - xyz[Q["N"]]))
            oca = angle(X["O"], X["C"], xyz[Q["N"]])
            omega = torsion(xyz[methyl], X["C"], xyz[Q["N"]], xyz[Q["CA"]])
            hn = (torsion(X["O"], X["C"], xyz[Q["N"]], xyz[Q["H"]])
                  if "H" in Q else None)
            links = [(R["C"], Q["N"]), (R["C"], R["O"]), (R["C"], methyl)]
            on_n = [h for h in ("H", "H1", "H2", "H3") if h in Q]
            check(on_n == ["H"], "%s: the residue after it has %s on its "
                  "nitrogen, not one H" % (tag, on_n))
        check(1.30 <= bond <= 1.36, "%s: the peptide bond is %.3f A" % (tag, bond))
        if built:
            check(120.0 <= oca <= 125.0, "%s: O-C-N is %.1f degrees" % (tag, oca))
        check(abs(omega) > 150.0, "%s: omega is %.1f degrees, not trans" % (tag, omega))
        if hn is not None:
            check(abs(hn) > 150.0, "%s: the amide H is at %.1f degrees to the "
                  "carbonyl O, not trans to it" % (tag, hn))
        for nm, i in R.items():
            if nm.startswith("H") and nm != "H":
                links.append((methyl, i))
        for p, q in links:
            nb[p].add(q)
            nb[q].add(p)
    if built:
        for i in sorted(cap):
            if not heavy[i]:
                continue
            near = {i} | nb[i]
            for j in list(nb[i]):
                near |= nb[j]
            d = np.linalg.norm(xyz - xyz[i], axis=1)
            d[~heavy] = np.inf
            d[sorted(near)] = np.inf
            j = int(np.argmin(d))
            check(d[j] >= 2.5, "%s %d %s in %s is %.2f A from %s %d %s"
                  % (at[i]["res"], at[i]["num"], at[i]["name"], where, d[j],
                     at[j]["res"], at[j]["num"], at[j]["name"]))
    return seen


if not os.path.isfile(FIXTURE):
    print("skip  no fragment at %s" % FIXTURE)
    sys.exit(77)

# ---- the torsion convention ----------------------------------------------------
# Set out by hand: looking from b to c, a points right and d is turned t
# degrees clockwise from it, which IUPAC calls +t.
a, b, c = np.array([1.0, 0.0, 0.0]), np.zeros(3), np.array([0.0, 1.0, 0.0])
for t in (0.0, 60.0, 90.0, 120.0, 180.0, -60.0, -90.0, -150.0):
    r = math.radians(t)
    d = c + np.array([math.cos(r), 0.0, -math.sin(r)])
    got = prepare._dihedral(a, b, c, d)
    check(abs(wrap(got - t)) < 1e-6,
          "a torsion of %+.0f set out by hand measures %+.1f" % (t, got))
    check(abs(wrap(torsion(a, b, c, d) - t)) < 1e-6,
          "this test's own torsion reads %+.0f set out by hand as %+.1f"
          % (t, torsion(a, b, c, d)))
# what _place builds is what _dihedral measures
rng = np.random.default_rng(7)
for t in (-170.0, -100.0, -45.0, 0.0, 30.0, 116.0, 180.0):
    a, b, c = rng.normal(size=(3, 3)) * 1.5
    d = prepare._place(a, b, c, 1.33, 116.2, t)
    got = prepare._dihedral(a, b, c, d)
    check(abs(wrap(got - t)) < 1e-6,
          "_place built a torsion of %+.0f and _dihedral measures %+.1f" % (t, got))
# and the chain as measured: every peptide bond in the fragment is trans,
# which a torsion read as 180 less than it is calls cis, and an L-amino
# acid's phi is negative almost everywhere, which a torsion of the wrong
# sign reads as positive
at = atoms(FIXTURE)
res = {}
for x in at:
    res.setdefault(x["num"], {})[x["name"]] = x["xyz"]
omegas, phis = [], []
for p in sorted(res):
    q = p + 1
    if q not in res:
        continue
    omegas.append(prepare._dihedral(res[p]["CA"], res[p]["C"], res[q]["N"],
                                    res[q]["CA"]))
    if "CB" in res[q]:
        phis.append(prepare._dihedral(res[p]["C"], res[q]["N"], res[q]["CA"],
                                      res[q]["C"]))
cis = [round(w) for w in omegas if abs(w) < 150.0]
check(len(omegas) > 20 and not cis,
      "the fragment's peptide bonds are trans, and %d of %d measure cis: %s"
      % (len(cis), len(omegas), cis))
neg = sum(1 for x in phis if x < 0)
check(len(phis) > 20 and neg >= 0.8 * len(phis),
      "%d of the fragment's %d non-glycine phi angles measure negative"
      % (neg, len(phis)))

tmp = tempfile.mkdtemp(prefix="capgeo_", dir=os.environ.get("TMPDIR") or None)
try:
    # ---- the caps, as cap_breaks writes them --------------------------------------
    broken = os.path.join(tmp, "two_breaks.pdb")
    with open(FIXTURE) as src, open(broken, "w") as dst:
        for line in src:
            if line.startswith("ATOM") and int(line[22:26]) in (182, 183):
                continue
            dst.write(line)
    gaps = prepare.backbone_gaps(broken)
    check([(g[0][1], g[1][1]) for g in gaps] == [(181, 184), (190, 195)],
          "the fragment's breaks are 181-184 and 190-195: %s" % gaps)
    capped = prepare.cap_breaks(broken, os.path.join(tmp, "capped.pdb"), gaps, [],
                                log=lambda *x: None)
    check(capped is not None, "cap_breaks capped nothing")
    if capped:
        seen = caps(capped, "cap_breaks' output")
        check(seen == {"NME": 2, "ACE": 2}, "two NMEs and two ACEs: %s" % seen)

    # ---- and the whole preparation, minimisation included -------------------------
    missing = []
    if not prepare.find_tool("pdb2pqr", "pdb2pqr30"):
        missing.append("pdb2pqr")
    if importlib.util.find_spec("openmm") is None:
        missing.append("OpenMM")
    if missing:
        for f in fails:
            print("FAIL  " + f)
        if fails:
            sys.exit(1)
        print("skip  %d check(s) passed, but the preparation was not run: no %s "
              "here" % (ran, " or ".join(missing)))
        sys.exit(77)
    out = prepare.run(broken, os.path.join(tmp, "prep"), ph=7.4, chain="A",
                      fill_loops=False, minimise_structure=True,
                      protonation="propka", log=lambda *x: None)
    before = os.path.join(tmp, "prep", "capped.pdb")
    check(os.path.isfile(before), "the preparation capped nothing")
    if os.path.isfile(before):
        seen = caps(before, "the preparation, before minimisation")
        check(seen == {"NME": 2, "ACE": 2},
              "two NMEs and two ACEs before minimising: %s" % seen)
    plog = open(os.path.join(tmp, "prep", "preparation.log")).read()
    check("minimisation: done" in plog, "the minimisation did not run:\n" + plog)
    check(out is not None, "the preparation wrote nothing")
    if out:
        seen = caps(out, "the minimised receptor", built=False)
        check(seen == {"NME": 2, "ACE": 2},
              "two NMEs and two ACEs after minimising: %s" % seen)
finally:
    shutil.rmtree(tmp, ignore_errors=True)

for f in fails:
    print("FAIL  " + f)
print("%d check(s), %d failed" % (ran, len(fails)))
sys.exit(1 if fails else 0)
