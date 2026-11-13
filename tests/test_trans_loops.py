"""Rebuilt loops come out as proteins: every peptide bond trans, every residue
of the L hand, nothing clashing, and the same loop every time.

By Aadhil Haq

Lengths and backbone angles alone passed loops with cis and twisted peptide
bonds (7A3R chain A: 274-275 at -2 degrees, 194-195 at -95) and residues of the
D hand, and a cis bond that is not before a proline stays cis through any
simulation that follows. What is checked:
- junction_geometry flags a cis bond, a twisted one, a residue of the D hand,
  an allo Thr and a clash, and passes a measured stretch and a cis bond before
  a proline;
- a rebuilt residue of the D hand is set right before the relaxation;
- PDBFixer's own loop on the fragment of 7A3R chain A (residues 170-200,
  191-194 missing) is flagged, and the relaxed loop passes: every peptide bond
  in or beside it within 30 degrees of trans, C-N 1.30-1.37 A, CA-CA
  3.70-3.95 A, every rebuilt residue of the L hand, no contact under 2.5 A;
- the relaxation holds the bonds trans and the residues L while it runs, and
  says so;
- it is the same every time: two relaxations byte for byte, and the attempts
  run side by side give what one after another give;
- the atoms it simulates feel the forces they would in the whole structure;
- a loop that fails the checks on every seed stops the preparation.

    python3 tests/test_trans_loops.py
"""
import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_pilot._vendor import prepare  # noqa: E402

FIXTURE = os.path.join(HERE, "data", "7a3r_loop_and_tail.pdb")
fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


def atoms(path, keep=None):
    out = []
    for line in open(path):
        if line.startswith(("ATOM", "HETATM")) and (keep is None or int(line[22:26]) in keep):
            out.append([line[:30], [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                        line[54:].rstrip("\n")])
    return out


def write(rows, path):
    with open(path, "w") as fh:
        for head, xyz, tail in rows:
            fh.write(head + "%8.3f%8.3f%8.3f" % tuple(xyz) + tail + "\n")
    return path


def turned(rows, after, angle):
    """rows with every atom from residue `after` on turned by `angle` degrees
    about the C(after - 1) to N(after) bond: the twist of that bond changes by
    that much."""
    import numpy as np
    xyz = {(int(h[22:26]), h[12:16].strip()): np.array(p) for h, p, t in rows}
    c, n = xyz[after - 1, "C"], xyz[after, "N"]
    k = (n - c) / np.linalg.norm(n - c)
    a = math.radians(angle)
    out = []
    for h, p, t in rows:
        p = np.array(p)
        if int(h[22:26]) >= after:
            v = p - n
            p = n + v * math.cos(a) + np.cross(k, v) * math.sin(a) + k * np.dot(k, v) * (1 - math.cos(a))
        out.append([h, list(p), t])
    return out


def omegas(path, lo, hi):
    """{'A190-191': omega} for each bond with a residue in lo..hi on either side."""
    import numpy as np
    at, order = {}, []
    for h, p, t in atoms(path):
        r = int(h[22:26])
        if not order or order[-1] != r:
            order.append(r)
        at[r, h[12:16].strip()] = np.array(p)
    out = {}
    for a, b in zip(order, order[1:]):
        if lo <= a <= hi or lo <= b <= hi:
            if (a, "C") in at and (b, "N") in at:
                out["A%d-%d" % (a, b)] = prepare._dihedral(at[a, "CA"], at[a, "C"], at[b, "N"], at[b, "CA"])
    return out


work = tempfile.mkdtemp(prefix="test_trans_loops_")
built5 = [("A", str(n), "") for n in range(179, 184)]
frag = atoms(FIXTURE, set(range(176, 187)))

# ---- the checks, on geometry made to fail them: no PDBFixer needed ---------
ok = prepare.junction_geometry(write(frag, os.path.join(work, "measured.pdb")), built5)
check(ok == [], "a measured stretch (176-186) passes every check: %s" % ok)

cis = prepare.junction_geometry(write(turned(frag, 181, 180.0), os.path.join(work, "cis.pdb")), built5)
check(any(f.startswith("A180-181:") and "omega" in f for f in cis),
      "a cis peptide bond is flagged, by its twist: %s" % cis)
twist = prepare.junction_geometry(write(turned(frag, 181, 90.0), os.path.join(work, "twist.pdb")), built5)
check(any(f.startswith("A180-181:") for f in twist), "a bond twisted 90 degrees is flagged: %s" % twist)

pro = atoms(FIXTURE, set(range(184, 191)))
pro_cis = write(turned(pro, 187, 180.0), os.path.join(work, "cis_pro.pdb"))
w = omegas(pro_cis, 187, 187)["A186-187"]
got = prepare.junction_geometry(pro_cis, [("A", "187", "")])
check(abs(w) < 30 and not any(f.startswith("A186-187:") for f in got),
      "a cis bond before a proline passes (omega %.0f): %s" % (w, got))

mirror = [[h, [-p[0], p[1], p[2]], t] for h, p, t in frag]
d = prepare.junction_geometry(write(mirror, os.path.join(work, "mirror.pdb")), built5)
check(all(any(f.startswith("A%d " % n) and "of the D hand" in f for f in d) for n in (180, 181, 182, 183))
      and not any("A179 " in f for f in d),
      "in the mirror image every rebuilt residue with a CB is of the D hand, and the glycine has none: %s" % d)

allo = [list(r) for r in frag]
og1 = next(r for r in allo if r[0][17:20] == "THR" and int(r[0][22:26]) == 180 and r[0][12:16].strip() == "OG1")
cg2 = next(r for r in allo if r[0][17:20] == "THR" and int(r[0][22:26]) == 180 and r[0][12:16].strip() == "CG2")
og1[1], cg2[1] = cg2[1], og1[1]
a = prepare.junction_geometry(write(allo, os.path.join(work, "allo.pdb")), built5)
check(any("A180 THR" in f and "allo" in f for f in a), "an allo Thr is flagged: %s" % a)

clash = [list(r) for r in frag]
ca176 = next(r for r in clash if int(r[0][22:26]) == 176 and r[0][12:16].strip() == "CA")[1]
cd1 = next(r for r in clash if int(r[0][22:26]) == 181 and r[0][12:16].strip() == "CD1")
cd1[1] = [ca176[0] + 2.0, ca176[1], ca176[2]]
c = prepare.junction_geometry(write(clash, os.path.join(work, "clash.pdb")), built5)
check(any("A181 CD1" in f and "nearer than 2.5 A" in f for f in c), "a clash is flagged: %s" % c)

# ---- the rest needs PDBFixer and OpenMM ------------------------------------
try:
    import openmm  # noqa: F401
    import pdbfixer  # noqa: F401
except ImportError:
    if fails:
        print("FAIL test_trans_loops: %d of %d\n  " % (len(fails), ran) + "\n  ".join(fails))
        sys.exit(1)
    print("skip  PDBFixer or OpenMM is not installed here; the rest was checked (%d checks)" % ran)
    sys.exit(77)

# a rebuilt residue of the D hand is set right before the relaxation, CB where an L residue's would be
import numpy as np  # noqa: E402
from openmm import app, unit  # noqa: E402

mirrored = app.PDBFile(os.path.join(work, "mirror.pdb"))
xyz = np.array(mirrored.positions.value_in_unit(unit.nanometer))
check(prepare._l_hand(mirrored.topology, xyz, set(built5)) == 4, "the four residues of the D hand are set right")
with open(os.path.join(work, "righted.pdb"), "w") as fh:
    app.PDBFile.writeFile(mirrored.topology, xyz * 10.0, fh, keepIds=True)
r = prepare.junction_geometry(os.path.join(work, "righted.pdb"), built5)
check(not any(f.startswith(("A180 ", "A181 ", "A182 ", "A183 ")) for f in r),
      "set right, each is of the L hand with its CA undistorted: %s" % r)

steps, built = [], []
fixed = prepare.fix_structure(FIXTURE, os.path.join(work, "rebuilt.pdb"), steps, fill_loops=True,
                              built=built)
raw = prepare.junction_geometry(fixed, built)
check(raw, "PDBFixer's own loop is flagged before it is relaxed (it was not, so nothing is checked)")

relaxed = prepare.relax_built(fixed, os.path.join(work, "relaxed.pdb"), steps, built)
check(prepare.junction_geometry(relaxed, built) == [],
      "the relaxed loop passes every check: %s" % prepare.junction_geometry(relaxed, built))
w = omegas(relaxed, 191, 194)
check(len(w) == 5 and all(abs(x) > 150 for x in w.values()),
      "every peptide bond in or beside the loop is trans: %s" % w)
text = str(steps[-1])
check("peptide bond(s) held trans" in text and "of the L hand" in text,
      "the relaxation says it held the bonds trans and the residues L: %s" % text[-300:])
check("5 peptide bond(s) held trans" in text, "the five bonds in and beside the loop were held: %s" % text)

# the same every time, and side by side as one after another
again = prepare.relax_built(fixed, os.path.join(work, "again.pdb"), [], built)
check(open(again).read() == open(relaxed).read(), "two relaxations of one loop are the same, byte for byte")


real = prepare._attempts


class _InTurn(real):
    def __init__(self, folder, n):
        real.__init__(self, folder, n)
        self.cpus = 1


prepare._attempts = _InTurn
try:
    in_turn = prepare.relax_built(fixed, os.path.join(work, "in_turn.pdb"), [], built)
finally:
    prepare._attempts = real
check(open(in_turn).read() == open(relaxed).read(),
      "the attempts run in turn give what they give side by side")

# the atoms simulated feel what they would in the whole structure
pdb = app.PDBFile(fixed)
ff = app.ForceField("amber14-all.xml")
m = app.Modeller(pdb.topology, pdb.positions)
with prepare._seeded_random(1), prepare._one_cpu_thread() as cpu:
    m.addHydrogens(ff, pH=7.0, platform=cpu)
system = ff.createSystem(m.topology, nonbondedMethod=app.CutoffNonPeriodic,
                         nonbondedCutoff=1.2 * unit.nanometer, constraints=None)
made = set(built)
held = prepare._hold_chemistry(system, m.topology, made)
check(held == (5, 4), "five bonds held trans and four residues L: %s" % (held,))
moving = [a.index for a in m.topology.atoms()
          if (a.residue.chain.id, a.residue.id, "") in made]
whole = m.positions.value_in_unit(unit.nanometer)
sub, keep = prepare._neighbourhood(system, whole, moving, 0.0 + prepare.NEIGHBOURHOOD_NM)


def forces(s, xyz):
    ctx = openmm.Context(s, openmm.VerletIntegrator(0.001), openmm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(xyz)
    return ctx.getState(getForces=True).getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer)


f_all = forces(system, whole)
f_sub = forces(sub, [whole[i] for i in keep])
where = {old: new for new, old in enumerate(keep)}
err = max(np.linalg.norm(f_all[i] - f_sub[where[i]]) / max(1.0, np.linalg.norm(f_all[i])) for i in moving)
check(err < 1e-5, "the forces on the moving atoms are the whole structure's (worst %.1e)" % err)

# a loop that fails the checks on every seed stops the preparation
real_jg, real_schedule = prepare.junction_geometry, prepare.RELAX_SCHEDULE
prepare.junction_geometry = lambda path, b: ["A191-192: C-N 1.33 A, CA-C-N 116, C-N-CA 122, omega 2, "
                                             "CA-CA 2.93 A"]
prepare.RELAX_SCHEDULE = ((10.0, 0.05), (5.0, 0.05))
rdir = os.path.join(work, "run")
os.makedirs(rdir)
try:
    prepare.run(FIXTURE, rdir, fill_loops=True, protonation="reduce", reduce_exe=None,
                minimise_structure=False, log=lambda *a: None)
    check(False, "a loop with a cis bond on every seed stops the preparation")
except prepare.PrepError as e:
    check("omega 2" in str(e) and "Untick" in str(e),
          "the stop names the bond that failed and the setting: %s" % str(e)[-200:])
finally:
    prepare.junction_geometry, prepare.RELAX_SCHEDULE = real_jg, real_schedule
log_text = open(os.path.join(rdir, "preparation.log")).read()
check("building them again from seed 2" in log_text and "building them again from seed 3" in log_text,
      "each seed's loops were tried before it stopped")

if fails:
    print("FAIL test_trans_loops: %d of %d\n  " % (len(fails), ran) + "\n  ".join(fails))
    sys.exit(1)
print("ok    test_trans_loops: %d checks" % ran)
