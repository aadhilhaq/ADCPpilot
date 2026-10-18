"""Missing residues: which are built, how they are relaxed, and what stops a
preparation.

By Aadhil Haq

What is checked, on a fragment of 7A3R chain A (residues 170-200, with 191-194
missing as deposited, and SEQRES running on to 210):
- the four-residue internal gap is built when loops are filled, and nothing
  when they are not;
- the ten residues after the last one present are NOT built: PDBFixer reports
  them, being in SEQRES, and building them laid 7A3R's 38-residue tail out as a
  straight strand 177 A long with stretched bonds;
- the built loop is relaxed before protonation, every junction with the
  measured chain ends up with a sane C-N bond and backbone angles, and the
  relaxation says how far the stems moved;
- the minimisation leaves built residues free rather than holding them where
  PDBFixer guessed them;
- a structure's own hydrogens are taken off before protonation (an HXT stopped
  PDB2PQR outright), but a kept ligand's are left;
- junction_geometry flags a closure PDBFixer left at 42 degrees.

    python3 tests/test_loops.py
"""
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock._vendor import prepare  # noqa: E402

FIXTURE = os.path.join(HERE, "data", "7a3r_loop_and_tail.pdb")
fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


def residues(path):
    seen = []
    for line in open(path):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            seen.append(int(line[22:26]))
    return seen


# ---- hydrogens as given: no PDBFixer needed --------------------------------
lines = [
    "ATOM      1  N   LEU A 248      48.000   7.000  35.000  1.00  0.00           N\n",
    "ATOM      2  HXT LEU A 248      48.595   7.632  35.069  1.00  0.00           H\n",
    "ATOM      3 1HD2 LEU A 248      48.100   7.100  35.100  1.00  0.00           H\n",
    "HETATM    4  H1  LIG B 301      10.000  10.000  10.000  1.00  0.00           H\n",
    "HETATM    5  H1  HOH W 401      20.000  20.000  20.000  1.00  0.00           H\n",
]
steps = []
kept = prepare.drop_hydrogens(lines, steps)
names = [l[12:16].strip() for l in kept]
check("HXT" not in names and "1HD2" not in names,
      "the polymer's own hydrogens, an HXT and a 1HD2 among them, are taken off: %s" % names)
check(any(l.startswith("HETATM") and l[17:20] == "LIG" for l in kept),
      "a kept ligand's hydrogen stays")
check(not any(l[17:20] == "HOH" for l in kept), "a water's hydrogen goes")
check(steps and "removed 3" in str(steps[-1]), "the log says how many went: %s" % (steps[-1] if steps else ""))

# ---- junction_geometry on a closure left at 42 degrees ----------------------
bad = os.path.join(tempfile.mkdtemp(), "bad.pdb")
with open(bad, "w") as fh:
    # residue 1 measured; residue 2 built, its N folded back onto 1's CA
    fh.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")
    fh.write("ATOM      2  C   ALA A   1       1.520   0.000   0.000  1.00  0.00           C\n")
    fh.write("ATOM      3  N   ALA A   2       0.600   0.990   0.000  1.00  0.00           N\n")
    fh.write("ATOM      4  CA  ALA A   2       1.000   2.400   0.000  1.00  0.00           C\n")
out = prepare.junction_geometry(bad, [("A", "2", "")])
check(len(out) == 1 and "A1-2" in out[0], "a junction folded back on itself is flagged: %s" % out)
check(prepare.junction_geometry(bad, []) == [],
      "junctions between measured residues are not this check's business")

# ---- the rest needs PDBFixer and OpenMM ------------------------------------
try:
    import openmm  # noqa: F401
    import pdbfixer  # noqa: F401
except ImportError:
    if fails:
        print("FAIL test_loops: %d of %d\n  " % (len(fails), ran) + "\n  ".join(fails))
        sys.exit(1)
    print("skip  PDBFixer or OpenMM is not installed here; the rest was checked (%d checks)" % ran)
    sys.exit(77)

work = tempfile.mkdtemp(prefix="test_loops_")

# filled
steps, built = [], []
fixed = prepare.fix_structure(FIXTURE, os.path.join(work, "rebuilt.pdb"), steps,
                              fill_loops=True, built=built)
got = residues(fixed)
text = str(steps[-1])
check([b[1] for b in built] == ["191", "192", "193", "194"],
      "the four-residue internal gap is what is built: %s" % built)
check(max(got) == 200, "nothing after 200, the last residue present: up to %d" % max(got))
check("10 at the C terminus" in text and "not built" in text,
      "the log says the ten-residue tail was left out: %s" % text)

relaxed = prepare.relax_built(fixed, os.path.join(work, "relaxed.pdb"), steps, built)
check(prepare.junction_geometry(relaxed, built) == [],
      "every junction is in shape after relaxation: %s" % prepare.junction_geometry(relaxed, built))
text = str(steps[-1])
check(text.startswith("  rebuilt loops relaxed: done") and "stem residue(s)" in text
      and "A RMSD" in text, "the relaxation says how far the stems moved: %s" % text)

# not filled
steps2, built2 = [], []
prepare.fix_structure(FIXTURE, os.path.join(work, "asis.pdb"), steps2, fill_loops=False,
                      built=built2)
check(built2 == [], "nothing is built when loops are not filled: %s" % built2)
check("LEFT 4 missing residue(s)" in str(steps2[-1]),
      "the log says the internal gap was left alone: %s" % steps2[-1])

# the minimisation leaves built residues free. In a preparation it runs after
# protonation, so the fragment is given hydrogens first, by OpenMM's templates.
from openmm import app  # noqa: E402

pdb = app.PDBFile(relaxed)
mod = app.Modeller(pdb.topology, pdb.positions)
mod.addHydrogens(app.ForceField("amber14-all.xml"), pH=7.0)
protonated = os.path.join(work, "protonated.pdb")
with open(protonated, "w") as fh:
    app.PDBFile.writeFile(mod.topology, mod.positions, fh, keepIds=True)
steps3 = []
mini = prepare.minimise(protonated, os.path.join(work, "mini.pdb"), steps3, free=built)
text = str(steps3[-1])
check(mini is not None and "atom(s) of rebuilt loops left free" in text,
      "the minimisation does not hold the rebuilt loop: %s" % text)

# the same seed builds the same loop, atom for atom
def built_coords(path):
    return [line[30:54] for line in open(path) if line.startswith("ATOM") and 191 <= int(line[22:26]) <= 194]


again = []
prepare.fix_structure(FIXTURE, os.path.join(work, "again.pdb"), [], fill_loops=True, built=again)
check(built_coords(os.path.join(work, "again.pdb")) == built_coords(fixed),
      "the same seed builds the same loop, every coordinate")

# the message that stops a preparation names the tool's own setting
real_jg = prepare.junction_geometry
prepare.junction_geometry = lambda path, b: ["A190-191: C-N 1.80 A, CA-C-N 60, C-N-CA 60"]
prepare.LOOPS_SETTING, saved = "a tool's own name for it", prepare.LOOPS_SETTING
try:
    prepare.relax_built(fixed, os.path.join(work, "never.pdb"), [], built)
    check(False, "a loop that will not fit stops the preparation")
except prepare.PrepError as e:
    check("Untick 'a tool's own name for it'" in str(e),
          "the stop message names the tool's own setting: %s" % str(e)[-90:])
prepare.junction_geometry, prepare.LOOPS_SETTING = real_jg, saved

# run() builds the loops again from the next seed when they cannot be made to fit
real_relax, seeds_tried = prepare.relax_built, []


def failing_twice(path, out_path, steps, b, log=print):
    seeds_tried.append(len(seeds_tried) + 1)
    if len(seeds_tried) < 3:
        raise prepare.PrepError("did not fit")
    return path


prepare.relax_built = failing_twice
rdir = os.path.join(work, "run")
os.makedirs(rdir)
try:
    prepare.run(FIXTURE, rdir, fill_loops=True, protonation="reduce", reduce_exe=None,
                minimise_structure=False, log=lambda *a: None)
except prepare.PrepError:
    pass                     # no reduce here, so no hydrogens: what matters is the log
prepare.relax_built = real_relax
log_text = open(os.path.join(rdir, "preparation.log")).read()
check(seeds_tried == [1, 2, 3], "relaxation was tried on each seed's loops in turn: %s" % seeds_tried)
check("building them again from seed 2" in log_text and "building them again from seed 3" in log_text,
      "the log says the loops were built again, and from which seed")

if fails:
    print("FAIL test_loops: %d of %d\n  " % (len(fails), ran) + "\n  ".join(fails))
    sys.exit(1)
print("ok    test_loops: %d checks" % ran)
