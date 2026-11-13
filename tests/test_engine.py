"""The parts that fail silently: target contents and the MM-GBSA split.

By Aadhil Haq

Both produce numbers when they are wrong. A stage that crashes announces
itself; a stage that returns a plausible energy computed on the wrong
molecules does not.

The last section checks that the stages call the shared engine in
adcp_pilot/_vendor rather than local reimplementations, which would
require the same fixes a second time.
"""
import json
import os
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from adcp_pilot import peptide, report, structure, target, verify   # noqa: E402
from adcp_pilot._vendor import mmgbsa, poseviewer                   # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


WORK = os.path.join(HERE, "_work")
os.makedirs(WORK, exist_ok=True)


def fake_trg(path, prefix, n_maps=30, rec_atoms=100):
    with zipfile.ZipFile(path, "w") as z:
        for i in range(n_maps):
            z.writestr("%s/rigidReceptor.t%d.map" % (prefix, i), "grid\n")
        z.writestr("%s/rigidReceptor.pdbqt" % prefix,
                   "".join("ATOM  %5d  CA  ALA A%4d       0.000   0.000"
                           "   0.000  0.00  0.00     0.000 C\n" % (i, i)
                           for i in range(rec_atoms)))
        z.writestr("%s/%s.cfg" % (prefix, prefix),
                   "[agfrCmdOptions]\nboxMode = user 1.0 2.0 3.0 30.0 30.0 30.0\n")


rec = os.path.join(WORK, "rec.pdbqt")
with open(rec, "w") as fh:
    for i in range(100):
        fh.write("ATOM  %5d  CA  ALA A%4d       0.000   0.000   0.000"
                 "  0.00  0.00     0.000 C\n" % (i, i))

# --- a good target passes --------------------------------------------------
good = os.path.join(WORK, "receptor.trg")
fake_trg(good, "receptor")
check(target.check(good, rec) is None,
      "a well-formed target was rejected: %s" % target.check(good, rec))
check(target.box_of(good) is not None, "the box could not be read back")

# --- and every way one goes wrong is caught --------------------------------
wrong_name = os.path.join(WORK, "receptor2.trg")
fake_trg(wrong_name, "somethingelse")
why = target.check(wrong_name, rec)
check(why and "folder inside" in why,
      "a target whose internal folder does not match its name was accepted "
      "(this is the one ADCP reports only as 'cannot locate map file')")

few = os.path.join(WORK, "few.trg")
fake_trg(few, "few", n_maps=3)
check(target.check(few, rec), "a target with three maps was accepted")

mismatch = os.path.join(WORK, "mismatch.trg")
fake_trg(mismatch, "mismatch", rec_atoms=55)
why = target.check(mismatch, rec)
check(why and "different structure" in why,
      "a target built from another receptor was accepted")

trunc = os.path.join(WORK, "trunc.trg")
with open(good, "rb") as a, open(trunc, "wb") as b:
    b.write(a.read(200))
check(target.check(trunc, rec), "a truncated target was accepted")
check(target.check(os.path.join(WORK, "nope.trg"), rec) == "does not exist",
      "a missing target was not reported as missing")

# --- the MM-GBSA masks come from the chains, in file order ----------------
com = os.path.join(WORK, "com.pdb")
with open(com, "w") as fh:
    # the peptide FIRST, which is what PDB2QR's chain sorting produces when
    # the receptor's letter sorts after the ligand's
    for i in range(1, 6):
        fh.write("ATOM  %5d  CA  GLY B%4d       0.000   0.000   0.000\n" % (i, i))
    for i in range(1, 21):
        fh.write("ATOM  %5d  CA  ALA E%4d       0.000   0.000   0.000\n"
                 % (i + 5, i))
rmask, lmask, nrec, nlig, how = mmgbsa.masks_by_chain(com, "B")
check(lmask == ":1-5",
      "the ligand mask is %r; the peptide is first in this file, so it is "
      "residues 1-5 and the old 'receptor then the rest' rule would have "
      "said :21-25" % lmask)
check(rmask == ":6-25", "the receptor mask is %r, expected :6-25" % rmask)
check((nrec, nlig) == (20, 5), "counted %d/%d residues" % (nrec, nlig))
# the ASL is what goes in the report, so it has to name the two chains the
# masks were actually built from
check("chain.name E" in how and "chain.name B" in how,
      "the description does not name the chains: %r" % how)

# --- verify() catches a split in the wrong place --------------------------
goodout = os.path.join(WORK, "good.dat")
open(goodout, "w").write(
    "Differences (Complex, Receptor, Ligand):\n"
    "BOND    0.0000\nANGLE   0.0000\nDIHED  -0.0000\n"
    "1-4 VDW 0.0000\n1-4 EEL 0.0000\nDELTA TOTAL  -42.1234\n")
ok, bad = verify.mmgbsa_split(goodout)
check(ok and not bad, "a clean result was reported as badly split")

badout = os.path.join(WORK, "bad.dat")
open(badout, "w").write(
    "Differences (Complex, Receptor, Ligand):\n"
    "BOND    0.0018\nANGLE   0.0576\nDIHED   7.3298\n"
    "1-4 VDW 0.9903\n1-4 EEL 36.8830\nDELTA TOTAL  -22.0424\n")
ok, bad = verify.mmgbsa_split(badout)
check(not ok, "the bonded terms from a wrongly split complex were accepted")
check(len(bad) >= 3, "only %d bonded term(s) flagged" % len(bad))

# --- the box clamps --------------------------------------------------------
atoms = [structure.Atom(name="CA", resname="ALA", chain="A", resnum=i,
                        x=float(i), y=0.0, z=0.0, element="C")
         for i in range(1, 4)]
_, size, _ = structure.box_for(atoms, [("A", 1), ("A", 2)], minimum=20.0,
                               maximum=30.0)
check(size[0] == 20.0, "two adjacent residues gave a %s box, not the floor"
      % size[0])
big = [structure.Atom(name="CA", resname="ALA", chain="A", resnum=i,
                      x=float(i * 20), y=0.0, z=0.0, element="C")
       for i in range(1, 6)]
_, size, _ = structure.box_for(big, [("A", 1), ("A", 5)], maximum=30.0)
check(size[0] == 30.0, "residues 80 A apart gave a %s box, not the ceiling"
      % size[0])

# --- the box drawn is the box docked into ---------------------------------
# The page draws structure.box_for; the run docks into receptor.box. One
# rule for both: along each axis, the span plus padding, clamped.
import types                                                    # noqa: E402
from adcp_pilot._vendor import receptor as vrec                   # noqa: E402
flat = [structure.Atom(name="CA", resname="ALA", chain="A", resnum=i,
                       x=x, y=y, z=z, element="C")
        for i, (x, y, z) in enumerate([(0.0, 0.0, 0.0), (10.0, 30.0, 2.0),
                                       (4.0, 12.0, 1.0)], 1)]
keys = [("A", 1), ("A", 2), ("A", 3)]
c1, s1, _ = structure.box_for(flat, keys, padding=8.0, minimum=20.0,
                              maximum=30.0)
stub = types.SimpleNamespace(
    hotspot_keys=keys, hotspots=[],
    atoms_of=lambda k: [types.SimpleNamespace(xyz=(a.x, a.y, a.z))
                        for a in flat])
c2, s2 = vrec.box(stub, types.SimpleNamespace(box_padding=8.0, box_min=20.0,
                                              box_max=30.0))
check([round(v, 6) for v in s1] == [round(float(v), 6) for v in s2]
      == [26.0, 30.0, 20.0],
      "the page draws a %s box and the run docks into %s"
      % (list(s1), [float(v) for v in s2]))
check(all(abs(a - float(b)) < 1e-6 for a, b in zip(c1, c2)),
      "the drawn and the docked box have different centers")

# --- the report leads with the mean of the top few ------------------------
row = report.summarise("AAAA", [{"score": -20.0}, {"score": -10.0},
                                {"score": -9.0}], top_n=2)
check(row["adcp_top2"] == -15.0,
      "the top-2 mean is %s, expected -15.0" % row["adcp_top2"])
check(row["adcp_best"] == -20.0, "the best score was not kept alongside it")
rows, key = report.rank([
    report.summarise("A", [{"score": -5.0}], mm=[-10.0], top_n=1),
    report.summarise("B", [{"score": -9.0}], mm=[-30.0], top_n=1)], top_n=1)
check(key.startswith("mmgbsa"),
      "the ranking used %s even though MM-GBSA was run" % key)
check(rows[0]["peptide"] == "B", "the better MM-GBSA did not rank first")

# --- the pose is laid on a clean topology, by key, not by file order -----
#
# ADCP writes the backbone through the chain and then comes back for the side
# chains, so a ranked pose lists a residue's atoms in two places. Copying
# those lines out in file order produced a PDB with duplicate residues, a
# topology twice the size it should be, and a PDB2PQR refusal: after the
# docking had already succeeded. The pose is therefore not copied: the
# peptide is built from the sequence and each atom is moved to the coordinate
# that matches it by (residue number, atom name), so the order in the file
# cannot matter.
SEQ = "AGGA"
tmpl = poseviewer.template(SEQ, capped=False)
check(tmpl is not None, "the peptide template for %s could not be built" % SEQ)
if tmpl is not None:
    atoms = [(a.GetMonomerInfo().GetResidueNumber(),
              a.GetMonomerInfo().GetName().strip(),
              a.GetMonomerInfo().GetResidueName().strip())
             for a in tmpl.GetAtoms()]
    # ADCP writes every backbone atom for every residue first, then every
    # side chain, so a residue's atoms are not contiguous in its output.
    back = [(r, n, rn) for r, n, rn in atoms if n in ("N", "CA", "C", "O")]
    side = [(r, n, rn) for r, n, rn in atoms if n not in ("N", "CA", "C", "O")]
    order = back + side
    import numpy as _np                                          # noqa: E402
    pose = dict(xyz=_np.array([(float(r), float(i), 0.0)
                               for i, (r, _n, _rn) in enumerate(order)]),
                resnum=_np.array([r for r, _n, _rn in order]),
                resname=[rn for _r, _n, rn in order],
                name=[n for _r, n, _rn in order])
    mol = poseviewer.transplant(SEQ, pose, capped=False)
    check(mol is not None,
          "an interleaved pose in ADCP's own order was refused")
    if mol is not None:
        check(mol.GetNumAtoms() == len(atoms),
              "the transplanted peptide has %d atoms, the template has %d: "
              "atoms were duplicated or dropped"
              % (mol.GetNumAtoms(), len(atoms)))
        conf = mol.GetConformer()
        want = {(r, n): (float(r), float(i), 0.0)
                for i, (r, n, _rn) in enumerate(order)}
        wrong = []
        for a in mol.GetAtoms():
            info = a.GetMonomerInfo()
            k = (info.GetResidueNumber(), info.GetName().strip())
            got = conf.GetAtomPosition(a.GetIdx())
            if k in want and abs(got.y - want[k][1]) > 1e-4:
                wrong.append(k)
        check(not wrong,
              "%d atom(s) took the coordinate of a different line, which is "
              "what copying in file order does: %s" % (len(wrong), wrong[:4]))
        seen = {}
        for a in mol.GetAtoms():
            info = a.GetMonomerInfo()
            seen.setdefault(info.GetResidueNumber(), set()).add(
                info.GetName().strip())
        check(len(seen) == len(SEQ),
              "the peptide came out with %d residues, not %d: the duplicate"
              "-residue bug" % (len(seen), len(SEQ)))

    # and a pose with an atom missing is refused, not scored
    short = dict(xyz=pose["xyz"][2:], resnum=pose["resnum"][2:],
                 resname=pose["resname"][2:], name=pose["name"][2:])
    check(poseviewer.transplant(SEQ, short, capped=False) is None,
          "a pose missing backbone atoms was accepted; the energy would be "
          "computed on a molecule with atoms left at template positions")

# --- the stages call the vendored pipeline, in its order -----------------
#
# A locally written receptor preparation or complex builder would run and
# produce files, and would only be found out when its energies turned out to
# be wrong. So the call sites are checked.
import inspect                                                   # noqa: E402
from adcp_pilot import job as jobmod, pipeline                    # noqa: E402

prep_src = inspect.getsource(pipeline.stage_prepare)
# The dock stage hands the rescoring to _stage_mmgbsa, which takes its poses
# from _mm_inputs, so the three are read as one, in the order they run.
_mm_src = inspect.getsource(pipeline._stage_mmgbsa)
check("_stage_mmgbsa(" in inspect.getsource(pipeline.stage_dock)
      and _mm_src.index("_mm_inputs(") < _mm_src.index("mmgbsa.run("),
      "stage_dock no longer rescores through _stage_mmgbsa, or it scores "
      "before it has taken its poses from _mm_inputs")
dock_src = (inspect.getsource(pipeline.stage_dock)
            + inspect.getsource(pipeline._mm_inputs) + _mm_src)

for call, src, why in [
        ("receptor.load_for_box(", prep_src,
         "the chains, the validated hotspots and the box"),
        ("prep.run(", prep_src,
         "the full Amber preparation: altloc, hetatms, disulfides, missing "
         "atoms, protonation, restrained minimization"),
        ("receptor.prepare_pdbqt(", prep_src, "the receptor as PDBQT"),
        ("target.build(", prep_src, "the grid maps, built and then checked"),
        ("adcp11.dock(", dock_src, "the docking search"),
        ("adcp.read_pose(", dock_src, "the all-atom ranked pose"),
        ("poseviewer.transplant(", dock_src,
         "the pose laid on a clean topology by key"),
        ("mmgbsa.run(", dock_src,
         "the rescoring: the batch entry, which caches per pose and "
         "runs them several at a time")]:
    check(call in src,
          "%s no longer calls %s: %s is the engine's and is not to be "
          "written again here" % ("stage_prepare" if src is prep_src
                                  else "stage_dock", call.rstrip("("), why))

check(prep_src.index("receptor.load_for_box(") < prep_src.index("prep.run(")
      < prep_src.index("receptor.prepare_pdbqt(")
      < prep_src.index("target.build("),
      "stage_prepare runs the preparation steps out of order; the box has to "
      "be sized from the structure as read, the PDBQT has to come from the "
      "prepared structure, and the maps from the PDBQT")
check(dock_src.index("adcp11.dock(") < dock_src.index("poseviewer.transplant(")
      < dock_src.index("mmgbsa.run("),
      "stage_dock rescores before it docks, or builds the complex from "
      "something other than the pose")
check("minimise_structure=bool(j.get(\"minimise\"" in prep_src,
      "stage_prepare does not pass the minimization through; protonation "
      "places a full set of hydrogens and nothing would relax them")
check(jobmod.default()["minimise"] is True,
      "a run created with no flags does not minimize")
check("verify.mmgbsa_split(" in dock_src,
      "an MM-GBSA value is kept without reading back whether the complex "
      "was split at the interface")

# nothing in the package may shadow a vendored module with a local rewrite
import os as _os                                                 # noqa: E402
_pkg = _os.path.join(_os.path.dirname(HERE), "adcp_pilot")
_local = {f[:-3] for f in _os.listdir(_pkg) if f.endswith(".py")}
_vend = {f[:-3] for f in _os.listdir(_os.path.join(_pkg, "_vendor"))
         if f.endswith(".py") and f != "__init__.py"}
check(not (_local & _vend),
      "these modules exist both here and in _vendor: %s. Two copies of the "
      "same stage is how the three failed runs happened."
      % sorted(_local & _vend))

# --- the affinity ruler runs one way only ---------------------------------
#
# Converting an MM-GBSA score into a Kd is not meaningful: an MM-GBSA
# interaction energy has no configurational entropy term, so it is not the
# quantity in dG = RT ln(Kd), and its magnitude follows the interface and the
# receptor more than the affinity. The absence of that conversion is asserted
# so that it is not added later.
from adcp_pilot import affinity                                   # noqa: E402
from adcp_pilot._vendor import mmgbsa as VMM                      # noqa: E402

check(hasattr(VMM, "dg_for_kd"),
      "the ruler is gone; there is now no way to answer 'what score means "
      "nanomolar'")
for gone in ("kd_for_dg", "predicted_kd", "kd_of", "affinity_of"):
    check(not hasattr(VMM, gone),
          "mmgbsa.%s exists. The inverse of this ruler cannot be computed "
          "from an MM-GBSA score: see the comment above this check" % gone)
check(abs(VMM.dg_for_kd(1e-9) - (-12.28)) < 0.01,
      "1 nM is not -12.28 kcal/mol; the constant or the temperature moved")
_window = VMM.dg_for_kd(1e-9) - VMM.dg_for_kd(1e-12)
check(abs(_window - 4.1) < 0.05,
      "the picomolar-to-nanomolar window is %.2f, not 4.1" % _window)

_lines = []
affinity.table(out=_lines.append)
_text = "\n".join(_lines)
check("4.1 kcal/mol" in _text, "the ruler does not print the window")
check("never a Kd" in _text or "not a Kd" in _text,
      "the ruler prints a table of energies without saying why no column "
      "reports a Kd, which is the only reason it exists")

# the spread reaches the table, and the caveat reaches the reader
_rows = [report.summarise("AAAA", [{"score": -20.0}], mm=[-30.0, -12.0],
                          top_n=2)]
check(_rows[0]["mmgbsa_spread"] == 18.0,
      "mmgbsa_spread is %r, expected 18.0" % _rows[0].get("mmgbsa_spread"))
check(report.summarise("A", [{"score": -1.0}], mm=[-5.0], top_n=2)
      ["mmgbsa_spread"] is None,
      "a single pose was given a spread; there is nothing to spread")
_tbl = report.table(_rows, top_n=2)
check("spread" in _tbl, "the spread column is not in the printed table")
# The results say what their columns are and nothing about what they are not.
# A note about Kd under a table with no Kd in it answered a discussion, not the
# table; the affinity command, which is about Kd, is where that belongs.
check("Kd" not in _tbl and "entropy" not in _tbl,
      "the printed table talks about quantities it does not contain")

# --- two peptides with one sequence would overwrite each other ------------
#
# A docking is written to docking/<sequence>/ and its result to
# results/<sequence>.json. Two entries carrying the same sequence therefore
# share a working directory AND a result file: the second overwrites the
# first, `outstanding` counts both as finished once one has run, and the
# table shows two rows from one calculation. That is not untidy, it is
# wrong, and nothing downstream would say so.
from adcp_pilot import job as _job                               # noqa: E402

_base = dict(_job.default(), receptor=__file__, hotspots="A:1", poses=5,
             mmgbsa=False)


def _problems(peps):
    return _job.validate(dict(_base, peptides=peps))


_dup = _problems([{"name": "a", "sequence": "AKSLW"},
                  {"name": "b", "sequence": "akslw"}])
check(any("same sequence" in p for p in _dup),
      "two entries with the same sequence were accepted; one would "
      "overwrite the other's results: %s" % _dup)
check(any("AKSLW" in p for p in _dup),
      "the refusal does not say WHICH sequence: %s" % _dup)

_dupname = _problems([{"name": "Test", "sequence": "AKSLW"},
                      {"name": "test", "sequence": "PLLSLL"}])
check(any("named" in p for p in _dupname),
      "two peptides with the same name were accepted. The name labels the "
      "row in the results table: %s" % _dupname)

_clean = _problems([{"name": "a", "sequence": "AKSLW"},
                    {"name": "b", "sequence": "PLLSLL"}])
check(not [p for p in _clean if "same sequence" in p or "named" in p],
      "two different peptides were refused: %s" % _clean)

# a blank name is filled in by position, and stays put between runs: the
# name reaches the results table and a log line, and one that moved would
# make two reports of the same job disagree
_named = _job.peptide_list(dict(_base, peptides=[
    {"name": "", "sequence": "AKSLW"},
    {"name": "   ", "sequence": "PLLSLL"},
    {"name": "mine", "sequence": "HYGSGSTW"}]))
check([n for n, _q in _named] == ["peptide_1", "peptide_2", "mine"],
      "blank names were not filled in by position: %s"
      % [n for n, _q in _named])
check(_job.peptide_list(dict(_base, peptides=[{"name": "", "sequence": "AK"}]))
      == _job.peptide_list(dict(_base, peptides=[{"name": "", "sequence": "AK"}])),
      "the same job file produces different names on two reads")

# blank names are NOT a refusal: they have an obvious right answer
check(not [p for p in _problems([{"name": "", "sequence": "AKSLW"}])
           if "named" in p],
      "a blank name was treated as a duplicate")

# --- the queue says which run each job belongs to -------------------------
# Three rows called adcp_prep, adcp_pilot and adcp_report tell a user with
# four runs queued nothing at all: not which receptor, not which peptides,
# not which of their own runs is the one that is stuck.
from adcp_pilot import pipeline as _pl                            # noqa: E402

check(_pl.job_name("/runs/DENV1_7A3R_screen", {}, "dock")
      == "DENV1_7A3R_screen-dock",
      "the job name does not carry the run directory's name")
check(_pl.job_name("/runs/r1", {"label": "PD-L1 alpha helix v2"}, "prep")
      == "PD-L1-alpha-helix-prep",
      "a typed label does not reach the job name, or is not made safe for "
      "SLURM: %r" % _pl.job_name("/runs/r1", {"label": "PD-L1 alpha helix v2"},
                                 "prep"))
check(not _pl.job_name("/runs/a-very-long-run-name-indeed", {},
                       "report").startswith("a-very-long-run-name-indeed"),
      "a long run name is not shortened, so squeue's narrow column shows "
      "nothing useful")
check("--" not in _pl.job_name("/runs/r1", {"label": "PD-L1 alpha helix v2"},
                               "dock"),
      "the length cut left a separator against the stage suffix")
# "adcp" is job_name's fallback and "dock" the stage: the two together are
# not the package's name, and the rename to adcp-pilot does not touch them.
check(_pl.job_name("/runs/   ", {"label": "   "}, "dock") == "adcp-dock",
      "a blank run name produces a nameless job rather than a fallback")
_src = inspect.getsource(_pl)
check('name="adcp_prep"' not in _src and 'name="adcp_pilot"' not in _src,
      "a stage still submits under the fixed name")

# --- the preparation estimate is not six times the measurement ------------
# The floor was sized for a TYPICAL receptor and then the per-atom and
# per-box terms were added on top of it, so a typical receptor was counted
# twice: 27-30 minutes shown against 4 min 36 s taken. An estimate no one
# believes is not a safety margin. What protects a run from the wall is the
# two-hour floor, not a padded stage.
from adcp_pilot import sizing as _sz                               # noqa: E402

# Measured with the ADCP 1.1 suite, the whole prepare stage, protonated by
# reduce, maps with AutoSite 1.1 and no gradient. The four-minute floor and
# reduce's minute are sized for a slow node, and they are most of a small
# receptor's estimate, so the bound is wider there.
for _atoms, _kinds, _box, _meas, _most, _what in (
        (705, 0, (30, 30, 30), 1.25, 6.0, "MDM2, 705 atoms, 30 A box"),
        (5912, 1, (26.6, 25.5, 29.1), 4.3, 3.5,
         "the 7A3R dimer, 5,912 atoms, one hetero group")):
    _e = _sz.prep_minutes(atoms=_atoms, hetatm_kinds=_kinds, box=_box,
                          gradient=False, autosite="1.1", protonation="reduce")
    check(_e >= _meas,
          "the prep estimate for %s is %.1f min, UNDER the %.1f min it took: "
          "the job would be killed at the wall" % (_what, _e, _meas))
    check(_e <= _meas * _most,
          "the prep estimate for %s is %.1f min against %.1f measured (%.1fx). "
          "Over-booking is free, but an estimate this far out stops being "
          "read" % (_what, _e, _meas, _e / _meas))

# the floor is a floor for the smallest receptor, not for a typical one
check(_sz.prep_minutes(atoms=0, box=None) == _sz.PREP_FLOOR_MINUTES,
      "the floor does not apply to a receptor of unknown size")
check(_sz.prep_minutes(atoms=5884, box=(20, 20, 20))
      > _sz.prep_minutes(atoms=1793, box=(20, 20, 20)),
      "a bigger receptor is not estimated to take longer, so the per-atom "
      "term is being swallowed by the floor")

# and the thing that actually protects the run is still there
check(_sz.FLOOR_MINUTES >= 120,
      "the two-hour walltime floor is gone; that, not a padded stage "
      "estimate, is what stops a job being killed at the wall")


# --- the docking model is waves, and it matches every measurement ---------
#
# It was core-minutes, K x replicas x Msteps / cores, fitted to ONE run,
# and it cannot be right: halving the replicas below the core count does not
# halve the wall time, because those replicas already ran side by side. It
# told a user two minutes for a job that ran twenty.
#
# Every figure below is ADCP 1.1.21 docking into MDM2 on a compute node,
# timed around adcp11.dock: the whole docking, from starting ADCP to the
# poses split into their files (sizing.py has the table).
TIMED = (
    (8, 250_000, 8, 6, 41.1, "8 x 250k on 8 cores, 6-mer"),
    (8, 1_000_000, 8, 6, 72.6, "8 x 1M on 8 cores, 6-mer"),
    (8, 2_500_000, 8, 6, 132.1, "8 x 2.5M on 8 cores, 6-mer"),
    (8, 250_000, 8, 9, 79.2, "8 x 250k on 8 cores, 9-mer"),
    (8, 1_000_000, 8, 9, 96.1, "8 x 1M on 8 cores, 9-mer"),
    (8, 2_500_000, 8, 9, 200.6, "8 x 2.5M on 8 cores, 9-mer"),
    (16, 1_000_000, 8, 9, 202.0, "16 x 1M on 8 cores, 9-mer, two waves"),
    (8, 250_000, 8, 13, 65.8, "8 x 250k on 8 cores, 13-mer"),
    (8, 1_000_000, 8, 13, 163.5, "8 x 1M on 8 cores, 13-mer"),
    (8, 2_500_000, 8, 13, 307.7, "8 x 2.5M on 8 cores, 13-mer"))
for _reps, _steps, _cores, _len, _meas, _what in TIMED:
    _m = _sz.dock_minutes(_len, _reps, _steps, _cores) * 60.0
    check(_m >= _meas,
          "%s is estimated at %.0f s, UNDER the %.0f s it took" % (_what, _m, _meas))
    check(_m <= _meas * 3.0,
          "%s is estimated at %.0f s against %.0f measured (%.1fx)"
          % (_what, _m, _meas, _m / _meas))

# below the core count, replicas are free: that is the intent
check(_sz.dock_minutes(6, 4, 1_000_000, 16)
      == _sz.dock_minutes(6, 16, 1_000_000, 16),
      "4 replicas and 16 replicas on 16 cores are estimated differently, so "
      "the model is still dividing rather than counting waves")
check(_sz.dock_minutes(6, 17, 1_000_000, 16)
      > _sz.dock_minutes(6, 16, 1_000_000, 16),
      "the 17th replica starts a second wave and must cost one")
check(_sz.dock_minutes(6, 32, 1_000_000, 16)
      > _sz.dock_minutes(6, 32, 1_000_000, 32),
      "more cores must finish the same replicas sooner")
check(not hasattr(_sz, "K"),
      "the core-minute constant is still here; two models will drift apart")

# and the estimate names the jobs the queue will show, not internal stages
_p = _sz.suggest(["AKLLWH"], replicas=10, steps=2_500_000, cores=48,
                 mmgbsa=True, mm_poses=5, atoms=2837, box=(30, 30, 30))
check("dock+mm" in _p["why"],
      "the estimate still quotes docking separately from the rescoring, "
      "while SLURM runs them as one job, which is what made a 2-minute "
      "docking look like a 20-minute one")
check("one job" in _p["why"] or "as one" in _p["why"],
      "nothing says the two run as a single task")
_q = _sz.suggest(["AKLLWH"], replicas=10, steps=2_500_000, cores=48,
                 mmgbsa=False, atoms=2837, box=(30, 30, 30))
check("dock+mm" not in _q["why"],
      "with rescoring off the task is still described as dock+mm")


# --- the estimate scales with peptide length, at every length -------------
#
# It counted only residues past TEN, so a 4-mer and a 10-mer were costed the
# same, and ADCP's own guidance is that the problem scales at every length:
# "The longer the peptide is, the more difficult is the docking problem"
# (ccsb.scripps.edu/adcp/documentation).
check(_sz.REFERENCE_LENGTH <= 6,
      "the length factor is anchored at %d residues, so everything shorter "
      "is costed identically" % _sz.REFERENCE_LENGTH)
_a = _sz.dock_minutes(6, 50, 5_000_000, 48)
_b = _sz.dock_minutes(10, 50, 5_000_000, 48)
_c = _sz.dock_minutes(20, 50, 5_000_000, 48)
check(_a < _b < _c,
      "a longer peptide is not estimated to take longer: 6->%.1f 10->%.1f "
      "20->%.1f min" % (_a, _b, _c))

# and a longer peptide, at the same settings, is never estimated under what
# it took: the per-residue factor carries the measurements above.
_by_len = {}
for _reps, _steps, _cores, _len, _meas, _what in TIMED:
    if (_reps, _steps) == (8, 2_500_000):
        _by_len[_len] = _meas
check(sorted(_by_len) == [6, 9, 13] and _by_len[6] < _by_len[9] < _by_len[13],
      "the measurements no longer cover three lengths: %s" % _by_len)


# --- ADCP 1.1: sequence syntax, output layout, OpenMM table, target files ---
import shutil as _sh                                            # noqa: E402
import tempfile as _tf                                          # noqa: E402
import zipfile as _zf                                           # noqa: E402
from adcp_pilot import adcp11, agfr_target, job as _job, peptide as _pep  # noqa: E402
from adcp_pilot import report as _rep                            # noqa: E402

_s, _p = _pep.check("n<BHD>pis&dvd", library={"BHD"})
check(not _p and len(_pep.residues(_s)) == 7,
      "an ADCP 1.1 sequence with a non-standard and a D residue was refused "
      "or miscounted: %r %r" % (_s, [x.message for x in _p]))
_s, _p = _pep.check("n<BHD>pis&dvd")
check(not _p, "ADCP 1.1's syntax needs naming before it is accepted: %r"
      % [x.message for x in _p])
_s, _p = _pep.check("n<ZZZ>pisdv", library={"BHD"})
check(any(x.kind == "unknown_nst" for x in _p),
      "a side chain missing from ADCP's library was accepted")
check(_job.key("NPISDVD") == "NPISDVD" and _job.key("N<BHD>PIS&DVD") == "N_BHD_PIS_DVD",
      "folder names changed for ordinary peptides, or kept <, > and &")
check(adcp11.adcp_sequence("N<BHD>PIS&DVD") == "n<BHD>pis&dvd",
      "the sequence passed to ADCP 1.1 lower-cased a library name")

_wd = _tf.mkdtemp()
for _f in ("out_out.pdb", "out_summary.dlg"):
    _sh.copy(os.path.join(HERE, "data", "adcp11", _f), _wd)
check(adcp11.split_poses(_wd) == 5
      and all(os.path.isfile(os.path.join(_wd, "out_ranked_%d.pdb" % k))
              for k in range(1, 6)),
      "ADCP 1.1's multi-model file was not split into one pose per file")
_m = adcp11.read_modes(_wd)
check([x["rank"] for x in _m] == [1, 2, 3, 4, 5] and _m[0]["affinity"] == -10.5,
      "ADCP 1.1's mode table was not read: %r" % _m[:2])
_o = adcp11.read_openmm(_wd)
check(len(_o) == 3 and _o[0]["rank_adcp"] == 3
      and _o[0]["e_complex_receptor"] == -262.4,
      "the OpenMM re-ranking table was not read: %r" % _o[:1])
_row = _rep.summarise("N<BHD>PIS&DVD", [{"score": -10.5}, {"score": -10.0}],
                      openmm=_o, engine="ADCP 1.1.21")
check(_row["length"] == 7 and _row["omm_best"] == -262.4
      and _row["omm_best_pose"] == 3 and _row["engine"] == "ADCP 1.1.21",
      "the report row miscounts residues or drops the OpenMM result: %r" % _row)

# A small target file with the layout agfr writes.
_t = _tf.mkdtemp()
_trg = os.path.join(_t, "mytarget.trg")
with _zf.ZipFile(_trg, "w") as _z:
    for _k in range(20):
        _z.writestr("mytarget/rigidReceptor.T%d.map" % _k, "0\n")
    _z.writestr("mytarget/rigidReceptor.pdbqt",
                "ATOM      1  N   ALA A   1      1.000   2.000   3.000  1.00  0.00    -0.300 N \n"
                "ATOM      2  CA  ALA A   1      2.000   2.000   3.000  1.00  0.00     0.100 C \n"
                "ATOM      3  OG  SER A   2      3.000   2.000   3.000  1.00  0.00    -0.400 OA\n")
    _z.writestr("mytarget/rigidReceptor.gpf", "npts 80 60 40 # grid\nspacing 0.375\n"
                "gridcenter 1.0 2.0 3.0\n")
_run = _tf.mkdtemp()
_trg2, _prep, _q, _c, _sz = agfr_target.import_target(_trg, _run, log=lambda *a: None)
check(os.path.basename(_trg2) == "mytarget.trg" and _c == [1.0, 2.0, 3.0]
      and _sz == [30.0, 22.5, 15.0],
      "a target file lost its name or its box: %r %r %r" % (_trg2, _c, _sz))
_pdb = open(_prep).read()
check(_pdb.count("ATOM") == 3 and " O\n" in _pdb.replace("  O\n", " O\n"),
      "the receptor from a target file was not written as a PDB with elements")
_j = _job.default(target_file=_trg, peptides=[{"name": "p", "sequence": "NPISDVD"}])
check(not [x for x in _job.validate(_j) if "receptor" in x.lower() or "site" in x.lower()],
      "a run from a target file still asks for a receptor or a site: %r"
      % _job.validate(_j))
if adcp11.find():
    _j = _job.default(receptor=_prep, hotspots="A:1",
                      peptides=[{"name": "p", "sequence": "N<ZZZ>PIS&DVD"}])
    check(any("ZZZ" in x for x in _job.validate(_j)),
          "a run accepted a side chain the installed ADCP 1.1 does not have")

# --- the Modifications dialog's list, and OpenMM parameter sets ------------
_mods = json.load(open(os.path.join(os.path.dirname(HERE), "adcp_pilot", "static",
                                    "modifications.json")))["modifications"]
check(len(_mods) > 100 and all(m["parent"] in _pep.STANDARD and m["coarse"]
                               in _pep.STANDARD and m["code"] and m["name"]
                               for m in _mods),
      "the modification list is empty or has an entry without a parent, a "
      "coarse residue, a code or a name")
_exe = adcp11.find()
if _exe:
    for _seq, _want in (("R&WGPGAF<KYN>GS", "swiss"), ("RWGPGAF<0A1>GS", "sannerlab"),
                        ("F<KYN>WGPGAS", None), ("RWGPGAWGS", "")):
        check(adcp11.openmm_set_for(_exe, _seq) == _want,
              "the OpenMM parameter set for %s is %r, expected %r"
              % (_seq, adcp11.openmm_set_for(_exe, _seq), _want))
    _lib = adcp11.library(_exe)
    check(all(m["code"] in _lib for m in _mods),
          "the modification list offers side chains ADCP 1.1 does not have")
else:
    print("skip  ADCP 1.1 is not installed here; its parameter checks were skipped")

# --- nothing in adcp_pilot/ shadows a module of the engine -------------------
# adcp_pilot/_vendor holds the engine. A module of the same name directly in
# adcp_pilot/ is easily imported in its place, and replaces a stage silently.
_vend = {f for f in os.listdir(os.path.join(os.path.dirname(HERE), "adcp_pilot", "_vendor"))
         if f.endswith(".py") and f != "__init__.py"}
_top = {f for f in os.listdir(os.path.join(os.path.dirname(HERE), "adcp_pilot"))
        if f.endswith(".py")}
check(not (_vend & _top), "these modules exist both in adcp_pilot/ and in _vendor/: %s"
      % sorted(_vend & _top))

# --- a failed OpenMM minimization keeps the docking --------------------------
# ADCP writes its poses, then exits 1 when OpenMM cannot use a residue's
# template (ADCP 1.1.21 does for FVAL, FPG2, FPG3, FTR and FT6). The poses
# are kept, the peptide is marked done so a resumed run does not dock it
# again, and the log says its poses are not minimized.
import subprocess as _sp                                            # noqa: E402
_wd2 = _tf.mkdtemp()
_trg = os.path.join(_wd2, "receptor.trg")
open(_trg, "wb").write(b"x")


def _fake_run(cmd, wd, exe, append=False):
    _sh.copy(os.path.join(HERE, "data", "adcp11", "out_out.pdb"), wd)
    with open(os.path.join(wd, "out_summary.dlg"), "w") as fh:
        fh.writelines(x for x in open(os.path.join(HERE, "data", "adcp11", "out_summary.dlg"))
                      if not x.startswith("OMM"))
    return _sp.CompletedProcess(cmd, 1)


_saved = (adcp11._run, adcp11.command, adcp11.openmm_set_for)
adcp11._run = _fake_run
adcp11.command = lambda exe: ["adcp"]
adcp11.openmm_set_for = lambda exe, seq: "sannerlab"
_said = []
try:
    _modes, _wdd = adcp11.dock("GAA<FVAL>AG", _trg, _wd2, "/nowhere/bin/adcp", key="fval",
                               minimize=1, log=_said.append)
finally:
    adcp11._run, adcp11.command, adcp11.openmm_set_for = _saved
check(_modes and os.path.isfile(os.path.join(_wdd, ".adcp_finished")),
      "a docking whose OpenMM minimization failed was not kept as done")
check(any("minimization failed" in x for x in _said),
      "a failed OpenMM minimization was not reported: %r" % _said)
_sh.rmtree(_wd2, ignore_errors=True)

# --- MM-GBSA uses the radii GB-OBC(II) was fitted with ------------------------
# igb=5 goes with mbondi2; tleap's default is mbondi. A run scored before the
# radii were set keeps mbondi, so it is never part one set and part another.
import tempfile as _tf                                             # noqa: E402
_leap = mmgbsa.LEAP.format(complex_pdb="c.pdb", receptor_pdb="r.pdb",
                           ligand_pdb="l.pdb", prefix="X", radii=mmgbsa.RADII)
check(mmgbsa.RADII == "mbondi2" and "set default PBRadii mbondi2" in _leap
      and _leap.index("PBRadii") < _leap.index("loadpdb"),
      "tleap does not build the complex with mbondi2 radii: %r" % _leap)
with _tf.TemporaryDirectory() as _t:
    _new = os.path.join(_t, "run1", "mmgbsa")
    check(mmgbsa.radii_for(_new) == "mbondi2"
          and open(os.path.join(_new, "RADII")).read().strip() == "mbondi2",
          "a new run's MM-GBSA folder is not given and does not record mbondi2")
    _old = os.path.join(_t, "run2", "mmgbsa")
    os.makedirs(_old)
    with open(os.path.join(_old, "ETFSDLWKLLPEN.dg"), "w") as _fh:
        _fh.write("# igb 5 saltcon 0.150 intdiel 1\n1 -40.0\n")
    check(mmgbsa.radii_for(_old) == "mbondi",
          "a run scored before the radii were set does not keep mbondi")
    _dg = os.path.join(_t, "X.dg")
    mmgbsa._write_cache(_dg, {1: -40.0}, 1.0, "mbondi2")
    check(mmgbsa._read_cache(_dg, 1.0, "mbondi2") == {1: -40.0}
          and mmgbsa._read_cache(_dg, 1.0, "mbondi") == {},
          "cached energies are reused under different radii")
check("kept.append(ranks[r - 1]" in dock_src,
      "MM-GBSA energies are recorded by their position among the scored "
      "poses, not by ADCP rank, so a skipped pose shifts every later one")
check('out["mmgbsa_radii"]' in dock_src,
      "a peptide's result does not record the radii its MM-GBSA used")

# --------------------------------------------------------------------------
# Pose viewer files. A finished run produced tables and nothing anybody could
# open: the poses were on disk only as ADCP's own out_ranked_*.pdb inside
# each peptide's docking directory, with no receptor beside them.
print("the pose viewer stage")
from adcp_pilot import pipeline as _pl                          # noqa: E402
from adcp_pilot import bundle as _bundle                        # noqa: E402
import inspect as _inspect                                     # noqa: E402

check(callable(getattr(_pl, "stage_poseviewer", None)),
      "there is a pose viewer stage")
_pvsrc = _inspect.getsource(_pl.stage_poseviewer)
# Poses are rebuilt, never copied: ADCP writes the peptide alone, ordered for
# docking rather than for reading.
check("transplant(" in _pvsrc,
      "poses are transplanted onto a freshly built peptide, not copied")
check("with_hydrogens(" in _pvsrc, "and carry their hydrogens")
# The receptor written beside them must be the one the poses were docked
# into, or the complex shows a pose against the wrong protonation.
check("prepared_pdb" in _pvsrc and "write_view_pdb(" in _pvsrc,
      "the receptor written beside them is the prepared structure")
check('"L"' in _pvsrc and '"A"' in _pvsrc,
      "receptor and peptide land on separate chains")
for _w in ("write_pdb(", "write_sdf(", "write_mae(", "write_maegz("):
    check(_w in _pvsrc, "the stage writes %s" % _w.rstrip("("))
# The sequence identifies the file. Names default to "peptide 1".
check("stem" in _pvsrc, "file names carry the sequence, not just the name")
# A run that produced a table must not be failed by a viewer file.
_rpt = _inspect.getsource(_pl.stage_report)
check("stage_poseviewer(" in _rpt, "the report stage calls it")
check("could not be written" in _rpt,
      "and a failure there does not fail the run")
# And they have to reach the person: a results download that omits them
# leaves the run as tables and logs again.
for _p in ("poseviewer/*.pdb", "poseviewer/*.mae", "poseviewer/*.maegz",
           "poseviewer/*.sdf"):
    check(_p in _bundle.RESULTS, "a results download carries %s" % _p)


# --------------------------------------------------------------------------
# More than one MM-GBSA engine on a run, and rescoring a finished one.
print("listing the engines does not import OpenMM")
# The page's first request lists every engine. Importing OpenMM to learn
# whether it is installed reads hundreds of files, and on a busy shared
# filesystem that outlasted the API test's 30 second limit. Finding the
# package answers the question without loading it.
import subprocess as _sp3                                       # noqa: E402
_code = ("import sys; sys.path.insert(0, %r); "
         "from adcp_pilot._vendor import gbengines as G; G.catalogue(); "
         "print('openmm' in sys.modules)" % os.path.dirname(HERE))
_out = _sp3.run([sys.executable, "-c", _code], capture_output=True,
                text=True).stdout.split()
check(_out == ["False"],
      "listing the engines imported OpenMM (%s)" % _out)

print("several MM-GBSA engines on one run")
from adcp_pilot import report as _rep                            # noqa: E402
_r = _rep.summarise("LLKWHY", [{"score": -14.0}], mm=[-25.0], name="P1",
                    engine="ADCP 1.1.21", mmgbsa_engine="amber",
                    mm_by_engine={"amber": [-25.0], "openmm": [-26.3]})
# The docking engine and the MM-GBSA engine are different questions. Reading
# one for the other labeled every energy in the .mae files "ADCP 1.1.21".
check(_r["engine"] == "ADCP 1.1.21" and _r["mmgbsa_engine"] == "amber",
      "a row keeps the docking engine and the MM-GBSA engine apart")
check(_r.get("mmgbsa_openmm_best") == -26.3 and _r.get("mmgbsa_amber_best") == -25.0,
      "each engine that scored a peptide has its own columns")
_old = _rep.summarise("LLKWHY", [{"score": -14.0}], mm=[-25.0])
check(_old["mmgbsa_engine"] == "amber",
      "a run from before engines were selectable reads as AmberTools")
_pv = _inspect.getsource(_pl.stage_poseviewer)
check('r.get("mmgbsa_engine")' in _pv and 'props["MMGBSA_engine"] = str(r["engine"])' not in _pv,
      "the pose viewer labels the MM-GBSA engine from mmgbsa_engine, not engine")

# Rescoring leaves the numbers a run is RANKED on exactly as they were: the
# engines are on different scales, so re-ranking on another would change a
# finished table's order on a number that was never its basis.
_rs = _inspect.getsource(_pl.stage_rescore)
check('open(f, "w")' not in _rs and "json.dump(" not in _rs,
      "a rescoring does not rewrite the run's result files")
check('engine=eng.name' in _rs, "a rescoring runs the engine it was asked for")
check('1 if eng.name == "schrodinger"' in _rs,
      "Prime rescores one pose at a time, so one column is one force field")
check(callable(getattr(_pl, "submit_rescore", None)),
      "a rescoring can be sent to the queue")
from adcp_pilot import worker as _wk                              # noqa: E402
check('startswith("rescore:")' in _inspect.getsource(_wk.main),
      "the worker runs a rescore stage")

# The page's shared engine values must be at the top level of app.js. They
# were first declared inside the tools callback, where the results view could
# not see them, and every run threw on opening. A syntax check passes that.
_js = open(os.path.join(os.path.dirname(HERE), "adcp_pilot", "static",
                        "app.js"), encoding="utf-8").read()
for _name in ("const GBLABEL",):
    _line = [l for l in _js.splitlines() if l.startswith(_name)]
    check(bool(_line), "%s is declared at the top level of app.js" % _name)
# Peptide residues are matched by position: the residue map numbers the
# peptide on from the receptor (B252 ...), the viewer keeps ADCP's numbers.
check("inView[at]" in _js and "focusResidue" in _js,
      "a peptide residue is found in the viewer by its position, not its label")


# Prime's headline dG_Bind subtracts the strain of a separately relaxed free
# peptide, which grew from 9 to 59 kcal/mol between a 5 and a 15 residue
# peptide and turned -62.2 into -2.9. The column must read the no-strain
# value, the quantity the single-trajectory engines report, and never fall
# back to the headline one when it is absent.
import csv as _csv, tempfile as _tf                              # noqa: E402
_d = _tf.mkdtemp()
_p = os.path.join(_d, "x-out.csv")
with open(_p, "w", newline="") as _fh:
    _w = _csv.writer(_fh)
    _w.writerow(["title", "r_psp_MMGBSA_dG_Bind", "r_psp_MMGBSA_dG_Bind(NS)",
                 "r_psp_Lig_Strain_Energy"])
    _w.writerow(["x", "-2.89", "-62.22", "59.32"])
_got = VMM._read_prime_csv(_p)
check(_got == -62.22, "Prime is read as its no-strain binding energy")
with open(_p, "w", newline="") as _fh:
    _w = _csv.writer(_fh)
    _w.writerow(["title", "r_psp_MMGBSA_dG_Bind"])
    _w.writerow(["x", "-2.89"])
_got2 = VMM._read_prime_csv(_p)
check(_got2 is None, "and gives nothing rather than the headline value when that is absent")


if fails:
    print("test_engine FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  targets are checked by content, the complex is split by chain, "
      "the stages call the vendored pipeline, and the ruler runs one "
      "way, %d checks" % ran)

