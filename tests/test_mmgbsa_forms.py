#!/usr/bin/env python3
"""MM-GBSA for D residues and for closed rings, and what an old job still means.

By Aadhil Haq

What is checked, without Amber: the folder name a D peptide's MM-GBSA gets,
and that it can never be an L peptide's; that a D residue reaches tleap under
its L name with its own coordinates; the closure test and its two cutoffs,
for a head to tail ring and for a disulfide with two and with three
cysteines; the tleap input for each kind of ring; that a ring with no closed
pose is refused with the reason and never built; that Prime refuses both;
and that a job file from before this change keeps skipping D residues and
rings, whether its peptides are docked now or rescored later.

With Amber: that ff14SB cannot tell a D peptide from its mirror image (the
torsion phases in a built topology, and the energy of a D peptide against
its exact mirror image as an L peptide, before and after a minimization);
that tleap alone puts the hydrogens of a D residue in the wrong place, and
mirror_hydrogens puts them right; and that both kinds of ring build, with
the closing bond in the topology, no termini on a head to tail ring and an
integer charge.

    python3 tests/test_mmgbsa_forms.py

Exit 0 passes, 77 skips (no Amber here, after everything else passed), and
anything else fails.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import numpy as np                                            # noqa: E402

from adcp_pilot import job as jobmod, mmforms, pipeline, report  # noqa: E402
from adcp_pilot._vendor import amberfiles as af                # noqa: E402
from adcp_pilot._vendor import library, mmgbsa, poseviewer, prep  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)
        print("  FAIL " + why)


TMP = tempfile.mkdtemp(prefix="adcp_mmforms_")
os.environ["HOME"] = TMP


def pose_of(seq, mirror=False, shift=(0.0, 0.0, 0.0)):
    """A pose dict, as adcp.read_pose gives one, from RDKit's build of the L
    sequence: heavy atoms, ADCP's numbering from 1, and reflected through
    the yz plane when asked, which makes every residue its D form."""
    m = library.build_3d(seq, capped=False)
    from rdkit import Chem
    m = Chem.RemoveHs(m)
    conf = m.GetConformer()
    xyz, resnum, resname, name = [], [], [], []
    for at in m.GetAtoms():
        p = conf.GetAtomPosition(at.GetIdx())
        info = at.GetMonomerInfo()
        x = np.array([p.x, p.y, p.z]) + np.array(shift)
        if mirror:
            x[0] = -x[0]
        xyz.append(x)
        resnum.append(info.GetResidueNumber())
        resname.append(info.GetResidueName().strip())
        name.append(info.GetName().strip())
    return {"xyz": np.array(xyz), "resnum": np.array(resnum),
            "resname": resname, "name": name}


def ca_sign(pose, k):
    """+1 for an L alpha carbon, -1 for a D one, in a pose dict."""
    a = {n: pose["xyz"][i] for i, n in enumerate(pose["name"])
         if int(pose["resnum"][i]) == k}
    return np.sign(af.signed_volume(a["N"], a["CA"], a["C"], a["CB"]))


def pdb_lines(pose, chain, first=1, drop=()):
    """ATOM lines for a pose dict, residues numbered from `first`."""
    out, n = [], 0
    for i, nm in enumerate(pose["name"]):
        if nm in drop:
            continue
        n += 1
        x, y, z = pose["xyz"][i]
        r = int(pose["resnum"][i]) - 1 + first
        el = nm[0]
        out.append("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                   "          %2s\n" % (n, nm if len(nm) == 4 else " " + nm,
                                        pose["resname"][i], chain, r, x, y,
                                        z, el))
    return out


# --- names and the job flag ---------------------------------------------------
check(mmforms.name("ETFSDLWKLLPEN") == "ETFSDLWKLLPEN",
      "a standard peptide's MM-GBSA folder keeps the name every run used")
check(mmforms.name("&T&NWY&A") == "_T_NWY_A",
      "a D residue's & becomes _ in the folder name (%s)" % mmforms.name("&T&NWY&A"))
check(mmforms.name("&TNWYA").lower() != mmforms.name("TNWYA").lower()
      and mmforms.name("D&T").lower() != mmforms.name("&DT").lower(),
      "a D peptide never shares a folder with an L one, even ignoring case")
check(mmforms.letters("ETF&SDLW&KLL") == "ETFSDLWKLL"
      and mmforms.d_positions("ETF&SDLW&KLL") == [4, 8],
      "the topology is built from the L letters, with the D positions noted")
check(jobmod.default().get(mmforms.JOB_KEY) is True,
      "a new run rescores D residues and closed rings")
# The shipped example as a job file from before the setting: without the key.
old_example = json.load(open(os.path.join(ROOT, "adcp_pilot", "example_run",
                                          "job.json")))
old_example.pop(mmforms.JOB_KEY, None)
check(not mmforms.enabled(old_example),
      "the shipped example, as made before, reads as a run that skips them")

# --- what is not scored, and why --------------------------------------------
new = jobmod.default(cyclic="backbone")
old = dict(new)
old.pop(mmforms.JOB_KEY)
check(mmforms.unscoreable("ETFSDLWKLLPEN", old) == "cyclic peptide",
      "an old ring run still skips its rings, in its own words")
check(mmforms.unscoreable("&TNWYANLEKLLR", dict(old, cyclic="none"))
      == "non-standard or D residues",
      "an old run still skips D residues, in its own words")
check(mmforms.unscoreable("ETFSDLWKLLPEN", new) is None
      and mmforms.unscoreable("&TNWYANLEKLLR", new) is None,
      "a new run scores rings and D residues")
check(mmforms.unscoreable("N<BHD>PIS&DVD", new) == "non-standard residues",
      "non-standard side chains are still not scored")
check(mmforms.engine_refusal("amber", "&TNWY", "backbone") is None
      and mmforms.engine_refusal("openmm", "&TNWY", "disulfide") is None,
      "AmberTools and OpenMM, which read tleap's topology, score both")
r_prime = mmforms.engine_refusal("schrodinger", "ETFSDLWKLLPEN", "backbone")
check(r_prime and "Prime" in r_prime and "ring" in r_prime,
      "Prime refuses a ring, with the reason (%s)" % r_prime)
check(mmforms.engine_refusal("schrodinger", "&TNWY", "none")
      and mmforms.engine_refusal("schrodinger", "TNWYA", "none") is None,
      "and a D peptide, and nothing else")

# --- the closure test -------------------------------------------------------
check(mmforms.BACKBONE_CLOSED == 2.5 and mmforms.DISULFIDE_CLOSED == 3.5,
      "the cutoffs are the ones the documentation states (%s, %s)"
      % (mmforms.BACKBONE_CLOSED, mmforms.DISULFIDE_CLOSED))
docs = open(os.path.join(ROOT, "adcp_pilot", "docs", "pages",
                         "rescoring.html"), encoding="utf-8").read()
check("2.5 &Aring;" in docs and "3.5 &Aring;" in docs,
      "and the rescoring page states both")


def ring_pose(gap, n=6):
    """A pose of n residues whose first N and last C are `gap` apart."""
    names, resn, xyz = [], [], []
    for k in range(1, n + 1):
        for nm in ("N", "CA", "C", "O"):
            names.append(nm)
            resn.append(k)
            xyz.append([3.8 * k, {"N": 0, "CA": 1, "C": 2, "O": 3}[nm], 0.0])
    xyz = np.array(xyz, dtype=float)
    first_n = 0
    last_c = len(names) - 2
    xyz[last_c] = xyz[first_n] + np.array([gap, 0.0, 0.0])
    return {"xyz": xyz, "resnum": np.array(resn), "resname": ["ALA"] * len(names),
            "name": names}


for gap, want in ((1.33, True), (2.5, True), (2.51, False), (4.3, False)):
    c = mmforms.closure(ring_pose(gap), "AAAAAA", "backbone")
    check(c["closed"] is want,
          "head to tail ends %.2f A apart are %s" % (gap, "closed" if want
                                                      else "open"))
c = mmforms.closure(ring_pose(4.3), "AAAAAA", "backbone")
check(c["why"] == "ring not closed in this pose: ends 4.30 A apart",
      "an open pose says how far apart its ends are (%s)" % c["why"])
check(mmforms.closure(ring_pose(1.4), "AAAAAA", "backbone")["pairs"] == [[1, 6]],
      "a closed head to tail ring joins residue 1 to the last")


def ss_pose(sg):
    """Cysteines at the positions given, with their SG at the points given."""
    names, resn, xyz, rn = [], [], [], []
    for k in range(1, 9):
        for nm in ("N", "CA", "C", "O") + (("CB", "SG") if k in sg else ()):
            names.append(nm)
            resn.append(k)
            rn.append("CYS" if k in sg else "ALA")
            xyz.append(sg[k] if nm == "SG" else [3.8 * k, 0.0, 0.0])
    return {"xyz": np.array(xyz, dtype=float), "resnum": np.array(resn),
            "resname": rn, "name": names}


c = mmforms.closure(ss_pose({1: [0, 0, 0], 8: [2.1, 0, 0]}), "CAAAAAAC",
                    "disulfide")
check(c["closed"] and c["pairs"] == [[1, 8]] and c["gap"] == 2.1,
      "sulfurs 2.1 A apart close a disulfide (%s)" % c)
c = mmforms.closure(ss_pose({1: [0, 0, 0], 8: [3.5, 0, 0]}), "CAAAAAAC",
                    "disulfide")
check(c["closed"], "sulfurs 3.5 A apart close it (%s)" % c)
c = mmforms.closure(ss_pose({1: [0, 0, 0], 8: [3.51, 0, 0]}), "CAAAAAAC",
                    "disulfide")
check(not c["closed"] and c["why"] == "ring not closed in this pose: sulfurs "
      "3.51 A apart", "sulfurs 3.51 A apart do not (%s)" % c["why"])
c = mmforms.closure(ss_pose({1: [0, 0, 0], 4: [2.0, 0, 0], 8: [2.3, 0, 0]}),
                    "CAACAAAC", "disulfide")
check(c["closed"] and c["pairs"] == [[4, 8]],
      "with three cysteines the closest pair is joined, each sulfur once (%s)"
      % c["pairs"])
check(mmforms.closure(ss_pose({1: [0, 0, 0], 8: [2.1, 0, 0]}), "CAAAAAA&C",
                      "disulfide")["closed"],
      "a D cysteine can close a disulfide")

# --- a D residue reaches tleap under its L name ------------------------------
pL = pose_of("TWIKN")
pD = pose_of("TWIKN", mirror=True)
check(all(ca_sign(pL, k) > 0 for k in range(1, 6))
      and all(ca_sign(pD, k) < 0 for k in range(1, 6)),
      "the reflected pose is the all D peptide")
molD = poseviewer.transplant("&T&W&I&K&N", pD)
check(molD is not None, "a D pose is laid on the L template atom for atom")
if molD is not None:
    mmforms.mark(molD, "&T&W&I&K&N", None)
    spec = mmgbsa.peptide_spec(molD)
    check(spec["d"] == [1, 2, 3, 4, 5] and spec["ring"] is None,
          "the pose carries its D positions to the topology builder (%s)" % spec)
    rn = [a.GetMonomerInfo().GetResidueName().strip()
          for a in molD.GetAtoms() if a.GetMonomerInfo().GetName().strip() == "CA"]
    check(rn == ["THR", "TRP", "ILE", "LYS", "ASN"],
          "D residues are named as their L forms, as ADCP 1.1 writes them (%s)" % rn)
    cip = [a.GetProp("_CIPCode") for a in molD.GetAtoms()
           if a.GetMonomerInfo().GetName().strip() == "CA"
           and a.HasProp("_CIPCode")]
    check(cip == ["R", "R", "R", "R", "R"],
          "and their stereo tags are read from the coordinates: D-Thr, D-Trp, "
          "D-Ile, D-Lys and D-Asn are R at the alpha carbon (%s)" % cip)
check(poseviewer.transplant("N<BHD>PIS", pD) is None,
      "a non-standard residue has no template, and gets None, not an error")
molL = poseviewer.transplant("TWIKN", pL)
check(molL is not None and not molL.HasProp(mmgbsa.D_PROP)
      and mmgbsa.peptide_spec(mmforms.mark(molL, "TWIKN", None))["d"] == [],
      "an L peptide carries nothing and is built as it always was")

# --- the tleap input for each ring -------------------------------------------
wd = os.path.join(TMP, "ringtext")
os.makedirs(wd)
recl = pdb_lines(pose_of("AAAAA", shift=(30, 0, 0)), "A", first=1)
pep = pose_of("CAGWC")
pepl = pdb_lines(pep, "B", first=6)
com = os.path.join(wd, "x_com_forleap.pdb")
open(com, "w").writelines(recl + ["TER\n"] + pepl + ["TER\n", "END\n"])
text = mmgbsa._ring_leap(com, {"kind": "backbone", "pairs": [[1, 5]]}, wd, "x",
                         "mbondi2")
check("loadPdbUsingSeq x_pep_forleap.pdb { CYS ALA GLY TRP CYS }" in text
      and "bond pep.1.N pep.5.C" in text and "combine { rec pep }" in text,
      "a head to tail ring is read with internal templates and joined N to C:\n%s"
      % text)
check(not any(" OXT " in l for l in open(com) if l[21:22] == "B"),
      "and the peptide has no OXT")
check(open(os.path.join(wd, "x_pep_forleap.pdb")).readline()[22:26] == "   1",
      "its peptide is numbered from 1 for the bond to name")
open(com, "w").writelines(recl + ["TER\n"] + pepl + ["TER\n", "END\n"])
text = mmgbsa._ring_leap(com, {"kind": "disulfide", "pairs": [[1, 5]]}, wd,
                         "x", "mbondi2")
names = [l[17:20] for l in open(com) if l.startswith("ATOM") and l[21] == "B"
         and l[12:16].strip() == "CA"]
check("bond pep.1.SG pep.5.SG" in text and "loadpdb x_pep_forleap.pdb" in text
      and names[0] == "CYX" and names[-1] == "CYX",
      "a disulfide names its cysteines CYX and joins the sulfurs (%s)" % names)
again = mmgbsa._ring_leap(com, {"kind": "disulfide", "pairs": [[1, 5]]}, wd,
                          "x", "mbondi2")
check(again == text, "and making the ring twice changes nothing")
bad = None
try:
    mmgbsa._ring_leap(com, {"kind": "backbone", "pairs": [[1, 3]]}, wd, "x",
                      "mbondi2")
except ValueError as exc:
    bad = str(exc)
check(bad and "joins residue 1 to 5" in bad,
      "a head to tail ring that does not join the ends is refused (%s)" % bad)

# --- an open ring is refused, with the reason --------------------------------
run = os.path.join(TMP, "openring")
pdir = os.path.join(run, "docking", "ETFSDLWKL")
os.makedirs(pdir)
far = pose_of("ETFSDLWKL")
for k in (1, 2, 3):
    with open(os.path.join(pdir, "out_ranked_%d.pdb" % k), "w") as fh:
        fh.writelines(pdb_lines(far, "A"))
j = jobmod.default(cyclic="backbone", mmgbsa=True, mmgbsa_poses=3)
out = {"modes": [{"rank": 1}, {"rank": 2}, {"rank": 3}]}
logs = []
pipeline._stage_mmgbsa(run, j, "ETFSDLWKL", out, pdir, 3, 1,
                       pipeline.docking_options(j), log=logs.append)
check(out.get("mmgbsa") == [] and out.get("mmgbsa_skipped")
      == "ring not closed in any of its top 3 poses",
      "a ring with no closed pose is not rescored, and says why (%s)"
      % out.get("mmgbsa_skipped"))
ps = out.get("mmgbsa_poses_skipped") or []
check(len(ps) == 3 and all(re.match(r"ring not closed in this pose: ends "
                                    r"\d+\.\d\d A apart$", p["why"]) for p in ps),
      "each pose says how far apart its ends are (%s)" % ps)
check(not os.path.isdir(os.path.join(run, "mmgbsa")),
      "and nothing was built for it, linear or otherwise")
row = report.summarise("ETFSDLWKL", [{"score": -10.0}], [],
                       mm_skipped=out["mmgbsa_skipped"], name="ring")
row2 = report.summarise("TSFAEYWNLLS", [{"score": -9.0}], [-40.0, -38.0],
                        name="lin")
txt = report.table([row, row2])
check("Unranked (no mmgbsa_top5): ring (ring not closed in any of its top 3 "
      "poses)." in txt, "results.txt names the peptide and its reason:\n%s"
      % txt.splitlines()[-1])

# --- an old job, resumed --------------------------------------------------------
old_j = dict(j)
old_j.pop(mmforms.JOB_KEY)
out = {"modes": [{"rank": 1}]}
logs = []
pipeline._stage_mmgbsa(run, old_j, "ETFSDLWKL", out, pdir, 3, 1,
                       pipeline.docking_options(old_j), log=logs.append)
check(out == {"modes": [{"rank": 1}], "mmgbsa": [],
              "mmgbsa_skipped": "cyclic peptide"},
      "an old ring run resumed records what it always did (%s)" % out)
check(any("the rescoring builds a linear peptide" in l for l in logs),
      "and logs what it always did (%s)" % logs)
out = {"modes": [{"rank": 1}]}
pipeline._stage_mmgbsa(run, dict(old_j, cyclic="none"), "&TNWYANLEKLLR", out,
                       pdir, 3, 1, {"cyclic": "none"}, log=logs.append)
check(out.get("mmgbsa_skipped") == "non-standard or D residues",
      "an old run resumed still skips a D peptide (%s)" % out)
ex = os.path.join(ROOT, "adcp_pilot", "example_run")
r_ex = json.load(open(os.path.join(ex, "results", "ETFSDLWKLLPEN.json")))
t_ex = pipeline.energy_terms(ex, r_ex)
check(t_ex.get("mm_total") is not None,
      "an old run's MM-GBSA folders are still found under the names they have "
      "(%s)" % t_ex)

# --- with Amber ----------------------------------------------------------------
home = mmgbsa.amber_home()
sander = prep.find_tool("sander")
if not home or not sander:
    print("skip  AmberTools is not here; the rest was checked (%d checks)" % ran)
    if fails:
        print("test_mmgbsa_forms FAILED (%d of %d)" % (len(fails), ran))
        sys.exit(1)
    sys.exit(77)
env = dict(os.environ, AMBERHOME=home)


def tleap(wd, text):
    f = os.path.join(wd, "t.leap")
    open(f, "w").write(text)
    subprocess.run([os.path.join(home, "bin", "tleap"), "-f", f], cwd=wd,
                   capture_output=True, text=True, env=env)


SP = """single point
&cntrl
  imin=1, maxcyc=1, ncyc=1, ntb=0, igb=5, cut=999.0, ntpr=1,
/
"""
MIN = """minimization
&cntrl
  imin=1, maxcyc=300, ncyc=150, ntb=0, igb=5, cut=999.0, ntpr=300,
/
"""


def sander_energy(wd, parm, crd, tag, inp=SP):
    """(total, {term: value}) of the last energy block sander prints."""
    open(os.path.join(wd, tag + ".in"), "w").write(inp)
    subprocess.run([sander, "-O", "-i", tag + ".in", "-p", parm, "-c", crd,
                    "-o", tag + ".out", "-r", tag + ".rst7"], cwd=wd,
                   capture_output=True, text=True, env=env)
    txt = open(os.path.join(wd, tag + ".out")).read()
    if "FINAL RESULTS" in txt:
        txt = txt[txt.index("FINAL RESULTS"):]
    i = txt.index("   NSTEP       ENERGY")
    block = txt[i:i + 700]
    terms = {k.strip(): float(v) for k, v in
             re.findall(r"([A-Z0-9\-]+(?: [A-Z]+)?)\s*=\s*(-?\d+\.\d+)", block)}
    return float(block.splitlines()[1].split()[1]), terms


def build(tag, pose, d):
    """The peptide alone, built the way MM-GBSA builds it: chain B, tleap's
    LEAP input, and the D residues' hydrogens placed by mirror_hydrogens."""
    wd = os.path.join(TMP, tag)
    os.makedirs(wd)
    com = os.path.join(wd, tag + "_com_forleap.pdb")
    # Without its OXT, as ADCP writes a pose: tleap adds it, and on a D
    # residue mirror_hydrogens places it.
    open(com, "w").writelines(pdb_lines(pose, "B", drop=("OXT",))
                              + ["TER\n", "END\n"])
    tleap(wd, mmgbsa.LEAP.format(complex_pdb=os.path.basename(com),
                                 prefix=tag, radii="mbondi2",
                                 receptor_pdb="", ligand_pdb=""))
    if d:
        shutil.copy(os.path.join(wd, tag + "_com.rst7"),
                    os.path.join(wd, tag + "_tleap_only.rst7"))
        got = mmgbsa.mirror_hydrogens(wd, tag, com, {"ring": None, "d": d}, 0,
                                      home, env, "mbondi2", log=print)
        check(got and got["atoms_placed"] > 0,
              "mirror_hydrogens placed the D residues' hydrogens (%s)" % got)
    return wd, tag + "_com.parm7", tag + "_com.rst7"


wdL, parmL, crdL = build("pepL", pL, [])
wdD, parmD, crdD = build("pepD", pD, [1, 2, 3, 4, 5])
P = af.read_parm7(os.path.join(wdL, parmL))
check(af.dihedral_phases(P) == [0.0, 180.0],
      "every torsion in the built topology has a phase of 0 or 180 degrees (%s)"
      % af.dihedral_phases(P))
check("CMAP" not in open(os.path.join(wdL, parmL)).read(),
      "and the topology has no CMAP term")
eL, tL = sander_energy(wdL, parmL, crdL, "spL")
eD, tD = sander_energy(wdD, parmD, crdD, "spD")
e0, t0 = sander_energy(wdD, parmD, "pepD_tleap_only.rst7", "sp0")
print("  GB energy alone: L %.4f, its mirror image as D %.4f, D with tleap's "
      "own hydrogens %.4f kcal/mol" % (eL, eD, e0))
check(abs(eL - eD) < 0.01 and all(abs(tL[k] - tD[k]) < 0.01 for k in tL),
      "a D peptide and its mirror image as an L peptide have the same GB "
      "energy, term by term (L %.4f, D %.4f)" % (eL, eD))
check(abs(e0 - eL) > 1.0,
      "while tleap's own hydrogens on a D residue would not (%.1f against %.1f)"
      % (e0, eL))
mL, _ = sander_energy(wdL, parmL, crdL, "minL", MIN)
mD, _ = sander_energy(wdD, parmD, crdD, "minD", MIN)
print("  after 300 steps of minimization: L %.4f, D %.4f kcal/mol" % (mL, mD))
check(abs(mL - mD) < 0.05,
      "and so do the two minimized the same way (L %.4f, D %.4f)" % (mL, mD))
xD = af.read_rst7(os.path.join(wdD, "minD.rst7"))
hand = af.handedness(af.stereo_centres(af.read_parm7(os.path.join(wdD, parmD)),
                                       [1, 2, 3, 4, 5]), xD)
check(all(h["volume"] < 0 for h in hand if h["centre"] == "CA")
      and all(h["h_opposite"] for h in hand),
      "every D center is still D after the minimization, its hydrogen opposite "
      "(%s)" % hand)
check(len([h for h in hand if h["centre"] == "CB"]) == 2,
      "the beta carbons of Thr and Ile are checked as well")

# Both rings, built by tleap from what _ring_leap writes.
for kind, pairs in (("backbone", [[1, 5]]), ("disulfide", [[1, 5]])):
    wdr = os.path.join(TMP, "ring_" + kind)
    os.makedirs(wdr)
    com = os.path.join(wdr, "r_com_forleap.pdb")
    open(com, "w").writelines(recl + ["TER\n"] + pepl + ["TER\n", "END\n"])
    tleap(wdr, mmgbsa._ring_leap(com, {"kind": kind, "pairs": pairs}, wdr, "r",
                                 "mbondi2"))
    pp = os.path.join(wdr, "r_com.parm7")
    check(os.path.isfile(pp), "tleap builds the %s ring" % kind)
    if not os.path.isfile(pp):
        continue
    PR = af.read_parm7(pp)
    xr = af.read_rst7(os.path.join(wdr, "r_com.rst7"))
    res = af.residues(PR)
    a, b = res[5], res[9]
    if kind == "backbone":
        t = af.bond_term(PR, xr, a[2]["N"], b[2]["C"])
        check(t["present"] and t["r0"] == 1.335,
              "the head to tail ring has its C-N bond (%s)" % t)
        check("OXT" not in b[2] and "H1" not in a[2] and "H" in a[2]
              and [x[1] for x in res[5:]] == ["CYS", "ALA", "GLY", "TRP", "CYS"],
              "and no termini: %s" % [(x[1], sorted(x[2])) for x in (a, b)])
    else:
        t = af.bond_term(PR, xr, a[2]["SG"], b[2]["SG"])
        check(t["present"] and abs(t["r0"] - 2.038) < 1e-3,
              "the disulfide has its S-S bond (%s)" % t)
        check(a[1] == "CYX" and b[1] == "CYX" and "HG" not in a[2]
              and "H1" in a[2] and "OXT" in b[2],
              "between two CYX, which here keep the peptide's termini "
              "(%s %s, %s %s)" % (a[1], sorted(a[2]), b[1], sorted(b[2])))
    q = sum(af.read_parm7(pp, ("CHARGE",))["CHARGE"]) / 18.2223
    check(abs(q - round(q)) < 1e-3, "the %s ring has an integer charge (%.4f)"
          % (kind, q))

if fails:
    print("test_mmgbsa_forms FAILED (%d of %d)" % (len(fails), ran))
    sys.exit(1)
print("ok  D residues and closed rings in MM-GBSA, old jobs kept, %d checks"
      % ran)
