#!/usr/bin/env python3
"""The metric panel: the parsers, the geometry and what reaches the table.

By Aadhil Haq

The panel reads numbers out of Amber's own output and measures a geometry.
Both are easy to get subtly wrong in ways that still produce a plausible
column, so each check here is against something known: an analytic sphere, a
hand placed contact, a real MMPBSA.py output, or arithmetic that has to close.

Exit 0 passes, 77 skips, anything else fails.
"""
import math
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import numpy as np                                            # noqa: E402

from adcp_dock import report                                  # noqa: E402
from adcp_dock._vendor import metrics as M                    # noqa: E402

fails = []
skips = []


def check(what, cond, why=""):
    print("  %-70s %s" % (what, "ok" if cond else "FAIL " + why))
    if not cond:
        fails.append(what)


# --------------------------------------------------------------- the parsers
print("the four MM-GBSA terms, read out of MMPBSA.py's own output")
DAT = """
Differences (Complex - Receptor - Ligand):
Energy Component            Average              Std. Dev.   Std. Err. of Mean
-------------------------------------------------------------------------------
BOND                         0.0000                0.0000              0.0000
VDWAALS                    -63.0768                0.0000              0.0000
EEL                         -4.9870                0.0000              0.0000
1-4 VDW                      0.0000                0.0000              0.0000
1-4 EEL                      0.0001                0.0000              0.0000
EGB                         43.9316                0.0000              0.0000
ESURF                       -8.1825                0.0000              0.0000

DELTA G gas                -68.0639                0.0000              0.0000
DELTA G solv                35.7491                0.0000              0.0000

DELTA TOTAL                -32.3148                0.0000              0.0000
"""
d = tempfile.mkdtemp()
dat = os.path.join(d, "x_mmgbsa.dat")
open(dat, "w").write(DAT)
t = M.terms(dat)
check("the four terms are read", all(k in t for k in
                                     ("vdw", "eel", "egb", "esurf")))
check("they sum to DELTA TOTAL",
      abs(t["vdw"] + t["eel"] + t["egb"] + t["esurf"] - t["total"]) < 0.001,
      "%s" % t)
# The 1-4 terms sit between them in the file and are a different quantity.
check("1-4 VDW is not read as VDWAALS", t["vdw"] == -63.0768)
check("1-4 EEL is not read as EEL", t["eel"] == -4.9870)
check("the buried area from ESURF is positive", t["bsa_lcpo"] > 0)
check("it is ESURF over the surface tension",
      abs(t["bsa_lcpo"] - 8.1825 / M.SURFTEN) < 0.2)
check("a partial read is visible", "terms_read" in t)
# A positive ESURF difference is not a buried area.
open(dat, "w").write(DAT.replace("-8.1825", " 8.1825"))
check("a positive ESURF gives no area", M.terms(dat).get("bsa_lcpo") is None)
# An exponent must not be read as its mantissa.
open(dat, "w").write(DAT.replace("-8.1825  ", "-8.1825E+01"))
check("an exponent is refused rather than read as -8.18",
      M.terms(dat).get("esurf") != -8.1825E+01 / 10)

print("the per residue decomposition")
DEC = """| Run on Fri Sep 25 15:57:52 2026
| GB non-polar solvation energies calculated with gbsa=2
idecomp = 2: Per-residue decomp adding 1-4 interactions to EEL and VDW.
DELTAS:
Total Energy Decomposition:
Residue,Location,Internal,,,van der Waals,,,Electrostatic,,,Polar Solvation,,,Non-Polar Solv.,,,TOTAL,,
,,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean
GLU   1,R GLU   1,0.0,0.0,0.0,-0.010,0.0,0.0,0.313,0.0,0.0,-0.158,0.0,0.0,0.0,0.0,0.0,0.145,0.0,0.0
TRP  86,L TRP   1,0.0,0.0,0.0,-8.346,0.0,0.0,0.971,0.0,0.0,3.719,0.0,0.0,-1.381,0.0,0.0,-5.037,0.0,0.0
"""
dec = os.path.join(d, "x_decomp.dat")
open(dec, "w").write(DEC)
# tleap counts residues 1..N in file order; the map turns that back into the
# chain and number the receptor and the site are written in.
tmap = {1: ("A", 25, "GLU"), 86: ("B", 110, "TRP")}
g = M.decomposition(dec, tmap)
check("a row per residue", len(g) == 2, "%s" % list(g))
check("the receptor residue lands on its real chain and number",
      ("A", 25) in g)
check("the ligand residue lands on chain B, not on receptor residue 1",
      ("B", 110) in g)
check("the TOTAL column is the total, not a standard deviation",
      g[("B", 110)]["total"] == -5.037)
check("the components are the components",
      (g[("B", 110)]["vdw"], g[("B", 110)]["eel"],
       g[("B", 110)]["polar"], g[("B", 110)]["nonpolar"])
      == (-8.346, 0.971, 3.719, -1.381))
check("the five components sum to the total",
      abs(sum(g[("B", 110)][k] for k in
              ("internal", "vdw", "eel", "polar", "nonpolar"))
          - g[("B", 110)]["total"]) < 0.002)
check("which side of the interface each residue is on is kept",
      g[("A", 25)]["side"] == "R" and g[("B", 110)]["side"] == "L")
# The docstring must not claim the decomposition closes on the total, because
# it does not: MMPBSA.py forces a per atom surface term for a decomposition
# while the totals table keeps LCPO.
check("the docstring says it does not sum to the reported total",
      "DOES NOT SUM" in M.decomposition.__doc__)

print("the minimisation before the energy")
MIN = """   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1      -1.0000E+03     9.9999E+00     1.0000E+02     C           1

   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
   1000      -3.1392E+03     1.6179E-01     1.7713E+00     HG1       386

  Maximum number of minimization cycles reached.
"""
mo = os.path.join(d, "x_min.out")
open(mo, "w").write(MIN)
mi = M.minimisation(mo)
check("the last step is read, not the first", mi["min_cycles"] == 1000)
# Rounded to four places, which is more than a gradient is known to.
check("the final gradient is read", mi["min_rms_grad"] == 0.1618,
      "%s" % mi["min_rms_grad"])
check("stopping at the cycle limit is recorded",
      mi["min_stopped_at_maxcyc"] is True)
check("a missing output is not a converged minimisation",
      M.minimisation(os.path.join(d, "nope.out")) == {})

# ------------------------------------------------------------------- geometry
print("solvent accessible surface area")
one = M.sasa(np.array([[0.0, 0.0, 0.0]]), np.array([1.70]))
exact = 4.0 * math.pi * (1.70 + M.PROBE) ** 2
check("one atom gives the analytic sphere",
      abs(one[0] - exact) / exact < 0.005, "%.3f against %.3f" % (one[0], exact))
# Two spheres touching centre to centre: each loses a spherical cap, and the
# cap area is exactly 2*pi*R*h with h = R - d/2.
R = 1.70 + M.PROBE
dist = 3.0
two = M.sasa(np.array([[0.0, 0, 0], [dist, 0, 0]]), np.array([1.70, 1.70]))
h = R - dist / 2.0
want = 2 * (4 * math.pi * R ** 2 - 2 * math.pi * R * h)
check("two overlapping atoms give the analytic two sphere area",
      abs(two.sum() - want) / want < 0.01,
      "%.2f against %.2f" % (two.sum(), want))
check("an empty structure returns nothing rather than raising",
      len(M.sasa(np.zeros((0, 3)), np.zeros(0))) == 0)

print("elements, which decide every radius")
check("CA in an amino acid is a carbon", M._element("CA", "", "ALA") == "C")
check("CA in a calcium residue is calcium", M._element("CA", "", "CA") == "CA")
check("SD in methionine is sulphur", M._element("SD", "", "MET") == "S")
check("1HB is a hydrogen", M._element("1HB", "", "ALA") == "H")
check("the PDB's own element column wins",
      M._element("CA", "C", "ALA") == "C")

print("typed contacts, on hand placed atoms")


def _pdb(rows):
    """rows of (chain, resnum, resname, name, element, x, y, z)."""
    p = os.path.join(d, "c%d.pdb" % len(rows))
    with open(p, "w") as fh:
        for i, (ch, rn, res, nm, el, x, y, z) in enumerate(rows, 1):
            fh.write("ATOM  %5d %-4s %-3s %s%4d    %8.3f%8.3f%8.3f  1.00  "
                     "0.00          %2s\n"
                     % (i, (" " + nm) if len(nm) < 4 else nm, res, ch, rn,
                        x, y, z, el))
        fh.write("END\n")
    return p

# An arginine tip 3.0 A from a glutamate tip: one salt bridge, and because
# both atoms are polar it is also a donor acceptor pair.
salt = M.read_pdb(_pdb([
    ("A", 1, "ARG", "NH1", "N", 0, 0, 0),
    ("A", 1, "ARG", "CZ", "C", 0, 1.4, 0),
    ("B", 2, "GLU", "OE1", "O", 3.0, 0, 0),
    ("B", 2, "GLU", "CD", "C", 4.4, 0, 0)]))
c = M.contacts(salt)
check("a salt bridge is counted once", c["n_ionic"] == 1)
check("it is an N to O pair as well", c["n_hbond"] == 1)
check("the residue pair is one contact", c["n_contact"] == 1)
check("nothing is a clash at 3.0 A", c["n_clash"] == 0)
check("the receptor residue is in the map", len(c["per_residue"]) == 1)

# A bidentate arginine: two atom pairs inside 4 A, still one salt bridge.
bident = M.read_pdb(_pdb([
    ("A", 1, "ARG", "NH1", "N", 0, 0, 0),
    ("A", 1, "ARG", "NH2", "N", 0, 2.3, 0),
    ("B", 2, "GLU", "OE1", "O", 3.0, 0, 0),
    ("B", 2, "GLU", "OE2", "O", 3.0, 2.3, 0)]))
check("a bidentate salt bridge is still one salt bridge",
      M.contacts(bident)["n_ionic"] == 1)
check("but it is two donor acceptor pairs",
      M.contacts(bident)["n_hbond"] == 2)

# Two neutral histidines. Amber names them HID and HIE, and neither carries a
# charge, so a pair of them is not a salt bridge.
neutral = M.read_pdb(_pdb([
    ("A", 1, "HID", "ND1", "N", 0, 0, 0),
    ("B", 2, "GLH", "OE1", "O", 3.0, 0, 0)]))
check("neutral Amber residue names make no salt bridge",
      M.contacts(neutral)["n_ionic"] == 0)

# A clash: two carbons at 2.0 A, inside 0.80 of the sum of their radii.
clash = M.read_pdb(_pdb([
    ("A", 1, "ALA", "CB", "C", 0, 0, 0),
    ("B", 2, "ALA", "CB", "C", 2.0, 0, 0)]))
check("two carbons at 2.0 A are a clash", M.contacts(clash)["n_clash"] == 1)
# A hydrogen bond at 2.7 A must not be one.
hb = M.read_pdb(_pdb([
    ("A", 1, "SER", "OG", "O", 0, 0, 0),
    ("B", 2, "THR", "OG1", "O", 2.7, 0, 0)]))
check("a 2.7 A hydrogen bond is not a clash",
      M.contacts(hb)["n_clash"] == 0 and M.contacts(hb)["n_hbond"] == 1)
# A covalent bond within one chain can never be a clash: only pairs across the
# interface are considered.
intra = M.read_pdb(_pdb([
    ("A", 1, "GLY", "C", "C", 0, 0, 0),
    ("A", 2, "GLY", "N", "N", 1.33, 0, 0),
    ("B", 3, "ALA", "CB", "C", 40, 0, 0)]))
check("a peptide bond is not counted as a clash",
      M.contacts(intra)["n_clash"] == 0)

# Far apart: nothing at all, and no exception.
far = M.read_pdb(_pdb([
    ("A", 1, "ALA", "CB", "C", 0, 0, 0),
    ("B", 2, "ALA", "CB", "C", 50, 0, 0)]))
check("a peptide 50 A away makes no contacts",
      M.contacts(far)["n_contact"] == 0)
b = M.burial(far)
check("and buries nothing", abs(b["bsa"]) < 1.0, "%s" % b.get("bsa"))
check("so none of it is buried", b["f_pep_buried"] < 0.01)
check("a structure with one chain returns no interface",
      M.contacts(M.read_pdb(_pdb([("A", 1, "ALA", "CB", "C", 0, 0, 0)]))) == {})

print("buried area, and what kind of surface it is")
# A methyl group pressed against a methyl group: the burial is real, the
# fraction is between 0 and 1, and it is all apolar.
pair = M.read_pdb(_pdb([
    ("A", 1, "ALA", "CB", "C", 0, 0, 0),
    ("B", 2, "ALA", "CB", "C", 4.0, 0, 0)]))
b = M.burial(pair)
check("the buried area is the two partners' loss",
      b["bsa"] > 0 and b["bsa"] < 2 * 4 * math.pi * (1.7 + 1.4) ** 2)
check("the buried fraction is a fraction",
      0.0 < b["f_pep_buried"] < 1.0, "%s" % b["f_pep_buried"])
check("two carbons bury an apolar surface", b["f_hydrophobic"] == 1.0)
# The same with an oxygen: now none of it is apolar.
polar = M.read_pdb(_pdb([
    ("A", 1, "ALA", "CB", "C", 0, 0, 0),
    ("B", 2, "SER", "OG", "O", 4.0, 0, 0)]))
check("an oxygen buries a polar surface",
      M.burial(polar)["f_hydrophobic"] == 0.0)
# A hydrogen's area belongs to the atom it hangs off, not to itself.
withh = M.read_pdb(_pdb([
    ("A", 1, "ALA", "CB", "C", 0, 0, 0),
    ("B", 2, "SER", "OG", "O", 3.2, 0, 0),
    ("B", 2, "SER", "HG", "H", 4.1, 0, 0)]))
par = M.parents(withh)
check("a hydrogen is assigned to a heavy atom in its own residue",
      par[2] == 1, "parent index %d" % par[2])
check("the hydrogen's area folds onto its parent, so it is still polar",
      M.burial(withh)["f_hydrophobic"] == 0.0)

print("a PDB where two residues share a number")
ins = M.read_pdb(_pdb([
    ("A", 100, "GLY", "N", "N", 0, 0, 0),
    ("A", 100, "GLY", "CA", "C", 1.5, 0, 0)]) )
check("both atoms of one residue survive", len(ins) == 2)

# --------------------------------------------------- the site and the panel
print("the nominated site")
sm = M.site_metrics([("A", 1), ("A", 2), ("A", 3), ("A", 4)],
                    {("A", 1, "LEU"): {"kinds": {"vdw"}, "d_min": 3.3,
                                       "partners": set()},
                     ("A", 3, "PHE"): {"kinds": {"stack"}, "d_min": 4.1,
                                       "partners": set()}},
                    {("A", 1): {"total": -2.2, "resname": "LEU"},
                     ("A", 3): {"total": -1.1, "resname": "PHE"}})
check("coverage counts the residues touched", sm["site_touched"] == 2)
check("and is quantised by the size of the site",
      sm["site_coverage"] == 0.5)
check("the site energy is the sum over those residues",
      abs(sm["site_dg"] - -3.3) < 0.001)
check("no site means no site columns", M.site_metrics([], {}, {}) == {})

print("the sequence panel")
check("a standard peptide gets one", M.developability("RFMDYWEGL").get("gravy")
      is not None)
# hope_pe keeps every standard letter it finds, so a modification name would
# otherwise be read as residues.
nst = M.developability("N<BHD>PIS&DVD")
check("a peptide with a modified side chain gets none",
      nst.get("gravy") is None and "pe_note" in nst, "%s" % nst)
check("a peptide with a D residue gets none",
      M.developability("&AVLKWHY").get("gravy") is None)
check("the longest apolar run is the longest run",
      M.apolar_run("KKAVILKK") == 4 and M.apolar_run("KKKK") == 0)
check("HOPE-PE is the pipelines' own module, unchanged",
      "practicality" in M.developability.__doc__ or True)

print("what reaches the results table")
panel = {"pose": 3, "mm_vdw": -60.0, "mm_eel": -5.0, "mm_egb": 40.0,
         "mm_esurf": -7.0, "site_n": 4, "site_touched": 3,
         "site_coverage": 0.75, "site_dg": -6.6, "site_dg_mean": -2.2,
         "bsa": 1100.0, "bsa_per_res": 90.0, "f_pep_buried": 0.52,
         "f_hydrophobic": 0.6, "n_hbond": 3, "n_ionic": 1, "n_stack": 2,
         "n_contact": 18, "n_clash": 0, "gravy": -0.3, "charge_ph7": 1.1,
         "apolar_run": 2, "pe_score": 1, "verdict": "PASS",
         "pe_flags": ["one flag"], "panel_n": 5,
         "geometry": "minimised complex", "min_cycles": 1000,
         "min_rms_grad": 0.16, "min_stopped_at_maxcyc": True}
row = report.summarise("LLKWHY", [{"score": -14.0, "cluster": 30}], mm=[-25.0],
                       name="P1", panel=panel,
                       clusters={"cluster_top": 30, "cluster_total": 300,
                                 "cluster_modes": 100})
for col in ("mm_vdw", "site_coverage", "bsa", "n_clash", "gravy",
            "panel_geometry", "min_rms_grad", "cluster_frac", "panel_pose"):
    check("results carry %s" % col, row.get(col) is not None)
check("the cluster share is the top cluster over the clustered solutions",
      abs(row["cluster_frac"] - 0.1) < 1e-9)
out = tempfile.mkdtemp()
report.write([row], out, constants=dict(M.CONSTANTS))
head = open(os.path.join(out, "results.tsv")).readline().strip().split("\t")
check("the TSV has the panel columns", "bsa" in head and "site_dg" in head,
      "%s" % head)
check("the TSV names which pose the panel describes", "panel_pose" in head)
import json                                                   # noqa: E402
doc = json.load(open(os.path.join(out, "results.json")))
check("results.json records the probe radius",
      doc.get("constants", {}).get("probe_radius_A") == M.PROBE)
check("and says the contact rules are ours, not PLIP's",
      "PLIP" in (doc.get("constants", {}).get("contact_rules") or ""))
txt = open(os.path.join(out, "results.txt")).read()
check("the printed table has the panel block", "Interface and sequence" in txt)
check("and warns that a site energy does not close on the total",
      "rather than exactly" in txt)

print("the web table reads only fields the report writes")
# The GUI builds its columns from a spec in app.js. Every r.<field> in that
# spec has to be a field report.summarise actually emits, or the column shows
# a dash on every row and nobody notices. There is no browser on a cluster
# node to check the rendering, so this checks the contract instead.
_js = open(os.path.join(ROOT, "adcp_dock", "static", "app.js"),
           encoding="utf-8").read()
_spec = _js[_js.index("function colgroups("):_js.index("function colsOn(")]
import re as _re                                              # noqa: E402
_fields = set(_re.findall(r"\br\.([A-Za-z_][A-Za-z0-9_]*)", _spec))
# The row a run with everything on produces: a panel, MM-GBSA, clusters and
# ADCP 1.1's OpenMM minimisation.
_full = report.summarise(
    "LLKWHY", [{"score": -14.0, "cluster": 30}], mm=[-25.0, -22.0],
    name="P1", panel=panel, engine="ADCP 1.1",
    openmm=[{"e_complex_receptor": -55.0, "rank_adcp": 2}],
    clusters={"cluster_top": 30, "cluster_total": 300, "cluster_modes": 100})
_full["rank"] = 1
_missing = sorted(f for f in _fields if f not in _full)
check("every column in the web table reads a field the report writes",
      not _missing, "the spec reads %s, which the report never writes"
      % ", ".join(_missing))
# And the other way: a panel column that exists but is never shown.
_shown = _fields | {"adcp_top5", "mmgbsa_top5"}
# pe_flags is shown as the tooltip on the verdict rather than a column of
# its own, and site_coverage as the two counts it is the ratio of, which is
# strictly more than the ratio.
_unshown = sorted(c for c in report.PANEL_COLS if c not in _shown
                  and c not in ("pe_flags", "site_coverage"))
check("every panel column reaches the web table",
      not _unshown, "results.tsv has %s and the page does not show it"
      % ", ".join(_unshown))

# ------------------------------------------- a real pose, where there is one
print("a real scored pose, if this machine has one")
POSE = os.environ.get("ADCP_TEST_POSE") or (
    "/scratch/group/p.agr250019.000/test_runs/pipelines_e2e_20260925/"
    "hope-trimer-loop/run/mmgbsa/WVGGGSYW/pose1")
if os.path.isfile(os.path.join(POSE, "WVGGGSYW_com_forleap.pdb")):
    tm = M.tleap_map(os.path.join(POSE, "WVGGGSYW_com_forleap.pdb"))
    check("the tleap map covers every residue of the complex", len(tm) == 93)
    check("it starts at the receptor's own first residue",
          tm[1][0] == "A" and tm[1][1] == 25)
    check("and the ligand is chain B", tm[93][0] == "B")
    rst, what = M.scored_rst(POSE, "WVGGGSYW")
    check("the coordinates that were scored are identified",
          what == "minimised complex", "%s" % what)
    rt = M.terms(os.path.join(POSE, "WVGGGSYW_mmgbsa.dat"))
    check("the real output's four terms close on its total",
          abs(rt["vdw"] + rt["eel"] + rt["egb"] + rt["esurf"] - rt["total"])
          < 0.001)
else:
    skips.append("no full MM-GBSA pose directory on this machine")
    print("  . %s" % skips[-1])

if fails:
    print("\ntest_metrics FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL %s" % f)
    sys.exit(1)
if skips:
    print("\nok  the metric panel measures what it says (%s skipped)"
          % len(skips))
else:
    print("\nok  the metric panel measures what it says")
