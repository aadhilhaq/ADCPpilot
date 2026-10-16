"""Choosing AutoSite's pockets before submitting, the box clipped to its maximum,
a ligand read from the structure ADapt was opened on, and WebGL.

By Aadhil Haq

What is checked, without docking: the quick receptor AutoSite is given (the
chains kept, a ligand chain left out, waters and unkept hetero groups out, one
conformer); AutoSite's table and its pockets' points read back; the box around
chosen pockets by the rule every box follows; a target whose translation points
are replaced and nothing else touched; the prepare stage's chosen points; the
site check saying when a box is clipped and whether residues fall outside it;
a ligand read from the structure ADapt was opened on, marked when it is in
ADapt's receptor and taken out of it before the maps are built; the job's
checks with ADapt; and the page's WebGL check and pocket controls. With the
ADCP 1.1 suite installed it also runs AutoSite on the example receptor.

    python3 tests/test_pockets.py
"""
import io
import json
import os
import shutil
import sys
import tempfile
import types
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import (adcp11, job as jobmod, pipeline, pockets,  # noqa: E402
                       server, structure)

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)
        print("  FAIL " + why)


TMP = tempfile.mkdtemp(prefix="adcp_pockets_")
EX = os.path.join(ROOT, "adcp_dock", "example_run", "receptor", "1YCR.pdb")
SITE = [("A", n) for n in (54, 58, 61, 62, 67, 72, 96, 99)]

# --- the quick receptor -------------------------------------------------------
pdb = os.path.join(TMP, "mix.pdb")
with open(pdb, "w") as fh:
    fh.write(
        "ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00           N\n"
        "ATOM      2  CA AALA A   1      11.639   6.071  -5.147  0.50  0.00           C\n"
        "ATOM      3  CA BALA A   1      11.739   6.171  -5.247  0.50  0.00           C\n"
        "ATOM      4  N   GLY B   1       1.104   2.134  -1.504  1.00  0.00           N\n"
        "HETATM    5  O   HOH A 201       5.000   5.000   5.000  1.00  0.00           O\n"
        "HETATM    6 ZN    ZN A 202       6.000   6.000   6.000  1.00  0.00          ZN\n"
        "END\n")
t = pockets.receptor_text(pdb)
check("HOH" not in t and " ZN " not in t, "waters and unkept hetero groups reach AutoSite")
check(t.count("CA ") == 1 and " CA AALA" in t, "more than one conformer reaches AutoSite")
check("GLY B" in t, "with no chains named, a chain is left out")
t = pockets.receptor_text(pdb, drop=["B"])
check("GLY B" not in t, "the ligand chain the box is built around reaches AutoSite")
t = pockets.receptor_text(pdb, keep_het=True)
check(" ZN " in t and "HOH" not in t, "a kept hetero group is lost, or a water kept")
t = pockets.receptor_text(pdb, chains=["A"])
check("GLY B" not in t and "ALA A" in t, "the chains kept are not the ones AutoSite sees")
try:
    pockets.receptor_text(pdb, chains=["Q"])
    check(False, "a receptor with nothing left is not refused")
except pockets.PocketError:
    pass
k1 = pockets.key_of(path=pdb, chains=["A"], centre=[1, 2, 3])
k2 = pockets.key_of(path=pdb, chains=["A"], centre=[1, 2, 4])
check(k1 != k2 and k1 == pockets.key_of(path=pdb, chains=["A"], centre=[1, 2, 3]),
      "the cache key does not follow the box")

# --- AutoSite's table, and its pockets' points -----------------------------------
asd = os.path.join(TMP, "as")
os.makedirs(asd)
with open(os.path.join(asd, "rec_summary.csv"), "w") as fh:
    fh.write("Cluster # | Energy | #points | Radius ofN,e,v,rg,epv,buriedness,v*buriedness^2/rg\n"
             "1,-80.516000,201.000000,4.150268,-0.400577,0.830396,152.098280\n"
             "2,-43.179000,125.000000,3.224645,-0.345432,0.806020,90.454212\n")
rows = pockets._summary(os.path.join(asd, "rec_summary.csv"))
check([r["rank"] for r in rows] == [1, 2] and rows[0]["points"] == 201
      and abs(rows[0]["score"] - 152.10) < 0.01 and abs(rows[1]["buriedness"] - 0.806) < 0.001,
      "AutoSite's table is not read as it ranks the pockets: %s" % rows)
with open(os.path.join(asd, "rec_cl_001.pdb"), "w") as fh:
    fh.write("ATOM      1   C   FP     1      22.060 -22.430  -1.190  1.0   0.364    0.001 C\n"
             "ATOM      2   C   FP     1      24.060 -20.430  -3.190  1.0   0.488    0.001 C\n")
check(pockets._points(os.path.join(asd, "rec_cl_001.pdb"))
      == [[22.06, -22.43, -1.19], [24.06, -20.43, -3.19]], "a pocket's points are misread")
check(len(pockets._thin([[i, 0, 0] for i in range(5000)], 1500)) == 1500,
      "the page is sent every point of a large pocket")

# --- the box around chosen pockets -----------------------------------------------------
c, s, clip = pockets.box_around([[0, 0, 0], [10, 4, 2]], padding=8, minimum=20, maximum=30)
check(c == [5.0, 2.0, 1.0] and s == [26.0, 20.0, 20.0] and not clip,
      "the box around the pockets is not span plus padding, held between the limits: %s %s"
      % (c, s))
c, s, clip = pockets.box_around([[0, 0, 0], [40, 4, 2]], padding=8, minimum=20, maximum=30)
check(s[0] == 30.0 and clip and clip[0][0] == "x", "a box held to the maximum is not said")

# --- the chosen points replace the target's, and nothing else changes ---------------
trg = os.path.join(TMP, "receptor.trg")
import numpy                                                    # noqa: E402
old = io.BytesIO()
numpy.save(old, numpy.zeros((7, 3)))
with zipfile.ZipFile(trg, "w") as z:
    z.writestr("receptor/translationPoints.npy", old.getvalue())
    z.writestr("receptor/rigidReceptor.C.map", "map body")
    z.writestr("receptor/data.pkl", b"\x80\x02}q\x00.")
npy = pockets.save_points([[1, 2, 3], [4, 5, 6]], os.path.join(TMP, "target", "site_points.npy"))
n = pockets.swap_points(trg, npy)
with zipfile.ZipFile(trg) as z:
    got = numpy.load(io.BytesIO(z.read("receptor/translationPoints.npy")))
    others = {i: z.read(i) for i in z.namelist() if not i.endswith(".npy")}
check(n == 2 and got.shape == (2, 3) and got.dtype == numpy.float64
      and got.tolist() == [[1, 2, 3], [4, 5, 6]], "the target's points are not the chosen ones")
check(others == {"receptor/rigidReceptor.C.map": b"map body",
                 "receptor/data.pkl": b"\x80\x02}q\x00."},
      "replacing the points changed something else in the target")
with zipfile.ZipFile(os.path.join(TMP, "bare.trg"), "w") as z:
    z.writestr("receptor/rigidReceptor.C.map", "x")
try:
    pockets.swap_points(os.path.join(TMP, "bare.trg"), npy)
    check(False, "a target without points is not refused")
except pockets.PocketError:
    pass

# --- the prepare stage's chosen points ---------------------------------------------------
rd = os.path.join(TMP, "run")
os.makedirs(os.path.join(rd, "target"))
pockets.save_points([[1, 1, 1]], os.path.join(rd, "target", "site_points.npy"))
j = {"site_points": os.path.join("target", "site_points.npy")}
check(pipeline._chosen_points(rd, j) == [[1.0, 1.0, 1.0]], "the run does not find its points")
check(pipeline._chosen_points(rd, {}) is None, "a run with no choice reads points")
try:
    pipeline._chosen_points(rd, {"site_points": "target/gone.npy"})
    check(False, "missing points are docked without them")
except SystemExit:
    pass
d = jobmod.default()
check(d.get("site_pockets") == {} and d.get("site_points") == "" and d.get("site_box") is None
      and d.get("pocket_score") == "peptide" and d.get("site_ligand_from") == "",
      "a new job does not carry the pocket and ligand-source settings")

# --- the box clipped to its maximum ----------------------------------------------------------
pts = [(0, 0, 0), (18, 10, 12)]
c = structure.clipped(pts, 8.0, 30.0)
check([a for a, _s, _p in c] == ["x"] and c[0][1] == 18.0 and c[0][2] == 6.0,
      "a site whose padding is trimmed is not reported with the padding left: %s" % c)
c = structure.clipped([(0, 0, 0), (41, 3, 3)], 8.0, 30.0)
check(c and c[0][1] > 30, "a site wider than the box is not reported as such")
check(structure.clipped([(0, 0, 0), (5, 5, 5)], 8.0, 30.0) == [],
      "a small site is reported clipped")

# --- a ligand read from the structure ADapt was opened on -------------------------------
if os.path.isfile(EX):
    recA = os.path.join(TMP, "adapt_rec_H.pdb")
    with open(EX) as src, open(recA, "w") as out:
        for line in src:
            if line.startswith(("ATOM", "HETATM")) and line[21:22] == "B":
                continue
            out.write(line)
    app = server.App(os.path.join(TMP, "runs"))
    body = {"receptor": recA, "site_mode": "ligand", "site_ligand": "B",
            "box_padding": 8, "box_min": 20, "box_max": 30}
    d = app.site(body)
    check(d.get("error"), "a ligand ADapt left out was found in ADapt's receptor")
    d = app.site(dict(body, site_ligand_from=EX))
    check(not d.get("error") and d.get("centre") and d.get("ligand_in_receptor") is False,
          "a ligand ADapt left out cannot set the box from the structure ADapt opened: %s" % d)
    check(d.get("contacts"), "the ligand's contacts are not measured")
    shutil.copy(EX, os.path.join(TMP, "orig.pdb"))
    d = app.site(dict(body, receptor=EX, site_ligand_from=os.path.join(TMP, "orig.pdb")))
    check(not d.get("error") and d.get("ligand_in_receptor") is True,
          "a ligand ADapt kept in its receptor is not marked: %s" % d)
    # The job's own checks follow the ligand's source, and with ADapt the
    # chains ticked on the page do not decide.
    j = jobmod.default(receptor=recA, site_mode="ligand", site_ligand="B",
                       site_ligand_from=EX, chains=["A", "B"],
                       adapt={"folder": "/x", "save": 1})
    check(not jobmod.site_problems(j), "the job refuses a ligand read from ADapt's original: %s"
          % jobmod.site_problems(j))
    j["site_ligand_from"] = ""
    check(jobmod.site_problems(j), "the job accepts a ligand that is not in the receptor read")
    # Taken out of ADapt's receptor, in both files, and said.
    q = os.path.join(TMP, "rec.pdbqt")
    shutil.copy(EX, q)
    p2 = os.path.join(TMP, "prepared.pdb")
    shutil.copy(EX, p2)
    said = []
    n = pipeline._drop_ligand("B", (q, p2), TMP, log=said.append)
    left = open(q).read()
    check(n > 0 and not any(l[21:22] == "B" for l in left.splitlines()
                            if l.startswith(("ATOM", "HETATM"))),
          "the ligand chain stays in ADapt's receptor")
    check(any("WARNING" in x and "taken out" in x for x in said)
          and "ligand:" in open(os.path.join(TMP, "preparation.log")).read(),
          "taking the ligand out is not said in the log and preparation.log")
    check(pipeline._drop_ligand("B", (q,), TMP, log=said.append) == 0,
          "a receptor without the ligand is changed")
    # The site check reports clipping per axis, with the padding left.
    d = app.site({"receptor": EX, "site_mode": "residues",
                  "hotspots": ", ".join("%s:%d" % x for x in SITE),
                  "box_padding": 8, "box_min": 20, "box_max": 30})
    check(all(c["span"] <= 30 for c in d.get("clipped", [])) and "box_padding" in d,
          "the MDM2 site is said to be wider than the box: %s" % d.get("clipped"))
    d = app.site({"receptor": EX, "site_mode": "residues", "hotspots": "A:25-109",
                  "box_padding": 8, "box_min": 20, "box_max": 30})
    check(any(c["span"] > 30 for c in d.get("clipped", [])),
          "a whole domain is not said to be wider than the box")
else:
    print("  (no example receptor; the ADapt and clipping checks were not run)")

# --- AutoSite itself, where the ADCP 1.1 suite is installed -----------------------------
roots = [r for r in (os.environ.get("ADCP_TEST_TOOLS", ""),
                     os.path.dirname(os.environ.get("ADFRSUITE_HOME", "").rstrip("/")))
         if r and os.path.isdir(r)]
exe = adcp11.find(roots)
if exe and pockets.autosite_script(exe) and os.path.isfile(EX):
    found = pockets.find(EX, os.path.join(TMP, "cache"), exe, chains=["A"],
                         centre=[31.06, -22.43, -7.19], size=[30, 30, 30])
    check(found["pockets"] and found["pockets"][0]["rank"] == 1
          and len(found["pockets"][0]["xyz"]) == found["pockets"][0]["points"],
          "AutoSite found no pocket in the MDM2 groove, or its points were lost")
    p1 = found["pockets"][0]["centre"]
    d2 = sum((p1[i] - [31.06, -22.43, -7.19][i]) ** 2 for i in range(3)) ** 0.5
    check(d2 < 12, "the best pocket is %.1f A from the p53 groove" % d2)
    again = pockets.find(EX, os.path.join(TMP, "cache"), exe, chains=["A"],
                         centre=[31.06, -22.43, -7.19], size=[30, 30, 30])
    check(again["key"] == found["key"], "asking again is not answered from the cache")
    sel, pts, _f = pockets.chosen(os.path.join(TMP, "cache"), found["key"], [1])
    check(sel[0]["rank"] == 1 and len(pts) == found["pockets"][0]["points"],
          "the chosen pocket's points are not the ones found")
else:
    print("  (ADCP 1.1 suite not found; AutoSite itself was not run)")

# --- the page ------------------------------------------------------------------------------
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
page = open(os.path.join(ROOT, "adcp_dock", "page.py")).read()
for want, why in (("function webglOK()", "the page does not ask the browser for WebGL"),
                  ("viewerIn($('viewD')", "the site view is made without the check"),
                  ("viewerIn($('viewR')", "the pose view is made without the check"),
                  ("noGL($('viewM'))", "Mol* is started without the check"),
                  ("get.webgl.org", "the message has no way to test WebGL"),
                  ("troubleshooting.html#webgl", "the message does not lead to the fix"),
                  ("adaptnogl", "ADapt's view is not covered"),
                  ("function clippedNote", "the clipped box is not said"),
                  ("function findPockets", "pockets cannot be found from the page"),
                  ("pocketsOut()", "the chosen pockets are not sent with the run"),
                  ("clearPockets();", "pockets outlive the site they were found for"),
                  ("ligandFrom()", "the ligand's source is not sent")):
    check(want in appjs, why)
for want, why in (("id=pocketfind", "no Find pockets button"),
                  ("id=pocketscore", "no choice of how pockets are ranked"),
                  ("id=adaptnogl", "no place for ADapt's WebGL note")):
    check(want in page, why)
docs = open(os.path.join(ROOT, "adcp_dock", "docs", "pages", "troubleshooting.html")).read()
check('id="webgl"' in docs and "--use-angle=swiftshader --enable-unsafe-swiftshader" in docs,
      "the troubleshooting page has no WebGL section with the software rendering switches")
site = open(os.path.join(ROOT, "adcp_dock", "docs", "pages", "site.html")).read()
check('id="pockets"' in site and "lig.</b> column" in site,
      "the site page does not describe choosing pockets or ligands with ADapt")

shutil.rmtree(TMP, ignore_errors=True)
if fails:
    print("test_pockets FAILED (%d of %d)" % (len(fails), ran))
    sys.exit(1)
print("ok  pockets found and chosen before submitting, boxes clipped, ADapt's ligands "
      "and WebGL, %d checks" % ran)
