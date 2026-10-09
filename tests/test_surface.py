"""Molecular surfaces: the MSMS call, its inputs, its cache and its fallback.

By Aadhil Haq

MSMS ships in ADFRsuite's bin folder and is called as the installed program.
What is checked: the radii it is given are pdb_to_xyzr's, read from
ADFRsuite's own copy of that table; the triangulation it writes comes back
whole and decodes to the counts it states; the second request for the same
structure is read from the cache; and every way it can fail (no msms, an msms
that writes nothing, a structure too large) comes back as a reason, which the
server turns into the viewer drawing its own surface.

    python3 tests/test_surface.py
"""
import base64
import json
import os
import struct
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import server, structure, surface, tools as toolmod  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


TMP = tempfile.mkdtemp(prefix="adcp_surface_")
os.environ["HOME"] = TMP
REC = os.path.join(ROOT, "adcp_dock", "example_run", "receptor", "1YCR.pdb")
# Where the tools are, for the checks that run them: ADCP_TEST_TOOLS, or else the folder
# holding ADFRsuite. Without either those checks are skipped.
ROOTS = [r for r in (os.environ.get("ADCP_TEST_TOOLS", ""),
                     os.path.dirname(os.environ.get("ADFRSUITE_HOME", "").rstrip("/")))
         if r and os.path.isdir(r)]
atoms = [a for a in structure.read(REC) if a.chain == "A"]

# --- the density follows the size -------------------------------------------
check(surface.density_for(700) == 3.0 and surface.density_for(6000) == 1.0
      and surface.density_for(60000) == 0.3,
      "the density does not fall as the structure grows: %s"
      % [surface.density_for(n) for n in (700, 6000, 60000)])

# --- hydrogens are left out, as pdb_to_xyzr leaves them ---------------------
xyzr, keys, which = surface.spheres(atoms)
heavy = [a for a in atoms if not surface._is_h(a)]
check(len(xyzr.strip().splitlines()) == len(heavy) == len(which),
      "one sphere per heavy atom, none for hydrogens")
check("A:25" in keys and "A:96" in keys,
      "each sphere is tied to its residue (%s ...)" % keys[:3])

# --- the fallbacks, every one a reason ---------------------------------------
for exe, why in ((None, "not installed"),
                 (os.path.join(TMP, "no_such_msms"), "not installed")):
    try:
        surface.compute(atoms, exe)
        check(False, "a missing msms gives a surface")
    except surface.SurfaceError as exc:
        check(why in str(exc), "a missing msms says so (%s)" % exc)

mute = os.path.join(TMP, "bin", "msms")
os.makedirs(os.path.dirname(mute))
open(mute, "w").write("#!/bin/sh\necho 'MSMS 2.6.1 started'\necho "
                      "'cannot triangulate' >&2\nexit 3\n")
os.chmod(mute, 0o755)
try:
    surface.compute(atoms, mute)
    check(False, "an msms that wrote nothing gives a surface")
except surface.SurfaceError as exc:
    check("wrote no surface" in str(exc) and "cannot triangulate" in str(exc),
          "an msms that wrote nothing is reported with its last words (%s)" % exc)

was = surface.MAX_ATOMS
surface.MAX_ATOMS = 10
try:
    surface.compute(atoms, mute)
    check(False, "a structure over the limit gives a surface")
except surface.SurfaceError as exc:
    check("more than this server triangulates" in str(exc),
          "a structure over the limit is refused before msms runs (%s)" % exc)
surface.MAX_ATOMS = was

# The server turns any of these into the viewer's own surface, with the reason.
runs = os.path.join(TMP, "runs")
app = server.App(runs)
for k in ("ADFRSUITE_HOME",):
    os.environ.pop(k, None)
app.tools_roots = []
_was_find = toolmod.adfr_program
toolmod.adfr_program = lambda name, tools=None, roots=(): None
got = app.surface({"path": [REC], "chains": ["A"]})
toolmod.adfr_program = _was_find
check(got.get("source") == "fallback" and "not installed" in got.get("why", ""),
      "without msms the page is told to draw its own surface, and why (%s)" % got)
got = app.surface({"path": [os.path.join(TMP, "nothing.pdb")]})
check(got.get("source") == "fallback", "a structure that is not there is a fallback")
got = app.surface({"dir": [os.path.join(TMP, "not_a_run")], "what": ["receptor"]})
check(got.get("source") == "fallback" and "no such run" in got.get("why", ""),
      "a run that is not there is a fallback (%s)" % got)

# --- the page draws the mesh, and its own surface when there is none ------------
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
page = open(os.path.join(ROOT, "adcp_dock", "page.py")).read()
check("viewer.addCustom(spec)" in appjs, "the page draws MSMS's triangles as a mesh")
check(appjs.count("$3Dmol.SurfaceType.VDW") >= 2,
      "the site view and the pose viewer both fall back to 3Dmol's own surface")
check("Surface drawn by the viewer" in appjs and "surface drawn by the viewer" in appjs,
      "a fallback is said under the view, quietly")
check("id=sitesurf" in page and "id=surfnote" in page and "id=surfmsg" in page,
      "the site view has its surface tick, and both views a line for the note")
check("let RSURF = 2" in appjs and 'class="pill on" id=rsurf' in page,
      "a pose's page draws the receptor's surface from the start, see-through")

# --- the real msms, where ADFRsuite has one ------------------------------------
exe = toolmod.discover(roots=ROOTS).adfr_program("msms")
if not exe:
    print("skip  msms is not installed here; the fallbacks were checked (%d checks)"
          % ran)
    if fails:
        for f in fails:
            print("  FAIL " + f)
        sys.exit(1)
    sys.exit(77)

table_path = surface.radii_file(surface.adfr_home(exe))
check(table_path and table_path.endswith(os.path.join("MolKit", "radii_patterns.py")),
      "pdb_to_xyzr's table is found in ADFRsuite (%s)" % table_path)
table = surface.read_radii_table(table_path)
# pdb_to_xyzr's united atom radii: a methyl carbon 2.00, a carbonyl oxygen
# 1.40, a carbonyl carbon 1.74, an Arg NH 1.80, a zinc 1.40.
for res, name, want in (("ALA", "CB", 2.00), ("GLY", "O", 1.40), ("GLY", "C", 1.74),
                        ("ARG", "NH1", 1.80), ("ZN", "ZN", 1.40)):
    r = surface.radius_of(res, name, "", table)
    check(abs(r - want) < 1e-6, "%s %s takes radius %.2f, not pdb_to_xyzr's %.2f"
          % (res, name, r, want))
check(surface.version(exe).startswith("MSMS"), "msms names itself (%r)" % surface.version(exe))

cache = os.path.join(TMP, "cache")
got = surface.cached(atoms, exe, cache)
check(got["source"] == "msms" and not got["cached"], "the first request runs msms")
check(got["probe"] == 1.5 and got["density"] == 3.0,
      "a 1.5 A probe, and the density for a small receptor (%s, %s)"
      % (got["probe"], got["density"]))
nv, nf = got["vertices"], got["triangles"]
verts = struct.unpack("<%df" % (3 * nv), base64.b64decode(got["verts"]))
norms = struct.unpack("<%df" % (3 * nv), base64.b64decode(got["norms"]))
faces = struct.unpack("<%dI" % (3 * nf), base64.b64decode(got["faces"]))
res = struct.unpack("<%dI" % nv, base64.b64decode(got["res"]))
check(nv > 5000 and nf > 2 * nv - 100,
      "MDM2's surface is a closed mesh of thousands of triangles (%d, %d)" % (nv, nf))
check(max(faces) < nv and min(faces) >= 0, "every triangle names vertices that exist")
check(all(abs(sum(norms[i + k] ** 2 for k in range(3)) - 1) < 0.02
          for i in range(0, 3 * nv, 3 * 97)), "the normals are unit vectors")
check(max(res) < len(got["residues"]), "every vertex names a residue that exists")
xs = [a.x for a in atoms]
check(min(xs) - 6 < min(verts[0::3]) and max(verts[0::3]) < max(xs) + 6,
      "the surface lies over the atoms, not elsewhere")
check("pdb_to_xyzr" in got["radii"], "the answer says which radii it used")
again = surface.cached(atoms, exe, cache)
check(again["cached"] and again["vertices"] == nv, "the second request is read from the cache")
check(len(os.listdir(cache)) == 1, "one cache file per structure")

# Through the server, for the site view and for a run.
os.environ["ADFRSUITE_HOME"] = surface.adfr_home(exe)
app = server.App(runs, tools_roots=ROOTS)
got = app.surface({"path": [REC], "chains": ["A"]})
check(got.get("source") == "msms" and got.get("vertices") == nv,
      "the site view's surface of chain A is the same surface (%s)"
      % {k: got.get(k) for k in ("source", "vertices", "why")})
both = app.surface({"path": [REC], "chains": ["A,B"]})
check(both.get("vertices", 0) != nv, "the chains asked for are the chains surfaced")
check(os.path.isdir(os.path.join(runs, "_surfaces")),
      "the server keeps its cache in the runs folder")
json.dumps(got)

if fails:
    print("test_surface FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  MSMS surfaces, their radii, cache and fallback, %d checks" % ran)
