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

# --- the file is read whole, whatever wrote it ----------------------------------
# Past 99,999 atoms OpenMM and others write serials in hexadecimal, and past
# 9,999 residues OpenMM writes the residue number so too; an atom whose serial
# could not be read used to be dropped, so a large receptor lost its last
# atoms from everything read from it, the surface included. Coordinates may
# run together in their columns. One conformer per residue is surfaced.
odd = os.path.join(TMP, "odd.pdb")


def pdb_line(rec, serial, name, alt, res, chain, num, icode, xyz, el):
    """One fixed-column line, the fields given as the file would hold them."""
    return "%-6s%5s %-4s%1s%3s %1s%4s%1s   %s  1.00  0.00          %2s\n" % (
        rec, serial, name, alt, res, chain, num, icode, xyz, el)


def xyz(x, y, z):
    return "%8.3f%8.3f%8.3f" % (x, y, z)


open(odd, "w").write(
    "MODEL 1\n"
    + pdb_line("ATOM", "99999", " CA", "", "ALA", "A", "1", "", xyz(10, 10, 10), "C")
    + pdb_line("ATOM", "186a0", " CB", "", "ALA", "A", "1", "", xyz(11.5, 10, 10), "C")
    + pdb_line("ATOM", "*****", " N", "", "ALA", "A", "2", "", xyz(13, 10, 10), "N")
    + pdb_line("ATOM", "186a2", " CA", "", "GLY", "A", "2710", "",
               "-100.123-200.456-300.789", "C")
    + pdb_line("ATOM", "186a3", " CA", "", "GLY", "A", "27aa", "", xyz(14, 10, 10), "C")
    + pdb_line("ATOM", "10", " CA", "A", "SER", "A", "5", "", xyz(20, 20, 20), "C")
    + pdb_line("ATOM", "11", " CA", "B", "SER", "A", "5", "", xyz(20.4, 20, 20), "C")
    + pdb_line("ATOM", "12", " OG", "A", "SER", "A", "5", "", xyz(21, 20, 20), "O")
    + pdb_line("ATOM", "13", " OG", "B", "SER", "A", "5", "", xyz(21.4, 20, 20), "O")
    + pdb_line("ATOM", "14", " CA", "", "ALA", "A", "6", "A", xyz(25, 20, 20), "C")
    + pdb_line("HETATM", "15", "ZN", "", " ZN", "A", "401", "", xyz(30, 20, 20), "ZN")
    + "ENDMDL\n")
got_odd = structure.read(odd)
check(len(got_odd) == 11,
      "every atom of a file with hexadecimal serials and residue numbers is read "
      "(%d of 11)" % len(got_odd))
check(any(a.resnum == 0x27aa for a in got_odd),
      "a residue number past 9,999 written in hexadecimal is read")
wide = [a for a in got_odd if a.resnum == 2710]
check(wide and (wide[0].x, wide[0].y, wide[0].z) == (-100.123, -200.456, -300.789),
      "coordinates that run together in their columns are read from the columns")
check([a.altloc for a in got_odd if a.resnum == 5] == ["A", "B", "A", "B"],
      "the alternate location of each atom is kept")
check(any(a.icode == "A" for a in got_odd) and any(a.hetatm for a in got_odd),
      "insertion codes and hetero atoms are read")
rows_odd, keys_odd, _w = surface.sphere_rows(got_odd)
check(len(rows_odd) == 9 and "A:5" in keys_odd and "A:6A" in keys_odd,
      "one conformer per residue is surfaced, the first (%d spheres, %s)"
      % (len(rows_odd), keys_odd))
check(any(abs(r[0] - 20.0) < 1e-6 for r in rows_odd)
      and not any(abs(r[0] - 20.4) < 1e-6 for r in rows_odd),
      "and it is the A conformer, the one the viewers draw")
# The site view draws the structure as /api/rawstructure writes it back out:
# one conformer, the one surfaced, in columns that hold however many atoms.
one = list(structure.one_conformer(got_odd))
check(len(one) == 9 and not any(a.altloc == "B" for a in one),
      "the conformer drawn is the conformer surfaced")
back = os.path.join(TMP, "back.pdb")
big = structure.Atom(serial=0, name="CA", resname="ALA", chain="A", resnum=12345,
                     icode="", x=-100.123, y=-200.456, z=-300.789, element="C",
                     hetatm=False, altloc="")
structure.write_pdb(one + [big], back)
again_odd = structure.read(back)
check(len(again_odd) == 10 and (again_odd[-1].x, again_odd[-1].y, again_odd[-1].z)
      == (-100.123, -200.456, -300.789),
      "a residue number past 9,999 does not push the coordinates out of their columns "
      "(%s)" % [(a.resnum, a.x, a.y, a.z) for a in again_odd[-1:]])
check(structure._fixed(123456, 5) == "1e240" and structure._fixed(99999, 5) == "99999",
      "serials past 99,999 are written in hexadecimal, as OpenMM writes them")
shifted = os.path.join(TMP, "shifted.pdb")
open(shifted, "w").write(
    pdb_line("ATOM", "1", " CA", "", "ALA", "A", "1", "", "-1000.123   5.000   6.000", "C"))
got_sh = structure.read(shifted)
check(len(got_sh) == 1 and got_sh[0].y == 5.0,
      "a line whose writer let a wide number push the columns along is still read "
      "(%s)" % [(a.x, a.y, a.z) for a in got_sh])

# --- pieces: atoms whose surfaces join -------------------------------------------
# Two atoms join when the probe cannot pass between them: closer than their two
# radii and the probe's diameter.
p2 = 2 * surface.PROBE
near = [(0.0, 0.0, 0.0, 1.8), (1.8 + 1.8 + p2 - 0.01, 0.0, 0.0, 1.8)]
far = [(0.0, 0.0, 0.0, 1.8), (1.8 + 1.8 + p2 + 0.01, 0.0, 0.0, 1.8)]
check(len(surface.pieces(near)) == 1, "two atoms the probe cannot pass between are one piece")
check(len(surface.pieces(far)) == 2, "two atoms the probe can pass between are two pieces")
# A chain whose neighbours in the file are apart, joined only across the grid.
zig = [(i * 2.0, (i % 2) * 30.0, 0.0, 1.7) for i in range(20)]
zig_parts = surface.pieces(zig)
check(len(zig_parts) == 2 and sorted(len(p) for p in zig_parts) == [10, 10],
      "atoms are joined across the grid, not only to their neighbours in the file "
      "(%s)" % [len(p) for p in zig_parts])

# The test structure that shows the bug: MDM2 and a copy of it 45 A away, as
# two chains that do not touch.
apart = os.path.join(TMP, "two_apart.pdb")
with open(apart, "w") as fh:
    for L in open(REC):
        if L.startswith("ATOM") and L[21] == "A":
            fh.write(L)
    fh.write("TER\n")
    for L in open(REC):
        if L.startswith("ATOM") and L[21] == "A":
            fh.write(L[:21] + "C" + L[22:30] + "%8.3f" % (float(L[30:38]) + 45.0)
                     + L[38:])
    fh.write("END\n")
apart_atoms = structure.read(apart)
rows_ap, _k, _w = surface.sphere_rows(apart_atoms)
parts_ap = surface.pieces(rows_ap)
check(len(parts_ap) == 2 and [len(p) for p in parts_ap] == [705, 705],
      "two chains apart are two pieces of 705 atoms (%s)" % [len(p) for p in parts_ap])
check(len(surface.pieces(surface.sphere_rows(atoms)[0])) == 1,
      "MDM2 alone is one piece")

# --- the cache is versioned ------------------------------------------------------
# A surface cached before every piece had one is not served again.
check(surface.CACHE_VERSION >= 4, "the cache version was raised with the pieces (%s)"
      % surface.CACHE_VERSION)
_was_v = surface.CACHE_VERSION
surface.CACHE_VERSION = 3
old_key = surface.cache_key(apart_atoms, "/no/msms")
surface.CACHE_VERSION = _was_v
check(old_key != surface.cache_key(apart_atoms, "/no/msms"),
      "a surface cached under the old version has another key")
stale = os.path.join(TMP, "stale_cache")
os.makedirs(stale)
json.dump({"source": "msms", "vertices": 1, "stale": True},
          open(os.path.join(stale, old_key + ".json"), "w"))
try:
    surface.cached(apart_atoms, "/no/msms", stale)
    check(False, "a missing msms gave a surface")
except surface.SurfaceError:
    check(True, "a stale surface is not served")

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

# --- every piece has its outer surface, and no cavity is drawn --------------------


def mesh_parts(answer):
    """The answer's mesh split into its connected parts: [(vertex ids, signed
    volume)]. The outer surface of a piece encloses its atoms, so its volume
    is positive; a cavity's surface faces into the cavity, so its is
    negative."""
    nv_, nf_ = answer["vertices"], answer["triangles"]
    vv = struct.unpack("<%df" % (3 * nv_), base64.b64decode(answer["verts"]))
    ff = struct.unpack("<%dI" % (3 * nf_), base64.b64decode(answer["faces"]))
    parent = list(range(nv_))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for t in range(0, len(ff), 3):
        a, b, c = find(ff[t]), find(ff[t + 1]), find(ff[t + 2])
        parent[a] = b
        parent[find(c)] = b
    vol, members = {}, {}
    for t in range(0, len(ff), 3):
        i, j, k = ff[t], ff[t + 1], ff[t + 2]
        ax, ay, az = vv[3 * i:3 * i + 3]
        bx, by, bz = vv[3 * j:3 * j + 3]
        cx, cy, cz = vv[3 * k:3 * k + 3]
        r = find(i)
        vol[r] = vol.get(r, 0.0) + (ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx)
                                    + az * (bx * cy - by * cx)) / 6.0
    for i in range(nv_):
        members.setdefault(find(i), []).append(i)
    return [(members[r], vol.get(r, 0.0)) for r in members], vv


sep = surface.compute(apart_atoms, exe)
one_a = surface.compute([a for a in apart_atoms if a.chain == "A"], exe)
check(sep["pieces"] == 2 and sep["vertices"] == 2 * one_a["vertices"],
      "two chains apart each have their surface, the same as one alone "
      "(%s pieces, %s vertices, %s for one)"
      % (sep["pieces"], sep["vertices"], one_a["vertices"]))
# What the old single run gave: the outer surface of one piece only.
xyzr_ap, _k, _w = surface.spheres(apart_atoms, table)
old_run = surface.run_msms(exe, xyzr_ap, tempfile.mkdtemp(dir=TMP), density=3.0)
check(len(old_run[0]) // 3 == one_a["vertices"],
      "one msms run over both chains surfaces only one of them, which was the bug "
      "(%d vertices)" % (len(old_run[0]) // 3))
parts_sep, vv_sep = mesh_parts(sep)
res_sep = struct.unpack("<%dI" % sep["vertices"], base64.b64decode(sep["res"]))
by_chain = {}
for ids, _v in parts_sep:
    by_chain.setdefault(frozenset(sep["residues"][res_sep[i]][0] for i in ids), []).append(len(ids))
check(sorted("".join(sorted(k)) for k in by_chain) == ["A", "C"],
      "each part of the mesh lies over one chain, and both chains have one (%s)"
      % {"".join(sorted(k)): v for k, v in by_chain.items()})
xs_sep = vv_sep[0::3]
check(min(xs_sep) < min(xs) and max(xs_sep) > max(xs) + 40,
      "the surface reaches from the first chain to the copy 45 A away")
check(all(v > 0 for _ids, v in parts_sep),
      "every part is an outer surface (%s)" % [round(v) for _i, v in parts_sep])

# 1YCR, MDM2 with the p53 peptide, holds a sealed cavity MSMS finds with
# -all_components; it is not drawn.
whole = structure.read(REC)
both_ans = surface.compute(whole, exe)
parts_both, _vv = mesh_parts(both_ans)
check(both_ans["pieces"] == 1 and len(parts_both) == 1 and parts_both[0][1] > 0,
      "1YCR has one outer surface and no cavity (%d parts, volumes %s)"
      % (len(parts_both), [round(v) for _i, v in parts_both]))

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
