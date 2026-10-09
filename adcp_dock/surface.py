"""Molecular surfaces for the viewers, computed by MSMS.

By Aadhil Haq

MSMS, by Michel Sanner, triangulates the solvent excluded surface: the
surface a probe sphere of 1.5 Angstrom traces as it rolls over the atoms.
The msms program ships in ADFRsuite's bin folder. It is called here as the
installed program, never copied, because its licence allows non-commercial
scientific research only (THIRD_PARTY_NOTICES.md).

MSMS's own inputs are used. Each atom is a sphere with the radius
pdb_to_xyzr, MSMS's own converter, would give it. ADFRsuite carries that
table as MolKit's radii_patterns.py, "array of radii taken from
pdb_to_xyzr", and it is read from the installed tree as data rather than
imported or copied. Hydrogens are left out and the heavy atoms take the
table's united atom radii, which is what pdb_to_xyzr does by default, so a
structure with hydrogens and the same structure without them give the same
surface. Where the table cannot be read, an element table of this module's
own stands in, and the answer says so.

The answer is the triangulation itself: vertices, normals and triangles,
packed as little-endian binary and base64 encoded, so a surface of fifty
thousand vertices costs the page about a megabyte rather than several. Each
vertex also carries the residue of the atom it lies over, so a viewer can
colour a patch of the surface by residue.

Where msms is absent, refuses, runs out of time or memory, or the structure
is larger than is sensible to triangulate on a shared login node, the answer
says why and the viewer draws its own surface instead.
"""
import ast
import base64
import hashlib
import json
import os
import re
import shutil
import struct
import subprocess
import tempfile
import threading
import time

#: The probe, in Angstrom: a water molecule, and MSMS's own default.
PROBE = 1.5

#: Above this many heavy atoms no surface is computed and the viewer draws
#: its own. MSMS on GroEL with GroES (1AON, 58,674 atoms) at the density
#: below took a few seconds and well under the memory limit on a compute
#: node; this leaves room above that for a slower login node.
MAX_ATOMS = 120000

#: Wall seconds and address space MSMS may use. It runs on the machine that
#: serves the page, which is usually a shared login node, where long or
#: large processes get the account flagged. GroEL with GroES peaked at about
#: 250 MB and 5 s, so these leave room without letting one surface hog it.
TIMEOUT = 30
MEMORY = 1 << 30

#: Each answer is cached by what it was computed from, so a structure is
#: triangulated once whoever asks for it next.
CACHE_VERSION = 3

#: One MSMS at a time, so two people opening large structures at once do
#: not both hold their memory on the same login node.
_RUN_LOCK = threading.Lock()


class SurfaceError(Exception):
    """Why no surface could be computed. The page shows it, quietly, and
    draws its own surface instead."""


def density_for(n_atoms):
    """Vertices per square Angstrom, by size: fine for a small receptor,
    coarser as it grows, so a small one is smooth and a large complex stays
    a size a browser draws without trouble. MSMS's own default is 1.0.

    Measured with ADFRsuite's msms 2.6.1 on a compute node: MDM2 (1YCR
    chain A, 705 heavy atoms) at 3.0 gives 11,488 vertices in 0.15 s; the
    dengue envelope dimer (7A3R, 5,912) at 1.0 gives 38,530 in 0.4 s; GroEL
    with GroES (1AON, 58,870) at 0.5 gives 234,229 in 4.5 s.
    """
    n = int(n_atoms or 0)
    if n <= 2000:
        return 3.0
    if n <= 5000:
        return 2.0
    if n <= 15000:
        return 1.0
    if n <= 40000:
        return 0.5
    return 0.3


# --- radii -----------------------------------------------------------------

#: Element radii for a heavy atom when ADFRsuite's table cannot be read, in
#: Angstrom, as united atoms: a carbon or nitrogen counts the hydrogens it
#: carries, as pdb_to_xyzr's own united radii do.
ELEMENT_RADII = {"C": 1.90, "N": 1.70, "O": 1.50, "S": 1.80, "P": 1.80,
                 "SE": 2.00, "ZN": 1.40, "MG": 1.60, "CA": 1.97, "FE": 1.30,
                 "MN": 1.24, "CU": 1.40, "CO": 1.25, "NI": 1.25, "CD": 1.49,
                 "NA": 1.80, "K": 2.20, "CL": 1.75, "BR": 1.85, "I": 1.98,
                 "F": 1.47}
DEFAULT_RADIUS = 1.80

_TABLES = {}


def adfr_home(msms_exe):
    """ADFRsuite's folder, from the msms program inside its bin folder."""
    if not msms_exe:
        return None
    return os.path.dirname(os.path.dirname(os.path.realpath(msms_exe)))


def radii_file(home):
    """MolKit's copy of pdb_to_xyzr's table in an ADFRsuite folder, or None."""
    if not home:
        return None
    for sub in ("CCSBpckgs", os.path.join("lib", "python2.7", "site-packages")):
        p = os.path.join(home, sub, "MolKit", "radii_patterns.py")
        if os.path.isfile(p):
            return p
    return None


def read_radii_table(path):
    """(patterns, radii) from radii_patterns.py, read as data.

    patterns is [(residue regex, atom regex, type number)] in the table's
    order, the first match winning, as in pdb_to_xyzr; radii is {type
    number: (covalent, explicit, united)}. Only the two literal assignments
    are evaluated, with ast.literal_eval, so nothing in the file is run.
    """
    if path in _TABLES:
        return _TABLES[path]
    import warnings
    with warnings.catch_warnings():
        # The file writes regular expressions in plain strings ("\\*"),
        # which a newer Python warns about while reading them the same way.
        warnings.simplefilter("ignore")
        tree = ast.parse(open(path).read(), path)
    found = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 \
                and isinstance(node.targets[0], ast.Name) \
                and node.targets[0].id in ("AtomRadiiPatterns", "AAradii"):
            found[node.targets[0].id] = ast.literal_eval(node.value)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pats = [(re.compile(r), re.compile(a), int(n))
                for r, a, n in found["AtomRadiiPatterns"]]
    radii = {int(k): (float(v[0]), float(v[1]), float(v[2]))
             for k, v in found["AAradii"].items()}
    _TABLES[path] = (pats, radii)
    return pats, radii


def radius_of(resname, name, element, table=None):
    """An atom's united atom radius: from pdb_to_xyzr's table when there is
    one and it names the atom, from the element otherwise."""
    if table is not None:
        pats, radii = table
        res = (resname or "XXX").strip() or "XXX"
        nm = (name or "").strip()
        for rp, ap, num in pats:
            if ap.search(nm) and rp.search(res):
                r = radii.get(num)
                if r:
                    return r[2]
                break
    el = (element or "").strip().upper() or (name or "").strip()[:1].upper()
    return ELEMENT_RADII.get(el, DEFAULT_RADIUS)


def _is_h(a):
    e = (getattr(a, "element", "") or "").strip().upper()
    if e:
        return e in ("H", "D")
    n = (getattr(a, "name", "") or "").strip().lstrip("0123456789")
    return n[:1] in ("H", "D")


def spheres(atoms, table=None):
    """(xyzr text, residue keys, index of each sphere's residue) for the
    heavy atoms. Hydrogens are left out, as pdb_to_xyzr leaves them out by
    default: the united atom radii already count them."""
    lines, keys, which, seen = [], [], [], {}
    for a in atoms:
        if _is_h(a):
            continue
        r = radius_of(a.resname, a.name, a.element, table)
        lines.append("%.3f %.3f %.3f %.2f" % (a.x, a.y, a.z, r))
        k = "%s:%s%s" % ((a.chain or "").strip() or "_", a.resnum,
                         (a.icode or "").strip())
        if k not in seen:
            seen[k] = len(keys)
            keys.append(k)
        which.append(seen[k])
    return "\n".join(lines) + "\n", keys, which


# --- msms itself -------------------------------------------------------------

def version(exe):
    """The first line msms prints about itself, or ''. It names its
    version and build."""
    if not exe:
        return ""
    try:
        p = subprocess.run([exe], capture_output=True, text=True, timeout=20,
                           stdin=subprocess.DEVNULL)
    except (OSError, subprocess.SubprocessError):
        return ""
    for line in (p.stdout + "\n" + p.stderr).splitlines():
        if "MSMS" in line and re.search(r"\d", line):
            return line.strip()
    return ""


def _limits():
    """Run in the child before msms starts: cap its address space."""
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_AS, (MEMORY, MEMORY))
    except (ImportError, ValueError, OSError):
        pass


def run_msms(exe, xyzr, workdir, probe=PROBE, density=1.0, timeout=TIMEOUT):
    """(vertices, normals, closest sphere per vertex, faces) from msms, or
    SurfaceError. Indices are 0-based."""
    inp = os.path.join(workdir, "in.xyzr")
    out = os.path.join(workdir, "surface")
    with open(inp, "w") as fh:
        fh.write(xyzr)
    cmd = [exe, "-if", inp, "-of", out, "-probe_radius", "%g" % probe,
           "-density", "%g" % density, "-no_area"]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout,
                           stdin=subprocess.DEVNULL, cwd=workdir,
                           preexec_fn=_limits if os.name == "posix" else None)
    except subprocess.TimeoutExpired:
        raise SurfaceError("MSMS took longer than %d s" % timeout)
    except OSError as exc:
        raise SurfaceError("msms could not be started (%s)" % exc)
    vert, face = out + ".vert", out + ".face"
    if not (os.path.isfile(vert) and os.path.isfile(face)):
        said = [l.strip() for l in (p.stdout + "\n" + p.stderr).splitlines()
                if l.strip()]
        why = said[-1] if said else "no output"
        if p.returncode and p.returncode < 0:
            why = "stopped by signal %d; %s" % (-p.returncode, why)
        raise SurfaceError("MSMS wrote no surface (%s)" % why[:160])
    v, n, atom = [], [], []
    for line in open(vert):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 8:
            continue                       # the counts line
        v.extend((float(f[0]), float(f[1]), float(f[2])))
        n.extend((float(f[3]), float(f[4]), float(f[5])))
        atom.append(max(0, int(f[7]) - 1))
    faces = []
    for line in open(face):
        if line.startswith("#"):
            continue
        f = line.split()
        if len(f) < 5:
            continue
        # In MSMS's own order: 3Dmol lights a triangle from the side its
        # winding faces, and MSMS winds them facing out, the way its normals
        # point.
        faces.extend((int(f[0]) - 1, int(f[1]) - 1, int(f[2]) - 1))
    if not v or not faces:
        raise SurfaceError("MSMS wrote an empty surface")
    return v, n, atom, faces


def _b64(fmt, values):
    return base64.b64encode(struct.pack("<%d%s" % (len(values), fmt),
                                        *values)).decode("ascii")


def compute(atoms, exe, probe=PROBE, density=None, timeout=TIMEOUT,
            workdir=None):
    """The surface of these atoms as a dict the page can draw.

    Raises SurfaceError with a reason a person can read.
    """
    if not exe or not os.path.isfile(exe):
        raise SurfaceError("MSMS is not installed here (ADFRsuite's bin/msms "
                           "was not found)")
    table_path = radii_file(adfr_home(exe))
    table = None
    if table_path:
        try:
            table = read_radii_table(table_path)
        except Exception:                                      # noqa: BLE001
            table = None
    xyzr, keys, which = spheres(atoms, table)
    n = len(which)
    if not n:
        raise SurfaceError("there are no heavy atoms to put a surface on")
    if n > MAX_ATOMS:
        raise SurfaceError("%d atoms is more than this server triangulates "
                           "(%d)" % (n, MAX_ATOMS))
    dens = float(density or density_for(n))
    tmp = tempfile.mkdtemp(prefix="msms_", dir=workdir)
    t0 = time.time()
    try:
        with _RUN_LOCK:
            v, nrm, atom, faces = run_msms(exe, xyzr, tmp, probe, dens,
                                           timeout)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    secs = time.time() - t0
    res = [which[a] if a < len(which) else 0 for a in atom]
    return {"source": "msms", "probe": probe, "density": dens,
            "atoms": n, "vertices": len(v) // 3, "triangles": len(faces) // 3,
            "seconds": round(secs, 2),
            "radii": ("pdb_to_xyzr united atom radii, read from %s"
                      % os.path.relpath(table_path, adfr_home(exe))
                      if table else "element radii (pdb_to_xyzr's table "
                      "was not found)"),
            "encoding": "base64 little-endian: verts and norms float32 x,y,z; "
                        "faces uint32 0-based; res uint32 index into residues",
            "verts": _b64("f", v), "norms": _b64("f", nrm),
            "faces": _b64("I", faces), "res": _b64("I", res),
            "residues": keys}


def cache_key(atoms, exe, probe=PROBE, density=None):
    """What a cached surface was computed from: the spheres, the probe, the
    density and the msms program."""
    h = hashlib.sha1()
    for a in atoms:
        if _is_h(a):
            continue
        h.update(("%s %s %s %s %.3f %.3f %.3f\n"
                  % (a.chain, a.resnum, a.resname, a.name, a.x, a.y, a.z))
                 .encode())
    try:
        st = os.stat(exe)
        h.update(("%s %d %d" % (exe, st.st_size, int(st.st_mtime))).encode())
    except (OSError, TypeError):
        h.update(b"no msms")
    h.update(("%g %s %d" % (probe, density, CACHE_VERSION)).encode())
    return h.hexdigest()


def cached(atoms, exe, cache_dir, probe=PROBE, density=None, log=None):
    """compute(), read from cache_dir when it was computed before, and
    written there when it was not. The cache is an optimisation: a folder
    that cannot be written is not an error."""
    key = cache_key(atoms, exe, probe, density)
    path = os.path.join(cache_dir, key + ".json") if cache_dir else None
    if path and os.path.isfile(path):
        try:
            got = json.load(open(path))
            got["cached"] = True
            return got
        except (OSError, ValueError):
            pass
    got = compute(atoms, exe, probe, density)
    if path:
        try:
            os.makedirs(cache_dir, exist_ok=True)
            tmp = path + ".part"
            with open(tmp, "w") as fh:
                json.dump(got, fh)
            os.replace(tmp, path)
        except OSError:
            pass
    got["cached"] = False
    return got
