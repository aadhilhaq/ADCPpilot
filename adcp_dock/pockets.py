"""AutoSite's pockets, found before anything is submitted, and docked into as chosen.

AGFR runs AutoSite inside every box it builds maps over, and with a box set by
residues or a ligand it merges every pocket it finds there into the points
the search starts from (agfr's -p all, its default for any box but the whole
receptor). Over the whole receptor it takes the best scoring pocket. Either
way the pockets were only seen after the run, in its target file.

Here AutoSite runs while the site is being set: on a quick receptor PDBQT
(the chains kept, waters out, hydrogens added by prepare_receptor), inside
the box, or over the whole receptor for a box built around a pocket. That
takes seconds. The pockets come back ranked, with their points, and the ones
chosen are what the run docks into: their points replace AutoSite's own in
the target file once agfr has built it (swap_points), so the pocket seen is
the pocket searched, even though the run's receptor is prepared more
thoroughly than this one.

agfr itself has a "-p user" mode that reads points from a .npy file, but its
option parser only accepts the mode names, so the file can never be given on
its command line. Replacing translationPoints.npy in the target is the same
thing done after the fact, and ADCP reads that file and nothing else for its
starting points.
"""
import csv
import glob
import hashlib
import io
import json
import os
import shutil
import subprocess
import time
import zipfile

from . import adcp11

#: The most AutoSite is given. Measured on a 98 residue receptor: 98 s inside
#: a 26 A box, 13 s over the whole receptor. A box is the slower of the two,
#: not the quicker, because the grid is finer inside it, so a large receptor
#: can reach this limit; find() says so rather than failing quietly.
TIMEOUT = 240
#: Its address space, capped as the surfaces' msms is, since it runs beside
#: the page's server.
MEMORY = 4 << 30
#: The OpenBLAS threads prepare_receptor and AutoSite are given. The suite's
#: OpenBLAS starts one per core and reserves memory for each, so under MEMORY
#: prepare_receptor wrote no receptor from 32 cores up (measured on a cluster
#: compute node: 16 threads work, 32 and 48 fail; also reported on a 128-core
#: workstation).
#: One thread finds the same pockets in the same time as four.
BLAS_THREADS = 4
#: Points drawn per pocket in the page's view; the run is given all of them.
DRAW = 1500
#: The residues AutoSite's quick receptor leaves out.
WATER = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"}


class PocketError(Exception):
    pass


def autosite_script(exe):
    """AutoSite's AS.py inside the ADCP 1.1 environment `exe` belongs to."""
    if not exe:
        return None
    found = glob.glob(os.path.join(adcp11.home_of(exe), "lib", "python3*",
                                   "site-packages", "AutoSite", "bin", "AS.py"))
    return found[0] if found else None


def _limits():
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_AS, (MEMORY, MEMORY))
    except (ImportError, ValueError, OSError):
        pass
    try:
        os.nice(10)
    except OSError:
        pass


def receptor_text(path, chains=(), drop=(), keep_het=False):
    """The structure as AutoSite is given it: the chains kept (all when none
    are named) less any in `drop` (a ligand chain the box is built around is
    not part of the receptor), protein atoms, hetero groups only when the
    receptor keeps them and never waters, and one conformer."""
    keep, drop = set(chains or ()), set(drop or ())
    out = []
    for line in open(path, errors="replace"):
        rec = line[:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        ch = line[21:22]
        if (keep and ch not in keep) or ch in drop:
            continue
        if line[16:17] not in (" ", "A"):
            continue
        if rec == "HETATM" and (not keep_het or line[17:20].strip() in WATER):
            continue
        out.append(line.rstrip("\n")[:66].ljust(66) + "\n")
    if not out:
        raise PocketError("nothing is left of the receptor once its chains "
                          "are chosen")
    return "".join(out) + "END\n"


def key_of(**parts):
    """The cache key: the structure as it is on disk, and every choice."""
    src = parts.get("path") or ""
    try:
        st = os.stat(src)
        parts["stamp"] = [st.st_size, int(st.st_mtime)]
    except OSError:
        pass
    blob = json.dumps(parts, sort_keys=True, default=str).encode()
    return hashlib.sha1(blob).hexdigest()[:16]


def _points(pdb):
    pts = []
    for line in open(pdb, errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            try:
                pts.append([float(line[30:38]), float(line[38:46]),
                            float(line[46:54])])
            except ValueError:
                continue
    return pts


def _thin(pts, most=DRAW):
    if len(pts) <= most:
        return pts
    step = len(pts) / float(most)
    return [pts[int(i * step)] for i in range(most)]


def _summary(csvfile):
    """AutoSite's ranked table: rank, energy, points, radius of gyration,
    energy per point, buriedness and score, in the order it ranks them."""
    rows = []
    with open(csvfile, newline="") as fh:
        for r in list(csv.reader(fh))[1:]:
            if len(r) < 7:
                continue
            try:
                rows.append({"rank": int(float(r[0])), "energy": round(float(r[1]), 2),
                             "points": int(float(r[2])), "radius": round(float(r[3]), 2),
                             "energy_per_point": round(float(r[4]), 3),
                             "buriedness": round(float(r[5]), 3),
                             "score": round(float(r[6]), 2)})
            except ValueError:
                continue
    return rows


def find(path, cache, exe, chains=(), drop=(), keep_het=False, centre=None,
         size=None, peptide=True, log=None):
    """The pockets AutoSite finds, ranked, each with its points.

    `centre` and `size` confine it to a box; without them it searches the
    whole receptor. `peptide` ranks them with AutoSite's peptide score
    (E * buriedness^1.5) rather than its general one. Cached under `cache`
    by everything that changes the answer, so asking again is instant.
    """
    script = autosite_script(exe)
    if not script:
        from . import mac10
        if mac10.engaged():
            raise PocketError("finding pockets first needs AutoSite 1.1, from "
                              "the ADCP 1.1 suite, which Scripps does not "
                              "build for macOS; on a Mac, AutoSite 1.0 places "
                              "the points inside the box as the maps are "
                              "built")
        raise PocketError("AutoSite 1.1 is not installed here; it comes with "
                          "the ADCP 1.1 suite")
    if not path or not os.path.isfile(path):
        raise PocketError("choose a receptor first")
    key = key_of(path=os.path.abspath(path), chains=sorted(chains or ()),
                 drop=sorted(drop or ()), het=bool(keep_het),
                 centre=[round(c, 2) for c in centre] if centre else None,
                 size=[round(s, 1) for s in size] if size else None,
                 peptide=bool(peptide))
    folder = os.path.join(cache, key)
    done = os.path.join(folder, "pockets.json")
    if os.path.isfile(done):
        try:
            return json.load(open(done))
        except ValueError:
            pass
    shutil.rmtree(folder, ignore_errors=True)
    os.makedirs(folder)
    home = adcp11.home_of(exe)
    env = adcp11.agfr_env(exe)
    # Here rather than in agfr_env: the targets and the docking run without
    # the memory cap, and are not to be held to a few threads for it.
    env["OPENBLAS_NUM_THREADS"] = str(BLAS_THREADS)
    t0 = time.time()
    if path.lower().endswith(".pdbqt"):
        shutil.copyfile(path, os.path.join(folder, "rec.pdbqt"))
    else:
        with open(os.path.join(folder, "rec.pdb"), "w") as fh:
            fh.write(receptor_text(path, chains, drop, keep_het))
        prep = [os.path.join(home, "bin", "prepare_receptor"), "-r", "rec.pdb",
                "-o", "rec.pdbqt", "-A", "hydrogens"]
        _run(prep, folder, env, "prepare_receptor")
        if not os.path.isfile(os.path.join(folder, "rec.pdbqt")):
            raise PocketError("prepare_receptor wrote no receptor; its output "
                              "is in %s" % os.path.join(folder, "prepare_receptor.log"))
    os.makedirs(os.path.join(folder, "as"))
    cmd = [os.path.join(home, "bin", "python"), script, "-r", "rec.pdbqt",
           "-o", "as"]
    if centre and size:
        # Written with = so that a center starting with a minus sign is not
        # read as an option; AutoSite evaluates the text as a tuple.
        cmd += ["--boxcenter=(%.3f,%.3f,%.3f)" % tuple(centre),
                "--boxdim=(%.1f,%.1f,%.1f)" % tuple(size)]
    if peptide:
        cmd += ["-p"]
    _run(cmd, folder, env, "autosite")
    summ = glob.glob(os.path.join(folder, "as", "*_summary.csv"))
    rows = _summary(summ[0]) if summ else []
    pockets = []
    for r in rows:
        f = glob.glob(os.path.join(folder, "as", "*_cl_%03d.pdb" % r["rank"]))
        pts = _points(f[0]) if f else []
        if not pts:
            continue
        lo = [min(p[i] for p in pts) for i in range(3)]
        hi = [max(p[i] for p in pts) for i in range(3)]
        r.update(centre=[round((lo[i] + hi[i]) / 2.0, 3) for i in range(3)],
                 span=[round(hi[i] - lo[i], 1) for i in range(3)],
                 xyz=pts)
        pockets.append(r)
    out = {"key": key, "pockets": pockets, "seconds": round(time.time() - t0, 1),
           "scored": "peptide" if peptide else "general",
           "box": {"centre": list(centre), "size": list(size)} if centre and size else None}
    with open(done, "w") as fh:
        json.dump(out, fh)
    if log:
        log("AutoSite found %d pocket(s) in %.1f s" % (len(pockets), out["seconds"]))
    return out


def _run(cmd, cwd, env, name):
    logf = os.path.join(cwd, "%s.log" % name)
    with open(logf, "w") as fh:
        fh.write(" ".join(cmd) + "\n")
        fh.flush()
        try:
            subprocess.run(cmd, cwd=cwd, env=env, stdout=fh, stderr=subprocess.STDOUT,
                           timeout=TIMEOUT, preexec_fn=_limits)
        except subprocess.TimeoutExpired:
            raise PocketError("%s ran for more than %d s and was stopped"
                              % (name, TIMEOUT))


def for_page(found):
    """What the page draws: every pocket with its points thinned."""
    out = dict(found)
    out["pockets"] = [dict(p, xyz=_thin(p.get("xyz") or [])) for p in found["pockets"]]
    return out


def chosen(cache, key, ranks):
    """The points of the pockets chosen, by rank, from a search already made."""
    done = os.path.join(cache, key, "pockets.json")
    if not os.path.isfile(done):
        raise PocketError("the pockets found earlier are no longer there; "
                          "find them again")
    found = json.load(open(done))
    want = set(int(r) for r in ranks)
    sel = [p for p in found["pockets"] if p["rank"] in want]
    if not sel:
        raise PocketError("none of the pockets chosen were found")
    pts = [x for p in sel for x in p["xyz"]]
    return sel, pts, found


def box_around(points, padding=8.0, minimum=20.0, maximum=30.0):
    """(center, size, clipped axes) around points, by the rule every box here
    follows: the span plus the padding on each side, held between the
    limits."""
    lo = [min(p[i] for p in points) for i in range(3)]
    hi = [max(p[i] for p in points) for i in range(3)]
    want = [hi[i] - lo[i] + 2 * padding for i in range(3)]
    size = [round(max(minimum, min(maximum, w)), 1) for w in want]
    centre = [round((lo[i] + hi[i]) / 2.0, 3) for i in range(3)]
    clipped = [("xyz"[i], round(want[i], 1)) for i in range(3)
               if want[i] > maximum + 0.05]
    return centre, size, clipped


def save_points(points, path):
    """The points as agfr writes them: an (N, 3) float64 array in a .npy."""
    import numpy
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    numpy.save(path, numpy.asarray(points, dtype="float64").reshape(-1, 3))
    return path


def swap_points(trg, npy):
    """The target with its translation points replaced by these: every other
    entry copied as it was. Returns how many points it now holds."""
    import numpy
    pts = numpy.load(npy)
    buf = io.BytesIO()
    numpy.save(buf, pts.astype("float64").reshape(-1, 3))
    tmp = trg + ".swap"
    found = False
    with zipfile.ZipFile(trg) as src, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as dst:
        for item in src.infolist():
            if item.filename.endswith("/translationPoints.npy") or \
                    item.filename == "translationPoints.npy":
                dst.writestr(item, buf.getvalue())
                found = True
            else:
                dst.writestr(item, src.read(item.filename))
    if not found:
        os.remove(tmp)
        raise PocketError("the target has no translationPoints.npy to replace")
    os.replace(tmp, trg)
    return int(pts.reshape(-1, 3).shape[0])
