"""Building the grid maps ADCP docks into, and verifying they are usable.

By Aadhil Haq

A .trg file is a zip archive: one grid map per atom type, plus the receptor
it was built from and the box it was built over. agfr writes it, and every
later stage trusts it.

Three failure modes produce no error of their own:

* a full disk truncates the archive, and a builder that only checks for the
  file's existence returns the truncated one on the next run;
* it is built under one name and used under another, so ADCP unpacks it and
  finds no maps; it warns once per atom type and then docks against an
  empty grid, returning scores;
* a stale copy in a working directory outlives the target it came from.

In each case ADCP still runs and still returns numbers, so the archive's
contents are checked before anything docks.

The ADCP 1.1 suite's agfr builds it, found from the same installation as
ADCP 1.1 (adcp11.find_agfr). It calls the suite's AutoGrid (4.2.7 of 2020),
and writes a constriction map as well, which ADCP does not read.

Two settings change what it builds, so both are spelled out on every command
here rather than left to the program's defaults: the AutoSite that places
the translation points (-asv 1.0, the program's default, or -asv 1.1), and
the receptor gradient, added only when asked (-g).
"""
import os
import re
import subprocess
import time
import zipfile

MIN_MAPS = 20

#: The AutoDock atom types a peptide needs a map for, by one-letter code.
#: Every residue carries the backbone set; the entries below add what its
#: side chain brings. agfr writes a map per type and ADCP loads one map for
#: each type its peptide uses, so a target built for another ligand, with
#: agfr's -m ligand or a hand-picked list in AGFRGUI, can be missing one.
#: ADCP then stops with "Missing gridmap file" after the job has queued,
#: which is what this check is here to prevent.
BACKBONE_TYPES = frozenset(("C", "N", "OA", "HD"))
SIDE_CHAIN_TYPES = {
    "A": frozenset(), "G": frozenset(), "P": frozenset(),
    "I": frozenset(), "L": frozenset(), "V": frozenset(),
    "R": frozenset(("N",)), "K": frozenset(("N",)),
    "N": frozenset(("N", "OA")), "Q": frozenset(("N", "OA")),
    "D": frozenset(("OA",)), "E": frozenset(("OA",)),
    "S": frozenset(("OA",)), "T": frozenset(("OA",)),
    "C": frozenset(("SA",)), "M": frozenset(("SA",)),
    "H": frozenset(("A", "N", "NA")),
    "F": frozenset(("A",)), "Y": frozenset(("A", "OA")),
    "W": frozenset(("A", "N")),
}


def types_needed(sequences):
    """The AutoDock types the maps must cover for these peptides, and the
    non-canonical side chains whose types cannot be told from the sequence.

    Returns (types, unknown). A side chain in angle brackets comes from
    ADCP's rotamer library and can carry a type no canonical residue does,
    so it is reported rather than guessed at.
    """
    types, unknown = set(), []
    for seq in sequences or ():
        rest = str(seq or "")
        while rest:
            m = re.match(r"&?([A-Za-z])(?:<([^>]*)>)?", rest)
            if not m:
                rest = rest[1:]
                continue
            code = m.group(1).upper()
            if m.group(2):
                unknown.append(m.group(2))
            types |= BACKBONE_TYPES
            types |= SIDE_CHAIN_TYPES.get(code, frozenset())
            rest = rest[m.end():]
    return types, sorted(set(unknown))


def types_in(names):
    """The atom types a target's map filenames carry."""
    got = set()
    for n in names:
        m = re.search(r"rigidReceptor\.([A-Za-z]{1,2})\.map$", n)
        if m:
            got.add(m.group(1))
    return got

#: Which agfr builds the maps: the ADCP 1.1 suite's. "auto", which new runs
#: said before, means it too.
AGFR_CHOICES = ("1.1", "auto")
#: Which AutoSite places the translation points and finds the pockets.
AUTOSITE_CHOICES = ("1.1", "1.0")
#: Where the box comes from.
SITE_MODES = ("residues", "ligand", "pocket")
#: Which pocket an AutoSite box is built around: AutoSite's best scoring
#: one (agfr -p best), or every pocket it finds merged into one box (-p all).
POCKET_CHOICES = ("best", "all")
#: The spacing the docking maps are built at, agfr's own default.
SPACING = 0.375


class TargetError(Exception):
    pass


def check(trg, receptor_pdbqt=None, sequences=None):
    """None if the target is usable, else a sentence saying what is wrong.

    With `sequences`, the maps are also checked against the atom types those
    peptides need, so a target that cannot dock them is refused here rather
    than by ADCP after the job has queued.
    """
    if not os.path.isfile(trg):
        return "does not exist"
    if not zipfile.is_zipfile(trg):
        return "not a zip archive (agfr did not finish, or the disk was full)"
    z = zipfile.ZipFile(trg)
    names = z.namelist()
    if not names:
        return "empty archive"
    maps = [n for n in names if n.endswith(".map")]
    if len(maps) < MIN_MAPS:
        return "only %d grid maps, expected at least %d" % (len(maps), MIN_MAPS)
    recs = [n for n in names if n.endswith("rigidReceptor.pdbqt")]
    if not recs:
        return "no receptor inside the archive"
    if receptor_pdbqt and os.path.isfile(receptor_pdbqt):
        got = sum(1 for L in z.read(recs[0]).decode("utf8", "replace").splitlines()
                  if L[:4] == "ATOM")
        want = sum(1 for L in open(receptor_pdbqt) if L[:4] == "ATOM")
        if got != want:
            return ("the receptor inside has %d atoms but %s has %d; this "
                    "target was built from a different structure"
                    % (got, os.path.basename(receptor_pdbqt), want))
    # ADCP unpacks to ./tmp_<job>/<target basename>/, so the folder
    # inside the archive must carry the target's own name or no map is
    # found.
    prefix = names[0].split("/")[0]
    stem = os.path.basename(trg)[:-4]
    if prefix != stem:
        return ("the folder inside the archive is %r but the target is named "
                "%r; ADCP looks for the second and finds nothing" % (prefix, stem))
    if sequences:
        want, unknown = types_needed(sequences)
        missing = sorted(want - types_in(names))
        if missing:
            return ("no grid map for the atom type(s) %s, which these peptides "
                    "need; this target was built for a different ligand. Build "
                    "the maps again over the same box with every atom type"
                    % ", ".join(missing))
    return None


def box_of(trg):
    """The box agfr recorded, as its own config states it."""
    try:
        z = zipfile.ZipFile(trg)
        for n in z.namelist():
            if n.endswith(".cfg"):
                for L in z.read(n).decode("utf8", "replace").splitlines():
                    if L.startswith("boxMode"):
                        return L.split("=", 1)[1].strip()
    except (OSError, zipfile.BadZipFile):
        pass
    return None


class Maps:
    """How a run's maps are built: by the ADCP 1.1 suite's AGFR, with which
    AutoSite, and whether the receptor gradient is added.

    `exe` is ADCP 1.1's executable, from which its AGFR is found.
    """

    agfr = "1.1"

    def __init__(self, exe=None, autosite="1.1", gradient=False):
        self.exe = exe
        self.autosite = str(autosite)
        self.gradient = bool(gradient)

    def runner(self, tools=None):
        """The program to start: the suite's own Python on its runAGFR.py."""
        from . import adcp11
        return adcp11.agfr_command(self.exe)

    def env(self):
        from . import adcp11
        return adcp11.agfr_env(self.exe)

    def flags(self):
        """Both settings, given explicitly whatever the program's defaults."""
        return ["-asv", self.autosite] + (["-g"] if self.gradient else [])

    def describe(self):
        return ("AGFR 1.1, AutoSite %s, receptor gradient %s"
                % (self.autosite, "on" if self.gradient else "off"))

    def record(self):
        return {"agfr": self.agfr, "autosite": self.autosite,
                "receptor_gradient": self.gradient}


def maps_for(j):
    """The Maps a job asks for.

    A job file from before the maps could be chosen has no agfr setting:
    its maps were built by an older AGFR with that program's own defaults,
    AutoSite 1.0 and the receptor gradient on, so they are built again with
    the same two, which comes closest to the maps it was docked into. A job
    file that records its settings is built with what it records.
    """
    from . import adcp11, mac10
    if mac10.engaged():
        # On a Mac, ADFRsuite 1.0's agfr with AutoSite 1.0, the gradient as
        # the run asks (mac10.maps).
        return mac10.maps(gradient=bool(j.get("receptor_gradient", True)))
    exe, script = adcp11.find_agfr(j.get("tools_roots") or [])
    if not (exe and script):
        raise TargetError("AGFR comes with the ADCP 1.1 suite, which is not "
                          "installed here")
    if "agfr" not in j:
        return Maps(exe=exe, autosite="1.0", gradient=True)
    return Maps(exe=exe, autosite=str(j.get("autosite") or "1.0"),
                gradient=bool(j.get("receptor_gradient", True)))


def command(maps, tools, receptor_pdbqt, name, centre=None, size=None,
            pocket=None, padding=None):
    """The agfr command line for a box given by its center and edges
    (-b user), or for a box around AutoSite's pocket (-b fill -p pocket)."""
    cmd = maps.runner(tools) + ["-r", os.path.abspath(receptor_pdbqt),
                                "-o", name]
    if pocket:
        cmd += ["-b", "fill", "-p", pocket]
        if padding is not None:
            cmd += ["-P", "%g" % float(padding)]
    else:
        cmd += ["-b", "user", "%.3f" % centre[0], "%.3f" % centre[1],
                "%.3f" % centre[2], "%.1f" % size[0], "%.1f" % size[1],
                "%.1f" % size[2]]
    return cmd + maps.flags()


def _run_agfr(cmd, outdir, name, env=None):
    """Run agfr where its output belongs. Returns (target or None, seconds)."""
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    with open(os.path.join(outdir, "%s_agfr.out" % name), "w") as fh:
        fh.write(" ".join(cmd) + "\n")
        fh.flush()
        subprocess.run(cmd, cwd=outdir, stdout=fh, stderr=subprocess.STDOUT,
                       env=env)
    trg = os.path.join(outdir, "%s.trg" % name)
    return (trg if os.path.isfile(trg) else None), round(time.time() - t0, 1)


def build(receptor_pdbqt, centre, size, outdir, tools, name="receptor",
          log=print, maps=None):
    """agfr over the box, then the check. Raises TargetError on failure.

    An existing .trg is returned untouched, which makes a resume cheap and a
    truncated target permanent. Validating on the way out covers both: a
    good target short-circuits, a bad one is deleted and rebuilt once.

    `maps` is a Maps, maps_for(job).
    """
    trg = os.path.join(outdir, "%s.trg" % name)
    why = check(trg, receptor_pdbqt) if os.path.isfile(trg) else "missing"
    attempts = 0
    while why and attempts < 2:
        if os.path.isfile(trg):
            log("  the grid target is unusable (%s); building it again" % why)
            try:
                os.remove(trg)
            except OSError:
                pass
        attempts += 1
        # agfr's output is kept beside the target.
        got, seconds = _run_agfr(
            command(maps, tools, receptor_pdbqt, name, centre, size),
            outdir, name, maps.env())
        if got is None:
            raise TargetError(
                "agfr built no target; its output is in %s"
                % os.path.join(outdir, "%s_agfr.out" % name))
        trg = got
        why = check(trg, receptor_pdbqt)
        if not why:
            log("  maps built in %.0f s" % seconds)
    if why:
        raise TargetError("%s: %s" % (os.path.basename(trg), why))
    log("  grid target %s: %d maps, box %s; %s"
        % (os.path.basename(trg), maps_in(trg), box_of(trg), maps.describe()))
    return trg


# AutoSite's pocket table as agfr prints it once it has chosen its pockets:
#     pocket|  energy | # of |Rad. of | energy |   bns    | score
#     number|         |points|gyration|per vol.|buriedness|v*b^2/rg
#          1   -175.24   513    6.44     -0.34      0.78      48.85
_POCKET_ROW = re.compile(r"^\s+(\d+)\s+(-?\d+\.\d+)\s+(\d+)\s+(-?\d+\.\d+)\s+"
                         r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$")


def pocket_table(text):
    """[{rank, energy, points, radius, energy_per_point, buriedness, score}]
    from agfr's output, the table it printed after "found N pocket(s)". The
    tables AutoSite 1.1 prints while it scans its cutoffs come before that
    line and are not read."""
    out, on = [], False
    for line in text.splitlines():
        if re.search(r"found \d+ pocket\(s\)", line):
            out, on = [], True
            continue
        if not on:
            continue
        m = _POCKET_ROW.match(line)
        if m:
            out.append({"rank": int(m.group(1)), "energy": float(m.group(2)),
                        "points": int(m.group(3)),
                        "radius": float(m.group(4)),
                        "energy_per_point": float(m.group(5)),
                        "buriedness": float(m.group(6)),
                        "score": float(m.group(7))})
        elif out and line.strip() and not line.strip().startswith("-"):
            on = False
    return out


def find_pocket(receptor_pdbqt, outdir, tools, maps, pocket="best",
                padding=8.0, log=print):
    """Where AutoSite puts the box: agfr -b fill, run for its box alone.

    agfr -b fill runs AutoSite over the whole receptor, takes the pocket
    asked for and sizes the box around its points with the padding on each
    side. It then builds the maps at the 1 A spacing AutoSite searched at
    (runAGFR.py sets the box with setBox(['receptor'], self.padding, 1.0)
    before AutoSite runs and setBoxForCoords(..., self.spacing) after, so the
    -s spacing is lost), and maps that coarse are not docked into here: the
    box is kept and the maps are built over it again at 0.375 A, as they are
    for every other kind of box.

    Returns {"centre", "size", "points", "pockets", "command", "seconds"}.
    """
    sub = os.path.join(outdir, "pocket")
    name = "pocket"
    # The gradient changes nothing AutoSite or the box sees, and costs time.
    # Built the same way as the run's own maps (on a Mac, mac10's).
    quick = type(maps)(exe=maps.exe, autosite=maps.autosite, gradient=False)
    cmd = command(quick, tools, receptor_pdbqt, name, pocket=pocket,
                  padding=padding)
    trg = os.path.join(sub, "%s.trg" % name)
    seconds = None
    if not (os.path.isfile(trg) and zipfile.is_zipfile(trg)):
        trg, seconds = _run_agfr(cmd, sub, name, quick.env())
        if trg is None:
            raise TargetError("AutoSite found no pocket to build a box around; "
                              "agfr's output is in %s"
                              % os.path.join(sub, "%s_agfr.out" % name))
    from . import agfr_target
    centre, size = agfr_target.read_box(trg)
    try:
        text = open(os.path.join(sub, "%s_agfr.out" % name),
                    errors="replace").read()
    except OSError:
        text = ""
    pts = translation_points(trg)
    table = pocket_table(text)
    log("  AutoSite %s found %d pocket(s); the box is around %s: center "
        "%.2f %.2f %.2f, %.1f x %.1f x %.1f A"
        % (maps.autosite, len(table),
           "the best scoring one" if pocket == "best" else "all of them",
           centre[0], centre[1], centre[2], size[0], size[1], size[2]))
    return {"centre": [float(c) for c in centre],
            "size": [float(s) for s in size], "points": pts,
            "pockets": table, "command": cmd, "seconds": seconds}


def translation_points(trg):
    """The AutoSite points inside a target, as [[x, y, z]], or [].

    Read from translationPoints.npy without unpickling anything: a target
    file can come from anywhere."""
    import io
    try:
        import numpy as np
        z = zipfile.ZipFile(trg)
        name = next((n for n in z.namelist()
                     if n.endswith("translationPoints.npy")), None)
        if not name:
            return []
        arr = np.load(io.BytesIO(z.read(name)), allow_pickle=False)
        return [[round(float(v), 3) for v in p] for p in arr.reshape(-1, 3)]
    except Exception:                                           # noqa: BLE001
        return []


def maps_in(trg):
    try:
        return sum(1 for n in zipfile.ZipFile(trg).namelist()
                   if n.endswith(".map"))
    except (OSError, zipfile.BadZipFile):
        return 0
