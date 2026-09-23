"""Building the grid maps ADCP docks into, and verifying they are usable.

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
"""
import os
import zipfile

from ._vendor import adcp

MIN_MAPS = 20


class TargetError(Exception):
    pass


def check(trg, receptor_pdbqt=None):
    """None if the target is usable, else a sentence saying what is wrong."""
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


def build(receptor_pdbqt, centre, size, outdir, tools, name="receptor",
          log=print):
    """The vendored builder, then the check. Raises TargetError on failure.

    build_target returns an existing .trg untouched, which makes a resume cheap
    and a truncated target permanent. Validating on the way out covers both: a
    good target short-circuits, a bad one is deleted and rebuilt once.
    """
    trg = adcp.build_target(receptor_pdbqt, centre, size, outdir, tools,
                            name=name)
    why = check(trg, receptor_pdbqt)
    if why:
        log("  the grid target is unusable (%s); building it again" % why)
        try:
            os.remove(trg)
        except OSError:
            pass
        trg = adcp.build_target(receptor_pdbqt, centre, size, outdir, tools,
                                name=name)
        why = check(trg, receptor_pdbqt)
        if why:
            raise TargetError("%s: %s" % (os.path.basename(trg), why))
    log("  grid target %s: %d maps, box %s"
        % (os.path.basename(trg), maps_in(trg), box_of(trg)))
    return trg


def maps_in(trg):
    try:
        return sum(1 for n in zipfile.ZipFile(trg).namelist()
                   if n.endswith(".map"))
    except (OSError, zipfile.BadZipFile):
        return 0
