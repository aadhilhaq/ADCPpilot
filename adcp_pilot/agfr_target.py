"""A run that starts from a target file made elsewhere, such as in AGFRGUI.

By Aadhil Haq

AGFR (and its interface, AGFRGUI) writes a .trg file: a zip archive holding
the grid maps, the receptor they were computed from, and the grid parameters.
A run given one skips this package's receptor preparation and map building.
The receptor for the viewer and MM-GBSA, and the box, are read from the file.
"""
import os
import re
import shutil
import zipfile

from . import target

# AutoDock atom types to elements, for writing the receptor as a PDB.
_AD2EL = {"A": "C", "C": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O",
          "OS": "O", "SA": "S", "S": "S", "HD": "H", "HS": "H", "H": "H",
          "P": "P", "F": "F", "Cl": "Cl", "CL": "Cl", "Br": "Br", "BR": "Br",
          "I": "I", "Mg": "Mg", "MG": "Mg", "Ca": "Ca", "CA": "Ca", "Mn": "Mn",
          "MN": "Mn", "Fe": "Fe", "FE": "Fe", "Zn": "Zn", "ZN": "Zn"}


#: The suffix AGFR and AGFRGUI give a target file. The page's file browser
#: lists these, and only these, when it is choosing a target.
SUFFIX = ".trg"


def is_target(path):
    """Is this a target file, by name? The contents are target.check's
    business."""
    return str(path).lower().endswith(SUFFIX)


def preview(src):
    """(receptor as PDB text, center, size) from a target file, for the page.

    Nothing is copied or unpacked to disk. The file is held to the same check
    import_target applies, so a target the page shows is one the run will
    take, and one the run would refuse is refused here for the same reason.
    Raises ValueError carrying that reason.
    """
    why = target.check(src)
    if why:
        raise ValueError(why)
    z = zipfile.ZipFile(src)
    rec = next(n for n in z.namelist() if n.endswith("rigidReceptor.pdbqt"))
    text = z.read(rec).decode("utf8", "replace")
    centre, size = read_box(src)
    return pdbqt_to_pdb(text), centre, size


def read_box(trg):
    """(center, size) in Angstrom from the grid parameter file in the target."""
    z = zipfile.ZipFile(trg)
    gpf = next((n for n in z.namelist() if n.endswith(".gpf")), None)
    if not gpf:
        raise ValueError("no grid parameter file (.gpf) inside %s" % trg)
    centre = npts = spacing = None
    for line in z.read(gpf).decode("utf8", "replace").splitlines():
        parts = line.split("#")[0].split()
        if not parts:
            continue
        if parts[0] == "gridcenter" and len(parts) >= 4:
            centre = [float(x) for x in parts[1:4]]
        elif parts[0] == "npts" and len(parts) >= 4:
            npts = [int(x) for x in parts[1:4]]
        elif parts[0] == "spacing" and len(parts) >= 2:
            spacing = float(parts[1])
    if centre is None or npts is None or spacing is None:
        raise ValueError("the grid parameter file in %s has no box" % trg)
    return centre, [n * spacing for n in npts]


def pdbqt_to_pdb(text):
    """The receptor from a PDBQT as a PDB: coordinates kept, charges and
    AutoDock types replaced by element symbols."""
    out = []
    for line in text.splitlines():
        if line[:6] in ("ATOM  ", "HETATM"):
            # The AutoDock type is the last field after the charge; read by
            # field rather than column, since writers differ by a column.
            tail = line[66:].split()
            ad = tail[-1] if len(tail) >= 2 else ""
            el = _AD2EL.get(ad, (ad or line[12:16].strip()[:1])[:2])
            out.append("%-66s          %2s" % (line[:66], el.upper() if len(el) == 1
                                               else el))
        elif line.startswith("TER"):
            out.append("TER")
    out.append("END")
    return "\n".join(out) + "\n"


def import_target(src, rundir, log=print):
    """Copy a target file into the run and unpack what the run needs from it.

    Returns (trg, prepared_pdb, pdbqt, center, size). The file keeps its own
    name: ADCP looks for maps under the folder named after the file.
    """
    why = target.check(src)
    if why:
        raise ValueError("the target file %s is unusable: %s"
                         % (os.path.basename(src), why))
    tdir = os.path.join(rundir, "target")
    rdir = os.path.join(rundir, "receptor")
    os.makedirs(tdir, exist_ok=True)
    os.makedirs(rdir, exist_ok=True)
    trg = os.path.join(tdir, os.path.basename(src))
    if os.path.abspath(trg) != os.path.abspath(src):
        shutil.copy(src, trg)
    z = zipfile.ZipFile(trg)
    rec = next(n for n in z.namelist() if n.endswith("rigidReceptor.pdbqt"))
    text = z.read(rec).decode("utf8", "replace")
    pdbqt = os.path.join(rdir, "receptor.pdbqt")
    open(pdbqt, "w").write(text)
    prepared = os.path.join(rdir, "prepared.pdb")
    open(prepared, "w").write(pdbqt_to_pdb(text))
    centre, size = read_box(trg)
    with open(os.path.join(rdir, "preparation.log"), "w") as fh:
        fh.write("Receptor preparation\n====================\n\n"
                 "in  : %s\n\n"
                 "  target file: done; receptor, maps and box taken from %s, "
                 "made outside this workflow (for example in AGFRGUI). This "
                 "package's receptor preparation was not run.\n"
                 % (src, os.path.basename(src)))
    log("  target file %s: box center %.2f %.2f %.2f, size %.1f %.1f %.1f A"
        % ((os.path.basename(src),) + tuple(centre) + tuple(size)))
    return trg, prepared, pdbqt, centre, size
