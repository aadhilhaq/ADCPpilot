"""Reading a receptor, and turning chosen residues into a docking box.

By Aadhil Haq

Deliberately small: enough PDB and mmCIF to list the chains, name the
residues, and find where a set of them sits in space. Anything more belongs
to a structure library, and a tool people install on a cluster is better off
without one it does not need.
"""
import gzip
import math
import os
import re

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "ASH": "D", "GLH": "E", "LYN": "K",
}
NUCLEIC = {"A", "C", "G", "U", "T", "DA", "DC", "DG", "DT", "DU"}
WATER = {"HOH", "WAT", "DOD", "H2O"}


class Atom:
    __slots__ = ("serial", "name", "resname", "chain", "resnum", "icode",
                 "x", "y", "z", "element", "hetatm", "altloc")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    @property
    def key(self):
        return (self.chain, self.resnum, self.icode or "")

    def is_hydrogen(self):
        e = (self.element or "").strip().upper()
        if e:
            return e == "H" or e == "D"
        n = (self.name or "").strip()
        return bool(n) and n[0] in "HD" and not n[:2] in ("HG", "HO", "HF")


def _open(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


# What this module can read, in one place. The file browser filter, the mmCIF
# test and /api/rawstructure all read these lists, so they cannot disagree.
#
# gzip is handled by _open for every format, so each plain suffix has a .gz
# twin.
CIF_SUFFIXES = (".cif", ".mmcif")
PDB_SUFFIXES = (".pdb", ".ent")
READABLE = tuple(s + g for s in PDB_SUFFIXES + CIF_SUFFIXES
                 for g in ("", ".gz"))


def is_cif(path):
    """Is this an mmCIF, by name? Compressed or not."""
    low = path.lower()
    return any(low.endswith(s + g) for s in CIF_SUFFIXES for g in ("", ".gz"))


def is_structure(path):
    """Is this a file this module will read? By name only."""
    return path.lower().endswith(READABLE)


def read(path, keep_water=False, keep_hetatm=True, model=1):
    """[Atom] from a PDB or mmCIF file. Only the first model, by default:
    an NMR ensemble docked as though it were one structure is twenty copies
    of the receptor in the same box."""
    if is_cif(path):
        atoms = _read_cif(path)
    else:
        atoms = _read_pdb(path, model)
    out = []
    for a in atoms:
        if a.resname in WATER and not keep_water:
            continue
        if a.hetatm and not keep_hetatm and a.resname not in AA3:
            continue
        out.append(a)
    return out


def _read_pdb(path, model):
    out, seen_model, skipping = [], False, False
    with _open(path) as fh:
        for L in fh:
            rec = L[:6]
            if rec == "MODEL ":
                seen_model = True
                # The number belongs in columns 11 to 14, and some writers
                # put it anywhere after the record name.
                try:
                    number = int(L[10:14])
                except ValueError:
                    try:
                        number = int((L[6:].split() or ["0"])[0])
                    except ValueError:
                        number = model
                skipping = number != model
                continue
            if rec == "ENDMDL" and seen_model:
                skipping = True
                continue
            if skipping or rec not in ("ATOM  ", "HETATM"):
                continue
            try:
                resnum = _number(L[22:26])
                x, y, z = _coords(L)
            except ValueError:
                continue
            out.append(Atom(serial=_serial(L[6:11]), name=L[12:16].strip(),
                            resname=L[17:20].strip(), chain=L[21].strip() or "A",
                            resnum=resnum, icode=L[26:27].strip(),
                            altloc=L[16:17].strip(), x=x, y=y, z=z,
                            element=L[76:78].strip(), hetatm=rec == "HETATM"))
    return out


def _serial(text):
    """An atom serial, which nothing here relies on, so it never costs the
    atom. Past 99,999 atoms OpenMM, VMD and others write it in hexadecimal
    ("186a0"), and some write stars. An atom whose serial could not be read
    used to be dropped, which left a large receptor's last atoms out of
    everything read from it, its surface included."""
    s = text.strip()
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return int(s, 16)
    except ValueError:
        return 0


def _number(text):
    """A residue number: decimal, or hexadecimal past 9,999 residues, which
    is how OpenMM writes one it cannot keep."""
    s = text.strip()
    try:
        return int(s)
    except ValueError:
        return int(s, 16)


_FLOAT = re.compile(r"-?\d+\.\d+")


def _coords(L):
    """x, y and z from their fixed columns, which hold three numbers even
    when they run together ("-100.123-200.456-300.789"). A line whose writer
    let a wide number push the columns along is read by the numbers in it
    instead, rather than lost."""
    try:
        return float(L[30:38]), float(L[38:46]), float(L[46:54])
    except ValueError:
        pass
    got = _FLOAT.findall(L[30:66])
    if len(got) < 3:
        raise ValueError("no coordinates in %r" % L[30:54])
    return float(got[0]), float(got[1]), float(got[2])


def _read_cif(path):
    """The atom_site loop, and only that."""
    cols, rows, in_loop, in_header = [], [], False, False
    with _open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("loop_"):
                if cols and rows:
                    break
                in_loop, in_header, cols = True, True, []
                continue
            if in_loop and s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1].split()[0])
                in_header = True
                continue
            if in_loop and cols and in_header and not s.startswith("_"):
                in_header = False
            if cols and not in_header:
                if not s or s.startswith("#") or s.startswith("loop_"):
                    if rows:
                        break
                    in_loop, cols = False, []
                    continue
                rows.append(s.split())
    if not cols or not rows:
        return []
    ix = {c: i for i, c in enumerate(cols)}

    def g(r, *names, **kw):
        for n in names:
            if n in ix and ix[n] < len(r):
                v = r[ix[n]]
                return kw.get("default") if v in (".", "?") else v
        return kw.get("default")

    out = []
    for r in rows:
        if len(r) < len(cols):
            continue
        try:
            out.append(Atom(
                serial=int(g(r, "id", default="0") or 0),
                name=g(r, "label_atom_id", "auth_atom_id", default="").strip('"'),
                resname=g(r, "label_comp_id", "auth_comp_id", default=""),
                chain=(g(r, "auth_asym_id", "label_asym_id", default="A") or "A"),
                resnum=int(float(g(r, "auth_seq_id", "label_seq_id", default="0") or 0)),
                icode=(g(r, "pdbx_PDB_ins_code", default="") or "").strip(),
                altloc=(g(r, "label_alt_id", default="") or "").strip(),
                x=float(g(r, "Cartn_x")), y=float(g(r, "Cartn_y")),
                z=float(g(r, "Cartn_z")),
                element=g(r, "type_symbol", default=""),
                hetatm=(g(r, "group_PDB", default="ATOM") == "HETATM")))
        except (TypeError, ValueError):
            continue
    return out


def chains(atoms):
    """[{id, kind, residues, first, last, sequence}] in file order."""
    order, by = [], {}
    for a in atoms:
        if a.chain not in by:
            by[a.chain] = {}
            order.append(a.chain)
        by[a.chain].setdefault(a.key, a.resname)
    out = []
    for c in order:
        keys = sorted(by[c], key=lambda k: (k[1], k[2]))
        names = [by[c][k] for k in keys]
        prot = sum(1 for n in names if n in AA3)
        nuc = sum(1 for n in names if n in NUCLEIC)
        kind = ("protein" if prot >= max(1, len(names) // 2)
                else "nucleic" if nuc >= max(1, len(names) // 2) else "other")
        out.append({"id": c, "kind": kind, "residues": len(keys),
                    "first": keys[0][1], "last": keys[-1][1],
                    "sequence": "".join(AA3.get(n, "X") for n in names)
                    if kind == "protein" else ""})
    return out


def residues(atoms, chain=None):
    """[{chain, resnum, icode, resname, code}] for a chain, or for all."""
    seen, out = set(), []
    for a in atoms:
        if chain and a.chain != chain:
            continue
        if a.key in seen:
            continue
        seen.add(a.key)
        out.append({"chain": a.chain, "resnum": a.resnum, "icode": a.icode or "",
                    "resname": a.resname, "code": AA3.get(a.resname, "X")})
    return out


def parse_hotspots(text, default_chain=None):
    """'A:41, A:45, 78, A:98-112' -> [(chain, resnum)].

    Delegated to the pipeline's own hotspot parser rather than reimplemented,
    so the form and the run always read a site the same way, ranges included.
    """
    from ._vendor.receptor import parse_hotspots as _parse
    return _parse(text, default_chain or "")


def centre_of(atoms, hotspots):
    """The middle of the chosen residues' heavy atoms, and how many were
    found. A hotspot the structure does not contain is the commonest setup
    mistake there is, so the count is returned rather than assumed."""
    want = {(c, int(n)) for c, n in hotspots}
    pts = [(a.x, a.y, a.z) for a in atoms
           if not a.is_hydrogen()
           and ((a.chain, a.resnum) in want
                or (None, a.resnum) in want)]
    if not pts:
        return None, 0, sorted({r for _, r in want})
    n = len(pts)
    ctr = (sum(p[0] for p in pts) / n, sum(p[1] for p in pts) / n,
           sum(p[2] for p in pts) / n)
    hit = {(a.chain, a.resnum) for a in atoms
           if (a.chain, a.resnum) in want or (None, a.resnum) in want}
    missed = sorted({r for c, r in want} - {r for c, r in hit})
    return ctr, len(hit), missed


def box_for(atoms, hotspots, padding=8.0, minimum=20.0, maximum=30.0):
    """(center, size, note), by the rule the run itself uses (receptor.box):
    along each axis, the span of the chosen residues plus padding on each
    side, clamped. Never smaller than `minimum`, because a handful of
    residues in a line would otherwise give a box a peptide cannot fit in,
    and never larger than `maximum`, because map building and the search
    both scale with the volume and a box the size of the protein is a blind
    docking run. The page draws this box, so it must be the one docked into.
    """
    ctr, found, missed = centre_of(atoms, hotspots)
    if ctr is None:
        return None, None, "none of the chosen residues are in the structure"
    want = {(c, int(n)) for c, n in hotspots}
    pts = [(a.x, a.y, a.z) for a in atoms
           if not a.is_hydrogen()
           and ((a.chain, a.resnum) in want or (None, a.resnum) in want)]
    spans = [max(p[i] for p in pts) - min(p[i] for p in pts) for i in range(3)]
    size = tuple(max(minimum, min(maximum, s + 2 * padding)) for s in spans)
    # The box's size is not repeated here: the page and the run's log both
    # state it on the next line, from the same numbers.
    note = ("%d residue(s) found, spanning %.1f \u00c5"
            % (found, max(spans)))
    if missed:
        note += "; NOT in the structure: " + ", ".join(str(m) for m in missed)
    return ctr, size, note


def clipped(points, padding, maximum):
    """[(axis, span, padding left)] for each axis along which the site and
    its padding are wider than the largest box allowed, so the box is held
    to the maximum there. Where the span itself is wider than the maximum,
    part of the site lies outside the box and cannot be docked into; where
    only the padding is cut, every residue is still inside it."""
    out = []
    if not points:
        return out
    for i, axis in enumerate("xyz"):
        span = max(p[i] for p in points) - min(p[i] for p in points)
        if span + 2 * padding > maximum + 0.05:
            out.append((axis, round(span, 1), round(max(0.0, (maximum - span) / 2.0), 1)))
    return out


def distance(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


# ---- a box from a ligand or peptide in the structure --------------------
#
# AGFR can place its box around a ligand (-b ligand), and AGFR 1.1 can pick
# that ligand out of a PDB itself (--pdbLigand). Only on a PDB it prepares
# itself, though: --pdbLigand belongs to the path where AGFR protonates the
# receptor with reduce and charges it, which would replace this tool's
# preparation. And a peptide picked that way stays in the receptor: a
# protein chain is flagged as receptor before the selection flags it as the
# ligand, and on 1YCR with --pdbLigand "chain B" the receptor in the target
# held chain B's 132 atoms, so the maps were built with the peptide still in
# its site. So the box is worked out here by AGFR's own rule for -b ligand
# (setBoxForCoords: the middle of the ligand's extent, and the extent plus
# the padding on each side) and given to agfr as a box, and the ligand is
# kept out of the receptor here.

#: Hetero groups that are crystallization additives (AGFR's own list, from
#: recFromPDB.PDBProfiler._additives). Offered, but after the ligands.
ADDITIVES = {"P4G", "EOH", "1PE", "DMS", "ACT", "ACY", "PEG", "PGE", "BME",
             "FMT", "BUD", "GOL", "MRY", "BU3", "PO4", "MLI", "PDO", "PG4",
             "SO4", "15P", "CO3", "PE8", "PE4", "EDO"}
#: A protein chain this short is offered as a peptide that could set the box.
PEPTIDE_MAX = 50
#: Residues with a heavy atom this close to the ligand, or to an AutoSite
#: pocket's points, are the site its box is built around.
CONTACT = 4.0


def ligands(atoms):
    """What in a structure could set the box: [{id, kind, label, atoms}].

    A peptide is named by its chain, like B; a hetero group as
    CHAIN:NAME:NUMBER, like A:NUT:201. Waters are left out and additives
    come last.
    """
    groups, order = {}, []
    for a in atoms:
        if not a.hetatm or a.resname in AA3 or a.resname in WATER:
            continue
        k = (a.chain, a.resname, a.resnum, a.icode or "")
        if k not in groups:
            groups[k] = 0
            order.append(k)
        if not a.is_hydrogen():
            groups[k] += 1
    out = []
    for ch in chains(atoms):
        if ch["kind"] != "protein" or ch["residues"] > PEPTIDE_MAX:
            continue
        seq = ch["sequence"]
        out.append({"id": ch["id"], "kind": "peptide",
                    "atoms": sum(1 for a in atoms if a.chain == ch["id"]
                                 and not a.is_hydrogen()),
                    "label": "chain %s, %d residue%s%s"
                    % (ch["id"], ch["residues"],
                       "" if ch["residues"] == 1 else "s",
                       " (%s)" % seq if seq and len(seq) <= 30 else "")})
    het = []
    for k in order:
        c, name, num, ic = k
        kind = ("additive" if name in ADDITIVES
                else "ion" if groups[k] == 1 else "ligand")
        het.append({"id": "%s:%s:%d%s" % (c, name, num, ic), "kind": kind,
                    "atoms": groups[k],
                    "label": "%s %s:%d%s, %d atom%s%s"
                    % (name, c, num, ic, groups[k],
                       "" if groups[k] == 1 else "s",
                       ", an additive" if kind == "additive" else "")})
    het.sort(key=lambda r: r["kind"] == "additive")
    return out + het


def parse_ligand(spec):
    """("chain", "B") or ("het", (chain, name, number, icode)), or None."""
    s = str(spec or "").strip()
    if not s:
        return None
    parts = s.split(":")
    if len(parts) == 1 and len(s) <= 4:
        return ("chain", s)
    if len(parts) == 3:
        c, name, num = (p.strip() for p in parts)
        ic = ""
        if num and num[-1].isalpha():
            num, ic = num[:-1], num[-1]
        try:
            return ("het", (c, name.upper(), int(num), ic))
        except ValueError:
            return None
    return None


def ligand_atoms(atoms, spec):
    """The heavy atoms of the ligand `spec` names, or []."""
    p = parse_ligand(spec)
    if not p:
        return []
    if p[0] == "chain":
        return [a for a in atoms if a.chain == p[1] and not a.is_hydrogen()]
    c, name, num, ic = p[1]
    return [a for a in atoms if a.chain == c and a.resname == name
            and a.resnum == num and (a.icode or "") == ic
            and not a.is_hydrogen()]


def ligand_problems(atoms, spec, chains_keep=(), keep_hetatms=False):
    """[str] why `spec` cannot set the box for a receptor of `chains_keep`
    (empty: every chain but a peptide ligand's own)."""
    p = parse_ligand(spec)
    if not p:
        return ["%r does not name a ligand. Write a peptide's chain, like B, "
                "or a hetero group as CHAIN:NAME:NUMBER, like A:NUT:201."
                % str(spec)]
    if not ligand_atoms(atoms, spec):
        have = ", ".join(r["id"] for r in ligands(atoms)[:12]) or "nothing"
        return ["%s is not in the structure. What it holds that could set "
                "the box: %s." % (spec, have)]
    keep = [c for c in (chains_keep or []) if c]
    if p[0] == "chain":
        if p[1] in keep:
            return ["Chain %s sets the box and is also one of the receptor's "
                    "chains, so the peptide would stay in the site it marks. "
                    "Untick chain %s among the chains to keep." % (p[1], p[1])]
        return []
    chain = p[1][0]
    if keep_hetatms and (not keep or chain in keep):
        return ["Keep heteroatoms is on, so %s would stay in the receptor and "
                "fill the site its box is built around. Turn Keep heteroatoms "
                "off, or place the box another way." % spec]
    return []


def ligand_box(lig, padding=8.0, minimum=20.0, maximum=30.0):
    """(center, size, note) for a box around a ligand, by AGFR's rule for
    -b ligand: centered on the middle of its extent, the extent plus the
    padding on each side. Then clamped, as a box around residues is."""
    if not lig:
        return None, None, "the ligand has no atoms in this structure"
    lo = [min(getattr(a, k) for a in lig) for k in "xyz"]
    hi = [max(getattr(a, k) for a in lig) for k in "xyz"]
    centre = tuple((lo[i] + hi[i]) / 2.0 for i in range(3))
    spans = [hi[i] - lo[i] for i in range(3)]
    size = tuple(max(minimum, min(maximum, s + 2 * padding)) for s in spans)
    note = ("%d ligand atoms, spanning %.1f \u00c5" % (len(lig), max(spans)))
    return centre, size, note


def _near(atoms, pts, cutoff, chains_keep, skip=()):
    keep = set(chains_keep or [])
    c2 = cutoff * cutoff
    found = set()
    for a in atoms:
        if id(a) in skip or a.is_hydrogen() or a.resname in WATER:
            continue
        if keep and a.chain not in keep:
            continue
        if a.hetatm and a.resname not in AA3:
            continue
        for p in pts:
            if ((a.x - p[0]) ** 2 + (a.y - p[1]) ** 2
                    + (a.z - p[2]) ** 2) <= c2:
                found.add((a.chain, a.resnum))
                break
    return sorted(found)


def contacts(atoms, lig, cutoff=CONTACT, chains_keep=None):
    """[(chain, resnum)] of the receptor residues with a heavy atom within
    `cutoff` of the ligand: the site its box is built around, which the
    results then measure each pose against."""
    return _near(atoms, [(a.x, a.y, a.z) for a in lig], cutoff, chains_keep,
                 skip={id(a) for a in lig})


def near_points(atoms, points, cutoff=CONTACT, chains_keep=None):
    """[(chain, resnum)] of the residues with a heavy atom within `cutoff`
    of any of `points`: the residues lining an AutoSite pocket."""
    return _near(atoms, points, cutoff, chains_keep)


def one_conformer(atoms):
    """The atoms with one alternate location per residue: those without one,
    and those of the first location the residue lists, normally A, which is
    the one 3Dmol draws. The site view is drawn from these and the surface
    is computed from them, so the two agree."""
    first = {}
    for a in atoms:
        alt = (getattr(a, "altloc", "") or "").strip()
        if not alt:
            yield a
            continue
        k = (a.chain, a.resnum, a.icode or "")
        if first.setdefault(k, alt) == alt:
            yield a


def _fixed(n, width):
    """A number that fits its columns: decimal while it can, then
    hexadecimal, as OpenMM writes them, so a number past the columns' reach
    never pushes the coordinates along."""
    if -(10 ** (width - 1)) < n < 10 ** width:
        return "%*d" % (width, n)
    return "%*x" % (width, n % (16 ** width))


def write_pdb(atoms, path, chains_keep=None):
    """Back out as PDB, which is what every tool downstream reads."""
    keep = set(chains_keep or [])
    n = 0
    with open(path, "w") as fh:
        last = None
        for a in atoms:
            if keep and a.chain not in keep:
                continue
            if last is not None and a.chain != last:
                fh.write("TER\n")
            n += 1
            fh.write("%-6s%5s %-4s %3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f"
                     "          %2s\n"
                     % ("HETATM" if a.hetatm else "ATOM", _fixed(n, 5),
                        (" " + a.name) if len(a.name) < 4 else a.name,
                        a.resname[:3], a.chain[:1], _fixed(a.resnum, 4),
                        a.icode or " ",
                        a.x, a.y, a.z, 1.0, 0.0, (a.element or "")[:2]))
            last = a.chain
        fh.write("TER\nEND\n")
    return n
