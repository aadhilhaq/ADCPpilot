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
                 "x", "y", "z", "element", "hetatm")

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
                skipping = int(L[10:14]) != model
                continue
            if rec == "ENDMDL" and seen_model:
                skipping = True
                continue
            if skipping or rec not in ("ATOM  ", "HETATM"):
                continue
            try:
                out.append(Atom(serial=int(L[6:11]), name=L[12:16].strip(),
                                resname=L[17:20].strip(), chain=L[21].strip() or "A",
                                resnum=int(L[22:26]), icode=L[26].strip(),
                                x=float(L[30:38]), y=float(L[38:46]),
                                z=float(L[46:54]),
                                element=L[76:78].strip(), hetatm=rec == "HETATM"))
            except ValueError:
                continue
    return out


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
    """(centre, size, note), by the rule the run itself uses (receptor.box):
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
    note = ("%d residue(s) found, span %.1f A, box %.1f x %.1f x %.1f A"
            % ((found, max(spans)) + size))
    if missed:
        note += "; NOT in the structure: " + ", ".join(str(m) for m in missed)
    return ctr, size, note


def distance(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


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
            fh.write("%-6s%5d %-4s %3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f"
                     "          %2s\n"
                     % ("HETATM" if a.hetatm else "ATOM", n,
                        (" " + a.name) if len(a.name) < 4 else a.name,
                        a.resname[:3], a.chain[:1], a.resnum, a.icode or " ",
                        a.x, a.y, a.z, 1.0, 0.0, (a.element or "")[:2]))
            last = a.chain
        fh.write("TER\nEND\n")
    return n
