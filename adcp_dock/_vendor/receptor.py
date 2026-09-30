"""Receptor loading (.pdb / .cif / .mae), docking box construction and PDBQT preparation.

By Aadhil Haq
"""
from __future__ import annotations
import os
import re, shutil, subprocess
from dataclasses import dataclass, field
import numpy as np

BACKBONE = {"N", "CA", "C", "O", "OXT"}

def is_hydrogen(name: str, element: str = "") -> bool:
    """Robust hydrogen test. PDB and PDBQT name hydrogens '1HD2' as well as 'HD2', so a bare name[0]=='H' test misses roughly 5% of them.
    """
    e = (element or "").strip()
    if e:
        return e.upper().startswith("H")
    return name.strip().lstrip("0123456789").upper().startswith("H")

def _serial(n):
    """A PDB atom serial that fits in columns 7 to 11.

    %5d prints six characters at 100,000, and every later field then moves one
    column right: the atom name lands in altLoc and the chain letter in resSeq.
    Writing hydrogens roughly doubles the atom count, so a large complex can
    pass that limit. Hybrid-36 is what the PDB itself uses beyond it.
    """
    if n < 100000:
        return "%5d" % n
    n -= 100000
    digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    out = ""
    for _ in range(5):
        n, r = divmod(n, 36)
        out = digits[r] + out
    return "A" + out[1:] if n == 0 else "*****"


def pdb_line(serial: int, atom, chain=None) -> str:
    """One ATOM record. chain=None uses the atom's own chain.

    Writing two chains under one letter merges them into one molecule: tleap
    then reads a continuous chain and bonds the C-terminus of the first to the
    N-terminus of the second, inventing a peptide bond across an interface.
    """
    """A correctly column-aligned PDB ATOM record.

    Columns matter: 13-16 name, 17 altLoc, 18-20 resName, 22 chain, 23-26 resSeq,
    31-54 coordinates, 77-78 element. Getting 17 wrong makes every residue look
    like it has an alternate location, which breaks prepare_receptor.
    """
    el = (atom.element or atom.name[0]).strip()
    nm = atom.name
    name = nm.ljust(4) if len(nm) >= 4 else (" " + nm).ljust(4)
    return (f"ATOM  {_serial(serial)} {name}"
            f" {atom.resname:>3} {(chain if chain is not None else (getattr(atom, 'chain', '') or 'A'))}{atom.resnum:4d}    "
            f"{atom.xyz[0]:8.3f}{atom.xyz[1]:8.3f}{atom.xyz[2]:8.3f}"
            f"{1.00:6.2f}{0.00:6.2f}          {el:>2}\n")

# Residue numbers run to about 10,000 in the largest structures; 100,000
# leaves room and keeps the packed code readable when it is printed.
_CHAIN_STRIDE = 100000


@dataclass
class Atom:
    resnum: int
    resname: str
    name: str
    element: str
    xyz: np.ndarray
    # Which chain this atom came from.
    #
    # Defaulted so every existing construction of an Atom (the Maestro
    # reader, the tests, anything that builds one positionally) still works
    # and means "the one chain there was". It stops being cosmetic the moment
    # two chains are read at once: residue 56 of chain A and residue 56 of
    # chain B are different residues, and everything downstream matched
    # hotspots on the number alone.
    chain: str = ""

    @property
    def key(self) -> str:
        """Chain-qualified residue identity, e.g. "A:56"."""
        return "%s:%d" % (self.chain, self.resnum)
    @property
    def is_backbone(self) -> bool:
        return self.name in BACKBONE

@dataclass
class Receptor:
    path: str
    chain: str                      # the PRIMARY chain: bare numbers mean this
    atoms: list[Atom]
    hotspots: list[int]             # residue numbers, primary chain unless keyed
    shell: list[int]
    # Every chain actually loaded, primary first. One entry for a
    # single-chain run, which is what it was before.
    chains: list[str] = field(default_factory=list)
    # Hotspots as (chain, resnum). `hotspots` above keeps the bare numbers so
    # the thirty-odd places that only ever had one chain go on working; this
    # is what anything chain-aware should read.
    hotspot_keys: list[tuple] = field(default_factory=list)

    def __post_init__(self):
        if not self.chains:
            self.chains = [self.chain] if self.chain else []
        if not self.hotspot_keys:
            self.hotspot_keys = [(self.chain, r) for r in self.hotspots]

    @property
    def multi_chain(self) -> bool:
        return len(self.chains) > 1

    def chain_index(self) -> dict:
        """chain -> position, so a (chain, resnum) pair can become one int."""
        return {c: i for i, c in enumerate(self.chains)}

    def residue_codes(self) -> "np.ndarray":
        """One integer per atom identifying its residue AND its chain.

        numpy membership tests on the residue number alone are wrong the
        moment two chains are loaded: both chains' residue 56 match. This
        packs the chain into the same integer so np.isin stays fast and
        starts being right.
        """
        idx = self.chain_index()
        return np.array([idx.get(a.chain, 0) * _CHAIN_STRIDE + a.resnum
                         for a in self.atoms], dtype=np.int64)

    def codes_for(self, keys) -> "np.ndarray":
        """The same integers for a list of (chain, resnum) pairs."""
        idx = self.chain_index()
        return np.array([idx.get(c, 0) * _CHAIN_STRIDE + r for c, r in keys],
                        dtype=np.int64)

    @property
    def coords(self) -> np.ndarray:
        return np.array([a.xyz for a in self.atoms])

    def atoms_of(self, resnums) -> list[Atom]:
        """Atoms of the given residues.

        Accepts bare numbers, which match those residues on any loaded chain, or
        (chain, resnum) pairs, which match exactly one residue. Chain-aware code
        should pass pairs: on a two-chain receptor a bare 56 matches residue 56 of
        both chains, and a docking box centred on both sits between them, in the
        solvent.
        """
        s = set(resnums)
        if s and all(isinstance(x, tuple) for x in s):
            return [a for a in self.atoms if (a.chain, a.resnum) in s]
        return [a for a in self.atoms if a.resnum in s]

def hetatms(path: str, chain: str) -> dict:
    """Non-protein residues present in the file.

    These are not carried into the receptor: waters, ions, metals and cofactors
    are dropped. On a target with a catalytic metal or a structural ion that
    matters, which is why it is reported.
    """
    from collections import Counter
    c = Counter()
    if path.lower().endswith((".cif", ".mmcif")):
        # Reported for mmCIF as well as PDB: the CIF reader drops the
        # same residues, and returning nothing here meant no warning
        # on that format.
        try:
            from Bio.PDB.MMCIFParser import MMCIFParser
            st = MMCIFParser(QUIET=True).get_structure("r", path)
            model = next(iter(st))
        except Exception:
            return {}
        for ch in model:
            if chain and ch.id != chain:
                continue
            for res in ch:
                if res.id[0] == " ":
                    continue                 # a residue of the chain itself
                nm = res.get_resname().strip()
                if nm != "HOH":
                    c[nm] += len([a for a in res])
        return dict(c)
    for l in open(path):
        if l.startswith("HETATM") and l[21] in (chain, " "):
            nm = l[17:20].strip()
            if nm != "HOH":
                c[nm] += 1
    return dict(c)


def chain_breaks(atoms):
    """Consecutive residue numbers with a gap. Unmodelled loops leave artificial
    cavities that docking will happily fill."""
    rn = sorted({a.resnum for a in atoms})
    return [(rn[i], rn[i + 1]) for i in range(len(rn) - 1) if rn[i + 1] - rn[i] > 1]


# Residues whose number alone does not identify them, per (file, chain).
# Filled by the PDB reader; read before a hotspot is chosen by number.
_INSERTION_CODES = {}


def insertion_codes(path, chain):
    """[(resnum, code)] for residues this chain distinguishes by a letter."""
    return list(_INSERTION_CODES.get((os.path.abspath(str(path)), chain), []))


def chains_wanted(chain) -> list:
    """The chains a config's `chain` field asks for, in order.

    "A" is one chain and stays one chain. "A,B" is two. The first is the
    PRIMARY: a hotspot written as a bare number belongs to it, which is what
    every config written before this meant and still means.
    """
    if isinstance(chain, (list, tuple)):
        return [str(c).strip() for c in chain if str(c).strip()]
    return [c.strip() for c in str(chain or "").split(",") if c.strip()]


def parse_hotspots(spec, primary: str) -> list:
    """[(chain, resnum)] from "56,123" or "A:56,B:45" or a mix of both.

    A bare number means the primary chain. That is not a convenience: it is
    what every existing config says, and those configs must not change
    meaning because the parser learned a new syntax.
    """
    out = []
    items = spec if isinstance(spec, (list, tuple)) else str(spec or "").split(",")
    for item in items:
        s = str(item).strip()
        if not s:
            continue
        if ":" in s:
            c, _, n = s.partition(":")
            c = c.strip() or primary
        else:
            c, n = primary, s
        n = str(n).strip()
        # A range such as "98-112", read the same way as the
        # configuration parser reads it, so a site means the same
        # thing wherever the text comes from.
        m = re.match(r"^(-?\d+)\s*[-\u2013]\s*(-?\d+)$", n)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            if hi < lo:
                lo, hi = hi, lo
            out.extend((c, r) for r in range(lo, hi + 1))
            continue
        try:
            out.append((c, int(n)))
        except ValueError:
            raise ValueError(
                "hotspot %r is not a residue number, a range like 98-112, or "
                "a chain and a residue number written as CHAIN:NUMBER" % s)
    return out


def _read_pdb(path: str, chain) -> list[Atom]:
    """Heavy atoms of one chain, from the first model only.

    An NMR structure holds many models of the same protein one after another
    in the same file. Reading them all would stack every model on top of the
    others, making the clash test meaningless, the contact counts many times
    too large, and the docking box the union of every model. _read_cif below
    also takes one model, so both formats give the same structure.
    """
    # A blank chain column is a chain. Many files have one, for example
    # anything written by a modelling tool, and pick.chains() reports it
    # as "_" so that it can appear in a menu; "_" must therefore match a
    # blank column.
    #
    # No chain named means every chain, as the interface says ("Chains to
    # keep; none ticked means all").
    _want = {(" " if c == "_" else c) for c in chains_wanted(chain)}
    _all = not _want
    out = []
    icodes = set()
    for l in open(path):
        if l.startswith("ENDMDL"):
            break                            # the first model is the receptor
        if not l.startswith("ATOM"):
            continue
        if not _all and l[21] not in _want:
            continue
        if l[16] not in (" ", "A"):          # keep one alternate location only
            continue
        name = l[12:16].strip()
        elem = l[76:78].strip()
        if is_hydrogen(name, elem):
            continue
        try:
            xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        except ValueError:
            continue
        # Column 27 is the insertion code, so 15A and 15B, or 100C to
        # 100J, share a residue number. Selecting that number as a
        # hotspot selects all of them at once, and the box centres on
        # their combined centroid. The number is kept as it is so
        # nothing downstream changes, and the codes are collected so
        # the ambiguity can be reported.
        ic = l[26] if len(l) > 26 else " "
        if ic.strip():
            icodes.add((int(l[22:26]), ic))
        out.append(Atom(int(l[22:26]), l[17:20].strip(), name, elem, xyz,
                        "_" if l[21] == " " else l[21]))
    _INSERTION_CODES[os.path.abspath(path), str(chain)] = sorted(icodes)
    return out

def _read_cif(path: str, chain) -> list[Atom]:
    """Heavy atoms of the wanted chains, from the first model.

    Every atom carries its chain, which matters as soon as a receptor is more
    than one chain.
    """
    from Bio.PDB.MMCIFParser import MMCIFParser
    want = chains_wanted(chain)
    st = MMCIFParser(QUIET=True).get_structure("r", path)
    model = next(iter(st))
    have = [c.id for c in model]
    # Same rule as the PDB reader: nothing asked for is everything.
    if not want:
        want = list(have)
    absent = [c for c in want if c not in have]
    if absent:
        raise ValueError("chain(s) %s not in %s; found %s"
                         % (", ".join(absent), path, ", ".join(have)))
    out, icodes = [], set()
    for cid in want:
        for res in model[cid]:
            if res.id[0] != " ":       # skip waters / heteroatoms
                continue
            if str(res.id[2]).strip():
                icodes.add((res.id[1], str(res.id[2])))
            for at in res:
                elem = (at.element or "").strip()
                if is_hydrogen(at.get_name(), elem):
                    continue
                out.append(Atom(res.id[1], res.get_resname().strip(),
                                at.get_name().strip(), elem,
                                np.array(at.coord, float), cid))
    _INSERTION_CODES[os.path.abspath(path), str(chain)] = sorted(icodes)
    return out

def load(cfg) -> Receptor:
    ext = os.path.splitext(cfg.receptor)[1].lower()
    # Maestro files come straight out of Schrodinger; RDKit reads them,
    # so no Schrodinger installation or licence is needed here.
    if ext in (".mae", ".maegz"):
        from .read_maestro import read_maestro
        atoms = read_maestro(cfg.receptor, cfg.chain, Atom, is_hydrogen)
    else:
        atoms = (_read_cif if ext in (".cif", ".mmcif") else _read_pdb)(cfg.receptor, cfg.chain)
    if not atoms:
        raise ValueError(
            f"no atoms read for chain '{cfg.chain}' from {cfg.receptor}")
    want = chains_wanted(cfg.chain)
    # The readers above keep every chain when none was named, so record
    # which chains that turned out to be. Otherwise `want` stays empty,
    # and everything below (primary, got, missing_chains, the Receptor's
    # own .chains, and the chain list handed to the prepare stage)
    # behaves as though no chain was read.
    if not want:
        want = sorted({a.chain for a in atoms})
    primary = want[0] if want else ""
    got = [c for c in want if any(a.chain == c for a in atoms)] or want
    missing_chains = [c for c in want if c not in got]
    if missing_chains:
        raise ValueError(
            f"chain(s) {', '.join(missing_chains)} are not in {cfg.receptor}; "
            f"found {', '.join(sorted({a.chain for a in atoms})) or 'none'}")

    # Hotspots by (chain, residue), so residue 56 of chain A and residue 56
    # of chain B are told apart. A bare number still means the primary chain.
    keys = parse_hotspots(getattr(cfg, "hotspots", []) or [], primary)
    present = {(a.chain, a.resnum) for a in atoms}
    absent = [k for k in keys if k not in present]
    if absent:
        spans = {c: (min(a.resnum for a in atoms if a.chain == c),
                     max(a.resnum for a in atoms if a.chain == c))
                 for c in got}
        raise ValueError(
            "hotspot residues absent: %s\n%s"
            % (", ".join("%s:%d" % k for k in absent),
               "; ".join("chain %s spans %d-%d" % (c, lo, hi)
                         for c, (lo, hi) in sorted(spans.items()))))

    keyset = set(keys)
    hot = [a for a in atoms if (a.chain, a.resnum) in keyset]
    hxyz = np.array([a.xyz for a in hot])
    allxyz = np.array([a.xyz for a in atoms])
    d = np.linalg.norm(allxyz[:, None, :] - hxyz[None, :, :], axis=-1).min(axis=1)
    shell_keys = sorted({(a.chain, a.resnum) for a, dd in zip(atoms, d)
                         if dd <= cfg.shell_radius} - keyset)
    # The bare-number lists stay what they were for a single-chain run, which
    # is what almost everything downstream still reads.
    rec = Receptor(cfg.receptor, primary, atoms,
                   [r for _c, r in keys], [r for _c, r in shell_keys],
                   chains=got, hotspot_keys=keys)
    rec.shell_keys = shell_keys
    # What the shell was built with, so the graded nearness ramp in
    # contacts.py uses the same number instead of its own constant.
    rec.shell_radius = float(cfg.shell_radius)
    return rec

def box(rec: Receptor, cfg) -> tuple[np.ndarray, np.ndarray]:
    """Return (center, size) in Angstrom, centred on the hotspot centroid."""
    h = np.array([a.xyz for a in rec.atoms_of(rec.hotspot_keys or rec.hotspots)])
    center = h.mean(axis=0)
    extent = h.max(axis=0) - h.min(axis=0)
    size = np.clip(extent + 2 * cfg.box_padding, cfg.box_min, cfg.box_max)
    return center, size


def chains_in_box(path: str, center, size, exclude=(), margin: float = 0.0,
                  min_atoms: int = 0):
    """Chains with heavy atoms inside the docking box -> {chain: n_atoms}.

    `min_atoms` drops chains that only graze the box: a rectangular box can clip
    a residue or two of a chain far from the site, which is not the same as a
    chain that forms part of the site.

    A box built around the hotspots of a multi-chain structure often reaches
    into a neighbouring chain, because an interface groove has two sides.
    Docking into that groove with only one side loaded places the peptide
    against a hole: it scores well because half the wall is missing, and the
    rescoring measures the same missing wall.

    Reads the whole file rather than the loaded Receptor, because the point is
    to find chains that were not loaded.
    """
    lo = np.asarray(center) - np.asarray(size) / 2.0 - margin
    hi = np.asarray(center) + np.asarray(size) / 2.0 + margin
    found = {}
    ext = os.path.splitext(str(path))[1].lower()
    if ext in (".cif", ".mmcif"):
        from Bio.PDB.MMCIFParser import MMCIFParser
        st = MMCIFParser(QUIET=True).get_structure("r", path)
        model = next(iter(st))
        for ch in model:
            if ch.id in exclude:
                continue
            n = 0
            for res in ch:
                for at in res:
                    if is_hydrogen(at.get_name(), at.element or ""):
                        continue
                    xyz = np.asarray(at.get_coord(), float)
                    if np.all(xyz >= lo) and np.all(xyz <= hi):
                        n += 1
            if n >= max(1, min_atoms):
                found[ch.id] = n
        return found
    for l in open(path, errors="replace"):
        if l.startswith("ENDMDL"):
            break
        if not l.startswith("ATOM"):
            continue
        c = "_" if l[21] == " " else l[21]
        if c in exclude:
            continue
        if l[16] not in (" ", "A"):
            continue
        if is_hydrogen(l[12:16].strip(), l[76:78].strip()):
            continue
        try:
            xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        except ValueError:
            continue
        if np.all(xyz >= lo) and np.all(xyz <= hi):
            found[c] = found.get(c, 0) + 1
    return {c: n for c, n in found.items() if n >= max(1, min_atoms)}


def load_for_box(cfg, log=print):
    """Load the receptor, then decide what else the box reaches into.

    Two passes, because the box is built from the hotspots, and only then can
    it say which other chains it covers.

    What to do about those chains depends on the structure, so it is not
    decided here. When the site is a groove between two chains of one
    assembly, leaving the second chain out means docking against a hole. When
    the other chain is the partner a designed peptide should displace, loading
    it would dock every peptide against that partner. The default reports and
    does nothing, and `box_chains` decides:

      "ignore": report what was found and load nothing more.
      "include": load every chain that is really in the box.
      "B,C": load exactly these, whatever the box says.

    "Really in the box" has a floor, because a rectangular box can clip a
    residue or two of a distant chain. box_chains_min_atoms is that floor,
    about two residues by default; a chain named explicitly is loaded whatever
    its count.
    """
    rec = load(cfg)
    centre, size = box(rec, cfg)
    want = str(getattr(cfg, "box_chains", "ignore") or "ignore").strip()
    floor = int(getattr(cfg, "box_chains_min_atoms", 20) or 0)
    found = chains_in_box(cfg.receptor, centre, size, exclude=set(rec.chains),
                          min_atoms=floor)
    if not found:
        return rec, centre, size
    listed = ", ".join("%s (%d atoms)" % (c, n)
                       for c, n in sorted(found.items(), key=lambda x: -x[1]))
    if want.lower() == "ignore":
        log("  the docking box also reaches chain(s) %s, which are NOT loaded."
            % listed)
        log("    set box_chains = include to dock against them, or name them "
            "(box_chains = %s). Leave it alone if they are the partner this "
            "peptide is meant to replace."
            % ",".join(sorted(found)))
        return rec, centre, size
    if want.lower() == "include":
        extra = sorted(found)
    else:
        extra = [c for c in chains_wanted(want) if c not in rec.chains]
        unseen = [c for c in extra if c not in found]
        if unseen:
            log("  box_chains names chain(s) %s, which have no atoms in the "
                "box; loading them anyway because they were asked for."
                % ", ".join(unseen))
    if not extra:
        return rec, centre, size
    import copy as _copy
    c2 = _copy.copy(cfg)
    c2.chain = ",".join(list(rec.chains) + extra)
    rec2 = load(c2)
    log("  loaded chain(s) %s as well: the box covers %s"
        % (", ".join(extra), listed))
    # The box is NOT recomputed. It is defined by the hotspots, and widening
    # it because more atoms are now in view would let it grow every pass.
    return rec2, centre, size


def hotspots_outside_box(rec: Receptor, center, size):
    """Hotspot residues with atoms outside the search box, and by how far.

    The box is clipped to box_max, and the result is compared against the
    residues it is meant to contain.
    """
    center = np.asarray(center, dtype=float)
    size = np.asarray(size, dtype=float)
    lo, hi = center - size / 2.0, center + size / 2.0
    keys = rec.hotspot_keys or [(rec.chain, r) for r in rec.hotspots]
    out = []
    for c, rn in keys:
        xyz = np.array([a.xyz for a in rec.atoms_of([(c, rn)])], dtype=float)
        if not len(xyz):
            continue
        past = np.maximum(lo - xyz, xyz - hi).max(axis=1)
        worst = float(past.max())
        if worst > 0:
            out.append((residue_label(rec, c, rn), round(worst, 1)))
    out.sort(key=lambda x: -x[1])
    return out


def residue_label(rec: Receptor, chain: str, resnum: int) -> str:
    """How a residue should be written back to whoever named it."""
    return str(resnum) if chain == rec.chain else "%s:%d" % (chain, resnum)

def _write_chain_pdb(rec: Receptor, workdir: str) -> str:
    """The loaded chains as a plain PDB, with no preparation applied.

    Each atom is written under its own chain, with a TER between chains. Two
    chains written under one letter would merge into a single molecule (see
    pdb_line). This is the fallback path, used when preparation is off or has
    failed, so it must be correct on its own.
    """
    path = os.path.join(workdir, "receptor_chain.pdb")
    with open(path, "w") as fh:
        last = None
        for i, a in enumerate(rec.atoms, 1):
            if last is not None and a.chain != last:
                fh.write("TER\n")
            fh.write(pdb_line(i, a))      # None: the atom's own chain
            last = a.chain
        fh.write("END\n")
    return path


def _run_prepare_receptor(prep: str, src: str, out: str, add_h: bool):
    """One prepare_receptor call. Leaves no stale output behind on failure."""
    if os.path.exists(out):
        os.remove(out)
    cmd = [prep, "-r", src, "-o", out] + (["-A", "hydrogens"] if add_h else [])
    return subprocess.run(cmd, capture_output=True, text=True)


def count_hydrogens(pdb_path: str) -> int:
    """Hydrogens in a PDB/PDBQT file, decided the way is_hydrogen() decides."""
    n = 0
    with open(pdb_path, errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")) and \
                    is_hydrogen(line[12:16], line[76:78]):
                n += 1
    return n


def prepare_pdbqt(rec: Receptor, workdir: str, tools,
                  prepared_pdb: str | None = None, log=print) -> str:
    """Write receptor.pdbqt using ADFRsuite prepare_receptor.

    If `prepared_pdb` carries hydrogens (the output of PDB2PQR does)
    prepare_receptor is NOT asked to add them, because that structure already
    has a complete, pH-appropriate set and geometric hydrogens would be worse.

    The hydrogens are counted rather than assumed. prep.run() returns
    prepared.pdb whether or not pdb2pqr ran: on failure it falls back to the
    un-protonated structure and returns it under exactly the same name, so the
    file's existence says nothing about whether it has been protonated.
    Withholding -A then wrote a receptor.pdbqt with no hydrogens on it, which
    prepare_receptor accepts in silence and AGFR rejects several minutes later
    with "It seems that hydrogen atoms are missing in the receptor": the
    symptom, reported by a different program, long after the cause, and with
    no mention of protonation at all.
    """
    os.makedirs(workdir, exist_ok=True)
    if prepared_pdb and os.path.isfile(prepared_pdb):
        src = prepared_pdb
        add_h = count_hydrogens(src) == 0
        if add_h:
            log("WARNING   the prepared structure has no hydrogens on it, so "
                "pdb2pqr did not run.")
            log("          Why is in receptor/preparation.log, on the "
                "'protonation' line.")
            log("          prepare_receptor will add them geometrically so the "
                "run can go on, but")
            log("          there are then no pH-dependent titration states: "
                "the run's pH is not")
            log("          reflected in the receptor. Fix the preparation and "
                "re-run if that matters.")
    else:
        src = _write_chain_pdb(rec, workdir)
        add_h = True
    out = os.path.join(workdir, "receptor.pdbqt")
    prep = tools.prepare_receptor
    if prep is None:
        raise RuntimeError(
            "prepare_receptor not found. Run ./install.sh to install ADFRsuite.")

    p = _run_prepare_receptor(prep, src, out, add_h)

    # Everything above can succeed and still produce a receptor nothing can
    # dock against, so the only question worth asking is asked directly. If
    # the answer is wrong, try the one other structure available (the chain
    # as it was read, without any preparation) because a prepared PDB that
    # defeats prepare_receptor is usually malformed in a way the plain
    # extract is not. Recovering here is what keeps a run going; stopping
    # would only trade AGFR's confusing failure for a clearer one.
    if not os.path.isfile(out) or count_hydrogens(out) == 0:
        if src != prepared_pdb:
            raise RuntimeError(
                "prepare_receptor could not produce a receptor with hydrogens "
                f"on it from {src}.\n" + (p.stdout or "") + (p.stderr or ""))
        log("WARNING   prepare_receptor could not add hydrogens to the "
            "prepared structure either.")
        log("          Falling back to the chain as it was read, with "
            "hydrogens added geometrically.")
        log("          The preparation (rebuilt side chains, capped termini, "
            "minimisation) is not")
        log("          in this receptor. Treat the scores as provisional and "
            "fix the preparation.")
        src = _write_chain_pdb(rec, workdir)
        p = _run_prepare_receptor(prep, src, out, add_h=True)

    if not os.path.isfile(out):
        raise RuntimeError(f"prepare_receptor failed:\n{p.stdout}\n{p.stderr}")
    if count_hydrogens(out) == 0:
        raise RuntimeError(
            "receptor.pdbqt was written with no hydrogens on it, from the "
            "prepared structure and from the raw chain alike. AGFR refuses "
            'such a receptor with "It seems that hydrogen atoms are missing '
            'in the receptor", several minutes into the grid step, so this '
            "stops here instead.\n"
            f"  written from : {src}\n"
            "Look at receptor/preparation.log for the protonation step, and "
            "at prepare_receptor's own output below.\n"
            + (p.stdout or "") + (p.stderr or ""))
    return out


# ------------------------------------------------- the prepared structure,
#                                                    everywhere downstream

def prepared_in(rundir: str) -> str | None:
    """The prepared receptor a run left in ``<rundir>/receptor/``, or None."""
    for name in ("prepared.pdb", "minimised.pdb", "protonated.pdb"):
        p = os.path.join(rundir, "receptor", name)
        if os.path.isfile(p) and os.path.getsize(p) > 0:
            return p
    return None


def rebind_prepared(rec: Receptor, cfg, prepared_pdb: str | None,
                    log=print) -> Receptor:
    """`rec` with the prepared structure's atoms in place of the input file's.

    The grid is built from the prepared structure and the poses come back in
    its frame, so everything that measures a pose against the receptor must
    read the same structure. The input file differs in two ways:

      * heavy atoms: preparation rebuilds side chains that the deposited
        structure truncates, and contacts counted against the truncated side
        chains undercount those residues;
      * coordinates: the restrained minimisation moves every atom slightly.

    The box is not recomputed: the caller already has the centre and size the
    grid was built with, and those must not move. Only the atoms change.

    Returns `rec` unchanged, with a reason logged, when the prepared file
    cannot stand in for it: no preparation ran, a hotspot is missing from it,
    or the structure moved far enough that the existing box no longer
    describes it.
    """
    if not prepared_pdb or not os.path.isfile(prepared_pdb):
        return rec
    try:
        import copy as _copy
        c = _copy.copy(cfg)
        c.receptor = prepared_pdb
        c.chain = ",".join(rec.chains) if rec.chains else cfg.chain
        got = load(c)
    except Exception as e:
        log(f"  the prepared structure could not be re-read ({e}); the input "
            f"structure is used downstream")
        return rec

    was = {(a.chain, a.resnum, a.name): a for a in rec.atoms}
    now = {(a.chain, a.resnum, a.name): a for a in got.atoms}
    lost = set(was) - set(now)
    if lost:
        # Preparation adds; it does not remove. Losing atoms means the two
        # files are not the same receptor (a different chain selection, a
        # renumbering) and every residue index downstream would be wrong.
        by_res = sorted({(c_, r) for c_, r, _n in lost})
        log(f"  the prepared structure is missing {len(lost)} atom(s) the "
            f"input has, over {len(by_res)} residue(s); it is not the same "
            f"structure, so the input is used downstream")
        return rec

    shared = [(was[k], now[k]) for k in was]
    moved = max((float(np.linalg.norm(a.xyz - b.xyz)) for a, b in shared),
                default=0.0)
    gained = len(now) - len(was)
    rec2 = Receptor(prepared_pdb, rec.chain, got.atoms,
                    list(rec.hotspots), list(rec.shell),
                    chains=list(rec.chains),
                    hotspot_keys=list(rec.hotspot_keys or []))
    rec2.shell_keys = list(getattr(rec, "shell_keys", None) or [])
    rec2.shell_radius = getattr(rec, "shell_radius", cfg.shell_radius)
    log(f"  downstream stages use the prepared receptor: "
        f"{gained} heavy atom(s) rebuilt, "
        f"largest atom shift {moved:.2f} A")
    return rec2


def write_view_pdb(path: str, prepared_pdb: str | None, rec: Receptor,
                   chain: str = "A", keep_hydrogens: bool = True) -> int:
    """The receptor as a viewer should see it. Returns the atom count.

    Written from the prepared structure when there is one, hydrogens
    included, because that is the structure the poses were docked into and
    the one whose protonation states the energies were computed with. The
    heavy-atom Receptor is the fallback, and it is the whole file that changes
    - complex_<seq>.pdb, receptor.pdb and every .mae built from it.
    """
    n = 0
    if prepared_pdb and os.path.isfile(prepared_pdb):
        want = {(" " if c == "_" else c) for c in rec.chains} or None
        with open(path, "w") as out:
            for l in open(prepared_pdb, errors="replace"):
                if l.startswith("ENDMDL"):
                    break
                if not l.startswith(("ATOM", "HETATM")):
                    continue
                if want is not None and l[21] not in want:
                    continue
                if not keep_hydrogens and is_hydrogen(l[12:16], l[76:78]):
                    continue
                n += 1
                # The record name is kept as it was: an ion or a kept water is
                # a HETATM, and calling it an ATOM tells a viewer to bond it
                # into the chain it happens to sit next to.
                out.write("%-6s%s %s%s%s%s"
                          % (l[:6].strip(), _serial(n), l[12:21], chain,
                             l[22:26], l[26:].rstrip("\n") + "\n"))
            out.write("END\n")
        if n:
            return n
    with open(path, "w") as out:
        for i, at in enumerate(rec.atoms, 1):
            out.write(pdb_line(i, at, chain))
            n = i
        out.write("END\n")
    return n
