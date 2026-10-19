"""Protein preparation, in the sense of Schrodinger's Protein Preparation Wizard.

By Aadhil Haq

pdb2pqr alone adds hydrogens and assigns titration states. A crystal
structure needs more before it can be docked against:

* one alternate location, not two overlapping copies of a side chain
* selenomethionine read as methionine
* waters removed, except any requested near the site
* other heteroatoms removed unless wanted, since an arbitrary cofactor
  cannot be parameterised here
* missing side-chain atoms rebuilt (or, if asked, left missing in the
  result), and missing residues if requested
* termini capped, so a chain break is not read as a charged terminus
* disulfides found and named, so they are not reduced by accident
* hydrogens and titration states at a stated pH
* a restrained minimisation, so the result is a structure the force field
  accepts

Every step is optional and every step is logged. A structure whose ligand
was removed, or which kept a clash, without saying so is worse than one
left alone, so the log is part of the output.

Not done: ligand protomer enumeration. Loop rebuilding is off by default,
because a modelled loop cannot be told from a measured one afterwards.
"""
from __future__ import annotations

import json
import math
import os
import shutil
import subprocess

# Residues that are not really heteroatoms: modified amino acids a structure
# uses in place of the ordinary one.
STANDIN = {"MSE": "MET", "SEP": "SER", "TPO": "THR", "PTR": "TYR",
           "CSO": "CYS", "MLY": "LYS", "HYP": "PRO"}

WATERS = {"HOH", "WAT", "H2O", "DOD", "TIP", "SOL"}

# Ions worth keeping when heteroatoms are dropped: structural, chargeable, and
# understood by the force fields involved.
KEEP_IONS = {"ZN", "MG", "CA", "MN", "FE", "CU", "NA", "K", "CL"}


class Step:
    """What one stage of preparation did, so the log can be read back."""

    def __init__(self, name, done, detail=""):
        self.name, self.done, self.detail = name, bool(done), detail

    def __str__(self):
        mark = "done" if self.done else "skipped"
        return f"  {self.name}: {mark}{('; ' + self.detail) if self.detail else ''}"


class PrepError(RuntimeError):
    """A preparation that cannot be trusted. Raised rather than logged, so the
    run stops here instead of docking against a receptor the log has already
    said is wrong."""


def _is_hydrogen(name, element=""):
    """Hydrogen test matching receptor.is_hydrogen: PDB names hydrogens
    '1HD2' as well as 'HD2', so name[0] == 'H' misses about one in twenty."""
    e = (element or "").strip()
    if e:
        return e.upper().startswith("H")
    return name.strip().lstrip("0123456789").upper().startswith("H")


def as_pdb_lines(path, workdir, log=print):
    """The structure as PDB lines, converting an mmCIF if that is what it is.

    Everything below reads fixed-column PDB: byte offsets for the chain, the
    residue number, the coordinates. An mmCIF is whitespace-delimited and read
    that way it is not merely wrong, it raises: 2,121 atoms came back as 18
    residues and then a coordinate parsed as "? -8.584".
    """
    if not str(path).lower().endswith((".cif", ".mmcif")):
        return open(path, errors="replace").read().splitlines()
    try:
        from openmm.app import PDBxFile, PDBFile
    except ImportError:
        log("  this is an mmCIF and openmm is not installed to convert it; "
            "save it as PDB first")
        return []
    out = os.path.join(workdir, "converted.pdb")
    try:
        cif = PDBxFile(path)
        with open(out, "w") as fh:
            PDBFile.writeFile(cif.topology, cif.positions, fh, keepIds=True)
    except Exception as e:
        log(f"  could not read {path} as mmCIF: {e}")
        return []
    log(f"  read as mmCIF and converted to PDB")
    return open(out, errors="replace").read().splitlines()


def looks_like_trace(lines):
    """Whether this is a backbone trace rather than a full structure.

    One atom per residue means every side chain is missing. Rebuilding them
    is a different job from preparation, it is modelling, and it takes long
    enough that starting it by accident costs a run.
    """
    res, atoms = set(), 0
    for l in lines:
        if not l.startswith("ATOM"):
            continue
        atoms += 1
        res.add((l[21], l[22:27]))
    if not res:
        return False, ""
    per = atoms / float(len(res))
    if per < 2.0:
        return True, (f"this is an alpha-carbon trace: {atoms} atoms for "
                      f"{len(res)} residues, about one each")
    return False, ""


def _atoms(lines):
    return [l for l in lines if l.startswith(("ATOM", "HETATM"))]


def _write(path, lines):
    """The lines, with exactly one END, however many the input arrived with.

    Most PDB files already end with END, and a second one is not harmless:
    pdb2pqr flushes the residue it is building when it meets END, so a second
    flush passes an empty atom list on and pdb2pqr fails with

        IndexError: list index out of range

    without naming END, the file or the residue. PDBFixer, when present,
    rewrites the file with a single END, so the problem appears only where
    PDBFixer is absent.

    ENDMDL is left alone: it closes a model, not the file.
    """
    body = [l for l in lines if l.strip() != "END"]
    with open(path, "w") as fh:
        for l in body:
            fh.write(l if l.endswith("\n") else l + "\n")
        fh.write("END\n")
    return path


def pick_altloc(lines, steps):
    """One conformer per atom: the highest occupancy, ties to the first.

    A structure with alternate locations has two copies of some side chains a
    fraction of an angstrom apart. Docking against both is docking against a
    protein that does not exist.
    """
    best, out, dropped = {}, [], 0
    for l in lines:
        if not l.startswith(("ATOM", "HETATM")):
            out.append(l)
            continue
        alt = l[16]
        if alt == " ":
            out.append(l)
            continue
        key = (l[21], l[22:27], l[12:16])
        try:
            occ = float(l[54:60])
        except ValueError:
            occ = 0.0
        if key not in best or occ > best[key][0]:
            best[key] = (occ, l)
    if best:
        keep = {id(v[1]) for v in best.values()}
        rebuilt, seen = [], set()
        for l in lines:
            if l.startswith(("ATOM", "HETATM")) and l[16] != " ":
                if id(l) in keep and l[12:16] + l[21] + l[22:27] not in seen:
                    seen.add(l[12:16] + l[21] + l[22:27])
                    rebuilt.append(l[:16] + " " + l[17:])
                else:
                    dropped += 1
            else:
                rebuilt.append(l)
        out = rebuilt
    steps.append(Step("alternate locations", dropped > 0,
                      f"kept the highest occupancy, dropped {dropped} atom(s)"
                      if dropped else "none present"))
    return out


def standardise(lines, steps):
    """Modified residues read as the residue they stand in for."""
    changed = {}
    out = []
    for l in lines:
        if l.startswith("HETATM") and l[17:20].strip() in STANDIN:
            was = l[17:20].strip()
            now = STANDIN[was]
            changed[was] = changed.get(was, 0) + 1
            # a stand-in is part of the chain, so it becomes ATOM
            l = "ATOM  " + l[6:17] + f"{now:>3}" + l[20:]
        out.append(l)
    steps.append(Step("modified residues", bool(changed),
                      ", ".join(f"{k}->{STANDIN[k]} x{v}"
                                for k, v in sorted(changed.items()))
                      or "none present"))
    return out


# Records that carry no coordinates and no connectivity: bibliography,
# database cross-references, provenance. Nothing downstream reads them, and
# their free-text fields in fixed columns are where files most often break: a
# UniProt entry name longer than the twelve columns DBREF allows spills into
# the next numeric field, which pdb2pqr then fails to read as an integer.
# SEQRES is kept: PDBFixer uses it to find missing residues, and without it
# every gap looks intentional.
METADATA = (
    "HEADER", "OBSLTE", "TITLE ", "SPLIT ", "CAVEAT", "COMPND", "SOURCE",
    "KEYWDS", "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT", "SPRSDE",
    "JRNL  ", "REMARK", "DBREF ", "DBREF1", "DBREF2", "SEQADV", "HETNAM",
    "HETSYN", "FORMUL", "SITE  ",
)


def drop_metadata(lines, steps):
    """Bibliographic and cross-reference records out.

    The file written from here is a working structure, not an archival one.
    """
    out = [l for l in lines if not l.startswith(METADATA)]
    n = len(lines) - len(out)
    steps.append(Step("header records", bool(n),
                      f"removed {n} bibliographic/cross-reference line(s) "
                      f"that no later step reads and that parsers trip over"
                      if n else "none present"))
    return out


def drop_waters(lines, steps, keep_within=0.0, near=None):
    """Waters out, except any within keep_within angstroms of the site."""
    import math
    kept, dropped = 0, 0
    out = []
    for l in lines:
        if l.startswith("HETATM") and l[17:20].strip() in WATERS:
            if keep_within > 0 and near:
                try:
                    x, y, z = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
                except ValueError:
                    dropped += 1
                    continue
                if any(math.dist((x, y, z), p) <= keep_within for p in near):
                    kept += 1
                    out.append(l)
                    continue
            dropped += 1
            continue
        out.append(l)
    steps.append(Step("waters", dropped > 0,
                      f"removed {dropped}" + (f", kept {kept} within "
                      f"{keep_within:g} A of the site" if kept else "")))
    return out


def drop_hetatms(lines, steps, keep_ions=True):
    """Ligands and cofactors out: nothing here can parameterise them."""
    dropped, kept = {}, {}
    out = []
    for l in lines:
        if l.startswith("HETATM"):
            name = l[17:20].strip()
            if keep_ions and name in KEEP_IONS:
                kept[name] = kept.get(name, 0) + 1
                out.append(l)
                continue
            dropped[name] = dropped.get(name, 0) + 1
            continue
        out.append(l)
    detail = ", ".join(f"{k} x{v}" for k, v in sorted(dropped.items())) or "none"
    if kept:
        detail += " (kept " + ", ".join(sorted(kept)) + ")"
    steps.append(Step("heteroatoms", bool(dropped), "removed " + detail))
    return out


def drop_hydrogens(lines, steps):
    """The hydrogens a structure arrives with, taken off.

    Both ways this protonates, reduce and PDB2PQR with PROPKA, build every
    hydrogen again from the heavy atoms. One already there is not neutral:
    a protonated C terminus (HXT) stops PDB2PQR outright ("Unable to debump
    biomolecule"), which left a receptor that had been prepared elsewhere
    with no hydrogens and no minimisation, docked anyway."""
    keep, n = [], 0
    for line in lines:
        # the polymer's and the waters' only: a kept ligand's hydrogens are
        # not rebuilt by either route, so they stay as given
        polymer = line.startswith("ATOM") or line[17:20].strip() in ("HOH", "WAT", "DOD")
        if polymer and line.startswith(("ATOM", "HETATM")) and _is_hydrogen(line[12:16], line[76:78]):
            n += 1
            continue
        keep.append(line)
    steps.append(Step("hydrogens as given", bool(n),
                      f"removed {n}; protonation below builds them again from the heavy atoms"
                      if n else "none present"))
    return keep


def disulfides(lines, steps, cutoff=2.5):
    """Find SG-SG pairs and say so. Naming is left to whatever consumes this.

    Reported rather than renamed: pdb2pqr works this out itself, and a CYX
    written here would be a residue name it does not know.
    """
    import math
    sg = []
    for l in lines:
        if not l.startswith("ATOM") or l[12:16].strip() != "SG":
            continue
        try:
            sg.append((l[21], l[22:27].strip(),
                       (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
        except (ValueError, IndexError):
            continue        # a line this cannot read is not a disulfide
    pairs = []
    for i in range(len(sg)):
        for j in range(i + 1, len(sg)):
            if math.dist(sg[i][2], sg[j][2]) <= cutoff:
                pairs.append(f"{sg[i][0]}{sg[i][1]}-{sg[j][0]}{sg[j][1]}")
    steps.append(Step("disulfides", bool(pairs),
                      ", ".join(pairs) if pairs else "none within "
                      f"{cutoff:g} A"))
    return pairs


# The heavy atoms of an amino acid's backbone, with the C terminal oxygen.
# Every other heavy atom an amino acid carries is its side chain.
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

# The names an amino acid can carry between here and the prepared file: the
# standard twenty, and the names pdb2pqr and the Amber force fields give to
# protonation states and disulphides.
AMINO_ACIDS = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM", "ASH", "GLH",
    "LYN", "ARN", "TYM"})


#: The longest gap rebuilt as a loop, as HOPE-MD's preparation has it
#: (model_loops_max). A longer one is a domain's worth of guessing.
MAX_LOOP = 12
#: PDBFixer's seeds for building loops, tried in turn. A seed makes the loop it
#: builds the same every time (left unseeded, 7A3R's loops came out differently
#: run to run and needed one to four relaxation attempts); a further seed is
#: tried only when no relaxation can make the first seed's loops fit.
LOOP_SEEDS = (1, 2, 3)
#: What the tool calling this names the setting that rebuilds loops, for the
#: message that says how to dock without them. Each tool sets its own.
LOOPS_SETTING = "Rebuild missing loops"


def _bonded_across(fixer, before, after, cutoff=2.0):
    """Whether two residues either side of a gap in the numbering are in fact
    joined: C to N for amino acids, O3' to P for nucleotides. Then the gap is
    a skipped number, not missing residues."""
    from openmm import unit
    pairs = (("C", "N"), ("O3'", "P"))
    for a_name, b_name in pairs:
        a = next((x for x in before.atoms() if x.name == a_name), None)
        b = next((x for x in after.atoms() if x.name == b_name), None)
        if a is not None and b is not None:
            p = fixer.positions[a.index].value_in_unit(unit.angstrom)
            q = fixer.positions[b.index].value_in_unit(unit.angstrom)
            return math.dist(tuple(p), tuple(q)) < cutoff
    return False


def fix_structure(path, out_path, steps, fill_loops=False, cap_termini=True,
                  log=print, rebuild_sidechains=True, model_missing=True,
                  built=None, max_loop=MAX_LOOP, rebuild_termini=False, seed=LOOP_SEEDS[0]):
    """Missing side-chain atoms rebuilt, and optionally missing residues.

    PDBFixer's findMissingResidues does NOT find only gaps between residues
    that are present. It aligns the chain with SEQRES, so it also returns
    every residue missing before the first one present and after the last:
    the disordered termini a crystal does not resolve. Built, those come out
    as a straight strand laid across solvent with stretched bonds; 7A3R's
    38-residue C-terminal tail spanned 177 A, more than a chain that long can
    reach at full stretch. So what is rebuilt is, as in HOPE-MD's preparation:
    internal gaps only (termini only if rebuild_termini), no longer than
    max_loop, and never a gap in the numbering whose two sides are bonded.

    Filling is off by default: a rebuilt loop is a prediction, and once
    written it is indistinguishable from measured coordinates. `built`, a
    list, receives the (chain, number, insertion code) of every residue made,
    so later steps can relax them and not hold them to positions nobody
    measured.

    With rebuild_sidechains off the same atoms are still built here when
    model_missing is on, because pdb2pqr completes any residue it is given
    whatever this step does, and the minimisation's force field has no
    template for an incomplete one. strip_modelled_side_chains takes them
    out again once both have run; the log says so instead of saying they
    were rebuilt.

    With model_missing off as well, the missing side-chain atoms are not
    built at all. That is right only when nothing later needs whole
    residues: reduce protonates a residue as it stands, so a receptor
    protonated with reduce and not minimised never needs them.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        steps.append(Step("missing atoms", False,
                          "pdbfixer is not installed, so nothing was rebuilt"))
        return None
    try:
        fixer = PDBFixer(filename=path)
        fixer.findMissingResidues()
        gaps = dict(fixer.missingResidues)
        chains = list(fixer.topology.chains())
        buildable, left = {}, []
        for key, names in gaps.items():
            ci, pos = key
            res = list(chains[ci].residues())
            cid = chains[ci].id
            if pos == 0 or pos == len(res):
                if not rebuild_termini:
                    left.append(f"{len(names)} at the {'N' if pos == 0 else 'C'} terminus of chain "
                                f"{cid}, disordered and not built")
                    continue
            elif _bonded_across(fixer, res[pos - 1], res[pos]):
                left.append(f"chain {cid} {res[pos - 1].id}-{res[pos].id}, a gap in the "
                            f"numbering only: the two residues are bonded")
                continue
            if len(names) > max_loop:
                left.append(f"{len(names)} in chain {cid}, longer than {max_loop}")
                continue
            buildable[key] = names
        fixer.missingResidues = buildable if fill_loops else {}
        present = {(r.chain.id, str(r.id).strip(), (r.insertionCode or "").strip())
                   for r in fixer.topology.residues()}
        fixer.findMissingAtoms()
        n_side = sum(len(v) for v in fixer.missingAtoms.values())
        # Of those, the ones on an amino acid's side chain. missingAtoms also
        # holds missing backbone atoms and nucleic acid atoms, which are
        # rebuilt whatever rebuild_sidechains says.
        n_sc = sum(1 for res, atoms in fixer.missingAtoms.items()
                   if res.name in AMINO_ACIDS
                   for a in atoms if a.name not in BACKBONE_ATOMS)
        n_term = sum(len(v) for v in fixer.missingTerminals.values())
        if not cap_termini:
            fixer.missingTerminals = {}
            n_term = 0
        if not rebuild_sidechains and not model_missing:
            # Only what is not an amino acid's side chain: backbone atoms
            # and nucleic acid atoms are still completed.
            for res in list(fixer.missingAtoms):
                keep = [a for a in fixer.missingAtoms[res]
                        if res.name not in AMINO_ACIDS
                        or a.name in BACKBONE_ATOMS]
                if keep:
                    fixer.missingAtoms[res] = keep
                else:
                    del fixer.missingAtoms[res]
        # PDBFixer's own minimisation on one thread of the CPU platform, so the
        # same loop comes out on every node (a GPU node's platform, or another
        # number of threads, sums the forces in another order)
        with _one_cpu_thread() as cpu:
            fixer.platform = cpu
            try:
                fixer.addMissingAtoms(seed=seed)     # the same loop every time
            except TypeError:                        # a PDBFixer too old to take one
                fixer.addMissingAtoms()
        made = [(r.chain.id, str(r.id).strip(), (r.insertionCode or "").strip())
                for r in fixer.topology.residues()
                if (r.chain.id, str(r.id).strip(), (r.insertionCode or "").strip()) not in present]
        if built is not None:
            built.extend(made)
        with open(out_path, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh,
                              keepIds=True)
    except Exception as e:
        steps.append(Step("missing atoms", False, f"pdbfixer failed: {e}"))
        return None

    if rebuild_sidechains:
        detail = f"rebuilt {n_side} side-chain atom(s)"
    elif not model_missing:
        detail = (f"left {n_sc} missing side-chain atom(s) out, as the "
                  f"structure has them; rebuilt {n_side - n_sc} other "
                  f"missing atom(s)")
    else:
        detail = (f"modelled {n_sc} missing side-chain atom(s) only so "
                  f"protonation and minimisation can run; they are taken out "
                  f"again at the end. Rebuilt {n_side - n_sc} other missing "
                  f"atom(s)")
    if n_term:
        detail += f", completed {n_term} terminal atom(s)"
    if buildable and not fill_loops:
        where = ", ".join(f"{len(v)} residue(s) after {k[1]}"
                          for k, v in list(buildable.items())[:4])
        detail += (f"; LEFT {sum(len(v) for v in buildable.values())} missing "
                   f"residue(s) alone ({where}); a rebuilt loop is a guess")
    elif buildable and fill_loops:
        detail += (f"; built {len(made)} missing residue(s) in {len(buildable)} loop(s), "
                   f"relaxed before protonation")
    if left:
        detail += "; not built: " + "; ".join(left[:6]) + (" ..." if len(left) > 6 else "")
    steps.append(Step("missing atoms", True, detail))
    return out_path


#: Measured residues either side of a rebuilt loop that may move a little
#: (on a spring) so the loop can close: the loop's stems.
FLANK = 2
#: The relaxation attempts, in order: the stems' spring (kcal/mol/A^2) and the
#: warm spell (ps). The first that leaves every peptide bond in shape is kept.
RELAX_SCHEDULE = ((10.0, 2.0), (5.0, 5.0), (2.0, 10.0), (1.0, 10.0), (0.5, 20.0))

#: What every peptide bond in or beside a rebuilt loop must measure after
#: relaxation, as ADapt's rebuilding checks its loops: the C-N bond (A), the
#: two backbone angles across it (degrees; ideal 1.33 A, 116 and 122), CA to
#: CA (A; 3.80 across a trans bond in crystal structures, a little more once
#: minimised in ff14SB) and the twist, omega (CA-C-N-CA), within `trans`
#: degrees of 180. A cis bond, within `trans` of 0, passes only before a
#: proline, with CA to CA in `ca_ca_cis`. Every rebuilt residue with a CB is
#: of the L hand: its improper N-C-CA-CB within `hand` degrees of HAND, and a
#: Thr or Ile has its CB of the natural hand as well. And no heavy atom of a
#: rebuilt residue is nearer than `contact` A to an amino acid's heavy atom
#: more than three bonds away (a hydrogen bond is 2.6 A and up).
JUNCTION = {"bond": (1.30, 1.37), "ca_c_n": (105.0, 130.0), "c_n_ca": (110.0, 135.0),
            "ca_ca": (3.70, 3.95), "ca_ca_cis": (2.70, 3.05), "trans": 30.0, "hand": 25.0,
            "contact": 2.5}
#: The improper N-C-CA-CB of an L residue, in degrees (114 to 129 in measured
#: structures), and its N-C-CA-HA (-119, within 1.4 on the 361 residues of a
#: prepared 7A3R chain A). Held on HA as well as CB, a CA stays tetrahedral:
#: held on CB alone, one in a loop often went nearly flat (improper 160 to
#: 165 degrees, N-CA-C over 120) and stayed so through every attempt.
HAND = 122.0
HAND_HA = -119.0
#: The chiral CB of Thr and Ile: (X - CB) x (Y - CB) . (CA - CB) is positive in
#: the natural form and negative in the allo form.
CHIRAL_CB = {"THR": ("OG1", "CG2"), "ILE": ("CG1", "CG2")}
#: The restraints that hold a rebuilt loop's chemistry while it relaxes and in
#: the minimisation that ends the preparation. Both are flat-bottomed, so a
#: loop inside their wells feels neither: each peptide bond in or beside it
#: trans (omega within TRANS_SLACK degrees of 180), and each rebuilt residue
#: and stem with a CB of the L hand (its impropers within HAND_SLACK degrees
#: of HAND and HAND_HA), at RESTRAINT_K kJ/mol/rad^2 outside the well.
TRANS_SLACK = 15.0
HAND_SLACK = 20.0
RESTRAINT_K = 400.0
#: How far from a moving atom (rebuilt, or a stem) the relaxation keeps the
#: structure, in nm: the 1.2 nm cutoff and 0.8 nm of room for the loop to move
#: in. Beyond it an atom cannot reach a moving one, and leaving it out of the
#: simulation changes no force on one.
NEIGHBOURHOOD_NM = 2.0


def _angle(a, b, c):
    v1 = [a[i] - b[i] for i in range(3)]
    v2 = [c[i] - b[i] for i in range(3)]
    n = math.sqrt(sum(x * x for x in v1)) * math.sqrt(sum(x * x for x in v2))
    return math.degrees(math.acos(max(-1.0, min(1.0, sum(x * y for x, y in zip(v1, v2)) / n))))


def _torsion(a, b, c, d):
    """The a-b-c-d torsion in degrees, from plain coordinates (_dihedral's
    sense: 180 is trans)."""
    import numpy as np
    return _dihedral(*(np.asarray(p, dtype=float) for p in (a, b, c, d)))


def junction_geometry(path, built):
    """Every peptide bond in or beside a rebuilt loop, and every rebuilt
    residue, checked against JUNCTION. Returns what is out of range, as text,
    one line per bond, residue or contact.

    Bond lengths alone are not enough: PDBFixer closed 7A3R's loops at
    1.31-1.36 A with the CA-C-N angle at 42 degrees, the N folded back onto
    the CA. Lengths and angles together are not enough either: they passed
    7A3R's loops with bonds cis (274-275 at -2 degrees) and twisted (194-195 at
    -95), and a cis bond that is not before a proline stays cis through any
    simulation that follows. So each bond's twist and CA to CA are checked, and
    each rebuilt residue's hand, and how near its atoms come to the rest."""
    import numpy as np
    at, order, names, heavy = {}, [], {}, []
    for line in open(path, errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            k = (line[21], line[22:26].strip(), line[26].strip())
            if not order or order[-1] != k:
                order.append(k)
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            at[k + (line[12:16].strip(),)] = xyz
            names[k] = line[17:20].strip()
            if names[k] in AMINO_ACIDS and not _is_hydrogen(line[12:16], line[76:78]):
                heavy.append((k, line[12:16].strip(), xyz))
    made = set(built)
    bad, links = [], []
    for a, b in zip(order, order[1:]):
        if a[0] != b[0] or not (a in made or b in made):
            continue
        ca1, c1, n2, ca2 = (at.get(a + ("CA",)), at.get(a + ("C",)),
                            at.get(b + ("N",)), at.get(b + ("CA",)))
        if None in (ca1, c1, n2, ca2):
            continue
        d = math.dist(c1, n2)
        if d > 2.5:
            continue                         # a chain break left open, not a junction
        links.append((a, b))
        t1, t2 = _angle(ca1, c1, n2), _angle(c1, n2, ca2)
        w, dca = _torsion(ca1, c1, n2, ca2), math.dist(ca1, ca2)
        cis = abs(w) <= JUNCTION["trans"] and names[b] == "PRO"
        lo, hi = JUNCTION["ca_ca_cis" if cis else "ca_ca"]
        if not (JUNCTION["bond"][0] <= d <= JUNCTION["bond"][1]
                and JUNCTION["ca_c_n"][0] <= t1 <= JUNCTION["ca_c_n"][1]
                and JUNCTION["c_n_ca"][0] <= t2 <= JUNCTION["c_n_ca"][1]
                and (cis or abs(w) >= 180.0 - JUNCTION["trans"])
                and lo <= dca <= hi):
            bad.append(f"{a[0]}{a[1]}-{b[1]}: C-N {d:.2f} A, CA-C-N {t1:.0f}, C-N-CA {t2:.0f}, "
                       f"omega {w:.0f}, CA-CA {dca:.2f} A")
    # The rebuilt residues, and their stems (which move with them): the hand of
    # every one, and how far a rebuilt one's CA is from the ideal
    moved = set()
    for i, k in enumerate(order):
        if k in made:
            moved.update(r for r in order[max(0, i - FLANK):i + FLANK + 1] if r[0] == k[0])
    for k in order:
        if k not in moved:
            continue
        n, c, ca, cb = (at.get(k + (x,)) for x in ("N", "C", "CA", "CB"))
        if None not in (n, c, ca, cb):
            imp = _torsion(n, c, ca, cb)
            if imp <= 0:
                bad.append(f"{k[0]}{k[1]} {names[k]}: of the D hand (improper N-C-CA-CB {imp:.0f})")
            elif k in made and abs(imp - HAND) > JUNCTION["hand"]:
                bad.append(f"{k[0]}{k[1]} {names[k]}: CA distorted (improper N-C-CA-CB {imp:.0f}, "
                           f"{HAND - JUNCTION['hand']:.0f} to {HAND + JUNCTION['hand']:.0f} allowed)")
        sub = CHIRAL_CB.get(names[k])
        if sub and None not in (ca, at.get(k + ("CB",)), at.get(k + (sub[0],)), at.get(k + (sub[1],))):
            o = np.asarray(at[k + ("CB",)])
            x, y = (np.asarray(at[k + (s,)]) - o for s in sub)
            if float(np.dot(np.cross(x, y), np.asarray(ca) - o)) <= 0:
                bad.append(f"{k[0]}{k[1]} {names[k]}: its CB of the allo hand")
    near = _closest_contact(heavy, made, links)
    if near and near[0] < JUNCTION["contact"]:
        bad.append(f"{near[1]}: {near[0]:.2f} A apart, nearer than {JUNCTION['contact']:g} A")
    return bad


def _closest_contact(heavy, made, links):
    """(distance, 'A193 CG .. A160 O') of the closest approach of a rebuilt
    residue's heavy atom to an amino acid's heavy atom more than three bonds
    away, or None. heavy: [(residue, name, xyz)]; links: the peptide bonds
    between residues, as (residue, next residue). Bonds are read from the
    distances: within a residue under 1.9 A, C to N along `links`, and SG to SG
    under 2.3 A (a disulphide)."""
    import numpy as np
    new = [i for i, h in enumerate(heavy) if h[0] in made]
    if not new or len(heavy) < 2:
        return None
    xyz = np.array([h[2] for h in heavy], dtype=float)
    close = []
    for i in new:
        d = np.sqrt(((xyz - xyz[i]) ** 2).sum(1))
        d[i] = np.inf
        for j in np.nonzero(d < 4.0)[0]:
            close.append((float(d[j]), i, int(j)))
    if not close:
        return None
    # the bonds among the atoms that came close, and their neighbours
    where = sorted({j for _d, i, j in close} | {i for _d, i, j in close})
    pocket = np.array([xyz[j] for j in where])
    region = [j for j in range(len(heavy))
              if float(np.min(((pocket - xyz[j]) ** 2).sum(1))) < 36.0]
    by_res, named, bonded = {}, {}, {j: set() for j in region}
    for j in region:
        by_res.setdefault(heavy[j][0], []).append(j)
        named[heavy[j][0], heavy[j][1]] = j

    def bond(p, q):
        bonded[p].add(q)
        bonded[q].add(p)
    for js in by_res.values():
        for p in js:
            for q in js:
                if p < q and math.dist(heavy[p][2], heavy[q][2]) < 1.9:
                    bond(p, q)
    for a, b in links:
        if (a, "C") in named and (b, "N") in named:
            bond(named[a, "C"], named[b, "N"])
    sg = [j for j in region if heavy[j][1] == "SG"]
    for p in sg:
        for q in sg:
            if p < q and math.dist(heavy[p][2], heavy[q][2]) < 2.3:
                bond(p, q)
    best = None
    for dist, i, j in sorted(close):
        if i not in bonded or j not in bonded:
            continue
        seen, edge = {i}, {i}
        for _ in range(3):
            edge = {y for x in edge for y in bonded.get(x, ())} - seen
            seen |= edge
        if j in seen:
            continue
        a, b = heavy[i], heavy[j]
        best = (dist, f"{a[0][0]}{a[0][1]}{a[0][2]} {a[1]} .. {b[0][0]}{b[0][1]}{b[0][2]} {b[1]}")
        break
    return best


class _one_cpu_thread:
    """OpenMM's CPU platform, its contexts made with one thread while this
    lasts (the default it had is put back after). Summed over several threads,
    in whatever order they finish, the same relaxation came out differently on
    two runs, and so did the hydrogens Modeller places; on one thread both are
    the same every time."""

    def __enter__(self):
        import openmm
        self.cpu = openmm.Platform.getPlatformByName("CPU")
        self.was = self.cpu.getPropertyDefaultValue("Threads")
        self.cpu.setPropertyDefaultValue("Threads", "1")
        return self.cpu

    def __exit__(self, *exc):
        self.cpu.setPropertyDefaultValue("Threads", self.was)
        return False


class _seeded_random:
    """Python's random module seeded for what runs inside, and its state put
    back after. Modeller.addHydrogens lays each new hydrogen at a random place
    and then minimises it, so unseeded the same structure came out of the
    relaxation two ways on two runs."""

    def __init__(self, seed):
        self.seed = seed

    def __enter__(self):
        import random
        self.state = random.getstate()
        random.seed(self.seed)

    def __exit__(self, *exc):
        import random
        random.setstate(self.state)
        return False


def _neighbourhood(system, xyz_nm, moving, radius_nm):
    """(a copy of `system` holding only the atoms within radius_nm of a moving
    one, the indices it kept, in order), or None when the system has a force
    this does not copy. Bonded terms, restraints and nonbonded exceptions are
    kept where all their atoms are; the masses and the nonbonded method are
    the system's own. With the radius past the cutoff, every force on a moving
    atom is what it was in the whole system (measured on 7A3R's loops: the
    same to 1 part in 10^5), and an attempt on one loop of a large protein
    simulates a fraction of it."""
    import numpy as np
    import openmm
    xyz = np.asarray(xyz_nm, dtype=float)
    near = np.zeros(len(xyz), dtype=bool)
    mv = np.asarray(sorted(moving), dtype=int)
    for chunk in np.array_split(mv, max(1, len(mv) // 256)):
        if len(chunk):
            near |= (((xyz[:, None, :] - xyz[None, chunk, :]) ** 2).sum(-1)
                     < radius_nm ** 2).any(axis=1)
    keep = [int(i) for i in np.nonzero(near)[0]]
    new = {old: i for i, old in enumerate(keep)}
    sub = openmm.System()
    for old in keep:
        sub.addParticle(system.getParticleMass(old))
    for f in system.getForces():
        if isinstance(f, openmm.HarmonicBondForce):
            g = openmm.HarmonicBondForce()
            for t in range(f.getNumBonds()):
                a, b, r0, k = f.getBondParameters(t)
                if a in new and b in new:
                    g.addBond(new[a], new[b], r0, k)
        elif isinstance(f, openmm.HarmonicAngleForce):
            g = openmm.HarmonicAngleForce()
            for t in range(f.getNumAngles()):
                a, b, c, th, k = f.getAngleParameters(t)
                if a in new and b in new and c in new:
                    g.addAngle(new[a], new[b], new[c], th, k)
        elif isinstance(f, openmm.PeriodicTorsionForce):
            g = openmm.PeriodicTorsionForce()
            for t in range(f.getNumTorsions()):
                a, b, c, d, per, ph, k = f.getTorsionParameters(t)
                if a in new and b in new and c in new and d in new:
                    g.addTorsion(new[a], new[b], new[c], new[d], per, ph, k)
        elif isinstance(f, openmm.CustomTorsionForce):
            g = openmm.CustomTorsionForce(f.getEnergyFunction())
            for t in range(f.getNumGlobalParameters()):
                g.addGlobalParameter(f.getGlobalParameterName(t),
                                     f.getGlobalParameterDefaultValue(t))
            for t in range(f.getNumPerTorsionParameters()):
                g.addPerTorsionParameter(f.getPerTorsionParameterName(t))
            for t in range(f.getNumTorsions()):
                a, b, c, d, params = f.getTorsionParameters(t)
                if a in new and b in new and c in new and d in new:
                    g.addTorsion(new[a], new[b], new[c], new[d], params)
        elif isinstance(f, openmm.NonbondedForce):
            g = openmm.NonbondedForce()
            g.setNonbondedMethod(f.getNonbondedMethod())
            g.setCutoffDistance(f.getCutoffDistance())
            g.setReactionFieldDielectric(f.getReactionFieldDielectric())
            g.setUseSwitchingFunction(f.getUseSwitchingFunction())
            g.setSwitchingDistance(f.getSwitchingDistance())
            g.setUseDispersionCorrection(f.getUseDispersionCorrection())
            for old in keep:
                g.addParticle(*f.getParticleParameters(old))
            for t in range(f.getNumExceptions()):
                a, b, q, sig, eps = f.getExceptionParameters(t)
                if a in new and b in new:
                    g.addException(new[a], new[b], q, sig, eps)
        elif isinstance(f, openmm.CustomExternalForce):
            g = openmm.CustomExternalForce(f.getEnergyFunction())
            for t in range(f.getNumGlobalParameters()):
                g.addGlobalParameter(f.getGlobalParameterName(t),
                                     f.getGlobalParameterDefaultValue(t))
            for t in range(f.getNumPerParticleParameters()):
                g.addPerParticleParameter(f.getPerParticleParameterName(t))
            for t in range(f.getNumParticles()):
                a, params = f.getParticleParameters(t)
                if a in new:
                    g.addParticle(new[a], params)
        elif isinstance(f, openmm.CMMotionRemover):
            g = openmm.CMMotionRemover(f.getFrequency())
        else:
            return None
        sub.addForce(g)
    return sub, keep


def _l_hand(topology, xyz, made):
    """xyz (nm, an n x 3 array) with every rebuilt residue's side chain set on
    the L hand before the relaxation, as ADapt sets its loops before refining
    them; returns how many residues were set right. PDBFixer's soft dynamics
    leave some CA of the D hand (seven of the sixteen it builds in 7A3R chain
    A), and a hydrogen placed on such a CA takes the place CB belongs in. So a
    residue whose improper N-C-CA-CB is out of JUNCTION's range has its side
    chain turned rigidly about CA until CB lies where an L residue's would,
    from N, CA and C (the constants place CB with ideal angles; not for Pro,
    whose ring a turn would open at N), and a Thr or Ile whose CB has the allo
    hand has two of its substituents swapped (CD1 going with CG1)."""
    import numpy as np

    def key(r):
        return (r.chain.id, str(r.id).strip(), (r.insertionCode or "").strip())
    righted = 0
    for r in topology.residues():
        if key(r) not in made:
            continue
        at = {a.name: a.index for a in r.atoms()
              if a.element is not None and a.element.symbol != "H"}
        if not all(x in at for x in ("N", "CA", "C", "CB")):
            continue
        n, ca, c, cb = (xyz[at[x]].copy() for x in ("N", "CA", "C", "CB"))
        imp = _dihedral(n, c, ca, cb)
        turned = False
        if r.name != "PRO" and (imp <= 0 or abs(imp - HAND) > JUNCTION["hand"]):
            b, d = (ca - n) * 10.0, (c - ca) * 10.0       # in A: the constants hold there, not in nm
            want = -0.58273431 * np.cross(b, d) + 0.56802827 * b - 0.54067466 * d
            have = cb - ca
            axis = np.cross(have, want)
            sin, cos = float(np.linalg.norm(axis)), float(np.dot(have, want))
            if sin < 1e-9 * float(np.linalg.norm(have) * np.linalg.norm(want)):
                axis = np.cross(have, b)         # CB straight opposite: any axis across it
                sin = 0.0
            axis = axis / float(np.linalg.norm(axis))
            angle = math.atan2(sin, cos)
            if abs(angle) > 1e-6:
                kx = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]],
                               [-axis[1], axis[0], 0]])
                turn = np.eye(3) + math.sin(angle) * kx + (1 - math.cos(angle)) * kx.dot(kx)
                side = [i for name, i in at.items() if name not in BACKBONE_ATOMS]
                xyz[side] = (xyz[side] - ca).dot(turn.T) + ca
                turned = True
        sub = CHIRAL_CB.get(r.name)
        if sub and all(s in at for s in sub):
            o = xyz[at["CB"]]
            x, y = (xyz[at[s]] - o for s in sub)
            if float(np.dot(np.cross(x, y), xyz[at["CA"]] - o)) <= 0:
                px, py = xyz[at[sub[0]]].copy(), xyz[at[sub[1]]].copy()
                xyz[at[sub[0]]], xyz[at[sub[1]]] = py, px
                if "CD1" in at:
                    xyz[at["CD1"]] = xyz[at["CD1"]] + py - px
                turned = True
        righted += turned
    return righted


def _ha_opposite(topology, xyz, made):
    """xyz (nm) with the HA of every rebuilt residue with a CB put opposite its
    CA's three other substituents, 0.109 nm from CA, where the L hand has it.
    Modeller lays a new hydrogen at random and minimises it, and on a CA that
    PDBFixer left nearly flat it can settle on CB's side. With HA held as well
    (HAND_HA), 7A3R chain A's loops passed on the second attempt with HA placed
    here first, and on the third without. Returns how many were placed."""
    import numpy as np

    def key(r):
        return (r.chain.id, str(r.id).strip(), (r.insertionCode or "").strip())
    placed = 0
    for r in topology.residues():
        if key(r) not in made:
            continue
        at = {a.name: a.index for a in r.atoms()}
        if not all(x in at for x in ("N", "CA", "C", "CB", "HA")):
            continue
        ca = xyz[at["CA"]].copy()
        u = sum((xyz[at[x]] - ca) / np.linalg.norm(xyz[at[x]] - ca) for x in ("N", "C", "CB"))
        if np.linalg.norm(u) > 1e-6:
            xyz[at["HA"]] = ca - 0.109 * u / np.linalg.norm(u)
            placed += 1
    return placed


def _hold_chemistry(system, topology, made, stems=()):
    """The restraints that keep a rebuilt loop's chemistry while it relaxes,
    added to `system`: every peptide bond in or beside it trans, and every
    rebuilt residue with a CB, and every stem, of the L hand, its CB and its
    HA each held to their side (TRANS_SLACK, HAND_SLACK, HAND, HAND_HA).
    PDBFixer's soft force field leaves loops with bonds anywhere
    from cis to trans, and Amber's force field keeps whatever it is given: its
    barrier holds a cis bond cis. Unheld, a stem on its spring was turned to
    the D hand (7A3R's A:275). Returns (bonds, residues) held."""
    import openmm

    def key(r):
        return (r.chain.id, str(r.id).strip(), (r.insertionCode or "").strip())

    trans = openmm.CustomTorsionForce(
        "0.5*k_trans*max(0, d-w_trans)^2; d=min(dt, 6.283185307179586-dt); "
        "dt=abs(theta-3.141592653589793)")
    trans.addGlobalParameter("k_trans", RESTRAINT_K)
    trans.addGlobalParameter("w_trans", math.radians(TRANS_SLACK))
    hand = openmm.CustomTorsionForce(
        "0.5*k_hand*max(0, abs(d)-w_hand)^2; d=atan2(sin(theta-t0), cos(theta-t0))")
    hand.addGlobalParameter("k_hand", RESTRAINT_K)
    hand.addGlobalParameter("w_hand", math.radians(HAND_SLACK))
    hand.addPerTorsionParameter("t0")
    bonded = {frozenset((a.index, b.index)) for a, b in topology.bonds()}
    res = list(topology.residues())

    def atom(r, name):
        return next((a.index for a in r.atoms() if a.name == name), None)
    for a, b in zip(res, res[1:]):
        if a.chain.index != b.chain.index or not (key(a) in made or key(b) in made):
            continue
        idx = (atom(a, "CA"), atom(a, "C"), atom(b, "N"), atom(b, "CA"))
        if None not in idx and frozenset((idx[1], idx[2])) in bonded:
            trans.addTorsion(*idx, [])
    held = 0
    for r in res:
        if key(r) in made or key(r) in stems:
            idx = (atom(r, "N"), atom(r, "C"), atom(r, "CA"), atom(r, "CB"))
            if None not in idx:
                hand.addTorsion(*idx, [math.radians(HAND)])
                held += 1
                ha = atom(r, "HA")
                if ha is not None:
                    hand.addTorsion(idx[0], idx[1], idx[2], ha, [math.radians(HAND_HA)])
    system.addForce(trans)
    system.addForce(hand)
    return trans.getNumTorsions(), held


def relax_built(path, out_path, steps, built, log=print):
    """The residues fix_structure made, relaxed before anything protonates them.

    Every measured atom is held exactly where it was measured (mass zero, so
    the minimiser cannot move it); only built residues move, with hydrogens
    from OpenMM's templates so the force field covers them. Done before
    protonation because reduce cannot place an amide hydrogen on a junction
    PDBFixer closed badly, and the minimisation that would have repaired it
    then fails for the want of that hydrogen. While it relaxes, every peptide
    bond in or beside the loop is held trans and every rebuilt residue of the
    L hand (_hold_chemistry); afterwards every one of those bonds and residues
    is checked (junction_geometry), and a loop that still does not fit stops
    the preparation.

    The same input relaxes the same way every time: the hydrogens are placed
    from a seed, each attempt's integrator is seeded before its context is
    made (OpenMM reads the seed then and never again), and each attempt runs
    on one thread, with only the atoms the moving ones can feel
    (_neighbourhood). On several threads the CPU platform sums its forces in a
    different order each run, and the dynamics carried the difference into
    which loop came out."""
    if not built:
        return path
    import openmm
    from openmm import app, unit
    made = set(built)

    def key(r):
        return (r.chain.id, str(r.id).strip(), (r.insertionCode or "").strip())

    try:
        pdb = app.PDBFile(path)
        ff = app.ForceField("amber14-all.xml")
        m = app.Modeller(pdb.topology, pdb.positions)
        other = [r for r in m.topology.residues() if r.name not in AMINO_ACIDS]
        if other:
            m.delete(other)              # relaxation only; the file keeps them
        import numpy as np
        xyz = np.array(m.positions.value_in_unit(unit.nanometer))
        righted = _l_hand(m.topology, xyz, made)
        m.positions = [openmm.Vec3(*p) for p in xyz] * unit.nanometer
        # Modeller lays each new hydrogen at random before minimising it
        with _seeded_random(LOOP_SEEDS[0]), _one_cpu_thread() as cpu:
            m.addHydrogens(ff, pH=7.0, platform=cpu)
        xyz = np.array(m.positions.value_in_unit(unit.nanometer))
        _ha_opposite(m.topology, xyz, made)
        m.positions = [openmm.Vec3(*p) for p in xyz] * unit.nanometer
        system = ff.createSystem(m.topology, nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=1.2 * unit.nanometer, constraints=None)
        # A loop held rigidly at both ends often cannot close: 7A3R's stayed
        # at C-N-CA 92 degrees. As loop modelling does, the stems (FLANK
        # measured residues either side) may give a little, on a spring to
        # where they were measured; everything further away does not move.
        res = list(m.topology.residues())
        idx = {key(r): i for i, r in enumerate(res)}
        stems = set()
        for k in made:
            i = idx.get(k)
            if i is None:
                continue
            for j in range(i - FLANK, i + FLANK + 1):
                if 0 <= j < len(res) and res[j].chain.id == res[i].chain.id and key(res[j]) not in made:
                    stems.add(key(res[j]))
        spring = openmm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        spring.addGlobalParameter("k", 10.0 * 418.4)        # 10 kcal/mol/A^2 in kJ/mol/nm^2
        for p in ("x0", "y0", "z0"):
            spring.addPerParticleParameter(p)
        for a in m.topology.atoms():
            k = key(a.residue)
            if k in made:
                continue
            if k in stems and a.element is not None and a.element.symbol != "H":
                spring.addParticle(a.index, m.positions[a.index].value_in_unit(unit.nanometer))
            elif k not in stems:
                system.setParticleMass(a.index, 0.0)
        system.addForce(spring)
        held = _hold_chemistry(system, m.topology, made, stems)
        start = {a.index: m.positions[a.index].value_in_unit(unit.angstrom)
                 for a in m.topology.atoms() if key(a.residue) in stems}
        moving = [a.index for a in m.topology.atoms()
                  if key(a.residue) in made or key(a.residue) in stems]
        whole = m.positions.value_in_unit(unit.nanometer)
        cut = _neighbourhood(system, whole, moving, NEIGHBOURHOOD_NM)
        sub, keep = cut if cut else (system, list(range(system.getNumParticles())))
        pos0 = [tuple(whole[i]) for i in keep]
        xml = openmm.XmlSerializer.serialize(sub)
        ctx = openmm.Context(sub, openmm.VerletIntegrator(0.001),
                             openmm.Platform.getPlatformByName("CPU"), {"Threads": "1"})
        ctx.setPositions([openmm.Vec3(*p) for p in pos0] * unit.nanometer)
        e0 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilocalories_per_mole)
        del ctx
    except Exception as e:  # noqa: BLE001 - any failure here leaves invented atoms unrelaxed
        raise PrepError(f"the {len(made)} rebuilt residue(s) could not be relaxed ({str(e)[:200]}). "
                        f"Untick '{LOOPS_SETTING}' to dock the structure as measured.") from e

    def write(pos):
        new = {}
        for a in m.topology.atoms():
            if (key(a.residue) in made or key(a.residue) in stems) \
                    and a.element is not None and a.element.symbol != "H":
                new[key(a.residue) + (a.name,)] = pos[a.index]
        out, moved = [], 0
        for line in open(path, errors="replace"):
            if line.startswith(("ATOM", "HETATM")):
                k = (line[21], line[22:26].strip(), line[26].strip(), line[12:16].strip())
                if k in new:
                    x, y, z = new[k]
                    line = line[:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + line[54:]
                    moved += 1
            out.append(line)
        with open(out_path, "w") as fh:
            fh.writelines(out)
        return moved

    # Each attempt starts again from PDBFixer's loop: minimise, a warm spell to
    # leave its closure, minimise again. Later attempts loosen the stems and
    # stay warm longer. They are independent, so they run side by side, one
    # process and one thread each, and the first in RELAX_SCHEDULE's order
    # whose every bond and residue passes is kept: the answer a run of them one
    # after another would give, in the time of the slowest one it needs.
    tried = []
    whole_a = [tuple(p) for p in m.positions.value_in_unit(unit.angstrom)]
    jobs = [(xml, pos0, n, k_spring, warm_ps)
            for n, (k_spring, warm_ps) in enumerate(RELAX_SCHEDULE, 1)]
    folder = os.path.join(os.path.dirname(os.path.abspath(out_path)), "relax_attempts")
    with _attempts(folder, len(jobs)) as run:
        failed = []
        for n, k_spring, e1, got in run(jobs):
            if isinstance(got, str):
                # an attempt that blew up (one on 7A3R chain A ended "Particle
                # coordinate is NaN") is one that did not pass: the next goes on
                failed.append(f"attempt {n}: {got[:200]}")
                continue
            pos = list(whole_a)
            for j, i in enumerate(keep):
                pos[i] = tuple(got[j])
            stem_rmsd = math.sqrt(sum(math.dist(tuple(pos[i]), tuple(p)) ** 2
                                      for i, p in start.items())
                                  / max(1, len(start))) if start else 0.0
            moved = write(pos)
            bad = junction_geometry(out_path, built)
            tried.append((len(bad), bad, k_spring))
            if not bad:
                break
        else:
            if not tried:
                raise PrepError(f"the {len(made)} rebuilt residue(s) could not be relaxed "
                                f"({failed[0]}). Untick '{LOOPS_SETTING}' to dock the structure "
                                f"as measured.")
            fewest = min(tried, key=lambda t: t[0])
            raise PrepError(f"{fewest[0]} fault(s) in the peptide bonds and residues of the rebuilt "
                            f"loops remain after {len(RELAX_SCHEDULE)} relaxation attempts "
                            f"({'; '.join(fewest[1][:4])}). Untick '{LOOPS_SETTING}' to dock the "
                            f"structure as measured.")
    steps.append(Step("rebuilt loops relaxed", True,
                      f"{len(made)} residue(s) free; {len(stems)} stem residue(s) either side on a "
                      f"{k_spring:g} kcal/mol/A^2 spring (attempt {n} of {len(RELAX_SCHEDULE)}), "
                      f"moved {stem_rmsd:.2f} A RMSD; every other measured atom held; "
                      f"{righted} rebuilt residue(s) set on the L hand first; {held[0]} peptide "
                      f"bond(s) held trans and {held[1]} residue(s) of the L hand while it "
                      f"relaxed; {moved} heavy atoms moved; energy of the {len(keep)} atoms "
                      f"simulated {e0:,.0f} -> {e1:,.0f} kcal/mol; every peptide bond in or "
                      f"beside the loops trans, C-N {JUNCTION['bond'][0]}-{JUNCTION['bond'][1]} A, "
                      f"CA-CA {JUNCTION['ca_ca'][0]}-{JUNCTION['ca_ca'][1]} A and backbone "
                      f"angles in range, every rebuilt residue of the L hand, none nearer than "
                      f"{JUNCTION['contact']:g} A to another atom"))
    return out_path


def _relax_attempt(job):
    """One attempt of RELAX_SCHEDULE: (attempt, spring, energy in kcal/mol,
    the positions in A of the system's atoms, or the error as text). Minimise,
    a warm spell to leave PDBFixer's closure, minimise again, on one thread.
    The integrator is seeded before its context is made, which is when OpenMM
    reads the seed."""
    xml, pos0, n, k_spring, warm_ps = job
    try:
        import openmm
        from openmm import unit
        system = openmm.XmlSerializer.deserialize(xml)
        integ = openmm.LangevinMiddleIntegrator(300 * unit.kelvin, 5.0 / unit.picosecond,
                                                0.0005 * unit.picoseconds)
        integ.setRandomNumberSeed(n)
        ctx = openmm.Context(system, integ, openmm.Platform.getPlatformByName("CPU"),
                             {"Threads": "1"})
        ctx.setPositions([openmm.Vec3(*p) for p in pos0] * unit.nanometer)
        ctx.setParameter("k", k_spring * 418.4)
        openmm.LocalEnergyMinimizer.minimize(ctx, 1.0, 5000)
        ctx.setVelocitiesToTemperature(300 * unit.kelvin, n)
        integ.step(int(round(warm_ps / 0.0005)))
        openmm.LocalEnergyMinimizer.minimize(ctx, 0.5, 10000)
        state = ctx.getState(getEnergy=True, getPositions=True)
        return (n, k_spring, state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole),
                [tuple(p) for p in state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)])
    except Exception as e:  # noqa: BLE001 - said by the caller, which stops the preparation
        return (n, k_spring, 0.0, f"{type(e).__name__}: {e}")


def _relax_main(argv):
    """One attempt run as a process of its own (_attempts): argv is the folder
    holding system.xml and start.json, the attempt's number, its spring and its
    warm spell; the result goes to attempt_<n>.json in that folder."""
    folder, n, k_spring, warm_ps = argv[0], int(argv[1]), float(argv[2]), float(argv[3])
    with open(os.path.join(folder, "system.xml")) as fh:
        xml = fh.read()
    with open(os.path.join(folder, "start.json")) as fh:
        pos0 = json.load(fh)
    n, k_spring, e1, got = _relax_attempt((xml, pos0, n, k_spring, warm_ps))
    with open(os.path.join(folder, f"attempt_{n}.json"), "w") as fh:
        json.dump({"e1": e1, "pos": got}, fh)


#: How a process of its own runs one attempt: this file, loaded by its path,
#: since each tool carries it in a package of its own name.
_ATTEMPT = ("import importlib.util, sys; "
            "s = importlib.util.spec_from_file_location('_relax_attempt', sys.argv[1]); "
            "m = importlib.util.module_from_spec(s); s.loader.exec_module(m); "
            "m._relax_main(sys.argv[2:])")


class _attempts:
    """Runs _relax_attempt over the jobs, as many at once as this process has
    CPUs, each in a process of its own on one thread, and yields the results in
    the jobs' order; leaving the block stops whatever is still running. The
    processes are started afresh, not forked: OpenMM's threads in this one
    would be copied half-held, and a script that called the preparation
    without a main guard would be run again by multiprocessing's spawn. With
    one CPU, or where a process cannot start, the attempts run here in turn,
    which gives the same results."""

    def __init__(self, folder, n):
        self.folder, self.n, self.procs = folder, n, {}
        try:
            self.cpus = len(os.sched_getaffinity(0))
        except (AttributeError, OSError):
            self.cpus = os.cpu_count() or 1

    def __enter__(self):
        return self.run

    def _start(self, job):
        import sys
        xml, pos0, n, k_spring, warm_ps = job
        env = dict(os.environ, OPENMM_CPU_THREADS="1")
        self.procs[n] = subprocess.Popen(
            [sys.executable, "-c", _ATTEMPT, os.path.abspath(__file__), self.folder, str(n),
             repr(k_spring), repr(warm_ps)],
            env=env, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)

    def run(self, jobs):
        if self.cpus <= 1 or self.n <= 1:
            for job in jobs:
                yield _relax_attempt(job)
            return
        try:
            os.makedirs(self.folder, exist_ok=True)
            with open(os.path.join(self.folder, "system.xml"), "w") as fh:
                fh.write(jobs[0][0])
            with open(os.path.join(self.folder, "start.json"), "w") as fh:
                json.dump(jobs[0][1], fh)
            width = min(self.cpus, len(jobs))
            for job in jobs[:width]:
                self._start(job)
        except OSError:
            for job in jobs:
                yield _relax_attempt(job)
            return
        for k, job in enumerate(jobs):
            n, k_spring = job[2], job[3]
            ahead = k + width - 1
            if ahead < len(jobs) and jobs[ahead][2] not in self.procs:
                self._start(jobs[ahead])
            proc = self.procs[n]
            _, err = proc.communicate()
            got = os.path.join(self.folder, f"attempt_{n}.json")
            if proc.returncode != 0 or not os.path.isfile(got):
                yield (n, k_spring, 0.0,
                       f"attempt {n} ended without a result: {(err or '').strip()[-300:]}")
                continue
            with open(got) as fh:
                res = json.load(fh)
            yield (n, k_spring, res["e1"], res["pos"])

    def __exit__(self, *exc):
        for proc in self.procs.values():
            if proc.poll() is None:
                proc.kill()
                proc.wait()
        shutil.rmtree(self.folder, ignore_errors=True)
        return False


def strip_modelled_side_chains(measured, path, out_path, steps):
    """Side chains left as the structure had them.

    `measured` is the structure before anything was rebuilt (tidied.pdb) and
    `path` the structure after protonation and minimisation. On every amino
    acid that `measured` has, a side-chain heavy atom it lacks is taken out,
    with the hydrogens on it. The backbone, the C terminal oxygen, residues
    the structure did not have at all (a rebuilt loop, a cap) and anything
    that is not an amino acid are left as they are.

    An atom counts as measured if `measured` has one of that name in that
    residue, or has a heavy atom of that residue within 0.9 A of it. The
    second covers a measured atom the force field renamed: the restrained
    minimisation moves a heavy atom a fraction of that, and a modelled atom
    sits at least a bond length, 1.4 A, from every measured one.
    """
    import math

    def heavy(l):
        return not _is_hydrogen(l[12:16], l[76:78])

    def key(l):
        return (l[21], l[22:27])

    def xyz(l):
        return (float(l[30:38]), float(l[38:46]), float(l[46:54]))

    had = {}
    for l in open(measured, errors="replace"):
        if (l.startswith(("ATOM", "HETATM")) and heavy(l)
                and l[17:20].strip() in AMINO_ACIDS):
            try:
                had.setdefault(key(l), []).append((l[12:16].strip(), xyz(l)))
            except ValueError:
                continue

    lines = open(path, errors="replace").read().splitlines()
    by_res = {}
    for i, l in enumerate(lines):
        if l.startswith(("ATOM", "HETATM")):
            by_res.setdefault(key(l), []).append(i)

    drop, where, n_heavy, n_h = set(), [], 0, 0
    for k, idx in by_res.items():
        resname = lines[idx[0]][17:20].strip()
        if k not in had or resname not in AMINO_ACIDS:
            continue
        names = {n for n, _p in had[k]}
        places = [p for _n, p in had[k]]
        try:
            atoms = [(i, xyz(lines[i]), heavy(lines[i])) for i in idx]
        except ValueError:
            continue                    # a line this cannot read; left alone
        gone = set()
        for i, p, is_heavy in atoms:
            name = lines[i][12:16].strip()
            if (not is_heavy or name in BACKBONE_ATOMS or name in names
                    or any(math.dist(p, q) < 0.9 for q in places)):
                continue
            gone.add(i)
        if not gone:
            continue
        # A hydrogen goes with the heavy atom it sits on, its nearest.
        heavies = [(i, p) for i, p, is_heavy in atoms if is_heavy]
        for i, p, is_heavy in atoms:
            if not is_heavy and min(heavies,
                                    key=lambda t: math.dist(p, t[1]))[0] in gone:
                drop.add(i)
                n_h += 1
        drop |= gone
        n_heavy += len(gone)
        where.append("%s:%s %s" % (k[0], k[1].strip(), resname))

    if not drop:
        steps.append(Step("side chains", True, "left as the structure has "
                          "them; none had atoms missing"))
        return None

    # CONECT records naming an atom that is gone would point at nothing.
    serials = set()
    for i in drop:
        try:
            serials.add(int(lines[i][6:11]))
        except ValueError:
            pass

    def stale(l):
        if not l.startswith("CONECT"):
            return False
        for a in range(6, 31, 5):
            try:
                if int(l[a:a + 5]) in serials:
                    return True
            except ValueError:
                continue
        return False

    _write(out_path, [l for i, l in enumerate(lines)
                      if i not in drop and not stale(l)])
    steps.append(Step(
        "side chains", True,
        f"left as the structure has them: {n_heavy} side-chain atom(s) it "
        f"lacks, and the {n_h} hydrogen(s) on them, were modelled only for "
        f"protonation and minimisation and are not in the prepared receptor; "
        f"{len(where)} residue(s): {', '.join(where[:8])}"
        + (", ..." if len(where) > 8 else "")))
    return out_path


# OpenMM and Amber do not name a cap's atoms the same way, and the file goes
# to tleap next. OpenMM writes the methyl of an NME as C with H1, H2 and H3
# on it; Amber calls them CH3 and HH31 to HH33, and tleap answers an unknown
# atom name by dropping the atom and rebuilding it. Renaming here keeps the
# hydrogens that were just minimised instead of having them guessed again.
CAP_NAMES = {"ACE": {"H1": "HH31", "H2": "HH32", "H3": "HH33"},
             "NME": {"C": "CH3", "H1": "HH31", "H2": "HH32", "H3": "HH33"},
             "NMA": {"C": "CH3", "H1": "HH31", "H2": "HH32", "H3": "HH33"}}


def amber_cap_names(path):
    """Cap atoms renamed from OpenMM's spelling to Amber's, in place."""
    out, changed = [], 0
    for l in open(path, errors="replace"):
        if l[:6].strip() in ("ATOM", "HETATM"):
            want = CAP_NAMES.get(l[17:20].strip(), {}).get(l[12:16].strip())
            if want:
                # A FOUR-CHARACTER NAME FILLS THE FIELD. Padding it like a
                # short one and then trimming to four wrote HH31 as HH3.
                field = want if len(want) == 4 else (" " + want).ljust(4)
                l = l[:12] + field + l[16:]
                changed += 1
        out.append(l)
    if changed:
        with open(path, "w") as fh:
            fh.writelines(out)
    return changed


def minimise(path, out_path, steps, restraint=10.0, max_iterations=500,
             log=print, free=()):
    """A restrained minimisation, the way the wizard's last step is.

    Heavy atoms are held with a spring rather than frozen, so hydrogens and
    clashes relax while the structure that was measured stays where it was
    measured. Reported as the RMSD it moved: if that is large, something was
    wrong before minimising and is worth looking at rather than accepting.
    """
    try:
        import openmm
        from openmm import app, unit
    except ImportError:
        steps.append(Step("minimisation", False,
                          "openmm is not installed, so nothing was minimised"))
        return None
    # The same chain-break problem affects OpenMM. It bonds residue i's C
    # to residue i+1's N for every consecutive pair within a chain, by
    # name rather than by distance, so a structure with an unmodelled
    # loop gets a long bond across the gap. The residues on either side,
    # which pdb2pqr has correctly given terminal atoms, then match no
    # template, and the error names a residue without mentioning the gap.
    # Naming the fragments as separate chains avoids this, and the
    # original chain letters are restored before anything downstream
    # reads the file.
    mapping = None
    gaps = backbone_gaps(path)
    if gaps:
        split = os.path.join(os.path.dirname(out_path), "split_for_min.pdb")
        mapping = split_at_gaps(path, split, gaps)
        if mapping:
            path = split
            steps.append(Step(
                "minimisation setup", True,
                "%d backbone break(s) named as separate chains, so the force "
                "field is not asked to bond across them" % len(gaps)))
    try:
        pdb = app.PDBFile(path)
        # tip3p IS WHERE THE IONS ARE. amber14-all.xml carries protein,
        # nucleic acid and lipid templates and no monovalent ions, so a
        # receptor keeping the sodium a crystal came with stopped the
        # minimisation with "No template found for residue 104 (NA). The
        # residue contains Na atoms, which are not supported by any template
        # in the force field" - and keep_ions is on by default, so that is
        # most structures. Nothing else changes: with implicit solvent there
        # are no waters to parameterise unless the caller kept some, and
        # then TIP3P is the right answer for them too.
        ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml",
                            "implicit/gbn2.xml")
        # The caps need their hydrogens, and the force field places
        # them: an NME is -NH-CH3, and placing those by hand would be
        # four further guesses.
        topology, positions = pdb.topology, pdb.positions
        if any(r.name in CAP_RESIDUES for r in topology.residues()):
            _m = app.Modeller(topology, positions)
            # seeded and on one thread: Modeller lays each new hydrogen at
            # random before minimising it
            with _seeded_random(LOOP_SEEDS[0]), _one_cpu_thread() as cpu:
                _m.addHydrogens(ff, platform=cpu)
            topology, positions = _m.topology, _m.positions
        # The caps are placed on the right bonds at the right lengths and
        # nothing more, so they are the one thing this does NOT hold still:
        # a restraint on a guessed position is a guess that cannot correct
        # itself. Everything measured stays where it was measured.
        try:
            # A cutoff, not NoCutoff. Generalised-Born with every pair
            # considered is quadratic in atom count, and measured it was 391
            # of the 393 seconds a preparation took: the rest of the work was
            # 2.3. At 1.2 nm the pairs that matter are still there and the
            # ones beyond it were contributing rounding.
            system = ff.createSystem(
                topology, nonbondedMethod=app.CutoffNonPeriodic,
                nonbondedCutoff=1.2 * unit.nanometer,
                constraints=app.HBonds)
        except ValueError as e:
            steps.append(Step("minimisation", False,
                              f"the force field does not cover this structure "
                              f"({str(e)[:120]}); it was left as it was"))
            return None

        # A spring per heavy atom, to its measured position.
        force = openmm.CustomExternalForce(
            "k*periodicdistance(x, y, z, x0, y0, z0)^2")
        force.addGlobalParameter(
            "k", restraint * unit.kilocalories_per_mole / unit.angstroms ** 2)
        for p in ("x0", "y0", "z0"):
            force.addPerParticleParameter(p)
        heavy, free_caps, free_built = 0, 0, 0
        # Built residues are guesses too, so by the rule above they are not
        # held either: holding them kept 6B88's rebuilt loops where PDBFixer
        # put them (6.4 million kcal/mol at the start, 0.15 A moved).
        guessed = set(free)
        for atom in topology.atoms():
            if atom.element is None or atom.element.symbol == "H":
                continue
            if atom.residue.name in CAP_RESIDUES:
                free_caps += 1          # placed by geometry; let it relax
                continue
            r = atom.residue
            if (r.chain.id, str(r.id).strip(), (r.insertionCode or "").strip()) in guessed:
                free_built += 1         # built by fix_structure; let it relax
                continue
            force.addParticle(atom.index, positions[atom.index])
            heavy += 1
        system.addForce(force)
        # Free is not unheld: the rebuilt loops keep the chemistry their
        # relaxation gave them, by the restraints it used (flat-bottomed, so a
        # loop inside their wells feels none of them). Released here, a
        # strained loop twisted a peptide bond 40 degrees from trans (7A3R's
        # 274-275).
        if guessed:
            _hold_chemistry(system, topology, guessed)

        integrator = openmm.LangevinIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        # One thread, so the same structure minimises to the same coordinates
        # every time: on several, the CPU platform sums its forces in whatever
        # order the threads finish.
        sim = app.Simulation(topology, system, integrator,
                             openmm.Platform.getPlatformByName("CPU"), {"Threads": "1"})
        sim.context.setPositions(positions)
        before = sim.context.getState(getEnergy=True).getPotentialEnergy()
        sim.minimizeEnergy(maxIterations=max_iterations)
        state = sim.context.getState(getPositions=True, getEnergy=True)
        after = state.getPotentialEnergy()

        import math
        p0 = positions.value_in_unit(unit.angstrom)
        p1 = state.getPositions().value_in_unit(unit.angstrom)
        moved = [sum((p1[i][k] - p0[i][k]) ** 2 for k in range(3))
                 for i, a in enumerate(topology.atoms())
                 if a.element is not None and a.element.symbol != "H"]
        rmsd = math.sqrt(sum(moved) / len(moved)) if moved else 0.0

        with open(out_path, "w") as fh:
            app.PDBFile.writeFile(sim.topology, state.getPositions(), fh,
                                  keepIds=True)
        amber_cap_names(out_path)
        if mapping:
            restore_chains(out_path, mapping)
    except Exception as e:
        steps.append(Step("minimisation", False,
                          f"failed: {type(e).__name__}: {str(e)[:140]}"))
        return None

    kc = unit.kilocalories_per_mole
    steps.append(Step(
        "minimisation", True,
        f"{heavy} heavy atoms restrained at {restraint:g} kcal/mol/A^2"
        + (f" ({free_caps} cap atom(s) left free to settle)" if free_caps
           else "")
        + (f" ({free_built} atom(s) of rebuilt loops left free, as no one measured them)"
           if free_built else "") + "; "
        f"energy {before.value_in_unit(kc):,.0f} -> {after.value_in_unit(kc):,.0f} "
        f"kcal/mol; heavy atoms moved {rmsd:.2f} A RMSD"))
    if rmsd > 1.0:
        steps.append(Step(
            "minimisation warning", True,
            f"{rmsd:.2f} A is a long way for a restrained minimisation. The "
            f"structure disagreed with the force field somewhere, which is worth "
            f"looking at before trusting a docking against it."))
    return out_path


def find_tool(*names):
    """An executable, from PATH or from this interpreter's own bin directory.

    shutil.which searches only PATH, so a job that runs the environment's
    python without putting its bin on PATH would find no pdb2pqr and skip
    protonation. sys.executable is <env>/bin/python, so its directory holds
    the rest of the environment's programs.
    """
    import sys
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    here = os.path.dirname(os.path.abspath(sys.executable))
    for name in names:
        candidate = os.path.join(here, name)
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


def _pdb2pqr(exe, path, out_path, ph, forcefield, propka, nodebump=False):
    """One pdb2pqr attempt. Returns (ok, why) and never raises.

    `propka` off drops the pH-dependent titration states and keeps the force
    field's standard ones, which loses fidelity and is reported.
    """
    pqr = os.path.splitext(out_path)[0] + ".pqr"
    cmd = [exe, f"--ff={forcefield}", "--keep-chain",
           f"--pdb-output={out_path}", path, pqr]
    if nodebump:
        # The last resort for a structure with a gap that could not
        # be split. Debumping resolves steric clashes, so turning it
        # off loses something, but less than having no hydrogens at
        # all.
        cmd.insert(1, "--nodebump")
    if propka:
        cmd[2:2] = [f"--with-ph={ph}", "--titration-state-method=propka"]
    if os.path.exists(out_path):
        os.remove(out_path)
    p = subprocess.run(cmd, capture_output=True, text=True)
    if os.path.isfile(out_path):
        return True, ""
    # The reason is in pdb2pqr's own CRITICAL lines, not in the traceback
    # that follows them, whose last lines only say "RuntimeError".
    err = (p.stderr or "").strip()
    sout = (p.stdout or "").strip()
    # Keep the whole thing. pdb2pqr writes its own log beside the output, but
    # only while it is alive to write it: on the run that prompted this, that
    # log simply stops mid-sentence and the last word is in the captured
    # output instead. Record the exit status too: it is the difference
    # between a program that raised and one that was killed, and no amount of
    # reading the text will tell you which.
    try:
        with open(os.path.splitext(out_path)[0] + ".pdb2pqr.txt", "w") as fh:
            fh.write(" ".join(cmd) + "\n\nexit status: %s\n\n" % p.returncode
                     + "--- stdout ---\n" + sout + "\n\n--- stderr ---\n" + err)
    except OSError:
        pass

    def _tail(text, n=4):
        return [l.strip() for l in text.splitlines() if l.strip()][-n:]

    lines = [l.strip() for l in (sout + "\n" + err).splitlines() if l.strip()]
    # CRITICAL is pdb2pqr's fatal level, so it is used when present.
    # Otherwise the last thing said is used, rather than searching for a
    # level marker.
    #
    # ERROR lines are not used: pdb2pqr logs ERROR for a malformed header
    # it then skips, so an ERROR line can point at a harmless record
    # while the real cause appears below it, and a Python traceback's
    # last line ("KeyError: ...") carries no marker at all.
    crit = [l for l in lines if l.startswith("CRITICAL")]
    why = crit[:3] or _tail(err) or _tail(sout)
    detail = " / ".join(why)[:500] if why else "no output at all"
    # A negative return code is a signal: -9 is the out-of-memory killer, -11
    # a segfault. Neither leaves anything in a log, and both look like "it
    # just stopped" unless the number is reported.
    if p.returncode is not None and p.returncode < 0:
        detail = ("killed by signal %d (nothing it logged says so); "
                  % -p.returncode) + detail
    else:
        detail = "exit %s: " % p.returncode + detail
    return False, detail


CHAIN_POOL = ("ABCDEFGHIJKLMNOPQRSTUVWXYZ"
              "abcdefghijklmnopqrstuvwxyz0123456789")


# How far apart a residue's C and the next residue's N must be before they
# are treated as not bonded. A peptide bond is 1.33 A.
#
# Low-resolution structures stretch real bonds: in one 3.6 A crystal
# structure the C to N distances ran up to 1.7 A, with nothing beyond. A real
# break cannot be shorter than about 3 A, since one missing residue puts the
# next N some 4.5 A from the previous C. So 2.5 A separates the two with room
# on both sides.
GAP_TOLERANCE = 2.5


def backbone_gaps(path, tol=GAP_TOLERANCE):
    """Consecutive residues whose C-N distance is too long to be a bond.

    Returns [((chain, resnum, icode), (chain, resnum, icode), distance)].

    This is the measurement that predicts whether protonation will run at
    all. pdb2pqr treats a backbone gap as an incomplete structure and gives
    up on the WHOLE file (not on the fragment, on the file) so a receptor
    with one unmodelled loop anywhere in it came out with no hydrogens, no
    titration states and, because minimisation needs hydrogens, no
    minimisation either.
    """
    import math
    per, keys = {}, []
    for l in open(path, errors="replace"):
        if l.startswith("ENDMDL"):
            break
        if l[:6].strip() not in ("ATOM", "HETATM"):
            continue
        nm = l[12:16].strip()
        if nm not in ("C", "N"):
            continue
        try:
            xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            key = (l[21], int(l[22:26]), l[26])
        except ValueError:
            continue
        if key not in per:
            keys.append(key)
        per.setdefault(key, {})[nm] = xyz
    # IN THE ORDER THE FILE HOLDS THEM, which is the order the chain runs in.
    #
    # This sorted the residues by number, and a residue number is not a
    # position: thrombin is numbered the chymotrypsin way, so chain H runs
    # 16, 17, ... 221, 221A, 222 ... and its light chain restarts at 1.
    # Sorting put residues side by side that the chain does not, and each
    # such pair came back as a break several angstroms wide. Six of the
    # seven breaks reported on 4DII were that, and the caps written for them
    # landed on the amide nitrogen of the residue that really follows: the
    # complex then had six pairs of atoms 0.07 to 0.11 A apart, and
    # mmgbsa.write_complex refused every one of them.
    out = []
    for a, b in zip(keys, keys[1:]):
        if a[0] != b[0]:
            continue
        c, n = per[a].get("C"), per[b].get("N")
        if c is None or n is None:
            continue
        d = math.dist(c, n)
        if d > tol:
            out.append((a, b, d))
    return out


def split_at_gaps(path, out_path, gaps):
    """Each backbone fragment written as its own chain. Returns {new: old}.

    A gap in the backbone is a chain break. Told that the fragments are
    separate chains, pdb2pqr protonates each of them and caps its ends, which
    is what the deposited file already implies. Given the gaps as they are, it
    refuses the whole file. --nodebump also avoids the refusal, but only by
    switching off the step that resolves steric clashes.

    The letters are temporary: restore_chains puts the originals back before
    anything downstream reads the file, because everything downstream matches
    residues on (chain, number).
    """
    cut = {g[0] for g in gaps}
    used = set()
    for l in open(path, errors="replace"):
        if l[:6].strip() in ("ATOM", "HETATM") and len(l) > 21:
            used.add(l[21])
    free = [c for c in CHAIN_POOL if c not in used]
    # Per chain, and only the first model: a single previous-residue
    # marker for the whole file assigns the new chain letter to whatever
    # chain follows the cut residue in file order.
    mapping, out = {}, []
    cur, prev = {}, {}
    for l in open(path, errors="replace"):
        if l.startswith("ENDMDL"):
            break
        if l[:6].strip() not in ("ATOM", "HETATM"):
            if l.startswith(("END", "TER")):
                out.append(l)
            continue
        try:
            chain = l[21]
            key = (chain, int(l[22:26]), l[26])
        except ValueError:
            out.append(l)
            continue
        if prev.get(chain) is not None and prev[chain] in cut \
                and key != prev[chain]:
            if not free:
                return None            # nothing left to rename with
            cur[chain] = free.pop(0)
            out.append("TER\n")
        prev[chain] = key
        letter = cur.get(chain, chain)
        mapping[letter] = chain
        out.append(l[:21] + letter + l[22:])
    if not any(v != k for k, v in mapping.items()):
        return None
    with open(out_path, "w") as fh:
        fh.writelines(out)
        if not out or not out[-1].startswith("END"):
            fh.write("END\n")
    return mapping


def restore_chains(path, mapping, out_path=None):
    """Put the original chain letters back. Returns the path written.

    The PQR written beside the PDB receives the same treatment, so that the
    two files agree about which chain every residue is on.
    """
    out_path = out_path or path
    for src in (path, os.path.splitext(path)[0] + ".pqr"):
        if not os.path.isfile(src):
            continue
        dst = out_path if src == path else os.path.splitext(out_path)[0] + ".pqr"
        rows = []
        for l in open(src, errors="replace"):
            if l[:6].strip() in ("ATOM", "HETATM") and len(l) > 21:
                l = l[:21] + mapping.get(l[21], l[21]) + l[22:]
            rows.append(l)
        with open(dst, "w") as fh:
            fh.writelines(rows)
    return out_path


#: The caps this adds at a chain break. Not restrained during minimisation,
#: and not counted as part of the measured structure.
CAP_RESIDUES = ("ACE", "NME", "NMA")


def _unit(v):
    import numpy as np
    n = float(np.linalg.norm(v))
    return v / n if n > 1e-6 else v


def _dihedral(a, b, c, d):
    """The a-b-c-d torsion in degrees, in the sense _place builds one.

    IUPAC's: 0 with a and d on the same side of b-c (cis), 180 on opposite
    sides (trans), positive when a must turn clockwise, looking from b to c,
    to eclipse d. It measured 180 less than that, and cap_breaks places an
    NME's nitrogen from this torsion with _place, so the nitrogen went to the
    mirror of where it belongs: on 7A3R as close as 0.65 A to the carbonyl
    oxygen it should be trans to, with O-C-N as small as 29 degrees where
    about 123 is right.
    """
    import math
    import numpy as np
    b1, b2, b3 = b - a, c - b, d - c
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return math.degrees(math.atan2(
        float(np.linalg.norm(b2)) * float(np.dot(b1, n2)),
        float(np.dot(n1, n2))))


def _place(a, b, c, bond, angle, dihedral):
    """A fourth atom, at `bond` from c, `angle` to b-c and `dihedral` to a-b.

    The standard construction, and the reason the caps are built with it:
    they were placed by adding scaled unit vectors instead, and for the ACE
    every one of those vectors reduced to the same axis. The carbonyl oxygen
    and the methyl carbon came out collinear with the nitrogen and 0.279 A
    apart, in every ACE of every structure this has ever capped. Angles and
    torsions cannot collapse that way.
    """
    import math
    import numpy as np
    angle, dihedral = math.radians(angle), math.radians(dihedral)
    bc = _unit(c - b)
    n = _unit(np.cross(b - a, bc))
    m = np.cross(n, bc)
    d = np.array([-bond * math.cos(angle),
                  bond * math.sin(angle) * math.cos(dihedral),
                  bond * math.sin(angle) * math.sin(dihedral)])
    return c + d[0] * bc + d[1] * m + d[2] * n


def _methyl(root, one, two, bond=1.090, angle=109.5):
    """Three hydrogens on `two`, staggered about the one-two bond."""
    return [_place(root, one, two, bond, angle, t) for t in (60.0, 180.0,
                                                             300.0)]


def _room(x, atoms, more=()):
    """How far x is from the nearest of `atoms` (an n x 3 array) and `more`."""
    import numpy as np
    near = [float(np.linalg.norm(atoms - x, axis=1).min())] if len(atoms) else []
    near += [float(np.linalg.norm(p - x)) for p in more]
    return min(near) if near else float("inf")


def cap_breaks(path, out_path, gaps, steps, log=print, amide_h=False):
    """ACE and NME on the ends a chain break creates. Returns out_path or None.

    A break is not a terminus. Splitting the chain so that pdb2pqr will
    protonate it leaves two artificial ends, and charging them as real termini
    puts a charge pair into the structure that the protein does not have.

    The residue after each break is given the amide hydrogen the cap makes
    it need, and the terminal hydrogens the protonation gave it are dropped:
    a capped residue is an ordinary amide, with one hydrogen on its
    nitrogen. `amide_h` is kept for the caller that knows reduce leaves none
    there at all, and no longer changes what is written.
    """
    import numpy as np
    if not gaps:
        return None
    rows = [l for l in open(path, errors="replace")]
    # index the atoms of the residues on each side of every break
    want = set()
    for a, b, _d in gaps:
        want.add(a)
        want.add(b)
    at = {}
    for l in rows:
        if l[:6].strip() not in ("ATOM", "HETATM"):
            continue
        try:
            key = (l[21], int(l[22:26]), l[26])
        except ValueError:
            continue
        if key in want:
            at.setdefault(key, {})[l[12:16].strip()] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])

    used = {}
    for l in rows:
        if l[:6].strip() in ("ATOM", "HETATM"):
            try:
                used.setdefault(l[21], set()).add(int(l[22:26]))
            except ValueError:
                pass

    def free_number(chain, near):
        taken = used.get(chain, set())
        for n in range(near, near + 200):
            if n not in taken:
                taken.add(n)
                return n
        return near

    def rec(serial, name, resname, chain, resnum, xyz, element):
        nm = (" " + name).ljust(4) if len(name) < 4 else name
        return ("ATOM  %5d %s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                "          %2s\n"
                % (serial, nm, resname, chain, resnum,
                   xyz[0], xyz[1], xyz[2], element))

    # Every heavy atom that stays, and whose it is, for the room an ACE is
    # turned to find below. The terminal oxygen of each residue before a
    # break goes (see below), and stands where the NME's nitrogen does.
    ends = {a for a, _b, _d in gaps}
    owner, others = [], []
    for l in rows:
        if l[:6].strip() not in ("ATOM", "HETATM") or _is_hydrogen(
                l[12:16], l[76:78]):
            continue
        try:
            xyz = [float(l[30:38]), float(l[38:46]), float(l[46:54])]
            key = (l[21], int(l[22:26]), l[26])
        except ValueError:
            continue
        if key in ends and l[12:16].strip() in ("OXT", "OT1", "OT2"):
            continue
        owner.append((key, l[12:16].strip()))
        others.append(xyz)
    others = np.array(others, dtype=float).reshape(-1, 3)
    placed = []       # the heavy atoms of the caps, as they are built

    add = {}          # (chain, resnum, icode) -> [lines to put after it]
    pre = {}          # (chain, resnum, icode) -> [lines to put before it]
    add_h = {}        # (chain, resnum, icode) -> [its amide hydrogen]
    n_cap = 0
    # The NMEs first: each stands where the residue's own carbonyl puts it,
    # so it has no freedom, and every ACE is then turned clear of all of them.
    nme_at = {}
    for a, b, _d in gaps:
        A = at.get(a, {})
        if {"N", "CA", "C", "O"} <= set(A):
            # The N is trans to the carbonyl O about CA-C, in its plane.
            psi = _dihedral(A["N"], A["CA"], A["C"], A["O"])
            npos = _place(A["N"], A["CA"], A["C"], 1.329, 116.2, psi + 180.0)
            mpos = _place(A["O"], A["C"], npos, 1.449, 121.5, 0.0)
            nme_at[a] = (npos, mpos)
            placed += [npos, mpos]
    # WHOLE RESIDUES, NAMED AS THE FORCE FIELDS NAME THEM. The caps used to
    # be two or three heavy atoms with no hydrogens, and the methyl carbon of
    # the NME was called C, which is not an atom NME has. OpenMM could not
    # match either template and the minimisation stopped on the first one:
    #     No template found for residue 0 (ACE). The set of atoms is similar
    #     to NGLY, but is missing 1 H atom and 1 N atom.
    # so every receptor with an unmodelled break came out unminimised, and
    # the step that says so is one line in a log nobody reads twice.
    for a, b, _d in gaps:
        A, B = at.get(a, {}), at.get(b, {})
        # NME on the fragment that ENDS at `a`: the amide, then its methyl.
        if a in nme_at:
            npos, mpos = nme_at[a]
            hpos = _place(A["O"], A["C"], npos, 1.010, 119.0, 180.0)
            num = free_number(a[0], a[1] + 1)
            rows_out = [rec(0, "N", "NME", a[0], num, npos, "N"),
                        rec(0, "H", "NME", a[0], num, hpos, "H"),
                        rec(0, "CH3", "NME", a[0], num, mpos, "C")]
            for i, h in enumerate(_methyl(A["C"], npos, mpos), 1):
                rows_out.append(rec(0, "HH3%d" % i, "NME", a[0], num, h, "H"))
            add.setdefault(a, []).extend(rows_out)
            n_cap += 1
        # ACE on the fragment that BEGINS at `b`: carbonyl C, its O, methyl.
        if {"N", "CA", "C"} <= set(B):
            # Where the carbonyl stands about b's N-CA bond, b's phi, was
            # never measured: the residue the cap stands for is missing. It
            # is built extended, as it always was, unless that brings one of
            # its heavy atoms within 2.8 A of another (b's N and CA aside,
            # which it is bonded to or one bond beyond), and then at the
            # first turn that clears them all, or failing that the one that
            # leaves the most room. Extended put the ACE's oxygen 2.1 A from
            # the next residue's own side chain twice on 7A3R, and 1.2 A from
            # the methyl of the NME across a two-residue gap.
            keep = np.array([not (k == b and n in ("N", "CA"))
                             for k, n in owner], dtype=bool)
            near = others[keep]
            best = None
            for phi in (180.0, -150.0, -120.0, -90.0, -60.0, -30.0, 0.0,
                        30.0, 60.0, 90.0, 120.0, 150.0):
                cpos = _place(B["C"], B["CA"], B["N"], 1.335, 121.7, phi)
                opos = _place(B["CA"], B["N"], cpos, 1.229, 122.7, 0.0)
                mpos = _place(B["CA"], B["N"], cpos, 1.508, 116.6, 180.0)
                room = min(_room(x, near, placed) for x in (cpos, opos, mpos))
                if best is None or room > best[0]:
                    best = (room, cpos, opos, mpos)
                if room >= 2.8:
                    break
            cpos, opos, mpos = best[1:]
            placed += [cpos, opos, mpos]
            num = free_number(b[0], max(1, b[1] - 1))
            rows_out = [rec(0, "C", "ACE", b[0], num, cpos, "C"),
                        rec(0, "O", "ACE", b[0], num, opos, "O"),
                        rec(0, "CH3", "ACE", b[0], num, mpos, "C")]
            for i, h in enumerate(_methyl(B["N"], cpos, mpos), 1):
                rows_out.append(rec(0, "HH3%d" % i, "ACE", b[0], num, h, "H"))
            pre.setdefault(b, []).extend(rows_out)
            n_cap += 1

    if not n_cap:
        return None

    # The cap replaces the terminal group rather than sitting beside it.
    # pdb2pqr treats each fragment end as a real end and completes it: an
    # OXT on the last residue, and H2 and H3 on the first. A cap added
    # without removing those gives a residue with both a C-terminal
    # carboxylate and a peptide bond onward, which no force field has a
    # template for.
    drop = set()
    for a, b, _d in gaps:
        for nm in ("OXT", "HXT", "OT1", "OT2"):
            drop.add((a, nm))
        # ALL THREE, ALWAYS. A capped residue is an ordinary amide: one
        # hydrogen on the nitrogen, named H. H1 was kept here unless the
        # reduce path asked otherwise, and pdb2pqr puts H1 where the
        # PRECEDING CARBONYL CARBON goes - so the cap's C landed 0.42 A from
        # it on thrombin, which write_complex then refused as two atoms at
        # the same place. The amide hydrogen below replaces all three.
        for nm in ("H1", "H2", "H3"):
            drop.add((b, nm))
        # And the one this pdb2pqr calls H, which is the first of the three
        # (H, H2, H3), not an amide hydrogen: kept, it stood where the cap's
        # carbonyl carbon goes, 0.41 A from it on 7A3R, or out of the amide
        # plane on every other break.
        if b in pre:
            drop.add((b, "H"))

    # THE AMIDE HYDROGEN, placed from the cap's carbonyl carbon, which is
    # the nitrogen's new neighbour, in the plane of the amide and trans to
    # the carbonyl oxygen. Written for every capped residue, whichever
    # protonation ran: the terminal hydrogens it put there have just been
    # dropped, so without this the nitrogen has none. A proline's nitrogen
    # carries none.
    if True:
        for a, b, _d in gaps:
            B = at.get(b, {})
            if not {"N", "CA"} <= set(B):
                continue
            ace = [c for c in pre.get(b, []) if c[12:16].strip() == "C"]
            if not ace:
                continue
            resname = next((l[17:20] for l in rows
                            if l[:6].strip() in ("ATOM", "HETATM")
                            and l[21] == b[0] and l[22:26].strip() == str(b[1])
                            and l[26] == b[2]), "UNK")
            if resname.strip() == "PRO":
                continue
            cpos = np.array([float(ace[0][30:38]), float(ace[0][38:46]),
                             float(ace[0][46:54])])
            hpos = B["N"] - 1.01 * _unit(_unit(cpos - B["N"])
                                         + _unit(B["CA"] - B["N"]))
            add_h.setdefault(b, []).append(
                rec(0, "H", resname.strip(), b[0], b[1], hpos, "H"))

    out, serial, seen = [], 0, None
    for l in rows:
        if l[:6].strip() in ("ATOM", "HETATM"):
            try:
                _k = (l[21], int(l[22:26]), l[26])
            except ValueError:
                _k = None
            if _k is not None and (_k, l[12:16].strip()) in drop:
                continue
        if l[:6].strip() not in ("ATOM", "HETATM"):
            if l.startswith(("TER", "END")):
                if seen in add:
                    for c in add.pop(seen):
                        serial += 1
                        out.append(c[:6] + "%5d" % serial + c[11:])
                out.append(l)
            continue
        try:
            key = (l[21], int(l[22:26]), l[26])
        except ValueError:
            out.append(l)
            continue
        if seen is not None and key != seen and seen in add:
            for c in add.pop(seen):
                serial += 1
                out.append(c[:6] + "%5d" % serial + c[11:])
        if key != seen and key in pre:
            for c in pre.pop(key):
                serial += 1
                out.append(c[:6] + "%5d" % serial + c[11:])
        seen = key
        serial += 1
        out.append(l[:6] + "%5d" % serial + l[11:])
        if key in add_h and l[12:16].strip() == "N":
            for c in add_h.pop(key):
                serial += 1
                out.append(c[:6] + "%5d" % serial + c[11:])
    for leftover in (add, pre, add_h):
        for _k, lines in leftover.items():
            for c in lines:
                serial += 1
                out.append(c[:6] + "%5d" % serial + c[11:])
    if not out or not out[-1].startswith("END"):
        out.append("END\n")
    with open(out_path, "w") as fh:
        fh.writelines(out)
    steps.append(Step(
        "capping breaks", True,
        "%d break end(s) capped with ACE/NME. A break is not a terminus, and "
        "leaving it to be protonated as one puts an NH3+ and a COO- at every "
        "hole in the chain, charges the protein does not have. The real "
        "termini are left charged, because those are" % n_cap))
    return out_path


def protonate(path, out_path, steps, ph=7.4, forcefield="AMBER", log=print):
    """Hydrogens and titration states at a stated pH, via pdb2pqr and PROPKA.

    PROPKA is tried first and plain force-field protonation second, because the
    two fail for different reasons.
    """
    exe = find_tool("pdb2pqr", "pdb2pqr30")
    if not exe:
        steps.append(Step("protonation", False,
                          "pdb2pqr is not installed, so no hydrogens were "
                          "added and no titration states assigned"))
        return None

    # A backbone gap is a chain break, and pdb2pqr has to be told so. It
    # refuses a structure whose backbone jumps ("Found gap in biomolecule
    # structure"), and it refuses the whole file, which would leave the
    # receptor with no hydrogens, no titration states and, since
    # minimisation needs hydrogens, no minimisation. Naming the fragments
    # as chains costs nothing. --nodebump also avoids the refusal by
    # switching off clash resolution; it remains the fallback when no
    # chain letters are left to rename with.
    src, mapping = path, None
    gaps = backbone_gaps(path)
    if gaps:
        split = os.path.join(os.path.dirname(out_path), "split_at_gaps.pdb")
        mapping = split_at_gaps(path, split, gaps)
        where = ", ".join("%d-%d (%.1f A)" % (a[1], b[1], d)
                          for a, b, d in gaps[:4])
        if mapping:
            src = split
            steps.append(Step(
                "backbone gaps", True,
                "%d unmodelled break(s) (%s%s) written as separate chains "
                "so pdb2pqr protonates the whole structure instead of "
                "refusing it. The fragment ends become termini, which is what "
                "the file already implies; the original chain letters go back "
                "on afterwards" % (len(gaps), where,
                                   ", ..." if len(gaps) > 4 else "")))
        else:
            steps.append(Step(
                "backbone gaps", False,
                "%d unmodelled break(s) (%s) and no free chain letter to "
                "name the fragments with; pdb2pqr will be asked not to debump "
                "instead" % (len(gaps), where)))

    ok, why = _pdb2pqr(exe, src, out_path, ph, forcefield, propka=True)
    if not ok and gaps and not mapping:
        ok, why = _pdb2pqr(exe, src, out_path, ph, forcefield, propka=True,
                           nodebump=True)
    if ok and mapping:
        restore_chains(out_path, mapping)
    fallback = False
    if not ok:
        steps.append(Step("protonation", False,
                          f"pdb2pqr with PROPKA failed: {why}"))
        if "integral" in why or "charge" in why:
            steps.append(Step(
                "protonation note", True,
                "a non-integer charge usually means a group the force field "
                "has no parameters for, such as ACE or NME caps from an earlier "
                "preparation, or a modified residue. Prepare the original "
                "structure rather than a prepared one, or keep_hetatms=False "
                "to drop what is not protein."))
        log("  protonation with PROPKA failed; retrying without it")
        ok, why2 = _pdb2pqr(exe, src, out_path, ph, forcefield,
                            propka=False)
        if not ok and gaps:
            ok, why2 = _pdb2pqr(exe, src, out_path, ph, forcefield,
                                propka=False, nodebump=True)
        if ok and mapping:
            restore_chains(out_path, mapping)
        fallback = ok
        if not ok:
            steps.append(Step("protonation retry", False,
                              f"pdb2pqr without PROPKA also failed: {why2}"))
            steps.append(Step(
                "protonation consequence", True,
                "WITHOUT HYDROGENS the minimisation cannot run either. The "
                "docking will add them geometrically and the run's pH will "
                "not be reflected in the receptor."))
            return None

    n_h = sum(1 for l in open(out_path, errors="replace")
              if l.startswith(("ATOM", "HETATM"))
              and _is_hydrogen(l[12:16], l[76:78]))
    if n_h == 0:
        # pdb2pqr exited cleanly and wrote a file with nothing on it. Treated
        # as a failure, because the only thing that matters here is whether
        # there are hydrogens, and there are not.
        steps.append(Step("protonation", False,
                          "pdb2pqr wrote a structure with no hydrogens on it"))
        return None
    if fallback:
        steps.append(Step(
            "protonation retry", True,
            f"{forcefield} force field WITHOUT PROPKA; {n_h} hydrogens added. "
            f"Titration states are the force field's standard ones, so pH "
            f"{ph} is not reflected in this structure; histidines in "
            f"particular may be in the wrong state."))
    else:
        steps.append(Step("protonation", True,
                          f"{forcefield} force field at pH {ph} with PROPKA "
                          f"titration states; {n_h} hydrogens added"))
    return out_path


#: How the receptor is protonated, by the names job.json and the page use:
#: reduce, as ADFR's own preparation and ADapt do it, or PDB2PQR with PROPKA
#: at a chosen pH.
PROTONATION = ("reduce", "propka")
#: What a job file written before the choice existed was protonated with.
#: It has no protonation key, and must still be prepared the same way, so
#: that an old run reproduces.
PROTONATION_LEGACY = "propka"


def protonation_of(job_or_value):
    """'reduce' or 'propka' for a job dict or a stored value. Absent or
    unknown reads as the legacy method, never as the new default."""
    v = job_or_value.get("protonation") if isinstance(job_or_value, dict) \
        else job_or_value
    v = str(v or "").strip().lower()
    return v if v in PROTONATION else PROTONATION_LEGACY


def reduce_version(exe):
    """reduce's own version line, such as 'reduce.3.23.130521', or ''."""
    if not exe:
        return ""
    try:
        p = subprocess.run([exe, "-version"], capture_output=True, text=True,
                           timeout=30)
    except (OSError, subprocess.SubprocessError):
        return ""
    for line in (p.stdout + "\n" + p.stderr).splitlines():
        if line.strip().startswith("reduce."):
            return line.strip()
    return ""


def reduce_het_dict(exe):
    """The het dictionary that ships beside reduce in the ADCP 1.1 suite,
    which is what ADFR's own protonation passes with -DB, or None."""
    if not exe:
        return None
    for cand in (exe + "_wwPDB_het_dict.txt",
                 os.path.join(os.path.dirname(exe),
                              "reduce_wwPDB_het_dict.txt")):
        if os.path.isfile(cand):
            return cand
    return None


def his_states(path):
    """{(chain, number with insertion code): 'HID' | 'HIE' | 'HIP' | 'HIS'}
    for every histidine, read from the hydrogens on its ring nitrogens."""
    hs = {}
    for l in open(path, errors="replace"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        if l[17:20] not in ("HIS", "HID", "HIE", "HIP"):
            continue
        hs.setdefault((l[21], l[22:27].strip()), set()).add(
            l[12:16].strip().upper())
    out = {}
    for k, names in hs.items():
        d, e = "HD1" in names, "HE2" in names
        out[k] = "HIP" if d and e else "HID" if d else "HIE" if e else "HIS"
    return out


def reduce_decisions(text):
    """reduce's USER MOD lines for His, Asn and Gln, as
    [(chain, number, residue, decision, turned over)], for example
    ('A', '149', 'HIS', 'no HD1', True)."""
    out = []
    for l in text.splitlines():
        if not l.startswith("USER  MOD"):
            continue
        parts = l.split(":")
        if len(parts) < 4 or len(parts[1]) < 3:
            continue
        field = parts[1][2:].split()
        if len(field) < 2 or field[1] not in ("HIS", "ASN", "GLN"):
            continue
        what = parts[2].strip()
        turned = what.startswith("FLIP")
        out.append((parts[1][1], field[0], field[1],
                    what[4:].strip() if turned else what, turned))
    return out


def _first_numbers(path):
    """The number of the first residue of each chain, by chain."""
    firsts = {}
    for l in open(path, errors="replace"):
        if l.startswith(("ATOM", "HETATM")) and l[21] not in firsts:
            try:
                firsts[l[21]] = int(l[22:26])
            except ValueError:
                continue
    return firsts


def _count_hydrogens(path):
    """Hydrogens in a PDB file, decided the way _is_hydrogen decides."""
    n = 0
    try:
        for l in open(path, errors="replace"):
            if l.startswith(("ATOM", "HETATM")) and _is_hydrogen(l[12:16],
                                                                 l[76:78]):
                n += 1
    except OSError:
        return 0
    return n


def _drop_hydrogens(path, out_path):
    """Every hydrogen line taken out, heavy atoms untouched. Count removed."""
    n = 0
    with open(path, errors="replace") as fh, open(out_path, "w") as out:
        for l in fh:
            if l.startswith(("ATOM", "HETATM")) and _is_hydrogen(l[12:16],
                                                                 l[76:78]):
                n += 1
                continue
            out.write(l)
    return n


def _trim_hydrogens(path, out_path, steps, exe, log=print):
    """`path` with the hydrogens it arrived with removed, or `path` itself.

    reduce -Trim is reduce's own way of doing this and the first half of the
    two-pass order it documents for an already-protonated file. It is run only
    when there is something to trim, so a deposited structure with no
    hydrogens takes the same path it always did.

    The result is checked rather than trusted. If -Trim could not run, or left
    hydrogens behind, the remainder is dropped here by name instead: the point
    of the trim is that no hydrogen reduce did not place can reach the force
    field, and a trim that half worked would leave exactly the stray atom this
    exists to remove.
    """
    had = _count_hydrogens(path)
    if not had:
        return path
    work = os.path.dirname(os.path.abspath(out_path))
    os.makedirs(work, exist_ok=True)      # run() has made it; a direct caller may not
    trimmed = os.path.join(work, "trimmed.pdb")
    why = ""
    try:
        p = subprocess.run([exe, "-Trim", path], capture_output=True,
                           text=True, timeout=3600)
        text = p.stdout or ""
        if text.count("\n") > 1 and _atoms(text.splitlines()):
            with open(trimmed, "w") as fh:
                fh.write(text)
        else:
            why = "reduce -Trim wrote no structure"
    except (OSError, subprocess.SubprocessError) as e:
        why = "reduce -Trim could not be run: %s" % e
    left = _count_hydrogens(trimmed) if os.path.isfile(trimmed) else None
    if why or left is None or left:
        # By name, so the guarantee holds whatever reduce did.
        n = _drop_hydrogens(path, trimmed)
        steps.append(Step(
            "trim hydrogens", True,
            "%d hydrogen(s) the structure arrived with were removed by name "
            "so reduce places every one of them%s"
            % (n, (" (%s)" % why) if why
               else (" (reduce -Trim left %d behind)" % left) if left
               else "")))
        log("  the structure's own hydrogens were removed; reduce places "
            "every hydrogen from here")
        return trimmed
    steps.append(Step(
        "trim hydrogens", True,
        "reduce -Trim removed the %d hydrogen(s) the structure arrived with, "
        "so reduce places every one of them" % had))
    log("  reduce -Trim removed the structure's own %d hydrogen(s); reduce "
        "places every hydrogen from here" % had)
    return trimmed


def protonate_reduce(path, out_path, steps, exe=None, log=print):
    """Hydrogens placed by reduce, the way ADFR's own preparation and ADapt
    place them.

    reduce -BUILD adds every hydrogen, turns Asn and Gln amides and His rings
    over where the hydrogen bonds say the file has them the wrong way round,
    and chooses which histidine nitrogens carry a proton from the same
    network. It takes no pH: Asp and Glu stay charged, Lys and Arg
    positive, and a histidine is given two protons only when its
    surroundings call for it. On the protein structures in the PDB solved
    by neutron diffraction it placed histidine protons as deposited in 55%
    of 1,457 cases, where pdb2pqr with PROPKA placed 31%.

    The flags are ADFR's (-BUILD with the het dictionary beside reduce),
    and one more: -Nterm set to the first residue number of the chains, so
    a chain that starts at residue 25 is charged as a terminus there, as
    pdb2pqr treats it. Without it reduce puts no hydrogen on that nitrogen
    at all, and the minimisation's force field has no template for the
    residue.

    HYDROGENS THE STRUCTURE ARRIVES WITH ARE TRIMMED FIRST, with reduce's own
    -Trim, which is the two-pass order reduce documents for a file that is
    already protonated. -BUILD adds what is missing and leaves what it finds,
    and what it finds was placed by whoever wrote the file. One such atom is
    enough to lose the minimisation: 7A3R arrives with 2,859 hydrogens, HXT
    among them, which is a proton on the C-terminal carboxylate. That makes
    the terminus neutral, Amber's templates have only the charged form, and
    createSystem then refuses the whole structure for want of a template for
    one residue - measured, and measured again with that single atom removed,
    where it is accepted. pdb2pqr never meets this because it rebuilds every
    hydrogen and decides the terminus itself; after the trim, so does reduce.

    A structure that arrives with no hydrogens is not trimmed and is
    protonated exactly as it was before any of this.

    Returns out_path, or None when reduce could not be run or wrote no
    hydrogens, with the reason in the step.
    """
    if not exe or not os.path.isfile(exe):
        steps.append(Step("reduce", False,
                          "reduce was not found in the ADCP 1.1 suite"))
        return None
    het = reduce_het_dict(exe)
    path = _trim_hydrogens(path, out_path, steps, exe, log)
    nterm = max([1] + list(_first_numbers(path).values()))
    cmd = [exe, "-BUILD", "-Nterm%d" % nterm]
    if het:
        cmd += ["-DB", het]
    cmd.append(path)
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except (OSError, subprocess.SubprocessError) as e:
        steps.append(Step("reduce", False, f"reduce could not be run: {e}"))
        return None
    # reduce's exit status is not a verdict: it returns 1 on structures it
    # protonated completely. What it wrote is.
    text = p.stdout or ""
    try:
        with open(os.path.splitext(out_path)[0] + ".reduce.txt", "w") as fh:
            fh.write(" ".join(cmd) + "\n\nexit status: %s\n\n--- stderr ---\n%s"
                     % (p.returncode, p.stderr or ""))
    except OSError:
        pass
    if "USER  MOD reduce" not in text[:4000]:
        said = [l.strip() for l in (p.stderr or "").splitlines() if l.strip()]
        steps.append(Step("reduce", False, "reduce wrote no structure (%s)"
                          % (" / ".join(said[-2:]) or "no output")[:300]))
        return None
    with open(out_path, "w") as fh:
        fh.write(text)
    n_h = sum(1 for l in text.splitlines()
              if l.startswith(("ATOM", "HETATM"))
              and _is_hydrogen(l[12:16], l[76:78]))
    if n_h == 0:
        steps.append(Step("reduce", False,
                          "reduce wrote a structure with no hydrogens on it"))
        return None
    his = his_states(out_path)
    counts = {}
    for v in his.values():
        counts[v] = counts.get(v, 0) + 1
    flips = {}
    for _c, _n, res, _d, flipped in reduce_decisions(text):
        if flipped:
            flips[res] = flips.get(res, 0) + 1
    ver = reduce_version(exe) or "reduce"
    detail = (f"{ver} -BUILD, as ADFR and ADapt run it, with -Nterm{nterm} "
              f"so each chain's first residue is a terminus; {n_h} hydrogens "
              f"added")
    if his:
        detail += "; histidines " + ", ".join(
            f"{k} x{counts[k]}" for k in ("HID", "HIE", "HIP", "HIS")
            if counts.get(k))
    detail += ("; turned over " + ", ".join(f"{v} {k.capitalize()}"
                                            for k, v in sorted(flips.items()))
               if flips else "; nothing turned over")
    steps.append(Step("protonation", True, detail))
    return out_path


def keep_chain(lines, chain, steps):
    """Only the chain(s) that will be docked against.

    `chain` may name several, comma-separated, which is what a run whose
    docking box spans an interface needs: the groove has two sides and
    preparing one of them leaves the peptide pressed against a hole.

    The run uses one chain. Preparing all of them rebuilds side chains, adds
    hydrogens and minimises nineteen copies of a protein that nothing will
    look at: on a twenty-chain cryo-EM structure that is the difference
    between a job that finishes and one that does not.
    """
    if not chain:
        return lines
    want = [c.strip() for c in str(chain).split(",") if c.strip()]
    seen = {l[21] for l in lines if l.startswith(("ATOM", "HETATM"))
            and len(l) > 21}
    absent = [c for c in want if c not in seen]
    if absent:
        steps.append(Step("chain", False,
                          f"chain {', '.join(absent)} is not in this file; found "
                          f"{', '.join(sorted(seen)) or 'none'}. Left as it is."))
        return lines
    keepset = set(want)
    def _keep(l):
        if l.startswith(("ATOM", "HETATM")):
            return len(l) > 21 and l[21] in keepset
        # Chain-scoped annotation for a chain being removed: SEQRES
        # names its chain in column 12, TER in column 22. SEQRES is
        # kept, since PDBFixer reads it to tell a real gap from an
        # intentional truncation.
        if l.startswith("SEQRES") and len(l) > 11:
            return l[11] in keepset
        if l.startswith("TER") and len(l) > 21 and l[21].strip():
            return l[21] in keepset
        return True

    kept = [l for l in lines if _keep(l)]
    others = sorted(seen - keepset)
    steps.append(Step("chain", bool(others),
                      f"kept chain{'s' if len(want) > 1 else ''} "
                      f"{', '.join(want)}"
                      + (f", dropped {', '.join(others)}" if others
                         else " (the only one)")))
    return kept


def _prep_key(pdb_in, ph, forcefield, chain, keep_waters_within,
              keep_hetatms, keep_ions, fill_loops, cap_termini,
              minimise_structure, restraint, rebuild_sidechains=True,
              protonation=PROTONATION_LEGACY):
    """Everything about a preparation that changes its result.

    The input file is identified by size and modification time as well as path,
    so editing the receptor between attempts invalidates the reuse. Hashing the
    contents would be better and is not worth reading a 14,000-atom file for on
    every start.
    """
    try:
        st = os.stat(pdb_in)
        stamp = [int(st.st_size), int(st.st_mtime)]
    except OSError:
        stamp = None
    key = {"receptor": os.path.abspath(pdb_in), "stamp": stamp,
           "ph": round(float(ph), 3), "forcefield": str(forcefield),
           "chain": str(chain), "waters": float(keep_waters_within),
           "hetatms": bool(keep_hetatms), "ions": bool(keep_ions),
           "loops": bool(fill_loops), "termini": bool(cap_termini),
           "minimise": bool(minimise_structure),
           "restraint": float(restraint)}
    # Recorded only when off. A preparation made before the setting existed
    # rebuilt side chains and has no such entry, and it must still be
    # recognised as made with the default.
    if not rebuild_sidechains:
        key["sidechains"] = False
    # The same for the protonation: a preparation made before it could be
    # chosen used PDB2PQR with PROPKA and has no entry. reduce takes no pH,
    # so with it the pH is not part of what the result depends on.
    if protonation_of(protonation) != PROTONATION_LEGACY:
        key["protonation"] = protonation_of(protonation)
        key["ph"] = None
    return key


def run(pdb_in, workdir, ph=7.4, forcefield="AMBER", log=print, chain="",
        keep_waters_within=0.0, site_points=None, keep_hetatms=False,
        keep_ions=True, fill_loops=False, cap_termini=True,
        minimise_structure=True, restraint=10.0, rebuild_sidechains=True,
        protonation=PROTONATION_LEGACY, reduce_exe=None):
    """Prepare a receptor. Returns the prepared PDB, or None if nothing ran.

    The order matters and is the wizard's: tidy the file first, rebuild what
    is missing, then protonate, then minimise. Protonating before rebuilding
    would put hydrogens on atoms that are about to move; minimising before
    protonating would minimise a structure with no hydrogens on it.

    `protonation` is "reduce" (reduce_exe, the ADCP 1.1 suite's) or "propka"
    (PDB2PQR with PROPKA at `ph`, what every run made before the choice
    existed used). Where reduce cannot run, PDB2PQR with PROPKA protonates
    instead and the log says so.

    Missing side chains with rebuild_sidechains off: PDB2PQR completes every
    residue it is given and the minimisation needs whole residues, so they
    are modelled for those two and taken out at the end. reduce protonates a
    residue as it stands, so with reduce they are modelled only when the
    minimisation will run, and not at all when it will not.
    """
    protonation = protonation_of(protonation)
    os.makedirs(workdir, exist_ok=True)

    # Already prepared by an earlier attempt at this run? A job killed at
    # its wall clock limit is resumed, and the resumed job would
    # otherwise repeat the whole preparation.
    want = _prep_key(pdb_in, ph, forcefield, chain, keep_waters_within,
                     keep_hetatms, keep_ions, fill_loops, cap_termini,
                     minimise_structure, restraint, rebuild_sidechains,
                     protonation)
    out = os.path.join(workdir, "prepared.pdb")
    keyf = os.path.join(workdir, "prepared.json")
    if os.path.isfile(out) and os.path.getsize(out) > 0:
        try:
            had = json.load(open(keyf))
        except (OSError, ValueError):
            had = None
        if had == want:
            log(f"  receptor already prepared for these settings, reusing "
                f"{os.path.basename(out)}")
            return out
        if had is not None:
            changed = [k for k in want if had.get(k) != want.get(k)]
            log("  re-preparing the receptor: "
                + ", ".join(changed) + " changed since it was prepared")

    steps = []
    log("  preparing the receptor")

    try:
        lines = as_pdb_lines(pdb_in, workdir, log)
    except OSError as e:
        log(f"  cannot read {pdb_in}: {e}")
        return None
    if not lines:
        return None

    n_before = len(_atoms(lines))
    lines = drop_metadata(lines, steps)
    lines = keep_chain(lines, chain, steps)

    # An alpha-carbon trace stops here rather than being rebuilt.
    # PDBFixer would try, since every side chain of every residue is
    # missing, and on a 5,768-residue trace it was still running after
    # twenty minutes. Refusing at once is better than a job that reaches
    # its time limit having done nothing else.
    ca_only, why = looks_like_trace(lines)
    if ca_only:
        for s in steps:
            log(str(s))
        log("  " + why)
        log("  Preparation stopped. Rebuilding every side chain is modelling "
            "rather than preparation, and it takes long enough to lose a run.")
        log("  What to do instead, in the order worth trying:")
        log("    1. Look for another structure of the same protein. A trace "
            "is usually one deposition among several, and the others often "
            "have side chains.")
        log("    2. Take the AlphaFold model: alphafold.ebi.ac.uk has a "
            "full-atom prediction for most UniProt entries, ready to use.")
        log("    3. Predict one yourself with ColabFold, which is already "
            "installed here for co-folding.")
        log("  A predicted model is a prediction: its backbone will not match "
            "the measured one exactly. Where the measurement matters, prefer "
            "another experimental structure over a model.")
        return None

    lines = pick_altloc(lines, steps)
    lines = standardise(lines, steps)
    lines = drop_hydrogens(lines, steps)
    lines = drop_waters(lines, steps, keep_waters_within, site_points)
    if not keep_hetatms:
        lines = drop_hetatms(lines, steps, keep_ions)
    disulfides(lines, steps)

    tidy = _write(os.path.join(workdir, "tidied.pdb"), lines)
    current = tidy

    # Missing side-chain atoms that will not be kept are modelled only for
    # what cannot work without them (see above).
    modelled = (not rebuild_sidechains
                and (protonation != "reduce" or minimise_structure))
    # Loops from PDBFixer's first seed; only when no relaxation can make them
    # fit are they built again from the next, and the log says so.
    source = current
    for k, seed in enumerate(LOOP_SEEDS):
        built = []
        fixed = fix_structure(source, os.path.join(workdir, "rebuilt.pdb"), steps,
                              fill_loops=fill_loops, cap_termini=cap_termini,
                              log=log, rebuild_sidechains=rebuild_sidechains,
                              model_missing=rebuild_sidechains or modelled, built=built,
                              seed=seed)
        current = fixed or source
        if not built:
            break
        try:
            current = relax_built(current, os.path.join(workdir, "relaxed.pdb"), steps,
                                  built, log=log)
            break
        except PrepError as e:
            if k + 1 < len(LOOP_SEEDS):
                steps.append(Step("rebuilt loops relaxed", False,
                                  f"the loops PDBFixer built from seed {seed} could not be made "
                                  f"to fit; building them again from seed {LOOP_SEEDS[k + 1]}"))
                continue
            steps.append(Step("rebuilt loops relaxed", False, str(e)))
            _write_prep_log(workdir, pdb_in, out, protonation, ph, forcefield, steps,
                            n_before, None)
            raise

    by_reduce = False
    if protonation == "reduce":
        prot = protonate_reduce(current,
                                os.path.join(workdir, "protonated.pdb"),
                                steps, exe=reduce_exe, log=log)
        by_reduce = prot is not None
        if not by_reduce:
            log("  reduce could not protonate the receptor; PDB2PQR with "
                "PROPKA does instead")
            steps.append(Step("protonation fallback", True,
                              f"PDB2PQR with PROPKA at pH {ph}, because "
                              f"reduce could not run"))
            prot = protonate(current, os.path.join(workdir, "protonated.pdb"),
                             steps, ph=ph, forcefield=forcefield, log=log)
            # pdb2pqr completes every residue it is given, so what it built
            # comes out at the end, as it would have had it been chosen.
            modelled = modelled or not rebuild_sidechains
    else:
        prot = protonate(current, os.path.join(workdir, "protonated.pdb"),
                         steps, ph=ph, forcefield=forcefield, log=log)
    current = prot or current

    # The caps go on here: after pdb2pqr, which has no parameters for ACE
    # or NME, and before OpenMM and tleap, which both have them.
    if cap_termini:
        _g = backbone_gaps(current)
        if _g:
            _capped = cap_breaks(current,
                                 os.path.join(workdir, "capped.pdb"), _g,
                                 steps, log=log, amide_h=by_reduce)
            current = _capped or current

    if minimise_structure:
        mini = minimise(current, os.path.join(workdir, "minimised.pdb"), steps,
                        restraint=restraint, log=log, free=built)
        current = mini or current
    else:
        steps.append(Step("minimisation", False, "not asked for"))

    # The rebuilt loops as they leave: the relaxation's checks again, on the
    # structure that is written, so nothing after it (protonation, the
    # minimisation) has twisted a bond or turned a residue unseen.
    if built:
        bad = junction_geometry(current, built)
        if bad:
            why = (f"{len(bad)} fault(s) in the peptide bonds and residues of the rebuilt loops "
                   f"after protonation and minimisation ({'; '.join(bad[:4])}). Untick "
                   f"'{LOOPS_SETTING}' to dock the structure as measured.")
            steps.append(Step("rebuilt loops checked", False, why))
            _write_prep_log(workdir, pdb_in, out, protonation, ph, forcefield, steps,
                            n_before, None)
            raise PrepError(why)
        steps.append(Step("rebuilt loops checked", True,
                          f"as written: every peptide bond in or beside the {len(built)} rebuilt "
                          f"residue(s) trans, every one of them of the L hand, none nearer than "
                          f"{JUNCTION['contact']:g} A to another atom"))

    # Last, so protonation and minimisation ran on whole residues, which
    # neither can do without: what was modelled for them comes out here.
    if modelled:
        bare = strip_modelled_side_chains(
            tidy, current, os.path.join(workdir, "sidechains_as_measured.pdb"),
            steps)
        current = bare or current

    # `out` was set above, where the reuse check needed it.
    if current != out:
        shutil.copy(current, out)

    for s in steps:
        log(str(s))
    n_after = len(_atoms(open(out, errors="replace").read().splitlines()))
    log(f"  {n_before} atoms in, {n_after} out -> {out}")

    # prepared.pdb is returned whether or not protonation ran, and the
    # path alone does not tell the caller which. A receptor without
    # hydrogens then fails later in agfr, with a message that points at
    # the grid step. The hydrogens are counted rather than inferred from
    # the protonation step's status, so a protonation that succeeded but
    # added nothing is caught too.
    n_h = sum(1 for line in open(out, errors="replace")
              if line.startswith(("ATOM", "HETATM"))
              and _is_hydrogen(line[12:16], line[76:78]))
    _write_prep_log(workdir, pdb_in, out, protonation, ph, forcefield, steps,
                    n_before, n_after)
    if n_h == 0:
        # Stopped, not warned: a banner in a log nobody opens let a receptor
        # with no hydrogens be docked, at a pH it never had.
        raise PrepError("the prepared receptor has no hydrogens: protonation failed (the "
                        "protonation line in preparation.log says why). Docking it would "
                        "ignore the run's pH, so the run stops here.")
    # Record what this preparation was made from, so the next attempt at this
    # run can tell "already done" from "done differently".
    try:
        with open(keyf, "w") as fh:
            json.dump(want, fh, indent=1)
    except OSError:
        pass                      # reuse is an optimisation, not the result
    return out


def _write_prep_log(workdir, pdb_in, out, protonation, ph, forcefield, steps,
                    n_before, n_after):
    """preparation.log, which the results page reads. Written on the way out of
    a failed preparation too, so the step that stopped it is on the page."""
    with open(os.path.join(workdir, "preparation.log"), "w") as fh:
        fh.write("Protein preparation\n===================\n\n")
        fh.write(f"in  : {pdb_in}\nout : {out}\n"
                 + (f"protonation: reduce   force field: {forcefield}\n\n"
                    if protonation == "reduce" else
                    f"pH  : {ph}   force field: {forcefield}\n\n"))
        for s in steps:
            fh.write(str(s) + "\n")
        fh.write(f"\n{n_before} atoms in, "
                 + (f"{n_after} out\n" if n_after is not None else "stopped before an output\n"))


def summarise(prepared_pdb, log=print):
    """What came out, in one line, for a caller that only wants the shape."""
    n_atom = n_h = 0
    res = set()
    try:
        for l in open(prepared_pdb, errors="replace"):
            if not l.startswith(("ATOM", "HETATM")):
                continue
            n_atom += 1
            if l[76:78].strip() == "H" or l[12:16].strip().startswith("H"):
                n_h += 1
            res.add((l[21], l[22:27]))
    except OSError:
        return
    log(f"  prepared: {len(res)} residues, {n_atom} atoms ({n_h} hydrogens)")
