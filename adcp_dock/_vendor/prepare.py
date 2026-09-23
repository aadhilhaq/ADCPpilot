"""Protein preparation, in the sense of Schrodinger's Protein Preparation Wizard.

pdb2pqr alone adds hydrogens and assigns titration states. A crystal
structure needs more before it can be docked against:

* one alternate location, not two overlapping copies of a side chain
* selenomethionine read as methionine
* waters removed, except any requested near the site
* other heteroatoms removed unless wanted, since an arbitrary cofactor
  cannot be parameterised here
* missing side-chain atoms rebuilt, and missing residues if requested
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


def fix_structure(path, out_path, steps, fill_loops=False, cap_termini=True,
                  log=print):
    """Missing side-chain atoms rebuilt, and optionally missing residues.

    PDBFixer's findMissingResidues finds gaps between residues that are
    present. Filling them is off by default: a rebuilt loop is a prediction,
    and once written it is indistinguishable from measured coordinates.
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
        if not fill_loops:
            fixer.missingResidues = {}
        fixer.findMissingAtoms()
        n_side = sum(len(v) for v in fixer.missingAtoms.values())
        n_term = sum(len(v) for v in fixer.missingTerminals.values())
        if not cap_termini:
            fixer.missingTerminals = {}
            n_term = 0
        fixer.addMissingAtoms()
        with open(out_path, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh,
                              keepIds=True)
    except Exception as e:
        steps.append(Step("missing atoms", False, f"pdbfixer failed: {e}"))
        return None

    detail = f"rebuilt {n_side} side-chain atom(s)"
    if n_term:
        detail += f", completed {n_term} terminal atom(s)"
    if gaps and not fill_loops:
        where = ", ".join(f"{len(v)} residue(s) after {k[1]}"
                          for k, v in list(gaps.items())[:4])
        detail += (f"; LEFT {sum(len(v) for v in gaps.values())} missing "
                   f"residue(s) alone ({where}); a rebuilt loop is a guess")
    elif gaps and fill_loops:
        detail += f"; built {sum(len(v) for v in gaps.values())} missing residue(s)"
    steps.append(Step("missing atoms", True, detail))
    return out_path


def minimise(path, out_path, steps, restraint=10.0, max_iterations=500,
             log=print):
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
        ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
        # The caps need their hydrogens, and the force field places
        # them: an NME is -NH-CH3, and placing those by hand would be
        # four further guesses.
        topology, positions = pdb.topology, pdb.positions
        if any(r.name in CAP_RESIDUES for r in topology.residues()):
            _m = app.Modeller(topology, positions)
            _m.addHydrogens(ff)
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
        heavy, free_caps = 0, 0
        for atom in topology.atoms():
            if atom.element is None or atom.element.symbol == "H":
                continue
            if atom.residue.name in CAP_RESIDUES:
                free_caps += 1          # placed by geometry; let it relax
                continue
            force.addParticle(atom.index, positions[atom.index])
            heavy += 1
        system.addForce(force)

        integrator = openmm.LangevinIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        sim = app.Simulation(topology, system, integrator,
                             openmm.Platform.getPlatformByName("CPU"))
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
           else "") + "; "
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
    per = {}
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
        per.setdefault(key, {})[nm] = xyz
    keys = sorted(per, key=lambda k: (k[0], k[1], k[2]))
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


def cap_breaks(path, out_path, gaps, steps, log=print):
    """ACE and NME on the ends a chain break creates. Returns out_path or None.

    A break is not a terminus. Splitting the chain so that pdb2pqr will
    protonate it leaves two artificial ends, and charging them as real termini
    puts a charge pair into the structure that the protein does not have.
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

    add = {}          # (chain, resnum, icode) -> [lines to put after it]
    pre = {}          # (chain, resnum, icode) -> [lines to put before it]
    n_cap = 0
    for a, b, _d in gaps:
        A, B = at.get(a, {}), at.get(b, {})
        # NME on the fragment that ENDS at `a`: N then its methyl.
        if {"CA", "C", "O"} <= set(A):
            d = -(_unit(A["CA"] - A["C"]) + _unit(A["O"] - A["C"]))
            npos = A["C"] + 1.329 * _unit(d)
            cpos = npos + 1.449 * _unit(_unit(npos - A["C"])
                                        + 0.5 * _unit(npos - A["O"]))
            num = free_number(a[0], a[1] + 1)
            add.setdefault(a, []).extend([
                rec(0, "N", "NME", a[0], num, npos, "N"),
                rec(0, "C", "NME", a[0], num, cpos, "C")])
            n_cap += 1
        # ACE on the fragment that BEGINS at `b`: carbonyl C, its O, methyl.
        if {"N", "CA"} <= set(B):
            d = -(_unit(B["CA"] - B["N"]))
            cpos = B["N"] + 1.335 * _unit(d)
            opos = cpos + 1.229 * _unit(_unit(cpos - B["N"])
                                        + 0.6 * _unit(B["CA"] - B["N"]))
            mpos = cpos + 1.508 * _unit(_unit(cpos - B["N"])
                                        - 0.9 * _unit(opos - cpos))
            num = free_number(b[0], max(1, b[1] - 1))
            pre.setdefault(b, []).extend([
                rec(0, "C", "ACE", b[0], num, cpos, "C"),
                rec(0, "O", "ACE", b[0], num, opos, "O"),
                rec(0, "CH3", "ACE", b[0], num, mpos, "C")])
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
        for nm in ("H2", "H3"):
            drop.add((b, nm))

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
    for leftover in (add, pre):
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
              minimise_structure, restraint):
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
    return {"receptor": os.path.abspath(pdb_in), "stamp": stamp,
            "ph": round(float(ph), 3), "forcefield": str(forcefield),
            "chain": str(chain), "waters": float(keep_waters_within),
            "hetatms": bool(keep_hetatms), "ions": bool(keep_ions),
            "loops": bool(fill_loops), "termini": bool(cap_termini),
            "minimise": bool(minimise_structure),
            "restraint": float(restraint)}


def run(pdb_in, workdir, ph=7.4, forcefield="AMBER", log=print, chain="",
        keep_waters_within=0.0, site_points=None, keep_hetatms=False,
        keep_ions=True, fill_loops=False, cap_termini=True,
        minimise_structure=True, restraint=10.0):
    """Prepare a receptor. Returns the prepared PDB, or None if nothing ran.

    The order matters and is the wizard's: tidy the file first, rebuild what
    is missing, then protonate, then minimise. Protonating before rebuilding
    would put hydrogens on atoms that are about to move; minimising before
    protonating would minimise a structure with no hydrogens on it.
    """
    os.makedirs(workdir, exist_ok=True)

    # Already prepared by an earlier attempt at this run? A job killed at
    # its wall clock limit is resumed, and the resumed job would
    # otherwise repeat the whole preparation.
    want = _prep_key(pdb_in, ph, forcefield, chain, keep_waters_within,
                     keep_hetatms, keep_ions, fill_loops, cap_termini,
                     minimise_structure, restraint)
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
    lines = drop_waters(lines, steps, keep_waters_within, site_points)
    if not keep_hetatms:
        lines = drop_hetatms(lines, steps, keep_ions)
    disulfides(lines, steps)

    tidy = _write(os.path.join(workdir, "tidied.pdb"), lines)
    current = tidy

    fixed = fix_structure(current, os.path.join(workdir, "rebuilt.pdb"), steps,
                          fill_loops=fill_loops, cap_termini=cap_termini,
                          log=log)
    current = fixed or current

    prot = protonate(current, os.path.join(workdir, "protonated.pdb"), steps,
                     ph=ph, forcefield=forcefield, log=log)
    current = prot or current

    # The caps go on here: after pdb2pqr, which has no parameters for ACE
    # or NME, and before OpenMM and tleap, which both have them.
    if cap_termini:
        _g = backbone_gaps(current)
        if _g:
            _capped = cap_breaks(current,
                                 os.path.join(workdir, "capped.pdb"), _g,
                                 steps, log=log)
            current = _capped or current

    if minimise_structure:
        mini = minimise(current, os.path.join(workdir, "minimised.pdb"), steps,
                        restraint=restraint, log=log)
        current = mini or current
    else:
        steps.append(Step("minimisation", False, "not asked for"))

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
    if n_h == 0:
        bar = "  " + "!" * 68
        log("")
        log(bar)
        log("  NO HYDROGENS ON THE PREPARED STRUCTURE")
        log("  The protonation line above says why. Docking from here uses")
        log("  hydrogens placed geometrically and IGNORES the run's pH, or")
        log("  fails outright in the grid step. Worth fixing before trusting")
        log("  any score this run produces.")
        log(bar)
        log("")

    with open(os.path.join(workdir, "preparation.log"), "w") as fh:
        fh.write("Protein preparation\n===================\n\n")
        fh.write(f"in  : {pdb_in}\nout : {out}\n"
                 f"pH  : {ph}   force field: {forcefield}\n\n")
        for s in steps:
            fh.write(str(s) + "\n")
        fh.write(f"\n{n_before} atoms in, {n_after} out\n")
    # Record what this preparation was made from, so the next attempt at this
    # run can tell "already done" from "done differently".
    try:
        with open(keyf, "w") as fh:
            json.dump(want, fh, indent=1)
    except OSError:
        pass                      # reuse is an optimisation, not the result
    return out


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
