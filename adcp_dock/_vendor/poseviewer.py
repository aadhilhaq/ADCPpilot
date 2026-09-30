"""Clean, viewer-ready pose files.

A docked pose cannot be handed to a molecular viewer as it is. Its atom
order follows the docking program's internal order rather than the residue
sequence, so a residue's atoms are not contiguous, and it carries no CONECT
records. Viewers then split one peptide into several ligands with duplicated
residues.

The fix is a coordinate transplant: build the peptide from its sequence so
the topology, atom order and bonds are correct, then copy the docked
coordinates onto it by matching (residue number, atom name). The geometry
does not change, only the representation.

Two writers are provided. The PDB writer is always available and writes
ordered atoms with explicit CONECT records. The Maestro writer is used when
$SCHRODINGER is set and writes a .maegz pose viewer.
"""
from __future__ import annotations
import os, subprocess, sys, tempfile
import numpy as np
from rdkit import Chem
from . import library, receptor, maewriter

CAP_RES = {"ACE", "NME", "NMA"}


def template(seq: str, capped: bool = True, ph: float | None = None) -> Chem.Mol | None:
    """Correctly built, correctly ordered peptide (heavy atoms only).

    `ph` sets the side-chain formal charges as the docking set them, so that a
    peptide written for a viewer carries the same charges it was docked with.
    """
    m = library.build_3d(seq, capped=capped, ph=ph)
    return Chem.RemoveHs(m) if m is not None else None


def with_hydrogens(mol: Chem.Mol) -> Chem.Mol:
    """The same molecule with hydrogens on it, placed geometrically.

    The transplant is heavy-atom only, because that is what a docked pose
    gives and what the matching needs. A viewer wants the hydrogens, and so
    does anyone checking that the charges are the ones the run used, since a
    protonated amine is three hydrogens rather than a plus sign somewhere.
    Returns the input unchanged if RDKit cannot place them.
    """
    if mol is None:
        return None
    try:
        # addResidueInfo, or the new hydrogens carry no PDBResidueInfo at all
        # and every writer downstream falls back to its "unknown" branch: 32
        # of the 61 atoms land in a residue called LIG, numbered 1, and the
        # peptide a viewer draws has its hydrogens in a separate molecule.
        return Chem.AddHs(mol, addCoords=True, addResidueInfo=True)
    except Exception:                                          # noqa: BLE001
        return mol


def transplant(seq: str, pose: dict, capped: bool | None = None,
               ph: float | None = None) -> Chem.Mol | None:
    """Docked coordinates on a clean topology, matched by (resnum, atom name).

    Vina poses come from the capped library ligand. ADCP writes uncapped peptides
    numbered from 1, so the template must match; when `capped` is None it is
    inferred from whether the pose contains ACE/NME residues.
    """
    if capped is None:
        capped = any(str(r).strip().upper() in CAP_RES for r in pose.get("resname", []))
    tmpl = template(seq, capped=capped, ph=ph)
    if tmpl is None:
        return None
    coords = {}
    for i, nm in enumerate(pose["name"]):
        coords[(int(pose["resnum"][i]), nm.strip())] = pose["xyz"][i]
    conf = tmpl.GetConformer()
    missing = []
    for at in tmpl.GetAtoms():
        info = at.GetMonomerInfo()
        key = (info.GetResidueNumber(), info.GetName().strip())
        if key in coords:
            x, y, z = coords[key]
            conf.SetAtomPosition(at.GetIdx(), (float(x), float(y), float(z)))
        else:
            missing.append(at.GetIdx())
    if missing:
        # A pose that cannot be matched atom-for-atom is not trustworthy; refuse
        # it rather than emit a peptide with atoms left at template positions.
        # The one tolerated exception is a terminal OXT, which several programs
        # (ADCP among them) simply do not write.
        names = {tmpl.GetAtomWithIdx(i).GetMonomerInfo().GetName().strip()
                 for i in missing}
        if names != {"OXT"}:
            return None
        ed = Chem.RWMol(tmpl)
        # The carbon that loses the OXT must not gain a hydrogen: the
        # C-terminal carbon is written C(=O)O in the SMILES and
        # carries noImplicit=False, so removing its OXT would add an
        # implicit hydrogen and change the molecule.
        for i in sorted(missing, reverse=True):
            for nb in ed.GetAtomWithIdx(i).GetNeighbors():
                nb.SetNoImplicit(True)
                nb.SetNumExplicitHs(0)
            ed.RemoveAtom(i)
        tmpl = ed.GetMol()
        Chem.SanitizeMol(tmpl)
    return tmpl


def _pdb_records(mol: Chem.Mol, chain: str, start_serial: int):
    """Ordered ATOM records plus CONECT records for one molecule."""
    lines, serial = [], start_serial
    idx2serial = {}
    for at in mol.GetAtoms():
        info = at.GetMonomerInfo()
        pos = mol.GetConformer().GetAtomPosition(at.GetIdx())
        a = receptor.Atom(info.GetResidueNumber(), info.GetResidueName().strip(),
                          info.GetName().strip(), at.GetSymbol(),
                          np.array([pos.x, pos.y, pos.z]))
        lines.append(receptor.pdb_line(serial, a, chain))
        idx2serial[at.GetIdx()] = serial
        serial += 1
    conect = []
    for b in mol.GetBonds():
        i, j = idx2serial[b.GetBeginAtomIdx()], idx2serial[b.GetEndAtomIdx()]
        conect.append(f"CONECT{i:5d}{j:5d}\n")
        conect.append(f"CONECT{j:5d}{i:5d}\n")
    return lines, conect, serial


def write_pdb(entries: list[tuple[str, Chem.Mol]], path: str,
              rec=None, chain: str = "L") -> int:
    """One MODEL per pose. CONECT records go inside each MODEL so the bonding is
    unambiguous when several models share serial numbers."""
    n = 0
    with open(path, "w") as fh:
        for model, (title, mol) in enumerate(entries, 1):
            if mol is None:
                continue
            fh.write(f"MODEL     {model:4d}\n")
            fh.write(f"REMARK   1 TITLE {title}\n")
            serial = 1
            if rec is not None:
                for a in rec.atoms:
                    fh.write(receptor.pdb_line(serial, a, "A"))
                    serial += 1
                fh.write("TER\n")
            lines, conect, serial = _pdb_records(mol, chain, serial)
            fh.writelines(lines)
            fh.write("TER\n")
            fh.writelines(conect)
            fh.write("ENDMDL\n")
            n += 1
        fh.write("END\n")
    return n


def write_sdf(entries: list[tuple[str, Chem.Mol]], path: str) -> int:
    """Explicit bond blocks; useful for cheminformatics tools."""
    w = Chem.SDWriter(path)
    n = 0
    for title, mol in entries:
        if mol is None:
            continue
        mol = Chem.Mol(mol)
        mol.SetProp("_Name", title)
        w.write(mol)
        n += 1
    w.close()
    return n


def write_maegz(entries, path: str, rec_pdb: str | None = None,
                props: dict | None = None) -> bool:
    """Maestro pose viewer, if a Schrodinger install can be found. Optional.

    Found rather than demanded: requiring $SCHRODINGER meant this quietly
    wrote nothing inside a batch job, where nothing has loaded the module,
    on clusters that have the whole suite installed.
    """
    from . import gbengines
    schro = gbengines.find_schrodinger()
    if not schro or not os.path.isfile(os.path.join(schro, "run")):
        return False
    with tempfile.TemporaryDirectory() as td:
        sdf = os.path.join(td, "ligands.sdf")
        if not write_sdf(entries, sdf):
            return False
        script = os.path.join(td, "conv.py")
        with open(script, "w") as fh:
            fh.write(
                "import sys\n"
                "from schrodinger import structure\n"
                "sdf, out, rec = sys.argv[1], sys.argv[2], sys.argv[3]\n"
                "w = structure.StructureWriter(out)\n"
                "if rec != 'NONE':\n"
                "    r = next(structure.StructureReader(rec))\n"
                "    r.title = 'receptor'\n"
                "    w.append(r)\n"
                "for st in structure.StructureReader(sdf):\n"
                "    w.append(st)\n"
                "w.close()\n")
        cmd = [os.path.join(schro, "run"), "python3", script, sdf, path,
               rec_pdb or "NONE"]
        subprocess.run(cmd, capture_output=True, text=True)
    return os.path.isfile(path)


def receptor_mol(pdb_path: str) -> Chem.Mol | None:
    """Receptor as an RDKit molecule with residue information intact, so it can
    be written as the first entry of a pose viewer."""
    m = Chem.MolFromPDBFile(pdb_path, removeHs=False, sanitize=False)
    if m is None:
        return None
    try:
        Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL ^
                         Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    except Exception:
        pass
    return m


def write_mae(entries, path: str, rec_pdb: str | None = None) -> int:
    """Maestro pose viewer, receptor included as the first entry.

    Written directly rather than converted, so no Schrodinger licence is needed.
    A PDB carries residues but no entry title Maestro will read; an SDF carries a
    title but no residues. The .mae carries both, plus the scores as properties.
    """
    rec = receptor_mol(rec_pdb) if rec_pdb else None
    return maewriter.write(entries, path, receptor_mol=rec)
