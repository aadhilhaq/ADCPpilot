"""Write Maestro .mae files directly.

By Aadhil Haq

Neither PDB nor SDF can carry everything a pose viewer needs:

  PDB   has residues and atom names, but no per-entry title that Maestro reads
        (REMARK lines are ignored), and no place for scores.
  SDF   carries a title and properties, but has no concept of residues, so every
        peptide arrives as UNK 900 with no sequence.

The .mae format carries all of it (title, residue and atom names, chains, bonds
and arbitrary typed properties) and is plain text, so writing it needs no
Schrodinger licence.
"""
from __future__ import annotations
from rdkit import Chem

_ORDER = {Chem.BondType.SINGLE: 1, Chem.BondType.DOUBLE: 2,
          Chem.BondType.TRIPLE: 3, Chem.BondType.AROMATIC: 1}

HEADER = ('{\n  s_m_m2io_version\n  :::\n  2.0.0\n}\n\n')


def _q(s: str) -> str:
    return '"' + str(s).replace('"', "'") + '"'


def _prop_block(props: dict) -> tuple[list[str], list[str]]:
    """Maestro types properties by name prefix: s_ string, r_ real, i_ integer."""
    names, values = [], []
    for k, v in props.items():
        if isinstance(v, bool):
            continue
        if isinstance(v, int):
            key = k if k.startswith("i_") else f"i_user_{k}"
            names.append(key); values.append(str(v))
        elif isinstance(v, float):
            key = k if k.startswith("r_") else f"r_user_{k}"
            names.append(key); values.append(f"{v:.4f}")
        else:
            key = k if k.startswith("s_") else f"s_user_{k}"
            names.append(key); values.append(_q(v))
    return names, values


def ct_block(mol: Chem.Mol, title: str, props: dict | None = None,
             default_chain: str = "L") -> str:
    """One structure (an f_m_ct block)."""
    conf = mol.GetConformer()
    props = dict(props or {})
    names, values = _prop_block(props)

    out = ["f_m_ct {", "  s_m_title"]
    out += [f"  {n}" for n in names]
    out += ["  :::", f"  {_q(title)}"]
    out += [f"  {v}" for v in values]

    out.append(f"  m_atom[{mol.GetNumAtoms()}] {{")
    out += ["    # First column is atom index #",
            "    r_m_x_coord", "    r_m_y_coord", "    r_m_z_coord",
            "    i_m_residue_number", "    s_m_pdb_residue_name",
            "    s_m_pdb_atom_name", "    i_m_atomic_number",
            # i_m_formal_charge, or Maestro reads every atom as neutral.
            #
            # It was harmless while the file held heavy atoms only: a lysine
            # NZ with one bond and no charge is an ordinary neutral amine to
            # look at. With the hydrogens now written, that same NZ has FOUR
            # bonds and no charge: a tetravalent neutral nitrogen, which is
            # not a thing, and Maestro says so on load.
            "    i_m_formal_charge",
            "    s_m_chain_name", "    :::"]
    for i, at in enumerate(mol.GetAtoms(), 1):
        p = conf.GetAtomPosition(at.GetIdx())
        info = at.GetMonomerInfo()
        if info is not None:
            resnum = info.GetResidueNumber()
            resname = info.GetResidueName().strip()[:4]
            aname = info.GetName().strip()[:4]
            chain = (info.GetChainId() or default_chain).strip() or default_chain
        else:
            resnum, resname, aname, chain = 1, "UNK", at.GetSymbol(), default_chain
        # Maestro expects the PDB atom-name field padded to four characters
        aname = aname.ljust(4) if len(aname) >= 4 else (" " + aname).ljust(4)
        out.append(f"    {i} {p.x:.6f} {p.y:.6f} {p.z:.6f} {resnum} "
                   f"{_q(resname.ljust(4))} {_q(aname)} {at.GetAtomicNum()} "
                   f"{at.GetFormalCharge()} {_q(chain)}")
    out += ["    :::", "  }"]

    bonds = list(mol.GetBonds())
    if bonds:
        out.append(f"  m_bond[{len(bonds)}] {{")
        out += ["    # First column is bond index #",
                "    i_m_from", "    i_m_to", "    i_m_order", "    :::"]
        for j, b in enumerate(bonds, 1):
            out.append(f"    {j} {b.GetBeginAtomIdx() + 1} {b.GetEndAtomIdx() + 1} "
                       f"{_ORDER.get(b.GetBondType(), 1)}")
        out += ["    :::", "  }"]
    out.append("}")
    return "\n".join(out) + "\n\n"


def write(entries: list[tuple[str, Chem.Mol, dict]], path: str,
          receptor_mol: Chem.Mol | None = None,
          receptor_title: str = "receptor") -> int:
    """Write a pose viewer. The receptor, when given, is the first entry so the
    file opens as a complex rather than needing a separate import."""
    n = 0
    with open(path, "w") as fh:
        fh.write(HEADER)
        if receptor_mol is not None:
            fh.write(ct_block(receptor_mol, receptor_title, {}, default_chain="A"))
            n += 1
        for title, mol, props in entries:
            if mol is None:
                continue
            fh.write(ct_block(mol, title, props))
            n += 1
    return n
