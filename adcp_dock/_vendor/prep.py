"""Protein preparation. The steps live in prepare.py; this is the name callers use.

By Aadhil Haq

prep.run() hands the work to prepare, which follows the steps of a
preparation wizard: alternate locations, modified residues, waters,
cofactors, missing side chains, termini, disulfides, protonation at a
stated pH, and a restrained minimisation. run() and summarise() keep the
same interface for existing callers.
"""
from __future__ import annotations

from .prepare import (  # noqa: F401
    run, summarise, looks_like_trace, keep_chain, find_tool,
)


def available() -> str | None:
    """Path to pdb2pqr, or None.

    Only pdb2pqr is asked about here, because it is the one step whose absence
    changes the answer rather than the polish: without it there are no
    hydrogens and no titration states. PDBFixer and OpenMM missing means less
    is done, and prepare.py says which.
    """
    return find_tool("pdb2pqr", "pdb2pqr30")
