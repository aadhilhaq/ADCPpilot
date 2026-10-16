# Third-party notices

The BSD 3-Clause License in [LICENSE](LICENSE) covers the code in this repository
written for GUI to ADCP Docking. It does not cover, and cannot relicense, the
programs and libraries below. Each keeps its own licence.

## Programs run as separate executables (not included in this repository)

The package starts these as separate processes. Users obtain and install them
under their own terms.

| Program | Owner | Licence | How it is obtained |
|---|---|---|---|
| ADCP (adcp) | Scripps Research Institute | GNU LGPL (version 2.1 or later, per runADCP.py) | Downloaded from Scripps by the user's installer: in the ADCP 1.1 suite, and in ADFRsuite 1.0 where `--with-adfr10` adds it or it is already installed |
| agfr, prepare_receptor, AutoGrid4, AutoSite | Scripps Research Institute | ADFRsuite terms, mostly LGPL; see https://ccsb.scripps.edu/adfr/downloads/ and the MGLTOOLS row below, which covers packages these programs load | As above |
| MSMS program and library (mslib), and the UT packages | Scripps Research Institute; University of Texas at Austin | Free for academic use only; commercial use needs a licence from their owners. The MSMS terms read: "You are free to use the program (MSMS) ONLY FOR NON-COMMERCIAL SCIENTIFIC RESEARCH. This license is issued to you as an individual. For institution wide use or commercial use of the software you will have to follow the official registration procedure." They ask that published work using MSMS cite Sanner, Olson and Spehner, Biopolymers 38:305-320, 1996. | These are the two components ADFRsuite's `LICENSE.txt` marks academic-only. The **msms program** (`bin/msms`) is called by this package as the program installed with ADFRsuite, never copied into this repository or distributed with it: the interface's server runs it on request to triangulate the molecular surfaces the site view and the pose viewer draw (`adcp_dock/surface.py`). It is given the radii from MolKit's copy of `pdb_to_xyzr`'s table, read from the installed ADFRsuite as data. Where it is absent the viewers draw their own surfaces, so nothing else depends on it. The library (mslib) and the UT packages are imported by nothing this package runs: mslib is imported only by ADFR's map-fixing utility, its graphical interface, DejaVu2's Vision interface and the Pmv viewer, and the UT packages only by graphical and volume-reading modules. The installer now passes `-l`, but `-l` was added after the installed trees were built: on those trees `mslibACA` was renamed to `mslib`, and `bin/msms` ships and runs. |
| MolKit, mglutil, Support, DejaVu2, PyBabel and Volume, in ADFRsuite's `CCSBpckgs` | Scripps Research Institute (Dr Sanner) | MGLTOOLS Software License Agreement, carried by each of these packages and not mentioned in ADFRsuite's `LICENSE.txt`: "The programs received by you will be used only for NON COMMERCIAL purposes. This license is issued to you as an individual." Clause 3 adds: "You will not sell or otherwise distribute commercially these programs or derivatives to any other party, whether with or without consideration." Clause 2 counts scientific research in an academic environment and teaching as non-commercial; commercial use needs an arrangement with Dr Sanner at Scripps Research. | On the live run path, not merely unpacked. `prepare_receptor` imports MolKit directly and loads mglutil, Support and PyBabel with it; `agfr` loads Support, Volume and mglutil. Receptor preparation and map generation therefore run under these terms. The ADCP 1.1 suite carries its own copies of these packages, without licence files; where ADFRsuite 1.0 is not installed, its `prepare_receptor` and `agfr` load those copies, which should be treated under the same terms. AutoDockTools carries a licence file too, but it is the MSMS and MSLIB agreement, which states that the rest of AutoDockTools is freely available for commercial and academic use. |
| reduce 3.23 (`bin/reduce` and its `reduce_wwPDB_het_dict.txt`) | J. Michael Word and the Richardson laboratory, Duke University | Distributed with ADFRsuite under ADFRsuite's terms (above); reduce prints "Copyright 1997-2013, J. Michael Word" and asks that work using it cite Word et al., J Mol Biol 285:1735-1747, 1999 | Called as the program installed with ADFRsuite (the same one ADFR's own preparation and ADapt run), never copied: the default receptor protonation, and MM-GBSA's complexes in such runs |
| ADapt | Michel Sanner | Not yet licensed by its author; no licence has been chosen | A separate program, cloned by the user, started by the interface on request and reached through its `/adapt/` proxy. Nothing of it is in this repository |
| PDB2PQR | APBS/PDB2PQR developers | BSD-3-Clause | pip |
| PROPKA (run by PDB2PQR) | PROPKA developers | LGPL-2.1 | pip, with PDB2PQR |
| ADCP 1.1 environment (optional): adcp, adfr, autosite and their libraries | Scripps Research Institute | ADCP under GNU LGPL; the environment also includes the MSMS library (mslib) and the UT packages, free for academic use only | Installed by the user with `adcp-dock install-adcp11`, which follows Scripps' own script |
| AmberTools: tleap, sander, MMPBSA.py | Amber developers | Several licences (GPL-3.0-or-later, LGPL-3.0-or-later, BSD, MIT, by component) | conda, only when MM-GBSA is wanted |

## Python libraries imported by the package (not included in this repository)

| Library | Licence | Needed for |
|---|---|---|
| NumPy | BSD-3-Clause | always |
| RDKit | BSD-3-Clause | always |
| OpenMM | MIT and LGPL (by component) | receptor minimization (optional) |
| PDBFixer | MIT | receptor preparation (optional) |
| Biopython | Biopython License / BSD-3-Clause | reading mmCIF structures (optional) |
| paramiko | LGPL-2.1 | the desktop launcher only |

## Included in this repository

| Component | Version | Licence | Files |
|---|---|---|---|
| 3Dmol.js | 2.5.5 | BSD-3-Clause (includes code from GLmol, Three.js and jQuery) | `adcp_dock/static/3Dmol-min.js`, licence in `3Dmol.LICENSE.txt` and `3Dmol-min.js.LICENSE.txt` |
| Mol* viewer | 4.9.0 | MIT (bundles third-party code listed in its licence file) | `adcp_dock/static/molstar.js`, `molstar.css`, licence in `molstar.LICENSE.txt` and `molstar.js.LICENSE.txt` |

## Data

| Data | Source | Terms | File |
|---|---|---|---|
| Codes, coarse residues and structures of ADCP 1.1's non-standard residues | Derived by `scripts/modifications_table.py` from ADCP 1.1's side-chain library and OpenMM templates, which the Scripps Research Institute built from the SwissSidechain database (Gfeller, Michielin and Zoete, SIB Swiss Institute of Bioinformatics) | Scripps states that it has permission to distribute this library and that SwissSidechain's terms are moving to CC BY 4.0 (ADCP 1.1 documentation) | `adcp_dock/static/modifications.json` |
| Residue names | RCSB Chemical Component Dictionary; PubChem | CC0 1.0; public domain | As above |

## Desktop launcher builds

The Windows and macOS launcher downloads are built with PyInstaller and contain
a Python interpreter, Tcl/Tk, paramiko (LGPL-2.1) and its dependencies
(cryptography, bcrypt, PyNaCl), each under its own licence.
