# Third-party notices

The BSD 3-Clause License in [LICENSE](LICENSE) covers the code in this repository
written for GUI to ADCP Docking. It does not cover, and cannot relicense, the
programs and libraries below. Each keeps its own licence.

## Programs run as separate executables (not included in this repository)

The package starts these as separate processes. Users obtain and install them
under their own terms.

| Program | Owner | Licence | How it is obtained |
|---|---|---|---|
| ADCP (adcp) and AGFR (agfr, with autogrid4), the ADCP 1.1 suite's `adcp` and `adfr` packages | Scripps Research Institute (Dr Sanner's laboratory) | Their sources (`runADCP.py`, `runAGFR.py`) state the GNU LGPL, version 2.1 or later; the packages' metadata names the MIT licence, and their `LICENSE` files hold the MIT text | Installed by the user's installer (`adcp_dock/adcp11_install.py`, which `install.sh`, `adcp-dock local` and the Colab notebook run), following Scripps' own script: each Scripps package is downloaded from ccsb.scripps.edu and checked against a pinned SHA-256 |
| AutoSite, and the suite's libraries agfr loads (ADFRcc, bhtree, mglutil, MolKit2, ProDy's Scripps build, Support, Volume) | sannerlab | BSD 2-Clause, "Copyright (c) 2022, sannerlab", in each package's `LICENSE` | As above. agfr loads these and none of the packages in the next row |
| prepare_receptor (AutoDockTools), and the MolKit and PyBabel it loads | Michel F. Sanner and the Scripps Research Institute | MGLTOOLS Software License Agreement, in each of these packages' `LICENSE`: "The programs received by you will be used only for NON COMMERCIAL purposes. This license is issued to you as an individual." Clause 2 counts scientific research in an academic environment and teaching as non-commercial; clause 3 adds: "You will not sell or otherwise distribute commercially these programs or derivatives to any other party, whether with or without consideration." Commercial use needs an arrangement with Dr Sanner at Scripps Research. | On the live run path: every receptor this package prepares is written as PDBQT by prepare_receptor, which imports AutoDockTools, MolKit and PyBabel (with bhtree and mglutil, BSD 2-Clause). A receptor from ADapt or from a target file is not |
| The MSMS library (mslib, mslibCOM) | sannerlab | BSD 2-Clause, "Copyright (c) 2022, sannerlab", in the package's `LICENSE`. Published work that shows an MSMS surface cites Sanner, Olson and Spehner, Biopolymers 38:305-320, 1996 | Run by the interface's server on request, with the suite's own Python through `adcp_dock/msms_mslib.py`, to triangulate the surfaces the site view and the pose viewer draw (`adcp_dock/surface.py`); never copied into this repository. The atoms' radii come from MolKit's copy of `pdb_to_xyzr`'s table (`radii_patterns.py`), read from the installed suite as data, and MolKit is under the MGLTOOLS licence (the row above). Where the suite is absent the viewers draw their own surfaces, so nothing else depends on it |
| reduce 3.23 (`bin/reduce` and its `reduce_wwPDB_het_dict.txt`) | J. Michael Word and the Richardson laboratory, Duke University | Shipped inside the ADCP 1.1 suite's `adfr` package (the first row); reduce prints "Copyright 1997-2013, J. Michael Word" and asks that work using it cite Word et al., J Mol Biol 285:1735-1747, 1999 | Called as the program installed with the suite (the same one ADFR's own preparation and ADapt run), never copied: the default receptor protonation, and MM-GBSA's complexes in such runs |
| The suite's other Scripps packages (the UT packages, DejaVu2, opengltk, PmvApp and the rest of the viewers' libraries) | sannerlab | BSD 2-Clause, "Copyright (c) 2022, sannerlab", in each package's `LICENSE` | Installed with the suite; nothing this package runs imports them |
| ADapt | Michel Sanner | Not yet licensed by its author; no licence has been chosen | A separate program, cloned by the user, started by the interface on request and reached through its `/adapt/` proxy. Nothing of it is in this repository |
| PDB2PQR | APBS/PDB2PQR developers | BSD-3-Clause | pip |
| PROPKA (run by PDB2PQR) | PROPKA developers | LGPL-2.1 | pip, with PDB2PQR |
| The ADCP 1.1 suite's environment: Python 3.7, OpenMM, Open Babel, PDBFixer, ParmEd and NumPy from conda-forge; Pillow, PySide2, pybel, Biopython and pypdb from PyPI | Their developers | Each under its own licence | Installed with the suite, as Scripps' own script installs them |
| AmberTools: tleap, sander, MMPBSA.py | Amber developers | Several licences (GPL-3.0-or-later, LGPL-3.0-or-later, BSD, MIT, by component) | conda, only when MM-GBSA is wanted |

### On a Mac only: ADFRsuite 1.0

The Scripps Research Institute builds no ADCP 1.1 for macOS, so on macOS, and
nowhere else, `adcp-dock local` installs ADFRsuite 1.0 and docks with it
(`adcp_dock/mac10.py`). It downloads `ADFRsuite_x86_64Darwin_1.0.tar.gz` from
ccsb.scripps.edu, checks it against a pinned SHA-256 and runs Scripps' own
`install.sh`. Nothing of it is in this repository.

| Program | Owner | Licence | How it is used |
|---|---|---|---|
| ADCP 1.0 (adcp), agfr with AutoGrid4 and AutoSite 1.0, and prepare_receptor | Scripps Research Institute | ADFRsuite's `LICENSE.txt`: "Unless otherwise specified, this software is distributed under the LGPL 2 open source license"; `runADCP.py` states the GNU LGPL, version 2.1 or later | Run as installed: they dock, build the maps and write the receptor as PDBQT |
| MolKit, mglutil, Support, DejaVu2, PyBabel and Volume, in its `CCSBpckgs` | Scripps Research Institute (Dr Sanner) | The MGLTOOLS Software License Agreement, in each of these packages' `LICENSE`, quoted in the table above: non-commercial use only | On the run path: `prepare_receptor` imports MolKit and loads mglutil, Support and PyBabel; agfr loads Support, Volume and mglutil. AutoDockTools' own `LICENSE` is the MSMS and MSLIB agreement |
| The MSMS program and library (`bin/msms`, mslibACA, mslibCOM), and the UT packages | Scripps Research Institute; University of Texas at Austin | Free for academic use only, as ADFRsuite's `LICENSE.txt` marks them; commercial use needs a licence from their owners. Published work that shows an MSMS surface cites Sanner, Olson and Spehner, Biopolymers 38:305-320, 1996 | The installer passes `-l`, which leaves Scripps' academic-or-commercial question about the library unanswered, so neither library is put in place and nothing this package runs imports it or the UT packages. The server would run the msms program, as installed, for the viewers' surfaces, but the macOS build is 32-bit, which macOS 10.15 and later do not run, so the viewers draw their own |
| reduce (`bin/reduce` and its `reduce_wwPDB_het_dict.txt`) | J. Michael Word and the Richardson laboratory, Duke University | Distributed with ADFRsuite under its terms; reduce asks that work using it cite Word et al., J Mol Biol 285:1735-1747, 1999 | Run as installed: the default receptor protonation on a Mac |
| ADFRsuite's Python 2.7, Open Babel, NumPy, ProDy and its other bundled libraries | Their developers | Each under its own licence | Installed with ADFRsuite; its programs run on them |

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
