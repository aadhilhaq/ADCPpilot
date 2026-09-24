# Third-party notices

The MIT License in [LICENSE](LICENSE) covers the code in this repository
written for GUI to ADCP Docking. It does not cover, and cannot relicense, the
programs and libraries below. Each keeps its own licence.

## Programs run as separate executables (not included in this repository)

The package starts these as separate processes. Users obtain and install them
under their own terms.

| Program | Owner | Licence | How it is obtained |
|---|---|---|---|
| ADCP (adcp) | Scripps Research Institute | GNU LGPL (version 2.1 or later, per runADCP.py) | Downloaded from Scripps by the user's installer, in ADFRsuite 1.0 |
| agfr, prepare_receptor, AutoGrid4, AutoSite | Scripps Research Institute | ADFRsuite terms, generally LGPL; see https://ccsb.scripps.edu/adfr/downloads/ | As above |
| MSMS, MSLIB (in the ADFRsuite download) | Scripps Research Institute | Free for academic use; commercial use needs a licence | Not used and not installed: the installers pass the ADFRsuite installer's `-l` switch, which leaves them out |
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

## Desktop launcher builds

The Windows and macOS launcher downloads are built with PyInstaller and contain
a Python interpreter, Tcl/Tk, paramiko (LGPL-2.1) and its dependencies
(cryptography, bcrypt, PyNaCl), each under its own licence.
