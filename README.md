<p align="center">
  <img src="brand/banner.png" alt="GUI to ADCP Docking" width="100%">
</p>

GUI to ADCP Docking is a browser interface and command line for peptide docking with [AutoDock CrankPep (ADCP)](https://ccsb.scripps.edu/adcp/). Give it a receptor structure and peptide sequences: it prepares the receptor, builds the grid maps, docks each peptide, can rescore the poses with MM-GBSA, and reports ranked results with the structures. The calculations run on a Slurm cluster, on your own computer, or in Google Colab.

Developed by Aadhil Haq, HOPE Lab, under the supervision of Dr. Sandun Fernando. ADCP and ADFRsuite are developed by the Scripps Research Institute; this project is independent of it.

## Features

* Receptor preparation: protonation at a chosen pH, missing side chains, a restrained minimization, and a docking box built around the residues you choose.
* ADCP 1.0 and ADCP 1.1. With ADCP 1.1: non-standard and D amino acids, chosen from a peptide's menu with a drawing of each; OpenMM minimization of the top poses; runs that start from an AGFR or AGFRGUI target file.
* Optional MM-GBSA rescoring with AmberTools.
* A ranked results table, a 3D pose viewer, and the citations for the tools each run used.
* An interface panel beside the scores: the four MM-GBSA energy terms, buried surface area and how much of the peptide is buried, typed contacts and clashes, per residue energies from Amber's own decomposition, how much of the nominated site the pose touches, and the sequence liabilities. Every number is measured on the geometry its energy was computed on, and the run records the cutoffs and radii it used.
* On Slurm, one array task per peptide; a run that stopped can be resumed.
* A desktop launcher for Windows and macOS that signs in to a cluster and opens the interface.

## Getting started

| Where the calculations run | How to install |
|---|---|
| A Slurm cluster | [INSTALL.md, section 1](INSTALL.md#1-on-a-slurm-cluster) |
| Your own computer (Linux, macOS, or Windows through WSL) | [INSTALL.md, section 2](INSTALL.md#2-on-your-own-computer) |
| Google Colab | [Open the notebook](https://colab.research.google.com/github/aadhilhaq/ADCP_docking/blob/main/notebooks/ADCP_docking_colab.ipynb); nothing is installed on your computer |

On a cluster, once the environment in INSTALL.md is set up, in a group or scratch folder rather than your home directory:

```bash
git clone https://github.com/aadhilhaq/ADCP_docking.git
cd ADCP_docking
pip install --no-cache-dir -e ".[local]"
python -m adcp_dock.adfr --dest tools/docking/ADFRsuite
export ADFRSUITE_HOME=$PWD/tools/docking/ADFRsuite
adcp-dock check
adcp-dock serve --port 0
```

The server prints an address with an access token. Open it through an SSH tunnel ([INSTALL.md, 1.5](INSTALL.md#15-start-the-page)) or with the [desktop launcher](INSTALL.md#5-the-desktop-launcher).

On your own computer:

```bash
python3 -m venv ~/.adcp-dock/venv
~/.adcp-dock/venv/bin/pip install "adcp-dock[local] @ git+https://github.com/aadhilhaq/ADCP_docking"
~/.adcp-dock/venv/bin/adcp-dock local
```

ADCP 1.1 is an optional extra install: `adcp-dock install-adcp11` ([INSTALL.md, section 4](INSTALL.md#4-adcp-11-optional)).

`adcp-dock example` adds a finished run, MDM2 with three peptides, to open before docking anything. The Colab notebook includes it.

## Command line

```bash
adcp-dock new runs/trial -r 1YCR --chain A \
    -s A:54,A:62,A:67,A:96,A:99 -p "p53: ETFSDLWKLLP" --poses 10 --mmgbsa
adcp-dock run runs/trial
adcp-dock status runs/trial
adcp-dock report runs/trial
adcp-dock panel runs/trial      # only for a run docked before 0.1.17
```

`-r` takes a structure file or a PDB ID, `-s` the residues the docking box covers, `-p` a peptide (repeat it for more) or `-f` a peptide file. Running `adcp-dock run` again docks only the peptides that have no result yet.

A run docked before version 0.1.17 has no interface panel in its results, because the panel is measured while each peptide is rescored. The run's page shows its MM-GBSA breakdown straight away, read from those files, and offers a **Measure them** button for the interface columns; `adcp-dock panel` does the same from the command line. Either one measures the panel from the files the MM-GBSA already left, submits that to Slurm, and rebuilds the report when it finishes. Nothing is docked or rescored again; the only new calculation is the per residue decomposition, which takes seconds per pose. A run whose intermediate files were deleted keeps its four energy terms and gets no interface columns, and the report says so.

## Peptides

Sequences of 5 to 30 residues in one-letter codes, one per line, optionally named (`p53: ETFSDLWKLLP`) or as FASTA. Lines starting with `#` are comments. With ADCP 1.1, write a non-standard residue as its coarse one-letter code followed by its library name in angle brackets, and a D residue with `&` before its letter (`N<BHD>PIS&DVD`), or choose them in the interface.

## Results

Each run writes `results.tsv`, `results.json` and `results.txt`, the ranked poses, a `metrics/` folder, and its logs. Peptides are ranked by the mean of up to five top poses: by MM-GBSA where it ran, otherwise by the ADCP score. Both are ranking aids, not binding free energies. Look at the poses, and make sure the search settings suit your system. The logs record every preparation decision and any step that could not run.

`metrics/` holds the interface panel: one file per peptide with every pose that was measured, and one residue map listing every residue of the peptide and each receptor residue at the interface, what kind of contact it makes, and its own share of the energy. The same map is shown under the 3D viewer when a peptide is opened. Every energy in the results and in the panel comes from AmberTools: tleap builds the topologies, sander minimizes the complex with the receptor restrained, and MMPBSA.py computes the energy and its per residue decomposition. This tool also runs OpenMM, but only in two places that produce no reported energy: the receptor's restrained minimization during preparation, and ADCP 1.1's own minimization of its top poses, whose energies are printed beside the scores and never ranked by. `results.json` names which of those a given run used. `results.json` carries the constants the panel was measured with, so a number in it can be compared with a number from somewhere else. Two limits are worth knowing before the columns are quoted: the contact counts test the distance between heavy atoms only, not hydrogen positions or ring orientation; and a per residue energy comes from a decomposition that agrees with the reported total to within a few per cent rather than exactly, and always slightly more negative, because Amber uses a per atom surface term in a decomposition where the total uses LCPO. The energy terms themselves agree to within 0.05 kcal/mol.

## Requirements

Python 3.9 or newer with NumPy and RDKit, ADFRsuite (installed by `python -m adcp_dock.adfr`, or on the first `adcp-dock local`), and PDB2PQR (installed with the package). Optional: AmberTools for MM-GBSA; OpenMM and PDBFixer for receptor repair and minimization; Biopython for mmCIF files; ADCP 1.1. `adcp-dock check` lists what it finds.

## Citation

This interface drives other programs, and each has its own paper. Cite the ones your run used; `adcp-dock cite <run>` and the `HOW_TO_CITE.txt`, `CITATION.bib` and `CITATION.ris` files written into every run list exactly those. ADCP and ADFRsuite are products of the Scripps Research Institute. The Amber developers ask for the AmberTools paper and the force field paper ([ambermd.org/CiteAmber.php](https://ambermd.org/CiteAmber.php)).

| Tool | Cite for | Reference | DOI |
|---|---|---|---|
| ADCP | every run | Zhang Y, Sanner MF. AutoDock CrankPep: combining folding and docking to predict protein-peptide complexes. Bioinformatics 35(24):5121-5127, 2019. | [10.1093/bioinformatics/btz459](https://doi.org/10.1093/bioinformatics/btz459) |
| ADFRsuite (agfr, prepare_receptor) | every run | Ravindranath PA, Forli S, Goodsell DS, Olson AJ, Sanner MF. AutoDockFR: advances in protein-ligand docking with explicitly specified binding site flexibility. PLoS Comput Biol 11(12):e1004586, 2015. | [10.1371/journal.pcbi.1004586](https://doi.org/10.1371/journal.pcbi.1004586) |
| AutoGrid4 | every run | Morris GM, Huey R, Lindstrom W, Sanner MF, Belew RK, Goodsell DS, Olson AJ. AutoDock4 and AutoDockTools4: automated docking with selective receptor flexibility. J Comput Chem 30(16):2785-2791, 2009. | [10.1002/jcc.21256](https://doi.org/10.1002/jcc.21256) |
| AutoSite | every run | Ravindranath PA, Sanner MF. AutoSite: an automated approach for pseudo-ligands prediction—from ligand-binding sites identification to predicting key ligand atoms. Bioinformatics 32(20):3142-3149, 2016. | [10.1093/bioinformatics/btw367](https://doi.org/10.1093/bioinformatics/btw367) |
| PDB2PQR | receptor preparation, and MM-GBSA | Jurrus E, Engel D, Star K, et al. Improvements to the APBS biomolecular solvation software suite. Protein Sci 27(1):112-128, 2018. | [10.1002/pro.3280](https://doi.org/10.1002/pro.3280) |
| PROPKA | receptor preparation, and MM-GBSA | Olsson MHM, Sondergaard CR, Rostkowski M, Jensen JH. PROPKA3: consistent treatment of internal and surface residues in empirical pKa predictions. J Chem Theory Comput 7(2):525-537, 2011. | [10.1021/ct100578z](https://doi.org/10.1021/ct100578z) |
| PDBFixer | receptor preparation (default) | Eastman P, et al. OpenMM 7: rapid development of high performance algorithms for molecular dynamics. PLoS Comput Biol 13(7):e1005659, 2017.  (PDBFixer is distributed with OpenMM.) | [10.1371/journal.pcbi.1005659](https://doi.org/10.1371/journal.pcbi.1005659) |
| OpenMM | receptor minimization (default) | Eastman P, et al. OpenMM 7: rapid development of high performance algorithms for molecular dynamics. PLoS Comput Biol 13(7):e1005659, 2017. | [10.1371/journal.pcbi.1005659](https://doi.org/10.1371/journal.pcbi.1005659) |
| AmberTools (tleap, sander) | runs with MM-GBSA | Case DA, et al. AmberTools. J Chem Inf Model 63(20):6183-6191, 2023. | [10.1021/acs.jcim.3c01153](https://doi.org/10.1021/acs.jcim.3c01153) |
| ff14SB force field | receptor minimization (default), and MM-GBSA | Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, Simmerling C. ff14SB: improving the accuracy of protein side chain and backbone parameters from ff99SB. J Chem Theory Comput 11(8):3696-3713, 2015. | [10.1021/acs.jctc.5b00255](https://doi.org/10.1021/acs.jctc.5b00255) |
| MMPBSA.py | runs with MM-GBSA | Miller BR III, McGee TD Jr, Swails JM, Homeyer N, Gohlke H, Roitberg AE. MMPBSA.py: an efficient program for end-state free energy calculations. J Chem Theory Comput 8(9):3314-3321, 2012. | [10.1021/ct300418h](https://doi.org/10.1021/ct300418h) |
| GB-OBC(II) solvation model, with its mbondi2 radii | runs with MM-GBSA (the default model) | Onufriev A, Bashford D, Case DA. Exploring protein native states and large-scale conformational changes with a modified generalized born model. Proteins 55(2):383-394, 2004. | [10.1002/prot.20033](https://doi.org/10.1002/prot.20033) |
| SwissSidechain | runs with non-standard residues (ADCP 1.1) | Gfeller D, Michielin O, Zoete V. SwissSidechain: a molecular and structural database of non-natural sidechains. Nucleic Acids Res 41(D1):D327-D332, 2013. | [10.1093/nar/gks991](https://doi.org/10.1093/nar/gks991) |
| GAFF | OpenMM minimization of peptides with non-standard residues | Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and testing of a general amber force field. J Comput Chem 25(9):1157-1174, 2004. | [10.1002/jcc.20035](https://doi.org/10.1002/jcc.20035) |
| Amber ff99SB force field | ADCP 1.1 poses minimized with OpenMM | Hornak V, Abel R, Okur A, Strockbine B, Roitberg A, Simmerling C. Comparison of multiple Amber force fields and development of improved protein backbone parameters. Proteins 65(3):712-725, 2006. | [10.1002/prot.21123](https://doi.org/10.1002/prot.21123) |
| GBn2 implicit solvent model | receptor minimization (default), ADCP 1.1 poses minimized in implicit solvent, and MM-GBSA run with GBn2 | Nguyen H, Roe DR, Simmerling C. Improved generalized Born solvent model parameters for protein simulations. J Chem Theory Comput 9(4):2020-2034, 2013. | [10.1021/ct3010485](https://doi.org/10.1021/ct3010485) |
| 3Dmol.js | figures made in the 3Dmol.js viewer | Rego N, Koes D. 3Dmol.js: molecular visualization with WebGL. Bioinformatics 31(8):1322-1324, 2015. | [10.1093/bioinformatics/btu829](https://doi.org/10.1093/bioinformatics/btu829) |
| Mol* | figures made in the Mol* viewer | Sehnal D, Bittrich S, Deshpande M, et al. Mol* Viewer: modern web app for 3D visualization and analysis of large biomolecular structures. Nucleic Acids Res 49(W1):W431-W437, 2021. | [10.1093/nar/gkab314](https://doi.org/10.1093/nar/gkab314) |
| NumPy | every run | Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with NumPy. Nature 585(7825):357-362, 2020. | [10.1038/s41586-020-2649-2](https://doi.org/10.1038/s41586-020-2649-2) |
| RDKit | every run | RDKit: Open-source cheminformatics. https://www.rdkit.org | [10.5281/zenodo.591637](https://doi.org/10.5281/zenodo.591637) |

```bash
adcp-dock cite                 # all supported tools
adcp-dock cite <run>           # tools recorded for one run
adcp-dock cite --bibtex        # use --ris for EndNote or Zotero
```

For this interface, use GitHub's **Cite this repository** button or [CITATION.cff](CITATION.cff), and record the software version used.

## Support

Report problems through [GitHub Issues](https://github.com/aadhilhaq/ADCP_docking/issues) or **Send feedback** in the interface, with the version, the settings and the relevant log. To run the tests, see [CONTRIBUTING.md](CONTRIBUTING.md).

## License

| What | License | Text |
|---|---|---|
| This package: all code and documents in this repository, except the rows below | BSD 3-Clause | [LICENSE](LICENSE) |
| 3Dmol.js, bundled structure viewer | BSD 3-Clause | `adcp_dock/static/3Dmol.LICENSE.txt` |
| Mol*, bundled structure viewer | MIT | `adcp_dock/static/molstar.LICENSE.txt` |
| Table of ADCP 1.1's non-standard residues (codes, structures) | Derived from ADCP 1.1's SwissSidechain library | [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) |
| ADCP, ADFRsuite, AmberTools, PDB2PQR and the other programs it runs | Their own licences; not included in this repository | [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) |

No Scripps code is included in this repository: ADCP and ADFRsuite are downloaded from the Scripps Research Institute by each user. Some of their components are free for academic use only; THIRD_PARTY_NOTICES.md says which, and which of them this package installs.
