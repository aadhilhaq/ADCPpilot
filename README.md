# GUI to ADCP Docking

GUI to ADCP Docking provides a browser interface and command-line workflow for peptide docking with [AutoDock CrankPep (ADCP)](https://ccsb.scripps.edu/adcp/). It connects receptor preparation, grid generation, docking, optional MM-GBSA rescoring, and reporting. Calculations can run on a Slurm cluster, a local computer, or Google Colab.

Provide a receptor structure and peptide sequences to generate predicted poses and a ranked results table. ADCP performs the docking calculations. This package manages the surrounding workflow.

Developed by Aadhil Haq, HOPE Lab, under the supervision of Dr. Sandun Fernando. This project is independent of the Scripps Research Institute, which developed ADCP and ADFRsuite.

## Workflow

```text
Receptor structure and peptide sequences
    → Receptor preparation
    → Receptor conversion and grid generation
    → ADCP docking for each peptide
    → Optional MM-GBSA rescoring
    → Structures, ranked results, and citation files
```

On Slurm, a preparation job is followed by an array with one task per peptide. Each task includes optional rescoring. The reporting job runs after the array terminates, so it can collect partial results when a task fails.

## Getting started

Choose where the calculations will run. [INSTALL.md](INSTALL.md) provides the full setup instructions, dependency requirements, and troubleshooting steps.

| Environment | Setup |
|---|---|
| Slurm cluster | Install the package and scientific tools on the cluster. Connect through SSH or the desktop launcher. |
| Local computer | Use Linux, macOS, or Windows through WSL, with the required scientific tools. |
| Google Colab | Open the [notebook](https://colab.research.google.com/github/aadhilhaq/ADCP_docking/blob/main/notebooks/ADCP_docking_colab.ipynb). Nothing is installed on your computer. |

### Slurm cluster

After setting up the environment described in [INSTALL.md](INSTALL.md#1-on-a-slurm-cluster):

```bash
git clone https://github.com/aadhilhaq/ADCP_docking.git
cd ADCP_docking
pip install -e ".[local]"
python -m adcp_dock.adfr --dest ~/tools/ADFRsuite
export ADFRSUITE_HOME=~/tools/ADFRsuite
adcp-dock check
adcp-dock serve --port 0
```

The server prints a local address and access token. Use SSH port forwarding to open that address from your computer, or use the optional [desktop launcher](https://github.com/aadhilhaq/ADCP_docking/releases). See [connection instructions](INSTALL.md#15-start-the-page).

Computational stages run as batch jobs. Follow your cluster's policy for hosting the interface server; an interactive-job option is described in the installation guide.

### Local computer

```bash
python3 -m venv ~/.adcp-dock/venv
~/.adcp-dock/venv/bin/pip install "adcp-dock[local] @ git+https://github.com/aadhilhaq/ADCP_docking"
~/.adcp-dock/venv/bin/adcp-dock local
```

The first launch installs ADFRsuite if needed. MM-GBSA requires a separate AmberTools installation. Platform requirements are listed in [INSTALL.md](INSTALL.md#2-on-your-own-computer).

### Google Colab

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aadhilhaq/ADCP_docking/blob/main/notebooks/ADCP_docking_colab.ipynb)

The notebook runs the same interface on a Colab machine. Step 1 installs the package and downloads ADFRsuite, which takes about five minutes. Step 2 opens the interface inside the notebook, or in a separate browser tab. The calculations run on the Colab machine, so the available cores, memory, and session length depend on your Colab plan. Colab deletes the session's files when it ends; download finished runs from the interface, or with step 3. MM-GBSA rescoring is not available in Colab, because AmberTools is not installed there. See [INSTALL.md](INSTALL.md#3-in-google-colab).

Local and Colab runs are useful for learning the workflow and checking inputs. Reduced search settings do not establish sampling convergence. Choose production settings for the target, peptide, and available resources.

### Command line

Example commands, after configuring the required tools:

```bash
adcp-dock new runs/trial -r 1YCR --chain A \
    -s A:54,A:62,A:67,A:96,A:99 -p "p53: ETFSDLWKLLP" --poses 10 --mmgbsa
adcp-dock run runs/trial
adcp-dock status runs/trial
adcp-dock report runs/trial
```

`-r` accepts a structure path or PDB ID. `-s` specifies residues used to define the docking box. Use `-p` for individual peptides or `-f` for a peptide file. Remove `--mmgbsa` if rescoring is not required. Review the search and resource settings before submitting a run.

Running `adcp-dock run` again submits peptides without completed results.

## Requirements

Python 3.9 or newer, NumPy, and RDKit are required.

| Tool | Source | Purpose |
|---|---|---|
| `adcp`, `agfr`, `prepare_receptor` | ADFRsuite | Docking, grid generation, and receptor conversion |
| `pdb2pqr` | PDB2PQR | Protonation at the selected pH |
| `tleap`, `sander`, `MMPBSA.py` | AmberTools | Optional MM-GBSA rescoring |
| `sbatch` | Slurm | Cluster execution |

PDBFixer, OpenMM, and Biopython support structural repair, restrained minimization, and mmCIF input. The preparation log reports unavailable steps.

Run `adcp-dock check` to inspect the installation. The ADFRsuite installer checks the downloaded archive against a pinned checksum. Fetching a receptor by PDB ID requires internet access; downloaded entries are cached.

## Peptide input

Use one-letter amino acid codes for sequences of 5 to 30 residues. The parser rejects B, J, O, U, X, and Z.

Supported formats are a sequence, a named sequence, or FASTA:

```text
ETFSDLWKLLP
```

```text
p53: ETFSDLWKLLP
```

```fasta
>p53
ETFSDLWKLLP
```

Lines beginning with `#` are comments. Use a colon to separate a name from its sequence. Spaces within sequence lines are removed.

## Results and interpretation

The report writes `results.tsv`, `results.json`, and `results.txt`. The default ranking uses the mean of up to five top poses. If MM-GBSA summaries are available, the report uses that ranking criterion; otherwise it uses ADCP scores. Peptides without a value for the selected criterion remain in the report without a rank.

The report also includes individual score summaries, pose counts, and available measures of score variation. Variation across poses does not establish statistical significance, uncertainty in binding affinity, or convergence of the search. Inspect the structures and sampling settings alongside the scores.

Use MM-GBSA as an approximate ranking aid within a consistently prepared run. These values are not experimentally measured binding free energies or dissociation constants.

## Workflow checks

The workflow checks grid-target contents before docking, identifies receptor and peptide selections by chain for rescoring, and checks bonded-energy cancellation in single-trajectory MM-GBSA calculations. Poses that fail the rescoring consistency check are excluded from those results.

Preparation decisions and incomplete steps are recorded in the logs, and the results page reports when the receptor preparation was incomplete. Review these records before interpreting a ranking. Passing workflow checks does not validate the predicted binding pose or affinity.

## Tests and contributions

Run the test suite on an allocated compute node:

```bash
srun --pty -c 8 -t 30 bash tests/run_all.sh
```

Add the account and partition required by your cluster. Test scripts report pass, skip, or failure; skipped checks have not been verified. See [CONTRIBUTING.md](CONTRIBUTING.md) for individual tests, contribution guidance, and the shared engine modules.

Report problems through [GitHub Issues](https://github.com/aadhilhaq/ADCP_docking/issues) or **Send feedback** in the interface. Include the software revision, settings, relevant log output, and steps to reproduce the problem.

## Citation

Cite ADCP when reporting docking calculations:

> Zhang Y, Sanner MF. AutoDock CrankPep: combining folding and docking to predict protein-peptide complexes. *Bioinformatics*. 2019;35(24):5121-5127. [DOI](https://doi.org/10.1093/bioinformatics/btz459).

Cite ADFRsuite and the preparation, rescoring, and visualization tools used in the analysis. The package can generate references for a specific run:

```bash
adcp-dock cite                 # all supported tools
adcp-dock cite <run>           # tools recorded for one run
adcp-dock cite --bibtex        # use --ris for EndNote or Zotero
```

Run outputs include `HOW_TO_CITE.txt`, `CITATION.bib`, and `CITATION.ris`. For this interface, use GitHub's **Cite this repository** button or [CITATION.cff](CITATION.cff), and record the software version used.

## License

This package is distributed under the [MIT license](LICENSE). External scientific programs retain their own licenses and terms of use. Check the terms for the installed versions of ADCP, ADFRsuite, AmberTools, PDB2PQR, and other dependencies.
