# Installing GUI to ADCP Docking

This guide covers installation, dependencies, and access to the browser interface for GUI to ADCP Docking.

Choose the environment where the calculations will run:

| where | what you do | section |
|---|---|---|
| **a Slurm cluster** | install once on the cluster, open the page from your laptop | [1](#1-on-a-slurm-cluster) |
| **your own computer** (Linux, macOS, or Windows through WSL) | install in a Python environment | [2](#2-on-your-own-computer) |
| **Google Colab** | install and run in the notebook runtime | [3](#3-in-google-colab) |

The desktop launcher ([4](#4-the-desktop-launcher)) is an optional convenience
for the first of those: it signs in, starts the tool on the cluster and opens
the page for you.

In Slurm mode, preparation, docking, and rescoring run as batch jobs. The interface server reads files and submits jobs. Follow your cluster policy when choosing where to run the server; section 1.6 describes an interactive-job option.

---

## 1. On a Slurm cluster

Cluster execution targets **Slurm on x86-64 Linux**. Account, partition, module, and network settings depend on the site. PBS, LSF, and SGE are not supported.

### 1.1 What it needs

| tool | from | needed for |
|---|---|---|
| `adcp`, `agfr`, `prepare_receptor` | ADFRsuite 1.0 (Scripps) | docking, grid maps, PDBQT (**required**) |
| `pdb2pqr` | pip / conda | protonation at a chosen pH (**required**) |
| `tleap`, `sander`, `MMPBSA.py` | AmberTools (conda) | MM-GBSA only (optional) |
| `sbatch`, `squeue`, `sacct` | Slurm | already on the cluster |

Python 3.9 or newer with **numpy** and **RDKit**, plus **pdbfixer**,
**OpenMM** and **Biopython**, which rebuild missing side chains, run the
restrained minimisation and read mmCIF.

### 1.2 Install

Run these commands in a location permitted by your cluster's installation policy. Module names vary by site.

```bash
module load Anaconda3                  # or Miniforge, or your own conda
conda create -n adcp -c conda-forge python=3.11 rdkit pdbfixer openmm \
      biopython ambertools -y          # drop ambertools if you will not rescore
conda activate adcp

git clone https://github.com/aadhilhaq/ADCP_docking.git
cd ADCP_docking
pip install -e ".[local]"              # the package, and pdb2pqr

python -m adcp_dock.adfr --dest ~/tools/ADFRsuite
```

The last command downloads ADFRsuite from Scripps, checks the archive against a pinned SHA-256 checksum, and installs it. The installer includes a certificate-error fallback; the checksum check remains required. A checksum verifies the expected archive contents, while TLS verification addresses the server connection.

**For a whole group**, install once in a shared directory and have everyone
use it: `--dest /path/to/group/tools/ADFRsuite`, with the repository cloned
next to it. Each person runs their own server from the shared install.

### 1.3 Point it at the tools

```bash
export ADFRSUITE_HOME=~/tools/ADFRsuite
export AMBERHOME=$CONDA_PREFIX          # only if you installed ambertools
adcp-dock check
```

`adcp-dock check` lists detected tools and missing dependencies without submitting jobs. Run it before starting a calculation. Put the two `export` lines in `~/.bashrc`, or in an `activate.sh` beside the checkout, which the desktop launcher sources before it starts the server.

### 1.4 Tell it about your cluster's rules

Set the required scheduler options below the run name and output folder. Leave them blank only if your cluster accepts its default settings.

* **Account**: many clusters refuse a job that does not name an allocation.
  `sacctmgr show assoc user=$USER format=account` lists yours.
* **Other Slurm options**: anything else your site requires, written as
  `--name=value`, for example `--qos=normal`, `--partition=short` or
  `--constraint=icelake`.

* **Acknowledgement**: if your computing facility asks to be acknowledged
  in publications, set `ADCP_ACKNOWLEDGEMENT` to its text, for example in
  `activate.sh`. Every run records it in its citation files. The launcher's
  cluster presets set it for their own facility.

The account and the other options are added to every job a run submits. From the command line, use
`--account` and `--sbatch-extra` on `adcp-dock new`. Write `--sbatch-extra`
with an `=`, as in `--sbatch-extra="--qos=normal --partition=short"`, so the
argument parser reads the quoted text as its value. Options the tool sets
itself (name, log file, cores, memory, time) cannot be overridden there, and
nothing that could end the line or run a command is accepted.

If the scheduler still refuses a job, the page explains the reason in plain
words, for example an invalid account, a walltime longer than the partition
allows, or a partition that is down.

### 1.5 Start the page

```bash
adcp-dock serve --port 0
```

`--port 0` lets the operating system pick a free port, which is what stops
two people on the same login node colliding. It prints the address, with a
one-time token:

```
GUI to ADCP Docking is at

    http://127.0.0.1:41823/?t=Xq4v...
```

The server binds to 127.0.0.1 and requires an access token. Keep the token private. To connect from your computer, forward the displayed port and open the token-bearing address in your browser:

```bash
ssh -N -L 41823:127.0.0.1:41823 you@login.your-cluster.edu
```

Or let the launcher do all of that (section 4): choose **Other (SLURM)** and
give it your cluster's login address.

### 1.6 If your cluster does not allow servers on login nodes

Some sites stop long-running processes on login nodes, however small. Run
the server inside a small interactive job instead. It still submits the
docking as separate jobs:

```bash
srun --pty -c 1 --mem=2G -t 12:00:00 bash      # add --account etc. if needed
conda activate adcp
hostname                                       # the node, for the tunnel
adcp-dock serve --port 0
```

Then tunnel through the login node to that compute node:

```bash
ssh -N -J you@login.your-cluster.edu -L 41823:127.0.0.1:41823 you@<node>
```

`-J` works where the cluster lets you SSH to a node you have a job on, when permitted by your site. If yours has **Open OnDemand**, an interactive desktop
session is the simplest route of all: run `adcp-dock serve` there and open
the address in that desktop's own browser.

A server on a compute node usually has **no internet**, so fetching a PDB id
from the RCSB fails there. Upload the structure instead, or browse to a file
already on the cluster.

### 1.7 Without the page at all

The main workflow is also available from the command line:

```bash
adcp-dock new $SCRATCH/adcp_runs/trial \
    -r 1YCR --chain A \
    -s A:54,A:62,A:67,A:96,A:99 \
    -f peptides.txt \
    --poses 10 --mmgbsa --cpus 48 --time 08:00:00
adcp-dock run    $SCRATCH/adcp_runs/trial
adcp-dock status $SCRATCH/adcp_runs/trial
adcp-dock report $SCRATCH/adcp_runs/trial
```

`-r` takes a path or a PDB id; `-s` names the residues the box covers; `-f`
reads a peptide file (`name: SEQUENCE` per line, or FASTA), or use `-p
SEQUENCE`, repeated. `run` submits three jobs that wait on each other:
preparation, one array task per peptide, then the report. Run it again after
a failure and only the peptides with no result are docked again.

### 1.8 Things that differ between clusters

| Setting | Possible problem | Action |
|---|---|---|
| **home quota** | runs are large | write runs under `$SCRATCH`; the default already does when `$SCRATCH` is set |
| **array limits** | a long peptide list is refused | check the maximum array size and concurrent-task limit; `--throttle` controls concurrency, not total array size |
| **walltime limits** | a job waits indefinitely with `PartitionTimeLimit` | the page shows it; lower the walltime or choose a partition that allows it |
| **cores per node** | a task asks for more than a node has | set *Cores per task* to your node size or less |
| **module names** | `module load Anaconda3` fails | `module avail` for your site's conda |
| **compute nodes without internet** | PDB ids fail from a job | fetch from the login node, or upload the file |

---

## 2. On your own computer

For a workstation or a laptop with no cluster: **Linux**, **macOS**, or
**Windows through WSL**. The work runs on that machine, one step after
another.

```bash
python3 -m venv ~/.adcp-dock/venv
~/.adcp-dock/venv/bin/pip install "adcp-dock[local] @ git+https://github.com/aadhilhaq/ADCP_docking"
~/.adcp-dock/venv/bin/adcp-dock local
```

Use a virtual environment to keep package dependencies separate from the system Python. On Ubuntu, `sudo apt install python3-venv` first if `venv` is missing.
The desktop launcher's **This computer** option looks for it in exactly this
place, so after these three lines you can start it from the launcher too.

The first `adcp-dock local` downloads and installs ADFRsuite from Scripps
(about 100 MB, checked against its pinned checksum) into
`~/.adcp-dock/ADFRsuite`, then opens the page in your browser. Subsequent launches reuse the installation. Runs are written to `~/adcp_runs`.

* **macOS**: ADFRsuite is built for Intel. On Apple silicon it runs through
  Rosetta 2; if you have never needed Rosetta, run
  `softwareupdate --install-rosetta` once first.
* **Windows**: there is no native build of what this needs, so use WSL.
  Install Ubuntu from the Microsoft Store, open it, and run the three lines
  above inside it. WSL shares `localhost` with Windows, so open the address
  it prints in your usual Windows browser.
* **Linux**: x86-64, which is what ADFRsuite is built for.

**Review the search settings.** Local mode displays the available cores and starts with reduced search settings. These settings are useful for checking the workflow, but do not establish convergence or suitability for a research result. Set the sampling budget for your system and available resources.

**MM-GBSA needs AmberTools**, which is not on pip:
`conda install -c conda-forge ambertools`. Make sure the AmberTools programs are available to the environment that runs `adcp-dock local`, and confirm their paths with `adcp-dock check`. Docking works without it.

`adcp-dock local` refuses to run on a machine that has `sbatch`, because that
is almost always a cluster login node. Use `adcp-dock serve` there; it sends
the work to Slurm. `--force` overrides the check on a workstation that
merely has the Slurm client installed.

---

## 3. In Google Colab

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/aadhilhaq/ADCP_docking/blob/main/notebooks/ADCP_docking_colab.ipynb)

The notebook runs the same interface on a Google Colab machine. Nothing is installed on your own computer. Press ▶ on each step in order.

1. **Install.** Installs the package, PDB2PQR, OpenMM and PDBFixer, then downloads ADFRsuite from the Scripps Research Institute and checks it against its pinned checksum. With **add_MM_GBSA** ticked (the default), it also installs AmberTools from conda-forge into its own folder, using micromamba, for MM-GBSA rescoring. It takes about five minutes, a few more with AmberTools, and ends with a single line such as *Installed in 6 min 10 s, with MM-GBSA*.
2. **Open GUI to ADCP Docking.** Starts the interface on the Colab machine and shows it inside the notebook, across the full width of the page. **Open in a separate tab** shows the same interface in its own browser tab. The list of runs includes a finished example, MDM2 with three peptides docked and rescored with MM-GBSA, to open and explore before docking anything.
3. **Download a run.** Optional. Finished runs can be downloaded from the interface; this step zips a run and downloads it instead.

Things to know before you start:

* The calculations run on the Colab machine, one step after another. The number of cores, the memory, and the session length depend on your Colab plan. The interface shows the cores it found and starts with small search settings.
* Colab deletes the session's files when the session ends. Download the runs you want to keep. **Runtime → Disconnect and delete runtime** removes every run at once; the example returns the next time step 2 runs.
* Each person who opens the notebook gets a separate Colab machine, so nobody else sees your runs.
* MM-GBSA rescoring needs AmberTools, which step 1 installs when **add_MM_GBSA** is ticked. Without it, docking and the ADCP scores work normally and the interface disables the MM-GBSA option. On two cores MM-GBSA is slow: allow several minutes per pose.
* ADFRsuite includes components that are free for academic use only (see the end of this guide).

Use the notebook to explore the workflow and inspect inputs. Reduced search settings are not a validated production protocol.

---

## 3a. ADCP 1.1 (optional)

ADCP 1.1 is the current version of AutoDock CrankPep. It adds non-standard and D amino acids, all-atom poses in one file, and an optional OpenMM minimization and re-ranking of the top poses. The Scripps Research Institute distributes it as its own Python 3.7 environment. Install it, on a cluster login node, a workstation or in Colab, with:

```bash
adcp-dock install-adcp11                      # into ~/.adcp-dock/adcpsuite
adcp-dock install-adcp11 --dest /path/adcpsuite
```

This follows Scripps' own installation script: it creates the environment with micromamba (downloaded if needed), installs OpenMM, Open Babel and PDBFixer from conda-forge, and fetches the Scripps packages from ccsb.scripps.edu, checking each against a pinned SHA-256. It takes about ten minutes and 3 GB. Scripps' own instructions are at <https://ccsb.scripps.edu/adcpv11/downloads/>.

The interface finds ADCP 1.1 at `$ADCPSUITE_HOME`, in a folder named `adcpsuite` beside ADFRsuite or in a tools folder given with `--tools-root`, at `~/.adcp-dock/adcpsuite`, or at `/content/adcpsuite` in Colab. ADFRsuite 1.0 is still needed, for receptor preparation and the grid maps.

Each run chooses its engine under **Docking engine**: *ADCP, newest installed* (the default) uses 1.1 when it is found and 1.0 otherwise. On the command line, `adcp-dock new ... --engine 1.1`. With ADCP 1.1:

* **Non-standard and D amino acids.** Write a non-standard residue as the one-letter code of the standard residue whose coarse model it uses, followed by its library name in angle brackets, and a D residue with `&` before its letter: `N<BHD>PIS&DVD` is seven residues. The names are those of ADCP's side-chain libraries (the SwissSideChain set). MM-GBSA rescoring is skipped for these peptides, since its force field has no parameters for them; use the OpenMM minimization instead.
* **OpenMM minimization.** *Poses to minimize with OpenMM* (`--omm-minimize N`) minimizes the top N poses in implicit solvent and reports the lowest E(complex) minus E(receptor) beside the scores. The minimized complexes and their Amber topologies stay in the run's docking folder.

A run can also start from a target file made in AGFR or AGFRGUI: give its path under *or a target file* (`--target-file`). The receptor, maps and box then come from that file, and this package's receptor preparation is not run. The file's name must not be changed, since ADCP finds the maps by it.

Scripps' ADCP 1.1 environment includes the MSMS library (mslib) and the UT packages, which are free for academic use only; check their terms before commercial use.

---

## 4. The desktop launcher

The launcher connects to a cluster over SSH, supports two-factor
authentication, starts `adcp-dock serve` on the login node, forwards a port
and opens your browser. Closing it stops the server.

### 4.1 Get it

* **Windows and macOS**: download from the repository's **Releases** page.
  Windows is a folder: keep it together and run the `.exe` inside. macOS is
  an Apple silicon `.app`.
* **Any machine with Python**: `pip install "adcp-dock[launcher] @
  git+https://github.com/aadhilhaq/ADCP_docking"`, then `adcp-dock-launcher`.

The bundles are not yet code-signed, so the first run shows a warning:
Windows says *Windows protected your PC* (**More info → Run anyway**); macOS
refuses to open it (**right-click → Open → Open**). Signing removes those
prompts.

### 4.2 What it asks for

| field | what to put |
|---|---|
| **Cluster** | a preset, or **Other (SLURM)** for any Slurm cluster |
| **Address** | your cluster's login node |
| **Username** | your cluster username |
| **Tool directory** | where the repository is cloned on the cluster |
| **ADFRsuite / Amber roots** | blank if `ADFRSUITE_HOME` is set there |

Everything but the password is remembered for next time. It signs in with a
password and a second factor (keyboard-interactive) or with your SSH key,
whichever your cluster uses.

---

## 5. If something does not work

| what you see | what it means |
|---|---|
| `N tool(s) missing` from `check` | ADFRsuite or Amber is not reachable: set `ADFRSUITE_HOME` / `AMBERHOME`, or pass `--tools-root` |
| `sbatch was not found` | Slurm is not available in this environment; load the required cluster environment or use `adcp-dock local` on a workstation |
| the downloaded ADFRsuite *does not match its pinned checksum* | the archive differs from the pinned version; inspect the installer message before retrying |
| `the account on the job is not one you can use` | this cluster needs an **Account** (1.4) |
| `the grid target is unusable` | the `.trg` was truncated or built from another structure. It is rebuilt once automatically; twice means agfr is failing, and its output is in the prepare log |
| `tleap did not build topologies` | the preparation log names the residue. It is usually a HETATM kept by *keep heteroatoms*, or a gap that *rebuild missing loops* would close |
| `the docking box also reaches chain(s) ...` | the box covers a chain that is not loaded. If it is the far wall of the groove, include it; if it is the partner your peptide should displace, leave it out |
| `pose N discarded: the bonded terms did not cancel` | the energy decomposition failed a consistency check; inspect the chain selections, topology, and rescoring log |
| the queue rejects the array | check the scheduler message for total array-size and concurrency limits; reduce the peptide list or adjust *Max tasks at once* as appropriate |

---

## 6. Uninstall

```bash
pip uninstall adcp-dock
rm -rf ADCP_docking ~/.adcp-dock
```

The commands above remove the named installation directories. Run outputs stored elsewhere remain available; move any outputs saved inside those directories before removing them.

---

Use `adcp-dock cite <run>` to list the tools recorded for a particular run, with
DOIs, and every finished run carries `HOW_TO_CITE.txt`, `CITATION.bib` and
`CITATION.ris`. ADFRsuite bundles components that are free for academic use
only (MSMS, the UT packages, STRIDE); commercial users should contact the
Scripps Research Institute.
