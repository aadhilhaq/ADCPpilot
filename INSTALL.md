# Installing GUI to ADCP Docking

This guide is written for someone who has used a terminal only a few times.
Every command is in a box: copy it, paste it into the terminal, press Enter,
and wait for it to finish before the next one.

Choose where the calculations will run:

| where | what you do | section |
|---|---|---|
| **a Slurm cluster** | install once on the cluster, open the page from your laptop | [1](#1-on-a-slurm-cluster) |
| **your own computer** (Linux, macOS, or Windows through WSL) | install in a Python environment | [2](#2-on-your-own-computer) |
| **Google Colab** | install and run in the notebook runtime | [3](#3-in-google-colab) |

The desktop launcher ([5](#5-the-desktop-launcher)) is an optional convenience
for the first of those: it signs in, starts the tool on the cluster and opens
the page for you.

---

## 1. On a Slurm cluster

Cluster execution targets **Slurm on x86-64 Linux**. PBS, LSF and SGE are not
supported.

### 1.1 Before you start

**Login nodes and compute nodes.** When you connect to a cluster with `ssh`,
you arrive on a *login node*: a computer shared by everyone who is logged in,
meant for small tasks such as copying files, downloading and submitting work.
The calculations run on *compute nodes*, which you reach only through the
scheduler, Slurm. Compute nodes usually cannot download anything, so the
install runs on a login node. Docking never does: adcp-dock sends every
calculation to Slurm for you.

**What you need:**

* an account on the cluster, and its login address (for example
  `login.your-cluster.edu`; your cluster's documentation gives it);
* a terminal on your own computer: **Terminal** on macOS, **PowerShell** on
  Windows 10 or 11, or any terminal on Linux. All of them have `ssh`;
* a folder on the cluster with room for about **150,000 files and 10 GB**.
  Clusters limit how many files you may keep as well as how much space they
  take, and this install is many small files: about 150,000 while it is
  being installed, and about 93,000 once it has deleted what it downloaded.
  Your home folder is usually too small: most clusters allow only a few
  thousand files there. Use your scratch or project folder. On many clusters
  `echo $SCRATCH` prints yours, and the cluster's quota command (often
  `showquota`, or `quota -s`) shows how much room it has.

Nothing else needs to be installed first, and nothing needs administrator
rights.

### 1.2 Install

**1. Connect to the cluster** from the terminal on your own computer. Use your
own username and your cluster's login address:

```bash
ssh you@login.your-cluster.edu
```

It asks for your password, and on many clusters for a second factor. Nothing
appears on the screen while you type the password; that is normal.

**2. Go to a folder with room** and download the code (under a minute):

```bash
cd $SCRATCH
git clone https://github.com/aadhilhaq/ADCP_docking.git
cd ADCP_docking
```

If `cd $SCRATCH` says *No such file or directory*, your cluster names the
folder differently; use the path its documentation gives, here and below.

**3. See what it would do** (a few seconds; it changes nothing):

```bash
./install.sh --check $SCRATCH/adcp
```

It lists what it found on the cluster, what it would add, and whether the
folder has room. If it says **NOT ENOUGH ROOM**, choose another folder before
going on.

**4. Install** (about 7 to 10 minutes when it builds everything, under a minute
when the cluster already has it all):

```bash
./install.sh $SCRATCH/adcp
```

`$SCRATCH/adcp` is where the installed software goes; it is created if it is
not there. The script prints what it is doing as it goes, never asks a
question, and writes nothing into your home folder. In order, it:

1. looks for what the cluster already has (below);
2. builds a Python environment with RDKit, AmberTools, PDB2PQR, OpenMM,
   PDBFixer and Biopython, from conda-forge;
3. downloads ADFRsuite from the Scripps Research Institute, checks it against
   a pinned checksum, and installs it;
4. installs ADCP 1.1 beside it with Scripps' own recipe
   ([section 4](#4-adcp-11-optional)); `--no-adcp11` leaves it out;
5. fetches ADapt beside it, once ADapt's repository is public
   ([ADapt](#adapt-optional));
6. deletes the packages it downloaded, and writes `activate.sh` in the
   checkout, which records where everything is, then checks the result.

Measured on a clean account on the Grace cluster at Texas A&M, the whole of
it took 7 and 10 minutes on two days, most of it the two environments, with
about 1 GB of downloads; a slower network or a busy login node takes longer.
It leaves about 95,000 files. It uses your
cluster's own conda when it is recent, and otherwise fetches micromamba, a
single file, into the folder; no conda module is needed.

**What success looks like.** The output ends with a summary and the next
commands, for example (your paths and times will differ):

```
== Summary ==
  Python environment  /scratch/you/adcp/envs/adcp
  AmberTools          in that environment
  ADFRsuite           /scratch/you/adcp/tools/docking/ADFRsuite
  ADCP 1.1            /scratch/you/adcp/tools/docking/adcpsuite
  ADapt               not installed (optional; ./update.sh fetches it once it is public)

  Finished in 6 min 58 s.
```

If something could not be installed, the summary also lists what is turned
off because of it, for example *MM-GBSA rescoring (needs AmberTools)*.

**5. Check it:**

```bash
source activate.sh
adcp-dock check
```

Every program should show a path, the Python packages a version, and the last
line should read `sbatch: found`. Docking needs `adcp`, `agfr`,
`prepare_receptor` and `pdb2pqr`; MM-GBSA rescoring needs the four Amber
programs.

#### Messages you may see, which are normal

* `Security Warning: This transaction includes executing package scripts`:
  micromamba says this whenever a package runs a small setup script of its
  own; the packages come from conda-forge, and nothing needs doing.
* `Executing post-link script for package 'ocl-icd-system'`: one of those
  scripts, which lets OpenMM find the system's OpenCL drivers.
* `ADapt: not installed: ... is not public yet`: ADCP docking works without
  it, and `./update.sh` fetches it once it is public.

#### If it stops

When something goes wrong, the script stops with a line such as
`install.sh STOPPED while building the environment`, followed by what usually
causes it. Nothing is lost: **run the same command again** and it carries on
from what is already there. The same is true if your connection drops while
it runs: log in again, `cd` to the checkout and run the same line.

| what you see | what to do |
|---|---|
| `NOT ENOUGH ROOM` | choose a folder with more room and give that instead, or free some files; nothing was changed |
| `Could not resolve host`, `Connection timed out`, or a download that failed | you are probably on a compute node; log out, log in to a login node and run it again |
| `Disk quota exceeded` or `No space left on device` | the folder ran out of room part way; free some, or start again with a folder that has more |
| `No Python environment ... and no folder was given` | add the folder: `./install.sh $SCRATCH/adcp` |
| `could not be added without changing what is already in ...` | an environment you pointed at (below) could not take a package without changing others, so it was left as it was; the summary says what that turns off |
| `ADapt: not installed: ... is not public yet` | nothing to do; ADCP docking works without it, and `./update.sh` fetches it later |
| `ADFRsuite did not install properly` | follow the two commands it prints: read its log, delete the half installed folder, and run the install again |
| `Run this with bash` or `do not source it` | type it exactly as shown: `./install.sh $SCRATCH/adcp` |

If it stops in the same place twice, open an issue at
<https://github.com/aadhilhaq/ADCP_docking/issues> with the last 30 lines of
the output.

### 1.3 Every time you use it

Each new terminal starts without the tools on its path. Load them with:

```bash
source /path/to/ADCP_docking/activate.sh
```

using the folder you cloned into (`pwd` inside it prints the full path).
`activate.sh` sets `ADFRSUITE_HOME`, `AMBERHOME`, `ADCP_PYTHON` and `PATH`.
The desktop launcher sources it before it starts the server, and every job a
run submits names the same programs itself, so a cluster that starts its jobs
with a clean environment still finds them.

`adcp-dock check` lists what it found without submitting anything. Run it
whenever something seems missing.

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

On the login node, after `source activate.sh`:

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

The server binds to 127.0.0.1 and requires that token, so only you can open
it. Keep the token private. To open it from your own computer, open a
**second** terminal there and forward the port it printed (41823 here):

```bash
ssh -N -L 41823:127.0.0.1:41823 you@login.your-cluster.edu
```

This one prints nothing and keeps running; that is the tunnel. Then paste the
whole address, token included, into your browser. Many clusters have several
login nodes behind one address, and the tunnel must reach the one running the
server. If the page does not open, your cluster's documentation gives the
address of each login node; use the same one in both terminals.

To stop, press Ctrl-C in both terminals. Runs already submitted keep going in
Slurm, and the page shows them again the next time it starts.

Or let the launcher do all of that (section 5): choose **Other (SLURM)** and
give it your cluster's login address.

### 1.6 If your cluster does not allow servers on login nodes

Some sites stop long-running processes on login nodes, however small. Run
the server inside a small interactive job instead. It still submits the
docking as separate jobs:

```bash
srun --pty -c 1 --mem=2G -t 12:00:00 bash      # add --account etc. if needed
source /path/to/ADCP_docking/activate.sh
hostname                                       # the node, for the tunnel
adcp-dock serve --port 0
```

Then tunnel through the login node to that compute node:

```bash
ssh -N -J you@login.your-cluster.edu -L 41823:127.0.0.1:41823 you@<node>
```

`-J` works where the cluster lets you SSH to a node you have a job on. If
yours has **Open OnDemand**, an interactive desktop session is the simplest
route of all: run `adcp-dock serve` there and open the address in that
desktop's own browser.

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
reads a peptide file (`name: SEQUENCE` per line, FASTA, or columns of name,
sequence and an optional ring), or use `-p SEQUENCE`, repeated. `run` submits three jobs that wait on each other:
preparation, one array task per peptide, then the report. Run it again after
a failure and only the peptides with no result are docked again.

### 1.8 Things that differ between clusters

| Setting | Possible problem | Action |
|---|---|---|
| **file allowance** | the install or a run stops with *Disk quota exceeded* | install and write runs under scratch or project space; `install.sh` checks for room before it starts |
| **array limits** | a long peptide list is refused | check the maximum array size and concurrent-task limit; `--throttle` controls concurrency, not total array size |
| **walltime limits** | a job waits indefinitely with `PartitionTimeLimit` | the page shows it; lower the walltime or choose a partition that allows it |
| **cores per node** | a task asks for more than a node has | *Cores per task* starts at what a node of the default partition has, up to 48; a partition with smaller nodes needs it lowered |
| **conda** | there is no conda module, or an old one | nothing: `install.sh` fetches micromamba when it finds no recent conda |
| **AmberTools module** | the cluster provides AmberTools | load it before installing (1.10); `install.sh` then uses it |
| **job environment** | the cluster starts jobs with a clean environment | nothing: each job script names ADFRsuite and AmberTools itself |
| **compute nodes without internet** | PDB ids fail from a job | fetch from the login node, or upload the file |

### 1.9 Update later

On a login node, in the folder you cloned into (under a minute):

```bash
cd /path/to/ADCP_docking
git pull; ./update.sh
```

`update.sh` moves the checkout to the newest version, brings ADapt up to date
beside ADFRsuite ([ADapt](#adapt-optional)), and ends with the same list as
`adcp-dock check`. It never builds an environment; when ADCP 1.1 is missing it
gives the command that adds it. A server that was already running keeps the
old code until you stop it with Ctrl-C and start it again.

### 1.10 Using what your cluster already has

`install.sh` looks for what it needs before it builds anything, and adds only
what is missing.

* **AmberTools as a module.** Many clusters provide AmberTools as a module.
  Load it before running `install.sh`, and the environment is built without
  its own copy:

  ```bash
  module spider AmberTools          # lists the versions, and what to load first
  module load AmberTools/22.3       # the name your cluster gives
  ./install.sh $SCRATCH/adcp
  ```

  It is used once its four programs (`tleap`, `sander`, `MMPBSA.py`,
  `ante-MMPBSA.py`) are there and `MMPBSA.py` runs. `install.sh` then writes
  small scripts in `tools/amber/bin` that carry the module's settings, so you
  do not need to load the module again later, and the module's Python
  settings never reach adcp-dock's own Python. If the module is updated or
  removed, run `./install.sh` again. The same goes for an AmberTools that is
  on your `PATH` or named by `AMBERHOME`.

* **A Python environment you already have.** Point at it with `ADCP_ENV`:

  ```bash
  ADCP_ENV=/path/to/env ./install.sh
  ```

  It needs Python 3.9 or newer. Whatever it lacks of numpy, RDKit, PDB2PQR,
  OpenMM, PDBFixer, Biopython and AmberTools is added to it, with conda for a
  conda environment and pip for any other, and every package already there is
  held at its version, so nothing another user relies on changes. If you
  cannot write to it, it is used as it is, and the script says what is
  missing and what that turns off (for example OpenMM: the receptor's
  restrained minimisation). If it lacks numpy, RDKit or PDB2PQR and cannot be
  given them, it is not used, and `install.sh` builds an environment in the
  folder you give instead. The adcp-dock code is never installed into an
  environment `install.sh` did not build; the `adcp-dock` command runs it
  from the checkout.

* **For a whole group.** Run it once in a shared folder and have everyone use
  that checkout and its `activate.sh`; each person runs their own server.
  `ADCP_SHARED=/path/to/group/folder` makes `install.sh` also look in that
  folder's `envs/adcp` and `tools/docking`.

#### By hand

The same thing step by step, for a site where you would rather build the
environment yourself. Everything stays inside the checkout, including the
conda environment and the download caches; module names vary by site.

```bash
git clone https://github.com/aadhilhaq/ADCP_docking.git
cd ADCP_docking
export CONDA_PKGS_DIRS=$PWD/.conda_pkgs PIP_CACHE_DIR=$PWD/.pip_cache

module load Anaconda3                  # or Miniforge, or your own conda
conda create -p ./env -c conda-forge python=3.11 rdkit pdbfixer openmm \
      biopython ambertools -y          # drop ambertools if you will not rescore
conda activate ./env
pip install -e ".[local]"              # the package, and pdb2pqr

python -m adcp_dock.adfr --dest tools/docking/ADFRsuite
```

The last command downloads ADFRsuite from Scripps, checks the archive against
a pinned SHA-256 checksum, and installs it. Then set the two variables
yourself in each new terminal:

```bash
export ADFRSUITE_HOME=/path/to/ADCP_docking/tools/docking/ADFRsuite
export AMBERHOME=$CONDA_PREFIX          # only if you installed ambertools
```

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
* Some components of ADFRsuite and of ADCP 1.1 are free for academic use only; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Use the notebook to explore the workflow and inspect inputs. Reduced search settings are not a validated production protocol.

---

## 4. ADCP 1.1 (optional)

ADCP 1.1 is the current version of AutoDock CrankPep. It adds non-standard and D amino acids, all-atom poses in one file, and an optional OpenMM minimization and re-ranking of the top poses. The Scripps Research Institute distributes it as its own Python 3.7 environment. On a cluster, `install.sh` installs it beside ADFRsuite unless given `--no-adcp11`. Otherwise, or to add it later, install it on a cluster login node, a workstation or in Colab with:

```bash
adcp-dock install-adcp11                      # beside ADFRsuite, in the tools folder
adcp-dock install-adcp11 --dest /path/adcpsuite
```

This follows Scripps' own installation script: it creates the environment with micromamba (downloaded if needed), installs OpenMM, Open Babel and PDBFixer from conda-forge, and fetches the Scripps packages from ccsb.scripps.edu, checking each against a pinned SHA-256. It takes a few minutes and 3 GB, and deletes the packages it downloaded once it is done. Scripps' own instructions are at <https://ccsb.scripps.edu/adcpv11/downloads/>.

Without `--dest` it goes beside ADFRsuite (`$ADFRSUITE_HOME/../adcpsuite`), in the installation's own tools folder, and on a computer without ADFRsuite set up yet into `~/.adcp-dock/adcpsuite`. The interface finds ADCP 1.1 at `$ADCPSUITE_HOME`, in a folder named `adcpsuite` beside ADFRsuite or in a tools folder given with `--tools-root`, at `~/.adcp-dock/adcpsuite`, or at `/content/adcpsuite` in Colab. ADFRsuite 1.0 is still needed, for `prepare_receptor` and for ADCP 1.0. The suite carries its own AGFR, found from the same installation, and the grid maps are built with it wherever the suite is installed (`--agfr auto`, the default; `--agfr 1.0` asks for ADFRsuite's).

Each run chooses its engine under **Docking engine**: *ADCP, newest installed* (the default) uses 1.1 when it is found and 1.0 otherwise. On the command line, `adcp-dock new ... --engine 1.1`. With ADCP 1.1:

* **Non-standard and D amino acids.** Write a non-standard residue as the one-letter code of the standard residue whose coarse model it uses, followed by its library name in angle brackets, and a D residue with `&` before its letter: `N<BHD>PIS&DVD` is seven residues. The names are those of ADCP 1.1's library, built by the Scripps Research Institute from the SwissSidechain database: 209 non-natural side chains, each in L and D form, among them phosphorylated serine, threonine and tyrosine, oxidised methionine and tryptophan, citrulline and nitrotyrosine, and many synthetic residues. The 197 that ADCP gives a coarse potential can be written in this syntax, and the page offers those. Common modifications such as lysine acetylation or methylation are not in it. ADCP builds and docks every residue in the library, but its accuracy on peptides containing them has not yet been benchmarked in a publication, so treat such results as exploratory. Runs that use them also cite SwissSidechain, and GAFF when they are minimized. MM-GBSA rescoring is skipped for these peptides, since its force field has no parameters for them; use the OpenMM minimization instead.
* **Choosing them in the page.** Open a peptide's &#8942; menu and choose *Non-standard and D amino acids*. Click a residue, then pick its D form or one of the non-standard residues ADCP has for it; the list offers only what the installed ADCP 1.1 contains. Each choice is drawn below the list as the residue sits in the chain, with the charge ADCP models it with, so variants such as SEP and SEP2 (charge &minus;1 and &minus;2) can be told apart. A D residue is drawn as the mirror image of the L one, which is how ADCP builds it. Saving writes the sequence in the syntax above.
* **OpenMM minimization.** *Poses to minimize with OpenMM* (`--omm-minimize N`) minimizes the top N poses in implicit solvent and reports the lowest E(complex) minus E(receptor) beside the scores. ADCP uses the Amber ff99SB force field with the GBn2 solvent model, holds receptor residues away from the peptide in place, and runs at most five minimizer iterations, its default, so this relaxes each docked complex rather than refining it. The minimized complexes and their Amber topologies stay in the run's docking folder. It covers nearly every residue the dialog offers. ADCP 1.1.21 has no parameters for IYR, PBF, PHI and TYI, and its parameters for FVAL, FPG2, FPG3, FTR and FT6 fail in OpenMM; a peptide containing one of these is still docked, its poses are left unminimized, and the run's log says so.

Every run, on either engine, also takes these settings; each is named the same in the page and on the command line:

* **Where the box goes** (`--site-mode`): around residues you pick (`-s`, the default), around a ligand or peptide in the structure (`--site-mode ligand --site-ligand B`, or a hetero group as `A:NUT:201`), or around the pocket AutoSite finds once the receptor is prepared (`--site-mode pocket`, with `--site-pocket best` or `all`).
* **The grid maps**: which AGFR builds them (`--agfr`), which AutoSite places the points the search starts from (`--autosite`, 1.1 by default), and the receptor gradient, off unless asked for (`--receptor-gradient`).
* **Peptides**: linear, or a ring head to tail or through a disulfide (`--cyclic none|backbone|disulfide` for the run; a disulfide needs two cysteines in each peptide docked that way, and one with fewer is refused), and the share of replicas that start as a helix (`--helical 25`, the default; 0 starts every replica as a coil, as runs did before). Each peptide can have its own ring: on its row in the page, or in a third column of the `-f` file (`p53,ETFSDLWKLLPEN,head to tail`, or a header row `name,sequence,ring`), with `linear`, `head to tail` or `disulfide`; a peptide without one takes `--cyclic`.
* **Output**: the ranked poses written per peptide (`--poses`, up to 10 with ADCP 1.0 and 100 with 1.1), and every replica's trajectory, packed into one archive per peptide (`--keep-trajectories`).

A run can also start from a target file made in AGFR or AGFRGUI: give its path under *or a target file* (`--target-file`). The receptor, maps and box then come from that file, and this package's receptor preparation is not run. The file's name must not be changed, since ADCP finds the maps by it.

Scripps' ADCP 1.1 environment includes the MSMS library (mslib) and the UT packages, which are free for academic use only; check their terms before commercial use. A 1.1 run still uses ADFRsuite 1.0's `prepare_receptor`, and the packages it loads there carry non-commercial terms of their own that ADFRsuite's `LICENSE.txt` does not list; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

### ADapt (optional)

ADapt, by Michel Sanner, shows what a PDB file holds before you dock into it: each biological assembly, alternate locations, side chains with missing atoms, missing residues, ligands, cofactors and metals. It prepares the receptor the way AGFR does, with AGFR's own code. It is a separate project, not included in this package, and it runs on the ADCP 1.1 environment above, so it is used only where ADCP 1.1 is installed.

#### How it is fetched

`install.sh` and `update.sh` fetch it from <https://github.com/michelsanner/adapt> into a folder named `adapt` beside ADFRsuite, where the interface looks for it. They fetch it anonymously: git runs so that it never asks for a username or password and never uses a saved credential, so the same thing happens for everyone, on every cluster.

* **While the repository is private**, as it is today, the fetch fails at once and the script says so in one line, then carries on:

  ```
  ADapt: not installed: https://github.com/michelsanner/adapt.git is not public yet. ADCP docking works without it; ./update.sh fetches it once it is public, or fetch it by hand with access (INSTALL.md, ADapt)
  ```

  Where ADapt is already installed, the line is instead `ADapt: 1005e87, kept as it is: https://github.com/michelsanner/adapt.git is not public yet`. If you have access now, clone it yourself beside ADFRsuite, over whichever address your access works with (ssh as here, or https with a token):

  ```bash
  source activate.sh
  git clone git@github.com:michelsanner/adapt.git "$(dirname "$ADFRSUITE_HOME")/adapt"
  ```

  `update.sh` never uses that address afterwards: it fetches from the public one.

* **Once it is public**, the next `./update.sh`, or `./install.sh` on a new installation, clones it. Nothing else needs doing.

`install.sh` fetches ADapt only where there is no `adapt` folder yet, and leaves one that is already there to `update.sh`. After updating the checkout, `update.sh` looks at that folder:

* **Missing**: it tries the anonymous clone above.
* **A git checkout**: it fetches the newest `master` from the https address, never through the remote the checkout has configured (on a cluster whose ssh setup gives github.com another repository's deploy key, that remote fails), and moves the checkout forward only when that is a plain fast forward. A checkout with local changes to its files, on another branch, or with commits of its own is left exactly as it is, and the line says why.
* **Not a git checkout**, for example one copied without its `.git` folder: it is left as it is, and the line says it cannot be updated automatically.

A checkout whose commit is no longer in the published history, as happens if the history is rewritten when the repository is made public, is left as it is too. In that case, and for a folder that is not a git checkout, move the folder aside (`mv adapt adapt.old` beside ADFRsuite) and the next `./update.sh` fetches a fresh checkout in its place.

After moving ADapt forward, `update.sh` checks that the new version still has `bin/adapt` and takes `--base-path`, which this interface needs to embed it. If it does not, ADapt is moved back to the commit it was on and the line says so, so a change to ADapt never leaves an installation without a working ADapt. The line always says which commit ADapt is at and why; the next one says whether ADCP 1.1 is installed and, when it is not, gives the command that adds it.

`ADAPT_REPO_URL` fetches from another address instead, a fork for example, and `ADAPT_REPO_BRANCH` follows another branch.

#### Where the interface finds it

At `$ADAPT_HOME`, in a folder named `adapt` beside ADFRsuite, in the installation's own `tools/docking`, beside this checkout or in a tools folder given with `--tools-root`, or at `~/.adcp-dock/adapt`. It needs ADapt from commit 6dae90b on, which takes `--base-path` and writes `_rec_H.pdb`. An older checkout is not offered, and the Receptor card says how to update it; `./update.sh` does.

#### In the page

When the interface finds ADapt and ADCP 1.1, ADapt opens by itself in a frame inside the Receptor card, on every structure loaded there: a PDB id fetched, a file on the cluster or one from your own computer. Opening it changes nothing. The line above the frame always says which receptor will be docked:

* **Docking will use this page's preparation.** Until you change anything in ADapt, the structure is prepared as the Receptor card's settings say, exactly as without ADapt.
* **Docking will use the receptor as set in ADapt.** From the first change there, an unticked chain or the waters taken out for example, ADapt decides what stays in the receptor: the assembly, hetero groups such as a catalytic zinc, metal charges, alternate locations and modified residues. The card's preparation settings are greyed out.

You do not need ADapt's own save button. When you submit, the page asks ADapt to write the receptor as it is set there (its `api/save`, the call its *receptor PDBQT* button makes), and reads back what ADapt wrote in its session folder: the receptor PDBQT, its `_rec_H.pdb` and `session.json`. The PDBQT is docked as ADapt wrote it, and this package's preparation does not run. The run copies the three files into its `receptor` folder and records the session, the save and their checksums. *Use ADapt's receptor now* does the same at once, so the site can be checked on the receptor that will be docked, and *Use this page's preparation instead* goes back. *Open full screen* shows the same session in a tab of its own; on a phone the card offers *Open in ADapt*, which opens it in a tab. From the command line: `adcp-dock new RUN --adapt-session FOLDER [--adapt-save N] -s ...`. ADapt keeps its sessions in `adapt_sessions` inside the runs folder, and it stops when the interface stops. The built-in Docs page *A receptor prepared in ADapt* shows each step.

ADapt runs beside the interface on the same machine, started with `--base-path /adapt`, and the interface's server passes `/adapt/` to it unchanged, as ADapt's README describes, so it reaches your browser over the same connection as the interface. Nothing else needs forwarding, whether you use the desktop launcher or your own `ssh -L`. Its first start imports AGFR's libraries and can take a minute or two on a busy login node, so the interface starts it in the background as soon as the page opens, and the card says what it is waiting for. To run ADapt on its own: `ADCPSUITE=/path/to/adcpsuite ADAPT_OUTDIR=/path/to/folder adapt/bin/adapt serve`.

### reduce and MSMS

The receptor is protonated by default with reduce, the program ADFR's own preparation and ADapt use, taken from ADFRsuite's `bin` folder (`reduce` 3.23 there). Checked against the protein structures in the PDB solved by neutron diffraction, where the hydrogens are placed from the data, reduce matched the deposited histidines more often than pdb2pqr (55% against 31% of 1,457 histidines) and turned fewer rings, and fewer Asn and Gln amides, the wrong way. What reduce does not do is pH: it never protonates an Asp or Glu and rarely puts hydrogens on both nitrogens of a histidine, whatever the pH. PDB2PQR with PROPKA at a chosen pH stays as the other choice (`--protonation propka --ph 5.5`), and a job file from before the choice keeps it, so old runs reproduce.

The site view and the pose viewer draw molecular surfaces computed by MSMS, the `bin/msms` program in ADFRsuite, which the interface's server calls on request. MSMS is free for non-commercial scientific research only, so it is called where ADFRsuite installed it and is not part of this package (see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)). Without it the viewers draw their own surfaces. `adcp-dock check` lists both programs.

---

## 5. The desktop launcher

The launcher does from your laptop what section 1.5 does by hand: it connects
to the cluster over SSH (two-factor sign in included), starts
`adcp-dock serve` on the login node, forwards a port and opens your browser.
Closing it stops the server; runs already submitted keep going in Slurm.

Install on the cluster first (section 1.2). The launcher only starts what is
installed there.

### 5.1 Get it

* **Windows and macOS**: download from the repository's **Releases** page.
  Windows is a folder: keep it together and run the `.exe` inside. macOS is
  an Apple silicon `.app`.
* **Any machine with Python**: `pip install "adcp-dock[launcher] @
  git+https://github.com/aadhilhaq/ADCP_docking"`, then `adcp-dock-launcher`.

The bundles are not yet code-signed, so the first run shows a warning:
Windows says *Windows protected your PC* (**More info → Run anyway**); macOS
refuses to open it (**right-click → Open → Open**). Signing removes those
prompts.

### 5.2 What it asks for

| field | what to put |
|---|---|
| **Cluster** | a preset, or **Other (SLURM)** for any Slurm cluster |
| **Address** | your cluster's login address, the one you use with `ssh` |
| **Username** | your cluster username |
| **Tool directory** | the folder you cloned into on the cluster, for example `/scratch/you/ADCP_docking`; `pwd` inside it prints it |
| **ADFRsuite / Amber roots** | leave blank: `activate.sh` in that folder names them |

Then press connect. It asks for your password, and for a second factor if
your cluster uses one, or signs in with your SSH key. The first start on a
busy login node can take a minute; the browser opens by itself when the page
is ready. Everything but the password is remembered for next time.

---

## 6. If something does not work

For problems while installing, see [1.2, If it stops](#if-it-stops).

| what you see | what it means |
|---|---|
| `adcp-dock: command not found` | this terminal has not loaded the tools yet: `source /path/to/ADCP_docking/activate.sh` (1.3) |
| `N tool(s) missing` from `check` | something the install recorded is not reachable any more; run `./install.sh` again in the checkout, which finds it or says what to do |
| `sbatch was not found` | Slurm is not available in this environment; load the required cluster environment or use `adcp-dock local` on a workstation |
| the downloaded ADFRsuite *does not match its pinned checksum* | the archive differs from the pinned version, usually an incomplete download; delete it as the message says and run the install again |
| `the account on the job is not one you can use` | this cluster needs an **Account** (1.4) |
| `the grid target is unusable` | the `.trg` was truncated or built from another structure. It is rebuilt once automatically; twice means agfr is failing, and its output is in the prepare log |
| `tleap did not build topologies` | the preparation log names the residue. It is usually a HETATM kept by *keep heteroatoms*, or a gap that *rebuild missing loops* would close |
| `the docking box also reaches chain(s) ...` | the box covers a chain that is not loaded. If it is the far wall of the groove, include it; if it is the partner your peptide should displace, leave it out |
| `pose N discarded: the bonded terms did not cancel` | the energy decomposition failed a consistency check; inspect the chain selections, topology, and rescoring log |
| the queue rejects the array | check the scheduler message for total array-size and concurrency limits; reduce the peptide list or adjust *Max tasks at once* as appropriate |

---

## 7. Uninstall

On a cluster, delete the folder you gave `install.sh` and the checkout:

```bash
rm -rf $SCRATCH/adcp $SCRATCH/ADCP_docking
```

If you pointed `install.sh` at an environment of your own (`ADCP_ENV`), the
packages it added there stay; remove them with that environment's own conda or
pip if you no longer want them.

On your own computer (section 2):

```bash
rm -rf ~/.adcp-dock
```

Runs are kept wherever you wrote them (by default `$SCRATCH/adcp_runs` on a
cluster and `~/adcp_runs` on your own computer). Move any you want to keep out
of the folders above before deleting them.
