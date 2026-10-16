"""Command line, for people who would rather not use the web page.

By Aadhil Haq

    adcp-dock check                              what can be found on this machine
    adcp-dock new  RUNDIR -r rec.pdb -s A:41,A:45 -p RWGPGAW -p KDGGVDG
    adcp-dock run  RUNDIR                        submit it to SLURM
    adcp-dock status RUNDIR
    adcp-dock report RUNDIR
    adcp-dock serve --port 8765                  the web page
    adcp-dock affinity [RUNDIR]                  what a target Kd is worth
    adcp-dock example                            a finished run to look at
"""
import argparse
import json
import os
import sys

from . import credit, fetch, job as jobmod, peptide, pipeline, slurm
from . import tools as toolmod


def _peptides(values, path, **kw):
    out = []
    # -p reads the same forms as a peptide file and the page, so that `-p
    # "p53: ETFSDLWKLLP"` names the peptide rather than being read as a
    # sequence containing a colon.
    bad = []
    for v in values or []:
        for r in peptide.parse_rows(v, **kw):
            if r["problems"]:
                bad.append("  %s: %s" % (v, r["problems"][0].message))
            p = {"name": r["name"] if r["named"] else "",
                 "sequence": r["sequence"]}
            if r["cyclic"]:
                p["cyclic"] = r["cyclic"]
            out.append(p)
    if bad:
        raise SystemExit("-p has peptides that cannot be docked:\n%s"
                         % "\n".join(bad))
    if path:
        # Parse problems are reported rather than discarded: a file
        # that does not parse would otherwise become a job whose
        # peptides were its comment lines, caught later by validation
        # and far from its cause. The ring column is read as the page
        # reads it, and a row that leaves it empty is docked as --cyclic
        # says.
        bad = []
        for r in peptide.parse_rows(open(path).read(), **kw):
            if r["problems"]:
                bad.append("  %s: %s" % (r["name"], " ".join(
                    x.message for x in r["problems"])))
            p = {"name": r["name"], "sequence": r["sequence"]}
            if r["cyclic"]:
                p["cyclic"] = r["cyclic"]
            out.append(p)
        if bad:
            raise SystemExit("%s has peptides that cannot be docked:\n%s"
                             % (path, "\n".join(bad)))
    for i, p in enumerate(out, 1):
        p["name"] = p["name"] or "peptide_%d" % i
    return out


def _local(a):
    """`adcp-dock local`: set this machine up, then open the interface.

    First use installs ADFRsuite (from Scripps, checksum-verified) and reports
    anything else that is missing; afterwards it starts in seconds. The work
    runs here rather than being submitted to SLURM.
    """
    import threading
    import webbrowser
    from . import adfr, slurm, server

    # Not on a cluster login node. A machine with sbatch that is not
    # inside a job is almost always one. On a cluster, `adcp-dock serve`
    # submits the work to SLURM instead.
    if slurm.available() and not os.environ.get("SLURM_JOB_ID") and not a.force:
        raise SystemExit(
            "This machine has sbatch, so it is probably a cluster login node, "
            "and dockings must not run on one.\n\n"
            "  On a cluster:  adcp-dock serve     (the work goes to SLURM)\n"
            "  Inside a job:  srun --pty -c 8 -t 4:00:00 adcp-dock local\n"
            "  Not a cluster: adcp-dock local --force")

    missing = []
    for mod, pkg in (("pdb2pqr", "pdb2pqr"), ("openmm", "openmm"),
                     ("pdbfixer", "pdbfixer"), ("rdkit", "rdkit")):
        try:
            __import__(mod)
        except ImportError:
            missing.append(pkg)
    if missing:
        raise SystemExit(
            "Missing: %s. Install everything this needs in one go:\n\n"
            "  pip install \"adcp-dock[local] @ "
            "git+https://github.com/aadhilhaq/ADCP_docking\"\n"
            % ", ".join(missing))

    tools = os.path.abspath(os.path.expanduser(a.tools))
    adfr.install(tools)                       # returns at once if installed
    os.environ["ADFRSUITE_HOME"] = tools
    os.environ["PATH"] = os.path.join(tools, "bin") + os.pathsep \
        + os.environ.get("PATH", "")
    os.environ["ADCP_EXECUTOR"] = "local"
    # The browser is on this machine, so the page may open a run's folder in
    # the file manager.
    os.environ["ADCP_DESKTOP"] = "1"

    runs = os.path.abspath(os.path.expanduser(a.runs))
    os.makedirs(runs, exist_ok=True)
    httpd = server.build("127.0.0.1", a.port, runs)
    url = "http://127.0.0.1:%d/?t=%s" % (httpd.server_address[1],
                                         httpd.token)
    from . import APP
    print("\n%s is running on this machine, %d core(s).\n\n    %s\n\n"
          "Runs are written to %s. Ctrl-C stops it; a run in progress stops "
          "with it.\nMM-GBSA needs AmberTools (conda install -c conda-forge "
          "ambertools); docking does not.\n"
          % (APP, slurm.usable_cores(), url, runs))
    if not a.no_browser:
        # Under WSL the browser is on the Windows side and localhost
        # is shared, so the address works there even when this cannot
        # open it.
        threading.Timer(1.0, lambda: webbrowser.open(url)).start()
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    return 0


def default_runs():
    """Where runs go when none is given.

    $SCRATCH on a cluster that provides it: a home directory usually has a
    small file allowance and a docking run writes tens of thousands of files.
    """
    for var in ("ADCP_RUNS", "SCRATCH"):
        root = os.environ.get(var)
        if root and os.path.isdir(root):
            return (root if var == "ADCP_RUNS"
                    else os.path.join(root, "adcp_runs"))
    return os.path.expanduser("~/adcp_runs")


def main(argv=None):
    ap = argparse.ArgumentParser(prog="adcp-dock", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")

    c = sub.add_parser("check", help="report which executables were found")
    c.add_argument("--tools-root", action="append", default=[])

    n = sub.add_parser("new", help="create a run directory")
    n.add_argument("rundir")
    n.add_argument("-r", "--receptor",
                   help="a path, or a four-character PDB id to fetch")
    n.add_argument("-s", "--site",
                   help="residues the box should cover, e.g. A:41,A:45,A:78")
    n.add_argument("--site-mode", choices=("residues", "ligand", "pocket"),
                   default="residues",
                   help="where the box comes from: the residues given with "
                        "-s (the default); a ligand or peptide in the "
                        "structure, named with --site-ligand; or the pocket "
                        "AutoSite finds (--site-pocket)")
    n.add_argument("--site-ligand", default="",
                   help="with --site-mode ligand: a peptide's chain, like B, "
                        "or a hetero group as CHAIN:NAME:NUMBER, like "
                        "A:NUT:201. The box is built around it and it is "
                        "kept out of the receptor")
    n.add_argument("--site-pocket", choices=("best", "all"), default="best",
                   help="with --site-mode pocket: AutoSite's best scoring "
                        "pocket (the default), or all the pockets it finds")
    n.add_argument("--agfr", choices=("auto", "1.1", "1.0"), default="auto",
                   help="which AGFR builds the grid maps: auto uses the "
                        "ADCP 1.1 suite's where it is installed")
    n.add_argument("--autosite", choices=("1.1", "1.0"), default="1.1",
                   help="which AutoSite places the translation points and "
                        "finds the pockets (1.1 by default)")
    n.add_argument("--receptor-gradient", action="store_true",
                   help="add the receptor gradient to the maps (AGFR's -g). "
                        "Off by default: it replaces the values inside the "
                        "receptor with a slope toward its surface, costs "
                        "time, and ADCP's authors found it did not help")
    n.add_argument("--target-file", default="",
                   help="a .trg target made elsewhere (for example in AGFRGUI): "
                        "its receptor, maps and box are used, and -r and -s "
                        "are not needed")
    n.add_argument("--engine", choices=("auto", "1.0", "1.1"), default="auto",
                   help="which ADCP: auto uses 1.1 where it is installed")
    n.add_argument("--omm-minimize", type=int, default=0,
                   help="ADCP 1.1: minimize this many top poses with OpenMM")
    n.add_argument("--omm-env", choices=("implicit", "vacuum"),
                   default="implicit")
    n.add_argument("-p", "--peptide", action="append", default=[])
    n.add_argument("-f", "--peptides-file",
                   help="a peptide file: one sequence per line, name: "
                        "SEQUENCE, FASTA, or columns (sequence, name and an "
                        "optional ring: linear, head to tail or disulfide), "
                        "with or without a header row")
    n.add_argument("--chain", action="append", default=[])
    n.add_argument("--protonation", choices=("reduce", "propka"),
                   default="reduce",
                   help="Protonate with: reduce (the default), ADFRsuite's, "
                        "as ADFR and ADapt protonate, which places histidine "
                        "protons from the hydrogen bonds and takes no pH; or "
                        "propka, PDB2PQR with PROPKA at the pH given with "
                        "--ph")
    n.add_argument("--ph", type=float, default=7.4,
                   help="the pH for --protonation propka (7.4 by default); "
                        "reduce takes none")
    n.add_argument("--adapt-session", default="",
                   help="Receptor from ADapt: a session folder ADapt wrote "
                        "(it holds session.json). Its receptor PDBQT is "
                        "docked as ADapt wrote it and its _rec_H.pdb is the "
                        "prepared receptor, so -r and this package's "
                        "preparation are not used")
    n.add_argument("--adapt-save", type=int, default=0,
                   help="with --adapt-session: which receptor save to use "
                        "(the newest by default)")
    n.add_argument("--no-minimise", dest="minimise",
                   action="store_false", default=True,
                   help="skip the restrained minimisation after "
                        "protonation (it is on by default)")
    n.add_argument("--keep-hetatms", action="store_true",
                   help="keep cofactors and other heteroatoms. Off by "
                        "default: they are dropped unless they are part of "
                        "the site being docked into")
    n.add_argument("--no-sidechains", dest="rebuild_sidechains",
                   action="store_false", default=True,
                   help="leave side chains as the structure has them: atoms "
                        "it lacks are not in the prepared receptor. Missing "
                        "side chain atoms are rebuilt by default")
    n.add_argument("--fill-loops", action="store_true",
                   help="rebuild whole missing loops as well. Off by "
                        "default, because a rebuilt loop is a guess and the "
                        "docking cannot tell it from a measurement")
    n.add_argument("--box-chains", default="ignore",
                   help="what to do about a chain the box reaches into but "
                        "that was not loaded: 'ignore' (report it), "
                        "'include', or a list like 'B,C'")
    n.add_argument("--poses", type=int, default=10,
                   help="ranked poses written per peptide: up to 10 with "
                        "ADCP 1.0 and up to 100 with ADCP 1.1 (adcp -m)")
    n.add_argument("--cyclic", choices=("none", "backbone", "disulfide"),
                   default="none",
                   help="dock every peptide linear (the default), cyclic "
                        "head to tail (adcp -cyc), or cyclic through a "
                        "disulfide between two cysteines (adcp -cys). A "
                        "peptide whose row in the -f file has a ring column "
                        "(linear, head to tail or disulfide) is docked with "
                        "that ring instead")
    n.add_argument("--helical", type=int, default=25,
                   help="the per cent of replicas that start as a helix, "
                        "the rest as coil (adcp -p). 25 by default; 0 starts "
                        "every replica as coil")
    n.add_argument("--keep-trajectories", action="store_true",
                   help="keep every replica's trajectory, packed into one "
                        "archive per peptide (adcp -k)")
    n.add_argument("--replicas", type=int, default=50)
    n.add_argument("--steps", type=int, default=2_500_000)
    n.add_argument("--mmgbsa", action="store_true")
    n.add_argument("--mmgbsa-poses", type=int, default=5)
    n.add_argument("--mmgbsa-model", choices=("obc2", "gbn2"), default="obc2",
                   help="the implicit solvent model for MM-GBSA: obc2, "
                        "GB-OBC(II), the default; or gbn2, GBn2 with its own "
                        "radii. Both run in AmberTools. Compare peptides only "
                        "within one model")
    n.add_argument("--mmgbsa-engine",
                   choices=("amber", "openmm", "schrodinger"),
                   default="amber",
                   help="which program computes the MM-GBSA energy: amber, "
                        "AmberTools, the default and the reference; openmm, "
                        "the same force field and topology solved again; or "
                        "schrodinger, Prime with OPLS and VSGB 2.0. Each "
                        "writes its own column. They are not comparable by "
                        "value, only by rank within one column")
    n.add_argument("--top-n", type=int, default=5)
    n.add_argument("--max-length", type=int, default=peptide.MAX_LENGTH)
    # No default partition. "cpu" is one cluster's name for it and not
    # another's, so defaulting it had runs refused before they started.
    # Empty means the scheduler's own default.
    n.add_argument("--partition", default="")
    n.add_argument("--account", default="",
                   help="the allocation to charge, where the cluster "
                        "requires one")
    n.add_argument("--sbatch-extra", default="",
                   help='other SLURM options, written with = so they are not read '
                        'as options of this command: '
                        '--sbatch-extra="--qos=normal --partition=short"')
    # 48, or what a node of the default partition has when that is fewer:
    # sbatch refuses a task asking for more cores than any node holds.
    n.add_argument("--cpus", type=int, default=None)
    n.add_argument("--time", default="08:00:00")
    n.add_argument("--throttle", type=int, default=0)
    n.add_argument("--tools-root", action="append", default=[])

    r = sub.add_parser("run", help="submit the run to SLURM")
    r.add_argument("rundir")
    r.add_argument("--restart", action="store_true",
                   help="submit every peptide again, including the ones "
                        "that already have a result")

    for name, helptext in (("status", "what the jobs are doing"),
                           ("report", "gather the results")):
        s = sub.add_parser(name, help=helptext)
        s.add_argument("rundir")

    # For a run docked before the interface panel existed. It measures the
    # panel from the files the MM-GBSA already left and docks nothing again.
    pn = sub.add_parser("panel",
                        help="measure the interface panel for a finished run "
                             "that has none, then rebuild its report")
    pn.add_argument("rundir")
    pn.add_argument("--here", action="store_true",
                    help="measure it on this machine instead of submitting it "
                         "to SLURM. It reads topologies and runs MMPBSA.py "
                         "once per pose, so this is for a workstation, not a "
                         "cluster login node.")

    rs = sub.add_parser("rescore",
                        help="score a finished run's poses again with another "
                             "MM-GBSA engine, then rebuild its report")
    rs.add_argument("rundir")
    rs.add_argument("--engine", required=True,
                    choices=("amber", "openmm", "schrodinger"),
                    help="the engine to score with. The run's own MM-GBSA "
                         "values, which it is ranked on, are left as they "
                         "are; this one's numbers appear beside them as "
                         "their own columns")
    rs.add_argument("--here", action="store_true",
                    help="run it on this machine instead of submitting it to "
                         "SLURM. It minimises every pose, so this is for a "
                         "workstation, not a cluster login node.")
    p = sub.add_parser("affinity",
                       help="what a target Kd is worth in kcal/mol")
    p.add_argument("rundir", nargs="?",
                   help="a finished run, to compare its own pose spread "
                        "against the window")
    p.add_argument("--temperature", type=float, default=298.15)
    p.add_argument("--spread", type=float, default=None)

    # Citation is a command rather than a README section: a method
    # section is written from the output directory months later, and the
    # tools omitted are consistently the docking engine and the
    # preparation steps.
    ct = sub.add_parser("cite", help="what to cite, for this tool and "
                                     "everything it runs")
    ct.add_argument("rundir", nargs="?",
                    help="a run directory: lists only what that run used")
    ct.add_argument("--verify", action="store_true",
                    help="resolve every DOI against Crossref (login node)")
    ct.add_argument("--bibtex", action="store_true", help="emit BibTeX")
    ct.add_argument("--ris", action="store_true",
                    help="emit RIS, which EndNote/Zotero/Mendeley import")

    lc = sub.add_parser(
        "local", help="set up and open the interface on this machine (a "
                      "laptop, a workstation or WSL), with no cluster")
    lc.add_argument("--tools", default="~/.adcp-dock/ADFRsuite",
                    help="where ADFRsuite is, or is installed on first use")
    lc.add_argument("--runs", default="~/adcp_runs")
    lc.add_argument("--port", type=int, default=8765)
    lc.add_argument("--no-browser", action="store_true",
                    help="print the address instead of opening it")
    lc.add_argument("--force", action="store_true",
                    help="run here even though this machine has sbatch")

    ia = sub.add_parser("install-adcp11",
                        help="install ADCP 1.1 (Scripps' conda recipe) so "
                             "runs can use it")
    ia.add_argument("--dest", default=None,
                    help="where to create it (default: beside ADFRsuite, in "
                         "the installation's tools folder)")

    ex = sub.add_parser("example",
                        help="copy a finished example run into the runs "
                             "folder, to open in the interface")
    ex.add_argument("--runs", default=None,
                    help="the runs folder, with the same default as serve")
    ex.add_argument("--force", action="store_true",
                    help="replace an example that is already there")

    s = sub.add_parser("serve", help="the web page")
    s.add_argument("--port", type=int, default=8765)
    s.add_argument("--host", default="127.0.0.1")
    s.add_argument("--runs", default=None,
                   help="where new runs go. Defaults to $SCRATCH/adcp_runs "
                        "on a cluster that has $SCRATCH and ~/adcp_runs "
                        "otherwise, and the page can move it afterwards")
    s.add_argument("--tools-root", action="append", default=[])
    s.add_argument("--open", action="store_true", help="print the URL only")

    a = ap.parse_args(argv)
    if not a.cmd:
        ap.print_help()
        return 2
    # The job scripts change into the run folder and then read job.json by
    # this path, so a relative one would be looked for inside itself.
    if getattr(a, "rundir", None):
        a.rundir = os.path.abspath(a.rundir)

    if a.cmd == "cite":
        from . import cite as citemod
        if a.verify:
            return citemod.verify()
        job = None
        if a.rundir:
            # No `import os` here: an import anywhere inside
            # main() makes the name local to all of main(), so
            # every other branch that uses the module-level os
            # raises UnboundLocalError.
            jp = os.path.join(a.rundir, "job.json")
            if not os.path.isfile(jp):
                print("no job.json in %s; showing everything instead\n"
                      % a.rundir)
            else:
                job = json.load(open(jp))
        if a.bibtex:
            print(citemod.bibtex(job), end="")
            return 0
        if a.ris:
            print(citemod.ris(job), end="")
            return 0
        citemod.text(job)
        return 0

    if a.cmd == "check":
        print(credit())
        print()
        T = toolmod.discover(roots=a.tools_root)
        missing = 0
        for name, path in T.report():
            print("  %-16s %s" % (name, path or "NOT FOUND  (%s)"
                                  % toolmod.NEEDED[name]))
            missing += 0 if path else 1
        for name, why in sorted(toolmod.ADFR_PROGRAMS.items()):
            path = T.adfr_program(name)
            if name == "msms" and not path:
                # The ADCP 1.1 suite has no msms program but MSMS's library.
                from . import surface as _surf
                path = _surf.describe(_surf.msms_command(T))
            print("  %-16s %s" % (name, path or "not found, optional (%s)"
                                  % why))
        # The Python packages, in the Python running this. The optional ones
        # turn a step off when they are missing, so they are named here.
        print("\nPython %s" % sys.executable)
        import importlib
        for mod, label, why in (
                ("numpy", "numpy", ""), ("rdkit", "RDKit", ""),
                ("openmm", "OpenMM", "the receptor's restrained minimisation"),
                ("pdbfixer", "PDBFixer", "rebuilding missing atoms and side chains"),
                ("Bio", "Biopython", "reading mmCIF files")):
            try:
                m = importlib.import_module(mod)
                state = getattr(m, "__version__", "") or "found"
                if mod == "openmm":
                    state = getattr(getattr(m, "version", None), "version", state)
            except Exception:                                   # noqa: BLE001
                state = ("not installed, optional (%s)" % why if why
                         else "NOT INSTALLED (required)")
                missing += 0 if why else 1
            print("  %-16s %s" % (label, state))
        print("\nsbatch: %s" % ("found" if slurm.available() else "NOT FOUND"))
        if missing:
            print("\n%d tool(s) missing. Set ADFRSUITE_HOME and AMBERHOME, or "
                  "pass --tools-root; or run ./install.sh again, which finds "
                  "them, and a cluster's AmberTools module once it is loaded."
                  % missing)
        return 1 if missing else 0

    if a.cmd == "new":
        # A four-character code is fetched from the RCSB and cached
        # beside the run, so the same id reused on a node with no
        # outbound route still works. Anything else is treated as a
        # path.
        adapt_rec = None
        if a.adapt_session:
            from . import adapt as adaptmod
            try:
                adapt_rec = adaptmod.read_save(a.adapt_session,
                                               a.adapt_save or None)
            except adaptmod.NoSave as exc:
                print("--adapt-session: %s" % exc)
                return 2
            a.receptor = adapt_rec["pdb"]
            print("receptor from ADapt: save %d (%s) of %s"
                  % (adapt_rec["save"], adapt_rec["tab"] or "tab not recorded",
                     adapt_rec["folder"]))
        if not a.target_file and not (
                a.receptor and (a.site or a.site_mode != "residues")):
            print("give -r and -s, or -r and --site-mode ligand or pocket, "
                  "or --adapt-session and -s, or --target-file")
            return 2
        if a.site_mode == "ligand" and not a.site_ligand \
                and not a.target_file:
            print("--site-mode ligand needs --site-ligand: the peptide's "
                  "chain, like B, or a hetero group, like A:NUT:201")
            return 2
        rec = a.receptor or ""
        if rec and fetch.looks_like_code(rec):
            cache = os.path.join(os.path.dirname(os.path.abspath(a.rundir)),
                                 ".pdb_cache")
            try:
                rec = fetch.fetch(rec, cache)
            except fetch.FetchError as exc:
                print("could not fetch %s: %s" % (a.receptor, exc))
                return 1
            print("fetched %s -> %s" % (a.receptor.upper(), rec))
        roots = [os.path.abspath(r) for r in a.tools_root]
        engine, exe = jobmod.resolve_engine({"adcp_engine": a.engine,
                                             "tools_roots": roots})
        kw = {}
        if engine == "1.1":
            from . import adcp11
            kw = {"extended": True, "library": adcp11.library(exe)}
        j = jobmod.default(
            receptor=os.path.abspath(rec) if rec else "", hotspots=a.site or "",
            target_file=os.path.abspath(a.target_file) if a.target_file else "",
            adcp_engine=a.engine, omm_minimize=a.omm_minimize,
            omm_env=a.omm_env,
            chains=a.chain, ph=a.ph, minimise=a.minimise, poses=a.poses,
            replicas=a.replicas, steps=a.steps, mmgbsa=a.mmgbsa,
            mmgbsa_poses=a.mmgbsa_poses, mmgbsa_model=a.mmgbsa_model,
            mmgbsa_engine=a.mmgbsa_engine,
            top_n=a.top_n,
            max_length=a.max_length, partition=a.partition,
            cpus=a.cpus if a.cpus else slurm.cores_for_cluster(48),
            time=a.time, throttle=a.throttle,
            account=a.account, sbatch_extra=a.sbatch_extra,
            keep_hetatms=a.keep_hetatms, fill_loops=a.fill_loops,
            rebuild_sidechains=a.rebuild_sidechains,
            box_chains=a.box_chains,
            site_mode=a.site_mode, site_ligand=a.site_ligand,
            site_pocket=a.site_pocket, agfr=a.agfr, autosite=a.autosite,
            receptor_gradient=a.receptor_gradient, cyclic=a.cyclic,
            helical=a.helical, keep_trajectories=a.keep_trajectories,
            protonation=a.protonation,
            adapt=({"folder": adapt_rec["folder"], "save": adapt_rec["save"]}
                   if adapt_rec else {}),
            tools_roots=roots,
            peptides=_peptides(a.peptide, a.peptides_file, **kw))
        if adapt_rec:
            jobmod.take_adapt(j)
        # A ring the peptide cannot close is refused here, before any job
        # file is written, with the reason: ADCP itself would dock it as a
        # linear peptide and say nothing. Each peptide with its own ring,
        # from the file's ring column, or with --cyclic's.
        bad = []
        for p in j["peptides"]:
            ring = jobmod.ring_of(j, p)
            pr = (peptide.cyclic_problem(p["sequence"], ring)
                  if ring in peptide.CYCLIC
                  else peptide.ring_problem(p.get("cyclic")))
            if pr:
                bad.append("  %s (%s), %s: %s"
                           % (p["name"], p["sequence"],
                              peptide.RING_AS.get(ring, ring), pr.message))
        if bad:
            raise SystemExit("%d peptide(s) cannot be docked with the ring "
                             "asked for (--cyclic %s, or the file's ring "
                             "column):\n%s"
                             % (len(bad), a.cyclic, "\n".join(bad)))
        problems = jobmod.validate(j)
        jobmod.save(a.rundir, j)
        print("wrote %s" % jobmod.path_of(a.rundir))
        for p in problems:
            print("  ! %s" % p)
        return 1 if problems else 0

    if a.cmd == "run":
        # A run that already has results is resumed rather than
        # restarted: resubmitting the whole array would repeat the
        # docking, which is the expensive half.
        todo = pipeline.outstanding(a.rundir)
        j = jobmod.load(a.rundir)
        n = len(jobmod.peptide_list(j))
        if a.restart or len(todo) == n:
            pipeline.submit(a.rundir)
        else:
            print("%d of %d peptide(s) already have a result; submitting "
                  "only the rest" % (n - len(todo), n))
            pipeline.resume(a.rundir)
        return 0

    if a.cmd == "status":
        j = jobmod.load(a.rundir)
        print("state: %s" % j.get("state"))
        for k, v in (j.get("jobs") or {}).items():
            print("  %-8s %-10s %s" % (k, v, slurm.state(v, a.rundir) or "(gone)"))
        return 0

    if a.cmd == "report":
        pipeline.stage_report(a.rundir)
        return 0

    if a.cmd == "rescore":
        if a.here:
            pipeline.stage_rescore(a.rundir, a.engine)
            pipeline.stage_report(a.rundir)
        else:
            got = pipeline.submit_rescore(a.rundir, a.engine)
            print("submitted: rescore %s, report %s"
                  % (got.get("rescore"), got.get("report")))
        return 0
    if a.cmd == "panel":
        if a.here:
            pipeline.stage_panel(a.rundir)
            pipeline.stage_report(a.rundir)
        else:
            pipeline.submit_panel(a.rundir)
        return 0

    if a.cmd == "affinity":
        from . import affinity
        spread = a.spread
        if spread is None and a.rundir:
            spread = affinity.spread_of(a.rundir)
        affinity.table(a.temperature, spread)
        return 0

    if a.cmd == "install-adcp11":
        from . import adcp11_install
        adcp11_install.install(a.dest or adcp11_install.default_dest())
        return 0

    if a.cmd == "example":
        from . import example
        path = example.install(a.runs or default_runs(), force=a.force)
        print("the example run is in %s; open the interface and choose it "
              "from the list of runs" % path)
        return 0

    if a.cmd in ("local", "serve"):
        # Line-buffered, whoever is reading. Both commands print the
        # address and then block serving the page; with stdout on a
        # pipe (the launcher, a notebook, a log file) Python would
        # hold that line in its buffer indefinitely.
        try:
            sys.stdout.reconfigure(line_buffering=True)
            sys.stderr.reconfigure(line_buffering=True)
        except (AttributeError, ValueError):
            pass

    if a.cmd == "local":
        return _local(a)

    if a.cmd == "serve":
        from . import server
        server.serve(host=a.host, port=a.port, runs=a.runs or default_runs(),
                     tools_roots=a.tools_root)
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
