"""Command line, for people who would rather not use the web page.

    adcp-dock check                              what can be found on this machine
    adcp-dock new  RUNDIR -r rec.pdb -s A:41,A:45 -p RWGPGAW -p KDGGVDG
    adcp-dock run  RUNDIR                        submit it to SLURM
    adcp-dock status RUNDIR
    adcp-dock report RUNDIR
    adcp-dock serve --port 8765                  the web page
    adcp-dock affinity [RUNDIR]                  what a target Kd is worth
"""
import argparse
import json
import os
import sys

from . import credit, fetch, job as jobmod, peptide, pipeline, slurm
from . import tools as toolmod


def _peptides(values, path):
    out = []
    # -p reads the same forms as a peptide file and the page, so that `-p
    # "p53: ETFSDLWKLLP"` names the peptide rather than being read as a
    # sequence containing a colon.
    bad = []
    for v in values or []:
        for name, seq, probs in peptide.parse_many(v):
            if probs:
                bad.append("  %s: %s" % (v, probs[0].message))
            named = ":" in v or v.lstrip().startswith(">")
            out.append({"name": name if named else "", "sequence": seq})
    if bad:
        raise SystemExit("-p has peptides that cannot be docked:\n%s"
                         % "\n".join(bad))
    if path:
        # Parse problems are reported rather than discarded: a file
        # that does not parse would otherwise become a job whose
        # peptides were its comment lines, caught later by validation
        # and far from its cause.
        bad = []
        for name, seq, probs in peptide.parse_many(open(path).read()):
            if probs:
                bad.append("  %s: %s" % (name, probs[0].message))
            out.append({"name": name, "sequence": seq})
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
    n.add_argument("-r", "--receptor", required=True,
                   help="a path, or a four-character PDB id to fetch")
    n.add_argument("-s", "--site", required=True,
                   help="residues the box should cover, e.g. A:41,A:45,A:78")
    n.add_argument("-p", "--peptide", action="append", default=[])
    n.add_argument("-f", "--peptides-file")
    n.add_argument("--chain", action="append", default=[])
    n.add_argument("--ph", type=float, default=7.4)
    n.add_argument("--no-minimise", dest="minimise",
                   action="store_false", default=True,
                   help="skip the restrained minimisation after "
                        "protonation (it is on by default)")
    n.add_argument("--keep-hetatms", action="store_true",
                   help="keep cofactors and other heteroatoms. Off by "
                        "default: they are dropped unless they are part of "
                        "the site being docked into")
    n.add_argument("--fill-loops", action="store_true",
                   help="rebuild whole missing loops, not only missing side "
                        "chains. Off by default, because a rebuilt loop is a "
                        "guess and the docking cannot tell it from a "
                        "measurement")
    n.add_argument("--box-chains", default="ignore",
                   help="what to do about a chain the box reaches into but "
                        "that was not loaded: 'ignore' (report it), "
                        "'include', or a list like 'B,C'")
    n.add_argument("--poses", type=int, default=10)
    n.add_argument("--replicas", type=int, default=50)
    n.add_argument("--steps", type=int, default=2_500_000)
    n.add_argument("--mmgbsa", action="store_true")
    n.add_argument("--mmgbsa-poses", type=int, default=5)
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
    n.add_argument("--cpus", type=int, default=48)
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
    ct.add_argument("--how", action="store_true",
                    help="how to make this tool itself citable (DOI, JOSS)")
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

    if a.cmd == "cite":
        from . import cite as citemod
        if a.verify:
            return citemod.verify()
        if a.how:
            citemod.how()
            return 0
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
        print("\nsbatch: %s" % ("found" if slurm.available() else "NOT FOUND"))
        if missing:
            print("\n%d tool(s) missing. Set ADFRSUITE_HOME and AMBERHOME, or "
                  "pass --tools-root." % missing)
        return 1 if missing else 0

    if a.cmd == "new":
        # A four-character code is fetched from the RCSB and cached
        # beside the run, so the same id reused on a node with no
        # outbound route still works. Anything else is treated as a
        # path.
        rec = a.receptor
        if fetch.looks_like_code(rec):
            cache = os.path.join(os.path.dirname(os.path.abspath(a.rundir)),
                                 ".pdb_cache")
            try:
                rec = fetch.fetch(rec, cache)
            except fetch.FetchError as exc:
                print("could not fetch %s: %s" % (a.receptor, exc))
                return 1
            print("fetched %s -> %s" % (a.receptor.upper(), rec))
        j = jobmod.default(
            receptor=os.path.abspath(rec), hotspots=a.site,
            chains=a.chain, ph=a.ph, minimise=a.minimise, poses=a.poses,
            replicas=a.replicas, steps=a.steps, mmgbsa=a.mmgbsa,
            mmgbsa_poses=a.mmgbsa_poses, top_n=a.top_n,
            max_length=a.max_length, partition=a.partition, cpus=a.cpus,
            time=a.time, throttle=a.throttle,
            account=a.account, sbatch_extra=a.sbatch_extra,
            keep_hetatms=a.keep_hetatms, fill_loops=a.fill_loops,
            box_chains=a.box_chains,
            tools_roots=[os.path.abspath(r) for r in a.tools_root],
            peptides=_peptides(a.peptide, a.peptides_file))
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

    if a.cmd == "affinity":
        from . import affinity
        spread = a.spread
        if spread is None and a.rundir:
            spread = affinity.spread_of(a.rundir)
        affinity.table(a.temperature, spread)
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
