"""A run, on disk: one directory, one job.json, everything else beside it.

By Aadhil Haq

A plain file rather than a database, so that a run survives the web app
being closed, can be read by anyone asked to look at it later, and can be
resumed by re-running the same command.
"""
import json
import os
import time

VERSION = 1
NAME = "job.json"
# ADCP builds side chains, and so writes an all-atom ranked pose, for its
# ten best clusters only (clusterADCP.py), so a run cannot keep more.
MAX_POSES = 10


def _max_length():
    from . import peptide as pep
    return pep.MAX_LENGTH


def default(**kw):
    j = {
        "version": VERSION,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "label": "",
        # the computing facility's required acknowledgement, if any
        "acknowledgement": os.environ.get("ADCP_ACKNOWLEDGEMENT", "").strip(),
        # receptor
        "receptor": "",          # path to the input structure
        "chains": [],            # chains to keep; empty means all
        "ph": 7.4,
        # Full Amber preparation is the default rather than an extra.
        # Protonating at a chosen pH sets the charges, and minimising
        # afterwards relaxes the hydrogens just placed along with the
        # clashes that come with them.
        "minimise": True,
        "minimise_steps": 500,
        "restraint": 10.0,       # kcal/mol/A^2 on the heavy atoms
        "forcefield": "AMBER",
        "keep_hetatms": False,   # cofactors: out unless they are the site
        "fill_loops": False,     # rebuilding whole missing loops is a guess
        # site
        "hotspots": "",          # "A:41, A:45" or "41,45"
        "shell_radius": 10.0,    # what counts as "near the site"
        "box_padding": 8.0,
        "box_min": 20.0,
        "box_max": 30.0,
        "spacing": 0.375,
        # A box built around an interface groove reaches into the
        # chain on the far side. Including it means docking against a
        # complete wall; omitting it means docking against a gap.
        # Which is correct depends on whether that chain is part of
        # the site or the partner being displaced.
        "box_chains": "ignore",
        "box_chains_min_atoms": 20,
        # peptides
        "peptides": [],          # [{name, sequence}]
        # The same limit the page and the peptide checks use, so a peptide
        # the page accepts is not refused at submission.
        "max_length": _max_length(),
        # docking
        "replicas": 50,
        "steps": 2_500_000,
        "poses": 10,             # how many ranked poses to keep per peptide
        "seed": None,
        # Which ADCP: "auto" uses 1.1 where it is installed and 1.0
        # (ADFRsuite) otherwise; "1.0" or "1.1" asks for one.
        "adcp_engine": "auto",
        # ADCP 1.1 only: minimize this many top poses with OpenMM after
        # docking, in "implicit" solvent or "vacuum". 0 does not minimize.
        "omm_minimize": 0,
        "omm_env": "implicit",
        # A target file (.trg) made elsewhere, for example in AGFRGUI. When
        # set, the receptor, maps and box come from it and this package's
        # receptor preparation is not run.
        "target_file": "",
        # rescoring
        "mmgbsa": False,
        # The generalised Born model: "obc2", GB-OBC(II), what every run
        # before it used, or "gbn2", Amber's GBn2 with its own radii. Both
        # run in AmberTools. A run keeps the model its MM-GBSA started with.
        "mmgbsa_model": "obc2",
        # Which program computes the energy. AmberTools is the default
        # because every number this project has published came from it, so a
        # run that does not ask for anything else stays comparable with them.
        "mmgbsa_engine": "amber",
        "mmgbsa_poses": 5,
        # How the poses of one peptide become a single number. The
        # table reports the mean of the top few either way, with the
        # best, the median, the count and the spread.
        "mmgbsa_aggregate": "mean",
        "mm_workers": 8,       # rescore this many of the ranked poses
        "top_n": 5,              # the mean reported is over this many
        # Resources. Blank deliberately: asked of the cluster at
        # submit time. Defaulting the partition fails on any cluster
        # that names its partitions differently.
        "partition": "",
        "cpus": 48,
        "mem": "32G",
        "time": "02:00:00",
        "mm_cpus": 8,
        "mm_mem": "32G",
        "mm_time": "12:00:00",
        "throttle": 0,           # max array tasks at once; 0 = no limit
        "tools_roots": [],
        "sbatch_extra": "",      # --qos=... --partition=... for this site
        "account": "",           # --account, where the cluster needs one
        # filled in as the run proceeds
        "state": "new",
        "jobs": {},
        "site": {},
    }
    j.update(kw)
    return j


def path_of(rundir):
    return os.path.join(rundir, NAME)


def save(rundir, j):
    os.makedirs(rundir, exist_ok=True)
    tmp = path_of(rundir) + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(j, fh, indent=1)
    os.replace(tmp, path_of(rundir))
    return path_of(rundir)


def load(rundir):
    p = path_of(rundir) if os.path.isdir(rundir) else rundir
    with open(p) as fh:
        return json.load(fh)


def resolve_engine(j):
    """("1.0" | "1.1" | "missing", ADCP 1.1 executable or None) for a job."""
    from . import adcp11
    want = str(j.get("adcp_engine") or "auto")
    if want == "1.0":
        return "1.0", None
    exe = adcp11.find(j.get("tools_roots") or [])
    if exe:
        return "1.1", exe
    return ("missing" if want == "1.1" else "1.0"), None


def key(sequence):
    """The name a peptide's docking folder and result file get.

    The sequence itself for ordinary peptides, so existing runs keep their
    names; ADCP 1.1 sequences carry <, > and &, which are replaced.
    """
    import re
    return sequence if sequence.isalpha() else re.sub(r"[^A-Za-z0-9]+", "_",
                                                      sequence).strip("_")


#: The engine names the registry knows. Read from it rather than written out
#: again, so this cannot drift from what can actually run.
def _gb_engine_names():
    try:
        from ._vendor import gbengines
        return tuple(gbengines.ORDER)
    except Exception:                                          # noqa: BLE001
        return ("amber", "openmm", "schrodinger")


GB_ENGINES = _gb_engine_names()


def validate(j):
    """[str] of things that would stop the run, in the order worth fixing."""
    from . import peptide as pep
    from . import slurm as _sl
    out = []
    # The site's own options are checked first: a job script the
    # scheduler rejects fails the run at its first stage, far from the
    # field where the option was entered.
    out.extend(_sl.site_options(j.get("account"), j.get("sbatch_extra"))[1])
    trg = str(j.get("target_file") or "").strip()
    if trg:
        from . import target as _t
        why = _t.check(trg) if os.path.isfile(trg) else "does not exist"
        if why:
            out.append("The target file %s is unusable: %s" % (trg, why))
    else:
        if not j.get("receptor"):
            out.append("No receptor structure chosen.")
        elif not os.path.isfile(j["receptor"]):
            out.append("The receptor file does not exist: %s" % j["receptor"])
        if not str(j.get("hotspots") or "").strip():
            out.append("No site chosen. Pick the residues the box should cover.")
    engine, exe = resolve_engine(j)
    if engine == "missing":
        out.append("ADCP 1.1 was asked for but is not installed here. Install "
                   "it (see INSTALL.md) or choose ADCP 1.0.")
    lib = None
    if engine == "1.1" and exe:
        from . import adcp11 as _a
        lib = _a.library(exe)
    if int(j.get("omm_minimize") or 0) and engine != "1.1":
        out.append("OpenMM minimization of poses needs ADCP 1.1.")
    if not j.get("peptides"):
        out.append("No peptides to dock.")
    for i, p in enumerate(j.get("peptides") or [], 1):
        seq, probs = pep.check(p.get("sequence", ""), extended=engine == "1.1",
                               library=lib,
                               max_length=j.get("max_length", pep.MAX_LENGTH))
        for pr in probs:
            out.append("Peptide %d (%s): %s"
                       % (i, p.get("name") or seq[:10] or "?", pr.message))
    # Duplicate sequences are fatal rather than untidy. A docking is
    # written to docking/<sequence>/ and its result to
    # results/<sequence>.json, so two entries with the same sequence
    # share a working directory and a result file.
    seqs, names = {}, {}
    for i, p in enumerate(j.get("peptides") or [], 1):
        seq = pep.clean(p.get("sequence", ""))
        label = (p.get("name") or "").strip()
        if seq:
            if seq in seqs:
                out.append(
                    "Peptides %d and %d are the same sequence (%s). Two "
                    "entries with one sequence would share a working "
                    "directory and a result file, so one would overwrite "
                    "the other. Remove one, or change the sequence."
                    % (seqs[seq], i, seq))
            else:
                seqs[seq] = i
        if label:
            key = label.lower()
            if key in names:
                out.append(
                    "Peptides %d and %d are both named %r. The name is what "
                    "labels the row in the results table, and two rows with "
                    "one name cannot be told apart."
                    % (names[key], i, label))
            else:
                names[key] = i
    if int(j.get("poses", 0)) < 1:
        out.append("Ask for at least one pose per peptide.")
    elif int(j.get("poses", 0)) > MAX_POSES:
        out.append("ADCP writes at most %d ranked poses per peptide; ask for "
                   "%d or fewer." % (MAX_POSES, MAX_POSES))
    if j.get("mmgbsa_model", "obc2") not in ("obc2", "gbn2"):
        out.append("The MM-GBSA model must be obc2 (GB-OBC(II)) or gbn2 "
                   "(GBn2), not %r." % j.get("mmgbsa_model"))
    if j.get("mmgbsa_engine", "amber") not in GB_ENGINES:
        out.append("The MM-GBSA engine must be one of %s, not %r."
                   % (", ".join(GB_ENGINES), j.get("mmgbsa_engine")))
    if j.get("mmgbsa") and int(j.get("mmgbsa_poses", 0)) < 1:
        out.append("MM-GBSA is on but no poses were chosen to rescore.")
    if j.get("mmgbsa") and int(j.get("mmgbsa_poses", 0)) > int(j.get("poses", 0)):
        out.append("MM-GBSA asks for %d poses but only %d are being kept."
                   % (j["mmgbsa_poses"], j["poses"]))
    return out


def peptide_list(j):
    """[(name, sequence)] cleaned, in the order given.

    A blank name becomes peptide_N by position, which is stable for a given job
    file: the name appears in the results table, and one that moved between
    runs would make two reports of the same job disagree.

    Duplicates are not resolved here; validate() refuses them.
    """
    from . import peptide as pep
    out = []
    for i, p in enumerate(j.get("peptides") or [], 1):
        seq = pep.clean(p.get("sequence", ""))
        out.append(((p.get("name") or "").strip() or "peptide_%d" % i, seq))
    return out
