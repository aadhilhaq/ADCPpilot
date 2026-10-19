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
# ADCP 1.1 writes every pose it reports with its side chains, up to -m, and
# its own default for -m is 100. ADCP 1.0, on a Mac, writes ten (mac10).
MAX_POSES = 100


def max_poses(engine):
    """The most ranked poses the engine writes per peptide."""
    from . import mac10
    return mac10.MAX_POSES if engine == mac10.ENGINE else MAX_POSES


#: The share of replicas that start as a helix in a new run, in per cent.
HELICAL_DEFAULT = 25
#: What a job file's adcp_engine may say. Every run docks with ADCP 1.1:
#: "auto", which new runs said before, means it too.
ENGINES = ("1.1", "auto")


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
        # Side chain atoms the structure lacks are rebuilt by default, which
        # is what every run made before this setting existed did, so a job
        # file without it reads as on. Off, the receptor that is docked into
        # keeps its side chains as the structure has them.
        "rebuild_sidechains": True,
        "fill_loops": False,     # rebuilding whole missing loops is a guess
        # How the receptor is protonated: "reduce", the ADCP 1.1 suite's, as
        # ADFR's own preparation and ADapt do it, or "propka", PDB2PQR with
        # PROPKA at the pH above, which the pH applies to. A job file without the
        # setting was protonated with PDB2PQR and PROPKA, and still is.
        "protonation": "reduce",
        # A receptor prepared in ADapt: its session folder and the receptor
        # save to use ({"folder": ..., "save": n}, with what was read back
        # from it). Its PDBQT is docked as ADapt wrote it, its _rec_H.pdb is
        # the prepared receptor, and this package's preparation is not run.
        # Empty for a receptor prepared here.
        "adapt": {},
        # site
        # Where the box comes from. "residues" builds it around the residues
        # below, which is what every run made before the choice existed did,
        # so a job file without it reads as that. "ligand" builds it around a
        # ligand or peptide in the structure, named in site_ligand as a chain
        # ("B") or as CHAIN:NAME:NUMBER ("A:NUT:201"). "pocket" builds it
        # around the pocket AutoSite finds: the best scoring one, or every
        # one it finds merged (site_pocket).
        "site_mode": "residues",
        "site_ligand": "",
        # The structure the ligand is read from, when it is not the receptor
        # docked into: with ADapt's receptor, the structure ADapt was opened
        # on, which still holds a ligand ADapt left out of its receptor.
        "site_ligand_from": "",
        "site_pocket": "best",
        # Pockets chosen on the page (pockets.py): the ranks of AutoSite's
        # pockets the search starts from, their points saved in the run as
        # site_points, and for a box built around a pocket, that box.
        "site_pockets": {},
        # How AutoSite ranks the pockets the page finds: its peptide score
        # (energy x buriedness^1.5) or its general one.
        "pocket_score": "peptide",
        "site_points": "",
        "site_box": None,
        "hotspots": "",          # "A:41, A:45" or "41,45"
        "shell_radius": 10.0,    # what counts as "near the site"
        "box_padding": 8.0,
        "box_min": 20.0,
        "box_max": 30.0,
        "spacing": 0.375,
        # The grid maps, built by the ADCP 1.1 suite's AGFR ("auto", which
        # new runs said before, means it too). A job file without the setting
        # was built by an older AGFR with its defaults, AutoSite 1.0 and the
        # receptor gradient on, and its maps are built again with those
        # (target.maps_for). New runs say all three, so a rebuilt target
        # always matches the one it replaces.
        "agfr": "1.1",
        "autosite": "1.1",
        "receptor_gradient": False,
        # A box built around an interface groove reaches into the
        # chain on the far side. Including it means docking against a
        # complete wall; omitting it means docking against a gap.
        # Which is correct depends on whether that chain is part of
        # the site or the partner being displaced.
        "box_chains": "ignore",
        "box_chains_min_atoms": 20,
        # peptides
        # [{name, sequence}], and "cyclic" on a peptide docked as its own
        # ring rather than as the run's (ring_of). A peptide without it, as
        # in every job file written before a peptide could carry one, is
        # docked as the run's "cyclic" below says.
        "peptides": [],
        # The same limit the page and the peptide checks use, so a peptide
        # the page accepts is not refused at submission.
        "max_length": _max_length(),
        # Linear ("none"), or closed head to tail ("backbone", ADCP's -cyc)
        # or through a disulfide between two cysteines ("disulfide", -cys).
        # The run's choice, which every peptide takes unless it has its own.
        "cyclic": "none",
        # The per cent of replicas that start as a helix, the rest as coil
        # (ADCP's -p). A job file without it was docked with every replica
        # starting as coil, so it reads as 0.
        "helical": HELICAL_DEFAULT,
        # docking
        "replicas": 50,
        "steps": 2_500_000,
        "poses": 10,             # how many ranked poses to keep per peptide
        "seed": None,
        # Every replica's trajectory, packed into one archive per peptide
        # (ADCP's -k). Off, they are not kept. A job file without it leaves
        # whatever ADCP left, as runs made before it did.
        "keep_trajectories": False,
        # Which ADCP: 1.1, the only one (ENGINES).
        "adcp_engine": "1.1",
        # Minimize this many top poses with ADCP's own OpenMM after docking,
        # in "implicit" solvent or "vacuum". 0 does not minimize.
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
        # because runs made before other engines could be chosen used it, so a
        # run that does not ask for anything else stays comparable with them.
        "mmgbsa_engine": "amber",
        "mmgbsa_poses": 5,
        # D residues and rings are rescored: a D residue with its L form's
        # template and the docked handedness, a ring only on the poses that
        # close it (mmforms). A job file without it was made when both were
        # skipped, and keeps skipping them, so one run never ranks peptides
        # scored under two rules.
        "mmgbsa_d_rings": True,
        # How the poses of one peptide become a single number. The
        # table reports the mean of the top few either way, with the
        # best, the median, the count and the spread.
        "mmgbsa_aggregate": "mean",
        "mm_workers": 8,       # rescore this many of the ranked poses
        "top_n": 5,              # the mean reported is over this many
        # Which order the results table is written in. Off keeps the
        # single-column order every run before this one used: MM-GBSA
        # where it ran, ADCP's own score where it did not. On orders
        # the table on the combined metric, which also takes in the
        # interface the panel measured. The combined score and its
        # components are columns either way, and results_combined.tsv
        # is written beside results.tsv either way, so this changes
        # the order of one file and never what the run reports.
        "rank_combined": False,
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
    """("1.1", the ADCP 1.1 executable) for a job, or ("missing", None)
    where the ADCP 1.1 suite is not installed.

    Every job docks with ADCP 1.1. A job file from before it was the only
    engine may name another in adcp_engine; what such a run still has to
    dock, on a resume, is docked with 1.1 as well.

    The one exception is a Mac, where Scripps builds no ADCP 1.1: there a
    run made on the Mac docks with ADFRsuite 1.0 (mac10.ENGINE, its folder),
    and one made elsewhere cannot dock ("missing"). A run made on a Mac
    cannot dock anywhere else ("mac-only"): its peptides would be docked by
    two engines.
    """
    from . import adcp11, mac10
    want = str(j.get("adcp_engine") or "1.1")
    if mac10.engaged():
        return ((mac10.ENGINE, mac10.home()) if want == mac10.ENGINE
                else ("missing", None))
    if want == mac10.ENGINE:
        return "mac-only", None
    exe = adcp11.find(j.get("tools_roots") or [])
    return ("1.1", exe) if exe else ("missing", None)


def for_this_machine(j):
    """A new job as this machine docks it: on a Mac, ADCP 1.0, named in
    adcp_engine and engine, with its maps from ADFRsuite 1.0's agfr and
    AutoSite 1.0; elsewhere left as it is (ADCP 1.1)."""
    from . import mac10
    if mac10.engaged():
        j["adcp_engine"] = mac10.ENGINE
        j["engine"] = mac10.LABEL
        # Its maps: ADFRsuite 1.0's agfr, which has AutoSite 1.0 alone.
        j["agfr"] = mac10.ENGINE
        j["autosite"] = "1.0"
    return j


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


def site_problems(j):
    """[str] why the job's site cannot place a box."""
    from . import structure, target
    mode = str(j.get("site_mode") or "residues")
    if mode not in target.SITE_MODES:
        return ["The site must come from %s, not %r."
                % (", ".join(target.SITE_MODES), mode)]
    if mode == "residues":
        if not str(j.get("hotspots") or "").strip():
            return ["No site chosen. Pick the residues the box should cover."]
        return []
    if mode == "pocket":
        pocket = str(j.get("site_pocket") or "best")
        if pocket not in target.POCKET_CHOICES:
            return ["The pocket must be %s, not %r."
                    % (" or ".join(target.POCKET_CHOICES), pocket)]
        return []
    spec = str(j.get("site_ligand") or "").strip()
    if not spec:
        return ["No ligand chosen. Pick the ligand or peptide in the structure "
                "the box should be built around."]
    rec = j.get("site_ligand_from") or j.get("receptor") or ""
    if not os.path.isfile(rec):
        return []
    try:
        atoms = structure.read(rec)
    except Exception as exc:                                    # noqa: BLE001
        return ["The receptor could not be read: %s" % exc]
    # With ADapt's receptor, ADapt chose what it holds: a ligand it kept is
    # taken out of it before the maps are built (pipeline._drop_ligand), so
    # neither the chains ticked here nor Keep heteroatoms decide that.
    adapt = bool((j.get("adapt") or {}).get("folder"))
    return structure.ligand_problems(atoms, spec,
                                     [] if adapt else (j.get("chains") or []),
                                     False if adapt else bool(j.get("keep_hetatms")))


def protonation_of(j):
    """'reduce' or 'propka' for a job: what it asks for, or what a job file
    written before the choice existed was protonated with."""
    from ._vendor import prepare
    if (j.get("adapt") or {}).get("folder"):
        return "reduce"
    return prepare.protonation_of(j)


def adapt_problems(j):
    """[str] why the job's ADapt receptor cannot be used."""
    a = j.get("adapt") or {}
    if not a.get("folder"):
        return []
    from . import adapt as adaptmod
    try:
        got = adaptmod.read_save(a["folder"], a.get("save"))
    except adaptmod.NoSave as exc:
        return ["The receptor from ADapt cannot be used: %s." % str(exc).rstrip(".")]
    out = []
    for k in ("pdbqt_md5", "pdb_md5"):
        if a.get(k) and a[k] != got.get(k):
            out.append("The receptor ADapt saved (%s, save %s) has changed since "
                       "this run was set up. Read it back again."
                       % (a["folder"], got.get("save")))
            break
    return out


def take_adapt(j):
    """A receptor from ADapt, read back and recorded in the job: its save,
    the md5 of its two files so a later change is caught, and the
    _rec_H.pdb as the receptor the site and the viewers read. Returns [str]
    of problems, empty when it is fine or there is none."""
    a = j.get("adapt") or {}
    if not a.get("folder"):
        j["adapt"] = {}
        return []
    from . import adapt as adaptmod
    try:
        got = adaptmod.read_save(a["folder"], a.get("save"))
    except adaptmod.NoSave as exc:
        return ["The receptor from ADapt cannot be used: %s"
                % str(exc).rstrip(".")]
    j["adapt"] = {k: got[k] for k in (
        "folder", "save", "tab", "pdbid", "input", "pdbqt", "pdb",
        "pdbqt_md5", "pdb_md5", "reduce", "biomolecule", "hetero",
        "adapt_version", "adapt_commit")}
    if a.get("session"):
        j["adapt"]["session"] = a["session"]
    j["receptor"] = got["pdb"]
    j["protonation"] = "reduce"
    return []


def map_problems(j):
    """[str] why the grid maps cannot be built as the job asks."""
    from . import mac10, target
    out = []
    agfr = str(j.get("agfr") or "1.1")
    if agfr == mac10.ENGINE and mac10.engaged():
        pass                    # ADFRsuite 1.0's agfr, on a Mac (mac10)
    elif agfr not in target.AGFR_CHOICES:
        out.append("The maps are built by AGFR 1.1, from the ADCP 1.1 suite: "
                   "agfr must be %s, not %r."
                   % (" or ".join(target.AGFR_CHOICES), agfr))
    autosite = str(j.get("autosite") or "1.0")
    if autosite not in target.AUTOSITE_CHOICES:
        out.append("AutoSite must be %s, not %r."
                   % (" or ".join(target.AUTOSITE_CHOICES), autosite))
    return out


def ring_of(j, p):
    """How one of the job's peptides is docked: "none", "backbone" or
    "disulfide".

    The peptide's own "cyclic" when it has one, and the run's otherwise,
    which is all a job file written before a peptide could carry its own
    ever has. A peptide's value is read as a person writes it ("head to
    tail", "linear"), and one that names no ring is returned as written, for
    validate() to refuse with the peptide's name.
    """
    from . import peptide as pep
    own = p.get("cyclic") if isinstance(p, dict) else None
    if own is not None and str(own).strip():
        got = pep.ring_choice(own)
        if got:
            return got
        if got is None:
            return str(own)
    return str((j or {}).get("cyclic") or "none")


def rings(j):
    """[ring] for the job's peptides, in the order peptide_list gives them."""
    return [ring_of(j, p) for p in (j.get("peptides") or [])]


def ring_for(j, seq):
    """The ring the peptide with this cleaned sequence is docked with.

    Results are kept by sequence, and a run never holds one sequence twice
    (validate), so the sequence is enough to find the peptide. One the job
    no longer lists is given the run's ring.
    """
    from . import peptide as pep
    for p in (j or {}).get("peptides") or []:
        if pep.clean(p.get("sequence", "")) == seq:
            return ring_of(j, p)
    return str((j or {}).get("cyclic") or "none")


def tidy_peptides(peptides):
    """The peptides as a job file keeps them: a peptide's own ring written
    the one way the code reads it ("none", "backbone", "disulfide"), and
    left off where it is empty or "as the run", so the run's choice applies.
    A value that names no ring is kept as written, so validate() can say
    which peptide it is on and what it says."""
    from . import peptide as pep
    out = []
    for p in peptides or []:
        q = dict(p) if isinstance(p, dict) else {"sequence": str(p)}
        own = q.get("cyclic")
        got = pep.ring_choice(own) if own is not None else ""
        if got == "":
            q.pop("cyclic", None)
        elif got:
            q["cyclic"] = got
        out.append(q)
    return out


def peptide_option_problems(j):
    """[str] why the job's cyclic or helical setting cannot be used, the
    run's or a peptide's own."""
    from . import peptide as pep
    out = []
    try:
        h = float(j.get("helical", 0))
        if not 0 <= h <= 100:
            raise ValueError
    except (TypeError, ValueError):
        out.append("The share of replicas that start as a helix must be a "
                   "number from 0 to 100, not %r." % j.get("helical"))
    mode = str(j.get("cyclic") or "none")
    if mode not in pep.CYCLIC:
        out.append("Cyclic must be %s, not %r." % (", ".join(pep.CYCLIC), mode))
        return out
    for i, p in enumerate(j.get("peptides") or [], 1):
        seq = pep.clean(p.get("sequence", ""))
        ring = ring_of(j, p)
        pr = (pep.cyclic_problem(seq, ring) if ring in pep.CYCLIC
              else pep.ring_problem(p.get("cyclic")))
        if pr:
            out.append("Peptide %d (%s): %s"
                       % (i, p.get("name") or seq[:10] or "?", pr.message))
    return out


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
        out.extend(adapt_problems(j))
        prot = j.get("protonation")
        if prot is not None and str(prot).strip().lower() not in ("reduce",
                                                                  "propka"):
            out.append("Protonation must be reduce or propka (PDB2PQR with "
                       "PROPKA at the pH given), not %r." % prot)
        out.extend(site_problems(j))
        out.extend(map_problems(j))
    from . import mac10
    engine, exe = resolve_engine(j)
    mac = engine == mac10.ENGINE
    if mac10.engaged() and not mac:
        out.append("This Mac docks with ADCP 1.0, because Scripps builds no "
                   "ADCP 1.1 for macOS, and this run asks for ADCP 1.1. "
                   "Dock it on Linux, a cluster or in Colab.")
    elif engine == "mac-only":
        out.append("This run was made on a Mac, to dock with ADCP 1.0; finish "
                   "it there, or dock its peptides in a new run here.")
    elif str(j.get("adcp_engine") or "1.1") not in ENGINES and not mac:
        out.append("The docking engine is ADCP 1.1: adcp_engine must be %s, "
                   "not %r." % (" or ".join(ENGINES), j.get("adcp_engine")))
    elif engine == "missing":
        out.append("ADCP 1.1 is not installed here. ./install.sh installs "
                   "it; INSTALL.md says how.")
    out.extend(peptide_option_problems(j))
    lib = None
    if exe and not mac:
        from . import adcp11 as _a
        lib = _a.library(exe)
    if mac and int(j.get("omm_minimize") or 0):
        out.append("Minimising poses with OpenMM needs ADCP 1.1, which "
                   "Scripps does not build for macOS.")
    if not j.get("peptides"):
        out.append("No peptides to dock.")
    for i, p in enumerate(j.get("peptides") or [], 1):
        seq, probs = pep.check(p.get("sequence", ""), library=lib,
                               max_length=j.get("max_length", pep.MAX_LENGTH),
                               extended=not mac)
        for pr in probs:
            out.append("Peptide %d (%s): %s"
                       % (i, p.get("name") or seq[:10] or "?", pr.message))
    # Duplicate sequences are fatal rather than untidy. A docking is
    # written to docking/<sequence>/ and its result to
    # results/<sequence>.json, so two entries with the same sequence
    # share a working directory and a result file.
    seqs, names = {}, {}
    peps = j.get("peptides") or []
    for i, p in enumerate(peps, 1):
        seq = pep.clean(p.get("sequence", ""))
        label = (p.get("name") or "").strip()
        if seq:
            if seq in seqs:
                first, second = (ring_of(j, peps[seqs[seq] - 1]),
                                 ring_of(j, p))
                if first != second:
                    # One sequence as two molecules, a chain and a ring, is
                    # a comparison worth making, but not inside one run.
                    out.append(
                        "Peptides %d and %d are the same sequence (%s), one "
                        "docked %s and the other %s. A run docks each "
                        "sequence once, since its working directory and "
                        "result file are named by the sequence; dock the "
                        "second in a run of its own to compare the two."
                        % (seqs[seq], i, seq,
                           pep.RING_AS.get(first, first),
                           pep.RING_AS.get(second, second)))
                else:
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
    most = max_poses(engine)
    if int(j.get("poses", 0)) < 1:
        out.append("Ask for at least one pose per peptide.")
    elif int(j.get("poses", 0)) > most:
        out.append("%s writes at most %d ranked poses per peptide; ask for "
                   "%d or fewer." % (mac10.LABEL if mac else "ADCP", most, most))
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
