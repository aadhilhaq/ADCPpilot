"""Estimating how long a run will take, and what to request from the queue.

By Aadhil Haq

Calibrated against ADCP 1.1.21 on a compute node: peptides of 6, 9 and 13
residues, one wave of 8 replicas on 8 cores at 250,000, 1,000,000 and
2,500,000 steps, and two waves of a 9-mer (the numbers are below, and in
tests/test_engine.py, which checks the estimate against each).

Everything rounds up. A job killed at the wall clock limit has produced
nothing and must start again, which costs more than the idle cores an
over-booking wastes.
"""
import math

# Wall seconds for one wave of replicas, at REFERENCE_LENGTH residues: a
# fixed part (starting ADCP's Python, unpacking the maps, clustering, building
# the poses' side chains) and a part per million steps. Measured with ADCP
# 1.1.21 on 8 cores, 8 replicas of a 6-mer (FSDLWK) into MDM2:
#
#     250,000 steps     41.1 s
#     1,000,000 steps   72.6 s
#     2,500,000 steps  132.1 s
#
# which is 31 s and 40.4 s per million steps.
REPLICA_SETUP_SECONDS = 31.0
REPLICA_SECONDS_PER_MSTEP = 41.0

# How much longer each residue past REFERENCE_LENGTH makes it. Measured the
# same way, at 2,500,000 steps: a 9-mer (RFMDYWEGL) 200.6 s, 1.5 times the
# 6-mer, and a 13-mer (ETFSDLWKLLPEN) 307.7 s, 2.3 times; 0.17 and 0.19 per
# residue, booked at 0.2.
PER_RESIDUE = 0.2

# The length at which the per-residue factor is 1.0.
#
# The factor applies at every length rather than only above ten residues,
# following ADCP's guidance that the difficulty scales throughout.
REFERENCE_LENGTH = 6

# Wall seconds for one pose's MM-GBSA against a receptor of MMGBSA_REF_ATOMS.
#
# tleap, one restrained sander minimisation and one MMPBSA.py call, in
# series; none is faster with more cores, only running several poses at once
# is. Cost scales with receptor size.
BOX_PER_RESIDUE = 1.14
BOX_INTERCEPT = 8.4
BOX_RESIDUAL_RMS = 2.7


def box_for_length(length):
    """(typical, safe) cube edge in A for a bound peptide of `length`."""
    n = max(1, int(length))
    typical = BOX_PER_RESIDUE * n + BOX_INTERCEPT
    return typical, typical + 2.0 * BOX_RESIDUAL_RMS


def box_warning(length, box_size, box_max=None):
    """None if the box holds that peptide, else a sentence saying so.

    `box_size` is the run's actual box, as one edge or three.
    """
    # box_size may be a numpy array, which has no truth value: `is None`
    # plus an explicit length is the check that works for a list, a
    # tuple, an array or a scalar.
    if not length or box_size is None:
        return None
    try:
        vals = [float(x) for x in box_size]
    except TypeError:
        vals = [float(box_size)]
    if not vals:
        return None
    edge = min(vals)
    typical, safe = box_for_length(length)
    if edge >= safe:
        return None
    # Below the typical line is the serious case: not "might not fit" but
    # "the median bound peptide of this length did not fit".
    worst = edge < typical
    msg = ("A %d-residue peptide needs about %.0f %s of box and this run has "
           "%.0f %s. ADCP searches inside the box and nowhere else, so what "
           "comes back will be folded to fit, and it will still score and "
           "still rank."
           % (int(length), typical if worst else safe,
              "\u00c5", edge, "\u00c5"))
    if not worst:
        msg = ("A %d-residue peptide needs %.0f to %.0f %s of box depending "
               "on whether it binds as a helix or extended, and this run has "
               "%.0f %s. An extended pose may not fit."
               % (int(length), typical, safe, "\u00c5", edge, "\u00c5"))
    if box_max is not None and float(box_max) < safe:
        msg += (" Raise 'Max box' to at least %.0f %s: the box is clamped "
                "to it." % (safe, "\u00c5"))
    return msg


MMGBSA_REF_SECONDS = 150.0
MMGBSA_REF_ATOMS = 895
MMGBSA_EXPONENT = 1.5


def mmgbsa_seconds(atoms=0):
    """Wall seconds for one pose against a receptor of `atoms` heavy atoms."""
    n = max(1, int(atoms or MMGBSA_REF_ATOMS))
    return MMGBSA_REF_SECONDS * (n / float(MMGBSA_REF_ATOMS)) ** MMGBSA_EXPONENT


# Poses are rescored concurrently (min(cores, poses, mm_workers) at a time)
# so wall time is the number of waves, not the number of poses.
MM_WORKERS_DEFAULT = 8

# The prepare stage: pdb2pqr, the restrained minimisation, prepare_receptor
# and agfr, as its own job with its own walltime. The floor grows with
# receptor size rather than being flat, and is set for the smallest receptor.
PREP_FLOOR_MINUTES = 4.0

# Per thousand heavy atoms, on top of the floor. Measured with the ADCP 1.1
# suite (reduce, AutoSite 1.1, no gradient), the whole stage took 1.25 min
# for MDM2 (1YCR chain A, 705 atoms, a 30 A box) and 4.3 min for the 7A3R
# dimer (5,912 atoms, one hetero group): 0.59 min per 1,000 atoms between
# them. This is about double that.
PREP_PER_1K_ATOMS = 1.0

# Per distinct non-water heteroatom residue type. Each must be typed by
# tleap, and each is an opportunity for preparation to need a second pass.
PREP_PER_HETATM_KIND = 2.0

# Rebuilding whole missing loops is a different amount of work from
# completing a side chain, so it is a multiplier rather than an addition.
PREP_FILL_LOOPS = 1.5

# The maps, measured on the MDM2 example's receptor (868 atoms) and its 30 A
# box, the ADCP 1.1 suite's agfr alone, one at a time on a compute node, on
# two days:
#
#     AutoSite 1.0, no receptor gradient (the defaults)     13.6 s   15.5 s
#     AutoSite 1.1 (-asv 1.1)                                16.3 s
#     the receptor gradient on (-g)                          36.9 s   68.4 s
#     -b fill (the AutoSite pocket)                           7.4 s
#
# The first agfr on a node whose file cache is cold takes up to a minute
# more, starting the suite's Python (71.5 s for AutoSite 1.1).
#
# agfr's grid, per 1,000 cubic Angstrom of box beyond a 20 A cube: about
# 0.009 min measured, booked at five times that.
GRID_PER_1K_BOX = 0.05
# The receptor gradient: 0.014 to 0.033 min per 1,000 measured, over the
# whole box, not beyond 20 A: it is paid on every grid point the receptor
# fills.
GRADIENT_PER_1K_BOX = 0.05
# AutoSite 1.1 scans its cutoffs where 1.0 runs once: 2.7 s more on MDM2.
AUTOSITE11_MINUTES = 0.25
AUTOSITE11_PER_1K_BOX = 0.01
# An AutoSite pocket: agfr -b fill searches the whole receptor at 1 A
# before the box exists, which grows with the receptor rather than the box.
# 7.4 s on MDM2; booked at 0.5 min plus 0.5 per 1,000 atoms.
POCKET_MINUTES = 0.5
POCKET_PER_1K_ATOMS = 0.5


# The protonation, on top of the per atom slope above, which was fitted to
# runs protonated by PDB2PQR with PROPKA. Measured on a compute node:
#
#                         MDM2 (706 atoms)   7A3R (5,912 atoms, 11 breaks)
#     PDB2PQR + PROPKA         1.2 s              4.1 s
#     reduce -BUILD            0.2 s              0.6 s
#
# (reduce 3.23, byte for byte the program in the ADCP 1.1 suite's bin folder)
#
# Booked at several times that, and reduce with a minute more: its search of
# a large network of Asn, Gln and His flips stops at 600 s a network.
PROTONATE_PER_1K_ATOMS = {"propka": 0.1, "reduce": 0.05}
PROTONATE_MINUTES = {"propka": 0.0, "reduce": 1.0}


def prep_minutes(atoms=0, hetatm_kinds=0, fill_loops=False, box=None,
                 gradient=True, autosite="1.0", pocket=False,
                 protonation="propka", prepared=True):
    """What to book for preparation, given what is known about the receptor.

    The defaults are what a job file written before these settings existed
    means: AutoSite 1.0 with the receptor gradient on (target.maps_for), and
    PDB2PQR with PROPKA. A box from an AutoSite pocket is not known until
    the stage runs, so `box` is then the largest it may be.

    `protonation` is "reduce" or "propka". `prepared` off is a receptor
    prepared elsewhere (in ADapt, or one that came prepared): it is copied
    in, so only the maps cost anything.
    """
    m = PREP_FLOOR_MINUTES
    if prepared:
        m += (max(0, int(atoms)) / 1000.0) * PREP_PER_1K_ATOMS
        m += max(0, int(hetatm_kinds)) * PREP_PER_HETATM_KIND
        how = protonation if protonation in PROTONATE_PER_1K_ATOMS else "propka"
        m += PROTONATE_MINUTES[how]
        m += (max(0, int(atoms)) / 1000.0) * PROTONATE_PER_1K_ATOMS[how]
    vol = 0.0
    if box:
        try:
            vol = float(box[0]) * float(box[1]) * float(box[2])
        except (TypeError, ValueError, IndexError):
            vol = 0.0
    m += max(0.0, (vol - 8000.0) / 1000.0) * GRID_PER_1K_BOX
    if gradient:
        m += (vol / 1000.0) * GRADIENT_PER_1K_BOX
    if str(autosite) == "1.1":
        m += AUTOSITE11_MINUTES + (vol / 1000.0) * AUTOSITE11_PER_1K_BOX
    if pocket:
        m += POCKET_MINUTES + (max(0, int(atoms)) / 1000.0) * POCKET_PER_1K_ATOMS
    if fill_loops:
        m *= PREP_FILL_LOOPS
    return max(PREP_FLOOR_MINUTES, m)


# Docking options that change what one peptide costs, measured with 8
# replicas of 300,000 steps on 4 cores (ADCP 1.1.21, MDM2):
#
#     linear                          77.2 s of search
#     backbone cyclic (-cyc)         122.3 s   x 1.58
#     disulfide (-cys), 2 cysteines   88.3 s   x 1.14 (another sequence)
#     -p 25                           77.2 s   x 1.00
#     every replica helical           82.3 s   x 1.07
#
# Booked with room on top of each.
CYCLIC_FACTOR = {"none": 1.0, "backbone": 1.8, "disulfide": 1.3}
# A share of helical starts is the same search from other starting shapes.
HELICAL_FACTOR = 1.1
# More poses: ADCP writes each one out with its side chains, and the run
# splits them into files. 90 more poses took 3.3 s more on MDM2.
POSE_SECONDS = 0.2
# Keeping the trajectories: packing each replica's trajectory and logs into
# one archive. The whole working folder of a 50 replica, 2,500,000 step run
# (194 files, 231 MB with the unpacked maps) packed in 10.5 to 11.4 s; the
# archive leaves the maps out.
PACK_SECONDS_PER_REPLICA = 0.4


def dock_option_minutes(replicas, poses=10, keep=False):
    """What the docking options add to one peptide's task, beyond the
    search itself: writing more poses, and packing the trajectories."""
    m = max(0, int(poses) - 10) * POSE_SECONDS / 60.0
    if keep:
        m += max(1, int(replicas)) * PACK_SECONDS_PER_REPLICA / 60.0
    return m


#: Kept for anything that still names them.
PREP_MINUTES = PREP_FLOOR_MINUTES
MMGBSA_SECONDS = MMGBSA_REF_SECONDS
#: on top of twice the estimate
BUFFER_MINUTES = 30.0
CORES_DEFAULT = 48


def dock_minutes(length, replicas, steps, cores, option_factor=1.0):
    """Wall minutes for one peptide's search.

    Waves, not core-minutes: `cores` replicas run side by side, so what
    matters is how many rounds are needed. Below the core count, adding
    replicas costs no additional wall time. `option_factor` is what a
    cyclic peptide or helical starts cost on top (search_factor).
    """
    cores = max(1, int(cores or 1))
    replicas = max(1, int(replicas))
    factor = 1.0 + PER_RESIDUE * max(0, int(length) - REFERENCE_LENGTH)
    factor *= max(1.0, float(option_factor or 1.0))
    per_replica = (REPLICA_SETUP_SECONDS
                   + REPLICA_SECONDS_PER_MSTEP * (float(steps) / 1e6)) * factor
    waves = -(-replicas // cores)                 # ceil
    # The clustering at the end is serial and does not scale with
    # anything measured here.
    return waves * per_replica / 60.0 * 1.08 + 0.2


def mmgbsa_minutes(poses, cores=0, workers=MM_WORKERS_DEFAULT, atoms=0):
    """Wall minutes to rescore `poses`, several at a time.

    They run concurrently, as many as there are workers, so the wall time is
    the number of waves rather than the number of poses.
    """
    poses = max(0, int(poses))
    if not poses:
        return 0.0
    n = max(1, min(int(cores or workers), poses, int(workers)))
    waves = -(-poses // n)                       # ceil
    return waves * mmgbsa_seconds(atoms) / 60.0


def per_task(length, replicas, steps, cores, mmgbsa=False, mm_poses=0,
             atoms=0):
    """The longest one array task will run, in minutes."""
    t = dock_minutes(length, replicas, steps, cores)
    if mmgbsa:
        t += mmgbsa_minutes(mm_poses, cores, atoms=atoms)
    return t


# Every booking is twice the estimate and half an hour more (BUFFER_MINUTES):
# a job killed at the wall has produced nothing and has already spent its
# queue time, while a cluster charges for the time a job uses.
SAFETY = 2.0

# No job is booked for less than two hours. The model accounts for the
# search and the rescoring, but not for a cold filesystem, a slow module
# load, or a loaded node.
FLOOR_MINUTES = 120.0

# The prepare job is never booked for less than four hours, what every
# prepare job was given before it was priced: AGFR's run time grows with the
# box and the receptor in ways the model only roughly follows.
PREP_FLOOR_MINUTES_BOOKED = 240.0


def book_prep(minutes):
    """The prepare job's booking in minutes: as _book, with a four hour floor."""
    return max(PREP_FLOOR_MINUTES_BOOKED, _book(minutes))

# Gathering the results: reading each peptide's JSON, ranking, writing three
# files. Seconds in practice, but it is a queued job.
REPORT_MINUTES = 2.0
# What the report job books, whatever the estimate: the submission asks for
# this, so it is what the page adds up when it states the total booked.
REPORT_WALLTIME = "01:00:00"
REPORT_BOOKED_MINUTES = 60.0


def search_factor(cyclic="none", helical=0):
    """How much longer one replica's search runs for these options."""
    f = CYCLIC_FACTOR.get(str(cyclic or "none"), max(CYCLIC_FACTOR.values()))
    try:
        if float(helical or 0) > 0:
            f *= HELICAL_FACTOR
    except (TypeError, ValueError):
        pass
    return f


def plan(peptides, replicas, steps, cores=None, mmgbsa=False, mm_poses=5,
         throttle=0, atoms=0, hetatm_kinds=0, fill_loops=False, box=None,
         options=None):
    """Every stage, added up, end to end.

    The estimate is cumulative: a run is three queued jobs in a chain, and
    what is wanted is when there will be an answer, not how long one array
    task takes.

    `options` holds the run's settings that change the cost: autosite,
    receptor_gradient, site_mode, cyclic, helical, poses,
    keep_trajectories, box_max, protonation, from_adapt, receptor_prepared.
    Missing ones cost what they did before the settings existed. "rings",
    when given, is each peptide's own ring in the order of `peptides`; a
    peptide without one is priced with the run's "cyclic".
    """
    o = dict(options or {})
    lengths = [len(p) if isinstance(p, str) else int(p) for p in peptides]
    given = list(o.get("rings") or [])
    pairs = [(x, given[i] if i < len(given) and given[i] else o.get("cyclic"))
             for i, x in enumerate(lengths) if x > 0] or [(14, o.get("cyclic"))]
    lengths = [x for x, _r in pairs]
    cores = int(cores or CORES_DEFAULT)
    # More cores than replicas is cores that cannot be given work.
    cores = max(1, min(cores, max(4, int(replicas))))
    longest = max(lengths)
    n = len(lengths)

    pocket = o.get("site_mode") == "pocket"
    if pocket:
        # The pocket's box is not known until AutoSite has run; the maps
        # are booked for the largest it may be.
        edge = float(o.get("box_max") or 30.0)
        prep_box = (edge, edge, edge)
    else:
        prep_box = box
    prep = prep_minutes(atoms, hetatm_kinds, fill_loops, prep_box,
                        gradient=bool(o.get("receptor_gradient", True)),
                        autosite=o.get("autosite", "1.0"), pocket=pocket,
                        protonation=o.get("protonation", "propka"),
                        prepared=not (o.get("from_adapt")
                                      or o.get("receptor_prepared")))
    # Each peptide is one array task, priced with its own length and ring.
    # Every task books the same walltime, so the booking is the costliest
    # task's: the longest peptide, unless a shorter one's ring makes it the
    # slower search.
    helical = o.get("helical", 0)
    costs = [(dock_minutes(x, replicas, steps, cores, search_factor(r, helical)),
              x, r) for x, r in pairs]
    worst = max(costs, key=lambda c: (c[0], c[1]))
    factor = search_factor(worst[2], helical)
    dock = worst[0]
    extra = dock_option_minutes(replicas, o.get("poses", 10),
                                bool(o.get("keep_trajectories")))
    dock += extra
    mm = mmgbsa_minutes(mm_poses, cores, atoms=atoms) if mmgbsa else 0.0
    task = dock + mm
    at_once = int(throttle) if int(throttle or 0) > 0 else n
    waves = max(1, -(-n // max(1, at_once)))
    array = task * waves
    total = prep + array + REPORT_MINUTES
    rings_differ = len({str(r or "none") for _x, r in pairs}) > 1
    return {
        "n": n, "longest": longest, "cores": cores, "waves": waves,
        "prep_minutes": round(prep, 1),
        "dock_minutes": round(dock, 1),
        "mmgbsa_minutes": round(mm, 1),
        "task_minutes": round(task, 1),
        "array_minutes": round(array, 1),
        "report_minutes": round(REPORT_MINUTES, 1),
        "total_minutes": round(total, 1),
        "options_minutes": round(extra, 2),
        "search_factor": round(factor, 2),
        # The task the booking was sized for, and whether the peptides are
        # docked with more than one kind of ring, which is when that task
        # need not be the longest peptide's.
        "sized_for": {"length": worst[1], "ring": str(worst[2] or "none")},
        "rings_differ": rings_differ,
        "core_hours": round(
            sum(c + extra
                + (mmgbsa_minutes(mm_poses, cores, atoms=atoms)
                   if mmgbsa else 0.0)
                for c, _x, _r in costs) * cores / 60.0, 1),
        # A run whose box cannot hold its own longest peptide costs
        # as much as one that can, and returns a folded pose that
        # still scores.
        "box_needed": [round(v, 1) for v in box_for_length(longest)],
        "box_warning": box_warning(longest, box),
    }


def _book(minutes):
    """A booking for one job: twice the estimate and the buffer, never under
    the floor."""
    return max(FLOOR_MINUTES, SAFETY * float(minutes) + BUFFER_MINUTES)


def book(minutes):
    """The walltime to ask for a job estimated at `minutes`, as hh:mm:ss."""
    return hhmmss(_book(minutes))


def suggest(peptides, replicas, steps, cores=None, mmgbsa=False, mm_poses=5,
            throttle=0, buffer_minutes=None, atoms=0, hetatm_kinds=0,
            fill_loops=False, box=None, options=None, prep_booked=None):
    """What to put in the form, and what to report to the user.

    `buffer_minutes` is accepted and ignored: the buffer is BUFFER_MINUTES.
    `prep_booked` is what the prepare job will book, in minutes, when the
    caller has worked it out the way the submission does (pipeline.prep_time
    prices the maps for the largest box the run allows); without it, the
    estimate's own booking is used.

    The total booked is the sum of what the three jobs ask the queue for,
    end to end: the prepare job, each wave of docking tasks, and the report.
    It is not the chain's estimate booked as if it were one job, which is a
    number no job asks for.
    """
    p = plan(peptides, replicas, steps, cores, mmgbsa, mm_poses, throttle,
             atoms, hetatm_kinds, fill_loops, box, options)
    asked = _book(p["task_minutes"])
    prep = (float(prep_booked) if prep_booked is not None
            else book_prep(p["prep_minutes"]))
    # Each docking task books the walltime the form is given, whole minutes,
    # and a throttled array runs its tasks in waves, one after another.
    end_to_end = (prep + math.ceil(asked) * p["waves"]
                  + REPORT_BOOKED_MINUTES)
    # Written as the jobs that appear in the queue rather than as
    # internal stages: the docking and the rescoring are one SLURM task,
    # so quoting them separately invites comparison with a queue row that
    # covers both.
    if p["mmgbsa_minutes"]:
        task = ("one dock+mm task of %s (docking %s, then MM-GBSA %s; "
                "Slurm runs them as one job)"
                % (human(p["dock_minutes"] + p["mmgbsa_minutes"]),
                   human(p["dock_minutes"]), human(p["mmgbsa_minutes"])))
    else:
        task = "one dock task of %s" % human(p["dock_minutes"])
    stages = ["A prepare job of %s" % human(p["prep_minutes"]), task,
              "a report job of %s" % human(p["report_minutes"])]
    if p.get("rings_differ"):
        # The peptides are docked with different rings, and a ring makes
        # the search slower, so the task sized for is named by its ring too.
        from . import peptide as _pep
        ring = p["sized_for"]["ring"]
        whom = ("for the peptide that takes longest to dock (%d residues, %s)"
                % (p["sized_for"]["length"],
                   _pep.RING_NAMES.get(ring, ring)))
    else:
        whom = "for the longest peptide (%d residues)" % p["longest"]
    out = dict(p)
    out.update({
        "estimate_minutes": p["task_minutes"],
        "walltime_minutes": int(math.ceil(asked)),
        "walltime": hhmmss(asked),
        "prep_walltime": hhmmss(prep),
        "end_to_end": hhmmss(end_to_end),
        # What each job books, for the page to say what the total is made of.
        "booked": {"prepare": hhmmss(prep), "dock": hhmmss(asked),
                   "waves": p["waves"],
                   "report": hhmmss(REPORT_BOOKED_MINUTES),
                   "total": hhmmss(end_to_end)},
        "total_core_hours": p["core_hours"],
        "why": (
            "%s, %s. "
            "%d peptide%s in %s, so about %s end to end; "
            "each job is booked for twice its own estimate and %s more, and "
            "never less than %s (the prepare job %s), because a job "
            "killed at the wall has produced nothing and has already spent "
            "its queue time."
            % (", then ".join(stages), whom, p["n"],
               "" if p["n"] == 1 else "s",
               "one wave" if p["waves"] == 1
               else "%d waves" % p["waves"],
               human(p["total_minutes"]), human(BUFFER_MINUTES),
               human(FLOOR_MINUTES), human(PREP_FLOOR_MINUTES_BOOKED))),
    })
    return out


def receptor_cost(path, chains=()):
    """(heavy atoms, distinct non-water heteroatom kinds) in a structure.

    Counted over the chains that will be kept rather than the whole file: a
    structure is routinely several copies of a complex, and costing a run by
    the file means costing work it will never do.
    """
    want = {c.strip() for c in (chains or []) if c and c.strip()}
    atoms, kinds = 0, set()
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                if line.startswith("ENDMDL"):
                    break
                if line[:4] not in ("ATOM", "HETA"):
                    continue
                if line[16] not in (" ", "A"):     # one altloc only
                    continue
                chain = line[21].strip() or "_"
                if want and chain not in want:
                    continue
                res = line[17:20].strip()
                if line.startswith("HETATM"):
                    if res in ("HOH", "WAT", "DOD", "D2O"):
                        continue
                    kinds.add(res)
                if (line[76:78].strip() or res[:1]).upper() == "H":
                    continue
                atoms += 1
    except (OSError, TypeError):
        return 0, 0
    return atoms, len(kinds)


def hhmmss(minutes):
    m = max(1, int(math.ceil(float(minutes))))
    return "%02d:%02d:00" % (m // 60, m % 60)


def parse_hhmmss(text):
    """'08:00:00' or '8:30' or '90' (minutes) -> minutes."""
    t = str(text or "").strip()
    if not t:
        return 0
    if ":" not in t:
        try:
            return int(float(t))
        except ValueError:
            return 0
    bits = [int(x or 0) for x in t.split(":")[:3]]
    while len(bits) < 3:
        bits.append(0)
    return bits[0] * 60 + bits[1] + (1 if bits[2] else 0)


def human(minutes):
    """A duration as a user would say it. "2 h 00 min" is not."""
    m = int(round(float(minutes)))
    if m < 60:
        return "%d min" % m
    h, r = m // 60, m % 60
    if not r:
        return "%d hour" % h if h == 1 else "%d hours" % h
    return "%d h %d min" % (h, r)
