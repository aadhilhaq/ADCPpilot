"""Estimating how long a run will take, and what to request from the queue.

By Aadhil Haq

Calibrated against measured runs: an 11-residue peptide at 100 replicas by
5 million steps took about eight minutes on 48 cores, and a 15-residue one
at the same settings about twelve.

Everything rounds up. A job killed at the wall clock limit has produced
nothing and must start again, which costs more than the idle cores an
over-booking wastes.
"""
import math

# Core-minutes per (replica x million steps) of ADCP search, at ten residues.
#
# Two measurements disagree by about a factor of two, so the larger is used.
REPLICA_SETUP_SECONDS = 4.8
REPLICA_SECONDS_PER_MSTEP = 34.0

PER_RESIDUE = 0.06

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
        msg = ("A %d-residue peptide needs %.0f-%.0f %s of box depending on "
               "whether it binds as a helix or extended, and this run has "
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

# Per thousand heavy atoms, on top of the floor. The measured slope is about
# 0.5 min per 1,000 atoms; this is double that.
PREP_PER_1K_ATOMS = 1.0

# Per distinct non-water heteroatom residue type. Each must be typed by
# tleap, and each is an opportunity for preparation to need a second pass.
PREP_PER_HETATM_KIND = 2.0

# Rebuilding whole missing loops is a different amount of work from
# completing a side chain, so it is a multiplier rather than an addition.
PREP_FILL_LOOPS = 1.5

# agfr's grid, per 1,000 cubic Angstrom of box beyond a 20 A cube.
PREP_PER_1K_BOX = 0.15


def prep_minutes(atoms=0, hetatm_kinds=0, fill_loops=False, box=None):
    """What to book for preparation, given what is known about the receptor."""
    m = PREP_FLOOR_MINUTES
    m += (max(0, int(atoms)) / 1000.0) * PREP_PER_1K_ATOMS
    m += max(0, int(hetatm_kinds)) * PREP_PER_HETATM_KIND
    if box:
        try:
            vol = float(box[0]) * float(box[1]) * float(box[2])
        except (TypeError, ValueError, IndexError):
            vol = 0.0
        m += max(0.0, (vol - 8000.0) / 1000.0) * PREP_PER_1K_BOX
    if fill_loops:
        m *= PREP_FILL_LOOPS
    return max(PREP_FLOOR_MINUTES, m)


#: Kept for anything that still names them.
PREP_MINUTES = PREP_FLOOR_MINUTES
MMGBSA_SECONDS = MMGBSA_REF_SECONDS
#: on top of twice the estimate
BUFFER_MINUTES = 30.0
CORES_DEFAULT = 48


def dock_minutes(length, replicas, steps, cores):
    """Wall minutes for one peptide's search.

    Waves, not core-minutes: `cores` replicas run side by side, so what
    matters is how many rounds are needed. Below the core count, adding
    replicas costs no additional wall time.
    """
    cores = max(1, int(cores or 1))
    replicas = max(1, int(replicas))
    factor = 1.0 + PER_RESIDUE * max(0, int(length) - REFERENCE_LENGTH)
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


# The safety factor applied to every booking.
SAFETY = 1.5

# No job is booked for less than two hours, with the safety factor applied on
# top for anything larger. The model accounts for the search and the
# rescoring, but not for a cold filesystem, a slow module load, or a loaded
# login node.
FLOOR_MINUTES = 120.0

# Gathering the results: reading each peptide's JSON, ranking, writing three
# files. Seconds in practice, but it is a queued job.
REPORT_MINUTES = 2.0


def plan(peptides, replicas, steps, cores=None, mmgbsa=False, mm_poses=5,
         throttle=0, atoms=0, hetatm_kinds=0, fill_loops=False, box=None):
    """Every stage, added up, end to end.

    The estimate is cumulative: a run is three queued jobs in a chain, and
    what is wanted is when there will be an answer, not how long one array
    task takes.
    """
    lengths = [len(p) if isinstance(p, str) else int(p) for p in peptides]
    lengths = [n for n in lengths if n > 0] or [14]
    cores = int(cores or CORES_DEFAULT)
    # More cores than replicas is cores that cannot be given work.
    cores = max(1, min(cores, max(4, int(replicas))))
    longest = max(lengths)
    n = len(lengths)

    prep = prep_minutes(atoms, hetatm_kinds, fill_loops, box)
    dock = dock_minutes(longest, replicas, steps, cores)
    mm = mmgbsa_minutes(mm_poses, cores, atoms=atoms) if mmgbsa else 0.0
    task = dock + mm
    at_once = int(throttle) if int(throttle or 0) > 0 else n
    waves = max(1, -(-n // max(1, at_once)))
    array = task * waves
    total = prep + array + REPORT_MINUTES
    return {
        "n": n, "longest": longest, "cores": cores, "waves": waves,
        "prep_minutes": round(prep, 1),
        "dock_minutes": round(dock, 1),
        "mmgbsa_minutes": round(mm, 1),
        "task_minutes": round(task, 1),
        "array_minutes": round(array, 1),
        "report_minutes": round(REPORT_MINUTES, 1),
        "total_minutes": round(total, 1),
        "core_hours": round(
            sum(dock_minutes(x, replicas, steps, cores)
                + (mmgbsa_minutes(mm_poses, cores, atoms=atoms)
                   if mmgbsa else 0.0)
                for x in lengths) * cores / 60.0, 1),
        # A run whose box cannot hold its own longest peptide costs
        # as much as one that can, and returns a folded pose that
        # still scores.
        "box_needed": [round(v, 1) for v in box_for_length(longest)],
        "box_warning": box_warning(longest, box),
    }


def _book(minutes):
    """A booking for one job: the safety factor, never under the floor."""
    return max(FLOOR_MINUTES, SAFETY * float(minutes))


def suggest(peptides, replicas, steps, cores=None, mmgbsa=False, mm_poses=5,
            throttle=0, buffer_minutes=None, atoms=0, hetatm_kinds=0,
            fill_loops=False, box=None):
    """What to put in the form, and what to report to the user.

    `buffer_minutes` is accepted and ignored: the safety margin is a
    multiplier, and the floor does what a fixed addition used to.
    """
    p = plan(peptides, replicas, steps, cores, mmgbsa, mm_poses, throttle,
             atoms, hetatm_kinds, fill_loops, box)
    asked = _book(p["task_minutes"])
    prep = _book(p["prep_minutes"])
    end_to_end = _book(p["total_minutes"])
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
    stages = ["a prepare job of %s" % human(p["prep_minutes"]), task,
              "a report job of %s" % human(p["report_minutes"])]
    out = dict(p)
    out.update({
        "estimate_minutes": p["task_minutes"],
        "walltime_minutes": int(math.ceil(asked)),
        "walltime": hhmmss(asked),
        "prep_walltime": hhmmss(prep),
        "end_to_end": hhmmss(end_to_end),
        "total_core_hours": p["core_hours"],
        "why": (
            "%s, for the longest peptide (%d residues). "
            "%d peptide%s in %s, so about %s end to end; "
            "no job is booked for less than %s, and anything bigger gets "
            "%.1fx its own stage, because a job killed at the wall has produced "
            "nothing and has already spent its queue time."
            % (", then ".join(stages), p["longest"], p["n"],
               "" if p["n"] == 1 else "s",
               "one wave" if p["waves"] == 1
               else "%d waves" % p["waves"],
               human(p["total_minutes"]), human(FLOOR_MINUTES), SAFETY)),
    })
    return out


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
