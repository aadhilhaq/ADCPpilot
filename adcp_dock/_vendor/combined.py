"""The combined ranking: the measurements the panel already makes, folded
into one order instead of sitting in the table being ignored.

By: Aadhil Haq

``rank.joint`` ranks on the docking score, the site contacts and MM-GBSA.
``panel`` then measures buried area, hydrogen bonds, salt bridges, stacking,
clashes, site coverage and the per residue decomposition, and says in its own
docstring that none of it changes the ranking. This module is where it does.

Nothing here replaces ``rank.joint``. It writes a second order beside the
first, for the same reason HOPE-PE does: a number that has been reported must
not change underneath the person who reported it.

Why these terms and not the obvious ones
----------------------------------------
Measured on the authors' own runs before anything was designed: three
receptors, 6B88, CLas SPase1 and a DENV envelope dimer with F-loop 12-mers;
21 constructs and 105 scored poses in all. Mean
Spearman rho within a run, so the figures are not reporting the difference
between the receptors:

                    vs docking score   vs MM-GBSA   vs length
    bsa                   -0.65          -0.51        +0.64
    n_contact             -0.63          -0.38        +0.36
    bsa_per_res           +0.05          -0.15        -0.54
    n_hbond               -0.37          -0.38        +0.02
    n_hbond + n_ionic     -0.44          -0.45        +0.04
    site_dg               -0.06          +0.27        +0.03

Buried area, contact count and the docking score are one quantity, not three.
They agree with each other at rho 0.83 (bsa against n_contact) and each tracks
peptide length. Adding all three z-scores would rank the library by size three
times over and call the result a consensus. So the size axis is counted once,
by the energy, and the interface terms are the ones that are left when size is
taken out of them:

  bsa_per_res   buried area per peptide residue. Rho +0.05 against the
                docking score: dividing by residue count removes exactly the
                bulk the score already carries, and what is left is how well
                the peptide fits rather than how much of it there is. This is
                not the ligand efficiency ``rank`` refuses - that divides the
                ENERGY by atom count and so ranks glycine-rich peptides whose
                contacts are backbone. Buried area per residue does not: on
                6B88 the glycine-rich WHAGGAGGGAGG has the lowest value in the
                run at 82 square angstroms per residue and MMNFGWH the
                highest at 120.

  n_hbond,      counted raw, NOT per residue. The instinct is to normalize
  n_ionic       every count the way bsa is normalized, and the measurements
                say not to: a hydrogen bond count is already independent of
                length (rho +0.02), so dividing by length would not remove a
                size bias, it would introduce one in the other direction
                (hydrogen bonds per residue: rho -0.39 against length). The
                two are summed into one term because they agree with each
                other at rho +0.42 and are the same kind of evidence.

  site_dg       the site residues' share of the per residue decomposition -
                "strong per-residue contributions", which is what the
                decomposition was computed for. The most independent term in
                the set: rho -0.06 against the docking score and +0.03 against
                length. Its +0.27 against MM-GBSA is the DISAGREEING sign for
                two energies, which is the point of having it - a construct can
                reach a good total energy away from the residues the run was
                aimed at, and total MM-GBSA cannot tell you that it did.

  site_coverage how much of the site is touched at all. Independent of
                everything (rho within 0.11 of zero), but 1.0 for every
                construct in two of the three runs measured, so on most runs
                it is constant and drops out of its own accord - a z-score of
                a constant column is zeros. It is kept because a run where it
                is not constant is a run where it matters.

Clashes are not a term. See CLASH_PENALTY.

Weights are trust, not an exchange rate
---------------------------------------
No weight here is fitted, and none should be presented as validated. A
designed construct that has not been tested has no affinity to fit against;
and the benchmark set of solved complexes that is available - peptides whose
bound answer is on the deposited coordinates - runs to eleven, nowhere near
enough to fit seven weights without fitting the noise instead.

So the terms are converted to a position within the run before they are
weighted, exactly as ``panel.AGREE`` argues, and a weight says HOW MUCH THAT
MEASURE IS TRUSTED and nothing else. The relative order follows the
literature ``panel`` already cites: MM-GBSA the best of the cheap predictors
of affinity but modest in absolute terms, buried area a real correlate of
binding across interfaces, an interface hydrogen bond worth roughly 0.5 to 1.5
kcal/mol, salt bridges similar and highly dependent on burial.

Two questions follow, and they are different. Whether the weights matter at
all is answered here: ``stability`` sweeps them and reports how much the
order moves. Whether this ordering beats the single column it sits beside is
NOT answered here and cannot be - that needs those solved complexes, asking
whether the score puts a known binder above decoys and how far down the list
you must go to find it. Until somebody runs that, this metric is better
argued for than the single column, not better shown.
"""
from __future__ import annotations

import numpy as np

from . import rank

# A term's value comes from the first of these keys the row carries, so one
# code path serves the pipelines' final_constructs.tsv and the docking
# console's results.tsv without either having to rename its columns.
#
# direction is which way is better and is the one thing in this module that
# silently inverts the whole ranking if it is wrong: "low" for an energy or a
# score, "high" for a count, an area or a fraction. There is a test for every
# row of this table.
#
# (name, keys, direction, default weight, what it is)
TERMS = [
    ("mmgbsa", ("mmgbsa", "mmgbsa_top5"), "low", None,
     "MM-GBSA interaction energy"),
    ("score", ("affinity", "adcp_top5", "score"), "low", None,
     "docking score"),
    ("fit", ("bsa_per_res",), "high", 0.6,
     "buried area per peptide residue"),
    ("specific", ("n_hbond",), "high", 0.4,
     "hydrogen bonds and salt bridges"),
    ("site", ("site_dg",), "low", 0.5,
     "the site residues' share of the per residue decomposition"),
    ("coverage", ("site_coverage",), "high", 0.3,
     "how much of the site is touched"),
    ("cluster", ("cluster_frac",), "high", 0.25,
     "share of the search's solutions in the top pose's cluster"),
]

# Salt bridges ride along with the hydrogen bonds rather than carrying a term.
# They are the same kind of evidence, they agree with hydrogen bonds at rho
# +0.42, and they are far rarer - 0 to 3 per pose across every run measured -
# so a term of their own would spend a whole axis on a count that is usually
# zero or one.
SPECIFIC_ALSO = ("n_ionic",)

# The energy axis, and why it is one axis rather than two.
#
# The docking score and MM-GBSA agree at rho +0.66 across the three runs: they
# are two estimates of the same quantity, one cheap and one considered. Giving
# each a term of its own would weight that quantity twice, which is the fault
# this module exists to avoid, so when a run has both they are averaged into
# one term - both visible, both contributing, the axis counted once. A run
# with only one of them uses that one, on the same axis and at the same
# weight, because the axis is "how strongly is it predicted to bind" and not
# "which program said so".
ENERGY_WEIGHT = 1.0

# Clashes do not raise anything, and they are not a gate either.
#
# A gate would remove a construct from the ranking for a measurement that is
# usually a geometry artifact of a rescoring minimization rather than evidence
# the design is bad, and removing it is the one thing that cannot be undone by
# the person reading the table. So it is a penalty, taken from the combined
# score and never added to it, and the count stays in the table beside it.
#
# It has not been exercised. n_clash is 0 for every pose of every construct in
# all three runs measured, so this penalty has never changed an order on real
# data and is tested only against constructed rows. Said here rather than
# discovered later.
CLASH_PENALTY = 0.5

# Below this many constructs a z-score is not worth computing: with three rows
# every term is -1, 0, +1 whatever the values were, and the combined score
# reports the arithmetic rather than the measurements.
MIN_ROWS = 4


def _num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x else x


def _value(row, keys):
    for k in keys:
        if k in row:
            v = _num(row[k])
            if v is not None:
                return v
    return None


def weights(**over):
    """The default weights, with any of them overridden by name."""
    w = {"energy": ENERGY_WEIGHT, "clash": CLASH_PENALTY}
    for name, _keys, _dir, default, _what in TERMS:
        if default is not None:
            w[name] = default
    bad = set(over) - set(w)
    if bad:
        raise KeyError("no such weight: %s" % ", ".join(sorted(bad)))
    w.update({k: float(v) for k, v in over.items()})
    return w


def _column(rows, keys, direction):
    """One term's values, or None when the term is off for this run.

    All-or-nothing, which ``rank.joint`` establishes and is right about: a
    z-score over a subset puts the measured and the unmeasured rows on
    different scales, so one row without a value turns the term off for
    everybody rather than being filled with a zero. A fabricated zero would
    make "not measured" and "measured as touching nothing" the same row.
    """
    vals = [_value(r, keys) for r in rows]
    if any(v is None for v in vals):
        return None
    v = np.array(vals, float)
    # The sign convention, in one place. Everything below this line is
    # "higher is better", so an energy is negated here and nowhere else.
    return -v if direction == "low" else v


def score(rows, w=None, min_rows=MIN_ROWS):
    """Add the combined score and its components to each row, best first.

    Returns ``(rows, info)``. Each row gains:

        combined_z        the score, higher being better
        combined_terms    {term: its contribution}, so the person can see WHY
                          a construct ranked where it did rather than being
                          handed one opaque number
        combined_off      the terms that were not available for every row and
                          so were switched off for the whole run

    Rows that could not be scored at all keep ``combined_z`` of None and sort
    after the ones that could. They have not been shown to be worse; they have
    not been measured.
    """
    w = w or weights()
    info = {"n": len(rows), "used": {}, "off": [], "constant": [],
            "weights": dict(w), "model": {}}
    if not rows:
        return [], info
    if len(rows) < min_rows:
        # Not a failure and not an error: the run is too small for a position
        # within it to mean anything, so the combined score is not computed
        # and the caller keeps whatever order it had.
        info["off"] = [t[0] for t in TERMS]
        info["too_few"] = True
        for r in rows:
            r["combined_z"] = None
            r["combined_terms"] = {}
            r["combined_off"] = "run too small for a combined score"
        return list(rows), info

    total = np.zeros(len(rows), float)
    parts = {}
    # Where this run's z-scores came from, kept so one construct's other poses
    # can be placed on the same scale afterwards. Re-deriving the scale from a
    # different set of poses would move every construct at once and report the
    # whole run's variation as if it were one construct's.
    model = info["model"]

    def fitted(name, col, keys):
        model[name] = {"mean": float(col.mean()), "sd": float(col.std()),
                       "keys": list(keys)}
        return rank.zscore(col)

    # The energy axis: the mean of whichever of the two estimates the run has.
    energy = []
    for name in ("mmgbsa", "score"):
        keys = next(t[1] for t in TERMS if t[0] == name)
        col = _column(rows, keys, "low")
        if col is None:
            info["off"].append(name)
            continue
        z = fitted(name, col, keys)
        if not z.any():
            info["constant"].append(name)
            model.pop(name, None)
            continue
        energy.append((name, z))
    if energy:
        z = sum(z for _n, z in energy) / len(energy)
        parts["energy"] = w["energy"] * z
        total = total + parts["energy"]
        info["used"]["energy"] = [n for n, _z in energy]

    for name, keys, direction, default, _what in TERMS:
        if name in ("mmgbsa", "score"):
            continue
        col = _column(rows, keys, direction)
        if col is None:
            info["off"].append(name)
            continue
        made_of = list(keys)
        if name == "specific":
            # Salt bridges join the hydrogen bonds only when every row has
            # them too; the same all-or-nothing rule, applied to the second
            # half of a term rather than to a term.
            extra = _column(rows, SPECIFIC_ALSO, "high")
            if extra is not None:
                col = col + extra
                made_of = ["n_hbond"] + list(SPECIFIC_ALSO)
            else:
                made_of = ["n_hbond"]
        z = fitted(name, col, keys)
        if not z.any():
            # A column with no spread - site coverage of 1.0 for everybody,
            # which is the usual case - contributes nothing, and saying so is
            # better than letting a weight sit in the table doing nothing.
            info["constant"].append(name)
            model.pop(name, None)
            continue
        parts[name] = w.get(name, 0.0) * z
        total = total + parts[name]
        info["used"][name] = made_of

    # The penalty, last, so it can only ever take away.
    clash = _column(rows, ("n_clash",), "high")
    if clash is not None and clash.any():
        pen = w["clash"] * rank.zscore(clash)
        # Only the rows that actually clash are touched. A z-score would hand
        # a bonus to every row below the mean, and a construct with no clashes
        # has earned nothing by not having any - it has merely not been
        # penalized.
        pen = np.where(clash > 0, np.maximum(pen, 0.0), 0.0)
        parts["clash"] = -pen
        total = total - pen
        info["used"]["clash"] = ["n_clash"]

    scored = bool(parts)
    for i, r in enumerate(rows):
        r["combined_terms"] = {k: float(v[i]) for k, v in parts.items()}
        r["combined_z"] = float(total[i]) if scored else None
        r["combined_off"] = ", ".join(info["off"])
    out = sorted(rows, key=lambda r: (0 if r.get("combined_z") is not None
                                      else 1,
                                      -(r.get("combined_z") or 0.0)))
    for i, r in enumerate(out, 1):
        r["rank_combined"] = i if r.get("combined_z") is not None else None
    return out, info


def _pose_rows(row, poses):
    """One row per measured pose of a construct, carrying that pose's values."""
    out = []
    for p in poses or []:
        d = dict(row)
        for k, v in (p or {}).items():
            if v is not None:
                d[k] = v
        # The pose's own energy, under the name the construct row uses.
        if p.get("dg") is not None:
            d["mmgbsa"] = p["dg"]
            d["mmgbsa_top5"] = p["dg"]
        out.append(d)
    return out


def spans(rows, poses_of, w=None, info=None):
    """sequence -> (lowest, highest) combined score over its own poses.

    What ``rank.overlapping`` does for MM-GBSA, for the combined score. A
    construct's row describes ONE pose chosen out of several, and the poses
    are tens of kcal/mol and hundreds of square angstroms apart; the order
    between two constructs whose poses cover each other is an artifact of
    reducing several poses to one number, not something the run established.

    The run's scale is held fixed while one construct's poses are swapped
    through it. Re-deriving the z-scores from each slab of poses instead would
    move every construct at once, and the range would then report how much the
    WHOLE RUN varies rather than how much this construct does - a much wider
    interval, and one that says nothing about the pair being compared.

    The criterion stays assumption-free: ranges overlap or they do not.
    Nothing is assumed about the distribution and there is nothing to tune.
    """
    if info is None:
        _scored, info = score([dict(r) for r in rows], w)
    model = info.get("model") or {}
    if not model:
        return {}
    w = w or info.get("weights") or weights()
    energy_used = info.get("used", {}).get("energy") or []

    def one(pose_row):
        """A pose's combined score on the run's own scale."""
        parts, got = {}, False
        for name, m in model.items():
            if not m["sd"]:
                continue
            direction = next(t[2] for t in TERMS if t[0] == name)
            v = _value(pose_row, m["keys"])
            if v is None:
                return None
            if name == "specific":
                extra = _value(pose_row, SPECIFIC_ALSO)
                if "n_ionic" in (info["used"].get("specific") or []):
                    if extra is None:
                        return None
                    v = v + extra
            if direction == "low":
                v = -v
            parts[name] = (v - m["mean"]) / m["sd"]
            got = True
        if not got:
            return None
        total = 0.0
        if energy_used:
            total += w["energy"] * (sum(parts[n] for n in energy_used)
                                    / len(energy_used))
        for name, z in parts.items():
            if name not in energy_used:
                total += w.get(name, 0.0) * z
        return total

    out = {}
    for r in rows:
        key = r.get("sequence") or r.get("peptide")
        vals = [v for v in (one(p) for p in _pose_rows(r, poses_of(r)))
                if v is not None]
        # One pose says nothing about spread, exactly as ``rank.overlapping``
        # requires two before it will group anything.
        if len(vals) > 1:
            out[key] = (min(vals), max(vals))
    return out


def unseparated(rows, pose_spans):
    """Neighbors in the combined order whose ranges cover each other.

    ``rank.overlapping`` with the combined score's own ranges: the mechanism
    is that function's, this only feeds it a different number.
    """
    def values_of(r):
        s = pose_spans.get(r.get("sequence") or r.get("peptide"))
        return list(s) if s else []

    return rank.overlapping(rows, None, key="sequence", value="combined_z",
                            values_of=values_of)


def covers(rows, pose_spans):
    """sequence -> the constructs whose range actually covers its own.

    ``unseparated`` groups a CHAIN: A overlaps B and B overlaps C puts all
    three in one group, which is the right shape for the caution but says
    more than the measurements do if it is read pairwise - A and C may be
    cleanly apart. The per construct column has to be true pairwise, because
    somebody comparing two rows will read it that way, so it is computed
    pairwise here rather than taken from the group.
    """
    out = {}
    for r in rows:
        k = r.get("sequence") or r.get("peptide")
        s = pose_spans.get(k)
        if not s:
            continue
        hit = []
        for other in rows:
            ok = other.get("sequence") or other.get("peptide")
            o = pose_spans.get(ok)
            if ok == k or not o:
                continue
            if not (o[0] > s[1] or o[1] < s[0]):
                hit.append(ok)
        out[k] = hit
    return out


def say(info, log=print):
    """What the combined score was actually computed from, in the log.

    A weighted score that does not say which terms were in it is a number
    nobody can check.
    """
    if info.get("too_few"):
        log("  combined ranking: %d construct(s) is too few for a position "
            "within the run to mean anything, so the order is unchanged"
            % info["n"])
        return
    if not info.get("used"):
        log("  combined ranking: nothing measured for every construct, so "
            "the order is unchanged")
        return
    w = info["weights"]
    log("  combined ranking, from %d term(s):" % len(info["used"]))
    for name in ("energy",) + tuple(t[0] for t in TERMS) + ("clash",):
        if name not in info["used"]:
            continue
        what = next((t[4] for t in TERMS if t[0] == name), None)
        if name == "energy":
            what = "binding energy (%s)" % " and ".join(info["used"]["energy"])
        elif name == "clash":
            what = "steric clashes, as a penalty"
        log("      %-9s weight %-5s %s"
            % (name, "%g" % w.get(name, 0.0), what))
    if info.get("constant"):
        log("      not ranked on (the same for every construct in this run): "
            "%s" % ", ".join(info["constant"]))
    if info.get("off"):
        log("      not ranked on (missing for at least one construct, so "
            "switched off for all): %s" % ", ".join(info["off"]))
    log("      The weights are relative trust, not fitted: there is no "
        "outcome data here to fit them on.")


def top(rows, n, w=None, pose_spans=None, log=None, with_ties=False):
    """The best ``n`` on the combined score, and what the cut does not settle.

    ``rank.with_ties`` hands back the near-ties as extra rows, which is right
    for a table somebody reads. This cut is usually a shortlist for molecular
    dynamics, where every extra construct is days of compute, so the default
    returns exactly ``n`` and NAMES what sits against the line instead of
    quietly lengthening the list. ``with_ties=True`` restores the other
    behavior for a caller that wants the longer table.

    Either way the near-ties are established rather than assumed: a construct
    joins the cut only if its range over its own poses covers the cut row's,
    so the tolerance is measured and not a number somebody typed.
    """
    n = max(0, int(n))
    for r in rows:
        r.setdefault("tied_with", "")
    if n <= 0 or len(rows) <= n:
        return list(rows)
    kept = list(rows[:n])
    cut = rows[n - 1]
    cs = (pose_spans or {}).get(cut.get("sequence") or cut.get("peptide"))
    near = []
    if cs:
        for r in rows[n:]:
            s = (pose_spans or {}).get(r.get("sequence") or r.get("peptide"))
            if not s or s[0] > cs[1] or s[1] < cs[0]:
                continue
            r["tied_with"] = str(n)
            near.append(r)
    if log and near:
        log("  the top %d is a cut somebody asked for, and the measurements "
            "do not place it there: %d construct(s) below it are not "
            "separated from rank %d, because scored over their own poses "
            "their combined scores cover rank %d's."
            % (n, len(near), n, n))
        for r in near:
            log("            %s" % (r.get("sequence") or r.get("peptide")))
        log("            More poses narrows this - it does not remove it.")
    return kept + near if with_ties else kept


# How far each weight is swept when the ordering's sensitivity is measured.
#
# Wide on purpose. A sweep over a narrow band around the chosen weights would
# find the ordering stable whatever the weights were worth, which is the
# answer the question is asking for and so is no answer at all. The range runs
# from a term switched off entirely to a term worth twice the energy, which
# covers every weighting anybody could argue for from this evidence.
SWEEP = (0.0, 0.25, 0.5, 1.0, 1.5, 2.0)


def _order(rows):
    return [r.get("sequence") or r.get("peptide") for r in rows
            if r.get("combined_z") is not None]


def _spearman(a, b):
    """Rank correlation between two orderings of the same constructs."""
    pos = {k: i for i, k in enumerate(b)}
    common = [k for k in a if k in pos]
    n = len(common)
    if n < 3:
        return None
    x = [a.index(k) for k in common]
    y = [pos[k] for k in common]
    mx, my = sum(x) / n, sum(y) / n
    num = sum((p - mx) * (q - my) for p, q in zip(x, y))
    dx = sum((p - mx) ** 2 for p in x) ** 0.5
    dy = sum((q - my) ** 2 for q in y) ** 0.5
    return (num / (dx * dy)) if dx and dy else None


def stability(rows, cuts=(5, 10, 20), sweep=SWEEP, log=None):
    """How much the order actually depends on the weights.

    The weights are not fitted - the module docstring says why - so none of
    them can be presented as validated. What CAN be established without
    fitting anything is whether the choice matters at all: every weight is
    swept over ``sweep`` one at a time, and the resulting order is compared
    with the order the defaults give.

    This measures the ordering's sensitivity to its own weights. It is not a
    check that the ordering is right, and cannot be read as one: a metric
    can be perfectly insensitive to its weights and still rank badly.

    If the top of the table barely moves, the weights are not worth arguing
    about and that is a real finding. If it churns, the honest presentation is
    the overlap rather than an order, and that is a more important finding.

    Returns {"baseline": [...], "runs": [...], "worst": {...}}, where each run
    carries the weight that was changed, the rank correlation against the
    baseline, and how many of each top-X survived.
    """
    base_rows, _info = score([dict(r) for r in rows])
    base = _order(base_rows)
    out = []
    names = ["energy"] + [t[0] for t in TERMS if t[3] is not None] + ["clash"]
    for name in names:
        for v in sweep:
            w = weights(**{name: v})
            got, _i = score([dict(r) for r in rows], w)
            order = _order(got)
            row = {"weight": name, "value": v, "rho": _spearman(base, order)}
            for c in cuts:
                if c <= len(base):
                    row["top%d" % c] = len(set(base[:c]) & set(order[:c]))
            out.append(row)
    live = [r for r in out if r["rho"] is not None]
    worst = min(live, key=lambda r: r["rho"]) if live else {}
    if log:
        log("  stability of the combined order against its weights "
            "(%d construct(s))" % len(base))
        for r in out:
            cells = " ".join("top%d %d/%d" % (c, r["top%d" % c], c)
                             for c in cuts if ("top%d" % c) in r)
            log("      %-9s = %-4g  rho %s   %s"
                % (r["weight"], r["value"],
                   "   -" if r["rho"] is None else "%5.2f" % r["rho"], cells))
        if worst:
            log("      worst case: %s = %g, rho %.2f"
                % (worst["weight"], worst["value"], worst["rho"]))
    return {"baseline": base, "runs": out, "worst": worst}


# ----------------------------------------------------------- reading a run
# A finished run, whichever tool made it. Tools write different tables: some
# put the panel in metric_panel.tsv beside final_constructs.tsv, this one puts
# everything in one results.tsv, and a person with a run directory should not
# have to know which.

#: (the ranked table, the panel beside it or None, what the combined table is
#: called). Order matters: the first that exists is the one used.
LAYOUTS = [
    ("final_constructs.tsv", "metric_panel.tsv",
     "final_constructs_combined.tsv"),
    ("results.tsv", None, "results_combined.tsv"),
]


def _read_tsv(path):
    import csv
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_run(rundir):
    """(rows, original order, the name to write beside it), or (None, ...).

    Rows carry 'sequence' whichever column the table calls it, so everything
    downstream has one name to use.
    """
    import os
    for ranked, panel, out in LAYOUTS:
        p = os.path.join(rundir, ranked)
        if not os.path.isfile(p):
            continue
        rows = _read_tsv(p)
        if not rows:
            return [], [], out
        side = {}
        if panel and os.path.isfile(os.path.join(rundir, panel)):
            for r in _read_tsv(os.path.join(rundir, panel)):
                s = (r.get("sequence") or "").strip()
                if s:
                    side[s] = r
        for r in rows:
            s = (r.get("sequence") or r.get("peptide") or "").strip()
            r["sequence"] = s
            # The ranked table wins where both carry a column: it is the one
            # place those numbers are published, and two files disagreeing
            # about a construct's MM-GBSA is worse than either being wrong.
            for k, v in (side.get(s) or {}).items():
                if v not in ("", None) and r.get(k) in ("", None):
                    r[k] = v
        return rows, [r["sequence"] for r in rows], out
    return None, [], None


def poses_from(rundir):
    """A reader for the per pose measurements the panel wrote, for ``spans``.

    Returns a function, and an empty list for any construct whose poses were
    never written - a run from before the panel existed simply has no ranges
    and therefore no claim about what is separated.
    """
    import json
    import os

    def poses_of(row):
        p = os.path.join(rundir, "metrics", "%s.json" % row.get("sequence"))
        if not os.path.isfile(p):
            return []
        try:
            return json.load(open(p)).get("poses") or []
        except (OSError, ValueError):
            return []
    return poses_of


def write_table(rows, path, info, original=None, spans_by=None):
    """The combined order, as a table beside the one the run already has.

    Every term that went into the score gets its own column. A single number
    nobody can take apart is not a ranking anybody can argue with, and the
    person reading has to be able to see that a construct is third because of
    its buried area rather than because of its energy.
    """
    # Rows first, path second, as the other writers in this engine have it.
    # Said out loud because the two are easy to swap and the swap
    # fails deep inside with an AttributeError about a string, which says
    # nothing about what the caller got wrong.
    if isinstance(rows, str) or not isinstance(path, str):
        raise TypeError("write_table(rows, path, ...): the rows come first")
    terms = sorted({k for r in rows for k in (r.get("combined_terms") or {})})
    was = {s: i for i, s in enumerate(original or [], 1)}
    # The measurements themselves, under whichever name this run's table gave
    # them, so the columns beside the z-scores are the numbers a reader can
    # check against the run's own file rather than a renamed copy of them.
    raw = []
    for _n, keys, _d, _w, _what in TERMS:
        for k in keys:
            if any(r.get(k) not in (None, "") for r in rows):
                raw.append(k)
    for k in SPECIFIC_ALSO + ("n_clash",):
        if any(r.get(k) not in (None, "") for r in rows):
            raw.append(k)
    cols = (["rank_combined", "rank_before", "moved", "sequence",
             "combined_z"]
            + ["z_" + t for t in terms] + raw
            + ["combined_lo", "combined_hi", "not_separated"])
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i, r in enumerate(rows, 1):
            before = was.get(r.get("sequence"))
            sp = (spans_by or {}).get(r.get("sequence"))
            cell = {
                "rank_combined": "" if r.get("combined_z") is None else i,
                "rank_before": before or "",
                "moved": "" if not before else ("0" if before == i
                                                else "%+d" % (before - i)),
                "sequence": r.get("sequence", ""),
                "combined_z": ("" if r.get("combined_z") is None
                               else "%.3f" % r["combined_z"]),
                "combined_lo": "" if not sp else "%.3f" % sp[0],
                "combined_hi": "" if not sp else "%.3f" % sp[1],
                "not_separated": r.get("not_separated", ""),
            }
            for t in terms:
                v = (r.get("combined_terms") or {}).get(t)
                cell["z_" + t] = "" if v is None else "%+.3f" % v
            for k in raw:
                cell[k] = "" if r.get(k) in (None, "") else str(r[k])
            fh.write("\t".join(str(cell.get(c, "")) for c in cols) + "\n")
    return path


def reorder(rundir, n=0, w=None, out=None, log=print, with_ties=False):
    """Rank a finished run on the combined metric, beside its own table.

    Nothing the run already wrote is read back out and overwritten. Numbers
    from these runs have been quoted in reports and in conversation, and a
    table that silently says something different the next time somebody opens
    it is worse than no table.
    """
    import os
    rows, original, name = read_run(rundir)
    if rows is None:
        log("no ranked table in %s" % rundir)
        return None
    if not rows:
        log("%s has an empty table, so there is nothing to reorder" % rundir)
        return None
    done, info = score(rows, w)
    say(info, log)
    if len(info.get("used") or {}) < 2:
        # One axis is not a combination. Writing the energy order out under
        # the name "combined" would hand somebody a second file that looks
        # like a second opinion and is the first one renamed.
        log("  only one measurement is available for every construct, so "
            "there is nothing to combine: no table is written.")
        log("  the metric panel measures the interface of %s from the files "
            "MM-GBSA already left." % rundir)
        return None
    sp = spans(done, poses_from(rundir), w, info)
    groups = unseparated(done, sp) if sp else []
    pairs = covers(done, sp) if sp else {}
    for r in done:
        r["not_separated"] = ", ".join(pairs.get(r.get("sequence")) or [])
    path = out or os.path.join(rundir, name)
    write_table(done, path, info, original, sp)
    log("  combined ranking -> %s" % path)
    log("  the run's own table is unchanged.")
    moved = sum(1 for a, b in zip(original,
                                  [r.get("sequence") for r in done]) if a != b)
    log("  %d of %d construct(s) are in a different position from the run's "
        "own order." % (moved, len(original)))
    if groups:
        log("  CAUTION %d construct(s) fall in %d chain(s) of neighbors the "
            "combined score does not separate: scored over their own poses "
            "each one's score covers the next one's, so a step of one place "
            "inside a chain is not something the run showed."
            % (sum(len(g) for g in groups), len(groups)))
        for g in groups:
            log("            %s" % " > ".join(x.get("sequence") for x in g))
        log("            Ends of a chain may still be apart - the "
            "not_separated column says which pairs actually overlap.")
    if n:
        picked = top(done, n, w, sp, log, with_ties)
        log("  top %d on the combined score:" % n)
        for i, r in enumerate(picked, 1):
            log("      %2d  %-24s %7.3f%s"
                % (i, r.get("sequence"), r.get("combined_z") or 0.0,
                   "   (not separated from the cut)" if r.get("tied_with")
                   else ""))
    return path


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        description="Rank a finished run on the combined metric: the binding "
                    "energy together with the interface measurements the "
                    "panel already made. Writes a table BESIDE the run's own; "
                    "nothing already written is changed.")
    ap.add_argument("rundir")
    ap.add_argument("--top", type=int, default=0, metavar="X",
                    help="also print the best X constructs")
    ap.add_argument("--with-ties", action="store_true",
                    help="include constructs the measurements cannot "
                         "separate from the cut in the top X")
    ap.add_argument("--out", help="write the table here instead of beside "
                                  "the run")
    ap.add_argument("--stability", action="store_true",
                    help="sweep the weights and report how much the order "
                         "depends on them")
    for name, _keys, _d, default, what in TERMS:
        if default is not None:
            ap.add_argument("--w-" + name, type=float, default=None,
                            metavar="W", help="weight on %s (default %g)"
                                              % (what, default))
    ap.add_argument("--w-energy", type=float, default=None, metavar="W",
                    help="weight on the binding energy (default %g)"
                         % ENERGY_WEIGHT)
    ap.add_argument("--w-clash", type=float, default=None, metavar="W",
                    help="penalty on steric clashes (default %g)"
                         % CLASH_PENALTY)
    a = ap.parse_args(argv)
    over = {k[2:]: v for k, v in vars(a).items()
            if k.startswith("w_") and v is not None}
    w = weights(**over)
    if a.stability:
        rows, _o, _n = read_run(a.rundir)
        if not rows:
            print("no ranked table in %s" % a.rundir)
            return 1
        stability(rows, log=print)
        return 0
    return 0 if reorder(a.rundir, a.top, w, a.out,
                        with_ties=a.with_ties) else 1


if __name__ == "__main__":
    raise SystemExit(main())


def apply_order(path, order, rank_col="rank", keep_as="rank_binding"):
    """Put an existing table into a new order without rebuilding it.

    The rows are moved and the rank column renumbered; every other column is
    carried across as the characters that were already in the file. A table
    rewritten from parsed values loses whatever the writer knew about
    formatting - how many decimals a column had, which blanks meant "not
    measured" - and this one is read by people and by three other modules.
    The position the row held before is kept in `keep_as`, so the previous
    ordering is still in the file rather than only in a file beside it.

    Only ever used on a table this run has just written. An older run's table
    is never touched: its numbers have been reported.
    """
    with open(path) as fh:
        lines = fh.read().splitlines()
    if len(lines) < 2:
        return path
    head = lines[0].split("\t")
    if keep_as in head:
        return path
    # Some tools call the column "sequence" and this one calls it "peptide".
    # Looking for one name only made this a no-op on half the runs it is for, and a no-op is the worst failure available here: the
    # table is left in the old order with nothing said about it.
    seq_at = next((head.index(c) for c in ("sequence", "peptide")
                   if c in head), None)
    if seq_at is None:
        return path
    rank_at = head.index(rank_col) if rank_col in head else None
    body = {}
    for line in lines[1:]:
        if line.strip():
            body[line.split("\t")[seq_at]] = line.split("\t")
    # Anything the new order does not mention keeps its place at the end,
    # rather than being dropped from a table it was already in.
    rest = [s for s in body if s not in order]
    out = [lines[0] + "\t" + keep_as]
    for i, s in enumerate(list(order) + rest, 1):
        cells = body.get(s)
        if cells is None:
            continue
        was = cells[rank_at] if rank_at is not None else ""
        if rank_at is not None:
            cells = list(cells)
            cells[rank_at] = str(i)
        out.append("\t".join(cells) + "\t" + was)
    with open(path, "w") as fh:
        fh.write("\n".join(out) + "\n")
    return path
