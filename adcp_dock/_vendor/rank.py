"""Ranking. Two objectives are maximised jointly: binding affinity (lowest Vina
score) and sequence-specific engagement (most side-chain contacts to the site).

These do not conflict - on a fixed-length library they correlate positively - so
a plain equal-weight z-score sum is appropriate. Ligand efficiency is deliberately
NOT used: on a library where every member has the same backbone and residue count,
dividing by heavy-atom count removes the very quantity being selected, and in
practice it ranks glycine-rich peptides whose contacts are almost entirely backbone.
"""
from __future__ import annotations
import numpy as np

def zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, float)
    s = v.std()
    return (v - v.mean()) / s if s > 0 else np.zeros_like(v)

# What a contact beside the hotspot is worth, relative to the hotspot itself.
#
# A neighbour of a hotspot is part of the same epitope patch and a peptide
# touching it IS engaging the site - so it counts, but never as much as the
# residue somebody actually asked for. sc_near is already graded from 1 at the
# hotspot to 0 at the edge of the shell, so this is the weight on that graded
# quantity rather than on a flat ring.
NEAR_WEIGHT = 0.5

# And what it is worth to reach two hotspots rather than one twice.
#
# A peptide spanning two hotspots is a different molecule from one perched on
# a single one: it engages more of the interface, and the residues it touches
# on the way are the rest of that interface rather than incidental neighbours.
# Nothing scored this at all - a pose touching one hotspot twenty times and a
# pose touching two hotspots ten times each were indistinguishable.
BRIDGE_WEIGHT = 0.35


def joint(rows: list[dict], w_score: float = 1.0, w_contact: float = 1.0,
          shell_weight: float = 0.5, w_mmgbsa: float = 1.0,
          near_weight: float = NEAR_WEIGHT,
          bridge_weight: float = BRIDGE_WEIGHT) -> list[dict]:
    """rows need 'score', 'sc_hot', 'sc_shell'. Adds 'z' and sorts best first.

    sc_near and bridge are used when every row carries them - older runs and
    anything that scored its poses before they existed keep the flat shell
    term exactly as it was.
    """
    if not rows:
        return []
    sc = np.array([r["score"] for r in rows], float)
    # A row whose contacts were never measured cannot contribute a contact
    # term, and must not be given one. The pose could not be read, so `m` was
    # empty; filling the missing terms with 0 would make "not measured" and
    # "measured as touching nothing" the same row, and would let a fabricated
    # zero be z-scored against real values. They are None instead, and one
    # None turns the term off for everybody - the same all-or-nothing rule the
    # graded and bridging terms already use, and for the same reason: a
    # z-score over a subset puts scored and unscored rows on different scales.
    measured = all(r.get("sc_hot") is not None and r.get("sc_shell") is not None
                   for r in rows)
    graded = measured and all(r.get("sc_near") is not None for r in rows)
    if not measured:
        z = w_score * zscore(-sc)
        for r, zz in zip(rows, z):
            r["z"] = float(zz)
        return sorted(rows, key=lambda r: -r["z"])
    if graded:
        # The graded neighbourhood replaces the flat one rather than adding to
        # it: they measure the same thing, one of them badly.
        ct = np.array([r["sc_hot"] + near_weight * float(r["sc_near"])
                       for r in rows], float)
    else:
        ct = np.array([r["sc_hot"] + shell_weight * r["sc_shell"]
                       for r in rows], float)
    z = w_score * zscore(-sc) + w_contact * zscore(ct)
    if measured and all(r.get("bridge") is not None for r in rows):
        br = np.array([float(r.get("bridge") or 0) for r in rows], float)
        if br.any():
            z = z + bridge_weight * zscore(br)
    # MM-GBSA is the more physical estimate of binding, so when it is
    # available for every row it belongs in the ranking rather than
    # sitting in the table being ignored. Only when every row has one:
    # it is computed for the top n_mmgbsa constructs, and ranking a
    # partly-scored set on it would compare unlike with unlike.
    gb = [r.get("mmgbsa") for r in rows]
    scored = [g for g in gb if g is not None]
    # The z term only when EVERY row has one: a z-score over a subset puts
    # the scored and the unscored rows on different scales, and adding it to
    # some rows and not others would rank them against each other on it.
    if scored and len(scored) == len(gb):
        z = z + w_mmgbsa * zscore(-np.array(gb, float))
    for r, zz in zip(rows, z):
        r["z"] = float(zz)
    if scored:
        # Most negative MM-GBSA first, z breaking ties.
        #
        # Requiring EVERY row to have one before ranking on it threw the
        # whole ordering away whenever a single construct was missing a
        # value - and one is missing whenever n_mmgbsa is less than
        # n_constructs, which is the ordinary way to save an hour a
        # construct on a large receptor, or whenever one MM-GBSA fails.
        # The table then reverted to z-ranking with the MM-GBSA column
        # sitting in it being ignored, silently: runs/colabrun ranks
        # -40.36 above -45.17 for exactly this reason.
        #
        # Rows without a value sort after the ones with it rather than
        # among them. They have not been shown to be worse - they have not
        # been measured - and interleaving them on a z that means something
        # different would assert an order the run never established.
        return sorted(rows, key=lambda r: (
            0 if r.get("mmgbsa") is not None else 1,
            r["mmgbsa"] if r.get("mmgbsa") is not None else 0.0,
            -r["z"]))
    return sorted(rows, key=lambda r: -r["z"])


# Side-chain heavy atoms. Vina's score grows with buried surface, so on a
# fixed-length library the raw score largely reports how big the side chains
# are: Trp beats Leu almost regardless of the pocket. Dividing by that gives
# binding per atom, which is what shape complementarity looks like.
#
# Reported alongside the raw score, never instead of it: efficiency on its own
# favours small glycine-rich peptides whose contacts are mostly backbone, which
# is why this module did not use it for the primary ranking. Two columns, both
# visible, is the honest presentation.
_SC_HEAVY = {"G": 0, "A": 1, "S": 2, "C": 2, "T": 3, "P": 3, "V": 3,
             "N": 4, "D": 4, "I": 4, "L": 4, "M": 4, "E": 5, "Q": 5,
             "K": 5, "H": 6, "F": 7, "R": 7, "Y": 8, "W": 10}


def overlapping(rows, detail, key="sequence", value="mmgbsa",
                values_of=None):
    """Neighbours in the ranking whose pose energies overlap, as groups.

    A construct's MM-GBSA is one number reduced from several poses, and those
    poses can be tens of kcal/mol apart. When two constructs' pose energies
    cover each other, the order between them is an artefact of the reduction,
    not something the calculation showed - and the table prints them one above
    the other with no hint of it.

    That is not a fault in this code and cannot be fixed by computing harder.
    Two independent MMPBSA calculations of the same system are documented to
    differ by more than 10 kcal/mol for small molecules and by up to 43 for
    flexible peptide ligands such as peptide-MHC complexes (Wan & Coveney,
    Phil. Trans. R. Soc. A 379:20200082, 2021). Measured here on this
    pipeline's own output: the same three constructs re-docked on the same grid
    moved by up to 27.3 kcal/mol and came back in the opposite order.

    The criterion is deliberately assumption-free - ranges overlap, or they do
    not. No standard error, no distribution assumed, nothing to tune. It
    understates rather than overstates how much is separated, which is the
    right direction for a caveat.

    ``detail`` is what mmgbsa.run fills: {sequence: {"values": [...]}}.
    Returns a list of groups, each a list of the row dicts that cannot be
    ordered against each other. Constructs with fewer than two poses are never
    grouped, because one pose says nothing about spread.
    """
    def span(r):
        if values_of is not None:
            v = [float(x) for x in (values_of(r) or []) if x is not None]
        else:
            d = (detail or {}).get(r.get(key)) or {}
            v = [float(x) for x in (d.get("values") or []) if x is not None]
        return (min(v), max(v)) if len(v) > 1 else None

    groups, run = [], []
    for r in rows:
        s_r = span(r)
        if s_r is None or r.get(value) is None:
            if len(run) > 1:
                groups.append([x for x, _ in run])
            run = []
            continue
        if run and not (s_r[0] > run[-1][1][1] or s_r[1] < run[-1][1][0]):
            run.append((r, s_r))
        else:
            if len(run) > 1:
                groups.append([x for x, _ in run])
            run = [(r, s_r)]
    if len(run) > 1:
        groups.append([x for x, _ in run])
    return groups


def say_overlaps(rows, detail, log=print, key="sequence", value="mmgbsa",
                 values_of=None, what="pose energies"):
    """Log which neighbours the measurements do not separate."""
    groups = overlapping(rows, detail, key, value, values_of)
    if not groups:
        return groups
    n = sum(len(g) for g in groups)
    log(f"  CAUTION {n} construct(s) are not separated by {value}: "
        f"their {what} overlap, so the order within each group below is "
        f"not something the calculation showed.")
    for g in groups:
        lo = min(float(x[value]) for x in g)
        hi = max(float(x[value]) for x in g)
        log(f"            {', '.join(str(x.get(key)) for x in g)} "
            f"({hi:.1f} to {lo:.1f} kcal/mol)")
    log("            More poses narrows this - it does not remove it. "
        "Re-run to see which orderings hold.")
    return groups


def sc_heavy(seq: str) -> int:
    """Side-chain heavy atoms in a sequence."""
    return sum(_SC_HEAVY.get(c, 4) for c in (seq or "").upper())


def efficiency(rows, score_key="affinity", seq_key="sequence"):
    """Add sc_heavy, efficiency and rank_efficiency. Row order is unchanged."""
    for r in rows:
        n = sc_heavy(r.get(seq_key, ""))
        sc = r.get(score_key)
        r["sc_heavy"] = n
        r["efficiency"] = (float(sc) / n) if (n and sc is not None) else None
    have = [r for r in rows if r.get("efficiency") is not None]
    for i, r in enumerate(sorted(have, key=lambda x: x["efficiency"]), 1):
        r["rank_efficiency"] = i
    return rows


# How close two constructs have to be before the table refuses to separate them.
#
# Both measures, not either: two designs that agree on the docking score and
# disagree on MM-GBSA are not tied, they are two different answers, and the
# table should keep the order it has rather than pretend the question is open.
TIE_AFFINITY = 2.0        # kcal/mol, ADCP
TIE_MMGBSA = 2.0          # kcal/mol, MM-GBSA


def with_ties(rows, n, tie_affinity=TIE_AFFINITY, tie_mmgbsa=TIE_MMGBSA,
              log=None):
    """The top ``n``, plus anything the measurements cannot separate from them.

    Asking for ten designs and being handed exactly ten hides the thing worth
    knowing: that the eleventh was indistinguishable from the tenth. The cut is
    a number somebody typed, and the measurements do not respect it - so where
    a row below the line is within tolerance of the last row above it on BOTH
    the docking score and MM-GBSA, it comes too, marked, and the person reading
    decides what to do about it.

    That is the opposite of padding the table. Only rows that reach the cut's
    own scores are kept, the reason is written in `tied_with`, and a table with
    no near-ties is exactly ``n`` rows long.

    ``rows`` must already be in the order the run reports. Returns the kept
    rows; each carries `tied_with` - "" for the ones that made the cut on their
    own, and the rank they are tied with for the rest.
    """
    n = max(0, int(n))
    for r in rows:
        r.setdefault("tied_with", "")
    if n <= 0 or len(rows) <= n:
        return list(rows)

    cut = rows[n - 1]
    ca, cg = cut.get("affinity"), cut.get("mmgbsa")
    kept = list(rows[:n])
    for i, r in enumerate(rows[n:], start=n + 1):
        a, g = r.get("affinity"), r.get("mmgbsa")
        if a is None or ca is None or abs(float(a) - float(ca)) > tie_affinity:
            continue
        # MM-GBSA has to agree too, and a row that has none cannot be shown to
        # be tied - it is simply unmeasured, which is not the same thing.
        if g is None or cg is None or abs(float(g) - float(cg)) > tie_mmgbsa:
            continue
        r["tied_with"] = str(n)
        kept.append(r)
    if log and len(kept) > n:
        log(f"  {len(kept) - n} construct(s) below the top {n} are within "
            f"{tie_affinity:g} kcal/mol of rank {n} on the docking score AND "
            f"{tie_mmgbsa:g} on MM-GBSA, so they are reported too - the cut "
            f"is a number that was asked for, not one the measurements "
            f"support.")
        for r in kept[n:]:
            log(f"            {r.get('sequence', '?')} "
                f"(affinity {r.get('affinity')}, mmgbsa {r.get('mmgbsa')})")
    return kept
