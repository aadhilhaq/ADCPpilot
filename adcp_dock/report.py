"""Turning a run into readable, sortable numbers.

By Aadhil Haq

The headline figure for a peptide is the mean of its top few poses rather
than its single best. A best score is the tail of a distribution the search
drew many times and rewards the luckiest replica; the mean of the top few
is steadier. Both are reported, since comparison with other tools needs the
best as well.
"""
import json
import os
import statistics as st

TOP_N = 5


#: The MM-GBSA engines a row can carry values from, in column order. Each
#: gets its own columns so their numbers are never read against one another:
#: they are on different scales, AmberTools and Prime 25 kcal/mol apart on the
#: same pose.
GB_ENGINES = ("amber", "openmm", "schrodinger")
GB_LABELS = {"amber": "AmberTools", "openmm": "OpenMM",
             "schrodinger": "Prime"}


def summarise(peptide, modes, mm=None, top_n=TOP_N, name="", openmm=None,
              engine="", panel=None, clusters=None, mmgbsa_engine="",
              mm_by_engine=None):
    """One row. `mm` is [dG] for the poses that were rescored, best-first.

    The name is carried from the result file rather than derived from the
    sequence: it is what the user typed beside the sequence, and it is what
    distinguishes two variants of the same construct.
    """
    scores = [m["score"] for m in modes or []]
    top = scores[:top_n]
    row = {
        "name": name or "",
        "peptide": peptide,
        # Residues, not characters: an ADCP 1.1 sequence such as N<BHD>PIS&DVD
        # is seven residues.
        "length": len(_residues(peptide)),
        "engine": engine or "",
        # The DOCKING engine is "engine" above. This is the program that
        # computed the MM-GBSA numbers, which is a different question with a
        # different answer, and reading one for the other labels an energy
        # with the name of a docking program. Runs from before engines were
        # selectable all used AmberTools.
        "mmgbsa_engine": (mmgbsa_engine or ("amber" if mm else "")),
        "poses": len(scores),
        "adcp_best": scores[0] if scores else None,
        "adcp_top%d" % top_n: round(st.mean(top), 2) if top else None,
        "adcp_sd": round(st.pstdev(top), 2) if len(top) > 1 else None,
    }
    if mm:
        vals = [v for v in mm if v is not None]
        mtop = vals[:top_n]
        row.update({
            "mmgbsa_n": len(vals),
            "mmgbsa_best": round(min(vals), 2) if vals else None,
            "mmgbsa_top%d" % top_n: round(st.mean(mtop), 2) if mtop else None,
            "mmgbsa_median": round(st.median(mtop), 2) if mtop else None,
            "mmgbsa_sd": round(st.pstdev(mtop), 2) if len(mtop) > 1 else None,
            # Most positive minus most negative across the scored
            # poses. It shows how much the poses disagree; it is
            # not a statistical test of the ranking.
            "mmgbsa_spread": (round(max(vals) - min(vals), 2)
                              if len(vals) > 1 else None),
        })
    # Every engine that has scored this peptide, each in its own columns.
    # These come from the caches on disk, so an engine run on a past job
    # afterwards appears here without the result file having to change.
    for eng in GB_ENGINES:
        vals = [v for v in ((mm_by_engine or {}).get(eng) or [])
                if v is not None]
        if not vals:
            continue
        mtop = vals[:top_n]
        row["mmgbsa_%s_n" % eng] = len(vals)
        row["mmgbsa_%s_best" % eng] = round(min(vals), 2)
        row["mmgbsa_%s_top%d" % (eng, top_n)] = round(st.mean(mtop), 2)
    if openmm:
        # ADCP 1.1's OpenMM minimization: E(complex) minus E(receptor) after
        # minimizing, per pose. Reported beside the scores, not ranked by.
        best = min(openmm, key=lambda o: o["e_complex_receptor"])
        row.update({"omm_n": len(openmm),
                    "omm_best": round(best["e_complex_receptor"], 1),
                    "omm_best_pose": best["rank_adcp"]})
    if clusters and clusters.get("cluster_total"):
        # How many of the search's clustered solutions fall in the top pose's
        # cluster, and what share of the clustered solutions that is. ADCP's
        # poses are cluster representatives, so this is the search's own
        # measure of convergence, and how much the poses differ from one
        # another is not. ADCP prints at most 100 modes, so the share is an
        # upper bound: the solutions in the modes it did not print are missing
        # from the denominator.
        row["cluster_top"] = clusters.get("cluster_top")
        row["cluster_modes"] = clusters.get("cluster_modes")
        row["cluster_frac"] = round(clusters["cluster_top"]
                                    / float(clusters["cluster_total"]), 3)
    if panel:
        row.update({c: panel.get(k) for c, k in PANEL_KEYS})
    return row


# The panel's columns, as (column, key in the per pose metrics). The energy
# terms, the site and the interface describe ONE pose: the lowest energy pose
# among those the panel could measure, or the top ADCP pose when there is no
# MM-GBSA. panel_pose says which it was, panel_n how many were measured, and
# panel_geometry which coordinates they were measured on.
PANEL_KEYS = [
    ("panel_pose", "pose"),
    ("mm_vdw", "mm_vdw"), ("mm_eel", "mm_eel"), ("mm_egb", "mm_egb"),
    ("mm_esurf", "mm_esurf"),
    ("site_n", "site_n"), ("site_touched", "site_touched"),
    ("site_coverage", "site_coverage"),
    ("site_dg", "site_dg"), ("site_dg_mean", "site_dg_mean"),
    ("bsa", "bsa"), ("bsa_per_res", "bsa_per_res"), ("bsa_sd", "bsa_sd"),
    ("bsa_basis", "bsa_basis"),
    ("f_pep_buried", "f_pep_buried"), ("f_hydrophobic", "f_hydrophobic"),
    ("n_hbond", "n_hbond"), ("n_ionic", "n_ionic"), ("n_stack", "n_stack"),
    ("n_contact", "n_contact"), ("n_clash", "n_clash"),
    ("gravy", "gravy"), ("charge_ph7", "charge_ph7"),
    ("apolar_run", "apolar_run"), ("pe_score", "pe_score"),
    ("pe_verdict", "verdict"), ("pe_flags", "pe_flags"),
    ("panel_geometry", "geometry"),
    # What the minimisation before the energy did. sander stops at maxcyc
    # rather than at a gradient for almost every pose, which is expected of a
    # restrained rescoring minimisation and is not the same as converging, so
    # both facts are published instead of a verdict.
    ("panel_n", "panel_n"),
    ("min_cycles", "min_cycles"), ("min_rms_grad", "min_rms_grad"),
    ("min_stopped_at_maxcyc", "min_stopped_at_maxcyc"),
]
PANEL_COLS = [c for c, _k in PANEL_KEYS]


def _residues(seq):
    from . import peptide
    return peptide.residues(seq)


def combined_scores(rows):
    """Attach the combined score and its components to every row.

    Computed whichever order the table is printed in. The panel measures
    buried area, hydrogen bonds, salt bridges, site coverage and the per
    residue decomposition for every pose, and until now none of it reached
    the ranking; these columns are the record of what it said. Returns the
    info block, or None when nothing could be scored.
    """
    from ._vendor import combined
    for r in rows:
        r.setdefault("sequence", r.get("peptide"))
    try:
        done, info = combined.score(rows)
    except Exception:                                           # noqa: BLE001
        # A second ordering is a convenience. A run that has produced a table
        # must not be failed because an extra column could not be computed.
        return None
    if len(info.get("used") or {}) < 2:
        # One axis is not a combination. A run with no panel can still have
        # its energy z-scored, and the result would be the single-column
        # order under a second name - a file and three columns implying a
        # second opinion that was never formed. Nothing is published.
        for r in rows:
            for k in ("combined_z", "combined_terms", "combined_off"):
                r.pop(k, None)
        return None
    order = {id(r): i for i, r in enumerate(done, 1)}
    for r in rows:
        r["rank_combined"] = (order[id(r)] if r.get("combined_z") is not None
                              else None)
        for k, v in (r.get("combined_terms") or {}).items():
            r["z_" + k] = round(v, 3)
        if r.get("combined_z") is not None:
            r["combined_z"] = round(r["combined_z"], 3)
    return info


def rank(rows, by=None, top_n=TOP_N, on_combined=False):
    """Most negative first. MM-GBSA decides when it was run, because that is
    the more considered estimate; ADCP's own score when it was not.

    `on_combined` orders on the combined metric instead: the same energy
    together with the interface measurements the panel already made. Off by
    default, so a run reads the way every run before it did. The combined
    columns are written either way - both orders are always in the file, and
    this only decides which one the rows are sorted into.
    """
    if on_combined and any(r.get("combined_z") is not None for r in rows):
        key = "combined_z"
        have = [r for r in rows if r.get(key) is not None]
        gone = [r for r in rows if r.get(key) is None]
        # Higher is better here, the opposite of every other column in this
        # table, because it is a z-score and not an energy.
        have.sort(key=lambda r: -r[key])
    else:
        key = by or ("mmgbsa_top%d" % top_n
                     if any(r.get("mmgbsa_top%d" % top_n) is not None
                            for r in rows)
                     else "adcp_top%d" % top_n)
        have = [r for r in rows if r.get(key) is not None]
        gone = [r for r in rows if r.get(key) is None]
        have.sort(key=lambda r: r[key])
    for i, r in enumerate(have, 1):
        r["rank"] = i
    for r in gone:
        r["rank"] = None
    return have + gone, key


def table(rows, top_n=TOP_N, on_combined=False):
    """Fixed-width, for a terminal or a log."""
    rows, key = rank(rows, top_n=top_n, on_combined=on_combined)
    mm = any(r.get("mmgbsa_n") for r in rows)
    omm = any(r.get("omm_n") for r in rows)
    head = ["#", "name", "peptide", "len", "poses", "ADCP best",
            "ADCP top%d" % top_n, "sd"]
    if mm:
        head += ["MM n", "MM best", "MM top%d" % top_n, "MM med", "sd",
                 "spread"]
    if omm:
        head += ["OpenMM best"]
    # The combined score sits at the end whether or not the table is ordered
    # on it. Printed where a run is ordered the old way it is the only hint
    # that a second answer exists; printed where it is the sort key it is the
    # column somebody will want to check the order against.
    comb = any(r.get("combined_z") is not None for r in rows)
    if comb:
        head += ["combined"]
    # Column widths are set from the content: a fixed name column either
    # truncates long names or leaves a gap on runs where nothing was
    # named.
    w = [max(4, len(head[0])),
         max([len(str(r.get("name") or "-")) for r in rows] + [len(head[1])])
         + 1 if rows else 10,
         (max(len(str(r["peptide"])) for r in rows) + 1) if rows else 10]
    out = ["ranked by %s, %s first"
           % (key, "highest" if key == "combined_z" else "most negative"), ""]
    fmt = ("%-4s %-*s %-*s %4s %5s %10s %11s %6s"
           + ("  %5s %8s %10s %8s %6s %7s" if mm else "")
           + ("  %11s" if omm else "")
           + ("  %8s" if comb else ""))

    def cell(v, nd=2):
        return "-" if v is None else ("%.*f" % (nd, v) if isinstance(v, float)
                                      else str(v))
    out.append(fmt % tuple([head[0], w[1], head[1], w[2], head[2]]
                           + head[3:]))
    for r in rows:
        vals = [cell(r.get("rank"), 0), w[1], (r.get("name") or "-"),
                w[2], r["peptide"], r["length"],
                r["poses"], cell(r.get("adcp_best"), 1),
                cell(r.get("adcp_top%d" % top_n), 2), cell(r.get("adcp_sd"), 2)]
        if mm:
            vals += [cell(r.get("mmgbsa_n"), 0), cell(r.get("mmgbsa_best"), 2),
                     cell(r.get("mmgbsa_top%d" % top_n), 2),
                     cell(r.get("mmgbsa_median"), 2), cell(r.get("mmgbsa_sd"), 2),
                     cell(r.get("mmgbsa_spread"), 2)]
        if omm:
            vals += [("%.1f (pose %s)" % (r["omm_best"], r["omm_best_pose"]))
                     if r.get("omm_n") else "-"]
        if comb:
            vals += [cell(r.get("combined_z"), 2)]
        out.append(fmt % tuple(vals))
    # An interface, not merely a pose: a run whose geometry was not kept still
    # has its energy terms, and they belong in the file rather than in a
    # printed block of dashes.
    if any(r.get("bsa") is not None for r in rows):
        out += [""] + panel_table(rows)
    out += _geometry_notes(rows)
    unranked = [r for r in rows if r.get("rank") is None]
    if unranked:
        out.append("")
        out.append("Unranked (no %s): %s. Peptides with non-standard or D "
                   "residues have no MM-GBSA value." % (key, ", ".join(
                       str(r.get("name") or r["peptide"]) for r in unranked)))
    return "\n".join(out)


def _geometry_notes(rows):
    """Two different cautions, which it is tempting to write as one.

    A pose scored without minimising is a caution about its ENERGY: that
    number is a single point on a complex the docking left where it was, and
    it does not belong in a column beside the others. A geometry that was not
    kept is a caution about the PANEL only, and says nothing at all about the
    energy beside it.
    """
    out = []
    raw = [r for r in rows
           if "without minimisation" in (r.get("panel_geometry") or "")]
    gone = [r for r in rows if "not kept" in (r.get("panel_geometry") or "")]
    if raw:
        out += ["", "Scored on a complex that was not minimised, so these "
                    "energies are not comparable with the",
                "others: %s." % ", ".join(str(r.get("name") or r["peptide"])
                                          for r in raw[:6])]
    if gone:
        out += ["", "The geometry these energies were computed on is not kept "
                    "with this run, so only the",
                "energy terms could be read and the interface columns are "
                "empty."]
    return out


def _tsv(v):
    """A cell. A list (the developability flags) is joined with semicolons,
    since a tab or a newline inside it would break the row."""
    if v is None:
        return ""
    if isinstance(v, (list, tuple)):
        return "; ".join(str(x) for x in v)
    return str(v)


def panel_table(rows):
    """The interface and sequence panel, one line per peptide in rank order.

    Printed as a block of its own under the ranking: the ranking is what a
    run is read for, and twenty more columns beside it would push the scores
    off a terminal.
    """
    head = ("#", "name", "pose", "BSA", "buried", "apolar", "HB", "SB",
            "stack", "clash", "site", "site dG", "GRAVY", "charge", "PE")
    wn = max([len(str(r.get("name") or r["peptide"])) for r in rows]
             + [len(head[1])]) + 1
    fmt = ("%-4s %-*s %4s %7s %6s %6s %3s %3s %5s %5s %5s %8s %6s %6s %5s")

    def c(v, nd=2):
        if v is None or v == "":
            return "-"
        return ("%.*f" % (nd, v)) if isinstance(v, float) else str(v)
    out = ["Interface and sequence panel (the pose each row describes is "
           "in the pose column)", ""]
    out.append(fmt % ((head[0], wn, head[1]) + head[2:]))
    for r in rows:
        if r.get("panel_pose") is None:
            continue
        site = ("%d/%d" % (r["site_touched"], r["site_n"])
                if r.get("site_n") and r.get("site_touched") is not None
                else "-")
        out.append(fmt % (c(r.get("rank"), 0), wn,
                          r.get("name") or r["peptide"],
                          c(r.get("panel_pose")), c(r.get("bsa"), 0),
                          c(r.get("f_pep_buried"), 2),
                          c(r.get("f_hydrophobic"), 2),
                          c(r.get("n_hbond")), c(r.get("n_ionic")),
                          c(r.get("n_stack")), c(r.get("n_clash")), site,
                          c(r.get("site_dg"), 1), c(r.get("gravy"), 2),
                          c(r.get("charge_ph7"), 1),
                          c(r.get("pe_verdict"))))
    # What site dG is, and why the breakdown does not add up exactly to the
    # MM-GBSA value, is explained in the documentation (results.html, "The
    # limits on these columns"). results.txt is read by the people the run
    # was for, and carries their numbers, not notes on method.
    return out


def write(rows, outdir, top_n=TOP_N, constants=None, on_combined=False,
          log=None):
    """results.tsv, results.json and a printable table. TSV because that is
    what goes into a spreadsheet without anyone being asked to parse JSON.

    `on_combined` orders results.tsv on the combined metric. Either way both
    orders are in the folder: the combined score and its components are
    columns of results.tsv, and results_combined.tsv is written beside it in
    combined order. A person reading one file is never without the other, and
    nothing that was reported from an earlier run changes meaning because a
    setting was switched.
    """
    os.makedirs(outdir, exist_ok=True)
    rows = list(rows)
    info = combined_scores(rows)
    if log and info is not None:
        from ._vendor import combined as _c
        _c.say(info, log)
    rows, key = rank(rows, top_n=top_n, on_combined=on_combined)
    cols = ["rank", "name", "peptide", "length", "poses", "adcp_best",
            "adcp_top%d" % top_n, "adcp_sd"]
    if any(r.get("combined_z") is not None for r in rows):
        # The components, not only the total: a construct that ranked third
        # did so for a reason, and the reason is these columns.
        cols += ["rank_combined", "combined_z"]
        cols += sorted({c for r in rows for c in r if c.startswith("z_")})
    if any(r.get("mmgbsa_n") for r in rows):
        cols += ["mmgbsa_n", "mmgbsa_best", "mmgbsa_top%d" % top_n,
                 "mmgbsa_median", "mmgbsa_sd", "mmgbsa_spread"]
    # One block per engine that has values anywhere in the run. Shown only
    # when more than one engine is present: with one, these repeat the
    # mmgbsa columns above under another name.
    present = [e for e in GB_ENGINES
               if any(r.get("mmgbsa_%s_n" % e) for r in rows)]
    if len(present) > 1:
        for e in present:
            cols += ["mmgbsa_%s_n" % e, "mmgbsa_%s_best" % e,
                     "mmgbsa_%s_top%d" % (e, top_n)]
    if any(r.get("mmgbsa_engine") for r in rows):
        cols += ["mmgbsa_engine"]
    if any(r.get("omm_n") for r in rows):
        cols += ["omm_n", "omm_best", "omm_best_pose"]
    if any(r.get("cluster_top") is not None for r in rows):
        cols += ["cluster_top", "cluster_frac", "cluster_modes"]
    if any(r.get("panel_pose") is not None for r in rows):
        cols += [c for c in PANEL_COLS
                 if any(r.get(c) not in (None, "", []) for r in rows)]
    if any(r.get("engine") for r in rows):
        cols += ["engine"]
    tsv = os.path.join(outdir, "results.tsv")
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(_tsv(r.get(c)) for c in cols) + "\n")
    # The other order, as a file of its own, so switching the setting is
    # never the only way to see it and a run carries both answers.
    if any(r.get("combined_z") is not None for r in rows):
        other = sorted([r for r in rows if r.get("combined_z") is not None],
                       key=lambda r: -r["combined_z"])
        other += [r for r in rows if r.get("combined_z") is None]
        with open(os.path.join(outdir, "results_combined.tsv"), "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in other:
                fh.write("\t".join(_tsv(r.get(c)) for c in cols) + "\n")
    doc = {"ranked_by": key, "top_n": top_n, "rows": rows}
    if constants:
        # What every panel number was measured with, so a reader of the file
        # alone can tell a probe radius from a cutoff.
        doc["constants"] = constants
    with open(os.path.join(outdir, "results.json"), "w") as fh:
        json.dump(doc, fh, indent=1)
    with open(os.path.join(outdir, "results.txt"), "w") as fh:
        fh.write(table(rows, top_n=top_n, on_combined=on_combined) + "\n")
    return tsv
