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


def summarise(peptide, modes, mm=None, top_n=TOP_N, name="", openmm=None,
              engine="", panel=None, clusters=None):
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


def rank(rows, by=None, top_n=TOP_N):
    """Most negative first. MM-GBSA decides when it was run, because that is
    the more considered estimate; ADCP's own score when it was not."""
    key = by or ("mmgbsa_top%d" % top_n
                 if any(r.get("mmgbsa_top%d" % top_n) is not None for r in rows)
                 else "adcp_top%d" % top_n)
    have = [r for r in rows if r.get(key) is not None]
    gone = [r for r in rows if r.get(key) is None]
    have.sort(key=lambda r: r[key])
    for i, r in enumerate(have, 1):
        r["rank"] = i
    for r in gone:
        r["rank"] = None
    return have + gone, key


def table(rows, top_n=TOP_N):
    """Fixed-width, for a terminal or a log."""
    rows, key = rank(rows, top_n=top_n)
    mm = any(r.get("mmgbsa_n") for r in rows)
    omm = any(r.get("omm_n") for r in rows)
    head = ["#", "name", "peptide", "len", "poses", "ADCP best",
            "ADCP top%d" % top_n, "sd"]
    if mm:
        head += ["MM n", "MM best", "MM top%d" % top_n, "MM med", "sd",
                 "spread"]
    if omm:
        head += ["OpenMM best"]
    # Column widths are set from the content: a fixed name column either
    # truncates long names or leaves a gap on runs where nothing was
    # named.
    w = [max(4, len(head[0])),
         max([len(str(r.get("name") or "-")) for r in rows] + [len(head[1])])
         + 1 if rows else 10,
         (max(len(str(r["peptide"])) for r in rows) + 1) if rows else 10]
    out = ["ranked by %s, most negative first" % key, ""]
    fmt = ("%-4s %-*s %-*s %4s %5s %10s %11s %6s"
           + ("  %5s %8s %10s %8s %6s %7s" if mm else "")
           + ("  %11s" if omm else ""))

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
        out.append(fmt % tuple(vals))
    if any(r.get("panel_pose") is not None for r in rows):
        out += [""] + panel_table(rows)
    unranked = [r for r in rows if r.get("rank") is None]
    if unranked:
        out.append("")
        out.append("Unranked (no %s): %s. Peptides with non-standard or D "
                   "residues have no MM-GBSA value." % (key, ", ".join(
                       str(r.get("name") or r["peptide"]) for r in unranked)))
    if mm:
        widest = max((r.get("mmgbsa_spread") or 0) for r in rows)
        out.append("")
        out.append("MM-GBSA is an ordering within this run against this "
                   "receptor.")
        out.append("It is not a Kd and not a binding free energy: there is "
                   "no entropy term")
        out.append("in it. The whole picomolar-to-nanomolar range is 4.1 "
                   "kcal/mol; the widest")
        out.append("pose spread")
        if widest:
            out[-1] += (" in this run is %.1f. `adcp-dock affinity` prints "
                        "the ruler." % widest)
        else:
            out[-1] += (" is what the spread column would show. "
                        "`adcp-dock affinity` prints the ruler.")
    return "\n".join(out)


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
    odd = [r for r in rows if (r.get("panel_geometry") or "minimised complex")
           != "minimised complex"]
    if odd:
        out.append("")
        out.append("Measured on something other than a minimised complex: %s. "
                   "The geometry column says which." % ", ".join(
                       str(r.get("name") or r["peptide"]) for r in odd[:6]))
    if any(r.get("site_dg") is not None for r in rows):
        out.append("")
        out.append("A site energy is MMPBSA.py's decomposition of that pose, "
                   "and sums to within a per cent or two of its total rather "
                   "than exactly: a decomposition uses a per atom surface "
                   "term where the total uses LCPO. The residue map is in "
                   "metrics/.")
    if any(r.get("min_stopped_at_maxcyc") for r in rows):
        out.append("")
        out.append("The minimisation before the energy stopped at its cycle "
                   "limit rather than at a gradient, which is usual for a "
                   "restrained rescoring minimisation. The final gradient is "
                   "in min_rms_grad in results.tsv.")
    return out


def write(rows, outdir, top_n=TOP_N, constants=None):
    """results.tsv, results.json and a printable table. TSV because that is
    what goes into a spreadsheet without anyone being asked to parse JSON."""
    os.makedirs(outdir, exist_ok=True)
    rows, key = rank(list(rows), top_n=top_n)
    cols = ["rank", "name", "peptide", "length", "poses", "adcp_best",
            "adcp_top%d" % top_n, "adcp_sd"]
    if any(r.get("mmgbsa_n") for r in rows):
        cols += ["mmgbsa_n", "mmgbsa_best", "mmgbsa_top%d" % top_n,
                 "mmgbsa_median", "mmgbsa_sd", "mmgbsa_spread"]
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
    doc = {"ranked_by": key, "top_n": top_n, "rows": rows}
    if constants:
        # What every panel number was measured with, so a reader of the file
        # alone can tell a probe radius from a cutoff.
        doc["constants"] = constants
    with open(os.path.join(outdir, "results.json"), "w") as fh:
        json.dump(doc, fh, indent=1)
    with open(os.path.join(outdir, "results.txt"), "w") as fh:
        fh.write(table(rows, top_n=top_n) + "\n")
    return tsv
