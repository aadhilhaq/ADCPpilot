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
              engine=""):
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
    return row


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


def write(rows, outdir, top_n=TOP_N):
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
    if any(r.get("engine") for r in rows):
        cols += ["engine"]
    tsv = os.path.join(outdir, "results.tsv")
    with open(tsv, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c))
                               for c in cols) + "\n")
    with open(os.path.join(outdir, "results.json"), "w") as fh:
        json.dump({"ranked_by": key, "top_n": top_n, "rows": rows}, fh, indent=1)
    with open(os.path.join(outdir, "results.txt"), "w") as fh:
        fh.write(table(rows, top_n=top_n) + "\n")
    return tsv
