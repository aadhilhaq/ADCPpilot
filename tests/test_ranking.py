#!/usr/bin/env python3
"""Which order the results table comes out in, and what is in it either way.

By Aadhil Haq

The table can be ordered on the single column it has always used or on the
combined metric. The thing worth pinning is not which order is better - that
cannot be shown from anything on this cluster - but that switching the flag
never costs the reader the other answer, and never changes a number. Both
orders are in the run folder whichever way the flag is set, and the rows
carry the same values in both.

The metric itself is tested where it lives, in HOPE-pipelines'
tests/test_combined.py. This covers only what this repository adds.

Exit 0 passes, 77 skips, anything else fails.
"""
import csv
import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import report                                  # noqa: E402

FAIL = []


def check(what, ok, detail=""):
    print(("  ok   " if ok else "  FAIL ") + what
          + (("  " + detail) if detail and not ok else ""))
    if not ok:
        FAIL.append(what)


def rows():
    """Five peptides whose two orderings differ.

    MM-GBSA runs one way and the interface the other, which is the case the
    flag exists for: where they agree the flag cannot be seen to do anything.
    """
    out = []
    for i, (gb, adcp, bsa_pr, hb, dg) in enumerate([
            (-70.0, -24.0, 80.0, 4, -2.0),
            (-60.0, -22.0, 120.0, 18, -14.0),
            (-50.0, -20.0, 110.0, 14, -11.0),
            (-40.0, -18.0, 100.0, 10, -8.0),
            (-30.0, -16.0, 90.0, 6, -5.0)]):
        out.append({"peptide": "PEP%d" % i, "name": "p%d" % i, "length": 7,
                    "poses": 5, "adcp_best": adcp - 1, "adcp_top5": adcp,
                    "adcp_sd": 0.5, "mmgbsa_n": 5, "mmgbsa_best": gb - 5,
                    "mmgbsa_top5": gb, "mmgbsa_median": gb, "mmgbsa_sd": 2.0,
                    "mmgbsa_spread": 10.0, "panel_pose": 1,
                    "bsa_per_res": bsa_pr, "n_hbond": hb, "n_ionic": 1,
                    "n_clash": 0, "site_dg": dg, "site_coverage": 1.0,
                    "bsa": bsa_pr * 7})
    return out


print("the two orderings")
single, key = report.rank(rows())
check("without the flag the table is ordered as it always was",
      key == "mmgbsa_top5"
      and [r["peptide"] for r in single] == ["PEP0", "PEP1", "PEP2", "PEP3",
                                             "PEP4"],
      "by %s: %s" % (key, [r["peptide"] for r in single]))

r = rows()
report.combined_scores(r)
comb, key = report.rank(r, on_combined=True)
check("with the flag the table is ordered on the combined score",
      key == "combined_z", "by %s" % key)
check("the combined order is highest first, not most negative first",
      all(comb[i]["combined_z"] >= comb[i + 1]["combined_z"]
          for i in range(len(comb) - 1)),
      "%s" % [r["combined_z"] for r in comb])
check("the two orderings actually differ on these rows",
      [x["peptide"] for x in comb] != [x["peptide"] for x in single],
      "%s" % [x["peptide"] for x in comb])

# The flag must not be able to change a measurement, only the order.
by_pep = {x["peptide"]: x for x in single}
check("no measured value changes with the flag",
      all(by_pep[x["peptide"]]["mmgbsa_top5"] == x["mmgbsa_top5"]
          and by_pep[x["peptide"]]["bsa_per_res"] == x["bsa_per_res"]
          for x in comb))

print("\nwhat reaches the files")
for flag in (False, True):
    with tempfile.TemporaryDirectory() as d:
        report.write(rows(), d, on_combined=flag)
        have = sorted(os.listdir(d))
        check("rank_combined=%s writes both tables" % flag,
              "results.tsv" in have and "results_combined.tsv" in have,
              "%s" % have)
        main = list(csv.DictReader(open(os.path.join(d, "results.tsv"),
                                        newline=""), delimiter="\t"))
        other = list(csv.DictReader(
            open(os.path.join(d, "results_combined.tsv"), newline=""),
            delimiter="\t"))
        check("rank_combined=%s: the combined score is a column either way"
              % flag,
              all(x.get("combined_z") for x in main))
        check("rank_combined=%s: the components are columns, not just the "
              "total" % flag,
              any(c.startswith("z_") for c in main[0]))
        check("rank_combined=%s: the two files hold the same peptides" % flag,
              sorted(x["peptide"] for x in main)
              == sorted(x["peptide"] for x in other))
        order_main = [x["peptide"] for x in main]
        order_other = [x["peptide"] for x in other]
        check("rank_combined=%s: results_combined.tsv is in combined order"
              % flag,
              order_other == [x["peptide"] for x in comb],
              "%s" % order_other)
        if flag:
            check("rank_combined=True: results.tsv is in combined order too",
                  order_main == order_other)
        else:
            check("rank_combined=False: results.tsv keeps the old order",
                  order_main == [x["peptide"] for x in single])
        txt = open(os.path.join(d, "results.txt")).read()
        check("rank_combined=%s: the printed table says what it is ranked by"
              % flag,
              ("ranked by combined_z" if flag
               else "ranked by mmgbsa_top5") in txt,
              txt.splitlines()[0] if txt else "")

print("\na run with nothing to combine")
bare = [{"peptide": "P%d" % i, "name": "", "length": 7, "poses": 3,
         "adcp_best": -20.0 + i, "adcp_top5": -19.0 + i, "adcp_sd": 0.3}
        for i in range(4)]
with tempfile.TemporaryDirectory() as d:
    report.write(bare, d, on_combined=True)
    have = sorted(os.listdir(d))
    check("a run with no panel writes no combined table",
          "results_combined.tsv" not in have, "%s" % have)
    got = list(csv.DictReader(open(os.path.join(d, "results.tsv"),
                                   newline=""), delimiter="\t"))
    check("and falls back to the single column rather than failing",
          [x["peptide"] for x in got] == ["P0", "P1", "P2", "P3"],
          "%s" % [x["peptide"] for x in got])

with tempfile.TemporaryDirectory() as d:
    report.write([], d, on_combined=True)
    check("an empty run still writes a table",
          os.path.isfile(os.path.join(d, "results.tsv")))

print()
if FAIL:
    print("%d FAILED:" % len(FAIL))
    for f in FAIL:
        print("  " + f)
    raise SystemExit(1)
print("all checks passed")
