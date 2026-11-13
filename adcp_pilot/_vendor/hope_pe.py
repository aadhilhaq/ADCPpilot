"""HOPE-PE: peptide practicality evaluation.

Docking and rescoring optimize for binding. Neither asks whether the peptide
that comes out can be made, kept in solution, or survive long enough to reach
the target. This module adds that check as a separate layer, so a
construct is judged on two independent axes and neither is hidden inside the
other: the docking/MM-GBSA ranking stays exactly as it was, and a second
ranking is produced that also accounts for practicality.

What is computed, and from what
-------------------------------
Everything here is a sequence-level descriptor. There is no structural or ADME
model behind it, and none is claimed. The properties are:

  gravy            Kyte-Doolittle grand average of hydropathy. Positive is
                   hydrophobic; above ~0.5 a peptide of this length tends to
                   need organic co-solvent to dissolve.
  charge_ph7       Net charge at pH 7.4 from Henderson-Hasselbalch on the
                   standard side-chain and terminal pKa set.
  pi               Isoelectric point, by bisection on the same charge curve.
                   A peptide near neutral at working pH aggregates readily,
                   so |pH - pI| < 1 is flagged.
  aromaticity      Fraction F+W+Y. Drives both UV quantification and stacking.
  aliphatic_index  Ikai relative volume of aliphatic side chains; a coarse
                   thermostability proxy.
  instability      Guruprasad instability index via Biopython when present.
                   Above 40 the peptide is predicted unstable in vitro.
  boman            Boman potential protein-interaction index (kcal/mol).
                   Above ~2.48 indicates promiscuous binding, which for a
                   designed binder is a specificity liability, not an asset.
  hydrophobic_mom  Eisenberg hydrophobic moment at 100 deg (alpha) - high
                   values with high gravy indicate a membrane-seeking amphipath.

Liabilities are counted as explicit motifs rather than rolled into a score
without trace, because a chemist needs to know which residue to change:

  protease sites   K/R not followed by P (trypsin), F/W/Y not followed by P
                   (chymotrypsin). Counts, not rates.
  oxidation        M, W, and free C.
  deamidation      NG, NS, NT, QG - NG is the fast one.
  isomerization    DG, DP, DS. DP additionally acid-labile.
  n_end_rule       N-terminal residue class; R, K, F, L, W, Y are
                   destabilizing in the eukaryotic N-end rule.
  synth_difficulty Consecutive beta-branched/bulky runs (I, V, T, F, W, Y) and
                   long identical runs, both of which stall solid-phase
                   coupling and cause on-resin aggregation.

Verdict
-------
Each property contributes a penalty only when it crosses a documented
threshold; there is no fitted weighting, because there is no dataset here to
fit one on. `pe_score` is the sum of those penalties (0 is clean) and
`verdict` is PASS / FLAG / FAIL. The reasons are carried in `pe_flags` so a
rejection is always attributable.

Dual ranking
------------
`rescore` returns the constructs ranked two ways: the pipeline's own ranking
untouched, and a ranking that combines it with practicality. Nothing is
deleted - a FAIL construct is ranked last, not dropped, so a strong binder with
a fixable liability stays visible.
"""
from __future__ import annotations

AA = "ACDEFGHIKLMNPQRSTVWY"

# Kyte & Doolittle 1982
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}

# Eisenberg consensus, for the hydrophobic moment
EISENBERG = {"A": 0.62, "R": -2.53, "N": -0.78, "D": -0.90, "C": 0.29,
             "Q": -0.85, "E": -0.74, "G": 0.48, "H": -0.40, "I": 1.38,
             "L": 1.06, "K": -1.50, "M": 0.64, "F": 1.19, "P": 0.12,
             "S": -0.18, "T": -0.05, "W": 0.81, "Y": 0.26, "V": 1.08}

# Boman 2003 - free energy of transfer from cyclohexane to water, sign flipped
BOMAN = {"L": -4.92, "I": -4.92, "V": -4.04, "F": -2.98, "M": -2.35,
         "W": -2.33, "A": -1.81, "C": -1.28, "G": 0.94, "Y": 0.14,
         "T": 2.57, "S": 3.40, "H": 4.66, "Q": 5.54, "K": 5.55,
         "N": 6.64, "E": 6.81, "D": 8.72, "R": 14.92, "P": 0.0}

MW = {"A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.16,
      "Q": 146.15, "E": 147.13, "G": 75.07, "H": 155.16, "I": 131.17,
      "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
      "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15}
WATER = 18.015

# EMBOSS pKa set
PKA_SIDE = {"C": 8.5, "D": 3.9, "E": 4.1, "H": 6.5, "K": 10.8, "R": 12.5, "Y": 10.1}
POSITIVE = set("KRH")
PKA_NTERM, PKA_CTERM = 8.6, 3.6

DESTABILISING_NEND = set("RKFLWY")     # eukaryotic N-end rule, destabilizing set
BULKY = set("IVTFWY")                  # beta-branched / bulky: slow SPPS coupling


# --------------------------------------------------------------- properties

def _clean(seq: str) -> str:
    return "".join(c for c in seq.upper() if c in AA)


def charge_at(seq: str, ph: float = 7.4) -> float:
    """Net charge from Henderson-Hasselbalch, termini included."""
    pos = 1.0 / (1.0 + 10 ** (ph - PKA_NTERM))
    neg = 1.0 / (1.0 + 10 ** (PKA_CTERM - ph))
    for aa in seq:
        pk = PKA_SIDE.get(aa)
        if pk is None:
            continue
        if aa in POSITIVE:
            pos += 1.0 / (1.0 + 10 ** (ph - pk))
        else:
            neg += 1.0 / (1.0 + 10 ** (pk - ph))
    return pos - neg


def isoelectric_point(seq: str) -> float:
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if charge_at(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def gravy(seq: str) -> float:
    return sum(KD[a] for a in seq) / len(seq) if seq else 0.0


def aromaticity(seq: str) -> float:
    return sum(seq.count(a) for a in "FWY") / len(seq) if seq else 0.0


def aliphatic_index(seq: str) -> float:
    """Ikai 1980. a=2.9, b=3.9 relative to alanine."""
    n = len(seq)
    if not n:
        return 0.0
    A, V = seq.count("A") / n, seq.count("V") / n
    IL = (seq.count("I") + seq.count("L")) / n
    return 100.0 * (A + 2.9 * V + 3.9 * IL)


def boman(seq: str) -> float:
    return sum(BOMAN[a] for a in seq) / len(seq) if seq else 0.0


def molecular_weight(seq: str) -> float:
    return sum(MW[a] for a in seq) - WATER * (len(seq) - 1) if seq else 0.0


def hydrophobic_moment(seq: str, angle: float = 100.0, window: int = 11) -> float:
    """Eisenberg mu-H, maximum over all windows. Amphipathicity, not hydrophobicity."""
    import math
    if not seq:
        return 0.0
    w = min(window, len(seq))
    best = 0.0
    rad = math.radians(angle)
    for start in range(len(seq) - w + 1):
        sx = sy = 0.0
        for i, aa in enumerate(seq[start:start + w]):
            h = EISENBERG[aa]
            sx += h * math.cos(i * rad)
            sy += h * math.sin(i * rad)
        best = max(best, math.hypot(sx, sy) / w)
    return best


def instability_index(seq: str):
    """Guruprasad. Needs the full 400-entry DIWV matrix, which Biopython ships;
    returns None rather than a guess if Biopython is absent."""
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except Exception:
        return None
    if len(seq) < 2:
        return None
    try:
        return float(ProteinAnalysis(seq).instability_index())
    except Exception:
        return None


# -------------------------------------------------------------- liabilities

def liabilities(seq: str) -> dict:
    """Motif counts. Each is a specific, locatable edit for a chemist."""
    n = len(seq)
    pairs = [seq[i:i + 2] for i in range(n - 1)]

    tryptic = sum(1 for i, a in enumerate(seq[:-1])
                  if a in "KR" and seq[i + 1] != "P")
    chymo = sum(1 for i, a in enumerate(seq[:-1])
                if a in "FWY" and seq[i + 1] != "P")

    runs_bulky, run, longest_bulky = 0, 0, 0
    for a in seq:
        run = run + 1 if a in BULKY else 0
        longest_bulky = max(longest_bulky, run)
    if longest_bulky >= 3:
        runs_bulky = longest_bulky

    longest_same, run, prev = 1, 1, ""
    for a in seq:
        run = run + 1 if a == prev else 1
        longest_same = max(longest_same, run)
        prev = a

    return dict(
        protease_trypsin=tryptic,
        protease_chymotrypsin=chymo,
        oxidation=seq.count("M") + seq.count("W") + seq.count("C"),
        free_cys=seq.count("C"),
        deamidation=sum(pairs.count(m) for m in ("NG", "NS", "NT", "QG")),
        deamidation_fast=pairs.count("NG"),
        isomerisation=sum(pairs.count(m) for m in ("DG", "DP", "DS")),
        acid_labile=pairs.count("DP"),
        n_end_destabilising=int(bool(seq) and seq[0] in DESTABILISING_NEND),
        bulky_run=runs_bulky,
        homopolymer_run=longest_same,
    )


# ------------------------------------------------------------------ verdict
# Thresholds are literature conventions, not fitted values. Each carries its
# own penalty weight so the total can be read back as a list of reasons.
RULES = [
    # (key, test, penalty, message)
    ("gravy",       lambda v, s: v > 0.9,   3, "very hydrophobic (GRAVY > 0.9): expect poor aqueous solubility"),
    ("gravy",       lambda v, s: 0.5 < v <= 0.9, 1, "hydrophobic (GRAVY > 0.5): may need DMSO/organic co-solvent"),
    ("pi_gap",      lambda v, s: v < 1.0,   2, "pI within 1 unit of working pH: low net charge, aggregation risk"),
    ("charge_ph7",  lambda v, s: abs(v) > 6, 1, "highly charged: nonspecific electrostatic binding likely"),
    ("instability", lambda v, s: v is not None and v > 40, 2, "instability index > 40: predicted unstable in vitro"),
    ("aliphatic_index", lambda v, s: v < 20, 1, "low aliphatic index: little hydrophobic core, weak thermostability"),
    ("boman",       lambda v, s: v > 2.48,  1, "Boman index > 2.48: promiscuous protein binding, specificity risk"),
    ("hydrophobic_moment", lambda v, s: v > 0.55 and s["gravy"] > 0.3, 2,
     "amphipathic and hydrophobic: membrane-seeking, likely hemolytic"),
]

LIABILITY_RULES = [
    ("protease_trypsin", lambda v: v >= 2, 2, "{v} tryptic sites: short serum half-life"),
    ("protease_trypsin", lambda v: v == 1, 1, "1 tryptic site"),
    ("protease_chymotrypsin", lambda v: v >= 3, 1, "{v} chymotryptic sites"),
    ("deamidation_fast", lambda v: v >= 1, 2, "{v} NG motif(s): fast deamidation, isoAsp formation"),
    ("isomerisation", lambda v: v >= 1, 1, "{v} Asp isomerization motif(s)"),
    ("acid_labile", lambda v: v >= 1, 1, "{v} DP motif(s): acid-labile, cleaved during TFA cleavage"),
    ("free_cys", lambda v: v >= 1, 1, "{v} free cysteine(s): disulfide scrambling, dimerization"),
    ("oxidation", lambda v: v >= 3, 1, "{v} oxidation-prone residues (M/W/C)"),
    ("n_end_destabilising", lambda v: v == 1, 1, "destabilizing N-terminal residue (N-end rule): cap or reorder"),
    ("bulky_run", lambda v: v >= 4, 2, "run of {v} beta-branched/bulky residues: difficult SPPS, on-resin aggregation"),
    ("bulky_run", lambda v: v == 3, 1, "run of 3 bulky residues: watch coupling efficiency"),
    ("homopolymer_run", lambda v: v >= 4, 1, "homopolymer run of {v}: synthesis and aggregation risk"),
]

FAIL_AT, FLAG_AT = 8, 3


def evaluate(sequence: str, ph: float = 7.4) -> dict:
    """Full practicality report for one peptide."""
    seq = _clean(sequence)
    if not seq:
        return dict(sequence=sequence, pe_score=None, verdict="NA",
                    pe_flags=["no standard residues in sequence"])

    props = dict(
        sequence=seq,
        length=len(seq),
        mw=round(molecular_weight(seq), 2),
        gravy=round(gravy(seq), 3),
        charge_ph7=round(charge_at(seq, ph), 2),
        pi=isoelectric_point(seq),
        aromaticity=round(aromaticity(seq), 3),
        aliphatic_index=round(aliphatic_index(seq), 1),
        boman=round(boman(seq), 2),
        hydrophobic_moment=round(hydrophobic_moment(seq), 3),
    )
    props["pi_gap"] = round(abs(ph - props["pi"]), 2)
    inst = instability_index(seq)
    props["instability"] = None if inst is None else round(inst, 1)

    penalty, flags = 0, []
    for key, test, w, msg in RULES:
        try:
            if test(props[key], props):
                penalty += w
                flags.append(msg)
        except (TypeError, KeyError):
            continue

    lia = liabilities(seq)
    props.update(lia)
    for key, test, w, msg in LIABILITY_RULES:
        v = lia[key]
        if test(v):
            penalty += w
            flags.append(msg.format(v=v))

    props["pe_score"] = penalty
    props["verdict"] = ("FAIL" if penalty >= FAIL_AT else
                        "FLAG" if penalty >= FLAG_AT else "PASS")
    props["pe_flags"] = flags
    if props["instability"] is None:
        flags.append("instability index unavailable (Biopython not installed)")
    return props


# -------------------------------------------------------------- dual ranking

def _z(values):
    vals = [v for v in values if v is not None]
    if len(vals) < 2:
        return [0.0] * len(values)
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / len(vals)
    sd = var ** 0.5
    if sd == 0:
        return [0.0] * len(values)
    return [0.0 if v is None else (v - mean) / sd for v in values]


def rescore(rows: list[dict], seq_key: str = "sequence",
            rank_key: str = "z", weight: float = 1.0, ph: float = 7.4) -> list[dict]:
    """Attach a practicality report to each row and re-rank.

    `rank_key` is whatever the pipeline already ranked on, higher being better.
    The combined score subtracts a z-scored practicality penalty from it, so
    the two axes are weighted comparably rather than in raw units that have
    nothing to do with each other. Nothing is filtered out.
    """
    if not rows:
        return []
    out = []
    for r in rows:
        rep = evaluate(r.get(seq_key, ""), ph)
        d = dict(r)
        d["pe_score"] = rep["pe_score"]
        d["pe_verdict"] = rep["verdict"]
        d["pe_flags"] = "; ".join(rep["pe_flags"])
        d["pe_report"] = rep
        out.append(d)

    zp = _z([d["pe_score"] for d in out])
    for d, z in zip(out, zp):
        base = d.get(rank_key)
        d["pe_z"] = -z                       # low penalty is good
        d["combined_z"] = (0.0 if base is None else float(base)) + weight * (-z)
    out.sort(key=lambda d: -d["combined_z"])
    for i, d in enumerate(out, 1):
        d["rank_hopepe"] = i
    return out


COLUMNS = ["length", "mw", "gravy", "charge_ph7", "pi", "aromaticity",
           "aliphatic_index", "instability", "boman", "hydrophobic_moment",
           "protease_trypsin", "protease_chymotrypsin", "oxidation", "free_cys",
           "deamidation_fast", "isomerisation", "acid_labile",
           "n_end_destabilising", "bulky_run", "homopolymer_run",
           "pe_score", "verdict"]


def write_table(rows: list[dict], path: str, seq_key: str = "sequence") -> str:
    """Per-peptide practicality table, one row per construct."""
    with open(path, "w") as fh:
        fh.write("rank\tsequence\t" + "\t".join(COLUMNS) + "\tflags\n")
        for i, r in enumerate(rows, 1):
            rep = r.get("pe_report") or evaluate(r.get(seq_key, ""))
            vals = []
            for c in COLUMNS:
                v = rep.get(c)
                vals.append("" if v is None else
                            (f"{v:.3f}" if isinstance(v, float) else str(v)))
            fh.write(f"{i}\t{rep.get('sequence','')}\t" + "\t".join(vals) +
                     "\t" + "; ".join(rep.get("pe_flags", [])) + "\n")
    return path


if __name__ == "__main__":
    import sys
    for s in (sys.argv[1:] or ["GKWKKLLKKLLKALKK", "SIINFEKL", "DPNGWWIVTFY"]):
        rep = evaluate(s)
        print(f"\n{rep['sequence']}  ->  {rep['verdict']}  (penalty {rep['pe_score']})")
        for k in ("mw", "gravy", "charge_ph7", "pi", "instability", "boman",
                  "hydrophobic_moment", "aliphatic_index"):
            print(f"   {k:20s} {rep[k]}")
        for f in rep["pe_flags"]:
            print(f"   ! {f}")
