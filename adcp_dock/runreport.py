"""The report of a run: what it says, for the lab's report page.

By Aadhil Haq

The report page (GET /report) is laid out by reportkit, the module the lab's
tools share (static/reportkit.js and reportkit.css); this decides what an ADCP
docking run's report says and sends it as the document reportkit reads
(GET /api/report). In order:

  Summary              which peptide ranks first, by how much, and whether the
                       run separates it from the ones after it
  The receptor and     how the receptor was prepared, or where it came from,
  the site             and the site and box the maps were built over
  The peptides ranked  the ranking, a chart of it, and the interface panel
  Each peptide         its poses and their scores, pictures of the top ones,
                       the interaction diagram and residue map of the pose the
                       panel describes, its ring, its sequence flags
  Views, Methods,      the views added from the viewer, a methods section
  Appendix             written from what the run did, with its references, and
                       the preparation log, the maps' log and every setting

Every caption states the result from the run's own numbers, then only what is
needed to read it. Every number is one the run wrote, or one the results page
works out the same way: a pose's contacts are classified by the rules the
viewer draws them with (static/app.js, classifyContacts), on the same complex
(server.structure_text).
"""
import csv
import html
import json
import os
import re
import statistics
import time

from . import job as jobmod

# ---- the contact rules, the viewer's own (app.js classifyContacts) --------------
HBOND_A, IONIC_A, STACK_A, VDW_A = 3.5, 4.0, 5.0, 4.0
POLAR = {"N", "O", "S"}
BASIC = {"ARG": ("NE", "NH1", "NH2"), "LYS": ("NZ",), "HIS": ("ND1", "NE2")}
ACIDIC = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
RING = {"PHE", "TYR", "TRP", "HIS"}
#: Most specific first: a residue pair is drawn as the most specific kind any of
#: its atom pairs makes.
SPECIFIC = ("ionic", "stack", "hbond", "vdw")
KIND_COLOUR = {"hbond": "#2B83BA", "ionic": "#F0C814", "stack": "#8CB366", "vdw": "#808080"}
KIND_WORD = {"hbond": "hydrogen bond", "ionic": "salt bridge", "stack": "stacking",
             "vdw": "van der Waals"}

#: Amber's names for a residue's states, read as the residue itself.
AMBER_NAMES = {"HIE": "His", "HID": "His", "HIP": "His", "HSD": "His", "HSE": "His",
               "HSP": "His", "CYX": "Cys", "CYM": "Cys", "ASH": "Asp", "GLH": "Glu",
               "LYN": "Lys"}
GB_LABELS = {"amber": "AmberTools", "openmm": "OpenMM", "schrodinger": "Prime"}
SMALL = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"]

#: Where the views added from the viewer are kept, inside the run.
FOLDER = "report"
IMAGE_MAX = 12 << 20


# ---- writing numbers and words --------------------------------------------------
def e(s):
    return html.escape(str(s if s is not None else ""), quote=True)


def num(v, d=1, sign=False):
    """A number as the report prints it: fixed decimals, a true minus sign,
    thousands grouped from 10,000 up, as reportkit.js does."""
    if v is None:
        return ""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return e(v)
    s = "%.*f" % (d, v)
    if abs(v) >= 10000:
        a, _, b = s.lstrip("-").partition(".")
        s = ("-" if v < 0 else "") + "{:,}".format(int(a)) + ("." + b if b else "")
    if sign and v > 0:
        s = "+" + s
    return s.replace("-", "\u2212")


def pm(m, sd, d=1):
    if m is None:
        return ""
    return num(m, d) + ("" if sd is None else " \u00b1 " + num(sd, d))


def count(n, word, plural=None):
    """'one pose', '12 poses': small numbers as words, as prose has them."""
    w = SMALL[n] if 0 <= n < len(SMALL) else str(n)
    return "%s %s" % (w, word if n == 1 else (plural or word + "s"))


def words(n):
    return SMALL[n] if 0 <= n < len(SMALL) else str(n)


def and_list(items):
    items = [str(i) for i in items if i not in (None, "")]
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def res_name(resn, n):
    """'TRP', 4 -> 'Trp4'; Amber's state names as the residue itself."""
    r = str(resn or "?").upper()
    r = AMBER_NAMES.get(r, r)
    return "%s%s" % (r[:1].upper() + r[1:].lower() if len(r) == 3 else r, n)


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


# ---- the contacts of a pose, as the viewer classifies them -----------------------
def _atoms(text):
    out = []
    for L in text.splitlines():
        if L[:4] != "ATOM" and L[:6] != "HETATM":
            continue
        el = (L[76:78].strip() or L[12:16].strip()[:1]).upper()
        if el == "H":
            continue
        try:
            out.append((L[12:16].strip(), L[17:20].strip(), L[21:22], int(L[22:26]), el,
                        float(L[30:38]), float(L[38:46]), float(L[46:54])))
        except ValueError:
            continue
    return out


def peptide_chain(text):
    """The peptide is the short chain, as app.js peptideChain reads it."""
    n = {}
    for L in text.splitlines():
        if L[:4] == "ATOM" or L[:6] == "HETATM":
            n[L[21:22]] = n.get(L[21:22], 0) + 1
    return min(n, key=lambda c: n[c]) if n else "Z"


def classify(text, pep=None):
    """{"pairs": one per residue pair, "receptor": the receptor residues in
    contact, "peptide": the peptide residues, "chain": the peptide's chain}.

    The rules and the bookkeeping are app.js classifyContacts', line for line:
    heavy atoms within 5.0 A; van der Waals within 4.0, a hydrogen bond between
    N, O or S within 3.5, a salt bridge between a basic and an acidic side-chain
    tip within 4.0, stacking between carbons of two ring residues within 5.0;
    one line per residue pair, at its closest atom pair, as its most specific
    kind."""
    import numpy as np
    pep = pep or peptide_chain(text)
    at = _atoms(text)
    rec = [a for a in at if a[2] != pep]
    lig = [a for a in at if a[2] == pep]
    out = {"pairs": [], "receptor": [], "peptide": [], "chain": pep}
    if not rec or not lig:
        return out
    R = np.array([a[5:8] for a in rec])
    Lg = np.array([a[5:8] for a in lig])
    d = np.sqrt(((R[:, None, :] - Lg[None, :, :]) ** 2).sum(-1))
    pairs, byres, seen = [], {}, {}
    for i, k in zip(*np.nonzero(d <= STACK_A)):
        a, b, dist = rec[i], lig[k], float(d[i, k])
        kinds = set()
        if dist <= VDW_A:
            kinds.add("vdw")
        if dist <= HBOND_A and a[4] in POLAR and b[4] in POLAR:
            kinds.add("hbond")
        if dist <= IONIC_A and ((a[0] in BASIC.get(a[1], ()) and b[0] in ACIDIC.get(b[1], ()))
                                or (a[0] in ACIDIC.get(a[1], ()) and b[0] in BASIC.get(b[1], ()))):
            kinds.add("ionic")
        if a[1] in RING and b[1] in RING and a[4] == "C" and b[4] == "C":
            kinds.add("stack")
        if not kinds:
            continue
        kind = next(s for s in SPECIFIC if s in kinds)
        key = (a[2], a[3], b[2], b[3])
        have = seen.get(key)
        if have:
            if dist < have["d"]:
                have["d"] = dist
            if SPECIFIC.index(kind) < SPECIFIC.index(have["kind"]):
                have["kind"] = kind
        else:
            seen[key] = {"chain": a[2], "resi": a[3], "resn": a[1], "bchain": b[2],
                         "bresi": b[3], "bresn": b[1], "d": dist, "kind": kind}
            pairs.append(seen[key])
        r = byres.setdefault((a[2], a[3]), {"chain": a[2], "resi": a[3], "resn": a[1],
                                            "kinds": set(), "d": dist})
        r["kinds"].add(kind)
        r["d"] = min(r["d"], dist)
    out["pairs"] = pairs
    out["receptor"] = sorted(byres.values(), key=lambda r: r["d"])
    out["peptide"] = sorted({(p["bresi"], p["bresn"]) for p in pairs})
    return out


def lid_svg(pairs, on=("ionic", "stack", "hbond"), band=8):
    """The interaction diagram, as the viewer draws it (app.js drawLid): the
    peptide along the middle, each receptor residue above or below the peptide
    residue it is closest to, lines coloured by kind, with a key of the kinds
    drawn. A long peptide is drawn in bands of `band` residues, one under
    another, so the labels stay legible on a page; a line to a residue in
    another band runs to the edge of its own. None when there is nothing to
    draw."""
    import math
    pairs = [p for p in pairs if p["kind"] in on]
    if not pairs:
        return None
    pep = {}
    for p in pairs:
        pep.setdefault((p["bchain"], p["bresi"]), {"resi": p["bresi"], "resn": p["bresn"]})
    peps = sorted(pep.items(), key=lambda kv: kv[1]["resi"])
    at = {k: i for i, (k, _v) in enumerate(peps)}
    rec = {}
    for p in pairs:
        r = rec.setdefault((p["chain"], p["resi"]), {"resi": p["resi"], "resn": p["resn"],
                                                     "d": 1e9, "at": 0, "links": []})
        r["links"].append(p)
        if p["d"] < r["d"]:
            r["d"] = p["d"]
            r["at"] = at.get((p["bchain"], p["bresi"]), 0)
    recs = sorted(rec.values(), key=lambda r: (r["at"], r["d"]))
    longest = max([6] + [len(res_name(r["resn"], r["resi"])) for r in recs]
                  + [len(res_name(v["resn"], v["resi"])) for _k, v in peps])
    RX, RY, LANE = max(32, int(math.ceil(longest * 4.2)) + 8), 14, 46
    BOXW = max(64, int(math.ceil(longest * 8.6)) + 12)
    SLOT = max(96, 2 * RX + 18, BOXW + 18)
    PAD = RX + 10

    def lane(lst):
        used, deepest = set(), -1
        for r in lst:
            L = 0
            while (L, r["at"]) in used:
                L += 1
            used.add((L, r["at"]))
            r["lane"] = L
            deepest = max(deepest, L)
        return deepest
    cols = min(band, len(peps))
    W = PAD * 2 + max(cols, 1) * SLOT
    lines, nodes = [], []
    y = 8
    for lo in range(0, len(peps), band):
        hi = min(len(peps), lo + band)
        mine = [r for r in recs if lo <= r["at"] < hi]
        above, below = mine[0::2], mine[1::2]
        da, db = lane(above), lane(below)
        midY = y + PAD + (da + 1) * LANE + RY

        def px(i, lo=lo):
            return PAD + (i - lo) * SLOT + SLOT / 2.0
        for r in above:
            r["x"], r["y"] = px(r["at"]), midY - 30 - (r["lane"] + 1) * LANE
        for r in below:
            r["x"], r["y"] = px(r["at"]), midY + 30 + (r["lane"] + 1) * LANE
        for r in mine:
            for p in r["links"]:
                i = at.get((p["bchain"], p["bresi"]))
                if i is None:
                    continue
                x2 = px(i) if lo <= i < hi else (PAD if i < lo else W - PAD)
                lines.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1.8" '
                             'stroke-dasharray="4 3"><title>%s to %s, %.2f \u00c5</title></line>'
                             % (r["x"], r["y"], x2, midY, KIND_COLOUR[p["kind"]],
                                e(res_name(r["resn"], r["resi"])), e(res_name(p["bresn"], p["bresi"])), p["d"]))
        nodes.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#D4763A" stroke-width="1" opacity="0.35"/>'
                     % (PAD, midY, PAD + (hi - lo) * SLOT, midY))
        for i in range(lo, hi):
            v = peps[i][1]
            x = px(i)
            nodes.append('<rect x="%.1f" y="%.1f" width="%d" height="32" rx="7" fill="#fdf2e9" stroke="#D4763A" '
                         'stroke-width="1.3"/><text x="%.1f" y="%.1f" text-anchor="middle" font-size="14" '
                         'font-weight="600" fill="#8a4620">%s</text>'
                         % (x - BOXW / 2.0, midY - 16, BOXW, x, midY + 5, e(res_name(v["resn"], v["resi"]))))
        for r in mine:
            kind = next(k for k in SPECIFIC if any(p["kind"] == k for p in r["links"]))
            nodes.append('<ellipse cx="%.1f" cy="%.1f" rx="%d" ry="%d" fill="#ffffff" stroke="%s" stroke-width="1.8"/>'
                         '<text x="%.1f" y="%.1f" text-anchor="middle" font-size="13" fill="#20262b">%s</text>'
                         % (r["x"], r["y"], RX, RY, KIND_COLOUR[kind], r["x"], r["y"] + 4.5,
                            e(res_name(r["resn"], r["resi"]))))
        y = midY + 30 + (db + 1) * LANE + RY + 18
    key, x = [], PAD
    for k in SPECIFIC:
        if not any(p["kind"] == k for p in pairs):
            continue
        key.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2" stroke-dasharray="4 3"/>'
                   '<text x="%.1f" y="%.1f" font-size="13" fill="#20262b">%s</text>'
                   % (x, y + 8, x + 26, y + 8, KIND_COLOUR[k], x + 32, y + 12.5, KIND_WORD[k]))
        x += 32 + 7.4 * len(KIND_WORD[k]) + 26
    H, Wk = y + 26, max(W, int(x) + 10)
    return ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" '
            'font-family="Helvetica, Arial, sans-serif" role="img" aria-label="interaction diagram">'
            '<rect width="100%%" height="100%%" fill="#ffffff"/>%s%s%s</svg>'
            % (Wk, H, Wk, H, "".join(lines), "".join(nodes), "".join(key)))


def pair_sentence(pairs, name, receptor, pose):
    """What a pose's diagram shows, in one sentence, from its pairs: the salt
    bridges and stacked pairs named, peptide residue first, and the
    hydrogen-bonded pairs counted."""
    ionic = [p for p in pairs if p["kind"] == "ionic"]
    stack = [p for p in pairs if p["kind"] == "stack"]
    hb = [p for p in pairs if p["kind"] == "hbond"]
    nrec = len({(p["chain"], p["resi"]) for p in pairs})

    def grouped(lst, verb, most=4):
        by = {}
        for p in sorted(lst, key=lambda p: p["d"]):
            by.setdefault((p["bresi"], p["bresn"]), []).append(res_name(p["resn"], p["resi"]))
        parts = ["%s %s %s" % (res_name(bn, bi), verb, and_list(rs)) for (bi, bn), rs in list(by.items())[:most]]
        return "; ".join(parts) + ("; and others" if len(by) > most else "")
    bits = []
    if ionic:
        bits.append("%s (%s)" % ("a salt bridge" if len(ionic) == 1 else count(len(ionic), "salt bridge"),
                                 grouped(ionic, "to")))
    if stack:
        bits.append("%s (%s)" % (count(len(stack), "stacked pair"), grouped(stack, "on")))
    if hb:
        bits.append(count(len(hb), "hydrogen-bonded pair"))
    if not bits:
        return ("In pose %s, %s touches %s through van der Waals contacts alone, with %s"
                % (pose, e(name), e(receptor), count(nrec, "receptor residue")))
    return ("In pose %s, %s holds %s by %s, through %s"
            % (pose, e(name), e(receptor), and_list(bits), count(nrec, "receptor residue")))


# ---- reading the run ------------------------------------------------------------
def _json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def receptor_name(j):
    """'6B88 chain A': the structure's code or file name, and its chains."""
    a = j.get("adapt") or {}
    base = a.get("pdbid") or os.path.basename(str(j.get("receptor") or "")) or "the receptor"
    base = re.sub(r"(\.(pdb|ent|cif|mmcif|pdbqt))?(\.gz)?$", "", base, flags=re.I)
    ch = [c for c in (j.get("chains") or ((j.get("site") or {}).get("chains")) or []) if c]
    if not ch:
        return base
    return "%s chain%s %s" % (base, "s" if len(ch) > 1 else "", and_list(ch))


def engine_label(j, rows):
    """Which ADCP docked the run, as the report names it."""
    from . import mac10
    if j.get("adcp_engine") == mac10.ENGINE or any(r.get("engine") == mac10.LABEL for r in rows):
        return mac10.LABEL
    if (j.get("converted") or {}).get("from") == "ADCP 1.0" or any(r.get("engine") == "ADCP 1.0" for r in rows):
        return "ADCP 1.0"
    got = next((r.get("engine") for r in rows if str(r.get("engine") or "").startswith("ADCP")), "")
    return got or "ADCP 1.1"


def _who(rundir):
    try:
        import pwd
        return pwd.getpwuid(os.stat(rundir).st_uid).pw_name
    except (ImportError, KeyError, OSError):
        return ""


def prep_steps(rundir):
    """[(step, outcome, detail)] from receptor/preparation.log, and its last
    line ('818 atoms in, 1449 out')."""
    path = os.path.join(rundir, "receptor", "preparation.log")
    steps, tail = [], ""
    try:
        lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    except OSError:
        return steps, tail, None
    for L in lines:
        m = re.match(r"^\s{2,}([a-z][a-z \-]+?):\s*(done|skipped|failed|not done)\b\s*(?:[;\-]\s*)?(.*)$", L)
        if m:
            steps.append((m.group(1), m.group(2), m.group(3).strip()))
        elif re.search(r"\d+ atoms in, \d+ out", L):
            tail = L.strip()
    return steps, tail, "\n".join(lines)


def prep_facts(steps, tail, j):
    """The facts of a preparation, in a few words each, for a caption."""
    by = {s: d for s, o, d in steps if o == "done"}
    out = []
    m = re.search(r"(\d+) atoms in, (\d+) out", tail or "")
    if m:
        out.append("%s atoms in, %s out" % m.groups())
    ch = by.get("chain", "")
    m = re.search(r"dropped ([A-Za-z0-9, ]+)", ch)
    if m:
        out.append("chain%s %s dropped" % ("s" if "," in m.group(1) else "", m.group(1).strip()))
    m = re.search(r"removed (\d+)", by.get("waters", ""))
    if m and int(m.group(1)):
        out.append("%s waters removed" % m.group(1))
    het = by.get("heteroatoms", "")
    if het and not het.startswith("removed none"):
        out.append(het.replace("removed ", "") + " removed" if het.startswith("removed") else het)
    miss = by.get("missing atoms", "")
    m = re.search(r"built (\d+) missing residue", miss)
    if m and int(m.group(1)):
        out.append("%s missing residues built" % m.group(1))
    m = re.search(r"rebuilt (\d+) side-chain atom", miss)
    if m and int(m.group(1)):
        out.append("%s side-chain atoms rebuilt" % m.group(1))
    prot = by.get("protonation", "")
    m = re.search(r"(\d+) hydrogens added", prot)
    if prot:
        how = ("by reduce" if "reduce" in prot else
               "by PDB2PQR with PROPKA at pH %s" % j.get("ph", 7.4) if "PROPKA" in prot else "")
        out.append(("%s hydrogens added %s" % (m.group(1), how) if m else "protonated %s" % how).strip())
    m = re.search(r"moved ([\d.]+) A RMSD", by.get("minimisation", ""))
    if m:
        out.append("heavy atoms moved %s \u00c5 in the restrained minimisation" % m.group(1))
    return out


def residue_map(rundir, seq):
    """(peptide rows, receptor rows, the pose they describe) from metrics/."""
    key = jobmod.key(seq)
    path = os.path.join(rundir, "metrics", "%s_residues.tsv" % key)
    if not os.path.isfile(path):
        return None, None, None
    pep, rec = [], []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            r.setdefault("side", "receptor")
            for k in ("dg", "d_min"):
                r[k] = _f(r.get(k)) if r.get(k) not in ("", None) else None
            (pep if r["side"] == "peptide" else rec).append(r)
    pep.sort(key=lambda r: int(r.get("resnum") or 0))
    for i, r in enumerate(pep, 1):
        r["pos"] = i
    pose = (_json(os.path.join(rundir, "metrics", "%s.json" % key)) or {}).get("best_pose")
    return pep, rec, pose


# ---- the report -----------------------------------------------------------------
OPTIONS = [
    {"key": "peptides", "label": "Peptides", "type": "number", "least": 1, "most": 50,
     "help": "How many peptides, from the top of the ranking, get a section of their own"},
    {"key": "poses", "label": "Poses drawn", "type": "number", "default": 3, "least": 1, "most": 5,
     "help": "How many of each peptide's top poses are drawn"},
    {"key": "diagrams", "label": "Diagrams for", "type": "select", "default": "panel",
     "choices": [["panel", "the pose the panel describes"], ["all", "every pose drawn"]],
     "help": "Which poses get an interaction diagram"},
    {"key": "vdw", "label": "van der Waals in diagrams", "type": "check", "default": False,
     "help": "Draw van der Waals contacts in the diagrams as well as the specific ones"},
]


def options(raw, total):
    """The options as asked for, held to what they may be."""
    out = {}
    for o in OPTIONS:
        d = o.get("default", min(total, 10) or 1)
        v = raw.get(o["key"], d)
        if o["type"] == "number":
            try:
                v = int(float(v))
            except (TypeError, ValueError):
                v = d
            v = max(o["least"], min(o["most"], v))
        elif o["type"] == "check":
            v = str(v).lower() in ("1", "true", "yes", "on")
        elif o["type"] == "select":
            v = v if v in [c[0] for c in o["choices"]] else d
        out[o["key"]] = v
    return out


class Report:
    """One run's report. `app` is the server's App: the report reads the run
    through the same calls the results page makes."""

    def __init__(self, app, rundir, raw_options=None):
        self.app, self.rundir = app, rundir
        self.st = app.run_status(rundir)
        if self.st.get("error"):
            raise ValueError(self.st["error"])
        self.j = self.st["job"]
        res = self.st.get("results") or {}
        self.rows = [r for r in (res.get("rows") or []) if isinstance(r, dict)]
        self.top_n = int(res.get("top_n") or self.j.get("top_n") or 5)
        self.key = res.get("ranked_by") or ("mmgbsa_top%d" % self.top_n
                                            if any(r.get("mmgbsa_top%d" % self.top_n) is not None
                                                   for r in self.rows) else "adcp_top%d" % self.top_n)
        self.by_mm = self.key.startswith("mmgbsa")
        # The poses a peptide's MM-GBSA mean is over: the top_n, or fewer when
        # the run rescored fewer.
        try:
            self.n_mm = max(1, min(self.top_n, int(self.j.get("mmgbsa_poses") or self.top_n)))
        except (TypeError, ValueError):
            self.n_mm = self.top_n
        self.constants = res.get("constants") or {}
        self.total = len(self.j.get("peptides") or [])
        self.opt = options(raw_options or {}, max(1, len(self.rows) or self.total))
        self.rec = receptor_name(self.j)
        self.engine = engine_label(self.j, self.rows)
        self.mm_engine = next((r.get("mmgbsa_engine") for r in self.rows if r.get("mmgbsa_engine")),
                              "amber" if any(r.get("mmgbsa_n") for r in self.rows) else "")
        self.ranked = [r for r in self.rows if r.get("rank") is not None]
        self.results = {}
        for r in self.rows:
            self.results[r["peptide"]] = _json(os.path.join(
                rundir, "results", "%s.json" % jobmod.key(r["peptide"]))) or {}
        self._ix = {}

    # -- small readings -----------------------------------------------------------
    def name(self, r):
        return r.get("name") or r.get("peptide")

    def mean(self, r):
        return r.get(self.key)

    def sd(self, r):
        return r.get("mmgbsa_sd") if self.by_mm else r.get("adcp_sd")

    def n_scored(self, r):
        return min(self.n_mm, r.get("mmgbsa_n") or self.n_mm) if self.by_mm else min(self.top_n, r.get("poses") or 0)

    def n_mean(self):
        """How many poses the ranking's mean is over, as a number."""
        return self.n_mm if self.by_mm else self.top_n

    def measure(self):
        return "MM-GBSA" if self.by_mm else "ADCP score"

    def within(self, lead):
        """The ranked peptides whose range over their poses reaches the
        first's: mean minus SD no higher than the first's mean plus SD."""
        m0, s0 = self.mean(lead), self.sd(lead) or 0.0
        out = []
        for r in self.ranked[1:]:
            m, s = self.mean(r), self.sd(r) or 0.0
            if m is not None and m0 is not None and m - s <= m0 + s0:
                out.append(r)
        return out

    def poses_of(self, r):
        """[{rank, score, cluster, mm, skipped, panel}] for one peptide."""
        res = self.results.get(r["peptide"]) or {}
        mm = dict(zip(res.get("mmgbsa_ranks") or [], res.get("mmgbsa") or []))
        skipped = {s.get("pose"): s.get("why") for s in res.get("mmgbsa_poses_skipped") or []}
        panel = {p.get("pose"): p for p in res.get("panel_poses") or [] if isinstance(p, dict)}
        out = []
        for m in res.get("modes") or []:
            k = m.get("rank")
            out.append({"rank": k, "score": m.get("score"), "cluster": m.get("cluster"),
                        "mm": mm.get(k), "skipped": skipped.get(k), "panel": panel.get(k)})
        return out

    def contacts(self, seq, rank):
        """The classified contacts of one pose, on the complex the viewer
        draws, kept while the files are unchanged."""
        pose = os.path.join(self.rundir, "docking", jobmod.key(seq), "out_ranked_%d.pdb" % int(rank))
        rec = (self.j.get("receptor_info") or {}).get("prepared_pdb") or ""
        try:
            stamp = (os.path.getmtime(pose), os.path.getmtime(rec))
        except OSError:
            return None
        k = (self.rundir, seq, int(rank))
        hit = _CONTACTS.get(k)
        if hit and hit[0] == stamp:
            return hit[1]
        text, err = self.app.structure_text(self.rundir, seq, rank, "complex")
        if err:
            return None
        got = classify(text)
        if len(_CONTACTS) > 500:
            _CONTACTS.clear()
        _CONTACTS[k] = (stamp, got)
        return got

    def structure(self, seq, rank, contacts):
        """What reportkit draws a pose with: the complex, seen side on, the
        receptor grey, the peptide orange, the residues in contact as sticks."""
        pep = (contacts or {}).get("chain") or "Z"
        rec_ch = sorted({p["chain"] for p in (contacts or {}).get("pairs") or []}) or \
            list(self.j.get("chains") or ["A"])
        layers = [{"sel": {}, "style": {"cartoon": {"color": "#b9c3cf"}}},
                  {"sel": {"chain": pep}, "style": {"cartoon": {"color": "#e8955a"},
                                                    "stick": {"radius": 0.17, "colorscheme": "orangeCarbon"}}}]
        by = {}
        for r in (contacts or {}).get("receptor") or []:
            by.setdefault(r["chain"], []).append(r["resi"])
        for ch, ns in sorted(by.items()):
            layers.append({"sel": {"chain": ch, "resi": sorted(ns)}, "add": True,
                           "style": {"stick": {"radius": 0.17, "colorscheme": "grayCarbon"}}})
        layers.append({"sel": {"elem": "H"}, "style": {}})
        return {"url": "/api/structure?" + _q(dir=self.rundir, peptide=seq, rank=rank, what="complex"),
                "format": "pdb", "side_on": [rec_ch, [pep]], "layers": layers}

    # -- the document -------------------------------------------------------------
    def build(self):
        j = self.j
        from . import APP, __version__
        created = str(j.get("created") or "")
        who = _who(self.rundir)
        doc = {
            "kicker": "ADCP docking report",
            "title": j.get("label") or os.path.basename(self.rundir.rstrip(os.sep)),
            "subtitle": self.subtitle(),
            "byline": "%s %s, report written %s; run created %s%s" % (
                APP, __version__, time.strftime("%Y-%m-%d"), created[:16] or "on a date it did not record",
                " by " + who if who else ""),
            "path": self.rundir,
            "options": self.option_list(),
            "notes": self.notes(),
            "sections": [s for s in (self.summary(), self.receptor(), self.ranking(), self.peptides()) if s],
            "views": self.views(),
            "views_hint": "None yet. In the results, show a pose and press Add to report on the viewer's bar.",
            "appendix": self.appendix(),
        }
        doc["methods"], doc["references"] = self.methods()
        return doc

    def option_list(self):
        out = []
        for o in OPTIONS:
            o = dict(o, value=self.opt[o["key"]])
            if o["key"] == "peptides":
                o["most"] = max(1, len(self.rows) or self.total)
                o["default"] = min(o["most"], 10)
            out.append(o)
        return out

    def subtitle(self):
        n = self.total
        s = "%s against %s, docked with %s" % (e(self.rec), count(n, "peptide"), e(self.engine))
        if self.mm_engine:
            s += " and rescored with MM-GBSA (%s)" % GB_LABELS.get(self.mm_engine, self.mm_engine)
        return s

    def notes(self):
        from . import report as rmod
        out = []
        note = rmod.engine_note(self.rows, self.j)
        if note:
            out.append({"kind": "p", "tone": "warn", "html": e(note)})
        for w in self.st.get("prep_warnings") or []:
            out.append({"kind": "p", "tone": "warn", "html": e(w)})
        return out

    # -- 1 summary ----------------------------------------------------------------
    def summary(self):
        j, rows = self.j, self.ranked
        bits = []
        mm_n = j.get("mmgbsa_poses")
        first = "%s %s docked against %s with %s" % (
            count(self.total, "peptide").capitalize(), "was" if self.total == 1 else "were",
            e(self.rec), e(self.engine))
        if self.mm_engine:
            first += (", and %s top %s rescored by MM-GBSA (%s)" % (
                "its" if self.total == 1 else "their", count(int(mm_n or self.top_n), "pose"),
                GB_LABELS.get(self.mm_engine, self.mm_engine)))
        bits.append(first + ".")
        if not rows:
            bits.append("None has a result yet: %s finished docking." % (
                "no peptide has" if not self.st.get("done") else count(self.st.get("done"), "peptide") + " have"))
        else:
            lead = rows[0]
            n = self.n_scored(lead) or self.top_n
            unit = " kcal/mol"
            if len(rows) == 1:
                bits.append("<b>%s</b> scores %s%s by %s, the mean over its top %s." % (
                    e(self.name(lead)), pm(self.mean(lead), self.sd(lead)), unit, self.measure(),
                    count(n, "pose")))
            else:
                second = rows[1]
                bits.append("<b>%s</b> (%s) has the most favourable %s, %s%s over its top %s, ahead of %s (%s)." % (
                    e(self.name(lead)), e(lead["peptide"]),
                    "MM-GBSA energy" if self.by_mm else "ADCP score", pm(self.mean(lead), self.sd(lead)), unit,
                    count(n, "pose"), e(self.name(second)), pm(self.mean(second), self.sd(second))))
                near = self.within(lead)
                if near:
                    bits.append("%s %s within one standard deviation of it, so the run does not separate %s." % (
                        and_list(e(self.name(r)) for r in near), "lies" if len(near) == 1 else "lie",
                        "the two" if len(near) == 1 else "them"))
                else:
                    bits.append("No other peptide comes within one standard deviation of it.")
                if self.by_mm:
                    k = "adcp_top%d" % self.top_n
                    by_adcp = sorted([r for r in rows if r.get(k) is not None], key=lambda r: r[k])
                    if by_adcp and by_adcp[0] is lead:
                        bits.append("ADCP's own score ranks it first too (%s)." % num(lead.get(k)))
                    elif by_adcp:
                        pos = by_adcp.index(lead) + 1 if lead in by_adcp else None
                        bits.append("By ADCP score alone %s ranks first (%s)%s." % (
                            e(self.name(by_adcp[0])), num(by_adcp[0].get(k)),
                            ", and %s is %s" % (e(self.name(lead)), _ordinal(pos)) if pos else ""))
            if lead.get("bsa") is not None:
                hb, sb = lead.get("n_hbond"), lead.get("n_ionic")
                site = ""
                if lead.get("site_n"):
                    site = (", touches all %s site residues" % words(lead["site_n"])
                            if lead.get("site_touched") == lead.get("site_n") else
                            ", touches %s of the %s site residues" % (lead.get("site_touched"), lead["site_n"]))
                bits.append("Its pose %s, the one the interface panel describes, buries %s \u00c5\u00b2%s and makes %s and %s." % (
                    lead.get("panel_pose"), num(lead["bsa"], 0), site,
                    count(int(hb or 0), "hydrogen bond"), count(int(sb or 0), "salt bridge")))
        unranked = [r for r in self.rows if r.get("rank") is None]
        if unranked:
            bits.append("%s %s no %s%s." % (
                and_list(e(self.name(r)) for r in unranked), "has" if len(unranked) == 1 else "have",
                self.measure(), " (%s)" % "; ".join(sorted({e(r["mmgbsa_skipped"]) for r in unranked
                                                            if r.get("mmgbsa_skipped")}))
                if any(r.get("mmgbsa_skipped") for r in unranked) else ""))
        missing = self.st.get("outstanding") or []
        if missing:
            bits.append("%s %s not finished docking." % (
                and_list(e(x.get("name") or x.get("sequence")) for x in missing[:6])
                + (" and %d more" % (len(missing) - 6) if len(missing) > 6 else ""),
                "has" if len(missing) == 1 else "have"))
        blocks = [{"kind": "p", "html": " ".join(bits)}]
        blocks.append({"kind": "p", "tone": "caveat", "html":
                       "ADCP scores and MM-GBSA energies order these peptides against this receptor; neither is "
                       "a binding free energy, and MM-GBSA has no entropy term."})
        return {"id": "summary", "heading": "Summary", "label": "summary", "blocks": blocks}

    # -- 2 the receptor and the site ------------------------------------------------
    def receptor(self):
        j = self.j
        blocks = []
        a = j.get("adapt") or {}
        steps, tail, _log = prep_steps(self.rundir)
        if str(j.get("target_file") or "").strip():
            blocks.append({"kind": "p", "html": "The run started from the target file <code>%s</code>: the "
                           "receptor was prepared and its maps built outside this tool, and docked as the file "
                           "gave them." % e(j["target_file"])})
        elif a.get("folder"):
            his = ", ".join("%s %s" % (k, v) for k, v in sorted((a.get("his") or {}).items())) if a.get("his") else ""
            blocks.append({"kind": "p", "html": (
                "The receptor was prepared in ADapt%s (session <code>%s</code>, save %s%s): %s%s, its hydrogens "
                "placed by reduce%s%s. This run docked it as ADapt saved it." % (
                    " " + e(a.get("adapt_version")) if a.get("adapt_version") else "", e(a.get("folder")),
                    e(a.get("save")), ", " + e(a.get("tab")) if a.get("tab") else "",
                    e(a.get("pdbid") or a.get("input") or "the structure"),
                    ", biomolecule " + e(a.get("biomolecule")) if a.get("biomolecule") else "",
                    " " + e(a.get("reduce")) if a.get("reduce") else "",
                    "; hetero groups kept: " + e(", ".join(a.get("hetero") or [])) if a.get("hetero") else ""))
                + (" Histidines: %s." % e(his) if his else "")})
        elif j.get("receptor_prepared"):
            prot = (j.get("receptor_info") or {}).get("protonation") or {}
            blocks.append({"kind": "p", "html": "%s came prepared (%s, pH %s) and was docked as it was given: it "
                           "was not prepared again here." % (e(self.rec), e(prot.get("by") or "before this run"),
                                                             e(j.get("ph", 7.4)))})
        if steps:
            facts = prep_facts(steps, tail, j)
            blocks.append({"kind": "table", "id": "prep", "caption": "%s went into the docking %s." % (
                e(self.rec), ("as this tool prepared it: " + e(and_list(facts))) if facts else
                "as this tool prepared it, step by step") + " Each step is listed in the order it ran.",
                "columns": [{"key": "step", "label": "Step"}, {"key": "out", "label": "Outcome"},
                            {"key": "detail", "label": "What it did"}],
                "rows": [{"step": s.capitalize(), "out": o, "detail": d} for s, o, d in steps],
                "note": "The whole log is in the appendix and in <code>receptor/preparation.log</code>."})
        site = j.get("site") or {}
        mode, hs, how = self.site_how()
        centre, size = site.get("centre"), site.get("size")
        if centre and size:
            dims = " \u00d7 ".join(num(v) for v in size)
            pockets = _pockets(self.rundir)
            cap = "The docking box, %s \u00c5, is %s%s." % (
                dims, how, "; AutoSite found %s in it" % count(pockets, "pocket") if pockets else "")
            maps = j.get("maps") or {}
            rows = [{"k": "Site from", "v": {"residues": "residues", "ligand": "a ligand", "pocket": "an AutoSite pocket",
                                             "target": "a target file"}.get(mode, mode or "residues")},
                    {"k": "Centre (\u00c5)", "v": ", ".join(num(v, 2) for v in centre)},
                    {"k": "Size (\u00c5)", "v": dims},
                    {"k": "Grid spacing (\u00c5)", "v": num(j.get("spacing", 0.375), 3)}]
            if hs:
                rows.append({"k": "Site residues", "v": ", ".join(hs)})
            if maps or j.get("agfr"):
                rows.append({"k": "Maps", "v": "AGFR %s%s, AutoSite %s%s" % (
                    e(maps.get("agfr") or j.get("agfr") or ""),
                    " (%s)" % e(maps["agfr_version"]) if maps.get("agfr_version") else "",
                    e(maps.get("autosite") or j.get("autosite") or ""),
                    ", receptor gradient on" if (maps.get("receptor_gradient") or j.get("receptor_gradient"))
                    else "")})
            if site.get("outside_box"):
                rows.append({"k": "Outside the box", "v": ", ".join("%s%s" % tuple(x[:2]) if isinstance(x, (list, tuple))
                                                                    else str(x) for x in site["outside_box"])})
            blocks.append({"kind": "table", "id": "site", "caption": cap,
                           "columns": [{"key": "k", "label": "What"}, {"key": "v", "label": "Value"}], "rows": rows})
            rec = (j.get("receptor_info") or {}).get("prepared_pdb")
            if rec and os.path.isfile(rec):
                layers = [{"sel": {}, "style": {"cartoon": {"color": "#b9c3cf"}}}]
                by = {}
                for c, n in site.get("hotspots") or []:
                    by.setdefault(c, []).append(int(n))
                for c, ns in sorted(by.items()):
                    layers.append({"sel": {"chain": c, "resi": ns}, "add": True,
                                   "style": {"stick": {"radius": 0.2, "colorscheme": "orangeCarbon"}}})
                layers.append({"sel": {"elem": "H"}, "style": {}})
                blocks.append({"kind": "figure", "id": "sitefig", "caption": (
                    "The box (blue) on %s%s; receptor as a grey cartoon." % (
                        e(self.rec), ", with %s as orange sticks" % _few(hs, "site residue") if hs else "")),
                    "structure": {"url": "/api/structure?" + _q(dir=self.rundir, what="receptor"),
                                  "format": "pdb", "layers": layers,
                                  "shapes": [{"box": {"center": dict(zip("xyz", centre)),
                                                      "dimensions": dict(zip("whd", size)),
                                                      "color": "#1a56db", "wireframe": True}}]}})
        return {"id": "receptor", "heading": "The receptor and the site", "label": "receptor", "blocks": blocks}

    def site_how(self):
        """(where the box came from, the site residues, how the box was placed
        in words)."""
        j = self.j
        site = j.get("site") or {}
        mode = "target" if str(j.get("target_file") or "").strip() else str(j.get("site_mode") or (
            "residues" if site.get("hotspots") else ""))
        hs = ["%s%s" % (c, n) for c, n in site.get("hotspots") or []]
        how = {"residues": "built around %s" % _few(hs, "site residue"),
               "ligand": "built around the ligand %s (%s)%s" % (
                   e(site.get("ligand") or j.get("site_ligand")), e(site.get("ligand_note") or ""),
                   ", which touches %s" % _few(hs, "residue") if hs else ""),
               "pocket": "placed by AutoSite on %s pocket%s, lined by %s" % (
                   "the chosen" if site.get("pocket") == "chosen" else "the best" if site.get("pocket") == "best"
                   else "every", "" if site.get("pocket") == "best" else "s", _few(hs, "residue")),
               "target": "read from the target file"}.get(mode, "built around %s" % _few(hs, "residue"))
        return mode, hs, how

    # -- 3 the peptides ranked ------------------------------------------------------
    def ranking(self):
        rows, tn = self.rows, self.top_n
        if not rows:
            return {"id": "ranking", "heading": "The peptides ranked", "label": "ranking",
                    "blocks": [{"kind": "p", "tone": "caveat", "html": "No peptide has a result yet."}]}
        lead = self.ranked[0] if self.ranked else None
        rings = any((r.get("ring") or "linear") != "linear" for r in rows)
        mm = any(r.get("mmgbsa_n") for r in rows)
        engs = [x for x in ("amber", "openmm", "schrodinger") if any(r.get("mmgbsa_%s_n" % x) for r in rows)]
        cols = [{"key": "rank", "label": "#", "num": True, "empty": "\u2013"},
                {"key": "name", "label": "Peptide"}, {"key": "peptide", "label": "Sequence"},
                {"key": "length", "label": "Length", "num": True}]
        if rings:
            cols.append({"key": "ring", "label": "Ring", "help": "How the peptide was docked: linear, closed head "
                         "to tail, or closed through a disulfide. A ring is rescored only on the poses that close it."})
        cols += [{"key": "adcp_best", "label": "ADCP best", "num": True, "digits": 1},
                 {"key": "adcp", "label": "ADCP top %d" % tn, "num": True,
                  "help": "Mean \u00b1 SD of the ADCP scores of the top %d poses, kcal/mol" % tn}]
        if mm:
            cols += [{"key": "mmgbsa_n", "label": "Rescored", "num": True,
                      "help": "How many poses MM-GBSA rescored"},
                     {"key": "mmgbsa_best", "label": "MM-GBSA best", "num": True, "digits": 1},
                     {"key": "mm", "label": "MM-GBSA top %d" % self.n_mm, "num": True,
                      "help": "Mean \u00b1 SD of the MM-GBSA energies of the top %d poses rescored, kcal/mol"
                      % self.n_mm}]
        if len(engs) > 1:
            cols += [{"key": "eng_" + x, "label": "%s top %d" % (GB_LABELS[x], tn), "num": True, "digits": 1,
                      "help": "%s's mean over the top %d poses; each engine is on its own scale" % (GB_LABELS[x], tn)}
                     for x in engs]
        if any(r.get("cluster_top") is not None for r in rows):
            cols.append({"key": "cluster", "label": "Top cluster", "num": True,
                         "help": "How many of ADCP's clustered solutions fall in the top pose's cluster, and "
                         "their share; an upper bound, as ADCP lists at most 100 clusters"})
        trows = []
        near = {id(r) for r in self.within(lead)} if lead else set()
        for r in rows:
            t = dict(r)
            t["adcp"] = {"html": pm(r.get("adcp_top%d" % tn), r.get("adcp_sd"), 2)}
            t["mm"] = {"html": pm(r.get("mmgbsa_top%d" % tn), r.get("mmgbsa_sd"))}
            for x in engs:
                t["eng_" + x] = r.get("mmgbsa_%s_top%d" % (x, tn))
            if r.get("cluster_top") is not None:
                t["cluster"] = "%s (%s %%)" % (r["cluster_top"], num(100 * (r.get("cluster_frac") or 0), 0))
            t["_class"] = "rk-mark" if r is lead else ("" if id(r) in near else "")
            trows.append(t)
        blocks = [{"kind": "table", "id": "ranking", "caption": self.ranking_caption(), "columns": cols,
                   "rows": trows,
                   "note": self.unranked_note()}]
        fig = self.ranking_chart()
        if fig:
            blocks.append(fig)
        panel = self.panel_table()
        if panel:
            blocks.append(panel)
        return {"id": "ranking", "heading": "The peptides ranked", "label": "ranking", "blocks": blocks}

    def ranking_caption(self):
        rows = self.ranked
        tn = self.top_n
        unit = "kcal/mol"
        how = ("rows are in order of the mean MM-GBSA energy of the top %s rescored, most negative first; "
               "\u00b1 is the SD over those poses, and every energy and score is in %s, lower being better."
               % (count(self.n_mm, "pose"), unit) if self.by_mm else
               "rows are in order of the mean ADCP score of the top %s, most negative first; \u00b1 is the SD "
               "over those poses, in %s, lower being better." % (count(tn, "pose"), unit))
        if not rows:
            return "No peptide could be ranked: none has a %s. %s" % (self.measure(), how[:1].upper() + how[1:])
        lead = rows[0]
        if len(rows) == 1:
            return "%s is the only peptide ranked, at %s %s by %s; %s" % (
                e(self.name(lead)), pm(self.mean(lead), self.sd(lead)), unit, self.measure(), how)
        second = rows[1]
        gap = (self.mean(second) or 0) - (self.mean(lead) or 0)
        return "%s ranks first by %s, %s %s, %s ahead of %s (%s); %s" % (
            e(self.name(lead)), self.measure(), pm(self.mean(lead), self.sd(lead)), unit, num(gap),
            e(self.name(second)), pm(self.mean(second), self.sd(second)), how)

    def unranked_note(self):
        un = [r for r in self.rows if r.get("rank") is None]
        if not un:
            return ""
        return "Not ranked, having no %s: %s." % (self.measure(), "; ".join(
            e(self.name(r)) + (" (%s)" % e(r["mmgbsa_skipped"]) if r.get("mmgbsa_skipped") else "") for r in un))

    def ranking_chart(self):
        rows = self.ranked
        if len(rows) < 2:
            return None
        lead = rows[0]
        near = {id(r) for r in self.within(lead)}
        tn = self.top_n
        charts = []
        if self.by_mm:
            charts.append({"axis": "MM-GBSA, mean of the top %d poses (kcal/mol)" % self.n_mm, "digits": 1, "points": [
                {"name": self.name(r), "v": r.get("mmgbsa_top%d" % tn), "sd": r.get("mmgbsa_sd"),
                 "dim": not (r is lead or id(r) in near)} for r in rows if r.get("mmgbsa_top%d" % tn) is not None]})
        adcp = [r for r in rows if r.get("adcp_top%d" % tn) is not None]
        adcp_near = near if not self.by_mm else None
        charts.append({"axis": "ADCP score, mean of the top %d poses (kcal/mol)" % tn, "digits": 2, "points": [
            {"name": self.name(r), "v": r.get("adcp_top%d" % tn), "sd": r.get("adcp_sd"),
             "dim": (not (r is lead or id(r) in adcp_near)) if adcp_near is not None else False} for r in adcp]})
        inside = [lead] + [r for r in rows if id(r) in near]
        outside = [r for r in rows if r is not lead and id(r) not in near]
        if outside and len(inside) > 1:
            head = ("%s cannot be told apart by %s, their ranges over the top %s overlapping %s's (dark bars); "
                    "%s %s less favourably (pale)" % (
                        and_list(e(self.name(r)) for r in inside), self.measure(), count(self.n_mean(), "pose"),
                        e(self.name(lead)),
                        and_list(e(self.name(r)) for r in outside), "scores" if len(outside) == 1 else "score"))
        elif outside:
            head = "%s stands clear of every other peptide by %s, by more than either's SD (dark bar)" % (
                e(self.name(lead)), self.measure())
        else:
            head = "No peptide stands clear of %s by %s: every range over the top %s overlaps its own" % (
                e(self.name(lead)), self.measure(), count(self.n_mean(), "pose"))
        return {"kind": "figure", "id": "rankchart", "bars": charts, "caption": head + (
            ". %s; bars are means over the top poses (%s), error bars one SD, lower is better." % (
                "Top, MM-GBSA; bottom, ADCP score" if self.by_mm else "ADCP score",
                "%s for MM-GBSA, %s for ADCP" % (words(self.n_mm), words(tn)) if self.by_mm else words(tn)))}

    def panel_table(self):
        rows = [r for r in self.rows if r.get("bsa") is not None]
        if not rows:
            return None
        cols = [{"key": "rank", "label": "#", "num": True, "empty": "\u2013"}, {"key": "name", "label": "Peptide"},
                {"key": "panel_pose", "label": "Pose", "num": True,
                 "help": "Which ADCP pose the row describes: the one with the lowest MM-GBSA energy"},
                {"key": "bsa", "label": "BSA (\u00c5\u00b2)", "num": True, "digits": 0,
                 "help": "Buried surface area, both partners together"},
                {"key": "f_pep_buried", "label": "Buried", "num": True, "digits": 2,
                 "help": "The share of the peptide's own surface the receptor covers"},
                {"key": "f_hydrophobic", "label": "Apolar", "num": True, "digits": 2,
                 "help": "The share of the peptide's buried area that is carbon or sulphur"},
                {"key": "n_hbond", "label": "HB", "num": True, "help": "N or O atom pairs within 3.5 \u00c5"},
                {"key": "n_ionic", "label": "SB", "num": True, "help": "Salt-bridged residue pairs"},
                {"key": "n_stack", "label": "Stack", "num": True, "help": "Aromatic residue pairs within 5.0 \u00c5"},
                {"key": "n_clash", "label": "Clash", "num": True,
                 "help": "Heavy atom pairs closer than 0.80 of the sum of their radii"}]
        if any(r.get("site_n") for r in rows):
            cols += [{"key": "site", "label": "Site", "num": True,
                      "help": "Site residues the pose touches, of how many"},
                     {"key": "site_dg", "label": "Site \u0394G", "num": True, "digits": 1,
                      "help": "The site residues' share of the MM-GBSA energy, kcal/mol"}]
        if any(r.get("pe_verdict") for r in rows):
            cols.append({"key": "pe", "label": "HOPE-PE", "help": "Sequence liabilities: verdict and penalty"})
        out = []
        for r in rows:
            t = dict(r)
            if r.get("site_n"):
                t["site"] = "%s/%s" % (r.get("site_touched"), r.get("site_n"))
            if r.get("pe_verdict"):
                t["pe"] = "%s (%s)" % (r["pe_verdict"], r.get("pe_score"))
            out.append(t)
        big = max(rows, key=lambda r: r.get("bsa") or 0)
        hb = max(rows, key=lambda r: r.get("n_hbond") or 0)
        if len(rows) == 1:
            one = rows[0]
            site = ("; it touches %s of the %s site residues" % (
                "all" if one.get("site_touched") == one.get("site_n") else one.get("site_touched"), one["site_n"])
                if one.get("site_n") else "")
            return {"kind": "table", "id": "panel", "columns": cols, "rows": out, "caption": (
                "%s's pose %s buries %s \u00c5\u00b2 and makes %s and %s%s. The pose with the lowest MM-GBSA energy, "
                "measured on its minimised complex; HB counts atom pairs, SB and Stack residue pairs." % (
                    e(self.name(one)), one.get("panel_pose"), num(one["bsa"], 0),
                    count(int(one.get("n_hbond") or 0), "hydrogen bond"),
                    count(int(one.get("n_ionic") or 0), "salt bridge"), site))}
        if big is hb:
            lead = "%s's pose buries the most surface, %s \u00c5\u00b2, and makes the most hydrogen bonds, %s" % (
                e(self.name(big)), num(big["bsa"], 0), hb.get("n_hbond"))
        else:
            lead = "%s's pose buries the most surface, %s \u00c5\u00b2, and %s's makes the most hydrogen bonds, %s" % (
                e(self.name(big)), num(big["bsa"], 0), e(self.name(hb)), hb.get("n_hbond"))
        site = ""
        sited = [r for r in rows if r.get("site_n")]
        if sited:
            full = [r for r in sited if r.get("site_touched") == r.get("site_n")]
            if len(full) == len(sited):
                site = "; every peptide touches all %s site residues" % words(sited[0]["site_n"]) if len(sited) > 1 \
                    else "; it touches all %s site residues" % words(sited[0]["site_n"])
            elif len({(r.get("site_touched"), r.get("site_n")) for r in sited}) == 1:
                site = "; every peptide touches %s of the %s site residues" % (
                    sited[0].get("site_touched"), sited[0].get("site_n"))
            else:
                lo = min(sited, key=lambda r: r.get("site_touched") or 0)
                hi = max(sited, key=lambda r: r.get("site_touched") or 0)
                site = "; site coverage runs from %s of %s (%s) to %s of %s (%s)" % (
                    lo.get("site_touched"), lo.get("site_n"), e(self.name(lo)), hi.get("site_touched"),
                    hi.get("site_n"), e(self.name(hi)))
        clash = [r for r in rows if r.get("n_clash")]
        clash_s = ("; %s clash" % and_list(e(self.name(r)) for r in clash)) if clash else ""
        return {"kind": "table", "id": "panel", "columns": cols, "rows": out, "caption": (
            "%s%s%s. Each row is one pose, the one with the lowest MM-GBSA energy (Pose), measured on the "
            "minimised complex; HB counts atom pairs, SB and Stack residue pairs." % (lead, site, clash_s))}

    # -- 4 each peptide -------------------------------------------------------------
    def peptides(self):
        rows = self.rows[:self.opt["peptides"]]
        subs = [self.peptide(r) for r in rows]
        if not subs:
            return None
        left = len(self.rows) - len(rows)
        blocks = []
        if left > 0:
            blocks.append({"kind": "p", "tone": "caveat", "html": "The top %s of %d are shown; raise Peptides on the "
                           "bar to show the rest." % (count(len(rows), "peptide"), len(self.rows))})
        return {"id": "peptides", "heading": "Each peptide", "label": "each peptide", "blocks": blocks,
                "sections": subs}

    def peptide(self, r):
        seq, name = r["peptide"], self.name(r)
        res = self.results.get(seq) or {}
        poses = self.poses_of(r)
        blocks = []
        # what it is, where it ranks, how it was docked
        intro = []
        if r.get("rank") is not None:
            intro.append("Rank %s of %d by %s, %s kcal/mol over its top %s." % (
                r["rank"], len(self.ranked), self.measure(), pm(self.mean(r), self.sd(r)),
                count(self.n_scored(r) or self.top_n, "pose")))
        else:
            intro.append("Not ranked: no %s%s." % (self.measure(), " (%s)" % e(r["mmgbsa_skipped"])
                                                    if r.get("mmgbsa_skipped") else ""))
        ring = r.get("ring") or "linear"
        from . import peptide as pepmod
        mode = res.get("cyclic") or jobmod.ring_for(self.j, seq)
        intro.append("%s residues, docked %s by %s; ADCP kept %s." % (
            r.get("length"), pepmod.RING_AS.get(mode, "as " + ring), e(r.get("engine") or self.engine),
            count(len(poses), "pose")))
        if r.get("cluster_top") is not None and res.get("cluster_total"):
            intro.append("Its top pose stands for %s of the %s solutions ADCP clustered into its %s modes (%s %%), "
                         "an upper bound since ADCP lists at most 100 clusters." % (
                             r["cluster_top"], res["cluster_total"], r.get("cluster_modes") or res.get("cluster_modes"),
                             num(100.0 * r["cluster_top"] / res["cluster_total"], 0)))
        blocks.append({"kind": "p", "html": " ".join(intro)})
        skipped = [p for p in poses if p["skipped"]]
        if ring != "linear" and skipped:
            blocks.append({"kind": "p", "html": "MM-GBSA rescored only the poses that close the ring; %s did not: %s." % (
                and_list("pose %s" % p["rank"] for p in skipped), "; ".join(
                    "pose %s, %s" % (p["rank"], e(p["skipped"])) for p in skipped))})
        elif skipped:
            blocks.append({"kind": "p", "html": "Poses not rescored: %s." % "; ".join(
                "pose %s, %s" % (p["rank"], e(p["skipped"])) for p in skipped)})
        if res.get("mmgbsa_skipped") and not r.get("mmgbsa_n"):
            blocks.append({"kind": "p", "html": "No MM-GBSA energy: %s." % e(res["mmgbsa_skipped"])})
        if poses:
            blocks.append(self.poses_table(r, poses))
        # the top poses, drawn
        shown = [p for p in poses][:self.opt["poses"]]
        ix = {p["rank"]: self.contacts(seq, p["rank"]) for p in shown}
        if shown:
            panels = []
            for k, p in enumerate(shown):
                c = ix.get(p["rank"])
                panels.append({"structure": self.structure(seq, p["rank"], c), "label": "(%s) pose %s, ADCP %s%s" % (
                    "abcde"[k], p["rank"], num(p["score"]), ", MM-GBSA %s" % num(p["mm"]) if p["mm"] is not None else "")})
            touch = [len((ix.get(p["rank"]) or {}).get("receptor") or []) for p in shown]
            blocks.append({"kind": "figure", "id": "poses_" + jobmod.key(seq), "panels": panels,
                           "columns": min(3, len(panels)), "caption": (
                               "%s's top %s on %s, scoring %s kcal/mol by ADCP, %s in contact with %s receptor "
                               "residues. Receptor grey, peptide orange, residues in contact (heavy atoms within "
                               "4 \u00c5, or ring carbons within 5 \u00c5) as sticks, each seen side on across the "
                               "interface." % (
                                   e(name), count(len(shown), "pose") if len(shown) > 1 else "pose",
                                   e(self.rec), and_list(num(p["score"]) for p in shown),
                                   "are" if len(shown) > 1 else "is", and_list(str(t) for t in touch)))})
        # the energy of the pose the panel describes, in its parts
        pp = r.get("panel_pose")
        if r.get("mm_vdw") is not None and pp is not None:
            dg = next((p["mm"] for p in poses if p["rank"] == pp), None)
            blocks.append({"kind": "p", "html": "Pose %s's MM-GBSA energy%s is van der Waals %s, electrostatic %s, polar "
                           "solvation %s and nonpolar solvation %s kcal/mol." % (
                               pp, ", %s," % num(dg) if dg is not None else "", num(r["mm_vdw"]), num(r["mm_eel"]),
                               num(r["mm_egb"], 1, True), num(r["mm_esurf"]))})
        # interaction diagrams
        which = [pp] if pp is not None else ([poses[0]["rank"]] if poses else [])
        if self.opt["diagrams"] == "all":
            which = list(dict.fromkeys(which + [p["rank"] for p in shown]))
        on = ("ionic", "stack", "hbond", "vdw") if self.opt["vdw"] else ("ionic", "stack", "hbond")
        for k in which:
            c = ix.get(k) if k in ix else self.contacts(seq, k)
            fig = self.lid_figure(name, k, c, on, "lid_%s_%s" % (jobmod.key(seq), k))
            if fig:
                blocks.append(fig)
        blocks += self.residue_tables(r)
        if r.get("omm_n"):
            blocks.append({"kind": "p", "html": "ADCP 1.1's OpenMM minimisation of %s: the lowest E(complex) "
                           "\u2212 E(receptor), %s kcal/mol, is pose %s's." % (
                               count(int(r["omm_n"]), "pose"), num(r.get("omm_best")), r.get("omm_best_pose"))})
        if r.get("pe_verdict"):
            flags = r.get("pe_flags") or []
            blocks.append({"kind": "p", "html": "HOPE-PE %s the sequence (%s, penalty %s)%s" % (
                "passes" if r["pe_verdict"] == "PASS" else "flags" if r["pe_verdict"] == "FLAG" else "fails",
                e(r["pe_verdict"]), e(r.get("pe_score")), ":" if flags else ".")})
            if flags:
                blocks.append({"kind": "list", "items": [e(f) for f in flags]})
        return {"id": "p_" + jobmod.key(seq), "heading": "%s (%s)" % (name, seq), "newpage": True, "blocks": blocks}

    def poses_table(self, r, poses):
        name = self.name(r)
        has_mm = any(p["mm"] is not None for p in poses)
        has_cl = any(p["cluster"] for p in poses)
        has_pan = any(p["panel"] for p in poses)
        cols = [{"key": "rank", "label": "Pose", "num": True},
                {"key": "score", "label": "ADCP", "num": True, "digits": 1, "help": "ADCP's score, kcal/mol"}]
        if has_cl:
            cols.append({"key": "cluster", "label": "Cluster", "num": True,
                         "help": "How many of ADCP's solutions this pose's cluster holds"})
        if has_mm or any(p["skipped"] for p in poses):
            cols.append({"key": "mm", "label": "MM-GBSA", "num": True, "digits": 1,
                         "help": "MM-GBSA \u0394G, kcal/mol", "empty": ""})
        if has_pan:
            cols += [{"key": "bsa", "label": "BSA (\u00c5\u00b2)", "num": True, "digits": 0},
                     {"key": "hb", "label": "HB", "num": True}, {"key": "sb", "label": "SB", "num": True}]
            if any((p["panel"] or {}).get("site_n") for p in poses):
                cols.append({"key": "site", "label": "Site", "num": True})
        rows = []
        for p in poses:
            pan = p["panel"] or {}
            t = {"rank": p["rank"], "score": p["score"], "cluster": p["cluster"],
                 "mm": p["mm"] if p["mm"] is not None else ({"html": "not rescored"} if p["skipped"] else None),
                 "bsa": pan.get("bsa"), "hb": pan.get("n_hbond"), "sb": pan.get("n_ionic"),
                 "site": "%s/%s" % (pan.get("site_touched"), pan.get("site_n")) if pan.get("site_n") else None}
            if p["rank"] == r.get("panel_pose"):
                t["_class"] = "rk-mark"
            rows.append(t)
        scored = [p for p in poses if p["mm"] is not None]
        first = poses[0]
        if scored:
            best = min(scored, key=lambda p: p["mm"])
            if best is first:
                lead = "Pose 1 is best by both ADCP (%s) and MM-GBSA (%s kcal/mol)" % (num(first["score"]), num(best["mm"]))
                rest = [p for p in scored if p is not best]
                if rest:
                    lead += "; the other poses rescored lie between %s and %s" % (
                        num(min(p["mm"] for p in rest)), num(max(p["mm"] for p in rest)))
            else:
                lead = ("Pose %s has the most favourable MM-GBSA energy, %s kcal/mol, though ADCP ranks it %s (%s); "
                        "pose 1, ADCP's best (%s), scores %s" % (
                            best["rank"], num(best["mm"]), _ordinal(best["rank"]), num(best["score"]),
                            num(first["score"]), num(first["mm"]) if first["mm"] is not None else "nothing, not rescored"))
        elif len(poses) > 1:
            lead = "Pose 1 has the best ADCP score, %s kcal/mol, %s ahead of pose 2" % (
                num(first["score"]), num((poses[1]["score"] or 0) - (first["score"] or 0)))
        else:
            lead = "ADCP kept one pose, scoring %s kcal/mol" % num(first["score"])
        return {"kind": "table", "id": "poses_t_" + jobmod.key(r["peptide"]), "columns": cols, "rows": rows,
                "caption": "%s. %s's poses in ADCP's order; scores and energies in kcal/mol, lower is better%s." % (
                    lead, e(name), "; the marked row is the pose the interface panel describes, HB counting atom "
                    "pairs and SB residue pairs" if has_pan else "")}

    def lid_figure(self, name, rank, c, on, fid):
        if c is None:
            return {"kind": "p", "tone": "caveat", "html": "Pose %s's structure could not be read, so it has no "
                    "interaction diagram." % rank}
        pairs = c.get("pairs") or []
        svg = lid_svg(pairs, on)
        shown = [p for p in pairs if p["kind"] in on]
        if not svg:
            return {"kind": "p", "tone": "caveat", "html": "Pose %s makes no %s with the receptor, so it has no "
                    "interaction diagram." % (rank, "contact" if "vdw" in on else
                                              "hydrogen bond, salt bridge or stacking")}
        return {"kind": "figure", "id": fid, "svg": svg, "caption": "%s; lines are coloured by kind (key below), %s, "
                "and each receptor residue sits beside the peptide residue it is closest to, on the docked pose." % (
                    pair_sentence(shown, name, self.rec, rank),
                    "van der Waals contacts drawn too" if "vdw" in on else "van der Waals contacts left out")}

    def rec_res(self, resn, chain, n):
        """A receptor residue as the report names it: 'Glu82', with its chain
        ('Glu A82') when the receptor has more than one."""
        multi = len(self.j.get("chains") or []) > 1
        return res_name(resn, ("%s%s" % (chain, n)) if multi and chain else n)

    def residue_tables(self, r):
        pep, rec, pose = residue_map(self.rundir, r["peptide"])
        if pep is None:
            return []
        name, key = self.name(r), jobmod.key(r["peptide"])
        out = []

        def kinds(s):
            got = (s or "").split(",")
            return ", ".join(KIND_WORD[k] for k in SPECIFIC if k in got) + (", clash" if "clash" in got else "")
        cols = [{"key": "res", "label": "Residue"}, {"key": "kinds", "label": "Contacts"},
                {"key": "d", "label": "Closest (\u00c5)", "num": True, "digits": 1,
                 "help": "The closest heavy-atom distance to the other partner"},
                {"key": "dg", "label": "Energy (kcal/mol)", "num": True, "digits": 2,
                 "help": "The residue's share of the MM-GBSA energy; negative helps binding"}]
        if pep:
            # The map is measured on the minimised complex, as the energy is.
            have = [x for x in pep if x["dg"] is not None]
            best = min(have, key=lambda x: x["dg"]) if have else None
            worst = max(have, key=lambda x: x["dg"]) if have else None
            rows = [{"res": res_name(x.get("resname"), x["pos"]), "kinds": kinds(x.get("kinds")),
                     "d": x["d_min"], "dg": x["dg"], "_class": "rk-mark" if x is best else ""} for x in pep]
            if best:
                cap = "Of %s's residues, %s contributes most to pose %s's MM-GBSA energy, %s kcal/mol" % (
                    e(name), res_name(best.get("resname"), best["pos"]), pose, num(best["dg"], 2))
                if worst is not None and worst["dg"] > 0.5:
                    cap += ", and %s costs the most, %s" % (res_name(worst.get("resname"), worst["pos"]),
                                                            num(worst["dg"], 2, True))
            else:
                cap = "%s's residues in pose %s, none with a share of the energy" % (e(name), pose)
            out.append({"kind": "table", "id": "pres_" + key, "columns": cols, "rows": rows,
                        "caption": cap + ". Every peptide residue, with its contacts, its closest heavy-atom distance "
                        "to the receptor and its share of the energy (negative helps); the shares add up to within "
                        "a few per cent of the total."})
        if rec:
            rec = sorted(rec, key=lambda x: (x["dg"] is None, x["dg"] if x["dg"] is not None else 0))
            best = rec[0] if rec[0]["dg"] is not None else None
            site = [x for x in rec if str(x.get("site")) == "1"]
            rows = [{"res": {"html": e(self.rec_res(x.get("resname"), x.get("chain"), x.get("resnum")))
                             + (" <b>site</b>" if str(x.get("site")) == "1" else "")},
                     "kinds": kinds(x.get("kinds")), "d": x["d_min"], "dg": x["dg"],
                     "_class": "rk-mark" if x is best else ""} for x in rec]
            if best:
                cap = "On the receptor, %s contributes most, %s kcal/mol, of the %d residues pose %s touches once minimised" % (
                    e(self.rec_res(best.get("resname"), best.get("chain"), best.get("resnum"))),
                    num(best["dg"], 2), len(rec), pose)
            else:
                cap = "Pose %s touches %d receptor residues" % (pose, len(rec))
            if site:
                cap += ", %s of them site residues (marked)" % words(len(site))
            out.append({"kind": "table", "id": "rres_" + key, "columns": cols, "rows": rows,
                        "caption": cap + ". Most favourable first; columns as in {{table:pres_%s}}." % key})
        return out

    # -- the views added from the viewer --------------------------------------------
    def views(self):
        out = []
        for it in load_items(self.rundir):
            seq, rank = it.get("peptide") or "", it.get("rank")
            r = next((x for x in self.rows if x["peptide"] == seq), None)
            name = self.name(r) if r else seq
            pose = next((p for p in self.poses_of(r) if p["rank"] == rank), None) if r else None
            c = self.contacts(seq, rank) if seq and rank else None
            bits = []
            if pose:
                bits.append("ADCP %s" % num(pose["score"]))
                if pose["mm"] is not None:
                    bits.append("MM-GBSA %s kcal/mol" % num(pose["mm"]))
            cap = "Pose %s of %s as drawn in the viewer (%s)%s%s." % (
                e(rank), e(name), {"molstar": "Mol*", "3dmol": "3Dmol"}.get(it.get("viewer"), "3Dmol"),
                ": " + ", ".join(bits) if bits else "",
                ", touching %s" % count(len(c.get("receptor") or []), "receptor residue") if c else "")
            v = {"id": it.get("id"), "caption": cap, "own": it.get("own") or "",
                 "png": "/api/report/image?" + _q(rundir=self.rundir, id=it.get("id")), "blocks": []}
            if c:
                on = ("ionic", "stack", "hbond", "vdw") if self.opt["vdw"] else ("ionic", "stack", "hbond")
                fig = self.lid_figure(name, rank, c, on, "vlid_%s" % it.get("id"))
                if fig:
                    v["blocks"].append(fig)
            out.append(v)
        return out

    # -- methods and references -----------------------------------------------------
    def methods(self):
        """(the methods as HTML with [n] citations, the references they cite).

        Written from what the run did: its preparation log, its job file, the
        engine its rows name. References are numbered in the order first cited;
        two keys with one paper (PDBFixer and OpenMM) share its number."""
        from . import APP, __version__, cite
        j = self.j
        keys = cite.used_by(j) | {"3dmol"}
        refs = {k: "%s doi:%s" % (ref, doi) for k, _w, _d, ref, doi in cite.REFERENCES}
        refs["tool"] = cite.self_citation()
        texts = []

        def c(*ks):
            nums = []
            for k in ks:
                if k != "tool" and (k not in keys or k not in refs):
                    continue
                if refs[k] not in texts:
                    texts.append(refs[k])
                nums.append(texts.index(refs[k]) + 1)
            return " [%s]" % ", ".join(str(n) for n in sorted(set(nums))) if nums else ""
        paras = []
        a = j.get("adapt") or {}
        prot = "reduce" if a.get("folder") else jobmod.protonation_of(j)
        rinfo = (j.get("receptor_info") or {}).get("protonation") or {}
        steps, tail, _l = prep_steps(self.rundir)
        done = {s_: d for s_, o, d in steps if o == "done"}
        tool = "%s %s%s" % (APP, __version__, c("tool"))
        if str(j.get("target_file") or "").strip():
            t = "Docking used a target file prepared outside this tool with AGFR%s, read in by %s." % (c("adfr"), tool)
        elif a.get("folder"):
            t = ("The receptor, %s, was prepared in ADapt, its hydrogens placed by reduce%s, and docked as ADapt "
                 "saved it; the docking was set up and analysed with %s." % (e(self.rec), c("reduce"), tool))
        elif j.get("receptor_prepared"):
            t = ("The receptor, %s, came prepared and was docked as given; the docking was set up and analysed "
                 "with %s." % (e(self.rec), tool))
        else:
            did = []
            if "waters" in done and not re.search(r"removed 0\b", done["waters"]):
                did.append("its waters removed")
            if "heteroatoms" in done:
                did.append("its heteroatoms removed")
            if "alternate locations" in done:
                did.append("alternate locations reduced to the most occupied")
            if "missing atoms" in done or not steps:
                built = re.search(r"built (\d+) missing residue", done.get("missing atoms", ""))
                did.append("missing atoms%s rebuilt with PDBFixer%s" % (" and residues" if built else "", c("pdbfixer")))
            ver = str(rinfo.get("version") or "").replace("reduce.", "")
            did.append("hydrogens added with reduce%s%s" % (" " + e(ver) if ver else "", c("reduce"))
                       if prot == "reduce" else
                       "hydrogens added by PDB2PQR%s with PROPKA%s titration states" % (c("pdb2pqr"), c("propka")))
            t = "The receptor, %s, was prepared with %s at pH %s: %s." % (e(self.rec), tool, e(j.get("ph", 7.4)),
                                                                        and_list(did))
            if j.get("minimise", True) and ("minimisation" in done or not steps):
                t += (" It was then minimised in OpenMM%s with ff14SB%s and the GBn2 implicit solvent%s, heavy atoms "
                      "restrained at %s kcal/mol/\u00c5\u00b2, for at most %s steps." % (
                          c("openmm"), c("ff14sb"), c("gbn2"), num(j.get("restraint", 10.0), 0),
                          e(j.get("minimise_steps", 500))))
        paras.append(("Receptor", t))
        site = j.get("site") or {}
        if site.get("centre") and site.get("size"):
            mode, hs, how = self.site_how()
            maps = j.get("maps") or {}
            t = ("Affinity maps were computed with AGFR%s%s and AutoGrid4%s over a %s \u00c5 box at %s \u00c5 "
                 "spacing, %s, and AutoSite%s placed the points the search starts from." % (
                     " " + e(maps.get("agfr")) if maps.get("agfr") else "", c("adfr"), c("autogrid"),
                     " \u00d7 ".join(num(v) for v in site["size"]), num(j.get("spacing", 0.375), 3), how,
                     c("autosite")))
            paras.append(("Site", t))
        from . import peptide as pepmod
        by_ring = {}
        for (nm, _sq), p in zip(jobmod.peptide_list(j), j.get("peptides") or []):
            by_ring.setdefault(jobmod.ring_of(j, p), []).append(nm)
        ringed = sorted(x for x in by_ring if x != "none")
        if ringed and len(by_ring) == 1:
            rtext = "; every peptide was docked %s" % pepmod.RING_AS.get(ringed[0], ringed[0])
        else:
            rtext = "".join("; %s %s docked %s" % (and_list(e(n) for n in by_ring[x]),
                                                   "was" if len(by_ring[x]) == 1 else "were",
                                                   pepmod.RING_AS.get(x, x)) for x in ringed)
        t = ("Each peptide was docked with %s%s (AutoDock CrankPep), %s replicas of %s Monte Carlo steps, keeping "
             "the top %s, each the representative of a cluster of solutions%s." % (
                 e(self.engine), c("adcp"), "{:,}".format(int(j.get("replicas") or 0)),
                 "{:,}".format(int(j.get("steps") or 0)), count(int(j.get("poses") or 0), "pose"), rtext))
        if j.get("helical"):
            t += " %s %% of the replicas started as a helix." % e(j.get("helical"))
        if (j.get("converted") or {}).get("from") == "ADCP 1.0":
            t += (" The peptides were docked with ADCP 1.0 (ADFRsuite 1.0) before this tool moved to ADCP 1.1; the "
                  "run was converted on %s, its poses and scores kept as 1.0 wrote them." % e(j["converted"].get("date")))
        if any(r.get("engine") == "ADCP 1.0 (macOS)" for r in self.rows):
            t += " ADCP 1.0 ran on a Mac, as Scripps builds no ADCP 1.1 for macOS."
        if int(j.get("omm_minimize") or 0):
            t += (" ADCP minimised the top %s poses in OpenMM%s with Amber ff99SB%s%s." % (
                e(j.get("omm_minimize")), c("openmm"), c("ff99sb"),
                " and GBn2%s" % c("gbn2") if (j.get("omm_env") or "implicit") == "implicit" else ""))
        if any(pepmod.nonstandard((p.get("sequence") or "").upper()) for p in j.get("peptides") or []):
            t += " Non-standard residues came from ADCP 1.1's SwissSidechain library%s." % c("swisssidechain")
        paras.append(("Docking", t))
        if self.mm_engine:
            model = j.get("mmgbsa_model") or "obc2"
            if self.mm_engine == "amber":
                t = ("The top %s of each peptide were rescored by MM-GBSA in AmberTools%s: tleap built each complex "
                     "with ff14SB%s, sander minimised it with the receptor restrained, and MMPBSA.py%s took the "
                     "energy of that one structure in %s, with no entropy term. A peptide's energy is the mean over "
                     "its top %s rescored." % (
                         count(int(j.get("mmgbsa_poses") or 0), "pose"), c("amber"), c("ff14sb"), c("mmpbsa"),
                         "GBn2 (igb=8)%s" % c("gbn2") if model == "gbn2" else
                         "GB-OBC(II) (igb=5)%s with mbondi2 radii and 0.15 M salt" % c("gbmodel"),
                         count(self.n_mm, "pose")))
            else:
                t = ("The top %s of each peptide were rescored by MM-GBSA with %s on the geometry sander "
                     "minimised%s, with no entropy term; its energies are on its own scale." % (
                         count(int(j.get("mmgbsa_poses") or 0), "pose"),
                         GB_LABELS.get(self.mm_engine, self.mm_engine), c("amber")))
            if any(r.get("bsa") is not None for r in self.rows):
                t += (" Each peptide's lowest-energy pose was measured on its minimised complex: buried surface by "
                      "Shrake and Rupley with Bondi radii and a 1.4 \u00c5 probe, and each residue's share of the "
                      "energy by MMPBSA.py's decomposition of the same structure.")
            paras.append(("Rescoring", t))
        paras.append(("Contacts", "Contacts in the interaction diagrams were classified by distance between heavy "
                      "atoms on the docked pose: a hydrogen bond, nitrogen, oxygen or sulphur atoms within 3.5 \u00c5; "
                      "a salt bridge, a basic and an acidic side-chain tip within 4.0 \u00c5; stacking, carbons of two "
                      "aromatic residues (Phe, Tyr, Trp, His) within 5.0 \u00c5; van der Waals, any other heavy atoms "
                      "within 4.0 \u00c5. Angles were not assessed."))
        paras.append(("Figures", "Structures were drawn with 3Dmol.js%s." % c("3dmol")))
        return ("".join("<p><b>%s.</b> %s</p>" % (h, t) for h, t in paras),
                [{"n": i, "text": t} for i, t in enumerate(texts, 1)])

    # -- the appendix ---------------------------------------------------------------
    def software(self):
        """[[what, html]]: what wrote the run and what it ran, from the run's own files,
        who ran it and when, and where anything in it came from."""
        from . import APP, __version__, page
        j = self.j
        rows = []
        wb = j.get("written_by") or {}
        if wb.get("version"):
            rows.append(["Written by", "%s %s%s" % (e(wb.get("app") or APP), e(wb["version"]),
                                                   ", commit %s" % e(wb["commit"]) if wb.get("commit") else "")])
        else:
            try:
                cite = open(os.path.join(self.rundir, "HOW_TO_CITE.txt"), errors="replace").read()
            except OSError:
                cite = ""
            m = re.search(r"%s v([\w.]+)" % re.escape(APP), cite)
            rows.append(["Written by", ("%s %s, as the run's HOW_TO_CITE.txt names it; the commit was not recorded"
                                        % (e(APP), e(m.group(1)))) if m else "not recorded in the run"])
        rows.append(["This report", "%s %s, commit %s" % (e(APP), e(__version__), e(page._build()))])
        who = _who(self.rundir)
        rows.append(["Run", "created %s%s; jobs %s" % (e(j.get("created") or "at a time not recorded"),
                                                     " by " + e(who) if who else "",
                                                     e(", ".join("%s %s" % kv for kv in (j.get("jobs") or {}).items()))
                                                     or "none")])
        rows.append(["Docking", e(", ".join(sorted({r.get("engine") for r in self.rows if r.get("engine")})) or self.engine)])
        maps = j.get("maps") or {}
        if maps or j.get("agfr"):
            rows.append(["Maps", "AGFR %s%s, AutoSite %s" % (
                e(maps.get("agfr") or j.get("agfr") or ""), " (%s)" % e(maps["agfr_version"]) if maps.get("agfr_version")
                else "", e(maps.get("autosite") or j.get("autosite") or ""))])
        prot = (j.get("receptor_info") or {}).get("protonation") or {}
        if prot.get("version") or prot.get("method"):
            rows.append(["Hydrogens", "%s%s" % (e(prot.get("method") or ""), " " + e(prot["version"]) if prot.get("version")
                                                else "")])
        if self.mm_engine:
            rows.append(["MM-GBSA", "%s, %s" % (GB_LABELS.get(self.mm_engine, self.mm_engine),
                                                 e(j.get("mmgbsa_model") or "obc2"))])
        env = {}
        for name in ("02_dock.sbatch", "01_prepare.sbatch"):
            try:
                for L in open(os.path.join(self.rundir, name), errors="replace"):
                    m = re.match(r"^export (AMBERHOME|ADCPSUITE_HOME|ADFRSUITE_HOME|ADAPT_HOME)=(\S+)", L.strip())
                    if m:
                        env.setdefault(m.group(1), m.group(2))
            except OSError:
                continue
        if env:
            rows.append(["Environment", "; ".join("%s <code>%s</code>" % (k, e(v)) for k, v in sorted(env.items()))])
        conv, a = j.get("converted") or {}, j.get("adapt") or {}
        came = []
        if conv.get("from"):
            came.append("docked with %s and converted on %s" % (e(conv["from"]), e(conv.get("date") or "")))
        if a.get("folder"):
            came.append("receptor from ADapt, session <code>%s</code> save %s" % (e(a["folder"]), e(a.get("save"))))
        if j.get("receptor_prepared"):
            came.append("receptor came prepared (%s)" % e(prot.get("by") or "before this run"))
        if str(j.get("target_file") or "").strip():
            came.append("maps from the target file <code>%s</code>" % e(j["target_file"]))
        rows.append(["Provenance", "; ".join(came) if came else "everything here was prepared and docked by this run"])
        return rows


    def appendix(self):
        j = self.j
        out = []
        _s, _t, log = prep_steps(self.rundir)
        if log:
            out.append({"heading": "The receptor's preparation log", "blocks": [{"kind": "pre", "text": log}]})
        tl = os.path.join(self.rundir, "target", "receptor.log")
        if os.path.isfile(tl):
            try:
                txt = open(tl, encoding="utf-8", errors="replace").read()
            except OSError:
                txt = ""
            if txt.strip():
                out.append({"heading": "The maps: AGFR's log", "blocks": [{"kind": "pre", "text": txt[:20000]}]})
        rows = [["Receptor file", e(j.get("receptor") or "")], ["Chains", e(", ".join(j.get("chains") or []))],
                ["pH", e(j.get("ph"))], ["Protonation", e(jobmod.protonation_of(j))],
                ["Minimisation", "%s steps, heavy atoms restrained at %s kcal/mol/\u00c5\u00b2" % (
                    e(j.get("minimise_steps")), e(j.get("restraint"))) if j.get("minimise", True) else "none"],
                ["Site from", e(j.get("site_mode") or ("residues" if j.get("hotspots") else ""))],
                ["Engine", e(self.engine)], ["Replicas", "{:,}".format(int(j.get("replicas") or 0))],
                ["Steps per replica", "{:,}".format(int(j.get("steps") or 0))], ["Poses kept", e(j.get("poses"))],
                ["Ring", e(", ".join(sorted({jobmod.ring_of(j, p) for p in j.get("peptides") or []})) or "none")],
                ["Random seed", e(j.get("seed") if j.get("seed") not in (None, "") else "ADCP's own")],
                ["MM-GBSA", ("%s, %s, %s poses, %s" % (
                    GB_LABELS.get(j.get("mmgbsa_engine") or "amber", j.get("mmgbsa_engine")),
                    e(j.get("mmgbsa_model") or "obc2"), e(j.get("mmgbsa_poses")), e(j.get("mmgbsa_aggregate") or "mean"))
                    if j.get("mmgbsa") else "no")],
                ["Ranked by", e(self.key)],
                ["Jobs", e(", ".join("%s %s" % kv for kv in (j.get("jobs") or {}).items()))],
                ["Created", e(j.get("created"))], ["Run folder", "<code>%s</code>" % e(self.rundir)]]
        out.append({"heading": "Settings", "blocks": [{"kind": "facts", "rows": rows}]})
        out.append({"heading": "Software and provenance", "blocks": [{"kind": "facts", "rows": self.software()}]})
        conv = j.get("converted") or {}
        if conv:
            out.append({"heading": "Conversion from ADCP 1.0", "blocks": [{"kind": "facts", "rows": [
                ["From", e(conv.get("from"))], ["Converted", e(conv.get("date"))], ["By", e(conv.get("by"))],
                ["What was done", e(conv.get("note"))]]}]})
        a = j.get("adapt") or {}
        if a.get("folder"):
            out.append({"heading": "The ADapt session", "blocks": [{"kind": "facts", "rows": [
                [k, e(v)] for k, v in (("Session folder", a.get("folder")), ("Save", a.get("save")),
                                       ("Assembly tab", a.get("tab")), ("PDB", a.get("pdbid")),
                                       ("Input", a.get("input")), ("Biomolecule", a.get("biomolecule")),
                                       ("reduce", a.get("reduce")),
                                       ("Hetero groups kept", ", ".join(a.get("hetero") or [])),
                                       ("ADapt", "%s %s" % (a.get("adapt_version") or "", a.get("adapt_commit") or "")))
                if v not in (None, "", [])]}]})
        return out


_CONTACTS = {}


def _q(**kw):
    import urllib.parse
    return urllib.parse.urlencode({k: v for k, v in kw.items() if v is not None})


def _ordinal(n):
    if n is None:
        return ""
    names = {1: "first", 2: "second", 3: "third", 4: "fourth", 5: "fifth", 6: "sixth", 7: "seventh",
             8: "eighth", 9: "ninth", 10: "tenth"}
    return names.get(int(n), "%dth" % int(n))


def _few(items, word):
    items = list(items)
    if not items:
        return "no %ss" % word
    if len(items) <= 6:
        return "%s %s" % (word + ("s" if len(items) > 1 else ""), and_list(items))
    return "%d %ss (%s and others)" % (len(items), word, ", ".join(items[:4]))


def _pockets(rundir):
    try:
        txt = open(os.path.join(rundir, "target", "receptor.log"), errors="replace").read()
    except OSError:
        return 0
    m = re.search(r"found (\d+) pocket", txt)
    return int(m.group(1)) if m else 0


# ---- the views kept with the run ------------------------------------------------
def _folder(rundir):
    return os.path.join(rundir, FOLDER)


def load_items(rundir):
    try:
        with open(os.path.join(_folder(rundir), "items.json")) as fh:
            got = json.load(fh)
        return got if isinstance(got, list) else []
    except (OSError, ValueError):
        return []


def save_items(rundir, items):
    folder = _folder(rundir)
    os.makedirs(folder, exist_ok=True)
    tmp = os.path.join(folder, ".items.json.tmp")
    with open(tmp, "w") as fh:
        json.dump(items, fh, indent=1)
    os.replace(tmp, os.path.join(folder, "items.json"))


def add_item(rundir, body):
    """A view from the viewer: its picture and the pose it shows."""
    import base64
    data = str(body.get("image") or "")
    if not data.startswith("data:image/png;base64,"):
        return {"error": "the view came without a picture"}
    try:
        png = base64.b64decode(data.split(",", 1)[1], validate=True)
    except ValueError:
        return {"error": "the picture could not be read"}
    if not png.startswith(b"\x89PNG") or len(png) > IMAGE_MAX:
        return {"error": "the picture is not a PNG, or is too large"}
    seq = str(body.get("peptide") or "")
    if seq and not os.path.isdir(os.path.join(rundir, "docking", jobmod.key(seq))):
        return {"error": "this run has no peptide %s" % seq}
    try:
        rank = int(body.get("rank") or 1)
    except (TypeError, ValueError):
        return {"error": "no such pose"}
    items = load_items(rundir)
    n = 1 + max([int(x.get("id", 0)) for x in items] or [0])
    folder = _folder(rundir)
    os.makedirs(folder, exist_ok=True)
    name = "view_%03d.png" % n
    with open(os.path.join(folder, name), "wb") as fh:
        fh.write(png)
    item = {"id": n, "peptide": seq, "rank": rank, "image": FOLDER + "/" + name,
            "viewer": "molstar" if body.get("viewer") == "molstar" else "3dmol",
            "own": str(body.get("own") or "")[:2000], "added": time.strftime("%Y-%m-%d %H:%M")}
    items.append(item)
    save_items(rundir, items)
    return {"item": item, "count": len(items)}


def edit_item(rundir, body):
    """Change a view's own words, move it, or take it out."""
    try:
        n = int(body.get("id"))
    except (TypeError, ValueError):
        return {"error": "no such view"}
    items = load_items(rundir)
    hit = next((x for x in items if int(x.get("id", 0)) == n), None)
    if hit is None:
        return {"error": "no such view"}
    if body.get("remove"):
        items.remove(hit)
        try:
            os.remove(os.path.join(rundir, hit["image"]))
        except OSError:
            pass
    else:
        if "own" in body:
            hit["own"] = str(body.get("own") or "")[:2000]
        if body.get("move") in (-1, 1):
            k = items.index(hit)
            m = k + int(body["move"])
            if 0 <= m < len(items):
                items[k], items[m] = items[m], items[k]
    save_items(rundir, items)
    return {"ok": True}


def item_image(rundir, n):
    hit = next((x for x in load_items(rundir) if str(x.get("id")) == str(n)), None)
    if not hit:
        return None
    path = os.path.realpath(os.path.join(rundir, hit.get("image") or ""))
    if not path.startswith(os.path.realpath(_folder(rundir)) + os.sep):
        return None
    try:
        with open(path, "rb") as fh:
            return fh.read()
    except OSError:
        return None
