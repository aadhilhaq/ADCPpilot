"""The report of a run: what it says, and the routes that serve it.

By Aadhil Haq

The report page is laid out by reportkit; what it says is runreport's, built
from the run (GET /api/report). What is checked, on the bundled example run and
on copies of it changed into the other kinds of run:

  - the document is the contract's: its sections, its blocks, its options;
  - every caption states a result from the run's own numbers, never a title,
    and the numbers are the ones the run wrote (results.json);
  - a pose's contacts are classified exactly as the viewer classifies them
    (app.js classifyContacts, run under Node where Node is found);
  - a run converted from ADCP 1.0, one docked on a Mac with ADCP 1.0, ring
    peptides, a receptor from ADapt, a target file, a run without MM-GBSA and
    an unfinished one each say what they are;
  - the options change what is shown;
  - views are added from the viewer, captioned, moved and taken out, and kept
    in the run's report folder; the page is served with its token only.

    python3 tests/test_report.py
"""
import html as _html
import json
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile
import threading
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zlib

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_pilot import example, mac10, page, runreport, server  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


TMP = tempfile.mkdtemp(prefix="adcp_report_")
os.environ["HOME"] = TMP
RUNS = os.path.join(TMP, "runs")
EX = example.install(RUNS)
app = server.App(RUNS)


def blocks_of(doc):
    """Every block in the document, with the heading it is under."""
    out = []

    def walk(sec):
        for b in sec.get("blocks") or []:
            out.append((sec.get("heading"), b))
        for c in sec.get("sections") or []:
            walk(c)
    for s in doc.get("sections") or []:
        walk(s)
    for v in doc.get("views") or []:
        for b in v.get("blocks") or []:
            out.append(("views", b))
    for a in doc.get("appendix") or []:
        for b in a.get("blocks") or []:
            out.append((a.get("heading"), b))
    return out


def text(s):
    return _html.unescape(re.sub(r"<[^>]+>", "", s or ""))


def first_sentence(cap):
    """Up to the first full stop that ends a sentence (not one in a number)."""
    m = re.search(r"\.(\s|$)", cap)
    return cap[:m.start()] if m else cap


#: What a caption may not open with: a title in place of a result.
TITLES = re.compile(r"^(Table|Figure|Interaction diagram|The peptides ranked|The top (poses|clusters)|"
                    r"Each peptide|Residue by residue|The interface panel|The receptor's preparation|"
                    r"Ranking of|Scores of)\b")


def well_formed(doc, where):
    """The contract's document: what reportkit can lay out."""
    for k in ("title", "subtitle", "byline", "options", "sections", "methods", "references", "appendix", "views"):
        check(k in doc, "%s: the report has no %r" % (where, k))
    names = {r.get("name") for r in ((json.load(open(os.path.join(doc["path"], "results.json"))).get("rows"))
                                     if os.path.isfile(os.path.join(doc["path"], "results.json")) else [])}
    for heading, b in blocks_of(doc):
        kind = b.get("kind")
        check(kind in ("p", "list", "table", "figure", "facts", "pre", "heading"),
              "%s, %s: a block of unknown kind %r" % (where, heading, kind))
        if kind == "table":
            check(all(c.get("key") and "label" in c for c in b.get("columns") or []),
                  "%s, %s: a column without a key or label" % (where, heading))
            check(isinstance(b.get("rows"), list), "%s, %s: a table without rows" % (where, heading))
        if kind == "figure":
            check(any(b.get(k) for k in ("svg", "png", "structure", "bars", "panels")),
                  "%s, %s: a figure with nothing to draw" % (where, heading))
            for p in b.get("panels") or []:
                check(any(p.get(k) for k in ("svg", "png", "structure", "bars")),
                      "%s, %s: a panel with nothing to draw" % (where, heading))
            if b.get("svg"):
                try:
                    ET.fromstring(b["svg"])
                    ok = True
                except ET.ParseError as exc:
                    ok = str(exc)
                check(ok is True, "%s, %s: the figure's SVG does not parse: %s" % (where, heading, ok))
        if kind in ("table", "figure") and "caption" in b:
            cap = text(b["caption"])
            check(cap.strip(), "%s, %s: an empty caption" % (where, heading))
            check(not TITLES.match(cap), "%s, %s: a caption opens with a title: %r" % (where, heading, cap[:80]))
            lead = first_sentence(cap)
            check(re.search(r"\d", lead) or any(n and n in lead for n in names),
                  "%s, %s: a caption's first sentence states no result: %r" % (where, heading, lead[:120]))
    json.dumps(doc)


# --- the example run, as it is ------------------------------------------------------
d = app.report_doc(EX, {})
check(not d.get("error"), "the example's report could not be written: %s" % d.get("error"))
res = json.load(open(os.path.join(EX, "results.json")))
rows = res["rows"]
if not d.get("error"):
    well_formed(d, "example")
    check(d["kicker"] == "ADCP docking report", "the report is not named an ADCP docking report")
    check(d["title"] == json.load(open(os.path.join(EX, "job.json")))["label"], "the title is not the run's label")
    check([o["key"] for o in d["options"]] == ["peptides", "poses", "diagrams", "vdw"],
          "the bar's options are %s" % [o["key"] for o in d["options"]])
    check([s["id"] for s in d["sections"]] == ["summary", "receptor", "ranking", "peptides"],
          "the sections are %s" % [s["id"] for s in d["sections"]])
    check(not d["notes"], "an ADCP 1.1 run carries a note above its report: %s" % d["notes"])
    lead = rows[0]
    summ = d["sections"][0]["blocks"][0]["html"]
    check("<b>%s</b>" % lead["name"] in summ and runreport.pm(lead["mmgbsa_top5"], lead["mmgbsa_sd"]) in summ,
          "the summary does not lead with %s at its own mean: %s" % (lead["name"], summ[:300]))
    rank = next(b for _h, b in blocks_of(d) if b.get("id") == "ranking")
    check(rank["caption"].startswith("%s ranks first by MM-GBSA, %s" % (
        lead["name"], runreport.pm(lead["mmgbsa_top5"], lead["mmgbsa_sd"]))),
        "the ranking's caption does not state the first peptide's own numbers: %s" % rank["caption"][:200])
    check([r.get("mmgbsa_best") for r in rank["rows"]] == [r.get("mmgbsa_best") for r in rows]
          and [r.get("adcp_best") for r in rank["rows"]] == [r.get("adcp_best") for r in rows],
          "the ranking's numbers are not results.json's")
    peps = d["sections"][3]["sections"]
    check(len(peps) == len(rows) and all(p.get("newpage") for p in peps),
          "each peptide does not have a page of its own: %d" % len(peps))
    ids = [b.get("id") for _h, b in blocks_of(d)]
    for want in ("prep", "site", "sitefig", "ranking", "rankchart", "panel"):
        check(want in ids, "the report has no %s block" % want)
    check(any(i and i.startswith("lid_") for i in ids), "no peptide has an interaction diagram")
    check(any(i and i.startswith("pres_") for i in ids), "no peptide has its residues one by one")
    figs = [b for _h, b in blocks_of(d) if b.get("kind") == "figure" and b.get("panels")]
    check(figs and all(len(f["panels"]) == 3 for f in figs), "the top three poses are not drawn by default")
    for f in figs:
        for p in f["panels"]:
            q = urllib.parse.parse_qs(urllib.parse.urlparse(p["structure"]["url"]).query)
            got, err = app.structure_text(q["dir"][0], q["peptide"][0], int(q["rank"][0]), q["what"][0])
            check(err is None and "ATOM" in (got or ""), "a pose picture asks for a structure the server cannot "
                  "give: %s" % err)
    # the methods, with their references numbered as cited
    m, refs = d["methods"], d["references"]
    nums = sorted({int(n) for g in re.findall(r"\[([\d, ]+)\]", m) for n in g.split(",")})
    check(nums == list(range(1, len(refs) + 1)), "the methods cite %s of %d references" % (nums, len(refs)))
    check(len({r["text"] for r in refs}) == len(refs), "a reference is listed twice")
    check("ADCP 1.1" in m and "reduce" in m and "AmberTools" in m, "the methods do not say what the run did")
    check([r["n"] for r in refs] == list(range(1, len(refs) + 1)), "the references are not numbered in order")
    heads = [a["heading"] for a in d["appendix"]]
    check("Settings" in heads and "The receptor's preparation log" in heads and "Software and provenance" in heads,
          "the appendix holds %s" % heads)
    soft = dict(next(a for a in d["appendix"] if a["heading"] == "Software and provenance")["blocks"][0]["rows"])
    check("ADCPpilot" in soft.get("Written by", "") and "ADCP 1.1" in soft.get("Docking", "")
          and "AGFR 1.1" in soft.get("Maps", "") and "prepared and docked by this run" in soft.get("Provenance", ""),
          "the appendix does not say what wrote the run and what it ran: %s" % soft)

# --- the contacts are the viewer's ------------------------------------------------
pdb, err = app.structure_text(EX, rows[0]["peptide"], 1, "complex")
check(err is None, "the example's first pose could not be read: %s" % err)
mine = runreport.classify(pdb or "")
check(mine["chain"] == "Z" and len(mine["pairs"]) > 10, "the first pose has %d contacts" % len(mine["pairs"]))
node = shutil.which("node") or "/sw/eb/sw/nodejs/20.13.1-GCCcore-13.3.0/bin/node"
if os.path.isfile(node) and pdb:
    appjs = open(os.path.join(ROOT, "adcp_pilot", "static", "app.js")).read()
    consts = [L for L in appjs.splitlines() if re.match(r"const (HBOND_A|POLAR|BASIC|ACIDIC|RING|SPECIFIC)\b", L)]
    start = appjs.index("function classifyContacts(pdb, pep){")
    end = appjs.index("\n}\n", start) + 3
    probe = os.path.join(TMP, "probe.js")
    with open(probe, "w") as fh:
        fh.write("\n".join(consts) + "\n" + appjs[start:end] + "\nconst fs = require('fs');\n"
                 "console.log(JSON.stringify(classifyContacts(fs.readFileSync(process.argv[2], 'utf8'), 'Z')));\n")
    with open(os.path.join(TMP, "complex.pdb"), "w") as fh:
        fh.write(pdb)
    env = dict(os.environ, LD_LIBRARY_PATH="/sw/eb/sw/GCCcore/13.3.0/lib64:" + os.environ.get("LD_LIBRARY_PATH", ""))
    r = subprocess.run([node, probe, os.path.join(TMP, "complex.pdb")], capture_output=True, text=True,
                       timeout=120, env=env)
    check(r.returncode == 0, "app.js's classifyContacts did not run under Node: %s" % r.stderr[-300:])
    if r.returncode == 0:
        js = json.loads(r.stdout)
        a = sorted((p["chain"], p["resi"], p["bresi"], p["kind"], round(p["d"], 4)) for p in js["pairs"])
        b = sorted((p["chain"], p["resi"], p["bresi"], p["kind"], round(p["d"], 4)) for p in mine["pairs"])
        check(a == b, "the report's contacts are not the viewer's: %d against %d pairs, first differences %s"
              % (len(b), len(a), [x for x in a if x not in b][:3] + [x for x in b if x not in a][:3]))
        check(sorted((x["chain"], x["resi"]) for x in js["res"]) == sorted((x["chain"], x["resi"]) for x in mine["receptor"]),
              "the report's receptor residues in contact are not the viewer's")
else:
    print("note  Node was not found, so the contacts were not held against app.js")
svg = runreport.lid_svg(mine["pairs"])
check(svg and ET.fromstring(svg) is not None and "Trp7" in svg, "the interaction diagram does not draw p53's Trp7")
check("hydrogen bond" in svg, "the diagram has no key")


# --- the other kinds of run -----------------------------------------------------
def variant(name, change_job=None, change_rows=None, change_result=None):
    """A copy of the example, changed into another kind of run."""
    runs = os.path.join(TMP, name)
    rd = example.install(runs)
    jp = os.path.join(rd, "job.json")
    j = json.load(open(jp))
    if change_job:
        change_job(j)
    json.dump(j, open(jp, "w"), indent=1)
    rp = os.path.join(rd, "results.json")
    r = json.load(open(rp))
    if change_rows:
        change_rows(r)
    json.dump(r, open(rp, "w"), indent=1)
    if change_result:
        change_result(rd)
    return rd, app.report_doc(rd, {})


def conv(j):
    j["converted"] = {"from": "ADCP 1.0", "date": "2026-10-18", "by": "convert_adcp10_runs.py",
                      "note": "Docked with ADCP 1.0; its poses and scores are 1.0's, unchanged."}


def conv_rows(r):
    for x in r["rows"]:
        x["engine"] = "ADCP 1.0"


rd, dv = variant("converted", conv, conv_rows)
check(not dv.get("error"), "a converted run's report: %s" % dv.get("error"))
if not dv.get("error"):
    well_formed(dv, "converted")
    check(dv["notes"] and dv["notes"][0]["tone"] == "warn" and "ADCP 1.0" in dv["notes"][0]["html"]
          and "not comparable" in dv["notes"][0]["html"], "a converted run does not say so above its report")
    check("docked with ADCP 1.0" in dv["subtitle"], "a converted run's subtitle names %s" % dv["subtitle"])
    check("converted on 2026-10-18" in dv["methods"], "a converted run's methods do not say it was converted")
    check("Conversion from ADCP 1.0" in [a["heading"] for a in dv["appendix"]], "no conversion in the appendix")
    soft = dict(next(a for a in dv["appendix"] if a["heading"] == "Software and provenance")["blocks"][0]["rows"])
    check("converted on 2026-10-18" in soft.get("Provenance", ""), "a converted run's provenance: %s" % soft)


def mac(j):
    j["adcp_engine"] = mac10.ENGINE


def mac_rows(r):
    for x in r["rows"]:
        x["engine"] = mac10.LABEL


rd, dv = variant("mac", mac, mac_rows)
check(not dv.get("error"), "a Mac run's report: %s" % dv.get("error"))
if not dv.get("error"):
    check(dv["notes"] and dv["notes"][0]["html"] == _html.escape(mac10.NOTE, quote=True),
          "a Mac run does not say so above its report: %s" % dv["notes"])
    check(mac10.LABEL in dv["subtitle"] and "ran on a Mac" in dv["methods"], "a Mac run's subtitle and methods")

PMI = "TSFAEYWNLLS"


def ring(j):
    for p in j["peptides"]:
        if p["sequence"] == PMI:
            p["cyclic"] = "backbone"


def ring_rows(r):
    for x in r["rows"]:
        if x["peptide"] == PMI:
            x["ring"] = "head to tail"


def ring_result(rd):
    p = os.path.join(rd, "results", PMI + ".json")
    x = json.load(open(p))
    x["cyclic"] = "backbone"
    x["mmgbsa_poses_skipped"] = [{"pose": 3, "why": "ring not closed: its ends are 6.1 \u00c5 apart", "gap": 6.1}]
    json.dump(x, open(p, "w"))


rd, dv = variant("ring", ring, ring_rows, ring_result)
check(not dv.get("error"), "a ring run's report: %s" % dv.get("error"))
if not dv.get("error"):
    well_formed(dv, "ring")
    rank = next(b for _h, b in blocks_of(dv) if b.get("id") == "ranking")
    check(any(c["key"] == "ring" for c in rank["columns"]), "the ranking has no ring column once a peptide is a ring")
    sec = next(s for s in dv["sections"][3]["sections"] if PMI in s["heading"])
    said = " ".join(b.get("html", "") for b in sec["blocks"] if b["kind"] == "p")
    check("as a ring head to tail" in said and "close the ring" in said and "pose 3" in said,
          "a ring peptide's page does not say how it was docked and which poses were open: %s" % said[:400])
    check("PMI was docked as a ring head to tail" in dv["methods"], "the methods do not name the ring peptide")


def adapt(j):
    j["adapt"] = {"folder": "/nowhere/adapt_session", "save": 2, "tab": "biomolecule 1", "pdbid": "1YCR",
                  "biomolecule": "1", "reduce": "reduce.3.23", "hetero": ["ZN"], "adapt_version": "ADapt 0.9",
                  "adapt_commit": "6dae90b"}


rd, dv = variant("adapt", adapt)
check(not dv.get("error"), "an ADapt run's report: %s" % dv.get("error"))
if not dv.get("error"):
    recp = " ".join(b.get("html", "") for b in dv["sections"][1]["blocks"] if b["kind"] == "p")
    check("prepared in ADapt" in recp and "save 2" in recp and "ZN" in recp,
          "the receptor's section does not say it came from ADapt: %s" % recp[:300])
    check("The ADapt session" in [a["heading"] for a in dv["appendix"]], "no ADapt session in the appendix")
    check("prepared in ADapt" in dv["methods"], "the methods do not say ADapt prepared the receptor")


def target(j):
    j["target_file"] = "/somewhere/receptor.trg"


rd, dv = variant("target", target)
if not dv.get("error"):
    recp = " ".join(b.get("html", "") for b in dv["sections"][1]["blocks"] if b["kind"] == "p")
    check("target file" in recp and "target file" in dv["methods"], "a run from a target file does not say so")
check(not dv.get("error"), "a target-file run's report: %s" % dv.get("error"))


def no_mm_job(j):
    j["mmgbsa"] = False


def no_mm_rows(r):
    r["ranked_by"] = "adcp_top5"
    for x in r["rows"]:
        for k in list(x):
            if k.startswith("mmgbsa") or k.startswith("mm_") or k in ("bsa", "panel_pose", "site_dg", "n_hbond"):
                x.pop(k)
    r["rows"].sort(key=lambda x: x["adcp_top5"])
    for i, x in enumerate(r["rows"], 1):
        x["rank"] = i


def no_mm_result(rd):
    for f in os.listdir(os.path.join(rd, "results")):
        p = os.path.join(rd, "results", f)
        x = json.load(open(p))
        for k in [k for k in x if k.startswith("mmgbsa") or k.startswith("panel")]:
            x.pop(k)
        json.dump(x, open(p, "w"))
    shutil.rmtree(os.path.join(rd, "mmgbsa"), ignore_errors=True)
    shutil.rmtree(os.path.join(rd, "metrics"), ignore_errors=True)


rd, dv = variant("adcponly", no_mm_job, no_mm_rows, no_mm_result)
check(not dv.get("error"), "a run without MM-GBSA: %s" % dv.get("error"))
if not dv.get("error"):
    well_formed(dv, "ADCP only")
    rank = next(b for _h, b in blocks_of(dv) if b.get("id") == "ranking")
    check("ranks first by ADCP score" in rank["caption"], "a run without MM-GBSA is not ranked by ADCP: %s"
          % rank["caption"][:120])
    check("Rescoring" not in dv["methods"] and "MM-GBSA" not in dv["subtitle"], "a run without MM-GBSA says it had it")


def unfinished_rows(r):
    r["rows"] = r["rows"][:2]


def unfinished_result(rd):
    os.remove(os.path.join(rd, "results", "RFMDYWEGL.json"))


rd, dv = variant("unfinished", None, unfinished_rows, unfinished_result)
check(not dv.get("error"), "an unfinished run's report: %s" % dv.get("error"))
if not dv.get("error"):
    summ = dv["sections"][0]["blocks"][0]["html"]
    check("12/1 core has not finished docking" in summ, "an unfinished run does not say what is missing: %s"
          % summ[-200:])

check("error" in app.report_doc(os.path.join(TMP, "no_such_run"), {}), "a folder that is no run has a report")
from adcp_pilot import job as _jobmod, __version__  # noqa: E402
_wb = _jobmod.default(written_by={"version": "0.0.1"})["written_by"]
check(_wb.get("version") == __version__ and _wb.get("app"), "a new job does not record the tool that wrote it: %s" % _wb)

# --- the options -----------------------------------------------------------------
d1 = app.report_doc(EX, {"peptides": "1", "poses": "1", "diagrams": "all", "vdw": "1"})
check(not d1.get("error"), "the report with options: %s" % d1.get("error"))
if not d1.get("error"):
    check(len(d1["sections"][3]["sections"]) == 1, "Peptides 1 does not show one peptide")
    figs = [b for _h, b in blocks_of(d1) if b.get("kind") == "figure" and b.get("panels")]
    check(figs and len(figs[0]["panels"]) == 1, "Poses drawn 1 does not draw one pose")
    lids = [b for _h, b in blocks_of(d1) if (b.get("id") or "").startswith("lid_")]
    check(lids and all("van der Waals contacts drawn too" in b["caption"] for b in lids),
          "van der Waals in diagrams does not say so in the caption")
    check({o["key"]: o["value"] for o in d1["options"]} == {"peptides": 1, "poses": 1, "diagrams": "all", "vdw": True},
          "the options are not echoed as used: %s" % {o["key"]: o["value"] for o in d1["options"]})
d2 = app.report_doc(EX, {"peptides": "999", "poses": "-4"})
check(not d2.get("error") and {o["key"]: o["value"] for o in d2["options"]}["poses"] == 1
      and len(d2["sections"][3]["sections"]) == len(rows), "options out of range are not held to their range")

# --- the routes, and the views added from the viewer -------------------------------
httpd = server.build(host="127.0.0.1", port=0, runs=RUNS)
PORT, TOKEN = httpd.server_address[1], httpd.app.token
threading.Thread(target=httpd.serve_forever, daemon=True).start()


def get(path, raw=False, token=True, **q):
    if token:
        q["t"] = TOKEN
    with urllib.request.urlopen("http://127.0.0.1:%d%s?%s" % (PORT, path, urllib.parse.urlencode(q)),
                                timeout=60) as r:
        body = r.read()
        return body if raw else json.loads(body.decode())


def post(path, body):
    req = urllib.request.Request("http://127.0.0.1:%d%s?t=%s" % (PORT, path, TOKEN), data=json.dumps(body).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode())


def png(w=4, h=3):
    raw = b"".join(b"\x00" + b"\xff\xff\xff" * w for _ in range(h))

    def chunk(t, data):
        return struct.pack(">I", len(data)) + t + data + struct.pack(">I", zlib.crc32(t + data) & 0xffffffff)
    return (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(raw)) + chunk(b"IEND", b""))


import base64  # noqa: E402
URI = "data:image/png;base64," + base64.b64encode(png()).decode()
try:
    get("/report", raw=True, token=False, rundir=EX)
    check(False, "the report page opened without the token")
except urllib.error.HTTPError as exc:
    check(exc.code == 403, "the report page without the token answered %d" % exc.code)
pg = get("/report", raw=True, rundir=EX).decode()
check("ReportKit.start(" in pg and TOKEN in pg, "the report page does not start reportkit with the token")
check(get("/static/reportkit.js", raw=True, token=False) == open(os.path.join(ROOT, "adcp_pilot", "static",
                                                                           "reportkit.js"), "rb").read(),
      "the server does not serve reportkit.js as it is")
doc = get("/api/report", rundir=EX, **{"opt.poses": "2"})
check(not doc.get("error") and doc["views"] == [] and {o["key"]: o["value"] for o in doc["options"]}["poses"] == 2,
      "GET /api/report with an option: %s" % doc.get("error"))
check("error" in post("/api/report/add", {"rundir": EX, "peptide": "NOSUCH", "rank": 1, "image": URI}),
      "a view of a peptide the run does not have is kept")
check("error" in post("/api/report/add", {"rundir": EX, "peptide": rows[0]["peptide"], "rank": 1,
                                          "image": "data:image/png;base64," + base64.b64encode(b"GIF89a").decode()}),
      "a picture that is no PNG is kept")
a1 = post("/api/report/add", {"rundir": EX, "peptide": rows[0]["peptide"], "rank": 2, "image": URI, "viewer": "molstar"})
check(a1.get("count") == 1 and os.path.isfile(os.path.join(EX, "report", "view_001.png")),
      "the view is not kept in the run's report folder: %s" % a1)
a2 = post("/api/report/add", {"rundir": EX, "peptide": rows[1]["peptide"], "rank": 1, "image": URI})
check(a2.get("count") == 2, "a second view is not kept: %s" % a2)
doc = get("/api/report", rundir=EX)
vs = doc.get("views") or []
check(len(vs) == 2 and vs[0]["caption"].startswith("Pose 2 of %s as drawn in the viewer (Mol*)" % rows[0]["name"])
      and "receptor residue" in vs[0]["caption"], "a view's caption does not say what it shows: %s"
      % (vs[0]["caption"] if vs else vs))
check(vs and vs[0]["blocks"] and vs[0]["blocks"][0].get("svg"), "a view has no interaction diagram beside it")
img = get("/api/report/image", raw=True, rundir=EX, id=1)
check(img.startswith(b"\x89PNG"), "the view's picture is not served")
e1 = post("/api/report/edit", {"rundir": EX, "id": 1, "own": "The helix lies in the cleft."})
check([v.get("own") for v in e1.get("views") or []][:1] == ["The helix lies in the cleft."],
      "a view's own words are not kept: %s" % e1)
e2 = post("/api/report/edit", {"rundir": EX, "id": 2, "move": -1, "opt": {"vdw": True}})
check([v["id"] for v in e2.get("views") or []] == [2, 1], "a view is not moved up: %s" % e2)
e3 = post("/api/report/edit", {"rundir": EX, "id": 1, "remove": True})
check([v["id"] for v in e3.get("views") or []] == [2] and not os.path.exists(os.path.join(EX, "report",
                                                                                          "view_001.png")),
      "a view taken out is not removed with its picture: %s" % e3)
check("error" in post("/api/report/edit", {"rundir": EX, "id": 99, "own": "x"}), "an edit of no view succeeded")
check("3Dmol" in get("/api/report", rundir=EX)["methods"], "the methods do not cite the viewer")
# the problem report keeps its own route
check(post("/api/report", {}).get("error") == "no run given", "POST /api/report is no longer the problem report")

# --- the page offers it ----------------------------------------------------------
appjs = open(os.path.join(ROOT, "adcp_pilot", "static", "app.js")).read()
check("reportLine(dir, rows)" in appjs and "Report (PDF)" in appjs and "/api/report/add" in appjs,
      "the run's page does not offer the report, or the viewer cannot add to it")
check("id=raddrep" in page.render("tok"), "the viewer's bar has no add to report")

httpd.shutdown()
shutil.rmtree(TMP, ignore_errors=True)
if fails:
    print("test_report FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the report of a run, its captions, its kinds of run and its views, %d checks" % ran)
