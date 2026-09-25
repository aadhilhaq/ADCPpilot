"""Rewrite the references section of CITATION.cff from the tool's own
citation list (adcp_dock/cite.py), with full metadata from Crossref.

By Aadhil Haq

    python scripts/cff_references.py

Run it after adding a reference to cite.py; tests/test_api.py checks that
every reference there is in CITATION.cff.
"""
import json
import os
import re
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from adcp_dock import cite  # noqa: E402

NOTES = {"adcp": "The docking engine. Scripps Research Institute.",
         "adfr": "ADFRsuite: agfr and prepare_receptor. Scripps Research Institute.",
         "autogrid": "AutoGrid4, which agfr uses to compute the affinity maps.",
         "autosite": "AutoSite, run by agfr on every grid target.",
         "pdb2pqr": "Receptor protonation at a chosen pH.",
         "propka": "The pKa prediction PDB2PQR uses.",
         "pdbfixer": "Rebuilding missing atoms (PDBFixer is part of OpenMM).",
         "openmm": "Restrained receptor minimization; ADCP 1.1's pose minimization.",
         "amber": "AmberTools (tleap, sander), when MM-GBSA is run.",
         "ff14sb": "The ff14SB force field, when MM-GBSA is run.",
         "mmpbsa": "MMPBSA.py, when MM-GBSA is run.",
         "gbmodel": "The GB-OBC(II) solvation model, when MM-GBSA is run.",
         "3dmol": "The 3Dmol.js structure viewer.",
         "molstar": "The Mol* structure viewer.",
         "numpy": "Array arithmetic.",
         "rdkit": "Peptide topologies."}


def q(text):
    return json.dumps(re.sub(r"\s+", " ", text).strip(), ensure_ascii=False)


def crossref(doi):
    req = urllib.request.Request("https://api.crossref.org/works/" + doi,
                                 headers={"User-Agent": "cff_references.py"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)["message"]


def clean(t):
    t = re.sub(r"\s*<(sub|sup)>\s*([^<]*?)\s*</\1>", r"\2", t or "")
    return re.sub(r"\s+([:;,])", r"\1", re.sub(r"\s+", " ", re.sub(r"<[^>]+>", "", t))).strip()


out, seen = ["references:"], set()
for key, when, did, ref, doi in cite.REFERENCES:
    if doi in seen:
        continue
    seen.add(doi)
    if doi.startswith("10.5281/zenodo"):
        out += ["  - type: software", "    title: %s" % q("RDKit: Open-source cheminformatics"),
                "    authors:", "      - name: %s" % q("The RDKit developers"),
                "    url: \"https://www.rdkit.org\"", "    doi: %s" % q(doi),
                "    notes: %s" % q(NOTES.get(key, did))]
        continue
    m = crossref(doi)
    out += ["  - type: article", "    title: %s" % q(clean(m["title"][0])), "    authors:"]
    for a in m.get("author", []):
        if a.get("family"):
            out.append("      - family-names: %s" % q(a["family"]))
            if a.get("given"):
                out.append("        given-names: %s" % q(a["given"]))
        elif a.get("name"):
            out.append("      - name: %s" % q(a["name"]))
    out.append("    journal: %s" % q(clean(m["container-title"][0])))
    if m.get("volume"):
        out.append("    volume: %s" % q(m["volume"]))
    if m.get("issue"):
        out.append("    issue: %s" % q(m["issue"]))
    page = m.get("page") or m.get("article-number") or ""
    if page:
        sp, _, ep = page.partition("-")
        out.append("    start: %s" % q(sp))
        if ep:
            out.append("    end: %s" % q(ep))
    out.append("    year: %d" % m["issued"]["date-parts"][0][0])
    out.append("    doi: %s" % q(m["DOI"].lower()))
    notes = NOTES.get(key, did)
    if key == "pdbfixer":
        notes = "OpenMM, including PDBFixer: receptor minimization, missing atoms and ADCP 1.1's pose minimization."
    out.append("    notes: %s" % q(notes))
p = os.path.join(ROOT, "CITATION.cff")
s = open(p, encoding="utf-8").read()
s = s[:s.index("references:")] + "\n".join(out) + "\n"
open(p, "w", encoding="utf-8").write(s)
print("wrote %d references" % len(seen))
