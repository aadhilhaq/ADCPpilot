"""Build adcp_dock/static/modifications.json, the list the Modifications
dialog offers, from ADCP 1.1's side-chain library and the PDB's Chemical
Component Dictionary.

    python scripts/modifications_table.py /path/to/adcpsuite

For each non-standard side chain in ADCP's Swiss library: its code, its name
and parent residue from the RCSB chemical component entry (where the code is
one), ADCP's coarse residue for it, and the code of its D form.
"""
import concurrent.futures as cf
import glob
import json
import os
import sys
import urllib.request

THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
         "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
         "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
         "TYR": "Y", "VAL": "V"}
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def chemcomp(code):
    try:
        with urllib.request.urlopen(
                "https://data.rcsb.org/rest/v1/core/chemcomp/" + code,
                timeout=30) as r:
            m = json.load(r)["chem_comp"]
        return m.get("name"), m.get("mon_nstd_parent_comp_id") or []
    except Exception:                                           # noqa: BLE001
        return None, []


rot = glob.glob(os.path.join(sys.argv[1], "lib", "python3*", "site-packages",
                             "ADCP", "data", "rotamers"))[0]
coarse = {}
for line in open(os.path.join(rot, "swiss.lib")):
    if line.startswith("rotamer "):
        p = line.split()
        coarse[p[1]] = p[2].upper()
rows = []
for line in open(os.path.join(rot, "SWISS_ROTAMER_LIST.dat")):
    if line.startswith("#") or not line.split():
        continue
    p = line.split()
    rows.append({"code": p[0], "d": p[3] if p[2] == "Y" else None,
                 "coarse": coarse.get(p[0])})
with cf.ThreadPoolExecutor(8) as ex:
    meta = dict(zip([r["code"] for r in rows],
                    ex.map(chemcomp, [r["code"] for r in rows])))
out = []
for r in rows:
    name, parents = meta[r["code"]]
    parent = next((THREE[x] for x in parents if x in THREE), None)
    # A side chain ADCP gives no coarse residue cannot be written in its
    # sequence syntax, so it is not offered.
    if not r["coarse"] or r["coarse"] == "O":
        continue
    if name and name.isupper():
        name = name.lower()          # the dictionary writes some names in capitals
    out.append({"code": r["code"], "name": name or r["code"],
                "parent": parent or r["coarse"], "coarse": r["coarse"],
                "d": r["d"]})
out.sort(key=lambda x: (x["parent"], x["name"].lower()))
dest = os.path.join(ROOT, "adcp_dock", "static", "modifications.json")
json.dump({"source": "ADCP 1.1 Swiss side-chain library; names and parents "
                     "from the RCSB Chemical Component Dictionary",
           "modifications": out}, open(dest, "w"), indent=1)
print("%d modifications, %d named from the RCSB" % (
    len(out), sum(1 for x in out if x["name"] != x["code"])))
