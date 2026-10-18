#!/usr/bin/env python3
"""The complex MM-GBSA builds keeps the peptide apart from the receptor,
whatever chains the receptor is on.

By Aadhil Haq

The peptide was always stamped chain B. With a receptor on chain B as well
(3Q47 and 1ELW in the benchmark), the reduce step read the receptor back as
part of the peptide, the TER between them was lost, and tleap bonded the
receptor's last residue to the peptide's first, 40 A away, and stopped:
every pose of every such run had no MM-GBSA. The interface panel measured
the same receptor as all peptide and found no interface.

What is checked: the peptide gets a chain the receptor does not use; a TER
stands between receptor and peptide, between receptor chains, and at a break
in a receptor chain, and nowhere inside the peptide; reduce's output is read
back as the peptide alone; a prepared complex left by the old assembly is
made again rather than reused; the split into receptor and ligand for
MMPBSA.py; the ring and D residue builds; the panel; and, with Amber, that
tleap builds the topology for a receptor on chain B, on chains A and B, and
on a blank chain, protonated by reduce or by PDB2PQR, linear, as a D peptide
and as both kinds of ring.

    python3 tests/test_chains.py

Exit 0 passes, 77 skips (no Amber here, after everything else passed), and
anything else fails.
"""
import json
import os
import shutil
import stat
import sys
import tempfile
import types

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import job as jobmod, pipeline                     # noqa: E402
from adcp_dock._vendor import adcp, mmgbsa, poseviewer, prep      # noqa: E402
from adcp_dock._vendor import receptor as recmod                  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)
        print("  FAIL " + why)


TMP = tempfile.mkdtemp(prefix="adcp_chains_")
EX = os.path.join(ROOT, "adcp_dock", "example_run")
REC = os.path.join(EX, "receptor", "1YCR.pdb")
P53 = "ETFSDLWKLLPEN"
POSE = adcp.read_pose(os.path.join(EX, "docking", P53), 1)


def relabel(path, out, mapping, drop=()):
    """A copy of a PDB file with its chains renamed, residues in `drop`
    (chain, number) left out, and only ATOM records."""
    with open(out, "w") as fh:
        for l in open(path):
            if not l.startswith("ATOM"):
                continue
            ch = l[21]
            if ch not in mapping or (ch, int(l[22:26])) in drop:
                continue
            fh.write(l[:21] + mapping[ch] + l[22:])
        fh.write("END\n")
    return out


def rec_of(path):
    """The receptor as MM-GBSA is handed it: heavy atoms with their chains."""
    return types.SimpleNamespace(atoms=recmod._read_pdb(path, ""))


def atoms_of(path):
    return [l for l in open(path) if l.startswith(("ATOM", "HETATM"))]


def blocks(path):
    """The file as tleap reads it: runs of atoms between TER records, each
    as [(chain, residue number)]."""
    out, cur = [], []
    for l in open(path):
        if l.startswith("TER") or l.startswith("END"):
            if cur:
                out.append(cur)
            cur = []
        elif l.startswith(("ATOM", "HETATM")):
            key = (l[21], int(l[22:26]))
            if not cur or cur[-1] != key:
                cur.append(key)
    if cur:
        out.append(cur)
    return out


mol = poseviewer.transplant(P53, POSE, ph=7.4)
check(mol is not None, "the example's p53 pose could not be laid on a "
      "topology, so nothing below means anything")

# --- the letter -------------------------------------------------------------
pc = getattr(mmgbsa, "peptide_chain", None)
check(pc is not None and pc({"A"}) == "B" and pc({"A", "C"}) == "B",
      "a receptor that does not use B no longer gives the peptide B")
check(pc is not None and pc({"B"}) == "C" and pc({"A", "B"}) == "C"
      and pc({"B", "C", "D"}) == "E" and pc({"_"}) == "B",
      "a receptor on B does not give the peptide a letter of its own")

# --- a receptor on chain B: the complex as written --------------------------
recB_pdb = relabel(REC, os.path.join(TMP, "recB.pdb"), {"A": "B"})
recB = rec_of(recB_pdb)
n_rec_atoms = len(recB.atoms)
n_rec_res = len({(a.chain, a.resnum) for a in recB.atoms})
wd = os.path.join(TMP, "write")
os.makedirs(wd)
_r, _l, com = mmgbsa._write_pdbs(recB, mol, wd, "x")
at = atoms_of(com)
rec_chains = {l[21] for l in at[:n_rec_atoms]}
pep_chains = {l[21] for l in at[n_rec_atoms:]}
check(len(at) > n_rec_atoms and rec_chains == {"B"} and len(pep_chains) == 1
      and not (pep_chains & rec_chains),
      "the peptide shares a chain with the receptor in the complex written "
      "(receptor %s, peptide %s)" % (sorted(rec_chains), sorted(pep_chains)))
b = blocks(com)
check(len(b) == 2 and len(b[0]) == n_rec_res and len(b[1]) == len(P53),
      "the complex written is not the receptor and the peptide, each ended "
      "by a TER (%s)" % [len(x) for x in b])


# --- reduce's output read back as the peptide alone -------------------------
fake_reduce = os.path.join(TMP, "reduce")
open(fake_reduce, "w").write('#!/bin/sh\n'
                             'for a in "$@"; do f="$a"; done\n'
                             'echo "USER  MOD reduce.3.23 stand in"\n'
                             'cat "$f"\n')
os.chmod(fake_reduce, os.stat(fake_reduce).st_mode | stat.S_IEXEC)
red = {"method": "reduce", "states": {}, "reduce": fake_reduce}
clean = mmgbsa.prep_complex(com, wd, "x", 7.4, log=lambda *a: None,
                            protonation=red)
at2 = atoms_of(clean)
pep_chain = at2[-1][21]
pep_res = {(l[21], l[22:27]) for l in at2 if l[21] == pep_chain}
check(len(at2) == len(at) and len(pep_res) == len(P53),
      "after reduce the peptide's chain holds %d residue(s), not the %d of "
      "the peptide, and the complex %d atoms, not %d"
      % (len(pep_res), len(P53), len(at2), len(at)))
b = blocks(clean)
check(len(b) == 2 and len(b[0]) == n_rec_res and len(b[1]) == len(P53)
      and b[0][-1][0] != b[1][0][0],
      "after reduce no TER stands between receptor and peptide, so tleap "
      "would bond the receptor's last residue to the peptide's first (%s)"
      % [len(x) for x in b])
rm, lm, nr, nl, _asl = mmgbsa.masks_by_chain(clean, pep_chain)
check(nr == n_rec_res and nl == len(P53),
      "the split for MMPBSA.py gives the receptor %d residue(s) and the "
      "ligand %d" % (nr, nl))

# A prepared complex the old assembly left is made again, not reused.
stale = os.path.join(TMP, "stale")
os.makedirs(stale)
_r, _l, com_s = mmgbsa._write_pdbs(recB, mol, stale, "s")
merged = [l[:21] + "B" + l[22:] for l in atoms_of(com_s)]
open(os.path.join(stale, "s_com_forleap.pdb"), "w").writelines(
    merged + ["TER\n", "TER\n", "END\n"])
clean_s = mmgbsa.prep_complex(com_s, stale, "s", 7.4, log=lambda *a: None,
                              protonation=red)
b = blocks(clean_s)
check(len(b) == 2 and len(b[1]) == len(P53),
      "a prepared complex left by the old assembly, receptor and peptide "
      "merged on chain B, was reused (%s)" % [len(x) for x in b])
again = mmgbsa.prep_complex(com_s, stale, "s", 7.4, log=lambda *a: None,
                            protonation=red)
check(again == clean_s and blocks(again) == blocks(clean_s),
      "a complex prepared from this very file is not kept for a resumed run")

# --- several receptor chains, a blank chain, and a broken chain ------------
two = os.path.join(TMP, "recAB.pdb")
lines = [l for l in open(REC) if l.startswith("ATOM") and l[21] == "A"]
with open(two, "w") as fh:
    fh.writelines(lines)
    for l in lines:          # a second copy, 60 A away, as chain B
        fh.write(l[:21] + "B" + l[22:30] + "%8.3f" % (float(l[30:38]) + 60)
                 + l[38:])
    fh.write("END\n")
recAB = rec_of(two)
wd2 = os.path.join(TMP, "two")
os.makedirs(wd2)
_r, _l, com2 = mmgbsa._write_pdbs(recAB, mol, wd2, "y")
clean2 = mmgbsa.prep_complex(com2, wd2, "y", 7.4, log=lambda *a: None,
                             protonation=red)
b = blocks(clean2)
chains = [{c for c, _n in x} for x in b]
check(len(b) == 3 and chains[0] == {"A"} and chains[1] == {"B"}
      and len(chains[2]) == 1 and not (chains[2] & {"A", "B"})
      and len(b[2]) == len(P53),
      "a receptor on chains A and B and the peptide are not three molecules "
      "with a letter each (%s)" % chains)
blank = relabel(REC, os.path.join(TMP, "recblank.pdb"), {"A": " "})
recX = rec_of(blank)
wd3 = os.path.join(TMP, "blank")
os.makedirs(wd3)
_r, _l, com3 = mmgbsa._write_pdbs(recX, mol, wd3, "z")
b = blocks(com3)
check(len(b) == 2 and b[1][0][0] not in {c for c, _n in b[0]},
      "a receptor with no chain letter and the peptide are not apart (%s)"
      % [sorted({c for c, _n in x}) for x in b])
gap = relabel(REC, os.path.join(TMP, "recgap.pdb"), {"A": "B"},
              drop={("A", n) for n in (60, 61, 62)})
recG = rec_of(gap)
wd4 = os.path.join(TMP, "gap")
os.makedirs(wd4)
_r, _l, com4 = mmgbsa._write_pdbs(recG, mol, wd4, "g")
b = blocks(com4)
check(len(b) == 3 and b[0][-1][1] == 59 and b[1][0][1] == 63
      and len(b[2]) == len(P53),
      "a receptor chain broken between 59 and 63 is not two pieces, each "
      "ended by a TER (%s)" % [(x[0], x[-1]) for x in b])
# Running it again on what it wrote changes nothing.
before = open(com4).read()
mmgbsa.separate_chains(com4, b[2][0][0]) if hasattr(
    mmgbsa, "separate_chains") else None
check(open(com4).read() == before, "separating the chains twice changes the "
      "file")

# --- a ring and a D peptide find the peptide by its chain ------------------
ring_dir = os.path.join(TMP, "ring")
os.makedirs(ring_dir)
shutil.copy(clean, os.path.join(ring_dir, "r_com_forleap.pdb"))
text = mmgbsa._ring_leap(os.path.join(ring_dir, "r_com_forleap.pdb"),
                         {"kind": "backbone", "pairs": [[1, len(P53)]]},
                         ring_dir, "r", "mbondi2")
rec_side = atoms_of(os.path.join(ring_dir, "r_rec_forleap.pdb"))
pep_side = atoms_of(os.path.join(ring_dir, "r_pep_forleap.pdb"))
check(len(rec_side) == n_rec_atoms
      and len({l[22:27] for l in pep_side}) == len(P53)
      and "bond pep.1.N pep.%d.C" % len(P53) in text,
      "a ring's receptor and peptide files are not the receptor and the "
      "peptide (%d receptor atoms of %d, %d peptide residues)"
      % (len(rec_side), n_rec_atoms, len({l[22:27] for l in pep_side})))
pl = mmgbsa._peptide_lines(open(clean).readlines(), mirror=True)
check(len({l[22:27] for l in pl if l.startswith("ATOM")}) == len(P53),
      "the reflected peptide a D residue's hydrogens are built from is not "
      "the peptide alone")

# --- the interface panel, with no MM-GBSA ----------------------------------
run = os.path.join(TMP, "panelrun")
pdir = os.path.join(run, "docking", P53)
os.makedirs(pdir)
shutil.copy(os.path.join(EX, "docking", P53, "out_ranked_1.pdb"), pdir)
j = jobmod.default(receptor=recB_pdb, chains=["B"], hotspots="B:54,B:58,B:62",
                   mmgbsa_poses=1, protonation="propka")
out = {"modes": [{"rank": 1, "score": -20.0}]}
pipeline._panel_geometry(run, j, P53, out, pdir, log=lambda *a: None)
bsa = (out.get("panel") or {}).get("bsa")
check(bsa is not None and bsa > 300,
      "the interface of p53 on MDM2 relabelled as chain B measures %s A2 of "
      "buried surface; on chain A it is over 1000" % bsa)

# --- with Amber: tleap builds every one of them -----------------------------
home = mmgbsa.amber_home()
if not home or not prep.find_tool("tleap"):
    shutil.rmtree(TMP, ignore_errors=True)
    if fails:
        print("test_chains FAILED (%d of %d)" % (len(fails), ran))
        sys.exit(1)
    print("skip  AmberTools is not here; the rest was checked (%d checks)" % ran)
    sys.exit(77)


def built(tag, rec, m, protonation=None, ph=7.4):
    """_topologies on one complex; (topo, the leap log's text)."""
    rec = types.SimpleNamespace(atoms=rec.atoms)
    if protonation:
        rec.protonation = protonation
    d = os.path.join(TMP, "leap_" + tag)
    os.makedirs(d)
    logs = []
    topo = mmgbsa._topologies(rec, m, d, tag, log=logs.append, prep_ph=ph)
    lg = os.path.join(d, "leap.log")
    return topo, (open(lg).read() if os.path.isfile(lg) else ""), logs


real_reduce = None
for root in (os.environ.get("ADCP_TEST_TOOLS", ""),
             os.environ.get("ADCPSUITE_HOME", "")):
    for cand in (os.path.join(root, "bin", "reduce"),
                 os.path.join(root, "adcpsuite", "bin", "reduce")):
        if root and os.path.isfile(cand):
            real_reduce = real_reduce or cand
red_real = dict(red, reduce=real_reduce or fake_reduce)
cases = [("b_reduce", recB, mol, red_real),
         ("ab_reduce", recAB, mol, red_real),
         ("blank_reduce", recX, mol, red_real),
         ("gap_reduce", recG, mol, red_real)]
if prep.find_tool("pdb2pqr", "pdb2pqr30"):
    cases.append(("b_propka", recB, mol, None))
pose_d = dict(POSE, xyz=np.array(POSE["xyz"]) * np.array([-1.0, 1.0, 1.0]))
molD = poseviewer.transplant("&" + "&".join(P53), pose_d, ph=7.4)
if molD is not None:
    from adcp_dock import mmforms
    mmforms.mark(molD, "&" + "&".join(P53), None)
    cases.append(("b_dpeptide", recB, molD, red_real))
for tag, rec, m, protonation in cases:
    topo, leaplog, logs = built(tag, rec, m, protonation)
    nres = len({(a.chain, a.resnum) for a in rec.atoms})
    check(topo is not None,
          "%s: no topology was built (%s)" % (tag, [l for l in logs
                                                   if "skipped" in l][-1:]))
    check("There is a bond of" not in leaplog,
          "%s: tleap bonded across a gap: %s" % (tag, [
              l for l in leaplog.splitlines() if "bond of" in l][:2]))
    if topo is not None:
        check(topo.n_lig_res == len(P53) and topo.n_rec_res == nres,
              "%s: the topology splits into %d receptor and %d peptide "
              "residue(s), not %d and %d" % (tag, topo.n_rec_res,
                                             topo.n_lig_res, nres, len(P53)))
        rec_ch = {(a.chain or "A") for a in rec.atoms}
        lc = getattr(topo, "lig_chain", None)
        of = getattr(mmgbsa, "peptide_chain_of", lambda *a: None)
        check(lc is not None and lc not in rec_ch
              and of(topo.workdir, tag) == lc,
              "%s: the peptide's chain %s is the receptor's, or not recorded"
              % (tag, lc))

# Both rings, built by tleap for a peptide beside a receptor on chain B.
from adcp_dock._vendor import library                             # noqa: E402
from rdkit import Chem                                            # noqa: E402


def toy(seq, chain, first, shift):
    m = Chem.RemoveHs(library.build_3d(seq, capped=False))
    conf, out, n = m.GetConformer(), [], 0
    for a in m.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        info = a.GetMonomerInfo()
        nm = info.GetName().strip()
        n += 1
        out.append("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                   "          %2s\n" % (n, nm if len(nm) == 4 else " " + nm,
                                        info.GetResidueName().strip(), chain,
                                        info.GetResidueNumber() - 1 + first,
                                        p.x + shift, p.y, p.z, nm[0]))
    return out


for kind, pairs in (("backbone", [[1, 5]]), ("disulfide", [[1, 5]])):
    d = os.path.join(TMP, "toy_" + kind)
    os.makedirs(d)
    c = os.path.join(d, "t_com_forleap.pdb")
    open(c, "w").writelines(toy("AAAAA", "B", 1, 30.0) + ["TER\n"]
                            + toy("CAGWC", "C", 6, 0.0) + ["TER\n", "END\n"])
    leapin = mmgbsa._ring_leap(c, {"kind": kind, "pairs": pairs}, d, "t",
                               "mbondi2")
    open(os.path.join(d, "t.leap"), "w").write(leapin)
    import subprocess
    subprocess.run([os.path.join(home, "bin", "tleap"), "-f", "t.leap"],
                   cwd=d, capture_output=True, text=True,
                   env=dict(os.environ, AMBERHOME=home))
    check(os.path.isfile(os.path.join(d, "t_com.parm7")),
          "tleap does not build the %s ring beside a receptor on chain B"
          % kind)

shutil.rmtree(TMP, ignore_errors=True)
if fails:
    print("test_chains FAILED (%d of %d)" % (len(fails), ran))
    sys.exit(1)
print("ok  the peptide has a chain of its own and TERs where the molecules "
      "end, for every receptor chain, %d checks" % ran)
