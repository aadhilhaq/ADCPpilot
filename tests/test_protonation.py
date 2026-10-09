"""Protonation: reduce by default, PROPKA at a pH when chosen, and what an old
job file means. And the histidine names MM-GBSA hands tleap under reduce.

By Aadhil Haq

What is checked: a new job asks for reduce and an old one, with no
protonation key, still means PDB2PQR with PROPKA, through job.json, the
command line, the preparation's reuse key, the citations and the batch
notes; the walltime prices both; reduce, ADFRsuite's, is run with ADFR's
flags and charges a chain that starts at residue 25 as a terminus there, so
the minimisation accepts it; a residue after a chain break is given the
amide hydrogen its cap calls for; and the complex MM-GBSA builds takes the
receptor's histidines from where reduce put the protons (HID, HIE, HIP),
keeps an amide reduce turned over as it is in the coordinates, and has
reduce choose the peptide's.

    python3 tests/test_protonation.py
"""
import contextlib
import io
import json
import math
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import cite, cli, job as jobmod, pipeline, sizing  # noqa: E402
from adcp_dock import tools as toolmod  # noqa: E402
from adcp_dock._vendor import mmgbsa, prepare  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


TMP = tempfile.mkdtemp(prefix="adcp_protonation_")
os.environ["HOME"] = TMP
REC = os.path.join(ROOT, "adcp_dock", "example_run", "receptor", "1YCR.pdb")
ROOTS = [r for r in ("/scratch/group/p.agr250019.000/tools/docking",)
         if os.path.isdir(r)]

# --- the default, and what an old job file means ------------------------------
check(jobmod.default()["protonation"] == "reduce", "a new run asks for reduce")
check(jobmod.protonation_of({}) == "propka",
      "a job file without the setting was protonated with PDB2PQR and PROPKA, "
      "and still is")
check(jobmod.protonation_of({"protonation": "reduce"}) == "reduce"
      and jobmod.protonation_of({"protonation": "propka"}) == "propka",
      "each choice reads as itself")
check(jobmod.protonation_of({"protonation": "propka",
                             "adapt": {"folder": "/x"}}) == "reduce",
      "a receptor from ADapt was protonated by reduce there, whatever the box said")
old_example = json.load(open(os.path.join(ROOT, "adcp_dock", "example_run", "job.json")))
check("protonation" not in old_example and jobmod.protonation_of(old_example) == "propka",
      "the shipped example, made before the choice, reads as PDB2PQR with PROPKA")
bad = jobmod.default(receptor=REC, hotspots="A:96", protonation="guess",
                     peptides=[{"name": "p", "sequence": "ETFSDLWKLLPEN"}])
check(any("Protonation must be reduce or propka" in p for p in jobmod.validate(bad)),
      "an unknown protonation is refused before anything is submitted")

# --- the reuse key: an old preparation is still recognised --------------------
k_old = prepare._prep_key(REC, 7.4, "AMBER", "A", 0.0, False, False, False, True,
                          True, 10.0)
k_prop = prepare._prep_key(REC, 7.4, "AMBER", "A", 0.0, False, False, False, True,
                           True, 10.0, protonation="propka")
k_red = prepare._prep_key(REC, 7.4, "AMBER", "A", 0.0, False, False, False, True,
                          True, 10.0, protonation="reduce")
k_red55 = prepare._prep_key(REC, 5.5, "AMBER", "A", 0.0, False, False, False, True,
                            True, 10.0, protonation="reduce")
check(k_old == k_prop and "protonation" not in k_prop,
      "a preparation made before the choice is recognised as PROPKA's, and reused")
check(k_red != k_prop and k_red.get("protonation") == "reduce",
      "a reduce preparation is not mistaken for a PROPKA one")
check(k_red == k_red55, "reduce takes no pH, so the pH does not force a new one")

# --- the command line, under the page's names -----------------------------------
def run_cli(*argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            rc = cli.main(list(argv))
        except SystemExit as e:
            rc = e.code
    return rc, out.getvalue()


for flag, want in (([], "reduce"), (["--protonation", "propka", "--ph", "5.5"], "propka")):
    rd = os.path.join(TMP, "cli_%s" % want)
    rc, out = run_cli("new", rd, "-r", REC, "-s", "A:96", "--chain", "A",
                      "-p", "ETFSDLWKLLPEN", *flag)
    j = jobmod.load(rd)
    check(j.get("protonation") == want,
          "adcp-dock new %s writes protonation %r (%r)" % (flag, want, j.get("protonation")))
check(jobmod.load(os.path.join(TMP, "cli_propka"))["ph"] == 5.5, "--ph reaches the job")
rc, out = run_cli("new", os.path.join(TMP, "cli_bad"), "-r", REC, "-s", "A:96",
                  "-p", "ETFSDLWKLLPEN", "--protonation", "maybe")
check(rc == 2, "an unknown --protonation is refused by the command line")
helptext = io.StringIO()
with contextlib.redirect_stdout(helptext):
    try:
        cli.main(["new", "--help"])
    except SystemExit:
        pass
check("--protonation" in helptext.getvalue() and "--adapt-session" in helptext.getvalue(),
      "the command line offers the protonation and ADapt's receptor")

# --- what each run says it did ---------------------------------------------------
j_red = jobmod.default(receptor=REC, hotspots="A:96", chains=["A"], mmgbsa=True,
                       peptides=[{"name": "p", "sequence": "ETFSDLWKLLPEN"}])
j_old = dict(j_red)
del j_old["protonation"]
notes_red = "\n".join(pipeline.batch_notes(j_red, "prepare") + pipeline.batch_notes(j_red, "dock"))
notes_old = "\n".join(pipeline.batch_notes(j_old, "prepare") + pipeline.batch_notes(j_old, "dock"))
check("protonation: reduce" in notes_red and "reduce protonates the peptide" in notes_red,
      "the batch notes say reduce, for the receptor and for MM-GBSA's complexes")
check("PDB2PQR with PROPKA at pH 7.4" in notes_old,
      "an old job's batch notes say PDB2PQR with PROPKA at its pH")
k_r, k_o = cite.used_by(j_red), cite.used_by(j_old)
check("reduce" in k_r and "pdb2pqr" not in k_r and "propka" not in k_r,
      "a reduce run cites reduce and not PDB2PQR or PROPKA (%s)" % sorted(k_r))
check({"pdb2pqr", "propka"} <= k_o and "reduce" not in k_o,
      "an old run still cites PDB2PQR and PROPKA (%s)" % sorted(k_o))
fell = dict(j_red, receptor_info={"protonation": {"method": "propka"}})
check({"pdb2pqr", "propka"} <= cite.used_by(fell),
      "a reduce run that fell back to PDB2PQR cites what actually ran")
ad = dict(j_red, adapt={"folder": "/x", "save": 1})
k_a = cite.used_by(ad)
check("reduce" in k_a and "pdbfixer" not in k_a and "openmm" not in k_a,
      "a receptor from ADapt cites reduce and none of this package's preparation")
check(any(k == "reduce" and "10.1006/jmbi.1998.2401" == d for k, _w, _x, _r, d in cite.REFERENCES),
      "reduce's paper is Word et al. 1999")

# --- the walltime prices both, and neither at zero ----------------------------------
base = sizing.prep_minutes(5000, 0, protonation="propka")
red = sizing.prep_minutes(5000, 0, protonation="reduce")
check(sizing.PROTONATE_PER_1K_ATOMS["reduce"] > 0 and sizing.PROTONATE_PER_1K_ATOMS["propka"] > 0
      and red > sizing.PREP_FLOOR_MINUTES and base > sizing.PREP_FLOOR_MINUTES,
      "both protonations are priced, neither at zero (%s, %s)" % (base, red))
check(sizing.prep_minutes(5000, 0, prepared=False) < base,
      "a receptor prepared in ADapt is copied, so it costs only the maps")
plan_r = sizing.plan([13], 50, 2500000, atoms=5000, options={"protonation": "reduce"})
plan_p = sizing.plan([13], 50, 2500000, atoms=5000, options={"protonation": "propka"})
check(plan_r["prep_minutes"] != plan_p["prep_minutes"],
      "the page's estimate follows the choice")

# --- the amide hydrogen a cap calls for -------------------------------------------------
gap = os.path.join(TMP, "gap.pdb")
open(gap, "w").write(
    "ATOM      1  N   ALA A  10       0.000   0.000   0.000  1.00  0.00           N\n"
    "ATOM      2  CA  ALA A  10       1.458   0.000   0.000  1.00  0.00           C\n"
    "ATOM      3  C   ALA A  10       2.009   1.420   0.000  1.00  0.00           C\n"
    "ATOM      4  O   ALA A  10       1.251   2.390   0.000  1.00  0.00           O\n"
    "ATOM      5  N   GLY A  14      12.000   3.000   0.000  1.00  0.00           N\n"
    "ATOM      6  CA  GLY A  14      13.400   3.300   0.300  1.00  0.00           C\n"
    "ATOM      7  C   GLY A  14      14.200   2.100   0.800  1.00  0.00           C\n"
    "ATOM      8  O   GLY A  14      13.700   0.980   0.900  1.00  0.00           O\n"
    "END\n")
g = prepare.backbone_gaps(gap)
for amide, want in ((False, 0), (True, 1)):
    out = prepare.cap_breaks(gap, os.path.join(TMP, "capped_%s.pdb" % amide), g, [],
                             amide_h=amide)
    hs = [l for l in open(out) if l[17:20] == "GLY" and l[12:16].strip() == "H"]
    check(len(hs) == want, "after a break, %s amide hydrogen with amide_h=%s"
          % (want, amide))
    if hs:
        n = [l for l in open(out) if l[17:20] == "GLY" and l[12:16].strip() == "N"][0]
        xyz = lambda l: [float(l[30:38]), float(l[38:46]), float(l[46:54])]  # noqa: E731
        check(abs(math.dist(xyz(hs[0]), xyz(n)) - 1.01) < 0.02,
              "and it sits 1.01 A from the nitrogen")

# --- the histidine names tleap is given -------------------------------------------------
def atom(serial, name, res, chain, num, x, y, z, el):
    nm = (" " + name).ljust(4) if len(name) < 4 else name
    return ("ATOM  %5d %s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, nm, res, chain, num, x, y, z, el))


ring = [("CG", 0, 0, 0, "C"), ("ND1", 1.2, 0.6, 0, "N"), ("CD2", -0.4, 1.3, 0, "C"),
        ("CE1", 1.0, 1.9, 0, "C"), ("NE2", 0.4, 2.3, 0, "N")]
prepared_lines, n = [], 0
for num, hs in ((73, ("HD1",)), (96, ("HE2",)), (120, ("HD1", "HE2"))):
    for nm, x, y, z, el in ring:
        n += 1
        prepared_lines.append(atom(n, nm, "HIS", "A", num, x + num, y, z, el))
    for h in hs:
        n += 1
        prepared_lines.append(atom(n, h, "HIS", "A", num, num + 2.0, 0.5, 0, "H"))
# A cysteine in a disulphide has no thiol hydrogen; a free one has.
for num, hs in ((30, ()), (31, ("HG",))):
    for nm, el in (("CB", "C"), ("SG", "S")):
        n += 1
        prepared_lines.append(atom(n, nm, "CYS", "A", num, num, 5.0, 0, el))
    for h in hs:
        n += 1
        prepared_lines.append(atom(n, h, "CYS", "A", num, num, 6.0, 0, "H"))
# An asparagine reduce turned over: its OD1 and ND2 as reduce left them.
for nm, x, el in (("CG", 0.0, "C"), ("OD1", 1.2, "O"), ("ND2", -1.1, "N")):
    n += 1
    prepared_lines.append(atom(n, nm, "ASN", "A", 79, 79 + x, 9.0, 0, el))
n += 1
prepared_lines.append(atom(n, "HD21", "ASN", "A", 79, 77.5, 9.5, 0, "H"))
prepared = os.path.join(TMP, "prepared.pdb")
open(prepared, "w").write("".join(prepared_lines) + "END\n")

states = mmgbsa.residue_states(prepared)
check(states.get(("A", 73)) == "HID" and states.get(("A", 96)) == "HIE"
      and states.get(("A", 120)) == "HIP",
      "HID, HIE and HIP from where the protons are (%s)" % states)
check(states.get(("A", 30)) == "CYX" and ("A", 31) not in states,
      "a cysteine without its thiol hydrogen is CYX, one with it stays CYS")
check(("A", 79) not in states, "an asparagine keeps its name")
heavy = os.path.join(TMP, "heavy.pdb")
open(heavy, "w").write("".join(l for l in prepared_lines if l[76:78].strip() != "H"))
check(mmgbsa.residue_states(heavy) == {},
      "a file with no hydrogens says nothing about any residue's state")

# The complex as _write_pdbs writes it: receptor heavy atoms named HIS, the
# peptide as chain B: here one whole histidine, MDM2's His96 moved clear.
com_lines = [l for l in prepared_lines if l[76:78].strip() != "H"]
pep = []
for i, l in enumerate(x for x in open(REC) if x.startswith("ATOM")
                      and x[21] == "A" and x[22:26].strip() == "96"):
    pep.append(atom(900 + i, l[12:16].strip(), "HIS", "B", 200,
                    float(l[30:38]) + 40, float(l[38:46]), float(l[46:54]),
                    l[76:78].strip() or l[12:16].strip()[0]))
com = os.path.join(TMP, "p_com.pdb")
open(com, "w").write("".join(com_lines) + "TER\n" + "".join(pep) + "TER\nEND\n")
wd = os.path.join(TMP, "mm")
os.makedirs(wd)
logged = []
clean = mmgbsa.prep_complex(com, wd, "p", 7.4, log=logged.append,
                            protonation={"method": "reduce", "states": states,
                                         "reduce": None})
names = {}
for l in open(clean):
    if l.startswith(("ATOM", "HETATM")):
        names[(l[21], int(l[22:26]))] = l[17:20]
check(names.get(("A", 73)) == "HID" and names.get(("A", 96)) == "HIE"
      and names.get(("A", 120)) == "HIP" and names.get(("A", 30)) == "CYX"
      and names.get(("A", 31)) == "CYS",
      "tleap is given the receptor's states as reduce placed them (%s)"
      % {k: v for k, v in names.items() if k[0] == "A"})
od1 = [l for l in open(clean) if l[17:20] == "ASN" and l[12:16].strip() == "OD1"][0]
check(abs(float(od1[30:38]) - 80.2) < 1e-3,
      "an amide reduce turned over stays turned in the coordinates tleap reads")
check(not any(l[12:16].strip().startswith("H") for l in open(clean)
              if l.startswith(("ATOM", "HETATM"))),
      "every hydrogen goes, and tleap builds them to match the names")
check(any("reduce was not found" in str(x) for x in logged),
      "without reduce the peptide keeps tleap's standard states, and the log says so")
legacy = mmgbsa.prep_complex.__code__.co_varnames
check("protonation" in legacy, "prep_complex takes the receptor's protonation")

# The pipeline attaches the states only under reduce.
class _Rec:
    pass


r1, r2 = _Rec(), _Rec()
pipeline._carry_protonation(r1, j_old, prepared, log=lambda *a: None)
pipeline._carry_protonation(r2, j_red, prepared, log=lambda *a: None)
check(not hasattr(r1, "protonation"),
      "a run protonated with PDB2PQR is rescored as it always was")
check(getattr(r2, "protonation", {}).get("states") == states,
      "a reduce run hands MM-GBSA the prepared receptor's states")
r3 = _Rec()
pipeline._carry_protonation(r3, dict(j_red, receptor_info={"protonation": {"method": "propka"}}),
                            prepared, log=lambda *a: None)
check(not hasattr(r3, "protonation"),
      "a reduce run that fell back to PDB2PQR is rescored the PDB2PQR way")

# --- reduce itself, where ADFRsuite has it --------------------------------------------
exe = toolmod.discover(roots=ROOTS).adfr_program("reduce")
if not exe:
    print("skip  reduce is not installed here; the rest was checked (%d checks)" % ran)
    if fails:
        for f in fails:
            print("  FAIL " + f)
        sys.exit(1)
    sys.exit(77)

check(os.path.basename(os.path.dirname(exe)) == "bin"
      and os.path.isfile(os.path.join(os.path.dirname(exe), "prepare_receptor"))
      or "adcpsuite" in exe,
      "reduce is ADFRsuite's own, beside prepare_receptor (%s)" % exe)
check(prepare.reduce_version(exe).startswith("reduce."),
      "its version is read (%r)" % prepare.reduce_version(exe))

# Chain A of 1YCR, tidied as the preparation tidies it.
wd2 = os.path.join(TMP, "prep")
os.makedirs(wd2)
st = []
lines = prepare.as_pdb_lines(REC, wd2)
lines = prepare.keep_chain(prepare.drop_metadata(lines, st), "A", st)
lines = prepare.drop_hetatms(prepare.drop_waters(lines, st), st, keep_ions=False)
tidy = prepare._write(os.path.join(wd2, "tidied.pdb"), lines)
out = prepare.protonate_reduce(tidy, os.path.join(wd2, "protonated.pdb"), st, exe=exe)
check(out is not None, "reduce protonated MDM2 (%s)" % [str(s) for s in st])
if out:
    first = [l[12:16].strip() for l in open(out)
             if l.startswith("ATOM") and l[21] == "A" and int(l[22:26]) == 25]
    check({"H1", "H2", "H3"} <= set(first),
          "a chain that starts at residue 25 is charged as a terminus there (%s)" % first)
    his = prepare.his_states(out)
    check(his == {("A", "73"): "HIE", ("A", "96"): "HIE"},
          "reduce gives MDM2's His73 and His96 their epsilon protons (%s)" % his)
    line = [str(s) for s in st if "protonation: done" in str(s)]
    check(line and "-BUILD" in line[0] and "HIE x2" in line[0],
          "the log line names reduce, its flags and the histidines (%s)" % line)
    check(any("USER  MOD" in l for l in open(out)),
          "reduce's own record of its decisions stays at the top of the file")
    check(os.path.isfile(os.path.join(wd2, "protonated.reduce.txt")),
          "reduce's command and what it said are kept beside it")

    # The whole preparation with reduce, minimisation included, where OpenMM
    # and PDBFixer are installed.
    try:
        import openmm  # noqa: F401
        import pdbfixer  # noqa: F401
        have_omm = True
    except ImportError:
        have_omm = False
    if have_omm:
        wd3 = os.path.join(TMP, "full")
        log = []
        got = prepare.run(REC, wd3, chain="A", protonation="reduce", reduce_exe=exe,
                          log=log.append)
        text = open(os.path.join(wd3, "preparation.log")).read()
        check(got and "protonation: done; reduce." in text,
              "the preparation protonates with reduce when asked")
        check("minimisation: done" in text,
              "and the minimisation accepts reduce's hydrogens (%s)"
              % [l for l in text.splitlines() if "minimisation" in l])
        check(prepare.his_states(got) == {("A", "73"): "HIE", ("A", "96"): "HIE"},
              "the histidines survive the minimisation as reduce placed them")
        bare = os.path.join(TMP, "bare")
        got2 = prepare.run(REC, bare, chain="A", protonation="reduce", reduce_exe=exe,
                           minimise_structure=False, rebuild_sidechains=False,
                           log=lambda *a: None)
        text2 = open(os.path.join(bare, "preparation.log")).read()
        check("left 0 missing side-chain atom(s) out" in text2
              and "side chains:" not in text2,
              "with reduce and no minimisation nothing is modelled, so nothing is "
              "stripped (%s)" % [l for l in text2.splitlines() if "missing" in l])
        broken = prepare.run(REC, os.path.join(TMP, "noreduce"), chain="A",
                             protonation="reduce", reduce_exe=os.path.join(TMP, "none"),
                             minimise_structure=False, log=lambda *a: None)
        text3 = open(os.path.join(TMP, "noreduce", "preparation.log")).read()
        check("reduce: skipped" in text3 and "protonation fallback: done" in text3,
              "without reduce the receptor is protonated by PDB2PQR, and the log says so")
    else:
        print("note  OpenMM or PDBFixer is not here; the full preparation was not run")

    # The peptide's histidine in a complex is chosen by reduce.
    wd4 = os.path.join(TMP, "mm2")
    os.makedirs(wd4)
    clean2 = mmgbsa.prep_complex(com, wd4, "q", 7.4, log=lambda *a: None,
                                 protonation={"method": "reduce", "states": states,
                                              "reduce": exe})
    pep_names = {l[17:20] for l in open(clean2) if l.startswith("ATOM") and l[21] == "B"}
    check(pep_names and pep_names <= {"HID", "HIE", "HIP"},
          "reduce chooses the peptide's histidine in the complex (%s)" % pep_names)
    check(os.path.isfile(os.path.join(wd4, "q_com_reduce.pdb")),
          "and what reduce wrote for the complex is kept")

if fails:
    print("test_protonation FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  reduce by default, PROPKA when chosen, old jobs kept, His names for "
      "tleap, %d checks" % ran)
