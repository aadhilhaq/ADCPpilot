"""A receptor that came prepared: docked as it is, its states kept, refused without hydrogens.

By Aadhil Haq

What is checked: --receptor-prepared reaches job.json with the pH and the chains; a receptor so
marked with no hydrogens is refused when the run is made and again by the prepare stage; the
prepared file's chains are copied line for line (names, numbers, hydrogens), its caps taken
whatever record they were written as and its waters left out; a prepared.json beside it is read
for who prepared it and at what pH; nothing of this package's preparation runs; prepare_receptor
gives the Amber names (HID, HIE, HIP, CYX, ASH, GLH, LYN) and the ACE and NME caps the AutoDock
types and charges it gives the standard residues; the PDBQT and the record keep the states and
the pH given; and MM-GBSA keeps the receptor's states while protonating the peptide at the run's
pH. The receptor is MDM2, chain A of 1YCR, prepared by HOPE Labs' preparation card (HOPE-MD's
Prepare, PDB2PQR with PROPKA at pH 6.5): tests/data/1ycr_A_prepared.pdb.

    python3 tests/test_receptor_prepared.py
"""
import contextlib
import io
import json
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import cite, cli, job as jobmod, pipeline, prepared, server  # noqa: E402
from adcp_dock import tools as toolmod  # noqa: E402
from adcp_dock._vendor import mmgbsa, prepare, receptor  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


TMP = tempfile.mkdtemp(prefix="adcp_rprep_")
os.environ["HOME"] = TMP
PREPARED = os.path.join(HERE, "data", "1ycr_A_prepared.pdb")
CRYSTAL = os.path.join(ROOT, "adcp_dock", "example_run", "receptor", "1YCR.pdb")
SITE = "A:54,A:58,A:62,A:67,A:72,A:93,A:96,A:99"
ROOTS = [r for r in (os.environ.get("ADCP_TEST_TOOLS", ""),
                     os.environ.get("ADCPSUITE_HOME", ""))
         if r and os.path.isdir(r)]
STD = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "CYX": "CYS", "ASH": "ASP", "GLH": "GLU",
       "LYN": "LYS"}


def residues(path, chains=None):
    """{(chain, number with insertion code): name} of the ATOM records."""
    out = {}
    for l in open(path, errors="replace"):
        if l.startswith(("ATOM", "HETATM")) and (chains is None or l[21] in chains):
            out.setdefault((l[21], l[22:27].strip()), l[17:20].strip())
    return out


def tally(names):
    got = {}
    for v in names.values():
        if v in prepared.STATE_NAMES or v in prepared.CAPS:
            got[v] = got.get(v, 0) + 1
    return got


def without_h(src, dst):
    with open(src) as fh, open(dst, "w") as out:
        for l in fh:
            if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() == "H":
                continue
            out.write(l)
    return dst


def renamed(src, dst, table):
    with open(src) as fh, open(dst, "w") as out:
        for l in fh:
            if l.startswith(("ATOM", "HETATM")) and l[17:20] in table:
                l = l[:17] + table[l[17:20]] + l[20:]
            out.write(l)
    return dst


def run_cli(argv):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        try:
            rc = cli.main(argv)
        except SystemExit as exc:
            rc = exc.code
    return rc, buf.getvalue()


BEFORE = residues(PREPARED)
check(tally(BEFORE) == {"HID": 2}, "the test receptor is the one described: %s" % tally(BEFORE))

# --- the setting, and the command line ------------------------------------------------------
check(jobmod.default()["receptor_prepared"] is False, "a new run prepares its receptor")
rd = os.path.join(TMP, "run")
rc, said = run_cli(["new", rd, "-r", PREPARED, "-s", SITE, "--chain", "A", "--ph", "6.5",
                    "--receptor-prepared", "-p", "ETFSDLWKLLPE"])
j = json.load(open(jobmod.path_of(rd))) if os.path.isfile(jobmod.path_of(rd)) else {}
check(rc == 0 and j.get("receptor_prepared") is True and j.get("ph") == 6.5
      and j.get("chains") == ["A"],
      "--receptor-prepared reaches job.json with the pH and the chain: rc %s, %s"
      % (rc, {k: j.get(k) for k in ("receptor_prepared", "ph", "chains")}))
bare = without_h(PREPARED, os.path.join(TMP, "bare.pdb"))
rc2, said2 = run_cli(["new", os.path.join(TMP, "run_bare"), "-r", bare, "-s", SITE,
                      "--chain", "A", "--receptor-prepared", "-p", "ETFSDLWKLLPE"])
check(rc2 == 1 and "carry no hydrogens" in said2 and "--receptor-prepared" in said2,
      "a receptor marked as prepared with no hydrogens is refused when the run is made: rc %s, %r"
      % (rc2, said2[-300:]))
rc3, _s3 = run_cli(["new", os.path.join(TMP, "run_plain"), "-r", bare, "-s", SITE,
                    "--chain", "A", "-p", "ETFSDLWKLLPE"])
check(rc3 == 0, "while the same file not marked as prepared is a run like any other: rc %s" % rc3)
check(any("ADapt" in p for p in jobmod.prepared_problems(
          {"receptor_prepared": True, "receptor": PREPARED, "adapt": {"folder": "/x"}})),
      "a receptor from ADapt cannot also be one that came prepared")
check(prepared.problems(os.path.join(TMP, "x.cif")) == [],
      "a file that is not there is left to the check that says so")
cif = os.path.join(TMP, "prepared.cif")
shutil.copy(PREPARED, cif)
check(any("read as a PDB" in p for p in prepared.problems(cif)),
      "a prepared receptor is read as PDB, and an mmCIF is refused with the reason")

# --- taking it in: line for line ------------------------------------------------------------
with_extra = os.path.join(TMP, "extra.pdb")
caps_cut = os.path.join(TMP, "cut.pdb")
with open(PREPARED) as fh, open(caps_cut, "w") as out:
    for l in fh:
        if l.startswith("ATOM") and 45 <= int(l[22:26]) <= 47:
            continue
        out.write(l)
capped = os.path.join(TMP, "capped.pdb")
prepare.cap_breaks(caps_cut, capped, [(("A", 44, " "), ("A", 48, " "), 0.0)], [],
                   log=lambda *a: None)
with open(capped) as fh, open(with_extra, "w") as out:
    for l in fh:
        if l.startswith("ATOM") and l[17:20] in ("ACE", "NME"):
            l = "HETATM" + l[6:]
        if l.startswith("END"):
            out.write("HETATM 9001  O   HOH A 901      30.000  30.000  30.000  1.00  0.00"
                      "           O\n")
            out.write("HETATM 9002 ZN    ZN A 902      31.000  30.000  30.000  1.00  0.00"
                      "          ZN\n")
        out.write(l)
took_to = os.path.join(TMP, "took.pdb")
took = prepared.take(with_extra, ["A"], took_to)
after = residues(took_to)
check(tally(after) == {"HID": 2, "ACE": 1, "NME": 1},
      "the states and the caps are what the file had: %s" % tally(after))
check(all(l.startswith("ATOM") for l in open(took_to)
          if l[17:20] in ("ACE", "NME")),
      "a cap written as HETATM is taken in as ATOM, which the receptor's readers read")
check(took["waters_dropped"] == 1 and took["groups_dropped"] == ["ZN A:902"]
      and "HOH" not in open(took_to).read(),
      "the water is left out, and the zinc with it unless hetero groups are kept: %s"
      % took["groups_dropped"])
kept = prepared.take(with_extra, ["A"], os.path.join(TMP, "kept.pdb"), keep_hetatms=True)
check(kept["hetero_kept"] == 1, "with hetero groups kept the zinc stays: %s" % kept)
plain = prepared.take(PREPARED, ["A"], os.path.join(TMP, "plain.pdb"))
src_atoms = [l[:66] for l in open(PREPARED) if l.startswith("ATOM")]
out_atoms = [l[:66] for l in open(os.path.join(TMP, "plain.pdb")) if l.startswith("ATOM")]
check(src_atoms == out_atoms and plain["hydrogens"] == sum(
          1 for l in open(PREPARED) if l.startswith("ATOM") and l[76:78].strip() == "H"),
      "every ATOM line is copied as it was, hydrogens and all: %d of %d"
      % (len(out_atoms), len(src_atoms)))
try:
    prepared.take(bare, ["A"], os.path.join(TMP, "never.pdb"))
    refused = ""
except prepared.NotPrepared as exc:
    refused = str(exc)
check("carry no hydrogens" in refused, "take refuses a receptor without hydrogens: %r" % refused)

# --- who prepared it, from the record beside it ----------------------------------------------
card = os.path.join(TMP, "card")
os.makedirs(os.path.join(card, "mdm2"))
mine = os.path.join(card, "mdm2", "prepared.pdb")
shutil.copy(PREPARED, mine)
with open(os.path.join(card, "prepared.json"), "w") as fh:
    json.dump({"prepared": 1, "ph": 6.5, "prep": {"protonation": "propka"},
               "structures": [{"name": "mdm2", "path": mine, "molecules": [
                   {"protonation": {"method": "propka", "ph": 6.5,
                                    "states": {"HID": 2}}}]}]}, fh)
said = prepared.provenance(mine)
check(said.get("ph") == 6.5 and said.get("method") == "propka" and "HOPE-MD" in said.get("by", ""),
      "the preparation card's record says who prepared it and at what pH: %s" % said)
check(prepared.provenance(PREPARED).get("record") is None,
      "a file with no record beside it says nothing about who prepared it")
other = os.path.join(TMP, "card2")
os.makedirs(os.path.join(other, "x"))
shutil.copy(PREPARED, os.path.join(other, "x", "prepared.pdb"))
with open(os.path.join(other, "prepared.json"), "w") as fh:
    json.dump({"ph": 5.0, "structures": [{"path": "/somewhere/else.pdb"}]}, fh)
check(prepared.provenance(os.path.join(other, "x", "prepared.pdb")).get("ph") is None,
      "a record that lists other files is not read as this one's")

# --- what is written about it -----------------------------------------------------------------
notes = "\n".join(pipeline.batch_notes(dict(j, mmgbsa=True), "prepare")
                  + pipeline.batch_notes(dict(j, mmgbsa=True), "dock"))
check("came prepared" in notes and "pH 6.5" in notes and "PROPKA protonates the peptide" in notes,
      "the job script says the receptor came prepared and how the peptide is protonated:\n%s"
      % notes)
keys = cite.used_by(dict(j, mmgbsa=False))
check(not keys & {"pdbfixer", "openmm", "pdb2pqr", "propka", "reduce"},
      "a run that prepared nothing cites no preparation program: %s"
      % sorted(keys & {"pdbfixer", "openmm", "pdb2pqr", "propka", "reduce"}))
check({"pdb2pqr", "propka"} <= cite.used_by(dict(j, mmgbsa=True)),
      "MM-GBSA protonates the peptide with PDB2PQR and PROPKA, and cites them")
check(pipeline.prep_time(j) <= pipeline.prep_time(dict(j, receptor_prepared=False)),
      "the prepare job books no preparation: %s against %s"
      % (pipeline.prep_time(j), pipeline.prep_time(dict(j, receptor_prepared=False))))

# --- MM-GBSA: the receptor's states kept, the peptide at the run's pH ------------------------
states = mmgbsa.residue_states(PREPARED, keep_named=True)
check(sorted(states.values()) == ["HID", "HID"],
      "the states read from the prepared receptor are its own names: %s" % states)
check(mmgbsa.residue_states(PREPARED) == {},
      "read as before, a residue named HID gives no state, which is why its name is kept")
heavy = [l for l in open(PREPARED) if l.startswith("ATOM") and l[76:78].strip() != "H"]
pep = [l[:21] + "B" + l[22:] for l in open(CRYSTAL)
       if l.startswith("ATOM") and l[21] == "B" and l[76:78].strip() != "H"]
com = os.path.join(TMP, "com.pdb")
open(com, "w").write("".join(heavy) + "TER\n" + "".join(pep) + "TER\nEND\n")
if prepare.find_tool("pdb2pqr", "pdb2pqr30"):
    wd = os.path.join(TMP, "mm")
    os.makedirs(wd)
    clean = mmgbsa.prep_complex(com, wd, "p", 3.0, log=lambda *a: None,
                                protonation={"method": "prepared", "states": states, "ph": 3.0,
                                             "source": PREPARED})
    got = residues(clean)
    rec_names = {k: v for k, v in got.items() if k[0] == "A"}
    pep_names = {k: v for k, v in got.items() if k[0] == "B"}
    check(rec_names == {k: v for k, v in BEFORE.items() if k[0] == "A"},
          "at pH 3 the receptor keeps every name it came prepared with: %s"
          % {k: (v, rec_names.get(k)) for k, v in BEFORE.items() if rec_names.get(k) != v})
    check({"GLH", "ASH"} & set(pep_names.values()),
          "while the peptide is protonated at the run's pH: %s" % sorted(set(pep_names.values())))
    check(not any(l[76:78].strip() == "H" for l in open(clean) if l.startswith("ATOM")),
          "and every hydrogen goes, for tleap to build from the names")
    prime = os.path.join(wd, "p_com_prepped.pdb")
    hd1 = {(l[21], l[22:26].strip()) for l in open(prime) if l.startswith("ATOM")
           and l[21] == "A" and l[17:20] == "HIS" and l[12:16].strip() == "HD1"}
    he2 = {(l[21], l[22:26].strip()) for l in open(prime) if l.startswith("ATOM")
           and l[21] == "A" and l[17:20] == "HIS" and l[12:16].strip() == "HE2"}
    check(len(hd1) == 2 and not he2 and any(l[21] == "B" for l in open(prime)
                                            if l.startswith("ATOM")),
          "Prime is given the receptor's own hydrogens in standard names, each HID with "
          "its proton on ND1 alone, beside the peptide: HD1 %s, HE2 %s" % (sorted(hd1), sorted(he2)))
    stamp = json.load(open(os.path.join(wd, "p_com_forleap.json")))
    check(stamp["how"]["method"] == "prepared" and stamp["how"]["ph"] == 3.0,
          "the kept complex is stamped with how it was protonated: %s" % stamp["how"])
    # A capped receptor through tleap: the caps' heavy atoms named as ff14SB names them.
    cheavy = [l for l in open(capped) if l.startswith("ATOM") and l[76:78].strip() != "H"]
    ccom = os.path.join(TMP, "ccom.pdb")
    open(ccom, "w").write("".join(cheavy) + "TER\n" + "".join(pep) + "TER\nEND\n")
    wd2 = os.path.join(TMP, "mmcap")
    os.makedirs(wd2)
    cstates = mmgbsa.residue_states(capped, keep_named=True)
    cclean = mmgbsa.prep_complex(ccom, wd2, "c", 6.5, log=lambda *a: None,
                                 protonation={"method": "prepared", "states": cstates,
                                              "source": capped})
    caps = {(l[17:20], l[12:16].strip()) for l in open(cclean)
            if l.startswith("ATOM") and l[17:20] in ("ACE", "NME")}
    check(caps == {("ACE", "C"), ("ACE", "O"), ("ACE", "CH3"), ("NME", "N"), ("NME", "C")},
          "the caps reach tleap in ff14SB's names, NME's methyl as C: %s" % sorted(caps))
    home = mmgbsa.amber_home()
    if home:
        import subprocess
        leap = os.path.join(wd2, "c.leap")
        open(leap, "w").write(mmgbsa.LEAP.format(radii="mbondi2",
                                                 complex_pdb=os.path.basename(cclean),
                                                 prefix="c"))
        subprocess.run([os.path.join(home, "bin", "tleap"), "-f", leap], cwd=wd2,
                       capture_output=True, text=True, env=dict(os.environ, AMBERHOME=home))
        check(os.path.isfile(os.path.join(wd2, "c_com.parm7"))
              and os.path.getsize(os.path.join(wd2, "c_com.parm7")) > 0,
              "and tleap builds the capped complex")
    else:
        print("skip  tleap on the capped complex needs Amber")
else:
    print("skip  MM-GBSA's peptide protonation needs PDB2PQR")

# --- the page and the log ---------------------------------------------------------------------
log_rd = os.path.join(TMP, "logged")
os.makedirs(os.path.join(log_rd, "receptor"))
with open(os.path.join(log_rd, "receptor", "preparation.log"), "w") as fh:
    fh.write("  receptor came prepared: done; x\n  protonation: done; as the receptor was "
             "prepared, 743 hydrogens kept\n")
check(server.prep_warnings(log_rd, j) == [],
      "a receptor that came prepared is not flagged as one whose preparation failed")

# --- the prepare stage, with the ADCP 1.1 suite ----------------------------------------------
T = toolmod.discover(roots=ROOTS)
if T.prepare_receptor is None or not T.suite_program("agfr"):
    print("skip  the prepare stage needs the ADCP 1.1 suite (ADCPSUITE_HOME)")
else:
    # prepare_receptor types the Amber names as it types the standard ones.
    for name, path in (("prepared", PREPARED), ("capped", capped)):
        std = renamed(path, os.path.join(TMP, name + "_std.pdb"), STD)
        types = []
        for src in (path, std):
            out = os.path.join(TMP, os.path.basename(src) + ".%d.pdbqt" % len(types))
            receptor._run_prepare_receptor(T.prepare_receptor, src, out, add_h=False)
            types.append({(l[21], l[22:27].strip(), l[12:16].strip()): (l[77:79].strip(),
                                                                        l[70:76].strip())
                          for l in open(out) if l.startswith("ATOM")})
        check(types[0] and types[0] == types[1],
              "%s: prepare_receptor gives the Amber names the types and charges of the standard "
              "ones (%d atoms, %d differ)"
              % (name, len(types[0]), sum(1 for k in types[0] if types[0][k] != types[1].get(k))))
        if name == "capped":
            cap = {k: v[0] for k, v in types[0].items() if k[1] in ("45", "47")}
            check(cap.get(("A", "45", "N")) == "N" and cap.get(("A", "47", "O")) == "OA"
                  and cap.get(("A", "45", "H")) == "HD",
                  "and the caps are typed as the amide they are: %s" % cap)
    called = []
    real = pipeline.prep.run
    pipeline.prep.run = lambda *a, **kw: called.append(1)
    try:
        jj = jobmod.load(rd)
        trg = pipeline.stage_prepare(rd, jj, log=lambda *a: None)
        jj = jobmod.load(rd)
        info = jj.get("receptor_info") or {}
        prot = info.get("protonation") or {}
        check(not called, "this package's preparation is not run")
        check(prot.get("method") == "prepared" and prot.get("ph") == 6.5
              and prot.get("states") == {"HID": 2},
              "the run records that the receptor came prepared, its states and the pH given: %s"
              % prot)
        plog = open(os.path.join(rd, "receptor", "preparation.log")).read()
        check("receptor came prepared: done" in plog and "pH  : 6.5" in plog,
              "and so does its preparation log:\n%s" % plog)
        check(residues(info.get("prepared_pdb", "")) == BEFORE,
              "the prepared receptor the run holds is the file it was given, name for name")
        q = residues(info.get("pdbqt", ""))
        check(tally(q) == tally(BEFORE),
              "the PDBQT docked into keeps the states: %s against %s" % (tally(q), tally(BEFORE)))
        hd1 = sum(1 for l in open(info["pdbqt"]) if l.startswith("ATOM")
                  and l[17:20] == "HID" and l[12:16].strip() == "HD1")
        check(hd1 == 2, "with each HID's ring proton on ND1: %d" % hd1)
        check(sorted(tuple(x) for x in jj["site"]["hotspots"])
              == sorted(("A", int(s.split(":")[1])) for s in SITE.split(",")),
              "the site residues keep their numbers: %s" % jj["site"]["hotspots"])
        check(trg and os.path.isfile(trg), "and the maps are built: %s" % trg)
        # The record reaches the results' constants.
        check(pipeline._constants(jj, rd).get("protonation", {}).get("method") == "prepared",
              "the results' constants say the receptor came prepared")
        # The same run made from the receptor without hydrogens, by hand.
        jb = dict(jj, receptor=bare)
        rb = os.path.join(TMP, "bare_run")
        try:
            pipeline.stage_prepare(rb, jb, log=lambda *a: None)
            stopped = ""
        except SystemExit as exc:
            stopped = str(exc)
        check("carry no hydrogens" in stopped and not called,
              "the prepare stage refuses it too, and prepares nothing: %r" % stopped[:200])
    finally:
        pipeline.prep.run = real

shutil.rmtree(TMP, ignore_errors=True)
if fails:
    print("test_receptor_prepared FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  a receptor that came prepared is docked as it is, %d checks" % ran)
