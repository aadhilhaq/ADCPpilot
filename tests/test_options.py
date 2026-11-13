"""The docking and map options: the commands each one builds, what is
refused, and what a run records and leaves on disk.

By Aadhil Haq

Nothing here docks or builds maps. ADCP and agfr are replaced by small
scripts that write down how they were called, so what is checked is the
command line each setting produces, the refusals made before anything is
submitted, the packing of trajectories and a resume after it, the box from
a ligand, and the round trip of every new setting through job.json, the
command line and the batch scripts.
"""
import contextlib
import io
import json
import os
import stat
import sys
import tarfile
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_pilot import (adcp11, cli, job as jobmod, peptide, pipeline,  # noqa: E402
                       sizing, slurm, structure, target, trajectories)

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


def run_cli(*argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            rc = cli.main(list(argv))
        except SystemExit as e:
            rc = e.code
        except Exception as e:                                  # noqa: BLE001
            rc = "CRASH %s: %s" % (type(e).__name__, e)
    return rc, out.getvalue()


TMP = tempfile.mkdtemp(prefix="adcp_options_")
REC = os.environ.get("ADCP_TEST_RECEPTOR") or os.path.join(ROOT, "examples", "receptor.pdb")
if not os.path.isfile(REC):
    # examples/receptor.pdb is the user's own and ignored by git, so a fresh clone has none.
    # The package's own example ships MDM2 (chain A) with p53's peptide (chain B), a site in A
    # and a ligand in B, which is what these checks ask of a receptor.
    REC = os.path.join(ROOT, "adcp_pilot", "example_run", "receptor", "1YCR.pdb")
# Where the ADCP 1.1 suite is, for the checks that need it: ADCP_TEST_TOOLS, or else
# ADCPSUITE_HOME. Without either they are skipped.
ROOTS = [r for r in (os.environ.get("ADCP_TEST_TOOLS", ""),
                     os.environ.get("ADCPSUITE_HOME", ""))
         if r and os.path.isdir(r)]


# --- % helical: what -p is given, and how many replicas start helical ------
#
# runADCP.py rewrites the sequence's case replica by replica only when
# 0 < p < 100; at 0 and 100 the case of the sequence decides.
seq, extra = adcp11.start_arguments("ETFSN<BHD>LW", 0)
check(seq == "etfsn<BHD>lw" and extra == [],
      "0%% helical should be the old behavior, lower case and no -p; got "
      "%r %r" % (seq, extra))
seq, extra = adcp11.start_arguments("ETFSN<BHD>LW", 100)
check(seq == "ETFSN<BHD>LW" and extra == [],
      "100%% helical should be upper case and no -p (ADCP ignores -p 100); "
      "got %r %r" % (seq, extra))
seq, extra = adcp11.start_arguments("ETFSN<BHD>LW", 25)
check(extra == ["-p", "25"] and "<BHD>" in seq,
      "25%% helical should pass -p 25 and keep the library name as written; "
      "got %r %r" % (seq, extra))
seq, extra = adcp11.start_arguments("etfsn<bhd>lw".upper(), 0)
check("<BHD>" in seq, "the library name's case was changed: %r" % seq)
check(adcp11.helical_replicas(50, 25) == 12
      and adcp11.helical_replicas(100, 25) == 25
      and adcp11.helical_replicas(50, 0) == 0
      and adcp11.helical_replicas(50, 100) == 50,
      "the count of helical starts does not follow runADCP.py's "
      "jobNum <= nbRuns * p / 100")

# --- the ADCP 1.1 command line for each option ------------------------------
args = adcp11.dock_arguments("ETFSDLWKLLPEN", "receptor.trg", 50, 2500000,
                             60, cpus=48, helical=25, cyclic="backbone",
                             keep=True)
for flag in (["-m", "60"], ["-p", "25"], ["-cyc"], ["-k"], ["-c", "48"],
             ["-N", "50"], ["-n", "2500000"]):
    i = args.index(flag[0]) if flag[0] in args else -1
    check(i >= 0 and args[i:i + len(flag)] == flag,
          "ADCP 1.1 was not given %s: %s" % (" ".join(flag), args))
args = adcp11.dock_arguments("CETFSDLWKLLPC", "receptor.trg", 50, 2500000,
                             10, cyclic="disulfide")
check("-cys" in args and "-cyc" not in args and "-k" not in args
      and "-p" not in args,
      "a disulfide ring, 0%% helical and no trajectories built %s" % args)
args = adcp11.dock_arguments("ETFSDLWKLLPEN", "receptor.trg", 50, 2500000, 10)
check(not {"-cys", "-cyc", "-k", "-p"} & set(args)
      and args[args.index("-s") + 1] == "etfsdlwkllpen",
      "a run with every new option at its old setting should build the old "
      "command, lower case and nothing added: %s" % args)


def fake(path, body):
    with open(path, "w") as fh:
        fh.write("#!/bin/sh\n" + body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IEXEC)
    return path


bin_dir = os.path.join(TMP, "bin")
os.makedirs(bin_dir)
fake(os.path.join(bin_dir, "adcp"), "exit 0\n")
trg = os.path.join(TMP, "receptor.trg")
open(trg, "w").write("not a real target")

# --- packing the trajectories, and a resume after it ------------------------
# ADCP 1.1 keeps out/ with -k: run_<n> files, out.pdb and the unpacked target.
wd11 = os.path.join(TMP, "dock11", "ETFSDLWKLLPEN")
os.makedirs(os.path.join(wd11, "out", "receptor"))
for i in (1, 2):
    for x in ("pdb", "out", "err"):
        open(os.path.join(wd11, "out", "run_%d.%s" % (i, x)), "w").write("x")
open(os.path.join(wd11, "out", "out.pdb"), "w").write("all of them")
open(os.path.join(wd11, "out", "receptor", "rigidReceptor.C.map"),
     "w").write("map")
open(os.path.join(wd11, "out_ranked_1.pdb"), "w").write("ATOM\nEND\n")
open(os.path.join(wd11, ".adcp_finished"), "w").write("done\n")
# A pack that was stopped part way left a partial archive behind.
open(trajectories.archive_path(wd11, "ETFSDLWKLLPEN") + ".part", "w").write("x")
arc = trajectories.finish(wd11, "ETFSDLWKLLPEN", True, log=lambda *a: None)
names = sorted(tarfile.open(arc).getnames()) if arc else []
check(names == sorted("out/run_%d.%s" % (i, x) for i in (1, 2)
                      for x in ("pdb", "out", "err")),
      "ADCP 1.1's archive holds %s; it should leave out out.pdb and the "
      "unpacked maps" % names)
check(not os.path.isdir(os.path.join(wd11, "out"))
      and os.path.isfile(os.path.join(wd11, trajectories.PACKED)),
      "the working folder was not removed after it was packed")
# A resume of this peptide must not dock it again: the ranked pose and the
# marker are still there, and the archive is what is kept.
called = []
_real_run = adcp11._run
adcp11._run = lambda cmd, wd, exe, append=False: called.append(cmd)
try:
    adcp11.dock("ETFSDLWKLLPEN", trg, os.path.join(TMP, "dock11"),
                os.path.join(TMP, "bin", "adcp"), key="ETFSDLWKLLPEN",
                keep=True, log=lambda *a: None)
finally:
    adcp11._run = _real_run
check(not called, "a peptide already docked and packed was docked again")
check(os.path.isfile(arc), "a resume removed the archive")
check(trajectories.finish(wd11, "ETFSDLWKLLPEN", True, log=lambda *a: None) == arc,
      "packing again does not find the archive already made")
# Off: a working folder left behind is removed; absent: left as it was.
wd = os.path.join(TMP, "loose")
os.makedirs(os.path.join(wd, "out"))
for n in (os.path.join("out", "run_1.pdb"), "out_ranked_1.pdb"):
    open(os.path.join(wd, n), "w").write("x")
trajectories.finish(wd, "K", None, log=lambda *a: None)
check(sorted(os.listdir(wd)) == ["out", "out_ranked_1.pdb"],
      "a job file without the setting had its trajectories touched")
trajectories.finish(wd, "K", False, log=lambda *a: None)
check(sorted(os.listdir(wd)) == ["out_ranked_1.pdb"],
      "with the setting off the working folder was not removed: %s"
      % os.listdir(wd))

# --- cyclic: refused before anything is submitted ---------------------------
check(peptide.cyclic_problem("ETFSDLWKLLPEN", "disulfide") is not None,
      "a disulfide ring with no cysteine was accepted")
check(peptide.cyclic_problem("ETFSDLWKLLPEC", "disulfide") is not None,
      "a disulfide ring with one cysteine was accepted")
check(peptide.cyclic_problem("CETFSDLWKLLPC", "disulfide") is None,
      "a disulfide ring with two cysteines was refused")
check(peptide.cyclic_problem("ETFSDLWKLLPEN", "backbone") is None,
      "a head to tail ring was refused")
check(peptide.cysteines("C<HCS>TFSDLWKLLPC") == 1,
      "a non-standard residue on the cysteine model was counted as a "
      "cysteine; its sulfur is not ADCP's SG")
msg = peptide.cyclic_problem("ETFSDLWKLLPEC", "disulfide").message
check("two cysteines" in msg and "one" in msg,
      "the refusal does not say why: %r" % msg)
probs = jobmod.validate(jobmod.default(
    receptor=REC, hotspots="A:54", cyclic="disulfide",
    peptides=[{"name": "p53", "sequence": "ETFSDLWKLLPEN"}]))
check(any("two cysteines" in p for p in probs),
      "validate does not refuse a disulfide ring without cysteines: %s"
      % probs)
rd = os.path.join(TMP, "ring")
rc, out = run_cli("new", rd, "-r", REC, "--chain", "A", "-s", "A:54",
                  "-p", "p53: ETFSDLWKLLPEN", "--cyclic", "disulfide")
# A SystemExit carrying a sentence is what the shell prints and exits 1 on.
said = (rc if isinstance(rc, str) else "") + out
check(rc not in (0, None) and "two cysteines" in said
      and "p53" in said and not os.path.isfile(os.path.join(rd, "job.json")),
      "the command line wrote a job for a disulfide ring that cannot close, "
      "or did not say why (rc %r): %s" % (rc, out[-300:]))

# --- the box from a ligand ---------------------------------------------------
atoms = structure.read(REC)
ligs = structure.ligands(atoms)
check([l["id"] for l in ligs] == ["B"] and ligs[0]["kind"] == "peptide",
      "1YCR's p53 peptide is not offered as what the box can be built "
      "around: %s" % ligs)
lig = structure.ligand_atoms(atoms, "B")
c, s, note = structure.ligand_box(lig, padding=4.0, minimum=0, maximum=100)
xs = [a.x for a in lig]
check(abs(c[0] - (min(xs) + max(xs)) / 2) < 1e-6
      and abs(s[0] - (max(xs) - min(xs) + 8.0)) < 1e-6,
      "the ligand box does not follow AGFR's rule for -b ligand")
c2, s2, _n = structure.ligand_box(lig, 8.0, 20.0, 30.0)
check(all(20.0 <= v <= 30.0 for v in s2),
      "the ligand box is not held between the minimum and maximum: %s"
      % (s2,))
contacts = structure.contacts(atoms, lig, chains_keep=["A"])
check(16 <= len(contacts) <= 18 and ("A", 58) in contacts,
      "the residues p53 touches in 1YCR are %d (%s); 17 expected"
      % (len(contacts), contacts[:5]))
check(structure.ligand_problems(atoms, "B", ["A", "B"], False),
      "a peptide ticked as one of the receptor's chains was accepted as "
      "the ligand that sets the box")
check(not structure.ligand_problems(atoms, "B", [], False),
      "a peptide with no chains ticked was refused")
check(structure.ligand_problems(atoms, "A:XYZ:9", [], False),
      "a ligand that is not in the structure was accepted")
check(pipeline.receptor_chains({"receptor": REC, "chains": [],
                                "site_mode": "ligand",
                                "site_ligand": "B"}) == ["A"],
      "with no chains ticked, the receptor kept the peptide that sets the box")
j = jobmod.default(receptor=REC, site_mode="ligand", site_ligand="B",
                   peptides=[{"name": "p53", "sequence": "ETFSDLWKLLPEN"}])
check(not jobmod.validate(j) or all("ADCP 1.1" in p for p in jobmod.validate(j)),
      "a box from chain B of 1YCR was refused: %s" % jobmod.validate(j))
check(not [p for p in jobmod.validate(dict(j, site_mode="pocket"))
           if "site" in p.lower() or "residue" in p.lower()],
      "a box from an AutoSite pocket asked for residues")
check([p for p in jobmod.validate(dict(j, site_mode="residues",
                                       hotspots=""))
       if "residues" in p],
      "a box from residues with no residues was accepted")

# --- the AGFR command for each map setting ----------------------------------
# Both settings are given explicitly, whatever the program's defaults.
check(target.Maps(None, "1.1", False).flags() == ["-asv", "1.1"]
      and target.Maps(None, "1.0", True).flags() == ["-asv", "1.0", "-g"],
      "AGFR's flags are wrong for the gradient or AutoSite: %s, %s"
      % (target.Maps(None, "1.1", False).flags(),
         target.Maps(None, "1.0", True).flags()))
check(target.AGFR_CHOICES == ("1.1", "auto") and "1.0" in target.AUTOSITE_CHOICES,
      "the maps are not AGFR 1.1's alone, or AutoSite 1.0 went with the older "
      "AGFR: %s %s" % (target.AGFR_CHOICES, target.AUTOSITE_CHOICES))
T = types.SimpleNamespace()
m = target.Maps(None, "1.1", False)
m.runner = lambda tools=None: ["agfr"]
cmd = target.command(m, T, "/r/receptor.pdbqt", "receptor",
                     (1.0, 2.0, 3.0), (30.0, 30.0, 30.0))
check(cmd == ["agfr", "-r", "/r/receptor.pdbqt", "-o", "receptor", "-b",
              "user", "1.000", "2.000", "3.000", "30.0", "30.0", "30.0",
              "-asv", "1.1"],
      "the map command changed: %s" % cmd)
cmd = target.command(m, T, "/r/r.pdbqt", "pocket", pocket="best", padding=8)
check(cmd[cmd.index("-b"):] == ["-b", "fill", "-p", "best", "-P", "8",
                                "-asv", "1.1"],
      "the pocket command is %s" % cmd)
# The pocket table agfr prints, read back.
table = target.pocket_table(
    "Scanning at: -0.36 -0.79 -0.6\n"
    "    1    -12.00    10    1.00    -0.10     0.500     1.00\n"
    "    found 2 pocket(s)\n\n"
    "    pocket|  energy | # of |Rad. of | energy |   bns    | score  \n"
    "    number|         |points|gyration|per vol.|buriedness|v*b^2/rg\n"
    "    ------+---------+------+--------+--------+----------+---------\n"
    "        1   -175.24   513    6.44     -0.34      0.78      48.85\n"
    "        2    -99.24   309    4.57     -0.32      0.75      37.65\n"
    "    using best score cluster with 513 points\n")
check([r["points"] for r in table] == [513, 309] and table[0]["score"] == 48.85,
      "the pocket table was read as %s" % table)
exe, script = adcp11.find_agfr(ROOTS)
if exe:
    m11 = target.maps_for(jobmod.default(tools_roots=ROOTS))
    check(m11.agfr == "1.1" and m11.autosite == "1.1" and not m11.gradient,
          "a new job does not build its maps with AGFR 1.1, AutoSite 1.1 "
          "and no gradient: %s" % m11.describe())
    cmd = target.command(m11, T, "/r/r.pdbqt", "receptor", (0, 0, 0),
                         (20, 20, 20))
    check(cmd[0].endswith("/bin/python") and cmd[1].endswith("runAGFR.py")
          and "-g" not in cmd and cmd[-2:] == ["-asv", "1.1"],
          "AGFR is not run as the suite's Python on runAGFR.py: %s" % cmd)
    # A job file from before the maps could be chosen, or one that named the
    # older AGFR: built again with what comes closest to its own maps.
    legacy = target.maps_for({"tools_roots": ROOTS})
    check((legacy.autosite, legacy.gradient) == ("1.0", True)
          and legacy.flags() == ["-asv", "1.0", "-g"],
          "a job file without the map settings is not built with AutoSite 1.0 "
          "and the gradient on: %s" % legacy.describe())
    named = target.maps_for({"tools_roots": ROOTS, "adcp_engine": "1.0",
                             "agfr": "1.0", "autosite": "1.1",
                             "receptor_gradient": False})
    check(named.flags() == ["-asv", "1.1"] and named.agfr == "1.1",
          "a job file that named the older AGFR and recorded its settings is "
          "not built with them by AGFR 1.1: %s" % named.describe())
else:
    print("  . the ADCP 1.1 suite is not here (ADCP_TEST_TOOLS); its AGFR's "
          "commands were not checked")

# --- the job record, the command line and the batch scripts -----------------
d = jobmod.default()
for k, v in (("site_mode", "residues"), ("site_ligand", ""),
             ("site_pocket", "best"), ("agfr", "1.1"), ("autosite", "1.1"),
             ("receptor_gradient", False), ("cyclic", "none"),
             ("helical", 25), ("keep_trajectories", False), ("poses", 10)):
    check(d.get(k) == v, "a new job's %s is %r, not %r" % (k, d.get(k), v))
old = {k: v for k, v in d.items()
       if k not in ("helical", "cyclic", "keep_trajectories")}
o = pipeline.docking_options(old)
# "seed": None is what every run did before the setting was honoured: ADCP
# draws its own seed. The old job file therefore still docks as it did.
check(o == {"helical": 0, "cyclic": "none", "seed": None, "keep": None},
      "a job file from before the options does not dock as it did: %s" % o)
rd = os.path.join(TMP, "cli_opts")
for opt in ("--engine", "--agfr"):
    rc, out = run_cli("new", os.path.join(TMP, "cli_old"), "-r", REC,
                      "-s", "A:54", opt, "1.0", "-p", "p53: ETFSDLWKLLPEN")
    check(rc == 2 and "invalid choice" in out
          and not os.path.isdir(os.path.join(TMP, "cli_old")),
          "%s 1.0 was accepted (rc %r)" % (opt, rc))
rc, out = run_cli("new", rd, "-r", REC, "--site-mode", "ligand",
                  "--site-ligand", "B", "--agfr", "auto", "--autosite", "1.0",
                  "--receptor-gradient", "--cyclic", "backbone",
                  "--helical", "40", "--keep-trajectories", "--poses", "50",
                  "--engine", "auto", "-p", "p53: ETFSDLWKLLPEN")
jf = os.path.join(rd, "job.json")
j = json.load(open(jf)) if os.path.isfile(jf) else {}
for k, v in (("site_mode", "ligand"), ("site_ligand", "B"), ("agfr", "auto"),
             ("autosite", "1.0"), ("receptor_gradient", True),
             ("cyclic", "backbone"), ("helical", 40),
             ("keep_trajectories", True), ("poses", 50)):
    check(j.get(k) == v, "--%s did not reach job.json: %r"
          % (k.replace("_", "-"), j.get(k)))
check(not any("poses" in p for p in jobmod.validate(j)) if j else False,
      "50 poses were refused")
check(any("at most 100" in p for p in jobmod.validate(dict(j, poses=101))),
      "101 poses were accepted")
check(any("adcp_engine must be" in p
          for p in jobmod.validate(dict(j, adcp_engine="1.0"))),
      "a new job naming another ADCP was accepted")
j2 = jobmod.load(jf) if j else {}
check(j2 == j, "job.json did not read back as it was written")
notes = "\n".join(pipeline.batch_notes(j, "prepare")
                  + pipeline.batch_notes(j, "dock")) if j else ""
for want in ("around the ligand B", "AGFR 1.1", "AutoSite 1.0", "ADCP 1.1",
             "receptor gradient on", "40% of replicas start as a helix",
             "adcp -p 40", "backbone cyclic", "50 pose(s)",
             "trajectories kept"):
    check(want in notes, "the batch scripts do not say %r:\n%s"
          % (want, notes))
# A run folder typed as a relative path reaches the job scripts absolute: they
# change into the run folder before they read job.json by that path.
if j:
    seen = []
    _submit = pipeline.submit
    pipeline.submit = lambda rundir, *a, **k: seen.append(rundir)
    here = os.getcwd()
    os.chdir(TMP)
    try:
        run_cli("run", "cli_opts")
    finally:
        os.chdir(here)
        pipeline.submit = _submit
    check([os.path.realpath(x) for x in seen]
          == [os.path.realpath(os.path.join(TMP, "cli_opts"))]
          and all(os.path.isabs(x) for x in seen),
          "adcp-pilot run passed the run folder on as %s; the job scripts "
          "need it absolute" % seen)
script = os.path.join(TMP, "02_dock.sbatch")
slurm.write_script(script, name="t", cpus=2, time="01:00:00", workdir=TMP,
                   logdir=TMP, python="python", jobfile=jf, mem="1G",
                   indices="0", notes=["line one\nrm -rf /", "two"])
text = open(script).read()
check("# line one rm -rf /" in text and "\nrm -rf" not in text,
      "a setting's text could end a comment line in the batch script")

# --- what each new option costs ---------------------------------------------
base = dict(agfr="1.1", autosite="1.0", receptor_gradient=False,
            site_mode="residues", cyclic="none", helical=0, poses=10,
            keep_trajectories=False)
plain = sizing.plan([13], 50, 2500000, 48, options=base)
for k, v in (("receptor_gradient", True), ("autosite", "1.1"),
             ("site_mode", "pocket")):
    p = sizing.plan([13], 50, 2500000, 48, atoms=1500, box=(30, 30, 30),
                    options=dict(base, **{k: v}))
    q = sizing.plan([13], 50, 2500000, 48, atoms=1500, box=(30, 30, 30),
                    options=base)
    check(p["prep_minutes"] > q["prep_minutes"],
          "%s=%r costs the prepare stage nothing" % (k, v))
for k, v in (("cyclic", "backbone"), ("cyclic", "disulfide"),
             ("keep_trajectories", True), ("poses", 100)):
    p = sizing.plan([13], 50, 2500000, 48, options=dict(base, **{k: v}))
    check(p["task_minutes"] > plain["task_minutes"],
          "%s=%r costs the docking task nothing" % (k, v))
check(sizing._book(60) >= 2 * 60 + 30 and sizing._book(1) >= 120,
      "a booking is not twice the estimate and a buffer, at least two hours")

if fails:
    print("test_options FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the new options build their commands, refuse what cannot run, "
      "pack and resume, and round trip, %d checks" % ran)
