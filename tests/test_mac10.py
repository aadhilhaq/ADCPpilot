"""ADCP 1.0 on a Mac, and nowhere else.

By Aadhil Haq

On macOS, where the Scripps Research Institute builds no ADCP 1.1, `adcp-dock
local` installs ADFRsuite 1.0 and docks with it (adcp_dock/mac10.py). No Mac
is needed here: the system mac10 sees is set by hand, and ADFRsuite is a
folder of stand-in programs that write what ADCP 1.0 writes. What is checked:

  * the gate: on any system but macOS nothing reaches ADCP 1.0, even with
    ADFRsuite 1.0 installed and named and its programs on PATH;
  * on a Mac: the tools, reduce and msms come from ADFRsuite 1.0; maps from
    its agfr with AutoSite 1.0 and the gradient as the run asks; a new run
    says "ADCP 1.0 (macOS)"; what ADCP 1.0 cannot do is refused;
  * the docking: the command ADCP 1.0 is given, and straight after it the
    folder in the 1.1 layout (modes.json in the form the 1.1 reader takes,
    the loose replica files packed or removed), resumed without docking
    again, and a failure said;
  * a whole dock task and its report, which say what the engine was and that
    the scores are not comparable with ADCP 1.1's;
  * the install: the macOS tarball, checked against its pinned SHA-256, and
    refused on any other system; `adcp-dock local` choosing it on a Mac
    only, and `adcp-dock check` naming it there and only there.

The real programs on a real Mac are tests/mac10_live.py's, which the manual
mac10 workflow runs on GitHub's macOS runner.

    python3 tests/test_mac10.py
"""
import contextlib
import hashlib
import io
import json
import os
import shutil
import stat
import sys
import tarfile
import tempfile
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import (adcp11, adcp11_install, cli, job as jobmod,  # noqa: E402
                       mac10, pipeline, server, target, tools)

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)
        print("  FAIL " + why)


TMP = tempfile.mkdtemp(prefix="adcp_mac10_")
P53 = "ETFSDLWKLLPEN"
TABLE = """Performing search (8 ADCP runs with 1000 steps each) ...
mode |  affinity  | clust. | ref. | clust. | rmsd | best | energy | best |
     | (kcal/mol) | rmsd   | rmsd |  size  | avg. | rmsd |  avg.  | run  |
-----+------------+--------+------+--------+------+------+--------+------+
   1        -21.1     0.0   999.0      21   999.0 999.0   -19.1    075
   2        -19.2    21.3   999.0      10   999.0 999.0    -9.5    582
   3        -18.7    12.2   999.0       1      NA      NA    1142
"""
FAKE_ADCP = r'''#!/bin/sh
printf "%s\n" "$@" > argv.txt
echo run >> runs.txt
if [ -f "$FAIL_MARK" ]; then
  echo "Building protein from sequence: nwr (3 amino acids, 1 chains)."
  for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do echo "$i ENDED WITH ERROR"; done
  exit 0
fi
cat "$TABLE_FILE"
for k in 1 2 3; do
  printf "ATOM      1  N   GLU A   1       %d.000   0.000   0.000  1.00  0.00           N\nEND\n" $k > out_ranked_$k.pdb
done
for i in 1 2 3; do echo "traj $i" > out_$i.pdb; echo "log $i" > out_$i.out; done
echo all > out.pdb; echo 1 > constrains; echo r > ramaprob.data; echo t > tmp.pdb
'''
table_file = os.path.join(TMP, "table.txt")
open(table_file, "w").write(TABLE)
os.environ["TABLE_FILE"] = table_file
os.environ["FAIL_MARK"] = os.path.join(TMP, "fail_now")


def make_suite(path):
    """A stand-in ADFRsuite 1.0: its programs, each writing down how it was
    called, and adcp writing what ADCP 1.0 writes."""
    os.makedirs(os.path.join(path, "bin"))
    for name in ("adcp", "agfr", "prepare_receptor", "reduce", "msms"):
        p = os.path.join(path, "bin", name)
        open(p, "w").write(FAKE_ADCP if name == "adcp" else
                           '#!/bin/sh\necho "%s $*" >> "%s"\n'
                           % (name, os.path.join(TMP, "called.txt")))
        os.chmod(p, os.stat(p).st_mode | stat.S_IEXEC)
    return path


def check_output():
    """What `adcp-dock check` prints."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        try:
            cli.main(["check"])
        except SystemExit:
            pass
    return out.getvalue()


ADFR = make_suite(os.path.join(TMP, "ADFRsuite"))
os.environ[mac10.ENV] = ADFR
_real_system = mac10._system


def as_system(name):
    mac10._system = lambda: name


# --- on Linux, nothing reaches ADCP 1.0 ------------------------------------------
# ADFRsuite 1.0 installed, named in the variable `adcp-dock local` sets, and its
# programs first on PATH: none of it is used, and every function that would
# run it fails the check if it is reached.
as_system("Linux")
reached = []
_saved = {n: getattr(mac10, n) for n in ("dock", "install", "maps", "program",
                                         "to_layout", "mode_table")}
for n in _saved:
    setattr(mac10, n, (lambda name: lambda *a, **k: reached.append(name))(n))
_path = os.environ["PATH"]
try:
    check(not mac10.is_mac() and mac10.home() is None and not mac10.engaged(),
          "the Mac route is engaged on Linux")
    os.environ["PATH"] = os.path.join(ADFR, "bin") + os.pathsep + _path
    T = tools.discover()
    check(not any((T.get(n) or "").startswith(ADFR)
                  for n in ("adcp", "agfr", "prepare_receptor")),
          "on Linux a docking program was taken from ADFRsuite 1.0: %s"
          % [T.get(n) for n in ("adcp", "agfr", "prepare_receptor")])
    check(not (tools.suite_program("reduce") or "").startswith(ADFR),
          "on Linux reduce was taken from ADFRsuite 1.0")
    # An msms on the PATH is the last of MSMS's places on any system, so
    # ADFRsuite's is off it here: what is checked is that its folder,
    # named for the Mac, is not one of them.
    os.environ["PATH"] = _path
    check(ADFR not in str(tools.msms_command(T)),
          "on Linux MSMS was taken from ADFRsuite 1.0: %s" % tools.msms_command(T))
    j11 = jobmod.default(peptides=[{"name": "p53", "sequence": P53}])
    check(jobmod.resolve_engine(j11)[0] in ("1.1", "missing"),
          "a run on Linux does not dock with ADCP 1.1")
    jmac = dict(j11, adcp_engine=mac10.ENGINE)
    check(jobmod.resolve_engine(jmac) == ("mac-only", None)
          and any("made on a Mac" in p for p in jobmod.validate(jmac)),
          "a run made on a Mac is docked, or not refused, on Linux")
    check(jobmod.for_this_machine(dict(j11)).get("adcp_engine") == "1.1",
          "a new run on Linux is not an ADCP 1.1 run")
    try:
        m = target.maps_for({"tools_roots": []})
        check(getattr(m, "agfr", None) == "1.1",
              "the maps on Linux are not AGFR 1.1's: %r" % m)
    except target.TargetError:
        pass                    # no ADCP 1.1 suite here: no maps at all
    rd = os.path.join(TMP, "linux_mac_run")
    os.makedirs(os.path.join(rd, "target"))
    open(os.path.join(rd, "target", "receptor.trg"), "w").write("t")
    jobmod.save(rd, dict(jmac, target=os.path.join(rd, "target", "receptor.trg")))
    _check = pipeline.target.check
    pipeline.target.check = lambda trg: None
    try:
        pipeline.stage_dock(rd, 0, log=lambda *a: None)
        stopped = ""
    except SystemExit as exc:
        stopped = str(exc)
    finally:
        pipeline.target.check = _check
    check("made on a Mac" in stopped,
          "a run made on a Mac was docked on Linux (%r)" % stopped)
    said = check_output()
    check("ADFRsuite" not in said and "ADCP 1.0" not in said,
          "adcp-dock check on Linux names ADFRsuite 1.0")
    check(not reached, "on Linux these were reached: %s" % reached)
    try:
        _saved["install"](os.path.join(TMP, "never"))
        refused = ""
    except SystemExit as exc:
        refused = str(exc)
    check("only on a Mac" in refused and not os.path.exists(os.path.join(TMP, "never")),
          "ADFRsuite 1.0 would install on Linux (%r)" % refused)
finally:
    for n, f in _saved.items():
        setattr(mac10, n, f)
    os.environ["PATH"] = _path

# --- on a Mac: the tools, the run, the maps ------------------------------------------
as_system("Darwin")
check(mac10.engaged() and mac10.home() == os.path.abspath(ADFR),
      "a Mac with ADFRsuite 1.0 named does not dock with it")
T = tools.discover()
for name in ("adcp", "agfr", "prepare_receptor"):
    check(T.get(name) == os.path.join(ADFR, "bin", name),
          "on a Mac %s is not ADFRsuite 1.0's: %r" % (name, T.get(name)))
check(tools.suite_program("reduce") == os.path.join(ADFR, "bin", "reduce"),
      "on a Mac reduce is not ADFRsuite 1.0's")
check(tools.msms_command(T) == os.path.join(ADFR, "bin", "msms"),
      "on a Mac MSMS is not ADFRsuite 1.0's msms program: %r"
      % tools.msms_command(T))
j = jobmod.for_this_machine(jobmod.default(
    replicas=8, steps=1000, poses=3, cpus=2, helical=25, cyclic="backbone",
    peptides=[{"name": "p53", "sequence": P53}]))
check(j["adcp_engine"] == mac10.ENGINE and j["engine"] == "ADCP 1.0 (macOS)",
      "a new run on a Mac does not say it docks with ADCP 1.0: %s"
      % {k: j.get(k) for k in ("adcp_engine", "engine")})
check(jobmod.resolve_engine(j) == (mac10.ENGINE, os.path.abspath(ADFR)),
      "a Mac run does not resolve to ADFRsuite 1.0")
check(any("ADCP 1.1" in p and "macOS" in p
          for p in jobmod.validate(dict(j, adcp_engine="1.1"))),
      "a run asking for ADCP 1.1 is not refused on a Mac")
probs = jobmod.validate(dict(j, peptides=[{"name": "nst",
                                          "sequence": "FSN<BHD>LW&KLL"}]))
check(any("macOS" in p and "Non-standard" in p for p in probs),
      "a non-standard peptide is not refused on a Mac: %s" % probs)
check(any("OpenMM" in p for p in jobmod.validate(dict(j, omm_minimize=2))),
      "OpenMM minimization is not refused on a Mac")
check(any("at most 10" in p for p in jobmod.validate(dict(j, poses=11)))
      and not any("at most" in p for p in jobmod.validate(dict(j, poses=10))),
      "a Mac run's poses are not held to ADCP 1.0's ten")
m = target.maps_for(j)
check(m.autosite == "1.0" and m.flags() == ["-ng"]
      and target.maps_for(dict(j, receptor_gradient=True)).flags() == []
      and "AGFR 1.0 (macOS)" in m.describe(),
      "the maps on a Mac are not AGFR 1.0's, AutoSite 1.0, with the gradient "
      "as the run asks: %s %s" % (m.describe(), m.flags()))
cmd = target.command(m, T, "/r/receptor.pdbqt", "receptor", (1, 2, 3),
                     (30, 30, 30))
check(cmd[0] == os.path.join(ADFR, "bin", "agfr") and cmd[-1] == "-ng"
      and "-asv" not in cmd,
      "the agfr command on a Mac is %s" % cmd)
said = check_output()
check("ADFRsuite 1.0" in said and mac10.NOTE in said and ADFR in said
      and "ADCP 1.1 suite" not in said,
      "adcp-dock check on a Mac does not say it docks with ADFRsuite 1.0")
check("ADCP 1.0 (macOS)" in "\n".join(pipeline.batch_notes(j, "dock")),
      "the dock script's notes do not name ADCP 1.0 (macOS)")

# --- the docking, and the 1.1 layout straight after --------------------------------
trg = os.path.join(TMP, "receptor.trg")
open(trg, "w").write("a target")
logs = []
modes, wd = mac10.dock(P53, trg, os.path.join(TMP, "dock"), key=P53,
                       replicas=8, steps=1000, cpus=2, helical=25,
                       cyclic="backbone", keep=False, log=logs.append)
argv = open(os.path.join(wd, "argv.txt")).read().split()
check(argv[argv.index("-t") + 1] == "receptor.trg"
      and argv[argv.index("-s") + 1] == P53 and argv[argv.index("-p") + 1] == "25"
      and "-cyc" in argv and argv[argv.index("-c") + 1] == "2"
      and argv[argv.index("-N") + 1] == "8",
      "ADCP 1.0 was not given the target, the sequence, -p 25, -cyc and the "
      "cores: %s" % argv)
want = [{"rank": 1, "affinity": -21.1, "cluster": 21},
        {"rank": 2, "affinity": -19.2, "cluster": 10},
        {"rank": 3, "affinity": -18.7, "cluster": 1}]
js = os.path.join(wd, "modes.json")
check(os.path.isfile(js) and json.load(open(js)) == want and modes == want,
      "modes.json is not ADCP 1.0's mode table: %s"
      % (open(js).read() if os.path.isfile(js) else "missing"))
check(open(js).read() == json.dumps(want, indent=1),
      "modes.json is not written as the 1.1 docking writes it")
check(adcp11.read_modes(wd) == want, "the 1.1 reader does not read it")
left = sorted(os.listdir(wd))
check(left == sorted([".adcp_finished", "adcp.log", "argv.txt", "modes.json",
                      "out_ranked_1.pdb", "out_ranked_2.pdb", "out_ranked_3.pdb",
                      "receptor.trg", "runs.txt"]),
      "the folder after the docking is not the 1.1 layout: %s" % left)
modes2, _wd = mac10.dock(P53, trg, os.path.join(TMP, "dock"), key=P53,
                         replicas=8, steps=1000, log=lambda *a: None)
check(open(os.path.join(wd, "runs.txt")).read().count("run") == 1
      and modes2 == want, "a finished peptide was docked again on a resume")
# Kept: the replicas in one archive, read back, under out/ as 1.1 keeps them.
modes, wdk = mac10.dock(P53, trg, os.path.join(TMP, "dock_keep"), key=P53,
                        replicas=8, steps=1000, keep=True, log=lambda *a: None)
arc = os.path.join(wdk, "%s_trajectories.tar.gz" % P53)
names = sorted(tarfile.open(arc).getnames()) if os.path.isfile(arc) else []
check(names == sorted("out/out_%d.%s" % (i, x) for i in (1, 2, 3)
                      for x in ("pdb", "out")),
      "the kept replicas are not in one archive: %s" % names)
check(not [n for n in os.listdir(wdk) if n.startswith("out_")
           and not n.startswith("out_ranked")] and "out.pdb" not in os.listdir(wdk)
      and "ramaprob.data" not in os.listdir(wdk),
      "loose replica files were left: %s" % sorted(os.listdir(wdk)))
# A docking ADCP could not do: said, and not marked finished.
open(os.environ["FAIL_MARK"], "w").write("x")
logs = []
modes, wdf = mac10.dock("NWR", trg, os.path.join(TMP, "dock_fail"), key="NWR",
                        replicas=20, steps=1000, log=logs.append)
os.remove(os.environ["FAIL_MARK"])
check(modes == [] and not os.path.isfile(os.path.join(wdf, ".adcp_finished"))
      and any("could not dock this at all" in l and "every replica failed" in l
              for l in logs),
      "a failed docking is not said, or is marked finished: %s" % logs)

# --- a whole dock task and its report --------------------------------------------
rd = os.path.join(TMP, "macrun")
os.makedirs(os.path.join(rd, "target"))
trg2 = os.path.join(rd, "target", "receptor.trg")
open(trg2, "w").write("a target")
jobmod.save(rd, dict(j, target=trg2, mmgbsa=False))
_check = pipeline.target.check
pipeline.target.check = lambda trg: None
try:
    lines = []
    pipeline.stage_dock(rd, 0, log=lines.append)
finally:
    pipeline.target.check = _check
res = json.load(open(os.path.join(rd, "results", P53 + ".json")))
check(res.get("engine") == "ADCP 1.0 (macOS)"
      and [m["score"] for m in res["modes"]] == [-21.1, -19.2, -18.7]
      and res.get("cluster_top") == 21,
      "the result does not name ADCP 1.0 (macOS) or lost its modes: %s"
      % {k: res.get(k) for k in ("engine", "modes", "cluster_top")})
check(any("ADCP 1.0 (macOS)" in l for l in lines),
      "the task's log does not name the engine: %s" % lines[:2])
pipeline.stage_report(rd, jobmod.load(rd), log=lambda *a: None)
txt = open(os.path.join(rd, "results.txt")).read()
rj = json.load(open(os.path.join(rd, "results.json")))
check(txt.startswith(mac10.NOTE) and rj.get("engine_note") == mac10.NOTE
      and all(r.get("engine") == "ADCP 1.0 (macOS)" for r in rj["rows"]),
      "the report does not say the run docked with ADCP 1.0 on a Mac, or a "
      "row does not name it")
check("not comparable with runs docked with ADCP 1.1" in mac10.NOTE
      and "no ADCP 1.1 for macOS" in mac10.NOTE, "the note does not say why")
app = server.App(TMP)
st = app.run_status(rd)
check(st.get("engine_note") == mac10.NOTE, "the run's page is not told")
check(app.tools().get("mac10", {}).get("engaged") is True,
      "the page is not told this Mac docks with ADCP 1.0")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("d.engine_note" in appjs and "d.mac10.note" in appjs,
      "the page does not show the note")

# --- a run docked with ADCP 1.0 before the move to 1.1, and converted since --------
from adcp_dock import report as reportmod  # noqa: E402

converted = {"converted": {"date": "2026-10-18", "from": "ADCP 1.0"}}
said = reportmod.engine_note([{"engine": "ADCP 1.0"}], converted)
check(said.startswith("Docked with ADCP 1.0, before this tool moved to ADCP 1.1")
      and "converted on 2026-10-18" in said and "not comparable" in said,
      "a converted 1.0 run's page does not say its scores are 1.0's: %r" % said)
check(reportmod.engine_note([{"engine": "ADCP 1.0"}]).startswith("Docked with ADCP 1.0"),
      "the rows alone do not say a run was docked with ADCP 1.0")
check(reportmod.engine_note([{"engine": "ADCP 1.1.21"}], {}) == "",
      "an ADCP 1.1 run is given a note")
check(reportmod.engine_note([], {"adcp_engine": mac10.ENGINE}) == mac10.NOTE,
      "a Mac run's job alone does not give the Mac note")

# --- the install: the macOS tarball, checked ---------------------------------------
# A stand-in tarball whose install.sh makes the three programs in -d.
src = os.path.join(TMP, "src", "ADFRsuite_x86_64Darwin_1.0")
os.makedirs(src)
open(os.path.join(src, "install.sh"), "w").write(
    '#!/bin/sh\nwhile [ $# -gt 0 ]; do [ "$1" = -d ] && D=$2; shift; done\n'
    'mkdir -p "$D/bin"; for b in adcp agfr prepare_receptor; do '
    'printf "#!/bin/sh\\n" > "$D/bin/$b"; chmod 755 "$D/bin/$b"; done\n')
tgz = os.path.join(TMP, "built.tar.gz")
with tarfile.open(tgz, "w:gz") as tf:
    tf.add(src, arcname="ADFRsuite_x86_64Darwin_1.0")
sha = hashlib.sha256(open(tgz, "rb").read()).hexdigest()
_fetch, _build = adcp11_install._fetch, mac10.BUILD
adcp11_install._fetch = lambda url, path, log: shutil.copy(tgz, path)
try:
    mac10.BUILD = (_build[0], _build[1], sha, 1)
    dest = os.path.join(TMP, "inst", "ADFRsuite")
    got = mac10.install(dest, log=lambda *a: None)
    check(got == dest and mac10.installed(dest),
          "ADFRsuite 1.0 did not install from its tarball")
    mac10.BUILD = (_build[0], _build[1], "0" * 64, 1)
    dest2 = os.path.join(TMP, "inst2", "ADFRsuite")
    try:
        mac10.install(dest2, log=lambda *a: None)
        said = ""
    except SystemExit as exc:
        said = str(exc)
    check("pinned checksum" in said and not os.path.exists(dest2)
          and not os.path.exists(os.path.join(TMP, "inst2", _build[1])),
          "a tarball that does not match its checksum was used or kept (%r)" % said)
finally:
    adcp11_install._fetch, mac10.BUILD = _fetch, _build
check(_build[1] == "ADFRsuite_x86_64Darwin_1.0.tar.gz" and len(_build[2]) == 64,
      "the macOS build is not the pinned Darwin tarball")

# --- `adcp-dock local` installs it on a Mac, and only there ------------------------
asked = []


def stop(dest=None, log=print):
    asked.append(dest)
    raise SystemExit("installed")


# The Python modules `adcp-dock local` asks for first, stood in for where
# this environment lacks one: what is checked comes after them.
for mod in ("pdb2pqr", "openmm", "pdbfixer", "rdkit"):
    try:
        __import__(mod)
    except ImportError:
        sys.modules[mod] = types.ModuleType(mod)
RUNS = os.path.join(TMP, "runs")
_install, _avail = mac10.install, cli.slurm.available
mac10.install, cli.slurm.available = stop, (lambda: False)
try:
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        try:
            cli.main(["local", "--no-browser", "--runs", RUNS])
            said = ""
        except SystemExit as exc:
            said = str(exc)
    check(said == "installed" and asked == [mac10.DEFAULT],
          "adcp-dock local on a Mac does not install ADFRsuite 1.0 into %s "
          "(%r, %s)" % (mac10.DEFAULT, said, asked))
    as_system("Linux")
    asked.clear()
    # A Linux machine with the ADCP 1.1 suite: the server is reached, and
    # ADFRsuite 1.0 is not.
    _find, _at, _build_srv = adcp11.find, adcp11.at, server.build
    adcp11.find = lambda roots=(): os.path.join(TMP, "adcpsuite", "bin", "adcp")
    adcp11.at = lambda home: None

    def no_server(*a, **k):
        raise SystemExit("server")

    server.build = no_server
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            try:
                cli.main(["local", "--no-browser", "--runs", RUNS])
                said = ""
            except SystemExit as exc:
                said = str(exc)
    finally:
        adcp11.find, adcp11.at, server.build = _find, _at, _build_srv
    check(said == "server" and not asked,
          "adcp-dock local on Linux reached ADFRsuite 1.0 (%r, %s)" % (said, asked))
finally:
    mac10.install, cli.slurm.available = _install, _avail
    mac10._system = _real_system

if fails:
    print("test_mac10 FAILED (%d of %d)" % (len(fails), ran))
    sys.exit(1)
print("ok  ADCP 1.0 on a Mac and nowhere else: the gate, the tools, the maps, "
      "the docking in the 1.1 layout, the report and the install, %d checks"
      % ran)
