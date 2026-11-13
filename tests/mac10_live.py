"""ADCP 1.0 on a real Mac: ADFRsuite 1.0 installed, and a tiny docking.

By Aadhil Haq

Run by the mac10 workflow (.github/workflows/mac10.yml), by hand, on a macOS
runner, and nowhere else: on any other system it stops with 77, a skip. It
is not in tests/run_all.sh, whose tests need no Mac.

    python3 tests/mac10_live.py [folder]      (default ~/.adcp-pilot/ADFRsuite)

What it does, on the Mac it runs on:

  * installs ADFRsuite 1.0 for macOS as `adcp-pilot local` does (mac10.install),
    the tarball checked against its pinned SHA-256;
  * builds grid maps over the example receptor with its agfr (AutoSite 1.0,
    no receptor gradient), and checks the target as a run does;
  * docks one short peptide with its ADCP 1.0, a few replicas of a few
    steps, and checks the folder is in the ADCP 1.1 layout afterwards:
    modes.json, the ranked poses and no loose replica files;
  * says whether ADFRsuite's msms program runs on this macOS.
"""
import json
import os
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_pilot import mac10, surface, target, tools  # noqa: E402

if not mac10.is_mac():
    print("skip: this checks ADCP 1.0 on a Mac, and this is %s"
          % mac10._system())
    sys.exit(77)

dest = os.path.abspath(os.path.expanduser(
    sys.argv[1] if len(sys.argv) > 1 else mac10.DEFAULT))
t0 = time.time()
mac10.install(dest)
os.environ[mac10.ENV] = dest
print("installed in %.0f s: %s" % (time.time() - t0, dest))
assert mac10.engaged() and mac10.home() == dest, "ADFRsuite 1.0 is not engaged"
T = tools.discover()
for name in ("adcp", "agfr", "prepare_receptor"):
    assert T.get(name) == os.path.join(dest, "bin", name), (name, T.get(name))
print("reduce:", tools.suite_program("reduce"))

work = tempfile.mkdtemp(prefix="adcp_mac10_live_")
rec = os.path.join(ROOT, "adcp_pilot", "example_run", "receptor",
                   "receptor.pdbqt")
site = json.load(open(os.path.join(ROOT, "adcp_pilot", "example_run",
                                   "job.json")))["site"]
t0 = time.time()
trg = target.build(rec, site["centre"], [26.0, 26.0, 26.0],
                   os.path.join(work, "target"), T,
                   maps=mac10.maps(gradient=False))
print("maps in %.0f s: %s" % (time.time() - t0, trg))

t0 = time.time()
modes, wd = mac10.dock("ETFSDLWKLLP", trg, os.path.join(work, "docking"),
                       replicas=4, steps=100000, cpus=2)
print("docked in %.0f s: %s" % (time.time() - t0, modes[:3]))
left = sorted(os.listdir(wd))
print("folder:", left)
assert modes and os.path.isfile(os.path.join(wd, "modes.json")), \
    "no modes.json: %s" % open(os.path.join(wd, "adcp.log")).read()[-2000:]
assert os.path.isfile(os.path.join(wd, "out_ranked_1.pdb")), "no ranked pose"
assert not [n for n in left if mac10._LOOSE.match(n) or n in mac10._SCRATCH], \
    "loose replica files were left"

exe = tools.msms_command(T)
said = surface.version(exe)
print("msms: %s: %s" % (exe, said or "does not run on this macOS (the "
                        "viewers draw their own surfaces)"))
print("ok  ADFRsuite 1.0 installed, maps built and a peptide docked with ADCP "
      "1.0 on this Mac, in the ADCP 1.1 layout")
