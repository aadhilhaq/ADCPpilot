"""An installation with the ADCP 1.1 suite, and nothing else for the docking.

By Aadhil Haq

install.sh installs the ADCP 1.1 suite, and everything the docking needs
comes from it. Checked here, with nothing else in sight (no tools roots, a
bare PATH): the required tools are found in the suite, prepare_receptor and
reduce are the suite's, MSMS is its library and draws the surface the msms
program drew (tests/data/msms_1YCR.json, made with the program and kept as
the reference), the job scripts name the suite without putting its bin
folder (and its python) on PATH, prepare_receptor is started with the
suite's own Python, the suite is found beside a checkout with nothing named,
and install.sh decides what to install as it says it does.

    python3 tests/test_adcp11_only.py

Without the ADCP 1.1 suite installed (ADCP_TEST_TOOLS or ADCPSUITE_HOME),
nothing can be checked and it exits 77.
"""
import base64
import hashlib
import json
import os
import shutil
import struct
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_pilot import adcp11, slurm, surface, tools  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)
        print("  FAIL " + why)


roots = [r for r in (os.environ.get("ADCP_TEST_TOOLS", ""),
                     os.environ.get("ADCPSUITE_HOME", ""))
         if r and os.path.isdir(r)]
exe = adcp11.find(roots)
if not exe:
    print("skip  the ADCP 1.1 suite is not installed here (ADCP_TEST_TOOLS or "
          "ADCPSUITE_HOME)")
    sys.exit(77)
SUITE = adcp11.home_of(exe)
EX = os.path.join(ROOT, "adcp_pilot", "example_run", "receptor", "1YCR.pdb")
REF = json.load(open(os.path.join(HERE, "data", "msms_1YCR.json")))

saved = dict(os.environ)
os.environ["ADCPSUITE_HOME"] = SUITE
here_bin = os.path.dirname(os.path.abspath(sys.executable))
os.environ["PATH"] = os.pathsep.join([here_bin, "/usr/bin", "/bin"])
try:
    T = tools.discover()
    check(T.missing_required() == [] or T.missing_required() == ["pdb2pqr"],
          "the suite alone misses required tools: %s" % T.missing_required())
    for name in ("adcp", "agfr", "prepare_receptor"):
        check(T.get(name) == os.path.join(SUITE, "bin", name),
              "%s is not the suite's: %r" % (name, T.get(name)))
    pr = T.get("prepare_receptor") or ""
    check(T.suite_program("reduce") == os.path.join(SUITE, "bin", "reduce"),
          "reduce is not the suite's: %r" % T.suite_program("reduce"))
    m = tools.msms_command(T)
    check(isinstance(m, list) and m[0] == os.path.join(SUITE, "bin", "python")
          and m[1].endswith("msms_mslib.py"),
          "MSMS is not the suite's library: %r" % m)
    check("ADCP 1.1 suite" in surface.describe(m), "the page is not told where MSMS is")
    check(tools.suite_home() == SUITE, "the suite is not found by name")
    # Another copy of the programs on PATH is never taken for the suite's.
    fake = tempfile.mkdtemp(prefix="adcp11only_path_")
    for name in ("adcp", "agfr", "prepare_receptor", "reduce"):
        with open(os.path.join(fake, name), "w") as fh:
            fh.write("#!/bin/sh\nexit 0\n")
        os.chmod(os.path.join(fake, name), 0o755)
    os.environ["PATH"] = fake + os.pathsep + os.environ["PATH"]
    T2 = tools.discover()
    check(all(T2.get(n) == os.path.join(SUITE, "bin", n)
              for n in ("adcp", "agfr", "prepare_receptor"))
          and T2.suite_program("reduce") == os.path.join(SUITE, "bin", "reduce"),
          "a program on PATH was taken for the suite's: %s"
          % [T2.get(n) for n in ("adcp", "agfr", "prepare_receptor")])
    os.environ["PATH"] = os.pathsep.join([here_bin, "/usr/bin", "/bin"])
    shutil.rmtree(fake, ignore_errors=True)
    # The job scripts name the suite, and keep its bin folder off PATH.
    lines = slurm.tool_exports().splitlines()
    check(any(l.startswith("export ADCPSUITE_HOME=") for l in lines),
          "a job script does not name the ADCP 1.1 suite: %s" % lines)
    check(not any(l.startswith("export PATH=") and os.path.join(SUITE, "bin") in l
                  for l in lines), "a job script puts the suite's python on PATH")
    # prepare_receptor, started as the receptor module starts it, with the
    # suite's own Python whatever is first on PATH.
    from adcp_pilot._vendor import receptor as rec
    tmp = tempfile.mkdtemp(prefix="adcp11only_")
    src = os.path.join(tmp, "recA.pdb")
    with open(EX) as fh, open(src, "w") as out:
        out.writelines(l for l in fh if l.startswith("ATOM") and l[21:22] == "A")
        out.write("END\n")
    p = rec._run_prepare_receptor(pr, src, os.path.join(tmp, "recA.pdbqt"), add_h=True)
    made = os.path.join(tmp, "recA.pdbqt")
    check(p.returncode == 0 and os.path.isfile(made) and os.path.getsize(made) > 10000,
          "the suite's prepare_receptor did not write a receptor: %s %s"
          % (p.returncode, (p.stdout + p.stderr)[-300:]))
    # The library draws the surface the msms program drew, vertex for vertex:
    # the same points, triangles and residue per point, compared through a
    # digest of them in a fixed order (the library lists them in another
    # order, which nothing that draws them minds).
    atoms = surface.read_text(open(EX).read(), "pdb")
    lib = surface.compute(atoms, m)

    def dec(v, f):
        r = base64.b64decode(v)
        return struct.unpack("<%d%s" % (len(r) // 4, f), r)

    def shape(d):
        v = dec(d["verts"], "f")
        pts = [tuple(round(x, 2) for x in v[i:i + 3]) for i in range(0, len(v), 3)]
        f = dec(d["faces"], "I")
        tris = sorted(tuple(sorted(pts[f[i + k]] for k in range(3)))
                      for i in range(0, len(f), 3))
        return sorted(zip(pts, dec(d["res"], "I"))), tris

    check(lib.get("vertices") == REF["vertices"]
          and lib.get("triangles") == REF["triangles"]
          and lib.get("density") == REF["density"],
          "mslib's surface of 1YCR has %s vertices and %s triangles; the msms "
          "program's had %d and %d" % (lib.get("vertices"), lib.get("triangles"),
                                       REF["vertices"], REF["triangles"]))
    got = hashlib.sha256(json.dumps(shape(lib)).encode()).hexdigest()
    check(got == REF["shape_sha256"],
          "mslib's surface of 1YCR is not the msms program's, point for point "
          "(%s, %s)" % (got[:12], REF["shape_sha256"][:12]))
    shutil.rmtree(tmp, ignore_errors=True)

    # Found with nothing named: beside a checkout, as install.sh leaves it
    # on a cluster (<folder>/ADCPpilot and <folder>/tools/docking/adcpsuite).
    os.environ.pop("ADCPSUITE_HOME", None)
    side = tempfile.mkdtemp(prefix="adcp11only_side_")
    os.makedirs(os.path.join(side, "tools", "docking"))
    os.symlink(SUITE, os.path.join(side, "tools", "docking", "adcpsuite"))
    os.makedirs(os.path.join(side, "ADCPpilot"))
    os.environ["ADCP_HOME"] = os.path.join(side, "ADCPpilot")
    found = adcp11.find()
    check(found and os.path.realpath(adcp11.home_of(found)) == os.path.realpath(SUITE),
          "the suite beside the installation is not found with nothing named: %r"
          % found)
    shutil.rmtree(side, ignore_errors=True)
finally:
    os.environ.clear()
    os.environ.update(saved)

# install.sh: what it would install into an empty folder. --check changes
# nothing.
empty = tempfile.mkdtemp(prefix="adcp11only_prefix_")
os.rmdir(empty)
env = {"PATH": "/usr/bin:/bin", "HOME": os.environ.get("HOME", "/tmp")}
# A checkout on the lab's cluster sits beside tools/docking, whose suite and
# ADapt are found where a test says there are none; a copy in a folder with
# nothing beside it is what "not installed" means here.
ALONE = os.path.join(tempfile.mkdtemp(prefix="adcp11only_alone_"), "ADCPpilot")
shutil.copytree(ROOT, ALONE, ignore=shutil.ignore_patterns(
    ".git", "tests", "build", "*.egg-info", "__pycache__", "notebooks", "tools",
    "activate.sh", "example_run", "_runs", "_work"))


def would(*opts, extra=None):
    p = subprocess.run(["bash", os.path.join(ALONE, "install.sh"), "--check"]
                       + list(opts), capture_output=True, text=True,
                       env=dict(env, **(extra or {})), timeout=600)
    return " ".join((p.stdout + p.stderr).split())


out = would(empty)
check("the ADCP 1.1 suite (adcp, agfr, prepare_receptor, reduce and the MSMS "
      "library), into %s" % os.path.join(empty, "tools", "docking", "adcpsuite") in out,
      "install.sh does not install the ADCP 1.1 suite into the folder given:\n"
      + out[-600:])
out = would(empty, extra={"ADCPSUITE_HOME": SUITE})
check(("USING the ADCP 1.1 suite at %s" % os.path.realpath(SUITE)) in out
      and "MSMS library), into" not in out,
      "install.sh does not adopt the suite ADCPSUITE_HOME names:\n" + out[-600:])
for gone in ("--with-adfr10", "--no-adcp11"):
    p = subprocess.run(["bash", os.path.join(ROOT, "install.sh"), "--check", gone],
                       capture_output=True, text=True, env=env, timeout=120)
    check(p.returncode == 2 and "there is no option" in p.stderr,
          "install.sh still takes %s" % gone)
check(not os.path.exists(empty), "install.sh --check made the folder")
src = open(os.path.join(ROOT, "install.sh")).read()
check("exit 1" in src.split("The ADCP 1.1 suite did not install")[1][:600],
      "install.sh carries on when the ADCP 1.1 suite fails to install")
check('echo "export ADCPSUITE_HOME=' in src and 'echo "export ADAPT_HOME=' in src,
      "activate.sh does not name the ADCP 1.1 suite and ADapt")
appjs = open(os.path.join(ROOT, "adcp_pilot", "static", "app.js")).read()
check("ADCP 1.1 is not installed here,</b>" in appjs
      and "so nothing can dock" in appjs,
      "the page does not say plainly when the ADCP 1.1 suite is missing")

if fails:
    print("test_adcp11_only FAILED (%d of %d)" % (len(fails), ran))
    sys.exit(1)
print("ok  an installation with the ADCP 1.1 suite: tools, surfaces, job "
      "scripts and install.sh, %d checks" % ran)
