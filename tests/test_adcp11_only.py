"""An installation with ADCP 1.1 and without ADFRsuite 1.0.

By Aadhil Haq

install.sh installs the ADCP 1.1 suite alone unless given --with-adfr10, so
everything ADFRsuite 1.0 used to supply has to come from the suite. Checked
here, with ADFRsuite 1.0 out of sight (no ADFRSUITE_HOME, nothing of it on
PATH, no tools roots): the required tools are found, prepare_receptor and
reduce are the suite's, MSMS is its library and draws the same surface as the
msms program, the job scripts name the suite without putting its bin folder
(and its python) on PATH, prepare_receptor is started with the suite's own
Python, and install.sh decides what to install as it says it does.

    python3 tests/test_adcp11_only.py

Without the ADCP 1.1 suite installed, the checks that need it are skipped
(exit 77 when nothing at all could be checked).
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import adcp11, slurm, structure, surface, tools  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)
        print("  FAIL " + why)


roots = [r for r in (os.environ.get("ADCP_TEST_TOOLS", ""),
                     os.path.dirname(os.environ.get("ADFRSUITE_HOME", "").rstrip("/")))
         if r and os.path.isdir(r)]
exe = adcp11.find(roots)
if not exe:
    print("skip  the ADCP 1.1 suite is not installed here (ADCP_TEST_TOOLS)")
    sys.exit(77)
SUITE = adcp11.home_of(exe)
ADFR10 = os.path.join(os.path.dirname(SUITE), "ADFRsuite")
MSMS_PROG = os.path.join(ADFR10, "bin", "msms")
EX = os.path.join(ROOT, "adcp_dock", "example_run", "receptor", "1YCR.pdb")

saved = dict(os.environ)
os.environ.pop("ADFRSUITE_HOME", None)
os.environ["ADCPSUITE_HOME"] = SUITE
here_bin = os.path.dirname(os.path.abspath(sys.executable))
os.environ["PATH"] = os.pathsep.join([here_bin, "/usr/bin", "/bin"])
try:
    T = tools.discover()
    check(T.missing_required() == [] or T.missing_required() == ["pdb2pqr"],
          "a 1.1-only installation misses required tools: %s" % T.missing_required())
    check(not T.get("adcp") and not T.get("agfr"),
          "bin/adcp or bin/agfr is taken for ADFRsuite 1.0's programs")
    pr = T.get("prepare_receptor") or ""
    check(pr == os.path.join(SUITE, "bin", "prepare_receptor"),
          "prepare_receptor is not the suite's: %r" % pr)
    check(T.adfr_program("reduce") == os.path.join(SUITE, "bin", "reduce"),
          "reduce is not the suite's: %r" % T.adfr_program("reduce"))
    m = tools.msms_command(T)
    check(isinstance(m, list) and m[0] == os.path.join(SUITE, "bin", "python")
          and m[1].endswith("msms_mslib.py"),
          "MSMS is not the suite's library: %r" % m)
    check("ADCP 1.1 suite" in surface.describe(m), "the page is not told where MSMS is")
    check(tools.suite11() == SUITE, "the suite is not found by name")
    # The job scripts name the suite, and keep its bin folder off PATH.
    lines = slurm.tool_exports().splitlines()
    check(any(l.startswith("export ADCPSUITE_HOME=") for l in lines),
          "a job script does not name the ADCP 1.1 suite: %s" % lines)
    check(not any(l.startswith("export PATH=") and os.path.join(SUITE, "bin") in l
                  for l in lines), "a job script puts the suite's python on PATH")
    # prepare_receptor, started as the receptor module starts it, with the
    # suite's own Python whatever is first on PATH.
    if os.path.isfile(EX):
        from adcp_dock._vendor import receptor as rec
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
        # The library draws the surface the program draws, vertex for vertex.
        atoms = surface.read_text(open(EX).read(), "pdb")
        lib = surface.compute(atoms, m)
        check(lib.get("vertices") and lib.get("faces"), "mslib drew no surface")
        if os.path.isfile(MSMS_PROG):
            # The same points, triangles and residue per point; the library
            # lists them in another order, which nothing that draws them minds.
            import base64
            import struct

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

            prog = surface.compute(atoms, MSMS_PROG)
            check(prog["vertices"] == lib["vertices"] and shape(prog) == shape(lib),
                  "mslib's surface differs from the msms program's")
        else:
            print("  (ADFRsuite 1.0's msms is not here to compare with)")
        shutil.rmtree(tmp, ignore_errors=True)
finally:
    os.environ.clear()
    os.environ.update(saved)

# install.sh: what it would install into an empty folder, by option. --check
# changes nothing.
empty = tempfile.mkdtemp(prefix="adcp11only_prefix_")
os.rmdir(empty)
env = {"PATH": "/usr/bin:/bin", "HOME": os.environ.get("HOME", "/tmp")}


def would(*opts):
    p = subprocess.run(["bash", os.path.join(ROOT, "install.sh"), "--check"]
                       + list(opts) + [empty], capture_output=True, text=True,
                       env=env, timeout=600)
    return p.stdout + p.stderr


out = would()
check("ADFRsuite 1.0 is not needed" in out and "ADCP 1.1, into" in out,
      "by default install.sh does not install ADCP 1.1 alone:\n" + out[-600:])
out = would("--with-adfr10")
check("ADFRsuite 1.0 (adcp" in out and "ADCP 1.1, into" in out,
      "--with-adfr10 does not add ADFRsuite 1.0 beside ADCP 1.1")
out = would("--no-adcp11")
check("ADFRsuite 1.0 (adcp" in out and "--no-adcp11 leaves it out" in out,
      "--no-adcp11 does not install ADFRsuite 1.0 alone")
check(not os.path.exists(empty), "install.sh --check made the folder")
src = open(os.path.join(ROOT, "install.sh")).read()
check("exit 1" in src.split("ADCP 1.1 did not install")[1][:900],
      "install.sh carries on when ADCP 1.1 fails to install")
check('echo "export ADCPSUITE_HOME=' in src and 'echo "export ADAPT_HOME=' in src,
      "activate.sh does not name ADCP 1.1 and ADapt")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("ADCP 1.1 is not installed here,</b>" in appjs,
      "the page does not say plainly when runs use ADCP 1.0")

if fails:
    print("test_adcp11_only FAILED (%d of %d)" % (len(fails), ran))
    sys.exit(1)
print("ok  an installation with ADCP 1.1 alone: tools, surfaces, job scripts and "
      "install.sh, %d checks" % ran)
