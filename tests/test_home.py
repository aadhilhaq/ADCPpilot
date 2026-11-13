"""Nothing this tool runs writes into the home folder.

By Aadhil Haq

A cluster's home folder often holds few files (10,000 is common), and three things
here would otherwise fill it: the ADCP 1.1 suite's own Python, whose modules
make ~/.adfrsuite, ~/.data/pybel and ~/.config as they load and know no
other place; matplotlib's font list, in ~/.cache; and the CUDA driver's
compiled kernels, in ~/.nv. Checked here: the suite is given a home folder
of each person's own (adcp11.private_home), and so is ADapt, which runs on
the suite's Python; every job script and a run on this machine keep their
caches in the run folder; the installer's micromamba is given a HOME beside
the installation. With the suite installed, its Python is started for real
with HOME pointing at an empty folder, which must still be empty after.

    python3 tests/test_home.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import adapt, adcp11, slurm  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)
        print("  FAIL " + why)


def files_in(folder):
    out = []
    for d, dirs, names in os.walk(folder):
        out += [os.path.relpath(os.path.join(d, n), folder) for n in dirs + names]
    return sorted(out)


tmp = tempfile.mkdtemp(prefix="adcp_home_")
saved = dict(os.environ)
try:
    user = __import__("getpass").getuser()
    os.environ.pop("HOPE_CACHE", None)

    # --- where the suite's home folder is -------------------------------------
    lab = os.path.join(tmp, "lab")
    suite = os.path.join(lab, "tools", "docking", "adcpsuite")
    os.makedirs(os.path.join(suite, "bin"))
    got = adcp11.private_home(suite)
    check(got == os.path.join(lab, ".cache", user, "adcpsuite-home"),
          "a suite in <lab>/tools/docking is given <lab>/.cache/<user>/"
          "adcpsuite-home, not %s" % got)
    check(os.path.isdir(got), "the suite's home folder is made")
    elsewhere = os.path.join(tmp, "own", "adcpsuite")
    os.makedirs(elsewhere)
    got = adcp11.private_home(elsewhere)
    check(got == os.path.join(tmp, "own", ".cache", user, "adcpsuite-home"),
          "a suite anywhere else keeps its home beside it, not %s" % got)
    os.environ["HOPE_CACHE"] = os.path.join(tmp, "named")
    got = adcp11.private_home(suite)
    check(got == os.path.join(tmp, "named", "adcpsuite-home"),
          "HOPE_CACHE names the folder: %s" % got)
    os.environ.pop("HOPE_CACHE")
    real = os.path.realpath(os.path.expanduser("~"))
    check(not os.path.realpath(adcp11.private_home(suite)).startswith(real + os.sep),
          "the suite's home folder is inside the real home folder")

    # --- the suite's programs, and ADapt, are given it -----------------------
    exe = os.path.join(suite, "bin", "adcp")
    env = adcp11._env(exe)
    check(env["HOME"] == adcp11.private_home(suite),
          "ADCP and AGFR are not started with the suite's own home folder")
    check(adcp11.agfr_env(exe)["HOME"] == env["HOME"],
          "AGFR's environment differs from ADCP's in HOME")
    seen = {}

    class Stop(Exception):
        pass

    def fake_popen(argv, **kw):
        seen.update(kw.get("env") or {})
        raise Stop()

    real_popen = adapt.subprocess.Popen
    adapt.subprocess.Popen = fake_popen
    try:
        adapt.Runner().start(os.path.join(tmp, "adapt"), suite,
                             os.path.join(tmp, "adapt_out"), port=0)
    except Stop:
        pass
    finally:
        adapt.subprocess.Popen = real_popen
    check(seen.get("HOME") == adcp11.private_home(suite),
          "ADapt is not started with the suite's own home folder: %r"
          % seen.get("HOME"))

    # --- job scripts and runs on this machine --------------------------------
    run = os.path.join(tmp, "run")
    os.makedirs(run)
    lines = slurm.cache_exports(run)
    for var in ("XDG_CACHE_HOME", "MPLCONFIGDIR", "CUDA_CACHE_PATH"):
        check(("export %s=%s" % (var, os.path.join(run, ".cache"))) in lines,
              "the job script does not put %s in the run folder" % var)
    from adcp_dock import pipeline  # noqa: E402
    common = pipeline._scripts(run, {"partition": "short"})
    check(lines in common["env"], "the job scripts carry no cache lines")
    text = slurm.SINGLE.format(name="x", partition="", cpus=1, mem="1G",
                               time="00:10:00", extra="", notes="", stage="report",
                               **{k: v for k, v in common.items()
                                  if k in ("workdir", "logdir", "python",
                                           "jobfile", "env")})
    check("export MPLCONFIGDIR=" in text and "export XDG_CACHE_HOME=" in text,
          "a rendered job script keeps matplotlib's cache in the home folder")

    # --- the installer -------------------------------------------------------
    src = open(os.path.join(ROOT, "adcp_dock", "adcp11_install.py")).read()
    check("HOME=solver_home" in src,
          "the ADCP 1.1 installer runs micromamba with the real HOME, which "
          "records the environment in ~/.conda/environments.txt")
    sh = open(os.path.join(ROOT, "install.sh")).read()
    check('HOME="$CACHE_BASE/home"' in sh and "export XDG_CACHE_HOME=" in sh,
          "install.sh runs its solver with the real HOME or XDG cache")

    # --- the suite itself, for real ------------------------------------------
    os.environ.clear()
    os.environ.update(saved)
    roots = [r for r in (os.environ.get("ADCP_TEST_TOOLS", ""),
                         os.environ.get("ADCPSUITE_HOME", ""))
             if r and os.path.isdir(r)]
    found = adcp11.find(roots)
    if not found:
        print("note  the ADCP 1.1 suite is not installed here, so its Python "
              "was not started")
    else:
        sentinel = os.path.join(tmp, "sentinel_home")
        os.makedirs(sentinel)
        os.environ["HOME"] = sentinel
        os.environ["HOPE_CACHE"] = os.path.join(tmp, "cache")
        env = adcp11._env(found)
        p = subprocess.run([os.path.join(adcp11.home_of(found), "bin", "python"),
                            "-c", "import MolKit2, pybel; print('loaded')"],
                           env=env, capture_output=True, text=True, timeout=600,
                           cwd=tmp)
        check("loaded" in p.stdout, "the suite's Python could not load MolKit2 "
              "and pybel: %s" % (p.stdout + p.stderr)[-400:])
        check(files_in(sentinel) == [],
              "the suite's Python wrote into HOME: %s" % files_in(sentinel))
        check(os.path.isdir(os.path.join(env["HOME"], ".adfrsuite")),
              "the suite's settings did not go to its own home folder")
finally:
    os.environ.clear()
    os.environ.update(saved)
    shutil.rmtree(tmp, ignore_errors=True)

if fails:
    print("test_home FAILED (%d of %d)" % (len(fails), ran))
    sys.exit(1)
print("ok  nothing is written into the home folder: the suite's Python, ADapt, "
      "job scripts and the installer, %d checks" % ran)
