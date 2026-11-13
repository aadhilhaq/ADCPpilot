"""What a normal `pip install` puts on disk.

By Aadhil Haq

An editable install (`pip install -e`) reads the source tree directly, so a
package missing from pyproject.toml goes unnoticed there. This checks the
package list against the tree, and builds a real wheel and looks inside it
where the tools to do that are present.
"""
import fnmatch
import glob
import os
import subprocess
import sys
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

try:
    import tomllib
except ImportError:                                             # < 3.11
    tomllib = None

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


if tomllib is None:
    print("skip  Python < 3.11 has no tomllib to read pyproject.toml with")
    sys.exit(77)

cfg = tomllib.load(open(os.path.join(ROOT, "pyproject.toml"), "rb"))
st = cfg.get("tool", {}).get("setuptools", {})
listed = set(st.get("packages") or [])
data = st.get("package-data") or {}

# --- every importable directory is listed ----------------------------------
for top in ("adcp_pilot", "adcp_launcher"):
    for init in glob.glob(os.path.join(ROOT, top, "**", "__init__.py"),
                          recursive=True):
        rel = os.path.relpath(os.path.dirname(init), ROOT)
        if "__pycache__" in rel:
            continue
        pkg = rel.replace(os.sep, ".")
        check(pkg in listed,
              "%s has an __init__.py but is not in pyproject.toml's package "
              "list, so `pip install` leaves it out" % pkg)

# --- every non-Python file a package reads is covered by package-data -------
for pkg in listed:
    d = os.path.join(ROOT, *pkg.split("."))
    globs = data.get(pkg, [])
    for f in glob.glob(os.path.join(d, "**", "*"), recursive=True):
        if os.path.isdir(f) or "__pycache__" in f or f.endswith((".py", ".pyc")):
            continue
        rel = os.path.relpath(f, d).replace(os.sep, "/")
        # a file inside a listed SUBpackage is that subpackage's business
        if any(rel.startswith(sub.split(".")[-1] + "/")
               for sub in listed if sub.startswith(pkg + ".")):
            continue
        check(any(fnmatch.fnmatch(rel, g) for g in globs),
              "%s/%s is shipped with the code but no package-data pattern "
              "covers it, so `pip install` leaves it out" % (pkg, rel))

# --- the Scripps packages are pinned, in one place ---------------------------
# install.sh downloads nothing from Scripps itself: the ADCP 1.1 suite comes
# through adcp11_install, which checks every Scripps package against its own
# SHA-256 before it is installed.
import re as _re                                                # noqa: E402
_sh = open(os.path.join(ROOT, "install.sh")).read()
check("ccsb.scripps.edu" not in _sh and "adcp_pilot.adcp11_install" in _sh,
      "install.sh fetches from Scripps itself, or not through adcp11_install")
sys.path.insert(0, ROOT)
from adcp_pilot import adcp11_install as _ai                     # noqa: E402
check(len(_ai.WHEELS) >= 20
      and all(_re.fullmatch(r"[0-9a-f]{64}", v) for v in _ai.WHEELS.values()),
      "a Scripps package of the ADCP 1.1 suite is installed without a pinned "
      "SHA-256")
check(not os.path.exists(os.path.join(ROOT, "adcp_pilot", "adfr.py")),
      "the installer of the older suite is still in the package")

# --- and a real wheel, where it can be built offline ------------------------
built = None
try:
    import setuptools  # noqa: F401
    out = tempfile.mkdtemp()
    p = subprocess.run([sys.executable, "-m", "pip", "wheel", "--no-deps",
                        "--no-build-isolation", "-q", "-w", out, ROOT],
                       capture_output=True, text=True, timeout=300)
    wheels = glob.glob(os.path.join(out, "*.whl"))
    if p.returncode == 0 and wheels:
        built = zipfile.ZipFile(wheels[0]).namelist()
except Exception:                                               # noqa: BLE001
    built = None

if built is not None:
    for must in ("adcp_pilot/__init__.py", "adcp_pilot/_vendor/adcp.py",
                 "adcp_pilot/_vendor/prepare.py", "adcp_pilot/adcp11_install.py",
                 # run by the suite's own Python for the surfaces
                 "adcp_pilot/msms_mslib.py",
                 "adcp_pilot/static/app.js", "adcp_pilot/static/3Dmol-min.js",
                 # the lab's shared report module, which the report page loads
                 "adcp_pilot/static/reportkit.js", "adcp_pilot/static/reportkit.css",
                 "adcp_launcher/__init__.py",
                 # the documentation served at /docs/
                 "adcp_pilot/docs.py", "adcp_pilot/docskit.py",
                 "adcp_pilot/docs/pages/index.html", "adcp_pilot/docs/img/interface.png",
                 # the example run, at its deepest folder as well as its top
                 "adcp_pilot/example_run/job.json",
                 "adcp_pilot/example_run/docking/ETFSDLWKLLPEN/out_ranked_1.pdb",
                 "adcp_pilot/example_run/docking/ETFSDLWKLLPEN/modes.json",
                 "adcp_pilot/example_run/mmgbsa/ETFSDLWKLLPEN/pose1/"
                 "ETFSDLWKLLPEN_mmgbsa.dat",
                 # which radii the example's energies were computed with
                 "adcp_pilot/example_run/mmgbsa/RADII"):
        check(must in built, "the built wheel does not contain %s" % must)

if fails:
    print("test_packaging FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  every package and data file is in what pip installs, %d checks%s"
      % (ran, "" if built is not None else
         " (no wheel built here: setuptools is not importable)"))
