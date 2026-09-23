"""What a normal `pip install` puts on disk.

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
for top in ("adcp_dock", "adcp_launcher"):
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
    for must in ("adcp_dock/__init__.py", "adcp_dock/_vendor/adcp.py",
                 "adcp_dock/_vendor/prepare.py", "adcp_dock/adfr.py",
                 "adcp_dock/static/app.js", "adcp_dock/static/3Dmol-min.js",
                 "adcp_launcher/__init__.py"):
        check(must in built, "the built wheel does not contain %s" % must)

if fails:
    print("test_packaging FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  every package and data file is in what pip installs, %d checks%s"
      % (ran, "" if built is not None else
         " (no wheel built here: setuptools is not importable)"))
