"""install.sh with environments it did not build, and AmberTools from a module.

By Aadhil Haq

Stand-ins for everything: an environment whose python answers imports from
a list of installed modules (and fails every import while a module's
PYTHONPATH is set, as a real conda Python does), a micromamba that records
what it was asked and "installs" by adding to that list, an AmberTools whose
MMPBSA.py runs only with its module's PYTHONPATH, and an ADCP 1.1 suite, so
nothing is downloaded. Against those:

  * a module's PYTHONPATH does not make install.sh think numpy is missing;
  * AmberTools from PATH or a module is used, through wrappers that work
    with the module unloaded, and no AmberTools is added to the environment;
  * what a writable environment lacks is added beside what is there
    (--freeze-installed), and only that;
  * one you cannot write to is left alone, with what is missing and what
    that turns off; one that lacks RDKit and cannot be given it is not used;
  * nothing is pip-installed into an environment install.sh did not build;
    the adcp-pilot command comes from the checkout instead;
  * activate.sh replaces PYTHONPATH, and there is no /opt/shared default;
  * install.sh and update.sh refuse to be sourced or run by sh;
  * update.sh reinstalls only into an environment that imports this checkout.

    python3 tests/test_install.py
"""
import os
import shutil
import stat
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


if not shutil.which("bash") or not shutil.which("git"):
    print("skip  bash or git is not available")
    sys.exit(77)

tmp = tempfile.mkdtemp(prefix="adcp_install_")
home = os.path.join(tmp, "home")
os.makedirs(home)


def write(path, text, mode=0o755):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, mode)


# --- a checkout: the scripts, and the package linked in -----------------------
dep = os.path.join(tmp, "ADCP_docking")
os.makedirs(os.path.join(dep, "scripts"))
for f in ("install.sh", "update.sh"):
    shutil.copy(os.path.join(ROOT, f), dep)
shutil.copy(os.path.join(ROOT, "scripts", "adapt_fetch.sh"), os.path.join(dep, "scripts"))
os.symlink(os.path.join(ROOT, "adcp_pilot"), os.path.join(dep, "adcp_pilot"))

# --- an ADCP 1.1 suite, so nothing is installed there ---------------------------
tools = os.path.join(tmp, "tools", "docking")
suite = os.path.join(tools, "adcpsuite")
write(os.path.join(suite, "bin", "adcp"), "#!/bin/sh\nexit 0\n")
write(os.path.join(suite, "lib", "python3.7", "site-packages", "ADCP", "runADCP.py"), "", 0o644)

# --- a stand-in Python environment ----------------------------------------------
# bin/python answers `-c "import X ..."` from installed.txt beside it, and fails
# every import while PYTHONPATH holds the module's site folder, which is what a
# conda Python does when a module's Python 3.9 packages come first.
MODULE_SITE = os.path.join(tmp, "module", "lib", "python3.9", "site-packages")
FAKE_PY = r'''#!/bin/bash
ENV="$(cd "$(dirname "$0")/.." && pwd)"
case ":${PYTHONPATH:-}:" in *":%(site)s:"*) BROKEN=1 ;; *) BROKEN=0 ;; esac
if [ "$1" = -V ]; then echo "Python 3.11.9"; exit 0; fi
if [ "$1" = -c ]; then
  set -- $2
  mod=${2%%[;,]*}
  case "$mod" in
    sys) exit 0 ;;
    sysconfig) echo "$ENV/lib/python3.11/site-packages"; exit 0 ;;
  esac
  [ "$BROKEN" = 1 ] && { echo "ImportError: numpy.core._multiarray_umath" >&2; exit 1; }
  grep -qx "$mod" "$ENV/installed.txt" && exit 0
  echo "ModuleNotFoundError: No module named '$mod'" >&2; exit 1
fi
if [ "$1" = -m ] && [ "$2" = pip ]; then
  echo "pip $*" >> "$ENV/pip.log"
  echo "PIP_NO_CACHE_DIR=${PIP_NO_CACHE_DIR:-} XDG_CACHE_HOME=${XDG_CACHE_HOME:-} HOME=$HOME" >> "$ENV/pip.env"
  if [ "$3" = install ]; then
    for a in "$@"; do case "$a" in
      openmm|pdbfixer|pdb2pqr|rdkit) echo "$a" >> "$ENV/installed.txt" ;;
      biopython) echo Bio >> "$ENV/installed.txt" ;;
    esac; done
    for a in "$@"; do [ "$a" = pdb2pqr ] && printf '#!/bin/sh\nexit 0\n' > "$ENV/bin/pdb2pqr" && chmod 755 "$ENV/bin/pdb2pqr"; done
  fi
  exit 0
fi
exit 0
''' % {"site": MODULE_SITE}


def make_env(name, modules, conda=True, amber=False, pdb2pqr=True):
    e = os.path.join(tmp, name)
    write(os.path.join(e, "bin", "python"), FAKE_PY)
    os.symlink("python", os.path.join(e, "bin", "python3"))
    os.makedirs(os.path.join(e, "lib", "python3.11", "site-packages"))
    write(os.path.join(e, "installed.txt"), "".join(m + "\n" for m in modules), 0o644)
    if conda:
        os.makedirs(os.path.join(e, "conda-meta"))
    else:
        write(os.path.join(e, "pyvenv.cfg"), "home = /usr/bin\n", 0o644)
    if pdb2pqr:
        write(os.path.join(e, "bin", "pdb2pqr"), "#!/bin/sh\nexit 0\n")
    if amber:
        for b in ("tleap", "sander", "MMPBSA.py", "ante-MMPBSA.py"):
            write(os.path.join(e, "bin", b), "#!/bin/sh\nexit 0\n")
    return e


# --- a micromamba that records what it is asked -------------------------------
solver_log = os.path.join(tmp, "solver.log")
write(os.path.join(dep, "tools", "micromamba"), r'''#!/bin/bash
echo "$*" >> "%(log)s"
[ "$1" = --version ] && { echo 2.0.0; exit 0; }
if [ "$1" = install ]; then
  env=""; prev=""
  for a in "$@"; do [ "$prev" = -p ] && env=$a; prev=$a; done
  # What a real one leaves: the packages it fetched, in the cache it is given.
  mkdir -p "$CONDA_PKGS_DIRS/somepkg-1.0/info"
  echo x > "$CONDA_PKGS_DIRS/somepkg-1.0/info/index.json"
  [ -f "$env/softlink" ] && ln -sf "$CONDA_PKGS_DIRS/somepkg-1.0/info/index.json" "$env/linked"
  for a in "$@"; do case "$a" in
    openmm|pdbfixer|pdb2pqr|rdkit) echo "$a" >> "$env/installed.txt" ;;
    biopython) echo Bio >> "$env/installed.txt" ;;
    ambertools) for b in tleap sander MMPBSA.py ante-MMPBSA.py; do
                  printf '#!/bin/sh\nexit 0\n' > "$env/bin/$b"; chmod 755 "$env/bin/$b"; done ;;
  esac; done
fi
exit 0
''' % {"log": solver_log})

# --- AmberTools as a module provides it -------------------------------------------
# MMPBSA.py runs only with the module's PYTHONPATH, as a module's does.
amber_home = os.path.join(tmp, "module", "AmberTools")
for b in ("tleap", "sander", "ante-MMPBSA.py", "ambpdb"):
    write(os.path.join(amber_home, "bin", b), "#!/bin/sh\necho \"$0 AMBERHOME=$AMBERHOME\"\n")
write(os.path.join(amber_home, "bin", "README"), "not a program\n", 0o644)
write(os.path.join(amber_home, "bin", "MMPBSA.py"),
      '#!/bin/sh\ncase ":${PYTHONPATH:-}:" in *":%s:"*) echo "usage: MMPBSA.py"; exit 0 ;; esac\n'
      'echo "ImportError: No module named MMPBSA_mods" >&2; exit 1\n' % MODULE_SITE)

base_env = {k: v for k, v in os.environ.items()
            if not k.startswith(("ADCP", "ADAPT", "AMBER", "CONDA", "MAMBA",
                                 "PYTHON", "LOADEDMODULES", "SLURM"))}
base_env.update(HOME=home, PATH="/usr/bin:/bin", ADCP_IGNORE_ROOM="1",
                ADCPSUITE_HOME=suite, ADAPT_REPO_URL=os.path.join(tmp, "no-adapt.git"),
                GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1")
module_env = dict(base_env, AMBERHOME=amber_home, PYTHONPATH=MODULE_SITE,
                  PATH=os.path.join(amber_home, "bin") + ":/usr/bin:/bin",
                  LD_LIBRARY_PATH=os.path.join(amber_home, "lib"),
                  LOADEDMODULES="GCC/11.2.0:AmberTools/22.3")


def install(env_vars, *args):
    p = subprocess.run(["bash", os.path.join(dep, "install.sh"), *args], cwd=dep,
                       env=env_vars, stdin=subprocess.DEVNULL, capture_output=True,
                       text=True, timeout=300)
    return p.returncode, " ".join((p.stdout + p.stderr).split())


def reset():
    for p in (os.path.join(dep, "activate.sh"), solver_log):
        if os.path.exists(p):
            os.remove(p)
    for d in ("amber", "bin", "conda_shared"):
        shutil.rmtree(os.path.join(dep, "tools", d), ignore_errors=True)


def installed(e):
    return open(os.path.join(e, "installed.txt")).read().split()


def activate():
    p = os.path.join(dep, "activate.sh")
    return open(p).read() if os.path.exists(p) else ""


# --- 1. a module's PYTHONPATH, and its AmberTools -------------------------------
e1 = make_env("env_noomm", ["numpy", "rdkit", "Bio"])
reset()
rc, out = install(dict(module_env, ADCP_ENV=e1), "--check")
check(rc == 0 and ("FOUND %s" % e1) in out and "lacks: numpy" not in out,
      "a module's PYTHONPATH made the environment look broken (exit %d):\n%s" % (rc, out[-900:]))
check("USING the cluster's own, from module AmberTools/22.3" in out,
      "AmberTools from the loaded module was not used:\n" + out[-900:])
check("OpenMM, PDBFixer, added to" in out and "AmberTools, added" not in out,
      "--check did not plan to add only OpenMM and PDBFixer:\n" + out[-900:])
check(not os.path.exists(solver_log) and not os.path.exists(os.path.join(dep, "activate.sh")),
      "install.sh --check changed something")

rc, out = install(dict(module_env, ADCP_ENV=e1))
check(rc == 0, "the install with a module's AmberTools failed (exit %d):\n%s" % (rc, out[-1500:]))
log = open(solver_log).read() if os.path.exists(solver_log) else ""
check("--freeze-installed" in log and ("-p %s" % e1) in log and "openmm" in log
      and "pdbfixer" in log and "ambertools" not in log and "create" not in log,
      "what was asked of the solver is not an addition of OpenMM and PDBFixer alone:\n" + log)
check(set(installed(e1)) >= {"openmm", "pdbfixer"}, "OpenMM and PDBFixer were not added")
check(not os.path.exists(os.path.join(e1, "pip.log")),
      "pip was run in an environment install.sh did not build")
wrapper = os.path.join(dep, "tools", "amber", "bin", "MMPBSA.py")
check(os.path.isfile(wrapper), "no AmberTools wrappers were written")
if os.path.isfile(wrapper):
    p = subprocess.run([wrapper, "--help"], env=dict(base_env), capture_output=True,
                       text=True, timeout=60)
    check(p.returncode == 0 and "usage" in p.stdout,
          "the MMPBSA.py wrapper does not run with the module unloaded:\n" + p.stdout + p.stderr)
    p = subprocess.run([os.path.join(dep, "tools", "amber", "bin", "tleap")],
                       env=dict(base_env), capture_output=True, text=True, timeout=60)
    check(("AMBERHOME=%s" % amber_home) in p.stdout,
          "the tleap wrapper does not give tleap its own AMBERHOME:\n" + p.stdout)
    # The rescoring runs more than the four programs checked for: the panel
    # writes the scored geometry with $AMBERHOME/bin/ambpdb.
    wbin = os.path.join(dep, "tools", "amber", "bin")
    check(os.access(os.path.join(wbin, "ambpdb"), os.X_OK)
          and not os.path.exists(os.path.join(wbin, "README")),
          "the wrappers do not cover every program in the module's bin folder: %r"
          % sorted(os.listdir(wbin)))
act = activate() if os.path.exists(os.path.join(dep, "activate.sh")) else ""
check('export AMBERHOME="%s"' % os.path.join(dep, "tools", "amber") in act,
      "activate.sh does not point AMBERHOME at the wrappers:\n" + act)
check('export PYTHONPATH="$ADCP_HOME"' in act and "${PYTHONPATH" not in act,
      "activate.sh adds to PYTHONPATH instead of replacing it:\n" + act)
check(os.access(os.path.join(dep, "tools", "bin", "adcp-pilot"), os.X_OK)
      and os.path.join(dep, "tools", "bin") in act,
      "the adcp-pilot command was not provided from the checkout:\n" + act)
check(os.listdir(home) == [], "install.sh wrote into the home directory: %r" % os.listdir(home))

# Run again with the module unloaded: the wrappers are found and kept.
rc, out = install(dict(base_env, ADCP_ENV=e1,
                       PATH=os.path.join(dep, "tools", "amber", "bin") + ":/usr/bin:/bin"),
                  "--check")
check(rc == 0 and "the wrappers an earlier run wrote" in out,
      "a second run does not recognize its own AmberTools wrappers:\n" + out[-900:])

# --- 2. AmberTools in neither: added to a writable conda environment -------------
e2 = make_env("env_noamber", ["numpy", "rdkit", "openmm", "pdbfixer", "Bio"])
reset()
rc, out = install(dict(base_env, ADCP_ENV=e2))
own_pkgs = os.path.join(dep, "tools", "conda_shared", "pkgs")
check("removed the download cache it was built from" in out and not os.path.exists(own_pkgs),
      "the download cache install.sh made was not removed after adding to the environment:\n"
      + out[-900:])
log = open(solver_log).read() if os.path.exists(solver_log) else ""
check(rc == 0 and "ambertools" in log and "--freeze-installed" in log,
      "AmberTools was not added to a writable environment that lacked it (exit %d):\n%s\n%s"
      % (rc, log, out[-900:]))
check('export AMBERHOME="%s"' % e2 in (activate() if rc == 0 else ""),
      "activate.sh does not name the AmberTools just added")

# The cache is kept when the environment links into it, and one named by
# whoever ran this is never touched.
e2b = make_env("env_softlinked", ["numpy", "rdkit", "openmm", "pdbfixer", "Bio"])
write(os.path.join(e2b, "softlink"), "", 0o644)
reset()
rc, out = install(dict(base_env, ADCP_ENV=e2b))
check(rc == 0 and "the download cache is kept" in out and os.path.isdir(own_pkgs)
      and os.path.isfile(os.path.join(e2b, "linked"))
      and open(os.path.join(e2b, "linked")).read() == "x\n",
      "a cache the environment links into was removed (exit %d):\n%s" % (rc, out[-900:]))
users = os.path.join(tmp, "users_conda_pkgs")
write(os.path.join(users, "theirs-1.0", "info", "index.json"), "theirs\n", 0o644)
e2c = make_env("env_userscache", ["numpy", "rdkit", "openmm", "pdbfixer", "Bio"])
reset()
rc, out = install(dict(base_env, ADCP_ENV=e2c, CONDA_PKGS_DIRS=users))
check(rc == 0 and os.path.isfile(os.path.join(users, "theirs-1.0", "info", "index.json"))
      and "removed the download cache" not in out,
      "a package cache named by the person running install.sh was cleaned (exit %d):\n%s"
      % (rc, out[-900:]))

# --- 3. one you cannot write to ----------------------------------------------------
e3 = make_env("env_readonly", ["numpy", "rdkit"])
reset()
for root, dirs, files in os.walk(e3):
    for n in dirs + files:
        p = os.path.join(root, n)
        if not os.path.islink(p):
            os.chmod(p, os.stat(p).st_mode & ~0o222)
os.chmod(e3, os.stat(e3).st_mode & ~0o222)
rc, out = install(dict(base_env, ADCP_ENV=e3))
check(rc == 0 and "which you cannot write to" in out
      and "the receptor's restrained minimization is off" in out
      and "MM-GBSA rescoring is off" in out,
      "a read-only environment was not reported with what is off (exit %d):\n%s"
      % (rc, out[-1500:]))
check(not os.path.exists(solver_log), "the solver was run on a read-only environment")
check(installed(e3) == ["numpy", "rdkit"], "a read-only environment was changed")
check("No AmberTools was found, so MM-GBSA rescoring is off" in activate(),
      "activate.sh does not say MM-GBSA is off")

# --- 4. lacking RDKit, read-only, and no folder: not used --------------------------
e4 = make_env("env_nordkit", ["numpy"])
os.chmod(os.path.join(e4, "conda-meta"), 0o555)
reset()
rc, out = install(dict(base_env, ADCP_ENV=e4))
check(rc == 1 and "which adcp-pilot cannot run without" in out and "./install.sh" in out
      and "with a folder that has room" in out,
      "an environment without RDKit that cannot be given it was not set aside (exit %d):\n%s"
      % (rc, out[-1200:]))

# --- 5. a pip environment (a venv) lacking pdb2pqr and OpenMM ------------------------
e5 = make_env("venv", ["numpy", "rdkit", "pdbfixer", "Bio"], conda=False, pdb2pqr=False)
reset()
rc, out = install(dict(module_env, ADCP_ENV=e5))
pip_log = open(os.path.join(e5, "pip.log")).read() if os.path.exists(os.path.join(e5, "pip.log")) else ""
check(rc == 0 and "pdb2pqr" in pip_log and "openmm" in pip_log and " -c " in pip_log
      and "ambertools" not in pip_log,
      "a venv was not given pdb2pqr and OpenMM with pip, held by constraints (exit %d):\n%s\n%s"
      % (rc, pip_log, out[-1200:]))
check(not os.path.exists(solver_log), "conda was used on a venv")
pip_env = open(os.path.join(e5, "pip.env")).read() if os.path.exists(os.path.join(e5, "pip.env")) else ""
check(pip_env and all("PIP_NO_CACHE_DIR=1 " in l and ("XDG_CACHE_HOME=%s" % home) not in l
                      and "XDG_CACHE_HOME= " not in l for l in pip_env.splitlines()),
      "pip was run without PIP_NO_CACHE_DIR=1 and a cache folder outside home:\n" + pip_env)
# And a real pip, given the same settings, keeps no cache anywhere.
p = subprocess.run([sys.executable, "-m", "pip", "cache", "dir"], capture_output=True, text=True,
                   env=dict(base_env, PIP_NO_CACHE_DIR="1"), timeout=120)
check("disabled" in (p.stdout + p.stderr).lower() or p.returncode != 0,
      "pip still names a cache with PIP_NO_CACHE_DIR=1: " + (p.stdout + p.stderr))
check(os.listdir(home) == [], "something wrote into the home directory: %r" % os.listdir(home))

# --- the script itself ---------------------------------------------------------------
sh = open(os.path.join(ROOT, "install.sh")).read()
check("/opt/shared" not in sh, "install.sh still looks in /opt/shared by default")
for script in ("install.sh", "update.sh"):
    p = subprocess.run(["bash", "-c", "source %s; echo after=$?" % os.path.join(dep, script)],
                       cwd=dep, env=dict(base_env), capture_output=True, text=True, timeout=60)
    check("do not source it" in p.stdout + p.stderr and "after=2" in p.stdout,
          "sourcing %s was not refused:\n%s" % (script, p.stdout + p.stderr))
    if shutil.which("sh") and os.path.realpath(shutil.which("sh")) != os.path.realpath(
            shutil.which("bash")):
        p = subprocess.run(["sh", os.path.join(dep, script)], cwd=dep, env=dict(base_env),
                           capture_output=True, text=True, timeout=60)
        check(p.returncode == 2 and "Run this with bash" in p.stderr,
              "running %s with sh was not refused:\n%s" % (script, p.stderr))
rc, out = install(dict(base_env), "--bogus")
check(rc == 2 and "there is no option --bogus" in out and "Usage" in out,
      "an unknown option was not explained:\n" + out[-400:])

# --- update.sh leaves an environment it did not build alone -------------------------
subprocess.run(["git", "init", "-q", dep], env=base_env, check=True)
reset()
write(os.path.join(dep, "activate.sh"),
      'export ADCP_HOME="%s"\nexport ADCP_PYTHON="%s"\n' % (dep, os.path.join(e1, "bin", "python3")),
      0o644)
if os.path.exists(os.path.join(e1, "pip.log")):
    os.remove(os.path.join(e1, "pip.log"))
p = subprocess.run(["bash", os.path.join(dep, "update.sh")], cwd=dep, env=dict(base_env),
                   stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=300)
said = p.stdout + p.stderr
check("nothing to reinstall" in said and not os.path.exists(os.path.join(e1, "pip.log")),
      "update.sh ran pip in an environment that does not import this checkout:\n" + said[-900:])

for root, dirs, files in os.walk(tmp):
    for n in dirs:
        try:
            os.chmod(os.path.join(root, n), 0o755)
        except OSError:
            pass
os.chmod(tmp, 0o755)
shutil.rmtree(tmp, ignore_errors=True)

if fails:
    print("test_install FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  install.sh adds only what is missing, uses a module's AmberTools through "
      "wrappers, and leaves environments it did not build alone, %d checks" % ran)
