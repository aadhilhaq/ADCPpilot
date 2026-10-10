"""install.sh and update.sh bring in ADapt, anonymously, and never stop for it.

By Aadhil Haq

A local bare repository stands in for GitHub (ADAPT_REPO_URL points at it),
and a small web server that answers every request with "sign in first" stands
in for the repository while it is private. Against those:

  * update.sh clones ADapt beside ADFRsuite when it is missing and public,
    where the interface's own adapt.find() then finds it;
  * while it is private, or cannot be reached, both scripts say so in one
    line, create nothing and finish normally;
  * update.sh moves a clean checkout forward, refuses one with local changes,
    on another branch or with a history rewritten upstream, and moves one
    back when the newest commit is not an ADapt this interface can embed;
  * a folder that is not a git checkout, or is inside another repository, is
    left alone;
  * no git command asks for anything, even with a terminal to ask on, and a
    stored credential helper or askpass program is never run;
  * install.sh fetches it into a fresh install, leaves an existing folder to
    update.sh, and --check says what an install would do.

    python3 tests/test_adapt_fetch.py
"""
import http.server
import os
import re
import select
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


GIT = shutil.which("git")
if not GIT or not shutil.which("bash"):
    print("skip  git or bash is not available")
    sys.exit(77)

tmp = tempfile.mkdtemp(prefix="adcp_adapt_fetch_")
home = os.path.join(tmp, "home")
os.makedirs(home)

# Recorders: each one, if it is ever run, leaves a line in its own file.
marks = {n: os.path.join(tmp, n + ".used") for n in ("helper", "askpass", "sshaskpass")}
for name, mark in marks.items():
    path = os.path.join(tmp, name)
    with open(path, "w") as fh:
        fh.write('#!/bin/sh\necho "$0 $*" >> "%s"\necho username=someone\necho password=secret\n'
                 % mark)
    os.chmod(path, 0o755)

# A git configuration of the test's own, with a credential helper in it, as a
# person who saved a GitHub token has; both the general entry and one for the
# stand-in server's address.
gitconfig = os.path.join(tmp, "gitconfig")

# Nothing from the environment this runs in may point at a real installation.
base_env = {k: v for k, v in os.environ.items()
            if k not in ("ADFRSUITE_HOME", "ADCP_HOME", "ADAPT_HOME", "ADCPSUITE_HOME",
                         "ADAPT_REPO_URL", "ADAPT_REPO_BRANCH", "GIT_TERMINAL_PROMPT",
                         "GIT_ASKPASS", "SSH_ASKPASS", "ADCP_PYTHON", "ADCP_ENV",
                         "PYTHONPATH", "GIT_DIR", "GIT_WORK_TREE")}
# Nor may PATH: update.sh puts ADapt beside the ADFRsuite whose programs are
# on it when activate.sh names none.
base_env["PATH"] = os.pathsep.join(
    d for d in base_env.get("PATH", "").split(os.pathsep)
    if d and not os.path.exists(os.path.join(d, "prepare_receptor")))
base_env.update(HOME=home, GIT_CONFIG_GLOBAL=gitconfig, GIT_CONFIG_NOSYSTEM="1",
                GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@example.org",
                GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@example.org",
                GIT_ASKPASS=os.path.join(tmp, "askpass"),
                SSH_ASKPASS=os.path.join(tmp, "sshaskpass"))


def git(cwd, *args):
    return subprocess.run([GIT, *args], cwd=cwd, env=base_env, capture_output=True,
                          text=True, check=True).stdout.strip()


# --- the stand-in for a private repository --------------------------------
class SignInFirst(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        self.send_response(401)
        self.send_header("WWW-Authenticate", 'Basic realm="GitHub"')
        self.send_header("Content-Length", "0")
        self.end_headers()

    do_POST = do_GET

    def log_message(self, *a):
        pass


server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), SignInFirst)
threading.Thread(target=server.serve_forever, daemon=True).start()
PRIVATE = "http://127.0.0.1:%d/michelsanner/adapt.git" % server.server_address[1]
UNREACHABLE = "http://127.0.0.1:1/michelsanner/adapt.git"

with open(gitconfig, "w") as fh:
    fh.write('[init]\n\tdefaultBranch = master\n'
             '[credential]\n\thelper = %s\n'
             '[credential "http://127.0.0.1:%d"]\n\thelper = %s\n'
             % (os.path.join(tmp, "helper"), server.server_address[1],
                os.path.join(tmp, "helper")))

# --- the stand-in for GitHub ----------------------------------------------
upstream = os.path.join(tmp, "github", "adapt.git")
os.makedirs(os.path.dirname(upstream))
git(tmp, "init", "-q", "--bare", "-b", "master", upstream)
author = os.path.join(tmp, "author")
git(tmp, "clone", "-q", upstream, author)
os.makedirs(os.path.join(author, "bin"))
os.makedirs(os.path.join(author, "adapt"))


def publish(message, embeddable=True, extra="", force=False):
    """A new commit on master upstream; its short id."""
    with open(os.path.join(author, "bin", "adapt"), "w") as fh:
        fh.write("#!/bin/sh\nexec python -m adapt \"$@\"\n")
    with open(os.path.join(author, "adapt", "__main__.py"), "w") as fh:
        fh.write("import argparse\nap = argparse.ArgumentParser()\n"
                 + ('ap.add_argument("--base-path")\n' if embeddable else "")
                 + "# %s\n%s" % (message, extra))
    git(author, "add", "-A")
    git(author, "commit", "-q", "-m", message)
    git(author, "push", "-q", *(["--force"] if force else []), "origin", "HEAD:master")
    return git(author, "rev-parse", "--short", "HEAD")


c1 = publish("one")

# --- a checkout of ADCP_docking, installed beside a stand-in ADFRsuite -----
fakepy = os.path.join(tmp, "python")
with open(fakepy, "w") as fh:
    fh.write("#!/bin/sh\nexit 0\n")
os.chmod(fakepy, 0o755)


def deployment(name):
    """A checkout with this update.sh and its helper, cloned from a remote of
    its own so its pull works, and an activate.sh naming a tools folder."""
    remote = os.path.join(tmp, name + "_origin.git")
    git(tmp, "init", "-q", "--bare", "-b", "main", remote)
    src = os.path.join(tmp, name + "_src")
    git(tmp, "clone", "-q", remote, src)
    os.makedirs(os.path.join(src, "scripts"))
    shutil.copy(os.path.join(ROOT, "update.sh"), src)
    shutil.copy(os.path.join(ROOT, "install.sh"), src)
    shutil.copy(os.path.join(ROOT, "scripts", "adapt_fetch.sh"), os.path.join(src, "scripts"))
    with open(os.path.join(src, ".gitignore"), "w") as fh:
        fh.write("activate.sh\ntools/\n")
    git(src, "add", "-A")
    git(src, "commit", "-q", "-m", "deploy")
    git(src, "push", "-q", "origin", "HEAD:main")
    dep = os.path.join(tmp, name)
    git(tmp, "clone", "-q", remote, dep)
    tools = os.path.join(tmp, name + "_tools", "docking")
    adfr = os.path.join(tools, "ADFRsuite")
    os.makedirs(os.path.join(adfr, "bin"))
    for b in ("adcp", "agfr", "prepare_receptor"):
        p = os.path.join(adfr, "bin", b)
        open(p, "w").write("#!/bin/sh\nexit 0\n")
        os.chmod(p, 0o755)
    with open(os.path.join(dep, "activate.sh"), "w") as fh:
        fh.write('export ADCP_HOME="%s"\nexport ADCP_PYTHON="%s"\n'
                 'export ADFRSUITE_HOME="%s"\n' % (dep, fakepy, adfr))
    return dep, tools


def update(dep, url, extra=None):
    env = dict(base_env, ADAPT_REPO_URL=url, **(extra or {}))
    p = subprocess.run(["bash", os.path.join(dep, "update.sh")], cwd=dep, env=env,
                       stdin=subprocess.DEVNULL, capture_output=True, text=True,
                       timeout=300)
    return p.returncode, p.stdout + p.stderr


def adapt_line(out):
    plain = re.sub(r"\x1b\[[0-9;]*m", "", out)
    return next((l.strip() for l in plain.splitlines() if l.strip().startswith("ADapt:")),
                "(no ADapt line)\n" + plain[-900:])


def head(path):
    return git(path, "rev-parse", "--short", "HEAD")


def leftovers(tools):
    return [n for n in os.listdir(tools) if n.startswith(".adapt.fetching")] \
        if os.path.isdir(tools) else []


dep, tools = deployment("deploy")
adapt_dir = os.path.join(tools, "adapt")

# --- private, then unreachable: one calm line, nothing created ------------
rc, out = update(dep, PRIVATE)
check(rc == 0, "update.sh failed while ADapt is private (exit %d):\n%s" % (rc, out[-900:]))
check("not public yet" in adapt_line(out) and "works without it" in adapt_line(out),
      "a private ADapt is not reported calmly: %s" % adapt_line(out))
check(not os.path.exists(adapt_dir) and not leftovers(tools),
      "a failed fetch left a folder behind: %s" % os.listdir(tools))
check("Username for" not in out and "Password for" not in out,
      "update.sh asked for a username while ADapt is private:\n" + out[-600:])
rc, out = update(dep, UNREACHABLE)
check(rc == 0 and "could not be reached" in adapt_line(out),
      "an unreachable ADapt is not reported calmly (exit %d): %s" % (rc, adapt_line(out)))
check(not os.path.exists(adapt_dir), "an unreachable ADapt left a folder behind")

# --- public: cloned where the interface finds it --------------------------
rc, out = update(dep, upstream)
check(rc == 0 and os.path.isfile(os.path.join(adapt_dir, "bin", "adapt")),
      "update.sh did not clone ADapt when it is public (exit %d):\n%s" % (rc, out[-900:]))
check(os.path.isdir(adapt_dir) and head(adapt_dir) == c1
      and ("fetched at %s" % c1) in adapt_line(out),
      "the fresh clone is not at the upstream commit, or not said: %s" % adapt_line(out))
check(not leftovers(tools), "the clone left its temporary folder behind")
found = subprocess.run(
    [sys.executable, "-c", "import sys; sys.path.insert(0, sys.argv[1]); "
     "from adcp_dock import adapt; print(adapt.find() or '')", ROOT],
    env=dict(base_env, ADFRSUITE_HOME=os.path.join(tools, "ADFRsuite"), ADCP_HOME=dep),
    capture_output=True, text=True).stdout.strip()
check(os.path.realpath(found) == os.path.realpath(adapt_dir),
      "the interface does not find the ADapt update.sh fetched (it finds %r)" % found)

# --- nothing new upstream ----------------------------------------------------
rc, out = update(dep, upstream)
check(rc == 0 and "the newest on master" in adapt_line(out) and head(adapt_dir) == c1,
      "an up to date ADapt is not reported as such: %s" % adapt_line(out))

# --- private again, with a checkout in place: kept, and said why -----------
rc, out = update(dep, PRIVATE)
check(rc == 0 and head(adapt_dir) == c1 and "kept as it is" in adapt_line(out)
      and "not public yet" in adapt_line(out),
      "an installed ADapt is not kept as it is while private: %s" % adapt_line(out))

# --- fast forward ------------------------------------------------------------
c2 = publish("two")
rc, out = update(dep, upstream)
check(rc == 0 and head(adapt_dir) == c2 and ("%s -> %s" % (c1, c2)) in adapt_line(out),
      "update.sh did not move ADapt forward from %s to %s: %s" % (c1, c2, adapt_line(out)))

# --- the newest commit is not an ADapt this interface can embed: undone -----
c3 = publish("three, which cannot be embedded", embeddable=False)
rc, out = update(dep, upstream)
check(rc == 0 and head(adapt_dir) == c2,
      "update.sh left ADapt on %s, which cannot be embedded (it is at %s)"
      % (c3, head(adapt_dir)))
check(("kept at %s" % c2) in adapt_line(out) and "undone" in adapt_line(out),
      "the rollback is not said: %s" % adapt_line(out))
check(git(adapt_dir, "status", "--porcelain") == "",
      "the rollback left the checkout with changes")
c4 = publish("four, embeddable again")
rc, out = update(dep, upstream)
check(head(adapt_dir) == c4, "after a rollback, a later good commit was not taken (at %s, "
      "upstream %s): %s" % (head(adapt_dir), c4, adapt_line(out)))

# --- local changes are never touched -----------------------------------------
main_py = os.path.join(adapt_dir, "adapt", "__main__.py")
with open(main_py, "a") as fh:
    fh.write("# a local edit\n")
c5 = publish("five")
rc, out = update(dep, upstream)
check(rc == 0 and head(adapt_dir) == c4 and "a local edit" in open(main_py).read(),
      "update.sh moved or overwrote an ADapt with local changes (at %s)" % head(adapt_dir))
check("local changes" in adapt_line(out),
      "a checkout with local changes was left without saying why: %s" % adapt_line(out))
git(adapt_dir, "checkout", "-q", "--", ".")

# --- another branch is left alone ---------------------------------------------
git(adapt_dir, "checkout", "-q", "-b", "mine")
rc, out = update(dep, upstream)
check(rc == 0 and head(adapt_dir) == c4 and "on branch mine" in adapt_line(out),
      "a checkout on another branch was not left alone: %s" % adapt_line(out))
git(adapt_dir, "checkout", "-q", "master")
rc, out = update(dep, upstream)
check(head(adapt_dir) == c5, "back on master, ADapt was not updated: %s" % adapt_line(out))

# --- the checkout's own remote is never used ----------------------------------
# As where a checkout came from a tarball and its origin is an ssh address
# that logs in with another repository's key.
git(adapt_dir, "remote", "set-url", "origin", "git@github.invalid:michelsanner/adapt.git")
c6 = publish("six")
rc, out = update(dep, upstream)
check(head(adapt_dir) == c6, "with an unusable origin, ADapt was not updated from the "
      "https address: %s" % adapt_line(out))

# --- a folder that is not a git checkout is left alone ------------------------
dep2, tools2 = deployment("plain")
plain = os.path.join(tools2, "adapt")
shutil.copytree(adapt_dir, plain, ignore=shutil.ignore_patterns(".git"))
before = sorted(os.listdir(plain))
rc, out = update(dep2, upstream)
check(rc == 0 and "not a git checkout" in adapt_line(out) and sorted(os.listdir(plain)) == before
      and not os.path.exists(os.path.join(plain, ".git")),
      "a plain ADapt folder was not left alone: %s" % adapt_line(out))

# ... nor one inside another repository, whose own history must not be touched:
# no ADFRsuite named, so ADapt lives in the checkout's own tools folder.
dep3, _ = deployment("inside")
with open(os.path.join(dep3, "activate.sh"), "w") as fh:
    fh.write('export ADCP_HOME="%s"\nexport ADCP_PYTHON="%s"\n' % (dep3, fakepy))
inside = os.path.join(dep3, "tools", "docking", "adapt")
shutil.copytree(plain, inside)
before_head = git(dep3, "rev-parse", "HEAD")
rc, out = update(dep3, upstream)
check(rc == 0 and "not a git checkout" in adapt_line(out)
      and git(dep3, "rev-parse", "HEAD") == before_head
      and git(dep3, "status", "--porcelain") == "",
      "an ADapt folder inside the ADCP_docking checkout was treated as a checkout: %s"
      % adapt_line(out))

# --- no ADFRsuite named in activate.sh: the one on PATH ------------------------
dep5, tools5 = deployment("onpath")
with open(os.path.join(dep5, "activate.sh"), "w") as fh:
    fh.write('export ADCP_HOME="%s"\nexport ADCP_PYTHON="%s"\n' % (dep5, fakepy))
rc, out = update(dep5, upstream, {"PATH": os.path.join(tools5, "ADFRsuite", "bin")
                                  + os.pathsep + base_env["PATH"]})
check(rc == 0 and os.path.isdir(os.path.join(tools5, "adapt", ".git")),
      "with ADFRsuite only on PATH, ADapt was not fetched beside it: %s" % adapt_line(out))

# --- ADCP 1.1, which ADapt runs on: found where the interface finds it -------
suite = os.path.join(tools, "adcpsuite")
os.makedirs(os.path.join(suite, "bin"))
os.makedirs(os.path.join(suite, "lib", "python3.7", "site-packages", "ADCP"))
open(os.path.join(suite, "bin", "adcp"), "w").close()
open(os.path.join(suite, "lib", "python3.7", "site-packages", "ADCP", "runADCP.py"), "w").close()
rc, out = update(dep, upstream)
check(("ADCP 1.1: %s" % suite) in out, "update.sh does not report the ADCP 1.1 it has:\n"
      + out[-600:])
exe = subprocess.run(
    [sys.executable, "-c", "import sys; sys.path.insert(0, sys.argv[1]); "
     "from adcp_dock import adcp11; print(adcp11.find() or '')", ROOT],
    env=dict(base_env, ADFRSUITE_HOME=os.path.join(tools, "ADFRsuite"), ADCP_HOME=dep),
    capture_output=True, text=True).stdout.strip()
check(os.path.realpath(exe) == os.path.realpath(os.path.join(suite, "bin", "adcp")),
      "update.sh and the interface disagree on ADCP 1.1 (%r)" % exe)
shutil.rmtree(suite)
rc, out = update(dep, upstream)
check(rc == 0 and "ADCP 1.1 is not installed" in out
      and ("install-adcp11 --dest %s" % suite) in out,
      "a missing ADCP 1.1 is not reported with the command that adds it:\n" + out[-600:])


# --- never a prompt, even with a terminal to prompt on --------------------------
def on_terminal(argv, env, seconds=120):
    """Runs argv with a terminal of its own. (output, finished): finished is
    False when it was still running after `seconds`, which is what a prompt
    waiting for an answer looks like."""
    import pty
    pid, fd = pty.fork()
    if pid == 0:
        try:
            os.execvpe(argv[0], argv, env)
        finally:
            os._exit(127)
    out, end = b"", time.time() + seconds
    while time.time() < end:
        r, _, _ = select.select([fd], [], [], 0.5)
        if r:
            try:
                chunk = os.read(fd, 4096)
            except OSError:
                chunk = b""
            if not chunk:
                break
            out += chunk
        elif os.waitpid(pid, os.WNOHANG)[0]:
            break
    else:
        os.kill(pid, signal.SIGKILL)
        os.waitpid(pid, 0)
        os.close(fd)
        return out.decode(errors="replace"), False
    try:
        os.waitpid(pid, 0)
    except ChildProcessError:
        pass
    os.close(fd)
    return out.decode(errors="replace"), True


try:
    import pty  # noqa: F401
    # Without the anonymous settings, git on a terminal does ask: proof the
    # harness would see a prompt. The askpass programs are removed so that it
    # asks on the terminal itself.
    plain_env = {k: v for k, v in base_env.items() if k not in ("GIT_ASKPASS", "SSH_ASKPASS")}
    plain_env["GIT_CONFIG_GLOBAL"] = os.devnull
    said, finished = on_terminal([GIT, "ls-remote", PRIVATE], plain_env, seconds=15)
    can_tell = (not finished) or "Username for" in said
except (ImportError, OSError):
    can_tell = False
if can_tell:
    shutil.rmtree(adapt_dir)
    said, finished = on_terminal(["bash", os.path.join(dep, "update.sh")],
                                 dict(base_env, ADAPT_REPO_URL=PRIVATE))
    check(finished and "Username for" not in said and "Password for" not in said,
          "update.sh waited at a prompt for ADapt on a terminal:\n" + said[-600:])
    check("not public yet" in said, "on a terminal, a private ADapt is not reported:\n"
          + said[-600:])
    git(tmp, "clone", "-q", upstream, adapt_dir)
    said, finished = on_terminal(["bash", os.path.join(dep, "update.sh")],
                                 dict(base_env, ADAPT_REPO_URL=PRIVATE))
    check(finished and "Username for" not in said,
          "update.sh waited at a prompt updating ADapt on a terminal:\n" + said[-600:])
else:
    print("note  no terminal could be made here; the prompt checks ran without one")

# --- install.sh ----------------------------------------------------------------
# A stand-in environment: its python answers the import checks and pip, and runs
# anything else with this Python, so the install's own steps run for real but
# install nothing into a real environment.
envdir = os.path.join(tmp, "env")
os.makedirs(os.path.join(envdir, "bin"))
with open(os.path.join(envdir, "bin", "python"), "w") as fh:
    fh.write('#!/bin/sh\ncase "$1 $2" in\n  "-c import "*) exit 0 ;;\n  "-m pip") exit 0 ;;\n'
             'esac\nexec "%s" "$@"\n' % sys.executable)
os.chmod(os.path.join(envdir, "bin", "python"), 0o755)
for b in ("python3", "tleap", "sander", "MMPBSA.py", "ante-MMPBSA.py", "pdb2pqr"):
    p = os.path.join(envdir, "bin", b)
    if b == "python3":
        os.symlink("python", p)
    else:
        open(p, "w").write("#!/bin/sh\nexit 0\n")
        os.chmod(p, 0o755)


def install(dep, tools, url, *args):
    env = dict(base_env, ADAPT_REPO_URL=url, ADCP_ENV=envdir,
               ADFRSUITE_HOME=os.path.join(tools, "ADFRsuite"), ADCP_SHARED=tmp)
    p = subprocess.run(["bash", os.path.join(dep, "install.sh"), *args], cwd=dep, env=env,
                       stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=300)
    # install.sh wraps its long lines; checked here as one line of text.
    return p.returncode, " ".join((p.stdout + p.stderr).split())


dep4, tools4 = deployment("fresh")
os.remove(os.path.join(dep4, "activate.sh"))
os.symlink(os.path.join(ROOT, "adcp_dock"), os.path.join(dep4, "adcp_dock"))
adapt4 = os.path.join(tools4, "adapt")

rc, out = install(dep4, tools4, upstream, "--check")
check(rc == 0 and ("ADapt, into %s, from %s" % (adapt4, upstream)) in out,
      "install.sh --check does not say it would fetch a public ADapt (exit %d):\n%s"
      % (rc, out[-900:]))
check(("ADCP 1.1, into %s" % os.path.join(tools4, "adcpsuite")) in out and "--no-adcp11" in out,
      "install.sh --check does not say it would install ADCP 1.1:\n" + out[-900:])
check(not os.path.exists(adapt4) and not os.path.exists(os.path.join(dep4, "activate.sh")),
      "install.sh --check changed something")
rc, out = install(dep4, tools4, PRIVATE, "--check", "--no-adcp11")
check(rc == 0 and "would not be fetched today" in out and "not public yet" in out,
      "install.sh --check does not say a private ADapt would be skipped:\n" + out[-900:])
check("--no-adcp11 leaves it out" in out and ("ADCP 1.1, into") not in out,
      "install.sh --check --no-adcp11 still says it would install ADCP 1.1:\n" + out[-900:])

# The install itself, with ADCP 1.1 already there so nothing is downloaded.
suite4 = os.path.join(tools4, "adcpsuite")
os.makedirs(os.path.join(suite4, "bin"))
os.makedirs(os.path.join(suite4, "lib", "python3.7", "site-packages", "ADCP"))
open(os.path.join(suite4, "bin", "adcp"), "w").close()
open(os.path.join(suite4, "lib", "python3.7", "site-packages", "ADCP", "runADCP.py"), "w").close()
rc, out = install(dep4, tools4, PRIVATE)
check(rc == 0 and "not public yet" in out and "works without it" in out,
      "install.sh did not carry on calmly while ADapt is private (exit %d):\n%s"
      % (rc, out[-1200:]))
check(not os.path.exists(adapt4) and os.path.isfile(os.path.join(dep4, "activate.sh")),
      "install.sh while ADapt is private: a folder was left, or the install did not finish")
check(("USING ADCP 1.1 at %s" % suite4) in out,
      "install.sh does not find the ADCP 1.1 beside ADFRsuite:\n" + out[-900:])
rc, out = install(dep4, tools4, upstream)
check(rc == 0 and os.path.isdir(adapt4) and head(adapt4) == c6 and "fetched at" in out,
      "install.sh did not fetch a public ADapt (exit %d):\n%s" % (rc, out[-1200:]))
publish("seven")
rc, out = install(dep4, tools4, upstream)
check(rc == 0 and head(adapt4) == c6 and "./update.sh keeps it up to date" in out,
      "install.sh did not leave an existing ADapt to update.sh (exit %d):\n%s"
      % (rc, out[-900:]))

# --- the published history rewritten, as can happen when a repository is ----
# made public: the checkout is left alone, and once moved aside a fresh one
# is fetched in its place.
was = head(adapt_dir)
git(author, "checkout", "-q", "--orphan", "published")
c_pub = publish("published afresh", force=True)
rc, out = update(dep, upstream)
check(rc == 0 and head(adapt_dir) == was and "separate ways" in adapt_line(out)
      and "move it aside" in adapt_line(out),
      "a checkout whose history was rewritten upstream was not left alone, with what "
      "to do: %s" % adapt_line(out))
os.rename(adapt_dir, adapt_dir + ".old")
rc, out = update(dep, upstream)
check(rc == 0 and os.path.isdir(adapt_dir) and head(adapt_dir) == c_pub,
      "after moving the old checkout aside, the rewritten history was not fetched: %s"
      % adapt_line(out))

# --- no helper, no askpass program, ever ------------------------------------
for name, mark in marks.items():
    check(not os.path.exists(mark), "the %s was run: %s" % (
        name, open(mark).read().strip() if os.path.exists(mark) else ""))
# and that the recorder would have caught it: plain git uses the helper
subprocess.run([GIT, "ls-remote", PRIVATE], env=dict(base_env, GIT_TERMINAL_PROMPT="0"),
               stdin=subprocess.DEVNULL, capture_output=True, timeout=60)
check(os.path.exists(marks["helper"]) or os.path.exists(marks["askpass"]),
      "the recorders did not notice plain git asking for a credential, so this test "
      "proves nothing about the scripts")

server.shutdown()
shutil.rmtree(tmp, ignore_errors=True)

if fails:
    print("test_adapt_fetch FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  ADapt is cloned when public, skipped quietly while private, moved forward only "
      "when safe and back when it would break, never prompts and never uses a stored "
      "credential, %d checks" % ran)
