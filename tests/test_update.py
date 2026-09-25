"""update.sh against a history rewritten upstream.

By Aadhil Haq

The case that matters is an https remote whose credential works: the fetch
succeeds and only the fast-forward fails. That is a rewrite, not a sign-in
problem, and the checkout must end up on the new history. A checkout with
local edits to tracked files is left alone.
"""
import os
import shutil
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


GIT = shutil.which("git")
if not GIT or not shutil.which("bash"):
    print("skip  git or bash is not available")
    sys.exit(77)

tmp = tempfile.mkdtemp(prefix="adcp_update_")
env = dict(os.environ, GIT_CONFIG_GLOBAL=os.devnull, GIT_CONFIG_NOSYSTEM="1",
           GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@example.org",
           GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@example.org")


def git(cwd, *args):
    return subprocess.run([GIT, *args], cwd=cwd, env=env, capture_output=True,
                          text=True, check=True).stdout.strip()


# a remote, a checkout of it, then the remote's history rewritten
bare = os.path.join(tmp, "remote.git")
git(tmp, "init", "-q", "--bare", "-b", "main", bare)
work = os.path.join(tmp, "work")
git(tmp, "clone", "-q", bare, work)
shutil.copy(os.path.join(ROOT, "update.sh"), work)
open(os.path.join(work, "README"), "w").write("one\n")
git(work, "add", ".")
git(work, "commit", "-q", "-m", "one")
git(work, "push", "-q", "origin", "HEAD:main")
deploys = {}
for name in ("clean", "edited"):
    deploys[name] = os.path.join(tmp, name)
    git(tmp, "clone", "-q", bare, deploys[name])
git(work, "checkout", "-q", "--orphan", "rewritten")
open(os.path.join(work, "README"), "w").write("two\n")
git(work, "add", ".")
git(work, "commit", "-q", "-m", "two")
git(work, "push", "-q", "--force", "origin", "rewritten:main")
new_head = git(work, "rev-parse", "--short", "HEAD")

# the remote reports itself as https, as Grace's does; nothing else changes
shim = os.path.join(tmp, "bin")
os.makedirs(shim)
with open(os.path.join(shim, "git"), "w") as fh:
    fh.write('#!/bin/sh\nif [ "$1" = remote ] && [ "$2" = get-url ]; then\n'
             '  echo https://github.com/example/repo.git; exit 0\nfi\n'
             'exec "%s" "$@"\n' % GIT)
fakepy = os.path.join(tmp, "python")
with open(fakepy, "w") as fh:
    fh.write("#!/bin/sh\nexit 0\n")
for f in (os.path.join(shim, "git"), fakepy):
    os.chmod(f, 0o755)
run_env = dict(env, PATH=shim + os.pathsep + env["PATH"], ADCP_PYTHON=fakepy)

open(os.path.join(deploys["edited"], "README"), "a").write("a local edit\n")
out = {}
for name, path in deploys.items():
    p = subprocess.run(["bash", os.path.join(path, "update.sh")], cwd=path, env=run_env,
                       capture_output=True, text=True, timeout=120)
    out[name] = p.stdout + p.stderr

check(git(deploys["clean"], "rev-parse", "--short", "HEAD") == new_head,
      "a clean checkout was not moved to the rewritten history:\n" + out["clean"][-600:])
check("rewritten upstream" in out["clean"] and "deploy key" not in out["clean"],
      "a rewrite over https was reported as a sign-in problem:\n" + out["clean"][-600:])
check(git(deploys["edited"], "rev-parse", "--short", "HEAD") != new_head
      and "a local edit" in open(os.path.join(deploys["edited"], "README")).read(),
      "a checkout with local edits was reset, and the edits discarded")
check("local edits" in out["edited"],
      "a checkout with local edits was left without saying why:\n" + out["edited"][-600:])
shutil.rmtree(tmp, ignore_errors=True)

if fails:
    print("test_update FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  update.sh follows a history rewritten upstream, over https too, %d checks" % ran)
