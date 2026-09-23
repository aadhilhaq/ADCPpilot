#!/usr/bin/env bash
# Update this checkout to the latest version, and say what it now has.
#
#     cd /path/to/ADCP_docking && git pull && ./update.sh
#
# The pull in front is not redundant: it is what brings in a newer copy of
# THIS script. Running it afterwards is what applies whatever that copy does.
#
# What it does NOT do: build an environment. install.sh owns that, it only
# has to run once, and it is the thing that needs a login node and a network.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

say()  { printf '\n\033[1m[update]\033[0m %s\n' "$1"; }
note() { printf '  %s\n' "$1"; }
warn() { printf '  \033[33m%s\033[0m\n' "$1"; }

[ -d .git ] || { warn "$HERE is not a git checkout. Use install.sh."; exit 1; }

# --- pull -------------------------------------------------------------------
# Never prompt. An https remote would ask for a GitHub username and password,
# which looks like a hang and fails anyway, since GitHub no longer accepts
# passwords. A cluster pulls with a read-only deploy key.
say "pulling"
URL=$(git remote get-url origin 2>/dev/null)
BEFORE=$(git rev-parse --short HEAD 2>/dev/null)
# Try it first: an https remote with a cached credential pulls normally. The
# prompt is disabled, so the worst outcome is a clean failure, which is when
# the deploy-key advice is printed.
OUT=$(GIT_TERMINAL_PROMPT=0 GIT_SSH_COMMAND="ssh -o BatchMode=yes" \
      git pull --ff-only 2>&1); RC=$?
case "$RC:$URL" in
  [1-9]*:https://*)
    warn "could not pull, and the remote is https so it cannot ask for a"
    warn "password here. This cluster pulls with a deploy key:"
    note "  git -C $HERE remote set-url origin \
$(printf '%s' "$URL" | sed -e 's|^https://[^/@]*@|https://|' -e 's|^https://\([^/]*\)/|git@\1:|')"
    warn 'If that says "Permission denied (publickey)", this repository has'
    warn "no deploy key yet: add this machine's public key to it on GitHub,"
    warn "read-only."
    note "$(printf '%s' "$OUT" | head -2)" ;;
  *)
    AFTER=$(git rev-parse --short HEAD 2>/dev/null)
    # A failed pull leaves the commit where it was, as a pull with nothing to
    # fetch does, so the exit status is what tells the two apart.
    if [ "$RC" -ne 0 ] && printf '%s' "$OUT" | grep -qi \
         "diverge\|unrelated histories\|not possible to fast-forward"; then
      # The history was rewritten upstream. --ff-only cannot cross that, and
      # this checkout is a deployment rather than a place to commit, so it
      # takes what the remote has. It refuses if tracked files have local
      # edits, so local changes are never discarded without asking.
      if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
        warn "the upstream history was rewritten, and this checkout has"
        warn "local edits to tracked files. Not discarding them. Save them,"
        warn "then:  git -C $HERE fetch origin && git -C $HERE reset --hard origin/main"
        git status --short --untracked-files=no | head -5
      else
        note "upstream history was rewritten; resetting to origin/main"
        GIT_TERMINAL_PROMPT=0 GIT_SSH_COMMAND="ssh -o BatchMode=yes" \
          git fetch -q origin && git reset -q --hard origin/main && RC=0
        AFTER=$(git rev-parse --short HEAD 2>/dev/null)
        [ "$RC" -eq 0 ] && note "now at $AFTER"
      fi
    fi
    if [ "$RC" -ne 0 ]; then
      warn "NOT UPDATED (git exited $RC), still at $AFTER"
      note "$(printf '%s' "$OUT" | head -3)"
    elif [ "$BEFORE" = "$AFTER" ]; then
      note "already at $AFTER"
    else
      note "$BEFORE -> $AFTER ($(git rev-list --count "$BEFORE..$AFTER") new commit(s))"
    fi ;;
esac

# --- reinstall the package, not the environment -----------------------------
# The code is installed into a shared environment with --no-deps, so this is
# a few files and nothing else in that environment moves. Without it, a pull
# updates the checkout while every job keeps importing the installed copy.
say "reinstalling adcp-dock into the environment it was installed in"
PY=""
# shellcheck disable=SC1091
[ -f "$HERE/activate.sh" ] && . "$HERE/activate.sh" >/dev/null 2>&1
for c in "${ADCP_PYTHON:-}" "$HERE/../envs/hope/bin/python" \
         "$(command -v python3 2>/dev/null)"; do
  [ -n "${c:-}" ] && [ -x "${c:-}" ] && { PY="$c"; break; }
done
if [ -z "$PY" ]; then
  warn "no python found. Run ./install.sh first."
elif "$PY" -m pip install --quiet --no-deps -e "$HERE" 2>&1 | sed 's/^/  /'; then
  note "done"
else
  warn "could not reinstall. On a login node:  $PY -m pip install --no-deps -e $HERE"
fi

# --- clear stale bytecode ---------------------------------------------------
# A .pyc whose source moved is read in preference to the source on some
# filesystems, and the symptom is an update that changed nothing.
say "clearing cached bytecode"
find "$HERE" -name '__pycache__' -type d -prune -exec rm -rf {} + 2>/dev/null
note "done"

# --- say what this install actually has now ---------------------------------
if [ -x "$HERE/check.sh" ]; then
  say "what this install has"
  "$HERE/check.sh" "$HERE"
fi

# A running server keeps the code it started with: page.py is imported once,
# so the old stylesheet and app.js stay until it is restarted, which looks
# like a failed update. Running servers on this node are therefore reported.
RUNNING=$(pgrep -af "adcp[-_]dock.*serve" 2>/dev/null | head -3)
# Only this node: each user starts their own server, and ssh between login
# nodes cannot be assumed. The page itself says when it is older than the
# install.
if [ -n "$RUNNING" ]; then
  printf '\n\033[1;33m[update] A SERVER IS STILL RUNNING THE OLD CODE\033[0m\n'
  printf '%s\n' "$RUNNING" | sed 's/^/  /'
  printf '  Nothing you see in the browser will change until it is restarted.\n'
  printf '  Ctrl-C that process, start it again, then hard-reload the page\n'
  printf '  with Ctrl-Shift-R. The footer shows the build it is serving.\n\n'
fi

cat <<MSG
Nothing already running has any of this yet. The server holds the code it
started with, in memory, until it is restarted:

  Ctrl-C the 'adcp-dock serve' process and start it again, then hard-reload
  the browser (Ctrl-Shift-R).

MSG
