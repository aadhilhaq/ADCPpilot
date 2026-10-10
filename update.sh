#!/usr/bin/env bash
# Update this checkout to the latest version, and say what it now has.
#
# By Aadhil Haq
#
#     cd /path/to/ADCP_docking && git pull; ./update.sh
#
# The pull in front is not redundant: it is what brings in a newer copy of
# THIS script. Running it afterwards is what applies whatever that copy does.
# The two are joined with ";", not "&&": when the published history has been
# rewritten, a plain pull refuses, and this script is what moves the checkout
# onto the new history.
#
# It also brings ADapt, Michel Sanner's separate tool, up to date beside
# ADFRsuite, or fetches it there when it is missing and its repository is
# public. That step never stops or fails the update.
#
# What it does NOT do: build an environment. install.sh owns that, it only
# has to run once, and it is the thing that needs a login node and a network.
# So a missing ADCP 1.1 is reported with the command that adds it.

# Run, not sourced, and by bash. These lines are read by any shell.
if [ -z "${BASH_VERSION:-}" ]; then
  echo "Run this with bash:  ./update.sh   (or: bash update.sh)" >&2
  return 2 2>/dev/null || exit 2
fi
if [ "${BASH_SOURCE[0]}" != "$0" ]; then
  echo "Run this script, do not source it:  ./update.sh" >&2
  return 2
fi
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
# Never prompt: an https remote with a cached credential pulls normally, and
# without one the pull fails at once instead of waiting for a password. The
# checkout's own ssh command is kept when it has one (core.sshCommand): a
# cluster whose ~/.ssh/config gives github.com another repository's deploy
# key pins this repository's key that way, and a plain "ssh" would offer the
# other key and be told the repository does not exist.
SSH_CMD=$(git config --get core.sshCommand 2>/dev/null || true)
OUT=$(GIT_TERMINAL_PROMPT=0 GIT_SSH_COMMAND="${SSH_CMD:-ssh} -o BatchMode=yes" \
      git pull --ff-only 2>&1); RC=$?
REWRITTEN=""
if [ "$RC" -ne 0 ] && printf '%s' "$OUT" | grep -qi \
     "diverg\|unrelated histories\|not possible to fast-forward"; then
  # The fetch worked and the history was rewritten upstream, which --ff-only
  # cannot cross. This checkout is a deployment rather than a place to
  # commit, so it takes what the remote has, unless tracked files have local
  # edits, which are never discarded without asking.
  if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    warn "the upstream history was rewritten, and this checkout has"
    warn "local edits to tracked files. Not discarding them. Save them,"
    warn "then:  git -C $HERE fetch origin && git -C $HERE reset --hard origin/main"
    git status --short --untracked-files=no | head -5
  else
    note "upstream history was rewritten; resetting to origin/main"
    git reset -q --hard origin/main && RC=0 && REWRITTEN=1
  fi
elif [ "$RC" -ne 0 ] && [ "${URL#https://}" != "$URL" ]; then
  # Checked only after the rewrite: an https remote whose fetch worked but
  # could not fast-forward is not a sign-in problem.
  warn "could not pull from $URL."
  warn "Run this on a login node: compute nodes usually have no internet."
  warn "If it fails there too and the message below asks for a username,"
  warn "the repository needs a sign in, which this never waits for. A cluster"
  warn "with access pulls over ssh instead, with a read-only deploy key:"
  note "  git -C $HERE remote set-url origin \
$(printf '%s' "$URL" | sed -e 's|^https://[^/@]*@|https://|' -e 's|^https://\([^/]*\)/|git@\1:|')"
fi
AFTER=$(git rev-parse --short HEAD 2>/dev/null)
if [ "$RC" -ne 0 ]; then
  warn "NOT UPDATED (git exited $RC), still at $AFTER"
  note "$(printf '%s' "$OUT" | head -3)"
elif [ -n "$REWRITTEN" ]; then
  note "$BEFORE -> $AFTER (the history was rewritten upstream)"
elif [ "$BEFORE" = "$AFTER" ]; then
  note "already at $AFTER"
else
  note "$BEFORE -> $AFTER ($(git rev-list --count "$BEFORE..$AFTER") new commit(s))"
fi

# --- reinstall the package, not the environment -----------------------------
# Only where install.sh installed it: an environment it built, which imports
# adcp_dock from this checkout. That is a few files with --no-deps, and
# nothing else in the environment moves. An environment it did not build is
# never written to; activate.sh puts this checkout first on PYTHONPATH there,
# so the pull is already what runs. Never the system python3, whose pip would
# fall back to installing into the home directory.
say "the adcp-dock command"
PY=""
# shellcheck disable=SC1091
[ -f "$HERE/activate.sh" ] && . "$HERE/activate.sh" >/dev/null 2>&1
for c in "${ADCP_PYTHON:-}" "$HERE/../envs/adcp/bin/python"; do
  [ -n "${c:-}" ] && [ -x "${c:-}" ] && { PY="$c"; break; }
done
if [ -z "$PY" ]; then
  warn "no environment is recorded here (activate.sh is missing or names none)."
  warn "Run ./install.sh <folder> once, on a login node; INSTALL.md says how."
else
  # Asked from /, so that the checkout itself is not what answers.
  AT=$(cd / && env -u PYTHONPATH -u PYTHONHOME "$PY" -c \
       "import adcp_dock, os; print(os.path.dirname(os.path.dirname(os.path.abspath(adcp_dock.__file__))))" \
       2>/dev/null)
  if [ -n "$AT" ] && [ "$AT" -ef "$HERE" ]; then
    if env -u PYTHONPATH "$PY" -m pip install --quiet --no-deps --no-cache-dir \
         --disable-pip-version-check -e "$HERE" 2>&1 | sed 's/^/  /'; then
      note "reinstalled into $(dirname "$(dirname "$PY")")"
    else
      warn "could not reinstall. On a login node:  $PY -m pip install --no-deps -e $HERE"
    fi
  else
    note "runs from this checkout with $PY (nothing to reinstall)"
  fi
fi

# --- clear stale bytecode ---------------------------------------------------
# A .pyc whose source moved is read in preference to the source on some
# filesystems, and the symptom is an update that changed nothing.
say "clearing cached bytecode"
find "$HERE" -name '__pycache__' -type d -prune -exec rm -rf {} + 2>/dev/null
note "done"

# --- ADapt, and the ADCP 1.1 it runs on ---------------------------------------
# ADapt is Michel Sanner's own repository, kept beside ADFRsuite, where the
# interface finds it. It is fetched anonymously from its https address, never
# through the remote its checkout has configured, and only ever moved forward
# along its branch; scripts/adapt_fetch.sh says what is left alone and why.
# Nothing here can stop this script, prompt, or make it fail: ADCP docking
# works without ADapt.
say "ADapt, and the ADCP 1.1 it runs on"
if [ -f "$HERE/scripts/adapt_fetch.sh" ]; then
  # shellcheck source=scripts/adapt_fetch.sh
  . "$HERE/scripts/adapt_fetch.sh"
  # The ADFRsuite this install uses: the one activate.sh names, else the one
  # whose programs are on PATH, as the interface would find them.
  ADFR_AT="${ADFRSUITE_HOME:-}"
  if [ -z "$ADFR_AT" ] && command -v prepare_receptor >/dev/null 2>&1; then
    ADFR_AT=$(dirname "$(dirname "$(command -v prepare_receptor)")")
  fi
  ADAPT_DIR=$(adapt_dest "$ADFR_AT" "${ADCP_HOME:-}")
  if [ -z "$ADAPT_DIR" ]; then
    note "ADapt: not fetched, since no ADFRsuite was found to put it beside (activate.sh names none, and none is on PATH; run ./install.sh)"
  else
    if [ -e "$ADAPT_DIR" ]; then
      adapt_update "$ADAPT_DIR"
    else
      adapt_clone "$ADAPT_DIR"
    fi
    case "$ADAPT_STATE" in
      changed|branch|detached|ahead|diverged|blocked|rolledback|broken|unreadable|fetched_old)
        warn "ADapt: $ADAPT_LINE" ;;
      *) note "ADapt: $ADAPT_LINE" ;;
    esac
    USED=$(adapt_found "$HERE")
    if [ -n "$USED" ] && [ "$(cd "$USED" && pwd -P)" != "$(cd "$ADAPT_DIR" 2>/dev/null && pwd -P)" ]; then
      note "  the interface uses the ADapt in $USED instead, which this leaves alone"
    fi
  fi
  SUITE=$(adcp11_found)
  if [ -n "$SUITE" ]; then
    note "ADCP 1.1: $SUITE"
  elif [ -n "$ADAPT_DIR" ]; then
    warn "ADCP 1.1 is not installed, and ADapt runs on it. To add it, on a login node (ten minutes, 3 GB):"
    note "  ${PY:-python3} -m adcp_dock.cli install-adcp11 --dest $(dirname "$ADAPT_DIR")/adcpsuite"
  else
    warn "ADCP 1.1 is not installed, and ADapt runs on it. ./install.sh adds both."
  fi
else
  note "scripts/adapt_fetch.sh is missing from this checkout; ADapt left as it is"
fi

# --- what this install has now ----------------------------------------------
# The checkout is installed in editable mode, so it is what runs. A copy left
# in site-packages by an older, non-editable install would shadow it, and a
# pull would then change files that nothing reads.
if [ -n "$PY" ]; then
  say "what this install has"
  WHERE=$("$PY" -c "import adcp_dock,os;print(os.path.dirname(os.path.abspath(adcp_dock.__file__)))" 2>/dev/null)
  if [ "$WHERE" = "$HERE/adcp_dock" ]; then
    note "adcp_dock $("$PY" -c 'import adcp_dock;print(adcp_dock.__version__)' 2>/dev/null), from this checkout"
  else
    warn "python imports adcp_dock from ${WHERE:-nowhere}, not from this checkout."
    warn "Run ./install.sh again (no folder needed), which records this checkout in activate.sh."
  fi
  "$PY" -m adcp_dock.cli check 2>&1 | sed 's/^/  /'
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
