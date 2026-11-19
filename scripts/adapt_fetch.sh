# ADapt, fetched and kept up to date beside the ADCP 1.1 suite, and the suite
# found where the interface finds it. Sourced by install.sh and update.sh; it
# is not a program of its own.
#
# By Aadhil Haq
#
# ADapt, by Michel Sanner, is its own repository. The interface finds it in a
# folder named adapt beside the ADCP 1.1 suite, among other places
# (adcp_pilot/adapt.py), and these functions put it there and keep it current.
#
# Every git command that reaches the repository is anonymous: it never
# prompts, never runs an askpass program and never uses a stored credential.
# So while the repository is private the attempt fails at once and quietly,
# and once it is public it works for anyone, on any cluster. A credential
# helper is not merely unused but harmful here: after a refused sign-in git
# tells the helper to erase the credential it offered, which would delete a
# person's saved GitHub token. And the remote an existing checkout has
# configured is never used: on a cluster whose ssh config gives github.com
# another repository's deploy key, a git@github.com: address fails.
#
# Nothing here exits, and every function returns 0, so a caller running under
# set -e carries on: ADCP docking works without ADapt, and ADapt never stops
# an install or an update. Each function leaves its answer in ADAPT_LINE, one
# line saying what ADapt is now at and why, and ADAPT_STATE, one word.
#
#     ADAPT_REPO_URL      where ADapt is fetched from (a fork, or a test's
#                         local repository); the default is Michel Sanner's
#     ADAPT_REPO_BRANCH   the branch followed, master by default
#     ADAPT_GIT_SECONDS   how long one fetch may take, 300 by default

ADAPT_URL="${ADAPT_REPO_URL:-https://github.com/sannerlab/adapt.git}"
ADAPT_BRANCH="${ADAPT_REPO_BRANCH:-master}"
ADAPT_LINE=""
ADAPT_STATE=""

# git, anonymous and unable to wait for anyone. The askpass programs are
# removed as well as the terminal prompt, since git asks them first even when
# terminal prompts are off. A machine whose network drops packets rather than
# refusing them would otherwise hold an update for minutes, hence the limit.
adapt_git() {
  local limit=""
  command -v timeout >/dev/null 2>&1 && limit="timeout ${ADAPT_GIT_SECONDS:-300}"
  env -u GIT_ASKPASS -u SSH_ASKPASS GIT_TERMINAL_PROMPT=0 GCM_INTERACTIVE=never \
      GIT_SSH_COMMAND="${GIT_SSH_COMMAND:-ssh} -o BatchMode=yes" \
      $limit git -c credential.helper= -c "credential.$ADAPT_URL.helper=" \
      -c core.askPass= "$@"
}

# The checks adcp_pilot/adapt.py makes: is_adapt() and embeddable().
adapt_is() { [ -f "$1/bin/adapt" ] && [ -f "$1/adapt/__main__.py" ]; }
adapt_embeds() { grep -q -e '--base-path' "$1/adapt/__main__.py" 2>/dev/null; }

# What kind of failure git's output and exit status describe: private (the
# repository asks for a sign-in), offline, slow, owner (the folder belongs to
# another user) or other.
adapt_kind() {  # $1 = git's output, $2 = its exit status
  if [ "${2:-1}" = 124 ]; then printf slow; return 0; fi
  case "$1" in
    *"dubious ownership"*) printf owner ;;
    *"could not read Username"*|*"could not read Password"*|*"Authentication failed"*|\
    *"terminal prompts disabled"*|*"Repository not found"*|*"requested URL returned error: 40"[134]*)
      printf private ;;
    *"Could not resolve"*|*"Failed to connect"*|*"Couldn't connect"*|*"timed out"*|\
    *"Network is unreachable"*|*"Connection refused"*|*"No route to host"*)
      printf offline ;;
    *) printf other ;;
  esac
}

# Why a fetch failed, in a few words.
adapt_why() {  # $1 = git's output, $2 = its exit status
  case "$(adapt_kind "$1" "${2:-1}")" in
    slow)    printf '%s did not answer within %s seconds' "$ADAPT_URL" "${ADAPT_GIT_SECONDS:-300}" ;;
    owner)   printf 'it belongs to another user, so git will not work in it for you' ;;
    private) printf '%s is not public yet' "$ADAPT_URL" ;;
    offline) printf '%s could not be reached from this machine' "$ADAPT_URL" ;;
    *)       printf '%s could not be fetched (git said: %s)' "$ADAPT_URL" \
               "$(printf '%s\n' "$1" | sed '/^[[:space:]]*$/d' | tail -1)" ;;
  esac
}

# The folder that holds the ADCP 1.1 suite, where the interface looks for
# ADapt, or the installation's own tools folder when that one cannot be
# written (the interface looks there too, through ADCP_HOME). A folder that
# is not there yet counts as writable: an install is about to create it.
adapt_tools() {  # $1 = the ADCP 1.1 suite (may be empty), $2 = the installation (ADCP_HOME)
  local side="" own=""
  [ -n "${1:-}" ] && side="$(dirname "$1")"
  [ -n "${2:-}" ] && own="$2/tools/docking"
  if [ -n "$side" ] && { [ -w "$side" ] || [ ! -e "$side" ] || [ -z "$own" ]; }; then
    printf '%s' "$side"
  elif [ -n "$own" ]; then
    printf '%s' "$own"
  fi
  return 0
}

# The folder that holds ADapt for this installation, whether or not it is
# there yet. An existing folder named adapt or ADapt beside the suite or in
# the installation's tools folder is the one kept up to date; otherwise it is
# fetched into adapt_tools.
adapt_dest() {  # $1 = the ADCP 1.1 suite (may be empty), $2 = the installation (ADCP_HOME)
  local d places=() tools
  [ -n "${1:-}" ] && places+=("$(dirname "$1")/adapt" "$(dirname "$1")/ADapt")
  [ -n "${2:-}" ] && places+=("$2/tools/docking/adapt" "$2/tools/docking/ADapt")
  for d in ${places[@]+"${places[@]}"}; do
    if [ -e "$d" ]; then printf '%s' "$d"; return 0; fi
  done
  tools=$(adapt_tools "${1:-}" "${2:-}")
  [ -n "$tools" ] && printf '%s' "$tools/adapt"
  return 0
}

# Whether a folder holds the ADCP 1.1 suite, as the interface decides it
# (adcp11.at): its bin/adcp and ADCP's runADCP.py.
adcp11_is() {
  [ -n "${1:-}" ] && [ -f "$1/bin/adcp" ] || return 1
  ls "$1"/lib/python3*/site-packages/ADCP/runADCP.py >/dev/null 2>&1
}

# The ADCP 1.1 suite, where the interface finds it, in its order
# (adcp11.candidates() without --tools-root): ADCPSUITE_HOME; inside and
# beside the installation (ADCP_HOME) and this checkout, as the folder
# itself, adcpsuite or tools/docking/adcpsuite; then the folder `adcp-pilot
# local` installs it in, and Colab's. The first that holds it, or nothing.
adcp11_found() {  # $1 = the checkout this interface runs from (optional)
  local d h places=() homes=()
  [ -n "${ADCPSUITE_HOME:-}" ] && places+=("$ADCPSUITE_HOME")
  [ -n "${ADCP_HOME:-}" ] && homes+=("$ADCP_HOME")
  [ -n "${1:-}" ] && homes+=("$1")
  for h in ${homes[@]+"${homes[@]}"}; do
    places+=("$h" "$h/adcpsuite" "$h/tools/docking/adcpsuite")
  done
  for h in ${homes[@]+"${homes[@]}"}; do
    places+=("$(dirname "$h")" "$(dirname "$h")/adcpsuite" \
             "$(dirname "$h")/tools/docking/adcpsuite")
  done
  places+=("$HOME/.adcp-pilot/adcpsuite" /content/adcpsuite)
  for d in "${places[@]}"; do
    if adcp11_is "$d"; then printf '%s' "$(cd "$d" && pwd -P)"; return 0; fi
  done
  return 0
}

# Where the interface finds ADapt, in its order (adapt.candidates() without
# --tools-root): the first that is ADapt, or nothing.
adapt_found() {  # $1 = the checkout this interface runs from
  local d h suite places=() homes=()
  [ -n "${ADAPT_HOME:-}" ] && places+=("$ADAPT_HOME")
  suite=$(adcp11_found "$1")
  [ -n "$suite" ] && places+=("$(dirname "$suite")/adapt" "$(dirname "$suite")/ADapt")
  [ -n "${ADCP_HOME:-}" ] && homes+=("$ADCP_HOME")
  homes+=("$1")
  for h in "${homes[@]}"; do
    places+=("$h/tools/docking/adapt" "$h/tools/docking/ADapt" \
             "$(dirname "$h")/tools/docking/adapt" "$(dirname "$h")/tools/docking/ADapt")
  done
  places+=("$(dirname "$1")/adapt" "$(dirname "$1")/ADapt" "$HOME/.adcp-pilot/adapt")
  for d in "${places[@]}"; do
    adapt_is "$d" && { printf '%s' "$d"; return 0; }
  done
  return 0
}

# Its commit, short; nothing for a folder that is not a checkout of its own,
# which would otherwise answer with the commit of a repository around it.
adapt_commit() {
  [ -e "$1/.git" ] && git -C "$1" rev-parse --short HEAD 2>/dev/null
  return 0
}

# Whether the repository can be read anonymously right now. Changes nothing.
adapt_public() {
  local out rc=0 ADAPT_GIT_SECONDS="${ADAPT_GIT_SECONDS:-30}"
  out=$(adapt_git ls-remote --heads "$ADAPT_URL" "$ADAPT_BRANCH" 2>&1) || rc=$?
  if [ "$rc" = 0 ] && [ -n "$out" ]; then
    ADAPT_STATE=public; ADAPT_LINE="$ADAPT_URL is public"
  elif [ "$rc" = 0 ]; then
    ADAPT_STATE=nobranch; ADAPT_LINE="$ADAPT_URL has no branch $ADAPT_BRANCH"
  else
    ADAPT_STATE=unreachable; ADAPT_LINE="$(adapt_why "$out" "$rc")"
  fi
  return 0
}

# A new checkout at $1, which is not there yet. It is cloned beside it under
# another name and moved into place only once it is known to be ADapt, so a
# failed or interrupted clone never leaves a folder that would later be taken
# for an installation.
adapt_clone() {  # $1 = the folder
  local dest="$1" tmp out rc=0 at
  tmp="$(dirname "$dest")/.adapt.fetching.$$"
  if ! mkdir -p "$(dirname "$dest")" 2>/dev/null; then
    ADAPT_STATE=failed
    ADAPT_LINE="not installed: $(dirname "$dest") could not be created"
    return 0
  fi
  # What an earlier, interrupted attempt left, then this one.
  rm -rf "$(dirname "$dest")"/.adapt.fetching.*
  out=$(adapt_git clone -q --branch "$ADAPT_BRANCH" "$ADAPT_URL" "$tmp" 2>&1) || rc=$?
  if [ "$rc" != 0 ]; then
    rm -rf "$tmp"
    ADAPT_STATE=$(adapt_kind "$out" "$rc")
    if [ "$ADAPT_STATE" = private ]; then
      ADAPT_LINE="not installed: $(adapt_why "$out" "$rc"). ADCP docking works without it; ./update.sh fetches it once it is public, or fetch it by hand with access (INSTALL.md, ADapt)"
    else
      ADAPT_LINE="not installed: $(adapt_why "$out" "$rc"). ADCP docking works without it; ./update.sh fetches it when this machine can reach it, or fetch it by hand (INSTALL.md, ADapt)"
    fi
    return 0
  fi
  if ! adapt_is "$tmp"; then
    rm -rf "$tmp"
    ADAPT_STATE=failed
    ADAPT_LINE="not installed: what $ADAPT_URL holds is not ADapt (it has no bin/adapt), so it was not kept"
    return 0
  fi
  if ! mv "$tmp" "$dest" 2>/dev/null; then
    rm -rf "$tmp"
    ADAPT_STATE=failed
    ADAPT_LINE="not installed: it was fetched but could not be put at $dest"
    return 0
  fi
  at=$(adapt_commit "$dest")
  if adapt_embeds "$dest"; then
    ADAPT_STATE=fetched
    ADAPT_LINE="fetched at $at into $dest"
  else
    ADAPT_STATE=fetched_old
    ADAPT_LINE="fetched at $at into $dest, but this version is older than the interface needs (no --base-path), so it is not offered until a newer one is fetched"
  fi
  return 0
}

# An existing folder at $1, brought up to the newest commit on the branch
# when that is safe, and otherwise left exactly as it is. "Safe" means: it is
# a git checkout of its own, on the branch followed, with no local changes
# to its files, and the newest commit only adds to what it has. After moving
# forward it must still be an ADapt this interface can embed; if not, it is
# moved back to the commit it was on, so a change upstream never leaves an
# installation without a working ADapt.
adapt_update() {  # $1 = the folder
  local home="$1" top out rc=0 old was new branch n tag=""
  if [ ! -e "$home/.git" ]; then
    ADAPT_STATE=notgit
    ADAPT_LINE="in $home, which is not a git checkout, so it cannot be updated automatically; left as it is (move it aside and ./update.sh fetches a checkout in its place)"
    return 0
  fi
  # A folder inside another repository answers with that repository's top,
  # so the top is compared with the folder itself.
  top=$(git -C "$home" rev-parse --show-toplevel 2>&1) || rc=$?
  if [ "$rc" != 0 ]; then
    ADAPT_STATE=unreadable
    if [ "$(adapt_kind "$top" "$rc")" = owner ]; then
      ADAPT_LINE="in $home, not updated: it belongs to another user, so git will not work in it for you; its owner's ./update.sh updates it"
    else
      ADAPT_LINE="in $home, not updated: git cannot read it ($(printf '%s\n' "$top" | sed '/^[[:space:]]*$/d' | tail -1)); left as it is"
    fi
    return 0
  fi
  if [ "$(cd "$top" && pwd -P)" != "$(cd "$home" && pwd -P)" ]; then
    ADAPT_STATE=notgit
    ADAPT_LINE="in $home, which is not a git checkout of its own, so it cannot be updated automatically; left as it is"
    return 0
  fi
  old=$(adapt_commit "$home")
  was=$(git -C "$home" rev-parse HEAD 2>/dev/null || true)
  adapt_embeds "$home" || tag=" (older than this interface needs)"
  branch=$(git -C "$home" symbolic-ref -q --short HEAD 2>/dev/null || true)
  if [ -z "$branch" ]; then
    ADAPT_STATE=detached
    ADAPT_LINE="$old$tag, not updated: $home is not on a branch, so it is left as it is"
    return 0
  fi
  if [ "$branch" != "$ADAPT_BRANCH" ]; then
    ADAPT_STATE=branch
    ADAPT_LINE="$old$tag on branch $branch, not updated: only $ADAPT_BRANCH is updated automatically, so it is left as it is"
    return 0
  fi
  if [ -n "$(git -C "$home" --no-optional-locks status --porcelain --untracked-files=no 2>/dev/null)" ]; then
    ADAPT_STATE=changed
    ADAPT_LINE="$old$tag, not updated: $home has local changes to its files, and they are never overwritten; left as it is"
    return 0
  fi
  rc=0
  out=$(adapt_git -C "$home" fetch -q --no-tags "$ADAPT_URL" "refs/heads/$ADAPT_BRANCH" 2>&1) || rc=$?
  if [ "$rc" != 0 ]; then
    ADAPT_STATE=$(adapt_kind "$out" "$rc")
    ADAPT_LINE="$old$tag, kept as it is: $(adapt_why "$out" "$rc")"
    return 0
  fi
  new=$(git -C "$home" rev-parse --short FETCH_HEAD 2>/dev/null || true)
  if [ -z "$new" ] || [ "$was" = "$(git -C "$home" rev-parse FETCH_HEAD 2>/dev/null)" ]; then
    ADAPT_STATE=current
    ADAPT_LINE="$old$tag, the newest on $ADAPT_BRANCH"
    return 0
  fi
  if git -C "$home" merge-base --is-ancestor FETCH_HEAD HEAD 2>/dev/null; then
    ADAPT_STATE=ahead
    ADAPT_LINE="$old$tag, not updated: it has commits of its own beyond $ADAPT_BRANCH on $ADAPT_URL; left as it is"
    return 0
  fi
  if ! git -C "$home" merge-base --is-ancestor HEAD FETCH_HEAD 2>/dev/null; then
    ADAPT_STATE=diverged
    ADAPT_LINE="$old$tag, not updated: its history and $ADAPT_BRANCH on $ADAPT_URL have gone separate ways; left as it is (move it aside and ./update.sh fetches a fresh checkout)"
    return 0
  fi
  n=$(git -C "$home" rev-list --count HEAD..FETCH_HEAD 2>/dev/null || true)
  rc=0
  out=$(git -C "$home" merge -q --ff-only FETCH_HEAD 2>&1) || rc=$?
  if [ "$rc" != 0 ]; then
    # --ff-only checks before it writes, so a refusal changed nothing.
    ADAPT_STATE=blocked
    ADAPT_LINE="$old$tag, not updated: git could not move it to $new ($(printf '%s\n' "$out" | sed '/^[[:space:]]*$/d' | head -1)); left as it was"
    return 0
  fi
  if adapt_is "$home" && adapt_embeds "$home"; then
    ADAPT_STATE=updated
    ADAPT_LINE="$old -> $new ($n new commit(s))"
    return 0
  fi
  # Back to where it was. Nothing local can be lost: it had no local changes.
  if git -C "$home" reset -q --hard "$was" >/dev/null 2>&1 && \
     [ "$(git -C "$home" rev-parse HEAD 2>/dev/null)" = "$was" ]; then
    ADAPT_STATE=rolledback
    ADAPT_LINE="kept at $old$tag: the newest version, $new, does not work with this interface (it has no bin/adapt, or no --base-path), so the update was undone"
  else
    ADAPT_STATE=broken
    ADAPT_LINE="moved to $new, which does not work with this interface, and could not be moved back to $old: run  git -C $home reset --hard $old"
  fi
  return 0
}
