#!/usr/bin/env bash
# Install GUI to ADCP Docking on a cluster, reusing whatever is already there.
#
# By Aadhil Haq
#
# The tools it needs (ADFRsuite for the docking, AmberTools and pdb2pqr for
# the rescoring, numpy and RDKit for everything else) are often already
# present in a shared scientific environment. Another copy of AmberTools is
# about 60,000 files and 15 to 40 minutes of work, and on a shared /scratch
# the file count matters more than the size. So this looks for them first,
# in the places a shared install puts them, and builds only what is missing.
#
#     ./install.sh                      look around, install what is missing
#     ./install.sh <shared-prefix> install under this prefix
#     ./install.sh --check              report and change nothing
#
# Run it on a login node: it downloads files, and compute nodes often have no
# outbound network.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHECK_ONLY=0
PREFIX=""
for a in "$@"; do
  case "$a" in
    --check) CHECK_ONLY=1 ;;
    -*) echo "unknown option: $a" >&2; exit 2 ;;
    *)  PREFIX="$a" ;;
  esac
done

say()  { printf '%s\n' "$*"; }
head_() { printf '\n== %s ==\n' "$*"; }

# --- where a shared install keeps things ----------------------------------
# Ordered by how likely they are to be the right answer on a given machine.
ENV_CANDIDATES=(
  "${ADCP_ENV:-}"
  ${ADCP_SHARED:-/opt/shared}/envs/adcp
)
ADFR_CANDIDATES=(
  "${ADFRSUITE_HOME:-}"
  ${ADCP_SHARED:-/opt/shared}/tools/docking/ADFRsuite
  ${ADCP_SHARED:-/opt/shared}/third_party/ADFRsuite
)

# --- what has to exist, and what each one is for --------------------------
PY_MODULES=(numpy rdkit)
PY_OPTIONAL=(pdbfixer openmm Bio)
AMBER_BINS=(tleap sander MMPBSA.py ante-MMPBSA.py pdb2pqr)
ADFR_BINS=(adcp agfr prepare_receptor)

found_env=""
found_adfr=""

head_ "Looking for an environment that already has what this needs"
for e in "${ENV_CANDIDATES[@]}"; do
  [ -n "$e" ] || continue
  [ -x "$e/bin/python" ] || continue
  missing=""
  for m in "${PY_MODULES[@]}"; do
    "$e/bin/python" -c "import $m" >/dev/null 2>&1 || missing="$missing $m"
  done
  for b in "${AMBER_BINS[@]}"; do
    [ -x "$e/bin/$b" ] || missing="$missing $b"
  done
  if [ -z "$missing" ]; then
    say "  USING $e"
    say "        it already has numpy, RDKit, AmberTools and pdb2pqr"
    found_env="$e"
    break
  else
    say "  $e"
    say "        missing:$missing"
  fi
done
[ -n "$found_env" ] || say "  nothing complete found"

head_ "Looking for ADFRsuite"
for d in "${ADFR_CANDIDATES[@]}"; do
  [ -n "$d" ] || continue
  ok=1
  for b in "${ADFR_BINS[@]}"; do
    [ -x "$d/bin/$b" ] || ok=0
  done
  if [ "$ok" = 1 ]; then
    say "  USING $d"
    found_adfr="$d"
    break
  fi
  [ -d "$d" ] && say "  $d exists but is missing one of: ${ADFR_BINS[*]}"
done
if [ -z "$found_adfr" ]; then
  if command -v adcp >/dev/null 2>&1 && command -v agfr >/dev/null 2>&1; then
    found_adfr="$(dirname "$(dirname "$(command -v adcp)")")"
    say "  USING $found_adfr (from PATH)"
  else
    say "  not found"
  fi
fi

# --- decide ---------------------------------------------------------------
head_ "What this machine still needs"
need_env=0; need_adfr=0
[ -n "$found_env" ]  || { need_env=1;  say "  an environment with numpy, RDKit, AmberTools and pdb2pqr"; }
[ -n "$found_adfr" ] || { need_adfr=1; say "  ADFRsuite (adcp, agfr, prepare_receptor)"; }
[ "$need_env$need_adfr" = "00" ] && say "  nothing. Everything it needs is already on this machine."

if [ "$CHECK_ONLY" = 1 ]; then
  say ""
  say "--check: nothing was changed."
  exit 0
fi

# --- the one file that sets everything up ---------------------------------
# Written as soon as the environment is known, and again at the end. The
# launcher sources it to find a Python new enough to run this package, so it
# must exist even if installing ADFRsuite fails partway.
ACT="$HERE/activate.sh"
write_activate() {
  # $1 = the environment, $2 = ADFRsuite (may be empty).
  #
  # Written line by line rather than with ${2:+...} inside a heredoc, whose
  # quotes would need escaping twice (for the heredoc and for the shell that
  # later sources the file); a wrong escape puts literal quote marks into
  # the variable's value.
  local env_dir="$1" adfr="$2"
  {
    echo "# Written by install.sh on $(date +%Y-%m-%d). Source it, do not run it."
    echo "#"
    echo "#     source $ACT"
    echo "#"
    echo "# Nothing here is specific to one person: it is read by whoever is"
    echo "# using the shared checkout, and by the desktop launcher every time"
    echo "# it connects, which is why ADCP_PYTHON matters: a cluster's system"
    echo "# python3 may be too old to parse this package."
    echo "export ADCP_HOME=\"$HERE\""
    if [ -n "$env_dir" ]; then
      echo "export ADCP_PYTHON=\"$env_dir/bin/python3\""
      echo "export AMBERHOME=\"$env_dir\""
      echo "export PATH=\"$env_dir/bin:\$PATH\""
    fi
    if [ -n "$adfr" ]; then
      echo "export ADFRSUITE_HOME=\"$adfr\""
      echo "export PATH=\"$adfr/bin:\$PATH\""
    fi
    echo "export PYTHONPATH=\"\$ADCP_HOME:\${PYTHONPATH:-}\""
  } > "$ACT"
}

# --- build only what is missing -------------------------------------------
if [ "$need_env" = 1 ]; then
  if [ -z "$PREFIX" ]; then
    cat >&2 <<'MSG'

No usable environment was found and no prefix was given, so this does not
know where to put one.

    ./install.sh <shared-prefix>

THE PREFIX IS NOT OPTIONAL on a shared cluster. Without one a conda
environment lands under your home directory, whose file allowance is small,
and the build dies partway with "Disk quota exceeded" having already spent
twenty minutes. envs/adcp is appended to whatever is given.

If a suitable environment already exists somewhere this did not look, point
at it instead and nothing is built at all:

    ADCP_ENV=/path/to/env ./install.sh

MSG
    exit 1
  fi
  ENVDIR="$PREFIX/envs/adcp"
  head_ "Building $ENVDIR"
  command -v conda >/dev/null 2>&1 || {
    echo "conda is not on PATH. 'module load Anaconda3' first." >&2; exit 1; }
  # Into the shared package cache, not ~/.conda, for the same reason the
  # prefix is not optional. Given to this command only, so no setting is
  # written into the home directory either.
  mkdir -p "$PREFIX/conda_shared/pkgs"
  export CONDA_PKGS_DIRS="$PREFIX/conda_shared/pkgs"
  say "  this is the slow part: 15-40 minutes, almost all of it AmberTools"
  conda create -y -p "$ENVDIR" -c conda-forge \
      python=3.11 "numpy<2" rdkit ambertools pdb2pqr \
      pdbfixer openmm biopython
  found_env="$ENVDIR"
fi

# The environment is settled now, whether it was found or just built. Write
# activate.sh here so that even if the next step fails the launcher can
# still find a usable Python.
if [ -n "$found_env" ]; then
  write_activate "$found_env" "$found_adfr"
  say ""
  say "Wrote $ACT (the launcher reads this to find the right Python)"
fi

if [ "$need_adfr" = 1 ]; then
  TOOLS="${PREFIX:-$HERE}/tools/docking"
  head_ "Installing ADFRsuite into $TOOLS"
  mkdir -p "$TOOLS"
  cd "$TOOLS"
  URL=https://ccsb.scripps.edu/adfr/download/1038/
  TAR=ADFRsuite_x86_64Linux_1.0.tar.gz
  if [ ! -s "$TAR" ]; then
    say "  downloading from Scripps (about 100 MB)"
    curl -fSL -o "$TAR" "$URL" || {
      cat >&2 <<'MSG'
  Could not download ADFRsuite.

  This is a login-node step: compute nodes here have no outbound network. If
  you are on one, that is the whole explanation.

  Otherwise fetch it by hand from https://ccsb.scripps.edu/adfr/downloads/,
  drop ADFRsuite_x86_64Linux_1.0.tar.gz in this directory, and run this
  script again; it will use the file that is already here.
MSG
      exit 1; }
  else
    say "  $TAR is already here, using it"
  fi

  # The same pinned checksum as adcp_dock/adfr.py: the archive is run as an
  # installer, so a truncated or substituted download must stop here.
  SHA=66539864ffb7ba87728d52ee50dde78586f551c7ca96aa455dae761875faf7e4
  got=$(sha256sum "$TAR" 2>/dev/null | cut -d' ' -f1)
  [ -n "$got" ] || got=$(shasum -a 256 "$TAR" | cut -d' ' -f1)
  if [ "$got" != "$SHA" ]; then
    echo "  $TAR does not match its published checksum (got $got)." >&2
    echo "  Delete it and run this script again." >&2
    exit 1
  fi
  say "  checksum verified"

  # What the archive unpacks into, read from the archive rather than assumed:
  # a wrong guess `cd`s into a directory that is not there and, under
  # `set -e`, ends the script with no message, immediately after a 100 MB
  # download.
  say "  unpacking"
  SRC=$(tar tzf "$TAR" | head -1 | cut -d/ -f1)
  [ -n "$SRC" ] || { echo "  $TAR is not a readable tar.gz; delete it and retry" >&2; exit 1; }
  tar xzf "$TAR"
  [ -d "$SRC" ] || { echo "  $TAR did not unpack into $SRC" >&2; exit 1; }
  cd "$SRC"

  say "  running the Scripps installer (a few minutes)"
  LOG="$TOOLS/adfr_install.log"
  # No licence answer is given on the user's behalf. The Scripps installer
  # asks "ACADEMIC INSTALLATION (Y/N)?" about MSMS, a licence declaration
  # only the installing user can make. This tool never calls MSMS, so -l
  # (the installer's own switch) skips the question and the MSMS library is
  # never activated; stdin is /dev/null so no prompt can block the install.
  if ! ./install.sh -d "$TOOLS/ADFRsuite" -c 0 -l < /dev/null > "$LOG" 2>&1; then
    say "  the installer returned non-zero; checking what it produced anyway"
  fi

  cd "$HERE"
  missing=""
  for b in "${ADFR_BINS[@]}"; do
    [ -x "$TOOLS/ADFRsuite/bin/$b" ] || missing="$missing $b"
  done
  if [ -n "$missing" ]; then
    cat >&2 <<MSG

  ADFRsuite did not install properly. Missing:$missing

  Its own log is at $LOG; the last few lines usually say why:

      tail -20 $LOG

  A part-installed tree is worse than none, because the next run would find
  it and believe it. Remove it before retrying:

      rm -rf $TOOLS/ADFRsuite

MSG
    exit 1
  fi
  found_adfr="$TOOLS/ADFRsuite"
  say "  installed at $found_adfr"
  say "  verified: ${ADFR_BINS[*]} are all present and executable"
fi

# --- install this package into that environment ---------------------------
head_ "Installing adcp-dock into $found_env"
# --no-deps, because the environment may be shared: a plain `pip install -e .`
# could upgrade packages other users depend on (numpy, for example). The
# dependencies are checked above instead; if numpy or RDKit were missing,
# this script would already have built an environment of its own.
if ! "$found_env/bin/python" -m pip install --quiet --no-deps --no-cache-dir -e "$HERE"; then
  cat >&2 <<MSG

  pip could not install into $found_env.

  If that is because the directory is not writable by you, it is a shared
  environment and that is correct. Two ways on:

      export PYTHONPATH=$HERE:\$PYTHONPATH      use it in place
      ./install.sh \$SCRATCH                     build one of your own

  The first needs nothing from anyone and is what activate.sh sets up
  anyway; adcp-dock runs perfectly well from the checkout.
MSG
  say "  continuing without installing the console script"
else
  say "  done, with --no-deps: nothing in that shared environment was"
  say "  upgraded, so it stays usable for other users"
fi

# --- and again, now that ADFRsuite is known -------------------------------
write_activate "$found_env" "$found_adfr"
say ""
say "Wrote $ACT"

# --- prove it ------------------------------------------------------------
head_ "Checking"
# shellcheck disable=SC1090
source "$ACT"
"$found_env/bin/python" -m adcp_dock.cli check || true

cat <<EOF

Next:

    source $ACT
    adcp-dock check                 what it found on this machine
    adcp-dock serve --runs \$SCRATCH/adcp_runs

Never dock on a login node. Everything heavy goes to SLURM, and adcp-dock
submits it there for you.
EOF
