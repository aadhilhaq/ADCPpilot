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
# Beside ADFRsuite it then puts ADCP 1.1, with the Scripps Research
# Institute's own recipe (adcp_dock/adcp11_install.py), and ADapt, Michel
# Sanner's tool for choosing what stays in a receptor, which runs on ADCP 1.1.
# ADapt is fetched anonymously from its public repository; while that
# repository is private this says so in one line and carries on, and
# ./update.sh fetches it once it is public. ADapt never stops an install.
#
#     ./install.sh                      look around, install what is missing
#     ./install.sh <shared-prefix> install under this prefix
#     ./install.sh --check              report and change nothing
#     ./install.sh --no-adcp11 ...      leave out ADCP 1.1 (ADapt needs it)
#
# Run it on a login node: it downloads files, and compute nodes often have no
# outbound network. Nothing else is needed first: a cluster's own conda is
# used when there is one, and micromamba (a single file) is fetched into the
# prefix when there is not. Run it again after a failure and it carries on
# from what is already there. Nothing is written to your home directory.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHECK_ONLY=0
NO_ADCP11=0
PREFIX=""
for a in "$@"; do
  case "$a" in
    --check) CHECK_ONLY=1 ;;
    --no-adcp11) NO_ADCP11=1 ;;
    -*) echo "unknown option: $a" >&2; exit 2 ;;
    *)  PREFIX="$a" ;;
  esac
done

say()  { printf '%s\n' "$*"; }
head_() { printf '\n== %s ==\n' "$*"; }
# A long line, wrapped and indented under its heading.
wrap() { printf '%s\n' "$*" | fold -s -w 74 | sed 's/^/  /'; }

# Finding and fetching ADapt, and finding ADCP 1.1, the way the interface
# finds them. The same file is what update.sh uses.
# shellcheck source=scripts/adapt_fetch.sh
. "$HERE/scripts/adapt_fetch.sh"

# --- where a shared install keeps things ----------------------------------
# Ordered by how likely they are to be the right answer on a given machine.
# The prefix's own come first after anything named outright, so running this
# again with the same prefix finds what the first run built.
ENV_CANDIDATES=(
  "${ADCP_ENV:-}"
  ${PREFIX:+"$PREFIX/envs/adcp"}
  ${ADCP_SHARED:-/opt/shared}/envs/adcp
)
ADFR_CANDIDATES=(
  "${ADFRSUITE_HOME:-}"
  ${PREFIX:+"$PREFIX/tools/docking/ADFRsuite"}
  "$HERE/tools/docking/ADFRsuite"
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

# Both go beside ADFRsuite, where the interface looks for them: beside the one
# found, or the one this is about to install. ADCP_HOME is this checkout, as
# activate.sh will say.
adfr_at="${found_adfr:-${PREFIX:-$HERE}/tools/docking/ADFRsuite}"
TOOLS_AT=$(adapt_tools "$adfr_at" "$HERE")
suite_found() { ADFRSUITE_HOME="$adfr_at" ADCP_HOME="$HERE" adcp11_found; }
head_ "Looking for ADCP 1.1 and ADapt"
found_suite=$(suite_found)
if [ -n "$found_suite" ]; then
  say "  USING ADCP 1.1 at $found_suite"
else
  say "  ADCP 1.1: not found"
fi
ADAPT_DIR=$(adapt_dest "$adfr_at" "$HERE")
adapt_at=$(adapt_commit "$ADAPT_DIR")
if adapt_is "$ADAPT_DIR"; then
  if adapt_embeds "$ADAPT_DIR"; then
    say "  ADapt at $ADAPT_DIR${adapt_at:+ ($adapt_at)}"
  else
    say "  ADapt at $ADAPT_DIR${adapt_at:+ ($adapt_at)}, older than this"
    say "        interface needs; ./update.sh brings it up to date"
  fi
elif [ -e "$ADAPT_DIR" ]; then
  say "  $ADAPT_DIR exists but is not ADapt; it is left as it is"
else
  say "  ADapt: not found"
fi

# --- decide ---------------------------------------------------------------
head_ "What this machine still needs"
need_env=0; need_adfr=0; need_suite=0; need_adapt=0
[ -n "$found_env" ]  || { need_env=1;  say "  an environment with numpy, RDKit, AmberTools and pdb2pqr"; }
[ -n "$found_adfr" ] || { need_adfr=1; say "  ADFRsuite (adcp, agfr, prepare_receptor)"; }
if [ -z "$found_suite" ]; then
  if [ "$NO_ADCP11" = 1 ]; then
    say "  ADCP 1.1 is not here, and --no-adcp11 leaves it out. ADapt, the"
    say "        non-standard and D residues and AGFR 1.1 need it."
  else
    need_suite=1
    say "  ADCP 1.1, into $TOOLS_AT/adcpsuite (about ten minutes"
    say "        and 3 GB; --no-adcp11 leaves it out)"
  fi
fi
if [ ! -e "$ADAPT_DIR" ]; then
  need_adapt=1
  if [ "$CHECK_ONLY" = 1 ]; then
    # Asked, not assumed: --check says what an install would do today.
    adapt_public
    if [ "$ADAPT_STATE" = public ]; then
      say "  ADapt, into $ADAPT_DIR, from $ADAPT_URL"
    else
      wrap "ADapt would not be fetched today: $ADAPT_LINE. ADCP docking works without it, and ./update.sh fetches it once it can."
    fi
  else
    say "  ADapt, into $ADAPT_DIR, if its repository is public"
  fi
fi
[ "$need_env$need_adfr$need_suite$need_adapt" = "0000" ] && say "  nothing. Everything it needs is already on this machine."

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
  # Whatever makes conda environments here: the cluster's own conda or
  # mamba when one is on PATH, else micromamba, a single static file, fetched
  # into the prefix once. So a cluster with no conda module, or one named
  # something this script has never heard of, needs nothing done first.
  SOLVER=""
  for c in "${CONDA_EXE:-}" mamba conda micromamba "$PREFIX/tools/micromamba"; do
    [ -n "$c" ] || continue
    if command -v "$c" >/dev/null 2>&1; then SOLVER="$(command -v "$c")"; break; fi
  done
  if [ -z "$SOLVER" ]; then
    case "$(uname -m)" in
      x86_64) MM=linux-64 ;; aarch64) MM=linux-aarch64 ;; ppc64le) MM=linux-ppc64le ;;
      *) echo "no conda here, and micromamba is not built for $(uname -m)" >&2; exit 1 ;;
    esac
    say "  no conda on PATH; fetching micromamba into $PREFIX/tools"
    mkdir -p "$PREFIX/tools"
    MMURL="https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-$MM"
    if command -v curl >/dev/null 2>&1; then
      curl -fsSL --retry 3 -o "$PREFIX/tools/micromamba.part" "$MMURL"
    else
      wget -q -O "$PREFIX/tools/micromamba.part" "$MMURL"
    fi || { echo "  could not download micromamba from $MMURL" >&2; exit 1; }
    mv "$PREFIX/tools/micromamba.part" "$PREFIX/tools/micromamba"
    chmod 775 "$PREFIX/tools/micromamba"
    SOLVER="$PREFIX/tools/micromamba"
  fi
  say "  using $SOLVER"
  # Into the prefix's package cache, not ~/.conda, for the same reason the
  # prefix is not optional. HOME itself points into the prefix for these
  # commands only: conda keeps a list of every environment it has made in
  # ~/.conda/environments.txt, and nothing about this install belongs there.
  # A cache already named by whoever ran this (HOPE Labs' Set up names its own)
  # is used rather than a second one beside it.
  export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$PREFIX/conda_shared/pkgs}"
  export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$PREFIX/conda_shared/mamba}"
  mkdir -p "$CONDA_PKGS_DIRS" "$PREFIX/conda_shared/home"
  SPECS=(python=3.11 "numpy<2" rdkit ambertools pdb2pqr pdbfixer openmm biopython)
  say "  this is the slow part: 15-40 minutes, almost all of it AmberTools"
  if [ -d "$ENVDIR/conda-meta" ]; then
    # An environment an earlier run began and did not finish: completed
    # in place rather than made again, so nothing already there is fetched
    # twice.
    say "  $ENVDIR was begun by an earlier run; completing it"
    HOME="$PREFIX/conda_shared/home" "$SOLVER" install -y -p "$ENVDIR" \
        -c conda-forge --override-channels "${SPECS[@]}"
  else
    HOME="$PREFIX/conda_shared/home" "$SOLVER" create -y -p "$ENVDIR" \
        -c conda-forge --override-channels "${SPECS[@]}"
  fi
  # The packages fetched to build it are not needed once it is built, and
  # they are tens of thousands of files.
  HOME="$PREFIX/conda_shared/home" "$SOLVER" clean -a -y >/dev/null 2>&1 || true
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
    if command -v curl >/dev/null 2>&1; then
      curl -fSL --retry 3 -o "$TAR" "$URL"
    else
      wget -q -O "$TAR" "$URL"
    fi || {
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
  # tar is still listing when head has its line and closes the pipe, so tar
  # dies of SIGPIPE; under pipefail that would end the script here, so its
  # status is not the pipeline's.
  SRC=$( { tar tzf "$TAR" 2>/dev/null || true; } | head -1 | cut -d/ -f1)
  [ -n "$SRC" ] || { echo "  $TAR is not a readable tar.gz; delete it and retry" >&2; exit 1; }
  tar xzf "$TAR"
  [ -d "$SRC" ] || { echo "  $TAR did not unpack into $SRC" >&2; exit 1; }
  cd "$SRC"

  say "  running the Scripps installer (a few minutes)"
  LOG="$TOOLS/adfr_install.log"
  # No licence answer is given on the user's behalf. The Scripps installer
  # asks "ACADEMIC INSTALLATION (Y/N)?" about MSMS, a licence declaration
  # only the installing user can make, and renames mslibACA or mslibCOM to
  # mslib, the importable MSMS library, once it is answered. This tool never
  # imports MSMS, so -l (the installer's own switch) skips the question and
  # with it that rename. It therefore governs a tree built here from scratch
  # and says nothing about one already at this path: unpacking over a tree
  # built without -l leaves the library where it was, and the msms program
  # is in bin/msms either way. What this tree holds is reported below.
  # stdin is /dev/null so no prompt can block the install.
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
  # Read back, rather than assumed from -l: CCSBpckgs/mslib exists only once
  # the academic-or-commercial question has been answered and the rename it
  # drives has happened.
  if [ -d "$found_adfr/CCSBpckgs/mslib" ]; then
    say "  the MSMS library is in place at CCSBpckgs/mslib, so its academic"
    say "  or commercial terms were accepted for this tree; this script passes"
    say "  -l, which skips that question, so it was left by an earlier install"
    say "  or put there by hand. Nothing adcp-dock runs imports it."
  elif [ -d "$found_adfr/CCSBpckgs/mslibACA" ]; then
    say "  the MSMS library was not activated: CCSBpckgs still holds mslibACA"
    say "  and mslibCOM under those names, so nothing can import mslib."
  else
    say "  this tree holds no MSMS library to import."
  fi
  if [ -f "$found_adfr/bin/msms" ]; then
    say "  the msms program itself is in bin/msms, as in every installation"
  fi
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

# --- ADCP 1.1, beside ADFRsuite --------------------------------------------
# With the package's own installer, which follows the Scripps Research
# Institute's recipe and checks every Scripps package against a pinned
# checksum. Optional for docking, which ADCP 1.0 does, so a failure here is
# reported and the install carries on; running this script again retries it.
adfr_at="$found_adfr"
TOOLS_AT=$(adapt_tools "$found_adfr" "$HERE")
if [ "$need_suite" = 1 ]; then
  SUITE_DEST="$TOOLS_AT/adcpsuite"
  head_ "Installing ADCP 1.1 into $SUITE_DEST"
  say "  the Scripps Research Institute's recipe: about ten minutes and 3 GB"
  # HOME points beside it for this one command, as for conda above: nothing
  # of it belongs in the home directory. A micromamba fetched above is
  # reused rather than fetched again.
  mkdir -p "$TOOLS_AT/.adcpsuite_home"
  if HOME="$TOOLS_AT/.adcpsuite_home" ADFRSUITE_HOME="$found_adfr" ADCP_HOME="$HERE" \
     PATH="${PREFIX:+$PREFIX/tools:}$PATH" PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}" \
     "$found_env/bin/python" -u -m adcp_dock.adcp11_install --dest "$SUITE_DEST" 2>&1 \
     | sed 's/^/  /'; then :; fi
  rmdir "$TOOLS_AT/.adcpsuite_home" 2>/dev/null || true
  found_suite=$(suite_found)
  if [ -n "$found_suite" ]; then
    say "  verified: ADCP 1.1 is at $found_suite"
  else
    cat >&2 <<MSG

  ADCP 1.1 did not install; the lines above say why. Docking with ADCP 1.0
  works without it. Run this script again to retry, or just this step:

      $found_env/bin/python -m adcp_dock.cli install-adcp11 --dest $SUITE_DEST

MSG
  fi
fi

# --- ADapt, beside it -------------------------------------------------------
# Fetched anonymously (scripts/adapt_fetch.sh): no prompt and no stored
# credential, so while the repository is private this is one calm line. A
# folder already there is left for ./update.sh, which knows when it is safe
# to move one forward.
head_ "ADapt"
ADAPT_DIR=$(adapt_dest "$found_adfr" "$HERE")
if [ -z "$ADAPT_DIR" ]; then
  say "  not fetched: there is no ADFRsuite folder to put it beside"
elif [ -e "$ADAPT_DIR" ]; then
  adapt_at=$(adapt_commit "$ADAPT_DIR")
  if adapt_is "$ADAPT_DIR"; then
    wrap "already at $ADAPT_DIR${adapt_at:+ ($adapt_at)}; ./update.sh keeps it up to date"
  else
    wrap "$ADAPT_DIR is already there and is not ADapt; it is left as it is"
  fi
else
  adapt_clone "$ADAPT_DIR"
  wrap "$ADAPT_LINE"
fi
if [ -z "$found_suite" ] && adapt_is "${ADAPT_DIR:-/nonexistent}"; then
  wrap "ADapt runs on ADCP 1.1, which is not installed, so the interface will not offer it until ADCP 1.1 is added."
fi

# --- prove it ------------------------------------------------------------
head_ "Checking"
# shellcheck disable=SC1090
source "$ACT"
"$found_env/bin/python" -m adcp_dock.cli check || true
say ""
say "  ADCP 1.1         ${found_suite:-not installed (optional; ADapt needs it)}"
if adapt_is "${ADAPT_DIR:-/nonexistent}"; then
  adapt_at=$(adapt_commit "$ADAPT_DIR")
  if adapt_embeds "$ADAPT_DIR"; then
    say "  ADapt            $ADAPT_DIR${adapt_at:+ ($adapt_at)}"
  else
    say "  ADapt            $ADAPT_DIR${adapt_at:+ ($adapt_at)}, older than the interface needs"
  fi
else
  say "  ADapt            not installed (optional; ./update.sh fetches it once it is public)"
fi

cat <<EOF

Next:

    source $ACT
    adcp-dock check                 what it found on this machine
    adcp-dock serve --runs \$SCRATCH/adcp_runs

Never dock on a login node. Everything heavy goes to SLURM, and adcp-dock
submits it there for you.
EOF
