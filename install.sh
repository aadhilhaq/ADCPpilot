#!/usr/bin/env bash
# Install GUI to ADCP Docking on a cluster, reusing whatever is already there.
#
# By Aadhil Haq
#
# The tools it needs (ADFRsuite for the docking, AmberTools and pdb2pqr for
# the rescoring, numpy and RDKit for everything else, and OpenMM, PDBFixer
# and Biopython for preparing the receptor) are often already on a cluster,
# in a shared scientific environment or as a module. A new environment with
# AmberTools is about 49,000 files, and on a shared /scratch the file count
# matters more than the size. So this looks for them first, adds only what
# is missing, and builds only what it cannot find.
#
# Beside ADFRsuite it then puts ADCP 1.1, with the Scripps Research
# Institute's own recipe (adcp_dock/adcp11_install.py), and ADapt, Michel
# Sanner's tool for choosing what stays in a receptor, which runs on ADCP 1.1.
# ADapt is fetched anonymously from its public repository; while that
# repository is private this says so in one line and carries on, and
# ./update.sh fetches it once it is public. ADapt never stops an install.
#
#     ./install.sh <folder>             install into this folder
#     ./install.sh                      use what is already here
#     ./install.sh --check [<folder>]   report what it would do, change nothing
#     ./install.sh --no-adcp11 ...      leave out ADCP 1.1 (ADapt needs it)
#
# Where it looks:
#   an environment: ADCP_ENV, <folder>/envs/adcp, the one an earlier run
#     wrote into activate.sh, and $ADCP_SHARED/envs/adcp when that is set;
#   AmberTools: in that environment, then wherever AMBERHOME or PATH points,
#     which is where `module load AmberTools` puts a cluster's own;
#   ADFRsuite: ADFRSUITE_HOME, activate.sh, <folder>/tools/docking,
#     this checkout's tools/docking, $ADCP_SHARED, then PATH.
#
# An environment it did not build is never rebuilt or upgraded. A package it
# lacks is added beside what is there, with everything already installed held
# where it is, and only when you can write to it; otherwise this says what is
# missing and what that turns off.
#
# Run it on a login node: it downloads files, and compute nodes often have no
# outbound network. Nothing else is needed first: a cluster's own conda is
# used when it is recent, and micromamba (a single file) is fetched into the
# folder when it is not. Run it again after a failure and it carries on from
# what is already there. It never asks a question, and nothing is written to
# your home directory.

# Run, not sourced, and by bash. These lines are read by any shell.
if [ -z "${BASH_VERSION:-}" ]; then
  echo "Run this with bash:  ./install.sh <folder>   (or: bash install.sh <folder>)" >&2
  return 2 2>/dev/null || exit 2
fi
if [ "${BASH_SOURCE[0]}" != "$0" ]; then
  echo "Run this script, do not source it:  ./install.sh <folder>" >&2
  return 2
fi
set -Eeuo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHECK_ONLY=0
NO_ADCP11=0
PREFIX=""

usage() {
  cat <<'EOF'
Usage:
    ./install.sh <folder>             install into this folder
    ./install.sh                      use what is already here
    ./install.sh --check [<folder>]   report what it would do, change nothing
    ./install.sh --no-adcp11 ...      leave out ADCP 1.1 (ADapt needs it)

Give a folder with room for about 150,000 files while it installs (about
93,000 stay), in your scratch or project space rather than your home
directory. INSTALL.md explains each step.
EOF
}

for a in "$@"; do
  case "$a" in
    --check) CHECK_ONLY=1 ;;
    --no-adcp11) NO_ADCP11=1 ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "install.sh: there is no option $a" >&2; usage >&2; exit 2 ;;
    *)
      if [ -n "$PREFIX" ]; then
        echo "install.sh: give one folder, not two ($PREFIX and $a)" >&2; exit 2
      fi
      PREFIX="$a" ;;
  esac
done
# A relative folder is made absolute, since activate.sh is read from
# anywhere: by the launcher, by every job, and by you in a new terminal.
if [ -n "$PREFIX" ]; then
  case "$PREFIX" in /*) ;; *) PREFIX="$PWD/$PREFIX" ;; esac
  while [ "${PREFIX%/}" != "$PREFIX" ] && [ "$PREFIX" != / ]; do PREFIX="${PREFIX%/}"; done
fi

say()  { printf '%s\n' "$*"; }
head_() { printf '\n== %s ==\n' "$*"; }
# A long line, wrapped and indented under its heading. Only between words: a
# long path stays whole, so it can still be copied.
wrap() {
  printf '%s\n' "$*" | awk -v w=74 '{
    line = ""; n = split($0, a, " ")
    for (i = 1; i <= n; i++) {
      if (line == "") line = a[i]
      else if (length(line) + 1 + length(a[i]) <= w) line = line " " a[i]
      else { print "  " line; line = a[i] }
    }
    if (line != "") print "  " line
  }'
}
mins() { local s=$(( $(date +%s) - $1 )); printf '%d min %02d s' $((s / 60)) $((s % 60)); }

# Every Python this runs is started without PYTHONPATH or PYTHONHOME. A
# cluster module (AmberTools is one) points PYTHONPATH at its own Python's
# packages, and an environment's Python then imports those instead of its
# own and fails, numpy first.
py_clean() { env -u PYTHONPATH -u PYTHONHOME "$@"; }

# --- when something fails ---------------------------------------------------
# Under `set -e` a failed command ends the script; this says where, and what
# to do, instead of leaving a bare error from conda or curl as the last word.
T0=$(date +%s)
STEP="starting"
on_error() {
  local rc=$1 line=$2
  # Inside $(...) a failure is the function's own business: the ADapt and
  # room helpers return 0 whatever their commands do, and are relied on to.
  [ "${BASH_SUBSHELL:-0}" -gt 0 ] && return 0
  trap - ERR
  {
    printf '\ninstall.sh STOPPED while %s (exit %s, line %s).\n\n' "$STEP" "$rc" "$line"
    wrap "The lines just above say why. Nothing is lost: run the same command again and it carries on from what is already there."
    say ""
    say "  If it stops in the same place again:"
    wrap "* \"Could not resolve host\", \"Connection timed out\" or a failed download: run it on a login node, since compute nodes usually have no internet."
    wrap "* \"Disk quota exceeded\" or \"No space left on device\": the folder has run out of room or of files. Your cluster's quota command (often showquota, or quota -s) shows how much you have; choose a folder with more room and run it again there."
    wrap "* anything else: copy the last 30 lines of this output into an issue at https://github.com/aadhilhaq/ADCP_docking/issues"
  } >&2
  exit "$rc"
}
trap 'on_error $? $LINENO' ERR

# Finding and fetching ADapt, and finding ADCP 1.1, the way the interface
# finds them. The same file is what update.sh uses.
# shellcheck source=scripts/adapt_fetch.sh
. "$HERE/scripts/adapt_fetch.sh"

# --- what it needs, and what each part is for --------------------------------
PY_REQUIRED=(numpy rdkit)
PY_OPTIONAL=(openmm pdbfixer Bio)
AMBER_BINS=(tleap sander MMPBSA.py ante-MMPBSA.py)
ADFR_BINS=(adcp agfr prepare_receptor)

name_of() {
  case "$1" in
    rdkit) echo RDKit ;; openmm) echo OpenMM ;; pdbfixer) echo PDBFixer ;;
    Bio) echo Biopython ;; amber) echo AmberTools ;; *) echo "$1" ;;
  esac
}
what_for() {
  case "$1" in
    numpy|rdkit) echo "everything: adcp-dock does not start without it" ;;
    pdb2pqr) echo "protonation with PROPKA at a chosen pH, and the fallback when reduce cannot run" ;;
    openmm) echo "the receptor's restrained minimisation" ;;
    pdbfixer) echo "rebuilding missing atoms and side chains in the receptor" ;;
    Bio) echo "reading mmCIF structure files" ;;
    amber) echo "MM-GBSA rescoring" ;;
  esac
}
conda_spec() {
  case "$1" in
    numpy) echo "numpy<2" ;; Bio) echo biopython ;; amber) echo ambertools ;; *) echo "$1" ;;
  esac
}
pip_spec() {  # nothing for AmberTools, which pip does not have
  case "$1" in
    numpy) echo "numpy<2" ;; Bio) echo biopython ;; amber) ;; *) echo "$1" ;;
  esac
}
names() { local n out=""; for n in "$@"; do out="$out, $(name_of "$n")"; done; printf '%s' "${out#, }"; }

# How many files and how much space each part needs while it is installed,
# and what stays once its download cache is removed. Measured on a clean
# install on Grace (Texas A&M, October 2026): the environment 48,927 files,
# ADFRsuite 11,272, ADCP 1.1 31,116, and their caches 50,579 more, 143,472
# in all; the peak is rounded up to 150,000. So a folder without room is
# found before anything is downloaded rather than part way through.
ROOM_ENV_FILES=90000;    ROOM_ENV_GB=5;    KEEP_ENV_FILES=49000;    KEEP_ENV_GB=3
ROOM_TOPUP_FILES=10000;  ROOM_TOPUP_GB=1;  KEEP_TOPUP_FILES=5000;   KEEP_TOPUP_GB=1
ROOM_ADFR_FILES=12000;   ROOM_ADFR_GB=1;   KEEP_ADFR_FILES=11500;   KEEP_ADFR_GB=1
ROOM_SUITE_FILES=48000;  ROOM_SUITE_GB=4;  KEEP_SUITE_FILES=31500;  KEEP_SUITE_GB=3

# --- this computer ------------------------------------------------------------
head_ "This computer"
say "  $(uname -s) $(uname -m), $(hostname 2>/dev/null || echo 'name unknown')"
if [ "$(uname -s)" != Linux ] || [ "$(uname -m)" != x86_64 ]; then
  wrap "This script is for a Linux cluster on x86-64, which is what ADFRsuite is built for. On your own Mac or Windows computer, follow INSTALL.md section 2 instead."
  [ "$CHECK_ONLY" = 1 ] || exit 1
fi
lacking=""
for c in tar gzip; do command -v "$c" >/dev/null 2>&1 || lacking="$lacking $c"; done
command -v curl >/dev/null 2>&1 || command -v wget >/dev/null 2>&1 || lacking="$lacking curl"
command -v sha256sum >/dev/null 2>&1 || command -v shasum >/dev/null 2>&1 || lacking="$lacking sha256sum"
if [ -n "$lacking" ]; then
  wrap "These standard programs are not on this machine:$lacking. Every Linux cluster has them, so this is probably not a login node; log in to one and run this again there. If it is, ask your cluster's support to add them."
  [ "$CHECK_ONLY" = 1 ] || exit 1
fi
command -v git >/dev/null 2>&1 || wrap "git is not here, so ADapt cannot be fetched. Everything else works without it."
if [ -n "${SLURM_JOB_ID:-}" ]; then
  wrap "This is running inside Slurm job $SLURM_JOB_ID. Compute nodes often cannot download anything; if a download fails, run it again on a login node."
fi
if [ -n "$PREFIX" ] && [ -n "${HOME:-}" ] && [ "$HOME" != / ]; then
  case "$PREFIX/" in
    "${HOME%/}"/*)
      wrap "$PREFIX is in your home directory. Most clusters allow only a few thousand files there, and this can need about 150,000 while it installs; a folder in your scratch or project space is usually the right place." ;;
  esac
fi

# --- what an earlier run recorded --------------------------------------------
ACT="$HERE/activate.sh"
recorded() {  # the value an earlier activate.sh gave one variable
  [ -f "$ACT" ] || return 0
  sed -n "s/^export $1=\"\(.*\)\"\$/\1/p" "$ACT" 2>/dev/null | head -n 1 || true
}
PREV_ENV=""
prev_py=$(recorded ADCP_PYTHON)
[ -n "$prev_py" ] && PREV_ENV="$(dirname "$(dirname "$prev_py")")"
PREV_ADFR=$(recorded ADFRSUITE_HOME)

# --- an environment ------------------------------------------------------------
# env_inspect DIR sets: E_OK (1 when its Python is 3.9 or newer), E_WHY (why
# not), E_REQ and E_OPT (what it lacks), E_AMBER (1 when AmberTools is in it),
# E_KIND (conda, venv or other) and E_WRITABLE (1 when you can add to it).
env_inspect() {
  local e=$1 m b site
  E_OK=0; E_WHY=""; E_REQ=""; E_OPT=""; E_AMBER=0; E_KIND=other; E_WRITABLE=0
  if [ ! -e "$e" ]; then E_WHY="not there"; return 0; fi
  if [ ! -x "$e/bin/python" ]; then E_WHY="has no bin/python"; return 0; fi
  if ! py_clean "$e/bin/python" -c "import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)" \
       >/dev/null 2>&1; then
    E_WHY="its Python does not start, or is older than 3.9"; return 0
  fi
  E_OK=1
  for m in "${PY_REQUIRED[@]}"; do
    py_clean "$e/bin/python" -c "import $m" >/dev/null 2>&1 || E_REQ="$E_REQ $m"
  done
  [ -x "$e/bin/pdb2pqr" ] || E_REQ="$E_REQ pdb2pqr"
  for m in "${PY_OPTIONAL[@]}"; do
    py_clean "$e/bin/python" -c "import $m" >/dev/null 2>&1 || E_OPT="$E_OPT $m"
  done
  E_AMBER=1
  for b in "${AMBER_BINS[@]}"; do [ -x "$e/bin/$b" ] || E_AMBER=0; done
  if [ -d "$e/conda-meta" ]; then E_KIND=conda
  elif [ -f "$e/pyvenv.cfg" ]; then E_KIND=venv
  fi
  site=$(py_clean "$e/bin/python" -c "import sysconfig; print(sysconfig.get_paths()['purelib'])" \
         2>/dev/null || true)
  if [ -w "$e" ] && [ -w "$e/bin" ] && { [ -z "$site" ] || [ -w "$site" ]; } \
     && { [ "$E_KIND" != conda ] || [ -w "$e/conda-meta" ]; }; then
    E_WRITABLE=1
  fi
  return 0
}

ENV_CANDIDATES=()
for e in "${ADCP_ENV:-}" "${PREFIX:+$PREFIX/envs/adcp}" "$PREV_ENV" \
         "${ADCP_SHARED:+$ADCP_SHARED/envs/adcp}"; do
  [ -n "$e" ] || continue
  case " ${ENV_CANDIDATES[*]-} " in *" $e "*) continue ;; esac
  ENV_CANDIDATES+=("$e")
done

head_ "Looking for a Python environment"
found_env=""
for e in ${ENV_CANDIDATES[@]+"${ENV_CANDIDATES[@]}"}; do
  env_inspect "$e"
  if [ "$E_OK" = 1 ]; then
    found_env="$e"
    break
  fi
  say "  $e: $E_WHY"
done
env_req=""; env_opt=""; env_amber=0; env_kind=""; env_writable=0
if [ -n "$found_env" ]; then
  env_req=$E_REQ; env_opt=$E_OPT; env_amber=$E_AMBER; env_kind=$E_KIND; env_writable=$E_WRITABLE
  say "  FOUND $found_env"
  case "$env_kind" in
    conda) kind_text="a conda environment" ;; venv) kind_text="a Python venv" ;;
    *) kind_text="a Python installation" ;;
  esac
  say "        $(py_clean "$found_env/bin/python" -V 2>&1 | head -n 1 || true), $kind_text that you $([ "$env_writable" = 1 ] && echo can || echo cannot) add to"
  # shellcheck disable=SC2086
  [ -z "$env_req$env_opt" ] || say "        it lacks: $(names $env_req $env_opt)"
else
  [ ${#ENV_CANDIDATES[@]} -gt 0 ] || say "  none named (ADCP_ENV is not set, and no folder was given)"
  say "  nothing usable found"
fi

# --- AmberTools -----------------------------------------------------------------
# In the environment is simplest. Otherwise a cluster's own: whatever AMBERHOME
# or PATH points at, which is where `module load AmberTools` puts it. That one
# is used only once it is complete and its MMPBSA.py runs here, and it is then
# reached through small wrapper scripts that carry the settings the module
# made (its PYTHONPATH and libraries), so it works later without the module
# loaded, and those settings never reach adcp-dock's own Python.
AMB_DIR=""; AMB_HOME=""; AMB_HOW=""; AMB_WHY=""; AMB_WRAPPED=0
amber_external() {
  local d b ok cands=() limit=""
  command -v timeout >/dev/null 2>&1 && limit="timeout 180"
  [ -n "${AMBERHOME:-}" ] && cands+=("$AMBERHOME/bin")
  if command -v tleap >/dev/null 2>&1; then cands+=("$(dirname "$(command -v tleap)")"); fi
  for d in ${cands[@]+"${cands[@]}"}; do
    [ -n "$found_env" ] && [ "$d" -ef "$found_env/bin" ] && continue
    ok=1
    for b in "${AMBER_BINS[@]}"; do [ -x "$d/$b" ] || ok=0; done
    if [ "$ok" = 0 ]; then
      AMB_WHY="$d has only part of AmberTools (it needs ${AMBER_BINS[*]})"; continue
    fi
    if grep -q "adcp-dock amber wrapper" "$d/tleap" 2>/dev/null; then
      AMB_DIR="$d"; AMB_WRAPPED=1
      AMB_HOME=$(sed -n "s/^AMBERHOME='\(.*\)'\$/\1/p" "$d/tleap" | head -n 1 || true)
      AMB_HOW="the wrappers an earlier run wrote, for $AMB_HOME"
      return 0
    fi
    if ! (cd / && $limit "$d/MMPBSA.py" --help >/dev/null 2>&1); then
      AMB_WHY="$d/MMPBSA.py does not run here, so that AmberTools is not used"; continue
    fi
    if command -v ldd >/dev/null 2>&1 && ldd "$d/sander" 2>/dev/null | grep -q "not found"; then
      AMB_WHY="$d/sander cannot find all its libraries, so that AmberTools is not used"; continue
    fi
    AMB_DIR="$d"
    AMB_HOME="$(dirname "$d")"
    if [ -n "${AMBERHOME:-}" ] && [ "$AMBERHOME/bin" -ef "$d" ]; then AMB_HOME="$AMBERHOME"; fi
    local mods
    mods=$(printf '%s' "${LOADEDMODULES:-}" | tr ':' '\n' | grep -i amber | tr '\n' ' ' || true)
    if [ -n "$mods" ]; then
      AMB_HOW="the cluster's own, from module ${mods% }"
    else
      AMB_HOW="the one on PATH"
    fi
    return 0
  done
  return 0
}

head_ "Looking for AmberTools (for MM-GBSA)"
AMBER_SRC=none        # env, external, add (to the environment), build (in a new one) or none
if [ -n "$found_env" ] && [ "$env_amber" = 1 ]; then
  AMBER_SRC=env
  say "  USING the AmberTools in $found_env"
else
  amber_external
  if [ -n "$AMB_DIR" ]; then
    AMBER_SRC=external
    say "  USING $AMB_HOW"
    say "        $AMB_DIR"
  else
    [ -z "$AMB_WHY" ] || wrap "$AMB_WHY"
    say "  not found in the environment, on PATH or in a loaded module"
  fi
fi

# --- ADFRsuite ------------------------------------------------------------------
ADFR_CANDIDATES=()
for d in "${ADFRSUITE_HOME:-}" "$PREV_ADFR" "${PREFIX:+$PREFIX/tools/docking/ADFRsuite}" \
         "$HERE/tools/docking/ADFRsuite" \
         "${ADCP_SHARED:+$ADCP_SHARED/tools/docking/ADFRsuite}" \
         "${ADCP_SHARED:+$ADCP_SHARED/third_party/ADFRsuite}"; do
  [ -n "$d" ] || continue
  case " ${ADFR_CANDIDATES[*]-} " in *" $d "*) continue ;; esac
  ADFR_CANDIDATES+=("$d")
done
found_adfr=""
head_ "Looking for ADFRsuite"
for d in ${ADFR_CANDIDATES[@]+"${ADFR_CANDIDATES[@]}"}; do
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

# --- decide -------------------------------------------------------------------
# What the environment found lacks, and whether it can be added there.
need_env=0; need_adfr=0; need_suite=0; need_adapt=0
TOPUP=()        # what to add to the environment found
CANNOT=()       # what it lacks and will go on lacking
if [ -n "$found_env" ]; then
  # shellcheck disable=SC2206
  lack=($env_req $env_opt)
  [ "$AMBER_SRC" = none ] && lack+=(amber)
  for m in ${lack[@]+"${lack[@]}"}; do
    if [ "$env_writable" = 1 ] && { [ "$env_kind" = conda ] || [ -n "$(pip_spec "$m")" ]; }; then
      TOPUP+=("$m")
    else
      CANNOT+=("$m")
    fi
  done
  # Without numpy, RDKit or pdb2pqr nothing runs. An environment that lacks
  # one and cannot be given it is not used; a new one is built instead.
  for m in ${CANNOT[@]+"${CANNOT[@]}"}; do
    case " $env_req " in
      *" $m "*)
        wrap "$found_env lacks $(name_of "$m"), which adcp-dock cannot run without, and this cannot add it there$([ "$env_writable" = 1 ] || echo ' (you cannot write to it)'). It is not used; a new environment is built instead."
        found_env=""; TOPUP=(); CANNOT=()
        [ "$AMBER_SRC" = env ] && AMBER_SRC=none
        break ;;
    esac
  done
fi
if [ -z "$found_env" ]; then
  need_env=1
  [ "$AMBER_SRC" = external ] || AMBER_SRC=build
fi
for m in ${TOPUP[@]+"${TOPUP[@]}"}; do [ "$m" = amber ] && AMBER_SRC=add; done

head_ "What this machine still needs"
if [ "$need_env" = 1 ]; then
  if [ -n "$PREFIX" ]; then
    say "  a Python environment, built in $PREFIX/envs/adcp"
  else
    say "  a Python environment; give a folder to build it in: ./install.sh <folder>"
  fi
  if [ "$AMBER_SRC" = external ]; then
    say "        without AmberTools, since the cluster's own is used"
  else
    say "        with AmberTools"
  fi
fi
if [ ${#TOPUP[@]} -gt 0 ]; then
  say "  $(names "${TOPUP[@]}"), added to $found_env"
  say "        beside what is there: nothing already installed is changed"
fi
for m in ${CANNOT[@]+"${CANNOT[@]}"}; do
  if [ "$env_writable" = 1 ]; then
    why="which pip cannot add it to"
  else
    why="which you cannot write to"
  fi
  if [ "$m" = amber ]; then
    wrap "AmberTools is missing from $found_env, $why, so $(what_for amber) is off. To turn it on, load your cluster's AmberTools module (module spider AmberTools names it) and run this again; or run this with a folder of your own (./install.sh <folder>) to build an environment that has it."
  else
    wrap "$(name_of "$m") is missing from $found_env, $why, so $(what_for "$m") is off. Whoever looks after that environment can add $(conda_spec "$m"); or run this with a folder of your own (./install.sh <folder>) to build one that has it."
  fi
done
[ -n "$found_adfr" ] || { need_adfr=1; say "  ADFRsuite (adcp, agfr, prepare_receptor), downloaded from Scripps"; }
if [ -z "$found_suite" ]; then
  if [ "$NO_ADCP11" = 1 ]; then
    say "  ADCP 1.1 is not here, and --no-adcp11 leaves it out. ADapt, the"
    say "        non-standard and D residues and AGFR 1.1 need it."
  else
    need_suite=1
    wrap "ADCP 1.1, into $TOOLS_AT/adcpsuite (a few minutes and 3 GB; --no-adcp11 leaves it out)"
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
if [ "$need_env$need_adfr$need_suite$need_adapt${#TOPUP[@]}" = "00000" ]; then
  if [ ${#CANNOT[@]} -eq 0 ]; then
    say "  nothing. Everything it needs is already on this machine."
  else
    say "  nothing it can add. Everything else it needs is already on this machine."
  fi
fi

# --- room -----------------------------------------------------------------------
# How many more files and kilobytes the folder at $1 may hold, as "FILES KB",
# with "-" for a number this cannot tell. The tighter of the filesystem's free
# space and any Lustre quota (user, group or project) on it. Best effort: a
# quota this cannot read is not a reason to stop.
room_left() {
  local d=$1 files="" kb="" t mnt line f k id kind who
  while [ ! -e "$d" ] && [ "$d" != / ]; do d=$(dirname "$d"); done
  # The filesystem itself.
  line=$(df -P -k "$d" 2>/dev/null | awk 'NR==2 {print $4}' || true)
  case "$line" in ''|*[!0-9]*) ;; *) kb=$line ;; esac
  line=$(df -P -i "$d" 2>/dev/null | awk 'NR==2 {print $4}' || true)
  case "$line" in ''|*[!0-9]*) ;; *) files=$line ;; esac
  # Lustre quotas, which df does not see.
  t=$(stat -f -c %T "$d" 2>/dev/null || true)
  if [ "$t" = lustre ] && command -v lfs >/dev/null 2>&1; then
    mnt=$(df -P "$d" 2>/dev/null | awk 'NR==2 {print $6}' || true)
    id=$(lsattr -pd "$d" 2>/dev/null | awk '{print $1}' || true)
    for kind in u g p; do
      case $kind in
        u) who=$(id -un 2>/dev/null || true) ;;
        g) who=$(stat -c %G "$d" 2>/dev/null || true) ;;
        p) who=$id; case "$who" in ''|0|*[!0-9]*) continue ;; esac ;;
      esac
      [ -n "$who" ] && [ -n "$mnt" ] || continue
      # -q: no headings. The numbers follow the filesystem's name, on the
      # same line or the next: used, quota, limit, grace for space in KB,
      # then the same for files. A "*" marks one that is over.
      # shellcheck disable=SC2046
      set -- $(lfs quota -q -$kind "$who" "$mnt" 2>/dev/null | tr -d '*' | tr -s ' \n' '  ' || true)
      [ $# -ge 9 ] || continue
      shift
      local kused=$1 kq=$2 klim=$3 fused=$5 fq=$6 flim=$7
      [ "$klim" != 0 ] 2>/dev/null || klim=$kq
      [ "$flim" != 0 ] 2>/dev/null || flim=$fq
      case "$klim$kused" in *[!0-9]*) ;; *)
        if [ "$klim" -gt 0 ]; then k=$((klim - kused)); [ $k -ge 0 ] || k=0
          if [ -z "$kb" ] || [ "$k" -lt "$kb" ]; then kb=$k; fi; fi ;;
      esac
      case "$flim$fused" in *[!0-9]*) ;; *)
        if [ "$flim" -gt 0 ]; then f=$((flim - fused)); [ $f -ge 0 ] || f=0
          if [ -z "$files" ] || [ "$f" -lt "$files" ]; then files=$f; fi; fi ;;
      esac
    done
  fi
  printf '%s %s' "${files:--}" "${kb:--}"
}
thousands() { printf '%s' "$1" | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta'; }

need_files=0; need_gb=0; keep_files=0; keep_gb=0
add_room() {  # PART: adds what that part needs while installing, and keeps
  local f g kf kg
  f=ROOM_$1_FILES; g=ROOM_$1_GB; kf=KEEP_$1_FILES; kg=KEEP_$1_GB
  need_files=$((need_files + ${!f})); need_gb=$((need_gb + ${!g}))
  keep_files=$((keep_files + ${!kf})); keep_gb=$((keep_gb + ${!kg}))
}
[ "$need_env" = 1 ] && add_room ENV
[ ${#TOPUP[@]} -gt 0 ] && add_room TOPUP
[ "$need_adfr" = 1 ] && add_room ADFR
[ "$need_suite" = 1 ] && add_room SUITE
ROOM_OK=1
if [ "$need_files" -gt 0 ]; then
  where="${PREFIX:-$HERE}"
  head_ "Room in $where"
  read -r have_files have_kb <<<"$(room_left "$where")"
  say "  needs about $(thousands $need_files) files and $need_gb GB while it installs;"
  say "  about $(thousands $keep_files) files and $keep_gb GB stay once its download cache is removed"
  if [ "$have_files" = - ] && [ "$have_kb" = - ]; then
    wrap "How much room you have there could not be read. Your cluster's quota command (often showquota, or quota -s) shows it; most clusters' documentation names it."
  else
    say "  you have $([ "$have_files" = - ] && echo 'an unknown number of' || thousands "$have_files") files and $([ "$have_kb" = - ] && echo 'an unknown amount' || echo "$((have_kb / 1048576)) GB") left there"
    if { [ "$have_files" != - ] && [ "$have_files" -lt "$need_files" ]; } || \
       { [ "$have_kb" != - ] && [ "$have_kb" -lt $((need_gb * 1048576)) ]; }; then
      ROOM_OK=0
      say ""
      say "  NOT ENOUGH ROOM."
      wrap "An install that runs out part way leaves a half built environment, and while your allowance is full your other jobs cannot write their files either. Choose a folder with more room, often a project or group space, and give that instead; or free some room (your cluster's quota command shows what you use), or ask your cluster's support for a larger file allowance. If you are sure there is room and this reading is wrong, ADCP_IGNORE_ROOM=1 ./install.sh ... goes ahead anyway."
    fi
  fi
fi

if [ "$CHECK_ONLY" = 1 ]; then
  say ""
  say "--check: nothing was changed."
  exit 0
fi

if [ "$need_env" = 1 ] && [ -z "$PREFIX" ]; then
  example='/path/with/room/adcp'
  [ -n "${SCRATCH:-}" ] && example="$SCRATCH/adcp"
  cat >&2 <<MSG

No Python environment with what adcp-dock needs was found, and no folder
was given to build one in. Run it again with a folder that has room:

    ./install.sh $example

Use your scratch or project space, not your home directory: most clusters
allow only a few thousand files in home, and this needs about 150,000
while it installs.
envs/adcp is added to the folder you give.

If a suitable environment already exists somewhere this did not look,
point at it instead and nothing is built:

    ADCP_ENV=/path/to/env ./install.sh

MSG
  exit 1
fi
if [ "$ROOM_OK" = 0 ] && [ "${ADCP_IGNORE_ROOM:-0}" != 1 ]; then
  say ""
  say "Stopped before changing anything."
  exit 1
fi

# --- where downloads, caches and helper programs go ----------------------------
# In the folder given, or in this checkout's tools folder (ignored by git),
# never in the home directory.
TOOLS_DIR="${PREFIX:+$PREFIX/tools}"; TOOLS_DIR="${TOOLS_DIR:-$HERE/tools}"
CACHE_BASE="${PREFIX:+$PREFIX/conda_shared}"; CACHE_BASE="${CACHE_BASE:-$HERE/tools/conda_shared}"
# Nothing is written into the home directory: pip keeps no cache (its own
# builds in a subprocess included, which do not see --no-cache-dir), and
# anything else that keeps one keeps it here. ADCP 1.1's installer and every
# pip run below inherit these.
export PIP_NO_CACHE_DIR=1 PIP_DISABLE_PIP_VERSION_CHECK=1
export XDG_CACHE_HOME="$CACHE_BASE/cache" MPLCONFIGDIR="$CACHE_BASE/cache/matplotlib"

download() {  # download URL FILE
  if command -v curl >/dev/null 2>&1; then
    curl -fSL --retry 3 --connect-timeout 30 -o "$2" "$1"
  else
    wget -q -O "$2" "$1"
  fi
}

# Whatever makes conda environments here: micromamba or mamba when one is on
# PATH, or the cluster's own conda when it is recent enough to solve quickly
# (23.10 or newer); otherwise micromamba, a single static file, fetched into
# the tools folder once. So a cluster with no conda, or an old one, needs
# nothing done first.
SOLVER=""
newer_or_same() { [ "$(printf '%s\n%s\n' "$2" "$1" | sort -V | head -n 1)" = "$2" ]; }
find_solver() {
  local c p v
  [ -z "$SOLVER" ] || return 0
  for c in "$TOOLS_DIR/micromamba" micromamba mamba "${CONDA_EXE:-}" conda; do
    [ -n "$c" ] || continue
    p=$(command -v "$c" 2>/dev/null) || continue
    [ -x "$p" ] || continue
    if [ "$(basename "$p")" = conda ]; then
      v=$(py_clean "$p" --version 2>/dev/null | awk '{print $2}' || true)
      if [ -z "$v" ] || ! newer_or_same "$v" 23.10; then
        say "  $p is conda ${v:-of an unknown version}, which solves slowly; micromamba is used instead"
        continue
      fi
    fi
    SOLVER=$p
    return 0
  done
  local mm url
  case "$(uname -m)" in
    x86_64) mm=linux-64 ;; aarch64) mm=linux-aarch64 ;; ppc64le) mm=linux-ppc64le ;;
    *) echo "  there is no conda here, and micromamba is not built for $(uname -m)" >&2; exit 1 ;;
  esac
  STEP="downloading micromamba"
  say "  no recent conda on PATH; fetching micromamba into $TOOLS_DIR"
  mkdir -p "$TOOLS_DIR"
  url="https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-$mm"
  if ! download "$url" "$TOOLS_DIR/micromamba.part"; then
    rm -f "$TOOLS_DIR/micromamba.part"
    wrap "Could not download micromamba from $url. Run this on a login node, which can reach the internet, and run it again." >&2
    exit 1
  fi
  mv "$TOOLS_DIR/micromamba.part" "$TOOLS_DIR/micromamba"
  chmod 775 "$TOOLS_DIR/micromamba"
  SOLVER="$TOOLS_DIR/micromamba"
}

# The solver, with its package cache in the folder rather than ~/.conda, and
# HOME pointing beside it for these commands only: conda keeps a list of every
# environment it has made in ~/.conda/environments.txt, and nothing about
# this install belongs there. A cache already named by whoever ran this is
# used rather than a second one beside it, and is left as it is afterwards.
OWN_CACHE=0
[ -n "${CONDA_PKGS_DIRS:-}" ] || OWN_CACHE=1
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$CACHE_BASE/pkgs}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$CACHE_BASE/mamba}"
# A login node is shared, and many clusters allow each person only a few
# cores across all of them; micromamba would otherwise use a thread for every
# core it can see. ADCP 1.1's installer inherits the same limits.
export MAMBA_DOWNLOAD_THREADS="${MAMBA_DOWNLOAD_THREADS:-4}"
export MAMBA_EXTRACT_THREADS="${MAMBA_EXTRACT_THREADS:-4}"
solver() {
  mkdir -p "$CONDA_PKGS_DIRS" "$CACHE_BASE/home"
  py_clean HOME="$CACHE_BASE/home" "$SOLVER" "$@"
}
# The packages fetched to build an environment are not needed once it is
# built, and they are tens of thousands of files: 50,579 on a clean install.
# Only the cache this script made is removed, never one a cluster or a person
# set up, and only when the environment holds no symbolic link into it. The
# environment's files are hard links to the cache's, or copies, so they stay
# whole when the cache's names for them go. A clean command is not used:
# conda's keeps every package an environment still links to, and both read
# every cache their configuration names, which can include someone else's.
remove_own_cache() {  # $1 = the environment just built or added to
  local pk="$CONDA_PKGS_DIRS" n
  [ "$OWN_CACHE" = 1 ] && [ -d "$pk" ] || return 0
  case "$pk" in "$CACHE_BASE"/*) ;; *) return 0 ;; esac
  if [ -n "$(find "$1" -type l -lname "$pk/*" -print -quit 2>/dev/null)" ]; then
    say "  the download cache is kept: $1 links into it"
    return 0
  fi
  n=$(find "$pk" 2>/dev/null | wc -l)
  rm -rf "$pk"
  say "  removed the download cache it was built from ($(thousands "$n") files)"
}

# --- add what the environment found lacks ---------------------------------------
# Beside what is there: --freeze-installed keeps every installed package at
# the version it has, and pip is held to the same with a constraints file, so
# nothing another user of a shared environment relies on moves. If that
# cannot be done, nothing is changed and this says so.
# All together first, which is one solve; when that cannot be done, one at a
# time, so that one package that cannot be added does not keep out the rest.
add_conda() {  # names...
  local specs=() m
  for m in "$@"; do specs+=("$(conda_spec "$m")"); done
  solver install -y -p "$found_env" --freeze-installed -c conda-forge --override-channels \
    "${specs[@]}"
}
add_pip() {  # names...
  local specs=() m cons rc=0
  for m in "$@"; do [ -z "$(pip_spec "$m")" ] || specs+=("$(pip_spec "$m")"); done
  [ ${#specs[@]} -gt 0 ] || return 1
  cons=$(mktemp "${TMPDIR:-/tmp}/adcp-constraints.XXXXXX")
  py_clean "$found_env/bin/python" -m pip freeze --disable-pip-version-check 2>/dev/null \
    | grep -v -e '^-e ' -e ' @ ' -e '^#' > "$cons" || true
  (cd / && py_clean "$found_env/bin/python" -m pip install --no-cache-dir --no-input \
     --disable-pip-version-check -c "$cons" "${specs[@]}") || rc=$?
  rm -f "$cons"
  return $rc
}
if [ ${#TOPUP[@]} -gt 0 ]; then
  head_ "Adding $(names "${TOPUP[@]}") to $found_env"
  T1=$(date +%s)
  how=add_pip
  if [ "$env_kind" = conda ]; then
    find_solver
    say "  using $SOLVER"
    how=add_conda
  else
    say "  using pip, holding every package already there at its version"
  fi
  STEP="adding $(names "${TOPUP[@]}") to $found_env"
  failed=()
  if ! "$how" "${TOPUP[@]}"; then
    if [ ${#TOPUP[@]} -gt 1 ]; then
      say ""
      say "  Not all at once; trying them one at a time."
      for m in "${TOPUP[@]}"; do
        say "  -> $(name_of "$m")"
        "$how" "$m" || failed+=("$m")
      done
    else
      failed=("${TOPUP[@]}")
    fi
  fi
  [ "$how" = add_conda ] && remove_own_cache "$found_env"
  STEP="checking $found_env"
  env_inspect "$found_env"
  env_req=$E_REQ; env_opt=$E_OPT; env_amber=$E_AMBER
  if [ ${#failed[@]} -eq 0 ]; then
    say "  added, in $(mins "$T1")"
  else
    wrap "$(names "${failed[@]}") could not be added without changing what is already in $found_env, so nothing was changed for it. The lines above say what stood in the way."
  fi
  [ "$AMBER_SRC" = add ] && { [ "$env_amber" = 1 ] && AMBER_SRC=env || AMBER_SRC=none; }
  if [ -n "$env_req" ]; then
    # shellcheck disable=SC2086
    cat >&2 <<MSG

  $found_env still lacks $(names $env_req), and adcp-dock cannot run
  without it. Either build an environment of its own, which leaves that
  one alone:

      ./install.sh /path/with/room/adcp

  or add the missing packages to that environment yourself and run this
  again.

MSG
    exit 1
  fi
fi

# --- build an environment, when none was usable ----------------------------------
BUILT_HERE=0
if [ "$need_env" = 1 ]; then
  ENVDIR="$PREFIX/envs/adcp"
  head_ "Building $ENVDIR"
  T1=$(date +%s)
  find_solver
  say "  using $SOLVER"
  SPECS=(python=3.11 "numpy<2" rdkit pdb2pqr pdbfixer openmm biopython)
  [ "$AMBER_SRC" = external ] || SPECS+=(ambertools)
  say "  this is the slow part: a few minutes on a fast network, most of it AmberTools"
  STEP="building the environment in $ENVDIR"
  if [ -d "$ENVDIR/conda-meta" ]; then
    # An environment an earlier run began and did not finish: completed in
    # place rather than made again, so nothing already there is fetched twice.
    say "  $ENVDIR was begun by an earlier run; completing it"
    solver install -y -p "$ENVDIR" -c conda-forge --override-channels "${SPECS[@]}"
  else
    solver create -y -p "$ENVDIR" -c conda-forge --override-channels "${SPECS[@]}"
  fi
  remove_own_cache "$ENVDIR"
  echo "Built by install.sh in $HERE on $(date +%Y-%m-%d)" > "$ENVDIR/.adcp-dock-env"
  found_env="$ENVDIR"
  env_inspect "$found_env"
  env_req=$E_REQ; env_opt=$E_OPT; env_amber=$E_AMBER; env_kind=$E_KIND; env_writable=$E_WRITABLE
  [ "$AMBER_SRC" = build ] && { [ "$env_amber" = 1 ] && AMBER_SRC=env || AMBER_SRC=none; }
  say "  built in $(mins "$T1")"
fi
[ -f "$found_env/.adcp-dock-env" ] && BUILT_HERE=1
# One an earlier version of this script built carries no mark, but has this
# checkout installed in it already. Asked from /, so that the checkout itself
# is not what answers.
if [ "$BUILT_HERE" = 0 ]; then
  at=$(cd / && py_clean "$found_env/bin/python" -c \
       "import adcp_dock, os; print(os.path.dirname(os.path.dirname(os.path.abspath(adcp_dock.__file__))))" \
       2>/dev/null || true)
  [ -n "$at" ] && [ "$at" -ef "$HERE" ] && BUILT_HERE=1
fi

# --- AmberTools from outside the environment: wrappers ---------------------------
AMB_ROOT="$TOOLS_DIR/amber"
sq() { printf "'%s'" "$(printf '%s' "$1" | sed "s/'/'\\\\''/g")"; }
# One for every program in that AmberTools' bin folder, not only the four
# checked for: the rescoring also runs ambpdb, and finds each program as
# $AMBERHOME/bin/<name>. Written into a new folder that then replaces the old
# one, so a program the module no longer has does not linger.
write_amber_wrappers() {
  local b f p new="$AMB_ROOT/bin.part"
  rm -rf "$new"
  mkdir -p "$new"
  for p in "$AMB_DIR"/*; do
    [ -f "$p" ] && [ -x "$p" ] || continue
    b=$(basename "$p")
    f="$new/$b"
    {
      echo '#!/bin/sh'
      echo "# adcp-dock amber wrapper, written by install.sh on $(date +%Y-%m-%d)."
      echo "# Runs $b from the AmberTools in $AMB_HOME with the settings it had"
      echo "# when install.sh found it (a cluster module makes them), so it works"
      echo "# without that module loaded and they never reach adcp-dock's own"
      echo "# Python. Run ./install.sh again to write it afresh."
      echo "AMBERHOME=$(sq "$AMB_HOME")"
      echo "PATH=$(sq "$AMB_DIR:$PATH")"':"$PATH"'
      if [ -n "${LD_LIBRARY_PATH:-}" ]; then
        echo "LD_LIBRARY_PATH=$(sq "$LD_LIBRARY_PATH")"'${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}'
        echo "export LD_LIBRARY_PATH"
      fi
      if [ -n "${PERL5LIB:-}" ]; then echo "PERL5LIB=$(sq "$PERL5LIB"); export PERL5LIB"; fi
      if [ -n "${PYTHONPATH:-}" ]; then
        echo "PYTHONPATH=$(sq "$PYTHONPATH"); export PYTHONPATH"
      else
        echo "unset PYTHONPATH"
      fi
      echo "unset PYTHONHOME"
      echo "export AMBERHOME PATH"
      echo "exec $(sq "$AMB_DIR/$b")"' "$@"'
    } > "$f"
    chmod 755 "$f"
  done
  rm -rf "$AMB_ROOT/bin"
  mv "$new" "$AMB_ROOT/bin"
}
if [ "$AMBER_SRC" = external ] && [ "$AMB_WRAPPED" = 0 ]; then
  STEP="writing the AmberTools wrappers"
  write_amber_wrappers
  head_ "AmberTools"
  wrap "The cluster's own AmberTools in $AMB_HOME is used through $(ls "$AMB_ROOT/bin" | wc -l) small scripts in $AMB_ROOT/bin, one for each of its programs, which carry the settings its module made. Nothing needs loading later; if that module is updated or removed, run ./install.sh again."
fi
if [ "$AMBER_SRC" = external ] && [ "$AMB_WRAPPED" = 1 ]; then AMB_ROOT="$(dirname "$AMB_DIR")"; fi

# --- the one file that sets everything up --------------------------------------
# Written as soon as the environment is known, and again at the end. The
# launcher sources it to find a Python new enough to run this package, so it
# must exist even if installing ADFRsuite fails partway.
write_activate() {
  # $1 = the environment, $2 = ADFRsuite (may be empty).
  #
  # Written line by line rather than with ${2:+...} inside a heredoc, whose
  # quotes would need escaping twice (for the heredoc and for the shell that
  # later sources the file); a wrong escape puts literal quote marks into
  # the variable's value.
  local env_dir="$1" adfr="$2"
  {
    echo "# Written by install.sh on $(date +%Y-%m-%d). Source it, do not run it:"
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
      echo "export PATH=\"$env_dir/bin:\$PATH\""
    fi
    case "$AMBER_SRC" in
      env)
        echo "export AMBERHOME=\"$env_dir\"" ;;
      external)
        echo "# AmberTools: the cluster's own, in $AMB_HOME, through wrappers"
        echo "# that carry its module's settings."
        echo "export AMBERHOME=\"$AMB_ROOT\""
        echo "export PATH=\"$AMB_ROOT/bin:\$PATH\"" ;;
      *)
        echo "# No AmberTools was found, so MM-GBSA rescoring is off. Load your"
        echo "# cluster's AmberTools module and run ./install.sh again to add it." ;;
    esac
    if [ -n "$adfr" ]; then
      echo "export ADFRSUITE_HOME=\"$adfr\""
      echo "export PATH=\"$adfr/bin:\$PATH\""
    fi
    if [ -x "$HERE/tools/bin/adcp-dock" ]; then
      echo "export PATH=\"$HERE/tools/bin:\$PATH\""
    fi
    echo "# Replaced rather than added to: a module's Python packages on it"
    echo "# (AmberTools adds some) stop this environment's Python importing numpy."
    echo "export PYTHONPATH=\"\$ADCP_HOME\""
  } > "$ACT.part"
  mv "$ACT.part" "$ACT"
}
STEP="writing activate.sh"
write_activate "$found_env" "$found_adfr"
say ""
say "Wrote $ACT (the launcher reads this to find the right Python)"

# --- ADFRsuite ---------------------------------------------------------------------
if [ "$need_adfr" = 1 ]; then
  TOOLS="${PREFIX:-$HERE}/tools/docking"
  head_ "Installing ADFRsuite into $TOOLS"
  T1=$(date +%s)
  STEP="installing ADFRsuite"
  mkdir -p "$TOOLS"
  cd "$TOOLS"
  URL=https://ccsb.scripps.edu/adfr/download/1038/
  TAR=ADFRsuite_x86_64Linux_1.0.tar.gz
  if [ ! -s "$TAR" ]; then
    say "  downloading from Scripps (about 100 MB)"
    STEP="downloading ADFRsuite"
    if ! download "$URL" "$TAR.part"; then
      rm -f "$TAR.part"
      cat >&2 <<MSG
  Could not download ADFRsuite from $URL

  This is a login-node step: compute nodes usually have no outbound
  network. If you are on one, that is the whole explanation.

  Otherwise fetch it by hand from https://ccsb.scripps.edu/adfr/downloads/,
  put ADFRsuite_x86_64Linux_1.0.tar.gz in $TOOLS,
  and run this script again; it uses the file that is already there.
MSG
      exit 1
    fi
    mv "$TAR.part" "$TAR"
  else
    say "  $TAR is already here, using it"
  fi

  # The same pinned checksum as adcp_dock/adfr.py: the archive is run as an
  # installer, so a truncated or substituted download must stop here.
  SHA=66539864ffb7ba87728d52ee50dde78586f551c7ca96aa455dae761875faf7e4
  got=$(sha256sum "$TAR" 2>/dev/null | cut -d' ' -f1 || true)
  [ -n "$got" ] || got=$(shasum -a 256 "$TAR" | cut -d' ' -f1)
  if [ "$got" != "$SHA" ]; then
    echo "  $TOOLS/$TAR does not match its published checksum (got $got)." >&2
    echo "  It is probably an incomplete download. Delete it and run this script again:" >&2
    echo "      rm $TOOLS/$TAR" >&2
    exit 1
  fi
  say "  checksum verified"

  # What the archive unpacks into, read from the archive rather than assumed:
  # a wrong guess `cd`s into a directory that is not there and, under
  # `set -e`, ends the script with no message, immediately after a 100 MB
  # download.
  say "  unpacking"
  STEP="unpacking ADFRsuite"
  # tar is still listing when head has its line and closes the pipe, so tar
  # dies of SIGPIPE; under pipefail that would end the script here, so its
  # status is not the pipeline's.
  SRC=$( { tar tzf "$TAR" 2>/dev/null || true; } | head -1 | cut -d/ -f1)
  [ -n "$SRC" ] || { echo "  $TAR is not a readable tar.gz; delete it (rm $TOOLS/$TAR) and run this again" >&2; exit 1; }
  tar xzf "$TAR"
  [ -d "$SRC" ] || { echo "  $TAR did not unpack into $SRC; delete it (rm $TOOLS/$TAR) and run this again" >&2; exit 1; }
  cd "$SRC"

  say "  running the Scripps installer (a few minutes)"
  STEP="running the ADFRsuite installer"
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
  it and believe it. Remove it, then run this script again:

      rm -rf $TOOLS/ADFRsuite

MSG
    exit 1
  fi
  found_adfr="$TOOLS/ADFRsuite"
  say "  installed at $found_adfr, in $(mins "$T1")"
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

# --- this package, and the adcp-dock command -------------------------------------
# Into an environment this script built, with pip and --no-deps. One it did
# not build is never written to for this: pip would point every user of a
# shared environment at this checkout. The command then comes from a small
# script in this checkout's tools folder, which runs the same thing.
head_ "The adcp-dock command"
STEP="installing adcp-dock"
pip_ok=0
if [ "$BUILT_HERE" = 1 ] && [ "$env_writable" = 1 ]; then
  if py_clean "$found_env/bin/python" -m pip install --quiet --no-deps --no-cache-dir \
       --disable-pip-version-check -e "$HERE"; then
    pip_ok=1
    rm -f "$HERE/tools/bin/adcp-dock"
    say "  installed into $found_env, which this script built"
  fi
fi
if [ "$pip_ok" = 0 ]; then
  mkdir -p "$HERE/tools/bin"
  {
    echo '#!/bin/sh'
    echo "# Written by install.sh: adcp-dock from this checkout, with the"
    echo "# environment's Python. Run ./install.sh again to write it afresh."
    echo "PYTHONPATH=$(sq "$HERE"); export PYTHONPATH"
    echo "unset PYTHONHOME"
    echo "exec $(sq "$found_env/bin/python3") -m adcp_dock.cli"' "$@"'
  } > "$HERE/tools/bin/adcp-dock.part"
  chmod 755 "$HERE/tools/bin/adcp-dock.part"
  mv "$HERE/tools/bin/adcp-dock.part" "$HERE/tools/bin/adcp-dock"
  wrap "$HERE/tools/bin/adcp-dock runs this checkout with $found_env/bin/python3. Nothing was written into that environment."
fi

# --- and again, now that ADFRsuite is known -------------------------------
STEP="writing activate.sh"
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
  T1=$(date +%s)
  STEP="installing ADCP 1.1"
  say "  the Scripps Research Institute's recipe: a few minutes and 3 GB"
  # HOME points beside it for this one command, as for conda above: nothing
  # of it belongs in the home directory. A micromamba fetched above is
  # reused rather than fetched again.
  mkdir -p "$TOOLS_AT/.adcpsuite_home"
  if HOME="$TOOLS_AT/.adcpsuite_home" ADFRSUITE_HOME="$found_adfr" ADCP_HOME="$HERE" \
     PATH="$TOOLS_DIR:$PATH" PYTHONPATH="$HERE" \
     "$found_env/bin/python" -u -m adcp_dock.adcp11_install --dest "$SUITE_DEST" 2>&1 \
     | sed 's/^/  /'; then :; fi
  rmdir "$TOOLS_AT/.adcpsuite_home" 2>/dev/null || true
  found_suite=$(suite_found)
  if [ -n "$found_suite" ]; then
    say "  verified: ADCP 1.1 is at $found_suite, in $(mins "$T1")"
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
STEP="fetching ADapt"
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
elif command -v git >/dev/null 2>&1; then
  adapt_clone "$ADAPT_DIR"
  wrap "$ADAPT_LINE"
else
  wrap "not fetched: git is not installed here. ADCP docking works without it."
fi
if [ -z "$found_suite" ] && adapt_is "${ADAPT_DIR:-/nonexistent}"; then
  wrap "ADapt runs on ADCP 1.1, which is not installed, so the interface will not offer it until ADCP 1.1 is added."
fi

# --- prove it ------------------------------------------------------------
head_ "Checking"
STEP="checking the install"
# shellcheck disable=SC1090
source "$ACT"
py_clean PYTHONPATH="$HERE" "$found_env/bin/python" -m adcp_dock.cli check || true

head_ "Summary"
off=()
# shellcheck disable=SC2206
for m in $env_req $env_opt; do off+=("$m"); done
[ "$AMBER_SRC" = env ] || [ "$AMBER_SRC" = external ] || off+=(amber)
say "  Python environment  $found_env"
case "$AMBER_SRC" in
  env) say "  AmberTools          in that environment" ;;
  external) say "  AmberTools          the cluster's own, $AMB_HOME" ;;
  *) say "  AmberTools          not found" ;;
esac
say "  ADFRsuite           $found_adfr"
say "  ADCP 1.1            ${found_suite:-not installed (optional; ADapt needs it)}"
if adapt_is "${ADAPT_DIR:-/nonexistent}"; then
  adapt_at=$(adapt_commit "$ADAPT_DIR")
  if adapt_embeds "$ADAPT_DIR"; then
    say "  ADapt               $ADAPT_DIR${adapt_at:+ ($adapt_at)}"
  else
    say "  ADapt               $ADAPT_DIR${adapt_at:+ ($adapt_at)}, older than the interface needs"
  fi
else
  say "  ADapt               not installed (optional; ./update.sh fetches it once it is public)"
fi
if [ ${#off[@]} -gt 0 ]; then
  say ""
  say "  Turned off, because something is missing:"
  for m in "${off[@]}"; do say "    $(what_for "$m") (needs $(name_of "$m"))"; done
fi
say ""
say "  Finished in $(mins "$T0")."

cat <<EOF

Next, in this terminal or any new one:

    source $ACT
    adcp-dock check                 lists what it found
    adcp-dock serve --port 0        starts the page; INSTALL.md 1.5 says how to open it

Or open it from your own computer with the desktop launcher (INSTALL.md 5).
Never dock on a login node: adcp-dock sends every calculation to Slurm.
EOF
