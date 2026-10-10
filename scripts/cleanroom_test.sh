#!/usr/bin/env bash
# The whole install of GUI to ADCP Docking, as a new user meets it, on any
# Slurm cluster, with a short summary to send back.
#
# By Aadhil Haq
#
#     bash scripts/cleanroom_test.sh <folder>            on a login node
#     bash scripts/cleanroom_test.sh --report <folder>   once its jobs finish
#
# In <folder>, which must not be inside this checkout, it clones the
# repository (from GitHub when it can be read anonymously, otherwise from this
# checkout), then runs install.sh --check and install.sh as a new user would:
# with nothing in the environment but HOME (an empty folder of its own), PATH
# (/usr/bin and /bin, plus wherever sbatch is), USER and TERM. Everything runs
# at the lowest priority, and micromamba is held to 4 download and 4 extract
# threads, for login nodes that allow only a few cores. It then runs
# ./update.sh, and submits a small docking run (receptor preparation with
# reduce and OpenMM, ADCP, MM-GBSA) through adcp-dock run and the test suite
# as one more job, all with the cluster's default Slurm account.
#
# Each step's exit code, time and the files it made, the quota before and
# after, the install's warnings and the job ids go to <folder>/logs/SUMMARY.txt,
# which is short enough to paste into a message. --report adds how the jobs
# ended. It leaves about 112,000 files: delete <folder> afterwards.
#
# Nothing it does touches this checkout, the home directory or anything
# outside <folder>, apart from the jobs it submits.

if [ -z "${BASH_VERSION:-}" ]; then
  echo "Run this with bash:  bash scripts/cleanroom_test.sh <folder>" >&2
  return 2 2>/dev/null || exit 2
fi
if [ "${BASH_SOURCE[0]}" != "$0" ]; then
  echo "Run this script, do not source it:  bash scripts/cleanroom_test.sh <folder>" >&2
  return 2
fi
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
GITHUB=https://github.com/aadhilhaq/ADCP_docking.git

usage() {
  cat <<'EOF'
Usage:
    bash scripts/cleanroom_test.sh <folder>            run the clean-room install
    bash scripts/cleanroom_test.sh --report <folder>   add how its jobs ended

<folder> must be outside this checkout and empty (or one this script began),
with room for about 112,000 files.
EOF
}

MODE=run
CR=""
for a in "$@"; do
  case "$a" in
    --report) MODE=report ;;
    -h|--help) usage; exit 0 ;;
    -*) echo "cleanroom_test.sh: there is no option $a" >&2; usage >&2; exit 2 ;;
    *) [ -z "$CR" ] || { echo "cleanroom_test.sh: give one folder" >&2; exit 2; }; CR="$a" ;;
  esac
done
[ -n "$CR" ] || { usage >&2; exit 2; }
case "$CR" in /*) ;; *) CR="$PWD/$CR" ;; esac
CR="$(realpath -m "$CR")"
case "$CR/" in
  "$REPO"/*)
    echo "cleanroom_test.sh: $CR is inside this checkout ($REPO)." >&2
    echo "Give a folder somewhere else, for example beside it:  $(dirname "$REPO")/cleanroom_$USER" >&2
    exit 2 ;;
esac

L="$CR/logs"
S="$L/SUMMARY.txt"
STATE="$L/state.env"
say() { printf '%s\n' "$*"; }
sum() { printf '%s\n' "$*" >> "$S"; }
files() { find "$@" 2>/dev/null | wc -l; }
# The folder's files, leaving out the logs this script writes itself.
made() { find "$CR" -path "$L" -prune -o -print 2>/dev/null | wc -l; }
cleanup_line() { say "When you have sent SUMMARY.txt, delete the test (about 112,000 files):  rm -rf $CR"; }

# The quota of the filesystem <folder> is on, however this cluster reports it.
quota_lines() {
  local d="$CR" mnt id
  while [ ! -e "$d" ] && [ "$d" != / ]; do d=$(dirname "$d"); done
  if command -v showquota >/dev/null 2>&1; then
    showquota 2>&1 | head -8
  elif command -v lfs >/dev/null 2>&1 && [ "$(stat -f -c %T "$d" 2>/dev/null)" = lustre ]; then
    mnt=$(df -P "$d" | awk 'NR==2 {print $6}')
    lfs quota -h -u "$USER" "$mnt" 2>&1 | tail -n +2 | head -3
    lfs quota -h -g "$(stat -c %G "$d")" "$mnt" 2>&1 | tail -n +2 | head -3
    id=$(lsattr -pd "$d" 2>/dev/null | awk '{print $1}')
    case "$id" in ''|0|*[!0-9]*) ;; *) lfs quota -h -p "$id" "$mnt" 2>&1 | tail -n +2 | head -3 ;; esac
  elif command -v quota >/dev/null 2>&1; then
    quota -s 2>&1 | head -8
  fi
  df -h --output=target,size,used,avail,pcent "$d" 2>/dev/null | tail -1 || df -h "$d" | tail -1
  df --output=target,itotal,iused,iavail,ipcent "$d" 2>/dev/null | tail -1 || df -i "$d" | tail -1
}

# ============================================================================
if [ "$MODE" = report ]; then
  [ -f "$STATE" ] || { echo "No clean-room test in $CR (no $STATE)." >&2; exit 1; }
  # shellcheck disable=SC1090
  . "$STATE"
  {
    echo
    echo "== Report, $(date '+%Y-%m-%d %H:%M')"
    if [ -z "${JOBS:-}" ]; then
      echo "no jobs were submitted (see the steps above)"
    elif command -v sacct >/dev/null 2>&1; then
      echo "jobs:"
      sacct -j "$(echo "$JOBS" | tr ' ' ',')" -X -n -P -o JobID,JobName,State,Elapsed,ExitCode 2>&1 \
        | awk -F'|' '{printf "  %-14s %-26s %-11s %-10s exit %s\n", $1, $2, $3, $4, $5}' | head -20
    else
      echo "jobs (sacct is not available; squeue shows only the unfinished):"
      squeue -j "$(echo "$JOBS" | tr ' ' ',')" -h 2>&1 | head -10
    fi
    if [ -n "${RUN_DIR:-}" ]; then
      if [ -f "$RUN_DIR/results.txt" ]; then
        echo "results.txt: present"
        head -6 "$RUN_DIR/results.txt" | sed 's/^/  /'
        grep -A3 "^Interface and sequence panel" "$RUN_DIR/results.txt" | sed 's/^/  /' | head -4
      else
        echo "results.txt: not there (yet). The last lines of each job log:"
        for f in "$RUN_DIR"/logs/*.out; do
          [ -f "$f" ] || continue
          echo "  $(basename "$f"):"
          tail -3 "$f" | sed 's/^/    /'
        done
      fi
    fi
    gate=$(ls -t "$L"/gate_*.log 2>/dev/null | head -1)
    if [ -n "$gate" ]; then
      echo "test suite ($(basename "$gate")):"
      grep -E "passed|FAILED|### gate exit|SELFTEST|check exit" "$gate" | sed 's/^/  /' | head -12
    else
      echo "test suite: no log yet"
    fi
    echo "update.sh: exit ${UPDATE_RC:-not run}"
  } >> "$S"
  if [ -n "${JOBS:-}" ] && command -v squeue >/dev/null 2>&1 \
     && [ -n "$(squeue -j "$(echo "$JOBS" | tr ' ' ',')" -h 2>/dev/null)" ]; then
    sum "Some jobs have not finished; run --report again when they have."
  fi
  sum "Cleanup:  rm -rf $CR"
  tail -n 40 "$S"
  echo
  say "Send this file back:  $S"
  cleanup_line
  exit 0
fi

# ============================================================================
# --- the folder ----------------------------------------------------------------
if [ -e "$CR" ] && [ -n "$(ls -A "$CR" 2>/dev/null)" ] && [ ! -f "$STATE" ]; then
  echo "$CR already holds other files, and was not begun by this script;" >&2
  echo "they are left alone. Give a folder that does not exist yet, for example" >&2
  echo "    $CR/cleanroom_$USER" >&2
  exit 2
fi
mkdir -p "$CR/home" "$L" || { echo "cannot create $CR" >&2; exit 1; }
cleanup_line
say ""
T0=$(date +%s)

# A new user's environment: nothing but these. sbatch's own folder is added to
# PATH, since some clusters keep it outside /usr/bin.
SBIN=""
command -v sbatch >/dev/null 2>&1 && SBIN=":$(dirname "$(command -v sbatch)")"
LOW="nice -n 19"
command -v ionice >/dev/null 2>&1 && LOW="$LOW ionice -c3"
# git never waits for a password: a private repository fails at once.
clean() {
  env -i HOME="$CR/home" PATH="/usr/bin:/bin$SBIN" USER="$USER" TERM="${TERM:-xterm}" \
      MAMBA_DOWNLOAD_THREADS=4 MAMBA_EXTRACT_THREADS=4 GIT_TERMINAL_PROMPT=0 $LOW "$@"
}

: > "$L/steps.tsv"
# step NAME LOGFILE command...: runs it, and records exit, seconds and files made.
step() {
  local name=$1 log=$2 t n0 n1 rc
  shift 2
  n0=$(made); t=$(date +%s)
  say "[$(date +%T)] $name"
  "$@" > "$log" 2>&1
  rc=$?
  n1=$(made)
  printf '%s\t%s\t%s\t%s\n' "$name" "$rc" "$(( $(date +%s) - t ))" "$(( n1 - n0 ))" >> "$L/steps.tsv"
  say "           exit $rc, $(( $(date +%s) - t )) s, $(( n1 - n0 )) files"
  return $rc
}
save_state() {
  {
    printf 'JOBS=%q\n' "${JOBS:-}"
    printf 'RUN_DIR=%q\n' "${RUN_DIR:-}"
    printf 'UPDATE_RC=%q\n' "${UPDATE_RC:-}"
    printf 'INSTALL_RC=%q\n' "${INSTALL_RC:-}"
  } > "$STATE"
}
JOBS=""; RUN_DIR=""; UPDATE_RC=""; INSTALL_RC=""
save_state

quota_lines > "$L/quota_before.txt" 2>&1

# --- clone ------------------------------------------------------------------------
SOURCE=""
if [ -d "$CR/ADCP_docking/.git" ]; then
  SOURCE="already there (an earlier run of this script)"
else
  if step "clone from GitHub" "$L/01_clone.txt" clean git clone -q "$GITHUB" "$CR/ADCP_docking"; then
    SOURCE="$GITHUB"
  else
    rm -rf "$CR/ADCP_docking"
    say "           GitHub cannot be read anonymously (still private?); cloning this checkout"
    step "clone from this checkout" "$L/01_clone_local.txt" clean \
      git clone -q "file://$REPO" "$CR/ADCP_docking"
    SOURCE="file://$REPO at $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null)"
  fi
fi
[ -d "$CR/ADCP_docking/.git" ] || { say "the clone failed; see $L"; }
cd "$CR/ADCP_docking" 2>/dev/null || cd "$CR"

# --- install ------------------------------------------------------------------------
if [ -x "$CR/ADCP_docking/install.sh" ]; then
  step "install.sh --check" "$L/02_check.txt" clean ./install.sh --check "$CR/prefix"
  step "install.sh" "$L/03_install.txt" clean ./install.sh "$CR/prefix"
  INSTALL_RC=$?
  save_state
fi

if [ "${INSTALL_RC:-1}" = 0 ]; then
  # --- update.sh, from the same account --------------------------------------------
  step "git pull; ./update.sh" "$L/04_update.txt" clean bash -c 'git pull; ./update.sh'
  UPDATE_RC=$?
  save_state

  # --- a small docking run, through the normal path ------------------------------------
  mkdir -p "$CR/runs"
  RUN_DIR="$CR/runs/trial"
  step "adcp-dock new and run" "$L/05_run_submit.txt" clean bash -c "
    source '$CR/ADCP_docking/activate.sh' && cd '$CR/runs' &&
    adcp-dock new trial -r 1YCR --chain A -s A:54,A:62,A:67,A:96,A:99 \
      -p 'p53: ETFSDLWKLLP' -p 'pdi: LTFEHYWAQLTS' --replicas 4 --steps 200000 \
      --poses 5 --omm-minimize 2 --mmgbsa --mmgbsa-poses 2 --cpus 4 --time 01:00:00 &&
    adcp-dock run trial"
  JOBS=$(grep -oE '^(prepare|dock|report): [0-9]+' "$L/05_run_submit.txt" | awk '{print $2}' | tr '\n' ' ')
  save_state

  # --- the test suite, on a compute node -----------------------------------------------
  cat > "$CR/gate.sbatch" <<EOF
#!/bin/bash
#SBATCH -J adcp-cleanroom-gate
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=24G
#SBATCH -t 02:00:00
#SBATCH -o $L/gate_%j.log
# Written by scripts/cleanroom_test.sh: the clean-room install's own check,
# the launcher's selftest and the test suite, in an empty environment.
W=$CR/work_\$SLURM_JOB_ID
mkdir -p \$W/home \$W/tmp
exec env -i HOME=\$W/home USER=\$USER TERM=xterm PATH=/usr/bin:/bin$SBIN TMPDIR=\$W/tmp \\
  SLURM_JOB_ID=\$SLURM_JOB_ID SLURM_CPUS_PER_TASK=\${SLURM_CPUS_PER_TASK:-8} \\
  ADCP_TEST_TOOLS=$CR/prefix/tools/docking PYTHONUNBUFFERED=1 bash -c '
cd $CR/ADCP_docking || exit 1
echo "node \$(hostname)  job \$SLURM_JOB_ID  \$(git log -1 --format=%h)"
source ./activate.sh
adcp-dock check; echo "check exit=\$?"
"\$ADCP_PYTHON" packaging/launcher.py --selftest 2>&1 | tail -2
bash tests/run_all.sh
echo "### gate exit=\$?"
'
EOF
  say "[$(date +%T)] submitting the test suite"
  out=$(cd "$CR" && clean sbatch "$CR/gate.sbatch" 2>&1)
  printf '%s\n' "$out" > "$L/06_gate_submit.txt"
  gid=$(printf '%s\n' "$out" | grep -oE 'Submitted batch job [0-9]+' | awk '{print $4}')
  [ -n "$gid" ] && JOBS="$JOBS$gid"
  if [ -n "$gid" ]; then
    say "           job $gid"
  else
    say "           not submitted: $(printf '%s' "$out" | tail -1)"
  fi
  save_state
fi

quota_lines > "$L/quota_after.txt" 2>&1

# --- SUMMARY.txt ------------------------------------------------------------------------
{
  echo "ADCP_docking clean-room test, $(date '+%Y-%m-%d %H:%M')"
  echo "host $(hostname), $(uname -sm), user $USER"
  echo "folder $CR"
  echo "code $SOURCE ($(git -C "$CR/ADCP_docking" log -1 --format='%h %ad' --date=short 2>/dev/null))"
  if [ "${INSTALL_RC:-}" = 0 ]; then
    echo "RESULT: install.sh finished"
  else
    echo "RESULT: install.sh STOPPED (exit ${INSTALL_RC:-not run}); its reasons are below"
  fi
  echo
  echo "== Steps (exit, seconds, files made)"
  awk -F'\t' '{printf "  %-26s exit %-3s %6s s %8s files\n", $1, $2, $3, $4}' "$L/steps.tsv"
  if [ -d "$CR/prefix" ]; then
    echo "  files now: environment $(files "$CR/prefix/envs/adcp"), ADFRsuite $(files "$CR/prefix/tools/docking/ADFRsuite"), ADCP 1.1 $(files "$CR/prefix/tools/docking/adcpsuite"), caches $(files "$CR/prefix/conda_shared" "$CR/prefix/tools/docking/mamba_root"), whole folder $(files "$CR")"
  fi
  home=$(ls -A "$CR/home" | tr '\n' ' ')
  echo "  written to the empty home: ${home:-nothing}"
  echo "  total $(( $(date +%s) - T0 )) s"
  echo
  echo "== Quota before"
  sed 's/^/  /' "$L/quota_before.txt"
  echo "== Quota after"
  sed 's/^/  /' "$L/quota_after.txt"
  echo
  if [ -f "$L/03_install.txt" ]; then
    echo "== install.sh: room"
    sed -n '/^== Room/,/^$/p' "$L/03_install.txt" | sed -n '2,4p'
    echo "== install.sh: warnings and errors (at most 40 lines)"
    # Leaving out micromamba's list of packages, and the lines that only say
    # a new account has nothing installed yet.
    grep -n -i -E 'warn|error|stopped|not enough|could not|cannot|failed|turned off|not installed|not found|missing|refused' \
      "$L/03_install.txt" | grep -v -E '^[0-9]+:[[:space:]]+[+-] |gpg-error|━' \
      | grep -v -E '^[0-9]+:[[:space:]]+(not found|ADCP 1.1: not found|ADapt: not found|not found in the environment, on PATH or in a loaded module)$' \
      | cut -c1-160 | head -40
    echo "== install.sh: its summary"
    sed -n '/== Summary ==/,/Finished in/p' "$L/03_install.txt" | head -16
    tail -3 "$L/03_install.txt" | sed '/^[[:space:]]*$/d; s/^/  last: /'
  fi
  [ -f "$L/04_update.txt" ] && echo "== update.sh: exit $UPDATE_RC; $(grep -c 'NOT FOUND' "$L/04_update.txt") NOT FOUND lines"
  echo
  echo "== Jobs: ${JOBS:-none}"
  [ -f "$L/05_run_submit.txt" ] && grep -E '^(prepare|dock|report):|rror|refused|denied' "$L/05_run_submit.txt" | head -6 | sed 's/^/  /'
  [ -f "$L/06_gate_submit.txt" ] && sed 's/^/  gate: /' "$L/06_gate_submit.txt" | tail -2
  echo
  [ -n "$JOBS" ] && echo "When the jobs finish:  bash $REPO/scripts/cleanroom_test.sh --report $CR"
  echo "Cleanup:  rm -rf $CR"
} > "$S"

say ""
say "Every step's output is in $L."
if [ -n "$JOBS" ]; then
  say "When the jobs ($JOBS) have finished, add how they ended:"
  say "    bash $REPO/scripts/cleanroom_test.sh --report $CR"
fi
cleanup_line
say "Send this file back:  $S"
