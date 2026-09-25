#!/usr/bin/env bash
# The whole test suite, in one command.
#
# By Aadhil Haq
#
# Run it on a compute node, adding the account and partition your cluster
# requires:
#
#     srun --pty -c 8 -t 30 bash tests/run_all.sh
#
# The checks that matter need ADFRsuite and Amber, so the suite runs on a
# cluster rather than on GitHub's runners. It does not replace a real
# docking run.
#
# Each test is a script: exit 0 passes, 77 skips, and anything else fails
# with its output attached. A skip is reported as a skip, never as a pass.
set -uo pipefail
cd "$(dirname "$0")/.." || exit 1
REPO=$PWD

# Not on a login node: the suite does real work, and login nodes are shared.
# ALLOW_LOGIN=1 is for re-running one failure by hand.
if [ -z "${SLURM_JOB_ID:-}" ] && [ "${ALLOW_LOGIN:-0}" != "1" ]; then
  echo "This is a login node. Ask for a compute node first:"
  echo
  echo "    srun --pty -p cpu -c 8 -t 30 bash tests/run_all.sh"
  echo
  exit 2
fi

PY="${ADCP_PYTHON:-}"
for c in "$PY" "$REPO/../envs/hope/bin/python" "$(command -v python3)"; do
  [ -n "${c:-}" ] && [ -x "${c:-}" ] && { PY="$c"; break; }
done
[ -x "${PY:-}" ] || { echo "no python found"; exit 2; }

export PYTHONPATH="$REPO:${PYTHONPATH:-}"
# Unbuffered, so a suite that hangs still shows which check it hung on.
export PYTHONUNBUFFERED=1
# shellcheck disable=SC1091
[ -f "$REPO/activate.sh" ] && . "$REPO/activate.sh" >/dev/null 2>&1

echo "python  $PY"
echo "commit  $(git -C "$REPO" rev-parse --short HEAD 2>/dev/null)"
echo "node    $(hostname)  job ${SLURM_JOB_ID:-none}"
echo

PASS=0; FAIL=0; SKIP=0; FAILED=(); SKIPPED=()
run(){  # run <label> <command...>
  local label=$1; shift
  printf '%-54s' "$label"
  local log rc; log=$(mktemp)
  "$@" >"$log" 2>&1; rc=$?
  if [ $rc -eq 0 ]; then
    printf 'ok\n'; PASS=$((PASS+1))
  elif [ $rc -eq 77 ]; then
    printf 'skip\n'; SKIP=$((SKIP+1)); SKIPPED+=("$label")
    sed 's/^/    . /' "$log" | tail -3
  else
    printf 'FAILED (%d)\n' "$rc"; FAIL=$((FAIL+1)); FAILED+=("$label")
    sed 's/^/    | /' "$log" | tail -25
  fi
  rm -f "$log"
}

# How many test files this checkout has, so a partial checkout shows in the
# summary.
N=$(ls "$REPO"/tests/test_*.py 2>/dev/null | wc -l)
echo "$N test files in this checkout"
echo
for t in "$REPO"/tests/test_*.py; do
  run "$(basename "$t")" "$PY" "$t"
done

# Checks that would otherwise exist only in the workflow file, so that the
# cluster runs them too.
run "the launcher reaches the end of --selftest" \
    bash -c '"$1" "$2" --selftest | grep -q "^SELFTEST OK$"' _ "$PY" "$REPO/packaging/launcher.py"

echo
printf '%d passed' "$PASS"
[ "$SKIP" -gt 0 ] && printf ', %d skipped' "$SKIP"
if [ "$FAIL" -eq 0 ]; then
  printf '. All checks passed.\n'
else
  printf ', %d FAILED:\n' "$FAIL"
  printf '  %s\n' "${FAILED[@]}"
  printf '\nFix these here and re-run. Do not push a red suite.\n'
fi
[ "$SKIP" -gt 0 ] && { printf 'skipped:\n'; printf '  %s\n' "${SKIPPED[@]}"; }
exit $(( FAIL > 0 ))
