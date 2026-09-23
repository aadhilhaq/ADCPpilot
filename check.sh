#!/usr/bin/env bash
# Which fixes does this installation actually have?
#
#     ./check.sh                     this install
#     ./check.sh /path/to/another    somewhere else
#
# A commit hash is not an answer to "did the update land". Each line names
# something that went wrong, so the answer is a column of YES and NO.
H="${1:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
say(){ printf '  %-46s %s\n' "$1" "$2"; }
has(){ grep -q "$2" "$H/$3" 2>/dev/null && say "$1" "YES" || say "$1" "NO"; }
echo "checking $H"
echo "  commit $(git -C "$H" rev-parse --short HEAD 2>/dev/null)"
echo

# --- the receptor a peptide is docked into ---------------------------------
# The preparation code is shared with the HOPE peptide pipeline. These lines
# say whether the copy in this checkout is the current one.
has "backbone gaps do not stop protonation"    'def backbone_gaps'   adcp_dock/_vendor/prepare.py
has "minimisation not asked to bond a gap"     'split_for_min'       adcp_dock/_vendor/prepare.py
has "pdb2pqr found without PATH set up"        'def find_tool'       adcp_dock/_vendor/prepare.py
has "tools found from the env, not just PATH"  'sys.executable'      adcp_dock/tools.py
has "every stage reads the prepared receptor"  'def rebind_prepared' adcp_dock/_vendor/receptor.py
has "viewer peptide carries docked charges"    'def with_hydrogens'  adcp_dock/_vendor/poseviewer.py
has "MM-GBSA scores the prepared receptor"     'def _rec_for_mmgbsa' adcp_dock/pipeline.py
has "the page serves the prepared receptor"    'prepared_pdb'        adcp_dock/server.py

# --- can this machine actually run a docking? ------------------------------
echo
# Which adcp_dock is being imported. update.sh installs in editable mode, so
# the checkout should be what runs; a non-editable copy left in site-packages
# by an older install would shadow it, and every pull would then update files
# that nothing reads.
echo "  the package python actually imports:"
for _c in "${ADCP_PYTHON:-}" "$H/../envs/hope/bin/python" \
          "$(command -v python3 2>/dev/null)"; do
  [ -n "${_c:-}" ] && [ -x "${_c:-}" ] || continue
  _where=$("$_c" -c "import adcp_dock,os;print(os.path.dirname(os.path.abspath(adcp_dock.__file__)))" 2>/dev/null)
  if [ -z "$_where" ]; then
    printf "    %-46s not installed\n" "$_c"
  elif [ "$_where" = "$H/adcp_dock" ]; then
    printf "    %s\n" "$_where"
    printf "    ^ this checkout, so a pull reaches the page\n"
  else
    printf "    \033[1;33m%s\033[0m\n" "$_where"
    printf "    ^ NOT this checkout (%s). A pull here changes nothing that\n" "$H"
    printf "      is served. Reinstall:  %s -m pip install --no-deps -e %s\n" "$_c" "$H"
  fi
  break
done
echo

echo "  the tools a job needs:"
_env=""
for _c in "${ADCP_PYTHON:-}" "$H/../envs/hope/bin/python" \
          "$(command -v python3 2>/dev/null)"; do
  [ -n "${_c:-}" ] && [ -x "${_c:-}" ] && { _env="$_c"; break; }
done
if [ -z "$_env" ]; then
  echo "    no python found to ask"
else
  echo "    $_env"
  # activate.sh is what install.sh wrote and what every job sources: it sets
  # ADFRSUITE_HOME and AMBERHOME, so the probe reads it too.
  # shellcheck disable=SC1091
  [ -f "$H/activate.sh" ] && . "$H/activate.sh" >/dev/null 2>&1
  # Asked through the package's own finder, so this reports exactly what a
  # job would find, including the environment's bin, which is where pdb2pqr
  # is and where PATH alone does not look.
  PYTHONPATH="$H:${PYTHONPATH:-}" "$_env" - <<'PY'
try:
    from adcp_dock import tools
except Exception as exc:                                      # noqa: BLE001
    print("    adcp_dock is not importable here (%s)" % exc)
    raise SystemExit(0)
t = tools.discover()
found = dict(t.report())
for n, what in (("adcp", "the docking search"),
                ("agfr", "the grid maps"),
                ("prepare_receptor", "the receptor as PDBQT"),
                ("pdb2pqr", "hydrogens and titration states"),
                ("tleap", "Amber topologies for MM-GBSA"),
                ("sander", "the minimisation before MM-GBSA"),
                ("MMPBSA.py", "the MM-GBSA itself"),
                ("ante-MMPBSA.py", "splitting the complex topology")):
    print("    %-44s %s" % ("%s (%s)" % (n, what), found.get(n) or "NO"))
if not all(found.get(n) for n in ("adcp", "agfr", "prepare_receptor")):
    print("\n    ADFRsuite is missing. On a login node:  ./install.sh")
PY
fi
echo
