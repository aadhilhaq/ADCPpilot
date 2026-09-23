"""The command line, as a user types it.

Drives each subcommand that submits nothing (new, status, cite, check)
in-process, and checks what it leaves on disk. It also checks that no
module imported at the top of cli.py is imported again inside main(),
which would make that name local to every branch of main().
"""
import contextlib
import io
import json
import os
import re
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import cli                                      # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


def run(*argv):
    """cli.main(argv) -> (exit status, what it printed). A SystemExit is an
    exit status, not a crash; any other exception is the crash being tested
    for, and is reported rather than raised."""
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(out):
        try:
            rc = cli.main(list(argv))
        except SystemExit as e:
            rc = e.code
        except Exception as e:                                  # noqa: BLE001
            rc = "CRASH %s: %s" % (type(e).__name__, e)
    return rc, out.getvalue()


# --- no module is shadowed inside main() -------------------------------------
#
# A module the file imports at the top must never be imported again inside
# main(), because that makes the name local to every branch of it.
src = open(os.path.join(ROOT, "adcp_dock", "cli.py")).read()
body = src[src.index("def main("):]
body = body[:re.search(r"\n(?=def |\Z)", body[1:]).start() + 1] \
    if re.search(r"\n(?=def )", body[1:]) else body
top = set(re.findall(r"^import (\w+)", src, re.M))
local = set(re.findall(r"^\s+import (\w+)\s*$", body, re.M))
check(not (top & local),
      "main() re-imports %s, which makes it local to every branch and "
      "breaks the ones that use it first" % sorted(top & local))

# --- new, then the commands that read what new wrote ------------------------
REC = os.environ.get("ADCP_TEST_RECEPTOR",
                     os.path.join(ROOT, "examples", "receptor.pdb"))
if not os.path.isfile(REC):
    # A receptor is needed only because `new` validates it. The smallest
    # thing it accepts will do; this is not a docking.
    REC = os.path.join(tempfile.mkdtemp(), "tiny.pdb")
    lines = []
    for i, (x, y, z) in enumerate([(0, 0, 0), (1.5, 0, 0), (3, 0, 0),
                                   (4.5, 0, 0)], 1):
        lines.append("ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00  0.00"
                     "           C" % (i, i, x, y, z))
    open(REC, "w").write("\n".join(lines) + "\nEND\n")

rd = os.path.join(tempfile.mkdtemp(), "trial")
rc, out = run("new", rd, "-r", REC, "--chain", "A", "-s", "A:1,A:2",
              "-p", "p53: ETFSDLWKLLP", "-p", "RWGPGAW",
              "--account", "abc123", "--sbatch-extra=--qos=normal")
check(not str(rc).startswith("CRASH"), "adcp-dock new crashed: %s" % rc)
check(os.path.isfile(os.path.join(rd, "job.json")),
      "adcp-dock new wrote no job.json: %s" % out[-300:])
if os.path.isfile(os.path.join(rd, "job.json")):
    j = json.load(open(os.path.join(rd, "job.json")))
    names = [p.get("name") for p in j.get("peptides", [])]
    check(names[:1] == ["p53"], "-p 'p53: ...' did not name the peptide: %r"
          % names)
    check(j.get("account") == "abc123" and j.get("sbatch_extra") ==
          "--qos=normal", "--account / --sbatch-extra were not recorded")
    check(j.get("partition", "") == "",
          "the partition defaulted to %r instead of SLURM's own" %
          j.get("partition"))

for sub in (("status", rd), ("cite", "--rundir", rd), ("cite", "--bibtex"),
            ("cite", "--ris"), ("check",)):
    rc, out = run(*sub)
    check(not str(rc).startswith("CRASH"),
          "adcp-dock %s crashed: %s" % (" ".join(sub[:1]), rc))

rc, out = run("new", os.path.join(tempfile.mkdtemp(), "bad"), "-r", REC,
              "-s", "A:1", "-p", "alpha SWQWWWQPQRDCR")
check(rc not in (0, None) and "colon" in (out + str(rc)),
      "a name glued to a sequence on -p was not refused with the fix: %r"
      % out[-200:])

if fails:
    print("test_cli FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  new, status, cite and check run as typed, %d checks" % ran)
