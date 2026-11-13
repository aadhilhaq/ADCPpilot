"""The command line, as a user types it.

By Aadhil Haq

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

from adcp_pilot import cli                                      # noqa: E402

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
src = open(os.path.join(ROOT, "adcp_pilot", "cli.py")).read()
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
check(not str(rc).startswith("CRASH"), "adcp-pilot new crashed: %s" % rc)
check(os.path.isfile(os.path.join(rd, "job.json")),
      "adcp-pilot new wrote no job.json: %s" % out[-300:])
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
          "adcp-pilot %s crashed: %s" % (" ".join(sub[:1]), rc))

rc, out = run("new", os.path.join(tempfile.mkdtemp(), "bad"), "-r", REC,
              "-s", "A:1", "-p", "alpha SWQWWWQPQRDCR")
check(rc not in (0, None) and "colon" in (out + str(rc)),
      "a name glued to a sequence on -p was not refused with the fix: %r"
      % out[-200:])

# --- the example run: installed, complete, and read by the page -------------
#
# The copy must point at its own folder, with no placeholder left, and every
# panel the page opens on a finished run must find what it reads.
runs = tempfile.mkdtemp()
rc, out = run("example", "--runs", runs)
ex = os.path.join(runs, "example_MDM2")
check(rc in (0, None) and os.path.isfile(os.path.join(ex, "job.json")),
      "adcp-pilot example did not install the example: %s %s" % (rc, out[-300:]))
if os.path.isfile(os.path.join(ex, "job.json")):
    for base, _dirs, names in os.walk(ex):
        for n in names:
            with open(os.path.join(base, n), "rb") as fh:
                check(b"@RUN@" not in fh.read(),
                      "the example's %s still holds the placeholder" % n)
    from adcp_pilot import server
    app = server.App(runs)
    mon = {r["rundir"]: r for r in app.monitor()["runs"]}
    check(ex in mon and mon[ex]["dot"] == "done",
          "the example is not listed as a finished run: %r"
          % {k: v.get("dot") for k, v in mon.items()})
    st = app.run_status(ex)
    rows = (st.get("results") or {}).get("rows") or []
    check([r.get("name") for r in rows if r.get("name")] and len(rows) == 3
          and not st.get("outstanding") and not st.get("prep_warnings"),
          "the example's results are incomplete: %d rows, outstanding %r, "
          "warnings %r" % (len(rows), st.get("outstanding"),
                           st.get("prep_warnings")))
    j = json.load(open(os.path.join(ex, "job.json")))
    for pep in j.get("peptides", []):
        seq = pep["sequence"]
        poses = app.pose_list(ex, seq)["poses"]
        check(len(poses) >= 5 and all(p.get("score") is not None
                                      for p in poses)
              and any(p.get("mmgbsa") is not None for p in poses),
              "the example has no scored poses with MM-GBSA for %s" % seq)
        text, err = app.structure_text(ex, seq, 1)
        check(text and not err and "ATOM" in text,
              "the example's top pose of %s does not open: %s" % (seq, err))
    rc, out = run("cite", ex)
    check(rc in (0, None) and "Amber" in out,
          "the example's citation list does not name AmberTools: %s"
          % out[-300:])
    # A second install leaves the first alone; --force replaces it.
    open(os.path.join(ex, "mine.txt"), "w").write("kept")
    run("example", "--runs", runs)
    check(os.path.isfile(os.path.join(ex, "mine.txt")),
          "a second adcp-pilot example replaced the example already there")
    run("example", "--runs", runs, "--force")
    check(not os.path.isfile(os.path.join(ex, "mine.txt"))
          and os.path.isfile(os.path.join(ex, "job.json")),
          "adcp-pilot example --force did not replace the example")

# Where runs go: $ADCP_RUNS when it, or the folder it would go in, exists (the
# launcher suggests one for a cluster it knows), else $SCRATCH, else home.
_tmp = tempfile.mkdtemp(prefix="adcp_runs_")
_env = {k: os.environ.get(k) for k in ("ADCP_RUNS", "SCRATCH")}
os.environ["SCRATCH"] = _tmp
os.environ["ADCP_RUNS"] = os.path.join(_tmp, "mine", "adcp_runs_v2")
check(cli.default_runs() == os.path.join(_tmp, "adcp_runs"),
      "a suggested runs folder whose parent does not exist was taken: %s" % cli.default_runs())
os.makedirs(os.path.join(_tmp, "mine"))
check(cli.default_runs() == os.path.join(_tmp, "mine", "adcp_runs_v2"),
      "a suggested runs folder whose parent exists was not taken: %s" % cli.default_runs())
for _k, _v in _env.items():
    if _v is None:
        os.environ.pop(_k, None)
    else:
        os.environ[_k] = _v

if fails:
    print("test_cli FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  new, status, cite and check run as typed, %d checks" % ran)
