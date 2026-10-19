"""reportkit, the lab's report module: one module, the same file in every tool.

By Aadhil Haq

reportkit.js and reportkit.css lay out every lab tool's report: the bar, the
zoom, Save as PDF, the numbering of tables and figures, the print style. Each
tool carries a copy, and the copies are one module, so they are compared byte
for byte with the other lab checkouts named in ADCP_SIBLING_CHECKOUTS (folders
separated as in PATH), as tests/test_surface.py does for surface.py. A folder
may hold the module itself (the staging copy) or be a checkout, whose copies
are <checkout>/<package>/static/reportkit.*.

Also checked: the report page loads the module and starts it; the module holds
what the contract (its README) says it does; and, where Node is found, the
module parses and formats a number as the server's captions do.

    python3 tests/test_reportkit.py
"""
import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import page, runreport  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


STATIC = os.path.join(ROOT, "adcp_dock", "static")
NAMES = ("reportkit.js", "reportkit.css")
mine = {}
for n in NAMES:
    p = os.path.join(STATIC, n)
    check(os.path.isfile(p), "adcp_dock/static/%s is missing" % n)
    mine[n] = open(p, "rb").read() if os.path.isfile(p) else b""

# --- the copies are one module ------------------------------------------------
_raw = os.environ.get("ADCP_SIBLING_CHECKOUTS", "")
if not _raw.strip():
    print("note  ADCP_SIBLING_CHECKOUTS is not set, so reportkit was not compared "
          "with other tools' copies")
else:
    seen = 0
    for where in (p.strip() for p in _raw.split(os.pathsep)):
        if not where or os.path.realpath(where) == os.path.realpath(ROOT):
            continue
        for n in NAMES:
            copies = sorted(set(glob.glob(os.path.join(where, "*", "static", n))
                                + glob.glob(os.path.join(where, n))))
            if not copies:
                # Not a failure: a tool that has no report yet carries none.
                print("note  %s carries no %s" % (where, n))
            for other in copies:
                seen += 1
                with open(other, "rb") as fh:
                    check(fh.read() == mine[n], "adcp_dock/static/%s differs from %s; the copies are one "
                          "module" % (n, other))
    print("note  reportkit compared with %d copies" % seen)

# --- what the module is -------------------------------------------------------
js, css = mine["reportkit.js"].decode("ascii", "replace"), mine["reportkit.css"].decode("ascii", "replace")
check(all(ord(c) < 128 for c in mine["reportkit.js"].decode("utf-8", "replace")),
      "reportkit.js is not plain ASCII: a tool may serve it without a charset")
for what, needle in (("its entry point", "window.ReportKit"), ("Save as PDF", "Save as PDF"),
                     ("the zoom, held to 50 to 200 %", "Math.max(0.5, Math.min(2"),
                     ("the copy of the methods", "Copy the methods and references"),
                     ("the table numbers", "<b>Table "), ("the figure numbers", "<b>Figure "),
                     ("references between blocks", "(table|figure):"),
                     ("structures drawn in 3Dmol", "$3Dmol.createViewer"),
                     ("the views' own words", "contenteditable")):
    check(needle in js, "reportkit.js has no %s (%r)" % (what, needle))
for what, needle in (("a print style", "@media print"), ("A4 pages", "size:A4"),
                     ("the zoom undone in print", "zoom:1!important"),
                     ("tables kept whole", "rk-keep{break-inside:avoid"),
                     ("figures kept whole", "break-inside:avoid"),
                     ("the bar left out of the PDF", "#rk-bar,.rk-noprint{display:none")):
    check(needle in css, "reportkit.css has no %s (%r)" % (what, needle))

# --- the report page loads it --------------------------------------------------
html = page.report_page("TOK", "/some/run")
check('href="/static/reportkit.css' in html and 'src="/static/reportkit.js' in html,
      "the report page does not load reportkit")
check("ReportKit.start(" in html and '"TOK"' in html and '"/some/run"' in html,
      "the report page does not start reportkit with its token and run")
check(html.index("3Dmol-min.js") < html.index("reportkit.js"),
      "3Dmol is not loaded before reportkit, which draws the structures with it")
check("<script>" not in page.report_page("x", "</script><script>alert(1)").split("ReportKit.start(")[1][:200],
      "a run folder's name can close the report page's script")

# --- where Node is found, the module parses and numbers as the server does -----
node = shutil.which("node") or "/sw/eb/sw/nodejs/20.13.1-GCCcore-13.3.0/bin/node"
env = dict(os.environ, LD_LIBRARY_PATH="/sw/eb/sw/GCCcore/13.3.0/lib64:" + os.environ.get("LD_LIBRARY_PATH", ""))
if os.path.isfile(node):
    probe = ("global.window = {}; global.document = {getElementById: () => null};"
             "require(%r); const K = window.ReportKit;"
             "console.log(JSON.stringify([K.num(-69.0123, 1), K.num(2092, 0), K.num(12345.6, 1),"
             " K.num(-0.04, 1), K.num(null, 1), K.num(3, 2)]));" % os.path.join(STATIC, "reportkit.js"))
    r = subprocess.run([node, "-e", probe], capture_output=True, text=True, timeout=60, env=env)
    check(r.returncode == 0, "reportkit.js does not load under Node: %s" % r.stderr[-400:])
    if r.returncode == 0:
        import json
        got = json.loads(r.stdout.strip().splitlines()[-1])
        want = [runreport.num(-69.0123, 1), runreport.num(2092, 0), runreport.num(12345.6, 1),
                runreport.num(-0.04, 1), "", runreport.num(3, 2)]
        check(got == want, "reportkit.js and runreport.num write numbers differently: %s, %s" % (got, want))
else:
    print("note  Node was not found, so reportkit.js was not run")

if fails:
    print("test_reportkit FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  reportkit, the shared report module, %d checks" % ran)
