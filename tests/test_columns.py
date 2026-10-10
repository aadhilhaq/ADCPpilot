"""Column widths: every table of runs, folders, results and residues can have
its columns resized, and the widths are remembered.

By Aadhil Haq

A grip on the right edge of each heading is dragged with a mouse, a pen or a
finger, moved with the arrow keys, or double clicked to fit its column. The
widths are kept in the browser per table and per column, a Reset column
widths button beside the table forgets them, and a browser that refuses
storage still resizes. What is checked here is the storage, run under Node
where Node can be found (ADCP_NODE, or node on the PATH), and the wiring, by
reading the page; the screenshot script drags the grips in a browser.

    python3 tests/test_columns.py
"""
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import page  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


html = page.render("TOKEN")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
css = html[html.index("<style>"):html.index("</style>")]
ccss = css[css.index("---- column widths"):css.index("---- a phone")]
cols_js = appjs[appjs.index("/* ===================== column widths"):
                appjs.index("/* ===================== the step guide")]

# --- the grip ------------------------------------------------------------------------
check("function colResize(t, key, opt)" in appjs, "there is no column resizer")
for attr in ("g.tabIndex = 0", "'role', 'separator'", "'aria-orientation', 'vertical'",
             "'aria-label', 'Resize the '", "'aria-valuenow'"):
    check(attr in cols_js, "the grip lacks %s" % attr)
# Pointer events, so a mouse, a pen and a finger all work, and a drag that
# leaves the grip is still followed.
check("addEventListener('pointerdown'" in cols_js and "setPointerCapture" in cols_js
      and "'pointercancel'" in cols_js and "mousedown" not in cols_js,
      "the grip does not use pointer events")
check("touch-action:none" in ccss and "cursor:col-resize" in ccss,
      "the grip does not take a finger's drag or show the resize cursor")
check("th:hover>.colgrip::after" in ccss and ".colgrip:focus-visible" in ccss,
      "the grip is not seen on hover, or not seen with the keyboard")
check("addEventListener('dblclick'" in cols_js and "fit(i)" in cols_js,
      "double clicking a grip does not fit its column")
check("COLW_STEP = 10" in cols_js and "'ArrowLeft'" in cols_js and "'ArrowRight'" in cols_js,
      "the arrow keys do not move a column ten pixels")
check("COLW_MIN = 48" in cols_js and "data-min" in cols_js,
      "a column can be dragged to nothing")
# A fixed layout, with the widths on <col>, once anything is resized: the
# results table's first row is its group headings, which span several
# columns, and would otherwise decide the widths.
check("t.classList.add('colfixed')" in cols_js and "createElement('colgroup')" in cols_js
      and "table.colfixed{table-layout:fixed}" in ccss,
      "a resized table does not switch to a fixed layout on its own columns")
check("text-overflow:ellipsis" in ccss and "td.title = text" in cols_js,
      "text cut short has no ellipsis, or no way to read the rest")
check(":where(th.colsz){position:relative}" in ccss,
      "the grip's place would stop a sticky heading sticking")
check("body.colmoving" in ccss and "user-select:none" in ccss,
      "dragging a grip selects the table's text")

# --- every table whose columns can be moved ----------------------------------------
for key, where in (("'runs'", "listRuns"), ("'folders'", "showFolderView"),
                   ("'results'", "openRun"), ("'res-peptide'", "showResidues"),
                   ("'res-receptor'", "showResidues")):
    i = appjs.index("function %s(" % where) if where != "openRun" \
        else appjs.index("async function openRun(")
    body = appjs[i:appjs.index("\n}\n", i)]
    check("colResize(" in body and key in body, "the %s table cannot be resized" % key)
# The pinned columns of the results table follow a resize.
check("after: () => pinResults(el)" in appjs,
      "the results table's pinned columns stay where a resize left them")
check("data-ck=\"' + esc(c.ck) + '\"" in appjs and "ck: g.key + ':' + (c.h || 'act')" in appjs,
      "the results columns are not kept by a key of their own")
check("data-ck=name data-min=120" in appjs,
      "the Folders name column has no room kept for long run names")

# --- the reset ----------------------------------------------------------------------
check("'Reset column widths'" in cols_js and "'ghost sm colreset hide'" in cols_js,
      "there is no Reset column widths button, or it shows before anything is resized")
check("rb.classList.toggle('hide', !Object.keys(colwFor(key)).length)" in cols_js,
      "the reset button does not show only after a resize")
check("colwReset(colwStore(), key);" in cols_js and "unfix();" in cols_js,
      "the reset does not forget the widths and lay the table out again")

# --- storage, every read and write guarded ---------------------------------------------
for line in cols_js.splitlines():
    if "localStorage" in line:
        check("try {" in line, "localStorage is touched outside a try: %s" % line.strip())
check("COLW_HERE" in cols_js, "where storage is refused, a redraw loses the widths")

NODE = os.environ.get("ADCP_NODE") or shutil.which("node")
if NODE and os.path.isfile(NODE if os.sep in NODE else shutil.which(NODE) or ""):
    src = "\n".join(re.search(r"const COLW_KEY = [^;]*;", appjs).group(0).splitlines()) + "\n"
    for f in ("colwAll", "colwGet", "colwClamp", "colwPut", "colwReset"):
        i = appjs.index("function %s(" % f)
        src += appjs[i:appjs.index("\n}\n", i) + 3]
    prog = src + r"""
const mem = () => { const m = new Map(); return {
  getItem: k => m.has(k) ? m.get(k) : null, setItem: (k, v) => m.set(k, String(v)),
  _raw: () => m.get(COLW_KEY)}; };
const out = {};
let s = mem();
out.empty = colwGet(s, 'folders');
out.put1 = colwPut(s, 'folders', 'name', 412.6);
out.put2 = colwPut(s, 'folders', 'size', 3);            // under the minimum
out.put3 = colwPut(s, 'results', 'docking:poses', 99999); // over the maximum
out.put4 = colwPut(s, 'folders', 'actions', 40, 90);     // its own minimum
out.folders = colwGet(s, 'folders');
out.results = colwGet(s, 'results');
out.raw = JSON.parse(s._raw());
out.reset = colwReset(s, 'folders');
out.afterReset = {folders: colwGet(s, 'folders'), results: colwGet(s, 'results')};
// A storage that refuses every read and write: nothing throws.
const refuses = {getItem: () => { throw new Error('denied'); },
                 setItem: () => { throw new Error('denied'); }};
out.refusedGet = colwGet(refuses, 'folders');
out.refusedPut = colwPut(refuses, 'folders', 'name', 300);
out.refusedReset = colwReset(refuses, 'folders');
// No storage at all, and storage holding something else.
out.none = colwGet(null, 'folders');
const junk = mem(); junk.setItem(COLW_KEY, '{not json');
out.junk = colwGet(junk, 'folders');
const arr = mem(); arr.setItem(COLW_KEY, '[1,2]');
out.arr = colwAll(arr);
const odd = mem(); odd.setItem(COLW_KEY, JSON.stringify({folders: {name: -5, size: 'x', when: 120}}));
out.odd = colwGet(odd, 'folders');
out.clamp = [colwClamp(10), colwClamp(10, 90), colwClamp(5000), colwClamp('abc')];
process.stdout.write(JSON.stringify(out));
"""
    r = subprocess.run([NODE, "-e", prog], capture_output=True, text=True, timeout=60)
    got = json.loads(r.stdout) if r.returncode == 0 and r.stdout else None
    check(got is not None, "the width storage could not be run under Node: %s" % r.stderr[-400:])
    if got:
        check(got["empty"] == {}, "a table never resized has widths: %r" % got["empty"])
        check(got["put1"] and got["put2"] and got["put3"] and got["put4"],
              "a width was not kept: %r" % got)
        check(got["folders"] == {"name": 413, "size": 48, "actions": 90},
              "the folders widths read back as %r" % got["folders"])
        check(got["results"] == {"docking:poses": 1200},
              "a width past the maximum was not held to it: %r" % got["results"])
        check(set(got["raw"]) == {"folders", "results"},
              "the widths are not kept per table: %r" % got["raw"])
        check(got["reset"] and got["afterReset"] == {"folders": {},
                                                    "results": {"docking:poses": 1200}},
              "reset did not forget that table alone: %r" % got["afterReset"])
        check(got["refusedGet"] == {} and got["refusedPut"] is False
              and got["refusedReset"] is False,
              "a storage that refuses breaks the resizer: %r" % got)
        check(got["none"] == {} and got["junk"] == {} and got["arr"] == {},
              "no storage, or storage holding something else, is not read as none")
        check(got["odd"] == {"when": 120}, "nonsense widths are kept: %r" % got["odd"])
        check(got["clamp"] == [48, 90, 1200, 48], "the limits are %r" % got["clamp"])
else:
    print("note  no Node here (set ADCP_NODE), so the width storage was checked by "
          "reading it only; the screenshot script drags the grips in a browser")

# --- the docs say so ------------------------------------------------------------------
for name in ("runs.html", "results.html"):
    text = open(os.path.join(ROOT, "adcp_dock", "docs", "pages", name)).read()
    check("Reset column widths" in text and "resize" in text.lower(),
          "the %s page does not say the columns can be resized" % name)

if fails:
    print("test_columns FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  every table's columns can be resized, kept and reset, %d checks" % ran)
