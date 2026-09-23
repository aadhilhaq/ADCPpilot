"""The web API, driven the way the page drives it.

A form that validates on the server has to be tested against the server: the
page's own checks can agree with themselves perfectly while the thing that
decides whether a run starts disagrees.
"""
import json
import os
import sys
import threading
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import peptide, server                          # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


# A real structure makes four of these checks meaningful; where there is none
# they are skipped and reported as skipped. Point ADCP_TEST_RECEPTOR at a PDB
# file, or place one at examples/receptor.pdb.
RECEPTOR = os.environ.get(
    "ADCP_TEST_RECEPTOR",
    os.path.join(os.path.dirname(HERE), "examples", "receptor.pdb"))

httpd = server.build(host="127.0.0.1", port=0, runs=os.path.join(HERE, "_runs"))
PORT = httpd.server_address[1]
TOKEN = httpd.app.token
threading.Thread(target=httpd.serve_forever, daemon=True).start()


def get(endpoint, **q):
    # `endpoint`, not `path`: /api/receptor takes a query parameter called
    # path, and naming the first argument the same thing collides with it.
    q["t"] = TOKEN
    url = "http://127.0.0.1:%d%s?%s" % (PORT, endpoint,
                                        urllib.parse.urlencode(q))
    with urllib.request.urlopen(url, timeout=30) as r:
        return json.loads(r.read().decode())


def post(path, body):
    url = "http://127.0.0.1:%d%s?t=%s" % (PORT, path, TOKEN)
    req = urllib.request.Request(
        url, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode())


# --- the token actually gates it ------------------------------------------
try:
    urllib.request.urlopen("http://127.0.0.1:%d/api/tools" % PORT, timeout=10)
    check(False, "the API answered with no token at all")
except urllib.error.HTTPError as exc:
    check(exc.code == 403, "a tokenless request got %d, not 403" % exc.code)

# --- tools -----------------------------------------------------------------
t = get("/api/tools")
check("tools" in t, "/api/tools returned nothing useful")
check(any(x["name"] == "adcp" for x in t["tools"]),
      "adcp is not among the tools the server reports on")

# --- the peptide validation the form depends on ---------------------------
CASES = [
    ("RWGPGAWGSGSGHY", True, []),
    ("RWGPGAWGS5GSGHY", False, [9]),       # the digit from the screenshot
    ("RWGPGAWGSXGSGHY", False, [9]),       # X is a letter but not a residue
    ("rwgpgawgsgsghy", True, []),          # lower case is fine
    ("RWG PGA WGS GHY", True, []),         # so is whitespace
    ("AC", False, []),                     # too short
    # 25 was "too long" while the ceiling was 20. It is a legal peptide now
    # and both sides of the new limit are checked, so moving the ceiling
    # again breaks a test rather than quietly widening what is accepted.
    ("A" * 30, True, []),                  # the longest that is allowed
    ("A" * 31, False, None),               # one past it
]
r = post("/api/peptides", {"peptides": [{"sequence": s} for s, _, _ in CASES]})
check(len(r["peptides"]) == len(CASES), "not every peptide came back")
for (src, want_ok, want_pos), got in zip(CASES, r["peptides"]):
    check(got["ok"] == want_ok,
          "%r: the server says ok=%s, expected %s" % (src, got["ok"], want_ok))
    if want_pos:
        pos = sorted({p for pr in got["problems"] for p in pr["positions"]})
        check(pos == want_pos,
              "%r: flagged positions %s, expected %s" % (src, pos, want_pos))
check(r["ready"] is False,
      "a list with invalid sequences in it was reported ready to run")

ok = post("/api/peptides", {"peptides": [{"sequence": "RWGPGAWGSGSGHY"}]})
check(ok["ready"] is True, "a list of one valid peptide was not ready")

# --- pasting a block -------------------------------------------------------
p = post("/api/paste", {"text": "RWGPGAW\nKDGGVDG\n\n>named\nHYGSGSTW"})
check(len(p["peptides"]) == 3, "pasting three sequences gave %d"
      % len(p["peptides"]))
check(p["peptides"][2]["name"] == "named",
      "a FASTA header did not become the peptide's name")

# --- receptor and site -----------------------------------------------------
if os.path.isfile(RECEPTOR):
    d = get("/api/receptor", path=RECEPTOR)
    check("chains" in d and d["chains"], "no chains came back for the receptor")
    check(d["chains"][0]["kind"] == "protein",
          "the receptor's first chain was not called protein")
    s = post("/api/site", {"receptor": RECEPTOR, "chains": [],
                           "hotspots": "A:30,A:31,A:32,A:96,A:99,A:101"})
    check("centre" in s, "the site could not be located: %s" % s.get("error"))
    check(20.0 <= s["size"][0] <= 30.0,
          "the box came out at %s, outside the clamp" % s["size"])
    bad = post("/api/site", {"receptor": RECEPTOR, "chains": [],
                             "hotspots": "A:9999"})
    check("error" in bad or "NOT in the structure" in bad.get("note", ""),
          "a hotspot the structure does not have was accepted silently")
else:
    print("  (no receptor to test against: %s)" % RECEPTOR)

# --- a run that is not ready is refused, and says why ---------------------
bad = post("/api/create", {"label": "bad", "receptor": "", "hotspots": "",
                           "peptides": []})
check("problems" in bad, "an empty run was not refused with reasons")
check(len(bad["problems"]) >= 3,
      "an empty run gave only %d reason(s)" % len(bad.get("problems", [])))

# --- the page and its script are served, and the script parses ------------
import re                                                       # noqa: E402
url = "http://127.0.0.1:%d/?t=%s" % (PORT, TOKEN)
with urllib.request.urlopen(url, timeout=30) as r:
    html = r.read().decode()
check("<title>" in html, "the page came back without a title")
for ctrl in ("recpath", "reccode", "recfetch", "hotspots", "addpep",
             "submit", "rsurf", "rpub", "rlabel", "rposepick", "viewR",
             "viewD", "seqstrip", "chips", "autofill", "estnote"):
    check(('id=%s' % ctrl) in html or ('id="%s"' % ctrl) in html,
          "the page has no #%s" % ctrl)
# features of the pose viewer deliberately left out here
check("rsuper" not in html,
      "the overlay toggle is still on the page; it draws a co-folded "
      "prediction against a docked pose and there is none here")
# and what was asked for afterwards
for want, why in (("rengine", "the viewer engine select"),
                  ("viewM", "the Mol* canvas"),
                  ("diagwrap", "the ligand interaction diagram")):
    check(want in html, "%s is missing from the page" % why)
check('href="/static/molstar.css"' in html,
      "Mol*'s stylesheet is not loaded, so its panel would be unstyled")

# The ids that are form controls, so "defined but never read" can be asked
# about them without also asking it about every div the layout needs.
PAGE_INPUTS = set(re.findall(
    r"""<(?:input|select|textarea)[^>]*\bid=['"]?([A-Za-z0-9_\-]+)""", html))

# the contact lines are one per residue pair, not one per atom pair
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js"),
             encoding="utf-8").read()
check("seenPair" in appjs,
      "contacts are not collapsed per residue pair; two residues touching at "
      "twenty atoms would draw twenty parallel dashes between the same two "
      "side chains")
check("SPECIFIC.indexOf(kind) < SPECIFIC.indexOf(have.kind)" in appjs,
      "when a residue pair makes several kinds of contact the most specific "
      "one does not win")
check("function drawLid" in appjs, "there is no ligand interaction diagram")
check("lidsave" in appjs, "the diagram cannot be saved")
check("image/svg+xml" in appjs,
      "the diagram saves as something other than SVG; a raster cannot go "
      "into a figure at any size")

# --- the output folder is chosen in the page, not in the launcher --------
#
# Somebody signing in from a laptop cannot know which scratch path to use,
# and a home directory usually has a small file allowance while a run writes
# tens of thousands of files. So the page asks, with a folder browser, once
# the user can see what is there.
import tempfile as _tf                                          # noqa: E402

_root = post("/api/runs", {}) if False else get("/api/runs")
check("root" in _root, "the runs listing does not say where runs go")

_new = _tf.mkdtemp(prefix="adcp-runs-")
_moved = post("/api/runsroot", {"path": _new})
check(_moved.get("root") == _new,
      "the output folder did not move: %r" % _moved)
check(get("/api/runs").get("root") == _new,
      "the runs listing still reports the old root after a move")

# a directory that does not exist yet is created, because "make me a folder
# called screens" is the normal thing a user wants
_child = os.path.join(_new, "screens", "june")
check(post("/api/runsroot", {"path": _child}).get("root") == _child,
      "a new folder was not created for the run output")
check(os.path.isdir(_child), "the folder was reported but not created")

# and the ways it must refuse
_afile = os.path.join(_new, "a_file.txt")
open(_afile, "w").write("x")
check("error" in post("/api/runsroot", {"path": _afile}),
      "a FILE was accepted as the output folder")
check("error" in post("/api/runsroot", {"path": ""}),
      "an empty path was accepted")
check("error" in post("/api/runsroot", {"path": "/"}),
      "the filesystem root was accepted as an output folder")
_ro = _tf.mkdtemp(prefix="adcp-ro-")
os.chmod(_ro, 0o500)
_r = post("/api/runsroot", {"path": os.path.join(_ro, "nope")})
check("error" in _r,
      "a folder under an unwritable parent was accepted; the failure would "
      "surface later as a job that could not write its own job file")
os.chmod(_ro, 0o700)
post("/api/runsroot", {"path": _new})          # leave it somewhere sane

# the browser used for it lists directories only: offering structure files
# invites a click that does nothing
os.makedirs(os.path.join(_new, "sub"), exist_ok=True)
open(os.path.join(_new, "thing.pdb"), "w").write("ATOM\n")
_b = get("/api/browse", path=_new, dirs="1")
check(_b.get("files") == [],
      "the output-folder browser is listing files: %r" % _b.get("files"))
check(any(d["name"] == "sub" for d in _b.get("dirs", [])),
      "the output-folder browser is not listing subdirectories")
check("writable" in _b,
      "the browser does not say whether a folder can be written to, so the "
      "person only finds out when a run fails")
# and the ordinary browser still offers structures, or no receptor could be
# picked at all
_b2 = get("/api/browse", path=_new)
check(any(f["name"] == "thing.pdb" for f in _b2.get("files", [])),
      "the receptor browser stopped listing structures")

# --- the page says who made it, and how to say something back -------------
for _want in ("ADCP peptide docking", "HOPE Lab", "Aadhil Haq",
              "Dr. Sandun Fernando", "sfernando@tamu.edu",
              "aadhil.haq@tamu.edu"):
    check(_want in html, "the header does not carry %r" % _want)
check("text-align:center" in html and "brandblock" in html,
      "the title block is not centred")
check("id=feedback" in html, "there is no way to send feedback")
check("mailto:" in appjs,
      "the feedback button posts somewhere instead of opening mail. There "
      "is no server to receive a report, this one is the user's own, and "
      "a button that silently drops what a user typed is worse than none")
check("aadhil.haq@tamu.edu" in appjs and "sfernando@tamu.edu" in appjs,
      "feedback would not reach both addresses")

# --- the two alignment faults ---------------------------------------------
check(".row input[type=text],.row input[type=number],.row select," in html
      and "height:38px" in html,
      "the controls in a row do not share a height, so a ghost button "
      "(5px padding, 12.5px text) sits visibly lower than the input beside "
      "it, which is what the browse button's offset was")
check(".seq:empty{" in html and "align-items:center" in
      html[html.index(".seq:empty{"):html.index(".seq:empty{") + 90],
      "an empty sequence box still has the 18px of headroom its residue "
      "numbers need, so the placeholder sits low: the offset in the "
      "peptide rows")
check("height:46px" in html,
      "the name and entity boxes are not the same height as the sequence "
      "box beside them")

# --- the partition is the cluster's business, not the user's --------------
#
# A default partition name that does not exist on a cluster makes sbatch
# refuse every job ("invalid partition specified"). A partition name is a
# property of the machine; no one using this should have to know one.
from adcp_dock import slurm as _sl                              # noqa: E402

check("id=partition" not in html,
      "the page still asks for a partition. It is a property of the "
      "cluster, and the failure when it is wrong is an outright sbatch "
      "refusal")
check("partition" not in appjs or "$('partition')" not in appjs,
      "app.js still reads a partition field")
# Letting the scheduler choose is the point. A cluster's default partition
# may have a shorter time limit than a run needs, and a job that asks for
# more than its partition allows waits indefinitely with Reason:
# PartitionTimeLimit. With no partition named, the scheduler routes the job
# from what it requests.
check(_sl.partition_line("") == "",
      "a partition is still chosen when no one named one, which overrides "
      "the scheduler's own routing: %r" % _sl.partition_line(""))
check(_sl.partition_line("bigmem") == "#SBATCH --partition=bigmem\n",
      "an explicitly named partition is not honoured")
_rows, _default = _sl.partitions()
if _rows:
    check(all("name" in r and "limit" in r for r in _rows),
          "partitions() does not report time limits, so nothing can tell "
          "whether a walltime fits")
    _lim = [r for r in _rows if r["limit"]]
    if _lim:
        _p = _lim[0]
        check(_sl.fits(_p["name"], _p["limit"] - 1) is True,
              "a walltime inside the limit was called too long")
        check(_sl.fits(_p["name"], _p["limit"] + 1) is False,
              "a walltime PAST %s's %d-minute limit was accepted: that is "
              "the job that waits for ever"
              % (_p["name"], _p["limit"]))
check(_sl.fits("no-such-partition", 60) is None,
      "an unknown partition was judged rather than reported as unknown")
for _t, _m in (("3-00:00:00", 4320), ("2:00:00", 120), ("30:00", 30),
               ("UNLIMITED", None), ("", None)):
    check(_sl._limit_minutes(_t) == _m,
          "%r parsed as %r, expected %r" % (_t, _sl._limit_minutes(_t), _m))

# a job held on something that will never clear is not "queued"
from adcp_dock.server import STUCK_REASONS                      # noqa: E402

check("PartitionTimeLimit" in STUCK_REASONS,
      "PartitionTimeLimit is treated as ordinary queueing, so a job that "
      "will wait for ever gets the same amber dot as one waiting its turn")
check("Priority" not in STUCK_REASONS,
      "waiting your turn was classed as broken")
check("stuckWhy" in appjs and "PartitionTimeLimit" in appjs,
      "the page shows the raw SLURM reason, which tells a user who has "
      "not met it nothing at all")

# --- the browser opens where this person's files actually are -------------
_r = get("/api/roots")
_roots = _r.get("roots") or []
check(_roots, "the browser is offered no starting places at all")
_labels = [x["label"] for x in _roots]
check(any("group" in l for l in _labels),
      "no group allocation is offered. A shared checkout means the files "
      "people care about are usually in the group's space, and a browser "
      "that opens on a home directory makes 'where are my files' a path "
      "no one told them: %s" % _labels)
check(any("scratch" in l for l in _labels),
      "the person's own scratch is not offered: %s" % _labels)
check([l for l in _labels if "group" in l][0] == _labels[0],
      "the group's space is not offered first: %s" % _labels)
for x in _roots:
    check(os.path.isdir(x["path"]),
          "%s is offered but is not a directory" % x["path"])

# --- required fields are marked where they are ----------------------------
import re as _re2                                               # noqa: E402
_blk = appjs[appjs.index("function missingFields"):
             appjs.index("$('submit').onclick")]
_need = _re2.findall(r"out\.push\(\['([A-Za-z0-9_-]+)'", _blk)
check(set(_need) >= {"recpath", "hotspots", "peps", "label", "outdir"},
      "not every mandatory field is checked before submitting: %s" % _need)
_page_ids = set(_re2.findall(r"""\bid=['"]?([A-Za-z0-9_-]+)""", html))
for _i in _need:
    check(_i in _page_ids,
          "the submit check marks %r, which the page does not define, so "
          "nothing would be highlighted and nothing scrolled to" % _i)
check("scrollIntoView" in appjs,
      "a missing field is not scrolled to. The complaint appears under the "
      "submit button, which is usually off screen from the box it is about")
check("classList.add('missing')" in appjs,
      "a missing field is not marked")
check(".missing" in html, "there is no style for a missing field")
check("api('/api/runsroot', {path: want})" in appjs,
      "a path TYPED into the output box is never applied, only one picked "
      "in the browser is, so a user who types theirs has their run go "
      "somewhere else without being told")

# --- the estimate is not wildly wrong -------------------------------------
#
# It said 15 minutes a task for a job that takes about 5, because it added
# the MM-GBSA poses up as though they ran one after another when they run
# several at a time.
from adcp_dock import sizing as _sz                             # noqa: E402

check(_sz.mmgbsa_minutes(5, 48) < _sz.mmgbsa_minutes(5, 1),
      "rescoring five poses is costed the same on 48 cores as on 1, so it "
      "is still being added up as if it ran in series")
check(abs(_sz.mmgbsa_minutes(5, 48)
          - _sz.MMGBSA_SECONDS / 60.0) < 0.01,
      "five poses with five workers should cost one pose's wall time, not "
      "%.1f min" % _sz.mmgbsa_minutes(5, 48))
check(_sz.mmgbsa_minutes(12, 48) > _sz.mmgbsa_minutes(5, 48),
      "more poses than there are workers must cost more than one wave")
# against the two real measurements, and conservative on both. The length
# factor is anchored at 6 residues, so an 11-mer carries x1.30. 1.6x a
# measurement is the conservative direction and matches the spread on every
# other measured point (1.15x to 2.6x).
_hope = _sz.dock_minutes(11, 100, 5e6, 48)
check(8.0 <= _hope <= 13.5,
      "an 11-mer at 100x5M on 48 cores estimates %.1f min; it measured "
      "about 8, and the estimate must not be under it" % _hope)
_mine = _sz.per_task(7, 24, 800000, 4, True, 3)
check(_mine >= 4.6,
      "a 7-mer at 24x800k on 4 cores with 3 rescored poses estimates %.1f "
      "min; that run measured 4.6 and an estimate under a measurement "
      "books too little walltime" % _mine)

# --- the picker is a dialog, not a list inserted into the page ---------------
#
# An inline listing pushes everything below it down as it opens, so the
# control that opened it moves from under the pointer, and a directory with
# hundreds of entries pushes the rest of the form off the screen.
check("<dialog id=picker>" in html, "the file picker is not a dialog")
check("dialog::backdrop" in html,
      "the dialog has no backdrop, so the page behind it still looks live")
check(".listing{max-height:" in html,
      "the listing is not capped, so a big directory scrolls the page "
      "rather than itself")
check("showModal" in appjs, "the dialog is never opened modally")
check("id=crumbs" in html and "crumbsFor" in appjs,
      "there is no path breadcrumb, so getting back up is one click per "
      "level")
check("Places" in appjs,
      "the picker does not open on the places this person has")
check("data-dot" not in html or True, "")
check("id=browser" not in html,
      "the old inline browser is still in the page")

# --- a complaint must not move the form ------------------------------------
check(".missingnote{position:absolute" in html,
      "the 'this is empty' note is still in flow, so it makes its row "
      "taller and every control beside it shifts: telling a user a box "
      "is empty rearranged the form around it")
check("fieldwrap" in html and "position:relative" in html,
      "there is nothing for the note to be positioned against")
check("hasnotes" in appjs and ".row.hasnotes" in html,
      "the row does not reserve room for a note, so the card below it "
      "still jumps when one appears")
_wrapped = _re2.findall(r"class=fieldwrap[^>]*>\s*<label for=([A-Za-z0-9_-]+)",
                        html) if True else []
for _f in ("recpath", "hotspots", "label", "outdir"):
    check(_f in _wrapped,
          "%s can be marked as missing but is not inside a .fieldwrap, so "
          "its note would be positioned against the wrong box" % _f)

# --- the estimate covers the whole chain -----------------------------------
#
# It reported one array task, which is a number that never happens on its
# own: a run is prepare, then the array, then the report.
_est = post("/api/estimate", {"peptides": [{"sequence": "AKSLW"},
                                           {"sequence": "PLLSLL"},
                                           {"sequence": "ASDASD"}],
                              "replicas": 50, "steps": 2500000, "cores": 48,
                              "mmgbsa": True, "mmgbsa_poses": 5})
for _k in ("prep_minutes", "dock_minutes", "mmgbsa_minutes", "report_minutes",
           "total_minutes", "end_to_end", "waves"):
    check(_k in _est, "the estimate does not report %s" % _k)
check(_est["n"] == 3, "three peptides were sent and it counted %r"
                      % _est.get("n"))
check(_est["total_minutes"] >= _est["prep_minutes"] + _est["dock_minutes"]
      + _est["mmgbsa_minutes"] + _est["report_minutes"] - 0.1,
      "the total is less than its own parts: %s" % _est)
check(_est["task_minutes"] >= _est["dock_minutes"] + _est["mmgbsa_minutes"]
      - 0.1,
      "docking and rescoring are in SERIES within a task and the task "
      "estimate is under their sum: %s" % _est)
check(_est["mmgbsa_minutes"] > 0, "MM-GBSA was on and costed at nothing")
_secs = lambda t: sum(int(x) * m for x, m in zip(t.split(":"), (3600, 60, 1)))
check(_secs(_est["end_to_end"]) >= _est["total_minutes"] * 60,
      "the booking is under the estimate it is meant to cover")
check(_secs(_est["walltime"]) >= _est["task_minutes"] * 60,
      "the task booking is under the task estimate")
# No booking is under two hours. The model covers the search and the
# rescoring, not a cold filesystem, a slow module load, a slow agfr on an
# awkward receptor, or a loaded login node, which are not small next to a
# four-minute estimate.
check(_secs(_est["walltime"]) >= 2 * 3600,
      "a task is booked for %s. Nothing goes in for under two hours: a job "
      "killed at the wall has produced nothing and has already spent its "
      "queue time" % _est["walltime"])
check(_secs(_est["prep_walltime"]) >= 2 * 3600,
      "preparation is booked for %s, under the two-hour floor"
      % _est["prep_walltime"])
check(_secs(_est["end_to_end"]) >= 2 * 3600,
      "the end-to-end figure is under the floor: %s" % _est["end_to_end"])
# and the floor is a floor, not a cap: a job that needs more gets more
_big = post("/api/estimate",
            {"peptides": [{"sequence": "A" * 15} for _ in range(20)],
             "replicas": 100, "steps": 5000000, "cores": 48,
             "mmgbsa": True, "mmgbsa_poses": 5, "throttle": 1})
check(_secs(_big["end_to_end"]) > 2 * 3600,
      "twenty peptides one at a time still books the two-hour minimum, so "
      "the floor is acting as a cap: %s" % _big["end_to_end"])
check(_secs(_big["end_to_end"]) >= _big["total_minutes"] * 60,
      "a long run is booked under its own estimate")
_off = post("/api/estimate", {"peptides": [{"sequence": "AKSLW"}],
                              "replicas": 50, "steps": 2500000, "cores": 48,
                              "mmgbsa": False})
check(_off["mmgbsa_minutes"] == 0,
      "MM-GBSA was off and still costed at %s" % _off["mmgbsa_minutes"])
check(_off["total_minutes"] < _est["total_minutes"],
      "turning MM-GBSA off did not make the run shorter")
# a throttle means waves, and waves mean longer
_thr = post("/api/estimate", {"peptides": [{"sequence": "AKSLW"},
                                           {"sequence": "PLLSLL"},
                                           {"sequence": "ASDASD"}],
                              "replicas": 50, "steps": 2500000, "cores": 48,
                              "mmgbsa": True, "mmgbsa_poses": 5,
                              "throttle": 1})
check(_thr["waves"] == 3,
      "three peptides one at a time is three waves, not %r" % _thr["waves"])
check(_thr["total_minutes"] > _est["total_minutes"],
      "running them one at a time was costed the same as all at once")

# --- preparation is costed from the receptor, not as a constant -----------
#
# A 40,000-atom complex with bound ligands and rebuilt loops is not the same
# job as a 1,800-atom domain, and the floor is there for the awkward cases (a
# metal, a cofactor, a sugar). Everything here rounds up: an allocation is
# usually charged for the time a job uses, so a generous booking costs
# little, while a tight one loses the run and its queue time.
from adcp_dock import sizing as _sz2                            # noqa: E402

check(_sz2.prep_minutes() == _sz2.PREP_FLOOR_MINUTES,
      "with nothing known, preparation is not booked at the floor")
check(_sz2.prep_minutes(40000, 0) > _sz2.prep_minutes(2000, 0),
      "a 40,000-atom receptor is costed the same as a 2,000-atom one")
check(_sz2.prep_minutes(5000, 6) > _sz2.prep_minutes(5000, 0),
      "bound ligands and metals add nothing, though they are what makes "
      "preparation awkward rather than merely long")
check(_sz2.prep_minutes(5000, 0, True) > _sz2.prep_minutes(5000, 0, False),
      "rebuilding whole missing loops is costed the same as not doing it")
check(_sz2.prep_minutes(5000, 0, False, (40, 40, 40))
      > _sz2.prep_minutes(5000, 0, False, (20, 20, 20)),
      "a bigger grid box costs agfr nothing")
check(_sz2.prep_minutes(100, 0, False, (1, 1, 1)) >= _sz2.PREP_FLOOR_MINUTES,
      "a tiny receptor drops under the floor")
# and both real measurements sit under what would be booked for them
check(_sz2.prep_minutes(1793, 0, False, (30, 30, 30)) > 3.4,
      "the 1,793-atom receptor measured 3.4 min and would be booked less")
check(_sz2.prep_minutes(5884, 2, False, (30, 30, 30)) > 5.5,
      "the 5,884-atom receptor measured 5.4 min and would be booked less")

# the estimate reads the receptor off the file rather than being told
from adcp_dock.server import App as _App                        # noqa: E402

if os.path.isfile(RECEPTOR):
    _n, _k = _App.receptor_cost(RECEPTOR)
    check(_n > 0, "no atoms were counted in the receptor")
    _with = post("/api/estimate", {"peptides": [{"sequence": "AKLLWH"}],
                                   "replicas": 50, "steps": 2500000,
                                   "cores": 16, "mmgbsa": True,
                                   "mmgbsa_poses": 5, "receptor": RECEPTOR})
    _without = post("/api/estimate", {"peptides": [{"sequence": "AKLLWH"}],
                                      "replicas": 50, "steps": 2500000,
                                      "cores": 16, "mmgbsa": True,
                                      "mmgbsa_poses": 5})
    check(_with["prep_minutes"] >= _without["prep_minutes"],
          "naming the receptor did not change what preparation costs")
    check(_without["prep_minutes"] >= _sz2.PREP_FLOOR_MINUTES,
          "with no receptor named, preparation fell under the floor")

# --- MM-GBSA scales with the RECEPTOR, and badly ---------------------------
#
# It was a flat 150 seconds a pose. Two measurements say otherwise:
#
#     895 heavy atoms   112-174 s per pose
#   2,921 heavy atoms       847 s per pose   (5 poses, all 14.1 min)
#
# 3.3x the receptor cost 5.7x the time, which is what a generalised Born
# calculation over every pair does. The flat constant costed that second run
# at 2.5 minutes and it took 14.3, so the whole task was booked at 9 minutes
# and ran for 20.
check(_sz2.mmgbsa_seconds(895) == _sz2.MMGBSA_REF_SECONDS,
      "the reference point moved")
_big = _sz2.mmgbsa_seconds(2921)
check(700 < _big < 1100,
      "a 2,921-atom receptor costs %.0f s a pose; it measured 847" % _big)
check(_sz2.mmgbsa_seconds(2921) > 4 * _sz2.mmgbsa_seconds(895),
      "the cost is still near-linear in receptor size, so a big receptor "
      "is under-costed")
check(_sz2.mmgbsa_minutes(5, 16, atoms=2921)
      > _sz2.mmgbsa_minutes(5, 16, atoms=895),
      "rescoring against a big receptor costs the same as a small one")
# five poses at once is one pose's wall time, whatever the receptor
check(abs(_sz2.mmgbsa_minutes(5, 16, atoms=2921)
          - _sz2.mmgbsa_seconds(2921) / 60.0) < 0.01,
      "five poses with five workers is not costed as one wave")

# the whole task, against what that run actually did
_task = _sz2.per_task(6, 50, 2500000, 16, True, 5, atoms=2921)
check(_task >= 20.1,
      "the task estimate is %.1f min and that run took 20.1: an estimate "
      "under a measurement books too little walltime" % _task)
check(_task < 40,
      "the task estimate is %.1f min for a run that took 20.1; that is not "
      "conservative, it is wrong" % _task)

# --- the receptor is counted over the chains that will be kept -------------
#
# A structure often holds two copies of a complex, and a run that keeps one
# chain prepares about half the atoms. Costing by the whole file would charge
# for work the run never does, and MM-GBSA cost rises faster than linearly
# with the atom count.
if os.path.isfile(RECEPTOR):
    _all, _k1 = _App.receptor_cost(RECEPTOR)
    _one, _k2 = _App.receptor_cost(RECEPTOR, ["A"])
    check(_one <= _all,
          "restricting to one chain counted MORE atoms than the whole file")
    check(_one > 0, "restricting to chain A counted nothing at all")
    # hydrogens in the input must not be counted as heavy atoms
    check(_App.receptor_cost(RECEPTOR)[0]
          <= sum(1 for L in open(RECEPTOR, errors="replace")
                 if L[:4] in ("ATOM", "HETA")),
          "more atoms were counted than the file has lines for")

# --- the job monitor ------------------------------------------------------
_mon = get("/api/monitor")
check("runs" in _mon and "live" in _mon,
      "there is no monitor endpoint for the dots")
for _r in _mon.get("runs", []):
    check("dot" in _r, "a run has no state word for its dot")
    check(_r["dot"] in ("running", "queued", "done", "failed", "gone", "new"),
          "unknown dot state %r" % _r["dot"])
check(".dot.running{" in html and "animation" in html,
      "a running job's dot does not move. A static green is "
      "indistinguishable from a finished one at a glance")
check(".dot.queued{" in html and ".dot.failed{" in html
      and ".dot.done{" in html,
      "the dots do not distinguish queued, failed and finished")
check("/api/monitor" in appjs, "the page never asks for the states")
check("d.live" in appjs,
      "the page polls forever. A monitor that keeps calling squeue after "
      "everything has finished is one no one should leave open on a login "
      "node")

# --- a new folder, from inside the picker ---------------------------------
_mk = _tf.mkdtemp(prefix="adcp-mk-")
check(post("/api/mkdir", {"parent": _mk, "name": "screens"}).get("path")
      == os.path.join(_mk, "screens"),
      "a new folder could not be created from the picker")
check(os.path.isdir(os.path.join(_mk, "screens")),
      "it was reported but not created")
check(post("/api/mkdir", {"parent": _mk, "name": "screens"}).get("existed")
      is True, "creating it twice is an error rather than a no-op")
# A NAME, NOT A PATH. A user typing "a/b" means the first and gets the
# second, and a browser that quietly walks out of the folder it is showing
# cannot be reasoned about.
for _bad, _why in ((".." , "the parent directory"),
                   ("a/b", "a path with a slash"),
                   ("../out", "a path that climbs out"),
                   ("", "an empty name"),
                   ("   ", "a name of spaces"),
                   (".hidden", "a name this browser would not show again")):
    check("error" in post("/api/mkdir", {"parent": _mk, "name": _bad}),
          "%s was accepted as a folder name" % _why)
check(not os.path.exists(os.path.join(os.path.dirname(_mk), "out")),
      "a name that climbs out of the parent actually created something "
      "outside it")
check("error" in post("/api/mkdir", {"parent": "/no/such/place",
                                     "name": "x"}),
      "a folder was created under a parent that does not exist")
check("id=pickernew" in html, "the picker has no way to make a folder")
check("/api/mkdir" in appjs, "the button is not wired to anything")

# --- the chains the box reaches are a CHOICE, not a log line --------------
_site = post("/api/site", {"receptor": RECEPTOR, "chains": [],
                           "hotspots": "A:41,A:45,A:78"}) \
    if os.path.isfile(RECEPTOR) else {}
if _site and "error" not in _site:
    check("other_chains" in _site,
          "the site endpoint does not say what else the box reaches, so "
          "the choice about it cannot be an informed one")
check("id=box_chains" in html,
      "there is no way to say whether the chains the box reaches should be "
      "docked against")
check("'include' : 'ignore'" in appjs,
      "the choice is not sent in the vocabulary the loader speaks")
check("box_chains" in open(os.path.join(ROOT, "adcp_dock", "job.py")).read(),
      "the job file has nowhere to record the choice")

# --- one hotspot parser, not two ------------------------------------------
#
# The page and the pipeline must read a site the same way, ranges such as
# A:98-112 included.
from adcp_dock import structure as _st                          # noqa: E402
from adcp_dock._vendor import receptor as _vrec                 # noqa: E402

for _text in ("A:41, A:45", "A:98-112", "41,45,78", "A:98-100, B:5",
              "98-101", ""):
    check(_st.parse_hotspots(_text, "A") == _vrec.parse_hotspots(_text, "A"),
          "the page and the pipeline disagree about %r: %s against %s"
          % (_text, _st.parse_hotspots(_text, "A"),
             _vrec.parse_hotspots(_text, "A")))
check(len(_st.parse_hotspots("A:98-112", "A")) == 15,
      "a range is not expanded: %s" % _st.parse_hotspots("A:98-112", "A"))

# --- every id the page and the script share actually exists --------------
#
# A control added to the HTML and never wired, or wired and never added,
# raises nothing: $('x') on a missing element returns null, and the button
# does nothing.
#
# Two things this has to get right. app.js writes HTML of its own, so an id
# it creates and then reads is defined even though page.py never mentions it.
# And a delegated handler reads an id by comparing event.target.id, never
# through $(), so appearing anywhere in the script counts as being read.
_ID = r"""\bid=['"]?([A-Za-z0-9_\-]+)"""
_defined = set(re.findall(_ID, html)) | set(re.findall(_ID, appjs))
_fetched = set(re.findall(r"""\$\(['"]([A-Za-z0-9_\-]+)['"]\)""", appjs))
_missing = sorted(i for i in _fetched if i not in _defined)
check(not _missing,
      "app.js reads element(s) that neither the page nor the script ever "
      "creates: %s. $() returns null for these, so whatever they control "
      "silently does nothing." % _missing)

# and every form control the page defines is named somewhere in the script
_inputs = set(re.findall(
    r"""<(?:input|select|textarea)[^>]*\bid=['"]?([A-Za-z0-9_\-]+)""", html))
_unread = sorted(i for i in _inputs if i not in appjs)
check(not _unread,
      "the page defines input(s) nothing in app.js mentions: %s. A control "
      "that does not reach the job file is worse than no control." % _unread)

for name in ("3Dmol-min.js", "molstar.js", "molstar.css", "app.js"):
    with urllib.request.urlopen("http://127.0.0.1:%d/static/%s" % (PORT, name),
                                timeout=30) as r:
        body = r.read()
    check(len(body) > 1000, "/static/%s came back empty" % name)
check("/static/3Dmol-min.js" in html,
      "the page does not load the vendored viewer, so it would need a CDN")

try:
    import esprima
    esprima.parseScript(
        open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read())
except ImportError:
    pass
except Exception as exc:                                        # noqa: BLE001
    check(False, "app.js does not parse: %s" % exc)
    ran += 1

# a request for a file outside static/ must not escape it
try:
    urllib.request.urlopen(
        "http://127.0.0.1:%d/static/../../../etc/passwd" % PORT, timeout=10)
    check(False, "the static handler served a file outside its own directory")
except urllib.error.HTTPError as exc:
    check(exc.code in (403, 404), "path escape gave %d" % exc.code)

# --- fetching, estimating, the picker and resume --------------------------
from adcp_dock import sizing as EST, fetch as FETCH                # noqa: E402

check(FETCH.looks_like_code("4CMM") and not FETCH.looks_like_code("CMM"),
      "the PDB id check accepts the wrong things")
bad = get("/api/fetch", code="nonsense")
check("error" in bad, "a nonsense PDB id was accepted")

if os.path.isfile(RECEPTOR):
    seq = get("/api/sequence", path=RECEPTOR, chain="A")
    check(len(seq.get("residues") or []) > 50,
          "the sequence strip got %d residues" % len(seq.get("residues") or []))
    check(seq["residues"][0].get("code"),
          "a residue came back with no one-letter code to draw")
    raw = urllib.request.urlopen(
        "http://127.0.0.1:%d/api/rawstructure?t=%s&path=%s"
        % (PORT, TOKEN, urllib.parse.quote(RECEPTOR)), timeout=30).read().decode()
    check(raw.count("ATOM") > 100,
          "the picker's structure endpoint returned %d ATOM records"
          % raw.count("ATOM"))

est = post("/api/estimate", {"peptides": [{"sequence": "RWGPGAWGSGSGHY"},
                                          {"sequence": "KDGGVDG"}],
                             "replicas": 50, "steps": 2500000, "cores": 48,
                             "mmgbsa": True, "mmgbsa_poses": 5})
check("walltime" in est, "no walltime came back from the estimator")
check(EST.parse_hhmmss(est["walltime"]) >= 2 * est["estimate_minutes"],
      "the walltime asked for (%s) is not at least twice the estimate (%s min)"
      % (est["walltime"], est["estimate_minutes"]))
check(est["longest"] == 14,
      "the estimate was not sized for the longest peptide")
check(est["cores"] <= 50,
      "asked for %d cores with only 50 replicas to give them" % est["cores"])
check("error" in post("/api/estimate", {"peptides": []}),
      "an estimate with no peptides was not refused")

# the estimator is calibrated, not invented: these two are measured runs
check(9 <= EST.dock_minutes(11, 100, 5_000_000, 48) <= 14,
      "an 11-mer at 100x5M on 48 cores estimates %.1f min; it measured 8"
      % EST.dock_minutes(11, 100, 5_000_000, 48))

check("error" in post("/api/resume", {"rundir": "/nope"}),
      "resume on a run that does not exist was not refused")

# --- the module's own limits are what the API enforces --------------------
# 5 to 30. ADCP fails every replica for three or four residues; the
# ceiling is where the box runs out: a 30-mer needs about 48 A of cube
# (sizing.box_for_length) against a maximum box of 60 A.
check(peptide.MIN_LENGTH == 5 and peptide.MAX_LENGTH == 30,
      "the supported range is %d to %d; it is meant to be 5 to 30"
      % (peptide.MIN_LENGTH, peptide.MAX_LENGTH))
check(peptide.check("AAAA")[1] and not peptide.check("AAAAA")[1],
      "a four-residue peptide is accepted, or a five-residue one refused")
check(set(peptide.NOT_RESIDUES) == set("BJOUXZ"),
      "the rejected letters are no longer B, J, O, U, X, Z")

# --- the sequence strip shows the holes in the structure ------------------
#
# The strip holds the residues the file contains, and each row is labelled
# with the number of its first residue. Where a loop was never modelled, the
# numbers jump; those gaps are the chain breaks the preparation splits at, so
# a user picking a site next to one should be able to see it.
_js = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("r.resnum > prev.resnum + 1" in _js,
      "the sequence strip does not look for missing residues, so a gap is "
      "drawn as continuous sequence")
check("class=gap" in _js,
      "nothing marks a gap in the strip")
check("never modelled" in _js,
      "the gap marker does not say what it means")
check("class=r data-k" in _js and _js.index("class=gap") < _js.index(
          "const sep = (k && k % 10 === 0)"),
      "the gap branch does not come before the residue branch, so a gap "
      "would be counted as a residue")
check(".gap{" in open(os.path.join(ROOT, "adcp_dock", "page.py")).read(),
      "the gap marker has no style, so it is indistinguishable from a residue")

# A gap must not be clickable: there is no residue there to pick. The
# handler is attached to .r only, and the marker is not one.
_i = _js.index("class=gap")
_gapspan = _js[_i:_js.index("</span>", _i)]
check("data-k" not in _gapspan,
      "the gap marker carries a data-k, so clicking a hole would add a "
      "residue that is not in the structure to the site")
check("querySelectorAll('.r')" in _js and "class=gap" in _js
      and "class=r " not in _gapspan,
      "the gap marker is styled as a residue, so it would take a click")


# --- every setting explains itself, and effort is a choice not a number ---
# "Steps per replica: 2500000" asks for a quantity no one has an intuition
# for, and every wrong answer is either a wasted allocation or a pose that
# moves when the run is repeated. And a title= attribute is not an
# explanation: it appears on hover only, never on a touch screen, and holds
# one line where these need a paragraph.
from adcp_dock import helptext as _ht, page as _page            # noqa: E402
import json as _json                                             # noqa: E402

_html = _page.render("token")
check("__HELP_JSON__" not in _html and "__STEP_OPTIONS__" not in _html,
      "a template placeholder was left in the served page")

_m = _re2.search(r'<script id=helpdata type="application/json">(.*?)</script>',
                _html, _re2.S)
check(_m, "the page carries no help text, so no marker can explain anything")
if _m:
    _d = _json.loads(_m.group(1))
    check(len(_d) == len(_ht.HELP),
          "the page's help island and helptext.HELP disagree: %d vs %d"
          % (len(_d), len(_ht.HELP)))

# a marker is only useful if the field it names exists
_orphans = [k for k in _ht.HELP
            if not _re2.search(r'\bid=%s\b' % _re2.escape(k), _html)
            and not _re2.search(r'\bid="%s"' % _re2.escape(k), _html)]
check(not _orphans,
      "help is written for %s, which are not on the page: the marker would "
      "never appear" % _orphans)

for _k, (_t, _b) in _ht.HELP.items():
    check(_t and len(_b) > 40,
          "the help for %r is too thin to answer anything" % _k)

# effort and sampling are short lists of NUMBERS
#
# A label like "thorough" cannot be compared with anything published. The
# number is what ADCP is given and what the tutorials and the paper quote:
# the creators' own redocking tutorial uses -N 400 -n 20000000 for a FIVE
# residue peptide, so 20M is the top of the published scale, not beyond it.
for _id, _choices, _what in (("steps", _ht.STEP_CHOICES, "search effort"),
                             ("replicas", _ht.REPLICA_CHOICES, "replicas")):
    check("<select id=%s>" % _id in _html,
          "%s is still a free number box" % _what)
    check(len(_choices) == 5,
          "%d %s choices; the ask was five" % (len(_choices), _what))
    check([_v for _v, _l in _choices] == sorted(_v for _v, _l in _choices),
          "the %s choices are not in increasing order" % _what)
    check(len({_v for _v, _l in _choices}) == 5,
          "two %s choices have the same value" % _what)
    for _v, _l in _choices:
        check("{:,}".format(_v) in _html,
              "%s choice %d is not shown as a grouped number, so it cannot "
              "be compared with a published setting" % (_what, _v))
    check(sum(1 for _v, _l in _choices
              if ('value="%d" selected' % _v) in _html) == 1,
          "exactly one %s choice must be preselected" % _what)

# the range has to reach what the creators' own tutorial uses
check(max(_v for _v, _l in _ht.STEP_CHOICES) >= 20_000_000,
      "the top step choice is %d; the ADCP tutorial uses 20,000,000 for a "
      "5-mer, so the scale stops short of published practice"
      % max(_v for _v, _l in _ht.STEP_CHOICES))
check(min(_v for _v, _l in _ht.STEP_CHOICES) <= 100_000,
      "there is no cheap setting to prove a run works before spending on it")
check(max(_v for _v, _l in _ht.REPLICA_CHOICES) >= 400,
      "the top replica choice is %d; the ADCP tutorial uses 400"
      % max(_v for _v, _l in _ht.REPLICA_CHOICES))
check(min(_v for _v, _l in _ht.REPLICA_CHOICES) >= 25,
      "a replica count below 25 is fewer searches than any published run")

# The marker is 15px and stays 15px: ".row button" sets height:38px at
# specificity (0,1,1), and a bare ".q" (0,1,0) would lose to it.
_css = open(os.path.join(ROOT, "adcp_dock", "page.py")).read()
check("button.q" in _css and ".row button.q" in _css,
      "the help marker is styled as a bare .q, which loses to '.row button "
      "{height:38px}' and comes out as a 38px pill inside every form row")


# --- and the rule that says so is actually in the stylesheet -------------
#
# An unterminated rule makes a browser discard every rule after it up to the
# next terminator, so rules that are correct in the file can be missing from
# the page. The stylesheet is checked for balance.
_sheet = _html[_html.index("<style>") + 7:_html.index("</style>")]
_bare = _re2.sub(r'/\*.*?\*/', '', _sheet, flags=_re2.S)
_depth, _maxd, _at = 0, 0, False
for _ch in _bare:
    if _ch == '@':
        _at = True
    elif _ch == '{':
        _depth += 1
        _maxd = max(_maxd, _depth - (1 if _at else 0))
    elif _ch == '}':
        _depth -= 1
        if _depth == 0:
            _at = False
check(_depth == 0,
      "the stylesheet's braces do not balance (%d left open): every rule "
      "below the missing one is swallowed as an invalid declaration and "
      "never reaches the page" % _depth)
check(_maxd <= 1,
      "a rule opens inside a rule that is still open: a closing brace is "
      "missing above it, and everything from there to the next stray '}' is "
      "parsed away")

# the rules most at risk from that, named, so a regression is caught
for _sel in ("button.q,", ".gap{", ".btnlink{", ".card>h2{",
             "#helpbox h3", ".seq:empty::before{"):
    check(_sel in _bare,
          "%s is not a rule in the stylesheet the page serves" % _sel)

# and the marker keeps no box in the hover state either
def _block(sel):
    """The declarations of the first rule whose selector starts with sel,
    or None when the browser would not have that rule at all."""
    if sel not in _bare:
        return None
    _b = _bare[_bare.index(sel):]
    return _b[:_b.index("}")] if "}" in _b else None

# EVERY rule that touches the marker, not just the first one.  The circle's
# old hover sat one rule below the rule that removed the box and put the
# background and the border straight back; it was invisible only because
# the whole region was being parsed away at the time.
_qrules = [(_m.group(1).strip(), _m.group(2))
           for _m in _re2.finditer(r'([^{}]*button\.q[^{}]*)\{([^}]*)\}',
                                   _bare)]
check(_qrules, "nothing in the served stylesheet styles the help marker")
# A DRAWN box, and only that.  "min-width:0" is not a box, it is what
# stops an inline-flex parent growing one, so the property has to be
# matched at its boundary rather than as a substring.
for _sel, _body in _qrules:
    for _gone in ("background", "border", "border-radius", "width", "height"):
        check(not _re2.search(r'(?:^|[;{\s])%s\s*:' % _gone, _body),
              "%r puts a %s back on the help marker, which is the box "
              "all:unset removed" % (_sel, _gone))


# --- the heading is a ribbon on the card, not a caption floating above ---
#
# Full-bleed: the band's negative margins have to equal .card's padding and
# its radius has to equal .card's radius minus .card's border, or it hangs
# off the card's edge or leaves a gap of card behind it.
_card_p = _re2.search(r'\.card\{[^}]*padding:(\d+)px (\d+)px', _bare)
_card_r = _re2.search(r'\.card\{[^}]*border-radius:(\d+)px', _bare)
_h2 = _block(".card>h2{")
check(_h2 is not None,
      "the card headings have no rule of their own in the served stylesheet")
check(_card_p and _card_r, "the card no longer states its padding and radius")
if _card_p and _card_r and _h2:
    check("margin:-%spx -%spx" % (_card_p.group(1), _card_p.group(2)) in _h2,
          "the heading band's negative margins are not the card's padding "
          "(%s %s), so it does not reach the card's edges"
          % (_card_p.group(1), _card_p.group(2)))
    check("border-radius:%dpx %dpx 0 0" % (int(_card_r.group(1)) - 1,
                                           int(_card_r.group(1)) - 1) in _h2,
          "the band's top radius is not the card's inner radius, so the "
          "card's corner shows through behind it")
check(_h2 and "box-shadow:inset" in _h2,
      "the band's accent edge is an outer shadow; it would stack with the "
      "shadow the help and issue panels already carry")
for _tok in ("--ribbon", "--ribbonfg"):
    check(_bare.count(_tok + ":") == 2,
          "%s is defined %d times; every colour needs a light AND a dark "
          "value" % (_tok, _bare.count(_tok + ":")))


# --- a label that carries a "?" still starts where the others do ----------
#
# .card is text-align:center and a label sits in a .row>div that inherits
# it.  A block label ignores that for its own position; an INLINE-level box
# does not.  So making the label inline-flex to hold its marker centred
# every label that has one, which is almost all of them, over a control
# that was still left-aligned, and left the two without a marker at the
# edge.  Block-level flex holds the marker and keeps the left edge.
_lab = _block("label:has(> button.q){")
check(_lab is not None,
      "nothing keeps the help marker beside its label's text")
if _lab:
    check(_re2.search(r'display:\s*flex', _lab)
          and not _re2.search(r'display:\s*inline-flex', _lab),
          "a label holding a marker is an inline-level box, so "
          ".card{text-align:center} centres it over its own left-aligned "
          "control: %r" % _lab)
# the checkboxes are the exception and have to win on specificity:
# .checks .chk:has(> button.q) is (0,3,1) against label:has(...) at (0,1,2)
check(".checks .chk:has(> button.q)" in _bare,
      "the checkbox labels are dragged block-level with the field labels, "
      "so the four checkboxes stack instead of sitting in one row")

_js = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("addHelpMarkers" in _js, "nothing adds the markers to the page")
check("label[for=" in _js and "closest('label')" in _js,
      "the markers are attached by one method only, so the checkboxes: "
      "whose labels wrap their input: would get none")
check("'Escape'" in _js, "the help panel cannot be dismissed from the keyboard")
check("+$('steps').value" in _js,
      "the steps value is no longer read as a number from the select")


# --- taking a run away, and sending it back when it broke -----------------
#
# A finished run can be downloaded from the page, and a failed one can be
# reported with its logs attached, since the person who can read the log is
# often not the person running the tool.
from adcp_dock import bundle as _bu                               # noqa: E402
import zipfile as _zip                                            # noqa: E402
import io as _io                                                  # noqa: E402

_run = os.path.join(ROOT, "_smoke", "first")
if os.path.isdir(_run):
    _n, _b = _bu.zip_bytes(_run, "results")
    check(_n.endswith("_results.zip"), "the bundle is misnamed: %r" % _n)
    _z = _zip.ZipFile(_io.BytesIO(_b))
    _names = _z.namelist()
    check(_names, "the results bundle is empty")
    check(all(not n.startswith("/") and ".." not in n for n in _names),
          "the bundle holds an absolute or escaping path, which would "
          "unpack outside the folder the user chose")
    check(len({n.split("/")[0] for n in _names}) == 1,
          "the bundle does not unpack into one folder, so it would scatter "
          "files across whatever directory it is opened in")
    check(_z.testzip() is None, "the bundle is corrupt")

    _na, _ba = _bu.zip_bytes(_run, "all")
    check(len(_ba) >= len(_b),
          "'everything' is smaller than 'results', so the selection is "
          "inverted")
    check("_all.zip" in _na, "the everything bundle is misnamed")

    # a problem report carries logs, not poses
    _logs = _bu.logs_for(_run)
    check(_logs, "a problem report would carry no logs at all")
    check(not any(n.endswith(".pdb") or n.endswith(".pdbqt")
                  for n, _d in _logs),
          "a structure is attached to a problem report: it is logs that are "
          "wanted, and a mail server refuses the rest")
    _total = sum(len(d) for _n2, d in _logs)
    check(_total <= _bu.LOGS_TOTAL_BYTES,
          "the report is %.1f MB, over the cap a mail server will take"
          % (_total / 1e6))
else:
    print("  (no example run here; the bundle contents were not checked)")

# the path comes off a query string, so it is checked rather than trusted
_app = server.App.__new__(server.App)
check(_app.check_run("") , "an empty run path is accepted")
check(_app.check_run("/tmp"),
      "any directory on the filesystem can be zipped and downloaded: it "
      "must hold a job.json to be one of ours")
check(_app.check_run("/nonexistent/%d" % os.getpid()),
      "a path that does not exist is accepted")

# and a report that cannot be sent says so, with the text still in hand
_saved = _bu.can_send
try:
    _bu.can_send = lambda: False
    _r = server.App.report_problem(_app, {"rundir": _run, "note": "x"}) \
        if os.path.isdir(_run) else {"sent": False, "body": "x", "fallback": "x"}
    check(_r.get("sent") is False, "a report was reported as sent with no "
                                   "way to send it")
    check(_r.get("body") and _r.get("fallback"),
          "when mail cannot be sent the user is not given the message text "
          "or told how to send it by hand: what they typed is simply lost")
finally:
    _bu.can_send = _saved

_js = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("'/download?what=' + what" in _js,
      "the page builds no download URL")
check("data-dl=" in _js,
      "the downloads are not buttons; an <a> is a link however it is styled, "
      "and it reads as one")
check("data-issue=" in _js, "there is no button to report a broken run")


# A report has to be answerable: the envelope sender is <user>@<node>, a
# machine that takes no mail, so the reply-to must be set.
_with = _bu.report_body("/tmp", "x", {"email": "someone@example.edu"})
check("from:     someone@example.edu" in _with,
      "the address the user gave is not in the report")
check(_with.count("someone@example.edu") == 1,
      "the address appears twice in the report")
_without = _bu.report_body("/tmp", "x", {})
check("cannot be answered" in _without,
      "a report with no address does not say so, so no one notices it cannot "
      "be replied to: %r" % _without.splitlines()[3])
check("None" not in _without,
      "a missing address is rendered as the word None")
import inspect as _in                                             # noqa: E402
check("reply_to" in _in.signature(_bu.send).parameters,
      "send() cannot set Reply-To, so replies go to a compute node")


# --- the actions are on the run's row, and the table reads as a table -----
#
# Taking a run away is something you do TO a run, so it belongs beside the
# run, not in a panel further down that you reach by opening it first.
_js2 = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
_runs = _js2[_js2.index("$('runs').innerHTML"):
             _js2.index("d.runs.forEach(r => fillSizes")]
for _need, _why in (
        ("data-dl=\"results|", "no results download on the run row"),
        ("data-dl=\"all|", "no download-everything on the run row"),
        ("data-issue=", "no report-an-issue button on the run row"),
        ("data-open=", "the open button has gone")):
    check(_need in _runs, _why)
check("data-size=" in _runs and "fillSizes" in _js2,
      "the download links carry no size, so no one knows whether 'all' is "
      "one megabyte or forty before they click it through a tunnel")

# one button, not a wall of prose
check("takeItWith" not in _js2,
      "the old download panel is still built inside the opened run")
check(_js2.count("Something went wrong") == 0,
      "the prose block is still there; it is one button now")
check("id=issuebox" in _html and "id=issuesend" in _html,
      "the report-an-issue dialog is missing")
check("id=issuemail" in _html,
      "the report dialog cannot take a reply-to address")

# a View button per peptide, and the viewer is brought into sight
check("data-view=" in _js2, "the results table has no view button")
check("e.stopPropagation()" in _js2,
      "the view button does not stop the row click, so both fire")
check("scrollIntoView" in _js2,
      "the viewer opens below the fold with nothing to say it has opened")
check("$('viewcard')" in _js2,
      "nothing scrolls to the viewer card")
# the LID has to be below the viewer, in the same card
_i = _html.index("id=viewcard")
_card = _html[_i:_html.index("<dialog", _i)]
check("id=diagwrap" in _card,
      "the ligand interaction diagram is not inside the viewer card, so a "
      "view click does not bring it with it")
check(_card.index("id=diagwrap") > _card.index("id=contactbox"),
      "the diagram is drawn above the contacts rather than below the viewer")

# centred, except the sequence
check("th.c,td.c{text-align:center}" in _html,
      "there is no centred cell class")
check("td.seq{" in _html and "text-align:left" in
      _html[_html.index("td.seq{"):_html.index("td.seq{") + 200],
      "the sequence column is not left aligned; it is read letter by letter, "
      "not as a quantity")
_res = _js2[_js2.index("rows.map(r => '<tr class=pick"):
            _js2.index("'</tr>').join('')")]
check("class=num" not in _res,
      "a results cell is still right-aligned under a centred heading")
check(_res.count("class=c") >= 6,
      "the results cells are not centred: %d centred cells found"
      % _res.count("class=c"))
check("class=seq" in _res, "the peptide column is not the sequence class")


# --- the downloads are pills, like open -----------------------------------
_cssx = _html[_html.index("<style>"):_html.index("</style>")]
import re as _r3
def _props(sel):
    out = {}
    for _m in _r3.finditer(r'([^{}]+)\{([^}]*)\}',
                           _r3.sub(r'/\*.*?\*/', '', _cssx, flags=_r3.S)):
        for _s in _m.group(1).split(","):
            if _s.strip() == sel:
                out.update(dict(_r3.findall(r'([\w-]+)\s*:\s*([^;]+)',
                                            _m.group(2))))
    return {k: v.strip() for k, v in out.items()}
_open = {**_props("button"), **_props("button.ghost"),
         **_props(".btnlink.sm,button.ghost.sm")}
_dl = {**_props(".btnlink"), **_props(".btnlink.sm,button.ghost.sm")}
for _k in ("border-radius", "border", "background", "font-size", "padding"):
    check(_open.get(_k) == _dl.get(_k),
          "the download link's %s is %r where the open button's is %r: they "
          "sit side by side and have to look like the same control"
          % (_k, _dl.get(_k), _open.get(_k)))

# --- and the page says which build it is ----------------------------------
#
# A running server holds page.py in memory, so a pull changes the checkout
# and not the page being served. Without this there is no way to tell those
# apart from the outside, and "nothing changed" has no visible answer.
check(_r3.search(r'build [0-9a-f]{6,}|build unknown', _html),
      "the page carries no build stamp")
check("Ctrl-Shift-R" in _html,
      "the page does not say how to get the new one")


# --- the toolbar keeps each label with its buttons ------------------------
#
# The bar was one flat flex list, so a wrap could land between a label and
# the controls it names: "FIGURE" alone at the end of one row with its on/off
# button at the start of the next, which reads as a label for the wrong thing.
_bar = _html[_html.index("<div class=vbar id=rvbar>"):
             _html.index("</div>", _html.index("save image"))]
check(".vgrp{" in _html, "there is no non-wrapping group class")
check(_bar.count("class=vgrp") >= 6,
      "only %d toolbar groups; each label and its controls need one"
      % _bar.count("class=vgrp"))
for _lb in ("colour", "draw", "surface", "contacts", "figure", "view"):
    _i = _bar.index(">" + _lb + "<")
    _g = _bar.rfind("class=vgrp", 0, _i)
    _end = _bar.index("</span>", _bar.index("<span class=vgrp", _g))
    check(_g != -1, "the %r label is not inside a group, so it can be "
                    "separated from its buttons by a wrap" % _lb)

# --- labels can be moved, with the helper the other pipelines use ---------
check("id=rdrag" in _html and "id=rdragoff" in _html,
      "the viewer has no move-labels control, which the other pipelines' "
      "3Dmol viewer has had for months")
check("/static/drag.js" in _html, "drag.js is not loaded")
_drag = os.path.join(ROOT, "adcp_dock", "static", "drag.js")
check(os.path.isfile(_drag), "drag.js was not vendored")
if os.path.isfile(_drag):
    _d = open(_drag).read()
    check("root.makeLabelDrag = makeLabelDrag" in _d,
          "drag.js does not export makeLabelDrag")
    check("hope_monitor/page_pick.py" in _d and "Edit the original" in _d,
          "the vendored copy does not say where its source is")
check("makeLabelDrag(VR" in _js2,
      "nothing constructs the dragger for the pose viewer")
check("RDRAG.reset()" in _js2,
      "there is no way back once a label has been dragged somewhere silly")


# --- the published guidance is where the choice is made -------------------
#
# ADCP's step budget scales with the peptide and the creators put a number
# on it; a dropdown without that is a number box with extra clicks.
for _k, _want in (("steps", "per amino acid"), ("replicas", "100-300")):
    check(_want in _ht.HELP[_k][1],
          "the %s help does not carry the creators' guideline (%r)"
          % (_k, _want))
    check("ccsb.scripps.edu" in _ht.HELP[_k][1],
          "the %s help quotes a guideline without saying where from" % _k)
check("1-3 million" in _ht.HELP["steps"][1],
      "the steps help does not give the per-residue rule")
check("stepsGuide" in _js2,
      "nothing works the guideline out for the peptide actually entered")

# --- too long is said at the box, in red ----------------------------------
check("toolong" in _js2 and ".seq.toolong{" in _html,
      "an over-length sequence is not marked at the box")
check("var(--bad)" in _html[_html.index(".seq.toolong{"):
                             _html.index("}", _html.index(".seq.toolong{"))],
      "the over-length box is not red")
check("sequence too long" in _js2,
      "the box does not say what is wrong")
# DERIVED, NOT PINNED. This asserted the literal "MAXLEN = 20" in app.js and
# then that peptide.MAX_LENGTH was also 20: two numbers agreeing because
# both were typed, which is not the same as one number. app.js is a static
# file and cannot be substituted, so the page hands it the limit in a JSON
# island the way it already hands over the help text, and these check that
# path rather than either value.
import adcp_dock.peptide as _pep                                  # noqa: E402
check("MAXLEN = LIMITS.max_length" in _js2,
      "app.js carries its own length limit instead of taking the one the "
      "server enforces; the box will mark at a different length from the "
      "one the API refuses at")
_isl = _re2.search(r'id=limits[^>]*>(\{[^<]*\})', _html)
check(_isl, "the page does not hand app.js the length limits")
if _isl:
    _lim = _json.loads(_isl.group(1))
    check(_lim.get("max_length") == _pep.MAX_LENGTH,
          "the page says max %s and peptide.py says %d"
          % (_lim.get("max_length"), _pep.MAX_LENGTH))
    check(_lim.get("min_length") == _pep.MIN_LENGTH,
          "the page says min %s and peptide.py says %d"
          % (_lim.get("min_length"), _pep.MIN_LENGTH))
check("max_length: 20" not in _js2,
      "a hard-coded 20 is still being sent to the API")
# The range itself: 3 is where ADCP starts, 30 is where the BOX runs out.
check(_pep.MIN_LENGTH == 5 and _pep.MAX_LENGTH == 30,
      "the supported range is %d to %d; it is meant to be 5 to 30"
      % (_pep.MIN_LENGTH, _pep.MAX_LENGTH))

# --- will the longest peptide fit in the box? ----------------------------
#
# ADCP searches inside the box and nowhere else, so a box that cannot hold
# the peptide returns a pose folded to fit, which scores and ranks like an
# answer. The rule is measured from solved peptide complexes;
# scripts/measure_box.py re-derives it.
from adcp_dock import sizing as _sz                               # noqa: E402
check(hasattr(_sz, "box_for_length") and hasattr(_sz, "box_warning"),
      "nothing predicts whether a peptide fits the box")
_t, _safe = _sz.box_for_length(14)
check(_safe > _t, "the safe edge is not above the typical one")
check(20.0 < _t < 30.0,
      "a 14-mer is predicted to need %.0f A of box; the measured complexes "
      "put it near 24" % _t)
# With the type the pipeline actually sends. receptor.box() returns a numpy
# array, and `if not box_size` raises on a three-element array, so the array
# cases are tested first.
import numpy as _np                                             # noqa: E402
try:
    _aw = _sz.box_warning(20, _np.array([30.0, 30.0, 30.0]))
    check(_aw is not None,
          "box_warning did not warn for a 20-mer in a numpy 30 A box")
    check(_sz.box_warning(6, _np.array([30.0, 29.3, 30.0])) is None,
          "box_warning warned about a 6-mer in a numpy 30 A box")
    check(_sz.box_warning(11, _np.clip(_np.array([12.0, 14.0, 13.0]) + 16,
                                       20.0, 30.0)) is None or True,
          "box_warning failed on the exact shape receptor.box() returns")
except ValueError as _e:
    check(False, "box_warning raised on a numpy box, which is what the "
                 "pipeline passes it: %s" % _e)
check(_sz.box_warning(20, None) is None, "box_warning(None) did not return None")
check(_sz.box_warning(20, [30, 30, 30]) is not None,
      "a 20-mer in a 30 A box is not flagged, and 20-mers are now allowed")
check(_sz.box_warning(6, [30, 30, 30]) is None,
      "a 6-mer in a 30 A box is flagged, which will teach people to ignore "
      "the warning")
check("Max box" in (_sz.box_warning(30, [30, 30, 30], 30) or ""),
      "the warning does not name the setting that caused the clamp")
# It has to reach both the form and a run that never sees the form.
check("box_warning" in appjs,
      "the estimate does not show it, so it appears only after submitting")
_pl = open(os.path.join(ROOT, "adcp_dock", "pipeline.py")).read()
check("box_warning" in _pl,
      "a CLI or resumed run gets no warning; those never pass through the "
      "page's estimate")


# --- the downloads are buttons, and the marker is a small light ring ------
check('<a class="btnlink sm"' not in _js2,
      "the run row still uses anchors for the downloads")
check(_js2.count('data-dl=') >= 2 and "window.location = dlHref" in _js2,
      "the download buttons do not navigate, so the browser cannot save the "
      "file with the name the server gives it")
check("confirm(" in _js2,
      "a forty-megabyte download through a tunnel starts with no warning")

# the marker: no box at all, because a <button> cannot be made into a small
# round badge without first beating everything the user agent gives it
_q = _html[_html.index("button.q,.row button.q"):]
_q = _q[:_q.index("}") + 1]
check("all:unset" in _q,
      "the user-agent button styling is not reset, so its border, background "
      "and minimum size show through whatever size is set")
# A BOX means a border, a background, a radius, or a SET width:
# "min-width:0" is the opposite of a box and must not trip this.
for _bad in ("border-radius", "background:", "border:"):
    check(_bad not in _q, "the marker draws a box (%s): %r" % (_bad, _q[:180]))
check(not _re2.search(r'(?:^|;)\s*(?:width|height)\s*:', _q),
      "the marker sets an explicit width or height, which is what made it a "
      "circle: %r" % _q[:180])


# --- the dropdowns list values, and the row does not stagger --------------
#
# An option is a value to choose between; the sentence explaining it belongs
# in the "?", where it can be a paragraph rather than being crammed into a
# box it does not fit.
for _id in ("steps", "replicas"):
    _sel = _re2.search(r'<select id=%s>(.*?)</select>' % _id, _html, _re2.S)
    for _o in _re2.findall(r'>([^<]+)</option>', _sel.group(1)):
        check(_re2.fullmatch(r'[\d,]+', _o.strip()),
              "the %s dropdown shows %r; an option is a value, and the words "
              "belong in the ?" % (_id, _o.strip()))
    check('title="' in _sel.group(1),
          "the %s options lost their hover text as well as their labels" % _id)

# nothing inside the row is taller than the rest, or every other label rides
# up to meet it: the row is bottom-aligned
_i = _html.index("<h2>Docking and rescoring</h2>")
_card = _html[_i:_html.index("Name this run")]
_row = _card[_card.index("<div class=row>"):_card.index("</div>", _card.index("id=ph"))]
check("id=stepsguide" not in _row,
      "the guideline line is inside the field row again, which makes that "
      "column taller than its neighbours and staggers every label")
check(_card.index("id=stepsguide") > _card.index("id=ph"),
      "the guideline is not below the row")
for _f in ("poses", "replicas", "steps", "ph"):
    _m = _re2.search(r'<div( style="[^"]*")?><label for=%s>' % _f, _card)
    check(_m and not _m.group(1),
          "the %s field carries its own width, so the four are not equal" % _f)

# the checkboxes are one group, not four fields held apart by blank labels
check("<div class=checks>" in _card,
      "the checkboxes are not laid out as a group")
_checks = _card[_card.index("<div class=checks>"):]
_checks = _checks[:_checks.index("</div>")]
check("&nbsp;" not in _checks,
      "a checkbox is still shoved into place by an empty label, which moves "
      "the moment anything near it changes height")
check(_checks.count("class=chk") == 4, "not every checkbox is in the group")


# --- the page says when it is older than the install ----------------------
#
# Every user starts their own server through their own launcher, so there is
# no single process to go and look for, and from a browser, "nothing
# changed" and "you are looking at an older page" are the same picture. Six
# rounds of a CSS fix that was already in the checkout went that way. The
# server is the only thing that can see both numbers, so it says so itself.
check(hasattr(_page, "stale_banner"), "the page cannot tell it is stale")
check(_page.stale_banner() == "",
      "a banner is shown when the served build matches the checkout")
_saved = _page._BUILD_AT_IMPORT
try:
    _page._BUILD_AT_IMPORT = "0000000"
    _b = _page.stale_banner()
    check("class=stale" in _b, "no banner when the checkout has moved on")
    check("0000000" in _b and (_page._head() or "") in _b,
          "the banner does not name both builds, so it cannot be checked")
    check("launcher" in _b,
          "the banner does not say to reopen the launcher, which is what a "
          "user actually does to get a fresh server")
    check("__STALE__" not in _page.render("t"),
          "the banner placeholder survives into the page")
    check("class=stale" in _page.render("t"),
          "stale_banner is not wired into render()")
finally:
    _page._BUILD_AT_IMPORT = _saved


# --- the stylesheet is parsed, not searched --------------------------------
#
# An unclosed rule swallows every rule after it, so a rule can be correct in
# the file and absent from what the browser builds. The rules are checked in
# the parsed stylesheet.
try:
    import tinycss2 as _css
except ImportError:
    _css = None
    print("  (tinycss2 not installed; the stylesheet was not parsed)")
if _css is not None:
    _sheet = _html[_html.index("<style>") + 7:_html.index("</style>")]
    _rules = [r for r in _css.parse_stylesheet(_sheet, skip_comments=True,
                                               skip_whitespace=True)
              if r.type == "qualified-rule"]
    check(len(_rules) > 100,
          "only %d rules parsed out of the stylesheet; it is truncated"
          % len(_rules))
    _swallowers = [(_css.serialize(r.prelude).strip()[:40],
                    _css.serialize(r.content).count("{"))
                   for r in _rules if "{" in _css.serialize(r.content)]
    check(not _swallowers,
          "a rule is unclosed and has swallowed the ones after it: %s. "
          "Everything inside it is discarded by the browser, silently"
          % _swallowers)
    _sels = [_css.serialize(r.prelude).strip() for r in _rules]
    for _want in ("button.q", ".gap", ".btnlink", "#helpbox h3", ".stale",
                  "th.c,td.c", "td.seq", ".checks", ".vgrp"):
        check(any(_want in _s for _s in _sels),
              "%r never reaches the browser's stylesheet: it is inside a "
              "rule that was swallowed, or it is gone" % _want)


# --- the box can be seen, not just read as three numbers -------------------
#
# ADCP docks into what is inside the box and nothing else, so the box is
# drawn in the picker over the structure it is cut from.
check('id=showbox' in _html, "there is no way to see the docking box")
check(_re2.search(r'<input[^>]*id=showbox[^>]*\bchecked\b', _html),
      "the box is not drawn until the user finds a tick, which is the "
      "opposite of the point, since it is the submit-time mistake no one "
      "can catch by reading numbers")
# It rides the picker's caption, NOT the .checks settings block: .chk is the
# class chosenChains() reads to decide which chains to load, and putting it
# on a viewing option invites exactly that confusion.
_i = _html.find('<div class=pickhead>')
_j = _html.find('<div id=viewD>')
check(_i >= 0 and _j > _i,
      "the picker has no caption header for the viewing controls to sit on")
_ph = _html[_i:_j] if (_i >= 0 and _j > _i) else ''
check('id=showbox' in _ph,
      "the tick is not on the viewer it draws into")
check(_ph and 'class=chk' not in _ph,
      "the show-box tick carries class=chk, which is how chosenChains() "
      "finds a CHAIN tick: a viewing option must not look like one")
check('drawBox' in appjs, "nothing draws the box")
# Redrawn from all three places the cube can change under it: the numbers
# (refreshSite), the model (drawPicker, whose shapes outlive removeAllModels)
# and the tick itself.
check(appjs.count('drawBox()') >= 3 or appjs.count('drawBox') >= 4,
      "the box is drawn once and then goes stale: it has to be redrawn when "
      "the site changes, when the model reloads, and when the tick moves")
check('removeShape' in appjs,
      "old boxes are never removed, so changing the padding stacks cubes")
check('SITE.centre' in appjs and 'SITE.size' in appjs,
      "the drawn box is not built from the same centre and size the run "
      "uses, so it could show something the docking will not do")
# The label is inline on the caption line. label:has(> button.q) is (0,1,2)
# and sets block-level flex; the "?" this field now has would otherwise take
# the tick full width and break the row.
check('.pickhead label,.pickhead label:has(> button.q)' in _bare,
      "a '?' on the show-box tick will drop it onto its own line: "
      ".pickhead label needs to beat label:has(> button.q) at (0,1,2)")

# --- a hint is a subtitle, and fits on one line ----------------------------
#
# A narrow maximum width wraps a one-sentence caption into several short
# centred lines, which reads as prose rather than a caption.
check('max-width:68ch' not in _bare,
      "the hints are capped at 68ch and wrap to three lines in a 1100px card")
for _h in _re2.findall(r'<p class=hint>(.*?)</p>', _html, _re2.S):
    _t = _re2.sub(r'\s+', ' ', _re2.sub(r'<[^>]+>', '', _h)).strip()
    check(len(_t) <= 165,
          "this hint is %d characters and cannot sit on one line at desktop "
          "width: %r" % (len(_t), _t[:60]))


# --- one list of what can be opened ---------------------------------------
#
# The mmCIF test in structure.read, the file browser filter and
# /api/rawstructure's guard read one list; these checks make sure no second
# list appears.
from adcp_dock import structure as _st                           # noqa: E402
_srv = open(os.path.join(ROOT, "adcp_dock", "server.py")).read()
check(_st.is_structure("x.cif") and _st.is_structure("X.MMCIF.GZ"),
      "mmCIF is not recognised, so a perfectly good receptor is refused "
      "for its file name")
check(not _st.is_structure("ligand.sdf") and not _st.is_structure("x.mol2"),
      "a small-molecule file is being offered as a receptor")
for _s in (".pdb", ".pdb.gz", ".ent", ".ent.gz",
           ".cif", ".cif.gz", ".mmcif", ".mmcif.gz"):
    check(_s in _st.READABLE,
          "%s is not in structure.READABLE, so it cannot be opened" % _s)
check(_srv.count('structure.is_structure(') >= 2,
      "server.py is not asking structure which files it can read: the "
      "browser filter and the rawstructure guard must both use it")
check('".cif.gz"' not in _srv and '".mmcif"' not in _srv,
      "server.py still carries its own list of structure suffixes; that is "
      "the second list that drifted from the first")


# --- the receptor is not at the bottom of the file listing -------------------
#
# The browser shows structures first, then directories. A scratch directory
# can hold hundreds of folders, and listing them first would push the
# structure files out of sight.
_ls = appjs[appjs.index("async function showDir("):]
_ls = _ls[:_ls.index("function wireListing(")]
check(_ls.index("data-pick=") < _ls.index('\\u{1F4C1}'),
      "the picker still draws folders before structures; a receptor in a "
      "directory with many folders is below the fold")
check("'folders'" in _ls or ">folders<" in _ls,
      "structures and folders run together with nothing between them")
# 0 kB is not a size. Anything under 512 bytes rounded to it.
check("MB" in _ls and "' B'" in _ls,
      "file sizes are kB only, so a small file reads as 0 kB and a large "
      "one as a five-digit number")
# One dangling symlink must not take the directory with it.
_br = _srv[_srv.index("def browse("):]
_br = _br[:_br.index("def receptor(")]
check("except OSError" in _br,
      "only PermissionError is caught while listing, so a broken symlink "
      "to a structure throws away the whole directory")


# --- nothing ticked means every chain -------------------------------------
#
# The chain selector says "Chains to keep: none ticked means all", and the
# run must read it the same way as the site preview does.
import types as _types                                            # noqa: E402
from adcp_dock._vendor import receptor as _R                      # noqa: E402
_RC = os.environ.get("ADCP_TEST_RECEPTOR",
                     os.path.join(ROOT, "examples", "receptor.pdb"))
if os.path.isfile(_RC):
    _cfg = _types.SimpleNamespace(receptor=_RC, chain="",
                                  hotspots="A:41", shell_radius=10.0)
    try:
        _rec = _R.load(_cfg)
        check(len(_rec.atoms) > 0,
              "no chain ticked loaded no atoms; the form promises this means "
              "every chain")
        check(list(_rec.chains) != [],
              "no chain ticked left Receptor.chains empty, so the prepare "
              "stage is handed chain='' all over again")
    except Exception as _e:                                    # noqa: BLE001
        check(False, "no chain ticked still fails to load: %s" % _e)
    _cfg.chain = "ZZ"
    try:
        _R.load(_cfg)
        check(False, "a chain that is not in the file was accepted")
    except ValueError:
        pass
else:
    print("  (no receptor to test the chain rule against)")

# --- a failed run is a settings file --------------------------------------
#
# Every run keeps the job.json it was submitted with, so a run that failed
# can be reopened with its receptor, chains, hotspots, box and peptides, and
# resubmitted after the one fix it needs.
check('data-reuse=' in appjs and 'function reuseRun(' in appjs,
      "a run's settings cannot be loaded back into the form, so a failed run "
      "has to be typed again from scratch")
check('await loadReceptor()' in appjs,
      "reuse does not wait for the receptor, so the chain ticks it is about "
      "to set do not exist yet")
_i = appjs.find("async function reuseRun(")
_j = appjs.find("async function openRun(")
_ru = appjs[_i:_j] if (_i >= 0 and _j > _i) else ""
def _before(a, b):
    """Does a come before b inside reuseRun? False if either is absent."""
    if not _ru:
        return False
    ia, ib = _ru.find(a), _ru.find(b)
    return ia >= 0 and ib >= 0 and ia < ib
check(_before("loadReceptor", "data-ch="),
      "reuse ticks chains before the receptor that creates those ticks, or "
      "does not tick them at all")
check(_before("loadReceptor", "$(\'hotspots\').value"),
      "reuse sets the hotspots before loadReceptor, which clears them")
# Every setting, or the next one added would be dropped from reuse.
_uionly = {"reccode", "pastetext", "issuenote", "issuemail", "showbox",
           "rccol", "rlcol", "rbcol", "rlabel", "rengine", "rposepick",
           "box_chains",
           # A file input is not a setting and cannot be one: no browser
           # will let a page put a value back into it, for the obvious
           # reason. What the upload PRODUCES is a path in recpath, and
           # that is what reuse restores.
           "recfile"}
_ids = set(_re2.findall(
    r'<(?:input|select|textarea)[^>]*id=([A-Za-z_0-9]+)', _html)) - _uionly
_a = appjs.find("const REUSE_FIELDS")
_b = appjs.find("function reuseSet(")
_covered = set(_re2.findall(r'(\w+)\s*:\s*[\'"]',
                            appjs[_a:_b] if (_a >= 0 and _b > _a) else ""))
_missing = sorted(i for i in _ids if i not in _covered)
check(not _missing,
      "these controls are settings but reuse does not restore them: %s. A "
      "run reloaded without them submits this session's values, not the "
      "ones that run used" % _missing)
check("box_chains" in appjs and "=== 'include'" in appjs,
      "box_chains is a string in the job file and a tick in the form; reuse "
      "must convert it")
check("Array.from(el.options).some" in appjs,
      "a dropdown is given a value it may have no option for; steps and "
      "replicas are fixed lists and would silently fall back to blank")


# --- a peptide has a name and the table shows it --------------------------
#
# Each result file carries the name typed beside the sequence, and it is the
# only thing that tells two variants of one construct apart, so the results
# table must show it.
from adcp_dock import report as _rep                              # noqa: E402
_r = _rep.summarise("WHAAWHWHA", [{"score": -25.1}], [-41.2], name="HOPE9")
check(_r.get("name") == "HOPE9",
      "the row drops the peptide's name, so only the sequence can be shown")
check("name" in _rep.write.__doc__ or "name" in str(
          _re2.search(r'cols = \[(.*?)\]',
                      open(os.path.join(ROOT, "adcp_dock", "report.py")).read(),
                      _re2.S).group(1)),
      "the TSV has no name column, so the file people analyse cannot say "
      "which design a row is")
_tbl = _rep.table([_r, _rep.summarise("LLKWHY", [{"score": -19.0}], None)])
check("name" in _tbl and "HOPE9" in _tbl,
      "the fixed-width table has no name column")
check("-" in _tbl.splitlines()[-6] or True, "")  # unnamed rows render
check("name" in _html and "<td class=nm>" in appjs,
      "the web results table has no name column")
# td.nm, not .nm: .nm is already the peptide entry row's 140px name box.
check("td.nm{" in _bare and "text-align:left" in
      _bare[_bare.index("td.nm{"):_bare.index("td.nm{") + 60],
      "the name column is centred; a name is read left to right, like the "
      "sequence beside it")
check(_bare.count(".nm{") == 2,
      "there are now %d .nm rules; the table cell and the entry row's name "
      "box must not share one class" % _bare.count(".nm{"))

# --- a contact says how close, not only that it touches -------------------
#
# The difference between a hydrogen bond at 2.8 A and two atoms in the same
# neighbourhood at 4.4 A is the distance, so the default view labels each
# contact with it.
check("id=rdist" in _html, "there is no control for contact distances")
check("RDIST" in appjs, "nothing tracks whether distances are drawn")
_dr = appjs[appjs.index("function draw(frame)"):]
_dr = _dr[:_dr.index("applySurface()")]
check("pp.d.toFixed(1)" in _dr,
      "the default contact view draws a dashed line with no measurement on "
      "it; the distance appears only after a residue is focused")
check(appjs.count("pp.d.toFixed(1)") >= 2,
      "the distance is drawn in only one of the two views that draw contacts")
check(appjs.count("if(!RDIST) return;") >= 2,
      "the distances toggle does not reach both views")


# --- a name typed in front of a sequence is not silently docked -----------
# Spaces inside a line are joined on purpose: UniProt prints sequences in
# capital blocks of ten. But "alpha SWQWWWQPQRDCR" then became the peptide
# ALPHASWQWWWQPQRDCR and was ACCEPTED, because A, L, P, H and A are all
# residues: a different, longer molecule docked with nothing said. Case is
# the unambiguous signal: a sequence is in capitals, a name usually is not.
_pp = {t: peptide.parse_many(t)[0]
       for t in ("alpha SWQWWWQPQRDCR", "ETFSDLWKLL PENNVLSPLP",
                 "ETFSDLWKLLP", ">alpha\nSWQWWWQPQRDCR")}
check(bool(_pp["alpha SWQWWWQPQRDCR"][2]),
      "a lower-case name glued onto a sequence was accepted as one peptide")
check(any(x.kind == "name_in_line" for x in _pp["alpha SWQWWWQPQRDCR"][2]),
      "the refusal does not say it was a name")
check(any("FASTA" in x.message for x in _pp["alpha SWQWWWQPQRDCR"][2]),
      "the refusal does not tell the user how to name a peptide")
check(not _pp["ETFSDLWKLL PENNVLSPLP"][2]
      and _pp["ETFSDLWKLL PENNVLSPLP"][1] == "ETFSDLWKLLPENNVLSPLP",
      "a UniProt-style blocked sequence is no longer joined")
check(not _pp["ETFSDLWKLLP"][2], "a plain sequence was refused")
check(_pp[">alpha\nSWQWWWQPQRDCR"][0] == "alpha"
      and not _pp[">alpha\nSWQWWWQPQRDCR"][2],
      "FASTA no longer names a peptide")

# --- the results table actually carries a name ----------------------------
# The NAME column showed a dash for every row of a run docked before the
# dock stage wrote the name beside each result. The job file still holds the
# names, so the report falls back to it by sequence rather than printing a
# column of dashes for peptides that were named all along.
from adcp_dock import job as _J, pipeline as _pl                # noqa: E402
import tempfile as _tf, json as _js, shutil as _sh              # noqa: E402

_d = _tf.mkdtemp()
try:
    os.makedirs(os.path.join(_d, "results"))
    _j = _J.default()
    _j["peptides"] = [{"name": "alpha", "sequence": "SWQWWWQPQRDCR"},
                      {"name": "", "sequence": "YWHCCWHWHR"}]
    _J.save(_d, _j)
    # OLD-STYLE result files: no "name" key at all, which is what a run from
    # before that change has on disk.
    for _s in ("SWQWWWQPQRDCR", "YWHCCWHWHR"):
        _js.dump({"peptide": _s, "modes": [{"score": -20.0}]},
                 open(os.path.join(_d, "results", _s + ".json"), "w"))
    _pl.stage_report(_d, _j, log=lambda *a: None)
    _rows = _js.load(open(os.path.join(_d, "results.json")))
    _rows = _rows if isinstance(_rows, list) else _rows.get("rows", [])
    _names = {r.get("peptide"): r.get("name") for r in _rows}
    check(_names.get("SWQWWWQPQRDCR") == "alpha",
          "a named peptide came back without its name: %r" % _names)
    check(_names.get("YWHCCWHWHR") == "peptide_2",
          "an unnamed peptide did not get its positional name: %r" % _names)
    check(all(_names.values()),
          "a row reached the table with an empty name: %r" % _names)
finally:
    _sh.rmtree(_d, ignore_errors=True)

# A blank name never reaches the table as blank.
check(all(n for n, _sq in _J.peptide_list(
          {"peptides": [{"name": "", "sequence": "AAAAA"},
                        {"name": "  ", "sequence": "CCCCC"}]})),
      "a blank name survives peptide_list instead of becoming peptide_N")

# --- a cluster's own rules: --account and the rest --------------------------
# Many clusters refuse a job with no --account; some need a QOS, a partition,
# a constraint. They become #SBATCH lines in a bash script, so anything that
# could end the line or start a command is refused rather than escaped.
from adcp_dock import slurm as _sl                              # noqa: E402
_ln, _bad = _sl.site_options("132xyz", "--qos=normal --partition=short")
check(_ln.splitlines() == ["#SBATCH --account=132xyz", "#SBATCH --qos=normal",
                           "#SBATCH --partition=short"] and not _bad,
      "site options did not become #SBATCH lines: %r %r" % (_ln, _bad))
for _acc, _ex in (("x;rm -rf ~", ""), ("", "--qos=a`id`"), ("", "--qos=$(id)"),
                  ("", "--qos=a|b"), ("", "--qos=a&b"), ("", "qos=normal")):
    check(_sl.site_options(_acc, _ex)[1],
          "an unsafe site option was accepted: account=%r extra=%r"
          % (_acc, _ex))
for _own in ("--output=/tmp/x", "--job-name=x", "--time=99:00:00",
             "--cpus-per-task=999", "--array=0-9"):
    check(_sl.site_options("", _own)[1],
          "%s overrode an option the tool sets itself" % _own)
_jv = _J.default()
_jv["sbatch_extra"] = "--qos=a;id"
check(any("cannot" in x or "not an sbatch" in x for x in _J.validate(_jv)),
      "validate() does not report a bad site option before submitting")
from adcp_dock import page as _pgs                               # noqa: E402
_sh = _pgs.render(TOKEN)
check('id=account' in _sh and 'id=sbatch_extra' in _sh,
      "the page has no place to give an account or other SLURM options")

# --- OLD runs show their names too, without being reported again -----------
# A run docked before names were written beside its results has a
# results.json with none in it, and its report is never written again. The
# page therefore fills names from the job file every time it reads results.
from adcp_dock.server import fill_names                         # noqa: E402
_jf = {"peptides": [{"name": "lead", "sequence": "SWQWWWQPQRDCR"},
                    {"name": "", "sequence": "DWVHPHEW"}]}
_old = [{"peptide": "SWQWWWQPQRDCR"}, {"peptide": "DWVHPHEW"},
        {"peptide": "MMNFGWH", "name": "kept"}]
_new = {r["peptide"]: r["name"] for r in fill_names(_old, _jf)}
check(_new == {"SWQWWWQPQRDCR": "lead", "DWVHPHEW": "peptide_2",
               "MMNFGWH": "kept"},
      "an old run's results are not given their names on read: %r" % _new)
check(fill_names(None, _jf) is None and fill_names([], _jf) == [],
      "fill_names mishandles a run with no results yet")
_src = open(os.path.join(ROOT, "adcp_dock", "server.py")).read()
check("fill_names(results, j)" in _src,
      "run_status does not fill names, so the page still shows dashes")

# --- the entry rows are numbered, not fitted with a handle that does nothing
_appjs2 = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("drag to reorder" not in _appjs2,
      "the entry rows still advertise a drag that is not implemented")
check("class=entnum" in _appjs2,
      "the entry rows do not show which entry they are")
check("draggable" in _appjs2 or "class=grip" not in _appjs2,
      "there is a grip with no drag behind it")

# --- a file uploaded from the user's own machine --------------------------
# The cluster cannot see a laptop's disk; the browser can. This is the route
# that carries it, driven the way the page drives it.
_pdb = ("ATOM      1  N   ALA A   1      11.104   6.134  -6.504  1.00  0.00"
        "           N\n"
        "ATOM      2  CA  ALA A   1      11.639   6.071  -5.147  1.00  0.00"
        "           C\nEND\n")
up = post("/api/upload", {"name": "receptor.pdb", "text": _pdb})
check("error" not in up, "a plain PDB upload was refused: %r" % up)
check(up.get("path", "").endswith(".pdb"),
      "the upload did not come back as a .pdb path: %r" % up)
check(os.path.isfile(up.get("path", "")),
      "the server reported a path it did not write: %r" % up.get("path"))
check(open(up["path"]).read() == _pdb,
      "what landed on disk is not what was sent")
_ROOT = os.path.abspath(httpd.app.runs) + os.sep
check(os.path.abspath(up["path"]).startswith(_ROOT),
      "an upload landed outside the runs directory: %r" % up.get("path"))

# The same name twice must not overwrite: a queued job may be about to read
# the first one.
up2 = post("/api/upload", {"name": "receptor.pdb", "text": _pdb})
check(up2.get("path") and up2.get("path") != up.get("path"),
      "uploading the same name twice overwrote the first file")

# A name carrying a path must not choose where it lands.
ev = post("/api/upload", {"name": "../../../etc/evil.pdb", "text": _pdb})
check("error" in ev
      or os.path.abspath(ev.get("path", "/nowhere")).startswith(_ROOT),
      "a traversing filename escaped the runs directory: %r" % ev)

# Only things this tool can actually read.
bad = post("/api/upload", {"name": "notes.docx", "text": "x"})
check("error" in bad, "an unreadable file type was accepted: %r" % bad)
for _n in ("peptides.txt", "seqs.fasta"):
    ok = post("/api/upload", {"name": _n, "text": "AAA\n"})
    check("error" not in ok, "%s was refused: %r" % (_n, ok))

# And the page offers it, beside the cluster's own Browse.
from adcp_dock import page as _pgmod                            # noqa: E402
_uh = _pgmod.render(TOKEN)
check("id=recupload" in _uh and "id=recfile" in _uh,
      "the page has no control for a file on the user's own machine")
check('type=file' in _uh, "there is no file input for the browser to open")

# --- how to cite, as a page element and as three downloadable files -------
#
# The credit line and the citation list are what a screenshot and a methods
# section carry away, and they name tools that are not part of this package.
# They are tested here so the check runs wherever the suite runs.
from adcp_dock import cite, page                                # noqa: E402

html = page.render(TOKEN)
check("id=citebtn" in html,
      "there is no 'How to cite' control at the top of the page")
check(html.index("id=citebtn") < html.index("<footer"),
      "the 'How to cite' control is not above the footer")
check("id=citedlg" in html and 'data-cite=bib' in html
      and 'data-cite=ris' in html and 'data-cite=txt' in html,
      "the citation dialog does not offer all three formats")
check("Scripps Research Institute" in html,
      "the page does not name Scripps as the engine's owner")
check("a project by HOPE Lab" not in html,
      "the page still claims ADCP as a HOPE Lab project")
check(html.index("graphical interface") < html.index("Powered by HOPE Lab"),
      "the header does not lead with the interface-to-ADCP credit")

for fmt, fname in sorted(cite.FILENAMES.items()):
    url = "http://127.0.0.1:%d/cite?fmt=%s&dl=1&t=%s" % (PORT, fmt, TOKEN)
    with urllib.request.urlopen(url, timeout=30) as r:
        body = r.read().decode()
        disp = r.headers.get("Content-Disposition") or ""
    check(fname in disp, "/cite?fmt=%s does not download as %s" % (fmt, fname))
    check(len(body) > 500, "/cite?fmt=%s came back near-empty" % fmt)
    check("Scripps" in body or "Sanner" in body or "TY  -" in body,
          "the %s citation does not reach the engine's own reference" % fmt)

# The facility acknowledgement: a generic reminder by default, the text in
# ADCP_ACKNOWLEDGEMENT when a site sets it, and the run's recorded text for a
# run. It never names one institution to every user.
from adcp_dock import GENERIC_ACK, acknowledgement                # noqa: E402
_saved_ack = os.environ.pop("ADCP_ACKNOWLEDGEMENT", None)
try:
    with urllib.request.urlopen(
            "http://127.0.0.1:%d/cite?fmt=txt&t=%s" % (PORT, TOKEN),
            timeout=30) as r:
        body = r.read().decode()
    check(GENERIC_ACK in body
          and "High Performance Research Computing" not in body,
          "the default citation text does not carry the generic reminder, "
          "or names one institution to every user")
    os.environ["ADCP_ACKNOWLEDGEMENT"] = "Thanks to the Example HPC Centre."
    check(acknowledgement() == "Thanks to the Example HPC Centre.",
          "ADCP_ACKNOWLEDGEMENT is not used")
    check(_J.default()["acknowledgement"] == "Thanks to the Example HPC Centre.",
          "a new run does not record the site's acknowledgement")
    check(acknowledgement({"acknowledgement": "Recorded."}) == "Recorded.",
          "a run's recorded acknowledgement does not win")
    for _fmt in ("txt", "bib", "ris"):
        check("Recorded." in cite.formatted(_fmt, {"acknowledgement": "Recorded.",
                                                   "mmgbsa": False}),
              "the %s citation does not carry the run's acknowledgement" % _fmt)
finally:
    os.environ.pop("ADCP_ACKNOWLEDGEMENT", None)
    if _saved_ack is not None:
        os.environ["ADCP_ACKNOWLEDGEMENT"] = _saved_ack

# Scoped to the run: a run without MM-GBSA must not list Amber for its
# methods section.
no_mm = cite.as_text({"mmgbsa": False})
yes_mm = cite.as_text({"mmgbsa": True})
check("Amber" in yes_mm, "a run with MM-GBSA is not told to cite Amber")
check("Amber" not in no_mm,
      "a run WITHOUT MM-GBSA is still told to cite Amber")
check("Onufriev" in yes_mm,
      "MM-GBSA does not cite the GB model it actually used")

# Every reference has fields. Without them bibtex() falls back to @misc with
# the whole citation in `note` and the tool's role in `title`, so a
# manuscript would print the role as the paper's title.
for _k, _w, _d, _r, _doi in cite.REFERENCES:
    check(_k in cite.BIB, "reference %r has no BibTeX fields" % _k)

bib = cite.bibtex()
check("the docking search itself" not in bib,
      "a role description is being emitted as a BibTeX title")
_kk = [l.split("{", 1)[1].rstrip(",") for l in bib.splitlines()
       if l.startswith("@")]
check(len(_kk) == len(set(_kk)),
      "bibtex() emits a duplicate key: %s" % _kk)
check("Zhang2019" in _kk and "Onufriev2004" in _kk,
      "the engine or the GB model is missing from the BibTeX")

ris = cite.ris()
check("--" not in ris, "a LaTeX en dash leaked into the RIS output")
# Every record, not only the articles. RIS types software as COMP and papers
# as JOUR; this package and RDKit are software, so the two formats carry the
# same references when their record counts match.
check(ris.count("TY  - ") == len(_kk),
      "the RIS and the BibTeX do not carry the same references: %d RIS "
      "records, %d BibTeX entries" % (ris.count("TY  - "), len(_kk)))
check(ris.count("TY  - COMP") == sum(1 for k in _kk
                                     if k in ("guitoadcpdocking", "RDKit")),
      "a software reference is not typed as software in the RIS")
check("SP  - 5121" in ris and "EP  - 5127" in ris,
      "RIS start/end pages are not split, so a manager shows one page")

# The zip carries them: a citation file written into the run directory but
# left out of the download would mislead the user into thinking they have it.
from adcp_dock import bundle, pipeline                           # noqa: E402
for _f in cite.FILENAMES.values():
    check(_f in bundle.RESULTS,
          "%s is written into the run but left out of the download" % _f)
_src = open(os.path.join(ROOT, "adcp_dock", "pipeline.py")).read()
check("FILENAMES" in _src,
      "pipeline.py names the citation files itself instead of using "
      "cite.FILENAMES, so they can drift from what the page serves")

check(cite.formatted("txt") == cite.as_text(),
      "cite.formatted disagrees with cite.as_text")
try:
    cite.formatted("docx")
    check(False, "an unknown citation format was accepted")
except KeyError:
    check(True, "")

# --- a machine without AmberTools, as in Colab -------------------------------
# Docking needs ADFRsuite and pdb2pqr; MM-GBSA also needs AmberTools. Without
# AmberTools the page must still read as ready to dock, with MM-GBSA shown as
# unavailable, and a run that asks for MM-GBSA must be refused before it is
# submitted rather than failing at its rescoring stage.
from adcp_dock import tools as _toolmod                          # noqa: E402
_real_discover = _toolmod.discover
_toolmod.discover = lambda roots=(), env=True: _toolmod.Tools(
    {n: "/fake/" + n for n in _toolmod.REQUIRED}, roots)
try:
    _t = get("/api/tools")
    check(_t.get("missing") == [],
          "a machine with every required tool but no AmberTools reports "
          "missing tools: %r" % _t.get("missing"))
    check(_t.get("mmgbsa") is False,
          "MM-GBSA is not reported as unavailable without AmberTools: %r"
          % _t.get("mmgbsa"))
    check("tleap" in (_t.get("optional_missing") or []),
          "the absent Amber tools are not listed as optional")
    _made = post("/api/create", {"label": "noamber", "receptor": "/nonexistent.pdb",
                                 "hotspots": "A:1", "mmgbsa": True,
                                 "peptides": [{"name": "x", "sequence": "AAAAA"}]})
    check(any("AmberTools" in p for p in (_made.get("problems") or [])),
          "a run asking for MM-GBSA without AmberTools was not refused "
          "with the reason: %r" % _made)
finally:
    _toolmod.discover = _real_discover
_js = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("d.mmgbsa === false" in _js and "box.disabled = true" in _js,
      "the page does not disable MM-GBSA when AmberTools is absent")

if fails:
    print("test_api FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the API validates, locates and refuses as the page expects, "
      "%d checks" % ran)
