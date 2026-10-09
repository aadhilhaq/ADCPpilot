"""The Folders view of the runs: the listing, the files it sends, and every
way of asking it for something outside the runs folder.

By Aadhil Haq

The view walks the runs folder over HTTP, so it is tested over HTTP, with the
token, exactly as the page asks. A runs folder is built here with a run in
it, a subfolder, files of each kind the page shows, and links that lead out:
a folder link, a file link, and a link to a folder inside. Then the server is
asked for paths outside the runs folder in every form a query string can
carry one, and has to refuse each.

    python3 tests/test_folders.py
"""
import json
import os
import shutil
import sys
import tempfile
import threading
import urllib.error
import urllib.parse
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import server                                      # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


# --- a runs folder, and things outside it ----------------------------------
TMP = tempfile.mkdtemp(prefix="adcp_folders_")
RUNS = os.path.join(TMP, "runs")
RUN = os.path.join(RUNS, "p53_screen_20261008_101500")
POSE_DIR = os.path.join(RUN, "docking", "ETFSDLWKLLPEN")
os.makedirs(POSE_DIR)
os.makedirs(os.path.join(RUN, "logs"))
os.makedirs(os.path.join(RUNS, "notes"))
with open(os.path.join(RUN, "job.json"), "w") as fh:
    json.dump({"label": "p53 screen", "peptides": []}, fh)
with open(os.path.join(RUN, "results.tsv"), "w") as fh:
    fh.write("rank\tpeptide\n1\tETFSDLWKLLPEN\n")
POSE = ("ATOM      1  CA  GLU A   1      11.104  13.207   2.100  1.00  0.00"
        "           C\nEND\n")
with open(os.path.join(POSE_DIR, "out_ranked_1.pdb"), "w") as fh:
    fh.write(POSE)
with open(os.path.join(RUN, "logs", "dock_0.out"), "w") as fh:
    fh.write("docking ETFSDLWKLLPEN\n")
with open(os.path.join(RUN, "receptor.trg"), "wb") as fh:
    fh.write(b"PK\x03\x04 a zip, which the page does not show")
# A page someone saved into the runs folder: it must never be served as a
# page, or it would run inside this one.
with open(os.path.join(RUNS, "notes", "page.html"), "w") as fh:
    fh.write("<script>alert(1)</script>")
with open(os.path.join(RUNS, "notes", ".hidden"), "w") as fh:
    fh.write("hidden")

OUTSIDE = os.path.join(TMP, "outside")
os.makedirs(OUTSIDE)
SECRET = os.path.join(OUTSIDE, "secret.txt")
with open(SECRET, "w") as fh:
    fh.write("not for the page")
os.symlink(OUTSIDE, os.path.join(RUNS, "escape"))              # a folder link out
os.symlink(SECRET, os.path.join(RUN, "leak.txt"))              # a file link out
os.symlink(RUN, os.path.join(RUNS, "latest"))                  # a link that stays in

httpd = server.build(host="127.0.0.1", port=0, runs=RUNS)
PORT = httpd.server_address[1]
TOKEN = httpd.app.token
threading.Thread(target=httpd.serve_forever, daemon=True).start()
# A cluster that sets a web proxy must not have 127.0.0.1 sent to it.
OPENER = urllib.request.build_opener(urllib.request.ProxyHandler({}))


def raw(endpoint, token=True, **q):
    """(status, headers, body bytes) for a GET."""
    if token:
        q["t"] = TOKEN
    url = "http://127.0.0.1:%d%s?%s" % (PORT, endpoint,
                                        urllib.parse.urlencode(q))
    try:
        with OPENER.open(url, timeout=30) as r:
            return r.status, r.headers, r.read()
    except urllib.error.HTTPError as exc:
        return exc.code, exc.headers, exc.read()


def folder(path=""):
    code, _h, body = raw("/api/folder", path=path)
    check(code == 200, "/api/folder answered %s for %r" % (code, path))
    return json.loads(body.decode())


# --- the token gates it, as it gates every other endpoint -------------------
for ep in ("/api/folder", "/api/runfile"):
    code, _h, _b = raw(ep, token=False, path="")
    check(code == 403, "%s answered %s with no token, not 403" % (ep, code))
code, _h, _b = raw("/api/folder", token=False, path="", t="not-the-token")
check(code == 403, "/api/folder answered %s to a wrong token" % code)

# --- the top of the runs folder ---------------------------------------------
top = folder("")
check(top.get("rel") == "" and top.get("parent") is None,
      "the top of the runs folder has a parent to go up to: %r" % top)
check(top.get("root") == RUNS, "the root reported is %r" % top.get("root"))
names = [d["name"] for d in top.get("dirs", [])]
check("escape" not in names,
      "a link leading out of the runs folder is listed: %r" % names)
check(os.path.basename(RUN) in names and "notes" in names,
      "the folders in the runs folder are not all listed: %r" % names)
check("latest" in names,
      "a link to a folder inside the runs folder was left out: %r" % names)
run_row = next((d for d in top.get("dirs", [])
                if d["name"] == os.path.basename(RUN)), {})
check(run_row.get("run") == "p53 screen",
      "a run folder is not marked as a run, with its label: %r" % run_row)
check(run_row.get("path") == RUN,
      "a run folder's path is %r, not the path the runs table uses"
      % run_row.get("path"))
check(isinstance(run_row.get("mtime"), (int, float)),
      "a folder has no modification time: %r" % run_row)
check(run_row.get("items") == len([n for n in os.listdir(RUN)
                                   if not n.startswith(".")]),
      "a folder's item count is %r" % run_row.get("items"))
notes_row = next((d for d in top["dirs"] if d["name"] == "notes"), {})
check(notes_row.get("run") is None, "a plain folder is marked as a run")
check(top.get("crumbs") == [{"name": "runs", "rel": ""}],
      "the trail at the top is %r" % top.get("crumbs"))

# --- inside a run -----------------------------------------------------------
inrun = folder(os.path.basename(RUN))
check(inrun.get("run") == "p53 screen",
      "a run's own folder does not say it is a run: %r" % inrun.get("run"))
check(inrun.get("path") == RUN, "a run's folder path is %r" % inrun.get("path"))
files = {f["name"]: f for f in inrun.get("files", [])}
check("leak.txt" not in files,
      "a file link leading out of the runs folder is listed: %r" % sorted(files))
check(files.get("results.tsv", {}).get("view") == "report",
      "a run's results file does not open as its results: %r"
      % files.get("results.tsv"))
check(files.get("job.json", {}).get("view") == "text",
      "job.json is not shown as text: %r" % files.get("job.json"))
check(files.get("receptor.trg", {}).get("view") == "",
      "a target file (a zip) is offered to be shown: %r"
      % files.get("receptor.trg"))
check(all(isinstance(f.get("size"), int) and isinstance(f.get("mtime"),
                                                         (int, float))
          for f in files.values()),
      "a file has no size or no modification time: %r" % files)
check([c["name"] for c in inrun["crumbs"]] == ["runs", os.path.basename(RUN)],
      "the trail inside a run is %r" % inrun["crumbs"])
check(inrun.get("parent") == "",
      "going up from a run does not lead to the top: %r" % inrun.get("parent"))

# --- a subfolder of a run ---------------------------------------------------
sub = folder(os.path.join(os.path.basename(RUN), "docking", "ETFSDLWKLLPEN"))
check(sub.get("run") is None, "a run's subfolder is called a run")
pose = next((f for f in sub.get("files", []) if f["name"] == "out_ranked_1.pdb"),
            {})
check(pose.get("view") == "structure",
      "a pose is not shown in the 3D viewer: %r" % pose)
check(sub.get("parent") == os.path.join(os.path.basename(RUN), "docking"),
      "going up from a subfolder leads to %r" % sub.get("parent"))
check(len(sub.get("crumbs", [])) == 4, "the trail has %r" % sub.get("crumbs"))
# a results file in a folder that is not a run is only text
notes = folder("notes")
check(not any(f["name"] == ".hidden" for f in notes.get("files", [])),
      "a hidden file is listed")

# The link that stays inside resolves to the run it points at.
via = folder("latest")
check(via.get("rel") == os.path.basename(RUN) and via.get("run") == "p53 screen",
      "a link inside the runs folder does not lead to the run: %r"
      % {k: via.get(k) for k in ("rel", "run", "error")})

# --- paths outside the runs folder, in every form ----------------------------
TRIES = [
    "..", "../", "../outside", "../outside/secret.txt", "../../..",
    os.path.basename(RUN) + "/../../outside",
    os.path.basename(RUN) + "/../..",
    "./../outside", "notes/../../outside",
    "escape", "escape/", "escape/secret.txt",
    os.path.basename(RUN) + "/leak.txt",
    OUTSIDE, SECRET, TMP, "/", "/etc", "/etc/passwd",
    RUNS + "/../outside", RUNS + "/escape",
    "//etc/passwd", "a\x00b", "notes/page.html\x00.txt",
]
for p in TRIES:
    got = folder(p)
    check("error" in got and not got.get("dirs") and not got.get("files"),
          "/api/folder listed something for %r: %r" % (p, got))
    code, h, body = raw("/api/runfile", path=p)
    check(code in (403, 404) and b"not for the page" not in body,
          "/api/runfile sent %s for %r: %r" % (code, p, body[:80]))
    code, h, body = raw("/api/runfile", path=p, dl="1")
    check(code in (403, 404) and b"not for the page" not in body,
          "/api/runfile with dl=1 sent %s for %r" % (code, p))
# The refusal says why, and names the runs folder, not the path asked for.
got = folder("../outside")
check("runs folder" in (got.get("error") or ""),
      "a path outside is refused without saying why: %r" % got)
code, _h, body = raw("/api/runfile", path="../outside/secret.txt")
check(code == 403, "a file outside is answered %s, not 403" % code)
# Encoded forms arrive decoded once by the query string; a second encoding
# is a name with % in it, which is not there.
for enc in ("%2e%2e", "%2e%2e%2foutside", "..%2foutside", "..%252f"):
    url = ("http://127.0.0.1:%d/api/folder?path=%s&t=%s"
           % (PORT, enc, TOKEN))
    with OPENER.open(url, timeout=30) as r:
        got = json.loads(r.read().decode())
    check("error" in got and not got.get("dirs"),
          "an encoded path %r listed something: %r" % (enc, got))
# A folder that is not there, and a file asked for as a folder.
check("error" in folder("no_such_folder"), "a missing folder listed something")
check("error" in folder(os.path.join(os.path.basename(RUN), "job.json")),
      "a file was listed as a folder")

# --- the files it sends -------------------------------------------------------
rel_pose = os.path.join(os.path.basename(RUN), "docking", "ETFSDLWKLLPEN",
                        "out_ranked_1.pdb")
code, h, body = raw("/api/runfile", path=rel_pose)
check(code == 200 and body.decode() == POSE,
      "a pose was not sent as it is on disk (%s)" % code)
check(h.get("Content-Type", "").startswith("text/plain"),
      "a file to show was sent as %r" % h.get("Content-Type"))
check(h.get("X-Content-Type-Options") == "nosniff",
      "a file to show can be sniffed into something else")
check("sandbox" in (h.get("Content-Security-Policy") or ""),
      "a file to show is not sandboxed")
check("attachment" not in (h.get("Content-Disposition") or ""),
      "a file to show was sent as a download")
code, h, body = raw("/api/runfile", path=rel_pose, dl="1")
check(code == 200 and "attachment" in (h.get("Content-Disposition") or "")
      and "out_ranked_1.pdb" in (h.get("Content-Disposition") or ""),
      "a download does not come as an attachment named for the file: %r"
      % h.get("Content-Disposition"))
check(h.get("Content-Type") == "application/octet-stream",
      "a download was sent as %r" % h.get("Content-Type"))
# An absolute path inside the runs folder is the same file.
code, _h, body = raw("/api/runfile", path=os.path.join(RUNS, rel_pose))
check(code == 200 and body.decode() == POSE,
      "an absolute path inside the runs folder was refused (%s)" % code)
# A page in the runs folder is text, never a page.
code, h, body = raw("/api/runfile", path="notes/page.html")
check(code == 200 and "html" not in h.get("Content-Type", ""),
      "an HTML file in the runs folder was served as %r"
      % h.get("Content-Type"))
# A folder is not a file.
code, _h, _b = raw("/api/runfile", path=os.path.basename(RUN))
check(code == 404, "a folder was sent as a file (%s)" % code)
# A file too large to show is refused for showing, and still downloads.
_was = server.App.VIEW_MAX
server.App.VIEW_MAX = 10
try:
    code, _h, body = raw("/api/runfile", path=rel_pose)
    check(code == 413, "a file over the limit was shown (%s)" % code)
    code, _h, body = raw("/api/runfile", path=rel_pose, dl="1")
    check(code == 200 and body.decode() == POSE,
          "a file over the showing limit could not be downloaded (%s)" % code)
    big = folder(os.path.join(os.path.basename(RUN), "docking",
                              "ETFSDLWKLLPEN"))
    check(all(f["view"] == "" for f in big["files"]),
          "a file over the limit is still offered to be shown")
finally:
    server.App.VIEW_MAX = _was
# A name that needs quoting in the header is made safe.
odd = os.path.join(RUNS, "notes", 'two "words" Å.txt')
with open(odd, "w") as fh:
    fh.write("x")
code, h, _b = raw("/api/runfile", path=os.path.relpath(odd, RUNS), dl="1")
cd = h.get("Content-Disposition") or ""
check(code == 200 and '"words"' not in cd and "filename*=UTF-8''" in cd,
      "a file name with quotes is not made safe in the header: %r" % cd)

# --- a long folder is cut, and says so -----------------------------------------
_was = server.App.FOLDER_MAX
server.App.FOLDER_MAX = 2
try:
    cut = folder(os.path.basename(RUN))
    check(len(cut["dirs"]) + len(cut["files"]) == 2 and cut.get("more", 0) > 0,
          "a folder over the limit was not cut, or did not say how many more: "
          "%r" % {k: cut.get(k) for k in ("more",)})
finally:
    server.App.FOLDER_MAX = _was

# --- the picker is unchanged -----------------------------------------------------
# /api/browse never lists every file: the Folders view asks for that through
# its own endpoint, which is held to the runs folder.
code, _h, body = raw("/api/browse", path=RUN, kind="all")
check("error" in json.loads(body.decode()),
      "/api/browse lists every file when asked: %r" % body[:120])
code, _h, body = raw("/api/browse", path=POSE_DIR)
pick = json.loads(body.decode())
check([f["name"] for f in pick.get("files", [])] == ["out_ranked_1.pdb"],
      "the receptor picker lists %r" % [f["name"] for f in pick.get("files", [])])
code, _h, body = raw("/api/browse", path=RUN)
check(not json.loads(body.decode()).get("files"),
      "the receptor picker lists files that are not structures")

# --- a runs folder reached through a link ---------------------------------------
LINKED = os.path.join(TMP, "runs_link")
os.symlink(RUNS, LINKED)
app2 = server.App(LINKED)
got = app2.runs_folder("")
check(not got.get("error") and got.get("root") == LINKED,
      "a runs folder reached through a link cannot be listed: %r"
      % got.get("error"))
inside = next((d for d in got.get("dirs", [])
               if d["name"] == os.path.basename(RUN)), {})
check(inside.get("path") == os.path.join(LINKED, os.path.basename(RUN)),
      "a run under a linked runs folder is named %r" % inside.get("path"))
check("error" in app2.runs_folder("../outside"),
      "a linked runs folder lets a path climb out of it")
check(app2.run_file(os.path.join(LINKED, "escape", "secret.txt"))[0] is None,
      "a linked runs folder sends a file through a link that leads out")

# --- the page offers both views ----------------------------------------------------
from adcp_dock import page                                         # noqa: E402
html = page.render("TOKEN")
for want in ('id=tab_runs', 'id=tab_folders', 'role=tablist',
             'id=folderspanel', 'id=runspanel', 'id=filedlg'):
    check(want in html, "the page has no %s" % want)
check(">All runs<" in html and ">Folders<" in html,
      "the two tabs are not called All runs and Folders")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
for want in ("/api/folder?path=", "/api/runfile?path=", "data-openrun",
             "function viewFile", "aria-selected"):
    check(want in appjs, "app.js has no %s" % want)

httpd.shutdown()
shutil.rmtree(TMP, ignore_errors=True)

if fails:
    print("test_folders FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the Folders view lists the runs folder, sends its files, and "
      "refuses every path out of it, %d checks" % ran)
