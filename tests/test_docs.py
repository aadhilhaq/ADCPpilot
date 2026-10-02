"""The built-in documentation: the pages, the pictures, and how they are served.

By Aadhil Haq

What is checked here is everything that can go wrong between the sidebar, the
page files and the interface: a page in the nav with no file, a link or a help
key pointing at a section that does not exist, a picture a page shows that was
never taken, a dash used as punctuation, and the routes themselves, which must
serve the docs without a token and refuse anything else.

    python3 tests/test_docs.py
"""
import os
import re
import struct
import sys
import threading
import urllib.error
import urllib.request
from http.server import ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import docs, docskit, page, server  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


# --- the pages themselves ---------------------------------------------------
# docskit's own audit: every page in the sidebar exists and has a heading and a
# lead, every link and every help key resolves to a section that is really
# there, every picture is in img/, and no page uses a dash as punctuation.
bad = docs.audit()
check(not bad, "the documentation does not audit clean:\n      "
      + "\n      ".join(bad[:12]))

# --- the copyright is LICENSE's, word for word ------------------------------
# The bar along the bottom of the interface and of every page states it, so it
# is held against the licence rather than typed twice and left to drift.
line3 = open(os.path.join(ROOT, "LICENSE")).read().splitlines()[2].strip()
check(docs.SITE.copyright == line3,
      "the docs' copyright is %r and line 3 of LICENSE is %r"
      % (docs.SITE.copyright, line3))

# --- every picture is a real PNG, large enough to read ----------------------
def png_width(path):
    """The width in pixels, from the IHDR chunk, or None if it is not a PNG."""
    with open(path, "rb") as fh:
        head = fh.read(24)
    if len(head) < 24 or head[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    return struct.unpack(">I", head[16:20])[0]


used = set()
for slug in docs.SITE.order:
    body = docskit.fragment(docs.SITE, slug) or ""
    used.update(re.findall(r'src="(img/[^"]+)"', body))
check(used, "no page shows a picture, so the screenshots are not being used")
for rel in sorted(used):
    full = os.path.join(docs.HERE, rel)
    w = png_width(full) if os.path.isfile(full) else None
    check(w is not None, "%s is missing or is not a PNG" % rel)
    check(w is None or w > 600,
          "%s is only %s px wide; a screenshot of the interface is wider" % (rel, w))

# --- the help box carries a docs link for every key -------------------------
# The "?" panel's "Read more in the Docs" is built from the help JSON the page
# embeds, so a key missing from it loses its link silently.
html = page.render("TOKEN-FOR-THE-TEST")
blob = re.search(r'<script id=helpdata type="application/json">(.*?)</script>',
                 html, re.S)
check(blob is not None, "the page has no help JSON for app.js to read")
if blob:
    import json
    help_json = json.loads(blob.group(1))
    for key, link in sorted(docs.HELP_LINKS.items()):
        got = (help_json.get(key) or {}).get("docs")
        check(got == link,
              "the help box for %s carries %r, not %r" % (key, got, link))

# --- the header button and the bar along the bottom -------------------------
check('id=docsbtn' in html and '/docs/index.html' in html,
      "the header has no Docs button pointing at /docs/index.html")
check('id=helpdocs' in html, "the help box has no 'Read more in the Docs' link")
# Last in the body, after the form and the footer. What follows it is only
# the stylesheet and the scripts, and app.js has to come after the bar or it
# cannot wire the bar's own link.
after = html.split('class="creditbar"')[-1]
check('class="creditbar"' in html and "</main>" in html.split('class="creditbar"')[0],
      "the credit bar does not come after the page's main content")
check(not re.search(r"<(?!/?(div|a|link|script)\b)[a-z]", after.split("</div>", 1)[-1]),
      "something other than the stylesheet and the scripts follows the credit bar")
check("app.js" in after, "app.js is loaded before the bar it has to wire")
check(docs.SITE.copyright.replace("&", "&amp;") in html,
      "the credit bar does not state the copyright")
check('id=citebar' in html or 'id="citebar"' in html,
      "the bar's How to cite has no id for app.js to wire to the citation dialog")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
check("$('citebar')" in appjs, "app.js no longer wires the bar's How to cite")
check("$('helpdocs')" in appjs, "app.js no longer fills in the help box's docs link")

# --- served, without a token ------------------------------------------------
app = server.App(os.path.join(HERE, "_work", "docs_runs"))
httpd = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
httpd.app = app
httpd.token = app.token
threading.Thread(target=httpd.serve_forever, daemon=True).start()
BASE = "http://127.0.0.1:%d" % httpd.server_address[1]
# A cluster that sets a web proxy must not have 127.0.0.1 sent to it.
os.environ["http_proxy"] = os.environ["HTTP_PROXY"] = "http://127.0.0.1:9"


def get(path, redirect=True):
    opener = urllib.request.build_opener(urllib.request.ProxyHandler({}))
    if not redirect:
        class NoRedirect(urllib.request.HTTPRedirectHandler):
            def redirect_request(self, *a, **k):
                return None
        opener = urllib.request.build_opener(urllib.request.ProxyHandler({}),
                                            NoRedirect)
    try:
        with opener.open(BASE + path, timeout=20) as r:
            return r.status, r.read(), r.headers
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read(), exc.headers


code, body, h = get("/docs/", redirect=False)
check(code == 302 and h.get("Location") == "/docs/index.html",
      "/docs/ does not redirect to the overview (%s %s)" % (code, h.get("Location")))
code, body, h = get("/docs")
check(code == 200 and b"<h1>" in body,
      "/docs does not reach the overview (%s)" % code)
code, body, h = get("/docs/index.html")
check(code == 200 and h.get("Content-Type", "").startswith("text/html"),
      "/docs/index.html is not served as HTML (%s %s)" % (code, h.get("Content-Type")))
check(b"GUI to ADCP Docking" in body and b'href="style.css' in body,
      "the page served is not the documentation's layout")
check(app.token.encode() not in body, "the docs carry the session token")
code, body, h = get("/docs/style.css")
check(code == 200 and h.get("Content-Type", "").startswith("text/css"),
      "/docs/style.css is not served as CSS (%s %s)" % (code, h.get("Content-Type")))
check(b"--accent:#1a56db" in body,
      "the stylesheet does not carry the interface's own colours")
one_img = sorted(used)[0]
code, body, h = get("/docs/" + one_img)
check(code == 200 and h.get("Content-Type") == "image/png",
      "/docs/%s is not served as a PNG (%s %s)" % (one_img, code, h.get("Content-Type")))
code, body, h = get("/docs/nope.html")
check(code == 404, "a page that does not exist is not a 404 (%s)" % code)
code, body, h = get("/docs/../server.py")
check(code == 404, "a path out of the docs folder is not refused (%s)" % code)
httpd.shutdown()
for k in ("http_proxy", "HTTP_PROXY"):
    os.environ.pop(k, None)

# --- the kit is the same file in every repository ---------------------------
# Four interfaces share docskit.py and are kept byte identical. The other
# copies are only present in a checkout that has them, so this is skipped
# rather than failed where they are absent.
SIBLING = "/scratch/group/p.agr250019.000/HOPE-pipelines/hope_monitor/docskit.py"
if os.path.isfile(SIBLING):
    mine = open(os.path.join(ROOT, "adcp_dock", "docskit.py"), "rb").read()
    check(mine == open(SIBLING, "rb").read(),
          "adcp_dock/docskit.py differs from %s; the four copies are kept "
          "identical" % SIBLING)
else:
    print("note  %s is not here, so the copies were not compared" % SIBLING)

if fails:
    print("test_docs FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the documentation audits clean, its pictures are real and it is "
      "served without a token, %d checks" % ran)
