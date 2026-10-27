"""The icon is one Windows will draw.

By Aadhil Haq

PIL's ICO writer compresses every size as PNG, and Explorer falls back to
the generic executable icon rather than read a PNG-compressed 16x16 image.
The small sizes are therefore written as uncompressed BMP, and this fails
the build if that stops being true.

    python3 tests/test_brand.py
"""
import os
import struct
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BRAND = os.path.join(ROOT, "brand")

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


# The mark is drawn by a script, so the build can always make one. If it has
# not been run, run it: an icon missing from a checkout is not a failure, it
# is a file no one generated yet.
ICO = os.path.join(BRAND, "adcp.ico")
PNG = os.path.join(BRAND, "adcp_logo_512.png")
if not (os.path.isfile(ICO) and os.path.isfile(PNG)):
    try:
        subprocess.run([sys.executable, os.path.join(BRAND, "make_brand.py")],
                       check=True, capture_output=True)
    except Exception as exc:                                    # noqa: BLE001
        print("the icon could not be drawn (%s); Pillow is needed only for "
              "this. Skipped." % exc)
        raise SystemExit(77)

check(os.path.isfile(ICO), "there is no adcp.ico")
check(os.path.isfile(PNG), "there is no 512px PNG for the macOS build")

with open(ICO, "rb") as fh:
    blob = fh.read()
reserved, kind, count = struct.unpack("<HHH", blob[:6])
check(reserved == 0 and kind == 1,
      "adcp.ico does not start with an icon header")
check(count >= 5, "the icon carries %d size(s); Windows picks from several "
                  "and falls back when the one it wants is absent" % count)

sizes, png_small = [], []
for i in range(count):
    w, h, _c, _r, _p, _bpp, length, offset = struct.unpack(
        "<BBBBHHII", blob[6 + 16 * i:6 + 16 * (i + 1)])
    n = w or 256
    sizes.append(n)
    data = blob[offset:offset + length]
    check(len(data) == length,
          "the %dx%d entry claims %d bytes and the file holds %d"
          % (n, n, length, len(data)))
    is_png = data[:8] == b"\x89PNG\r\n\x1a\n"
    if n <= 48 and is_png:
        png_small.append(n)

check(16 in sizes, "there is no 16x16; that is the one the taskbar and the "
                   "title bar use")
check(32 in sizes, "there is no 32x32")
check(not png_small,
      "the %s entr%s PNG-compressed. Explorer will not read a "
      "PNG-compressed icon at that size and draws the generic executable "
      "icon instead: the build succeeds and the thing on the desktop looks "
      "anonymous." % (", ".join("%dx%d" % (n, n) for n in png_small),
                      "y is" if len(png_small) == 1 else "ies are"))

# and the PNG the macOS build turns into an .icns is square and big enough
try:
    from PIL import Image
    im = Image.open(PNG)
    check(im.size == (512, 512),
          "the logo is %s; iconutil wants 512x512 to make every size from"
          % (im.size,))
    check(im.mode in ("RGBA", "RGB"), "the logo is in mode %r" % im.mode)
except ImportError:
    pass

# one place says who made this, and the icon script is part of the build
check(os.path.isfile(os.path.join(BRAND, "make_brand.py")),
      "the icon is committed with no script that makes it, so it cannot be "
      "changed by anyone who does not have the original")

# --- the tab wears the same mark as the launcher --------------------------
#
# The served page had no <link rel=icon> at all, so a browser with several
# tabs open showed this one blank among them, while the launcher that
# started it carries a logo. Inlined as data URIs: no route to add, no second
# file to ship, and no cache to go stale against the launcher's copy.
sys.path.insert(0, ROOT)
from adcp_dock import page as _page, brand as _brand              # noqa: E402
import base64 as _b64                                            # noqa: E402
import re as _re                                                 # noqa: E402

_html = _page.render("token")
check("__FAVICON__" not in _html and "__FAVICON_PNG__" not in _html,
      "a favicon placeholder was left in the served page")
_m = _re.search(r'<link rel="icon" type="image/svg\+xml" href="data:image/'
                r'svg\+xml;base64,([^"]+)"', _html)
check(_m, "the page has no SVG favicon, so the tab shows a blank page icon")
if _m:
    _svg = _b64.b64decode(_m.group(1)).decode()
    check(_svg == _brand.LOGO_SVG,
          "the tab icon is not the launcher's mark, so the two disagree "
          "about what this program looks like")
_p = _re.search(r'type="image/png"[^>]*href="data:image/png;base64,([^"]+)"',
                _html)
check(_p, "no PNG fallback, so Safari before 16 shows nothing")
if _p:
    _raw = _b64.b64decode(_p.group(1))
    check(_raw[:8] == b"\x89PNG\r\n\x1a\n", "the PNG fallback is not a PNG")
    check(_raw == _brand.logo_png(),
          "the PNG fallback is not the launcher's mark")


# --- the help marker matches the related tools -------------------------------
#
# A 15px circle with a 10px "?", the same as the related pipelines use, so it
# reads as something to press rather than as part of the label.
_css = open(os.path.join(ROOT, "adcp_dock", "page.py")).read()
_q = _css[_css.index(".q{"):_css.index(".q:hover")]
# No box at all. A <button> arrives with a user-agent border, background,
# padding and minimum size, each of which overrides a declared size.
# all:unset removes them, leaving a "?" in the accent color at the size of
# the text beside it.
check("all:unset" in _q,
      "the marker does not reset the user-agent button styling, so its "
      "border, background and minimum size show through whatever else is "
      "set: %r" % _q[:200])
for _gone in ("border-radius", "background:", "border:"):
    check(_gone not in _q,
          "the marker still draws a box (%s), which is the thing that was "
          "wrong with it" % _gone)
# a SET width, not min-width:0, which is the opposite of a box
import re as _re5
check(not _re5.search(r'(?:^|;)\s*(?:width|height)\s*:', _q),
      "the marker sets an explicit width or height: %r" % _q[:180])
check("font-size:.92em" in _q,
      "the marker is not sized against the text it sits beside")
check("color:var(--accent)" in _q,
      "nothing marks the ? as something you can press")


if fails:
    print("test_brand FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the icon is one Explorer can read at every size, %d checks" % ran)
