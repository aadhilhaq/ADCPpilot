"""Draw the README banner from the same palette as the application mark.

By Aadhil Haq

Committed as a script and not only as a PNG so the banner can be redrawn
when the wording changes. `python brand/make_banner.py` rewrites
brand/banner.png in place.
"""
import os
import sys

from PIL import Image, ImageDraw, ImageFont

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from make_brand import DEEP, TEAL, WARM, PAPER, draw            # noqa: E402
# The name and tagline come from the package, not a second copy here.
sys.path.insert(0, os.path.dirname(HERE))
from adcp_pilot.brand import NAME as BRAND_NAME, TAGLINE as BRAND_TAGLINE  # noqa: E402,E501

MUTED = "#9fc4c4"

# The README banner.
STANDARD = dict(W=1280, H=340, title=58, body=21, small=19, x=248,
                y_title=96, y_tag=176, y_rule=252, y_small=274,
                rule=(96, 4), plate=156, plate_x=64, glyph=120, radius=34,
                amp=34, dot=9, line=5, motif_min=940, motif_pad=48,
                motif_gap=72, name="banner.png")


def _font(name, size):
    """DejaVu, from matplotlib's copy or the system, whichever is present."""
    here = []
    try:
        import matplotlib
        here.append(os.path.join(os.path.dirname(matplotlib.__file__),
                                 "mpl-data", "fonts", "ttf"))
    except Exception:                                           # noqa: BLE001
        pass
    here += ["/usr/share/fonts/dejavu", "/usr/share/fonts/truetype/dejavu",
             "/Library/Fonts", "C:\\Windows\\Fonts"]
    for d in here:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return ImageFont.truetype(p, size)
    # A banner is a nicety. Falling back to the bitmap font gives something
    # ugly rather than nothing, and never fails a build over a missing font.
    return ImageFont.load_default()


def _mix(a, b, t):
    a = tuple(int(a[i:i + 2], 16) for i in (1, 3, 5))
    b = tuple(int(b[i:i + 2], 16) for i in (1, 3, 5))
    return tuple(int(a[i] + (b[i] - a[i]) * t) for i in range(3))


def banner(path=None, L=STANDARD):
    W, H = L["W"], L["H"]
    im = Image.new("RGB", (W, H), DEEP)
    d = ImageDraw.Draw(im)

    # A diagonal wash rather than a flat fill, so the mark on the left sits
    # on the darkest part and the text keeps its contrast across the width.
    for x in range(W):
        d.line([(x, 0), (x, H)], fill=_mix(DEEP, "#0b2e37", x / float(W)))

    f_title = _font("DejaVuSans-Bold.ttf", L["title"])
    f_body = _font("DejaVuSans.ttf", L["body"])
    f_small = _font("DejaVuSans.ttf", L["small"])

    x = L["x"]
    lines = [
        (L["y_title"], BRAND_NAME, f_title, PAPER),
        # One line, not two. Broken after "Institute", the first half would
        # sit under the title and read as though the GUI were Scripps'.
        (L["y_tag"], BRAND_TAGLINE, f_body, MUTED),
        (L["y_small"], "Peptide docking, from receptor to ranked table.",
         f_small, MUTED),
    ]

    # The motif goes beside the bottom line, which is the short one. Its
    # start is measured from the text, so it never runs through a sentence
    # when the wording changes.
    last_y, last_t, last_f, _ = lines[-1]
    free = x + int(d.textlength(last_t, font=last_f)) + L["motif_gap"]
    start = max(free, L["motif_min"])
    if W - start > 5 * L["motif_pad"]:
        span = W - start - L["motif_pad"]
        step = span / 6.0
        base = last_y + L["dot"] * 0.67
        # Where the zigzag runs under the tagline, its peaks stop short of
        # the tagline's descenders.
        amp = L["amp"]
        tag_end = x + d.textlength(BRAND_TAGLINE, font=f_body)
        if start < tag_end:
            tag_bottom = L["y_tag"] + L["body"] * 1.25
            amp = min(amp, base - L["dot"] - tag_bottom - 6)
        pts = [(start + i * step, base - (amp if i % 2 else 0))
               for i in range(7)]
        d.line(pts, fill=_mix(DEEP, TEAL, 0.42), width=L["line"],
               joint="curve")
        for i, (px, py) in enumerate(pts):
            r = L["dot"]
            d.ellipse([px - r, py - r, px + r, py + r],
                      fill=_mix(DEEP, WARM if i % 2 == 0 else TEAL, 0.55))

    # The mark, on a panel a shade lighter than the ground. The app icon
    # carries its own background; here the panel is the only one.
    plate, px0, g = L["plate"], L["plate_x"], L["glyph"]
    py = (H - plate) // 2
    d.rounded_rectangle([px0, py, px0 + plate, py + plate],
                        radius=L["radius"], fill=_mix(DEEP, TEAL, 0.18))
    mark = draw(512, plate=False).resize((g, g), Image.LANCZOS)
    im.paste(mark, (px0 + (plate - g) // 2, py + (plate - g) // 2), mark)

    for ly, text, font, colour in lines:
        d.text((x, ly), text, font=font, fill=colour)

    # A short warm rule under the wordmark: an accent, not a divider.
    rw, rh = L["rule"]
    d.rectangle([x, L["y_rule"], x + rw, L["y_rule"] + rh], fill=WARM)

    path = path or os.path.join(HERE, L["name"])
    im.save(path, optimize=True)
    return path


if __name__ == "__main__":
    print(banner())
