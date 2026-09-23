#!/usr/bin/env python3
"""The mark: one SVG, every raster size, and the .ico.

Drawn in code rather than committed as an image, so the mark can be
changed by editing code.

The mark shows a peptide entering a pocket: a teal C-shaped receptor
groove with four residues on a diagonal running into its mouth. It must be
recognisable at 16 pixels in a taskbar.

    python3 brand/make_brand.py
"""
import base64
import os
import struct

from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))

DEEP = "#113f4a"
TEAL = "#1a8a8a"
WARM = "#d67a3a"
PAPER = "#f7f9f9"

LOGO_SVG = (
    '<svg width="62" height="62" viewBox="0 0 62 62" '
    'xmlns="http://www.w3.org/2000/svg">'
    '<defs><linearGradient id="ag" x1="0" y1="0" x2="1" y2="1">'
    '<stop offset="0%" stop-color="#1a8a8a"/>'
    '<stop offset="100%" stop-color="#113f4a"/></linearGradient></defs>'
    '<rect x="1" y="1" width="60" height="60" rx="13" fill="#113f4a"/>'
    '<path d="M46 14 A20 20 0 1 0 46 48" fill="none" '
    'stroke="url(#ag)" stroke-width="7" stroke-linecap="round"/>'
    '<path d="M22 40 L30 35 L38 30 L46 25" fill="none" stroke="#f7f9f9" '
    'stroke-width="2.6" stroke-linecap="round"/>'
    '<circle cx="22" cy="40" r="4.2" fill="#d67a3a"/>'
    '<circle cx="30" cy="35" r="4.2" fill="#f7f9f9"/>'
    '<circle cx="38" cy="30" r="4.2" fill="#d67a3a"/>'
    '<circle cx="46" cy="25" r="4.2" fill="#f7f9f9"/>'
    '</svg>'
)


def draw(size=512, plate=True):
    """The mark. `plate=False` gives the glyph alone, on transparency, for
    placing on a surface that already provides its own panel."""
    im = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    s = size / 62.0
    def S(v):
        return v * s
    if plate:
        d.rounded_rectangle([S(1), S(1), S(61), S(61)], radius=int(S(13)),
                            fill=DEEP)
    # the receptor groove, an open C
    d.arc([S(11), S(11), S(51), S(51)], start=315, end=225,
          fill=TEAL, width=max(1, int(S(7))))
    # the peptide running into it
    pts = [(22, 40), (30, 35), (38, 30), (46, 25)]
    d.line([(S(x), S(y)) for x, y in pts], fill=PAPER,
           width=max(1, int(S(2.6))), joint="curve")
    for i, (x, y) in enumerate(pts):
        r = S(4.2)
        d.ellipse([S(x) - r, S(y) - r, S(x) + r, S(y) + r],
                  fill=WARM if i % 2 == 0 else PAPER)
    return im


def write_ico(im, path, sizes=(16, 24, 32, 48, 64, 128, 256)):
    """A .ico Explorer can actually read at every size.

    PIL's own ICO writer compresses every size as PNG, and Windows falls
    back to the generic executable icon rather than read a PNG-compressed
    16x16. The small sizes go in as uncompressed BMP; tests/test_brand.py
    fails the build if that stops being true.
    """
    entries, blobs, offset = [], [], 6 + 16 * len(sizes)
    for n in sizes:
        r = im.resize((n, n), Image.LANCZOS)
        if n >= 64:
            import io
            b = io.BytesIO()
            r.save(b, format="PNG")
            data = b.getvalue()
        else:
            data = _bmp(r)
        entries.append((n, len(data), offset))
        blobs.append(data)
        offset += len(data)
    with open(path, "wb") as fh:
        fh.write(struct.pack("<HHH", 0, 1, len(sizes)))
        for n, length, off in entries:
            fh.write(struct.pack("<BBBBHHII", n % 256, n % 256, 0, 0,
                                 1, 32, length, off))
        for b in blobs:
            fh.write(b)


def _bmp(im):
    n = im.size[0]
    px = im.convert("RGBA").load()
    rows = []
    for y in range(n - 1, -1, -1):
        row = bytearray()
        for x in range(n):
            r, g, b, a = px[x, y]
            row += bytes((b, g, r, a))
        rows.append(bytes(row))
    xor = b"".join(rows)
    stride = ((n + 31) // 32) * 4
    header = struct.pack("<IiiHHIIiiII", 40, n, n * 2, 1, 32, 0,
                         len(xor) + stride * n, 0, 0, 0, 0)
    return header + xor + b"\x00" * (stride * n)


def main():
    im = draw(512)
    with open(os.path.join(HERE, "adcp_logo.svg"), "w") as fh:
        fh.write(LOGO_SVG)
    for n in (16, 24, 32, 48, 64, 128, 256, 512):
        im.resize((n, n), Image.LANCZOS).save(
            os.path.join(HERE, "adcp_logo_%d.png" % n))
    write_ico(im, os.path.join(HERE, "adcp.ico"))
    # The 128px mark, embedded in brand.py as base64 so a BUNDLED launcher
    # carries it without needing brand/ to travel alongside.
    with open(os.path.join(HERE, "adcp_logo_128.png"), "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode()
    print("wrote svg, 8 PNG sizes and adcp.ico in %s" % HERE)
    print("128px base64 is %d chars; paste into adcp_dock/brand.py "
          "if it changed" % len(b64))
    with open(os.path.join(HERE, "_logo_128.b64"), "w") as fh:
        for i in range(0, len(b64), 76):
            fh.write(b64[i:i + 76] + "\n")
    print("also written to brand/_logo_128.b64")


if __name__ == "__main__":
    main()
