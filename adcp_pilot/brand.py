"""Name, tagline and palette, defined once.

By Aadhil Haq

The launcher, the page, the CLI and the README read these rather than
carrying their own copies.

The mark travels as bytes: a bundled launcher has no checkout beside it, so
the 128-pixel PNG is embedded here in base64. Everything else in brand/
(the .ico for the Windows build, the .icns for macOS, the other sizes) is
for the build rather than the running program.
"""
from __future__ import annotations

import base64
import os
import subprocess

# The name states what this is: ADCP is the Scripps Research Institute's, and
# a pilot is what this package provides, a way to fly a docking study from one
# page. Everything that displays a name reads these three.
NAME = "ADCPpilot"
TAGLINE = "ADCP by the Scripps Research Institute; the interface by HOPE Lab"
LONG_NAME = "A graphical interface to AutoDock CrankPep, end to end on a cluster"

LAB = "HOPE Lab"
AUTHORS = ("Aadhil Haq",)
PI = "Dr. Sandun Fernando"

#: The palette the mark and the web page share.
DEEP = "#113f4a"
TEAL = "#1a8a8a"
WARM = "#d67a3a"
PAPER = "#f7f9f9"

LOGO_SVG = """<svg width="62" height="62" viewBox="0 0 62 62" xmlns="http://www.w3.org/2000/svg"><defs><linearGradient id="ag" x1="0" y1="0" x2="1" y2="1"><stop offset="0%" stop-color="#1a8a8a"/><stop offset="100%" stop-color="#113f4a"/></linearGradient></defs><rect x="1" y="1" width="60" height="60" rx="13" fill="#113f4a"/><path d="M46 14 A20 20 0 1 0 46 48" fill="none" stroke="url(#ag)" stroke-width="7" stroke-linecap="round"/><path d="M22 40 L30 35 L38 30 L46 25" fill="none" stroke="#f7f9f9" stroke-width="2.6" stroke-linecap="round"/><circle cx="22" cy="40" r="4.2" fill="#d67a3a"/><circle cx="30" cy="35" r="4.2" fill="#f7f9f9"/><circle cx="38" cy="30" r="4.2" fill="#d67a3a"/><circle cx="46" cy="25" r="4.2" fill="#f7f9f9"/></svg>"""

_LOGO_PNG_B64 = (
    "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"
    "SUVORK5CYII="
)


def logo_png():
    """The 128-pixel mark, as PNG bytes."""
    return base64.b64decode("".join(_LOGO_PNG_B64.split()))


def asset(name):
    """A path inside brand/, for the files that are not embedded."""
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(here, "brand", name)


def version():
    """The current commit, short.

    A bundled launcher has no git and no checkout, so the build writes the
    commit into the bundle and the entry script puts it here. Falling back to
    git covers running from source; "unknown" covers a tarball.
    """
    stamped = os.environ.get("ADCP_BUILD")
    if stamped:
        return stamped.split()[0]
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for candidate in (here, os.path.join(here, "..")):
        try:
            out = subprocess.run(
                ["git", "-C", candidate, "rev-parse", "--short", "HEAD"],
                capture_output=True, text=True, timeout=10)
            v = (out.stdout or "").strip()
            if v:
                return v
        except Exception:                                      # noqa: BLE001
            pass
    return "unknown"


def credit(sep=" \u00b7 "):
    from . import __version__
    return sep.join(["%s %s" % (NAME, __version__),
                     "a project by %s" % LAB,
                     ", ".join(AUTHORS),
                     "PI %s" % PI])
