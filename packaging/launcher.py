"""The script PyInstaller bundles.

By Aadhil Haq

A bundled executable runs its entry script as a top-level script with no
package around it, so this imports ABSOLUTELY and does nothing else.
Pointing PyInstaller at adcp_launcher/__main__.py instead produces a build
that fails on the first double-click with "attempted relative import with no
known parent package", which is the one failure a windowed executable
cannot report, because a windowed executable has no stdout.
"""

import os
import sys


def _stamp():
    """Publish the commit this bundle was built from.

    A bundle has no git to ask, so the workflow writes BUILD_VERSION into it
    and this puts it where brand.version() looks. Without it the window says
    "build unknown", which is unhelpful when a user wants to know whether
    they have the latest version.

    sys._MEIPASS is where PyInstaller unpacks a one-file bundle; a folder
    build has the file beside the executable instead, so both are tried.
    """
    base = getattr(sys, "_MEIPASS", os.path.dirname(os.path.abspath(__file__)))
    for candidate in (os.path.join(base, "BUILD_VERSION"),
                      os.path.join(os.path.dirname(base), "BUILD_VERSION"),
                      os.path.join(os.path.dirname(os.path.dirname(base)),
                                   "BUILD_VERSION")):
        try:
            with open(candidate) as fh:
                text = fh.read().strip()
            if text:
                os.environ.setdefault("ADCP_BUILD", text)
                return
        except OSError:
            continue


_stamp()

from adcp_launcher.entry import main                            # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
