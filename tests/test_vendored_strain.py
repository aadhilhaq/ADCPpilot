"""The peptide strain in adcp_pilot/_vendor is hope_core's, byte for byte.

By Aadhil Haq

adcp_pilot/_vendor/strain.py is the module HOPE-pipelines maintains as
hope_core/strain.py, carried here unchanged, which is the only reason the two
tools' strain columns can be compared. Unlike the vendored preparation it
keeps the master's own docstrings, so the two are compared as bytes: any
difference fails, naming the first line where they part.

hope_core/strain.py is looked for in a HOPE-pipelines checkout beside this one
or beside any folder above it, and HOPE_PIPELINES_CHECKOUT names another,
which is then the only place looked, as tests/test_vendored_prepare.py does.
Where none holds it the test is a skip naming the places it looked, never a
pass.

    python3 tests/test_vendored_strain.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MINE = os.path.join(ROOT, "adcp_pilot", "_vendor", "strain.py")


def checkouts():
    """The HOPE-pipelines checkouts to look in, in order."""
    named = os.environ.get("HOPE_PIPELINES_CHECKOUT", "").strip()
    if named:
        return [named]
    out, folder = [], os.path.dirname(ROOT)
    while True:
        out.append(os.path.join(folder, "HOPE-pipelines"))
        if os.path.dirname(folder) == folder:
            return out
        folder = os.path.dirname(folder)


places = checkouts()
master = next((m for m in (os.path.join(c, "hope_core", "strain.py") for c in places)
               if os.path.isfile(m)), None)
if not master:
    print("skip  adcp_pilot/_vendor/strain.py was compared with nothing: no "
          "HOPE-pipelines checkout holding hope_core/strain.py is %s (looked at %s)"
          % ("where HOPE_PIPELINES_CHECKOUT says" if os.environ.get(
              "HOPE_PIPELINES_CHECKOUT", "").strip()
             else "beside this one or beside a folder above it", ", ".join(places)))
    sys.exit(77)

mine, theirs = open(MINE, "rb").read(), open(master, "rb").read()
if mine != theirs:
    a, b = mine.splitlines(), theirs.splitlines()
    at = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
    print("test_vendored_strain FAILED")
    print("  FAIL adcp_pilot/_vendor/strain.py is not %s: they part at line %d "
          "(%d lines here, %d there). Copy it again rather than edit it here"
          % (master, at + 1, len(a), len(b)))
    sys.exit(1)
print("ok  adcp_pilot/_vendor/strain.py is %s byte for byte, %d bytes"
      % (master, len(mine)))
