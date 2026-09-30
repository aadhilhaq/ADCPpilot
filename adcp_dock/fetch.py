"""Fetching a structure from the RCSB, and caching it.

By Aadhil Haq

Cached under the runs directory rather than re-downloaded, because compute
nodes frequently have no outbound network: a run that succeeded on the
login node would otherwise fail on the node that executes it.
"""
import gzip
import os
import re
import urllib.error
import urllib.request

CODE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
#: PDB first, then mmCIF: some entries are too large for the PDB format and
#: only exist as CIF, and asking for the one that is not there is a 404 rather
#: than an error worth reporting.
SOURCES = [
    ("pdb", "https://files.rcsb.org/download/%s.pdb"),
    ("cif", "https://files.rcsb.org/download/%s.cif"),
]
TIMEOUT = 30


class FetchError(Exception):
    pass


def looks_like_code(text):
    return bool(CODE.match((text or "").strip()))


def cached(code, cachedir):
    code = code.strip().upper()
    for ext in ("pdb", "cif"):
        p = os.path.join(cachedir, "%s.%s" % (code, ext))
        if os.path.isfile(p) and os.path.getsize(p) > 200:
            return p
    return None


def fetch(code, cachedir, timeout=TIMEOUT, log=print):
    """Returns the path to the downloaded file. Cached hits do not touch the
    network at all, which is what makes this safe to call from a batch job."""
    code = (code or "").strip().upper()
    if not looks_like_code(code):
        raise FetchError("%r is not a PDB id. They are four characters, the "
                         "first a digit, such as 4CMM, 1FZV or 7LRI." % code)
    os.makedirs(cachedir, exist_ok=True)
    hit = cached(code, cachedir)
    if hit:
        log("  %s is already here: %s" % (code, hit))
        return hit
    problems = []
    for ext, url in SOURCES:
        out = os.path.join(cachedir, "%s.%s" % (code, ext))
        try:
            log("  fetching %s" % (url % code))
            with urllib.request.urlopen(url % code, timeout=timeout) as r:
                data = r.read()
        except urllib.error.HTTPError as exc:
            problems.append("%s: HTTP %s" % (ext, exc.code))
            continue
        except Exception as exc:                                # noqa: BLE001
            problems.append("%s: %s" % (ext, exc))
            continue
        if data[:2] == b"\x1f\x8b":
            data = gzip.decompress(data)
        if len(data) < 200:
            problems.append("%s: the file came back empty" % ext)
            continue
        tmp = out + ".part"
        with open(tmp, "wb") as fh:
            fh.write(data)
        os.replace(tmp, out)
        log("  %s: %d kB" % (os.path.basename(out), len(data) // 1024))
        return out
    raise FetchError("could not fetch %s from the RCSB (%s). If this machine "
                     "has no route to the internet, download the file "
                     "elsewhere and give the path instead."
                     % (code, "; ".join(problems)))
