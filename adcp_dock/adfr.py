"""Download, verify and install ADFRsuite, the Scripps Research Institute's
docking suite, which provides adcp, agfr and prepare_receptor.

By Aadhil Haq

    python -m adcp_dock.adfr --dest ~/tools/ADFRsuite

Standard library only: this runs before anything else is installed, on a
fresh machine.

DOWNLOAD VERIFICATION. ccsb.scripps.edu serves its certificate without the
intermediate that chains it to a trusted root, so scripted downloads fail
certificate verification even though browsers, which fetch the missing
intermediate themselves, succeed. The download is therefore verified by
content: each tarball's SHA-256 is pinned below. A verified connection is
tried first; if the chain fails, the file is fetched without it and must
still match the pinned hash exactly, or it is deleted and the install stops.

The pins were established independently of the download that produced them:
the Linux tarball's adcp binary and Python sources are byte-identical to an
ADFRsuite installed from Scripps six weeks earlier, and the macOS tarball's
platform-independent sources are byte-identical to the Linux tarball's.

LICENSING. ADFRsuite is distributed mostly under the LGPL 2, but bundles
components that are free for academic use only: MSMS (Scripps), the UT
packages (UT Austin) and STRIDE. Every installation unpacks them, but the
programs this package runs do not use them. The MSMS library is activated
only when the installing user answers an academic-or-commercial question,
and this installer does not answer it on their behalf: it passes -l, the
installer's own switch for skipping that question, so it stays inactive.
Scripps asks that ADFRsuite be downloaded from their site rather than
redistributed, which is what this does.
"""
import argparse
import hashlib
import os
import platform
import shutil
import ssl
import subprocess
import sys
import tarfile
import urllib.request

#: platform -> (download URL, file name, SHA-256, bytes)
BUILDS = {
    "linux-x86_64": (
        "https://ccsb.scripps.edu/adfr/download/1038/",
        "ADFRsuite_x86_64Linux_1.0.tar.gz",
        "66539864ffb7ba87728d52ee50dde78586f551c7ca96aa455dae761875faf7e4",
        103496717),
    "darwin-x86_64": (
        "https://ccsb.scripps.edu/adfr/download/1033/",
        "ADFRsuite_x86_64Darwin_1.0.tar.gz",
        "9556680809c76859b159b7709fa1a6e84716e2cb41a0ab8bd72aecbb7960b413",
        90475481),
}
BINARIES = ("adcp", "agfr", "prepare_receptor")
PAGE = "https://ccsb.scripps.edu/adfr/downloads/"
LICENCE_NOTE = (
    "ADFRsuite is the Scripps Research Institute's. It is mostly LGPL 2, but "
    "bundles MSMS, the UT packages and STRIDE, which are free for academic "
    "use only; commercial use needs licences from their owners. See "
    "LICENSE.txt in the installed suite.")


def which_build():
    """The build for this machine, and a note if it will run translated."""
    sysname = platform.system().lower()
    mach = platform.machine().lower()
    if sysname == "linux" and mach in ("x86_64", "amd64"):
        return "linux-x86_64", ""
    if sysname == "darwin":
        note = ""
        if mach in ("arm64", "aarch64"):
            note = ("This Mac has Apple silicon. ADFRsuite is built for Intel "
                    "and runs through Rosetta 2; if it is not installed, run "
                    "`softwareupdate --install-rosetta` once.")
        return "darwin-x86_64", note
    if sysname == "linux":
        raise SystemExit("ADFRsuite is built for x86-64 Linux; this machine "
                         "is %s. There is no build for it." % mach)
    raise SystemExit(
        "No ADFRsuite build for %s here. On Windows, use WSL (Ubuntu) and run "
        "this inside it; WSL runs the Linux build." % platform.system())


def installed(dest):
    return all(os.path.isfile(os.path.join(dest, "bin", b)) for b in BINARIES)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _fetch(url, path, log):
    """Verified TLS first; the unverified fallback exists only because the
    caller checks the SHA-256 afterwards and discards anything that fails."""
    def get(ctx):
        req = urllib.request.Request(url, headers={"User-Agent": "adcp-dock"})
        with urllib.request.urlopen(req, context=ctx, timeout=120) as r, \
                open(path + ".part", "wb") as out:
            shutil.copyfileobj(r, out, 1 << 20)
        os.replace(path + ".part", path)
    try:
        get(ssl.create_default_context())
        return "verified connection"
    except (ssl.SSLError, urllib.error.URLError) as exc:
        if "CERTIFICATE" not in str(exc).upper():
            raise
        log("  the server's certificate chain is incomplete (a known "
            "misconfiguration at ccsb.scripps.edu); fetching without it and "
            "checking the file's SHA-256 instead")
        get(ssl._create_unverified_context())                  # noqa: SLF001
        return "unverified connection, SHA-256 checked"


def download(key, cache, log=print):
    url, name, sha, size = BUILDS[key]
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, name)
    if os.path.isfile(path) and _sha256(path) == sha:
        log("  %s is already here and matches its checksum" % name)
        return path
    log("  downloading %s from the Scripps Research Institute (%.0f MB)"
        % (name, size / 1e6))
    try:
        how = _fetch(url, path, log)
    except Exception as exc:                                    # noqa: BLE001
        raise SystemExit(
            "Could not download ADFRsuite (%s).\n\nGet %s by hand from\n  %s\n"
            "put it in %s, and run this again; it will use the file it finds."
            % (exc, name, PAGE, cache))
    got = _sha256(path)
    if got != sha:
        os.remove(path)
        raise SystemExit(
            "The downloaded %s does not match its pinned checksum and has been "
            "deleted.\n  expected %s\n  got      %s\nThis means the file was "
            "damaged in transit or is not the file Scripps published. Retry; "
            "if it persists, download it by hand from %s." % (name, sha, got,
                                                              PAGE))
    log("  checksum matches (%s)" % how)
    return path


def install(dest, cache=None, log=print):
    """ADFRsuite into `dest`. Returns `dest`. Safe to call again."""
    dest = os.path.abspath(os.path.expanduser(dest))
    if installed(dest):
        log("ADFRsuite is already installed in %s" % dest)
        return dest
    key, note = which_build()
    if note:
        log(note)
    cache = cache or os.path.dirname(dest)
    tar = download(key, cache, log)

    log("  unpacking")
    with tarfile.open(tar) as t:
        top = t.getnames()[0].split("/")[0]
        # Restrict extraction to paths inside the archive's own
        # folder: tarfile otherwise writes wherever a member name
        # says, including "../" and absolute paths.
        for m in t.getmembers():
            p = os.path.normpath(m.name)
            if p.startswith("..") or os.path.isabs(p):
                raise SystemExit("refusing %s: a member escapes the archive"
                                 % m.name)
        # Use the 'tar' filter where Python provides it (3.12+): it
        # refuses the same escapes and preserves the executable bits
        # the installer needs. Older versions rely on the check
        # above.
        if hasattr(tarfile, "tar_filter"):
            t.extractall(cache, filter="tar")
        else:
            t.extractall(cache)
    src = os.path.join(cache, top)

    # The academic-or-commercial question about MSMS is skipped with -l
    # rather than answered: it is a licence declaration only the installing
    # user can make, and nothing in this package calls MSMS.
    # stdin is /dev/null so no prompt can block the install.
    log("  running the Scripps installer (a few minutes)")
    logf = os.path.join(cache, "adfr_install.log")
    with open(logf, "w") as fh:
        subprocess.run(["bash", "./install.sh", "-d", dest, "-c", "0", "-l"],
                       cwd=src, stdin=subprocess.DEVNULL, stdout=fh,
                       stderr=subprocess.STDOUT)
    missing = [b for b in BINARIES
               if not os.path.isfile(os.path.join(dest, "bin", b))]
    if missing:
        raise SystemExit("ADFRsuite did not install: missing %s. The "
                         "installer's own log is %s; its last lines usually "
                         "say why. Remove %s before retrying: a half-installed "
                         "suite is worse than none, because the next run "
                         "would find it and trust it."
                         % (", ".join(missing), logf, dest))
    shutil.rmtree(src, ignore_errors=True)
    log("ADFRsuite installed in %s" % dest)
    log(LICENCE_NOTE)
    log("The MSMS library was not activated: that needs you to state academic "
        "or commercial use, and this tool does not use it. To activate it, run "
        "ADFRsuite's own installer yourself, without -l.")
    return dest


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m adcp_dock.adfr",
                                 description=__doc__.split("\n\n")[0])
    ap.add_argument("--dest", default="~/tools/ADFRsuite",
                    help="where to install it (default ~/tools/ADFRsuite)")
    ap.add_argument("--cache", default=None,
                    help="where the tarball goes (default: beside --dest)")
    ap.add_argument("--download-only", action="store_true")
    a = ap.parse_args(argv)
    if a.download_only:
        key, _n = which_build()
        download(key, a.cache or os.path.dirname(
            os.path.abspath(os.path.expanduser(a.dest))))
        return 0
    install(a.dest, a.cache)
    print("\nexport ADFRSUITE_HOME=%s"
          % os.path.abspath(os.path.expanduser(a.dest)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
