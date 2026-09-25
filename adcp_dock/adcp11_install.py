"""Install ADCP 1.1 the way the Scripps Research Institute's own script does
(adcpsuite_micromamba.sh, October 2024), into one folder.

By Aadhil Haq

    python -m adcp_dock.adcp11_install [--dest DIR]      (or adcp-dock install-adcp11)

ADCP 1.1 is a separate Python 3.7 environment with its own OpenMM, Open Babel
and PDBFixer, plus the Scripps packages (adcp, adfr and their libraries). It
is created with micromamba, which is downloaded if it is not already here.
The Scripps packages are fetched from ccsb.scripps.edu and each is checked
against the SHA-256 below before it is installed.
"""
import argparse
import os
import shutil
import subprocess
import sys

from .adfr import _fetch, _sha256

BASE = "https://ccsb.scripps.edu/mamba/pip/py37"
DOWNLOADS = "https://ccsb.scripps.edu/adcpv11/downloads/"
#: On a computer without ADFRsuite set up yet, or without a checkout.
DEFAULT = os.path.expanduser("~/.adcp-dock/adcpsuite")
# What Scripps' script installs from conda-forge and from PyPI, with its pins.
CONDA = ["python=3.7", "setuptools=59.8.0", "numpy=1.21.6", "openbabel=2.4.1",
         "openmm=7.6.0", "parmed=3.4.3", "mmtf-python=1.1.3", "pdbfixer=1.8.1"]
PYPI = ["Pillow==9.5.0", "pyside2==5.15.2.1", "pybel==0.15.5",
        "biopython==1.81", "pypdb==2.3", "wget"]
#: path under BASE -> SHA-256 of the wheel Scripps published.
WHEELS = {
    "adcp/adcp-0.0.25": "5dd6ec8cf25824c495c03cda7f08df5c72237acfb40a050c9957682332c813d6",
    "adfr/adfr-0.0.22": "9f99868fa9900f9d1f6e44523856c1eabc8f4ad5ae53dc9bd6c52962a1f51e88",
    "adfrcc/adfrcc-0.0.1": "b4a5d56233e13a4a0290915f29ff5f46818ebe3fec56c31d7fad743b4f0d677b",
    "appframework/appframework-0.0.3": "7ee7fee27c338967e5981a7dd1b7757fc1433ca3c9407f31e2d9284d4c13285d",
    "autosite/autosite-0.0.4": "f88eeec78a7a9eba96d88b30084f9caa6f77378efeaf312a18c72562c9699d46",
    "bhtree/bhtree-0.0.1": "235d7bcc8fbbb1dbba73df52f62d0ac7882a11d79dc17b135b973c9f174327fb",
    "dejavu2/dejavu2-0.0.4": "3f9612c9015870515e627439333955fe34dfe402761cb63f0a517e845d553fbc",
    "geomutils/geomutils-0.0.1": "ca5f53509ccba026d1e390570d864c22c99e995bbf28bbdde8dc32e7cbfbdd5b",
    "gle/gle-0.0.1": "1411dcbd39ec61752b0681afe0c3715bae2791c9307c9d686abfbbe525f62207",
    "mglkey/mglkey-0.0.1": "1281932e4d58e9731ac3b36b5dc3786d492184cbbe6ba957f4f396b3d370a239",
    "mglutil/mglutil-0.0.4": "64f97431e48c6230043811bf8668eaa62f40e3bf7807a5bb1f5cba5185fca9c6",
    "molkit2/molkit2-0.0.5": "227180881ff8740440f2f8c7d532860cc651e2e29a01d4f92b856a68c0f23d91",
    "mslib/mslib-0.0.1": "43dcdd0e7619ffb0e49a3ba1978e0b575ba2416c4e87b03a970e010fcda669ec",
    "mslibcom/mslibcom-0.0.1": "7c3fbfbf11931fae088eef47a41fe8901f783a0d2143369b6a9d38e2ecf7cb22",
    "opengltk/opengltk-0.0.1": "8d36f73b5f0e0c660e463aa6c98cf840553eb72fe7354fcdaa14d8a44a8e3cbf",
    "pmvapp/pmvapp-0.0.6": "00a52227350335ea55ef53798efff612f87ec79847174cc9a544a05143c0aa14",
    "prody/prody-0.0.4": "b1326b7b8bd145e64185ec8adc083ff4210a4ef93c2bb36fb779bc573fc0aa49",
    "pyglf/pyglf-0.0.1": "9c89addc66ba42074e2d10a3db4039d90109679a4ce6da7fd2bfb09e04434635",
    "support/support-0.0.2": "a2bed345be1707874739b5f28204ab565504cd2571f540b8cdf56f99e0b4c9d5",
    "utpackages/utpackages-0.0.1": "e1276d5074548c45ad07af43f0e0460e83fa95ed61ce45d15d7837c300aa8b01",
    "volume/volume-0.0.2": "87616a7e4e6e277fe73047f5e7986c7c589a004ad284a7a5f04370a617a5980c",
    "legacy/autodocktools/autodocktools-0.0.1": "0afddc7bba1590677e99cd863c353e00960bea8d4a03e6b523485b440f95901f",
    "legacy/molkit/molkit-0.0.2": "47fae4f93986e51a85455c7921cdaa7e9d7d39c61e425c9922d2b14932f24295",
    "legacy/pybabel/pybabel-0.0.1": "993dea44d34306611d4b028826d2e1ed8be5fd04ff10bcc028966346ae11664d",
}
# Installed first and one at a time, then the rest without dependencies, as
# in Scripps' script.
FIRST = ("prody", "autosite", "mslib")


def _run(cmd, log, **kw):
    log("  $ " + " ".join(cmd[:6]) + (" ..." if len(cmd) > 6 else ""))
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if p.returncode:
        raise SystemExit("this step failed:\n  %s\n%s"
                         % (" ".join(cmd), (p.stdout + p.stderr)[-3000:]))
    return p


def micromamba(where, log=print):
    """A micromamba executable: one on PATH, or one downloaded into `where`."""
    found = shutil.which("micromamba")
    if found and os.path.isfile(found):
        return found
    exe = os.path.join(where, "bin", "micromamba")
    if not os.path.isfile(exe):
        os.makedirs(where, exist_ok=True)
        log("  downloading micromamba")
        arch = "linux-64" if sys.platform.startswith("linux") else (
            "osx-arm64" if os.uname().machine == "arm64" else "osx-64")
        _run(["bash", "-c", "curl -Ls https://micro.mamba.pm/api/micromamba/%s/"
              "latest | tar -xj -C %s bin/micromamba" % (arch, where)], log)
    return exe


def default_dest():
    """Where ADCP 1.1 goes when no --dest is given: beside ADFRsuite, in the
    installation's own tools folder, where the interface looks for it first.
    Not the home directory, whose file allowance on a cluster an environment
    of this size can use up."""
    adfr = os.environ.get("ADFRSUITE_HOME")
    if adfr:
        return os.path.join(os.path.dirname(os.path.abspath(adfr)), "adcpsuite")
    home = os.environ.get("ADCP_HOME")
    if home:
        return os.path.join(home, "tools", "docking", "adcpsuite")
    return DEFAULT


def install(dest=None, log=print):
    """Create the ADCP 1.1 environment at `dest`; returns its adcp executable."""
    from . import adcp11
    dest = os.path.abspath(os.path.expanduser(dest or default_dest()))
    if adcp11.find([os.path.dirname(dest)]) or os.path.isfile(
            os.path.join(dest, "bin", "adcp")):
        log("ADCP 1.1 is already installed at %s" % dest)
        return os.path.join(dest, "bin", "adcp")
    base = os.path.dirname(dest)
    mm = micromamba(os.path.join(base, "micromamba"), log)
    # Downloads and caches stay beside the installation, not in home.
    env = dict(os.environ, MAMBA_ROOT_PREFIX=os.path.join(base, "mamba_root"),
               XDG_CACHE_HOME=os.path.join(base, "cache"),
               PIP_CACHE_DIR=os.path.join(base, "cache", "pip"))
    log("creating the ADCP 1.1 environment at %s (Python 3.7, conda-forge)" % dest)
    _run([mm, "create", "-y", "-q", "-p", dest, "-c", "conda-forge", "--no-rc"]
         + CONDA, log, env=env)
    py = os.path.join(dest, "bin", "python")
    _run([py, "-m", "pip", "install", "-q"] + PYPI, log, env=env)
    cache = os.path.join(base, "adcpsuite_pkgs")
    os.makedirs(cache, exist_ok=True)
    log("downloading the ADCP 1.1 packages from the Scripps Research Institute")
    for path, sha in WHEELS.items():
        name = os.path.basename(path) + "-py3-none-any.whl"
        target = os.path.join(cache, name)
        if not (os.path.isfile(target) and _sha256(target) == sha):
            _fetch("%s/%s-py3-none-any.whl" % (BASE, path), target, log)
            got = _sha256(target)
            if got != sha:
                os.remove(target)
                raise SystemExit("%s does not match its pinned checksum and "
                                 "was deleted (got %s). Retry, or see %s."
                                 % (name, got, DOWNLOADS))
    wheels = {os.path.basename(p).split("-")[0]: os.path.join(
        cache, os.path.basename(p) + "-py3-none-any.whl") for p in WHEELS}
    for name in FIRST:
        _run([py, "-m", "pip", "install", "-q", "--no-index",
              "--find-links=" + cache, wheels[name]], log, env=env)
    rest = [w for n, w in wheels.items() if n not in FIRST]
    _run([py, "-m", "pip", "install", "-q", "--no-cache-dir", "--no-dependencies",
          "--no-index", "--find-links=" + cache] + rest, log, env=env)
    exe = adcp11.find([base]) or os.path.join(dest, "bin", "adcp")
    log("installed ADCP %s at %s" % (adcp11.version(exe) or "1.1", dest))
    return exe


def main(argv=None):
    ap = argparse.ArgumentParser(prog="adcp-dock install-adcp11",
                                 description=__doc__.split("\n\n")[0])
    ap.add_argument("--dest", default=None,
                    help="where to create it (default: beside ADFRsuite, "
                         "or %s)" % DEFAULT)
    a = ap.parse_args(argv)
    install(a.dest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
