"""ADCP 1.0 on a Mac: the one place it still runs.

By Aadhil Haq

The Scripps Research Institute builds the ADCP 1.1 suite for Linux only, so
`adcp-dock local` on macOS (the desktop launcher's "This computer") docks
with ADFRsuite 1.0 instead: its adcp, agfr (AutoSite 1.0), prepare_receptor,
reduce and msms. Everywhere else, Linux, Colab and every cluster, runs ADCP
1.1, and nothing here is reached: every entry point asks engaged(), which is
False on any system but macOS.

Straight after each docking the peptide's folder is put into the layout the
rest of the package reads, the ADCP 1.1 one: modes.json written from the
mode table in adcp.log, the ranked poses ADCP 1.0 writes itself
(out_ranked_<n>.pdb), and its loose replica files packed into one archive or
removed. The reader of ADCP 1.0's mode table lives in this module and is
used nowhere else.

A run docked this way says so: its job file and every result name the
engine "ADCP 1.0 (macOS)", and the run's page and report say that its scores
are not comparable with runs docked with ADCP 1.1.

    python -m adcp_dock.mac10 [--dest ~/.adcp-dock/ADFRsuite]   (macOS only)

DOWNLOAD VERIFICATION. ccsb.scripps.edu has served its certificate without
the intermediate that chains it to a trusted root, so the tarball is
verified by content: its SHA-256 is pinned below, and a download that does
not match it is deleted and the install stops.

LICENSING. ADFRsuite is mostly LGPL 2, but bundles MSMS and the UT packages,
free for academic use only, and six packages in CCSBpckgs (MolKit, mglutil,
Support, DejaVu2, PyBabel and Volume) carry the MGLTOOLS license, for
non-commercial use only; prepare_receptor and agfr load five of the six.
The installer passes -l, which skips Scripps' academic-or-commercial
question about the MSMS library, an answer only the installing user can
give; the msms program is in bin/msms either way, but the macOS build is
32-bit (ppc and i386), which macOS 10.15 and later do not run, so there the
viewers draw their own surfaces. THIRD_PARTY_NOTICES.md sets this out.
"""
import argparse
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tarfile

#: How a run docked here names its engine, in job.json (adcp_engine) and in
#: every result (engine).
ENGINE = "1.0 (macOS)"
LABEL = "ADCP 1.0 (macOS)"
#: What the run's page and its report say about such a run.
NOTE = ("Docked on a Mac with ADCP 1.0 (ADFRsuite 1.0), because the Scripps "
        "Research Institute builds no ADCP 1.1 for macOS. Its scores are not "
        "comparable with runs docked with ADCP 1.1.")
#: Where `adcp-dock local` puts ADFRsuite 1.0 on a Mac, and the variable it
#: names it in for the server and its tasks.
DEFAULT = os.path.expanduser("~/.adcp-dock/ADFRsuite")
ENV = "ADCP_MACOS_ADFRSUITE"
#: ADFRsuite 1.0 for macOS: download URL, file name, SHA-256, bytes. Intel
#: binaries, which run on Apple silicon through Rosetta 2.
BUILD = ("https://ccsb.scripps.edu/adfr/download/1033/",
         "ADFRsuite_x86_64Darwin_1.0.tar.gz",
         "9556680809c76859b159b7709fa1a6e84716e2cb41a0ab8bd72aecbb7960b413",
         90475481)
PAGE = "https://ccsb.scripps.edu/adfr/downloads/"
BINARIES = ("adcp", "agfr", "prepare_receptor")
#: ADCP 1.0 builds side chains, and so writes an all-atom ranked pose, for
#: its ten best clusters only.
MAX_POSES = 10


def _system():
    return platform.system()


def is_mac():
    return _system() == "Darwin"


def installed(dest):
    return bool(dest) and all(os.path.isfile(os.path.join(dest, "bin", b))
                              for b in BINARIES)


def home():
    """ADFRsuite 1.0's folder on a Mac where it is installed, else None.
    Always None on any other system."""
    if not is_mac():
        return None
    for d in (os.environ.get(ENV), DEFAULT):
        if d and installed(d):
            return os.path.abspath(d)
    return None


def engaged():
    """Whether this machine docks with ADCP 1.0: a Mac with ADFRsuite 1.0."""
    return home() is not None


def program(name):
    """One of ADFRsuite 1.0's programs on a Mac, or None."""
    h = home()
    p = os.path.join(h, "bin", name) if h else None
    return p if p and os.path.isfile(p) and os.access(p, os.X_OK) else None


# --- installing it -------------------------------------------------------------

def install(dest=None, log=print):
    """ADFRsuite 1.0 into `dest` on a Mac. Returns `dest`. Safe to call
    again. Refuses on any other system."""
    if not is_mac():
        raise SystemExit("ADFRsuite 1.0 is installed only on a Mac, where the "
                         "ADCP 1.1 suite does not run; this is %s."
                         % _system())
    from .adcp11_install import _fetch, _sha256
    dest = os.path.abspath(os.path.expanduser(dest or DEFAULT))
    if installed(dest):
        log("ADFRsuite 1.0 is already installed in %s" % dest)
        return dest
    if platform.machine().lower() in ("arm64", "aarch64"):
        log("This Mac has Apple silicon. ADFRsuite is built for Intel and runs "
            "through Rosetta 2; if it is not installed, run "
            "`softwareupdate --install-rosetta` once.")
    url, name, sha, size = BUILD
    cache = os.path.dirname(dest)
    os.makedirs(cache, exist_ok=True)
    tar = os.path.join(cache, name)
    if not (os.path.isfile(tar) and _sha256(tar) == sha):
        log("  downloading %s from the Scripps Research Institute (%.0f MB)"
            % (name, size / 1e6))
        try:
            _fetch(url, tar, log)
        except Exception as exc:                                # noqa: BLE001
            raise SystemExit("Could not download ADFRsuite (%s). Get %s by "
                             "hand from %s, put it in %s, and run this again."
                             % (exc, name, PAGE, cache))
        got = _sha256(tar)
        if got != sha:
            os.remove(tar)
            raise SystemExit("The downloaded %s does not match its pinned "
                             "checksum and was deleted (got %s). Retry, or "
                             "download it by hand from %s." % (name, got, PAGE))
    log("  unpacking")
    with tarfile.open(tar) as t:
        top = t.getnames()[0].split("/")[0]
        for m in t.getmembers():
            p = os.path.normpath(m.name)
            if p.startswith("..") or os.path.isabs(p):
                raise SystemExit("refusing %s: a member escapes the archive"
                                 % m.name)
        if hasattr(tarfile, "tar_filter"):
            t.extractall(cache, filter="tar")
        else:
            t.extractall(cache)
    src = os.path.join(cache, top)
    # -l skips Scripps' academic-or-commercial question about the MSMS
    # library rather than answering it; stdin is /dev/null so nothing waits.
    log("  running the Scripps installer (a few minutes)")
    logf = os.path.join(cache, "adfr_install.log")
    with open(logf, "w") as fh:
        subprocess.run(["bash", "./install.sh", "-d", dest, "-c", "0", "-l"],
                       cwd=src, stdin=subprocess.DEVNULL, stdout=fh,
                       stderr=subprocess.STDOUT)
    if not installed(dest):
        raise SystemExit("ADFRsuite did not install; its own log is %s. Remove "
                         "%s before retrying." % (logf, dest))
    shutil.rmtree(src, ignore_errors=True)
    log("ADFRsuite 1.0 installed in %s" % dest)
    return dest


# --- the grid maps -------------------------------------------------------------

def maps(gradient=True):
    """How the maps are built here: ADFRsuite 1.0's agfr, AutoSite 1.0, and
    the receptor gradient as the run asks (agfr 1.0 adds it unless -ng)."""
    from . import target

    class Maps10(target.Maps):
        agfr = "1.0 (macOS)"

        def runner(self, tools=None):
            return [program("agfr")]

        def env(self):
            return None

        def flags(self):
            return [] if self.gradient else ["-ng"]

        def describe(self):
            return ("AGFR 1.0 (macOS), AutoSite 1.0, receptor gradient %s"
                    % ("on" if self.gradient else "off"))

    return Maps10(exe=None, autosite="1.0", gradient=gradient)


# --- docking, and the 1.1 layout straight after ----------------------------------

#: ADCP 1.0's mode table: mode, affinity, clust. rmsd, ref. rmsd, clust. size,
#: rmsd avg., best rmsd, energy avg., best run.
_HEADER = re.compile(r"^mode \|\s+affinity\s+\|\s+clust\.\s+\|\s+ref\.\s+\|\s+clust\.")
_ROW = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+\S+\s+\S+\s+(\d+)")
#: The loose replica files ADCP 1.0 leaves in the peptide's folder, and its
#: scratch beside them.
_LOOSE = re.compile(r"^out(_\d+)?\.(pdb|out)$")
_SCRATCH = ("constrains", "ramaprob.data", "tmp.pdb")


def mode_table(adcp_log):
    """[{rank, affinity, cluster}] from ADCP 1.0's adcp.log, as modes.json
    holds it, or [] when it holds no finished table."""
    seen, rows, ranks = False, [], set()
    try:
        fh = open(adcp_log, errors="replace")
    except OSError:
        return []
    with fh:
        for line in fh:
            if _HEADER.match(line):
                seen = True
                continue
            m = _ROW.match(line) if seen else None
            if not m:
                continue
            rank = int(m.group(1))
            if rank in ranks:          # a second copy of the table
                break
            ranks.add(rank)
            rows.append({"rank": rank, "affinity": float(m.group(2)),
                         "cluster": int(m.group(3))})
    return rows


def _loose(wd):
    return sorted(n for n in os.listdir(wd)
                  if (_LOOSE.match(n) or n in _SCRATCH)
                  and not n.startswith("out_ranked")
                  and os.path.isfile(os.path.join(wd, n)))


def failure_reason(wd):
    """What ADCP 1.0 said went wrong, in one line, or ''. Read before its
    replica logs are removed."""
    try:
        text = open(os.path.join(wd, "adcp.log"), errors="replace").read()
    except OSError:
        text = ""
    n_err = text.count("ENDED WITH ERROR")
    if n_err:
        return ("every replica failed" if n_err >= 20
                else "%d replica(s) failed" % n_err)
    for name in sorted(n for n in os.listdir(wd) if _LOOSE.match(n)
                       and n.endswith(".out"))[:3]:
        for line in open(os.path.join(wd, name), errors="replace"):
            if line.strip().startswith("ERROR!"):
                return line.strip()[len("ERROR!"):].strip()
    return ""


def to_layout(wd, key, keep, log=print):
    """A finished ADCP 1.0 docking folder put into the ADCP 1.1 layout:
    modes.json from adcp.log, and the loose replica files packed into
    <key>_trajectories.tar.gz (keep) or removed. Returns the modes."""
    from . import trajectories
    table = mode_table(os.path.join(wd, "adcp.log"))
    tmp = os.path.join(wd, "modes.json.part")
    with open(tmp, "w") as fh:
        json.dump(table, fh, indent=1)
    os.replace(tmp, os.path.join(wd, "modes.json"))
    names = _loose(wd)
    if keep and names:
        arc = trajectories.archive_path(wd, key)
        members = [n for n in names if _LOOSE.match(n) and n != "out.pdb"]
        sizes = {n: os.path.getsize(os.path.join(wd, n)) for n in members}
        part = arc + ".part"
        with tarfile.open(part, "w:gz") as tf:
            for n in members:
                tf.add(os.path.join(wd, n),
                       arcname=os.path.join(trajectories.WORKING, n))
        with tarfile.open(part, "r:gz") as tf:
            held = {m.name: m.size for m in tf.getmembers() if m.isfile()}
        if any(held.get(os.path.join(trajectories.WORKING, n)) != s
               for n, s in sizes.items()):
            os.remove(part)
            log("      the trajectory archive did not read back whole; the "
                "replica files are left as they are")
            return table
        os.replace(part, arc)
        log("  packed %d trajectory file(s) into %s" % (len(members),
                                                        os.path.basename(arc)))
    for n in names:
        try:
            os.remove(os.path.join(wd, n))
        except OSError:
            pass
    return table


def dock(sequence, trg, workdir, key=None, replicas=100, steps=5_000_000,
         cpus=0, log=print, helical=0, cyclic="none", keep=False, seed=None):
    """Run ADFRsuite 1.0's adcp for one peptide, then put its folder into the
    1.1 layout. Returns (modes, the peptide's docking folder)."""
    from . import adcp11
    exe = program("adcp")
    if not exe:
        raise SystemExit("ADCP 1.0 runs only on a Mac with ADFRsuite 1.0 "
                         "installed (adcp-dock local installs it)")
    wd = os.path.join(workdir, key or sequence)
    os.makedirs(wd, exist_ok=True)
    # ADCP 1.0 looks for the maps under ./tmp_<job>/<target name>/, so the
    # target is copied beside the job and passed by its bare name.
    local_trg = os.path.join(wd, os.path.basename(trg))
    src = os.path.abspath(trg)
    if (not os.path.isfile(local_trg)
            or os.path.getsize(local_trg) != os.path.getsize(src)):
        shutil.copy(src, local_trg)
    marker = os.path.join(wd, ".adcp_finished")
    if (os.path.isfile(marker) and os.path.isfile(os.path.join(wd, "modes.json"))
            and os.path.isfile(os.path.join(wd, "out_ranked_1.pdb"))):
        return adcp11.read_modes(wd), wd
    seq_arg, start = adcp11.start_arguments(sequence, helical)
    cmd = [exe, "-t", os.path.basename(trg), "-s", seq_arg, "-N", str(replicas),
           "-n", str(steps), "-O", "-o", "out"] + start
    cmd += adcp11.cyclic_flags(cyclic)
    if cpus:
        cmd += ["-c", str(cpus)]
    if seed is not None:
        cmd += ["-S", str(int(seed))]
    log("  adcp: %s" % " ".join(cmd[1:]))
    env = dict(os.environ)
    env.pop("OMP_NUM_THREADS", None)          # ADCP deadlocks if this is set
    with open(os.path.join(wd, "adcp.log"), "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=wd,
                           env=env)
    if not os.path.isfile(os.path.join(wd, "out_ranked_1.pdb")):
        why = failure_reason(wd)
        log("      %s: %s" % (sequence, "ADCP could not dock this at all: " + why
                              if why else "ADCP ranked no pose"))
    modes = to_layout(wd, key or sequence, keep, log=log)
    if r.returncode == 0 and modes:
        with open(marker, "w") as fh:
            fh.write("ADCP 1.0 (macOS): %d replicas, %d steps\n"
                     % (replicas, steps))
    return adcp11.read_modes(wd), wd


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m adcp_dock.mac10",
                                 description=__doc__.split("\n\n")[0])
    ap.add_argument("--dest", default=DEFAULT,
                    help="where to install it (default %s)" % DEFAULT)
    a = ap.parse_args(argv)
    install(a.dest)
    print("\nexport %s=%s" % (ENV, os.path.abspath(os.path.expanduser(a.dest))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
