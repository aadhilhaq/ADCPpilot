"""Each replica's trajectory, kept as one archive per peptide or not at all.

By Aadhil Haq

ADCP 1.1 searches in a working folder inside the peptide's docking folder
(out/, from -w . -o out) and deletes it when the docking ends, unless it is
given -k. ADCP 1.0 has no such switch and leaves every replica's trajectory
and log loose in the docking folder (out_<n>.pdb and out_<n>.out, with out.pdb
holding them all again).

Either way that is hundreds of files per peptide, and on a cluster with a
quota on files rather than bytes a screen of fifty peptides can fill it. So
what a run keeps is decided by one setting, keep_trajectories:

  on   each replica's trajectory and logs are packed into
       <key>_trajectories.tar.gz in the peptide's docking folder, and the
       loose files are removed (see members() for the two copies left out);
  off  ADCP 1.1 removes its own working folder; ADCP 1.0's loose trajectories
       are removed here;
  absent (a job file written before the setting existed): nothing is touched,
       which is what those runs did.

The archive is written under another name and renamed into place, then a
marker is written, then the loose files are removed, so a task killed at any
point is finished by the next attempt without losing anything and without
packing a half-removed folder.
"""
import os
import re
import shutil
import tarfile

#: ADCP 1.1's working folder, inside the peptide's docking folder.
WORKING = "out"
#: ADCP 1.0's loose trajectory files.
_LOOSE10 = re.compile(r"^out(_\d+)?\.(pdb|out)$")
#: What else ADCP 1.0 leaves beside them, and nothing reads again: the
#: restraint file it writes, its copy of the Ramachandran table (a megabyte
#: per peptide) and the last pose its clustering built side chains on.
_SCRATCH10 = ("constrains", "ramaprob.data", "tmp.pdb")
#: Written once the archive is complete; the loose files may then go.
PACKED = ".trajectories_packed"


def archive_name(key):
    return "%s_trajectories.tar.gz" % key


def archive_path(wd, key):
    return os.path.join(wd, archive_name(key))


def loose(wd, engine):
    """The trajectory files and folders a docking left in `wd`, relative to
    it. The ranked poses, ADCP's log and anything this package wrote are
    never among them."""
    if not os.path.isdir(wd):
        return []
    if engine == "1.1":
        return [WORKING] if os.path.isdir(os.path.join(wd, WORKING)) else []
    return sorted(n for n in os.listdir(wd)
                  if (_LOOSE10.match(n) or n in _SCRATCH10)
                  and not n.startswith("out_ranked")
                  and os.path.isfile(os.path.join(wd, n)))


def members(wd, engine):
    """[(path on disk, name in the archive)] for what an archive keeps: each
    replica's trajectory and its logs, under out/ whichever ADCP wrote them.

    Two things in the working folder are left out because they are copies:
    out.pdb, every replica's trajectory again in one file, and, for ADCP 1.1,
    the target it unpacked there (about 100 MB of maps for a 30 A box, the
    same maps as the run's own .trg). runADCP.py makes both again when it is
    asked to recluster (-r): it unpacks the target before it looks at -r,
    and concatenates the replicas after.
    """
    out = []
    if engine == "1.1":
        top = os.path.join(wd, WORKING)
        if os.path.isdir(top):
            for n in sorted(os.listdir(top)):
                p = os.path.join(top, n)
                if os.path.isfile(p) and n != "out.pdb":
                    out.append((p, os.path.join(WORKING, n)))
        return out
    for n in loose(wd, engine):
        if _LOOSE10.match(n) and n != "out.pdb":
            out.append((os.path.join(wd, n), os.path.join(WORKING, n)))
    return out


def _remove(wd, names):
    for n in names:
        p = os.path.join(wd, n)
        if os.path.isdir(p) and not os.path.islink(p):
            shutil.rmtree(p, ignore_errors=True)
        else:
            try:
                os.remove(p)
            except OSError:
                pass


def finish(wd, key, engine, keep, log=print):
    """Pack, remove or leave a finished peptide's trajectories.

    Returns the archive's path when one was written or is already there,
    and None otherwise. `keep` is the job's keep_trajectories: True, False,
    or None for a job file that predates it.
    """
    if keep is None:
        return None
    names = loose(wd, engine)
    arc = archive_path(wd, key)
    marker = os.path.join(wd, PACKED)
    if not keep:
        if names:
            _remove(wd, names)
            log("  removed %d trajectory file(s); the run was asked not to "
                "keep them" % _count(names, wd))
        return arc if os.path.isfile(arc) else None
    if os.path.isfile(marker) and os.path.isfile(arc):
        # Packed by an earlier attempt that stopped before it had removed
        # everything: what is left is already in the archive.
        _remove(wd, names)
        return arc
    if not names:
        return arc if os.path.isfile(arc) else None
    keep_these = members(wd, engine)
    part = arc + ".part"
    with tarfile.open(part, "w:gz") as tf:
        for path, arcname in keep_these:
            tf.add(path, arcname=arcname)
    os.replace(part, arc)
    with open(marker, "w") as fh:
        fh.write("%d file(s) packed into %s\n"
                 % (len(keep_these), os.path.basename(arc)))
    n = _count(names, wd)
    _remove(wd, names)
    others = n - len(keep_these)
    log("  packed %d trajectory file(s) into %s (%.1f MB) and removed them%s"
        % (len(keep_these), os.path.basename(arc), os.path.getsize(arc) / 1e6,
           ", with %d other file(s) ADCP left beside them" % others
           if others > 0 else ""))
    return arc


def forget(wd):
    """Before a docking runs again in `wd`: the archive left by an earlier
    docking no longer describes what is about to be written."""
    for n in (PACKED,):
        try:
            os.remove(os.path.join(wd, n))
        except OSError:
            pass


def _count(names, wd):
    n = 0
    for name in names:
        p = os.path.join(wd, name)
        if os.path.isdir(p):
            n += sum(len(f) for _b, _d, f in os.walk(p))
        else:
            n += 1
    return n


def kept(wd, key):
    """The archive for a peptide, or None."""
    p = archive_path(wd, key)
    return p if os.path.isfile(p) else None
