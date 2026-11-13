"""Packaging a run: as a zip to download, or as mail when it failed.

By Aadhil Haq

A run finishes on a cluster and its results stay there. This produces a zip
of what is worth keeping, and a much smaller bundle of logs and settings
for a problem report.
"""

from __future__ import annotations

import fnmatch
import io
import os
import time
import zipfile

#: The files worth keeping from a finished run, in the order a reader wants
#: them. Globs, matched against the path relative to the run directory.
RESULTS = (
    # Citation files first, so that a download carries what to cite with
    # it: the method section is written from this folder.
    "HOW_TO_CITE.txt", "CITATION.bib", "CITATION.ris",
    "job.json",
    "results.tsv", "results.csv", "*.tsv", "*.csv",
    "results.json",
    # The interface panel and its residue maps. results.tsv holds one row per
    # peptide; these hold every pose that was measured and every residue at
    # the interface, which is what a figure is drawn from.
    "metrics/*.json", "metrics/*_residues.tsv",
    "report/*",
    "receptor/prepared.pdb", "receptor/preparation.log",
    # The structures somebody opens. Without these a download is tables and
    # logs, and the poses are only ADCP's own out_ranked_*.pdb with no
    # receptor beside them.
    "poseviewer/*.pdb", "poseviewer/*.sdf",
    "poseviewer/*.mae", "poseviewer/*.maegz",
    "docking/*/out_ranked_*.pdb", "docking/*/modes.json",
    "docking/*/adcp.log",
    "mmgbsa/*/*/FINAL_RESULTS*", "mmgbsa/*/*.csv",
    # The MM-GBSA output the energy terms are read from, and the per residue
    # decomposition, so the panel can be checked against its own source.
    "mmgbsa/*/*/*_mmgbsa.dat", "mmgbsa/*/*/*_decomp.dat",
    "mmgbsa/*/*/*_mmpbsa.in", "mmgbsa/*/*/*_min.out",
    "mmgbsa/*/*/*_scored_on.txt", "mmgbsa/RADII",
    "logs/*.out", "logs/*.err",
    "*.sbatch",
)

# What a problem report needs, which is different: no poses, no structures,
# nothing large; the logs, the settings, and what the scheduler reported.
LOGS = (
    "job.json",
    "logs/*.out", "logs/*.err",
    "*.sbatch",
    "receptor/preparation.log",
    "receptor/protonated.pdb2pqr.txt",
    "docking/*/adcp.log",
    "mmgbsa/*/*/*.log",
    "target/*.log",
)

# A single file larger than this is truncated into a report rather than
# dropped: the end of a log is the part that says what went wrong.
LOG_TAIL_BYTES = 200_000

# A whole problem report. Beyond this, mail servers begin refusing.
LOGS_TOTAL_BYTES = 6_000_000


def _matches(rel, patterns):
    return any(fnmatch.fnmatch(rel, p) for p in patterns)


def walk(rundir, patterns=None):
    """(relative path, absolute path) for the files a bundle should hold.

    `patterns` of None means everything, which is what "download all" is.
    """
    out = []
    rundir = os.path.abspath(rundir)
    for base, dirs, names in os.walk(rundir):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        for n in sorted(names):
            if n.startswith("."):
                continue
            path = os.path.join(base, n)
            if os.path.islink(path) or not os.path.isfile(path):
                continue
            rel = os.path.relpath(path, rundir)
            if patterns is None or _matches(rel, patterns):
                out.append((rel, path))
    return out


def zip_bytes(rundir, what="results", name=None):
    """The bundle, in memory. Returns (filename, bytes).

    In memory because the alternative is writing a zip into the run directory
    and then deciding when to delete it.
    """
    patterns = None if what == "all" else RESULTS
    files = walk(rundir, patterns)
    label = name or os.path.basename(os.path.normpath(rundir)) or "run"
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for rel, path in files:
            try:
                z.write(path, os.path.join(label, rel))
            except OSError:
                continue                    # a file that vanished mid-walk
        if not files:
            z.writestr(os.path.join(label, "EMPTY.txt"),
                       "This run directory held no files matching %r.\n"
                       % what)
    return "%s_%s.zip" % (label, what), buf.getvalue()


def _tail(path, cap=LOG_TAIL_BYTES):
    """The last `cap` bytes, with a line saying so. None if unreadable."""
    try:
        size = os.path.getsize(path)
        with open(path, "rb") as fh:
            if size > cap:
                fh.seek(size - cap)
                head = (b"[ earlier %d bytes not included; this is the last "
                        b"%d ]\n" % (size - cap, cap))
                return head + fh.read()
            return fh.read()
    except OSError:
        return None


def logs_for(rundir):
    """[(name, bytes)] for a problem report: the logs, capped and tailed."""
    out, total = [], 0
    for rel, path in walk(rundir, LOGS):
        data = _tail(path)
        if data is None:
            continue
        if total + len(data) > LOGS_TOTAL_BYTES:
            continue                        # skip this one, try the rest
        total += len(data)
        out.append((rel.replace(os.sep, "_"), data))
    return out


def sendmail():
    """The local mail transfer agent, or None."""
    from shutil import which
    for p in ("/usr/sbin/sendmail", "/usr/lib/sendmail"):
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    return which("sendmail")


def can_send():
    return bool(sendmail())


def send(to, subject, body, attachments=(), sender=None, cc=None,
         reply_to=None):
    """Post a problem report. Returns (ok, what happened).

    Never raises: a report that cannot be sent must say so on the page, so the
    user can download the bundle and attach it instead.
    """
    import subprocess
    from email.message import EmailMessage

    exe = sendmail()
    if not exe:
        return False, ("this machine has no mail transfer agent, so nothing "
                       "can be posted from it")
    msg = EmailMessage()
    msg["To"] = to
    if cc:
        msg["Cc"] = cc
    msg["From"] = sender or "%s@%s" % (os.environ.get("USER", "adcp"),
                                       os.uname().nodename)
    # Where a reply goes. The envelope sender is <user>@<compute-node>,
    # which takes no mail and identifies no one, so the reply-to is set
    # explicitly.
    if reply_to:
        msg["Reply-To"] = reply_to
    msg["Subject"] = subject
    msg.set_content(body)
    for name, data in attachments:
        msg.add_attachment(data, maintype="application", subtype="octet-stream",
                           filename=name)
    try:
        p = subprocess.run([exe, "-t", "-oi"], input=msg.as_bytes(),
                           capture_output=True, timeout=120)
    except (OSError, subprocess.SubprocessError) as exc:
        return False, "the mail command failed: %s" % exc
    if p.returncode != 0:
        return False, ("the mail command exited %d: %s"
                       % (p.returncode,
                          (p.stderr or b"").decode("utf-8", "replace")[:300]))
    return True, "sent to %s" % to


def report_body(rundir, note, extra=None):
    """The message, with the details that are easily omitted."""
    extra = dict(extra or {})
    # Taken out of `extra` rather than printed alongside it, where it
    # appeared twice.
    who = (extra.pop("email", "") or "").strip()
    lines = [note.strip() or "(no note given)", "",
             "--- this part is automatic ---",
             ("from:     %s" % who) if who else
             "from:     no address given, so this report cannot be answered",
             "run:      %s" % rundir,
             "when:     %s" % time.strftime("%Y-%m-%d %H:%M:%S"),
             "host:     %s" % os.uname().nodename,
             "user:     %s" % os.environ.get("USER", "?")]
    for k, v in sorted(extra.items()):
        lines.append("%-9s %s" % (k + ":", v))
    files = logs_for(rundir)
    lines += ["", "attached: %d file(s), %s"
              % (len(files), ", ".join(n for n, _d in files[:12]) or "none")]
    return "\n".join(lines) + "\n"


def sizes(rundir):
    """{"results": (files, bytes), "all": (files, bytes)} for one run.

    Asked before anything is clicked: "download everything" on a run of
    several gigabytes over a tunneled connection is worth knowing in advance.
    """
    n_all = b_all = n_res = b_res = 0
    for rel, path in walk(rundir):
        try:
            size = os.path.getsize(path)
        except OSError:
            continue
        n_all += 1
        b_all += size
        if _matches(rel, RESULTS):
            n_res += 1
            b_res += size
    return {"results": (n_res, b_res), "all": (n_all, b_all)}


def human_size(n):
    """Bytes in human-readable form. GB above a thousand megabytes."""
    n = int(n or 0)
    if n >= 1e9:
        return "%.1f GB" % (n / 1e9)
    if n < 1e6:
        return "<1 MB" if n else "empty"
    return "%d MB" % round(n / 1e6)
