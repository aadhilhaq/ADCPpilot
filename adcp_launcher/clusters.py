"""What differs between one cluster and the next, as a table rather than in code.

Selecting a site fills in the rest: a user should need to know their
account and password, not a group allocation path.

"Other (SLURM)" covers any other cluster; "This computer" runs the work
locally instead.
"""
from __future__ import annotations

import os
import shlex

# Where "Report a problem" goes. One place, so that a fork has a single line
# to change.
FEEDBACK_TO = os.environ.get("ADCP_FEEDBACK_TO") or "aadhil.haq@tamu.edu"

OTHER = "Other (SLURM)"
#: No cluster at all: `adcp-dock local` on this machine (or in WSL).
THIS = "This computer"

# Texas A&M HPRC asks for this acknowledgement in publications. The two
# presets below pass it to the server, so each run records it.
TAMU_ACK = ("The authors acknowledge the Texas A&M High Performance Research "
            "Computing (HPRC) facility for providing the computational "
            "resources used in this work. URL: https://hprc.tamu.edu")

SITES = {
    # Presets for two public HPC systems. The install path is where the
    # repository was cloned on that cluster; it is remembered after the
    # first connection, so it only has to be typed once.
    "Grace (TAMU)": {
        "host": "grace.hprc.tamu.edu",
        "auth": "password",
        "install": "~/ADCP_docking",
        "runs": "/scratch/user/$USER/adcp_runs",
        "hint": "Your cluster username",
        "label": "Username",
        "modules": ["Anaconda3/2024.02-1", "Anaconda3", "Miniconda3"],
        "note": "Password and a second-factor approval on every login.",
        "ack": TAMU_ACK,
    },
    "ACES (TAMU)": {
        "host": "aces.hprc.tamu.edu",
        "auth": "password",
        "install": "~/ADCP_docking",
        "runs": "/scratch/user/$USER/adcp_runs",
        "hint": "Your cluster username",
        "label": "Username",
        "modules": ["Anaconda3/2024.02-1", "Anaconda3", "Miniconda3"],
        "note": "Password and a second-factor approval on every login.",
        "ack": TAMU_ACK,
    },
    OTHER: {
        "host": "",
        "auth": "password",
        "install": "~/ADCP_docking",
        "runs": "~/adcp_runs",
        "hint": "Whatever you log in with",
        "label": "Username",
        "modules": ["Anaconda3", "Miniconda3"],
        "note": "Any cluster with SLURM and the tools installed.",
    },
    THIS: {
        "host": "",
        "auth": "none",
        "install": "",
        "runs": "~/adcp_runs",
        "hint": "",
        "label": "Username",
        "modules": [],
        "local": True,
        "note": "No cluster: the docking runs on this computer (on Windows, "
                "inside WSL). The boxes above are not used.",
    },
}


def sites():
    """The names, with the real clusters first and Other last."""
    return [n for n in SITES if n not in (OTHER, THIS)] + [OTHER, THIS]


def _fill(template, user):
    """A path template, only once there is something to fill it with.

    "/scratch/user//adcp_runs" is worse than an empty box, because it looks
    like an answer.
    """
    if not template:
        return ""
    if "$USER" in template and not user:
        return ""
    return template.replace("$USER", user or "")


def defaults_for(site_name, user=""):
    """Every field the form can fill in by itself, for this site and account.

    This is the whole of "only your account has to be typed". The window
    calls it when the site changes and again when the account is typed, and
    writes whatever comes back into the boxes the user has not edited.
    """
    site = SITES.get(site_name) or {}
    return {
        "host": site.get("host", ""),
        "note": site.get("note", ""),
        "hint": site.get("hint", ""),
        "label": site.get("label", "Username"),
        "auth": site.get("auth", "password"),
        "feedback": FEEDBACK_TO,
        "install": _fill(site.get("install", ""), user),
        "runs": _fill(site.get("runs", ""), user),
    }


# PORTS
#
# The remote port defaults to 0, meaning the server picks a free one and
# prints it, which is the only choice that cannot collide on a shared login
# node. The local port is derived from the username so that it is stable for
# one person and unlikely to collide between people.
AUTO = 0
PORT_BASE = 9100
PORT_SPAN = 800  # 9100 to 9899


def suggested_port(user):
    """A stable port for this person, in this tool's own band.

    Derived from the account name: uncoordinated, repeatable, and still
    editable. Mixed with the tool's name so that the same person gets a
    different number from another tool's launcher.
    """
    seed = (user or "adcp") + "|adcp-docking"
    return PORT_BASE + (sum(ord(c) * (i + 1) for i, c in enumerate(seed))
                        % PORT_SPAN)


def parse_port(text):
    """What the port box means. Blank or 'auto' -> let the server choose."""
    text = (text or "").strip().lower()
    if not text or text in ("auto", "0", "any"):
        return AUTO
    if text.isdigit() and 1024 <= int(text) <= 65535:
        return int(text)
    return AUTO


def remote_path(path):
    """A path for the far end's shell, with ~ still able to expand.

    shlex.quote('~/ADCP_docking') gives '~/ADCP_docking': in single quotes,
    which is exactly what stops the remote shell expanding the tilde. The
    launcher then ran `cd '~/ADCP_docking'` and the cluster answered "No
    such file or directory", naming a path that does exist. Quoting is still
    required for everything after it, because a path can contain a space.
    """
    path = (path or "").strip()
    if path == "~":
        return '"$HOME"'
    if path.startswith("~/"):
        return '"$HOME"/' + shlex.quote(path[2:])
    return shlex.quote(path)
