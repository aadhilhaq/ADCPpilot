"""What the launcher's cluster box offers, as a table rather than in code.

By Aadhil Haq

"Slurm cluster" is any cluster with Slurm: the person types its login
address once and the launcher remembers it. "This computer" runs the work
locally instead.
"""
from __future__ import annotations

import os
import shlex

# Where "Report a problem" goes. One place, so that a fork has a single line
# to change.
FEEDBACK_TO = os.environ.get("ADCP_FEEDBACK_TO") or "aadhil.haq@tamu.edu"

#: Any cluster with Slurm, reached at the address the person types.
OTHER = "Slurm cluster"
#: No cluster at all: `adcp-dock local` on this machine (or in WSL).
THIS = "This computer"

SITES = {
    OTHER: {
        "host": "",
        "auth": "password",
        "install": "~/ADCP_docking",
        "runs": "~/adcp_runs",
        "hint": "Whatever you log in with",
        "label": "Username",
        "modules": ["Anaconda3", "Miniconda3"],
        "note": "Any cluster with Slurm. Install it there once: clone the repository and "
                "run ./install.sh with a folder that has room (INSTALL.md, 1.2); then type "
                "the cluster's login address and the clone's path above. Both are "
                "remembered.",
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
    """The names, the cluster first and this computer last."""
    return [OTHER, THIS]


def site_of(saved):
    """The site a saved settings file names, or the cluster when it names
    one this version does not offer, so that a remembered address and path
    still come back."""
    return saved if saved in SITES else OTHER


def _fill(template, user):
    """A path template, only once there is something to fill it with.

    "/home//adcp_runs" is worse than an empty box, because it looks like an
    answer.
    """
    if not template:
        return ""
    if "$USER" in template and not user:
        return ""
    return template.replace("$USER", user or "")


def defaults_for(site_name, user=""):
    """Every field the form can fill in by itself, for this site and account.

    The window calls it when the site changes and again when the account is
    typed, and writes whatever comes back into the boxes the person has not
    edited.
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
