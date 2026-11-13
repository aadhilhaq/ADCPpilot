"""What the launcher does when started, however it was started.

By Aadhil Haq

Kept out of __main__.py so there is one implementation whether it is reached
by `python -m adcp_launcher`, by the console script, or by a bundled
executable running a script that imports it.
"""
import sys


def main(argv=None):
    argv = sys.argv if argv is None else argv

    if "--version" in argv:
        from . import credit
        print(credit())
        return 0

    # --selftest is for the build. A windowed executable cannot be driven
    # without a display, but it can be asked to prove that it starts, that
    # what it needs came along in the bundle, and that its imports resolve.
    if "--selftest" in argv:
        from . import session, tunnel, credit                  # noqa: F401
        line = session.command("/x/ADCP_docking", runs="~/r")
        assert "adcp_pilot.cli serve" in line, line
        assert "--host 127.0.0.1" in line, line
        assert tunnel.AUTH_PATIENCE >= 120, "the two-factor timeout was shortened"
        assert session.URL.search("http://127.0.0.1:8123/?t=abcDEF_-"), \
            "the server's URL line would not be recognized"
        # The window needs a display and paramiko needs installing,
        # neither of which is true on a cluster. Both are reported
        # rather than required, so this check can run anywhere.
        notes = []
        try:
            import paramiko
            notes.append("paramiko " + paramiko.__version__)
        except ImportError:
            notes.append("paramiko NOT INSTALLED: the launcher cannot "
                         "connect until it is")
        try:
            from . import app
            assert hasattr(app, "App") and hasattr(app, "main"), "no window"
            notes.append("window ok")
        except ImportError as exc:
            notes.append("no tkinter here (%s)" % exc)
        # A stable last line, so that a build can confirm the
        # launcher started without matching the product name, which
        # changes.
        print("%s\n  %s\nSELFTEST OK" % (credit(), "\n  ".join(notes)))
        return 0

    from .app import main as run
    run()
    return 0
