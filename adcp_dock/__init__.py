"""GUI to ADCP Docking: peptide docking with ADCP, on a SLURM cluster or locally.

By Aadhil Haq

A receptor structure and a list of peptide sequences go in; prepared
structures, docked poses, scores and optional MM-GBSA energies come out.
Computational stages run as batch jobs, or in order on the local machine
when local execution is selected.

Main modules: pipeline (the stages), server and page (the web interface),
cli (the command line), slurm (scheduling and local execution), cite
(citations), adcp11 and adcp11_install (the ADCP 1.1 suite, and installing
it). The preparation, docking and rescoring engine is in _vendor.
"""

__version__ = "0.1.21"

# Authorship and attribution, defined once and read by the launcher, the
# page, the CLI and the README. ADCP and the suite it comes in are products
# of the Scripps Research Institute; this package is an interface to them.
APP = "GUI to ADCP Docking"
ENGINE = "ADCP"
ENGINE_OWNER = "Scripps Research Institute"
ENGINE_URL = "https://ccsb.scripps.edu/adcp/"
BLURB = ("Peptide docking tool to be run in an interactive manner for "
         "users")
LAB = "HOPE Lab"
AUTHORS = ("Aadhil Haq", "Sandun Fernando")
PI = "Dr. Sandun Fernando"
POWERED_BY = "Powered by HOPE Lab"
# A computing facility may require an acknowledgment in publications.
# Which one applies depends on where the work ran, so it is configuration:
# ADCP_ACKNOWLEDGEMENT holds the facility's text, and each run records it.
GENERIC_ACK = ("Acknowledge the computing facility that ran this work, if "
               "its policy requires it.")


def acknowledgement(job=None):
    """The acknowledgment recorded for a run, else the one set in
    ADCP_ACKNOWLEDGEMENT, else a generic reminder."""
    import os
    text = ((job or {}).get("acknowledgement")
            or os.environ.get("ADCP_ACKNOWLEDGEMENT", "")).strip()
    return text or GENERIC_ACK
REPO = "https://github.com/aadhilhaq/ADCP_docking"


# The attribution sentence. engine_credit() adds the disclaimer.
ENGINE_CREDIT = ("A graphical interface to %s, a product of the %s. %s."
                 % (ENGINE, ENGINE_OWNER, POWERED_BY))


def credit(sep=" · "):
    """The credit line shown in the interface: the attribution sentence alone.
    """
    return ENGINE_CREDIT


# The disclaimer carried by the page footer.
DISCLAIMER = "Not affiliated with or endorsed by Scripps."


def engine_credit():
    """The attribution sentence and disclaimer, for wherever the tool is shown or published.
    """
    return ENGINE_CREDIT + " " + DISCLAIMER
