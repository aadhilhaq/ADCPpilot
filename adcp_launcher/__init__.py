"""The desktop launcher: open the docking tool from a laptop.

By Aadhil Haq

It makes an SSH connection to the cluster, starts the tool's own web server
there, forwards a local port to it and opens a browser. Nothing about the
tool has to be installed on the laptop except this launcher and paramiko.
"""

__version__ = "0.1.17"
# The same wording as the tool: ADCP and ADFRsuite belong to the Scripps
# Research Institute, and this package launches a GUI that drives them. The
# launcher is installed on computers that may not have adcp_dock, so it
# repeats these lines instead of importing them, and a test checks that the
# two copies match.
APP = "GUI to ADCP Docking"
ENGINE = "ADCP"
ENGINE_OWNER = "Scripps Research Institute"
BLURB = ("Peptide docking tool to be run in an interactive manner for "
         "users")
LAB = "HOPE Lab"
AUTHORS = ("Aadhil Haq",)
PI = "Dr. Sandun Fernando"
POWERED_BY = "Powered by HOPE Lab"


# The attribution sentence, word for word as the page's header carries it.
# What this is comes first; whose interface it is comes after.
ENGINE_CREDIT = ("A graphical interface to %s, a product of the %s. %s."
                 % (ENGINE, ENGINE_OWNER, POWERED_BY))

#: What to put in a paper. The launcher cannot import adcp_dock, it is
#: installed on laptops that have no copy of it, so it carries the
#: short form and points at the tool for the full list with DOIs.
HOW_TO_CITE = (
    "How to cite: ADCP (Zhang & Sanner 2019, "
    "doi:10.1093/bioinformatics/btz459) and this interface "
    "(github.com/aadhilhaq/ADCP_docking). Every finished run downloads "
    "with HOW_TO_CITE.txt, CITATION.bib and CITATION.ris listing what "
    "that run actually used, with DOIs.")



def credit(sep=" · "):
    """The credit line: the sentence, and nothing appended to it.

    It names the engine and its owner before the interface, because that is
    the order the claim has to be made in, and it stops there. Authors, PI
    and version have their own places (the page header, the launcher's own
    header); run onto the end of this they turn a statement of provenance
    into a byline, which is the thing it was rewritten to stop being.
    """
    return ENGINE_CREDIT
