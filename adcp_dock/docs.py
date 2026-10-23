"""The built-in documentation, served at /docs/ by the interface.

By Aadhil Haq

The layout, the stylesheet and the credit bar come from docskit, a module
shared with other tools by the same authors. This module holds only what is
this interface's own: the sidebar, the pages each "?" leads to, the
colours of this interface and the credit line.

A page is an HTML fragment in adcp_dock/docs/pages/<slug>.html: an <h1>, a
lead paragraph, then sections under <h2>. Screenshots are in
adcp_dock/docs/img/, taken of the running interface by
tools/docs_screenshots.py, so they show what the page really looks like.
"""
import os

from . import APP, __version__
from . import docskit

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "docs")

#: The sidebar, in order: (slug, title, section).
NAV = [
    ("index", "Overview", "Start here"),
    ("getting-started", "Getting started", "Start here"),
    ("first-run", "Your first run", "Start here"),

    ("receptor", "The receptor", "Preparing a run"),
    ("site", "The site and the box", "Preparing a run"),
    ("peptides", "The peptides", "Preparing a run"),
    ("docking", "Docking settings", "Preparing a run"),
    ("rescoring", "Rescoring with MM-GBSA", "Preparing a run"),
    ("submitting", "Resources and submitting", "Preparing a run"),

    ("runs", "Following your runs", "Results"),
    ("results", "Reading the results", "Results"),

    ("case-mdm2", "MDM2 and the p53 peptide", "Cases"),
    ("case-two-chains", "A site between two chains", "Cases"),
    ("case-histidines", "Histidines at the site", "Cases"),
    ("case-modified", "Modified and D amino acids", "Cases"),

    ("tools", "Tools", "Reference"),
    ("troubleshooting", "When something goes wrong", "Reference"),
    ("citing", "How to cite", "Reference"),
]

#: Where each "?" on the page leads: the help key in helptext.HELP, then the
#: page and the section that explains that setting at length.
HELP_LINKS = {
    "reccode": "receptor.html#loading",
    "recpath": "receptor.html#loading",
    "ph": "receptor.html#protonation",
    "protonation": "receptor.html#protonation",
    "sitesurf": "site.html#surface",
    "minimise": "receptor.html#minimisation",
    "keep_hetatms": "receptor.html#heteroatoms",
    "rebuild_sidechains": "receptor.html#side-chains",
    "fill_loops": "receptor.html#loops",
    "targetfile": "receptor.html#target-file",
    "adaptframe": "receptor.html#adapt",
    "hotspots": "site.html#residues",
    "site_mode": "site.html#where-the-box-goes",
    "site_ligand": "site.html#ligand",
    "site_pocket": "site.html#pocket",
    "autosite": "site.html#autosite",
    "receptor_gradient": "site.html#gradient",
    "cyclic": "peptides.html#cyclic",
    "helical": "peptides.html#helical",
    "keep_trajectories": "docking.html#trajectories",
    "keep_files": "runs.html#tidy",
    "pad": "site.html#box-size",
    "bmin": "site.html#box-size",
    "bmax": "site.html#box-size",
    "showbox": "site.html#seeing-the-box",
    "box_chains": "site.html#other-chains",
    "pastetext": "peptides.html#formats",
    "poses": "docking.html#poses",
    "replicas": "docking.html#replicas",
    "steps": "docking.html#steps",
    "mmgbsa": "rescoring.html#what-it-does",
    "mmposes": "rescoring.html#poses-to-rescore",
    "topn": "rescoring.html#averaging",
    "cpus": "submitting.html#cores",
    "walltime": "submitting.html#walltime",
    "throttle": "submitting.html#throttle",
    "label": "submitting.html#naming",
    "outdir": "submitting.html#output-folder",
    "account": "submitting.html#account",
    "sbatch_extra": "submitting.html#slurm-options",
}

#: Line 3 of LICENSE, word for word. tests/test_docs.py holds the two
#: together, because the bar along the bottom of every page states it.
COPYRIGHT = ("Copyright (c) 2026, Aadhil Haq and Sandun Fernando, "
             "HOPE Lab, Texas A&M University")

#: Who does the work, for the bar along the bottom. A bar has room only for
#: the programs every run depends on and the two viewers a figure comes from;
#: the full list, scoped to one run, is behind How to cite.
CREDITS = [
    'docking by <a href="https://ccsb.scripps.edu/adcp/" target="_blank"'
    ' rel="noopener">ADCP</a>, AGFR and AutoSite (Scripps Research)',
    "receptor hydrogens by reduce (Word et al., J Mol Biol 1999), or "
    "PDB2PQR with PROPKA when chosen",
    "MM-GBSA by AmberTools",
    "viewers 3Dmol.js (Rego &amp; Koes, Bioinformatics 2015) and Mol* "
    "(Sehnal et al., Nucleic Acids Res 2021), surfaces by MSMS (Sanner et "
    "al., Biopolymers 1996)",
]

#: The interface's own colours, from the :root block in page.py, under the
#: names the kit uses. Everything else (warnings, code, marks) keeps the
#: kit's own values.
LIGHT = {"bg": "#f7f8fb", "card": "#fff", "line": "#dfe3ee", "fg": "#1b2333",
         "dim": "#66718a", "accent": "#1a56db", "accent2": "#e8efff"}
DARK = {"bg": "#0f1219", "card": "#161b26", "line": "#28303f", "fg": "#e7ecf5",
        "dim": "#94a0b8", "accent": "#5b8cf0", "accent2": "#1d2942"}

#: Everything the shared layout needs to know about these docs.
SITE = docskit.Site(
    app=APP,
    tagline="a graphical interface to ADCP",
    version=__version__,
    copyright=COPYRIGHT,
    nav=NAV,
    root=HERE,
    help_links=HELP_LINKS,
    credits=CREDITS,
    # The key app.js stores the light or dark choice under, so the docs open
    # in the theme the interface is in.
    theme_key="adcp.theme",
    light=LIGHT,
    dark=DARK,
    shots="tools/docs_screenshots.py",
    # The bar's How to cite opens the interface's own citation dialog rather
    # than a page, so it is a link to nowhere that app.js intercepts.
    cite_href="#",
)


def render(slug):
    """One page in the shared layout, or None when there is no such page."""
    return docskit.render(SITE, slug)


def respond(path):
    """(status, content type, body, location) for a request under /docs/."""
    return docskit.respond(SITE, path)


def audit():
    """Everything wrong with the pages, as a list of sentences."""
    return docskit.audit(SITE)
