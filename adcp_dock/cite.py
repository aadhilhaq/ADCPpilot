"""`adcp-dock cite`: what to cite, for this tool and everything it runs.

By Aadhil Haq

A method section is usually written months after the run, from the output
directory rather than the repository, so the citation list is generated
into each run alongside the results. It is scoped to what that run actually
used: a run without MM-GBSA is not given an Amber citation to include.

Every DOI here is checked against Crossref or DataCite by `--verify`.
"""
from __future__ import annotations

from . import peptide
from . import (APP, AUTHORS, ENGINE, ENGINE_OWNER, ENGINE_URL, LAB,
               PI, REPO, __version__, acknowledgement)

# (key, what it did, citation, doi). `when` is the condition under which it
# must be cited: "always", or the setting that enables it.
REFERENCES = [
    ("adcp", "always", "the docking search itself",
     "Zhang Y, Sanner MF. AutoDock CrankPep: combining folding and "
     "docking to predict protein-peptide complexes. Bioinformatics "
     "35(24):5121-5127, 2019.",
     "10.1093/bioinformatics/btz459"),
    ("adfr", "always", "grid maps (agfr) and receptor preparation "
     "(prepare_receptor)",
     "Ravindranath PA, Forli S, Goodsell DS, Olson AJ, Sanner MF. "
     "AutoDockFR: advances in protein-ligand docking with explicitly "
     "specified binding site flexibility. PLoS Comput Biol 11(12):e1004586, "
     "2015.",
     "10.1371/journal.pcbi.1004586"),
    # agfr computes the affinity maps with AutoGrid4 (ADFR/utils/MakeGrids.py).
    ("autogrid", "always", "the affinity maps agfr computes with AutoGrid4",
     "Morris GM, Huey R, Lindstrom W, Sanner MF, Belew RK, Goodsell DS, "
     "Olson AJ. AutoDock4 and AutoDockTools4: automated docking with "
     "selective receptor flexibility. J Comput Chem 30(16):2785-2791, 2009.",
     "10.1002/jcc.21256"),
    # AutoSite is run by agfr on every target to place the translational
    # points ADCP searches from, and agfr prints its own request to cite
    # both papers.
    ("autosite", "always", "the translational points agfr computes for "
     "every grid target",
     "Ravindranath PA, Sanner MF. AutoSite: an automated approach for "
     "pseudo-ligands prediction—from ligand-binding sites identification "
     "to predicting key ligand atoms. Bioinformatics 32(20):3142-3149, 2016.",
     "10.1093/bioinformatics/btw367"),
    ("reduce", "receptor preparation with reduce (the default), a receptor "
     "from ADapt, or MM-GBSA in such a run",
     "the receptor's hydrogens: the histidine protons and the Asn, Gln and "
     "His flips",
     "Word JM, Lovell SC, Richardson JS, Richardson DC. Asparagine and "
     "glutamine: using hydrogen atom contacts in the choice of side-chain "
     "amide orientation. J Mol Biol 285(4):1735-1747, 1999.",
     "10.1006/jmbi.1998.2401"),
    ("pdb2pqr", "receptor preparation with PROPKA, or MM-GBSA in such a run",
     "protonation states and titration at the chosen pH",
     "Jurrus E, Engel D, Star K, et al. Improvements to the APBS "
     "biomolecular solvation software suite. Protein Sci 27(1):112-128, 2018.",
     "10.1002/pro.3280"),
    ("propka", "receptor preparation with PROPKA, or MM-GBSA in such a run",
     "the pKa prediction pdb2pqr uses",
     "Olsson MHM, Sondergaard CR, Rostkowski M, Jensen JH. PROPKA3: "
     "consistent treatment of internal and surface residues in empirical "
     "pKa predictions. J Chem Theory Comput 7(2):525-537, 2011.",
     "10.1021/ct100578z"),
    ("pdbfixer", "receptor preparation, or ADCP 1.1's minimization",
     "rebuilding missing side chains and atoms",
     "Eastman P, et al. OpenMM 7: rapid development of high performance "
     "algorithms for molecular dynamics. PLoS Comput Biol 13(7):e1005659, "
     "2017.  (PDBFixer is distributed with OpenMM.)",
     "10.1371/journal.pcbi.1005659"),
    ("openmm", "the receptor's minimization, or ADCP 1.1's",
     "the restrained energy minimisation",
     "Eastman P, et al. OpenMM 7: rapid development of high performance "
     "algorithms for molecular dynamics. PLoS Comput Biol 13(7):e1005659, "
     "2017.",
     "10.1371/journal.pcbi.1005659"),
    ("amber", "mmgbsa", "tleap, sander and the force field used for MM-GBSA",
     "Case DA, et al. AmberTools. J Chem Inf Model 63(20):6183-6191, 2023.",
     "10.1021/acs.jcim.3c01153"),
    # Asked for by https://ambermd.org/CiteAmber.php with the force field.
    # The receptor's minimization uses it too (OpenMM's amber14-all.xml).
    ("ff14sb", "the receptor's minimization, or MM-GBSA",
     "the ff14SB protein force field of the receptor's minimization and of "
     "MM-GBSA",
     "Maier JA, Martinez C, Kasavajhala K, Wickstrom L, Hauser KE, "
     "Simmerling C. ff14SB: improving the accuracy of protein side chain and "
     "backbone parameters from ff99SB. J Chem Theory Comput 11(8):3696-3713, "
     "2015.",
     "10.1021/acs.jctc.5b00255"),
    ("mmpbsa", "mmgbsa", "the MM-GBSA calculation itself",
     "Miller BR III, McGee TD Jr, Swails JM, Homeyer N, Gohlke H, Roitberg "
     "AE. MMPBSA.py: an efficient program for end-state free energy "
     "calculations. J Chem Theory Comput 8(9):3314-3321, 2012.",
     "10.1021/ct300418h"),
    ("gbmodel", "mmgbsa",
     "the generalised Born model itself (igb=5, GB-OBC II, 0.15 M salt)",
     "Onufriev A, Bashford D, Case DA. Exploring protein native states and "
     "large-scale conformational changes with a modified generalized born "
     "model. Proteins 55(2):383-394, 2004.",
     "10.1002/prot.20033"),
    # ADCP 1.1's non-standard residues and their rotamers come from the
    # SwissSidechain database, and its OpenMM parameters for them are GAFF
    # (parmchk2 output, GAFF atom types), as its documentation describes.
    ("swisssidechain", "peptides with non-standard residues (ADCP 1.1)",
     "the non-standard side chains and rotamers in ADCP 1.1's library",
     "Gfeller D, Michielin O, Zoete V. SwissSidechain: a molecular and "
     "structural database of non-natural sidechains. Nucleic Acids Res "
     "41(D1):D327-D332, 2013.",
     "10.1093/nar/gks991"),
    ("gaff", "OpenMM minimization of peptides with non-standard residues",
     "the GAFF parameters ADCP 1.1 minimizes non-standard residues with",
     "Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. Development and "
     "testing of a general amber force field. J Comput Chem 25(9):1157-1174, "
     "2004.",
     "10.1002/jcc.20035"),
    # ADCP 1.1 minimizes with OpenMM using Amber ff99SB (amber99sb.xml) and,
    # in implicit solvent, the GBn2 model (gbn2.xml): openMMmethods.py. The
    # receptor's minimization uses GBn2 as well (implicit/gbn2.xml).
    ("ff99sb", "ADCP 1.1 poses minimized with OpenMM",
     "the Amber ff99SB force field ADCP 1.1 minimizes poses with",
     "Hornak V, Abel R, Okur A, Strockbine B, Roitberg A, Simmerling C. "
     "Comparison of multiple Amber force fields and development of improved "
     "protein backbone parameters. Proteins 65(3):712-725, 2006.",
     "10.1002/prot.21123"),
    ("gbn2", "the receptor's minimization, ADCP 1.1's in implicit solvent, "
     "or MM-GBSA run with GBn2",
     "the GBn2 implicit solvent model",
     "Nguyen H, Roe DR, Simmerling C. Improved generalized Born solvent model "
     "parameters for protein simulations. J Chem Theory Comput "
     "9(4):2020-2034, 2013.",
     "10.1021/ct3010485"),
    ("msms", "if a figure with a surface from the viewer is used",
     "the molecular surfaces the viewers draw, computed by MSMS",
     "Sanner MF, Olson AJ, Spehner JC. Reduced surface: an efficient way to "
     "compute molecular surfaces. Biopolymers 38(3):305-320, 1996.",
     "10.1002/(SICI)1097-0282(199603)38:3<305::AID-BIP4>3.0.CO;2-Y"),
    ("3dmol", "if a figure from the viewer is used",
     "the interactive structure viewer",
     "Rego N, Koes D. 3Dmol.js: molecular visualization with WebGL. "
     "Bioinformatics 31(8):1322-1324, 2015.",
     "10.1093/bioinformatics/btu829"),
    # Mol* is offered beside 3Dmol.js as a viewer, so a figure made with
    # it needs its own citation.
    ("molstar", "if a figure from the Mol* viewer is used",
     "the Mol* structure viewer",
     "Sehnal D, Bittrich S, Deshpande M, et al. Mol* Viewer: modern web app "
     "for 3D visualization and analysis of large biomolecular structures. "
     "Nucleic Acids Res 49(W1):W431-W437, 2021.",
     "10.1093/nar/gkab314"),
    ("numpy", "always", "the array arithmetic under the box and the poses",
     "Harris CR, Millman KJ, van der Walt SJ, et al. Array programming with "
     "NumPy. Nature 585(7825):357-362, 2020.",
     "10.1038/s41586-020-2649-2"),
    ("rdkit", "always", "the peptide topologies every pose is laid onto",
     "RDKit: Open-source cheminformatics. https://www.rdkit.org",
     "10.5281/zenodo.591637"),
]

# DOIs registered with DataCite rather than Crossref, mostly Zenodo software
# archives. Crossref's API does not know them, so verify() asks DataCite.
DATACITE = {"10.5281/zenodo.591637"}

# The same references as structured fields. The prose strings above are what
# a person reads; these are what a reference manager and a typesetter need.
# They are kept alongside rather than parsed out of the prose, because
# splitting an author list out of a prose citation produces entries that look
# correct and are not.
BIB = {
    "adcp": dict(
        bibkey="Zhang2019",
        author="Zhang, Yuqi and Sanner, Michel F.",
        title="AutoDock CrankPep: combining folding and docking to predict "
              "protein--peptide complexes",
        journal="Bioinformatics", year=2019, volume="35", number="24",
        pages="5121--5127"),
    "adfr": dict(
        bibkey="Ravindranath2015",
        author="Ravindranath, Pradeep Anand and Forli, Stefano and "
               "Goodsell, David S. and Olson, Arthur J. and "
               "Sanner, Michel F.",
        title="AutoDockFR: advances in protein--ligand docking with "
              "explicitly specified binding site flexibility",
        journal="PLoS Computational Biology", year=2015,
        volume="11", number="12", pages="e1004586"),
    "autosite": dict(
        bibkey="Ravindranath2016",
        author="Ravindranath, Pradeep Anand and Sanner, Michel F.",
        title="AutoSite: an automated approach for pseudo-ligands prediction"
              "--from ligand-binding sites identification to predicting key "
              "ligand atoms",
        journal="Bioinformatics", year=2016, volume="32", number="20",
        pages="3142--3149"),
    "pdb2pqr": dict(
        bibkey="Jurrus2018",
        author="Jurrus, Elizabeth and Engel, Dave and Star, Keith and "
               "others",
        title="Improvements to the APBS biomolecular solvation software "
              "suite",
        journal="Protein Science", year=2018, volume="27", number="1",
        pages="112--128"),
    "reduce": dict(
        bibkey="Word1999",
        author="Word, J. Michael and Lovell, Simon C. and Richardson, Jane S. "
               "and Richardson, David C.",
        title="Asparagine and glutamine: using hydrogen atom contacts in the "
              "choice of side-chain amide orientation",
        journal="Journal of Molecular Biology", year=1999, volume="285",
        number="4", pages="1735--1747"),
    "msms": dict(
        bibkey="Sanner1996",
        author="Sanner, Michel F. and Olson, Arthur J. and Spehner, "
               "Jean-Claude",
        title="Reduced surface: an efficient way to compute molecular "
              "surfaces",
        journal="Biopolymers", year=1996, volume="38", number="3",
        pages="305--320"),
    "propka": dict(
        bibkey="Olsson2011",
        author="Olsson, Mats H. M. and S{\\o}ndergaard, Chresten R. and "
               "Rostkowski, Michal and Jensen, Jan H.",
        title="PROPKA3: consistent treatment of internal and surface "
              "residues in empirical pKa predictions",
        journal="Journal of Chemical Theory and Computation", year=2011,
        volume="7", number="2", pages="525--537"),
    "pdbfixer": dict(
        bibkey="Eastman2017",
        author="Eastman, Peter and others",
        title="OpenMM 7: rapid development of high performance algorithms "
              "for molecular dynamics",
        journal="PLoS Computational Biology", year=2017, volume="13",
        number="7", pages="e1005659",
        note="PDBFixer is distributed with OpenMM."),
    "openmm": dict(
        bibkey="Eastman2017",
        author="Eastman, Peter and others",
        title="OpenMM 7: rapid development of high performance algorithms "
              "for molecular dynamics",
        journal="PLoS Computational Biology", year=2017, volume="13",
        number="7", pages="e1005659"),
    "amber": dict(
        bibkey="Case2023",
        author="Case, David A. and others",
        title="AmberTools",
        journal="Journal of Chemical Information and Modeling", year=2023,
        volume="63", number="20", pages="6183--6191"),
    "mmpbsa": dict(
        bibkey="Miller2012",
        author="Miller III, Bill R. and McGee Jr., T. Dwight and "
               "Swails, Jason M. and Homeyer, Nadine and Gohlke, Holger "
               "and Roitberg, Adrian E.",
        title="MMPBSA.py: an efficient program for end-state free energy "
              "calculations",
        journal="Journal of Chemical Theory and Computation", year=2012,
        volume="8", number="9", pages="3314--3321"),
    "autogrid": dict(
        bibkey="Morris2009",
        author="Morris, Garrett M. and Huey, Ruth and Lindstrom, William and "
               "Sanner, Michel F. and Belew, Richard K. and Goodsell, David S. "
               "and Olson, Arthur J.",
        title="AutoDock4 and AutoDockTools4: automated docking with selective "
              "receptor flexibility",
        journal="Journal of Computational Chemistry", year=2009, volume="30",
        number="16", pages="2785--2791"),
    "ff14sb": dict(
        bibkey="Maier2015",
        author="Maier, James A. and Martinez, Carmenza and Kasavajhala, "
               "Koushik and Wickstrom, Lauren and Hauser, Kevin E. and "
               "Simmerling, Carlos",
        title="ff14SB: improving the accuracy of protein side chain and "
              "backbone parameters from ff99SB",
        journal="Journal of Chemical Theory and Computation", year=2015,
        volume="11", number="8", pages="3696--3713"),
    "gbmodel": dict(
        bibkey="Onufriev2004",
        author="Onufriev, Alexey and Bashford, Donald and Case, David A.",
        title="Exploring protein native states and large-scale "
              "conformational changes with a modified generalized born "
              "model",
        journal="Proteins: Structure, Function, and Bioinformatics",
        year=2004, volume="55", number="2", pages="383--394"),
    # RDKit has no paper; its citable record is a Zenodo archive, so it
    # is emitted as @software with a publisher and URL.
    "rdkit": dict(
        bibkey="RDKit", kind="software",
        author="Landrum, Greg and others",
        title="RDKit: Open-source cheminformatics",
        publisher="Zenodo", year=2026, url="https://www.rdkit.org"),
    "molstar": dict(
        bibkey="Sehnal2021",
        author="Sehnal, David and Bittrich, Sebastian and Deshpande, Mandar "
               "and others",
        title="Mol* Viewer: modern web app for 3D visualization and analysis "
              "of large biomolecular structures",
        journal="Nucleic Acids Research", year=2021, volume="49",
        number="W1", pages="W431--W437"),
    "numpy": dict(
        bibkey="Harris2020",
        author="Harris, Charles R. and Millman, K. Jarrod and van der Walt, "
               "St{\\'e}fan J. and others",
        title="Array programming with NumPy",
        journal="Nature", year=2020, volume="585", number="7825",
        pages="357--362"),
    "swisssidechain": dict(
        bibkey="Gfeller2013",
        author="Gfeller, David and Michielin, Olivier and Zoete, Vincent",
        title="SwissSidechain: a molecular and structural database of "
              "non-natural sidechains",
        journal="Nucleic Acids Research", year=2013, volume="41",
        number="D1", pages="D327--D332"),
    "gaff": dict(
        bibkey="Wang2004",
        author="Wang, Junmei and Wolf, Romain M. and Caldwell, James W. and "
               "Kollman, Peter A. and Case, David A.",
        title="Development and testing of a general amber force field",
        journal="Journal of Computational Chemistry", year=2004,
        volume="25", number="9", pages="1157--1174"),
    "ff99sb": dict(
        bibkey="Hornak2006",
        author="Hornak, Viktor and Abel, Robert and Okur, Asim and Strockbine, "
               "Bentley and Roitberg, Adrian and Simmerling, Carlos",
        title="Comparison of multiple Amber force fields and development of "
              "improved protein backbone parameters",
        journal="Proteins: Structure, Function, and Bioinformatics", year=2006,
        volume="65", number="3", pages="712--725"),
    "gbn2": dict(
        bibkey="Nguyen2013",
        author="Nguyen, Hai and Roe, Daniel R. and Simmerling, Carlos",
        title="Improved generalized Born solvent model parameters for protein "
              "simulations",
        journal="Journal of Chemical Theory and Computation", year=2013,
        volume="9", number="4", pages="2020--2034"),
    "3dmol": dict(
        bibkey="Rego2015",
        author="Rego, Nicholas and Koes, David",
        title="3Dmol.js: molecular visualization with WebGL",
        journal="Bioinformatics", year=2015, volume="31", number="8",
        pages="1322--1324"),
}

# The papers every published result must cite, whatever else the run used:
# the docking search, the grid maps it docks into (agfr, from the AutoDockFR
# paper), and the AutoSite points agfr places on every target. Each is named
# here only by its key and by the name it is known by; the reference itself is
# the entry above, so the page's header and the citation list cannot name
# different papers.
MUST_CITE = (("adcp", ENGINE), ("adfr", "AutoDockFR / AGFR"),
             ("autosite", "AutoSite"))


def _authors_short(author):
    """'Zhang and Sanner', or 'Ravindranath et al.', from a BibTeX author
    field."""
    names = [a.strip() for a in author.split(" and ") if a.strip()]
    family = [n.split(",")[0].strip() for n in names]
    if len(family) > 2 or (family and family[-1] == "others"):
        return family[0] + " et al."
    return " and ".join(family)


def must_cite():
    """[(name, short reference, full reference, doi)] for MUST_CITE, read
    from REFERENCES and BIB."""
    full = {k: (ref, doi) for k, _w, _d, ref, doi in REFERENCES}
    out = []
    for key, name in MUST_CITE:
        b = BIB[key]
        ref, doi = full[key]
        out.append((name, "%s, %s" % (_authors_short(b["author"]), b["year"]),
                     ref, doi))
    return out


def must_cite_sentence(html=False):
    """The sentence the page's header and the citation list both carry.

    With `html`, escaped for the page, each program's name in bold with its
    DOI on hover, and each reference kept on one line, so a narrow header
    breaks between papers rather than inside "et al.". The words are the
    same either way.
    """
    import html as _html
    e = _html.escape if html else (lambda s: s)
    parts = []
    for name, short, _ref, doi in must_cite():
        if html:
            parts.append('<b title="doi:%s">%s</b> <span class=ref>(%s)</span>'
                         % (e(doi), e(name), e(short)))
        else:
            parts.append("%s (%s)" % (name, short))
    listed = ", ".join(parts[:-1]) + " and " + parts[-1]
    return ("Results published using this tool must cite the %s papers."
            % listed)


# This package's DOI, once a release is archived. Left blank rather than
# filled with a placeholder.
SELF_DOI = ""


def self_citation(year=2026):
    who = ", ".join(AUTHORS)
    doi = ("  doi:%s" % SELF_DOI) if SELF_DOI else ""
    return ("%s (%s). %s v%s: a browser and command-line workflow for "
            "peptide docking with %s. %s, Texas A&M University.%s"
            % (who, year, APP, __version__, ENGINE, LAB, doi))


def used_by(job):
    """The keys a given run actually used, from its job.json."""
    keys = {k for k, when, *_ in REFERENCES if when == "always"}
    seqs = [(p.get("sequence") or "").upper() for p in job.get("peptides") or []]
    # Which program placed the receptor's hydrogens: what the prepare stage
    # recorded where it ran (reduce replaced by PDB2PQR when reduce could
    # not run), else what the job asks for. A job file written before the
    # choice existed was protonated with PDB2PQR and PROPKA.
    from_adapt = bool((job.get("adapt") or {}).get("folder"))
    done = ((job.get("receptor_info") or {}).get("protonation") or {})
    if from_adapt:
        prot = "reduce"
    else:
        from . import job as _jobmod
        prot = done.get("method") or _jobmod.protonation_of(job)
    by = {"reduce"} if prot == "reduce" else {"pdb2pqr", "propka"}
    # A run given a target file skips this package's preparation, and so
    # does a receptor from ADapt, which reduce protonated there.
    if from_adapt:
        keys |= by
    elif not str(job.get("target_file") or "").strip():
        # PDBFixer checks every receptor for missing atoms, and the
        # minimization is OpenMM with ff14SB and GBn2 (prepare.minimise).
        keys |= by | {"pdbfixer"}
        if job.get("minimise", True):
            keys |= {"openmm", "ff14sb", "gbn2"}
    # MM-GBSA runs only on peptides without non-standard or D residues, and
    # protonates each complex first: with PDB2PQR and PROPKA, or under
    # reduce with the receptor's own states and reduce for the peptide.
    if job.get("mmgbsa") and (not seqs or any(not peptide.is_extended(s) for s in seqs)):
        keys |= {"amber", "ff14sb", "mmpbsa"} | by
        # The paper for the generalised Born model the energies were computed
        # in: OBC(II) by default, GBn2 when the run asked for it.
        keys.add("gbn2" if job.get("mmgbsa_model") == "gbn2" else "gbmodel")
    omm = int(job.get("omm_minimize") or 0)
    if omm:
        # ADCP 1.1 completes each pose with PDBFixer before minimizing it.
        keys |= {"openmm", "pdbfixer", "ff99sb"}
        if (job.get("omm_env") or "implicit") == "implicit":
            keys.add("gbn2")
    if any(peptide.nonstandard((p.get("sequence") or "").upper())
           for p in job.get("peptides") or []):
        keys.add("swisssidechain")
        if int(job.get("omm_minimize") or 0):
            keys.add("gaff")
    return keys


def text(job=None, out=print):
    """The citation block. With a job, only what that run actually used."""
    keys = used_by(job) if job is not None else None
    out("HOW TO CITE\n")
    out("This interface:")
    out("  " + self_citation())
    out("")
    out("%s is a product of the %s and is not ours. It must be cited "
        "whatever else is:" % (ENGINE, ENGINE_OWNER))
    out("  %s" % ENGINE_URL)
    out("")
    # The sentence the page's header carries, built from the same entries.
    out(must_cite_sentence())
    out("")
    if keys is not None:
        out("For this run, from its job.json:")
    else:
        out("Everything this tool can run. Cite what your run used:")
    for key, when, did, ref, doi in REFERENCES:
        if keys is not None and key not in keys:
            continue
        mark = "  " if when == "always" else "  "
        out("%s- %s" % (mark, did))
        out("      %s" % ref)
        out("      doi:%s" % doi)
        if when != "always":
            out("      (only if: %s)" % when)
    out("")
    out("ACKNOWLEDGEMENT")
    out("  " + acknowledgement(job))
    out("")
    out("PI: %s" % PI)


def as_text(job=None):
    """text() as a string, for the page and anything else that serves it."""
    buf = []
    text(job, out=lambda t="": buf.append(t))
    return "\n".join(buf) + "\n"


def formatted(fmt, job=None):
    """One of the three downloadable forms. Raises KeyError on anything else."""
    return {"txt": as_text, "bib": bibtex, "ris": ris}[fmt](job)


FILENAMES = {"txt": "HOW_TO_CITE.txt",
             "bib": "CITATION.bib",
             "ris": "CITATION.ris"}


def verify(out=print, timeout=20):
    """Resolve every DOI against Crossref and print what it resolves to.

    Needs outbound network. Run it after editing REFERENCES: a citation that
    looks right but resolves to a different paper is the failure that survives
    proofreading.
    """
    import json
    import urllib.request
    bad = 0
    out("%-10s %-30s %s" % ("key", "doi", "Crossref"))
    seen = {}
    for key, _when, _did, _ref, doi in REFERENCES:
        # One request per DOI. PDBFixer and OpenMM share a paper, and
        # requesting the same DOI twice in quick succession is
        # rate-limited, which reported a correct DOI as broken.
        if doi in seen:
            out("%-10s %-30s same DOI as %s; %s" % (
                key, doi, seen[doi], "not re-checked"))
            continue
        seen[doi] = key
        if doi in DATACITE:
            try:
                r = urllib.request.urlopen("https://api.datacite.org/dois/"
                                           + doi, timeout=timeout)
                a = json.loads(r.read().decode())["data"]["attributes"]
                out("%-10s %-30s %s (DataCite)" % (key, doi,
                                                   a["titles"][0]["title"][:52]))
            except Exception as exc:                            # noqa: BLE001
                bad += 1
                out("%-10s %-30s DID NOT RESOLVE on DataCite (%s)"
                    % (key, doi, type(exc).__name__))
            continue
        try:
            r = urllib.request.urlopen("https://api.crossref.org/works/" + doi,
                                       timeout=timeout)
            m = json.loads(r.read().decode())["message"]
            title = (m.get("title") or [""])[0]
            import re as _re
            title = _re.sub(r"<[^>]+>", "", title).replace("\n", " ")
            title = " ".join(title.split())
            # The print issue where there is one. Crossref's
            # `issued` is the earliest date it holds, which for
            # several journals is online-first: a year before
            # the issue that is conventionally cited.
            dates = (m.get("published-print") or m.get("issued") or {})
            year = (dates.get("date-parts") or [[None]])[0][0]
            first = (m.get("author") or [{}])[0].get("family", "?")
            out("%-10s %-30s %s, %s (%s)" % (key, doi, first, title[:52], year))
        except Exception as exc:                                # noqa: BLE001
            bad += 1
            out("%-10s %-30s DID NOT RESOLVE (%s)"
                % (key, doi, type(exc).__name__))
    out("")
    out("%d of %d resolved." % (len(REFERENCES) - bad, len(REFERENCES)))
    return 1 if bad else 0


def _bibkey(key, ref):
    """The key a manuscript cites. From BIB where there is one, so that a
    \\cite in the paper and the .bib a user downloads agree."""
    if key in BIB:
        return BIB[key]["bibkey"]
    first = "".join(c for c in ref.split()[0] if c.isalnum()).lower()
    year = "".join(c for c in ref if c.isdigit())[-4:] or "0000"
    return "%s%s_%s" % (first or "ref", year, key)


def bibtex(job=None, year=2026):
    """Every citation for a run as BibTeX, ours included.

    Emitted rather than hand-kept: a .bib that is edited separately from
    REFERENCES drifts from it, and the drift is invisible until someone
    submits.
    """
    keys = used_by(job) if job is not None else None
    out = ["% Generated by `adcp-dock cite --bibtex`. Do not edit by hand;",
           "% regenerate instead, or it will drift from the tool.", ""]
    out.append("@software{%s," % "guitoadcpdocking")
    out.append("  author       = {%s}," % " and ".join(AUTHORS))
    out.append("  title        = {{%s: a browser and command-line "
               "workflow for peptide docking with %s}}," % (APP, ENGINE))
    out.append("  year         = {%d}," % year)
    out.append("  version      = {%s}," % __version__)
    out.append("  organization = {{%s, Texas A&M University}}," % LAB)
    out.append("  url          = {%s}," % REPO)
    out.append("  note         = {No DOI yet; cite the repository until a "
               "release is archived}")
    out.append("}")
    out.append("")
    seen = set()
    for key, when, did, ref, doi in REFERENCES:
        if keys is not None and key not in keys:
            continue
        bk = _bibkey(key, ref)
        # pdbfixer and openmm are the same paper; emitting it twice
        # produces a duplicate key.
        if bk in seen:
            continue
        seen.add(bk)
        b = BIB.get(key)
        if not b:
            # No structured fields for this one. The prose goes
            # in `note` rather than being split up, because
            # guessing at an author/title boundary produces an
            # entry that looks correct and is not.
            out += ["@misc{%s," % bk,
                    "  doi   = {%s}," % doi,
                    "  note  = {%s}," % ref.replace("{", "").replace("}", ""),
                    "  title = {{%s}}" % did, "}", ""]
            continue
        if b.get("kind") == "software":
            out += ["@software{%s," % bk,
                    "  author    = {%s}," % b["author"],
                    "  title     = {{%s}}," % b["title"],
                    "  publisher = {%s}," % b["publisher"],
                    "  year      = {%d}," % b["year"],
                    "  url       = {%s}," % b["url"],
                    "  doi       = {%s}" % doi, "}", ""]
            continue
        out.append("@article{%s," % bk)
        out.append("  author  = {%s}," % b["author"])
        out.append("  title   = {{%s}}," % b["title"])
        out.append("  journal = {%s}," % b["journal"])
        out.append("  year    = {%d}," % b["year"])
        out.append("  volume  = {%s}," % b["volume"])
        if b.get("number"):
            out.append("  number  = {%s}," % b["number"])
        out.append("  pages   = {%s}," % b["pages"])
        if b.get("note"):
            out.append("  note    = {%s}," % b["note"])
        out.append("  doi     = {%s}" % doi)
        out.append("}")
        out.append("")
    out.append("% Acknowledgement to include in the manuscript:")
    for line in _wrap(acknowledgement(job), 74):
        out.append("%% %s" % line)
    return "\n".join(out) + "\n"


def ris(job=None, year=2026):
    """The same list as RIS, which EndNote, Zotero and Mendeley all import."""
    keys = used_by(job) if job is not None else None
    out = []
    out += ["TY  - COMP",
            "TI  - %s: a browser and command-line workflow for peptide "
            "docking with %s" % (APP, ENGINE)]
    for a in AUTHORS:
        bits = a.split()
        out.append("AU  - %s, %s" % (bits[-1], " ".join(bits[:-1])))
    out += ["PY  - %d" % year, "UR  - %s" % REPO,
            "PB  - %s, Texas A&M University" % LAB, "ER  - ", ""]
    seen = set()
    for key, when, did, ref, doi in REFERENCES:
        if keys is not None and key not in keys:
            continue
        b = BIB.get(key)
        if not b:
            out += ["TY  - JOUR", "TI  - %s" % did, "N1  - %s" % ref,
                    "DO  - %s" % doi, "ER  - ", ""]
            continue
        # pdbfixer and openmm are one paper; importing it twice
        # leaves a duplicate for the user to merge.
        if b["bibkey"] in seen:
            continue
        seen.add(b["bibkey"])
        # "--" is a BibTeX en dash and is printed literally by a
        # reference manager.
        if b.get("kind") == "software":
            out += ["TY  - COMP", "TI  - %s" % b["title"]]
            out += ["AU  - %s" % a.strip() for a in b["author"].split(" and ")]
            out += ["PB  - %s" % b["publisher"], "PY  - %d" % b["year"],
                    "UR  - %s" % b["url"], "DO  - %s" % doi,
                    "N1  - used for: %s" % did, "ER  - ", ""]
            continue
        out += ["TY  - JOUR", "TI  - %s" % b["title"].replace("--", "-")]
        # RIS takes one AU per author, "Family, Given", by which
        # the BIB strings are already written.
        for a in b["author"].split(" and "):
            out.append("AU  - %s" % a.strip())
        out += ["JO  - %s" % b["journal"], "PY  - %d" % b["year"],
                "VL  - %s" % b["volume"]]
        if b.get("number"):
            out.append("IS  - %s" % b["number"])
        # SP and EP separately: a reference manager given
        # "5121--5127" in SP prints it as the start page.
        pages = b["pages"].replace("--", "-")
        if "-" in pages:
            sp, ep = pages.split("-", 1)
            out += ["SP  - %s" % sp, "EP  - %s" % ep]
        else:
            out.append("SP  - %s" % pages)
        out += ["DO  - %s" % doi, "N1  - used for: %s" % did, "ER  - ", ""]
    out += ["", "Acknowledgement to include in the manuscript:",
            acknowledgement(job)]
    return "\n".join(out) + "\n"


def _wrap(text, width):
    words, line, out = text.split(), "", []
    for w in words:
        if len(line) + len(w) + 1 > width:
            out.append(line); line = w
        else:
            line = (line + " " + w).strip()
    if line:
        out.append(line)
    return out
