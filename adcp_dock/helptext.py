"""What each setting on the page means, in one place.

By Aadhil Haq

Written for a reader who has not seen the code: each entry says what the
setting does, what changes if it is moved, and where published guidance
exists, what that guidance is.
"""

HELP = {
    # --- receptor and site -------------------------------------------------
    "reccode": (
        "PDB identifier",
        "Four characters, like 7A3R. The structure is fetched from the RCSB "
        "and cached, so each id is downloaded once per machine. Use the file "
        "box instead when the structure is unpublished, or when you have "
        "already prepared it yourself."),
    "recpath": (
        "A structure already on this machine",
        "A .pdb, .cif or .ent file, gzipped or not. Nothing is copied: the "
        "path is read where it is, so it must be somewhere the compute nodes "
        "can also see. A file under your home or scratch directory usually "
        "is; one on your laptop is not. Use \"From this computer\" to upload "
        "a file from your own machine.\n\n"
        "A .pdbqt is not taken: it is a receptor already prepared, and this "
        "tool prepares its own. To dock into one you prepared yourself, make "
        "a target file from it in AGFR or AGFRGUI and give that instead."),
    "targetfile": (
        "A target file from AGFR or AGFRGUI",
        "A .trg file holds a prepared receptor, the grid maps computed from "
        "it and the box they cover. Given one, the run docks into it as it "
        "is: the receptor preparation here does not run, and the box is the "
        "file's own, so no site residues are needed.\n\n"
        "Load reads the file, checks it the way the run will, and shows its "
        "receptor and box in the site view. A file that is truncated, holds "
        "too few maps or was renamed after it was made is refused with the "
        "reason. Do not rename a target: ADCP finds its maps by the file's "
        "name."),
    "adaptframe": (
        "ADapt, on the structure you loaded",
        "ADapt, by Michel Sanner, opens here on every structure you load. It "
        "shows what the file holds: each biological assembly, the hetero "
        "groups, metals and their charges, alternate locations, incomplete "
        "side chains and modified residues. Tick or untick what should stay "
        "in the receptor in its tree, and look at the result in its view.\n\n"
        "Opening it changes nothing: docking uses this page's own preparation "
        "until you change something in ADapt. From the first change the line "
        "above the frame says that docking will use the receptor as set in "
        "ADapt, and the preparation settings are greyed out, since ADapt and "
        "AGFR decide them. You do not need ADapt's save button: when you "
        "submit, ADapt writes the receptor as it is set there, and it is read "
        "back and recorded with the run. Use ADapt's receptor now writes it "
        "at once, so the site can be checked on it.\n\n"
        "Use this page's preparation instead goes back; only a later change "
        "in ADapt brings its receptor back. Hide folds the frame away and "
        "Open full screen shows the same session in a tab of its own. On a "
        "phone a button opens it in a tab instead."),
    "hotspots": (
        "The residues that define the site",
        "CHAIN:NUMBER, like A:98, or bare numbers when the receptor is one "
        "chain. A range works too: A:98-112.\n\n"
        "These do two things. They place the docking box, which is built "
        "around them and then padded; and poses are checked against them "
        "afterwards, so a peptide that lands beside your residues rather than "
        "on them is reported as such.\n\n"
        "Pick the residues you believe a binder must touch, not the whole "
        "pocket. The box grows to cover them, and a box that is too big "
        "spreads the search more thinly everywhere."),
    "site_mode": (
        "Where the box goes",
        "Residues you pick: the box spans them plus the padding. Every run "
        "made before this choice existed placed its box this way.\n\n"
        "A ligand or peptide in the structure: the box spans it plus the "
        "padding, as AGFR places a box around a ligand, and the ligand is "
        "kept out of the receptor. The residues it touches, within 4 \u00c5, "
        "become the site every pose is measured against.\n\n"
        "A pocket AutoSite finds: once the receptor is prepared, on the "
        "compute node, AutoSite searches the whole receptor and the box goes "
        "around its best scoring pocket, or around all of them. The residues "
        "lining that pocket become the site. Nothing can be drawn before the "
        "run, because the pocket is found on the prepared receptor; load the "
        "run's target file afterwards to see it.\n\n"
        "However it is placed, the box is held between the minimum and the "
        "maximum."),
    "site_ligand": (
        "The ligand or peptide the box is built around",
        "What the loaded structure holds that a box could be built around: "
        "peptides, the protein chains of up to 50 residues, named by their "
        "chain; and hetero groups, named CHAIN:NAME:NUMBER, like A:NUT:201. "
        "Crystallisation additives come last and waters are not listed.\n\n"
        "The ligand is never docked into. A peptide's chain is left out of "
        "the receptor, and if you tick it among the chains to keep the run "
        "is refused. A hetero group goes with the other heteroatoms; with "
        "Keep heteroatoms on it would stay and fill its own site, so that is "
        "refused too.\n\n"
        "AGFR 1.1 can pick a ligand out of a PDB itself (--pdbLigand), but "
        "only on a receptor it prepares itself, and it keeps a peptide chosen "
        "that way in the receptor. So the box is worked out here, by AGFR's "
        "own rule, and given to AGFR as a box."),
    "site_pocket": (
        "Which pocket",
        "The best scoring pocket is AutoSite's first, ranked by how many "
        "points it has, how strongly they bind and how buried they are "
        "(agfr -p best).\n\n"
        "Every pocket puts one box around all the pockets AutoSite finds "
        "(agfr -p all). On a whole protein that is close to docking blind, "
        "and the maximum box still applies.\n\n"
        "agfr builds the maps for a pocket at 1 \u00c5 spacing, too coarse to dock "
        "into, so its box is kept and the maps are built over it again at "
        "0.375 \u00c5, as for any other box."),
    "agfr": (
        "Which AGFR builds the maps",
        "AGFR builds the grid maps over the box that the peptide is docked "
        "into. AGFR 1.1 comes with the ADCP 1.1 suite and is used wherever "
        "that suite is installed; AGFR 1.0 comes with ADFRsuite.\n\n"
        "Over the same receptor and box the two give the same box, the same "
        "receptor and the same starting points. Their maps differ slightly, "
        "because each runs its own build of AutoGrid: by at most 0.15 "
        "kcal/mol where a peptide atom can sit, and by up to 1.2 in the "
        "hydrogen bond donor map. ADCP 1.0 and ADCP 1.1 dock into maps from "
        "either."),
    "autosite": (
        "Which AutoSite",
        "AutoSite marks the places in the box where a ligand is likely to "
        "bind. ADCP's search starts from them, and for a box from a pocket "
        "they are the pocket.\n\n"
        "AutoSite 1.0 keeps only the grid points that bind most strongly, so "
        "it tends to miss the weakly bound parts of a large ligand such as a "
        "peptide. 1.1 relaxes its cutoffs and scans several of them, which is "
        "why the ADCP authors recommend it for peptides.\n\n"
        "It changes where the search starts, not the maps. On the MDM2 "
        "example the two placed 246 and 380 points, and on a short search "
        "that can change which pose comes out best."),
    "receptor_gradient": (
        "The receptor gradient",
        "Inside the receptor no atom can be, and the maps hold very high "
        "values there. The gradient replaces them with values that slope "
        "toward the receptor's surface, so a peptide atom caught inside is "
        "pushed out (AGFR's -g). Where a peptide atom can actually sit the "
        "maps are unchanged.\n\n"
        "It is off unless you tick it: the ADCP authors found it expensive "
        "and of no real help to the search. ADFRsuite's AGFR 1.0 adds it "
        "unless told not to, and on the MDM2 example it took 2.5 minutes "
        "there against 14 seconds without; AGFR 1.1 took 23 seconds more."),
    "pad": (
        "Padding around what marks the site",
        "How far the box extends beyond the residues, the ligand or the "
        "pocket, in angstroms (\u00c5), on each side. 8 \u00c5 leaves room for a peptide to "
        "lie across the site. With less than about 5 \u00c5, a longer peptide has "
        "nowhere to put its ends."),
    "bmin": (
        "Smallest the box may be",
        "A minimum, in Angstrom. Two or three hotspot residues close together "
        "would otherwise give a box barely larger than they are, which a "
        "peptide cannot enter from any direction."),
    "bmax": (
        "Largest the box may be",
        "A maximum, in Angstrom. The search is spread over the whole box, so "
        "a wider box does not cost more time; it spreads the search more "
        "thinly. If your residues span more than this, they are probably two "
        "sites and are better docked separately."),

    # --- peptides ----------------------------------------------------------
    "pastetext": (
        "The peptides to dock",
        "One per line. Three forms are read:\n\n"
        "    ETFSDLWKLLP            named by its sequence\n"
        "    p53: ETFSDLWKLLP       named p53\n"
        "    >p53                   FASTA, named p53\n"
        "    ETFSDLWKLLP\n\n"
        "Lines starting with # are ignored.\n\n"
        "Name a peptide with a colon, not a space. Spaces inside a sequence "
        "are joined, because sequences are often printed in blocks of ten, "
        "so 'alpha SWQWWWQPQRDCR' would otherwise be read as one peptide, "
        "ALPHASWQWWWQPQRDCR. The page refuses that and says why.\n\n"
        "Standard amino acids only, in one-letter code, 5 to 30 residues. "
        "On a cluster each peptide becomes its own Slurm task, so twenty "
        "peptides are twenty tasks that queue and run independently, and "
        "one failing does not lose the others."),

    # --- docking -----------------------------------------------------------
    "poses": (
        "Poses kept per peptide",
        "How many of ADCP's ranked poses to write and keep (ADCP's -m). They "
        "are what you look at afterwards and what MM-GBSA rescores, so this "
        "sets the size of the answer rather than the size of the search.\n\n"
        "Ten is the default and enough to see whether the peptide's poses "
        "agree on where it binds. ADCP 1.0 writes at most ten; ADCP 1.1 "
        "writes up to 100, one for each cluster its search found, best "
        "first. More poses cost a few seconds to write, not a longer "
        "search."),
    "cyclic": (
        "A chain or a ring",
        "A linear chain is how every run before this setting was docked.\n\n"
        "A ring head to tail joins the first residue's nitrogen to the last "
        "residue's carbon (ADCP's -cyc). Each replica's search then takes "
        "about one and a half times as long.\n\n"
        "A ring through a disulfide joins two cysteines (ADCP's -cys), so "
        "every peptide needs at least two C, and one with fewer is refused "
        "before anything is submitted: ADCP itself would dock it as a linear "
        "chain and say nothing.\n\n"
        "Both are restraints ADCP's search pulls on, not bonds it starts "
        "with, so check that the ends, or the two sulfurs, meet in the poses "
        "you use. A ring is not rescored with MM-GBSA here, which builds a "
        "linear peptide. ADCP 1.1's OpenMM minimization can be used "
        "instead: it adds the ring's closing bond, and for a disulfide it "
        "does so only where the two sulfurs are within 3.5 \u00c5 of each "
        "other."),
    "helical": (
        "Replicas that start helical",
        "Each replica of ADCP's search starts from a shape: a helix or a "
        "coil. This is the share that start as a helix (ADCP's -p). ADCP "
        "gives replicas 1 to N x p / 100, rounded down, a helical start and "
        "the rest a coil: at 25% with 50 replicas, replicas 1 to 12.\n\n"
        "ADCP decides this from the case of the sequence's letters, upper "
        "case for helix and lower case for coil, and -p rewrites them replica "
        "by replica, so how a sequence is typed makes no difference. The "
        "library names of non-standard residues, inside angle brackets, are "
        "passed as written.\n\n"
        "At 0 every replica starts as a coil, which is how every run before "
        "this setting was docked; at 100 every one starts as a helix. ADCP's "
        "-p does nothing at 0 or 100, so those two are given to ADCP as an all "
        "lower case or all upper case sequence instead.\n\n"
        "The search is as long either way; only where it starts changes."),
    "keep_trajectories": (
        "Keep every replica's trajectory",
        "Each replica of ADCP's search writes its trajectory: every pose it "
        "accepted on the way. ADCP 1.1 deletes them when the docking ends "
        "unless told to keep them (its -k); ADCP 1.0 always leaves them, and "
        "they are deleted here.\n\n"
        "Ticked, they are packed into one archive per peptide, beside its "
        "poses: docking/<peptide>/<peptide>_trajectories.tar.gz. About 150 "
        "files per peptide become one. Two copies are left out of it: the "
        "maps ADCP unpacked from the target, and out.pdb, which is every "
        "trajectory again in one file. ADCP makes both again when the "
        "archive is unpacked in the peptide's folder and asked to recluster "
        "(-r).\n\n"
        "Packing adds a few seconds per peptide, and the archive is tens of "
        "megabytes, so leave this off unless you will look at the "
        "trajectories."),
    "replicas": (
        "Independent searches",
        "ADCP runs this many searches from different random starts and keeps "
        "the best. They are the main protection against one search settling "
        "in the wrong place.\n\n"
        "The authors' guideline is 100 to 300:\n\n"
        "    \"A guideline should be running 100-300 replicas with a maximum "
        "Monte Carlo steps of 1-3 million per amino acid in the peptide.\"\n"
        "    (ccsb.scripps.edu/adcp/documentation)\n\n"
        "The v1.1 paper used 300, the original paper 80, and the redocking "
        "tutorial suggests 100. ADCP's own default is 50, below all of "
        "them.\n\n"
        "The cost is not linear. Replicas run across the cores the task has, "
        "so up to the core count extra replicas add little time: 100 "
        "replicas on 48 cores run in three waves, the same as 144 would. "
        "Fewer replicas than cores leaves cores idle."),

    "steps": (
        "Search effort per replica",
        "How long each search runs. ADCP's guidance is that this scales "
        "with the length of the peptide:\n\n"
        "    \"A guideline should be running 100-300 replicas with a maximum "
        "Monte Carlo steps of 1-3 million per amino acid in the peptide "
        "(10-30 million steps for a 10-mer peptide).\"\n"
        "    (ccsb.scripps.edu/adcp/documentation)\n\n"
        "So pick roughly your peptide's length in residues times 1 to 3 "
        "million:\n\n"
        "    5 residues     5,000,000 to 15,000,000\n"
        "    8 residues     8,000,000 to 24,000,000\n"
        "    10 residues   10,000,000 to 30,000,000\n"
        "    15 residues   15,000,000 to 45,000,000\n\n"
        "The redocking tutorial says the same: \"we suggest running 100 "
        "replicas with 1 million steps per amino acid\". The v1.1 paper used "
        "300 replicas at 1 million per residue; the original paper used 80 "
        "at 3 million per residue.\n\n"
        "It is a maximum, not a fixed cost: a search that converges early "
        "stops. Too few steps show up as poses that disagree with each other "
        "and a best score that changes when the run is repeated."),

    "protonation": (
        "How the receptor's hydrogens are placed",
        "reduce, the default, is the program ADFR's own preparation and "
        "ADapt use. It adds every hydrogen, turns Asn and Gln amides and "
        "His rings over where the hydrogen bonds say the file has them the "
        "wrong way round, and puts each histidine's protons where the "
        "hydrogen bonds around it call for them. Checked against protein "
        "structures solved by neutron diffraction, where the hydrogens are "
        "seen, it placed histidine protons as deposited in 55% of 1,457 "
        "cases. It takes no pH: aspartates and glutamates stay charged, "
        "lysines and arginines positive.\n\n"
        "PROPKA at a pH predicts each titratable residue's pKa and sets its "
        "state at the pH you give, with PDB2PQR placing the hydrogens. Use "
        "it when the pH matters to your question, an endosomal target at "
        "5.5 for example. On the same neutron structures it placed 31% of "
        "the histidines as deposited.\n\n"
        "MM-GBSA keeps the receptor's choice: under reduce its complexes are "
        "named from the receptor's own histidines, and under PROPKA they are "
        "protonated again at the same pH. Runs made before this choice "
        "existed used PDB2PQR with PROPKA, and reusing their settings keeps "
        "it."),
    "ph": (
        "pH the receptor is protonated at",
        "Used only when the receptor is protonated with PROPKA. PROPKA "
        "decides each titratable residue's state at this pH, which sets the "
        "charges the docking sees: whether a histidine is neutral or "
        "positive, whether an aspartate is charged.\n\n"
        "7.4 is physiological. If your assay runs at a different pH, use "
        "that: an endosomal target at 5.5 has a different surface. The same "
        "pH is used again when MM-GBSA builds its complexes. reduce takes "
        "no pH, so with reduce chosen this box is greyed out."),
    "minimise": (
        "Minimise after protonation",
        "A restrained minimisation: heavy atoms are held near their measured "
        "positions while the newly placed hydrogens and any clashes "
        "relax.\n\n"
        "Leave it on. Hydrogens are placed by geometry and have not been "
        "relaxed against the force field, so docking without this docks "
        "into strain that is not really there. It takes under a minute for "
        "a few hundred residues."),
    "keep_hetatms": (
        "Keep heteroatoms",
        "Cofactors, ligands, metals and other non-protein atoms are removed "
        "unless you tick this.\n\n"
        "Tick it when what is bound there is part of the site, such as a "
        "structural zinc or a haem your peptide must sit against. Leave it "
        "off for a crystallisation additive or a ligand you intend to "
        "displace. The minimisation's force field has no parameters for "
        "most of what this keeps, metal ions included, so it is then "
        "skipped, and the results page says so."),
    "rebuild_sidechains": (
        "Rebuild missing side chains",
        "Many structures lack the ends of some side chains, usually on the "
        "surface where they were never resolved. On, which is the default, "
        "PDBFixer rebuilds the missing atoms and they are docked against "
        "like the rest.\n\n"
        "Off, the receptor that is docked into keeps its side chains as the "
        "structure has them. With reduce and no minimisation the missing "
        "atoms are never built: reduce places hydrogens on a residue as it "
        "stands. The minimisation needs whole residues, and PDB2PQR "
        "completes every residue it is given, so with either of those the "
        "missing atoms are modelled for a moment and then taken out with "
        "the hydrogens on them. The backbone, the terminal oxygen and any "
        "loop you ask to be rebuilt are completed either way.\n\n"
        "A rebuilt side chain is placed from a template, not measured. Turn "
        "this off when one would sit in your site and you would rather dock "
        "against the gap. MM-GBSA still needs whole residues, so it "
        "completes them in its own complex."),
    "fill_loops": (
        "Rebuild missing loops",
        "Whole missing residues, the loops the structure never resolved. "
        "Missing side chain atoms are a separate setting, just before "
        "this one.\n\n"
        "A rebuilt loop is a model, and once written it looks like a "
        "measured one. Turn this on only when a gap sits against your site "
        "and leaving a hole there would be worse than filling it. A rebuilt "
        "loop is built whole, side chains included, whichever way the side "
        "chain setting is."),
    "mmgbsa": (
        "Rescore with MM-GBSA",
        "A physics-based rescoring of the docked poses: the complex is "
        "protonated with the peptide bound, minimised, and its interaction "
        "energy computed.\n\n"
        "It separates peptides that the docking score ranks together, and it "
        "is the slowest step in the run, taking roughly as long again as the "
        "docking. It gives an ordering within one run against one receptor, "
        "not a binding constant. It needs AmberTools."),
    "mmposes": (
        "Poses to rescore",
        "MM-GBSA runs on this many of each peptide's poses. The cost grows "
        "in proportion, and it is the largest cost of the run.\n\n"
        "More than one is worth having: a single pose can be lucky, and the "
        "spread across poses shows whether the peptide has settled on one "
        "way of binding."),
    "topn": (
        "Average the best few",
        "The reported MM-GBSA value is the mean of this many best pose "
        "energies rather than the single best, which makes it less sensitive "
        "to one outlier. Set it to 1 to report the best pose alone."),

    # --- resources ---------------------------------------------------------
    "cpus": (
        "Cores per task",
        "Each peptide runs as one task with this many cores. ADCP spreads its "
        "replicas across them, so more cores are faster until there are more "
        "cores than replicas.\n\n"
        "Asking for a whole node usually means a longer wait in the queue; "
        "half a node often starts and finishes sooner."),
    "walltime": (
        "Time to book",
        "The scheduler stops the job at this limit, finished or not, so it "
        "is better to book generously: clusters usually charge for the time "
        "used, not the time booked. The estimate beside this box already "
        "includes a margin, and the line under it adds up what every job of "
        "the run books: the prepare job, the docking tasks and the report."),
    "throttle": (
        "Most tasks at once",
        "Limits how many peptides run at the same time, however many are "
        "queued. Use it on a busy cluster, or to keep a large screen from "
        "taking the whole allocation. Empty means no limit."),
    "label": (
        "What to call this run",
        "It names the run directory and the Slurm jobs, so it is what you "
        "see in the queue when several runs are going. Use something that "
        "tells runs apart, such as the receptor and what you are trying, "
        "rather than 'test2'."),
    "outdir": (
        "Where the run is written",
        "Everything the run produces goes here: the prepared receptor, the "
        "poses, the scores and the report. On a cluster it must be visible "
        "to the compute nodes, which a scratch directory is. Runs are never "
        "overwritten: the run directory is the name you give followed by the "
        "date and time, so a second run with the same name sits beside the "
        "first rather than on top of it."),
    "sitesurf": (
        "Show the receptor's surface",
        "Draws the molecular surface of the chains you keep over the "
        "structure, with the residues you picked coloured blue on it, so a "
        "groove shows as a groove and the patch the box is built around "
        "shows with it.\n\n"
        "The surface is computed by MSMS, Michel Sanner's program, which "
        "ships with ADFRsuite: the surface a probe of 1.5 \u00c5, a water "
        "molecule, traces as it rolls over the atoms. Where MSMS is not "
        "installed or cannot triangulate the structure, the viewer draws "
        "its own surface instead and says so under the view.\n\n"
        "This changes only what you see; it has no effect on the docking."),
    "showbox": (
        "Show the docking box",
        "Draws the box in the structure above, as a cage around the residues "
        "you picked.\n\n"
        "ADCP docks into what is inside this box and nothing else, so it is "
        "worth checking before you submit: whether it covers the groove, "
        "whether it reaches a neighbouring chain, and how much of it is "
        "solvent. If it misses part of the site, add a residue there or "
        "raise the padding; if it covers far too much, lower the padding or "
        "the maximum. It is drawn from the same numbers the run will use, "
        "and it updates as you change them.\n\n"
        "This changes only what you see; it has no effect on the docking."),
    "account": (
        "The allocation to charge",
        "Leave it empty on clusters that charge your default account.\n\n"
        "Many clusters refuse a job that does not name one, which the page "
        "reports as 'the account on the job is not one you can use'. Enter "
        "the account name exactly as your cluster gives it (on Slurm, "
        "`sacctmgr show assoc user=$USER` lists them). It is added to every "
        "job this run submits."),
    "sbatch_extra": (
        "Anything else this cluster requires",
        "Slurm options your site needs that this form does not ask about, "
        "each as --name=value, separated by spaces:\n\n"
        "    --qos=normal\n"
        "    --partition=short\n"
        "    --constraint=icelake\n\n"
        "They are added to every job this run submits. The options the tool "
        "sets itself (the name, the log file, cores, memory and time) cannot "
        "be overridden here, because the page reads them back. Nothing that "
        "could end the line or run a command is accepted."),
    "box_chains": (
        "Other chains the box reaches",
        "The box built around your residues also covers part of another "
        "chain.\n\n"
        "Tick this when that chain is the far wall of the same groove: "
        "without it the peptide is docked against a gap where the wall "
        "should be. Leave it unticked when that chain is the partner your "
        "peptide should displace, since docking against it would dock "
        "against what you want removed.\n\n"
        "Ticking loads the whole chain, not only the part inside the box."),
}

# Search effort, as five numeric choices: the number is what ADCP is given
# and what tutorials and papers quote, whereas a label such as "thorough"
# cannot be compared with anything published.
STEP_CHOICES = [
    (100_000, "a smoke test only, far below any published setting"),
    (500_000, "a rough look; under one residue's worth"),
    (2_500_000, "ADCP's default: about 2 residues at 1 million each"),
    (10_000_000, "a 10-mer at the guideline's lower bound"),
    (20_000_000, "a 10-mer at 2 million per residue, or a 20-mer at 1 million"),
]

# Independent searches. ADCP runs these across the cores the task holds, so
# below the core count they add little wall time, and they are the main
# defence against one search settling in the wrong place.
REPLICA_CHOICES = [
    (25, "below every published setting; a quick look only"),
    (50, "ADCP's default, below the 100 to 300 guideline"),
    (100, "the guideline's lower bound, and the tutorial's suggestion"),
    (200, "within the guideline; for a difficult site"),
    (400, "above the v1.1 paper's 300"),
]
