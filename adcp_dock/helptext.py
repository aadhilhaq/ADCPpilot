"""What each setting on the page means, in one place.

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
        "a file from your own machine."),
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
    "pad": (
        "Padding around the residues",
        "How far the box extends beyond the hotspot residues, in Angstrom. "
        "8 A leaves room for a peptide to lie across the site. With less "
        "than about 5 A, a longer peptide has nowhere to put its ends."),
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
        "How many ranked poses ADCP writes out. They are what you look at "
        "afterwards and what MM-GBSA rescores, so this sets the size of the "
        "answer rather than the size of the search.\n\n"
        "Ten is enough to see whether the peptide's poses agree on where it "
        "binds. Keeping more costs disk space, not time."),
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

    "ph": (
        "pH the receptor is prepared at",
        "PROPKA decides each titratable residue's state at this pH, which "
        "sets the charges the docking sees: whether a histidine is neutral or "
        "positive, whether an aspartate is charged.\n\n"
        "7.4 is physiological. If your assay runs at a different pH, use "
        "that: an endosomal target at 5.5 has a different surface."),
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
        "displace. An unusual cofactor may not have force-field parameters, "
        "so keeping one can make the preparation fail."),
    "fill_loops": (
        "Rebuild missing loops",
        "Missing side chains are always rebuilt; this setting is about whole "
        "missing residues, the loops the structure never resolved.\n\n"
        "A rebuilt loop is a model, and once written it looks like a "
        "measured one. Turn this on only when a gap sits against your site "
        "and leaving a hole there would be worse than filling it."),
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
        "includes a margin."),
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
        "overwritten; a name already in use is refused."),
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
