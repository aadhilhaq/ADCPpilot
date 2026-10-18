"""What ADCP accepts as a peptide, and why a given string does not qualify.

By Aadhil Haq

Problems are returned as structured objects (which characters, at which
positions, and what is wrong with the length) so a form can highlight the
offending letters rather than printing "invalid sequence".

ALPHABET. ADCP builds peptides from a rotamer library keyed on the twenty
standard amino acids. B, J, O, U, X and Z are excluded: B and Z are the
ambiguity codes (Asx, Glx), J is Leu-or-Ile, X is unspecified, O is
pyrrolysine and U is selenocysteine. None has a rotamer ADCP can build.
"""
import re

STANDARD = "ACDEFGHIKLMNPQRSTVWY"
#: Letters that are letters but not residues ADCP can build.
NOT_RESIDUES = "BJOUXZ"
# ADCP cannot build a chain shorter than five residues: with three or four,
# every replica fails. Measured with the same residues at each length.
MIN_LENGTH = 5
# The upper bound is set by the docking box rather than by the search: a
# 30-mer needs about 48 A of cube (sizing.box_for_length) against a maximum
# box of 60 A. Above roughly 20 residues the cost estimate, not this limit,
# is the number to read before submitting.
MAX_LENGTH = 30

_WS = re.compile(r"\s+")

# ADCP 1.1 also takes non-standard and D amino acids: an optional "&" for the
# D form, the one-letter code of the standard residue whose coarse potential
# is used, and the non-standard side chain's library name in angle brackets.
# "N<BHD>PIS&DVD" is seven residues, the first an asparagine derivative and
# the fifth a D-aspartate.
EXTENDED = re.compile(r"&?[A-Z](?:<[A-Z0-9]{1,4}>)?")


def is_extended(seq):
    """Whether a sequence uses the ADCP 1.1 syntax."""
    return "<" in seq or "&" in seq


def residues(seq):
    """The residues of a cleaned sequence, one string each ("N<BHD>", "&D")."""
    return EXTENDED.findall(seq) if is_extended(seq) else list(seq)


def nonstandard(seq):
    """The side-chain library names used in a sequence."""
    return re.findall(r"<([A-Z0-9]{1,4})>", seq)


class Problem:
    """One thing wrong with a sequence.

    `positions` are 0-based indices into the CLEANED sequence, which is what
    a form should be highlighting: the person did not type the whitespace
    that was stripped out.
    """

    def __init__(self, kind, message, positions=()):
        self.kind = kind
        self.message = message
        self.positions = list(positions)

    def __repr__(self):                                    # pragma: no cover
        return "Problem(%r, %r, %r)" % (self.kind, self.message, self.positions)

    def as_dict(self):
        return {"kind": self.kind, "message": self.message,
                "positions": self.positions}


# A peptide is docked as a linear chain unless the run says otherwise. ADCP
# can also close it into a ring: through the backbone, head to tail (-cyc),
# or through a disulfide between two of its cysteines (-cys). The run makes
# one choice, which every peptide takes unless the peptide carries its own
# (job.ring_of), so one run can dock some peptides linear and others as
# rings.
CYCLIC = ("none", "backbone", "disulfide")
#: What each choice is called in the logs and the batch scripts.
CYCLIC_LABELS = {"none": "linear", "backbone": "backbone cyclic, head to tail",
                 "disulfide": "cyclic through a disulfide"}
#: What each is called beside a peptide: on its row in the page, in the
#: results table, and in the ring column of a peptide file.
RING_NAMES = {"none": "linear", "backbone": "head to tail",
              "disulfide": "disulfide"}
#: The same in a sentence: "docked as a ring head to tail".
RING_AS = {"none": "as a linear chain", "backbone": "as a ring head to tail",
           "disulfide": "as a ring through a disulfide"}

# What people write for each choice, in a peptide file's ring column or a
# form. Read with case, dashes and underscores ignored, so "Head-to-tail"
# and "head_to_tail" are "head to tail". ADCP's own flags are here too, as
# cyc and cys, since someone used to ADCP will write those. An empty cell,
# "as the run" or "default" leaves the peptide to the run's own choice.
_RING_WORDS = {
    "linear": "none", "none": "none",
    "head to tail": "backbone", "backbone": "backbone", "cyclic": "backbone",
    "backbone cyclic": "backbone", "cyc": "backbone",
    "disulfide": "disulfide", "disulphide": "disulfide", "cys": "disulfide",
    "": "", "as the run": "", "default": "",
}


def _ring_text(text):
    return " ".join(str(text if text is not None else "").lower()
                    .replace("-", " ").replace("_", " ").split())


def ring_choice(text):
    """The ring a cell or a form names: "none", "backbone" or "disulfide",
    "" for a blank (the run's own choice), and None for anything else."""
    return _RING_WORDS.get(_ring_text(text))


def ring_problem(text):
    """Why `text` is not a ring a peptide can be given, as a Problem."""
    return Problem(
        "ring",
        "The ring is given as '%s', which is not one of the choices: linear "
        "(or none), head to tail (or backbone, or cyclic), or disulfide. "
        "Leave it empty to dock this peptide as the run is set."
        % str(text if text is not None else "").strip())


def cysteines(seq):
    """How many residues of a cleaned sequence a disulfide can join: C, and
    its D form &C.

    A non-standard residue built on the cysteine coarse model, such as
    C<HCS>, is not counted. Its side chain is the library's, and whether it
    has a thiol where cysteine's is depends on the residue.
    """
    return sum(1 for r in residues(seq) if r in ("C", "&C"))


def cyclic_problem(seq, mode):
    """A Problem when `mode` cannot close this peptide, else None.

    ADCP does not check this itself. Given -cys and a peptide with one
    cysteine or none, ADCP 1.1.21 docked it and reported poses with no
    warning, so a run would describe as cyclic a peptide that was
    never closed. It is refused here instead, before anything is submitted.
    """
    mode = str(mode or "none")
    if mode not in CYCLIC:
        return Problem("cyclic", "The cyclic choice must be one of %s, not %r."
                       % (", ".join(CYCLIC), mode))
    if mode == "disulfide":
        n = cysteines(seq)
        if n < 2:
            return Problem(
                "cyclic",
                "A disulfide ring needs two cysteines and this peptide has %s. "
                "Add the two cysteines the ring closes through, or dock it as "
                "a linear chain or a ring head to tail."
                % ("none" if n == 0 else "one"))
    return None


def clean(text):
    """Upper-cased, with whitespace and a FASTA header removed.

    People paste FASTA, and they paste sequences wrapped over several lines
    with numbers down the side. Dropping the header and the whitespace is not
    being clever: refusing a pasted FASTA record would only teach them to
    delete the first line by hand every time.
    """
    if not text:
        return ""
    lines = [l for l in str(text).splitlines() if not l.startswith(">")]
    return _WS.sub("", "".join(lines)).upper()


def check(text, min_length=MIN_LENGTH, max_length=MAX_LENGTH, library=None):
    """(cleaned, [Problem]). An empty problem list means ADCP will take it.

    ADCP 1.1's syntax for non-standard and D residues is accepted, and each
    side-chain name is checked against `library` when one is given.
    """
    seq = clean(text)
    out = []
    if not seq:
        return seq, [Problem("empty", "Enter a peptide sequence.")]
    if is_extended(seq):
        return seq, _check_extended(seq, min_length, max_length, library)

    bad_letters, not_letters = [], []
    for i, c in enumerate(seq):
        if c in STANDARD:
            continue
        (bad_letters if c in NOT_RESIDUES else not_letters).append(i)

    if not_letters:
        out.append(Problem(
            "not_a_letter",
            "The sequence must be one-letter amino acid codes. "
            + _names(seq, not_letters) + " is not one.",
            not_letters))
    if bad_letters:
        out.append(Problem(
            "ambiguous",
            "The sequence must consist only of letters except "
            + ", ".join(NOT_RESIDUES) + ". "
            + _names(seq, bad_letters)
            + " cannot be built as a residue.",
            bad_letters))
    if len(seq) < min_length:
        out.append(Problem(
            "too_short",
            "ADCP needs at least %d residues; this has %d."
            % (min_length, len(seq))))
    if len(seq) > max_length:
        out.append(Problem(
            "too_long",
            "This has %d residues and the limit is %d. Longer peptides can be "
            "docked by raising the limit, but the search cost climbs steeply "
            "with length." % (len(seq), max_length),
            list(range(max_length, len(seq)))))
    return seq, out


def _check_extended(seq, min_length, max_length, library):
    out = []
    # Whatever the pattern does not account for is a mistake in the syntax.
    covered = set()
    for m in EXTENDED.finditer(seq):
        covered.update(range(m.start(), m.end()))
    stray = [i for i in range(len(seq)) if i not in covered]
    if stray:
        out.append(Problem(
            "syntax",
            "Write a non-standard residue as a one-letter code followed by its "
            "library name in angle brackets, such as N<BHD>, and a D residue "
            "with & before its letter, such as &D. " + _names(seq, stray)
            + " does not fit.", stray))
    for m in EXTENDED.finditer(seq):
        letter = m.group(0).lstrip("&")[0]
        if letter not in STANDARD:
            i = m.start() + (1 if m.group(0).startswith("&") else 0)
            out.append(Problem("ambiguous", "%s is not one of the twenty "
                               "standard residues ADCP uses for the coarse "
                               "model." % letter, [i]))
    if library is not None:
        unknown = sorted({n for n in nonstandard(seq) if n not in library})
        if unknown:
            out.append(Problem(
                "unknown_nst",
                "%s %s not in ADCP's side-chain library."
                % (", ".join(unknown), "is" if len(unknown) == 1 else "are")))
    n = len(residues(seq))
    if n < min_length:
        out.append(Problem("too_short", "ADCP needs at least %d residues; "
                           "this has %d." % (min_length, n)))
    if n > max_length:
        out.append(Problem("too_long", "This has %d residues and the limit "
                           "is %d." % (n, max_length)))
    return out


def _names(seq, idx):
    """'K at position 4' / 'X at positions 2, 9': 1-based, as a person counts."""
    chars = sorted({seq[i] for i in idx})
    where = ", ".join(str(i + 1) for i in idx[:8]) + ("..." if len(idx) > 8 else "")
    return "%s at position%s %s" % ("/".join(chars),
                                    "" if len(idx) == 1 else "s", where)


def ok(text, **kw):
    seq, problems = check(text, **kw)
    return seq if not problems else None


def looks_like_a_name(line):
    """Is this line a NAME followed by a sequence, rather than a sequence?

    Spaces within a line are joined deliberately, because sequences are often
    printed in upper-case blocks of ten: "ETFSDLWKLL PENNVLSPLP" is one
    peptide. The consequence is that "alpha SWQWWWQPQRDCR" would also be read
    as one peptide, ALPHASWQWWWQPQRDCR, since A, L, P, H and A are all
    residues: silently docking a different, longer molecule.

    Case distinguishes the two: sequences are written in capitals, names
    usually are not. A first word containing a lower-case letter, followed by
    words that are entirely capitals, is treated as a name and refused with
    the FASTA or colon form to use instead. An all-capital name cannot be
    distinguished from a blocked sequence and is still joined.
    """
    words = str(line or "").split()
    if len(words) < 2:
        return False
    head, rest = words[0], words[1:]
    return (any(c.islower() for c in head)
            and all(w.isupper() and w.isalpha() for w in rest))


# Column names a header row might use, so a CSV that says which way round it
# is does not have to be guessed at. The ring column is optional.
_NAME_WORDS = ("name", "label", "id", "peptide_name", "construct")
_SEQ_WORDS = ("sequence", "seq", "peptide", "aa", "residues")
_RING_COLUMN = ("ring", "cyclic", "cyclisation", "cyclization", "closure",
                "topology")


def _reads_as_sequence(text, **kw):
    """True when this field is a usable sequence on its own.

    Used to work out which field of a row is the sequence, so both
    `SEQUENCE,name` and `name,SEQUENCE` work without being told which.
    """
    t = (text or "").strip()
    if not t:
        return False
    seq, probs = check(t, **kw)
    return bool(seq) and not probs


def _fields(line):
    """(separator, fields) for a line split on a tab, a comma or a
    semicolon, the first of those the line holds; (None, None) for a line
    with none. Empty fields are kept, so the column a header row names is
    still the right one in a row that leaves a cell empty."""
    for sep in ("\t", ",", ";"):
        if sep in line:
            return sep, [x.strip() for x in line.split(sep)]
    return None, None


def _split_row(line):
    """The filled fields of a line split on a tab, a comma or a semicolon,
    when there are two or three of them.

    A colon is handled separately and earlier, because `name: SEQUENCE` is
    the documented one line form. Two fields are a sequence and its name, or
    a sequence and its ring; three are a sequence, a name and a ring. Four or
    more is a spreadsheet row this function should not guess at.
    """
    _sep, parts = _fields(line)
    if not parts:
        return None
    parts = [x for x in parts if x != ""]
    return parts if len(parts) in (2, 3) else None


def _header(line):
    """Where each column is, when this line names its columns.

    {"sep": separator, "sequence": i, and "name": i, "ring": i or both}, or
    None. Every filled field must be a column's name, and one of them the
    sequence's.
    """
    sep, parts = _fields(line)
    if not parts:
        return None
    cols = {}
    for i, x in enumerate(parts):
        w = x.lower()
        if not w:
            continue
        role = ("sequence" if w in _SEQ_WORDS else "name" if w in _NAME_WORDS
                else "ring" if w in _RING_COLUMN else None)
        if role is None or role in cols:
            return None
        cols[role] = i
    if "sequence" not in cols or len(cols) < 2:
        return None
    cols["sep"] = sep
    return cols


def _seq_and_name(pair, header_hint, **kw):
    """(sequence index, name index, Problem or None) for two fields, one of
    which is the sequence; None when neither reads as one."""
    first = _reads_as_sequence(pair[0], **kw)
    second = _reads_as_sequence(pair[1], **kw)
    if first != second:
        # exactly one of them is a sequence, so there is nothing to decide
        return ((0, 1) if first else (1, 0)) + (None,)
    if not (first or second):
        return None
    # Both spell out in amino acid letters, which ordinary words do:
    # "native" is N A T I V E. Case settles most of them, because a sequence
    # is written in capitals and a label usually is not.
    up = [x.isupper() for x in pair]
    if up[0] != up[1]:
        return ((0, 1) if up[0] else (1, 0)) + (None,)
    # Still tied, so the first is taken as the sequence and the row says so
    # rather than choosing silently.
    return 0, 1, Problem(
        "ambiguous_columns",
        "Both columns read as sequences, so the first was taken as the "
        "sequence and '%s' as the name. Put a header row, %s, to say which is "
        "which." % (pair[1], header_hint))


def _with_ring(cell, seq, probs):
    """(ring, problems) for a row whose ring cell says `cell`. The ring is ""
    for an empty cell, which leaves the peptide to the run's own choice, and
    for a cell that names no ring, which is a problem said in the row."""
    got = ring_choice(cell)
    if got is None:
        return "", list(probs) + [ring_problem(cell)]
    if got and seq:
        pr = cyclic_problem(seq, got)
        if pr:
            probs = list(probs) + [pr]
    return got, probs


def _three(parts, **kw):
    """(sequence index, name index, ring index, Problem or None) for a row
    of three fields and no header, or None.

    The ring is the field that names one, the last such when more than one
    could; the other two are a sequence and its name. When no field names a
    ring, the field left over once a sequence and a name are found is taken
    as the ring the row meant, and is said to be wrong rather than dropped.
    """
    for want_ring in (True, False):
        for r in (2, 1, 0):
            if want_ring and ring_choice(parts[r]) is None:
                continue
            rest = [i for i in (0, 1, 2) if i != r]
            got = _seq_and_name([parts[rest[0]], parts[rest[1]]],
                                "'name,sequence,ring'", **kw)
            if got:
                si, ni, pr = got
                return rest[si], rest[ni], r, pr
    return None


def parse_rows(text, **kw):
    """One sequence per line, `name: SEQUENCE`, a row of columns, or FASTA.

    A row separated by a tab, a comma or a semicolon holds the sequence and
    its name, in either order: the field that reads as a sequence is taken
    as the sequence and the other as the name. A third field is the
    peptide's ring, written as linear, head to tail or disulfide (or one of
    the other names ring_choice reads); the ring may also stand alone beside
    the sequence when the row has no name. A header row naming the columns
    (sequence, name, ring) settles it for the whole list, which is what a
    spreadsheet export carries anyway.

    Returns [{"name", "sequence", "cyclic", "named", "problems"}]. "cyclic"
    is "none", "backbone" or "disulfide" for a row that names its ring, and
    "" for one that does not, which is docked as the run is set. "named" is
    whether the text gave the name, rather than it being the sequence's
    first letters. Blank lines and
    `#` comments are dropped; duplicates are kept, since two identical
    peptides in a list is a mistake worth reporting rather than silently
    fixing.

    `name: SEQUENCE` is the one-line way to name a peptide. A colon cannot be
    part of a sequence, so unlike a space it cannot be mistaken for one.
    """
    out = []
    name = None
    buf = []
    head = None

    def add(given, seq, ring, probs, fallback=None):
        out.append({"name": given or fallback or seq[:12],
                    "sequence": seq, "cyclic": ring, "named": bool(given),
                    "problems": list(probs)})

    def flush():
        if buf or name:
            seq, probs = check("".join(buf), **kw)
            add(name, seq, "", probs,
                fallback=seq[:12] or "peptide %d" % (len(out) + 1))

    for line in str(text or "").splitlines():
        if line.lstrip().startswith("#"):
            continue
        if line.startswith(">"):
            flush()
            name, buf[:] = line[1:].strip() or None, []
            continue
        if not line.strip():
            if name is None:
                continue
            flush()
            name, buf[:] = None, []
            continue
        # A header row fixes the columns for everything after it, so a list
        # whose names happen to be spellable in amino acid letters is not
        # guessed at row by row.
        if name is None and not out and head is None:
            got = _header(line)
            if got:
                head = got
                continue
        if name is None and head is not None and head["sep"] in line:
            cells = [x.strip() for x in line.split(head["sep"])]
            if not any(cells):
                continue                  # a spreadsheet's empty row

            def cell(role):
                i = head.get(role)
                return cells[i] if i is not None and i < len(cells) else ""
            seq, probs = check(cell("sequence"), **kw)
            ring, probs = _with_ring(cell("ring"), seq, probs)
            add(cell("name"), seq, ring, probs)
            continue
        if name is None:
            parts = _split_row(line)
            if parts and len(parts) == 2:
                # A sequence and its ring, with no name: the other field
                # names a ring and the first reads as the sequence.
                rings = [i for i in (0, 1) if ring_choice(parts[i])]
                if len(rings) == 1 and _reads_as_sequence(
                        parts[1 - rings[0]], **kw):
                    seq, probs = check(parts[1 - rings[0]], **kw)
                    ring, probs = _with_ring(parts[rings[0]], seq, probs)
                    add(None, seq, ring, probs)
                    continue
                got = _seq_and_name(parts, "'sequence,name' or "
                                    "'name,sequence'", **kw)
                if got:
                    si, ni, pr = got
                    seq, probs = check(parts[si], **kw)
                    add(parts[ni], seq, "",
                        list(probs) + ([pr] if pr else []))
                    continue
                # neither is a sequence: fall through and report on the line
                # as a whole rather than inventing a split
            elif parts:
                got = _three(parts, **kw)
                if got:
                    si, ni, ri, pr = got
                    seq, probs = check(parts[si], **kw)
                    ring, probs = _with_ring(parts[ri], seq, probs)
                    add(parts[ni], seq, ring,
                        list(probs) + ([pr] if pr else []))
                    continue
        if name is None and ":" in line:
            label, _c, rest = line.partition(":")
            if label.strip() and rest.strip():
                seq, probs = check(rest, **kw)
                add(label.strip(), seq, "", probs)
                continue
        if name is None:
            seq, probs = check(line, **kw)
            if looks_like_a_name(line):
                probs = list(probs) + [Problem(
                    "name_in_line",
                    "This looks like a name followed by a sequence, and it "
                    "would be docked as one longer peptide (%s). To name a "
                    "peptide, put a colon after the name ('%s: %s') or use "
                    "FASTA." % (seq, line.split()[0], " ".join(line.split()[1:])))]
            add(None, seq, "", probs, fallback=seq[:12] or line.strip()[:12])
        else:
            buf.append(line)
    if name is not None:
        flush()
    return out


def parse_many(text, **kw):
    """parse_rows as [(name, seq, problems)], for a caller that has no use
    for the ring column."""
    return [(r["name"], r["sequence"], r["problems"])
            for r in parse_rows(text, **kw)]
