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
# the fifth a D-aspartate. ADCP 1.0 has no such syntax.
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
# or through a disulfide between two of its cysteines (-cys). One choice per
# run, applied to every peptide in it.
CYCLIC = ("none", "backbone", "disulfide")
#: What each choice is called in the logs and the batch scripts.
CYCLIC_LABELS = {"none": "linear", "backbone": "backbone cyclic, head to tail",
                 "disulfide": "cyclic through a disulfide"}


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
    cysteine or none, ADCP 1.1.21 and 1.0 both docked it and reported poses
    with no warning, so a run would describe as cyclic a peptide that was
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
                "linear or backbone cyclic." % ("none" if n == 0 else "one"))
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


def check(text, min_length=MIN_LENGTH, max_length=MAX_LENGTH, extended=False,
          library=None):
    """(cleaned, [Problem]). An empty problem list means ADCP will take it.

    With `extended` (ADCP 1.1), the non-standard syntax is accepted, and each
    side-chain name is checked against `library` when one is given.
    """
    seq = clean(text)
    out = []
    if not seq:
        return seq, [Problem("empty", "Enter a peptide sequence.")]
    if is_extended(seq):
        return seq, _check_extended(seq, min_length, max_length, extended,
                                    library)

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


def _check_extended(seq, min_length, max_length, extended, library):
    if not extended:
        return [Problem("needs_adcp11",
                        "Non-standard (<XXX>) and D (&) amino acids need "
                        "ADCP 1.1; this run uses ADCP 1.0.")]
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
# is does not have to be guessed at.
_NAME_WORDS = ("name", "label", "id", "peptide_name", "construct")
_SEQ_WORDS = ("sequence", "seq", "peptide", "aa", "residues")


def _reads_as_sequence(text, **kw):
    """True when this field is a usable sequence on its own.

    Used to work out which half of a two column row is the sequence, so both
    `SEQUENCE,name` and `name,SEQUENCE` work without being told which.
    """
    t = (text or "").strip()
    if not t:
        return False
    seq, probs = check(t, **kw)
    return bool(seq) and not probs


def _split_pair(line):
    """A line split on a tab, a comma or a semicolon, if it is in two parts.

    A colon is handled separately and earlier, because `name: SEQUENCE` is
    the documented one line form. Only a clean two field split counts: three
    or more fields is a spreadsheet row this function should not guess at.
    """
    for sep in ("\t", ",", ";"):
        if sep in line:
            parts = [x.strip() for x in line.split(sep)]
            parts = [x for x in parts if x != ""]
            if len(parts) == 2:
                return parts
            return None
    return None


def _header_order(line):
    """(sequence index, name index) if this line names its columns."""
    parts = _split_pair(line)
    if not parts:
        return None
    low = [x.strip().lower() for x in parts]
    if low[0] in _SEQ_WORDS and low[1] in _NAME_WORDS:
        return 0, 1
    if low[0] in _NAME_WORDS and low[1] in _SEQ_WORDS:
        return 1, 0
    return None


def parse_many(text, **kw):
    """One sequence per line, `name: SEQUENCE`, a two column row, or FASTA.

    A two column row separated by a tab, a comma or a semicolon works in
    either order: the field that reads as a sequence is taken as the
    sequence and the other as the name. A header row naming the columns
    settles it for the whole list, which is what a spreadsheet export
    carries anyway.

    Returns [(name, seq, problems)]. Blank lines and `#` comments are dropped;
    duplicates are kept, since two identical peptides in a list is a mistake
    worth reporting rather than silently fixing.

    `name: SEQUENCE` is the one-line way to name a peptide. A colon cannot be
    part of a sequence, so unlike a space it cannot be mistaken for one.
    """
    out = []
    name = None
    buf = []
    order = None

    def flush():
        if buf or name:
            seq, probs = check("".join(buf), **kw)
            out.append((name or seq[:12] or "peptide %d" % (len(out) + 1),
                        seq, probs))

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
        # A header row fixes the order for everything after it, so a list
        # whose names happen to be spellable in amino acid letters is not
        # guessed at row by row.
        if name is None and not out and order is None:
            got = _header_order(line)
            if got:
                order = got
                continue
        if name is None:
            pair = _split_pair(line)
            if pair:
                if order is not None:
                    si, ni = order
                    seq, probs = check(pair[si], **kw)
                    out.append((pair[ni].strip() or seq[:12], seq, probs))
                    continue
                first = _reads_as_sequence(pair[0], **kw)
                second = _reads_as_sequence(pair[1], **kw)
                if first != second:
                    # exactly one of them is a sequence, so there is nothing
                    # to decide
                    si, ni = (0, 1) if first else (1, 0)
                    seq, probs = check(pair[si], **kw)
                    out.append((pair[ni].strip() or seq[:12], seq, probs))
                    continue
                if first and second:
                    # Both spell out in amino acid letters, which ordinary
                    # words do: "native" is N A T I V E. Case settles most
                    # of them, because a sequence is written in capitals and
                    # a label usually is not.
                    up = [x.isupper() for x in pair]
                    if up[0] != up[1]:
                        si, ni = (0, 1) if up[0] else (1, 0)
                        seq, probs = check(pair[si], **kw)
                        out.append((pair[ni].strip() or seq[:12], seq, probs))
                        continue
                    # Still tied, so the first column is taken as the
                    # sequence and the row says so rather than choosing
                    # silently.
                    seq, probs = check(pair[0], **kw)
                    probs = list(probs) + [Problem(
                        "ambiguous_columns",
                        "Both columns read as sequences, so the first was "
                        "taken as the sequence and '%s' as the name. Put a "
                        "header row, 'sequence,name' or 'name,sequence', to "
                        "say which is which." % pair[1])]
                    out.append((pair[1].strip(), seq, probs))
                    continue
                # neither is a sequence: fall through and report on the line
                # as a whole rather than inventing a split
        if name is None and ":" in line:
            label, _c, rest = line.partition(":")
            if label.strip() and rest.strip():
                seq, probs = check(rest, **kw)
                out.append((label.strip(), seq, probs))
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
            out.append((seq[:12] or line.strip()[:12], seq, probs))
        else:
            buf.append(line)
    if name is not None:
        flush()
    return out
