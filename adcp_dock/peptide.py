"""What ADCP accepts as a peptide, and why a given string does not qualify.

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


def check(text, min_length=MIN_LENGTH, max_length=MAX_LENGTH):
    """(cleaned, [Problem]). An empty problem list means ADCP will take it."""
    seq = clean(text)
    out = []
    if not seq:
        return seq, [Problem("empty", "Enter a peptide sequence.")]

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


def parse_many(text, **kw):
    """One sequence per line, `name: SEQUENCE`, or a multi-record FASTA.

    Returns [(name, seq, problems)]. Blank lines and `#` comments are dropped;
    duplicates are kept, since two identical peptides in a list is a mistake
    worth reporting rather than silently fixing.

    `name: SEQUENCE` is the one-line way to name a peptide. A colon cannot be
    part of a sequence, so unlike a space it cannot be mistaken for one.
    """
    out = []
    name = None
    buf = []

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
