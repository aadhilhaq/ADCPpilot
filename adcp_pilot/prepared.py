"""A receptor that arrives prepared: checked, taken in as it is, and recorded.

By Aadhil Haq

A receptor prepared elsewhere - in a HOPE Labs flow, by HOPE-MD's Prepare at the flow's pH - is
docked with the protonation it was given. Preparing it again here would throw that away: the
preparation strips every hydrogen and protonates again, so the histidines, the disulfides and the
neutral acids and bases it decided would be decided a second time, by another program, at
whatever pH this run happened to have.

So with receptor_prepared on, nothing of this package's preparation runs. What is still done is
the part that is not a decision: the chains asked for are taken out, the waters left behind, and
the hydrogens counted, because a receptor said to be prepared and carrying none is not prepared.
It is refused, with the reason, rather than protonated quietly.

What a prepared receptor looks like: hydrogens on it, the Amber names that carry a state (HID, HIE,
HIP, CYX, CYM, ASH, GLH, LYN), its own chain letters and numbering, and ACE and NME caps where a
break was capped. prepare_receptor types every one of those as it types the standard residue
(tests/test_receptor_prepared.py compares the two), so the names are passed to it unchanged.
"""
import json
import os

#: The residue names that carry a protonation state or a disulfide, as the Amber force fields
#: and PDB2PQR name them. Counted for the record; nothing here renames them.
STATE_NAMES = ("HID", "HIE", "HIP", "CYX", "CYM", "ASH", "GLH", "LYN")
#: The caps a preparation puts on the ends a chain break leaves. NMA is Maestro's name for NME.
CAPS = ("ACE", "NME", "NMA")
WATERS = ("HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL")

FLAG = "--receptor-prepared"


class NotPrepared(ValueError):
    """A receptor said to be prepared that is not."""


def _is_h(line):
    element = line[76:78].strip().upper() if len(line) > 76 else ""
    if element:
        return element in ("H", "D")
    name = line[12:16].strip()
    return name.lstrip("0123456789")[:1].upper() in ("H", "D")


def _chain(line):
    return line[21] if line[21] != " " else "_"


def _wanted(chains):
    return {c.strip() for c in (chains or []) if str(c).strip()}


def polymer_lines(path, chains=()):
    """The receptor's lines from the first model: every ATOM record of the chains asked for (all
    of them when none is), and the caps whatever record they were written as."""
    want = _wanted(chains)
    out = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if want and _chain(line) not in want:
                continue
            if line.startswith("HETATM") and line[17:20].strip() not in CAPS:
                continue
            if line[16] not in (" ", "A"):
                continue
            out.append(line)
    return out


def hydrogens(path, chains=()):
    """{chain: (atoms, hydrogens)} over the receptor's polymer, chain by chain."""
    got = {}
    for line in polymer_lines(path, chains):
        n, h = got.get(_chain(line), (0, 0))
        got[_chain(line)] = (n + 1, h + (1 if _is_h(line) else 0))
    return got


def problems(path, chains=()):
    """[str] why this receptor cannot be docked as a prepared one; empty when it can."""
    if not path or not os.path.isfile(path):
        return []
    if os.path.splitext(path)[1].lower() not in (".pdb", ".ent"):
        return ["A prepared receptor is read as a PDB file, and %s is not one. Write the "
                "prepared receptor as PDB, or dock this file without %s so it is prepared "
                "here." % (os.path.basename(path), FLAG)]
    try:
        got = hydrogens(path, chains)
    except OSError as exc:
        return ["The prepared receptor could not be read: %s" % exc]
    if not got:
        return ["The prepared receptor %s has no atoms on chain(s) %s."
                % (os.path.basename(path), ", ".join(sorted(_wanted(chains))) or "any")]
    bare = sorted(c for c, (_n, h) in got.items() if h == 0)
    if bare:
        return ["The receptor is marked as prepared (%s), but chain(s) %s of %s carry no "
                "hydrogens. A prepared receptor holds the hydrogens its preparation placed, and "
                "without them its protonation states are not known. Prepare it first (HOPE-MD's "
                "Prepare does it at the flow's pH), or dock it without %s so it is prepared here."
                % (FLAG, ", ".join(bare), os.path.basename(path), FLAG)]
    return []


def states(path, chains=()):
    """{name: count} of the residues whose names carry a state, and of the caps."""
    seen, out = set(), {}
    for line in polymer_lines(path, chains):
        key = (_chain(line), line[22:27])
        name = line[17:20].strip()
        if key in seen:
            continue
        seen.add(key)
        if name in STATE_NAMES or name in CAPS:
            out[name] = out.get(name, 0) + 1
    return out


def provenance(path):
    """What the file, or a record beside it, says about how it was prepared: {} when nothing does.

    HOPE Labs' preparation card leaves prepared.json in its own folder, one level above each
    structure it prepared, naming the pH and how the molecules were protonated. A file of the
    same name with .json added, beside it, is read the same way, and so are REMARK lines in the
    file that say who prepared it.
    """
    real = os.path.realpath(path)
    said = {}
    remarks = []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    break
                if line.startswith("REMARK") and "prepar" in line.lower():
                    remarks.append(line[6:].strip())
    except OSError:
        pass
    if remarks:
        said["remarks"] = remarks[:6]
    here = os.path.dirname(real)
    for cand in (real + ".json", os.path.splitext(real)[0] + ".json",
                 os.path.join(here, "prepared.json"),
                 os.path.join(os.path.dirname(here), "prepared.json")):
        try:
            rec = json.load(open(cand, encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if not isinstance(rec, dict):
            continue
        mine = None
        for s in rec.get("structures") or []:
            if isinstance(s, dict) and s.get("path") and \
                    os.path.realpath(s["path"]) == real:
                mine = s
                break
        if rec.get("structures") is not None and mine is None:
            continue
        prep = rec.get("prep") or {}
        said["record"] = cand
        if rec.get("ph") is not None:
            said["ph"] = float(rec["ph"])
        method = prep.get("protonation") or ""
        mols = (mine or {}).get("molecules") or []
        methods = sorted({(m.get("protonation") or {}).get("method") or "" for m in mols} - {""})
        said["method"] = ", ".join(methods) or method
        said["by"] = ("HOPE-MD's Prepare, through HOPE Labs' preparation card"
                      if mine is not None else rec.get("by") or "")
        if mine is not None:
            said["structure"] = mine.get("name") or ""
        break
    return said


def describe(info):
    """The one line a log or a report says about a receptor that came prepared."""
    by = info.get("by") or "a preparation done before this run"
    made = info.get("method")
    ph = info.get("ph")
    given = info.get("ph_given")
    bits = [by]
    if made:
        bits.append("protonated with %s" % ("PDB2PQR and PROPKA" if made == "propka" else made))
    if ph is not None:
        bits.append("at pH %g" % ph)
    elif given is not None:
        bits.append("at pH %g, as this run was told" % given)
    text = ", ".join(bits)
    st = info.get("states") or {}
    if st:
        text += "; " + ", ".join("%s x%d" % kv for kv in sorted(st.items()))
    return text


def take(path, chains, out, keep_hetatms=False):
    """Write the prepared receptor's chains to `out`, as they are. Returns {what was taken}.

    Lines are copied rather than rewritten, so every name, number, insertion code and
    coordinate is the preparation's own. A cap written as HETATM becomes ATOM, which is what
    the readers downstream take a receptor's residues from. Waters are left out; the other
    hetero groups go with them unless the run keeps them, as for a receptor prepared here.
    Refuses (NotPrepared) a receptor with no hydrogens on a chain.
    """
    bad = problems(path, chains)
    if bad:
        raise NotPrepared(bad[0])
    want = _wanted(chains)
    lines, hetero, dropped = [], [], {"waters": 0, "groups": set()}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if want and _chain(line) not in want:
                continue
            if line[16] not in (" ", "A"):
                continue
            name = line[17:20].strip()
            if line.startswith("HETATM"):
                if name in CAPS:
                    line = "ATOM  " + line[6:]
                elif name in WATERS:
                    dropped["waters"] += 1
                    continue
                elif keep_hetatms:
                    hetero.append(line)
                    continue
                else:
                    dropped["groups"].add("%s %s:%s" % (name, _chain(line),
                                                         line[22:27].strip()))
                    continue
            lines.append(line if line.endswith("\n") else line + "\n")
    with open(out, "w") as fh:
        last = None
        for line in lines:
            if last is not None and line[21] != last:
                fh.write("TER\n")
            fh.write(line)
            last = line[21]
        fh.write("TER\n")
        for line in hetero:
            fh.write(line if line.endswith("\n") else line + "\n")
        fh.write("END\n")
    got = hydrogens(out)
    took = {"atoms": sum(n for n, _h in got.values()),
            "hydrogens": sum(h for _n, h in got.values()),
            "chains": sorted(got), "states": states(out),
            "waters_dropped": dropped["waters"],
            "groups_dropped": sorted(dropped["groups"]),
            "hetero_kept": len(hetero)}
    return took
