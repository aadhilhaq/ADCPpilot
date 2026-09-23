#!/usr/bin/env python3
"""Keep adcp_dock/_vendor in step with the HOPE peptide pipeline's hope_core.

    python scripts/refresh_vendor.py --audit                  (no checkout needed)
    python scripts/refresh_vendor.py --from <HOPE checkout> --check
    python scripts/refresh_vendor.py --from <HOPE checkout>

The vendored modules carry hope_core's code unchanged. Their comments are
edited for this repository, so the comparison is by code: comments and
docstrings are ignored, and a module counts as changed only when its code
differs. A refresh copies only the modules whose code changed upstream;
review the comments of any file it copies.

estimate.py is an excerpt: only the names listed in EXCERPTS are kept, and
only those names are compared.

The modules import each other relatively (`from .receptor import
is_hydrogen`), so no import needs rewriting.
"""
import argparse
import ast
import os
import shutil
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VENDOR = os.path.join(HERE, "adcp_dock", "_vendor")

#: What is taken, and what each one is for.
MODULES = {
    "prepare.py": "the full receptor preparation, through to minimisation",
    "prep.py": "the entry point the pipeline calls the preparation through",
    "receptor.py": "chains, the site shell, the box, and the PDBQT",
    "adcp.py": "building the grid target and running the docking search",
    "mmgbsa.py": "the rescoring, with chain-ordered masks",
    "poseviewer.py": "a docked pose laid on a clean topology, matched by key",
    "maewriter.py": "Maestro output, which poseviewer writes through",
    "library.py": "peptides built from sequence",
    "estimate.py": "how many dockings to run at once (an excerpt)",
}

#: Modules kept in part. Only these top-level names are compared.
EXCERPTS = {
    "estimate.py": ("CORES_PER_DOCKING", "plan_concurrency"),
}

#: Deliberately not taken.
SKIPPED = {
    "tools.py": "adcp_dock/tools.py also finds tleap, sander and MMPBSA.py",
    "pairs.py": "pose pairing belongs to a growth loop this tool does not run",
    "sites.py": "site selection for the same loop",
}


def _code(tree, names=None):
    """The module's code with docstrings removed, optionally only `names`."""
    body = []
    for node in tree.body:
        if names is not None:
            named = (getattr(node, "name", None) in names
                     or (isinstance(node, ast.Assign)
                         and any(getattr(t, "id", None) in names
                                 for t in node.targets)))
            if not named:
                continue
        body.append(node)
    for node in ast.walk(ast.Module(body=body, type_ignores=[])):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef,
                             ast.ClassDef)) and node.body:
            first = node.body[0]
            if (isinstance(first, ast.Expr)
                    and isinstance(getattr(first, "value", None), ast.Constant)
                    and isinstance(first.value.value, str)):
                node.body[0] = ast.Pass()
    if names is None and body and isinstance(body[0], ast.Expr) \
            and isinstance(getattr(body[0], "value", None), ast.Constant):
        body = body[1:]
    return ast.dump(ast.Module(body=body, type_ignores=[]))


def same_code(upstream, local, names=None):
    a = ast.parse(open(upstream, encoding="utf-8").read())
    b = ast.parse(open(local, encoding="utf-8").read())
    return _code(a, names) == _code(b, names)


def audit():
    """The vendored directory holds exactly what MODULES lists.

    A listed module that is missing fails at import, at the far end of a
    queued job. A module that is present but not listed is left behind by
    the next refresh while everything around it is updated.
    """
    here = sorted(f for f in os.listdir(VENDOR)
                  if f.endswith(".py") and f != "__init__.py")
    missing = [n for n in sorted(MODULES) if n not in here]
    extra = [n for n in here if n not in MODULES]
    for n in sorted(MODULES):
        print("  %-9s %-16s %s" % ("present" if n in here else "MISSING",
                                   n, MODULES[n]))
    for n in extra:
        print("  %-9s %-16s not on the list" % ("EXTRA", n))

    # Nothing outside _vendor may shadow a vendored module.
    pkg = os.path.join(HERE, "adcp_dock")
    local = {f for f in os.listdir(pkg) if f.endswith(".py")}
    shadow = sorted(local & set(MODULES))

    bad = 0
    if missing:
        print("\n%d module(s) are listed and not present." % len(missing))
        bad += 1
    if extra:
        print("\n%d module(s) are present and not listed." % len(extra))
        bad += 1
    if shadow:
        print("\nthese exist both in adcp_dock/ and in _vendor/: %s"
              % ", ".join(shadow))
        bad += 1
    if not bad:
        print("\n%d vendored module(s), all accounted for." % len(MODULES))
    return 1 if bad else 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--from", dest="src",
                    help="a HOPE pipeline checkout, or its hope_core/")
    ap.add_argument("--audit", action="store_true",
                    help="check the vendored directory against the list; "
                         "needs no HOPE checkout")
    ap.add_argument("--check", action="store_true",
                    help="report code differences and change nothing")
    a = ap.parse_args(argv)

    if a.audit:
        return audit()
    if not a.src:
        print("give --from <HOPE checkout>, or --audit")
        return 2

    src = a.src
    if os.path.isdir(os.path.join(src, "hope_core")):
        src = os.path.join(src, "hope_core")
    if not os.path.isdir(src):
        print("no such directory: %s" % src)
        return 2

    same, changed, missing = [], [], []
    for name in sorted(MODULES):
        s, d = os.path.join(src, name), os.path.join(VENDOR, name)
        if not os.path.isfile(s):
            missing.append(name)
            continue
        if os.path.isfile(d) and same_code(s, d, EXCERPTS.get(name)):
            same.append(name)
            continue
        changed.append(name)
        if not a.check and name not in EXCERPTS:
            shutil.copy2(s, d)

    for n in same:
        print("  same      %s" % n)
    for n in changed:
        verb = "differs  " if a.check or n in EXCERPTS else "updated  "
        print("  %s %-16s %s" % (verb, n, MODULES[n]))
    for n in missing:
        print("  MISSING   %-16s not in %s" % (n, src))
    for n, why in sorted(SKIPPED.items()):
        print("  skipped   %-16s %s" % (n, why))

    if missing:
        print("\n%d module(s) are not in that checkout." % len(missing))
        return 1
    excerpt_changed = [n for n in changed if n in EXCERPTS]
    if excerpt_changed:
        print("\nupdate by hand: %s (only %s are kept)"
              % (", ".join(excerpt_changed),
                 "; ".join(", ".join(EXCERPTS[n]) for n in excerpt_changed)))
    if a.check and changed:
        print("\n%d vendored module(s) differ in code from %s. Run without "
              "--check to take the newer code, then re-run the tests."
              % (len(changed), src))
        return 1
    if changed and not a.check:
        print("\nCopied files carry the upstream comments; review them.")
    print("\n%d vendored, %d unchanged." % (len(MODULES), len(same)))
    return 1 if excerpt_changed else 0


if __name__ == "__main__":
    sys.exit(main())
