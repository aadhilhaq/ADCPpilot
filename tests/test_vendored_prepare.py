"""The receptor preparation in adcp_dock/_vendor is hope_core's code.

By Aadhil Haq

adcp_dock/_vendor/prepare.py is the preparation HOPE-pipelines maintains as
hope_core/prepare.py, carried here with this tool's own docstrings and
comments. So the two are compared as code rather than as text: both are
parsed, every docstring is left out (a parse holds no comments), and what is
left must match statement for statement. A difference fails, naming the
first top-level definition where the two part, so a fix made in one and not
in the other is found here rather than in a run.

hope_core/prepare.py is looked for in a HOPE-pipelines checkout beside this
one or beside any folder above it, which is how an installation keeps its
tools side by side. HOPE_PIPELINES_CHECKOUT names another checkout, and is
then the only place looked. Where none is found the test is a skip naming
the places it looked, never a pass.

    python3 tests/test_vendored_prepare.py
"""
import ast
import difflib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MINE = os.path.join(ROOT, "adcp_dock", "_vendor", "prepare.py")

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


def checkouts():
    """The HOPE-pipelines checkouts to look in, in order."""
    named = os.environ.get("HOPE_PIPELINES_CHECKOUT", "").strip()
    if named:
        return [named]
    out, folder = [], os.path.dirname(ROOT)
    while True:
        out.append(os.path.join(folder, "HOPE-pipelines"))
        if os.path.dirname(folder) == folder:
            return out
        folder = os.path.dirname(folder)


def name_of(stmt):
    """A top-level statement as a reader would look for it in the file."""
    if isinstance(stmt, (ast.FunctionDef, ast.AsyncFunctionDef)):
        return "def %s()" % stmt.name
    if isinstance(stmt, ast.ClassDef):
        return "class %s" % stmt.name
    if isinstance(stmt, ast.Assign):
        return " = ".join(ast.unparse(t) for t in stmt.targets) + " ="
    if isinstance(stmt, (ast.AnnAssign, ast.AugAssign)):
        return ast.unparse(stmt.target) + " ="
    return "`%s`" % ast.unparse(stmt).splitlines()[0][:60]


def code_of(path):
    """(name, line, code) for each top-level statement of a module, the code
    being its parse with every docstring left out."""
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), filename=path)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                             ast.AsyncFunctionDef)):
            body = node.body
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                node.body = body[1:]
    return [(name_of(s), s.lineno, ast.dump(s)) for s in tree.body]


places = checkouts()
master = next((m for m in (os.path.join(c, "hope_core", "prepare.py")
                           for c in places) if os.path.isfile(m)), None)
if not master:
    if os.environ.get("HOPE_PIPELINES_CHECKOUT", "").strip():
        print("skip  adcp_dock/_vendor/prepare.py was compared with nothing: "
              "%s, the checkout HOPE_PIPELINES_CHECKOUT names, holds no "
              "hope_core/prepare.py" % places[0])
    else:
        print("skip  adcp_dock/_vendor/prepare.py was compared with nothing: "
              "no HOPE-pipelines checkout holding hope_core/prepare.py is "
              "beside this one or beside a folder above it (looked at %s), "
              "and HOPE_PIPELINES_CHECKOUT, which names one elsewhere, is "
              "not set" % ", ".join(places))
    sys.exit(77)

mine, theirs = code_of(MINE), code_of(master)
ops = [op for op in difflib.SequenceMatcher(
    None, [c for _n, _l, c in mine], [c for _n, _l, c in theirs],
    autojunk=False).get_opcodes() if op[0] != "equal"]
where, more = "", 0
if ops:
    tag, i1, i2, j1, j2 = ops[0]
    if tag == "delete":
        where = "%s (line %d) is here and not there" % mine[i1][:2]
    elif tag == "insert":
        where = "%s (line %d there) is not here" % theirs[j1][:2]
    elif mine[i1][0] == theirs[j1][0]:
        where = "%s differs (line %d here, line %d there)" % (
            mine[i1][0], mine[i1][1], theirs[j1][1])
    else:
        where = "%s (line %d) stands where it has %s (line %d)" % (
            mine[i1][:2] + theirs[j1][:2])
    more = sum(max(b - a, d - c) for _t, a, b, c, d in ops) - 1
check(not ops,
      "adcp_dock/_vendor/prepare.py is not the code of %s: %s%s. Copy it "
      "again and put this tool's docstrings back, rather than edit it here"
      % (master, where,
         "; %d more top-level statement(s) differ after it" % more if more
         else ""))

if fails:
    print("test_vendored_prepare FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("note  compared with %s, docstrings and comments left out" % master)
print("ok  adcp_dock/_vendor/prepare.py is hope_core's code, all %d top-level "
      "statements" % len(mine))
