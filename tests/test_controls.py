"""Every control a person must press to make progress is plainly a button.

By Aadhil Haq

The page has one hierarchy of buttons. A filled primary is the one action a
card's step is waiting for; an outlined secondary, in the accent, at least
36 px tall with 14 px text, is any other action someone needs; a ghost,
outlined in the ink, is a utility such as close, copy or refresh; and a link
is only for what nobody needs in order to go on. A label, such as the "run"
beside a folder, is not shaped like a button at all. This reads the page and
its script to hold every control to that, and checks that the colors that
carry the words are readable in both themes.

    python3 tests/test_controls.py
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import page  # noqa: E402

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


html = page.render("TOKEN")
appjs = open(os.path.join(ROOT, "adcp_dock", "static", "app.js")).read()
css = html[html.index("<style>"):html.index("</style>")]


def props(sel):
    """The declarations a selector gets from the rules that name it alone."""
    out = {}
    for m in re.finditer(r"([^{}]+)\{([^}]*)\}", re.sub(r"/\*.*?\*/", "", css, flags=re.S)):
        if sel in [x.strip() for x in m.group(1).split(",")]:
            out.update(dict(re.findall(r"([\w-]+)\s*:\s*([^;]+)", m.group(2))))
    return {k: v.strip() for k, v in out.items()}


def classes_of(id_):
    m = re.search(r"<(?:button|a)\b[^>]*\bid=%s\b[^>]*>" % re.escape(id_), html)
    if not m:
        return None
    c = re.search(r'class=(?:"([^"]*)"|([\w-]+))', m.group(0))
    return set((c.group(1) or c.group(2)).split()) if c else set()


# --- the hierarchy itself -------------------------------------------------------------
sec = props("button.secondary")
check(sec.get("min-height") == "36px" and sec.get("font-size") == "14px",
      "a secondary button is not at least 36 px tall at 14 px: %r" % sec)
check("var(--accent)" in sec.get("border", "") and sec.get("color") == "var(--accent)",
      "a secondary button is not outlined in the accent: %r" % sec)
pri = props("button.primary")
check(pri.get("background") == "var(--accent)" and pri.get("color") == "var(--on-accent)",
      "a primary button is not filled, with the ink the theme gives a filled control: %r" % pri)
ghost = props("button.ghost")
check(ghost.get("color") == "var(--fg)",
      "a ghost button is still pale gray text: %r" % ghost)
check("var(--pill-line)" in props("button").get("border", ""),
      "a plain button's edge is the faint line color")
check(props(".btnlink").get("color") == "var(--fg)", "a link styled as a button is pale")
check("button:focus-visible" in css and ".ld-retry:focus-visible" in css,
      "the buttons do not show the keyboard's focus")

# --- one primary per card, for its step's main action -----------------------------------
MAIN = {"receptor": "recfetch", "site": "sitecheck", "peptides": "addpep",
        "docking": "autofill", "submit": "submit"}
# Each card up to whatever comes next: another card, a floating panel or a
# dialog, none of which belongs to it.
cards = re.split(r"(?=<div class=card|<div class=\"card|<div id=issuebox|<div id=helpbox|<dialog)",
                 html)
for step, main in MAIN.items():
    body = next((c for c in cards if c.startswith("<div class=card id=step-%s " % step)), "")
    shown = re.findall(r'<button class=(?:primary|"primary(?! hide)[^"]*") id=(\w+)', body)
    check(shown == [main], "the %s card's filled buttons are %r, not just %s"
          % (step, shown, main))

# --- what someone needs is never a ghost ------------------------------------------------
NEEDED = ["recfetch", "recload", "recbrowse", "recupload", "trgload", "trgbrowse",
          "sitecheck", "addpep", "pastepep", "pepfilebtn", "pasteadd", "autofill",
          "outbrowse", "submit", "adaptnow", "adaptdrop", "adaptopen", "adaptupdate",
          "adaptread", "adaptback", "adaptfull", "runsbrowse", "measuresopen", "msgo",
          "pickeruse", "modsave",
          "filedl", "foldercopy", "fbcopy", "fbmail", "fbissue"]
for id_ in NEEDED:
    c = classes_of(id_)
    check(c is not None and ("primary" in c or "secondary" in c),
          "%s, which someone needs, is %r rather than a primary or a secondary"
          % (id_, c))
# And the ones the script draws: open a run, open it from Folders, view a
# peptide, add what a run lacks, run the missing ones, try again.
for frag in ("'<button class=secondary data-open=\"'",
             "'<button class=secondary data-openrun=\"'",
             "'<button class=secondary data-view=\"'",
             "'.</span><button class=secondary type=button id=msadd data-ms=\"'",
             "'<button class=secondary id=resume",
             "'<button class=secondary data-openrun=\"' + esc(d.path)",
             "'<button type=button class=\"secondary ld-retry hide\">Try again</button>'",
             "'<button type=button class=\"primary ld-retry hide\">Try again</button></div>'"):
    check(frag.strip("'") in appjs, "the script draws %s as something else" % frag)
for frag in ('class="ghost sm" data-openrun', 'class=ghost data-open=',
             'class="ghost sm" data-view=', 'class="ghost sm" id=msadd',
             'class=ghost id=measuresopen'):
    check(frag not in appjs, "a needed action is still a ghost: %s" % frag)

# --- a label is not a button ------------------------------------------------------------
check("'<span class=tag title=\"this folder is a run\">run</span>'" in appjs
      and 'class=pill title="this folder is a run"' not in appjs,
      "the run label in Folders is still shaped like a button")
tag = props(".tag")
check(tag.get("cursor") == "default" and "border" not in tag,
      "a label has a pointer or a border, like a button: %r" % tag)
check("<span class=tag>' + esc(c.kind)" in appjs,
      "a chain's description is still shaped like a button")

# --- the viewer's toolbar ----------------------------------------------------------------
vp = props(".vbar .pill")
check(vp.get("color") == "var(--fg)" and "var(--pill-line)" in vp.get("border", ""),
      "the viewer's toolbar is still pale: %r" % vp)

# --- readable in both themes -------------------------------------------------------------
HEX = r"--([\w-]+):(#[0-9a-fA-F]{6}\b|#[0-9a-fA-F]{3}\b)"
LIGHT = dict(re.findall(HEX, css[:css.index("@media")]))
DARK = dict(re.findall(HEX, page.DARK_VARS))


def lum(h):
    if len(h) == 4:
        h = "#" + "".join(c * 2 for c in h[1:])
    r, g, b = (int(h[i:i + 2], 16) / 255 for i in (1, 3, 5))
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4  # noqa: E731
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b)


def contrast(a, b):
    la, lb = lum(a), lum(b)
    return (max(la, lb) + 0.05) / (min(la, lb) + 0.05)


for theme, t in (("light", LIGHT), ("dark", DARK)):
    for ink, ground in (("accent", "card"), ("on-accent", "accent"), ("fg", "card"),
                        ("fg", "accent2"), ("dim", "card"), ("bad", "card")):
        c = contrast(t[ink], t[ground])
        check(c >= 4.5, "in the %s theme --%s on --%s is %.2f:1, under 4.5"
              % (theme, ink, ground, c))

if fails:
    print("test_controls FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  every needed control is plainly a button, one primary per card, %d checks" % ran)
