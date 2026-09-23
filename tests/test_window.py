"""The launcher window builds, and the site table fills it in.

Where there is a display (the Windows and macOS build runners have one)
the real window is constructed and driven. Where there is not, the checks
that need no display still run and this file says which half was skipped.

The source checks here cover only what a headless machine cannot
exercise; where a display exists they are never the whole test.
"""
import inspect
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from adcp_dock import brand                                     # noqa: E402
from adcp_launcher import clusters                              # noqa: E402

fails, ran, skipped = [], 0, 0


def check(cond, why):
    global ran
    ran += 1
    if not cond:
        fails.append(why)


# --- the brand is real, and its own -----------------------------------------
png = brand.logo_png()
check(png[:8] == b"\x89PNG\r\n\x1a\n",
      "brand.logo_png() is not a PNG, so the window would open with no mark")
check(len(png) > 2000, "the embedded mark is %d bytes; that is not a 128px "
                       "image" % len(png))
check(brand.NAME and brand.TAGLINE, "the brand has no name or tagline")
check("<svg" in brand.LOGO_SVG and "</svg>" in brand.LOGO_SVG,
      "the vector mark is not well formed")
check(brand.version(), "brand.version() returned nothing; the window shows "
                       "this and must not be blank")

# Three launchers sit side by side in one taskbar. This one's mark must not
# be a sibling's or the peptide pipeline's.
for other in ("hope_aptamer", "hope_core"):
    try:
        mod = __import__("%s.brand" % other, fromlist=["brand"])
    except ImportError:
        continue
    try:
        check(mod.logo_png() != png,
              "this launcher's mark is byte-identical to %s's: at 16 "
              "pixels in a taskbar they would be the same icon" % other)
        check(mod.NAME != brand.NAME,
              "this launcher has the same name as %s's" % other)
    except Exception:                                           # noqa: BLE001
        pass

# --- the three launchers cannot tread on each other -------------------------
from adcp_launcher import app as appmod                         # noqa: E402

check(".adcp" in os.path.basename(appmod.CONFIG),
      "the settings file is %r, which is not this tool's own. Three "
      "launchers sharing one file means each overwrites the others' "
      "cluster, account and paths" % appmod.CONFIG)
for other in (".hope_launcher.json", ".hope_aptamer_launcher.json"):
    check(os.path.basename(appmod.CONFIG) != other,
          "the settings file collides with %s" % other)

# The port band is this tool's own, and a sibling's formula is the one to
# stay clear of: 8000 + (sum of the account's characters mod 1000).
for user in ("aadhil.haq", "sandun", "a", "zzzzzzzz"):
    mine = clusters.suggested_port(user)
    theirs = 8000 + (sum(ord(c) for c in user) % 1000)
    check(mine != theirs,
          "%r gets port %d from this launcher and %d from a sibling "
          "launcher's formula; one person running both would collide"
          % (user, mine, theirs))
    check(clusters.PORT_BASE <= mine < clusters.PORT_BASE + clusters.PORT_SPAN,
          "%d is outside this tool's declared band" % mine)
check(clusters.PORT_BASE >= 9000,
      "the port band starts at %d, inside a sibling launcher's 8000 to 8999"
      % clusters.PORT_BASE)

# and the default is auto, which is the only thing that cannot collide
check(clusters.parse_port("") == 0 and clusters.parse_port("auto") == 0,
      "a blank or 'auto' port box does not mean 'let the server choose'")
check(clusters.parse_port("9443") == 9443, "a typed port is not honoured")
check(clusters.parse_port("80") == 0,
      "a privileged port was accepted; it cannot be bound and the failure "
      "would look like the server dying")
check(clusters.parse_port("nonsense") == 0,
      "nonsense in the port box is not falling back to auto")

# --- the window is driven by the table, not by hard-coded strings -----------
src = inspect.getsource(appmod.App)
check("clusters.defaults_for" in src or "_apply_site" in src,
      "the form does not consult the site table, so picking a cluster would "
      "fill in nothing")
check("brand.logo_png" in src or "logo_png" in src,
      "the window never loads the mark")
check("_touched" in src,
      "nothing tracks which boxes the person typed into, so changing the "
      "cluster would silently overwrite their own answers")

# --- the runs directory is NOT a box in the launcher ------------------------
check(not hasattr(appmod.App, "runs") and "self.runs" not in src,
      "the launcher still has a runs directory field. A user signing in "
      "from a laptop cannot know which scratch path they are meant to use, "
      "and a home directory is usually too small for a docking run")
check('"", ' in src or 'install.get().strip(), ""' in src,
      "the launcher is still sending a runs directory to the far end")

# --- the banner is aligned by arithmetic, not by eye --------------------------
# The header is drawn on a Canvas from banner_layout, which is tested here
# because Tk cannot run on the machine this is developed on. The things that
# go wrong in a banner are exactly the arithmetic: ragged left edges, text
# sitting high or low against the mark, and decoration drawn through words.
from adcp_launcher.app import banner_layout, Banner, BANNER    # noqa: E402

HT, WD = [30, 17, 15], [260, 520, 110]
for W in (680, 760, 1200):
    L = banner_layout(W, Banner.H, HT, WD)
    x1, y1, x2, y2 = L["panel"]
    check(abs((y1 + y2) / 2.0 - Banner.H / 2.0) < 0.51,
          "at %dpx the mark's panel is not centred vertically" % W)
    top = L["text_ys"][0]
    bottom = L["text_ys"][-1] + HT[-1]
    check(abs((top + bottom) / 2.0 - Banner.H / 2.0) < 0.51,
          "at %dpx the text block is not centred on the mark's centre line "
          "(top %.1f, bottom %.1f)" % (W, top, bottom))
    check(0 <= top and bottom <= Banner.H,
          "at %dpx the text runs outside the banner" % W)
    check(L["text_x"] > x2, "at %dpx the text overlaps the mark" % W)
    if L["motif"]:
        check(L["motif"][0] >= L["text_x"] + max(WD),
              "at %dpx the peptide motif starts inside the widest line" % W)
        check(L["motif"][1] <= W, "at %dpx the motif runs off the edge" % W)
# Too narrow for the motif: it is left out, not squeezed through the text.
L = banner_layout(680, Banner.H, HT, [260, 560, 110])
check(L["motif"] is None or L["motif"][0] >= L["text_x"] + 560,
      "with no room for it, the motif was squeezed in over the tagline")
check(banner_layout(1200, Banner.H, HT, WD)["motif"] is not None,
      "on a wide window there is room for the motif and it was left out")
check(all(BANNER[k].startswith("#") and len(BANNER[k]) == 7 for k in BANNER),
      "a banner colour is not a #rrggbb value")
check(BANNER["deep"] == brand.DEEP and BANNER["teal"] == brand.TEAL
      and BANNER["warm"] == brand.WARM,
      "the launcher banner's palette has drifted from the brand's")

# --- every dialog is placed over the window that raised it -------------------
#
# A dialog placed away from the launcher is easy to miss, and typing then
# goes into whichever window has focus. The placement arithmetic is tested
# directly because the Tk part needs a display. A parent that has not been
# laid out reports width 1 and rootx 0, which must not be trusted.
from adcp_launcher.app import centre_on                         # noqa: E402

SW, SH = 1920, 1080
DW, DH = 300, 140

x, y = centre_on(720, 600, 165, 170, DW, DH, SW, SH)
check(x == 165 + (720 - DW) // 2,
      "a dialog over a real parent is not horizontally centred on it: x=%d"
      % x)
check(170 < y < 170 + 600,
      "the dialog is not inside the parent vertically: y=%d" % y)

# A parent that has not been laid out reports width 1 and rootx 0. Trusting
# it gives 0 + (1 - 300)//2 = -150, clamped to 0: the top-left corner of the
# screen.
x, y = centre_on(1, 1, 0, 0, DW, DH, SW, SH)
check((x, y) != (0, 0),
      "a parent reporting 1x1 at (0,0) still puts the dialog in the "
      "top-left corner. That is the bug, not the fix")
check(abs(x - (SW - DW) // 2) <= 1,
      "the fallback is not the centre of the screen: x=%d" % x)

# a parent smaller than the dialog cannot be a real parent either
x, y = centre_on(200, 100, 900, 800, DW, DH, SW, SH)
check(abs(x - (SW - DW) // 2) <= 1,
      "a parent smaller than the dialog was trusted: x=%d" % x)

# and it never lands where the title bar cannot be reached
for pw, ph, px, py in ((720, 600, 1800, 1000), (720, 600, -400, -300),
                       (1, 1, 0, 0), (720, 600, 165, 170)):
    x, y = centre_on(pw, ph, px, py, DW, DH, SW, SH)
    check(0 <= x <= SW - DW,
          "x=%d is off the screen for parent %dx%d at (%d,%d)"
          % (x, pw, ph, px, py))
    check(0 <= y <= SH - DH,
          "y=%d is off the screen for parent %dx%d at (%d,%d)"
          % (y, pw, ph, px, py))

# a dialog bigger than the screen is pinned to the origin rather than given
# a negative coordinate
x, y = centre_on(1, 1, 0, 0, 2400, 140, SW, SH)
check(x == 0, "an over-wide dialog was given x=%d" % x)

psrc = inspect.getsource(appmod.place_over)
check("update_idletasks" in psrc,
      "place_over measures before the layout has run. Tk reports width 1 "
      "until then")
check("parent.update_idletasks" in psrc,
      "place_over updates the dialog but not the PARENT, and the parent is "
      "as capable of reporting width 1 as the dialog is, which is the "
      "measurement that was lying")
check("deiconify" in psrc,
      "the dialog is not shown by place_over, so it must have been mapped "
      "before it had a position, and on Windows a mapped Toplevel keeps "
      "where the window manager first put it")
check("except Exception" in psrc and "log" in psrc,
      "a failure to place the dialog is still swallowed silently; two "
      "rounds went on a dialog that was not being positioned at all and "
      "said nothing")

asrc = inspect.getsource(appmod.App.ask)
check("withdraw()" in asrc,
      "the dialog is created visible. It is mapped the moment it exists, "
      "and geometry set afterwards does not move it on Windows")
check(asrc.index("withdraw()") < asrc.index("place_over"),
      "the dialog is hidden after it is positioned, which is backwards")
check("place_over" in asrc, "the sign-in dialog is not placed")
check("focus_force" in asrc,
      "the sign-in dialog does not take the keyboard, so a password can be "
      "typed into whatever was focused instead")
check(asrc.index("place_over") < asrc.index("grab_set"),
      "the grab is taken before the dialog is placed and shown")

# --- a settings file from an older version must not outrank the site table --
check("SCHEMA" in dir(appmod),
      "there is no settings schema version, so a remembered path written "
      "before the site table existed goes on winning forever")
check('cfg.get("schema") == SCHEMA' in src,
      "the saved settings are trusted without checking which version wrote "
      "them, which is how the Tool directory stayed at ~/ADCP_docking on a "
      "cluster whose checkout is in a shared allocation")

# --- where there is a display, build the real window ------------------
if not os.environ.get("DISPLAY") and sys.platform not in ("win32", "darwin"):
    skipped = 1
    print("  (no display; the window is NOT built here. It is built on the "
          "Windows and macOS build runners before any bundle is uploaded.)")
else:
    import tkinter as tk
    root = tk.Tk()
    root.withdraw()
    a = appmod.App(root)
    root.update_idletasks()
    check(a.site.get() in clusters.sites(),
          "the cluster box opened on %r, which is not a site" % a.site.get())
    check(a.host.get() == clusters.SITES[a.site.get()]["host"],
          "the address was not filled from the site table")
    check(a.install.get() == clusters.SITES[a.site.get()]["install"],
          "the tool directory was not filled in: %r" % a.install.get())
    # typing the account fills the runs directory
    a.user.set("someone")
    root.update_idletasks()
    check(a.user_label.cget("text")
          == clusters.SITES[a.site.get()].get("label", "Username"),
          "the account label did not come from the site table")
    # a box the person edits is never overwritten by a later site change
    a.install.set("/my/own/checkout")
    a.site.set([s for s in clusters.sites() if s != a.site.get()
                and s in clusters.SITES][0])
    a._on_site()
    root.update_idletasks()
    check(a.install.get() == "/my/own/checkout",
          "changing the cluster overwrote a path the person had typed")
    # and a real dialog lands on the launcher, where a display exists to
    # check it.
    root.deiconify()
    root.geometry("700x560+180+150")
    root.update_idletasks()
    dlg = tk.Toplevel(root)
    dlg.withdraw()
    tk.Label(dlg, text="Password:").pack(padx=12, pady=8)
    tk.Entry(dlg, width=40).pack(padx=12, pady=8)
    appmod.place_over(dlg, root)
    root.update_idletasks()
    dlg.update_idletasks()
    dx, dy = dlg.winfo_rootx(), dlg.winfo_rooty()
    dw, dh = dlg.winfo_width(), dlg.winfo_height()
    px, py = root.winfo_rootx(), root.winfo_rooty()
    pw, ph = root.winfo_width(), root.winfo_height()
    check((dx, dy) != (0, 0),
          "the dialog opened at the top-left corner of the screen, which is "
          "the bug this file exists for")
    check(px - 40 <= dx and dx + dw <= px + pw + 40,
          "the dialog is at x=%d..%d; the launcher spans %d..%d"
          % (dx, dx + dw, px, px + pw))
    check(py - 40 <= dy and dy <= py + ph,
          "the dialog is at y=%d; the launcher spans %d..%d"
          % (dy, py, py + ph))
    dlg.destroy()
    root.destroy()

if fails:
    print("test_window FAILED (%d)" % len(fails))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  the launcher has its own mark, its own settings and its own port "
      "band, %d checks%s" % (ran, "; the window itself was skipped"
                             if skipped else ""))
