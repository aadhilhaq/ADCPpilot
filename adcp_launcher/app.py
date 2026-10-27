"""The window. Tkinter, because it is in the standard library and a launcher
that needs a toolkit installed before it can ask for a password is not a
launcher.

By Aadhil Haq

The shape is deliberately one screen: where to connect, who you are, where
the tool lives, and a button. Everything else is a log, because the useful
thing when a connection fails is what the far end said, not a spinner.
"""
import json
import os
import platform
import queue
import threading
import webbrowser

import tkinter as tk
from tkinter import messagebox, ttk

from . import (APP, BLURB, LAB, PI, AUTHORS, __version__, credit,
               HOW_TO_CITE)
from . import clusters
from . import session as sessionmod, tunnel

try:
    from adcp_dock import brand
except ImportError:                                            # noqa: BLE001
    brand = None

# This launcher's own settings file, not shared with any sibling launcher:
# several may be open at once, and one settings file between them would have
# each overwriting the others' cluster and username.
CONFIG = os.path.join(os.path.expanduser("~"), ".adcp_docking_launcher.json")


#: The banner's palette; the same four as the README banner, the page and
#: the icon. Duplicated from adcp_dock.brand because the window must still
#: draw if that import fails.
BANNER = dict(deep="#113f4a", deeper="#0b2e37", teal="#1a8a8a",
              warm="#d67a3a", paper="#f7f9f9", muted="#9fc4c4",
              faint="#6f9797")


def banner_layout(width, height, line_heights, line_widths, panel=80,
                  gutter=18, gap=20, line_gap=3, motif_min=140,
                  motif_clear=44):
    """Where everything in the header banner goes, as plain arithmetic.

    Kept apart from Tk for the same reason centre_on is: the drawing cannot
    run without a display, and the arithmetic is where alignment goes
    wrong. What it guarantees:

    - every text line starts at the same x, so the left edges align;
    - the text block is centered vertically as a BLOCK, on the same center
      line as the mark's panel, rather than each line being placed by eye;
    - the peptide motif starts after the widest line plus a clear margin,
      and is left out entirely when there is not room for it; a backbone
      drawn through the tagline is worse than none.
    """
    px = gutter
    py = (height - panel) / 2.0
    tx = px + panel + gap
    block = sum(line_heights) + line_gap * (len(line_heights) - 1)
    y = (height - block) / 2.0
    ys = []
    for lh in line_heights:
        ys.append(y)
        y += lh + line_gap
    start = tx + max(line_widths or [0]) + motif_clear
    stop = width - gutter
    motif = (start, stop) if stop - start >= motif_min else None
    return {"panel": (px, py, px + panel, py + panel),
            "text_x": tx, "text_ys": ys, "motif": motif}


def centre_on(pw, ph, px, py, w, h, sw, sh):
    """Where a w x h dialog goes over a pw x ph parent at (px, py).

    A pure function, because the Tk part cannot be exercised without a display
    and the arithmetic is where placement goes wrong. Returns a position on
    screen, clamped so the dialog is never placed off it.
    """
    if pw > w and ph > h:
        x = px + (pw - w) // 2
        y = py + (ph - h) // 3
    else:
        x = (sw - w) // 2
        y = (sh - h) // 3
    # Never off the top or left, and never so far right or down that the
    # title bar leaves the screen: a dialog you cannot drag back is worse
    # than one in the corner.
    x = max(0, min(x, max(0, sw - w)))
    y = max(0, min(y, max(0, sh - h)))
    return x, y


def place_over(win, parent, log=None):
    """Put a dialog over the window that raised it, and show it there.

    The window is created withdrawn and deiconified here, so it is never mapped
    at a position that was not chosen.
    """
    try:
        if parent is not None:
            parent.update_idletasks()
        win.update_idletasks()
        w = max(win.winfo_width(), win.winfo_reqwidth())
        h = max(win.winfo_height(), win.winfo_reqheight())
        pw = ph = px = py = 0
        if parent is not None and parent.winfo_viewable():
            pw, ph = parent.winfo_width(), parent.winfo_height()
            px, py = parent.winfo_rootx(), parent.winfo_rooty()
        x, y = centre_on(pw, ph, px, py, w, h,
                         win.winfo_screenwidth(), win.winfo_screenheight())
        # Size AND position, so the window manager has nothing left to
        # decide, and only then is it shown.
        win.geometry("%dx%d+%d+%d" % (w, h, x, y))
        win.deiconify()
        win.lift()
    except Exception as exc:                                    # noqa: BLE001
        if log:
            log("could not place the dialog (%s: %s); it will open wherever "
                "the window manager puts it" % (type(exc).__name__, exc))
        try:
            win.deiconify()
        except Exception:                                       # noqa: BLE001
            pass


#: Bumped when a default changes in a way that makes a REMEMBERED answer
#: wrong. A settings file from an older schema keeps the cluster and the
#: account, those are the person's, and lets the site table fill the paths
#: again. Without this, a default an older version wrote went on outranking
#: the table's, and the launcher kept pointing at a folder with nothing in
#: it.
SCHEMA = 2


def load():
    try:
        with open(CONFIG) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return {}


def save(d):
    try:
        with open(CONFIG, "w") as fh:
            json.dump(dict(d, schema=SCHEMA), fh, indent=1)
    except OSError:
        pass


class Banner(tk.Canvas):
    """The header, drawn as the README banner is: the mark on its panel, the
    name, the tagline and the build, and a peptide backbone to the right.

    Drawn, not pasted. A PNG banner scaled to the window blurs its text on a
    high-DPI screen and cannot reflow; canvas text is set in the system font
    at the system's size, stays sharp, and the layout is recomputed on every
    resize. The arithmetic lives in banner_layout, where it can be tested.
    """

    H = 108

    def __init__(self, parent, mark, name, tagline, build):
        super().__init__(parent, height=self.H, highlightthickness=0, bd=0,
                         bg=BANNER["deep"])
        self.mark, self.lines = mark, (name, tagline, build)
        import tkinter.font as tkfont
        fam = tkfont.nametofont("TkDefaultFont").actual()["family"]
        self.fonts = (tkfont.Font(family=fam, size=17, weight="bold"),
                      tkfont.Font(family=fam, size=10),
                      tkfont.Font(family=fam, size=9))
        self.bind("<Configure>", self._draw)

    @staticmethod
    def _mix(a, b, t):
        a = [int(a[i:i + 2], 16) for i in (1, 3, 5)]
        b = [int(b[i:i + 2], 16) for i in (1, 3, 5)]
        return "#%02x%02x%02x" % tuple(int(a[k] + (b[k] - a[k]) * t)
                                       for k in range(3))

    def _round(self, x1, y1, x2, y2, r, **kw):
        pts = [x1 + r, y1, x2 - r, y1, x2, y1, x2, y1 + r, x2, y2 - r,
               x2, y2, x2 - r, y2, x1 + r, y2, x1, y2, x1, y2 - r,
               x1, y1 + r, x1, y1]
        return self.create_polygon(pts, smooth=True, **kw)

    def _draw(self, _e=None):
        self.delete("all")
        w, h = max(self.winfo_width(), 1), self.H
        # The wash: darkest on the right, so the mark and the name sit on the
        # lighter side and keep their contrast.
        for x in range(0, w, 4):
            self.create_rectangle(x, 0, x + 4, h, width=0, fill=self._mix(
                BANNER["deep"], BANNER["deeper"], x / float(w)))
        widths = [f.measure(t) for f, t in zip(self.fonts, self.lines)]
        heights = [f.metrics("linespace") for f in self.fonts]
        L = banner_layout(w, h, heights, widths)

        if L["motif"]:
            a, b = L["motif"]
            n, mid = 7, h / 2.0
            step = (b - a) / float(n - 1)
            pts = [(a + i * step, mid + (13 if i % 2 else -13))
                   for i in range(n)]
            self.create_line(*[c for p in pts for c in p], width=4,
                             smooth=True, capstyle="round",
                             fill=self._mix(BANNER["deep"], BANNER["teal"], .42))
            for i, (x, y) in enumerate(pts):
                col = BANNER["warm"] if i % 2 == 0 else BANNER["teal"]
                self.create_oval(x - 8, y - 8, x + 8, y + 8, width=0,
                                 fill=self._mix(BANNER["deep"], col, .55))

        x1, y1, x2, y2 = L["panel"]
        self._round(x1, y1, x2, y2, 18, width=0,
                    fill=self._mix(BANNER["deep"], BANNER["teal"], .18))
        if self.mark is not None:
            self.create_image((x1 + x2) / 2.0, (y1 + y2) / 2.0,
                              image=self.mark)

        colours = (BANNER["paper"], BANNER["muted"], BANNER["faint"])
        for text, font, y, col in zip(self.lines, self.fonts, L["text_ys"],
                                      colours):
            self.create_text(L["text_x"], y, text=text, font=font, fill=col,
                             anchor="nw")


class App:
    def __init__(self, root):
        self.root = root
        self.q = queue.Queue()
        # Every line, plain and technical alike. The technical ones
        # are kept even while hidden, because they are wanted once a
        # connection has already failed.
        self.lines = []
        self.t = None
        self.session = None
        cfg = load()

        root.title("%s launcher" % (brand.NAME if brand else APP))
        root.minsize(680, 560)
        pad = dict(padx=10, pady=5)

        # The mark travels in adcp_dock.brand as bytes, so a bundled
        # launcher carries it without brand/ having to travel alongside. A
        # missing icon is never fatal: the window opening without one is far
        # better than a launcher that will not open because it cannot draw
        # itself.
        self.icon = None
        if brand is not None:
            try:
                self.icon = tk.PhotoImage(data=brand.logo_png())
                root.iconphoto(True, self.icon)
            except Exception:                                  # noqa: BLE001
                self.icon = None

        # THE BANNER, edge to edge above the form. The mark is the embedded
        # 128px logo halved: subsample, not a second asset, so there is one
        # source for the mark. Everything else in it is drawn.
        self.mark = None
        if self.icon is not None:
            try:
                self.mark = self.icon.subsample(2, 2)
            except Exception:                                  # noqa: BLE001
                self.mark = None
        self.banner = Banner(
            root, self.mark,
            brand.NAME if brand else APP,
            brand.TAGLINE if brand else BLURB,
            "build %s" % (brand.version() if brand else __version__))
        self.banner.pack(fill="x", pady=(0, 6))

        form = ttk.LabelFrame(root, text="Connection")
        form.pack(fill="x", **pad)

        # Selecting a site fills in what the table knows. On a cluster
        # the address and the username are typed once and remembered. A
        # saved site this version does not offer becomes the cluster,
        # with its address and path kept.
        self.site = tk.StringVar(value=clusters.site_of(
            cfg.get("site", clusters.sites()[0])))
        ttk.Label(form, text="Cluster").grid(row=0, column=0, sticky="e",
                                             padx=6, pady=4)
        self.site_box = ttk.Combobox(form, textvariable=self.site,
                                     values=clusters.sites(),
                                     state="readonly")
        self.site_box.grid(row=0, column=1, sticky="ew", padx=6, pady=4)
        self.site_box.bind("<<ComboboxSelected>>", self._on_site)

        self.host = self._row(form, 1, "Address", cfg.get("host", ""))
        self.user_label = ttk.Label(form, text="Username")
        self.user_label.grid(row=2, column=0, sticky="e", padx=6, pady=4)
        self.user = tk.StringVar(value=cfg.get("user", ""))
        e = ttk.Entry(form, textvariable=self.user)
        e.grid(row=2, column=1, sticky="ew", padx=6, pady=4)
        self.user.trace_add("write", lambda *a: self._on_user())

        self.install = self._row(form, 3, "Tool directory",
                                 cfg.get("install", ""))
        # No runs directory here, deliberately: a user signing in
        # from a laptop has no way to know which of a cluster's
        # scratch paths to use, and a field that must be filled
        # before anything works is one that will be filled wrongly.
        self.port = self._row(form, 4, "Port on the cluster",
                              cfg.get("port", "auto"))
        ttk.Label(form, text="auto lets the server pick one that is free, "
                             "which is what stops two people colliding",
                  foreground="#5f7070").grid(row=5, column=1, sticky="w",
                                             padx=6, pady=(0, 4))
        self.tools = self._row(form, 6, "ADCP 1.1 suite / Amber roots",
                               cfg.get("tools", ""))
        ttk.Label(form, text="comma separated; leave blank once activate.sh "
                             "exists on the cluster",
                  foreground="#5f7070").grid(row=7, column=1, sticky="w",
                                             padx=6, pady=(0, 4))
        self.note = ttk.Label(form, text="", foreground="#5f7070",
                              wraplength=560, justify="left")
        self.note.grid(row=8, column=1, sticky="w", padx=6, pady=(0, 6))
        # Which fields the user has taken over. A field they typed
        # into is never rewritten by a later site change.
        self._touched = set()
        for name in ("host", "install"):
            var = getattr(self, name)
            var.trace_add("write", self._mark(name, var))
        self._apply_site(first=True, cfg=cfg)

        form.columnconfigure(1, weight=1)

        bar = ttk.Frame(root)
        bar.pack(fill="x", **pad)
        self.go = ttk.Button(bar, text="Connect and open", command=self.start)
        self.go.pack(side="left")
        self.open_btn = ttk.Button(bar, text="Open in browser",
                                   command=self.open_browser, state="disabled")
        self.open_btn.pack(side="left", padx=6)
        self.stop_btn = ttk.Button(bar, text="Disconnect", command=self.close,
                                   state="disabled")
        self.stop_btn.pack(side="left")
        self.status = ttk.Label(bar, text="not connected", foreground="#5f7070")
        self.status.pack(side="right")

        # Plain progress lines, with the technical detail behind a
        # switch. The ssh command, the host key check and the server
        # output help whoever is debugging a connection; they are
        # kept, but hidden by default.
        box = ttk.LabelFrame(root, text="What is happening")
        box.pack(fill="both", expand=True, **pad)
        self.log = tk.Text(box, height=15, wrap="word",
                           font=("TkFixedFont", 10))
        self.log.pack(fill="both", expand=True, side="left")
        sb = ttk.Scrollbar(box, command=self.log.yview)
        sb.pack(fill="y", side="right")
        self.log.configure(yscrollcommand=sb.set, state="disabled")

        self.verbose = tk.BooleanVar(value=False)
        ttk.Checkbutton(root, text="Show technical details",
                        variable=self.verbose, command=self._replay
                        ).pack(anchor="w", padx=12)

        # The credit line along the bottom: the line a screenshot
        # carries with it.
        ttk.Label(root, text=credit(), foreground="#5f7070", wraplength=690,
                  justify="left").pack(anchor="w", padx=10, pady=(4, 8))

        root.protocol("WM_DELETE_WINDOW", self.quit)
        self.pump()

    # ---- the site table drives the form ---------------------------------
    def _mark(self, name, var):
        """Remember that a box was typed into, so nothing overwrites it."""
        def cb(*_a):
            if getattr(self, "_filling", False):
                return
            self._touched.add(name)
        return cb

    def _set(self, name, value):
        """Fill a box the person has not taken over."""
        if name in self._touched:
            return
        self._filling = True
        try:
            getattr(self, name).set(value)
        finally:
            self._filling = False

    def _apply_site(self, first=False, cfg=None):
        d = clusters.defaults_for(self.site.get(), self.user.get().strip())
        cfg = cfg or {}
        if first:
            # A saved answer counts as the user's own, so a
            # remembered path is not replaced by the table's,
            # but only if it was saved by this version, since an
            # older settings file holds a default that was never
            # a deliberate choice.
            if cfg.get("schema") == SCHEMA:
                for name in ("host", "install"):
                    if cfg.get(name):
                        self._touched.add(name)
        self._set("host", d["host"])
        self._set("install", d["install"])
        self.user_label.configure(text=d["label"])
        note = d["note"]
        if d["hint"]:
            note = "%s  %s" % (note, d["hint"])
        self.note.configure(text=note.strip())

    def _on_site(self, _event=None):
        # A different cluster means different paths, so the fields
        # the user has not edited are refilled. _apply_site already
        # leaves edited fields alone, so nothing else is needed here.
        self._apply_site()

    def _on_user(self):
        self._apply_site()

    def suggested_port(self):
        return clusters.suggested_port(self.user.get().strip())

    def _row(self, parent, r, label, value, combo=None):
        ttk.Label(parent, text=label).grid(row=r, column=0, sticky="e",
                                           padx=6, pady=4)
        var = tk.StringVar(value=value)
        if combo:
            w = ttk.Combobox(parent, textvariable=var, values=list(combo))
        else:
            w = ttk.Entry(parent, textvariable=var)
        w.grid(row=r, column=1, sticky="ew", padx=6, pady=4)
        return var

    # ---- the log and the queue ------------------------------------------
    def say(self, text):
        """A line anyone should read."""
        self.q.put(("log", str(text)))

    def detail(self, text):
        """Machinery. Shown only with the switch on, but always KEPT, so a
        failure can print it without asking the user to reconnect with the
        box ticked, which is the one moment they cannot reproduce."""
        for ln in str(text).splitlines():
            if ln.strip():
                self.q.put(("detail", ln.rstrip()))

    def from_cluster(self, text):
        """The server's own output. Two lines of it are worth reading."""
        for ln in str(text).splitlines():
            if not ln.strip():
                continue
            if "runs are kept in" in ln:
                self.say("  " + ln.strip())
            else:
                self.detail(ln)

    def _replay(self):
        """Redraw the log when the details switch is flipped."""
        self.log.configure(state="normal")
        self.log.delete("1.0", "end")
        for kind, line in self.lines:
            if kind == "detail" and not self.verbose.get():
                continue
            self.log.insert("end", ("    " + line if kind == "detail"
                                    else line) + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")

    def show_details(self):
        """Turn the switch on and redraw: used when something failed."""
        self.verbose.set(True)
        self._replay()

    def pump(self):
        try:
            while True:
                kind, payload = self.q.get_nowait()
                if kind in ("log", "detail"):
                    self.lines.append((kind, payload.rstrip()))
                    if kind == "detail" and not self.verbose.get():
                        continue
                    self.log.configure(state="normal")
                    self.log.insert(
                        "end", ("    " if kind == "detail" else "")
                        + payload.rstrip() + "\n")
                    self.log.see("end")
                    self.log.configure(state="disabled")
                elif kind == "status":
                    self.status.configure(text=payload)
                elif kind == "ready":
                    self.open_btn.configure(state="normal")
                    self.stop_btn.configure(state="normal")
                    self.go.configure(state="normal")
                    self.open_browser()
                elif kind == "details":
                    self.show_details()
                elif kind == "failed":
                    self.go.configure(state="normal")
                    messagebox.showerror(APP, payload)
        except queue.Empty:
            pass
        self.root.after(120, self.pump)

    # ---- the two-factor dialog, on the main thread -------------------------
    def ask(self, title, instructions, prompts):
        """One dialog per challenge, however many fields it carries.

        This is called from the worker thread, so the window is built on the
        main thread and the worker waits on an event: a Tk widget touched
        from another thread is a crash on some platforms and a hang on the
        rest.
        """
        done = threading.Event()
        answers = []

        def build():
            win = tk.Toplevel(self.root)
            # Hidden until it has a position. A Toplevel is
            # mapped as soon as it is created, and on Windows it
            # then stays where the window manager first put it,
            # the top-left corner of the screen, whatever
            # geometry is set afterwards.
            win.withdraw()
            win.title(title or "%s: sign in" % APP)
            win.transient(self.root)
            if instructions and instructions.strip():
                ttk.Label(win, text=instructions.strip(), wraplength=460,
                          justify="left").pack(anchor="w", padx=12, pady=(12, 4))
            entries = []
            for text, echo in prompts:
                ttk.Label(win, text=text.strip() or "Response",
                          wraplength=460, justify="left").pack(
                              anchor="w", padx=12, pady=(8, 2))
                e = ttk.Entry(win, width=44, show="" if echo else "•")
                e.pack(fill="x", padx=12)
                entries.append(e)
            if entries:
                entries[0].focus_set()
            row = ttk.Frame(win)
            row.pack(fill="x", padx=12, pady=12)

            def ok(*_):
                answers.extend(e.get() for e in entries)
                win.destroy()
                done.set()

            def cancel():
                win.destroy()
                done.set()

            ttk.Button(row, text="Continue", command=ok).pack(side="right")
            ttk.Button(row, text="Cancel", command=cancel).pack(side="right",
                                                                padx=6)
            win.bind("<Return>", ok)
            win.bind("<Escape>", lambda *_a: cancel())
            win.protocol("WM_DELETE_WINDOW", cancel)
            # Over the launcher, sized, and holding the keyboard.
            # A dialog in the far corner of the screen is easy to
            # type past, and a password typed into another window
            # is lost.
            win.resizable(False, False)
            place_over(win, self.root, log=self.say)
            win.grab_set()
            win.attributes("-topmost", True)
            win.after(300, lambda: win.attributes("-topmost", False))
            if entries:
                entries[0].focus_force()

        self.root.after(0, build)
        done.wait()
        return answers or [""] * len(prompts)

    # ---- the work -------------------------------------------------------
    def start(self):
        host, user = self.host.get().strip(), self.user.get().strip()
        if self.site.get() == clusters.THIS:
            self.go.configure(state="disabled")
            self.q.put(("status", "starting here"))
            threading.Thread(target=self.work_here, daemon=True).start()
            return
        if not host or not user:
            messagebox.showwarning(APP, "The cluster and your username are both "
                                        "needed.")
            return
        save({"site": self.site.get(), "port": self.port.get(),
              "host": host, "user": user, "install": self.install.get(),
              "tools": self.tools.get()})
        self.go.configure(state="disabled")
        self.q.put(("status", "connecting"))
        threading.Thread(target=self.work, args=(host, user), daemon=True
                         ).start()

    def work_here(self):
        """This computer: start `adcp-dock local`, read its address, open it."""
        from . import localrun
        try:
            self.say("Starting on this computer...")
            if platform.system() == "Darwin":
                self.say("The first start installs ADFRsuite 1.0 from the "
                         "Scripps Research Institute, which takes a few "
                         "minutes. On a Mac the docking is ADCP 1.0: Scripps "
                         "builds no ADCP 1.1 for macOS, and its scores are not "
                         "comparable with ADCP 1.1 runs.")
            else:
                self.say("The first start installs the ADCP 1.1 suite with "
                         "the Scripps Research Institute's own recipe, which "
                         "takes a few minutes and about 3 GB.")
            self.session = localrun.start(on_output=self.detail)
            self.say("Ready. Opening it in your browser.")
            self.say("  If it does not open, paste this in: %s"
                     % self.session.url)
            self.say("")
            self.say(HOW_TO_CITE)
            self.q.put(("status", "running here"))
            self.q.put(("ready", ""))
        except Exception as exc:                                # noqa: BLE001
            self.say("")
            self.say(str(exc))
            self.q.put(("details", ""))
            self.q.put(("status", "not running"))
            self.q.put(("failed", str(exc).split("\n\n")[0]))

    def work(self, host, user):
        try:
            self.say("Connecting to %s as %s..." % (host, user))
            self.say("Complete your cluster's second sign-in step, if it "
                     "has one, when it is asked for.")
            self.t = tunnel.connect(host, user, self.ask, log=self.detail)
            self.say("Signed in.")
            self.q.put(("status", "starting the tool"))
            self.say("Starting the tool on the login node...")
            roots = [r.strip() for r in self.tools.get().split(",") if r.strip()]
            port = clusters.parse_port(self.port.get())
            if port:
                self.detail("asking the cluster for port %d" % port)
            # A facility's acknowledgment is not set here: the server
            # reads ADCP_ACKNOWLEDGEMENT from activate.sh on the cluster,
            # which the command sources first. A named cluster's table only
            # suggests where runs go and the Slurm account its jobs need,
            # where activate.sh has not said.
            suggest = clusters.server_env(self.site.get(), user)
            if suggest:
                self.detail("suggested to the server: %s" % ", ".join(
                    "%s=%s" % kv for kv in sorted(suggest.items())))
            self.session = sessionmod.start(
                self.t, self.install.get().strip(), "",
                port=port, tools_roots=roots, on_output=self.from_cluster,
                forward=tunnel.forward, free_port=tunnel.free_port,
                defaults=suggest)
            self.say("Ready. Opening it in your browser.")
            # The link stays on screen because "Open in browser" opens the
            # DEFAULT browser, and on a managed laptop that is not always the
            # one the person is signed in to.
            self.say("  If it does not open, paste this in: %s"
                     % self.session.url)
            self.say("")
            # Not the credit sentence: the label along the bottom
            # of the window already shows it, where it stays
            # visible instead of scrolling away.
            self.say(HOW_TO_CITE)
            self.q.put(("status", "connected"))
            self.q.put(("ready", ""))
        except Exception as exc:                                # noqa: BLE001
            self.say("")
            self.say("Could not connect: %s" % exc)
            # The detail is shown without being asked for:
            # telling a user to enable it and try again asks
            # them to reproduce a failure.
            self.say("The technical detail is below; send it with a report.")
            self.q.put(("details", ""))
            self.q.put(("status", "not connected"))
            self.q.put(("failed", str(exc)))

    def open_browser(self):
        if self.session:
            webbrowser.open(self.session.url)

    def close(self):
        if self.session:
            self.session.stop()
            self.session = None
        if self.t:
            try:
                self.t.close()
            except Exception:                                   # noqa: BLE001
                pass
            self.t = None
        self.open_btn.configure(state="disabled")
        self.stop_btn.configure(state="disabled")
        self.q.put(("status", "not connected"))
        self.say("Disconnected.")

    def quit(self):
        self.close()
        self.root.destroy()


def main():
    root = tk.Tk()
    try:
        ttk.Style().theme_use("clam")
    except Exception:                                           # noqa: BLE001
        pass
    App(root)
    root.mainloop()
