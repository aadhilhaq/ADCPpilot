#!/usr/bin/env python3
"""Screenshots for the documentation, taken of the running interface.

By Aadhil Haq

    python tools/docs_screenshots.py --python /path/to/env/bin/python \\
        [--tools-root DIR] [--out adcp_dock/docs/img]

Starts the interface with a fresh runs folder holding the shipped example
run, drives the page the way a person would (loads the example's receptor,
sets its site, pastes its peptides, opens a help box, opens the finished run,
opens ADapt when it is installed) and saves each view the docs show. Run it
again after the page changes, so the pictures never describe an older page.

Each view is taken inside its own step: one that fails is reported by name at
the end, with the others still taken, because a missing picture in a page is
easier to act on than a traceback from the first click that moved.

Needs Playwright with Chromium in the Python running this script; the
interface itself runs in the Python given with --python. Run it on a compute
node: it starts a server, a browser and possibly ADapt.
"""
import argparse
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import traceback

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLE = os.path.join(HERE, "adcp_dock", "example_run")

#: Every picture the pages use. Reported as missing when a step did not
#: produce it, so the docs and this script cannot drift apart silently.
WANT = ["interface.png", "receptor-loaded.png", "site-residues.png",
        "site-box.png", "site-other-chain.png", "peptides.png",
        "peptides-refused.png", "docking-settings.png", "help-box.png",
        "runs.png", "results-table.png", "pose-viewer.png", "cite.png",
        "adapt.png"]


def start_server(python, runs, roots):
    cmd = [python, "-u", "-m", "adcp_dock.cli", "serve", "--host", "127.0.0.1",
           "--port", "0", "--runs", runs]
    for r in roots:
        cmd += ["--tools-root", r]
    srv = subprocess.Popen(cmd, cwd=HERE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                           text=True, start_new_session=True)
    t0 = time.time()
    while time.time() - t0 < 180:
        line = srv.stdout.readline()
        if not line:
            break
        m = re.search(r"(http://127\.0\.0\.1:\d+/\?t=\S+)", line)
        if m:
            return srv, m.group(1)
    os.killpg(srv.pid, signal.SIGTERM)
    raise SystemExit("the interface did not start")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--python", required=True, help="the Python the interface runs in")
    ap.add_argument("--tools-root", action="append", default=[])
    ap.add_argument("--out", default=os.path.join(HERE, "adcp_dock", "docs", "img"))
    a = ap.parse_args()
    from playwright.sync_api import sync_playwright

    job = json.load(open(os.path.join(EXAMPLE, "job.json")))
    work = tempfile.mkdtemp(prefix="docshots_")
    runs = os.path.join(work, "runs")
    subprocess.run([a.python, "-m", "adcp_dock.cli", "example", "--runs", runs],
                   cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    receptor = os.path.join(work, "1YCR.pdb")
    shutil.copy(os.path.join(EXAMPLE, "receptor", "1YCR.pdb"), receptor)
    os.makedirs(a.out, exist_ok=True)
    srv, url = start_server(a.python, runs, a.tools_root)
    print("interface at", url, flush=True)
    made, errors, failed = [], [], []

    def shot(target, name):
        path = os.path.join(a.out, name)
        target.screenshot(path=path)
        made.append(name)
        print("  saved %-24s %6d bytes" % (name, os.path.getsize(path)), flush=True)

    def shot_clip(selectors, name, pad=10, least=0):
        """The region covering these elements, padded, and never narrower
        than `least`: a picture of one narrow panel is unreadable once a page
        scales it to the width of a figure."""
        boxes = [pg.locator(s).bounding_box() for s in selectors]
        boxes = [b for b in boxes if b]
        if not boxes:
            raise RuntimeError("none of %s is on the page" % (selectors,))
        x0 = min(b["x"] for b in boxes) - pad
        y0 = min(b["y"] for b in boxes) - pad
        x1 = max(b["x"] + b["width"] for b in boxes) + pad
        y1 = max(b["y"] + b["height"] for b in boxes) + pad
        if least and x1 - x0 < least:
            x0 = max(0.0, x1 - least)
        vw = pg.viewport_size["width"]
        vh = pg.viewport_size["height"]
        clip = {"x": max(0.0, x0), "y": max(0.0, y0),
                "width": min(x1, vw) - max(0.0, x0),
                "height": min(y1, vh) - max(0.0, y0)}
        path = os.path.join(a.out, name)
        pg.screenshot(path=path, clip=clip)
        made.append(name)
        print("  saved %-24s %6d bytes" % (name, os.path.getsize(path)), flush=True)

    def step(label, fn):
        """One view. A step that raises is named and the rest still run."""
        try:
            fn()
        except Exception as exc:                                # noqa: BLE001
            failed.append("%s: %s" % (label, str(exc).splitlines()[0][:160]))
            print("  FAILED %s: %s" % (label, exc), flush=True)
            traceback.print_exc(limit=1)

    try:
        with sync_playwright() as pw:
            br = pw.chromium.launch(args=["--use-gl=swiftshader", "--enable-unsafe-swiftshader"])
            ctx = br.new_context(viewport={"width": 1360, "height": 900}, color_scheme="light")
            ctx.add_init_script("try{localStorage.setItem('adcp.theme','light')}catch(e){}")
            pg = ctx.new_page()
            pg.on("pageerror", lambda e: errors.append(str(e)[:200]))
            pg.goto(url, wait_until="load", timeout=120000)
            pg.wait_for_timeout(2500)
            card = lambda sel: pg.locator(sel).locator(
                "xpath=ancestor-or-self::div[contains(@class,'card')][1]")

            # the page as it opens
            def opening():
                pg.evaluate("window.scrollTo(0, 0)")
                shot(pg, "interface.png")
            step("interface", opening)

            # the receptor, loaded from a file on this machine
            def load_receptor():
                pg.fill("#recpath", receptor)
                pg.click("#recload")
                pg.wait_for_function(
                    "() => /atoms/.test(document.getElementById('recinfo').textContent)",
                    timeout=120000)
                for ch in (job.get("chains") or []):
                    pg.check('#chainbox .chk[data-ch="%s"]' % ch)
                pg.wait_for_timeout(1500)
                shot(card("#reccode"), "receptor-loaded.png")
            step("receptor", load_receptor)

            # the site: the example's residues, the box, and what else it reaches
            def site():
                pg.fill("#hotspots", job["hotspots"])
                pg.fill("#pad", str(job.get("box_padding", 8)))
                pg.fill("#bmin", str(job.get("box_min", 20)))
                pg.fill("#bmax", str(job.get("box_max", 30)))
                pg.click("#sitecheck")
                pg.wait_for_function(
                    "() => /box/.test(document.getElementById('sitenote').textContent)",
                    timeout=60000)
                pg.wait_for_timeout(2000)
                shot(card("#hotspots"), "site-residues.png")
            step("site", site)

            def box():
                pg.check("#showbox")
                pg.locator("#sitepick").scroll_into_view_if_needed()
                # The viewer draws the structure, then the cage, in software.
                pg.wait_for_timeout(4000)
                # Both ways of picking a residue, side by side: the viewer
                # alone is narrower than the pages show a figure.
                shot(pg.locator("#sitepick"), "site-box.png")
            step("box", box)

            def other_chain():
                pg.wait_for_selector("#otherchains:not(.hide)", timeout=30000)
                # The notice, with the box summary above it, so the picture
                # says which box reached that chain.
                shot_clip(["#sitenote", "#otherchains"], "site-other-chain.png")
            step("other chain", other_chain)

            # the peptides, pasted, and one line the page refuses
            def peptides():
                pg.click("#pastepep")
                text = "\n".join("%s: %s" % (p["name"], p["sequence"])
                                 for p in job["peptides"])
                pg.fill("#pastetext", text)
                pg.click("#pasteadd")
                pg.wait_for_timeout(1500)
                shot(card("#peps"), "peptides.png")
            step("peptides", peptides)

            def refused():
                # Two lines the page will not add: a name written with a space,
                # which would otherwise be read as one longer sequence, and a
                # peptide shorter than ADCP can build. Nothing is added, the
                # text stays in the box, and the message line names each line
                # and its problem. The verdict comes back from the server, so
                # the wait is on the message.
                pg.click("#pastepep")
                pg.fill("#pastetext", "alpha SWQWWWQPQRDCR\ntiny: AAA")
                pg.click("#pasteadd")
                pg.wait_for_function(
                    "() => /nothing added/.test("
                    "document.getElementById('pepmsg').textContent)", timeout=30000)
                pg.wait_for_timeout(800)
                shot(card("#peps"), "peptides-refused.png")
            step("refused", refused)

            # the docking and rescoring settings, and one help box
            def settings():
                card("#poses").scroll_into_view_if_needed()
                pg.wait_for_timeout(500)
                shot(card("#poses"), "docking-settings.png")
            step("settings", settings)

            def helpbox():
                # The "?" buttons are added by the page to each label that has
                # an entry, so the one for a setting is inside its own label.
                pg.click('label[for="replicas"] button.q')
                pg.wait_for_selector("#helpbox:not(.hide)", timeout=30000)
                pg.wait_for_timeout(600)
                # The panel with the page behind it: on its own it is 420 px
                # wide, which is narrower than a figure here is shown.
                shot_clip(["#helpbox"], "help-box.png", pad=24, least=760)
                pg.click("#helpclose")
            step("help box", helpbox)

            # the runs, and the finished example opened
            def runs_card():
                c = card("#runs")
                c.scroll_into_view_if_needed()
                pg.wait_for_selector("#runs table", timeout=60000)
                pg.wait_for_timeout(1500)
                shot(c, "runs.png")
            step("runs", runs_card)

            def results():
                pg.click("#runs button[data-open]")
                # The table is drawn once the report has been read back.
                pg.wait_for_selector("#rundetail table tr.pick", timeout=90000)
                pg.wait_for_timeout(2500)
                pg.locator("#rundetail").scroll_into_view_if_needed()
                shot(pg.locator("#rundetail"), "results-table.png")
            step("results", results)

            def viewer():
                # The viewer is hidden until a peptide is chosen, which is what
                # the view button on a row does.
                pg.locator("#rundetail button[data-view]").first.click()
                pg.wait_for_selector("#viewcard:not(.hide)", timeout=60000)
                pg.locator("#viewcard").scroll_into_view_if_needed()
                # The pose, its contacts and the residue map, drawn in software.
                pg.wait_for_timeout(8000)
                shot(pg.locator("#viewcard"), "pose-viewer.png")
            step("pose viewer", viewer)

            # how to cite
            def cite():
                pg.evaluate("window.scrollTo(0, 0)")
                pg.click("#citebtn")
                pg.wait_for_selector("dialog[open] #citetext", timeout=60000)
                pg.wait_for_timeout(1200)
                shot(pg.locator("dialog[open]").first, "cite.png")
                pg.keyboard.press("Escape")
            step("cite", cite)

            # ADapt, through the interface, when it is installed here
            def adapt():
                if not pg.locator("#adaptopen").is_visible():
                    print("  ADapt is not offered here, so no adapt.png",
                          flush=True)
                    return
                with ctx.expect_page(timeout=60000) as popup:
                    pg.click("#adaptopen")
                tab = popup.value
                tab.wait_for_url("**/adapt/?t=*", timeout=330000)
                tab.wait_for_load_state("load")
                tab.set_input_files("#file", receptor)
                tab.wait_for_function(
                    "() => !document.body.innerText.includes('No structure loaded')",
                    timeout=240000)
                tab.wait_for_timeout(5000)
                shot(tab, "adapt.png")
            step("adapt", adapt)
            br.close()
    finally:
        try:
            os.killpg(srv.pid, signal.SIGTERM)
        except OSError:
            pass
        shutil.rmtree(work, ignore_errors=True)
    print("\nmade:", " ".join(made))
    print("page errors:", errors[:5] or "none")
    print("steps that failed:", failed or "none")
    missing = [n for n in WANT if n not in made]
    print("missing:", missing or "none")
    return 1 if missing or errors or failed else 0


if __name__ == "__main__":
    sys.exit(main())
