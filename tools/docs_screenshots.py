#!/usr/bin/env python3
"""Screenshots for the documentation, taken of the running interface.

By Aadhil Haq

    python tools/docs_screenshots.py --python /path/to/env/bin/python \\
        [--tools-root DIR] [--out adcp_dock/docs/img] [--work DIR]

The pictures show the folders the run is in, so --work chooses them: the
folder named there holds the receptor and an `adcp_runs` beside it, and its
path is what appears in the Browse and Output folder boxes.

Starts the interface with a fresh runs folder holding the shipped example
run, drives the page the way a person would (loads the example's receptor,
sets its site, pastes its peptides, opens a help box, opens the finished run,
and, when ADapt is installed, shows it starting and then in the Receptor
card, opens it full screen, unticks the bound peptide in it and has ADapt write
the receptor)
and saves each view the docs show. Run it
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
        "site-box.png", "site-other-chain.png", "site-ligand.png",
        "site-maps.png", "peptides.png", "peptides-refused.png",
        "peptide-rings.png", "peptide-ring-refused.png",
        "peptide-options.png", "docking-settings.png", "help-box.png",
        "runs.png", "results-table.png", "pose-viewer.png", "cite.png",
        "adapt.png", "guide.png", "submit.png", "runs-folders.png",
        "file-viewer.png", "site-surface.png", "receptor-adapt.png",
        "adapt-inline.png", "adapt-loading.png", "add-measures.png"]


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
    ap.add_argument("--work", default="",
                    help="the folder the receptor and the runs folder are "
                         "made in, whose path the pictures show. The default "
                         "is a temporary folder, removed afterwards; a folder "
                         "named here is left as it is.")
    a = ap.parse_args()
    from playwright.sync_api import sync_playwright

    job = json.load(open(os.path.join(EXAMPLE, "job.json")))
    work = a.work or tempfile.mkdtemp(prefix="docshots_")
    temporary = not a.work
    os.makedirs(work, exist_ok=True)
    runs = os.path.join(work, "adcp_runs")
    subprocess.run([a.python, "-m", "adcp_dock.cli", "example", "--runs", runs],
                   cwd=HERE, check=True, stdout=subprocess.DEVNULL)
    receptor = os.path.join(work, "1YCR.pdb")
    shutil.copy(os.path.join(EXAMPLE, "receptor", "1YCR.pdb"), receptor)
    os.makedirs(a.out, exist_ok=True)
    srv, url = start_server(a.python, runs, a.tools_root)
    print("interface at", url, flush=True)
    made, errors, failed = [], [], []

    def bar(shown):
        """The bar along the bottom is fixed to the window, so in a picture
        of one part of the page it would lie across that part's foot. It is
        hidden for those, and kept in pictures of the whole window."""
        pg.evaluate("document.querySelectorAll('.creditbar').forEach("
                    "b => b.style.visibility = %s)" % ("''" if shown else "'hidden'"))

    def shot(target, name):
        path = os.path.join(a.out, name)
        whole = target is pg
        if not whole:
            bar(False)
        try:
            target.screenshot(path=path)
        finally:
            if not whole:
                bar(True)
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
        bar(False)
        try:
            pg.screenshot(path=path, clip=clip)
        finally:
            bar(True)
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
            # The steps across the top follow the page down. In a picture of
            # one card that would lay the strip over the card's top, so the
            # strip is held at the top of the page for every picture but its
            # own.
            sticky = lambda on: pg.evaluate(
                "document.getElementById('guide').style.position = %s"
                % ("''" if on else "'static'"))
            sticky(False)

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
                # ADapt starting on it, where ADapt is installed: the loading
                # card in the space its frame will fill, taken while one step
                # is ticked and the next is under way. A start already warm
                # can take only seconds, so it is taken the moment that holds.
                if pg.locator("#adaptframewrap .ld-stage").count():
                    pg.evaluate("document.getElementById('adaptsec')"
                                ".scrollIntoView({block: 'start'})")
                    pg.wait_for_function(
                        "() => !!document.querySelector('#adaptframewrap .ld-step.ld-ok')"
                        " && !!document.querySelector('#adaptframewrap .ld-step.ld-now')",
                        timeout=120000, polling="raf")
                    shot(pg.locator("#adaptsec"), "adapt-loading.png")
                pg.wait_for_timeout(1500)
                # ADapt opens in the card on the structure by itself, where it
                # is installed. Folded for this picture, which is about the
                # chains and the preparation; its own pictures come later.
                if pg.locator("#adaptpanel").is_visible():
                    pg.evaluate("adaptFold(true)")
                    pg.wait_for_timeout(500)
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

            # the receptor's surface, from MSMS, over the same view
            def site_surface():
                pg.check("#sitesurf")
                pg.wait_for_function(
                    "() => !document.querySelector('#surfnote .ld-inline')",
                    timeout=180000)
                pg.locator("#sitepick").scroll_into_view_if_needed()
                pg.wait_for_timeout(5000)
                shot(pg.locator("#sitepick"), "site-surface.png")
                pg.uncheck("#sitesurf")
                pg.wait_for_timeout(1500)
            step("site surface", site_surface)

            def other_chain():
                pg.wait_for_selector("#otherchains:not(.hide)", timeout=30000)
                pg.evaluate("document.getElementById('sitenote')"
                            ".scrollIntoView({block: 'center'})")
                pg.wait_for_timeout(400)
                # The notice, with the box summary above it, so the picture
                # says which box reached that chain.
                shot_clip(["#sitenote", "#otherchains"], "site-other-chain.png")
            step("other chain", other_chain)

            # the box around the crystal peptide instead, then back again
            def site_ligand():
                pg.select_option("#site_mode", "ligand")
                pg.wait_for_function(
                    "() => document.getElementById('site_ligand').value === 'B'",
                    timeout=30000)
                pg.wait_for_function(
                    "() => /touch it/.test("
                    "document.getElementById('sitenote').textContent)",
                    timeout=60000)
                card("#site_mode").scroll_into_view_if_needed()
                # The viewer redraws the peptide in sticks and the new box.
                pg.wait_for_timeout(4000)
                shot(card("#site_mode"), "site-ligand.png")
            step("site ligand", site_ligand)

            def site_maps():
                pg.locator("#mapsgroup").scroll_into_view_if_needed()
                pg.wait_for_timeout(500)
                shot_clip(["#mapsgroup"], "site-maps.png", pad=14, least=760)
            step("maps", site_maps)

            def residues_again():
                pg.select_option("#site_mode", "residues")
                pg.wait_for_function(
                    "() => /box/.test(document.getElementById('sitenote').textContent)",
                    timeout=60000)
                pg.wait_for_timeout(1500)
            step("residues again", residues_again)

            # the steps across the top: the receptor and the site done, the
            # peptides to do now, with the Peptides card lit under the strip
            def guide():
                pg.wait_for_function(
                    "() => document.querySelector('#guide li[data-step=site]')"
                    ".classList.contains('done') && document.querySelector("
                    "'#guide li[data-step=peptides]').classList.contains('now')",
                    timeout=30000)
                sticky(True)
                pg.evaluate("document.getElementById('step-peptides')"
                            ".scrollIntoView({block: 'start'})")
                pg.wait_for_timeout(800)
                shot_clip(["#guide", "#step-peptides .hint"], "guide.png", pad=10)
                sticky(False)
            step("guide", guide)

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

            # one run, three kinds of peptide: a fourth peptide pasted with
            # its ring as a third column, and p53 set head to tail on its row
            def rings():
                # The paste box is still open from the refused paste.
                pg.evaluate("document.getElementById('pastebox')"
                            ".classList.remove('hide')")
                pg.fill("#pastetext", "ss_PMI,CTSFAEYWNLLC,disulfide")
                pg.click("#pasteadd")
                pg.wait_for_function(
                    "() => document.querySelectorAll('#peps .ent').length >= 4",
                    timeout=30000)
                pg.select_option("#peps .ent[data-i='0'] select.pr", "backbone")
                pg.wait_for_function(
                    "() => /rows differ/.test("
                    "document.getElementById('pepoptsnote').textContent)",
                    timeout=30000)
                pg.wait_for_timeout(1200)
                shot(card("#peps"), "peptide-rings.png")
            step("rings", rings)

            # a ring the peptide cannot close: 12/1 has no cysteine
            def ring_refused():
                row = "#peps .ent[data-i='2']"
                pg.select_option(row + " select.pr", "disulfide")
                pg.wait_for_function(
                    "() => /two cysteines/.test(document.querySelector("
                    "\"#peps .ent[data-i='2'] .seqmsg\").textContent)",
                    timeout=30000)
                pg.wait_for_timeout(600)
                shot_clip([row], "peptide-ring-refused.png", pad=8, least=760)
                # Back to three peptides docked as the run, for the pictures
                # after this one.
                pg.select_option(row + " select.pr", "")
                pg.select_option("#peps .ent[data-i='0'] select.pr", "")
                pg.click("#peps button[data-rm='3']")
                pg.wait_for_timeout(800)
            step("ring refused", ring_refused)

            # how every peptide is built and started
            def peptide_options():
                # Centered, so the labels above the fields are on screen too.
                pg.evaluate("document.getElementById('cyclic')"
                            ".scrollIntoView({block: 'center'})")
                pg.wait_for_timeout(500)
                shot_clip(['label[for="cyclic"]', "#cyclic",
                           'label[for="helical"]', "#helical", "#pepoptsnote"],
                          "peptide-options.png", pad=8, least=760)
            step("peptide options", peptide_options)

            # the docking and rescoring settings, and one help box
            def settings():
                card("#poses").scroll_into_view_if_needed()
                pg.wait_for_timeout(500)
                shot(card("#poses"), "docking-settings.png")
            step("settings", settings)

            # the last step: the run's name, its folder, and the button
            def submit_card():
                card("#submit").scroll_into_view_if_needed()
                pg.wait_for_timeout(500)
                shot(card("#submit"), "submit.png")
            step("submit card", submit_card)

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

            # Add measures, from the Runs card, with every measure and every
            # run that lacks one chosen.
            def add_measures():
                pg.click("#measuresopen")
                pg.wait_for_function(
                    "() => !/Reading/.test(document.getElementById('mslead')"
                    ".textContent)", timeout=60000)
                for box in ("#msmall", "#msrall"):
                    if pg.locator(box).is_enabled() and not pg.locator(box).is_checked():
                        pg.locator(box).click()
                pg.wait_for_timeout(600)
                shot(pg.locator("#measuresdlg"), "add-measures.png")
                pg.keyboard.press("Escape")
            step("add measures", add_measures)

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

            # the runs folder walked: inside the example run, then one pose
            # drawn from its docking folder
            def folders():
                pg.click("#tab_folders")
                pg.wait_for_selector("#flist .ftable", timeout=60000)
                pg.locator("#flist tr:has(.tag) button.fgo").first.click()
                pg.wait_for_function(
                    "() => /This folder is the run/.test("
                    "document.getElementById('fhere').textContent)", timeout=60000)
                # The tabs, the trail and the run's folder, without the run
                # opened below it by the step before.
                pg.evaluate("document.querySelector('.tabs')"
                            ".scrollIntoView({block: 'start'})")
                pg.wait_for_timeout(800)
                shot_clip([".tabs", "#flist"], "runs-folders.png", pad=14,
                          least=760)
            step("folders", folders)

            def file_viewer():
                pg.locator("#flist button.fgo", has_text="docking").first.click()
                at = ("() => (document.querySelector('#fcrumbs b') || {})"
                      ".textContent")
                pg.wait_for_function(at + " === 'docking'", timeout=60000)
                pg.locator("#flist tbody button.fgo").nth(1).click()
                pg.wait_for_function(at + " !== 'docking'", timeout=60000)
                pg.wait_for_selector("#flist button[data-kind=structure]",
                                     timeout=60000)
                pg.locator("#flist button.ghost[data-kind=structure]").first.click()
                pg.wait_for_selector("#filedlg[open] #file3d:not(.hide)",
                                     timeout=60000)
                # drawn in software
                pg.wait_for_timeout(5000)
                shot(pg.locator("#filedlg"), "file-viewer.png")
                pg.click("#fileclose")
                pg.click("#tab_runs")
            step("file viewer", file_viewer)

            # how to cite
            def cite():
                pg.evaluate("window.scrollTo(0, 0)")
                pg.click("#citebtn")
                pg.wait_for_selector("dialog[open] #citetext", timeout=60000)
                pg.wait_for_timeout(1200)
                shot(pg.locator("dialog[open]").first, "cite.png")
                pg.keyboard.press("Escape")
            step("cite", cite)

            # ADapt, in the Receptor card, when it is installed here
            def adapt():
                if not pg.locator("#adaptsec").is_visible() \
                        or not pg.evaluate("ADAPT_OK"):
                    print("  ADapt is not offered here, so no ADapt pictures",
                          flush=True)
                    return
                pg.evaluate("adaptFold(false)")
                # Opened on the receptor when it was loaded: its frame shows
                # the session once ADapt has started and read the file.
                pg.wait_for_function(
                    "() => { try { const d = document.getElementById('adaptframe')"
                    ".contentDocument; return !!d && d.querySelectorAll('#tree li').length > 0"
                    " && !!d.querySelector('#viewer canvas'); }"
                    " catch(e){ return false; } }", timeout=330000)
                pg.wait_for_function(
                    "() => /this page.s preparation/.test("
                    "document.getElementById('adaptsays').textContent)", timeout=60000)
                pg.wait_for_timeout(5000)
                shot(pg.locator("#adaptsec"), "adapt-inline.png")
                # The same session full screen, in a tab of its own.
                with ctx.expect_page(timeout=60000) as popup:
                    pg.click("#adaptfull")
                tab = popup.value
                tab.wait_for_load_state("load")
                tab.wait_for_function(
                    "() => document.querySelectorAll('#tree li').length > 0",
                    timeout=240000)
                tab.wait_for_timeout(5000)
                shot(tab, "adapt.png")
                tab.close()
                pg.bring_to_front()
                # The bound p53 peptide, chain B, taken out of the receptor in
                # ADapt's tree: the card switches to ADapt's receptor.
                frame = pg.frame_locator("#adaptframe")
                # A chain's own row has no box; the protein under it has.
                frame.locator("#tree").evaluate(
                    "t => { const li = [...t.querySelectorAll('li.kind-chain')].find("
                    "x => x.querySelector('.lab').textContent.trim() === 'B');"
                    " const cat = [...li.querySelectorAll('li.kind-category')].find("
                    "x => /protein/.test(x.querySelector('.lab').textContent));"
                    " (cat || li).querySelector('input.cb[data-col=rec]').dataset.shot = 'B'; }")
                frame.locator('input.cb[data-shot="B"]').click()
                pg.wait_for_function(
                    "() => /as set in ADapt/.test("
                    "document.getElementById('adaptsays').textContent)", timeout=60000)
                # Written now and read back, as a submission would.
                pg.click("#adaptnow")
                pg.wait_for_selector("#adaptuse:not(.hide)", timeout=300000)
                pg.wait_for_function(
                    "() => /atoms/.test(document.getElementById('recinfo').textContent)",
                    timeout=120000)
                pg.evaluate("adaptFold(true)")
                pg.wait_for_timeout(2000)
                pg.evaluate("document.getElementById('step-receptor')"
                            ".scrollIntoView({block: 'start'})")
                pg.wait_for_timeout(800)
                shot(card("#reccode"), "receptor-adapt.png")
            step("adapt", adapt)
            br.close()
    finally:
        try:
            os.killpg(srv.pid, signal.SIGTERM)
        except OSError:
            pass
        if temporary:
            shutil.rmtree(work, ignore_errors=True)
        else:
            print("the run and the receptor are in", work)
    print("\nmade:", " ".join(made))
    print("page errors:", errors[:5] or "none")
    print("steps that failed:", failed or "none")
    missing = [n for n in WANT if n not in made]
    print("missing:", missing or "none")
    return 1 if missing or errors or failed else 0


if __name__ == "__main__":
    sys.exit(main())
