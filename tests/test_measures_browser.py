"""Add measures in a real browser: a headless Chrome walks the dialog the way a
person does.

By Aadhil Haq

The interface is served from this process on a runs folder of made up runs,
with the scheduler stood in for, so Add queues nothing real. Chrome opens the
page and: opens Add measures from the Runs card, checks the focus is in it,
ticks Select all for the measures and for the runs, chooses a range with a
shift-click, lists Every run and back, finds a run by name, presses Add and
reads what was queued, then Done. It opens a run that lacks the strain, reads
the line over its table, opens the dialog from it on that run, and closes it
with Esc; and opens a run that has the strain, whose table carries the
Peptide strain column with the reason for an empty cell on hover.

It needs Chrome, chromedriver and selenium, named by ADCP_CHROME,
ADCP_CHROMEDRIVER and ADCP_SELENIUM (a folder to import selenium from, if it
is not installed); without them it is a skip. ADCP_CHROME_PROFILE is a folder
on the node's own disk for Chrome's profile. ADCP_SHOTS names a folder for
four screenshots: the dialog with selections, the line over a run's table, a
table with the strain filled in, and the dialog after Add. Run it on a compute
node: it starts a server and a browser.

    python3 tests/test_measures_browser.py
"""
import json
import os
import shutil
import sys
import tempfile
import threading
import time
from http.server import ThreadingHTTPServer

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
if os.environ.get("ADCP_SELENIUM"):
    sys.path.insert(0, os.environ["ADCP_SELENIUM"])

CHROME = os.environ.get("ADCP_CHROME", "")
DRIVER = os.environ.get("ADCP_CHROMEDRIVER", "")
if not (CHROME and os.path.isfile(CHROME) and DRIVER and os.path.isfile(DRIVER)):
    print("skip: Chrome or chromedriver is not named (ADCP_CHROME, ADCP_CHROMEDRIVER)")
    sys.exit(77)
try:
    from selenium import webdriver
    from selenium.webdriver.chrome.options import Options
    from selenium.webdriver.chrome.service import Service
    from selenium.webdriver.common.by import By
    from selenium.webdriver.common.action_chains import ActionChains
    from selenium.webdriver.common.keys import Keys
    from adcp_dock import job as jobmod, report, server, slurm
    from adcp_dock._vendor import gbengines, strain as strainmod
except ImportError as exc:
    print("skip: %s" % exc)
    sys.exit(77)

fails, ran = [], 0


def check(cond, why):
    global ran
    ran += 1
    print("   %s %s" % ("ok  " if cond else "FAIL", why), flush=True)
    if not cond:
        fails.append(why)


WORK = tempfile.mkdtemp(prefix="adcp_msbrowser_")
RUNS = os.path.join(WORK, "runs")
os.makedirs(RUNS)
SHOTS = os.environ.get("ADCP_SHOTS", "")
FILES = ("_com.parm7", "_rec.parm7", "_lig.parm7", "_min.rst7")


def make(name, peptides, created="2026-09-22 13:02:24", ran=None, report_rows=True,
         strain_vals=None):
    """A finished run of [(sequence, scored poses, has strain, has panel)]."""
    rd = os.path.join(RUNS, name)
    j = jobmod.default(mmgbsa=True, label=name, cpus=24,
                       peptides=[{"name": "pep%d" % i, "sequence": s}
                                 for i, (s, *_r) in enumerate(peptides, 1)])
    j.update(created=created, state="done")
    if ran:
        j["mmgbsa_engines_run"] = ran
    jobmod.save(rd, j)
    os.makedirs(os.path.join(rd, "results"))
    rows = []
    for i, (seq, n, has_strain, has_panel) in enumerate(peptides, 1):
        mm = [round(-30.0 - 3.1 * k - i, 4) for k in range(n)]
        r = {"name": "pep%d" % i, "peptide": seq, "engine": "ADCP 1.1.21",
             "modes": [{"rank": k, "score": -15.0 - k} for k in range(1, n + 1)],
             "mmgbsa": mm, "mmgbsa_ranks": list(range(1, n + 1)), "mmgbsa_engine": "amber"}
        if has_strain:
            vals = (strain_vals or {}).get(seq) or [20.0 + k for k in range(n)]
            r.update(strain=[{"pose": k + 1, "strain": v,
                              "why": None if v is not None else
                              "sander did not finish the relaxation"}
                             for k, v in enumerate(vals)],
                     strain_engine="amber", strain_version=strainmod.VERSION)
        if has_panel:
            r["panel"] = {"pose": 1, "bsa": 800.0 + 10 * i, "n_hbond": 4, "n_ionic": 1,
                          "n_clash": 0, "bsa_per_res": 120.0}
        json.dump(r, open(os.path.join(rd, "results", "%s.json" % seq), "w"), indent=1)
        for k in range(1, n + 1):
            pd = os.path.join(rd, "mmgbsa", seq, "pose%d" % k)
            os.makedirs(pd)
            for t in FILES:
                open(os.path.join(pd, seq + t), "w").write("x\n")
        rows.append(report.summarise(seq, r["modes"], mm, name=r["name"],
                                     engine=r["engine"], mmgbsa_engine="amber",
                                     panel=r.get("panel"), mm_ranks=r["mmgbsa_ranks"],
                                     strain=r.get("strain"), strain_engine="amber"))
    if report_rows and rows:
        report.write(rows, rd)
    return rd


make("mdm2_round2_20260922_181436",
     [("RFMDYWEGL", 5, False, True), ("TSFAEYWNLLS", 5, False, True)])
lacking = make("p53_variants_20260922_154712",
               [("ETFSDLWKLLPEN", 5, False, False), ("LTFEHYWAQLTS", 5, False, False)])
make("short_check_20261018_202009", [("CRWYPC", 3, False, True)], created="2026-10-18 20:20:09")
make("openmm_check_20261009_091813", [("SQETFSDLWKLL", 5, False, True)],
     created="2026-10-09 09:18:13", ran=["amber", "openmm"])
strained = make("test3_20260922_130224",
                [("LLKWHY", 5, True, True), ("WHWHLLPF", 2, True, True)],
                strain_vals={"LLKWHY": [41.2, 30.5, 22.07, 24.8, 19.31],
                             "WHWHLLPF": [27.4, None]})
make("test6_20260922_131946", [("GSWHLLP", 2, True, True)])
blocked = os.path.join(RUNS, "maps_only_20260922_175057")
jobmod.save(blocked, dict(jobmod.default(label="maps only"), created="2026-09-22 17:50:57"))
os.makedirs(os.path.join(blocked, "results"))
json.dump({"peptide": "KKKK", "mmgbsa": []}, open(os.path.join(blocked, "results", "KKKK.json"), "w"))

# The scheduler, stood in for: every job is given an id and nothing runs.
asked = []
ids = iter(range(20370101, 20379999))
slurm.available = lambda: True
slurm.local_mode = lambda: False
slurm.state = lambda jid, rundir=None: {}
slurm.within_limit = lambda t, partition="", extra="": t
slurm.cancel = lambda jid: ""
slurm.why_pending = lambda jid, rundir=None: ""


def fake_submit(script, depends_on=None, kind="afterok"):
    jid = str(next(ids))
    asked.append((os.path.basename(os.path.dirname(script)), os.path.basename(script), jid))
    return jid


slurm.submit = fake_submit
gbengines.OPENMM.available = lambda: True
gbengines.SCHRODINGER.available = lambda: True
gbengines.SCHRODINGER.plan = lambda timeout=30, max_age=None: {"readable": False}

app = server.App(RUNS)
httpd = ThreadingHTTPServer(("127.0.0.1", 0), server.Handler)
httpd.app = app
httpd.token = app.token
threading.Thread(target=httpd.serve_forever, daemon=True).start()
URL = "http://127.0.0.1:%d/?t=%s" % (httpd.server_address[1], app.token)

def chrome(attempt):
    """A headless Chrome, on a fresh profile."""
    opts = Options()
    opts.binary_location = CHROME
    # A profile on a node's own disk where one is named: made on a shared file
    # system it takes Chrome about a minute, and chromedriver gives up first.
    prof = os.path.join(os.environ.get("ADCP_CHROME_PROFILE") or WORK, "profile%d" % attempt)
    for a in ("--headless=new", "--no-sandbox", "--disable-dev-shm-usage",
              "--window-size=1500,1100", "--no-first-run", "--hide-scrollbars",
              "--user-data-dir=" + prof):
        opts.add_argument(a)
    opts.set_capability("goog:loggingPrefs", {"browser": "ALL"})
    return webdriver.Chrome(service=Service(DRIVER), options=opts)


# Chrome now and then fails to start on a busy node (its DevToolsActivePort is
# not written in time); that is the browser, not the page, so it is tried again.
d = None
for attempt in range(3):
    try:
        d = chrome(attempt)
        break
    except Exception as exc:                                    # noqa: BLE001
        print("note  Chrome did not start (%s), trying again" % str(exc).splitlines()[0],
              flush=True)
        time.sleep(3)
if d is None:
    httpd.shutdown()
    shutil.rmtree(WORK, ignore_errors=True)
    print("test_measures_browser FAILED: Chrome would not start three times")
    sys.exit(1)
js = d.execute_script
print("browser:", d.capabilities.get("browserName"), d.capabilities.get("browserVersion"),
      flush=True)


def wait(cond, most=30):
    t0 = time.time()
    while time.time() - t0 < most:
        try:
            if cond():
                return True
        except Exception:                                       # noqa: BLE001
            pass
        time.sleep(0.25)
    return False


def lead():
    return d.find_element(By.ID, "mslead").text


def shot(name, el=None):
    if not SHOTS:
        return
    os.makedirs(SHOTS, exist_ok=True)
    path = os.path.join(SHOTS, name)
    js("document.querySelectorAll('.creditbar').forEach(b => b.style.visibility = 'hidden')")
    if el is not None:
        el.screenshot(path)
    else:
        d.save_screenshot(path)
    js("document.querySelectorAll('.creditbar').forEach(b => b.style.visibility = '')")
    print("   shot", path, flush=True)


def rows():
    return d.find_elements(By.CSS_SELECTOR, "#msruns tbody tr[data-i]")


try:
    d.get(URL)
    check(wait(lambda: d.find_elements(By.CSS_SELECTOR, "#runs table")),
          "the page lists the runs")
    # --- inside one run --------------------------------------------------------------
    js("openRun(arguments[0])", lacking)
    check(wait(lambda: os.path.basename(lacking) in d.find_element(By.ID, "rundetail").text
               and d.find_elements(By.ID, "msline")), "a run that lacks the strain has the line")
    line = d.find_element(By.ID, "msline")
    check(line.text.startswith("This run has no peptide strain or interface measures. "
                               "One short job on a CPU node adds both.")
          and "Add them" in line.text, "the line says what the run lacks: %r" % line.text)
    js("arguments[0].scrollIntoView({block: 'center'})", line)
    time.sleep(0.4)
    shot("adcp_measures_2_line.png", d.find_element(By.ID, "rundetail"))
    d.find_element(By.ID, "msadd").click()
    check(wait(lambda: js("return document.getElementById('measuresdlg').open")
               and lead().startswith("Peptide strain and Interface measures for 1 run")),
          "Add them opens the dialog on this run with those measures ticked: %r" % lead())
    sel = [e.text for e in d.find_elements(By.CSS_SELECTOR, "#msruns tbody tr.sel td.msrun")]
    check(sel == [os.path.basename(lacking)], "and that run chosen: %s" % sel)
    d.switch_to.active_element.send_keys(Keys.ESCAPE)
    check(wait(lambda: not js("return document.getElementById('measuresdlg').open")),
          "Esc closes the dialog")

    # --- a run with the strain -------------------------------------------------------
    js("openRun(arguments[0])", strained)
    check(wait(lambda: os.path.basename(strained) in d.find_element(By.ID, "rundetail").text
               and d.find_elements(By.CSS_SELECTOR, ".rtable table.rt")),
          "the run with the strain shows its table")
    heads = [e.text for e in d.find_elements(By.CSS_SELECTOR, ".rtable tr.heads th")]
    check("PEPTIDE STRAIN" in [h.upper() for h in heads],
          "the table has a Peptide strain column: %s" % heads)
    th = [e for e in d.find_elements(By.CSS_SELECTOR, ".rtable tr.heads th")
          if e.text.upper() == "PEPTIDE STRAIN"]
    check(th and th[0].get_attribute("title").startswith("What the peptide pays to hold"),
          "its heading carries the help")
    txt = d.find_element(By.CSS_SELECTOR, ".rtable").text
    empty = d.find_elements(By.CSS_SELECTOR, ".rtable span[title='sander did not finish the relaxation']")
    check("19.31" in txt and empty,
          "the lowest energy pose's strain is shown, and an empty cell gives its reason on hover")
    check(not d.find_elements(By.ID, "msline"), "a run that has it all has no line over its table")
    js("arguments[0].scrollIntoView({block: 'center'})", d.find_element(By.CSS_SELECTOR, ".rtable"))
    time.sleep(0.4)
    shot("adcp_measures_3_table.png", d.find_element(By.ID, "rundetail"))
    # --- opened from the Runs card ------------------------------------------------
    btn = d.find_element(By.ID, "measuresopen")
    js("arguments[0].scrollIntoView({block: 'center'})", btn)
    btn.click()
    check(wait(lambda: "Reading" not in lead() and rows()), "the dialog opens and reads the runs")
    check(js("return document.getElementById('measuresdlg').open"), "it is open, as a dialog")
    check(lead().startswith("Now choose the runs."),
          "with the light measures ticked, the line says to choose the runs: %r" % lead())
    act = js("const a = document.activeElement; return a ? a.dataset.key || a.id : ''")
    check(act == "m:strain", "the focus is on the first measure: %r" % act)
    counts = [e.text for e in d.find_elements(By.CSS_SELECTOR, "#msmeasures .mscount")]
    check(counts == ["missing in 4", "missing in 1", "missing in 5", "missing in 6"],
          "each measure says how many runs lack it: %s" % counts)
    # The keyboard: Space unticks the focused measure, and the focus stays on it.
    d.switch_to.active_element.send_keys(Keys.SPACE)
    check(wait(lambda: not js("return document.querySelector('[data-key=\"m:strain\"]').checked"))
          and js("return document.activeElement.dataset.key") == "m:strain",
          "Space toggles a measure and keeps the focus on it")
    mall = d.find_element(By.ID, "msmall")
    check(js("return arguments[0].indeterminate", mall), "Select all shows some ticked")
    d.find_element(By.CSS_SELECTOR, "label.msall").click()
    ticked = js("return [...document.querySelectorAll('#msmeasures input')]"
                ".filter(b => b.checked).length")
    check(ticked == 4 and mall.is_selected() and not js("return arguments[0].indeterminate", mall),
          "Select all ticks all four measures: %s" % ticked)
    n_listed = len(rows())
    rall = d.find_element(By.ID, "msrall")
    check(d.find_element(By.ID, "msralllab").text == "Select all %d listed" % n_listed,
          "the runs' Select all says how many it takes: %r"
          % d.find_element(By.ID, "msralllab").text)
    rall.click()
    sel = len(d.find_elements(By.CSS_SELECTOR, "#msruns tbody tr.sel"))
    check(sel == n_listed and lead().startswith("Peptide strain, Interface measures,"),
          "Select all chooses every run listed and the line says what Add will queue: "
          "%d of %d, %r" % (sel, n_listed, lead()))
    rall.click()
    check(not d.find_elements(By.CSS_SELECTOR, "#msruns tbody tr.sel"),
          "Select all again chooses none")
    r = rows()
    r[1].find_element(By.CSS_SELECTOR, "td.msrun").click()
    ActionChains(d).key_down(Keys.SHIFT).click(rows()[3].find_element(
        By.CSS_SELECTOR, "td.msrun")).key_up(Keys.SHIFT).perform()
    sel = [e.get_attribute("data-i") for e in d.find_elements(By.CSS_SELECTOR, "#msruns tbody tr.sel")]
    check(sel == ["1", "2", "3"], "a shift-click chooses the range between two runs: %s" % sel)
    check(js("return arguments[0].indeterminate", rall), "the runs' Select all shows some chosen")
    tally = d.find_element(By.ID, "mstally").text
    check(tally.startswith("3 of %d runs chosen" % n_listed), "the tally counts them: %r" % tally)
    d.find_element(By.ID, "msfall").click()
    every = rows()
    blocked_rows = d.find_elements(By.CSS_SELECTOR, "#msruns tbody tr.blocked")
    check(len(every) == 7 and any("no poses were scored" in b.text for b in blocked_rows)
          and d.find_element(By.ID, "msfall").get_attribute("aria-pressed") == "true",
          "Every run lists the runs with nothing to add too, the blocked one with its "
          "reason: %d rows" % len(every))
    chk = blocked_rows[0].find_element(By.CSS_SELECTOR, "input") if blocked_rows else None
    check(chk is not None and not chk.is_enabled(), "a blocked run cannot be chosen")
    d.find_element(By.ID, "msfmissing").click()
    find = d.find_element(By.ID, "msfind")
    find.send_keys("p53")
    check(wait(lambda: len(rows()) == 1) and "p53" in rows()[0].text
          and "not listed" in d.find_element(By.ID, "mstally").text,
          "Find a run narrows the list and the tally says how many are not listed")
    find.send_keys(Keys.BACKSPACE * 3)
    check(wait(lambda: len(rows()) == n_listed), "clearing the search lists them again")
    shot("adcp_measures_1_dialog.png", d.find_element(By.ID, "measuresdlg"))
    go = d.find_element(By.ID, "msgo")
    check(go.text == "Add to 3 runs" and go.is_enabled(), "Add says how many runs: %r" % go.text)
    go.click()
    check(wait(lambda: lead().startswith("Queued")), "Add queues the jobs: %r" % lead())
    queued = d.find_element(By.ID, "msqueued").text
    check("job 2037" in queued and "then its report" in queued and go.text == "Done"
          and not d.find_element(By.ID, "mscancel").is_displayed(),
          "after Add the dialog lists what was queued for each run, and Done: %r" % queued[:300])
    scripts = sorted({s for _r, s, _j in asked})
    check("04_measures.sbatch" in scripts and "03_report.sbatch" in scripts
          and "06_tidy.sbatch" in scripts,
          "the stood in scheduler was given the short job, the report and the tidy: %s" % scripts)
    shot("adcp_measures_4_queued.png", d.find_element(By.ID, "measuresdlg"))
    go.click()
    check(wait(lambda: not js("return document.getElementById('measuresdlg').open")),
          "Done closes the dialog")

    errs = [e for e in (d.get_log("browser") if hasattr(d, "get_log") else [])
            if e.get("level") == "SEVERE" and "favicon" not in e.get("message", "")]
    check(not errs, "the page logged no errors: %s" % errs[:3])
finally:
    d.quit()
    httpd.shutdown()
    shutil.rmtree(WORK, ignore_errors=True)

if fails:
    print("test_measures_browser FAILED (%d of %d)" % (len(fails), ran))
    for f in fails:
        print("  FAIL " + f)
    sys.exit(1)
print("ok  Add measures works in Chrome as a person uses it, %d checks" % ran)
