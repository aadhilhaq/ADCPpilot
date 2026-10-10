"""The web page: one file, no build step, no network access at run time.

By Aadhil Haq

The peptide list is laid out as a row per sequence, with the entry number,
the type, and a sequence box that numbers its own residues.

Served by the same process that performs the work, so the page cannot be a
version ahead of the server it is talking to.
"""
import os

from . import docskit

PAGE = r"""<!doctype html><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>ADCP peptide docking</title>
<!-- The same mark the launcher carries. It is inlined as a data URI, so it
     needs no route, no second file and no cache. -->
<link rel="icon" type="image/svg+xml" href="__FAVICON__">
<!-- Safari takes SVG favicons only from version 16, and falls back to the
     last icon it understands rather than the first. -->
<link rel="icon" type="image/png" sizes="128x128" href="__FAVICON_PNG__">
<link rel="apple-touch-icon" href="__FAVICON_PNG__">
<script>/* ?embed=1: shown inside a notebook frame */
if (new URLSearchParams(location.search).get("embed"))
  document.documentElement.classList.add("embed");
/* The light or dark choice made with the switch in the header, restored
   before the page is drawn. Without one, the system setting decides. */
try { const t = localStorage.getItem("adcp.theme");
      if (t === "light" || t === "dark") document.documentElement.dataset.theme = t;
} catch (e) {}</script>
<style>
:root{
  color-scheme:light;
  --bg:#f7f8fb; --card:#fff; --panel:#eceef6; --line:#dfe3ee; --fg:#1b2333;
  --dim:#66718a; --accent:#1a56db; --accent2:#e8efff;
  --hydro:#b4632a; --polar:#1f3a6e; --bad:#c0392b; --warn:#c2571a;
  --good:#1a7f5a; --soft:#f2f4fa;
  /* the heading band across the top of a card, and the ink on it */
  --ribbon:#eef1f9; --ribbonfg:#2b3550;
  /* A pill's own border, which has to read against --card in both themes,
     and the ink that goes on a filled one. */
  --pill-line:#c3cbdd; --on-accent:#ffffff;
}
/* Dark: the system's setting, unless the switch in the header chose light;
   and always when the switch chose dark. */
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){__DARK_VARS__}}
:root[data-theme=dark]{__DARK_VARS__}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:14px/1.5 -apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
header{position:relative;background:var(--card);border-bottom:1px solid var(--line);
  padding:18px 24px 14px;text-align:center}
/* Docs in the corner, where a documentation link is looked for, and out of
   the centred block so it keeps its place whatever the title block holds. */
.docscorner{position:absolute;top:12px;right:16px}
.brandblock{max-width:760px;margin:0 auto}
.brandblock h1{margin:0 0 6px}
.brandblock .sub{display:block;line-height:1.55}
.brandblock a{color:var(--accent);text-decoration:none}
.brandblock a:hover{text-decoration:underline}
.headbar{margin-top:10px;display:flex;gap:12px;align-items:center;
  justify-content:center;flex-wrap:wrap}
h1{font-size:20px;margin:0;font-weight:650;letter-spacing:-.01em}
.sub{color:var(--dim);font-size:12.5px}
main{max-width:1180px;margin:0 auto;padding:20px 24px 80px}
/* On a wide window, and in a frame that fills one, the form has room
   for more than a 1180px column. The structure viewer and the results
   table are the two things that use it. */
@media (min-width:1500px){main{max-width:1440px}}
@media (min-width:1900px){main{max-width:1720px}}
/* Embedded in a notebook: the notebook shows the banner and the credit line
   directly above the frame, so the header keeps only its buttons, and the
   form uses the whole width of the frame. */
.embed header{padding:8px 18px;text-align:left}
.embed .brandblock{max-width:none}
.embed .brandblock h1,.embed .brandblock>.sub{display:none}
/* Except what a paper must cite: the notebook's banner names ADCP but not
   the papers, and this is the one line a methods section depends on. */
.embed .brandblock>.sub.mustcite{display:block;text-align:left;margin-bottom:6px}
.mustcite .ref{white-space:nowrap}
.embed .headbar{margin-top:0;justify-content:flex-end}
.embed main{max-width:none;padding:14px 18px 60px}
/* Centred by default. Sequences, the residue strip, logs and the contact
   list keep their left edge, because reading those depends on where the
   line starts. */
.card{background:var(--card);border:1px solid var(--line);border-radius:14px;
  padding:18px 20px;margin-bottom:18px;text-align:center}
.card .row{justify-content:center}
.card label{text-align:left}
.seq,.seq2,pre,#contactbox,#browser,#chainbox,.ent{text-align:left}
.seq2 .ln,#contactbox{justify-content:flex-start}
/* ---- a card's heading is part of the card ----------------------------
   A ribbon across the card's head: a band that fills the card from border to
   border, carries the card's top radius and has a 2px accent edge.

   The band is full-bleed, so its negative margins must equal .card's padding
   and its radius must equal .card's radius minus its border:
       .card padding 18px 20px          margin:-18px -20px
       .card radius 14px, border 1px    border-radius:13px 13px 0 0
   Change one and change the other.

   The accent edge is an inset shadow, so it follows the rounded corners and
   adds nothing beneath the card that could stack with the box-shadow of
   #helpbox and #issuebox.

   .card>h2 rather than a bare h2, so the band applies only to a card's own
   heading. The help and issue panels head themselves with an h3 and are
   unaffected. */
.card>h2{
  margin:-18px -20px 14px;
  padding:12px 20px 11px;
  border-radius:13px 13px 0 0;
  background:var(--ribbon);
  border-bottom:1px solid var(--line);
  box-shadow:inset 0 2px 0 var(--accent);
  font-size:12.5px;font-weight:700;text-transform:uppercase;
  letter-spacing:.1em;color:var(--ribbonfg);text-align:center}
/* The line under a heading is a subtitle rather than body text. */
/* A subtitle, not a paragraph: one sentence under a heading, kept on one line
   at desktop width and wrapped on a phone. */
.hint{color:var(--dim);font-size:12.5px;line-height:1.55;margin:0 auto 14px}

/* The picker is a dialog rather than a list inside the page. An inline
   listing pushes everything below it down as it opens, so the control that
   opened it moves from under the pointer, and a large directory pushes the
   rest of the form off the screen. */
dialog{border:1px solid var(--line);border-radius:12px;padding:0;
  width:min(680px,92vw);background:var(--card);color:var(--fg)}
dialog::backdrop{background:rgba(0,0,0,.45)}
.dlg-h{padding:12px 16px;border-bottom:1px solid var(--line);
  display:flex;gap:10px;align-items:baseline;flex-wrap:wrap}
.dlg-h b{font-size:14px}
.crumbs{font-size:12px;color:var(--dim);word-break:break-all;flex:1 1 100%}
.crumbs button{background:0;border:0;color:var(--accent);cursor:pointer;
  padding:0 2px;font:inherit}
/* Every scrolling list has a maximum height, for the same reason. */
.listing{max-height:56vh;overflow:auto;padding:6px 0}
.listing button{display:flex;width:100%;gap:10px;background:0;border:0;
  color:inherit;font:inherit;text-align:left;padding:7px 16px;cursor:pointer;
  align-items:center;border-radius:0}
.listing button:hover{background:var(--soft)}
.listing .sz{margin-left:auto;color:var(--dim);font-size:12px}
.listing .place{font-weight:600}
.dlg-f{padding:10px 16px;border-top:1px solid var(--line);display:flex;
  gap:10px;align-items:center;justify-content:space-between;flex-wrap:wrap}
.dlg-f>span{flex-wrap:wrap}

/* A job's state at a glance: a running job's dot moves, because a static one
   looks finished. */
.dot{width:9px;height:9px;border-radius:50%;flex:0 0 auto;
  background:var(--dim);display:inline-block}
.dot.running{background:var(--good);animation:pulse 1.1s ease-in-out infinite}
.dot.queued{background:#c8922f}
.dot.done{background:var(--good)}
.dot.failed{background:var(--warn)}
.dot.gone{background:var(--dim)}
@keyframes pulse{0%,100%{opacity:1}50%{opacity:.25}}
.jobrow{display:flex;align-items:center;gap:9px;padding:5px 2px;
  font-size:12px;border-bottom:1px solid var(--line)}
.jobrow:last-child{border-bottom:0}
.jobrow .jname{flex:0 0 96px;color:var(--dim)}
.jobrow .jid{font-variant-numeric:tabular-nums;flex:0 0 110px}
.jobrow .jst{flex:1 1 auto}
/* A required box left empty is marked where it is, because a message at the
   bottom of the form is often out of sight of the box. */
/* The theme button, matching the launcher's: an icon that says which way
   the switch goes, and a one word label. */
#themebtn{display:inline-flex;align-items:center;gap:6px}
#themebtn svg{width:13px;height:13px;fill:currentColor;flex:0 0 auto}
#themebtn .moon{display:none}
:root[data-theme=dark] #themebtn .moon{display:inline}
:root[data-theme=dark] #themebtn .sun{display:none}
@media (prefers-color-scheme:dark){
  :root:not([data-theme=light]) #themebtn .moon{display:inline}
  :root:not([data-theme=light]) #themebtn .sun{display:none}}
.pickfilter{display:flex;align-items:center;gap:9px;padding:8px 14px;
  border-bottom:1px solid var(--line)}
.pickfilter input{flex:1 1 auto;min-width:0}
.missing,.missing:focus{border-color:var(--warn) !important;
  outline:2px solid var(--warn)}
/* A container is not an input. #peps has no border for a border colour to
   change and no height when it is empty, so an outline traced nothing and
   the page said something was missing while showing nothing in red. */
.missingbox{border:1px dashed var(--warn);border-radius:6px;
  min-height:34px;padding:6px 9px;outline:none !important}
/* The note must not move anything. In flow it would make its row taller and
   shift every control beside it, so it is taken out of flow and placed under
   its box. */
.fieldwrap{position:relative}
.missingnote{position:absolute;left:0;top:100%;margin-top:3px;
  color:var(--warn);font-size:11.5px;line-height:1.3;pointer-events:none}
/* and the row leaves room for one, always, so nothing moves when it
   appears or goes away */
.row.hasnotes{margin-bottom:26px}
/* The label keeps its "?" beside the text rather than letting it wrap onto a
   line of its own, with a real gap rather than a margin on the marker. */
label{display:block;font-size:11.5px;color:var(--dim);margin-bottom:4px}
/* Flex, not inline-flex. .card sets text-align:center and a label inherits it
   through .row>div. A block label still starts at the left edge of its
   control, but an inline-level box is centred by text-align, so inline-flex
   would move every label with a "?" into the middle of its column.
   Block-level flex keeps the label full width with its content at the start,
   and still stops the marker wrapping onto its own line.
       label                    (0,0,1)  display:block
       label:has(> button.q)    (0,1,2)  wins: :has() takes the specificity
                                         of its argument, "> button.q" (0,1,1)
   The checkboxes are the exception: their labels wrap the input and sit
   inline, several to a row, so they keep inline-flex and the 7px gap .checks
   sets. (0,1,0)+(0,1,0)+(0,1,1) = (0,3,1) beats (0,1,2). */
label:has(> button.q){display:flex;align-items:center;gap:.3em;
  justify-content:flex-start}
.checks .chk:has(> button.q){display:inline-flex;gap:7px}
/* :has() is recent; without it the marker still sits inline beside the text.
   Nothing depends on the flex. */
input[type=text],input[type=number],select,textarea{font:inherit;
  padding:8px 11px;border:1px solid var(--line);border-radius:9px;
  background:var(--card);color:var(--fg);width:100%}
input:focus,select:focus,textarea:focus{outline:2px solid var(--accent2);
  border-color:var(--accent)}
textarea{font-family:ui-monospace,Menlo,monospace;font-size:12.5px;min-height:76px}
.row{display:flex;gap:12px;flex-wrap:wrap;align-items:flex-end}
/* A row of checkboxes is not a row of fields: they are laid out as one group,
   without spacer labels to push them level with inputs. */
.checks{display:flex;flex-wrap:wrap;gap:10px 22px;align-items:center;
  margin-top:14px}
.checks .chk{display:inline-flex;align-items:center;gap:7px;margin:0;
  font-size:12.5px;color:var(--fg);cursor:pointer;line-height:1.2}
.checks .chk input{width:auto;margin:0;flex:0 0 auto}
/* The receptor's preparation, as one group under its chains. A panel rather
   than a heading of its own: it belongs to the Receptor card, and a second
   ribbon inside a card would read as a second card. */
.prepgroup{margin-top:12px;padding:10px 14px 12px;background:var(--soft);
  border:1px solid var(--line);border-radius:10px;text-align:left}
.grouphead{font-size:10.5px;text-transform:uppercase;letter-spacing:.06em;
  color:var(--dim);font-weight:650}
.grouphead span{text-transform:none;letter-spacing:0;font-weight:400;
  margin-left:6px}
.prepgroup .checks{margin-top:8px}
/* The pH is a number inside a step, so it is sized to its four characters
   rather than to a form column. */
.prepsteps input[type=number]{width:74px;padding:4px 8px}
/* A choice inside a step is as wide as its words, not a form column, and
   may wrap on a phone rather than push the panel wider than the screen. */
.prepsteps select{width:auto;max-width:100%;padding:4px 8px}
.prepsteps .chk{flex-wrap:wrap;max-width:100%}
.prepgroup .why:empty{display:none}
#prepnote,#mapsnote{margin-top:8px}
/* The box switch sits on the viewer it draws into. The site row reports the
   box as a centre and edge lengths, which cannot be judged against a groove
   by eye, so the box is drawn in the picker over the structure, and the tick
   that hides it goes on the viewer's caption line, beside what it controls.
   Not .checks: that is a block of settings with its own 14px top margin and
   22px column gap, and this is one tick at the end of a caption. */
.pickhead{display:flex;align-items:baseline;justify-content:space-between;
  gap:10px 16px;flex-wrap:wrap;margin-bottom:5px}
.pickhead .why{flex:1 1 210px;margin:0;text-align:left}
/* Both selectors, because label:has(> button.q) is (0,1,2) and would
   otherwise win the display property and take the tick full width the
   moment this field is given a "?". .pickhead label:has(...) is (0,2,2). */
.pickhead label,.pickhead label:has(> button.q){
  display:inline-flex;align-items:center;gap:7px;margin:0;flex:0 0 auto;
  font-size:12px;color:var(--fg);cursor:pointer;line-height:1.2;
  white-space:nowrap}
.pickhead label input{width:auto;margin:0;flex:0 0 auto}
/* The two ticks that change what the view draws, side by side, wrapping
   under each other on a phone rather than past the edge. */
.pickticks{display:inline-flex;gap:6px 16px;flex-wrap:wrap;flex:0 1 auto}
#surfnote{margin-top:4px;font-size:11.5px;text-align:left}
#surfnote:empty{display:none}
/* What stays a choice of the structure file and not of ADapt's receptor is
   greyed out while ADapt's is used, as is the pH while reduce is chosen. */
.prepsteps .chk.off{opacity:.45}
#adaptuse .bar{margin-top:8px}
/* ADapt in the Receptor card: the line saying which receptor is docked, its
   buttons, and the frame, as tall as ADapt's tree and view need. */
#adaptsays b{font-weight:600}
#adaptacts{justify-content:flex-start;flex-wrap:wrap;gap:8px}
#adaptacts:not(:has(button:not(.hide))){display:none}
.adaptbar{display:flex;align-items:center;gap:6px 12px;flex-wrap:wrap;
  margin-top:10px}
.adaptlab{margin:0;font-weight:600;font-size:12.5px;color:var(--fg)}
.adaptbar .why{flex:1 1 200px;margin:0}
/* Always holding its place, so the buttons stay at the right end. */
.adaptbar #adaptwait{display:block}
#adaptfull.hide{display:none}
#adaptframewrap{margin-top:8px;border:1px solid var(--line);border-radius:10px;
  overflow:hidden;background:var(--card)}
#adaptframe{display:block;width:100%;height:min(78vh,760px);min-height:520px;
  border:0;background:var(--card);transition:opacity .5s ease}
#adaptframewrap{position:relative}
#adaptframewrap.ld-under #adaptframe{opacity:0}
/* Folded, the stage is folded with the frame and the line above says it. */
#adaptpanel:not(.isfolded) #adaptwait .ld-inline{display:none}
/* Folded, the frame keeps its size out of sight, so ADapt's 3D view is not
   drawn into a box of no size and need not be drawn again. */
#adaptframewrap.folded{height:0;border:0;margin:0;visibility:hidden}
#adaptsec.adapton{border-color:var(--accent)}
/* A button in a row is laid out by .row button, which outranks .hide. */
.row button.hide{display:none}
#adaptwhat b{font-weight:600}
#adaptwhat code{word-break:break-all}
.row>div{flex:1 1 150px;min-width:0}
/* Every control in a row has the same height. Otherwise a .ghost button (5px
   padding, 12.5px text) is shorter than the input beside it, and bottom
   alignment leaves it visibly low. A shared height lines up a row of mixed
   controls at both edges. */
.row input[type=text],.row input[type=number],.row select,
.row button,.row .fixedval{height:38px;box-sizing:border-box}
.row button{display:inline-flex;align-items:center;justify-content:center}
.row .stepper button{height:38px}
/* ---- the Docking and rescoring card: four parts on a grid -------------
   Search, Engine, Rescoring and Resources, each under a small heading in the
   style of the Receptor card's groups, with a hairline between parts rather
   than a panel round each, which inside a card would read as a second card.

   A part's row is a grid: three equal columns for the search and the
   engine, four for the rescoring and the resources. Of the four, the middle
   two hold the menus and are a little wider than the outer two, which hold
   numbers and the button; at 1000 px equal columns cut off "AmberTools, the
   reference". Every field starts on its column's edge, and the row is
   bottom aligned, so a label that wraps rises and its control stays level
   with the others.

   A caption is a grid item of its own, in the column of the field it
   describes. Inside the field it would make that field taller than its
   neighbours and lift every other label to meet it. The row gap is the
   space between rows of fields once they wrap; a caption takes most of it
   back, so it sits close under its field.

   Narrower, four columns become two, and the engine's two fields share a
   line. On a phone every part is one column, the fields in the order they
   are written, so each caption comes straight after its own field. */
.dpart{text-align:left}
.dpart+.dpart{margin-top:16px;padding-top:14px;border-top:1px solid var(--line)}
.dpart>.grouphead{margin:0 0 9px}
.dpart>.row{display:grid;gap:12px;align-items:end}
.dpart.cols3>.row{grid-template-columns:repeat(3,minmax(0,1fr))}
.dpart.cols4>.row{grid-template-columns:minmax(0,.85fr) minmax(0,1.15fr)
  minmax(0,1.15fr) minmax(0,.85fr)}
.dpart .dcap{align-self:start;margin-top:-7px}
.dpart>.why{margin-top:7px}
.dpart>.why:empty{display:none}
.dpart>.checks{margin-top:10px}
.dpart>.grouphead+.checks{margin-top:0}
.dpart>.checks+.row{margin-top:10px}
#enginenote{grid-row:2;grid-column:1}
#mmcap{grid-row:2;grid-column:2/span 2}
#autofill{width:100%}
/* Rescoring off: its settings stay where they are, greyed and disabled, as
   the receptor's preparation greys a setting that does not apply. */
#mmopts.off{opacity:.45}
@media (max-width:940px){
  .dpart.cols4>.row,.dpart.dengine>.row{
    grid-template-columns:repeat(2,minmax(0,1fr))}
  #mmcap{grid-row:3;grid-column:1}
}
@media (max-width:600px){
  .dpart.cols3>.row,.dpart.cols4>.row,.dpart.dengine>.row{
    grid-template-columns:minmax(0,1fr)}
  #enginenote,#mmcap{grid-row:auto;grid-column:auto}
}
/* ---- buttons: one hierarchy -------------------------------------------
   Four kinds, so a person new to the page can tell at a glance what to press:
     primary    filled in the accent: the one action a card's step waits for
     secondary  outlined in the accent, 36px tall at 14px: an action someone
                needs, that is not the card's main one
     ghost      outlined in the ink: utilities, such as close, copy, refresh,
                and the extra actions on a row
     a link     only for what nobody needs in order to carry on
   None is pale: a needed action never sits in grey text on a faint line.
   Every one shows the keyboard's focus. */
button{font:inherit;padding:8px 16px;border-radius:999px;cursor:pointer;
  border:1px solid var(--pill-line);background:var(--card);color:var(--fg)}
button:hover{background:var(--soft)}
button:focus-visible,a.btnlink:focus-visible{outline:2px solid var(--accent);
  outline-offset:2px}
button.primary{background:var(--accent);border-color:var(--accent);
  color:var(--on-accent);font-weight:600;padding:10px 26px}
button.primary:hover{background:var(--accent);filter:brightness(1.08)}
button.secondary,a.btnlink.secondary{min-height:36px;font-size:14px;
  font-weight:600;padding:6px 16px;color:var(--accent);background:var(--card);
  border:1.5px solid var(--accent);line-height:1.3}
button.secondary:hover,a.btnlink.secondary:hover{background:var(--accent2);
  color:var(--accent)}
a.btnlink.secondary{display:inline-flex;align-items:center}
/* In a table a button keeps its words on one line. */
td button.secondary{white-space:nowrap}
button.ghost{font-size:13px;padding:5px 13px;color:var(--fg);min-height:32px}
button.ghost.danger{color:var(--bad);
  border-color:color-mix(in srgb,var(--bad) 55%,var(--pill-line))}
/* The peptide menu and the non-standard residue dialog, laid out as AlphaFold
   Server's: a menu under the row's dots, then a dialog in which a residue is
   clicked to add a row with a type chosen from that residue's options. */
.pepmenu{position:absolute;z-index:40;background:var(--card);border:1px solid var(--line);
  border-radius:10px;box-shadow:0 6px 24px rgba(0,0,0,.14);padding:6px 0;min-width:250px}
.pepmenu button{display:flex;gap:12px;align-items:center;width:100%;border:0;border-radius:0;
  background:none;padding:9px 16px;font-size:13.5px;color:var(--fg);text-align:left;cursor:pointer}
.pepmenu button:hover{background:var(--soft)}
.pepmenu button[disabled]{color:var(--dim);cursor:default;background:none}
.pepmenu .mi{width:18px;text-align:center;font-size:15px}
#moddlg{width:min(760px,94vw)}
#modseq{font:14px ui-monospace,Menlo,Consolas,monospace;letter-spacing:.12em;line-height:2.2;
  padding:10px 12px;border:1px solid var(--line);border-radius:8px;margin:10px 16px;word-break:break-all}
#modseq .r{cursor:pointer;padding:2px 1px;border-radius:3px}
#modseq .r:hover{background:var(--soft)}
#modseq .r.on,.seq .mod{background:#f6e85a;color:#1b1b1b;border-radius:3px}
#modseq .grp{position:relative;margin-right:.9em}
#modseq .grp::before{content:attr(data-n);position:absolute;right:0;top:-1.3em;font-size:10px;
  color:var(--dim);letter-spacing:0}
.modrow{display:flex;gap:10px;align-items:center;background:var(--soft);border-radius:10px;
  padding:10px 12px;margin:8px 16px}
.modrow .chip{background:#f6e85a;color:#1b1b1b;border-radius:8px;padding:10px 12px;
  font:13px ui-monospace,Menlo,Consolas,monospace;min-width:54px}
.modrow select{flex:1 1 auto;min-width:0}
.modnote{margin:10px 16px 0;padding:10px 12px;border-radius:10px;background:var(--soft);
  font-size:12.5px;color:var(--fg)}
/* Only the middle scrolls, so Save stays in view however many residues
   are modified. */
#moddlg[open]{display:flex;flex-direction:column;max-height:min(90vh,980px)}
.modbody{overflow:auto;flex:1 1 auto;min-height:0;padding-bottom:4px}
/* The structure of each choice, drawn by the server from ADCP's own
   residue definitions: a residue as it sits in a chain, with wavy bonds
   where the chain continues. */
#modviews{margin:2px 16px 10px}
#modviews:empty{display:none}
.mvhead{font-size:11.5px;letter-spacing:.06em;text-transform:uppercase;color:var(--dim);
  margin:10px 0 6px}
.mvgrid{display:grid;grid-template-columns:repeat(auto-fill,minmax(210px,1fr));gap:10px}
.mv{margin:0;border:1px solid var(--line);border-radius:10px;background:var(--card);
  padding:8px 10px 10px;display:flex;flex-direction:column;gap:4px;min-width:0}
.mv figcaption{display:flex;gap:8px;align-items:baseline;font-size:12.5px;line-height:1.35;
  min-width:0}
.mv .pos{flex:none;font:12px ui-monospace,Menlo,Consolas,monospace;background:#f6e85a;
  color:#1b1b1b;border-radius:5px;padding:1px 6px}
.mv .mvname{overflow-wrap:anywhere}
.mv img{display:block;width:100%;height:auto;aspect-ratio:320/240}
.mv .meta{color:var(--dim);font-size:11.5px;overflow-wrap:anywhere}
.mv.nodraw img{display:none}
button.icon{border:0;background:none;color:var(--dim);padding:4px 7px;
  font-size:15px;line-height:1}
button:disabled{opacity:.45;cursor:default}
.why{color:var(--dim);font-size:12.5px}
.bad{color:var(--bad)} .good{color:var(--good)} .warn{color:var(--warn)}

/* ---- the peptide rows --------------------------------------------------- */
/* One row per peptide. Every column has its label above its control and every
   control is the same height, so the three line up on one baseline. */
.ent{display:flex;gap:12px;align-items:flex-start;background:var(--panel);
  border-radius:12px;padding:12px 14px;margin-bottom:10px;text-align:left}
.ent .entnum{color:var(--dim);font-size:12px;font-weight:600;min-width:18px;text-align:right;font-variant-numeric:tabular-nums;user-select:none}
.ent .col{display:flex;flex-direction:column;min-width:0}
.ent .fixed{flex:0 0 130px}
.ent .nm{flex:0 0 140px}
/* The row's own ring: as wide as its longest choice, "head to tail". */
.ent .ring{flex:0 0 132px}
.ent .sqwrap{flex:1 1 auto}
.ent .tools{display:flex;gap:2px;align-self:center}
.ent label{margin-bottom:4px;height:14px;line-height:14px}
/* The controls share one height, so their tops and bottoms line up.
   46px rather than 38: the sequence box carries its residue numbers above the
   letters, so it is taller than a plain input. */
.ent .fixedval,.ent .nm input,.ent .ring select{height:46px}
/* A ring the peptide cannot close is marked on the choice as well as under
   the sequence, since either may be the one to change. */
.ent .ring select.err{border-color:var(--warn)}
.ent .fixedval{background:var(--card);border:1px solid var(--line);
  border-radius:9px;padding:0 11px;display:flex;align-items:center;
  font-size:13.5px;color:var(--fg)}

/* the sequence editor: letters in tens, each group numbered above itself */
.seq{background:var(--card);border:1px solid var(--line);border-radius:9px;
  padding:18px 12px 7px;font-family:ui-monospace,Menlo,monospace;
  font-size:13.5px;line-height:1.9;min-height:46px;box-sizing:border-box;
  overflow-wrap:anywhere;white-space:pre-wrap;text-align:left}
/* An empty box has no residue numbers, so the 18px of headroom they need
   would only push the placeholder down. It is centred instead, and the box
   keeps the same height either way, so nothing jumps when the first letter is
   typed. */
.seq:empty{padding:0 12px;display:flex;align-items:center}
.seq:focus{outline:2px solid var(--accent2);border-color:var(--accent)}
.seq.err{border-color:var(--warn)}
/* Too long is a hard stop, not a warning: ADCP will not build a peptide
   longer than the maximum length, so the box turns the same red a wrong
   letter gets, since both mean the run will not start. */
.seq.toolong{border-color:var(--bad);border-width:2px;
  background:rgba(192,57,43,.05)}
.seqmsg.stop{color:var(--bad);font-weight:600}
/* Every rule here must be closed. A missing closing brace makes a browser
   discard everything after it up to the next terminator, so rules that are
   correct in the file would be absent from the page. tests/test_api.py parses
   the stylesheet to check this. */
.seq:empty::before{content:attr(data-ph);color:var(--dim);
  font-family:-apple-system,"Segoe UI",Roboto,sans-serif;font-size:13px}
/* A gap in the structure, shown where it is: the letters are what the file
   contains, so the row number alone would look like a numbering error. */
.gap{color:var(--warn,#d67a3a);letter-spacing:1px;padding:0 3px;
  cursor:help;font-weight:700}
/* The "?" beside a label. Quiet until wanted, and a real button so it is
   reachable by keyboard and works on a touch screen, which title= does not. */
/* A 15px circle at 10px, the same marker the related tools use. A bare 12px
   "?" read as part of the label rather than as something to press. */
/* button.q rather than .q: ".row button" sets height:38px at specificity
   (0,1,1), and a bare ".q" (0,1,0) would lose to it and stretch every marker
   in a form row. Matching that specificity later in the sheet makes the rule
   apply. Minimum and maximum sizes too, because inline-flex can grow a box to
   its line-height. */
/* No circle. A <button> arrives with a user-agent border, background, padding
   and minimum size, each of which overrides a declared size:
       https://css-tricks.com/overriding-default-button-styles/
       https://ishadeed.com/article/styling-the-good-old-button/

   So there is no box: a "?" in the accent colour at the size of the text
   beside it, which changes colour under the pointer. */
button.q,.row button.q,.listing button.q{
  all:unset;
  font:inherit;font-size:.92em;font-weight:700;
  color:var(--accent);cursor:pointer;
  margin-left:.35em;padding:0 .15em;
  vertical-align:baseline;line-height:inherit;
  /* a flex item defaults to min-width:auto; one glyph either way,
     but stated so it can never be the thing that grows. */
  min-width:0;min-height:0;flex:0 0 auto;
  display:inline}
button.q:hover,button.q:focus{color:var(--fg);text-decoration:underline}
button.q:focus-visible{outline:1px dotted var(--accent);outline-offset:2px}
/* No hover background or border here: the marker has no box in any state. */
/* Both panels head themselves with an h3, set to one size here rather than
   left to the user agent (about 21px). */
#helpbox h3,#issuebox h3{font-size:15px}
/* A link that looks like a button: a download has to be a real <a href>, or
   the browser cannot save it. It looks exactly like the ghost button beside
   it. */
.btnlink{display:inline-block;font:inherit;font-size:13px;
  padding:5px 13px;border-radius:999px;border:1px solid var(--pill-line);
  background:var(--card);color:var(--fg);text-decoration:none;
  cursor:pointer;vertical-align:middle;line-height:1.35}
.btnlink:hover{background:var(--soft);color:var(--fg);text-decoration:none}
.grp{position:relative;display:inline-block;margin-right:1.15em;
  letter-spacing:.19em}
.grp::after{content:attr(data-n);position:absolute;right:.19em;top:-1.35em;
  font-size:10px;letter-spacing:0;color:var(--dim);font-weight:500;
  white-space:nowrap}
.seq .h{color:var(--hydro)} .seq .p{color:var(--polar)}
/* A non-standard or D residue is one residue in the ruler, so it reads as
   one block: its characters close together rather than spaced like the
   letters around it, with the same gap after it as after any letter. */
.seq .mod{letter-spacing:0;padding:0 2px;margin-right:.19em}
.seq .x{color:var(--bad);font-weight:700;background:rgba(192,57,43,.13);
  border-radius:2px}
.seqmsg{font-size:12.5px;margin-top:6px;color:var(--warn)}
.seqlen{font-size:11.5px;color:var(--dim);margin-top:5px}
.bar{display:flex;gap:12px;align-items:center;margin-top:14px;
  justify-content:center}
/* The buttons under the peptide list keep their width, and a long message
   beside them, such as why a list was refused, wraps instead of squeezing
   them into a column of words. */
#step-peptides .bar>button{flex:0 0 auto;white-space:nowrap}
/* On a phone a peptide's row has no room for its columns side by side: the
   type, the name and the ring share a line, the sequence takes the next one
   whole, and the buttons under the list wrap rather than run off the
   screen. */
@media (max-width:600px){
  .ent{flex-wrap:wrap}
  .ent .nm,.ent .fixed,.ent .ring{flex:1 1 110px}
  .ent .sqwrap{flex:1 1 100%;order:3}
  .ent .tools{order:2}
  .bar{flex-wrap:wrap}
}
.centre{display:flex;justify-content:center;margin-top:18px}

/* ---- results and the viewer ----------------------------------------- */
.sitegrid{display:grid;grid-template-columns:minmax(0,1.15fr) minmax(0,1fr);
  gap:14px;margin-top:12px}
@media (max-width:820px){.sitegrid{grid-template-columns:minmax(0,1fr)}}
#viewD{width:100%;height:330px;border:1px solid var(--line);border-radius:10px;
  background:#fff;position:relative}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]) #viewD{background:#0b0e14}}
:root[data-theme=dark] #viewD{background:#0b0e14}
.seq2{font-family:ui-monospace,Menlo,monospace;font-size:12.5px;
  max-height:330px;overflow:auto;border:1px solid var(--line);
  border-radius:10px;padding:9px 11px;background:var(--card);line-height:1.85}
.seq2 .ln{display:flex;gap:9px;align-items:baseline;white-space:pre}
.seq2 .no{color:var(--dim);font-size:10.5px;min-width:3.8em;text-align:right;
  user-select:none}
.seq2 .r{cursor:pointer;border-radius:3px;padding:0 .5px}
/* A sentence in the strip, where the residues would be, is prose: it is not
   set in the strip's monospace. */
.seq2 .why{font-family:-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.seq2 .r:hover{background:var(--accent2)}
.seq2 .r.on{background:var(--accent);color:#fff;font-weight:700}
.chip{display:inline-block;font-size:11.5px;padding:3px 10px;border-radius:999px;
  background:var(--accent2);color:var(--accent);border:1px solid transparent;
  margin:0 6px 6px 0;cursor:pointer;font-family:ui-monospace,Menlo,monospace}
.chip:hover{border-color:var(--accent)}
.stepper{display:flex;align-items:stretch;gap:0}
.stepper input{border-radius:0;text-align:center;
  -moz-appearance:textfield}
.stepper input::-webkit-outer-spin-button,
.stepper input::-webkit-inner-spin-button{-webkit-appearance:none;margin:0}
.stepper .st{border-radius:9px 0 0 9px;padding:8px 11px;font-size:15px;
  line-height:1;color:var(--dim);border-right:0}
.stepper .st:last-child{border-radius:0 9px 9px 0;border-right:1px solid
  var(--line);border-left:0}
table{border-collapse:collapse;width:100%;font-size:13px}
th,td{text-align:left;padding:6px 10px;border-bottom:1px solid var(--line)}
th{font-size:11px;text-transform:uppercase;letter-spacing:.05em;color:var(--dim)}
td.num{text-align:right;font-variant-numeric:tabular-nums}
/* Centred, except names and sequences. Numbers read as a column when their
   digits line up under the heading; a sequence is read letter by letter, so
   it stays left-aligned and monospaced. */
th.c,td.c{text-align:center}
td.c{font-variant-numeric:tabular-nums}
td.seq{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;
  text-align:left;letter-spacing:.02em}
/* A name is read left to right, like td.seq, but it is words rather than
   letters counted one by one, so it is not monospaced. td.nm rather than .nm:
   .nm is already the name box in the peptide entry row. */
td.nm{text-align:left;white-space:nowrap}
/* The results table. It is as wide as its columns need and scrolls sideways
   inside its own box, never squeezing a column to fit the window. Both heading
   rows stay put while it scrolls down, and the three columns that say WHICH
   peptide a row is stay put while it scrolls across. Its cells have classes of
   their own: the sequence editor's .seq box is not a table cell style.
   Separate borders, because a sticky cell under border-collapse leaves its
   border behind when it scrolls. */
.rtable{max-height:64vh;overflow:auto;border:1px solid var(--line);
  border-radius:10px}
.rtable table.rt{width:max-content;min-width:100%;border-collapse:separate;
  border-spacing:0;font-size:12.5px}
.rtable table.rt th,.rtable table.rt td{white-space:nowrap;padding:7px 11px;
  border-bottom:1px solid var(--line);background:var(--card)}
.rtable table.rt thead th{position:sticky;top:0;z-index:3}
.rtable table.rt tr.groups th{font-size:10px;text-align:center;
  padding-top:6px;padding-bottom:5px}
.rtable table.rt .pin{position:sticky;z-index:2}
.rtable table.rt thead th.pin{z-index:5}
.rtable table.rt .gsep{border-left:1px solid var(--line)}
.rtable table.rt .pin-last{box-shadow:inset -1px 0 0 var(--line)}
.rtable table.rt td.rseq{font-family:ui-monospace,SFMono-Regular,Menlo,
  monospace;letter-spacing:.02em;text-align:left}
.rtable table.rt td.rnm{text-align:left}
.rtable table.rt tbody tr.pick:hover td{background:var(--soft)}
.rtable table.rt tbody tr.on td{background:var(--accent2)}
/* Group switches, so a reader after the scores is not made to scroll past the
   whole panel to reach them. */
.colpick{display:flex;gap:13px;flex-wrap:wrap;align-items:center;
  margin:9px 0 6px;font-size:12px;color:var(--dim)}
.colpick label{display:inline-flex;gap:5px;align-items:center;cursor:pointer}
/* Residue by residue, under the viewer: the peptide's residues beside the
   receptor's, each table scrolling in its own box. */
.restitle{font-size:14px;margin:4px 0 8px}
.resgrid{display:grid;grid-template-columns:repeat(auto-fit,minmax(300px,1fr));
  gap:14px}
.reshead{font-size:12px;color:var(--dim);margin:0 0 5px}
.restbl{max-height:360px;overflow:auto;border:1px solid var(--line);
  border-radius:10px}
.restbl table{width:100%;border-collapse:separate;border-spacing:0;
  font-size:12.5px}
.restbl th{position:sticky;top:0;background:var(--card);z-index:1}
.restbl th,.restbl td{padding:6px 10px;white-space:nowrap;
  border-bottom:1px solid var(--line)}
.restbl td.resnm{text-align:left}
.restbl tr.site td{background:var(--soft)}
/* Every residue opens the 3D view on itself and what it touches. The
   pointer and the hover say so; the chosen one stays marked, so it is clear
   which residue the view is showing. */
.restbl tr.resgo{cursor:pointer}
.restbl tr.resgo:hover td{background:var(--soft)}
.restbl tr.resgo.on td{background:var(--soft);
  box-shadow:inset 3px 0 0 var(--accent)}
.restbl tr.resgo.on td+td{box-shadow:none}
.panelnote{margin:9px 0 2px;padding:10px 12px;border:1px solid var(--line);
  border-radius:10px;font-size:13px}
/* The extra actions on a row: smaller, but still plainly buttons. */
.btnlink.sm,button.ghost.sm{font-size:13px;padding:4px 12px;min-height:32px}
/* A label on a row, such as "run" or "site": it says what something is and is
   not pressed, so it has no border and no pointer. */
.tag{display:inline-block;font-size:11.5px;font-weight:600;letter-spacing:.01em;
  padding:0 7px;border-radius:5px;margin-left:6px;
  background:color-mix(in srgb,var(--accent) 14%,var(--card));color:var(--fg);
  vertical-align:1px;cursor:default;line-height:1.5}
.sz{color:var(--dim);font-size:10.5px;margin-left:3px}
#issuebox textarea,#issuebox input{width:100%}
tr.pick{cursor:pointer} tr.pick:hover{background:var(--soft)}
tr.on{background:var(--accent2)}
/* A pill is a control, so it needs a surface to sit on. With only a border
   and no background it had nothing to show on a dark card, where the line
   colour is close to the card's own, and a selected one had no style at all
   outside the viewer bar: .pill.on was written as .vbar .pill.on, so the
   contact pills could be clicked but never looked clicked. */
.pill{display:inline-block;font-size:11px;padding:2px 9px;border-radius:999px;
  border:1px solid var(--pill-line);background:var(--card);color:var(--fg);
  cursor:pointer}
.pill:hover{border-color:var(--dim)}
/* Selected: filled, and outlined as well, so it reads as chosen without
   relying on the fill alone. */
.pill.on{background:var(--accent);border-color:var(--accent);
  color:var(--on-accent);font-weight:650;
  box-shadow:0 0 0 2px var(--accent2)}
.pill.on:hover{border-color:var(--accent)}
/* The interaction types carry their colour as a dot rather than as the
   text colour. Yellow and mid grey on a dark card are unreadable, and the
   dot keeps the colour that matches the drawing without spending it on the
   label. */
.kdot{display:inline-block;width:9px;height:9px;border-radius:50%;
  margin-right:5px;vertical-align:-1px;border:1px solid rgba(0,0,0,.25)}
.list{max-height:240px;overflow:auto;border:1px solid var(--line);
  border-radius:10px}
.list button{display:block;width:100%;text-align:left;border:0;border-radius:0;
  border-bottom:1px solid var(--line);background:none;padding:6px 11px;
  font-size:12.5px}
pre{background:var(--soft);padding:11px;border-radius:9px;overflow:auto;
  font-size:12px;max-height:240px}
.hide{display:none}
/* The bar along the bottom: version, copyright and who does the work. The
   same style as the documentation's, from the shared kit, so the two match
   wherever either changes. */
__CREDITBAR_CSS__
a.btnlink.ghost{font-size:12.5px;padding:5px 12px;color:var(--dim)}
footer{color:var(--dim);font-size:12px;text-align:center;padding:18px 0 4px;
  border-top:1px solid var(--line);margin-top:22px}
/* the viewer's own bar: pills, grouped and labelled */
/* A group does not wrap inside itself, so a label is never separated from the
   buttons it names. */
.vgrp{display:inline-flex;gap:6px;align-items:center;flex:0 0 auto}
/* On a phone a group wider than the screen wraps inside itself, its label
   still first, rather than push the page sideways. */
@media (max-width:600px){.vgrp{flex-wrap:wrap;max-width:100%;flex:0 1 auto}}
.vbar{display:flex;flex-wrap:wrap;gap:6px;align-items:center;
  justify-content:center;
  padding:9px 10px;background:var(--soft);border:1px solid var(--line);
  border-radius:10px;margin-bottom:10px}
.vbar .lb{font-size:10.5px;text-transform:uppercase;letter-spacing:.06em;
  color:var(--dim);font-weight:650;margin-right:1px}
.vbar .pill{font-size:12px;padding:4px 12px;border-radius:999px;
  border:1px solid var(--pill-line);background:var(--card);color:var(--fg);
  cursor:pointer;min-height:28px}
.vbar .pill:hover{border-color:var(--accent)}
.vbar .pill:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.vbar .pill.on{background:var(--accent);border-color:var(--accent);
  color:var(--on-accent);font-weight:600}
.vsep{width:1px;align-self:stretch;background:var(--line);margin:0 4px}
.cpick{display:inline-flex;align-items:center;gap:5px;font-size:11px;
  color:var(--dim)}
.cpick input{width:26px;height:22px;padding:0;border-radius:6px;
  border:1px solid var(--line);background:none}
#viewR,#viewM{width:100%;height:min(58vh,520px);border:1px solid var(--line);
  border-radius:10px;position:relative;background:#fff}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]) :is(#viewR,#viewM){background:#0b0e14}}
:root[data-theme=dark] :is(#viewR,#viewM){background:#0b0e14}
/* the ligand interaction diagram */
#diagwrap h3{font-size:12px;margin:0 0 8px;text-transform:uppercase;
  letter-spacing:.07em;color:var(--dim);font-weight:650;display:flex;
  gap:10px;align-items:center;justify-content:center}
.diagbox{overflow:auto;border:1px solid var(--line);border-radius:10px;
  padding:10px;background:var(--card)}
.diagbox svg{display:block;margin:0 auto}
/* The page says when it is out of date. Each user starts their own server,
   and in a browser an unchanged page and an older page look the same. The
   server can see both the version it started with and the version on disk, so
   it reports the difference. */
.stale{background:var(--warn);color:#241a0c;border-radius:12px;
  padding:12px 16px;margin:0 0 16px;text-align:left;font-size:13px;
  line-height:1.6}
.stale .why2{opacity:.75;font-size:12px}
.stale code{background:rgba(0,0,0,.14);padding:1px 6px;border-radius:5px;
  font-family:ui-monospace,Menlo,monospace}

/* ---- the step guide --------------------------------------------------
   Five numbers across the top, one per step to a run, in the order the
   cards below take them. A number fills green, with a tick in place of the
   digit, when the page's own checks say its step is done, so the state is
   not carried by colour alone; the step to do next is filled in the accent
   colour. Each number is a button that scrolls to its card.

   It stays in view while the page scrolls, on a screen with room for it.
   On a phone it stays at the top: there it would take a tenth of the screen
   from the form for as long as the form is used.

   The line between two numbers is drawn from the centre of one to the
   centre of the next, short of both circles: 50% is the centre of a step,
   and 20px is the circle's radius (15px) and a gap. */
/* The band behind the strip is 6px wider than the strip on each side, so
   the 3px ring of a lit card scrolling under it does not show beside it. */
.guide{position:sticky;top:0;z-index:30;margin:-8px -6px 12px;
  padding:8px 6px;background:var(--bg)}
.guide ol{list-style:none;margin:0;padding:10px 12px 9px;display:flex;
  background:var(--card);border:1px solid var(--line);border-radius:14px}
.guide li{flex:1 1 0;min-width:0;position:relative;text-align:center}
.guide li+li::before{content:"";position:absolute;top:15px;height:2px;
  left:calc(-50% + 20px);right:calc(50% + 20px);border-radius:1px;
  background:var(--line)}
.guide li.done+li::before{background:var(--good)}
/* No user agent box: the circle is the button's face. */
.guide button{all:unset;box-sizing:border-box;width:100%;cursor:pointer;
  display:flex;flex-direction:column;align-items:center;gap:4px;
  padding:0 2px 1px;border-radius:10px}
.guide button:hover .gnum{border-color:var(--accent)}
.guide button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
.gnum{width:30px;height:30px;border-radius:50%;flex:0 0 auto;
  display:inline-flex;align-items:center;justify-content:center;
  font-size:13px;font-weight:700;font-variant-numeric:tabular-nums;
  border:2px solid var(--pill-line);background:var(--card);color:var(--dim);
  transition:background-color .25s,border-color .25s,color .25s,
    box-shadow .25s}
.gtick{display:none;width:17px;height:17px;fill:currentColor}
.glab{font-size:12px;line-height:1.25;color:var(--dim);
  overflow-wrap:anywhere}
.guide li.now .gnum{background:var(--accent);border-color:var(--accent);
  color:var(--on-accent);box-shadow:0 0 0 4px var(--accent2)}
.guide li.now .glab{color:var(--fg);font-weight:650}
.guide li.done .gnum{background:var(--good);border-color:var(--good);
  color:var(--on-accent)}
.guide li.done .gn{display:none}
.guide li.done .gtick{display:block}
.guide li.done .glab{color:var(--good);font-weight:600}
/* Busy: a turning arc round the number while ADapt starts after a structure
   is loaded, and two words under the label. Nothing is blocked by it. */
.gnum{position:relative}
.guide li.busy .gnum::after{content:"";position:absolute;inset:-6px;
  border-radius:50%;border:2px solid transparent;border-top-color:var(--accent);
  border-right-color:var(--accent);animation:ld-turn 1s linear infinite}
.guide li.busy .glab::after{content:"ADapt starting";display:block;
  font-size:10.5px;font-weight:600;color:var(--accent)}
/* Words for a screen reader, which cannot see the fill or the tick. */
.vh{position:absolute;width:1px;height:1px;margin:-1px;padding:0;
  overflow:hidden;clip:rect(0 0 0 0);white-space:nowrap;border:0}
@media (max-width:600px),(max-height:560px){.guide{position:static}}
@media (max-width:600px){
  .guide ol{padding:8px 4px 7px}
  .gnum{width:26px;height:26px;font-size:12px}
  .gtick{width:15px;height:15px}
  .glab{font-size:10.5px}
  .guide li+li::before{top:13px;left:calc(-50% + 17px);right:calc(50% + 17px)}
}

/* ---- the cards follow the guide --------------------------------------
   The card of the step to do next is lit. The cards after it are grey, and
   a finished card is green with a tick where its number was. The number in
   a card's heading is a counter over the cards that are steps, so it cannot
   disagree with the strip's order.

   Grey is a cue and nothing more. Nothing in a grey card is disabled or
   made unclickable: a target file fills two steps at once, and a reused
   run fills four, and neither should have to be walked through. A grey
   card comes back to full colour under the pointer or the keyboard.

   scroll-margin-top keeps a card's heading clear of the strip when a step
   number scrolls to it. */
main{counter-reset:step}
.card[data-step]{counter-increment:step;scroll-margin-top:96px;
  transition:border-color .25s,box-shadow .25s}
.card[data-step]>*{transition:opacity .25s,filter .25s}
.card[data-step]>h2::before{content:counter(step);display:inline-flex;
  align-items:center;justify-content:center;width:19px;height:19px;
  margin-right:8px;border-radius:50%;border:1.5px solid currentColor;
  box-sizing:border-box;font-size:10.5px;letter-spacing:0;
  vertical-align:1px}
.card.step-now{border-color:var(--accent);box-shadow:0 0 0 3px var(--accent2)}
.card.step-now>h2{box-shadow:inset 0 3px 0 var(--accent);color:var(--fg)}
.card.step-now>h2::before{background:var(--accent);border-color:var(--accent);
  color:var(--on-accent)}
.card.step-done{border-color:var(--good);
  border-color:color-mix(in srgb,var(--good) 55%,var(--line))}
.card.step-done>h2{box-shadow:inset 0 2px 0 var(--good);
  background:color-mix(in srgb,var(--good) 9%,var(--ribbon))}
.card.step-done>h2::before{content:"\2713";background:var(--good);
  border-color:var(--good);color:var(--on-accent);font-size:11.5px}
.card.step-later>h2{box-shadow:inset 0 2px 0 var(--pill-line);
  color:var(--dim)}
.card.step-later>:not(h2){filter:grayscale(1);opacity:.8}
.card.step-later:hover>:not(h2),
.card.step-later:focus-within>:not(h2){filter:none;opacity:1}
@media (prefers-reduced-motion:reduce){
  .gnum,.card[data-step],.card[data-step]>*{transition:none}
  .guide li.busy .gnum::after{animation:none;border-style:dashed;
    border-color:var(--accent)}
}

/* ---- the runs: two views ---------------------------------------------- */
.tabs{display:flex;gap:6px;justify-content:center;flex-wrap:wrap;
  margin:0 0 12px}
.tabs button{font-size:13px;padding:6px 16px;color:var(--fg);min-height:34px}
.tabs button[aria-selected=true]{background:var(--accent);
  border-color:var(--accent);color:var(--on-accent);font-weight:600}
/* A table of runs is wider than a phone, and scrolls in its own box rather
   than pushing the page sideways. Its row buttons wrap onto a second line
   before the table grows past its card, and a date is never broken at its
   hyphens. */
#runs{overflow-x:auto}
#runs td.acts{white-space:normal;line-height:2.1;min-width:16em}
#runs td.why.c{white-space:nowrap}
/* A path has no spaces to break at. The runs folder beside Browse may break
   anywhere; an open run's title and folder break only where a word cannot
   fit, so "completed" beside a job is never split. */
#runsroot{overflow-wrap:anywhere;min-width:0}
#rundetail{overflow-wrap:break-word;min-width:0}
@media (max-width:600px){
  .jobrow .jname{flex-basis:64px}
  .jobrow .jid{flex-basis:76px}
}
/* The Folders view: where it is, as a trail back to the runs folder, then
   the folder's contents. Paths are long, so every part may wrap. */
.fcrumbs{display:flex;flex-wrap:wrap;gap:2px 3px;align-items:baseline;
  font-size:12.5px;color:var(--dim);margin:0 0 9px;text-align:left;
  overflow-wrap:anywhere}
.fcrumbs button{all:unset;cursor:pointer;color:var(--accent);
  padding:1px 4px;border-radius:5px}
.fcrumbs button:hover{background:var(--soft);text-decoration:underline}
.fcrumbs button:focus-visible{outline:2px solid var(--accent)}
.fcrumbs b{color:var(--fg);padding:1px 4px}
.fhere{margin:0 0 10px;text-align:left}
.ftable{overflow:auto;max-height:62vh;border:1px solid var(--line);
  border-radius:10px}
.ftable table{font-size:12.5px}
.ftable th{position:sticky;top:0;z-index:1;background:var(--card)}
.ftable td{vertical-align:middle}
.ftable td.fname{overflow-wrap:anywhere;min-width:11em}
.ftable .fgo{all:unset;cursor:pointer;display:inline-flex;gap:7px;
  align-items:center;color:var(--fg)}
.ftable .fgo:hover{color:var(--accent);text-decoration:underline}
.ftable .fgo:focus-visible{outline:2px solid var(--accent);border-radius:4px}
.ftable svg.fic{width:15px;height:15px;flex:0 0 auto;fill:var(--dim)}
.ftable svg.fic.dir{fill:var(--accent)}
.ftable td.fsz,.ftable td.fmod{white-space:nowrap;color:var(--dim);
  font-variant-numeric:tabular-nums}
.ftable td.facts{white-space:nowrap;text-align:right}
.ftable .pill{margin-left:6px;vertical-align:1px}
.fmeta{display:none;color:var(--dim);font-size:11px;margin-top:2px}
@media (max-width:600px){
  .ftable th.fmod,.ftable td.fmod{display:none}
  .fmeta{display:block}
  .ftable td.facts{white-space:normal}
  .ftable td.facts button{margin:2px 0}
}
/* A file from the runs folder, shown in the page: a structure in the same
   3D viewer the page uses, anything else as text. */
#filedlg{width:min(980px,94vw)}
#filedlg[open]{display:flex;flex-direction:column;max-height:92vh}
#file3d{height:min(62vh,560px);position:relative;background:#fff;
  flex:0 0 auto}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]) #file3d{background:#0b0e14}}
:root[data-theme=dark] #file3d{background:#0b0e14}
#filetext{margin:0;border-radius:0;max-height:none;flex:1 1 auto;
  min-height:0;white-space:pre;overflow:auto}
#filedlg .dlg-f{flex-wrap:wrap}

/* ---- waiting ------------------------------------------------------------
   One component for every wait longer than a moment (app.js, Wait), in two
   sizes. The stage fills the space something will appear in: a quiet outline
   of what is coming, and over it a card with a small molecule drawing itself,
   a title, why to wait, the real steps ticking as they happen, the time so
   far, and a line that it is fine to carry on. The inline size is one line,
   for a wait beside a button. Both turn into an error with Try again.

   The colours are the page's tokens, so both themes follow. With less motion
   asked for, the shimmer and the drawing stop and everything is shown
   still. */
.ld-host{position:relative}
.ld-stage{position:absolute;inset:0;z-index:4;display:flex;align-items:flex-start;
  justify-content:center;padding:clamp(18px,7%,64px) 16px 18px;background:var(--card);
  overflow:auto;transition:opacity .5s ease}
.ld-stage.ld-leave{opacity:0;pointer-events:none}
/* ADapt's layout, sketched: its tree beside its view. */
.ld-skel{position:absolute;inset:0;display:grid;
  grid-template-columns:minmax(0,34%) minmax(0,1fr);gap:14px;padding:16px;
  opacity:.9}
.ld-sk-tree{display:flex;flex-direction:column;gap:11px;padding-top:6px}
.ld-sk-tree i{display:block;height:11px;border-radius:6px;width:78%}
.ld-sk-tree i:nth-child(3n+2){width:62%;margin-left:16px}
.ld-sk-tree i:nth-child(3n){width:48%;margin-left:32px}
.ld-sk-tree i:first-child{height:15px;width:44%;margin-bottom:6px}
.ld-sk-view{border-radius:12px;position:relative;overflow:hidden}
.ld-sk-tree i,.ld-sk-view{background:linear-gradient(100deg,var(--soft) 30%,
  color-mix(in srgb,var(--line) 70%,var(--soft)) 50%,var(--soft) 70%);
  background-size:300% 100%;animation:ld-shimmer 1.8s linear infinite}
.ld-sk-view b{position:absolute;left:50%;top:50%;width:46%;aspect-ratio:1;
  transform:translate(-50%,-50%);border-radius:50%;
  border:2px dashed color-mix(in srgb,var(--line) 80%,transparent)}
@keyframes ld-shimmer{0%{background-position:100% 0}100%{background-position:0 0}}
.ld-card{position:relative;width:min(440px,100%);background:var(--card);
  border:1px solid var(--line);border-radius:16px;padding:22px 24px 18px;
  box-shadow:0 18px 50px rgba(15,25,45,.16),0 2px 6px rgba(15,25,45,.06);
  text-align:center;color:var(--fg)}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]) .ld-card{
  box-shadow:0 18px 50px rgba(0,0,0,.55)}}
:root[data-theme=dark] .ld-card{box-shadow:0 18px 50px rgba(0,0,0,.55)}
.ld-title{margin:6px 0 4px;font-size:18px;font-weight:650;letter-spacing:-.01em}
.ld-why{margin:0 auto 14px;color:var(--dim);font-size:13.5px;line-height:1.5;
  max-width:34em}
/* The mark: a ring of six atoms whose bonds draw in one after another, then
   an aromatic ring inside, and one atom going round. */
.ld-mark{width:62px;height:62px;display:block;margin:0 auto;overflow:visible}
.ld-mark .ld-bond,.ld-mark .ld-arom{fill:none;stroke:var(--accent);
  stroke-width:2.6;stroke-linecap:round;stroke-dasharray:1;stroke-dashoffset:0}
.ld-mark .ld-arom{stroke-width:1.8;opacity:.55}
.ld-mark .ld-atom{fill:var(--accent);stroke:var(--card);stroke-width:1.6;
  transform-box:fill-box;transform-origin:center}
.ld-mark .ld-orbit{fill:none;stroke:var(--accent2);stroke-width:1.4}
.ld-mark .ld-base{fill:none;stroke:color-mix(in srgb,var(--accent) 22%,var(--card));
  stroke-width:2.6;stroke-linejoin:round}
.ld-mark .ld-e{fill:var(--good)}
.ld-mark .ld-spin{transform-origin:32px 32px;animation:ld-turn 3s linear infinite}
.ld-mark .ld-bond{animation:ld-draw 2.8s ease-in-out infinite both}
.ld-mark .ld-atom{animation:ld-pop 2.8s ease-in-out infinite both}
.ld-mark .ld-arom{animation:ld-ring 2.8s ease-in-out infinite both}
.ld-mark .b2,.ld-mark .a2{animation-delay:.16s}
.ld-mark .b3,.ld-mark .a3{animation-delay:.32s}
.ld-mark .b4,.ld-mark .a4{animation-delay:.48s}
.ld-mark .b5,.ld-mark .a5{animation-delay:.64s}
.ld-mark .b6,.ld-mark .a6{animation-delay:.8s}
@keyframes ld-draw{0%{stroke-dashoffset:1;opacity:.35}22%{stroke-dashoffset:0;opacity:1}
  78%{stroke-dashoffset:0;opacity:1}100%{stroke-dashoffset:0;opacity:.35}}
@keyframes ld-pop{0%{transform:scale(.6);opacity:.5}16%{transform:scale(1.18);opacity:1}
  26%{transform:scale(1)}78%{transform:scale(1);opacity:1}100%{transform:scale(.6);opacity:.5}}
@keyframes ld-ring{0%,30%{stroke-dashoffset:1}55%,80%{stroke-dashoffset:0}100%{stroke-dashoffset:1}}
@keyframes ld-turn{to{transform:rotate(360deg)}}
/* The steps, each with its state in its icon and its words: waiting, under
   way, done (a tick, and when), or where it stopped (a cross). */
.ld-steps{list-style:none;margin:0 auto;padding:12px 14px;max-width:360px;
  text-align:left;background:color-mix(in srgb,var(--soft) 55%,var(--card));
  border:1px solid var(--line);
  border-radius:12px;display:flex;flex-direction:column;gap:9px}
.ld-step{display:flex;align-items:center;gap:10px;font-size:13.5px;
  color:var(--dim);line-height:1.35}
.ld-ic{flex:0 0 auto;width:18px;height:18px;border-radius:50%;
  border:2px solid var(--pill-line);display:inline-flex;align-items:center;
  justify-content:center;box-sizing:border-box}
.ld-ic svg{width:12px;height:12px;fill:var(--on-accent);display:none}
.ld-step .ld-lab{flex:1 1 auto;min-width:0}
.ld-at{flex:0 0 auto;font-size:12px;color:var(--dim);
  font-variant-numeric:tabular-nums}
.ld-step.ld-now{color:var(--fg);font-weight:600}
.ld-step.ld-now .ld-ic{border-color:var(--accent2);border-top-color:var(--accent);
  border-right-color:var(--accent);animation:ld-turn .9s linear infinite}
.ld-step.ld-ok{color:var(--fg)}
.ld-step.ld-ok .ld-ic{background:var(--good);border-color:var(--good)}
.ld-step.ld-ok .ld-ic .ld-tick,.ld-step.ld-no .ld-ic .ld-x{display:block}
.ld-step.ld-no{color:var(--bad);font-weight:600}
.ld-step.ld-no .ld-ic{background:var(--bad);border-color:var(--bad)}
.ld-meta{margin:12px 0 0;font-size:12.5px;color:var(--dim)}
.ld-time,.ld-at{white-space:nowrap}
.ld-time{font-variant-numeric:tabular-nums;font-weight:600;color:var(--fg)}
.ld-slow{margin:10px auto 0;max-width:360px;text-align:left;font-size:13px;
  line-height:1.45;padding:9px 12px;border-radius:10px;color:var(--fg);
  background:color-mix(in srgb,var(--accent2) 75%,var(--card));
  border-left:3px solid var(--accent)}
.ld-keep{margin:10px 0 0;font-size:12.5px;color:var(--dim)}
.ld-keep svg{width:13px;height:13px;fill:currentColor;vertical-align:-2px;
  margin-right:4px}
/* Gone wrong: the card says so in red, in plain words, with Try again. */
.ld-errbox{margin:2px 0 14px}
.ld-errtext{margin:0 auto 6px;max-width:34em;font-size:13.5px;line-height:1.5}
.ld-detail{margin:0 auto 12px;max-width:36em;font-size:12px;color:var(--dim);
  overflow-wrap:anywhere;text-align:left;padding:7px 10px;border-radius:8px;
  background:color-mix(in srgb,var(--soft) 55%,var(--card));
  border:1px solid var(--line)}
.ld-card.ld-bad{border-color:color-mix(in srgb,var(--bad) 55%,var(--line))}
.ld-card.ld-bad .ld-title{color:var(--bad)}
.ld-badmark{display:none;width:52px;height:52px;margin:0 auto;fill:var(--bad)}
.ld-card.ld-bad .ld-mark{display:none}
.ld-card.ld-bad .ld-badmark{display:block}
.ld-stage.ld-halt .ld-sk-tree i,.ld-stage.ld-halt .ld-sk-view{animation:none;
  opacity:.6}
.ld-retry:focus-visible{outline:3px solid var(--accent);outline-offset:3px}
.ld-card.ld-fine .ld-title{color:var(--good)}
/* Ready: the ring stands complete and still. */
.ld-card.ld-fine .ld-mark *{animation:none}
.ld-card.ld-fine .ld-mark .ld-bond,.ld-card.ld-fine .ld-mark .ld-arom{
  stroke-dashoffset:0;opacity:1}
.ld-card.ld-fine .ld-mark .ld-arom{opacity:.55}
/* A whole tab of its own (Open full screen while ADapt starts). */
.ld-page{position:fixed;inset:0;background:var(--bg)}
.ld-page .ld-stage{background:var(--bg);align-items:center}
/* The inline size: one line beside a button, for a shorter wait. */
.ld-inline{display:inline-flex;align-items:center;gap:8px;max-width:100%;
  padding:5px 13px 5px 9px;border-radius:14px;background:var(--accent2);
  border:1px solid color-mix(in srgb,var(--accent) 35%,transparent);
  color:var(--fg);font-size:13px;font-weight:600;line-height:1.35;
  text-align:left;vertical-align:middle}
.ld-inline .ld-time{font-weight:400;color:var(--dim);flex:0 0 auto}
.ld-inline .ld-lab{min-width:0}
.ld-mini{width:18px;height:18px;flex:0 0 auto;overflow:visible}
.ld-mini circle{fill:var(--accent);animation:ld-blink 1.2s linear infinite}
.ld-mini .m2{animation-delay:.2s}.ld-mini .m3{animation-delay:.4s}
.ld-mini .m4{animation-delay:.6s}.ld-mini .m5{animation-delay:.8s}
.ld-mini .m6{animation-delay:1s}
.ld-mini path{fill:none;stroke:var(--accent);stroke-width:1.2;opacity:.35}
@keyframes ld-blink{0%,100%{opacity:.2}30%{opacity:1}}
.ld-inline.ld-bad{background:color-mix(in srgb,var(--bad) 9%,var(--card));
  border-color:color-mix(in srgb,var(--bad) 45%,transparent);flex-wrap:wrap;
  font-weight:400}
.ld-inline.ld-bad .ld-lab{font-weight:600;color:var(--bad)}
.ld-inline .ld-x2{width:16px;height:16px;fill:var(--bad);flex:0 0 auto;
  display:none}
.ld-inline.ld-bad .ld-x2{display:block}
.ld-inline.ld-bad .ld-mini{display:none}
.ld-inline .ld-errtext{color:var(--fg);flex:1 1 220px;min-width:0}
.ld-inline .ld-retry{margin-left:2px;min-height:32px;font-size:13px;
  padding:3px 13px}
/* The heading's badge and the strip's step while ADapt starts. */
.cardbusy{display:inline-flex;align-items:center;gap:6px;margin-left:10px;
  padding:2px 10px 2px 6px;border-radius:999px;background:var(--accent2);
  color:var(--fg);font-size:11.5px;font-weight:650;letter-spacing:0;
  text-transform:none;vertical-align:1px;border:0;cursor:pointer;
  font-family:inherit;line-height:1.5}
.cardbusy .ld-mini{width:14px;height:14px}
.cardbusy.ld-bad{background:color-mix(in srgb,var(--bad) 12%,var(--card))}
.cardbusy .ld-x2{display:block;width:14px;height:14px;fill:var(--bad)}
.cardbusy:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
@media (prefers-reduced-motion:reduce){
  .ld-sk-tree i,.ld-sk-view,.ld-mark .ld-bond,.ld-mark .ld-atom,
  .ld-mark .ld-arom,.ld-mark .ld-spin,.ld-step.ld-now .ld-ic,
  .ld-mini circle{animation:none}
  .ld-mark .ld-bond,.ld-mark .ld-arom{stroke-dashoffset:0}
  .ld-step.ld-now .ld-ic{border-color:var(--accent)}
  .ld-mini circle{opacity:.75}
  .ld-stage{transition:none}
  #adaptframe{transition:none}
}
@media (max-width:600px){
  .ld-skel{grid-template-columns:minmax(0,1fr)}
  .ld-sk-tree{display:none}
  .ld-card{padding:18px 16px 14px}
  .ld-title{font-size:16.5px}
}

/* ---- ADapt's receptor, as a status card --------------------------------
   The headline is what docking will use; the facts are chips; the long
   details fold away; and the one action the state calls for is a filled
   button, with the way back beside it as a secondary. */
.aucard{margin-top:10px;background:var(--card);border:1px solid var(--line);
  border-left:4px solid var(--accent);border-radius:12px;padding:12px 14px;
  text-align:left}
.aucard.au-current{border-left-color:var(--good)}
.aucard.au-newer,.aucard.au-pending{border-left-color:var(--warn)}
.auhead{display:flex;align-items:flex-start;gap:11px;flex-wrap:wrap}
.auicon{flex:0 0 auto;width:34px;height:34px;border-radius:10px;
  background:var(--accent2);color:var(--accent);display:inline-flex;
  align-items:center;justify-content:center}
.auicon svg{width:21px;height:21px;fill:currentColor}
.aumain{flex:1 1 220px;min-width:0}
.autitle{font-size:15px;font-weight:650;color:var(--fg);line-height:1.3}
.ausub{font-size:12.5px;color:var(--dim);margin-top:2px;overflow-wrap:anywhere}
.austate{flex:0 0 auto;display:inline-flex;align-items:center;gap:6px;
  font-size:12.5px;font-weight:650;padding:4px 11px;border-radius:999px;
  align-self:center;color:var(--fg)}
.austate svg{width:15px;height:15px;flex:0 0 auto}
.austate:empty{display:none}
.austate.ok{background:color-mix(in srgb,var(--good) 14%,var(--card))}
.austate.ok svg{fill:var(--good)}
.austate.todo{background:color-mix(in srgb,var(--warn) 16%,var(--card))}
.austate.todo svg{fill:var(--warn)}
.austate.info{background:var(--accent2)}
.auchips{display:flex;flex-wrap:wrap;gap:6px;margin-top:10px}
.auchips:empty{display:none}
.auchip{display:inline-flex;align-items:baseline;gap:6px;font-size:12.5px;
  padding:3px 10px;border-radius:8px;background:var(--card);
  border:1px solid var(--pill-line);color:var(--fg);max-width:100%;
  overflow-wrap:anywhere}
.auchip b{font-weight:600;color:var(--dim);font-size:11px;
  text-transform:uppercase;letter-spacing:.04em}
.aunote{margin-top:10px;font-size:13px;line-height:1.5;color:var(--fg)}
.aunote:empty{display:none}
.auacts{display:flex;flex-wrap:wrap;gap:8px 10px;margin-top:12px;
  align-items:center}
.aumore{margin-top:10px;font-size:12.5px}
.aumore summary{cursor:pointer;color:var(--accent);font-weight:600;
  display:inline-flex;align-items:center;gap:5px;padding:2px 0;list-style:none}
.aumore summary::-webkit-details-marker{display:none}
.aumore summary::before{content:"";width:0;height:0;border-style:solid;
  border-width:4.5px 0 4.5px 7px;border-color:transparent transparent transparent
  currentColor;transition:transform .15s}
.aumore[open] summary::before{transform:rotate(90deg)}
@media (prefers-reduced-motion:reduce){.aumore summary::before{transition:none}}
.aumore summary:focus-visible{outline:2px solid var(--accent);outline-offset:2px;
  border-radius:4px}
.aumore .why{margin-top:6px;line-height:1.6}

/* ---- column widths -----------------------------------------------------
   Every table of runs, results and residues can have its columns resized
   (app.js, colResize): a grip on the right edge of each heading, dragged with
   a mouse, a pen or a finger, moved with the arrow keys once focused, or
   double clicked to fit the column to its widest entry. Once a table has a
   resized column its layout is fixed, so the other columns keep their
   widths, and text too long for its column ends in an ellipsis with the
   whole of it on hover. */
/* No weight of its own, so a heading that sticks keeps sticking. */
:where(th.colsz){position:relative}
.colgrip{position:absolute;top:0;right:0;bottom:0;width:11px;z-index:6;
  cursor:col-resize;touch-action:none;display:flex;justify-content:center;
  align-items:center;outline:none}
.colgrip::after{content:"";width:2px;height:62%;min-height:12px;
  border-radius:1px;background:transparent;transition:background-color .15s}
th:hover>.colgrip::after{background:var(--pill-line)}
.colgrip:hover::after,.colgrip.on::after{background:var(--accent);width:3px}
.colgrip:focus-visible::after{background:var(--accent);width:3px}
.colgrip:focus-visible{box-shadow:inset 0 0 0 2px var(--accent);border-radius:3px}
table.colfixed{table-layout:fixed}
table.colfixed th,table.colfixed td{overflow:hidden;text-overflow:ellipsis}
table.colfixed td.fname{white-space:nowrap;overflow-wrap:normal}
table.colfixed td.fname .fgo{max-width:100%;min-width:0}
table.colfixed td.fname .fgo span{overflow:hidden;text-overflow:ellipsis;
  white-space:nowrap;min-width:0}
table.colfixed td.acts,table.colfixed td.facts{white-space:normal}
body.colmoving,body.colmoving *{cursor:col-resize !important;
  user-select:none !important}
.colreset{margin:6px 0 0;display:block;margin-left:auto}
.colreset.hide{display:none}
@media (prefers-reduced-motion:reduce){.colgrip::after{transition:none}}

/* ---- a phone ----------------------------------------------------------
   The residue strip lays a row of thirty residues on one line, wider than a
   phone. There each row wraps at its groups of ten, beside its number,
   instead of running out of its box on the right. */
@media (max-width:600px){
  .seq2 .ln>span:last-child{white-space:normal;min-width:0}
  /* The Docs button sits in the corner; on a phone the title reaches the
     corner, so the header starts lower and the two never overlap. */
  header{padding-top:46px}
}
</style>

__STALE__
<header>
  <a class="ghost btnlink docscorner" id=docsbtn href="/docs/index.html" target=_blank
     rel=noopener title="How to use this page, with worked examples">Docs</a>
  <div class=brandblock>
    <h1>GUI to ADCP Docking</h1>
    <!-- What belongs to whom, in the header: ADCP and ADFRsuite are the
         Scripps Research Institute's, and this is a GUI and workflow that
         drives them. A screenshot carries the header with it. -->
    <div class=sub>A graphical interface to
      <a href="https://ccsb.scripps.edu/adcp/" target="_blank"
         rel="noopener">ADCP</a>, a product of the
      <b>Scripps Research Institute</b>. <b>Powered by HOPE Lab.</b></div>
    <!-- What a paper has to cite, where it is seen before anything is run.
         The sentence is cite.must_cite_sentence(), read from the same
         entries as the How to cite list, so the two cannot name different
         papers. The link opens that list. -->
    <div class="sub mustcite">__MUST_CITE__
      <a href="#" id=citehead>How to cite</a></div>
    <div class=sub>Authors: Aadhil Haq &middot; PI Dr. Sandun Fernando</div>
    <div class=sub>
      <a href="mailto:sfernando@tamu.edu">sfernando@tamu.edu</a> &middot;
      <a href="mailto:aadhil.haq@tamu.edu">aadhil.haq@tamu.edu</a>
    </div>
    <div class=headbar>
      <!-- How to cite, at the top, where it is looked for when a methods
           section is being written. The engine and the preparation tools are
           the citations most often left out. -->
      <button class=ghost id=citebtn>How to cite</button>
      <button class=ghost id=feedback>Send feedback</button>
      <button class=ghost id=themebtn aria-pressed=false
              title="Switch between light and dark"><svg class=sun viewBox="0 0 24 24" aria-hidden=true><path d="M12 17a5 5 0 1 1 0-10 5 5 0 0 1 0 10zm0-13a1 1 0 0 1-1-1V2a1 1 0 0 1 2 0v1a1 1 0 0 1-1 1zm0 18a1 1 0 0 1-1-1v-1a1 1 0 0 1 2 0v1a1 1 0 0 1-1 1zM4 13H3a1 1 0 0 1 0-2h1a1 1 0 0 1 0 2zm17 0h-1a1 1 0 0 1 0-2h1a1 1 0 0 1 0 2zM6.3 7.7 5.6 7a1 1 0 0 1 1.4-1.4l.7.7a1 1 0 0 1-1.4 1.4zm11 11-.7-.7a1 1 0 0 1 1.4-1.4l.7.7a1 1 0 0 1-1.4 1.4zm.7-11a1 1 0 0 1-.7-1.7l.7-.7A1 1 0 0 1 19.4 7l-.7.7a1 1 0 0 1-.7.3zm-11 11a1 1 0 0 1-.7-1.7l.7-.7a1 1 0 0 1 1.4 1.4l-.7.7a1 1 0 0 1-.7.3z"/></svg><svg class=moon viewBox="0 0 24 24" aria-hidden=true><path d="M12.3 22a10 10 0 0 1-1.4-19.9 1 1 0 0 1 1 1.6 7.5 7.5 0 0 0 9.4 11.3 1 1 0 0 1 1.3 1.4A10 10 0 0 1 12.3 22z"/></svg><span id=themelab>Dark</span></button>
      <span class=sub id=toolnote></span>
    </div>
  </div>
</header>

<main>
<!-- The steps to a run, in the order the cards below take them. app.js
     (guide) marks each done, to do now or still to come, from the checks the
     page already makes, and lights the cards to match. Each number scrolls
     to its card. The tick drawn in a finished number is one path, inline,
     so it needs no font that may be missing. -->
<!-- The strip and its five cards share one box, so the strip follows the
     page down only while one of those cards is in view, and leaves with the
     last of them instead of covering the runs and their results. -->
<div class=stepsflow>
<nav class=guide id=guide aria-label="Steps to a run">
  <ol>
    <li data-step=receptor><button type=button data-goto=receptor><span
      class=gnum aria-hidden=true><span class=gn>1</span><svg class=gtick
      viewBox="0 0 24 24"><path d="M9.5 16.6 4.9 12l-1.6 1.6 6.2 6.2L21 8.3l-1.6-1.6z"/></svg></span><span
      class=glab>Receptor</span><span class="vh gstate"></span></button></li>
    <li data-step=site><button type=button data-goto=site><span
      class=gnum aria-hidden=true><span class=gn>2</span><svg class=gtick
      viewBox="0 0 24 24"><path d="M9.5 16.6 4.9 12l-1.6 1.6 6.2 6.2L21 8.3l-1.6-1.6z"/></svg></span><span
      class=glab>Site</span><span class="vh gstate"></span></button></li>
    <li data-step=peptides><button type=button data-goto=peptides><span
      class=gnum aria-hidden=true><span class=gn>3</span><svg class=gtick
      viewBox="0 0 24 24"><path d="M9.5 16.6 4.9 12l-1.6 1.6 6.2 6.2L21 8.3l-1.6-1.6z"/></svg></span><span
      class=glab>Peptides</span><span class="vh gstate"></span></button></li>
    <li data-step=docking><button type=button data-goto=docking><span
      class=gnum aria-hidden=true><span class=gn>4</span><svg class=gtick
      viewBox="0 0 24 24"><path d="M9.5 16.6 4.9 12l-1.6 1.6 6.2 6.2L21 8.3l-1.6-1.6z"/></svg></span><span
      class=glab>Docking</span><span class="vh gstate"></span></button></li>
    <li data-step=submit><button type=button data-goto=submit><span
      class=gnum aria-hidden=true><span class=gn>5</span><svg class=gtick
      viewBox="0 0 24 24"><path d="M9.5 16.6 4.9 12l-1.6 1.6 6.2 6.2L21 8.3l-1.6-1.6z"/></svg></span><span
      class=glab>Submit</span><span class="vh gstate"></span></button></li>
  </ol>
</nav>

<div class=card id=step-receptor data-step=receptor>
  <h2>Receptor</h2>
  <p class=hint>The protein the peptides will be docked into, and how it is
    prepared before anything is docked into it.</p>
  <div class=row>
    <div style="flex:0 0 150px">
      <label for=reccode>PDB id</label>
      <input type=text id=reccode placeholder="4CMM" maxlength=4
             style="text-transform:uppercase"></div>
    <div style="flex:0 0 auto"><button class=primary id=recfetch>Fetch from RCSB</button></div>
    <div class=why style="flex:0 0 auto;padding-bottom:9px">or</div>
    <div class=fieldwrap style="flex:3 1 280px">
      <label for=recpath>a file already on this machine</label>
      <input type=text id=recpath
             placeholder="/path/to/receptor.pdb  (or .cif, .mmcif, .ent, .gz)"></div>
    <div style="flex:0 0 auto"><button class=secondary id=recload>Load</button></div>
    <div style="flex:0 0 auto"><button class=secondary id=recbrowse>Browse&hellip;</button></div>
    <!-- From the user's own computer. "Browse" walks the cluster's
         filesystem; this reads a file on the user's machine and sends it
         through the connection that is already open. -->
    <!-- No .pdbqt. A PDBQT is a receptor someone has already prepared, and
         this workflow prepares its receptor from a PDB or mmCIF; nothing
         here builds maps from a PDBQT as it stands. The way to dock into a
         receptor prepared elsewhere is a target file, below. -->
    <div style="flex:0 0 auto">
      <button class=secondary id=recupload>From this computer&hellip;</button>
      <input type=file id=recfile style="display:none"
             accept=".pdb,.ent,.cif,.mmcif,.txt"></div>
  </div>
  <div id=recinfo class=why style="margin-top:9px"></div>
  <div id=chainbox class="hide" style="margin-top:10px"></div>
  <!-- The preparation is the receptor's, so its choices are here rather than
       among the docking settings, beside the chains they are applied to and
       in the order prepare.run applies them: heteroatoms are dropped, missing
       atoms are rebuilt, hydrogens go on at the pH, then the minimisation. -->
  <div class=prepgroup id=prepgroup>
    <div class=grouphead>Receptor preparation
      <span>in the order it is done</span></div>
    <div class="checks prepsteps">
      <label class=chk><input type=checkbox id=keep_hetatms>
        <span>Keep heteroatoms</span></label>
      <!-- On by default, which is what every run before this setting did. -->
      <label class=chk><input type=checkbox id=rebuild_sidechains checked>
        <span>Rebuild missing side chains</span></label>
      <label class=chk><input type=checkbox id=fill_loops>
        <span>Rebuild missing loops</span></label>
      <!-- reduce by default, as ADFR's own preparation and ADapt protonate;
           PROPKA when the pH matters, and the pH box is for it alone. A run
           made before the choice existed used PROPKA, and reusing its
           settings sets it back. -->
      <label class=chk for=protonation><span>Protonate with</span>
        <select id=protonation>
          <option value=reduce selected>reduce</option>
          <option value=propka>PROPKA at a pH</option>
        </select></label>
      <label class=chk for=ph id=phlabel><span>pH</span>
        <input type=number id=ph value="7.4" step="0.1" min="0" max="14"></label>
      <label class=chk><input type=checkbox id=minimise checked>
        <span>Minimise the receptor after protonation</span></label>
    </div>
    <div class=why id=prepnote></div>
  </div>
  <!-- A target made in AGFR or AGFRGUI replaces the receptor preparation
       and the site: the run takes its receptor, maps and box from it. -->
  <div class=row style="margin-top:12px">
    <div class=fieldwrap style="flex:3 1 280px">
      <label for=targetfile>or a target file made in AGFR or AGFRGUI (.trg)</label>
      <input type=text id=targetfile
             placeholder="/path/to/target.trg  (optional: its receptor, maps and box are used)"></div>
    <div style="flex:0 0 auto"><button class=secondary id=trgload>Load</button></div>
    <div style="flex:0 0 auto"><button class=secondary id=trgbrowse>Browse&hellip;</button></div>
  </div>
  <div id=targetnote class=why style="margin-top:9px"></div>
  <!-- ADapt, a separate tool, shown only when the server can find it. It
       shows what a structure file holds and prepares the receptor the way
       AGFR does. Once a structure is loaded it opens here on it by itself,
       in a frame through this server's /adapt/, and the line above the
       frame always says which receptor will be docked: this page's own
       preparation until anything is changed in ADapt, then the receptor as
       set there, which ADapt writes when the run is submitted. -->
  <div id=adaptsec class="prepgroup adaptsec hide">
    <div class=grouphead>ADapt <span>what the structure holds, and what stays in the receptor</span></div>
    <div class=why id=adaptsays style="margin-top:6px">ADapt, a separate tool by
      Michel Sanner, opens here on the structure once one is loaded above. It
      shows each assembly, the hetero groups, metals, alternate locations and
      modified residues, and lets you choose what stays in the receptor.</div>
    <div class=bar id=adaptacts style="margin-top:8px">
      <button class="primary hide" id=adaptnow type=button>Use ADapt&rsquo;s receptor now</button>
      <button class="secondary hide" id=adaptdrop type=button>Use this page&rsquo;s preparation instead</button>
    </div>
    <!-- The receptor ADapt wrote, once it is read back: a status card with
         the one thing that matters as its headline, the key facts as chips,
         the long details folded away, and only the action the state calls
         for. app.js (adaptCard) fills it. Above the frame, where it is seen
         without scrolling past ADapt. -->
    <div id=adaptuse class="aucard hide" role=region aria-labelledby=autitle>
      <div class=auhead>
        <span class=auicon aria-hidden=true><svg viewBox="0 0 24 24"><path
          d="M7 7.5l10 3M7 7.5l3 10M17 10.5l-7 7" fill=none stroke=currentColor
          stroke-width=1.8 stroke-linecap=round /><circle cx=7 cy=7.5 r=3 /><circle
          cx=17 cy=10.5 r=3 /><circle cx=10 cy=17.5 r=3 /></svg></span>
        <div class=aumain>
          <div class=autitle id=autitle>Docking will use ADapt&rsquo;s receptor</div>
          <div class=ausub id=ausub></div>
        </div>
        <span class=austate id=austate role=status aria-live=polite></span>
      </div>
      <div class=auchips id=auchips></div>
      <div class=aunote id=aunote></div>
      <div class=auacts id=auacts>
        <button class="primary hide" id=adaptupdate type=button>Use ADapt&rsquo;s current receptor</button>
        <button class="primary hide" id=adaptread type=button>Read the newer save</button>
        <button class=secondary id=adaptback type=button>Use this page&rsquo;s preparation instead</button>
      </div>
      <details class=aumore id=aumore><summary>Details</summary>
        <div class=why id=adaptwhat></div></details>
    </div>
    <!-- The frame: tall enough for ADapt's tree beside its 3D view. It
         may not navigate this page (no allow-top-navigation). -->
    <div id=adaptpanel class=hide>
      <div class=adaptbar>
        <label for=adaptframe class=adaptlab>ADapt on <span id=adaptwhat2>this structure</span></label>
        <span class=why id=adaptwait></span>
        <a class="btnlink secondary hide" id=adaptfull target=_blank rel=noopener
           href="#">Open full screen</a>
        <button class=ghost id=adaptfold type=button aria-expanded=true
                aria-controls=adaptframewrap>Hide</button>
      </div>
      <!-- While ADapt starts, the space the frame will fill holds the
           loading stage (app.js, Wait): a quiet outline of ADapt's tree and
           view, and over it a card whose steps tick as they happen. -->
      <div id=adaptframewrap class=hide>
        <iframe id=adaptframe title="ADapt" loading=eager
          sandbox="allow-scripts allow-same-origin allow-forms allow-downloads allow-popups allow-popups-to-escape-sandbox"
          referrerpolicy=same-origin></iframe>
      </div>
    </div>
    <!-- On a phone the frame is too narrow for ADapt's tree and view, so a
         button opens it in a tab of its own on the same session. -->
    <div id=adaptrow class="row hide" style="margin-top:10px;align-items:center">
      <div style="flex:0 0 auto"><button class=secondary id=adaptopen type=button>Open in ADapt</button></div>
      <div class=why id=adaptnote style="flex:1 1 220px">Opens ADapt in a tab of its
        own, on the structure above.</div>
      <div id=adaptrowwait style="flex:1 1 100%"></div>
    </div>
  </div>
</div>

<div class=card id=step-site data-step=site>
  <h2>Site</h2>
  <p class=hint>Where the docking box goes: around residues you pick, a ligand
    or peptide in the structure, or a pocket AutoSite finds, padded and held
    between the limits.</p>
  <!-- Three ways to place the box. Residues is what every run before the
       choice did. A ligand or peptide is listed from the structure loaded
       above. An AutoSite pocket is found on the compute node once the
       receptor is prepared, so it has no box to draw until then. -->
  <div class=row>
    <div style="flex:0 0 250px"><label for=site_mode>Build the box around</label>
      <select id=site_mode>
        <option value=residues selected>residues you pick</option>
        <option value=ligand>a ligand or peptide</option>
        <option value=pocket>an AutoSite pocket</option>
      </select></div>
    <div id=ligandfield class=hide style="flex:2 1 280px">
      <label for=site_ligand>Ligand or peptide</label>
      <select id=site_ligand><option value="">load a receptor first</option></select></div>
    <div id=pocketfield class=hide style="flex:2 1 280px">
      <label for=site_pocket>Pocket</label>
      <select id=site_pocket>
        <option value=best selected>the best scoring pocket</option>
        <option value=all>every pocket it finds, in one box</option>
      </select></div>
  </div>
  <div class=row style="margin-top:10px">
    <div class=fieldwrap id=hotfield style="flex:3 1 300px">
      <label for=hotspots>Residues (CHAIN:NUMBER, or bare numbers)</label>
      <input type=text id=hotspots placeholder="A:41, A:45, A:78"></div>
    <div><label for=pad>Padding (&#8491;)</label>
      <input type=number id=pad value="8" step="0.5" min="0" max="20"></div>
    <div><label for=bmin>Min box</label>
      <input type=number id=bmin value="20" step="1" min="10" max="40"></div>
    <div><label for=bmax>Max box</label>
      <input type=number id=bmax value="30" step="1" min="10" max="60"></div>
    <div style="flex:0 0 auto"><button class=primary id=sitecheck>Check site</button></div>
  </div>
  <div id=sitenote class=why style="margin-top:9px"></div>
  <!-- A box around an interface groove often covers the chain on the far
       side. Docking without it places the peptide against a gap; docking with
       it may dock against the partner the peptide should displace. The user
       decides. -->
  <div id=otherchains class="why hide" style="margin-top:9px">
    <label class=why style="display:flex;gap:8px;align-items:flex-start">
      <input type=checkbox id=box_chains style="width:auto;margin-top:2px">
      <span id=otherchainstext></span>
    </label>
  </div>
  <!-- How the grid maps over the box are built: three settings of one
       program, agfr, in one panel like the receptor's preparation. -->
  <div class=prepgroup id=mapsgroup>
    <div class=grouphead>Grid maps <span>what the box is filled with</span></div>
    <div class="checks prepsteps">
      <label class=chk for=agfr><span>Built by</span>
        <select id=agfr>
          <option value=auto selected>AGFR, newest installed</option>
          <option value="1.1">AGFR 1.1, from the ADCP 1.1 suite</option>
          <option value="1.0">AGFR 1.0, from ADFRsuite</option>
        </select></label>
      <label class=chk for=autosite><span>AutoSite</span>
        <select id=autosite>
          <option value="1.1" selected>1.1</option>
          <option value="1.0">1.0</option>
        </select></label>
      <label class=chk><input type=checkbox id=receptor_gradient>
        <span>Add the receptor gradient</span></label>
    </div>
    <div class=why id=mapsnote></div>
  </div>
  <div id=chips style="margin-top:9px"></div>
  <div class=sitegrid id=sitepick>
    <div>
      <div class=pickhead>
        <span class=why id=pickhint>Click a residue in the structure to add
          it to the site. Click it again to drop it.</span>
        <span class=pickticks>
          <label><input type=checkbox id=sitesurf>Show the surface</label>
          <label><input type=checkbox id=showbox checked>Show the
            docking box</label></span>
      </div>
      <div id=viewD></div>
      <div class=why id=surfnote></div>
    </div>
    <div>
      <div class=pickhead>
        <span class=why><span id=seqhint>Or click it in the sequence.</span>
          <span id=seqchain></span></span>
      </div>
      <div id=seqstrip class=seq2></div>
    </div>
  </div>
</div>

<div class=card id=step-peptides data-step=peptides>
  <h2>Peptides</h2>
  <p class=hint>One-letter codes; pasted FASTA, spaces and case are fine.
    For non-standard and D amino acids (ADCP 1.1), use a peptide's &#8942; menu.</p>
  <div id=peps></div>
  <div class=bar>
    <button class=primary id=addpep>+ Add peptide</button>
    <button class=secondary id=pastepep>Paste a list&hellip;</button>
    <!-- The server has always accepted .txt, .csv, .fasta, .fa and .seq
         peptide lists; there was simply no control to send one, so a list
         could only be pasted. The file is read in the browser, the same way
         the receptor is, because the page is served through a tunnel and the
         cluster cannot see the machine the browser runs on. -->
    <button class=secondary id=pepfilebtn>Upload a list&hellip;</button>
    <input type=file id=pepfile style="display:none"
           accept=".txt,.csv,.fasta,.fa,.seq,text/plain">
    <span style="flex:1"></span>
    <!-- Its own element. #pepcount is rewritten by the validation debounce
         160 ms after any change, so a message put there disappears. -->
    <span class=why id=pepmsg></span>
    <span class=why id=pepcount></span>
  </div>
  <div id=pastebox class="hide" style="margin-top:11px">
    <label for=pastetext>One per line: <code>SEQUENCE</code>, <code>name: SEQUENCE</code>, FASTA, or columns with a ring, such as <code>p53,ETFSDLWKLLPEN,head to tail</code></label>
    <textarea id=pastetext placeholder="RWGPGAWGSGSGHY&#10;KDGGVDGPGDS"></textarea>
    <div style="margin-top:8px"><button class=secondary id=pasteadd>Add these</button></div>
  </div>
  <!-- How every peptide in the run is built and started: as a chain or a
       ring, and how many replicas begin as a helix. Peptide settings rather
       than search settings, so they sit with the peptides. The ring here is
       the run's: a row left "as the run" takes it, and a row may choose its
       own. -->
  <div class=row style="margin-top:12px">
    <div style="flex:0 0 250px"><label for=cyclic>The run docks each peptide as</label>
      <select id=cyclic>
        <option value=none selected>a linear chain</option>
        <option value=backbone>a ring, head to tail</option>
        <option value=disulfide>a ring through a disulfide</option>
      </select></div>
    <div style="flex:0 0 210px"><label for=helical>Replicas that start helical (%)</label>
      <input type=number id=helical value="25" min="0" max="100" step="5"></div>
    <div class=why id=pepoptsnote style="flex:1 1 240px;align-self:center;
         text-align:left"></div>
  </div>
</div>

<!-- One panel rather than a tooltip per field: a title= attribute appears
     only on hover, never on a touch screen, and cannot hold a paragraph.
     Clicking any "?" opens this. -->
<!-- Opened by the button on each run's row, where the user is when they want
     it. -->
<div id=issuebox class="card hide"
     style="position:fixed;left:50%;top:12%;transform:translateX(-50%);
            width:min(560px,92vw);z-index:60;
            box-shadow:0 10px 40px rgba(0,0,0,.3)">
  <button class=ghost id=issueclose
          style="float:right;margin:-4px -4px 0 8px">close</button>
  <h3 style="margin:0 0 4px">Report an issue</h3>
  <div class=why style="margin-bottom:9px">This sends the logs (not the
    poses) with the settings and what Slurm said.
    <div id=issuerun style="font-family:ui-monospace,monospace;font-size:12px;
      margin-top:4px"></div></div>
  <textarea id=issuenote rows=4
            placeholder="what you did, and what happened instead"></textarea>
  <input id=issuemail style="margin-top:8px"
         placeholder="your email, so we can reply (optional)">
  <div style="margin-top:10px"><button class=primary id=issuesend>Send the
    logs</button></div>
  <div class=why id=issuemsg style="margin-top:9px"></div>
</div>

<div id=helpbox class="card hide"
     style="position:fixed;right:18px;bottom:18px;max-width:420px;z-index:50;
            box-shadow:0 8px 30px rgba(0,0,0,.28)">
  <button class=ghost id=helpclose
          style="float:right;margin:-4px -4px 0 8px">close</button>
  <h3 id=helptitle style="margin:0 0 8px"></h3>
  <div id=helpbody class=why style="white-space:pre-wrap;line-height:1.55"></div>
  <a id=helpdocs class=hide href="#" target=_blank rel=noopener
     style="display:inline-block;margin-top:10px;font-size:13px">Read more in the Docs</a>
</div>

<div class=card id=step-docking data-step=docking>
  <h2>Docking and rescoring</h2>
  <!-- Four parts, each under a small heading: the search, the engine, the
       rescoring and the resources. The fields of a part sit on one grid (CSS
       at .dpart), and a caption sits under the field it describes, in that
       field's column. -->
  <div class="dpart cols3">
    <div class=grouphead>Search</div>
    <!-- Three equal fields, with nothing taller than the rest, so the labels
         line up. The step guideline is its own line below the row. The
         preparation pH was here and is now with the rest of the receptor's
         preparation, in the Receptor card. -->
    <div class=row>
      <div><label for=poses>Poses per peptide</label>
        <input type=number id=poses value="10" min="1" max="100"></div>
      <div><label for=replicas>Search replicas</label>
        __REPLICA_OPTIONS__</div>
      <div><label for=steps>Steps per replica</label>
        __STEP_OPTIONS__</div>
    </div>
    <div class=why id=stepsguide></div>
    <!-- Off by default: every replica's trajectory is hundreds of files and
         tens of megabytes per peptide. On, they are kept as one archive. -->
    <div class=checks>
      <label class=chk><input type=checkbox id=keep_trajectories>
        <span>Keep every replica's trajectory</span></label>
    </div>
  </div>
  <!-- Which ADCP. 1.1 adds non-standard and D amino acids (N<BHD>, &D in a
       sequence) and an OpenMM minimization of the top poses. The version
       found is the engine menu's caption, under it in its column. -->
  <div class="dpart cols3 dengine">
    <div class=grouphead>Engine</div>
    <div class=row>
      <div><label for=adcp_engine>Docking engine</label>
        <select id=adcp_engine>
          <option value=auto>ADCP, newest installed</option>
          <option value="1.1">ADCP 1.1</option>
          <option value="1.0">ADCP 1.0</option>
        </select></div>
      <div class="why dcap" id=enginenote></div>
      <div><label for=omm_minimize>Poses to minimize with OpenMM</label>
        <input type=number id=omm_minimize value="0" min="0" max="10"></div>
    </div>
  </div>
  <!-- The rescoring is decided after the docking, so its switch stays in
       this card, at the head of its part. The three receptor preparation
       switches that shared its row are in the Receptor card. -->
  <div class="dpart cols4">
    <div class=grouphead>Rescoring</div>
    <div class=checks>
      <!-- On by default. The docking score orders poses of the same
           peptide well and is a poor way to compare different peptides,
           which is what a screen is for, so a run without the rescoring
           produces a table that cannot be ranked on. It is switched off
           and disabled below if AmberTools is not installed. -->
      <label class=chk><input type=checkbox id=mmgbsa checked>
        <span>Rescore with MM-GBSA</span></label>
    </div>
    <!-- Unticked, these stay in place, greyed and disabled (app.js,
         mmOpts), so nothing below them moves. -->
    <div id=mmopts class=row>
      <div><label for=mmposes>Poses to rescore</label>
        <input type=number id=mmposes value="5" min="1" max="10"></div>
      <div><label for=mmmodel title="The implicit solvent model the MM-GBSA energy is computed in. GB-OBC(II) is the default and what every earlier run used; GBn2 is a newer model with its own atomic radii. Both run in AmberTools and both give per residue energies. The two give different numbers, so compare peptides only within one model.">Solvation model</label>
        <select id=mmmodel>
          <option value="obc2" selected>GB-OBC(II), the default</option>
          <option value="gbn2">GBn2</option>
        </select></div>
      <div><label for=mmengine title="Which program computes the MM-GBSA energy. AmberTools is the reference and every earlier number came from it. OpenMM scores the same pose, on the same geometry sander minimised, with its own generalised Born code: it lands within about 1.3 kcal/mol of AmberTools and is a check on the solver rather than a different estimate. Schrodinger Prime uses OPLS with VSGB 2.0, a different model entirely, on its own scale. Each writes its own column and none is compared with another by value.">MM-GBSA engine</label>
        <select id=mmengine>
          <option value="amber" selected>AmberTools, the reference</option>
          <option value="openmm">OpenMM, a second opinion on the same pose</option>
          <option value="schrodinger">Schrodinger Prime, OPLS</option>
        </select></div>
      <!-- The force field, the solvent model and the radii the chosen engine
           computes with, under the solvation and engine menus. -->
      <div class="why dcap" id=mmcap><span class=note id=mmenginenote></span></div>
      <div><label for=topn>Report the mean of the top</label>
        <input type=number id=topn value="5" min="1" max="50"></div>
    </div>
  </div>
  <div class="dpart cols4">
    <div class=grouphead>Resources</div>
    <div class=row>
      <div><label for=cpus>Cores per task</label>
        <div class=stepper>
          <button class=st data-step="cpus:-4" type=button>&minus;</button>
          <input type=number id=cpus value="48" min="1" max="256" step="4">
          <button class=st data-step="cpus:4" type=button>+</button></div></div>
      <div><label for=walltime>Walltime</label>
        <div class=stepper>
          <button class=st data-step="walltime:-30" type=button>&minus;</button>
          <input type=text id=walltime value="02:00:00">
          <button class=st data-step="walltime:30" type=button>+</button></div></div>
      <div><label for=throttle>Max tasks at once</label>
        <div class=stepper>
          <button class=st data-step="throttle:-1" type=button>&minus;</button>
          <input type=number id=throttle value="0" min="0" max="200">
          <button class=st data-step="throttle:1" type=button>+</button></div></div>
      <div><button class=primary id=autofill title="Work out how long this will take and ask
the queue for twice that and half an hour more, and at least two hours">Estimate &amp; fill</button></div>
    </div>
    <!-- Under the row the button fills: what the estimate assumed, and what
         the jobs will book. -->
    <div id=estnote class=why></div>
  </div>
</div>

<!-- The last step: what the run is called, where it is written, the
     cluster's own options, and the button. A card of its own, so the guide's
     fifth step has a card to light. -->
<div class=card id=step-submit data-step=submit>
  <h2>Submit</h2>
  <p class=hint>Name the run, choose where it is written, and send it to the
    queue.</p>
  <!-- Where the run goes: the output folder and the run name are the last
       decisions, so they sit beside the submit button. There is no partition
       field: a partition is a property of the cluster, and a default name
       fails on any cluster that does not have it. -->
  <div class=row>
    <div class=fieldwrap style="flex:1 1 220px"><label for=label>Name this run *</label>
      <input type=text id=label placeholder="e.g. p53_screen"></div>
    <div class=fieldwrap style="flex:2 1 320px"><label for=outdir>Output folder *</label>
      <input type=text id=outdir placeholder="where the run directory goes">
      </div>
    <div style="flex:0 0 auto"><label>&nbsp;</label>
      <button class=secondary id=outbrowse>Browse&hellip;</button></div>
  </div>
  <div id=outmsg class=why style="margin-top:6px"></div>
  <!-- This cluster's own rules. Many clusters refuse a job without an
       --account, and some need a QOS, a constraint or a particular partition.
       They are properties of the site, typed once, and carried into every job
       script the run writes. Leave them empty where the cluster needs none. -->
  <div class=row style="margin-top:10px">
    <div class=fieldwrap style="flex:1 1 200px"><label for=account>Account
      <button class=q type=button data-help=account>?</button></label>
      <input type=text id=account placeholder="only if your cluster needs one"></div>
    <div class=fieldwrap style="flex:2 1 320px"><label for=sbatch_extra>Other
      SLURM options <button class=q type=button data-help=sbatch_extra>?</button></label>
      <input type=text id=sbatch_extra
             placeholder="e.g. --qos=normal --partition=short"></div>
  </div>

  <div class=centre><button class=primary id=submit>Submit to the queue</button></div>
  <div id=submsg class=why style="margin-top:11px;text-align:center"></div>
</div>
</div><!-- .stepsflow -->

<div class=card>
  <h2>Runs</h2>
  <!-- Two views of the same runs folder: every run in one table, or the
       folder itself, walked folder by folder down to the files. Tabs in the
       WAI-ARIA pattern: the arrow keys move between them. -->
  <div class=tabs role=tablist aria-label="How the runs are shown">
    <button type=button role=tab id=tab_runs aria-selected=true
            aria-controls=runspanel data-tab=runs>All runs</button>
    <button type=button role=tab id=tab_folders aria-selected=false
            aria-controls=folderspanel data-tab=folders tabindex=-1>Folders</button>
  </div>
  <div class=bar style="margin-top:0;margin-bottom:10px">
    <button class=ghost id=refresh>Refresh</button>
    <!-- Runs made elsewhere, or by someone else, are reached by pointing this
         at their folder. Without it the only runs visible were the ones this
         session happened to write. -->
    <button class=secondary id=runsbrowse>Browse&hellip;</button>
    <span class=why id=runsroot></span>
  </div>
  <div id=runspanel role=tabpanel aria-labelledby=tab_runs>
    <div id=runs></div>
  </div>
  <!-- The runs folder as folders and files. Only what is inside the runs
       folder is listed: the server holds every path to it. A folder that is
       a run opens as one; a structure, a log or a table opens in the page;
       anything else downloads. -->
  <div id=folderspanel role=tabpanel aria-labelledby=tab_folders class=hide>
    <nav class=fcrumbs id=fcrumbs aria-label="Where in the runs folder"></nav>
    <div id=fhere class=fhere></div>
    <div id=flist></div>
  </div>
  <div id=rundetail style="margin-top:14px"></div>
</div>

<div class="card hide" id=viewcard>
  <h2>Pose</h2>
  <div class=vbar id=rvbar>
    <span class=vgrp>
    <span class=lb>colour</span>
    <button class="pill on" data-rcol=element>element</button>
    <button class=pill data-rcol=chain>chain</button>
    <button class=pill data-rcol=plain>plain</button>
    <button class=pill data-rcol=custom>custom</button>
    <label class=cpick title="receptor colour">rec
      <input type=color id=rccol value="#4A6FA5"></label>
    <label class=cpick title="peptide colour">pep
      <input type=color id=rlcol value="#D4763A"></label>
    <label class=cpick title="colour of everything outside the close-up">bg
      <input type=color id=rbcol value="#DFE3E6"></label>
    </span>
    <span class=vsep></span>
    <span class=vgrp>
    <span class=lb>draw</span>
    <button class="pill on" data-rdraw=cartoon>cartoon</button>
    <button class=pill data-rdraw=trace>trace</button>
    <button class=pill data-rdraw=sticks>sticks</button>
    </span>
    <span class=vsep></span>
    <span class=vgrp>
    <span class=lb>surface</span>
    <button class="pill on" id=rsurf title="The receptor's molecular surface, computed by MSMS, in the receptor's colour. See-through to start with, so the peptide stays in view. Cycles mesh, off, solid, see-through.">rec see-through</button>
    <button class=pill id=rsurfl title="The peptide's molecular surface, computed by MSMS, in the peptide's colour. Cycles solid, see-through, mesh, off.">pep off</button>
    </span>
    <span class=vsep></span>
    <span class=vgrp>
    <span class=lb>contacts</span>
    <button class="pill on" id=rcontact>on</button>
    <span class=lb>distances</span>
    <button class="pill on" id=rdist title="The measured heavy-atom distance on every contact drawn, in angstroms. A dashed line says two residues touch; the number says how closely, which is the difference between a hydrogen bond and two atoms in the same neighbourhood.">on</button>
    <span class=lb>figure</span>
    <button class=pill id=rpub title="Active-site view for a paper: pale ribbon, contacting residues and the peptide as sticks, bonds as black dashes, labels as plain text.">off</button>
    <select id=rlabel style="width:auto;padding:4px 9px;font-size:11.5px">
      <option value=auto>labels: auto</option>
      <option value=boxed>labels: boxed</option>
      <option value=plain>labels: plain</option>
      <option value=none>labels: off</option>
    </select>
    <button class=pill id=rdrag title="Drag the residue labels out of the way. A line keeps each one joined to its atom, and it is drawn into a saved image too.">move labels</button>
    <button class=pill id=rdragoff title="Put every label back on its atom">labels back</button>
    </span>
    <span class=vsep></span>
    <span class=vgrp>
    <span class=lb>view</span>
    <button class=pill id=rreset>reset</button>
    <button class=pill id=rspin>spin</button>
    <button class=pill id=rshot>save image</button>
    </span>
    <span class=vsep></span>
    <span class=vgrp>
    <span class=lb>engine</span>
    <select id=rengine style="width:auto;padding:4px 9px;font-size:11.5px"
            title="Which viewer draws the pose. Both are bundled; Mol* is the
larger of the two and takes a moment the first time.">
      <option value=3dmol>3Dmol</option>
      <option value=molstar>Mol*</option>
    </select>
    </span>
    <span class=vsep></span>
    <span class=vgrp>
    <span class=lb>pose</span>
    <select id=rposepick style="width:auto;max-width:210px;padding:4px 9px;
      font-size:11.5px"></select>
    </span>
</div>
  <div id=viewR></div>
  <div id=viewM class=hide></div>
  <div class=why id=viewmsg style="margin-top:7px"></div>
  <div class=why id=surfmsg style="margin-top:2px;font-size:11.5px"></div>
  <div id=contactbox style="margin-top:10px"></div>
  <div id=diagwrap style="margin-top:14px"></div>
  <div id=residuebox style="margin-top:16px"></div>
</div>
<dialog id=citedlg>
  <div class=dlg-h><b>How to cite</b>
    <span class=why id=citescope></span></div>
  <div class=listing><pre id=citetext
    style="white-space:pre-wrap;margin:0;font-size:12px"></pre></div>
  <div class=dlg-f>
    <span class=why id=citenote>Download for your reference manager:</span>
    <span style="display:flex;gap:8px;align-items:center">
      <button class=secondary data-cite=txt>.txt</button>
      <button class=secondary data-cite=bib>BibTeX</button>
      <button class=secondary data-cite=ris>RIS / EndNote</button>
      <button class=secondary id=citecopy>Copy</button>
      <button class=ghost id=citeclose>Close</button>
    </span>
  </div>
</dialog>
<!-- Feedback as a dialog rather than a jump to mailto: a mail link opened
     from inside a notebook's frame is blocked without a word. The addresses
     and the details to paste are always shown, so a report can be sent
     whatever the browser allows. -->
<div class="pepmenu hide" id=pepmenu role=menu>
  <button id=pm_mods role=menuitem><span class=mi><svg width="16" height="16" viewBox="0 0 24 24" aria-hidden="true"><path fill="currentColor" d="M11 5h2v6h6v2h-6v6h-2v-6H5v-2h6z"/></svg></span>Non-standard and D amino acids</button>
  <button id=pm_dup role=menuitem><span class=mi><svg width="16" height="16" viewBox="0 0 24 24" aria-hidden="true"><path fill="currentColor" d="M8 3h11v13h-2V5H8V3zM5 7h10v14H5V7zm2 2v10h6V9H7z"/></svg></span>Duplicate</button>
  <button id=pm_rm role=menuitem><span class=mi><svg width="16" height="16" viewBox="0 0 24 24" aria-hidden="true"><path fill="currentColor" d="M9 3h6l1 2h4v2H4V5h4l1-2zm-3 6h12l-1 12H7L6 9zm4 2v8h2v-8h-2zm4 0v8h2v-8h-2z"/></svg></span>Delete</button>
</div>
<dialog id=moddlg>
  <div class=dlg-h><b>Non-standard and D amino acids</b></div>
  <div class=modbody>
    <div class=modnote>Click a residue, then choose its D form or a non-standard
      residue from ADCP 1.1's library. The library comes from SwissSidechain:
      phosphorylated, oxidised and halogenated residues among many synthetic
      ones. Saving changes the sequence here; it is not typed.</div>
    <div id=modseq></div>
    <div id=modrows></div>
    <div id=modviews></div>
  </div>
  <div class=dlg-f>
    <span class=why id=modmsg></span>
    <span style="display:flex;gap:8px;align-items:center">
      <button class=ghost id=modcancel>Cancel</button>
      <button class=primary id=modsave disabled>Save</button>
    </span>
  </div>
</dialog>
<!-- A file from the runs folder, opened from the Folders view: a structure
     in the 3D viewer, a log or a table as text. -->
<dialog id=filedlg aria-labelledby=filetitle>
  <div class=dlg-h><b id=filetitle>File</b>
    <span class=crumbs id=filepath></span></div>
  <div id=file3d class=hide></div>
  <pre id=filetext class=hide></pre>
  <div class=dlg-f>
    <span class=why id=filenote></span>
    <span style="display:flex;gap:8px;align-items:center">
      <a class="btnlink secondary" id=filedl href="#">Download</a>
      <button class=ghost id=fileclose>Close</button>
    </span>
  </div>
</dialog>
<dialog id=folderdlg>
  <div class=dlg-h><b>Run folder</b>
    <span class=why id=folderhint></span></div>
  <div class=listing><pre id=foldertext
    style="white-space:pre-wrap;word-break:break-all;margin:0;font-size:12.5px"></pre></div>
  <div class=dlg-f>
    <span class=why id=foldernote></span>
    <span style="display:flex;gap:8px;align-items:center">
      <button class=secondary id=foldercopy>Copy path</button>
      <button class=ghost id=folderclose>Close</button>
    </span>
  </div>
</dialog>
<dialog id=fbdlg>
  <div class=dlg-h><b>Send feedback</b>
    <span class=why>Write to aadhil.haq@tamu.edu, copying sfernando@tamu.edu,
      or open an issue on GitHub. Paste the details below into the message:
      they say which version and machine the report is about.</span></div>
  <div class=listing><textarea id=fbtext readonly rows=7 aria-label="Details
    to include in a report" style="width:100%;box-sizing:border-box;border:0;
    background:transparent;resize:vertical;font:12px ui-monospace,Menlo,Consolas,monospace"></textarea></div>
  <div class=dlg-f>
    <span class=why id=fbnote></span>
    <span style="display:flex;gap:8px;align-items:center;flex-wrap:wrap">
      <button class=secondary id=fbcopy>Copy details</button>
      <a class="btnlink secondary" id=fbmail target=_blank rel=noopener>Open email</a>
      <a class="btnlink secondary" id=fbissue target=_blank rel=noopener
         href="https://github.com/aadhilhaq/ADCP_docking/issues/new">Open a GitHub issue</a>
      <button class=ghost id=fbclose>Close</button>
    </span>
  </div>
</dialog>
<dialog id=picker>
  <div class=dlg-h><b id=pickertitle>Choose a file</b>
    <div class=crumbs id=crumbs></div></div>
  <!-- A scratch directory holds hundreds of folders and scrolling to one by
       eye is the slow part of this dialog. Typing filters what is listed,
       here rather than on the server, because the listing is already in the
       page and a round trip per keystroke would be slower than the scroll. -->
  <div class=pickfilter>
    <input type=search id=pickfind autocomplete=off
           placeholder="Type to filter this folder…"
           aria-label="Filter the list below">
    <span class=why id=pickfindnote></span>
  </div>
  <div class=listing id=listing></div>
  <div class=dlg-f>
    <span class=why id=pickernote></span>
    <span style="display:flex;gap:8px;align-items:center">
      <button class=ghost id=pickernew>+ New folder</button>
      <button class=primary id=pickeruse>Use this folder</button>
      <button class=ghost id=pickerclose>Close</button>
    </span>
  </div>
</dialog>

<!-- The footer carries the disclaimer; the credit sentence is in the header,
     so it is not repeated here. -->
<footer><span class=why>__DISCLAIMER__</span>
  <br><span class=why>How to cite this tool and everything it runs:
    <a href="#" id=citefoot>open the citation list</a>, or
    <code>adcp-dock cite</code>. Every download carries
    HOW_TO_CITE.txt, CITATION.bib and CITATION.ris for that run.</span>
  <br><span class=why>build __BUILD__. If this is
  not the commit you just pulled, the server is still running the old code:
  stop it, start it again, then reload with Ctrl-Shift-R.</span></footer>
</main>
__CREDITBAR__

<link rel="stylesheet" href="/static/molstar.css">
<script src="/static/3Dmol-min.js"></script>
<script src="/static/molstar.js"></script>
<!-- The limits, from the module that enforces them. app.js is a static file
     and cannot be substituted, so the page passes peptide.py's limits in this
     data island, the same way as the help text. -->
<script id=limits type="application/json">__LIMITS__</script>
<script id=helpdata type="application/json">__HELP_JSON__</script>
<script src="/static/drag.js"></script>
<script src="/static/app.js"></script>
"""

# The dark colour tokens, used twice in the page: for the system setting and
# for the switch in the header.
DARK_VARS = """
  color-scheme:dark;
  --bg:#0f1219; --card:#161b26; --panel:#1b2130; --line:#28303f; --fg:#e7ecf5;
  --dim:#94a0b8; --accent:#5b8cf0; --accent2:#1d2942;
  --hydro:#e0996a; --polar:#8fb3e8; --bad:#e2796a; --warn:#e0a06a;
  --good:#4cbf94; --soft:#141a24;
  --ribbon:#1b2231; --ribbonfg:#cdd8ec;
  /* The dark accent is light enough that white on it is about 3:1, under
     what small text needs, so a filled pill takes dark ink instead. */
  --pill-line:#3d4757; --on-accent:#0f1219;
"""
PAGE = PAGE.replace("__DARK_VARS__", DARK_VARS)
# The credit bar's style, as the documentation defines it.
PAGE = PAGE.replace("__CREDITBAR_CSS__", docskit.CREDITBAR_CSS.strip())



# The commit this process started with, read once at import. A running server
# holds page.py in memory, so this is the code actually being served however
# far the checkout moves afterwards.
_BUILD_AT_IMPORT = None


def _head():
    """The checkout's commit right now. None when it cannot be read."""
    import subprocess
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        r = subprocess.run(["git", "-C", here, "rev-parse", "--short", "HEAD"],
                           capture_output=True, text=True, timeout=5)
        return r.stdout.strip() or None
    except Exception:                                          # noqa: BLE001
        return None


def _build():
    """The build being SERVED, fixed at import whatever the checkout does."""
    global _BUILD_AT_IMPORT
    if _BUILD_AT_IMPORT is None:
        _BUILD_AT_IMPORT = _head() or "unknown"
    return _BUILD_AT_IMPORT


def _host():
    import socket
    try:
        return socket.gethostname().split(".")[0]
    except Exception:                                          # noqa: BLE001
        return "the machine it is running on"


def stale_banner():
    """Reported when the checkout has moved past what is being served.

    The server is the only thing that can see both the build it started with
    and the checkout on disk now.
    """
    served, head = _build(), _head()
    if not head or served == "unknown" or served == head:
        return ""
    return ('<div class=stale role=status>'
            '<b>This page is older than the version you have installed.</b>'
            '<br>It is served by a process that started at <code>'
            + _esc(served) + '</code>; the checkout is now at <code>'
            + _esc(head) + '</code>, so nothing pulled since is on this page.'
            '<br><b>Close the launcher and open it again</b>: that '
            'starts a fresh server with the code you pulled. Then reload '
            'this page with Ctrl-Shift-R.'
            '<br><span class=why2>Running it by hand instead? The process is '
            'on <code>' + _esc(_host()) + '</code>.</span></div>')


def _esc(v):
    return (str(v).replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def _number(n):
    """12500000 -> "12,500,000". The number, grouped, because that is what
    ADCP is given and what a paper quotes."""
    return "{:,}".format(int(n))


def choice_options(choices, selected, element_id):
    """A <select> of (value, hint) pairs, labelled with the NUMBER first."""
    out = ['<select id=%s>' % element_id]
    for value, hint in choices:
            # The number alone: a dropdown is a list of values to
            # choose between, and
        # the sentence explaining each belongs in the help panel,
        # where it can be a paragraph. The hint remains as the
        # option's title, for a hover.
        out.append('<option value="%d"%s title="%s">%s</option>'
                   % (value, " selected" if value == selected else "",
                      _esc(hint), _esc(_number(value))))
    out.append('</select>')
    return "".join(out)


def _local_pick(which):
    """Smaller starting values when the work runs on the local machine.

    The cluster defaults are sized for a compute node; on two cores they take
    the better part of an hour per peptide. Local execution starts at the
    smallest settings that are still a real search, and the dropdown still
    offers every other choice.
    """
    from . import slurm
    if not slurm.local_mode():
        return {}
    return {"selected": {"steps": 500_000, "replicas": 25}[which]}


def step_options(selected=2_500_000):
    """The search-effort choices, as a <select> of numbers.

    A free number box asks for a quantity few users have an intuition for, and
    a wrong answer is either a wasted allocation or a pose that moves when the
    run is repeated.
    """
    from . import helptext
    return choice_options(helptext.STEP_CHOICES, selected, "steps")


def replica_options(selected=50):
    """The same, for the number of independent searches."""
    from . import helptext
    return choice_options(helptext.REPLICA_CHOICES, selected, "replicas")


def favicon_href():
    """The launcher's logo, as a data URI for the tab icon."""
    import base64
    from . import brand
    raw = brand.LOGO_SVG.encode("utf-8")
    return "data:image/svg+xml;base64," + base64.b64encode(raw).decode("ascii")


def favicon_png_href():
    """The same mark as a PNG, for a browser that will not take the SVG."""
    import base64
    from . import brand
    return ("data:image/png;base64,"
            + base64.b64encode(brand.logo_png()).decode("ascii"))


def render(token):
    from . import credit
    import json
    from . import helptext
    from . import docs
    # "<" written as its escape, so a help text that names a path such as
    # docking/<peptide>/ can never read as a tag, or close the script it is
    # carried in.
    blob = json.dumps({k: {"title": t, "body": b, "docs": docs.HELP_LINKS.get(k, "")}
                       for k, (t, b) in helptext.HELP.items()}).replace(
                           "<", "\\u003c")
    # The bar the documentation draws along the bottom of every page, drawn
    # along the bottom of this one. Its How to cite opens the citation dialog
    # rather than a page, so app.js wires it by this id.
    bar = docskit.creditbar_html(docs.SITE, "/docs/index.html", cite_id="citebar")
    from . import peptide as _pep
    limits = json.dumps({"min_length": _pep.MIN_LENGTH,
                         "max_length": _pep.MAX_LENGTH})
    from . import DISCLAIMER
    from . import cite
    return (PAGE.replace("__DISCLAIMER__", DISCLAIMER)
                .replace("__MUST_CITE__", cite.must_cite_sentence(html=True))
                .replace("__LIMITS__", limits)
                .replace("__STEP_OPTIONS__", step_options(**_local_pick("steps")))
                .replace("__REPLICA_OPTIONS__",
                         replica_options(**_local_pick("replicas")))
                .replace("__HELP_JSON__", blob)
                .replace("__CREDITBAR__", bar)
                .replace("__BUILD__", _build())
                .replace("__STALE__", stale_banner())
                .replace("__FAVICON__", favicon_href())
                .replace("__FAVICON_PNG__", favicon_png_href()))
