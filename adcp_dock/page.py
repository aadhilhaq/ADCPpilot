"""The web page: one file, no build step, no network access at run time.

By Aadhil Haq

The peptide list is laid out as a row per sequence, with the entry number,
the type, and a sequence box that numbers its own residues.

Served by the same process that performs the work, so the page cannot be a
version ahead of the server it is talking to.
"""
import os

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
  document.documentElement.classList.add("embed");</script>
<style>
:root{
  --bg:#f7f8fb; --card:#fff; --panel:#eceef6; --line:#dfe3ee; --fg:#1b2333;
  --dim:#66718a; --accent:#1a56db; --accent2:#e8efff;
  --hydro:#b4632a; --polar:#1f3a6e; --bad:#c0392b; --warn:#c2571a;
  --good:#1a7f5a; --soft:#f2f4fa;
  /* the heading band across the top of a card, and the ink on it */
  --ribbon:#eef1f9; --ribbonfg:#2b3550;
}
@media (prefers-color-scheme:dark){:root{
  --bg:#0f1219; --card:#161b26; --panel:#1b2130; --line:#28303f; --fg:#e7ecf5;
  --dim:#94a0b8; --accent:#5b8cf0; --accent2:#1d2942;
  --hydro:#e0996a; --polar:#8fb3e8; --bad:#e2796a; --warn:#e0a06a;
  --good:#4cbf94; --soft:#141a24;
  --ribbon:#1b2231; --ribbonfg:#cdd8ec;
}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:14px/1.5 -apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
header{background:var(--card);border-bottom:1px solid var(--line);
  padding:18px 24px 14px;text-align:center}
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
  gap:10px;align-items:center;justify-content:space-between}

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
.missing,.missing:focus{border-color:var(--warn) !important;
  outline:2px solid var(--warn)}
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
.row>div{flex:1 1 150px;min-width:0}
/* Every control in a row has the same height. Otherwise a .ghost button (5px
   padding, 12.5px text) is shorter than the input beside it, and bottom
   alignment leaves it visibly low. A shared height lines up a row of mixed
   controls at both edges. */
.row input[type=text],.row input[type=number],.row select,
.row button,.row .fixedval{height:38px;box-sizing:border-box}
.row button{display:inline-flex;align-items:center;justify-content:center}
.row .stepper button{height:38px}
button{font:inherit;padding:8px 16px;border-radius:999px;cursor:pointer;
  border:1px solid var(--line);background:var(--card);color:var(--fg)}
button:hover{background:var(--soft)}
button.primary{background:var(--accent);border-color:var(--accent);color:#fff;
  font-weight:600;padding:10px 26px}
button.primary:hover{filter:brightness(1.08)}
button.ghost{font-size:12.5px;padding:5px 12px;color:var(--dim)}
button.ghost.danger{color:var(--bad)}
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
.ent .sqwrap{flex:1 1 auto}
.ent .tools{display:flex;gap:2px;align-self:center}
.ent label{margin-bottom:4px;height:14px;line-height:14px}
/* The three controls share one height, so their tops and bottoms line up.
   46px rather than 38: the sequence box carries its residue numbers above the
   letters, so it is taller than a plain input. */
.ent .fixedval,.ent .nm input{height:46px}
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
   the browser cannot save it. */
.btnlink{display:inline-block;font:inherit;font-size:12.5px;
  padding:5px 12px;border-radius:999px;border:1px solid var(--line);
  background:var(--card);color:var(--dim);text-decoration:none;
  cursor:pointer;vertical-align:middle;line-height:1.35}
.btnlink:hover{background:var(--soft);color:var(--fg);text-decoration:none}
.grp{position:relative;display:inline-block;margin-right:1.15em;
  letter-spacing:.19em}
.grp::after{content:attr(data-n);position:absolute;right:.19em;top:-1.35em;
  font-size:10px;letter-spacing:0;color:var(--dim);font-weight:500}
.seq .h{color:var(--hydro)} .seq .p{color:var(--polar)}
.seq .x{color:var(--bad);font-weight:700;background:rgba(192,57,43,.13);
  border-radius:2px}
.seqmsg{font-size:12.5px;margin-top:6px;color:var(--warn)}
.seqlen{font-size:11.5px;color:var(--dim);margin-top:5px}
.bar{display:flex;gap:12px;align-items:center;margin-top:14px;
  justify-content:center}
.centre{display:flex;justify-content:center;margin-top:18px}

/* ---- results and the viewer ----------------------------------------- */
.sitegrid{display:grid;grid-template-columns:minmax(0,1.15fr) minmax(0,1fr);
  gap:14px;margin-top:12px}
@media (max-width:820px){.sitegrid{grid-template-columns:minmax(0,1fr)}}
#viewD{width:100%;height:330px;border:1px solid var(--line);border-radius:10px;
  background:#fff;position:relative}
@media (prefers-color-scheme:dark){#viewD{background:#0b0e14}}
.seq2{font-family:ui-monospace,Menlo,monospace;font-size:12.5px;
  max-height:330px;overflow:auto;border:1px solid var(--line);
  border-radius:10px;padding:9px 11px;background:var(--card);line-height:1.85}
.seq2 .ln{display:flex;gap:9px;align-items:baseline;white-space:pre}
.seq2 .no{color:var(--dim);font-size:10.5px;min-width:3.8em;text-align:right;
  user-select:none}
.seq2 .r{cursor:pointer;border-radius:3px;padding:0 .5px}
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
/* Row actions have to fit on one line beside the run. */
.btnlink.sm,button.ghost.sm{font-size:11.5px;padding:3px 10px}
.sz{color:var(--dim);font-size:10.5px;margin-left:3px}
#issuebox textarea,#issuebox input{width:100%}
tr.pick{cursor:pointer} tr.pick:hover{background:var(--soft)}
tr.on{background:var(--accent2)}
.pill{display:inline-block;font-size:11px;padding:2px 9px;border-radius:999px;
  border:1px solid var(--line);color:var(--dim)}
.list{max-height:240px;overflow:auto;border:1px solid var(--line);
  border-radius:10px}
.list button{display:block;width:100%;text-align:left;border:0;border-radius:0;
  border-bottom:1px solid var(--line);background:none;padding:6px 11px;
  font-size:12.5px}
pre{background:var(--soft);padding:11px;border-radius:9px;overflow:auto;
  font-size:12px;max-height:240px}
.hide{display:none}
footer{color:var(--dim);font-size:12px;text-align:center;padding:18px 0 4px;
  border-top:1px solid var(--line);margin-top:22px}
/* the viewer's own bar: pills, grouped and labelled */
/* A group does not wrap inside itself, so a label is never separated from the
   buttons it names. */
.vgrp{display:inline-flex;gap:6px;align-items:center;flex:0 0 auto}
.vbar{display:flex;flex-wrap:wrap;gap:6px;align-items:center;
  justify-content:center;
  padding:9px 10px;background:var(--soft);border:1px solid var(--line);
  border-radius:10px;margin-bottom:10px}
.vbar .lb{font-size:10.5px;text-transform:uppercase;letter-spacing:.06em;
  color:var(--dim);font-weight:650;margin-right:1px}
.vbar .pill{font-size:11.5px;padding:4px 11px;border-radius:999px;
  border:1px solid var(--line);background:var(--card);color:var(--dim);
  cursor:pointer}
.vbar .pill.on{background:var(--accent);border-color:var(--accent);color:#fff;
  font-weight:600}
.vsep{width:1px;align-self:stretch;background:var(--line);margin:0 4px}
.cpick{display:inline-flex;align-items:center;gap:5px;font-size:11px;
  color:var(--dim)}
.cpick input{width:26px;height:22px;padding:0;border-radius:6px;
  border:1px solid var(--line);background:none}
#viewR,#viewM{width:100%;height:min(58vh,520px);border:1px solid var(--line);
  border-radius:10px;position:relative;background:#fff}
@media (prefers-color-scheme:dark){#viewR,#viewM{background:#0b0e14}}
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
</style>

__STALE__
<header>
  <div class=brandblock>
    <h1>GUI to ADCP Docking</h1>
    <!-- What belongs to whom, in the header: ADCP and ADFRsuite are the
         Scripps Research Institute's, and this is a GUI and workflow that
         drives them. A screenshot carries the header with it. -->
    <div class=sub>A graphical interface to
      <a href="https://ccsb.scripps.edu/adcp/" target="_blank"
         rel="noopener">ADCP</a>, a product of the
      <b>Scripps Research Institute</b>. <b>Powered by HOPE Lab.</b></div>
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
      <span class=sub id=toolnote></span>
    </div>
  </div>
</header>

<main>
<div class=card>
  <h2>Receptor</h2>
  <p class=hint>The protein the peptides will be docked into. It is protonated
    at the pH you choose before anything else happens.</p>
  <div class=row>
    <div style="flex:0 0 150px">
      <label for=reccode>PDB id</label>
      <input type=text id=reccode placeholder="4CMM" maxlength=4
             style="text-transform:uppercase"></div>
    <div style="flex:0 0 auto"><button id=recfetch>Fetch from RCSB</button></div>
    <div class=why style="flex:0 0 auto;padding-bottom:9px">or</div>
    <div class=fieldwrap style="flex:3 1 280px">
      <label for=recpath>a file already on this machine</label>
      <input type=text id=recpath
             placeholder="/path/to/receptor.pdb  (or .cif, .mmcif, .ent, .gz)"></div>
    <div style="flex:0 0 auto"><button id=recload>Load</button></div>
    <div style="flex:0 0 auto"><button class=ghost id=recbrowse>Browse&hellip;</button></div>
    <!-- From the user's own computer. "Browse" walks the cluster's
         filesystem; this reads a file on the user's machine and sends it
         through the connection that is already open. -->
    <div style="flex:0 0 auto">
      <button class=ghost id=recupload>From this computer&hellip;</button>
      <input type=file id=recfile style="display:none"
             accept=".pdb,.ent,.cif,.mmcif,.pdbqt,.txt"></div>
  </div>
  <div id=recinfo class=why style="margin-top:9px"></div>
  <div id=chainbox class="hide" style="margin-top:10px"></div>
  <!-- A target made in AGFR or AGFRGUI replaces the receptor preparation
       and the site: the run takes its receptor, maps and box from it. -->
  <div class=row style="margin-top:10px">
    <div class=fieldwrap style="flex:3 1 280px">
      <label for=targetfile>or a target file made in AGFR or AGFRGUI (.trg)</label>
      <input type=text id=targetfile
             placeholder="/path/to/target.trg  (optional: its receptor, maps and box are used)"></div>
  </div>
</div>

<div class=card>
  <h2>Site</h2>
  <p class=hint>The residues the docking box should cover. Along each axis the box
    spans them plus the padding, clamped between the minimum and maximum.</p>
  <div class=row>
    <div class=fieldwrap style="flex:3 1 300px">
      <label for=hotspots>Residues (CHAIN:NUMBER, or bare numbers)</label>
      <input type=text id=hotspots placeholder="A:41, A:45, A:78"></div>
    <div><label for=pad>Padding (&#8491;)</label>
      <input type=number id=pad value="8" step="0.5" min="0" max="20"></div>
    <div><label for=bmin>Min box</label>
      <input type=number id=bmin value="20" step="1" min="10" max="40"></div>
    <div><label for=bmax>Max box</label>
      <input type=number id=bmax value="30" step="1" min="10" max="60"></div>
    <div style="flex:0 0 auto"><button id=sitecheck>Check site</button></div>
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
  <div id=chips style="margin-top:9px"></div>
  <div class=sitegrid id=sitepick>
    <div>
      <div class=pickhead>
        <span class=why>Click a residue in the structure to add it to the
          site. Click it again to drop it.</span>
        <label><input type=checkbox id=showbox checked>Show the
          docking box</label>
      </div>
      <div id=viewD></div>
    </div>
    <div>
      <div class=pickhead>
        <span class=why>Or click it in the sequence.
          <span id=seqchain></span></span>
      </div>
      <div id=seqstrip class=seq2></div>
    </div>
  </div>
</div>

<div class=card>
  <h2>Peptides</h2>
  <p class=hint>One-letter codes; pasted FASTA, spaces and case are fine.
    For non-standard and D amino acids (ADCP 1.1), use a peptide's &#8942; menu.</p>
  <div id=peps></div>
  <div class=bar>
    <button id=addpep>+ Add peptide</button>
    <button class=ghost id=pastepep>Paste a list&hellip;</button>
    <span style="flex:1"></span>
    <span class=why id=pepcount></span>
  </div>
  <div id=pastebox class="hide" style="margin-top:11px">
    <label for=pastetext>One per line: <code>SEQUENCE</code>, <code>name: SEQUENCE</code>, or FASTA</label>
    <textarea id=pastetext placeholder="RWGPGAWGSGSGHY&#10;KDGGVDGPGDS"></textarea>
    <div style="margin-top:8px"><button class=ghost id=pasteadd>Add these</button></div>
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
</div>

<div class=card>
  <h2>Docking and rescoring</h2>
  <!-- Four equal fields, with nothing taller than the rest, so the labels
       line up. The step guideline is its own line below the row. -->
  <div class=row>
    <div><label for=poses>Poses per peptide</label>
      <input type=number id=poses value="10" min="1" max="10"></div>
    <div><label for=replicas>Search replicas</label>
      __REPLICA_OPTIONS__</div>
    <div><label for=steps>Steps per replica</label>
      __STEP_OPTIONS__</div>
    <div><label for=ph>Preparation pH</label>
      <input type=number id=ph value="7.4" step="0.1" min="0" max="14"></div>
  </div>
  <div class=why id=stepsguide style="margin-top:7px"></div>
  <!-- Which ADCP. 1.1 adds non-standard and D amino acids (N<BHD>, &D in a
       sequence) and an OpenMM minimization of the top poses. -->
  <div class=row style="margin-top:12px">
    <div><label for=adcp_engine>Docking engine</label>
      <select id=adcp_engine>
        <option value=auto>ADCP, newest installed</option>
        <option value="1.1">ADCP 1.1</option>
        <option value="1.0">ADCP 1.0</option>
      </select></div>
    <div><label for=omm_minimize>Poses to minimize with OpenMM</label>
      <input type=number id=omm_minimize value="0" min="0" max="10"></div>
    <div class=why id=enginenote style="flex:2 1 260px;align-self:center"></div>
  </div>
  <div class=row style="margin-top:12px">
    <!-- The checkboxes are their own row, evenly spaced: a group of four of
         the same kind, laid out as one. -->
    <div class=checks>
      <label class=chk><input type=checkbox id=minimise checked>
        <span>Minimise the receptor after protonation</span></label>
      <label class=chk><input type=checkbox id=keep_hetatms>
        <span>Keep heteroatoms</span></label>
      <label class=chk><input type=checkbox id=fill_loops>
        <span>Rebuild missing loops</span></label>
      <label class=chk><input type=checkbox id=mmgbsa>
        <span>Rescore with MM-GBSA</span></label>
    </div>
    <div id=mmopts class="row hide" style="flex:1 1 auto">
      <div><label for=mmposes>Poses to rescore</label>
        <input type=number id=mmposes value="5" min="1" max="10"></div>
      <div><label for=topn>Report the mean of the top</label>
        <input type=number id=topn value="5" min="1" max="50"></div>
    </div>
  </div>
  <div class=row style="margin-top:12px">
    <div style="flex:0 0 168px"><label for=cpus>Cores per task</label>
      <div class=stepper>
        <button class=st data-step="cpus:-4" type=button>&minus;</button>
        <input type=number id=cpus value="48" min="1" max="256" step="4">
        <button class=st data-step="cpus:4" type=button>+</button></div></div>
    <div style="flex:0 0 190px"><label for=walltime>Walltime</label>
      <div class=stepper>
        <button class=st data-step="walltime:-30" type=button>&minus;</button>
        <input type=text id=walltime value="02:00:00">
        <button class=st data-step="walltime:30" type=button>+</button></div></div>
    <div style="flex:0 0 150px"><label for=throttle>Max tasks at once</label>
      <div class=stepper>
        <button class=st data-step="throttle:-1" type=button>&minus;</button>
        <input type=number id=throttle value="0" min="0" max="200">
        <button class=st data-step="throttle:1" type=button>+</button></div></div>
    <div style="flex:0 0 auto"><label>&nbsp;</label>
      <button id=autofill title="Work out how long this will take and ask
the queue for twice that plus half an hour">Estimate &amp; fill</button></div>
  </div>

  <!-- Where the run goes: the output folder and the run name are the last
       decisions, so they sit beside the submit button. There is no partition
       field: a partition is a property of the cluster, and a default name
       fails on any cluster that does not have it. -->
  <div class=row style="margin-top:12px">
    <div class=fieldwrap style="flex:1 1 220px"><label for=label>Name this run *</label>
      <input type=text id=label placeholder="e.g. SIRPa_screen"></div>
    <div class=fieldwrap style="flex:2 1 320px"><label for=outdir>Output folder *</label>
      <input type=text id=outdir placeholder="where the run directory goes">
      </div>
    <div style="flex:0 0 auto"><label>&nbsp;</label>
      <button class=ghost id=outbrowse>Browse&hellip;</button></div>
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

  <div id=estnote class=why style="margin-top:9px"></div>
  <div class=centre><button class=primary id=submit>Submit to the queue</button></div>
  <div id=submsg class=why style="margin-top:11px;text-align:center"></div>
</div>

<div class=card>
  <h2>Runs</h2>
  <div class=bar style="margin-top:0;margin-bottom:10px">
    <button class=ghost id=refresh>Refresh</button>
    <span class=why id=runsroot></span>
  </div>
  <div id=runs></div>
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
    <button class=pill id=rsurf title="A surface over the receptor, in the receptor's colour. Cycles solid, see-through, mesh, off.">rec off</button>
    <button class=pill id=rsurfl title="A surface over the peptide, in the peptide's colour. Cycles solid, see-through, mesh, off.">pep off</button>
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
  <div id=contactbox style="margin-top:10px"></div>
  <div id=diagwrap style="margin-top:14px"></div>
</div>
<dialog id=citedlg>
  <div class=dlg-h><b>How to cite</b>
    <span class=why id=citescope></span></div>
  <div class=listing><pre id=citetext
    style="white-space:pre-wrap;margin:0;font-size:12px"></pre></div>
  <div class=dlg-f>
    <span class=why id=citenote>Download for your reference manager:</span>
    <span style="display:flex;gap:8px;align-items:center">
      <button class=ghost data-cite=txt>.txt</button>
      <button class=ghost data-cite=bib>BibTeX</button>
      <button class=ghost data-cite=ris>RIS / EndNote</button>
      <button class=ghost id=citecopy>Copy</button>
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
      <button id=modsave disabled>Save</button>
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
      <button class=ghost id=foldercopy>Copy path</button>
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
      <button class=ghost id=fbcopy>Copy details</button>
      <a class=btnlink id=fbmail target=_blank rel=noopener>Open email</a>
      <a class=btnlink id=fbissue target=_blank rel=noopener
         href="https://github.com/aadhilhaq/ADCP_docking/issues/new">Open a GitHub issue</a>
      <button class=ghost id=fbclose>Close</button>
    </span>
  </div>
</dialog>
<dialog id=picker>
  <div class=dlg-h><b id=pickertitle>Choose a file</b>
    <div class=crumbs id=crumbs></div></div>
  <div class=listing id=listing></div>
  <div class=dlg-f>
    <span class=why id=pickernote></span>
    <span style="display:flex;gap:8px;align-items:center">
      <button class=ghost id=pickernew>+ New folder</button>
      <button class=ghost id=pickeruse>Use this folder</button>
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
    blob = json.dumps({k: {"title": t, "body": b}
                       for k, (t, b) in helptext.HELP.items()})
    from . import peptide as _pep
    limits = json.dumps({"min_length": _pep.MIN_LENGTH,
                         "max_length": _pep.MAX_LENGTH})
    from . import DISCLAIMER
    return (PAGE.replace("__DISCLAIMER__", DISCLAIMER)
                .replace("__LIMITS__", limits)
                .replace("__STEP_OPTIONS__", step_options(**_local_pick("steps")))
                .replace("__REPLICA_OPTIONS__",
                         replica_options(**_local_pick("replicas")))
                .replace("__HELP_JSON__", blob)
                .replace("__BUILD__", _build())
                .replace("__STALE__", stale_banner())
                .replace("__FAVICON__", favicon_href())
                .replace("__FAVICON_PNG__", favicon_png_href()))
