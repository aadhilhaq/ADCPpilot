/* The page's behaviour. Plain JavaScript, no framework and no build step.
   Served by the same process that does the work, so it always matches the
   server it talks to. */
'use strict';
const T = new URLSearchParams(location.search).get('t') || '';
const $ = i => document.getElementById(i);
const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g,
  c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c]));
const fmt = (v, n) => (v == null ? '-' : Number(v).toFixed(n));

async function api(path, body){
  const url = path + (path.includes('?') ? '&' : '?') + 't=' + encodeURIComponent(T);
  const opt = body ? {method:'POST', headers:{'Content-Type':'application/json'},
                      body: JSON.stringify(body)} : {};
  const r = await fetch(url, opt);
  const ct = r.headers.get('Content-Type') || '';
  if(ct.indexOf('json') < 0) return {text: await r.text(), status: r.status};
  try { return await r.json(); }
  catch(e){ return {error:'the server did not answer with JSON (' + r.status + ')'}; }
}

api('/api/tools').then(d => {
  if(d.error) return;
  const miss = d.missing || [];
  const cores = d.cores || 1;
  $('toolnote').innerHTML = miss.length
    ? '<span class=bad>missing: ' + miss.map(esc).join(', ') + '</span>'
    : (d.local ? '<span class=good>ready</span> &middot; runs on this '
                 + 'machine, ' + cores + ' core' + (cores === 1 ? '' : 's')
                 + '; keep searches small'
      : d.sbatch ? '<span class=good>ready</span>'
                 : '<span class=bad>sbatch not found</span>')
    + (d.mmgbsa === false ? ' &middot; MM-GBSA needs AmberTools, '
                            + 'which is not installed here' : '');
  // Without AmberTools the rescoring cannot run, so it cannot be chosen.
  if (d.mmgbsa === false) {
    const box = $('mmgbsa');
    box.checked = false;
    box.disabled = true;
    box.parentNode.title = 'MM-GBSA needs AmberTools (tleap, sander, MMPBSA.py)';
    $('mmopts').classList.add('hide');
  }
});

/* ===================== receptor ====================================== */
let RECEPTOR = null;
$('recload').onclick = loadReceptor;
$('recpath').onkeydown = e => { if(e.key === 'Enter') loadReceptor(); };

/* From the user's own computer. The page is served from the cluster and
   reached through the launcher's tunnel, so the browser is on the laptop and
   the cluster cannot see its disk. The browser reads the file and posts its
   text; the server writes it under <runs>/uploads and returns a cluster path,
   and from there it is an ordinary receptor. Nothing is installed on the
   laptop: this is the file input every browser has. */
const recmsg = (t, bad) => {
  $('recinfo').innerHTML = bad ? '<span class=bad>' + esc(t) + '</span>' : esc(t);
};
$('recupload').onclick = () => $('recfile').click();
$('recfile').onchange = async () => {
  const f = $('recfile').files && $('recfile').files[0];
  if(!f) return;
  const MB = 64;
  if(f.size > MB * 1024 * 1024){
    recmsg('That file is ' + (f.size / 1e6).toFixed(0) + ' MB; the limit is '
      + MB + ' MB. Copy it to the cluster and use Browse instead.', true);
    $('recfile').value = '';
    return;
  }
  recmsg('sending ' + f.name + '…');
  let text;
  try { text = await f.text(); }
  catch(e){ recmsg('could not read that file: ' + e, true);
            $('recfile').value = ''; return; }
  const d = await api('/api/upload', {name: f.name, text: text});
  /* Cleared either way, so choosing the same file again still fires a change
     event, and a failed upload can be retried. */
  $('recfile').value = '';
  if(d.error){ recmsg(d.error, true); return; }
  $('recpath').value = d.path;
  loadReceptor();
};

async function loadReceptor(){
  const p = $('recpath').value.trim();
  if(!p) return;
  $('recinfo').textContent = 'reading…';
  const d = await api('/api/receptor?path=' + encodeURIComponent(p));
  if(d.error){ $('recinfo').innerHTML = '<span class=bad>' + esc(d.error) + '</span>';
               return; }
  gotReceptor(d);
}

function gotReceptor(d){
  RECEPTOR = d;
  $('recinfo').innerHTML = (d.fetched ? 'fetched: ' : '') + esc(d.atoms) + ' atoms, ' + d.chains.length + ' chain(s)';
  $('chainbox').classList.remove('hide');
  $('chainbox').innerHTML = '<label>Chains to keep (none ticked means all)</label>'
    + d.chains.map(c =>
      '<label class=why style="display:inline-flex;gap:6px;margin:0 16px 4px 0">'
      + '<input type=checkbox class=chk data-ch="' + esc(c.id) + '" style="width:auto">'
      + esc(c.id) + ' <span class=pill>' + esc(c.kind) + ' · ' + c.residues
      + ' res ' + c.first + ' to ' + c.last + '</span></label>').join('');
  $('chainbox').querySelectorAll('.chk').forEach(x =>
    x.onchange = () => drawPicker());
  HOT = [];
  drawChips();
  drawPicker();
}
$('recfetch').onclick = fetchCode;
$('reccode').onkeydown = e => { if(e.key === 'Enter') fetchCode(); };
async function fetchCode(){
  const c = $('reccode').value.trim().toUpperCase();
  if(!c) return;
  $('recinfo').textContent = 'fetching ' + c + ' from the RCSB...';
  const d = await api('/api/fetch?code=' + encodeURIComponent(c));
  if(d.error){ $('recinfo').innerHTML = '<span class=bad>' + esc(d.error)
    + '</span>'; return; }
  $('recpath').value = d.path;
  gotReceptor(d);
}
/* ===================== the picker ==================================== */
/* One dialog for both jobs: choosing a receptor and choosing where runs go. It
   opens on "Places": the group allocation first, then the user's scratch, then
   home, because starting in a home directory makes the user's files hard to
   find. */
let ROOTS = null, PICK = null, RUNS_ROOT = '';

async function roots(){
  if(ROOTS) return ROOTS;
  const d = await api('/api/roots');
  ROOTS = (d && d.roots) || [];
  return ROOTS;
}

function openPicker(opts){
  PICK = opts;                       // {dirs, title, start, onPick}
  $('pickertitle').textContent = opts.title;
  $('pickeruse').style.display = opts.dirs ? '' : 'none';
  $('pickernew').style.display = opts.dirs ? '' : 'none';
  $('pickernote').textContent = '';
  if(opts.start) showDir(opts.start); else showPlaces();
  const d = $('picker');
  if(d.showModal) d.showModal(); else d.setAttribute('open', '');
}
function closePicker(){
  const d = $('picker');
  if(d.close) d.close(); else d.removeAttribute('open');
}
$('pickerclose').onclick = closePicker;
$('picker').addEventListener('cancel', closePicker);

async function showPlaces(){
  await roots();
  $('crumbs').innerHTML = '<b>Places</b>';
  $('listing').innerHTML = (ROOTS || []).map(r =>
    '<button class=place data-go="' + esc(r.path) + '">\u{1F4CC} '
    + esc(r.label) + '<span class=sz>' + esc(r.path)
    + '</span></button>').join('')
    || '<div class=why style="padding:12px 16px">nowhere to start</div>';
  wireListing();
}

function crumbsFor(path){
  /* Places / scratch / user / name, every part clickable, so going back up is
     one click. */
  const parts = String(path || '').split('/').filter(Boolean);
  let at = '';
  const bits = ['<button data-go="">Places</button>'];
  parts.forEach(p => { at += '/' + p;
    bits.push('<button data-go="' + esc(at) + '">' + esc(p) + '</button>'); });
  return bits.join(' / ');
}

async function showDir(path){
  if(!path) return showPlaces();
  const q = PICK && PICK.dirs ? '/api/browse?dirs=1&path=' : '/api/browse?path=';
  const d = await api(q + encodeURIComponent(path));
  if(d.error){
    $('listing').innerHTML = '<div class=why style="padding:12px 16px">'
      + esc(d.error) + '</div>';
    return;
  }
  PICK.at = d.path;
  $('crumbs').innerHTML = crumbsFor(d.path);
  $('pickernote').innerHTML = d.writable ? ''
    : '<span class=bad>not writable by you</span>';
  /* Structures before folders, since the user is looking for a receptor. A
     scratch directory can hold hundreds of folders, and listing them first
     would push the structure files out of sight. Sizes under 512 bytes would
     round to 0 kB, so they scale with the size. */
  const fs = d.files || [];
  const size = n => n >= 1048576 ? (n / 1048576).toFixed(1) + ' MB'
                  : n >= 1024    ? Math.round(n / 1024) + ' kB'
                                 : n + ' B';
  $('listing').innerHTML =
      '<button data-go="' + esc(d.parent) + '">\u2191 up</button>'
    + fs.map(x => '<button data-pick="' + esc(x.path)
        + '">\u{1F4C4} ' + esc(x.name) + '<span class=sz>'
        + size(x.size) + '</span></button>').join('')
    + (fs.length && d.dirs.length
       ? '<div class=why style="padding:7px 16px 4px">folders</div>' : '')
    + d.dirs.map(x => '<button data-go="' + esc(x.path) + '">\u{1F4C1} '
        + esc(x.name) + '</button>').join('')
    + ((d.dirs.length || fs.length) ? ''
       : '<div class=why style="padding:12px 16px">nothing here</div>');
  wireListing();
}

function wireListing(){
  $('listing').querySelectorAll('button[data-go]').forEach(x =>
    x.onclick = () => showDir(x.getAttribute('data-go')));
  $('listing').querySelectorAll('button[data-pick]').forEach(x =>
    x.onclick = () => { const v = x.getAttribute('data-pick');
                        closePicker(); PICK.onPick(v); });
  $('crumbs').querySelectorAll('button[data-go]').forEach(x =>
    x.onclick = () => showDir(x.getAttribute('data-go')));
}

/* A folder that does not exist yet is often what is wanted, for example a new
   screen beside earlier ones. It is created where the browser is and then
   opened, so the next click is "Use this folder". */
$('pickernew').onclick = async () => {
  if(!PICK || !PICK.at) return;
  const name = prompt('Name for the new folder, inside\n' + PICK.at);
  if(name === null) return;
  const d = await api('/api/mkdir', {parent: PICK.at, name: name});
  if(d.error){
    $('pickernote').innerHTML = '<span class=bad>' + esc(d.error) + '</span>';
    return;
  }
  $('pickernote').textContent = d.existed
    ? 'that folder was already here' : 'created';
  showDir(d.path);
};

$('pickeruse').onclick = () => {
  if(!PICK || !PICK.at) return;
  const v = PICK.at;
  closePicker();
  PICK.onPick(v);
};

$('recbrowse').onclick = () => openPicker({
  dirs: false, title: 'Choose a receptor',
  start: $('recpath').value.trim(),
  onPick: v => { $('recpath').value = v; clearMissingOn('recpath');
                 loadReceptor(); }});

$('outbrowse').onclick = () => openPicker({
  dirs: true, title: 'Choose where runs go',
  start: $('outdir').value.trim() || RUNS_ROOT,
  onPick: v => useRuns(v)});

async function useRuns(path){
  const d = await api('/api/runsroot', {path: path});
  if(d.error){
    $('outmsg').innerHTML = '<span class=bad>' + esc(d.error) + '</span>';
    return;
  }
  RUNS_ROOT = d.root;
  $('outdir').value = d.root;
  clearMissingOn('outdir');
  if($('runsroot')) $('runsroot').textContent = d.root;
  $('outmsg').textContent = 'New runs go to ' + d.root
    + '. Runs already started stay where they are.';
  listRuns();
}

const chosenChains = () => Array.from(document.querySelectorAll('.chk:checked'))
  .map(x => x.getAttribute('data-ch'));

/* ===================== site ========================================== */
/* A hotspot is a chain and a residue number, held as "A:56". A bare number
   would land on whichever chain carries that number. */
let HOT = [], VD = null, SEQRES = [], SITE = null;
const keyOf = (c, r) => c + ':' + r;

function activeChains(){
  const picked = chosenChains();
  if(picked.length) return picked;
  return ((RECEPTOR && RECEPTOR.chains) || [])
    .filter(c => c.kind === 'protein').map(c => c.id);
}
function toggleHot(k){
  const i = HOT.indexOf(k);
  if(i < 0) HOT.push(k); else HOT.splice(i, 1);
  HOT.sort((a, b) => a.split(':')[0].localeCompare(b.split(':')[0])
                     || (+a.split(':')[1]) - (+b.split(':')[1]));
  $('hotspots').value = HOT.join(', ');
  drawChips(); markSeq(); stylePicker(); siteSoon();
}
function drawChips(){
  $('chips').innerHTML = HOT.length
    ? HOT.map(k => '<span class=chip data-k="' + esc(k) + '">' + esc(k)
        + ' ✕</span>').join('')
      + '<button class=ghost id=hotclear style="margin-left:4px">clear</button>'
    : '<span class=why>no residues picked yet</span>';
  $('chips').querySelectorAll('.chip').forEach(c =>
    c.onclick = () => toggleHot(c.getAttribute('data-k')));
  const cl = $('hotclear');
  if(cl) cl.onclick = () => { HOT = []; $('hotspots').value = '';
    drawChips(); markSeq(); stylePicker(); siteSoon(); };
}
/* Typing in the box stays authoritative: a list of residue numbers pasted in
   must not be ignored because it was not clicked. */
$('hotspots').oninput = () => {
  HOT = $('hotspots').value.split(/[,;\n]/).map(x => x.trim()).filter(Boolean)
    .map(x => x.indexOf(':') >= 0 ? x : ((activeChains()[0] || 'A') + ':' + x));
  drawChips(); markSeq(); stylePicker(); siteSoon();
};

async function drawPicker(){
  if(!RECEPTOR) return;
  const chs = activeChains();
  $('seqchain').textContent = chs.length > 1 ? '(' + chs.join(', ') + ')'
    : (chs[0] ? '(chain ' + chs[0] + ')' : '');
  SEQRES = [];
  for(const c of chs){
    const d = await api('/api/sequence?path='
      + encodeURIComponent(RECEPTOR.path) + '&chain=' + encodeURIComponent(c));
    (d.residues || []).forEach(r => SEQRES.push(r));
  }
  drawSeqStrip();
  if(typeof $3Dmol === 'undefined'){
    $('viewD').innerHTML = '<div class=why style="padding:14px">the 3D viewer '
      + 'did not load; pick residues from the sequence instead</div>';
    return;
  }
  const r = await fetch('/api/rawstructure?t=' + encodeURIComponent(T)
    + '&path=' + encodeURIComponent(RECEPTOR.path));
  const pdb = await r.text();
  if(!pdb || pdb.indexOf('ATOM') < 0){
    $('viewD').innerHTML = '<div class=why style="padding:14px">'
      + esc(pdb.slice(0, 200)) + '</div>';
    return;
  }
  if(!VD) VD = $3Dmol.createViewer($('viewD'), {backgroundColor:'white'});
  VD.removeAllModels();
  VD.addModel(pdb, 'pdb');
  stylePicker();
  /* Only the chains in use are clickable, and the click carries the chain it
     landed on. */
  VD.setClickable({chain: chs}, true, a => {
    if(a && a.resi != null)
      toggleHot(keyOf(a.chain || chs[0], parseInt(a.resi, 10)));
  });
  VD.zoomTo({chain: chs});
  VD.render();
  drawBox();
}
function stylePicker(){
  if(!VD) return;
  const chs = activeChains();
  VD.setStyle({}, {cartoon:{color:'#cfd8e3'}});
  VD.setStyle({chain: chs}, {cartoon:{color:'spectrum'}});
  const by = {};
  HOT.forEach(k => { const b = k.split(':'); (by[b[0]] = by[b[0]] || []).push(+b[1]); });
  Object.keys(by).forEach(c => VD.setStyle({chain:c, resi:by[c]},
    {cartoon:{color:'#1a56db'}, stick:{radius:0.22, color:'#1a56db'}}));
  VD.render();
}
/* The box, drawn where it is decided.

   ADCP docks into what is inside this box and nothing else. Drawn over the
   structure, it shows at a glance whether it covers the groove, whether it
   reaches the chain behind it, and how much of it is solvent, so the padding
   and limits above can be aimed rather than guessed.

   Two shapes: a wireframe alone disappears against a cartoon, and a solid
   alone hides the residues the box is checked against. A face at 0.12 opacity
   gives the box depth and the edges show where it stops. Both are tracked,
   because removeAllShapes() would also remove anything else the viewer is
   given. */
let BOXSHAPES = [];
function drawBox(){
  if(!VD) return;
  BOXSHAPES.forEach(sh => { try { VD.removeShape(sh); } catch(e){} });
  BOXSHAPES = [];
  const tick = $('showbox');
  if(tick && tick.checked && SITE && SITE.centre && SITE.size){
    const c = SITE.centre, z = SITE.size;
    const at = {center:{x:+c[0], y:+c[1], z:+c[2]},
                dimensions:{w:+z[0], h:+z[1], d:+z[2]}};
    // the solid first, so the edges draw over their own face
    BOXSHAPES.push(VD.addBox(Object.assign({}, at,
      {color:'#1a56db', opacity:0.12})));
    BOXSHAPES.push(VD.addBox(Object.assign({}, at,
      {color:'#1a56db', wireframe:true})));
  }
  VD.render();
}

function drawSeqStrip(){
  if(!SEQRES.length){
    $('seqstrip').innerHTML = '<span class=why>load a receptor first</span>';
    return;
  }
  // What is not there is shown too.
  //
  // This strip holds the residues the structure contains, and each row is
  // labelled with the number of its first residue. Where a loop was never
  // modelled, the numbers jump (rows starting 30, 31, 33), which would read as
  // a numbering fault rather than a gap in the structure. Those gaps are the
  // backbone breaks the preparation splits the chain at, so a user picking a
  // site should see that a residue next to their hotspot is missing.
  const rows = [];
  let cur = [], n = 0;
  SEQRES.forEach((r, i) => {
    const prev = i ? SEQRES[i - 1] : null;
    if(prev && prev.chain === r.chain && r.resnum > prev.resnum + 1){
      cur.push({gap: r.resnum - prev.resnum - 1,
                from: prev.resnum + 1, to: r.resnum - 1});
    }
    cur.push(r); n++;
    if(n === 30){ rows.push(cur); cur = []; n = 0; }
  });
  if(cur.length) rows.push(cur);

  const out = [];
  rows.forEach(row => {
    const first = row.find(x => !x.gap) || row[0];
    let k = 0;
    const body = row.map(x => {
      if(x.gap)
        return '<span class=gap title="' + x.gap + ' residue'
          + (x.gap === 1 ? '' : 's') + ' missing from the structure ('
          + x.from + (x.gap === 1 ? '' : ' to ' + x.to)
          + ' were never modelled)">\u00b7\u00b7\u00b7</span>';
      const sep = (k && k % 10 === 0) ? ' ' : '';
      k++;
      return sep + '<span class=r data-k="' + esc(x.chain + ':' + x.resnum)
        + '" title="' + esc(x.resname + ' ' + x.chain + ':' + x.resnum)
        + '">' + esc(x.code) + '</span>';
    }).join('');
    out.push('<div class=ln><span class=no>' + esc(first.chain) + ':'
      + first.resnum + '</span><span>' + body + '</span></div>');
  });
  $('seqstrip').innerHTML = out.join('');
  $('seqstrip').querySelectorAll('.r').forEach(el =>
    el.onclick = () => toggleHot(el.getAttribute('data-k')));
  markSeq();
}
function markSeq(){
  const on = {};
  HOT.forEach(k => { on[k] = 1; });
  $('seqstrip').querySelectorAll('.r').forEach(el =>
    el.classList.toggle('on', !!on[el.getAttribute('data-k')]));
}

// The site is checked as the residues are typed, not only when a button is
// pressed, because the question it answers ("the box also reaches chain B,
// include it?") must be seen before submitting. Debounced, because each check
// reads the structure on the server.
let SITE_T = null;
function siteSoon(){
  clearTimeout(SITE_T);
  SITE_T = setTimeout(() => { refreshSite(); }, 400);
}
async function refreshSite(){
  const d = await api('/api/site', {
    receptor: $('recpath').value.trim(), chains: chosenChains(),
    hotspots: $('hotspots').value, box_padding:+$('pad').value,
    box_min:+$('bmin').value, box_max:+$('bmax').value});
  SITE = d.error ? null : d;
  $('sitenote').innerHTML = d.error
    ? '<span class=bad>' + esc(d.error) + '</span>'
    : esc(d.note) + '<br>centre ' + d.centre.join(', ') + '; box '
      + d.size.map(x => x.toFixed(1)).join(' × ') + ' Å';
  // What else the box reaches, offered as a choice. Unticked, the docking
  // loads only the chains chosen above, so the peptide is placed against a gap
  // where the far wall of the groove should be; ticked, those chains are
  // loaded too. Neither is right in general: on a binder-design target the
  // other chain is the one being displaced, so the user decides.
  const others = (d.other_chains || []);
  const box = $('otherchains');
  if(!others.length){ box.classList.add('hide'); }
  else {
    box.classList.remove('hide');
    // Both numbers, because they answer different questions: how far into the
    // box the chain reaches, and how much is loaded if it is included.
    // Including a chain loads all of it, since half a chain is not a molecule.
    const list = others.map(o => 'chain ' + o.chain + ' (' + o.atoms
      + ' atoms in the box, ' + (o.chain_atoms || o.atoms)
      + ' in the chain)').join(', ');
    $('otherchainstext').innerHTML =
      '<b>The box also reaches ' + esc(list) + '.</b><br>'
      + 'Tick to load and dock against the whole of '
      + (others.length > 1 ? 'those chains' : 'that chain')
      + '; right when they are the far wall of the same groove. '
      + 'Leave it unticked when they are the partner this peptide is meant '
      + 'to displace; the docking then sees only the chains chosen above.';
  }
  // The numbers changed, so the box is redrawn.
  drawBox();
}
$('sitecheck').onclick = refreshSite;
/* Redrawing only. The box the run uses is computed on the server from the
   residues, the padding and the limits. */
if($('showbox')) $('showbox').onchange = drawBox;
// The box changes when its padding or limits change, and so do the chains it
// reaches, so those are checked again. The chain ticks decide which chains
// count as "other", which is the same question.
['pad', 'bmin', 'bmax'].forEach(id => {
  const el = $(id);
  if(el) el.addEventListener('input', siteSoon);
});
document.addEventListener('change', e => {
  if(e.target && e.target.classList && e.target.classList.contains('chk'))
    siteSoon();
});

/* ===================== peptides ====================================== */
/* The sequence editor is contenteditable, so residues can be numbered and
   coloured as they are typed. The numbers live in a ::after pseudo-element
   keyed off a data attribute; as real text they would end up in textContent
   and become part of the sequence. */
const HYDRO = new Set(['A','V','L','I','M','F','W','C','P']);
const STD = new Set('ACDEFGHIKLMNPQRSTVWY'.split(''));
let PEPS = [{name:'', sequence:''}];
/* The longest peptide ADCP will build. peptide.py holds the same number and
   the server refuses anything longer; this lets the box say so while the
   sequence is typed. */
/* Read from the page, which takes it from peptide.py, so the page and the
   server always agree. */
const LIMITS = (() => {
  try { return JSON.parse($('limits').textContent); }
  catch(e){ return {min_length: 5, max_length: 30}; }
})();
const MAXLEN = LIMITS.max_length;
const MINLEN = LIMITS.min_length;

function renderSeq(el, text, bad){
  if(!text){ el.innerHTML = ''; return; }
  const out = [];
  for(let i = 0; i < text.length; i += 10){
    const chunk = text.slice(i, i + 10);
    out.push('<span class=grp data-n="' + Math.min(i + 10, text.length) + '">'
      + Array.from(chunk).map((c, k) => {
          const j = i + k;
          const cls = bad.has(j) ? 'x' : (HYDRO.has(c) ? 'h' : 'p');
          return '<span class=' + cls + '>' + esc(c) + '</span>';
        }).join('') + '</span>');
  }
  el.innerHTML = out.join('');
}
/* The caret position, counted in characters of the plain text: the only
   measure that survives the box being re-rendered. */
function caretOffset(el){
  const s = window.getSelection();
  if(!s || !s.rangeCount || !el.contains(s.anchorNode)) return null;
  const r = s.getRangeAt(0).cloneRange();
  r.selectNodeContents(el);
  r.setEnd(s.getRangeAt(0).endContainer, s.getRangeAt(0).endOffset);
  return r.toString().length;
}
function setCaret(el, off){
  if(off == null) return;
  const walk = document.createTreeWalker(el, NodeFilter.SHOW_TEXT);
  let n, seen = 0;
  while((n = walk.nextNode())){
    if(seen + n.length >= off){
      const r = document.createRange();
      r.setStart(n, Math.max(0, off - seen));
      r.collapse(true);
      const s = window.getSelection();
      s.removeAllRanges(); s.addRange(r);
      return;
    }
    seen += n.length;
  }
  const r = document.createRange();
  r.selectNodeContents(el); r.collapse(false);
  const s = window.getSelection(); s.removeAllRanges(); s.addRange(r);
}

function drawPeps(){
  $('peps').innerHTML = PEPS.map((p, i) =>
    '<div class=ent data-i="' + i + '">'
    /* The entry number, which is what a user looks for when a message names
       peptide_3. Entries are not reordered: each peptide is an independent
       job, and the results table is ranked by score rather than by entry
       order. */
    + '<span class=entnum title="entry ' + (i + 1) + '">' + (i + 1)
    +   '</span>'
    + '<div class="col fixed"><label>Entity type</label>'
    +   '<div class=fixedval>Peptide</div></div>'
    + '<div class="col nm"><label>Name</label>'
    +   '<input type=text class=pn value="' + esc(p.name) + '" placeholder="peptide '
    +   (i + 1) + '"></div>'
    + '<div class="col sqwrap"><label>Sequence</label>'
    +   '<div class="seq ps" contenteditable="plaintext-only" spellcheck=false '
    +     'data-ph="Type or paste a peptide sequence"></div>'
    +   '<div class=seqmsg></div><div class=seqlen></div></div>'
    + '<div class=tools>'
    +   '<button class=icon data-dup="' + i + '" title="Duplicate">⋮</button>'
    +   '<button class=icon data-rm="' + i + '" title="Remove">✕</button>'
    + '</div></div>').join('');

  $('peps').querySelectorAll('.ps').forEach((el, i) => {
    renderSeq(el, PEPS[i].sequence, new Set());
    el.addEventListener('input', () => {
      const off = caretOffset(el);
      const raw = el.textContent.replace(/\s+/g, '').toUpperCase();
      PEPS[i].sequence = raw;
      const bad = new Set();
      Array.from(raw).forEach((c, k) => { if(!STD.has(c)) bad.add(k); });
      renderSeq(el, raw, bad);
      setCaret(el, off == null ? raw.length : Math.min(off, raw.length));
      // A sequence that is too long is flagged at the box, not at submission,
      // since ADCP will not build a peptide longer than MAXLEN residues.
      const over = raw.length > MAXLEN;
      el.classList.toggle('toolong', over);
      const m = el.parentElement ? el.parentElement.querySelector('.seqmsg') : null;
      if(m){
        m.classList.toggle('stop', over);
        if(over) m.textContent = 'sequence too long: ' + raw.length
                               + ' residues, and ADCP folds at most ' + MAXLEN;
        else if(m.classList.contains('stop')) m.textContent = '';
      }
      validate();
      // The guideline is per residue, so it changes as the sequence is typed.
      if(window.stepsGuide) stepsGuide();
    });
    el.addEventListener('paste', e => {
      e.preventDefault();
      const t = (e.clipboardData || window.clipboardData).getData('text');
      document.execCommand('insertText', false,
        t.split('\n').filter(l => l[0] !== '>').join('').replace(/\s+/g, ''));
    });
  });
  $('peps').querySelectorAll('.pn').forEach((el, i) => {
    el.oninput = () => { PEPS[i].name = el.value; };
  });
  $('peps').querySelectorAll('button[data-rm]').forEach(b => {
    b.onclick = () => { PEPS.splice(+b.getAttribute('data-rm'), 1);
      if(!PEPS.length) PEPS.push({name:'', sequence:''});
      drawPeps(); };
  });
  $('peps').querySelectorAll('button[data-dup]').forEach(b => {
    b.onclick = () => { const i = +b.getAttribute('data-dup');
      PEPS.splice(i + 1, 0, {name:'', sequence: PEPS[i].sequence});
      drawPeps(); };
  });
  validate();
}
$('addpep').onclick = () => { PEPS.push({name:'', sequence:''}); drawPeps();
  const b = $('peps').querySelectorAll('.ps'); b[b.length - 1].focus(); };

let VT = null;
function validate(){
  clearTimeout(VT);
  VT = setTimeout(async () => {
    const d = await api('/api/peptides', {peptides: PEPS, max_length: MAXLEN});
    if(d.error) return;
    d.peptides.forEach((p, i) => {
      const row = $('peps').querySelector('.ent[data-i="' + i + '"]');
      if(!row) return;
      const box = row.querySelector('.ps'), msg = row.querySelector('.seqmsg'),
            len = row.querySelector('.seqlen');
      box.classList.toggle('err', !p.ok && p.sequence.length > 0);
      len.textContent = p.sequence ? p.sequence.length + ' residues' : '';
      msg.textContent = (p.ok || !p.sequence) ? ''
        : p.problems.map(x => x.message).join(' ');
      /* Re-mark from the server's own verdict, which decides whether the run
         starts; the local check is only for responsiveness. */
      if(!p.ok && p.sequence && document.activeElement !== box){
        const bad = new Set([].concat.apply([],
          p.problems.map(x => x.positions || [])));
        renderSeq(box, p.sequence, bad);
      }
    });
    const n = d.peptides.filter(p => p.ok && p.sequence).length;
    $('pepcount').textContent = n + ' valid peptide' + (n === 1 ? '' : 's');
    $('submit').disabled = !d.ready;
  }, 160);
}
$('pastepep').onclick = () => $('pastebox').classList.toggle('hide');
$('pasteadd').onclick = async () => {
  const d = await api('/api/paste', {text: $('pastetext').value, max_length: MAXLEN});
  if(d.error || !d.peptides) return;
  if(PEPS.length === 1 && !PEPS[0].sequence) PEPS = [];
  d.peptides.forEach(p => PEPS.push({name: p.name, sequence: p.sequence}));
  $('pastetext').value = '';
  $('pastebox').classList.add('hide');
  drawPeps();
};
$('mmgbsa').onchange = () =>
  $('mmopts').classList.toggle('hide', !$('mmgbsa').checked);

/* ===================== estimate and steppers ========================= */
$('autofill').onclick = async () => {
  const d = await api('/api/estimate', {peptides: PEPS,
    replicas:+$('replicas').value, steps:+$('steps').value,
    cores:+$('cpus').value, mmgbsa: $('mmgbsa').checked,
    mmgbsa_poses:+$('mmposes').value, throttle:+$('throttle').value,
    // What the receptor costs to prepare, as well as what the peptides cost to
    // dock. A large glycosylated structure with rebuilt loops is a different
    // job from a 1,800-atom domain.
    receptor: $('recpath').value.trim(),
    // The chains that will actually be kept: a structure often holds two
    // copies of a complex, and costing the whole file would charge for work
    // the run never does.
    chains: chosenChains(),
    fill_loops: $('fill_loops').checked,
    box: (SITE && SITE.size) || null});
  if(d.error){
    $('estnote').innerHTML = '<span class=bad>' + esc(d.error) + '</span>';
    return;
  }
  $('cpus').value = d.cores;
  $('walltime').value = d.walltime;
  $('estnote').innerHTML = esc(d.why)
    + '<br>about ' + esc(d.end_to_end) + ' of queue time booked in total, '
    + d.total_core_hours + ' core-hours of work.'
    // The box is part of the estimate. A run whose box cannot hold its own
    // longest peptide costs as much as one that can, and returns a pose folded
    // to fit. It is shown here so it can be fixed while the padding and
    // maximum box are on screen.
    + (d.box_warning
       ? '<br><span class=bad>\u26a0 ' + esc(d.box_warning) + '</span>'
       : '');
};
/* The up and down buttons. Walltime steps in half hours and never goes below
   one minute. */
document.addEventListener('click', e => {
  const b = e.target.closest ? e.target.closest('.st') : null;
  if(!b) return;
  const bits = b.getAttribute('data-step').split(':');
  const el = $(bits[0]), d = +bits[1];
  if(bits[0] === 'walltime'){
    el.value = minToHhmm(Math.max(1, hhmmToMin(el.value) + d));
  } else {
    const lo = +(el.min || 0), hi = +(el.max || 1e9);
    el.value = Math.max(lo, Math.min(hi, (+el.value || 0) + d));
  }
});
function hhmmToMin(t){
  const b = String(t || '').split(':').map(x => parseInt(x, 10) || 0);
  while(b.length < 3) b.push(0);
  return b[0] * 60 + b[1] + (b[2] ? 1 : 0);
}
function minToHhmm(m){
  const h = Math.floor(m / 60), r = m % 60;
  return (h < 10 ? '0' : '') + h + ':' + (r < 10 ? '0' : '') + r + ':00';
}

/* ===================== submit ======================================== */
/* A required box left empty is marked and scrolled to, because a message under
   the submit button is often out of sight of the box it is about. */
function clearMissing(){
  document.querySelectorAll('.missing').forEach(e =>
    e.classList.remove('missing'));
  document.querySelectorAll('.missingnote').forEach(e => e.remove());
  document.querySelectorAll('.row.hasnotes').forEach(e =>
    e.classList.remove('hasnotes'));
}
function clearMissingOn(id){
  const e = $(id);
  if(!e) return;
  e.classList.remove('missing');
  const w = e.closest('.fieldwrap') || e.parentNode;
  const n = w && w.querySelector('.missingnote');
  if(n) n.remove();
}
function markMissing(id, why){
  const e = $(id);
  if(!e) return;
  e.classList.add('missing');
  // The note is placed out of flow, under the box it belongs to. In flow it
  // would make its row taller and shift every control beside it; the row keeps
  // room for one instead, so nothing moves.
  const wrap = e.closest('.fieldwrap');
  const note = document.createElement('div');
  note.className = 'missingnote';
  note.textContent = why;
  (wrap || e.parentNode || e).appendChild(note);
  const row = e.closest('.row');
  if(row) row.classList.add('hasnotes');
  return e;
}
function missingFields(){
  const out = [];
  if(!$('recpath').value.trim())
    out.push(['recpath', 'Choose a receptor, or fetch one by PDB id.']);
  if(!$('hotspots').value.trim())
    out.push(['hotspots', 'Pick the residues the docking box should cover.']);
  if(!PEPS.length)
    out.push(['peps', 'Add at least one peptide.']);
  if(!$('label').value.trim())
    out.push(['label', 'Give the run a name: it is what the folder and the '
                     + 'results table are called.']);
  if(!$('outdir').value.trim())
    out.push(['outdir', 'Choose where the run directory goes.']);
  return out;
}
/* Clear the complaint as soon as the box is filled, rather than at the next
   submit. */
['recpath','hotspots','label','outdir'].forEach(id => {
  const e = $(id);
  if(e) e.addEventListener('input', () => {
    if(e.value.trim()) clearMissingOn(id);
  });
});

$('submit').onclick = async () => {
  clearMissing();
  const missing = missingFields();
  if(missing.length){
    missing.forEach(([id, why]) => markMissing(id, why));
    const first = $(missing[0][0]);
    if(first){
      first.scrollIntoView({behavior: 'smooth', block: 'center'});
      if(first.focus) try { first.focus({preventScroll: true}); } catch(e){}
    }
    $('submsg').innerHTML = '<span class=bad>'
      + missing.length + (missing.length === 1 ? ' thing is' : ' things are')
      + ' still needed</span>';
    return;
  }
  $('submsg').textContent = 'submitting…';
  // A path typed into the box has to work as well as one chosen in the
  // browser. It is applied here rather than on every keystroke, so a
  // half-typed path never becomes the output folder.
  const want = $('outdir').value.trim();
  if(want && want !== RUNS_ROOT){
    const moved = await api('/api/runsroot', {path: want});
    if(moved.error){
      $('outdir').classList.add('missing');
      $('submsg').innerHTML = '<span class=bad>' + esc(moved.error)
        + '</span>';
      $('outdir').scrollIntoView({behavior:'smooth', block:'center'});
      return;
    }
    RUNS_ROOT = moved.root;
    $('outdir').value = moved.root;
    if($('runsroot')) $('runsroot').textContent = moved.root;
  }
  const d = await api('/api/create', {
    label: $('label').value || 'run', receptor: $('recpath').value.trim(),
    chains: chosenChains(), hotspots: $('hotspots').value, ph:+$('ph').value,
    minimise: $('minimise').checked, box_padding:+$('pad').value,
    keep_hetatms: $('keep_hetatms').checked,
    fill_loops: $('fill_loops').checked,
    box_min:+$('bmin').value, box_max:+$('bmax').value, peptides: PEPS,
    box_chains: $('box_chains').checked ? 'include' : 'ignore',
    poses:+$('poses').value, replicas:+$('replicas').value, steps:+$('steps').value,
    mmgbsa: $('mmgbsa').checked, mmgbsa_poses:+$('mmposes').value,
    top_n:+$('topn').value,
    cpus:+$('cpus').value, time: $('walltime').value,
    throttle:+$('throttle').value,
    account: $('account').value.trim(),
    sbatch_extra: $('sbatch_extra').value.trim()});
  if(d.error){
    $('submsg').innerHTML = '<span class=bad>' + esc(d.error) + '</span>'
      + (d.problems ? '<ul style="text-align:left">'
          + d.problems.map(x => '<li>' + esc(x) + '</li>').join('') + '</ul>' : '');
    return;
  }
  $('submsg').innerHTML = '<span class=good>submitted</span>: '
    + esc(d.rundir) + '<br>' + Object.entries(d.jobs)
        .map(([k, v]) => esc(k) + ' ' + esc(v)).join(' · ');
  listRuns();
};

/* ===================== how to cite ==================================== */
/* Scoped to the open run when there is one: a run without MM-GBSA must not be
   given an Amber citation, since a methods section that claims a calculation
   the run never did is worse than a missing reference. With no run open, it
   lists everything the tool can run and says so. */
async function showCite(){
  const scoped = RUN ? '?run=' + encodeURIComponent(RUN) : '';
  $('citescope').textContent = RUN
    ? 'what this run actually used' : 'everything this tool can run';
  $('citetext').textContent = 'loading...';
  $('citenote').textContent = 'Download for your reference manager:';
  $('citedlg').showModal();
  const d = await api('/cite' + scoped);
  $('citetext').textContent = d.text != null
    ? d.text : (d.error || 'could not load the citation list');
}
$('citebtn').onclick = showCite;
if($('citefoot')) $('citefoot').onclick = e => { e.preventDefault(); showCite(); };
$('citeclose').onclick = () => $('citedlg').close();
$('citecopy').onclick = () => copyText($('citetext').textContent, $('citetext'),
  $('citenote'));
document.querySelectorAll('[data-cite]').forEach(b => {
  b.onclick = () => {
    /* A real navigation, not a blob: the server sets the filename and the
       browser's own download handling takes it from there. */
    const u = '/cite?fmt=' + b.dataset.cite + '&dl=1&t=' + encodeURIComponent(T)
      + (RUN ? '&run=' + encodeURIComponent(RUN) : '');
    location.href = u;
  };
});

/* ===================== feedback ====================================== */
/* Email or a GitHub issue, from a dialog. There is no central server to
   receive a report, and a button that silently drops what a user wrote is
   worse than none. The dialog shows the addresses and the details every
   report needs (build, runs folder, browser), so a report can be sent even
   where the browser blocks mail links, as it does inside a notebook frame. */
async function copyText(text, node, note){
  try {
    await navigator.clipboard.writeText(text);
    note.textContent = 'Copied.';
    return;
  } catch(e){}
  /* Where the clipboard API is refused, select the text and ask the browser
     to copy the selection; failing that, the text is left selected. */
  try {
    if(node && node.select){ node.focus(); node.select(); }
    else if(node){
      const r = document.createRange(); r.selectNodeContents(node);
      const s = window.getSelection(); s.removeAllRanges(); s.addRange(r);
    }
    if(node && document.execCommand('copy')){ note.textContent = 'Copied.'; return; }
  } catch(e){}
  note.textContent = 'Selected: press Ctrl+C (Cmd+C on a Mac).';
}
$('feedback').onclick = () => {
  const details = [
    'page: ' + location.host,
    'build: ' + (($('toolnote') || {}).textContent || 'unknown'),
    'runs: ' + (RUNS_ROOT || 'not set'),
    'browser: ' + navigator.userAgent].join('\n');
  $('fbtext').value = details;
  $('fbnote').textContent = '';
  const body = '\n\n--- please leave the lines below ---\n' + details;
  $('fbmail').href = 'mailto:aadhil.haq@tamu.edu?cc=sfernando@tamu.edu'
    + '&subject=' + encodeURIComponent('GUI to ADCP Docking: feedback')
    + '&body=' + encodeURIComponent(body);
  $('fbissue').href = 'https://github.com/aadhilhaq/ADCP_docking/issues/new?body='
    + encodeURIComponent(body);
  $('fbdlg').showModal();
};
$('fbcopy').onclick = () => copyText($('fbtext').value, $('fbtext'), $('fbnote'));
$('fbclose').onclick = () => $('fbdlg').close();

/* ===================== runs and the job monitor ====================== */
/* A dot per run, and a dot per job under the open run. A running job's dot
   pulses, because a static green looks the same as a finished one.
     pulsing green  something is running
     amber          queued, waiting for the scheduler
     green          finished
     orange         a job failed, timed out or ran out of memory
     grey           cancelled, or nothing submitted yet
   The page polls only while something is live, so it does not keep calling
   squeue on a login node after everything has finished. */
const DOT_WORDS = {running:'running', queued:'queued', done:'finished',
                   failed:'failed', gone:'cancelled', new:'not submitted'};
function dot(kind){
  return '<span class="dot ' + esc(kind || 'gone') + '" title="'
    + esc(DOT_WORDS[kind] || kind || '') + '"></span>';
}
/* A pending reason in plain words: "PartitionTimeLimit" is accurate but means
   little to a user who has not met it. */
const STUCK_WHY = {
  PartitionTimeLimit: 'the walltime asked for is longer than this partition '
    + 'allows, so it will wait indefinitely. Lower the walltime, or name a '
    + 'partition that permits it',
  PartitionNodeLimit: 'more nodes than this partition allows',
  QOSMaxWallDurationPerJobLimit: 'longer than your QOS permits per job',
  InvalidAccount: 'the account on the job is not one you can use',
  InvalidQOS: 'the QOS on the job is not one you can use',
  PartitionDown: 'that partition is down',
  PartitionInactive: 'that partition is not accepting jobs',
  ReqNodeNotAvail: 'no node matching what was asked for is available',
  BadConstraints: 'the constraints asked for cannot be satisfied'
};
function stuckWhy(reason){
  return STUCK_WHY[reason] || (reason + ': this will not clear on its own');
}
let MONITOR = null;
async function monitor(){
  const d = await api('/api/monitor');
  if(d && d.runs) paintDots(d.runs);
  if(MONITOR) clearTimeout(MONITOR);
  // Only while something is live, and not too often: squeue on a busy login
  // node is not free.
  if(d && d.live) MONITOR = setTimeout(monitor, 15000);
}
function paintDots(rows){
  rows.forEach(r => {
    const cell = document.querySelector('[data-dot="' + cssEsc(r.rundir) + '"]');
    if(cell) cell.innerHTML = dot(r.dot) + ' <span class=why>'
      + esc(DOT_WORDS[r.dot] || r.dot) + '</span>';
    const panel = document.querySelector(
      '[data-jobs="' + cssEsc(r.rundir) + '"]');
    if(panel) panel.innerHTML = (r.stages || []).map(st =>
      '<div class=jobrow>' + dot(st.dot)
      + '<span class=jname>' + esc(st.stage) + '</span>'
      + '<span class=jid>' + esc(st.id) + '</span>'
      + '<span class=jst>' + esc(Object.entries(st.counts || {})
          .map(([k, n]) => n + ' ' + k.toLowerCase()).join(', ')
          || DOT_WORDS[st.dot] || '')
      + (st.why ? (st.stuck
          ? '<span class=bad>: ' + esc(stuckWhy(st.why)) + '</span>'
          : '<span class=why>: ' + esc(st.why) + '</span>') : '')
      + '</span></div>').join('')
      || '<div class=why>nothing has been submitted for this run yet</div>';
  });
}
/* An attribute selector needs its value escaped, and a run directory is a path
   full of dots and slashes. */
function cssEsc(v){
  return String(v).replace(/["\\]/g, '\\$&');
}

$('refresh').onclick = () => { listRuns(); monitor(); };
async function listRuns(){
  const d = await api('/api/runs');
  if(d.error) return;
  $('runsroot').textContent = d.root;
  RUNS_ROOT = d.root;
  // The output box starts where the server is already writing, so it is never
  // empty on a first visit.
  if(!$('outdir').value.trim()) $('outdir').value = d.root;
  // The actions are on the row: taking a run away is something done to that
  // run. Sizes are filled in after the table is drawn, because each needs a
  // walk of the directory and the table should not wait.
  $('runs').innerHTML = d.runs.length
    ? '<table><tr><th>run</th><th class=c>created</th><th class=c>peptides</th>'
      + '<th class=c>MM-GBSA</th><th class=c>state</th>'
      + '<th class=c>actions</th></tr>'
      + d.runs.map(r => '<tr><td>' + esc(r.label) + '</td><td class="why c">'
          + esc(r.created) + '</td><td class=c>' + r.peptides + '</td><td class=c>'
          + (r.mmgbsa ? 'yes' : 'no') + '</td><td class=c data-dot="'
          + esc(r.rundir) + '">' + esc(r.state)
          + '</td><td class=c style="white-space:nowrap">'
          + '<button class=ghost data-open="' + esc(r.rundir) + '">open</button> '
          // A failed run is a settings file. Every run keeps the job.json it
          // was submitted with, so a run that failed can be reopened with its
          // receptor, chains, hotspots, box and peptides, and resubmitted
          // after the one fix it needs.
          + '<button class="ghost sm" data-reuse="' + esc(r.rundir)
          + '" title="load this run\'s settings back into the form">'
          + 'reuse settings</button> '
          // Buttons, not anchors. An <a> takes the link colour and underline
          // from any rule that outranks its class and reads as a link; these
          // sit in a row of buttons and behave as buttons, and the navigation
          // is one line of script.
          + '<button class="ghost sm" data-dl="results|' + esc(r.rundir)
          + '" title="tables, report, prepared receptor and ranked poses">'
          + 'results<span class=sz data-size="' + esc(r.rundir) + ':results"></span></button> '
          + '<button class="ghost sm" data-dl="all|' + esc(r.rundir)
          + '" title="the whole run directory, grid maps included">'
          + 'all<span class=sz data-size="' + esc(r.rundir) + ':all"></span></button> '
          + '<button class=ghost data-issue="' + esc(r.rundir)
          + '" title="send us the logs for this run">report an issue</button> '
          + '<button class="ghost sm" data-folder="' + esc(r.rundir)
          + '" title="where this run is on disk">folder</button> '
          + '<button class="ghost sm danger" data-del="' + esc(r.rundir)
          + '" data-label="' + esc(r.label) + '" title="delete this run and '
          + 'everything in its folder">delete</button>'
          + '</td></tr>').join('') + '</table>'
    : '<div class=why>no runs yet</div>';
  $('runs').querySelectorAll('button[data-open]').forEach(b =>
    b.onclick = () => openRun(b.getAttribute('data-open')));
  $('runs').querySelectorAll('button[data-issue]').forEach(b =>
    b.onclick = () => openIssue(b.getAttribute('data-issue')));
  $('runs').querySelectorAll('button[data-reuse]').forEach(b =>
    b.onclick = () => reuseRun(b.getAttribute('data-reuse')));
  $('runs').querySelectorAll('button[data-folder]').forEach(b =>
    b.onclick = () => showFolder(b.getAttribute('data-folder')));
  $('runs').querySelectorAll('button[data-del]').forEach(b =>
    b.onclick = () => deleteRun(b.getAttribute('data-del'),
                                b.getAttribute('data-label')));
  // A plain navigation, so the browser saves the file under the name the
  // server gives and shows its own progress. Fetching it into memory first
  // would hold a large run in the tab.
  $('runs').querySelectorAll('button[data-dl]').forEach(b =>
    b.onclick = () => {
      const [what, dir] = b.getAttribute('data-dl').split('|');
      const s = SIZES[dir];
      if(s && s[what] && s[what].bytes > 2e8 &&
         !confirm('That is ' + s[what].human + ' across ' + s[what].files
                  + ' files, through the SSH tunnel.\n\nContinue?')) return;
      window.location = dlHref(dir, what);
    });
  d.runs.forEach(r => fillSizes(r.rundir));
  monitor();
}

let WATCH = null, RUN = null;
/* ===================== a run's settings, back in the form ============ */
/* The order matters here, so this is not a loop over the job file. The chain
   ticks do not exist until the receptor has been read (loadReceptor builds
   them from the chains it finds), so ticking waits for that, and the site
   check waits for the ticks, or it would size the box against the wrong
   chains. */
const REUSE_FIELDS = {
  recpath:'receptor', hotspots:'hotspots', label:'label', outdir:'outdir',
  pad:'box_padding', bmin:'box_min', bmax:'box_max',
  ph:'ph', poses:'poses', replicas:'replicas', steps:'steps',
  cpus:'cpus', walltime:'time', throttle:'throttle',
  mmposes:'mmgbsa_poses', topn:'top_n',
  account:'account', sbatch_extra:'sbatch_extra'
};
const REUSE_TICKS = {minimise:'minimise', keep_hetatms:'keep_hetatms',
                     fill_loops:'fill_loops', mmgbsa:'mmgbsa'};

function reuseSet(id, v){
  const el = $(id);
  if(!el || v === undefined || v === null) return;
  /* A dropdown cannot hold a value it has no option for: steps and replicas
     are fixed choices, and assigning a value that is not listed leaves the
     select empty. A run whose setting is no longer in the list keeps its own
     number, marked, rather than being rounded to something no one chose. */
  if(el.tagName === 'SELECT'
     && !Array.from(el.options).some(o => o.value === String(v))){
    const o = document.createElement('option');
    o.value = String(v);
    o.textContent = Number(v).toLocaleString() + ' (from this run)';
    el.appendChild(o);
  }
  el.value = v;
}

async function reuseRun(dir){
  const d = await api('/api/run?dir=' + encodeURIComponent(dir));
  if(d.error || !d.job){
    $('outmsg').textContent = 'could not read that run: '
      + (d.error || 'no settings were saved with it');
    return;
  }
  const j = d.job;
  for(const id in REUSE_FIELDS) reuseSet(id, j[REUSE_FIELDS[id]]);
  for(const id in REUSE_TICKS){
    const el = $(id);
    if(el && j[REUSE_TICKS[id]] !== undefined) el.checked = !!j[REUSE_TICKS[id]];
  }
  if($('box_chains')) $('box_chains').checked = j.box_chains === 'include';

  PEPS = (j.peptides && j.peptides.length)
    ? j.peptides.map(p => ({name: p.name || '', sequence: p.sequence || ''}))
    : [{name:'', sequence:''}];
  drawPeps();

  // The receptor first, because the chain ticks are drawn from it.
  if(j.receptor){
    clearMissingOn('recpath');
    await loadReceptor();
    (j.chains || []).forEach(c => {
      const el = document.querySelector('.chk[data-ch="' + c + '"]');
      if(el) el.checked = true;
    });
  }
  // Then the hotspots, which loadReceptor clears, and only then the box.
  if(j.hotspots !== undefined){
    $('hotspots').value = j.hotspots;
    HOT = String(j.hotspots).split(/[,;\n]/).map(x => x.trim()).filter(Boolean)
      .map(x => x.indexOf(':') >= 0 ? x : ((activeChains()[0] || 'A') + ':' + x));
    drawChips(); markSeq(); stylePicker();
  }
  await refreshSite();
  drawBox();
  $('outmsg').textContent = 'Settings loaded from ' + (j.label || dir)
    + '. Nothing has been submitted; change what went wrong and create the run.';
  $('recpath').scrollIntoView({behavior:'smooth', block:'start'});
}

async function openRun(dir){
  clearInterval(WATCH);
  RUN = dir;
  const show = async () => {
    const d = await api('/api/run?dir=' + encodeURIComponent(dir));
    const el = $('rundetail');
    if(d.error){ el.innerHTML = '<span class=bad>' + esc(d.error) + '</span>';
                 return; }
    monitor();
    const rows = (d.results && d.results.rows) || [];
    const mm = rows.some(r => r.mmgbsa_n);
    const tn = (d.results && d.results.top_n) || 5;
    el.innerHTML = '<div style="font-weight:600">' + esc(d.job.label || dir)
      + '</div><div class=why>' + esc(dir) + ' · ' + esc(d.job.state)
      + ' · ' + d.done + ' of ' + d.total + ' peptide(s) finished'
      + (d.job.site && d.job.site.note ? '<br>site: ' + esc(d.job.site.note) : '')
      + '</div>'
      /* What the preparation did not do, shown above the numbers it affects. A
         receptor docked without protonation still produces a results table,
         and nothing in the table says so. */
      + ((d.prep_warnings || []).length
        ? '<div class=bad style="margin:8px 0">' + d.prep_warnings.map(esc)
            .join('<br>') + '</div>' : '')
      // One row per queued job, with a dot. Filled in by the monitor, so there
      // is one squeue call for the whole page rather than one per panel.
      + '<div class=jobs data-jobs="' + esc(dir) + '" style="margin:8px 0">'
      + Object.entries(d.jobs || {}).map(([k, v]) => '<div class=jobrow>'
          + '<span class="dot queued"></span>'
          + '<span class=jname>' + esc(k) + '</span>'
          + '<span class=jid>' + esc(v.id) + '</span>'
          + '<span class=jst>' + esc(Object.keys(v.state || {}).join(', ')
              || 'checking') + '</span></div>').join('') + '</div>'
      + (rows.length
        ? '<table><tr><th class=c>#</th><th>name</th><th>peptide</th>'
          + '<th class=c>len</th>'
          + '<th class=c>poses</th>'
          + '<th class=c>ADCP best</th><th class=c>ADCP top' + tn + '</th>'
          + (mm ? '<th class=c>MM-GBSA best</th><th class=c>MM-GBSA top' + tn
                  + '</th>'
                  + '<th class=c>median</th><th class=c title="Most positive minus most '
                  + 'negative score across the scored poses. It shows how '
                  + 'much the poses disagree; it is not a statistical test '
                  + 'of the ranking.">n / '
                  + 'spread</th>' : '') + '<th class=c></th></tr>'
          // Centred, except the sequence. Numbers read as a column when their
          // digits line up under the heading. The peptide stays left-aligned
          // and monospaced, because it is read letter by letter.
          + rows.map(r => '<tr class=pick data-pep="' + esc(r.peptide) + '">'
              + '<td class=c>' + (r.rank == null ? '-' : r.rank) + '</td>'
              // The name is left-aligned for the same reason as the sequence:
              // both are read, not compared down a column.
              + '<td class=nm>' + esc(r.name || '-') + '</td>'
              + '<td class=seq>' + esc(r.peptide)
              + '</td><td class=c>' + r.length + '</td><td class=c>' + r.poses
              + '</td><td class=c>' + fmt(r.adcp_best, 1) + '</td><td class=c>'
              + fmt(r['adcp_top' + tn], 2) + '</td>'
              + (mm ? '<td class=c>' + fmt(r.mmgbsa_best, 2) + '</td>'
                      + '<td class=c>' + fmt(r['mmgbsa_top' + tn], 2) + '</td>'
                      + '<td class=c>' + fmt(r.mmgbsa_median, 2) + '</td>'
                      + '<td class=c>' + (r.mmgbsa_n || '-') + ' / '
                      + fmt(r.mmgbsa_spread, 1) + '</td>' : '')
              + '<td class=c><button class="ghost sm" data-view="'
              + esc(r.peptide) + '">view</button></td>'
              + '</tr>').join('') + '</table>'
          + (mm ? '<div class=why style="margin-top:6px">MM-GBSA is an '
                  + 'ordering within this run against this receptor. It is '
                  + 'not a K<sub>D</sub> and not a binding free energy: '
                  + 'there is no entropy term in it. The whole picomolar-to-'
                  + 'nanomolar range is 4.1 kcal/mol, which is narrower than '
                  + 'the pose spread of a single peptide here. '
                  + '<code>adcp-dock affinity</code> prints the ruler.</div>'
             : '')
        : '<div class=why>no results yet</div>')
      + ((d.outstanding && d.outstanding.length)
        ? '<div style="margin-top:11px;padding:11px;border:1px solid '
          + 'var(--line);border-radius:10px"><b>' + d.outstanding.length
          + ' peptide(s) have no result</b> <span class=why>'
          + esc(d.outstanding.slice(0, 6).map(x => x.sequence).join(', '))
          + (d.outstanding.length > 6 ? '...' : '') + '</span><br>'
          + '<button id=resume style="margin-top:8px">Run the missing ones '
          + 'again</button> <span class=why>the ones that finished are kept; '
          + 'only these go back to the queue</span>'
          + '<div class=why id=resumemsg style="margin-top:7px"></div></div>'
        : '')
      + (d.log_tail ? '<pre>' + esc(d.log_tail) + '</pre>' : '');
    const rb = $('resume');
    if(rb) rb.onclick = async () => {
      $('resumemsg').textContent = 'resubmitting...';
      const r = await api('/api/resume', {rundir: dir});
      $('resumemsg').innerHTML = r.error
        ? '<span class=bad>' + esc(r.error) + '</span>'
        : '<span class=good>' + r.todo + ' peptide(s) resubmitted</span>'
          + ', array ' + esc(r.indices)
          + (r.prepared_again ? ' (the receptor is being prepared again)' : '');
      show();
    };
    const pickRow = (tr, pep) => {
      el.querySelectorAll('tr.pick').forEach(x => x.classList.remove('on'));
      tr.classList.add('on');
      showPeptide(dir, pep);
      // The viewer is built below the table, so it is scrolled into view.
      setTimeout(() => {
        const v = $('viewcard');
        if(v && v.scrollIntoView)
          v.scrollIntoView({behavior: 'smooth', block: 'start'});
      }, 60);
    };
    el.querySelectorAll('tr.pick').forEach(tr => {
      tr.onclick = () => pickRow(tr, tr.getAttribute('data-pep'));
      const b = tr.querySelector('button[data-view]');
      if(b) b.onclick = e => { e.stopPropagation();
                               pickRow(tr, b.getAttribute('data-view')); };
    });
  };
  await show();
  WATCH = setInterval(show, 15000);
}

/* ===================== the 3D panel ================================== */
/* Adapted from the HOPE peptide pipeline's results viewer, without the
   controls that have nothing to act on here: the overlay of a co-folded
   prediction, and a choice of viewer engine.

   Note: writing the second engine's name with its asterisk followed by a slash
   inside a block comment would end the comment early. */
let VR = null, RDRAG = null, MODEL = null, PEPCHAIN = null, POSES = [], WHICH = null;
let RCOLOR = 'element', RDRAW = 'cartoon', RCONTACT = true, RPUB = false;
/* The distance on the dashed line. A contact labelled only with its two
   residue names says that they touch but not how closely, which is the
   difference between a hydrogen bond at 2.8 A and two atoms in the same
   neighbourhood at 4.4 A. */
let RDIST = true;
let RSURF = 0, LSURF = 0, RLABEL = 'auto', RSPIN = false;
let RFOCUS = [], RFOCUSFRAMED = false, RPUBFRAMED = false;
const RCUSTOM = {rec:'#4A6FA5', lig:'#D4763A', bg:'#DFE3E6'};
const SURFNAME = ['off', 'solid', 'see-through', 'mesh'];
const PUB_BOND = '#111111';
const KINDS = [['hbond','H-bond','#2B83BA'], ['ionic','Salt bridge','#F0C814'],
               ['stack','Ring stacking','#8CB366'], ['vdw','van der Waals','#808080']];
const SPECIFIC = ['ionic','stack','hbond','vdw'];
const KINDON = {hbond:1, ionic:1, stack:1, vdw:0};
const HBOND_A = 3.5, IONIC_A = 4.0, STACK_A = 5.0, VDW_A = 4.0;
const POLAR = {N:1, O:1, S:1};
const BASIC = {ARG:['NE','NH1','NH2'], LYS:['NZ'], HIS:['ND1','NE2']};
const ACIDIC = {ASP:['OD1','OD2'], GLU:['OE1','OE2']};
const RING = {PHE:1, TYR:1, TRP:1, HIS:1};
const colourOfKind = k => (KINDS.find(x => x[0] === k) || [0,0,'#888'])[2];

function carbonScheme(hex){
  try {
    const base = ($3Dmol.elementColors && $3Dmol.elementColors.defaultColors) || {};
    const map = Object.assign({}, base); map.C = hex;
    return {prop:'elem', map: map};
  } catch(e){ return 'grayCarbon'; }
}
function labelBox(col, extra){
  return Object.assign({fontSize:11, fontColor:'#14202b', backgroundColor:'#ffffff',
    backgroundOpacity:0.92, borderThickness:0.8, borderColor: col || '#9aa6b1',
    inFront:true}, extra || {});
}
function labelBare(col, extra){
  return Object.assign({fontSize:13, fontColor: col || '#14202b',
    showBackground:false, inFront:true}, extra || {});
}
function labelFor(bare, col, extra){
  const m = (RLABEL === 'auto') ? (bare ? 'plain' : 'boxed') : RLABEL;
  if(m === 'none') return null;
  return m === 'plain' ? labelBare(col, extra) : labelBox(col, extra);
}
function surfStyle(mode, hex){
  if(mode === 1) return {opacity:0.92, color:hex};
  if(mode === 2) return {opacity:0.50, color:hex};
  return {opacity:0.95, color:hex, wireframe:true, linewidth:0.5};
}
const focusKey = (c, r) => c + ':' + r;

async function showPeptide(dir, pep){
  $('viewcard').classList.remove('hide');
  WHICH = {dir: dir, pep: pep};
  const d = await api('/api/poses?dir=' + encodeURIComponent(dir)
                      + '&peptide=' + encodeURIComponent(pep));
  POSES = (d && d.poses) || [];
  $('rposepick').innerHTML = POSES.map(p =>
    '<option value="' + p.rank + '">pose ' + p.rank
    + (p.score != null ? '  ' + p.score.toFixed(1) : '')
    + (p.mmgbsa != null ? '  MM ' + p.mmgbsa.toFixed(1) : '') + '</option>').join('');
  if(!POSES.length){ $('viewmsg').textContent = 'no poses for ' + pep; return; }
  loadPose(POSES[0].rank);
}
$('rposepick').onchange = () => loadPose(+$('rposepick').value);

async function loadPose(rank){
  if(!WHICH) return;
  $('viewmsg').textContent = 'loading pose ' + rank + '…';
  const r = await api('/api/structure?dir=' + encodeURIComponent(WHICH.dir)
    + '&peptide=' + encodeURIComponent(WHICH.pep) + '&rank=' + rank
    + '&what=complex');
  if(r.error){ $('viewmsg').innerHTML = '<span class=bad>' + esc(r.error)
    + '</span>'; return; }
  const pdb = r.text != null ? r.text : r;
  if(typeof $3Dmol === 'undefined'){
    $('viewmsg').innerHTML = '<span class=bad>the 3D viewer did not load</span>';
    return;
  }
  showEngine();
  if(engineIs('molstar')){
    CONTACTS = classifyContacts(pdb, peptideChain(pdb));
    PEPCHAIN = peptideChain(pdb);
    drawContactList(RCONTACT ? CONTACTS.res : []);
    drawLid(RCONTACT ? CONTACTS : null);
    await drawMolstar(pdb);
    $('viewmsg').textContent = 'pose ' + rank + ' in Mol* - peptide is chain '
      + PEPCHAIN + ', ' + CONTACTS.res.length + ' receptor residue(s) in contact';
    return;
  }
  if(!VR){
    VR = $3Dmol.createViewer($('viewR'), {backgroundColor:'white'});
    // Label dragging, using the same helper as the HOPE peptide pipeline's
    // viewer. A figure author moves labels clear of the structure before
    // saving; the helper keeps a dragged label inside the canvas and draws its
    // leader line into the saved PNG, which pngURI() does not include because
    // the line is SVG.
    try {
      if(window.makeLabelDrag) RDRAG = window.makeLabelDrag(VR, $('viewR'), {});
    } catch(e){ RDRAG = null; }
  }
  VR.removeAllModels(); VR.removeAllShapes(); VR.removeAllLabels();
  try { VR.removeAllSurfaces(); } catch(e){}
  MODEL = VR.addModel(pdb, 'pdb');
  PEPCHAIN = peptideChain(pdb);
  RFOCUS = []; RFOCUSFRAMED = false; RPUBFRAMED = false;
  CONTACTS = classifyContacts(pdb, PEPCHAIN);
  draw(true);
  $('viewmsg').textContent = 'pose ' + rank + ': peptide is chain '
    + PEPCHAIN + ', ' + CONTACTS.res.length + ' receptor residue(s) in contact';
}

/* The peptide is the short chain, read from the file rather than assumed: the
   server writes it as chain Z, but a complex assembled another way still
   works. */
function peptideChain(pdb){
  const n = {};
  pdb.split('\n').forEach(L => {
    if(L.slice(0,4) !== 'ATOM' && L.slice(0,6) !== 'HETATM') return;
    const c = L[21]; n[c] = (n[c] || 0) + 1;
  });
  const ks = Object.keys(n);
  return ks.sort((a,b) => n[a] - n[b])[0] || 'Z';
}

let CONTACTS = {res:[], pairs:[]};
function classifyContacts(pdb, pep){
  const at = [];
  pdb.split('\n').forEach(L => {
    if(L.slice(0,4) !== 'ATOM' && L.slice(0,6) !== 'HETATM') return;
    const el = (L.slice(76,78).trim() || L.slice(12,16).trim()[0] || '').toUpperCase();
    if(el === 'H') return;
    at.push({name:L.slice(12,16).trim(), resn:L.slice(17,20).trim(), ch:L[21],
             resi:parseInt(L.slice(22,26),10), el:el,
             x:+L.slice(30,38), y:+L.slice(38,46), z:+L.slice(46,54)});
  });
  const lig = at.filter(a => a.ch === pep), rec = at.filter(a => a.ch !== pep);
  const pairs = [], byres = {}, seenPair = {};
  for(const a of rec) for(const b of lig){
    const d = Math.hypot(a.x-b.x, a.y-b.y, a.z-b.z);
    if(d > STACK_A) continue;
    const k = {};
    if(d <= VDW_A) k.vdw = 1;
    if(d <= HBOND_A && POLAR[a.el] && POLAR[b.el]) k.hbond = 1;
    if(d <= IONIC_A
       && ((BASIC[a.resn]||[]).indexOf(a.name) >= 0 && (ACIDIC[b.resn]||[]).indexOf(b.name) >= 0
        || (ACIDIC[a.resn]||[]).indexOf(a.name) >= 0 && (BASIC[b.resn]||[]).indexOf(b.name) >= 0))
      k.ionic = 1;
    if(d <= STACK_A && RING[a.resn] && RING[b.resn] && a.el === 'C' && b.el === 'C')
      k.stack = 1;
    if(!Object.keys(k).length) continue;
    const kind = SPECIFIC.find(s => k[s]) || 'vdw';
    /* One line per residue pair, not per atom pair. Two residues in contact
       touch at several atoms (a tryptophan against a tyrosine can give twenty
       qualifying pairs), and drawing each would fan parallel dashes between
       the same two side chains. The closest atom pair gives the distance and
       the most specific kind gives the colour. */
    const pk = a.ch + ':' + a.resi + '|' + b.ch + ':' + b.resi;
    const have = seenPair[pk];
    if (have){
      if (d < have.d){ have.d = d;
        have.ax = a.x; have.ay = a.y; have.az = a.z;
        have.bx = b.x; have.by = b.y; have.bz = b.z; }
      if (SPECIFIC.indexOf(kind) < SPECIFIC.indexOf(have.kind))
        have.kind = kind;
    } else {
      const rec1 = {chain:a.ch, resi:a.resi, resn:a.resn, ax:a.x, ay:a.y,
                    az:a.z, bchain:b.ch, bresi:b.resi, bresn:b.resn,
                    bx:b.x, by:b.y, bz:b.z, d:d, kind:kind};
      seenPair[pk] = rec1;
      pairs.push(rec1);
    }
    const rk = a.ch + ':' + a.resi;
    if(!byres[rk]) byres[rk] = {chain:a.ch, resi:a.resi, resn:a.resn, kinds:{},
                                d:d, x:a.x, y:a.y, z:a.z};
    byres[rk].kinds[kind] = 1;
    if(d < byres[rk].d){ byres[rk].d = d; byres[rk].x = a.x; byres[rk].y = a.y;
                         byres[rk].z = a.z; }
  }
  return {res: Object.values(byres).sort((p,q) => p.d - q.d), pairs: pairs};
}
const kindOf = r => SPECIFIC.find(k => r.kinds[k]) || 'vdw';

function baseStyle(){
  if(RDRAW === 'sticks') return {stick:{radius:0.16}};
  if(RDRAW === 'trace')  return {cartoon:{style:'trace', thickness:0.3}};
  return {cartoon:{}};
}
function colourFor(isPep){
  if(RCOLOR === 'custom') return {color: isPep ? RCUSTOM.lig : RCUSTOM.rec};
  if(RCOLOR === 'plain')  return {color:'#9fb3b6'};
  if(RCOLOR === 'chain')  return {colorscheme:'chain'};
  return {colorscheme:'default'};
}

function draw(frame){
  if(!VR || MODEL === null) return;
  VR.removeAllShapes(); VR.removeAllLabels();
  try { VR.removeAllSurfaces(); } catch(e){}
  const bs = baseStyle();
  VR.setStyle({model:MODEL},
    Object.assign({}, bs, {cartoon: Object.assign({}, bs.cartoon || {},
      colourFor(false))}));
  VR.setStyle({model:MODEL, chain:PEPCHAIN},
    {stick:Object.assign({radius:0.2}, colourFor(true))});

  const shown = CONTACTS.res.filter(r => KINDS.some(([k]) => KINDON[k] && r.kinds[k]));
  const focused = applyFocus();
  const pub = (!focused && RPUB) ? applyPublication(shown) : false;
  if(RCONTACT && !focused && !pub){
    shown.forEach(r => {
      const col = colourOfKind(kindOf(r));
      VR.setStyle({model:MODEL, chain:r.chain, resi:r.resi},
                  {stick:{radius:0.17, color:col}, cartoon:{color:col}});
      const lb = labelFor(false, col, {fontSize:10, position:{x:r.x,y:r.y,z:r.z}});
      if(lb) VR.addLabel(r.resn + r.resi, lb);
    });
    CONTACTS.pairs.forEach(pp => {
      if(!KINDON[pp.kind]) return;
      const col = colourOfKind(pp.kind);
      VR.addCylinder({start:{x:pp.ax,y:pp.ay,z:pp.az}, end:{x:pp.bx,y:pp.by,z:pp.bz},
        radius: pp.kind === 'vdw' ? 0.035 : 0.055, dashed:true, dashLength:0.25,
        gapLength:0.22, fromCap:1, toCap:1, color: col});
      // On the midpoint of the line it measures, in the contact's own colour,
      // to one decimal place: the precision a 3 A structure supports.
      if(!RDIST) return;
      const dl = labelFor(false, col, {fontSize:9,
        position:{x:(pp.ax+pp.bx)/2, y:(pp.ay+pp.by)/2, z:(pp.az+pp.bz)/2}});
      if(dl) VR.addLabel(pp.d.toFixed(1) + ' \u00c5', dl);
    });
  }
  applySurface();
  drawContactList(RCONTACT ? CONTACTS.res : []);
  drawLid(RCONTACT ? CONTACTS : null);
  if(frame && !RFOCUS.length && !RPUB) VR.zoomTo({model:MODEL});
  VR.render();
}

function applySurface(){
  if(!VR || MODEL === null) return;
  const draw1 = (mode, sel, hex) => {
    if(!mode) return;
    try { VR.addSurface($3Dmol.SurfaceType.VDW, surfStyle(mode, hex), sel); }
    catch(e){}
  };
  draw1(RSURF, {model:MODEL, chain:PEPCHAIN, invert:true}, RCUSTOM.rec);
  draw1(LSURF, {model:MODEL, chain:PEPCHAIN}, RCUSTOM.lig);
}

function applyFocus(){
  if(!RFOCUS.length || !VR || MODEL === null) return false;
  const want = {};
  RFOCUS.forEach(f => { want[focusKey(f.chain, f.resi)] = 1; });
  const mine = CONTACTS.pairs.filter(pp => want[focusKey(pp.chain, pp.resi)]);
  if(!mine.length){ RFOCUS = []; return false; }
  VR.setStyle({model:MODEL}, {cartoon:{color:RCUSTOM.bg, opacity:0.55}});
  RFOCUS.forEach(f => VR.setStyle({model:MODEL, chain:f.chain, resi:f.resi},
    {stick:{radius:0.2, colorscheme:carbonScheme(RCUSTOM.rec)}}));
  const far = {};
  mine.forEach(pp => { far[pp.bchain + ':' + pp.bresi] = pp; });
  Object.keys(far).forEach(k => {
    const pp = far[k];
    VR.setStyle({model:MODEL, chain:pp.bchain, resi:pp.bresi},
      {stick:{radius:0.2, colorscheme:carbonScheme(RCUSTOM.lig)}});
    const lb = labelFor(false, RCUSTOM.lig, {position:{x:pp.bx,y:pp.by,z:pp.bz}});
    if(lb) VR.addLabel(pp.bresn + pp.bresi, lb);
  });
  const seen = {};
  mine.forEach(pp => {
    const k = focusKey(pp.chain, pp.resi);
    if(seen[k]) return; seen[k] = 1;
    const rl = labelFor(false, RCUSTOM.rec, {position:{x:pp.ax,y:pp.ay,z:pp.az}});
    if(rl) VR.addLabel(pp.resn + pp.resi, rl);
  });
  mine.forEach(pp => {
    if(!KINDON[pp.kind]) return;
    const col = colourOfKind(pp.kind);
    VR.addCylinder({start:{x:pp.ax,y:pp.ay,z:pp.az}, end:{x:pp.bx,y:pp.by,z:pp.bz},
      radius:0.07, dashed:true, dashLength:0.25, gapLength:0.22, fromCap:1,
      toCap:1, color:col});
    if(!RDIST) return;
    const dl = labelFor(false, col, {fontSize:9,
      position:{x:(pp.ax+pp.bx)/2, y:(pp.ay+pp.by)/2, z:(pp.az+pp.bz)/2}});
    if(dl) VR.addLabel(pp.d.toFixed(1) + ' \u00c5', dl);
  });
  if(!RFOCUSFRAMED){
    RFOCUSFRAMED = true;
    const by = {};
    const add = (c, r) => { (by[c] = by[c] || []).push(r); };
    RFOCUS.forEach(f => add(f.chain, f.resi));
    Object.keys(far).forEach(k => add(far[k].bchain, far[k].bresi));
    const parts = Object.keys(by).map(c => ({model:MODEL, chain:c, resi:by[c]}));
    try { VR.zoomTo(parts.length === 1 ? parts[0] : {or: parts}); }
    catch(e){ try { VR.zoomTo({model:MODEL, chain:Object.keys(by)}); } catch(e2){} }
    VR.zoom(0.6);
  }
  return true;
}

function applyPublication(shown){
  if(!VR || MODEL === null) return false;
  VR.setStyle({model:MODEL}, {cartoon:{color:RCUSTOM.bg, opacity:0.65}});
  VR.setStyle({model:MODEL, chain:PEPCHAIN},
    {stick:{radius:0.20, colorscheme:carbonScheme(RCUSTOM.lig)}});
  shown.forEach(r => {
    VR.setStyle({model:MODEL, chain:r.chain, resi:r.resi},
      {stick:{radius:0.16, colorscheme:carbonScheme(RCUSTOM.rec)},
       cartoon:{color:RCUSTOM.bg, opacity:0.65}});
    const lb = labelFor(true, RCUSTOM.rec, {position:{x:r.x,y:r.y,z:r.z}});
    if(lb) VR.addLabel(r.resn + r.resi, lb);
  });
  const named = {};
  CONTACTS.pairs.forEach(pp => {
    if(!KINDON[pp.kind]) return;
    const k = pp.bchain + ':' + pp.bresi;
    if(named[k]) return; named[k] = 1;
    const lb = labelFor(true, RCUSTOM.lig, {position:{x:pp.bx,y:pp.by,z:pp.bz}});
    if(lb) VR.addLabel(pp.bresn + pp.bresi, lb);
  });
  CONTACTS.pairs.forEach(pp => {
    if(!KINDON[pp.kind]) return;
    VR.addCylinder({start:{x:pp.ax,y:pp.ay,z:pp.az}, end:{x:pp.bx,y:pp.by,z:pp.bz},
      radius:0.032, dashed:true, dashLength:0.30, gapLength:0.26, fromCap:1,
      toCap:1, color:PUB_BOND});
  });
  if(shown.length && !RPUBFRAMED){
    RPUBFRAMED = true;
    try { VR.zoomTo({model:MODEL, chain:shown.map(r => r.chain),
                     resi:shown.map(r => r.resi)}); } catch(e){}
  }
  return true;
}

function drawContactList(res){
  const box = $('contactbox');
  if(!box) return;
  if(!RCONTACT || !res.length){ box.innerHTML = ''; return; }
  const on = res.filter(r => KINDS.some(([k]) => KINDON[k] && r.kinds[k]));
  const count = k => res.filter(r => r.kinds[k]).length;
  box.innerHTML = '<div style="margin-bottom:6px">' + KINDS.map(([k, lab, col]) =>
      '<label class=why style="color:' + col + ';margin-right:14px">'
      + '<input type=checkbox data-kind="' + k + '"' + (KINDON[k] ? ' checked' : '')
      + ' style="width:auto"> ' + esc(lab) + ' (' + count(k) + ')</label>').join('')
    + '</div>'
    + (on.length ? '<div>' + on.map(r =>
        '<button class="pill ctcrow' + (RFOCUS.some(f => f.chain === r.chain
          && f.resi === r.resi) ? ' on' : '') + '" data-go="' + r.chain + ':'
        + r.resi + '" style="margin:0 5px 5px 0">' + esc(r.resn) + r.resi
        + ' · ' + r.d.toFixed(1) + ' Å</button>').join('') + '</div>'
      : '<span class=why>nothing in the kinds that are ticked</span>');
  box.querySelectorAll('input[data-kind]').forEach(b => b.onchange = () => {
    KINDON[b.getAttribute('data-kind')] = b.checked ? 1 : 0; draw(false); });
  box.querySelectorAll('.ctcrow').forEach(b => b.onclick = () => {
    const [c, r] = b.getAttribute('data-go').split(':');
    const key = focusKey(c, +r);
    const at = RFOCUS.findIndex(f => focusKey(f.chain, f.resi) === key);
    if(at >= 0) RFOCUS.splice(at, 1); else RFOCUS.push({chain:c, resi:+r});
    RFOCUSFRAMED = false;
    draw(false);
  });
}

/* ===================== the ligand interaction diagram ================ */
/* The peptide along the middle, each receptor residue above or below the
   peptide residue it is closest to, with lines coloured by the kind of
   interaction. Flat, labelled and printable, for a figure.

   Lanes are assigned before the canvas is sized. Two receptor residues nearest
   the same peptide residue need two lanes, and a height estimated from an
   average could clip the second outside the viewBox. */
function bestKind(links){
  for (const k of SPECIFIC) if (links.some(l => l.kind === k)) return k;
  return 'vdw';
}
function drawLid(cls){
  const box = $('diagwrap');
  if(!box) return;
  const pairs = ((cls && cls.pairs) || []).filter(pp => KINDON[pp.kind]);
  if(!pairs.length){
    box.innerHTML = '<h3>Ligand interaction diagram</h3><div class=why>'
      + 'Nothing to draw: no contacts in the kinds that are ticked.</div>';
    return;
  }
  const pep = {};
  pairs.forEach(pp => {
    const k = pp.bchain + ':' + pp.bresi;
    if(!pep[k]) pep[k] = {chain:pp.bchain, resi:pp.bresi, resn:pp.bresn};
  });
  const pepList = Object.keys(pep).map(k => pep[k])
                        .sort((a, b) => a.resi - b.resi);
  const pepAt = {};
  pepList.forEach((r, i) => { pepAt[r.chain + ':' + r.resi] = i; });

  const rec = {};
  pairs.forEach(pp => {
    const k = pp.chain + ':' + pp.resi;
    const e = rec[k] || (rec[k] = {chain:pp.chain, resi:pp.resi, resn:pp.resn,
                                   d:1e9, at:0, links:[]});
    e.links.push(pp);
    if(pp.d < e.d){
      e.d = pp.d;
      const i = pepAt[pp.bchain + ':' + pp.bresi];
      e.at = (i === undefined) ? 0 : i;
    }
  });
  const recList = Object.keys(rec).map(k => rec[k])
                        .sort((a, b) => a.at - b.at || a.d - b.d);

  const SLOT = 92, LANE = 44, RX = 31, RY = 13, PAD = RX + 8;
  const above = [], below = [];
  recList.forEach((r, i) => { (i % 2 ? below : above).push(r); });
  const lane = list => {
    const used = {};
    let deepest = 0;
    list.forEach(r => {
      let L = 0;
      while(used[L + ':' + r.at]) L++;
      used[L + ':' + r.at] = 1;
      r.lane = L;
      if(L > deepest) deepest = L;
    });
    return deepest;
  };
  const deepA = lane(above), deepB = lane(below);
  const W = PAD * 2 + Math.max(pepList.length, 1) * SLOT;
  const midY = PAD + (deepA + 1) * LANE + RY;
  const H = midY + 30 + (deepB + 1) * LANE + RY + PAD - RX + 8;
  const px = i => PAD + i * SLOT + SLOT / 2;
  above.forEach(r => { r.x = px(r.at); r.y = midY - 30 - (r.lane + 1) * LANE; });
  below.forEach(r => { r.x = px(r.at); r.y = midY + 30 + (r.lane + 1) * LANE; });

  const out = ['<svg viewBox="0 0 ' + W + ' ' + H + '" width="' + W
    + '" height="' + H + '" role=img aria-label="ligand interaction diagram">'];
  recList.forEach(r => r.links.forEach(pp => {
    const i = pepAt[pp.bchain + ':' + pp.bresi];
    if(i === undefined) return;
    out.push('<line x1="' + r.x.toFixed(1) + '" y1="' + r.y.toFixed(1)
      + '" x2="' + px(i).toFixed(1) + '" y2="' + midY + '" stroke="'
      + colourOfKind(pp.kind) + '" stroke-width="1.6" stroke-dasharray="4 3" '
      + 'opacity="0.85"><title>' + esc((r.resn || '?') + r.resi + ' to '
      + (pp.bresn || '?') + pp.bresi + '  ' + pp.d.toFixed(2) + ' A')
      + '</title></line>');
  }));
  pepList.forEach((r, i) => {
    const x = px(i);
    out.push('<rect x="' + (x - 30) + '" y="' + (midY - 15) + '" width="60" '
      + 'height="30" rx="7" fill="#fdf2e9" stroke="#D4763A" '
      + 'stroke-width="1.3"></rect>');
    out.push('<text x="' + x + '" y="' + (midY + 4.5) + '" text-anchor=middle '
      + 'font-size="12" font-weight="600" fill="#8a4620">'
      + esc((r.resn || '?') + r.resi) + '</text>');
  });
  out.push('<line x1="' + PAD + '" y1="' + midY + '" x2="' + (W - PAD)
    + '" y2="' + midY + '" stroke="#D4763A" stroke-width="1" opacity="0.35">'
    + '</line>');
  recList.forEach(r => {
    const col = colourOfKind(bestKind(r.links));
    out.push('<ellipse cx="' + r.x.toFixed(1) + '" cy="' + r.y.toFixed(1)
      + '" rx="' + RX + '" ry="' + RY + '" fill="#ffffff" stroke="' + col
      + '" stroke-width="1.6"></ellipse>');
    out.push('<text x="' + r.x.toFixed(1) + '" y="' + (r.y + 4).toFixed(1)
      + '" text-anchor=middle font-size="11.5" fill="#20262b">'
      + esc((r.resn || '?') + r.resi) + '</text>');
  });
  out.push('</svg>');
  box.innerHTML = '<h3>Ligand interaction diagram<button class=pill id=lidsave'
    + ' title="Save it as an SVG, which goes into a figure at any size">save'
    + '</button></h3><div class=diagbox>' + out.join('') + '</div>'
    + '<div class=why style="margin-top:7px">The peptide runs along the '
    + 'middle; each receptor residue sits over the one it is closest to. '
    + 'Hover a line for the distance. The ticks above filter this too.</div>';
}
/* An SVG rather than a screenshot: it is vector, so it goes into a figure at
   any size, and it is what the panel already is. */
document.addEventListener('click', e => {
  const t = e.target.closest ? e.target.closest('#lidsave') : null;
  if(!t) return;
  const svg = document.querySelector('#diagwrap svg');
  if(!svg) return;
  const copy = svg.cloneNode(true);
  copy.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  /* A white background: a transparent one takes on whatever is behind it in a
     document, and the labels are dark. */
  const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  bg.setAttribute('x', '0'); bg.setAttribute('y', '0');
  bg.setAttribute('width', '100%'); bg.setAttribute('height', '100%');
  bg.setAttribute('fill', '#ffffff');
  copy.insertBefore(bg, copy.firstChild);
  const txt = '<?xml version="1.0" encoding="UTF-8"?>\n'
    + new XMLSerializer().serializeToString(copy);
  const a = document.createElement('a');
  a.href = URL.createObjectURL(new Blob([txt], {type:'image/svg+xml'}));
  a.download = ((WHICH ? WHICH.pep : 'pose')) + '_interactions.svg';
  a.click();
});

/* ===================== the second engine ============================= */
/* Mol* as well as 3Dmol. Both are bundled, so the choice costs nothing at run
   time; Mol* is larger and takes a moment the first time it is used. It draws
   the same complex and the same contacts with a different renderer. */
let ENGINE = '3dmol', MV = null, MVBUSY = false;
try {
  const saved = localStorage.getItem('adcp.viewer');
  if(saved === 'molstar' || saved === '3dmol') ENGINE = saved;
} catch(e){}

function engineIs(name){ return ENGINE === name; }
function showEngine(){
  const mol = engineIs('molstar');
  $('viewR').classList.toggle('hide', mol);
  $('viewM').classList.toggle('hide', !mol);
  const sel = $('rengine');
  if(sel && sel.value !== ENGINE) sel.value = ENGINE;
  /* controls that belong to 3Dmol say so rather than doing nothing */
  ['rsurf','rsurfl','rpub','rlabel','rcontact'].forEach(id => {
    const el = $(id);
    if(el) el.style.opacity = mol ? 0.45 : 1;
  });
}
async function drawMolstar(pdb){
  if(typeof molstar === 'undefined'){
    $('viewmsg').innerHTML = '<span class=bad>Mol* did not load</span>';
    return;
  }
  if(MVBUSY) return;
  MVBUSY = true;
  try {
    if(!MV){
      MV = await molstar.Viewer.create($('viewM'), {
        layoutIsExpanded:false, layoutShowControls:false,
        layoutShowSequence:true, layoutShowLog:false,
        viewportShowExpand:true, collapseLeftPanel:true});
    }
    await MV.plugin.clear();
    await MV.loadStructureFromData(pdb, 'pdb');
  } catch(err){
    $('viewmsg').innerHTML = '<span class=bad>Mol* could not draw this: '
      + esc(err) + '</span>';
  } finally { MVBUSY = false; }
}
document.addEventListener('change', e => {
  if(!e.target || e.target.id !== 'rengine') return;
  ENGINE = e.target.value;
  try { localStorage.setItem('adcp.viewer', ENGINE); } catch(err){}
  showEngine();
  if(WHICH) loadPose(+($('rposepick').value || 1));
});

/* the toolbar */
document.addEventListener('click', e => {
  const b = e.target.closest ? e.target.closest(
    '#rvbar .pill[data-rcol],#rvbar .pill[data-rdraw],#rsurf,#rsurfl,'
    + '#rcontact,#rpub,#rreset,#rspin,#rshot') : null;
  if(!b || !VR) return;
  if(b.hasAttribute('data-rcol')){
    RCOLOR = b.getAttribute('data-rcol');
    document.querySelectorAll('#rvbar .pill[data-rcol]').forEach(x =>
      x.classList.toggle('on', x === b));
  } else if(b.hasAttribute('data-rdraw')){
    RDRAW = b.getAttribute('data-rdraw');
    document.querySelectorAll('#rvbar .pill[data-rdraw]').forEach(x =>
      x.classList.toggle('on', x === b));
  } else if(b.id === 'rsurf' || b.id === 'rsurfl'){
    const lg = b.id === 'rsurfl';
    const v = ((lg ? LSURF : RSURF) + 1) % 4;
    if(lg) LSURF = v; else RSURF = v;
    b.classList.toggle('on', !!v);
    b.textContent = (lg ? 'pep ' : 'rec ') + SURFNAME[v];
  } else if(b.id === 'rcontact'){
    RCONTACT = !RCONTACT;
    b.classList.toggle('on', RCONTACT);
    b.textContent = RCONTACT ? 'on' : 'off';
    if(!RCONTACT){ RFOCUS = []; RPUB = false; RPUBFRAMED = false;
      $('rpub').classList.remove('on'); $('rpub').textContent = 'off'; }
  } else if(b.id === 'rdist'){
    RDIST = !RDIST;
    b.classList.toggle('on', RDIST);
    b.textContent = RDIST ? 'on' : 'off';
  } else if(b.id === 'rpub'){
    RPUB = !RPUB; RPUBFRAMED = false;
    b.classList.toggle('on', RPUB);
    b.textContent = RPUB ? 'on' : 'off';
    if(RPUB && !RCONTACT){ RCONTACT = true;
      $('rcontact').classList.add('on'); $('rcontact').textContent = 'on'; }
  } else if(b.id === 'rreset'){
    RFOCUS = []; RPUB = false; RPUBFRAMED = false; RSPIN = false;
    $('rpub').classList.remove('on'); $('rpub').textContent = 'off';
    $('rspin').classList.remove('on');
    try { VR.spin(false); } catch(err){}
    draw(true); return;
  } else if(b.id === 'rspin'){
    RSPIN = !RSPIN; b.classList.toggle('on', RSPIN);
    try { VR.spin(RSPIN ? 'y' : false); } catch(err){}
    return;
  } else if(b.id === 'rshot'){
    try {
      const a = document.createElement('a');
      a.href = VR.pngURI();
      a.download = (WHICH ? WHICH.pep : 'pose') + '.png';
      a.click();
    } catch(err){}
    return;
  }
  draw(false);
});
document.addEventListener('input', e => {
  const t = e.target;
  if(!t || (t.id !== 'rccol' && t.id !== 'rlcol' && t.id !== 'rbcol')) return;
  if(t.id === 'rccol') RCUSTOM.rec = t.value;
  else if(t.id === 'rbcol') RCUSTOM.bg = t.value;
  else RCUSTOM.lig = t.value;
  if(t.id !== 'rbcol'){
    RCOLOR = 'custom';
    document.querySelectorAll('#rvbar .pill[data-rcol]').forEach(x =>
      x.classList.toggle('on', x.getAttribute('data-rcol') === 'custom'));
  }
  draw(false);
});
document.addEventListener('change', e => {
  if(!e.target || e.target.id !== 'rlabel') return;
  RLABEL = e.target.value;
  draw(false);
});

drawPeps();
listRuns();
showEngine();


/* ===================== what each setting means ======================= */
/* A "?" beside every label that has an entry, and one panel that shows it.
   Added here rather than written into the markup, so a field cannot be moved,
   renamed or added and lose its explanation: the list is the JSON island, and
   each entry gets a marker wherever its label is. */
const HELP = (() => {
  try { return JSON.parse($('helpdata').textContent); }
  catch(e){ return {}; }
})();

function showHelp(key){
  const h = HELP[key];
  if(!h) return;
  $('helptitle').textContent = h.title;
  $('helpbody').textContent = h.body;
  $('helpbox').classList.remove('hide');
}
function hideHelp(){ $('helpbox').classList.add('hide'); }

function addHelpMarkers(){
  Object.keys(HELP).forEach(key => {
    // The label that names this field, however it is written: for= on a label,
    // or a label wrapping the input (the checkboxes use the latter).
    let lab = document.querySelector('label[for="' + key + '"]');
    if(!lab){
      const el = $(key);
      lab = el && el.closest ? el.closest('label') : null;
    }
    if(!lab || lab.querySelector('.q')) return;
    const b = document.createElement('button');
    b.type = 'button'; b.className = 'q'; b.textContent = '?';
    b.title = 'what is this?';
    b.setAttribute('aria-label', 'what does ' + HELP[key].title + ' mean?');
    b.onclick = e => { e.preventDefault(); e.stopPropagation(); showHelp(key); };
    lab.appendChild(b);
  });
}

$('helpclose').onclick = hideHelp;
document.addEventListener('keydown', e => { if(e.key === 'Escape') hideHelp(); });
addHelpMarkers();

/* ============ taking a run away, and reporting one that failed ========== */
/* On the run's row, where the user is when they want it. */

function dlHref(dir, what){
  return '/download?what=' + what + '&run=' + encodeURIComponent(dir)
    + '&t=' + encodeURIComponent(T);
}

/* The size before the click: downloading everything from a finished run can be
   tens of megabytes through an SSH tunnel. Requested per row after the table
   is drawn, so the page does not wait on each directory walk. */
const SIZES = {};
async function fillSizes(dir){
  if(SIZES[dir] === undefined){
    SIZES[dir] = null;                      // in flight; do not ask twice
    const d = await api('/api/sizes?run=' + encodeURIComponent(dir));
    SIZES[dir] = (d && !d.error) ? d : false;
  }
  const d = SIZES[dir];
  if(!d) return;
  ['results', 'all'].forEach(what => {
    const el = document.querySelector(
      '[data-size="' + cssEsc(dir + ':' + what) + '"]');
    if(el && d[what]) el.textContent = ' ' + d[what].human;
  });
}

/* Where the run is. Opened in the file manager when the interface runs on
   this computer; otherwise the folder is on another machine (a cluster, a
   Colab session), so its path is shown with how to reach it. */
async function showFolder(dir){
  const r = await api('/api/open_folder', {rundir: dir});
  if(r.error){ alert(r.error); return; }
  $('foldertext').textContent = r.path;
  $('folderhint').textContent = r.opened
    ? 'Opened in your file manager.' : (r.hint || '');
  $('foldernote').textContent = '';
  $('folderdlg').showModal();
}
$('foldercopy').onclick = () => copyText($('foldertext').textContent,
  $('foldertext'), $('foldernote'));
$('folderclose').onclick = () => $('folderdlg').close();

/* Deleting is permanent, so the question names the folder and its size, and
   the server refuses while any of the run's jobs is queued or running. */
async function deleteRun(dir, label){
  const s = SIZES[dir];
  const size = (s && s.all) ? ' (' + s.all.human + ' in ' + s.all.files
    + ' files)' : '';
  if(!confirm('Delete the run "' + (label || dir) + '"?\n\n' + dir + size
              + '\n\nThe folder and everything in it are removed. This '
              + 'cannot be undone.')) return;
  const r = await api('/api/delete', {rundir: dir});
  if(r.error){ alert('Not deleted: ' + r.error); return; }
  delete SIZES[dir];
  if(RUN === dir){
    RUN = null;
    clearInterval(WATCH);
    $('rundetail').innerHTML = '';
    $('viewcard').classList.add('hide');
  }
  listRuns();
}

/* One button; the form appears when it is wanted. */
function openIssue(dir){
  const box = $('issuebox');
  $('issuerun').textContent = dir;
  $('issuemsg').textContent = '';
  $('issuenote').value = '';
  try { $('issuemail').value = localStorage.getItem('adcp_email') || ''; }
  catch(e){}
  box.classList.remove('hide');
  $('issuenote').focus();
  $('issuesend').onclick = async () => {
    const mail = $('issuemail').value.trim();
    try { if(mail) localStorage.setItem('adcp_email', mail); } catch(e){}
    $('issuesend').disabled = true;
    $('issuemsg').textContent = 'sending...';
    const r = await api('/api/report', {rundir: dir, email: mail,
      note: $('issuenote').value,
      extra: {browser: navigator.userAgent, page: location.host}});
    $('issuesend').disabled = false;
    if(r.sent){
      $('issuemsg').innerHTML = '<span class=good>Thank you. '
        + esc(r.note) + '</span>';
      return;
    }
    // Not silently: the note stays on screen and explains how to send it by
    // hand, since a button that drops what a user wrote is worse than none.
    $('issuemsg').innerHTML = '<span class=bad>'
      + esc(r.error || 'could not send') + '</span><br>'
      + esc(r.fallback || '') + '<br><a class=btnlink href="'
      + dlHref(dir, 'results') + '">download the logs and send them</a>';
  };
}
if($('issueclose')) $('issueclose').onclick = () =>
  $('issuebox').classList.add('hide');

/* ===================== moving the labels ============================= */
/* Delegated, because the toolbar is rebuilt whenever a pose is drawn, and a
   handler bound to the old button would be lost with it. */
document.addEventListener('click', e => {
  const t = e.target && e.target.closest
          ? e.target.closest('#rdrag,#rdragoff') : null;
  if(!t) return;
  if(!RDRAG){
    const m = $('viewmsg');
    if(m) m.textContent = 'label dragging could not start; the 3D viewer '
                        + 'itself may not have loaded';
    return;
  }
  if(t.id === 'rdragoff'){
    // A way back to the default layout without reloading the run.
    RDRAG.reset();
    const b = $('rdrag');
    if(b) b.classList.remove('on');
    return;
  }
  const on = !RDRAG.isEnabled();
  if(on) RDRAG.enable(); else RDRAG.disable();
  t.classList.toggle('on', on);
});

/* ============ what the guideline says for this peptide ================== */
/* ADCP's step budget should scale with the peptide (the guideline is 1 to 3
   million steps per amino acid). It is worked out here from the longest
   sequence entered, and the page says when the chosen number is below the
   published minimum. */
window.stepsGuide = stepsGuide;
function stepsGuide(){
  const el = $('stepsguide');
  if(!el) return;
  const longest = (PEPS || []).reduce(
    (n, p) => Math.max(n, (p.sequence || '').length), 0);
  if(!longest){
    el.textContent = 'ADCP suggests 1-3 million steps per residue. '
                   + 'Add a peptide and this says what that comes to.';
    return;
  }
  const lo = longest * 1e6, hi = longest * 3e6;
  const got = +$('steps').value;
  const f = n => n.toLocaleString();
  el.innerHTML = 'Longest peptide is <b>' + longest + '</b> residues, so the '
    + 'guideline is <b>' + f(lo) + '</b> to <b>' + f(hi) + '</b> steps'
    + (got < lo
       ? ': <span class=bad>' + f(got) + ' is below that</span>. An '
         + 'under-sampled run does not fail; it returns a ranked pose list '
         + 'that looks like a good one.'
       : ': <span class=good>' + f(got) + ' is in range</span>.');
}
['steps', 'replicas'].forEach(id => {
  const e = $(id);
  if(e) e.addEventListener('change', stepsGuide);
});
stepsGuide();
