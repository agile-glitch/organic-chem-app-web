/* MolView — the Molecule tab: a 2D drawing and a 3D model of one molecule, with statistics for whatever the pointer
   is over (atoms and bonds) and a side panel for the whole molecule. This file only joins the pieces; the facts and
   every caveat's wording come from js/molinfo.js (MolInfo), so all views say the same thing.

     window.MolView.show(textOrGraph, {view: '2d'|'3d', source}) → Promise<{ok, error?}>
         text: a name, SMILES, CAS number or MOL/SDF text; graph: a Chem drawing (e.g. the Sketcher's).

   Pieces (every index below is a RECORD index, shared by 2D, 3D and the panel: the molecule's atoms 0..n-1, then
   one H per implicit H, grouped by parent in ascending parent order, as RDKit's add_hs numbers them):
   - input → RDKit molecule: Chem.searchMolecule resolves names, CAS numbers and library entries (its `via`,
     `matched` and `warn` are shown, because a misspelt name silently becomes the closest compound); a SMILES goes to
     RDKit directly (it keeps isotopes, radicals and bracket H); a SMILES RDKit refuses is an error with RDKit's own
     reason, never replaced by another reader's guess (the app's lenient reader would drop bracket H counts and show
     a different compound). A name's drawing goes through a molfile. MOL/SDF text (a file, or pasted into the box)
     goes to RDKit as it is: its first line is the title, even when blank, so it is never trimmed (except one extra
     blank line before the header of text RDKit refuses, dropped with a note). A file that is not MOL/SDF gets a
     file error, not a failed name lookup.
   - salts and mixtures: RDKit get_frags; each component is its own molecule (largest first by default).
   - MolInfo.analyse(R, M) → the record: atoms 0..n-1 = M's atoms, then one H per implicit H grouped by parent. It runs
     before the 2D layout, because RDKit's InChI reads E/Z from 2D coordinates and a layout made here must not decide it.
   - 2D: M's RDKit molblock (sketch coordinates when there are any, else RDKit's layout; CoordGen for macrocycles) at
     RDKit's 1.5 bond length → Mol2DView; with "Show H" the H come from RDKit add_hs_in_place and are matched to the
     record's H by parent. Both are checked against the record (elements, bond list) before anything is drawn.
   - 3D: Mol3DEngine (OpenChemLib, MMFF94) builds from the same molblock; its topology (elements, H parents, bond
     set, charges) is checked against the record and a model that does not match is never shown. PubChem's 3D
     conformer (PubChem3D, only on a click) arrives already mapped to record order and checked there. A 3D MOL/SDF
     file with every H also gives its own coordinates as a conformer ("File": as given, no energy), mapped and checked
     the same way. Above 999 atoms (the V2000 molfile limit) no 3D model is built.
   - three.js and OpenChemLib load the first time the 3D view is used (RDKitLoad.script), not with the page.
   Accuracy rules applied here: a value that is not reliable is labelled or not shown (MolInfo decides per value);
   an unoptimised geometry shows no measured numbers; a name read as the "closest" name is flagged; a double bond whose
   geometry a name leaves open is passed to RDKit as "either" (a molfile's 2D drawing would otherwise fix it as E
   or Z); stereochemistry the input leaves open is named as the 3D model has it (read back from its coordinates);
   PubChem and file geometry is labelled as such; every 3D value comes from the conformer on screen. */
(() => {
  'use strict';
  const C = window.Chem, MI = window.MolInfo;
  const $ = id => document.getElementById(id);
  const page = $('page-mol');
  if (!page || !C || !MI) return;

  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const MINUS = '−';
  const num = (v, dp) => (v < 0 ? MINUS + (-v).toFixed(dp) : v.toFixed(dp));
  const signed = (v, dp) => { const r = Number(v.toFixed(dp)); return r === 0 ? (0).toFixed(dp) : (r > 0 ? '+' : MINUS) + Math.abs(r).toFixed(dp); };
  const EL_NUM = new WeakMap();   // per-element atom numbers (C1, C2, N1 …), as MolInfo names atoms
  const elNum = rec => { if (!EL_NUM.has(rec.atoms)) { const c = {}; EL_NUM.set(rec.atoms, rec.atoms.map(a => (c[a.el] = (c[a.el] || 0) + 1))); } return EL_NUM.get(rec.atoms); };
  const atomName = (rec, i) => rec.atoms[i].el + elNum(rec)[i];
  const THREE_SRC = 'vendor/three/three.js';
  const MAX_3D_ATOMS = 999;      // a V2000 molfile (the 3D builder's input and every 3D export) holds at most 999 atoms
  const NAMER_MAX_HEAVY = 100;   // the app's namer takes seconds (and more) above this size

  /* how each kind of value is marked (MolInfo's row kinds) */
  const KIND = {
    exact: { tag: '', title: 'follows exactly from the structure' },
    convention: { tag: 'conv.', title: 'a convention or definition (bookkeeping), not a measurement' },
    model: { tag: 'model', title: "a model's estimate, not a measured value" },
    table: { tag: 'table', title: 'a tabulated average or constant, not a value for this molecule' },
    '3d': { tag: '3D', title: 'measured on the 3D model shown (this one conformer)' },
    unchecked: { tag: 'unchecked', title: "the app's namer wrote this name, but it could not be read back to confirm it describes exactly this structure" },
  };
  const LABEL_TEXT = { none: '', element: 'element', index: 'atom number', cip: 'R/S', hybridization: 'hybridization',
    formal: 'formal charge', partial: 'Gasteiger–Marsili partial charge (model)', oxidation: 'oxidation state', lonepairs: 'lone pairs',
    electrons: 'electrons — each lone pair as two dots (the bonds are the bonding pairs)',
    pka: 'pKa of each acidic site and pKaH of each basic site, in THIS molecule',
    steric: 'steric crowding at the reacting atom of each fragment, in THIS molecule (open / moderately hindered / hindered / very hindered)',
    electronic: 'electronic effects at the reacting atom of each fragment, in THIS molecule (partial charge, ring Σσ, conjugation)',
    taft: 'Taft parameters of the groups on the reacting atom of each fragment, in THIS molecule, shown as Σσ*/Es: Σσ* (polar, > 0 pulls electrons, < 0 pushes) / Es of the bulkiest group (steric, more negative = bulkier); CH3 = 0 for both, – = no value' };

  /* ================================================================ state */
  const S = {
    R: null, inited: false,
    exportBW: false,       // the export dialog's last choices: the picture in black and white (true) or in colour,
    exportTransparent: false,   // and a transparent background (true) or a white one
    textbook2d: true,      // the 2D drawing in textbook style (all black) or with coloured heteroatoms
    gen: 0,                // the molecule shown: changes only when a new input has been read (build results are
    req: 0,                // tied to it); req counts Show requests, so an older one still waiting for RDKit gives way
    fileBlock: null,       // a 3D MOL/SDF input's first record, as read (for the "File" conformer)
    closed: new Set(),     // side-panel sections the user collapsed (by title), kept across re-renders
    info: null,            // how the input was read: {source, via, matched, warn, notes: [{kind, html}]}
    mols: [],              // RDKit objects owned by the current input (deleted when the next one is shown)
    full: null, comps: [], cs: null, combined: null,
    view: '2d', showH2d: false, showH3d: true, labels: 'none', style: 'ballstick', surface: 'none', spin: false,
    measureN: 0, picks: [], measures: [], sel: null, group: null, notesOpen: false,
    sort: { atoms: { key: 'i', dir: 1 }, bonds: { key: 'k', dir: 1 } },
  };
  let v2 = null, v2info = null;               // Mol2DView and what it currently shows {cs, withH, N, M}
  let v3 = null, v3cs = null;                 // Mol3DView and the component whose topology it holds
  let threeLoading = null, surfaceResult = null, contextLost = false;

  /* ================================================================ small DOM helpers */
  const show = (node, on) => { if (node) node.hidden = !on; };
  function download(name, blob) {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 1500);
  }
  const slug = s => (String(s || 'molecule').normalize('NFKD').replace(/[^\w\-]+/g, '_').replace(/^_+|_+$/g, '').slice(0, 48) || 'molecule');
  /* Export file names: "SMILES (common name, IUPAC name)". A file name cannot hold / \ : * ? " < > | on Windows, and
     SMILES use / and \ for E/Z, so those get look-alike characters (∕ ⧵ ꞉ ✱ ？ ＂ ‹ › ǀ): the name reads as the SMILES
     but is not one to paste back. Names are only added when the app has them; the whole name is kept under 180
     characters (Windows paths are limited to 260). */
  const FILE_SAFE = { '/': '∕', '\\': '⧵', ':': '꞉', '*': '✱', '?': '？', '"': '＂', '<': '‹', '>': '›', '|': 'ǀ' };
  const fileSafe = s => String(s).replace(/[\/\\:*?"<>|]/g, ch => FILE_SAFE[ch]).replace(/[\x00-\x1f]/g, '').replace(/[. ]+$/g, '');
  function exportName(cs) {
    const smi = (cs.record && cs.record.smiles && cs.record.smiles.isomeric) || '';
    const nm = cs.names || {}, names = [nm.common, nm.iupacChecked === false ? null : nm.iupac].filter(Boolean);   // an unchecked name only appears where it is tagged
    let name = smi ? fileSafe(smi) : slug(titleOf(cs));
    if (names.length) name += ' (' + fileSafe(names.join(', ')) + ')';
    if (name.length > 180) {
      name = name.slice(0, 177);
      const open = (name.match(/\(/g) || []).length, close = (name.match(/\)/g) || []).length;
      if (open > close) name += '…)';
    }
    return name || 'molecule';
  }
  function setStatus(text, busy) {
    const st = $('molStatus');
    if (!text) { st.hidden = true; st.innerHTML = ''; return; }
    st.hidden = false;
    st.innerHTML = (busy ? '<span class="mv-spin" aria-hidden="true"></span>' : '') + esc(text);
  }
  // the status line of the molecule still shown (after a request that did not replace it)
  function restoreStatus() {
    const B = S.cs && S.cs.b3d;
    if (S.view === '3d' && B) setStatus(B.statusText || '', !B.done);
    else setStatus('');
  }
  // let the browser paint (a status line) before a long synchronous step
  const nextPaint = () => new Promise(res => {
    let done = false;
    const go = () => { if (!done) { done = true; res(); } };
    try { requestAnimationFrame(() => setTimeout(go, 0)); } catch (e) {}
    setTimeout(go, 60);                            // rAF does not run in a hidden tab
  });

  /* ================================================================ input → RDKit molecule */
  function smallestRing(g, bond) {                   // ring size through `bond` in a Chem graph, 0 if not in a ring
    const adj = new Map(g.atoms.map(a => [a.id, []]));
    for (const b of g.bonds) if (b !== bond) { adj.get(b.a).push(b.b); adj.get(b.b).push(b.a); }
    const dist = new Map([[bond.a, 0]]), q = [bond.a];
    while (q.length) {
      const x = q.shift();
      if (x === bond.b) return dist.get(x) + 1;
      for (const y of adj.get(x)) if (!dist.has(y)) { dist.set(y, dist.get(x) + 1); q.push(y); }
    }
    return 0;
  }
  /* A Chem drawing as a molfile RDKit reads exactly: radicals as M  RAD (Chem.toMolfile writes none), and a double
     bond whose geometry the name left open (ezUnspec) marked "either" (3) — otherwise RDKit would read E or Z from
     the 2D layout (checked: "but-2-ene" came out as (E) without this). Ring double bonds in rings < 8 have no E/Z. */
  function graphMolfile(g) {
    const lines = C.toMolfile(g, '').split('\n'), na = g.atoms.length;
    const deg = new Map(g.atoms.map(a => [a.id, 0]));
    g.bonds.forEach(b => { deg.set(b.a, deg.get(b.a) + 1); deg.set(b.b, deg.get(b.b) + 1); });
    g.bonds.forEach((b, k) => {
      if (b.order !== 2 || !b.ezUnspec || b.stereo || deg.get(b.a) < 2 || deg.get(b.b) < 2) return;
      const rs = smallestRing(g, b);
      if (rs && rs < 8) return;
      const L = lines[4 + na + k];
      if (L) lines[4 + na + k] = L.slice(0, 9) + '  3' + L.slice(12);
    });
    const rad = [];
    g.atoms.forEach((a, i) => { if (a.radical) rad.push('M  RAD  1' + String(i + 1).padStart(4) + '   2'); });
    let t = lines.join('\n');
    if (rad.length) t = t.replace('M  END', rad.join('\n') + '\nM  END');
    return t;
  }
  // RDKit.js keeps a molblock's H atoms unless told otherwise (a SMILES loses them anyway); the record wants them
  // implicit, except isotopic H (D, T), which RDKit keeps as atoms
  function rdMol(text) {
    const R = S.R;
    let m = null;
    try { m = C.looksLikeMolfile(text) ? R.get_mol(text, JSON.stringify({ removeHs: true })) : R.get_mol(text); } catch (e) { m = null; }
    if (m && ((m.is_valid && !m.is_valid()) || !m.get_num_atoms())) { m.delete(); m = null; }
    return m;
  }
  const looksLikeSmiles = t => /^[A-Za-z0-9@+\-\[\]()=#$%\/\\.:*]+$/.test(t);
  /* SDF: the first record; MOL: as is. Header lines 1-3 are never touched: line 1 is the title even when it is blank
     (RDKit, Open Babel and most exporters write a blank title), and removing it shifts the counts line. The split
     consumes the $$$$ line with its line end, so a later record's blank title would survive too. */
  function firstRecord(text) {
    const parts = String(text).split(/\r?\n\$\$\$\$[^\n]*(?:\r?\n|$)/).filter(p => /M {2}END|V2000|V3000/.test(p));
    return { block: parts[0] || String(text), count: parts.length };
  }
  /* RDKit's own reason for refusing a SMILES: its error log, captured for this one parse (RDKit.js logging stays
     off otherwise). Only the useful lines, with a plain-words version of the common valence and aromaticity errors. */
  function rdkitReason(text) {
    const R = S.R;
    let h = null, log = '';
    try { h = R.set_log_capture ? R.set_log_capture('rdApp.error') : null; } catch (e) { h = null; }
    try {
      if (h) h.clear_buffer();
      let m = null;
      try { m = R.get_mol(text); } catch (e) { m = null; }
      if (m) m.delete();
      log = h ? String(h.get_buffer() || '') : '';
    } finally { if (h) { try { h.delete(); } catch (e) {} } }
    const lines = [...new Set(log.split(/\r?\n/).map(l => l.replace(/^\[\d\d:\d\d:\d\d\]\s*/, '').trim())
      .filter(l => l && !/^~*\^$/.test(l) && l !== text && !/^SMILES Parse Error: Failed parsing SMILES/.test(l)))];
    const out = [];
    for (const l of lines) {
      let x;
      if ((x = /Explicit valence for atom # (\d+) ([A-Za-z]+), (\d+), is greater than permitted/.exec(l))) {
        out.push(`atom ${+x[1] + 1} (${x[2]}) would have ${x[3]} bonds (counting its H), more than ${x[2]} can have with the charge written` +
          (/[+-]/.test(text) ? '' : ' (for an ion, write the charge in the brackets, e.g. [NH4+])'));
      } else if ((x = /Can't kekulize mol\.\s*Unkekulized atoms:\s*([\d ]+)/.exec(l))) {
        out.push(`the aromatic ring at atoms ${x[1].trim().split(/\s+/).map(v => +v + 1).join(', ')} cannot be given alternating double bonds (an aromatic N that carries H is written [nH])`);
      } else out.push(l.replace(/^SMILES Parse Error:\s*/, '').replace(/\s*:\s*$/, ''));   // "…around position 10:" introduces the echo removed above
    }
    return out.join('; ');
  }
  function smilesRefusal(text) {
    const why = rdkitReason(text);
    return `RDKit cannot read “${text}” as SMILES` + (why ? ': ' + why : '') + '. Nothing else is shown in its place.';
  }

  // a molfile's counts line (V2000 or V3000): where line 4 of a MOL record must be
  const isCountsLine = l => /^\s*\d+\s+\d+.*V[23]000\s*$/.test(String(l || ''));
  // → {mol, info, block?} or throws Error(user-readable message)
  function readInput(input, source) {
    const info = { source: source || 'typed', via: '', matched: '', warn: '', notes: [] };
    if (input && typeof input === 'object' && Array.isArray(input.atoms)) {
      if (!input.atoms.length) throw new Error('The drawing is empty.');
      const m = rdMol(graphMolfile(input));
      if (!m) throw new Error('RDKit could not read this drawing (check valences and charges).');
      info.via = source === 'sketcher' ? 'the Sketcher drawing' : 'a drawing';
      return { mol: m, info };
    }
    const raw = String(input || '').replace(/^﻿/, '');
    const text = raw.trim();                        // names, SMILES and CAS numbers only: a molfile keeps its lines
    if (!text) throw new Error(source === 'file' ? 'This file is empty (it has no text, or only blank lines), so there is no molecule to read.'
      : 'Type a SMILES, common name or IUPAC name first (or open a MOL/SDF file).');
    if (C.looksLikeMolfile(text)) {
      let { block, count } = firstRecord(raw.replace(/\s+$/, ''));
      let m = rdMol(block);
      // One extra blank line before the header (as copied from some web pages) moves every header line down one, so
      // RDKit refuses the text. Only then, and only when the counts line sits one line below its place, the text is
      // read once more without that first line, and the page says so. (Otherwise the first line is never removed:
      // it is the title, even when blank.)
      if (!m) {
        const L = block.split(/\r?\n/);
        if (L.length > 5 && !L[0].trim() && isCountsLine(L[4]) && !isCountsLine(L[3])) {
          const alt = L.slice(1).join('\n');
          m = rdMol(alt);
          if (m) {
            block = alt;
            info.notes.push({ kind: 'warn', html: `The ${source === 'file' ? 'file' : 'MOL text'} starts with an extra blank line before its three header lines, so it was read without that line.` });
          }
        }
      }
      if (!m) throw new Error(source === 'file' ? 'RDKit could not read this MOL/SDF file.' : 'RDKit could not read this MOL text.');
      info.via = source === 'file' ? 'MOL/SDF file' : source === 'pasted' ? 'pasted MOL text' : 'MOL text';
      if (count > 1) info.notes.push({ kind: 'info', html: `The SDF holds ${count} records; the first one is shown.` });
      return { mol: m, info, block };
    }
    // A file that is not MOL/SDF: its text is not a name to look up (the resolver's "not recognised as a name …"
    // would read as if one had been typed). A file of one line (a SMILES or a name saved as text) is read as typed.
    if (source === 'file') {
      const lines = text.split(/\r?\n/).map(l => l.trim()).filter(Boolean);
      const why = 'This file is not a MOL or SDF file (it has no V2000 or V3000 connection table and no “M\u00a0\u00a0END” line)';   // no-break spaces: HTML shows both
      if (lines.length !== 1) throw new Error(why + '.');
      try { return readText(lines[0], info); } catch (e) { throw new Error(`${why}, and its one line is not a structure the app can read: ${e.message}`); }
    }
    return readText(text, info);
  }
  // a name, SMILES or CAS number (typed, or the one line of a text file) → {mol, info}
  function readText(text, info) {
    // a long SMILES that RDKit reads: the app's name resolver is not needed (it is slow on long text, and for a
    // SMILES its reading is discarded anyway). Short text still goes to the resolver first: "CO" or "PCC" can be
    // a name or acronym as well as a SMILES, and the note below says so.
    if (text.length > 60 && looksLikeSmiles(text)) {
      const m = rdMol(text);
      if (m) { info.via = 'SMILES'; return { mol: m, info }; }
    }
    let r = null, err = null;
    try { r = C.searchMolecule(text); } catch (e) { err = e; }
    if (r && r.via === 'SMILES') {
      const m = rdMol(text);
      if (m) { info.via = 'SMILES'; return { mol: m, info }; }
      // RDKit refuses it (e.g. [NH4]: an impossible valence). No structure honestly fits, so none is shown: the
      // app's lenient reader would drop the bracket H counts and show a different compound (ammonia).
      throw new Error(smilesRefusal(text));
    }
    if (r && r.graph && r.graph.atoms.length) {
      const m = rdMol(graphMolfile(r.graph));
      if (!m) throw new Error(`"${text}" was read as ${r.matched || r.via}, but RDKit could not read that structure.`);
      info.via = r.via; info.matched = r.matched || ''; info.warn = r.warn || '';
      // the same text may also be a valid SMILES for something else: say so rather than pick silently
      if (looksLikeSmiles(text)) {
        const alt = rdMol(text);
        if (alt) {
          const a = alt.get_smiles(), b = m.get_smiles();
          alt.delete();
          if (a !== b) info.notes.push({ kind: 'warn', html: `“${esc(text)}” was read as ${esc(r.via)} (${esc(r.matched || text)}). Read as a SMILES string it would be <code>${esc(a)}</code>, a different structure.` });
        }
      }
      return { mol: m, info };
    }
    // the app's resolver failed: RDKit may still read it as SMILES
    const m = rdMol(text);
    if (m) { info.via = 'SMILES (RDKit)'; return { mol: m, info }; }
    // text with SMILES-only characters: RDKit's reason says more than the name resolver's
    if (looksLikeSmiles(text) && /[\[\]()=#@\/\%]|[A-Za-z]\d/.test(text)) throw new Error(smilesRefusal(text));
    // two identifiers in one box ("103-30-0 trans-stilbene"): the box reads one at a time
    if (/\s/.test(text) && text.split(/\s+/).some(w => /^\d{2,7}-\d\d-\d$/.test(w)))
      throw new Error(`“${text}” holds more than one identifier. Type one at a time: a CAS number, a name or a SMILES.`);
    // a program error inside the app's name reader (a TypeError, not a refusal): say that, not the raw JS message
    if (err && !(err instanceof TypeError || err instanceof RangeError || err instanceof ReferenceError)) throw new Error(err.message);
    if (err) throw new Error(`The app's name reader could not handle “${text}” (it stopped with an internal error). Try the SMILES instead.`);
    throw new Error(`Could not read “${text}”.`);
  }

  /* ================================================================ the molecule and its components */
  function freeMols() {
    for (const cs of S.comps) if (cs && cs.b3d && cs.b3d.job) { try { cs.b3d.job.cancel(); } catch (e) {} }
    for (const m of S.mols) { try { m.delete(); } catch (e) {} }
    S.mols = []; S.full = null; S.comps = []; S.cs = null; S.combined = null;
  }

  async function showMolecule(input, opts) {
    opts = opts || {};
    init();
    const req = ++S.req;
    setStatus('Starting RDKit…', true);
    try { S.R = await window.RDKitLoad(); } catch (e) {
      if (req === S.req) { restoreStatus(); showError('RDKit could not start: ' + e.message); }
      return { ok: false, error: 'RDKit could not start: ' + e.message };
    }
    if (req !== S.req) return { ok: false, error: 'superseded by a newer request' };
    setStatus('Reading the input…', true);           // name resolution and large molecules take a moment
    await nextPaint();
    if (req !== S.req) return { ok: false, error: 'superseded by a newer request' };
    let got;
    try { got = readInput(input, opts.source); } catch (e) {
      // a name or CAS number the curated library does not know: the saved PubChem lookups, then the broad offline
      // library (js/pclib.js), before giving up
      const text = typeof input === 'string' ? input.trim() : '';
      if (lookupable(text) && window.PCLib && !opts.noFallback) {
        const r = await window.PCLib.resolve(text);
        if (req !== S.req) return { ok: false, error: 'superseded by a newer request' };
        if (r.ok) { restoreStatus(); return loadResolved(r.smiles, r.note); }
        if (r.choices) { restoreStatus(); showChoices(r, text); return { ok: false, error: r.error }; }
      }
      // the molecule on screen stays, with its 3D build and its status line: only this request failed
      restoreStatus();
      showError(e.message, text || null);
      return { ok: false, error: e.message };
    }
    // components (salts, mixtures): each is analysed on its own
    let frags = null, fragAtoms = null;
    const typed = typedSmiles(input, got.mol);
    try {
      const f = got.mol.get_frags();
      try { fragAtoms = JSON.parse(f.mappings).fragsMolAtomMapping; } catch (e) { fragAtoms = null; }
      const ml = f.molList, n = ml.size();
      if (n > 1) { frags = []; for (let i = 0; i < n; i++) frags.push(ml.at(i)); }
      else if (n === 1) { const c = ml.at(0); c.delete(); }
      ml.delete();
    } catch (e) { frags = null; }
    const mols = frags || [got.mol];
    const heavy = mols.map(countHeavy);
    let best = 0;
    heavy.forEach((h, i) => { if (h > heavy[best]) best = i; });
    // a large molecule takes seconds to lay out and analyse (one synchronous step): say so first, and give way to a
    // newer request made meanwhile before anything of the molecule shown is replaced (so a later failure of that
    // request still leaves the old molecule whole)
    if (heavy[best] > 150) {
      setStatus('Laying out and analysing a large molecule…', true);
      await nextPaint();
      if (req !== S.req) {
        for (const m of (frags || []).concat(got.mol)) { try { m.delete(); } catch (e) {} }
        return { ok: false, error: 'superseded by a newer request' };
      }
    }
    ++S.gen;                                       // from here on the molecule shown changes
    freeMols();
    hideCard();
    S.measureN = 0;                                // measuring belongs to one molecule
    S.info = got.info;
    S.fileBlock = got.block || null;
    S.full = got.mol; S.mols.push(got.mol);
    if (frags) S.mols.push(...frags);
    $('molSide').scrollTop = 0;                    // a new molecule starts at its name, not mid-list
    S.typed = typed && (!frags || fragAtoms) ? typed : null;
    S.comps = mols.map((m, i) => ({ i, mol: m, heavy: heavy[i], record: null, b3d: null, typedAtoms: frags ? (fragAtoms ? fragAtoms[i] : null) : null }));
    if (frags) { try { S.combined = MI.combined(S.R, S.full); } catch (e) { S.combined = null; } }
    if (opts.view === '2d' || opts.view === '3d') S.view = opts.view;
    const ok = activate(best);
    fillComponentSelect();
    setStatus('');
    renderResolved();
    if (!ok) {
      const error = S.cs && S.cs.error || 'could not analyse the molecule';
      showError(error);
      return { ok: false, error };
    }
    if (S.view === '3d') {
      const B = S.cs.b3d;
      if (!v3) setStatus('Loading the 3D viewer (three.js)…', true);
      else setStatus(B ? B.statusText || '' : '', !!B && !B.done);
    }
    return { ok: true };
  }
  /* The typed SMILES and where each of its atoms is written, when the molecule was read straight from that text:
     RDKit numbers the atoms in the order they are typed, so the k-th atom token is atom k. Anything else (a name, a
     drawing, a file, or a SMILES whose explicit [H] atoms RDKit removed) gives null, and the canonical SMILES is used. */
  const SMILES_ATOM = /\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnops*]/g;
  function typedSmiles(input, mol) {
    const text = typeof input === 'string' ? input.trim() : '';
    if (!text || !looksLikeSmiles(text)) return null;
    const at = [];
    text.replace(SMILES_ATOM, (tok, off) => { at.push([off, off + tok.length]); return tok; });
    if (at.length !== mol.get_num_atoms()) return null;
    let m = null;
    try { m = S.R.get_mol(text); return m && m.get_smiles() === mol.get_smiles() ? { text, at } : null; }
    catch (e) { return null; } finally { if (m) m.delete(); }
  }
  function countHeavy(m) { try { return JSON.parse(m.get_descriptors()).NumHeavyAtoms; } catch (e) { return 0; } }

  function prepareComponent(cs) {
    const M = cs.mol;
    try {
      const hc = M.has_coords ? M.has_coords() : 0;
      cs.from3D = hc === 3;
      // analysed BEFORE any layout made here: RDKit's InChI reads double-bond geometry from 2D coordinates, so a
      // layout (which has to draw an open C=C one way or the other) would decide it — checked: CC=CC got
      // (E)-but-2-ene's InChIKey after set_new_coords(). Coordinates the input brought (a drawing, a file) stay.
      cs.record = MI.analyse(S.R, M);
      cs.pka = window.MoleculePka ? window.MoleculePka.forMol(S.R, M) : null;   // js/reactions.js: measured, estimated or typical
      cs.effects = window.MoleculeEffects ? window.MoleculeEffects.forMol(S.R, M, cs.record) : null;   // steric / electronic per fragment
      if (S.typed) {   // "Location in SMILES" shows the SMILES as typed: this component's atom j is typed atom typedAtoms[j]
        const at = cs.record.atoms.slice(0, cs.record.n).map((x, j) => S.typed.at[cs.typedAtoms ? cs.typedAtoms[j] : j] || null);
        if (at.every(Boolean)) cs.record.smilesLoc = { text: S.typed.text, at };
      }
      if (cs.from3D && S.fileBlock && S.comps.length === 1) cs.fileConf = fileConformer(cs, S.fileBlock);
      if (hc !== 2) layout2D(cs);                    // no coordinates, or a 3D file: RDKit's 2D layout for the drawing
      // RDKit's bond length (1.5): names, drawings and CoordGen come at 1.0, where add_hs_in_place piles H onto
      // ring centres and other atoms (checked: cyclohexanol 6 overlapping pairs at 1.0, none at 1.5). Scale only.
      try { M.normalize_depiction(0, -1); } catch (e) {}
      cs.mb2d = M.get_molblock();
      cs.draw = drawModel(cs, false);                // checked against the record; throws on a mismatch
    } catch (e) {
      cs.error = 'The molecule could not be analysed: ' + (e && e.message || e);
      cs.record = null;
      return false;
    }
    try { cs.pub = window.PubChem3D ? window.PubChem3D.eligible(S.R, M) : { ok: false, reason: 'PubChem3D is not loaded.' }; }
    catch (e) { cs.pub = { ok: false, reason: 'PubChem eligibility could not be checked.' }; }
    const r = cs.record;
    cs.formula = r.formulaText || r.formula || '';
    return true;
  }

  /* RDKit's default 2D layout; CoordGen for macrocycles (a ring of 8 or more atoms) or when the default layout
     stretches a bond past 1.3 × the median (vancomycin: 1.71 × by default, 1.14 × with CoordGen). CoordGen's layout
     is kept only if it reads back with the same stereochemistry. */
  function layout2D(cs) { layoutMol(S.R, cs.mol, cs.record.rings.some(ring => ring.length >= 8)); }
  function layoutMol(R, M, macrocycle) {
    M.set_new_coords();
    if (!macrocycle && bondSpread(M.get_molblock()) <= 1.3) return;
    const smi = M.get_smiles();
    let same = false;
    try {
      M.set_new_coords(true);
      const back = R.get_mol(M.get_molblock(), JSON.stringify({ removeHs: true }));
      same = !!back && back.get_smiles() === smi;
      if (back) back.delete();
    } catch (e) { same = false; }
    if (!same) M.set_new_coords();
  }
  /* This tab's 2D drawing of a SMILES, for other pages (the Reactions page) to draw a molecule the same way: the same
     layout steps and the same final normalisation as prepareComponent. RDKit numbers the atoms in the order the SMILES
     writes them. → {atoms: [{x, y, el}], bonds: [{a, b, order, flag}]} (flag 1 = wedge, 6 = hash, narrow end at a) or null */
  function depict2D(R, smiles) {
    let m = null;
    try { m = R.get_mol(smiles); } catch (e) { m = null; }
    if (!m) return null;
    try {
      if ((m.is_valid && !m.is_valid()) || !m.get_num_atoms()) return null;
      const J = JSON.parse(m.get_json()), ext = (J.molecules[0].extensions || []).find(e => e.name === 'rdkitRepresentation') || {};
      layoutMol(R, m, (ext.atomRings || []).some(ring => ring.length >= 8));
      try { m.normalize_depiction(0, -1); } catch (e) {}
      return parseMolblock(m.get_molblock());
    } catch (e) { return null; } finally { m.delete(); }
  }
  function bondSpread(mb) {                          // longest bond / median bond of a 2D molblock
    const P = parseMolblock(mb), L = P.bonds.map(b => Math.hypot(P.atoms[b.a].x - P.atoms[b.b].x, P.atoms[b.a].y - P.atoms[b.b].y)).sort((a, b) => a - b);
    if (L.length < 3) return 1;
    const med = L[L.length >> 1];
    return med > 0 ? L[L.length - 1] / med : 1;
  }

  /* A 3D MOL/SDF file's own coordinates as a conformer, in record order (UX: a user who opens a crystal or computed
     structure measures that structure, not a rebuilt one). Only when the file has every H: its heavy atoms are
     found by their coordinates (RDKit copied them into M, whose H it removed), each record H takes an H bonded to
     its parent in the file, and the file's bond set must equal the record's — the same check as the 3D builder's
     topology (checkTopology). → {coords} or {error} */
  function fileConformer(cs, block) {
    const r = cs.record;
    let F, P;
    try { F = parseMolblock(block); P = parseMolblock(cs.mol.get_molblock()); } catch (e) { return { error: 'the file could not be read for its coordinates' }; }
    if (r.N > MAX_3D_ATOMS) return { error: `more than ${MAX_3D_ATOMS} atoms` };
    if (P.atoms.length !== r.n) return { error: 'its atoms could not be matched' };
    const isH = a => a.el === 'H' || a.el === 'D' || a.el === 'T';
    const sameEl = (a, b) => a.el === b.el || (isH(a) && isH(b));
    const used = new Int8Array(F.atoms.length), toFile = new Int32Array(r.N).fill(-1);
    for (let i = 0; i < r.n; i++) {
      const p = P.atoms[i];
      let hit = -1;
      for (let j = 0; j < F.atoms.length; j++) {
        const f = F.atoms[j];
        if (used[j] || !sameEl(f, p) || Math.abs(f.x - p.x) > 2e-3 || Math.abs(f.y - p.y) > 2e-3 || Math.abs(f.z - p.z) > 2e-3) continue;
        if (hit >= 0) return { error: 'two of its atoms share a position' };
        hit = j;
      }
      if (hit < 0) return { error: 'its atoms could not be matched' };
      toFile[i] = hit; used[hit] = 1;
    }
    const nbr = F.atoms.map(() => []);
    for (const b of F.bonds) { nbr[b.a].push(b.b); nbr[b.b].push(b.a); }
    for (let i = r.n; i < r.N; i++) {
      const fp = toFile[r.atoms[i].parent];
      const h = nbr[fp].filter(j => !used[j] && F.atoms[j].el === 'H').sort((a, b) => a - b)[0];
      if (h === undefined) return { error: 'the file does not give every hydrogen' };
      toFile[i] = h; used[h] = 1;
    }
    if (used.some(u => !u)) return { error: 'the file has atoms the structure read from it does not' };
    const fromFile = new Int32Array(F.atoms.length);
    toFile.forEach((j, i) => { fromFile[j] = i; });
    const res = {
      N: r.N, elements: Array.from(toFile, (j, i) => (i < r.n ? r.atoms[i].el : F.atoms[j].el)),
      parents: Array.from(toFile, (j, i) => (i < r.n ? r.atoms[i].parent : fromFile[nbr[j][0]])),
      bonds: F.bonds.map(b => [fromFile[b.a], fromFile[b.b]]),
    };
    for (let i = r.n; i < r.N; i++) if (nbr[toFile[i]].length !== 1) return { error: 'a hydrogen in the file has more than one bond' };
    const why = checkTopology(res, r);
    if (why) return { error: 'internal atom-mapping check failed (' + why + ')' };
    const coords = new Float64Array(3 * r.N);
    toFile.forEach((j, i) => { const a = F.atoms[j]; coords[3 * i] = a.x; coords[3 * i + 1] = a.y; coords[3 * i + 2] = a.z; });
    return { coords };
  }

  function activate(ci) {
    const cs = S.comps[ci];
    if (!cs) return false;
    if (S.cs && S.cs !== cs && S.cs.b3d && S.cs.b3d.job && !S.cs.b3d.done) {
      try { S.cs.b3d.job.cancel(); } catch (e) {}  // switching component: stop the other build (its conformers are kept)
    }
    S.cs = cs;
    S.sel = null; S.group = null; S.picks = []; S.measures = [];
    hideCard();
    if (!cs.record && !cs.error) prepareComponent(cs);
    if (cs.error) {                                 // nothing of the previous molecule may stay on screen
      if (v2) { v2.setMolecule(null); v2info = null; }
      if (v3) { v3cs = null; shownCoords = null; $('mol3d').hidden = true; }
      setStatus('');
      $('molLegend').hidden = true;
      $('molEmpty').innerHTML = `<div><h2>No molecule shown</h2><p class="mv-err">${esc(cs.error)}</p></div>`;
      $('molEmpty').hidden = false;
      renderAll();
      showError(cs.error);
      return false;
    }
    $('molEmpty').hidden = true;
    render2D();
    renderLabels();
    const gen = S.gen;
    if (!cs.names && cs.heavy > NAMER_MAX_HEAVY) {  // the namer takes seconds to minutes on large molecules
      cs.names = { iupac: null, iupacNote: `not attempted: with more than ${NAMER_MAX_HEAVY} heavy atoms the molecule is too large for the app's namer`, common: null, cas: null };
    }
    if (!cs.names) {
      setTimeout(() => {                            // the namer can take ~0.2 s: after the first paint
        if (gen !== S.gen || cs !== S.cs || cs.names) return;
        try { cs.names = MI.names(S.R, cs.mol); } catch (e) { cs.names = { iupac: null, iupacNote: "the app's namer failed on this structure", common: null, cas: null }; }
        renderHead(); renderSections();
      }, 30);
    }
    if (S.view === '3d') enter3D();
    applyView();
    return true;
  }

  function fillComponentSelect() {
    const wrap = $('molCompWrap'), sel = $('molComp');
    if (S.comps.length < 2) { wrap.hidden = true; sel.innerHTML = ''; return; }
    const comb = S.combined && S.combined.components;
    sel.innerHTML = S.comps.map((cs, i) => {
      let f = cs.formula;
      if (!f && comb) { const c = comb[i]; f = c ? c.formulaText : ''; }
      if (!f) { try { f = cs.mol.get_smiles(); } catch (e) { f = '?'; } }
      return `<option value="${i}">${i + 1}: ${esc(f)} (${cs.heavy} heavy atom${cs.heavy === 1 ? '' : 's'})</option>`;
    }).join('');
    sel.value = String(S.cs ? S.cs.i : 0);
    wrap.hidden = false;
  }

  /* ================================================================ 2D */
  /* atoms {x, y, z, el} and bonds {a, b, order, flag (V2000 stereo code)} of a molblock. V3000 too: RDKit writes it
     above 999 atoms, and reading only the V2000 counts line made such molecules fail the mapping check. */
  function parseMolblock(mb) {
    const L = String(mb).split(/\r\n|\r|\n/);
    if (/V3000\s*$/.test(L[3] || '')) return parseV3000(L);
    const na = parseInt(L[3].slice(0, 3), 10), nb = parseInt(L[3].slice(3, 6), 10);
    const atoms = [], bonds = [];
    for (let i = 0; i < na; i++) { const s = L[4 + i]; atoms.push({ x: parseFloat(s.slice(0, 10)), y: parseFloat(s.slice(10, 20)), z: parseFloat(s.slice(20, 30)) || 0, el: s.slice(31, 34).trim() }); }
    for (let k = 0; k < nb; k++) {
      const s = L[4 + na + k];
      bonds.push({ a: parseInt(s.slice(0, 3), 10) - 1, b: parseInt(s.slice(3, 6), 10) - 1, order: parseInt(s.slice(6, 9), 10), flag: parseInt(s.slice(9, 12), 10) || 0 });
    }
    return { atoms, bonds };
  }
  function parseV3000(L) {
    const v = [];
    let buf = '';
    for (let i = 4; i < L.length; i++) {
      const l = L[i];
      if (/^M {2}END/.test(l)) break;
      if (!l.startsWith('M  V30 ')) continue;
      const t = buf + l.slice(7);
      if (/-\s*$/.test(t)) { buf = t.replace(/-\s*$/, ''); continue; }   // continuation line
      buf = ''; v.push(t.trim());
    }
    const atoms = [], bonds = [];
    let blk = '';
    for (const l of v) {
      if (/^BEGIN (ATOM|BOND)/.test(l)) { blk = l.slice(6); continue; }
      if (/^END /.test(l)) { blk = ''; continue; }
      const f = l.split(/\s+/);
      if (blk === 'ATOM') atoms.push({ x: parseFloat(f[2]), y: parseFloat(f[3]), z: parseFloat(f[4]) || 0, el: f[1] });
      else if (blk === 'BOND') {
        const order = +f[1], cfg = +((/CFG=(\d)/.exec(l) || [])[1] || 0);
        // V3000 CFG → the V2000 codes used below: single 1 wedge (1), 2 either (4), 3 hash (6); double 2 either (3)
        const flag = order === 1 ? ({ 1: 1, 2: 4, 3: 6 }[cfg] || 0) : order === 2 && cfg === 2 ? 3 : 0;
        bonds.push({ a: +f[2] - 1, b: +f[3] - 1, order, flag });
      }
    }
    return { atoms, bonds };
  }
  const STEREO_FLAG = { 1: 'wedge', 6: 'hash', 4: 'wavy' };
  const pairKey = (a, b) => (a < b ? a + '-' + b : b + '-' + a);
  const mapFail = what => { throw new Error('internal atom-mapping check failed (' + what + '), so it is not shown'); };

  /* Mol2DView input in record numbering. Without H: M's molblock (atoms 0..n-1 = record atoms). With H: RDKit's
     add_hs_in_place places the H; they are matched to the record's H by parent (H on one atom are interchangeable). */
  function drawModel(cs, withH) {
    const r = cs.record, n = r.n;
    let mb = cs.mb2d;
    if (withH) {
      const c = cs.mol.copy();
      try {
        c.add_hs_in_place(); mb = c.get_molblock();
        // safety net: an H on top of another atom → lay the H copy out again, H included (CoordGen first: crowded
        // macrocycles such as vancomycin keep 4-6 such H with the other layouts, none with CoordGen), and keep the
        // least crowded layout, the first one on a tie (paclitaxel: 6 such pairs as placed, 2 with CoordGen)
        let crowd = crowdedH(parseMolblock(mb));
        for (const cg of [true, false]) {
          if (!crowd) break;
          try {
            c.set_new_coords(cg); try { c.normalize_depiction(0, -1); } catch (e) {}
            const t = c.get_molblock(), n = crowdedH(parseMolblock(t));
            if (n < crowd) { mb = t; crowd = n; }
          } catch (e) {}
        }
      } finally { c.delete(); }
    }
    const P = parseMolblock(mb);
    const nA = withH ? r.N : n;
    if (P.atoms.length !== nA) mapFail(`2D: ${P.atoms.length} atoms, the record has ${nA}`);
    // a SMILES wildcard * (atomic number 0) comes back from RDKit's molfile writer as R (or A, Q, R#): the same atom
    const sameEl = (mb, rec) => mb === rec || (rec === '*' && /^(R#?|A|Q|\*)$/.test(mb));
    for (let i = 0; i < n; i++) if (!sameEl(P.atoms[i].el, r.atoms[i].el)) mapFail(`2D: atom ${i + 1} is ${P.atoms[i].el}, the record has ${r.atoms[i].el}`);
    for (let i = 0; i < n; i++) if (r.atoms[i].el === '*') P.atoms[i].el = '*';   // draw it as * like the input
    const toRec = new Int32Array(P.atoms.length).fill(-1);
    for (let i = 0; i < n; i++) toRec[i] = i;
    if (withH) {
      const hOfParent = new Map();                 // record: parent → its appended H, in order
      for (let i = n; i < r.N; i++) { const p = r.atoms[i].parent; if (!hOfParent.has(p)) hOfParent.set(p, []); hOfParent.get(p).push(i); }
      const used = new Map();
      for (const b of P.bonds) {
        for (const [h, p] of [[b.a, b.b], [b.b, b.a]]) {
          if (h < n || p >= n) continue;
          if (P.atoms[h].el !== 'H') mapFail('2D: an added atom is not H');
          const list = hOfParent.get(p) || [], t = used.get(p) || 0;
          if (t >= list.length) mapFail(`2D: too many H on atom ${p + 1}`);
          toRec[h] = list[t]; used.set(p, t + 1);
        }
      }
      if ([...toRec].some(v => v < 0)) mapFail('2D: an H has no parent');
    }
    // record bonds 0..m0-1 are M's bonds (both ends < n); the X–H bonds of the appended H follow
    const recBond = new Map(r.bonds.map(b => [pairKey(b.a, b.b), b.index]));
    const nb = withH ? r.bonds.length : r.bonds.filter(b => b.a < n && b.b < n).length;
    if (P.bonds.length !== nb) mapFail(`2D: ${P.bonds.length} bonds, the record has ${nb}`);
    const atoms = new Array(nA), bonds = new Array(nb);
    P.atoms.forEach((p, j) => {
      const i = toRec[j], a = r.atoms[i];
      const hCount = withH ? 0 : r.atoms.reduce((s, x) => s + (x.parent === i ? 1 : 0), 0);
      atoms[i] = { el: a.el, x: p.x, y: p.y, charge: a.charge, isotope: a.isotope, hCount, radical: a.radical };
    });
    for (const b of P.bonds) {
      const ra = toRec[b.a], rb = toRec[b.b], k = recBond.get(pairKey(ra, rb));
      if (k === undefined || k >= nb) mapFail(`2D: bond ${b.a + 1}–${b.b + 1} is not in the record`);
      const rb0 = r.bonds[k];
      let stereo = null;
      if (b.order === 2 && b.flag === 3) stereo = 'either';
      else if (b.order === 1 && STEREO_FLAG[b.flag]) stereo = STEREO_FLAG[b.flag];
      bonds[k] = { a: rb0.a, b: rb0.b, order: rb0.order, aromatic: rb0.aromatic, stereo, narrow: stereo === 'wedge' || stereo === 'hash' ? ra : undefined };
    }
    for (let k = 0; k < nb; k++) if (!bonds[k]) mapFail('2D: a bond is missing');
    return { atoms, bonds, rings: r.rings, N: nA, M: nb, withH };
  }

  // how crowded the H of a 2D layout are: the number of atom pairs, at least one of them H, closer than 0.3 × the
  // mean bond length (0 = no H drawn on top of another atom)
  function crowdedH(P) {
    if (!P.bonds.length) return 0;
    const A = P.atoms, d = (i, j) => Math.hypot(A[i].x - A[j].x, A[i].y - A[j].y);
    const cut = 0.3 * P.bonds.reduce((s, b) => s + d(b.a, b.b), 0) / P.bonds.length;
    let n = 0;
    for (let i = 0; i < A.length; i++) {
      if (A[i].el !== 'H') continue;
      for (let j = 0; j < A.length; j++) if (j !== i && (A[j].el !== 'H' || j > i) && d(i, j) < cut) n++;
    }
    return n;
  }

  function render2D() {
    const cs = S.cs;
    if (!cs || !cs.record) return;
    if (!v2) v2 = window.Mol2DView.create($('mol2d'), { maxBondPx: 64, textbook: S.textbook2d, onHover: (h, x, y) => onHover(h, x, y, '2d'), onClick: (h, e) => onClick(h, e, '2d') });
    let model = cs.draw;
    if (S.showH2d) {
      try { cs.drawH = cs.drawH || drawModel(cs, true); model = cs.drawH; }
      catch (e) { cs.drawHError = e.message; model = cs.draw; }
    }
    v2.setMolecule({ atoms: model.atoms, bonds: model.bonds, rings: model.rings });
    v2info = { cs, withH: model.withH, N: model.N, M: model.M };
    applyLabels2D();
    applyHighlights();
  }

  /* ================================================================ 3D: loading, building, conformers */
  function ensureThree() {
    if (window.THREE && window.THREE.WebGLRenderer) return Promise.resolve();
    if (!threeLoading) {
      threeLoading = window.RDKitLoad.script(THREE_SRC).then(() => { if (!window.THREE) throw new Error('three.js did not load'); });
      threeLoading.catch(() => { threeLoading = null; });
    }
    return threeLoading;
  }
  function topology3D(r) {
    return {
      atoms: r.atoms.map(a => ({ el: a.el === 'H' && a.isotope === 2 ? 'D' : a.el === 'H' && a.isotope === 3 ? 'T' : a.el, z: a.z, rvdw: a.rvdw, rcov: a.rcov, isH: a.appended })),
      bonds: r.bonds.map(b => ({ a: b.a, b: b.b, order: b.order, aromatic: b.aromatic })),
      rings: r.rings,
    };
  }
  function create3D() {
    if (v3) return;
    const box = $('mol3d');
    box.hidden = false;
    v3 = window.Mol3DView.create(box, {
      onHover: (h, x, y) => onHover(h, x, y, '3d'), onClick: (h, e) => onClick(h, e, '3d'),
      onSpinChange: on => { S.spin = on; syncControls(); },
      onSurface: (res, err) => onSurfaceResult(res, err),
      onContextLost: () => { contextLost = true; renderBanners(); },
      onContextRestored: () => { contextLost = false; renderBanners(); },
    });
    v3.setStyle(S.style);
    v3.setSpin(S.spin);
    applyPickMode();
  }
  // Measure picks atoms only: a click on a bond (with Show H, often a C–H bond in front of a methyl carbon) or just
  // beside an atom then picks the nearest atom instead of being ignored
  function applyPickMode() { if (v3 && v3.setPickMode) v3.setPickMode(S.measureN ? 'atoms' : 'all'); }
  // put the active component into the 3D view (topology, labels, highlights, surface) if it is not there already
  function bind3D(cs) {
    if (v3cs === cs) return;
    v3.setMolecule(topology3D(cs.record));
    v3cs = cs; shownCoords = null;
    v3.setShowH(S.showH3d);
    surfaceResult = null;
    renderLegend();
    const c = currentConf(cs);
    // the values (cs.geo) always come from the coordinate array on screen
    if (c) { v3.setCoords(c.coords, { keepView: false }); shownCoords = c.coords; cs.geo = geometryOf(cs, c); }
    applySurface();
    applyLabels3D();
    applyMeasures();
    applyHighlights();
  }

  async function enter3D() {
    const cs = S.cs, gen = S.gen;
    if (!cs || !cs.record) return;
    ensureB3D(cs);                                  // the file's own conformer, if any, is there before the first draw
    if (!v3) {
      setStatus('Loading the 3D viewer (three.js)…', true);
      try { await ensureThree(); } catch (e) {
        if (gen === S.gen) { setStatus(''); addBanner3D(cs, 'error', 'The 3D viewer could not load: ' + esc(e.message)); }
        return;
      }
      if (gen !== S.gen || S.view !== '3d' || cs !== S.cs) { setStatus(''); return; }
      try { create3D(); } catch (e) { setStatus(''); addBanner3D(cs, 'error', 'The 3D viewer could not start (WebGL): ' + esc(e.message)); return; }
    }
    if (S.view !== '3d' || cs !== S.cs) return;
    bind3D(cs);
    // the conformer picked while the view showed something else (PubChem answering in 2D): draw it now, so the
    // model on screen and every 3D value come from the same coordinates
    const cur = currentConf(cs);
    if (cur && shownCoords !== cur.coords) displayConformer(cs, false);
    const B = cs.b3d;
    // no build yet, or one cancelled before its first conformer (e.g. by switching component): build now.
    // B.done is false from startBuild's first line, so a 2D↔3D toggle while OpenChemLib loads starts no second build.
    if (B.done && !B.job && !B.error && !B.confs.some(c => c.source === 'OpenChemLib')) startBuild(cs, null);
    else setStatus(B.statusText || '', !B.done);
    renderAll();
  }

  // done: true = no build running (a record made by a banner, PubChem or the file conformer must not look busy,
  // or no build would ever start for it); startBuild sets it false
  function newB3D() { return { job: null, result: null, confs: [], cur: null, userPicked: false, done: true, error: null, banners: [], pub: null, pubBusy: false, checked: false, statusText: '' }; }
  function ensureB3D(cs) {
    const B = cs.b3d || (cs.b3d = newB3D());
    const F = cs.fileConf;
    if (F && F.coords && !B.confs.some(c => c.source === 'file')) {
      B.confs.unshift({ key: 'file', source: 'file', coords: F.coords, energy: NaN, relEnergy: NaN, optimised: true, flaggedBonds: [], warnings: [] });
      if (!B.cur) { B.cur = 'file'; B.userPicked = true; }   // the file's structure is what was opened: shown first
    }
    return B;
  }

  // engine topology == record: elements, H parents, charges and the bond set (the record numbering above); '' or the reason
  function checkTopology(res, r) {
    if (res.N !== r.N) return `the 3D model has ${res.N} atoms, the record ${r.N}`;
    for (let i = 0; i < r.N; i++) {
      const a = r.atoms[i], e = res.elements[i];
      const el = e === 'D' || e === 'T' ? 'H' : e;
      if (el !== a.el) return `atom ${i + 1} is ${e} in the 3D model, ${a.el} in the record`;
      if (res.parents && (res.parents[i] | 0) !== a.parent) return `H ${i + 1} sits on a different atom`;
      if (res.charges && (res.charges[i] | 0) !== a.charge) return `atom ${i + 1} has a different charge`;
      if (res.isotopes && a.isotope && (res.isotopes[i] | 0) !== a.isotope) return `atom ${i + 1} has a different isotope`;
    }
    if (!res.bonds || res.bonds.length !== r.bonds.length) return 'the bond count differs';
    const set = new Set(r.bonds.map(b => pairKey(b.a, b.b)));
    for (const b of res.bonds) if (!set.has(pairKey(b[0], b[1]))) return `bond ${b[0] + 1}–${b[1] + 1} is not in the record`;
    return '';
  }

  function startBuild(cs, previous) {
    const B = ensureB3D(cs);
    const wild = cs.record.atoms.filter(a => a.el === '*').length;
    if (cs.record.N > MAX_3D_ATOMS || wild) {      // refused up front, with the reason (not a failed mapping check)
      B.done = true; B.job = null; B.statusText = '';
      B.error = wild
        ? `the structure contains a wildcard atom (*), which stands for "any atom" and has no size, bonds or force-field type, so no real 3D shape exists. Replace * with an actual group (e.g. C for methyl) to get a 3D model; the 2D drawing and the values from the structure are unaffected`
        : `this molecule has ${cs.record.N} atoms with its hydrogens; 3D models are built for at most ${MAX_3D_ATOMS} atoms (the limit of the V2000 molfile the 3D builder reads and the 3D exports write). The 2D drawing and the values from the structure are unaffected`;
      if (cs === S.cs) { restoreStatus(); renderAll(); }
      return;
    }
    B.done = false; B.error = null;
    B.statusText = previous ? 'Searching for more conformers…' : 'Loading OpenChemLib…';
    if (cs === S.cs && S.view === '3d') setStatus(B.statusText, true);
    renderControls();
    const E = window.Mol3DEngine;
    // results count while this component belongs to the molecule shown (a new molecule makes new components; a
    // failed Show leaves them, and their build, in place)
    const live = () => S.comps.includes(cs);
    const go = () => {
      if (!live()) return;
      if (cs !== S.cs) { B.done = true; B.statusText = ''; return; }   // switched away while loading: build on return
      if (!previous) B.statusText = 'Building the 3D model (OpenChemLib conformer generator, then MMFF94)…';
      if (cs === S.cs && S.view === '3d') setStatus(B.statusText, true);
      const job = E.build(cs.mb2d, { previous: previous || undefined, onProgress: r => { if (live() && B.job === job) takeResult(cs, r, false); } });
      B.job = job;
      job.promise.then(r => { if (live() && B.job === job) takeResult(cs, r, true); });
    };
    E.ready().then(mode => {
      if (mode === 'main' && !S.mainNoted) {
        S.mainNoted = true;
        const why = (E.status() || {}).workerError;
        addBanner3D(cs, 'info', 'The 3D engine runs on the page itself (no background worker' + (why ? ': ' + esc(why) : '') + '), so the page may pause while a model is built.');
      }
      go();
    }, e => {
      if (!live()) return;
      B.done = true; B.error = 'The 3D engine (OpenChemLib) could not load: ' + (e && e.message || e);
      B.statusText = '';
      if (cs === S.cs) { setStatus(''); renderAll(); }
    });
  }

  function takeResult(cs, r, final) {
    const B = cs.b3d;
    if (!r) return;
    if (final) { B.done = true; B.job = null; }
    if (!r.ok) {
      if (final && !r.cancelled) {
        if (!B.confs.some(c => c.source === 'OpenChemLib')) B.error = r.error || 'the 3D model could not be built';
        else B.error = null;
        B.statusText = '';
      }
      if (final && r.cancelled) B.statusText = '';
      if (cs === S.cs) { if (S.view === '3d') setStatus(''); renderAll(); }
      return;
    }
    if (!B.checked) {
      const why = checkTopology(r, cs.record);
      if (why) {
        B.error = 'internal atom-mapping check failed (' + why + '), so the 3D model is not shown';
        if (B.job) { try { B.job.cancel(); } catch (e) {} }
        B.job = null; B.done = true; B.statusText = '';
        if (cs === S.cs) { setStatus(''); renderAll(); }
        return;
      }
      B.checked = true;
    }
    B.result = r;
    const n = r.conformers.length;
    // keep the same conformer object (and coordinate array) for a conformer already known, so the view is not
    // rebuilt at every progress report and its read-back stereo stays cached
    const old = new Map(B.confs.map(c => [c.key, c]));
    B.confs = B.confs.filter(c => c.source === 'file').concat(r.conformers.map((c, t) => {
      const was = old.get('ocl:' + c.id);
      return Object.assign(was || { stereoRead: null }, {
        key: 'ocl:' + c.id, id: c.id, source: 'OpenChemLib', coords: was ? was.coords : c.coords, energy: c.energy, relEnergy: c.relEnergy,
        optimised: c.optimised, flaggedBonds: c.flaggedBonds || [], warnings: c.warnings || [], rank: t + 1, of: n, note: c.note || null,
      });
    }), B.confs.filter(c => c.source === 'PubChem'));
    B.engineWarnings = r.warnings || [];
    B.search = r.search || null;
    // the conformer shown may have been merged into a duplicate with a lower energy: follow it
    if (B.cur && B.cur.startsWith('ocl:') && !B.confs.some(c => c.key === B.cur)) {
      const id = +B.cur.slice(4), m = (r.merged || []).find(x => x[0] === id);
      B.cur = m ? 'ocl:' + m[1] : null;
    }
    const low = B.confs.find(c => c.source === 'OpenChemLib');
    const first = !B.cur;
    if (first) B.cur = low ? low.key : B.confs.length ? B.confs[0].key : null;
    // when the search ends, show the lowest-energy conformer unless the user picked one (or PubChem's / the file's)
    if (final && !B.userPicked && low && (B.cur || '').startsWith('ocl:')) B.cur = low.key;
    if (final) B.statusText = '';
    else B.statusText = `Searching for more conformers… ${n} found so far`;
    if (cs !== S.cs) return;
    if (S.view === '3d') setStatus(B.statusText, !final);
    displayConformer(cs, first);
    renderAll();
  }

  function currentConf(cs) {
    const B = cs && cs.b3d;
    if (!B || !B.cur) return null;
    return B.confs.find(c => c.key === B.cur) || null;
  }
  /* A lone atom or ion ([Na+], [Cl-]): MMFF94 can type it, but there is no geometry to optimise and its "energy" is
     0 by definition, so no energy is shown and it is not called minimised. The engine says so (result.singleAtom and
     the conformer's note); a component of one atom with no H is the same case if an engine does not. → '' or the note */
  const SINGLE_ATOM = 'single atom: no geometry to optimise';
  function singleAtomNote(cs, c) {
    if (!cs || !c || c.source !== 'OpenChemLib') return '';
    const B = cs.b3d;
    return (B && B.result && B.result.singleAtom) || (cs.record && cs.record.N === 1) ? String(c.note || SINGLE_ATOM) : '';
  }
  const hasEnergy = (cs, c) => !!c && c.optimised && Number.isFinite(c.energy) && !singleAtomNote(cs, c);
  // A lone atom's geometry is a point whatever MMFF94 makes of it ([2H]: no atom type for a free H atom), so it is not
  // "approximate": it gets the single-atom note only, not the not-optimised banner and "—" values
  const geometryOK = (cs, c) => !!c && (c.optimised || !!singleAtomNote(cs, c));
  /* Conformers come out of the generator in unrelated orientations. So that switching conformer keeps the view
     meaningful, a conformer about to be shown is turned onto the one on screen (least-squares fit of the heavy atoms,
     Horn's quaternion method: always a proper rotation, never a mirror image, so every length, angle, dihedral and
     configuration is unchanged). The coordinates are changed in place, so exports and measurements use what is shown. */
  function eigenSym4(A) {                          // cyclic Jacobi; → eigenvector of the largest eigenvalue
    const a = A.map(r => r.slice()), V = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    const scale = A.reduce((s, r) => s + r.reduce((t, x) => t + x * x, 0), 0) || 1;
    for (let sweep = 0; sweep < 80; sweep++) {
      let off = 0;
      for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) off += a[p][q] * a[p][q];
      if (off <= 1e-26 * scale) break;
      for (let p = 0; p < 4; p++) for (let q = p + 1; q < 4; q++) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const th = (a[q][q] - a[p][p]) / (2 * a[p][q]), t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < 4; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = c * x - s * y; a[k][q] = s * x + c * y; }
        for (let k = 0; k < 4; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = c * x - s * y; a[q][k] = s * x + c * y; }
        for (let k = 0; k < 4; k++) { const x = V[k][p], y = V[k][q]; V[k][p] = c * x - s * y; V[k][q] = s * x + c * y; }
      }
    }
    let best = 0;
    for (let j = 1; j < 4; j++) if (a[j][j] > a[best][best]) best = j;
    return [V[0][best], V[1][best], V[2][best], V[3][best]];
  }
  function alignOnto(X, ref, idx) {                // rotate + translate X (in place) onto ref, fitting atoms idx
    const m = idx.length;
    if (!m) return;
    const cx = [0, 0, 0], cr = [0, 0, 0];
    for (const i of idx) for (let d = 0; d < 3; d++) { cx[d] += X[3 * i + d] / m; cr[d] += ref[3 * i + d] / m; }
    const S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];     // S[u][v] = Σ x_u · ref_v (moving × reference)
    for (const i of idx) for (let u = 0; u < 3; u++) for (let v = 0; v < 3; v++) S[u][v] += (X[3 * i + u] - cx[u]) * (ref[3 * i + v] - cr[v]);
    const [[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]] = S;
    const q = m < 2 ? [1, 0, 0, 0] : eigenSym4([
      [xx + yy + zz, yz - zy, zx - xz, xy - yx],
      [yz - zy, xx - yy - zz, xy + yx, zx + xz],
      [zx - xz, xy + yx, -xx + yy - zz, yz + zy],
      [xy - yx, zx + xz, yz + zy, -xx - yy + zz]]);
    const [q0, q1, q2, q3] = q;
    const R = [[q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)],
      [2 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1)],
      [2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3]];
    for (let i = 0; i < X.length / 3; i++) {
      const p = [X[3 * i] - cx[0], X[3 * i + 1] - cx[1], X[3 * i + 2] - cx[2]];
      for (let d = 0; d < 3; d++) X[3 * i + d] = R[d][0] * p[0] + R[d][1] * p[1] + R[d][2] * p[2] + cr[d];
    }
  }
  let shownCoords = null;
  function displayConformer(cs, frame) {
    const c = currentConf(cs);
    if (c && v3 && v3cs === cs && !frame && shownCoords && shownCoords !== c.coords && shownCoords.length === c.coords.length) {
      const r = cs.record, heavy = r.atoms.filter(a => !a.appended).map(a => a.index);
      alignOnto(c.coords, shownCoords, heavy.length >= 3 ? heavy : r.atoms.map(a => a.index));
    }
    cs.geo = c ? geometryOf(cs, c) : null;
    if (!v3 || v3cs !== cs || !c) return;
    if (shownCoords !== c.coords) {
      v3.setCoords(c.coords, { keepView: !frame });
      shownCoords = c.coords;
      applyMeasures();
    }
  }
  function geometryOf(cs, c) {
    try { return MI.geometry(cs.record, c.coords, { optimised: geometryOK(cs, c), flaggedBonds: c.flaggedBonds }); }
    catch (e) { return null; }
  }
  function pickConformer(key) {
    const cs = S.cs;
    if (!cs || !cs.b3d) return;
    cs.b3d.cur = key; cs.b3d.userPicked = true;
    displayConformer(cs, false);
    renderAll();
  }

  /* ---------- PubChem 3D (only on a click) ---------- */
  async function loadPubChem() {
    const cs = S.cs, P = window.PubChem3D;
    if (!cs || !P) return;
    const B = ensureB3D(cs);
    if (!cs.pub || !cs.pub.ok) { addBanner3D(cs, 'info', esc(cs.pub ? cs.pub.reason : 'PubChem 3D is not available.')); return; }
    if (P.busyFor() > 0) { addBanner3D(cs, 'warn', esc(P.TEXT.busy)); renderControls(); return; }
    B.pubBusy = true;
    B.banners = B.banners.filter(b => !b.pubchem);
    renderControls();
    if (S.view === '3d') setStatus('Asking PubChem for its 3D model…', true);
    const Mc = cs.mol.copy();                     // PubChem3D reads it after the network: keep our own copy
    let res;
    try { res = await P.load(S.R, Mc, cs.record); } catch (e) { res = { ok: false, kind: 'error', error: 'PubChem 3D failed: ' + (e && e.message || e) }; }
    finally { try { Mc.delete(); } catch (e) {} }
    B.pubBusy = false;
    if (!S.comps.includes(cs)) return;
    if (cs === S.cs && S.view === '3d') setStatus(B.statusText || '', !B.done && !!B.job);
    if (!res.ok) {
      const mm = res.mismatch;
      let html = esc(res.error);
      if (mm) {
        html = esc(mm.message) + ` <a href="${esc(mm.url)}" target="_blank" rel="noopener">Open on PubChem</a>`;
        if (mm.pubchemSmiles) html += ` <button class="mv-link" data-act="pubform" data-smiles="${esc(mm.pubchemSmiles)}" data-cid="${esc(res.cid)}">Load PubChem's form instead</button>`;
        if (mm.pubchemSmilesNote) html += `<br><small>${esc(mm.pubchemSmilesNote)}</small>`;
      }
      B.banners.push({ kind: res.kind === 'mismatch' || res.kind === 'identity' ? 'warn' : res.kind === 'ineligible' || res.kind === 'not-found' || res.kind === 'no-3d' ? 'info' : 'error', html, pubchem: true });
      if (cs === S.cs) renderAll();
      return;
    }
    // the coordinates are in record order (PubChem3D checked the mapping); check the size once more here
    if (!res.coords || res.coords.length !== 3 * cs.record.N) {
      B.banners.push({ kind: 'error', html: 'Internal atom-mapping check failed (PubChem coordinates do not fit this molecule), so PubChem\'s 3D model is not used.', pubchem: true });
      if (cs === S.cs) renderAll();
      return;
    }
    B.pub = res;
    B.confs = B.confs.filter(c => c.source !== 'PubChem').concat([{
      key: 'pubchem', source: 'PubChem', coords: res.coords, energy: NaN, relEnergy: NaN, optimised: true, flaggedBonds: [],
      warnings: [], cid: res.cid, conformerId: res.conformerId,
    }]);
    // PubChem's own stereo banner: the page's single stereo banner names the configuration read back from the
    // model shown, PubChem's included (kept here for its mixed-stereoisomers note)
    if (res.identity && res.identity.banner) B.banners.push({ kind: 'warn', html: esc(res.identity.banner), pubchem: true, stereo: true });
    B.cur = 'pubchem'; B.userPicked = true;
    // every path: the values follow the conformer now current (even for a component not shown, or in 2D);
    // displayConformer below recomputes them after turning the model onto the one on screen
    cs.geo = geometryOf(cs, currentConf(cs));
    if (cs !== S.cs) return;
    if (S.view === '3d' && v3 && v3cs === cs) displayConformer(cs, false);
    renderAll();
  }

  /* ================================================================ hover card, selection, highlights */
  const card = $('molCard');
  let cardKey = '';
  function hideCard() { card.hidden = true; cardKey = ''; }
  function rowsFor(hit, full) {
    const cs = S.cs, r = cs && cs.record;
    if (!r || !hit) return [];
    const geo = S.view === '3d' ? cs.geo : null;
    let rows = hit.type === 'atom' ? MI.atomRows(r, hit.index, geo) : MI.bondRows(r, hit.index, geo);
    if (!full) rows = rows.filter(x => x.main);
    const c = S.view === '3d' ? currentConf(cs) : null;
    rows = rows.map(x => geometrySource(x, c));
    // the engine flags atoms MMFF94 has no types for (radical, carbocation, carbanion): measured values that involve
    // such an atom — or any atom of its conjugated system, whose bonds MMFF94 then also gets wrong — are marked
    // approximate, with the engine's reason
    const fl = flaggedAtoms(cs, c);
    if (fl.size) {
      const touched = hit.type === 'atom' ? [hit.index] : [r.bonds[hit.index].a, r.bonds[hit.index].b];
      const why = [...new Set(touched.map(i => fl.get(i)).filter(Boolean))];
      if (why.length) rows = rows.map(x => (/^(angles|length|dihedral)$/.test(x.key) && /\d/.test(x.value) && !/flagged/.test(x.value)
        ? Object.assign({}, x, { value: x.value + ' (approximate)', note: why.join(' ') + (x.note ? ' ' + x.note : '') }) : x));
    }
    return rows;
  }
  /* atom → the engine's reason, for the atoms it flags and every atom conjugated with them (a charge or radical on
     one atom of a delocalised system: MMFF94 mistreats the whole system — tropylium's seven C–C all come out at a
     double-bond length, not only the two at the atom the charge is drawn on) */
  function flaggedAtoms(cs, c) {
    const m = new Map(), B = cs && cs.b3d;
    if (!c || c.source !== 'OpenChemLib' || !B || !B.result || !B.result.atomFlags) return m;
    const r = cs.record, adj = r.atoms.map(() => []);
    for (const b of r.bonds) if (b.conjugated || b.aromatic) { adj[b.a].push(b.b); adj[b.b].push(b.a); }
    for (const f of B.result.atomFlags) {
      const q = [f.atom];
      m.set(f.atom, f.text);
      while (q.length) for (const y of adj[q.shift()] || []) if (!m.has(y)) { m.set(y, f.text); q.push(y); }
    }
    return m;
  }
  const PUBCHEM_GEOMETRY = 'PubChem3D computed geometry of this conformer (OMEGA, MMFF94s; not an energy minimum), not an experimental structure';
  const FILE_GEOMETRY = 'geometry from the opened file, used as given (not minimised or otherwise changed by this app)';
  // MolInfo words 3D notes for the app's own models (MMFF94): name the real source for PubChem's and the file's
  function geometrySource(x, c) {
    const t = c && c.source === 'PubChem' ? PUBCHEM_GEOMETRY : c && c.source === 'file' ? FILE_GEOMETRY : null;
    return t && x.note && x.note.includes(MI.TEXT.geometry3d) ? Object.assign({}, x, { note: x.note.split(MI.TEXT.geometry3d).join(t) }) : x;
  }
  function tagHTML(kind) {
    const K = KIND[kind];
    return K && K.tag ? `<span class="mv-tag mv-tag-${kind === '3d' ? 'three' : kind}" title="${esc(K.title)}">${esc(K.tag)}</span>` : '';
  }
  function onHover(hit, x, y, from) {
    if (from !== S.view) return;
    if (!hit || !S.cs || !S.cs.record) { hideCard(); return; }
    const key = hit.type + hit.index + '|' + S.view + '|' + (S.view === '3d' ? (S.cs.b3d && S.cs.b3d.cur) : '');
    if (key !== cardKey) {
      cardKey = key;
      const rows = rowsFor(hit, false);
      if (!rows.length) { hideCard(); return; }
      const head = rows[0], rest = rows.slice(1);
      card.innerHTML = `<div class="mv-card-h">${esc(head.value)}</div>` +
        rest.map(rw => `<div class="mv-card-r"><span class="mv-k">${esc(rw.label)}</span><span class="mv-v${rw.mark ? ' mv-mono' : ''}">${valueHTML(rw)}${tagHTML(rw.kind)}</span></div>`).join('') +
        `<div class="mv-card-f">${S.measureN && S.view === '3d' ? 'click to pick this atom for the measurement' : 'click to pin every value, with notes and sources'}</div>`;
    }
    card.hidden = false;
    const w = card.offsetWidth, h = card.offsetHeight, vw = window.innerWidth, vh = window.innerHeight;
    let left = x + 18, top = y + 18;
    if (left + w > vw - 8) left = Math.max(8, x - w - 18);
    if (top + h > vh - 8) top = Math.max(8, vh - h - 8);
    card.style.left = left + 'px'; card.style.top = top + 'px';
  }
  function onClick(hit, e, from) {
    if (from !== S.view || !S.cs || !S.cs.record) return;
    if (S.view === '3d' && S.measureN) {
      if (hit && hit.type === 'atom') addPick(hit.index);
      return;
    }
    const hadGroup = !!S.group;
    if (!hit) { S.sel = null; S.group = null; }
    else S.sel = { type: hit.type, index: hit.index };
    applyHighlights();
    renderSelected();
    renderTables();                                 // the table row of the selection
    if (hadGroup && !S.group) renderSections();     // the group row loses its "selected" mark with its halo
    if (S.sel) scrollPanelTo($('molSel'));
  }
  // scroll the side panel (never the page) so `node` is visible: in the stacked layout (narrow windows) the panel is
  // not a scroll container, and scrolling the page would move the model off screen at every click
  function scrollPanelTo(node) {
    const side = $('molSide');
    if (!node || !side) return;
    const oy = getComputedStyle(side).overflowY;
    if (!/auto|scroll/.test(oy) || side.scrollHeight <= side.clientHeight) return;
    const nr = node.getBoundingClientRect(), sr = side.getBoundingClientRect();
    if (nr.top < sr.top) side.scrollTop += nr.top - sr.top - 4;
    else if (nr.bottom > sr.bottom) side.scrollTop += Math.min(nr.bottom - sr.bottom + 4, nr.top - sr.top - 4);
  }
  function selSpec(s) { return s ? (s.type === 'atom' ? { atoms: [s.index], bonds: [] } : { atoms: [], bonds: [s.index] }) : null; }
  function applyHighlights() {
    const sel = selSpec(S.sel);
    const grp = S.picks.length ? { atoms: S.picks.slice(), bonds: [] } : S.group;
    if (S.view === '2d' && v2 && v2info && v2info.cs === S.cs) {
      // same colours as the 3D halos: selection orange, group green (hover stays blue in both)
      const f = (s, color) => s ? { atoms: s.atoms.filter(i => i < v2info.N), bonds: s.bonds.filter(k => k < v2info.M), color } : null;
      v2.highlight('select', f(sel, '#ff8c00')); v2.highlight('group', f(grp, '#16a34a'));
    }
    if (S.view === '3d' && v3 && v3cs === S.cs) { v3.highlight('select', sel); v3.highlight('group', grp); }
  }
  function groupSpec(atoms) {
    const set = new Set(atoms), r = S.cs.record;
    return { atoms: atoms.slice(), bonds: r.bonds.filter(b => set.has(b.a) && set.has(b.b)).map(b => b.index) };
  }
  function toggleGroup(spec, id) {
    S.group = S.group && S.group.id === id ? null : Object.assign(spec, { id });
    applyHighlights();
    renderSections();
  }

  /* ================================================================ labels */
  /* lone pairs per atom, for the 'electrons' labels: an atom whose count could not be worked out (an unusual
     electron count) gets none, and the foot line says how many */
  function lonePairCounts(record) {
    return record.atoms.map(a => (a.isH || a.lonePairs === null ? 0 : a.lonePairs));
  }
  function unknownLonePairs(record) {
    return record.atoms.filter(a => !a.isH && a.lonePairs === null).length;
  }
  /* 'pKa' labels: the acidic atom of each acidic site "pKa 7.15", the basic atom of each basic site "pKaH 4.6" (sites,
     values and sources from js/reactions.js moleculePkas: measured IUPAC value, Hammett estimate or typical value) */
  function pkaLabels(cs) {
    const N = cs.record.atoms.length, out = new Array(N).fill('');
    if (!cs.pka) return out;
    for (const s_ of cs.pka.acid) for (const a of s_.atoms || [s_.atom]) if (a < N) out[a] = 'pKa ' + String(s_.pKa).replace('-', '−');
    for (const s_ of cs.pka.base) for (const a of s_.atoms || [s_.atom]) if (a < N) out[a] = (out[a] ? out[a] + ' / ' : '') + 'pKaH ' + String(s_.pKa).replace('-', '−');
    return out;
  }
  /* 'Steric' / 'Electronic' labels: on each fragment's reacting atom (js/reactions.js effectsOfMol) */
  function effectLabels(cs, kind) {
    const N = cs.record.atoms.length, out = new Array(N).fill('');
    for (const e of cs.effects || []) {
      if (e.fragment === 'F009' || e.atom >= N || out[e.atom]) continue;
      if (kind === 'taft') {
        const t = e.taft || {}, sg = v => (v >= 0 ? '+' : '−') + Math.abs(v).toFixed(2);
        out[e.atom] = t.sigmaStar == null && t.EsMin == null ? '' : `${t.sigmaStar != null ? sg(t.sigmaStar) : '–'}/${t.EsMin != null ? sg(t.EsMin) : '–'}`;   // σ*/Es: short, so it fits beside a crowded atom
        continue;
      }
      out[e.atom] = kind === 'steric' ? e.steric.cls : (e.electronic.charge != null ? `δ ${e.electronic.charge >= 0 ? '+' : '−'}${Math.abs(e.electronic.charge).toFixed(2)}` : '')
        + (kind === 'electronic' && e.electronic.sigma != null ? ` Σσ ${e.electronic.sigma >= 0 ? '+' : '−'}${Math.abs(e.electronic.sigma).toFixed(2)}` : '');
    }
    return out;
  }
  function applyLabels2D() {
    if (!v2 || !v2info || !S.cs || !S.cs.record) return;
    const r = S.cs.record;
    v2.setElectrons(S.labels === 'electrons' ? lonePairCounts(r).slice(0, v2info.N) : null);
    if (S.labels === 'none' || S.labels === 'electrons') { v2.setAnnotations(null); return; }
    const all = S.labels === 'pka' ? pkaLabels(S.cs) : S.labels === 'steric' || S.labels === 'electronic' || S.labels === 'taft' ? effectLabels(S.cs, S.labels) : MI.labelsFor(r, S.labels);
    v2.setAnnotations(all.slice(0, v2info.N));
  }
  function applyLabels3D() {
    if (!v3 || v3cs !== S.cs || !S.cs) return;
    const r = S.cs.record;
    v3.setElectrons(S.labels === 'electrons' ? lonePairCounts(r) : null);
    v3.setLabels(S.labels === 'none' || S.labels === 'electrons' ? null : S.labels === 'pka' ? pkaLabels(S.cs)
      : S.labels === 'steric' || S.labels === 'electronic' || S.labels === 'taft' ? effectLabels(S.cs, S.labels) : MI.labelsFor(r, S.labels));
  }
  function renderLabels() { applyLabels2D(); applyLabels3D(); renderFoot(); }

  /* ================================================================ surface and legend */
  function potentialCharges(r) {
    if (r.atoms.some(a => a.gasteiger === null || !Number.isFinite(a.gasteiger))) return null;
    return Float64Array.from(r.atoms.map(a => a.gasteiger));
  }
  function applySurface() {
    if (!v3 || v3cs !== S.cs) return;
    surfaceResult = null;
    renderLegend();
    if (S.surface === 'none') { v3.setSurface(null); return; }
    const cfg = { mode: S.surface === 'potential' ? 'potential' : 'plain', opacity: 0.75 };
    if (cfg.mode === 'potential') {
      const q = potentialCharges(S.cs.record);
      if (!q) { S.surface = 'plain'; cfg.mode = 'plain'; addBanner3D(S.cs, 'info', 'No potential surface: ' + esc(S.cs.record.gasteigerNote || 'partial charges are not available') + '.'); syncControls(); }
      else cfg.charges = q;
    }
    const cs = S.cs;
    v3.setSurface(cfg).then(res => { if (cs === S.cs) onSurfaceResult(res, null); }, err => { if (cs === S.cs) onSurfaceResult(null, err); });
  }
  function onSurfaceResult(res, err) {
    if (err) {
      surfaceResult = null;
      addBanner3D(S.cs, 'warn', 'The surface was not drawn: ' + esc(err.message || err));
      renderLegend();
      return;
    }
    if (res) surfaceResult = res;
    renderLegend();
  }
  /* What the potential legend says — one place for the page's legend and the exported PNG, so they cannot drift
     apart. null when no potential surface is shown. */
  function legendInfo() {
    const res = surfaceResult;
    if (S.view !== '3d' || S.surface !== 'potential' || !res || !S.cs) return null;
    const c = res.clamp, units = res.units || 'kcal/mol per e';
    const fmt = v => (Math.abs(v) >= 10 ? num(v, 0) : num(v, 1));
    const sgn = v => (Number(fmt(v).replace(MINUS, '-')) > 0 ? '+' : '') + fmt(v);   // both ends signed alike: "+141 to +146"
    const conf = currentConf(S.cs);
    if (conf && !geometryOK(S.cs, conf)) {         // an approximate geometry: colours only, no numbers
      return { title: 'Potential on the surface', clamp: c, minLabel: 'δ−', midLabel: '', maxLabel: 'δ+',
        scale: 'red δ− · blue δ+. The geometry is not force-field optimised, so no values are given.', note: MI.TEXT.esp };
    }
    if (!(c > 0)) return { title: `Potential on the surface (${units})`, clamp: 0, minLabel: '', midLabel: '', maxLabel: '', scale: 'The model potential is zero everywhere on this surface.', note: MI.TEXT.esp };
    return { title: `Potential on the surface (${units})`, clamp: c, minLabel: fmt(-c), midLabel: '0', maxLabel: '+' + fmt(c),
      scale: `red δ− · blue δ+ · colours clamped at ±${fmt(c)} (98th percentile of |V|)${res.fixedClamp ? ' (fixed)' : ''}; on this surface ${sgn(res.vmin)} to ${sgn(res.vmax)}`,
      note: MI.TEXT.esp };
  }
  function renderLegend() {
    const L = $('molLegend'), g = legendInfo();
    if (!g) { L.hidden = true; return; }
    const c = g.clamp, col = t => window.Mol3DView.espColor(t * c, c);
    // a strip under the model (not over it): title, colour bar with its ends, then the scale line and the method
    L.innerHTML = `<div class="mv-legend-t">${esc(g.title)}</div>` +
      (c > 0 ? `<div class="mv-legend-scale"><span>${esc(g.minLabel)}</span><span class="mv-legend-bar" style="background:linear-gradient(to right,${[-1, -0.5, 0, 0.5, 1].map(col).join(',')})"></span><span>${esc(g.maxLabel)}</span></div>` : '') +
      `<div class="mv-legend-s">${esc(dot(g.scale))} <span class="mv-legend-n">${esc(dot(g.note))}</span></div>`;
    L.hidden = false;
  }
  const dot = s => (/[.!?]$/.test(String(s)) ? String(s) : s + '.');

  /* ================================================================ measurements (3D) */
  /* A measurement on the conformer shown → {names, value (the panel's text), label (the model's), kind}. A value is
     never given silently where the model is not reliable: an approximate (not force-field optimised) model gives "—",
     as the hover rows do; a measurement that uses an atom of a bond the MMFF94 gap check flags (CO₂, SO₂, N₂O: the
     minimiser had no real parameters for it and stretched it) is marked "(model geometry not reliable here)", on the
     model too; one in a charged or radical system MMFF94 cannot type is marked approximate. */
  const NOT_RELIABLE = MI.TEXT.unreliableHere || '(model geometry not reliable here)';
  const APPROX_SYSTEM = '(approximate: MMFF94 cannot type an atom of this charged or radical system)';
  function gapAtoms(c) {                           // atoms of the bonds the MMFF94 gap check flags on conformer c
    const s = new Set();
    for (const f of (c && c.flaggedBonds) || []) { s.add(f[0]); s.add(f[1]); }
    return s;
  }
  function measureText(m, c) {
    const names = m.atoms.map(i => atomName(S.cs.record, i)).join('–');
    if (!c) return { names, value: '— (no 3D model)', label: '—' };
    if (!geometryOK(S.cs, c)) return { names, value: MI.TEXT.notOptimised, label: '—' };
    const X = c.coords, P = i => [X[3 * i], X[3 * i + 1], X[3 * i + 2]], a = m.atoms;
    const gap = gapAtoms(c), fl = flaggedAtoms(S.cs, c);
    const bad = a.some(i => gap.has(i)), approx = !bad && a.some(i => fl.has(i));
    const mark = bad ? ' ' + NOT_RELIABLE : approx ? ' ' + APPROX_SYSTEM : '', short = bad ? ' ' + NOT_RELIABLE : approx ? ' (approximate)' : '';
    let v, kind, extra = '';
    if (a.length === 2) { const p = P(a[0]), q = P(a[1]); v = Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]).toFixed(3) + ' Å'; kind = 'distance'; }
    else if (a.length === 3) { v = MI.util.angle(P(a[0]), P(a[1]), P(a[2])).toFixed(1) + '°'; kind = 'angle'; }
    else {
      const t = MI.util.dihedral(P(a[0]), P(a[1]), P(a[2]), P(a[3]));
      v = num(t, 1) + '°'; extra = ' (' + MI.util.klynePrelog(t).name + ')'; kind = 'dihedral';
    }
    return { names, value: v + extra + mark, label: v + short, kind, unreliable: bad };
  }
  function applyMeasures() {
    if (!v3 || v3cs !== S.cs) return;
    const c = currentConf(S.cs);
    v3.setMeasurements(c ? S.measures.map(m => ({ atoms: m.atoms, text: measureText(m, c).label })) : null);
  }
  function addPick(i) {
    if (S.picks[S.picks.length - 1] === i) return;
    if (S.picks.includes(i)) return;
    S.picks.push(i);
    if (S.picks.length >= S.measureN) { S.measures.push({ atoms: S.picks.slice() }); S.picks = []; applyMeasures(); }
    applyHighlights();
    renderMeasures(); renderFoot();
  }

  /* ================================================================ banners, resolved line, foot */
  // text that can be a name or a CAS number (not a SMILES, not MOL text)
  const lookupable = t => !!t && !/[\n=#\[\]@]/.test(t) && (/[a-z]{3}/i.test(t) || /^\d{2,7}-\d\d-\d$/.test(t));
  // a structure found by name elsewhere (offline library, PubChem): shown as SMILES, the label says where it came from
  function loadResolved(smiles, note) {
    return showMolecule(smiles, { source: 'pubchem', noFallback: true }).then(r => {
      if (r.ok && S.info) { S.info.label = note; renderResolved(); }
      return r;
    });
  }
  // several different compounds under one name or CAS number: the user chooses, nothing is picked for them
  function showChoices(res, text) {
    const box = $('molResolved');
    box.innerHTML = `<span class="mv-err">${esc(res.error)}</span> `;
    for (const c of res.choices) {
      const b = document.createElement('button');
      b.className = 'mv-pc-lookup'; b.textContent = c.title + (c.formula ? ` (${c.formula})` : '');
      b.title = c.smiles;
      b.addEventListener('click', () => loadResolved(c.smiles, `From ${res.src || 'PubChem'}: CID ${c.cid}, ${c.title} (chosen from ${res.choices.length} records under “${text}”${res.casNote || ''}).`));
      box.append(b, ' ');
    }
  }
  function showError(msg, lookupText) {
    const box = $('molResolved');
    box.innerHTML = `<span class="mv-err">${esc(msg)}</span>` + (S.cs && S.cs.record ? ' <span class="mv-dim">(the molecule shown below is unchanged)</span>' : '');
    // a name or CAS number no offline source knows: offer PubChem (only on a click; it goes over the internet)
    const P = window.PubChem3D;
    if (lookupable(lookupText) && P && P.lookupName) {
      const b = document.createElement('button');
      b.className = 'mv-pc-lookup'; b.textContent = 'Look it up on PubChem';
      b.title = 'Sends the text to PubChem (NCBI) over the internet; the answer is kept in this browser.';
      b.addEventListener('click', () => pubchemName(lookupText));
      box.append(' ', b);
    }
  }
  async function pubchemName(text) {
    const box = $('molResolved'), P = window.PubChem3D;
    box.innerHTML = '<span class="mv-dim">Asking PubChem…</span>';
    const res = await P.lookupName(text);
    if (res.ok) { $('molInput').value = text; loadResolved(res.smiles, res.note + (res.saved ? ' (from your saved lookups)' : '')); return; }
    if (res.choices && res.choices.length) { showChoices(Object.assign({ src: 'PubChem' }, res), text); return; }
    showError(res.error);
  }
  function renderResolved() {
    const box = $('molResolved'), I = S.info;
    if (!I) { box.innerHTML = ''; return; }
    if (I.label) { box.innerHTML = `<span class="mv-via">${esc(I.label)}</span>`; return; }
    if (I.source === 'sketcher' && !I.matched) { box.innerHTML = '<span class="mv-via">From the Sketcher drawing</span>'; return; }
    const how = { 'common name': 'a common name', 'molecule library': "the app's molecule library", 'CAS number': 'a CAS number', 'IUPAC name': 'an IUPAC name',
      'closest name': 'the closest known name', acronym: 'an acronym', 'AGILES code': 'an AGILES code' }[I.via];
    let h = `<span class="mv-via">` + (I.matched && how
      ? `Read as <b>${esc(I.matched)}</b> (from ${esc(how)})` : `Read as ${esc(I.via)}${I.matched && I.matched !== I.via ? ': <b>' + esc(I.matched) + '</b>' : ''}`) + '</span>';
    if (I.warn) h += ` <span class="mv-warn">${esc(I.warn)}</span>`;
    box.innerHTML = h;
  }
  function addBanner3D(cs, kind, html) {
    if (!cs) return;
    const B = ensureB3D(cs);
    if (!B.banners.some(b => b.html === html)) B.banners.push({ kind, html });
    if (cs === S.cs) renderBanners();
  }

  /* ---------- stereochemistry the input leaves open (one banner) ---------- */
  /* → {atoms, bonds} (record indices) whose configuration the input does not give: exactly the open stereo the
     record reports and the Stereochemistry section lists (MolInfo: RDKit's '?' centres, open cis/trans ring stereo,
     unspecified double bonds and axes), in 2D and 3D alike. Never the 3D builder's own list: it also names centres
     MolInfo rules out (ATP's phosphate P, a norbornene bridgehead: stereo.ignored), and the banner and the panel must
     say the same thing. */
  function openStereo(cs) {
    const r = cs.record, st = r.stereo;
    const ignored = new Set((st.ignored || []).map(x => x.index));
    const atoms = st.atoms.filter(s => s.label === '?' && !ignored.has(s.index)).map(s => s.index).sort((a, b) => a - b);
    const bonds = st.bonds.filter(s => s.label === 'unspecified').map(s => s.index).sort((a, b) => a - b);
    return { atoms, bonds };
  }
  const bondText = (r, k) => { const b = r.bonds[k]; return atomName(r, b.a) + (b.order === 2 && !b.aromatic ? '=' : '–') + atomName(r, b.b) + (b.order === 2 && !b.aromatic ? '' : ' (axis)'); };
  /* The configuration a conformer has, read back from its coordinates by RDKit (stereo from 3D, new CIP labeller),
     in record indices — cached on the conformer (its coordinates only ever turn rigidly, which keeps it). */
  function modelStereo(cs, c) {
    if (!c) return null;
    if (c.stereoRead !== undefined && c.stereoRead !== null) return c.stereoRead;
    let m = null, out = null;
    try {
      m = S.R.get_mol(window.Mol3DEngine.molfile(cs.record, c.coords, 'model'), JSON.stringify({ removeHs: false }));
      if (m) {
        const t = JSON.parse(m.get_stereo_tags() || '{}'), atoms = new Map(), bonds = new Map();
        for (const [i, lab] of t.CIP_atoms || []) atoms.set(i, String(lab).replace(/[()]/g, ''));
        for (const [a, b, lab] of t.CIP_bonds || []) bonds.set(pairKey(a, b), String(lab).replace(/[()]/g, ''));
        out = { atoms, bonds };
      }
    } catch (e) { out = null; } finally { if (m) { try { m.delete(); } catch (e) {} } }
    c.stereoRead = out || false;
    return out;
  }
  /* ---------- cis/trans of open ring stereo in the model shown ---------- */
  /* CIP names ring stereo only as r/s or R/S ("C4 (s, pseudoasymmetric)" for decalin), which does not say cis or
     trans. For the record's open ring stereo, the model's cis/trans is read from its coordinates:
     - a ring fusion a–b (both open fusion atoms, each with one exocyclic neighbour: H or an angular substituent): the
       dihedral X–a–b–Y is gauche for a cis fusion (decalin about 60°, less in five- and four-membered rings) and anti
       for a trans one (decalin about 180°): |τ| ≤ 90° cis, |τ| ≥ 130° trans;
     - ring substituents (an open ring atom with two ring bonds, one H and one other substituent): the side of the
       ring that (substituent − H) points to, against the ring's mean-plane (Newell) normal; the same side = cis.
       Open stereocentres of that kind are included too (1,2-dimethylcyclohexane: C2 (S), C7 (R) is the cis form)
       when every open centre of the ring is one and there are at most four, and so are open stereocentres at a ring
       fusion (hydrindane's C4 and C8, which RDKit names R/S rather than r/s).
     Nothing is said where the test is not clear-cut: 90° < |τ| < 130°, an (X − H) within 20° of the ring plane (that
     ring then gets no statement at all), rings of more than 8 atoms, or an atom without exactly one H.
     Checked on models of wedge-drawn cis and trans decalin, hydrindane, 1,2-, 1,3- and 1,4-dimethylcyclohexane,
     1,3-dimethylcyclopentane and 1,2-dimethylcyclopropane, and on the banner of each open structure against RDKit's
     reading of its model (molview_selftest.html).
     cand: [{index, ring: 'fusion'|'substituent'}] → [text] */
  const RING_FACE_MAX = 8, FACE_MIN = Math.sin(20 * Math.PI / 180), CIS_MAX = 90, TRANS_MIN = 130;
  const inRingOrder = (g, a, b) => { const ia = g.indexOf(a), ib = g.indexOf(b); return ia >= 0 && ib >= 0 && (Math.abs(ia - ib) === 1 || Math.abs(ia - ib) === g.length - 1); };
  const listText = xs => (xs.length > 1 ? xs.slice(0, -1).join(', ') + ' and ' + xs[xs.length - 1] : xs[0] || '');
  const byNumber = (p, q) => p.localeCompare(q, 'en', { numeric: true });
  function ringCisTrans(r, X, cand) {
    const P = i => [X[3 * i], X[3 * i + 1], X[3 * i + 2]];
    const ringNb = i => r.atoms[i].neighbours.filter(n => r.bonds[n.bond].inRing).map(n => n.atom);
    const exo = i => r.atoms[i].neighbours.filter(n => !r.bonds[n.bond].inRing).map(n => n.atom);
    const out = [];
    const fus = new Set(cand.filter(x => x.ring === 'fusion').map(x => x.index)), seen = new Set();
    for (const a of [...fus].sort((p, q) => p - q)) for (const b of ringNb(a)) {
      if (!fus.has(b) || seen.has(pairKey(a, b))) continue;
      seen.add(pairKey(a, b));
      const sizes = r.rings.filter(g => inRingOrder(g, a, b)).map(g => g.length).sort((p, q) => p - q);
      const xa = exo(a), xb = exo(b);
      if (sizes.length < 2 || sizes[1] > RING_FACE_MAX || xa.length !== 1 || xb.length !== 1) continue;
      const t = Math.abs(MI.util.dihedral(P(xa[0]), P(a), P(b), P(xb[0])));
      const rel = t <= CIS_MAX ? 'cis' : t >= TRANS_MIN ? 'trans' : '';
      if (rel) out.push(`the ring fusion ${atomName(r, Math.min(a, b))}–${atomName(r, Math.max(a, b))} is ${rel}`);
    }
    // ring substituents, ring by ring (the smallest ring through each atom)
    const byRing = new Map();
    for (const x of cand) {
      if (x.ring !== 'substituent') continue;
      const rs = r.rings.map((g, t) => [g, t]).filter(([g]) => g.includes(x.index)).sort((p, q) => p[0].length - q[0].length);
      if (!rs.length || rs[0][0].length > RING_FACE_MAX) continue;
      const t = rs[0][1];
      if (!byRing.has(t)) byRing.set(t, []);
      byRing.get(t).push(x.index);
    }
    for (const [t, atoms] of byRing) {
      if (atoms.length < 2) continue;
      const g = r.rings[t], pts = g.map(P), c = [0, 1, 2].map(d => pts.reduce((s, p) => s + p[d], 0) / pts.length);
      const n = [0, 0, 0];
      for (let k = 0; k < pts.length; k++) {           // Newell's normal of the ring polygon (area-weighted, any pucker)
        const u = pts[k].map((v, d) => v - c[d]), w = pts[(k + 1) % pts.length].map((v, d) => v - c[d]);
        n[0] += u[1] * w[2] - u[2] * w[1]; n[1] += u[2] * w[0] - u[0] * w[2]; n[2] += u[0] * w[1] - u[1] * w[0];
      }
      const nn = Math.hypot(n[0], n[1], n[2]);
      if (!(nn > 1e-6)) continue;
      const side = new Map();
      let clear = true;
      for (const i of atoms) {
        const e = exo(i), h = e.filter(j => r.atoms[j].el === 'H'), x = e.filter(j => r.atoms[j].el !== 'H');
        if (e.length !== 2 || h.length !== 1 || x.length !== 1) { clear = false; break; }
        const ph = P(h[0]), v = P(x[0]).map((q, d) => q - ph[d]);
        const f = (v[0] * n[0] + v[1] * n[1] + v[2] * n[2]) / (nn * Math.hypot(v[0], v[1], v[2]));
        if (!(Math.abs(f) >= FACE_MIN)) { clear = false; break; }
        side.set(i, f > 0 ? 1 : -1);
      }
      if (!clear) continue;
      const up = atoms.filter(i => side.get(i) > 0).map(i => atomName(r, i)).sort(byNumber);
      const down = atoms.filter(i => side.get(i) < 0).map(i => atomName(r, i)).sort(byNumber);
      if (atoms.length === 2) out.push(`${listText(up.concat(down).sort(byNumber))} are ` + (up.length === 2 || down.length === 2 ? 'cis (on the same side of the ring)' : 'trans (on opposite sides of the ring)'));
      else if (!up.length || !down.length) out.push(`${listText(up.length ? up : down)} are all on one side of the ring (all cis)`);
      else out.push(`${listText(up)} ${up.length === 1 ? 'is' : 'are'} on one side of the ring and ${listText(down)} on the other`);
    }
    return out;
  }
  // the record's open ring stereo that ringCisTrans reads (see there) → [{index, ring}]
  function ringCandidates(cs) {
    const r = cs.record, ignored = new Set((r.stereo.ignored || []).map(x => x.index));
    const open = r.stereo.atoms.filter(s => s.label === '?' && !ignored.has(s.index));
    const ringDeg = i => r.atoms[i].neighbours.filter(n => r.bonds[n.bond].inRing).length;
    const out = open.filter(s => s.kind === 'ring' && (s.ring === 'fusion' || s.ring === 'substituent')).map(s => ({ index: s.index, ring: s.ring }));
    // open stereocentres in a ring (1,2-dimethylcyclohexane): per ring, only when every open centre of that ring is a
    // ring substituent of the simple kind (two ring bonds, one H, one other substituent) and there are at most four
    const simple = i => { const e = r.atoms[i].neighbours.filter(n => !r.bonds[n.bond].inRing).map(n => n.atom);
      return ringDeg(i) === 2 && e.length === 2 && e.filter(j => r.atoms[j].el === 'H').length === 1; };
    const centres = open.filter(s => s.kind === 'centre' && ringDeg(s.index) > 0);
    // open stereocentres at a ring fusion (hydrindane: RDKit names them R/S, not r/s, so the record lists them as
    // centres): three ring bonds and one other neighbour; ringCisTrans pairs them by their shared bond
    for (const s of centres) {
      const exoN = r.atoms[s.index].neighbours.filter(n => !r.bonds[n.bond].inRing).length;
      if (ringDeg(s.index) === 3 && exoN === 1 && !out.some(x => x.index === s.index)) out.push({ index: s.index, ring: 'fusion' });
    }
    for (const g of r.rings) {
      const inG = centres.filter(s => g.includes(s.index));
      if (inG.length < 2 || inG.length > 4 || !inG.every(s => simple(s.index))) continue;
      for (const s of inG) if (!out.some(x => x.index === s.index)) out.push({ index: s.index, ring: 'substituent' });
    }
    return out;
  }

  /* ---------- open axial stereo (allene, hindered biaryl, spiro compound, alkylidene ring) ---------- */
  /* RDKit's stereo perception, and so the record, has no axial stereo elements, so the record's open list cannot
     name them; the panel says "not determined: … axial chirality possible". The 3D builder does report them
     (result.stereo.unspecified: the atoms or bonds OpenChemLib treats as stereo), and its model shows one form. Its
     entries that the record neither lists nor rules out (stereo.ignored) are grouped into one element per axis: an
     allene's centre → the allene; a single bond → the axis about it; an exocyclic C=C and the listed atoms of its ring
     (OpenChemLib lists both, "C2, C5=C6", for the one axis of an alkylidenecyclohexane) → that axis; ring atoms of a
     spiro compound → its spiro centre. Anything else is named as the builder names it. → [{name, axial}] */
  function axialOpen(cs) {
    const B = cs.b3d, st = B && B.result && B.result.stereo, r = cs.record;
    if (!st || !Array.isArray(st.unspecified) || !st.unspecified.length) return [];
    const listedA = new Set(r.stereo.atoms.map(s => s.index).concat((r.stereo.ignored || []).map(x => x.index)));
    const listedB = new Set(r.stereo.bonds.map(s => s.index));
    const bondOf = (a, b) => r.bonds.findIndex(x => (x.a === a && x.b === b) || (x.a === b && x.b === a));
    let atoms = [], bonds = [];
    for (const u of st.unspecified) {
      if (u.atom != null) { if (u.atom >= 0 && u.atom < r.N && !listedA.has(u.atom)) atoms.push(u.atom); }
      else if (u.bond) { const k = bondOf(u.bond[0], u.bond[1]); if (k >= 0 && !listedB.has(k)) bonds.push(k); }
    }
    atoms = [...new Set(atoms)].sort((p, q) => p - q); bonds = [...new Set(bonds)].sort((p, q) => p - q);
    if (!atoms.length && !bonds.length) return [];
    const out = [], usedA = new Set(), ringsOf = i => r.rings.filter(g => g.includes(i));
    const dbl = i => r.atoms[i].neighbours.filter(n => n.order === 2 && !r.bonds[n.bond].aromatic).map(n => n.atom);
    for (const k of bonds) {
      const b = r.bonds[k];
      if (b.order === 2 && !b.aromatic && !b.inRing) {
        const r0 = [b.a, b.b].find(x => ringsOf(x).length);
        const mates = r0 == null ? [] : atoms.filter(i => !usedA.has(i) && ringsOf(r0).some(g => g.includes(i)));
        if (mates.length) {
          mates.forEach(i => usedA.add(i));
          out.push({ name: `the axis through ${listText(mates.map(i => atomName(r, i)))} and ${atomName(r, b.a)}=${atomName(r, b.b)} (an alkylidene ring)`, axial: true });
        } else out.push({ name: `${atomName(r, b.a)}=${atomName(r, b.b)}`, axial: false });
      } else if (b.order === 1 && !b.aromatic) out.push({ name: `the ${atomName(r, b.a)}–${atomName(r, b.b)} axis`, axial: true });
      else out.push({ name: bondText(r, k), axial: false });
    }
    // spiro centres: an atom in two rings that share only it
    const spiro = [];
    for (let s = 0; s < r.n; s++) {
      const rs = ringsOf(s);
      if (rs.some((g, x) => rs.some((h, y) => y > x && g.filter(v => h.includes(v)).length === 1))) spiro.push(s);
    }
    const bySpiro = new Map();
    for (const i of atoms) {
      if (usedA.has(i)) continue;
      const d = dbl(i);
      if (d.length === 2) {                          // the centre of an allene (or a longer cumulene): walk to both ends
        const sides = d.map(start => {
          const part = [];
          let prev = i, cur = start;
          while (cur != null && cur !== i && !part.includes(cur)) {
            part.push(cur);
            const nx = dbl(cur).filter(x => x !== prev);
            if (dbl(cur).length !== 2 || nx.length !== 1) break;
            prev = cur; cur = nx[0];
          }
          return part;
        });
        const chain = sides[0].slice().reverse().concat([i], sides[1]);
        usedA.add(i);
        out.push({ name: `the ${chain.length === 3 ? 'allene' : 'cumulene'} ${chain.map(x => atomName(r, x)).join('=')}`, axial: true });
        continue;
      }
      const s = spiro.find(x => x !== i && ringsOf(i).some(g => g.includes(x)));
      if (s != null) { if (!bySpiro.has(s)) bySpiro.set(s, []); bySpiro.get(s).push(i); usedA.add(i); }
    }
    for (const [s, list] of bySpiro) out.push({ name: `the spiro centre ${atomName(r, s)} (the axis through ${listText(list.map(i => atomName(r, i)))})`, axial: true });
    for (const i of atoms) if (!usedA.has(i)) out.push({ name: atomName(r, i), axial: false });
    return out;
  }
  function axialBanner(cs, c) {
    if (!c) return '';
    const els = axialOpen(cs);
    if (!els.length) return '';
    const r = cs.record, names = listText(els.map(e => e.name));
    const who = c.source === 'PubChem' ? "PubChem's 3D model" : c.source === 'file' ? "the file's 3D structure" : 'the 3D model';
    const chose = c.source === 'OpenChemLib' ? ' (the builder chose it)' : '';
    // mirror images only when this axis is the molecule's one stereo element; otherwise the forms are diastereomers
    const alone = els.length === 1 && els[0].axial && !r.stereo.atoms.length && !r.stereo.bonds.length;
    return `The input does not specify the configuration of ${names}, so ${who} shows ` +
      (alone ? `one of its two mirror-image forms${chose}.` : `one of the possible configurations${els.length > 1 ? ' of each' : ''}${chose}.`);
  }

  /* the ways a structure comes from the app's name library (Chem.searchMolecule's via): such a name may stand for one
     stereoisomer the library does not hold, so the input did not leave the stereochemistry open, the library did.
     An IUPAC name, a SMILES, a file or a drawing states its own stereochemistry. */
  const NAME_VIAS = new Set(['molecule library', 'common name', 'CAS number', 'acronym', 'closest name']);
  function stereoBanner(cs) {
    const r = cs.record, is3 = S.view === '3d', c = is3 ? currentConf(cs) : null;
    const open = openStereo(cs);
    if (!open.atoms.length && !open.bonds.length) return '';
    const all = open.atoms.map(i => atomName(r, i)).concat(open.bonds.map(k => bondText(r, k)));
    const where = all.length > 12 ? `${all.length} places (${all.slice(0, 10).join(', ')} and ${all.length - 10} more)` : all.join(', ');
    const I = S.info, byName = !!(I && NAME_VIAS.has(I.via));
    const anySet = r.stereo.atoms.some(s => s.label !== '?') || r.stereo.bonds.some(s => s.label !== 'unspecified');
    const head = (byName ? `The name was resolved to a structure without stereochemistry at ${where} (the app's name library has none ${anySet ? 'there' : 'for it'}).`
      : `Stereochemistry not specified in the input at ${where}.`) + ' Values computed from the structure (formula, masses, counts, InChI) hold for any stereoisomer; ';
    if (!c) return head + 'a 3D model has to pick one.';
    const got = modelStereo(cs, c);
    const lab = (x, p) => (x === 'r' || x === 's' ? `${p} (${x}, pseudoasymmetric)` : `${p} (${x})`);
    const named = open.atoms.map(i => { const x = got && got.atoms.get(i); return x && x !== '?' ? lab(x, atomName(r, i)) : atomName(r, i) + ' (not named by CIP)'; })
      .concat(open.bonds.map(k => { const b = r.bonds[k], x = got && got.bonds.get(pairKey(b.a, b.b)); return x ? `${bondText(r, k)} (${x})` : bondText(r, k) + ' (one of its forms)'; }));
    const who = c.source === 'PubChem' ? "PubChem's 3D model" : c.source === 'file' ? "the file's 3D structure" : 'the 3D model (the builder picked one)';
    let t = head + `${who} and its measurements are for one of them: ${named.length > 12 ? named.slice(0, 10).join(', ') + ` and ${named.length - 10} more` : named.join(', ')}`;
    // cis/trans of the open ring stereo, where the model's geometry shows it unambiguously
    let ct = [];
    try { ct = ringCisTrans(r, c.coords, ringCandidates(cs)); } catch (e) { ct = []; }
    t += ct.length ? `; in this model ${ct.join('; ')}.` : '.';
    const B = cs.b3d;
    if (c.source === 'PubChem' && B && B.pub && B.pub.identity && B.pub.identity.mixedStereo) t += " PubChem's other conformers of this record mix stereoisomers, so only its default conformer is used.";
    return t;
  }
  /* the engine's own open-stereo list ("stereochemistry at … / the axial (atropisomeric) configuration about … is not
     specified in the input"): replaced by stereoBanner, which names the record's list and the model's form. Every
     other engine warning is shown, among them the cumulene one ("E/Z at the cumulene C2=C3=C4=C5 is not taken from
     the input…"): the record has no E/Z for a cumulene, so nothing else says that the model picked one. */
  const isStereoWarning = w => /^(stereochemistry at |the axial \(atropisomeric\) configuration)/i.test(w) && /not specified in the input/i.test(w);
  const sentence = s => { const t = String(s).trim(); return dot(t.charAt(0).toUpperCase() + t.slice(1)); };

  function renderBanners() {
    const box = $('molBanners'), cs = S.cs, out = [];
    const add = (kind, html) => out.push(`<div class="mv-banner mv-${kind}">${html}</div>`);
    if (S.info) for (const n of S.info.notes) add(n.kind, n.html);
    if (cs && cs.record) {
      if (S.comps.length > 1) add('info', `This input has ${S.comps.length} components (a salt or mixture); each is analysed on its own. Showing component ${cs.i + 1}: ${esc(cs.formula)}.` +
        (S.combined ? ` Whole input: ${esc(S.combined.formulaText || '')}, ${S.combined.mw.rounded.toFixed(2)} g/mol.` : '') +
        (cs.from3D ? " The 3D view rebuilds each component's model; the file's own coordinates are shown only for a file with one molecule." : ''));
      if (cs.from3D) {
        const F = cs.fileConf;
        add('info', "The file's 3D coordinates were used to read the stereochemistry; the 2D drawing is laid out by RDKit." +
          (F && F.coords ? ' In 3D the file\'s own coordinates are the conformer "File" (as given: not minimised, no energy); OpenChemLib\'s conformers are listed with it.'
            : F && F.error ? ` In 3D the model is rebuilt by OpenChemLib; the file's own coordinates are not shown (${esc(F.error)}).` : ''));
      }
      const st = stereoBanner(cs);
      if (st) add('warn', esc(st));
      if (S.view === '3d') { const ax = axialBanner(cs, currentConf(cs)); if (ax) add('warn', esc(ax)); }
      if (S.view === '2d' && S.showH2d && cs.drawHError) add('error', 'Hydrogens could not be drawn: ' + esc(cs.drawHError));
    }
    if (S.view === '3d' && cs && cs.b3d) {
      const B = cs.b3d, c = currentConf(cs);
      if (B.error) add('error', 'No 3D model from OpenChemLib: ' + esc(B.error));
      // a lone atom has only the single-atom note (foot, conformer list): its "geometry" is a point, so MMFF94 typing
      // failures ([2H]: no type for a free H atom) and radical flags say nothing about it
      if (c && c.source === 'OpenChemLib' && !singleAtomNote(cs, c)) {
        for (const w of new Set(B.engineWarnings || [])) if (!isStereoWarning(w)) add('warn', esc(sentence(w)));
        if (!c.optimised) {
          const why = c.warnings.map(w => w.replace(/^approximate geometry \(not force-field optimised\):\s*/, '')).join('; ');
          add('warn', esc(MI.TEXT.notOptimisedBanner.replace(/^a/, 'A')) + (why ? '. Reason: ' + esc(why) + '.' : '.'));
        }
        if (c.flaggedBonds.length) add('warn', flaggedBanner(cs, c));
      }
      for (const b of B.banners) if (!b.stereo) add(b.kind, b.html);
      if (contextLost) add('warn', 'The 3D graphics context was lost; the model comes back when the browser restores it.');
    }
    box.innerHTML = out.join('');
  }
  /* Bonds the MMFF94 gap check flags (CO₂, SO₂, N₂O: the minimiser had no real parameters and stretched them): the
     engine's reason for each, then what that means for the rest of the model, as the panel and Measure show it */
  function flaggedBanner(cs, c) {
    const r = cs.record, atoms = [...gapAtoms(c)].sort((p, q) => p - q).map(i => atomName(r, i));
    const byText = new Map();                        // bonds with the same reason (CO₂'s two C=O) share one sentence
    for (const f of c.flaggedBonds) { if (!byText.has(f[2])) byText.set(f[2], []); byText.get(f[2]).push(`${atomName(r, f[0])}–${atomName(r, f[1])}`); }
    const per = [...byText].map(([t, bs]) => (bs.length === 1 ? `Bond ${bs[0]}: ` : `Bonds ${listText(bs)}, each: `) + sentence(t)).join(' ');
    const g = cs.geo, withheld = !!(g && g.molecule && g.molecule.withheld), hidden = withheld || !molLevel3DShown(cs);
    const or = xs => (xs.length > 1 ? xs.slice(0, -1).join(', ') + ' or ' + xs[xs.length - 1] : xs[0] || '');
    return esc(per + ` The model geometry is not reliable near ${c.flaggedBonds.length === 1 ? 'this bond' : 'these bonds'}: bond angles and dihedrals ` +
      `that involve ${or(atoms)}, and any measurement that uses ${atoms.length === 1 ? 'it' : 'one of them'}, are marked "${NOT_RELIABLE.replace(/^[(]|[)]$/g, '')}"` +
      (hidden ? ', and no molecule-level 3D values are shown (longest extent, radius of gyration, principal moments, shape, van der Waals volume and area' +
        (withheld ? ', dipole verdict' : '') + ').' : '.'));
  }
  // does the panel show molecule-level numbers measured on the model (MolInfo leaves them out for a flagged model)?
  function molLevel3DShown(cs) {
    try {
      const secs = MI.moleculeSections(cs.record, cs.geo, extraFor(cs));
      return secs.some(sec => sec.rows.some(rw => /^(rg|pmi|shape|extent|volume|area)$/.test(rw.key) && /\d/.test(rw.value)));
    } catch (e) { return true; }
  }
  function renderFoot() {
    const box = $('molFoot'), cs = S.cs, parts = [];
    if (!cs || !cs.record) { box.innerHTML = ''; return; }
    const r = cs.record;
    if (S.view === '3d') {
      const c = currentConf(cs), B = cs.b3d;
      if (c && c.source === 'OpenChemLib') {
        const one = singleAtomNote(cs, c);
        parts.push(`Conformer ${c.rank} of ${c.of}` + (one ? ' · ' + esc(one) : hasEnergy(cs, c) ? ` · MMFF94 energy ${energyText(c)}` : ' · not force-field optimised') +
          ' · OpenChemLib' + (B && B.done && B.result && B.result.generated ? ` · ${B.result.generated} generated` : ''));
      } else if (c && c.source === 'file') {
        parts.push(`Conformer "File": the coordinates of the opened file, as given (not minimised; no energy)${S.info && S.info.fileName ? ' · ' + esc(S.info.fileName) : ''}`);
      } else if (c && c.source === 'PubChem' && B.pub) {
        const A = B.pub.attribution;
        parts.push(`${esc(A.text)} · <a href="${esc(A.url)}" target="_blank" rel="noopener">record</a> · <a href="${esc(A.policyUrl)}" target="_blank" rel="noopener">${esc(A.policyText)}</a>`);
        if (B.pub.energyNote) parts.push(esc(B.pub.energyNote));
      }
      if (S.measureN) parts.push(`<b>Measure:</b> click ${S.measureN} atoms (${S.picks.length} picked) — Esc cancels`);
    } else {
      parts.push('2D drawing · hover an atom or bond for its values, click to pin them · wheel zooms, drag pans, double-click fits');
    }
    if (S.labels !== 'none') {
      let t = 'Labels: ' + LABEL_TEXT[S.labels];
      if (S.labels === 'pka') {
        const P = cs.pka, F = ((window.REACTION_RULES || {}).fragments) || {};
        const item = (s_, kind) => `${esc([...new Set(s_.fragments)].map(f => F[f] ? F[f].name : f).join(' / '))} ${kind} ${esc(String(s_.pKa))} (${esc(s_.source)})`;
        const list = P ? P.acid.map(s_ => item(s_, 'pKa')).concat(P.base.map(s_ => item(s_, 'pKaH'))) : [];
        t += list.length ? ': ' + list.join(' · ') : ': no acidic or basic site found by the Fragments sheet';
      }
      if (S.labels === 'steric' || S.labels === 'electronic' || S.labels === 'taft') {
        const taftText = e => [e.taft && e.taft.electronicText, e.taft && e.taft.stericText].filter(Boolean).join('; ') || 'no group with a Taft value on it';
        const list = (cs.effects || []).filter(e => e.fragment !== 'F009').map(e => `${esc(e.name)} (atom ${e.atom + 1}): ${esc(S.labels === 'steric' ? e.steric.text : S.labels === 'taft' ? taftText(e) : e.electronic.text)}`);
        t += list.length ? ': ' + list.join(' · ') : ': no fragment of the Fragments sheet found';
      }
      if (S.labels === 'electrons') {
        const unk = unknownLonePairs(r);
        t += ' · ' + MI.TEXT.asDrawn + (unk ? ` · ${unk} atom${unk === 1 ? '' : 's'} with an unusual electron count get none` : '');
        if (S.view === '2d') t += ' · only atoms whose symbol is drawn (a plain carbon has no lone pair)';
      }
      if (S.labels === 'cip') {
        const ez = MI.bondLabelsFor(r, 'cip').map((l, k) => (l ? `${atomName(r, r.bonds[k].a)}=${atomName(r, r.bonds[k].b)} ${l === '?' ? 'E/Z not specified' : '(' + l + ')'}` : null)).filter(Boolean);
        if (ez.length) t += ' · double bonds: ' + ez.join(', ');
      }
      parts.push(esc(t));
    }
    box.innerHTML = parts.join('<br>');
  }
  // "52.34 kcal/mol (ΔE 0.04 vs the lowest found)" — ΔE to 2 dp as in the conformer table: at 1 dp a conformer 0.04
  // above the lowest read "ΔE 0.0", as if it were the lowest
  function energyText(c) {
    return `${num(c.energy, 2)} kcal/mol (` + (c.rank === 1 ? 'the lowest found' : `ΔE ${num(c.relEnergy, 2)} vs the lowest found`) + ')';
  }

  /* ================================================================ side panel */
  function rowHTML(rw, opts) {
    opts = opts || {};
    const note = rw.note ? `<div class="mv-note">${esc(rw.note)}</div>` : '';
    const mono = /^(smiles|smiles0|smilespos|inchi|inchikey)$/.test(rw.key) ? ' mv-mono' : '';
    return `<div class="mv-row${rw.note ? ' mv-has-note' : ''}${opts.click ? ' mv-click' : ''}${opts.on ? ' mv-on' : ''}"${opts.attrs || ''}>` +
      `<span class="mv-k">${esc(rw.label)}${rw.note ? '<span class="mv-i" title="show the note">i</span>' : ''}</span>` +
      `<span class="mv-v${mono}${copyable(rw) ? ` mv-copy" data-copy="${esc(copyValue(rw))}" title="${rw.mark ? 'Click to copy the underlined part' : 'Click to copy'}` : ''}">${valueHTML(rw)}${tagHTML(rw.kind)}</span>${note}</div>`;
  }
  // the naming rows (IUPAC and common name, CAS, SMILES, InChI, InChIKey): click the value to copy it
  const copyable = rw => /^(iupac|common|cas|formula|isoformula|smiles|smiles0|smilespos|inchi|inchikey)$/.test(rw.key) && rw.value && rw.value !== '—' && !/^not defined/.test(rw.value);
  // formulas are shown with sub/superscript digits (C₉H₈O₄); they are copied as plain text (C9H8O4) so they paste anywhere
  const SUBSUP = { '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4', '₅': '5', '₆': '6', '₇': '7', '₈': '8', '₉': '9',
    '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4', '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9', '⁺': '+', '⁻': '-' };
  const plainFormula = t => String(t).replace(/[₀-₉⁰¹²³⁴-⁹⁺⁻]/g, c => SUBSUP[c] || c);
  // "Location in SMILES" copies only what is underlined (several underlined pieces, as for a ring bond, joined in order)
  const copyValue = rw => {
    if (/formula$/.test(rw.key)) return plainFormula(rw.value);
    if (!rw.mark) return rw.value;
    const marks = typeof rw.mark[0] === 'number' ? [rw.mark] : rw.mark;
    return marks.slice().sort((p, q) => p[0] - q[0]).map(([a, b]) => String(rw.value).slice(a, b)).join('');
  };
  function copyText(text, done) {
    const fallback = () => {
      const ta = document.createElement('textarea');
      ta.value = text; ta.style.position = 'fixed'; ta.style.opacity = '0';
      document.body.appendChild(ta); ta.select();
      try { if (document.execCommand('copy')) done(); } catch (e) { /* nothing more to try */ }
      ta.remove();
    };
    if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(done, fallback); else fallback();
  }
  // a row's value; row.mark = [start, end) or a list of them underlines those parts ("Location in SMILES")
  function valueHTML(rw) {
    if (!rw.mark) return esc(rw.value);
    const s = String(rw.value), marks = typeof rw.mark[0] === 'number' ? [rw.mark] : rw.mark;
    let html = '', pos = 0;
    for (const [a, b] of marks.slice().sort((p, q) => p[0] - q[0])) {
      if (a < pos) continue;
      html += esc(s.slice(pos, a)) + '<u class="mv-mark">' + esc(s.slice(a, b)) + '</u>';
      pos = b;
    }
    return html + esc(s.slice(pos));
  }
  function titleOf(cs) {
    const nm = cs.names;
    return (nm && (nm.common || (nm.iupacChecked === false ? null : nm.iupac))) || cs.formula || 'molecule';
  }
  function renderHead() {
    const box = $('molHead'), cs = S.cs;
    if (!cs || !cs.record) { box.innerHTML = '<h2>No molecule</h2>'; return; }
    const r = cs.record, nm = cs.names;
    const title = titleOf(cs);
    const sub = [r.formulaText, r.mw ? r.mw.rounded.toFixed(2) + ' g/mol' : null].filter(Boolean).join(' · ');
    box.innerHTML = `<h2>${esc(title)}</h2><div class="mv-sub">${esc(sub)}${nm ? '' : ' <span class="mv-dim">(looking up names…)</span>'}</div>
      <div class="mv-legendline">Tags: ${['model', 'table', '3d', 'convention', 'unchecked'].map(k => tagHTML(k) + ' ' + esc(KIND[k].title.split(',')[0])).join(' · ')}; untagged values follow exactly from the structure.</div>
      <label class="mv-notes-t"><input type="checkbox" id="molNotesAll"${S.notesOpen ? ' checked' : ''}> show every note (method, source, caveat)</label>`;
    $('molNotesAll').addEventListener('change', e => { S.notesOpen = e.target.checked; $('molSide').classList.toggle('mv-notes-open', S.notesOpen); });
  }
  function renderSelected() {
    const box = $('molSel'), cs = S.cs;
    if (!cs || !cs.record || !S.sel) {
      box.innerHTML = `<div class="mv-sec-h">Selected atom or bond</div><p class="mv-hint">Hover over an atom or bond for its values; click it to pin them here with their notes.</p>`;
      return;
    }
    const rows = rowsFor(S.sel, true);
    box.innerHTML = `<div class="mv-sec-h">${esc(S.sel.type === 'atom' ? 'Selected atom' : 'Selected bond')} <button class="mv-x" data-act="unpin" title="Clear the selection">×</button></div>` +
      rows.map(rw => rowHTML(rw)).join('');
  }
  function extraFor(cs) {
    const ex = { names: cs.names || null };
    if (S.comps.length > 1) ex.component = { index: cs.i, count: S.comps.length, formula: cs.formula };
    if (S.combined) ex.combined = S.combined;
    if (S.view === '3d') {
      const c = currentConf(cs);
      const one = singleAtomNote(cs, c);
      if (c) ex.conformer = c.source === 'PubChem'
        ? { energy: null, relEnergy: null, source: `PubChem CID ${c.cid}${c.conformerId ? ', conformer ' + c.conformerId : ''} (computed by PubChem3D, not experimental)` }
        : c.source === 'file'
          ? { energy: null, relEnergy: null, source: 'the opened file (its own coordinates, as given: not minimised, no energy)' }
          : one
            ? { energy: null, relEnergy: null, source: `OpenChemLib (${one})` }
            : { energy: c.optimised ? c.energy : null, relEnergy: c.optimised ? c.relEnergy : null, source: `OpenChemLib, conformer ${c.rank} of ${c.of}` + (c.optimised ? ' (MMFF94-minimised)' : ' (not force-field optimised)') };
    }
    return ex;
  }
  function renderSections() {
    const box = $('molSections'), cs = S.cs;
    if (!cs || !cs.record) { box.innerHTML = ''; return; }
    const r = cs.record, geo = S.view === '3d' ? cs.geo : null;
    let secs;
    try { secs = MI.moleculeSections(r, geo, extraFor(cs)); } catch (e) { box.innerHTML = `<p class="mv-err">${esc(e.message)}</p>`; return; }
    // 3D notes name the geometry's real source (PubChem, the file), as the hover card does; ΔE to 2 dp as the table
    const c = S.view === '3d' ? currentConf(cs) : null;
    secs = secs.map(sec => ({ title: sec.title, rows: sec.rows.map(rw => {
      const x = geometrySource(rw, c);
      return x.key === 'energy' && c && c.source === 'OpenChemLib' && hasEnergy(cs, c) ? Object.assign({}, x, { value: energyText(c) }) : x;
    }) }));
    const html = secs.map(sec => {
      let body;
      if (sec.title === 'Functional groups' && r.groups.length) {
        body = sec.rows.map((rw, t) => rowHTML(rw, { click: true, on: S.group && S.group.id === 'g' + t, attrs: ` data-act="group" data-g="${t}" title="highlight this group"` })).join('');
      } else body = sec.rows.map(rw => rowHTML(rw)).join('');
      if (sec.title === 'Stereochemistry' && (r.stereo.atoms.length || r.stereo.bonds.length)) {
        const chips = r.stereo.atoms.map(s => `<button class="mv-chip${S.group && S.group.id === 'a' + s.index ? ' mv-on' : ''}" data-act="stereoA" data-i="${s.index}">${esc(atomName(r, s.index))} ${esc(s.label === '?' ? '(not specified)' : s.label === 'r' || s.label === 's' ? '(' + s.label + ', pseudoasymmetric)' : '(' + s.label + ')')}</button>`)
          .concat(r.stereo.bonds.map(s => { const b = r.bonds[s.index]; return `<button class="mv-chip${S.group && S.group.id === 'b' + s.index ? ' mv-on' : ''}" data-act="stereoB" data-k="${s.index}">${esc(atomName(r, b.a))}=${esc(atomName(r, b.b))} ${esc(s.label === 'unspecified' ? '(E/Z not specified)' : '(' + s.label + ')')}</button>`; }));
        body += `<div class="mv-chips" title="click to highlight">${chips.join('')}</div>`;
      }
      // a section the user collapsed stays collapsed through re-renders (chips, view switches, conformer updates)
      return `<details class="mv-sec"${S.closed.has(sec.title) ? '' : ' open'}><summary>${esc(sec.title)}</summary>${body}</details>`;
    }).join('');
    box.innerHTML = html;
  }
  function renderConfs() {
    const box = $('molConfs'), cs = S.cs;
    if (!cs || S.view !== '3d' || !cs.b3d || !cs.b3d.confs.length) { box.innerHTML = ''; return; }
    const B = cs.b3d, oc = B.confs.filter(c => c.source === 'OpenChemLib');
    let h = `<div class="mv-sec-h">Conformers <small>${B.done ? '' : 'searching…'}</small></div><table class="mv-tbl mv-conf-tbl"><thead><tr><th>#</th><th title="MMFF94 steric energy relative to the lowest conformer found">ΔE kcal/mol</th><th>MMFF94 energy</th><th>source</th></tr></thead><tbody>`;
    let single = '';
    for (const c of B.confs) {
      const on = c.key === B.cur ? ' class="mv-on"' : '', one = singleAtomNote(cs, c);
      if (one) single = one;
      if (c.source === 'OpenChemLib') h += `<tr data-act="conf" data-key="${esc(c.key)}"${on}><td>${c.rank}</td><td>${hasEnergy(cs, c) && Number.isFinite(c.relEnergy) ? num(c.relEnergy, 2) : '—'}</td><td${one ? ` title="${esc(one)}"` : ''}>${one ? 'none (single atom)' : hasEnergy(cs, c) ? num(c.energy, 2) : 'not optimised'}</td><td>OpenChemLib</td></tr>`;
      else if (c.source === 'file') h += `<tr data-act="conf" data-key="${esc(c.key)}"${on}><td>F</td><td>—</td><td title="the file's coordinates are shown as given; no energy is computed for them">not computed</td><td>the file (as given)</td></tr>`;
      else h += `<tr data-act="conf" data-key="${esc(c.key)}"${on}><td>P</td><td>—</td><td title="PubChem's energies are MMFF94s without electrostatics and not comparable">not comparable</td><td>PubChem CID ${esc(c.cid)}</td></tr>`;
    }
    h += '</tbody></table>';
    if (single) {                                   // one atom: no energy, nothing to search
      box.innerHTML = h + `<p class="mv-hint">${esc(sentence(single))} There is no energy to compare and no conformer search.</p>`;
      return;
    }
    const s = B.search;
    const stop = s && B.done ? { time: 'stopped at the time limit', exhausted: 'no further distinct torsion combinations', k: `all ${s.target} tries done`, cancelled: 'search stopped', untypable: 'MMFF94 could not type an atom, so no energies to compare', 'no-search': 'the search could not start for this structure' }[s.stoppedBy] : '';
    h += `<p class="mv-hint">${esc(MI.TEXT.mmffEnergy)}. ${oc.length} distinct conformer${oc.length === 1 ? '' : 's'} kept` +
      (s ? ` (${s.tried} tried${s.notConverged ? ', ' + s.notConverged + ' did not converge' : ''}${s.stereoDropped ? ', ' + s.stereoDropped + ' dropped for wrong stereo' : ''}${s.otherConfiguration ? ', ' + s.otherConfiguration + ' with another configuration at open stereo' : ''})` : '') +
      (stop ? '; ' + esc(stop) : '') + '. "Lowest" means the lowest found, not necessarily the global minimum. No Boltzmann populations are given.</p>';
    box.innerHTML = h;
  }
  function renderMeasures() {
    const box = $('molMeas'), cs = S.cs;
    if (!cs || S.view !== '3d' || (!S.measures.length && !S.measureN)) { box.innerHTML = ''; return; }
    const c = currentConf(cs);
    let h = `<div class="mv-sec-h">Measurements ${S.measures.length ? '<button class="mv-link" data-act="measclear">clear all</button>' : ''}</div>`;
    if (!S.measures.length) h += `<p class="mv-hint">Click ${S.measureN} atoms in the 3D model.</p>`;
    h += S.measures.map((m, t) => { const x = measureText(m, c); return `<div class="mv-row"><span class="mv-k">${esc(x.names)}</span><span class="mv-v">${esc(x.value)}${tagHTML('3d')} <button class="mv-x" data-act="measdel" data-t="${t}" title="remove">×</button></span></div>`; }).join('');
    if (S.measures.length) h += `<p class="mv-hint">${esc(c && c.source === 'PubChem' ? PUBCHEM_GEOMETRY : c && c.source === 'file' ? FILE_GEOMETRY : MI.TEXT.geometry3d)}; dihedral sign and names: IUPAC (Klyne–Prelog).</p>`;
    box.innerHTML = h;
  }

  /* atoms and bonds tables (sortable; a click selects) */
  // [key, header, header tooltip]: short headers so the tables fit the panel; the most used columns first
  const ATOM_COLS = [
    ['i', '#', 'atom number'], ['el', 'El.', 'element'], ['hyb', 'Hybrid.', 'hybridization'], ['cip', 'CIP', 'CIP label (R/S, r/s pseudoasymmetric, ? not specified)'],
    ['charge', 'Chg', 'formal charge (bookkeeping, for the Lewis structure as drawn)'], ['q', 'q (G–M)', 'Gasteiger–Marsili partial charge (empirical model), in e'],
    ['h', 'H', 'hydrogens on the atom'], ['lp', 'LP', 'lone pairs (for the Lewis structure as drawn)'], ['ox', 'Ox.', 'oxidation state'], ['arom', 'Ar.', 'aromatic'],
  ];
  const BOND_COLS = [['k', '#', 'bond number'], ['atoms', 'Bond', 'atoms joined'], ['type', 'Type', 'bond type'], ['len', 'Length Å', 'length on the 3D conformer shown'],
    ['pol', 'Polarity', 'polarity class from the electronegativity difference'], ['rot', 'Rot.', 'rotatable'], ['conj', 'Conj.', 'conjugated'], ['ring', 'Ring', 'in a ring']];
  function atomRowsTable(r) {
    return r.atoms.map((a, i) => ({
      i, key: { i, el: a.el, charge: a.charge, q: a.gasteiger ?? -99, hyb: a.hybridization.label || '', h: a.hCount, lp: a.lonePairs ?? -1, ox: a.oxidationState.value ?? -99, cip: a.cip || '', arom: a.aromatic ? 1 : 0 },
      cells: [atomName(r, i), a.el, a.isH ? '—' : a.hybridization.label ? MI.util.hybText(a.hybridization.label) : '—', a.cip || '',
        a.charge ? signed(a.charge, 0) : '0', a.gasteiger === null ? '—' : signed(a.gasteiger, 3), a.isH ? '' : String(a.hCount),
        a.lonePairs === null ? '—' : String(a.lonePairs), a.oxidationState.text, a.aromatic ? '✓' : ''],
    }));
  }
  function bondRowsTable(r, geo, fl) {
    return r.bonds.map((b, k) => {
      const len = geo && geo.bondLength && (geo.optimised || (geo.flagged && geo.flagged.has(k))) ? geo.bondLength[k] : null;
      const mark = geo && geo.flagged && geo.flagged.has(k) ? ' (flagged)' : fl && (fl.has(b.a) || fl.has(b.b)) ? ' (approx.)' : '';
      return {
        i: k, key: { k, atoms: b.a * 10000 + b.b, type: b.aromatic ? 1.5 : b.order, len: len ?? -1, pol: b.dEN ?? -1, rot: b.rotatable.yes ? 1 : 0, conj: b.conjugated ? 1 : 0, ring: b.inRing ? 1 : 0 },
        cells: [String(k + 1), atomName(r, b.a) + '–' + atomName(r, b.b), b.aromatic ? 'aromatic' : b.typeText, len !== null && Number.isFinite(len) ? len.toFixed(3) + mark : '—',
          b.polarity ? b.polarity.cls : '—', b.rotatable.yes ? '✓' : '', b.conjugated ? '✓' : '', b.inRing ? '✓' : ''],
      };
    });
  }
  function tableHTML(kind, cols, rows) {
    const st = S.sort[kind];
    rows.sort((x, y) => { const a = x.key[st.key], b = y.key[st.key]; return (a < b ? -1 : a > b ? 1 : x.i - y.i) * st.dir; });
    const selType = kind === 'atoms' ? 'atom' : 'bond';
    return `<table class="mv-tbl"><thead><tr>${cols.map(([k, t, tip]) => `<th data-act="sort" data-tbl="${kind}" data-key="${k}" class="${st.key === k ? 'mv-sorted' : ''}" title="${esc(tip ? tip + ' (click to sort)' : 'sort')}">${esc(t)}${st.key === k ? (st.dir > 0 ? ' ▲' : ' ▼') : ''}</th>`).join('')}</tr></thead><tbody>` +
      rows.map(rw => `<tr data-act="pick" data-type="${selType}" data-i="${rw.i}"${S.sel && S.sel.type === selType && S.sel.index === rw.i ? ' class="mv-on"' : ''}>${rw.cells.map(c => `<td>${esc(c)}</td>`).join('')}</tr>`).join('') + '</tbody></table>';
  }
  function renderTables() {
    const cs = S.cs;
    if (!cs || !cs.record) { $('molAtoms').innerHTML = ''; $('molBonds').innerHTML = ''; return; }
    const r = cs.record;
    if ($('molAtomsSec').open) $('molAtoms').innerHTML = tableHTML('atoms', ATOM_COLS, atomRowsTable(r)) +
      `<p class="mv-hint">Chg: formal charge; q (G–M): ${esc(MI.TEXT.gasteiger)}; LP: lone pairs; Ox.: oxidation state; Ar.: aromatic. Formal charge, lone pairs and oxidation state are ${esc(MI.TEXT.asDrawn)}. Click a row to select the atom.</p>`;
    const is3 = S.view === '3d';
    if ($('molBondsSec').open) $('molBonds').innerHTML = tableHTML('bonds', BOND_COLS, bondRowsTable(r, is3 ? cs.geo : null, is3 ? flaggedAtoms(cs, currentConf(cs)) : null)) +
      `<p class="mv-hint">Length: measured on the 3D conformer shown (3D view only; "flagged" = MMFF94 parameter gap, "approx." = in a charged or radical system MMFF94 cannot type, see the banner). Polarity: ${esc(MI.TEXT.polarityCutoffs)}. Rot.: rotatable, ${esc(MI.TEXT.rotatable)}. Conj.: conjugated.</p>`;
  }
  function renderMethods() {
    const box = $('molMethods'), cs = S.cs, MD = window.MolData;
    const pub = cs && cs.b3d && cs.b3d.pub;
    const li = (t, u) => `<li>${esc(t)}${u ? ` <a href="${esc(u)}" target="_blank" rel="noopener">${esc(u.replace(/^https?:\/\//, ''))}</a>` : ''}</li>`;
    let h = '<ul class="mv-src">';
    h += li('Structure handling, SMILES, InChI, descriptors (TPSA, Crippen logP/MR, Lipinski counts), CIP labels: RDKit.js 2026.03.6, BSD-3-Clause licence.', 'https://www.rdkit.org');
    h += li("Names, CAS numbers and the molecule library: the app's own resolver and namer (js/chem.js); names are shown only when they describe exactly this structure.");
    h += li('3D models: OpenChemLib 9.25.0 (BSD-3-Clause) conformer generator (torsions and bond lengths from COD crystal statistics), minimised with the MMFF94 force field (T. A. Halgren, J. Comput. Chem. 1996, 17, 490–519). Energies are MMFF94 steric energies, comparable only between conformers of one molecule.', 'https://github.com/cheminfo/openchemlib-js');
    h += li('3D rendering: three.js r185, MIT licence. Marching-cubes tables from three.js (Paul Bourke, Cory Gene Bloyd).', 'https://threejs.org');
    h += li('Partial charges: Gasteiger–Marsili PEOE (J. Gasteiger, M. Marsili, Tetrahedron 1980, 36, 3219–3228), ported from RDKit and checked against it; an empirical model.');
    h += li('Dipole moment: ' + String(MI.TEXT.dipole || '').trim());
    h += li('Stereochemistry the input leaves open: the values computed from the structure alone hold for every stereoisomer. The configuration of the 3D model shown at each open stereocentre and double bond is read back from its coordinates by RDKit (new CIP labeller) and named in the banner (R/S, r/s, E/Z), with cis or trans for open ring fusions and ring substituents where the geometry of the model decides it clearly (the dihedral across a ring fusion; the side of the ring each substituent is on). An open stereo axis (an allene, a hindered biaryl, a spiro compound, an alkylidene ring) is named in the banner, which says the model shows one of its forms but not which one: RDKit gives an axis no CIP label.');
    if (cs && cs.fileConf && cs.fileConf.coords) h += li('3D conformer "File": the coordinates of the opened MOL/SDF file, as given (not minimised, no energy); its atoms are matched to the structure read from the same file and checked like every 3D model.');
    h += li('Surface potential: ' + MI.TEXT.esp + '; V = 332.0637 Σ qᵢ/rᵢ kcal/mol per e.');
    // the sources behind the app's own rules (js/molinfo.js), each checked by tools/test_molinfo.js
    h += li('Molecule statistics (formula, counts, rings, degree of unsaturation): computed from the structure and checked field by field against RDKit (Python) on reference molecules.');
    h += li(`Functional groups: ${MI.TEXT.functionalGroups || "the app's own pattern table"}.`);
    h += li('Chirality: stereocentres and their R/S and E/Z labels from RDKit (new CIP labeller); a molecule is called chiral when it is not superimposable on its mirror image, tested on its canonical SMILES and InChI; meso as IUPAC defines it (an achiral member of a set of stereoisomers that also has chiral members).');
    h += li("Hybridization and conjugation: ports of RDKit's rules, checked against RDKit (Python) on 21,869 molecules. Lone pairs, VSEPR steric number and shape, σ/π counts: textbook conventions for the Lewis structure as drawn, checked against 489 hand-worked values.");
    h += li('Oxidation states: the electrons of each bond given to the more electronegative atom (Pauling electronegativities; aromatic bonds averaged over the Kekulé structures), for the Lewis structure as drawn. A value is flagged as borderline where a bond it depends on joins atoms within 0.05 in Pauling electronegativity, or where Allen electronegativities, which the IUPAC 2016 definition uses (P. Karen, P. McArdle, J. Takats, Pure Appl. Chem. 2016, 88, 831–839), order the two atoms the other way.');
    h += li('Bond polarity classes: ' + dot(MI.TEXT.polarityCutoffs) + ' Torsion names: Klyne–Prelog, as in the IUPAC Gold Book ("torsion angle").');
    h += '</ul>';
    if (MD && MD.SOURCES) {
      h += '<div class="mv-sec-h">Data tables (js/mol_data.js)</div><ul class="mv-src">' +
        Object.entries(MD.SOURCES).map(([k, v]) => `<li><b>${esc(k)}</b>: ${esc(v)}</li>`).join('') + '</ul>';
    }
    if (pub && window.PubChem3D) {
      h += '<div class="mv-sec-h">PubChem (3D coordinates shown)</div><ul class="mv-src">' +
        window.PubChem3D.CREDITS.map(c => li(c.text, c.url)).join('') +
        li(window.PubChem3D.TEXT.caveat) + '</ul>';
    }
    box.innerHTML = h;
  }
  function renderAll() {
    renderBanners(); renderFoot(); renderHead(); renderSelected(); renderSections(); renderConfs(); renderMeasures(); renderTables(); renderMethods(); renderControls();
  }

  /* ================================================================ controls */
  function exportOptions() {
    const cs = S.cs, has3D = !!(cs && currentConf(cs));
    if (S.view === '2d') return [['png2d', 'PNG image'], ['svg', 'SVG image'], ['mol2d', 'MOL file (2D)']];
    return [['png3d', 'PNG image', has3D], ['mol3d', 'MOL file (3D, with H)', has3D], ['xyz', 'XYZ file', has3D], ['sdf', 'SDF (all conformers)', !!(cs && cs.b3d && cs.b3d.confs.some(c => c.source === 'OpenChemLib'))]];
  }
  function renderControls() {
    const cs = S.cs, is3 = S.view === '3d', B = cs && cs.b3d;
    $('molView2d').setAttribute('aria-pressed', String(!is3)); $('molView3d').setAttribute('aria-pressed', String(is3));
    page.classList.toggle('mv-is3d', is3);
    syncControls();
    // conformer select
    const sel = $('molConf');
    if (B && B.confs.length) {
      sel.innerHTML = B.confs.map(c => `<option value="${esc(c.key)}">${c.source === 'PubChem' ? 'PubChem CID ' + esc(c.cid) : c.source === 'file' ? 'File (as given)'
        : c.rank + (singleAtomNote(cs, c) ? ' · single atom' : hasEnergy(cs, c) && Number.isFinite(c.relEnergy) ? (c.rank === 1 ? ' · lowest found' : ' · ΔE ' + num(c.relEnergy, 2) + ' kcal/mol') : ' · not optimised')}</option>`).join('') +
        (B.done ? '' : '<option disabled>searching…</option>');
      sel.value = B.cur || '';
      sel.disabled = false;
    } else { sel.innerHTML = `<option>${B && !B.done ? 'building…' : '—'}</option>`; sel.disabled = true; }
    $('molMore').disabled = !(B && B.done && B.result && B.result.ok && B.confs.some(c => c.source === 'OpenChemLib' && hasEnergy(cs, c)));
    const P = window.PubChem3D, pb = $('molPubchem'), busy = P ? P.busyFor() : 0;
    pb.disabled = !cs || !cs.record || !P || (B && B.pubBusy) || busy > 0 || !!(B && B.pub);
    pb.title = !P ? 'PubChem3D is not loaded' : B && B.pub ? 'PubChem\'s 3D model is loaded (see the conformer list)' : busy > 0 ? P.TEXT.busy : cs && cs.pub && !cs.pub.ok ? cs.pub.reason : P.TEXT.hint;
    pb.textContent = B && B.pubBusy ? 'Asking PubChem…' : P ? P.TEXT.button : 'Load PubChem 3D';
    if (busy > 0 && !renderControls.timer) renderControls.timer = setTimeout(() => { renderControls.timer = null; renderControls(); }, busy + 50);
    // export
    const fmt = $('molExportFmt'), prev = fmt.value, opts = exportOptions();
    if (opts.every(o => o[2] === false)) {         // nothing to export yet: say why instead of an empty box
      fmt.innerHTML = `<option value="" disabled selected>${B && !B.done ? '3D model building…' : B && B.error ? 'no 3D model' : 'no 3D model yet'}</option>`;
    } else {
      fmt.innerHTML = opts.map(([v, t, ok]) => `<option value="${v}"${ok === false ? ' disabled' : ''}>${esc(t)}</option>`).join('');
      if (opts.some(o => o[0] === prev && o[2] !== false)) fmt.value = prev;
      else fmt.value = opts.find(o => o[2] !== false)[0];
    }
    $('molExport').disabled = !cs || !cs.record || opts.every(o => o[2] === false);
    $('molReset').disabled = !cs || !cs.record;
  }
  function syncControls() {
    const is3 = S.view === '3d';
    const sh = $('molShowH'), on = is3 ? S.showH3d : S.showH2d;
    sh.setAttribute('aria-pressed', String(on));
    $('molSpin').setAttribute('aria-pressed', String(S.spin));
    $('molStyle').value = S.style; $('molLabels').value = S.labels; $('molDraw').value = S.textbook2d ? 'textbook' : 'color'; $('molSurface').value = S.surface; $('molMeasure').value = String(S.measureN);
    const q = S.cs && S.cs.record ? potentialCharges(S.cs.record) : null;
    const po = $('molSurface').querySelector('option[value="potential"]');
    po.disabled = !!(S.cs && S.cs.record) && !q;
    po.title = po.disabled ? (S.cs.record.gasteigerNote || 'partial charges are not available') : MI.TEXT.esp;
  }
  function applyView() {
    const is3 = S.view === '3d';
    $('mol2d').hidden = is3;
    if (v3) $('mol3d').hidden = !is3;
    else $('mol3d').hidden = true;
    hideCard();
    if (!is3) setStatus('');
    else if (S.cs && S.cs.b3d) setStatus(S.cs.b3d.statusText || '', !S.cs.b3d.done);
    renderLegend();
    applyHighlights();
    renderAll();
  }
  function setView(view) {
    if (view === S.view) return;
    S.view = view;
    if (view === '2d') { S.picks = []; if (S.spin && v3) { v3.setSpin(false); S.spin = false; } }
    applyView();
    if (view === '3d') enter3D();
    else if (v2 && v2info && v2info.cs !== S.cs) render2D();
  }

  /* ================================================================ exports */
  const IMAGE_FMTS = new Set(['png2d', 'svg', 'png3d']);
  /* The Export button first asks colour or black and white (one of the two). Only images have a colour; MOL, XYZ and
     SDF files are text, so they go straight out. The choice is remembered for the next export. */
  function askExport() {
    const cs = S.cs;
    if (!cs || !cs.record) return;
    const fmt = $('molExportFmt').value;
    if (!IMAGE_FMTS.has(fmt)) { doExport(false, false); return; }
    const dlg = $('molExportDlg');
    if (!dlg || typeof dlg.showModal !== 'function') { doExport(false, false); return; }
    const pick = dlg.querySelector(`input[name="molExportMode"][value="${S.exportBW ? 'bw' : 'color'}"]`);
    if (pick) pick.checked = true;
    const bgPick = dlg.querySelector(`input[name="molExportBg"][value="${S.exportTransparent ? 'transparent' : 'white'}"]`);
    if (bgPick) bgPick.checked = true;
    $('molExportDlgNote').textContent = 'Both save the picture exactly as shown (view, style, labels, surface). ' +
      'Black and white keeps the same picture with the colours turned to shades of grey. ' +
      'A transparent background lets the picture sit on any page or slide colour.';
    dlg.showModal();
  }
  /* an image blob → the same image in shades of grey (ITU-R BT.709 luminance, the usual sRGB→grey weights) */
  async function greyPNG(blob) {
    const bmp = await createImageBitmap(blob);
    const cv = document.createElement('canvas'); cv.width = bmp.width; cv.height = bmp.height;
    const ctx = cv.getContext('2d'); ctx.drawImage(bmp, 0, 0); bmp.close && bmp.close();
    const im = ctx.getImageData(0, 0, cv.width, cv.height), d = im.data;
    for (let i = 0; i < d.length; i += 4) { const y = Math.round(0.2126 * d[i] + 0.7152 * d[i + 1] + 0.0722 * d[i + 2]); d[i] = d[i + 1] = d[i + 2] = y; }
    ctx.putImageData(im, 0, 0);
    return new Promise((res, rej) => cv.toBlob(b => (b ? res(b) : rej(new Error('could not encode the PNG'))), 'image/png'));
  }
  /* an SVG text → the same drawing wrapped in a desaturating filter (feColorMatrix is plain SVG 1.1: every viewer) */
  function greySVG(text) {
    const doc = new DOMParser().parseFromString(text, 'image/svg+xml'), root = doc.documentElement;
    if (root.nodeName !== 'svg') return text;
    const NS = 'http://www.w3.org/2000/svg';
    const defs = doc.createElementNS(NS, 'defs'), f = doc.createElementNS(NS, 'filter'), m = doc.createElementNS(NS, 'feColorMatrix');
    f.setAttribute('id', 'mv-grey'); m.setAttribute('type', 'saturate'); m.setAttribute('values', '0');
    f.appendChild(m); defs.appendChild(f);
    const g = doc.createElementNS(NS, 'g'); g.setAttribute('filter', 'url(#mv-grey)');
    while (root.firstChild) g.appendChild(root.firstChild);
    root.appendChild(defs); root.appendChild(g);
    return new XMLSerializer().serializeToString(doc);
  }
  async function doExport(bw, transparent) {
    const cs = S.cs;
    if (!cs || !cs.record) return;
    const fmt = $('molExportFmt').value, base = exportName(cs), c = currentConf(cs), suffix = bw ? '_bw' : '';
    try {
      // atom labels that print a model's values (partial charges) carry the model's name in the image
      const cap2 = S.labels === 'partial' ? { caption: 'Atom labels: ' + MI.TEXT.gasteiger }
        : S.labels === 'electrons' ? { caption: 'Dots: lone pairs, ' + MI.TEXT.asDrawn } : {};
      if (transparent) cap2.transparent = true;
      if (fmt === 'png2d') { const b = await v2.exportPNG(2, cap2); download(base + suffix + '.png', bw ? await greyPNG(b) : b); }
      else if (fmt === 'svg') { const t = v2.exportSVG(cap2); download(base + suffix + '.svg', new Blob([bw ? greySVG(t) : t], { type: 'image/svg+xml' })); }
      else if (fmt === 'mol2d') download(base + '_2d.mol', new Blob([titleOf(cs) + cs.mb2d.slice(cs.mb2d.indexOf('\n'))], { type: 'chemical/x-mdl-molfile' }));
      else if (fmt === 'png3d' && v3 && c) {
        // with the potential surface on, the image carries the legend (scale and method), as the page does: a red/blue
        // map without them would read as a quantum-chemical ESP map
        // A potential that is zero everywhere (Cl2) has no colour scale: its words go in as a caption, with no bar, as
        // on the page. (Labels the view had to leave out are listed in the image by Mol3DView itself.)
        const g = legendInfo(), o = transparent ? { transparent: true } : {}, cap = [];
        if (g && g.clamp > 0) o.legend = { title: g.title, minLabel: g.minLabel, midLabel: g.midLabel, maxLabel: g.maxLabel, note: dot(g.scale) + ' ' + dot(g.note) };
        else if (g) cap.push(dot(g.title) + ' ' + dot(g.scale) + ' ' + dot(g.note));
        // atom labels that print a model's values carry the model's name too
        if (S.labels === 'partial') cap.push('Atom labels: ' + MI.TEXT.gasteiger);
        if (S.labels === 'electrons') cap.push('Dots: lone pairs, ' + MI.TEXT.asDrawn);
        if (cap.length) o.caption = cap;
        // in black and white the red/blue potential map becomes shades of grey: the legend words still say which end
        // is which, and the bar's ends are labelled, so the picture stays readable
        const b = await v3.exportPNG(2, o);
        download(base + '_3d' + suffix + '.png', bw ? await greyPNG(b) : b);
      } else if ((fmt === 'mol3d' || fmt === 'xyz') && c) {
        const E = window.Mol3DEngine;
        const title = c.source === 'PubChem' ? cs.b3d.pub.attribution.exportTitle
          : c.source === 'file' ? `${titleOf(cs)} coordinates from the opened file (as given)`
            : `${titleOf(cs)} OpenChemLib conformer ${c.rank}` + (singleAtomNote(cs, c) ? ', ' + singleAtomNote(cs, c) : hasEnergy(cs, c) ? `, MMFF94 ${c.energy.toFixed(2)} kcal/mol` : ', not optimised');
        const text = fmt === 'mol3d' ? E.molfile(cs.record, c.coords, title) : E.xyz(cs.record, c.coords, title);
        download(base + (c.source === 'PubChem' ? '_pubchem' : c.source === 'file' ? '_file' : '_conf' + c.rank) + (fmt === 'mol3d' ? '.mol' : '.xyz'), new Blob([text], { type: 'text/plain' }));
      } else if (fmt === 'sdf') {
        const confs = cs.b3d.confs.filter(x => x.source === 'OpenChemLib');
        download(base + '_conformers.sdf', new Blob([window.Mol3DEngine.sdf(cs.record, confs, titleOf(cs))], { type: 'chemical/x-mdl-sdfile' }));
      }
    } catch (e) { addBanner3D(cs, 'error', 'Export failed: ' + esc(e.message || e)); renderBanners(); }
  }

  /* ================================================================ wiring */
  function init() {
    if (S.inited) return;
    S.inited = true;
    window.RDKitLoad().catch(() => {});             // start RDKit early: the first Show is then quick
    $('molEmpty').innerHTML = `<div><h2>Molecule viewer</h2><p>Paste or type a SMILES, a common name or an IUPAC name above (a CAS number, MOL text or a MOL/SDF file also work) — then hover over atoms and bonds for their values.</p>
      <p class="mv-examples">Try: ${['aspirin', 'caffeine', 'L-alanine', 'cholesterol', '(E)-but-2-ene', '2-methylpropan-1-ol', 'CC(=O)[O-].[Na+]'].map(x => `<button class="mv-chip" data-example="${esc(x)}">${esc(x)}</button>`).join(' ')}</p></div>`;
    renderAll();
  }

  const input = $('molInput');
  // Show submits exactly what is in the box; an empty box replaces any earlier error with the reason
  const submit = text => {
    if (String(text || '').trim()) showMolecule(text, { source: 'typed' });
    else showError('Type a SMILES, common name or IUPAC name first (or open a MOL/SDF file).');
  };
  // the name list stays as a dropdown, but nothing is completed or pre-selected while typing: short SMILES ("CO",
  // "CI", "NO") would otherwise become names ("CO2", "citral", "non-1,3-diene"). Enter shows what was typed; a
  // suggestion needs ArrowDown or a click.
  if (window.Sketcher && window.Sketcher.attachSuggest) window.Sketcher.attachSuggest(input, $('molSuggest'), submit, undefined, { autofill: false });
  $('molShow').addEventListener('click', () => submit(input.value));
  // a single-line box drops line breaks, so MOL text pasted into it could never be read: take it from the clipboard
  // as it is and read it as a MOL file (placeholder and error texts no longer promise typed MOL text)
  input.addEventListener('paste', e => {
    const t = e.clipboardData && e.clipboardData.getData('text');
    if (!t || !/[\r\n]/.test(t.trim()) || !C.looksLikeMolfile(t)) return;
    e.preventDefault();
    input.value = '';
    showMolecule(t, { source: 'pasted' }).then(r => { if (r.ok && S.info) { S.info.label = 'Read from the pasted MOL text'; renderResolved(); } });
  });
  $('molFile').addEventListener('change', () => {
    const f = $('molFile').files[0];
    if (!f) return;
    f.text().then(t => {
      $('molFile').value = ''; input.value = '';
      showMolecule(t, { source: 'file' }).then(r => { if (r.ok && S.info) { S.info.label = 'Read from the file ' + f.name; S.info.fileName = f.name; renderResolved(); renderFoot(); } });
    });
  });
  $('molComp').addEventListener('change', e => { activate(+e.target.value); });
  $('molView2d').addEventListener('click', () => setView('2d'));
  $('molView3d').addEventListener('click', () => setView('3d'));
  $('molStyle').addEventListener('change', e => { S.style = e.target.value; if (v3) v3.setStyle(S.style); });
  $('molShowH').addEventListener('click', () => {
    if (S.view === '3d') { S.showH3d = !S.showH3d; if (v3) v3.setShowH(S.showH3d); }
    else { S.showH2d = !S.showH2d; render2D(); renderBanners(); }
    hideCard(); syncControls();
  });
  $('molLabels').addEventListener('change', e => { S.labels = e.target.value; renderLabels(); });
  $('molDraw').addEventListener('change', e => { S.textbook2d = e.target.value === 'textbook'; if (v2) v2.setStyle({ textbook: S.textbook2d }); });
  $('molSpin').addEventListener('click', () => { S.spin = !S.spin; if (v3) v3.setSpin(S.spin); syncControls(); });
  $('molSurface').addEventListener('change', e => { S.surface = e.target.value; applySurface(); });
  $('molMeasure').addEventListener('change', e => { S.measureN = +e.target.value; S.picks = []; applyPickMode(); hideCard(); applyHighlights(); renderMeasures(); renderFoot(); });
  $('molReset').addEventListener('click', () => { if (S.view === '3d') { if (v3) v3.resetView(); } else if (v2) v2.fit(); });
  $('molConf').addEventListener('change', e => pickConformer(e.target.value));
  $('molMore').addEventListener('click', () => { const cs = S.cs; if (cs && cs.b3d && cs.b3d.result) startBuild(cs, cs.b3d.result); });
  $('molPubchem').addEventListener('click', loadPubChem);
  $('molExport').addEventListener('click', askExport);
  {
    const dlg = $('molExportDlg');
    if (dlg) {
      // the export runs from the button (or Enter) itself, not from the dialog's close event: that event is queued as
      // a task and does not arrive while the page is hidden. Escape and Cancel just close: nothing is exported.
      const go = e => {
        e.preventDefault();
        const bw = dlg.querySelector('input[name="molExportMode"]:checked');
        S.exportBW = !!bw && bw.value === 'bw';
        const bg = dlg.querySelector('input[name="molExportBg"]:checked');
        S.exportTransparent = !!bg && bg.value === 'transparent';
        dlg.close();
        doExport(S.exportBW, S.exportTransparent);
      };
      $('molExportGo').addEventListener('click', go);
      dlg.querySelector('form').addEventListener('submit', go);
      $('molExportCancel').addEventListener('click', () => dlg.close());
      dlg.addEventListener('click', e => { if (e.target === dlg) dlg.close(); });   // a click on the backdrop cancels
    }
  }
  $('molAtomsSec').addEventListener('toggle', renderTables);
  $('molBondsSec').addEventListener('toggle', renderTables);
  // remember which molecule sections the user collapsed ('toggle' does not bubble: capture phase)
  $('molSections').addEventListener('toggle', e => {
    const d = e.target;
    if (!d || d.tagName !== 'DETAILS' || !d.classList.contains('mv-sec')) return;
    const t = d.querySelector('summary') ? d.querySelector('summary').textContent : '';
    if (d.open) S.closed.delete(t); else S.closed.add(t);
  }, true);

  // clicks inside the page that act on data-act attributes (panel rows, chips, banners)
  page.addEventListener('click', e => {
    const ex = e.target.closest('[data-example]');
    if (ex) { input.value = ex.dataset.example; submit(ex.dataset.example); return; }
    const cp = e.target.closest('.mv-copy');
    if (cp && $('molSide').contains(cp)) {           // a name or identifier: copy it (clicking the label still opens the note)
      copyText(cp.dataset.copy, () => {
        if (cp.dataset.busy) return;
        cp.dataset.busy = '1';
        const was = cp.innerHTML;
        cp.textContent = 'Copied';
        setTimeout(() => { cp.innerHTML = was; delete cp.dataset.busy; }, 900);
      });
      return;
    }
    const t = e.target.closest('[data-act]');
    if (!t || !page.contains(t)) {
      const noteRow = e.target.closest('.mv-has-note');
      if (noteRow && $('molSide').contains(noteRow)) noteRow.classList.toggle('mv-open');
      return;
    }
    const cs = S.cs, r = cs && cs.record;
    switch (t.dataset.act) {
      case 'group': { const G = r.groups[+t.dataset.g]; if (G) toggleGroup(groupSpec(G.atoms), 'g' + t.dataset.g); break; }
      case 'stereoA': toggleGroup({ atoms: [+t.dataset.i], bonds: [] }, 'a' + t.dataset.i); break;
      case 'stereoB': { const k = +t.dataset.k, b = r.bonds[k]; toggleGroup({ atoms: [b.a, b.b], bonds: [k] }, 'b' + k); break; }
      case 'unpin': S.sel = null; applyHighlights(); renderSelected(); renderTables(); break;
      case 'pick': S.sel = { type: t.dataset.type, index: +t.dataset.i }; applyHighlights(); renderSelected(); renderTables(); break;
      case 'sort': { const st = S.sort[t.dataset.tbl]; if (st.key === t.dataset.key) st.dir = -st.dir; else { st.key = t.dataset.key; st.dir = 1; } renderTables(); break; }
      case 'conf': pickConformer(t.dataset.key); break;
      case 'measdel': S.measures.splice(+t.dataset.t, 1); applyMeasures(); renderMeasures(); break;
      case 'measclear': S.measures = []; S.picks = []; applyMeasures(); applyHighlights(); renderMeasures(); renderFoot(); break;
      case 'pubform': {
        const smi = t.dataset.smiles, cid = t.dataset.cid;
        input.value = smi;
        showMolecule(smi, { view: '3d', source: 'pubchem' }).then(res => {
          if (res.ok && S.info) { S.info.label = `Showing PubChem's form of CID ${cid} in place of the structure you gave (SMILES ${smi})`; renderResolved(); }
        });
        break;
      }
    }
  });
  hideCard();
  $('molStage').addEventListener('pointerleave', hideCard);

  // keyboard: only while this tab is shown, never while typing in a field
  document.addEventListener('keydown', e => {
    if (!page.classList.contains('active')) return;
    const tg = e.target, tag = tg && tg.tagName;
    if (tag === 'INPUT' || tag === 'SELECT' || tag === 'TEXTAREA' || (tg && tg.isContentEditable)) return;
    if (e.key === 'Escape') {
      if (S.picks.length) { S.picks = []; applyHighlights(); renderFoot(); return; }
      if (S.sel || S.group) { S.sel = null; S.group = null; applyHighlights(); renderSelected(); renderSections(); renderTables(); }
    }
  });

  // start on the first visit to the tab (app.js opens the page named in the URL without a click)
  const tab = document.querySelector('.tab[data-page="mol"]');
  if (tab) tab.addEventListener('click', init);
  if (location.hash === '#mol') init();

  window.MolView = {
    show: (textOrGraph, opts) => showMolecule(textOrGraph, opts || {}),
    depict2D,
    // for tests (molview_selftest.html, browser checks): the page state and the two views; not for page code
    _debug: { state: S, views: () => ({ v2, v3 }), hover: (hit, x, y) => onHover(hit, x, y, S.view),
      ringCisTrans, ringCandidates, axialOpen, legendInfo, measureText: (atoms, c) => measureText({ atoms }, c || currentConf(S.cs)) },
  };
})();
