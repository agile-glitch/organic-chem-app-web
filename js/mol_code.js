/* MolCode — AGILES, the molecule code of the Molecule tab: one molecule with every atom label the tab shows, as an SDF
   record, and that record compressed into one line for pasting. (Written as ACM before 2026-10-15: ACM1: codes and ACM.*
   items are still read.)

   The SDF record is an ordinary MOL block (the 2D drawing, as laid out) followed by data items, so RDKit, Open Babel and
   any SDF reader open it as the molecule and keep the items as properties:

     >  <AGILES.format>         AGILES 1
     >  <AGILES.app>            the app version that wrote it ("1.1.0 (2026-10-15)")
     >  <AGILES.h_parents>      the parent of each H atom (nHeavy+1 … nHeavy+nH), as atom numbers, in one line
     >  <AGILES.name> … <AGILES.cas>   the names and identifiers (IDS below: common name, IUPAC name and whether it was
                                checked, SMILES, InChI, InChIKey, CAS number), one line each, an unknown one left out
     >  <AGILES.atoms>          "13 8": the heavy atoms (the MOL block's atoms 1..13, in its order), then the H atoms 14..21,
                                one per implicit H, grouped by parent in ascending parent order (as RDKit add_hs numbers them)
     >  <AGILES.labels.MODE>    one line per labelled atom: "<atom number> <label text as the tab shows it>"; an atom with no
                                label has no line, and an item with no lines says the mode labels nothing on this molecule
     >  <AGILES.labels.angles>  one line per bond angle: "<centre> <atom> <atom> <text>", measured on the 3D model
     >  <AGILES.note>           free text (e.g. why the bond angles are missing)

   MODE is each of MODES below (the Labels list without 'none', 'element' and 'index', which follow from the MOL block,
   and 'electrons', which draws the lone-pairs numbers as dots). Atom numbers are 1-based, as in the MOL block.

   The code for pasting is  AGILES1:<base64url of the record compressed with raw DEFLATE>  (AGILES, version 1).
   CompressionStream does the DEFLATE (every current browser, Node 18+).

     MolCode.build({molblock, title, app, nHeavy, nH, ids: {name, iupac, …}, labels: {mode: string[]}, angles: [{a,i,j,text}] | null, notes: []}) -> SDF text
     MolCode.parse(text) -> {ok, molblock, app, nHeavy, nH, ids, labels: {mode: string[nHeavy+nH]}, angles, notes, error}
     MolCode.hasData(text) -> true when the text is (or starts with) a record this file wrote
     MolCode.isCode(text) -> true for an 'AGILES1:' code (or an older 'ACM1:' one)
     MolCode.encode(sdf) -> Promise<code>      MolCode.decode(code) -> Promise<SDF text>                                    */
(function (root) {
  'use strict';
  const PREFIX = 'AGILES1:', OLD_PREFIX = 'ACM1:', FORMAT = 'AGILES 1', NS = 'AGILES.', OLD_NS = 'ACM.';
  /* the molecule's names and identifiers, as the tab's Identity section shows them (each item one line; an unknown one is
     left out): name (common name), iupac (the app's namer; iupac_checked 'yes' when its name reads back to exactly this
     structure, 'no' when it could not be checked), smiles (the SMILES as the molecule was given: typed, or the library
     entry's for a name; else the app's own writer), smiles_rdkit (RDKit's canonical SMILES: how the tab reads the
     molecule, and the order its drawing numbers the atoms in), inchi, inchikey, cas (from the app's library), formula
     (Hill order, plain digits, charge at the end: C9H8O4, C2H3O2-), isotopic_formula (only for an isotope-labelled molecule,
     isotopes after their element, subscript digits, as the tab writes it: CHD₃O), mw (molecular weight as the tab shows it: "180.16 g/mol"),
     monoisotopic ("180.0423 Da"; for an isotope-labelled molecule the exact mass of that isotopologue) */
  const IDS = ['name', 'iupac', 'iupac_checked', 'smiles', 'smiles_rdkit', 'inchi', 'inchikey', 'cas', 'formula', 'isotopic_formula', 'mw', 'monoisotopic'];
  const MODES = ['cip', 'hybridization', 'formal', 'partial', 'oxidation', 'lonepairs', 'pka', 'steric', 'electronic', 'taft'];

  const item = (name, lines) => `>  <${name}>\n${lines.filter(l => l !== '').join('\n')}${lines.some(l => l !== '') ? '\n' : ''}\n`;
  const oneLine = s => String(s).replace(/\s*[\r\n]+\s*/g, ' ').trim();   // a label never breaks the item it sits in

  function build(o) {
    let mb = String(o.molblock || '');
    if (!/M {2}END/.test(mb)) throw new Error('MolCode.build: the MOL block has no "M  END" line');
    mb = mb.slice(0, mb.indexOf('M  END') + 6) + '\n';
    if (o.title != null) mb = oneLine(o.title) + mb.slice(mb.indexOf('\n'));
    let out = mb + item(NS + 'format', [FORMAT]);
    if (o.app) out += item(NS + 'app', [oneLine(o.app)]);
    out += item(NS + 'atoms', [`${o.nHeavy | 0} ${o.nH | 0}`]);
    if (o.hParents && o.hParents.length) out += item(NS + 'h_parents', [o.hParents.map(x => x + 1).join(' ')]);
    const I = o.ids || {};
    for (const k of IDS) if (I[k] != null && oneLine(I[k]) !== '') out += item(NS + k, [oneLine(I[k])]);
    const L = o.labels || {};
    for (const m of MODES) {
      if (!L[m]) continue;
      out += item(NS + 'labels.' + m, L[m].map((t, i) => (t == null || oneLine(t) === '' ? '' : `${i + 1} ${oneLine(t)}`)));
    }
    if (o.angles) out += item(NS + 'labels.angles', o.angles.map(x => `${x.a + 1} ${x.i + 1} ${x.j + 1} ${oneLine(x.text)}`));
    for (const n of o.notes || []) out += item(NS + 'note', [oneLine(n)]);
    return out + '$$$$\n';
  }

  const hasData = text => /^>[^\n]*<(AGILES|ACM)\.format>/m.test(String(text || ''));

  function parse(text) {
    const t = String(text || '').replace(/\r\n?/g, '\n');
    const end = t.indexOf('M  END');
    if (end < 0) return { ok: false, error: 'This is not an SDF record (it has no "M  END" line).' };
    const molblock = t.slice(0, end + 6) + '\n';
    const rest = t.slice(end + 6).split(/\n\$\$\$\$/)[0];
    const items = [];
    let cur = null;
    for (const line of rest.split('\n')) {
      const h = /^>.*?<([^>]+)>/.exec(line);
      if (h) { cur = { name: h[1], lines: [] }; items.push(cur); continue; }
      if (!cur) continue;
      if (line.trim() === '') { cur = null; continue; }
      cur.lines.push(line);
    }
    const get = n => items.find(x => x.name === NS + n) || items.find(x => x.name === OLD_NS + n);   // AGILES.*, or ACM.* before the rename
    const fmt = get('format');
    if (!fmt) return { ok: false, error: 'This SDF record has no AGILES data (no AGILES.format item).' };
    const version = /(?:AGILES|AgileCHEM molecule code) (\d+)/.exec(fmt.lines.join(' '));
    if (!version) return { ok: false, error: `Unknown AGILES format “${fmt.lines.join(' ')}”.` };
    if (+version[1] > 1) return { ok: false, error: `This AGILES code is version ${version[1]}; this version of the app reads version 1. Update the app to read it.` };
    const atomsItem = get('atoms'), nums = atomsItem ? atomsItem.lines.join(' ').trim().split(/\s+/).map(Number) : [];
    if (nums.length !== 2 || !nums.every(Number.isInteger)) return { ok: false, error: 'The AGILES code has no atom count (AGILES.atoms).' };
    const [nHeavy, nH] = nums, N = nHeavy + nH;
    const labels = {};
    for (const m of MODES) {
      const it = get('labels.' + m);
      if (!it) continue;
      const arr = new Array(N).fill('');
      for (const line of it.lines) {
        const x = /^\s*(\d+)\s+(.*\S)\s*$/.exec(line);
        if (!x || +x[1] < 1 || +x[1] > N) return { ok: false, error: `The AGILES code's ${m} labels name an atom that is not there (“${line}”).` };
        arr[+x[1] - 1] = x[2];
      }
      labels[m] = arr;
    }
    let angles = null;
    const an = get('labels.angles');
    if (an) {
      angles = [];
      for (const line of an.lines) {
        const x = /^\s*(\d+)\s+(\d+)\s+(\d+)\s+(.*\S)\s*$/.exec(line);
        if (!x || [x[1], x[2], x[3]].some(v => +v < 1 || +v > N)) return { ok: false, error: `The AGILES code's bond angles name an atom that is not there (“${line}”).` };
        angles.push({ a: +x[1] - 1, i: +x[2] - 1, j: +x[3] - 1, text: x[4] });
      }
    }
    const hp = get('h_parents');
    const hParents = hp ? hp.lines.join(' ').trim().split(/\s+/).map(x => +x - 1) : null;
    const app = get('app'), ids = {};
    for (const k of IDS) { const it = get(k); if (it && it.lines.length) ids[k] = it.lines.join(' ').trim(); }
    return { ok: true, version: +version[1], molblock, app: app ? app.lines.join(' ') : '', nHeavy, nH, hParents, labels, angles, ids,
      notes: items.filter(x => x.name === NS + 'note' || x.name === OLD_NS + 'note').map(x => x.lines.join(' ')) };
  }

  const isCode = text => /^(AGILES1|ACM1):[A-Za-z0-9_-]+={0,2}$/.test(String(text || '').trim());

  async function pipe(bytes, stream) {
    const out = new Response(new Blob([bytes]).stream().pipeThrough(stream));
    return new Uint8Array(await out.arrayBuffer());
  }
  function toB64url(bytes) {
    let s = '';
    for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000));
    return btoa(s).replace(/\+/g, '-').replace(/\//g, '_').replace(/=+$/, '');
  }
  function fromB64url(s) {
    const b = atob(s.replace(/-/g, '+').replace(/_/g, '/') + '==='.slice((s.length + 3) % 4));
    const u = new Uint8Array(b.length);
    for (let i = 0; i < b.length; i++) u[i] = b.charCodeAt(i);
    return u;
  }
  async function encode(sdf) {
    if (typeof CompressionStream === 'undefined') throw new Error('This browser cannot compress (no CompressionStream); export the SDF file instead.');
    return PREFIX + toB64url(await pipe(new TextEncoder().encode(String(sdf)), new CompressionStream('deflate-raw')));
  }
  async function decode(code) {
    const t = String(code || '').trim();
    if (!isCode(t)) throw new Error('An AGILES code starts with ' + PREFIX + ' followed by letters, digits, - and _.');
    if (typeof DecompressionStream === 'undefined') throw new Error('This browser cannot decompress (no DecompressionStream), so it cannot read AGILES codes.');
    let bytes;
    try { bytes = await pipe(fromB64url(t.slice(t.indexOf(':') + 1)), new DecompressionStream('deflate-raw')); }
    catch (e) { throw new Error('This AGILES code is damaged or cut short (it could not be unpacked). Copy the whole code again.'); }
    const sdf = new TextDecoder().decode(bytes);
    if (!hasData(sdf)) throw new Error('This AGILES code unpacks to something that is not a molecule record.');
    return sdf;
  }

  const api = { PREFIX, OLD_PREFIX, FORMAT, MODES, IDS, build, parse, hasData, isCode, encode, decode };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.MolCode = api;
})(typeof window !== 'undefined' ? window : globalThis);
