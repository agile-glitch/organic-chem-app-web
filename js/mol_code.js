/* MolCode — the molecule code of the Molecule tab: one molecule with every atom label the tab shows, as an SDF record, and
   that record compressed into one line for pasting.

   The SDF record is an ordinary MOL block (the 2D drawing, as laid out) followed by data items, so RDKit, Open Babel and
   any SDF reader open it as the molecule and keep the items as properties:

     >  <ACM.format>            AgileCHEM molecule code 1
     >  <ACM.app>               the app version that wrote it ("1.0.9 (2026-10-13)")
     >  <ACM.atoms>             "13 8": the heavy atoms (the MOL block's atoms 1..13, in its order), then the H atoms 14..21,
                                one per implicit H, grouped by parent in ascending parent order (as RDKit add_hs numbers them)
     >  <ACM.labels.MODE>       one line per labelled atom: "<atom number> <label text as the tab shows it>"; an atom with no
                                label has no line, and an item with no lines says the mode labels nothing on this molecule
     >  <ACM.labels.angles>     one line per bond angle: "<centre> <atom> <atom> <text>", measured on the 3D model
     >  <ACM.note>              free text (e.g. why the bond angles are missing)

   MODE is each of MODES below (the Labels list without 'none', 'element' and 'index', which follow from the MOL block,
   and 'electrons', which draws the lone-pairs numbers as dots). Atom numbers are 1-based, as in the MOL block.

   The code for pasting is  ACM1:<base64url of the record compressed with raw DEFLATE>  ('ACM1' = AgileCHEM molecule,
   version 1). CompressionStream does the DEFLATE (every current browser, Node 18+).

     MolCode.build({molblock, title, app, nHeavy, nH, labels: {mode: string[]}, angles: [{a,i,j,text}] | null, notes: []}) -> SDF text
     MolCode.parse(text) -> {ok, molblock, app, nHeavy, nH, labels: {mode: string[nHeavy+nH]}, angles, notes, error}
     MolCode.hasData(text) -> true when the text is (or starts with) a record this file wrote
     MolCode.isCode(text) -> true for an 'ACM1:' code
     MolCode.encode(sdf) -> Promise<code>      MolCode.decode(code) -> Promise<SDF text>                                    */
(function (root) {
  'use strict';
  const PREFIX = 'ACM1:', FORMAT = 'AgileCHEM molecule code 1';
  const MODES = ['cip', 'hybridization', 'formal', 'partial', 'oxidation', 'lonepairs', 'pka', 'steric', 'electronic', 'taft'];

  const item = (name, lines) => `>  <${name}>\n${lines.filter(l => l !== '').join('\n')}${lines.some(l => l !== '') ? '\n' : ''}\n`;
  const oneLine = s => String(s).replace(/\s*[\r\n]+\s*/g, ' ').trim();   // a label never breaks the item it sits in

  function build(o) {
    let mb = String(o.molblock || '');
    if (!/M {2}END/.test(mb)) throw new Error('MolCode.build: the MOL block has no "M  END" line');
    mb = mb.slice(0, mb.indexOf('M  END') + 6) + '\n';
    if (o.title != null) mb = oneLine(o.title) + mb.slice(mb.indexOf('\n'));
    let out = mb + item('ACM.format', [FORMAT]);
    if (o.app) out += item('ACM.app', [oneLine(o.app)]);
    out += item('ACM.atoms', [`${o.nHeavy | 0} ${o.nH | 0}`]);
    const L = o.labels || {};
    for (const m of MODES) {
      if (!L[m]) continue;
      out += item('ACM.labels.' + m, L[m].map((t, i) => (t == null || oneLine(t) === '' ? '' : `${i + 1} ${oneLine(t)}`)));
    }
    if (o.angles) out += item('ACM.labels.angles', o.angles.map(x => `${x.a + 1} ${x.i + 1} ${x.j + 1} ${oneLine(x.text)}`));
    for (const n of o.notes || []) out += item('ACM.note', [oneLine(n)]);
    return out + '$$$$\n';
  }

  const hasData = text => /^>[^\n]*<ACM\.format>/m.test(String(text || ''));

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
    const get = n => items.find(x => x.name === n);
    const fmt = get('ACM.format');
    if (!fmt) return { ok: false, error: 'This SDF record has no molecule-code data (no ACM.format item).' };
    const version = /AgileCHEM molecule code (\d+)/.exec(fmt.lines.join(' '));
    if (!version) return { ok: false, error: `Unknown molecule-code format “${fmt.lines.join(' ')}”.` };
    if (+version[1] > 1) return { ok: false, error: `This molecule code is version ${version[1]}; this version of the app reads version 1. Update the app to read it.` };
    const atomsItem = get('ACM.atoms'), nums = atomsItem ? atomsItem.lines.join(' ').trim().split(/\s+/).map(Number) : [];
    if (nums.length !== 2 || !nums.every(Number.isInteger)) return { ok: false, error: 'The molecule code has no atom count (ACM.atoms).' };
    const [nHeavy, nH] = nums, N = nHeavy + nH;
    const labels = {};
    for (const m of MODES) {
      const it = get('ACM.labels.' + m);
      if (!it) continue;
      const arr = new Array(N).fill('');
      for (const line of it.lines) {
        const x = /^\s*(\d+)\s+(.*\S)\s*$/.exec(line);
        if (!x || +x[1] < 1 || +x[1] > N) return { ok: false, error: `The molecule code's ${m} labels name an atom that is not there (“${line}”).` };
        arr[+x[1] - 1] = x[2];
      }
      labels[m] = arr;
    }
    let angles = null;
    const an = get('ACM.labels.angles');
    if (an) {
      angles = [];
      for (const line of an.lines) {
        const x = /^\s*(\d+)\s+(\d+)\s+(\d+)\s+(.*\S)\s*$/.exec(line);
        if (!x || [x[1], x[2], x[3]].some(v => +v < 1 || +v > N)) return { ok: false, error: `The molecule code's bond angles name an atom that is not there (“${line}”).` };
        angles.push({ a: +x[1] - 1, i: +x[2] - 1, j: +x[3] - 1, text: x[4] });
      }
    }
    const app = get('ACM.app');
    return { ok: true, version: +version[1], molblock, app: app ? app.lines.join(' ') : '', nHeavy, nH, labels, angles,
      notes: items.filter(x => x.name === 'ACM.note').map(x => x.lines.join(' ')) };
  }

  const isCode = text => new RegExp('^' + PREFIX + '[A-Za-z0-9_-]+={0,2}$').test(String(text || '').trim());

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
    if (!isCode(t)) throw new Error('A molecule code starts with ' + PREFIX + ' followed by letters, digits, - and _.');
    if (typeof DecompressionStream === 'undefined') throw new Error('This browser cannot decompress (no DecompressionStream), so it cannot read molecule codes.');
    let bytes;
    try { bytes = await pipe(fromB64url(t.slice(PREFIX.length)), new DecompressionStream('deflate-raw')); }
    catch (e) { throw new Error('This molecule code is damaged or cut short (it could not be unpacked). Copy the whole code again.'); }
    const sdf = new TextDecoder().decode(bytes);
    if (!hasData(sdf)) throw new Error('This molecule code unpacks to something that is not a molecule record.');
    return sdf;
  }

  const api = { PREFIX, FORMAT, MODES, build, parse, hasData, isCode, encode, decode };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.MolCode = api;
})(typeof window !== 'undefined' ? window : globalThis);
