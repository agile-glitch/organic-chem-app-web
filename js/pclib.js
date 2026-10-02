/* PCLib — the broad offline compound library: every PubChem compound with a CAS number (built by
   tools/pubchem_lib.py into data/pclib/), looked up by CAS number or by PubChem's main name (the title).

   PCLib.lookup(text) → Promise<{ok:true, matches:[{cid, title, smiles, cas:[…]}], via:'CAS'|'name'}
                                | {ok:false, kind:'not-found'|'not-installed'|'error', error}>
   PCLib.info() → Promise<index | null>   counts and build date, null when data/pclib/ is not there
   PCLib.resolve(text) → Promise<{ok, smiles, cid, title, note, saved?} | {ok:false, choices?, error}>
       what the search boxes use when the curated library (js/chem.js) does not know the text: first the lookups
       saved from PubChem in this browser, then this library. Several different compounds under one CAS number or
       name come back as `choices`; nothing is picked for the user.

   The data are pieces of gzip + base64 in <script> files (a double-clicked page cannot read binary files), read
   one at a time when a lookup needs them and kept for the session: startup is unaffected. Hashing (32-bit FNV-1a over
   the UTF-8 bytes of the key) and the name key (lower case, spaces collapsed) match the build script. */
window.PCLib = (() => {
  'use strict';
  const BASE = 'data/pclib/';
  const CAS_RE = /^\d{2,7}-\d\d-\d$/;
  const scripts = new Map(), pieces = new Map();
  let indexP = null;

  function script(src) {
    if (!scripts.has(src)) {
      scripts.set(src, new Promise((res, rej) => {
        const s = document.createElement('script');
        s.src = src;
        s.onload = () => { s.remove(); res(); };
        s.onerror = () => { s.remove(); scripts.delete(src); rej(new Error('could not load ' + src)); };
        document.head.appendChild(s);
      }));
    }
    return scripts.get(src);
  }
  function info() {
    if (!indexP) indexP = script(BASE + 'index.js').then(() => window.PCLIB_INDEX || null, () => { indexP = null; return null; });
    return indexP;
  }
  function fnv1a(key) {
    let h = 0x811C9DC5;
    for (const b of new TextEncoder().encode(key)) { h ^= b; h = Math.imul(h, 0x01000193) >>> 0; }
    return h;
  }
  const nameKey = s => String(s).trim().toLowerCase().replace(/\s+/g, ' ');
  async function gunzip(b64) {
    const bin = atob(b64), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    return new Response(new Blob([u8]).stream().pipeThrough(new DecompressionStream('gzip'))).text();
  }
  /* one piece's lines, loaded once per session */
  function piece(name) {
    if (!pieces.has(name)) {
      pieces.set(name, (async () => {
        await script(BASE + name + '.js');
        const b64 = window.PCLIB && window.PCLIB[name];
        if (b64 == null) throw new Error(name + ' held no data');
        delete window.PCLIB[name];                    // the lines below are what is kept
        return b64 ? (await gunzip(b64)).split('\n') : [];
      })());
      pieces.get(name).catch(() => pieces.delete(name));
    }
    return pieces.get(name);
  }
  async function rows(cids, idx) {
    const out = [];
    for (const cid of cids) {
      let lo = 0, hi = idx.cidStarts.length - 1, p = -1;      // the last piece that starts at or before this CID
      while (lo <= hi) { const m = (lo + hi) >> 1; if (idx.cidStarts[m] <= cid) { p = m; lo = m + 1; } else hi = m - 1; }
      if (p < 0) continue;
      const L = await piece('cid_' + String(p).padStart(4, '0'));
      let a = 0, b = L.length - 1;                              // rows are in CID order
      while (a <= b) {
        const m = (a + b) >> 1, c = parseInt(L[m], 10);
        if (c === cid) {
          const [, title, smiles, cas] = L[m].split('\t');
          out.push({ cid, title: title || `CID ${cid}`, smiles, cas: cas ? cas.split(';') : [] });
          break;
        }
        if (c < cid) a = m + 1; else b = m - 1;
      }
    }
    return out;
  }
  async function lookup(text) {
    const q = String(text || '').trim();
    if (!q) return { ok: false, kind: 'not-found', error: 'Nothing to look up.' };
    try {
      const idx = await info();
      if (!idx) return { ok: false, kind: 'not-installed', error: 'The broad offline library (data/pclib/) is not installed.' };
      const isCas = CAS_RE.test(q), key = isCas ? q : nameKey(q);
      const kind = isCas ? 'cas' : 'name', n = isCas ? idx.nCas : idx.nName;
      const L = await piece(kind + '_' + String(fnv1a(key) % n).padStart(3, '0'));
      const cids = [...new Set(L.filter(l => l.startsWith(key + '\t')).map(l => parseInt(l.slice(key.length + 1), 10)))];
      if (!cids.length) return { ok: false, kind: 'not-found', error: `"${q}" is not in the offline PubChem library either.` };
      const matches = await rows(cids.slice(0, 12), idx);
      if (!matches.length) return { ok: false, kind: 'error', error: 'The offline library\'s record could not be read.' };
      return { ok: true, matches, via: isCas ? 'CAS' : 'name', more: Math.max(0, cids.length - 12) };
    } catch (e) {
      return { ok: false, kind: 'error', error: 'The offline library could not be read: ' + (e && e.message || e) };
    }
  }

  /* the search boxes' fallback, in order: lookups saved from PubChem in this browser, then this library */
  async function resolve(text) {
    const P = window.PubChem3D;
    const saved = P && P.savedName ? P.savedName(text) : null;
    if (saved && saved.ok) return Object.assign({}, saved, { note: saved.note + ' (from your saved PubChem lookups)', saved: true });
    const r = await lookup(text);
    if (!r.ok) return r;
    const idx = await info();
    const src = `the offline copy of PubChem (${idx && idx.built ? 'built ' + idx.built : 'PubChem bulk data'})`;
    const casNote = r.via === 'CAS' ? '; CAS numbers in PubChem are synonyms supplied by depositors' : '';
    if (r.matches.length === 1) {
      const m = r.matches[0];
      return { ok: true, smiles: m.smiles, cid: m.cid, title: m.title,
        note: `From ${src}: CID ${m.cid}, ${m.title}${casNote}.` };
    }
    return { ok: false, kind: 'ambiguous', error: `The offline copy of PubChem has ${r.matches.length}${r.more ? '+' : ''} different records under "${text}". Choose one:`,
      choices: r.matches.map(m => ({ cid: m.cid, title: m.title, smiles: m.smiles, formula: '' })), src, casNote };
  }

  return { lookup, info, resolve };
})();
