/* PubChem3D — "Load PubChem 3D" on the Molecule tab: PubChem's computed 3D conformer of the molecule shown, used only
   when PubChem's record is provably the same structure, with its atoms put into the app's (record) numbering.

   PubChem3D.eligible(R, M) → {ok, reason, smiles}              checks made before any request
   PubChem3D.load(R, M, record?) → Promise<result>              never rejects; result.ok false carries `error`
   PubChem3D.mapSdf(R, M, record?, sdfText, {cid}?) → result    the same checks on an SDF already in hand (no network)
   PubChem3D.busyFor() → ms until requests are allowed again (after an HTTP 503), 0 when free
   PubChem3D.manufacturing(R, M) → Promise   how the compound is made: PubChem's "Methods of Manufacturing" (HSDB)
                                            entries with their cited sources, for the Retrosynthesis page
   PubChem3D.identify(R, smiles), PubChem3D.papersBoth(cidA, cidB)   papers mentioning two compounds together, kept
                                            in this browser (see below); shown by js/litref.js
   PubChem3D.TEXT, PubChem3D.CREDITS                            wording for the button, hint and Methods & sources
   result = {ok, cid, conformerId, coords: Float64Array(3N) in record order, sdf, energy, energyNote,
             identity: {status: 'identical'|'stereo-undefined', banner?, assigned: [..], mixedStereo},
             attribution: {text, url, policyText, policyUrl, exportTitle}}
          | {ok: false, kind, error, cid?, mismatch?: {kind, message, pubchemSmiles, pubchemSmilesNote?, url}}
   `record` is MolInfo's record (atoms 0..n-1 = M's atoms, then one H per implicit H grouped by parent, ascending);
   without one, that same order is built from M.

   Policy and facts, from PubChem's and NCBI's own pages: usage limits pubchem.ncbi.nlm.nih.gov/docs/programmatic-access
   (≤5 requests/s) and /docs/dynamic-request-throttling (HTTP 503), NCBI policies www.ncbi.nlm.nih.gov/home/about/policies/
   (≤3 requests/s, disclaimer), PUG-REST pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial (meant for use from web pages),
   PubChem3D coverage pubchem.ncbi.nlm.nih.gov/docs/pubchem3d and Bolton et al. 2011 (cited in CREDITS below); the
   CORS behaviour (no preflight support) was tested against the live service:
   - Requests only on an explicit click (the caller's job). One request in flight app-wide and ≥350 ms between request
     starts (PubChem asks for ≤5 requests/s, NCBI ≤3/s); 20 s timeout; HTTP 503 means throttled: no requests for 60 s
     and never an automatic retry. No custom headers and no credentials, so both requests are CORS "simple requests"
     (PubChem does not answer preflights). SMILES → CID and CID → SDF are cached for the session.
   - PubChem computes 3D only for one covalent unit of ≤50 heavy atoms of H C N O F Si P S Cl Br I (plus ≤15 rotatable
     bonds and MMFF94s atom types, which PubChem's own 404 reports), so other input is refused before any request.
   - Identity: FixedH InChI of the query vs the SDF. Standard InChIKeys miss tautomers and zwitterions (2-hydroxypyridine
     comes back as 2-pyridone, CID 8871). Any difference outside the stereo layers is refused, with the reason.
   - Stereo: every configuration the user specified must be reproduced (CIP labels perceived from the SDF's 3D
     coordinates by RDKit's new CIP labeller, get_stereo_tags); configurations the user left open are taken from
     PubChem's conformer and named in a banner (for undefined stereo PubChem mixes stereoisomers across conformers).
   - Final check: the mapped coordinates, read back by RDKit in record order (with the user's open stereo marked
     "either"), must give exactly the query's canonical SMILES.
   - Energies are PubChem's (MMFF94s with the Coulomb terms removed); nothing is recomputed here.
   Data: PubChem (NCBI/NLM), in the public domain; the NCBI disclaimer must be linked wherever it is shown. */
window.PubChem3D = (() => {
  'use strict';
  const API = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug';
  const MIN_GAP_MS = 350, TIMEOUT_MS = 20000, BUSY_MS = 60000, MAX_HEAVY = 50, MAX_MATCHES = 10000;
  const ALLOWED = new Set([1, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53]);   // PubChem3D's elements
  const SYMBOLS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr ' +
    'Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re ' +
    'Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh ' +
    'Fl Mc Lv Ts Og').split(' ');
  const sym = z => SYMBOLS[z - 1] || (z === 0 ? '*' : '#' + z);
  // "−4.70": a typographic minus, 2 decimals (as the app's own energies); PubChem's SDF gives 4
  const fmt = (v, dp) => { const r = +v.toFixed(dp); return r < 0 ? '−' + (-r).toFixed(dp) : Math.abs(r).toFixed(dp); };   // no "−0.00"
  const recordUrl = cid => `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}#section=3D-Conformer`;
  const POLICY_URL = 'https://www.ncbi.nlm.nih.gov/home/about/policies/';

  const TEXT = {
    button: 'Load PubChem 3D',
    hint: 'Sends this structure (as SMILES) to PubChem (NCBI) over the internet.',
    offline: 'Offline — PubChem 3D needs an internet connection. The app\'s own 3D model (OpenChemLib, MMFF94) is still shown.',
    unreachable: 'Could not reach PubChem (no internet connection, or the request was blocked). The app\'s own 3D model (OpenChemLib, MMFF94) is still shown.',
    notFound: 'This exact structure is not in PubChem, so there is no PubChem 3D model.',
    no3d: cid => `PubChem has no computed 3D model for CID ${cid}. It computes 3D only for single-component structures with ≤50 heavy atoms, ≤15 rotatable bonds and only H C N O F Si P S Cl Br I.`,
    busy: 'PubChem is busy or limiting requests — try again in a minute.',
    timeout: 'PubChem did not answer.',
    rejected: 'PubChem could not read this structure (it rejected the SMILES), so there is no PubChem 3D model.',
    energyMethod: 'MMFF94s without electrostatics (PubChem)',
    source: 'computed by PubChem3D (OMEGA, MMFF94s), not experimental',
    caveat: 'PubChem\'s 3D conformers are computed, not experimental, and are not energy minima; PubChem notes that its method aims at shapes seen in protein–ligand complexes, which can differ from the shape in the gas phase.',
  };
  const CREDITS = [
    { text: 'Kim S, Chen J, Cheng T, et al. PubChem 2025 update. Nucleic Acids Res. 2025;53(D1):D1516–D1525.', url: 'https://doi.org/10.1093/nar/gkae1059' },
    { text: 'Bolton EE, et al. PubChem3D: a new resource for scientists. J Cheminform. 2011;3:32.', url: 'https://doi.org/10.1186/1758-2946-3-32' },
    { text: 'NCBI policies and disclaimer (PubChem data: NCBI/NLM, public domain)', url: POLICY_URL },
  ];

  class Fail extends Error { constructor(kind, message, extra) { super(message); this.kind = kind; Object.assign(this, extra || {}); } }
  const failure = (kind, error, extra) => Object.assign({ ok: false, kind, error }, extra || {});
  const pairKey = (a, b) => a < b ? a + '-' + b : b + '-' + a;

  /* ---- RDKit helpers (every molecule made here is deleted) ---- */
  function molJson(m) {
    const doc = JSON.parse(m.get_json()), mol = doc.molecules[0] || {}, dA = (doc.defaults || {}).atom || {}, dB = (doc.defaults || {}).bond || {};
    const atoms = (mol.atoms || []).map(a => ({ z: a.z ?? dA.z ?? 6, chg: a.chg ?? dA.chg ?? 0, isotope: a.isotope ?? dA.isotope ?? 0,
      nRad: a.nRad ?? dA.nRad ?? 0, hs: a.impHs ?? dA.impHs ?? 0 }));
    const bonds = (mol.bonds || []).map(b => ({ a: b.atoms[0], b: b.atoms[1], bo: b.bo ?? dB.bo ?? 1 }));
    const conf = (mol.conformers || []).find(c => c.dim === 3) || null;
    return { atoms, bonds, coords: conf ? conf.coords : null };
  }
  function stereoTags(m) {   // new CIP labeller: {atoms: Map index → 'R'|'S'|'r'|'s'|'?', bonds: Map 'a-b' → 'E'|'Z'}
    const t = JSON.parse(m.get_stereo_tags() || '{}'), atoms = new Map(), bonds = new Map();
    for (const [i, lab] of t.CIP_atoms || []) atoms.set(i, String(lab).replace(/[()]/g, ''));
    for (const [a, b, lab] of t.CIP_bonds || []) bonds.set(pairKey(a, b), String(lab).replace(/[()]/g, ''));
    return { atoms, bonds };
  }
  function getMol(R, text, removeHs) {
    let m = null;
    try { m = R.get_mol(text, JSON.stringify({ removeHs })); } catch (e) { m = null; }
    if (m && ((m.is_valid && !m.is_valid()) || !m.get_num_atoms())) { m.delete(); m = null; }
    return m;
  }
  const DEFINED = new Set(['R', 'S', 'r', 's']);

  /* ---- before any request ---- */
  function eligible(R, M) {
    let j, smiles = '';
    try { j = molJson(M); smiles = M.get_smiles(); } catch (e) { return { ok: false, reason: 'There is no structure to look up.' }; }
    if (!j.atoms.length || !smiles) return { ok: false, reason: 'There is no structure to look up.' };
    if (smiles.includes('.')) return { ok: false, reason: 'PubChem computes 3D models only for a single molecule, not for salts or mixtures. Choose one component.' };
    const bad = [...new Set(j.atoms.filter(a => !ALLOWED.has(a.z)).map(a => sym(a.z)))];
    if (bad.length) return { ok: false, reason: `PubChem computes 3D models only for molecules made of H C N O F Si P S Cl Br I (this one contains ${bad.join(', ')}).` };
    const heavy = j.atoms.filter(a => a.z > 1).length;
    if (heavy > MAX_HEAVY) return { ok: false, reason: `PubChem computes 3D models only for molecules with at most ${MAX_HEAVY} heavy atoms (this one has ${heavy}).` };
    /* InChI, the identity check below, does not record unpaired electrons (·CH2–CH2· and ethene share one InChI) */
    if (j.atoms.some(a => a.nRad)) return { ok: false, reason: 'PubChem 3D is not offered for radicals: the identity check (InChI) does not record unpaired electrons, so PubChem\'s record could not be confirmed to be the same species.' };
    return { ok: true, reason: '', smiles };
  }

  /* ---- network: one request at a time, spaced, with a timeout; a 503 stops everything for a minute ---- */
  const net = { chain: Promise.resolve(), lastStart: 0, busyUntil: 0 };
  const cidCache = new Map(), sdfCache = new Map(), pending = new Map();
  const sleep = ms => new Promise(r => setTimeout(r, ms));
  const isOffline = () => typeof navigator !== 'undefined' && navigator.onLine === false;
  const busyFor = () => Math.max(0, net.busyUntil - Date.now());

  function request(url, init) {
    const run = async () => {
      if (busyFor()) throw new Fail('busy', TEXT.busy);
      for (let wait; (wait = net.lastStart + MIN_GAP_MS - Date.now()) > 0; ) await sleep(wait);   // timers may fire early
      const ac = new AbortController(), timer = setTimeout(() => ac.abort(), TIMEOUT_MS);
      try {
        const sent = fetch(url, Object.assign({ credentials: 'omit', signal: ac.signal }, init || {}));
        // the time is taken once fetch has started the request (not before): then the next start is ≥ MIN_GAP_MS after
        // this one even when the millisecond clock ticks in between (the scripted-fetch test measured 349 ms)
        net.lastStart = Date.now();
        const res = await sent;
        const text = await res.text();   // the timeout covers the body too; a 503 body may be HTML
        if (res.status === 503) { net.busyUntil = Date.now() + BUSY_MS; throw new Fail('busy', TEXT.busy); }
        if (res.status === 504) throw new Fail('timeout', TEXT.timeout);
        return { status: res.status, text };
      } catch (e) {
        if (e instanceof Fail) throw e;
        if (ac.signal.aborted) throw new Fail('timeout', TEXT.timeout);
        throw new Fail('offline', isOffline() ? TEXT.offline : TEXT.unreachable);
      } finally { clearTimeout(timer); }
    };
    const p = net.chain.then(run, run);
    net.chain = p.catch(() => {});
    return p;
  }
  const cached = v => { if (v instanceof Fail) throw v; return v; };

  async function getCid(smiles) {
    if (cidCache.has(smiles)) return cached(cidCache.get(smiles));
    const r = await request(API + '/compound/smiles/cids/JSON', { method: 'POST', body: new URLSearchParams({ smiles }) });
    if (r.status === 400) { const f = new Fail('rejected', TEXT.rejected); cidCache.set(smiles, f); throw f; }
    if (r.status === 404) { cidCache.set(smiles, 0); return 0; }
    if (r.status !== 200) throw new Fail('http', `PubChem answered with an error (HTTP ${r.status}).`);
    let list = null;
    try { list = JSON.parse(r.text).IdentifierList.CID; } catch (e) { list = null; }
    if (!Array.isArray(list) || !list.length || !list.every(c => Number.isInteger(c) && c >= 0)) throw new Fail('http', 'PubChem\'s answer could not be read.');
    if (list.length > 1) throw new Fail('ambiguous', `PubChem matched this structure to ${list.length} records (CID ${list.slice(0, 4).join(', ')}${list.length > 4 ? ', …' : ''}), so no single 3D model can be chosen.`);
    cidCache.set(smiles, list[0]);
    return list[0];   // 0 = not in PubChem
  }
  async function getSdf(cid) {
    if (sdfCache.has(cid)) return cached(sdfCache.get(cid));
    const r = await request(`${API}/compound/cid/${cid}/SDF?record_type=3d`);
    if (r.status === 404) { const f = new Fail('no-3d', TEXT.no3d(cid)); sdfCache.set(cid, f); throw f; }
    if (r.status !== 200) throw new Fail('http', `PubChem answered with an error (HTTP ${r.status}).`);
    if (!/M {2}END/.test(r.text)) throw new Fail('parse', 'PubChem\'s 3D file could not be read.');
    sdfCache.set(cid, r.text);
    return r.text;
  }
  function fetchRecord(smiles) {   // a second click while the first is under way shares its requests
    if (!pending.has(smiles)) {
      const p = (async () => { const cid = await getCid(smiles); return { cid, sdf: cid ? await getSdf(cid) : null }; })();
      pending.set(smiles, p);
      p.then(() => pending.delete(smiles), () => pending.delete(smiles));
    }
    return pending.get(smiles);
  }

  async function load(R, M, record) {
    const el = eligible(R, M);
    if (!el.ok) return failure('ineligible', el.reason);
    if (busyFor()) return failure('busy', TEXT.busy);
    if (isOffline()) return failure('offline', TEXT.offline);
    let cid = null;
    try {
      const got = await fetchRecord(el.smiles);
      cid = got.cid;
      if (!cid) return failure('not-found', TEXT.notFound, { cid: 0 });
      return mapSdf(R, M, record, got.sdf, { cid });
    } catch (e) {
      if (e instanceof Fail) return failure(e.kind, e.message, cid ? { cid } : {});
      return failure('error', 'PubChem 3D failed: ' + (e && e.message || e), cid ? { cid } : {});
    }
  }

  /* ---- SDF text ---- */
  function splitSdf(text) {   // first record: {block (molfile), tags: {NAME: 'value'}}
    const rec = String(text).replace(/\r\n?/g, '\n').split(/^\${4}\s*$/m)[0];
    const end = rec.search(/^M {2}END\s*$/m);
    if (end < 0) return null;
    const block = rec.slice(0, end) + 'M  END\n', tags = {};
    let name = null;
    for (const line of rec.slice(end).split('\n').slice(1)) {   // "> <NAME>" then value lines up to a blank line
      const m = /^>.*<([^>]+)>/.exec(line);
      if (m) { name = m[1]; tags[name] = ''; } else if (!line.trim()) name = null; else if (name) tags[name] += (tags[name] ? '\n' : '') + line.trim();
    }
    return { block, tags };
  }
  /* bond stereo flags off: RDKit then reads every configuration from the coordinates (PubChem writes "either" for
     double bonds it leaves undefined, although its conformer has one geometry) */
  function coordinateStereoOnly(block) {
    const L = block.split('\n');
    if (!L[3] || !/V2000/.test(L[3])) return null;
    const nA = parseInt(L[3].slice(0, 3), 10), nB = parseInt(L[3].slice(3, 6), 10);
    if (!(nA > 0) || !(nB >= 0) || L.length < 4 + nA + nB) return null;
    for (let k = 4 + nA; k < 4 + nA + nB; k++) if (L[k].length >= 12) L[k] = L[k].slice(0, 9) + '  0' + L[k].slice(12);
    return L.join('\n');
  }

  /* ---- InChI layers: which part of the identity differs ---- */
  function inchiLayers(inchi) {
    const parts = inchi.split('/'), out = new Map();
    let sec = 'main';
    out.set('main:formula', parts[1] || '');
    for (let k = 2; k < parts.length; k++) {
      const p = parts[k], L = p[0], body = p.slice(1);
      if (L === 'i') { sec = sec.startsWith('fixed') ? 'fixed-iso' : 'main-iso'; out.set(sec + ':i', body); continue; }
      if (L === 'f') { sec = 'fixed'; out.set('fixed:formula', body); continue; }
      if (L === 'r') { sec = 'reconnected'; out.set('reconnected:formula', body); continue; }
      out.set(sec + ':' + L, body);
    }
    return out;
  }
  const STEREO_LAYERS = new Set(['b', 't', 'm', 's']);
  function identityDifference(qI, pI) {   // null when only stereo layers differ (or nothing)
    if (qI === pI) return null;
    const A = inchiLayers(qI), B = inchiLayers(pI);
    const diff = [...new Set([...A.keys(), ...B.keys()])].filter(k => !STEREO_LAYERS.has(k.split(':')[1]) && A.get(k) !== B.get(k));
    if (!diff.length) return null;
    const has = re => diff.some(k => re.test(k));
    if (has(/^main:(formula|c)$/)) return 'compound';
    if (has(/:(q|p)$/)) return 'charge';
    if (has(/-iso:/)) return 'isotope';
    return 'tautomer';   // same formula, skeleton and charge; the H sit elsewhere
  }
  const MISMATCH = {
    compound: 'a different compound from',
    charge: 'a different charge state (protonation) of',
    isotope: 'a different isotope pattern of',
    tautomer: 'a different tautomer of',
    zwitterion: 'a different charge-separated (zwitterionic) form of',
    stereo: 'a different stereoisomer of',
  };
  const mismatchText = (kind, cid) => `PubChem's record (CID ${cid}) is ${MISMATCH[kind]} your structure, so its 3D model is not used for your molecule.`;

  /* ---- the query in record order: heavy atoms = M's atoms, then the record's H grouped by parent ---- */
  function queryOf(M, record) {
    const j = molJson(M), n = j.atoms.length, hOf = Array.from({ length: n }, () => []);
    const bad = msg => { throw new Fail('internal', 'Internal atom-mapping check failed (' + msg + '), so PubChem\'s 3D model is not used.'); };
    const elNo = (() => { const c = {}; return j.atoms.map(a => { const e = sym(a.z); return (c[e] = (c[e] || 0) + 1); }); })();
    const nm = i => sym(j.atoms[i].z) + elNo[i];   // "C3": per element, as the app names atoms
    let N = n;
    if (record) {
      const atoms = record.atoms || [];
      N = atoms.length;
      if (record.n !== n) bad('record and molecule differ in atom count');
      if (record.N != null && record.N !== N) bad('record atom count');
      for (let i = 0; i < n; i++) {
        const a = atoms[i], z = a.z != null ? a.z : SYMBOLS.indexOf(a.el) + 1;
        if (z !== j.atoms[i].z) bad('element of ' + nm(i));
      }
      for (let i = n; i < N; i++) {
        const a = atoms[i], p = a.parent;
        if (!(a.isH || a.el === 'H' || a.z === 1) || !(p >= 0 && p < n)) bad('added hydrogen number ' + (i + 1));
        hOf[p].push(i);
      }
      for (let i = 0; i < n; i++) if (hOf[i].length !== j.atoms[i].hs) bad('hydrogen count of ' + nm(i));
    } else {
      for (let i = 0; i < n; i++) for (let k = 0; k < j.atoms[i].hs; k++) hOf[i].push(N++);
    }
    const bonds = j.bonds.map(b => [b.a, b.b, b.bo]);
    hOf.forEach((hs, p) => hs.forEach(h => bonds.push([p, h, 1])));
    if (record && Array.isArray(record.bonds)) {   // the record's bond set must be this one
      const want = new Set(bonds.map(([a, b]) => pairKey(a, b)));
      if (record.bonds.length !== bonds.length || !record.bonds.every(b => want.has(pairKey(b.a, b.b)))) bad('bond list');
    }

    /* stereo the user specified: M's CIP labels (the record's come from the same labeller and must agree); open
       stereo = '?' atoms and "either" double bonds, from either source */
    const tags = stereoTags(M), cip = new Array(n).fill(null), ez = new Map();
    for (const [i, lab] of tags.atoms) cip[i] = lab;
    for (const [k, lab] of tags.bonds) ez.set(k, lab);
    const mb = M.get_molblock().split('\n'), counts = /V2000/.test(mb[3] || '') ? mb[3] : '';
    const mA = parseInt(counts.slice(0, 3), 10) || 0, mB = parseInt(counts.slice(3, 6), 10) || 0;
    for (let k = 4 + mA; k < 4 + mA + mB && k < mb.length; k++) {   // undefined double bonds: order 2, stereo 3 ("either")
      if (+mb[k].slice(6, 9) === 2 && +mb[k].slice(9, 12) === 3) ez.set(pairKey(+mb[k].slice(0, 3) - 1, +mb[k].slice(3, 6) - 1), 'unspecified');
    }
    if (record && record.atoms) {
      for (let i = 0; i < n; i++) {
        if (!('cip' in record.atoms[i])) continue;
        const r = record.atoms[i].cip || null, t = cip[i];
        if ((DEFINED.has(r) || DEFINED.has(t)) && r !== t) bad('stereo label of ' + nm(i));
        if (r === '?') cip[i] = '?';
      }
    }
    if (record && Array.isArray(record.bonds)) {
      for (const b of record.bonds) {
        if (!b || !('ez' in b)) continue;
        const k = pairKey(b.a, b.b), r = b.ez || null, t = ez.get(k) || null;
        if ((r === 'E' || r === 'Z' || t === 'E' || t === 'Z') && r !== t) bad('E/Z label of ' + nm(b.a) + '=' + nm(b.b));
        if (r === 'unspecified') ez.set(k, r);
      }
    }
    return { j, n, N, hOf, bonds, cip, ez, smiles: M.get_smiles() };
  }

  /* heavy atoms by substructure match (all symmetry-equivalent matches are tried), then H by parent */
  function chooseMatch(G, g, gTags, q) {
    let matches;
    try { matches = JSON.parse(G.get_substruct_matches(q.M, JSON.stringify({ uniquify: false, maxMatches: MAX_MATCHES }))); } catch (e) { matches = null; }
    if (!Array.isArray(matches)) matches = [];
    const nG = g.atoms.length, adj = Array.from({ length: nG }, () => []), gBonds = new Set();
    for (const b of g.bonds) { adj[b.a].push(b.b); adj[b.b].push(b.a); gBonds.add(pairKey(b.a, b.b)); }
    let sameTopology = 0;
    for (const m of matches) {
      const heavy = m.atoms;
      if (!heavy || heavy.length !== q.n) continue;
      if (!heavy.every((k, i) => { const a = g.atoms[k], b = q.j.atoms[i]; return a.z === b.z && a.chg === b.chg && a.isotope === b.isotope; })) continue;
      const perm = new Int32Array(q.N).fill(-1), used = new Set(heavy);
      heavy.forEach((k, i) => { perm[i] = k; });
      let ok = nG === q.N;
      for (let i = 0; i < q.n && ok; i++) {   // the k-th H of a parent (SDF order) → the k-th added H of that parent
        const hs = adj[heavy[i]].filter(k => g.atoms[k].z === 1 && !used.has(k)).sort((a, b) => a - b);
        if (hs.length !== q.hOf[i].length || hs.some(k => g.atoms[k].isotope || g.atoms[k].chg)) { ok = false; break; }
        hs.forEach((k, t) => { perm[q.hOf[i][t]] = k; used.add(k); });
      }
      if (!ok || used.size !== nG || g.bonds.length !== q.bonds.length || !q.bonds.every(([a, b]) => gBonds.has(pairKey(perm[a], perm[b])))) continue;
      sameTopology++;
      if (!q.cip.every((lab, i) => !DEFINED.has(lab) || gTags.atoms.get(perm[i]) === lab)) continue;
      let ezOk = true;
      for (const [k, lab] of q.ez) {
        if (lab !== 'E' && lab !== 'Z') continue;
        const [a, b] = k.split('-').map(Number);
        if (gTags.bonds.get(pairKey(perm[a], perm[b])) !== lab) { ezOk = false; break; }
      }
      if (ezOk) return perm;
    }
    /* at the cap the right match may be among those not listed, so no stereo verdict is given then */
    if (matches.length >= MAX_MATCHES) throw new Fail('internal', 'The molecule is too symmetric to check the atom mapping, so PubChem\'s 3D model is not used.');
    if (sameTopology) throw new Fail('stereo', null);
    throw new Fail('internal', 'Internal atom-mapping check failed (PubChem\'s atoms could not be matched to yours), so PubChem\'s 3D model is not used.');
  }

  /* V2000 molfile of the query in record order with the mapped coordinates; the stereo the user left open is marked
     "either" (wavy bond from the atom; stereo 3 on the double bond) so RDKit reads back only what was specified */
  function recordMolfile(q, coords) {
    const p3 = v => String(v).padStart(3), p4 = v => String(v).padStart(4), f = v => v.toFixed(4).padStart(10);
    const el = i => i < q.n ? sym(q.j.atoms[i].z) : 'H';
    const bonds = q.bonds.map(([a, b, bo]) => ({ a, b, bo, st: 0 }));
    bonds.forEach(bd => { if (bd.bo === 2 && q.ez.get(pairKey(bd.a, bd.b)) === 'unspecified') bd.st = 3; });
    q.cip.forEach((lab, i) => {
      if (lab !== '?') return;
      const bd = bonds.find(x => x.bo === 1 && !x.st && x.b === q.hOf[i][0] && x.a === i) || bonds.find(x => x.bo === 1 && !x.st && (x.a === i || x.b === i));
      if (bd) { if (bd.b === i) { bd.b = bd.a; bd.a = i; } bd.st = 4; }
    });
    const L = ['', '  OrgChem           3D', '', p3(q.N) + p3(bonds.length) + '  0  0  1  0  0  0  0  0999 V2000'];
    for (let i = 0; i < q.N; i++) L.push(f(coords[3 * i]) + f(coords[3 * i + 1]) + f(coords[3 * i + 2]) + ' ' + el(i).padEnd(3) + ' 0' + '  0'.repeat(11));
    for (const bd of bonds) L.push(p3(bd.a + 1) + p3(bd.b + 1) + p3(bd.bo) + p3(bd.st) + '  0  0  0');
    const prop = (tag, list) => { for (let k = 0; k < list.length; k += 8) { const c = list.slice(k, k + 8); L.push('M  ' + tag + p3(c.length) + c.map(([i, v]) => p4(i + 1) + p4(v)).join('')); } };
    prop('CHG', q.j.atoms.map((a, i) => [i, a.chg]).filter(x => x[1]));
    prop('ISO', q.j.atoms.map((a, i) => [i, a.isotope]).filter(x => x[1]));
    L.push('M  END');
    return L.join('\n');
  }

  /* ---- parse, check identity, map, verify ---- */
  function mapSdf(R, M, record, sdfText, info) {
    const made = [], keep = m => { if (m) made.push(m); return m; };
    let cid = info && info.cid;
    try {
      const parsed = splitSdf(sdfText || '');
      if (!parsed) return failure('parse', 'PubChem\'s 3D file could not be read.', { cid });
      const tags = parsed.tags;
      if (!cid) cid = parseInt(tags.PUBCHEM_COMPOUND_CID, 10) || null;
      const geoBlock = coordinateStereoOnly(parsed.block);
      const G = geoBlock && keep(getMol(R, geoBlock, false));
      if (!G) return failure('parse', 'PubChem\'s 3D file could not be read.', { cid });
      const g = molJson(G);
      if (!g.coords || g.coords.length !== g.atoms.length) return failure('parse', 'PubChem\'s file has no 3D coordinates.', { cid });
      const q = queryOf(M, record);
      q.M = M;

      /* identity: FixedH InChI (tautomer, charge, isotopes); stereo is checked on the mapping below */
      let qI = '', gI = '';
      try { qI = M.get_inchi('-FixedH'); gI = G.get_inchi('-FixedH'); } catch (e) { qI = gI = ''; }
      if (!qI || !gI) return failure('identity', 'PubChem\'s record could not be compared with your structure (no InChI), so its 3D model is not used.', { cid });
      let kind = identityDifference(qI, gI);
      const mismatch = k => {
        const P = keep(getMol(R, parsed.block, true));
        const udef = (parseInt(tags.PUBCHEM_ATOM_UDEF_STEREO_COUNT, 10) || 0);
        const mm = { kind: k, message: mismatchText(k, cid), pubchemSmiles: P ? P.get_smiles() : null, url: recordUrl(cid) };
        if (udef && mm.pubchemSmiles) mm.pubchemSmilesNote = `PubChem's record leaves ${udef} stereocentre${udef > 1 ? 's' : ''} undefined; this SMILES takes ${udef > 1 ? 'their configurations' : 'its configuration'} from PubChem's conformer.`;
        return failure('mismatch', mm.message, { cid, mismatch: mm });
      };
      if (kind === 'tautomer') {   // a proton moved onto/off a charged site: say zwitterion, not tautomer
        const charged = js => js.atoms.filter(a => a.chg).length;
        if (charged(q.j) !== charged(g)) kind = 'zwitterion';
      }
      if (kind) return mismatch(kind);

      const gTags = stereoTags(G);
      let perm;
      try { perm = chooseMatch(G, g, gTags, q); } catch (e) {
        if (e instanceof Fail && e.kind === 'stereo') return mismatch('stereo');
        throw e;
      }
      const coords = new Float64Array(3 * q.N);
      for (let i = 0; i < q.N; i++) for (let d = 0; d < 3; d++) coords[3 * i + d] = g.coords[perm[i]][d];
      /* stereo the user left open that RDKit does not flag in the query (e.g. cis/trans across a ring, CC1CCC(C)CC1):
         unlabelled in the query but labelled in PubChem's conformer, so it is open too */
      q.cip.forEach((lab, i) => { const got = gTags.atoms.get(perm[i]); if (lab == null && got && got !== '?') q.cip[i] = '?'; });
      for (const b of q.j.bonds) {
        const k = pairKey(b.a, b.b);
        if (b.bo === 2 && !q.ez.has(k) && gTags.bonds.has(pairKey(perm[b.a], perm[b.b]))) q.ez.set(k, 'unspecified');
      }

      /* final check, independent of the steps above: RDKit reads the result back in record order */
      const O = keep(getMol(R, recordMolfile(q, coords), true));
      if (!O || O.get_smiles() !== q.smiles) {
        return failure('internal', 'Internal atom-mapping check failed (the mapped 3D model does not read back as your structure), so PubChem\'s 3D model is not used.', { cid });
      }

      /* stereo the user left open: name what PubChem's conformer has there, atoms named as elsewhere in the app ("C3":
         the element, then the 1-based record index) */
      const elNo = (() => { const c = {}; return q.j.atoms.map(a => { const e = sym(a.z); return (c[e] = (c[e] || 0) + 1); }); })(), hWritten = q.j.atoms.filter(a => sym(a.z) === 'H').length;
      const assigned = [], atomWords = [], bondWords = [], name = i => (i < q.n ? sym(q.j.atoms[i].z) + elNo[i] : 'H' + (hWritten + i - q.n + 1));
      q.cip.forEach((lab, i) => {
        if (lab !== '?') return;
        const got = gTags.atoms.get(perm[i]);
        if (!got || got === '?') return;
        assigned.push({ kind: 'atom', index: i, label: got });
        atomWords.push(`${name(i)} is (${got})${got === 'r' || got === 's' ? ', pseudoasymmetric' : ''}`);
      });
      for (const [k, lab] of q.ez) {
        if (lab !== 'unspecified') continue;
        const [a, b] = k.split('-').map(Number), got = gTags.bonds.get(pairKey(perm[a], perm[b]));
        if (!got) continue;
        assigned.push({ kind: 'bond', atoms: [a, b], label: got });
        bondWords.push(`${name(a)}=${name(b)} is (${got})`);
      }
      const udef = (parseInt(tags.PUBCHEM_ATOM_UDEF_STEREO_COUNT, 10) || 0) + (parseInt(tags.PUBCHEM_BOND_UDEF_STEREO_COUNT, 10) || 0);
      const identity = { status: assigned.length ? 'stereo-undefined' : 'identical', assigned, mixedStereo: udef > 0 };
      if (assigned.length) {
        identity.banner = 'Stereochemistry not fully specified in your structure, so PubChem\'s 3D model shows one stereoisomer.' +
          (atomWords.length ? ' In the model ' + atomWords.join('; ') + '.' : '') +
          (bondWords.length ? ' Double-bond geometry not specified; in the model ' + bondWords.join('; ') + '.' : '') +
          (udef > 0 ? ' PubChem\'s other conformers of this record mix stereoisomers, so only its default conformer is used.' : '');
      }

      const conformerId = tags.PUBCHEM_CONFORMER_ID || null;
      const e = parseFloat(tags.PUBCHEM_MMFF94_ENERGY), energy = Number.isFinite(e) ? e : null;
      const energyNote = energy == null ? 'PubChem gives no energy for this conformer.'
        : `PubChem's energy for this conformer: ${fmt(energy, 2)} kcal/mol, ${TEXT.energyMethod}. Not comparable with the app's MMFF94 energies; PubChem's conformers are not energy minima.`;
      const url = recordUrl(cid);
      return {
        ok: true, source: 'PubChem', cid, conformerId, coords, sdf: sdfText, energy, energyNote, note: TEXT.caveat, identity,
        attribution: {
          text: `3D coordinates: PubChem CID ${cid}${conformerId ? ' · conformer ' + conformerId : ''} · ${TEXT.source} · NCBI/NLM data, public domain`,
          url, policyText: 'NCBI policies & disclaimer', policyUrl: POLICY_URL,
          exportTitle: `PubChem CID ${cid}${conformerId ? ' conformer ' + conformerId : ''} (PubChem3D computed, not experimental)`,
        },
      };
    } catch (e) {
      if (e instanceof Fail) return failure(e.kind, e.message, { cid });
      return failure('error', 'PubChem 3D failed: ' + (e && e.message || e), { cid });
    } finally {
      for (const m of made) { try { m.delete(); } catch (e) {} }
    }
  }

  /* ---- how a compound is made: PubChem's "Methods of Manufacturing" section ----
     manufacturing(R, M) → Promise<{ok, cid, title, entries:[{text, cite, note}], source, url, policyUrl} | {ok:false, kind, error}>
     The entries come from the Hazardous Substances Data Bank (NLM), each with the publication HSDB quotes (mostly
     Ullmann's Encyclopedia of Industrial Chemistry or Kirk-Othmer). Two requests, through the same queue and rules as
     the 3D model: the SMILES → PubChem's CID, title and InChIKey, then the section. The record is used only when its
     standard InChIKey equals the one computed here for this structure (a lookup can land on another form). */
  const mfgCache = new Map();
  async function manufacturing(R, M) {
    let smiles = '', key = '';
    try { smiles = M.get_smiles(); key = R.get_inchikey_for_inchi(M.get_inchi()); } catch (e) { /* below */ }
    if (!smiles || !key) return failure('ineligible', 'There is no structure to look up.');
    if (smiles.includes('.')) return failure('ineligible', 'Look up one component at a time, not a salt or mixture.');
    if (mfgCache.has(smiles)) return mfgCache.get(smiles);
    if (busyFor()) return failure('busy', TEXT.busy);
    if (isOffline()) return failure('offline', 'Offline — this information comes from PubChem and needs an internet connection.');
    try {
      const r = await request(API + '/compound/smiles/property/Title,InChIKey/JSON', { method: 'POST', body: new URLSearchParams({ smiles }) });
      if (r.status === 400) return failure('rejected', 'PubChem could not read this structure.');
      if (r.status === 404) return remember(smiles, failure('not-found', 'This exact structure is not in PubChem.'));
      if (r.status !== 200) return failure('http', `PubChem answered with an error (HTTP ${r.status}).`);
      let p = null;
      try { p = JSON.parse(r.text).PropertyTable.Properties; } catch (e) { p = null; }
      if (!Array.isArray(p) || !p.length) return failure('http', 'PubChem\'s answer could not be read.');
      if (p.length > 1) return failure('ambiguous', `PubChem matched this structure to ${p.length} records, so none is shown.`);
      const { CID: cid, Title: title, InChIKey: theirs } = p[0];
      if (!cid) return remember(smiles, failure('not-found', 'This exact structure is not in PubChem.'));
      if (theirs !== key) {
        const same = theirs && theirs.slice(0, 14) === key.slice(0, 14);
        return remember(smiles, failure('mismatch', `PubChem's closest record (CID ${cid}, ${title || 'untitled'}) is ${same ? 'a different stereoisomer or charge/isotope form of' : 'a different compound from'} this structure, so its manufacturing notes are not shown.`, { cid }));
      }
      const v = await request(`https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/${cid}/JSON?heading=Methods+of+Manufacturing`);
      if (v.status === 404) return remember(smiles, failure('none', `PubChem (CID ${cid}, ${title}) has no "Methods of Manufacturing" section for this compound.`, { cid }));
      if (v.status !== 200) return failure('http', `PubChem answered with an error (HTTP ${v.status}).`);
      let rec = null;
      try { rec = JSON.parse(v.text).Record; } catch (e) { rec = null; }
      if (!rec) return failure('parse', 'PubChem\'s answer could not be read.');
      const refs = new Map((rec.Reference || []).map(x => [x.ReferenceNumber, x]));
      const entries = [], sources = new Set();
      (function walk(s) {
        for (const sec of s.Section || []) {
          for (const inf of sec.Information || []) {
            const text = ((inf.Value || {}).StringWithMarkup || []).map(x => x.String || '').join(' ').trim();
            if (!text) continue;
            const cite = (inf.Reference || [])[0] || ((inf.ExtendedReference || [])[0] || {}).Citation || '';
            const src = refs.get(inf.ReferenceNumber);
            if (src && src.SourceName) sources.add(src.SourceName);
            entries.push({ text, cite, note: inf.Description || '' });
          }
          walk(sec);
        }
      })(rec);
      if (!entries.length) return remember(smiles, failure('none', `PubChem (CID ${cid}, ${title}) has no manufacturing notes for this compound.`, { cid }));
      return remember(smiles, { ok: true, cid, title, entries, source: [...sources].join(', ') || 'PubChem',
        url: `https://pubchem.ncbi.nlm.nih.gov/compound/${cid}#section=Methods-of-Manufacturing`, policyUrl: POLICY_URL });
    } catch (e) {
      if (e instanceof Fail) return failure(e.kind, e.kind === 'offline' ? 'Could not reach PubChem (no internet connection, or the request was blocked).' : e.message);
      return failure('error', 'The PubChem lookup failed: ' + (e && e.message || e));
    }
  }
  function remember(smiles, v) { mfgCache.set(smiles, v); return v; }   // answers that will not change this session

  /* ---- literature: papers that mention two compounds together (PubChem's "Chemical Co-Occurrences in Literature")
     PubChem links each paper to the compounds it mentions. Two compounds in one paper does not mean the paper
     describes a reaction between them, and the callers say so. Answers are kept in this browser (localStorage), so a
     pair looked up once is there offline afterwards: the user's own growing record of the literature.
     identify(R, smiles) → Promise<{ok, cid, title, inchikey} | {ok:false, kind, error}>   exact standard InChIKey match only
     papersBoth(cidA, cidB, {refresh}) → Promise<{ok, count, articles:[{title, authors, journal, date, review, doi, pmid}],
                                                  fetched, saved} | {ok:false, kind, error}>   PubChem's 5 most relevant */
  const LS = {
    get(k) { try { const v = localStorage.getItem(k); return v ? JSON.parse(v) : null; } catch (e) { return null; } },
    set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { /* storage full or blocked: not kept */ } },
  };
  const idCache = new Map();
  async function identify(R, smiles) {
    let m = null, key = '', canon = '';
    try { m = R.get_mol(smiles); if (m && m.get_num_atoms()) { canon = m.get_smiles(); key = R.get_inchikey_for_inchi(m.get_inchi()); } }
    catch (e) { /* below */ } finally { if (m) m.delete(); }
    if (!key) return failure('ineligible', 'RDKit could not read this structure.');
    if (canon.includes('.')) return failure('ineligible', 'A salt or mixture: PubChem is asked about one compound at a time.');
    if (idCache.has(key)) return idCache.get(key);
    const saved = LS.get('pc-id:' + key);
    if (saved) { idCache.set(key, saved); return saved; }
    if (busyFor()) return failure('busy', TEXT.busy);
    const keep = v => { idCache.set(key, v); LS.set('pc-id:' + key, v); return v; };
    try {
      const r = await request(API + '/compound/smiles/property/Title,InChIKey/JSON', { method: 'POST', body: new URLSearchParams({ smiles: canon }) });
      if (r.status === 400) return failure('rejected', 'PubChem could not read this structure.');
      if (r.status === 404) return keep(failure('not-found', 'This exact structure is not in PubChem.'));
      if (r.status !== 200) return failure('http', `PubChem answered with an error (HTTP ${r.status}).`);
      let p = null;
      try { p = JSON.parse(r.text).PropertyTable.Properties; } catch (e) { p = null; }
      if (!Array.isArray(p) || p.length !== 1) return failure('http', 'PubChem\'s answer could not be read.');
      const { CID: cid, Title: title, InChIKey: theirs } = p[0];
      if (!cid) return keep(failure('not-found', 'This exact structure is not in PubChem.'));
      if (theirs !== key) return keep(failure('mismatch', `PubChem's closest record (CID ${cid}, ${title || 'untitled'}) is not exactly this structure, so it is not used.`, { cid }));
      return keep({ ok: true, cid, title: title || `CID ${cid}`, inchikey: key });
    } catch (e) {
      if (e instanceof Fail) return failure(e.kind, e.kind === 'offline' ? 'Could not reach PubChem (no internet connection, or the request was blocked).' : e.message);
      return failure('error', 'The PubChem lookup failed: ' + (e && e.message || e));
    }
  }
  async function papersBoth(cidA, cidB, opts) {
    /* (pc-pair2: answers saved under the old pc-pair: key asked one direction only and are ignored) */
    const a = Math.min(cidA, cidB), b = Math.max(cidA, cidB), k = `pc-pair2:${a}-${b}`;
    if (!(opts && opts.refresh)) { const s = LS.get(k); if (s) return Object.assign({}, s, { saved: true }); }
    if (a === b) return failure('same', 'Both are the same compound.');
    if (busyFor()) return failure('busy', TEXT.busy);
    try {
      /* The service answers from id_1's list of co-mentioned compounds, and each list keeps only that compound's most
         co-mentioned neighbours. So the direction matters: indole (798) → indoline (10328) says "No links found",
         while indoline → indole finds 884 papers, because indole is mentioned with so many compounds that indoline
         is not on its list. Ask one way, and if that finds nothing, the other. */
      const ask = async (x, y) => {
        const r = await request(`https://pubchem.ncbi.nlm.nih.gov/link_db/link_db_server.cgi?format=JSON&type=ChemicalNeighbor&operation=GetAllLinks&id_1=${x}&id_2=${y}&response_type=display`);
        let j = null;
        try { j = JSON.parse(r.text); } catch (e) { j = null; }
        // "No links found": the pair is not on x's list (not an error, and not proof that no paper mentions both)
        if (j && j.Fault && /no links found/i.test(j.Fault.Message || '')) j = { LinkDataSet: { LinkData: [] } };
        if (!j || r.status !== 200 || j.Fault) throw new Fail('http', 'PubChem\'s literature service answered with an error' + (j && j.Fault && j.Fault.Message ? ': ' + j.Fault.Message : ` (HTTP ${r.status}).`));
        return j;
      };
      let j = await ask(b, a);                             // the larger CID is usually the rarer compound: its list is likelier to hold the other
      if (!(((j.LinkDataSet || {}).LinkData || [])[0])) { if (busyFor()) return failure('busy', TEXT.busy); j = await ask(a, b); }
      const ev = (((j.LinkDataSet || {}).LinkData || [])[0] || {}).Evidence;
      const cn = ev && ev.ChemicalNeighbor;
      const articles = ((cn && cn.Article) || []).map(x => ({ title: x.Title || '', authors: x.Author || '', journal: x.Journal || '',
        date: x.PublicationDate || '', review: !!x.IsReview, doi: x.DOI || '', pmid: x.PMID || null }));
      const res = { ok: true, count: cn ? +cn.ArticleCount || articles.length : 0, articles, fetched: new Date().toISOString().slice(0, 10) };
      LS.set(k, res);
      return Object.assign({}, res, { saved: false });
    } catch (e) {
      if (e instanceof Fail) return failure(e.kind, e.kind === 'offline' ? 'Could not reach PubChem (no internet connection, or the request was blocked).' : e.message);
      return failure('error', 'The PubChem lookup failed: ' + (e && e.message || e));
    }
  }

  /* ---- a name or CAS number the app's offline library does not know → PubChem ----
     lookupName(text) → Promise<{ok, cid, title, smiles, formula, cas, url, note} | {ok:false, kind, error, choices?}>
     PubChem matches the text against its synonyms (names, CAS numbers, other identifiers). When it finds several
     records that are different compounds (a different connectivity), none is picked: the choices are returned so the
     user can choose. CAS numbers in PubChem are synonyms supplied by depositors, and the note says so. Answers are kept
     in this browser, so a name looked up once works offline afterwards. Only on an explicit click. */
  async function lookupName(text) {
    const q = String(text || '').trim();
    if (!q) return failure('ineligible', 'Nothing to look up.');
    const k = 'pc-name:' + q.toLowerCase();
    const saved = LS.get(k);
    if (saved) return Object.assign({}, saved, { saved: true });
    if (busyFor()) return failure('busy', TEXT.busy);
    if (isOffline()) return failure('offline', 'Offline — looking a name up on PubChem needs an internet connection.');
    const isCas = /^\d{2,7}-\d\d-\d$/.test(q);
    try {
      const r = await request(API + '/compound/name/property/Title,SMILES,MolecularFormula,InChIKey/JSON', { method: 'POST', body: new URLSearchParams({ name: q }) });
      if (r.status === 404) { const f = failure('not-found', `PubChem does not know "${q}" either.`); LS.set(k, f); return f; }
      if (r.status === 400) return failure('rejected', `PubChem could not use "${q}" as a name.`);
      if (r.status !== 200) return failure('http', `PubChem answered with an error (HTTP ${r.status}).`);
      let p = null;
      try { p = JSON.parse(r.text).PropertyTable.Properties; } catch (e) { p = null; }
      if (!Array.isArray(p) || !p.length) return failure('http', 'PubChem\'s answer could not be read.');
      const blocks = new Set(p.map(x => String(x.InChIKey || '').slice(0, 14)));
      if (blocks.size > 1) return failure('ambiguous', `PubChem has ${p.length} different compounds under "${q}". Choose one, or paste its SMILES.`,
        { choices: p.slice(0, 8).map(x => ({ cid: x.CID, title: x.Title || `CID ${x.CID}`, smiles: x.SMILES, formula: x.MolecularFormula })) });
      const x = p[0];
      if (!x.SMILES) return failure('http', 'PubChem\'s record has no structure.');
      const res = { ok: true, cid: x.CID, title: x.Title || `CID ${x.CID}`, smiles: x.SMILES, formula: x.MolecularFormula || '', cas: isCas ? q : '',
        url: `https://pubchem.ncbi.nlm.nih.gov/compound/${x.CID}`,
        note: `From PubChem: CID ${x.CID}, ${x.Title || ''}` + (isCas ? ' (CAS numbers in PubChem are synonyms supplied by depositors; check an important one with CAS Common Chemistry)' : '') +
          (p.length > 1 ? `; ${p.length} PubChem records match, all the same compound apart from stereo or form` : '') + '.' };
      LS.set(k, res);
      return Object.assign({}, res, { saved: false });
    } catch (e) {
      if (e instanceof Fail) return failure(e.kind, e.kind === 'offline' ? 'Could not reach PubChem (no internet connection, or the request was blocked).' : e.message);
      return failure('error', 'The PubChem lookup failed: ' + (e && e.message || e));
    }
  }

  // a name looked up on PubChem before, from this browser's saved record (no request): the search boxes use it first
  function savedName(text) {
    const s = LS.get('pc-name:' + String(text || '').trim().toLowerCase());
    return s && s.ok ? s : null;
  }

  return { eligible, load, mapSdf, busyFor, manufacturing, identify, papersBoth, lookupName, savedName, TEXT, CREDITS };
})();
