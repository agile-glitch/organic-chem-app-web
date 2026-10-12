/* Atlas — what is known about a functional group IN ITS SURROUNDINGS, learned from recorded reactions, for groups never met before.

   A fragment instance is one of the 205 named groups of the Fragments sheet at one place in one molecule. Its identity is the group plus the
   surroundings of its reacting atom at radius 1, 2 and 3 (tools/atlas_core.py: a Weisfeiler-Lehman hash of element, aromaticity, branching, H count,
   charge, ring membership and bond kinds), so the number of distinct fragments is not a list somebody wrote: it is every environment that occurs
   (millions are in the tables, billions are possible). tools/atlas_mine.py counted, over 478,000 atom-mapped patent reactions, how often each was
   TRANSFORMED by the reaction and how often it was left alone, per kind of reagent and per named reagent; data/atlas/ holds those counts.

   For a molecule and a reagent the page asks, for every group:  how likely is it that this reaction transforms it?  The answer backs off from the
   most specific context that has evidence to the coarser ones (p_k = (h_k + M p_(k-1)) / (n_k + M)), then a small logistic model fitted on held-out
   reactions combines the chains and the evidence (tools/atlas_eval.py: figures in data/atlas/meta.js). A context nobody recorded is still judged, from
   what its parts are like, and says so (evidence level 'none'), with the model's own confidence.

   Atlas.ready()                          -> Promise<meta | null>   (null: data/atlas is not there; the page works without)
   Atlas.describe(R, smiles)              -> [{F, fid, name, atom, e: [e0..e3], steric, charge, ring, strain, hyb, nAtoms}]   one per group instance
   Atlas.judge(R, subs, reagentLine)      -> [{smiles, fid, name, atom, p, level, n, h, chain[], confidence, ...}]  (all instances of all subs)
   Atlas.explain(j)                       -> one line of text for an instance judgement
   Atlas.mix32 / fnv32 / wl(record)       -> the hashes, for the parity test (tools/atlas_parity.js)                                                   */
(() => {
  'use strict';
  const SKIP_ROLES = new Set(['counter-ion / spectator ion', 'solvent', 'salt (work-up or additive)', 'other', 'drying agent']);
  const FAM = { F: 1, R0: 2, R1: 3, R2: 4, R3: 5, U1: 6, U2: 7, U3: 8, P0: 9, P1: 10 };

  // ---------------------------------------------------------------- hashing (tools/atlas_keys.py, tools/atlas_core.py)
  const mix32 = (a, b) => {
    let x = (a ^ Math.imul(b ^ 0x9E3779B9, 0x85EBCA6B)) | 0;
    x ^= x >>> 16; x = Math.imul(x, 0x85EBCA6B); x ^= x >>> 13; x = Math.imul(x, 0xC2B2AE35); x ^= x >>> 16;
    return x >>> 0;
  };
  const key = (fam, F, a, b) => mix32(mix32(mix32(FAM[fam], F), a || 0), b || 0);
  const enc = new TextEncoder();
  const fnv32 = s => { let h = 0x811c9dc5; for (const c of enc.encode(s)) h = Math.imul(h ^ c, 0x01000193); return h >>> 0; };

  /* environment hashes of every heavy atom at radius 0..3, from a MolInfo record (js/molinfo.js): the same refinement as atlas_core.wl_hashes */
  function wl(rec) {
    const A = rec.atoms, n = rec.n;
    const heavyNb = i => A[i].neighbours.filter(x => x.atom < n && !A[x.atom].isH);
    const inv0 = i => `${A[i].z},${A[i].aromatic ? 1 : 0},${heavyNb(i).length},${A[i].hCount},${A[i].charge},${(A[i].ringSizes || []).length ? 1 : 0}`;
    let cur = []; for (let i = 0; i < n; i++) cur.push(fnv32(inv0(i)));
    const out = [cur];
    for (let r = 1; r < 4; r++) {
      const nxt = [];
      for (let i = 0; i < n; i++) {
        const nb = heavyNb(i).map(x => `${x.aromatic ? 'A' : String(Math.trunc(x.order))}:${cur[x.atom]}`).sort();
        nxt.push(fnv32(`${cur[i]}[${nb.join(',')}]`));
      }
      cur = nxt; out.push(cur);
    }
    return out;
  }

  // ---------------------------------------------------------------- the data
  let metaP = null, META = null;
  const loadScript = src => new Promise((res, rej) => {
    const el = document.createElement('script'); el.src = src;
    el.onload = () => { el.remove(); res(); }; el.onerror = () => { el.remove(); rej(new Error('could not load ' + src)); };
    document.head.appendChild(el);
  });
  function ready() {
    if (!metaP) metaP = (async () => {
      try { await loadScript('data/atlas/meta.js'); META = window.ATLAS_META || null; } catch (e) { META = null; }
      return META;
    })();
    return metaP;
  }
  const shards = new Map();
  async function gunzipB64(b64) {
    const bin = atob(b64), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    return new Response(new Blob([u8]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
  }
  function parseShard(buf) {
    const dv = new DataView(buf), u8 = new Uint8Array(buf), count = dv.getUint32(0, true);
    let pos = 4;
    const rd = () => { let v = 0, m = 1, b; do { b = u8[pos++]; v += (b & 0x7f) * m; m *= 128; } while (b & 0x80); return v; };
    const k = new Float64Array(count), n = new Uint32Array(count), h = new Uint32Array(count);
    let prev = 0;
    for (let i = 0; i < count; i++) { prev += rd(); k[i] = prev; n[i] = rd(); h[i] = rd(); }
    return { k, n, h };
  }
  function shard(F) {
    const name = 'f/' + String(F).padStart(3, '0');
    if (!shards.has(name)) shards.set(name, (async () => {
      if (!META || !(F in META.shardEntries)) return null;
      try { await loadScript('data/atlas/' + name + '.js'); } catch (e) { return null; }
      const P = (window.ATLAS_PARTS || {})[name];
      if (!P) return null;
      delete window.ATLAS_PARTS[name];
      try { return parseShard(await gunzipB64(P.gz)); } catch (e) { shards.delete(name); return null; }
    })());
    return shards.get(name);
  }
  function get(sh, k) {                                    // [met, transformed] or [0, 0]
    if (!sh) return [0, 0];
    let lo = 0, hi = sh.k.length - 1;
    while (lo <= hi) {
      const mid = (lo + hi) >> 1, v = sh.k[mid];
      if (v === k) return [sh.n[mid], sh.h[mid]];
      if (v < k) lo = mid + 1; else hi = mid - 1;
    }
    return [0, 0];
  }

  // ---------------------------------------------------------------- the instances of a molecule
  const qcache = new Map();
  const qmol = (R, id, smarts) => { if (!qcache.has(id)) { let q = null; try { q = R.get_qmol(smarts); } catch (e) { q = null; } qcache.set(id, q); } return qcache.get(id); };
  function reactiveAtom(rec, atoms) {                     // js/reactions.js effectsOfMol pick; tools/atlas_core.reactive_atom
    const A = rec.atoms, n = rec.n;
    const heavy = i => A[i].neighbours.filter(x => x.atom < n && !A[x.atom].isH).map(x => x.atom);
    const bo = (i, j) => { const e = A[i].neighbours.find(x => x.atom === j); return e ? (e.aromatic ? 1.5 : e.order) : 0; };
    const isC = i => A[i].el === 'C', multi = (i, els) => heavy(i).some(j => els.includes(A[j].el) && bo(i, j) >= 2);
    return atoms.find(i => isC(i) && multi(i, ['O', 'N', 'S']))
      ?? atoms.find(i => isC(i) && heavy(i).some(j => ['F', 'Cl', 'Br', 'I'].includes(A[j].el)))
      ?? atoms.find(i => ['N', 'O', 'S'].includes(A[i].el) && A[i].hCount > 0)
      ?? atoms.find(i => ['N', 'O', 'S'].includes(A[i].el))
      ?? atoms.find(i => isC(i) && heavy(i).some(j => A[j].el === 'C' && bo(i, j) >= 2 && bo(i, j) !== 1.5))
      ?? atoms[0];
  }
  function stericScore(rec, x) {                          // tools/atlas_core.steric_score
    const A = rec.atoms, n = rec.n;
    const heavy = i => A[i].neighbours.filter(e => e.atom < n && !A[e.atom].isH).map(e => e.atom);
    const bo = (i, j) => { const e = A[i].neighbours.find(y => y.atom === j); return e ? (e.aromatic ? 1.5 : e.order) : 0; };
    const partner = heavy(x).filter(j => !(['O', 'S'].includes(A[j].el) && bo(x, j) === 2));
    let k = 0;
    for (const j of partner) {
      let kk = 0;
      if (A[j].aromatic) {
        const ring = (rec.rings || []).find(r => r.length === 6 && r.includes(j) && r.every(q => A[q].aromatic));
        if (ring) { const p = ring.indexOf(j); kk = [ring[(p + 1) % 6], ring[(p + 5) % 6]].filter((o, i, a) => a.indexOf(o) === i).reduce((t, o) => t + heavy(o).filter(q => !ring.includes(q)).length, 0); }
      } else kk = heavy(j).filter(q => q !== x).length;
      k = Math.max(k, kk);
    }
    return partner.length + k;
  }
  function hybridization(rec, x) {                        // tools/atlas_core.hybridization: 1 sp, 2 sp2 / aromatic, 3 sp3
    const a = rec.atoms[x];
    if (a.aromatic) return 2;
    const orders = a.neighbours.filter(e => e.atom < rec.n).map(e => e.order);
    if (orders.some(o => o === 3) || orders.filter(o => o === 2).length >= 2) return 1;
    return orders.some(o => o === 2) ? 2 : 3;
  }
  const smallestRing = (rec, x) => { const s = rec.atoms[x].ringSizes || []; return s.length ? Math.min(...s) : 0; };
  function angleStrain(rec, x) {
    const n = smallestRing(rec, x);
    if (!n) return 0;
    const ideal = { 3: 109.5, 2: 120, 1: 180 }[hybridization(rec, x)];
    return Math.abs((n - 2) * 180 / n - ideal);
  }

  /* the reacting atom and what it is bonded to, for people:  C(=O)(-O)(-c)  (lower case = aromatic) */
  function whereText(rec, x) {
    const A = rec.atoms, a = A[x], sym = i => (A[i].aromatic ? A[i].el.toLowerCase() : A[i].el);
    const nb = a.neighbours.filter(e => e.atom < rec.n && !A[e.atom].isH).map(e => (e.aromatic ? ':' : e.order === 2 ? '=' : e.order === 3 ? '#' : '-') + sym(e.atom));
    nb.sort((p, q) => (p[0] === '=' ? 0 : 1) - (q[0] === '=' ? 0 : 1) || p.localeCompare(q));
    return sym(x) + (a.hCount ? 'H' + (a.hCount > 1 ? a.hCount : '') : '') + nb.map(t => `(${t})`).join('');
  }

  const descCache = new Map();
  function describe(R, smiles) {
    if (descCache.has(smiles)) return descCache.get(smiles);
    const RR = window.REACTION_RULES || {}, ids = META ? META.fragments : Object.keys(RR.fragments || {}), out = [];
    let m = null;
    try {
      m = R.get_mol(smiles);
      if (m && window.MolInfo) {
        const rec = window.MolInfo.analyse(R, m), h = wl(rec);
        ids.forEach((id, F) => {
          const f = (RR.fragments || {})[id]; if (!f) return;
          const q = qmol(R, id, f.smarts); if (!q) return;
          let ms = []; try { const r = JSON.parse(m.get_substruct_matches(q) || '[]'); ms = Array.isArray(r) ? r : (r && r.atoms ? [r] : []); } catch (e) { ms = []; }   // one match comes back as an object, several as an array
          for (const mt of ms.slice(0, 8)) {
            const atoms = mt.atoms.filter(i => i < rec.n), x = reactiveAtom(rec, atoms);
            if (x == null) continue;
            const g = rec.atoms[x].gasteiger;
            out.push({ F, fid: id, name: f.name, atom: x, atoms, where: whereText(rec, x), e: [h[0][x], h[1][x], h[2][x], h[3][x]], steric: stericScore(rec, x),
              charge: g == null || !isFinite(g) ? null : g, ring: smallestRing(rec, x), strain: angleStrain(rec, x), hyb: hybridization(rec, x), nAtoms: atoms.length });
          }
        });
      }
    } catch (e) { /* no instances */ } finally { if (m) m.delete(); }
    descCache.set(smiles, out);
    return out;
  }

  // ---------------------------------------------------------------- reagent side of an instance
  const canonOf = (R, s) => { let m = null; try { m = R.get_mol(s); return m ? m.get_smiles() : s; } catch (e) { return s; } finally { if (m) m.delete(); } };
  const keyOf = (R, s) => { let m = null; try { m = R.get_mol(s); return m ? R.get_inchikey_for_inchi(m.get_inchi()) : ''; } catch (e) { return ''; } finally { if (m) m.delete(); } };
  const partHash = (R, s) => fnv32(keyOf(R, s) || s);
  const acidPka = new Map();
  function acidTable(R) {
    if (!acidPka.size) for (const e of (window.REACTION_RULES || {}).acidsBases || []) if (e.role === 'acid' && e.pKa != null) acidPka.set(canonOf(R, e.smiles), e.pKa);
    return acidPka;
  }
  /* the reagent side for the instance in molecule `own`: the reagent line's parts and the other starting materials that the Reagents sheet knows */
  function reagentSide(R, own, subs, reagentLine) {
    const RG = (window.REACTION_RULES || {}).reagents || {}, seen = new Set(), parts = [];
    const add = s => { if (seen.has(s)) return; seen.add(s); const r = RG[s]; if (r && !SKIP_ROLES.has(r.role)) parts.push({ smiles: s, role: r.role, pKaH: r.pKaH }); };
    String(reagentLine || '').split('.').filter(Boolean).forEach(p => add(canonOf(R, p)));
    subs.forEach(s => { if (s !== own) add(canonOf(R, s)); });
    parts.sort((a, b) => (a.smiles < b.smiles ? -1 : a.smiles > b.smiles ? 1 : 0));
    const chosen = parts.slice(0, 3);
    const roles = [...new Set(parts.map(p => p.role))].sort();
    return { parts: chosen, cls: roles.length ? roles.slice(0, 3).join('|') : '-' };
  }

  // ---------------------------------------------------------------- the judgement
  const shrink = (prev, n, h, m) => (h + m * prev) / (n + m);
  const logit = p => { p = Math.min(1 - 1e-4, Math.max(1e-4, p)); return Math.log(p / (1 - p)); };
  const sigmoid = z => 1 / (1 + Math.exp(-Math.max(-30, Math.min(30, z))));

  function features(R, inst, side, sh) {
    const M = META.m, minN = META.minN, F = inst.F, e = inst.e;
    const clsH = fnv32(side.cls);
    const look = (fam, a, b) => { const [n, h] = get(sh, key(fam, F, a, b)); return n >= minN[fam] ? [n, h] : [0, 0]; };
    const [nF, hF] = look('F'); const pF = shrink(META.base, nF, hF, 2.0);
    let [n, h] = look('R0', clsH); const pR0 = shrink(pF, n, h, M), nR0 = n;
    const chain = { F: [pF, nF], R0: [pR0, nR0] };
    let prev = pR0;
    for (const lv of [1, 2, 3]) { [n, h] = look('R' + lv, clsH, e[lv]); prev = shrink(prev, n, h, M); chain['R' + lv] = [prev, n, h]; }
    prev = pF;
    for (const lv of [1, 2, 3]) { [n, h] = look('U' + lv, e[lv]); prev = shrink(prev, n, h, M); chain['U' + lv] = [prev, n, h]; }
    let bestN = -1, pP0 = pR0, pP1 = pR0, nP0 = 0, nP1 = 0;
    for (const p of side.parts) {
      const ph = partHash(R, p.smiles);
      const [n0, h0] = look('P0', ph), [n1, h1] = look('P1', ph, e[1]);
      const p0 = shrink(pF, n0, h0, M), p1 = shrink(p0, n1, h1, M);
      if (n0 > bestN) { bestN = n0; pP0 = p0; pP1 = p1; nP0 = n0; nP1 = n1; }
    }
    chain.P0 = [pP0, nP0]; chain.P1 = [pP1, nP1];
    const v = {};
    for (const k of ['F', 'R0', 'R1', 'R2', 'R3', 'U1', 'U2', 'U3', 'P0', 'P1']) { v['logit_' + k] = logit(chain[k][0]); if (k !== 'F') v['logn_' + k] = Math.log1p(chain[k][1]); }
    // descriptors of the reacting atom (tools/atlas_eval.py; present only when the fitted model uses them)
    const q = inst.charge, RR = window.REACTION_RULES || {}, fr = (RR.fragments || {})[inst.fid] || {};
    Object.assign(v, {
      steric: Math.min(9, Math.max(0, inst.steric)), steric_sq: Math.min(9, Math.max(0, inst.steric)) ** 2 / 9,
      charge: q == null ? 0 : q, charge_missing: q == null ? 1 : 0, abs_charge: q == null ? 0 : Math.abs(q),
      ring3: +(inst.ring === 3), ring4: +(inst.ring === 4), ring5: +(inst.ring === 5), ring6: +(inst.ring === 6), ring7p: +(inst.ring >= 7), in_ring: +(inst.ring > 0),
      strain: inst.strain / 30, sp: +(inst.hyb === 1), sp2: +(inst.hyb === 2), match_atoms: Math.min(12, inst.nAtoms) / 12,
    });
    let pb = NaN, pa = NaN;
    const at = acidTable(R);
    for (const p of side.parts) {
      if (p.role === 'base' && p.pKaH != null && !(pb >= p.pKaH)) pb = p.pKaH;
      if (p.role === 'acid' && at.has(p.smiles) && !(pa <= at.get(p.smiles))) pa = at.get(p.smiles);
    }
    const dpb = pb - (fr.pKa == null ? NaN : fr.pKa), dpa = (fr.pKaH == null ? NaN : fr.pKaH) - pa;
    const clip = (x, lo, hi) => (isNaN(x) ? 0 : Math.max(lo, Math.min(hi, x))) / 10;
    Object.assign(v, { has_base: +!isNaN(pb), dpka_base: clip(dpb, -25, 35), has_acid: +!isNaN(pa), dpka_acid: clip(dpa, -25, 35) });
    return { v, chain };
  }
  function combine(v) {
    const M = META.model;
    let z = M.coef[0];
    M.names.forEach((nm, i) => { const x = v[nm]; if (x === undefined) throw new Error('atlas: the model wants ' + nm); z += (x - M.mean[i]) / M.std[i] * M.coef[i + 1]; });
    return sigmoid(z);
  }
  const LEVELS = [['R3', 'radius 3, this kind of reagent'], ['U3', 'radius 3, any reagent'], ['R2', 'radius 2, this kind of reagent'], ['U2', 'radius 2, any reagent'],
                  ['R1', 'radius 1, this kind of reagent'], ['U1', 'radius 1, any reagent'], ['R0', 'this group with this kind of reagent'], ['F', 'this group']];
  /* how well documented the answer is: the deepest context that has evidence */
  function confidenceOf(chain) {
    for (const [k, label] of LEVELS) { const c = chain[k]; if (c && c[1] > 0) return { level: k, label, n: c[1], h: c[2] == null ? null : c[2] }; }
    return { level: 'none', label: 'no recorded reaction', n: 0, h: null };
  }
  function band(conf) {
    if (conf.level === 'R3' || conf.level === 'U3') return conf.n >= 20 ? 'high' : 'medium';
    if (conf.level === 'R2' || conf.level === 'U2') return conf.n >= 20 ? 'high' : 'medium';
    if (conf.level === 'R1' || conf.level === 'U1') return conf.n >= 50 ? 'medium' : 'low';
    return 'low';
  }

  async function judge(R, subs, reagentLine) {
    await ready();
    if (!META || !R) return [];
    const out = [];
    for (const smi of subs) {
      const side = reagentSide(R, smi, subs, reagentLine);
      for (const inst of describe(R, smi)) {
        const sh = await shard(inst.F);
        const { v, chain } = features(R, inst, side, sh);
        let p; try { p = combine(v); } catch (e) { p = chain.R1[0]; }
        const conf = confidenceOf(chain);
        out.push({ smiles: smi, fid: inst.fid, name: inst.name, atom: inst.atom, p, chain, evidence: conf, confidence: band(conf), cls: side.cls,
                   parts: side.parts.map(x => x.smiles), steric: inst.steric, charge: inst.charge, ring: inst.ring, strain: inst.strain, hyb: inst.hyb });
      }
    }
    return out;
  }
  const pct = x => (x >= 0.995 ? '>99' : x < 0.005 ? '<1' : String(Math.round(x * 100))) + ' %';
  function explain(j) {
    const c = j.evidence;
    const ev = c.level === 'none' ? 'no recorded reaction of this kind: estimated from what its parts are like'
      : `${c.n.toLocaleString()} recorded at ${c.label}${c.h != null ? `, ${Math.round(100 * c.h / c.n)} % of them transformed` : ''}`;
    return `${j.name}: ${pct(j.p)} chance this reaction transforms it (${ev}; ${j.confidence} confidence)`;
  }

  window.Atlas = { ready, describe, judge, explain, mix32, fnv32, wl, key, features: (...a) => features(...a), combine, shard, get meta() { return META; }, band, confidenceOf };
})();
