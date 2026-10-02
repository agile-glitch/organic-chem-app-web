/* RXNKB — reaction prediction from the nearest recorded precedents, in the browser.

   The method (the same one tools/rxnkb.py measures on held-out reactions before anything is shipped):
     1. fingerprint your reactants (Morgan radius 2, 512 bits — RDKit.js, identical to the Python build);
     2. find the K recorded reactions whose reactants are most similar (Tanimoto);
     3. apply each of those reactions' TEMPLATES (the atoms that changed + their neighbours) to your molecules;
     4. rank the products by the similarity of the precedents behind each one.
   Nothing is predicted without precedent, every answer carries the reactions it came from, and the
   confidence shown is the measured hit rate of predictions that looked like this one on held-out data.

   Data: data/rxnkb/ (built by tools/rxnkb.py pack) — meta.js, templates_*.js, fp_*.js (fingerprints as base64:
   a double-clicked page cannot read binary files), prec_*.js (the precedents' details, read on demand).
   RDKit.js: vendor/rdkit/ (BSD-3), started by js/rdkit_load.js (the WebAssembly comes from RDKit_minimal.wasm.js
   for the same reason). */
window.RXNKB = (() => {
  'use strict';
  const S = { loading: null, ready: false, R: null, meta: null, T: [], TN: [], FP: null, POP: null, TID: null, AG: null, N: 0, rx: new Map(), prec: new Map() };
  const BASE = 'data/rxnkb/';

  const loadScript = src => new Promise((res, rej) => {
    const s = document.createElement('script'); s.src = src;
    s.onload = () => { s.remove(); res(); }; s.onerror = () => { s.remove(); rej(new Error('could not load ' + src)); };
    document.head.appendChild(s);
  });
  function bytesOf(b64) { const bin = atob(b64), u8 = new Uint8Array(bin.length); for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i); return u8; }
  /* text in the pack (templates, precedent details) is gzip inside base64; the browser gunzips natively */
  async function gunzipJson(b64) { const stream = new Blob([bytesOf(b64)]).stream().pipeThrough(new DecompressionStream('gzip')); return JSON.parse(await new Response(stream).text()); }
  async function gunzipBytes(b64) { const stream = new Blob([bytesOf(b64)]).stream().pipeThrough(new DecompressionStream('gzip')); return new Uint8Array(await new Response(stream).arrayBuffer()); }
  const pop32 = x => { x = x - ((x >>> 1) & 0x55555555); x = (x & 0x33333333) + ((x >>> 2) & 0x33333333); return Math.imul((x + (x >>> 4)) & 0x0F0F0F0F, 0x01010101) >>> 24; };

  /* ---- loading: once per session, with progress ---- */
  function load(onProgress) {
    if (S.ready) return Promise.resolve(S.meta);
    if (S.loading) return S.loading;
    const say = (msg, frac) => { try { onProgress && onProgress(msg, frac); } catch (e) {} };
    S.loading = (async () => {
      say('Reading the knowledge-base index…', 0);
      await loadScript(BASE + 'meta.js');
      const meta = window.RXNKB_META; if (!meta) throw new Error('data/rxnkb/meta.js held no data');
      say('Starting RDKit…', 0.03);
      /* one RDKit instance per page, shared with the Molecule tab (js/rdkit_load.js starts it from the embedded bytes) */
      S.R = await window.RDKitLoad();
      for (let k = 0; k < meta.template_files.length; k++) {
        say('Loading reaction templates…', 0.06 + 0.14 * k / meta.template_files.length);
        await loadScript(BASE + meta.template_files[k]);
        const part = window.RXNKB_T.pop(), body = part.gz ? await gunzipJson(part.gz) : part;
        for (let j = 0; j < body.smarts.length; j++) { S.T[part.start + j] = body.smarts[j]; S.TN[part.start + j] = body.n[j]; }
      }
      const N = meta.n, W = meta.fp_bits >>> 5;
      S.N = N; S.W = W; S.FP = new Uint32Array(N * W); S.AG = new Uint32Array(N * 2); S.TID = new Uint32Array(N); S.POP = new Uint16Array(N);
      for (let k = 0; k < meta.fp_files.length; k++) {
        say(`Loading ${N.toLocaleString()} recorded reactions…`, 0.2 + 0.78 * k / meta.fp_files.length);
        await loadScript(BASE + meta.fp_files[k]);
        const c = window.RXNKB_FP.pop(), raw = c.z ? gunzipBytes : async b => bytesOf(b);
        const [fp, ag, tid] = await Promise.all([raw(c.fp), raw(c.ag), raw(c.tid)]);
        S.FP.set(new Uint32Array(fp.buffer, fp.byteOffset, fp.byteLength >>> 2), c.start * W);
        S.AG.set(new Uint32Array(ag.buffer, ag.byteOffset, ag.byteLength >>> 2), c.start * 2);
        S.TID.set(new Uint32Array(tid.buffer, tid.byteOffset, tid.byteLength >>> 2), c.start);
      }
      for (let i = 0, o = 0; i < N; i++) { let p = 0; for (let w = 0; w < W; w++, o++) { const x = S.FP[o]; if (x) p += pop32(x); } S.POP[i] = p; }
      S.meta = meta; S.ready = true; say('Ready.', 1);
      return meta;
    })();
    S.loading.catch(() => { S.loading = null; });
    return S.loading;
  }

  /* ---- molecules ---- */
  const FP_OPTS = JSON.stringify({ radius: 2, nBits: 512 }), FLAT = JSON.stringify({ doIsomericSmiles: false });
  function mol(smiles) { let m = null; try { m = S.R.get_mol(smiles); } catch (e) { m = null; } if (m && m.is_valid && !m.is_valid()) { m.delete(); m = null; } return m; }
  function canonical(smiles) { const m = mol(smiles); if (!m) return null; const out = { iso: m.get_smiles(), flat: m.get_smiles(FLAT) }; m.delete(); return out; }
  function fingerprint(smiles) {
    const m = mol(smiles); if (!m) return null;
    const bits = m.get_morgan_fp(FP_OPTS); m.delete();
    /* always the 512-bit fingerprint (verified identical to the Python build); a 256-bit pack is that folded in half */
    const w = new Uint32Array(S.W), nb = S.W * 32;
    for (let i = 0; i < bits.length; i++) if (bits.charCodeAt(i) === 49) { const j = i % nb; w[j >>> 5] |= (1 << (j & 31)); }
    let pop = 0; for (let k = 0; k < S.W; k++) pop += pop32(w[k]);
    return { w, pop };
  }
  /* agents are hashed by their canonical SMILES (FNV-1a, as in tools/rxnkb.py) into 64 bits */
  function agentBits(smilesList) {
    const v = new Uint32Array(2);
    for (const raw of smilesList || []) for (const part of String(raw).split('.')) {
      const c = canonical(part); if (!c) continue;
      let h = 0x811C9DC5; for (let i = 0; i < c.iso.length; i++) { h ^= c.iso.charCodeAt(i); h = Math.imul(h, 0x01000193) >>> 0; }
      const b = h % 64; v[b >>> 5] |= (1 << (b & 31));
    }
    return v;
  }

  /* ---- the K nearest recorded reactions ---- */
  function nearest(q, K) {
    const { FP, POP, N, W } = S, qp = q.pop, qw = q.w;
    const idx = new Int32Array(K), sim = new Float32Array(K); let filled = 0, worst = 0;
    for (let i = 0, o = 0; i < N; i++, o += W) {
      const p = POP[i];
      if (filled === K) { const ub = p < qp ? p / qp : qp / p; if (ub <= worst) continue; }      // cannot beat the current K-th
      let c = 0; for (let w = 0; w < W; w++) { const x = FP[o + w] & qw[w]; if (x) c += pop32(x); }
      const s = c / (p + qp - c);
      if (filled < K) { let j = filled++; while (j > 0 && sim[j - 1] < s) { sim[j] = sim[j - 1]; idx[j] = idx[j - 1]; j--; } sim[j] = s; idx[j] = i; if (filled === K) worst = sim[K - 1]; }
      else if (s > worst) { let j = K - 1; while (j > 0 && sim[j - 1] < s) { sim[j] = sim[j - 1]; idx[j] = idx[j - 1]; j--; } sim[j] = s; idx[j] = i; worst = sim[K - 1]; }
    }
    return { idx: idx.subarray(0, filled), sim: sim.subarray(0, filled) };
  }

  /* ---- applying a template to the query molecules (every assignment of molecules to its reactant patterns) ---- */
  function permutations(n, k) { const out = []; const go = (cur, used) => { if (cur.length === k) { out.push(cur.slice()); return; } for (let i = 0; i < n; i++) if (!used[i]) { used[i] = true; cur.push(i); go(cur, used); cur.pop(); used[i] = false; } }; go([], []); return out; }
  function reaction(tid) {
    let rx = S.rx.get(tid);
    if (rx === undefined) {
      try { rx = S.R.get_rxn(S.T[tid]); } catch (e) { rx = null; }
      if (S.rx.size >= 3000) { const first = S.rx.keys().next().value; const old = S.rx.get(first); if (old) old.delete(); S.rx.delete(first); }
      S.rx.set(tid, rx || null);
    }
    return rx || null;
  }
  function applyTemplate(tid, mols, reactantFlats) {
    const rx = reaction(tid); if (!rx) return [];
    const nT = S.T[tid].split('>>')[0].split('.').length;
    if (nT > mols.length) return [];
    const found = new Map();
    for (const perm of permutations(mols.length, nT)) {
      const ml = new S.R.MolList(); perm.forEach(i => ml.append(mols[i]));
      let res = null; try { res = rx.run_reactants(ml, 30); } catch (e) { res = null; }
      ml.delete(); if (!res) continue;
      const n = res.size();
      for (let a = 0; a < n; a++) {
        const pl = res.get(a), sz = pl.size();
        for (let b = 0; b < sz; b++) { const pm = pl.at(b); let smi = null; try { smi = pm.get_smiles(); } catch (e) {} pm.delete();
          const c = smi && canonical(smi); if (c && !reactantFlats.has(c.flat) && !found.has(c.flat)) found.set(c.flat, c.iso); }
        pl.delete();
      }
      res.delete();
    }
    return [...found].map(([flat, iso]) => ({ flat, iso }));
  }

  /* ---- confidence: the measured hit rate of first answers that looked like this one ---- */
  function confidence(bestSim, share) {
    const m = S.meta, bin = (edges, v) => { let k = 0; for (let j = 0; j < edges.length - 1; j++) if (v >= edges[j]) k = j; return k; };
    const i = bin(m.sim_edges, bestSim), k = bin(m.share_edges, share), cell = m.calibration[i][k];
    if (cell && cell[0] != null) return { p: cell[0], n: cell[1] };
    let ok = 0, n = 0; m.calibration[i].forEach(c => { if (c && c[0] != null) { ok += c[0] * c[1]; n += c[1]; } });   // fall back to the whole similarity row
    return n ? { p: ok / n, n } : { p: null, n: 0 };               // no measurement: the page then does not rely on it
  }

  /* ---- predict ---- */
  function predict(query, opts = {}) {
    if (!S.ready) throw new Error('RXNKB is not loaded');
    const m = S.meta, K = opts.K || m.K, alpha = m.alpha, beta = opts.agents === false ? 0 : m.beta, t0 = performance.now();
    const smiles = (query.reactants || []).map(s => String(s).trim()).filter(Boolean);
    const mols = [], flats = new Set(), isos = [];
    for (const s of smiles) for (const part of s.split('.')) { const x = mol(part); if (!x) { mols.forEach(y => y.delete()); return { error: 'RDKit could not read "' + part + '"' }; } mols.push(x); const c = { iso: x.get_smiles(), flat: x.get_smiles(FLAT) }; flats.add(c.flat); isos.push(c.iso); }
    if (!mols.length) return { error: 'no reactants' };
    if (mols.length > 4) { mols.forEach(y => y.delete()); return { error: 'at most four reactant molecules' }; }
    const q = fingerprint(isos.slice().sort().join('.'));
    const qa = agentBits(query.agents), hasAgents = !!(qa[0] | qa[1]);
    const nb = nearest(q, K), tScan = performance.now();
    const byTid = new Map(), cand = new Map();
    for (let j = 0; j < nb.idx.length; j++) {
      const i = nb.idx[j], s = nb.sim[j], tid = S.TID[i];
      if (!byTid.has(tid)) byTid.set(tid, applyTemplate(tid, mols, flats));
      const prods = byTid.get(tid); if (!prods.length) continue;
      let w = Math.pow(s, alpha);
      if (beta && hasAgents) { const a0 = S.AG[2 * i], a1 = S.AG[2 * i + 1]; const inter = pop32(a0 & qa[0]) + pop32(a1 & qa[1]), uni = pop32(a0 | qa[0]) + pop32(a1 | qa[1]); w *= 1 + beta * (uni ? inter / uni : 0); }
      if (m.gamma) w *= Math.pow(1 + Math.log10(Math.max(1, S.TN[tid] || 1)), m.gamma);          // how common the transformation is overall
      if (m.split) w /= prods.length;
      for (const p of prods) {
        let c = cand.get(p.flat); if (!c) { c = { flat: p.flat, smiles: p.iso, score: 0, best: 0, ws: [], support: [], templates: new Set() }; cand.set(p.flat, c); }
        c.ws.push(w); if (s > c.best) c.best = s;
        c.support.push({ i, sim: s, tid }); c.templates.add(tid);
      }
    }
    mols.forEach(y => y.delete());
    /* as in tools/rxnkb.py rank(): max = the best precedent; topM = the mean of the best M (a missing one counts as 0,
       so a single precedent cannot decide on its own); sum = all of them */
    const topM = m.agg === 'max' ? 1 : /^top\d+$/.test(m.agg) ? +m.agg.slice(3) : 0;
    for (const c of cand.values()) {
      if (!topM) c.score = c.ws.reduce((t, w) => t + w, 0);
      else { const top = c.ws.sort((a, b) => b - a).slice(0, topM); c.score = top.reduce((t, w) => t + w, 0) / topM; }
      delete c.ws;
    }
    const list = [...cand.values()].sort((a, b) => b.score - a.score);
    const total = list.reduce((t, c) => t + c.score, 0) || 1;
    list.forEach(c => { c.share = c.score / total; c.confidence = confidence(c.best, c.share); c.templates = [...c.templates]; c.templateCount = c.templates.reduce((t, id) => t + (S.TN[id] || 0), 0); });
    return { candidates: list, neighbours: Array.from(nb.idx, (i, j) => ({ i, sim: nb.sim[j], tid: S.TID[i] })), bestSim: nb.sim.length ? nb.sim[0] : 0,
      reactants: isos, usedAgents: hasAgents, ms: { scan: tScan - t0, total: performance.now() - t0 } };
  }

  /* ---- the precedents' details, read on demand ---- */
  async function evidence(indices) {
    const m = S.meta, out = [];
    for (const i of indices) {
      const c = Math.floor(i / m.prec_chunk);
      if (!S.prec.has(c)) {
        await loadScript(BASE + 'prec_' + String(c).padStart(4, '0') + '.js');
        const raw = window.RXNKB_PREC[c]; delete window.RXNKB_PREC[c];
        S.prec.set(c, raw.gz ? Object.assign(await gunzipJson(raw.gz), { start: raw.start }) : raw);
        if (S.prec.size > 80) S.prec.delete(S.prec.keys().next().value);
      }
      const ch = S.prec.get(c), j = i - ch.start;
      out.push({ i, reactants: ch.r[j], product: ch.p[j], agents: ch.a[j] ? ch.a[j].split('.') : [], cite: ch.c[j], year: ch.y[j], yield: ch.yl[j] === '' ? null : +ch.yl[j], tid: S.TID[i] });
    }
    return out;
  }

  return { load, predict, evidence, canonical, get ready() { return S.ready; }, get meta() { return S.meta; }, get rdkit() { return S.R; }, template: tid => ({ smarts: S.T[tid], n: S.TN[tid] }) };
})();
