/* Relevance — which recorded reaction template does this reaction use? The learned tier's main predictor, in the browser.

   Built by tools/rel_train.py (numpy, from the 1.09 M learned reactions), measured on patents it never saw (tools/eval_top1.py), and
   shipped as data/learned/m/ (tools/rel_export.py). This file does the same steps as the Python (tools/rel_features.py,
   tools/rel_fuse.py, tools/rel_reactivity.py), bit for bit where it matters:

     1. features     the starting materials' Morgan bits, the reagent's (Morgan bits, one bit per part, its Reagents-sheet roles,
                     its fragments), the fragments of the starting materials (found once / at least twice), a few counts
     2. the model    bits -> 1024 -> 512 -> a probability for each of the ~39,000 common templates (+ "none of them")
     3. products     the most probable templates are run on the molecules (the page's runLearned), each once as learned and once
                     with every unchanged leaf atom stripped (data/learned/m/gen3.js); the probability behind each product is
                     the sum over the templates that give it (a template's share split over its different outcomes)
     4. ranking      a logistic ranker over the products: that mass, how many templates, the product's size, what the reaction made
                     and destroyed (Morgan difference), and which functional group it consumed against how readily the reagent
                     usually takes that group (the recorded reactivity order of the reagent: data/learned/m/prior.js)
     5. confidence   the first product's chance of being right, calibrated on held-out patents

   Relevance.load()                        -> Promise (once; ~30 MB of int8 weights)
   Relevance.predict(ctx)                  -> { candidates: [{ products, score, share, mass, templates, ... }], confidence, pOther, ... }
     ctx = { R (RDKit), subs, reagent, tmpl(t) -> {id, smarts, leaving}, run(smarts) -> [{products, used}], rc: Map(template # -> recorded reactions with this exact reagent), H: page helpers }
   No molecule leaves the browser. */
window.Relevance = (() => {
  'use strict';
  const BASE = 'data/learned/m/';
  const S = { meta: null, ready: null, w1: null, w2: null, w3: [], classes: null, gen3: null, fusion: null, prior: null, roleIdx: null, fragIdx: null, solvIdx: null, colMap: null };
  const TYPED = { int8: Int8Array, uint8: Uint8Array, int16: Int16Array, int32: Int32Array, uint32: Uint32Array, float32: Float32Array, float64: Float64Array };
  const EPS = 1e-5, FLAT = JSON.stringify({ doIsomericSmiles: false });

  const loadScript = src => new Promise((res, rej) => {
    const el = document.createElement('script'); el.src = src;
    el.onload = () => { el.remove(); res(); }; el.onerror = () => { el.remove(); rej(new Error('could not load ' + src)); };
    document.head.appendChild(el);
  });
  async function inflate(b64, gz) {
    const bin = atob(b64), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    if (!gz) return u8.buffer;
    return await new Response(new Blob([u8]).stream().pipeThrough(new DecompressionStream('gzip'))).arrayBuffer();
  }
  /* a part written by tools/rel_export.py write_part: named typed arrays packed in one (gzipped) base64 blob */
  async function part(name) {
    await loadScript(BASE + name.replace('/', '_') + '.js');
    const P = (window.LEARNED_PARTS || {})['m/' + name];
    if (!P) throw new Error('data/learned/m/' + name + '.js held no data');
    delete window.LEARNED_PARTS['m/' + name];
    const buf = await inflate(P.b64, P.gz), arrays = {};
    for (const [k, a] of Object.entries(P.arrays)) { const T = TYPED[a.dtype]; arrays[k] = new T(buf, a.offset, a.bytes / T.BYTES_PER_ELEMENT); }
    P.b64 = null;
    return Object.assign({}, P, { arrays });
  }
  /* a part that is plain gzipped JSON (the fusion weights, the reactivity table) */
  async function jsonPart(name) {
    await loadScript(BASE + name + '.js');
    const P = (window.LEARNED_PARTS || {})['m/' + name];
    if (!P) throw new Error('data/learned/m/' + name + '.js held no data');
    delete window.LEARNED_PARTS['m/' + name];
    return JSON.parse(new TextDecoder().decode(await inflate(P.gz, true)));
  }

  function load(onProgress) {
    if (S.ready) return S.ready;
    const say = m => { try { onProgress && onProgress(m); } catch (e) {} };
    S.ready = (async () => {
      say('Reading the relevance model…');
      await loadScript(BASE + 'meta.js');
      const meta = S.meta = window.RELEVANCE_META;
      if (!meta) throw new Error('data/learned/m/meta.js held no data');
      S.roleIdx = new Map(meta.roles.map((r, i) => [r, i]));
      S.fragIdx = new Map(meta.fragments.map((f, i) => [f, i]));
      S.solvIdx = new Map((meta.solvents || []).map((n, i) => [n, i]));
      if (meta.input_cols) { S.colMap = new Int32Array(meta.layout.bits).fill(-1); meta.input_cols.forEach((b, j) => { S.colMap[b] = j; }); }
      S.w1 = (await part('w1')).arrays; S.w2 = (await part('w2')).arrays;
      S.w3 = [];
      for (let k = 0; k < meta.chunks.length; k++) { say(`Reading the output layer (${k + 1}/${meta.chunks.length})…`); const p = await part('w3_' + k); S.w3.push({ start: p.start, n: p.n, q: p.arrays.q, s: p.arrays.s, b: p.arrays.b }); }
      S.classes = (await part('classes')).arrays.t;
      try { const g = await jsonPart('gen3'); S.gen3 = { u: g.u, i: Int32Array.from(g.i) }; } catch (e) { S.gen3 = null; }
      try { S.fusion = await jsonPart('fusion'); } catch (e) { S.fusion = null; }
      try { S.prior = await jsonPart('prior'); } catch (e) { S.prior = null; }
      say('Ready.');
      return meta;
    })();
    S.ready.catch(() => { S.ready = null; });
    return S.ready;
  }

  /* ---- features (tools/rel_features.py) ---- */
  function fnv1a32(text) { let h = 0x811C9DC5; for (const b of new TextEncoder().encode(text)) { h ^= b; h = Math.imul(h, 0x01000193) >>> 0; } return h >>> 0; }
  const MORGAN = {};
  function morganBits(R, smiles, nBits) {
    let m = null; try { m = R.get_mol(smiles); } catch (e) { m = null; }
    if (!m) return null;
    try {
      const s = m.get_morgan_fp(MORGAN[nBits] || (MORGAN[nBits] = JSON.stringify({ radius: 2, nBits }))), out = [];
      for (let i = 0; i < s.length; i++) if (s.charCodeAt(i) === 49) out.push(i);
      return out;
    } finally { m.delete(); }
  }
  const canonical = (R, smiles) => { let m = null; try { m = R.get_mol(smiles); return m ? m.get_smiles() : null; } catch (e) { return null; } finally { if (m) m.delete(); } };
  const flatOf = (R, smiles) => { let m = null; try { m = R.get_mol(smiles); return m ? m.get_smiles(FLAT) : smiles; } catch (e) { return smiles; } finally { if (m) m.delete(); } };
  const heavyOf = (R, smiles) => { let m = null; try { m = R.get_mol(smiles); return m ? m.get_num_atoms(true) : 0; } catch (e) { return 0; } finally { if (m) m.delete(); } };     // heavy atoms only, as the Python counts

  /* the set bits of the model's input vector; H = { fragmentCounts(R, list), fragmentsOf(R, smiles), keyOf(R, smiles), roles, solvents (names from the conditions box) } */
  function features(R, subs, reagent, H) {
    const L = S.meta.layout, nf = S.meta.fragments.length, bits = [];
    const sub = new Set();
    for (const s of subs) { const b = morganBits(R, s, S.meta.fp_bits.sub); if (b) b.forEach(i => sub.add(i)); }
    sub.forEach(i => bits.push(L.sub_fp[0] + i));
    const parts = [...new Set(String(reagent || '').split('.').filter(Boolean))].sort();
    const rb = new Set(), ids = new Set(), rl = new Set(), fr = new Set();
    for (const p of parts) {
      const b = morganBits(R, p, S.meta.fp_bits.reag); if (!b) continue;
      b.forEach(i => rb.add(i));
      ids.add(fnv1a32(H.keyOf(R, p) || p) % S.meta.fp_bits.id);
      const info = (H.roles || {})[canonical(R, p) || p] || (H.roles || {})[p];
      if (info) for (const role of [info.role, info.also]) if (S.roleIdx.has(role)) rl.add(S.roleIdx.get(role));
      for (const id of H.fragmentsOf(R, p)) if (S.fragIdx.has(id)) fr.add(S.fragIdx.get(id));
    }
    rb.forEach(i => bits.push(L.reag_fp[0] + i)); ids.forEach(i => bits.push(L.reag_id[0] + i)); rl.forEach(i => bits.push(L.roles[0] + i));
    const counts = H.fragmentCounts(R, subs);
    for (const [id, f] of S.fragIdx) { const c = counts[id] || 0; if (c >= 1) bits.push(L.frag_sub[0] + f); if (c >= 2) bits.push(L.frag_sub[0] + nf + f); }
    fr.forEach(f => bits.push(L.frag_reag[0] + f));
    if (L.solv) for (const name of (H.solvents || [])) if (S.solvIdx.has(name)) bits.push(L.solv[0] + S.solvIdx.get(name));      // the recorded solvents (the conditions box)
    const m0 = L.meta[0], ns = subs.length, np = parts.length;
    bits.push(m0 + Math.min(ns, 4) - 1);                                    // 1, 2, 3, 4+ starting materials
    if (reagent) bits.push(m0 + 4);                                         // a reagent is given
    bits.push(m0 + 5 + (np === 0 ? 0 : np === 1 ? 1 : np === 2 ? 2 : np <= 4 ? 3 : 4));
    return Int32Array.from(bits);
  }

  /* ---- the network (int8 weights, one scale per output unit) ---- */
  function forward(bits) {
    const [H1, H2] = S.meta.hidden, w1 = S.w1, w2 = S.w2, acc = new Int32Array(H1), map = S.colMap;
    for (let k = 0; k < bits.length; k++) {
      const i = map ? map[bits[k]] : bits[k]; if (i < 0) continue;
      const row = i * H1;
      for (let j = 0; j < H1; j++) acc[j] += w1.q[row + j];
    }
    const h1 = new Float32Array(H1);
    for (let j = 0; j < H1; j++) { const v = acc[j] * w1.s[j] + w1.b[j]; h1[j] = v > 0 ? v : 0; }
    const acc2 = new Float32Array(H2);
    for (let k = 0; k < H1; k++) { const a = h1[k]; if (a === 0) continue; const row = k * H2; for (let j = 0; j < H2; j++) acc2[j] += a * w2.q[row + j]; }
    const h2 = new Float32Array(H2);
    for (let j = 0; j < H2; j++) { const v = acc2[j] * w2.s[j] + w2.b[j]; h2[j] = v > 0 ? v : 0; }
    const C = S.meta.classes, logits = new Float32Array(C + 1);
    for (const ch of S.w3) {
      const n = ch.n, acc3 = new Float32Array(n);
      for (let k = 0; k < H2; k++) { const a = h2[k]; if (a === 0) continue; const row = k * n; for (let c = 0; c < n; c++) acc3[c] += a * ch.q[row + c]; }
      for (let c = 0; c < n; c++) logits[ch.start + c] = acc3[c] * ch.s[c] + ch.b[c];
    }
    let mx = -Infinity; for (let c = 0; c <= C; c++) if (logits[c] > mx) mx = logits[c];
    let z = 0; for (let c = 0; c <= C; c++) { logits[c] = Math.exp(logits[c] - mx); z += logits[c]; }
    for (let c = 0; c <= C; c++) logits[c] /= z;
    return logits;                                     // probabilities; the last one is "none of the common templates"
  }
  function topClasses(p, n) {
    const C = S.meta.classes, idx = new Uint32Array(C); for (let i = 0; i < C; i++) idx[i] = i;
    // partial selection: a heap-free threshold pass is enough for n ~ 100
    const out = Array.from(idx).sort((a, b) => p[b] - p[a]).slice(0, n);
    return out.map(c => ({ c, t: S.classes[c], p: p[c] }));
  }

  /* ---- the reactivity prior (tools/rel_reactivity.py Prior) ---- */
  function priorRates(exactKey, roleSig) {
    const P = S.prior; if (!P) return null;
    const K = P.k, nf = P.fragments.length;
    const key = (exactKey || '') + '|' + (roleSig == null ? '\u0000' : roleSig);
    if (!P._cache) P._cache = new Map();
    if (P._cache.has(key)) return P._cache.get(key);
    const p = new Float64Array(nf);
    const g = P.global;
    for (let f = 0; f < nf; f++) p[f] = (g.nc[f] + 1) / (g.np[f] + 2);
    const r = roleSig != null ? P.role[roleSig] : null;
    if (r) for (let f = 0; f < nf; f++) p[f] = (r.nc[f] + K * p[f]) / (r.np[f] + K);
    const e = exactKey ? P.exact[exactKey] : null;
    if (e) for (let f = 0; f < nf; f++) p[f] = (e.nc[f] + K * p[f]) / (e.np[f] + K);
    for (let f = 0; f < nf; f++) p[f] = Math.min(1 - 1e-3, Math.max(1e-3, p[f]));
    P._cache.set(key, p);
    return p;
  }
  function priorFeatures(sCounts, cCounts, rates) {
    let a = 0, b = 0, nc = 0, nfm = 0;
    for (let f = 0; f < sCounts.length; f++) {
      const s = sCounts[f], c = cCounts[f];
      if (s > c) { nc++; if (rates) a += Math.log(rates[f] / (1 - rates[f])); }
      else if (s >= 1 && rates) b += Math.log(1 - rates[f]);
      if (c > s) nfm++;
    }
    return [a, b, nc, nfm];
  }
  const roleSignature = (reagent, roles, skip) => {
    if (!reagent) return '';
    const sig = new Set();
    for (const p of reagent.split('.')) { const r = (roles || {})[p]; if (!r) return null; if (!skip.has(r.role)) sig.add(r.role); }
    return sig.size ? [...sig].sort().join('|') : null;
  };
  const SKIP_ROLES = new Set(['counter-ion / spectator ion', 'solvent', 'salt (work-up or additive)', 'other', 'drying agent']);

  /* piecewise-linear map through the knots (x increasing), flat outside them: numpy.interp */
  function interp(v, xs, ys) {
    if (v <= xs[0]) return ys[0];
    const n = xs.length; if (v >= xs[n - 1]) return ys[n - 1];
    let lo = 0, hi = n - 1; while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (xs[mid] <= v) lo = mid; else hi = mid; }
    const t = (v - xs[lo]) / (xs[hi] - xs[lo]);
    return ys[lo] + t * (ys[hi] - ys[lo]);
  }

  /* ---- products, ranking, confidence ---- */
  async function predict(ctx, opts = {}) {
    await load();
    const R = ctx.R, H = ctx.H, topn = opts.topn || 100, F = S.fusion;
    const subs = ctx.subs.flatMap(s => s.split('.')).filter(Boolean), hasReagent = !!ctx.reagent;
    const x = features(R, subs, ctx.reagent, H);
    const prob = forward(x), pOther = prob[S.meta.classes], top = topClasses(prob, topn);
    const topP = top.length ? top[0].p : 0;
    const cands = new Map();
    let tried = 0;
    const G3 = S.gen3;                                                                 // the same templates with every unchanged leaf atom stripped (tools/template_general.py level 3)
    const II = window.Interference, marks = !!(II && II.ready && II.mapsLoaded && ctx.runMarked);        // the interference layer (js/interference.js): templates run with their changed atoms marked
    const add = (outs, src, t, p, seen, T, grp) => {
      for (const o of outs) {
        const iso = o.products.slice().sort().join('.');
        const flat = o.products.map(s => flatOf(R, s)).sort().join('.');
        let c = cands.get(flat);
        if (!c) { c = { flat, sum: 0, max: 0, n: 0, sumG3: 0, rcMax: 0, variants: new Map() }; cands.set(flat, c); }
        const w = p / Math.max(1, outs.length);
        if (src === 'mlp') { c.sum += w; c.max = Math.max(c.max, p); c.n++; if (seen > c.rcMax) c.rcMax = seen; } else c.sumG3 += w;
        const v = c.variants.get(iso) || { iso, products: o.products, sum: 0, sumG3: 0, used: o.used, usedN: 0, best: null, marked: o.marked || null, grp };      // marked: the first producer's products with the changed atoms marked; grp: its relaxed template
        if (src === 'mlp') v.sum += w; else v.sumG3 += w;
        v.usedN = Math.max(v.usedN, o.used.length);
        if (!v.best || (src === 'mlp' && (v.best.relaxed || p > v.best.p))) v.best = { t, p, id: T.id, leaving: T.leaving || '', relaxed: src !== 'mlp' };      // the specific template speaks for a product when there is one
        c.variants.set(iso, v);
      }
    };
    for (const { c: cls, t, p } of top) {
      const T = await ctx.tmpl(t); if (!T) continue;
      tried++;
      const gi = G3 ? G3.i[cls] : -1, grp = gi >= 0 ? gi : -1 - cls;                                   // candidates of one relaxed template are the different sites it could react at
      if (marks) {
        add(ctx.runMarked(II.markTemplate(T.smarts, II.mapsOf('c', cls))), 'mlp', t, p, ctx.rc ? (ctx.rc.get(t) || 0) : 0, T, grp);
        if (gi >= 0) add(ctx.runMarked(II.markTemplate(G3.u[gi], II.mapsOf('g', gi))), 'g3', t, p, 0, T, grp);
      } else {
        add(ctx.run(T.smarts), 'mlp', t, p, ctx.rc ? (ctx.rc.get(t) || 0) : 0, T, grp);      // rc: how many recorded reactions used this exact reagent with this template
        if (gi >= 0) add(ctx.run(G3.u[gi]), 'g3', t, p, 0, T, grp);
      }
    }
    const list = [...cands.values()];
    for (const c of list) { c.rep = [...c.variants.values()].sort((a, b) => (b.sum + 1e-3 * b.sumG3) - (a.sum + 1e-3 * a.sumG3))[0]; }
    if (!list.length) return { candidates: [], pOther, topP, tried, confidence: null, subHeavy: Math.max(1, subs.reduce((t, s) => t + heavyOf(R, s), 0)) };
    list.sort((a, b) => b.sum - a.sum);
    const tot = list.reduce((t, c) => t + c.sum, 0) || 1;
    // features of the candidates: what the reaction made and destroyed, size, fragments consumed
    const subBits = new Set(); for (const s of subs) { const b = morganBits(R, s, S.meta.fp_bits.sub); if (b) b.forEach(i => subBits.add(i)); }
    const subHeavy = Math.max(1, subs.reduce((t, s) => t + heavyOf(R, s), 0));
    const subCounts = S.prior ? H.fragmentCounts(R, subs) : null;
    const sCounts = S.prior ? S.prior.fragments.map(f => subCounts[f] || 0) : null;
    const rates = S.prior ? priorRates(H.exactKey, H.roleSig) : null;
    const pctx = (rates && II && II.hasCtx) ? II.ctxRates(rates, sCounts) : null;                         // the reagent's rates adjusted for the other groups present (tools/int_context.py)
    const fragIds = S.prior ? S.prior.fragments : null;
    list.slice(0, 200).forEach(c => {
      const prods = c.rep.products, pb = new Set();
      for (const s of prods) { const b = morganBits(R, s, S.meta.fp_bits.sub); if (b) b.forEach(i => pb.add(i)); }
      c.formed = [...pb].filter(i => !subBits.has(i)); c.destroyed = [...subBits].filter(i => !pb.has(i));
      c.heavy = prods.reduce((t, s) => t + heavyOf(R, s), 0); c.nMols = prods.length;
      if (S.prior) {
        const pc = H.fragmentCounts(R, prods);
        c.cCounts = S.prior.fragments.map(f => pc[f] || 0);
        c.prior = priorFeatures(sCounts, c.cCounts, rates);
        if (pctx) c.ctx = priorFeatures(sCounts, c.cCounts, pctx);
        if (marks && c.rep.marked) {                                                                     // where the template reacted, and what the learned interference weights make of that site
          const consumed = new Set(); for (let f = 0; f < sCounts.length; f++) if (sCounts[f] > c.cCounts[f]) consumed.add(fragIds[f]);
          c.family = II.familyOf(consumed); c.flags = II.flags(R, c.rep.marked); c.intRaw = II.score(c.flags, c.family);
        }
      }
    });
    if (marks) {                                                                                         // the gap to the best site of the same relaxed template, and the log of the site's share among them
      const best = new Map(), tot = new Map();
      for (const c of list) if (c.flags) best.set(c.rep.grp, Math.max(best.has(c.rep.grp) ? best.get(c.rep.grp) : -1e9, c.intRaw));
      for (const c of list) if (c.flags) tot.set(c.rep.grp, (tot.get(c.rep.grp) || 0) + Math.exp(c.intRaw - best.get(c.rep.grp)));
      for (const c of list) {
        if (c.flags) { const mx = best.get(c.rep.grp); c.intRel = c.intRaw - mx; c.intLshare = c.intRaw - (mx + Math.log(tot.get(c.rep.grp))); }
        else { c.intRel = 0; c.intLshare = 0; }
      }
    }
    list.forEach((c, i) => {
      const share = c.sum / tot, size = Math.min(3, Math.max(-3, Math.log(Math.max(1, c.heavy || subHeavy) / subHeavy))), pr = c.prior || [0, 0, 0, 0];
      const hr = hasReagent ? 1 : 0;
      c.f = { log_sum: Math.log(c.sum + EPS), log_max: Math.log(c.max + EPS), log_n: Math.log1p(c.n), log_share: Math.log(share + EPS), is_first: i === 0 ? 1 : 0, log_size: size,
        extra_mols: (c.nMols || 1) - 1, sum_x_reagent: Math.log(c.sum + EPS) * hr, sum_x_other: Math.log(c.sum + EPS) * pOther, share_x_reagent: Math.log(share + EPS) * hr,
        prior_consumed: pr[0], prior_left: pr[1], prior_both: pr[0] + pr[1], n_consumed: pr[2], n_formed: pr[3],
        log_sum_g3: Math.log(c.sumG3 + EPS), only_gen: c.sum === 0 ? 1 : 0,
        frac_used: Math.min(1, (c.rep.usedN || 0) / Math.max(1, subs.length)),
        rc_log: Math.log1p(c.rcMax || 0), rc_seen: (hasReagent && c.rcMax > 0) ? 1 : 0,
        int_raw: c.intRaw || 0, int_rel: c.intRel || 0, int_lshare: c.intLshare || 0,
        ctx_consumed: c.ctx ? c.ctx[0] : 0, ctx_left: c.ctx ? c.ctx[1] : 0, ctx_both: c.ctx ? c.ctx[0] + c.ctx[1] : 0 };
    });
    // the fused score
    if (F) {
      const names = F.dense, nd = names.length;
      for (const c of list) {
        let s = 0;
        for (let j = 0; j < nd; j++) s += ((c.f[names[j]] || 0) - F.mu[j]) / F.sd[j] * F.w_dense[j];
        if (c.formed) { for (const b of c.formed) s += F.w_formed[b]; for (const b of c.destroyed) s += F.w_destroyed[b]; }
        c.score = s;
      }
      list.sort((a, b) => b.score - a.score);
    } else list.forEach(c => { c.score = Math.log(c.sum + EPS); });
    const mx = list[0].score; let z = 0; list.forEach(c => { c.e = Math.exp(c.score - mx); z += c.e; });
    list.forEach(c => { c.share = c.e / z; });
    let confidence = null;
    if (F && F.confidence) {
      const p1 = list[0].share, p2 = list.length > 1 ? list[1].share : 0;
      // the confidence model's inputs, by name (data/learned/m/fusion.js lists the names it was fitted with; tools/rel_fuse.py confidence_features)
      const have = {
        log_share: Math.log(p1 + 1e-4), log_margin: Math.log(p1 - p2 + 1e-4), log_mass: Math.log(list[0].sum + 1e-5), p_other: pOther, has_reagent: hasReagent ? 1 : 0,
        log_n_candidates: Math.log1p(list.length), top_template_p: topP,
        reagent_never_with_template: (hasReagent && !(list[0].rcMax > 0)) ? 1 : 0,       // the leading templates were never recorded with this exact reagent
        log_sub_heavy: Math.log(subHeavy) };
      const names = F.confidence_features || Object.keys(have), v = names.map(n => have[n] == null ? 0 : have[n]);
      const m = F.confidence; let zz = m.w[m.w.length - 1];
      for (let j = 0; j < v.length; j++) zz += (v[j] - m.mu[j]) / m.sd[j] * m.w[j];
      confidence = 1 / (1 + Math.exp(-zz));
      if (m.iso) confidence = interp(confidence, m.iso.x, m.iso.y);        // the monotone recalibration measured on held-out patents (tools/rel_fuse.py isotonic_knots)
    }
    return { candidates: list, confidence, pOther, topP, tried, features: x.length, subHeavy };
  }

  return { load, features, forward, topClasses, predict, canonical, heavyOf, roleSignature, SKIP_ROLES, get ready() { return !!(S.meta && S.w3.length); }, get meta() { return S.meta; }, get state() { return S; } };
})();
