/* Interference — everything that hinders (or helps) a template, tested on a candidate product.

   The knowledge is the sheet "Interferences" of reaction_rules.xlsx (exported by tools/rules_to_js.py into window.REACTION_RULES.interferences):
     SITE patterns   a SMARTS whose first atom is the atom the template changes: an ortho nitro group on the ring carbon that reacts, a branch on the
                     carbon that is attacked ... A learned template is run with its changed atoms marked (isotope 900 on the product side of the template,
                     markTemplate below), so every product says where the reaction happened, and the patterns are tested at the marked atoms.
     ADVISORIES      a reagent (a role, a fragment or a part as typed) and a group that is present: a Grignard and an alcohol, LiAlH4 and an ester ...
   How much each site pattern matters, and in which direction, is learned from the recorded reactions (tools/int_contests.py, tools/int_train.py) and shipped
   as data/learned/i/weights.js; the ranking in js/relevance.js uses the resulting score and this file explains it. Same logic as the Python
   (tools/int_marks.py, tools/interferences.py, tools/interference_defs.py); tools/int_check_js.py compares the two.

   Interference.markTemplate(smarts, maps)    the template with an isotope on the product atoms the template changes (maps = their atom-map numbers)
   Interference.flags(R, markedSmilesList)    keys of the site patterns that hold at a marked atom of the products
   Interference.familyOf(consumedIds)         the reaction family (the learned weights are kept per family)
   Interference.score(flags, family)          the learned log-score of a site (sum of the weights of its flags + the Hammett terms)
   Interference.explain(flags, family)        the flags with their learned weight and the sheet's text, strongest first
   Interference.advisories(ctx)               the reagent + group advisories that hold; ctx = { roles: Set, frags: Set, parts: Set, solventTypes: Set, present: Set, consumed: Set }
   No molecule leaves the browser. */
window.Interference = (() => {
  'use strict';
  const MARK = 900, BASE = 'data/learned/i/';
  const STRICT = /;(?:D\d|H\d|[+-]\d?)/;
  const S = { rows: null, site: null, adv: null, W: null, ctx: null, ready: null, q: new Map(), pre: new Map(), cache: new Map(), maps: null, mapsReady: null };

  /* ---- marking a template (tools/int_marks.py mark_product_side / mark_template) ---- */
  function markProductSide(p, markNew, only) {
    const out = [], n = p.length;
    let i = 0;
    while (i < n) {
      const ch = p[i];
      if (ch === '[') {
        let depth = 1, j = i + 1;
        while (j < n && depth) { if (p[j] === '[') depth++; else if (p[j] === ']') depth--; j++; }
        let body = p.slice(i + 1, j - 1);
        const mm = body.match(/:(\d+)$/);
        const changed = mm ? (only ? only.has(parseInt(mm[1], 10)) : STRICT.test(body)) : markNew;
        if (changed && !/^\d/.test(body)) body = MARK + body;
        out.push('[' + body + ']'); i = j;
      } else {
        const two = p.substr(i, 2);
        if (two === 'Cl' || two === 'Br' || ('BCNOPSFI'.includes(ch) && !(i + 1 < n && p[i + 1] === 'l' && ch === 'C')) || 'bcnops*'.includes(ch)) {
          const tok = (two === 'Cl' || two === 'Br') ? two : ch;
          out.push(markNew ? '[' + MARK + tok + ']' : tok); i += tok.length;
        } else { out.push(ch); i++; }
      }
    }
    return out.join('');
  }
  function markTemplate(smarts, maps, markNew = true) {
    const k = smarts.indexOf('>>');
    if (k < 0) return smarts;
    const only = maps == null ? null : (maps instanceof Set ? maps : new Set(maps));
    return smarts.slice(0, k) + '>>' + markProductSide(smarts.slice(k + 2), markNew, only);
  }
  const stripMarks = smi => smi.replace(new RegExp('\\[' + MARK, 'g'), '[');

  /* ---- loading ---- */
  const loadScript = src => new Promise((res, rej) => {
    const el = document.createElement('script'); el.src = src;
    el.onload = () => { el.remove(); res(); }; el.onerror = () => { el.remove(); rej(new Error('could not load ' + src)); };
    document.head.appendChild(el);
  });
  async function gunzipJson(b64) {
    const bin = atob(b64), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    return JSON.parse(await new Response(new Blob([u8]).stream().pipeThrough(new DecompressionStream('gzip'))).text());
  }
  async function part(name) {
    await loadScript(BASE + name + '.js');
    const P = (window.LEARNED_PARTS || {})['i/' + name];
    if (!P) throw new Error('data/learned/i/' + name + '.js held no data');
    delete window.LEARNED_PARTS['i/' + name];
    return gunzipJson(P.gz);
  }
  function rulesRows() { return (window.REACTION_RULES || {}).interferences || null; }
  /* the sheet's rows (always) and the learned weights (when installed) */
  function load() {
    if (S.ready) return S.ready;
    S.ready = (async () => {
      const rows = rulesRows();
      if (!rows) throw new Error('no Interferences sheet in data/reaction_rules.js');
      S.rows = rows;
      S.site = rows.filter(r => r.smarts);
      S.adv = rows.filter(r => !r.smarts && r.groups && r.groups.length && r.page && !['size', 'reagent', 'conditions'].includes(r.needs));
      S.byKey = new Map(rows.map(r => [r.key, r]));
      try { S.W = await part('weights'); } catch (e) { S.W = null; }
      try { S.ctx = await part('ctx'); } catch (e) { S.ctx = null; }
      if (S.W) { S.famIdx = new Map(S.W.families.map((f, i) => [f, i])); S.keyIdx = new Map(S.W.keys.map((k, i) => [k, i])); }
      return true;
    })();
    S.ready.catch(() => { S.ready = null; });
    return S.ready;
  }
  /* the changed atom-map numbers of every class template ('c') and of its relaxed form ('g'): data/learned/i/maps.js (tools/int_export.py) */
  function loadMaps() {
    if (S.mapsReady) return S.mapsReady;
    S.mapsReady = (async () => { S.maps = await part('maps'); return S.maps; })();
    S.mapsReady.catch(() => { S.mapsReady = null; });
    return S.mapsReady;
  }
  function mapsOf(kind, idx) {                     // (a template RDKit could not read has no changed atom: the empty set, as the Python gives)
    const m = S.maps && S.maps[kind] && S.maps[kind][idx];
    return new Set(m == null ? [] : String(m).split(',').filter(Boolean).map(Number));
  }

  /* ---- site patterns ---- */
  const markFirst = sm => '[' + MARK + sm.slice(1);                                       // [c;$(c:c-X)] -> [900c;$(c:c-X)]
  const ANCHOR_PRE = { 'aromatic C': '[900c]', 'aromatic n': '[900n]', 'sp3 C': '[900C;X4]', 'carbonyl C': '[900C;X3;$(C=[O,N,S])]', 'alkene C': '[900C;X3;$(C=C)]', 'nitrogen': '[900N]', 'oxygen': '[900O]', 'sulfur': '[900S]' };
  function qmol(R, sm) {
    if (!S.q.has(sm)) { let q = null; try { q = R.get_qmol(sm); } catch (e) { q = null; } S.q.set(sm, q); }
    return S.q.get(sm);
  }
  const matches = (m, q) => { if (!q) return false; try { return Object.keys(JSON.parse(m.get_substruct_match(q) || '{}')).length > 0; } catch (e) { return false; } };
  function flagsOfMol(R, m) {
    const out = [], okAnchor = new Map(), hasSub = new Map();
    for (const r of S.site) {
      const pre = ANCHOR_PRE[r.anchor];
      if (pre) {
        if (!okAnchor.has(r.anchor)) okAnchor.set(r.anchor, matches(m, qmol(R, pre)));
        if (!okAnchor.get(r.anchor)) continue;
      }
      if (r.prefilter) {                                                                  // the group a ring-substituent pattern needs must be in the molecule at all
        if (!hasSub.has(r.prefilter)) hasSub.set(r.prefilter, matches(m, qmol(R, r.prefilter)));
        if (!hasSub.get(r.prefilter)) continue;
      }
      if (matches(m, qmol(R, markFirst(r.smarts)))) out.push(r.key);
    }
    return out;
  }
  /* the keys of the site patterns that hold at a marked atom of some product (the marked product SMILES of one candidate) */
  function flags(R, markedList) {
    if (!S.site) return [];
    const got = new Set();
    for (const smi of markedList) {
      for (let part of String(smi).split('.')) {
        if (!part.includes(String(MARK))) continue;
        part = part.replace(/\[[23]H\]/g, '[H]');                                         // deuterium and tritium are atoms to RDKit and would count in D / X: read them as hydrogen
        if (S.cache.has(part)) { S.cache.get(part).forEach(k => got.add(k)); continue; }
        let m = null, fl = [];
        try { m = R.get_mol(part); if (m) fl = flagsOfMol(R, m); } catch (e) { fl = []; } finally { if (m) m.delete(); }
        if (S.cache.size > 50000) S.cache.clear();
        S.cache.set(part, fl); fl.forEach(k => got.add(k));
      }
    }
    return [...got].sort();
  }

  /* ---- families and learned weights ---- */
  function familyOf(consumed) {
    const F = (S.W && S.W.family_fragments) || [];
    for (const [name, ids] of F) if (ids.some(i => consumed.has(i))) return name;
    return 'other';
  }
  /* sum of the weights of the flags (global + family) and the Hammett terms: rho(position) x sum of sigma at that position */
  function score(fl, family) {
    const W = S.W; if (!W) return 0;
    const fi = S.famIdx.get(family); const fam = fi == null ? null : W.fam[family];
    let s = 0; const sg = [0, 0, 0];
    for (const k of fl) {
      const i = S.keyIdx.get(k); if (i == null) continue;
      s += W.g[i] + (fam ? fam[i] : 0);
      const t = W.sigma[k]; if (t) sg[t[0]] += t[1];
    }
    for (let p = 0; p < 3; p++) s += (W.rho_g[p] + (fam ? W.rho_f[family][p] : 0)) * sg[p];
    return s;
  }
  function explain(fl, family) {
    const W = S.W; if (!W) return [];
    const fam = W.fam[family];
    return fl.map(k => {
      const i = S.keyIdx.get(k), row = S.byKey.get(k) || {};
      let w = i == null ? 0 : W.g[i] + (fam ? fam[i] : 0);
      const t = W.sigma[k]; if (t) w += (W.rho_g[t[0]] + (fam ? W.rho_f[family][t[0]] : 0)) * t[1];             // a ring substituent also acts through its Hammett constant
      return { key: k, name: row.name || k, weight: w, effect: row.effect || '', why: row.why || '', kind: row.kind || '', anchor: row.anchor || '' };
    }).sort((a, b) => Math.abs(b.weight) - Math.abs(a.weight));
  }

  /* ---- chemoselectivity from the other groups present (tools/int_context.py): the reagent's consumption rate of each group, adjusted for its competitors ---- */
  function ctxRates(rates, sCounts) {
    const C = S.ctx; if (!C) return null;
    const nf = sCounts.length, p = new Float64Array(nf), present = [];
    for (let g = 0; g < nf; g++) if (sCounts[g] >= 1) present.push(g);
    for (let f = 0; f < nf; f++) {
      if (sCounts[f] < 1) { p[f] = rates[f]; continue; }
      let z = Math.log(rates[f] / (1 - rates[f])) + C.b[f] + (sCounts[f] >= 2 ? C.mult[f] : 0);
      const row = C.beta[f];
      if (row) for (const [g, v] of row) if (sCounts[g] >= 1) z += v;
      p[f] = Math.min(1 - 1e-3, Math.max(1e-3, 1 / (1 + Math.exp(-z))));
    }
    return p;
  }

  /* ---- advisories: a reagent and a group that is present ---- */
  const ELEM = new Map();
  function hasElement(parts, sym) {                                                      // some reagent part, as typed, has this element as a bracket atom ([Pd], [Pd+2] ...)
    if (!ELEM.has(sym)) ELEM.set(sym, new RegExp('\\[(?:\\d+)?' + sym.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '(?![a-z])'));
    const re = ELEM.get(sym);
    for (const p of parts) if (re.test(p)) return true;
    return false;
  }
  function tokenHolds(t, ctx) {
    if (t.startsWith('not:')) return !tokenHolds(t.slice(4), ctx);                       // "not:smiles:[Ni]": the reagent must NOT have this part
    const k = t.indexOf(':'), kind = t.slice(0, k), val = t.slice(k + 1);
    if (kind === 'role') return ctx.roles.has(val);
    if (kind === 'frag') return ctx.frags.has(val);
    if (kind === 'smiles') return ctx.parts.has(val);
    if (kind === 'elem') return hasElement(ctx.parts, val);
    if (kind === 'solvent') return !!ctx.solventTypes && ctx.solventTypes.has(val);
    return false;
  }
  /* ctx = { roles, frags (found in the reagent parts), parts (as typed), solventTypes, present (fragments of the starting materials), consumed (the ones the leading product consumed) } */
  function advisories(ctx) {
    if (!S.adv) return [];
    const out = [];
    for (const a of S.adv) {
      if (a.needs === 'solvent' && !(ctx.solventTypes && ctx.solventTypes.size)) continue;
      if (a.reagent.length && !a.reagent.some(alt => alt.every(t => tokenHolds(t, ctx)))) continue;      // no reagent side: any reagent (or none)
      const have = a.groupWhere === 'reagent' ? ctx.frags : a.groupWhere === 'both' ? new Set([...ctx.present, ...ctx.frags]) : ctx.present;
      const hit = a.groups.filter(g => have.has(g));
      if (hit.length < (a.minDistinct || 1)) continue;
      const by = hit.filter(g => !ctx.consumed.has(g)), tg = hit.filter(g => ctx.consumed.has(g));
      if (!(a.show === 'bystander' ? by : a.show === 'target' ? tg : hit).length) continue;
      out.push({ key: a.key, name: a.name, kind: a.kind, effectKind: a.effect_kind, effect: a.effect, why: a.why, groups: by, target: tg });
    }
    return out;
  }

  return { load, loadMaps, mapsOf, markTemplate, stripMarks, flags, familyOf, score, explain, advisories, ctxRates, MARK,
           get ready() { return !!S.site && !!S.W; }, get hasWeights() { return !!S.W; }, get hasCtx() { return !!S.ctx; }, get mapsLoaded() { return !!S.maps; }, get state() { return S; } };
})();
