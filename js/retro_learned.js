/* RetroLearned — retrosynthesis from recorded reactions.

   The taught disconnections of js/chem.js are textbook reactions. This adds what chemists actually recorded: the 39,000 most common reaction
   types of the 1.09 million patent reactions behind the Reactions page (data/learned/t/: each a reaction template, the changed atoms and their
   neighbours as SMARTS). A template read backwards says "a molecule like this can be made from precursors like those". All of them are run
   on the target; a disconnection is ranked by how often that kind of reaction was recorded, how much simpler its precursors are, and how many
   of them can be bought. For each one the page also has what was recorded for that kind of reaction (data/learned/c/: reagent lines, solvents,
   yield) and asks the Fragment Atlas (js/atlas.js) what else in the precursors those reagents would likely transform.

   Limits, stated where they matter: templates are matched on the target's structure alone, so a proposal is a candidate, not a proof that the
   forward reaction works on these molecules (check it in the Reactions tab: "Open in Reactions"); the popularity ranking favours routine
   chemistry (protections, hydrolyses) over clever chemistry; patent records carry no temperature or time.

   RetroLearned.disconnect(R, smiles, opts)  -> [{precursors:[{smiles, buyable, heavy}], recorded, templates:[ids], score, ...}]   one step, best first
   RetroLearned.route(R, smiles, opts)       -> a tree in the shape js/retro.js renders: {g, transform:{name, reagent, note}, children} | {g, basic:true}
   RetroLearned.mount(host, smiles)          -> adds the section to the Retrosynthesis page                                                              */
window.RetroLearned = (() => {
  'use strict';
  const C = window.Chem;
  const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const el = (tag, cls, html) => { const n = document.createElement(tag); if (cls) n.className = cls; if (html != null) n.innerHTML = html; return n; };
  const shortName = n => String(n || '').split(' (')[0];
  const fmtN = n => (n == null ? '' : n.toLocaleString());

  // ---------------------------------------------------------------- data (the same part files as the Reactions page: data/learned/)
  const loadScript = src => new Promise((res, rej) => {
    const s = document.createElement('script'); s.src = src;
    s.onload = () => { s.remove(); res(); }; s.onerror = () => { s.remove(); rej(new Error('could not load ' + src)); };
    document.head.appendChild(s);
  });
  async function gunzipJson(b64) {
    const bin = atob(b64), u8 = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
    return JSON.parse(await new Response(new Blob([u8]).stream().pipeThrough(new DecompressionStream('gzip'))).text());
  }
  const parts = new Map();
  function part(name) {                                       // 't/003', 'c/012'
    if (!parts.has(name)) parts.set(name, (async () => {
      try { await loadScript('data/learned/' + name + '.js'); } catch (e) { return null; }
      const P = (window.LEARNED_PARTS || {})[name];
      if (!P) return null;
      delete window.LEARNED_PARTS[name];
      try { return await gunzipJson(P.gz); } catch (e) { parts.delete(name); return null; }
    })());
    return parts.get(name);
  }
  let TEMPLATES = null;                                       // {id[], smarts[], recorded[]}: commonest first
  async function templates(n) {
    const per = 5000, files = Math.ceil(n / per);
    if (TEMPLATES && TEMPLATES.loaded >= files) return TEMPLATES;
    const T = { id: [], smarts: [], recorded: [], loaded: 0 };
    for (let k = 0; k < files; k++) {
      const b = await part('t/' + String(k).padStart(3, '0'));
      if (!b) break;
      T.id.push(...b.id); T.smarts.push(...b.smarts); T.recorded.push(...b.recorded); T.loaded = k + 1;
    }
    TEMPLATES = T;
    return T;
  }
  const tNumber = id => parseInt(String(id).slice(1), 10) - 1;
  async function recipesOf(id) {
    const t = tNumber(id), tab = await part('c/' + String(Math.floor(t / 2000)).padStart(3, '0'));
    return tab && tab[String(t)] ? tab[String(t)] : null;     // [n, median yield, recipes, solvents, consumed ids, formed ids]
  }

  // ---------------------------------------------------------------- one step backwards
  const breathe = () => new Promise(r => setTimeout(r, 0));
  const rxnCache = new Map();
  function reversed(R, T, i) {
    if (rxnCache.has(i)) return rxnCache.get(i);
    const [l, r] = T.smarts[i].split('>>');
    let rx = null;
    if (r && l && !r.includes('.')) { try { rx = R.get_rxn(r + '>>' + l); } catch (e) { rx = null; } }     // a template with two products cannot be run on one target
    rxnCache.set(i, rx);
    return rx;
  }
  const heavyOf = (R, smi) => { let m = null; try { m = R.get_mol(smi); return m ? m.get_num_atoms() : 0; } catch (e) { return 0; } finally { if (m) m.delete(); } };
  const canonOf = (R, smi) => { let m = null; try { m = R.get_mol(smi); return m ? m.get_smiles() : null; } catch (e) { return null; } finally { if (m) m.delete(); } };
  const buyCache = new Map();
  function buyable(smi) {
    if (buyCache.has(smi)) return buyCache.get(smi);
    let ok = false;
    try { const g = C.parseSmiles(smi); ok = !!g && (C.isBasicMaterial(g) || C.buyableMatches(g, 2).some(x => x.exact)); } catch (e) { ok = false; }
    buyCache.set(smi, ok);
    return ok;
  }
  async function disconnect(R, smiles, opts = {}) {
    const T = await templates(opts.templates || 15000);
    const target = canonOf(R, smiles);
    if (!target) return [];
    const th = heavyOf(R, target), tm = R.get_mol(target), found = new Map();
    try {
      for (let i = 0; i < T.smarts.length; i++) {
        if (i % 400 === 399) { if (opts.onProgress) opts.onProgress(i / T.smarts.length); await breathe(); }
        const rx = reversed(R, T, i); if (!rx) continue;
        let res = null, ml = null;
        try { ml = new R.MolList(); ml.append(tm); res = rx.run_reactants(ml, 30); } catch (e) { if (ml) ml.delete(); continue; }
        try {
          for (let k = 0; k < res.size(); k++) {
            const pl = res.get(k), smis = [];
            for (let b = 0; b < pl.size(); b++) { const pm = pl.at(b); try { smis.push(pm.get_smiles()); } catch (e) { /* unreadable */ } pm.delete(); }
            pl.delete();
            if (!smis.length || smis.length > 3) continue;
            const canon = smis.map(s => canonOf(R, s));
            if (canon.some(c => !c || c === target)) continue;
            const key = canon.slice().sort().join('.');
            let f = found.get(key);
            if (!f) { f = { key, precursors: canon.slice().sort(), recorded: 0, templates: [], best: i }; found.set(key, f); }
            if (!f.templates.includes(T.id[i])) { f.templates.push(T.id[i]); f.recorded += T.recorded[i]; if (T.recorded[i] > T.recorded[f.best]) f.best = i; }
          }
        } finally { res.delete(); ml.delete(); }
      }
    } finally { tm.delete(); }
    const out = [];
    for (const f of found.values()) {
      const pre = f.precursors.map(smiles => ({ smiles, buyable: buyable(smiles), heavy: heavyOf(R, smiles) }));
      const maxH = Math.max(...pre.map(p => p.heavy)), totH = pre.reduce((t, p) => t + p.heavy, 0);
      if (pre.some(p => p.heavy < 1)) continue;
      const bigger = Math.max(0, maxH - th), nBuy = pre.filter(p => p.buyable).length;
      const score = 0.6 * Math.log(1 + f.recorded) + 2.5 * nBuy / pre.length - (bigger > 0 ? 0.8 + 0.35 * Math.min(bigger, 8) : 0) - 0.04 * Math.max(0, totH - th - 3);
      out.push({ precursors: pre, recorded: f.recorded, templates: f.templates, template: T.id[f.best], score, targetHeavy: th, bigger });
    }
    out.sort((a, b) => b.score - a.score);
    return out.slice(0, opts.top || 12);
  }
  /* what was recorded for this disconnection's kind of reaction, and what the reagents would also do to the precursors */
  async function annotate(R, d, opts = {}) {
    const row = await recipesOf(d.template);
    d.kind = row ? { n: row[0], yield: row[1], recipes: row[2], solvents: row[3], consumed: row[4], formed: row[5] } : null;
    const RR = window.REACTION_RULES || {}, F = RR.fragments || {};
    const names = ids => String(ids || '').split('+').map(x => x.trim()).filter(Boolean).map(id => shortName((F[id] || {}).name || id).split(' / ')[0]).filter(n => !/C–H/.test(n));
    d.title = d.kind ? `${names(d.kind.consumed).join(' + ') || 'a bond-forming reaction'} → ${names(d.kind.formed).join(' + ') || 'a new bond'}` : 'a recorded reaction type';
    d.options = [];
    if (d.kind) {
      const targets = new Set(String(d.kind.consumed).split('+').map(x => x.trim()).filter(Boolean));
      const pre = d.precursors.map(p => p.smiles);
      try {                                                   // and the groups that really disappear between these precursors and the target
        const P = window.RxnPage;
        if (P && P.fragmentCounts) { const a = P.fragmentCounts(R, pre), b = P.fragmentCounts(R, [opts.target || '']); Object.keys(a).forEach(id => { if ((a[id] || 0) > (b[id] || 0)) targets.add(id); }); }
      } catch (e) { /* the usual ones only */ }
      // a record with no reagent line usually means the reagents were not captured: list real reagent lines first
      const recipes = d.kind.recipes.slice().sort((a, b) => (a[0] === '' ? 1 : 0) - (b[0] === '' ? 1 : 0) || b[1] - a[1]);
      for (const rc of recipes.slice(0, opts.recipes || 3)) {
        const opt = { reagent: rc[0], count: rc[1], source: rc[2], yield: rc[3], collateral: null };
        const sol = (d.kind.solvents || []).find(x => x[0]);
        opt.cond = sol ? sol[0] : '';
        if (window.Atlas && opts.atlas !== false) {
          try {
            await window.Atlas.ready();
            const J = await window.Atlas.judge(R, pre, rc[0]);
            const byAtom = new Map();
            J.forEach(j => { const k = j.smiles + '|' + j.atom; if (!byAtom.has(k)) byAtom.set(k, { names: [shortName(j.name)], fids: [j.fid], p: j.pb, n: j.nReagentB || 0, confidence: j.confidenceB }); else { const c = byAtom.get(k); c.n = Math.max(c.n, j.nReagentB || 0); c.fids.push(j.fid); if (j.pb > c.p) { c.p = j.pb; c.confidence = j.confidenceB; } if (!c.names.includes(shortName(j.name))) c.names.push(shortName(j.name)); } });
            const by = [...byAtom.values()].filter(x => !x.fids.some(id => targets.has(id)) && !x.fids.every(id => id === 'F009' || id === 'F108'));
            let keep = 1; by.forEach(x => { if (x.p >= 0.03) keep *= 1 - x.p; });
            const unknown = by.filter(x => x.p < 0.10 && x.n < 15);
            opt.collateral = { unknown: unknown.slice(0, 2).map(x => x.names[0].split(' / ')[0].slice(0, 40)), spared: by.length - unknown.length, risk: 1 - keep, risky: by.filter(x => x.p >= 0.10).sort((a, b) => b.p - a.p).slice(0, 2).map(x => ({ names: [x.names[0].split(' / ')[0].slice(0, 40)], p: x.p })) };
          } catch (e) { opt.collateral = null; }
        }
        d.options.push(opt);
      }
    }
    return d;
  }
  /* a reagent line for people: "sodium hydroxide", "HATU + DIPEA"; the SMILES stays in the tooltip */
  function reagentLabel(R, line) {
    if (!line) return 'no reagent recorded';
    let name = '';
    try { name = window.RxnPage ? window.RxnPage.reagentLabel(R, line) : ''; } catch (e) { name = ''; }
    name = name || line;
    return name.length > 70 ? name.slice(0, 68) + '…' : name;
  }

  // ---------------------------------------------------------------- a whole route
  function transformOf(R, d) {
    const k = d.kind, top = d.options && d.options[0];
    const bits = [];
    if (top) bits.push(top.reagent ? reagentLabel(R, top.reagent) : 'no reagent recorded');
    if (top && top.cond) bits.push(top.cond);
    return {
      name: d.title + ' (recorded)',
      reagent: bits.join(', ') + (top ? `  (recorded ${fmtN(top.count)} times)` : ''),
      note: `This kind of reaction (${d.template}) was recorded ${fmtN(d.recorded)} times in patents${k && k.yield != null ? `, median yield ${k.yield} %` : ''}. A candidate from the reversed template: check it on these molecules in the Reactions tab.`,
      code: null,
    };
  }
  async function route(R, smiles, opts = {}) {
    const depth = opts.depth == null ? 4 : opts.depth, seen = opts.seen || new Set();
    const g = C.parseSmiles(smiles);
    if (buyable(smiles)) return { g, basic: true };
    if (depth <= 0 || seen.has(smiles)) return { g };
    seen.add(smiles);
    const options = opts.first ? [opts.first] : (await disconnect(R, smiles, { templates: opts.templates || 8000, top: 6 })).filter(d => d.bigger <= 2);
    for (const d of options.slice(0, 3)) {
      await annotate(R, d, { atlas: false, recipes: 1, target: smiles });
      const kids = [];
      let dead = 0;
      for (const p of d.precursors) {
        const child = await route(R, p.smiles, { depth: depth - 1, seen: new Set(seen), templates: opts.templates });
        kids.push(child);
        if (!child.basic && !(child.children && child.children.length)) dead++;
      }
      if (dead === 0 || d === options[Math.min(2, options.length - 1)]) {
        const t = transformOf(R, d);
        t.code = `${d.precursors.map(p => p.smiles).join('.')}>${d.options[0] ? d.options[0].reagent : ''}>${smiles}` + (d.options[0] && d.options[0].cond ? ` | cond=${d.options[0].cond}` : '');
        return { g, transform: t, children: kids, recorded: d.recorded };
      }
    }
    return { g };
  }

  // ---------------------------------------------------------------- the page section
  const STYLE = `
    .rl-sec { margin:14px 0; } .rl-sec > summary { cursor:pointer; font-weight:bold; font-size:15px; }
    .rl-card { border:1px solid #e3e3e3; border-radius:8px; padding:8px 10px; margin:8px 0; background:#fff; }
    .rl-eq { display:flex; align-items:center; gap:8px; flex-wrap:wrap; } .rl-eq .rl-plus { color:#888; } .rl-eq .rl-arrow { color:#555; font-size:20px; }
    .rl-mol { text-align:center; font-size:11px; color:#666; max-width:190px; } .rl-mol svg { max-width:100%; height:auto; }
    .rl-tag { display:inline-block; font-size:11px; border-radius:3px; padding:0 5px; margin-left:4px; background:#eef3fa; color:#1f4e78; } .rl-tag.buy { background:#e8f5e9; color:#1e7d32; } .rl-tag.make { background:#fff4e0; color:#a15c00; }
    .rl-title { font-weight:bold; margin-top:4px; } .rl-sub { color:#666; font-size:12px; margin:2px 0; }
    .rl-opt { margin:3px 0; font-size:13px; } .rl-warn { color:#b3261e; font-size:12px; } .rl-mid { color:#a15c00; font-size:12px; } .rl-ok { color:#1e7d32; font-size:12px; }
    .rl-btns { margin-top:6px; display:flex; gap:8px; flex-wrap:wrap; } .rl-btns button, .rl-open { font:12px Arial,sans-serif; padding:3px 10px; border:1px solid #1f4e78; background:#fff; color:#1f4e78; border-radius:4px; cursor:pointer; }
    .rl-btns button:hover, .rl-open:hover { background:#1f4e78; color:#fff; }
    .rl-status { color:#777; font-size:13px; margin:6px 0; } .rl-progress { height:4px; background:#eee; border-radius:2px; overflow:hidden; margin:4px 0; } .rl-progress div { height:100%; width:0; background:#1f4e78; transition:width .15s; }`;
  function drawn(smiles, w, h) { try { return window.MolDraw.svg(C.parseSmiles(smiles), w, h); } catch (e) { return el('span'); } }
  const nameOf = smiles => { try { const g = C.parseSmiles(smiles); const n = C.displayNameFor(g); return n || C.formula(g); } catch (e) { return ''; } };

  function card(R, d, host, target, startRoute) {
    const c = el('div', 'rl-card');
    const eq = el('div', 'rl-eq');
    d.precursors.forEach((p, k) => {
      if (k) eq.appendChild(el('span', 'rl-plus', '+'));
      const m = el('div', 'rl-mol'); m.appendChild(drawn(p.smiles, 150, 90));
      m.appendChild(el('div', null, `${esc(nameOf(p.smiles))}<span class="rl-tag ${p.buyable ? 'buy' : 'make'}">${p.buyable ? 'buy' : 'make'}</span>`));
      eq.appendChild(m);
    });
    eq.appendChild(el('span', 'rl-arrow', '⟶'));
    const t = el('div', 'rl-mol'); t.appendChild(drawn(target, 150, 90)); t.appendChild(el('div', null, 'target')); eq.appendChild(t);
    c.appendChild(eq);
    c.appendChild(el('div', 'rl-title', esc(d.title)));
    c.appendChild(el('div', 'rl-sub', `recorded ${fmtN(d.recorded)} times${d.templates.length > 1 ? ` (${d.templates.length} templates)` : ''}${d.kind && d.kind.yield != null ? ` · median yield ${d.kind.yield} %` : ''}` +
      `${d.kind && d.kind.solvents && d.kind.solvents.length ? ' · solvents: ' + esc(d.kind.solvents.slice(0, 3).map(s => (s[0] || 'none stated') + ' ' + Math.round(100 * s[1] / d.kind.n) + ' %').join(', ')) : ''}`));
    d.options.forEach(o => {
      const col = o.collateral;
      const risk = col ? (col.risky.length ? col.risky.map(x => `<span class="${x.p >= 0.3 ? 'rl-warn' : 'rl-mid'}">${x.p >= 0.3 ? '⚠ likely also reacts' : 'may also react'}: ${esc(x.names[0])} (${Math.round(x.p * 100)} %)</span>`).join(' ') : (col.spared > 0 && !(col.unknown && col.unknown.length) ? '<span class="rl-ok">recorded as sparing the other groups</span>' : '') + (col.unknown && col.unknown.length ? ` <span class="rl-mid">no precedent for ${esc(col.unknown.join(', '))} with this reagent</span>` : '')) : '';
      c.appendChild(el('div', 'rl-opt', `<b title="${esc(o.reagent)}">${esc(reagentLabel(R, o.reagent))}</b> <span class="rl-sub">${fmtN(o.count)} recorded${o.yield != null ? ' · best yield ' + o.yield + ' %' : ''}${o.cond ? ' · in ' + esc(o.cond) : ''}</span> ${risk}`));
    });
    const b = el('div', 'rl-btns');
    const open = el('button', null, 'Open this step in Reactions');
    open.title = 'Fills the Reactions tab with the precursors, the most recorded reagent and solvent: check what it does to these molecules';
    open.addEventListener('click', () => {
      const top = d.options[0];
      const code = `${d.precursors.map(p => p.smiles).join('.')}>${top ? top.reagent : ''}>` + (top && top.cond ? ` | cond=${top.cond}` : '');
      if (window.Reactions) window.Reactions.load(code);
    });
    b.appendChild(open);
    const plan = el('button', null, 'Plan the rest from here');
    plan.title = 'Take this step first, then work backwards from the precursors you cannot buy';
    plan.addEventListener('click', () => startRoute(d));
    b.appendChild(plan);
    c.appendChild(b);
    return c;
  }

  async function mount(host, smiles) {
    if (!document.getElementById('rl-style')) { const st = document.createElement('style'); st.id = 'rl-style'; st.textContent = STYLE; document.head.appendChild(st); }
    const sec = el('details', 'rl-sec'); sec.open = true;
    sec.appendChild(el('summary', null, 'Disconnections from recorded reactions <small style="font-weight:400;color:#777">1.09 million patent reactions, run backwards</small>'));
    const status = el('div', 'rl-status', 'Starting…'), bar = el('div', 'rl-progress'), fill = el('div'); bar.appendChild(fill);
    const list = el('div'), routeHost = el('div');
    sec.append(status, bar, list, routeHost);
    host.appendChild(sec);
    let R = null;
    try { R = window.RDKitLoad ? await window.RDKitLoad() : null; } catch (e) { R = null; }
    if (!R) { status.textContent = 'RDKit could not start, so the recorded reactions cannot be run.'; bar.remove(); return; }
    const target = canonOf(R, smiles);
    if (!target) { status.textContent = 'This structure could not be read by RDKit.'; bar.remove(); return; }
    status.textContent = 'Running the 15,000 commonest recorded reaction types backwards on the target…';
    let ds = [];
    try { ds = await disconnect(R, target, { templates: 15000, top: 10, onProgress: f => { fill.style.width = Math.round(f * 100) + '%'; } }); } catch (e) { status.textContent = 'Something went wrong: ' + e.message; bar.remove(); return; }
    bar.remove();
    if (!ds.length) { status.textContent = 'No recorded reaction type fits this structure as a product (it may be very small, or an unusual skeleton).'; return; }
    status.textContent = `${ds.length} disconnections, best first (ranked by how often each kind was recorded, how much simpler the precursors are, and what can be bought). Looking up recorded reagents…`;
    for (const d of ds) await annotate(R, d, { atlas: true, recipes: 3, target });
    status.textContent = `${ds.length} disconnections, best first. A proposal is a candidate: open it in the Reactions tab to see what it does to these molecules.`;
    const startRoute = async first => {
      routeHost.innerHTML = '<div class="rl-status">Planning a route from recorded reactions (works backwards from every piece that cannot be bought; a few seconds per step)…</div>';
      let tree = null;
      try { tree = await route(R, target, { depth: 4, first }); } catch (e) { routeHost.innerHTML = `<div class="banner">The planner failed: ${esc(e.message)}</div>`; return; }
      routeHost.innerHTML = '';
      const inner = el('div');
      routeHost.appendChild(inner);
      if (!window.Retro) return;
      window.Retro.renderPlan(inner, C.parseSmiles(target), tree, {
        heading: 'A route built from recorded reactions', deadTag: 'no recorded route', recordedNote: false, depthNote: false,
        stepBadge: s => (s.t && s.t.code ? `<button type="button" class="rl-open" data-code="${esc(s.t.code)}" title="Fill the Reactions tab with this step">Open in Reactions</button> ` : ''),
      });
      inner.addEventListener('click', e => { const b = e.target.closest('.rl-open'); if (b && window.Reactions) window.Reactions.load(b.dataset.code); });
    };
    ds.forEach(d => list.appendChild(card(R, d, list, target, startRoute)));
    const auto = el('div', 'rl-btns'); const ab = el('button', null, 'Plan a whole route from recorded reactions');
    ab.addEventListener('click', () => startRoute(null)); auto.appendChild(ab); list.insertBefore(auto, list.firstChild);
  }

  return { disconnect, annotate, route, mount, templates, part };
})();
