/* Advisor — which reagents and conditions could you use for these starting materials?

   Two sources, shown separately because they are different kinds of evidence:

   Textbook rules  The 102 reactions of reaction_rules.xlsx whose reactant pattern fits the molecules: the reagent each uses, what it makes, and the
                   conditions it requires. These come with explanations, mechanisms and stereochemistry once run on the page.
   Recorded        The relevance model (js/relevance.js) says which kinds of recorded reaction these molecules take part in, with no reagent given;
                   data/learned/c/ (tools/cond_build.py) says what chemists used for each kind: the reagent lines, solvents and median yield of the
                   1.09 million recorded reactions. Patent records carry no temperature or time, so none is invented here.

   Every reagent idea is screened with the Fragment Atlas (js/atlas.js): for the molecules at hand, how likely is that reagent to transform a group the
   reaction is NOT meant for? A reagent that does the job and spares the rest ranks above one that does the job and also reduces the nitro group.

   Advisor.toggle() / open() / close()          the panel on the Reactions page (button "Suggest reagents")
   Advisor.suggest(state, onStatus)             -> { textbook: [...], recorded: [...], notes: [...] }   (the data behind the panel; used by the tests)   */
window.Advisor = (() => {
  'use strict';
  const RP = () => window.RxnPage;
  const esc = s => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const shortName = n => String(n || '').split(' (')[0];
  const pct = p => (p >= 0.995 ? '>99' : p < 0.005 ? '<1' : String(Math.round(p * 100))) + ' %';
  const arene = new Set(['F009', 'F108']);

  /* ---- the groups a reaction is for, and what a reagent would also do (Fragment Atlas) ---- */
  function namesOf(ids) {
    const RR = window.REACTION_RULES || {}, F = RR.fragments || {};
    const list = String(ids || '').split('+').map(x => x.trim()).filter(Boolean).map(id => ({ id, name: shortName((F[id] || {}).name || id) }));
    const keep = list.filter(x => !/C–H/.test(x.name));
    return (keep.length ? keep : list).map(x => x.name);
  }
  function consumedGroups(R, subs, products) {
    const RP_ = RP(), out = new Set(), RR = window.REACTION_RULES || {};
    try {
      const a = RP_.fragmentCounts(R, subs), b = RP_.fragmentCounts(R, products || []);
      Object.keys(a).forEach(id => { if ((a[id] || 0) > (b[id] || 0) && !/C–H/.test(((RR.fragments || {})[id] || {}).name || '')) out.add(id); });
    } catch (e) { /* none */ }
    return out;
  }
  /* what the reagent line would also do: the Atlas's chance for every group that the reaction is not for */
  async function collateral(R, subs, reagentLine, targets) {
    const AT = window.Atlas;
    if (!AT) return null;
    let J = []; try { J = await AT.judge(R, subs, reagentLine); } catch (e) { return null; }
    if (!J.length) return null;
    const byAtom = new Map();
    J.forEach(j => {
      const k = j.smiles + '|' + j.atom;
      if (!byAtom.has(k)) byAtom.set(k, { names: [shortName(j.name)], fids: [j.fid], p: j.p, evidence: j.evidence, confidence: j.confidence });
      else { const c = byAtom.get(k); if (!c.names.includes(shortName(j.name))) c.names.push(shortName(j.name)); c.fids.push(j.fid); if (j.p > c.p) { c.p = j.p; c.evidence = j.evidence; c.confidence = j.confidence; } }
    });
    const bystanders = [...byAtom.values()].filter(x => !x.fids.some(id => targets.has(id)) && !x.fids.every(id => arene.has(id)));
    const risky = bystanders.filter(x => x.p >= 0.15).sort((a, b) => b.p - a.p);
    let keep = 1; bystanders.forEach(x => { if (x.p >= 0.05) keep *= (1 - x.p); });
    // how likely the reagent is to transform the group the reaction is for (the best such group), judged the same way
    const tg = [...byAtom.values()].filter(x => x.fids.some(id => targets.has(id)));
    return { risk: 1 - keep, risky: risky.slice(0, 3), spared: bystanders.length - risky.length, targetP: tg.length ? Math.max(...tg.map(x => x.p)) : null };
  }

  /* ---- gathering ---- */
  async function textbook(R, state, status) {
    const P = RP();
    status('Checking the textbook rules…');
    const ideas = P.textbookIdeas(R, state.subs), out = [];
    for (const it of ideas) {
      const targets = new Set(); (it.rule.consumes || []).forEach(alt => alt.forEach(id => targets.add(id)));
      const options = [];
      for (const rg of it.reagents.slice(0, 6)) {
        const col = await collateral(R, state.subs, rg.smiles, targets);
        options.push({ label: rg.name, reagent: rg.smiles, collateral: col });
      }
      if (!it.reagents.length) options.push({ label: 'no reagent', reagent: '', collateral: null });
      const val = o => (o.collateral ? (o.collateral.targetP == null ? 0.12 : o.collateral.targetP) * (1 - 0.7 * o.collateral.risk) : 0.12);
      options.sort((a, b) => val(b) - val(a));
      out.push({ kind: 'textbook', rule: it.rule, products: it.products, alternatives: it.alternatives, cond: it.conditions, options,
                 plausible: options.length ? val(options[0]) : 0.12 });
    }
    return out;
  }
  async function recorded(R, state, status) {
    const P = RP(), notes = [];
    status('Loading the reaction model (about 30 MB, once)…');
    let rel = null; try { rel = await P.relevance(R, state.subs, ''); } catch (e) { rel = null; }
    if (!rel || !rel.candidates || !rel.candidates.length) {
      const heavy = state.subs.reduce((t, s) => t + ((s.match(/[A-Z]/g) || []).length), 0);
      return { items: [], notes: [heavy < 8 ? 'The recorded reactions are patent-sized molecules (20 or more heavy atoms is the norm); these are too small for the learned model, so only the textbook rules are offered.' : 'The learned model proposes nothing for these molecules.'] };
    }
    const items = [], seen = new Set();
    status('Looking up what was used for each kind of reaction…');
    for (const c of rel.candidates) {
      if (items.length >= 8 || c.share < 0.01) break;
      const best = c.rep && c.rep.best; if (!best) continue;
      const key = c.rep.products.slice().sort().join('.');
      if (seen.has(key)) continue; seen.add(key);
      const tab = await P.learnedPart('c/' + String(Math.floor(best.t / 2000)).padStart(3, '0'));
      const row = tab && tab[String(best.t)];
      const targets = consumedGroups(R, state.subs, c.rep.products);
      const options = [];
      for (const rc of (row ? row[2] : []).slice(0, 4)) {
        const line = rc[0], col = await collateral(R, state.subs, line, targets);
        const sol = (row[3] || []).find(x => x[0]) || null;
        options.push({ label: P.reagentLabel(R, line) || 'no reagent', reagent: line, count: rc[1], share: rc[1] / row[0], source: rc[2], yield: rc[3], collateral: col,
                       cond: sol ? sol[0] : '' });
      }
      options.sort((a, b) => (b.share * (1 - 0.7 * (b.collateral ? b.collateral.risk : 0))) - (a.share * (1 - 0.7 * (a.collateral ? a.collateral.risk : 0))));
      const top = options[0];
      items.push({ kind: 'recorded', key, first: items.length === 0 && c === rel.candidates[0], confidence: rel.confidence, products: c.rep.products, share: c.share, template: best.id, n: row ? row[0] : null, yield: row ? row[1] : null,
                   consumed: row ? namesOf(row[4]) : [], formed: row ? namesOf(row[5]) : [], solvents: row ? row[3] : [], options,
                   score: c.share * (top ? top.share * (1 - 0.7 * (top.collateral ? top.collateral.risk : 0)) : 0.01) });
    }
    items.sort((a, b) => b.score - a.score);
    return { items, notes, confidence: rel.confidence };
  }
  async function suggest(state, status) {
    const P = RP(), R = await P.loadRD();
    status = status || (() => {});
    if (!R) return { textbook: [], recorded: [], notes: ['RDKit could not start.'] };
    if (!state.subs.length) return { textbook: [], recorded: [], notes: ['Add the starting molecules on the left first.'] };
    const tb = await textbook(R, state, status);
    const rec = await recorded(R, state, status);
    // a recorded kind of reaction that gives the same product as a textbook rule is one idea: its recorded reagents join the rule's card
    const only = [];
    for (const it of rec.items) {
      // the same product(s), or the recorded product is one of the rule's products (hydrolysis gives the acid and the alcohol; the record keeps the acid)
      const t = tb.find(x => !x.recorded && (x.products.slice().sort().join('.') === it.key || it.products.every(p => x.products.includes(p))));
      if (t) { t.recorded = it; } else only.push(it);
    }
    tb.forEach(t => { t.score = t.recorded ? 1 + t.recorded.share : t.plausible; });
    tb.sort((a, b) => b.score - a.score);
    return { textbook: tb, recorded: only, notes: rec.notes };
  }

  // ---------------------------------------------------------------- the panel
  let panel = null, running = 0;
  const STYLE = `
    #rxn-advice { position:fixed; top:96px; right:12px; bottom:12px; width:560px; max-width:calc(100vw - 24px); box-sizing:border-box; background:#fff; border:1px solid #ccc; border-radius:8px;
                  box-shadow:0 6px 24px rgba(0,0,0,.16); z-index:9; display:flex; flex-direction:column; font:13px Arial,sans-serif; color:#333; }
    #rxn-advice[hidden] { display:none; }
    .adv-head { display:flex; align-items:center; gap:10px; padding:10px 12px 6px; } .adv-head b { flex:1; font-size:14px; }
    .adv-close { border:0; background:none; font-size:20px; line-height:1; cursor:pointer; color:#777; }
    .adv-status { padding:0 12px 6px; color:#777; min-height:18px; }
    .adv-list { flex:1; overflow:auto; padding:0 12px 12px; }
    .adv-sec { margin:12px 0 4px; font-size:11px; color:#888; text-transform:uppercase; letter-spacing:.04em; }
    .adv-note { color:#777; font-size:12px; margin:6px 0; }
    .adv-card { border:1px solid #e3e3e3; border-radius:8px; padding:8px 10px; margin:8px 0; }
    .adv-title { font-weight:bold; } .adv-tag { font-weight:normal; font-size:11px; background:#eef3fa; color:#1f4e78; border-radius:3px; padding:0 5px; margin-left:6px; }
    .adv-body { display:flex; gap:10px; align-items:center; margin-top:4px; } .adv-prod { display:flex; gap:2px; align-items:center; flex-shrink:0; max-width:46%; overflow:hidden; }
    .adv-sub { color:#777; font-size:12px; margin:2px 0; }
    .adv-opt { display:flex; align-items:baseline; gap:6px; flex-wrap:wrap; margin:4px 0; padding-top:4px; border-top:1px dotted #eee; }
    .adv-opt:first-of-type { border-top:0; }
    .adv-name { font-weight:bold; } .adv-n { color:#777; font-size:12px; }
    .adv-use { font:12px Arial,sans-serif; padding:2px 10px; border:1px solid #1f4e78; background:#fff; color:#1f4e78; border-radius:4px; cursor:pointer; margin-left:auto; }
    .adv-use:hover { background:#1f4e78; color:#fff; }
    .adv-warn { color:#b3261e; font-size:12px; width:100%; } .adv-mid { color:#a15c00; font-size:12px; width:100%; } .adv-ok { color:#1e7d32; font-size:12px; width:100%; }`;
  function ensure() {
    if (panel) return panel;
    const host = document.getElementById('page-reactions');
    if (!document.getElementById('adv-style')) { const st = document.createElement('style'); st.id = 'adv-style'; st.textContent = STYLE; document.head.appendChild(st); }
    panel = document.createElement('div'); panel.id = 'rxn-advice'; panel.hidden = true;
    panel.innerHTML = '<div class="adv-head"><b>Reagents and conditions for these molecules</b><button type="button" class="adv-refresh" title="Look again at the molecules on the page">Refresh</button><button type="button" class="adv-close" title="Close">×</button></div>' +
      '<div class="adv-status"></div><div class="adv-list"></div>';
    host.appendChild(panel);
    panel.querySelector('.adv-close').addEventListener('click', close);
    panel.querySelector('.adv-refresh').addEventListener('click', run);
    return panel;
  }
  const btn = () => document.getElementById('rxn-suggest');
  function open() { ensure().hidden = false; btn().classList.add('on'); run(); }
  function close() { if (panel) panel.hidden = true; btn().classList.remove('on'); running++; }
  function toggle() { if (panel && !panel.hidden) close(); else open(); }

  const optHtml = (o, i, kind) => {
    const col = o.collateral;
    let risk = '';
    if (col && col.risky.length) risk = col.risky.map((x, k) => `<div class="${x.p >= 0.5 ? 'adv-warn' : 'adv-mid'}">${x.p >= 0.5 ? '⚠ likely also reacts' : 'may also react'}: ${esc(x.names.join(' / '))} (${pct(x.p)}, ${esc(x.confidence)} confidence)</div>`).join('');
    else if (col) risk = `<div class="adv-ok">spares the other groups (each under 15 %)</div>`;
    const bits = [];
    if (o.count) bits.push(`${o.count.toLocaleString()} recorded${o.share ? ` (${pct(o.share)} of this kind)` : ''}`);
    if (o.yield != null) bits.push(`best yield ${o.yield} %${o.source ? ' · ' + esc(o.source) : ''}`);
    return `<div class="adv-opt"><span class="adv-name">${esc(o.label)}</span><span class="adv-n">${bits.join(' · ')}${o.cond ? (bits.length ? ' · ' : '') + (o.condLabel ? esc(o.condLabel) + ': ' : 'in ') + esc(o.cond) : ''}</span>` +
      `<button type="button" class="adv-use" data-k="${kind}" data-i="${i}" title="Put this reagent${o.cond ? ' and solvent' : ''} on the page and run it">Use ▸</button>${risk}</div>`;
  };
  function render(res, state) {
    const list = panel.querySelector('.adv-list');
    const thumbs = (products) => products.slice(0, 3).map(s => { const t = RP().thumb(s, 72); return t ? t.outerHTML : ''; }).join('<span class="adv-n">+</span>');
    let html = '';
    const likely = res.textbook.filter(it => it.recorded || it.plausible >= 0.15), unlikely = res.textbook.filter(it => !likely.includes(it));
    const card = (it, i) => {
        const rc = it.recorded, tag = rc && rc.first ? `<span class="adv-tag">the model's first choice${rc.confidence != null ? ' · about ' + pct(rc.confidence) + ' sure' : ''}</span>` : '';
        const have = new Set(it.options.map(o => o.reagent));
        const extra = rc ? rc.options.filter(o => !have.has(o.reagent)) : [];
        return `<div class="adv-card"><div class="adv-title">${esc(it.rule.name)}<span class="adv-tag">textbook ${esc(it.rule.id)}</span>${tag}</div>` +
          `<div class="adv-body"><div class="adv-prod">${thumbs(it.products)}</div><div style="flex:1;min-width:0">` +
          `<div class="adv-sub">Textbook conditions: ${esc(it.rule.conditions || '')}${it.alternatives > 1 ? ` · ${it.alternatives} places it could react: the page picks by reactivity` : ''}` +
          `${rc && rc.n ? ` · also recorded ${rc.n.toLocaleString()} times${rc.yield != null ? ', median yield ' + rc.yield + ' %' : ''}` : ''}</div>` +
          it.options.map((o, k) => optHtml(Object.assign({}, o, { cond: it.cond, condLabel: 'conditions' }), `${i}:${k}`, 't')).join('') +
          (extra.length ? '<div class="adv-sub" style="margin-top:6px">Also used in patents for this product:</div>' + extra.map((o, k) => optHtml(o, `${i}:${it.options.length + k}`, 't')).join('') : '') +
          '</div></div></div>';
    };
    if (likely.length) html += '<div class="adv-sec">Textbook reactions (explained, with mechanisms when run)</div>' + likely.map(it => card(it, res.textbook.indexOf(it))).join('');
    if (res.recorded.length) {
      html += '<div class="adv-sec">Recorded in patents (what chemists used for this kind of reaction)</div>';
      html += res.recorded.map((it, i) => {
        const title = it.consumed.length || it.formed.length ? `${esc(it.consumed.join(' + ') || 'a reaction')} → ${esc(it.formed.join(' + ') || 'a new bond')}` : 'a recorded reaction type';
        const first = it.first ? `<span class="adv-tag">the model's first choice${it.confidence != null ? ' · about ' + pct(it.confidence) + ' sure' : ''}</span>` : `<span class="adv-tag">${pct(it.share)} of the model's weight</span>`;
        return `<div class="adv-card"><div class="adv-title">${title}${first}</div>` +
          `<div class="adv-body"><div class="adv-prod">${thumbs(it.products)}</div><div style="flex:1;min-width:0">` +
          `<div class="adv-sub">${it.n ? `${it.n.toLocaleString()} recorded reactions of this kind` : ''}${it.yield != null ? ` · median yield ${it.yield} %` : ''}${it.solvents.length ? ' · solvents: ' + esc(it.solvents.slice(0, 3).map(s => (s[0] || 'none stated') + ' ' + Math.round(100 * s[1] / it.n) + ' %').join(', ')) : ''}</div>` +
          it.options.map((o, k) => optHtml(o, `${i}:${k}`, 'r')).join('') + '</div></div></div>';
      }).join('');
    }
    if (unlikely.length) html += `<details><summary class="adv-note" style="cursor:pointer">${unlikely.length} textbook reaction${unlikely.length > 1 ? 's' : ''} fit the pattern but are unlikely (recorded reactions rarely transform that group with these reagents)</summary>` + unlikely.map(it => card(it, res.textbook.indexOf(it))).join('') + '</details>';
    if (!html) html = '<div class="adv-note">Nothing to suggest for these molecules.</div>';
    if (res.notes.length) html += res.notes.map(n => `<div class="adv-note">${esc(n)}</div>`).join('');
    html += '<div class="adv-note">Patent records have no temperature or time: take those from the textbook rule, the literature or your own judgement. The percentages for other groups are rates over the conditions people used (Fragment Atlas).</div>';
    list.innerHTML = html;
    list.querySelectorAll('.adv-use').forEach(b => b.addEventListener('click', async () => {
      const [i, k] = b.dataset.i.split(':').map(Number);
      const it = (b.dataset.k === 't' ? res.textbook : res.recorded)[i];
      const all = b.dataset.k === 't' ? it.options.concat((it.recorded ? it.recorded.options : []).filter(x => !it.options.some(o => o.reagent === x.reagent))) : it.options;
      const o = all[k], fromRule = b.dataset.k === 't' && k < it.options.length;
      await RP().apply({ reagent: o.reagent, cond: fromRule ? it.cond : o.cond });
    }));
  }
  async function run() {
    const my = ++running, st = panel.querySelector('.adv-status'), list = panel.querySelector('.adv-list');
    list.innerHTML = '';
    const state = RP().state();
    if (!state.subs.length) { st.textContent = ''; list.innerHTML = '<div class="adv-note">Add the starting molecules on the left (the boxes, or type a reaction code), then press Refresh.</div>'; return; }
    let res;
    try { res = await suggest(state, t => { if (my === running) st.textContent = t; }); } catch (e) { if (my === running) st.textContent = 'Something went wrong: ' + e.message; return; }
    if (my !== running) return;                                       // the panel was closed or refreshed meanwhile
    st.textContent = `${res.textbook.length} textbook and ${res.recorded.length} recorded ideas for ${state.subs.length} molecule${state.subs.length > 1 ? 's' : ''}`;
    render(res, state);
  }

  return { toggle, open, close, suggest, run };
})();
