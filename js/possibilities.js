/* Possibilities page — "here is what I have; what can I make?"

   Two answers, both from the lifted engines:
     1. Direct reactions between the chemicals on the shelf: every molecule × every reagent (a written
        reagent like NaBH4, or another molecule on the shelf) through RXN.run.
     2. Common target molecules and how far the shelf gets you: each target's starting materials come
        from Chem.retroAnalyze (computed once in the background and cached in localStorage); the ones
        not on the shelf are what you would need to buy — grouped by 0, 1, 2, 3, 4+ extra purchases.
   Targets = the common-name table + the library's molecules that have a trivial name (aspirin,
   caffeine …), minus anything that is itself a starting material.                                   */
(() => {
  'use strict';
  const C = window.Chem, M = window.MolDraw, R = window.Reagents;
  const addBox = document.getElementById('possAdd'), chips = document.getElementById('possChips'), out = document.getElementById('possOut');
  const status = document.getElementById('possStatus');
  if (!addBox || !C || typeof RXN === 'undefined') return;

  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const el = (tag, cls, html) => { const n = document.createElement(tag); if (cls) n.className = cls; if (html != null) n.innerHTML = html; return n; };
  const nameOf = g => { try { return C.displayNameFor(g) || C.formula(g); } catch (e) { return '?'; } };
  const keyOf = g => { try { return C.canonicalKey(g); } catch (e) { return null; } };
  const smilesOf = g => { try { return C.toSmiles(g); } catch (e) { return ''; } };
  const store = { get(k) { try { return JSON.parse(localStorage.getItem(k)); } catch (e) { return null; } }, set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) {} } };

  /* ================= the shelf ================= */
  let shelf = store.get('poss.shelf') || [];          // [{text}]
  const parsed = new Map();                             // text -> {kind:'mol'|'reagent'|'bad', g, key, name, ruleId, label}
  function parseItem(text) {
    if (parsed.has(text)) return parsed.get(text);
    let item = null;
    try { const g = C.searchMolecule(text).graph; item = { kind: 'mol', g, key: keyOf(g), name: nameOf(g), smiles: smilesOf(g) }; } catch (e) {}
    if (!item) {
      window.REAGENT_SUBSTRATE = null;
      let r = null; try { r = R.resolve(text); } catch (e) {}
      if (r && r.id) { const rule = C.RULES.find(x => x.id === r.id); item = { kind: 'reagent', ruleId: r.id, name: rule ? rule.label : text }; }
    }
    if (!item) item = { kind: 'bad', name: text };
    parsed.set(text, item);
    return item;
  }
  function renderShelf() {
    chips.innerHTML = '';
    shelf.forEach((s, i) => {
      const it = parseItem(s.text);
      const chip = el('span', 'chip ' + it.kind);
      chip.innerHTML = esc(it.name) + (it.kind === 'mol' && it.name.toLowerCase() !== s.text.toLowerCase() ? ` <small>${esc(s.text)}</small>` : '') +
        (it.kind === 'reagent' ? ' <small>reagent</small>' : it.kind === 'bad' ? ' <small>not recognised</small>' : '');
      const x = el('button', null, '✕'); x.title = 'Remove';
      x.addEventListener('click', () => { shelf.splice(i, 1); save(); });
      chip.appendChild(x); chips.appendChild(chip);
    });
    if (!shelf.length) chips.innerHTML = '<span class="hint">Nothing on the shelf yet.</span>';
  }
  function add(text) {
    text = (text || '').trim(); if (!text) return;
    if (shelf.some(s => s.text.toLowerCase() === text.toLowerCase())) return;
    shelf.push({ text }); addBox.value = ''; save();
  }
  function save() { store.set('poss.shelf', shelf); renderShelf(); refresh(); }

  /* ================= targets: computed once, cached ================= */
  // the planner's rule list is part of the key: adding a disconnection (js/retro_ext.js) recomputes the cached plans
  const CACHE_KEY = 'poss.targets.v1.t' + ((window.Chem && window.Chem.RETRO_TRANSFORMS) || []).length + ((window.RetroExt && window.RetroExt.version) || '');
  let targets = null;                                   // [{name, smiles, key, tags, leaves:[{key,smiles,name}], nSteps, dead, route:[str]}]
  let building = false;
  function targetList() {
    const seen = new Set(), list = [];
    const push = (name, smiles, tags) => {
      let g; try { g = C.parseSmiles(smiles); } catch (e) { return; }
      const key = keyOf(g); if (!key || seen.has(key)) return;
      seen.add(key); list.push({ name, smiles, key, tags: tags || [] });
    };
    Object.entries(C.COMMON_NAMES || {}).forEach(([n, s]) => push(n, s, null));
    (C.MOL_LIBRARY || []).filter(e => e.alt && e.alt.length).forEach(e => push(e.alt[0] || e.name, e.smiles, e.tags));
    /* fill in tags for the common-name entries from the library where it has them */
    const byKey = new Map(); (C.MOL_LIBRARY || []).forEach(e => { try { byKey.set(keyOf(C.parseSmiles(e.smiles)), e.tags); } catch (x) {} });
    list.forEach(t => { if (!t.tags.length && byKey.has(t.key)) t.tags = byKey.get(t.key) || []; });
    return list;
  }
  function analyse(t) {
    const g = C.parseSmiles(t.smiles);
    const tree = C.retroAnalyze(g);
    if (tree.basic) return null;                                 // a starting material itself: not a target
    const leaves = new Map(); let nSteps = 0, dead = false; const route = [];
    (function walk(n) {
      if (n.basic) { const k = keyOf(n.g); if (k && !leaves.has(k)) leaves.set(k, { key: k, smiles: smilesOf(n.g), name: nameOf(n.g) }); return; }
      if (!n.children || !n.children.length) { dead = true; const k = keyOf(n.g); if (k && !leaves.has(k)) leaves.set(k, { key: k, smiles: smilesOf(n.g), name: nameOf(n.g), dead: true }); return; }
      n.children.forEach(walk);
      nSteps++;
      route.push((n.transform && n.transform.name) || '?');
    })(tree);
    return { leaves: [...leaves.values()], nSteps, dead, route };
  }
  function buildTargets(done) {
    const cached = store.get(CACHE_KEY);
    const list = targetList();
    if (cached && cached.n === list.length && cached.items) {
      targets = list.map(t => Object.assign(t, cached.items[t.key])).filter(t => t.leaves);
      done(); return;
    }
    if (building) return; building = true;
    const items = {}; let i = 0;
    const bar = el('div', 'progress'); const fill = el('div'); bar.appendChild(fill);
    status.textContent = 'Analysing ' + list.length + ' common molecules (once; the result is kept)…';
    status.appendChild(bar);
    (function step() {
      const end = Math.min(list.length, i + 6);
      for (; i < end; i++) { try { const a = analyse(list[i]); if (a) items[list[i].key] = a; } catch (e) {} }
      fill.style.width = Math.round(100 * i / list.length) + '%';
      if (i < list.length) { setTimeout(step, 0); return; }
      store.set(CACHE_KEY, { n: list.length, items });
      targets = list.map(t => Object.assign(t, items[t.key])).filter(t => t.leaves);
      building = false; status.textContent = ''; done();
    })();
  }

  /* ================= results ================= */
  let have = new Set();                                 // canonical keys of the molecules on the shelf
  function refresh() {
    if (!shelf.length) { out.innerHTML = '<p class="hint">Add the chemicals you have to see what they can make.</p>'; return; }
    if (!targets) { buildTargets(refresh); if (!targets) { out.innerHTML = '<p class="hint">Working out the starting materials of every common molecule — this happens once.</p>'; return; } }
    out.innerHTML = '';
    const items = shelf.map(s => parseItem(s.text));
    const mols = items.filter(i => i.kind === 'mol'), reagents = items.filter(i => i.kind === 'reagent');
    have = new Set(mols.map(m => m.key));
    const maxMissing = +document.getElementById('possMax').value;
    const filter = document.getElementById('possFilter').value.trim().toLowerCase();

    /* 1. direct reactions between what is on the shelf */
    out.appendChild(el('h2', null, 'Reactions between your chemicals <small>one step, major product</small>'));
    const pairs = el('div', 'pairs'); const seenProd = new Set(); let nPairs = 0;
    for (const sub of mols) {
      const partners = reagents.map(r => ({ id: r.ruleId, label: r.name, g: null }))
        .concat(mols.filter(m => m !== sub).map(m => { const rule = R.moleculeRule(m.g, m.name); return rule ? { id: rule.id, label: m.name, g: m.g } : null; }).filter(Boolean));
      for (const p of partners) {
        let res = null; try { window.REAGENT_SUBSTRATE = sub.g; res = RXN.run(sub.g, p.id, {}); } catch (e) {}
        if (!res || res.error || !res.major) continue;
        const k = keyOf(res.major.graph); if (!k || seenProd.has(k + '|' + sub.key + '|' + p.id)) continue; seenProd.add(k + '|' + sub.key + '|' + p.id);
        if (filter && !(res.major.name + ' ' + sub.name + ' ' + p.label).toLowerCase().includes(filter)) continue;
        nPairs++;
        const card = el('div', 'pair');
        card.appendChild(M.svg(res.major.graph, 110, 70));
        const txt = el('div', 'txt', `<b>${esc(res.major.name)}</b><div class="sub">${esc(sub.name)} + ${esc(p.label)}` +
          (res.outcomes.length > 1 ? ` · ${Math.round(res.major.share * 100)}% major` : '') + (res.major.stereo ? ' · ' + esc(res.major.stereo.kind) : '') + '</div>');
        const go = el('button', 'plan', 'Open in Reactions'); go.addEventListener('click', () => { document.querySelector('.tab[data-page="reactions"]').click(); window.Reactions?.set(p.g ? [sub.smiles, smilesOf(p.g)] : [sub.smiles], p.g ? '' : p.label); });
        txt.appendChild(go); card.appendChild(txt); pairs.appendChild(card);
      }
    }
    out.appendChild(nPairs ? pairs : el('p', 'hint', mols.length < 1 ? 'Add at least one molecule.' : 'None of the taught reactions runs between these ' + (mols.length + reagents.length) + ' items.'));

    /* 2. common targets, by how many purchases they need. A target counts only if at least one of its
       starting materials is already on the shelf — otherwise every one-ingredient molecule would be
       "one purchase away" whatever you own. */
    const limit = maxMissing >= 4 ? 99 : maxMissing;
    const rows = targets.map(t => {
      const missing = t.leaves.filter(l => !have.has(l.key));
      return { t, missing, n: missing.length };
    }).filter(r => r.n <= limit && r.n < r.t.leaves.length)
      .filter(r => !filter || (r.t.name + ' ' + r.t.tags.join(' ') + ' ' + r.t.route.join(' ')).toLowerCase().includes(filter));
    const groups = new Map();
    rows.forEach(r => { const k = Math.min(r.n, 4); if (!groups.has(k)) groups.set(k, []); groups.get(k).push(r); });
    const heading = { 0: 'You can make these now', 1: 'One more purchase away', 2: 'Two more purchases away', 3: 'Three more purchases away', 4: 'Four or more purchases away' };
    for (let k = 0; k <= 4; k++) {
      const rs = groups.get(k); if (!rs) { if (k <= maxMissing) out.appendChild(el('h2', null, heading[k] + ' <small>none</small>')); continue; }
      rs.sort((a, b) => a.t.nSteps - b.t.nSteps || a.t.name.localeCompare(b.t.name));
      out.appendChild(el('h2', null, heading[k] + ` <small>${rs.length} molecule${rs.length === 1 ? '' : 's'}</small>`));
      if (k === 0) { out.appendChild(targetCards(rs)); continue; }
      /* group by the exact set of purchases, biggest unlock first */
      const bySet = new Map();
      rs.forEach(r => { const id = r.missing.map(m => m.key).sort().join('|'); if (!bySet.has(id)) bySet.set(id, { missing: r.missing, rows: [] }); bySet.get(id).rows.push(r); });
      [...bySet.values()].sort((a, b) => b.rows.length - a.rows.length).slice(0, 40).forEach(grp => {
        const d = el('details', 'unlock'); if (grp.rows.length >= 3 || bySet.size <= 3) d.open = true;
        d.appendChild(el('summary', null, 'Buy ' + grp.missing.map(m => `<span class="need"><b>${esc(m.name)}</b></span>`).join(' + ') +
          ` <small>→ unlocks ${grp.rows.length} molecule${grp.rows.length === 1 ? '' : 's'}</small>`));
        d.appendChild(el('div', 'links', grp.missing.map(m => `${esc(m.name)}: ` + shopLinks(m.name)).join('<br>')));
        d.appendChild(targetCards(grp.rows));
        out.appendChild(d);
      });
    }
    if (mols.length && !rows.length) out.appendChild(el('p', 'hint', 'Nothing within that many purchases — allow more extra molecules above.'));
  }
  const shopLinks = name => {
    const q = encodeURIComponent(name);
    return `<a href="https://pubchem.ncbi.nlm.nih.gov/#query=${q}" target="_blank" rel="noopener">PubChem</a>` +
      `<a href="https://www.sigmaaldrich.com/US/en/search/${q}?focus=products" target="_blank" rel="noopener">Sigma-Aldrich</a>` +
      `<a href="https://www.fishersci.com/us/en/catalog/search/products?keyword=${q}" target="_blank" rel="noopener">Fisher</a>`;
  };
  function targetCards(rs) {
    const box = el('div', 'targets');
    rs.forEach(r => {
      const t = r.t;
      const card = el('div', 'target');
      let g = null; try { g = C.parseSmiles(t.smiles); } catch (e) {}
      if (g) card.appendChild(M.svg(g, 230, 90));
      card.appendChild(el('div', 'nm', esc(t.name)));
      card.appendChild(el('div', 'sub', esc(t.nSteps + ' step' + (t.nSteps === 1 ? '' : 's') + ' · ' + t.leaves.length + ' starting material' + (t.leaves.length === 1 ? '' : 's')) + (t.dead ? ' · one piece has no taught route' : '')));
      card.appendChild(el('div', 'need', 'From: ' + t.leaves.map(l => have.has(l.key) ? esc(l.name) : `<b>${esc(l.name)}</b>`).join(', ')));
      card.appendChild(el('div', 'steps', esc(t.route.slice().reverse().join(' → '))));
      const b = el('button', 'plan', 'Plan it'); b.addEventListener('click', () => { document.getElementById('retroTarget').value = t.smiles; document.querySelector('.tab[data-page="retro"]').click(); window.Retro.plan(t.smiles); });
      card.appendChild(b); box.appendChild(card);
    });
    return box;
  }
  /* ================= wiring ================= */
  const SAMPLE = ['phenol', 'acetic acid', 'toluene', 'benzaldehyde', 'ethanol', 'cyclohexene', 'NaBH4', 'HBr', 'Br2', 'KMnO4'];
  window.Sketcher.attachSuggest(addBox, document.getElementById('possSuggest'), add, () => {
    const reag = C.RULES.filter(r => !String(r.id).startsWith('mol:')).map(r => ({ label: r.label, formula: 'reagent', k: r.label.toLowerCase() }));
    return window.Sketcher.nameIndex().concat(reag);
  });
  document.getElementById('possAddBtn').addEventListener('click', () => add(addBox.value));
  document.getElementById('possSample').addEventListener('click', () => { SAMPLE.forEach(s => { if (!shelf.some(x => x.text === s)) shelf.push({ text: s }); }); save(); });
  document.getElementById('possClear').addEventListener('click', () => { shelf = []; save(); });
  document.getElementById('possMax').addEventListener('change', refresh);
  document.getElementById('possFilter').addEventListener('input', () => { clearTimeout(refresh._t); refresh._t = setTimeout(refresh, 250); });
  renderShelf();
  /* start the one-off analysis as soon as the tab is first shown */
  document.querySelector('.tab[data-page="possib"]').addEventListener('click', () => { if (!targets) buildTargets(refresh); refresh(); }, { once: false });
  if (location.hash === '#possib') { buildTargets(refresh); refresh(); }
  window.Possibilities = { refresh, add };
})();
