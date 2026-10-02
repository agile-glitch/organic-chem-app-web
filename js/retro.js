/* Retrosynthesis page — type a target, get a plan back to purchasable starting materials.

   Everything chemical comes from js/chem.js:
     searchMolecule(text)          name / SMILES / CAS → graph
     retroAnalyze(g)               the disconnection tree: {g, name, transform{name,reagent,note}, children} | {basic:true}
     isBasicMaterial(g)            the "you can buy this" rule
     buyableMatches(g, n)          catalogue entries that match a piece exactly (or contain it)
     retroAlternatives(g)          every other disconnection that fits the target
     retroSearchRoutes(g, opts)    a Monte-Carlo tree search for whole routes (slower; run after the page is up)
   The page turns the tree into: a verdict, the shopping list, the forward steps in order, the backward
   tree, the alternatives, and the searched routes.                                                     */
(() => {
  'use strict';
  const C = window.Chem, M = window.MolDraw;
  const input = document.getElementById('retroTarget');
  const out = document.getElementById('retroOut');
  if (!input || !C) return;

  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const el = (tag, cls, html) => { const n = document.createElement(tag); if (cls) n.className = cls; if (html != null) n.innerHTML = html; return n; };
  const nameOf = g => { try { const n = C.displayNameFor(g); if (n) return n; } catch (e) {} try { return C.formula(g); } catch (e) { return '?'; } };
  const iupacOf = g => { try { const n = C.iupacName(g); return C.nameFailed && C.nameFailed(n) ? '' : n; } catch (e) { return ''; } };
  const formulaOf = g => { try { return C.formula(g); } catch (e) { return ''; } };
  const smilesOf = g => { try { return C.toSmiles(g); } catch (e) { return ''; } };
  /* the SMILES shown on a card: RDKit's canonical form (aromatic rings in lowercase, as SMILES are usually written)
     once RDKit has started, else the app's own (Kekulé) SMILES, which describes the same structure */
  let RD = null;
  const smilesShown = g => {
    const s = smilesOf(g);
    if (!RD || !s) return s;
    let m = null;
    try { m = RD.get_mol(s); return m && m.get_num_atoms() ? m.get_smiles() : s; } catch (e) { return s; } finally { if (m) m.delete(); }
  };

  /* a picture + name + tag, used everywhere on the page */
  function molCard(g, tag, extra) {
    const card = el('div', 'molcard');
    card.appendChild(M.svg(g, 170, 100));
    const common = nameOf(g), sys = iupacOf(g);
    const smi = smilesShown(g);
    card.appendChild(el('div', 'nm', esc(common) + (sys && sys !== common ? `<small>${esc(sys)}</small>` : '') +
      (smi ? `<small class="smi" title="SMILES">${esc(smi)}</small>` : '') + `<small>${esc(formulaOf(g))}</small>`));
    if (tag) card.appendChild(el('span', 'tag ' + tag.cls, tag.text));
    if (extra) card.appendChild(extra);
    return card;
  }

  /* ---- how a starting material is made or obtained: PubChem's "Methods of Manufacturing" (from the Hazardous
     Substances Data Bank, each entry with the publication it cites). Asked only on a click; the text is PubChem's,
     shown with its sources, not written by this app. ---- */
  const mfgPanel = () => { const p = el('div', 'card mfg'); p.hidden = true; return p; };
  function mfgButton(g, panel, label) {
    const b = el('button', 'mfg-btn', esc(label));
    b.title = 'Industrial production and natural sources, from PubChem. Sends this structure (as SMILES) to PubChem (NCBI) over the internet.';
    b.addEventListener('click', () => howMade(g, panel, b));
    return b;
  }
  async function howMade(g, panel, btn) {
    const name = nameOf(g);
    panel.hidden = false;
    panel.innerHTML = `<h3>How ${esc(name)} is made</h3><p class="hint">Asking PubChem…</p>`;
    panel.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
    btn.disabled = true;
    let res, m = null;
    try {
      if (!window.PubChem3D || !window.PubChem3D.manufacturing) throw new Error('the PubChem module is not loaded');
      const R = RD || await window.RDKitLoad();
      m = R.get_mol(smilesOf(g));
      if (!m || !m.get_num_atoms()) throw new Error('RDKit could not read this structure');
      res = await window.PubChem3D.manufacturing(R, m);
    } catch (e) { res = { ok: false, error: 'The lookup failed: ' + (e && e.message || e) }; }
    finally { if (m) m.delete(); btn.disabled = false; }
    if (!res.ok) { panel.innerHTML = `<h3>How ${esc(name)} is made</h3><p class="hint">${esc(res.error)}</p>`; return; }
    const items = res.entries.map(x => `<li>${esc(x.text)}${x.cite ? `<div class="cite">Source quoted: ${esc(x.cite)}</div>` : ''}</li>`).join('');
    panel.innerHTML = `<h3>How ${esc(name)} is made <small>industrial production and natural sources</small></h3>` +
      `<ol class="mfg-list">${items}</ol>` +
      `<p class="hint">"Methods of Manufacturing", ${esc(res.source)}, via PubChem CID ${res.cid} (${esc(res.title)}) · ` +
      `<a href="${esc(res.url)}" target="_blank" rel="noopener">open the record</a> · NCBI/NLM data — ` +
      `<a href="${esc(res.policyUrl)}" target="_blank" rel="noopener">NCBI policies and disclaimer</a>. ` +
      `Descriptions of industrial processes, quoted as PubChem gives them; not checked by this app.</p>`;
  }

  /* ---- read the tree ---- */
  function walk(node, fn, depth = 0) { fn(node, depth); (node.children || []).forEach(c => walk(c, fn, depth + 1)); }
  function forwardSteps(node, steps = []) {
    if (!node.children || !node.children.length) return steps;
    node.children.forEach(c => forwardSteps(c, steps));      // make the pieces first
    steps.push({ product: node, reactants: node.children, t: node.transform || {} });
    return steps;
  }

  function plan(text) {
    text = (text || '').trim(); if (!text) return;
    const rh = document.getElementById('retroRouteOut'), RC0 = window.RetroChimeraPage;
    if (rh && !(RC0 && RC0.searching())) rh.innerHTML = '';     // routes belong to the target they were searched for
    let g;
    try { g = C.searchMolecule(text).graph; } catch (e) {
      out.innerHTML = `<div class="banner">Could not read "${esc(text)}": ${esc(e.message)}</div>`;
      // not in the curated library: the saved PubChem lookups and the broad offline library first (no internet)
      if (window.PCLib && !/[=#\[\]@]/.test(text) && (/[a-z]{3}/i.test(text) || /^\d{2,7}-\d\d-\d$/.test(text))) {
        window.PCLib.resolve(text).then(r => {
          if (r.ok) { document.getElementById('retroTarget').value = r.smiles; plan(r.smiles); return; }
          if (r.choices) out.firstChild.insertAdjacentHTML('beforeend', '<br>' + esc(r.error) + '<br>' +
            r.choices.map(c => `CID ${c.cid} ${esc(c.title)}: <code>${esc(c.smiles)}</code>`).join('<br>') + '<br>Paste the SMILES you mean.');
        });
      }
      // a name or CAS number the offline library does not know: offer PubChem (only on a click)
      const P = window.PubChem3D;
      if (P && P.lookupName && !/[=#\[\]@]/.test(text) && (/[a-z]{3}/i.test(text) || /^\d{2,7}-\d\d-\d$/.test(text))) {
        const b = el('button', 'mfg-btn', 'Look it up on PubChem');
        b.title = 'Sends the text to PubChem (NCBI) over the internet; the answer is kept in this browser.';
        b.addEventListener('click', async () => {
          b.disabled = true; b.textContent = 'Asking PubChem…';
          const res = await P.lookupName(text);
          if (res.ok) { document.getElementById('retroTarget').value = res.smiles; plan(res.smiles); return; }
          b.remove();
          out.firstChild.insertAdjacentHTML('beforeend', `<br>${esc(res.error)}` +
            (res.choices ? '<br>' + res.choices.map(c => `CID ${c.cid} ${esc(c.title)}: <code>${esc(c.smiles)}</code>`).join('<br>') : ''));
        });
        out.firstChild.append(' ', b);
      }
      return;
    }
    out.innerHTML = '<p class="searching">Working backwards…</p>';
    // start RDKit for the cards' SMILES (once per page; shared with the other tabs), waiting at most 5 s for it
    const rd = window.RDKitLoad ? window.RDKitLoad().then(R => { RD = R; }, () => {}) : Promise.resolve();
    Promise.race([rd, new Promise(r => setTimeout(r, 5000))]).then(() => setTimeout(() => show(g), 20));
  }

  function show(g) {
    out.innerHTML = '';
    let tree;
    try { tree = C.retroAnalyze(g); } catch (e) { out.innerHTML = `<div class="banner">The analysis failed: ${esc(e.message)}</div>`; return; }
    const r = renderPlan(out, g, tree);
    extras(g, tree);
    if (!r.basic) precedents(g, nameOf(g));
  }

  /* The verdict, the shopping list and the reactions in order, for any tree of the retroAnalyze shape
     {g, transform{name, reagent, note}, children} | {g, basic:true}. The Retrosynthesis page passes its own analysis;
     the RetroChimera tab passes routes it searched (js/retrochimera.js), with opts:
       deadTag          the tag on an end piece that is neither basic nor a catalogue item (default 'no taught route')
       deadVerdict(ns)  the banner for those pieces (default: the "honest verdict" below)
       stats            extra HTML lines for the verdict card
       stepBadge(s)     HTML put before a step's explanation (s = {product, reactants, t})
       arrowText(s)     the text over the arrow (null: the first clause of the reagent text, as here)
       depthNote        false to leave out the depth-limit banner
       heading          the verdict's title (default 'How to make …')
       recordedNote     false to leave out the "recorded syntheses … bottom of the page" line
     Returns {basic} so the caller knows whether the target itself is a starting material. */
  function renderPlan(host, g, tree, opts = {}) {
    const steps = forwardSteps(tree);
    const leaves = [], dead = [];
    /* an end piece with no taught route that is an exact item in the embedded catalogue is bought, not a dead end
       (the route search below counts it as in stock the same way) */
    const inCatalogue = n => { try { return C.buyableMatches(n.g, 2).some(x => x.exact); } catch (e) { return false; } };
    walk(tree, n => { if (n.basic) leaves.push(n); else if (!n.children || !n.children.length) { if (inCatalogue(n)) { n.catalogue = true; leaves.push(n); } else dead.push(n); } });
    const seen = new Set();
    const buyList = leaves.filter(n => { const k = smilesOf(n.g); if (seen.has(k)) return false; seen.add(k); return true; });
    let maxDepth = 0; walk(tree, (n, d) => { maxDepth = Math.max(maxDepth, d); });

    /* 1. verdict */
    const v = el('div', 'card verdict');
    v.appendChild(M.svg(g, 220, 130));
    const common = nameOf(g), sys = iupacOf(g);
    const txt = el('div');
    txt.appendChild(el('div', 'big', opts.heading || 'How to make ' + esc(common)));
    txt.appendChild(el('div', 'sub', esc([sys !== common ? sys : '', formulaOf(g), smilesOf(g)].filter(Boolean).join(' · '))));
    if (tree.basic) txt.appendChild(el('div', 'stat', 'This is itself a <b>starting material</b>: small and common enough to buy rather than make.'));
    else txt.appendChild(el('div', 'stat', `<b>${steps.length}</b> reaction${steps.length === 1 ? '' : 's'} · <b>${buyList.length}</b> thing${buyList.length === 1 ? '' : 's'} to buy` +
      (dead.length ? ` · <b style="color:#b3261e">${dead.length}</b> piece${dead.length === 1 ? '' : 's'} with no taught route` : '')));
    const P = window.Precedent, made = P ? P.exact(g) : [];
    if (P && opts.recordedNote !== false) txt.appendChild(el('div', 'stat', made.length ? `<b>${made.length}</b> recorded synthes${made.length === 1 ? 'is' : 'es'} in the reaction database — see the bottom of the page.` : 'No recorded synthesis in the reaction database (a sample; see the note at the bottom).'));
    if (opts.stats) txt.appendChild(el('div', 'stat', opts.stats));
    v.appendChild(txt);
    host.appendChild(v);
    if (dead.length) host.appendChild(el('div', 'banner', opts.deadVerdict ? opts.deadVerdict(dead) : '<b>Honest verdict:</b> the plan reaches ' + dead.map(n => esc(nameOf(n.g))).join(', ') +
      ' and knows no taught reaction that makes ' + (dead.length === 1 ? 'it' : 'them') + '. ' +
      (dead.length === 1 ? 'It is' : 'They are') + ' listed as something to buy; check a supplier catalogue.'));
    if (!tree.basic && opts.depthNote !== false && maxDepth >= 6) host.appendChild(el('div', 'banner', 'The analysis hit its depth limit (6 steps), so the deepest pieces may not be the simplest possible.'));
    if (tree.basic) {
      const panel = mfgPanel();
      txt.appendChild(mfgButton(g, panel, 'How is it made? (PubChem)'));
      host.appendChild(panel);
      return { basic: true };
    }

    /* 2. shopping list */
    host.appendChild(el('h2', null, 'What to buy <small>starting materials</small>'));
    const shop = el('div', 'shop'), shopPanel = mfgPanel();
    buyList.concat(dead).forEach(n => {
      let extra = null;
      try {
        const m = C.buyableMatches(n.g, 2).filter(x => x.exact);
        if (m.length && m[0].entry.cas) extra = el('div', 'cas', 'CAS ' + esc(m[0].entry.cas));
      } catch (e) {}
      const card = molCard(n.g, n.basic || n.catalogue ? { cls: 'buy', text: 'buy' } : { cls: 'dead', text: opts.deadTag || 'no taught route' }, extra);
      card.appendChild(mfgButton(n.g, shopPanel, 'How is it made?'));
      shop.appendChild(card);
    });
    host.appendChild(shop);
    host.appendChild(shopPanel);

    /* 3. the reactions, in the order you would run them */
    host.appendChild(el('h2', null, 'The reactions, in order <small>forward synthesis</small>'));
    const madeIn = new Map();
    steps.forEach((s, i) => {
      const card = el('div', 'card step');
      card.appendChild(el('div', 'num', String(i + 1)));
      const body = el('div');
      const eqn = el('div', 'eqn');
      s.reactants.forEach((r, k) => {
        if (k) eqn.appendChild(el('span', 'plus', '+'));
        const key = smilesOf(r.g);
        const tag = r.basic || r.catalogue ? { cls: 'buy', text: 'buy' } : madeIn.has(key) ? { cls: 'made', text: 'from step ' + madeIn.get(key) } : { cls: 'dead', text: 'buy' };
        eqn.appendChild(molCard(r.g, tag));
      });
      const arrow = el('div', 'arrow');
      const at = opts.arrowText ? opts.arrowText(s) : null;
      arrow.appendChild(el('div', null, esc(at != null ? at : shortReagent(s.t.reagent))));
      arrow.appendChild(el('div', 'ar', '⟶'));
      eqn.appendChild(arrow);
      const isTarget = s.product === tree;
      eqn.appendChild(molCard(s.product.g, isTarget ? { cls: 'target', text: 'target' } : { cls: 'made', text: 'step ' + (i + 1) }));
      madeIn.set(smilesOf(s.product.g), i + 1);
      body.appendChild(eqn);
      const nPrec = window.Precedent ? window.Precedent.exact(s.product.g).length : 0;
      body.appendChild(el('div', 'why', (opts.stepBadge ? opts.stepBadge(s) : '') + `<b>${esc(s.t.name || '')}.</b> ${esc(s.t.note || '')}` + (s.t.reagent ? `<br><b>Conditions:</b> ${esc(s.t.reagent)}` : '') +
        (nPrec ? `<br><b>Recorded:</b> ${nPrec} reaction${nPrec === 1 ? '' : 's'} in the database make${nPrec === 1 ? 's' : ''} this product.` : '')));
      // has this step been done? papers naming a reactant and the product (PubChem; not necessarily this reaction)
      if (window.LitRef) {
        const box = window.LitRef.panel(), prodName = nameOf(s.product.g), prodSmi = smilesShown(s.product.g);
        body.appendChild(window.LitRef.button('Has this been done? Papers mentioning both (PubChem)',
          () => s.reactants.map(r => ({ a: { smiles: smilesShown(r.g), name: nameOf(r.g) }, b: { smiles: prodSmi, name: prodName } })), box));
        body.appendChild(box);
      }
      card.appendChild(body);
      host.appendChild(card);
    });

    return { basic: false };
  }

  /* recorded syntheses of the target, from the reaction database */
  function precedents(g, name) {
    const P = window.Precedent; if (!P) return;
    const made = P.exact(g), skipped = P.skipped(g);
    out.appendChild(el('h2', null, 'Recorded syntheses <small>from the reaction database</small>'));
    const box = el('div');
    if (skipped) box.appendChild(el('p', 'hint', `${skipped} record${skipped === 1 ? '' : 's'} naming this molecule as the product ${skipped === 1 ? 'was' : 'were'} skipped as data errors (the recorded product is a solvent or a fragment of a much larger reactant).`));
    if (made.length) {
      box.appendChild(el('p', null, `<b>${made.length}</b> recorded reaction${made.length === 1 ? '' : 's'} make${made.length === 1 ? 's' : ''} <b>${esc(name)}</b>` + (made.length > 5 ? ' — the first 5:' : ':')));
      made.slice(0, 5).forEach(x => box.appendChild(P.card(x)));
    } else {
      box.appendChild(el('p', null, `No recorded reaction in the database makes <b>${esc(name)}</b>.`));
      const sim = P.similar(g, 3);
      if (sim.length) { box.appendChild(el('p', 'hint', 'Closest recorded products:')); sim.forEach(s => { const c = P.card(s.rows[0]); c.prepend(el('div', 'prec-sim', Math.round(s.sim * 100) + '% similar')); box.appendChild(c); }); }
    }
    box.appendChild(el('div', 'prec-note', P.NOTE));
    out.appendChild(box);
  }

  /* the first clause of a reagent description, for the arrow label */
  const shortReagent = r => { if (!r) return ''; const s = r.split(/[—;(]/)[0].trim(); return s.length > 48 ? s.slice(0, 46) + '…' : s; };

  /* 4. the backward tree, 5. alternatives at the target, 6. searched routes */
  function extras(g, tree) {
    if (!tree.basic) {
      const d = el('details'); d.appendChild(el('summary', null, 'Backward analysis: the disconnection tree'));
      const box = el('div', 'tree');
      const build = node => {
        const li = el('li');
        const row = el('div', 'node');
        row.appendChild(M.svg(node.g, 120, 72));
        const how = node.transform ? `<b>${esc(node.transform.name)}</b> — ${esc(node.transform.note || '')}` : (node.basic ? '<b>buy</b>' : node.catalogue ? '<b>buy</b> (catalogue item)' : '<b>no taught disconnection</b>');
        row.appendChild(el('div', 'how', `<div>${esc(nameOf(node.g))}</div>${how}`));
        li.appendChild(row);
        if (node.children && node.children.length) { const ul = el('ul'); node.children.forEach(c => ul.appendChild(build(c))); li.appendChild(ul); }
        return li;
      };
      const ul = el('ul'); ul.appendChild(build(tree)); box.appendChild(ul);
      d.appendChild(box); out.appendChild(d);
    }

    let alts = [];
    try { alts = C.retroAlternatives(g); } catch (e) {}
    if (alts.length > 1 || (alts.length === 1 && !tree.transform)) {
      const d = el('details'); d.appendChild(el('summary', null, `Other ways to disconnect the target (${alts.length})`));
      d.appendChild(el('ul', 'alts', alts.map(a => `<li><b>${esc(a.name)}</b> → ${a.pieces.map(esc).join(' + ')}<br><span style="color:#666">${esc(a.reagent || '')}</span></li>`).join('')));
      out.appendChild(d);
    }

    /* the route search is the slow part (a few seconds): run it after the page has painted */
    const d = el('details'); d.open = true;
    d.appendChild(el('summary', null, 'Route search <small style="font-weight:400;color:#777">Monte-Carlo tree search over all taught reactions</small>'));
    const body = el('div', 'searching', 'Searching for complete routes…'); d.appendChild(body); out.appendChild(d);
    setTimeout(() => {
      let routes = [];
      try { routes = C.retroSearchRoutes(g, { iterations: 60 }); } catch (e) { body.textContent = 'The search failed: ' + e.message; return; }
      if (!routes.length) { body.className = ''; body.innerHTML = '<span class="hint">No route was found in which every starting material is an exact catalogue item. ' +
        'The plan above still stands; its leaves are judged by the "small and common" rule rather than by catalogue hits.</span>'; return; }
      body.className = '';
      const ol = el('ol', 'routes');
      routes.slice(0, 8).forEach(r => {
        const stepsHtml = r.steps.map(s => `${esc(s.from.name)} ⇐ ${s.to.map(t => `<span class="${t.stock ? 'st' : ''}">${esc(t.name)}${t.stock ? ' ✓' : ''}</span>`).join(' + ')} <span class="sc">(${esc(s.name)})</span>`).join('<br>');
        ol.appendChild(el('li', null, `<span class="sc">${r.steps.length} step${r.steps.length === 1 ? '' : 's'} · score ${r.score.toFixed(2)}</span><br>${stepsHtml}`));
      });
      body.appendChild(ol);
      body.appendChild(el('div', 'hint', '✓ = an exact catalogue item. Score: 0.95 × fraction of pieces in stock + a small bonus for short routes.'));
    }, 30);
  }

  /* ---- wiring ---- */
  window.Sketcher.attachSuggest(input, document.getElementById('retroSuggest'), plan);
  document.getElementById('retroGo').addEventListener('click', () => plan(input.value));
  document.getElementById('retroFromSketch').addEventListener('click', () => {
    const s = window.Sketcher.smiles();
    if (!s) { out.innerHTML = '<div class="banner">Nothing is drawn on the Sketcher tab yet.</div>'; return; }
    input.value = s; plan(s);
  });
  /* Find whole routes: the same search, output and evidence labels as the RetroChimera tab (js/retrochimera.js owns
     it; this tab supplies its own target box, button and output area). One search runs at a time. */
  const routeBtn = document.getElementById('retroRoutes');
  if (routeBtn) routeBtn.addEventListener('click', () => {
    const RC = window.RetroChimeraPage, host = document.getElementById('retroRouteOut');
    if (!RC || !RC.findRoutes) { host.innerHTML = '<div class="banner">The RetroChimera route search is not loaded.</div>'; return; }
    if (RC.searching()) {                                      // one search at a time: this stops whichever is running
      RC.stopSearch();
      if (RC.searchingOn && RC.searchingOn() !== routeBtn)
        host.innerHTML = '<p class="hint">A route search started on the RetroChimera tab was still running; it has been stopped. Click <b>Find whole routes</b> again to search for this target.</p>';
      return;
    }
    RC.findRoutes({ host, btn: routeBtn, text: input.value, sug: document.getElementById('retroSuggest') });
  });

  window.Retro = { plan, renderPlan, useRDKit: R => { RD = R; } };
})();
