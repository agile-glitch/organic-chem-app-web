/* Validation — the reaction engine graded against a recorded-reaction set.

   Method (rebuilt from organic-chem-tool/docs/finding-the-next-gap.md and
   reaction-database-and-modern-reactions.md):
     For each recorded reaction, offer the engine every sensible way in:
       * the largest reactant with every taught reagent the recorded agents map to
         (Reagents.resolve on the agent's SMILES / name, plus a small structure map);
       * each ordered pair of reactants (A, B) through the molecule + molecule path
         (Reagents.moleculeRule(B) → RXN.run(A, 'mol:…')).
     Then one strict question, by canonical key: is the engine's MAJOR product the recorded product?
     Grading stops at the first route that is right (nothing beats "right"), which is only faster.
   Per reaction the verdict is one of:
     right        a route's major product is the recorded product
     ranked       the recorded product appears among a route's outcomes, but not as the major one
     wrong        the engine offered products, none of them the recorded one
     refused      every route refused with a reason (no live site, gate, no mapped reagent)
     data_error   the recorded product is a solvent or a fragment of a much larger reactant
     unreadable   a molecule the engine cannot parse (an element outside its valence table)
   Sets: the built-in database, plus any public set loaded through js/sources.js. Runs in chunks so the
   page stays responsive; each set's result is kept in localStorage and can be exported.            */
(() => {
  'use strict';
  const C = window.Chem, R = window.Reagents, P = window.Precedent, S = window.Sources, NAMES = window.NAMES || { solvents: {}, agents: {} };
  const out = document.getElementById('validOut'), status = document.getElementById('validStatus'), bar = document.querySelector('#validBar div');
  const setSel = document.getElementById('validSet'), srcBox = document.getElementById('validSources');
  if (!out || !C || !R || !P || typeof RXN === 'undefined') return;
  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const el = (tag, cls, html) => { const n = document.createElement(tag); if (cls) n.className = cls; if (html != null) n.innerHTML = html; return n; };
  const store = { get(k) { try { return JSON.parse(localStorage.getItem(k)); } catch (e) { return null; } }, set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) {} }, del(k) { try { localStorage.removeItem(k); } catch (e) {} } };

  /* ---- the set under test ---- */
  const DATA = { db: window.REACTIONS };
  let D = DATA.db, M = D.molecules, RX = D.reactions, setId = 'db';
  const modeSel = document.getElementById('validMode');
  const mode = () => (modeSel ? modeSel.value : 'rules');
  const KEY = () => 'bench.v2.' + setId + '.' + mode();
  function chooseSet(id) {
    if (running) stop();
    setId = id; D = DATA[id]; M = D.molecules; RX = D.reactions;
    graphCache.clear(); agentCache.clear();
    state = loadState();
    document.getElementById('validRun').textContent = state && state.i && state.i < RX.length ? 'Resume' : 'Run the benchmark';
    document.getElementById('validRun').disabled = false;
    bar.style.width = state && state.i ? Math.round(100 * state.i / RX.length) + '%' : '0';
    status.textContent = '';
    render();
  }
  function addSetOption(id, data, title) {
    DATA[id] = data;
    if (![...setSel.options].some(o => o.value === id)) { const o = document.createElement('option'); o.value = id; o.textContent = `${title} (${data.reactions.length.toLocaleString()})`; setSel.appendChild(o); }
  }
  function renderSources() {
    srcBox.innerHTML = '';
    (S ? S.SETS : []).forEach(set => {
      const on = S.loaded().includes(set.id);
      const row = el('div', 'srcrow');
      row.innerHTML = `<b>${esc(set.title)}</b> <span class="hint">[${esc(set.license)}]</span><div class="hint">${esc(set.note)}</div>`;
      if (on) { const d = S.get(set.id); row.appendChild(el('div', 'hint', `Loaded: ${d.reactions.length.toLocaleString()} reactions over ${d.molecules.length.toLocaleString()} molecules` + (d.unreadable ? ` (${d.unreadable} molecules the engine cannot read)` : '') + (d.citation ? ' · ' + esc(d.citation) : ''))); }
      else {
        const b = el('button', 'use', 'Load'); b.style.marginTop = '4px';
        b.addEventListener('click', () => { b.disabled = true; b.textContent = 'Loading…'; S.load(set.id).then(() => renderSources()).catch(e => { b.textContent = 'Load'; b.disabled = false; row.appendChild(el('div', 'hint', '⚠ ' + esc(e.message))); }); });
        row.appendChild(b);
      }
      srcBox.appendChild(row);
    });
    if (!S) srcBox.appendChild(el('p', 'hint', 'Sources are not available.'));
  }
  if (S) { S.loaded().forEach(id => addSetOption(id, S.get(id), S.SETS.find(s => s.id === id).title)); S.onLoad((id, data) => { addSetOption(id, data, S.SETS.find(s => s.id === id).title); renderSources(); }); }
  setSel.addEventListener('change', () => chooseSet(setSel.value));
  if (modeSel) modeSel.addEventListener('change', () => chooseSet(setId));          // each predictor keeps its own saved result

  /* ---- mapping a recorded agent to a taught reagent ---- */
  const AGENT_MAP = {                                    // structures the name/formula tables do not catch
    '[BH4-]': 'nabh4', '[Na+].[BH4-]': 'nabh4', '[AlH4-]': 'lialh4', '[Li+].[AlH4-]': 'lialh4', '[H][H]': 'h2', '[HH]': 'h2',
    'O=[Mn](=O)(=O)[O-]': 'kmno4', 'O=[Cr](=O)=O': 'jones', 'O=C(OO)c1cccc(Cl)c1': 'mcpba', 'OOC(=O)c1cccc(Cl)c1': 'mcpba',
    'O=[N+]([O-])O': 'nitration', 'BrBr': 'br2', 'ClCl': 'cl2', 'Cl': 'hcl', 'Br': 'hbr', 'I': 'hi', '[OH-]': 'naoh', '[Na+].[OH-]': 'naoh', '[K+].[OH-]': 'naoh',
    'CC(C)(C)[O-]': 'tbuok', 'CC[O-]': 'naoet', 'C[O-]': 'naoet', '[C-]#N': 'nacn', 'N#[C-]': 'nacn', 'O=S(Cl)Cl': 'soCl2', 'O=C(Cl)C(=O)Cl': 'soCl2',
    'ClP(Cl)Cl': 'pcl5', 'ClP(Cl)(Cl)(Cl)Cl': 'pcl5', 'P(Br)(Br)Br': 'pbr3', 'BrP(Br)Br': 'pbr3', 'O=N[O-]': 'nano2', '[O-]N=O': 'nano2',
    'O=C1CCC(=O)N1Br': 'allylbr', 'CC(C)C[AlH]CC(C)C': 'dibal', 'O=C(O)C(F)(F)F': 'boc_off', 'CC(=O)OC(C)=O': 'ac_on',
    'CC(C)(C)OC(=O)OC(=O)OC(C)(C)C': 'boc_on', 'CC(C)(C)[Si](C)(C)Cl': 'tbs_on', '[F-]': 'tbs_off', 'CCCC[N+](CCCC)(CCCC)CCCC.[F-]': 'tbs_off',
    'BrCc1ccccc1': 'bn_on', 'ClCc1ccccc1': 'bn_on', '[Fe]': 'nitrored', 'Cl[Sn]Cl': 'nitrored', '[Sn]': 'nitrored', '[Zn]': 'nitrored',
    'NN': 'wolffkishner', 'O=S(=O)(O)O': 'aqacid', 'OS(=O)(=O)O': 'aqacid', 'O=[Os](=O)(=O)=O': 'oso4', 'B': 'hydroboration', 'B1CCCC1': 'hydroboration',
    'CS(=O)(=O)Cl': 'ms_on', 'Cc1ccc(S(=O)(=O)Cl)cc1': 'ts_on', 'C1=COCCC1': 'thp_on', 'O=[Mn]=O': 'mno2', 'O=[Ag]': 'ag2o',
    'Cc1ccc(S(=O)(=O)O)cc1': 'aqacid', 'O=C(OOC(=O)c1ccccc1)c1ccccc1': 'radpoly', 'CC(C)(C#N)N=NC(C)(C)C#N': 'radpoly',
  };
  /* a bare metal in the record means a hydrogenation whose H₂ was not written down: it can be an alkene
     reduction, a nitro reduction or a hydrogenolysis, so all of those doors are offered */
  const METAL_H2 = ['h2', 'nitrored', 'cbz_off', 'bn_off'];
  ['[Pd]', '[Pt]', '[Ni]', '[Rh]', 'O=[Pt]=O', '[Pd].[C]', '[C].[Pd]'].forEach(s => { AGENT_MAP[s] = METAL_H2; });
  Object.assign(AGENT_MAP, { 'O=S(=O)(O)C(F)(F)F': 'boc_off', '[Li]CCCC': 'lda', 'CC(C)[N-]C(C)C': 'lda', 'O=C1CCC(=O)N1Cl': 'allylbr',
    'OO': 'mcpba', 'CC(=O)OO': 'mcpba', 'O=C(O)OO': 'mcpba', 'NO': 'oxime', '[N-]=[N+]=[N-]': 'azide' });
  const agentCache = new Map();
  function rulesForAgent(i) {
    if (agentCache.has(i)) return agentCache.get(i);
    const s = M[i].s; let ids = [];
    if (!NAMES.solvents[s]) {                              // a solvent is never the reagent
      if (AGENT_MAP[s]) ids = ids.concat(AGENT_MAP[s]).filter(id => C.RULES.some(r => r.id === id));
      else {
        window.REAGENT_SUBSTRATE = null;
        for (const text of [NAMES.agents[s], s].filter(Boolean)) {
          let r = null; try { r = R.resolve(text, true); } catch (e) {}
          if (r && r.id && !String(r.id).startsWith('mol:')) { ids.push(r.id); break; }
        }
      }
      ids = ids.filter(id => id !== 'hydration');          // water in the flask is not, by itself, a hydration
    }
    agentCache.set(i, ids); return ids;
  }
  const graphCache = new Map();
  const graphOf = i => { if (!graphCache.has(i)) { let g = null; try { g = C.parseSmiles(M[i].s); } catch (e) {} graphCache.set(i, g); } return graphCache.get(i); };
  const keyOf = g => { try { return C.libKeyHash(C.canonicalKey(g)); } catch (e) { return null; } };
  const heavy = i => (M[i].f.match(/[A-Z][a-z]?\d*/g) || []).reduce((n, t) => n + (t[0] === 'H' && !/^H[a-z]/.test(t) ? 0 : (parseInt(t.replace(/\D/g, ''), 10) || 1)), 0);

  /* ---- grade one recorded reaction ---- */
  function grade(x) {
    const want = M[x.p].k;
    if (!P.plausible(x)) return { v: 'data_error' };
    if (!want || x.r.some(i => !M[i].k)) return { v: 'unreadable' };
    const reactants = x.r.filter(i => graphOf(i));
    if (!reactants.length) return { v: 'unreadable' };
    const routes = [];                                     // cheapest first: a written reagent ~6 ms, a molecule pair ~110 ms
    /* The records file a reagent that contributes atoms (HNO₃, SOCl₂, NBS) under REACTANTS, and one that does
       not (a base, a catalyst) under agents — so taught reagents are looked for in both lists. The substrate is
       then a reactant that is not itself a reagent: the two largest are tried. */
    const ids = new Set(); x.a.concat(x.r).forEach(i => rulesForAgent(i).forEach(id => ids.add(id)));
    const organic = reactants.filter(i => !rulesForAgent(i).length);
    const subs = (organic.length ? organic : reactants).slice().sort((p, q) => heavy(q) - heavy(p)).slice(0, 2);
    const X = window.RXX;                                 // the extension families are more specific, so they go first
    for (const big of subs) for (const id of ids) {
      const rule = C.RULES.find(r => r.id === id);
      if (X) routes.push({ sub: big, via: id, label: (rule ? rule.label : id) + ' [ext]', run: () => X.runReagent(graphOf(big), id) });
      routes.push({ sub: big, via: id, label: rule ? rule.label : id, run: () => RXN.run(graphOf(big), id, {}) });
    }
    if (reactants.length >= 2) for (const a of reactants) for (const b of reactants) {
      if (a === b) continue;
      if (X) routes.push({ sub: a, via: 'mol:' + b, label: M[b].f + ' (as partner) [ext]', run: () => X.run(graphOf(a), graphOf(b)) });
      routes.push({ sub: a, via: 'mol:' + b, label: M[b].f + ' (as partner)', run: () => { const rule = R.moleculeRule(graphOf(b), M[b].s); return rule ? RXN.run(graphOf(a), rule.id, {}) : null; } });
    }
    if (!routes.length) return { v: 'refused', why: 'no reagent mapped and only one reactant' };
    /* `best` is the most favourable verdict over all routes (the harness's generous reading); `first` is the answer
       of the first route that produced a product — what a user would be shown, so the fair basis for combining */
    const rank = { right: 3, ranked: 2, wrong: 1, refused: 0 };
    let best = null, first = null;
    const rulesBefore = C.RULES.length;
    try {
      for (const rt of routes) {
        let res = null; try { res = rt.run(); } catch (e) { res = null; }
        let v = 'refused', got = null, why = '';
        if (res && !res.error && res.major) {
          got = res.major.smiles;
          if (keyOf(res.major.graph) === want) v = 'right';
          else if (res.outcomes.some(o => keyOf(o.graph) === want)) v = 'ranked';
          else v = 'wrong';
          if (!first) first = { v, smiles: got, label: rt.label, ext: !!res.ext, via: rt.via, sub: rt.sub };
        } else why = res && res.error ? res.error : 'no result';
        if (!best || rank[v] > rank[best.v]) best = { v, via: rt.via, label: rt.label, sub: rt.sub, got, why };
        if (v === 'right' && first) break;
      }
    } finally { C.RULES.length = rulesBefore; }              // drop the temporary mol: rules
    /* the verdict is the COMMITTED answer's; the generous reading is kept as vAny, and shown only as a note */
    if (first) return { v: first.v, vAny: best.v, via: first.via, label: first.label, sub: first.sub, got: first.smiles, why: '', first };
    return { ...best, vAny: best.v, first: null };
  }

  /* ---- grading the precedent-based predictor (js/rxnkb.js), alone or combined with the rules as on the page ----
     The recorded product is compared as RDKit canonical SMILES without stereo. `conf` keeps the confidence the
     predictor stated, so the run can check that "High" really is right more often than "Low". */
  /* The arbitration the Reactions page uses (js/reactions.js, decide()), repeated here so the benchmark measures
     exactly what a user sees:
       - no calibrated confidence, or below RELIABLE → the precedents are not relied on; the mechanistic answer, if any
       - the mechanistic model's committed answer equals the precedent's first → that answer ("agree")
       - it is among the precedent's top 3 → the mechanistic answer ("pick")
       - otherwise → the precedent's first                                                                        */
  const RELIABLE = 0.2;
  function gradeKB(x, combined) {
    const K = window.RXNKB;
    if (!P.plausible(x)) return { v: 'data_error' };
    const want = K.canonical(M[x.p].s);
    if (!want) return { v: 'unreadable' };
    const res = K.predict({ reactants: x.r.map(i => M[i].s), agents: x.a.map(i => M[i].s) });
    if (res.error) return { v: 'unreadable' };
    const c = res.candidates, top = c[0] || null;
    const conf = top && top.confidence.p != null ? top.confidence.p : null;
    const verdict = flat => (flat === want.flat ? 'right' : c.slice(0, 5).some(y => y.flat === want.flat) ? 'ranked' : 'wrong');
    if (!combined) {
      if (!top) return { v: 'refused', why: 'none of the most similar recorded reactions applies to these molecules', conf: null };
      return { v: verdict(top.flat), via: 'kb', label: 'precedents (closest ' + Math.round(100 * top.best) + ' %)', got: top.smiles, conf, sub: x.r[0] };
    }
    const g = grade(x), rf = g.first ? K.canonical(g.first.smiles) : null, rfFlat = rf ? rf.flat : null;
    const reliable = top && conf != null && conf >= RELIABLE;
    let answer = null, smiles = null, how = '';
    if (reliable && rfFlat && rfFlat === top.flat) { answer = top.flat; smiles = top.smiles; how = 'agree'; }
    else if (reliable && rfFlat && c.slice(0, 3).some(y => y.flat === rfFlat)) { answer = rfFlat; smiles = g.first.smiles; how = 'pick'; }
    else if (reliable) { answer = top.flat; smiles = top.smiles; how = 'precedent'; }
    else if (rfFlat) { answer = rfFlat; smiles = g.first.smiles; how = 'rules'; }
    if (!answer) return { v: 'refused', why: top ? 'precedents too distant to rely on, and the mechanistic model has no answer' : 'no precedent applies, and the mechanistic model has no answer', conf: null };
    return { v: verdict(answer), via: how === 'rules' || how === 'pick' ? g.via : 'kb', label: { agree: 'both methods agree', pick: 'mechanistic pick among the precedent top 3', precedent: 'precedents', rules: 'mechanistic model (no close precedent)' }[how],
      got: smiles, conf: how === 'precedent' ? conf : null, how, sub: g.sub != null ? g.sub : x.r[0] };
  }
  const gradeBy = x => (mode() === 'kb' ? gradeKB(x, false) : mode() === 'combined' ? gradeKB(x, true) : grade(x));

  /* ---- the run: chunks, progress, persistence ---- */
  /* results saved by earlier versions (key bench.v1.*) graded the rule engine generously — right if ANY route was
     right — so they are discarded rather than shown under the committed-answer grading used now */
  try { Object.keys(localStorage).filter(k => k.startsWith('bench.v1')).forEach(k => localStorage.removeItem(k)); } catch (e) {}
  function loadState() { return store.get(KEY()) || null; }
  let state = loadState();                                 // {i, rows: {id: {v, via, label, sub, got, why, conf}}}
  let running = false;
  async function run() {
    if (running) return; running = true;
    document.getElementById('validRun').disabled = true; document.getElementById('validStop').disabled = false;
    if (mode() !== 'rules') {
      try { await window.RXNKB.load(msg => { status.textContent = msg; }); }
      catch (e) { running = false; status.textContent = 'The precedent knowledge base is not installed (data/rxnkb/).'; document.getElementById('validRun').disabled = false; document.getElementById('validStop').disabled = true; return; }
    }
    if (!state) state = { i: 0, rows: {}, started: new Date().toISOString(), set: setId, mode: mode() };
    const t0 = performance.now(), i0 = state.i, myKey = KEY();
    (function step() {
      if (!running || myKey !== KEY()) return;
      const end = Math.min(RX.length, state.i + 4);
      for (; state.i < end; state.i++) {
        const x = RX[state.i]; let r;
        try { r = gradeBy(x); } catch (e) { r = { v: 'refused', why: 'harness error: ' + e.message }; }
        state.rows[x.id] = r;
      }
      bar.style.width = Math.round(100 * state.i / RX.length) + '%';
      const rate = (state.i - i0) / ((performance.now() - t0) / 1000);
      status.textContent = `${state.i.toLocaleString()} / ${RX.length.toLocaleString()} graded` + (rate > 0 ? ` · ~${Math.round((RX.length - state.i) / rate)} s left` : '');
      if (state.i % 200 === 0) { store.set(myKey, state); render(); }
      if (state.i < RX.length) setTimeout(step, 0);
      else { running = false; state.finished = new Date().toISOString(); store.set(myKey, state); status.textContent = 'Complete.'; document.getElementById('validStop').disabled = true; render(); }
    })();
  }
  function stop() { running = false; store.set(KEY(), state); document.getElementById('validRun').disabled = false; document.getElementById('validStop').disabled = true; document.getElementById('validRun').textContent = 'Resume'; render(); }

  /* ---- results ---- */
  const pct = (n, d) => (d ? (100 * n / d).toFixed(1) + ' %' : '—');
  function summarise(rows, filter) {
    const s = { n: 0, right: 0, ranked: 0, wrong: 0, refused: 0, data_error: 0, unreadable: 0 };
    for (const id in rows) { const x = RX[id]; if (!x || (filter && !filter(x))) continue; s.n++; s[rows[id].v]++; }
    s.gradable = s.n - s.data_error - s.unreadable;
    s.attempted = s.right + s.ranked + s.wrong;
    return s;
  }
  function render() {
    out.innerHTML = '';
    if (!state || !state.i) { out.appendChild(el('p', 'hint', 'Not run yet for this set.')); return; }
    const rows = state.rows, all = summarise(rows);
    const kpi = (lbl, big, sub) => `<div class="kpi"><div class="lbl">${lbl}</div><div class="big">${big}</div><div class="hint">${sub || ''}</div></div>`;
    out.appendChild(el('div', 'kpis',
      kpi('Major product right', pct(all.right, all.gradable), `${all.right.toLocaleString()} of ${all.gradable.toLocaleString()} gradable reactions`) +
      kpi('Right, of those attempted', pct(all.right, all.attempted), `${all.attempted.toLocaleString()} attempted`) +
      kpi('Engine offered a product', pct(all.attempted, all.gradable), `${all.refused.toLocaleString()} refused with a reason`) +
      kpi('Recorded product in its ranking', pct(all.right + all.ranked, all.gradable), `${all.ranked.toLocaleString()} ranked but not major`) +
      kpi('Set aside', (all.data_error + all.unreadable).toLocaleString(), `${all.data_error} data errors · ${all.unreadable} unreadable`)));
    out.appendChild(el('p', 'hint', `${state.i.toLocaleString()} of ${RX.length.toLocaleString()} reactions graded in ${esc(setSel.selectedOptions[0].textContent)}` + (state.finished ? ` · finished ${state.finished.slice(0, 16).replace('T', ' ')}` : ' (partial)') +
      '. "Right" is measured on gradable reactions, so the denominator includes every reaction the engine refused, and it grades the one answer the predictor commits to.'));
    const anyRight = Object.values(rows).filter(r => r.vAny === 'right').length;
    if (Object.values(rows).some(r => r.vAny)) out.appendChild(el('p', 'hint', `For comparison, the generous reading (right if ANY route the rule engine tried was right, which a user never sees): ${pct(anyRight, all.gradable)}.`));
    const hows = {}; Object.values(rows).forEach(r => { if (r.how) { hows[r.how] = hows[r.how] || [0, 0]; hows[r.how][0]++; hows[r.how][1] += r.v === 'right'; } });
    if (Object.keys(hows).length) out.appendChild(el('p', 'hint', 'How the combined answer was reached: ' + Object.entries(hows).map(([h, [n, ok]]) => `${esc({ agree: 'both methods agree', pick: 'mechanistic pick among the precedent top 3', precedent: 'precedents alone', rules: 'mechanistic alone (no close precedent)' }[h] || h)}: ${ok} of ${n} right (${pct(ok, n)})`).join(' · ') + '.'));
    const table = (title, groups) => {
      const t = el('table'); t.innerHTML = `<thead><tr><th>${title}</th><th class="num">gradable</th><th class="num">right</th><th class="num">right %</th><th class="num">ranked</th><th class="num">wrong</th><th class="num">refused</th><th class="num">right of attempted</th></tr></thead>`;
      const tb = el('tbody');
      groups.forEach(([name, s]) => { if (!s.gradable) return; tb.innerHTML += `<tr><td>${esc(name)}</td><td class="num">${s.gradable}</td><td class="num">${s.right}</td><td class="num">${pct(s.right, s.gradable)}</td><td class="num">${s.ranked}</td><td class="num">${s.wrong}</td><td class="num">${s.refused}</td><td class="num">${pct(s.right, s.attempted)}</td></tr>`; });
      t.appendChild(tb); return t;
    };
    /* does the stated confidence mean anything? first answers grouped by the label the predictor gave them */
    if (Object.values(rows).some(r => r.conf != null)) {
      out.appendChild(el('h2', null, 'Is the stated confidence honest? <small>first answers, by the label shown on the Reactions page</small>'));
      const bands = [['High (≥ 80 %)', 0.8, 1.01], ['Medium (50–80 %)', 0.5, 0.8], ['Low (25–50 %)', 0.25, 0.5], ['Very low (< 25 %)', 0, 0.25]];
      const t = el('table'); t.innerHTML = '<thead><tr><th>Stated</th><th class="num">answers</th><th class="num">right</th><th class="num">measured here</th></tr></thead>';
      const tb = el('tbody');
      bands.forEach(([name, lo, hi]) => { let n = 0, ok = 0; for (const id in rows) { const r = rows[id]; if (r.conf == null || r.conf < lo || r.conf >= hi || !['right', 'ranked', 'wrong'].includes(r.v)) continue; n++; ok += r.v === 'right'; }
        tb.innerHTML += `<tr><td>${name}</td><td class="num">${n}</td><td class="num">${ok}</td><td class="num">${pct(ok, n)}</td></tr>`; });
      t.appendChild(tb); out.appendChild(t);
    }
    if (D.sources.length > 1) { out.appendChild(el('h2', null, 'By source')); out.appendChild(table('Source', D.sources.map((s, i) => [s, summarise(rows, x => x.src === i)]))); }
    if (D.classes.length > 1) {
      out.appendChild(el('h2', null, 'By reaction class <small>sorted by how many are wrong or refused — the to-do list</small>'));
      const cls = D.classes.map((c, i) => [c, summarise(rows, x => x.cls === i)]).filter(([, s]) => s.gradable).sort((a, b) => (b[1].wrong + b[1].refused) - (a[1].wrong + a[1].refused));
      out.appendChild(table('Class', cls));
    }
    if (RX.some(x => x.year)) {
      out.appendChild(el('h2', null, 'By decade'));
      const decades = [...new Set(RX.map(x => x.year && Math.floor(x.year / 10) * 10).filter(Boolean))].sort();
      out.appendChild(table('Decade', decades.map(d => [d + 's', summarise(rows, x => x.year && Math.floor(x.year / 10) * 10 === d)])));
    }

    /* misses */
    out.appendChild(el('h2', null, 'The misses <small>read the reactions, not just the labels</small>'));
    const f = el('div', 'filters');
    f.innerHTML = `<select id="vfV"><option value="">wrong + ranked + refused</option><option value="wrong">wrong</option><option value="ranked">ranked, not major</option><option value="refused">refused</option><option value="right">right (spot-check)</option></select>
      ${D.classes.length > 1 ? `<select id="vfC"><option value="">all classes</option>${D.classes.map((c, i) => `<option value="${i}">${esc(c)}</option>`).join('')}</select>` : ''}
      ${D.sources.length > 1 ? `<select id="vfS"><option value="">all sources</option>${D.sources.map((c, i) => `<option value="${i}">${esc(c)}</option>`).join('')}</select>` : ''}
      <input id="vfQ" placeholder="filter by reason, reagent or SMILES">`;
    out.appendChild(f);
    const list = el('div'); out.appendChild(list);
    const draw = () => {
      const v = document.getElementById('vfV').value, c = document.getElementById('vfC') ? document.getElementById('vfC').value : '', s = document.getElementById('vfS') ? document.getElementById('vfS').value : '', q = document.getElementById('vfQ').value.toLowerCase();
      list.innerHTML = ''; let n = 0;
      for (const id in rows) {
        const r = rows[id], x = RX[id]; if (!x) continue;
        if (v ? r.v !== v : !['wrong', 'ranked', 'refused'].includes(r.v)) continue;
        if (c !== '' && x.cls !== +c) continue; if (s !== '' && x.src !== +s) continue;
        const text = ((r.why || '') + ' ' + (r.label || '') + ' ' + D.classes[x.cls] + ' ' + x.r.map(i => M[i].s).join(' ') + ' ' + M[x.p].s).toLowerCase();
        if (q && !text.includes(q)) continue;
        if (++n > 60) break;
        const m = el('div', 'miss');
        const tagCls = { right: 'buy', ranked: 'made', wrong: 'target', refused: 'dead' }[r.v];
        m.innerHTML = `<span class="tag ${tagCls}">${r.v}</span><b>${esc(D.classes[x.cls])}</b> <span class="hint">· ${esc(D.sources[x.src])} · #${x.id}${x.cite ? ' · ' + esc(x.cite) : ''}</span>` +
          `<div class="mono">${x.r.map(i => esc(M[i].s)).join(' + ')}${x.a.length ? '  [' + x.a.map(i => esc(NAMES.agents[M[i].s] || M[i].f)).join(', ') + ']' : ''} → <b>${esc(M[x.p].s)}</b></div>` +
          (r.v === 'refused' ? `<div class="hint">${esc(r.why || '')}${r.label ? ' — via ' + esc(r.label) : ''}</div>` :
            `<div class="hint">engine (via ${esc(r.label)}): <span class="mono">${esc(r.got || '')}</span></div>`);
        if (r.sub != null && r.via) {
          const b = el('button', 'use', 'Open in Reactions'); b.style.marginTop = '4px';
          b.addEventListener('click', () => { document.querySelector('.tab[data-page="reactions"]').click(); window.Reactions.set(x.r.map(i => M[i].s), r.via.startsWith('mol:') ? '' : String(r.label || '').replace(' [ext]', '')); });
          m.appendChild(b);
        }
        list.appendChild(m);
      }
      if (!n) list.appendChild(el('p', 'hint', 'Nothing matches.'));
      else if (n > 60) list.appendChild(el('p', 'hint', 'Showing the first 60.'));
    };
    ['vfV', 'vfC', 'vfS'].forEach(id => { const e = document.getElementById(id); if (e) e.addEventListener('change', draw); });
    document.getElementById('vfQ').addEventListener('input', draw);
    draw();
  }

  /* ---- wiring ---- */
  document.getElementById('validRun').addEventListener('click', run);
  document.getElementById('validStop').addEventListener('click', stop);
  document.getElementById('validReset').addEventListener('click', () => { if (running) stop(); state = null; store.del(KEY()); bar.style.width = '0'; status.textContent = ''; document.getElementById('validRun').textContent = 'Run the benchmark'; document.getElementById('validRun').disabled = false; render(); });
  document.getElementById('validExport').addEventListener('click', () => {
    if (!state) return;
    const blob = new Blob([JSON.stringify({ set: setId, method: 'see js/validation.js', graded: state.i, of: RX.length, started: state.started, finished: state.finished || null, summary: summarise(state.rows),
      by_source: D.sources.map((s, i) => [s, summarise(state.rows, x => x.src === i)]), by_class: D.classes.map((c, i) => [c, summarise(state.rows, x => x.cls === i)]), rows: state.rows }, null, 1)], { type: 'application/json' });
    const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'benchmark-' + setId + '.json'; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  });
  renderSources();
  chooseSet('db');
  window.Validation = { run, stop, grade, summarise, chooseSet };
})();
