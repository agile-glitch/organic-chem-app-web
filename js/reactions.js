/* Reactions page — reactants (+ reagents) → the predicted product, with the evidence for it.

   Two predictors, asked in this order and shown side by side:
     1. RXNKB (js/rxnkb.js)  the nearest recorded precedents: their transformation applied to your molecules,
        ranked by how similar those precedents are, with a confidence measured on held-out reactions and the
        patent reactions themselves as evidence. Covers whatever recurs in 40 years of patents.
     2. RXN + RXX (js/rxn_engine.js, js/rxn_ext.js)  the mechanistic rule engine: which group, which site, which
        way round, the stereochemistry and the reasons — for the reagents it has been taught.
   When both answer, the page says whether they agree. When neither has grounds, it says so and predicts nothing.
   Chem (js/chem.js) reads the names/SMILES; Reagents (js/reagents.js) reads a written reagent into a rule. */
(() => {
  'use strict';
  const C = window.Chem, M = window.MolDraw, R = window.Reagents, KB = window.RXNKB, EQ = window.Equation;
  const box = document.getElementById('rxnReactants'), B = document.getElementById('rxnB'), out = document.getElementById('rxnOut'), status = document.getElementById('rxnStatus');
  if (!box || !C || typeof RXN === 'undefined' || !R) return;   // RXN is a top-level const, not a window property

  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const el = (tag, cls, html) => { const n = document.createElement(tag); if (cls) n.className = cls; if (html != null) n.innerHTML = html; return n; };
  const nameOf = g => { try { return C.displayNameFor(g) || C.formula(g); } catch (e) { return '?'; } };
  const formulaOf = g => { try { return C.formula(g); } catch (e) { return ''; } };
  const smilesOf = g => { try { return C.toSmiles(g); } catch (e) { return ''; } };
  const pct = x => (x >= 0.995 ? '>99' : x < 0.005 ? '<1' : Math.round(x * 100)) + '%';
  const graphCache = new Map();
  const graphOfSmiles = s => { if (!graphCache.has(s)) { let g = null; try { g = C.parseSmiles(s); } catch (e) {} graphCache.set(s, g); } return graphCache.get(s); };
  const nameOfSmiles = s => { const g = graphOfSmiles(s); return g ? nameOf(g) : ''; };

  function molCard(g, caption, w = 170, h = 100) {
    const card = el('div', 'molcard');
    card.appendChild(M.svg(g, w, h));
    card.appendChild(el('div', 'nm', esc(caption != null ? caption : nameOf(g)) + `<small>${esc(formulaOf(g))}</small>`));
    return card;
  }
  /* a molecule picture from a SMILES: RDKit's depiction when it is loaded, the app's own otherwise */
  function pic(smiles, w, h, captionHtml) {
    const d = el('div', 'rdpic'); let svg = '';
    try { if (KB && KB.ready) { const m = KB.rdkit.get_mol(smiles); if (m) { svg = m.get_svg(w, h).replace(/<\?xml[^>]*\?>/, ''); m.delete(); } } } catch (e) { svg = ''; }
    if (svg) d.innerHTML = svg; else { const g = graphOfSmiles(smiles); if (g) d.appendChild(M.svg(g, w, h)); }
    if (captionHtml != null) d.appendChild(el('div', 'nm', captionHtml));
    return d;
  }

  /* ================= the form ================= */
  const rows = [];
  function addRow(value = '') {
    if (rows.length >= 4) return null;
    const wrap = el('div', 'rxn-reactant', '<span class="import-wrap"><input spellcheck="false" autocomplete="off"><div class="suggest" hidden></div></span><button class="rm" title="Remove this reactant">✕</button>');
    const input = wrap.querySelector('input'); input.value = value;
    const row = { wrap, input };
    window.Sketcher.attachSuggest(input, wrap.querySelector('.suggest'), () => predict());
    wrap.querySelector('.rm').addEventListener('click', () => { if (rows.length > 1) { rows.splice(rows.indexOf(row), 1); wrap.remove(); relabel(); } else input.value = ''; });
    rows.push(row); box.appendChild(wrap); relabel();
    return row;
  }
  function relabel() { rows.forEach((r, i) => { r.input.placeholder = i === 0 ? 'name, SMILES or CAS — e.g. benzaldehyde' : 'another reactant (optional)'; }); document.getElementById('rxnAdd').disabled = rows.length >= 4; }
  function setForm(reactants, reagentText) {
    rows.splice(0).forEach(r => r.wrap.remove());
    (reactants && reactants.length ? reactants : ['']).slice(0, 4).forEach(s => addRow(s));
    if (rows.length < 2) addRow('');
    B.value = reagentText || '';
  }
  function readEnv() {
    const T = document.getElementById('rxnT').value, pH = document.getElementById('rxnPH').value;
    return { T: T === '' ? null : +T + 273.15, pH: pH === '' ? null : +pH, P: 1, solvent: document.getElementById('rxnSolv').value,
      light: document.getElementById('rxnLight').checked, equiv: document.getElementById('rxnEquiv').value === 'excess' ? 'excess' : 1 };
  }

  /* ================= reagents written over the arrow → structures =================
     [label, SMILES, ...aliases]. Salts are written as their ions, the way the patent data records them. */
  const REAGENTS = [
    ['NaBH₄', '[Na+].[BH4-]', 'sodium borohydride'], ['LiAlH₄', '[Li+].[AlH4-]', 'lah', 'lithium aluminium hydride', 'lithium aluminum hydride'], ['NaH', '[H-].[Na+]', 'sodium hydride'],
    ['NaOH', '[Na+].[OH-]', 'sodium hydroxide'], ['KOH', '[K+].[OH-]', 'potassium hydroxide'], ['LiOH', '[Li+].[OH-]', 'lithium hydroxide'],
    ['K₂CO₃', 'O=C([O-])[O-].[K+].[K+]', 'potassium carbonate'], ['Na₂CO₃', 'O=C([O-])[O-].[Na+].[Na+]', 'sodium carbonate'], ['Cs₂CO₃', 'O=C([O-])[O-].[Cs+].[Cs+]', 'cesium carbonate', 'caesium carbonate'],
    ['NaHCO₃', 'O=C([O-])O.[Na+]', 'sodium bicarbonate'], ['Et₃N', 'CCN(CC)CC', 'tea', 'triethylamine', 'net3'], ['DIPEA', 'CCN(C(C)C)C(C)C', 'hunigs base', 'diea', 'ipr2net'],
    ['pyridine', 'c1ccncc1', 'py'], ['DMAP', 'CN(C)c1ccncc1'], ['DBU', 'C1CCC2=NCCCN2CC1'], ['imidazole', 'c1c[nH]cn1'],
    ['HCl', 'Cl', 'hydrochloric acid'], ['HBr', 'Br'], ['HI', 'I'], ['H₂SO₄', 'O=S(=O)(O)O', 'sulfuric acid'], ['HNO₃', 'O=[N+]([O-])O', 'nitric acid'], ['AcOH', 'CC(=O)O', 'acetic acid', 'hoac'],
    ['TFA', 'O=C(O)C(F)(F)F', 'trifluoroacetic acid'], ['TsOH', 'Cc1ccc(S(=O)(=O)O)cc1', 'ptsa', 'ptsoh'], ['H₂', '[H][H]', 'hydrogen'], ['Pd/C', '[Pd]', 'pd', 'pdc', 'palladium on carbon', 'palladium'],
    ['H₂ / Pd', '[H][H].[Pd]', 'h2/pd', 'h2/pdc', 'h2pd', 'h2pdc', 'h2,pd/c'], ['Pt', '[Pt]', 'pto2', 'platinum'], ['Raney Ni', '[Ni]', 'raney nickel', 'ni'],
    ['Pd(PPh₃)₄', 'c1ccc(P(c2ccccc2)c2ccccc2)cc1.[Pd]', 'tetrakis', 'pdpph34'], ['PPh₃', 'c1ccc(P(c2ccccc2)c2ccccc2)cc1', 'triphenylphosphine'], ['Pd(OAc)₂', 'CC(=O)[O-].CC(=O)[O-].[Pd+2]', 'pdoac2', 'palladium acetate'],
    ['CuI', '[Cu]I', 'copper iodide'], ['SOCl₂', 'O=S(Cl)Cl', 'thionyl chloride'], ['(COCl)₂', 'O=C(Cl)C(=O)Cl', 'oxalyl chloride', 'cocl2'], ['POCl₃', 'O=P(Cl)(Cl)Cl'], ['PBr₃', 'BrP(Br)Br'], ['PCl₅', 'ClP(Cl)(Cl)(Cl)Cl'],
    ['mCPBA', 'O=C(OO)c1cccc(Cl)c1', 'm-cpba'], ['H₂O₂', 'OO', 'hydrogen peroxide'], ['KMnO₄', 'O=[Mn](=O)(=O)[O-].[K+]', 'potassium permanganate'], ['CrO₃', 'O=[Cr](=O)=O', 'jones'],
    ['PCC', 'O=[Cr](=O)([O-])Cl.c1cc[nH+]cc1'], ['DMP', 'CC(=O)OI1(OC(C)=O)(OC(C)=O)OC(=O)c2ccccc21', 'dess-martin', 'dess martin periodinane'], ['MnO₂', 'O=[Mn]=O'],
    ['DCC', 'C(=NC1CCCCC1)=NC1CCCCC1'], ['EDC', 'CCN=C=NCCCN(C)C', 'edci', 'edac'], ['HOBt', 'On1nnc2ccccc21'], ['HATU', 'CN(C)C(On1nnc2cccnc21)=[N+](C)C.F[P-](F)(F)(F)(F)F'],
    ['NBS', 'O=C1CCC(=O)N1Br', 'n-bromosuccinimide'], ['NCS', 'O=C1CCC(=O)N1Cl'], ['AIBN', 'CC(C)(C#N)N=NC(C)(C)C#N'], ['Br₂', 'BrBr', 'bromine'], ['Cl₂', 'ClCl', 'chlorine'], ['I₂', 'II', 'iodine'],
    ['NaCN', '[C-]#N.[Na+]', 'sodium cyanide'], ['KCN', '[C-]#N.[K+]'], ['NaN₃', '[N-]=[N+]=[N-].[Na+]', 'sodium azide'], ['NH₃', 'N', 'ammonia'], ['NH₄Cl', '[Cl-].[NH4+]', 'ammonium chloride'],
    ['NH₂OH', 'NO', 'hydroxylamine'], ['N₂H₄', 'NN', 'hydrazine'], ['LDA', 'CC(C)[N-]C(C)C.[Li+]'], ['n-BuLi', '[Li]CCCC', 'buli', 'nbuli', 'butyllithium'], ['tBuOK', 'CC(C)(C)[O-].[K+]', 'kotbu', 'potassium tert-butoxide'],
    ['NaOEt', 'CC[O-].[Na+]', 'sodium ethoxide'], ['NaOMe', 'C[O-].[Na+]', 'sodium methoxide'], ['MeI', 'CI', 'iodomethane', 'methyl iodide'], ['BnBr', 'BrCc1ccccc1', 'benzyl bromide'],
    ['Boc₂O', 'CC(C)(C)OC(=O)OC(=O)OC(C)(C)C'], ['TBSCl', 'CC(C)(C)[Si](C)(C)Cl', 'tbdmscl'], ['TBAF', 'CCCC[N+](CCCC)(CCCC)CCCC.[F-]'], ['Ac₂O', 'CC(=O)OC(C)=O', 'acetic anhydride'],
    ['TsCl', 'Cc1ccc(S(=O)(=O)Cl)cc1', 'tosyl chloride'], ['MsCl', 'CS(=O)(=O)Cl', 'mesyl chloride'], ['AlCl₃', 'Cl[Al](Cl)Cl', 'aluminium chloride', 'aluminum chloride'], ['BF₃', 'FB(F)F'],
    ['Zn', '[Zn]', 'zinc'], ['Fe', '[Fe]', 'iron'], ['SnCl₂', 'Cl[Sn]Cl'], ['Mg', '[Mg]', 'magnesium'], ['NaBH₃CN', '[BH3-]C#N.[Na+]', 'sodium cyanoborohydride'], ['NaBH(OAc)₃', 'CC(=O)O[BH-](OC(C)=O)OC(C)=O.[Na+]', 'stab', 'sodium triacetoxyborohydride'],
    /* solvents */
    ['THF', 'C1CCOC1', 'tetrahydrofuran'], ['CH₂Cl₂', 'ClCCl', 'dcm', 'dichloromethane', 'methylene chloride'], ['DMF', 'CN(C)C=O'], ['DMSO', 'CS(C)=O'], ['MeOH', 'CO', 'methanol'], ['EtOH', 'CCO', 'ethanol'],
    ['water', 'O', 'h2o'], ['toluene', 'Cc1ccccc1'], ['acetone', 'CC(C)=O'], ['MeCN', 'CC#N', 'acetonitrile'], ['EtOAc', 'CCOC(C)=O', 'ethyl acetate'], ['Et₂O', 'CCOCC', 'ether', 'diethyl ether'],
    ['dioxane', 'C1COCCO1', '1,4-dioxane'], ['CHCl₃', 'ClC(Cl)Cl', 'chloroform'], ['benzene', 'c1ccccc1'], ['hexane', 'CCCCCC'], ['iPrOH', 'CC(C)O', 'ipa', 'isopropanol'], ['DME', 'COCCOC'],
    /* bare ions, for naming what the precedents record */
    ['Na⁺', '[Na+]'], ['K⁺', '[K+]'], ['Li⁺', '[Li+]'], ['Cs⁺', '[Cs+]'], ['Cl⁻', '[Cl-]'], ['Br⁻', '[Br-]'], ['I⁻', '[I-]'], ['OH⁻', '[OH-]'], ['CO₃²⁻', 'O=C([O-])[O-]'], ['HCO₃⁻', 'O=C([O-])O'],
    ['BH₄⁻', '[BH4-]'], ['AlH₄⁻', '[AlH4-]'], ['H⁻', '[H-]'], ['NH₄⁺', '[NH4+]'], ['AcO⁻', 'CC(=O)[O-]'], ['SO₄²⁻', 'O=S(=O)([O-])[O-]'], ['CN⁻', '[C-]#N'], ['N₃⁻', '[N-]=[N+]=[N-]'], ['F⁻', '[F-]'],
  ];
  window.REAGENT_LABELS = REAGENTS;                                 // [label, SMILES, …aliases]; also used by js/retrochimera.js
  const normKey = s => String(s).toLowerCase().replace(/[₀-₉]/g, d => '₀₁₂₃₄₅₆₇₈₉'.indexOf(d)).replace(/[⁺]/g, '+').replace(/[⁻]/g, '-').replace(/[\s.\-()\[\]]/g, '');
  const REAGENT_BY_KEY = new Map();
  REAGENTS.forEach(([label, smi, ...aliases]) => [label, ...aliases].forEach(a => { const k = normKey(a); if (!REAGENT_BY_KEY.has(k)) REAGENT_BY_KEY.set(k, { label, smiles: smi }); }));
  function agentsFromText(text) {
    const tokens = [];
    String(text || '').split(/[,;]|\s+then\s+|\s+and\s+/i).map(s => s.trim()).filter(Boolean).forEach(tok => {
      const tryOne = t => { const hit = REAGENT_BY_KEY.get(normKey(t)); if (hit) return hit; try { const g = C.searchMolecule(t).graph; const s = smilesOf(g); if (s) return { label: nameOf(g), smiles: s }; } catch (e) {} return null; };
      const whole = tryOne(tok);
      if (whole) { tokens.push({ text: tok, ...whole }); return; }
      const parts = tok.split('/').map(s => s.trim()).filter(Boolean);
      if (parts.length > 1) parts.forEach(p => { const h = tryOne(p); tokens.push(h ? { text: p, ...h } : { text: p, smiles: null }); });
      else tokens.push({ text: tok, smiles: null });
    });
    return { tokens, smiles: tokens.filter(t => t.smiles).map(t => t.smiles) };
  }
  let AGENT_LABEL = null;                                           // canonical SMILES → label, for the precedents' recorded agents
  function agentLabel(smiles) {
    if (!AGENT_LABEL && KB && KB.ready) { AGENT_LABEL = new Map(); REAGENTS.forEach(([label, smi]) => { if (smi.indexOf('.') < 0) { const c = KB.canonical(smi); if (c && !AGENT_LABEL.has(c.iso)) AGENT_LABEL.set(c.iso, label); } }); }
    return (AGENT_LABEL && AGENT_LABEL.get(smiles)) || null;
  }

  /* ================= reading the reactants ================= */
  const SPECTATOR_ION = /^\[(Li|Na|K|Cs|Rb|Mg|Ca|Zn|Cl|Br|I|F)[+-]\d?\]$/;
  async function readReactants() {
    const list = [], spectators = [];
    for (const r of rows) {
      const text = r.input.value.trim(); if (!text) continue;
      let g; try { g = C.searchMolecule(text).graph; } catch (e) {
        // not in the curated library: saved PubChem lookups, then the broad offline library (js/pclib.js)
        const lookupable = !/[=#\[\]@]/.test(text) && (/[a-z]{3}/i.test(text) || /^\d{2,7}-\d\d-\d$/.test(text));
        const res = lookupable && window.PCLib ? await window.PCLib.resolve(text) : null;
        if (res && res.ok) { try { g = C.parseSmiles(res.smiles); } catch (e2) { g = null; } }
        if (!g) {
          if (res && res.choices) return { error: `"${text}": ${res.error} ` + res.choices.map(c => `${c.title} = ${c.smiles}`).join(' · ') + ' (paste the SMILES you mean).' };
          return { error: `Could not read the reactant "${text}": ${e.message}` };
        }
      }
      let comps; try { comps = C.componentsOf(g); } catch (e) { comps = [g.atoms.map(a => a.id)]; }
      for (const ids of comps) {
        const set = new Set(ids);
        const sub = comps.length === 1 ? g : { atoms: g.atoms.filter(a => set.has(a.id)).map(a => ({ ...a })), bonds: g.bonds.filter(b => set.has(b.a) && set.has(b.b)).map(b => ({ ...b })), nextId: g.nextId };
        const smiles = smilesOf(sub); if (!smiles) continue;
        if (SPECTATOR_ION.test(smiles)) { spectators.push(smiles); continue; }          // a counter-ion is a reagent, not a reactant
        list.push({ text, g: sub, smiles, name: nameOf(sub) });
      }
    }
    return { list: list.slice(0, 4), spectators };
  }

  /* ================= the mechanistic model ================= */
  function mechanistic(reactants, reagentText, flatOf) {
    const env = readEnv(), tries = [];
    const graphs = reactants.map(r => r.g).slice().sort((a, b) => b.atoms.length - a.atoms.length);
    const ids = [];
    if (reagentText) {
      const texts = [reagentText].concat(reagentText.split(/[,;]|\s+then\s+/i).map(s => s.trim()).filter(Boolean));
      for (const t of texts) {
        window.REAGENT_SUBSTRATE = graphs[0];
        let r = null; try { r = R.resolve(t); } catch (e) { r = null; }
        if (r && r.id && !String(r.id).startsWith('mol:') && !ids.some(x => x.id === r.id)) ids.push({ id: r.id, via: r.via || 'name', text: t });
      }
    }
    for (const g of graphs) for (const rid of ids) {
      let res = null; try { res = window.RXX ? RXX.runReagent(g, rid.id) : null; } catch (e) { res = null; }
      if (!res) { try { res = RXN.run(g, rid.id, env); } catch (e) { res = null; } }
      if (res) tries.push({ sub: g, res, rr: { id: rid.id, via: rid.via + (res.ext ? '; the added family "' + res.rule.label + '"' : ''), text: rid.text } });
    }
    if (graphs.length >= 2) for (const a of graphs) for (const b of graphs) {
      if (a === b) continue;
      let res = null; try { res = window.RXX ? RXX.run(a, b) : null; } catch (e) { res = null; }
      if (!res) { try { const rule = R.moleculeRule(b, nameOf(b)); res = rule ? RXN.run(a, rule.id, env) : null; } catch (e) { res = null; } }
      if (res) tries.push({ sub: a, res, rr: { id: res.rule.id, via: 'the two molecules as partners' } });
    }
    /* `first` is the model's committed answer (the first route, in this order, that gives a product): the one the
       arbitration uses, as in the benchmark. The other routes are kept only to explain an agreeing product. */
    const good = tries.filter(t => !t.res.error && t.res.major);
    good.forEach(t => { try { t.flat = flatOf(t.res.major.smiles); } catch (e) { t.flat = null; } });
    return { first: good[0] || null, good, fail: good.length ? null : (tries[0] || null) };
  }

  /* ================= combining the two =================
     Measured with the same logic (Validation tab, "combined") on every 10th reaction of the built-in set — 1,434
     graded, none of whose products the precedent knowledge base learned from:
       agree  the mechanistic answer equals the precedents' first         right 235 of 304 (77 %)
       pick   it is 2nd or 3rd among the precedent candidates             its pick right 39 of 93 (42 %); the precedents' first 25 of 93 (27 %)
       alone  the precedents' first, with no mechanistic support          right 395 of 1,015 (39 %)
       closest precedent under 30 % similar                               right 0 of 8 — so nothing is asserted there */
  const MEASURED = { agree: [235, 304], pick: [39, 93], pickFirst: [25, 93] };
  function decide(cands, mech) {
    const top = cands[0] || null, conf = top && top.confidence.p != null ? top.confidence.p : null;
    const reliable = !!top && conf != null && conf >= RELIABLE;
    const rf = mech && mech.first ? mech.first : null;
    if (reliable && rf && rf.flat === top.flat) return { how: 'agree', flat: top.flat, smiles: top.smiles, cand: top, route: rf };
    const k = reliable && rf && rf.flat ? cands.slice(0, 3).findIndex(c => c.flat === rf.flat) : -1;
    if (k > 0) return { how: 'pick', flat: rf.flat, smiles: cands[k].smiles, cand: cands[k], rank: k + 1, route: rf };
    if (reliable) return { how: 'precedent', flat: top.flat, smiles: top.smiles, cand: top, route: (mech && mech.good.find(t => t.flat === top.flat)) || null };
    if (rf) return { how: 'rules', flat: rf.flat, smiles: rf.res.major.smiles, route: rf };
    return { how: 'none' };
  }

  /* ================= predict ================= */
  let runId = 0;
  async function predict() {
    const my = ++runId;
    const rd = await readReactants();
    if (rd.error) { out.innerHTML = `<div class="banner">${esc(rd.error)}</div>`; return; }
    if (!rd.list.length) { out.innerHTML = '<div class="banner">Enter at least one reactant.</div>'; return; }
    const reagentText = B.value.trim(), agents = agentsFromText(reagentText);
    out.innerHTML = '<p class="hint">Working…</p>';
    let kb = null, kbNote = '';
    if (KB) {
      try {
        await KB.load((msg, frac) => { if (my !== runId) return; status.innerHTML = esc(msg) + '<div class="progress"><div style="width:' + Math.round(100 * frac) + '%"></div></div>'; });
        if (my !== runId) return;
        status.textContent = '';
        kb = KB.predict({ reactants: rd.list.map(r => r.smiles), agents: agents.smiles.concat(rd.spectators) });
        if (kb.error) { kbNote = kb.error; kb = null; }
      } catch (e) { status.textContent = ''; kbNote = 'The precedent knowledge base is not installed (data/rxnkb/, built by tools/rxnkb.py), so only the mechanistic model is used.'; }
    }
    if (my !== runId) return;
    const flatOf = s => { const c = KB && KB.ready ? KB.canonical(s) : null; return c ? c.flat : null; };
    const mech = mechanistic(rd.list, reagentText, flatOf);
    render(my, rd, agents, kb, kbNote, mech);
  }

  const band = p => (p == null ? ['none', 'not rated'] : p >= 0.8 ? ['high', 'High confidence'] : p >= 0.5 ? ['medium', 'Medium confidence'] : p >= 0.25 ? ['low', 'Low confidence'] : ['low', 'Very low confidence']);
  const confBadge = c => band(c ? c.p : null);
  const RELIABLE = 0.2;                                             // below this measured hit rate — or with no measurement — the page does not rely on precedents

  function render(my, rd, agents, kb, kbNote, mech) {
    out.innerHTML = '';
    const reactantSmiles = rd.list.map(r => r.smiles);
    const cands = kb ? kb.candidates : [];
    const top = cands[0] || null;
    const dec = decide(cands, mech);
    window.Reactions.last = { how: dec.how, smiles: dec.smiles || null, p: top && top.confidence ? top.confidence.p : null,   // what was decided, for tools that grade the page
      rules: mech && mech.first ? mech.first.res.major.smiles : null, precedent: top ? top.smiles : null };
    const topReliable = dec.how !== 'rules' && dec.how !== 'none';
    const unknownAgents = agents.tokens.filter(t => !t.smiles).map(t => t.text);

    /* ---- 1. the answer ---- */
    const card = el('div', 'card');
    const answer = el('div', 'answer');
    const eqn = el('div', 'eqn');
    rd.list.forEach((r, i) => { if (i) eqn.appendChild(el('span', 'plus', '+')); eqn.appendChild(pic(r.smiles, 180, 110, esc(r.name) + `<small>${esc(formulaOf(r.g))}</small>`)); });
    const arrow = el('div', 'arrow');
    arrow.appendChild(el('div', null, esc(agents.tokens.map(t => t.label || t.text).join(', '))));
    arrow.appendChild(el('div', 'ar', '⟶'));
    eqn.appendChild(arrow);
    /* the headline: drawn by the rule engine when it made this product (it draws stereochemistry), else by RDKit */
    let headline = null;
    if (dec.how !== 'none') {
      const r = dec.route && dec.route.res && dec.route.flat === dec.flat ? dec.route.res.major : null;
      headline = { smiles: dec.smiles, name: (r && r.name) || nameOfSmiles(dec.smiles), graph: r ? r.graph : null, stereo: r ? r.stereo : null };
      if (headline.stereo && headline.stereo.forms && headline.stereo.forms.length > 1)
        headline.stereo.forms.forEach((f, i) => { if (i) eqn.appendChild(el('span', 'plus', '+')); eqn.appendChild(molCard(f.graph, headline.name + ' ' + f.label, 180, 110)); });
      else eqn.appendChild(headline.graph ? molCard(headline.graph, headline.name, 190, 110) : pic(headline.smiles, 200, 120, `<b>${esc(headline.name || '')}</b><small>${esc(formulaOf(graphOfSmiles(headline.smiles) || { atoms: [] }))}</small>`));
    } else eqn.appendChild(el('div', 'rdpic', '<div class="nm" style="color:#b3261e;width:190px">no prediction</div>'));
    answer.appendChild(eqn); card.appendChild(answer);

    const basis = el('div', 'basis'), rate = ([a, b]) => `${pct(a / b)} (${a} of ${b})`;
    const backing = c => `Backed by <b>${c.support.length}</b> of the ${kb.neighbours.length} most similar recorded reactions (closest ${pct(c.best)} similar)`;
    if (dec.how === 'agree') {
      const pA = MEASURED.agree[0] / MEASURED.agree[1], pP = top.confidence.p, [cls, label] = band(Math.max(pA, pP || 0));
      basis.innerHTML = `<span class="conf ${cls}">${label}</span> &nbsp;<b>Two independent methods agree.</b> ${backing(top)}, and the mechanistic model gives the same product. ` +
        `On reactions the program had never seen, answers where the two agreed were right ${rate(MEASURED.agree)}; first answers with this much precedent alone were right ${pct(pP)} (n = ${top.confidence.n}).`;
    } else if (dec.how === 'pick') {
      const [cls, label] = band(MEASURED.pick[0] / MEASURED.pick[1]);
      basis.innerHTML = `<span class="conf ${cls}">${label}</span> &nbsp;The precedents rank this product <b>no. ${dec.rank}</b>; the mechanistic model gives it as its answer, so it is put first. ` +
        `When that happened on reactions the program had never seen, the mechanistic pick was right ${rate(MEASURED.pick)}, against ${rate(MEASURED.pickFirst)} for the precedents' own first choice. ${backing(dec.cand)}.`;
    } else if (dec.how === 'precedent') {
      const [cls, label] = confBadge(top.confidence);
      basis.innerHTML = `<span class="conf ${cls}">${label}</span> &nbsp;${backing(top)}. On reactions the program had never seen, first answers with this much precedent behind them were right <b>${pct(top.confidence.p)}</b> of the time (n = ${top.confidence.n}).` +
        (mech && mech.first ? ' The mechanistic model gives a different product, not among the top three candidates; see below.' : '');
    } else if (dec.how === 'rules') {
      basis.innerHTML = `<span class="conf none">rule-based, unconfirmed</span> &nbsp;` + (top ? `The recorded precedents are too distant to rely on (closest ${pct(top.best)} similar${top.confidence.p != null ? '; candidates like this were right ' + pct(top.confidence.p) + ' of the time' : ', a range where testing found no reliable answers'}), so this is the mechanistic model's answer on its own. It has no recorded precedent behind it.` : 'No recorded precedent applies to these molecules, so this is the mechanistic model\'s answer on its own. It has no recorded precedent behind it.');
    } else {
      basis.innerHTML = '<span class="conf none">no prediction</span> &nbsp;' + (kb ? (top ? `The nearest recorded reactions are too unlike these molecules to assert a product (closest ${pct(top.best)} similar${top.confidence.p != null ? '; in testing, candidates this weak were right ' + pct(top.confidence.p) + ' of the time' : ''}). They are listed below for orientation only.` : 'None of the most similar recorded reactions has a transformation that applies to these molecules.') : '') +
        (mech && mech.fail && mech.fail.res.error ? ' The mechanistic model: ' + esc(mech.fail.res.error) : (!B.value.trim() && rd.list.length < 2 ? ' Add a reagent or a second reactant.' : ''));
    }
    card.appendChild(basis);
    const notes = [];
    if (kb) notes.push(`Compared with ${KB.meta.n.toLocaleString()} recorded reactions in ${Math.round(kb.ms.total)} ms` + (kb.usedAgents ? ', taking the reagents into account.' : '.'));
    if (kbNote) notes.push(kbNote);
    if (unknownAgents.length) notes.push('Not recognised, so ignored by the precedent search: ' + unknownAgents.join(', ') + '.');
    if (notes.length) card.appendChild(el('div', 'hint', esc(notes.join(' '))));
    if (headline) {
      const use = el('button', 'use', 'Use the product as the next starting material →'); use.style.marginTop = '8px';
      use.addEventListener('click', () => { setForm([headline.smiles], ''); window.scrollTo(0, 0); rows[1] ? rows[1].input.focus() : B.focus(); });
      card.appendChild(use);
    }
    out.appendChild(card);

    /* ---- 2. the candidates and their evidence ---- */
    if (cands.length) {
      out.appendChild(el('h2', null, (topReliable ? 'Candidates' : 'Weak candidates, for orientation only') + ` <small>${Math.min(5, cands.length)} of ${cands.length} · ranked by the similar precedent behind each</small>`));
      const strip = el('div', 'cands'), evBox = el('div');
      const firstShown = Math.max(0, cands.findIndex(c => c.flat === dec.flat));
      cands.slice(0, 5).forEach((c, k) => {
        const [cls, label] = confBadge(c.confidence);
        const cd = el('div', 'cand' + (k === firstShown ? ' on' : ''));
        cd.appendChild(pic(c.smiles, 226, 110));
        cd.appendChild(el('div', null, `<b>${esc(nameOfSmiles(c.smiles) || c.smiles)}</b>`));
        cd.appendChild(el('div', 'bar', `<div style="width:${Math.round(100 * c.share)}%"></div>`));
        cd.appendChild(el('div', 'hint', `${pct(c.share)} of the support · ${c.support.length} precedent${c.support.length === 1 ? '' : 's'} · closest ${pct(c.best)} <span class="conf ${cls}" style="margin-left:4px">${label.replace(' confidence', '')}</span>`));
        cd.addEventListener('click', () => { strip.querySelectorAll('.cand').forEach(x => x.classList.remove('on')); cd.classList.add('on'); evidence(evBox, c, my); });
        strip.appendChild(cd);
      });
      out.appendChild(strip);
      out.appendChild(el('h2', null, 'The evidence <small>the recorded reactions this candidate comes from</small>'));
      out.appendChild(evBox);
      evidence(evBox, cands[firstShown], my);
    }

    /* ---- 3. the mechanistic model: the route behind the answer if it has one, else its own first answer ---- */
    const shown = (dec.route && dec.route.flat === dec.flat ? dec.route : null) || (mech && (mech.first || mech.fail)) || null;
    if (shown) {
      out.appendChild(el('h2', null, 'Mechanistic model <small>which group, which site, which way round — for the reagents it has been taught</small>'));
      const note = dec.how === 'agree' ? ['yes', '✓ Its answer is the same as the precedent-based one: two independent methods give the same product.']
        : dec.how === 'pick' ? ['yes', `Its answer is candidate no. ${dec.rank} of the precedent list, and is shown first for that reason.`]
        : dec.how === 'precedent' && shown.res.major && shown.flat !== dec.flat ? ['no', 'Its answer differs from the precedent-based one and is not among the top three candidates. Check the evidence before relying on either.']
        : dec.how === 'rules' ? ['no', 'No close recorded precedent supports this answer; it rests on the rules alone.'] : null;
      if (note) out.appendChild(el('div', 'agree ' + note[0], note[1]));
      mechanisticView(shown);
    }

    /* ---- 4. balanced equation, and whether the product has been made ---- */
    if (headline) {
      out.appendChild(el('h2', null, 'Balanced equation <small>for the answer above</small>'));
      let eq = null;
      if (dec.route && dec.route.flat === dec.flat && dec.route.res.major) { try { eq = EQ.forRule(dec.route.sub, dec.route.res); } catch (e) { eq = null; } }
      if (!eq) eq = EQ.forSmiles(reactantSmiles, headline.smiles, headline.name);
      out.appendChild(eq || el('p', 'hint', 'The equation could not be written for this product.'));
      const pg = headline.graph || graphOfSmiles(headline.smiles);
      if (pg) precedents(rd.list[0].g, pg, headline.name || headline.smiles);
      // papers naming a reactant and the product (PubChem; on a click; not reaction-level evidence, and said so)
      if (window.LitRef) {
        out.appendChild(el('h2', null, 'Mentioned together in the literature <small>papers in PubChem that name a reactant and the product — not necessarily this reaction</small>'));
        const box = window.LitRef.panel();
        out.appendChild(window.LitRef.button('Look up papers (PubChem)', () => rd.list.map(r => ({ a: { smiles: r.smiles, name: r.name }, b: { smiles: headline.smiles, name: headline.name } })), box));
        out.appendChild(box);
      }
    }
    whatElse(rd.list[0].g);
  }

  /* the recorded reactions behind one candidate: the closest few in full, and what was used in all of them */
  async function evidence(container, cand, my) {
    container.innerHTML = '<p class="hint">Reading the supporting reactions…</p>';
    const sup = cand.support.slice().sort((a, b) => b.sim - a.sim);
    let recs = [];
    try { recs = await KB.evidence(sup.map(s => s.i)); } catch (e) { container.innerHTML = `<p class="hint">The precedent details could not be read: ${esc(e.message)}</p>`; return; }
    if (my !== runId) return;
    container.innerHTML = '';
    const counts = new Map(), yields = [], years = [];
    recs.forEach(r => { new Set(r.agents).forEach(a => counts.set(a, (counts.get(a) || 0) + 1)); if (r.yield != null) yields.push(r.yield); if (+r.year) years.push(+r.year); });
    const topAgents = [...counts.entries()].sort((a, b) => b[1] - a[1]).slice(0, 10);
    yields.sort((a, b) => a - b);
    const stats = el('div', 'condstats');
    stats.innerHTML = `<b>What was used</b> in these ${recs.length} reactions: ` + (topAgents.length ? topAgents.map(([s, n]) => `<span class="chipc" title="${esc(s)}">${esc(agentLabel(s) || s)} <small>${n}</small></span>`).join('') : 'no reagents recorded') +
      (yields.length ? `<br><b>Recorded yields:</b> median ${Math.round(yields[Math.floor(yields.length / 2)])} % (range ${Math.round(yields[0])}–${Math.round(yields[yields.length - 1])} %, n = ${yields.length})` : '') +
      (years.length ? ` · patents ${Math.min(...years)}–${Math.max(...years)}` : '');
    container.appendChild(stats);
    recs.slice(0, 6).forEach((r, k) => {
      const ev = el('div', 'evid'), eq = el('div', 'eqn');
      r.reactants.split('.').forEach((s, i) => { if (i) eq.appendChild(el('span', 'plus', '+')); eq.appendChild(pic(s, 130, 80)); });
      eq.appendChild(el('span', 'ar', '⟶')); eq.appendChild(pic(r.product, 150, 80));
      ev.appendChild(eq);
      const cite = window.Precedent && window.Precedent.citeHtml ? window.Precedent.citeHtml(r.cite) : esc(r.cite);
      ev.appendChild(el('div', 'meta', `<span class="simtag">${pct(sup[k].sim)} similar</span> · ` + (r.agents.length ? '<b>with</b> ' + r.agents.slice(0, 8).map(s => esc(agentLabel(s) || s)).join(', ') + ' · ' : '') +
        (r.yield != null ? `<b>yield</b> ${r.yield} % · ` : '') + cite + (r.year ? ' (' + esc(r.year) + ')' : '')));
      container.appendChild(ev);
    });
    if (recs.length > 6) container.appendChild(el('p', 'hint', `Showing the 6 closest of ${recs.length}. Source: ${esc(KB.meta.source)}.`));
    else container.appendChild(el('p', 'hint', 'Source: ' + esc(KB.meta.source) + '.'));
  }

  /* the rule engine's own products, shares, reasons and candidate sites */
  function mechanisticView(mech) {
    const { sub, res, rr } = mech;
    if (res.error) {
      out.appendChild(el('div', 'banner', '<b>No reaction from the mechanistic model.</b> ' + esc(res.error)));
      if (res.ignored && res.ignored.length) out.appendChild(el('p', 'hint', 'Groups present: ' + res.ignored.map(i => esc(i.type) + (i.why ? ' (' + esc(i.why) + ')' : '')).join('; ')));
      return;
    }
    out.appendChild(el('p', 'hint', `${esc(nameOf(sub))} with ${esc(res.rule.label)}` + (rr.via && rr.via !== 'name' ? ` (read via ${esc(rr.via)})` : '') + (R.PRED_COND[res.rule.id] ? ' · ' + esc(R.PRED_COND[res.rule.id]) : '') + '.'));
    const list = el('div', 'products');
    res.outcomes.forEach(o => {
      const p = el('div', 'product ' + o.tag);
      p.appendChild(el('div', 'head', `<span class="share">${pct(o.share)}</span><span class="tag ${o.tag === 'major' ? 'buy' : o.tag === 'minor' ? 'made' : 'dead'}">${o.tag}</span>`));
      p.appendChild(el('div', null, `<b>${esc(o.name)}</b> <span class="hint">${esc(o.formula)} · <code>${esc(o.smiles)}</code></span>`));
      const forms = el('div', 'forms');
      if (o.stereo && o.stereo.forms && o.stereo.forms.length) o.stereo.forms.forEach(f => forms.appendChild(molCard(f.graph, f.label + (o.stereo.forms.length > 1 ? ' · ' + pct(f.share) + ' of this product' : ''), 160, 100)));
      else forms.appendChild(molCard(o.graph, o.label || '', 160, 100));
      p.appendChild(forms);
      if (o.why || o.label) p.appendChild(el('div', 'why', (o.label ? `<b>${esc(o.label)}.</b> ` : '') + esc(o.why || '')));
      if (o.stereo) p.appendChild(el('div', 'stereo', `<b>${esc(o.stereo.kind)}</b> — ${esc(o.stereo.note)}`));
      else if (o.stereoHow) p.appendChild(el('div', 'stereo', esc(o.stereoHow)));
      const use = el('button', 'use', 'Use as the next starting material →');
      use.addEventListener('click', () => { setForm([o.smiles], ''); window.scrollTo(0, 0); });
      p.appendChild(use);
      list.appendChild(p);
    });
    out.appendChild(list);
    if (res.byproducts && res.byproducts.length) out.appendChild(el('div', 'byprod', 'By-products: ' + res.byproducts.map(b => esc(b.name || b.smiles || b)).join(', ')));
    const S = res.summary || {};
    const why = el('div', 'card whybox'); why.style.marginTop = '10px';
    if (S.chemo) why.appendChild(el('p', null, '<b>Which group:</b> ' + esc(S.chemo)));
    if (S.site) why.appendChild(el('p', null, '<b>Which site:</b> ' + esc(S.site)));
    if (S.regio) why.appendChild(el('p', null, '<b>Which way round:</b> ' + esc(S.regio)));
    if (S.stereo) why.appendChild(el('p', null, '<b>Stereochemistry:</b> ' + esc(S.stereo)));
    if (S.conditions) why.appendChild(el('p', null, '<b>Conditions:</b> ' + esc(S.conditions)));
    out.appendChild(why);
    if (res.sites && res.sites.length) {
      const d = el('details'); d.appendChild(el('summary', null, 'Every candidate site (' + res.sites.length + ')'));
      const t = el('table', 'sites');
      t.innerHTML = '<thead><tr><th>Site</th><th>Share</th><th>Reason</th></tr></thead><tbody>' +
        res.sites.map(s => `<tr class="${s.inert ? 'dead' : ''}"><td>${esc(s.label)}</td><td>${s.inert ? '—' : pct(s.share || 0)}</td><td>${esc(s.inert ? 'inert: ' + s.inert : s.reason || '')}</td></tr>`).join('') + '</tbody>';
      d.appendChild(t); out.appendChild(d);
    }
  }

  /* ---- has this product been made? exact-structure search over the recorded-reaction sets ---- */
  function precedents(sub, product, productName) {
    const P = window.Precedent; if (!P) return;
    const made = P.exact(product), used = P.from(sub), skipped = P.skipped(product);
    out.appendChild(el('h2', null, 'Has it been done? <small>this exact product in the recorded-reaction sets</small>'));
    const bx = el('div');
    if (skipped) bx.appendChild(el('p', 'hint', `${skipped} record${skipped === 1 ? '' : 's'} naming this molecule as the product ${skipped === 1 ? 'was' : 'were'} skipped: the recorded product is a solvent or a fragment of a much larger reactant, so it is not a synthesis of it.`));
    if (made.length) {
      bx.appendChild(el('p', null, `<b>${made.length}</b> recorded reaction${made.length === 1 ? '' : 's'} make${made.length === 1 ? 's' : ''} <b>${esc(productName)}</b>` + (made.length > 5 ? ' — the first 5:' : ':')));
      made.slice(0, 5).forEach(x => bx.appendChild(P.card(x)));
    } else {
      bx.appendChild(el('p', null, `No recorded reaction in these sets makes <b>${esc(productName)}</b>.`));
      const sim = P.similar(product, 4);
      if (sim.length) { bx.appendChild(el('p', 'hint', 'Closest recorded products (structural similarity):')); sim.forEach(s => { const c = P.card(s.rows[0]); c.prepend(el('div', 'prec-sim', Math.round(s.sim * 100) + '% similar')); bx.appendChild(c); }); }
    }
    if (used.length) {
      const d = el('details'); d.appendChild(el('summary', null, `${used.length} recorded reaction${used.length === 1 ? '' : 's'} start${used.length === 1 ? 's' : ''} from ${esc(nameOf(sub))}`));
      used.slice(0, 5).forEach(x => d.appendChild(P.card(x)));
      if (used.length > 5) d.appendChild(el('p', 'hint', 'Showing the first 5.'));
      bx.appendChild(d);
    }
    bx.appendChild(el('div', 'prec-note', P.NOTE));
    out.appendChild(bx);
  }

  /* what else this starting material would do: every taught reagent with a live site, on demand */
  function whatElse(sub) {
    const d = el('details'); d.appendChild(el('summary', null, `What else would ${esc(nameOf(sub))} react with? (mechanistic model)`));
    const body = el('div', 'hint', 'Click to find out — every taught reagent is tried.');
    const btn = el('button', 'use', 'Try every reagent'); btn.style.marginTop = '6px';
    btn.addEventListener('click', () => {
      body.textContent = 'Working…';
      setTimeout(() => {
        let hits = [];
        try { hits = RXN.reactiveReagents(sub, {}).filter(h => !String(h.rule.id).startsWith('mol:')); } catch (e) {}
        body.className = 'reactives'; body.innerHTML = '';
        if (!hits.length) { body.className = 'hint'; body.textContent = 'None of the taught reagents has a site on this molecule.'; return; }
        hits.forEach(h => { const b = el('button', null, esc(h.rule.label) + ' → ' + esc(h.res.major.name)); b.addEventListener('click', () => { setForm([smilesOf(sub)], h.rule.label); predict(); window.scrollTo(0, 0); }); body.appendChild(b); });
      }, 20);
      btn.remove();
    });
    d.appendChild(body); d.appendChild(btn); out.appendChild(d);
  }

  /* ================= wiring ================= */
  let REAG = null;
  const reagentIndex = () => {
    if (!REAG) {
      const seen = new Set(); REAG = [];
      REAGENTS.forEach(([label]) => { const k = label.toLowerCase(); if (!/[⁺⁻]$/.test(label) && !seen.has(k)) { seen.add(k); REAG.push({ label, formula: 'reagent', k: normKey(label) }); } });
      C.RULES.forEach(r => { if (String(r.id).startsWith('mol:')) return; const k = r.label.toLowerCase(); if (!seen.has(k)) { seen.add(k); REAG.push({ label: r.label, formula: 'taught reagent', k }); } });
    }
    return REAG;
  };
  window.Sketcher.attachSuggest(B, document.getElementById('rxnBSuggest'), () => predict(), reagentIndex);
  document.getElementById('rxnGo').addEventListener('click', predict);
  document.getElementById('rxnAdd').addEventListener('click', () => { const r = addRow(''); if (r) r.input.focus(); });
  document.getElementById('rxnFromSketch').addEventListener('click', () => {
    const s = window.Sketcher.smiles();
    if (!s) { out.innerHTML = '<div class="banner">Nothing is drawn on the Sketcher tab yet.</div>'; return; }
    setForm(s.split('.'), B.value); if (B.value.trim() || s.indexOf('.') > 0) predict(); else B.focus();
  });
  setForm([''], '');
  /* other pages hand a reaction over with Reactions.set([...reactants], 'reagent text') */
  window.Reactions = { predict, set: (reactants, reagentText) => { setForm(reactants, reagentText); return predict(); } };
})();
