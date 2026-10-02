/* Equation — the balanced chemical equation for a predicted reaction.

   Left: the reactants (+ the reagent as a molecule, where one is known: its SMILES from the reagent table,
   the partner molecule, or a net stand-in for oxidants and hydrides: [O], H₂, H–OH).
   Right: the product (+ the mechanistic engine's by-products). Whatever is still unbalanced is closed with
   small molecules (H₂O, HBr, NaBr …), one at a time or as a pair (water in, HOMgBr out). If it cannot be
   closed the card says what is left over instead of pretending.
     Equation.forRule(sub, res)                 a result from RXN.run / RXX.run
     Equation.forSmiles([reactants], product)   any reactants → product (the precedent-based predictor)      */
window.Equation = (() => {
  'use strict';
  const C = window.Chem, R = window.Reagents;
  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const el = (tag, cls, html) => { const n = document.createElement(tag); if (cls) n.className = cls; if (html != null) n.innerHTML = html; return n; };
  const nameOf = g => { try { return C.displayNameFor(g) || C.formula(g); } catch (e) { return '?'; } };
  const formulaOf = g => { try { return C.formula(g); } catch (e) { return ''; } };

  /* a species is {name, n, g?, t (atom tally per unit), f (formula)}; net reagents carry an explicit
     tally because a bare [O] would otherwise be given two implicit hydrogens by the valence rules */
  const species = (smiOrTally, name, extra = {}) => {
    if (typeof smiOrTally === 'string') { const g = C.parseSmiles(smiOrTally); return { g, t: C.atomTally([g]), f: formulaOf(g), name: name || nameOf(g), n: 1, ...extra }; }
    const f = Object.entries(smiOrTally).map(([e, k]) => e + (k > 1 ? k : '')).join('');
    return { t: smiOrTally, f, name: name || f, n: 1, ...extra };
  };
  const OX = { O: 1 }, H2 = '[H][H]';
  /* [source, label, formula?]: source = SMILES or an atom tally; null = a catalyst / nothing consumed.
     `net` reagents are written as what they contribute; their real by-products are left out of the sum. */
  const NET_REAGENT = {
    h2: [H2, 'H₂'], lindlar: [H2, 'H₂'], nabh4: [H2, 'H₂ (H⁻ from NaBH₄, then H⁺ in the work-up)'], lialh4: [H2, 'H₂ (H⁻ from LiAlH₄, then H⁺ in the work-up)'],
    dibal: [H2, 'H₂ (H⁻, then H⁺)'], nabh3cn: [H2, 'H₂ (H⁻, then H⁺)'], na_nh3: [H2, 'H₂ (2 e⁻ + 2 H⁺)'], nitrored: [H2, 'H₂ (or Fe / HCl)'],
    hydroboration: ['O', 'H–OH (BH₃, then H₂O₂ / NaOH)', 'H₂O'], alkynehydb: ['O', 'H–OH', 'H₂O'], oso4: ['OO', 'H₂O₂ (OsO₄ cat., NMO as the oxidant)'],
    pcc: [OX, '[O] (from Cr(VI))', '[O]'], kmno4: [OX, '[O] (from MnO₄⁻)', '[O]'], jones: [OX, '[O] (from Cr(VI))', '[O]'], dmp: [OX, '[O] (from I(V))', '[O]'],
    swern: [OX, '[O] (from DMSO)', '[O]'], tempo: [OX, '[O] (from NaOCl)', '[O]'], mno2: [OX, '[O] (from MnO₂)', '[O]'], ag2o: [OX, '[O] (from Ag₂O)', '[O]'], etard: [OX, '[O] (from CrO₂Cl₂)', '[O]'],
    naocl: [OX, '[O] (from OCl⁻)', '[O]'], ozonolysis: [{ O: 3 }, 'O₃ (then Me₂S)'], ozonolysis_ox: [{ O: 3 }, 'O₃ (then H₂O₂)'],
    hbr_peroxide: ['Br', 'HBr (ROOR initiator)'], aqacid: ['O', 'H₂O (H⁺ cat.)'], epox_h2o: ['O', 'H₂O (H⁺ cat.)'], alkynehyd: ['O', 'H₂O (Hg²⁺ cat.)'],
    hydration: ['O', 'H₂O (H⁺ cat.)'], epox_roh: ['CO', 'CH₃OH (H⁺ cat.)'], alkoxyhg: ['CO', 'CH₃OH (Hg(OAc)₂, then NaBH₄)'],
    clemmensen: [{ H: 4 }, '2 H₂ (Zn(Hg) / HCl)'], wolffkishner: [{ H: 4 }, '2 H₂ (N₂H₄ / KOH; N₂ leaves)'],
    allylbr: ['BrBr', 'Br₂ (delivered by NBS)'], nano2: ['O=N[O-]', 'HNO₂'],
    h2so4_dehydrate: [null, 'H₂SO₄ (catalyst)'], tautomer: [null, 'nothing — a proton moves'], thermal: [null, 'heat only'],
    acetal_off: ['O', 'H₂O (H⁺ cat.)'], cyclise: [null, 'nothing — intramolecular'],
    /* ionic reagents with their counter-ion, so the salt by-product balances */
    naoh: ['[Na+].[OH-]', 'NaOH', 'NaOH'], naoet: ['CC[O-].[Na+]', 'NaOEt', 'NaOEt'], tbuok: ['CC(C)(C)[O-].[K+]', 'tBuOK', 'tBuOK'], nacn: ['[C-]#N.[Na+]', 'NaCN', 'NaCN'],
  };
  /* short names for the reagents that come from the reagent table's SMILES */
  const REAGENT_SHORT = { hbr: ['HBr', 'HBr'], hcl: ['HCl', 'HCl'], hi: ['HI', 'HI'], br2: ['Br₂', 'Br₂'], cl2: ['Cl₂', 'Cl₂'], ar_br2: ['Br₂ (FeBr₃ cat.)', 'Br₂'],
    mcpba: ['mCPBA', 'C₇H₅ClO₃'], pbr3: ['PBr₃', 'PBr₃'], nitration: ['HNO₃ (H₂SO₄ cat.)', 'HNO₃'], sulfonation: ['SO₃ (in H₂SO₄)', 'SO₃'], naocl: ['NaOCl', 'NaOCl'],
    etard: ['CrO₂Cl₂', 'CrO₂Cl₂'], soCl2: ['SOCl₂', 'SOCl₂'] };
  /* small molecules used to close the balance: [smiles, name, formula, cost]; H₂ and O₂ are last resorts */
  const SMALL = [['O', 'water', 'H₂O'], ['Br', 'HBr', 'HBr'], ['Cl', 'HCl', 'HCl'], ['I', 'HI', 'HI'], ['F', 'HF', 'HF'], ['O=C=O', 'CO₂', 'CO₂'], ['N#N', 'N₂', 'N₂'], ['N', 'NH₃', 'NH₃'],
    ['CO', 'methanol', 'CH₃OH'], ['CCO', 'ethanol', 'C₂H₅OH'], ['CC(=O)O', 'acetic acid', 'CH₃CO₂H'], ['CC(C)(C)O', 'tert-butanol', 'C₄H₉OH'],
    ['[Na+].[Br-]', 'NaBr', 'NaBr'], ['[Na+].[Cl-]', 'NaCl', 'NaCl'], ['[Na+].[I-]', 'NaI', 'NaI'], ['[K+].[Br-]', 'KBr', 'KBr'], ['[K+].[Cl-]', 'KCl', 'KCl'], ['[K+].[I-]', 'KI', 'KI'],
    ['[Na+].[OH-]', 'NaOH', 'NaOH'], ['[Li+].[OH-]', 'LiOH', 'LiOH'], ['OS(=O)(=O)O', 'H₂SO₄', 'H₂SO₄'], ['O=[N+]([O-])O', 'HNO₃', 'HNO₃'], ['O=S=O', 'SO₂', 'SO₂'],
    ['O[Mg]Br', 'HOMgBr', 'HOMgBr'], ['O[Mg]Cl', 'HOMgCl', 'HOMgCl'], ['O[Mg]I', 'HOMgI', 'HOMgI'], ['OP(=O)(O)O', 'H₃PO₄', 'H₃PO₄'], ['OP(O)O', 'H₃PO₃', 'H₃PO₃'],
    ['OB(O)O', 'B(OH)₃', 'B(OH)₃'], ['OB(O)Br', 'B(OH)₂Br', 'B(OH)₂Br'], ['OB(O)Cl', 'B(OH)₂Cl', 'B(OH)₂Cl'], ['OB(O)I', 'B(OH)₂I', 'B(OH)₂I'],
    ['OO', 'H₂O₂', 'H₂O₂'], [H2, 'H₂', 'H₂', 2], ['O=O', 'O₂', 'O₂', 2]];
  const sub2 = s => String(s).replace(/(\d+)/g, m => m.replace(/\d/g, d => '₀₁₂₃₄₅₆₇₈₉'[d]));
  const sumT = xs => { const o = {}; xs.forEach(x => { for (const e in x.t) o[e] = (o[e] || 0) + x.t[e] * x.n; }); return o; };
  let smalls = null;
  const diffOf = (L, Rt) => {
    const a = sumT(L), b = sumT(Rt), d = {};
    new Set([...Object.keys(a), ...Object.keys(b)]).forEach(e => { const v = (a[e] || 0) - (b[e] || 0); if (v) d[e] = v; });
    return d;                                    // positive = surplus on the left, so needed on the right
  };
  const fits = (sm, d, sign) => {                // how many of `sm` fit into the gap on one side
    const se = Object.keys(sm.t);
    if (!se.every(e => Math.sign(d[e] || 0) === sign)) return 0;
    return Math.min(...se.map(e => Math.floor(Math.abs(d[e]) / sm.t[e])));
  };
  function close(L, Rt) {                        // → {L, Rt, cost} or null
    if (!smalls) smalls = SMALL.map(([s, name, f, cost]) => { try { const x = species(s, name); x.f = f; x.cost = cost || 1; return x; } catch (e) { return null; } }).filter(Boolean);
    L = L.slice(); Rt = Rt.slice(); let cost = 0;
    for (let round = 0; round < 3; round++) {
      const d = diffOf(L, Rt); if (!Object.keys(d).length) return { L, Rt, cost };
      let best = null;
      for (const sm of smalls) for (const sign of [1, -1]) {
        const n = fits(sm, d, sign); if (n < 1) continue;
        const covered = Object.keys(sm.t).reduce((s, e) => s + n * sm.t[e], 0);
        if (!best || covered > best.covered) best = { adds: [{ sm, n, sign }], covered, cost: sm.cost };
      }
      if (!best) {                               // a pair, one on each side, when nothing fits alone
        for (const a of smalls) for (let na = 1; na <= 2 && !best; na++) {
          const d2 = { ...d }; for (const e in a.t) d2[e] = (d2[e] || 0) + na * a.t[e];
          for (const b of smalls) { const nb = fits(b, d2, 1); if (nb < 1) continue;
            const rest = { ...d2 }; for (const e in b.t) rest[e] = (rest[e] || 0) - nb * b.t[e];
            if (Object.values(rest).every(v => !v)) { best = { adds: [{ sm: a, n: na, sign: -1 }, { sm: b, n: nb, sign: 1 }], covered: 0, cost: a.cost + b.cost }; break; } }
        }
      }
      if (!best) return null;
      best.adds.forEach(x => (x.sign > 0 ? Rt : L).push({ ...x.sm, n: x.n, added: true }));
      cost += best.cost;
    }
    return Object.keys(diffOf(L, Rt)).length ? null : { L, Rt, cost };
  }
  /* the first reactant and the product scale together (3 ROH + PBr₃ → 3 RBr + H₃PO₃); the cheapest closure
     wins — fewest added species, and a bigger stoichiometry only when it removes them */
  function balance(left, right, reagent) {
    let best = null;
    for (let s = 1; s <= 3; s++) {
      left[0].n = s; right[0].n = s;
      for (let coef = 1; coef <= 4; coef++) {
        if (reagent) reagent.n = coef; else if (coef > 1) break;
        const sol = close(left, right); if (!sol) continue;
        const cost = sol.cost + (s - 1) * 0.5 + (coef - 1) * 0.1;
        if (!best || cost < best.cost) best = { ...sol, cost, s, coef };
      }
    }
    if (best) { left[0].n = right[0].n = best.s; if (reagent) reagent.n = best.coef; return { ok: true, L: best.L, Rt: best.Rt }; }
    left[0].n = right[0].n = 1; if (reagent) reagent.n = 1;
    return { ok: false, L: left, Rt: right };
  }
  function card(left, right, { reagent = null, catalyst = null, skipped = [], note = '' } = {}) {
    const box = el('div', 'card');
    const { ok, L, Rt } = balance(left, right, reagent);
    const side = xs => xs.map(x => (x.n > 1 ? x.n + ' ' : '') + sub2(x.f)).join(' + ');
    const names = xs => xs.map(x => (x.n > 1 ? x.n + ' × ' : '') + x.name + (x.added ? ' (balancing)' : '')).join(' + ');
    box.appendChild(el('div', null, `<div style="font:500 17px system-ui,sans-serif">${esc(side(L))} ⟶ ${esc(side(Rt))}` +
      (catalyst ? ` <span class="hint" style="font-size:13px">[${esc(catalyst)}]</span>` : '') + '</div>' +
      `<div class="hint" style="margin-top:4px">${esc(names(L))} ⟶ ${esc(names(Rt))}</div>`));
    const d = diffOf(L, Rt), tally = sumT(L);
    box.appendChild(el('div', 'hint', ok
      ? '✓ Balanced: ' + Object.entries(tally).map(([e, n]) => e + ' ' + n).join(', ') + ' on each side.' +
        (Rt.some(x => x.added) || L.some(x => x.added) ? ' Species marked "balancing" are inferred from the atom count.' : '') +
        (reagent && reagent.net ? ' The reagent is written as its net contribution; the real reagent and work-up are in the conditions' + (skipped.length ? ' (its own by-products: ' + esc(skipped.join(', ')) + ')' : '') + '.' : '') + (note ? ' ' + esc(note) : '')
      : '⚠ Not fully balanced — left over: ' + Object.entries(d).map(([e, v]) => e + ' ' + (v > 0 ? '+' : '') + v).join(', ') +
        '. Only the organic product is modelled, not every inorganic by-product.'));
    return box;
  }

  function forRule(sub, res) {
    const id = res.rule.id; let reagent = null, catalyst = null;
    try {
      if (res.rule.eq) { reagent = species(res.rule.eq.smiles, res.rule.eq.label, { reagent: true }); reagent.f = res.rule.eq.f; }   // an added family says what to write
      else if (res.rule.molecule) reagent = species(C.toSmiles(res.rule.molecule), nameOf(res.rule.molecule), { reagent: true });
      else if (NET_REAGENT[id]) { const [src, lab, f] = NET_REAGENT[id]; if (src) { reagent = species(src, lab, { reagent: true, net: true }); if (f) reagent.f = f; } else catalyst = lab; }
      else { const t = R.RULE_REAGENTS[id]; if (t && t.smiles) { const sh = REAGENT_SHORT[id]; reagent = species(t.smiles, sh ? sh[0] : res.rule.label.split(' (')[0], { reagent: true }); if (sh) reagent.f = sh[1]; } }
    } catch (e) { reagent = null; }
    const left = [species(C.toSmiles(sub), nameOf(sub))]; if (reagent) left.push(reagent);
    const right = [species(C.toSmiles(res.major.graph), res.major.name)];
    /* the engine's by-products belong to the real reagent; when the reagent is written as a net
       contribution ([O], H₂) they would unbalance the sum, so they are named in the note instead */
    const skipped = [];
    (res.byproducts || []).forEach(b => {
      if (reagent && reagent.net) { if (b.name) skipped.push(b.name); return; }
      if (b.smiles) { try { right.push(species(b.smiles, b.name || null)); } catch (e) {} }
    });
    return card(left, right, { reagent, catalyst, skipped });
  }
  function forSmiles(reactants, product, productName) {
    try {
      const left = reactants.map(s => species(s)), right = [species(product, productName || null)];
      return card(left, right, { note: 'Reagents that do not end up in the product are not part of this sum.' });
    } catch (e) { return null; }
  }
  return { forRule, forSmiles };
})();
