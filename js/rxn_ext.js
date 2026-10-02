/* RXX — reaction families added on top of the lifted engine, chosen from the benchmark's misses.

   The lifted engine (js/rxn_engine.js) is left verbatim; this layer owns the families it lacks. Each
   family matches a substrate + partner (or a substrate + a written reagent), builds the product graph
   with the Chem helpers, and returns a result in the same shape as RXN.run so the Reactions page and
   the benchmark can treat it identically:
     { rule: {id, label, molecule?}, sub, sites, outcomes: [{graph, share, tag, name, formula, smiles, label, why, stereoHow}],
       major, byproducts, notes, summary: {chemo, site, regio, stereo, conditions}, ext: true }
   Families (why each is here is in its `why`):
     carbamate   R–N=C=O + R'–OH  → R–NH–C(=O)–OR'          urea with an amine
     cyanohydrin R2C=O + CN⁻      → R2C(OH)–C≡N              (NaCN / KCN, then acid)
     pyridazinone γ-keto acid/ester + N2H4 → 6-aryl-4,5-dihydropyridazin-3(2H)-one
     pyrazole    1,3-diketone + N2H4       → 3,5-disubstituted 1H-pyrazole (Knorr)
   RXX.run(A, B) tries every family on the pair (both orders); RXX.runReagent(A, ruleId) handles the
   written reagents (nacn → cyanohydrin, hydrazine → the cyclocondensations).                          */
(() => {
  'use strict';
  const C = window.Chem;
  if (!C) return;
  const atom = (g, id) => g.atoms.find(a => a.id === id);
  const nb = (g, id) => C.neighbors(g, id).filter(n => n.atom);
  const bond = (g, a, b) => C.bondBetween(g, a, b);
  const hOn = (g, a) => C.implicitH(g, a);
  const isCarbonylO = (g, o) => { const n = nb(g, o.id); return n.length === 1 && n[0].bond.order === 2 && n[0].atom.element === 'C'; };

  /* ---- graph helpers ---- */
  function merge(A, B) {                              // one graph holding both, B's ids shifted; returns [g, mapB]
    const g = C.clone(A); const off = Math.max(0, ...g.atoms.map(a => a.id)); const map = new Map();
    B.atoms.forEach(a => { map.set(a.id, a.id + off); g.atoms.push({ ...a, id: a.id + off, x: a.x + 200, y: a.y }); });
    B.bonds.forEach(b => g.bonds.push({ ...b, a: b.a + off, b: b.b + off, ...(b.narrow !== undefined ? { narrow: b.narrow + off } : {}) }));
    g.nextId = Math.max(...g.atoms.map(a => a.id)) + 1;
    return [g, map];
  }
  function cutFragment(g, keepId, cutId) {            // remove cutId and everything reachable from it without passing keepId
    const gone = new Set([cutId]); const stack = [cutId];
    while (stack.length) { const id = stack.pop(); nb(g, id).forEach(n => { if (n.atom.id !== keepId && !gone.has(n.atom.id)) { gone.add(n.atom.id); stack.push(n.atom.id); } }); }
    g.atoms = g.atoms.filter(a => !gone.has(a.id)); g.bonds = g.bonds.filter(b => !gone.has(b.a) && !gone.has(b.b));
  }
  const setOrder = (g, a, b, order) => { const x = bond(g, a, b); if (x) x.order = order; else g.bonds.push({ a, b, order }); };
  const finish = g => { try { C.layoutGraph(g); } catch (e) {} return g; };

  /* ---- site finders ---- */
  const isocyanates = g => g.atoms.filter(c => c.element === 'C' && nb(g, c.id).length === 2 &&
    nb(g, c.id).every(n => n.bond.order === 2) && nb(g, c.id).some(n => n.atom.element === 'N') && nb(g, c.id).some(n => n.atom.element === 'O'))
    .map(c => ({ c, n: nb(g, c.id).find(n => n.atom.element === 'N').atom, o: nb(g, c.id).find(n => n.atom.element === 'O').atom }));
  const alcohols = g => g.atoms.filter(o => o.element === 'O' && !o.charge && hOn(g, o) === 1 && nb(g, o.id).length === 1 &&
    nb(g, o.id)[0].atom.element === 'C' && !nb(g, nb(g, o.id)[0].atom.id).some(n => n.atom.element === 'O' && n.bond.order === 2));
  const amines = g => g.atoms.filter(n => n.element === 'N' && !n.charge && hOn(g, n) >= 1 && nb(g, n.id).every(x => x.bond.order === 1 && x.atom.element === 'C') &&
    !nb(g, n.id).some(x => nb(g, x.atom.id).some(y => y.atom.element === 'O' && y.bond.order === 2)));
  const carbonyls = g => g.atoms.filter(c => c.element === 'C' && nb(g, c.id).some(n => n.atom.element === 'O' && n.bond.order === 2 && nb(g, n.atom.id).length === 1))
    .map(c => { const o = nb(g, c.id).find(n => n.atom.element === 'O' && n.bond.order === 2).atom; const others = nb(g, c.id).filter(n => n.atom.id !== o.id);
      const kind = others.some(n => ['O', 'N', 'Cl', 'S'].includes(n.atom.element)) ? 'acyl' : others.length === 1 ? 'aldehyde' : 'ketone';
      return { c, o, kind, others: others.map(n => n.atom) }; });
  const isCyanide = g => g.atoms.length === 2 && g.atoms.some(a => a.element === 'C') && g.atoms.some(a => a.element === 'N') && g.bonds[0] && g.bonds[0].order === 3;
  const isHydrazine = g => g.atoms.length === 2 && g.atoms.every(a => a.element === 'N') && g.bonds[0] && g.bonds[0].order === 1;
  const acidOrEsterC = (g, c) => {                    // the acyl carbon of an acid or ester, with its leaving oxygen
    const o = nb(g, c.id).find(n => n.atom.element === 'O' && n.bond.order === 1 && (hOn(g, n.atom) === 1 || nb(g, n.atom.id).length === 2));
    return o ? o.atom : null;
  };
  const chainBetween = (g, a, b, len) => {            // a–CH2–…–b with `len` sp3 carbons between, or null
    const path = [a.id]; const walk = (id, depth, prev) => {
      if (depth === len) return bond(g, id, b.id) ? true : false;
      for (const n of nb(g, id)) { if (n.atom.id === prev || n.atom.element !== 'C' || n.bond.order !== 1) continue;
        if (nb(g, n.atom.id).some(x => x.bond.order > 1 && x.atom.id !== id)) continue; path.push(n.atom.id); if (walk(n.atom.id, depth + 1, id)) return true; path.pop(); }
      return false; };
    return walk(a.id, 0, -1) ? path.slice(1) : null;
  };

  /* ---- the families ---- */
  const FAMILIES = [
    { id: 'ext:carbamate', label: 'isocyanate + alcohol / amine (carbamate or urea)', needs: 'an isocyanate and an alcohol or amine',
      match(A, B) { const iso = isocyanates(A); if (!iso.length) return null; const oh = alcohols(B), nh = amines(B); if (!oh.length && !nh.length) return null; return { iso: iso[0], nuc: oh[0] || nh[0], kind: oh.length ? 'carbamate' : 'urea' }; },
      build(A, B, m) { const [g, map] = merge(A, B); setOrder(g, m.iso.c.id, m.iso.n.id, 1); setOrder(g, m.iso.c.id, map.get(m.nuc.id), 1); return finish(g); },
      why: m => `The isocyanate carbon is a strong electrophile (it sits between two heteroatoms, N=C=O); the ${m.kind === 'carbamate' ? 'alcohol' : 'amine'} adds to it and the proton moves to nitrogen. No catalyst is needed; a tertiary amine or dibutyltin dilaurate speeds it up. Product: a ${m.kind}.`,
      conditions: 'room temperature, aprotic solvent (CH₂Cl₂, THF, toluene); Et₃N or DABCO catalytic', stereo: 'No new stereocentre.' },
    { id: 'ext:cyanohydrin', label: 'cyanide addition to a carbonyl (cyanohydrin)', needs: 'an aldehyde or ketone and cyanide',
      match(A, B) { if (!isCyanide(B)) return null; const cs = carbonyls(A).filter(x => x.kind !== 'acyl'); if (!cs.length) return null; cs.sort((p, q) => (p.kind === 'aldehyde' ? 0 : 1) - (q.kind === 'aldehyde' ? 0 : 1)); return { site: cs[0], cn: B.atoms.find(a => a.element === 'C') }; },
      build(A, B, m) { const [g, map] = merge(A, B); setOrder(g, m.site.c.id, m.site.o.id, 1); const cn = atom(g, map.get(m.cn.id)); cn.charge = 0; setOrder(g, m.site.c.id, cn.id, 1); g.atoms.forEach(a => { if (a.element === 'N' && a.charge) a.charge = 0; }); return finish(g); },
      why: m => `Cyanide is a good carbon nucleophile and the C=O carbon is electrophilic; CN⁻ adds to the carbonyl carbon and the alkoxide is protonated on work-up. ${m.site.kind === 'aldehyde' ? 'An aldehyde is less hindered than a ketone and reacts first.' : 'Ketones react, more slowly than aldehydes; hindered ketones may not.'} The equilibrium favours the cyanohydrin for aldehydes and unhindered ketones.`,
      conditions: 'NaCN or KCN with a little acid (pH ≈ 8–10, HCN in situ), 0–25 °C; or TMSCN then aqueous acid', stereo: 'A new stereocentre is made from a flat carbonyl with no facial bias: the cyanohydrin is racemic.',
      eq: { smiles: 'C#N', label: 'HCN (from NaCN + acid)', f: 'HCN' } },
    { id: 'ext:pyridazinone', label: 'hydrazine + γ-keto acid / ester (4,5-dihydropyridazin-3(2H)-one)', needs: 'a 1,4-keto acid or ester and hydrazine',
      match(A, B) { if (!isHydrazine(B)) return null;
        for (const k of carbonyls(A).filter(x => x.kind === 'ketone' || x.kind === 'aldehyde')) for (const a of carbonyls(A).filter(x => x.kind === 'acyl')) {
          const lo = acidOrEsterC(A, a.c); if (!lo) continue; const chain = chainBetween(A, k.c, a.c, 2); if (chain) return { k, a, lo, chain }; }
        return null; },
      build(A, B, m) { const [g, map] = merge(A, B); const [n1, n2] = B.atoms.map(a => map.get(a.id));
        cutFragment(g, m.k.c.id, m.k.o.id); setOrder(g, m.k.c.id, n1, 2);          // hydrazone end
        cutFragment(g, m.a.c.id, m.lo.id); setOrder(g, m.a.c.id, n2, 1);          // hydrazide end closes the ring
        return finish(g); },
      why: () => 'Hydrazine has two nucleophilic nitrogens two atoms apart. One condenses with the ketone to a hydrazone; the other then attacks the acid or ester carbonyl in the same molecule, four atoms away, and closes a six-membered ring: a 4,5-dihydropyridazin-3(2H)-one. This is the standard route to that heterocycle (levosimendan-type cardiotonics are made this way). The open-chain hydrazide the engine used to stop at is the intermediate.',
      conditions: 'hydrazine hydrate, ethanol or acetic acid, reflux', stereo: 'C4 or C5 may become a stereocentre if substituted; unsubstituted, the product is achiral.',
      eq: { smiles: 'NN', label: 'hydrazine', f: 'N₂H₄' }, water: 2 },
    { id: 'ext:pyrazole', label: 'hydrazine + 1,3-diketone (Knorr pyrazole)', needs: 'a 1,3-dicarbonyl and hydrazine',
      match(A, B) { if (!isHydrazine(B)) return null; const ks = carbonyls(A).filter(x => x.kind === 'ketone' || x.kind === 'aldehyde');
        for (const p of ks) for (const q of ks) { if (p === q) continue; const chain = chainBetween(A, p.c, q.c, 1); if (chain) return { p, q, mid: chain[0] }; } return null; },
      build(A, B, m) { const [g, map] = merge(A, B); const [n1, n2] = B.atoms.map(a => map.get(a.id));
        cutFragment(g, m.p.c.id, m.p.o.id); cutFragment(g, m.q.c.id, m.q.o.id);
        setOrder(g, m.p.c.id, n1, 2); setOrder(g, m.q.c.id, n2, 1); setOrder(g, m.q.c.id, m.mid, 2); setOrder(g, m.p.c.id, m.mid, 1);
        return finish(g); },
      why: () => 'The Knorr pyrazole synthesis: each hydrazine nitrogen condenses with one carbonyl of the 1,3-diketone; the second condensation closes a five-membered ring, and loss of the second water gives the aromatic pyrazole. With two different R groups the two tautomers are the same compound (1H-pyrazole NH exchanges).',
      conditions: 'hydrazine hydrate, ethanol, reflux; acetic acid catalytic', stereo: 'Aromatic product, no stereocentre.',
      eq: { smiles: 'NN', label: 'hydrazine', f: 'N₂H₄' }, water: 2 },
  ];

  function result(fam, A, B, m) {
    let graph; try { graph = fam.build(A, B, m); } catch (e) { return null; }
    try { C.validateValences(graph); } catch (e) { return null; }
    const name = (() => { try { return C.displayNameFor(graph) || C.formula(graph); } catch (e) { return C.formula(graph); } })();
    const o = { graph, share: 1, tag: 'major', name, formula: C.formula(graph), smiles: C.toSmiles(graph), label: fam.label, why: fam.why(m), stereoHow: fam.stereo, stereo: null };
    try { if (graph.atoms.some(a => a.element === 'C' && C.isStereocenter(graph, a) && !A.atoms.some(x => x.id === a.id && C.isStereocenter(A, x)))) o.stereo = { kind: 'racemic', forms: [], note: fam.stereo }; } catch (e) {}
    return { rule: { id: fam.id, label: fam.label, molecule: B, eq: fam.eq || null }, ruleId: fam.id, sub: A, ext: true,
      sites: [{ label: fam.needs, share: 1, reason: 'the only site this family acts on' }], outcomes: [o], major: o,
      byproducts: fam.water ? Array(fam.water).fill({ smiles: 'O', name: 'water (condensation)' }) : [],
      notes: [], summary: { chemo: fam.label + '.', site: '', regio: fam.why(m), stereo: fam.stereo, conditions: 'Conditions: ' + fam.conditions + '.' }, env: { T: 298, TAuto: true, excess: false } };
  }
  /* try every family on the pair, both ways round; the first match wins (families are ordered specific → general) */
  function run(A, B) {
    if (!A || !B || !A.atoms || !B.atoms || !A.atoms.length || !B.atoms.length) return null;
    for (const [X, Y] of [[A, B], [B, A]]) for (const fam of FAMILIES) { let m = null; try { m = fam.match(X, Y); } catch (e) {} if (m) { const r = result(fam, X, Y, m); if (r) return r; } }
    return null;
  }
  /* a written reagent that one of these families understands */
  const REAGENT_MOL = { nacn: '[C-]#N', wolffkishner: 'NN', hydrazine: 'NN' };
  function runReagent(A, ruleId) {
    const smi = REAGENT_MOL[ruleId]; if (!smi) return null;
    let B; try { B = C.parseSmiles(smi); } catch (e) { return null; }
    return run(A, B);
  }
  window.RXX = { run, runReagent, FAMILIES };
})();
