/* RXN: the reaction selectivity engine — RXN.run(graph, ruleId, env).
   Lifted verbatim from organic-chem-tool/src/rxn_engine.js; it depends only on the Chem global (js/chem.js).
   Do not hand-edit; re-copy from the old tool if it changes. */

/* =====================================================================
   REACTION SELECTIVITY ENGINE (RXN) — why a reagent goes where it goes.

   The taught rules in Chem.RULES know HOW each reagent changes a group.
   This layer decides WHERE, when a molecule offers a choice:

     chemo-selectivity  — which KIND of group the reagent attacks and why the
                          others are left alone (NaBH4 takes the ketone, not
                          the ester; PCC takes the alcohol, not the alkene);
     site selectivity   — which ONE of several identical groups reacts first
                          (the more substituted alkene for HBr, the least
                          hindered one for BH3, the 1° alcohol for PBr3);
     regio-selectivity  — which way round a group adds (Markovnikov /
                          anti-Markovnikov, Zaitsev / Hofmann);
     stereo             — racemic pairs, syn/anti, inversion, meso;
     conditions         — temperature, pH, pressure, solvent, light and
                          reagent equivalents, which move the numbers.

   Every relative rate is a textbook-grade estimate, and the numbers that
   come from a measured series say so in the text they generate. The
   engine never invents a mechanism: each family below mirrors the rule
   card it belongs to, and the products it draws are built with the same
   graph helpers the rules use.
   ===================================================================== */
const RXN = (() => {
  const C = () => Chem;                       // the shared engine (late-bound)
  const nb = (g, id) => C().neighbors(g, id).filter(n => n.atom);
  const deg = (g, a) => C().carbonDegree(g, a);                      // carbon neighbours
  const hOn = (g, a) => C().implicitH(g, a);
  const atom = (g, id) => g.atoms.find(a => a.id === id);
  const bondOf = (g, a, b) => C().bondBetween(g, a, b);
  const clone = g => C().clone(g);
  const log10 = Math.log10;

  /* ---------- environment ---------- */
  const DEFAULT_ENV = { T: null, pH: null, P: 1, solvent: 'auto', light: false, equiv: 1 };
  /* what the reagent itself does to the medium — used when pH/solvent are 'auto' */
  const REGIME = {
    hbr: { pH: 0, solvent: 'nonpolar' }, hcl: { pH: 0, solvent: 'nonpolar' }, hi: { pH: 0, solvent: 'nonpolar' },
    hydration: { pH: 0.5, solvent: 'water' }, h2so4_dehydrate: { pH: -1, solvent: 'neat', T: 413 },
    hbr_peroxide: { pH: 1, solvent: 'nonpolar' }, br2: { pH: 7, solvent: 'nonpolar' }, cl2: { pH: 7, solvent: 'nonpolar' },
    mcpba: { pH: 6, solvent: 'nonpolar' }, oso4: { pH: 7, solvent: 'protic' }, hydroboration: { pH: 7, solvent: 'ether' },
    h2: { pH: 7, solvent: 'protic' }, pcc: { pH: 4, solvent: 'nonpolar' }, kmno4: { pH: 10, solvent: 'water' },
    pbr3: { pH: 4, solvent: 'ether', T: 273 }, naoh: { pH: 14, solvent: 'water' }, naoet: { pH: 15, solvent: 'protic', T: 328 },
    tbuok: { pH: 17, solvent: 'protic', T: 328 }, nacn: { pH: 11, solvent: 'aprotic' }, nabh4: { pH: 9, solvent: 'protic' },
    lialh4: { pH: 7, solvent: 'ether' },
    sulfonation: { pH: -2, solvent: 'neat', T: 313 }, nitration: { pH: -2, solvent: 'neat', T: 303 }, ar_br2: { pH: 3, solvent: 'nonpolar' },
    etard: { pH: 7, solvent: 'nonpolar', T: 273 }, naocl: { pH: 5, solvent: 'water', T: 293 },
    /* the modern catalog — each reagent brings the conditions it is actually run at, so a
       gate like "the Swern must stay below −60 °C" is not tripped by the room-temperature default */
    dmp: { pH: 6, solvent: 'nonpolar', T: 298 }, swern: { pH: 7, solvent: 'nonpolar', T: 195 },
    tempo: { pH: 8, solvent: 'water', T: 273 }, jones: { pH: 0, solvent: 'water', T: 293 },
    nano2: { pH: 1, solvent: 'water', T: 275 }, sandmeyer_cl: { pH: 1, solvent: 'water', T: 278 }, sandmeyer_br: { pH: 1, solvent: 'water', T: 278 },
    sandmeyer_cn: { pH: 7, solvent: 'water', T: 278 }, diazo_i: { pH: 1, solvent: 'water', T: 278 }, diazo_f: { pH: 1, solvent: 'water', T: 373 },
    diazo_oh: { pH: 1, solvent: 'water', T: 353 }, diazo_h: { pH: 1, solvent: 'water', T: 278 },
    soCl2: { pH: 7, solvent: 'nonpolar', T: 350 }, pcl5: { pH: 7, solvent: 'nonpolar', T: 298 },
    ch2n2: { pH: 7, solvent: 'ether', T: 273 },
    ar_cl2: { pH: 7, solvent: 'nonpolar', T: 298 }, ar_i2: { pH: 1, solvent: 'nonpolar', T: 298 },
    epox_h2o: { pH: 1, solvent: 'water', T: 298 }, epox_hx: { pH: 1, solvent: 'nonpolar', T: 273 }, epox_roh: { pH: 1, solvent: 'protic', T: 323 },
    clemmensen: { pH: 0, solvent: 'water', T: 373 }, wolffkishner: { pH: 14, solvent: 'protic', T: 453 },
    ozonolysis: { pH: 7, solvent: 'nonpolar', T: 195 }, ozonolysis_ox: { pH: 7, solvent: 'nonpolar', T: 195 },
    alkynehyd: { pH: 1, solvent: 'water', T: 333 }, alkynehydb: { pH: 7, solvent: 'ether', T: 298 },
    allylbr: { pH: 7, solvent: 'nonpolar', T: 350 }, aqacid: { pH: 1, solvent: 'water', T: 298 },
    ag2o: { pH: 8, solvent: 'water', T: 298 }, alkoxyhg: { pH: 7, solvent: 'protic', T: 298 },
    radpoly: { pH: 7, solvent: 'neat', T: 373 }, cationic: { pH: 1, solvent: 'nonpolar', T: 250 },
    anionic: { pH: 14, solvent: 'aprotic', T: 273 }, ziegler: { pH: 7, solvent: 'nonpolar', T: 333 },
    nabh3cn: { pH: 6, solvent: 'protic', T: 298 },
    tautomer: { pH: 7, solvent: 'protic', T: 298 }, lda: { pH: 14, solvent: 'ether', T: 195 },
    enolate_base: { pH: 13, solvent: 'protic', T: 298 }, thermal: { pH: 7, solvent: 'neat', T: 403 },
    mno2: { pH: 7, solvent: 'nonpolar', T: 298 }, dibal: { pH: 7, solvent: 'nonpolar', T: 195 },
    lindlar: { pH: 7, solvent: 'protic', T: 298 }, na_nh3: { pH: 11, solvent: 'aprotic', T: 195 },
    nitrored: { pH: 7, solvent: 'protic', T: 298 },
    boc_on: { pH: 9, solvent: 'nonpolar', T: 298 }, boc_off: { pH: 0, solvent: 'nonpolar', T: 298 },
    cbz_on: { pH: 9, solvent: 'water', T: 273 }, cbz_off: { pH: 7, solvent: 'protic', T: 298 },
    fmoc_on: { pH: 9, solvent: 'water', T: 273 }, fmoc_off: { pH: 11, solvent: 'aprotic', T: 298 },
    tbs_on: { pH: 8, solvent: 'aprotic', T: 298 }, tbs_off: { pH: 7, solvent: 'ether', T: 298 },
    bn_on: { pH: 14, solvent: 'aprotic', T: 273 }, bn_off: { pH: 7, solvent: 'protic', T: 298 },
    thp_on: { pH: 4, solvent: 'nonpolar', T: 298 }, ac_on: { pH: 8, solvent: 'nonpolar', T: 298 },
    ms_on: { pH: 8, solvent: 'nonpolar', T: 273 }, ts_on: { pH: 8, solvent: 'nonpolar', T: 273 },
    acetal_on: { pH: 3, solvent: 'nonpolar', T: 383 }, acetal_off: { pH: 1, solvent: 'water', T: 298 },
  };
  const SOLVENT_NAMES = { water: 'water', protic: 'a protic solvent (EtOH / MeOH)', aprotic: 'a polar aprotic solvent (DMSO / DMF)',
    nonpolar: 'a non-polar solvent (CH₂Cl₂ / hexane)', ether: 'dry ether / THF', neat: 'neat (no solvent)' };
  function resolveEnv(ruleId, env) {
    const e = Object.assign({}, DEFAULT_ENV, env || {});
    const reg = REGIME[ruleId] || { pH: 7, solvent: 'nonpolar' };
    e.pHAuto = (e.pH === null || e.pH === undefined || e.pH === 'auto');
    e.pHVal = e.pHAuto ? reg.pH : +e.pH;
    e.solventAuto = (!e.solvent || e.solvent === 'auto');
    e.solventVal = e.solventAuto ? reg.solvent : e.solvent;
    e.TAuto = (e.T === null || e.T === undefined || e.T === 'auto');
    e.T = e.TAuto ? (reg.T || 298) : (+e.T || 298); e.P = +e.P || 1;
    /* a plain "heat" / Δ means the reagent's usual heat — never LESS than the temperature its regime calls for */
    if (e.heatAuto && reg.T && e.T < reg.T) e.T = reg.T;
    e.excess = e.equiv === 'excess' || (+e.equiv >= 2);
    e.regime = reg;
    return e;
  }
  const degC = T => Math.round(T - 273.15);

  /* ---------- small structural helpers ---------- */
  function bondInRing(g, a, b) {                  // still connected with the bond removed?
    const seen = new Set([a]); const st = [a];
    while (st.length) {
      const id = st.pop();
      for (const n of nb(g, id)) {
        if (id === a && n.atom.id === b) continue;
        if (n.atom.id === b) return true;
        if (!seen.has(n.atom.id)) { seen.add(n.atom.id); st.push(n.atom.id); }
      }
    }
    return false;
  }
  function hasPiToCarbon(g, id, except) {          // this atom carries a C=C / C≡C / aromatic bond (not to `except`)
    return nb(g, id).some(n => n.atom.element === 'C' && n.bond.order > 1 && n.atom.id !== except);
  }
  function isAromatic(g, id) { try { return C().isAromaticCarbon(g, id); } catch (e) { return false; } }
  /* is this atom the heteroatom of an aromatic FIVE-membered ring (thiophene, pyrrole, furan)?
     If so its lone pair belongs to the ring's six π electrons and is not available — a
     thiophene is not oxidised to a sulfoxide, and a pyrrole nitrogen is not basic. */
  function aromFiveHetero(g, id) {
    try {
      const sys = C().ringSystemsOf(g, C().ringAtomsOf(g)).find(r => r.indexOf(id) >= 0);
      if (!sys) return false;
      const ns = nb(g, id);
      if (!(ns.length === 2 && ns.every(n => n.atom.element === 'C' && isAromatic(g, n.atom.id)))) return false;
      return C().ringsOfSystem(g, sys).some(r => r.length === 5 && r.indexOf(id) >= 0);
    } catch (e) { return false; }
  }
  /* a pyridine-type ring nitrogen: no hydrogen, no charge, two aromatic-carbon neighbours, and
     NOT the heteroatom of a five-ring (that would be a pyrrole) */
  function pyridineN(g, a) {
    if (a.element !== 'N' || a.charge || hOn(g, a) !== 0) return false;
    const ns = nb(g, a.id);
    return ns.length === 2 && ns.every(n => n.atom.element === 'C' && isAromatic(g, n.atom.id)) && !aromFiveHetero(g, a.id);
  }
  /* a sulfide: divalent neutral sulfur between two carbons, with no oxygen already on it */
  function sulfideAt(g, a) {
    if (a.element !== 'S' || a.charge) return false;
    const ns = nb(g, a.id);
    return ns.length === 2 && ns.every(n => n.bond.order === 1 && n.atom.element === 'C') &&
      !ns.some(n => n.atom.element === 'O') && !aromFiveHetero(g, a.id);
  }
  /* a tertiary aliphatic amine (an N-oxide precursor), not an amide */
  function tertAmine(g, a) {
    if (a.element !== 'N' || a.charge || hOn(g, a) !== 0) return false;
    const ns = nb(g, a.id);
    return ns.length === 3 && ns.every(n => n.atom.element === 'C' && n.bond.order === 1) &&
      !ns.some(n => nb(g, n.atom.id).some(m => m.atom.element === 'O' && m.bond.order === 2)) && !isAromatic(g, a.id);
  }
  function conjugatedTo(g, id, except) {           // 'aryl' | 'vinyl' | null — a π system next door
    for (const n of nb(g, id)) {
      if (n.atom.id === except || n.atom.element !== 'C') continue;
      if (isAromatic(g, n.atom.id)) return 'aryl';
      if (hasPiToCarbon(g, n.atom.id, id)) return 'vinyl';
    }
    return null;
  }
  function ewgNext(g, id, except) {                // a carbonyl / nitrile carbon next door (electron-withdrawing)
    return nb(g, id).some(n => n.atom.id !== except && n.atom.element === 'C' &&
      nb(g, n.atom.id).some(m => (m.atom.element === 'O' && m.bond.order === 2) || (m.atom.element === 'N' && m.bond.order === 3)));
  }
  function heteroNext(g, id, except) {             // O or N directly on this carbon (enol ether / enamine)
    return nb(g, id).some(n => n.atom.id !== except && (n.atom.element === 'O' || n.atom.element === 'N') && n.bond.order === 1);
  }
  function hindrance(g, id, except) {              // crowding around a carbon: branched neighbours
    let h = 0;
    nb(g, id).forEach(n => {
      if (n.atom.id === except || n.atom.element !== 'C') return;
      const d = deg(g, n.atom);
      if (d >= 4) h += 2; else if (d === 3) h += 1;
    });
    return h;
  }
  const clsName = d => d === 0 ? 'methyl' : d === 1 ? 'primary (1°)' : d === 2 ? 'secondary (2°)' : d === 3 ? 'tertiary (3°)' : 'quaternary';
  function nameOf(g) {
    try { const d = C().displayNameFor(g); if (d && !C().nameFailed(d) && d.length < 120 && !/ (only|when|because) /.test(d)) return d; } catch (e) {}
    try { const s = C().iupacName(g); if (s && !C().nameFailed(s) && s.length < 90 && !/ (only|when|because) /.test(s)) return s; } catch (e) {}
    try { const c = C().commonNameFor(g); if (c) return c; } catch (e) {}
    try { return C().formula(g); } catch (e) { return '?'; }
  }
  function keyOf(g) { try { return C().canonicalKey(g) + '|' + C().stereoSignature(g); } catch (e) { return Math.random() + ''; } }
  function molKey(g) { try { return C().canonicalKey(g); } catch (e) { return Math.random() + ''; } }
  const pct = x => (x > 0 && x < 0.005) ? '<1' : String(Math.round(x * 100));
  /* give a freshly added atom the id the PRODUCT uses for it, so the player can
     watch the same atom travel through the intermediate */
  function renameAtom(g, from, to) {
    if (from === to || g.atoms.some(a => a.id === to)) return;
    const a = atom(g, from); if (!a) return;
    a.id = to; g.bonds.forEach(b => { if (b.a === from) b.a = to; if (b.b === from) b.b = to; });
    if (g.nextId <= to) g.nextId = to + 1;
  }

  /* carbocation (or radical) stability at carbon X once it holds the + charge:
     a log10 scale where each class is worth about two powers of ten */
  function cationScore(g, X, partner) {
    const d = deg(g, X);                           // carbon neighbours in the SUBSTRATE (partner included)
    let s = d <= 1 ? 0 : d === 2 ? 2 : 4;          // 1° ~ 0, 2° ~ 2, 3° ~ 4
    const conj = conjugatedTo(g, X.id, partner ? partner.id : undefined);
    if (conj) s += 2.5;                            // allylic / benzylic resonance
    if (heteroNext(g, X.id, partner ? partner.id : undefined)) s += 3;   // oxocarbenium / iminium
    if (ewgNext(g, X.id, partner ? partner.id : undefined)) s -= 3;      // next to C=O: destabilised
    return s;
  }
  const catLabel = (g, X, partner) => {
    const d = deg(g, X);
    const conj = conjugatedTo(g, X.id, partner ? partner.id : undefined);
    return (d <= 1 ? 'primary' : d === 2 ? 'secondary' : 'tertiary') + (conj ? ' and ' + (conj === 'aryl' ? 'benzylic' : 'allylic') : '');
  };
  const share = (scores) => {                      // 10^score, normalised
    const w = scores.map(s => Math.pow(10, s));
    const t = w.reduce((a, b) => a + b, 0) || 1;
    return w.map(x => x / t);
  };

  /* ---------- stereo helpers ---------- */
  function mirror(g) { const m = clone(g); m.atoms.forEach(a => { a.x = -a.x; }); return m; }
  function centresOf(g) { return g.atoms.filter(a => a.element === 'C' && C().isStereocenter(g, a)).map(a => a.id); }
  function definedCentres(g) {
    return centresOf(g).filter(id => { const r = C().assignRS(g, id); return r && r.label; });
  }
  /* The two enantiomers of a product, when the reaction has just made a new
     stereocentre without any facial bias. Returns null when the product has
     no new centre, one form when it is meso. */
  function enantioForms(before, after, changed, how) {
    const newC = centresOf(after).filter(id => changed.has(id) && !(atom(before, id) && C().isStereocenter(before, atom(before, id))));
    if (!newC.length) return null;
    const chiralBefore = definedCentres(before).length > 0;
    const drawn = after.bonds.some(b => b.stereo && newC.includes(b.narrow));
    let A = clone(after);
    if (!drawn) {
      /* draw ONE configuration at each new centre (the mirror gives the other) */
      newC.forEach(id => { try { C().realizeStereo(A, id, () => { const r = C().assignRS(A, id); return !!(r && r.label); }); } catch (e) {} });
    }
    const B = mirror(A);
    let meso = false;
    try { meso = C().stereoSignature(A) === C().stereoSignature(B) && C().stereoSignature(A) !== ''; } catch (e) {}
    const lab = gg => { try { return newC.map(id => { const r = C().assignRS(gg, id); return r && r.label ? r.label : '?'; }).join(','); } catch (e) { return '?'; } };
    if (meso) return { kind: 'meso', centres: newC, forms: [{ graph: A, label: 'meso — one compound', share: 1 }],
      note: 'The two new stereocentres are mirror images of each other inside the same molecule (an internal mirror plane), so the “pair” is one meso compound, not two enantiomers.' };
    if (chiralBefore) return { kind: 'diastereomers', centres: newC,
      forms: [{ graph: A, label: '(' + lab(A) + ')', share: 0.5 }, { graph: mirrorNewOnly(A, newC), label: '(' + lab(mirrorNewOnly(A, newC)) + ')', share: 0.5 }],
      note: 'The starting material already has a stereocentre, so the two faces are NOT equivalent: the products are diastereomers, and the ratio depends on which face is less hindered (drawn 50:50 because that facial bias is not modelled).' };
    return { kind: 'racemic', centres: newC,
      forms: [{ graph: A, label: '(' + lab(A) + ')', share: 0.5 }, { graph: B, label: '(' + lab(B) + ')', share: 0.5 }],
      note: how || 'The new stereocentre is made through a flat intermediate or from an achiral starting material, so both faces are attacked equally: a 50:50 racemic mixture. One wedge alone would claim a single enantiomer, which is why both are drawn.' };
  }
  function mirrorNewOnly(g, ids) {                 // invert only the listed centres (swap wedge/dash they own)
    const m = clone(g);
    m.bonds.forEach(b => { if (b.stereo && ids.includes(b.narrow)) b.stereo = b.stereo === 'wedge' ? 'dash' : 'wedge'; });
    return m;
  }
  /* Draw a syn (same-face) relationship: the new substituent at c1 on a wedge,
     the new one at c2 on a wedge too. Ring/other marks are left alone. */
  function markSyn(g, c1, n1, c2, n2) {
    const b1 = bondOf(g, c1, n1), b2 = bondOf(g, c2, n2);
    if (b1) { b1.stereo = 'wedge'; b1.narrow = c1; }
    if (b2) { b2.stereo = 'wedge'; b2.narrow = c2; }
  }
  /* Put an ADDED hydrogen on the same face as a wedge substituent at the other
     carbon: since H is implicit, mark one existing substituent of c with a dash. */
  function markImplicitHSyn(g, c, avoid) {
    const cands = nb(g, c).filter(n => n.atom.id !== avoid && !n.bond.stereo);
    cands.sort((p, q) => (bondInRing(g, c, p.atom.id) ? 1 : 0) - (bondInRing(g, c, q.atom.id) ? 1 : 0));
    if (!cands.length) return false;
    cands[0].bond.stereo = 'dash'; cands[0].bond.narrow = c;
    return true;
  }
  /* SN2 inversion: rebuild the centre with the nucleophile where the leaving
     group was, on the OPPOSITE face. Works from the substrate's own parity. */
  function invertAt(gBefore, gAfter, cId, lgId, nuId) {
    try {
      const orderB = C().canonicalStereoOrder(gBefore, cId);
      if (orderB.length !== 4) return false;
      const parB = C().geometricParity(gBefore, cId, orderB);
      if (parB === null) return false;                       // substrate was flat: nothing defined to invert
      const orderA = orderB.map(o => (o.kind === 'atom' && o.id === lgId) ? { kind: 'atom', id: nuId } : o);
      gAfter.bonds.forEach(b => { if (b.narrow === cId) { delete b.stereo; delete b.narrow; } });
      return C().realizeStereo(gAfter, cId, () => C().geometricParity(gAfter, cId, orderA) === 1 - parB);
    } catch (e) { return false; }
  }
  function ezOf(g, a, b) { try { const r = C().assignEZ(g, bondOf(g, a, b)); return r ? r.label : null; } catch (e) { return null; } }
  function ezForms(g, a, b, eShare, why) {
    const cur = ezOf(g, a, b);
    if (!cur) return null;
    const other = clone(g); C().flipAlkeneEnd(other, a, b);
    const E = cur === 'E' ? g : other, Z = cur === 'E' ? other : g;
    return [{ graph: E, label: '(E) — trans', share: eShare, why: why || 'The E (trans) alkene keeps the two larger groups apart, so it is the more stable and the major geometric isomer.' },
            { graph: Z, label: '(Z) — cis', share: 1 - eShare, why: 'The Z (cis) alkene forms too, but the crowding between its groups makes it the minor isomer.' }];
  }

  /* ---------- SITE DESCRIPTORS ---------- */
  function alkeneSites(g) {
    const seenRing = new Set();
    let aromSys = null;
    return C().findGroups(g).filter(x => x.type === 'alkene').map(x => {
      const [A, B] = x.atoms;
      const sub = deg(g, A) - 1 + deg(g, B) - 1;          // alkyl groups on the C=C (0..4)
      const s = { kind: 'alkene', atoms: [A.id, B.id], A, B, bond: x.bond, sub,
        aromatic: !!x.aromatic, ring: bondInRing(g, A.id, B.id),
        ewg: ewgNext(g, A.id, B.id) || ewgNext(g, B.id, A.id),
        conj: conjugatedTo(g, A.id, B.id) || conjugatedTo(g, B.id, A.id),
        hind: hindrance(g, A.id, B.id) + hindrance(g, B.id, A.id),
        label: 'C' + A.id + '=C' + B.id };
      if (s.aromatic) {
        /* one site per aromatic ring system, not one per Kekulé bond */
        if (!aromSys) { try { aromSys = C().ringSystemsOf(g, g.atoms.filter(a => isAromatic(g, a.id)).map(a => a.id)); } catch (e) { aromSys = []; } }
        const sys = aromSys.find(r => r.includes(A.id));
        const key = sys ? sys.slice().sort((p, q) => p - q).join(',') : s.label;
        if (seenRing.has(key)) return null;
        seenRing.add(key);
        s.ringIds = sys || [A.id, B.id];
        s.label = 'aromatic ring (' + s.ringIds.length + ' atoms)';
      }
      return s;
    }).filter(Boolean);
  }
  function alkyneSites(g) {
    return C().findGroups(g).filter(x => x.type === 'alkyne').map(x => ({ kind: 'alkyne', atoms: [x.atoms[0].id, x.atoms[1].id], A: x.atoms[0], B: x.atoms[1], bond: x.bond,
      sub: deg(g, x.atoms[0]) - 1 + deg(g, x.atoms[1]) - 1, label: 'C' + x.atoms[0].id + '≡C' + x.atoms[1].id }));
  }
  function alcoholSites(g) {
    return C().findGroups(g).filter(x => x.type === 'alcohol').map(x => {
      const c = x.atoms[0], o = x.atoms[1];
      const d = deg(g, c);
      const phenol = isAromatic(g, c.id);
      const acid = !phenol && nb(g, c.id).some(n => n.atom.element === 'O' && n.bond.order === 2);
      const enol = !phenol && !acid && nb(g, c.id).some(n => n.bond.order > 1);
      return { kind: 'alcohol', atoms: [c.id, o.id], c, o, bond: x.bond, deg: d, cls: phenol ? 'phenol' : acid ? 'carboxylic acid O–H' : enol ? 'enol' : clsName(d),
        hind: hindrance(g, c.id), conj: conjugatedTo(g, c.id), phenol, acid, enol, hOnC: hOn(g, c), label: 'C' + c.id + '–OH' };
    });
  }
  function halideSites(g) {
    return C().findGroups(g).filter(x => x.type === 'halide').map(x => {
      const c = x.atoms[0], xa = x.atoms[1];
      const d = deg(g, c);
      const vinyl = nb(g, c.id).some(n => n.bond.order > 1), aryl = isAromatic(g, c.id);
      return { kind: 'halide', atoms: [c.id, xa.id], c, x: xa, bond: x.bond, deg: d, cls: aryl ? 'aryl' : vinyl ? 'vinyl' : clsName(d), el: xa.element,
        hind: hindrance(g, c.id), conj: conjugatedTo(g, c.id), betas: C().betaHydrogens(g, c.id), label: 'C' + c.id + '–' + xa.element };
    });
  }
  function carbonylSites(g) {
    const fgs = C().findFunctionalGroups(g);
    return C().findGroups(g).filter(x => x.type === 'carbonyl').map(x => {
      const c = x.atoms[0], o = x.atoms[1];
      const fg = fgs.find(f => f.atoms.includes(c.id) && f.atoms.includes(o.id));
      const kind = fg ? fg.name : 'carbonyl';
      const conjugated = nb(g, c.id).some(n => n.atom.element === 'C' && hasPiToCarbon(g, n.atom.id, c.id));
      return { kind: 'carbonyl', atoms: [c.id, o.id], c, o, bond: x.bond, cls: kind, hind: hindrance(g, c.id), conjugated, label: 'C' + c.id + '=O' };
    });
  }

  /* ---------- FAMILIES ---------- */
  const F = {};

  /* --- electrophilic addition of H–X or H–OH through a carbocation --- */
  F.addHX = {
    rules: ['hbr', 'hcl', 'hi', 'hydration'],
    needs: 'alkene', fallback: ['epoxOpen', 'etherCleave', 'esterHyd', 'alcHX', 'amideHyd'],
    sites(g, ruleId, env) {
      return alkeneSites(g).map(s => {
        if (s.aromatic) return Object.assign(s, { inert: 'aromatic C=C — adding across it would spend the ring’s ~150 kJ/mol of aromatic stabilisation, so benzene rings substitute rather than add' });
        const sA = cationScore(g, s.A, s.B), sB = cationScore(g, s.B, s.A);
        const best = Math.max(sA, sB);
        const X = sA >= sB ? s.A : s.B;
        s.score = best;
        s.reason = 'a ' + ['un', 'mono', 'di', 'tri', 'tetra'][Math.min(4, s.sub)] + 'substituted C=C — protonating it puts the + charge on a ' + catLabel(g, X, X === s.A ? s.B : s.A) + ' carbon' +
          (s.conj && !s.ewg ? ' (resonance-stabilised)' : '') +
          (s.ewg ? '; but the neighbouring C=O pulls electron density out of this π bond, which slows it badly' : '');
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const rad = ruleId === 'hbr' && env.light;   // hν switches HBr to the radical chain
      const Xel = ruleId === 'hydration' ? 'O' : ruleId === 'hbr' ? 'Br' : ruleId === 'hcl' ? 'Cl' : 'I';
      const sA = cationScore(g0, s.A, s.B), sB = cationScore(g0, s.B, s.A);
      const out = [];
      const build = (Xc, Hc, tag) => {
        const g = clone(g0); const b = bondOf(g, s.A.id, s.B.id); b.order = 1;
        const x = C().addAtomNear(g, Xel, atom(g, Xc.id));
        return { graph: g, changed: new Set([Xc.id, Hc.id, x.id]), added: [x.id] };
      };
      const dirs = [{ Xc: s.A, Hc: s.B, sc: sA }, { Xc: s.B, Hc: s.A, sc: sB }];
      dirs.sort((p, q) => q.sc - p.sc);
      const tie = Math.abs(dirs[0].sc - dirs[1].sc) < 1e-9;
      let fr = share(dirs.map(d => rad ? -d.sc : d.sc));     // radical: X goes to the carbon that leaves the MORE stable radical
      if (rad) { dirs.reverse(); fr = [fr[1], fr[0]]; }
      const xName = ruleId === 'hydration' ? 'OH' : Xel;
      dirs.forEach((d, i) => {
        const p = build(d.Xc, d.Hc);
        const lab = tie ? (i === 0 ? xName + ' on C' + d.Xc.id : xName + ' on C' + d.Xc.id + ' (equivalent)')
          : (i === 0 ? (rad ? 'anti-Markovnikov' : 'Markovnikov') + ' product' : (rad ? 'Markovnikov' : 'anti-Markovnikov') + ' product');
        p.label = lab; p.share = fr[i];
        p.why = tie ? 'Both alkene carbons carry the same number of hydrogens, so neither carbocation is better — the two orientations are equally likely.'
          : rad ? (i === 0 ? 'Light splits H–Br into radicals. Br• adds to the LESS substituted carbon so the unpaired electron lands on the ' + catLabel(g0, d.Hc, d.Xc) + ' carbon — the more stable radical. That reverses Markovnikov.' : 'Adding Br• to the other carbon would leave a less stable radical; only a trace goes this way.')
          : (i === 0 ? 'H⁺ adds to the carbon that already has more hydrogens, so the + charge lands on the ' + catLabel(g0, d.Xc, d.Hc) + ' carbon — the more stable carbocation — and ' + xName + ' finishes there (Markovnikov).'
                     : 'The other orientation would need a ' + catLabel(g0, d.Xc, d.Hc) + ' carbocation, which is far less stable — essentially none forms.');
        p.stereoHow = rad ? 'The carbon radical is essentially planar, so bromide is delivered to either face — racemic at that centre.' : 'The carbocation is flat (sp²), so the nucleophile attacks either face equally — a new stereocentre comes out racemic.';
        out.push(p);
      });
      /* conjugated diene: the first cation is allylic, so X can be captured at
         either end of the allyl system — 1,2- (kinetic) or 1,4- (thermodynamic) */
      if (!rad && ruleId !== 'hydration') {
        const main = dirs[0]; const Xc = main.Xc;
        const partner = nb(g0, Xc.id).find(n => n.atom.element === 'C' && n.bond.order === 1 && n.atom.id !== main.Hc.id &&
          nb(g0, n.atom.id).some(m => m.bond.order === 2 && m.atom.element === 'C' && m.atom.id !== Xc.id && !isAromatic(g0, m.atom.id)));
        if (partner) {
          const Cc = partner.atom;
          const Dd = nb(g0, Cc.id).find(m => m.bond.order === 2 && m.atom.element === 'C' && m.atom.id !== Xc.id).atom;
          const g = clone(g0); bondOf(g, s.A.id, s.B.id).order = 1; bondOf(g, Cc.id, Dd.id).order = 1; bondOf(g, Xc.id, Cc.id).order = 2;
          const x = C().addAtomNear(g, Xel, atom(g, Dd.id));
          const T = env.T || 298;
          const f12 = Math.max(0.15, Math.min(0.85, 0.8 - 0.6 * (T - 193) / 120));   // 80:20 at −80 °C → 20:80 at 40 °C
          const p12 = out[0]; const tot = p12.share;
          p12.share = tot * f12;
          p12.label = '1,2-addition product (kinetic)';
          p12.allyl = { Xc: Xc.id, Cc: Cc.id, Dd: Dd.id, mode: '1,2' };
          p12.why = 'H⁺ adds to the end carbon C' + main.Hc.id + ' of the conjugated diene, so the + charge lands on C' + Xc.id + ' next to the other C=C: an allylic carbocation, stabilised by resonance over C' + Xc.id + ' and C' + Dd.id + '. ' + xName + '⁻ captured at C' + Xc.id + ' (right next to the new C–H) gives the 1,2-adduct. It forms FASTER — the kinetic product — so it dominates cold (≈ 80% at −80 °C).';
          const p14 = { graph: g, changed: new Set([Xc.id, main.Hc.id, Cc.id, Dd.id, x.id]), added: [x.id], share: tot * (1 - f12),
            label: '1,4-addition product (thermodynamic)', allyl: { Xc: Xc.id, Cc: Cc.id, Dd: Dd.id, mode: '1,4' },
            why: 'The same allylic carbocation, captured by ' + xName + '⁻ at its far end C' + Dd.id + ': the π bond shifts to C' + Xc.id + '=C' + Cc.id + ' and the halogen lands four atoms from the new H — 1,4- (conjugate) addition. The double bond ends up more substituted, so this isomer is MORE STABLE — the thermodynamic product. Warm, the adducts equilibrate (the C–X bond re-ionises to the same cation) and it takes over (≈ 80% at 40 °C); the new C=C is mostly E.',
            stereoHow: 'The allyl cation is planar; ' + xName + '⁻ adds to either face, and the new C=C forms mainly as the E (trans) isomer.' };
          out.splice(1, 0, p14);
        }
      }
      /* rearrangement: hydride or methyl shift to a better cation, then X there */
      if (!rad && !(out[0].allyl)) {
        const main = dirs[0];
        const Xc = main.Xc;
        const shifts = [];
        nb(g0, Xc.id).forEach(n => {
          if (n.atom.id === main.Hc.id || n.atom.element !== 'C' || n.bond.order !== 1) return;
          const N = n.atom;
          /* after the shift N holds the charge: its degree counts Xc as a neighbour still */
          const gain = (deg(g0, N) <= 1 ? 0 : deg(g0, N) === 2 ? 2 : 4) + (conjugatedTo(g0, N.id, Xc.id) ? 2.5 : 0);
          const cur = cationScore(g0, Xc, main.Hc);
          if (hOn(g0, N) >= 1 && gain > cur + 0.5) shifts.push({ N, kind: 'hydride', gain });
          else if (hOn(g0, N) === 0 && deg(g0, N) === 4 && nb(g0, N.id).some(m => m.atom.element === 'C' && deg(g0, m.atom) === 1 && hOn(g0, m.atom) === 3) && gain - 1 > cur + 0.5)
            shifts.push({ N, kind: 'methyl', gain: gain - 1 });
        });
        if (shifts.length) {
          shifts.sort((p, q) => q.gain - p.gain);
          const sh = shifts[0];
          const g = clone(g0); bondOf(g, s.A.id, s.B.id).order = 1;
          let target = sh.N.id;
          if (sh.kind === 'methyl') {
            /* a methyl walks from N to Xc; the + (and then X) stays on N */
            const me = nb(g, sh.N.id).find(m => m.atom.element === 'C' && deg(g, m.atom) === 1 && hOn(g, m.atom) === 3);
            const mb = bondOf(g, sh.N.id, me.atom.id); if (mb.a === sh.N.id) mb.a = Xc.id; else mb.b = Xc.id;
          }
          const x = C().addAtomNear(g, Xel, atom(g, target));
          const gLab = sh.gain >= 4 ? 'tertiary' : 'secondary';
          const p = { graph: g, changed: new Set([Xc.id, target, x.id]), added: [x.id],
            label: 'rearranged product (' + (sh.kind === 'hydride' ? '1,2-hydride' : '1,2-methyl') + ' shift)',
            why: 'The first-formed ' + catLabel(g0, Xc, main.Hc) + ' carbocation sits next to a carbon that can hand over ' + (sh.kind === 'hydride' ? 'a hydride' : 'a methyl group') +
              '. That 1,2-shift moves the + charge to a ' + gLab + (conjugatedTo(g0, sh.N, Xc.id) ? ', resonance-stabilised' : '') + ' carbon — a more stable cation — and ' + xName + ' is captured THERE. With HCl on 3-methylbut-1-ene this rearranged halide is about 60% of the product, the unrearranged one about 40%.',
            stereoHow: 'The rearranged carbocation is flat too, so this product is also racemic at the new centre.' };
          /* the rearranged product takes ~60% of what the main orientation would have made */
          const mainOut = out[0];
          p.share = mainOut.share * 0.6; mainOut.share *= 0.4;
          mainOut.why += ' A 1,2-shift competes here (see the rearranged product), so this direct product is the smaller share.';
          out.unshift(p);
        }
      }
      return out;
    },
    gate(g, ruleId, env) {
      if (ruleId === 'hydration' && env.pHVal > 2.5)
        return 'Acid-catalysed hydration needs a genuinely acidic medium (H₃O⁺ from H₂SO₄, pH ≲ 1). At pH ' + env.pHVal.toFixed(1) + ' there is too little H₃O⁺ to protonate the alkene, so nothing happens.';
      if (ruleId !== 'hydration' && env.pHVal > 3)
        return 'H–' + (ruleId === 'hbr' ? 'Br' : ruleId === 'hcl' ? 'Cl' : 'I') + ' IS the acid here. Buffering the medium to pH ' + env.pHVal.toFixed(1) + ' neutralises it to the halide salt, which does not add to alkenes.';
      return null;
    },
    condNotes(ruleId, env) {
      const n = [];
      if (ruleId !== 'hydration' && env.T !== undefined && env.T < 233) n.push('Very cold (' + degC(env.T) + ' °C): with a conjugated diene the 1,2-adduct (kinetic) is trapped — nothing re-ionises, so the faster-formed product is what you isolate.');
      if (ruleId !== 'hydration' && env.T > 303) n.push('Warm (' + degC(env.T) + ' °C): with a conjugated diene the allylic C–X bond re-ionises, the adducts equilibrate and the 1,4-adduct (thermodynamic, more substituted C=C) wins.');
      if (ruleId === 'hbr' && env.light) n.push('hν is on: light splits H–Br homolytically, so the radical chain (anti-Markovnikov) runs instead of the ionic addition. This switch works for HBr only — for HCl the H-abstraction step is uphill, for HI the I• addition is.');
      if (ruleId === 'hcl' && env.light) n.push('hν has no effect on HCl: the chain would need Cl–H to give up its hydrogen to a carbon radical, which is endothermic, so the ionic Markovnikov path stays in charge.');
      if (env.T > 373) n.push('At ' + degC(env.T) + ' °C the addition is reversible enough that some alkene is re-formed by E1; the equilibrium still favours the adduct.');
      if (env.T < 273) n.push('Cold (' + degC(env.T) + ' °C) slows the reaction but sharpens the Markovnikov preference and discourages rearrangement.');
      return n;
    },
    byproducts: () => [],
    stages(g0, s, out, ruleId, env) {
      if (ruleId === 'hbr' && env.light) return null;
      const Xc = [...out.changed].map(id => atom(g0, id)).filter(a => a && a.element === 'C');
      const xAtom = atom(out.graph, out.added[0]);
      const xC = nb(out.graph, xAtom.id)[0].atom;
      const g = clone(g0); bondOf(g, s.A.id, s.B.id).order = 1;
      if (out.allyl) {
        const al = out.allyl;
        if (al.mode === '1,4') { bondOf(g, al.Cc, al.Dd).order = 1; bondOf(g, al.Xc, al.Cc).order = 2; atom(g, al.Dd).charge = 1; }
        else atom(g, al.Xc).charge = 1;
        return [{ label: 'allylic carbocation' + (al.mode === '1,4' ? ' (resonance form, + on C' + al.Dd + ')' : ' (+ on C' + al.Xc + ', shared with C' + al.Dd + ')'), graph: g }];
      }
      const cc = atom(g, xC.id); if (cc) cc.charge = 1;
      return [{ label: 'carbocation intermediate', graph: g }];
    },
  };

  /* --- radical anti-Markovnikov HBr --- */
  F.addRad = {
    rules: ['hbr_peroxide'], needs: 'alkene',
    sites(g, ruleId, env) {
      return alkeneSites(g).map(s => {
        if (s.aromatic) return Object.assign(s, { inert: 'aromatic C=C — a bromine radical does not add to a benzene ring; the aromatic stabilisation is worth far more than the new C–Br bond' });
        const sA = cationScore(g, s.A, s.B), sB = cationScore(g, s.B, s.A);
        s.score = Math.max(sA, sB) * 0.5;            // radicals care about substitution too, but less steeply
        s.reason = 'a ' + ['un', 'mono', 'di', 'tri', 'tetra'][Math.min(4, s.sub)] + 'substituted C=C — Br• adds so that the unpaired electron lands on a ' + catLabel(g, sA >= sB ? s.A : s.B, sA >= sB ? s.B : s.A) + ' carbon';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const e2 = Object.assign({}, env, { light: true });
      return F.addHX.variants(g0, s, 'hbr', e2);
    },
    gate: () => null,
    condNotes(ruleId, env) {
      const n = ['Peroxides (ROOR) — or light — split H–Br into radicals; the chain then adds Br• first, which is why the regiochemistry flips. Without an initiator this reagent is plain HBr (Markovnikov).'];
      if (env.T > 353) n.push('Heat helps the peroxide fall apart (its O–O bond is only ~150 kJ/mol), so the chain starts faster at ' + degC(env.T) + ' °C.');
      return n;
    },
    byproducts: () => [],
    stages(g0, s, out) {
      const xAtom = atom(out.graph, out.added[0]); const xC = nb(out.graph, xAtom.id)[0].atom;
      const other = xC.id === s.A.id ? s.B : s.A;
      const g = clone(g0); bondOf(g, s.A.id, s.B.id).order = 1;
      const br = C().addAtomNear(g, 'Br', atom(g, xC.id)); renameAtom(g, br.id, out.added[0]);
      const r = atom(g, other.id); if (r) r.radical = true;
      return [{ label: 'carbon radical intermediate', graph: g }];
    },
  };

  /* --- Br2 / Cl2 through the bridged halonium ion --- */
  F.addX2 = {
    rules: ['br2', 'cl2'], needs: 'alkene', fallback: ['alphaHal', 'hvz'],
    sites(g, ruleId, env) {
      return alkeneSites(g).map(s => {
        if (s.aromatic) return Object.assign(s, { inert: 'aromatic C=C — Br₂ alone does not add to a benzene ring (that needs FeBr₃ and gives SUBSTITUTION)' });
        s.score = 1.8 * s.sub - (s.ewg ? 3 : 0) + (s.conj && !s.ewg ? 0.8 : 0);
        s.reason = 'a ' + ['un', 'mono', 'di', 'tri', 'tetra'][Math.min(4, s.sub)] + 'substituted alkene — each alkyl group makes the π bond richer in electrons and roughly 50–100× faster at attacking Br₂ (ethene 1 : propene 61 : 2-methylpropene 5,400 : 2,3-dimethylbut-2-ene 920,000, bromination in methanol)' +
          (s.ewg ? '; the C=O next door pulls electrons out of it, which slows it badly' : '');
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const X = ruleId === 'br2' ? 'Br' : 'Cl';
      const out = [];
      const water = env.solventVal === 'water';
      const build = (nuc2) => {
        const g = clone(g0); bondOf(g, s.A.id, s.B.id).order = 1;
        const sA = cationScore(g0, s.A, s.B), sB = cationScore(g0, s.B, s.A);
        const more = sA >= sB ? s.A : s.B, less = more === s.A ? s.B : s.A;     // the carbon that carries more + in the bridged ion
        const c1 = atom(g, less.id), c2 = atom(g, more.id);
        const x1 = C().addAtomNear(g, X, c1);                     // halogen from the bridge on one carbon…
        const x2 = C().addAtomNear(g, nuc2, c2);                  // …nucleophile opens the bridge at the other, from the back
        const b1 = bondOf(g, c1.id, x1.id), b2 = bondOf(g, c2.id, x2.id);
        if (b1) { b1.stereo = 'wedge'; b1.narrow = c1.id; }
        if (b2) { b2.stereo = 'dash'; b2.narrow = c2.id; }
        return { graph: g, changed: new Set([c1.id, c2.id, x1.id, x2.id]), added: [x1.id, x2.id] };
      };
      const di = build(X);
      di.label = 'vicinal di' + (X === 'Br' ? 'bromide' : 'chloride') + ' (anti)';
      di.why = 'The π bond attacks ' + X + '₂ and a cyclic ' + (X === 'Br' ? 'bromonium' : 'chloronium') + ' ion bridges both carbons. ' + X + '⁻ can only reach the face opposite the bridge, so the two halogens end up ANTI (one wedge, one dash).';
      di.stereoHow = 'Anti addition fixes the RELATIVE configuration; with no chiral influence both mirror-image anti products form equally.';
      if (water) {
        const hh = build('O');
        hh.label = (X === 'Br' ? 'bromohydrin' : 'chlorohydrin') + ' (anti, Markovnikov OH)';
        hh.why = 'In water the solvent is present in huge excess, so H₂O — not ' + X + '⁻ — opens the ' + (X === 'Br' ? 'bromonium' : 'chloronium') + ' ion. It attacks the carbon that carries more of the positive charge (the more substituted one), so OH lands there, anti to the halogen.';
        hh.stereoHow = di.stereoHow;
        hh.share = 0.85; di.share = 0.15;
        di.why += ' In water this is the minor path — water out-competes the halide ion for the bridged intermediate.';
        out.push(hh, di);
      } else { di.share = 1; out.push(di); }
      return out;
    },
    gate: () => null,
    condNotes(ruleId, env) {
      const n = [];
      if (env.solventVal === 'water') n.push('Water as solvent: the halonium ion is opened by H₂O, so the halohydrin is the main product rather than the dihalide.');
      if (env.light) n.push('hν with ' + (ruleId === 'br2' ? 'Br₂' : 'Cl₂') + ' also starts a radical chain — with an alkane C–H available that is the substitution shown in Mechanism & Energy; here the ionic addition to the π bond is the faster process.');
      return n;
    },
    byproducts: (ruleId, env) => env.solventVal === 'water' ? [{ smiles: ruleId === 'br2' ? 'Br' : 'Cl', name: 'H' + (ruleId === 'br2' ? 'Br' : 'Cl') }] : [],
    stages(g0, s, out, ruleId) {
      const X = ruleId === 'br2' ? 'Br' : 'Cl';
      const g = clone(g0); bondOf(g, s.A.id, s.B.id).order = 1;
      const br = C().addAtomNear(g, X, atom(g, s.A.id)); br.charge = 1;
      g.bonds.push({ a: br.id, b: s.B.id, order: 1 });
      if (out && out.added && out.added.length) renameAtom(g, br.id, out.added[0]);
      return [{ label: (X === 'Br' ? 'bromonium' : 'chloronium') + ' ion — bridged intermediate', graph: g }];
    },
  };

  /* --- concerted syn additions: mCPBA, OsO4, hydroboration, H2 --- */
  F.epox = {
    rules: ['mcpba'], needs: 'alkene',
    sites(g) {
      /* a sulfide first, then a basic nitrogen, then the alkene — the order of nucleophilicity
         towards the peracid's electrophilic oxygen, which is what actually decides it */
      const extra = [];
      g.atoms.forEach(a => {
        if (sulfideAt(g, a)) extra.push({ kind: 'sulfide', s: a, score: 3.0, label: 'sulfide (S' + a.id + ')',
          reason: 'the sulfide — sulfur’s lone pair is a better nucleophile than any alkene here, so it is oxidised first, and one equivalent stops at the sulfoxide' });
        else if (pyridineN(g, a)) extra.push({ kind: 'noxide', n: a, pyridine: true, score: 2.0, label: 'ring nitrogen (N' + a.id + ')',
          reason: 'a pyridine-type ring nitrogen — oxidised to the N-oxide, which is usually the point: the N-oxide nitrates at C4 and opens C2 to nucleophiles' });
        else if (tertAmine(g, a)) extra.push({ kind: 'noxide', n: a, pyridine: false, score: 1.7, label: 'tertiary amine (N' + a.id + ')',
          reason: 'a tertiary amine — oxidised to the N-oxide' });
      });
      return extra.concat(alkeneSites(g).map(s => {
        if (s.aromatic) return Object.assign(s, { inert: 'aromatic C=C — peracids do not epoxidise benzene rings' });
        s.score = 1.4 * s.sub - (s.ewg ? 4 : 0);
        s.reason = 'a ' + ['un', 'mono', 'di', 'tri', 'tetra'][Math.min(4, s.sub)] + 'substituted alkene — peracid epoxidation is electrophilic, and each alkyl group speeds it ~20–30× (ethene 1 : monosubstituted ~24 : trisubstituted ~6,500)' +
          (s.ewg ? '; conjugation with C=O makes this π bond electron-poor, and peracids leave those alone (basic H₂O₂ is the reagent for enones)' : '');
        return s;
      }));
    },
    variants(g0, s) {
      if (s.kind === 'sulfide' || s.kind === 'noxide') {
        const g = clone(g0);
        const at = s.kind === 'sulfide' ? s.s.id : s.n.id;
        const o = C().addAtomNear(g, 'O', atom(g, at));    // this already bonds O to `at`, order 1
        if (s.kind === 'sulfide') {
          const bd = bondOf(g, at, o.id); if (bd) bd.order = 2;
          return [{ graph: g, changed: new Set([at, o.id]), added: [o.id], share: 1, label: 'sulfoxide',
            why: 'Sulfur’s lone pair takes the peracid’s electrophilic outer oxygen — a better nucleophile than any alkene in the molecule, which is why the sulfide goes first. One equivalent stops here: the sulfoxide sulfur is already electron-poor, so the sulfone needs a second equivalent, a stronger oxidant, or heat.',
            stereoHow: 'The sulfur becomes a stereocentre (two different carbons, an oxygen, a lone pair). A plain peracid attacks both lone pairs equally, so the sulfoxide is racemic — one enantiomer needs a chiral titanium or vanadium catalyst with a hydroperoxide.' }];
        }
        const nA = atom(g, at); if (nA) nA.charge = 1;
        const oA = atom(g, o.id); if (oA) oA.charge = -1;
        return [{ graph: g, changed: new Set([at, o.id]), added: [o.id], share: 1, label: s.pyridine ? 'pyridine N-oxide' : 'amine N-oxide',
          why: 'The nitrogen lone pair takes the outer oxygen of the peracid, giving a genuine N⁺–O⁻ dipole. ' +
            (s.pyridine ? 'For a pyridine this is usually a means and not an end: the N-oxide is nitrated far more easily than pyridine itself (and at C4), and C2 becomes open to nucleophiles — which is how most 2-substituted pyridines are made. PCl₃, P(OEt)₃ or H₂/Pd removes the oxide again afterwards.'
                        : 'On heating, a tertiary amine oxide with a β-hydrogen does the Cope elimination instead, giving the alkene and a hydroxylamine.'),
          stereoHow: null }];
      }
      const g = clone(g0); bondOf(g, s.A.id, s.B.id).order = 1;
      const c1 = atom(g, s.A.id), c2 = atom(g, s.B.id);
      const o = C().addAtomNear(g, 'O', c1); g.bonds.push({ a: o.id, b: c2.id, order: 1 });
      const b1 = bondOf(g, c1.id, o.id), b2 = bondOf(g, c2.id, o.id);
      if (b1) { b1.stereo = 'wedge'; b1.narrow = c1.id; }
      if (b2) { b2.stereo = 'wedge'; b2.narrow = c2.id; }
      return [{ graph: g, changed: new Set([c1.id, c2.id, o.id]), added: [o.id], share: 1, label: 'epoxide (syn)',
        why: 'One concerted step: the peracid’s outer oxygen is delivered to one face of the π bond while its O–H proton swings back to the carbonyl. Both new C–O bonds are made on the SAME face, so whatever was cis on the alkene stays cis on the epoxide.',
        stereoHow: 'Syn and stereospecific; from an achiral alkene both faces are attacked equally, so the epoxide is racemic.' }];
    },
    gate: () => null,
    condNotes(ruleId, env) { const n = []; if (env.T > 313) n.push('Peracids are run cold (0 °C to room temperature); at ' + degC(env.T) + ' °C mCPBA starts to decompose and, with a ketone present, Baeyer–Villiger oxidation competes.'); return n; },
    byproducts: () => [{ smiles: 'OC(=O)c1cccc(Cl)c1', name: '3-chlorobenzoic acid (the spent peracid)' }],
    stages: () => null,
  };
  F.diol = {
    rules: ['oso4'], needs: 'alkene',
    sites(g) {
      return alkeneSites(g).map(s => {
        if (s.aromatic) return Object.assign(s, { inert: 'aromatic C=C — OsO₄ leaves benzene rings alone' });
        s.score = 0.6 * s.sub - (s.ewg ? 2 : 0) - (s.sub === 4 ? 1.5 : 0) - 0.4 * s.hind;
        s.reason = 'a ' + ['un', 'mono', 'di', 'tri', 'tetra'][Math.min(4, s.sub)] + 'substituted alkene — OsO₄ is electrophilic (electron-rich alkenes react faster) but bulky, so a crowded tetrasubstituted C=C is slow' +
          (s.ewg ? '; the C=O next door deactivates it' : '');
        return s;
      });
    },
    variants(g0, s) {
      const g = clone(g0); bondOf(g, s.A.id, s.B.id).order = 1;
      const c1 = atom(g, s.A.id), c2 = atom(g, s.B.id);
      const o1 = C().addAtomNear(g, 'O', c1), o2 = C().addAtomNear(g, 'O', c2);
      markSyn(g, c1.id, o1.id, c2.id, o2.id);
      return [{ graph: g, changed: new Set([c1.id, c2.id, o1.id, o2.id]), added: [o1.id, o2.id], share: 1, label: 'syn-1,2-diol',
        why: 'OsO₄ adds to one face of the alkene in a single concerted [3+2] step to give a cyclic osmate ester; hydrolysis then releases the diol with BOTH oxygens on the same face (syn). NMO re-oxidises the osmium so only a catalytic amount is needed.',
        stereoHow: 'Syn addition fixes the relative configuration; from an achiral alkene the two mirror-image syn diols form equally.' }];
    },
    gate: () => null,
    condNotes(ruleId, env) { const n = []; if (env.T > 323) n.push('Dihydroxylation is normally run at 0–25 °C; heating with excess oxidant risks over-oxidation to the diketone / cleavage.'); return n; },
    byproducts: () => [{ smiles: 'CN1CCOCC1', name: 'N-methylmorpholine (from NMO)' }],
    stages: () => null,
  };
  F.hydrobor = {
    rules: ['hydroboration'], needs: 'alkene',
    sites(g) {
      return alkeneSites(g).map(s => {
        if (s.aromatic) return Object.assign(s, { inert: 'aromatic C=C — borane does not add to benzene rings' });
        s.score = -1.2 * s.sub - 0.5 * s.hind;
        s.reason = 'a ' + ['un', 'mono', 'di', 'tri', 'tetra'][Math.min(4, s.sub)] + 'substituted alkene — hydroboration is controlled by STERICS, so the least crowded C=C reacts first (terminal > 1,2-disubstituted ≈ 1,1-disubstituted > trisubstituted ≫ tetrasubstituted)';
        return s;
      });
    },
    variants(g0, s) {
      const dA = deg(g0, s.A), dB = deg(g0, s.B);
      const hA = hindrance(g0, s.A.id, s.B.id), hB = hindrance(g0, s.B.id, s.A.id);
      /* boron (then OH) goes to the LESS substituted / less hindered carbon */
      const less = (dA + 0.3 * hA) <= (dB + 0.3 * hB) ? s.A : s.B, more = less === s.A ? s.B : s.A;
      const tie = dA === dB && hA === hB;
      const diff = Math.abs(deg(g0, more) - deg(g0, less));
      const styr = conjugatedTo(g0, more.id, less.id) === 'aryl' && deg(g0, less) === 1;
      let frac = tie ? 0.5 : styr ? 0.80 : diff >= 2 ? 0.99 : (deg(g0, more) === 2 ? 0.94 : 0.98);
      const build = (OHc, Hc) => {
        const g = clone(g0); bondOf(g, s.A.id, s.B.id).order = 1;
        const o = C().addAtomNear(g, 'O', atom(g, OHc.id));
        const bo = bondOf(g, OHc.id, o.id); if (bo) { bo.stereo = 'wedge'; bo.narrow = OHc.id; }
        /* the H added at the other carbon is on the same face as the OH */
        if (C().isStereocenter(g, atom(g, Hc.id))) markImplicitHSyn(g, Hc.id, OHc.id);
        return { graph: g, changed: new Set([OHc.id, Hc.id, o.id]), added: [o.id] };
      };
      const p1 = build(less, more), p2 = build(more, less);
      p1.share = frac; p2.share = 1 - frac;
      p1.label = tie ? 'OH on C' + less.id : 'anti-Markovnikov alcohol';
      p2.label = tie ? 'OH on C' + more.id + ' (equivalent)' : 'Markovnikov alcohol';
      p1.why = tie ? 'Both ends of the C=C are equally substituted, so boron has no preference.' :
        'Boron is the electrophile and it is bulky: it bonds to the LESS hindered carbon, and the hydride goes to the more substituted one, in one four-centre step (no carbocation, no rearrangement). The oxidation then swaps B for OH with retention. ' +
        (styr ? 'On a styrene the benzylic position competes (about 80:20, H. C. Brown).' : 'H. C. Brown measured ' + pct(frac) + ':' + pct(1 - frac) + ' for this substitution pattern with BH₃·THF; 9-BBN pushes it to >99:1.');
      p2.why = tie ? 'The equivalent orientation.' : 'Putting boron on the more crowded carbon is the minor path — sterics and the small partial charge both work against it.';
      p1.stereoHow = p2.stereoHow = 'Syn addition: H and OH end up on the same face. From an achiral alkene both faces react equally, so a new stereocentre is racemic.';
      return [p1, p2];
    },
    gate: (g, ruleId, env) => env.solventVal === 'water' ? 'BH₃ is destroyed by water (it hydrolyses to boric acid and H₂) — hydroboration needs dry THF; the H₂O₂/NaOH is only added afterwards.' : null,
    condNotes(ruleId, env) { return env.T > 323 ? ['Hydroboration is fast at 0–25 °C; heating is unnecessary and lets the alkylborane isomerise along the chain.'] : []; },
    byproducts: () => [{ smiles: 'OB(O)O', name: 'boric acid (after the H₂O₂ / NaOH step)' }],
    stages: () => null,
  };
  F.h2 = {
    rules: ['h2', 'lindlar', 'na_nh3', 'nabh3cn'], needs: 'alkene-or-alkyne', structural: true,
    sites(g, ruleId, env) {
      const hard = env.P >= 50 && env.T >= 373;
      /* the two half-reductions of an alkyne care only about alkynes */
      if (ruleId === 'lindlar' || ruleId === 'na_nh3') {
        return alkyneSites(g).map(x => { x.score = 1; x.reason = 'an alkyne — this reagent stops at the alkene, and which alkene it gives is decided by the mechanism, not by the substrate'; return x; })
          .concat(alkeneSites(g).map(x => Object.assign(x, { inert: ruleId === 'lindlar' ? 'an alkene — the Lindlar catalyst is deliberately poisoned (Pb(OAc)₂ / quinoline) so that it is too weak to touch a C=C. That is the entire point of it' : 'an alkene — sodium in ammonia reduces only a triple bond (or an aromatic ring, in a Birch reduction); an isolated C=C has no orbital low enough to take the electron' })));
      }
      /* C=N is reducible too, and that is the second half of every reductive amination.
         NaBH3CN is the reagent that does ONLY this: mild enough to leave the ketone alone at
         pH 6-7 and still reduce the iminium, which is what makes a one-pot procedure work. */
      const imines = [];
      g.atoms.forEach(a => {
        if (a.element !== 'C') return;
        const dn = nb(g, a.id).find(n => n.atom.element === 'N' && n.bond.order === 2);
        if (!dn) return;
        imines.push({ kind: 'imine', c: a, n: dn.atom, score: ruleId === 'nabh3cn' ? 3 : 0.5,
          label: 'imine C=N (C' + a.id + ')',
          reason: 'the C=N — reduced to the amine' + (ruleId === 'nabh3cn' ? '; sodium cyanoborohydride is mild and tolerates a mildly acidic pH, so it reduces the iminium and leaves the ketone it came from alone' : ', on the catalyst surface, just as a C=C is') });
      });
      if (ruleId === 'nabh3cn') {
        return imines.concat(alkeneSites(g).map(x => Object.assign(x, { inert: 'an alkene — NaBH₃CN is a hydride donor, not a catalyst: it reduces only the polarised C=N (and a C=O slowly)' })));
      }
      const s1 = imines.concat(alkeneSites(g).map(s => {
        if (s.aromatic) {
          if (hard) { s.score = -3; s.reason = 'an aromatic ring — at ' + env.P + ' atm and ' + degC(env.T) + ' °C it can be hydrogenated, but all its π bonds go at once, after any ordinary C=C'; s.aromaticRing = true; return s; }
          return Object.assign(s, { inert: 'an aromatic ring — needs ~50–100 atm H₂ and heat (Rh or Ru is better than Pd); at ' + env.P + ' atm and ' + degC(env.T) + ' °C it survives. Set pressure ≥ 50 atm and T ≥ 100 °C to see it reduced' });
        }
        s.score = -1.0 * s.sub - 0.5 * s.hind;
        s.reason = 'a ' + ['un', 'mono', 'di', 'tri', 'tetra'][Math.min(4, s.sub)] + 'substituted alkene — the alkene must lie flat on the metal, so the least substituted, least hindered C=C is hydrogenated fastest';
        return s;
      }));
      const s2 = alkyneSites(g).map(s => { s.score = 1.0 - 0.6 * s.sub; s.reason = 'an alkyne — it binds the metal surface more strongly than an alkene, so it is reduced first (and, with enough H₂, all the way to the alkane)'; return s; });
      return s1.concat(s2);
    },
    variants(g0, s, ruleId, env) {
      /* the two half-reductions of an alkyne: same product skeleton, opposite geometry, and
         the reason for each is the mechanism rather than the substrate */
      if (ruleId === 'lindlar' || ruleId === 'na_nh3') {
        const g = clone(g0); const b = bondOf(g, s.A.id, s.B.id);
        if (!b || b.order !== 3) return [];
        b.order = 2;
        try { C().layoutGraph(g); } catch (e) {}
        const want = ruleId === 'lindlar' ? 'Z' : 'E';
        try {
          const bb = C().bondBetween(g, s.A.id, s.B.id);
          if (bb) delete bb.ezUnspec;
          const lab = C().assignEZ(g, bb);
          if (lab && lab.label && lab.label !== want) C().flipAlkeneEnd(g, s.A.id, s.B.id);
        } catch (e) {}
        const cis = ruleId === 'lindlar';
        return [{ graph: g, changed: new Set([s.A.id, s.B.id]), added: [], share: 1,
          label: cis ? 'cis (Z) alkene — Lindlar' : 'trans (E) alkene — dissolving metal',
          why: cis
            ? 'The Lindlar catalyst is palladium on calcium carbonate, deliberately poisoned with lead(II) acetate and quinoline. Poisoning it does two things: it weakens the surface enough that the alkene formed no longer sticks (so the reduction stops after one addition), while the alkyne — which binds far more strongly — still adsorbs. Both hydrogens are delivered to the face lying on the metal, so the product is the CIS alkene. Stopping at the alkene is the whole purpose; on ordinary Pd–C the same substrate goes straight to the alkane.'
            : 'Sodium (or lithium) dissolved in liquid ammonia at −78 °C is a solution of solvated electrons. One electron goes into the alkyne π* to give a radical anion, which is protonated by ammonia to a vinyl radical; the vinyl radical adopts the geometry with its two substituents apart (the trans radical is lower in energy and the two invert rapidly, so the population is what matters); a second electron gives the vinyl anion, which keeps that geometry because sp² carbanions invert slowly, and a second ammonia protonates it. The result is the TRANS alkene — the complement to Lindlar, and the reason both reagents are taught together.',
          stereoHow: cis ? 'Both hydrogens arrive on the same face from the metal surface: cis (Z).' : 'The trans vinyl radical/anion is the lower-energy geometry and it is retained: trans (E).' }];
      }
      if (s.kind === 'imine') {
        const g = clone(g0);
        const b = bondOf(g, s.c.id, s.n.id); if (b) b.order = 1;
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.c.id, s.n.id]), added: [], share: 1, label: 'amine (C=N reduced)',
          why: 'The C=N is reduced to a C–N single bond, and what you have is an AMINE. ' +
            (ruleId === 'nabh3cn'
              ? 'Sodium cyanoborohydride is the reagent that makes reductive amination a one-pot procedure. The cyano group pulls electron density off boron, so it is a far weaker hydride donor than NaBH₄ — weak enough to ignore a ketone at pH 6–7, but still able to reduce the much more electrophilic protonated IMINIUM as soon as it forms. So the carbonyl, the amine and the hydride can all sit in the flask together: the equilibrium keeps making a little iminium and the hydride keeps taking it out. Sodium triacetoxyborohydride does the same job without the cyanide.'
              : 'On a catalyst surface a C=N behaves much as a C=C does. Run with ammonia or an amine in the same pot this is REDUCTIVE AMINATION — the most useful way there is to make an amine, because the obvious alternative, alkylating an amine with a halide, cannot be stopped after one alkylation.') + ' ' +
            'Count the nitrogens and the class of the product follows: ammonia gives a primary amine, a primary amine a secondary one, a secondary amine a tertiary one. The nitrogen gains exactly one carbon.',
          stereoHow: 'The C=N is flat, so hydride arrives from either face and a new stereocentre comes out racemic. Doing it enantioselectively takes a chiral catalyst, which is an industry in itself.' }];
      }
      if (s.aromaticRing) {
        const g = clone(g0);
        const ids = new Set(s.ringIds || []);
        g.bonds.forEach(b => { if (ids.has(b.a) && ids.has(b.b)) b.order = 1; });
        return [{ graph: g, changed: ids, added: [], share: 1, label: 'fully hydrogenated ring', why: 'Under forcing conditions the aromatic ring is reduced. You cannot stop part-way: the first addition destroys the aromatic stabilisation that was protecting the other two π bonds, so they go immediately.', stereoHow: 'All hydrogens are delivered from the metal face — cis.' }];
      }
      const g = clone(g0); const b = bondOf(g, s.A.id, s.B.id);
      const was = b.order;
      b.order = 1;
      const c1 = atom(g, s.A.id), c2 = atom(g, s.B.id);
      /* syn: both new hydrogens on the metal face → mark a substituent at each new centre with a wedge */
      const newC = [c1, c2].filter(c => C().isStereocenter(g, c));
      if (newC.length === 2) {
        const pick = c => { const cands = nb(g, c.id).filter(n => n.atom.id !== c1.id && n.atom.id !== c2.id); cands.sort((p, q) => (bondInRing(g, c.id, p.atom.id) ? 1 : 0) - (bondInRing(g, c.id, q.atom.id) ? 1 : 0)); return cands[0]; };
        const n1 = pick(c1), n2 = pick(c2);
        if (n1 && n2) markSyn(g, c1.id, n1.atom.id, c2.id, n2.atom.id);
      }
      return [{ graph: g, changed: new Set([c1.id, c2.id]), added: [], share: 1,
        label: was === 3 ? 'alkane (both π bonds reduced)' : 'alkane (syn H₂)',
        why: 'H₂ dissociates on the palladium surface; the ' + (was === 3 ? 'alkyne' : 'alkene') + ' adsorbs flat on the same surface and the two hydrogens are delivered one after the other to the face touching the metal. ' + (was === 3 ? 'The cis-alkene formed first is itself reduced on the spot — stopping there needs a poisoned catalyst (Lindlar).' : ''),
        stereoHow: 'Syn: both new C–H bonds are on the same face (cis on a ring).' }];
    },
    gate: () => null,
    condNotes(ruleId, env) {
      const n = [];
      if (env.P > 1) n.push('H₂ pressure ' + env.P + ' atm: rate scales with the surface coverage of hydrogen, so higher pressure means faster reduction and reaches harder substrates (trisubstituted alkenes, and above ~50 atm with heat, aromatic rings).');
      if (env.P >= 50 && env.T >= 373) n.push('Forcing conditions (≥ 50 atm, ≥ 100 °C): aromatic rings are now reduced as well.');
      if (env.excess) n.push('Excess H₂: every C=C and C≡C is reduced, not just the fastest one.');
      return n;
    },
    byproducts: () => [],
    stages: () => null,
  };

  /* --- alcohol oxidation: PCC, KMnO4/CrO3 --- */
  /* BH₃·THF is normally the hydroboration reagent, but with no alkene in the flask it does the
     job medicinal chemists actually keep it for: it reduces an AMIDE to the amine (and a
     nitrile to the amine) while leaving an ester untouched — the reverse of LiAlH₄'s order. */
  F.amideRed = {
    rules: [], needs: 'amide',
    sites(g) {
      const out = [];
      carbonylSites(g).forEach(s => {
        if (s.cls !== 'amide') return;
        s.toAmine = true; s.score = -0.6 * s.hind;
        s.reason = 'amide — borane coordinates to the carbonyl OXYGEN first, which is why it reduces an amide faster than an ester (an ester’s oxygens are much poorer donors); the reduction then runs to the amine through the iminium';
        out.push(s);
      });
      g.atoms.forEach(a => {
        if (a.element !== 'C') return;
        const t = nb(g, a.id).find(n => n.atom.element === 'N' && n.bond.order === 3);
        if (t) out.push({ cls: 'nitrile', kind: 'nitrile', c: a, n: t.atom, hind: 0, toAmine: true, score: -0.9,
          label: 'nitrile (C' + a.id + '≡N)', reason: 'nitrile — borane reduces it to the primary amine, more slowly than an amide' });
      });
      return out;
    },
    variants(g0, s, ruleId) { return F.reduce.variants(g0, s, 'lialh4'); },
    byproducts: () => [],
    condNotes() { return ['BH₃·THF or BH₃·SMe₂, THF, 0 °C → reflux; then methanol or HCl to break up the amine–borane complex (this workup is not optional — without it the product stays as the borane adduct).',
      'The selectivity is the point: BH₃ reduces amide > nitrile > ester ≈ acid, while LiAlH₄ reduces ester > amide. With both an amide and an ester present, BH₃ is the reagent that touches only the amide.'] },
  };

  F.oxid = {
    rules: ['pcc', 'kmno4', 'naocl', 'dmp', 'swern', 'tempo', 'jones', 'mno2'], needs: 'alcohol', fallback: ['aldOx', 'benzylOx', 'haloform'],
    sites(g, ruleId, env) {
      return alcoholSites(g).map(s => {
        /* MnO₂ is the one oxidant that is picky about WHICH alcohol: it only takes an allylic
           or benzylic one, which is exactly what makes it useful in the middle of a synthesis */
        if (ruleId === 'mno2' && !s.conj && !s.phenol && !s.acid) return Object.assign(s, { inert: 'a plain saturated alcohol — activated MnO₂ only oxidises an ALLYLIC or BENZYLIC alcohol (that is the point of using it), so this one is left alone' });
        if (s.phenol) return Object.assign(s, { inert: 'phenol — the OH is on an aromatic carbon with no C–H to lose, so it is not oxidised at carbon' });
        if (s.acid) return Object.assign(s, { inert: 'a carboxylic acid O–H — that carbon is already at the top oxidation level (C=O plus O–H); there is nothing left to oxidise' });
        if (s.enol) return Object.assign(s, { inert: 'enol — it tautomerises to the carbonyl on its own rather than being oxidised' });
        if (s.hOnC === 0) return Object.assign(s, { inert: 'tertiary alcohol — the carbinol carbon has no hydrogen to remove, so it cannot be oxidised (it would need a C–C bond to break)' });
        s.score = (s.deg === 2 ? 0.4 : 0) + (s.conj ? 0.5 : 0) - 0.5 * s.hind;
        s.reason = s.cls + ' alcohol' + (s.conj ? ', ' + (s.conj === 'aryl' ? 'benzylic' : 'allylic') + ' (a weaker C–H, so faster)' : '') +
          (s.deg === 2 ? ' — for Cr(VI) the hydride transfer from a secondary carbinol C–H is a little faster than from a primary one' : '') +
          (s.hind ? ', but crowded around the carbinol carbon, which slows the bulky oxidant' : '');
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0); const c = atom(g, s.c.id);
      bondOf(g, s.c.id, s.o.id).order = 2;
      /* which oxidants run past the aldehyde to the acid: the ones used in water */
      const toAcid = (ruleId === 'kmno4' || ruleId === 'naocl' || ruleId === 'jones' || (ruleId === 'tempo' && env && env.excess)) && s.hOnC >= 2;
      let added = [];
      if (toAcid) { const o2 = C().addAtomNear(g, 'O', c); added = [o2.id]; }
      const what = toAcid ? 'carboxylic acid' : (s.deg <= 1 ? (s.hOnC >= 2 ? 'aldehyde' : 'ester / carbonyl') : 'ketone');
      return [{ graph: g, changed: new Set([s.c.id, s.o.id].concat(added)), added, share: 1, label: what,
        why: ruleId === 'dmp'
          ? 'The Dess–Martin periodinane — hypervalent iodine(V) — is the modern workhorse. The alcohol oxygen displaces one acetate from iodine, then the carbinol C–H is handed to a second acetate in a five-membered transition state: out comes the ' + (s.deg <= 1 ? 'aldehyde' : 'ketone') + ', acetic acid, and iodine(III). Room temperature, 30 minutes, neutral, no chromium, and it stops dead at the aldehyde because there is no water to hydrate it. Its one vice is cost — and that the iodine(V) reagent is shock-sensitive when dry.'
          : ruleId === 'swern'
          ? 'The Swern oxidation. DMSO and oxalyl chloride first make the chlorodimethylsulfonium salt (this is the step that must be done at −78 °C — above about −60 °C it decomposes); the alcohol attacks sulfur to give an alkoxysulfonium salt, triethylamine removes a proton from one S-methyl to make a sulfonium ylide, and the ylide plucks the carbinol C–H intramolecularly. The ' + (s.deg <= 1 ? 'aldehyde' : 'ketone') + ' is released with dimethyl sulfide — which is why a Swern is recognisable down the corridor.'
          : ruleId === 'tempo'
          ? (toAcid ? 'TEMPO (a stable nitroxyl radical) with bleach as the terminal oxidant: the real oxidant is the oxoammonium ion, which takes the alcohol to the aldehyde; in water the aldehyde hydrates and is oxidised again, so the product is the acid. Catalytic in TEMPO, cheap, and metal-free.' : 'TEMPO with bleach (or BAIB) as the terminal oxidant. The oxoammonium ion formed from TEMPO takes the hydride from the carbinol carbon; the hydroxylamine left behind is re-oxidised, so a few mol % of TEMPO turns over the whole reaction. Metal-free, and selective for a PRIMARY alcohol over a secondary one — which is unusual and useful.')
          : ruleId === 'jones'
          ? (toAcid ? 'Jones oxidation — CrO₃ in aqueous sulfuric acid. The chromate ester forms and collapses to the aldehyde, which hydrates in the water and is oxidised straight on to the carboxylic acid. Brutal but cheap; nothing acid-sensitive survives it.' : 'Jones oxidation (CrO₃ / H₂SO₄, acetone/water): the secondary alcohol goes to the ketone through the chromate ester, and stops there because that carbon has no hydrogen left.')
          : ruleId === 'mno2'
          ? 'Activated manganese dioxide, heterogeneous and mild, and selective in a way that matters: it oxidises only an ALLYLIC or BENZYLIC alcohol (the C–H is weaker there and the radical intermediate is stabilised), leaving ordinary saturated alcohols untouched. Stir it in CH₂Cl₂ and filter — this is how an allylic alcohol is turned into an enal without touching anything else.'
          : ruleId === 'naocl'
          ? (s.deg <= 1 ? 'Household bleach (NaOCl) in acetic acid makes HOCl, which chlorinates the alcohol oxygen; loss of HCl from the O–Cl (an E2-like step on the C–H) gives the aldehyde, and in water the hydrate is oxidised again to the acid. A cheap, chromium-free oxidant.' : 'Household bleach (NaOCl) in acetic acid makes HOCl; the alcohol oxygen is chlorinated to an alkyl hypochlorite (R₂CH–O–Cl), and base removes the carbinol C–H while chloride leaves — the ketone. No chromium, no heavy-metal waste: the green replacement for the Jones oxidation.')
          : ruleId === 'pcc'
          ? (s.deg <= 1 ? 'PCC is a mild Cr(VI) oxidant used in dry CH₂Cl₂: it removes the carbinol C–H and O–H to give the aldehyde and STOPS — without water the aldehyde cannot form the hydrate that a second oxidation would need.' : 'PCC takes the secondary alcohol to the ketone; with no second C–H on that carbon it can go no further.')
          : (toAcid ? 'Hot aqueous KMnO₄ (or Jones CrO₃) oxidises the primary alcohol to the aldehyde, the aldehyde hydrates in water, and the hydrate is oxidised again — straight through to the carboxylic acid.' : 'The secondary alcohol is oxidised to the ketone, which has no C–H left on that carbon and stops there.'),
        stereoHow: s.deg === 2 ? 'The carbinol carbon becomes sp² and loses any stereochemistry it had.' : null }];
    },
    gate(g, ruleId, env) {
      if (ruleId === 'swern' && env.T > 220 && !env.TAuto) return 'The Swern must be run at −78 °C. The chlorodimethylsulfonium intermediate decomposes above about −60 °C (to a species that chlorinates the alcohol instead), so a warm Swern gives the alkyl chloride, not the carbonyl.';
      if (ruleId === 'dmp' && env.solventVal === 'water') return 'The Dess–Martin periodinane is used in dry CH₂Cl₂. Water hydrolyses the iodine(V) reagent, and in water the aldehyde would be carried on to the acid — which is the opposite of why DMP is chosen.';
      if (ruleId === 'pcc' && env.solventVal === 'water') return 'PCC is used in DRY CH₂Cl₂. In water the aldehyde hydrates and is oxidised again, and the reagent itself decomposes — that is the KMnO₄/Jones outcome, not the PCC one.';
      return null;
    },
    condNotes(ruleId, env) {
      const n = [];
      if (ruleId === 'kmno4') {
        if (env.pHVal < 3) n.push('Acidic permanganate (pH ' + env.pHVal.toFixed(1) + ') is the strongest form: alkenes present would be CLEAVED, not just dihydroxylated.');
        if (env.T > 333) n.push('Hot KMnO₄ (' + degC(env.T) + ' °C) also cleaves any C=C in the molecule.');
      }
      if (env.excess) n.push('Excess oxidant: every oxidisable alcohol goes, not only the fastest one.');
      else n.push('One equivalent: the fastest alcohol is oxidised first; the others wait.');
      return n;
    },
    byproducts: (ruleId) => ruleId === 'pcc' ? [{ smiles: null, name: 'Cr(III) residues + pyridinium chloride' }] : ruleId === 'naocl' ? [{ smiles: '[Na+].[Cl-]', name: 'NaCl (the hypochlorite is reduced to chloride)' }, { smiles: 'O', name: 'water' }] : [{ smiles: null, name: 'MnO₂ (brown) / Mn²⁺' }],
    stages: () => null,
    crossReact(g, ruleId, env) {                   // KMnO4 also attacks alkenes
      if (ruleId !== 'kmno4') return [];
      const al = alkeneSites(g).filter(s => !s.aromatic);
      if (!al.length) return [];
      const s = al[0];
      const out = F.diol.variants(g, s)[0];
      out.label = 'syn-diol at the C=C (cold, dilute KMnO₄)';
      out.why = 'Permanganate is not selective for the alcohol: cold, dilute, basic KMnO₄ adds to the C=C first (syn dihydroxylation through a cyclic manganate ester), and hot or acidic KMnO₄ cleaves it. To oxidise only the alcohol in the presence of an alkene use PCC.';
      out.share = env.T <= 283 ? 0.6 : 0.25;
      out.crossTag = 'alkene';
      return [out];
    },
  };

  /* --- PBr3 --- */
  F.pbr3 = {
    rules: ['pbr3'], needs: 'alcohol',
    sites(g) {
      return alcoholSites(g).map(s => {
        if (s.phenol) return Object.assign(s, { inert: 'phenol — an aryl C–O cannot be displaced by bromide (no backside for SN2, no aryl cation)' });
        if (s.acid) return Object.assign(s, { inert: 'a carboxylic acid O–H — PBr₃ would turn it into the acyl bromide (R–CO–Br), a different reaction not drawn here' });
        if (s.enol) return Object.assign(s, { inert: 'enol — tautomerises to the carbonyl instead' });
        if (s.deg >= 3) return Object.assign(s, { inert: 'tertiary alcohol — the SN2-type displacement of the phosphite is blocked by the three groups; a 3° alcohol goes to the bromide through its carbocation with HBr instead' });
        s.score = (s.deg <= 1 ? 1 : 0) - 1.2 * s.hind + (s.conj ? 0.5 : 0);
        s.reason = s.cls + ' alcohol — bromide displaces the phosphite ester from the back, so the less crowded carbinol carbon reacts faster (1° > 2°)' + (s.hind ? '; branching next to it gets in the way' : '');
        return s;
      });
    },
    variants(g0, s) {
      const g = clone(g0);
      C().removeAtom(g, s.o.id);
      const br = C().addAtomNear(g, 'Br', atom(g, s.c.id));
      /* the new atom sits where the O was, and the centre inverts */
      const oB = atom(g0, s.o.id); const b = atom(g, br.id); b.x = oB.x; b.y = oB.y;
      let inverted = false;
      if (C().isStereocenter(g0, s.c)) inverted = invertAt(g0, g, s.c.id, s.o.id, br.id);
      return [{ graph: g, changed: new Set([s.c.id, br.id]), added: [br.id], removed: [s.o.id], share: 1, label: 'alkyl bromide (inversion)',
        why: 'The alcohol oxygen attacks phosphorus, turning the OH into a good leaving group (–OPBr₂); bromide then displaces it in one backside step. No carbocation, so no rearrangement.',
        stereoHow: inverted ? 'The carbon was a defined stereocentre and it INVERTS — the bromine is drawn on the opposite face from where the OH was.' : 'SN2-type displacement: a stereocentre here would invert.' }];
    },
    gate: (g, ruleId, env) => env.solventVal === 'water' ? 'PBr₃ hydrolyses instantly in water (to H₃PO₃ and HBr). It is used in dry ether at 0 °C.' : null,
    condNotes(ruleId, env) { return env.excess ? ['Excess PBr₃: every primary and secondary OH is converted; one PBr₃ can serve three alcohols.'] : ['One equivalent (of the three Br it carries): the least hindered OH reacts first.']; },
    byproducts: () => [{ smiles: 'OP(O)O', name: 'phosphorous acid, H₃PO₃' }],
    stages: () => null,
  };

  /* --- E1 dehydration --- */
  F.dehyd = {
    rules: ['h2so4_dehydrate'], needs: 'alcohol',
    sites(g, ruleId, env) {
      return alcoholSites(g).map(s => {
        if (s.phenol) return Object.assign(s, { inert: 'phenol — there is no sp³ carbinol carbon and no β-hydrogen; phenols do not dehydrate' });
        if (s.acid) return Object.assign(s, { inert: 'a carboxylic acid O–H — no carbinol carbon, so no carbocation and no alkene; acids survive hot acid' });
        if (s.enol) return Object.assign(s, { inert: 'enol — tautomerises to the carbonyl instead' });
        const betas = nb(g, s.c.id).filter(n => n.atom.element === 'C' && n.bond.order === 1 && hOn(g, n.atom) > 0);
        if (!betas.length) return Object.assign(s, { inert: 'no β-hydrogen next to this carbinol carbon, so no alkene can form from it' });
        const need = s.deg >= 3 ? 298 : s.deg === 2 ? 373 : 443;
        const cat = cationScore(g, s.c, null);
        s.needT = need; s.betas = betas;
        if (env.T < need - 15) return Object.assign(s, { inert: s.cls + ' alcohol — its ' + (s.deg >= 3 ? 'tertiary' : s.deg === 2 ? 'secondary' : 'primary') + ' carbocation only forms fast enough above about ' + degC(need) + ' °C; at ' + degC(env.T) + ' °C it just sits protonated' });
        s.score = cat;
        s.reason = s.cls + ' alcohol — after protonation water leaves to give a ' + catLabel(g, s.c, null) + ' carbocation; the more stable that cation, the lower the temperature needed (3° ~25–80 °C, 2° ~100 °C, 1° ~170 °C)';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const out = [];
      const build = (catId, betaId, chg) => {
        const g = clone(g0);
        C().removeAtom(g, s.o.id);
        if (chg) chg(g);
        bondOf(g, catId, betaId).order = 2;
        return g;
      };
      const alkeneSub = (g, a, b) => deg(g, atom(g, a)) - 1 + deg(g, atom(g, b)) - 1;
      const addForms = (catId, betas, baseShare, pre, whyPre) => {
        const cands = betas.map(bt => { const g = build(catId, bt.atom.id, pre); return { beta: bt, g, sub: alkeneSub(g, catId, bt.atom.id), ring: bondInRing(g, catId, bt.atom.id) }; });
        const fr = share(cands.map(c => 0.7 * c.sub));                                // ~83:17 per substitution step
        cands.forEach((c, i) => {
          const lab = ['ethene', 'monosubstituted', 'disubstituted', 'trisubstituted', 'tetrasubstituted'][Math.min(4, c.sub)] + ' alkene';
          const ez = ezForms(c.g, catId, c.beta.atom.id, 0.75);
          const zaitsev = i === 0 && cands.length > 1;
          const why = whyPre + (cands.length > 1 ? (zaitsev ? 'Zaitsev: losing the β-H that gives the MORE substituted (more stable) C=C is the major path — about 83:17 per extra alkyl group on the double bond.' : 'The less substituted alkene is the minor (Hofmann-type) product; the carbocation is only weakly selective, so a real amount still forms.') : 'Only one kind of β-hydrogen is available, so only this alkene can form.');
          if (ez) ez.forEach(f => out.push({ graph: f.graph, changed: new Set([catId, c.beta.atom.id]), added: [], removed: [s.o.id], beta: c.beta.atom.id, cat: catId, share: baseShare * fr[i] * f.share, label: lab + ' ' + f.label, why: why + ' ' + f.why, stereoHow: null }));
          else out.push({ graph: c.g, changed: new Set([catId, c.beta.atom.id]), added: [], removed: [s.o.id], beta: c.beta.atom.id, cat: catId, share: baseShare * fr[i], label: lab, why, stereoHow: null });
        });
      };
      /* rearrangement of the first cation? */
      let shift = null;
      const cur = cationScore(g0, s.c, null);
      nb(g0, s.c.id).forEach(n => {
        if (n.atom.element !== 'C' || n.bond.order !== 1) return;
        const N = n.atom;
        const gain = (deg(g0, N) <= 1 ? 0 : deg(g0, N) === 2 ? 2 : 4) + (conjugatedTo(g0, N.id, s.c.id) ? 2.5 : 0);
        if (hOn(g0, N) >= 1 && gain > cur + 0.5 && (!shift || gain > shift.gain)) shift = { N, gain, kind: 'hydride' };
      });
      if (shift) {
        const N = shift.N;
        const betas2 = nb(g0, N.id).filter(n => n.atom.element === 'C' && n.bond.order === 1 && n.atom.id !== s.c.id && hOn(g0, n.atom) > 0)
          .concat([{ atom: s.c, bond: bondOf(g0, s.c.id, N.id) }]);           // the original carbinol carbon now has an extra H
        addForms(N.id, betas2, 0.65, null, 'The ' + catLabel(g0, s.c, null) + ' carbocation formed first is next to a ' + (deg(g0, N) >= 3 ? 'tertiary' : 'secondary') + ' C–H: a 1,2-hydride shift moves the + charge there (a more stable cation) BEFORE the β-hydrogen is lost, so the alkene grows from the rearranged cation. ');
        addForms(s.c.id, s.betas, 0.35, null, 'Some of the first-formed cation loses a β-H before it can rearrange. ');
      } else {
        addForms(s.c.id, s.betas, 1, null, 'Protonation of the OH, loss of water to the ' + catLabel(g0, s.c, null) + ' carbocation, then a β-hydrogen is removed. ');
      }
      /* fold identical structures together (the same alkene reachable two ways) */
      const seen = new Map();
      out.forEach(o => { const k = keyOf(o.graph); if (seen.has(k)) { seen.get(k).share += o.share; } else seen.set(k, o); });
      return [...seen.values()];
    },
    gate(g, ruleId, env) {
      if (env.pHVal > 1.5) return 'Dehydration needs the OH protonated to make water the leaving group — that takes concentrated acid (pH < 1). At pH ' + env.pHVal.toFixed(1) + ' the alcohol is essentially all un-protonated, and OH⁻ itself never leaves.';
      return null;
    },
    condNotes(ruleId, env) {
      const n = ['Heat drives it: the elimination makes two molecules from one, so TΔS pushes the equilibrium toward the alkene and water at high temperature.'];
      if (env.T < 373) n.push('At ' + degC(env.T) + ' °C only tertiary alcohols dehydrate at a useful rate (secondary need ~100 °C, primary ~170 °C).');
      return n;
    },
    byproducts: () => [{ smiles: 'O', name: 'water' }],
    stages(g0, s, out) {
      const g = clone(g0); C().removeAtom(g, s.o.id);
      /* find which carbon holds the charge in this outcome: the alkene carbon that was the carbinol carbon or the shifted one */
      const ids = [...out.changed];
      const catId = ids.includes(s.c.id) && bondOf(out.graph, s.c.id, ids.find(i => i !== s.c.id)) ? s.c.id : ids[0];
      const c = atom(g, catId); if (c) c.charge = 1;
      return [{ label: 'carbocation (after loss of water)', graph: g }];
    },
  };

  /* --- substitution / elimination on alkyl halides --- */
  const NUC_KEY = { naoh: 'NaOH', naoet: 'NaOEt', tbuok: 'tBuOK', nacn: 'NaCN' };
  const NUC_TEXT = { naoh: 'hydroxide', naoet: 'ethoxide', tbuok: 'tert-butoxide', nacn: 'cyanide' };
  F.halide = {
    rules: ['naoh', 'naoet', 'tbuok', 'nacn'], needs: 'halide', fallback: ['saponify', 'amideHyd', 'acidBase', 'cannizzaro', 'enol'],
    sites(g, ruleId, env) {
      return halideSites(g).map(s => {
        if (s.cls === 'aryl') return Object.assign(s, { inert: 'aryl halide — the ring blocks backside attack and an aryl cation is hopeless; neither SN2 nor SN1 (nor E2) is available' });
        if (s.cls === 'vinyl') return Object.assign(s, { inert: 'vinyl halide — sp² carbon, no SN2 backside and no stable cation; inert to these reagents' });
        if (s.el === 'F') return Object.assign(s, { inert: 'C–F — fluoride is a terrible leaving group (the strongest C–halogen bond and the most basic halide); it does not leave' });
        const lg = s.el === 'I' ? 1.5 : s.el === 'Br' ? 1.0 : 0;
        const spec = C().NUCLEOPHILES[NUC_KEY[ruleId]];
        const sn2 = s.deg === 0 ? 1.5 : s.deg === 1 ? 0 : s.deg === 2 ? -1.7 : -6;
        const e2 = s.deg === 3 ? 1 : s.deg === 2 ? 0 : -1;
        const bulky = spec.bulky;
        const canE = s.betas.length > 0;
        const path = (ruleId === 'nacn') ? sn2 : (bulky ? (canE ? e2 : sn2) : (s.deg === 3 ? (canE ? e2 : -6) : Math.max(sn2, canE ? e2 - 0.5 : -9)));
        s.score = lg + path - 1.2 * s.hind + (s.conj ? 0.8 : 0);
        s.reason = s.cls + ' C–' + s.el + ' — ' + (s.el === 'I' ? 'iodide is the best leaving group (weakest C–X bond, least basic anion)' : s.el === 'Br' ? 'bromide is a good leaving group' : 'chloride is a fair leaving group, ~10–30× slower than bromide') +
          (s.deg <= 1 ? '; an unhindered carbon is wide open to backside attack' : s.deg === 2 ? '; a secondary carbon is the SN2/E2 borderline' : '; a tertiary carbon is closed to SN2, so only elimination (or SN1) is possible') +
          (s.hind ? '; branching next to it slows SN2 sharply (neopentyl-type crowding)' : '');
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const key = NUC_KEY[ruleId];
      const spec = C().NUCLEOPHILES[key];
      const heat = env.T >= 323;
      let mech = null; try { mech = C().chooseMechanism(g0, s.c.id, key, { heat }); } catch (e) {}
      const canE = s.betas.length > 0;
      /* base E:S odds by substrate class and reagent, then conditions */
      let eFrac;
      if (ruleId === 'nacn') eFrac = s.deg >= 3 ? 0.85 : s.deg === 2 ? 0.10 : 0.0;
      else if (ruleId === 'tbuok') eFrac = s.deg >= 3 ? 1.0 : s.deg === 2 ? 0.97 : 0.85;
      else eFrac = s.deg >= 3 ? 0.97 : s.deg === 2 ? 0.79 : 0.10;
      if (!canE) eFrac = 0;
      const odds0 = eFrac / Math.max(1e-6, 1 - eFrac);
      let odds = odds0 * Math.pow(2, (env.T - 298) / 30);            // heat favours elimination (~×2 per 30 °C)
      if (env.solventVal === 'aprotic') odds /= 3;                     // naked anion: SN2 speeds up more than E2
      if (env.solventVal === 'nonpolar') odds *= 1.5;
      if (eFrac >= 0.999) odds = 1e6; if (eFrac <= 0.001) odds = 0;
      eFrac = odds / (1 + odds);
      const out = [];
      /* substitution product */
      if (eFrac < 0.995) {
        const g = clone(g0);
        C().removeAtom(g, s.x.id);
        const xB = atom(g0, s.x.id);
        let nuId; const addedIds = [];
        if (ruleId === 'nacn') { const cn = C().addAtomNear(g, 'C', atom(g, s.c.id)); const n = C().addAtomNear(g, 'N', cn); bondOf(g, cn.id, n.id).order = 3; nuId = cn.id; addedIds.push(cn.id, n.id); }
        else if (ruleId === 'naoet') { const o = C().addAtomNear(g, 'O', atom(g, s.c.id)); const c1 = C().addAtomNear(g, 'C', o); const c2 = C().addAtomNear(g, 'C', c1); nuId = o.id; addedIds.push(o.id, c1.id, c2.id); }
        else if (ruleId === 'tbuok') { const o = C().addAtomNear(g, 'O', atom(g, s.c.id)); const cq = C().addAtomNear(g, 'C', o); const m1 = C().addAtomNear(g, 'C', cq); const m2 = C().addAtomNear(g, 'C', cq); const m3 = C().addAtomNear(g, 'C', cq); nuId = o.id; addedIds.push(o.id, cq.id, m1.id, m2.id, m3.id); }
        else { const o = C().addAtomNear(g, 'O', atom(g, s.c.id)); nuId = o.id; addedIds.push(o.id); }
        const nu = atom(g, nuId); nu.x = xB.x; nu.y = xB.y;
        let inverted = false;
        const sn1 = mech && mech.major === 'SN1';
        if (!sn1 && C().isStereocenter(g0, s.c)) inverted = invertAt(g0, g, s.c.id, s.x.id, nuId);
        const what = ruleId === 'nacn' ? 'nitrile (SN2, chain +1 carbon)' : ruleId === 'naoh' ? 'alcohol (SN2)' : ruleId === 'naoet' ? 'ethyl ether (SN2, Williamson)' : 'tert-butyl ether';
        out.push({ graph: g, changed: new Set([s.c.id, nuId]), added: addedIds, removed: [s.x.id], share: 1 - eFrac, label: what,
          why: (sn1 ? 'SN1: the leaving group departs on its own to a carbocation, which the nucleophile then captures from either face.' :
            NUC_TEXT[ruleId].charAt(0).toUpperCase() + NUC_TEXT[ruleId].slice(1) + ' attacks the ' + s.cls + ' carbon from the side opposite the ' + s.el + ' while the C–' + s.el + ' bond breaks — one step, no intermediate (SN2).') +
            (canE && eFrac > 0.05 ? ' It competes with elimination: the same anion can take a β-hydrogen instead.' : ''),
          stereoHow: sn1 ? 'A flat carbocation is captured from either face — racemic.' : (inverted ? 'INVERSION: the centre turns inside-out like an umbrella, so the nucleophile is drawn on the face opposite the one the halogen left.' : 'Backside attack inverts a stereocentre (this carbon has none defined).'),
          mech: sn1 ? 'SN1' : 'SN2' });
      }
      /* elimination products */
      if (canE && eFrac > 0.005) {
        const hof = ruleId === 'tbuok';
        const cands = s.betas.map(bt => {
          const g = clone(g0); C().removeAtom(g, s.x.id); bondOf(g, s.c.id, bt.carbon.id).order = 2;
          return { bt, g, sub: deg(g, atom(g, s.c.id)) - 1 + deg(g, atom(g, bt.carbon.id)) - 1 };
        });
        let fr;
        if (cands.length === 1) fr = [1];
        else {
          const subs = cands.map(c => c.sub);
          const maxS = Math.max(...subs), minS = Math.min(...subs);
          /* Brown's numbers: 3° halide NaOEt 71:29 Zaitsev, tBuOK 28:72; 2° halide NaOEt 81:19, tBuOK 47:53 */
          const zai = hof ? (s.deg >= 3 ? 0.28 : 0.47) : (s.deg >= 3 ? 0.71 : 0.81);
          fr = cands.map(c => c.sub === maxS ? zai : (1 - zai));
          const nMax = cands.filter(c => c.sub === maxS).length, nMin = cands.length - nMax;
          fr = cands.map(c => c.sub === maxS ? zai / nMax : (1 - zai) / Math.max(1, nMin));
          if (maxS === minS) fr = cands.map(() => 1 / cands.length);
        }
        const mechName = (mech && (mech.major === 'E1' || mech.major === 'E2')) ? mech.major : 'E2';
        cands.forEach((c, i) => {
          const lab = ['ethene', 'monosubstituted', 'disubstituted', 'trisubstituted', 'tetrasubstituted'][Math.min(4, c.sub)] + ' alkene';
          const isMajorRegio = c.sub === Math.max(...cands.map(x => x.sub));
          const why = (cands.length > 1
            ? (hof ? (isMajorRegio ? 'Zaitsev alkene — the more substituted C=C — but tert-butoxide is too bulky to reach the crowded internal β-H easily, so this is the MINOR alkene (H. C. Brown: 28% with a 3° bromide, 47% with 2-bromobutane).'
                                     : 'Hofmann alkene: the bulky base takes the most exposed β-hydrogen (the CH₃ / CH₂ end), giving the LESS substituted C=C as the major alkene (72% with a 3° bromide, 53% with 2-bromobutane).')
                   : (isMajorRegio ? 'Zaitsev alkene: a small base removes the β-H that gives the MORE substituted, more stable C=C (H. C. Brown: 71% with a 3° bromide, 81% with 2-bromobutane, NaOEt/EtOH).'
                                     : 'The less substituted alkene is the minor product with a small base — the transition state leading to it is less alkene-like and less stabilised.'))
            : 'Only one kind of β-hydrogen is available, so only this alkene can form.') +
            ' The β-H and the ' + s.el + ' must be anti-periplanar (180° apart) in the same plane for the orbitals to overlap as the C=C forms.';
          const ez = ezForms(c.g, s.c.id, c.bt.carbon.id, 0.67);
          if (ez) ez.forEach(f => out.push({ graph: f.graph, changed: new Set([s.c.id, c.bt.carbon.id]), added: [], removed: [s.x.id], beta: c.bt.carbon.id, share: eFrac * fr[i] * f.share, label: lab + ' ' + f.label, why: why + ' ' + f.why, stereoHow: null, mech: mechName }));
          else out.push({ graph: c.g, changed: new Set([s.c.id, c.bt.carbon.id]), added: [], removed: [s.x.id], beta: c.bt.carbon.id, share: eFrac * fr[i], label: lab, why, stereoHow: null, mech: mechName });
        });
      }
      const seen = new Map();
      out.forEach(o => { const k = keyOf(o.graph); if (seen.has(k)) seen.get(k).share += o.share; else seen.set(k, o); });
      return [...seen.values()].sort((p, q) => q.share - p.share);
    },
    gate(g, ruleId, env) {
      if (env.pHVal < 7) return NUC_TEXT[ruleId].charAt(0).toUpperCase() + NUC_TEXT[ruleId].slice(1) + ' is a base/anion: in an acidic medium (pH ' + env.pHVal.toFixed(1) + ') it is protonated to ' + (ruleId === 'nacn' ? 'HCN' : ruleId === 'naoh' ? 'water' : 'the alcohol') + ', which is neither a strong nucleophile nor a strong base. No reaction.';
      return null;
    },
    condNotes(ruleId, env) {
      const n = [];
      if (env.T >= 323) n.push('Heat (' + degC(env.T) + ' °C) favours elimination over substitution — it makes more molecules, so TΔS is on its side, and E2 has the higher activation energy.');
      if (env.solventVal === 'aprotic') n.push('Polar aprotic solvent (DMSO/DMF): the anion is not hydrogen-bonded, so it is far more nucleophilic — SN2 speeds up more than E2 and substitution gains share.');
      if (env.solventVal === 'protic' || env.solventVal === 'water') n.push('Protic solvent: hydrogen bonding cages the anion and slows SN2; with a tertiary or benzylic substrate it also opens the SN1/E1 door.');
      if (env.solventVal === 'nonpolar') n.push('Non-polar solvent: the ionic reagent barely dissolves — expect a slow, heterogeneous reaction.');
      return n;
    },
    byproducts: (ruleId, env, s) => [{ smiles: null, name: 'Na' + (s ? s.el : 'X') + ' (salt)' }].concat(
      ruleId === 'naoh' ? [{ smiles: 'O', name: 'water (from E2)' }] : ruleId === 'naoet' ? [{ smiles: 'CCO', name: 'ethanol (from E2)' }] : ruleId === 'tbuok' ? [{ smiles: 'CC(C)(C)O', name: 'tert-butanol (from E2)' }] : []),
    stages(g0, s, out) {
      if (out.mech === 'SN1' || out.mech === 'E1') { const g = clone(g0); C().removeAtom(g, s.x.id); const c = atom(g, s.c.id); if (c) c.charge = 1; return [{ label: 'carbocation', graph: g }]; }
      return null;
    },
  };

  /* --- hydride reductions --- */
  F.reduce = {
    rules: ['nabh4', 'lialh4', 'dibal'], needs: 'carbonyl',
    sites(g, ruleId, env) {
      const out = carbonylSites(g).map(s => {
        const strong = ruleId === 'lialh4' || ruleId === 'dibal';
        /* DIBAL exists for one job: stopping an ester at the aldehyde */
        if (ruleId === 'dibal' && s.cls === 'amide') return Object.assign(s, { inert: 'amide — DIBAL takes it to the ALDEHYDE too (through the same tetrahedral aluminium complex), a transformation not drawn here' });
        if (s.cls === 'amide' && ruleId === 'lialh4') { s.toAmine = true; s.score = -2.4 - 0.6 * s.hind;   // the LEAST reactive carbonyl: an ester or an acid in the same molecule goes first
          s.reason = 'amide — the slowest carbonyl of all to reduce (the nitrogen lone pair feeds the C=O and kills its electrophilicity), which is why it needs LiAlH₄ and not NaBH₄; and it goes all the way past the alcohol stage to the AMINE, because the nitrogen never leaves';
          return s; }
        if (s.cls === 'amide') return Object.assign(s, { inert: strong ? 'amide — DIBAL takes it to the aldehyde by design, not to the amine' : 'amide — the nitrogen lone pair donates into the C=O, making it far too weak an electrophile for NaBH₄' });
        if (s.cls === 'ester' && !strong) return Object.assign(s, { inert: 'ester — the OR oxygen’s lone pair feeds the carbonyl and makes it a poor electrophile; NaBH₄ is not a strong enough hydride donor (LiAlH₄ would reduce it)' });
        if (s.cls === 'carboxylic acid' && !strong) return Object.assign(s, { inert: 'carboxylic acid — NaBH₄ is quenched by the acidic O–H (H₂ evolves) and the carboxylate anion left behind is not electrophilic' });
        s.score = (s.cls === 'aldehyde' ? 1.0 : s.cls === 'acyl chloride' ? 1.5 : s.cls === 'ketone' ? 0 : s.cls === 'ester' ? -1.5 : s.cls === 'carboxylic acid' ? -2 : -0.5) - 0.6 * s.hind - (s.conjugated ? 0.4 : 0);
        s.reason = s.cls + (s.cls === 'aldehyde' ? ' — the least hindered, most electrophilic carbonyl (only one alkyl group donating into it), so it is reduced first' :
          s.cls === 'ketone' ? ' — two alkyl groups donate into the C=O and crowd it, so it is slower than an aldehyde' :
          s.cls === 'ester' ? ' — the OR lone pair makes it a weak electrophile; LiAlH₄ takes it only after any aldehyde/ketone, and removes the OR to give the primary alcohol' :
          s.cls === 'carboxylic acid' ? ' — first deprotonated by the hydride (H₂ evolves), then the carboxylate is reduced slowly to the primary alcohol' : '') +
          (s.hind ? '; crowded, which slows hydride delivery' : '') + (s.conjugated ? '; conjugated with C=C (some 1,4-reduction can compete)' : '');
        return s;
      });
      /* a nitrile is reduced to a PRIMARY AMINE by LiAlH₄ (two hydrides, through the imine
         anion) — and to the aldehyde by DIBAL, which is that reagent's other reason to exist */
      if (ruleId === 'lialh4' || ruleId === 'dibal') {
        g.atoms.forEach(a => {
          if (a.element !== 'C') return;
          const t = nb(g, a.id).find(n => n.atom.element === 'N' && n.bond.order === 3);
          if (!t) return;
          out.push({ cls: 'nitrile', kind: 'nitrile', c: a, n: t.atom, hind: 0, toAmine: ruleId === 'lialh4',
            score: -0.9, label: 'nitrile (C' + a.id + '≡N)',
            reason: ruleId === 'lialh4'
              ? 'nitrile — reduced past the imine stage to the primary amine; slower than a ketone, so a ketone in the same molecule goes first'
              : 'nitrile — DIBAL stops at the metalated imine, which hydrolyses to the ALDEHYDE on workup (the Stephen-type reduction)' });
        });
      }
      return out;
    },
    variants(g0, s, ruleId) {
      const g = clone(g0);
      if (s.cls === 'nitrile') {
        const bd = bondOf(g, s.c.id, s.n.id);
        if (s.toAmine) {
          bd.order = 1;
          try { C().layoutGraph(g); } catch (e) {}
          return [{ graph: g, changed: new Set([s.c.id, s.n.id]), added: [], share: 1, label: 'primary amine (nitrile reduced)',
            why: 'Two hydrides, one after the other. The first adds to the nitrile carbon and gives the metalated imine (C=N–[Al]); unlike an ester’s tetrahedral intermediate this cannot expel anything, so it simply takes a second hydride and the aqueous workup gives the primary amine. ' +
              'This is the standard way to lengthen a chain by one carbon AND end with an amine: alkyl halide → nitrile (NaCN) → amine (LiAlH₄). DIBAL at −78 °C stops the same reaction at the imine, which hydrolyses to the aldehyde instead.',
            stereoHow: null }];
        }
        bd.order = 2;
        const nAt = g.atoms.find(a => a.id === s.n.id);
        if (nAt) { nAt.element = 'O'; }
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.c.id, s.n.id]), added: [], share: 1, label: 'aldehyde (DIBAL stops at the imine)',
          why: 'One hydride gives the metalated imine, which is stable at −78 °C. There is nothing for it to expel, and no second hydride is delivered, so the acid workup hydrolyses C=N to C=O: the nitrile has become an aldehyde. LiAlH₄ under the same conditions would carry on to the primary amine.',
          stereoHow: null, byprod: [{ smiles: 'N', name: 'NH₃ (from hydrolysis of the imine)' }] }];
      }
      if (s.cls === 'amide' && s.toAmine) {
        /* the carbonyl oxygen leaves altogether and the C–N bond stays: C=O becomes CH2 */
        try { C().removeAtom(g, s.o.id); } catch (e) {}
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.c.id]), added: [], share: 1, label: 'amine (amide reduced, C=O → CH₂)',
          why: 'An amide is reduced to the AMINE, not to an alcohol, and the reason is what the tetrahedral intermediate can throw out. The first hydride gives an alkoxide with both a nitrogen and an O⁻ on the same carbon; aluminium takes the oxygen and what leaves is the OXIDE, giving an iminium ion. That is a far better electrophile than the amide ever was, so the second hydride is delivered at once and the nitrogen ends up on a CH₂. ' +
            'Nothing stops halfway: there is no way to isolate the hemiaminal. Note the selectivity — BH₃·THF reduces an amide and leaves an ESTER alone, which is exactly the opposite of LiAlH₄, and is the usual choice when both are present.',
          stereoHow: null }];
      }
      const c = atom(g, s.c.id);
      bondOf(g, s.c.id, s.o.id).order = 1;
      let added = [];
      let label = s.cls === 'aldehyde' ? 'primary alcohol' : s.cls === 'ketone' ? 'secondary alcohol' : 'alcohol';
      if (ruleId === 'dibal' && s.cls === 'ester') {
        /* one hydride only: the tetrahedral aluminium alkoxide is stable at −78 °C and does not
           collapse until the workup, so the aldehyde is what comes out */
        const lv = nb(g0, s.c.id).find(n => n.atom.id !== s.o.id && n.atom.element === 'O' && n.bond.order === 1);
        if (lv) {
          const half = cutFragment(g, s.c.id, lv.atom.id);
          bondOf(g, s.c.id, s.o.id).order = 2;          // it stops AT the aldehyde: C=O stays
          try { C().layoutGraph(g); } catch (e) {}
          return [{ graph: g, changed: new Set([s.c.id, s.o.id]), added: [], share: 1, label: 'aldehyde (DIBAL, one hydride only)',
            why: 'DIBAL-H (diisobutylaluminium hydride) at −78 °C delivers exactly ONE hydride to the ester. The tetrahedral intermediate that results is an aluminium alkoxide, and — unlike the lithium version from LiAlH₄ — it is stable at that temperature and does NOT collapse to expel the alkoxide. So no second hydride can be delivered, and when the aldehyde is finally released at the workup there is no reducing agent left to touch it. This is the standard way to get an aldehyde from an ester, and it works for exactly the same reason a Weinreb amide does: the intermediate is held together until the acid is added.',
            stereoHow: null, _byover: [{ graph: bmTidy(half), name: 'the alcohol (from the ester’s OR group)' }] }];
        }
      }
      if (s.cls === 'ester' || s.cls === 'carboxylic acid' || s.cls === 'acyl chloride') {
        /* the OR / OH / Cl leaves after the first hydride, a second hydride gives the primary alcohol */
        const lv = nb(g0, s.c.id).find(n => n.atom.id !== s.o.id && ((n.atom.element === 'O' && n.bond.order === 1) || C().HALOGENS.includes(n.atom.element)));
        if (lv) {
          const b = bondOf(g, s.c.id, lv.atom.id);
          g.bonds = g.bonds.filter(x => x !== b);
          if (s.cls === 'ester') { /* the alcohol fragment R'OH stays as a second molecule */ }
          else C().removeAtom(g, lv.atom.id);
        }
        label = 'primary alcohol' + (s.cls === 'ester' ? ' + the ester’s alcohol' : '');
      }
      return [{ graph: g, changed: new Set([s.c.id, s.o.id]), added, share: 1, label,
        why: (ruleId === 'nabh4' ? 'BH₄⁻ hands a hydride to the electrophilic carbonyl carbon; the π electrons move onto oxygen, and methanol (the solvent) protonates the alkoxide. ' : 'AlH₄⁻ delivers hydride to the carbonyl carbon; the aqueous workup afterwards protonates the alkoxide. ') +
          (s.cls === 'ester' ? 'The tetrahedral intermediate throws out the OR group to give an aldehyde, which is reduced again at once — so an ester takes two hydrides and gives two alcohols.' :
           s.cls === 'carboxylic acid' ? 'The acid is first deprotonated (H₂ bubbles off); the carboxylate is then reduced, slowly, through the aldehyde stage to the primary alcohol.' : ''),
        stereoHow: s.cls === 'ketone' ? 'The carbonyl is flat, so hydride arrives from either face equally — a new stereocentre is racemic.' : null }];
    },
    gate(g, ruleId, env) {
      if (ruleId === 'nabh4' && env.pHVal < 4) return 'NaBH₄ is destroyed by acid — below pH ~4 it hydrolyses to borate and H₂ faster than it can reduce anything. Keep the medium neutral or mildly basic (MeOH or EtOH).';
      if (ruleId === 'lialh4' && (env.solventVal === 'water' || env.solventVal === 'protic')) return 'LiAlH₄ reacts violently with water and alcohols (it is a hydride source, and O–H is an acid to it) — the reagent is destroyed before it touches the carbonyl. It must be used in DRY ether or THF; water comes only at the workup.';
      return null;
    },
    condNotes(ruleId, env) {
      const n = [];
      if (env.T <= 253) n.push('At ' + degC(env.T) + ' °C the aldehyde/ketone difference is exaggerated — cold NaBH₄ reduces an aldehyde while leaving a ketone almost untouched.');
      if (env.excess) n.push('Excess hydride: every reducible carbonyl is reduced.');
      else n.push('One equivalent of hydride: the most electrophilic, least hindered carbonyl is reduced first.');
      return n;
    },
    byproducts: (ruleId) => ruleId === 'nabh4' ? [{ smiles: 'OB(O)O', name: 'boric acid (as borate) after workup' }] : [{ smiles: null, name: 'Al(OH)₃ + LiOH on aqueous workup' }],
    stages(g0, s) { const g = clone(g0); bondOf(g, s.c.id, s.o.id).order = 1; const o = atom(g, s.o.id); if (o) o.charge = -1; return [{ label: 'alkoxide intermediate', graph: g }]; },
  };

  /* --- acid-catalysed ester hydrolysis: what aqueous H–X / H₃O⁺ does when there is no C=C --- */
  function esterSites(g) {
    return carbonylSites(g).filter(s => s.cls === 'ester').map(s => {
      const eo = nb(g, s.c.id).find(n => n.atom.element === 'O' && n.bond.order === 1 && n.atom.id !== s.o.id);
      if (!eo) return null;
      const r = nb(g, eo.atom.id).find(n => n.atom.id !== s.c.id && n.atom.element === 'C');
      if (!r) return null;
      const aryl = isAromatic(g, r.atom.id);
      s.kind = 'ester'; s.eo = eo.atom; s.r = r.atom; s.aryl = aryl; s.atoms = [s.c.id, s.o.id, eo.atom.id];
      s.score = (aryl ? 1.5 : 0) - 0.8 * s.hind;
      s.label = 'ester C' + s.c.id + '(=O)–O–C' + r.atom.id;
      s.reason = (aryl ? 'a phenyl ester — the phenoxide-like leaving group is good, so it hydrolyses faster than an alkyl ester' : 'an alkyl ester') + (s.hind ? '; branching next to the carbonyl slows the attack of water' : '');
      return s;
    }).filter(Boolean);
  }
  function splitOff(g, keepId) {              // keep the fragment containing keepId, return the other fragment as its own graph
    const seen = new Set([keepId]); const st = [keepId];
    while (st.length) { const id = st.pop(); nb(g, id).forEach(n => { if (!seen.has(n.atom.id)) { seen.add(n.atom.id); st.push(n.atom.id); } }); }
    const other = { atoms: g.atoms.filter(a => !seen.has(a.id)).map(a => Object.assign({}, a)), bonds: g.bonds.filter(b => !seen.has(b.a) && !seen.has(b.b)).map(b => Object.assign({}, b)), nextId: g.nextId };
    g.atoms = g.atoms.filter(a => seen.has(a.id)); g.bonds = g.bonds.filter(b => seen.has(b.a) && seen.has(b.b));
    return other;
  }
  F.esterHyd = {
    rules: [], needs: 'ester',
    sites(g) { return esterSites(g); },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      g.bonds = g.bonds.filter(b => !((b.a === s.c.id && b.b === s.eo.id) || (b.b === s.c.id && b.a === s.eo.id)));
      const oh = C().addAtomNear(g, 'O', atom(g, s.c.id));
      let alcohol = splitOff(g, s.c.id);
      /* the product shown is the bigger fragment; the other is listed as also formed */
      let main = g, changed = new Set([s.c.id, s.o.id, oh.id]), added = [oh.id], mainIs = 'acid';
      if (alcohol.atoms.length > g.atoms.length) { main = alcohol; alcohol = g; changed = new Set([s.eo.id]); added = []; mainIs = 'alcohol'; }
      s.other = alcohol; s.otherIs = mainIs === 'acid' ? (s.aryl ? 'phenol' : 'alcohol') : 'acid';
      const acidName = ruleId === 'hydration' ? 'H₂SO₄ / H₂O' : 'aqueous H–' + (ruleId === 'hbr' ? 'Br' : ruleId === 'hcl' ? 'Cl' : 'I');
      return [{ graph: main, changed, added, share: 1,
        label: mainIs === 'acid' ? 'carboxylic acid (ester hydrolysed) + ' + (s.aryl ? 'the phenol' : 'the alcohol') : (s.aryl ? 'phenol' : 'alcohol') + ' (ester hydrolysed) + the carboxylic acid',
        why: acidName + ' has no C=C to add to, so the acid does the only thing it can here: it catalyses hydrolysis of the ester. H₃O⁺ protonates the carbonyl oxygen (making the carbon a much better electrophile), water adds to give a tetrahedral intermediate, a proton hops onto the O–' + (s.aryl ? 'aryl' : 'alkyl') + ' oxygen, and that ' + (s.aryl ? 'phenol' : 'alcohol') + ' leaves; losing H⁺ from the carbonyl regenerates the catalyst. It is the reverse of Fischer esterification — an equilibrium pushed by the large excess of water.',
        stereoHow: null }];
    },
    gate(g, ruleId, env) {
      if (!env.solventAuto && !['water', 'protic'].includes(env.solventVal))
        return 'In ' + (SOLVENT_NAMES[env.solventVal] || env.solventVal) + ' (no water) H–X has no water to hydrolyse the ester with, and there is no C=C for it to add to — nothing happens. Switch the solvent to water (aqueous acid) to hydrolyse the ester.';
      return null;
    },
    condNotes(ruleId, env) {
      const n = ['No C=C here, so H–X behaves as an aqueous acid catalyst (H₃O⁺). Acid-catalysed ester hydrolysis is slow at room temperature — it is normally run at reflux for hours; every step is reversible, so a large excess of water drives it.'];
      if (env.T < 313) n.push('At ' + degC(env.T) + ' °C expect this to take many hours; heating to reflux (~100 °C) is the usual practice.');
      if (env.T >= 353) n.push('At ' + degC(env.T) + ' °C the hydrolysis runs at a practical rate.');
      return n;
    },
    byproducts(ruleId, env, site) {
      const out = [];
      if (site && site.other) { let nm = ''; try { nm = nameOf(site.other); } catch (e) {} let smi = null; try { smi = C().toSmiles(site.other); } catch (e) {} out.push({ smiles: smi, name: (nm || 'the ' + site.otherIs) + ' (the other half of the ester)', graph: site.other }); }
      return out;
    },
    stages(g0, s) {
      const st = [];
      const g1 = clone(g0); atom(g1, s.o.id).charge = 1; st.push({ label: 'protonated carbonyl (H₃O⁺ lends a proton to C=O)', graph: g1 });
      const g2 = clone(g0); bondOf(g2, s.c.id, s.o.id).order = 1; C().addAtomNear(g2, 'O', atom(g2, s.c.id)); st.push({ label: 'tetrahedral intermediate (water has added)', graph: g2 });
      /* the proton moves to the O–aryl / O–alkyl oxygen: that is what turns it into a leaving group */
      const g3 = clone(g2); atom(g3, s.eo.id).charge = 1; st.push({ label: 'proton on the leaving oxygen (O⁺–H) — ready to leave', graph: g3 });
      return st;
    },
  };



  /* --- H–X on an ALCOHOL: the O–H is protonated and water leaves — Sₙ1 for 3° (and 2°), Sₙ2 for 1° ---
     1-butanol + NaBr / H₂SO₄ → 1-bromobutane; tert-butyl alcohol + conc. HCl → tert-butyl chloride */
  F.alcHX = {
    rules: [], needs: 'alcohol',
    sites(g, ruleId, env) {
      if (ruleId === 'hydration') return [];
      return alcoholSites(g).map(s => {
        if (s.phenol) return Object.assign(s, { inert: 'phenol — an aryl C–O cannot be substituted (no Sₙ2 backside, no aryl cation)' });
        if (s.acid) return Object.assign(s, { inert: 'a carboxylic acid O–H, not an alcohol' });
        if (s.enol) return Object.assign(s, { inert: 'an enol O–H — it tautomerises instead' });
        s.mech = s.deg >= 3 ? 'SN1' : s.deg === 2 ? 'SN1' : 'SN2';
        s.score = s.deg >= 3 ? 2 : s.deg === 2 ? 0.8 : 0.4;
        s.reason = s.deg >= 3 ? 'a tertiary alcohol: protonate the O–H, water leaves, the 3° carbocation forms in seconds and the halide traps it (Sₙ1 — the Lucas test is instant)'
          : s.deg === 2 ? 'a secondary alcohol: Sₙ1 through a 2° carbocation (slow without heat or ZnCl₂; the cation may rearrange)'
          : 'a primary alcohol: the protonated O–H is displaced by the halide from the back (Sₙ2) — needs conc. HBr with H₂SO₄ and heat (HCl needs ZnCl₂)';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const X = ruleId === 'hbr' ? 'Br' : ruleId === 'hcl' ? 'Cl' : 'I';
      const g = clone(g0);
      const c = atom(g, s.c.id);
      C().removeAtom(g, s.o.id);
      const x = C().addAtomNear(g, X, c);
      const cls = s.deg >= 3 ? 'tertiary' : s.deg === 2 ? 'secondary' : 'primary';
      return [{ graph: g, changed: new Set([s.c.id, x.id]), added: [x.id], share: 1, mech: s.mech, mechKind: s.mech,
        label: cls + ' alkyl ' + (X === 'Br' ? 'bromide' : X === 'Cl' ? 'chloride' : 'iodide') + ' (' + s.mech + ' on the protonated alcohol)',
        why: 'An –OH is a hopeless leaving group, but the acid protonates it to –OH₂⁺, and WATER is an excellent leaving group. ' +
          (s.mech === 'SN1' ? 'Water leaves on its own to give the ' + cls + ' carbocation, which the halide ion traps — Sₙ1' + (s.deg === 2 ? ' (a 2° cation: slower, and it can rearrange to a more stable one)' : ' (a 3° cation forms in seconds — the Lucas test)') : 'A primary carbon cannot make a cation, so the halide ion pushes water out from the back in one step — Sₙ2; concentrated acid and heat are needed (NaBr / H₂SO₄ for the bromide, ZnCl₂ for the chloride)') + '. Overall the OH is replaced by ' + X + '.',
        stereoHow: s.mech === 'SN1' ? 'The flat carbocation is attacked from either face: a stereocentre at that carbon is racemised.' : 'Backside attack inverts the configuration at that carbon.' }];
    },
    gate() { return null; },
    condNotes(ruleId, env) {
      const n = ['H–X on an alcohol: the acid protonates the O–H first, then water leaves. Concentrated acid, no water to speak of.'];
      if (env.T < 313) n.push('At ' + degC(env.T) + ' °C only a tertiary alcohol reacts quickly; primary and secondary alcohols need heat (reflux with NaBr / H₂SO₄).');
      return n;
    },
    byproducts: () => [{ smiles: 'O', name: 'water (the leaving group)' }],
    stages(g0, s, out, ruleId, env) {
      const st = [];
      const g1 = clone(g0); atom(g1, s.o.id).charge = 1; st.push({ label: 'protonated alcohol (R–OH₂⁺): water is now the leaving group', graph: g1 });
      if (s.mech === 'SN1') { const g2 = clone(g0); C().removeAtom(g2, s.o.id); atom(g2, s.c.id).charge = 1; st.push({ label: 'carbocation (water has left)', graph: g2 }); }
      return st;
    },
  };
  /* --- NaOH on an ESTER: saponification — hydroxide adds to C=O, the alkoxide leaves, the acid is
     deprotonated at once (irreversible); acidifying the workup gives the carboxylic acid ---
     methyl salicylate + NaOH → salicylic acid (after H₃O⁺) + methanol */
  F.saponify = {
    rules: [], needs: 'ester',
    sites(g, ruleId, env) { return ruleId === 'naoh' ? esterSites(g) : []; },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      g.bonds = g.bonds.filter(b => !((b.a === s.c.id && b.b === s.eo.id) || (b.b === s.c.id && b.a === s.eo.id)));
      const oh = C().addAtomNear(g, 'O', atom(g, s.c.id));
      let alcohol = splitOff(g, s.c.id);
      let main = g, changed = new Set([s.c.id, s.o.id, oh.id]), added = [oh.id], mainIs = 'acid';
      if (alcohol.atoms.length > g.atoms.length) { main = alcohol; alcohol = g; changed = new Set([s.eo.id]); added = []; mainIs = 'alcohol'; }
      s.other = alcohol; s.otherIs = mainIs === 'acid' ? (s.aryl ? 'phenol' : 'alcohol') : 'acid';
      return [{ graph: main, changed, added, share: 1, mechKind: 'acylsub',
        label: mainIs === 'acid' ? 'carboxylic acid (after acidifying the carboxylate) + ' + (s.aryl ? 'the phenol' : 'the alcohol') : (s.aryl ? 'phenol' : 'alcohol') + ' (ester saponified) + the carboxylic acid',
        why: 'SAPONIFICATION. Hydroxide attacks the ester carbonyl carbon (tetrahedral intermediate), the alkoxide is expelled, and the carboxylic acid formed is deprotonated on the spot by the alkoxide / hydroxide — that last step is irreversible, which is why base hydrolysis goes to completion where acid hydrolysis is an equilibrium. In the flask the product is the CARBOXYLATE salt; acidifying the workup (HCl) gives the free acid drawn here.',
        stereoHow: null }];
    },
    gate(g, ruleId, env) {
      if (!env.solventAuto && !['water', 'protic'].includes(env.solventVal)) return 'Saponification needs water (or aqueous alcohol) — hydroxide has to be in solution.';
      return null;
    },
    condNotes(ruleId, env) {
      const n = ['No C–X here, so hydroxide does its other job: it hydrolyses the ester (saponification). One equivalent of NaOH is CONSUMED — it ends up as the carboxylate salt — so use at least one equivalent, usually more, and heat.'];
      if (env.T < 313) n.push('At ' + degC(env.T) + ' °C this is slow; refluxing in aqueous NaOH / ethanol for an hour is the usual practice.');
      return n;
    },
    byproducts(ruleId, env, site) {
      const out = [];
      if (site && site.other) { let nm = ''; try { nm = nameOf(site.other); } catch (e) {} let smi = null; try { smi = C().toSmiles(site.other); } catch (e) {} out.push({ smiles: smi, name: (nm || 'the ' + site.otherIs) + ' (the other half of the ester)', graph: site.other }); }
      return out;
    },
    stages(g0, s) {
      const st = [];
      const g2 = clone(g0); bondOf(g2, s.c.id, s.o.id).order = 1; atom(g2, s.o.id).charge = -1; C().addAtomNear(g2, 'O', atom(g2, s.c.id)); st.push({ label: 'tetrahedral intermediate (hydroxide has added; the alkoxide O⁻)', graph: g2 });
      return st;
    },
  };

  /* --- the haloform reaction: NaOCl (or Br₂ / I₂ + NaOH) on a METHYL KETONE — the CH₃ is trihalogenated,
     hydroxide cleaves the C–C bond, CHX₃ leaves: acetophenone + bleach → benzoic acid + chloroform --- */
  F.haloform = {
    rules: [], needs: 'carbonyl',
    sites(g, ruleId, env) {
      if (ruleId !== 'naocl') return [];
      return carbonylSites(g).map(s => {
        const me = nb(g, s.c.id).find(n => n.atom.element === 'C' && n.bond.order === 1 && hOn(g, n.atom) === 3);
        if (s.cls !== 'ketone' && s.cls !== 'aldehyde') return Object.assign(s, { inert: 'not a ketone / aldehyde' });
        if (!me) return Object.assign(s, { inert: 'no CH₃ on the carbonyl — only a METHYL ketone gives the haloform reaction (the CH₃ must be trihalogenated to make a leaving group)' });
        s.me = me.atom; s.score = 1; s.reason = 'a methyl ketone: the CH₃ next to the C=O is trihalogenated, then cleaved off as CHCl₃';
        s.label = 'methyl ketone C' + s.c.id + '(=O)–CH₃';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      C().removeAtom(g, s.me.id);
      const oh = C().addAtomNear(g, 'O', atom(g, s.c.id));
      return [{ graph: g, changed: new Set([s.c.id, oh.id]), added: [oh.id], share: 1, mechKind: 'acylsub',
        label: 'carboxylic acid (haloform reaction) + chloroform',
        why: 'THE HALOFORM REACTION (the iodoform test run preparatively with bleach). Hydroxide makes the enolate; the enolate takes Cl from HOCl; each chlorine makes the remaining α-hydrogens MORE acidic, so it does not stop until the CH₃ is a CCl₃. Now hydroxide adds to the C=O and the tetrahedral intermediate expels ⁻CCl₃ — a leaving group only because three chlorines stabilise it — which is protonated to chloroform. One carbon shorter: the methyl ketone has become the carboxylic acid (as its salt; acidify to get the acid).',
        stereoHow: null }];
    },
    gate(g, ruleId, env) { if (!env.solventAuto && !['water', 'protic'].includes(env.solventVal)) return 'The haloform reaction needs aqueous hypochlorite / hydroxide.'; return null; },
    condNotes: () => ['No alcohol to oxidise here, so the hypochlorite does its other classic job on the methyl ketone: the haloform reaction. Three equivalents of NaOCl are consumed (one per chlorine) plus one of base; keep it warm (50–60 °C) and work up with acid to free the carboxylic acid.'],
    byproducts: () => [{ smiles: 'ClC(Cl)Cl', name: 'chloroform (CHCl₃ — the haloform)' }, { smiles: '[Na+].[Cl-]', name: 'NaCl' }],
    stages(g0, s) {
      const st = [];
      const g1 = clone(g0); const me = atom(g1, s.me.id); [0, 1, 2].forEach(() => C().addAtomNear(g1, 'Cl', me)); try { C().layoutGraph(g1); } catch (e) {}
      st.push({ label: 'trichloromethyl ketone (the CH₃ has been chlorinated three times through the enolate)', graph: g1 });
      const g2 = clone(g1); bondOf(g2, s.c.id, s.o.id).order = 1; atom(g2, s.o.id).charge = -1; C().addAtomNear(g2, 'O', atom(g2, s.c.id)); try { C().layoutGraph(g2); } catch (e) {}
      st.push({ label: 'tetrahedral intermediate (hydroxide has added) — ⁻CCl₃ is about to leave', graph: g2 });
      return st;
    },
  };

  /* --- Cannizzaro: an aldehyde with NO α-hydrogen + concentrated NaOH disproportionates —
     one molecule is reduced to the alcohol, another oxidised to the carboxylate --- */
  F.amideHyd = {
    rules: [], needs: 'amide',
    sites(g, ruleId, env) {
      const out = [];
      /* Hydrolysis needs WATER or hydroxide. tert-Butoxide and ethoxide are bases, not
         hydroxide: tert-butoxide is far too bulky to add to an amide carbonyl, and neither
         brings any water with it. Cyanide is not a hydrolysing agent at all. Without this
         guard, KOtBu used as the base of a cross-coupling was read as hydrolysing an amide
         somewhere else in the molecule. */
      if (ruleId === 'tbuok' || ruleId === 'naoet' || ruleId === 'nacn') return out;
      g.atoms.forEach(a => {
        if (a.element !== 'C') return;
        const ns = nb(g, a.id);
        const dO = ns.find(n => n.atom.element === 'O' && n.bond.order === 2);
        const nAt = ns.find(n => n.atom.element === 'N' && n.bond.order === 1);
        if (!dO || !nAt) return;
        if (ns.some(n => n.atom.element === 'O' && n.bond.order === 1)) return;      // a carbamate, not an amide
        const nRing = (() => { try { return C().ringAtomsOf(g).indexOf(nAt.atom.id) >= 0 && C().ringAtomsOf(g).indexOf(a.id) >= 0; } catch (e) { return false; } })();
        out.push({ kind: 'amide', c: a, o: dO.atom, n: nAt.atom, lactam: nRing,
          score: (hOn(g, nAt.atom) >= 1 ? 0.4 : 0) - (nRing ? 0.3 : 0),
          label: 'amide (C' + a.id + '=O, N' + nAt.atom.id + ')',
          reason: 'the amide — hydrolysed, but slowly: the nitrogen lone pair donates into the C=O and makes it a poor electrophile, and the leaving group would be an amide anion. This is precisely why an amine is protected AS an amide, and why removing it needs hours of reflux in 6 M acid or base rather than a room-temperature stir' + (nRing ? '; and being a lactam, the ring opens' : '') });
      });
      return out;
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      const basic = ruleId === 'naoh' || ruleId === 'naoet';
      /* the C–N bond breaks: the acid (or its carboxylate) and the amine come apart */
      const amine = cutFragment(g, s.c.id, s.n.id);
      const o = C().addAtomNear(g, 'O', atom(g, s.c.id));
      if (basic) { const oa = atom(g, o.id); if (oa) oa.charge = -1; }
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.c.id, o.id]), added: [o.id], share: 1,
        label: basic ? 'carboxylate + the free amine' : 'carboxylic acid + the amine salt',
        why: 'Amide hydrolysis. Hydroxide (or water, with acid) adds to the carbonyl carbon; the tetrahedral intermediate then has to expel either OH⁻ or the amide anion, and the amide anion is the WORSE leaving group of the two — which is why most of the intermediate simply falls back to the amide and the reaction is slow. ' +
          (basic ? 'What drives it to completion in base is the last step, not the first: once the acid is released it is deprotonated to the carboxylate, and that is irreversible. So base hydrolysis needs a full equivalent of hydroxide, not a catalytic amount. '
                 : 'In acid it is the amine that is taken out of play, as its ammonium salt — the same trick from the other end. ') +
          'Either way: hours at reflux in 6 M acid or base. An amide that came off in five minutes at room temperature would be no use as a protecting group.',
        stereoHow: null, _byover: [{ graph: bmTidy(amine), name: basic ? 'the free amine' : 'the amine, as its ammonium salt' }] }];
    },
    byproducts: () => [],
    condNotes(ruleId, env) { return ['6 M HCl or 6 M NaOH, reflux several hours — an amide is far more robust than an ester and needs the harsher conditions.',
      'A milder alternative when the rest of the molecule will not stand it: an anilide comes off with acid in a few hours at 100 °C, and a trifluoroacetamide comes off with mild base at room temperature (the three fluorines make its carbonyl much more electrophilic).'] },
    stages: () => null,
  };

  F.cannizzaro = {
    rules: [], needs: 'carbonyl',
    sites(g, ruleId, env) {
      if (ruleId !== 'naoh') return [];
      return carbonylSites(g).map(s => {
        if (s.cls !== 'aldehyde') return Object.assign(s, { inert: 'not an aldehyde' });
        const alpha = nb(g, s.c.id).find(n => n.atom.element === 'C');
        if (alpha && hOn(g, alpha.atom) > 0 && !isAromatic(g, alpha.atom.id)) return Object.assign(s, { inert: 'this aldehyde HAS α-hydrogens — with hydroxide it goes into the aldol reaction, not the Cannizzaro' });
        s.score = 1; s.reason = 'an aldehyde with no α-hydrogen (no enolate possible): hydride transfer between two molecules — the Cannizzaro reaction'; s.label = 'aldehyde C' + s.c.id + '=O';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const red = clone(g0); bondOf(red, s.c.id, s.o.id).order = 1;
      const ox = clone(g0); C().addAtomNear(ox, 'O', atom(ox, s.c.id));
      const why = 'CANNIZZARO. Hydroxide adds to one aldehyde to give a tetrahedral alkoxide; that alkoxide hands its C–H hydrogen — as HYDRIDE — to the C=O of a second aldehyde. The donor becomes the carboxylic acid (deprotonated at once to the carboxylate) and the acceptor becomes the alkoxide, protonated to the primary alcohol at the workup. Half the aldehyde is oxidised, half reduced — a disproportionation, possible only because there are no α-hydrogens to make an enolate.';
      return [
        { graph: red, changed: new Set([s.c.id, s.o.id]), added: [], share: 0.5, mechKind: 'red', label: 'primary alcohol (the reduced half)', why, stereoHow: null },
        { graph: ox, changed: new Set([s.c.id]), added: [], share: 0.5, mechKind: 'ox', label: 'carboxylic acid (the oxidised half — as its salt until acidified)', why, stereoHow: null },
      ];
    },
    gate() { return null; },
    condNotes: () => ['No C–X and no ester, but an aldehyde without α-hydrogens: concentrated NaOH (50 %) makes it disproportionate — the Cannizzaro reaction. Expect a 1:1 mixture of the alcohol and the acid; separate them by extracting the alcohol from the basic solution, then acidifying to precipitate the acid.'],
    byproducts: () => [],
    stages(g0, s) {
      const g1 = clone(g0); bondOf(g1, s.c.id, s.o.id).order = 1; atom(g1, s.o.id).charge = -1; C().addAtomNear(g1, 'O', atom(g1, s.c.id)); try { C().layoutGraph(g1); } catch (e) {}
      return [{ label: 'tetrahedral alkoxide (hydroxide has added to one aldehyde) — the hydride donor', graph: g1 }];
    },
  };
  /* --- Étard oxidation: CrO₂Cl₂ takes a benzylic CH₃ / CH₂ to the aldehyde / ketone through the Étard complex --- */
  function benzylicSites(g) {
    const out = [];
    g.atoms.forEach(a => {
      if (a.element !== 'C' || isAromatic(g, a.id)) return;
      const ns = nb(g, a.id);
      const ring = ns.find(n => n.atom.element === 'C' && isAromatic(g, n.atom.id)); if (!ring) return;
      if (ns.some(n => n.bond.order > 1)) return;                                   // vinyl / carbonyl carbons are not benzylic C–H
      const h = hOn(g, a);
      const het = ns.find(n => n.atom.element !== 'C');
      const s = { kind: 'benzylic', atoms: [a.id, ring.atom.id], c: a, ring: ring.atom, h, label: (h === 3 ? 'benzylic CH₃' : h === 2 ? 'benzylic CH₂' : h === 1 ? 'benzylic CH' : 'quaternary benzylic C') + ' (C' + a.id + ')' };
      if (het) s.inert = 'this benzylic carbon already carries a heteroatom (' + het.atom.element + ') — that is a different functional group, not the alkyl C–H the Étard reagent attacks';
      else if (h === 0) s.inert = 'no benzylic C–H to abstract';
      else if (h === 1) s.inert = 'only one benzylic hydrogen: the Étard complex would hydrolyse to a tertiary alcohol, not a carbonyl — the reaction is not synthetically useful here';
      else { s.score = h === 3 ? 1 : 0.7; s.reason = h === 3 ? 'a methyl on the ring — the classic Étard substrate (toluene → benzaldehyde)' : 'a benzylic CH₂ — the hydrogens are weak (the radical/cation is delocalised into the ring), so chromyl chloride takes both and leaves a ketone'; }
      out.push(s);
    });
    return out;
  }
  F.etard = {
    rules: ['etard'], needs: 'benzylic',
    sites(g, ruleId, env) { return benzylicSites(g); },
    variants(g0, s, ruleId, env) {
      const g = clone(g0); const c = atom(g, s.c.id);
      /* the Étard complex, drawn before it is hydrolysed: ArCH(O–CrCl₂OH)₂ — two chromium(V) esters on the same carbon */
      const cx = clone(g0); const cc = atom(cx, s.c.id); const cxIds = [];
      for (let i = 0; i < 2; i++) {
        const o = C().addAtomNear(cx, 'O', cc); const cr = C().addAtomNear(cx, 'Cr', o); const cl1 = C().addAtomNear(cx, 'Cl', cr), cl2 = C().addAtomNear(cx, 'Cl', cr), oh = C().addAtomNear(cx, 'O', cr);
        cxIds.push(o.id, cr.id, cl1.id, cl2.id, oh.id);
      }
      try { C().layoutGraph(cx); } catch (e) {}
      const o2 = C().addAtomNear(g, 'O', c); bondOf(g, c.id, o2.id).order = 2;
      const what = s.h === 3 ? 'aromatic aldehyde' : 'aryl ketone';
      return [{ graph: g, changed: new Set([s.c.id, o2.id]), added: [o2.id], share: 1, label: what + ' (Étard oxidation)',
        why: 'Chromyl chloride (CrO₂Cl₂, 2 equivalents, in CS₂ or CCl₄ at 0 °C) attacks the weak benzylic C–H: an ene-type transfer of hydrogen to a Cr=O gives a chromium(V) ester, and a second CrO₂Cl₂ adds to the other benzylic hydrogen. The result is the ÉTARD COMPLEX — a brown, insoluble solid with two O–Cr groups on the same carbon — which is filtered off and hydrolysed with water to release the ' + (s.h === 3 ? 'aldehyde' : 'ketone') + '. Because the carbonyl is only formed at the hydrolysis, it is never exposed to the oxidant, so it is NOT over-oxidised to the acid (the reason this beats KMnO₄ for making ArCHO).',
        stereoHow: null,
        _stages: [{ label: 'Étard complex — ArCH(O–CrCl₂OH)₂ (brown solid)', graph: cx, changed: new Set([s.c.id].concat(cxIds)) }] }];
    },
    gate(g, ruleId, env) {
      if (env.solventVal === 'water' && !env.solventAuto) return 'Chromyl chloride is hydrolysed violently by water. The Étard reaction is run in a dry non-polar solvent (CS₂ or CCl₄); water is added only afterwards, to break up the Étard complex.';
      return null;
    },
    condNotes(ruleId, env) {
      const n = ['Two equivalents of CrO₂Cl₂ per benzylic carbon: both are needed to build the Étard complex.'];
      if (env.T > 300) n.push('Keep it cold (0 °C): warm chromyl chloride chlorinates the ring and the side chain as well.');
      if (env.excess) n.push('Excess reagent: every benzylic CH₃/CH₂ in the molecule is oxidised.'); else n.push('With one benzylic site the product is clean; with several, the methyl reacts first.');
      return n;
    },
    byproducts: () => [{ smiles: null, name: 'Cr(III) chloride / hydroxide (green, water-soluble) — from hydrolysing the Étard complex' }],
    stages(g0, s, out) { return out._stages || null; },
  };
  /* --- electrophilic aromatic substitution: sulfonation, nitration, bromination --- */
  const EAS = {
    sulfonation: { E: 'SO₃H', name: 'sulfonation', reagent: 'SO₃ / H₂SO₄', el: 'SO3', steric: 0.6, by: [{ smiles: 'O', name: 'water' }] },
    h2so4_dehydrate: { E: 'SO₃H', name: 'sulfonation', reagent: 'conc. H₂SO₄', el: 'SO3', steric: 0.6, by: [{ smiles: 'O', name: 'water' }] },
    nitration: { E: 'NO₂', name: 'nitration', reagent: 'HNO₃ / H₂SO₄', el: 'NO2+', steric: 0.15, by: [{ smiles: 'O', name: 'water' }] },
    ar_br2: { E: 'Br', name: 'bromination', reagent: 'Br₂ / FeBr₃', el: 'Br+', steric: 0.5, by: [{ smiles: 'Br', name: 'HBr' }] },
    ar_cl2: { E: 'Cl', name: 'chlorination', reagent: 'Cl₂ / FeCl₃ (or AlCl₃)', el: 'Cl+', steric: 0.45, by: [{ smiles: 'Cl', name: 'HCl' }] },
    ar_i2: { E: 'I', name: 'iodination', reagent: 'I₂ / HNO₃ (an oxidant is needed)', el: 'I+', steric: 0.6, by: [{ smiles: 'I', name: 'HI' }] },
  };
  /* how a group on the ring steers the next electrophile: v > 0 activates, v < 0 deactivates;
     dir 'op' = ortho/para director, 'm' = meta director */
  function directorOf(g, ringId, subId, ringSet) {
    const a = atom(g, subId), r = atom(g, ringId);
    const ns = nb(g, subId).filter(n => n.atom.id !== ringId);
    const el = a.element;
    if (el === 'O') {
      if (ns.some(n => n.atom.element === 'C' && nb(g, n.atom.id).some(m => m.atom.element === 'O' && m.bond.order === 2))) return { v: 1.0, dir: 'op', word: 'acyloxy (an ester O — its lone pair is shared with the C=O, so only a moderate activator)' };
      return { v: 2.6, dir: 'op', hetero: true, word: hOn(g, a) ? 'hydroxy (a strong activator — its lone pair feeds the ring)' : 'alkoxy (a strong activator — its lone pair feeds the ring)' };
    }
    if (el === 'N') {
      if (a.charge === 1) return { v: -3.2, dir: 'm', word: 'nitro (a strong deactivator, meta director)' };
      if (ns.some(n => n.atom.element === 'C' && nb(g, n.atom.id).some(m => m.atom.element === 'O' && m.bond.order === 2))) return { v: 1.4, dir: 'op', hetero: true, word: 'amido (a moderate activator)' };
      if (nb(g, subId).some(n => n.bond.order === 3)) return { v: -1.6, dir: 'm', word: 'cyano (a deactivator, meta director)' };
      return { v: 3.0, dir: 'op', hetero: true, word: 'amino (a very strong activator)' };
    }
    if (['F', 'Cl', 'Br', 'I'].includes(el)) return { v: -0.7, dir: 'op', hetero: true, word: 'halo (a deactivator that still directs ortho/para — its lone pair stabilises the σ-complex)' };
    if (el === 'S') return { v: -2.0, dir: 'm', word: 'sulfonyl (a deactivator, meta director)' };
    if (el === 'C') {
      const dblO = nb(g, subId).some(n => n.atom.element === 'O' && n.bond.order === 2);
      const tplN = nb(g, subId).some(n => n.atom.element === 'N' && n.bond.order === 3);
      if (dblO) return { v: -1.6, dir: 'm', word: 'carbonyl (an electron-withdrawing deactivator, meta director)' };
      if (tplN) return { v: -1.6, dir: 'm', word: 'cyano (a deactivator, meta director)' };
      const halos = ns.filter(n => ['F', 'Cl', 'Br'].includes(n.atom.element)).length;
      if (halos >= 3) return { v: -1.8, dir: 'm', word: 'trihalomethyl (a deactivator, meta director)' };
      if (isAromatic(g, subId)) return { v: 0.5, dir: 'op', word: 'aryl (a weak activator, ortho/para director by resonance)' };
      if (nb(g, subId).some(n => n.bond.order > 1)) return { v: 0.4, dir: 'op', word: 'alkenyl (a weak activator)' };
      return { v: 0.65, dir: 'op', word: 'alkyl (a weak activator — hyperconjugation feeds the ring; ortho/para director)', bulky: ns.length >= 2 };
    }
    return { v: 0, dir: 'op', word: el };
  }
  /* =====================================================================
     WHERE THE ELECTROPHILE GOES ON A FUSED ARENE

     On a single ring, ortho / meta / para against each substituent decides it.
     On a fused system that reasoning has nothing to say: naphthalene carries no
     substituent at all, yet bromination gives the 1-isomer nine times out of ten.
     What decides it is the σ-complex — whether the positive charge can be parked
     somewhere good while EVERY OTHER RING keeps its aromatic sextet.

     So this counts exactly that. Remove the attacked carbon, then pair the
     remaining atoms into π bonds, leaving either one carbon over to carry the
     charge or one pyrrole-type heteroatom pulled into a π bond (an iminium or
     oxocarbenium). A form is "sextet-preserving" when every ring that does not
     contain the attacked atom has all of its own atoms paired inside itself; a
     charge position reachable that way counts in full, one reachable only by
     breaking another ring counts for little. Each position then scores by how
     well that atom holds a positive charge.

       naphthalene   attack C1 → the charge can sit on C2 or C4 with the other
                     ring intact (two places); attack C2 → only C1 (one place),
                     which is the 9 : 1 α : β ratio
       indole        attack C3 → N1 takes the charge as an iminium, worth more
                     than any carbon; attack C2 → no sextet-preserving form at all
       quinoline     in the pyridine ring the charge cannot avoid the nitrogen, so
                     the benzo ring is attacked, at its own α positions C5 and C8
       2-naphthol    attack C1 → the charge can sit on the carbon carrying the OH

     Written for fused systems only. A single ring keeps the ortho/meta/para
     scoring below, which is tuned against measured isomer ratios. */
  function sigmaComplex(g, sys, rings) {
    const set = new Set(sys);
    const elOf = id => atom(g, id).element;
    const ringsWith = new Map(sys.map(id => [id, rings.filter(r => r.indexOf(id) >= 0)]));
    /* a heteroatom whose lone pair is part of the aromatic cloud (pyrrole, furan,
       the N–H of indole) may either sit out of the π bonding or be pulled into it
       as an onium; a pyridine-type N holds its lone pair in the ring plane and so
       must take a π bond — if it cannot get one, it is left holding the charge */
    const donates = id => {
      const e = elOf(id);
      if (e === 'C') return false;
      if (e === 'N') return nb(g, id).length >= 3 || hOn(g, atom(g, id)) >= 1;
      return (ringsWith.get(id) || []).some(r => r.length === 5);
    };
    const exo = id => nb(g, id).filter(n => !set.has(n.atom.id));
    /* –OH, –OR, –NH2, –NHC(=O)R all put a lone pair onto the ring; –NO2 does not, even though
       its nitrogen is single-bonded to the ring. Reading it as a donor is what made the app
       nitrate 1-nitronaphthalene at C4 instead of on the far ring. */
    const donorAt = n => (n.atom.element === 'O' || n.atom.element === 'N') && n.bond.order === 1 &&
      !n.atom.charge && !nb(g, n.atom.id).some(m => m.bond.order > 1);
    const weightOf = id => {
      const e = elOf(id);
      if (e !== 'C') return donates(id) ? (e === 'N' ? 3.5 : e === 'O' ? 1.9 : 1.5) : -2.5;
      const out = exo(id);
      if (out.some(donorAt)) return 3.0;
      if (out.some(n => n.bond.order > 1 || (n.atom.element === 'N' && n.atom.charge > 0))) return 0.15;
      if (out.some(n => ['F', 'Cl', 'Br', 'I'].includes(n.atom.element))) return 0.8;
      if (out.some(n => n.atom.element === 'C')) return 1.6;          // alkyl: hyperconjugation
      return 1.0;
    };
    const wordFor = id => {
      const e = elOf(id);
      if (e !== 'C') return donates(id) ? 'the ring ' + e + ' as an ' + (e === 'N' ? 'iminium' : 'oxocarbenium') : 'the ring ' + e + ' (which is why this position is so bad)';
      const out = exo(id);
      if (out.some(donorAt)) return 'the carbon carrying the ' + (out.find(n => n.atom.element === 'O' && donorAt(n)) ? '–O' : '–N') + ' group';
      if (out.some(n => n.bond.order > 1)) return 'a carbon carrying an electron-withdrawing group';
      return 'a ring carbon';
    };

    /* every place the charge can go, for an attack at atom i */
    function chargeSites(i) {
      const list = sys.filter(x => x !== i);
      const adj = new Map(list.map(x => [x, nb(g, x).filter(n => set.has(n.atom.id) && n.atom.id !== i).map(n => n.atom.id)]));
      const other = rings.filter(r => r.indexOf(i) < 0);
      const partner = new Map();                 // atom -> its π-bond partner, or null for a lone pair
      const found = new Map();                   // charge position -> how much of the rest survives
      /* A ring survives when every one of its atoms is paired inside it. Breaking a benzene ring
         costs almost everything; breaking a five-membered heteroaromatic ring costs much less,
         and is exactly what the useful σ-complexes do — the iminium form of the indole C3 adduct
         has no aromatic five-ring at all, and the N of carbazole reaches C3 and C6 that way. */
      const done = (carrier) => {
        if (carrier === null) return;            // no charge anywhere: not a σ-complex
        const intact = r => r.every(x => {
          if (x === carrier) return false;
          const p = partner.get(x);
          return p === null || (p !== undefined && r.indexOf(p) >= 0);
        });
        const broken = other.filter(r => !intact(r));
        /* a five-membered heteroaromatic ring is much less aromatic than benzene, so a form that
           sacrifices it is a real contributor — and more so when what it buys is an iminium on a
           pyrrole-type nitrogen, which is the strongest stabilisation available. Sacrificing a
           benzene ring (a quinoidal form) is worth very little. */
        const f = !broken.length ? 1 : !broken.every(r => r.length === 5) ? 0.12
          : (elOf(carrier) === 'N' && donates(carrier) ? 0.5 : 0.2);
        if (!found.has(carrier) || found.get(carrier).f < f) found.set(carrier, { f, partner: new Map(partner) });
      };
      const rec = (k, carrier) => {
        while (k < list.length && partner.has(list[k])) k++;
        if (k === list.length) { done(carrier); return; }
        const a = list[k];
        for (const b of adj.get(a)) {
          if (partner.has(b)) continue;
          const both = donates(a) && donates(b);
          if (both) continue;                    // two oniums would be a dication
          const onium = donates(a) ? a : donates(b) ? b : null;
          if (onium !== null && carrier !== null) continue;
          partner.set(a, b); partner.set(b, a);
          rec(k + 1, onium === null ? carrier : onium);
          partner.delete(a); partner.delete(b);
        }
        if (donates(a)) { partner.set(a, null); rec(k + 1, carrier); partner.delete(a); }
        else if (carrier === null) { partner.set(a, null); rec(k + 1, a); partner.delete(a); }
      };
      rec(0, null);
      return found;
    }

    /* Take this atom and one other atom of the same ring away: is what remains nothing but whole
       benzene rings? That is what makes the 9,10-positions of anthracene and phenanthrene
       special. Those two carbons carry the whole of the ring system's "extra" π system between
       them; addition across them (or substitution at one, with the charge on the other) leaves
       two complete aromatic rings either side, which is why both react there and not at an α
       position. Naphthalene has no such pair, which is why it has no meso position. */
    function mesoBond(i) {
      const six = rings.filter(r => r.length === 6);
      const mates = [...new Set(rings.filter(r => r.indexOf(i) >= 0).flatMap(r => r))].filter(x => x !== i);
      for (const j of mates) {
        const gone = new Set([i, j]);
        const rest = sys.filter(x => !gone.has(x));
        if (!rest.length) continue;
        const keep = six.filter(r => r.every(x => !gone.has(x)));
        if (keep.length < 2) continue;
        const covered = new Set(); keep.forEach(r => r.forEach(x => covered.add(x)));
        if (rest.every(x => covered.has(x))) return j;
      }
      return null;
    }

    return function score(i) {
      const found = chargeSites(i);
      if (!found.size) return null;
      let sc = 0; const clean = [], dirty = [];
      found.forEach((rec, p) => {
        /* a charge position that only exists by breaking a benzene ring is a real but minor
           contributor — the quinoidal forms. Counting those at anything like full weight made
           8-nitroquinoline beat 5-nitroquinoline tenfold; the two form about equally. */
        sc += weightOf(p) * rec.f;
        (rec.f === 1 ? clean : dirty).push(p);
      });
      const meso = mesoBond(i);
      if (meso !== null) sc += 1.6;
      const words = [...new Set(clean.map(wordFor))];
      const nOther = rings.length - (ringsWith.get(i) || []).length;
      let why;
      if (!clean.length) why = 'the σ-complex here cannot keep the other ring aromatic, whichever way the charge is drawn';
      else why = 'the σ-complex can put the + on ' + (clean.length === 1 ? words.join(' or ') : clean.length + ' positions (' + words.join(', ') + ')') +
        (nOther > 0 ? ' while the other ring' + (nOther > 1 ? 's stay' : ' stays') + ' aromatic' : '');
      if (meso !== null) why = 'a meso position — this carbon and one other hold the ring system\u2019s extra π system between them, so reacting here leaves two complete benzene rings — and ' + why;
      return { score: sc, why, clean: clean.length, meso: meso !== null };
    };
  }

  function areneSites(g, ruleId, env) {
    const spec = EAS[ruleId]; if (!spec) return [];
    let systems = [];
    try { systems = C().ringSystemsOf(g, g.atoms.filter(a => isAromatic(g, a.id)).map(a => a.id)); } catch (e) { systems = []; }
    /* a ring heteroatom is not always flagged aromatic (the N of indole, of quinoline), and
       without it the five-membered ring does not close, so the system reads as one ring and
       the σ-complex analysis never runs. Put back any ring heteroatom bridging the system —
       heteroatoms only, so the saturated carbons of tetralin stay out of it. */
    systems = systems.map(sys => {
      const set = new Set(sys);
      const add = g.atoms.filter(a => a.element !== 'C' && !set.has(a.id) &&
        nb(g, a.id).filter(n => set.has(n.atom.id)).length >= 2 &&
        nb(g, a.id).every(n => set.has(n.atom.id) || n.bond.order === 1)).map(a => a.id);
      return add.length ? sys.concat(add) : sys;
    });
    const out = [];
    systems.forEach(sys => {
      const set = new Set(sys);
      /* ring distances (ortho 1, meta 2, para 3) inside the system */
      const dist = (a, b) => { const seen = new Map([[a, 0]]); const q = [a]; while (q.length) { const id = q.shift(); if (id === b) return seen.get(id); nb(g, id).forEach(n => { if (set.has(n.atom.id) && !seen.has(n.atom.id)) { seen.set(n.atom.id, seen.get(id) + 1); q.push(n.atom.id); } }); } return 9; };
      const hetero = sys.filter(id => atom(g, id).element !== 'C');
      const subs = [];
      sys.forEach(id => nb(g, id).forEach(n => { if (!set.has(n.atom.id) && atom(g, id).element === 'C') subs.push({ ring: id, sub: n.atom.id, d: directorOf(g, id, n.atom.id, set) }); }));
      /* pyridine-type N: two ring bonds and NO hydrogen or substituent, so its lone pair is in
         the ring plane and the ring is deactivated. The N–H of pyrrole and indole also has two
         ring bonds but puts its lone pair into the π cloud, which activates the ring instead —
         counting it as a pyridine is what once made the app brominate pyrrole at C3. */
      const ringDeact = hetero.some(id => atom(g, id).element === 'N' && nb(g, id).length === 2 && hOn(g, atom(g, id)) === 0);
      let rings = [];
      try { rings = C().ringsOfSystem(g, sys); } catch (e) { rings = []; }
      const sigma = rings.length > 1 ? sigmaComplex(g, sys, rings) : null;
      sys.forEach(id => {
        const a = atom(g, id);
        if (a.element !== 'C' || hOn(g, a) < 1) return;
        let score = 0; const why = [];
        subs.forEach(su => {
          const r = dist(id, su.ring);
          const rel = r === 1 ? 'ortho' : r === 2 ? 'meta' : r === 3 ? 'para' : null;
          if (!rel) return;
          /* ortho pays a steric price; heteroatom donors and halogens pay more (their
             para isomer is also the better-solvated, thermodynamic one) */
          const st = spec.steric + (su.d.bulky ? 0.5 : 0) + (su.d.hetero ? 0.5 : 0);
          let c;
          if (su.d.dir === 'op' && su.d.v > 0) c = rel === 'para' ? su.d.v : rel === 'ortho' ? su.d.v - st : -0.4;
          else if (su.d.dir === 'op') c = rel === 'para' ? su.d.v : rel === 'ortho' ? su.d.v - st : su.d.v - 1.3;
          else c = rel === 'meta' ? su.d.v + 1.2 : rel === 'ortho' ? su.d.v - 0.3 : su.d.v;
          score += c; why.push(rel + ' to the ' + su.d.word.split(' (')[0]);
        });
        hetero.forEach(h => { const r = dist(id, h); const hel = atom(g, h).element;
          if (ringDeact && hel === 'N') { score += r === 2 ? -1.0 : -2.5; why.push(r === 2 ? 'β to the ring N (the least deactivated position of a pyridine)' : 'α/γ to the ring N (the σ-complex would put + next to N⁺)'); }
          else if (r === 1) { score += 1.2; why.push('next to the ring ' + hel + ' (C2 of a five-membered heteroaromatic — the σ-complex has the most resonance forms)'); } });
        /* on a fused system the σ-complex decides it, not ortho/meta/para */
        let sig = null;
        if (sigma) {
          sig = sigma(id);
          if (!sig) return;
          /* a pyridine-type N deactivates its OWN ring hard (the σ-complex would have to put +
             next to N⁺) and the ring fused to it only mildly — which is why quinoline nitrates
             on the benzo ring, at C5 and C8 equally. Scoring that by distance to the N, as a
             single ring would, made C8 look better than C5 for no reason. */
          const myRings = rings.filter(r => r.indexOf(id) >= 0);
          const deact = hetero.filter(h => atom(g, h).element === 'N' && nb(g, h).length === 2 && hOn(g, atom(g, h)) === 0)
            .reduce((t, h) => t + (myRings.some(r => r.indexOf(h) >= 0) ? -1.6 : -0.5), 0);
          score = sig.score + deact;
          if (deact < 0) why.push(myRings.some(r => r.some(x => atom(g, x).element === 'N' && nb(g, x).length === 2 && hOn(g, atom(g, x)) === 0))
            ? 'this ring carries a pyridine-type N, which deactivates it strongly' : 'the ring fused to it carries a pyridine-type N, a mild deactivation');
          why.length = 0; why.push(sig.why);
          /* Sulfonation is reversible, so hot and left alone it settles on the isomer that is
             most STABLE rather than the one that forms fastest. On naphthalene the 1-sulfonic
             acid is the fast one but its SO3H is crowded against the H at the 8-position (peri
             strain); above about 150 °C it comes off again and the 2-isomer is what is left. */
          const peri = nb(g, id).some(n => set.has(n.atom.id) && rings.filter(r => r.indexOf(n.atom.id) >= 0).length > 1);
          if (ruleId === 'sulfonation' && env.T >= 413 && peri) {
            score -= 1.8;
            why.push('hot enough for the reversible sulfonation to settle on the more stable isomer, and an SO₃H here is crowded against the peri hydrogen');
          }
        }
        const s = { kind: 'arene', atoms: [id], c: a, ringIds: sys, score, label: 'C' + id + ' (ring position)',
          reason: (why.length ? why.join(', ') : 'an unsubstituted ring position — every position is equivalent'), rel: why[0] || '' };
        if (sig) s.sigma = sig;
        s.subs = subs;
        if (ringDeact && hetero.length && score < -1.5 && !subs.some(su => su.d.v > 1)) s.slow = 'pyridine rings are strongly deactivated (N⁺–H under the acid); forcing conditions are needed';
        out.push(s);
      });
    });
    return out;
  }
  F.eas = {
    rules: ['sulfonation', 'nitration', 'ar_br2', 'ar_cl2', 'ar_i2'], needs: 'arene',
    sites(g, ruleId, env) { return areneSites(g, ruleId, env); },
    variants(g0, s, ruleId, env) {
      const spec = EAS[ruleId]; const g = clone(g0); const c = atom(g, s.c.id);
      const added = []; let changed = new Set([c.id]);
      if (spec.el === 'SO3') {
        const S = C().addAtomNear(g, 'S', c); added.push(S.id);
        const o1 = C().addAtomNear(g, 'O', S), o2 = C().addAtomNear(g, 'O', S), o3 = C().addAtomNear(g, 'O', S);
        bondOf(g, S.id, o1.id).order = 2; bondOf(g, S.id, o2.id).order = 2; added.push(o1.id, o2.id, o3.id);
      } else if (spec.el === 'NO2+') {
        const N = C().addAtomNear(g, 'N', c); N.charge = 1; added.push(N.id);
        const o1 = C().addAtomNear(g, 'O', N), o2 = C().addAtomNear(g, 'O', N);
        bondOf(g, N.id, o1.id).order = 2; o2.charge = -1; added.push(o1.id, o2.id);
      } else {
        const X = C().addAtomNear(g, spec.el === 'Cl+' ? 'Cl' : spec.el === 'I+' ? 'I' : 'Br', c); added.push(X.id);
      }
      added.forEach(id => changed.add(id));
      const relWord = s.rel ? s.rel.replace(/ \(.*$/, '') : '';
      const top = (s.subs || []).slice().sort((p, q) => Math.abs(q.d.v) - Math.abs(p.d.v))[0];
      return [{ graph: g, changed, added, share: 1,
        label: (relWord ? relWord.split(' to the ')[0] + '-' : '') + spec.name.replace('ation', 'ated') + ' ring' + (relWord ? ' (' + relWord + ')' : ''),
        why: 'Electrophilic aromatic substitution: the ring’s π electrons attack ' + spec.E.replace('H', '') + (spec.el === 'NO2+' ? '⁺ (nitronium, made by H₂SO₄ protonating HNO₃)' : spec.el === 'Br+' ? ' (Br₂ polarised by FeBr₃ — an electrophilic bromine)' : ' (sulfur trioxide, present in the concentrated acid)') +
          ' to give the σ-complex (arenium ion); losing the H⁺ from that carbon restores the aromatic ring. ' +
          (top ? 'The ' + top.d.word + ' decides the position: ' + s.reason + '.' : 'Every ring position is the same, so there is one product.'),
        stereoHow: null }];
    },
    gate(g, ruleId, env) {
      if (ruleId === 'sulfonation' && env.solventVal === 'water' && !env.solventAuto) return 'Dilute (aqueous) sulfuric acid does not sulfonate — the SO₃ concentration is nil. It needs concentrated or fuming H₂SO₄.';
      return null;
    },
    condNotes(ruleId, env) {
      const n = [];
      if (ruleId === 'sulfonation' || ruleId === 'h2so4_dehydrate') n.push('Sulfonation is reversible: dilute aqueous acid and heat take the SO₃H back off (a way to use it as a temporary blocking group). Hot conditions favour the thermodynamic para isomer.');
      if (ruleId === 'nitration') n.push('Keep it below 50 °C — nitration is exothermic, and an activated ring will be nitrated twice with excess nitrating mixture.');
      if (ruleId === 'ar_br2') n.push('Without the FeBr₃ Lewis acid, Br₂ does not touch a benzene ring (phenols and anilines are the exception — they brominate in plain Br₂ water).');
      if (env.excess) n.push('Excess reagent: an activated ring is substituted again, at the next best position.');
      return n;
    },
    byproducts: (ruleId) => (EAS[ruleId] || EAS.sulfonation).by,
    /* excess electrophile: the major product is substituted once more, at its own
       best position — unless the first group has deactivated the ring */
    excess(g0, ruleId, env, all, live) {
      const top = all.slice().sort((p, q) => q.share - p.share)[0]; if (!top) return;
      const s2 = areneSites(top.graph, ruleId, env).filter(x => !x.inert).sort((p, q) => q.score - p.score);
      if (!s2.length || s2[0].score < -0.6) return;
      const v = F.eas.variants(top.graph, s2[0], ruleId, env)[0]; if (!v) return;
      all.forEach(o => { o.share *= 0.35; });
      all.push({ graph: v.graph, changed: new Set([...top.changed, ...v.changed]), added: [...(top.added || []), ...v.added], share: 0.65, site: live[0], exhaustive: true,
        label: 'substituted twice (excess reagent) — ' + (s2[0].rel || 'second position'),
        why: 'With the electrophile in excess the ring is attacked again. The first group (' + (EAS[ruleId] || EAS.sulfonation).E + ') is ' + (s2[0].score < 0 ? 'deactivating, so the second substitution is slower and goes ' : 'still outweighed by the original director, so the second group goes ') + (s2[0].rel || 'to the next best position') + '.',
        stereoHow: null });
    },
    stages(g0, s, out, ruleId) {
      /* the σ-complex: the new group on the carbon (still holding its H), the ring’s
         aromaticity broken and the + charge on the neighbouring carbon */
      try {
        const g = clone(out.graph); const c = atom(g, s.c.id);
        const dbl = nb(g, c.id).find(n => n.bond.order === 2 && (s.ringIds || []).includes(n.atom.id));
        if (!dbl) return null;
        dbl.bond.order = 1; dbl.atom.charge = 1;
        return [{ label: 'σ-complex (arenium ion) — the ring has lost its aromaticity for a moment', graph: g }];
      } catch (e) { return null; }
    },
  };
  F.dehyd.fallback = 'eas';
  F.hydrobor.fallback = 'amideRed';

  /* =====================================================================
     MOLECULE + MOLECULE — "react with any molecule" (the common-reagents
     setting switched off). The typed reagent is a full structure; the
     engine looks at what the two molecules ARE (an acid and an amine, an
     alkyl halide and a thiol, a diene and a dienophile …) and runs the
     taught bimolecular reaction that fits. The reagent lives on the dynamic
     rule (rule.molecule); the substrate is the drawn molecule.
     ===================================================================== */
  const bmNb = (g, id) => nb(g, id);
  const bmArom = (g, id) => isAromatic(g, id);
  const bmDblO = (g, c) => bmNb(g, c).find(n => n.atom.element === 'O' && n.bond.order === 2 && bmNb(g, n.atom.id).length === 1);
  const bmIsOH = (g, n) => n.atom.element === 'O' && n.bond.order === 1 && hOn(g, n.atom) >= 1 && bmNb(g, n.atom.id).length === 1;
  const bmSp3 = (g, a) => a.element === 'C' && !bmArom(g, a.id) && bmNb(g, a.id).every(n => n.bond.order === 1);
  /* -------------------------------------------------------------------
     Diels-Alder regiochemistry: the ortho-para rule.
     Both numberings of a diene are enumerated as separate matches, so the
     two regiochemical answers both exist already - this is what chooses
     between them. Treat the diene as a polarised nucleophile and the
     dienophile as a polarised electrophile and let the ends that carry the
     most charge pair up. Pushing the lone pair through the diene shows
     which terminus is the nucleophilic one:
       a donor on C1 or C3  ->  C4 is nucleophilic
       a donor on C2 or C4  ->  C1 is nucleophilic
     and on the dienophile the carbon NOT carrying the withdrawing group is
     the electrophilic one. Pair those two and a 1-substituted diene comes
     out "ortho", a 2-substituted diene "para" - which is the rule as it is
     taught, arrived at from the charges rather than memorised.
     ------------------------------------------------------------------- */
  const DA_DONOR = { O: 2.2, N: 2.4, S: 1.8, F: 0.4, Cl: 0.3, Br: 0.3, I: 0.3 };
  const daPullAt = (g, c) => {
    let s = 0;
    bmNb(g, c.id).forEach(n => {
      const e = n.atom.element;
      if (e === 'N' && (n.atom.charge === 1 || n.bond.order === 3)) { s += 2.5; return; }            /* nitro, nitrile */
      if (e === 'C' && (bmDblO(g, n.atom.id) || bmNb(g, n.atom.id).some(m => m.atom.element === 'N' && m.bond.order === 3))) { s += 2.2; return; }
      if (e === 'S' && bmNb(g, n.atom.id).filter(m => m.atom.element === 'O' && m.bond.order === 2).length >= 2) { s += 2.2; return; }
      if (e === 'C' && bmArom(g, n.atom.id)) s += 0.3;
    });
    return s;
  };
  const daDonorAt = (g, c, skip) => {
    let s = 0;
    bmNb(g, c.id).forEach(n => {
      if (skip.indexOf(n.atom.id) >= 0) return;
      const e = n.atom.element;
      if (DA_DONOR[e] !== undefined) {
        if (n.bond.order > 1) { s -= 1.5; return; }                                                  /* already pulled into a pi bond */
        if (e === 'N' && n.atom.charge === 1) { s -= 2.5; return; }
        s += DA_DONOR[e]; return;
      }
      if (e === 'C') {
        if (bmDblO(g, n.atom.id) || bmNb(g, n.atom.id).some(m => m.atom.element === 'N' && m.bond.order === 3)) { s -= 2.0; return; }
        s += bmArom(g, n.atom.id) ? 0.8 : 0.9;
      }
    });
    return s;
  };
  /* > 0 when THIS pairing (c1 to a, c4 to b) is the one that is actually observed */
  function daRegio(X, d, Y, dp) {
    const skip = [d.c1.id, d.c2.id, d.c3.id, d.c4.id];
    const n1 = daDonorAt(X, d.c2, skip) + daDonorAt(X, d.c4, skip);   /* donors that make C1 nucleophilic */
    const n4 = daDonorAt(X, d.c1, skip) + daDonorAt(X, d.c3, skip);   /* donors that make C4 nucleophilic */
    const bias = (n1 - n4) * (daPullAt(Y, dp.b) - daPullAt(Y, dp.a));
    if (!bias) return 0;
    return Math.max(-0.9, Math.min(0.9, bias * 0.25));
  }
  /* how the substituents ended up related round the new ring, for the write-up */
  function daRelation(d, dp, X, Y) {
    const skip = [d.c1.id, d.c2.id, d.c3.id, d.c4.id];
    const on = c => daDonorAt(X, c, skip) > 0.5;
    const pulled = daPullAt(Y, dp.a) > daPullAt(Y, dp.b) ? 'a' : (daPullAt(Y, dp.b) > daPullAt(Y, dp.a) ? 'b' : null);
    if (!pulled) return null;
    /* ring order is c1 c2 c3 c4 b a, so the withdrawing group sits at position 6 (a) or 5 (b) */
    const wPos = pulled === 'a' ? 6 : 5;
    let dPos = null;
    if (on(d.c1)) dPos = 1; else if (on(d.c2)) dPos = 2; else if (on(d.c3)) dPos = 3; else if (on(d.c4)) dPos = 4;
    if (dPos === null) return null;
    let gap = Math.abs(wPos - dPos); if (gap > 3) gap = 6 - gap;
    return { rel: gap === 1 ? 'ortho (1,2)' : (gap === 2 ? 'meta (1,3)' : 'para (1,4)'), dPos: dPos, donorTerminal: dPos === 1 || dPos === 4 };
  }
  /* what a molecule can offer — every handle it carries, by role */
  const METALS = ['Mg', 'Li', 'Na', 'K', 'Zn'];
  /* How badly does this aromatic carbon want a nucleophile? Sₙ Ar goes through a Meisenheimer
     anion, and that anion is only bearable if the negative charge can land on an
     electron-withdrawing group or on a ring nitrogen — which it can only do from ortho or para.
     Null when nothing activates the position, which is why chlorobenzene is inert to amines
     and 2-chloropyridine and 4-chloronitrobenzene are not. */
  function bmSnArActivation(g, id) {
    let sys = null;
    try { sys = C().ringSystemsOf(g, C().ringAtomsOf(g)).find(r => r.indexOf(id) >= 0); } catch (e) { sys = null; }
    if (!sys) return null;
    const set = new Set(sys);
    const dist = (a, b) => { const seen = new Map([[a, 0]]); const q = [a]; while (q.length) { const x = q.shift(); if (x === b) return seen.get(x); bmNb(g, x).forEach(n => { if (set.has(n.atom.id) && !seen.has(n.atom.id)) { seen.set(n.atom.id, seen.get(x) + 1); q.push(n.atom.id); } }); } return 99; };
    const hits = [];
    sys.forEach(x => {
      const d = dist(id, x);
      if (d !== 1 && d !== 3) return;                                  // ortho or para only
      const a = atom(g, x);
      if (a.element === 'N' && bmNb(g, x).length === 2 && hOn(g, a) === 0) { hits.push({ w: d === 1 ? 2.4 : 2.0, word: 'the ring nitrogen ' + (d === 1 ? 'right next to it' : 'para to it') }); return; }
      bmNb(g, x).forEach(n => {
        if (set.has(n.atom.id)) return;
        const e = n.atom.element, where = d === 1 ? 'ortho' : 'para';
        if (e === 'N' && n.atom.charge === 1 && bmNb(g, n.atom.id).some(m => m.atom.element === 'O')) hits.push({ w: 2.5, word: 'the nitro group ' + where + ' to it' });
        else if (e === 'C' && n.bond.order === 1 && bmNb(g, n.atom.id).some(m => m.atom.element === 'N' && m.bond.order === 3)) hits.push({ w: 1.4, word: 'the nitrile ' + where + ' to it' });
        else if (e === 'C' && n.bond.order === 1 && bmDblO(g, n.atom.id)) hits.push({ w: 1.3, word: 'the carbonyl ' + where + ' to it' });
        else if (e === 'S' && bmNb(g, n.atom.id).filter(m => m.atom.element === 'O' && m.bond.order === 2).length === 2) hits.push({ w: 1.2, word: 'the sulfonyl group ' + where + ' to it' });
      });
    });
    if (!hits.length) return null;
    return { w: hits.reduce((t, h) => t + h.w, 0), words: [...new Set(hits.map(h => h.word))] };
  }
  /* how willingly this C–X bond adds to Pd(0): I > OTf ≈ Br ≫ Cl ≫ F, and a vinyl halide
     rather more than an aryl one */
  const bmOxAdd = h => (h.otf ? 0.55 : { I: 0.7, Br: 0.5, Cl: -0.3, F: -2.2 }[h.x.element] !== undefined ? { I: 0.7, Br: 0.5, Cl: -0.3, F: -2.2 }[h.x.element] : 0) + (h.aryl ? 0 : 0.1);
  /* a penalty for doing Sₙ Ar at this halide when a more palladium-reactive one is available
     somewhere else on the same molecule: I > Br >> Cl decides where the chemistry happens */
  const HAL_RANK = { I: 3, Br: 2, Cl: 1, F: 0 };
  function betterHalideElsewhere(X, h) {
    if (h.otf) return 0;
    const mine = HAL_RANK[h.x.element] === undefined ? 0 : HAL_RANK[h.x.element];
    let best = mine;
    X.atoms.forEach(a => {
      if (!bmArom(X, a.id)) return;
      bmNb(X, a.id).forEach(n => {
        const r = HAL_RANK[n.atom.element];
        if (r !== undefined && n.bond.order === 1 && a.id !== h.c.id) best = Math.max(best, r);
      });
    });
    return best > mine ? 1.4 : 0;
  }
  const bmXWord = h => h.otf ? 'triflate' : h.x.element === 'I' ? 'iodide' : h.x.element === 'Br' ? 'bromide' : h.x.element === 'Cl' ? 'chloride' : h.x.element + '-ide';

  function bmProfile(g) {
    const P = { acyl: [], sulfonyl: [], amine: [], alcohol: [], phenol: [], thiol: [], halide: [], carbonyl: [], enol: [], grignard: [], diene: [], dienophile: [], arene: [], carboxylate: [], cyanide: [], azide: [], alkoxide: [],
      cuprate: [], halideSp2: [], enone: [], co2: [], ester: [], nitrile: [], epoxide: [], acetylene: [], strongBase: [], acidicH: [],
      boron: [], stannane: [], azoleCH: [], peracid: [], sulfide: [], oxidN: [], activeCH: [], ylide: [], phosphonate: [],
      biNuc: [], biElec: [], nxs: [], thionate: [], amideC: [], nitrosating: [], metal0: [],
      enolate: [], enamine: [], esterC: [] };
    g.atoms.forEach(a => {
      const ns = bmNb(g, a.id);
      if (a.element === 'C') {
        const dO = bmDblO(g, a.id);
        if (dO) {
          const oh = ns.find(n => bmIsOH(g, n)), cl = ns.find(n => n.atom.element === 'Cl'), oc = ns.find(n => n.atom.element === 'O' && n.bond.order === 1 && !bmIsOH(g, n)), nN = ns.find(n => n.atom.element === 'N' && n.bond.order === 1);
          const cC = ns.filter(n => n.atom.element === 'C');
          if (oh) { P.acyl.push({ c: a, o: dO.atom, lg: oh.atom, kind: 'acid' }); P.acidicH.push({ at: oh.atom, pka: 5, what: 'carboxylic acid O–H' }); }
          else if (cl) P.acyl.push({ c: a, o: dO.atom, lg: cl.atom, kind: 'acyl chloride' });
          else if (oc) { const other = bmNb(g, oc.atom.id).find(n => n.atom.id !== a.id); const anh = other && other.atom.element === 'C' && bmDblO(g, other.atom.id); P.acyl.push({ c: a, o: dO.atom, lg: oc.atom, kind: anh ? 'anhydride' : 'ester', anhC: anh ? other.atom : null }); if (!anh) { P.ester.push({ c: a, o: dO.atom, or: oc.atom }); cC.forEach(n => { if (bmSp3(g, n.atom) && hOn(g, n.atom) >= 1) P.acidicH.push({ at: n.atom, pka: 25, what: 'α-C–H of the ester', enolate: true }); }); } }
          else if (nN) { /* an amide: not an acyl donor here */ }
          else if (a.charge === 0 && !cl && !(ns.length === 2 && ns.every(n => n.atom.element === 'O' && n.bond.order === 2))) {
            const kind = cC.length === 2 ? 'ketone' : 'aldehyde';
            P.carbonyl.push({ c: a, o: dO.atom, kind, conj: cC.some(n => hasPiToCarbon(g, n.atom.id, a.id)) });
            /* an enone: C=C conjugated to the C=O — the β-carbon takes a soft nucleophile */
            cC.forEach(n => { const beta = bmNb(g, n.atom.id).find(m => m.atom.element === 'C' && m.bond.order === 2 && m.atom.id !== a.id && !bmArom(g, m.atom.id)); if (beta && !bmArom(g, n.atom.id)) P.enone.push({ cC: a, o: dO.atom, alpha: n.atom, beta: beta.atom }); });
            /* α-hydrogens make it an enol(ate) nucleophile */
            cC.forEach(n => { if (bmSp3(g, n.atom) && hOn(g, n.atom) >= 1) P.enol.push({ c: a, o: dO.atom, alpha: n.atom, kind }); });
            cC.forEach(n => { if (bmSp3(g, n.atom) && hOn(g, n.atom) >= 1) P.acidicH.push({ at: n.atom, pka: 20, what: 'α-C–H of the ' + kind, enolate: true }); });
          }
          if (a.charge === -1 && !oh) { /* skip */ }
          if (ns.find(n => n.atom.element === 'O' && n.bond.order === 1 && n.atom.charge === -1)) P.carboxylate.push({ c: a });
        }
        if (a.charge === -1 && ns.some(n => n.atom.element === 'N' && n.bond.order === 3)) P.cyanide.push({ c: a });
        const mg = ns.find(n => METALS.includes(n.atom.element)); if (mg) { const kind = mg.atom.element === 'Mg' ? 'Grignard reagent' : mg.atom.element === 'Li' ? 'organolithium' : mg.atom.element === 'Zn' ? 'organozinc' : 'metal acetylide';
          const acet = ns.some(n => n.bond.order === 3);
          P.grignard.push({ c: a, mg: mg.atom, kind: acet ? 'metal acetylide' : kind, acet });
          if (mg.atom.element === 'Li' && !acet) P.strongBase.push({ c: a, m: mg.atom, kind: 'organolithium', conj: 'the alkane', maxPka: 45 });
          else if (mg.atom.element === 'Mg' && !acet) P.strongBase.push({ c: a, m: mg.atom, kind: 'Grignard reagent', conj: 'the alkane', maxPka: 32 });
          else if (acet && (mg.atom.element === 'Na' || mg.atom.element === 'Li')) P.strongBase.push({ c: a, m: mg.atom, kind: 'metal acetylide', conj: 'the terminal alkyne', maxPka: 22 }); }
        /* the transmetalation partners of a cross-coupling: a carbon on boron (boronic acid,
           pinacol ester, or a potassium trifluoroborate) or a carbon on tin (a stannane) */
        const B = ns.find(n => n.atom.element === 'B' && n.bond.order === 1);
        if (B) {
          const bn = bmNb(g, B.atom.id);
          const oxy = bn.filter(n => n.atom.element === 'O').length;
          const flu = bn.filter(n => n.atom.element === 'F').length;
          const esterified = bn.some(n => n.atom.element === 'O' && bmNb(g, n.atom.id).some(m => m.atom.element === 'C'));
          if (oxy >= 2 || flu >= 3) P.boron.push({ c: a, b: B.atom, kind: flu >= 3 ? 'trifluoroborate salt' : esterified ? 'boronate ester' : 'boronic acid', aryl: bmArom(g, a.id) });
        }
        const Sn = ns.find(n => n.atom.element === 'Sn' && n.bond.order === 1);
        if (Sn && bmNb(g, Sn.atom.id).filter(n => n.atom.element === 'C').length >= 4) P.stannane.push({ c: a, sn: Sn.atom, aryl: bmArom(g, a.id),
          sp2: bmArom(g, a.id) || ns.some(n => n.atom.element === 'C' && n.bond.order >= 2) });
        const cu = ns.find(n => n.atom.element === 'Cu'); if (cu) { const other = bmNb(g, cu.atom.id).find(n => n.atom.id !== a.id && n.atom.element === 'C'); if (other && a.id < other.atom.id) P.cuprate.push({ c: a, c2: other.atom, cu: cu.atom }); }
        if (dO && ns.filter(n => n.atom.element === 'O' && n.bond.order === 2).length === 2 && ns.length === 2) P.co2.push({ c: a });
        if (ns.some(n => n.atom.element === 'N' && n.bond.order === 3)) P.nitrile.push({ c: a, n: ns.find(n => n.atom.element === 'N' && n.bond.order === 3).atom });
        /* an epoxide carbon: O in a three-ring */
        ns.forEach(n => { if (n.atom.element === 'O' && n.bond.order === 1) { const oth = bmNb(g, n.atom.id).find(m => m.atom.id !== a.id && m.atom.element === 'C'); if (oth && bmNb(g, a.id).some(m => m.atom.id === oth.atom.id)) P.epoxide.push({ o: n.atom, c: a, c2: oth.atom, sub: deg(g, a) }); } });
        /* a terminal alkyne C–H (pKa ≈ 25): deprotonated by strong bases */
        if (ns.some(n => n.bond.order === 3 && n.atom.element === 'C') && hOn(g, a) === 1) { P.acetylene.push({ c: a }); P.acidicH.push({ at: a, pka: 25, what: 'terminal alkyne C–H' }); }
        const X = ns.find(n => ['Br', 'Cl', 'I', 'F'].includes(n.atom.element) && n.bond.order === 1);
        if (X && X.atom.element === 'F' && bmSp3(g, a)) { /* an alkyl fluoride does not do Sₙ2 */ }
        else if (X && bmSp3(g, a)) P.halide.push({ c: a, x: X.atom, deg: deg(g, a), hind: hindrance(g, a.id), benzylic: ns.some(n => n.atom.element === 'C' && bmArom(g, n.atom.id)) });
        /* A SULFONATE ESTER is a leaving group on sp³ carbon exactly as a halide is, and this is
           the whole reason mesylates, tosylates and triflates are made: an alcohol cannot be
           displaced (HO⁻ is far too basic to leave), so the O–H is turned into O–SO₂R and the
           same carbon becomes an ordinary Sₙ2 electrophile. Reactivity: OTf ≫ OMs ≈ OTs ≈ Br.
           An ARYL sulfonate is not included here — sp² carbon does not do Sₙ2 — and the aryl
           triflate already has its own handle for palladium a few lines below. */
        if (bmSp3(g, a)) {
          const OS = ns.find(n => n.atom.element === 'O' && n.bond.order === 1 && bmNb(g, n.atom.id).some(m => m.atom.element === 'S' &&
            bmNb(g, m.atom.id).filter(k => k.atom.element === 'O' && k.bond.order === 2).length === 2));
          if (OS) {
            const S = bmNb(g, OS.atom.id).find(m => m.atom.element === 'S').atom;
            const on = bmNb(g, S.id).find(k => k.atom.element === 'C');
            const cf3 = on && bmNb(g, on.atom.id).filter(f => f.atom.element === 'F').length === 3;
            const me = on && !bmArom(g, on.atom.id) && hOn(g, on.atom) === 3;
            const word = cf3 ? 'triflate (OTf)' : me ? 'mesylate (OMs)' : on && bmArom(g, on.atom.id) ? 'tosylate (OTs)' : 'sulfonate ester';
            P.halide.push({ c: a, x: OS.atom, s: S, sulf: true, lgWord: word, tf: !!cf3,
              deg: deg(g, a), hind: hindrance(g, a.id), benzylic: ns.some(n => n.atom.element === 'C' && bmArom(g, n.atom.id)) });
          }
        }
        if (X && !bmSp3(g, a) && (bmArom(g, a.id) || ns.some(n => n.atom.element === 'C' && n.bond.order === 2))) {
          /* Sₙ2 never happens on an sp² carbon, but two other things do: a palladium
             cross-coupling (any aryl or vinyl halide will do), and — when the ring is short of
             electrons in the right place — plain addition–elimination, Sₙ Ar. */
          P.halideSp2.push({ c: a, x: X.atom, aryl: bmArom(g, a.id), act: bmArom(g, a.id) ? bmSnArActivation(g, a.id) : null, otf: false });
        }
        /* an aryl triflate is an aryl halide as far as palladium is concerned */
        const OTf = ns.find(n => n.atom.element === 'O' && n.bond.order === 1 && bmNb(g, n.atom.id).some(m => m.atom.element === 'S' &&
          bmNb(g, m.atom.id).filter(k => k.atom.element === 'O' && k.bond.order === 2).length === 2 &&
          bmNb(g, m.atom.id).some(k => k.atom.element === 'C' && bmNb(g, k.atom.id).filter(f => f.atom.element === 'F').length === 3)));
        if (OTf && !bmSp3(g, a) && bmArom(g, a.id)) P.halideSp2.push({ c: a, x: OTf.atom, aryl: true, act: bmSnArActivation(g, a.id), otf: true });
        if (bmArom(g, a.id) && hOn(g, a) >= 1) P.arene.push({ c: a });
        /* An acidic C–H on an electron-poor heteroaromatic ring — the handle for direct
           C–H arylation. There is no organometallic partner in that reaction: palladium
           deprotonates the ring itself (concerted metalation–deprotonation, which is what the
           carboxylate base is for) and couples it to an aryl halide. Which C–H depends on how
           well the ring can carry the negative charge: flanked by two ring nitrogens is best,
           then next to one, and an electron-withdrawing group on the ring helps further. */
        if (bmArom(g, a.id) && hOn(g, a) === 1) {
          let sys = null;
          try { sys = C().ringSystemsOf(g, C().ringAtomsOf(g)).find(r => r.indexOf(a.id) >= 0); } catch (e) { sys = null; }
          if (sys) {
            const set = new Set(sys);
            const ringHet = sys.filter(x => atom(g, x).element !== 'C');
            if (ringHet.length) {
              const nbIn = ns.filter(n => set.has(n.atom.id));
              const flank = nbIn.filter(n => atom(g, n.atom.id).element !== 'C').length;
              const five = (() => { try { return C().ringsOfSystem(g, sys).some(r => r.length === 5 && r.indexOf(a.id) >= 0); } catch (e) { return false; } })();
              /* an EWG anywhere on this ring acidifies every C–H on it */
              const ewg = sys.some(x => bmNb(g, x).some(n => !set.has(n.atom.id) &&
                ((n.atom.element === 'C' && (bmDblO(g, n.atom.id) || bmNb(g, n.atom.id).some(m => m.atom.element === 'N' && m.bond.order === 3))) ||
                 (n.atom.element === 'N' && n.atom.charge === 1))));
              P.azoleCH.push({ c: a, flank, five, ewg, nHet: ringHet.length, sys });
            }
          }
        }
        /* alkenes: a dienophile (any non-aromatic C=C) and, chained, a diene */
        ns.forEach(n => { if (n.atom.element === 'C' && n.bond.order === 2 && !bmArom(g, a.id) && !bmArom(g, n.atom.id) && a.id < n.atom.id) {
          const ewg = [a, n.atom].some(x => bmNb(g, x.id).some(m => (m.atom.element === 'C' && bmDblO(g, m.atom.id)) || (m.atom.element === 'N' && m.bond.order === 3) ||
            (m.atom.element === 'N' && m.atom.charge === 1) ||
            (m.atom.element === 'S' && bmNb(g, m.atom.id).filter(q => q.atom.element === 'O' && q.bond.order === 2).length >= 2)));
          P.dienophile.push({ a, b: n.atom, ewg });
          /* a=b–c=d */
          bmNb(g, n.atom.id).forEach(m => { if (m.atom.element === 'C' && m.bond.order === 1 && m.atom.id !== a.id && !bmArom(g, m.atom.id)) bmNb(g, m.atom.id).forEach(q => { if (q.atom.element === 'C' && q.bond.order === 2 && q.atom.id !== n.atom.id && !bmArom(g, q.atom.id)) P.diene.push({ c1: a, c2: n.atom, c3: m.atom, c4: q.atom }); }); });
          bmNb(g, a.id).forEach(m => { if (m.atom.element === 'C' && m.bond.order === 1 && m.atom.id !== n.atom.id && !bmArom(g, m.atom.id)) bmNb(g, m.atom.id).forEach(q => { if (q.atom.element === 'C' && q.bond.order === 2 && q.atom.id !== a.id && !bmArom(g, q.atom.id)) P.diene.push({ c1: n.atom, c2: a, c3: m.atom, c4: q.atom }); }); });
        } });
      }
      if (a.element === 'O' && ns.length === 1 && hOn(g, a) >= 1 && !a.charge) {
        const c = ns[0].atom;
        if (c.element === 'C' && !bmDblO(g, c.id)) { if (bmArom(g, c.id)) { P.phenol.push({ o: a, c }); P.acidicH.push({ at: a, pka: 10, what: 'phenol O–H' }); } else if (bmSp3(g, c)) { P.alcohol.push({ o: a, c, deg: deg(g, c) }); P.acidicH.push({ at: a, pka: 16, what: 'alcohol O–H' }); } }
      }
      /* strong bases that are not carbon: amide ion, hydride */
      if (a.element === 'N' && a.charge === -1 && !ns.length) P.strongBase.push({ c: a, kind: 'sodium amide (NaNH₂)', conj: 'ammonia', maxPka: 33 });
      if (a.element === 'N' && a.charge === -1 && ns.length === 2 && ns.every(n => n.atom.element === 'C')) P.strongBase.push({ c: a, kind: 'lithium diisopropylamide (LDA)', conj: 'the secondary amine', maxPka: 33 });
      if (a.element === 'H' && a.charge === -1) P.strongBase.push({ c: a, kind: 'sodium hydride (NaH)', conj: 'hydrogen gas', maxPka: 30 });
      if (a.element === 'O' && a.charge === -1 && ns.length === 1 && ns[0].atom.element === 'C' && !bmDblO(g, ns[0].atom.id)) P.alkoxide.push({ o: a, c: ns[0].atom });
      if (a.element === 'S' && ns.length === 1 && hOn(g, a) >= 1 && ns[0].atom.element === 'C') P.thiol.push({ s: a, c: ns[0].atom });
      if (a.element === 'N' && !a.charge && hOn(g, a) >= 1 && ns.every(n => n.bond.order === 1) && !ns.some(n => n.atom.element === 'C' && bmDblO(g, n.atom.id)) && !ns.some(n => n.atom.element === 'S')) P.amine.push({ n: a, h: hOn(g, a), aryl: ns.some(n => bmArom(g, n.atom.id)), deg: ns.length });
      if (a.element === 'N' && !a.charge && hOn(g, a) >= 1 && ns.every(n => n.bond.order === 1) && !ns.some(n => n.atom.element === 'S')) { const amide = ns.some(n => n.atom.element === 'C' && bmDblO(g, n.atom.id)); P.acidicH.push({ at: a, pka: amide ? 17 : ns.some(n => bmArom(g, n.atom.id)) ? 30 : 36, what: amide ? 'amide N–H' : ns.some(n => bmArom(g, n.atom.id)) ? 'aniline N–H' : 'amine N–H' }); }
      if (a.element === 'N' && a.charge === -1 && ns.length === 1) P.azide.push({ n: a });
      if (a.element === 'S' && ns.length === 4 && ns.filter(n => n.atom.element === 'O' && n.bond.order === 2).length === 2) { const cl = ns.find(n => n.atom.element === 'Cl'); if (cl) P.sulfonyl.push({ s: a, cl: cl.atom }); }
      /* An ACTIVE METHYLENE: an sp³ C–H flanked by electron-withdrawing groups. One EWG gives a
         pKa near 20 (a ketone α-H); TWO give a pKa near 11–13 (malonate, cyanoacetate,
         malononitrile, nitroalkane), and that is the dividing line that matters — a doubly
         activated carbon is deprotonated by a mere amine or carbonate, so these condensations
         run without any strong base, and the alkene they give is conjugated to BOTH groups,
         which is why the water leaves on its own instead of needing to be driven off. */
      if (a.element === 'C' && bmSp3(g, a) && hOn(g, a) >= 1) {
        let w = 0, ringNH = false;
        ns.forEach(n => {
          const t = n.atom;
          if (t.element === 'C' && n.bond.order === 1) {
            if (bmDblO(g, t.id)) w += 1;
            else if (bmNb(g, t.id).some(m => m.atom.element === 'N' && m.bond.order === 3)) w += 1;
          } else if (t.element === 'N' && t.charge === 1 && bmNb(g, t.id).filter(m => m.atom.element === 'O').length === 2) w += 1;
          else if (t.element === 'S' && bmNb(g, t.id).filter(m => m.atom.element === 'O' && m.bond.order === 2).length === 2) w += 0.9;
          else if (t.element === 'N' && hOn(g, t) >= 1) ringNH = true;
        });
        if (w >= 0.9) P.activeCH.push({ c: a, w, ringNH, doubly: w >= 1.8 || (w >= 0.9 && ringNH) });
      }
      /* A WITTIG partner. Both forms count: the isolated ylide (Ph₃P=CH–R) and the phosphonium
         SALT it is made from, since the base that makes the ylide is part of the recipe. */
      if (a.element === 'C') {
        const pp = ns.find(n => n.atom.element === 'P');
        if (pp && bmNb(g, pp.atom.id).filter(n => n.atom.element === 'C' && n.atom.id !== a.id).length >= 3) {
          const isYlide = pp.bond.order === 2 || a.charge === -1;
          if (isYlide || hOn(g, a) >= 1) {
            /* an ester or nitrile on the same carbon STABILISES the ylide, and a stabilised
               ylide gives the E alkene where an unstabilised one gives Z */
            const stab = ns.some(n => n.atom.element === 'C' && (bmDblO(g, n.atom.id) || bmNb(g, n.atom.id).some(m => m.atom.element === 'N' && m.bond.order === 3)));
            P.ylide.push({ c: a, p: pp.atom, salt: !isYlide, stab });
          }
        }
        /* A PHOSPHONATE: (RO)₂P(=O)–CH₂–EWG, the Horner–Wadsworth–Emmons reagent. */
        if (pp && hOn(g, a) >= 1) {
          const pn = bmNb(g, pp.atom.id);
          const dbl = pn.some(n => n.atom.element === 'O' && n.bond.order === 2);
          const ors = pn.filter(n => n.atom.element === 'O' && n.bond.order === 1).length;
          if (dbl && ors >= 2) P.phosphonate.push({ c: a, p: pp.atom,
            stab: ns.some(n => n.atom.element === 'C' && (bmDblO(g, n.atom.id) || bmNb(g, n.atom.id).some(m => m.atom.element === 'N' && m.bond.order === 3))) });
        }
      }
      /* An N-HALO IMIDE (NBS, NCS, NIS) or sulfuryl chloride. These are the reagents that
         halogenate an electron-rich ring — a pyrrole, an indole, an anisole — where molecular
         bromine with a Lewis acid would be far too violent and would give polyhalogenation.
         They work because the N–X bond is weak and polarised, so the ring meets a small,
         well-behaved amount of X⁺ at any moment. */
      if (a.element === 'N' && !a.charge) {
        const x = ns.find(n => ['Br', 'Cl', 'I'].includes(n.atom.element) && n.bond.order === 1);
        const imide = ns.filter(n => n.atom.element === 'C' && bmDblO(g, n.atom.id)).length >= 1;
        if (x && imide) P.nxs.push({ n: a, x: x.atom, el: x.atom.element, word: 'N-' + (x.atom.element === 'Br' ? 'bromo' : x.atom.element === 'Cl' ? 'chloro' : 'iodo') + 'succinimide (N' + x.atom.element + 'S)' });
      }
      if (a.element === 'S' && ns.filter(n => n.atom.element === 'O' && n.bond.order === 2).length === 2 &&
          ns.filter(n => n.atom.element === 'Cl').length === 2)
        P.nxs.push({ n: a, x: ns.find(n => n.atom.element === 'Cl').atom, el: 'Cl', word: 'sulfuryl chloride (SO₂Cl₂)' });
      /* LAWESSON'S REAGENT or P₄S₁₀: a phosphorus with sulfur doubly bonded to it, which swaps
         a C=O for a C=S. */
      /* An ENOLATE, drawn as the anion: C=C–O⁻. It is one species with the charge spread over
         both ends, and the CARBON is where it reacts — which is the single most important fact
         in the chapter. Also the ENOL (C=C–OH), which does the same chemistry more slowly. */
      if (a.element === 'O' && (a.charge === -1 || hOn(g, a) === 1) && ns.length === 1 && ns[0].atom.element === 'C' && ns[0].bond.order === 1) {
        const cO = ns[0].atom;
        const dbl = bmNb(g, cO.id).find(n => n.atom.element === 'C' && n.bond.order === 2 && !bmArom(g, n.atom.id));
        if (dbl) P.enolate.push({ o: a, cO, c: dbl.atom, anion: a.charge === -1,
          sub: deg(g, dbl.atom), hind: hindrance(g, dbl.atom.id) });
      }
      /* An ENAMINE: N–C=C, with the nitrogen's lone pair pushing electron density onto the
         far carbon. That far carbon is nucleophilic, which is what makes the Stork reaction
         work — and an enamine is a neutral nucleophile, so no strong base is ever needed. */
      if (a.element === 'N' && !a.charge && hOn(g, a) === 0) {
        const vin = ns.find(n => n.atom.element === 'C' && n.bond.order === 1 && bmNb(g, n.atom.id).some(m => m.atom.element === 'C' && m.bond.order === 2 && !bmArom(g, m.atom.id)));
        if (vin) {
          const far = bmNb(g, vin.atom.id).find(m => m.atom.element === 'C' && m.bond.order === 2 && !bmArom(g, m.atom.id));
          if (far && !bmDblO(g, vin.atom.id)) P.enamine.push({ n: a, cN: vin.atom, c: far.atom });
        }
      }
      /* an ESTER carbonyl as an ACYL ELECTROPHILE for a Claisen condensation */
      if (a.element === 'C' && bmDblO(g, a.id)) {
        const or = ns.find(n => n.atom.element === 'O' && n.bond.order === 1 && hOn(g, n.atom) === 0 && bmNb(g, n.atom.id).some(m => m.atom.element === 'C' && m.atom.id !== a.id));
        if (or) P.esterC.push({ c: a, o: ns.find(n => n.atom.element === 'O' && n.bond.order === 2).atom, or: or.atom });
      }
      /* A NITROSATING AGENT: nitrite, nitrous acid, or an alkyl nitrite. All three become NO⁺
         in acid, and that is the electrophile the amine attacks. */
      if (a.element === 'N' && ns.length <= 2) {
        const dO = ns.find(n => n.atom.element === 'O' && n.bond.order === 2);
        const sO = ns.find(n => n.atom.element === 'O' && n.bond.order === 1);
        if (dO && (ns.length === 1 || (sO && !ns.some(n => n.atom.element === 'C'))) && a.charge !== 1)
          P.nitrosating.push({ n: a, o: dO.atom, kind: 'nitrite' });
        else if (dO && sO && bmNb(g, sO.atom.id).some(n => n.atom.element === 'C'))
          P.nitrosating.push({ n: a, o: dO.atom, or: sO.atom, kind: 'alkyl nitrite' });
      }
      /* an AMIDE carbonyl — the substrate thionation cares about most, and not part of P.acyl
         (which only lists acyl donors with a leaving group) */
      if (a.element === 'C' && bmDblO(g, a.id)) {
        const nn = ns.find(n => n.atom.element === 'N' && n.bond.order === 1);
        const ox = ns.find(n => n.atom.element === 'O' && n.bond.order === 1);
        if (nn && !ox) P.amideC.push({ c: a, o: ns.find(n => n.atom.element === 'O' && n.bond.order === 2).atom, n: nn.atom, kind: 'amide' });
      }
      if (a.element === 'P' && ns.some(n => n.atom.element === 'S' && n.bond.order === 2))
        P.thionate.push({ p: a, s: ns.find(n => n.atom.element === 'S' && n.bond.order === 2).atom });
      /* the three lone pairs a peracid competes for, in the order it takes them */
      if (sulfideAt(g, a)) P.sulfide.push({ s: a, aryl: ns.some(n => bmArom(g, n.atom.id)) });
      if (pyridineN(g, a)) P.oxidN.push({ n: a, pyridine: true });
      else if (tertAmine(g, a)) P.oxidN.push({ n: a, pyridine: false });
      /* A PERACID (R–CO–O–OH: mCPBA, peracetic acid) or hydrogen peroxide. The reactive part is
         the same in both — a weak O–O bond whose outer oxygen is electrophilic. */
      if (a.element === 'O' && ns.length === 2 && hOn(g, a) === 0) {
        const oo = ns.find(n => n.atom.element === 'O' && n.bond.order === 1 && hOn(g, n.atom) === 1 && bmNb(g, n.atom.id).length === 1);
        const cc = ns.find(n => n.atom.element === 'C' && bmDblO(g, n.atom.id));
        if (oo && cc) P.peracid.push({ o1: a, o2: oo.atom, c: cc.atom, kind: 'peracid' });
      }
      if (a.element === 'O' && hOn(g, a) === 1 && ns.length === 1 && ns[0].atom.element === 'O' &&
          hOn(g, ns[0].atom) === 1 && bmNb(g, ns[0].atom.id).length === 1 && a.id < ns[0].atom.id) {
        P.peracid.push({ o1: a, o2: ns[0].atom, c: null, kind: 'hydrogen peroxide' });
      }
    });
    /* A BARE METAL — magnesium, lithium or zinc written on its own ([Mg], Mg, [Mg+2] all
       count: a lone atom is the metal, whatever charge the SMILES happened to carry). This is
       the reagent that MAKES an organometallic: Mg inserts into a C–X bond to give the
       Grignard reagent, two lithiums give the organolithium and LiX. */
    if (g.atoms.length === 1 && ['Mg', 'Li', 'Zn'].includes(g.atoms[0].element)) P.metal0.push({ m: g.atoms[0], el: g.atoms[0].element });
    try { P.biNuc = biNucleophiles(g); } catch (e) { P.biNuc = []; }
    try { P.biElec = biElectrophiles(g); } catch (e) { P.biElec = []; }
    /* An oxidant is an oxidant. mCPBA's own aryl chloride, its carbonyl and its acidic O–H are
       spectators: the O–O bond reacts orders of magnitude faster than any of them, so once a
       peracid is recognised the rest of that molecule stops offering itself as a partner. */
    if (P.peracid.length) {
      ['acyl', 'halideSp2', 'halide', 'arene', 'azoleCH', 'dienophile', 'diene', 'amine', 'alcohol', 'phenol', 'carboxylate', 'acidicH', 'sulfide', 'oxidN'].forEach(k => { if (P[k]) P[k] = []; });
    }
    return P;
  }
  /* put B beside A in one graph; B's ids follow A's */
  function bmMerge(A, B) {
    const g = clone(A); const off = Math.max(g.nextId || 0, ...g.atoms.map(a => a.id + 1));
    const map = new Map();
    const ax = Math.max(...g.atoms.map(a => a.x)) + 140, bx0 = Math.min(...B.atoms.map(a => a.x));
    const ay = g.atoms.reduce((t, a) => t + a.y, 0) / g.atoms.length, by = B.atoms.reduce((t, a) => t + a.y, 0) / B.atoms.length;
    B.atoms.forEach(a => { map.set(a.id, off + a.id); g.atoms.push(Object.assign({}, a, { id: off + a.id, x: a.x - bx0 + ax, y: a.y - by + ay })); });
    B.bonds.forEach(b => g.bonds.push(Object.assign({}, b, { a: map.get(b.a), b: map.get(b.b) })));
    g.nextId = off + Math.max(...B.atoms.map(a => a.id)) + 1;
    return { g, map };
  }
  const bmBond = (g, a, b, order) => { g.bonds.push({ a, b, order: order || 1 }); };
  /* cut the bond a–b and hand back the fragment that went with b (for a triflate leaving) */
  function splitOff2(g, a, b) {
    g.bonds = g.bonds.filter(x => !((x.a === a && x.b === b) || (x.b === a && x.a === b)));
    return splitOff(g, a);
  }
  /* the three stages of a palladium cycle, drawn: Ar–Pd(II)–X, then the two organic groups on
     the same palladium, then the product. Pd is shown explicitly because that is the point. */
  function bmPdStages(prod, cAr, cPartner, m, idX, idY) {
    try {
      const st = [];
      const a = clone(prod);
      /* undo the coupling on a copy so palladium can be drawn between the two carbons */
      a.bonds = a.bonds.filter(b => !((b.a === cAr && b.b === cPartner) || (b.b === cAr && b.a === cPartner)));
      const pd = C().addAtomNear(a, 'Pd', a.atoms.find(x => x.id === cAr));
      bmBond(a, pd.id, cPartner, 1);
      bmTidy(a);
      st.push({ label: 'both groups on the same palladium(II) — cis, and about to couple', graph: a, changed: new Set([cAr, cPartner, pd.id]) });
      return st;
    } catch (e) { return null; }
  }
  const bmName = g => { try { return nameOf(g); } catch (e) { return C().formula(g); } };
  const bmTidy = g => { try { C().layoutGraph(g); } catch (e) {} return g; };

  /* every way the two molecules could react, best first */
  function bmMatches(A, B) {
    const PA = bmProfile(A), PB = bmProfile(B);
    const out = [];
    const both = (fnAB) => { fnAB(A, PA, B, PB, false); fnAB(B, PB, A, PA, true); };
    /* 1. acyl donor + O/N nucleophile → ester / amide */
    both((X, PX, Y, PY, swapped) => {
      PX.acyl.forEach(ac => {
        const nucs = [].concat(PY.amine.map(n => ({ kind: 'amide', at: n.n, deg: n.deg, aryl: n.aryl })), PY.alcohol.map(n => ({ kind: 'ester', at: n.o, deg: n.deg })), PY.phenol.map(n => ({ kind: 'ester', at: n.o, phenol: true })));
        nucs.forEach(nu => {
          if (ac.kind === 'ester' && nu.kind === 'ester') return;                 // transesterification: not taught
          if (ac.kind === 'acid' && nu.phenol) return;                            // a phenol needs the anhydride / acyl chloride
          if (nu.kind === 'amide' && nu.deg >= 3) return;                         // a tertiary amine has no N–H to lose
          const score = (ac.kind === 'acyl chloride' ? 3 : ac.kind === 'anhydride' ? 2.5 : ac.kind === 'acid' ? 1.5 : 1) + (nu.kind === 'amide' ? 1.1 : 0) - (nu.deg >= 3 ? 1 : 0) - (nu.phenol ? 0.4 : 0);
          out.push({ join: 'acyl', score, X, Y, swapped, ac, nu });
        });
      });
      PX.sulfonyl.forEach(su => {
        PY.amine.forEach(n => { if (n.deg >= 3) return; out.push({ join: 'sulfonamide', score: 3, X, Y, swapped, su, nu: { at: n.n, kind: 'amide' } }); });
        /* the alcohol: this is the reaction that MAKES a leaving group out of an alcohol */
        PY.alcohol.forEach(n => out.push({ join: 'sulfonamide', score: 2.6 - 0.3 * (n.deg >= 2 ? 1 : 0), X, Y, swapped, su, nu: { at: n.o, kind: 'sulfonate' } }));
        PY.phenol.forEach(n => out.push({ join: 'sulfonamide', score: 2.3, X, Y, swapped, su, nu: { at: n.o, kind: 'sulfonate', phenol: true } }));
      });
    });
    /* 1b. an aryl or vinyl halide + a coupling partner: the palladium-catalysed couplings, and
           Sₙ Ar for the cases where the ring is activated and no metal is needed at all */
    both((X, PX, Y, PY, swapped) => {
      PX.halideSp2.forEach(h => {
        const xr = bmOxAdd(h);
        PY.boron.forEach(b => out.push({ join: 'suzuki', score: 4.2 + xr, X, Y, swapped, h, p: b }));
        /* only the sp² group on tin transfers: the methyls of a trimethylaryl stannane are
           spectators, which is why PhSnMe₃ gives the biaryl and not toluene */
        PY.stannane.forEach(t => out.push({ join: 'stille', score: 3.4 + xr + (t.sp2 ? 0.8 : -2.2), X, Y, swapped, h, p: t }));
        PY.acetylene.forEach(k => out.push({ join: 'sonogashira', score: 3.9 + xr, X, Y, swapped, h, p: k }));
        PY.grignard.forEach(gr => { if (gr.acet) return; out.push({ join: gr.mg.element === 'Zn' ? 'negishi' : 'kumada', score: 3.2 + xr, X, Y, swapped, h, p: gr }); });
        PY.amine.forEach(n => {
          const nu = { at: n.n, deg: n.deg, aryl: n.aryl, kind: 'amine' };
          out.push({ join: 'buchwald', score: 3.6 + xr - (n.deg >= 3 ? 1.2 : 0), X, Y, swapped, h, p: nu });
          if (h.act && !h.otf) out.push({ join: 'snar', score: 2.4 + Math.min(2.4, h.act.w) - (n.deg >= 3 ? 1.2 : 0) - betterHalideElsewhere(X, h), X, Y, swapped, h, p: nu });
        });
        if (h.act && !h.otf) {
          PY.alkoxide.forEach(n => out.push({ join: 'snar', score: 2.7 + Math.min(2.4, h.act.w), X, Y, swapped, h, p: { at: n.o, kind: 'alkoxide' } }));
          PY.thiol.forEach(n => out.push({ join: 'snar', score: 2.7 + Math.min(2.4, h.act.w), X, Y, swapped, h, p: { at: n.s, kind: 'thiol' } }));
          PY.phenol.forEach(n => out.push({ join: 'snar', score: 2.1 + Math.min(2.4, h.act.w), X, Y, swapped, h, p: { at: n.o, kind: 'phenol' } }));
        }
        /* direct C–H arylation: no boron, no tin, no zinc — the heteroarene's own C–H */
        PY.azoleCH.forEach(k => {
          /* charge stabilisation at that carbon decides which C–H is taken. Two flanking ring
             nitrogens (C2 of an imidazole, oxazole or thiazole) is the classic answer; but when
             the ring also carries an electron-withdrawing group the other position competes and
             usually wins, which is what the recorded arylations of 4-substituted azoles show. */
          const acid = 1.4 * k.flank + (k.five ? 0.4 : 0) + 0.3 * (k.nHet - 1) + (k.ewg && k.flank < 2 ? 1.6 : 0);
          /* If the same heteroarene still has a free N–H, that nitrogen is what gets arylated:
             N-arylation is easier, and it is what the Buchwald conditions in the flask are for.
             C–H arylation is the reaction you run on an N-protected ring, with a carboxylate
             base. Ranking C–H above N–H cost 25 of the recorded AstraZeneca couplings. */
          const hasNH = PY.amine.length > 0;
          out.push({ join: 'charyl', score: 2.6 + xr + acid * 0.5 - (hasNH ? 1.6 : 0), X, Y, swapped, h, p: k, acid: acid });
        });
        /* Heck: the aryl group ends up on the LESS substituted end of the alkene and the new
           double bond comes out E — insertion, then syn β-hydride elimination */
        PY.dienophile.forEach(d => {
          if (hOn(Y, d.a) + hOn(Y, d.b) === 0) return;
          out.push({ join: 'heck', score: 3.1 + xr, X, Y, swapped, h, p: d });
        });
      });
    });
    /* 1c. a peracid (or H₂O₂) + something with a lone pair. WHICH lone pair goes first is the
           whole question, and nucleophilicity decides it, not what looks reactive:
           sulfide ≫ amine / pyridine N > electron-rich alkene. */
    both((X, PX, Y, PY, swapped) => {
      PX.peracid.forEach(px => {
        PY.sulfide.forEach(su => out.push({ join: 'oxidS', score: 4.4, X, Y, swapped, px, su }));
        PY.oxidN.forEach(nn => out.push({ join: 'oxidN', score: 3.4 - (nn.pyridine ? 0 : 0.3), X, Y, swapped, px, nn }));
        PY.dienophile.forEach(d => { if (d.ewg) return; out.push({ join: 'oxidC', score: 2.4, X, Y, swapped, px, d }); });
      });
    });
    /* 2. alkyl halide + nucleophile → Sₙ2 */
    both((X, PX, Y, PY, swapped) => {
      PX.halide.forEach(h => {
        if (h.deg >= 3) return;
        const nucs = [].concat(PY.amine.map(n => ({ kind: 'amine', at: n.n, name: 'N-alkylation' })), PY.thiol.map(n => ({ kind: 'thiol', at: n.s, name: 'S-alkylation (thioether)' })),
          PY.alkoxide.map(n => ({ kind: 'alkoxide', at: n.o, name: 'Williamson ether synthesis' })), PY.alcohol.map(n => ({ kind: 'alcohol', at: n.o, name: 'Williamson ether synthesis (after NaH deprotonates the alcohol)' })),
          PY.phenol.map(n => ({ kind: 'phenol', at: n.o, name: 'Williamson ether synthesis (K₂CO₃ deprotonates the phenol)' })), PY.cyanide.map(n => ({ kind: 'cyanide', at: n.c, name: 'cyanide substitution' })),
          PY.azide.map(n => ({ kind: 'azide', at: n.n, name: 'azide substitution' })), PY.carboxylate.map(n => ({ kind: 'carboxylate', at: n.c, name: 'carboxylate alkylation (ester)' })));
        nucs.forEach(nu => out.push({ join: 'sn2', score: (nu.kind === 'thiol' || nu.kind === 'alkoxide' || nu.kind === 'cyanide' ? 2.5 : nu.kind === 'amine' ? 2.2 : 1.6) - (h.deg === 2 ? 0.6 : 0) + (h.benzylic ? 0.4 : 0) + (h.tf ? 0.5 : 0), X, Y, swapped, h, nu }));
      });
    });
    /* 3. carbanion reagents: Grignard / organolithium / acetylide + carbonyl, CO₂, ester, nitrile, epoxide; Gilman 1,4 and coupling */
    both((X, PX, Y, PY, swapped) => {
      PX.grignard.forEach(gr => {
        PY.carbonyl.forEach(co => out.push({ join: 'grignard', score: 3 + (PY.enone.some(e => e.cC === co.c) ? -0.3 : 0), X, Y, swapped, gr, co }));
        PY.co2.forEach(co => out.push({ join: 'grignard-co2', score: 3, X, Y, swapped, gr, co }));
        PY.ester.forEach(es => out.push({ join: 'grignard-ester', score: 2.4, X, Y, swapped, gr, es }));
        PY.nitrile.forEach(ni => out.push({ join: 'grignard-nitrile', score: 2.2, X, Y, swapped, gr, ni }));
        PY.epoxide.forEach(ep => out.push({ join: 'grignard-epoxide', score: 2.3, X, Y, swapped, gr, ep }));
        if (gr.acet) PY.halide.filter(h => h.deg <= 1 && !h.sulf).forEach(h => out.push({ join: 'acetylide-sn2', score: 2.6, X, Y, swapped, gr, h }));
      });
      /* An epoxide with a real nucleophile and NO acid: a plain Sₙ2 on a strained ring, at the
         LESS hindered carbon — the opposite end from the acid-catalysed opening. */
      PX.epoxide.forEach(ep => {
        const nucs = [].concat(
          PY.alkoxide.map(n => ({ at: n.o, kind: 'alkoxide', word: 'the alkoxide' })),
          PY.thiol.map(n => ({ at: n.s, kind: 'thiolate', word: 'the thiolate' })),
          PY.cyanide.map(n => ({ at: n.c, kind: 'cyanide', word: 'cyanide' })),
          PY.azide.map(n => ({ at: n.n, kind: 'azide', word: 'azide' })),
          PY.amine.filter(n => n.deg <= 2).map(n => ({ at: n.n, kind: 'amine', word: 'the amine' })));
        nucs.forEach(nu => out.push({ join: 'epoxNu', score: 3.2 + (nu.kind === 'amine' ? -0.4 : 0), X, Y, swapped, ep, nu }));
      });
      PX.cuprate.forEach(cu => {
        PY.enone.forEach(en => out.push({ join: 'gilman-14', score: 3.2, X, Y, swapped, cu, en }));
        PY.halide.filter(h => h.deg <= 1).forEach(h => out.push({ join: 'gilman-couple', score: 2.6, X, Y, swapped, cu, h }));
        PY.halideSp2.forEach(h => out.push({ join: 'gilman-couple', score: 2.5, X, Y, swapped, cu, h, sp2: true }));
        PY.acyl.filter(ac => ac.kind === 'acyl chloride').forEach(ac => out.push({ join: 'gilman-acyl', score: 2.8, X, Y, swapped, cu, ac }));
        PY.carbonyl.forEach(co => { if (!PY.enone.some(e => e.cC === co.c)) out.push({ join: 'gilman-none', score: 0.5, X, Y, swapped, cu, co }); });
      });
      /* a strong base takes the most acidic proton */
      PX.strongBase.forEach(bs => { const ah = PY.acidicH.filter(a => a.pka <= (bs.maxPka || 40)).sort((p, q) => p.pka - q.pka)[0]; if (ah) out.push({ join: 'deprotonate', score: 1.8 + (ah.pka <= 5 ? 0.5 : 0) + (ah.pka >= 30 ? -0.6 : 0), X, Y, swapped, bs, ah }); });
    });
    /* 4. aldol: an enolisable carbonyl + a carbonyl electrophile */
    both((X, PX, Y, PY, swapped) => {
      PX.enol.forEach(en => PY.carbonyl.forEach(co => {
        if (co.c === en.c && X === Y) return;
        /* the electrophile should be the aldehyde, ideally one with no α-H (an aromatic aldehyde) */
        const coHasAlpha = PY.enol.some(e => e.c === co.c);
        const score = 1.2 + (co.kind === 'aldehyde' ? 0.8 : 0) + (!coHasAlpha ? 0.8 : 0) + (en.kind === 'ketone' ? 0.3 : 0);
        out.push({ join: 'aldol', score, X, Y, swapped, en, co });
      }));
    });
    /* 4b. Knoevenagel condensation, Michael addition, Wittig and Horner–Wadsworth–Emmons.
           All four make a bond from a stabilised carbanion, and all four are decided by the
           same question: how acidic is that C–H, and what does the carbanion meet. */
    both((X, PX, Y, PY, swapped) => {
      PX.activeCH.forEach(ac => {
        if (!ac.doubly) return;
        /* + an aldehyde or ketone: the water goes on its own, because the alkene that forms is
           conjugated to both withdrawing groups */
        PY.carbonyl.forEach(co => {
          const coHasAlpha = PY.activeCH.some(e => bmNb(Y, co.c.id).some(n => n.atom.id === e.c.id));
          out.push({ join: 'knoev', score: 3.0 + (co.kind === 'aldehyde' ? 0.9 : 0) + (coHasAlpha ? 0 : 0.4), X, Y, swapped, ac, co });
        });
        /* + an enone: 1,4-addition. This is THE Michael reaction — the classic version uses a
           malonate or a cyanoacetate exactly because a mere amine or carbonate deprotonates it. */
        PY.enone.forEach(en => out.push({ join: 'michael', score: 3.3, X, Y, swapped, ac, en }));
        /* + an alkyl halide: the malonic ester and acetoacetic ester syntheses start here */
        PY.halide.forEach(h => { if (h.deg >= 3) return; out.push({ join: 'acAlkyl', score: 3.1 + (h.benzylic ? 0.3 : 0) - (h.deg === 2 ? 0.7 : 0), X, Y, swapped, ac, h }); });
      });
      /* Wittig: a phosphorus ylide (or the salt plus base) + a carbonyl → an alkene, and the
         C=O oxygen ends up on phosphorus. */
      PX.ylide.forEach(y => PY.carbonyl.forEach(co => out.push({ join: 'wittig', score: 4.0 - (y.salt ? 0.3 : 0), X, Y, swapped, y, co })));
      PX.phosphonate.forEach(ph => PY.carbonyl.forEach(co => out.push({ join: 'hwe', score: 3.6 + (ph.stab ? 0.4 : -1.2), X, Y, swapped, ph, co })));
    });
    /* 4c. HETEROCYCLE-FORMING CONDENSATION: two nucleophiles a fixed distance apart meet one
           or two electrophilic carbons a fixed distance apart, and the ring closes. The only
           arithmetic is the ring size, and only five and six are allowed. */
    both((X, PX, Y, PY, swapped) => {
      PX.biNuc.forEach(bn => PY.biElec.forEach(be => {
        if (X === Y) return;
        /* ring size = the nucleophiles + whatever is between them + the carbons + whatever is
           between those. A span of −1 means "the same atom does both jobs". */
        const nAt = bn.span < 0 ? 1 : 2 + bn.span;
        const cAt = be.span < 0 ? 1 : 2 + be.span;
        const size = nAt + cAt;
        if (size !== 5 && size !== 6) return;
        if (bn.span < 0 && be.span < 0) return;               // one atom onto one carbon is not a ring
        /* a thiourea's sulfur takes the C–X carbon and its nitrogen the C=O: that is the
           regiochemistry of the Hantzsch thiazole synthesis, and it is not arbitrary */
        let score = 3.4 + (size === 5 ? 0.5 : 0)
          + (be.span >= 0 ? 1.0 : 0)                          // a real bis-electrophile: both carbons are used
          + (bn.arom ? 0.5 : 0)                               // a benzo-fused product is a strong driving force
          + (be.kind === 'orthoester' ? 0.6 : 0)
          + (bn.kind === 'thiourea' && be.kind === 'alpha-halo ketone' ? 0.8 : 0)
          + ((bn.kind === 'hydrazine' || bn.kind === 'hydroxylamine') && be.span === 1 ? 0.8 : 0)
          + ((bn.kind === 'amidine' || bn.kind === 'urea') && be.span === 1 ? 0.7 : 0)
          + (bn.kind === 'amine' && be.span === 2 ? 0.6 : 0);
        out.push({ join: 'hetcyc', score, X, Y, swapped, bn, be, size });
      }));
    });
    /* 4d. halogenation of an electron-rich ring with NBS / NCS / NIS or SO₂Cl₂, and thionation
           of a C=O with Lawesson's reagent. Both are one reagent meeting one site, but the
           reagent arrives as a molecule rather than as a chosen reagent id. */
    both((X, PX, Y, PY, swapped) => {
      PX.nxs.forEach(nx => { if (PY.arene.length) out.push({ join: 'arhal', score: 3.2, X, Y, swapped, nx }); });
      PX.thionate.forEach(th => {
        PY.acyl.forEach(ac => { if (ac.kind !== 'acid' && ac.kind !== 'acyl chloride' && ac.kind !== 'anhydride') out.push({ join: 'thion', score: 3.4, X, Y, swapped, th, ac }); });
        PY.carbonyl.forEach(co => out.push({ join: 'thion', score: 2.6, X, Y, swapped, th, co }));
        /* an amide is the classic substrate: amide >> ketone > ester, and an acid not at all */
        PY.amideC.forEach(am => out.push({ join: 'thion', score: 4.0, X, Y, swapped, th, ac: am }));
      });
    });
    /* 4h. STEP GROWTH: two molecules, each with TWO reactive ends, joining alternately. It is
           the same acyl substitution or addition you already know, repeated — and the reason
           the molecular weight climbs slowly is that any two pieces can join, so long chains
           only appear near the very end of the reaction. */
    both((X, PX, Y, PY, swapped) => {
      const twoAcyl = PX.acyl.filter(a => a.kind === 'acid' || a.kind === 'acyl chloride');
      const nucN = PY.amine.filter(n => n.deg <= 2);
      const nucO = PY.alcohol.concat(PY.phenol);
      const isocyanate = [];
      Y.atoms.forEach(a => { if (a.element !== 'N') return;
        const c = bmNb(Y, a.id).find(n => n.atom.element === 'C' && n.bond.order === 2 && bmNb(Y, n.atom.id).some(m => m.atom.element === 'O' && m.bond.order === 2));
        if (c) isocyanate.push({ n: a, c: c.atom }); });
      const bi = (arr) => arr.length >= 2;
      /* `elec` is the monomer with the two electrophilic ends and `nuc` the one with the two
         nucleophilic ends — named outright, so the builder never has to guess which is which */
      if (bi(twoAcyl) && bi(nucN)) out.push({ join: 'stepPoly', score: 4.6, X, Y, swapped, kind: 'amide', elec: X, nuc: Y });
      else if (bi(twoAcyl) && bi(nucO)) out.push({ join: 'stepPoly', score: 4.6, X, Y, swapped, kind: 'ester', elec: X, nuc: Y });
      if (bi(isocyanate) && bi(PX.alcohol.concat(PX.phenol))) out.push({ join: 'stepPoly', score: 4.6, X, Y, swapped, kind: 'urethane', elec: Y, nuc: X });
    });
    /* 4f. oxidative insertion of a metal into a C–X bond: the FORMATION of a Grignard reagent
           (or an organolithium). Alkyl, vinyl and aryl halides all do it; fluorides do not. */
    both((X, PX, Y, PY, swapped) => {
      PX.metal0.forEach(mt => {
        const acidic = PY.acidicH.filter(a => a.pka <= 25);       // O–H, N–H, COOH, terminal alkyne C–H
        const all = PY.halide.filter(h => !h.sulf).map(h => ({ c: h.c, x: h.x, sp2: false, deg: h.deg, benzylic: h.benzylic }))
          .concat(PY.halideSp2.filter(h => !h.otf).map(h => ({ c: h.c, x: h.x, sp2: true, aryl: h.aryl })));
        all.forEach(h => {
          if (h.x.element === 'F') return;
          const rank = { I: 0.6, Br: 0.4, Cl: 0 }[h.x.element] || 0;
          out.push({ join: 'metalInsert', score: 3.6 + rank - (acidic.length ? 3.0 : 0), X, Y, swapped, mt, h, acidic: acidic[0] || null });
        });
      });
    });
    /* 4e. diazotisation as a molecule + molecule reaction: an aryl amine meets a nitrite */
    both((X, PX, Y, PY, swapped) => {
      PX.nitrosating.forEach(ni => {
        let sites = [];
        try { sites = arylAmineSites(Y); } catch (e) { sites = []; }
        sites.forEach(am => {
          if (am.kind === 'tertiary') return;
          out.push({ join: 'diazotise', score: (am.kind === 'primary' && am.aryl ? 4.2 : am.kind === 'secondary' ? 3.0 : 1.0), X, Y, swapped, ni, am });
        });
      });
    });
    /* 4g. THE ENOLATE AND THE ENAMINE, meeting each of the four electrophiles the course uses:
           an alkyl halide (α-alkylation), a carbonyl (aldol), an enone (Michael), an ester
           (Claisen). One nucleophile, four reactions — which is the tidiest way to hold the
           whole of chapter 22 in your head. */
    both((X, PX, Y, PY, swapped) => {
      PX.enolate.forEach(en => {
        const boost = en.anion ? 0 : -1.4;                    // the neutral enol is slower
        PY.halide.forEach(h => { if (h.deg >= 3) return; out.push({ join: 'alkylate', score: 3.4 + boost + (h.benzylic ? 0.3 : 0) - (h.deg === 2 ? 0.8 : 0), X, Y, swapped, en, h }); });
        PY.carbonyl.forEach(co => out.push({ join: 'aldolE', score: 3.2 + boost + (co.kind === 'aldehyde' ? 0.6 : 0), X, Y, swapped, en, co }));
        PY.enone.forEach(e2 => out.push({ join: 'michaelE', score: 3.6 + boost, X, Y, swapped, en, e2 }));
        PY.esterC.forEach(es => out.push({ join: 'claisen', score: 3.0 + boost, X, Y, swapped, en, es }));
      });
      PX.enamine.forEach(em => {
        PY.halide.forEach(h => { if (h.deg >= 3) return; out.push({ join: 'stork', score: 3.0 + (h.benzylic ? 0.4 : 0), X, Y, swapped, em, h }); });
        PY.enone.forEach(e2 => out.push({ join: 'storkMichael', score: 3.2, X, Y, swapped, em, e2 }));
      });
      /* a SECONDARY amine and a ketone give the ENAMINE, not the imine — there is no second
         N–H to lose, so the molecule dehydrates towards carbon instead */
      PX.carbonyl.forEach(co => PY.amine.filter(n => n.deg === 2).forEach(n => {
        const hasAlpha = PX.enol.some(e => e.c === co.c) || alphaSites(X).some(a => a.c.id === co.c.id);
        if (!hasAlpha) return;
        out.push({ join: 'enamine', score: 2.2, X, Y, swapped, co, nu: { at: n.n } });
      }));
    });
    /* 5. carbonyl + 1° amine → imine */
    both((X, PX, Y, PY, swapped) => { PX.carbonyl.forEach(co => PY.amine.filter(n => n.h >= 2).forEach(n => out.push({ join: 'imine', score: 1.0, X, Y, swapped, co, nu: { at: n.n } }))); });
    /* 6. Diels–Alder */
    both((X, PX, Y, PY, swapped) => { PX.diene.forEach(d => PY.dienophile.forEach(dp => { if (X === Y) return;
      /* the regiochemistry term only ever PENALISES the wrong pairing: the best a
         Diels-Alder can score is unchanged, so it cannot elbow another family aside */
      out.push({ join: 'da', score: 2 + (dp.ewg ? 1 : 0) + Math.min(0, daRegio(X, d, Y, dp)), X, Y, swapped, d, dp }); })); });
    /* 7. Friedel–Crafts */
    both((X, PX, Y, PY, swapped) => {
      if (!PX.arene.length) return;
      /* Friedel–Crafts needs a ring that is not deactivated (no C=O / NO₂ / SO₃H / CN on it — the
         AlCl₃ complex will not attack such a ring) and no free amine (it ties up the Lewis acid);
         an O–H / N–H on the substrate is acylated first, so the ring reaction is only a minor path */
      let ringOk = true;
      try { const sites = areneSites(X, 'sulfonation', { T: 298, pHVal: 0 }).filter(s => !s.inert); ringOk = sites.length > 0 && !sites.every(s => (s.subs || []).some(su => su.d && su.d.dir === 'm' && su.d.v <= -1.5)); } catch (e) {}
      const donor = PY.acyl.some(ac => ac.kind === 'acyl chloride' || ac.kind === 'anhydride') || PY.halide.some(h => h.deg >= 2);
      if (!ringOk) { if (donor) out.push({ join: 'fc-none', score: -9, X, Y, swapped, why: 'the ring carries a strong electron-withdrawing group (C=O, NO₂, SO₃H, CN …), so it is DEACTIVATED — Friedel–Crafts acylation and alkylation fail on deactivated rings (the AlCl₃ complex cannot attack it)' }); return; }
      if (PX.amine.some(n => n.aryl)) { if (donor) out.push({ join: 'fc-none', score: -9, X, Y, swapped, why: 'the ring carries a free NH₂ — its lone pair ties up the Lewis acid, so Friedel–Crafts fails on anilines (the amine is acylated instead)' }); return; }
      const hetero = PX.alcohol.length + PX.phenol.length + PX.amine.length > 0 ? 1.6 : 0;
      PY.acyl.filter(ac => ac.kind === 'acyl chloride' || ac.kind === 'anhydride').forEach(ac => out.push({ join: 'fc-acyl', score: 2.2 - hetero, X, Y, swapped, ac }));
      PY.halide.filter(h => h.deg >= 2).forEach(h => out.push({ join: 'fc-alkyl', score: 1.4 - hetero, X, Y, swapped, h }));
    });
    out.sort((p, q) => q.score - p.score);
    return out;
  }

  /* build the product of one match; returns { graph, changed, added, byprod, label, why, mechKind, stages, eqn } */
  function bmBuild(m, env, sub) {
    const A = m.X, B = m.Y;                     // X carries the electrophile (acyl / halide / …), Y the partner
    /* the drawn molecule keeps its atom ids (the player follows atoms by id); the reagent's are offset */
    const subIsX = A === sub;
    const { g, map } = bmMerge(sub, subIsX ? B : A);
    const idX = subIsX ? (id => id) : (id => map.get(id)), idY = subIsX ? (id => map.get(id)) : (id => id);
    const by = [], changed = new Set();
    const heat = env && env.T >= 330;
    let label = '', why = '', mechKind = 'acylsub', eqn = null, stages = null, nameHint = null, stereoHow = null;
    if (m.join === 'acyl' || m.join === 'sulfonamide') {
      const ac = m.join === 'acyl' ? m.ac : null;
      const eC = idX(ac ? ac.c.id : m.su.s.id), nuAt = idY(m.nu.at.id);
      bmBond(g, eC, nuAt, 1);
      /* the tetrahedral intermediate as a stage — the leaving group still on, C=O opened to O⁻, N⁺ if an amine */
      try { const t = clone(g); const dO = t.atoms.find(a => a.id === (ac ? idX(ac.o.id) : -1)); if (dO) { const bd = t.bonds.find(b => (b.a === eC && b.b === dO.id) || (b.b === eC && b.a === dO.id)); if (bd) { bd.order = 1; dO.charge = -1; } if (m.nu.kind === 'amide' || !ac) { const na = t.atoms.find(a => a.id === nuAt); if (na) na.charge = 1; } bmTidy(t); stages = [{ label: 'tetrahedral intermediate (the nucleophile has added to the C=O)', graph: t }]; } } catch (e) {}
      /* the leaving group */
      if (ac) {
        const lg = idX(ac.lg.id);
        if (ac.kind === 'acid') { C().removeAtom(g, lg); by.push({ smiles: 'O', name: 'water' }); }
        else if (ac.kind === 'acyl chloride') { C().removeAtom(g, lg); by.push({ smiles: 'Cl', name: 'HCl' }); }
        else { g.bonds = g.bonds.filter(b => !((b.a === eC && b.b === lg) || (b.b === eC && b.a === lg))); const other = splitOff(g, eC); by.push({ graph: bmTidy(other), name: ac.kind === 'anhydride' ? 'the carboxylic acid (other half of the anhydride)' : 'the alcohol (from the ester)' }); }
      } else { C().removeAtom(g, idX(m.su.cl.id)); by.push({ smiles: 'Cl', name: 'HCl' }); }
      changed.add(eC); changed.add(nuAt);
      const isAmide = m.nu.kind === 'amide';
      label = !ac ? (isAmide ? 'sulfonamide' : m.nu.phenol ? 'aryl sulfonate ester' : 'sulfonate ester (a leaving group made from the alcohol)') : isAmide ? 'amide' : 'ester';
      const donor = !ac ? 'sulfonyl chloride' : ac.kind;
      why = (isAmide ? 'The amine nitrogen' : 'The ' + (m.nu.phenol ? 'phenol' : 'alcohol') + ' oxygen') + ' attacks the ' + (ac ? 'carbonyl carbon' : 'sulfur') + ' of the ' + donor + ' (nucleophilic ' + (ac ? 'acyl' : 'sulfonyl') + ' substitution): addition gives the tetrahedral intermediate, then the leaving group (' + by[0].name + ') is expelled and the ' + label + ' is left.' +
        (!ac && !isAmide ? ' The point of doing this is what comes NEXT: the alcohol oxygen could never leave on its own (hydroxide is far too basic), but as a sulfonate ester the same carbon is an ordinary Sₙ2 electrophile, and the sulfonate anion leaves about as readily as bromide. Pyridine or triethylamine is there to take up the HCl; keep it cold, because a benzylic or allylic sulfonate can ionise. Nothing about the carbon changes in this step, so a stereocentre keeps its configuration here — the inversion comes in the substitution afterwards.' : '') +
        (ac && ac.kind === 'acid' ? (isAmide ? ' With the FREE acid nothing useful happens on mixing: the amine simply takes the proton and you have the ammonium carboxylate salt, which needs 150–200 °C to condense. ' +
          'So in practice the acid is activated first, and this — amide formation from an acid and an amine — is the single most-run reaction in medicinal chemistry. The reagent choices: a carbodiimide (EDC·HCl with HOBt, or DCC — HOBt is there to intercept the O-acylisourea before it rearranges to the dead N-acylurea, and to suppress racemisation of an α-stereocentre); a uronium salt (HATU, HBTU) when the acid or the amine is hindered, which is the fastest and most reliable; CDI, which gives an isolable acyl imidazole; T3P on scale, because its by-products wash out; or the classical two-step route through the acyl chloride with SOCl₂ or oxalyl chloride and then the amine with a base. All of them do the same thing — turn the –OH into a better leaving group — and all go through the same tetrahedral intermediate drawn here.' :
          ' Fischer esterification: an acid catalyst and heat, with the water removed or the alcohol in excess, since every step is reversible.') : '');
      eqn = { R: by.filter(b => b.smiles).map(b => [b.smiles === 'O' ? 'H2O' : 'HCl', b.name + ' — the BYPRODUCT']), cat: ac && ac.kind === 'acid' ? (isAmide ? ['Δ', 'heat — drives off the water'] : ['H₂SO₄ cat., Δ', 'the acid catalyst and heat of a Fischer esterification']) : ['pyridine', 'a base to take up the HCl'] };
    } else if (m.join === 'suzuki' || m.join === 'stille' || m.join === 'negishi' || m.join === 'kumada' || m.join === 'sonogashira') {
      /* All four transmetalation couplings do the same three things to the drawing: the C–X
         bond goes, the partner's carbon–metal (or carbon–boron) bond goes, and the two carbons
         are joined. What differs is what carries the carbon, how forgiving the conditions are,
         and what is left over. */
      const cX = idX(m.h.c.id), xId = idX(m.h.x.id);
      const pC = idY(m.join === 'sonogashira' ? m.p.c.id : m.p.c.id);
      const carrier = m.join === 'suzuki' ? idY(m.p.b.id) : m.join === 'stille' ? idY(m.p.sn.id) : (m.join === 'negishi' || m.join === 'kumada') ? idY(m.p.mg.id) : null;
      if (m.h.otf) { const tf = splitOff2(g, cX, xId); by.push({ graph: bmTidy(tf), name: 'the triflate anion' }); }
      else { C().removeAtom(g, xId); }
      if (carrier !== null) g.bonds = g.bonds.filter(b => !((b.a === pC && b.b === carrier) || (b.b === pC && b.a === carrier)));
      bmBond(g, cX, pC, 1);
      if (carrier !== null) {
        const off = splitOff(g, cX);
        by.push({ graph: bmTidy(off), name: m.join === 'suzuki' ? 'the boron residue (boric acid, or the diol from a boronate ester — washed out with the base)' : m.join === 'stille' ? 'the tin residue R₃SnX — toxic, and the reason this coupling is avoided when a Suzuki will do' : 'the ' + (m.join === 'negishi' ? 'zinc' : 'magnesium') + ' halide salt' });
      } else {
        by.push({ smiles: m.h.otf ? null : m.h.x.element, name: 'H' + (m.h.otf ? 'OTf' : m.h.x.element) + ' — taken up by the amine base' });
      }
      changed.add(cX); changed.add(pC);
      mechKind = 'crosscouple';
      const named = { suzuki: 'Suzuki–Miyaura', stille: 'Stille', negishi: 'Negishi', kumada: 'Kumada', sonogashira: 'Sonogashira' }[m.join];
      const partWord = m.join === 'suzuki' ? m.p.kind : m.join === 'stille' ? 'organostannane' : m.join === 'negishi' ? 'organozinc' : m.join === 'kumada' ? 'Grignard reagent' : 'terminal alkyne';
      label = (m.join === 'sonogashira' ? 'internal alkyne' : m.p.aryl && m.h.aryl ? 'biaryl' : 'coupled product') + ' (' + named + ' coupling)';
      why = 'A palladium cross-coupling, and the catalytic cycle is the same three steps every time. ' +
        '(1) OXIDATIVE ADDITION: Pd(0) inserts into the C–' + (m.h.otf ? 'O' : m.h.x.element) + ' bond of the aryl ' + bmXWord(m.h) + ', giving a square-planar Ar–Pd(II)–' + (m.h.otf ? 'OTf' : m.h.x.element) + '. This is the slow step, and it is why the halide matters so much: iodide and triflate go easily, bromide readily, chloride only with a strong electron-rich ligand, fluoride not at all. ' +
        '(2) TRANSMETALATION: the ' + partWord + ' hands its carbon to palladium' + (m.join === 'suzuki' ? ' — but only after the base converts the neutral boronic acid into the boronate, which is why a Suzuki without base does nothing' : m.join === 'sonogashira' ? ' — copper(I) first makes the copper acetylide, which is the species that actually transfers' : '') + ', displacing the halide. ' +
        '(3) REDUCTIVE ELIMINATION: the two organic groups, now cis on the same palladium, couple and fall off as the product, regenerating Pd(0) for the next turn. Nothing is consumed but the two partners and the base.';
      eqn = { cat: m.join === 'suzuki' ? ['Pd(PPh₃)₄ (cat.), K₂CO₃', 'the palladium catalyst and the base that activates the boron']
        : m.join === 'sonogashira' ? ['Pd(PPh₃)₂Cl₂ / CuI (cat.), Et₃N', 'palladium, the copper co-catalyst and an amine base']
        : m.join === 'stille' ? ['Pd(PPh₃)₄ (cat.)', 'the palladium catalyst — no base needed, which is this coupling\u2019s one advantage']
        : ['Pd or Ni catalyst', 'the metal catalyst; keep it dry — the organometallic is basic'] };
      stages = bmPdStages(g, cX, pC, m, idX, idY);
    } else if (m.join === 'charyl') {
      const cX = idX(m.h.c.id), cH = idY(m.p.c.id);
      if (m.h.otf) { const tf = splitOff2(g, cX, idX(m.h.x.id)); by.push({ graph: bmTidy(tf), name: 'the triflate anion' }); }
      else C().removeAtom(g, idX(m.h.x.id));
      bmBond(g, cX, cH, 1);
      changed.add(cX); changed.add(cH);
      by.push({ smiles: null, name: 'the carboxylic acid (pivalic or acetic) regenerated, and the base\u2019s salt of H' + (m.h.otf ? 'OTf' : m.h.x.element) });
      mechKind = 'crosscouple';
      label = 'arylated heteroarene (direct C–H arylation)';
      why = 'Direct C–H arylation — a cross-coupling with only one partner prepared. There is no boronic acid, no stannane and no Grignard: the heteroarene is coupled through its own C–H bond, which is why this reaction matters so much industrially (it removes a whole step, and the waste that goes with it). ' +
        'Pd(0) adds into the aryl–' + bmXWord(m.h) + ' bond as usual. Then, instead of transmetalation, comes CONCERTED METALATION–DEPROTONATION: the carboxylate on palladium (pivalate or acetate — this is what the carboxylate base is actually for, and why a carbonate alone is much slower) reaches across and takes the heteroarene\u2019s C–H at the same moment the C–Pd bond forms, in one transition state. Reductive elimination then gives the biaryl. ' +
        'Which C–H is taken follows how well that carbon can carry negative charge: a position between two ring nitrogens goes first, then one next to a single nitrogen' + (m.p.ewg ? ', and an electron-withdrawing group on the ring shifts the balance to the position it acidifies — which is the outcome recorded for 4-substituted azoles' : '') + '.';
      eqn = { cat: ['Pd(OAc)₂ or PdCl(allyl) (cat.), PCy₃, KOPiv, DMA, 100–140 °C', 'palladium with a phosphine, and a carboxylate base — the carboxylate is a reagent in the mechanism, not just a proton sponge'] };
      stages = bmPdStages(g, cX, cH, m, idX, idY);
    } else if (m.join === 'buchwald') {
      const cX = idX(m.h.c.id), nAt = idY(m.p.at.id);
      if (m.h.otf) { const tf = splitOff2(g, cX, idX(m.h.x.id)); by.push({ graph: bmTidy(tf), name: 'the triflate anion' }); }
      else C().removeAtom(g, idX(m.h.x.id));
      bmBond(g, cX, nAt, 1);
      changed.add(cX); changed.add(nAt);
      by.push({ smiles: null, name: 'the base\u2019s salt (Cs₂CO₃ or NaOtBu takes up the H' + (m.h.otf ? 'OTf' : m.h.x.element) + ')' });
      mechKind = 'crosscouple';
      label = 'aryl amine (Buchwald–Hartwig amination)';
      why = 'Buchwald–Hartwig amination — the modern way to put a nitrogen on an aromatic ring, and the reaction this app sees more often than any other in the patent literature. Pd(0) adds into the C–' + (m.h.otf ? 'O' : m.h.x.element) + ' bond; the base deprotonates the amine on palladium to give an Ar–Pd(II)–NR₂ amido complex; reductive elimination forms the C–N bond and returns Pd(0). ' +
        'The ligand is the whole story: a bulky electron-rich biaryl phosphine (XPhos, RuPhos, BrettPhos) speeds up both the oxidative addition and the reductive elimination, and picking the wrong one is the usual reason a Buchwald fails. ' + (m.p.aryl ? 'An aniline is a poorer nucleophile than an alkylamine but couples well with the right ligand.' : m.p.deg >= 3 ? 'A tertiary amine has no N–H to lose, so it cannot couple at all.' : 'A primary or secondary alkylamine is the easy case.');
      eqn = { cat: ['Pd₂(dba)₃ / XPhos (cat.), Cs₂CO₃', 'palladium, a bulky biaryl phosphine ligand and a mild base'] };
      stages = bmPdStages(g, cX, nAt, m, idX, idY);
    } else if (m.join === 'heck') {
      /* the aryl group lands on the less substituted alkene carbon; the C=C reappears E */
      const cX = idX(m.h.c.id);
      const aH = hOn(m.Y, m.p.a), bH = hOn(m.Y, m.p.b);
      const endA = aH >= bH ? m.p.a : m.p.b, endB = aH >= bH ? m.p.b : m.p.a;
      const cEnd = idY(endA.id), cIn = idY(endB.id);
      if (m.h.otf) { const tf = splitOff2(g, cX, idX(m.h.x.id)); by.push({ graph: bmTidy(tf), name: 'the triflate anion' }); }
      else C().removeAtom(g, idX(m.h.x.id));
      bmBond(g, cX, cEnd, 1);
      const bd = g.bonds.find(b => (b.a === cEnd && b.b === cIn) || (b.b === cEnd && b.a === cIn));
      if (bd) bd.order = 1;
      /* β-hydride elimination puts the new double bond between the two former alkene carbons
         again — one hydrogen has moved from the far carbon to palladium and then to the base */
      if (bd) bd.order = 2;
      changed.add(cX); changed.add(cEnd); changed.add(cIn);
      by.push({ smiles: null, name: 'Et₃N·H' + (m.h.otf ? 'OTf' : m.h.x.element) + ' (the amine base takes up the acid)' });
      mechKind = 'crosscouple';
      label = 'arylated alkene, E (Heck reaction)';
      why = 'The Heck reaction. Pd(0) adds into the aryl–' + (m.h.otf ? 'triflate' : bmXWord(m.h)) + ' bond, then the alkene coordinates and the Ar–Pd unit adds ACROSS it syn — aryl to the less hindered carbon, palladium to the other. Now syn β-hydride elimination removes a hydrogen from the carbon next to palladium and regenerates a C=C; the amine base takes the H–Pd–X apart to give back Pd(0). ' +
        'Because both the addition and the elimination are syn and the intermediate rotates to put the bulky groups apart, the new double bond comes out E almost exclusively. No organometallic partner is needed at all, which is why the Heck is the coupling of choice when a plain alkene is what you have.';
      eqn = { cat: ['Pd(OAc)₂ / PPh₃ (cat.), Et₃N', 'palladium, phosphine and an amine base'] };
      stages = bmPdStages(g, cX, cEnd, m, idX, idY);
    } else if (m.join === 'snar') {
      const cX = idX(m.h.c.id), nuAt = idY(m.p.at.id);
      const xEl = m.h.x.element;
      /* the Meisenheimer intermediate: the nucleophile added, the ring dearomatised, the
         negative charge parked on the activating group */
      let meis = null;
      try {
        meis = clone(g);
        bmBond(meis, cX, nuAt, 1);
        const na = meis.atoms.find(a => a.id === nuAt); if (na && na.charge === -1) na.charge = 0; else if (na) na.charge = 1;
        const ring = meis.bonds.find(b => (b.a === cX || b.b === cX) && b.order === 2);
        if (ring) { ring.order = 1; const far = meis.atoms.find(a => a.id === (ring.a === cX ? ring.b : ring.a)); if (far) far.charge = -1; }
        bmTidy(meis);
      } catch (e) { meis = null; }
      C().removeAtom(g, idX(m.h.x.id));
      bmBond(g, cX, nuAt, 1);
      const nuA = g.atoms.find(a => a.id === nuAt);
      if (nuA && nuA.charge === -1) nuA.charge = 0;
      changed.add(cX); changed.add(nuAt);
      by.push({ smiles: xEl, name: 'H' + xEl + (m.p.kind === 'amine' ? ' (a second equivalent of the amine takes it up)' : '') });
      mechKind = 'SNAr';
      label = m.p.kind === 'amine' ? 'aryl amine (Sₙ Ar)' : m.p.kind === 'thiol' ? 'aryl thioether (Sₙ Ar)' : 'aryl ether (Sₙ Ar)';
      why = 'Nucleophilic aromatic substitution — no metal needed, because this ring is activated. The nucleophile adds straight onto the carbon carrying the ' + bmXWord(m.h) + ', which breaks the aromaticity and gives the MEISENHEIMER anion; ' +
        'that anion survives only because the negative charge can be delocalised onto ' + (m.h.act ? m.h.act.words.join(' and ') : 'the activating group') + '. The ring then pushes the ' + xEl + '⁻ out and the aromaticity comes back. ' +
        'Two consequences worth knowing: the leaving-group order is the OPPOSITE of Sₙ2 (F > Cl > Br > I, because what matters is how much the C–X bond polarises the ring, not how well X⁻ leaves), and an activating group meta to the halide does nothing at all — the charge cannot reach it.';
      eqn = { R: [[xEl === 'F' ? 'HF' : 'H' + xEl, by[0].name + ' — the BYPRODUCT']], cat: ['K₂CO₃ or Et₃N, DMSO/DMF, Δ', 'a base to mop up the acid, in a polar aprotic solvent'] };
      if (meis) stages = [{ label: 'Meisenheimer anion — the ring has lost its aromaticity and the charge sits on the activating group', graph: meis }];
    } else if (m.join === 'oxidS' || m.join === 'oxidN' || m.join === 'oxidC') {
      /* exactly one oxygen atom is transferred, and what is left of the oxidant is the
         carboxylic acid (from a peracid) or water (from H₂O₂) */
      const oGone = idX(m.px.o2.id), oStay = idX(m.px.o1.id);
      const oat = g.atoms.find(a => a.id === oGone);
      if (oat) { oat.charge = 0; oat.hAdd = 0; oat.hExp = 0; }
      if (m.join === 'oxidC') {
        const a1 = idY(m.d.a.id), b1 = idY(m.d.b.id);
        const bd = g.bonds.find(b => (b.a === a1 && b.b === b1) || (b.b === a1 && b.a === b1));
        if (bd) bd.order = 1;
        bmBond(g, oGone, a1, 1); bmBond(g, oGone, b1, 1);
        changed.add(a1); changed.add(b1); changed.add(oGone);
        label = 'epoxide'; mechKind = 'epox';
        why = 'One concerted step (the butterfly transition state): the alkene’s π electrons attack the peracid’s electrophilic OUTER oxygen while its O–H proton swings across to the carbonyl oxygen of the same molecule. Both new C–O bonds are made on the SAME face, so whatever was cis on the alkene stays cis on the epoxide, and the more electron-rich the alkene the faster it goes.';
      } else if (m.join === 'oxidS') {
        const sAt = idY(m.su.s.id);
        bmBond(g, oGone, sAt, 2);
        changed.add(sAt); changed.add(oGone);
        label = 'sulfoxide'; mechKind = 'oxid';
        why = 'Sulfur’s lone pair is the best nucleophile in the flask — better than an alkene’s π bond and better than an amine — so the sulfide takes the peracid’s oxygen first, whatever else is present. ' +
          'And it STOPS at the sulfoxide with one equivalent: the sulfoxide sulfur is already electron-poor and much slower to oxidise again, so reaching the sulfone means deliberately using two equivalents, a stronger oxidant, or heat. ' +
          'The sulfur becomes a stereocentre here (two different carbons, an oxygen and a lone pair); a plain peracid gives it racemic.';
      } else {
        const nId = idY(m.nn.n.id);
        bmBond(g, oGone, nId, 1);
        const nAt = g.atoms.find(a => a.id === nId);
        if (nAt) nAt.charge = 1;
        if (oat) oat.charge = -1;
        changed.add(nId); changed.add(oGone);
        label = m.nn.pyridine ? 'pyridine N-oxide' : 'amine N-oxide'; mechKind = 'oxid';
        why = 'The nitrogen lone pair takes the outer oxygen of the peracid, giving a real N⁺–O⁻ dipole (drawn with the charges, not with a dative arrow). ' +
          (m.nn.pyridine ? 'For a pyridine this is usually a means rather than an end: the N-oxide nitrates far more easily than pyridine itself and does so at C4, and C2 becomes open to nucleophiles — which is how most 2-substituted pyridines are actually made. PCl₃, P(OEt)₃ or H₂/Pd takes the oxide off again afterwards.'
                         : 'A tertiary amine gives the N-oxide cleanly. Heat it with a β-hydrogen present and it does the Cope elimination instead, giving the alkene and a hydroxylamine.');
      }
      /* only now is the spent oxidant cut away, with the transferred oxygen already bonded */
      const spent = cutFragment(g, oGone, oStay);
      by.push(Object.assign({ graph: bmTidy(spent) }, { name: m.px.kind === 'peracid'
        ? 'the carboxylic acid (mCPBA → 3-chlorobenzoic acid; peracetic acid → acetic acid)' : 'water' }));
      try { C().layoutGraph(g); } catch (e) {}
      eqn = { cat: [m.px.kind === 'peracid' ? 'CH₂Cl₂ or AcOH, 0 °C → rt' : 'AcOH or MeOH, rt',
        'one equivalent of oxidant — two would take a sulfide on to the sulfone'] };
    } else if (m.join === 'sn2') {
      const cX = idX(m.h.c.id), nuAt = idY(m.nu.at.id);
      const xEl = m.h.x.element;
      const lgWord = m.h.sulf ? m.h.lgWord : xEl;
      if (m.h.sulf) {
        /* the whole –OSO₂R group leaves as the sulfonate anion: one atom cannot simply be deleted */
        bmBond(g, cX, nuAt, 1);
        const frag = cutFragment(g, cX, idX(m.h.x.id));
        by.push({ graph: bmTidy(frag), name: 'the ' + m.h.lgWord + ' anion — the leaving group' });
      } else {
        C().removeAtom(g, idX(m.h.x.id));
        bmBond(g, cX, nuAt, 1);
      }
      const nuA = g.atoms.find(a => a.id === nuAt);
      if (nuA && nuA.charge === -1) nuA.charge = 0;               // the alkoxide / cyanide / azide is now neutral
      changed.add(cX); changed.add(nuAt);
      const salt = m.nu.kind === 'alkoxide' || m.nu.kind === 'cyanide' || m.nu.kind === 'azide' || m.nu.kind === 'carboxylate';
      if (!m.h.sulf) by.push(salt ? { smiles: '[Na+].[' + xEl + '-]', name: 'Na' + xEl } : { smiles: xEl, name: 'H' + xEl + (m.nu.kind === 'amine' ? ' (taken up by a second amine as the salt)' : '') });
      label = m.nu.name; mechKind = 'SN2';
      why = 'Sₙ2: the ' + (m.nu.kind === 'amine' ? 'nitrogen' : m.nu.kind === 'thiol' ? 'sulfur' : m.nu.kind === 'cyanide' ? 'cyanide carbon' : 'oxygen') + ' lone pair attacks the ' + clsName(m.h.deg) + ' carbon carrying the ' + lgWord + ' from the back side, and the ' + lgWord + ' leaves. ' +
        (m.h.sulf ? 'That carbon started as an ALCOHOL. Hydroxide is far too basic to be displaced, so the O–H was converted to a sulfonate ester first — which is the only purpose of MsCl, TsCl and Tf₂O. The sulfonate anion that leaves is stabilised by three oxygens' + (m.h.tf ? ' and by the trifluoromethyl group, which is what makes a triflate the best of the three by orders of magnitude' : '') + ', so it goes as readily as bromide. ' : '') +
        (m.h.deg === 2 ? 'A secondary halide is slower and a strong base would compete by E2. ' : '') + (m.nu.kind === 'alcohol' ? 'The alcohol itself is too weak a nucleophile — NaH makes the alkoxide first.' : m.nu.kind === 'amine' ? 'Use excess amine: the product amine is itself a nucleophile and would react again.' : '');
      eqn = { R: m.h.sulf ? [] : [[by[0].smiles === xEl ? 'H' + xEl : 'Na' + xEl, by[0].name + ' — the BYPRODUCT']], cat: m.nu.kind === 'alcohol' ? ['NaH', 'the base that makes the alkoxide'] : m.nu.kind === 'phenol' ? ['K₂CO₃', 'the base that makes the phenoxide'] : undefined };
    } else if (m.join === 'grignard') {
      const cC = idY(m.co.c.id), oId = idY(m.co.o.id), rC = idX(m.gr.c.id), mgId = idX(m.gr.mg.id);
      /* R leaves the magnesium and lands on the carbonyl carbon; the C=O becomes C–O⁻ → C–OH at the workup */
      bmBond(g, cC, rC, 1);
      g.bonds = g.bonds.filter(b => !((b.a === rC && b.b === mgId) || (b.b === rC && b.a === mgId)));
      const mEl = m.gr.mg.element, isMg = mEl === 'Mg';
      const mgFrag = splitOff(g, cC); by.push({ graph: mgFrag, name: isMg ? 'Mg(OH)X (magnesium salt, from the workup)' : mEl + 'OH (from the workup)' });
      const bd = g.bonds.find(b => (b.a === cC && b.b === oId) || (b.b === cC && b.a === oId)); if (bd) bd.order = 1;
      changed.add(cC); changed.add(rC); changed.add(oId);
      const addWord = isMg ? 'Grignard addition' : m.gr.kind + ' addition';
      label = (m.co.kind === 'aldehyde' ? (m.gr.acet ? '2° propargylic alcohol (' : '2° alcohol (') : (m.gr.acet ? '3° propargylic alcohol (' : '3° alcohol (')) + addWord + ')'; mechKind = 'grignard';
      why = 'The ' + (isMg ? 'Grignard' : m.gr.kind) + ' carbon (a carbanion equivalent) attacks the electrophilic carbonyl carbon; the π electrons move onto oxygen as the magnesium alkoxide, and the dilute-acid workup gives the alcohol. Everything must be dry until the workup.';
      eqn = { L: [['H2O', 'the acid workup — the PROTON SOURCE for the alkoxide']], cat: ['dry Et₂O; then H₃O⁺', 'anhydrous ether, then aqueous acid'] };
    } else if (m.join === 'epoxNu') {
      /* which carbon? With no acid there is no protonation and no positive charge to
         stabilise, so it is sterics alone: the LESS substituted carbon. */
      const eps = [m.ep];
      const a = m.ep.c, b = m.ep.c2;
      const subA = deg(m.X, a), subB = deg(m.X, b);
      const hit = subA <= subB ? a : b, other = subA <= subB ? b : a;
      const cHit = idX(hit.id), cOther = idX(other.id), oId = idX(m.ep.o.id);
      const nuAt = idY(m.nu.at.id);
      const bd = g.bonds.find(x => (x.a === cHit && x.b === oId) || (x.b === cHit && x.a === oId));
      if (bd) g.bonds = g.bonds.filter(x => x !== bd);
      bmBond(g, cHit, nuAt, 1);
      const nA = g.atoms.find(x => x.id === nuAt); if (nA && nA.charge === -1) nA.charge = 0;
      const oA = g.atoms.find(x => x.id === oId); if (oA) { oA.charge = 0; oA.hAdd = 0; }
      try { markAnti(g, cHit, nuAt, cOther, oId); } catch (e) {}
      changed.add(cHit); changed.add(oId); changed.add(nuAt);
      label = 'ring opened at the LESS hindered carbon (β-' + (m.nu.kind === 'amine' ? 'amino' : m.nu.kind === 'cyanide' ? 'hydroxy nitrile' : m.nu.kind === 'thiolate' ? 'hydroxy thioether' : m.nu.kind === 'azide' ? 'azido' : 'alkoxy') + ' alcohol)';
      mechKind = 'SN2';
      why = 'An epoxide opened by a real nucleophile, with no acid anywhere — and the regiochemistry is the OPPOSITE of the acid-catalysed opening, which is the point of the comparison. ' +
        'There is no protonated oxygen here and so no positive charge building on carbon; nothing rewards attacking the more substituted end. It is a plain Sₙ2, and Sₙ2 goes where there is room: the LESS hindered carbon. The oxygen leaves as an alkoxide and picks up a proton at the workup. ' +
        'Put the two side by side: acid → the more substituted carbon, base or a strong nucleophile → the less substituted one. And in both cases the nucleophile arrives from the face opposite the oxygen, so the two groups end up anti. ' +
        'What makes an epoxide worth the trouble is that a three-membered ring holds about 27 kcal/mol of strain, so it reacts under conditions no other ether would notice — and it delivers two functional groups on adjacent carbons in one step.';
      stereoHow = 'Backside attack: inversion at the carbon attacked, and the nucleophile and the OH end up anti. On a ring that means trans.';
      eqn = { L: [['H2O', 'the workup — the proton source for the alkoxide']], cat: ['the sodium salt of the nucleophile, then H₃O⁺', 'no acid: acid would send it to the other carbon'] };
    } else if (m.join === 'grignard-co2' || m.join === 'grignard-ester' || m.join === 'grignard-nitrile' || m.join === 'grignard-epoxide' || m.join === 'acetylide-sn2') {
      const rC = idX(m.gr.c.id), mId = idX(m.gr.mg.id);
      const dropMetal = keep => { g.bonds = g.bonds.filter(b => !((b.a === rC && b.b === mId) || (b.b === rC && b.a === mId))); const frag = splitOff(g, keep); by.push({ graph: frag, name: m.gr.mg.element === 'Mg' ? 'Mg(OH)X (magnesium salt, from the workup)' : m.gr.mg.element + 'OH (from the workup)' }); };
      const kindWord = m.gr.kind;
      if (m.join === 'grignard-co2') {
        const cC = idY(m.co.c.id); bmBond(g, cC, rC, 1); dropMetal(cC);
        const o1 = bmNb(g, cC).find(n => n.atom.element === 'O'); if (o1) { const bd = g.bonds.find(b => (b.a === cC && b.b === o1.atom.id) || (b.b === cC && b.a === o1.atom.id)); if (bd) bd.order = 1; }
        changed.add(cC); changed.add(rC); label = 'carboxylic acid (carboxylation with CO₂)'; mechKind = 'grignard';
        why = 'The ' + kindWord + ' attacks the electrophilic carbon of CO₂ (dry ice) to give the magnesium carboxylate; acid workup liberates the carboxylic acid — one carbon longer than the halide it came from.';
        eqn = { L: [['H2O', 'the acid workup — proton source']], cat: ['dry ice; then H₃O⁺', 'solid CO₂, then aqueous acid'] };
      } else if (m.join === 'grignard-ester') {
        const cC = idY(m.es.c.id), oId = idY(m.es.o.id), orId = idY(m.es.or.id);
        /* two equivalents: the first gives the ketone (which cannot be stopped at), the second the tertiary alcohol */
        bmBond(g, cC, rC, 1); dropMetal(cC);
        g.bonds = g.bonds.filter(b => !((b.a === cC && b.b === orId) || (b.b === cC && b.a === orId))); const alk = splitOff(g, cC); by.push({ graph: bmTidy(alk), name: 'the alcohol from the ester (leaves as the alkoxide)' });
        const R2 = bmMerge(g, m.X); const rC2 = R2.map.get(m.gr.c.id), m2 = R2.map.get(m.gr.mg.id); Object.assign(g, R2.g);
        bmBond(g, cC, rC2, 1); g.bonds = g.bonds.filter(b => !((b.a === rC2 && b.b === m2) || (b.b === rC2 && b.a === m2))); splitOff(g, cC);
        const bd = g.bonds.find(b => (b.a === cC && b.b === oId) || (b.b === cC && b.a === oId)); if (bd) bd.order = 1;
        [cC, rC, rC2, oId].forEach(id => changed.add(id)); label = '3° alcohol (two equivalents add to the ester)'; mechKind = 'grignard';
        why = 'The first ' + kindWord + ' adds to the ester carbonyl and the alkoxide leaves, giving a ketone — which is MORE reactive than the ester, so a second equivalent adds at once. Two R groups end up on the carbon: a tertiary alcohol after the workup.';
        eqn = { L: [['H2O', 'the acid workup']], cat: ['2 equiv; dry Et₂O; then H₃O⁺', 'two equivalents of the organometallic'] };
      } else if (m.join === 'grignard-nitrile') {
        const cC = idY(m.ni.c.id), nId = idY(m.ni.n.id);
        bmBond(g, cC, rC, 1); dropMetal(cC);
        /* the imine anion is hydrolysed at the workup: C≡N → C=O */
        C().removeAtom(g, nId); const oNew = C().addAtomNear(g, 'O', g.atoms.find(a => a.id === cC)); const bd = g.bonds.find(b => (b.a === cC && b.b === oNew.id) || (b.b === cC && b.a === oNew.id)); if (bd) bd.order = 2;
        by.push({ smiles: 'N', name: 'ammonia (from hydrolysing the imine)' });
        changed.add(cC); changed.add(rC); changed.add(oNew.id); label = 'ketone (addition to the nitrile, then hydrolysis)'; mechKind = 'grignard';
        why = 'The ' + kindWord + ' adds once to the C≡N to give a metal imine salt; it cannot add twice because the anion is not electrophilic. Aqueous acid at the workup hydrolyses the imine to the ketone.';
        eqn = { L: [['H2O', 'the hydrolysis', 2]], cat: ['dry Et₂O; then H₃O⁺', ''] };
      } else if (m.join === 'grignard-epoxide') {
        /* attack at the less substituted carbon; the oxygen stays on the other one */
        const ep = m.ep; const eps = bmProfile(m.Y).epoxide.filter(e => e.o === ep.o).sort((p, q) => p.sub - q.sub); const pick = eps[0] || ep;
        const cA = idY(pick.c.id), cB = idY(pick.c2.id), oId = idY(pick.o.id);
        g.bonds = g.bonds.filter(b => !((b.a === cA && b.b === oId) || (b.b === cA && b.a === oId)));
        bmBond(g, cA, rC, 1); dropMetal(cA);
        [cA, cB, oId, rC].forEach(id => changed.add(id)); label = 'alcohol (epoxide opened, two carbons from the new C–C bond)'; mechKind = 'SN2';
        why = 'The ' + kindWord + ' attacks the LESS hindered epoxide carbon from the back (Sₙ2), the strained ring springs open, and the alkoxide on the other carbon is protonated at the workup — the chain grows by two carbons.';
        eqn = { L: [['H2O', 'the acid workup']], cat: ['dry Et₂O; then H₃O⁺', ''] };
      } else {  /* acetylide-sn2 */
        const cX = idY(m.h.c.id); C().removeAtom(g, idY(m.h.x.id)); bmBond(g, cX, rC, 1);
        g.bonds = g.bonds.filter(b => !((b.a === rC && b.b === mId) || (b.b === rC && b.a === mId))); const salt = splitOff(g, cX); by.push({ graph: salt, name: m.gr.mg.element + m.h.x.element + ' (the salt)' });
        changed.add(cX); changed.add(rC); label = 'internal alkyne (acetylide alkylation)'; mechKind = 'SN2';
        why = 'The acetylide carbon is a carbanion nucleophile: it attacks the primary carbon carrying the halide from the back (Sₙ2). Only primary halides work — secondary and tertiary ones eliminate with a base this strong.';
        eqn = { R: [[m.gr.mg.element + m.h.x.element, 'the salt']], cat: ['liquid NH₃ or THF', ''] };
      }
    } else if (m.join === 'gilman-14' || m.join === 'gilman-couple' || m.join === 'gilman-acyl' || m.join === 'gilman-none') {
      if (m.join === 'gilman-none') return null;
      const rC = idX(m.cu.c.id), cuId = idX(m.cu.cu.id);
      const dropCu = keep => { g.bonds = g.bonds.filter(b => !((b.a === rC && b.b === cuId) || (b.b === rC && b.a === cuId))); const frag = splitOff(g, keep); by.push({ graph: frag, name: 'RCu (the other R stays on copper) + LiX' }); };
      if (m.join === 'gilman-14') {
        const cB = idY(m.en.beta.id), cA = idY(m.en.alpha.id);
        bmBond(g, cB, rC, 1); const bd = g.bonds.find(b => (b.a === cA && b.b === cB) || (b.b === cA && b.a === cB)); if (bd) bd.order = 1; dropCu(cB);
        [cB, cA, rC].forEach(id => changed.add(id)); label = 'β-alkylated ketone (conjugate 1,4-addition)'; mechKind = 'gilman';
        why = 'A Gilman reagent (lithium dialkylcuprate) is a SOFT carbon nucleophile: it adds to the β-carbon of the enone (1,4- or Michael addition), not to the C=O the way a Grignard would. The enolate that forms is protonated at the workup, so the C=O is back and R sits on the β-carbon.';
        eqn = { L: [['H2O', 'the workup']], cat: ['R₂CuLi, Et₂O, −78 °C; then H₃O⁺', 'the cuprate, cold, then aqueous workup'] };
      } else if (m.join === 'gilman-couple') {
        const cX = idY(m.h.c.id); C().removeAtom(g, idY(m.h.x.id)); bmBond(g, cX, rC, 1); dropCu(cX);
        changed.add(cX); changed.add(rC); label = 'C–C coupling (Corey–House)'; mechKind = 'gilman';
        why = m.sp2 ? 'One R of the cuprate replaces the halide on the ' + (m.h.aryl ? 'aromatic' : 'vinyl') + ' carbon — the Corey–House coupling. An sp² C–X cannot do Sₙ2 at all, and a Grignard would not couple here, but the cuprate does (with retention of the alkene geometry for a vinyl halide).' : 'One R of the cuprate replaces the halide on the primary carbon — the Corey–House coupling. Cuprates do this cleanly where a Grignard would give a mess (Wurtz coupling and elimination).';
        eqn = { R: [['LiX', 'lithium halide'], ['RCu', 'the second R stays on copper']], cat: ['Et₂O, 0 °C', ''] };
      } else {
        const cC = idY(m.ac.c.id); bmBond(g, cC, rC, 1); C().removeAtom(g, idY(m.ac.lg.id)); dropCu(cC);
        changed.add(cC); changed.add(rC); label = 'ketone (cuprate + acyl chloride)'; mechKind = 'acylsub';
        why = 'A cuprate is nucleophilic enough for an acyl chloride but not for the ketone it makes, so the reaction stops cleanly at the ketone — unlike a Grignard, which would add a second time.';
        eqn = { R: [['LiCl', 'lithium chloride'], ['RCu', 'the second R stays on copper']], cat: ['Et₂O, −78 °C', ''] };
      }
    } else if (m.join === 'deprotonate') {
      const hAt = idY(m.ah.at.id), bC = idX(m.bs.c.id);
      const a = g.atoms.find(x => x.id === hAt); if (a) a.charge = -1;
      /* the base leaves as its conjugate acid — drawn as a by-product */
      const baseFrag = splitOff(g, hAt); const mAt = baseFrag.atoms.find(x => ['Li', 'Na', 'K', 'Mg'].includes(x.element));
      if (mAt) { const hal = baseFrag.atoms.find(x => ['Cl', 'Br', 'I'].includes(x.element)); const metal = C().parseSmiles(hal ? '[' + mAt.element + '+]' : '[' + mAt.element + '+]'); by.push({ graph: metal, name: (hal ? mAt.element + hal.element + '⁺' : mAt.element + '⁺') + ' (counter-ion)' }); }
      if (m.bs.kind === 'organolithium' || m.bs.kind === 'Grignard reagent' || m.bs.kind === 'metal acetylide') {
        /* R–M + H–A → R–H + M–A: the conjugate acid is the reagent's own carbon skeleton with the metal (and its halide) stripped off */
        const alk = clone(baseFrag); const gone = new Set(alk.atoms.filter(x => METALS.includes(x.element) || ['Cl', 'Br', 'I'].includes(x.element)).map(x => x.id));
        alk.atoms = alk.atoms.filter(x => !gone.has(x.id)); alk.bonds = alk.bonds.filter(bd => !gone.has(bd.a) && !gone.has(bd.b));
        if (alk.atoms.length) by.push({ graph: bmTidy(alk), name: (m.bs.kind === 'metal acetylide' ? 'the terminal alkyne' : 'the alkane R–H') + ' (the base’s conjugate acid)' }); else by.push({ smiles: 'C', name: 'the alkane (the base’s conjugate acid)' });
      } else if (m.bs.kind.indexOf('LDA') >= 0) by.push({ smiles: 'CC(C)NC(C)C', name: 'diisopropylamine (the base’s conjugate acid)' });
      else by.push({ smiles: m.bs.kind.indexOf('amide') >= 0 ? 'N' : '[H][H]', name: m.bs.conj + ' (the base’s conjugate acid)' });
      changed.add(hAt); label = (m.ah.enolate ? 'enolate' : m.ah.what.replace(/ .*/, '') + ' anion') + ' (deprotonation by ' + m.bs.kind + ')'; mechKind = 'acidbase';
      /* name the anion from its parent: ethanol → ethoxide, phenol → phenoxide, ethanoic acid → ethanoate, propyne → propynide */
      try { const parent = bmName(m.Y); const cation = mAt ? (mAt.element === 'Li' ? 'lithium ' : mAt.element === 'Na' ? 'sodium ' : mAt.element === 'Mg' ? 'magnesium ' : 'potassium ') : '';
        let an = null;
        if (/ic acid$/.test(parent)) an = parent.replace(/ic acid$/, 'ate'); else if (/carboxylic acid$/.test(parent)) an = parent.replace(/carboxylic acid$/, 'carboxylate');
        else if (parent === 'phenol') an = 'phenoxide'; else if (/anol$/.test(parent)) an = parent.replace(/anol$/, 'oxide'); else if (/ol$/.test(parent)) an = parent.replace(/ol$/, 'olate'); else if (/yne$/.test(parent)) an = parent.replace(/yne$/, 'ynide');
        else if (/phenol$/.test(parent)) an = parent.replace(/phenol$/, 'phenoxide');
        else if (parent === 'aniline') an = 'anilide'; else if (/amine$/.test(parent)) an = parent.replace(/amine$/, 'amide'); else if (/aniline$/.test(parent)) an = parent.replace(/aniline$/, 'anilide');
        if (m.ah.enolate) an = 'enolate of ' + parent;
        if (an) nameHint = cation + an + (cation ? '' : ' ion'); } catch (e) {}
      why = 'An acid–base step: ' + m.bs.kind + ' is a far stronger base than the ' + m.ah.what + ' is an acid (pKa ≈ ' + m.ah.pka + '), so the proton is removed completely and the anion is left, ready to act as a nucleophile in the next step (with an alkyl halide, a carbonyl, an epoxide …).' + (m.ah.enolate ? ' The negative charge is drawn on carbon, but it is shared with the carbonyl oxygen — this is the ENOLATE, the nucleophile of the aldol, Claisen and alkylation reactions. A bulky base (LDA) at −78 °C takes the LESS hindered α-hydrogen (the kinetic enolate).' : '');
      eqn = { R: [[m.bs.conj === 'the alkane' ? 'RH' : m.bs.conj === 'the terminal alkyne' ? 'RC≡CH' : m.bs.conj === 'ammonia' ? 'NH3' : m.bs.conj === 'the secondary amine' ? 'iPr2NH' : 'H2', m.bs.conj]] };
    } else if (m.join === 'aldol') {
      const al = idX(m.en.alpha.id), cE = idY(m.co.c.id), oE = idY(m.co.o.id);
      bmBond(g, al, cE, 1); const bd = g.bonds.find(b => (b.a === cE && b.b === oE) || (b.b === cE && b.a === oE)); if (bd) bd.order = 1;
      changed.add(al); changed.add(cE); changed.add(oE);
      const aldolG = clone(g);
      /* a β-hydroxy carbonyl dehydrates when the alkene it gives would be CONJUGATED — to an
         aryl ring or a second carbonyl — because that conjugation is what pays for losing the
         alcohol. Heat does the same job for an ordinary aliphatic pair. */
      const conjAld = nb(g, cE).some(n => n.atom.id !== al && (isAromatic(g, n.atom.id) || bmDblO(g, n.atom.id)));
      if (heat || conjAld) { /* dehydrate to the enone */
        const bAl = g.bonds.find(b => (b.a === al && b.b === cE) || (b.b === al && b.a === cE)); if (bAl) bAl.order = 2;
        C().removeAtom(g, oE); by.push({ smiles: 'O', name: 'water' });
        label = 'α,β-unsaturated carbonyl (aldol condensation)'; stages = [{ label: 'the aldol (β-hydroxy carbonyl) — before it loses water', graph: bmTidy(aldolG) }];
      } else label = 'β-hydroxy carbonyl (aldol addition)';
      mechKind = 'aldol';
      if (conjAld && !heat) nameHint = null;
      why = 'Base removes an α-hydrogen from the ' + m.en.kind + ' to give the enolate, a carbon nucleophile; it attacks the carbonyl carbon of the ' + m.co.kind + ' and the alkoxide is protonated by the solvent. ' + ((heat || conjAld) ? 'The β-hydroxy carbonyl then loses water by E1cB — deprotonate the α-carbon, expel the hydroxide from the β one — and the C=C that results is CONJUGATED with the carbonyl' + (conjAld && !heat ? ' AND with the ring or second carbonyl next to it, which is why this one condenses without needing heat.' : ', which is what pays for it.') : 'Cold and with an aliphatic partner it stops at the β-hydroxy carbonyl; heating would dehydrate it to the enone.') + ' If both partners have α-hydrogens you get four products — a useful crossed aldol needs one partner with NO α-hydrogen, such as benzaldehyde.';
      eqn = { R: (heat || conjAld) ? [['H2O', 'water — lost in the dehydration']] : [], cat: ['NaOH cat.' + (heat ? ', Δ' : ', 0 °C'), 'the base that makes the enolate; regenerated'] };
    } else if (m.join === 'acAlkyl') {
      const al = idX(m.ac.c.id), cX = idY(m.h.c.id);
      if (m.h.sulf) { bmBond(g, al, cX, 1); const f = cutFragment(g, cX, idY(m.h.x.id)); by.push({ graph: bmTidy(f), name: 'the ' + m.h.lgWord + ' anion' }); }
      else { C().removeAtom(g, idY(m.h.x.id)); bmBond(g, al, cX, 1); by.push({ smiles: '[Na+].[' + m.h.x.element + '-]', name: 'Na' + m.h.x.element }); }
      changed.add(al); changed.add(cX);
      label = 'α-alkylated (the doubly activated carbon)'; mechKind = 'SN2';
      why = 'This is the alkylation that actually works, and it works because of the pKa. The carbon between two esters has a pKa near 13, so sodium ethoxide deprotonates it COMPLETELY — there is no ketone-strength base needed and none of the elimination that ruins a direct α-alkylation. The anion then does a clean Sₙ2 on the halide. ' +
        'It is step one of the MALONIC ESTER SYNTHESIS. Alkylate, then hydrolyse both esters to the diacid, then heat: a malonic acid loses CO₂ through a six-membered hydrogen-bonded transition state, and what is left is a mono-substituted acetic acid. The second ester group existed only to make this alkylation possible, and then it leaves. ' +
        'The same three steps with ethyl acetoacetate (the ACETOACETIC ESTER SYNTHESIS) give a methyl ketone instead. And with a DIhalide you alkylate twice on the same carbon and make a ring — which is how cyclobutane, cyclopentane and cyclohexane carboxylic acids are built. ' +
        'A second alkylation on the same carbon is possible and often wanted; it needs a second equivalent of base and is slower, because the carbon is now more hindered.';
      eqn = { cat: ['NaOEt (1 equiv), EtOH, then the halide; then H₃O⁺ and heat to decarboxylate', 'the ethoxide deprotonates the doubly activated carbon completely'] };
    } else if (m.join === 'knoev' || m.join === 'michael') {
      const al = idX(m.ac.c.id);
      if (m.join === 'knoev') {
        const cE = idY(m.co.c.id), oE = idY(m.co.o.id);
        bmBond(g, al, cE, 1);
        const addG = clone(g);                                  // the aldol-type adduct, as a stage
        const bd = g.bonds.find(b => (b.a === cE && b.b === oE) || (b.b === cE && b.a === oE));
        if (bd) bd.order = 1;
        try { const t = clone(g); bmTidy(t); stages = [{ label: 'the β-hydroxy adduct — it does not survive, because the alkene that follows is conjugated to both withdrawing groups', graph: t }]; } catch (e) {}
        const b2 = g.bonds.find(b => (b.a === al && b.b === cE) || (b.b === al && b.a === cE));
        if (b2) b2.order = 2;
        C().removeAtom(g, oE);
        by.push({ smiles: 'O', name: 'water' });
        changed.add(al); changed.add(cE);
        label = 'Knoevenagel condensation product (an alkene conjugated to both groups)'; mechKind = 'aldol';
        why = 'A Knoevenagel condensation. The carbon between the two withdrawing groups has a pKa around 11–13 — a hundred million times more acidic than a plain ketone α-H — so a secondary amine (piperidine, β-alanine) or potassium carbonate is base enough; nothing stronger is needed, and that is the practical point of the reaction. ' +
          'The carbanion adds to the ' + m.co.kind + ' carbon, and the β-hydroxy adduct then loses water WITHOUT heating, because the alkene it gives is conjugated to both groups at once. That conjugation is also why the reaction is not reversible in the way a plain aldol is. ' +
          'The amine catalyst often does a second job: it condenses with the aldehyde to an iminium, which is a better electrophile than the aldehyde itself.';
        eqn = { R: [['H2O', 'water — lost in the condensation']], cat: ['piperidine or β-alanine (cat.), EtOH or toluene; or K₂CO₃', 'a weak amine base — the doubly activated C–H needs nothing stronger'] };
      } else {
        const cB = idY(m.en.beta.id), cA = idY(m.en.alpha.id);
        bmBond(g, al, cB, 1);
        const bd = g.bonds.find(b => (b.a === cA && b.b === cB) || (b.b === cA && b.a === cB));
        if (bd) bd.order = 1;
        changed.add(al); changed.add(cB); changed.add(cA);
        label = '1,4-adduct (Michael addition)'; mechKind = 'gilman';
        why = 'The Michael addition. The stabilised carbanion is a SOFT nucleophile and adds to the β-carbon of the enone (1,4-, conjugate addition) rather than to the C=O — the opposite of what a Grignard does, and for a real reason: 1,4-addition is under thermodynamic control and forms the stronger C–C bond while keeping the C=O intact, and a soft nucleophile has time to find the softer site because its addition is reversible. ' +
          'The enolate that results picks up a proton from the solvent at the α-carbon. Because the nucleophile is doubly activated, a catalytic amount of a weak base (NaOEt, K₂CO₃, or an amine) is all that is required, and the same carbon can go on to close a ring — a Michael addition followed by an intramolecular aldol is the Robinson annulation.';
        eqn = { cat: ['NaOEt or K₂CO₃ (cat.), EtOH', 'a catalytic weak base; the doubly activated C–H needs no more'] };
      }
    } else if (m.join === 'wittig' || m.join === 'hwe') {
      const cE = idY(m.co.c.id), oE = idY(m.co.o.id);
      const al = idX(m.join === 'wittig' ? m.y.c.id : m.ph.c.id);
      const pId = idX(m.join === 'wittig' ? m.y.p.id : m.ph.p.id);
      /* the new C=C first, then phosphorus is cut away with the old carbonyl oxygen */
      const bd = g.bonds.find(b => (b.a === cE && b.b === oE) || (b.b === cE && b.a === oE));
      if (bd) g.bonds = g.bonds.filter(b => b !== bd);
      bmBond(g, al, cE, 2);
      const pFrag = cutFragment(g, al, pId);
      /* the ylide was drawn as a carbanion; in the alkene that carbon is neutral again */
      const alAt = g.atoms.find(a => a.id === al); if (alAt) { alAt.charge = 0; alAt.hAdd = 0; }
      try { C().removeAtom(g, oE); } catch (e) {}
      changed.add(al); changed.add(cE);
      const stab = m.join === 'hwe' ? true : m.y.stab;
      by.push(m.join === 'wittig'
        ? { smiles: 'O=P(c1ccccc1)(c1ccccc1)c1ccccc1', name: 'triphenylphosphine oxide — the driving force of the whole reaction' }
        : { name: 'the dialkyl phosphate salt — water-soluble, which is why HWE is easier to work up than a Wittig' });
      label = (m.join === 'wittig' ? 'alkene (Wittig reaction)' : 'alkene, E (Horner–Wadsworth–Emmons)');
      mechKind = 'wittig';
      why = (m.join === 'wittig'
        ? 'A Wittig reaction, and the one thing to understand about it is where the driving force comes from: the P=O bond of triphenylphosphine oxide is one of the strongest bonds in organic chemistry, and forming it pays for everything else. ' +
          (m.y.salt ? 'The phosphonium salt drawn here is the precursor — a base (n-BuLi, NaH, or NaOH for a stabilised one) removes the C–H next to phosphorus to give the ylide. ' : '') +
          'The ylide carbon attacks the carbonyl carbon while its oxygen closes onto phosphorus, giving the four-membered oxaphosphetane, which falls apart to the alkene and Ph₃P=O. ' +
          'The C=C ends up exactly where the C=O was — that certainty is why the Wittig is used instead of an elimination, which would leave the position of the double bond in doubt. '
        : 'The Horner–Wadsworth–Emmons reaction: the phosphonate version of the Wittig. The phosphonate C–H is acidic enough for NaH, LiHMDS or even DBU (the two P–O groups help), the anion adds to the carbonyl, and the same four-membered ring collapses to the alkene. ')
        + (m.join === 'hwe'
          ? 'Two practical advantages over the Wittig, and both matter on scale: the by-product is a water-soluble phosphate rather than triphenylphosphine oxide (which is a nuisance to remove by chromatography), and the selectivity is a clean E.'
          : (stab ? 'This ylide is STABILISED (an ester or nitrile on the ylide carbon), so it reacts slowly, reversibly, and gives the E alkene.'
                  : 'This ylide is UNSTABILISED, so it reacts fast and irreversibly through the cis-oxaphosphetane and gives the Z alkene — which is the reason to choose a Wittig over an HWE when Z is what is wanted.'));
      stereoHow = m.join === 'hwe' || stab
        ? 'E-selective: the reaction is reversible enough to pass through the lower-energy trans-oxaphosphetane.'
        : 'Z-selective: an unstabilised ylide adds irreversibly through the cis-oxaphosphetane, and the alkene keeps that geometry.';
      eqn = { cat: [m.join === 'wittig' ? (m.y.salt ? 'n-BuLi or NaH (to make the ylide), THF' : 'THF or CH₂Cl₂, rt') : 'NaH or LiHMDS or DBU, THF',
        m.join === 'wittig' ? 'the base that makes the ylide' : 'a mild base — the phosphonate C–H is acidic enough'] };
    } else if (m.join === 'hetcyc') {
      const bn = m.bn, be = m.be;
      const n1 = idX(bn.n1.id), n2 = idX(bn.n2.id);
      /* which nucleophile goes to which carbon. With a thiourea the SULFUR takes the carbon
         that carries the halide (soft nucleophile, soft electrophile) and the nitrogen takes
         the carbonyl — the Hantzsch regiochemistry. Otherwise the more nucleophilic nitrogen
         takes the more electrophilic carbon. */
      let cA = idY(be.c1.id), cB = idY(be.c2.id);
      if (bn.kind === 'thiourea' && be.kind === 'alpha-halo ketone') { cA = idY(be.c1.id); cB = idY(be.c2.id); }
      const lgOut = [];
      const solo = (cA === cB);                    // one carbon taking both nucleophiles
      const dropLG = (e, cid) => {
        if (!e) return;
        if (e.kind === 'halide') { try { C().removeAtom(g, idY(e.lg.id)); } catch (err) {} lgOut.push({ smiles: e.lg.element, name: 'H' + e.lg.element + ' — taken up by the base' }); return; }
        const acyl = e.kind === 'ester' || e.kind === 'acid' || e.kind === 'acyl chloride';
        if (e.lg) {
          if (e.kind === 'ester') { const frag = cutFragment(g, cid, idY(e.lg.id)); lgOut.push({ graph: bmTidy(frag), name: 'the alcohol (from the ester)' }); }
          else { try { C().removeAtom(g, idY(e.lg.id)); } catch (err) {} lgOut.push({ smiles: e.kind === 'acyl chloride' ? 'Cl' : 'O', name: e.kind === 'acyl chloride' ? 'HCl' : 'water' }); }
        }
        /* the C=O goes only when this carbon has to become part of the aromatic ring itself.
           An ester or acid carbon that takes one nucleophile becomes an AMIDE and keeps its
           oxygen — that is the 4-oxo group of a pyrimidin-4(3H)-one. */
        if (e.o && (solo || !acyl)) { try { C().removeAtom(g, idY(e.o.id)); } catch (err) {} lgOut.push({ smiles: 'O', name: 'water' }); }
      };
      /* the atoms between the two electrophilic carbons, read off the PARTNER molecule before
         anything is bonded (afterwards the shortest path would run through the new nitrogen) */
      let midIds = [];
      if (!solo && be.span >= 1) { try { midIds = bondRingPath(m.Y, be.c1.id, be.c2.id).map(idY); } catch (e) { midIds = []; } }
      /* a nucleophile double-bonded to its own carbon must give that bond up first */
      [bn.n1, bn.n2].forEach(nu => {
        if (!bn.mid) return;
        const bd = g.bonds.find(b => (b.a === idX(nu.id) && b.b === idX(bn.mid.id)) || (b.b === idX(nu.id) && b.a === idX(bn.mid.id)));
        if (bd && bd.order === 2) bd.order = 1;
      });
      /* bond first, then expel — the same order the cross-couplings had to learn */
      bmBond(g, n1, cA, 1);
      if (!(bn.span < 0 && cA === cB)) bmBond(g, n2, cB, 1);
      /* an orthoester keeps its C–H and loses every OR; a nitrile keeps its nitrogen */
      if (be.kind === 'orthoester') {
        be.ors.forEach(o => { const frag = cutFragment(g, cA, idY(o.id)); lgOut.push({ graph: bmTidy(frag), name: 'the alcohol (one of three from the orthoester)' }); });
      } else if (be.kind === 'mono-nitrile') {
        const nid = idY(be.nitrile.id);
        const bd = g.bonds.find(b => (b.a === cA && b.b === nid) || (b.b === cA && b.a === nid));
        if (bd) bd.order = 1;
      } else if (solo) {
        dropLG(be.e1, cA);
      } else {
        dropLG(be.e1, cA); dropLG(be.e2, cB);
      }
      /* the ring, in order, so it can be Kekulé-ed */
      const ring = [];
      ring.push(n1);
      if (bn.span === 2 && bn.path) { ring.push(idX(bn.path[0].id)); ring.push(idX(bn.path[1].id)); }
      else if (bn.span === 1 && bn.mid) ring.push(idX(bn.mid.id));
      if (n2 !== n1) ring.push(n2);
      if (cB !== cA) {
        ring.push(cB);
        midIds.slice().reverse().forEach(id => ring.push(id));
      }
      if (ring[ring.length - 1] !== cA) ring.push(cA);
      const uniq = ring.filter((v, i) => ring.indexOf(v) === i);
      let nPi = 0;
      if (uniq.length === m.size) nPi = aromatiseRing(g, uniq);
      lgOut.forEach(b => by.push(b));
      /* an ALDEHYDE (or a ketone) leaves the ring one oxidation state short of aromatic — the
         two extra hydrogens have to be taken by something, and in practice that is air, sodium
         metabisulfite, DDQ or a little iodine. Say so rather than pretending. */
      const needsOx = (be.kind === 'mono-aldehyde' || be.kind === 'mono-ketone');
      if (needsOx && nPi) by.push({ name: 'two hydrogens — taken by the oxidant (air, Na₂S₂O₅, DDQ or I₂), which is what makes the last step of this sequence an OXIDATIVE cyclisation' });
      uniq.forEach(id => changed.add(id));
      const nm = ringNameOf(g, uniq, bn.arom);
      label = nm + ' (condensation — the ring is built from both molecules)';
      mechKind = 'imine';
      why = 'Two molecules, one new ring. ' + bn.word.charAt(0).toUpperCase() + bn.word.slice(1) + ' meets ' + be.word + ', and each nucleophile adds to a carbon. ' +
        'The mechanism is nothing new — it is nucleophilic addition to a carbonyl, twice: add, then throw out ' +
        (lgOut.some(x => x.name && x.name.indexOf('water') >= 0) ? 'water' : 'the leaving group') + ', then do it again with the second nucleophile. ' +
        'What makes it go to completion is the RING: a ' + m.size + '-membered ring closes with almost no strain, and ' +
        (nPi ? 'the product is AROMATIC — every step before it is reversible, and aromatisation is the step that is not, which is why these condensations need nothing but heat and a trace of acid. '
             : 'the ring here cannot aromatise, so the product is the saturated heterocycle and the reaction stays reversible — keep removing the water. ') +
        (needsOx ? 'One honest caveat: starting from an aldehyde the immediate product is the non-aromatic ring (the aminal), which is TWO HYDROGENS away from the aromatic one. Something has to take them — air on standing, or sodium metabisulfite, DDQ or iodine if you want it in an hour rather than a day. Starting from the carboxylic acid, the ester or the orthoester instead, no oxidant is needed at all, because those carbons arrive already at the right oxidation state. ' : '') +
        (bn.kind === 'thiourea' ? 'The regiochemistry is set by which atom attacks which: sulfur is the softer, more polarisable nucleophile and takes the C–X carbon, while nitrogen takes the carbonyl. That is why a Hantzsch thiazole synthesis gives one isomer and not a mixture. ' : '') +
        (bn.kind === 'hydrazine' && be.span === 1 ? 'With an unsymmetrical 1,3-diketone the two possible pyrazoles are NOT formed equally: the more nucleophilic nitrogen of the hydrazine attacks the more electrophilic (less hindered, less conjugated) carbonyl first, and that first step decides the regiochemistry. ' : '');
      eqn = { cat: ['EtOH or AcOH, 80–120 °C (a trace of acid, water removed)', 'heat and a little acid — the ring closure is what drives it, not the reagent'] };
    } else if (m.join === 'arhal') {
      /* the ring position comes from the same directing-effect scoring the EAS family uses */
      const sites = areneSites(m.Y, 'ar_br2', env).filter(x => !x.inert).sort((p, q) => q.score - p.score);
      if (!sites.length) return null;
      const ringC = idY(sites[0].c.id);
      const xId = idX(m.nx.x.id);
      g.bonds = g.bonds.filter(b => b.a !== xId && b.b !== xId);
      const spent = cutFragment(g, xId, idX(m.nx.n.id));
      bmBond(g, xId, ringC, 1);
      const xa = g.atoms.find(a => a.id === xId); if (xa) { xa.charge = 0; xa.hAdd = 0; }
      by.push({ graph: bmTidy(spent), name: m.nx.word.indexOf('sulfuryl') >= 0 ? 'SO₂ and HCl' : 'succinimide (the spent reagent)' });
      changed.add(ringC); changed.add(xId);
      label = 'aryl ' + (m.nx.el === 'Br' ? 'bromide' : m.nx.el === 'Cl' ? 'chloride' : 'iodide') + ' (aromatic halogenation)';
      mechKind = 'eas';
      why = 'Electrophilic aromatic substitution, with ' + m.nx.word + ' as the halogen source. The ring attacks the halogen, the σ-complex (arenium ion) loses its hydrogen to the succinimide anion that has just left, and aromaticity is restored. ' +
        'Why not bromine and a Lewis acid? Because on an electron-rich ring — a pyrrole, an indole, an anisole, an aniline — Br₂/FeBr₃ is far too reactive: it halogenates two or three times over and oxidises what is left. The N–X bond of a succinimide is weak and polarised, so the ring only ever meets a trace of X⁺, and the reaction stops after one substitution. ' +
        'The position (' + (sites[0].label || 'the best position') + ') is decided by the groups already on the ring: ' + (sites[0].reason || 'the usual directing effects') + '.';
      eqn = { cat: ['DMF or MeCN or CCl₄, 0 °C → rt (in the dark)', 'a polar solvent helps; no Lewis acid, and no light — light would start a radical side-chain reaction instead'] };
    } else if (m.join === 'thion') {
      const ac = m.ac || m.co;
      const cId = idY(ac.c.id), oId = idY(ac.o.id);
      const oAt = g.atoms.find(a => a.id === oId);
      if (oAt) oAt.element = 'S';
      changed.add(cId); changed.add(oId);
      /* the phosphorus reagent leaves as its oxide; cut it away from the sulfur it gave up */
      try {
        const sGone = idX(m.th.s.id);
        const rest = cutFragment(g, idX(m.th.p.id), sGone);
        by.push({ graph: bmTidy(rest), name: 'the phosphorus reagent, now carrying the oxygen (P=S becomes P=O)' });
        const drop = new Set([idX(m.th.p.id)]);
        const st = [idX(m.th.p.id)];
        while (st.length) { const x = st.pop(); nb(g, x).forEach(n => { if (!drop.has(n.atom.id)) { drop.add(n.atom.id); st.push(n.atom.id); } }); }
        g.atoms = g.atoms.filter(a => !drop.has(a.id));
        g.bonds = g.bonds.filter(b => !drop.has(b.a) && !drop.has(b.b));
      } catch (e) {}
      try { C().layoutGraph(g); } catch (e) {}
      label = 'thioamide / thioketone (C=O → C=S)'; mechKind = 'acylsub';
      why = 'Thionation. Lawesson’s reagent (and P₄S₁₀, which is the older and cruder version) exists in solution as a small amount of a reactive dithiophosphine ylide; that adds across the C=O to give a four-membered ring with phosphorus, oxygen, carbon and sulfur, and the ring then collapses the other way — taking the oxygen onto phosphorus and leaving sulfur on carbon. The driving force is the P=O bond, exactly as in a Wittig. ' +
        'The selectivity is worth remembering: amide ≫ ketone > ester, and a carboxylic acid does not work at all. So in a molecule with both an amide and an ester, one equivalent touches only the amide. ' +
        'And this is done for a reason, not for novelty: a thioamide is the standard way into a thiazole (Hantzsch), and swapping one oxygen for sulfur changes hydrogen bonding, pKa and metabolism, which is why it appears so often in medicinal chemistry.';
      eqn = { cat: ['toluene or THF, 60–110 °C', 'heat; one equivalent of the reagent supplies two sulfurs'] };
    } else if (m.join === 'stepPoly') {
      /* Build the alternating chain by repeatedly joining a fresh copy of each monomer. Four
         links for an amide (the slides draw it acid-terminated), five for an ester, carbonate
         or urethane (drawn with the nucleophile at both ends) — enough to see the pattern. */
      const linkAmide = m.kind === 'amide';
      const acidG = m.elec, nucG = m.nuc;
      const chain = { atoms: [], bonds: [], nextId: 1 };
      const nUnits = linkAmide ? 4 : 5;
      let openEnd = null;                                 // the id on the chain waiting to join
      const usedEnds = [];
      for (let i = 0; i < nUnits; i++) {
        const isNuc = (i % 2 === 0);
        const piece = clone(isNuc ? nucG : acidG);
        const map = mergeGraph(chain, piece);
        /* the two ends of this monomer, in the pasted numbering */
        let ends = [];
        if (isNuc) {
          if (m.kind === 'urethane') {
            /* the diol's two O–H oxygens */
            piece.atoms.forEach(a => { if (a.element === 'O' && hOn(piece, a) === 1) ends.push(map.get(a.id)); });
          } else {
            piece.atoms.forEach(a => { if ((a.element === 'N' && hOn(piece, a) >= 1 && !nb(piece, a.id).some(n => bmDblO(piece, n.atom.id))) || (a.element === 'O' && hOn(piece, a) === 1 && !nb(piece, a.id).some(n => bmDblO(piece, n.atom.id)))) ends.push(map.get(a.id)); });
          }
        } else if (m.kind === 'urethane') {
          piece.atoms.forEach(a => { if (a.element === 'N' && nb(piece, a.id).some(n => n.atom.element === 'C' && n.bond.order === 2 && bmDblO(piece, n.atom.id))) ends.push(map.get(a.id)); });
        } else {
          piece.atoms.forEach(a => { if (a.element !== 'C' || !bmDblO(piece, a.id)) return;
            const lg = nb(piece, a.id).find(n => (n.atom.element === 'O' && n.bond.order === 1 && hOn(piece, n.atom) >= 1) || n.atom.element === 'Cl');
            if (lg) ends.push(map.get(a.id)); });
        }
        ends = ends.filter((v, k) => ends.indexOf(v) === k);
        if (ends.length < 1) break;
        if (openEnd !== null) {
          /* join: the nucleophile's heteroatom to the electrophile's carbon, and the
             electrophile loses its leaving group (water, HCl) */
          const nucId = isNuc ? ends[0] : openEnd;
          const acId = isNuc ? openEnd : ends[0];
          if (m.kind === 'urethane') {
            /* an isocyanate ADDS — nothing is lost at all, which is why polyurethane foams
               are made by simply mixing the two components */
            const nId = isNuc ? openEnd : ends[0];
            const oId = isNuc ? ends[0] : openEnd;
            const cAt = nb(chain, nId).find(n => n.atom.element === 'C' && n.bond.order === 2);
            if (cAt) { const bd = bondOf(chain, nId, cAt.atom.id); if (bd) bd.order = 1; bondOf2(chain, cAt.atom.id, oId, 1); }
          } else {
            const lg = nb(chain, acId).find(n => (n.atom.element === 'O' && n.bond.order === 1 && hOn(chain, n.atom) >= 1) || n.atom.element === 'Cl');
            if (lg) { try { cutFragment(chain, acId, lg.atom.id); } catch (e) {} }
            bondOf2(chain, acId, nucId, 1);
            const nA = chain.atoms.find(x => x.id === nucId); if (nA) nA.hAdd = 0;
          }
          ends = ends.filter(e => e !== ends[0]);
        }
        openEnd = ends.length ? ends[ends.length - 1] : null;
        if (openEnd === null) break;
      }
      /* whatever acyl chloride is left at the end of the chain meets water at the workup */
      chain.atoms.slice().forEach(a => {
        if (a.element !== 'Cl') return;
        const c = nb(chain, a.id).find(n => n.atom.element === 'C' && bmDblO(chain, n.atom.id));
        if (!c) return;
        a.element = 'O'; a.hAdd = 0;
      });
      try { C().layoutGraph(chain); } catch (e) {}
      g.atoms = chain.atoms; g.bonds = chain.bonds; g.nextId = chain.nextId;
      chain.atoms.forEach(a => changed.add(a.id));
      const linkWord = m.kind === 'amide' ? 'polyamide (nylon / aramid)' : m.kind === 'urethane' ? 'polyurethane' : 'polyester (or polycarbonate)';
      label = 'the ' + linkWord + ' — a few repeat units shown';
      mechKind = 'acylsub';
      why = 'Step growth, and the contrast with chain growth is worth holding on to. Here every molecule in the flask is reactive at BOTH ends, so any two pieces can join: two monomers, a monomer and a dimer, two dimers. The consequence is that the molecular weight climbs only at the very END of the reaction — at 95 % conversion the average chain is still only about 20 units long — which is why a polycondensation has to be driven almost to completion, with the water or HCl removed continuously. ' +
        (m.kind === 'amide'
          ? 'Each link is an amide, made exactly as any amide is: the amine attacks the acid (as its ammonium carboxylate salt on heating, or far more easily from the acid CHLORIDE) and water or HCl leaves. Nylon 6,6 is adipic acid with hexamethylenediamine — the two sixes are the carbon counts. Kevlar is the same chemistry with two aromatic monomers, and its strength comes from the flat, fully hydrogen-bonded sheets those rigid rings allow.'
          : m.kind === 'urethane'
            ? 'Each link is a carbamate, and notice what is NOT lost: an isocyanate is already at the right oxidation state, so the alcohol simply ADDS across the C=N. Nothing is expelled, no water has to be removed, and the two components can be mixed and poured — which is why polyurethane is the foam in your furniture. Add a little water on purpose and it makes CO₂ with some of the isocyanate, and that is the blowing agent.'
            : 'Each link is an ester, made by the same Fischer chemistry as any ester, with the water removed to drive it. PET is terephthalic acid with ethylene glycol. A polycarbonate is the same idea with phosgene or a carbonate, giving O–CO–O links and a rigid, transparent, tough material.') + ' ' +
        'The number of units drawn here is illustrative — a real chain is thousands long — but the LINKAGE and the alternation are exactly right.';
      eqn = { R: m.kind === 'urethane' ? [] : [['H2O', 'water (or HCl from the acid chloride) — removed continuously to drive the chain out']], cat: [m.kind === 'urethane' ? 'mix and pour; a tin or amine catalyst' : 'heat, with the water distilled out (or the diacid chloride with a base)', 'the condensation has to be driven nearly to completion before the chains get long'] };
    } else if (m.join === 'metalInsert') {
      const cId = idY(m.h.c.id), xId = idY(m.h.x.id), mId = idX(m.mt.m.id);
      const el = m.mt.el, xEl = m.h.x.element;
      const mAt = g.atoms.find(a => a.id === mId); if (mAt) mAt.charge = 0;
      g.bonds = g.bonds.filter(b => !((b.a === cId && b.b === xId) || (b.b === cId && b.a === xId)));
      bmBond(g, cId, mId, 1);
      if (el === 'Li') {
        /* R–X + 2 Li → R–Li + LiX: the halide leaves with the second lithium */
        try { C().removeAtom(g, xId); } catch (e) {}
        by.push({ smiles: '[Li+].[' + xEl + '-]', name: 'Li' + xEl + ' (from the second equivalent of lithium)' });
      } else {
        bmBond(g, mId, xId, 1);                        // R–Mg–X: the metal sits BETWEEN carbon and halogen
      }
      changed.add(cId); changed.add(mId);
      const word = el === 'Mg' ? 'Grignard reagent' : el === 'Li' ? 'organolithium' : 'organozinc';
      label = word + ' (R–' + (el === 'Mg' ? 'MgX' : el === 'Li' ? 'Li' : 'ZnX') + ' — the metal has inserted into the C–' + xEl + ' bond)';
      mechKind = 'grignard';
      /* the name: R-yl + magnesium bromide / lithium, for the cases where "R-yl" is not in doubt */
      try {
        const rh = clone(m.Y); C().removeAtom(rh, m.h.x.id);
        const parent = (C().displayNameFor ? C().displayNameFor(rh) : '') || '';
        const heavy = rh.atoms.filter(a => a.element !== 'H').length;
        let ryl = null;
        if (parent === 'benzene') ryl = 'phenyl';
        else if (/^cyclo[a-z]+ane$/.test(parent) && heavy === rh.atoms.length && rh.atoms.every(a => a.element === 'C')) ryl = parent.replace(/ane$/, 'yl');
        else if (/^(meth|eth|prop|but|pent|hex|hept|oct)ane$/.test(parent) && m.h.deg <= 1 && !m.h.sp2) ryl = parent.replace(/ane$/, 'yl');
        else if (parent === 'ethene' || parent === 'ethylene') ryl = 'vinyl';
        if (ryl) nameHint = ryl + (el === 'Mg' ? 'magnesium ' + ({ Cl: 'chloride', Br: 'bromide', I: 'iodide' }[xEl] || 'halide') : el === 'Li' ? 'lithium' : 'zinc ' + ({ Cl: 'chloride', Br: 'bromide', I: 'iodide' }[xEl] || 'halide'));
      } catch (e) {}
      why = 'Oxidative insertion. ' + (el === 'Mg' ? 'Magnesium metal inserts itself between carbon and halogen: the metal gives up two electrons, and what was an electrophilic carbon (C–X, polarised δ+ on carbon) becomes a NUCLEOPHILIC one (C–Mg, polarised δ− on carbon). That reversal of polarity — umpolung — is the entire point of making a Grignard reagent: it turns an alkyl halide into a carbanion equivalent that will add to a carbonyl, an epoxide, a nitrile or CO₂. ' : el === 'Li' ? 'Two lithium atoms are consumed for each C–X bond: one becomes the C–Li bond, the other leaves as the lithium halide. The carbon goes from δ+ to δ− — this is the same polarity reversal as a Grignard, but an organolithium is more reactive and more basic. ' : 'Zinc inserts into the C–X bond to give an organozinc, the mild organometallic used in Negishi couplings and Reformatsky reactions. ') +
        'It happens on the metal surface, through single-electron transfer to the C–X σ* orbital, and the surface has to be clean: a crystal of iodine or a little 1,2-dibromoethane is added to scour off the oxide layer, and the reaction often needs a moment of warming to start and then runs exothermically. ' +
        'Two things are non-negotiable. The ether or THF must be DRY — the Grignard reagent is a strong base and one molecule of water destroys one molecule of reagent, giving R–H. And the molecule must carry no acidic hydrogen of its own: an O–H, an N–H, a carboxylic acid or a terminal alkyne would protonate the organometallic the instant it formed.' +
        (m.h.sp2 ? (m.h.aryl ? ' An aryl halide works, though aryl chlorides need THF and heat; the reagent formed is an arylmagnesium halide.' : ' A vinyl halide works too, and the alkene geometry is kept.') : m.h.deg >= 3 ? ' A tertiary halide is fine here — this is not a substitution, so there is no steric problem at carbon.' : '') +
        (m.acidic ? ' In THIS molecule that condition is not met: it carries ' + m.acidic.what + ' (pKa ≈ ' + m.acidic.pka + '), so the organometallic would be quenched by its own substrate as fast as it forms. Protect that group first (a TBS ether, an acetal), or make the Grignard from a different piece.' : '');
      eqn = { cat: [(el === 'Mg' ? 'Mg turnings, dry Et₂O or THF, I₂ (cat.) to start' : el === 'Li' ? 'Li metal (2 equiv), dry Et₂O or pentane, 0 °C' : 'Zn dust, THF, 1,2-dibromoethane (cat.)'), 'anhydrous conditions throughout — water destroys the product'] };
    } else if (m.join === 'diazotise') {
      /* the whole nitrosating molecule is consumed; what stays on the nitrogen is one N */
      const nId = idY(m.am.n.id);
      const nAt = g.atoms.find(a => a.id === nId);
      const nNew = idX(m.ni.n.id);
      /* strip the oxygens off the nitrite's nitrogen and bond it to the amine */
      const oIds = nb(g, nNew).filter(x => x.atom.element === 'O').map(x => x.atom.id);
      oIds.forEach(o => { try { const f = cutFragment(g, nNew, o); by.push({ graph: bmTidy(f), name: 'water (or the alcohol, from an alkyl nitrite)' }); } catch (e) {} });
      bmBond(g, nId, nNew, 1);
      if (m.am.kind === 'secondary') {
        /* no second N–H: it stops at the nitrosamine, N–N=O */
        const oNew = C().addAtomNear(g, 'O', g.atoms.find(a => a.id === nNew));
        const bd = bondOf(g, nNew, oNew.id); if (bd) bd.order = 2;
        label = 'N-nitrosamine'; mechKind = 'acylsub';
        why = 'The amine attacks NO⁺ (which is what a nitrite becomes in acid) to give the N-nitroso compound. With a PRIMARY amine this same intermediate carries on — it tautomerises and loses water to give the diazonium — but a secondary amine has no second N–H, so it stops here. Worth knowing rather than memorising: N-nitrosamines are potent carcinogens, and the sartan and ranitidine recalls came from this reaction happening unintentionally when a secondary amine met a nitrite impurity.';
      } else {
        const bd = bondOf(g, nId, nNew); if (bd) bd.order = 3;
        if (nAt) nAt.charge = 1;
        label = 'aryl diazonium salt'; mechKind = 'acylsub';
        why = 'Diazotisation. The nitrite becomes nitrous acid and then the nitrosonium ion NO⁺; the amine nitrogen attacks it, the N-nitrosamine tautomerises to Ar–N=N–OH, and a last protonation lets water leave. What is left is the diazonium ion, with N₂ on the ring as the best leaving group aromatic chemistry has. ' +
          'Keep it at 0–5 °C and use it where you made it: warmer, it decomposes to the phenol and nitrogen gas, and the dry solid salts are shock-sensitive. From here the ring will take Cl, Br or CN with copper(I) (Sandmeyer), I with plain KI, F by heating the tetrafluoroborate, OH with warm water, or nothing at all — H₃PO₂ removes the group altogether, which is how a directing amine is used and then discarded.';
      }
      changed.add(nId); changed.add(nNew);
      try { C().layoutGraph(g); } catch (e) {}
      eqn = { cat: ['NaNO₂, aq. HCl, 0–5 °C', 'nitrous acid made in the flask; cold, always'] };
    } else if (m.join === 'alkylate' || m.join === 'aldolE' || m.join === 'michaelE' || m.join === 'claisen') {
      /* the enolate's CARBON is the nucleophile; the C=C collapses back to C=O */
      const cNu = idX(m.en.c.id), cCarbonyl = idX(m.en.cO.id), oId = idX(m.en.o.id);
      const restore = () => {
        const bd = g.bonds.find(b => (b.a === cNu && b.b === cCarbonyl) || (b.b === cNu && b.a === cCarbonyl));
        if (bd) bd.order = 1;
        const bo = g.bonds.find(b => (b.a === cCarbonyl && b.b === oId) || (b.b === cCarbonyl && b.a === oId));
        if (bo) bo.order = 2;
        const oa = g.atoms.find(a => a.id === oId); if (oa) { oa.charge = 0; oa.hAdd = 0; }
      };
      if (m.join === 'alkylate') {
        const cX = idY(m.h.c.id);
        if (m.h.sulf) { bmBond(g, cNu, cX, 1); const f = cutFragment(g, cX, idY(m.h.x.id)); by.push({ graph: bmTidy(f), name: 'the ' + m.h.lgWord + ' anion' }); }
        else { C().removeAtom(g, idY(m.h.x.id)); bmBond(g, cNu, cX, 1); by.push({ smiles: '[Na+].[' + m.h.x.element + '-]', name: 'Na' + m.h.x.element }); }
        restore();
        changed.add(cNu); changed.add(cX);
        label = 'α-alkylated carbonyl'; mechKind = 'SN2';
        why = 'α-Alkylation: the enolate carbon attacks the alkyl halide by Sₙ2, and the carbonyl comes back. In one step you have made a new C–C bond next to a carbonyl, which is why this looks like the obvious way to build a molecule. ' +
          'In practice it is fussy, and it is worth knowing why. The enolate is also a strong BASE, so with anything but a methyl, primary or benzylic halide it does E2 instead of substitution — a secondary halide gives mostly alkene and a tertiary one gives nothing else. The product ketone is also still acidic, so it can be deprotonated and alkylated again (polyalkylation), and with an unsymmetrical ketone you have to have chosen the right enolate to begin with. ' +
          'Those three problems are exactly what the malonic ester synthesis and the Stork enamine reaction were invented to get around.';
        eqn = { cat: ['LDA or NaOEt first, then the halide; THF', 'the base makes the enolate; the halide must be methyl, primary or benzylic'] };
      } else if (m.join === 'aldolE') {
        const cE = idY(m.co.c.id), oE = idY(m.co.o.id);
        bmBond(g, cNu, cE, 1);
        const bd = g.bonds.find(b => (b.a === cE && b.b === oE) || (b.b === cE && b.a === oE)); if (bd) bd.order = 1;
        const oe = g.atoms.find(a => a.id === oE); if (oe) oe.charge = 0;
        restore();
        changed.add(cNu); changed.add(cE); changed.add(oE);
        const addG = clone(g);
        /* does the aldol dehydrate? It does when the alkene that results is conjugated to
           something — an aryl ring or a second carbonyl — or when the mixture is heated. */
        const conj = nb(g, cE).some(n => n.atom.id !== cNu && (isAromatic(g, n.atom.id) || bmDblO(g, n.atom.id)));
        if (heat || conj) {
          const b2 = g.bonds.find(b => (b.a === cNu && b.b === cE) || (b.b === cNu && b.a === cE)); if (b2) b2.order = 2;
          try { C().removeAtom(g, oE); } catch (e) {}
          by.push({ smiles: 'O', name: 'water' });
          label = 'α,β-unsaturated carbonyl (aldol CONDENSATION)';
          try { const t = clone(addG); bmTidy(t); stages = [{ label: 'the aldol (β-hydroxy carbonyl) — before it loses water', graph: t }]; } catch (e) {}
        } else label = 'β-hydroxy carbonyl (aldol ADDITION)';
        mechKind = 'aldol';
        why = 'The aldol reaction: an enolate adds to a carbonyl. Count the carbons and you can see the design — the α-carbon of one molecule bonds to the carbonyl carbon of the other, so the product is a β-hydroxy carbonyl, and the new bond is always between those two positions. ' +
          ((heat || conj)
            ? 'It does not stop there. The β-hydroxy carbonyl still has an acidic α-hydrogen, and losing water from it gives a C=C that is CONJUGATED to the carbonyl — a genuinely more stable arrangement. So the reaction runs on to the enone, by E1cB: deprotonate the α-carbon, then expel hydroxide from the β one. ' +
              (conj ? 'Here the alkene ends up conjugated with a ring or a second carbonyl as well, which is why the dehydration happens even without heating.' : 'Heat is what carries it over.')
            : 'Cold and with an aliphatic partner it stops at the β-hydroxy carbonyl. Warm it and it dehydrates to the enone.') + ' ' +
          'One practical warning that the textbook makes into an exam question: if BOTH partners have α-hydrogens and both are electrophilic, you get four products. A useful crossed aldol needs one partner with no α-hydrogen (benzaldehyde, formaldehyde) or a preformed enolate added to the other partner in a separate step.';
        eqn = { R: (heat || conj) ? [['H2O', 'water — lost in the condensation']] : [], cat: ['NaOH or NaOEt (cat.)' + ((heat || conj) ? ', Δ' : ', 0 °C'), 'the base makes the enolate and is regenerated'] };
      } else if (m.join === 'michaelE') {
        const cB = idY(m.e2.beta.id), cA = idY(m.e2.alpha.id);
        bmBond(g, cNu, cB, 1);
        const bd = g.bonds.find(b => (b.a === cA && b.b === cB) || (b.b === cA && b.a === cB)); if (bd) bd.order = 1;
        restore();
        changed.add(cNu); changed.add(cB); changed.add(cA);
        label = '1,4-adduct (Michael addition)'; mechKind = 'gilman';
        why = 'The Michael addition: an enolate adds to the β-carbon of an enone, not to its carbonyl. Both sites are electrophilic — the resonance structures of an enone put positive character on the carbonyl carbon AND on the β-carbon — and which one is attacked depends on the nucleophile. ' +
          'A stabilised enolate is a SOFT, reversible nucleophile: it adds and retreats until it finds the site that gives the more stable product, and 1,4-addition wins because the C=O survives. A hard, irreversible nucleophile such as a Grignard or a hydride grabs the carbonyl instead (1,2-addition). That soft/hard split is the whole story of conjugate addition. ' +
          'A Michael addition followed by an intramolecular aldol condensation is the ROBINSON ANNULATION, which builds a new six-membered ring with an enone in it — the standard route into steroid skeletons.';
        eqn = { cat: ['NaOEt or K₂CO₃ (cat.), EtOH', 'catalytic base: a doubly activated donor needs no more'] };
      } else {
        /* CLAISEN: the enolate attacks an ester carbonyl and the alkoxide leaves */
        const cE = idY(m.es.c.id), orId = idY(m.es.or.id);
        bmBond(g, cNu, cE, 1);
        const f = cutFragment(g, cE, orId);
        by.push({ graph: bmTidy(f), name: 'the alkoxide (ethoxide) — and this is the step that must be irreversible' });
        restore();
        changed.add(cNu); changed.add(cE);
        label = 'β-keto ester (Claisen condensation)'; mechKind = 'acylsub';
        why = 'The Claisen condensation. An ester enolate attacks the carbonyl of a second ester, and the tetrahedral intermediate expels the alkoxide — so unlike an aldol this is an ACYL SUBSTITUTION, and the product is a β-keto ester. ' +
          'Now the part that is always asked about: every step so far is reversible and slightly uphill, so what makes the reaction go? The PRODUCT. A β-keto ester has a hydrogen between two carbonyls with a pKa of about 11, and the ethoxide base (pKa 16) deprotonates it completely. That final, irreversible deprotonation is what pulls the whole equilibrium over, and it is why a FULL equivalent of base is needed — not a catalytic amount — and why the acid workup at the end is not optional: it is what gives you back the neutral β-keto ester. ' +
          'It also explains a rule that otherwise looks arbitrary: the ester must have TWO α-hydrogens. With only one, the product has none between the carbonyls, the final deprotonation is impossible, and the reaction will not go. Run intramolecularly on a diester it is the DIECKMANN condensation and gives a cyclic β-keto ester.';
        eqn = { cat: ['NaOEt (1 full equiv), EtOH; then H₃O⁺', 'a full equivalent, because the product is deprotonated at the end — that is what drives it'] };
      }
    } else if (m.join === 'stork' || m.join === 'storkMichael') {
      const cNu = idX(m.em.c.id), cN = idX(m.em.cN.id), nId = idX(m.em.n.id);
      if (m.join === 'stork') {
        const cX = idY(m.h.c.id);
        C().removeAtom(g, idY(m.h.x.id));
        bmBond(g, cNu, cX, 1);
        changed.add(cNu); changed.add(cX);
        by.push({ smiles: m.h.x.element, name: 'H' + m.h.x.element });
      } else {
        const cB = idY(m.e2.beta.id), cA = idY(m.e2.alpha.id);
        bmBond(g, cNu, cB, 1);
        const bd = g.bonds.find(b => (b.a === cA && b.b === cB) || (b.b === cA && b.a === cB)); if (bd) bd.order = 1;
        changed.add(cNu); changed.add(cB);
      }
      /* the iminium that results is hydrolysed on workup, giving the ketone back */
      const bd2 = g.bonds.find(b => (b.a === cNu && b.b === cN) || (b.b === cNu && b.a === cN));
      if (bd2) bd2.order = 1;
      const f = cutFragment(g, cN, nId);
      const o = C().addAtomNear(g, 'O', g.atoms.find(a => a.id === cN));
      const bo = bondOf(g, cN, o.id); if (bo) bo.order = 2;
      by.push({ graph: bmTidy(f), name: 'the secondary amine, recovered at the hydrolysis' });
      changed.add(cN); changed.add(o.id);
      try { C().layoutGraph(g); } catch (e) {}
      label = 'α-' + (m.join === 'stork' ? 'alkylated' : 'functionalised') + ' ketone (Stork enamine ' + (m.join === 'stork' ? 'alkylation' : 'Michael addition') + ')';
      mechKind = 'SN2';
      why = 'The Stork enamine reaction, and it exists to solve the three problems of direct α-alkylation. An enamine is a NEUTRAL carbon nucleophile — the nitrogen lone pair pushes electron density onto the far carbon — so there is no strong base in the flask, nothing to cause E2 elimination, and nothing to deprotonate the product and alkylate it twice. ' +
        'Three steps: condense the ketone with a secondary amine (pyrrolidine or morpholine) to the enamine, alkylate or add to a Michael acceptor, then hydrolyse the iminium with aqueous acid to get the ketone back. The amine is recovered, so it is effectively a catalyst. ' +
        'It also chooses its own regiochemistry: an enamine forms towards the LESS substituted side, because the more substituted one is twisted out of conjugation by A-strain with the amine. So the Stork route alkylates the less hindered α-carbon, which is often the one you wanted.';
      eqn = { cat: ['the enamine first (R₂NH, cat. acid, −H₂O), then the halide, then H₃O⁺', 'the amine is recovered at the hydrolysis'] };
    } else if (m.join === 'enamine') {
      const cC = idX(m.co.c.id), oId = idX(m.co.o.id), nAt = idY(m.nu.at.id);
      /* find the alpha carbon that will carry the new C=C */
      let alphaId = null;
      try { const a0 = alphaSites(m.X).filter(a => a.c.id === m.co.c.id).sort((p, q) => p.sub - q.sub)[0]; if (a0) alphaId = idX(a0.alpha.id); } catch (e) {}
      C().removeAtom(g, oId);
      bmBond(g, cC, nAt, 1);
      if (alphaId !== null) { const ab = g.bonds.find(b => (b.a === cC && b.b === alphaId) || (b.b === cC && b.a === alphaId)); if (ab) ab.order = 2; changed.add(alphaId); }
      changed.add(cC); changed.add(nAt);
      by.push({ smiles: 'O', name: 'water' });
      label = 'enamine'; mechKind = 'imine';
      why = 'A SECONDARY amine cannot give an imine: after it adds to the carbonyl and water leaves, the nitrogen has no hydrogen left to lose, so the molecule takes a proton from the α-CARBON instead and the double bond ends up between the two carbons. That is an ENAMINE, N–C=C. ' +
        'The reason it matters is what that double bond does. The nitrogen lone pair is conjugated into it, so the far carbon carries real electron density and behaves as a neutral carbon nucleophile — no base required. That is the basis of the Stork reaction, and of a great deal of modern organocatalysis. ' +
        'The enamine forms towards the LESS substituted α-carbon, because the alternative is twisted out of conjugation, and mild acid catalysis plus removal of the water drives it. Aqueous acid reverses the whole thing and gives the ketone back.';
      eqn = { R: [['H2O', 'water — removed to drive the reaction']], cat: ['cat. TsOH, benzene, −H₂O (Dean–Stark)', 'a trace of acid, with the water taken out'] };
    } else if (m.join === 'imine') {
      const cC = idX(m.co.c.id), oId = idX(m.co.o.id), nAt = idY(m.nu.at.id);
      C().removeAtom(g, oId); bmBond(g, cC, nAt, 2); changed.add(cC); changed.add(nAt);
      by.push({ smiles: 'O', name: 'water' }); label = 'imine (Schiff base)'; mechKind = 'imine';
      why = 'The amine adds to the carbonyl (carbinolamine), then loses water under mild acid to give the C=N. Reversible — remove the water to drive it. With a hydride present (NaBH₃CN) this becomes a reductive amination.';
      eqn = { R: [['H2O', 'water — the BYPRODUCT']], cat: ['H⁺ cat., −H₂O', 'a trace of acid; water removed'] };
    } else if (m.join === 'da') {
      const d = m.d, dp = m.dp;
      const c1 = idX(d.c1.id), c2 = idX(d.c2.id), c3 = idX(d.c3.id), c4 = idX(d.c4.id), a = idY(dp.a.id), b = idY(dp.b.id);
      const setO = (p, q, o) => { const bd = g.bonds.find(x => (x.a === p && x.b === q) || (x.b === p && x.a === q)); if (bd) bd.order = o; };
      /* was the diene held in a ring? then the two ends are already bridged and the
         adduct is BICYCLIC - cyclopentadiene is the case everyone meets first. */
      let bridge = null;
      try {
        const seen = new Set([m.d.c2.id, m.d.c3.id]); const st = [[m.d.c1.id, 0]]; seen.add(m.d.c1.id);
        while (st.length) { const [x, dd] = st.shift();
          if (x === m.d.c4.id) { bridge = dd - 1; break; }
          nb(m.X, x).forEach(n => { if (!seen.has(n.atom.id)) { seen.add(n.atom.id); st.push([n.atom.id, dd + 1]); } }); }
      } catch (e) { bridge = null; }
      const regio = bridge === null ? daRelation(m.d, m.dp, m.X, m.Y) : null;
      setO(c1, c2, 1); setO(c2, c3, 2); setO(c3, c4, 1); setO(a, b, 1);
      bmBond(g, c1, a, 1); bmBond(g, c4, b, 1);
      [c1, c2, c3, c4, a, b].forEach(id => changed.add(id));
      const STEM = { 7: 'hept', 8: 'oct', 9: 'non', 10: 'dec', 11: 'undec', 12: 'dodec' };
      const cage = bridge !== null && bridge >= 1
        ? (bridge === 1 ? 'norbornene \u2014 bicyclo[2.2.1]hept-2-ene' : 'bicyclo[2.2.' + bridge + ']' + (STEM[6 + bridge] || 'cyclo') + '-2-ene')
        : null;
      label = (cage || 'cyclohexene') + ' (Diels–Alder adduct)' + (regio ? ' — ' + regio.rel : ''); mechKind = 'DA';
      why = 'One concerted step: the diene’s 4π electrons and the dienophile’s 2π electrons reorganise through a single cyclic transition state — two new σ-bonds at the ENDS of the diene (C1 and C4) and a new C=C across its middle. ' +
        (dp.ewg ? 'The electron-withdrawing group on the dienophile lowers its LUMO, which is what makes this fast.' : 'Without an electron-withdrawing group on the dienophile it needs heat or a Lewis acid.') +
        (cage ? ' The diene here is locked in a ring, so its two ends were already joined: the piece bridging them stays put and becomes the one-carbon bridge of the cage. That is why cyclopentadiene gives a BRIDGED bicyclic — a norbornene — and never two fused five-membered rings. Counting it out: bridgehead, two-carbon bridge carrying the new C=C, bridgehead, two-carbon bridge from the dienophile, and the original sp³ CH₂ over the top as C7.' : '') +
        (regio ? ' Regiochemistry: the ' + (regio.donorTerminal ? '1-substituted diene gives the "ortho" adduct' : '2-substituted diene gives the "para" adduct') +
          ' — here ' + regio.rel + '. Push the donor’s lone pair through the diene and one terminus carries the charge; the dienophile is polarised the other way, with the carbon AWAY from the withdrawing group left electron-poor. Those two ends find each other, and that is the whole of the ortho–para rule.' : '');
      stereoHow = 'Suprafacial on both partners, so it is stereospecific: whatever was cis on the dienophile stays cis in the ring, and the diene’s termini keep their geometry. ' +
        (dp.ewg ? 'Relative to the ring, the ENDO adduct is the one you isolate: in the transition state the dienophile sits UNDER the diene with its π system overlapping the diene’s C2–C3, which pushes the withdrawing group onto the crowded face. It is the faster product, not the more stable one — heat it long enough and the exo isomer takes over. ' : '') +
        (cage ? 'In the bicyclic adduct that means the substituent points TOWARDS the C=C bridge and AWAY from the one-carbon CH₂ bridge.' : '') +
        ' From achiral partners both faces react equally, so the product is racemic.';
      eqn = { over: 'Δ' };
    } else if (m.join === 'fc-acyl' || m.join === 'fc-alkyl') {
      /* the best ring position, by the same directing-effect scoring as the EAS family */
      const sites = areneSites(A, 'sulfonation', env).filter(s => !s.inert).sort((p, q) => q.score - p.score);
      if (!sites.length) return null;
      const ringC = idX(sites[0].c.id);
      let eC;
      if (m.join === 'fc-acyl') { eC = idY(m.ac.c.id); bmBond(g, ringC, eC, 1); const lg = idY(m.ac.lg.id); if (m.ac.kind === 'acyl chloride') { C().removeAtom(g, lg); by.push({ smiles: 'Cl', name: 'HCl' }); } else { g.bonds = g.bonds.filter(bd => !((bd.a === eC && bd.b === lg) || (bd.b === eC && bd.a === lg))); const other = splitOff(g, ringC); by.push({ graph: bmTidy(other), name: 'the carboxylic acid (other half of the anhydride)' }); } }
      else { eC = idY(m.h.c.id); bmBond(g, ringC, eC, 1); C().removeAtom(g, idY(m.h.x.id)); by.push({ smiles: m.h.x.element, name: 'H' + m.h.x.element }); }
      changed.add(ringC); changed.add(eC);
      label = m.join === 'fc-acyl' ? 'aryl ketone (Friedel–Crafts acylation, ' + (sites[0].rel || 'ring') + ')' : 'alkylarene (Friedel–Crafts alkylation, ' + (sites[0].rel || 'ring') + ')'; mechKind = 'eas';
      why = 'AlCl₃ makes the ' + (m.join === 'fc-acyl' ? 'acylium ion' : 'carbocation') + '; the ring’s π electrons attack it at ' + (sites[0].reason || 'the free position') + ', the σ-complex loses a proton and the ring is aromatic again.' + (m.join === 'fc-acyl' ? ' The ketone deactivates the ring, so it stops after one.' : ' The product is more reactive than the starting arene, so use excess arene.');
      eqn = { R: [[by[0].smiles === 'Cl' ? 'HCl' : by[0].smiles ? 'H' + by[0].smiles : 'RCOOH', 'the BYPRODUCT']], cat: ['AlCl₃', 'the Lewis acid catalyst'] };
    } else return null;
    bmTidy(g);
    /* A reagent written as a salt (a phosphonium chloride, sodium hydroxide) leaves its
       counter-ion in the graph, and a spectator ion is not part of the product. Anything left
       that contains none of the atoms this reaction touched is moved to the by-products. */
    try {
      const keep = new Set();
      const seed = [...changed][0];
      if (seed !== undefined) {
        const st = [seed]; keep.add(seed);
        while (st.length) { const x = st.pop(); nb(g, x).forEach(n => { if (!keep.has(n.atom.id)) { keep.add(n.atom.id); st.push(n.atom.id); } }); }
        const strays = g.atoms.filter(a => !keep.has(a.id)).map(a => a.id);
        if (strays.length && keep.size) {
          const set = new Set(strays);
          const frag = { atoms: g.atoms.filter(a => set.has(a.id)).map(a => Object.assign({}, a)),
            bonds: g.bonds.filter(b => set.has(b.a) && set.has(b.b)).map(b => Object.assign({}, b)), nextId: g.nextId };
          g.atoms = g.atoms.filter(a => keep.has(a.id));
          g.bonds = g.bonds.filter(b => keep.has(b.a) && keep.has(b.b));
          by.push({ graph: bmTidy(frag), name: frag.atoms.length === 1 && frag.atoms[0].charge ? 'the counter-ion (a spectator)' : 'a spectator — it takes no part in this reaction' });
        }
      }
    } catch (e) {}
    return { graph: g, changed, added: [], byprod: by, label, why, mechKind, stages, eqn, nameHint, stereoHow };
  }

  F.bimol = {
    rules: [], needs: 'any',
    sites(g, ruleId, env) {
      const rule = C().RULES.find(r => r.id === ruleId); const B = rule && rule.molecule; if (!B) return [];
      const ms = bmMatches(g, B);
      if (!ms.length) return [];
      /* one site per match kind — the best of each */
      const seen = new Set(); const out = [];
      const fn = ms.find(m => m.join === 'fc-none'); if (fn && !ms.some(m => m.join !== 'fc-none')) return [{ kind: 'join', join: 'fc-none', atoms: [g.atoms[0].id], c: g.atoms[0], label: 'Friedel–Crafts on this ring', inert: fn.why }];
      const gn = ms.find(m => m.join === 'gilman-none'); if (gn && !ms.some(m => m.join !== 'gilman-none')) return [{ kind: 'join', join: 'gilman-none', atoms: [g.atoms[0].id], c: g.atoms[0], label: 'cuprate + simple carbonyl', inert: 'a Gilman reagent (R₂CuLi) is a soft nucleophile — it does 1,4-addition to enones and couples with halides / acyl chlorides, but leaves an isolated C=O alone (use a Grignard for that)' }];
      const HANDLE0 = { acyl: m => [m.ac.c, m.nu.at], sulfonamide: m => [m.su.s, m.nu.at], sn2: m => [m.h.c, m.nu.at], grignard: m => [m.gr.c, m.co.c], aldol: m => [m.en.alpha, m.co.c], imine: m => [m.co.c, m.nu.at], da: m => [m.d.c1, m.dp.a], 'fc-acyl': m => [null, m.ac.c], 'fc-alkyl': m => [null, m.h.c],
        epoxNu: m => [m.nu.at, m.ep.c],
        'grignard-co2': m => [m.gr.c, m.co.c], 'grignard-ester': m => [m.gr.c, m.es.c], 'grignard-nitrile': m => [m.gr.c, m.ni.c], 'grignard-epoxide': m => [m.gr.c, m.ep.c], 'acetylide-sn2': m => [m.gr.c, m.h.c],
        'gilman-14': m => [m.cu.c, m.en.beta], 'gilman-couple': m => [m.cu.c, m.h.c], 'gilman-acyl': m => [m.cu.c, m.ac.c], 'gilman-none': m => [m.cu.c, m.co.c], 'fc-none': m => [null, null], deprotonate: m => [m.bs.c, m.ah.at],
        suzuki: m => [m.h.c, m.p.c], stille: m => [m.h.c, m.p.c], negishi: m => [m.h.c, m.p.c], kumada: m => [m.h.c, m.p.c],
        sonogashira: m => [m.h.c, m.p.c], buchwald: m => [m.h.c, m.p.at], heck: m => [m.h.c, m.p.a], snar: m => [m.h.c, m.p.at],
        charyl: m => [m.h.c, m.p.c],
        oxidS: m => [m.px.o2, m.su.s], oxidN: m => [m.px.o2, m.nn.n], oxidC: m => [m.px.o2, m.d.a],
        knoev: m => [m.ac.c, m.co.c], michael: m => [m.ac.c, m.en.beta], acAlkyl: m => [m.ac.c, m.h.c], wittig: m => [m.y.c, m.co.c], hwe: m => [m.ph.c, m.co.c],
        hetcyc: m => [m.bn.n1, m.be.c1],
        arhal: m => [m.nx.x, null], thion: m => [m.th.s, (m.ac || m.co).c], diazotise: m => [m.ni.n, m.am.n],
        metalInsert: m => [m.mt.m, m.h.c], stepPoly: m => [null, null],
        alkylate: m => [m.en.c, m.h.c], aldolE: m => [m.en.c, m.co.c], michaelE: m => [m.en.c, m.e2.beta],
        claisen: m => [m.en.c, m.es.c], stork: m => [m.em.c, m.h.c], storkMichael: m => [m.em.c, m.e2.beta],
        enamine: m => [m.co.c, m.nu.at] };
      const HANDLE = Object.assign({}, HANDLE0, { 'fc-acyl': m => [m.X.atoms.find(a => bmArom(m.X, a.id)) || m.X.atoms[0], m.ac.c], 'fc-alkyl': m => [m.X.atoms.find(a => bmArom(m.X, a.id)) || m.X.atoms[0], m.h.c] });
      const JOIN_WORD = { acyl: 'acyl substitution', sulfonamide: 'sulfonamide formation', sn2: 'Sₙ2 substitution', grignard: 'carbanion addition to the carbonyl', aldol: 'aldol reaction', imine: 'imine formation', da: 'Diels–Alder cycloaddition', 'fc-acyl': 'Friedel–Crafts acylation', 'fc-alkyl': 'Friedel–Crafts alkylation',
        epoxNu: 'epoxide opening by a nucleophile (Sₙ2, less hindered end)',
        'grignard-co2': 'carboxylation with CO₂', 'grignard-ester': 'double addition to the ester', 'grignard-nitrile': 'addition to the nitrile', 'grignard-epoxide': 'epoxide opening', 'acetylide-sn2': 'acetylide alkylation', 'gilman-14': 'conjugate (1,4-) addition', 'gilman-couple': 'Corey–House coupling', 'gilman-acyl': 'cuprate acylation', 'gilman-none': 'no reaction (a cuprate does not add to a simple C=O)', deprotonate: 'deprotonation',
        suzuki: 'Suzuki–Miyaura coupling', stille: 'Stille coupling', negishi: 'Negishi coupling', kumada: 'Kumada coupling',
        sonogashira: 'Sonogashira coupling', buchwald: 'Buchwald–Hartwig amination', heck: 'Heck reaction', snar: 'nucleophilic aromatic substitution (Sₙ Ar)',
        charyl: 'direct C–H arylation',
        oxidS: 'sulfide → sulfoxide', oxidN: 'N-oxidation', oxidC: 'epoxidation',
        knoev: 'Knoevenagel condensation', michael: 'Michael addition (conjugate 1,4-)', acAlkyl: 'α-alkylation of a doubly activated carbon', wittig: 'Wittig reaction', hwe: 'Horner–Wadsworth–Emmons reaction',
        hetcyc: 'heterocycle-forming condensation',
        arhal: 'aromatic halogenation', thion: 'thionation (C=O → C=S)', diazotise: 'diazotisation',
        metalInsert: 'oxidative insertion (Grignard / organolithium formation)', stepPoly: 'step-growth polymerisation',
        alkylate: 'α-alkylation', aldolE: 'aldol reaction', michaelE: 'Michael addition (conjugate 1,4-)',
        claisen: 'Claisen condensation', stork: 'Stork enamine alkylation', storkMichael: 'Stork enamine Michael addition',
        enamine: 'enamine formation' };
      ms.forEach(m => {
        if (m.join === 'fc-none') return;
        const [hx0, hy0] = HANDLE0[m.join](m);
        const key = m.join + ':' + (hx0 ? hx0.id : '-') + ':' + (hy0 ? hy0.id : '-') + ':' + (m.swapped ? 1 : 0);
        if (seen.has(key) || out.length >= 4) return; seen.add(key);
        const [hx, hy] = HANDLE[m.join](m);
        const aAtom = (m.X === g ? hx : hy) || g.atoms[0];
        const NUC_IS_X = { grignard: true, aldol: true, da: true, 'fc-acyl': true, 'fc-alkyl': true, 'grignard-co2': true, 'grignard-ester': true, 'grignard-nitrile': true, 'grignard-epoxide': true, 'acetylide-sn2': true, 'gilman-14': true, 'gilman-couple': true, 'gilman-acyl': true, deprotonate: true };
        out.push({ kind: 'join', join: m.join, m, atoms: [aAtom.id], c: aAtom, score: m.score, label: JOIN_WORD[m.join] || m.join,
          subHandle: aAtom, rgHandle: (m.X === g ? hy : hx), subIsNuc: (m.X === g) === !!NUC_IS_X[m.join], subIsX: m.X === g,
          reason: 'the two molecules carry a matching pair of handles for ' + (JOIN_WORD[m.join] || m.join) + (ms.filter(x => x.join === m.join).length > 1 ? ' (at more than one site — the most reactive one is used)' : '') });
      });
      return out;
    },
    variants(g0, s, ruleId, env) {
      const r = bmBuild(s.m, env, g0); if (!r) return [];
      s.byprod = r.byprod; s.mechKind = r.mechKind; s.eqn = r.eqn;
      return [{ graph: r.graph, changed: r.changed, added: r.added, share: 1, label: r.label, why: r.why, stereoHow: r.stereoHow || null, mechKind: r.mechKind, eqn: r.eqn, _stages: r.stages, nameHint: r.nameHint }];
    },
    byproducts: (ruleId, env, site) => (site && site.byprod) ? site.byprod.map(b => ({ smiles: b.smiles || null, graph: b.graph || null, name: b.name })) : [],
    condNotes(ruleId, env) { return ['Molecule + molecule: the engine matched the functional groups of the two structures and ran the taught reaction between them. Conditions (base, acid catalyst, heat) are the ones that reaction normally needs — see the reagent line.']; },
    stages(g0, s, out) { return out._stages || null; },
  };

  /* =====================================================================
     PROTECTING GROUPS

     Half of a real medicinal-chemistry sequence is protection and
     deprotection, and none of it was here. The chemistry is simple — hang a
     group on, take it off again — so what matters is saying WHICH conditions
     take WHICH group off and, above all, which groups survive each other.
     That is the whole reason a chemist chooses Boc over Cbz: orthogonality.
     ===================================================================== */

  /* graft a fragment written as SMILES onto an atom; returns the new atom ids */
  function graft(g, atId, smiles, attachIdx) {
    const frag = C().parseSmiles(smiles);
    const off = Math.max(g.nextId || 0, ...g.atoms.map(a => a.id + 1));
    const map = new Map();
    const host = atom(g, atId);
    frag.atoms.forEach(a => { map.set(a.id, off + a.id); g.atoms.push(Object.assign({}, a, { id: off + a.id, x: host.x + 42 + (a.x - frag.atoms[0].x) * 0.9, y: host.y + 30 + (a.y - frag.atoms[0].y) * 0.9 })); });
    frag.bonds.forEach(b => g.bonds.push(Object.assign({}, b, { a: map.get(b.a), b: map.get(b.b) })));
    g.nextId = off + Math.max(...frag.atoms.map(a => a.id)) + 1;
    const anchor = map.get(frag.atoms[attachIdx === undefined ? 0 : attachIdx].id);
    bondOf2(g, atId, anchor, 1);
    try { C().layoutGraph(g); } catch (e) {}
    return [...map.values()];
  }
  const bondOf2 = (g, a, b, order) => { g.bonds.push({ a, b, order: order || 1 }); };
  /* cut atom `at` away from `keep` and delete the whole fragment that leaves with it */
  function cutFragment(g, keep, at) {
    g.bonds = g.bonds.filter(b => !((b.a === keep && b.b === at) || (b.b === keep && b.a === at)));
    const gone = new Set([at]); const st = [at];
    while (st.length) { const x = st.pop(); nb(g, x).forEach(n => { if (!gone.has(n.atom.id)) { gone.add(n.atom.id); st.push(n.atom.id); } }); }
    const frag = { atoms: g.atoms.filter(a => gone.has(a.id)).map(a => Object.assign({}, a)), bonds: g.bonds.filter(b => gone.has(b.a) && gone.has(b.b)).map(b => Object.assign({}, b)), nextId: g.nextId };
    g.atoms = g.atoms.filter(a => !gone.has(a.id));
    g.bonds = g.bonds.filter(b => !gone.has(b.a) && !gone.has(b.b));
    return frag;
  }

  /* a nitrogen that can be protected: an N–H that is not already an amide or sulfonamide */
  function freeAmineSites(g) {
    const out = [];
    g.atoms.forEach(a => {
      if (a.element !== 'N' || a.charge || hOn(g, a) < 1) return;
      const ns = nb(g, a.id);
      if (ns.some(n => n.bond.order > 1)) return;
      if (ns.some(n => n.atom.element === 'C' && nb(g, n.atom.id).some(m => m.atom.element === 'O' && m.bond.order === 2))) return;   // amide
      if (ns.some(n => n.atom.element === 'S')) return;                                                                              // sulfonamide
      const aryl = ns.some(n => n.atom.element === 'C' && isAromatic(g, n.atom.id));
      out.push({ kind: 'amine', atoms: [a.id], c: a, n: a, deg: ns.length, aryl,
        label: (ns.length === 1 ? 'primary' : ns.length === 2 ? 'secondary' : 'tertiary') + ' amine N' + a.id + (aryl ? ' (an aniline — less nucleophilic)' : ''),
        score: (ns.length === 1 ? 0.4 : ns.length === 2 ? 0.2 : -3) - (aryl ? 0.5 : 0) });
    });
    return out;
  }
  /* a nitro group, and a C=N that a hydride can reduce */
  function nitroSites(g) {
    const out = [];
    g.atoms.forEach(a => {
      if (a.element !== 'N') return;
      /* `deg` counts CARBON neighbours, which is zero for a nitro oxygen — count bonds instead */
      const os = nb(g, a.id).filter(n => n.atom.element === 'O' && nb(g, n.atom.id).length === 1);
      const cAt = nb(g, a.id).find(n => n.atom.element === 'C');
      if (os.length === 2 && cAt) out.push({ kind: 'nitro', atoms: [a.id].concat(os.map(o => o.atom.id)), c: a, n: a, os: os.map(o => o.atom), cAt: cAt.atom,
        aryl: isAromatic(g, cAt.atom.id), score: isAromatic(g, cAt.atom.id) ? 0.3 : 0,
        label: 'nitro group on ' + (isAromatic(g, cAt.atom.id) ? 'the ring' : 'C' + cAt.atom.id) + ' (N' + a.id + ')' });
    });
    return out;
  }
  /* what is ALREADY protected, and what takes it off */
  const PG_KEY = {
    boc: { off: 'boc_off', word: 'Boc carbamate', by: 'TFA in CH₂Cl₂, or 4 M HCl in dioxane' },
    cbz: { off: 'cbz_off', word: 'Cbz carbamate', by: 'H₂ over Pd–C (hydrogenolysis)' },
    fmoc: { off: 'fmoc_off', word: 'Fmoc carbamate', by: '20 % piperidine in DMF' },
    tbs: { off: 'tbs_off', word: 'silyl ether', by: 'TBAF, or HF·pyridine, or dilute AcOH' },
    acetal: { off: 'acetal_off', word: '1,3-dioxolane acetal', by: 'dilute aqueous acid' },
    bn: { off: 'bn_off', word: 'benzyl ether', by: 'H₂ over Pd–C (hydrogenolysis)' },
  };
  function protectedSites(g, which) {
    const out = [];
    g.atoms.forEach(a => {
      const ns = nb(g, a.id);
      /* a carbamate: N–C(=O)–O–R. Which R decides which group it is. */
      if (a.element === 'C') {
        const dO = ns.find(n => n.atom.element === 'O' && n.bond.order === 2 && nb(g, n.atom.id).length === 1);
        const oR = ns.find(n => n.atom.element === 'O' && n.bond.order === 1);
        const nN = ns.find(n => n.atom.element === 'N');
        if (dO && oR && nN) {
          const r = nb(g, oR.atom.id).find(n => n.atom.id !== a.id);
          if (r && r.atom.element === 'C') {
            const rn = nb(g, r.atom.id);
            /* tert-butyl: three methyls on the carbon that holds the ester oxygen. `deg` counts
               CARBON neighbours only, so it is 3 here, not 4 — the O is not counted. */
            const tBu = hOn(g, r.atom) === 0 && rn.filter(n => n.atom.element === 'C' && hOn(g, n.atom) === 3).length === 3;
            /* benzyl (Cbz): a CH2 on an aromatic carbon. Fluorenylmethyl (Fmoc): a CH2 on the
               sp3 C9 of a fluorene — not aromatic itself, but attached to a big fused system. */
            const benzylic = hOn(g, r.atom) >= 1 && rn.some(n => n.atom.element === 'C' && isAromatic(g, n.atom.id));
            let fused = false;
            if (hOn(g, r.atom) >= 1) {
              try {
                const sys = C().ringSystemsOf(g, C().ringAtomsOf(g));
                fused = rn.some(n => n.atom.element === 'C' && sys.some(sy => sy.length >= 12 && sy.indexOf(n.atom.id) >= 0));
              } catch (e) {}
            }
            const kind = tBu ? 'boc' : fused ? 'fmoc' : benzylic ? 'cbz' : null;
            if (kind) out.push({ kind: 'protected', pg: kind, atoms: [a.id, oR.atom.id, r.atom.id], c: a, n: nN.atom, cO: a, oR: oR.atom, r: r.atom,
              label: PG_KEY[kind].word + ' on N' + nN.atom.id, score: 0 });
          }
        }
        /* an acetal: one carbon carrying two single-bonded oxygens, each on carbon */
        const oo = ns.filter(n => n.atom.element === 'O' && n.bond.order === 1 && nb(g, n.atom.id).some(m => m.atom.element === 'C' && m.atom.id !== a.id));
        if (oo.length === 2 && !ns.some(n => n.atom.element === 'O' && n.bond.order === 2))
          out.push({ kind: 'protected', pg: 'acetal', atoms: [a.id].concat(oo.map(o => o.atom.id)), c: a, os: oo.map(o => o.atom), label: '1,3-dioxolane acetal at C' + a.id, score: 0 });
      }
      /* a silyl ether, and a benzyl ether */
      if (a.element === 'O' && ns.length === 2) {
        const si = ns.find(n => n.atom.element === 'Si');
        const cSide = ns.find(n => n.atom.element === 'C');
        if (si && cSide) out.push({ kind: 'protected', pg: 'tbs', atoms: [a.id, si.atom.id], c: a, o: a, si: si.atom, cSide: cSide.atom, label: 'silyl ether on O' + a.id, score: 0 });
        const bz = ns.filter(n => n.atom.element === 'C' && hOn(g, n.atom) >= 2 && nb(g, n.atom.id).some(m => m.atom.element === 'C' && isAromatic(g, m.atom.id)));
        const other = ns.find(n => n.atom.element === 'C' && bz.indexOf(n) < 0);
        if (bz.length === 1 && other && !si) out.push({ kind: 'protected', pg: 'bn', atoms: [a.id, bz[0].atom.id], c: a, o: a, r: bz[0].atom, cSide: other.atom, label: 'benzyl ether on O' + a.id, score: 0 });
      }
    });
    return which ? out.filter(x => x.pg === which) : out;
  }

  const PROT = {
    boc_on: { on: 'amine', frag: 'C(=O)OC(C)(C)C', word: 'Boc-protected amine', reagent: 'Boc₂O (di-tert-butyl dicarbonate)',
      why: 'The amine nitrogen attacks one carbonyl of Boc₂O and tert-butyl carbonate leaves, taking a proton and falling apart to CO₂ and tert-butanol — which is why the reaction is irreversible and needs only a mild base. The nitrogen is now a carbamate: still there, but no longer nucleophilic or basic, so it will sit quietly through an alkylation, an acylation or a coupling.',
      by: [{ smiles: 'OC(C)(C)C', name: 'tert-butanol (with CO₂)' }] },
    cbz_on: { on: 'amine', frag: 'C(=O)OCc1ccccc1', word: 'Cbz-protected amine', reagent: 'CbzCl (benzyl chloroformate)',
      why: 'Benzyl chloroformate acylates the nitrogen; a base takes up the HCl. Cbz does the same job as Boc but comes off by hydrogenolysis rather than acid, which is exactly why both exist: they are orthogonal.',
      by: [{ smiles: 'Cl', name: 'HCl (taken up by the base)' }] },
    fmoc_on: { on: 'amine', frag: 'C(=O)OCC1c2ccccc2-c2ccccc21', word: 'Fmoc-protected amine', reagent: 'Fmoc-Cl or Fmoc-OSu',
      why: 'The third member of the orthogonal set: Fmoc comes off with a mild secondary amine (piperidine) and is untouched by acid or hydrogen, which is what makes solid-phase peptide synthesis work.',
      by: [{ smiles: 'Cl', name: 'HCl (taken up by the base)' }] },
    tbs_on: { on: 'alcohol', frag: '[Si](C)(C)C(C)(C)C', word: 'TBS silyl ether', reagent: 'TBSCl / imidazole (or TBSOTf, 2,6-lutidine)',
      why: 'The alcohol oxygen attacks silicon and chloride leaves. A silyl ether is not nucleophilic, not acidic and not oxidisable, so the alcohol is out of the way — and the Si–O bond is broken again by fluoride, which nothing else in the molecule cares about. Bulk sets the pace: a primary alcohol silylates much faster than a secondary, and a tertiary one may not react at all.',
      by: [{ smiles: 'Cl', name: 'HCl (taken up by imidazole)' }] },
    thp_on: { on: 'alcohol', frag: 'C1CCCCO1', word: 'THP ether', reagent: 'dihydropyran, cat. TsOH',
      why: 'Acid protonates the enol ether of dihydropyran, the alcohol adds to the oxocarbenium, and the result is a mixed acetal — stable to base, to hydride and to organometallics, and removed by dilute aqueous acid. Cheap, but it makes a new stereocentre, so the NMR gets messy.', by: [] },
    ac_on: { on: 'both', frag: 'C(C)=O', word: 'acetate ester', reagent: 'Ac₂O / pyridine (or AcCl)',
      why: 'Straight acylation of the alcohol. An acetate is the cheapest protection there is and comes off with K₂CO₃ in methanol or aqueous LiOH — but it is an ester, so it will not survive a hydride or a Grignard.',
      by: [{ smiles: 'CC(=O)O', name: 'acetic acid' }] },
    bn_on: { on: 'alcohol', frag: 'Cc1ccccc1', word: 'benzyl ether', reagent: 'BnBr, NaH (Williamson)',
      why: 'The alkoxide displaces bromide from benzyl bromide. A benzyl ether is about as robust as protection gets — acid, base, hydride and organometallics all leave it alone — and it comes off by hydrogenolysis over palladium, which is also its limitation: anything else in the molecule that hydrogenates will go too.',
      by: [{ smiles: 'Br', name: 'NaBr' }] },
    ms_on: { on: 'alcohol', frag: 'S(=O)(=O)C', word: 'mesylate', reagent: 'MsCl / Et₃N',
      why: 'Not a protection at all — an ACTIVATION. An alcohol is a hopeless leaving group; a mesylate is an excellent one, so this is the standard way to turn C–OH into something Sₙ2 will displace, with the configuration at that carbon untouched (no bond to it is broken here).',
      by: [{ smiles: 'Cl', name: 'Et₃N·HCl' }] },
    ts_on: { on: 'alcohol', frag: 'S(=O)(=O)c1ccc(C)cc1', word: 'tosylate', reagent: 'TsCl / pyridine',
      why: 'The same activation as a mesylate, with a crystalline tosyl group that is easier to handle and to see. The alcohol oxygen attacks sulfur and chloride leaves — again, no bond to the carbon is broken, so the stereocentre is retained until the nucleophile arrives and inverts it.',
      by: [{ smiles: 'Cl', name: 'pyridine·HCl' }] },
  };

  F.protect = {
    rules: Object.keys(PROT).concat(['acetal_on']), needs: 'protectable', structural: true,
    sites(g, ruleId) {
      if (ruleId === 'acetal_on') return carbonylSites(g).map(x =>
        (x.cls === 'ketone' || x.cls === 'aldehyde') ? Object.assign(x, { score: x.cls === 'aldehyde' ? 0.5 : 0 })
          : Object.assign(x, { inert: 'a ' + x.cls + ' — only a ketone or an aldehyde forms an acetal; an ester or amide has no electrophilic carbonyl left to add two alcohols to' }));
      const spec = PROT[ruleId]; if (!spec) return [];
      if (spec.on === 'amine') return freeAmineSites(g);
      const ohs = alcoholSites(g).map(x => x.phenol ? Object.assign(x, { score: (x.score || 0) - 0.3 }) : x)
        .map(x => x.acid ? Object.assign(x, { inert: 'a carboxylic acid, not an alcohol — this reagent is for an O–H on carbon' }) : x);
      /* an acylation is not fussy about which heteroatom, but it IS fussy about which is the
         better nucleophile: an amine goes long before a phenol, which is why 4-aminophenol and
         acetic anhydride give paracetamol and not the aryl acetate */
      if (spec.on === 'both') return freeAmineSites(g).map(x => Object.assign(x, { score: (x.score || 0) + 2.2 })).concat(ohs);
      return ohs;
    },
    variants(g0, s, ruleId) {
      if (ruleId === 'acetal_on') {
        /* the C=O becomes a 1,3-dioxolane: two oxygens on the old carbonyl carbon, joined by
           the glycol's two carbons */
        const g = clone(g0);
        const c = s.c.id;
        bondOf(g, c, s.o.id).order = 1;                       // the carbonyl O keeps its place as one acetal O
        const added = graft(g, s.o.id, 'CCO', 0);              // –O–CH₂CH₂–OH
        const tail = added[added.length - 1];                  // the new OH oxygen
        bondOf2(g, c, tail, 1);                                // close the five-membered ring
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([c, s.o.id].concat(added)), added, share: 1, label: '1,3-dioxolane (carbonyl protected)',
          why: 'Ethylene glycol and an acid catalyst, with the water removed as it forms (a Dean–Stark trap) to drive the equilibrium. Protonate the carbonyl, add one alcohol, lose water through the oxocarbenium, then close the ring with the second alcohol — a cyclic acetal, which is what makes it favourable: one molecule of diol, one ring, entropy on your side. ' +
            'The point of doing it: an acetal is completely inert to base, to hydride and to organometallics, so a ketone can be hidden while a Grignard is added somewhere else, then brought back with dilute aqueous acid. It is the classic answer to “this molecule has two carbonyls and I only want to touch one”.',
          stereoHow: 'The carbonyl carbon becomes sp³ but carries two identical-through-the-ring oxygens, so no new stereocentre is created.' }];
      }
      const spec = PROT[ruleId];
      const g = clone(g0);
      const at = (spec.on === 'amine' || s.kind === 'amine') ? s.n.id : s.o.id;
      const added = graft(g, at, spec.frag, 0);
      return [{ graph: g, changed: new Set([at].concat(added)), added, share: 1, label: spec.word,
        why: spec.why, stereoHow: (ruleId === 'ms_on' || ruleId === 'ts_on') ? 'No bond to the carbinol carbon is broken, so its configuration is untouched — the inversion comes later, when the nucleophile displaces the sulfonate.' : null }];
    },
    byproducts: (ruleId) => ruleId === 'acetal_on' ? [{ smiles: 'O', name: 'water — removed with a Dean–Stark trap to drive the equilibrium' }] : ((PROT[ruleId] || {}).by || []),
    condNotes(ruleId, env) {
      const n = [];
      if (ruleId === 'acetal_on') n.push('Ethylene glycol, cat. TsOH, benzene or toluene at reflux with a Dean–Stark trap. Every step is reversible: without removing the water there is no acetal.');
      if (ruleId === 'boc_on') n.push('Boc₂O with NaOH/dioxane or Et₃N/CH₂Cl₂; the by-products are CO₂ and tert-butanol, so nothing needs removing.');
      if (ruleId === 'tbs_on') n.push('TBSCl + imidazole in DMF for a primary or secondary alcohol; a hindered one needs TBSOTf with 2,6-lutidine.');
      if (ruleId === 'ms_on' || ruleId === 'ts_on') n.push('Keep it cold (0 °C) and use a hindered base: a warm mesylate with chloride around simply becomes the chloride.');
      n.push('Orthogonality is the point of choosing one group over another: Boc off with acid, Cbz off with H₂/Pd, Fmoc off with piperidine, a silyl ether off with fluoride. Each is blind to the other three.');
      return n;
    },
  };

  const DEPROT = {
    boc_off: { pg: 'boc', word: 'free amine (Boc removed)', reagent: 'TFA / CH₂Cl₂',
      why: 'Acid protonates the carbamate carbonyl and the tert-butyl group leaves as the tert-butyl cation — which is why this is so easy and so clean: the cation loses a proton to give isobutylene, the carbamic acid left behind falls apart to CO₂ and the free amine, and both by-products are gases. Nothing else in the molecule is touched by TFA at room temperature, unless it happens to be another acid-labile group.',
      by: [{ smiles: 'CC(C)=C', name: 'isobutylene (gas)' }, { smiles: 'O=C=O', name: 'CO₂ (gas)' }] },
    cbz_off: { pg: 'cbz', word: 'free amine (Cbz removed)', reagent: 'H₂, Pd–C',
      why: 'Hydrogenolysis: palladium cleaves the benzylic C–O bond, giving toluene and a carbamic acid that immediately loses CO₂. Neutral conditions throughout — the reason to pick Cbz when the molecule cannot see acid. Watch what else in the molecule hydrogenates.',
      by: [{ smiles: 'Cc1ccccc1', name: 'toluene' }, { smiles: 'O=C=O', name: 'CO₂' }] },
    fmoc_off: { pg: 'fmoc', word: 'free amine (Fmoc removed)', reagent: '20 % piperidine / DMF',
      why: 'A base removes the fluorenyl C9–H, and the aromatic dibenzofulvene that forms drags the carbamate off with it (E1cb); piperidine then traps the fulvene. Mild, fast, and completely blind to acid-labile groups, which is what makes Fmoc/tBu peptide synthesis possible.',
      by: [{ smiles: 'O=C=O', name: 'CO₂' }] },
    tbs_off: { pg: 'tbs', word: 'free alcohol (silyl removed)', reagent: 'TBAF / THF (or HF·pyridine)',
      why: 'Fluoride attacks silicon — the Si–F bond is the strongest single bond fluorine makes, which is what drives it — and the alkoxide picks up a proton on workup. Nothing else in an organic molecule competes for fluoride, so this is about as selective a deprotection as exists.',
      by: [] },
    acetal_off: { pg: 'acetal', word: 'ketone / aldehyde unmasked', reagent: 'dilute aqueous HCl (acetone/water)',
      why: 'Acetal hydrolysis, and every step is the reverse of how it was made: protonate an acetal oxygen, lose the alcohol to give an oxocarbenium, add water, and lose the second alcohol. Water in excess pushes it all the way back to the carbonyl. Base does nothing at all to an acetal — that is the whole reason it is a useful mask for a ketone.',
      by: [{ smiles: 'OCCO', name: 'ethylene glycol' }] },
    bn_off: { pg: 'bn', word: 'free alcohol (benzyl removed)', reagent: 'H₂, Pd–C',
      why: 'Hydrogenolysis of the benzylic C–O bond over palladium, giving toluene and the free alcohol. Neutral and reliable — but it will also reduce any alkene, alkyne, nitro group or azide in the molecule, so it is the wrong choice if one of those has to survive.',
      by: [{ smiles: 'Cc1ccccc1', name: 'toluene' }] },
  };

  F.deprotect = {
    rules: Object.keys(DEPROT), needs: 'protected', structural: true,
    sites(g, ruleId) {
      const spec = DEPROT[ruleId]; if (!spec) return [];
      return protectedSites(g, spec.pg).map(x => Object.assign(x, { score: 0 }));
    },
    variants(g0, s, ruleId) {
      const spec = DEPROT[ruleId];
      const g = clone(g0);
      const changed = new Set();
      if (spec.pg === 'acetal') {
        /* A cyclic acetal is a RING: cutting one C–O bond does not detach anything, because the
           other oxygen still holds the diol on. So break BOTH bonds to the carbon first, then
           delete whatever fell off, then give the carbon a fresh double-bonded oxygen. */
        const c = s.c.id;
        const oIds = s.os.map(o => o.id);
        oIds.forEach(o => { g.bonds = g.bonds.filter(b => !((b.a === c && b.b === o) || (b.b === c && b.a === o))); });
        /* everything now unreachable from the carbon goes */
        const keep = new Set([c]); const st = [c];
        while (st.length) { const x = st.pop(); nb(g, x).forEach(n => { if (!keep.has(n.atom.id)) { keep.add(n.atom.id); st.push(n.atom.id); } }); }
        g.atoms = g.atoms.filter(a => keep.has(a.id));
        g.bonds = g.bonds.filter(b => keep.has(b.a) && keep.has(b.b));
        const o2 = C().addAtomNear(g, 'O', atom(g, c));
        bondOf(g, c, o2.id).order = 2;
        changed.add(c); changed.add(o2.id);
      } else if (spec.pg === 'tbs') {
        cutFragment(g, s.o.id, s.si.id);
        changed.add(s.o.id);
      } else if (spec.pg === 'bn') {
        cutFragment(g, s.o.id, s.r.id);
        changed.add(s.o.id);
      } else {
        /* a carbamate: the whole C(=O)O–R unit leaves the nitrogen */
        cutFragment(g, s.n.id, s.cO.id);
        changed.add(s.n.id);
      }
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed, added: [], share: 1, label: spec.word, why: spec.why, stereoHow: null }];
    },
    byproducts: (ruleId) => (DEPROT[ruleId] || {}).by || [],
    condNotes(ruleId, env) {
      const n = [];
      if (ruleId === 'boc_off') n.push('25–50 % TFA in CH₂Cl₂ for 30 min at room temperature, or 4 M HCl in dioxane; the product is isolated as the salt and freebased if the next step needs the neutral amine.');
      if (ruleId === 'cbz_off' || ruleId === 'bn_off') n.push('1 atm H₂ over 10 % Pd–C in methanol. Anything else reducible in the molecule — an alkene, an alkyne, a nitro group, an aryl halide — goes as well.');
      if (ruleId === 'tbs_off') n.push('TBAF in THF, or HF·pyridine when the molecule cannot see a strong base; a primary TBS ether also comes off in dilute AcOH/water, which lets a primary be removed while a secondary stays on.');
      if (ruleId === 'acetal_off') n.push('Dilute aqueous acid with the water in excess — the same equilibrium that made the acetal, pushed the other way.');
      return n;
    },
  };

  /* nitro → amine: the standard way to get an aniline, since the ring is nitrated first */
  F.nitro = {
    rules: ['nitrored'], needs: 'nitro', structural: true,
    sites(g) { return nitroSites(g); },
    variants(g0, s) {
      const g = clone(g0);
      const n = atom(g, s.n.id);
      s.os.forEach(o => { try { C().removeAtom(g, o.id); } catch (e) {} });
      n.charge = 0;
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.n.id]), added: [], share: 1, label: s.aryl ? 'aniline (nitro reduced)' : 'primary amine (nitro reduced)',
        why: 'Six electrons and six protons take the nitro group all the way down to the amine, through the nitroso and the hydroxylamine — none of which is normally seen, because each is reduced faster than the one before. ' +
          'Three reagents do it and the choice is about what else is in the molecule: H₂ over Pd–C is cleanest but also reduces alkenes, alkynes and benzyl groups; Fe or Zn in aqueous HCl is old, cheap and leaves those alone; SnCl₂ is the mild choice that tolerates almost everything. ' +
          (s.aryl ? 'This is how nearly every aniline is made — nitrate the ring first (the nitro group is a meta director), then reduce, because aniline itself cannot be nitrated cleanly.' : ''),
        stereoHow: null }];
    },
    byproducts: () => [{ smiles: 'O', name: 'water (2 equivalents)' }],
    condNotes: () => ['H₂ (1–4 atm) over Pd–C or PtO₂; or Fe / Zn / Sn in aqueous HCl; or SnCl₂·2H₂O in ethanol. Iron in acid is the one to use when an aryl halide or an alkene has to survive.'],
  };

  /* =====================================================================
     CYCLISATION — the same joins, but within one molecule

     Everything in the molecule + molecule engine needs two molecules. A great
     deal of real chemistry does not: an amino ester closes to a lactam on
     heating, a hydroxy acid to a lactone, an amino alkyl halide to a
     pyrrolidine, a phenol with a chloride four atoms away to a benzofuran.
     In the external reaction set, "ring formation" was the largest family the
     engine simply could not express.

     The chemistry is the same nucleophilic substitution; what is new is that
     both halves are already tethered, so the only extra question is RING SIZE.
     Five wins, then six; three and seven are possible and slower; four is
     bad (both strain and a poor approach angle); eight and above needs high
     dilution or it oligomerises instead, and the app says so.
     ===================================================================== */
  const RING_PREF = { 3: -1.2, 4: -2.6, 5: 1.2, 6: 0.9, 7: -0.6, 8: -1.8, 9: -2.2, 10: -2.4 };
  const ringPref = n => (RING_PREF[n] !== undefined ? RING_PREF[n] : -2.8);
  function bondPath(g, a, b) {
    if (a === b) return 0;
    const seen = new Map([[a, 0]]); const q = [a];
    while (q.length) {
      const x = q.shift();
      for (const n of nb(g, x)) {
        if (seen.has(n.atom.id)) continue;
        seen.set(n.atom.id, seen.get(x) + 1);
        if (n.atom.id === b) return seen.get(n.atom.id);
        q.push(n.atom.id);
      }
    }
    return 99;
  }

  /* =====================================================================
     HETEROCYCLE-FORMING CONDENSATIONS — two molecules become a ring

     This is the largest single family of reactions in the medicinal-chemistry
     literature that a rule-per-reagent engine cannot express, because the
     product is not the substrate with a group changed: it is a new aromatic
     ring built from atoms contributed by both partners.

     Every one of them is the same reaction twice. One partner brings TWO
     nucleophilic atoms a fixed distance apart (an o-diamine, a hydrazine, an
     amidine, a thiourea); the other brings ONE or TWO electrophilic carbons a
     fixed distance apart (an aldehyde, an acid, a 1,3-diketone, an α-halo
     ketone). Each nucleophile adds to a carbon, each carbon throws out its
     leaving group as water or an alcohol, and the ring that closes is
     aromatic — which is what makes the whole sequence irreversible and is why
     these reactions work as well as they do in a flask with nothing but heat
     and a little acid.

     So the app does not need a rule for benzimidazoles and another for
     pyrazoles. It needs to recognise the two halves, add up the distances, and
     check that the ring comes out five- or six-membered.
     ===================================================================== */

  /* the two nucleophilic atoms, and how many atoms sit between them */
  function biNucleophiles(g) {
    const out = [];
    const freeN = a => a.element === 'N' && !a.charge && hOn(g, a) >= 1 &&
      !nb(g, a.id).some(n => n.atom.element === 'C' && bmDblO(g, n.atom.id) && nb(g, n.atom.id).filter(m => m.atom.element === 'N').length === 1);
    const oh = a => (a.element === 'O' || a.element === 'S') && !a.charge && hOn(g, a) === 1;
    const twoApart = (x, y) => {                   // x–A–B–y : exactly two atoms between
      for (const p of nb(g, x.id)) for (const q of nb(g, p.atom.id))
        if (q.atom.id !== x.id && nb(g, q.atom.id).some(r => r.atom.id === y.id) &&
            p.atom.element === 'C' && q.atom.element === 'C') return [p.atom, q.atom];
      return null;
    };
    g.atoms.forEach(a => {
      /* N–N (a hydrazine, an aryl hydrazine, a semicarbazide) and N–O (a hydroxylamine) */
      if (freeN(a)) nb(g, a.id).forEach(n => {
        if (n.bond.order !== 1) return;
        if (n.atom.element === 'N' && hOn(g, n.atom) >= 1 && a.id < n.atom.id)
          out.push({ n1: a, n2: n.atom, span: 0, kind: 'hydrazine', word: 'hydrazine (N–N)' });
        if (n.atom.element === 'O' && hOn(g, n.atom) === 1)
          out.push({ n1: n.atom, n2: a, span: 0, kind: 'hydroxylamine', word: 'hydroxylamine (H₂N–OH)' });
      });
      /* N–C(=N or =O or =S)–N : an amidine, a guanidine, a urea or a thiourea.
         For a thiourea the sulfur is the better nucleophile of the two, and that is
         exactly what decides the regiochemistry of a Hantzsch thiazole synthesis. */
      if (a.element === 'C') {
        const dbl = nb(g, a.id).find(n => n.bond.order === 2 && ['N', 'O', 'S'].includes(n.atom.element));
        const sng = nb(g, a.id).filter(n => n.bond.order === 1 && n.atom.element === 'N' && hOn(g, n.atom) >= 1);
        if (dbl && sng.length) {
          if (dbl.atom.element === 'N' && hOn(g, dbl.atom) >= 1)
            out.push({ n1: dbl.atom, n2: sng[0].atom, span: 1, mid: a, kind: 'amidine', word: 'amidine / guanidine (N–C=N)' });
          else if (dbl.atom.element === 'N')
            out.push({ n1: sng[0].atom, n2: sng.length > 1 ? sng[1].atom : dbl.atom, span: 1, mid: a, kind: 'amidine', word: 'amidine / guanidine (N–C=N)' });
          if (dbl.atom.element === 'S' && sng.length >= 1)
            out.push({ n1: dbl.atom, n2: sng[0].atom, span: 1, mid: a, kind: 'thiourea', word: 'thiourea / thioamide (N–C=S)' });
          if (dbl.atom.element === 'O' && sng.length >= 2)
            out.push({ n1: sng[0].atom, n2: sng[1].atom, span: 1, mid: a, kind: 'urea', word: 'urea (N–CO–N)' });
        }
      }
      /* 1,2 on a chain or a ring: an o-diamine, a 2-aminophenol, a 2-aminothiophenol */
      if (freeN(a)) g.atoms.forEach(bq => {
        if (bq.id <= a.id && freeN(bq)) return;
        if (!(freeN(bq) || oh(bq)) || bq.id === a.id) return;
        const path = twoApart(a, bq);
        if (!path) return;
        const kind = bq.element === 'N' ? 'diamine' : bq.element === 'O' ? 'aminophenol' : 'aminothiol';
        if (kind === 'diamine' && bq.id < a.id) return;                 // once per pair
        out.push({ n1: a, n2: bq, span: 2, path, kind,
          word: kind === 'diamine' ? 'a 1,2-diamine' : kind === 'aminophenol' ? 'a 2-aminophenol / amino alcohol' : 'a 2-aminothiophenol / amino thiol',
          arom: isAromatic(g, path[0].id) && isAromatic(g, path[1].id) });
      });
      /* one primary amine, used TWICE — this is how a Paal–Knorr pyrrole is made */
      if (freeN(a) && hOn(g, a) >= 2) out.push({ n1: a, n2: a, span: -1, kind: 'amine', word: 'a primary amine (both bonds come from the same nitrogen)' });
      /* ammonia, water and hydrogen sulfide do the same job for a furan or a thiophene */
      if ((a.element === 'O' || a.element === 'S') && hOn(g, a) === 2 && nb(g, a.id).length === 0)
        out.push({ n1: a, n2: a, span: -1, kind: a.element === 'O' ? 'water' : 'sulfide', word: a.element === 'O' ? 'water' : 'hydrogen sulfide' });
    });
    return out;
  }

  /* the electrophilic carbon or carbons, with what each has to throw out */
  function biElectrophiles(g) {
    const out = [];
    const acylLike = a => {
      if (a.element !== 'C') return null;
      const dO = nb(g, a.id).find(n => n.atom.element === 'O' && n.bond.order === 2);
      if (!dO) return null;
      const oh = nb(g, a.id).find(n => n.atom.element === 'O' && n.bond.order === 1 && hOn(g, n.atom) === 1);
      const or = nb(g, a.id).find(n => n.atom.element === 'O' && n.bond.order === 1 && hOn(g, n.atom) === 0);
      const cl = nb(g, a.id).find(n => n.atom.element === 'Cl');
      const h = hOn(g, a) >= 1;
      return { c: a, o: dO.atom, lg: oh ? oh.atom : or ? or.atom : cl ? cl.atom : null,
        kind: oh ? 'acid' : or ? 'ester' : cl ? 'acyl chloride' : h ? 'aldehyde' : 'ketone' };
    };
    g.atoms.forEach(a => {
      const ac = acylLike(a);
      /* ONE carbon that takes both nucleophiles: an aldehyde, an acid, an ester, a nitrile */
      if (ac) out.push({ c1: a, c2: a, span: -1, e1: ac, e2: ac, kind: 'mono-' + ac.kind, word: 'a single ' + ac.kind + ' carbon (it takes both bonds)' });
      if (a.element === 'C' && nb(g, a.id).some(n => n.atom.element === 'N' && n.bond.order === 3))
        out.push({ c1: a, c2: a, span: -1, nitrile: nb(g, a.id).find(n => n.bond.order === 3).atom, kind: 'mono-nitrile', word: 'a nitrile carbon' });
      /* an ORTHOESTER: HC(OEt)₃ — three leaving groups on one carbon, and the classic way to
         put a bare C–H at the 2-position of a benzimidazole or a purine */
      if (a.element === 'C' && bmSp3(g, a)) {
        const ors = nb(g, a.id).filter(n => n.atom.element === 'O' && n.bond.order === 1 && hOn(g, n.atom) === 0);
        if (ors.length >= 3) out.push({ c1: a, c2: a, span: -1, ors: ors.map(x => x.atom), kind: 'orthoester', word: 'an orthoester (three OR groups on one carbon)' });
      }
      /* TWO carbonyl carbons, 1,2- / 1,3- / 1,4- to each other */
      if (!ac) return;
      g.atoms.forEach(bq => {
        if (bq.id <= a.id) return;
        const ac2 = acylLike(bq);
        if (!ac2) return;
        const d = bondPath(g, a.id, bq.id);
        if (d < 1 || d > 3) return;
        /* the atoms between them must be carbons, and not part of an aromatic ring */
        if (isAromatic(g, a.id) || isAromatic(g, bq.id)) return;
        out.push({ c1: a, c2: bq, span: d - 1, e1: ac, e2: ac2, kind: (d + 1) + '-dicarbonyl',
          word: 'a 1,' + (d + 1) + '-dicarbonyl' });
      });
      /* an α-HALO KETONE: the carbonyl is one electrophile and the C–X carbon the other */
      nb(g, a.id).forEach(n => {
        if (n.atom.element !== 'C' || n.bond.order !== 1 || !bmSp3(g, n.atom)) return;
        const x = nb(g, n.atom.id).find(m => ['Cl', 'Br', 'I'].includes(m.atom.element));
        if (!x) return;
        out.push({ c1: n.atom, c2: a, span: 0, e1: { c: n.atom, lg: x.atom, kind: 'halide' }, e2: ac,
          kind: 'alpha-halo ketone', word: 'an α-halo ketone (C–X and C=O next to each other)' });
      });
    });
    return out;
  }

  /* Kekulé the new ring if it can be aromatic: find the largest set of ring bonds whose two
     ends both still have room for a π bond, no two sharing an atom. Six bonds at most, so a
     plain search is instant — and if nothing can be matched the ring is simply left saturated,
     which is the right answer for an imidazoline or an oxazolidinone. */
  /* the atoms strictly between a and b along the shortest path (used to walk a 1,4-dicarbonyl) */
  function bondRingPath(g, a, b) {
    const prev = new Map([[a, null]]); const q = [a];
    while (q.length) {
      const x = q.shift();
      for (const n of nb(g, x)) {
        if (prev.has(n.atom.id)) continue;
        prev.set(n.atom.id, x);
        if (n.atom.id === b) { const out = []; let cur = x; while (cur !== null && cur !== a) { out.unshift(cur); cur = prev.get(cur); } return out; }
        q.push(n.atom.id);
      }
    }
    return [];
  }

  function aromatiseRing(g, ring) {
    const hasPi = id => nb(g, id).some(n => n.bond.order >= 2);
    const room = id => {
      const a = atom(g, id);
      if (!a || a.charge || hasPi(id)) return false;
      const sig = nb(g, id).length;
      if (a.element === 'C') return sig <= 3;
      if (a.element === 'N') return sig <= 2;
      return false;                                  // O and S give the ring their lone pair instead
    };
    const bonds = [];
    for (let i = 0; i < ring.length; i++) {
      const p = ring[i], q = ring[(i + 1) % ring.length];
      if (room(p) && room(q) && bondOf(g, p, q)) bonds.push([p, q]);
    }
    /* a ring nitrogen next to a ring C=O is the AMIDE nitrogen: it keeps its hydrogen, and the
       double bond goes elsewhere. That is what makes a pyrimidin-4(3H)-one the 3H tautomer and
       not the 1H one, and the same for a pyrazol-5-one. */
    const amideN = id => atom(g, id).element === 'N' && nb(g, id).some(n => n.atom.element === 'C' &&
      nb(g, n.atom.id).some(m => m.atom.element === 'O' && m.bond.order === 2));
    const wOf = ([p, q]) => {
      const e = [atom(g, p).element, atom(g, q).element].sort().join('');
      if (amideN(p) || amideN(q)) return 0.35;
      if (e === 'CC' || e === 'CN') return 1.0;
      if (e === 'CS') return 0.5;
      if (e === 'NN') return 0.15;                 // an N=N in the ring means it is NOT aromatic
      if (e === 'NO' || e === 'OO') return 0.05;
      return 0.3;
    };
    const candSet = ring.filter(room);
    let best = [], bestScore = -1e9;
    const walk = (i, used, acc, w) => {
      /* A maximum matching is not automatically the right one. What makes the ring aromatic is
         that every ring CARBON is in a π bond; a nitrogen left over simply carries the N–H (and
         an oxygen or sulfur was never a candidate). So a carbon that could have taken a double
         bond and did not is a real penalty — it is the difference between a pyrimidin-4(3H)-one
         and its non-aromatic 4(5H)- tautomer. */
      let pen = 0;
      candSet.forEach(id => { if (!used.has(id) && atom(g, id).element === 'C') pen++; });
      const sc = acc.length * 10 + w - 4 * pen;
      if (sc > bestScore) { bestScore = sc; best = acc.slice(); }
      for (let k = i; k < bonds.length; k++) {
        const [p, q] = bonds[k];
        if (used.has(p) || used.has(q)) continue;
        used.add(p); used.add(q); acc.push(bonds[k]);
        walk(k + 1, used, acc, w + wOf(bonds[k]));
        acc.pop(); used.delete(p); used.delete(q);
      }
    };
    walk(0, new Set(), [], 0);
    best.forEach(([p, q]) => { const bd = bondOf(g, p, q); if (bd) bd.order = 2; });
    return best.length;
  }

  /* what the finished ring is called, from its size and what is in it */
  function ringNameOf(g, ring, fused) {
    const els = ring.map(id => atom(g, id).element);
    const n = els.filter(e => e === 'N').length, o = els.filter(e => e === 'O').length, sx = els.filter(e => e === 'S').length;
    const five = ring.length === 5;
    if (five && n === 2 && !o && !sx) return fused ? 'benzimidazole' : 'imidazole or pyrazole';
    if (five && n === 1 && o) return fused ? 'benzoxazole' : 'oxazole or isoxazole';
    if (five && n === 1 && sx) return fused ? 'benzothiazole' : 'thiazole';
    if (five && n === 1) return 'pyrrole';
    if (five && o) return 'furan';
    if (five && sx) return 'thiophene';
    if (five && n === 3) return 'triazole';
    if (ring.length === 6 && n === 2) return fused ? 'quinoxaline / quinazoline' : 'pyrimidine or pyrazine';
    if (ring.length === 6 && n === 1) return 'pyridine';
    return ring.length + '-membered ring';
  }

  /* Draw an ANTI (opposite-face) relationship: the new group at c1 on a wedge and the one at
     c2 on a dash. This is what epoxide opening, halogen addition and every backside attack on
     a ring actually gives, and it is the difference between a trans-diol and a cis one. */
  function markAnti(g, c1, n1, c2, n2) {
    const b1 = bondOf(g, c1, n1), b2 = bondOf(g, c2, n2);
    if (b1) { b1.stereo = 'wedge'; b1.narrow = c1; }
    if (b2) { b2.stereo = 'dash'; b2.narrow = c2; }
  }
  /* the connected pieces of a graph, as id sets */
  function componentsOf(g) {
    const seen = new Set(), out = [];
    g.atoms.forEach(a => {
      if (seen.has(a.id)) return;
      const set = new Set([a.id]); const st = [a.id]; seen.add(a.id);
      while (st.length) { const x = st.pop(); nb(g, x).forEach(n => { if (!seen.has(n.atom.id)) { seen.add(n.atom.id); set.add(n.atom.id); st.push(n.atom.id); } }); }
      out.push(set);
    });
    return out;
  }

  /* =====================================================================
     THE COURSE FAMILIES

     Everything below was added because the Organic II lecture slides teach it
     and the engine could not express it. Each one is written the way the
     lectures present it: what happens, what controls it, and the mistake a
     student actually makes.
     ===================================================================== */

  /* --- epoxide opening: the one reaction where acid and base give DIFFERENT ends --- */
  function epoxideSites(g) {
    const out = [];
    g.atoms.forEach(o => {
      if (o.element !== 'O' || o.charge || hOn(g, o) !== 0) return;
      const ns = nb(g, o.id).filter(n => n.atom.element === 'C');
      if (ns.length !== 2) return;
      const a = ns[0].atom, b = ns[1].atom;
      if (!nb(g, a.id).some(n => n.atom.id === b.id)) return;          // the three-ring
      out.push({ o, a, b, subA: deg(g, a), subB: deg(g, b),
        hindA: hindrance(g, a.id), hindB: hindrance(g, b.id) });
    });
    return out;
  }
  /* which carbon the nucleophile hits. Under ACID the ring oxygen is protonated first, the
     C–O bond stretches towards the carbon that can carry positive charge, and the nucleophile
     arrives at the MORE substituted carbon (Markovnikov-like). Under BASE or with a strong
     nucleophile there is no protonation and no charge to stabilise, so it is a plain Sₙ2 and
     the nucleophile takes the LESS hindered carbon. That single contrast is the most-tested
     idea in the epoxide chapter. */
  F.epoxOpen = {
    rules: ['epox_h2o', 'epox_hx', 'epox_roh'], needs: 'epoxide',
    sites(g, ruleId, env) {
      return epoxideSites(g).map(s => Object.assign(s, { kind: 'epoxide',
        label: 'epoxide (C' + s.a.id + '–O–C' + s.b.id + ')',
        score: 1, reason: 'the epoxide — a three-membered ring carries about 27 kcal/mol of strain, and that strain is what makes an ether, normally the most inert group in organic chemistry, open up under conditions no ordinary ether would notice' }));
    },
    variants(g0, s, ruleId, env) {
      const acidic = (env.pHVal === undefined || env.pHVal < 8);
      /* Which nucleophile? The rule that got here may be HBr, HCl or HI (through the addition
         family's fallback), water, or an alcohol — and the halide, not water, is what attacks. */
      const XOF = { hbr: 'Br', hcl: 'Cl', hi: 'I', epox_hx: 'Br' };
      const xEl = XOF[ruleId] || env.xEl || null;
      const nucFor = xEl ? { word: 'the halide', prod: 'halohydrin' }
        : (ruleId === 'epox_roh') ? { word: 'the alcohol', prod: 'β-alkoxy alcohol' }
        : { word: 'water', prod: 'anti-1,2-diol' };
      const out = [];
      /* both ends, scored: acid prefers the more substituted carbon, base the less hindered */
      [[s.a, s.b, s.subA, s.subB], [s.b, s.a, s.subB, s.subA]].forEach(([hit, other, subH, subO]) => {
        const g = clone(g0);
        const bd = bondOf(g, hit.id, s.o.id);
        if (!bd) return;
        g.bonds = g.bonds.filter(x => x !== bd);
        const oAt = atom(g, s.o.id);
        let added = [];
        if (xEl) {
          const x = C().addAtomNear(g, xEl, atom(g, hit.id)); added.push(x.id);
        } else if (ruleId === 'epox_roh') {
          const o2 = C().addAtomNear(g, 'O', atom(g, hit.id));
          const me = C().addAtomNear(g, 'C', o2); added.push(o2.id, me.id);
        } else {
          const o2 = C().addAtomNear(g, 'O', atom(g, hit.id)); added.push(o2.id);
        }
        /* anti: the nucleophile arrives on the face opposite the departing oxygen */
        try { markAnti(g, hit.id, added[0], other.id, s.o.id); } catch (e) {}
        try { C().layoutGraph(g); } catch (e) {}
        const score = acidic ? subH : -subH - (hit === s.a ? s.hindA : s.hindB);
        out.push({ graph: g, changed: new Set([hit.id, s.o.id].concat(added)), added, _score: score,
          label: nucFor.prod + ' — opened at the ' + (subH >= 2 ? 'MORE' : 'less') + ' substituted carbon (C' + hit.id + ')',
          why: (acidic
            ? 'Acid first, then the nucleophile. The ring oxygen is protonated, and that turns a poor leaving group into a good one; as the C–O bond begins to break, positive charge builds on carbon, so the bond that stretches is the one at the carbon that can BEST carry it — the more substituted one. The nucleophile then arrives there. The transition state is somewhere between Sₙ1 and Sₙ2: there is no free carbocation (the nucleophile is already arriving), but the regiochemistry is the one a carbocation would give. '
            : 'No acid, so no protonation: this is a plain Sₙ2 on a strained ring. The nucleophile attacks the LESS hindered carbon, from the back side, and the ring oxygen leaves as an alkoxide which the workup protonates. ') +
            'Either way the attack is from the face opposite the oxygen, so the two new groups end up ANTI — which is why cyclohexene oxide gives the trans-diol and never the cis one. ' +
            (nucFor.prod === 'anti-1,2-diol' ? 'Note that this is the same trans-diol you would get from an alkene with a peracid and then water, and the OPPOSITE relative configuration to the syn-diol from OsO₄.' : ''),
          stereoHow: 'Anti opening, and inversion at the carbon attacked. On a ring the product is trans; from a chiral epoxide the configuration at the attacked carbon is inverted and the other centre is untouched.' });
      });
      const tot = out.reduce((t, o) => t + Math.pow(10, o._score), 0);
      out.forEach(o => { o.share = Math.pow(10, o._score) / tot; delete o._score; });
      out.sort((a, b) => b.share - a.share);
      return out;
    },
    byproducts: () => [],
    gate: () => null,
    condNotes(ruleId, env) {
      return [ruleId === 'epox_hx' ? 'The hydrogen halide, cold, in an inert solvent — the halide is the nucleophile and the acid is what activates the ring.'
        : ruleId === 'epox_roh' ? 'The alcohol as the solvent with a trace of strong acid; warm. With SODIUM alkoxide instead (no acid) the other end is attacked.'
        : 'Dilute aqueous acid. With aqueous hydroxide instead, the same diol forms but at the less hindered carbon.',
        'The regiochemistry flips with the conditions, and that is the whole point of the reaction: acid → the more substituted carbon, base or a strong nucleophile → the less hindered one.'];
    },
    stages: () => null,
  };

  /* --- ether cleavage: the only way to break an ether, and it needs the strongest acids --- */
  F.etherCleave = {
    rules: [], needs: 'ether',
    sites(g, ruleId, env) {
      const out = [];
      g.atoms.forEach(o => {
        if (o.element !== 'O' || o.charge || hOn(g, o) !== 0) return;
        const cs = nb(g, o.id).filter(n => n.atom.element === 'C' && n.bond.order === 1);
        if (cs.length !== 2) return;
        if (cs.some(n => bmDblO(g, n.atom.id))) return;                    // an ester, not an ether
        const a = cs[0].atom, b = cs[1].atom;
        if (nb(g, a.id).some(n => n.atom.id === b.id)) return;             // an epoxide: its own family
        const arA = isAromatic(g, a.id), arB = isAromatic(g, b.id);
        if (arA && arB) return;                                            // diaryl ether: HX will not touch it
        out.push({ kind: 'ether', o, a, b, arA, arB,
          label: 'ether (C' + a.id + '–O–C' + b.id + ')',
          score: 0,
          reason: 'the ether oxygen — protonated by the strong acid, which is the only way to make this bond breakable at all' + (arA || arB ? '; one side is aromatic, and an aryl C–O bond cannot do Sₙ2 or Sₙ1, so the ALKYL side is the one that goes and a phenol is left' : '') });
      });
      return out;
    },
    variants(g0, s, ruleId, env) {
      const xEl = ruleId === 'hcl' ? 'Cl' : ruleId === 'hi' ? 'I' : 'Br';
      const out = [];
      /* which side leaves: the aryl side never does; otherwise the one that can carry charge
         better (3° > 2° > 1°), and a symmetric ether gives two of the same halide */
      const ends = [[s.a, s.b, s.arA], [s.b, s.a, s.arB]].filter(([hit, , ar]) => !ar);
      ends.forEach(([hit, keep]) => {
        const g = clone(g0);
        const frag = cutFragment(g, s.o.id, hit.id);
        const x = C().addAtomNear(frag, xEl, frag.atoms.find(a => a.id === hit.id));
        try { C().layoutGraph(frag); } catch (e) {}
        /* the remaining piece keeps the oxygen as an O–H (an alcohol, or a phenol) */
        try { C().layoutGraph(g); } catch (e) {}
        const dg = deg(g0, hit);
        out.push({ graph: frag, changed: new Set([hit.id, x.id]), added: [x.id], _score: dg + (isAromatic(g0, hit.id) ? -9 : 0),
          label: 'alkyl ' + (xEl === 'I' ? 'iodide' : xEl === 'Br' ? 'bromide' : 'chloride') + ' + the alcohol' + (s.arA || s.arB ? ' (a PHENOL)' : ''),
          why: 'Ethers are the inert group of organic chemistry — they are what you run reactions IN — and this is the one thing that breaks them. It takes concentrated HBr or HI, hot. ' +
            'Two steps: the acid protonates the oxygen (making it a decent leaving group), then the halide attacks carbon and the alcohol leaves. Which mechanism depends on the carbon: a methyl or primary carbon goes by Sₙ2 (so the halide attacks the less hindered side), while a tertiary or benzylic carbon ionises first and goes by Sₙ1. ' +
            (s.arA || s.arB ? 'Here one side is an aromatic ring, and an aryl C–O bond will do NEITHER — an aryl cation is impossibly high in energy and backside attack on an sp² carbon is geometrically impossible. So the alkyl side always goes, and the product is the phenol plus the alkyl halide. This is exactly how a methyl ether is used to protect a phenol and then removed. ' : '') +
            'HCl barely works (chloride is too poor a nucleophile) and HI works best. With EXCESS acid and a symmetrical dialkyl ether you end up with two molecules of the alkyl halide, because the alcohol released is itself converted.',
          stereoHow: dg >= 3 ? 'Through a carbocation, so a stereocentre there is scrambled.' : dg === 2 ? 'Sₙ2 at a secondary carbon: inversion.' : null,
          _byover: [{ graph: bmTidy(g), name: (s.arA || s.arB) ? 'the phenol' : 'the alcohol (with excess HX it becomes a second alkyl halide)' }] });
      });
      if (!out.length) return [];
      const tot = out.reduce((t, o) => t + Math.pow(10, o._score), 0);
      out.forEach(o => { o.share = Math.pow(10, o._score) / tot; delete o._score; });
      out.sort((a, b) => b.share - a.share);
      return out;
    },
    byproducts: () => [],
    condNotes() { return ['Concentrated HBr or HI, at reflux, usually in excess. HCl is too weak a nucleophile to be useful.',
      'This is how a methyl or benzyl ether protecting group comes off a phenol — and BBr₃ in dichloromethane at −78 °C does the same job under conditions the rest of a molecule can survive.'] },
    stages: () => null,
  };

  /* --- C=O all the way to CH2: Clemmensen and Wolff-Kishner --- */
  F.deoxy = {
    rules: ['clemmensen', 'wolffkishner'], needs: 'carbonyl',
    sites(g, ruleId, env) {
      return carbonylSites(g).map(s => {
        if (s.cls !== 'ketone' && s.cls !== 'aldehyde') return Object.assign(s, { inert: s.cls + ' — only a KETONE or an ALDEHYDE is taken down to CH₂; an ester or an amide is not' });
        s.score = (s.cls === 'aldehyde' ? 0.3 : 0) - 0.3 * s.hind;
        s.reason = s.cls + ' — reduced all the way to CH₂: the oxygen is removed altogether, not turned into an alcohol';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      try { C().removeAtom(g, s.o.id); } catch (e) {}
      try { C().layoutGraph(g); } catch (e) {}
      const wk = ruleId === 'wolffkishner';
      return [{ graph: g, changed: new Set([s.c.id]), added: [], share: 1,
        label: 'the alkane (C=O → CH₂' + (wk ? ', Wolff–Kishner' : ', Clemmensen' ) + ')',
        why: wk
          ? 'The Wolff–Kishner reduction. Hydrazine condenses with the ketone to give the HYDRAZONE (C=N–NH₂) — ordinary imine chemistry — and then strong base at 175–200 °C does something unusual: it deprotonates the terminal NH₂, the resulting anion shifts a hydride onto carbon, and nitrogen leaves as N₂ gas. Losing a molecule of nitrogen is what makes the whole thing irreversible. ' +
            'Choose it when the molecule cannot see acid: it runs in strong BASE. The Clemmensen does the same job in strong acid, so between the two you can deoxygenate almost anything.'
          : 'The Clemmensen reduction. Zinc amalgam in concentrated hydrochloric acid takes a ketone or aldehyde straight to the CH₂ — the oxygen is removed, not reduced to an alcohol. The mechanism is still argued over; it happens on the zinc surface and does not appear to go through the alcohol, because the alcohol itself is not reduced under these conditions. ' +
            'Choose it when the molecule can survive strong ACID but not strong base. The Wolff–Kishner is its base-tolerant twin. ' +
            'The classic use for both: a Friedel–Crafts ACYLATION followed by this reduction is how you put a straight primary alkyl chain on a ring, because a Friedel–Crafts alkylation would rearrange it.',
        stereoHow: null,
        byprod: wk ? [{ smiles: 'N#N', name: 'N₂ — leaves as a gas, which is what drives it' }, { smiles: 'O', name: 'water' }] : [{ smiles: 'O', name: 'water' }] }];
    },
    byproducts: () => [],
    gate(g, ruleId, env) {
      if (ruleId === 'wolffkishner' && env.T < 400 && !env.TAuto) return null;
      return null;
    },
    condNotes(ruleId, env) {
      return ruleId === 'wolffkishner'
        ? ['H₂NNH₂ then KOH (or KOtBu) in diethylene glycol, 175–200 °C — the high temperature is not optional, because the nitrogen has to be expelled.',
           'Runs in strong base, so acid-sensitive groups survive. Esters do not (they would be saponified).']
        : ['Zn(Hg) amalgam in concentrated HCl, at reflux.',
           'Runs in strong acid, so base-sensitive groups survive; acetals, silyl ethers and tert-butyl esters do not.'];
    },
    stages(g0, s, out) {
      if (!out || !out.label || out.label.indexOf('Wolff') < 0) return null;
      try {
        const g = clone(g0);
        const n1 = C().addAtomNear(g, 'N', atom(g, s.c.id));
        const n2 = C().addAtomNear(g, 'N', atom(g, n1.id));
        try { C().removeAtom(g, s.o.id); } catch (e) {}
        const bd = bondOf(g, s.c.id, n1.id); if (bd) bd.order = 2;
        try { C().layoutGraph(g); } catch (e) {}
        return [{ label: 'the hydrazone (C=N–NH₂) — ordinary imine chemistry, before the base does its work', graph: g }];
      } catch (e) { return null; }
    },
  };

  /* --- ozonolysis: the alkene is cut in half --- */
  F.ozone = {
    rules: ['ozonolysis', 'ozonolysis_ox'], needs: 'alkene',
    sites(g, ruleId, env) {
      return alkeneSites(g).map(s => {
        if (s.aromatic) return Object.assign(s, { inert: 'aromatic C=C — ozone does not cleave a benzene ring under these conditions' });
        s.score = 1;
        s.reason = 'the C=C — ozone cleaves it completely, and each carbon keeps one oxygen';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const ox = ruleId === 'ozonolysis_ox';
      const g = clone(g0);
      const A = s.A.id, B = s.B.id;
      const bd = bondOf(g, A, B);
      if (!bd) return [];
      /* cut the bond, put an oxygen on each carbon */
      const hA = hOn(g0, s.A), hB = hOn(g0, s.B);
      g.bonds = g.bonds.filter(x => x !== bd);
      const oA = C().addAtomNear(g, 'O', atom(g, A)); bondOf(g, A, oA.id).order = 2;
      const oB = C().addAtomNear(g, 'O', atom(g, B)); bondOf(g, B, oB.id).order = 2;
      const added = [oA.id, oB.id];
      /* an OXIDATIVE workup takes any aldehyde on to the acid */
      if (ox) {
        [[A, hA], [B, hB]].forEach(([c, h]) => { if (h >= 1) { const o3 = C().addAtomNear(g, 'O', atom(g, c)); added.push(o3.id); } });
      }
      try { C().layoutGraph(g); } catch (e) {}
      /* if the alkene was in a ring the molecule stays in one piece (a dicarbonyl); if not,
         it falls into two, and both halves are real products */
      let pieces = [];
      try { pieces = componentsOf(g); } catch (e) { pieces = []; }
      const halves = [];
      if (pieces.length > 1) {
        /* the alkene was not in a ring, so the molecule really is in two pieces now: the bigger
           one is the product on the page and the smaller is listed as the other half */
        pieces.sort((p, q) => q.size - p.size);
        const keep = pieces[0];
        pieces.slice(1).forEach(set => {
          const frag = { atoms: g.atoms.filter(a => set.has(a.id)).map(a => Object.assign({}, a)),
            bonds: g.bonds.filter(b => set.has(b.a) && set.has(b.b)).map(b => Object.assign({}, b)), nextId: g.nextId };
          halves.push({ graph: bmTidy(frag), name: 'the other half of the alkene' });
        });
        g.atoms = g.atoms.filter(a => keep.has(a.id));
        g.bonds = g.bonds.filter(b => keep.has(b.a) && keep.has(b.b));
        try { C().layoutGraph(g); } catch (e) {}
      }
      const words = (hA >= 1 || hB >= 1) ? (ox ? 'carboxylic acid(s)' : 'aldehyde(s) / ketone(s)') : 'two ketones';
      return [{ graph: g, changed: new Set([A, B].concat(added)), added, share: 1, _byover: halves.length ? halves : undefined,
        label: 'the C=C cleaved — ' + words + (ox ? ' (oxidative workup)' : ' (reductive workup)'),
        why: 'Ozonolysis cuts the double bond in half and hands each carbon an oxygen. Ozone adds across the π bond to give a molozonide, which falls apart and recombines into the ozonide; then the workup decides what you isolate. ' +
          (ox ? 'With an OXIDISING workup (H₂O₂) any carbon that had a hydrogen on it ends up as a CARBOXYLIC ACID, and a carbon with no hydrogen can only be a ketone.'
              : 'With a REDUCING workup — dimethyl sulfide, or zinc in acetic acid — you stop at the aldehyde or ketone. That is the point of Me₂S: it mops up the extra oxygen (becoming DMSO) so the aldehyde survives.') + ' ' +
          'This makes ozonolysis a structure-determination tool as much as a synthesis: cut the alkene, identify the two halves, and you know where the double bond was. ' +
          (hA === 0 && hB === 0 ? 'Both carbons here are fully substituted, so both halves are ketones and the workup makes no difference.' : '') +
          (pieces.length > 1 ? '' : ' The alkene was in a ring, so nothing falls off — the molecule opens out into a single chain with a carbonyl at each end.'),
        stereoHow: null,
        byprod: ox ? [] : [{ smiles: 'CS(C)=O', name: 'DMSO (from the Me₂S that took the extra oxygen)' }] }];
    },
    byproducts: () => [],
    condNotes(ruleId) {
      return ruleId === 'ozonolysis_ox'
        ? ['O₃ in CH₂Cl₂ at −78 °C, then H₂O₂ — the oxidative workup, which takes any aldehyde on to the acid.']
        : ['O₃ in CH₂Cl₂ at −78 °C, then Me₂S (or Zn / AcOH) — the reductive workup, which stops at the aldehyde.',
           'Cold and dilute: ozonides are explosive and are never isolated.'];
    },
    stages: () => null,
  };

  /* --- hydration of an alkyne: the two regiochemistries, and the enol you never see --- */
  F.alkyneHyd = {
    rules: ['alkynehyd', 'alkynehydb'], needs: 'alkyne',
    sites(g, ruleId, env) {
      const out = [];
      g.atoms.forEach(a => {
        nb(g, a.id).forEach(n => {
          if (n.bond.order !== 3 || n.atom.element !== 'C' || a.id >= n.atom.id) return;
          const hA = hOn(g, a), hB = hOn(g, n.atom);
          out.push({ kind: 'alkyne', a, b: n.atom, terminal: hA === 1 || hB === 1,
            label: 'alkyne (C' + a.id + '≡C' + n.atom.id + ')', score: 1,
            reason: (hA === 1 || hB === 1 ? 'a terminal alkyne — the two reagents give OPPOSITE ends, which is the whole reason both exist' : 'an internal alkyne — hydration gives a mixture unless the two ends are electronically different') });
        });
      });
      return out;
    },
    variants(g0, s, ruleId, env) {
      const mark = ruleId === 'alkynehyd';        // Markovnikov: HgSO4 / H2SO4 / H2O
      const g = clone(g0);
      const hA = hOn(g0, s.a), hB = hOn(g0, s.b);
      /* Markovnikov puts the oxygen on the MORE substituted carbon; hydroboration on the less */
      let put = s.a, keep = s.b;
      if (mark) { if (hA > hB) { put = s.b; keep = s.a; } }
      else { if (hA < hB) { put = s.b; keep = s.a; } }
      const bd = bondOf(g, s.a.id, s.b.id); if (bd) bd.order = 1;
      const o = C().addAtomNear(g, 'O', atom(g, put.id));
      bondOf(g, put.id, o.id).order = 2;
      try { C().layoutGraph(g); } catch (e) {}
      const aldehyde = hOn(g, atom(g, put.id)) >= 1;
      return [{ graph: g, changed: new Set([s.a.id, s.b.id, o.id]), added: [o.id], share: 1,
        label: (aldehyde ? 'aldehyde' : 'ketone') + ' — ' + (mark ? 'Markovnikov' : 'anti-Markovnikov') + ' hydration',
        why: 'Water adds across the triple bond, and the ENOL that forms immediately tautomerises to the carbonyl — which is why you never isolate the enol and why the product looks like an addition of water plus a shift. ' +
          (mark ? 'Mercury(II) makes it Markovnikov: the oxygen ends up on the MORE substituted carbon, so a terminal alkyne gives a METHYL KETONE. The mercury is there because water alone is far too poor a nucleophile for an alkyne; Hg²⁺ coordinates the π bond first.'
                : 'Hydroboration makes it anti-Markovnikov: boron, being the bigger and less electronegative atom, goes to the LESS substituted carbon, and the peroxide oxidation replaces it with OH there. So a terminal alkyne gives an ALDEHYDE. Disiamylborane (or 9-BBN) is used rather than BH₃ because a bulky borane adds only once and will not touch the vinylborane again.') + ' ' +
          'These two reagents on the same terminal alkyne give different compounds — a methyl ketone or an aldehyde — and that contrast is the reason both are taught.',
        stereoHow: null,
        byprod: mark ? [] : [{ name: 'the boron residue (oxidised by H₂O₂ / NaOH)' }] }];
    },
    byproducts: () => [],
    condNotes(ruleId) {
      return ruleId === 'alkynehyd'
        ? ['H₂SO₄ with a catalytic amount of HgSO₄, in water.', 'A terminal alkyne gives the methyl ketone; an internal one gives a mixture unless the two ends differ.']
        : ['Disiamylborane (Sia₂BH) or 9-BBN, then H₂O₂ / NaOH.', 'A bulky borane is essential: BH₃ would add a second time.'];
    },
    stages: (g0, s, out) => null,
  };

  /* --- allylic and benzylic bromination: the radical that is worth making --- */
  F.allylHal = {
    rules: ['allylbr'], needs: 'alkene-or-arene', structural: true,
    sites(g, ruleId, env) {
      const out = [];
      g.atoms.forEach(a => {
        if (a.element !== 'C' || !bmSp3(g, a) || hOn(g, a) < 1) return;
        const ns = nb(g, a.id);
        const allylic = ns.some(n => n.atom.element === 'C' && nb(g, n.atom.id).some(m => m.bond.order === 2 && m.atom.element === 'C' && !isAromatic(g, m.atom.id)));
        const benzylic = ns.some(n => n.atom.element === 'C' && isAromatic(g, n.atom.id));
        if (!allylic && !benzylic) return;
        const nH = hOn(g, a), dg = deg(g, a);
        out.push({ kind: 'radical', c: a, allylic, benzylic,
          label: (benzylic ? 'benzylic' : 'allylic') + ' C–H (C' + a.id + ', ' + (dg <= 1 ? 'primary' : dg === 2 ? 'secondary' : 'tertiary') + ')',
          score: (dg >= 3 ? 1.0 : dg === 2 ? 0.6 : 0) + (benzylic && allylic ? 0.5 : 0),
          reason: (benzylic ? 'a benzylic C–H: the radical left behind is delocalised into the ring, which is worth about 13 kcal/mol' : 'an allylic C–H: the radical is delocalised over three carbons') + ', so this bond breaks far more easily than any ordinary C–H — about 88 kcal/mol against 98 for a plain secondary one' });
      });
      return out;
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      const br = C().addAtomNear(g, 'Br', atom(g, s.c.id));
      try { C().layoutGraph(g); } catch (e) {}
      const out = [{ graph: g, changed: new Set([s.c.id, br.id]), added: [br.id], _score: 1,
        label: (s.benzylic ? 'benzylic' : 'allylic') + ' bromide',
        why: 'A radical chain, and the reason it is selective is the STRENGTH of the bond it breaks. A bromine atom abstracts the hydrogen from the allylic or benzylic position because the radical left behind is resonance-stabilised; that bond is roughly 10 kcal/mol weaker than an ordinary C–H, and with a bromine atom (which is only just reactive enough) a 10 kcal/mol difference is enormous. ' +
          'NBS is used rather than bromine itself, and not for convenience: it keeps the concentration of Br₂ in the flask very LOW, by reacting with the HBr produced to release just enough Br₂ to continue. Low bromine concentration is what suppresses ordinary electrophilic ADDITION to the double bond — with a bottle of Br₂ you would get the dibromide instead, and that is the mistake this reaction is designed to avoid. ' +
          (s.allylic ? 'One warning: the allylic radical is delocalised, so bromine can end up at EITHER end of it. With an unsymmetrical alkene that means two products, and the double bond may appear to have moved.' : ''),
        stereoHow: 'The radical is planar, so a new stereocentre here comes out racemic.',
        byprod: [{ name: 'succinimide (the spent NBS)' }] }];
      /* the allylic radical is delocalised: bromine can arrive at the far end too */
      if (s.allylic) {
        try {
          const ns = nb(g0, s.c.id).filter(n => n.atom.element === 'C');
          for (const n of ns) {
            const far = nb(g0, n.atom.id).find(m => m.bond.order === 2 && m.atom.element === 'C' && m.atom.id !== s.c.id && !isAromatic(g0, m.atom.id));
            if (!far) continue;
            const g2 = clone(g0);
            bondOf(g2, n.atom.id, far.atom.id).order = 1;
            bondOf(g2, s.c.id, n.atom.id).order = 2;
            const br2 = C().addAtomNear(g2, 'Br', atom(g2, far.atom.id));
            try { C().layoutGraph(g2); } catch (e) {}
            out.push({ graph: g2, changed: new Set([s.c.id, n.atom.id, far.atom.id, br2.id]), added: [br2.id], _score: 0.4,
              label: 'the other end of the allylic radical — bromide at C' + far.atom.id + ', double bond moved',
              why: 'The allylic radical is one delocalised species with the spin spread over three carbons, so bromine can be delivered to either end. This is the product from the far end, and the double bond has apparently walked along the chain. With a symmetrical radical the two are the same compound; with an unsymmetrical one you get both, and the more substituted (more stable) alkene usually predominates.',
              stereoHow: null, byprod: [{ name: 'succinimide (the spent NBS)' }] });
            break;
          }
        } catch (e) {}
      }
      const tot = out.reduce((t, o) => t + Math.pow(10, o._score), 0);
      out.forEach(o => { o.share = Math.pow(10, o._score) / tot; delete o._score; });
      out.sort((a, b) => b.share - a.share);
      return out;
    },
    byproducts: () => [],
    condNotes() { return ['NBS with a radical initiator (AIBN or a peroxide) or light, in CCl₄ or benzene, at reflux.',
      'Br₂ with light does the same thing, but only if the bromine is kept very dilute — otherwise it adds across the double bond instead.'] },
    stages: () => null,
  };

  /* --- benzylic oxidation: the whole side chain becomes one carboxyl --- */
  F.benzylOx = {
    rules: [], needs: 'benzylic',
    sites(g, ruleId, env) {
      const out = [];
      g.atoms.forEach(a => {
        if (a.element !== 'C' || !bmSp3(g, a)) return;
        const ar = nb(g, a.id).find(n => n.atom.element === 'C' && isAromatic(g, n.atom.id));
        if (!ar) return;
        if (hOn(g, a) === 0) return;                    // a quaternary benzylic carbon survives
        out.push({ kind: 'benzylic', c: a, ring: ar.atom, score: 1,
          label: 'benzylic carbon (C' + a.id + ')',
          reason: 'the benzylic carbon — the ring survives and the whole side chain is cut back to a single carboxyl carbon' });
      });
      return out;
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      /* everything on the benzylic carbon except the ring goes; the carbon becomes COOH */
      nb(g0, s.c.id).forEach(n => {
        if (n.atom.id === s.ring.id) return;
        try { const f = cutFragment(g, s.c.id, n.atom.id); } catch (e) {}
      });
      const o1 = C().addAtomNear(g, 'O', atom(g, s.c.id)); bondOf(g, s.c.id, o1.id).order = 2;
      const o2 = C().addAtomNear(g, 'O', atom(g, s.c.id));
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.c.id, o1.id, o2.id]), added: [o1.id, o2.id], share: 1,
        label: 'benzoic acid (the side chain oxidised back to COOH)',
        why: 'Hot permanganate (or chromic acid) oxidises an alkyl group on a benzene ring all the way back to a CARBOXYLIC ACID, and it does not matter how long that chain was: propylbenzene, butylbenzene and cumene all give benzoic acid. What is left attached to the ring is a single carbon bearing the COOH. ' +
          'It starts at the benzylic C–H, because that is the weakest bond in the molecule, and then runs away with itself — every carbon beyond the first is lost as CO₂. ' +
          'The one requirement is a benzylic HYDROGEN. tert-Butylbenzene has none, and it survives this treatment completely — which is the standard exam question. ' +
          'The aromatic ring itself is untouched; benzene rings are remarkably resistant to oxidation, which is part of what "aromatic stability" means in practice.',
        stereoHow: null }];
    },
    byproducts: () => [{ name: 'CO₂ (every carbon beyond the first) and MnO₂' }],
    condNotes() { return ['KMnO₄, hot aqueous base, then acid to protonate the carboxylate; or Na₂Cr₂O₇ / H₂SO₄, hot.',
      'A benzylic hydrogen is required: tert-butylbenzene does not react at all.'] },
    stages: () => null,
  };

  /* --- aqueous acid, the workup that is also a reaction --- */
  F.aqAcid = {
    rules: ['aqacid'], needs: 'any', structural: true,
    sites(g, ruleId, env) {
      const out = [];
      /* what dilute aqueous acid actually does, in the order it does it */
      g.atoms.forEach(a => {
        const ns = nb(g, a.id);
        if (a.element === 'C') {
          /* an acetal or ketal: two single-bonded oxygens on one sp3 carbon */
          const os = ns.filter(n => n.atom.element === 'O' && n.bond.order === 1);
          if (bmSp3(g, a) && os.length === 2) out.push({ kind: 'acetal', c: a, os: os.map(x => x.atom), score: 3, label: 'acetal (C' + a.id + ')', reason: 'an acetal — it exists only because acid was removed, and putting acid and water back takes it straight off' });
          /* an imine / oxime / hydrazone */
          const dn = ns.find(n => n.atom.element === 'N' && n.bond.order === 2);
          if (dn) out.push({ kind: 'imine', c: a, n: dn.atom, score: 2.5, label: 'imine C=N (C' + a.id + ')', reason: 'an imine — made by losing water, so water and acid take it apart again' });
          /* a nitrile */
          const tn = ns.find(n => n.atom.element === 'N' && n.bond.order === 3);
          if (tn) out.push({ kind: 'nitrile', c: a, n: tn.atom, score: 1, label: 'nitrile (C' + a.id + '≡N)', reason: 'a nitrile — hot aqueous acid hydrolyses it, through the amide, to the carboxylic acid' });
        }
      });
      /* an epoxide */
      epoxideSites(g).forEach(e => out.push(Object.assign({ kind: 'epoxideAq', score: 2 }, e, { label: 'epoxide', reason: 'an epoxide — dilute acid opens it to the anti-diol at the more substituted carbon' })));
      return out;
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      if (s.kind === 'acetal') {
        /* both C–O bonds go, the carbonyl comes back */
        const frags = [];
        s.os.forEach(o => { try { const f = cutFragment(g, s.c.id, o.id); if (f.atoms.length) frags.push(f); } catch (e) {} });
        /* a cyclic acetal is one ring: cutting the first bond detaches nothing */
        const keep = g.atoms.some(a => a.id === s.c.id);
        const o = C().addAtomNear(g, 'O', atom(g, s.c.id)); bondOf(g, s.c.id, o.id).order = 2;
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.c.id, o.id]), added: [o.id], share: 1,
          label: 'the ketone / aldehyde, unmasked',
          why: 'An acetal is a protected carbonyl, and this is the deprotection. Every step of its formation was reversible — the only thing that made the acetal is REMOVING water; put water back, with acid, and the equilibrium runs the other way. ' +
            'That reversibility is the whole design: an acetal is stable to base, to hydride, to Grignard reagents and to oxidants, and comes off with nothing more than dilute acid and water. Le Châtelier does the work in both directions.',
          stereoHow: null, _byover: frags.filter(f => f.atoms.length).map(f => ({ graph: bmTidy(f), name: 'the diol (or two alcohols) released' })) }];
      }
      if (s.kind === 'imine') {
        const frag = cutFragment(g, s.c.id, s.n.id);
        const o = C().addAtomNear(g, 'O', atom(g, s.c.id)); bondOf(g, s.c.id, o.id).order = 2;
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.c.id, o.id]), added: [o.id], share: 1,
          label: 'the carbonyl back, and the free amine',
          why: 'Imine formation is an equilibrium, and this is the other direction. Water adds to the C=N, and the tetrahedral carbinolamine that results collapses the other way — expelling the amine instead of water. Mild aqueous acid is all it takes, which is why imines are hydrolysed on workup unless you take care. ' +
            'The same reaction, run deliberately, is how a ketone protected as an imine (or an enamine) is recovered, and how the Wolff–Kishner’s hydrazone would come apart if the base did not get there first.',
          stereoHow: null, _byover: [{ graph: bmTidy(frag), name: 'the amine, as its ammonium salt in the acid' }] }];
      }
      if (s.kind === 'nitrile') {
        const bd = bondOf(g, s.c.id, s.n.id);
        try { C().removeAtom(g, s.n.id); } catch (e) {}
        const o1 = C().addAtomNear(g, 'O', atom(g, s.c.id)); bondOf(g, s.c.id, o1.id).order = 2;
        const o2 = C().addAtomNear(g, 'O', atom(g, s.c.id));
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.c.id, o1.id, o2.id]), added: [o1.id, o2.id], share: 1,
          label: 'carboxylic acid (nitrile hydrolysed)',
          why: 'Hot aqueous acid hydrolyses a nitrile to the carboxylic acid, and it goes through the AMIDE on the way — which means you can stop there if you want the amide, by using milder conditions and less water. ' +
            'This is the second half of the most useful chain-extension in the course: an alkyl halide plus cyanide gives the nitrile (one carbon longer), and this step turns it into the acid. Ammonia leaves as the ammonium ion.',
          stereoHow: null, byprod: [{ smiles: 'N', name: 'ammonia, as the ammonium salt' }] }];
      }
      if (s.kind === 'epoxideAq') {
        return F.epoxOpen.variants(g0, s, 'epox_h2o', Object.assign({}, env, { pHVal: 1 }));
      }
      return [];
    },
    byproducts: () => [],
    condNotes() { return ['Dilute aqueous acid (H₃O⁺), room temperature for an acetal or an imine, heat for a nitrile.',
      'This is the same thing written as "workup" at the end of a Grignard or a hydride reduction — there it only protonates an alkoxide, but on a molecule with an acetal or an imine it does real chemistry, and that is worth noticing before you write it down.'] },
    stages: () => null,
  };

  /* --- silver oxide: the mild way to an acid from an aldehyde --- */
  F.aldOx = {
    rules: ['ag2o'], needs: 'carbonyl',
    sites(g, ruleId, env) {
      return carbonylSites(g).map(s => {
        if (s.cls !== 'aldehyde') return Object.assign(s, { inert: s.cls + ' — silver(I) oxide oxidises an ALDEHYDE only; a ketone has no hydrogen on the carbonyl carbon to lose' });
        s.score = 1;
        s.reason = 'the aldehyde — oxidised to the acid under conditions so mild that an alkene, an alcohol or a sulfide in the same molecule is untouched';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      const o2 = C().addAtomNear(g, 'O', atom(g, s.c.id));
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.c.id, o2.id]), added: [o2.id], share: 1,
        label: 'carboxylic acid (mild oxidation)',
        why: 'Silver(I) oxide takes an aldehyde to the acid and stops. It is the mild end of the oxidation scale: an alkene, a secondary alcohol, a sulfide or an aromatic ring in the same molecule all survive, where permanganate or chromic acid would attack them. ' +
          'This is the same chemistry as the Tollens test — the silver is reduced to the metal, which is what plates out as the silver mirror — used as a preparation rather than as a test.',
        stereoHow: null, byprod: [{ name: 'silver metal (the Tollens mirror)' }] }];
    },
    byproducts: () => [],
    condNotes() { return ['Ag₂O in aqueous THF, or Tollens’ reagent (Ag(NH₃)₂⁺) — room temperature.',
      'Chemoselective: choose it when the molecule contains something a stronger oxidant would ruin.'] },
    stages: () => null,
  };

  /* --- alkoxymercuration: Markovnikov ether from an alkene, with no rearrangement --- */
  F.alkoxyHg = {
    rules: ['alkoxyhg'], needs: 'alkene',
    sites(g, ruleId, env) {
      return alkeneSites(g).map(s => {
        if (s.aromatic) return Object.assign(s, { inert: 'aromatic C=C — mercury(II) does not add to a benzene ring' });
        s.score = 1; s.reason = 'the C=C — mercury(II) bridges it, then the alcohol opens the bridge at the more substituted carbon';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      const hA = hOn(g0, s.A), hB = hOn(g0, s.B);
      const put = hA <= hB ? s.A : s.B;               // Markovnikov: the more substituted carbon
      const bd = bondOf(g, s.A.id, s.B.id); if (bd) bd.order = 1;
      const o = C().addAtomNear(g, 'O', atom(g, put.id));
      const me = C().addAtomNear(g, 'C', o);
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.A.id, s.B.id, o.id, me.id]), added: [o.id, me.id], share: 1,
        label: 'ether — Markovnikov, with no rearrangement',
        why: 'Alkoxymercuration–demercuration: the Markovnikov ether from an alkene and an alcohol, in two steps and with no carbocation anywhere. ' +
          'Mercury(II) acetate forms a bridged mercurinium ion across the double bond; the alcohol opens that bridge by attacking the carbon that carries more positive charge — the more substituted one — and sodium borohydride then replaces the mercury with hydrogen. ' +
          'Why bother, when acid and an alcohol would also add across an alkene? Because the acid route goes through a real carbocation, which rearranges and eliminates; this one does not. It is the same relationship as oxymercuration to acid-catalysed hydration.',
        stereoHow: 'Anti addition across the double bond (the alcohol attacks the face opposite mercury), though the demercuration step is radical and scrambles that centre.',
        byprod: [{ name: 'Hg and the borate salts (from the NaBH₄ step)' }] }];
    },
    byproducts: () => [],
    condNotes() { return ['Hg(OAc)₂ in the alcohol as solvent, then NaBH₄ — two separate operations in one flask.',
      'Choose it over acid-catalysed addition when a carbocation would rearrange.'] },
    stages: () => null,
  };

  /* =====================================================================
     ENOLS, ENOLATES AND THE ALPHA CARBON

     The whole second half of the course runs on one idea: a hydrogen on the
     carbon NEXT to a carbonyl is acidic (pKa about 20 for a ketone, 25 for an
     ester, 11 for a 1,3-diketone), and taking it off makes a CARBON
     nucleophile. Every condensation in chapter 22 — aldol, Claisen, Michael,
     malonic ester, Robinson — is that nucleophile meeting an electrophile.

     So the enol and the enolate are made first-class species here. You can ask
     the app for them directly, and once you have one it reacts.
     ===================================================================== */
  function alphaSites(g) {
    const out = [];
    g.atoms.forEach(c => {
      if (c.element !== 'C') return;
      const dO = nb(g, c.id).find(n => n.atom.element === 'O' && n.bond.order === 2);
      if (!dO) return;
      const ns = nb(g, c.id);
      const or = ns.find(n => n.atom.element === 'O' && n.bond.order === 1);
      const nn = ns.find(n => n.atom.element === 'N' && n.bond.order === 1);
      const cls = or ? (hOn(g, or.atom) ? 'carboxylic acid' : 'ester') : nn ? 'amide' : (hOn(g, c) >= 1 ? 'aldehyde' : 'ketone');
      ns.forEach(n => {
        if (n.atom.element !== 'C' || n.bond.order !== 1) return;
        if (!bmSp3(g, n.atom) || hOn(g, n.atom) < 1) return;
        /* how many OTHER carbonyls flank this same carbon? Two makes it a 1,3-dicarbonyl,
           pKa about 11, and a mere carbonate or amine is base enough for it. */
        const flank = nb(g, n.atom.id).filter(m => m.atom.element === 'C' && m.atom.id !== c.id && bmDblO(g, m.atom.id)).length;
        const pka = flank >= 1 ? 11 : cls === 'ester' ? 25 : cls === 'amide' ? 30 : 20;
        out.push({ kind: 'alpha', c, o: dO.atom, alpha: n.atom, cls, pka, flank,
          sub: deg(g, n.atom), hind: hindrance(g, n.atom.id), nH: hOn(g, n.atom) });
      });
    });
    return out;
  }

  F.enol = {
    rules: ['tautomer', 'lda', 'enolate_base'], needs: 'carbonyl',
    sites(g, ruleId, env) {
      const kinetic = ruleId === 'lda';
      return alphaSites(g).map(s => {
        /* KINETIC vs THERMODYNAMIC, the most-tested contrast in the chapter.
           LDA is bulky, used cold, and one full equivalent: it takes the most ACCESSIBLE
           proton and the enolate it makes cannot go back, so the LESS substituted side wins.
           An alkoxide is small, used warm and catalytically: every enolate is formed and
           re-formed until the mixture settles on the MORE substituted (more stable, more
           highly conjugated) one. */
        s.score = kinetic ? (-1.2 * (s.sub - 1) - 0.5 * s.hind + 0.3 * s.nH)
                          : (1.0 * (s.sub - 1) - 0.2 * s.hind)
                  + (s.flank ? 4 : 0);                       // a doubly activated C–H wins outright
        s.label = 'α-C–H on C' + s.alpha.id + ' (' + (s.sub <= 1 ? 'primary' : s.sub === 2 ? 'secondary' : 'tertiary') + ', pKa ≈ ' + s.pka + ')';
        s.reason = 'the α-hydrogen — pKa about ' + s.pka + (s.flank ? ', because this carbon sits BETWEEN two carbonyls and the anion is delocalised over both of them' : ', acidic only because the anion is delocalised onto the carbonyl oxygen') +
          (kinetic ? '; LDA is bulky and cold, so it takes the most accessible one — here the ' + (s.sub <= 1 ? 'less' : 'more') + ' substituted side'
                   : '; with an alkoxide the enolates equilibrate, so the more substituted and more stable one is what you get');
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      const bd = bondOf(g, s.c.id, s.o.id); if (bd) bd.order = 1;
      const ab = bondOf(g, s.c.id, s.alpha.id); if (ab) ab.order = 2;
      const oAt = atom(g, s.o.id);
      const anion = ruleId !== 'tautomer';
      if (oAt) oAt.charge = anion ? -1 : 0;
      try { C().layoutGraph(g); } catch (e) {}
      const kinetic = ruleId === 'lda';
      return [{ graph: g, changed: new Set([s.c.id, s.o.id, s.alpha.id]), added: [], share: 1,
        label: anion ? (kinetic ? 'the KINETIC enolate (less substituted)' : 'the THERMODYNAMIC enolate (more substituted)') : 'the enol',
        nameHint: null,
        why: anion
          ? 'The base removes the α-hydrogen and what is left is the ENOLATE — one anion with its charge spread over the α-carbon and the oxygen. That is the point: the oxygen carries most of the negative charge (it is the more electronegative atom), but the CARBON is where it reacts, because carbon is the softer, more polarisable end. ' +
            (kinetic
              ? 'LDA gives the KINETIC enolate. Three things make it so: lithium diisopropylamide is very bulky, so it reaches the least hindered proton; it is used at −78 °C, so nothing has the energy to equilibrate; and it is used in full stoichiometric amount, so once a proton is taken there is no acid left to put it back. Deprotonation is irreversible, and the enolate you get is the one that formed FASTEST — from the less substituted side. '
              : 'An alkoxide gives the THERMODYNAMIC enolate. The base is small, the temperature is ambient, and it is catalytic: there is always some ketone left to protonate an enolate back, so the enolates interconvert until the mixture settles on the most STABLE one — the more substituted enolate, whose C=C carries more alkyl groups, exactly as with any alkene. ') +
            'Choosing between them is how you control which side of an unsymmetrical ketone reacts, and it is the reason both reagents are in the chapter.'
          : 'Keto–enol tautomerism. A trace of acid or base moves one hydrogen from the α-carbon to the oxygen, and the C=O becomes C=C–OH. It is an EQUILIBRIUM and it lies far to the left: acetone is about 0.0002 % enol at equilibrium, and cyclohexanone not much more. ' +
            'The enol matters anyway, for two reasons. It is the nucleophile in acid-catalysed α-halogenation, and it is why an α-stereocentre RACEMISES on standing with a trace of acid or base — the enol is flat, so when the proton comes back it can arrive on either face. A 1,3-diketone is different: there the enol is conjugated to the second carbonyl and hydrogen-bonded to it, and the enol can be the major form (2,4-pentanedione is about 80 % enol).',
        stereoHow: anion ? null : 'The enol is flat at the α-carbon, so any stereocentre there is lost — which is exactly why a trace of acid or base racemises an α-substituted ketone.',
        byprod: anion ? [{ smiles: ruleId === 'lda' ? 'CC(C)NC(C)C' : 'CCO', name: ruleId === 'lda' ? 'diisopropylamine (the conjugate acid of LDA)' : 'ethanol (the conjugate acid of the ethoxide)' }] : [] }];
    },
    byproducts: () => [],
    gate(g, ruleId, env) {
      if (ruleId === 'lda' && env.T > 253 && !env.TAuto) return 'LDA is used at −78 °C. Warmer than about −20 °C the enolates equilibrate and you no longer have the kinetic one — at which point an alkoxide would have been the cheaper reagent.';
      return null;
    },
    condNotes(ruleId, env) {
      return ruleId === 'lda'
        ? ['LDA (1.0 equiv), THF, −78 °C. Made in the flask from diisopropylamine and n-BuLi.',
           'Bulky, cold and stoichiometric — the three conditions for a kinetic enolate.']
        : ruleId === 'enolate_base'
          ? ['NaOEt in ethanol (or NaOH in water), room temperature, catalytic.',
             'Small, warm and catalytic — the three conditions for the thermodynamic enolate.']
          : ['A trace of acid or base is all it takes; no reagent is consumed.',
             'The equilibrium lies far towards the ketone, but the enol is the reactive species in α-halogenation and the reason α-stereocentres racemise.'];
    },
    stages: () => null,
  };

  /* --- alpha halogenation, in acid and in base, and the haloform cleavage --- */
  F.alphaHal = {
    rules: [], needs: 'carbonyl',
    sites(g, ruleId, env) {
      const basic = (env.pHVal !== undefined && env.pHVal >= 9);
      return alphaSites(g).map(s => {
        s.score = (basic ? 0.4 * s.nH : -0.4 * (s.sub - 1)) + (s.flank ? 2 : 0);
        s.label = 'α-C–H on C' + s.alpha.id;
        s.reason = basic
          ? 'the α-hydrogen — in BASE the enolate is the nucleophile, and each halogen put on makes the next α-hydrogen MORE acidic, so base-promoted halogenation does not stop at one'
          : 'the α-hydrogen — in ACID the enol is the nucleophile, and each halogen put on makes the next enolisation slower, so acid-promoted halogenation stops cleanly at ONE';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const basic = (env.pHVal !== undefined && env.pHVal >= 9);
      const xEl = ruleId === 'cl2' ? 'Cl' : ruleId === 'i2' ? 'I' : 'Br';
      const g = clone(g0);
      const nH = hOn(g0, s.alpha);
      /* in acid: exactly one halogen. In base: all of them — and if that carbon is a METHYL
         next to the carbonyl, the trihalomethyl group then leaves (the haloform reaction). */
      const many = basic ? nH : 1;
      const added = [];
      for (let i = 0; i < many; i++) { const x = C().addAtomNear(g, xEl, atom(g, s.alpha.id)); added.push(x.id); }
      let haloform = false, frag = null;
      if (basic && nH === 3 && deg(g0, s.alpha) <= 1) {
        /* CX3 on a carbonyl is a leaving group: hydroxide adds to the C=O and the
           trihalomethyl anion is expelled, giving the carboxylate and the haloform */
        haloform = true;
        frag = cutFragment(g, s.c.id, s.alpha.id);
        const oh = C().addAtomNear(g, 'O', atom(g, s.c.id)); oh.charge = -1; added.push(oh.id);
      }
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.alpha.id, s.c.id].concat(added)), added, share: 1,
        label: haloform ? 'the carboxylate (haloform reaction — the methyl group left as CH' + xEl + '₃)'
             : basic ? 'the poly-α-halo ketone (base does not stop at one)' : 'the α-halo ketone (one halogen, cleanly)',
        why: basic
          ? 'In BASE the nucleophile is the enolate, and the trouble is what the first halogen does to the molecule: it WITHDRAWS electrons, so the remaining α-hydrogens become even more acidic and the second and third halogenations are faster than the first. You cannot stop at one. ' +
            (haloform ? 'And with a methyl ketone it goes further still. Once all three hydrogens are replaced, CX₃ is a good enough leaving group that hydroxide can add to the carbonyl and expel it — so a methyl ketone in excess halogen and base gives the CARBOXYLATE plus chloroform, bromoform or iodoform. With iodine that yellow iodoform precipitate is a classical test for a methyl ketone (and, historically, for acetone in diabetic urine).'
                       : 'For the same reason, if you want ONE halogen you must use acid.')
          : 'In ACID the nucleophile is the ENOL, and this is the reason the reaction is well behaved: the slow step is enolisation, and the halogen already installed makes the carbonyl harder to enolise again. So acid-promoted halogenation stops cleanly after one. ' +
            'Note that the rate does not depend on the halogen concentration at all — that observation is how the mechanism was worked out, because it means the halogen is not involved in the slow step.',
        stereoHow: null,
        byprod: haloform ? [{ graph: bmTidy(frag), name: 'the haloform, CH' + xEl + '₃' }] : [{ smiles: xEl, name: 'H' + xEl }] }];
    },
    byproducts: () => [],
    condNotes(ruleId, env) {
      return ['Br₂ or Cl₂ in acetic acid for ONE halogen; Br₂ or I₂ with hydroxide for all of them.',
        'A methyl ketone with excess halogen and base gives the carboxylate and the haloform — the iodoform test.'];
    },
    stages: () => null,
  };

  /* --- Hell-Volhard-Zelinsky: the only way to halogenate a carboxylic acid's alpha carbon --- */
  F.hvz = {
    rules: [], needs: 'carbonyl',
    sites(g, ruleId, env) {
      return alphaSites(g).filter(s => s.cls === 'carboxylic acid').map(s => {
        s.score = -0.3 * (s.sub - 1);
        s.label = 'α-C–H of the acid (C' + s.alpha.id + ')';
        s.reason = 'the α-hydrogen of a carboxylic acid — a carboxylic acid barely enolises on its own, so PBr₃ first turns it into the acyl bromide, which enolises easily';
        return s;
      });
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      const br = C().addAtomNear(g, 'Br', atom(g, s.alpha.id));
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.alpha.id, br.id]), added: [br.id], share: 1,
        label: 'α-bromo carboxylic acid (Hell–Volhard–Zelinsky)',
        why: 'A carboxylic acid will not α-halogenate directly: its enol content is negligible, because the OH group is already feeding electrons into the carbonyl. The trick is to convert it first. ' +
          'Phosphorus tribromide (or a catalytic amount of red phosphorus with Br₂) makes the ACYL BROMIDE, and an acyl bromide enolises readily; that enol is brominated, and the aqueous workup hydrolyses the acyl bromide back to the acid. So the net result is bromine on the α-carbon of the acid. ' +
          'It matters because the α-bromo acid is a crossroads: displace the bromide with ammonia and you have an α-amino acid, with hydroxide an α-hydroxy acid, or eliminate to get an α,β-unsaturated acid.',
        stereoHow: 'Through a flat enol, so a new α-stereocentre is racemic.',
        byprod: [{ smiles: 'Br', name: 'HBr' }] }];
    },
    byproducts: () => [],
    condNotes() { return ['Br₂ with PBr₃ (or a catalytic amount of red phosphorus), then water.',
      'The α-bromo acid is a springboard: ammonia gives the α-amino acid, hydroxide the α-hydroxy acid.'] },
    stages: () => null,
  };

  /* --- heat as a reagent: decarboxylation of a beta-keto acid or a malonic acid --- */
  F.thermal = {
    rules: ['thermal'], needs: 'carbonyl', structural: true,
    sites(g, ruleId, env) {
      const out = [];
      g.atoms.forEach(c => {
        /* a carboxylic acid whose BETA carbon carries another carbonyl */
        const dO = nb(g, c.id).find(n => n.atom.element === 'O' && n.bond.order === 2);
        const oh = nb(g, c.id).find(n => n.atom.element === 'O' && n.bond.order === 1 && hOn(g, n.atom) >= 1);
        if (c.element !== 'C' || !dO || !oh) return;
        const alpha = nb(g, c.id).find(n => n.atom.element === 'C' && n.bond.order === 1);
        if (!alpha) return;
        const beta = nb(g, alpha.atom.id).find(m => m.atom.element === 'C' && m.atom.id !== c.id && bmDblO(g, m.atom.id));
        if (!beta) return;
        out.push({ kind: 'betaketoacid', c, o: dO.atom, oh: oh.atom, alpha: alpha.atom, beta: beta.atom, score: 1,
          label: 'β-keto acid / malonic acid (C' + c.id + ')',
          reason: 'a carboxyl with another carbonyl two atoms away — this is the one arrangement that loses CO₂ on gentle heating, and the reason is the six-membered ring the molecule can make of its own hydrogen bond' });
      });
      return out;
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      /* the whole COOH leaves as CO2; the alpha carbon keeps a hydrogen */
      try { cutFragment(g, s.alpha.id, s.c.id); } catch (e) {}
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.alpha.id]), added: [], share: 1,
        label: 'decarboxylated — CO₂ lost',
        why: 'A β-keto acid (or a malonic acid) loses carbon dioxide when it is merely warmed, and almost no other carboxylic acid does. The reason is geometric: the acid’s O–H can hydrogen-bond to the OTHER carbonyl’s oxygen through a six-membered ring, and in that arrangement the C–C bond breaks and the proton transfers in one concerted motion. What comes off is CO₂, and what is left is the ENOL of the product, which tautomerises at once to the ketone or acid. ' +
          'This is the step that makes the malonic ester and acetoacetic ester syntheses work. You alkylate the doubly activated carbon (easy, because its pKa is 11), hydrolyse the esters to the diacid, and then simply heat: one carboxyl leaves as gas and you are left with the mono-substituted product. The second ester group was never meant to survive — it was there only to make the alkylation possible.',
        stereoHow: null,
        byprod: [{ smiles: 'O=C=O', name: 'CO₂ — leaves as a gas' }] }];
    },
    byproducts: () => [],
    condNotes() { return ['Just heat — 100–150 °C, often in the same pot as the preceding hydrolysis.',
      'Only a β-keto acid or a malonic acid does this. An ordinary carboxylic acid needs copper and quinoline at 200 °C or more.'] },
    stages: () => null,
  };

  /* =====================================================================
     POLYMERISATION

     A polymer is not a new kind of chemistry: it is one reaction you already
     know, repeated. Chain growth is addition to a C=C, over and over, and the
     only question is what carries the chain — a radical, a cation, an anion or
     a metal. Step growth is an acyl substitution, over and over, between two
     molecules that each have TWO ends.

     The app draws three repeat units. A real chain has thousands; three is
     enough to see the pattern and to check the regiochemistry, which is what
     matters in an exam: the units join HEAD TO TAIL, so the substituents end
     up on alternate carbons.
     ===================================================================== */
  /* paste a whole graph into another, returning the id map */
  function mergeGraph(g, frag) {
    const off = Math.max(g.nextId || 0, ...g.atoms.map(a => a.id + 1), 1);
    const map = new Map();
    frag.atoms.forEach(a => { map.set(a.id, off + a.id); g.atoms.push(Object.assign({}, a, { id: off + a.id, x: (a.x || 0) + 60, y: (a.y || 0) })); });
    frag.bonds.forEach(b => g.bonds.push(Object.assign({}, b, { a: map.get(b.a), b: map.get(b.b) })));
    g.nextId = off + Math.max(...frag.atoms.map(a => a.id)) + 1;
    return map;
  }
  /* the C=C of a vinyl monomer, with the more substituted carbon identified */
  function vinylOf(g) {
    let best = null;
    g.atoms.forEach(a => {
      if (a.element !== 'C') return;
      nb(g, a.id).forEach(n => {
        if (n.atom.element !== 'C' || n.bond.order !== 2 || a.id >= n.atom.id || isAromatic(g, a.id)) return;
        const hA = hOn(g, a), hB = hOn(g, n.atom);
        const head = hA >= hB ? a : n.atom;            // the CH2 end
        const tail = hA >= hB ? n.atom : a;            // the substituted end
        const cand = { head, tail, subs: 4 - hOn(g, tail) - 1 };
        if (!best || cand.subs > best.subs) best = cand;
      });
    });
    return best;
  }
  function chainGrow(mono, n, alt) {
    /* build H-[CH2-CHX]n-H, head to tail; `alt` is an optional second monomer */
    const units = [];
    for (let i = 0; i < n; i++) units.push(clone(alt && i % 2 === 1 ? alt : mono));
    const g = { atoms: [], bonds: [], nextId: 1 };
    let prevTail = null;
    units.forEach(u => {
      const v = vinylOf(u);
      if (!v) return;
      const bd = bondOf(u, v.head.id, v.tail.id); if (bd) bd.order = 1;
      const map = mergeGraph(g, u);
      const head = map.get(v.head.id), tail = map.get(v.tail.id);
      if (prevTail !== null) bondOf2(g, prevTail, head, 1);
      prevTail = tail;
    });
    try { C().layoutGraph(g); } catch (e) {}
    return g;
  }

  F.polymer = {
    rules: ['radpoly', 'cationic', 'anionic', 'ziegler'], needs: 'alkene', structural: true,
    sites(g, ruleId, env) {
      const v = vinylOf(g);
      if (!v) return [];
      /* which initiator suits which monomer is a real question, not a detail: the growing end
         has to be stabilised by the substituent, and the three mechanisms want opposite things */
      const tailNs = nb(g, v.tail.id);
      const ewg = tailNs.some(n => (n.atom.element === 'C' && (bmDblO(g, n.atom.id) || nb(g, n.atom.id).some(m => m.atom.element === 'N' && m.bond.order === 3))) || (n.atom.element === 'N' && n.atom.charge === 1));
      const edg = tailNs.some(n => (n.atom.element === 'O') || (n.atom.element === 'C' && isAromatic(g, n.atom.id)) || (n.atom.element === 'C' && hOn(g, n.atom) === 3));
      let ok = true, note = '';
      if (ruleId === 'anionic' && !ewg) { ok = false; note = 'anionic polymerisation needs an ELECTRON-POOR alkene: the chain end is a carbanion, and only an ester, nitrile or nitro group next to it can stabilise that. Styrene is borderline; ethylene and isobutylene will not go at all.'; }
      if (ruleId === 'cationic' && !edg) { ok = false; note = 'cationic polymerisation needs an ELECTRON-RICH alkene: the chain end is a carbocation, so it wants an alkyl, aryl or alkoxy group to stabilise it. Isobutylene and vinyl ethers are ideal; acrylonitrile is hopeless.'; }
      if (!ok) return [{ kind: 'vinyl', c: v.tail, inert: note, label: 'the C=C' }];
      return [{ kind: 'vinyl', v, c: v.tail, score: 1, label: 'the C=C (monomer)',
        reason: 'the C=C — it opens and joins to the next monomer, thousands of times over' + (ewg ? '; the withdrawing group here stabilises an anionic chain end' : edg ? '; the donating group here stabilises a cationic chain end' : '') }];
    },
    variants(g0, s, ruleId, env) {
      const g = chainGrow(g0, 3, null);
      const mech = { radpoly: 'radical', cationic: 'cationic', anionic: 'anionic', ziegler: 'Ziegler–Natta' }[ruleId] || 'radical';
      return [{ graph: g, changed: new Set(g.atoms.map(a => a.id)), added: [], share: 1,
        label: 'the polymer — three repeat units shown (' + mech + ' chain growth)',
        why: 'Chain growth: an addition to the C=C, repeated. Three steps that never change — INITIATION makes the reactive end, PROPAGATION adds monomer after monomer to it, TERMINATION destroys it — and the monomer is consumed steadily throughout, so the molecular weight is set by how long a chain survives before it terminates. ' +
          (mech === 'radical'
            ? 'A peroxide breaks homolytically on warming and the radical adds to the alkene; the new radical adds to the next monomer, and so on. Termination is two radicals meeting — either COMBINING to one long chain or DISPROPORTIONATING, where one takes a hydrogen from the other and you get one saturated end and one alkene end. Radical chains also transfer to themselves, which is why low-density polyethylene is BRANCHED and soft while the Ziegler–Natta product is linear and stiff.'
            : mech === 'cationic'
              ? 'A Lewis acid with a trace of water makes a proton, which adds to the alkene to give the more stable carbocation; that cation adds to the next monomer. It only works on electron-rich alkenes, because every chain end is a carbocation — isobutylene (giving butyl rubber) is the classic case.'
              : mech === 'anionic'
                ? 'An anionic initiator (butyllithium, or a trace of hydroxide) adds to the electron-poor alkene and the carbanion that results adds to the next monomer. The chain ends are anions and there is nothing for them to react with, so they stay alive — LIVING polymerisation — which is why this is the route to block copolymers and to very narrow molecular-weight distributions. It is also why cyanoacrylate glue sets on your fingers: the moisture is the initiator.'
                : 'A titanium(IV) chloride / alkylaluminium catalyst inserts each monomer into a Ti–C bond on a metal surface, and the surface controls the geometry. Two consequences: the chain is LINEAR (no branching, so the polymer crystallises, and that is high-density polyethylene) and it is STEREOREGULAR — isotactic polypropylene has every methyl group on the same side and is a useful solid, while the atactic polymer is a sticky mess. Ziegler and Natta shared the 1963 Nobel Prize for exactly this.') + ' ' +
          'Note the regiochemistry: HEAD TO TAIL. Each new monomer adds so that its CH₂ end joins the previous unit’s substituted carbon, because that puts the radical (or ion) on the carbon that can stabilise it. So the substituents end up on alternate carbons, never adjacent.',
        stereoHow: ruleId === 'ziegler' ? 'The catalyst surface sets the configuration of every stereocentre: isotactic (all the same side). A radical chain gives atactic — random — because the chain end is flat.' : 'The chain end is flat, so each new stereocentre is set at random: the polymer is ATACTIC, and that is why radical polypropylene is useless as a plastic.' }];
    },
    byproducts: () => [],
    condNotes(ruleId) {
      return ruleId === 'radpoly' ? ['A peroxide or AIBN (0.1–1 %), 60–150 °C, often under pressure.', 'Chain transfer gives branches; branches stop the polymer crystallising, which is what makes LDPE soft.']
        : ruleId === 'cationic' ? ['BF₃ with a trace of water, or H₂SO₄, at low temperature.', 'Only electron-rich alkenes: isobutylene, vinyl ethers, styrene.']
        : ruleId === 'anionic' ? ['n-BuLi or NaOH, in a dry aprotic solvent.', 'Only electron-poor alkenes. The chain ends stay alive, so a second monomer added later gives a block copolymer.']
        : ['TiCl₄ with Et₂AlCl or Et₃Al, a heterogeneous catalyst, mild temperature and pressure.', 'Linear and stereoregular — the reason polypropylene is a structural material at all.'];
    },
    stages: () => null,
  };

  /* --- turning a carboxylic acid into something reactive --- */
  F.acylCl = {
    rules: ['soCl2', 'pcl5', 'ch2n2'], needs: 'carbonyl',
    sites(g, ruleId, env) {
      const out = [];
      g.atoms.forEach(c => {
        if (c.element !== 'C') return;
        const dO = nb(g, c.id).find(n => n.atom.element === 'O' && n.bond.order === 2);
        const oh = nb(g, c.id).find(n => n.atom.element === 'O' && n.bond.order === 1 && hOn(g, n.atom) >= 1);
        if (!dO || !oh) return;
        out.push({ kind: 'acid', c, o: dO.atom, oh: oh.atom, score: 1,
          label: 'carboxylic acid (C' + c.id + ')',
          reason: ruleId === 'ch2n2' ? 'the carboxylic acid — its O–H is acidic enough to protonate diazomethane, and that is the whole mechanism'
            : 'the carboxylic acid — the OH is turned into a leaving group, which is the only way to make this carbon reactive towards a weak nucleophile' });
      });
      /* thionyl chloride also converts an ALCOHOL to an alkyl chloride, which is its other job */
      if (ruleId === 'soCl2') {
        alcoholSites(g).forEach(a => { a.kind = 'alcohol'; a.score = 0.4; a.label = 'alcohol (C' + a.c.id + ')';
          a.reason = 'the alcohol — SOCl₂ converts it to the alkyl chloride, with inversion; this is its other standard use'; out.push(a); });
      }
      return out;
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      if (s.kind === 'alcohol') {
        const o = s.o || nb(g0, s.c.id).find(n => n.atom.element === 'O' && hOn(g0, n.atom) >= 1);
        const oid = (s.o && s.o.id) !== undefined ? s.o.id : (o && o.atom ? o.atom.id : null);
        if (oid === null) return [];
        try { C().removeAtom(g, oid); } catch (e) {}
        const cl = C().addAtomNear(g, 'Cl', atom(g, s.c.id));
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.c.id, cl.id]), added: [cl.id], share: 1, label: 'alkyl chloride',
          why: 'Thionyl chloride turns the alcohol into an alkyl chloride. The oxygen attacks sulfur first, which converts the terrible leaving group (hydroxide) into an excellent one (chlorosulfite); chloride then displaces it. Both by-products are gases — SO₂ and HCl — so the reaction drives itself and the workup is trivial, which is why SOCl₂ is preferred to HCl for this.',
          stereoHow: 'Backside attack by chloride: INVERSION at that carbon (in pyridine; in ether an internal-return pathway can give retention, which is a classic exam subtlety).',
          byprod: [{ smiles: 'O=S=O', name: 'SO₂ (a gas)' }, { smiles: 'Cl', name: 'HCl (a gas)' }] }];
      }
      if (ruleId === 'ch2n2') {
        const me = C().addAtomNear(g, 'C', atom(g, s.oh.id));
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.oh.id, me.id]), added: [me.id], share: 1, label: 'methyl ester',
          why: 'Diazomethane methylates a carboxylic acid, and the mechanism is two steps with nothing else involved. The acid protonates diazomethane to give the methyldiazonium ion; the carboxylate then does an Sₙ2 on that methyl group and nitrogen leaves as N₂ gas. ' +
            'It is the mildest esterification there is — instant, at 0 °C, quantitative, with no acid catalyst and no water to remove — and it is selective for acids, leaving alcohols and phenols alone. The catch is the reagent: diazomethane is a toxic, explosive gas that has to be generated and used in solution behind a shield, which is why it belongs to research rather than to teaching labs. TMS-diazomethane is the safe modern stand-in.',
          stereoHow: null, byprod: [{ smiles: 'N#N', name: 'N₂ — leaves as a gas, which is what makes it irreversible' }] }];
      }
      /* the acyl chloride */
      try { C().removeAtom(g, s.oh.id); } catch (e) {}
      const cl = C().addAtomNear(g, 'Cl', atom(g, s.c.id));
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.c.id, cl.id]), added: [cl.id], share: 1, label: 'acyl chloride',
        why: 'The acyl chloride, and the reason for making one is the whole logic of chapter 21. A carboxylic acid is the LEAST reactive of the acid derivatives: its OH is a hopeless leaving group, and any decent nucleophile is basic enough to simply take the acidic proton instead, leaving an unreactive carboxylate. ' +
          'Converting it to the acyl chloride puts the carbon at the TOP of the reactivity series, and from there everything downhill is easy — an ester with an alcohol, an amide with an amine, an anhydride with a carboxylate, a ketone with a cuprate. ' +
          'Mechanistically SOCl₂ (or PCl₃, or PCl₅, or oxalyl chloride) does the same thing in each case: the acid oxygen attacks the electrophilic sulfur or phosphorus, making the OH into a good leaving group, and chloride then displaces it through the tetrahedral intermediate. Every by-product is a gas or is easily washed out, and the acyl chloride is normally used in the next step without purification — it hydrolyses in air.',
        stereoHow: null,
        byprod: ruleId === 'pcl5' ? [{ smiles: 'O=P(Cl)(Cl)Cl', name: 'POCl₃' }, { smiles: 'Cl', name: 'HCl' }]
          : [{ smiles: 'O=S=O', name: 'SO₂ (a gas)' }, { smiles: 'Cl', name: 'HCl (a gas)' }] }];
    },
    byproducts: () => [],
    condNotes(ruleId) {
      return ruleId === 'ch2n2' ? ['CH₂N₂ in ether, 0 °C — added until the yellow colour persists, which is the endpoint.',
          'Explosive and toxic: generated in solution and never distilled. TMS-diazomethane in methanol is the safe substitute.']
        : ['SOCl₂ (neat or in toluene), reflux; or PCl₃, PCl₅, or oxalyl chloride at room temperature.',
           'Used immediately in the next step — an acyl chloride hydrolyses in moist air.'];
    },
    stages: () => null,
  };

  /* --- the simplest reaction in the course, and the one most often written down wrong --- */
  F.acidBase = {
    rules: [], needs: 'any',
    sites(g, ruleId, env) {
      const out = [];
      g.atoms.forEach(a => {
        /* a carboxylic acid: pKa about 5, so ANY of hydroxide, bicarbonate or carbonate does it */
        if (a.element === 'O' && hOn(g, a) === 1) {
          const c = nb(g, a.id).find(n => n.atom.element === 'C');
          if (c && bmDblO(g, c.atom.id)) { out.push({ kind: 'acid', o: a, score: 3, label: 'carboxylic acid O–H (O' + a.id + ')',
            reason: 'the carboxylic acid — pKa about 5, so even sodium bicarbonate (pKa 6.4) deprotonates it, and that is the classical test that distinguishes an acid from a phenol' }); return; }
          if (c && isAromatic(g, c.atom.id)) out.push({ kind: 'phenol', o: a, score: 1.5, label: 'phenol O–H (O' + a.id + ')',
            reason: 'the phenol — pKa about 10: hydroxide deprotonates it but bicarbonate does not, which is how a phenol is separated from a carboxylic acid' });
        }
        /* an ammonium: hydroxide takes the proton back off */
        if (a.element === 'N' && a.charge === 1 && hOn(g, a) >= 1)
          out.push({ kind: 'ammonium', n: a, score: 2.5, label: 'ammonium N–H (N' + a.id + ')',
            reason: 'the ammonium ion — pKa about 10 for an alkylammonium, so hydroxide takes the proton off and frees the amine' });
      });
      return out;
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      if (s.kind === 'ammonium') {
        const n = atom(g, s.n.id); if (n) { n.charge = 0; n.hAdd = 0; }
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.n.id]), added: [], share: 1, label: 'the free amine',
          why: 'Hydroxide takes the proton off the ammonium ion and the free amine is released. ' +
            'This is not a footnote: it is how amines are PURIFIED and separated. Shake a mixture with dilute HCl and the amine goes into the water as its ammonium salt while everything neutral stays in the ether; make the aqueous layer basic again and the amine comes back out. The same trick in reverse separates a carboxylic acid from a neutral compound. ' +
            'It is also why amine drugs are sold as hydrochloride salts — the salt is a crystalline, water-soluble solid, while the free base is often an oil.',
          stereoHow: null, byprod: [{ smiles: 'O', name: 'water' }] }];
      }
      const o = atom(g, s.o.id);
      if (o) { o.charge = -1; o.hAdd = 0; }
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.o.id]), added: [], share: 1,
        label: s.kind === 'acid' ? 'the carboxylate salt' : 'the phenoxide salt',
        why: (s.kind === 'acid'
          ? 'The acid is deprotonated to the carboxylate. The pKa is about 5 — roughly ten powers of ten more acidic than an alcohol — and the whole reason is what happens AFTER the proton leaves: the negative charge is shared equally between two oxygens, which are identical by resonance, so the anion is genuinely stabilised rather than merely tolerated. An alkoxide has nowhere to put the charge. ' +
            'The practical consequence is the solubility test. A carboxylic acid dissolves in sodium bicarbonate solution with visible fizzing (the carbonic acid released decomposes to CO₂); a phenol does not, because at pKa 10 it is not acidic enough to protonate bicarbonate. That one observation separates the two classes in a teaching lab, and it is the basis of the acid–base extraction.'
          : 'The phenol is deprotonated to the phenoxide. At pKa 10 a phenol is about a million times more acidic than an ordinary alcohol, because the negative charge on the oxygen is delocalised into the ring — onto the two ortho positions and the para one. Any electron-withdrawing group on those positions helps further: 4-nitrophenol is pKa 7, and 2,4,6-trinitrophenol (picric acid) is pKa 0.4, which is why it is called an acid at all. ' +
            'But it is not acidic enough for bicarbonate, and that is the test: hydroxide dissolves a phenol, bicarbonate does not.'),
        stereoHow: null, byprod: [{ smiles: 'O', name: 'water' }] }];
    },
    byproducts: () => [],
    condNotes() { return ['NaOH, NaHCO₃ or Na₂CO₃ in water, room temperature — instantaneous.',
      'The classical separation: bicarbonate takes out the carboxylic acid, hydroxide then takes out the phenol, and what is left in the organic layer is neutral.'] },
    stages: () => null,
  };

  F.cyclise = {
    rules: ['cyclise'], needs: 'any', structural: true,
    sites(g) {
      let P = null;
      try { P = bmProfile(g); } catch (e) { return []; }
      const out = [];
      const add = (kind, nuAt, eC, lg, word, extra) => {
        const d = bondPath(g, nuAt.id, eC.id);
        if (d < 2 || d > 11) return;                        // 2 bonds apart is a 3-ring
        const size = d + 1;
        /* already in a ring together? then this bond would be a bridge, not a new ring */
        let already = false;
        try { already = C().ringSystemsOf(g, C().ringAtomsOf(g)).some(r => r.indexOf(nuAt.id) >= 0 && r.indexOf(eC.id) >= 0); } catch (e) {}
        if (already) return;
        out.push(Object.assign({ kind: 'cyclise', what: kind, atoms: [nuAt.id, eC.id], c: eC, nu: nuAt, lg: lg || null,
          size, score: ringPref(size) + (extra && extra.bonus ? extra.bonus : 0),
          label: word + ' → ' + size + '-membered ring (N' + nuAt.id + ' to C' + eC.id + ')',
          reason: size === 5 || size === 6 ? 'a ' + size + '-membered ring: the ends meet with almost no strain and the entropy cost is small, which is why five and six dominate' :
            size === 3 ? 'a 3-membered ring: strained, but the two ends are so close that it forms anyway when nothing better is available' :
            size === 4 ? 'a 4-membered ring: both strained AND a poor approach angle for backside attack — the worst of the options' :
            size === 7 ? 'a 7-membered ring: possible, but slower — the chain has to give up a lot of freedom to reach round' :
            'a ' + size + '-membered ring: at this size the two ends rarely find each other, and the reaction needs high dilution or it joins DIFFERENT molecules instead (oligomer, not ring)' }, extra || {}));
      };
      /* acyl + N–H or O–H in the same molecule: a lactam or a lactone */
      P.acyl.forEach(ac => {
        P.amine.forEach(n => add('lactam', n.n, ac.c, ac.lg, 'lactam closure', { ac, nuKind: 'amine', bonus: 0.5 }));
        P.alcohol.forEach(n => add('lactone', n.o, ac.c, ac.lg, 'lactone closure', { ac, nuKind: 'alcohol' }));
        P.phenol.forEach(n => add('lactone', n.o, ac.c, ac.lg, 'aryl lactone closure', { ac, nuKind: 'phenol' }));
      });
      /* alkyl halide + a nucleophile in the same molecule */
      P.halide.forEach(h => {
        if (h.deg >= 3) return;
        P.amine.forEach(n => add('cyclic amine', n.n, h.c, h.x, 'cyclic amine closure', { h, nuKind: 'amine', bonus: 0.4 }));
        P.alcohol.forEach(n => add('cyclic ether', n.o, h.c, h.x, 'cyclic ether closure', { h, nuKind: 'alcohol' }));
        P.phenol.forEach(n => add('cyclic ether', n.o, h.c, h.x, 'benzo-fused ether closure', { h, nuKind: 'phenol', bonus: 0.3 }));
        P.thiol.forEach(n => add('cyclic thioether', n.s, h.c, h.x, 'cyclic thioether closure', { h, nuKind: 'thiol', bonus: 0.3 }));
      });
      /* an activated aryl halide with the nucleophile on a tether: intramolecular Sₙ Ar */
      P.halideSp2.forEach(h => {
        if (!h.act || h.otf) return;
        P.amine.forEach(n => add('cyclic amine', n.n, h.c, h.x, 'intramolecular Sₙ Ar', { h, sp2: true, nuKind: 'amine', bonus: 0.4 }));
        P.phenol.forEach(n => add('cyclic ether', n.o, h.c, h.x, 'intramolecular Sₙ Ar', { h, sp2: true, nuKind: 'phenol' }));
        P.alcohol.forEach(n => add('cyclic ether', n.o, h.c, h.x, 'intramolecular Sₙ Ar', { h, sp2: true, nuKind: 'alcohol' }));
      });
      return out;
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      const by = [];
      const eC = s.c.id, nuAt = s.nu.id;
      let stages = null;
      if (s.ac) {
        /* the tetrahedral intermediate, then the leaving group goes */
        try {
          const t = clone(g);
          bondOf2(t, nuAt, eC, 1);
          const dO = t.atoms.find(a => a.id === s.ac.o.id);
          if (dO) { const bd = bondOf(t, eC, dO.id); if (bd) { bd.order = 1; dO.charge = -1; } }
          bmTidy(t);
          stages = [{ label: 'tetrahedral intermediate — the ring is closed, the leaving group not yet gone', graph: t }];
        } catch (e) {}
        bondOf2(g, nuAt, eC, 1);
        const kind = s.ac.kind;
        const lg = s.ac.lg.id;
        if (kind === 'acid') { try { C().removeAtom(g, lg); } catch (e) {} by.push({ smiles: 'O', name: 'water' }); }
        else if (kind === 'acyl chloride') { try { C().removeAtom(g, lg); } catch (e) {} by.push({ smiles: 'Cl', name: 'HCl' }); }
        else {
          const frag = cutFragment(g, eC, lg);
          by.push({ graph: bmTidy(frag), name: 'the alcohol (from the ester)' });
        }
      } else if (s.h && s.h.sulf) {
        bondOf2(g, nuAt, eC, 1);
        const frag = cutFragment(g, eC, s.lg.id);
        const a0 = g.atoms.find(x => x.id === nuAt);
        if (a0 && a0.charge === -1) a0.charge = 0;
        by.push({ graph: bmTidy(frag), name: 'the ' + s.h.lgWord + ' anion — the leaving group' });
      } else {
        const xEl = s.lg.element;
        bondOf2(g, nuAt, eC, 1);
        try { C().removeAtom(g, s.lg.id); } catch (e) {}
        const a = g.atoms.find(x => x.id === nuAt);
        if (a && a.charge === -1) a.charge = 0;
        by.push({ smiles: xEl === 'F' ? null : xEl, name: 'H' + xEl + ' — taken up by the base' });
      }
      try { C().layoutGraph(g); } catch (e) {}
      const ringWord = s.size + '-membered';
      return [{ graph: g, changed: new Set([eC, nuAt]), added: [], share: 1, byprod: by,
        label: s.what + ' (' + ringWord + ', intramolecular)',
        why: 'The same substitution the engine would do between two molecules, except that both halves are already in this one, so it closes a ring. ' +
          (s.ac ? 'The ' + s.nuKind + ' adds to the carbonyl carbon, the tetrahedral intermediate collapses and ' + by[0].name + ' leaves. '
                : s.sp2 ? 'The ' + s.nuKind + ' attacks the activated aromatic carbon (addition–elimination through the Meisenheimer anion) and the halide leaves. '
                : 'The ' + s.nuKind + ' displaces the halide by Sₙ2 from the back side. ') +
          'Ring size is what decides whether this happens at all: ' + s.reason + '. ' +
          (s.size >= 8 ? 'At this size the honest advice is high dilution and slow addition, and to expect polymer if you get it wrong.' :
           'The intramolecular version is much faster than the same reaction between two molecules, because the two ends are held together — the effective concentration of the nucleophile at that carbon is enormous.'),
        stereoHow: null, _stages: stages }];
    },
    byproducts: () => [],
    condNotes(ruleId, env) {
      return ['Intramolecular: what usually starts it is a mild base (K₂CO₃, NaH) for a halide, or heat — with a coupling agent (EDC, T3P) for a lactam from the free acid.',
        'Dilution matters here in a way it does not for two separate molecules: the more dilute, the more a ring beats a chain. For anything above a seven-membered ring, slow addition into a large volume is the difference between product and polymer.'];
    },
    stages(g0, s, out) { return out._stages || null; },
  };

  /* =====================================================================
     DIAZONIUM CHEMISTRY — the one place an aromatic ring accepts almost
     anything

     Electrophilic aromatic substitution can put a small set of groups on a
     ring, and nucleophilic aromatic substitution needs the ring to be
     electron-poor first. Neither can give you an aryl iodide, an aryl
     fluoride, an aryl nitrile, or a phenol from an ordinary benzene ring. The
     diazonium salt can give you all of them, from one starting material, and
     that is why it has survived in the syllabus and in the plant.

     Two steps, and the first is always the same. A primary AROMATIC amine and
     nitrous acid (made in the flask from NaNO2 and HCl, which is why it must
     be cold) give the diazonium ion, in which N2 sits on the ring as the best
     leaving group in aromatic chemistry. Then whatever nucleophile or copper
     salt is present decides what replaces it.

     Cold is not a detail: above about 5 degrees an aryl diazonium loses N2 on
     its own and the yield leaves as gas. An ALIPHATIC diazonium does that
     instantly at any temperature, which is why this chemistry is for aromatic
     amines only.
     ===================================================================== */
  function arylAmineSites(g) {
    const out = [];
    g.atoms.forEach(a => {
      if (a.element !== 'N' || a.charge) return;
      const ns = nb(g, a.id);
      const cs = ns.filter(n => n.atom.element === 'C' && n.bond.order === 1);
      const h = hOn(g, a);
      /* an amide or sulfonamide nitrogen is not basic enough to be diazotised */
      if (ns.some(n => n.atom.element === 'C' && bmDblO(g, n.atom.id)) ||
          ns.some(n => n.atom.element === 'S' && nb(g, n.atom.id).filter(m => m.atom.element === 'O' && m.bond.order === 2).length === 2)) return;
      const onRing = cs.filter(n => isAromatic(g, n.atom.id));
      if (h >= 2 && cs.length === 1) out.push({ n: a, c: cs[0].atom, aryl: onRing.length === 1, kind: 'primary' });
      else if (h === 1 && cs.length === 2) out.push({ n: a, c: cs[0].atom, aryl: onRing.length >= 1, kind: 'secondary' });
      else if (h === 0 && cs.length === 3 && onRing.length === 1) out.push({ n: a, c: onRing[0].atom, aryl: true, kind: 'tertiary' });
    });
    return out;
  }
  function diazoniumSites(g) {
    const out = [];
    g.atoms.forEach(a => {
      if (a.element !== 'N' || a.charge !== 1) return;
      const t = nb(g, a.id).find(n => n.atom.element === 'N' && n.bond.order === 3 && nb(g, n.atom.id).length === 1);
      const c = nb(g, a.id).find(n => n.atom.element === 'C');
      if (t && c) out.push({ n1: a, n2: t.atom, c: c.atom, aryl: isAromatic(g, c.atom.id) });
    });
    return out;
  }

  F.diazo = {
    rules: ['nano2'], needs: 'amine', structural: true,
    sites(g) {
      return arylAmineSites(g).map(s => {
        if (s.kind === 'tertiary') return Object.assign(s, { label: 'tertiary aromatic amine (N' + s.n.id + ')',
          inert: 'a tertiary aromatic amine has no N–H to lose, so it cannot become a diazonium — nitrous acid attacks its RING instead, para, to give the C-nitroso compound' });
        if (s.kind === 'secondary') return Object.assign(s, { label: 'secondary amine (N' + s.n.id + ')', score: 0.5, nitroso: true,
          reason: 'a secondary amine — one N–H is not enough for a diazonium, so it stops at the N-NITROSAMINE (and those are potent carcinogens: this is the reaction behind the nitrosamine limits in drug substances)' });
        if (!s.aryl) return Object.assign(s, { label: 'primary aliphatic amine (N' + s.n.id + ')',
          inert: 'a primary ALIPHATIC amine does form the diazonium, and then loses N₂ instantly at any temperature to give a carbocation — the products are an alcohol, alkenes and rearranged material in a mixture nobody wants. Only an aromatic diazonium is stable enough to use' });
        return Object.assign(s, { label: 'primary aromatic amine (N' + s.n.id + ')', score: 2,
          reason: 'a primary aromatic amine — the ring stabilises the diazonium ion by conjugation, which is the whole reason this chemistry exists' });
      });
    },
    variants(g0, s, ruleId, env) {
      const g = clone(g0);
      if (s.nitroso) {
        const o = C().addAtomNear(g, 'N', atom(g, s.n.id));
        const o2 = C().addAtomNear(g, 'O', atom(g, o.id));
        const bd = bondOf(g, o.id, o2.id); if (bd) bd.order = 2;
        try { C().layoutGraph(g); } catch (e) {}
        return [{ graph: g, changed: new Set([s.n.id, o.id, o2.id]), added: [o.id, o2.id], share: 1, label: 'N-nitrosamine',
          why: 'Nitrous acid nitrosates the nitrogen: the amine attacks the nitrosonium ion (NO⁺, which is what HNO₂ becomes in acid) and the N–N=O compound is left. With a PRIMARY amine the same intermediate goes on — it tautomerises to the diazohydroxide and loses water to give the diazonium. A secondary amine has no second N–H, so it stops here. ' +
            'This is not an academic point: N-nitrosamines are potent carcinogens, and the recalls of sartan and ranitidine medicines came from exactly this reaction happening where nobody intended it — a secondary amine and a nitrite impurity in the same process stream.',
          stereoHow: null }];
      }
      const nAt = atom(g, s.n.id);
      const nNew = C().addAtomNear(g, 'N', nAt);
      const bd = bondOf(g, s.n.id, nNew.id); if (bd) bd.order = 3;
      nAt.charge = 1;
      try { C().layoutGraph(g); } catch (e) {}
      return [{ graph: g, changed: new Set([s.n.id, nNew.id]), added: [nNew.id], share: 1, label: 'aryl diazonium salt',
        why: 'Diazotisation. NaNO₂ and HCl make nitrous acid in the flask, and in acid that becomes the nitrosonium ion NO⁺. The amine nitrogen attacks it; the N-nitrosamine that results tautomerises to the diazohydroxide Ar–N=N–OH, and one more protonation lets water leave. What is left is the diazonium ion, with N₂ sitting on the ring. ' +
          'Keep it between 0 and 5 °C. Warmer, and the salt decomposes on its own — N₂ leaves, the aryl cation grabs whatever water is there, and your yield walks out of the flask as gas and phenol. Diazonium salts are also not isolated dry on any scale: solid aryl diazonium salts are shock-sensitive explosives. Make it cold, use it cold, in the same vessel. ' +
          'Now the ring will accept things electrophilic substitution cannot deliver: Cl, Br and CN with a copper(I) salt (Sandmeyer), I with plain KI, F by heating the tetrafluoroborate (Balz–Schiemann), OH with warm water, H with hypophosphorous acid, or an azo dye with a phenol or another aniline.',
        stereoHow: null }];
    },
    byproducts: () => [{ smiles: 'O', name: 'water' }],
    gate(g, ruleId, env) {
      if (env.T > 288 && !env.TAuto) return 'Diazotisation must be run at 0–5 °C. At ' + degC(env.T) + ' °C the diazonium salt decomposes as fast as it forms — N₂ is lost, and what you isolate is the phenol (plus tar), not the salt you wanted.';
      return null;
    },
    condNotes(ruleId, env) {
      return ['NaNO₂ (1.05 equiv) added to the amine in aqueous HCl or H₂SO₄, kept at 0–5 °C with ice; the nitrous acid is generated in situ because it cannot be bought.',
        'Test for excess nitrous acid on starch–iodide paper — that is the classical endpoint, and it still works.',
        'Never isolate the dry salt: solid aryl diazonium salts are shock-sensitive. Use the cold solution as it stands.'];
    },
    stages: () => null,
  };

  const DIAZO_SUB = {
    sandmeyer_cl: { put: 'Cl', word: 'aryl chloride', reagent: 'CuCl / HCl', radical: true },
    sandmeyer_br: { put: 'Br', word: 'aryl bromide', reagent: 'CuBr / HBr', radical: true },
    sandmeyer_cn: { put: 'CN', word: 'aryl nitrile (benzonitrile type)', reagent: 'CuCN / KCN', radical: true },
    diazo_i:      { put: 'I',  word: 'aryl iodide', reagent: 'KI', radical: false },
    diazo_f:      { put: 'F',  word: 'aryl fluoride', reagent: 'HBF₄, then heat (Balz–Schiemann)', radical: false },
    diazo_oh:     { put: 'OH', word: 'phenol', reagent: 'H₂O, warm', radical: false },
    diazo_h:      { put: 'H',  word: 'the arene (the amine group removed altogether)', reagent: 'H₃PO₂ (hypophosphorous acid) or EtOH', radical: true },
  };

  F.diazoSub = {
    rules: Object.keys(DIAZO_SUB), needs: 'diazonium', structural: true,
    sites(g, ruleId) {
      return diazoniumSites(g).map(s => Object.assign(s, { score: s.aryl ? 2 : 0,
        label: 'diazonium (C' + s.c.id + '–N₂⁺)',
        reason: 'the diazonium group — N₂ is the best leaving group there is on an aromatic ring, because what leaves is a gas' }));
    },
    variants(g0, s, ruleId, env) {
      const spec = DIAZO_SUB[ruleId] || DIAZO_SUB.diazo_oh;
      const g = clone(g0);
      const cId = s.c.id;
      try { C().removeAtom(g, s.n2.id); } catch (e) {}
      try { C().removeAtom(g, s.n1.id); } catch (e) {}
      const added = [];
      if (spec.put === 'CN') {
        const c = C().addAtomNear(g, 'C', atom(g, cId));
        const n = C().addAtomNear(g, 'N', atom(g, c.id));
        const bd = bondOf(g, c.id, n.id); if (bd) bd.order = 3;
        added.push(c.id, n.id);
      } else if (spec.put === 'OH') {
        const o = C().addAtomNear(g, 'O', atom(g, cId));
        added.push(o.id);
      } else if (spec.put !== 'H') {
        const x = C().addAtomNear(g, spec.put, atom(g, cId));
        added.push(x.id);
      }
      try { C().layoutGraph(g); } catch (e) {}
      const mech = spec.radical
        ? 'This one is a RADICAL reaction, not a substitution. Copper(I) hands one electron to the diazonium ion; that breaks the C–N bond, N₂ leaves as a gas, and an ARYL RADICAL is left, which takes its halide (or cyanide) from copper(II) and regenerates copper(I). The evidence is in what it does NOT do: no rearrangement, and it works just as well on rings that would never support a cation. Hypophosphorous acid works the same way, donating a hydrogen atom to the aryl radical instead.'
        : spec.put === 'F'
          ? 'Balz–Schiemann: the diazonium tetrafluoroborate is isolated (this is the one diazonium salt that is stable enough to filter and dry) and then heated dry. It loses N₂ to give the aryl cation, which takes a fluoride from BF₄⁻ — and this is one of the very few practical routes to an aryl fluoride, because fluoride is far too poor a nucleophile for anything else.'
          : 'No copper needed here. N₂ leaves on warming to give the ARYL CATION — a genuinely high-energy species, which is why this needs heat where the Sandmeyer does not — and the nucleophile takes it from there.';
      return [{ graph: g, changed: new Set([cId].concat(added)), added, share: 1, label: spec.word,
        why: 'The diazonium group is replaced by ' + (spec.put === 'H' ? 'a hydrogen' : spec.put === 'OH' ? 'a hydroxyl' : spec.put === 'CN' ? 'a nitrile' : 'a ' + spec.put) + ', and the nitrogen leaves as N₂ gas. ' +
          'That is what makes every one of these reactions irreversible: the leaving group bubbles out of the flask. ' + mech + ' ' +
          (spec.put === 'H' ? 'Reductive deamination looks pointless until you need it: an NH₂ group is put on the ring precisely BECAUSE it directs the next substitution where you want it, and is then removed once its work is done. It is the classic way to make a 1,3,5-pattern that no direct route gives.' :
           spec.put === 'I' ? 'Aryl iodides are made this way and essentially no other way from a simple ring — and the iodide is then the best possible partner for a palladium coupling.' :
           spec.put === 'OH' ? 'Note that this is also the DECOMPOSITION reaction of a diazonium salt: it is what happens by accident whenever the ice bath is allowed to warm up.' : ''),
        stereoHow: null }];
    },
    byproducts: (g, s, ruleId) => [{ smiles: 'N#N', name: 'N₂ — leaves as a gas, which is what makes this irreversible' }],
    gate(g, ruleId, env) {
      if (ruleId === 'diazo_oh' && env.T < 313 && !env.TAuto) return 'The phenol needs WARMING (50–100 °C): the aryl cation only forms when the salt is made to decompose. Cold, the diazonium salt simply sits there.';
      return null;
    },
    condNotes(ruleId, env) {
      const spec = DIAZO_SUB[ruleId] || {};
      return [spec.reagent + ' — added to the cold diazonium solution' + (ruleId === 'diazo_oh' || ruleId === 'diazo_f' ? ', then warmed' : ''),
        'Nitrogen is evolved; the gas coming off is how you know it is working.'];
    },
    stages: () => null,
  };

  const FAMILY_OF = {};
  Object.keys(F).forEach(k => F[k].rules.forEach(r => { FAMILY_OF[r] = F[k]; }));

  /* what the reagent ignores — the chemo-selectivity sentence */
  const IGNORE = {
    alkene: { addHX: null, addRad: null, addX2: null, epox: null, diol: null, hydrobor: null, h2: null,
      oxid: 'PCC has no way to touch a C=C — it needs an O–H and a carbinol C–H to strip', pbr3: 'PBr₃ only reacts with O–H groups', dehyd: 'the acid can protonate a C=C too, but with nothing to add it just returns', halide: 'an alkoxide/cyanide is not electrophilic, so it has nothing to give a π bond', reduce: 'a hydride has nothing to attack in an electron-rich C=C — only polar C=O bonds' },
    alcohol: { addHX: 'H–X can also protonate an OH, but the alkene is the far better nucleophile and adds first', addX2: 'Br₂ does not react with an O–H', epox: 'peracids leave O–H alone', diol: 'OsO₄ ignores an existing OH', hydrobor: 'borane is quenched by O–H only slowly — the alkene wins', h2: 'H₂/Pd does not touch an alcohol', halide: 'the alcohol’s OH⁻ is a hopeless leaving group, so the base cannot substitute or eliminate it', reduce: 'an alcohol is already at the reduced level; hydride does nothing to it (NaBH₄ merely tolerates the O–H)', addRad: 'Br• takes an allylic/benzylic C–H long before an O–H' },
    halide: { addHX: 'the C–X stays — HX adds to π bonds only', addX2: null, oxid: 'a C–X survives the oxidant', reduce: 'NaBH₄ does not reduce alkyl halides (LiAlH₄ can, slowly)', dehyd: 'the halide is not displaced under these conditions', pbr3: null, epox: null, diol: null, hydrobor: null, h2: 'H₂/Pd can cleave a benzylic C–X, but an ordinary alkyl halide survives' },
    carbonyl: { addHX: 'HX does not add across a C=O (the reverse is favoured)', addX2: 'Br₂ does not add to a C=O (though it can α-brominate a ketone via the enol, slowly)', oxid: 'a ketone/ester is already oxidised and stops here; an aldehyde WOULD be oxidised further by KMnO₄', pbr3: 'PBr₃ ignores a C=O', dehyd: 'the acid protonates the carbonyl reversibly, but nothing leaves', halide: 'hydroxide/alkoxide would add to an aldehyde reversibly (hydrate/hemiacetal), but the substitution at C–X is what goes to completion', epox: 'mCPBA leaves an aldehyde/ester alone at 0 °C (a ketone can undergo Baeyer–Villiger if pushed)', diol: 'OsO₄ ignores C=O', hydrobor: 'BH₃ reduces a C=O only slowly; the alkene reacts first', h2: 'a C=O survives Pd/H₂ at 1 atm', addRad: 'a C=O is not attacked by Br•' },
    alkyne: { addHX: 'an alkyne adds HX too, but more slowly than the alkene (a vinyl cation is worse than an alkyl one)', addX2: 'Br₂ adds to alkynes more slowly than to alkenes', epox: 'peracids do not epoxidise alkynes', diol: 'OsO₄ is slow with alkynes', hydrobor: 'BH₃ hydroborates alkynes too (a second, slower site)', oxid: 'the oxidant leaves a C≡C alone (hot KMnO₄ would cleave it)', halide: null, reduce: null, pbr3: null, dehyd: null, addRad: null },
  };
  const NEED_WORD = { any: 'a matching pair of functional groups (an acid and an amine, a halide and a nucleophile, a diene and a dienophile …)', arene: 'an aromatic ring with a free C–H', ester: 'an ester', alkene: 'a C=C', 'alkene-or-alkyne': 'a C=C or C≡C', alcohol: 'an alcohol O–H on an sp³ carbon', halide: 'a C–X leaving group', carbonyl: 'a C=O', benzylic: 'a benzylic CH₃ or CH₂ (an alkyl carbon attached to an aromatic ring)' };

  /* ---------- the run ---------- */
  function run(g0, ruleId, envIn) {
    const rule = C().RULES.find(r => r.id === ruleId);
    let fam = FAMILY_OF[ruleId] || (String(ruleId).indexOf('mol:') === 0 ? F.bimol : null);
    if (!rule || !fam || !g0 || !g0.atoms || !g0.atoms.length) return null;
    const env = resolveEnv(ruleId, envIn);
    const res = { rule, ruleId, env, sub: g0, sites: [], outcomes: [], byproducts: [], notes: [], summary: {} };
    let gate = fam.gate ? fam.gate(g0, ruleId, env) : null;
    let sites = [];
    try { sites = fam.sites(g0, ruleId, env) || []; } catch (e) { sites = []; }
    let live = sites.filter(s => !s.inert);
    /* no live site for the reagent's own job? some reagents have a second job (aqueous HX → ester hydrolysis) */
    if (!live.length && fam.fallback) {
      const fbs = [].concat(fam.fallback).map(k => F[k]).filter(Boolean);
      for (const fb of fbs) {
        let s2 = [];
        try { s2 = fb.sites(g0, ruleId, env) || []; } catch (e) { s2 = []; }
        const l2 = s2.filter(s => !s.inert);
        if (l2.length) {
          const why = sites.length ? sites.filter(s => s.inert).map(s => s.label + ': ' + s.inert).join('; ') : 'it has no ' + (NEED_WORD[fam.needs] || fam.needs).replace(/^an? /, '');
          res.fallbackNote = rule.label + ' cannot do its usual job here (' + why + '), so it acts on the ' + fb.needs + ' instead.';
          fam = fb; sites = s2; live = l2; gate = fam.gate ? fam.gate(g0, ruleId, env) : null;
          break;
        }
      }
    }
    res.sites = sites;
    /* --- chemo sentence: which kinds of group are present, which one the reagent wants --- */
    const present = new Set(C().findGroups(g0).map(x => x.type));
    const famKey = Object.keys(F).find(k => F[k] === fam);
    const ignored = [...present].filter(t => !(fam.needs === t || (fam.needs === 'alkene-or-alkyne' && (t === 'alkene' || t === 'alkyne'))))
      .map(t => { const r = IGNORE[t] && IGNORE[t][famKey]; return { type: t, why: r === undefined ? null : r }; });
    res.ignored = ignored;
    if (gate) { res.error = gate; res.summary.conditions = gate; return res; }
    if (!live.length) {
      const inertWhy = sites.filter(s => s.inert).map(s => s.label + ': ' + s.inert);
      res.error = sites.length ? 'Nothing here that ' + rule.label + ' can attack — ' + inertWhy.join('; ') + '.'
        : 'This reagent needs ' + (NEED_WORD[fam.needs] || fam.needs) + ', and this molecule has none' + (ignored.length ? ' (it has ' + ignored.map(i => i.type).join(', ') + ', which ' + rule.label + ' leaves alone)' : '') + '.';
      return res;
    }
    /* --- site shares --- */
    const sc = live.map(s => s.score || 0);
    const fr = share(sc);
    live.forEach((s, i) => { s.share = fr[i]; });
    live.sort((p, q) => q.share - p.share);
    res.chosen = live[0];
    /* --- outcomes per site --- */
    const all = [];
    live.forEach(s => {
      let vs = [];
      try { vs = fam.variants(g0, s, ruleId, env) || []; } catch (e) { vs = []; }
      vs.forEach(v => { v.site = s; v.share = (v.share === undefined ? 1 : v.share) * s.share; all.push(v); });
    });
    if (fam.crossReact) { try { fam.crossReact(g0, ruleId, env).forEach(v => { v.site = null; all.forEach(o => { o.share *= (1 - v.share); }); all.push(v); }); } catch (e) {} }
    /* --- excess reagent: every live site reacts (main variant at each) --- */
    if (env.excess && fam.excess) { try { fam.excess(g0, ruleId, env, all, live); } catch (e) {} }
    if (env.excess && live.length > 1 && !fam.excess) {
      try {
        let g = g0; const changed = new Set(); let ok = true;
        for (const s of live) {
          const cur = fam.sites(g, ruleId, env).filter(x => !x.inert);
          const sIds = s.ringIds || s.atoms;
          const same = cur.find(x => x.atoms[0] === s.atoms[0] && x.atoms[1] === s.atoms[1]) || cur.find(x => (x.ringIds || x.atoms).some(id => sIds.includes(id)));
          if (!same) { ok = false; break; }
          const vs = fam.variants(g, same, ruleId, env) || [];
          if (!vs.length) { ok = false; break; }
          vs.sort((p, q) => (q.share || 0) - (p.share || 0));
          g = vs[0].graph; vs[0].changed.forEach(id => changed.add(id));
        }
        if (ok && g !== g0) {
          all.forEach(o => { o.share *= 0.12; });
          all.unshift({ graph: g, changed, added: [], share: 0.88, site: live[0], exhaustive: true, label: 'all ' + live.length + ' sites reacted (excess reagent)',
            why: 'With the reagent in excess nothing has to compete: every reactive site is consumed. The single-site products below are what you would isolate if the reagent ran out first.', stereoHow: null });
        }
      } catch (e) {}
    }
    /* --- merge identical products, sort, tag --- */
    const seen = new Map();
    all.forEach(o => { const k = keyOf(o.graph); if (seen.has(k)) { seen.get(k).share += o.share; } else seen.set(k, o); });
    let outs = [...seen.values()].filter(o => o.share > 0.002);
    const tot = outs.reduce((a, o) => a + o.share, 0) || 1;
    outs.forEach(o => { o.share /= tot; });
    outs.sort((p, q) => q.share - p.share);
    outs.forEach((o, i) => {
      o.tag = i === 0 ? 'major' : o.share >= 0.05 ? 'minor' : 'trace';
      o.name = nameOf(o.graph);
      if (o.nameHint) o.name = o.nameHint;
      try { o.formula = C().formula(o.graph); } catch (e) { o.formula = ''; }
      try { o.smiles = C().toSmiles(o.graph); } catch (e) { o.smiles = ''; }
      /* a wedge on a carbon that is not a stereocentre says nothing (the CH2–OH
         of a syn-diol, say) and only confuses the picture — drop it */
      try { o.graph.bonds.forEach(b => { if (b.stereo && b.narrow !== undefined) { const c = atom(o.graph, b.narrow); if (!c || !C().isStereocenter(o.graph, c)) { delete b.stereo; delete b.narrow; } } }); } catch (e) {}
      try { o.stereo = enantioForms(g0, o.graph, o.changed || new Set(), o.stereoHow); } catch (e) { o.stereo = null; }
      /* a racemate is named with RS/SR descriptors, not as one enantiomer */
      if (o.stereo && o.stereo.kind === 'racemic' && /\(\d+[a-z]?[RS](,\d+[a-z]?[RS])*\)/.test(o.name)) o.name = o.name.replace(/\((\d+[a-z]?[RS](?:,\d+[a-z]?[RS])*)\)/, (m, inner) => '(' + inner.replace(/(\d+[a-z]?)([RS])/g, (x, loc, d) => loc + (d === 'R' ? 'RS' : 'SR')) + ')');
      try { o.stages = fam.stages ? fam.stages(g0, o.site || live[0], o, ruleId, env) : null; } catch (e) { o.stages = null; }
    });
    res.outcomes = outs;
    res.major = outs[0];
    try { res.byproducts = fam.byproducts ? fam.byproducts(ruleId, env, live[0]) : []; } catch (e) { res.byproducts = []; }
    try { res.notes = fam.condNotes ? fam.condNotes(ruleId, env) : []; } catch (e) { res.notes = []; }
    if (res.fallbackNote) res.notes.unshift(res.fallbackNote);
    /* --- the summary sentences --- */
    const S = res.summary;
    const grpWord = { any: 'matching functional group', arene: 'aromatic ring', alkene: 'C=C', alcohol: 'alcohol', halide: 'alkyl halide', carbonyl: 'carbonyl', alkyne: 'C≡C', benzylic: 'benzylic C–H' };
    const kindsHere = [...present].map(t => grpWord[t] || t);
    const igText = ignored.filter(i => i.why).map(i => (grpWord[i.type] || i.type) + ' — ' + i.why);
    S.chemo = (kindsHere.length > 1 ? 'The molecule offers ' + kindsHere.join(', ') + '. ' : '') +
      rule.label + ' goes for the ' + (grpWord[fam.needs] || fam.needs) + (fam.needs === 'alkene-or-alkyne' ? ' (C=C / C≡C)' : '') +
      (igText.length ? '; the rest is left alone: ' + igText.join('; ') + '.' : '.');
    if (live.length > 1) {
      const top = live[0];
      S.site = 'There are ' + live.length + ' ' + (grpWord[fam.needs === 'alkene-or-alkyne' ? 'alkene' : fam.needs] || fam.needs) + ' sites' +
        (sites.length > live.length ? ' that can react (' + (sites.length - live.length) + ' more cannot: ' + sites.filter(s => s.inert).map(s => s.label + ' — ' + s.inert).join('; ') + ')' : '') +
        '. ' + top.label + ' reacts first (' + pct(top.share) + '%): ' + top.reason + '. ' +
        live.slice(1).map(s => s.label + ' (' + pct(s.share) + '%): ' + s.reason).join('; ') + '.' +
        (env.excess ? ' With excess reagent all of them react.' : ' With one equivalent, that ordering is what decides the product mix.');
    } else if (sites.length > 1) {
      S.site = 'Of the ' + sites.length + ' ' + (grpWord[fam.needs === 'alkene-or-alkyne' ? 'alkene' : fam.needs] || fam.needs) + ' sites only ' + live[0].label + ' can react: ' +
        sites.filter(s => s.inert).map(s => s.label + ' is ' + s.inert).join('; ') + '.';
    } else {
      S.site = 'One reactive site, ' + live[0].label + ': ' + live[0].reason + '.';
    }
    S.regio = outs.length > 1 ? outs.slice(0, 3).map(o => o.label + ' ' + pct(o.share) + '%').join(' · ') + '. ' + (outs[0].why || '') : (outs[0].why || '');
    const st = outs[0].stereo;
    S.stereo = st ? (st.kind === 'racemic' ? 'Racemic: ' : st.kind === 'meso' ? 'Meso: ' : st.kind === 'diastereomers' ? 'Diastereomers: ' : '') + st.note : (outs[0].stereoHow || '');
    S.conditions = [
      'Conditions used: ' + degC(env.T) + ' °C' + (env.TAuto ? ' (the reagent’s usual)' : '') + ', pH ' + (env.pHAuto ? env.pHVal + ' (set by the reagent itself)' : env.pHVal.toFixed(1) + ' (yours)') + ', ' + env.P + ' atm, ' + (SOLVENT_NAMES[env.solventVal] || env.solventVal) + (env.light ? ', hν on' : '') + ', ' + (env.excess ? 'excess reagent' : '1 equivalent') + '.'
    ].concat(res.notes).join(' ');
    return res;
  }

  /* the reagents that DO something to this molecule, for the React menu */
  function reactiveReagents(g, env) {
    const out = [];
    C().RULES.forEach(r => {
      let res = null; try { res = run(g, r.id, env); } catch (e) { res = null; }
      if (res && !res.error && res.major) out.push({ rule: r, res });
    });
    return out;
  }

  return { run, reactiveReagents, resolveEnv, DEFAULT_ENV, REGIME, SOLVENT_NAMES, families: F, familyOf: FAMILY_OF, mirror, pct, degC };
})();

