/* RetroExt — disconnections added to the Retrosynthesis planner's list (Chem.RETRO_TRANSFORMS) without editing
   js/chem.js, which is lifted verbatim from organic-chem-tool (as js/rxn_ext.js adds forward reaction families).
   Load it right after js/chem.js, so every page that plans (Retrosynthesis, Possibilities) sees the same list.

   Added:
   - Alkyne ← vicinal dibromide (double dehydrohalogenation). Undoing it puts one Br on each alkyne carbon and makes
     the C≡C single; the planner's own "Vicinal dibromide ← alkene" then takes the dibromide back to the alkene. That
     gives the textbook two-step route from an alkene: Br₂ adds anti, then two E2 eliminations of HBr with a strong
     base, e.g. trans-stilbene → meso-1,2-dibromo-1,2-diphenylethane → diphenylacetylene (KOH in triethylene glycol at
     about 190 °C, the usual teaching-lab route). Either diastereomer of the dibromide (meso or racemic) gives the same
     alkyne, so the dibromide is drawn without stereochemistry. A terminal alkyne needs a third equivalent of NaNH₂:
     the alkyne C–H is deprotonated, and the work-up (H₃O⁺) gives the alkyne back; the reagent text says so.
   - Directing effects for nitration, sulfonation and aromatic halogenation (see below the additions).
   - Indole ← indoline (dehydrogenation). Undoing it makes the indole's C2=C3 bond single. Offered only when that
     indoline can be bought: indoles are mostly built by ring-forming reactions (Fischer and others), and oxidising an
     indoline is the way when the indoline is on the shelf. Without that limit every indole-containing target
     (tryptophan, say) would be offered an indoline nobody sells. Reagents are the ones most recorded for it in the
     app's reaction knowledge base: MnO₂ in CH₂Cl₂ (5 of the 12 closest records), DDQ or chloranil. */
(() => {
  'use strict';
  const C = window.Chem;
  if (!C || !Array.isArray(C.RETRO_TRANSFORMS)) return;
  const byId = (g, id) => g.atoms.find(a => a.id === id);
  const bondsOf = (g, id) => g.bonds.filter(b => b.a === id || b.b === id);
  const hCount = (g, a) => {                      // implicit H on a neutral carbon: 4 − bond orders
    const used = bondsOf(g, a.id).reduce((s, b) => s + (b.order || 1), 0);
    return a.element === 'C' && !a.charge ? Math.max(0, 4 - used) : 0;
  };

  const ALKYNE_FROM_DIBROMIDE = {
    name: 'Alkyne ← vicinal dibromide (double E2)',
    find(g) {
      for (const b of g.bonds) {
        if (b.order !== 3) continue;
        const A = byId(g, b.a), B = byId(g, b.b);
        if (!A || !B || A.element !== 'C' || B.element !== 'C' || A.charge || B.charge) continue;
        // not in a ring (a C≡C in a small ring is not made this way) and not a cumulene end
        if (bondsOf(g, A.id).some(x => x !== b && x.order > 1) || bondsOf(g, B.id).some(x => x !== b && x.order > 1)) continue;
        const terminal = hCount(g, A) > 0 || hCount(g, B) > 0;
        const p = C.clone(g);
        const pb = p.bonds.find(x => (x.a === A.id && x.b === B.id) || (x.a === B.id && x.b === A.id));
        pb.order = 1;
        for (const at of [A, B]) {                  // one Br on each carbon, drawn beside it
          const id = p.nextId++;
          p.atoms.push({ id, element: 'Br', charge: 0, x: at.x, y: at.y - 48 });
          p.bonds.push({ a: at.id, b: id, order: 1 });
        }
        return { precursors: [p],
          reagent: terminal
            ? 'NaNH₂ (3 equiv) in liquid NH₃, then H₃O⁺ work-up — two E2 eliminations of HBr; the third equivalent deprotonates the terminal alkyne, and the work-up gives it back'
            : 'KOH in triethylene glycol at about 190 °C, or NaNH₂ in liquid NH₃ — two successive E2 eliminations of HBr',
          note: 'A C≡C comes from two HBr eliminations out of a vicinal dibromide (the first gives a vinyl bromide, the second, harder one the alkyne). The dibromide comes from the alkene plus Br₂; either diastereomer (meso or racemic) gives the same alkyne.' };
      }
      return null;
    },
  };

  const order = (g, x, y) => { const b = g.bonds.find(b => (b.a === x && b.b === y) || (b.a === y && b.b === x)); return b ? (b.order || 1) : 0; };
  const nbrs = (g, id) => bondsOf(g, id).map(b => (b.a === id ? b.b : b.a));
  /* a Kekulé benzene ring through the bond c–d: a 6-cycle of carbons with three double bonds */
  function benzeneThrough(g, c, d) {
    const walk = (path) => {
      const last = path[path.length - 1];
      if (path.length === 6) return order(g, last, path[0]) ? path : null;
      for (const n of nbrs(g, last)) {
        if (path.includes(n) || byId(g, n).element !== 'C') continue;
        const r = walk(path.concat(n)); if (r) return r;
      }
      return null;
    };
    const ring = walk([d, c]);
    if (!ring) return false;
    let dbl = 0; for (let i = 0; i < 6; i++) if (order(g, ring[i], ring[(i + 1) % 6]) === 2) dbl++;
    return dbl === 3;
  }
  const buyable = g => { try { if (C.isBasicMaterial(g)) return true; } catch (e) {} try { return C.buyableMatches(g, 2).some(x => x.exact); } catch (e) { return false; } };

  const INDOLE_FROM_INDOLINE = {
    name: 'Indole ← indoline oxidation (dehydrogenation)',
    find(g) {
      for (const N of g.atoms) {
        if (N.element !== 'N' || N.charge || nbrs(g, N.id).length > 3) continue;
        // the pyrrole ring read N–C2=C3–C3a…C7a–N, with C3a–C7a the bond shared with a benzene ring
        for (const a of nbrs(g, N.id)) for (const b of nbrs(g, a)) {
          if (b === N.id || order(g, N.id, a) !== 1 || order(g, a, b) !== 2) continue;
          const A = byId(g, a), B = byId(g, b);
          if (A.element !== 'C' || B.element !== 'C' || A.charge || B.charge) continue;
          for (const c of nbrs(g, b)) {
            if (c === a || order(g, b, c) !== 1 || byId(g, c).element !== 'C') continue;
            const d = nbrs(g, c).find(x => x !== b && x !== N.id && order(g, x, N.id) === 1 && byId(g, x).element === 'C');
            if (d === undefined || !benzeneThrough(g, c, d)) continue;
            const p = C.clone(g);
            p.bonds.find(x => (x.a === a && x.b === b) || (x.a === b && x.b === a)).order = 1;
            if (!buyable(p)) return null;
            return { precursors: [p],
              reagent: 'MnO₂ (excess) in CH₂Cl₂ at room temperature, or DDQ or chloranil in dioxane or toluene — removes 2 H from C2–C3',
              note: 'Indolines lose two hydrogens easily because the product is aromatic (the pyrrole ring joins the 10 π-electron indole system). Activated MnO₂, DDQ and chloranil are the oxidants most recorded for it; Pd/C on heating also works.' };
          }
        }
      }
      return null;
    },
  };

  for (const t of [ALKYNE_FROM_DIBROMIDE, INDOLE_FROM_INDOLINE]) if (!C.RETRO_TRANSFORMS.some(x => x.name === t.name)) C.RETRO_TRANSFORMS.push(t);

  /* ---- directing effects for the electrophilic aromatic substitutions in Chem's list ----
     'Nitration undone', 'Sulfonation undone' and 'Aryl halide undone' (js/chem.js) take the group off whatever ring
     position it sits on, so they offered ethyl 4-nitrobenzoate ← ethyl benzoate + HNO₃/H₂SO₄, which gives mostly the
     meta isomer (an ester is a meta director). Each is wrapped here: the step is kept only when the group's position
     fits the substituents already on the precursor's ring, by the textbook rules —
       an ortho/para director present (–OH, –OR, –OCOR, –NH₂/–NHR/–NHCOR, alkyl, aryl, vinyl, halogen, –SR): the new
         group must be ortho or para to at least one of them (activators and halogens outrank meta directors);
       otherwise meta directors only (C=O, C≡N, C=N, CF₃, NO₂, sulfonyl, N⁺): it must be meta to every one;
       an unsubstituted ring: any position.
     Ortho/para mixtures are not split further. Rings other than benzene are left as they were. */
  function benzeneRing(g, id) {                      // the six ring carbons in order, starting at id, or null
    const walk = path => {
      const last = path[path.length - 1];
      if (path.length === 6) {
        if (!order(g, last, path[0])) return null;
        let dbl = 0; for (let i = 0; i < 6; i++) if (order(g, path[i], path[(i + 1) % 6]) === 2) dbl++;
        return dbl === 3 ? path : null;
      }
      for (const n of nbrs(g, last)) {
        if (path.includes(n) || byId(g, n).element !== 'C') continue;
        const r = walk(path.concat(n)); if (r) return r;
      }
      return null;
    };
    return walk([id]);
  }
  function directs(g, ringAtom, x) {                 // 'op', 'm' or null for the substituent atom x on ringAtom
    const X = byId(g, x), el = X.element;
    if (['F', 'Cl', 'Br', 'I', 'O'].includes(el)) return 'op';
    if (el === 'N') return X.charge > 0 ? 'm' : 'op';          // nitro and ammonium N carry the + charge in these graphs
    if (el === 'S') return nbrs(g, x).filter(n => byId(g, n).element === 'O').length >= 2 ? 'm' : 'op';
    if (el === 'C') {
      const multiple = bondsOf(g, x).some(b => {
        const o = b.a === x ? b.b : b.a; if (o === ringAtom) return false;
        const e = byId(g, o).element;
        return (b.order === 2 && (e === 'O' || e === 'N')) || (b.order === 3 && e === 'N');
      });
      if (multiple || nbrs(g, x).filter(n => byId(g, n).element === 'F').length >= 3) return 'm';
      return 'op';                                    // alkyl, aryl, vinyl
    }
    return null;
  }
  /* a free amino group (–NH₂, –NHR, –NR₂; not an amide or sulfonamide N) on the ring: in HNO₃/H₂SO₄ or fuming H₂SO₄ it
     is protonated (–NH₃⁺ directs meta) and the ring is oxidised, and with Br₂ or Cl₂ the ring is substituted three
     times over (aniline → 2,4,6-tribromoaniline). The textbook way is to protect it as the acetanilide first, so the
     step is not offered on the free amine. */
  const freeAmine = (g, x) => byId(g, x).element === 'N' && !byId(g, x).charge &&
    !nbrs(g, x).some(n => bondsOf(g, n).some(b => b.order === 2 && ['O', 'S'].includes(byId(g, b.a === n ? b.b : b.a).element)));
  function positionFits(p, attach) {
    const ring = benzeneRing(p, attach); if (!ring) return true;
    const inRing = new Set(ring), op = [], m = [];
    if (ring.some(r => nbrs(p, r).some(x => !inRing.has(x) && freeAmine(p, x)))) return false;
    ring.forEach((r, i) => {
      if (r === attach) return;
      const dist = Math.min(i, 6 - i);                // 1 ortho, 2 meta, 3 para
      for (const x of nbrs(p, r)) {
        if (inRing.has(x)) continue;
        const d = directs(p, r, x);
        if (d === 'op') op.push(dist); else if (d === 'm') m.push(dist);
      }
    });
    if (op.length) return op.some(d => d === 1 || d === 3);
    if (m.length) return m.every(d => d === 2);
    return true;
  }
  const EAS = new Set(['Nitration undone', 'Sulfonation undone', 'Aryl halide undone']);
  for (const t of C.RETRO_TRANSFORMS) {
    if (!EAS.has(t.name) || t._directing) continue;
    const inner = t.find;
    t.find = function (g) {
      const r = inner.call(this, g);
      if (!r || !r.precursors || r.precursors.length !== 1) return r;
      const p = r.precursors[0], kept = new Set(p.atoms.map(a => a.id));
      const cut = g.bonds.find(b => kept.has(b.a) !== kept.has(b.b));     // the bond from the ring to the removed group
      if (!cut) return r;
      return positionFits(p, kept.has(cut.a) ? cut.a : cut.b) ? r : null;
    };
    t._directing = true;
  }

  /* ---- phenyl esters: phenol + acid ANHYDRIDE, not Fischer esterification ----
     chem.js's 'Ester disconnection' gives every ester as acid + alcohol by Fischer esterification. For an aryl ester
     (the ester O on a benzene ring) that is wrong: a phenol's O is a poor nucleophile (its lone pair is shared with the
     ring) and the equilibrium lies on the side of the acid and phenol, so aspirin is made from salicylic acid and
     acetic ANHYDRIDE with a few drops of H₂SO₄ or H₃PO₄ (or an acyl chloride with pyridine), not from acetic acid.
     The wrapper keeps the disconnection and turns the acid piece into its symmetric anhydride. */
  function anhydrideOf(acid, cId) {
    const oh = nbrs(acid, cId).find(x => byId(acid, x).element === 'O' && order(acid, cId, x) === 1 && nbrs(acid, x).length === 1);
    if (oh === undefined) return null;
    const a = C.clone(acid), map = new Map();
    for (const at of acid.atoms) {
      if (at.id === oh) continue;
      const id = a.nextId++; map.set(at.id, id);
      a.atoms.push(Object.assign({}, at, { id, x: (at.x || 0) + 60, y: at.y || 0 }));
    }
    for (const b of acid.bonds) if (b.a !== oh && b.b !== oh) a.bonds.push(Object.assign({}, b, { a: map.get(b.a), b: map.get(b.b) }));
    a.bonds.push({ a: oh, b: map.get(cId), order: 1 });      // R–C(=O)–O–C(=O)–R
    return a;
  }
  const ESTER = C.RETRO_TRANSFORMS.find(t => t.name === 'Ester disconnection');
  if (ESTER && !ESTER._aryl) {
    const inner = ESTER.find;
    ESTER.find = function (g) {
      const r = inner.call(this, g);
      if (!r || !r.precursors || r.precursors.length !== 2) return r;
      const [alcohol, acid] = r.precursors;
      /* the acyl carbon: in the acid piece, a C with =O and an –OH that has no other neighbour */
      const acyl = acid.atoms.find(x => x.element === 'C' && nbrs(acid, x.id).some(n => byId(acid, n).element === 'O' && order(acid, x.id, n) === 2) &&
        nbrs(acid, x.id).some(n => byId(acid, n).element === 'O' && order(acid, x.id, n) === 1 && nbrs(acid, n).length === 1));
      /* the ester O in the alcohol piece is now an –OH; is its carbon aromatic? */
      const phenolic = alcohol.atoms.some(o => o.element === 'O' && nbrs(alcohol, o.id).length === 1 &&
        (() => { const c = nbrs(alcohol, o.id)[0]; try { return byId(alcohol, c).element === 'C' && C.isAromaticCarbon(alcohol, c); } catch (e) { return false; } })() &&
        g.atoms.some(x => x.id === o.id) && nbrs(g, o.id).length === 2);
      if (!acyl || !phenolic) return r;
      const anh = anhydrideOf(acid, acyl.id);
      if (!anh) return r;
      return { precursors: [alcohol, anh],
        reagent: 'anhydride, H₂SO₄ or H₃PO₄ (a few drops), warm — or the acyl chloride with pyridine; a phenol is too weak a nucleophile for Fischer esterification',
        note: 'An aryl ester breaks back at the acyl C–O bond to the phenol and an acylating agent. The phenol O is a poor nucleophile (its lone pair is shared with the ring), so the acid itself will not do: the anhydride (or acyl chloride) acylates it, as in aspirin from salicylic acid and acetic anhydride.' };
    };
    ESTER._aryl = true;
  }

  /* `version` changes whenever what these additions return changes (js/possibilities.js keys its cache on it) */
  window.RetroExt = { transforms: [ALKYNE_FROM_DIBROMIDE, INDOLE_FROM_INDOLINE], version: 'eas2-aryl1', positionFits };
})();
