// "Open Mechanism in Draw", the Draw half, in a headless browser: Draw.showMechanism lays every workbook mechanism (the 52
// rules with a Mechanism sheet) out on the Draw canvas, and an engine-form step (SN2 of bromoethane by HO⁻) too.
//   python tools/serve.py   (another terminal)        node tools/browser_mechanism.js [http://localhost:8000/app.html] [screenshot dir]
// For every mechanism: no exception; every arrow token of the Mechanism sheet drawn (red curved arrows = tokens, fishhooks =
// the '~' tokens); every atom carries exactly its lone pairs and unpaired electrons as dots (plus the ones drawn at an
// arrow's tail on an atom that shows none, a metal), none on a bond or a label, each pair clear of the others; no map
// numbers; the step rows do not overlap. The arrows' ends: an arrow making a new bond between two far-apart atoms ends at
// the atom attacked, not in the empty middle; no two heads on one spot; an arrow from a bond starts on its middle; no
// head lands on another arrow's tail (or the electrons it leaves). Also the new arrow kinds (curved, fishhook) behave
// like the other arrows.
// Screenshots of a few mechanisms are saved (default: the system temp folder). The exit status is 1 when a check fails.
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const base = process.argv[2] || 'http://localhost:8000/app.html';
const shots = process.argv[3] || require('os').tmpdir();
let fails = 0, n = 0;
const check = (name, ok, detail) => { n++; if (!ok) { fails++; console.log('FAIL', name, detail === undefined ? '' : '\n    ' + JSON.stringify(detail).slice(0, 600)); } else console.log('ok  ', name); };
(async () => {
  const browser = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const ctx = await browser.newContext({ viewport: { width: 1500, height: 1000 }, acceptDownloads: true });
  const page = await ctx.newPage();
  const errors = []; page.on('pageerror', e => errors.push(e.message)); page.on('console', m => { if (m.type() === 'error') errors.push(m.text()); });
  await page.goto(base + '#draw', { waitUntil: 'load' });
  await page.waitForFunction(() => window.Draw && window.REACTION_RULES && window.Chem, null, { timeout: 30000 });
  check('Draw.showMechanism is a function', await page.evaluate(() => typeof window.Draw.showMechanism === 'function'));

  /* the page-side checks of one drawn mechanism (Draw.lastMechanism): returns {problems: [...], ...counts} */
  await page.evaluate(() => {
    window.__mechCheck = (expect) => {
      const D = window.Draw, C = window.Chem, L = D.lastMechanism, problems = [];
      const items = L.items, boxes = items.filter(i => i instanceof D.MoleculeBox), els = items.filter(i => i instanceof D.ElectronItem);
      const arrows = items.filter(i => i instanceof D.ArrowItem && (i.kind === 'curved' || i.kind === 'fishhook'));
      if (expect) {
        if (L.unresolved.length) problems.push('unresolved arrow tokens: ' + L.unresolved.join(' | '));
        if (L.tokens !== expect.tokens || L.arrows !== expect.tokens) problems.push(`arrows: ${L.arrows} drawn of ${L.tokens} tokens (sheet has ${expect.tokens})`);
        if (arrows.length !== expect.tokens) problems.push(`curved/fishhook items ${arrows.length} != ${expect.tokens} tokens`);
        if (arrows.filter(a => a.kind === 'fishhook').length !== expect.fish) problems.push(`fishhooks ${arrows.filter(a => a.kind === 'fishhook').length} != ${expect.fish}`);
      }
      if (arrows.some(a => a.color !== '#d62828')) problems.push('a mechanism arrow is not red');
      // electrons: each dot belongs to its nearest atom; every atom has exactly its lone pairs and unpaired electrons
      const atoms = boxes.flatMap(b => b.lay.atoms.map((a, i) => ({ b, a, i, P: { x: b.x + b.atomXY[i].x, y: b.y + b.atomXY[i].y }, lab: b.atomLabels[i], pairs: 0, singles: 0 })));
      const owner = new Map();
      for (const e of els) {
        let best = null; for (const t of atoms) { const d = Math.hypot(t.P.x - e.x, t.P.y - e.y); if (!best || d < best.d) best = { t, d }; }
        if (!best || best.d > 26) { problems.push(`a dot far from every atom (${Math.round(e.x)}, ${Math.round(e.y)})`); continue; }
        if (e.pair) best.t.pairs++; else best.t.singles++;
        owner.set(e, best.t);
      }
      let lpAtoms = 0;
      const extra = (t, pair) => (L.extras || []).filter(x => x.box === t.b && x.id === t.a.id && x.pair === pair).length;
      for (const t of atoms) {
        const want = C.lonePairCount(t.b.lay, t.a) + extra(t, true), rad = (t.a.radical === true ? 1 : Math.round(+t.a.radical || 0)) + extra(t, false);
        if (want && /^(O|N|S|F|Cl|Br|I)$/.test(t.a.element)) lpAtoms++;
        if (t.pairs !== want || t.singles !== rad) problems.push(`${t.a.element}${t.a.charge ? (t.a.charge > 0 ? '+' : '-') : ''} in ${t.b.smiles}: ${t.pairs} pairs/${t.singles} single, wants ${want}/${rad}`);
      }
      // no dot on a label or on a bond line (the bond as drawn: it stops at a label)
      const rectDist = (x, y, r) => Math.hypot(Math.max(r.x0 - x, 0, x - r.x1), Math.max(r.y0 - y, 0, y - r.y1));
      const segDist = (p, a, b) => { const vx = b.x - a.x, vy = b.y - a.y, l2 = vx * vx + vy * vy || 1, t = Math.max(0, Math.min(1, ((p.x - a.x) * vx + (p.y - a.y) * vy) / l2)); return Math.hypot(p.x - a.x - vx * t, p.y - a.y - vy * t); };
      const exitDist = (rects, ux, uy) => { let best = 0; for (const r of rects) { let t0 = -1e9, t1 = 1e9, hit = true; for (const [o0, o1, u] of [[r.x0, r.x1, ux], [r.y0, r.y1, uy]]) { if (Math.abs(u) < 1e-9) { if (o0 > 0 || o1 < 0) { hit = false; break; } continue; } let a = o0 / u, c = o1 / u; if (a > c) [a, c] = [c, a]; t0 = Math.max(t0, a); t1 = Math.min(t1, c); } if (hit && t1 >= Math.max(t0, 0)) best = Math.max(best, t1); } return best; };
      const segs = boxes.flatMap(b => b.bondPairs.map(([i, j]) => {
        const p = { x: b.x + b.atomXY[i].x, y: b.y + b.atomXY[i].y }, q = { x: b.x + b.atomXY[j].x, y: b.y + b.atomXY[j].y }, L = Math.hypot(q.x - p.x, q.y - p.y) || 1, ux = (q.x - p.x) / L, uy = (q.y - p.y) / L;
        const tp = b.atomLabels[i] ? exitDist(b.atomLabels[i].rects, ux, uy) + 1 : 0, tq = b.atomLabels[j] ? exitDist(b.atomLabels[j].rects, -ux, -uy) + 1 : 0;
        return [{ x: p.x + ux * tp, y: p.y + uy * tp }, { x: q.x - ux * tq, y: q.y - uy * tq }];
      }));
      let onLabel = 0, onBond = 0;
      for (const e of els) for (const o of (e.pair ? [3.6, -3.6] : [0])) {
        const d = { x: e.x + o * Math.cos(e.angle), y: e.y + o * Math.sin(e.angle) };
        if (atoms.some(t => t.lab && t.lab.rects.some(r => rectDist(d.x - t.P.x, d.y - t.P.y, r) < 2.4))) onLabel++;
        if (segs.some(([p, q]) => segDist(d, p, q) < 2.4 + 1)) onBond++;
      }
      if (onLabel) problems.push(onLabel + ' dot(s) on a label');
      if (onBond) problems.push(onBond + ' dot(s) on a bond');
      // the groups of dots on one atom read as groups: two groups' nearest dots further apart than a pair's own two (7.2 px)
      const dotsOf = e => (e.pair ? [3.6, -3.6] : [0]).map(o => ({ x: e.x + o * Math.cos(e.angle), y: e.y + o * Math.sin(e.angle) }));
      let merged = 0;
      for (let i = 0; i < els.length; i++) for (let j = i + 1; j < els.length; j++) {
        if (owner.get(els[i]) !== owner.get(els[j])) continue;
        const d = Math.min(...dotsOf(els[i]).flatMap(p => dotsOf(els[j]).map(q => Math.hypot(p.x - q.x, p.y - q.y))));
        if (d < 10) merged++;
      }
      if (merged) problems.push(merged + ' pair(s) of electron groups on one atom closer than 10 px (they read as one run of dots)');
      // the arrows' ends
      const near = (p, q, r) => Math.hypot(p.x - q.x, p.y - q.y) < r;
      for (const a of L.arrowInfo || []) {
        const it = a.item, tag = `step ${a.step} ${a.raw}`;
        if (a.tk === 'bond' && a.far && !a.mate && a.X && !near(it.p2, a.X, 22)) problems.push(`${tag}: a new bond between far-apart atoms, but the head is ${Math.round(Math.hypot(it.p2.x - a.X.x, it.p2.y - a.X.y))} px from the atom attacked`);
        if (a.tailLines) {
          const [l] = a.tailLines.map(([p, q]) => ({ p, q, d: segDist(it.p1, p, q) })).sort((x, y) => x.d - y.d);
          if (!l || l.d > 3.5) problems.push(`${tag}: the tail is ${l ? l.d.toFixed(1) : '?'} px off its bond's lines`);
          else { const L2 = Math.hypot(l.q.x - l.p.x, l.q.y - l.p.y) || 1, t = ((it.p1.x - l.p.x) * (l.q.x - l.p.x) + (it.p1.y - l.p.y) * (l.q.y - l.p.y)) / L2; if (Math.abs(t - L2 / 2) > 4) problems.push(`${tag}: the tail is ${Math.abs(t - L2 / 2).toFixed(1)} px from the middle of its bond`); }
        }
        const ctl = { x: (it.p1.x + it.p2.x) / 2 + 2 * it.bend.x, y: (it.p1.y + it.p2.y) / 2 + 2 * it.bend.y }, arrive = { x: it.p2.x - ctl.x, y: it.p2.y - ctl.y };
        const angle = (u, v) => Math.acos(Math.max(-1, Math.min(1, (u.x * v.x + u.y * v.y) / ((Math.hypot(u.x, u.y) * Math.hypot(v.x, v.y)) || 1)))) * 180 / Math.PI;
        // a lone pair pushed into a bond of its own atom: the head ends on the middle of that bond, coming across it
        if (a.lpOwn) {
          const dm = Math.hypot(it.p2.x - a.tgtMid.x, it.p2.y - a.tgtMid.y), dl = Math.min(...a.tgtLines.map(([p, q]) => segDist(it.p2, p, q)));
          const [p, q] = a.tgtLines[0], along = angle(arrive, { x: q.x - p.x, y: q.y - p.y }), lq = Math.hypot(q.x - p.x, q.y - p.y) || 1;
          const sl = Math.abs(((it.p2.x - a.tgtMid.x) * (q.x - p.x) + (it.p2.y - a.tgtMid.y) * (q.y - p.y)) / lq);   // along the bond, from its middle
          if (sl > 6 || dl > 5) problems.push(`${tag}: lone pair into its own bond, but the head is ${dm.toFixed(1)} px from the bond's middle (${dl.toFixed(1)} px off its line)`);
          if (Math.min(along, 180 - along) < 40) problems.push(`${tag}: lone pair into its own bond, but the head comes in along the bond (${Math.round(Math.min(along, 180 - along))}°)`);
          if (Math.hypot(it.p2.x - it.p1.x, it.p2.y - it.p1.y) < 16) problems.push(`${tag}: lone pair into its own bond: a stub of an arrow`);
        }
        // an arrow onto an atom (a lone pair, or a bond's electrons onto one of its own atoms) ends at its letters, pointing into it
        if ((a.tk === 'own' || a.tk === 'atom') && a.X && a.Xrects.length) {
          const d = Math.hypot(it.p2.x - a.X.x, it.p2.y - a.X.y), off = angle(arrive, { x: a.X.x - it.p2.x, y: a.X.y - it.p2.y });
          if (d > 17) problems.push(`${tag}: the head is ${d.toFixed(1)} px from the atom it ends at`);
          if (off > 50) problems.push(`${tag}: the head points ${Math.round(off)}° away from the atom it ends at`);
          // the arrowhead clear of the atom's charge
          if (a.Xcharged && a.Xrects[1]) {
            const v = { x: -arrive.x / (Math.hypot(arrive.x, arrive.y) || 1), y: -arrive.y / (Math.hypot(arrive.x, arrive.y) || 1) };
            const tri = [0, 0.45, -0.45].map(t => ({ x: it.p2.x + (t ? 9.5 : 0) * (v.x * Math.cos(t) - v.y * Math.sin(t)), y: it.p2.y + (t ? 9.5 : 0) * (v.x * Math.sin(t) + v.y * Math.cos(t)) }));
            const dq = Math.min(...tri.map(q => rectDist(q.x - a.X.x, q.y - a.X.y, a.Xrects[1])));
            if (dq < 1.5) problems.push(`${tag}: the arrowhead touches the atom's charge (${dq.toFixed(1)} px)`);
          }
        }
        // an arrow crosses no bond line other than those of the atoms it leaves from and goes to
        {
          const S = it.allSamples();
          const cr = (o, p, q) => (p.x - o.x) * (q.y - o.y) - (p.y - o.y) * (q.x - o.x);
          let x = 0;
          for (const [p, q] of segs) {
            const P0 = a.ends.find(e => Math.hypot(e.x - p.x, e.y - p.y) < 14), Q0 = a.ends.find(e => Math.hypot(e.x - q.x, e.y - q.y) < 14);
            if (P0 && Q0) continue;
            for (let i = 1; i < S.length; i++) if (cr(S[i - 1], S[i], p) * cr(S[i - 1], S[i], q) < 0 && cr(p, q, S[i - 1]) * cr(p, q, S[i]) < 0) x++;
          }
          if (x) problems.push(`${tag}: the arrow crosses ${x} bond line(s)`);
        }
        for (const b of L.arrowInfo) {
          if (b === a || b.step !== a.step) continue;
          // two arrows of one step that cross cannot be told apart (two fishhooks meeting at one new bond may touch)
          if (L.arrowInfo.indexOf(b) > L.arrowInfo.indexOf(a) && !(a.mate && b.mate)) {
            const P = a.item.allSamples(), Q = b.item.allSamples(), cr = (o, p, q) => (p.x - o.x) * (q.y - o.y) - (p.y - o.y) * (q.x - o.x);
            let x = 0; for (let i = 1; i < P.length; i++) for (let j = 1; j < Q.length; j++) if (cr(P[i - 1], P[i], Q[j - 1]) * cr(P[i - 1], P[i], Q[j]) < 0 && cr(Q[j - 1], Q[j], P[i - 1]) * cr(Q[j - 1], Q[j], P[i]) < 0) x++;
            if (x) problems.push(`${tag} and ${b.raw}: the arrows cross`);
          }
          if (near(it.p2, b.item.p2, 8) && L.arrowInfo.indexOf(b) > L.arrowInfo.indexOf(a)) problems.push(`${tag} and ${b.raw}: two heads on one spot`);
          if (near(it.p2, b.item.p1, 10) || (b.src && near(it.p2, b.src, 10))) problems.push(`${tag}: the head lands on the tail of ${b.raw}`);
        }
      }
      // map numbers never show: not in a box's SMILES, not in any text drawn
      if (boxes.some(b => /:\d+\]/.test(b.smiles))) problems.push('a map number in a box SMILES');
      if (items.some(i => /\[[^\]]*:\d+\]/.test(i.g.textContent || ''))) problems.push('a map number in the drawing');
      // rows: no two overlap; within a row no two molecule boxes overlap; the reaction arrow is clear of the curved ones
      const meet = (a, b) => a.x < b.x + b.w && b.x < a.x + a.w && a.y < b.y + b.h && b.y < a.y + a.h;
      for (let i = 0; i < L.rows.length; i++) for (let j = i + 1; j < L.rows.length; j++) if (meet(L.rows[i].box, L.rows[j].box)) problems.push(`rows ${L.rows[i].step} and ${L.rows[j].step} overlap`);
      for (const r of L.rows) {
        const bs = r.items.filter(i => i instanceof D.MoleculeBox).map(i => i.bbox());
        for (let i = 0; i < bs.length; i++) for (let j = i + 1; j < bs.length; j++) if (meet(bs[i], bs[j])) problems.push(`step ${r.step}: two molecule boxes overlap`);
        const fwd = r.items.find(i => i instanceof D.ArrowItem && i.kind === 'forward'), curls = r.items.filter(i => i instanceof D.ArrowItem && i.kind !== 'forward');
        if (fwd && curls.some(c => c.allSamples().some(q => q.x > fwd.p1.x - 2))) problems.push(`step ${r.step}: a curved arrow reaches the reaction arrow`);
        if (r.items.some(i => !L.items.includes(i))) problems.push('a row item is not in the drawing');
      }
      const head = items.filter(i => !L.rows.some(r => r.items.includes(i)));
      if (L.rows.length && head.some(h => L.rows.some(r => meet(h.bbox(), r.box)))) problems.push('the heading overlaps a step');
      // arrows through labels other than their own ends: counted, not failed (a crowded step may need it)
      let through = 0;
      for (const a of arrows) {
        const pts = a.allSamples().slice(5, -5);
        if (pts.some(q => atoms.some(t => t.lab && Math.hypot(t.P.x - a.p1.x, t.P.y - a.p1.y) > 24 && Math.hypot(t.P.x - a.p2.x, t.P.y - a.p2.y) > 24 && t.lab.rects.some(r => rectDist(q.x - t.P.x, q.y - t.P.y, r) < 1)))) through++;
      }
      return { problems, arrows: arrows.length, boxes: boxes.length, electrons: els.length, lpAtoms, through, top: L.top, bottom: L.bottom, mode: D.mode, selected: D.selection.length };
    };
  });

  /* every workbook mechanism, as the Reactions page sends it (plain JSON, onExample) */
  const ids = await page.evaluate(() => window.REACTION_RULES.rules.filter(r => (r.mechanism || []).length).map(r => r.id));
  check('52 rules with a mechanism', ids.length === 52, ids.length);
  let totalThrough = 0, totalArrows = 0, lpAtoms = 0;
  const shotIds = { R001: 'Fischer esterification', R004: 'radical HBr', R013: 'alkyne hydrogenation', R026: 'SN2', R029: 'E2', R009: 'epoxidation', R041: 'Wittig' };
  for (const id of ids) {
    const r = await page.evaluate(id => {
      const rule = window.REACTION_RULES.rules.find(x => x.id === id), txt = x => (x == null ? '' : String(x));
      const code = txt(rule.code).split(' | ')[0].split('>');
      const steps = rule.mechanism.map((st, i) => {
        const out = { step: i + 1, type: txt(st.type), from: txt(st.from), to: txt(st.to), made: txt(st.made), broken: txt(st.broken), note: txt(st.note), smiles: txt(st.smiles) };
        if (st.before) { out.before = String(st.before); out.arrows = txt(st.arrows); }
        return out;
      });
      const payload = structuredClone({ title: 'Mechanism for ' + rule.name, name: rule.name, ruleId: rule.id, status: rule.mechStatus || '', onExample: true,
        reaction: { reactants: (code[0] || '').split('.').filter(Boolean), reagent: code[1] || '', products: (code[2] || '').split('.').filter(Boolean), conditions: txt(rule.conditions) }, steps });
      const toks = steps.filter(s => s.before && s.arrows).flatMap(s => s.arrows.split(';').map(t => t.trim()).filter(Boolean));
      window.Draw.clearAll();
      let ok, err = null;
      try { ok = window.Draw.showMechanism(payload); } catch (e) { err = String(e && e.stack || e); }
      if (err) return { err };
      return Object.assign({ ok, scroll: window.Draw.scroll.scrollTop }, window.__mechCheck({ tokens: toks.length, fish: toks.filter(t => t.startsWith('~')).length }));
    }, id);
    if (r.err) { check(`${id}: no exception`, false, r.err); continue; }
    check(`${id}: drawn (${r.arrows} arrows, ${r.boxes} molecules, ${r.electrons} electron items)`, r.ok === true && !r.problems.length && r.mode === 'select' && r.selected === 0, r.problems);
    totalThrough += r.through; totalArrows += r.arrows; lpAtoms += r.lpAtoms;
    if (shotIds[id]) {
      await page.evaluate(() => { window.Draw.scroll.scrollTop = 0; window.Draw.scroll.scrollLeft = 0; });
      const sb = await page.locator('.dr-scroll svg').boundingBox();
      const w = await page.evaluate(() => Math.max(...window.Draw.lastMechanism.items.map(i => { const b = i.bbox(); return b.x + b.w; })));
      await page.setViewportSize({ width: 1500, height: Math.ceil(Math.max(1000, r.bottom + 120)) });
      await page.screenshot({ path: `${shots}/mechanism_${id}.png`, clip: { x: sb.x, y: sb.y, width: Math.min(1490 - sb.x, w + 20), height: r.bottom + 20 } });
      await page.setViewportSize({ width: 1500, height: 1000 });
      console.log(`      screenshot ${shots}/mechanism_${id}.png (${shotIds[id]})`);
    }
  }
  check('O, N, S and halogen atoms carry their lone pairs (counted exactly per atom above)', lpAtoms > 300, lpAtoms);
  console.log(`      ${totalThrough} of ${totalArrows} curved arrows pass over an atom label on the way (not a failure)`);

  /* an engine step built by hand: SN2 of bromoethane by HO⁻ (graph laid out, ids kept, gArrows lp→atom and bond→atom) */
  const eng = await page.evaluate(() => {
    const graph = { atoms: [
      { id: 1, element: 'C', charge: 0, radical: 0, x: 0, y: 0 }, { id: 2, element: 'C', charge: 0, radical: 0, x: 1.3, y: -0.75 },
      { id: 3, element: 'Br', charge: 0, radical: 0, x: 2.6, y: 0 }, { id: 7, element: 'O', charge: -1, radical: 0, x: -0.2, y: -2.2 } ],
      bonds: [{ a: 1, b: 2, order: 1 }, { a: 2, b: 3, order: 1 }] };
    const payload = { title: 'Mechanism for SN2 substitution', name: 'SN2 substitution', ruleId: 'SN2', status: '', onExample: false,
      reaction: { reactants: ['CCBr'], reagent: '[OH-]', products: ['CCO'], conditions: 'DMSO' },
      steps: [{ step: 1, type: 'SN2: the nucleophile attacks as the leaving group leaves', note: 'Backside attack.', smiles: 'CCO.[Br-]', graph,
        gArrows: [{ from: { kind: 'lp', atom: 7 }, to: { kind: 'atom', atom: 2 }, half: false, label: 'lone pair forms the new C–O bond' },
                  { from: { kind: 'bond', a: 2, b: 3 }, to: { kind: 'atom', atom: 3 }, half: false, label: 'C–Br electrons leave with Br' }] }] };
    window.Draw.clearAll();
    const ok = window.Draw.showMechanism(structuredClone(payload));
    const box = window.Draw.lastMechanism.items.find(i => i instanceof window.Draw.MoleculeBox);
    const lens = box.lay.bonds.map(b => { const p = box.lay.atoms.find(a => a.id === b.a), q = box.lay.atoms.find(a => a.id === b.b); return Math.hypot(p.x - q.x, p.y - q.y); });
    return Object.assign({ ok, ids: box.lay.atoms.map(a => a.id), lens, oxy: box.atomPoint(7) && box.atomPoint(2) ? Math.hypot(box.atomPoint(7).x - box.atomPoint(2).x, box.atomPoint(7).y - box.atomPoint(2).y) : 0 }, window.__mechCheck({ tokens: 2, fish: 0 }));
  });
  check('engine step: drawn with both arrows, HO⁻ and Br with their lone pairs', eng.ok === true && !eng.problems.length && eng.arrows === 2, eng.problems);
  check('engine step: setGraph keeps the atom ids and scales the bonds to Draw\'s 46 px', eng.ids.join() === '1,2,3,7' && eng.lens.every(l => Math.abs(l - 46) < 0.5), eng);
  check('engine step: the reagent stays where the engine parked it (scaled)', Math.abs(eng.oxy - Math.hypot(1.5, 1.45) / 1.5 * 46) < 1, eng.oxy);
  {
    const sb = await page.locator('.dr-scroll svg').boundingBox();
    await page.screenshot({ path: `${shots}/mechanism_engine_SN2.png`, clip: { x: sb.x, y: sb.y, width: 800, height: eng.bottom + 20 } });
    console.log(`      screenshot ${shots}/mechanism_engine_SN2.png (engine form)`);
  }

  /* a second mechanism goes under the first, and the canvas scrolls to it */
  const two = await page.evaluate(() => {
    const D = window.Draw, rule = window.REACTION_RULES.rules.find(r => r.id === 'R042');
    const pay = { title: 'Mechanism for ' + rule.name, status: rule.mechStatus, onExample: true, reaction: {}, steps: rule.mechanism.map((s, i) => ({ step: i + 1, type: s.type, note: s.note, smiles: s.smiles, before: s.before, arrows: s.arrows || '' })) };
    const firstBottom = D.lastMechanism.bottom, n0 = D.items.length;
    D.showMechanism(pay);
    const L = D.lastMechanism;
    D.showMechanism(pay);                       // and a third: the canvas grows to hold it
    const L3 = D.lastMechanism;
    return { firstBottom, top: L.top, n0, n: D.items.length, scroll: D.scroll.scrollTop, top3: L3.top, bottom3: L3.bottom, sceneH: D.sceneSize.h, svgH: +D.svg.getAttribute('height') };
  });
  check('a second mechanism starts under the first', two.top > two.firstBottom && two.n > two.n0, two);
  check('the canvas grows to hold a mechanism drawn below the others, and scrolls to it', two.sceneH >= two.bottom3 && two.svgH >= two.bottom3 && Math.abs(two.scroll - (two.top3 - 20)) < 2, two);

  /* the new arrow kinds behave like the other arrows: drawn with the toolbar, bent, selected, coloured, copied, exported */
  await page.evaluate(() => window.Draw.clearAll());
  await page.evaluate(() => { window.Draw.scroll.scrollTop = 0; window.Draw.scroll.scrollLeft = 0; });
  const svg = page.locator('.dr-scroll svg'), o = await svg.boundingBox();
  const at = (x, y) => [o.x + x, o.y + y];
  const drag = async (x1, y1, x2, y2) => { await page.mouse.move(...at(x1, y1)); await page.mouse.down(); await page.mouse.move(...at((x1 + x2) / 2, (y1 + y2) / 2), { steps: 3 }); await page.mouse.move(...at(x2, y2), { steps: 3 }); await page.mouse.up(); };
  check('no curved-arrow tool in the toolbar (the mechanisms still draw curved arrows)', await page.locator('.dr-bar button', { hasText: '↷' }).count() === 0);
  // each toolbar arrow is drawn straight and shows its three circles; the tool stays, and its circles still work
  for (const [k, sym, y] of [['forward', '→', 120], ['equilibrium', '⇌', 220], ['resonance', '↔', 320], ['fishhook', '⇀', 420]]) {
    await page.locator('.dr-bar button', { hasText: new RegExp('^' + sym + '$') }).click(); await drag(100, y, 260, y);
    const st = await page.evaluate(() => { const a = window.Draw.items[window.Draw.items.length - 1]; return { kind: a.kind, bent: a.bent(), circles: a.g.querySelectorAll('circle').length, sel: window.Draw.selection.length }; });
    check(`${k}: drawn straight, selected, three circles`, st.kind === k && !st.bent && st.circles === 3 && st.sel === 1, st);
    await drag(260, y, 320, y);
    const len = await page.evaluate(() => { const a = window.Draw.items[window.Draw.items.length - 1]; return Math.round(a.p2.x - a.p1.x); });
    check(`${k}: dragging the end circle (tool still on) stretches it`, Math.abs(len - 220) < 4 && await page.evaluate(n => window.Draw.items.length === n, [120, 220, 320, 420].indexOf(y) + 1), len);
    const hm = await page.evaluate(() => window.Draw.items[window.Draw.items.length - 1].handle());
    await drag(hm.x, hm.y, hm.x, hm.y - 40);
    check(`${k}: dragging the middle circle bends it`, await page.evaluate(() => window.Draw.items[window.Draw.items.length - 1].bent()));
  }
  const fish = await page.evaluate(() => { const a = window.Draw.items[3]; return { head: a.g.querySelectorAll('polygon').length, barb: a.g.querySelectorAll('line').length }; });
  check('the half-head arrow has one barb, no full head', fish.head === 0 && fish.barb === 1, fish);
  await page.keyboard.press('Escape');
  await page.evaluate(() => window.Draw.clearAll());
  await page.evaluate(() => { const D = window.Draw, a = new D.ArrowItem('curved', { x: 100, y: 200 }, { x: 260, y: 200 }); D.addItem(a); a.autoBend();
    const f = new D.ArrowItem('fishhook', { x: 100, y: 340 }, { x: 260, y: 340 }); D.addItem(f); f.bend = { x: 0, y: -48 }; f.render(); });
  await page.keyboard.press('Escape');
  const h = await page.evaluate(() => window.Draw.items[0].handle());
  await page.mouse.click(...at(h.x, h.y));
  check('clicking the curve selects the curved arrow', await page.evaluate(() => window.Draw.selection.length === 1 && window.Draw.selection[0].kind === 'curved'));
  await drag(h.x, h.y, h.x, h.y - 40);
  const bend = await page.evaluate(() => window.Draw.items[0].bend);
  check('dragging the middle bends it further', bend.y < -60, bend);
  await drag(260, 200, 300, 260);
  check('dragging an end moves it', await page.evaluate(() => { const a = window.Draw.items[0]; return Math.abs(a.p2.x - 300) < 3 && Math.abs(a.p2.y - 260) < 3; }));
  await page.locator('.dr-bar button', { hasText: 'Color' }).click(); await page.locator('.dr-color i[data-c="#e03131"]').first().click();
  check('Color recolours it', await page.evaluate(() => window.Draw.items[0].color === '#e03131' && window.Draw.items[0].g.querySelector('path').getAttribute('stroke') === null && window.Draw.items[0].g.querySelector('g').getAttribute('stroke') === '#e03131'));
  await page.mouse.click(...at(900, 700));
  await page.keyboard.press('Escape');
  const hb = await page.evaluate(() => window.Draw.items[1].handle());
  await page.mouse.click(...at(hb.x, hb.y));
  await page.keyboard.press('Control+c'); await page.keyboard.press('Control+v');
  const pasted = await page.evaluate(() => { const [a, , c] = window.Draw.items; const b = window.Draw.items[1]; return { n: window.Draw.items.length, kind: c && c.kind, bend: c && c.bend, orig: b.bend, sel: window.Draw.selection.length }; });
  check('copy and paste keeps a fishhook a fishhook, with its bend', pasted.n === 3 && pasted.kind === 'fishhook' && pasted.bend.x === pasted.orig.x && pasted.bend.y === pasted.orig.y && pasted.sel === 1, pasted);
  const [dl] = await Promise.all([page.waitForEvent('download'), page.locator('.dr-bar button', { hasText: 'Export PNG' }).click()]);
  const buf = require('fs').readFileSync(await dl.path());
  check('Export PNG with mechanism arrows saves a PNG', buf.slice(1, 4).toString() === 'PNG' && buf.length > 2000, buf.length);

  check('no page errors', errors.length === 0, errors);
  await browser.close();
  console.log(`${n - fails} of ${n} checks passed` + (fails ? ` — ${fails} FAILED` : ''));
  process.exit(fails ? 1 : 0);
})();
