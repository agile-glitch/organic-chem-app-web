/* Draw.showMechanism(payload): a reaction's mechanism laid out on the Draw canvas the way a textbook draws it. One row
   per step: the species at the start of the step with every lone pair (••) and unpaired electron (•) on their atoms, red
   curved arrows for the electrons that move (full head: a pair; fishhook: one electron), a reaction arrow, and what the
   step makes; the step's note under the row. Everything is an ordinary Draw item (molecule boxes, text, arrows,
   electrons), so it can be moved, recoloured, edited, copied and exported like anything drawn by hand.

   The payload is plain JSON (it crosses to the home screen by postMessage), from the Reactions page's "Open Mechanism in
   Draw" button (mechanismPayload in js/reactions.js):
     {title, name, ruleId, status, onExample, reaction: {reactants, reagent, products, conditions}, steps}
   Each step keeps its words (type, note, ...) and what it makes (smiles), and draws its start either from the workbook's
   atom-mapped 'before' SMILES with the Mechanism sheet's arrow notation ("5>1" lone pair on atom map 5 to atom 1,
   "1-2>3" bond 1–2 onto atom 3, "1-2>2-4" into a new bond 2–4, "~" in front = one electron), or from an engine graph
   {atoms, bonds} (Chem.buildMechanism, already laid out) with gArrows ({from: lp | bond, to: atom | bond, half}).

   The species are turned (30° steps) and, when they carry no stereo marks, mirrored so the atoms that bond face each
   other, an attack comes from the side away from the group that leaves, and the arrows do not cross; a chain whose two
   reacting ends point apart (a diene drawn s-trans) is turned about a single bond to bring them together.
   The arrows follow textbook rules: an arrow leaves the lone pair (or the unpaired electron) that faces its target, or the
   middle of a bond; it ends at the target atom's letter (never its charge or its H's), at the middle of a bond that gains
   a π bond, at the middle of a new bond whose two atoms are already close, and at the attacked atom when they are far
   apart (the middle of empty space between two molecules points at nothing). Two fishhooks making one bond end side by
   side, and two leaving one bond start on either side of it. Each arrow bends to whichever side keeps it clear of atoms,
   bonds, dots, the other arrows' ends and the step's other arrows (crossing one costs most).
   Draw.lastMechanism says what was drawn (for tools/browser_mechanism.js). */
(() => {
  'use strict';
  const D = window.Draw, C = window.Chem;
  if (!D || !C) return;

  const RED = '#d62828', TAU = Math.PI * 2, DEG = Math.PI / 180;
  const X0 = 40, Y0 = 40;                      // where a mechanism starts on an empty canvas
  const GAP = 18;                              // between two boxes of a row (a box already keeps 30 px round its molecule)
  const FWD_GAP = 6, FWD_LEN = 70;             // the black reaction arrow between the start and what the step makes
  const BLOCK_GAP = 28, NOTE_PX = 13, PLUS_PX = 20;
  const DOT_R = 2.4, PAIR_HALF = 3.6;          // ElectronItem: dot radius, and each dot of a pair 3.6 px off its centre
  const DOT_CLEAR = DOT_R + 2.2;               // a dot keeps this far from label ink and from bonds
  const GROUP_SEP = 13;                        // two electron groups on one atom: nearest dots this far apart (a pair's own: 7.2)
  const BOND = 46;                             // Draw's bond length
  const FAR = 1.35 * BOND;                     // a new bond longer than this: its arrow ends at the atom attacked
  const HEAD_GAP = 3.5, BARE_GAP = 9;          // an arrow's tip stops this far from the target's letters (a bare carbon: from its vertex)
  const FISH_GAP = 7;                          // two fishhooks making one bond: each tip this far to its own side of where it forms
  const SHARE_OFF = 2.5;                       // two arrows leaving one bond: their tails this far either side of its middle

  const norm = a => ((a % TAU) + TAU) % TAU;
  const angDist = (p, q) => { const d = Math.abs(norm(p) - norm(q)); return d > Math.PI ? TAU - d : d; };
  const mid = (p, q) => ({ x: (p.x + q.x) / 2, y: (p.y + q.y) / 2 });
  const dist = (p, q) => Math.hypot(p.x - q.x, p.y - q.y);
  const unit = (x, y) => { const l = Math.hypot(x, y) || 1; return { x: x / l, y: y / l }; };
  const mean = xs => xs.reduce((s, x) => s + x, 0) / (xs.length || 1);
  const along = (p, u, t) => ({ x: p.x + u.x * t, y: p.y + u.y * t });
  const segDist = (p, a, b) => {
    const vx = b.x - a.x, vy = b.y - a.y, l2 = vx * vx + vy * vy || 1;
    const t = Math.max(0, Math.min(1, ((p.x - a.x) * vx + (p.y - a.y) * vy) / l2));
    return Math.hypot(p.x - (a.x + vx * t), p.y - (a.y + vy * t));
  };
  const cross2 = (o, a, b) => (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x);
  const segsCross = (p1, p2, q1, q2) => cross2(p1, p2, q1) * cross2(p1, p2, q2) < 0 && cross2(q1, q2, p1) * cross2(q1, q2, p2) < 0;
  const rectDist = (x, y, r) => Math.hypot(Math.max(r.x0 - x, 0, x - r.x1), Math.max(r.y0 - y, 0, y - r.y1));
  /* how far a ray from the atom's centre runs before it leaves the label's ink (rects about the atom) */
  function exitDist(rects, ux, uy) {
    let best = 0;
    for (const r of rects) {
      let t0 = -Infinity, t1 = Infinity, hit = true;
      for (const [o0, o1, u] of [[r.x0, r.x1, ux], [r.y0, r.y1, uy]]) {
        if (Math.abs(u) < 1e-9) { if (o0 > 0 || o1 < 0) { hit = false; break; } continue; }
        let a = o0 / u, b = o1 / u; if (a > b) [a, b] = [b, a];
        t0 = Math.max(t0, a); t1 = Math.min(t1, b);
      }
      if (hit && t1 >= Math.max(t0, 0)) best = Math.max(best, t1);
    }
    return best;
  }
  const quad = (s, c, e, t) => { const u = 1 - t; return { x: u * u * s.x + 2 * u * t * c.x + t * t * e.x, y: u * u * s.y + 2 * u * t * c.y + t * t * e.y }; };
  const control = (s, e, bend) => { const m = mid(s, e); return { x: m.x + 2 * bend.x, y: m.y + 2 * bend.y }; };   // ArrowItem's: the curve passes through mid + bend
  const curve = (s, e, bend, n = 24) => { const c = control(s, e, bend), out = []; for (let i = 0; i <= n; i++) out.push(quad(s, c, e, i / n)); return out; };
  const crossings = (a, b) => { let n = 0; for (let i = 1; i < a.length; i++) for (let j = 1; j < b.length; j++) if (segsCross(a[i - 1], a[i], b[j - 1], b[j])) n++; return n; };

  function wrap(str, n = 100) {                // long text broken into lines of about n characters, at spaces
    const out = [];
    for (const para of String(str || '').split('\n')) {
      let line = '';
      for (const w of para.split(/\s+/).filter(Boolean)) { if (line && (line + ' ' + w).length > n) { out.push(line); line = w; } else line = line ? line + ' ' + w : w; }
      out.push(line);
    }
    return out.join('\n').trim();
  }
  function textItem(str, size, opts = {}) { const t = new D.TextBox(str); t.size = size; t.bold = !!opts.bold; t.color = opts.color || '#000'; t.render(); return t; }
  function plusItem(cx, cy) {                  // a '+' whose glyph is centred on (cx, cy)
    const t = textItem('+', PLUS_PX), b = t.bbox();
    t.x = cx - b.w / 2; t.y = cy - (4 + PLUS_PX * 0.95 - PLUS_PX * 0.34); t.render();
    return t;
  }
  /* a species for the heading: a formula for a small inorganic one (HBr, H₂O, CN⁻; a salt of them: NaCN), else its common
     name, else its IUPAC name, else the SMILES as given. The Reactions page hands over SMILES, and 'Br' (HBr) reads as a
     bromine atom */
  const SUB = '₀₁₂₃₄₅₆₇₈₉', SUP = { '+': '⁺', '-': '⁻', 2: '²', 3: '³' };
  const subDigits = f => f.replace(/\d/g, d => SUB[+d]);
  function readable(smi) {
    const s = String(smi || '').trim(); if (!s) return '';
    const parse = x => { try { const g = C.parseSmiles(x); return g && g.atoms.length ? g : null; } catch (e) { return null; } };
    const g = parse(s); if (!g) return s;
    const heavy = gr => gr.atoms.filter(a => a.element !== 'H');
    const tiny = gr => { const hv = heavy(gr); return hv.length <= 3 && (!hv.some(a => a.element === 'C') || !/H/.test(C.formula(gr))); };
    const chg = gr => gr.atoms.reduce((t, a) => t + (a.charge || 0), 0);
    const formula = (gr, withCharge) => { const q = chg(gr); return subDigits(C.formula(gr)) + (withCharge && q ? (Math.abs(q) > 1 ? SUP[Math.abs(q)] || '' : '') + SUP[q > 0 ? '+' : '-'] : ''); };
    const named = gr => {
      try { const n = C.commonNameFor(gr); if (n) return n; } catch (e) { /* none */ }
      try { const n = C.iupacName(gr); if (n && !(C.nameFailed && C.nameFailed(n)) && !/\+/.test(n)) return n; } catch (e) { /* none */ }
      return null;
    };
    const parts = s.split('.').filter(Boolean);
    if (parts.length === 1) return tiny(g) ? formula(g, true) : named(g) || s;
    const whole = (() => { try { return C.commonNameFor(g); } catch (e) { return null; } })();
    if (whole) return whole;
    const gs = parts.map(parse);
    if (gs.every(x => x && tiny(x))) return gs.slice().sort((p, q) => chg(q) - chg(p)).map(x => formula(x, false)).join('');   // a salt: NaCN
    return parts.map((p, k) => (gs[k] ? (tiny(gs[k]) ? formula(gs[k], true) : named(gs[k]) || p) : p)).join(' + ');
  }

  /* ---------- the species ---------- */
  const stripMaps = s => String(s || '').replace(/\[([^\]\[:]+):\d+\]/g, '[$1]');   // [CH2:2] → [CH2]: the drawing never shows map numbers
  /* a laid-out molecule for a SMILES, copied from the last one made for the same SMILES (a step's product is usually the
     next step's start). Without stereo marks the parse's own layout is used as it is: the box's setMolecule lays it out
     again keeping R/S and E/Z, and working those out takes over half a second for a triphenylphosphonium betaine, for
     nothing when the SMILES gives no stereo. null when the SMILES cannot be read */
  const made = new Map();
  function moleculeBox(smiles) {
    let t = made.get(smiles);
    if (t === undefined) {
      t = new D.MoleculeBox();
      if (/[@\/\\]/.test(smiles)) t.setMolecule(smiles);
      else { let g = null; try { g = C.parseSmiles(smiles); } catch (e) { g = null; } t.setGraph(g, smiles); }
      if (t.lay) { straighten(t); if (placeHs(t.lay)) t.refreshDisplay(); } else t = null;
      if (made.size > 300) made.clear();
      made.set(smiles, t);
    }
    if (!t) return null;
    const b = D.cloneItem(t); b.angle = 0; b.refreshDisplay();
    return b;
  }
  /* map number → {i: atom index, el: element}, read from the SMILES text: the parser numbers atoms in the order written */
  function mapIndex(piece) {
    const out = new Map(), re = /\[([^\]]+)\]|Br|Cl|[BCNOPSFIbcnops]/g;
    let m, k = 0;
    while ((m = re.exec(piece))) {
      const mm = m[1] && m[1].match(/^\d*([A-Z][a-z]?|[a-z]+)[^:]*:(\d+)$/);
      if (mm) out.set(+mm[2], { i: k, el: mm[1][0].toUpperCase() + mm[1].slice(1) });
      k++;
    }
    return out;
  }
  function smilesSpecies(piece) {              // {box, byMap: map number → atom id} or {text} when it cannot be read
    const clean = stripMaps(piece), box = moleculeBox(clean);
    if (!box) return { text: clean };
    const byMap = new Map(), idx = mapIndex(piece), at = box.lay.atoms;
    if ([...idx.values()].every(v => at[v.i] && at[v.i].element === v.el)) idx.forEach((v, n) => byMap.set(n, at[v.i].id));
    else {                                     // (not expected) the text and the parse disagree: draw the mapped parse itself
      let mapped = null; try { mapped = C.parseSmiles(piece); } catch (e) { mapped = null; }
      if (mapped) { box.setGraph(mapped); straighten(box); placeHs(box.lay); box.refreshDisplay(); mapped.atoms.forEach((a, i) => { if (a.mapNum && box.lay.atoms[i]) byMap.set(a.mapNum, box.lay.atoms[i].id); }); }
    }
    return { box, byMap };
  }
  /* Chem's layout draws C≡C, C≡N and C=C=C bent; a textbook draws them straight. Each such atom with two neighbours gets
     its smaller side swung round it until the two bonds are in line (a turn: wedges keep their meaning) */
  function straighten(box) {
    const lay = box.lay; if (!lay || lay.atoms.length < 3) return;
    const byId = new Map(lay.atoms.map(a => [a.id, a]));
    const nbrs = id => lay.bonds.filter(b => b.a === id || b.b === id).map(b => ({ id: b.a === id ? b.b : b.a, order: b.order }));
    const side = (from, block) => {                // the atoms reached from `from` without passing `block`
      const seen = new Set([from]), todo = [from];
      while (todo.length) { const id = todo.pop(); for (const n of nbrs(id)) if (n.id !== block && !seen.has(n.id)) { seen.add(n.id); todo.push(n.id); } }
      return seen;
    };
    let moved = false;
    for (const a of lay.atoms) {
      const ns = nbrs(a.id);
      if (ns.length !== 2 || !(ns.some(n => n.order === 3) || ns.every(n => n.order === 2))) continue;
      const sides = ns.map((n, k) => ({ n, other: ns[1 - k], set: side(n.id, a.id) })).filter(s => !s.set.has(s.other.id)).sort((p, q) => p.set.size - q.set.size);
      if (!sides.length) continue;                 // in a ring: leave it
      const s = sides[0], o = byId.get(s.other.id), m = byId.get(s.n.id);
      const delta = Math.atan2(a.y - o.y, a.x - o.x) - Math.atan2(m.y - a.y, m.x - a.x);
      if (Math.abs(Math.atan2(Math.sin(delta), Math.cos(delta))) < 1e-3) continue;
      const cs = Math.cos(delta), sn = Math.sin(delta);
      for (const id of s.set) { const p = byId.get(id), dx = p.x - a.x, dy = p.y - a.y; p.x = a.x + dx * cs - dy * sn; p.y = a.y + dx * sn + dy * cs; }
      moved = true;
    }
    if (!moved) return;
    // a molecule that is now one straight line (2-butyne, a nitrile) is laid flat
    const at = lay.atoms; let far = [at[0], at[1], -1];
    for (const p of at) for (const q of at) { const d = Math.hypot(p.x - q.x, p.y - q.y); if (d > far[2]) far = [p, q, d]; }
    const [p0, q0] = far, ux = (q0.x - p0.x) / far[2], uy = (q0.y - p0.y) / far[2];
    if (at.every(r => Math.abs((r.x - p0.x) * uy - (r.y - p0.y) * ux) < 2)) {
      const t = -Math.atan2(uy, ux), cs = Math.cos(t), sn = Math.sin(t);
      at.forEach(r => { const dx = r.x - p0.x, dy = r.y - p0.y; r.x = p0.x + dx * cs - dy * sn; r.y = p0.y + dx * sn + dy * cs; });
    }
    box.refreshDisplay();
  }
  /* an explicit H (an engine graph's, or [H] in a SMILES) as laid out can lie on another bond, crowd another atom, or
     continue a bond in a straight line (H–C2–C3 reading as one bond). Such an H goes into the widest free gap round its
     atom at bond length (with one other bond: 120° from it, a zigzag). An H on a wedge keeps its place. true: one moved */
  function placeHs(lay) {
    if (!lay || lay.atoms.length < 3) return false;
    const at = lay.atoms, byId = new Map(at.map(a => [a.id, a])), nb = new Map(at.map(a => [a.id, []]));
    lay.bonds.forEach(b => { if (nb.has(b.a) && nb.has(b.b)) { nb.get(b.a).push({ id: b.b, bond: b }); nb.get(b.b).push({ id: b.a, bond: b }); } });
    const lens = lay.bonds.map(b => { const p = byId.get(b.a), q = byId.get(b.b); return p && q ? Math.hypot(p.x - q.x, p.y - q.y) : 0; }).filter(x => x > 1e-6);
    const L = lens.length ? mean(lens) : BOND;
    const segs = () => lay.bonds.map(b => [byId.get(b.a), byId.get(b.b), b]).filter(([p, q]) => p && q);
    let moved = false;
    for (const h of at) {
      if (h.element !== 'H' || nb.get(h.id).length !== 1) continue;
      const { id: cid, bond } = nb.get(h.id)[0], c = byId.get(cid);
      if (!c || bond.stereo || c.element === 'H') continue;
      const others = nb.get(c.id).filter(n => n.id !== h.id).map(n => byId.get(n.id));
      if (!others.length) continue;
      const dirs = others.map(o => Math.atan2(o.y - c.y, o.x - c.x)), th = Math.atan2(h.y - c.y, h.x - c.x);
      const badAt = (x, y) => at.some(o => o !== h && o !== c && Math.hypot(o.x - x, o.y - y) < 0.6 * L) ||
        segs().some(([p, q, b]) => b !== bond && p !== c && q !== c && p !== h && q !== h && segDist({ x: (x + c.x) / 2, y: (y + c.y) / 2 }, p, q) < 0.25 * L);
      const onBond = Math.min(...dirs.map(d => angDist(d, th))) < 45 * DEG, inLine = others.length >= 2 && dirs.some(d => angDist(d, th) > 155 * DEG);
      if (!onBond && !inLine && !badAt(h.x, h.y)) continue;
      const cands = others.length === 1 ? [dirs[0] + 120 * DEG, dirs[0] - 120 * DEG] : [spread(dirs, 1)[0]];
      if (others.length >= 2) for (let k = 0; k < 24; k++) cands.push(k * 15 * DEG);
      const score = t => { const x = c.x + L * Math.cos(t), y = c.y + L * Math.sin(t); return Math.min(...dirs.map(d => angDist(d, t))) * 10 + (badAt(x, y) ? -100 : 0); };
      const best = cands.reduce((b, t) => (score(t) > score(b) + 1e-6 ? t : b), cands[0]);
      h.x = c.x + L * Math.cos(best); h.y = c.y + L * Math.sin(best); moved = true;
    }
    return moved;
  }
  /* a 1,2-shift's H goes in the gap round its carbon next to the C+ it moves to, so its arrow need not cross a bond */
  function shiftH(ar) {
    const h = ar.Y, c = same(ar.A, h) ? ar.B : ar.A, lay = h.box.lay;
    const at = id => lay.atoms.find(a => a.id === id), H = at(h.id), Cc = at(c.id), X = at(ar.X.id);
    if (!H || H.element !== 'H' || !Cc || !X || lay.bonds.filter(b => b.a === H.id || b.b === H.id).length !== 1) return false;
    const others = lay.bonds.filter(b => (b.a === Cc.id || b.b === Cc.id) && b.a !== H.id && b.b !== H.id).map(b => at(b.a === Cc.id ? b.b : b.a)).filter(Boolean);
    if (others.length < 2) return false;
    const L = Math.hypot(H.x - Cc.x, H.y - Cc.y) || BOND, dirs = others.map(o => Math.atan2(o.y - Cc.y, o.x - Cc.x)).map(norm).sort((p, q) => p - q);
    const want = Math.atan2(X.y - Cc.y, X.x - Cc.x);
    let best = null;
    dirs.forEach((d, k) => { const e = k + 1 < dirs.length ? dirs[k + 1] : dirs[0] + TAU, mid2 = (d + e) / 2; if (e - d < 50 * DEG) return; if (!best || angDist(mid2, want) < angDist(best, want)) best = mid2; });
    if (best == null) return false;
    H.x = Cc.x + L * Math.cos(best); H.y = Cc.y + L * Math.sin(best);
    return true;
  }
  function graphSpecies(graph) {
    const box = new D.MoleculeBox(); box.setGraph(graph);
    if (box.lay && placeHs(box.lay)) box.refreshDisplay();
    return box.lay ? { box, byId: true } : { text: '(structure)' };
  }
  const hasStereo = box => /[@\/\\]/.test(box.smiles || '') || box.lay.bonds.some(b => b.stereo);
  /* two atoms of one molecule that the step brings together (the ends of a diene, a peracid's H and its C=O) but the layout
     put far apart across a single bond (s-trans): one side is reflected across that bond, a turn about it that changes
     nothing about the molecule, when that brings them at least 10 px closer without crowding another atom */
  function bringTogether(box, pairs) {
    const lay = box.lay; if (!lay || !pairs.length || hasStereo(box)) return false;
    const byId = new Map(lay.atoms.map(a => [a.id, a])), adj = new Map(lay.atoms.map(a => [a.id, []]));
    lay.bonds.forEach(b => { if (adj.has(b.a) && adj.has(b.b)) { adj.get(b.a).push(b.b); adj.get(b.b).push(b.a); } });
    const reach = (from, block) => { const seen = new Set([from]), todo = [from]; while (todo.length) { const id = todo.pop(); for (const n of adj.get(id)) if (n !== block && !seen.has(n)) { seen.add(n); todo.push(n); } } return seen; };
    const path = (p, q) => { const prev = new Map([[p, null]]), todo = [p]; while (todo.length) { const id = todo.shift(); if (id === q) break; for (const n of adj.get(id)) if (!prev.has(n)) { prev.set(n, id); todo.push(n); } } if (!prev.has(q)) return null; const out = []; for (let v = q; v != null; v = prev.get(v)) out.unshift(v); return out; };
    const lens = lay.bonds.map(b => { const p = byId.get(b.a), q = byId.get(b.b); return p && q ? dist(p, q) : 0; }).filter(x => x > 1e-6), L = lens.length ? mean(lens) : BOND;
    let moved = false;
    for (const [p, q] of pairs) for (let round = 0; round < 3; round++) {
      const ids = path(p, q); if (!ids) break;
      const d0 = dist(byId.get(p), byId.get(q)); let best = null;
      for (let k = 0; k + 1 < ids.length; k++) {
        const u = ids[k], v = ids[k + 1], b = lay.bonds.find(x => (x.a === u && x.b === v) || (x.a === v && x.b === u));
        if (!b || b.order !== 1 || adj.get(u).length < 2 || adj.get(v).length < 2) continue;
        const side = reach(v, u); if (side.has(u)) continue;      // in a ring
        const set = side.size * 2 <= lay.atoms.length ? side : new Set(lay.atoms.map(a => a.id).filter(id => !side.has(id)));
        const U = byId.get(u), V = byId.get(v), ax = unit(V.x - U.x, V.y - U.y), out = new Map();
        for (const id of set) { const a = byId.get(id), dx = a.x - U.x, dy = a.y - U.y, t = dx * ax.x + dy * ax.y; out.set(id, { x: U.x + 2 * t * ax.x - dx, y: U.y + 2 * t * ax.y - dy }); }
        const at = id => out.get(id) || byId.get(id);
        if (lay.atoms.some(a => out.has(a.id) && lay.atoms.some(o => !out.has(o.id) && o.id !== u && o.id !== v && dist(at(a.id), o) < 0.55 * L))) continue;
        const d = dist(at(p), at(q));
        if (d < d0 - 0.22 * L && (!best || d < best.d)) best = { d, out };
      }
      if (!best) break;
      best.out.forEach((xy, id) => { const a = byId.get(id); a.x = xy.x; a.y = xy.y; });
      moved = true;
    }
    if (moved) box.refreshDisplay();
    return moved;
  }

  /* ---------- the arrows of a step: [{fish, A, B, Cc, Dd, raw}], each end {sp, box, id} ---------- */
  function resolveArrows(step, species) {
    const find = n => {
      for (let k = 0; k < species.length; k++) {
        const sp = species[k]; if (!sp.box) continue;
        const id = sp.byId ? (sp.box.atomIndex(n) >= 0 ? n : undefined) : sp.byMap.get(n);
        if (id !== undefined) return { sp: k, box: sp.box, id };
      }
      return null;
    };
    const raw = [];
    if (step.graph) {
      for (const ar of step.gArrows || []) {
        const f = ar.from || {}, t = ar.to || {};
        raw.push({ fish: !!ar.half, a: f.kind === 'bond' ? f.a : f.atom, b: f.kind === 'bond' ? f.b : null, c: t.kind === 'bond' ? t.a : t.atom, d: t.kind === 'bond' ? t.b : null, raw: ar.label || JSON.stringify(ar) });
      }
    } else {
      for (const tok of String(step.arrows || '').split(';').map(x => x.trim()).filter(Boolean)) {
        const m = tok.match(/^(~)?(\d+)(?:-(\d+))?>(\d+)(?:-(\d+))?$/);
        raw.push(m ? { fish: !!m[1], a: +m[2], b: m[3] ? +m[3] : null, c: +m[4], d: m[5] ? +m[5] : null, raw: tok } : { bad: true, raw: tok });
      }
    }
    const arrows = [], unresolved = [];
    for (const r of raw) {
      if (r.bad) { unresolved.push(r.raw); continue; }
      const A = find(r.a), B = r.b != null ? find(r.b) : null, Cc = find(r.c), Dd = r.d != null ? find(r.d) : null;
      if (!A || !Cc || (r.b != null && !B) || (r.d != null && !Dd)) { unresolved.push(r.raw); continue; }
      arrows.push(classify({ fish: r.fish, A, B, Cc, Dd, raw: r.raw }));
    }
    // two fishhooks making one bond end side by side; two arrows leaving one bond start on either side of it
    const key = (p, q) => [p.sp + ':' + p.id, q.sp + ':' + q.id].sort().join('|');
    arrows.forEach((ar, i) => arrows.forEach((o, j) => {
      if (j <= i) return;
      if (ar.fish && o.fish && ar.tk === 'bond' && o.tk === 'bond' && key(ar.Cc, ar.Dd) === key(o.Cc, o.Dd) && !ar.mate && !o.mate) { ar.mate = o; o.mate = ar; }
      if (ar.B && o.B && key(ar.A, ar.B) === key(o.A, o.B) && !ar.share && !o.share) { ar.share = o; o.share = ar; }
    }));
    return { arrows, tokens: raw.length, unresolved };
  }
  const same = (p, q) => !!(p && q && p.box === q.box && p.id === q.id);
  const P = r => r.box.atomPoint(r.id);
  const labelOf = r => r.box.atomLabels[r.box.atomIndex(r.id)] || null;
  const rectsOf = r => { const l = labelOf(r); return l ? l.rects : []; };
  const glyphOf = rects => (rects.length ? [rects[0]] : []);            // the element's letters (then the charge, then the H's)
  const bonded = (p, q) => p.box === q.box && p.box.lay.bonds.some(b => (b.a === p.id && b.b === q.id) || (b.a === q.id && b.b === p.id));
  /* what an arrow ends on: 'atom' (a lone pair or a bond's electrons going to an atom), 'own' (a bond's electrons onto one of
     its own atoms), 'pi' (into a bond that is there already), 'bond' (a new bond; X: the atom attacked, the one not giving
     the electrons, Y: its partner) */
  function classify(ar) {
    if (!ar.Dd) { ar.tk = ar.B && (same(ar.Cc, ar.A) || same(ar.Cc, ar.B)) ? 'own' : 'atom'; ar.X = ar.Cc; return ar; }
    if (bonded(ar.Cc, ar.Dd)) { ar.tk = 'pi'; return ar; }
    ar.tk = 'bond';
    const src = [ar.A, ar.B].filter(Boolean), out = [ar.Cc, ar.Dd].filter(r => !src.some(s => same(s, r)));
    if (out.length === 1) { ar.X = out[0]; ar.Y = same(ar.X, ar.Cc) ? ar.Dd : ar.Cc; }
    // a 1,2-shift (C–H moving to the C+ next door): the H and the C+ are about a bond apart, and the middle between them
    // lies on the C–C bond, so the arrow ends at the C+ as for a far bond
    if (ar.X && ar.B && src.some(r => same(r, ar.Y))) { const o = same(ar.A, ar.Y) ? ar.B : ar.A; ar.shift = bonded(o, ar.X); }
    return ar;
  }
  const isFar = ar => ar.tk === 'bond' && !!ar.X && (ar.shift || dist(P(ar.Cc), P(ar.Dd)) > FAR);
  const ends = ar => [ar.A, ar.B, ar.Cc, ar.Dd].filter(Boolean);
  /* where the electrons go, for turning the species and for the lone pair an arrow leaves from */
  function targetPoint(ar) {
    if (ar.tk === 'atom' || ar.tk === 'own') return P(ar.X);
    if (ar.tk === 'bond' && isFar(ar) && !ar.mate) return P(ar.X);
    return mid(P(ar.Cc), P(ar.Dd));
  }
  const sourcePoint = ar => (ar.B ? mid(P(ar.A), P(ar.B)) : P(ar.A));

  /* ---------- the row: species left to right ---------- */
  const molMid = box => (box.h - box.captionHeight()) / 2;
  /* each molecule centred on y = 0, or nudged up/down to meet the species an arrow comes from: returns where the row has
     got to and where the '+' signs go */
  function placeSpecies(species, arrows) {
    let x = 0; const plusAt = [];
    species.forEach((sp, k) => {
      if (k) { plusAt.push(x + GAP / 2 + 8); x += GAP + 16; }
      if (!sp.box) { const t = sp.item || (sp.item = textItem(sp.text, 15)), b = t.bbox(); t.x = x; t.y = -b.h / 2; t.render(); x += b.w; return; }
      const box = sp.box; box.x = x; box.y = -molMid(box);
      const mine = [], theirs = [];
      for (const ar of arrows) {
        const e = ends(ar);
        if (e.some(r => r.sp === k) && e.some(r => r.sp < k)) e.forEach(r => { if (r.sp === k) mine.push(r); else if (r.sp < k) theirs.push(r); });
      }
      if (mine.length && theirs.length) {
        const want = mean(theirs.map(r => P(r).y)), have = mean(mine.map(r => P(r).y));
        const lim = 0.3 * Math.max(box.h, ...theirs.map(r => r.box.h));
        box.y += Math.max(-lim, Math.min(lim, want - have));
      }
      x += box.w;
    });
    return { x, plusAt };
  }
  /* how good a placement of the species is (lower is better): how far the electrons travel between species (source to
     target) and across each new bond that joins two species, an attack from the same side as the group leaving the atom
     attacked (an SN2 is a backside attack), and arrows or new bonds that would cross */
  function layoutCost(arrows) {
    let c = 0; const segs = [];
    // an atom that bonds to (or gives electrons to) another molecule has that molecule in an open gap between its bonds
    // (the middle of a gap of 100° or more: opposite a lone bond, or between two of a C+'s three)
    const facing = (r, toward) => {
      const p = P(r), box = r.box, i = box.atomIndex(r.id), dirs = [];
      box.bondPairs.forEach(([a, b]) => { if (a !== i && b !== i) return; const q = box.atomXY[a === i ? b : a]; dirs.push(norm(Math.atan2(box.y + q.y - p.y, box.x + q.x - p.x))); });
      if (!dirs.length) return 0;
      dirs.sort((x, y) => x - y);
      // a gap that holds the charge ('+' of a C+ written in it) is no way in
      const lab = box.atomLabels[i], q = lab && (box.lay.atoms[i] || {}).charge && lab.rects[1], qa = q ? Math.atan2((q.y0 + q.y1) / 2, (q.x0 + q.x1) / 2) : null;
      const v = Math.atan2(toward.y - p.y, toward.x - p.x); let best = null;
      dirs.forEach((d, k) => { const e = k + 1 < dirs.length ? dirs[k + 1] : dirs[0] + TAU; if (e - d >= 100 * DEG && (qa == null || angDist((d + e) / 2, qa) > 40 * DEG)) { const a = angDist((d + e) / 2, v); if (best == null || a < best) best = a; } });
      return best == null ? 0 : (1 - Math.cos(best)) * 15;
    };
    for (const ar of arrows) {
      const S = sourcePoint(ar), T = targetPoint(ar), apart = ends(ar).some(r => r.sp !== ar.A.sp);
      if (apart) { c += dist(S, T); segs.push([S, T]); }
      if (ar.X && ar.X.sp !== ar.A.sp) c += facing(ar.X, S);
      if (!ar.B && ar.X && ar.X.sp !== ar.A.sp) c += facing(ar.A, P(ar.X));
      if (ar.Dd && ar.Cc.sp !== ar.Dd.sp) c += facing(ar.Cc, P(ar.Dd)) + facing(ar.Dd, P(ar.Cc));
      // one electron of bond Br–Cu going to a new C···Br bond: with Cu straight behind the Br the fishhook has to wrap
      // round the Br (and its lone pairs); a bond at the side gives it a clear path
      const pivot = ar.mate && ar.B && [ar.Cc, ar.Dd].find(r => same(r, ar.A) || same(r, ar.B));
      if (pivot) {
        const o = same(pivot, ar.A) ? ar.B : ar.A, nb = same(pivot, ar.Cc) ? ar.Dd : ar.Cc, p0 = P(pivot);
        const u = unit(P(o).x - p0.x, P(o).y - p0.y), v = unit(P(nb).x - p0.x, P(nb).y - p0.y);
        c += Math.max(0, -(u.x * v.x + u.y * v.y)) * 40;
      }
      if (ar.Dd && ar.Cc.sp !== ar.Dd.sp) { c += dist(P(ar.Cc), P(ar.Dd)); segs.push([P(ar.Cc), P(ar.Dd)]); }
      const X = ar.X; if (!X || X.sp === ar.A.sp) continue;
      for (const lg of arrows) {
        if (lg === ar || lg.tk !== 'own' || !(same(lg.A, X) || same(lg.B, X)) || same(lg.Cc, X)) continue;
        const lb = X.box.lay.bonds.find(b => (b.a === lg.A.id && b.b === lg.B.id) || (b.a === lg.B.id && b.b === lg.A.id));
        if (lb && lb.order > 1) continue;          // a C=O's π electrons going to the O: the O stays, nothing leaves
        const x = P(X), u = unit(S.x - x.x, S.y - x.y), v = unit(P(lg.Cc).x - x.x, P(lg.Cc).y - x.y);
        c += (1 + u.x * v.x + u.y * v.y) * 150;
      }
    }
    for (let i = 0; i < segs.length; i++) for (let j = i + 1; j < segs.length; j++) {
      const [p1, p2] = segs[i], [q1, q2] = segs[j];
      if ([q1, q2].some(q => dist(q, p1) < 2 || dist(q, p2) < 2)) continue;
      if (segsCross(p1, p2, q1, q2)) c += 90;
    }
    return c;
  }
  /* turn each species that an arrow joins to another (30° steps) and, without stereo marks, mirror it, for the lowest
     layoutCost, one species at a time until nothing improves; a turn costs a little, so a molecule stays as it is drawn
     elsewhere unless turning it helps */
  function orient(species, arrows) {
    const movable = species.map((sp, k) => k).filter(k => species[k].box && !species[k].byId && arrows.some(ar => ends(ar).some(r => r.sp === k) && ends(ar).some(r => r.sp !== k))).slice(0, 5);
    if (!movable.length) return;
    const poses = k => { const out = []; for (let a = -150; a <= 180; a += 30) for (const m of hasStereo(species[k].box) ? [false] : [false, true]) out.push({ a, m }); return out; };
    const apply = (k, pose) => {
      const sp = species[k], box = sp.box;
      if (!!sp.mirror !== pose.m) { box.lay.atoms.forEach(a => { a.x = -a.x; }); sp.mirror = pose.m; }
      if (box.angle !== pose.a || sp.mirror !== sp.drawnMirror) { box.angle = pose.a; sp.drawnMirror = sp.mirror; box.refreshDisplay(); }
    };
    const poseCost = pose => (pose.a === 0 ? 0 : pose.a === 180 ? 3 : 12) + (pose.m ? 4 : 0);
    const cur = new Map(movable.map(k => [k, { a: 0, m: false }]));
    const score = () => { placeSpecies(species, arrows); return layoutCost(arrows) + [...cur.values()].reduce((s, p) => s + poseCost(p), 0); };
    movable.forEach(k => apply(k, cur.get(k)));
    let best = score();
    for (let round = 0; round < 3; round++) {
      let better = false;
      for (const k of movable) {
        let pick = cur.get(k);
        for (const pose of poses(k)) {
          cur.set(k, pose); apply(k, pose);
          const c = score();
          if (c < best - 1e-6) { best = c; pick = pose; better = true; }
        }
        cur.set(k, pick); apply(k, pick);
      }
      if (!better) break;
    }
  }
  /* the pairs of atoms in one molecule that the step brings together: both ends of a new bond within it, or two of its
     atoms that both bond to one other molecule (a diene's ends) */
  function togetherPairs(species, arrows) {
    const out = new Map();
    const add = (p, q) => { if (!p.box || p.box !== q.box || p.id === q.id) return; const k = species.findIndex(s => s.box === p.box); (out.get(k) || out.set(k, []).get(k)).push([p.id, q.id]); };
    const news = arrows.filter(ar => ar.tk === 'bond');
    for (const ar of news) if (ar.Cc.sp === ar.Dd.sp) add(ar.Cc, ar.Dd);
    for (let i = 0; i < news.length; i++) for (let j = i + 1; j < news.length; j++) {
      const a = news[i], b = news[j];
      for (const [p, p2] of [[a.Cc, a.Dd], [a.Dd, a.Cc]]) for (const [q, q2] of [[b.Cc, b.Dd], [b.Dd, b.Cc]])
        if (p.sp === q.sp && p2.sp === q2.sp && p.sp !== p2.sp && !same(p, q)) add(p, q);
    }
    return out;
  }

  /* ---------- electrons ---------- */
  /* the bond lines of a box as drawn (each line of a double or triple bond, stopping at the labels), in scene coordinates,
     each with the bond it belongs to: [{P, Q, k}] */
  function inkLines(box) {
    const out = [], pts = box.atomXY, pairs = box.bondPairs;
    const own = (x1, y1, x2, y2) => { let best = -1, bd = 1e9; const m = { x: (x1 + x2) / 2, y: (y1 + y2) / 2 }; pairs.forEach(([i, j], k) => { const d = segDist(m, pts[i], pts[j]); if (d < bd) { bd = d; best = k; } }); return best; };
    box.g.querySelectorAll('line').forEach(l => {
      const x1 = +l.getAttribute('x1'), y1 = +l.getAttribute('y1'), x2 = +l.getAttribute('x2'), y2 = +l.getAttribute('y2');
      if (![x1, y1, x2, y2].every(Number.isFinite)) return;
      out.push({ P: { x: box.x + x1, y: box.y + y1 }, Q: { x: box.x + x2, y: box.y + y2 }, k: own(x1, y1, x2, y2) });
    });
    box.g.querySelectorAll('polygon').forEach(pg => {          // a wedge: its narrow end to the middle of its wide end
      const q = String(pg.getAttribute('points') || '').trim().split(/\s+/).map(s => s.split(',').map(Number));
      if (q.length !== 3) return;
      const a = { x: q[0][0], y: q[0][1] }, b = mid({ x: q[1][0], y: q[1][1] }, { x: q[2][0], y: q[2][1] });
      out.push({ P: { x: box.x + a.x, y: box.y + a.y }, Q: { x: box.x + b.x, y: box.y + b.y }, k: own(a.x, a.y, b.x, b.y) });
    });
    return out;
  }
  function atomEntries(boxes) {                // every atom with lone pairs or an unpaired electron: where it is and what is round it
    const out = new Map();
    boxes.forEach(box => {
      const ink = inkLines(box);
      box.lay.atoms.forEach((atom, i) => {
        let pairs = 0; try { pairs = C.lonePairCount(box.lay, atom) || 0; } catch (e) { pairs = 0; }
        const singles = atom.radical === true ? 1 : Math.max(0, Math.round(+atom.radical || 0));
        const p = box.atomXY[i]; if (!p) return;
        const at = { x: box.x + p.x, y: box.y + p.y }, lab = box.atomLabels[i], rects = lab ? lab.rects : [];
        const segs = ink.filter(l => segDist(at, l.P, l.Q) < 40).map(l => ({ a: { x: l.P.x - at.x, y: l.P.y - at.y }, b: { x: l.Q.x - at.x, y: l.Q.y - at.y } }));
        const occ = [], bondDirs = [];
        box.bondPairs.forEach(([a, b]) => { if (a !== i && b !== i) return; const o = box.atomXY[a === i ? b : a]; const th = Math.atan2(o.y - p.y, o.x - p.x); occ.push(th); bondDirs.push(th); });
        const hAng = lab ? (lab.hAngle != null ? lab.hAngle : lab.hSide ? (lab.hSide > 0 ? 0 : Math.PI) : null) : null;
        if (hAng != null) occ.push(hAng, hAng - 40 * DEG, hAng + 40 * DEG);   // the H's (OH, NH2, H2O): the dots go round the element, not beside the H
        const chg = atom.charge && rects[1], occQ = [];                      // the charge, at a corner of the element
        if (chg && bondDirs.length) occQ.push(Math.atan2((chg.y0 + chg.y1) / 2, (chg.x0 + chg.x1) / 2));
        const inside = [];                       // the inside of a ring is no place for a lone pair
        for (const ring of ringsAt(box, i)) {
          const c = { x: mean(ring.map(k => box.atomXY[k].x)), y: mean(ring.map(k => box.atomXY[k].y)) }, th = Math.atan2(c.y - p.y, c.x - p.x);
          const q1 = box.atomXY[ring[1]], q2 = box.atomXY[ring[ring.length - 1]];
          inside.push({ th, half: angDist(Math.atan2(q1.y - p.y, q1.x - p.x), Math.atan2(q2.y - p.y, q2.x - p.x)) / 2 });
          occ.push(th, th - 25 * DEG, th + 25 * DEG);
        }
        out.set(entryKey({ box, id: atom.id }, boxes), { box, i, id: atom.id, at, rects, segs, occ: occ.concat(occQ), occNoQ: occ, bondDirs, inside, bonded: bondDirs.length > 0 || hAng != null, pairs, singles, fixed: [], extra: [] });
      });
    });
    return out;
  }
  const entryKey = (r, boxes) => r.id + '@' + boxes.indexOf(r.box);
  function ringsAt(box, i) {                     // the smallest ring (atom indexes) through each pair of i's bonds, up to 8 atoms
    const adj = new Map(); box.bondPairs.forEach(([a, b]) => { (adj.get(a) || adj.set(a, []).get(a)).push(b); (adj.get(b) || adj.set(b, []).get(b)).push(a); });
    const ns = adj.get(i) || [], rings = [];
    for (let x = 0; x < ns.length; x++) for (let y = x + 1; y < ns.length; y++) {
      const prev = new Map([[ns[x], -1]]); let q = [ns[x]];
      for (let depth = 0; depth < 7 && q.length && !prev.has(ns[y]); depth++) {
        const nq = []; for (const u of q) for (const v of adj.get(u) || []) if (v !== i && !prev.has(v)) { prev.set(v, u); nq.push(v); }
        q = nq;
      }
      if (!prev.has(ns[y])) continue;
      const ring = [i]; for (let v = ns[y]; v !== -1; v = prev.get(v)) ring.push(v);
      rings.push(ring);
    }
    return rings;
  }
  /* how far out along th a dot (or a pair, set across th) clears the label ink and the bond lines round the atom */
  function clearRadius(e, th, pair) {
    const ux = Math.cos(th), uy = Math.sin(th), offs = pair ? [PAIR_HALF, -PAIR_HALF] : [0];
    for (let R = e.rects.length ? 7 : 9; R < 60; R += 0.5) {
      const ok = offs.every(o => {
        const x = R * ux - o * uy, y = R * uy + o * ux;
        return e.rects.every(r => rectDist(x, y, r) >= DOT_CLEAR) && e.segs.every(s => segDist({ x, y }, s.a, s.b) >= DOT_CLEAR);
      });
      if (ok) return R;
    }
    return 60;
  }
  /* `count` directions spread evenly over the free arcs between the occupied angles; an arc inside a ring counts for
     little, so electrons go there only when everything outside is taken */
  function spread(occ, count, inside = []) {
    if (!count) return [];
    if (!occ.length) return Array.from({ length: count }, (_, k) => k * TAU / count);
    const s = occ.map(norm).sort((a, b) => a - b);
    const gaps = s.map((a, k) => ({ from: a, size: (k + 1 < s.length ? s[k + 1] : s[0] + TAU) - a, n: 0 }));
    gaps.forEach(g => { g.w = inside.some(r => angDist(g.from + g.size / 2, r.th) < r.half) ? 0.15 : 1; });
    for (let k = 0; k < count; k++) gaps.reduce((b, g) => (g.w * g.size / (g.n + 2) > b.w * b.size / (b.n + 2) ? g : b)).n++;
    return gaps.flatMap(g => Array.from({ length: g.n }, (_, k) => g.from + g.size * (k + 1) / (g.n + 1)));
  }
  /* every way of sharing `count` groups among the free arcs between the occupied angles (an arc takes one per 60°), each
     arc's share spaced evenly in it (half a space at its ends: two pairs in one arc sit as far apart as they can); the
     arcs inside a ring only when the outside has no room */
  function spreads(occ, count, inside = []) {
    if (!occ.length) return [Array.from({ length: count }, (_, k) => k * TAU / count)];
    const s = occ.map(norm).sort((a, b) => a - b);
    const gaps = s.map((a, k) => ({ from: a, size: (k + 1 < s.length ? s[k + 1] : s[0] + TAU) - a })).filter(g => g.size > 1e-6);
    gaps.forEach(g => { g.ring = inside.some(r => angDist(g.from + g.size / 2, r.th) < r.half); g.cap = Math.floor(g.size / (60 * DEG)); });
    const room = gs => gs.reduce((t, g) => t + g.cap, 0);
    let open = gaps.filter(g => !g.ring && g.cap > 0);
    if (room(open) < count) open = gaps.filter(g => g.cap > 0);
    const out = [spread(occ, count, inside)];
    if (room(open) < count) return out;
    const rec = (k, left, acc) => {
      if (out.length > 60) return;
      if (k === open.length) { if (!left) out.push(acc.flatMap(([g, n]) => Array.from({ length: n }, (_, j) => g.from + g.size * (j + 0.5) / n))); return; }
      for (let n = Math.min(left, open[k].cap); n >= 0; n--) rec(k + 1, left - n, acc.concat([[open[k], n]]));
    };
    rec(0, count, []);
    return out;
  }
  /* the free direction nearest to `want` (at least 40° from what is occupied, and not pushed far out by a wide label such
     as NH2), for the electrons an arrow leaves from */
  function nearestFree(e, occ, want, pair) {
    let best = null;
    for (let k = 0; k < 72; k++) {
      const th = k * 5 * DEG, clear = occ.length ? Math.min(...occ.map(o => angDist(th, o))) : Math.PI;
      const sc = clear >= 40 * DEG ? angDist(th, want) + Math.max(0, clearRadius(e, th, pair) - 14) * 4 * DEG : 10 + (Math.PI - clear);
      if (!best || sc < best.sc) best = { th, sc };
    }
    return best.th;
  }
  const groupDots = (f, at) => { const cx = at.x + f.R * Math.cos(f.th), cy = at.y + f.R * Math.sin(f.th); return (f.pair ? [PAIR_HALF, -PAIR_HALF] : [0]).map(o => ({ x: cx - o * Math.sin(f.th), y: cy + o * Math.cos(f.th) })); };
  const groupGap = (f, g) => { const o = { x: 0, y: 0 }; let m = 1e9; for (const p of groupDots(f, o)) for (const q of groupDots(g, o)) m = Math.min(m, dist(p, q)); return m; };
  /* the groups of an atom kept apart, so each pair reads as a pair: a group that may still move (not one an arrow leaves)
     goes further out, up to 17 px, until its nearest dot is GROUP_SEP from every other group's. Returns how far they went */
  function separate(groups) {
    let pushed = 0;
    for (let it = 0; it < 40; it++) {
      let moved = false;
      for (let i = 0; i < groups.length; i++) for (let j = i + 1; j < groups.length; j++) {
        if (groupGap(groups[i], groups[j]) >= GROUP_SEP) continue;
        const mv = [groups[i], groups[j]].filter(f => f.rest && f.R < 17);
        mv.forEach(f => { f.R += 0.5; pushed += 0.5; }); if (mv.length) moved = true;
      }
      if (!moved) break;
    }
    return pushed;
  }
  function fixSources(entries, arrows, boxes) {  // phase 1: the electrons arrows leave from, each turned toward its target
    arrows.forEach(ar => {
      if (ar.B) return;
      const e = entries.get(entryKey(ar.A, boxes)); if (!e) return;
      const used = kind => e.fixed.filter(f => f.pair === kind).length;
      // a fishhook leaves the unpaired electron (or one electron of a pair); a full arrow leaves a lone pair
      const freeS = e.singles > used(false), freeP = e.pairs > used(true);
      let pair = ar.fish ? (freeS ? false : freeP ? true : null) : (freeP ? true : freeS ? false : null), extra = false;
      if (pair === null) {
        const own = e.fixed.find(f => f.pair === !ar.fish) || e.fixed[0];
        if (own) { ar.src = own; return; }
        // nothing drawn on the atom (a metal: Pd's d electrons, Na's 3s electron): the electrons the arrow moves are drawn at its tail
        pair = !ar.fish; extra = true;
      }
      const T = targetPoint(ar), want = Math.atan2(T.y - e.at.y, T.x - e.at.x);
      const th = e.bonded ? nearestFree(e, e.occ.concat(e.fixed.map(f => f.th)), want, pair) : want;
      const f = { pair, th, R: clearRadius(e, th, pair), source: true, extra };
      const others = e.fixed.slice();
      for (let k = 0; k < 40 && others.some(g => groupGap(f, g) < GROUP_SEP); k++) f.R += 0.5;   // clear of a pair another arrow leaves
      f.R0 = f.R;
      e.fixed.push(f); ar.src = f; ar.srcEntry = e;
      if (extra) e.extra.push(f);
    });
  }
  function placeRest(entries, arrowCurves) {     // phase 3: every other electron, clear of bonds, labels, the fixed ones and arrows
    const items = [];
    for (const e of entries.values()) {
      const nPairs = Math.max(0, e.pairs - e.fixed.filter(f => f.pair && !f.extra).length), nSingles = Math.max(0, e.singles - e.fixed.filter(f => !f.pair && !f.extra).length);
      const n = nPairs + nSingles;
      const near = [], close = [];               // arrows passing close to this atom: their directions, and their points
      for (const pts of arrowCurves) for (const q of pts) { const d = dist(q, e.at); if (d < 26 && d > 1) near.push(Math.atan2(q.y - e.at.y, q.x - e.at.x)); if (d < 40) close.push(q); }
      if (n) {
        // candidate sets of directions, each judged by how far out its groups have to sit to stay clear and apart, and how
        // near the arrows' lines they come (a crowded atom may put a pair beside an arrow rather than far out)
        const sets = [];
        if (!e.bonded && !e.fixed.length) {      // a lone atom or ion: a symmetric cross (or line), turned to sit best
          for (let phi = 0; phi < TAU / Math.max(n, 2) - 1e-9; phi += 15 * DEG) sets.push({ ds: Array.from({ length: n }, (_, k) => phi + k * TAU / n), bias: phi * 0.5 });
        } else {
          // the charge's corner taken as occupied, or not (a pair may sit by it, a little further out, to keep the pairs
          // apart); the arrows' directions taken as occupied, or not
          const fx = e.fixed.map(f => f.th);
          for (const [occ, bias] of [[e.occ.concat(fx, near), 0], [e.occNoQ.concat(fx, near), 2], [e.occ.concat(fx), 3], [e.occNoQ.concat(fx), 5]])
            for (const base of spreads(occ, n, e.inside)) for (const sh of [0, -12, 12]) sets.push({ ds: base.map(t => t + sh * DEG), bias: bias + Math.abs(sh) * 0.08 });
        }
        let best = null;
        for (const set of sets) {
          const groups = set.ds.map((th, k) => ({ pair: k < nPairs, th, R: clearRadius(e, th, k < nPairs), rest: true }));
          const all = e.fixed.concat(groups), pushed = separate(all);
          let sc = set.bias + groups.reduce((s, g) => s + g.R, 0) + pushed * 0.5;
          for (const g of groups) {
            for (const q of groupDots(g, e.at)) for (const a of close) { const d = dist(q, a); if (d < 8) sc += (8 - d) ** 2 * 1.5; }
            const bclear = e.bondDirs.length ? Math.min(...e.bondDirs.map(b => angDist(b, g.th))) : Math.PI;
            if (bclear < 35 * DEG) sc += (35 * DEG - bclear) / DEG * 0.8;
            sc += e.inside.some(r => angDist(g.th, r.th) < r.half) ? 25 : 0;
            for (const o of all) if (o !== g) { const d = groupGap(g, o); if (d < GROUP_SEP) sc += (GROUP_SEP - d) * 4 + (d < 10 ? 40 : 0); }
          }
          if (!best || sc < best.sc) best = { sc, groups };
        }
        e.fixed.push(...best.groups);
      }
    }
    for (const e of entries.values()) for (const f of e.fixed) {
      const x = e.at.x + f.R * Math.cos(f.th), y = e.at.y + f.R * Math.sin(f.th);
      const el = new D.ElectronItem(f.pair, { x, y }, f.th + Math.PI / 2); el.render(); f.item = el; items.push(el);
    }
    return items;
  }
  const dotPoints = f => (f.pair ? [PAIR_HALF, -PAIR_HALF] : [0]).map(o => ({ x: f.cx - o * Math.sin(f.th), y: f.cy + o * Math.cos(f.th) }));

  /* ---------- one arrow's geometry: start, end and bend ---------- */
  /* the middle of bond p–q as drawn (it stops at the labels), and its direction and normal */
  function bondMid(p, q) {
    const a = P(p), b = P(q), L = dist(a, b) || 1, u = { x: (b.x - a.x) / L, y: (b.y - a.y) / L };
    const ra = rectsOf(p), rb = rectsOf(q), ta = ra.length ? exitDist(ra, u.x, u.y) + 1 : 0, tb = rb.length ? exitDist(rb, -u.x, -u.y) + 1 : 0;
    return { m: along(a, u, (ta + L - tb) / 2), u, n: { x: -u.y, y: u.x } };
  }
  function bondOf(ctx, p, q) { return ctx.bonds.find(b => (same(b.p, p) && same(b.q, q)) || (same(b.p, q) && same(b.q, p))); }
  /* the side of a bond with more room beside its middle, away from the rest of the molecule (pushArrow's hook side); a
     lone diatomic, or no difference: up */
  function awaySide(ctx, p, q) {
    const { m, n } = bondMid(p, q);
    const crowd = sd => {
      const c = along(m, n, sd * 18); let k = 0;
      for (const o of ctx.atoms) { if (same(o.ref, p) || same(o.ref, q)) continue; const d = dist(o.P, c); if (d < 40) k += 40 - d; }
      for (const l of ctx.ink) { const d = segDist(c, l.P, l.Q); if (d < 14) k += 14 - d; }
      return k;
    };
    const a = crowd(1), b = crowd(-1);
    if (Math.abs(a - b) > 3) return a < b ? 1 : -1;
    return n.y > 0 || (Math.abs(n.y) < 1e-6 && n.x > 0) ? -1 : 1;
  }
  /* what is round a point where an arrow starts or ends: other labels, bonds, dots (the ones another arrow leaves from most),
     and the other arrows' ends (a head on another arrow's tail reads as a double-headed arrow) */
  function spotPen(p, ar, ctx, skip) {
    let k = 0;
    for (const o of ctx.atoms) {
      if (skip.atoms.some(r => same(r, o.ref))) continue;
      if (o.rects.length) { const d = Math.min(...o.rects.map(r => rectDist(p.x - o.P.x, p.y - o.P.y, r))); if (d < 6) k += (6 - d) ** 2 * 3; }
      else { const d = dist(p, o.P); if (d < 8) k += (8 - d) ** 2 * 2; }
    }
    for (const l of ctx.ink) { if (skip.bonds.includes(l.bond)) continue; const d = segDist(p, l.P, l.Q); if (d < 5) k += (5 - d) ** 2 * 4; }
    for (const f of ctx.dots) { if (f === ar.src) continue; for (const q of f.pts) { const d = dist(p, q); if (d < 9) k += (9 - d) ** 2 * 2; if (f.source && d < 15) k += (15 - d) ** 2 * 3; } }
    for (const o of ctx.placed) if (o.ar !== ar) for (const q of [o.p1, o.p2]) { const d = dist(p, q); if (d < 16) k += (16 - d) ** 2 * 3; }
    return k;
  }
  /* where an arrow may start: [{p, out (the direction it leaves electrons in), bondDir, side (a shared bond: which side)}] */
  function tailOptions(ar, ctx) {
    const a = P(ar.A);
    if (ar.B) {                                    // on the middle of the bond, on its outermost line on the side the arrow leaves to
      const { m, u, n } = bondMid(ar.A, ar.B), bondDir = Math.atan2(u.y, u.x), b = bondOf(ctx, ar.A, ar.B);
      const offs = b ? ctx.ink.filter(l => l.bond === b).map(l => { const q = mid(l.P, l.Q); return (q.x - m.x) * n.x + (q.y - m.y) * n.y; }) : [0];
      const taken = ar.share && ctx.placed.find(o => o.ar === ar.share), out = [];
      if (offs.length < 2 && !ar.share) return [{ p: m, bondDir, side: 0, n }];   // one line: either side
      // two arrows leaving one bond: each tail also a little toward its own target, so the two do not join into one S
      const T = targetPoint(ar), toward = ar.share ? Math.sign((T.x - m.x) * u.x + (T.y - m.y) * u.y) * SHARE_OFF : 0;
      for (const sd of [1, -1]) {
        if (taken && taken.side === sd) continue;
        const outer = Math.max(0, ...offs.map(o => o * sd));
        out.push({ p: along(along(m, n, sd * (outer + (ar.share ? SHARE_OFF : 0))), u, toward), bondDir, side: sd, n });
      }
      return out;
    }
    if (ar.src) {                                  // the electrons facing the target, or turned up to 80° round the atom (an arrow
      const f = ar.src, e = ar.srcEntry, out = [];  // that would cross another, or a bond, may leave from another side)
      if (f.th0 === undefined) f.th0 = f.th;
      const taken = e ? e.occ.concat(e.fixed.filter(g => g !== f).map(g => g.th)) : [];
      for (const d of [0, -40, 40, -80, 80]) {
        const th = f.th0 + d * DEG;
        if (d && (!e || taken.some(o => angDist(o, th) < 38 * DEG))) continue;
        const R = d && e ? clearRadius(e, th, f.pair) : f.R0 || f.R;
        if (d && e && e.fixed.some(g => g !== f && groupGap({ pair: f.pair, th, R }, g) < GROUP_SEP)) continue;
        out.push({ p: { x: a.x + (R + 4.5) * Math.cos(th), y: a.y + (R + 4.5) * Math.sin(th) }, out: th, src: { th, R }, pen: Math.abs(d) / 40 * 5 });
      }
      return out;
    }
    const T = targetPoint(ar), u = unit(T.x - a.x, T.y - a.y), r = (rectsOf(ar.A).length ? exitDist(rectsOf(ar.A), u.x, u.y) : 4) + 4;
    return [{ p: along(a, u, r), out: Math.atan2(u.y, u.x) }];
  }
  /* where an arrow may end, each with what that spot costs: [{p, pen}] */
  function headOptions(ar, ctx, tail) {
    const s = tail.p;
    const out = [];
    const atomEnds = (X, pref, prefW, hard) => {   // round atom X: at its letters (a bare carbon: by its vertex), clear of its charge, H's and bonds
      const c = P(X), all = rectsOf(X), g = glyphOf(all), rest = all.slice(1), charged = !!(X.box.lay.atoms[X.box.atomIndex(X.id)] || {}).charge, bdirs = ctx.bonds.filter(b => same(b.p, X) || same(b.q, X)).map(b => { const o = same(b.p, X) ? b.Q : b.P; return Math.atan2(o.y - c.y, o.x - c.x); });
      const skip = { atoms: [X], bonds: [] };
      for (let k = 0; k < 24; k++) {
        const th = k * 15 * DEG, w = { x: Math.cos(th), y: Math.sin(th) }, dev = angDist(th, pref);
        if (hard && dev > hard) continue;
        const p = along(c, w, g.length ? exitDist(g, w.x, w.y) + HEAD_GAP : BARE_GAP);
        let pen = prefW * (dev / (30 * DEG)) ** 2;
        // clear of its H's by 10 px, of its charge by 6 (the '+' of a C+ sits where an attack on it comes from)
        rest.forEach((r, j) => { const lim = j === 0 && charged ? 6 : 10, d = rectDist(p.x - c.x, p.y - c.y, r); if (d < lim) pen += (lim - d) ** 2 * 4; });
        for (const b of bdirs) { const a = angDist(th, b); if (a < 40 * DEG) pen += ((40 * DEG - a) / (10 * DEG)) ** 2 * 6; }
        out.push({ p, pen: pen + spotPen(p, ar, ctx, skip), c });
      }
    };
    if (ar.tk === 'atom') { const c = P(ar.X); atomEnds(ar.X, Math.atan2(s.y - c.y, s.x - c.x), 1); return out; }
    if (ar.tk === 'own') {
      const other = same(ar.X, ar.A) ? ar.B : ar.A, { n } = bondMid(ar.A, ar.B), c = P(ar.X), o = P(other);
      // beside the atom on the hook's side, a little toward the bond it comes from (pushArrow)
      const side = tail.side || ar.hookSide || 1, back = unit(o.x - c.x, o.y - c.y), w = unit(n.x * side + 0.45 * back.x, n.y * side + 0.45 * back.y);
      atomEnds(ar.X, Math.atan2(w.y, w.x), 2, 150 * DEG); return out;
    }
    if (ar.tk === 'pi') {                          // the middle of the bond, just off its lines on either side
      const { m, u, n } = bondMid(ar.Cc, ar.Dd), b = bondOf(ctx, ar.Cc, ar.Dd), mine = b ? ctx.ink.filter(l => l.bond === b) : [];
      const skip = { atoms: [ar.Cc, ar.Dd], bonds: b ? [b] : [] };
      for (const sd of [1, -1]) for (const t of [0, -5, 5]) {
        let off = 4.5; const base = along(m, u, t);
        while (off < 16 && mine.some(l => segDist(along(base, n, sd * off), l.P, l.Q) < 4)) off += 1;
        const p = along(base, n, sd * off), facing = (s.x - m.x) * n.x * sd + (s.y - m.y) * n.y * sd > 0 ? 0 : 3;
        out.push({ p, pen: facing + Math.abs(t) * 0.3 + spotPen(p, ar, ctx, skip) });
      }
      return out;
    }
    // a new bond
    const c = P(ar.Cc), d = P(ar.Dd), L = dist(c, d), u = unit(d.x - c.x, d.y - c.y);
    if (ar.mate) {                                 // two fishhooks: each tip FISH_GAP to its own side of where the bond forms
      const m = L > FAR ? mid(c, d) : nearMid(ar, ctx), mine = [ar.A, ar.B].filter(Boolean);
      let sd = mine.some(r => same(r, ar.Cc)) ? -1 : mine.some(r => same(r, ar.Dd)) ? 1 : 0;
      if (!sd) sd = dist(sourcePoint(ar), c) < dist(sourcePoint(ar), d) ? -1 : 1;
      if (sd === (ar.mate.fishSd || 0)) sd = -sd;
      ar.fishSd = sd;
      const p = along(m, u, sd * FISH_GAP);
      out.push({ p, pen: spotPen(p, ar, ctx, { atoms: [ar.Cc, ar.Dd], bonds: [] }) - matePen(p, ar, ctx) });
      return out;
    }
    if (ar.X && isFar(ar)) { const x = P(ar.X), y = P(ar.Y); atomEnds(ar.X, Math.atan2(y.y - x.y, y.x - x.x), 3, 65 * DEG); return out; }
    for (const t of [0.5, 0.42, 0.58, 0.34, 0.66]) {
      const p = { x: c.x + (d.x - c.x) * t, y: c.y + (d.y - c.y) * t };
      out.push({ p, pen: Math.abs(t - 0.5) * 20 + spotPen(p, ar, ctx, { atoms: [ar.Cc, ar.Dd], bonds: [] }) });
    }
    return out;
  }
  /* the mate's tip is meant to be near: what spotPen charged for it is given back */
  function matePen(p, ar, ctx) { const o = ctx.placed.find(x => x.ar === ar.mate); if (!o) return 0; const d = dist(p, o.p2); return d < 16 ? (16 - d) ** 2 * 3 : 0; }
  /* where a new bond between two close atoms forms: its middle, or a little toward one end when another atom is there */
  function nearMid(ar, ctx) {
    const c = P(ar.Cc), d = P(ar.Dd);
    for (const t of [0.5, 0.56, 0.44, 0.62, 0.38, 0.68, 0.32]) {
      const q = { x: c.x + (d.x - c.x) * t, y: c.y + (d.y - c.y) * t };
      if (ctx.atoms.every(o => same(o.ref, ar.Cc) || same(o.ref, ar.Dd) || dist(o.P, q) >= 18)) return q;
    }
    return mid(c, d);
  }
  function arrowGeometry(ar, ctx) {
    if (ar.tk === 'own' && !ar.hookSide) ar.hookSide = awaySide(ctx, ar.A, ar.B);
    const srcBond = ar.B ? bondOf(ctx, ar.A, ar.B) : null, tgtBond = ar.tk === 'pi' ? bondOf(ctx, ar.Cc, ar.Dd) : null;
    // what the curve may pass close to: the atom it ends at, or the ends of the bond it ends on (two fishhooks meeting in
    // the space between two far-apart atoms end at neither)
    const farMate = ar.mate && dist(P(ar.Cc), P(ar.Dd)) > FAR;
    const skip0 = { atoms: farMate ? [] : ar.tk === 'pi' || (ar.tk === 'bond' && !(ar.X && isFar(ar) && !ar.mate)) ? [ar.Cc, ar.Dd] : [ar.X || ar.Cc], bonds: [tgtBond].filter(Boolean) };
    const centre = (() => { const os = ctx.atoms.filter(o => o.sp === ar.A.sp || o.sp === ar.Cc.sp); return { x: mean(os.map(o => o.P.x)), y: mean(os.map(o => o.P.y)) }; })();
    let best = null;
    for (const tail of tailOptions(ar, ctx)) {
      // the source bond's line the tail is on is passed over; its other lines (a C=C's second line) are not
      const skip = Object.assign({}, skip0, { lines: new Set(srcBond ? ctx.ink.filter(l => l.bond === srcBond && segDist(tail.p, l.P, l.Q) < 5.5) : []) });
      const s = tail.p, heads = headOptions(ar, ctx, tail).sort((p, q) => p.pen - q.pen).slice(0, 8);
      let tailPen = 0;                             // a tail beside another arrow's head reads as one arrow with two heads
      for (const o of ctx.placed) if (o.ar !== ar) { const d = dist(s, o.p2); if (d < 16) tailPen += (16 - d) ** 2 * 3; }
      for (const head of heads) {
        // a double bond's arrow leaves from the line on the side its head is
        const across = tail.side && !ar.share && ((head.p.x - s.x) * tail.n.x + (head.p.y - s.y) * tail.n.y) * tail.side < -2 ? 40 : 0;
        const e = head.p, L0 = Math.max(dist(s, e), 1), ch = unit(e.x - s.x, e.y - s.y);
        const mags = [...new Set((L0 < 40 ? [7, 10, 13, 17, 22, 28] : [0.12, 0.2, 0.3, 0.42, 0.58, 0.78].map(f => Math.max(10, Math.min(95, f * L0)))).map(Math.round))];
        for (const side of [1, -1]) for (const mag of mags) {
          const bend = { x: -ch.y * side * mag, y: ch.x * side * mag }, cp = control(s, e, bend);
          let sc = head.pen + (tail.pen || 0) + tailPen + across + arrowCost(s, e, bend, ar, ctx, skip, tail) + 0.08 * mag;   // (a wide loop is longer)
          if (head.c) {                            // a head comes in toward the atom, not sliding along its letters
            const v = unit(e.x - cp.x, e.y - cp.y), r = unit(head.c.x - e.x, head.c.y - e.y), a = Math.acos(Math.max(-1, Math.min(1, v.x * r.x + v.y * r.y)));
            if (a > 55 * DEG) sc += ((a - 55 * DEG) / (10 * DEG)) ** 2 * 4;
          }
          sc += (bend.x * (centre.x - (s.x + e.x) / 2) + bend.y * (centre.y - (s.y + e.y) / 2)) > 0 ? (ar.tk === 'own' ? 6 : 3) : 0;   // rather away from the molecules
          if (tail.side && ar.share) {             // a shared bond's arrow leaves on its own side, not straight on from the other's start
            const toSide = (cp.x - s.x) * tail.n.x + (cp.y - s.y) * tail.n.y; if (toSide * tail.side < 0) sc += 60;
            const o = ctx.placed.find(x => x.ar === ar.share);
            if (o) { const v = unit(cp.x - s.x, cp.y - s.y), w = unit(control(o.p1, o.p2, o.bend).x - o.p1.x, control(o.p1, o.p2, o.bend).y - o.p1.y); if (v.x * w.x + v.y * w.y < -0.8) sc += 25; }
          }
          if (ar.tk === 'own' && !ar.share) { const { n } = bondMid(ar.A, ar.B); if (((cp.x - s.x) * n.x + (cp.y - s.y) * n.y) * ar.hookSide < 0) sc += 8; }
          if (ar.mate) {                           // the two fishhooks on one side of the new bond, tips side by side (not tip to tip)
            const c = P(ar.Cc), d = P(ar.Dd);
            if (Math.sign(cross2(c, d, cp)) !== ar.forceSide) sc += 200;
            const arrive = unit(e.x - cp.x, e.y - cp.y), lu = unit(d.x - c.x, d.y - c.y), a = Math.acos(Math.min(1, Math.abs(arrive.x * lu.x + arrive.y * lu.y)));
            if (a < 40 * DEG) sc += ((40 * DEG - a) / (10 * DEG)) ** 2 * 6;
          }
          if (!best || sc < best.sc) best = { sc, p1: s, p2: e, bend, side: tail.side || 0, src: tail.src };
        }
      }
    }
    return best;
  }
  function arrowCost(s, e, bend, ar, ctx, skip, tail) {
    const c = control(s, e, bend), pts = [];
    for (let i = 0; i <= 24; i++) pts.push(quad(s, c, e, i / 24));
    let k = 0;
    for (let i = 2; i <= 24; i++) {
      const q = pts[i];
      for (const o of ctx.atoms) {
        if (skip.atoms.some(r => same(r, o.ref))) {   // the target's own charge and H's still count, to the tip: the arrow ends at its letters
          for (const r of o.rects.slice(1)) { const d = rectDist(q.x - o.P.x, q.y - o.P.y, r); if (d < 5) k += (5 - d) ** 2 * 6; }
          continue;
        }
        if (i > 22) continue;
        if (o.rects.length) { const d = Math.min(...o.rects.map(r => rectDist(q.x - o.P.x, q.y - o.P.y, r))); if (d < 5) k += (5 - d) ** 2 * 2; }
        else { const d = dist(q, o.P); if (d < 9) k += (9 - d) ** 2; }
      }
      if (i > 22) continue;
      // (the line the tail sits on is passed over only as the arrow leaves it: one that comes back across it is charged)
      for (const l of ctx.ink) { if (skip.bonds.includes(l.bond) || (i < 6 && skip.lines && skip.lines.has(l))) continue; const d = segDist(q, l.P, l.Q); if (d < 6) k += (6 - d) ** 2 * (skip.lines && skip.lines.has(l) ? 4 : 1); }
      for (const f of ctx.dots) { if (f === ar.src) continue; for (const d0 of f.pts) { const d = dist(q, d0); if (d < 7) k += (7 - d) ** 2; } }
      for (const o of ctx.placed) { if (o.ar === ar) continue; for (const p of o.pts) { const d = dist(q, p); if (d < 8) k += (8 - d) ** 2 * 0.5; } }
      // room for the lone pairs still to be drawn round an atom the arrow neither leaves nor ends at
      for (const r of ctx.room) { if ((!ar.B && same(r.ref, ar.A)) || skip.atoms.some(x => same(x, r.ref))) continue; const d = dist(q, r.P); if (d < 26) k += (26 - d) ** 2 * 0.2 * r.n; }
    }
    for (const o of ctx.placed) if (o.ar !== ar && o.ar !== ar.mate) k += crossings(pts, o.pts) * 120;   // arrows that cross cannot be told apart
    const v0 = Math.atan2(c.y - s.y, c.x - s.x);
    if (tail.out != null) { const a = angDist(v0, tail.out); if (a > 75 * DEG) k += ((a - 75 * DEG) / (10 * DEG)) ** 2 * 20; }   // leave the electrons outward, not back over the atom
    if (tail.bondDir != null) { const a = Math.min(angDist(v0, tail.bondDir), angDist(v0, tail.bondDir + Math.PI)); if (a < 30 * DEG) k += ((30 * DEG - a) / (10 * DEG)) ** 2 * 20; }   // not along the bond it leaves
    return k;
  }
  /* every arrow placed in turn, then each placed again knowing where all the others went. Two fishhooks making one bond
     are placed together with the arrows that share their source bonds, on whichever side of the new bond suits them all */
  function placeArrows(arrows, ctx) {
    const drop = ar => { const k = ctx.placed.findIndex(o => o.ar === ar); if (k >= 0) ctx.placed.splice(k, 1); };
    const put = (ar, g) => {
      drop(ar); if (!g) return;
      ctx.placed.push({ ar, p1: g.p1, p2: g.p2, bend: g.bend, side: g.side, pts: curve(g.p1, g.p2, g.bend, 40) });
      if (g.src && ar.src && ar.srcEntry) {         // the electrons it leaves from go where the arrow starts
        const f = ar.src, e = ar.srcEntry; f.th = g.src.th; f.R = g.src.R;
        f.cx = e.at.x + f.R * Math.cos(f.th); f.cy = e.at.y + f.R * Math.sin(f.th); f.pts = dotPoints(f);
      }
    };
    const group = ar => [ar, ar.mate, ar.share, ar.mate.share].filter((x, i, all) => x && all.indexOf(x) === i);
    for (let pass = 0; pass < 2; pass++) for (const ar of arrows) {
      const lead = ar.mate ? (arrows.indexOf(ar.mate) < arrows.indexOf(ar) ? ar.mate : ar) : ar.share && ar.share.mate ? (ar.share.mate && arrows.indexOf(ar.share.mate) < arrows.indexOf(ar.share) ? ar.share.mate : ar.share) : null;
      if (!lead) { put(ar, arrowGeometry(ar, ctx)); continue; }
      if (lead !== ar) continue;                 // placed with its group
      const gr = group(ar); let best = null;
      for (const side of [1, -1]) {
        ar.forceSide = ar.mate.forceSide = side;
        gr.forEach(drop);
        let sc = 0; const gs = [];
        for (const x of gr) { const g = arrowGeometry(x, ctx); if (!g) { sc = Infinity; break; } put(x, g); gs.push([x, g]); sc += g.sc; }
        if (!best || sc < best.sc) best = { sc, gs, side };
      }
      gr.forEach(drop);
      ar.forceSide = ar.mate.forceSide = best.side;
      best.gs.forEach(([x, g]) => put(x, g));
    }
    return arrows.map(ar => ctx.placed.find(o => o.ar === ar));
  }

  /* ---------- one step: a row of items laid out about y = 0, then moved into place ---------- */
  function buildRow(step) {
    const row = [], species = [];
    if (step.graph) species.push(graphSpecies(step.graph));
    else if (step.before) String(step.before).split('.').filter(Boolean).forEach(piece => species.push(smilesSpecies(piece)));
    const res = (step.graph || step.before) ? resolveArrows(step, species) : { arrows: [], tokens: 0, unresolved: [] };
    for (const ar of res.arrows) if (ar.shift && shiftH(ar)) ar.Y.box.refreshDisplay();
    togetherPairs(species, res.arrows).forEach((pairs, k) => { if (species[k] && species[k].box && !species[k].byId) bringTogether(species[k].box, pairs); });
    orient(species, res.arrows);
    let { x, plusAt } = placeSpecies(species, res.arrows);
    species.forEach(sp => { if (sp.box) { sp.box.render(); row.push(sp.box); } else row.push(sp.item); });
    /* the electrons arrows leave from, the curved arrows, then every other electron */
    const startBoxes = species.filter(sp => sp.box).map(sp => sp.box);
    const entries = atomEntries(startBoxes);
    fixSources(entries, res.arrows, startBoxes);
    const ctx = { atoms: [], bonds: [], ink: [], dots: [], placed: [] };
    startBoxes.forEach(box => {
      const sp = species.findIndex(s => s.box === box), ink = inkLines(box);
      box.lay.atoms.forEach((a, i) => ctx.atoms.push({ ref: { box, id: a.id }, box, sp, P: { x: box.x + box.atomXY[i].x, y: box.y + box.atomXY[i].y }, rects: box.atomLabels[i] ? box.atomLabels[i].rects : [] }));
      box.bondPairs.forEach(([i, j], k) => {
        const b = { p: { box, id: box.lay.atoms[i].id }, q: { box, id: box.lay.atoms[j].id }, P: { x: box.x + box.atomXY[i].x, y: box.y + box.atomXY[i].y }, Q: { x: box.x + box.atomXY[j].x, y: box.y + box.atomXY[j].y } };
        ctx.bonds.push(b);
        ink.filter(l => l.k === k).forEach(l => ctx.ink.push({ P: l.P, Q: l.Q, bond: b }));
      });
    });
    for (const e of entries.values()) for (const f of e.fixed) { f.cx = e.at.x + f.R * Math.cos(f.th); f.cy = e.at.y + f.R * Math.sin(f.th); f.pts = dotPoints(f); ctx.dots.push(f); }
    ctx.room = [...entries.values()].map(e => ({ ref: { box: e.box, id: e.id }, P: e.at, n: Math.max(0, e.pairs + e.singles - e.fixed.filter(f => !f.extra).length) })).filter(r => r.n > 0);
    const placed = placeArrows(res.arrows, ctx), arrows = [];
    placed.forEach((g, k) => {
      if (!g) return;
      const ar = res.arrows[k], it = new D.ArrowItem(ar.fish ? 'fishhook' : 'curved', g.p1, g.p2); it.bend = g.bend; it.color = RED; it.render();
      arrows.push(it); ar.item = it;
    });
    const electrons = placeRest(entries, ctx.placed.map(o => o.pts));
    /* each '+' on the middle line, or as little above / below it as keeps it clear of the arrows crossing the gap */
    const dots = electrons.map(e => ({ x: e.x, y: e.y })), heads = ctx.placed.map(o => o.p2), allPts = ctx.placed.flatMap(o => o.pts);
    for (const px of plusAt) {
      let best = null;
      for (const dy of [0, -6, 6, -12, 12, -18, 18, -24, 24, -30, 30, -38, 38, -46, 46, -56, 56, -66, 66]) {
        const q = { x: px, y: dy };
        const clear = Math.min(14, ...allPts.map(p => dist(p, q) - 16), ...heads.map(p => dist(p, q) - 22), ...ctx.atoms.map(o => dist(o.P, q) - 18), ...dots.map(d => dist(d, q) - 12));
        if (!best || clear > best.clear + 1e-6) best = { q, clear };
        if (clear >= 0) break;
      }
      row.push(plusItem(best.q.x, best.q.y));
    }
    /* the reaction arrow, past the start and any arrow that swings out to its right, then what the step makes */
    if (step.smiles) {
      if (species.length) x = Math.max(x, ...allPts.map(p => p.x + 10)) + FWD_GAP;
      const fwd = new D.ArrowItem('forward', { x, y: 0 }, { x: x + FWD_LEN, y: 0 }); fwd.render(); row.push(fwd);
      x += FWD_LEN + FWD_GAP;
      const products = [];
      String(step.smiles).split('.').filter(Boolean).forEach((piece, k) => {
        if (k) { row.push(plusItem(x + GAP / 2 + 8, 0)); x += GAP + 16; }
        const box = moleculeBox(stripMaps(piece));
        if (!box) { const t = textItem(stripMaps(piece), 15), b = t.bbox(); t.x = x; t.y = -b.h / 2; t.render(); row.push(t); x += b.w; return; }
        box.x = x; box.y = -molMid(box); box.render(); row.push(box); products.push(box); x += box.w;
      });
      electrons.push(...placeRest(atomEntries(products), []));
    }
    row.push(...electrons, ...arrows);         // drawn in this order: molecules, text, dots, then the arrows on top
    // what was drawn, for the checks: each arrow with its ends (the electrons at its tail that no atom has of its own too)
    const info = res.arrows.filter(ar => ar.item).map(ar => {
      const b = ar.B ? bondOf(ctx, ar.A, ar.B) : null, far = isFar(ar);
      return { raw: ar.raw, item: ar.item, fish: ar.fish, tk: ar.tk, far, mate: !!ar.mate, share: !!ar.share,
        tailBond: ar.B ? bondMid(ar.A, ar.B).m : null, tailLines: b ? ctx.ink.filter(l => l.bond === b).map(l => [l.P, l.Q]) : null,
        src: ar.src ? { x: ar.src.cx, y: ar.src.cy } : null, X: ar.X ? P(ar.X) : null, Y: ar.Y ? P(ar.Y) : null,
        Xrects: ar.X ? rectsOf(ar.X) : [], newBond: ar.Dd ? [P(ar.Cc), P(ar.Dd)] : null };
    });
    const extras = [];
    for (const e of entries.values()) for (const f of e.extra) extras.push({ box: e.box, id: e.id, pair: f.pair, item: f.item });
    return { row, res, species, info, extras };
  }

  function union(items) {
    let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
    items.forEach(it => { const b = it.bbox(); x0 = Math.min(x0, b.x); y0 = Math.min(y0, b.y); x1 = Math.max(x1, b.x + b.w); y1 = Math.max(y1, b.y + b.h); });
    return { x: x0, y: y0, w: x1 - x0, h: y1 - y0 };
  }

  function showMechanism(payload) {
    payload = payload || {};
    const steps = Array.isArray(payload.steps) ? payload.steps : [];
    // under whatever is on the canvas already (an empty canvas starts at the top left)
    let y = D.items.length ? Math.max(...D.items.map(it => { const b = it.bbox(); return b.y + b.h; })) + 60 : Y0;
    const top = y, added = [], rows = [];
    const put = (it, list) => { (list || added).push(it); return it; };
    const line = (str, size, opts) => { const t = textItem(str, size, opts); t.x = X0; t.y = y; t.render(); y += t.bbox().h; return put(t); };
    line(payload.title || ('Mechanism for ' + (payload.name || 'this reaction')), 22, { bold: true });
    if (payload.status) line('Status: ' + payload.status, 13, { color: /^verified/.test(payload.status) ? '#2b8a3e' : '#9a3412' });
    if (payload.onExample) line("Drawn on the workbook's example molecules (Mechanism sheet); the arrows are the same for your molecules.", 13, { color: '#555' });
    const rx = payload.reaction || {};
    if ((rx.reactants || []).length || (rx.products || []).length) {
      const left = (rx.reactants || []).concat(rx.reagent ? [rx.reagent] : []).map(readable).filter(Boolean).join(' + ');
      const right = (rx.products || []).map(readable).filter(Boolean).join(' + ');
      line(wrap('Your reaction: ' + left + (rx.conditions ? ' (' + rx.conditions + ')' : '') + (right ? ' → ' + right : ''), 110), 14, { color: '#333' });
    }
    y += 18;
    let tokens = 0, drawn = 0; const unresolved = [], arrowInfo = [], extras = [];
    steps.forEach((st, k) => {
      const block = [], no = st.step != null ? st.step : k + 1;
      const title = textItem(`Step ${no}: ${st.type || ''}`.trim(), 16, { bold: true });
      title.x = X0; title.y = y; title.render(); block.push(title);
      let bottom = title.y + title.bbox().h;
      const { row, res, info, extras: ex } = buildRow(st);
      tokens += res.tokens; drawn += res.arrows.filter(ar => ar.item).length;
      res.unresolved.forEach(t => unresolved.push(`step ${no}: ${t}`));
      if (!st.graph && !res.tokens) {          // the Mechanism sheet has no arrows for it yet (said, so it does not look like a slip)
        const t = textItem('The workbook has no arrows for this step yet.', NOTE_PX, { color: '#777' });
        t.x = X0; t.y = bottom + 2; t.render(); block.push(t); bottom = t.y + t.bbox().h;
      }
      if (row.length) {                        // under the title; the first box at the left margin unless an arrow bulges out past it
        const u = union(row), dx = X0 + Math.max(0, -u.x), dy = bottom + 2 - u.y;
        row.forEach(it => it.translate(dx, dy));
        const mv = p => p && { x: p.x + dx, y: p.y + dy };
        info.forEach(a => { a.step = no; a.tailBond = mv(a.tailBond); a.tailLines = a.tailLines && a.tailLines.map(l => l.map(mv)); a.src = mv(a.src); a.X = mv(a.X); a.Y = mv(a.Y); a.newBond = a.newBond && a.newBond.map(mv); });
        block.push(...row); bottom = dy + u.y + u.h;
      }
      arrowInfo.push(...info); extras.push(...ex);
      if (st.note) {
        const note = textItem(wrap(st.note, 100), NOTE_PX, { color: '#444' });
        note.x = X0; note.y = bottom + 2; note.render(); block.push(note); bottom = note.y + note.bbox().h;
      }
      rows.push({ step: no, items: block, box: union(block) });
      added.push(...block);
      y = bottom + BLOCK_GAP;
    });
    if (!steps.length) line('There are no steps to draw for this mechanism.', 14, { color: '#555' });
    // the scene is 3200 × 2200; the longest workbook mechanism (6 steps) fits on an empty canvas, but one drawn under
    // others may not, so the scene grows downward instead of piling the new mechanism on top of what is there
    const bottom = union(added).y + union(added).h;
    D.ensureSceneHeight(bottom + 60);
    D.clearSelection();
    added.forEach(it => D.addItem(it));
    D.setMode('select');
    D.scrollToScene(0, Math.max(0, top - 20));
    D.lastMechanism = { title: payload.title || '', tokens, arrows: drawn, unresolved, rows, items: added, top, bottom, arrowInfo, extras };
    return true;
  }

  D.showMechanism = showMechanism;
})();
