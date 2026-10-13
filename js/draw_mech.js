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

   The arrows follow js/reactions.js pushArrow: an arrow leaves the lone pair (or the unpaired electron) that faces its
   target, or the middle of a bond; it stops short of the target atom's label, or ends at the middle of the new bond; a
   bond's electrons moving onto one of its own atoms curl out to the side away from the rest of the molecule. Each arrow
   bends to whichever side keeps it clear of atoms, bonds, dots, '+' signs and the step's other arrows.
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

  const norm = a => ((a % TAU) + TAU) % TAU;
  const angDist = (p, q) => { const d = Math.abs(norm(p) - norm(q)); return d > Math.PI ? TAU - d : d; };
  const mid = (p, q) => ({ x: (p.x + q.x) / 2, y: (p.y + q.y) / 2 });
  const dist = (p, q) => Math.hypot(p.x - q.x, p.y - q.y);
  const unit = (x, y) => { const l = Math.hypot(x, y) || 1; return { x: x / l, y: y / l }; };
  const mean = xs => xs.reduce((s, x) => s + x, 0) / (xs.length || 1);
  const segDist = (p, a, b) => {
    const vx = b.x - a.x, vy = b.y - a.y, l2 = vx * vx + vy * vy || 1;
    const t = Math.max(0, Math.min(1, ((p.x - a.x) * vx + (p.y - a.y) * vy) / l2));
    return Math.hypot(p.x - (a.x + vx * t), p.y - (a.y + vy * t));
  };
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
      if (t.lay) straighten(t); else t = null;
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
      if (mapped) { box.setGraph(mapped); straighten(box); mapped.atoms.forEach((a, i) => { if (a.mapNum && box.lay.atoms[i]) byMap.set(a.mapNum, box.lay.atoms[i].id); }); }
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
  function graphSpecies(graph) {
    const box = new D.MoleculeBox(); box.setGraph(graph);
    return box.lay ? { box, byId: true } : { text: '(structure)' };
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
      arrows.push({ fish: r.fish, A, B, Cc, Dd, raw: r.raw });
    }
    return { arrows, tokens: raw.length, unresolved };
  }
  const same = (p, q) => p && q && p.box === q.box && p.id === q.id;
  const P = r => r.box.atomPoint(r.id);
  const rectsOf = r => { const l = r.box.atomLabels[r.box.atomIndex(r.id)]; return l ? l.rects : []; };

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
        const ends = [ar.A, ar.B, ar.Cc, ar.Dd].filter(Boolean);
        if (ends.some(r => r.sp === k) && ends.some(r => r.sp < k)) ends.forEach(r => { if (r.sp === k) mine.push(r); else if (r.sp < k) theirs.push(r); });
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
  /* how far the electrons travel between species: from the source (atom or bond middle) to the target, and across each
     new bond that joins two species */
  function linkLength(arrows) {
    const pt = (a, b) => (b ? mid(P(a), P(b)) : P(a));
    let c = 0;
    for (const ar of arrows) {
      if (ar.A.sp !== ar.Cc.sp || (ar.Dd && ar.Dd.sp !== ar.A.sp)) c += dist(pt(ar.A, ar.B), pt(ar.Cc, ar.Dd));
      if (ar.Dd && ar.Cc.sp !== ar.Dd.sp) c += dist(P(ar.Cc), P(ar.Dd));
    }
    return c;
  }
  /* turn the species (0° or 180°: a turn keeps R/S) so the atoms the arrows join across species are as close as they can
     be: every combination is tried (a step has at most four species) */
  function orient(species, arrows) {
    const turn = species.map((sp, k) => k).filter(k => species[k].box && !species[k].byId && arrows.some(ar => [ar.A, ar.B, ar.Cc, ar.Dd].some(r => r && r.sp === k) && [ar.A, ar.B, ar.Cc, ar.Dd].some(r => r && r.sp !== k))).slice(0, 5);
    let best = null;
    for (let mask = 0; mask < (1 << turn.length); mask++) {
      turn.forEach((k, j) => { const ang = mask & (1 << j) ? 180 : 0, box = species[k].box; if (box.angle !== ang) { box.angle = ang; box.refreshDisplay(); } });
      placeSpecies(species, arrows);
      const c = linkLength(arrows) + mask * 1e-6;
      if (!best || c < best.c) best = { c, mask };
    }
    if (best) turn.forEach((k, j) => { const ang = best.mask & (1 << j) ? 180 : 0, box = species[k].box; if (box.angle !== ang) { box.angle = ang; box.refreshDisplay(); } });
  }

  /* ---------- electrons ---------- */
  function atomEntries(boxes) {                // every atom with lone pairs or an unpaired electron: where it is and what is round it
    const out = new Map();
    boxes.forEach(box => {
      box.lay.atoms.forEach((atom, i) => {
        let pairs = 0; try { pairs = C.lonePairCount(box.lay, atom) || 0; } catch (e) { pairs = 0; }
        const singles = atom.radical === true ? 1 : Math.max(0, Math.round(+atom.radical || 0));
        const p = box.atomXY[i]; if (!p) return;
        const at = { x: box.x + p.x, y: box.y + p.y }, lab = box.atomLabels[i], rects = lab ? lab.rects : [];
        const segs = [], occ = [];
        box.bondPairs.forEach(([a, b]) => { if (a !== i && b !== i) return; const o = box.atomXY[a === i ? b : a]; segs.push({ x: o.x - p.x, y: o.y - p.y }); occ.push(Math.atan2(o.y - p.y, o.x - p.x)); });
        const hSide = lab ? lab.hSide : 0;
        if (hSide) occ.push(hSide > 0 ? 0 : Math.PI);                       // the H label (OH, NH2, H2O)
        const chg = atom.charge && rects[1];                                 // the charge, up to the right
        if (chg && segs.length) occ.push(Math.atan2((chg.y0 + chg.y1) / 2, (chg.x0 + chg.x1) / 2));
        const inside = [];                       // the inside of a ring is no place for a lone pair
        for (const ring of ringsAt(box, i)) {
          const c = { x: mean(ring.map(k => box.atomXY[k].x)), y: mean(ring.map(k => box.atomXY[k].y)) }, th = Math.atan2(c.y - p.y, c.x - p.x);
          const q1 = box.atomXY[ring[1]], q2 = box.atomXY[ring[ring.length - 1]];
          inside.push({ th, half: angDist(Math.atan2(q1.y - p.y, q1.x - p.x), Math.atan2(q2.y - p.y, q2.x - p.x)) / 2 });
          occ.push(th, th - 25 * DEG, th + 25 * DEG);
        }
        out.set(entryKey({ box, id: atom.id }, boxes), { box, i, id: atom.id, at, rects, segs, occ, inside, bonded: segs.length > 0 || !!hSide, pairs, singles, fixed: [] });
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
  /* how far out along th a dot (or a pair, set across th) clears the label ink and the atom's own bonds */
  function clearRadius(e, th, pair) {
    const ux = Math.cos(th), uy = Math.sin(th), offs = pair ? [PAIR_HALF, -PAIR_HALF] : [0];
    for (let R = e.rects.length ? 7 : 9; R < 60; R += 0.5) {
      const ok = offs.every(o => {
        const x = R * ux - o * uy, y = R * uy + o * ux;
        return e.rects.every(r => rectDist(x, y, r) >= DOT_CLEAR) && e.segs.every(s => segDist({ x, y }, { x: 0, y: 0 }, s) >= DOT_CLEAR);
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
  function fixSources(entries, arrows, boxes) {  // phase 1: the electrons arrows leave from, each turned toward its target
    arrows.forEach(ar => {
      if (ar.B) return;
      const e = entries.get(entryKey(ar.A, boxes)); if (!e) return;
      const used = kind => e.fixed.filter(f => f.pair === kind).length;
      // a fishhook leaves the unpaired electron (or one electron of a pair); a full arrow leaves a lone pair
      const freeS = e.singles > used(false), freeP = e.pairs > used(true);
      const pair = ar.fish ? (freeS ? false : freeP ? true : null) : (freeP ? true : freeS ? false : null);
      if (pair === null) {                     // nothing free there (Na, Pd, a pair already used): share one, or start beside the atom
        const own = e.fixed.find(f => f.pair === !ar.fish) || e.fixed[0];
        if (own) ar.src = own;
        return;
      }
      const T = ar.Dd ? mid(P(ar.Cc), P(ar.Dd)) : P(ar.Cc), want = Math.atan2(T.y - e.at.y, T.x - e.at.x);
      const th = e.bonded ? nearestFree(e, e.occ.concat(e.fixed.map(f => f.th)), want, pair) : want;
      const f = { pair, th, R: clearRadius(e, th, pair) };
      e.fixed.push(f); ar.src = f; ar.srcEntry = e;
    });
  }
  function placeRest(entries, arrowCurves) {     // phase 3: every other electron, clear of bonds, labels, the fixed ones and arrows
    const items = [];
    for (const e of entries.values()) {
      const nPairs = Math.max(0, e.pairs - e.fixed.filter(f => f.pair).length), nSingles = Math.max(0, e.singles - e.fixed.filter(f => !f.pair).length);
      const n = nPairs + nSingles;
      const near = [];                           // arrows passing close to this atom take those directions
      for (const pts of arrowCurves) for (const q of pts) { const d = dist(q, e.at); if (d < 26 && d > 1) near.push(Math.atan2(q.y - e.at.y, q.x - e.at.x)); }
      let dirs;
      if (!n) dirs = [];
      else if (!e.bonded && !e.fixed.length) {     // a lone atom or ion: a symmetric cross (or line), turned to sit best
        let best = null;
        for (let phi = 0; phi < TAU / Math.max(n, 2) - 1e-9; phi += 15 * DEG) {
          const ds = Array.from({ length: n }, (_, k) => phi + k * TAU / n);
          const sc = ds.reduce((s, th, k) => s + clearRadius(e, th, k < nPairs) + near.filter(a => angDist(a, th) < 35 * DEG).length * 30, 0) + phi * 0.5;
          if (!best || sc < best.sc) best = { ds, sc };
        }
        dirs = best.ds;
      } else dirs = spread(e.occ.concat(e.fixed.map(f => f.th), near), n, e.inside);
      dirs.forEach((th, k) => { const pair = k < nPairs; e.fixed.push({ pair, th, R: clearRadius(e, th, pair), rest: true }); });
    }
    for (const e of entries.values()) for (const f of e.fixed) {
      const x = e.at.x + f.R * Math.cos(f.th), y = e.at.y + f.R * Math.sin(f.th);
      const el = new D.ElectronItem(f.pair, { x, y }, f.th + Math.PI / 2); el.render(); f.item = el; items.push(el);
    }
    return items;
  }
  const dotPoints = f => (f.pair ? [PAIR_HALF, -PAIR_HALF] : [0]).map(o => ({ x: f.cx - o * Math.sin(f.th), y: f.cy + o * Math.cos(f.th) }));

  /* ---------- one arrow's geometry: start, end and bend ---------- */
  function arrowGeometry(ar, ctx) {
    const a = P(ar.A);
    // a bond's electrons onto one of its own atoms: a short hook out to the side away from the rest of the molecule (pushArrow)
    if (ar.B && !ar.Dd && (same(ar.Cc, ar.A) || same(ar.Cc, ar.B))) {
      const b = P(ar.B), c = P(ar.Cc), s = mid(a, b), u = unit(b.x - a.x, b.y - a.y);
      let n = { x: -u.y, y: u.x };
      const near = ctx.atoms.filter(o => o.box === ar.A.box && !same(o.ref, ar.A) && !same(o.ref, ar.B) && dist(o.P, s) < 70);
      if (near.length) { const cx = mean(near.map(o => o.P.x)) - s.x, cy = mean(near.map(o => o.P.y)) - s.y; if (cx * n.x + cy * n.y > 0) n = { x: -n.x, y: -n.y }; }
      else if (n.y > 0 || (Math.abs(n.y) < 1e-6 && n.x > 0)) n = { x: -n.x, y: -n.y };   // a lone diatomic: the hook goes up
      const st = { x: s.x + n.x * 6, y: s.y + n.y * 6 };
      const w = unit(n.x + 0.45 * (s.x - c.x) / (dist(s, c) || 1), n.y + 0.45 * (s.y - c.y) / (dist(s, c) || 1)), rc = rectsOf(ar.Cc);
      const r = rc.length ? exitDist(rc, w.x, w.y) + 5 : 11;   // beside the atom, just clear of its label, on the hook's side
      const en = { x: c.x + w.x * r, y: c.y + w.y * r };
      const ch = unit(en.x - st.x, en.y - st.y); let m = { x: -ch.y, y: ch.x };
      if (m.x * n.x + m.y * n.y < 0) m = { x: -m.x, y: -m.y };
      return { p1: st, p2: en, bend: { x: m.x * 10, y: m.y * 10 } };
    }
    // the start: the electrons that face the target, the middle of the bond, or (no electrons there) beside the atom
    const T0 = ar.Dd ? bondTarget(ar, ctx) : P(ar.Cc);
    let s, out = null, bondDir = null;
    if (ar.B) { const b = P(ar.B); s = mid(a, b); bondDir = Math.atan2(b.y - a.y, b.x - a.x); }
    else if (ar.src) { out = ar.src.th; const R = ar.src.R + (ar.src.pair ? 4.5 : 4.5); s = { x: a.x + R * Math.cos(out), y: a.y + R * Math.sin(out) }; }
    else { const u = unit(T0.x - a.x, T0.y - a.y); out = Math.atan2(u.y, u.x); const r = (rectsOf(ar.A).length ? exitDist(rectsOf(ar.A), u.x, u.y) : 4) + 4; s = { x: a.x + u.x * r, y: a.y + u.y * r }; }
    const skipAtoms = ar.Dd ? [] : [ar.Cc], skipBonds = [[ar.A, ar.B], [ar.Cc, ar.Dd]].filter(([p, q]) => p && q);
    const centre = { x: mean(ctx.atoms.filter(o => o.sp === ar.A.sp || o.sp === ar.Cc.sp).map(o => o.P.x)), y: mean(ctx.atoms.filter(o => o.sp === ar.A.sp || o.sp === ar.Cc.sp).map(o => o.P.y)) };
    const L0 = Math.max(dist(s, T0), 1);
    const mags = [...new Set((L0 < 40 ? [10, 13, 17, 22, 28] : [0.2, 0.3, 0.42, 0.58, 0.78].map(f => Math.max(10, Math.min(85, f * L0)))).map(Math.round))];
    let best = null;
    for (const side of [1, -1]) for (const mag of mags) {
      let e = T0;
      const bendFor = endPt => { const u = unit(endPt.x - s.x, endPt.y - s.y); return { x: -u.y * side * mag, y: u.x * side * mag }; };
      let bend = bendFor(e);
      if (!ar.Dd) {                            // stop short of the target's label, along the way the arrow comes in
        const c = P(ar.Cc), ctl = control(s, c, bend), u = unit(ctl.x - c.x, ctl.y - c.y), rects = rectsOf(ar.Cc);
        const r = rects.length ? exitDist(rects, u.x, u.y) + 4 : 8;
        e = { x: c.x + u.x * r, y: c.y + u.y * r }; bend = bendFor(e);
      }
      const sc = arrowCost(s, e, bend, ar, ctx, skipAtoms, skipBonds, out, bondDir) + 0.03 * mag +
        (((bend.x * (centre.x - (s.x + e.x) / 2) + bend.y * (centre.y - (s.y + e.y) / 2)) > 0) ? 3 : 0);   // rather away from the molecules
      if (!best || sc < best.sc) best = { sc, p1: s, p2: e, bend };
    }
    return best;
  }
  /* where an arrow forming a new bond c–d ends: its middle, or a little toward one end when the middle falls on another
     atom (two species side by side can put a third atom there) */
  function bondTarget(ar, ctx) {
    const c = P(ar.Cc), d = P(ar.Dd);
    for (const t of [0.5, 0.56, 0.44, 0.62, 0.38, 0.68, 0.32]) {
      const q = { x: c.x + (d.x - c.x) * t, y: c.y + (d.y - c.y) * t };
      if (ctx.atoms.every(o => same(o.ref, ar.Cc) || same(o.ref, ar.Dd) || dist(o.P, q) >= 18)) return q;
    }
    return mid(c, d);
  }
  function arrowCost(s, e, bend, ar, ctx, skipAtoms, skipBonds, out, bondDir) {
    const c = control(s, e, bend);
    let k = 0;
    for (let i = 2; i <= 22; i++) {
      const q = quad(s, c, e, i / 24);
      for (const o of ctx.atoms) {
        if (skipAtoms.some(r => same(r, o.ref))) continue;
        if (o.rects.length) { const d = Math.min(...o.rects.map(r => rectDist(q.x - o.P.x, q.y - o.P.y, r))); if (d < 5) k += (5 - d) ** 2 * 2; }
        else { const d = dist(q, o.P); if (d < 9) k += (9 - d) ** 2; }
      }
      for (const b of ctx.bonds) {
        if (skipBonds.some(([p, r]) => (same(p, b.p) && same(r, b.q)) || (same(p, b.q) && same(r, b.p)))) continue;
        const cl = b.clear || 7, d = segDist(q, b.P, b.Q); if (d < cl) k += (cl - d) ** 2;
      }
      for (const f of ctx.dots) { if (f === ar.src) continue; for (const d0 of f.pts) { const d = dist(q, d0); if (d < 7) k += (7 - d) ** 2; } }
      for (const pts of ctx.curves) for (const p of pts) { const d = dist(q, p); if (d < 8) k += (8 - d) ** 2 * 0.5; }
    }
    const v0 = Math.atan2(c.y - s.y, c.x - s.x);
    if (out != null) { const a = angDist(v0, out); if (a > 75 * DEG) k += ((a - 75 * DEG) / (10 * DEG)) ** 2 * 20; }   // leave the electrons outward, not back over the atom
    if (bondDir != null) { const a = Math.min(angDist(v0, bondDir), angDist(v0, bondDir + Math.PI)); if (a < 30 * DEG) k += ((30 * DEG - a) / (10 * DEG)) ** 2 * 20; }   // not along the bond it leaves
    return k;
  }

  /* ---------- one step: a row of items laid out about y = 0, then moved into place ---------- */
  function buildRow(step) {
    const row = [], species = [];
    if (step.graph) species.push(graphSpecies(step.graph));
    else if (step.before) String(step.before).split('.').filter(Boolean).forEach(piece => species.push(smilesSpecies(piece)));
    const res = (step.graph || step.before) ? resolveArrows(step, species) : { arrows: [], tokens: 0, unresolved: [] };
    orient(species, res.arrows);
    let { x, plusAt } = placeSpecies(species, res.arrows);
    species.forEach(sp => { if (sp.box) { sp.box.render(); row.push(sp.box); } else row.push(sp.item); });
    /* the electrons arrows leave from, the curved arrows, then every other electron */
    const startBoxes = species.filter(sp => sp.box).map(sp => sp.box);
    const entries = atomEntries(startBoxes);
    fixSources(entries, res.arrows, startBoxes);
    const ctx = { atoms: [], bonds: [], dots: [], curves: [] };
    startBoxes.forEach(box => {
      const sp = species.findIndex(s => s.box === box);
      box.lay.atoms.forEach((a, i) => ctx.atoms.push({ ref: { box, id: a.id }, box, sp, P: { x: box.x + box.atomXY[i].x, y: box.y + box.atomXY[i].y }, rects: box.atomLabels[i] ? box.atomLabels[i].rects : [] }));
      box.bondPairs.forEach(([i, j]) => ctx.bonds.push({ p: { box, id: box.lay.atoms[i].id }, q: { box, id: box.lay.atoms[j].id }, P: { x: box.x + box.atomXY[i].x, y: box.y + box.atomXY[i].y }, Q: { x: box.x + box.atomXY[j].x, y: box.y + box.atomXY[j].y } }));
    });
    for (const e of entries.values()) for (const f of e.fixed) { f.cx = e.at.x + f.R * Math.cos(f.th); f.cy = e.at.y + f.R * Math.sin(f.th); f.pts = dotPoints(f); ctx.dots.push(f); }
    const arrows = [];
    for (const ar of res.arrows) {
      const g = arrowGeometry(ar, ctx);
      const it = new D.ArrowItem(ar.fish ? 'fishhook' : 'curved', g.p1, g.p2); it.bend = g.bend; it.color = RED; it.render();
      arrows.push(it); ctx.curves.push(curve(g.p1, g.p2, g.bend, 40));
    }
    const electrons = placeRest(entries, ctx.curves);
    /* each '+' on the middle line, or as little above / below it as keeps it clear of the arrows crossing the gap */
    const dots = electrons.map(e => ({ x: e.x, y: e.y }));
    for (const px of plusAt) {
      let best = null;
      for (const dy of [0, -6, 6, -12, 12, -18, 18, -24, 24, -30, 30, -38, 38, -46, 46, -56, 56]) {
        const q = { x: px, y: dy };
        const clear = Math.min(14, ...ctx.curves.flat().map(p => dist(p, q) - 14), ...ctx.atoms.map(o => dist(o.P, q) - 18), ...dots.map(d => dist(d, q) - 12));
        if (!best || clear > best.clear + 1e-6) best = { q, clear };
        if (clear >= 0) break;
      }
      row.push(plusItem(best.q.x, best.q.y));
    }
    /* the reaction arrow, past the start and any arrow that swings out to its right, then what the step makes */
    if (step.smiles) {
      if (species.length) x = Math.max(x, ...ctx.curves.flat().map(p => p.x + 10)) + FWD_GAP;
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
    return { row, res, species };
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
      const left = (rx.reactants || []).concat(rx.reagent ? [rx.reagent] : []).join(' + ');
      const right = (rx.products || []).join(' + ');
      line(wrap('Your reaction: ' + left + (rx.conditions ? ' (' + rx.conditions + ')' : '') + (right ? ' → ' + right : ''), 110), 14, { color: '#333' });
    }
    y += 18;
    let tokens = 0, drawn = 0; const unresolved = [];
    steps.forEach((st, k) => {
      const block = [];
      const title = textItem(`Step ${st.step != null ? st.step : k + 1}: ${st.type || ''}`.trim(), 16, { bold: true });
      title.x = X0; title.y = y; title.render(); block.push(title);
      const { row, res } = buildRow(st);
      tokens += res.tokens; drawn += res.arrows.length;
      res.unresolved.forEach(t => unresolved.push(`step ${st.step != null ? st.step : k + 1}: ${t}`));
      let bottom = title.y + title.bbox().h;
      if (row.length) {                        // under the title; the first box at the left margin unless an arrow bulges out past it
        const u = union(row), dx = X0 + Math.max(0, -u.x), dy = bottom + 2 - u.y;
        row.forEach(it => it.translate(dx, dy));
        block.push(...row); bottom = dy + u.y + u.h;
      }
      if (st.note) {
        const note = textItem(wrap(st.note, 100), NOTE_PX, { color: '#444' });
        note.x = X0; note.y = bottom + 2; note.render(); block.push(note); bottom = note.y + note.bbox().h;
      }
      rows.push({ step: st.step != null ? st.step : k + 1, items: block, box: union(block) });
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
    D.lastMechanism = { title: payload.title || '', tokens, arrows: drawn, unresolved, rows, items: added, top, bottom };
    return true;
  }

  D.showMechanism = showMechanism;
})();
