/* Sketcher — a 2D molecule editor on an SVG, with the tool set of the PubChem sketcher.

   The drawing IS a Chem graph (js/chem.js), so naming and SMILES need no conversion:
     atoms [{id, element, charge, radical?, x, y}]
     bonds [{a, b, order, stereo?: 'wedge'|'hash', narrow?: atomId}]   narrow = the stereocentre end

   Left palette
     select    drag an atom (or a marquee-selected group) to move it; Delete removes the selection
     rotate    drag to turn the selection (or the whole drawing) about its centre
     erase     click an atom or a bond
     element   click empty space → place; click an atom → relabel;
               drag → bonded atoms of that element: a short drag one, a longer drag a zigzag chain,
               starting from the atom pressed or from empty space; release on an atom to bond to it
     table     the periodic table; picks any element the namer knows
     chain     drag → zigzag carbon chain
     bonds     click a bond → set its type; drag from an atom → new carbon with that bond;
               click empty space → a fresh C–C; wedge/hash clicked twice flips the narrow end
     charge    click an atom → +1 / −1;  radical: click an atom → unpaired electron
     rings     click empty space → ring; click an atom → ring hanging off it; click a bond → fused ring
     templates a window of common ring systems and molecules to paste

   Fragment row (Et … tBu, Ph, CHO, CO2H, NO2, SO3H): click an atom to attach, empty space to place.
   Toolbar: new, undo/redo (Ctrl+Z / Ctrl+Y), clean, mirror, add/strip hydrogens, import (SMILES,
   IUPAC name or a MOL/SDF/SMILES file), export (SMILES / name / formula / MOL / SVG / PNG).
   Keyboard while hovering: an element letter relabels the atom (c n o s p f i h, l = Cl, b = Br);
   1 2 3 set a bond's order.                                                                      */
(() => {
  'use strict';
  const C = window.Chem;
  const svg = document.getElementById('sketchSvg');
  const info = document.getElementById('sketchInfo');
  const grid = document.getElementById('toolGrid');
  const fragBar = document.getElementById('fragBar');
  const bar = document.getElementById('sketchBar');
  const main = document.querySelector('.sketch-main');
  const importBox = document.getElementById('sketchImport');
  const fileBox = document.getElementById('sketchFile');
  const formatBox = document.getElementById('sketchFormat');
  if (!svg || !C) return;

  const NS = 'http://www.w3.org/2000/svg';
  const BOND = 48;                     // standard bond length, px (two grid cells)
  const SNAP = 24;                     // grid spacing for free-placed atoms
  const ATOM_R = 12, BOND_R = 7;       // hit radii
  const DEG30 = Math.PI / 6;

  let g = { atoms: [], bonds: [], nextId: 1 };
  const undoStack = [], redoStack = [];
  let tool = null;
  let selected = new Set();
  let hover = { atom: null, bond: null };
  let popup = null;                    // the open periodic-table / templates panel, if any

  /* ================= model ================= */
  const atom = id => g.atoms.find(a => a.id === id);
  const bondsOf = id => g.bonds.filter(b => b.a === id || b.b === id);
  const otherEnd = (b, id) => (b.a === id ? b.b : b.a);
  const neighboursIn = (G, id) => G.bonds.filter(b => b.a === id || b.b === id).map(b => G.atoms.find(a => a.id === otherEnd(b, id)));
  const neighbours = id => neighboursIn(g, id);
  const bondBetween = (p, q) => g.bonds.find(b => (b.a === p && b.b === q) || (b.a === q && b.b === p));

  function snapshot() { undoStack.push(JSON.stringify(g)); if (undoStack.length > 200) undoStack.shift(); redoStack.length = 0; }
  /* An edit that would give an atom more bonds than it can have is taken back, with a message, when the drawing was
     valid before it: one stray click with a bond tool on a full atom (an alkyne carbon, say) must not quietly turn
     the molecule into another compound. A drawing that is already invalid can still be edited, so it can be fixed. */
  let notice = '';                     // shown once in the readout below the canvas
  const valenceError = gr => { try { C.validateValences(gr); return ''; } catch (e) { return e.message; } };
  function guarded(edit) {
    const wasValid = !valenceError(g), depth = undoStack.length;
    edit();
    if (!wasValid || undoStack.length <= depth) return;
    const err = valenceError(g);
    if (!err) return;
    g = JSON.parse(undoStack[depth]);                  // back to the drawing before this edit
    undoStack.length = depth;
    selected.clear();
    hover = { atom: null, bond: null };               // those pointed into the drawing that was just replaced
    notice = 'Not done — that edit would leave an atom with too many bonds (' + err.replace(/\s*— that atom cannot exist as written\.?$/, '') +
      '). To change the bonding there, first make room: for example, turn a triple bond into a single bond.';
    render();
  }
  function undo() { if (!undoStack.length) return; redoStack.push(JSON.stringify(g)); g = JSON.parse(undoStack.pop()); selected.clear(); render(); }
  function redo() { if (!redoStack.length) return; undoStack.push(JSON.stringify(g)); g = JSON.parse(redoStack.pop()); selected.clear(); render(); }

  function addAtom(element, x, y) { const a = { id: g.nextId++, element, charge: 0, x, y }; g.atoms.push(a); return a; }
  /* a new end atom for a bond from `from`: an atom already at that spot is reused, never a second one stacked on it
     (dragging a bond along a ring bond and letting go just short of the neighbour snaps exactly onto the neighbour) */
  function atomAtOrNew(from, x, y) {
    const hit = atomAt(x, y);
    return hit && hit !== from ? hit : addAtom('C', x, y);
  }
  function addBond(p, q, order = 1, stereo = null) {
    if (p === q) return null;
    let b = bondBetween(p, q);
    if (!b) { b = { a: p, b: q, order }; g.bonds.push(b); }
    else b.order = order;
    if (stereo) { b.stereo = stereo; b.narrow = p; } else { delete b.stereo; delete b.narrow; }
    return b;
  }
  function removeAtom(id) {
    g.bonds = g.bonds.filter(b => b.a !== id && b.b !== id);
    g.atoms.splice(g.atoms.findIndex(a => a.id === id), 1);
    selected.delete(id);
  }
  function removeBond(b) { g.bonds.splice(g.bonds.indexOf(b), 1); }
  /* the atoms an edit-the-whole-thing action applies to: the selection, else everything */
  const targetAtoms = () => (selected.size > 1 ? [...selected].map(atom) : g.atoms);
  const centroid = atoms => ({ x: atoms.reduce((s, a) => s + a.x, 0) / atoms.length, y: atoms.reduce((s, a) => s + a.y, 0) / atoms.length });

  /* ================= geometry ================= */
  const dist = (x1, y1, x2, y2) => Math.hypot(x2 - x1, y2 - y1);
  const snapGrid = v => Math.round(v / SNAP) * SNAP;
  const polar = (x, y, ang, r) => ({ x: x + r * Math.cos(ang), y: y + r * Math.sin(ang) });

  function atomAt(x, y) {
    let best = null, bd = ATOM_R;
    for (const a of g.atoms) { const d = dist(x, y, a.x, a.y); if (d < bd) { bd = d; best = a; } }
    return best;
  }
  function bondAt(x, y) {
    let best = null, bd = BOND_R;
    for (const b of g.bonds) {
      const p = atom(b.a), q = atom(b.b);
      const L2 = (q.x - p.x) ** 2 + (q.y - p.y) ** 2 || 1;
      const t = Math.max(0.15, Math.min(0.85, ((x - p.x) * (q.x - p.x) + (y - p.y) * (q.y - p.y)) / L2));
      const d = dist(x, y, p.x + t * (q.x - p.x), p.y + t * (q.y - p.y));
      if (d < bd) { bd = d; best = b; }
    }
    return best;
  }

  /* end point of a bond dragged out of `a` towards (x, y): standard length, angle snapped to 30° */
  function snapDir(a, x, y) {
    const ang = Math.round(Math.atan2(y - a.y, x - a.x) / DEG30) * DEG30;
    return polar(a.x, a.y, ang, BOND);
  }
  /* the emptiest direction out of an atom, for a plain click */
  function freeAngleIn(G, a) {
    const nb = neighboursIn(G, a.id);
    if (!nb.length) return 0;
    if (nb.length === 1) {                       // continue a zigzag: 120° from the one neighbour
      const back = Math.atan2(nb[0].y - a.y, nb[0].x - a.x);
      const c1 = back + 4 * DEG30, c2 = back - 4 * DEG30;
      return Math.sin(c1) <= Math.sin(c2) ? c1 : c2;                // prefer going "up"
    }
    let best = 0, bestScore = -1;
    for (let k = 0; k < 24; k++) {
      const ang = k * DEG30 / 2;
      const p = polar(a.x, a.y, ang, BOND);
      const score = Math.min(...nb.map(n => dist(p.x, p.y, n.x, n.y)));
      if (score > bestScore) { bestScore = score; best = ang; }
    }
    return best;
  }
  const freeAngle = a => freeAngleIn(g, a);

  /* ---- rings: a regular n-gon; existing atoms within reach are reused so rings fuse ---- */
  function ringRadius(n) { return BOND / (2 * Math.sin(Math.PI / n)); }
  function buildRing(n, cx, cy, startAng, dir, aromatic) {
    const ids = [];
    for (let k = 0; k < n; k++) {
      const p = polar(cx, cy, startAng + dir * k * 2 * Math.PI / n, ringRadius(n));
      const hit = atomAt(p.x, p.y);
      ids.push(hit ? hit.id : addAtom('C', p.x, p.y).id);
    }
    let parity = 0;
    if (aromatic) { const shared = bondBetween(ids[0], ids[1]); if (shared && shared.order === 1) parity = 1; }
    for (let k = 0; k < n; k++) {
      const p = ids[k], q = ids[(k + 1) % n];
      if (bondBetween(p, q)) continue;
      const hasDouble = id => bondsOf(id).some(b => b.order === 2);
      const order = aromatic && k % 2 === parity && !hasDouble(p) && !hasDouble(q) ? 2 : 1;
      g.bonds.push({ a: p, b: q, order });
    }
  }
  function ringAtPoint(n, x, y, aromatic) { buildRing(n, x, y, n % 2 ? -Math.PI / 2 : Math.PI / n + Math.PI / 2, 1, aromatic); }
  function ringOnAtom(n, a, aromatic) {
    const ang = freeAngle(a);
    const c = polar(a.x, a.y, ang, ringRadius(n));
    buildRing(n, c.x, c.y, ang + Math.PI, 1, aromatic);
  }
  function ringOnBond(n, b, aromatic) {
    const p = atom(b.a), q = atom(b.b);
    const mx = (p.x + q.x) / 2, my = (p.y + q.y) / 2;
    const L = dist(p.x, p.y, q.x, q.y) || 1;
    let nx = -(q.y - p.y) / L, ny = (q.x - p.x) / L;
    const side = [...neighbours(p.id), ...neighbours(q.id)].filter(o => o !== p && o !== q)
      .reduce((s, o) => s + Math.sign((o.x - mx) * nx + (o.y - my) * ny), 0);
    if (side > 0) { nx = -nx; ny = -ny; }
    const apothem = ringRadius(n) * Math.cos(Math.PI / n);
    const cx = mx + nx * apothem, cy = my + ny * apothem;
    const angP = Math.atan2(p.y - cy, p.x - cx), angQ = Math.atan2(q.y - cy, q.x - cx);
    let d = angQ - angP; while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI;
    buildRing(n, cx, cy, angP, Math.sign(d) || 1, aromatic);
  }

  /* ---- chain: zigzag from a start point towards (x, y) ---- */
  function chainPoints(sx, sy, x, y) {
    const L = dist(sx, sy, x, y);
    const n = Math.max(1, Math.round(L / (BOND * Math.cos(DEG30))));
    const ang = Math.round(Math.atan2(y - sy, x - sx) / (DEG30 / 2)) * (DEG30 / 2);
    const pts = [{ x: sx, y: sy }];
    for (let k = 0; k < n; k++) { const prev = pts[pts.length - 1]; pts.push(polar(prev.x, prev.y, ang + (k % 2 ? -DEG30 : DEG30), BOND)); }
    return pts;
  }

  /* ---- element tool drag: a short drag makes one bond, a longer one a zigzag chain of that element.
     Starts at the atom under the pointer, or at a fresh atom where empty space was pressed. ---- */
  function elDrag(d, x, y) {
    const s = d.from || { x: snapGrid(d.x), y: snapGrid(d.y) };
    const target = atomAt(x, y);
    if (d.from && target && target !== d.from) return { snapTo: target };
    if (dist(x, y, s.x, s.y) <= BOND * 1.4) return { pts: [s, snapDir(s, x, y)] };
    return { pts: chainPoints(s.x, s.y, x, y) };
  }

  /* ---- fragments and templates: a parsed, laid-out graph pasted into the drawing ---- */
  const fragCache = new Map();
  function fragmentGraph(smiles) {
    if (fragCache.has(smiles)) return fragCache.get(smiles);
    let G = null;
    try {
      G = C.parseSmiles(smiles);
      if (G.bonds.length) C.layoutGraph(G);
      const lens = G.bonds.map(b => { const p = G.atoms.find(a => a.id === b.a), q = G.atoms.find(a => a.id === b.b); return dist(p.x, p.y, q.x, q.y); });
      const scale = lens.length ? BOND / (lens.reduce((s, l) => s + l, 0) / lens.length) : 1;
      const o = G.atoms[0];
      const ox = o.x, oy = o.y;
      G.atoms.forEach(a => { a.x = (a.x - ox) * scale; a.y = (a.y - oy) * scale; });   // attachment atom at the origin
    } catch (e) { G = null; }
    fragCache.set(smiles, G);
    return G;
  }
  /* paste `frag`: attached to `anchor` through the fragment's first atom, or free with its first atom at (x, y);
     `centred` puts the fragment's centroid at the point instead (templates) */
  function pasteFragment(frag, x, y, anchor, centred) {
    const first = frag.atoms[0];
    let rot = 0, tx = x, ty = y;
    if (anchor) {
      const at = polar(anchor.x, anchor.y, freeAngle(anchor), BOND);
      tx = at.x; ty = at.y;
      /* turn the fragment so its own free direction points back at the anchor */
      const own = freeAngleIn(frag, first);
      rot = Math.atan2(anchor.y - at.y, anchor.x - at.x) - own;
    } else if (centred) {
      const c = centroid(frag.atoms); tx = x - c.x; ty = y - c.y;
    }
    const cos = Math.cos(rot), sin = Math.sin(rot), map = new Map();
    for (const a of frag.atoms) {
      const n = addAtom(a.element, tx + a.x * cos - a.y * sin, ty + a.x * sin + a.y * cos);
      n.charge = a.charge || 0; map.set(a.id, n.id);
    }
    for (const b of frag.bonds) g.bonds.push({ a: map.get(b.a), b: map.get(b.b), order: b.order || 1 });
    if (anchor) addBond(anchor.id, map.get(first.id), 1);
  }

  /* ---- clean-up: Chem's layout, one molecule at a time, then the molecules side by side ----
     (laid out together, the engine folds a vertex of the second ring; alone, every ring is exact) */
  /* one molecule (its atom ids) laid out by Chem, every bond scaled to our standard length; returns its atoms and box */
  function layoutPiece(ids) {
    const set = new Set(ids);
    const sub = { atoms: g.atoms.filter(a => set.has(a.id)), bonds: g.bonds.filter(b => set.has(b.a) && set.has(b.b)), nextId: g.nextId };
    window.MolDraw.layoutKeepingStereo(sub);      // keeps R/S and E/Z; falls back to the old picture
    const lens = sub.bonds.map(b => dist(atom(b.a).x, atom(b.a).y, atom(b.b).x, atom(b.b).y)).filter(l => l > 0);
    const scale = lens.length ? BOND / (lens.reduce((s, l) => s + l, 0) / lens.length) : 1;
    sub.atoms.forEach(a => { a.x *= scale; a.y *= scale; });
    const xs = sub.atoms.map(a => a.x), ys = sub.atoms.map(a => a.y);
    return { atoms: sub.atoms, x0: Math.min(...xs), x1: Math.max(...xs), y0: Math.min(...ys), y1: Math.max(...ys) };
  }
  const componentsNow = () => { try { return C.componentsOf(g); } catch (e) { return [g.atoms.map(a => a.id)]; } };
  const ringCount = () => g.bonds.length - g.atoms.length + componentsNow().length;   // independent rings (cyclomatic number)
  /* Closing a ring by hand leaves it lopsided (a pentagon drawn on the 30° grid is not regular), so when an edit adds a
     ring, the molecule holding the newest bond is redrawn in the standard layout, centred where it was drawn. Other
     molecules on the canvas stay put. One undo gives the hand-drawn version back. */
  function snapRingMolecule() {
    const last = g.bonds[g.bonds.length - 1]; if (!last) return;
    const ids = componentsNow().find(c => c.includes(last.a)); if (!ids || ids.length < 3) return;
    const before = centroid(g.atoms.filter(a => ids.includes(a.id)));
    const p = layoutPiece(ids);
    const dx = before.x - (p.x0 + p.x1) / 2, dy = before.y - (p.y0 + p.y1) / 2;
    p.atoms.forEach(a => { a.x += dx; a.y += dy; });
  }
  /* the Auto: Tidy toggle (on by default; remembered in this browser) */
  let autoTidy = true;
  try { autoTidy = localStorage.getItem('sk-autotidy') !== 'off'; } catch (e) {}
  function setAutoTidy(on) {
    autoTidy = on;
    try { localStorage.setItem('sk-autotidy', on ? 'on' : 'off'); } catch (e) {}
    const b = bar.querySelector('[data-act="autotidy"]');
    if (b) { b.classList.toggle('on', on); b.setAttribute('aria-pressed', String(on)); }
  }
  function cleanLayout() {
    if (!g.atoms.length) return;
    const pieces = componentsNow().map(layoutPiece);
    /* line the pieces up left to right with a gap, then centre the row in the box */
    const GAP = BOND * 1.5;
    const total = pieces.reduce((s, p) => s + (p.x1 - p.x0), 0) + GAP * (pieces.length - 1);
    const box = svg.getBoundingClientRect();
    let cursor = box.width / 2 - total / 2;
    for (const p of pieces) {
      const dx = cursor - p.x0, dy = box.height / 2 - (p.y0 + p.y1) / 2;
      p.atoms.forEach(a => { a.x += dx; a.y += dy; });
      cursor += (p.x1 - p.x0) + GAP;
    }
  }
  function mirror(axis) {
    const atoms = targetAtoms(); if (!atoms.length) return;
    const c = centroid(atoms);
    for (const a of atoms) { if (axis === 'x') a.y = 2 * c.y - a.y; else a.x = 2 * c.x - a.x; }
  }
  function addHydrogens() {
    for (const a of g.atoms.slice()) {
      if (a.element === 'H') continue;
      const n = C.implicitH(g, a);
      for (let k = 0; k < n; k++) { const p = polar(a.x, a.y, freeAngle(a), BOND * 0.7); const h = addAtom('H', p.x, p.y); addBond(a.id, h.id); }
    }
  }
  function stripHydrogens() {
    for (const a of g.atoms.filter(x => x.element === 'H' && !x.charge && bondsOf(x.id).length <= 1)) removeAtom(a.id);
  }

  /* ================= rendering ================= */
  const el = (name, attrs = {}, text) => {
    const n = document.createElementNS(NS, name);
    for (const k in attrs) n.setAttribute(k, attrs[k]);
    if (text != null) n.textContent = text;
    return n;
  };
  const labelled = a => a.element !== 'C' || a.charge || a.radical || bondsOf(a.id).length === 0;
  const chargeText = a => (Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '') + (a.charge > 0 ? '+' : '−');

  function drawBond(b, cls) {
    const phantom = (id, pt) => ({ id, element: 'C', charge: 0, x: pt.x, y: pt.y });
    const p = b._start ? phantom(-2, b._start) : atom(b.a);
    const q = b._end ? phantom(-1, b._end) : atom(b.b);
    const dx = q.x - p.x, dy = q.y - p.y, L = Math.hypot(dx, dy) || 1;
    const ux = dx / L, uy = dy / L, nx = -uy, ny = ux;
    const trimP = p.id > 0 && labelled(p) ? 11 : 0, trimQ = q.id > 0 && labelled(q) ? 11 : 0;
    const x1 = p.x + ux * trimP, y1 = p.y + uy * trimP, x2 = q.x - ux * trimQ, y2 = q.y - uy * trimQ;
    const gEl = el('g', { class: 'bond' + (cls ? ' ' + cls : ''), 'data-bond': g.bonds.indexOf(b) });

    if (b.stereo && b.order === 1) {
      const flip = b.narrow === q.id;
      const sx = flip ? x2 : x1, sy = flip ? y2 : y1, ex = flip ? x1 : x2, ey = flip ? y1 : y2;
      if (b.stereo === 'wedge') {
        gEl.appendChild(el('polygon', { class: 'wedge', points: `${sx},${sy} ${ex + nx * 4},${ey + ny * 4} ${ex - nx * 4},${ey - ny * 4}` }));
      } else {
        const len = dist(sx, sy, ex, ey), steps = Math.max(3, Math.floor(len / 4));
        for (let k = 1; k <= steps; k++) {
          const t = k / steps, w = 0.8 + 3.4 * t, px = sx + (ex - sx) * t, py = sy + (ey - sy) * t;
          gEl.appendChild(el('line', { class: 'hash', x1: px + nx * w, y1: py + ny * w, x2: px - nx * w, y2: py - ny * w }));
        }
      }
      return gEl;
    }
    let offsets = { 1: [0], 2: [0, 6], 3: [-6, 0, 6] }[b.order] || [0];
    if (b.order === 2) {
      const others = [...(p.id > 0 ? neighbours(p.id) : []), ...(q.id > 0 ? neighbours(q.id) : [])].filter(o => o !== p && o !== q);
      if (!others.length) offsets = [-3, 3];
      else { const side = others.reduce((s, o) => s + Math.sign((o.x - p.x) * nx + (o.y - p.y) * ny), 0); offsets = [0, side < 0 ? -6 : 6]; }
    }
    for (const o of offsets) {
      const shrink = b.order === 2 && o !== 0 ? 5 : 0;
      gEl.appendChild(el('line', { x1: x1 + nx * o + ux * shrink, y1: y1 + ny * o + uy * shrink, x2: x2 + nx * o - ux * shrink, y2: y2 + ny * o - uy * shrink }));
    }
    return gEl;
  }

  /* is (x, y) on the H of this atom's label ("OH": the H right of the O)? The label is centred on the atom, 10 px a
     character, 3 px padding (drawAtom) */
  function onHydrogenOf(a, x) {
    if (!labelled(a) || a.element === 'H') return false;
    let h = 0; try { h = C.implicitH(g, a); } catch (e) { h = 0; }
    if (!h) return false;
    const text = a.element + 'H' + (h > 1 ? h : '') + (a.charge ? chargeText(a) : '') + (a.radical ? '•' : '');
    const left = a.x - (10 * text.length + 6) / 2 + 3;
    return x >= left + 10 * a.element.length;
  }
  function drawAtom(a) {
    const gEl = el('g', { class: 'atom' + (hover.atom === a ? ' hover' : ''), 'data-atom': a.id });
    if (selected.has(a.id)) gEl.appendChild(el('circle', { class: 'sel', cx: a.x, cy: a.y, r: 13 }));
    if (hover.atom === a && !labelled(a)) gEl.appendChild(el('circle', { class: 'atom-hover', cx: a.x, cy: a.y, r: 7 }));
    if (!labelled(a)) return gEl;
    const h = a.element === 'H' ? 0 : C.implicitH(g, a);
    const text = a.element + (h ? 'H' : '') + (h > 1 ? h : '') + (a.charge ? chargeText(a) : '') + (a.radical ? '•' : '');
    const w = 10 * text.length + 6;
    gEl.appendChild(el('rect', { class: 'label-bg', x: a.x - w / 2, y: a.y - 11, width: w, height: 22, rx: 3 }));
    const t = el('text', { class: 'label', x: a.x, y: a.y }, a.element);
    if (h) t.appendChild(el('tspan', { class: 'hcount' }, 'H' + (h > 1 ? h : '')));
    if (a.charge) t.appendChild(el('tspan', { class: 'charge' }, chargeText(a)));
    if (a.radical) t.appendChild(el('tspan', { class: 'charge' }, '•'));
    gEl.appendChild(t);
    return gEl;
  }

  function render(preview) {
    svg.replaceChildren();
    for (const b of g.bonds) svg.appendChild(drawBond(b, hover.bond === b ? 'hover' : ''));
    if (preview && preview.bond) svg.appendChild(drawBond(preview.bond, 'preview'));
    if (preview && preview.points) for (let k = 1; k < preview.points.length; k++) svg.appendChild(drawBond({ order: 1, _start: preview.points[k - 1], _end: preview.points[k] }, 'preview'));
    if (preview && preview.ghost) for (const p of preview.ghost) svg.appendChild(el('circle', { class: 'ghost', cx: p.x, cy: p.y, r: 2.5 }));
    for (const a of g.atoms) svg.appendChild(drawAtom(a));
    if (preview && preview.marquee) {
      const m = preview.marquee;
      svg.appendChild(el('rect', { class: 'marquee', x: Math.min(m.x1, m.x2), y: Math.min(m.y1, m.y2), width: Math.abs(m.x2 - m.x1), height: Math.abs(m.y2 - m.y1) }));
    }
    if (!preview) readout();
  }

  /* ================= readout: formula, SMILES, IUPAC ================= */
  const esc = s => String(s).replace(/[&<>]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;' }[c]));
  const escAttr = s => esc(s).replace(/"/g, '&quot;');
  /* A common name adds nothing when it is the systematic name again: the same apart from spacing, hyphens and case
     (benzene), or the same with its locants left out ("butanol" for butan-1-ol — which, without the locant, could
     be any of the butanols, so it is not shown as this molecule's name). Function declarations, so readout() can use
     them wherever it runs. */
  function plainName(s) { return String(s).toLowerCase().replace(/[\s\-–,'()]/g, ''); }
  function sameName(common, iupac) {
    const c = plainName(common), n = plainName(iupac);
    return c === n || c === plainName(String(iupac).replace(/\b\d+(,\d+)*-/g, '').replace(/-\d+(,\d+)*-/g, '-'));
  }
  function readout() {
    if (!g.atoms.length) {
      info.innerHTML = '<span class="hint">Pick a tool on the left and click in the box. Drag out of an atom to grow a bond; ' +
        'hover an atom and type an element letter to change it.</span>';
      return;
    }
    let err = '';
    try { C.validateValences(g); } catch (e) { err = e.message; }
    let stacked = 0;
    for (let i = 0; i < g.atoms.length; i++) for (let j = i + 1; j < g.atoms.length; j++)
      if (dist(g.atoms[i].x, g.atoms[i].y, g.atoms[j].x, g.atoms[j].y) < 6) stacked++;
    if (stacked && !notice) notice = `${stacked === 1 ? 'Two atoms lie' : stacked + ' pairs of atoms lie'} exactly on top of each other, so the drawing holds atoms you cannot see (the formula and SMILES include them). Press Tidy to pull them apart, or drag the atom away.`;
    let smiles = '', name = '', formula = '', common = '';
    try { formula = C.formula(g); } catch (e) {}
    try { smiles = C.toSmiles(g); } catch (e) {}
    if (!err) {
      try { name = C.iupacName(g); } catch (e) {}
      /* the everyday name, when the library knows one and it is not just the systematic name again
         (ethyl ethanoate is named "ethyl acetate"; butan-1-ol has no other name, so no second row) */
      try { common = C.commonNameFor(g) || ''; } catch (e) {}
      if (common && name && sameName(common, name)) common = '';
    }
    const row = (label, value, cls) => `<div class="ro-row" data-copy="${escAttr(value)}" title="Click to copy the ${label} to the clipboard" tabindex="0" role="button">` +
      `<b>${label}</b> <span class="${cls || ''}">${esc(value)}</span><span class="ro-copied" aria-live="polite"></span></div>`;
    info.innerHTML = (notice ? `<div class="err">${esc(notice)}</div>` : '') +
      row('Formula', formula) + row('SMILES', smiles, 'ro-smiles') +
      (err ? `<div><b>IUPAC</b> <span class="err">${esc(err)}</span></div>` : row('IUPAC', name)) +
      (common ? row('Common name', common) : '');
    notice = '';
  }
  /* click (or Enter / Space) a readout row to copy its value */
  async function copyRow(el) {
    const text = el.getAttribute('data-copy') || '';
    if (!text) return;
    const flash = el.querySelector('.ro-copied');
    let ok = false;
    try { await navigator.clipboard.writeText(text); ok = true; } catch (e) { ok = false; }
    if (!ok) {
      /* a page opened from a file (or a browser that refuses the clipboard API here) — select the value and use the
         old copy command; if that is refused too, the text is left selected so Ctrl+C works */
      const span = el.querySelector('span');
      try {
        const r = document.createRange(); r.selectNodeContents(span);
        const sel = getSelection(); sel.removeAllRanges(); sel.addRange(r);
        ok = document.execCommand('copy');
        if (ok) sel.removeAllRanges();
      } catch (e2) {}
    }
    if (flash) { flash.textContent = ok ? 'copied' : 'selected — press Ctrl+C'; clearTimeout(flash._t); flash._t = setTimeout(() => { flash.textContent = ''; }, 1600); }
  }
  info.addEventListener('click', e => { const r = e.target.closest('.ro-row'); if (r) copyRow(r); });
  info.addEventListener('keydown', e => {
    const r = e.target.closest('.ro-row');
    if (r && (e.key === 'Enter' || e.key === ' ')) { e.preventDefault(); copyRow(r); }
  });
  const currentSmiles = () => { try { return C.toSmiles(g); } catch (e) { return ''; } };
  const currentName = () => { try { C.validateValences(g); return C.iupacName(g); } catch (e) { return ''; } };

  /* ================= palette ================= */
  const ico = inner => `<svg viewBox="0 0 24 24">${inner}</svg>`;
  const poly = (n, aromatic) => {
    const r = 9.5, pts = [];
    for (let k = 0; k < n; k++) { const a = -Math.PI / 2 + k * 2 * Math.PI / n; pts.push(`${12 + r * Math.cos(a)},${12.5 + r * Math.sin(a)}`); }
    return ico(`<polygon points="${pts.join(' ')}"/>` + (aromatic ? '<circle cx="12" cy="12.5" r="5.5"/>' : ''));
  };
  const TOOLS = [
    { id: 'select', title: 'Select / move: drag atoms; drag empty space for a marquee; Delete removes', html: ico('<path d="M6 3l12 10h-7l-3 7z"/>') },
    { id: 'rotate', title: 'Rotate: drag to turn the selection (or the whole drawing)', html: ico('<path d="M19 12a7 7 0 1 1-2.5-5.4"/><path d="M19 4v4h-4"/>') },
    { id: 'ptable', title: 'Periodic table: pick any element', html: ico('<rect x="3" y="5" width="18" height="14" rx="1"/><path d="M3 10h18M3 15h18M8 5v14M13 5v14M18 5v14"/>') },
    { id: 'el:C', title: 'Carbon: click to place; drag to draw a chain', html: 'C' }, { id: 'el:N', title: 'Nitrogen: click to place; drag to draw a chain', html: 'N' },
    { id: 'el:O', title: 'Oxygen: click to place; drag to draw a chain', html: 'O' }, { id: 'el:S', title: 'Sulfur: click to place; drag to draw a chain', html: 'S' },
    { id: 'el:P', title: 'Phosphorus: click to place; drag to draw a chain', html: 'P' }, { id: 'el:F', title: 'Fluorine', html: 'F' },
    { id: 'el:Cl', title: 'Chlorine', html: 'Cl' }, { id: 'el:Br', title: 'Bromine', html: 'Br' },
    { id: 'el:I', title: 'Iodine', html: 'I' }, { id: 'el:H', title: 'Explicit hydrogen', html: 'H' },
    { id: 'chain', title: 'Chain: drag to draw a zigzag of carbons', html: ico('<polyline points="3,16 8,8 13,16 18,8 22,16"/>') },
    { id: 'bond:1', title: 'Single bond', html: ico('<line x1="4" y1="12" x2="20" y2="12"/>') },
    { id: 'bond:2', title: 'Double bond', html: ico('<line x1="4" y1="9.5" x2="20" y2="9.5"/><line x1="4" y1="14.5" x2="20" y2="14.5"/>') },
    { id: 'bond:3', title: 'Triple bond', html: ico('<line x1="4" y1="7.5" x2="20" y2="7.5"/><line x1="4" y1="12" x2="20" y2="12"/><line x1="4" y1="16.5" x2="20" y2="16.5"/>') },
    { id: 'wedge', title: 'Up-wedge bond (towards the viewer); click a wedge again to flip it', html: ico('<polygon points="4,12 20,7.5 20,16.5" style="fill:#111"/>') },
    { id: 'hash', title: 'Down-wedge (hashed) bond; click a hash again to flip it', html: ico('<line x1="6" y1="11" x2="6" y2="13"/><line x1="9.5" y1="10" x2="9.5" y2="14"/><line x1="13" y1="9" x2="13" y2="15"/><line x1="16.5" y1="8" x2="16.5" y2="16"/><line x1="20" y1="7" x2="20" y2="17"/>') },
    { id: 'charge:1', title: 'Positive charge: click an atom', html: ico('<circle cx="12" cy="12" r="8"/><path d="M12 8v8M8 12h8"/>') },
    { id: 'charge:-1', title: 'Negative charge: click an atom', html: ico('<circle cx="12" cy="12" r="8"/><path d="M8 12h8"/>') },
    { id: 'radical', title: 'Radical (unpaired electron): click an atom', html: ico('<circle cx="12" cy="12" r="8"/><circle cx="12" cy="12" r="1.6" style="fill:#111"/>') },
    { id: 'ring:6:1', title: 'Benzene', html: poly(6, true) },
    { id: 'ring:3', title: 'Cyclopropane', html: poly(3) }, { id: 'ring:4', title: 'Cyclobutane', html: poly(4) },
    { id: 'ring:5', title: 'Cyclopentane', html: poly(5) }, { id: 'ring:6', title: 'Cyclohexane', html: poly(6) },
    { id: 'ring:7', title: 'Cycloheptane', html: poly(7) }, { id: 'ring:8', title: 'Cyclooctane', html: poly(8) },
    { id: 'tpl', title: 'Templates: common ring systems and molecules to paste', html: ico('<rect x="3" y="3" width="7" height="7"/><rect x="14" y="3" width="7" height="7"/><rect x="3" y="14" width="7" height="7"/><rect x="14" y="14" width="7" height="7"/>') },
    { id: 'erase', title: 'Eraser: an atom removes it with its H (the whole OH on the O); the H of a label (the H in OH) removes that H as H⁺, leaving O⁻; a bond steps down one order per click (≡ → = → – → gone)', html: ico('<path d="M4 15l8-8 7 7-4 4H9z"/><path d="M9 18l-2-2"/>') },
  ];
  grid.innerHTML = TOOLS.map(t => `<button class="tool" data-id="${t.id}" title="${t.title}">${t.html}</button>`).join('');

  /* two rows of fragments: alkyl groups (as PubChem has them), then functional groups by name.
     The SMILES starts at the attachment atom. */
  const ALKYLS = [
    ['Et', 'CC', 'Ethyl'], ['nPr', 'CCC', 'n-Propyl'], ['iPr', 'C(C)C', 'Isopropyl'], ['nBu', 'CCCC', 'n-Butyl'],
    ['iBu', 'CC(C)C', 'Isobutyl'], ['sBu', 'C(C)CC', 'sec-Butyl'], ['tBu', 'C(C)(C)C', 'tert-Butyl'], ['Ph', 'c1ccccc1', 'Phenyl'],
    ['Bn', 'Cc1ccccc1', 'Benzyl'], ['Vinyl', 'C=C', 'Vinyl'], ['Allyl', 'CC=C', 'Allyl'],
  ];
  const GROUPS = [
    ['Alcohol', 'OH', 'O'], ['Ether', 'OMe', 'OC'], ['Thiol', 'SH', 'S'], ['Amine', 'NH<sub>2</sub>', 'N'],
    ['Aldehyde', 'CHO', 'C=O'], ['Ketone', 'C(=O)CH<sub>3</sub>', 'C(=O)C'], ['Carboxylic acid', 'CO<sub>2</sub>H', 'C(=O)O'],
    ['Ester', 'CO<sub>2</sub>Me', 'C(=O)OC'], ['Amide', 'CONH<sub>2</sub>', 'C(=O)N'], ['Acyl chloride', 'COCl', 'C(=O)Cl'],
    ['Nitrile', 'C≡N', 'C#N'], ['Nitro', 'NO<sub>2</sub>', '[N+](=O)[O-]'], ['Sulfonic acid', 'SO<sub>3</sub>H', 'S(=O)(=O)O'],
  ];
  fragBar.innerHTML =
    `<div class="row">${ALKYLS.map(([label, smi, name]) => `<button class="tool" data-id="frag:${smi}" title="${name}: click an atom to attach it, or empty space to place it">${label}</button>`).join('')}</div>` +
    `<div class="row">${GROUPS.map(([name, formula, smi]) => `<button class="tool" data-id="frag:${smi}" title="${name}: click an atom to attach it, or empty space to place it">${name}<small>${formula}</small></button>`).join('')}</div>`;

  const TEMPLATES = [
    ['Naphthalene', 'c1ccc2ccccc2c1'], ['Anthracene', 'c1ccc2cc3ccccc3cc2c1'], ['Phenanthrene', 'c1ccc2c(c1)ccc1ccccc12'],
    ['Biphenyl', 'c1ccc(cc1)-c1ccccc1'], ['Pyridine', 'c1ccncc1'], ['Pyrimidine', 'c1cncnc1'], ['Pyrrole', 'c1cc[nH]c1'],
    ['Furan', 'c1ccoc1'], ['Thiophene', 'c1ccsc1'], ['Imidazole', 'c1cnc[nH]1'], ['Indole', 'c1ccc2[nH]ccc2c1'],
    ['Quinoline', 'c1ccc2ncccc2c1'], ['Purine', 'c1ncc2nc[nH]c2n1'], ['Cyclohexene', 'C1CCC=CC1'], ['Cyclopentadiene', 'C1C=CC=C1'],
    ['Piperidine', 'C1CCNCC1'], ['Pyrrolidine', 'C1CCNC1'], ['Morpholine', 'C1COCCN1'], ['Tetrahydrofuran', 'C1CCOC1'],
    ['1,4-Dioxane', 'C1COCCO1'], ['Norbornane', 'C1CC2CCC1C2'], ['Adamantane', 'C1C2CC3CC1CC(C2)C3'], ['Decalin', 'C1CCC2CCCCC2C1'],
    ['Steroid nucleus', 'C1CCC2C(C1)CCC1C2CCC2CCCC12'], ['Toluene', 'Cc1ccccc1'], ['Phenol', 'Oc1ccccc1'], ['Aniline', 'Nc1ccccc1'],
    ['Benzoic acid', 'OC(=O)c1ccccc1'], ['Acetic acid', 'CC(=O)O'], ['Acetone', 'CC(C)=O'], ['Ethyl acetate', 'CCOC(C)=O'],
    ['Urea', 'NC(N)=O'], ['Glycine', 'NCC(=O)O'], ['Alanine', 'CC(N)C(=O)O'],
  ];

  function pickTool(id) {
    if (id === 'ptable') { openPeriodicTable(); return; }
    if (id === 'tpl') { openTemplates(); return; }
    closePopup();
    const [kind, ...rest] = id.split(':');
    tool = { id, kind, arg: rest.join(':'), arg2: rest[1] };
    if (kind === 'ring') { tool.arg = rest[0]; }
    document.querySelectorAll('.tool').forEach(b => b.classList.toggle('active', b.dataset.id === id || (b.dataset.id === 'ptable' && kind === 'el' && b.dataset.el === tool.arg)));
    svg.style.cursor = kind === 'select' ? 'default' : kind === 'rotate' ? 'grab' : kind === 'erase' ? 'not-allowed' : 'crosshair';
    if (kind !== 'select' && kind !== 'rotate') { selected.clear(); render(); }
  }
  document.addEventListener('click', e => {
    const b = e.target.closest('.tool'); if (!b) return;
    /* the periodic-table button remembers the last element picked: one click re-selects it,
       a second click (when it is already the tool) opens the table again */
    const remembered = b.dataset.id === 'ptable' && b.dataset.el;
    const alreadyOn = remembered && tool && tool.kind === 'el' && tool.arg === b.dataset.el;
    pickTool(remembered && !alreadyOn ? 'el:' + b.dataset.el : b.dataset.id);
  });
  pickTool('el:C');

  /* ---- popups ---- */
  function closePopup() { if (popup) { popup.remove(); popup = null; } }
  function openPopup(title, bodyHtml) {
    closePopup();
    popup = document.createElement('div');
    popup.className = 'popup';
    popup.innerHTML = `<h4>${title}<button type="button" title="Close">✕</button></h4>${bodyHtml}`;
    popup.querySelector('h4 button').addEventListener('click', closePopup);
    main.appendChild(popup);
  }
  document.addEventListener('pointerdown', e => { if (popup && !popup.contains(e.target) && !e.target.closest('.tool')) closePopup(); });

  /* the periodic table: [symbol, period, group]; lanthanides/actinides on two extra rows */
  const PT = 'H:1:1 He:1:18 Li:2:1 Be:2:2 B:2:13 C:2:14 N:2:15 O:2:16 F:2:17 Ne:2:18 Na:3:1 Mg:3:2 Al:3:13 Si:3:14 P:3:15 S:3:16 Cl:3:17 Ar:3:18 ' +
    'K:4:1 Ca:4:2 Sc:4:3 Ti:4:4 V:4:5 Cr:4:6 Mn:4:7 Fe:4:8 Co:4:9 Ni:4:10 Cu:4:11 Zn:4:12 Ga:4:13 Ge:4:14 As:4:15 Se:4:16 Br:4:17 Kr:4:18 ' +
    'Rb:5:1 Sr:5:2 Y:5:3 Zr:5:4 Nb:5:5 Mo:5:6 Tc:5:7 Ru:5:8 Rh:5:9 Pd:5:10 Ag:5:11 Cd:5:12 In:5:13 Sn:5:14 Sb:5:15 Te:5:16 I:5:17 Xe:5:18 ' +
    'Cs:6:1 Ba:6:2 Lu:6:3 Hf:6:4 Ta:6:5 W:6:6 Re:6:7 Os:6:8 Ir:6:9 Pt:6:10 Au:6:11 Hg:6:12 Tl:6:13 Pb:6:14 Bi:6:15 Po:6:16 At:6:17 Rn:6:18 ' +
    'Fr:7:1 Ra:7:2 Lr:7:3 ' +
    'La:9:3 Ce:9:4 Pr:9:5 Nd:9:6 Pm:9:7 Sm:9:8 Eu:9:9 Gd:9:10 Tb:9:11 Dy:9:12 Ho:9:13 Er:9:14 Tm:9:15 Yb:9:16 ' +
    'Ac:10:3 Th:10:4 Pa:10:5 U:10:6 Np:10:7 Pu:10:8 Am:10:9 Cm:10:10 Bk:10:11 Cf:10:12 Es:10:13 Fm:10:14 Md:10:15 No:10:16';
  function openPeriodicTable() {
    const cells = PT.split(' ').map(s => s.split(':')).map(([sym, p, gcol]) => {
      const ok = C.VALENCE[sym] !== undefined;
      return `<button data-sym="${sym}" style="grid-row:${p};grid-column:${gcol}" ${ok ? '' : 'disabled'} title="${ok ? sym : sym + ' — not in the namer\'s valence table'}">${sym}</button>`;
    }).join('');
    openPopup('Periodic table', `<div class="ptable">${cells}</div>`);
    popup.querySelector('.ptable').addEventListener('click', e => {
      const b = e.target.closest('button[data-sym]'); if (!b || b.disabled) return;
      const pt = grid.querySelector('[data-id="ptable"]');
      pt.dataset.el = b.dataset.sym; pt.innerHTML = b.dataset.sym;
      pickTool('el:' + b.dataset.sym);
    });
  }
  function openTemplates() {
    const items = TEMPLATES.filter(([, smi]) => fragmentGraph(smi)).map(([name, smi]) =>
      `<button data-smi="${smi}">${name}<small>${esc(smi)}</small></button>`).join('');
    openPopup('Templates — pick one, then click in the box to place it', `<div class="tpl-grid">${items}</div>`);
    popup.querySelector('.tpl-grid').addEventListener('click', e => {
      const b = e.target.closest('button[data-smi]'); if (!b) return;
      pickTool('paste:' + b.dataset.smi);
      grid.querySelector('[data-id="tpl"]').classList.add('active');
    });
  }

  /* ================= toolbar ================= */
  function loadGraph(parsed) {
    snapshot();
    g = { atoms: parsed.atoms.map(a => ({ id: a.id, element: a.element, charge: a.charge || 0, ...(a.radical ? { radical: true } : {}), x: a.x || 0, y: a.y || 0 })),
          bonds: parsed.bonds.map(b => ({ a: b.a, b: b.b, order: b.order || 1, ...(b.stereo ? { stereo: b.stereo, narrow: b.narrow } : {}) })),
          nextId: parsed.nextId || Math.max(0, ...parsed.atoms.map(a => a.id)) + 1 };
    selected.clear();
    cleanLayout();
    render();
    rdkitMacrocycleLayout();
  }
  /* Chem's layout draws rings of up to 16 atoms and leaves a bigger ring's closing bond as one long line across the
     page (a 38-membered macrolide came out with a bond 27 times the standard length). RDKit's layout is what the
     Molecule tab draws, so a molecule with such a ring is laid out by RDKit once it has started (a few seconds the
     first time) and the drawing is redrawn. It is used only when RDKit reads the new drawing as exactly the same
     molecule (configuration of every stereocentre and double bond included); otherwise the drawing stays as it was. */
  function hasMacrocycle(gr) {
    const adj = new Map(gr.atoms.map(a => [a.id, []]));
    gr.bonds.forEach(b => { adj.get(b.a).push(b.b); adj.get(b.b).push(b.a); });
    for (const b of gr.bonds) {                       // the smallest ring through each bond: the shortest path between its ends that avoids it
      const dist = new Map([[b.a, 0]]), queue = [b.a];
      while (queue.length) {
        const cur = queue.shift();
        if (cur === b.b) break;
        for (const nb of adj.get(cur)) {
          if (cur === b.a && nb === b.b) continue;
          if (!dist.has(nb)) { dist.set(nb, dist.get(cur) + 1); queue.push(nb); }
        }
      }
      if (dist.has(b.b) && dist.get(b.b) + 1 > 16) return true;
    }
    return false;
  }
  async function rdkitMacrocycleLayout() {
    if (!window.RDKitLoad || !g.atoms.length || !hasMacrocycle(g) || componentsNow().length !== 1) return;
    const mine = g, n = g.atoms.length, m = g.bonds.length;
    try {
      const R = await window.RDKitLoad();
      if (g !== mine || g.atoms.length !== n || g.bonds.length !== m) return;          // the drawing was replaced or edited meanwhile
      const mb = C.toMolfile(g, '');
      const ref = R.get_mol(mb);
      if (!ref) return;
      const want = ref.get_smiles(); ref.delete();
      const M = R.get_mol(mb);
      M.set_new_coords(true); M.normalize_depiction(0, -1);
      const mb2 = M.get_molblock(); M.delete();
      const back = R.get_mol(mb2), same = !!back && back.get_smiles() === want;
      if (back) back.delete();
      if (!same) return;
      const p = C.parseMolfile(mb2);
      if (p.atoms.length !== n) return;
      const idOf = new Map(p.atoms.map((a, i) => [a.id, g.atoms[i].id]));              // the molfile keeps the atom order
      const lens = p.bonds.map(b => { const u = p.atoms.find(a => a.id === b.a), v = p.atoms.find(a => a.id === b.b); return Math.hypot(u.x - v.x, u.y - v.y); }).filter(l => l > 0).sort((u, v) => u - v);
      const xs = p.atoms.map(a => a.x), ys = p.atoms.map(a => a.y), w = Math.max(...xs) - Math.min(...xs), h = Math.max(...ys) - Math.min(...ys);
      const box = svg.getBoundingClientRect();
      let k = BOND / (lens[lens.length >> 1] || 1);
      if (w * k > box.width * 0.92) k = box.width * 0.92 / w;                          // a big ring is drawn smaller to fit the canvas (not below half a bond)
      if (h * k > box.height * 0.92) k = Math.min(k, box.height * 0.92 / h);
      k = Math.max(k, 0.5 * BOND / (lens[lens.length >> 1] || 1));
      const cx = box.width / 2 - k * (Math.min(...xs) + Math.max(...xs)) / 2, cy = box.height / 2 - k * (Math.min(...ys) + Math.max(...ys)) / 2;
      p.atoms.forEach((a, i) => { g.atoms[i].x = a.x * k + cx; g.atoms[i].y = a.y * k + cy; });
      g.bonds.forEach(b => {                                                           // wedges and dashes as RDKit drew them
        delete b.stereo; delete b.narrow;
        const q = p.bonds.find(x => (idOf.get(x.a) === b.a && idOf.get(x.b) === b.b) || (idOf.get(x.a) === b.b && idOf.get(x.b) === b.a));
        if (q && q.stereo) { b.stereo = q.stereo; if (q.narrow != null) b.narrow = idOf.get(q.narrow); }
      });
      render();
    } catch (e) { /* RDKit unavailable or refused: the drawing stays as Chem laid it out */ }
  }
  function importText(text) {
    text = text.trim(); if (!text) return;
    closeSuggest();
    let parsed = null, error = '';
    if (C.looksLikeMolfile && C.looksLikeMolfile(text)) { try { parsed = C.parseMolfile(text); } catch (e) { error = e.message; } }
    /* searchMolecule knows common names (aspirin), the molecule library, CAS numbers, IUPAC names and SMILES */
    if (!parsed) { try { const r = C.searchMolecule(text); parsed = r && r.graph; } catch (e) { error = e.message; } }
    if (!parsed) { try { parsed = C.parseSmiles(text); } catch (e) { error = error || e.message; } }
    if (!parsed || !parsed.atoms.length) { try { parsed = C.parseIupacName(text); } catch (e) { error = error || e.message; } }
    if (!parsed || !parsed.atoms || !parsed.atoms.length) {
      if ((/[a-z]{3}/i.test(text) || /^\d{2,7}-\d\d-\d$/.test(text)) && !/[=#\[\]@]/.test(text)) { importFallback(text, error); return; }
      alert('Could not read "' + text.slice(0, 60) + '".\n' + error); return;
    }
    loadGraph(parsed);
  }
  /* a name or CAS number the curated library does not know: the saved PubChem lookups and the broad offline library
     first (js/pclib.js; no internet), then, asked first, PubChem itself */
  async function importFallback(text, error) {
    const L = window.PCLib;
    if (L) {
      info.innerHTML = '<span class="hint">Looking in the offline PubChem library…</span>';
      let res = await L.resolve(text);
      if (!res.ok && res.choices && res.choices.length) res = choose(res, text);
      if (res && res.ok) { importSmiles(res.smiles, res.note); return; }
      if (res === null) { render(); return; }                 // the choice was cancelled
    }
    const P = window.PubChem3D;
    if (P && P.lookupName && confirm('"' + text.slice(0, 60) + '" is not in the app\'s offline libraries.\n\nLook it up on PubChem? (This sends the text to PubChem, NCBI, over the internet.)')) {
      pubchemImport(text);
      return;
    }
    render();
    alert('Could not read "' + text.slice(0, 60) + '".\n' + error);
  }
  // several different compounds under one name or CAS number: the user picks one (null = cancelled)
  function choose(res, text) {
    const list = res.choices.map((c, i) => `${i + 1}. ${c.title}${c.formula ? ' (' + c.formula + ')' : ''}`).join('\n');
    const n = parseInt(prompt(res.error + '\n\n' + list + '\n\nType the number of the one you mean:', '1'), 10);
    const c = res.choices[n - 1];
    return c ? { ok: true, smiles: c.smiles, note: `From ${res.src || 'PubChem'}: CID ${c.cid}, ${c.title} (chosen from ${res.choices.length} records under "${text}"${res.casNote || ''}).` } : null;
  }
  function importSmiles(smiles, note) {
    let g2 = null;
    try { g2 = C.parseSmiles(smiles); } catch (e) { render(); alert('The structure found could not be read: ' + e.message); return; }
    loadGraph(g2);
    notice = note;
    render();
  }
  async function pubchemImport(text) {
    info.innerHTML = '<span class="hint">Asking PubChem…</span>';
    let res = await window.PubChem3D.lookupName(text);
    if (!res.ok && res.choices && res.choices.length) res = choose(Object.assign({ src: 'PubChem' }, res), text);
    if (!res) { render(); return; }
    if (!res.ok) { render(); alert(res.error); return; }
    importSmiles(res.smiles, res.note + (res.saved ? ' (from your saved lookups)' : ''));
  }

  /* ---- name suggestions under a text box: the molecule library's names, alternative names
     and common names, prefix matches first. attachSuggest() is shared with the other pages. ---- */
  let NAMES = null;
  function nameIndex() {
    if (NAMES) return NAMES;
    const seen = new Set(); NAMES = [];
    const add = (label, formula) => { const k = label.toLowerCase(); if (!seen.has(k)) { seen.add(k); NAMES.push({ label, formula: formula || '', k }); } };
    (C.MOL_LIBRARY || []).forEach(e => { add(e.name, e.formula); (e.alt || []).forEach(a => add(a, e.formula)); });
    Object.keys(C.COMMON_NAMES || {}).forEach(n => add(n, ''));
    return NAMES;
  }
  /* The Sketcher's search box: every name the molecule library knows (its names, mostly IUPAC, and its alternative
     and common names) plus the common-name list, each with its formula and CAS number where the library has one; and
     the CAS numbers themselves, so typing digits offers "50-78-2  aspirin". */
  let SEARCH = null;
  function searchIndex() {
    if (SEARCH) return SEARCH;
    const seen = new Set(), names = [], cas = [];
    const add = (label, formula, casNo) => {
      const k = label.toLowerCase(); if (seen.has(k)) return; seen.add(k);
      names.push({ label, formula: formula || '', cas: casNo || '', k });
    };
    (C.MOL_LIBRARY || []).forEach(e => {
      add(e.name, e.formula, e.cas); (e.alt || []).forEach(a => add(a, e.formula, e.cas));
      if (e.cas && !seen.has('cas:' + e.cas)) { seen.add('cas:' + e.cas); cas.push({ label: e.cas, name: e.name, formula: e.formula || '', k: e.cas }); }
    });
    Object.keys(C.COMMON_NAMES || {}).forEach(n => add(n, '', ''));
    SEARCH = { names, cas };
    return SEARCH;
  }
  /* opts.autofill === false: no inline completion and nothing pre-selected, so Enter submits exactly what was typed
     (the Molecule tab, where "CO" is a SMILES, not the start of "CO2"); the default keeps the behaviour below */
  function attachSuggest(input, box, onSubmit, indexFn, opts) {
    let idx = -1;
    const index = indexFn || nameIndex;
    const autofill = !opts || opts.autofill !== false;
    const searchMode = !!(opts && opts.search);        // the Sketcher: names incl. IUPAC names, and CAS numbers
    const close = () => { box.hidden = true; box.innerHTML = ''; idx = -1; };
    const move = d => {
      const items = [...box.children]; if (!items.length) return;
      idx = (idx + d + items.length) % items.length;
      items.forEach((it, i) => it.classList.toggle('on', i === idx));
      items[idx].scrollIntoView({ block: 'nearest' });
    };
    const suggest = e => {
      /* what the user actually typed: the part before any autofilled, still-selected tail */
      const typed = input.selectionStart < input.value.length && input.selectionEnd === input.value.length
        ? input.value.slice(0, input.selectionStart) : input.value;
      const q = typed.trim().toLowerCase();
      if (searchMode) { searchSuggest(typed.trim(), q, e); return; }
      if (q.length < 2 || /^[\[\]()=#@+\-\d]/.test(q)) { close(); return; }     // SMILES / CAS: no name list
      const starts = [], contains = [];
      for (const n of index()) {
        if (n.k.startsWith(q)) starts.push(n); else if (n.k.includes(q)) contains.push(n);
        if (starts.length >= 12) break;
      }
      starts.sort((a, b) => a.label.length - b.label.length);   // "prop" offers propane before propan-2-ol
      const hits = starts.concat(contains).slice(0, 12);
      if (!hits.length) { close(); return; }
      box.innerHTML = hits.map(h => `<div data-name="${esc(h.label)}">${esc(h.label)}<span>${esc(h.formula)}</span></div>`).join('');
      box.hidden = false; idx = -1;
      /* inline autofill while typing forward (not when deleting): complete to the best name and
         keep the added tail selected, so typing on replaces it and Enter accepts it */
      const typing = e && e.type === 'input' && e.inputType && !e.inputType.startsWith('delete');
      if (autofill && typing && starts.length && starts[0].k !== q) {
        input.value = typed + starts[0].label.slice(typed.length);
        input.setSelectionRange(typed.length, input.value.length);
        idx = 0; box.firstChild.classList.add('on');
      }
    };
    /* the Sketcher's mode. A CAS number (digits and hyphens) lists matching CAS numbers with their names; anything
       with a lowercase letter lists names (so IUPAC names such as "2-methylpropan-1-ol" or "(2S)-butan-2-ol" are
       found too); text with no lowercase letter is taken for a SMILES ("CCO", "CO": no list, never completed into a
       name). Completion fills in the best match and leaves the added part selected: keep typing to replace it, Enter
       to take it, the arrow keys to choose another. */
    function searchSuggest(raw, q, e) {
      const I = searchIndex(), isCas = /^\d[\d-]*$/.test(q);
      if (!isCas && (q.length < 2 || !/[a-z]/.test(raw) || /[\[\]=#@\\\/]/.test(raw))) { close(); return; }
      const starts = [], contains = [];
      for (const n of (isCas ? I.cas : I.names)) {
        if (n.k.startsWith(q)) starts.push(n); else if (!isCas && q.length >= 3 && n.k.includes(q)) contains.push(n);
        if (starts.length >= 300) break;
      }
      starts.sort((a, b) => a.label.length - b.label.length || a.label.localeCompare(b.label));
      const hits = starts.concat(contains).slice(0, 12);
      if (!hits.length) { close(); return; }
      // a name from the common-name list carries no formula or CAS: work them out for the rows shown (once each)
      for (const h of hits) if (!isCas && !h.formula && !h.done) {
        h.done = true;
        try { const g2 = C.searchMolecule(h.label).graph; h.formula = C.formula(g2) || ''; h.cas = (C.casFor && C.casFor(g2)) || ''; } catch (err) {}
      }
      box.innerHTML = hits.map(h => isCas
        ? `<div data-name="${esc(h.label)}">${esc(h.label)} <em>${esc(h.name)}</em><span>${esc(h.formula)}</span></div>`
        : `<div data-name="${esc(h.label)}">${esc(h.label)}<span>${esc([h.formula, h.cas ? 'CAS ' + h.cas : ''].filter(Boolean).join(' · '))}</span></div>`).join('');
      box.hidden = false; idx = -1;
      const typing = e && e.type === 'input' && e.inputType && !e.inputType.startsWith('delete');
      if (autofill && typing && starts.length && starts[0].k !== q) {
        input.value = raw + starts[0].label.slice(raw.length);
        input.setSelectionRange(raw.length, input.value.length);
        idx = 0; box.firstChild.classList.add('on');
      }
    }
    input.addEventListener('input', suggest);
    input.addEventListener('focus', suggest);
    input.addEventListener('blur', () => setTimeout(close, 150));
    input.addEventListener('keydown', e => {
      if (e.key === 'ArrowDown' && !box.hidden) { e.preventDefault(); move(1); }
      else if (e.key === 'ArrowUp' && !box.hidden) { e.preventDefault(); move(-1); }
      else if (e.key === 'Escape') close();
      else if (e.key === 'Enter') {
        const on = box.querySelector('.on');
        if (on) input.value = on.dataset.name;
        close(); onSubmit(input.value);
      }
    });
    box.addEventListener('pointerdown', e => {
      const it = e.target.closest('[data-name]'); if (!it) return;
      e.preventDefault(); input.value = it.dataset.name; close(); onSubmit(input.value);
    });
    return { close };
  }
  const importSuggest = attachSuggest(importBox, document.getElementById('sketchSuggest'), importText, null, { search: true });
  const closeSuggest = () => importSuggest.close();
  function download(name, blob) {
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob); a.download = name; a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }
  /* a standalone SVG of the drawing: cropped to the molecule, styles inlined */
  function svgText() {
    const saveHover = hover, saveSel = selected;
    hover = { atom: null, bond: null }; selected = new Set(); render({});
    const pad = 24;
    const xs = g.atoms.map(a => a.x), ys = g.atoms.map(a => a.y);
    const x0 = Math.min(...xs) - pad, y0 = Math.min(...ys) - pad, w = Math.max(...xs) - x0 + pad, h = Math.max(...ys) - y0 + pad;
    const clone = svg.cloneNode(true);
    clone.removeAttribute('id'); clone.removeAttribute('style');
    clone.setAttribute('xmlns', NS); clone.setAttribute('viewBox', `${x0} ${y0} ${w} ${h}`);
    clone.setAttribute('width', w); clone.setAttribute('height', h);
    const style = document.createElementNS(NS, 'style');
    style.textContent = '.bond line{stroke:#111;stroke-width:2;stroke-linecap:round}.wedge{fill:#111}.hash{stroke:#111;stroke-width:2}' +
      '.label{font:500 18px "Segoe UI",system-ui,sans-serif;text-anchor:middle;dominant-baseline:central;fill:#111}.label-bg{fill:#fff}' +
      '.hcount{font-size:13px}.charge{font-size:12px;baseline-shift:super}';
    clone.prepend(style);
    const bg = el('rect', { x: x0, y: y0, width: w, height: h, fill: '#fff' }); clone.insertBefore(bg, style.nextSibling);
    hover = saveHover; selected = saveSel; render();
    return { text: new XMLSerializer().serializeToString(clone), w, h };
  }
  function doExport() {
    if (!g.atoms.length) return;
    const fmt = formatBox.value;
    const txt = (name, s, type = 'text/plain') => download(name, new Blob([s], { type }));
    switch (fmt) {
      case 'smiles': txt('molecule.smi', currentSmiles() + '\n'); break;
      case 'name': txt('molecule-name.txt', currentName() + '\n'); break;
      case 'common': {
        // only a name the app's library or common-name list gives for exactly this structure
        let n = ''; try { n = (C.commonNameFor && C.commonNameFor(g)) || ((C.libEntryFor && C.libEntryFor(g)) || {}).name || ''; } catch (e) {}
        if (!n) { alert('The app knows no common name for this structure.'); break; }
        txt('molecule-common-name.txt', n + '\n'); break;
      }
      case 'cas': {
        let n = ''; try { n = C.casFor ? C.casFor(g) : ''; } catch (e) {}
        if (!n) { alert('This structure is not among the CAS numbers in the app\'s offline library.'); break; }
        txt('molecule-cas.txt', n + '\n'); break;
      }
      case 'formula': txt('molecule-formula.txt', C.formula(g) + '\n'); break;
      case 'mol': { try { txt('molecule.mol', C.toMolfile(g, 'sketch'), 'chemical/x-mdl-molfile'); } catch (e) { alert(e.message); } break; }
      case 'svg': txt('molecule.svg', svgText().text, 'image/svg+xml'); break;
      case 'png': {
        const { text, w, h } = svgText();
        const img = new Image();
        img.onload = () => {
          const cv = document.createElement('canvas'); cv.width = w * 2; cv.height = h * 2;
          const ctx = cv.getContext('2d'); ctx.scale(2, 2); ctx.drawImage(img, 0, 0);
          cv.toBlob(b => download('molecule.png', b), 'image/png');
        };
        img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(text);
        break;
      }
    }
  }
  setAutoTidy(autoTidy);                           // the button shows the remembered state
  bar.addEventListener('click', e => {
    const b = e.target.closest('button'); if (!b) return;
    const edit = fn => { if (!g.atoms.length) return; snapshot(); fn(); render(); };
    switch (b.dataset.act) {
      case 'undo': undo(); break;
      case 'redo': redo(); break;
      case 'autotidy': setAutoTidy(!autoTidy); if (autoTidy) edit(cleanLayout); break;
      case 'clear': if (g.atoms.length) { snapshot(); g = { atoms: [], bonds: [], nextId: 1 }; selected.clear(); render(); } break;
      case 'mirrorx': edit(() => mirror('x')); break;
      case 'mirrory': edit(() => mirror('y')); break;
      case 'addH': edit(addHydrogens); break;
      case 'stripH': edit(stripHydrogens); break;
      case 'import': importText(importBox.value); break;
      case 'export': doExport(); break;
      case 'copy': { const s = currentSmiles(); if (s) navigator.clipboard.writeText(s).then(() => flash(b, 'Copied')); break; }
    }
  });
  fileBox.addEventListener('change', () => {
    const f = fileBox.files[0]; if (!f) return;
    f.text().then(t => { importText(t); fileBox.value = ''; });
  });
  function flash(btn, msg) { const t = btn.textContent; btn.textContent = msg; setTimeout(() => (btn.textContent = t), 900); }

  /* ================= pointer handling ================= */
  let drag = null;
  const pos = e => { const r = svg.getBoundingClientRect(); return { x: e.clientX - r.left, y: e.clientY - r.top }; };

  svg.addEventListener('pointerdown', e => {
    if (e.button !== 0 || !tool) return;
    closePopup();
    const { x, y } = pos(e);
    const a = atomAt(x, y), b = a ? null : bondAt(x, y);
    const k = tool.kind;
    drag = { x, y, moved: false, from: a, bond: b };
    if (k === 'select') {
      if (a) {
        if (!selected.has(a.id)) { if (!e.shiftKey) selected.clear(); selected.add(a.id); }
        drag.group = [...selected].map(atom);
      } else if (b) {
        if (!e.shiftKey) selected.clear();
        selected.add(b.a); selected.add(b.b);
        drag.group = [...selected].map(atom);
      } else {
        if (!e.shiftKey) selected.clear();
        drag.marquee = true;
      }
      if (drag.group) drag.start = drag.group.map(p => ({ x: p.x, y: p.y }));
      render();
    } else if (k === 'rotate') {
      drag.group = targetAtoms();
      if (!drag.group.length) { drag = null; return; }
      drag.start = drag.group.map(p => ({ x: p.x, y: p.y }));
      drag.centre = centroid(drag.group);
      drag.angle0 = Math.atan2(y - drag.centre.y, x - drag.centre.x);
    }
    svg.setPointerCapture(e.pointerId);
  });

  svg.addEventListener('pointermove', e => {
    const { x, y } = pos(e);
    if (!drag) {
      const a = atomAt(x, y), b = a ? null : bondAt(x, y);
      if (a !== hover.atom || b !== hover.bond) { hover = { atom: a, bond: b }; render(); }
      return;
    }
    if (!drag.moved && dist(x, y, drag.x, drag.y) < 4) return;
    drag.moved = true;
    const k = tool.kind;
    if (k === 'select') {
      if (drag.marquee) render({ marquee: { x1: drag.x, y1: drag.y, x2: x, y2: y } });
      else {
        const dx = x - drag.x, dy = y - drag.y;
        drag.group.forEach((p, i) => { p.x = drag.start[i].x + dx; p.y = drag.start[i].y + dy; });
        if (!drag.snapped) { snapshot(); drag.snapped = true; }
        render({});
      }
    } else if (k === 'rotate') {
      let th = Math.atan2(y - drag.centre.y, x - drag.centre.x) - drag.angle0;
      th = Math.round(th / (Math.PI / 36)) * (Math.PI / 36);          // 5° steps
      const cos = Math.cos(th), sin = Math.sin(th), c = drag.centre;
      drag.group.forEach((p, i) => { const sx = drag.start[i].x - c.x, sy = drag.start[i].y - c.y; p.x = c.x + sx * cos - sy * sin; p.y = c.y + sx * sin + sy * cos; });
      if (!drag.snapped) { snapshot(); drag.snapped = true; }
      render({});
    } else if (k === 'chain') {
      const s = drag.from ? drag.from : { x: snapGrid(drag.x), y: snapGrid(drag.y) };
      const pts = chainPoints(s.x, s.y, x, y);
      render({ points: pts, ghost: pts.slice(1) });
    } else if (k === 'el') {
      const r = elDrag(drag, x, y);
      if (r.snapTo) render({ bond: { a: drag.from.id, order: 1, _end: { x: r.snapTo.x, y: r.snapTo.y } } });
      else render({ points: r.pts, ghost: r.pts.slice(1) });
    } else if ((k === 'bond' || k === 'wedge' || k === 'hash') && drag.from) {
      const target = atomAt(x, y);
      const end = target && target !== drag.from ? { x: target.x, y: target.y } : snapDir(drag.from, x, y);
      render({ bond: { a: drag.from.id, order: k === 'bond' ? +tool.arg : 1, _end: end } });
    }
  });

  svg.addEventListener('pointerup', e => {
    if (!drag) return;
    const { x, y } = pos(e);
    const d = drag; drag = null;
    const k = tool.kind;
    if (k === 'select') {
      if (d.marquee && d.moved) {
        const x1 = Math.min(d.x, x), x2 = Math.max(d.x, x), y1 = Math.min(d.y, y), y2 = Math.max(d.y, y);
        g.atoms.forEach(a => { if (a.x >= x1 && a.x <= x2 && a.y >= y1 && a.y <= y2) selected.add(a.id); });
      }
      render(); return;
    }
    if (k === 'rotate') { render(); return; }
    const rings = ringCount();
    guarded(() => editAt(d, x, y, k));
    if (autoTidy && k !== 'erase' && ringCount() > rings) {   // Auto: Tidy is on and a ring was just closed: snap it
      snapshot();
      try { snapRingMolecule(); } catch (e) { undoStack.pop(); }
      render();
    }
  });
  function editAt(d, x, y, k) {
    if (d.moved) {
      if (k === 'chain') {
        snapshot();
        const s = d.from;
        const pts = chainPoints(s ? s.x : snapGrid(d.x), s ? s.y : snapGrid(d.y), x, y);
        let prev = s || addAtom('C', pts[0].x, pts[0].y);
        for (let i = 1; i < pts.length; i++) {
          const hit = atomAt(pts[i].x, pts[i].y);
          const n = hit && hit !== prev ? hit : addAtom('C', pts[i].x, pts[i].y);
          addBond(prev.id, n.id, 1); prev = n;
        }
      } else if (k === 'el') {
        snapshot();
        const r = elDrag(d, x, y);
        if (r.snapTo) addBond(d.from.id, r.snapTo.id, 1);
        else {
          let prev = d.from || addAtom(tool.arg, r.pts[0].x, r.pts[0].y);
          for (let i = 1; i < r.pts.length; i++) {
            const hit = atomAt(r.pts[i].x, r.pts[i].y);
            const n = hit && hit !== prev ? hit : addAtom(tool.arg, r.pts[i].x, r.pts[i].y);
            addBond(prev.id, n.id, 1); prev = n;
          }
        }
      } else if ((k === 'bond' || k === 'wedge' || k === 'hash') && d.from) {
        snapshot();
        const target = atomAt(x, y);
        const order = k === 'bond' ? +tool.arg : 1;
        const stereo = k === 'wedge' || k === 'hash' ? k : null;
        if (target && target !== d.from) addBond(d.from.id, target.id, order, stereo);
        else { const p = snapDir(d.from, x, y); const n = atomAtOrNew(d.from, p.x, p.y); addBond(d.from.id, n.id, order, stereo); }
      }
      render(); return;
    }
    click(x, y, d.from, d.bond);
  }
  svg.addEventListener('pointerleave', () => { if (!drag && (hover.atom || hover.bond)) { hover = { atom: null, bond: null }; render(); } });

  function click(x, y, a, b) {
    const k = tool.kind;
    switch (k) {
      case 'el':
        snapshot();
        if (a) a.element = tool.arg; else addAtom(tool.arg, snapGrid(x), snapGrid(y));
        break;
      case 'bond': {
        const order = +tool.arg;
        snapshot();
        // the single-bond tool on an existing bond raises its order one step per click: – → = → ≡ (a triple stays triple)
        if (b) { b.order = order === 1 ? Math.min(3, b.order + 1) : order; if (b.order !== 1) { delete b.stereo; delete b.narrow; } }
        else if (a) { const p = polar(a.x, a.y, freeAngle(a), BOND); const n = atomAtOrNew(a, p.x, p.y); addBond(a.id, n.id, order); }
        else { const p = addAtom('C', snapGrid(x), snapGrid(y)); const q = addAtom('C', p.x + BOND, p.y); addBond(p.id, q.id, order); }
        break;
      }
      case 'wedge': case 'hash':
        snapshot();
        if (b) {
          if (b.stereo === k) b.narrow = otherEnd(b, b.narrow);
          else { b.order = 1; b.stereo = k; b.narrow = dist(x, y, atom(b.a).x, atom(b.a).y) < dist(x, y, atom(b.b).x, atom(b.b).y) ? b.a : b.b; }
        } else if (a) { const p = polar(a.x, a.y, freeAngle(a), BOND); const n = atomAtOrNew(a, p.x, p.y); addBond(a.id, n.id, 1, k); }
        break;
      case 'charge':
        if (a) { snapshot(); a.charge = (a.charge || 0) === +tool.arg ? 0 : +tool.arg; }
        break;
      case 'radical':
        if (a) { snapshot(); if (a.radical) delete a.radical; else a.radical = true; }
        break;
      case 'ring': {
        const n = +tool.arg, aromatic = tool.arg2 === '1';
        snapshot();
        if (b) ringOnBond(n, b, aromatic); else if (a) ringOnAtom(n, a, aromatic); else ringAtPoint(n, snapGrid(x), snapGrid(y), aromatic);
        break;
      }
      case 'frag': case 'paste': {
        const frag = fragmentGraph(tool.arg); if (!frag) break;
        snapshot();
        if (a) pasteFragment(frag, 0, 0, a);
        else pasteFragment(frag, snapGrid(x), snapGrid(y), null, k === 'paste');
        break;
      }
      case 'chain':
        snapshot();
        if (!a) addAtom('C', snapGrid(x), snapGrid(y));
        break;
      case 'erase':
        if (a && onHydrogenOf(a, x)) { snapshot(); a.charge = (a.charge || 0) - 1; }      // OH → O⁻: the H leaves as H⁺
        else if (a) { snapshot(); removeAtom(a.id); }                                    // the atom with its H (the whole OH)
        else if (b) {                                                                   // one bond order per click
          snapshot();
          if (b.order > 1) { b.order -= 1; delete b.ezUnspec; } else removeBond(b);
        }
        break;
    }
    hover = { atom: null, bond: null };
    render();
  }

  /* ================= keyboard ================= */
  const KEY_EL = { c: 'C', n: 'N', o: 'O', s: 'S', p: 'P', f: 'F', i: 'I', h: 'H', l: 'Cl', b: 'Br' };
  document.addEventListener('keydown', e => {
    if (e.target.tagName === 'INPUT' || e.target.tagName === 'SELECT') return;
    if (!document.getElementById('page-sketcher').classList.contains('active')) return;
    const key = e.key.toLowerCase();
    if (e.ctrlKey && key === 'z') { e.preventDefault(); e.shiftKey ? redo() : undo(); return; }
    if (e.ctrlKey && key === 'y') { e.preventDefault(); redo(); return; }
    if (e.ctrlKey && key === 'a') { e.preventDefault(); pickTool('select'); g.atoms.forEach(a => selected.add(a.id)); render(); return; }
    if (e.key === 'Escape') { drag = null; selected.clear(); closePopup(); render(); return; }
    if (e.key === 'Delete' || e.key === 'Backspace') {
      if (selected.size) { snapshot(); [...selected].forEach(removeAtom); render(); }
      else if (hover.atom) { snapshot(); removeAtom(hover.atom.id); hover.atom = null; render(); }
      else if (hover.bond) { snapshot(); removeBond(hover.bond); hover.bond = null; render(); }
      return;
    }
    if (hover.atom && KEY_EL[key] && !e.ctrlKey) { guarded(() => { snapshot(); hover.atom.element = KEY_EL[key]; render(); }); return; }
    if (hover.bond && /^[123]$/.test(e.key)) guarded(() => { snapshot(); hover.bond.order = +e.key; if (hover.bond.order > 1) delete hover.bond.stereo; render(); });
  });

  render();
  window.Sketcher = {
    get graph() { return g; },
    set graph(v) { snapshot(); g = v; selected.clear(); render(); },
    smiles: currentSmiles, name: currentName, load: loadGraph, render, attachSuggest, nameIndex,
  };
})();
