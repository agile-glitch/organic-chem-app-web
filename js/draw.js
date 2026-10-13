/* Draw — a canvas for reaction schemes made of molecule boxes, text, arrows and electron dots.
   A port of the Draw page of the AgileCHEM desktop program (draw_page.py, draw_items.py): the same tools, the same
   behaviour, and one tab (it is not numbered Draw1, Draw2: there is only one canvas).

   Tools     Select (rectangle over empty space, or over the blank part of a molecule box; Ctrl/Shift adds), Molecule box,
             Text (x₂ / x² for subscripts and superscripts), → ⇌ ↔ arrows, • electron, •• lone pair.
   Molecule  double-click a box: type a name, IUPAC name, CAS number or SMILES (suggestions as you type), or draw it in
             the Sketcher and send it back. Display options: hydrogens, common name, IUPAC name, SMILES, CAS.
             Selected box: drag a border to resize, drag just outside the molecule to turn it (Shift = 15° steps).
   Arrows    drag an end to move it (sticks straight near horizontal; Shift = 15° steps; Alt = free), drag the middle to bend.
   Colour    Color ▾ colours the selection: lassoed atoms/bonds, whole molecules, text, arrows, electrons.
   Keys      Delete / Backspace removes the selection, Esc returns to Select. Ctrl + wheel zooms. Export PNG at 4x.

   Molecules are laid out by Chem (js/chem.js, stereo kept by MolDraw.layoutKeepingStereo); names, IUPAC names and CAS
   numbers come from the same offline library as the Sketcher, so nothing here needs the internet. */
(() => {
  'use strict';
  const C = window.Chem, page = document.getElementById('page-draw');
  if (!C || !page) return;

  const NS = 'http://www.w3.org/2000/svg';
  const SCENE_W = 3200, SCENE_H = 2200;
  const BLUE = '#4a86e8', FONT = 'Arial, Helvetica, sans-serif';
  const EDGE = 7, MOVE_BOND = 16, MOVE_ATOM = 20, ROT_NEAR = 34, MIN_W = 90, MIN_H = 70, PAD = 30, CAPTION_PX = 14, TEXT_PX = 16;
  const store = { get(k) { try { return JSON.parse(localStorage.getItem(k)); } catch (e) { return null; } }, set(k, v) { try { localStorage.setItem(k, JSON.stringify(v)); } catch (e) { /* private window */ } } };
  const SHOW = Object.assign({ hydrogens: false, common: false, iupac: false, smiles: false, cas: false }, store.get('draw.show') || {});

  const h = (tag, attrs, parent) => {
    const n = document.createElementNS(NS, tag);
    for (const k in (attrs || {})) n.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(n);
    return n;
  };
  const esc = s => String(s ?? '').replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
  const dist = (a, b) => Math.hypot(a.x - b.x, a.y - b.y);
  const manhattan = (a, b) => Math.abs(a.x - b.x) + Math.abs(a.y - b.y);
  const segDist = (p, a, b) => {                                   // distance from p to the segment ab
    const vx = b.x - a.x, vy = b.y - a.y, l2 = vx * vx + vy * vy || 1;
    const t = Math.max(0, Math.min(1, ((p.x - a.x) * vx + (p.y - a.y) * vy) / l2));
    return Math.hypot(p.x - (a.x + vx * t), p.y - (a.y + vy * t));
  };
  const inRect = (p, r) => p.x >= r.x && p.x <= r.x + r.w && p.y >= r.y && p.y <= r.y + r.h;
  const rectsMeet = (a, b) => a.x <= b.x + b.w && b.x <= a.x + a.w && a.y <= b.y + b.h && b.y <= a.y + a.h;
  const normRect = (a, b) => ({ x: Math.min(a.x, b.x), y: Math.min(a.y, b.y), w: Math.abs(a.x - b.x), h: Math.abs(a.y - b.y) });
  let measureCtx = null;
  const textWidth = (s, px) => { if (!measureCtx) measureCtx = document.createElement('canvas').getContext('2d'); measureCtx.font = px + 'px ' + FONT; return measureCtx.measureText(s).width; };

  /* ================= items ================= */
  class Item {
    constructor() { this.g = h('g'); this.selected = false; this.scene = null; }
    setSelected(on) { if (on !== this.selected) { this.selected = on; this.hover = null; if (!on && this.selAtoms) { this.selAtoms = new Set(); this.selBonds = new Set(); } this.render(); } }
    press() { return null; }                                       // a drag of the item's own (resize, turn, handle), or null: a plain move
    hoverAt() { return null; }                                     // a CSS cursor, or null
    exportNode() { const was = this.selected; this.selected = false; this.render(); const n = this.g.cloneNode(true); this.selected = was; this.render(); return n; }
  }

  /* ---- a molecule in a box ---- */
  class MoleculeBox extends Item {
    constructor(w = 220, hgt = 160) {
      super();
      this.x = 0; this.y = 0; this.w = w; this.h = hgt;
      this.smiles = ''; this.graph = null; this.lay = null; this.angle = 0; this.autoFit = true;
      this.atomColors = new Map(); this.bondColors = new Map(); this.selAtoms = new Set(); this.selBonds = new Set();
      this.atomXY = []; this.bondPairs = []; this.hover = null; this.rotating = false; this.rotDelta = 0;
    }
    /* ---- the molecule: laid out once; turned by `angle` and centred in the box on every draw */
    setMolecule(smiles) {
      if (smiles !== this.smiles) { this.angle = 0; this.autoFit = true; }
      this.smiles = smiles || ''; this.atomColors = new Map(); this.bondColors = new Map(); this.clearLasso();
      this.graph = null; this.lay = null;
      if (this.smiles) { try { this.graph = C.parseSmiles(this.smiles); } catch (e) { this.graph = null; } }
      if (this.graph && this.graph.atoms.length) {
        const g = { atoms: this.graph.atoms.map(a => ({ ...a })), bonds: this.graph.bonds.map(b => ({ ...b })), nextId: this.graph.nextId };
        if (g.bonds.length && window.MolDraw) window.MolDraw.layoutKeepingStereo(g); else C.layoutGraph(g);
        this.lay = g;
      } else this.graph = null;
      this.refreshDisplay();
    }
    get canonical() { try { return this.graph ? C.toSmiles(this.graph) : this.smiles; } catch (e) { return this.smiles; } }
    info() {
      const g = this.graph, out = {};
      if (!g) return out;
      if (SHOW.common) { try { const n = C.commonNameFor ? C.commonNameFor(g) : C.displayNameFor(g); if (n) out.common = n; } catch (e) { /* none */ } }
      if (SHOW.iupac) { try { const n = C.iupacName(g); out.iupac = n && !/^\s*\(?(?:could not|cannot|unnamed|\?)/i.test(n) ? n : ''; } catch (e) { out.iupac = ''; } }
      if (SHOW.cas) { try { out.cas = C.casFor ? C.casFor(g) : null; } catch (e) { out.cas = null; } }
      return out;
    }
    captionLines() {
      if (!this.graph) return [];
      const inf = this.info(), lines = [];
      if (SHOW.common && inf.common) lines.push(inf.common);
      if (SHOW.iupac) lines.push(inf.iupac || 'IUPAC name not found');
      if (SHOW.smiles) lines.push(this.canonical);
      if (SHOW.cas) lines.push('CAS ' + (inf.cas || 'not found'));
      return lines;
    }
    wrapped(width) {                                               // caption lines broken anywhere to fit the width
      const out = [];
      for (const t of this.captionLines()) {
        let cur = '';
        for (const ch of t) { if (cur && textWidth(cur + ch, CAPTION_PX) > width) { out.push(cur); cur = ch; } else cur += ch; }
        out.push(cur);
      }
      return out;
    }
    captionHeight(width) { const n = this.wrapped(Math.max((width ?? this.w) - 16, 40)).length; return n ? 8 + n * (CAPTION_PX + 3) : 0; }
    /* the laid-out atoms, turned and centred: [{x, y}] in box coordinates, for a box of this size */
    place(width, height) {
      const lay = this.lay; if (!lay) return null;
      const a = this.angle * Math.PI / 180, cos = Math.cos(a), sin = Math.sin(a);
      const pts = lay.atoms.map(at => ({ x: at.x * cos - at.y * sin, y: at.x * sin + at.y * cos }));
      const xs = pts.map(p => p.x), ys = pts.map(p => p.y);
      const minx = Math.min(...xs), maxx = Math.max(...xs), miny = Math.min(...ys), maxy = Math.max(...ys);
      return { pts, w: maxx - minx, h: maxy - miny, cx: (minx + maxx) / 2, cy: (miny + maxy) / 2, to: (cx, cy) => pts.map(p => ({ x: cx + p.x - (minx + maxx) / 2, y: cy + p.y - (miny + maxy) / 2 })) };
    }
    neededSize() {
      const pl = this.place(); if (!pl) return [this.w, this.h];
      const w = Math.max(MIN_W, pl.w + 2 * PAD);
      return [w, Math.max(MIN_H, pl.h + 2 * PAD + this.captionHeight(w))];
    }
    resizeAboutCentre(w, hgt) {
      if (Math.abs(w - this.w) < 0.5 && Math.abs(hgt - this.h) < 0.5) return;
      this.x += (this.w - w) / 2; this.y += (this.h - hgt) / 2; this.w = w; this.h = hgt;
    }
    refreshDisplay() {                                             // fit the box to the molecule (or stop it clipping once sized by hand), redraw
      if (this.lay) {
        const [w, hgt] = this.neededSize();
        if (this.autoFit) this.resizeAboutCentre(w, hgt); else this.resizeAboutCentre(Math.max(this.w, w), Math.max(this.h, hgt));
      }
      this.render();
    }
    /* ---- selection of atoms and bonds with the rectangle */
    clearLasso() { if (this.selAtoms.size || this.selBonds.size) { this.selAtoms = new Set(); this.selBonds = new Set(); this.render(); } }
    lasso(rect, additive) {                                        // rect in scene coordinates
      if (!additive) { this.selAtoms = new Set(); this.selBonds = new Set(); }
      this.atomXY.forEach((p, i) => { if (inRect({ x: this.x + p.x, y: this.y + p.y }, rect)) this.selAtoms.add(i); });
      this.bondPairs.forEach(([i, j], k) => { const a = this.atomXY[i], b = this.atomXY[j]; if (inRect({ x: this.x + (a.x + b.x) / 2, y: this.y + (a.y + b.y) / 2 }, rect)) this.selBonds.add(k); });
      this.render();
      return this.selAtoms.size > 0 || this.selBonds.size > 0;
    }
    applyColor(color) {                                            // the lassoed atoms/bonds, or the whole molecule; null = back to black
      const picked = this.selAtoms.size > 0 || this.selBonds.size > 0;
      const atoms = picked ? this.selAtoms : new Set(this.atomXY.keys()), bonds = picked ? this.selBonds : new Set(this.bondPairs.keys());
      for (const [table, keys] of [[this.atomColors, atoms], [this.bondColors, bonds]]) for (const k of keys) { if (color == null) table.delete(k); else table.set(k, color); }
      this.render();
    }
    /* ---- geometry */
    bbox() { return { x: this.x, y: this.y, w: this.w, h: this.h }; }
    pad() { return this.selected ? EDGE + 6 : 3; }
    hit(p) {
      const m = this.pad(), top = this.selected ? 14 : 0;
      return inRect(p, { x: this.x - m, y: this.y - m - top, w: this.w + 2 * m, h: this.h + 2 * m + top });
    }
    nearDrawing(q) {                                               // [distance to the nearest bond, to the nearest atom], q in box coordinates
      let nb = 1e9, na = 1e9;
      this.atomXY.forEach(p => { na = Math.min(na, dist(p, q)); });
      this.bondPairs.forEach(([i, j]) => { nb = Math.min(nb, segDist(q, this.atomXY[i], this.atomXY[j])); });
      return [nb, na];
    }
    zone(p) {                                                      // null in the middle; {resize: [L, R, T, B]} on the border; {rotate} just outside the molecule
      if (!this.selected) return null;
      const x = p.x - this.x, y = p.y - this.y, w = this.w, hh = this.h;
      const onX = x >= -EDGE && x <= w + EDGE, onY = y >= -EDGE && y <= hh + EDGE;
      let left = Math.abs(x) <= EDGE && onY, right = Math.abs(x - w) <= EDGE && onY, top = Math.abs(y) <= EDGE && onX, bottom = Math.abs(y - hh) <= EDGE && onX;
      for (const [cx, cy, f] of [[0, 0, 'LT'], [w, 0, 'RT'], [0, hh, 'LB'], [w, hh, 'RB']]) {   // generous corners
        if (Math.abs(x - cx) <= 12 && Math.abs(y - cy) <= 12) { left = left || f.includes('L'); right = right || f.includes('R'); top = top || f.includes('T'); bottom = bottom || f.includes('B'); }
      }
      if (left || right || top || bottom) return { resize: [left, right, top, bottom] };
      if (this.atomXY.length) { const [nb, na] = this.nearDrawing({ x, y }); if (nb > MOVE_BOND && na > MOVE_ATOM && Math.min(nb, na - 4) <= ROT_NEAR) return { rotate: true }; }
      return null;
    }
    glyphHit(p) {                                                  // does this point grab the box: an atom, a bond, (selected) its border / rotate band, or an empty box
      if (this.zone(p)) return true;
      if (!this.atomXY.length) return true;
      const [nb, na] = this.nearDrawing({ x: p.x - this.x, y: p.y - this.y });
      return na <= MOVE_ATOM || nb <= MOVE_BOND;
    }
    hoverAt(p) {
      const z = this.zone(p), key = z ? JSON.stringify(z) : '';
      if (key !== (this.hover ? JSON.stringify(this.hover) : '')) { this.hover = z; this.render(); }
      if (!z) return null;
      if (z.rotate) return 'grab';
      const [l, r, t, b] = z.resize;
      return (l && t) || (r && b) ? 'nwse-resize' : (r && t) || (l && b) ? 'nesw-resize' : l || r ? 'ew-resize' : 'ns-resize';
    }
    leave() { if (this.hover) { this.hover = null; this.render(); } }
    press(p, e) {
      const z = this.zone(p); if (!z) return null;
      if (z.resize) {
        const x0 = this.x, y0 = this.y, w0 = this.w, h0 = this.h, [l, r, t, b] = z.resize;
        return { move: q => {
          const dx = q.x - p.x, dy = q.y - p.y;
          let x = x0, y = y0, w = w0, hh = h0;
          if (r) w = Math.max(MIN_W, w0 + dx);
          if (b) hh = Math.max(MIN_H, h0 + dy);
          if (l) { w = Math.max(MIN_W, w0 - dx); x = x0 + (w0 - w); }
          if (t) { hh = Math.max(MIN_H, h0 - dy); y = y0 + (h0 - hh); }
          this.x = x; this.y = y; this.w = w; this.h = hh; this.autoFit = false; this.render();
        }, up: () => {} };
      }
      /* turn: the molecule is turned as a picture about the box centre while the mouse is down, and laid out again on release */
      const before = [this.w, this.h];
      this.makeRoomToTurn();
      const centre = { x: this.x + this.w / 2, y: this.y + (this.h - this.captionHeight()) / 2 };
      const a0 = Math.atan2(p.y - centre.y, p.x - centre.x) * 180 / Math.PI, base = this.angle;
      this.rotating = true; this.rotDelta = 0; this.render();
      return { move: (q, ev) => {
        let delta = Math.atan2(q.y - centre.y, q.x - centre.x) * 180 / Math.PI - a0;
        if (ev.shiftKey) delta = Math.round((base + delta) / 15) * 15 - base;
        this.rotDelta = delta; this.render();
      }, up: () => {
        this.rotating = false; this.angle = (((base + this.rotDelta + 180) % 360) + 360) % 360 - 180; this.rotDelta = 0;
        if (!this.autoFit) this.resizeAboutCentre(before[0], before[1]);
        this.refreshDisplay();
      } };
    }
    makeRoomToTurn() {                                             // grow about the centre until the molecule fits at any angle
      if (!this.atomXY.length) return;
      const xs = this.atomXY.map(q => q.x), ys = this.atomXY.map(q => q.y), bx = (Math.min(...xs) + Math.max(...xs)) / 2, by = (Math.min(...ys) + Math.max(...ys)) / 2;
      const radius = Math.max(...this.atomXY.map(q => Math.hypot(q.x - bx, q.y - by))), side = 2 * (radius + PAD);
      this.resizeAboutCentre(Math.max(this.w, side), Math.max(this.h, side + this.captionHeight()));
      this.render();
    }
    /* ---- drawing */
    render() {
      const g = this.g; g.textContent = ''; g.setAttribute('transform', `translate(${this.x},${this.y})`);
      const capH = this.captionHeight(), hasMol = !!this.lay, turning = this.rotating && Math.abs(this.rotDelta) > 1e-6;
      this.atomXY = []; this.bondPairs = [];
      if (!hasMol) {
        h('rect', { x: 0, y: 0, width: this.w, height: this.h, rx: 4, fill: 'none', stroke: '#b4b4bc', 'stroke-width': 1.5, 'stroke-dasharray': '6 4' }, g);
        const t = h('text', { x: this.w / 2, y: this.h / 2, 'text-anchor': 'middle', 'dominant-baseline': 'middle', fill: '#9a9aa2', 'font-size': 13, 'font-family': FONT }, g);
        t.textContent = 'Double-click to add a molecule';
      } else {
        const pl = this.place(), cx = this.w / 2, cy = (this.h - capH) / 2, pts = pl.to(cx, cy);
        this.atomXY = pts; this.bondPairs = this.lay.bonds.map(b => [this.lay.atoms.findIndex(a => a.id === b.a), this.lay.atoms.findIndex(a => a.id === b.b)]);
        if (this.selected && ((this.hover && this.hover.rotate) || this.rotating)) {   // the band that turns the molecule when dragged: one even shape
          const band = h('g', { fill: BLUE, stroke: BLUE, 'stroke-linecap': 'round', opacity: 0.15 }, g);
          this.bondPairs.forEach(([i, j]) => h('line', { x1: pts[i].x, y1: pts[i].y, x2: pts[j].x, y2: pts[j].y, 'stroke-width': 2 * ROT_NEAR }, band));
          pts.forEach(q => h('circle', { cx: q.x, cy: q.y, r: ROT_NEAR, stroke: 'none' }, band));
        }
        const mol = h('g', {}, g);
        if (turning) mol.setAttribute('transform', `rotate(${this.rotDelta} ${cx} ${cy})`);
        const halo = h('g', { fill: 'rgba(74,134,232,.35)', stroke: 'rgba(74,134,232,.35)', 'stroke-linecap': 'round' }, mol);
        this.selBonds.forEach(k => { const [i, j] = this.bondPairs[k] || []; if (i != null) h('line', { x1: pts[i].x, y1: pts[i].y, x2: pts[j].x, y2: pts[j].y, 'stroke-width': 11 }, halo); });
        this.selAtoms.forEach(i => { if (pts[i]) h('circle', { cx: pts[i].x, cy: pts[i].y, r: 11, stroke: 'none' }, halo); });
        this.drawStructure(mol, pts);
        if (capH) {
          const t = h('text', { 'text-anchor': 'middle', fill: '#333', 'font-size': CAPTION_PX, 'font-family': FONT }, g), lines = this.wrapped(Math.max(this.w - 16, 40));
          lines.forEach((ln, k) => { const s = h('tspan', { x: this.w / 2, y: this.h - capH + 8 + (k + 1) * (CAPTION_PX + 3) - 4 }, t); s.textContent = ln; });
        }
      }
      if (this.selected) this.drawSelection(g);
    }
    drawStructure(parent, pts) {
      const lay = this.lay, atoms = lay.atoms, bonds = lay.bonds, idx = new Map(atoms.map((a, i) => [a.id, i]));
      const bondsOf = id => bonds.filter(b => b.a === id || b.b === id);
      const nbrs = id => bondsOf(id).map(b => atoms[idx.get(b.a === id ? b.b : b.a)]);
      const px = 17, LW = 2, hCount = a => (a.element === 'H' ? 0 : C.implicitH(lay, a));
      const labelled = a => a.element !== 'C' || a.charge || a.radical || bondsOf(a.id).length === 0 || (SHOW.hydrogens && hCount(a) > 0);
      const hWest = (a, i) => { const ns = nbrs(a.id); return ns.filter(o => pts[idx.get(o.id)].x > pts[i].x + 1).length > ns.filter(o => pts[idx.get(o.id)].x < pts[i].x - 1).length; };   // neighbours mostly to the east: the hydrogens go west
      const hasH = a => hCount(a) > 0 && (a.element !== 'C' || SHOW.hydrogens || bondsOf(a.id).length === 0);
      const box = (a, i) => {                                      // half-extents of the label, from the atom: the bond stops at its outline
        const ew = a.element.length * px * 0.62 / 2 + 2, hw = hasH(a) ? px * 0.62 + (hCount(a) > 1 ? px * 0.45 : 0) : 0, cw = a.charge ? px * 0.7 : 0;
        const west = hWest(a, i);
        return { l: ew + (west ? hw : 0), r: ew + (west ? 0 : hw) + cw, v: px * 0.58 };
      };
      const trimAlong = (a, i, ux, uy) => { const bx = box(a, i), ex = ux >= 0 ? bx.r : bx.l; return Math.min(Math.abs(ux) < 1e-6 ? 1e9 : ex / Math.abs(ux), Math.abs(uy) < 1e-6 ? 1e9 : bx.v / Math.abs(uy)) + 1; };
      bonds.forEach((b, k) => {
        const ia = idx.get(b.a), ib = idx.get(b.b), P = pts[ia], Q = pts[ib], A = atoms[ia], B = atoms[ib];
        const dx = Q.x - P.x, dy = Q.y - P.y, L = Math.hypot(dx, dy) || 1, ux = dx / L, uy = dy / L, nx = -uy, ny = ux;
        const tp = labelled(A) ? trimAlong(A, ia, ux, uy) : 0, tq = labelled(B) ? trimAlong(B, ib, -ux, -uy) : 0;
        const x1 = P.x + ux * tp, y1 = P.y + uy * tp, x2 = Q.x - ux * tq, y2 = Q.y - uy * tq, col = this.bondColors.get(k) || '#000';
        const grp = h('g', { stroke: col, fill: col, 'stroke-width': LW, 'stroke-linecap': 'round' }, parent);
        if (b.stereo && b.order === 1) {
          const flip = b.narrow === B.id, sx = flip ? x2 : x1, sy = flip ? y2 : y1, ex = flip ? x1 : x2, ey = flip ? y1 : y2, wd = 5.5;
          if (b.stereo === 'wedge') h('polygon', { points: `${sx},${sy} ${ex + nx * wd},${ey + ny * wd} ${ex - nx * wd},${ey - ny * wd}`, 'stroke-width': 1 }, grp);
          else for (let s = 1; s <= 7; s++) { const t = s / 7, ww = wd * t, qx = sx + (ex - sx) * t, qy = sy + (ey - sy) * t; h('line', { x1: qx + nx * ww, y1: qy + ny * ww, x2: qx - nx * ww, y2: qy - ny * ww, 'stroke-width': 1.6 }, grp); }
          return;
        }
        const off = 5;
        let offsets = { 1: [0], 2: [0, off * 1.7], 3: [-off * 1.7, 0, off * 1.7] }[b.order] || [0];
        if (b.order === 2) {
          const others = [...nbrs(A.id), ...nbrs(B.id)].filter(o => o !== A && o !== B);
          if (!others.length) offsets = [-off, off];
          else { const side = others.reduce((t, o) => { const q = pts[idx.get(o.id)]; return t + Math.sign((q.x - P.x) * nx + (q.y - P.y) * ny); }, 0); offsets = [0, side < 0 ? -off * 1.7 : off * 1.7]; }   // the second line goes on the side the ring is on
        }
        for (const o of offsets) { const sh = b.order === 2 && o !== 0 ? 4 : 0; h('line', { x1: x1 + nx * o + ux * sh, y1: y1 + ny * o + uy * sh, x2: x2 + nx * o - ux * sh, y2: y2 + ny * o - uy * sh }, grp); }
      });
      atoms.forEach((a, i) => {
        if (!labelled(a)) return;
        const P = pts[i], col = this.atomColors.get(i) || '#000';
        const hN = hCount(a), showHs = hasH(a), west = hWest(a, i);
        const ew = a.element.length * px * 0.62;
        h('circle', { cx: P.x, cy: P.y, r: px * 0.55, fill: '#fff' }, parent);
        const t = h('text', { x: P.x, y: P.y, 'text-anchor': 'middle', 'dominant-baseline': 'central', 'font-size': px, 'font-family': FONT, fill: col }, parent);
        t.textContent = a.element;
        if (a.charge) { const c = h('text', { x: P.x + ew / 2 + 1, y: P.y - px * 0.35, 'font-size': px * 0.66, 'font-family': FONT, fill: col, 'dominant-baseline': 'central' }, parent); c.textContent = (Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '') + (a.charge > 0 ? '+' : '−'); }
        if (showHs) {
          const hx = west ? P.x - ew / 2 - 1 : P.x + ew / 2 + 1 + (a.charge ? px * 0.7 : 0);
          const ht = h('text', { x: hx, y: P.y, 'text-anchor': west ? 'end' : 'start', 'dominant-baseline': 'central', 'font-size': px, 'font-family': FONT, fill: col }, parent);
          ht.textContent = 'H'; if (hN > 1) { const s = h('tspan', { 'font-size': px * 0.7, dy: px * 0.22 }, ht); s.textContent = String(hN); }
        }
      });
    }
    drawSelection(g) {
      if (this.rotating) {
        const t = h('text', { x: this.w / 2, y: -8, 'text-anchor': 'middle', fill: BLUE, 'font-size': 13, 'font-family': FONT }, g);
        t.textContent = Math.round((((this.angle + this.rotDelta + 180) % 360) + 360) % 360 - 180) + '°'; return;
      }
      h('rect', { x: 0, y: 0, width: this.w, height: this.h, fill: 'none', stroke: BLUE, 'stroke-width': 1.5, 'stroke-dasharray': '6 4' }, g);
      if (this.hover && this.hover.resize) {
        const [l, r, t, b] = this.hover.resize, ln = (x1, y1, x2, y2) => h('line', { x1, y1, x2, y2, stroke: BLUE, 'stroke-width': 3.5, 'stroke-linecap': 'round' }, g);
        if (l) ln(0, 0, 0, this.h); if (r) ln(this.w, 0, this.w, this.h); if (t) ln(0, 0, this.w, 0); if (b) ln(0, this.h, this.w, this.h);
      }
      for (const x of [0, this.w / 2, this.w]) for (const y of [0, this.h / 2, this.h]) if (!(x === this.w / 2 && y === this.h / 2)) h('rect', { x: x - 3.5, y: y - 3.5, width: 7, height: 7, fill: '#fff', stroke: BLUE, 'stroke-width': 1.2 }, g);
    }
    translate(dx, dy) { this.x += dx; this.y += dy; this.g.setAttribute('transform', `translate(${this.x},${this.y})`); }
    meetsRect(r) { return rectsMeet(this.bbox(), r); }
  }

  /* ---- a text box (rich: subscripts and superscripts) ---- */
  class TextBox extends Item {
    constructor(text = '') {
      super();
      this.x = 0; this.y = 0; this.runs = text ? [{ t: text, s: 0 }] : [];   // [{t, s: 0 normal | -1 sub | 1 super}], "\n" between lines
      this.color = '#000'; this.editing = false;
    }
    plain() { return this.runs.map(r => r.t).join(''); }
    html() {
      return this.runs.map(r => esc(r.t).replace(/\n/g, '<br>')).map((t, i) => { const s = this.runs[i].s; return s < 0 ? `<sub>${t}</sub>` : s > 0 ? `<sup>${t}</sup>` : t; }).join('');
    }
    setFromEditor(root) {
      const runs = [];
      const walk = (n, s) => {
        for (const c of n.childNodes) {
          if (c.nodeType === 3) { if (c.nodeValue) runs.push({ t: c.nodeValue.replace(/ /g, ' '), s }); }
          else if (c.nodeName === 'BR') runs.push({ t: '\n', s: 0 });
          else if (c.nodeName === 'SUB') walk(c, -1);
          else if (c.nodeName === 'SUP') walk(c, 1);
          else if (c.nodeName === 'DIV' || c.nodeName === 'P') { if (runs.length && !/\n$/.test(runs[runs.length - 1].t)) runs.push({ t: '\n', s: 0 }); walk(c, s); }
          else walk(c, s);
        }
      };
      walk(root, 0);
      while (runs.length && /^\n*$/.test(runs[runs.length - 1].t)) runs.pop();
      const merged = [];
      runs.forEach(r => { const l = merged[merged.length - 1]; if (l && l.s === r.s && r.t !== '\n' && l.t !== '\n') l.t += r.t; else merged.push({ ...r }); });
      this.runs = merged;
    }
    toggleScript(s) {                                              // the whole text, when the box is only selected
      const on = this.runs.length && this.runs.every(r => r.s === s);
      this.runs.forEach(r => { r.s = on ? 0 : s; });
      this.render();
    }
    render() {
      const g = this.g; g.textContent = ''; g.setAttribute('transform', `translate(${this.x},${this.y})`);
      const t = h('text', { x: 4, y: 4, 'font-size': TEXT_PX, 'font-family': FONT, fill: this.color, 'xml:space': 'preserve', style: 'white-space:pre' }, g);
      if (this.editing) t.setAttribute('visibility', 'hidden');
      let line = 0, off = 0, first = true;
      const lh = TEXT_PX * 1.25;
      this.runs.forEach(r => {
        r.t.split('\n').forEach((part, k) => {
          if (k > 0) { line++; first = true; off = 0; }
          if (part === '') return;                      // nothing to draw; the next text starts this line (an empty tspan would not move it)
          const sp = h('tspan', { 'font-size': r.s ? TEXT_PX * 0.7 : TEXT_PX }, t);
          const want = r.s < 0 ? TEXT_PX * 0.25 : r.s > 0 ? -TEXT_PX * 0.4 : 0;
          if (first) { sp.setAttribute('x', 4); sp.setAttribute('y', 4 + TEXT_PX * 0.95 + line * lh + want); first = false; }
          else if (want !== off) sp.setAttribute('dy', want - off);
          off = want; sp.textContent = part;
        });
      });
      if (this.selected && !this.editing) { const b = this.localBox(); h('rect', { x: b.x - 1, y: b.y - 1, width: b.w + 2, height: b.h + 2, fill: 'none', stroke: BLUE, 'stroke-width': 1.5, 'stroke-dasharray': '6 4' }, g); }
    }
    localBox() {
      const lines = this.plain().split('\n'), px = TEXT_PX;
      const w = Math.max(...lines.map(l => textWidth(l, px)), 16) + 8;
      return { x: 0, y: 0, w, h: Math.max(lines.length, 1) * px * 1.25 + 8 };
    }
    bbox() { const b = this.localBox(); return { x: this.x, y: this.y, w: b.w, h: b.h }; }
    hit(p) { return inRect(p, this.bbox()); }
    translate(dx, dy) { this.x += dx; this.y += dy; this.g.setAttribute('transform', `translate(${this.x},${this.y})`); }
    meetsRect(r) { return rectsMeet(this.bbox(), r); }
  }

  /* ---- arrows ---- */
  class StraightSnap {                                             // sticky snap to a straight arrow within 20° of horizontal
    constructor() { this.snapped = false; this.armed = true; }
    begin(anchor, point) { const straight = Math.abs(point.y - anchor.y) < 0.5; this.snapped = straight; this.armed = !straight && !this.inZone(anchor, point); }
    inZone(a, p) { const dx = p.x - a.x, dy = p.y - a.y; if (Math.hypot(dx, dy) < 1) return false; const ang = Math.abs(Math.atan2(dy, dx) * 180 / Math.PI); return Math.min(ang, 180 - ang) <= 20; }
    update(a, p) {
      if (!this.inZone(a, p)) { this.snapped = false; this.armed = true; return { x: p.x, y: p.y }; }
      const flat = { x: p.x, y: a.y };
      if (this.snapped) { if (Math.abs(p.y - a.y) > 14) { this.snapped = false; this.armed = false; return { x: p.x, y: p.y }; } return flat; }
      if (this.armed) { this.snapped = true; return flat; }
      return { x: p.x, y: p.y };
    }
  }
  const snapAngle = (a, p, step = 15) => { const dx = p.x - a.x, dy = p.y - a.y, len = Math.hypot(dx, dy), s = step * Math.PI / 180, ang = Math.round(Math.atan2(dy, dx) / s) * s; return { x: a.x + len * Math.cos(ang), y: a.y + len * Math.sin(ang) }; };

  class ArrowItem extends Item {
    constructor(kind, p1, p2) { super(); this.kind = kind; this.p1 = { ...p1 }; this.p2 = { ...p2 }; this.color = '#000'; this.bend = { x: 0, y: 0 }; this.straight = new StraightSnap(); }
    setEnd(p) { this.p2 = { ...p }; this.render(); }
    mid() { return { x: (this.p1.x + this.p2.x) / 2, y: (this.p1.y + this.p2.y) / 2 }; }
    handle() { const m = this.mid(); return { x: m.x + this.bend.x, y: m.y + this.bend.y }; }
    control() { const m = this.mid(); return { x: m.x + this.bend.x * 2, y: m.y + this.bend.y * 2 }; }
    bent() { return Math.abs(this.bend.x) + Math.abs(this.bend.y) > 0; }
    lines() {                                                      // [start, end, headAtEnd, headAtStart, barbSide]
      const len = dist(this.p1, this.p2); if (len < 1) return [];
      if (this.kind === 'forward') return [[this.p1, this.p2, true, false, 0]];
      if (this.kind === 'resonance') return [[this.p1, this.p2, true, true, 0]];
      const nx = -(this.p2.y - this.p1.y) / len * 4, ny = (this.p2.x - this.p1.x) / len * 4;
      return [[{ x: this.p1.x - nx, y: this.p1.y - ny }, { x: this.p2.x - nx, y: this.p2.y - ny }, true, false, 1], [{ x: this.p1.x + nx, y: this.p1.y + ny }, { x: this.p2.x + nx, y: this.p2.y + ny }, false, true, 1]];
    }
    pathD(s, e) { const c = this.control(), off = { x: s.x - this.p1.x, y: s.y - this.p1.y }; return this.bent() ? `M${s.x},${s.y} Q${c.x + off.x},${c.y + off.y} ${e.x},${e.y}` : `M${s.x},${s.y} L${e.x},${e.y}`; }
    samples(s, e) {                                                // points along the drawn line, for hit tests
      const c = this.control(), off = { x: s.x - this.p1.x, y: s.y - this.p1.y }, out = [];
      for (let i = 0; i <= 24; i++) {
        const t = i / 24;
        out.push(this.bent() ? { x: (1 - t) * (1 - t) * s.x + 2 * (1 - t) * t * (c.x + off.x) + t * t * e.x, y: (1 - t) * (1 - t) * s.y + 2 * (1 - t) * t * (c.y + off.y) + t * t * e.y } : { x: s.x + (e.x - s.x) * t, y: s.y + (e.y - s.y) * t });
      }
      return out;
    }
    allSamples() { return this.lines().flatMap(([s, e]) => this.samples(s, e)); }
    head(parent, tip, tail, side) {
      const L = dist(tail, tip); if (L < 1) return;
      const ux = (tip.x - tail.x) / L, uy = (tip.y - tail.y) / L, base = Math.atan2(uy, ux) + Math.PI, size = 11, spread = 0.45;
      const pts = (side === 0 ? [1, -1] : [side]).map(s => ({ x: tip.x + size * Math.cos(base + s * spread), y: tip.y + size * Math.sin(base + s * spread) }));
      if (side === 0) h('polygon', { points: `${tip.x},${tip.y} ${pts[0].x},${pts[0].y} ${pts[1].x},${pts[1].y}`, fill: this.color, stroke: this.color, 'stroke-width': 1, 'stroke-linejoin': 'round' }, parent);
      else h('line', { x1: tip.x, y1: tip.y, x2: pts[0].x, y2: pts[0].y }, parent);
    }
    tailFor(s, e, atEnd) { if (!this.bent()) return atEnd ? s : e; const c = this.control(); return { x: c.x + (s.x - this.p1.x), y: c.y + (s.y - this.p1.y) }; }
    render() {
      const g = this.g; g.textContent = '';
      const grp = h('g', { stroke: this.color, fill: 'none', 'stroke-width': 2, 'stroke-linecap': 'round', 'stroke-linejoin': 'round' }, g);
      for (const [s, e, headEnd, headStart, side] of this.lines()) {
        h('path', { d: this.pathD(s, e) }, grp);
        if (headEnd) this.head(grp, e, this.tailFor(s, e, true), side);
        if (headStart) this.head(grp, s, this.tailFor(s, e, false), side);
      }
      if (this.selected) for (const q of [this.p1, this.p2, this.handle()]) h('circle', { cx: q.x, cy: q.y, r: 5, fill: '#fff', stroke: BLUE, 'stroke-width': 1.5 }, g);
    }
    bbox() { const pts = [this.p1, this.p2, this.control()]; const xs = pts.map(q => q.x), ys = pts.map(q => q.y); return { x: Math.min(...xs) - 12, y: Math.min(...ys) - 12, w: Math.max(...xs) - Math.min(...xs) + 24, h: Math.max(...ys) - Math.min(...ys) + 24 }; }
    hit(p) {
      if (this.selected && dist(p, this.handle()) <= 8) return true;
      for (const [s, e] of this.lines()) { const pts = this.samples(s, e); for (let i = 1; i < pts.length; i++) if (segDist(p, pts[i - 1], pts[i]) <= 7) return true; }
      return false;
    }
    nearEnd(p) { for (const [n, q] of [[1, this.p1], [2, this.p2], [3, this.handle()]]) if (dist(p, q) <= 9) return n; return 0; }
    hoverAt(p) { return this.selected && this.nearEnd(p) ? 'crosshair' : null; }
    snapped(anchor, p, e) { if (e.shiftKey) return snapAngle(anchor, p); if (e.altKey) return { ...p }; return this.straight.update(anchor, p); }
    press(p) {
      if (!this.selected) return null;
      const end = this.nearEnd(p); if (!end) return null;
      if (end === 3) return { move: q => { const m = this.mid(), b = { x: q.x - m.x, y: q.y - m.y }; this.bend = Math.abs(b.x) + Math.abs(b.y) < 8 ? { x: 0, y: 0 } : b; this.render(); }, up: () => {} };
      const anchor = end === 1 ? this.p2 : this.p1; this.straight.begin(anchor, end === 1 ? this.p1 : this.p2);
      return { move: (q, e) => { const r = this.snapped(anchor, q, e); if (end === 1) this.p1 = r; else this.p2 = r; this.render(); }, up: () => {} };
    }
    translate(dx, dy) { this.p1.x += dx; this.p1.y += dy; this.p2.x += dx; this.p2.y += dy; this.render(); }
    meetsRect(r) { return this.allSamples().some(q => inRect(q, r)); }
  }

  /* ---- an electron (one dot) or a lone pair (two) ---- */
  class ElectronItem extends Item {
    constructor(pair, pos, angle = 0) { super(); this.pair = pair; this.x = pos.x; this.y = pos.y; this.angle = angle; this.color = '#000'; }
    ring() { return 11; }
    onRing(p) { const d = Math.hypot(p.x - this.x, p.y - this.y); return this.pair && this.selected && d >= 6 && d <= this.ring() + 6; }
    bbox() { const r = this.ring() + 8; return { x: this.x - r, y: this.y - r, w: 2 * r, h: 2 * r }; }
    hit(p) { const r = this.selected && this.pair ? this.ring() + 6 : 10; return Math.hypot(p.x - this.x, p.y - this.y) <= r; }
    hoverAt(p) { return this.onRing(p) ? 'grab' : null; }
    press(p) { return this.onRing(p) ? { move: q => { this.angle = Math.atan2(q.y - this.y, q.x - this.x); this.render(); }, up: () => {} } : null; }
    render() {
      const g = this.g; g.textContent = ''; g.setAttribute('transform', `translate(${this.x},${this.y})`);
      if (this.pair) { const dx = 3.6 * Math.cos(this.angle), dy = 3.6 * Math.sin(this.angle); h('circle', { cx: -dx, cy: -dy, r: 2.4, fill: this.color }, g); h('circle', { cx: dx, cy: dy, r: 2.4, fill: this.color }, g); }
      else h('circle', { cx: 0, cy: 0, r: 2.6, fill: this.color }, g);
      if (this.selected) h('circle', { cx: 0, cy: 0, r: this.pair ? this.ring() : 9, fill: 'none', stroke: BLUE, 'stroke-width': 1.5, 'stroke-dasharray': '4 3' }, g);
    }
    translate(dx, dy) { this.x += dx; this.y += dy; this.g.setAttribute('transform', `translate(${this.x},${this.y})`); }
    meetsRect(r) { return inRect({ x: this.x, y: this.y }, { x: r.x - 8, y: r.y - 8, w: r.w + 16, h: r.h + 16 }); }
  }

  /* ================= the page ================= */
  const MODES = [
    ['select', 'Select', 'Select and move things (Esc)'],
    ['molecule', 'Molecule box', 'Drag out a box, then import or draw a molecule in it'],
    ['text', 'Text', 'Click to add a text box'],
    ['arrow-forward', '→', 'Reaction arrow (clicks to straight near horizontal; Shift = 15° steps, Alt = no snapping)'],
    ['arrow-equilibrium', '⇌', 'Equilibrium arrow'],
    ['arrow-resonance', '↔', 'Resonance arrow'],
    ['electron', '•', 'Electron: click to place a single electron dot'],
    ['lonepair', '••', 'Lone pair: click to place two dots; drag away from the atom to turn them'],
  ];
  const SWATCHES = ['#000000', '#e03131', '#f08c00', '#2f9e44', '#0c8599', '#1971c2', '#6741d9', '#c2255c', '#862e9c', '#7f5539', '#868e96', '#ffffff'];

  const css = document.createElement('style');
  css.textContent = `
    #page-draw.active { display:flex; flex-direction:column; height:calc(100vh - 44px); }
    .dr-bar { display:flex; align-items:center; gap:6px; padding:8px; flex-wrap:wrap; }
    .dr-bar button, .dr-dlg button { background:#f4f4f6; border:1px solid #d6d6dc; border-radius:6px; padding:6px 12px; color:#222; font:13px system-ui,sans-serif; cursor:pointer; }
    .dr-bar button:hover, .dr-dlg button:hover { background:#ebebef; }
    .dr-bar button.on, .dr-dlg button.on { background:#e0ecff; border-color:#7aa7ff; }
    .dr-gap { width:10px; } .dr-fill { flex:1; }
    .dr-scroll { position:relative; flex:1; overflow:auto; background:#fff; border-top:1px solid #e3e3e3; }
    .dr-scroll svg { display:block; user-select:none; }
    .dr-edit { position:absolute; transform-origin:0 0; font:${TEXT_PX}px ${FONT}; line-height:1.25; padding:4px; min-width:16px; min-height:20px; outline:1.5px dashed ${BLUE}; background:#fff; white-space:pre; z-index:5; }
    .dr-color { position:absolute; z-index:30; background:#fff; border:1px solid #d6d6dc; border-radius:8px; padding:10px; box-shadow:0 6px 24px rgba(0,0,0,.18); width:210px; font:12px system-ui,sans-serif; color:#555; }
    .dr-color .sw { display:grid; grid-template-columns:repeat(6,1fr); gap:5px; margin:6px 0; }
    .dr-color .sw i { height:24px; border-radius:5px; border:1px solid #ccc; cursor:pointer; }
    .dr-color input[type=color] { width:100%; height:30px; border:1px solid #ccc; border-radius:5px; padding:0; background:#fff; }
    .dr-color button { width:100%; margin-top:6px; background:#f4f4f6; border:1px solid #d6d6dc; border-radius:6px; padding:5px; cursor:pointer; }
    .dr-modal { position:fixed; inset:0; background:rgba(0,0,0,.35); z-index:60; display:flex; align-items:center; justify-content:center; }
    .dr-dlg { background:#fff; border-radius:10px; width:min(840px,94vw); max-height:92vh; overflow:auto; padding:14px; box-shadow:0 10px 40px rgba(0,0,0,.3); font:13px system-ui,sans-serif; }
    .dr-dlg .row { display:flex; gap:6px; align-items:center; margin-bottom:8px; position:relative; flex-wrap:wrap; }
    .dr-dlg input[type=text] { flex:1; min-width:200px; padding:7px 9px; border:1px solid #d6d6dc; border-radius:6px; font:13px system-ui,sans-serif; }
    .dr-dlg input[type=text]:focus { border-color:#7aa7ff; outline:none; }
    .dr-dlg .status { color:#666; margin:2px 0 8px; min-height:18px; }
    .dr-dlg .preview { border:1px solid #e3e3e3; border-radius:8px; height:300px; display:flex; align-items:center; justify-content:center; background:#fff; margin-bottom:10px; overflow:hidden; }
    .dr-menu { position:absolute; z-index:40; min-width:150px; background:#fff; border:1px solid #d6d6dc; border-radius:8px; padding:4px; box-shadow:0 6px 24px rgba(0,0,0,.18); }
    .dr-menu button { display:flex; justify-content:space-between; gap:16px; width:100%; border:0; background:none; padding:6px 10px; border-radius:5px; font:13px system-ui,sans-serif; color:#222; cursor:pointer; text-align:left; }
    .dr-menu button:hover:not(:disabled) { background:#eef3ff; }
    .dr-menu button:disabled { color:#aaa; cursor:default; }
    .dr-menu kbd { font:12px system-ui,sans-serif; color:#888; }
    .dr-send { position:fixed; left:50%; bottom:18px; transform:translateX(-50%); z-index:70; background:#1b3a6b; color:#fff; border-radius:10px; padding:10px 14px; font:13px system-ui,sans-serif; display:flex; gap:10px; align-items:center; box-shadow:0 6px 24px rgba(0,0,0,.3); }
    .dr-send button { border:0; border-radius:6px; padding:6px 12px; cursor:pointer; font:13px system-ui,sans-serif; }
  `;
  document.head.appendChild(css);

  page.innerHTML = '';
  const bar = document.createElement('div'); bar.className = 'dr-bar'; page.appendChild(bar);
  const scroll = document.createElement('div'); scroll.className = 'dr-scroll'; page.appendChild(scroll);
  const svg = h('svg', { width: SCENE_W, height: SCENE_H, viewBox: `0 0 ${SCENE_W} ${SCENE_H}` });
  scroll.appendChild(svg);
  h('rect', { x: 0, y: 0, width: SCENE_W, height: SCENE_H, fill: '#fff' }, svg);
  const world = h('g', {}, svg), overlay = h('g', {}, svg);

  const items = [], sel = new Set();
  let mode = 'select', zoom = 1, drag = null, lassoRect = null, lassoEl = null, rubberEl = null;
  const modeButtons = {};
  const mk = (txt, title, fn, cls) => { const b = document.createElement('button'); b.textContent = txt; b.title = title || ''; if (cls) b.className = cls; b.addEventListener('mousedown', e => e.preventDefault()); b.addEventListener('click', fn); bar.appendChild(b); return b; };
  const gap = () => { const d = document.createElement('span'); d.className = 'dr-gap'; bar.appendChild(d); };
  MODES.forEach(([m, txt, tip]) => {
    modeButtons[m] = mk(txt, tip, () => setMode(m));
    if (m === 'text') {
      for (const [t2, tip2, s] of [['x₂', 'Subscript the selected text in a text box', -1], ['x²', 'Superscript the selected text in a text box', 1]]) mk(t2, tip2, () => toggleScript(s));
      for (const [t2, tip2, sym] of [['☀', 'Light: put ☀ in the text box being edited, or a new text box (e.g. over a reaction arrow)', '☀'],
                                     ['hν', 'hν (light, photon energy): put it in the text box being edited, or a new text box', 'hν']]) mk(t2, tip2, () => insertSymbol(sym));
    }
    if (['text', 'arrow-resonance', 'lonepair'].includes(m)) gap();
  });
  gap();
  mk('Delete', 'Delete the selected items (Delete key)', () => deleteSelected());
  const colorBtn = mk('Color ▾', 'Color the selection. With the Select tool, drag a rectangle around atoms, bonds, text or arrows first', () => toggleColorPanel());
  mk('Clear all', '', () => { items.slice().forEach(removeItem); setMode('select'); });
  const fill = document.createElement('span'); fill.className = 'dr-fill'; bar.appendChild(fill);
  mk('Export PNG', 'Save the scheme as a high-resolution PNG', () => exportPng());

  function setMode(m) {
    mode = m;
    Object.entries(modeButtons).forEach(([k, b]) => b.classList.toggle('on', k === m));
    svg.style.cursor = m === 'select' ? 'default' : 'crosshair';
  }
  function addItem(it) { items.push(it); it.scene = true; world.appendChild(it.g); it.render(); return it; }
  function removeItem(it) { const i = items.indexOf(it); if (i >= 0) items.splice(i, 1); sel.delete(it); it.g.remove(); if (editing && editing.item === it) editing = null; }
  function clearSelection() { sel.forEach(i => i.setSelected(false)); sel.clear(); }
  function select(it, on = true) { if (on) { sel.add(it); it.setSelected(true); } else { sel.delete(it); it.setSelected(false); } }
  function deleteSelected() { [...sel].forEach(removeItem); }
  function itemsAt(p) { return items.filter(i => i.hit(p)).reverse(); }
  const pt = e => { const m = svg.getScreenCTM().inverse(), q = new DOMPoint(e.clientX, e.clientY).matrixTransform(m); return { x: q.x, y: q.y }; };

  /* ---- colour */
  let colorPanel = null;
  const recent = () => (store.get('draw.recent') || []).filter(c => /^#[0-9a-f]{6}$/i.test(c)).slice(0, 10);
  function applyColor(color) {                                     // the selection; null = back to black (mechanism-free: arrows and electrons black)
    const list = [...sel];
    list.forEach(it => {
      if (it instanceof MoleculeBox) it.applyColor(color);
      else { it.color = color || '#000'; it.render(); }
    });
    return list.length > 0;
  }
  function pickColor(hex) { if (applyColor(hex)) { store.set('draw.recent', [hex, ...recent().filter(c => c !== hex)].slice(0, 10)); } renderColorPanel(); }
  function renderColorPanel() {
    if (!colorPanel) return;
    colorPanel.innerHTML = 'Colors<div class="sw">' + SWATCHES.map(c => `<i data-c="${c}" style="background:${c}" title="${c}"></i>`).join('') + '</div>' +
      '<input type="color" value="#1971c2" title="Any color"><div style="margin-top:6px">Recent</div><div class="sw">' + (recent().map(c => `<i data-c="${c}" style="background:${c}" title="${c}"></i>`).join('') || '<span style="grid-column:1/-1">none yet</span>') + '</div><button data-reset>Reset to black</button>';
    colorPanel.querySelectorAll('i[data-c]').forEach(i => i.addEventListener('click', () => pickColor(i.dataset.c)));
    colorPanel.querySelector('input[type=color]').addEventListener('change', e => pickColor(e.target.value));
    colorPanel.querySelector('[data-reset]').addEventListener('click', () => applyColor(null));
  }
  function toggleColorPanel() {
    if (colorPanel) { colorPanel.remove(); colorPanel = null; return; }
    colorPanel = document.createElement('div'); colorPanel.className = 'dr-color';
    const r = colorBtn.getBoundingClientRect(), pr = page.getBoundingClientRect();
    colorPanel.style.left = (r.left - pr.left) + 'px'; colorPanel.style.top = (r.bottom - pr.top + 2) + 'px';
    page.style.position = 'relative'; page.appendChild(colorPanel); renderColorPanel();
  }
  document.addEventListener('pointerdown', e => { if (colorPanel && !colorPanel.contains(e.target) && e.target !== colorBtn) { colorPanel.remove(); colorPanel = null; } }, true);

  /* ---- text editing: a contenteditable box over the canvas; the text is stored as runs and drawn in SVG */
  let editing = null;
  function startEditing(tb) {
    stopEditing(); sel.forEach(i => { if (i !== tb) select(i, false); });
    const ed = document.createElement('div'); ed.className = 'dr-edit'; ed.contentEditable = 'true'; ed.spellcheck = false;
    ed.style.left = (tb.x * zoom) + 'px'; ed.style.top = (tb.y * zoom) + 'px'; ed.style.transform = `scale(${zoom})`; ed.style.color = tb.color;
    ed.innerHTML = tb.html(); scroll.appendChild(ed);
    tb.editing = true; tb.render(); editing = { item: tb, el: ed };
    ed.addEventListener('keydown', e => {
      if (e.key === 'Escape') { e.preventDefault(); ed.blur(); }
      else if (e.key === 'Enter' && !e.shiftKey && !e.isComposing) { e.preventDefault(); ed.blur(); clearSelection(); }   // Enter finishes the text and deselects it; Shift+Enter is a new line
      e.stopPropagation();
    });
    ed.addEventListener('blur', () => { if (editing && editing.el === ed) stopEditing(); });
    ed.focus();
    const range = document.createRange(); range.selectNodeContents(ed); range.collapse(false); const s = getSelection(); s.removeAllRanges(); s.addRange(range);
  }
  function stopEditing() {
    if (!editing) return;
    const { item, el } = editing; editing = null;
    item.setFromEditor(el); item.editing = false; el.remove();
    if (!item.plain().trim()) removeItem(item); else item.render();
  }
  /* a symbol (☀, hν) at the caret of the text box being edited, else in a new text box in the middle of what is in view */
  function insertSymbol(sym) {
    if (editing) { editing.el.focus(); document.execCommand('insertText', false, sym); return; }
    const tb = new TextBox(sym);
    tb.x = (scroll.scrollLeft + scroll.clientWidth / 2) / zoom - 12; tb.y = (scroll.scrollTop + scroll.clientHeight / 2) / zoom - 12;
    addItem(tb); clearSelection(); select(tb); setMode('select');
  }
  function toggleScript(s) {
    if (editing) { document.execCommand(s < 0 ? 'subscript' : 'superscript'); return; }
    [...sel].filter(i => i instanceof TextBox).forEach(i => i.toggleScript(s));
  }

  /* ---- the pointer */
  const startsLasso = p => {
    for (const it of itemsAt(p)) { if (it instanceof MoleculeBox && !it.glyphHit(p)) continue; return false; }
    return true;
  };
  function updateHover(p) {
    let cur = null;
    for (const it of sel) { const c = it.hoverAt(p); if (c) { cur = c; break; } }
    svg.style.cursor = cur || (mode === 'select' ? 'default' : 'crosshair');
  }
  svg.addEventListener('pointermove', e => { if (!drag) updateHover(pt(e)); });
  svg.addEventListener('pointerleave', () => { sel.forEach(i => { if (i.leave) i.leave(); }); });

  svg.addEventListener('pointerdown', e => {
    if (e.button !== 0) return;
    stopEditing();
    const p = pt(e), additive = e.ctrlKey || e.shiftKey || e.metaKey;
    if (mode === 'select') {
      for (const it of sel) { const d = it.press(p, e); if (d) { start(d); return; } }       // a resize / turn / handle of a selected item
      if (startsLasso(p)) {
        lassoRect = { a: p, b: p, additive };
        start({ move: q => { lassoRect.b = q; if (manhattan(lassoRect.a, q) > 3) drawLasso(); }, up: () => finishLasso() });
        return;
      }
      const top = itemsAt(p).find(i => !(i instanceof MoleculeBox) || i.glyphHit(p));
      if (additive) select(top, !sel.has(top));
      else if (!sel.has(top)) { clearSelection(); select(top); }
      let last = p;
      start({ move: q => { const dx = q.x - last.x, dy = q.y - last.y; last = q; sel.forEach(i => i.translate(dx, dy)); }, up: () => {} });
      if (top instanceof MoleculeBox) top.clearLasso();                                        // grabbing the whole molecule drops any lasso picks
      return;
    }
    clearSelection();
    if (mode === 'molecule') {
      rubberEl = h('rect', { fill: 'none', stroke: BLUE, 'stroke-width': 1.5, 'stroke-dasharray': '6 4' }, overlay);
      const a = p;
      start({ move: q => { const r = normRect(a, q); ['x', 'y'].forEach(k => rubberEl.setAttribute(k, r[k])); rubberEl.setAttribute('width', r.w); rubberEl.setAttribute('height', r.h); },
              up: q => {
                rubberEl.remove(); rubberEl = null;
                let r = normRect(a, q); if (r.w < 30 || r.h < 30) r = { x: a.x - 110, y: a.y - 80, w: 220, h: 160 };
                const box = new MoleculeBox(r.w, r.h); box.x = r.x; box.y = r.y; addItem(box); select(box); setMode('select'); editMolecule(box);
              } });
    } else if (mode === 'text') {
      const tb = new TextBox(''); tb.x = p.x; tb.y = p.y; addItem(tb); setMode('select'); startEditing(tb); e.preventDefault();
    } else if (mode.startsWith('arrow-')) {
      const ar = addItem(new ArrowItem(mode.slice(6), p, p));
      start({ move: (q, ev) => ar.setEnd(ar.snapped(ar.p1, q, ev)), up: () => { if (manhattan(ar.p1, ar.p2) < 12) ar.setEnd({ x: ar.p1.x + 140, y: ar.p1.y }); } });   // the tool stays: Esc or Select to stop
    } else if (mode === 'electron' || mode === 'lonepair') {
      const el = addItem(new ElectronItem(mode === 'lonepair', p));
      start({ move: q => { if (manhattan(q, p) >= 10) { el.angle = Math.atan2(q.y - p.y, q.x - p.x) + Math.PI / 2; el.render(); } }, up: () => {} });   // the tool stays: Esc or Select to stop
    }
  });
  function start(d) {
    drag = d;
    const mv = e => { if (drag) drag.move(pt(e), e); };
    const up = e => { document.removeEventListener('pointermove', mv); document.removeEventListener('pointerup', up); const d2 = drag; drag = null; if (d2) d2.up(pt(e), e); };
    document.addEventListener('pointermove', mv); document.addEventListener('pointerup', up);
  }
  function drawLasso() {
    if (!lassoEl) lassoEl = h('rect', { fill: 'rgba(74,134,232,.14)', stroke: BLUE, 'stroke-width': 1.5, 'stroke-dasharray': '6 4' }, overlay);
    const r = normRect(lassoRect.a, lassoRect.b); ['x', 'y'].forEach(k => lassoEl.setAttribute(k, r[k])); lassoEl.setAttribute('width', r.w); lassoEl.setAttribute('height', r.h);
  }
  function finishLasso() {
    const { a, b, additive } = lassoRect; lassoRect = null; if (lassoEl) { lassoEl.remove(); lassoEl = null; }
    const r = normRect(a, b);
    if (r.w < 6 && r.h < 6) {                                      // a plain click: whatever is under it
      if (!additive) clearSelection();
      const top = itemsAt(a)[0]; if (top) select(top);
      return;
    }
    if (!additive) clearSelection();
    items.forEach(it => {
      if (it instanceof MoleculeBox) { if (it.lasso(r, additive)) select(it); }
      else if (it.meetsRect(r)) select(it);
    });
    items.forEach(it => { if (it instanceof MoleculeBox && !sel.has(it)) it.clearLasso(); });
  }
  svg.addEventListener('dblclick', e => {
    const top = itemsAt(pt(e))[0];
    if (top instanceof MoleculeBox) editMolecule(top); else if (top instanceof TextBox) { select(top); startEditing(top); }
  });
  scroll.addEventListener('wheel', e => {
    if (!e.ctrlKey) return;
    e.preventDefault();
    zoom = Math.max(0.3, Math.min(3, zoom * (e.deltaY < 0 ? 1.15 : 1 / 1.15)));
    svg.setAttribute('width', SCENE_W * zoom); svg.setAttribute('height', SCENE_H * zoom);
  }, { passive: false });
  /* ---- copy and paste (Ctrl+C / Ctrl+V, or right-click): whole items, kept on this page (molecules, text, arrows, electrons) */
  const TRANSIENT = new Set(['g', 'scene', 'selected', 'hover', 'editing', 'rotating']);
  function cloneItem(it) {
    const c = Object.create(Object.getPrototypeOf(it));
    for (const [k, v] of Object.entries(it)) {
      if (TRANSIENT.has(k)) continue;
      if (v instanceof StraightSnap) { c[k] = new StraightSnap(); continue; }
      try { c[k] = structuredClone(v); } catch (e) { c[k] = v; }            // (nothing here holds a DOM node or a function; just in case)
    }
    c.g = h('g'); c.selected = false; c.scene = null; c.hover = null;
    if ('editing' in it) c.editing = false;
    if ('rotating' in it) c.rotating = false;
    return c;
  }
  let clip = [], pasteN = 0;
  function copySelected() { if (!sel.size) return false; clip = items.filter(i => sel.has(i)).map(cloneItem); pasteN = 0; return true; }
  function paste(at) {                                             // at: where the copied group's top-left goes (else a step down-right of the original)
    if (!clip.length) return;
    const fresh = clip.map(cloneItem);
    let dx, dy;
    if (at) {
      const bs = clip.map(i => i.bbox()), x0 = Math.min(...bs.map(b => b.x)), y0 = Math.min(...bs.map(b => b.y));
      dx = at.x - x0; dy = at.y - y0;
    } else { pasteN++; dx = dy = 24 * pasteN; }
    clearSelection();
    fresh.forEach(it => { addItem(it); it.translate(dx, dy); it.render(); select(it); });
    setMode('select');
  }

  let menu = null;
  const closeMenu = () => { if (menu) { menu.remove(); menu = null; } };
  svg.addEventListener('contextmenu', e => {
    e.preventDefault(); stopEditing(); closeMenu();
    const p = pt(e), top = itemsAt(p).find(i => !(i instanceof MoleculeBox) || i.glyphHit(p));
    if (top && !sel.has(top)) { clearSelection(); select(top); }
    menu = document.createElement('div'); menu.className = 'dr-menu';
    const r = scroll.getBoundingClientRect();
    menu.style.left = (e.clientX - r.left + scroll.scrollLeft) + 'px'; menu.style.top = (e.clientY - r.top + scroll.scrollTop) + 'px';
    for (const [label, key, on, fn] of [['Copy', 'Ctrl+C', sel.size > 0, () => copySelected()], ['Paste', 'Ctrl+V', clip.length > 0, () => paste(p)],
                                        ['Delete', 'Del', sel.size > 0, () => deleteSelected()]]) {
      const b = document.createElement('button'); b.type = 'button'; b.disabled = !on;
      b.innerHTML = `<span>${label}</span><kbd>${key}</kbd>`;
      b.addEventListener('click', () => { closeMenu(); fn(); });
      menu.appendChild(b);
    }
    scroll.appendChild(menu);
  });
  document.addEventListener('pointerdown', e => { if (menu && !menu.contains(e.target)) closeMenu(); }, true);
  scroll.addEventListener('scroll', closeMenu);

  document.addEventListener('keydown', e => {
    if (!page.classList.contains('active') || editing || document.querySelector('.dr-modal')) return;
    if (/^(INPUT|TEXTAREA|SELECT)$/.test((e.target || {}).tagName) || (e.target && e.target.isContentEditable)) return;
    const mod = e.ctrlKey || e.metaKey, k = e.key.toLowerCase();
    if (mod && k === 'c') { if (copySelected()) e.preventDefault(); }
    else if (mod && k === 'x') { if (copySelected()) { e.preventDefault(); deleteSelected(); } }
    else if (mod && k === 'v') { if (clip.length) { e.preventDefault(); paste(); } }
    else if (e.key === 'Delete' || e.key === 'Backspace') { e.preventDefault(); deleteSelected(); }
    else if (e.key === 'Escape') { closeMenu(); setMode('select'); clearSelection(); }
  });

  /* ---- export */
  function exportPng() {
    if (!items.length) return;
    clearSelection();
    let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
    items.forEach(it => { const b = it.bbox(); x0 = Math.min(x0, b.x); y0 = Math.min(y0, b.y); x1 = Math.max(x1, b.x + b.w); y1 = Math.max(y1, b.y + b.h); });
    x0 -= 30; y0 -= 30; x1 += 30; y1 += 30;
    const w = x1 - x0, hgt = y1 - y0, S = 4;
    const out = h('svg', { xmlns: NS, width: Math.round(w * S), height: Math.round(hgt * S), viewBox: `${x0} ${y0} ${w} ${hgt}` });
    h('rect', { x: x0, y: y0, width: w, height: hgt, fill: '#fff' }, out);
    items.forEach(it => out.appendChild(it.exportNode()));
    const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(out)], { type: 'image/svg+xml;charset=utf-8' }));
    const img = new Image();
    img.onload = () => {
      const cv = document.createElement('canvas'); cv.width = Math.round(w * S); cv.height = Math.round(hgt * S);
      const ctx = cv.getContext('2d'); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, cv.width, cv.height); ctx.drawImage(img, 0, 0, cv.width, cv.height);
      URL.revokeObjectURL(url);
      cv.toBlob(b => { const a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = 'scheme.png'; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000); }, 'image/png');
    };
    img.src = url;
  }

  /* ---- the molecule dialog: type a name / CAS / SMILES, or draw it in the Sketcher and send it back */
  const refreshBoxes = () => items.forEach(i => { if (i instanceof MoleculeBox) i.refreshDisplay(); });
  function setOption(name, on) { SHOW[name] = !!on; store.set('draw.show', SHOW); refreshBoxes(); }
  function readMolecule(text) {                                    // → SMILES or '' (the same sources the Sketcher's import uses)
    text = text.trim(); if (!text) return '';
    let g = null;
    try { const r = C.searchMolecule(text); g = r && r.graph; } catch (e) { /* next */ }
    if (!g && window.IupacGen) { const s = window.IupacGen.smilesFor(text); if (s) { try { g = C.parseSmiles(s); } catch (e) { /* next */ } } }
    if (!g) { try { g = C.parseSmiles(text); } catch (e) { /* next */ } }
    if (!g || !g.atoms || !g.atoms.length) { try { g = C.parseIupacName(text); } catch (e) { /* none */ } }
    return g && g.atoms && g.atoms.length ? C.toSmiles(g) : '';
  }
  function editMolecule(box) {
    const modal = document.createElement('div'); modal.className = 'dr-modal';
    modal.innerHTML = `<div class="dr-dlg"><div class="row"><input type="text" placeholder="SMILES, IUPAC name or CAS number  (e.g. CCO, ethanol, 64-17-5)" spellcheck="false" autocomplete="off"><button data-import>Import</button><div class="suggest" hidden></div></div>
      <div class="status">Import a molecule, or draw one in the Sketcher.</div><div class="preview"></div>
      <div class="row"><button data-sketch title="Open this molecule in the Sketcher to draw or edit it; a button there sends it back">Draw in Sketcher</button><button data-clear>Clear</button><span class="dr-gap"></span><span data-opts class="row" style="margin:0"></span><span class="dr-fill"></span><button data-cancel>Cancel</button><button data-done style="background:#e0ecff;border-color:#7aa7ff">Done</button></div></div>`;
    document.body.appendChild(modal);
    const q = s => modal.querySelector(s), input = q('input'), status = q('.status'), previewBox = q('.preview');
    const original = { smiles: box.smiles }, pv = new MoleculeBox(); let current = box.smiles;
    const pvSvg = h('svg', { width: 520, height: 280, viewBox: '0 0 520 280' }); previewBox.appendChild(pvSvg);
    function redraw() {
      pv.setMolecule(current); pvSvg.textContent = '';
      pv.autoFit = true; pv.x = (520 - pv.w) / 2; pv.y = (280 - pv.h) / 2;
      const sc = Math.min(1, 500 / pv.w, 270 / pv.h); pvSvg.setAttribute('viewBox', `${-(520 / sc - 520) / 2} ${-(280 / sc - 280) / 2} ${520 / sc} ${280 / sc}`);
      pv.render(); pvSvg.appendChild(pv.g);
    }
    function use(smiles, note) { current = smiles; box.setMolecule(smiles); redraw(); status.textContent = note || (smiles ? 'Imported.' : 'Empty.'); }
    async function doImport(text) {
      text = (text || input.value).trim(); if (!text) return;
      let smi = readMolecule(text);
      if (!smi && window.PCLib && ((/[a-z]{3}/i.test(text) || /^\d{2,7}-\d\d-\d$/.test(text)) && !/[=#\[\]@]/.test(text))) {
        status.textContent = 'Looking in the offline PubChem library…';
        try { const res = await window.PCLib.resolve(text); if (res && res.ok) smi = res.smiles; } catch (e) { /* none */ }
      }
      if (smi) use(smi, 'Imported "' + text.slice(0, 60) + '".'); else status.textContent = 'Could not read "' + text.slice(0, 60) + '". Try a SMILES, an IUPAC name or a CAS number.';
    }
    if (window.Sketcher && window.Sketcher.attachSuggest) window.Sketcher.attachSuggest(input, q('.suggest'), doImport, null, { search: true });
    else input.addEventListener('keydown', e => { if (e.key === 'Enter') doImport(); });
    q('[data-import]').addEventListener('click', () => doImport());
    q('[data-clear]').addEventListener('click', () => { use('', 'Cleared.'); });
    const opts = q('[data-opts]');
    [['hydrogens', 'Show H', 'Label carbons with their hydrogens (CH3, CH2, CH)'], ['common', 'Common Name', 'Show the common name, if the molecule has one'], ['iupac', 'IUPAC', 'Show the IUPAC name of the molecule'], ['smiles', 'SMILES', 'Show the SMILES code of the molecule'], ['cas', 'CAS', 'Show the CAS number of the molecule (when the offline library has it)']].forEach(([n, t, tip]) => {
      const b = document.createElement('button'); b.textContent = t; b.title = tip; b.classList.toggle('on', SHOW[n]);
      b.addEventListener('click', () => { setOption(n, !SHOW[n]); b.classList.toggle('on', SHOW[n]); redraw(); }); opts.appendChild(b);
    });
    const close = () => { modal.remove(); };
    q('[data-cancel]').addEventListener('click', () => { if (current !== original.smiles && !original.smiles && !current) { /* nothing */ } close(); });
    q('[data-done]').addEventListener('click', close);
    q('[data-sketch]').addEventListener('click', () => {
      if (!window.Sketcher) return;
      try { window.Sketcher.graph = current ? C.parseSmiles(current) : { atoms: [], bonds: [], nextId: 1 }; } catch (e) { /* empty */ }
      close();
      document.querySelector('.tab[data-page="sketcher"]').click();
      const bar2 = document.createElement('div'); bar2.className = 'dr-send';
      bar2.innerHTML = '<span>Draw the molecule, then send it to the box on the Draw tab.</span><button data-send style="background:#fff;color:#1b3a6b">Send to Draw</button><button data-no style="background:transparent;color:#cfe0ff">Cancel</button>';
      document.body.appendChild(bar2);
      bar2.querySelector('[data-no]').addEventListener('click', () => { bar2.remove(); document.querySelector('.tab[data-page="draw"]').click(); });
      bar2.querySelector('[data-send]').addEventListener('click', () => {
        let smi = ''; try { smi = window.Sketcher.smiles(); } catch (e) { /* none */ }
        bar2.remove(); document.querySelector('.tab[data-page="draw"]').click();
        if (smi && items.includes(box)) box.setMolecule(smi);
      });
    });
    modal.addEventListener('pointerdown', e => { if (e.target === modal) close(); });
    redraw(); setTimeout(() => input.focus(), 30);
    if (box.smiles) input.value = box.smiles;
  }

  setMode('select');
  window.Draw = {
    items, MoleculeBox, TextBox, ArrowItem, ElectronItem, setMode, addItem, select, clearSelection, deleteSelected, applyColor, exportPng, SHOW,
    get mode() { return mode; }, get selection() { return [...sel]; }, get zoom() { return zoom; },
  };
})();
