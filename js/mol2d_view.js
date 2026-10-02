/* Mol2DView — the Molecule tab's 2D structure: a publication-style SVG drawing you can hover, click, zoom and pan.
   No chemistry happens here: every atom, bond order, H count, charge, ring and wedge is drawn exactly as given.

   const d = Mol2DView.create(container, { onHover(hit, clientX, clientY), onClick(hit, event), maxBondPx: 56, textbook })
       hit = {type:'atom', index} | {type:'bond', index} | null, with the indices given to setMolecule (record
       indices). onHover fires on every pointer move over an atom or bond (so a card can follow the pointer) and
       once with null on leaving; onClick fires on a click that did not pan (hit null = empty space).
   d.setMolecule({ atoms:[{el, x, y, charge, isotope, hCount, radical}],
                   bonds:[{a, b, order:1|2|3, stereo:null|'wedge'|'hash'|'wavy'|'either', aromatic}],
                   rings:number[][] })                                    null clears. Returns {warnings:[]}.
       x, y in molfile units with y pointing up (as a V2000 block has them); the drawing is scaled so the median
       bond is one standard bond. A wedge or hash starts (narrow end) at atom a, the stereocentre, as in a molfile
       bond line; a bond may carry `narrow: atomIndex` to say otherwise. 'wavy'/'either' on a single bond is drawn
       wavy, on a double bond crossed. Aromatic bonds are drawn with the Kekulé order they carry. H atoms passed in
       are drawn like any other atom; hCount is the number of H written into the atom's label, so pass 0 for an
       atom whose H are atoms. radical = number of unpaired electrons (true = 1), drawn as dots. rings (SSSR) put
       the second line of a ring double bond inside its ring. Clears highlights and annotations.
   d.setElectrons(number[] | null)        lone PAIRS per atom, drawn as pairs of dots around the atom's symbol
   d.setAnnotations(string[] | null)      small italic text beside atoms (R/S, charges, numbers …), in its own
                                          colour, clear of bonds, labels and the sub/superscript places of labels
                                          → {placed, overlapping, hOmitted}: a note with no free spot is left out,
                                          never drawn over anything, and so are notes on H when there are more than
                                          40 heavy atoms; the view says so in its bottom-left corner
   d.highlight('hover'|'select'|'group', {atoms, bonds, color?} | null)   translucent halos under the drawing, one
                                          set per kind; the view sets 'hover' itself unless opts.autoHover === false
   d.setStyle({textbook})  textbook true: all black, heavier lines, regular-weight labels (as printed in textbooks);
                                          false: heteroatoms coloured. Keeps the view, highlights and annotations.
   d.fit()       fit the drawing into the container (bond at most maxBondPx on screen); d.resize() re-measures
   d.exportSVG({highlights, caption, transparent}) → standalone SVG text (cropped, white background unless transparent; select/group halos only if
                                          asked; caption: a string or lines of text set under the drawing; when notes
                                          are left out, a sentence saying so is added to the caption)
   d.exportPNG(scale = 2, {highlights, caption, transparent}) → Promise<Blob>;   d.dispose()
   Extras: d.element (the <svg>); d.clientPosition('atom'|'bond', i) → {x, y}, the client pixel of an atom or a
   bond's middle (for tests).
   Mouse: wheel zooms about the pointer, drag pans, double-click on empty space fits. Touch: drag, pinch.
   Keyboard when focused: + and − zoom, 0 fits, arrow keys pan.

   Drawing conventions (the usual ones, as in IUPAC's 2008 recommendations on structure diagrams — Brecher,
   Pure Appl. Chem. 80, 277–410 — and in RDKit's drawer): H written on the side away from the bonds (OH / HO,
   NH2 / H2N, or above/below), charges as right superscripts (moved to another corner only when a bond points
   straight at that spot), mass numbers as left superscripts, carbon written only when it carries a charge,
   isotope or radical or has no bonds, the second line of a ring double bond inside the ring, a double bond to a
   terminal labelled atom centred, the narrow end of a wedge at the stereocentre. Which side the H goes follows
   RDKit (MolDraw2D, DrawMol getAtomOrientation): east or west unless the open direction is steeper than 70°,
   then above or below (never for one-bond atoms); isolated O, S, Se, Te, Po and halogens write the H first (H2O,
   HCl). Proportions follow RDKit's defaults (label font 0.4 × bond, second line of a double bond 0.15 × bond);
   line width and font follow the app's sketcher (2 px at a 48 px bond, Segoe UI 500). Colours: N blue, O red,
   P orange, S amber, Br brown, I violet as RDKit's 2D palette (S darkened to read on white); F and Cl green as
   in Jmol, darkened; other elements the Jmol colour from MolData darkened to at least 3.5:1 contrast on white.
   Bonds are coloured half and half by their atoms. Pieces of a salt closer than 0.3 bond are set apart in a row.
   Options: maxBondPx (fit cap, default 56), autoHover (default true).                                           */
(() => {
  'use strict';
  const NS = 'http://www.w3.org/2000/svg';

  /* ---- proportions, in drawing units (one standard bond = BL; at zoom 1 a unit is a CSS pixel) ---- */
  const BL = 48;                     // the sketcher's bond length
  const LW = 2;                      // bond line width (sketcher: 2 at 48)
  const FONT = 19;                   // element symbol, 0.4 × BL as RDKit
  const SMALL = 0.7 * FONT;          // H count, charge, mass number
  const NOTE = 0.6 * FONT;           // annotations
  const FAMILY = '"Segoe UI", system-ui, sans-serif';
  const WEIGHT = 500;
  const OFF = 0.15 * BL;             // distance between the lines of a double / triple bond
  const TRIM = 0.8 * OFF;            // an inner ring line is this much shorter at each end
  const WEDGE = 0.085 * BL;          // half width of the wide end of a wedge or hash
  const HASH_GAP = 0.1 * BL;
  const HASH_LW = 1.4;
  const WAVE = 0.07 * BL;            // amplitude of a wavy bond
  const PAD = 0.15 * FONT;           // gap between a label and a bond that stops at it
  const BASE = 0.36 * FONT;          // baseline below the atom centre, so capitals are centred on it
  const SUP = 0.02 * FONT;           // superscript baseline (charge, mass number)
  const SUB = 0.56 * FONT;           // subscript baseline (H count)
  const HLINE = 0.95 * FONT;         // H written above/below the symbol: distance between the two
  const DOT_R = 0.09 * FONT;         // radical dot, and each electron of a lone pair
  const ELECTRON_INK = '#2f4f8d';    // lone-pair dots (coloured drawing); black in the textbook style
  const HIT_ATOM = 0.3 * BL, HIT_BOND = 0.36 * BL;
  const HALO_ATOM = 0.33 * BL, HALO_BOND = 0.3 * BL;
  const INK = '#111';
  /* textbook style: every atom and bond in black, lines a little heavier and labels in regular weight, as structures
     are printed in organic chemistry textbooks; the default style colours heteroatoms (see COLORS) */
  const TEXTBOOK = { ink: '#000', lw: 2.6, hashLw: 1.8, weight: 400 };
  const NOTE_INK = '#9c2f86';        // annotations: a colour no element label uses (6.7:1 on white), in italics
  const NOTE_GAP = 0.2 * NOTE;       // clear space between two annotations
  const LABEL_GAP = 0.3 * FONT;      // clear space between an annotation and any atom label
  const TAN70 = Math.tan(70 * Math.PI / 180);
  const H_FIRST = new Set(['O', 'S', 'Se', 'Te', 'Po', 'F', 'Cl', 'Br', 'I', 'At']);
  const HALO = { group: { color: '#f2b600', opacity: 0.38 }, select: { color: '#2f6fd8', opacity: 0.3 }, hover: { color: '#2f6fd8', opacity: 0.2 } };

  /* RDKit palette hues (N 0000ff, O ff0000, P ff8000, S cccc00, F 33cccc, Cl 00cd00, Br 804c1a, I a11ff0),
     darkened for text on white; F and Cl green as in Jmol */
  const COLORS = { C: INK, H: INK, N: '#2342d8', O: '#e00000', S: '#b38600', P: '#e07000', F: '#3f9a1c',
    Cl: '#15a015', Br: '#9c3a1a', I: '#8b1fc0', B: '#c85f5f', Si: '#9a7048', Se: '#c06a00', Na: '#8a44d0',
    K: '#7a36c0', Li: '#9a58d8', Mg: '#4a8a00', Ca: '#3a8a00', '*': INK };

  const el = (name, attrs = {}, text) => {
    const n = document.createElementNS(NS, name);
    for (const k in attrs) n.setAttribute(k, attrs[k]);
    if (text != null) n.textContent = text;
    return n;
  };
  const r2 = v => Math.round(v * 100) / 100;                 // keep the SVG text short

  /* ---- colours ---- */
  function luminance(hex) {
    const f = v => { v /= 255; return v <= 0.03928 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); };
    const [r, g, b] = [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
    return 0.2126 * f(r) + 0.7152 * f(g) + 0.0722 * f(b);
  }
  /* darken a colour until it reads on white (contrast ratio ≥ 3.5) */
  function legible(hex) {
    if (!/^#[0-9a-f]{6}$/i.test(hex || '')) return INK;
    let [r, g, b] = [1, 3, 5].map(i => parseInt(hex.slice(i, i + 2), 16));
    const hx = () => '#' + [r, g, b].map(v => Math.round(v).toString(16).padStart(2, '0')).join('');
    for (let k = 0; k < 40 && 1.05 / (luminance(hx()) + 0.05) < 3.5; k++) { r *= 0.9; g *= 0.9; b *= 0.9; }
    return hx();
  }
  function elementColor(sym, given) {
    if (given) return legible(given);
    if (COLORS[sym]) return COLORS[sym];
    let c = INK;
    try {
      const e = window.MolData && window.MolData.element && window.MolData.element(sym);
      if (e && e.color) c = legible(e.color);
    } catch (err) { /* no table: ink */ }
    COLORS[sym] = c;
    return c;
  }

  /* ---- text widths: measured with a canvas in the same font, so labels can be laid out before the SVG is
     visible (getBBox needs a rendered element) ---- */
  let ctx2d = null;
  const widths = new Map();
  function textWidth(text, size) {
    let w = widths.get(text);
    if (w === undefined) {
      w = 0;
      try {
        if (!ctx2d) ctx2d = document.createElement('canvas').getContext('2d');
        ctx2d.font = `${WEIGHT} 100px ${FAMILY}`;
        w = ctx2d.measureText(text).width / 100;
      } catch (err) { w = 0; }
      if (!(w > 0)) w = 0.6 * text.length;                 // no canvas: an average glyph width
      widths.set(text, w);
    }
    return w * size;
  }

  /* ---- geometry ---- */
  const angDiff = (a, b) => { let d = Math.abs(a - b) % (2 * Math.PI); return d > Math.PI ? 2 * Math.PI - d : d; };
  const boxOf = (x0, y0, x1, y1) => ({ x0, y0, x1, y1 });
  const grow = (b, g) => ({ x0: b.x0 - g, y0: b.y0 - g, x1: b.x1 + g, y1: b.y1 + g });
  const shift = (b, dx, dy) => ({ x0: b.x0 + dx, y0: b.y0 + dy, x1: b.x1 + dx, y1: b.y1 + dy });
  const boxesMeet = (a, b) => a.x0 < b.x1 && b.x0 < a.x1 && a.y0 < b.y1 && b.y0 < a.y1;
  /* the parameter interval where the line p + t (q − p) is inside box b, or null */
  function slab(x1, y1, x2, y2, b) {
    let t0 = -Infinity, t1 = Infinity;
    for (const [p, d, lo, hi] of [[x1, x2 - x1, b.x0, b.x1], [y1, y2 - y1, b.y0, b.y1]]) {
      if (Math.abs(d) < 1e-12) { if (p < lo || p > hi) return null; continue; }
      const ta = (lo - p) / d, tb = (hi - p) / d;
      t0 = Math.max(t0, Math.min(ta, tb)); t1 = Math.min(t1, Math.max(ta, tb));
    }
    return t0 <= t1 ? [t0, t1] : null;
  }
  const segMeetsBox = (x1, y1, x2, y2, b) => { const s = slab(x1, y1, x2, y2, b); return !!s && s[0] <= 1 && s[1] >= 0; };
  /* how far along p→q the segment first gets clear of a group of boxes around p (0 when p is outside them) */
  function clearOf(x1, y1, x2, y2, boxes) {
    let t = 0;
    for (let pass = 0; pass < 6; pass++) {
      let moved = false;
      for (const b of boxes) {
        const s = slab(x1, y1, x2, y2, b);
        if (s && s[0] <= t + 1e-9 && s[1] > t + 1e-9) { t = s[1]; moved = true; }
      }
      if (!moved) break;
    }
    return t;
  }
  /* where the line p + s·d meets the segment e→m (s, t parameters), or null */
  function meet(px, py, dx, dy, ex, ey, wx, wy) {
    const den = dx * wy - dy * wx;
    if (Math.abs(den) < 1e-9) return null;
    const qx = ex - px, qy = ey - py;
    return { s: (qx * wy - qy * wx) / den, t: (qx * dy - qy * dx) / den };
  }

  /* =============================== the depiction =============================== */
  /* Builds everything that does not depend on the view: positions, labels, bond lines, bounding box. */
  function depict(mol, style) {
    const mono = !!(style && style.mono);
    const warnings = [];
    const src = mol.atoms || [];
    const atoms = src.map((a, i) => {
      const sym = String(a.el || a.element || '*');
      const rad = a.radical === true ? 1 : Math.max(0, Math.round(+a.radical || 0));
      return { i, el: sym, x: +a.x || 0, y: +a.y || 0, charge: Math.round(+a.charge || 0), isotope: Math.round(+a.isotope || 0),
        hCount: Math.max(0, Math.round(+a.hCount || 0)), radical: rad, color: mono ? TEXTBOOK.ink : elementColor(sym, a.color), nb: [] };
    });
    const N = atoms.length;
    const bonds = [];
    (mol.bonds || []).forEach((b, k) => {
      const a = b && b.a, c = b && b.b;
      if (!Number.isInteger(a) || !Number.isInteger(c) || a < 0 || c < 0 || a >= N || c >= N || a === c) {
        warnings.push(`bond ${k} skipped: atoms ${a}, ${c}`); return;
      }
      let stereo = b.stereo ? String(b.stereo).toLowerCase() : null;
      if (stereo === 'dash') stereo = 'hash';
      const narrow = b.narrow === c ? c : a;
      const order = b.order === 0 ? 0 : (+b.order || 1);
      bonds.push({ k, a, b: c, order, stereo, aromatic: !!b.aromatic, narrow });
      atoms[a].nb.push({ j: c, k }); atoms[c].nb.push({ j: a, k });
    });

    /* scale: the median bond becomes one standard bond; molfile y points up, SVG y down */
    const lens = bonds.map(b => Math.hypot(atoms[b.a].x - atoms[b.b].x, atoms[b.a].y - atoms[b.b].y)).filter(l => l > 1e-6).sort((p, q) => p - q);
    const med = lens.length ? lens[lens.length >> 1] : 1.5;
    const k = BL / med;
    for (const a of atoms) { a.X = a.x * k; a.Y = -a.y * k; }

    const labelled = a => a.el !== 'C' || a.charge !== 0 || a.isotope > 0 || a.radical > 0 || a.nb.length === 0;
    for (const a of atoms) { a.labelled = labelled(a); a.pieces = []; a.dots = []; a.boxes = []; }
    for (const a of atoms) if (a.labelled) layoutLabel(a, atoms);
    separateFragments(atoms, bonds);
    for (const a of atoms) finishLabel(a);

    /* rings: which bonds are in which ring, and ring centres, for placing the second line of a double bond */
    const ringsOf = new Map();
    const bkey = (p, q) => (p < q ? p + ',' + q : q + ',' + p);
    const rings = [];
    for (const r of (Array.isArray(mol.rings) ? mol.rings : [])) {
      if (!Array.isArray(r) || r.length < 3 || r.some(i => !Number.isInteger(i) || i < 0 || i >= N)) continue;
      const ring = { atoms: r, cx: 0, cy: 0, doubles: 0 };
      for (const i of r) { ring.cx += atoms[i].X / r.length; ring.cy += atoms[i].Y / r.length; }
      rings.push(ring);
      r.forEach((i, t) => { const key = bkey(i, r[(t + 1) % r.length]); if (!ringsOf.has(key)) ringsOf.set(key, []); ringsOf.get(key).push(ring); });
    }
    for (const b of bonds) {
      b.rings = ringsOf.get(bkey(b.a, b.b)) || [];
      if (b.order === 2) for (const ring of b.rings) ring.doubles++;
    }
    for (const b of bonds) drawBond(b, atoms);

    /* the drawing's bounding box */
    const bbox = { x0: Infinity, y0: Infinity, x1: -Infinity, y1: -Infinity };
    const take = (x, y, r = 0) => { bbox.x0 = Math.min(bbox.x0, x - r); bbox.y0 = Math.min(bbox.y0, y - r); bbox.x1 = Math.max(bbox.x1, x + r); bbox.y1 = Math.max(bbox.y1, y + r); };
    for (const a of atoms) {
      take(a.X, a.Y, LW);
      for (const p of a.pieces) { take(p.box.x0, p.box.y0); take(p.box.x1, p.box.y1); }
      for (const d of a.dots) take(d.x, d.y, DOT_R);
    }
    for (const b of bonds) {
      for (const s of b.segs) { take(s.x1, s.y1, LW); take(s.x2, s.y2, LW); }
      for (const p of b.polys) for (const q of p.pts) take(q[0], q[1], 0.5);
      for (const q of b.extentPts || []) take(q[0], q[1], 0);
    }
    if (!N) Object.assign(bbox, { x0: 0, y0: 0, x1: 0, y1: 0 });
    return { atoms, bonds, byK: new Map(bonds.map(b => [b.k, b])), bbox, warnings, notes: [] };
  }

  /* ---- labels ---- */
  /* which side the H (and dots) go: E W N S */
  function orientation(a, atoms) {
    if (!a.nb.length) return H_FIRST.has(a.el) ? 'W' : 'E';
    let sx = 0, sy = 0;
    const dirs = [];
    for (const { j } of a.nb) {
      const dx = atoms[j].X - a.X, dy = atoms[j].Y - a.Y, L = Math.hypot(dx, dy) || 1;
      sx -= dx / L; sy -= dy / L; dirs.push(Math.atan2(dy, dx));
    }
    let o;
    if (Math.hypot(sx, sy) < 0.1) o = 'E';                              // balanced: no preferred side
    else if (Math.abs(sy) > TAN70 * Math.abs(sx)) o = a.nb.length > 1 ? (sy < 0 ? 'N' : 'S') : 'E';
    else o = sx > 0 ? 'E' : 'W';
    /* the H must not sit on a bond (crowded atoms): then take the freest side */
    if (a.hCount) {
      const ANG = { E: 0, W: Math.PI, S: Math.PI / 2, N: -Math.PI / 2 };
      const clear = c => Math.min(...dirs.map(d => angDiff(d, ANG[c])));
      if (clear(o) < 0.6) {
        let best = o, bc = clear(o);
        for (const c of ['E', 'W', 'S', 'N']) {
          if (a.nb.length === 1 && (c === 'N' || c === 'S')) continue;
          if (clear(c) > bc + 1e-9) { best = c; bc = clear(c); }
        }
        o = best;
      }
    }
    return o;
  }

  const chargeText = c => (Math.abs(c) > 1 ? Math.abs(c) : '') + (c > 0 ? '+' : '−');

  /* Lays out one label around (0, 0): pieces {text, size, x (left), base, w, box}; the atom's symbol is centred
     on the atom. Positions are relative until finishLabel adds the atom position. */
  function layoutLabel(a, atoms) {
    const o = a.orient = orientation(a, atoms);
    const put = (text, size, x, base) => {
      const w = textWidth(text, size);
      const p = { text, size, x, base, w, box: boxOf(x, base - 0.72 * size, x + w, base) };
      a.pieces.push(p);
      return p;
    };
    const sym = put(a.el, FONT, -textWidth(a.el, FONT) / 2, BASE);
    let left = sym.x, right = sym.x + sym.w;
    if (a.isotope) { const t = String(a.isotope), w = textWidth(t, SMALL); put(t, SMALL, left - w, SUP); left -= w; }
    const hN = a.hCount, count = hN > 1 ? String(hN) : '';
    if (hN && o === 'E') {
      const h = put('H', FONT, right, BASE); right += h.w;
      if (count) { const n = put(count, SMALL, right, SUB); right += n.w; }
    } else if (hN && o === 'W') {
      if (count) { const w = textWidth(count, SMALL); put(count, SMALL, left - w, SUB); left -= w; }
      const w = textWidth('H', FONT); put('H', FONT, left - w, BASE); left -= w;
    } else if (hN) {
      const dy = o === 'S' ? HLINE : -HLINE, w = textWidth('H', FONT);
      put('H', FONT, -w / 2, BASE + dy);
      if (count) put(count, SMALL, w / 2, SUB + dy);
    }

    /* the bonds out of this atom, relative to it, for keeping the charge and dots off them */
    const segs = a.nb.map(({ j }) => [0, 0, atoms[j].X - a.X, atoms[j].Y - a.Y]);
    const free = box => !segs.some(s => segMeetsBox(s[0], s[1], s[2], s[3], grow(box, LW))) &&
      !a.pieces.some(p => boxesMeet(p.box, box));

    /* the charge: upper right of the label (IUPAC), unless a bond points straight at that spot (within 30°);
       then upper left, lower right, above, below, lower left. A bond passing near it is clipped around it. */
    if (a.charge) {
      const t = chargeText(a.charge), w = textWidth(t, SMALL), g = 0.04 * FONT;
      const symR = sym.x + sym.w, symL = sym.x;
      const dirs = segs.map(q => Math.atan2(q[3], q[2]));
      const cands = [[right + g, SUP], [left - g - w, SUP], [symR + g, SUB + 0.1 * FONT], [-w / 2, -0.42 * FONT],
        [-w / 2, 0.94 * FONT], [symL - g - w, SUB + 0.1 * FONT]];
      const boxAt = ([x, base]) => boxOf(x, base - 0.72 * SMALL, x + w, base);
      const ok = box => {
        const ang = Math.atan2((box.y0 + box.y1) / 2, (box.x0 + box.x1) / 2);
        return !dirs.some(d => angDiff(d, ang) < Math.PI / 6) && !a.pieces.some(p => boxesMeet(p.box, box));
      };
      const pick = cands.find(c => ok(boxAt(c))) || cands[0];
      put(t, SMALL, pick[0], pick[1]);
    }
    if (a.radical) {
      const n = Math.min(a.radical, 4), gap = 0.14 * FONT, step = 0.36 * FONT;
      const ext = () => a.pieces.reduce((e, p) => ({ x0: Math.min(e.x0, p.box.x0), x1: Math.max(e.x1, p.box.x1), y0: Math.min(e.y0, p.box.y0), y1: Math.max(e.y1, p.box.y1) }), { x0: Infinity, x1: -Infinity, y0: Infinity, y1: -Infinity });
      const e = ext();
      const along = (cx, cy, vertical) => Array.from({ length: n }, (_, t) => {
        const off = (t - (n - 1) / 2) * step;
        return vertical ? { x: cx, y: cy + off } : { x: cx + off, y: cy };
      });
      const cands = {
        E: along(e.x1 + gap + DOT_R, 0, true), W: along(e.x0 - gap - DOT_R, 0, true),
        N: along(0, sym.box.y0 - gap - DOT_R, false), S: along(0, sym.box.y1 + gap + DOT_R, false),
      };
      const order = [o, 'E', 'W', 'N', 'S'];
      const dotBox = ds => boxOf(Math.min(...ds.map(d => d.x)) - DOT_R, Math.min(...ds.map(d => d.y)) - DOT_R, Math.max(...ds.map(d => d.x)) + DOT_R, Math.max(...ds.map(d => d.y)) + DOT_R);
      const side = order.find(s => free(dotBox(cands[s]))) || o;
      a.dots = cands[side];
    }
  }

  /* relative → absolute; the boxes a bond stops at */
  function finishLabel(a) {
    for (const p of a.pieces) { p.x += a.X; p.base += a.Y; p.box = shift(p.box, a.X, a.Y); }
    for (const d of a.dots) { d.x += a.X; d.y += a.Y; }
    a.boxes = a.pieces.map(p => grow(p.box, PAD)).concat(a.dots.map(d => boxOf(d.x - DOT_R - PAD, d.y - DOT_R - PAD, d.x + DOT_R + PAD, d.y + DOT_R + PAD)));
    const all = a.pieces.map(p => p.box).concat(a.dots.map(d => boxOf(d.x - DOT_R, d.y - DOT_R, d.x + DOT_R, d.y + DOT_R)));
    a.extent = all.length ? all.reduce((e, b) => boxOf(Math.min(e.x0, b.x0), Math.min(e.y0, b.y0), Math.max(e.x1, b.x1), Math.max(e.y1, b.y1))) : null;
  }

  /* The pieces of a salt or mixture, when any two come closer than 0.3 bond (e.g. two ions the layout put 1 Å
     apart), are laid out again in a row, left to right, each moved whole: a translation changes nothing about
     the structure or its wedges. Left alone when nothing is that close. */
  function separateFragments(atoms, bonds) {
    const N = atoms.length;
    if (N < 2) return;
    const root = atoms.map((_, i) => i);
    const find = i => (root[i] === i ? i : (root[i] = find(root[i])));
    for (const b of bonds) root[find(b.a)] = find(b.b);
    const groups = new Map();
    atoms.forEach((a, i) => { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(a); });
    if (groups.size < 2) return;
    const frags = [...groups.values()].map(list => {
      const bx = boxOf(Infinity, Infinity, -Infinity, -Infinity);
      for (const a of list) {
        const bs = a.labelled ? a.pieces.map(p => shift(p.box, a.X, a.Y)) : [boxOf(a.X - LW, a.Y - LW, a.X + LW, a.Y + LW)];
        for (const b of bs) { bx.x0 = Math.min(bx.x0, b.x0); bx.y0 = Math.min(bx.y0, b.y0); bx.x1 = Math.max(bx.x1, b.x1); bx.y1 = Math.max(bx.y1, b.y1); }
      }
      return { list, box: bx };
    });
    const gap = 0.6 * BL, near = 0.3 * BL;                             // move them only when closer than `near`
    let clash = false;
    for (let p = 0; p < frags.length && !clash; p++) for (let q = p + 1; q < frags.length; q++) if (boxesMeet(grow(frags[p].box, near / 2), grow(frags[q].box, near / 2))) { clash = true; break; }
    if (!clash) return;
    frags.sort((p, q) => (p.box.x0 + p.box.x1) - (q.box.x0 + q.box.x1));
    const cy = (frags[0].box.y0 + frags[0].box.y1) / 2;
    let x = frags[0].box.x1 + gap;
    for (let t = 1; t < frags.length; t++) {
      const f = frags[t], dx = x - f.box.x0, dy = cy - (f.box.y0 + f.box.y1) / 2;
      for (const a of f.list) { a.X += dx; a.Y += dy; }
      x += f.box.x1 - f.box.x0 + gap;
    }
  }

  /* ---- bonds ---- */
  /* Where the second line of a double bond goes: 'ring' (inside the ring), 'side' (towards the substituents)
     or 'centre' (both lines either side of the atom–atom line). sign: which side for ring/side. */
  function doubleStyle(b, atoms, P, Q, nx, ny) {
    const side = o => Math.sign((atoms[o].X - P.X) * nx + (atoms[o].Y - P.Y) * ny);
    if (b.rings.length) {
      const ring = b.rings.slice().sort((r, s) => (s.doubles - r.doubles) || ((r.atoms.length === 6 ? 0 : 1) - (s.atoms.length === 6 ? 0 : 1)) || (r.atoms.length - s.atoms.length))[0];
      const sg = Math.sign((ring.cx - P.X) * nx + (ring.cy - P.Y) * ny);
      if (sg) return { mode: 'ring', sign: sg };
    }
    const others = P.nb.filter(n => n.j !== Q.i).map(n => n.j).concat(Q.nb.filter(n => n.j !== P.i).map(n => n.j));
    if ((P.nb.length === 1 && P.labelled) || (Q.nb.length === 1 && Q.labelled) || !others.length) return { mode: 'centre' };
    const sum = others.reduce((s, o) => s + side(o), 0);
    if (sum) return { mode: 'side', sign: Math.sign(sum) };
    if (P.nb.length === 1 || Q.nb.length === 1) return { mode: 'centre' };
    const first = others.map(side).find(s => s !== 0);
    return first ? { mode: 'side', sign: first } : { mode: 'centre' };
  }

  /* The drawn pieces of a bond: segs (lines, each {x1,y1,x2,y2,c1,c2,dash?}), polys (wedges), paths (wavy),
     and the hit line. Line ends stop at labels; colours split half and half. */
  function drawBond(b, atoms) {
    const P = atoms[b.a], Q = atoms[b.b];
    const L = Math.hypot(Q.X - P.X, Q.Y - P.Y) || 1e-9;
    const ux = (Q.X - P.X) / L, uy = (Q.Y - P.Y) / L, nx = -uy, ny = ux;
    b.segs = []; b.polys = []; b.paths = [];
    b.hit = { x1: P.X, y1: P.Y, x2: Q.X, y2: Q.Y };

    /* a line from p to q (both absolute), clipped at labelled ends; returns the clipped ends or null */
    const clip = (x1, y1, x2, y2) => {
      let t1 = P.labelled ? clearOf(x1, y1, x2, y2, P.boxes) : 0;
      let t2 = Q.labelled ? clearOf(x2, y2, x1, y1, Q.boxes) : 0;
      t1 = Math.min(t1, 0.49); t2 = Math.min(t2, 0.49);
      if (1 - t1 - t2 < 0.02) return null;
      const dx = x2 - x1, dy = y2 - y1;
      return { x1: x1 + dx * t1, y1: y1 + dy * t1, x2: x2 - dx * t2, y2: y2 - dy * t2 };
    };
    const line = (x1, y1, x2, y2, extra) => {
      const c = clip(x1, y1, x2, y2);
      if (c) b.segs.push(Object.assign(c, { c1: P.color, c2: Q.color }, extra || {}));
    };
    const off = (o, t1 = 0, t2 = 0) => [P.X + nx * o + ux * t1, P.Y + ny * o + uy * t1, Q.X + nx * o - ux * t2, Q.Y + ny * o - uy * t2];

    const st = b.stereo;
    if (b.order === 1 && (st === 'wedge' || st === 'hash')) {
      const Nw = b.narrow === b.b ? Q : P, Ww = Nw === P ? Q : P;
      const c = clip(P.X, P.Y, Q.X, Q.Y);
      if (!c) return;
      const [sx, sy, ex, ey] = Nw === P ? [c.x1, c.y1, c.x2, c.y2] : [c.x2, c.y2, c.x1, c.y1];
      const len = Math.hypot(ex - sx, ey - sy);
      const vx = nx, vy = ny;
      if (st === 'wedge') {
        const mx = (sx + ex) / 2, my = (sy + ey) / 2, hw = WEDGE / 2;
        if (Nw.color === Ww.color) b.polys.push({ pts: [[sx, sy], [ex + vx * WEDGE, ey + vy * WEDGE], [ex - vx * WEDGE, ey - vy * WEDGE]], color: Nw.color });
        else {
          b.polys.push({ pts: [[sx, sy], [mx + vx * hw, my + vy * hw], [mx - vx * hw, my - vy * hw]], color: Nw.color });
          b.polys.push({ pts: [[mx + vx * hw, my + vy * hw], [ex + vx * WEDGE, ey + vy * WEDGE], [ex - vx * WEDGE, ey - vy * WEDGE], [mx - vx * hw, my - vy * hw]], color: Ww.color });
        }
      } else {
        const steps = Math.max(3, Math.round(len / HASH_GAP));
        for (let s = 1; s <= steps; s++) {
          const t = s / steps, w = 0.8 + (WEDGE - 0.8) * t, px = sx + (ex - sx) * t, py = sy + (ey - sy) * t;
          b.segs.push({ x1: px + vx * w, y1: py + vy * w, x2: px - vx * w, y2: py - vy * w, c1: t <= 0.5 ? Nw.color : Ww.color, c2: t <= 0.5 ? Nw.color : Ww.color, width: HASH_LW, butt: true });
        }
      }
      return;
    }
    if (b.order === 1 && (st === 'wavy' || st === 'either')) {
      const c = clip(P.X, P.Y, Q.X, Q.Y);
      if (!c) return;
      const len = Math.hypot(c.x2 - c.x1, c.y2 - c.y1), n = Math.max(2, Math.round(len / (0.14 * BL)));
      const half = [[], []];
      for (let s = 0; s < n; s++) {
        const t0 = s / n, t1 = (s + 1) / n, tm = (t0 + t1) / 2, sg = s % 2 ? -1 : 1;
        const pt = t => [c.x1 + (c.x2 - c.x1) * t, c.y1 + (c.y2 - c.y1) * t];
        const [ax, ay] = pt(t0), [bx, by] = pt(t1), [mx, my] = pt(tm);
        half[s < n / 2 ? 0 : 1].push({ ax, ay, bx, by, cx: mx + nx * 2 * WAVE * sg, cy: my + ny * 2 * WAVE * sg });
      }
      const d = list => list.map((q, t) => (t ? '' : `M${r2(q.ax)} ${r2(q.ay)}`) + `Q${r2(q.cx)} ${r2(q.cy)} ${r2(q.bx)} ${r2(q.by)}`).join('');
      if (P.color === Q.color) b.paths.push({ d: d(half[0].concat(half[1])), color: P.color });
      else half.forEach((h, t) => h.length && b.paths.push({ d: d(h), color: t ? Q.color : P.color }));
      b.extentPts = half[0].concat(half[1]).map(q => [q.cx, q.cy]);
      return;
    }
    if (b.order === 2 && (st === 'either' || st === 'wavy')) {           // crossed: E/Z unknown
      const h = OFF / 2;
      line(P.X + nx * h, P.Y + ny * h, Q.X - nx * h, Q.Y - ny * h);
      line(P.X - nx * h, P.Y - ny * h, Q.X + nx * h, Q.Y + ny * h);
      return;
    }
    if (b.order === 2 || b.order === 1.5 || b.order === 4) {
      const aromaticOnly = b.order !== 2;                                // no Kekulé order given: dashed inner line
      const ds = doubleStyle(b, atoms, P, Q, nx, ny);
      if (ds.mode === 'centre' && !aromaticOnly) {
        const h = OFF / 2;
        for (const sg of [1, -1]) {
          let [x1, y1, x2, y2] = off(sg * h);
          if (!P.labelled) [x1, y1] = miter(P, Q.i, x1, y1, ux, uy, atoms);
          if (!Q.labelled) [x2, y2] = miter(Q, P.i, x2, y2, -ux, -uy, atoms);
          line(x1, y1, x2, y2);
        }
        return;
      }
      const sg = ds.sign || 1;
      line(P.X, P.Y, Q.X, Q.Y);
      const tP = P.nb.length > 1 ? TRIM : 0, tQ = Q.nb.length > 1 ? TRIM : 0;
      const [x1, y1, x2, y2] = off(sg * OFF, tP, tQ);
      line(x1, y1, x2, y2, aromaticOnly ? { dash: true } : null);
      b.hit = { x1: P.X + nx * sg * OFF / 2, y1: P.Y + ny * sg * OFF / 2, x2: Q.X + nx * sg * OFF / 2, y2: Q.Y + ny * sg * OFF / 2 };
      return;
    }
    if (b.order === 3) {
      for (const o of [-OFF, 0, OFF]) { const [x1, y1, x2, y2] = off(o); line(x1, y1, x2, y2); }
      return;
    }
    line(P.X, P.Y, Q.X, Q.Y, b.order === 0 || b.order !== 1 ? { dash: true } : null);
  }

  /* A centred double-bond line ending at an unlabelled atom E is extended or shortened to meet E's other bonds,
     so the join is clean (the "=O" of acetone meets both C–C lines). (x, y): the line's end at E; (dx, dy): the
     line's direction, pointing towards E. */
  function miter(E, otherEnd, x, y, dx, dy, atoms) {
    let best = null;
    for (const { j } of E.nb) {
      if (j === otherEnd) continue;
      const M = atoms[j], len = Math.hypot(M.X - E.X, M.Y - E.Y) || 1;
      const m = meet(x, y, dx, dy, E.X, E.Y, (M.X - E.X) / len, (M.Y - E.Y) / len);
      if (m && m.t > 0 && m.t < 0.5 * len && Math.abs(m.s) < 2 * OFF && (!best || Math.abs(m.s) < Math.abs(best.s))) best = m;
    }
    return best ? [x + dx * best.s, y + dy * best.s] : [x, y];
  }

  /* ---- annotations. Each goes to the free spot nearest its atom, preferring the side away from the atom's
     bonds. A small number set beside an atom label reads as part of it (a charge, an H count, a mass number: "OH"
     with a "2" at its top right reads OH²), so notes are drawn in italics in their own colour, kept a gap away
     from every label, never put level with a label's line of text just left or right of it (a sub- or
     superscript's place), and always clearly nearer its own atom than any other. A note with no such free spot
     (clear of bonds, labels and other notes too) is left out, never drawn anyway; so are the notes on H atoms
     when there are more than 40 heavy atoms (they crowd the heavy atoms' notes out). The hover card has every
     value. ---- */
  function isPlainH(a, atoms) { return a.el === 'H' && a.nb.length === 1 && atoms[a.nb[0].j].el !== 'H'; }
  /* every drawn bond line as a segment [x1, y1, x2, y2] (wedges by their outline, wavy bonds by their axis) */
  function bondSegs(scene) {
    const { atoms, bonds } = scene, segs = [];
    for (const b of bonds) {
      for (const s of b.segs) segs.push([s.x1, s.y1, s.x2, s.y2]);
      for (const p of b.polys) for (let t = 0; t < p.pts.length; t++) { const u = p.pts[t], v = p.pts[(t + 1) % p.pts.length]; segs.push([u[0], u[1], v[0], v[1]]); }
      if (b.paths.length) segs.push([atoms[b.a].X, atoms[b.a].Y, atoms[b.b].X, atoms[b.b].Y]);
    }
    return segs;
  }

  /* ---- electrons: each lone pair as two dots beside the atom's symbol, as in a Lewis structure (the bonds are the
     bonding pairs; a radical's unpaired electron is already a dot in the label). The pairs go on the sides of the
     symbol — above, below, left, right, or the diagonals — spread as far as they can from each other and from the
     atom's bonds and its H (a repulsion score, as VSEPR spreads electron pairs), with the four sides preferred, as
     textbooks draw them. A pair never goes over a bond, a label or another pair; one with no free side is left out
     and counted, like a note. Only atoms with a drawn symbol get pairs (a carbon with a lone pair is a carbanion and
     always carries its charge, so it is labelled). ---- */
  function placeElectrons(scene, counts) {
    const { atoms } = scene, marks = [], info = { pairs: 0, omitted: 0 };
    if (!counts) return { marks, info };
    const segs = bondSegs(scene);
    const SEP = 0.36 * FONT, GAP = 0.2 * FONT, R = DOT_R;
    const blocks = [];
    for (const a of atoms) {
      a.pieces.forEach((p, k) => blocks.push({ i: a.i, k, box: grow(p.box, 0.04 * FONT) }));
      for (const d of a.dots) blocks.push({ i: a.i, k: -1, box: boxOf(d.x - R, d.y - R, d.x + R, d.y + R) });
    }
    const taken = [];                                   // boxes of the pairs placed so far
    const DIRS = [];
    for (let s = 0; s < 8; s++) DIRS.push({ ang: s * Math.PI / 4, card: s % 2 === 0 });
    const energy = (u, v) => 1 / Math.max(1e-3, 1 - Math.cos(angDiff(u, v)));   // pairs push apart, like VSEPR
    counts.forEach((n, i) => {
      n = Math.round(+n || 0);
      if (n <= 0 || i >= atoms.length) return;
      const a = atoms[i];
      if (!a.labelled || !a.pieces.length) { info.omitted += n; return; }
      /* the symbol: the full-size piece nearest the atom point */
      let symK = 0, best = Infinity;
      a.pieces.forEach((p, k) => {
        if (p.size !== FONT) return;
        const d = Math.hypot((p.box.x0 + p.box.x1) / 2 - a.X, (p.box.y0 + p.box.y1) / 2 - a.Y);
        if (d < best) { best = d; symK = k; }
      });
      const sb = a.pieces[symK].box, cx = (sb.x0 + sb.x1) / 2, cy = (sb.y0 + sb.y1) / 2;
      /* what the pairs keep away from: the bonds, and the atom's own H (a condensed O–H is a bond too) */
      const busy = a.nb.map(({ j }) => Math.atan2(atoms[j].Y - a.Y, atoms[j].X - a.X));
      a.pieces.forEach((p, k) => { if (k !== symK && p.size === FONT) busy.push(Math.atan2((p.box.y0 + p.box.y1) / 2 - cy, (p.box.x0 + p.box.x1) / 2 - cx)); });
      /* the free sides */
      const free = [];
      for (const D of DIRS) {
        const dx = Math.cos(D.ang), dy = Math.sin(D.ang);
        const tx = Math.abs(dx) > 1e-6 ? ((dx > 0 ? sb.x1 : sb.x0) - cx) / dx : Infinity;
        const ty = Math.abs(dy) > 1e-6 ? ((dy > 0 ? sb.y1 : sb.y0) - cy) / dy : Infinity;
        const r = Math.min(Math.abs(tx), Math.abs(ty)) + GAP + R;
        const px = cx + dx * r, py = cy + dy * r, ox = -dy * SEP / 2, oy = dx * SEP / 2;
        const dots = [{ x: px + ox, y: py + oy }, { x: px - ox, y: py - oy }];
        const box = boxOf(Math.min(dots[0].x, dots[1].x) - R, Math.min(dots[0].y, dots[1].y) - R, Math.max(dots[0].x, dots[1].x) + R, Math.max(dots[0].y, dots[1].y) + R);
        let hard = busy.some(b => angDiff(b, D.ang) < Math.PI / 5);
        if (!hard) for (const sg of segs) if (segMeetsBox(sg[0], sg[1], sg[2], sg[3], grow(box, LW))) { hard = true; break; }
        if (!hard) for (const bl of blocks) if (!(bl.i === a.i && bl.k === symK) && boxesMeet(box, bl.box)) { hard = true; break; }
        if (!hard) for (const t of taken) if (boxesMeet(grow(box, 0.5 * R), t)) { hard = true; break; }
        if (!hard) for (const o of atoms) if (o !== a && Math.hypot(px - o.X, py - o.Y) < 0.3 * BL) { hard = true; break; }
        if (!hard) free.push({ ang: D.ang, card: D.card, dots, box });
      }
      /* the set of n free sides with the lowest repulsion (pairs with each other and with the bonds), sides preferred */
      const want = Math.min(n, free.length);
      let pick = null, pickE = Infinity;
      const choose = (start, set) => {
        if (set.length === want) {
          let E = 0;
          for (let x = 0; x < set.length; x++) {
            for (const b of busy) E += energy(set[x].ang, b);
            for (let y = x + 1; y < set.length; y++) E += energy(set[x].ang, set[y].ang);
            if (!set[x].card) E += 0.35;
          }
          if (E < pickE) { pickE = E; pick = set.slice(); }
          return;
        }
        for (let k = start; k < free.length; k++) { set.push(free[k]); choose(k + 1, set); set.pop(); }
      };
      if (want) choose(0, []);
      for (const p of pick || []) { marks.push({ i, dots: p.dots, box: p.box }); taken.push(p.box); info.pairs++; }
      info.omitted += n - want;
    });
    return { marks, info };
  }

  function placeNotes(scene, texts) {
    const { atoms } = scene;
    const notes = [], info = { placed: 0, overlapping: 0, hOmitted: 0 };
    if (!texts) return { notes, info };
    const segs = bondSegs(scene);
    const blocks = [], zones = [];
    for (const a of atoms) {
      for (const p of a.pieces) blocks.push(grow(p.box, LABEL_GAP));
      for (const d of a.dots) blocks.push(boxOf(d.x - DOT_R - LABEL_GAP, d.y - DOT_R - LABEL_GAP, d.x + DOT_R + LABEL_GAP, d.y + DOT_R + LABEL_GAP));
      if (!a.labelled) blocks.push(boxOf(a.X - 0.08 * BL, a.Y - 0.08 * BL, a.X + 0.08 * BL, a.Y + 0.08 * BL));
      /* the sub- and superscript places: level with each line of the label, just left and right of it */
      if (a.labelled && a.extent) {
        for (const p of a.pieces) {
          if (p.size !== FONT) continue;                                 // a line of the label: its symbol or H
          const y0 = p.box.y0 - 0.35 * FONT, y1 = p.box.y1 + 0.35 * FONT;
          zones.push(boxOf(a.extent.x0 - 1.2 * FONT, y0, a.extent.x0, y1), boxOf(a.extent.x1, y0, a.extent.x1 + 1.2 * FONT, y1));
        }
      }
    }
    let heavy = 0;
    for (const a of atoms) if (!isPlainH(a, atoms)) heavy++;
    const skipH = heavy > 40;
    const h = 0.74 * NOTE;
    texts.forEach((text, i) => {
      if (text == null || text === '' || i >= atoms.length) return;
      text = String(text);
      const a = atoms[i];
      if (skipH && isPlainH(a, atoms)) { info.hOmitted++; return; }
      const w = textWidth(text, NOTE);
      /* the open side: away from the bonds (straight down for an atom without bonds, or when they balance out) */
      let ox = 0, oy = 0;
      for (const { j } of a.nb) { const dx = atoms[j].X - a.X, dy = atoms[j].Y - a.Y, L = Math.hypot(dx, dy) || 1; ox -= dx / L; oy -= dy / L; }
      const open = Math.hypot(ox, oy) < 0.1 ? Math.PI / 2 : Math.atan2(oy, ox);
      let best = null;
      for (const ring of [0, 1, 2]) {
        for (let s = 0; s < 24; s++) {
          const ang = s * Math.PI / 12, dx = Math.cos(ang), dy = Math.sin(ang);
          /* start just outside the label (or the atom point) in this direction */
          let r0 = 0.2 * BL;
          if (a.extent) {
            const ex = Math.max(Math.abs(a.extent.x0 - a.X), Math.abs(a.extent.x1 - a.X)), ey = Math.max(Math.abs(a.extent.y0 - a.Y), Math.abs(a.extent.y1 - a.Y));
            const tx = Math.abs(dx) > 1e-6 ? ((dx > 0 ? a.extent.x1 : a.extent.x0) - a.X) / dx : Infinity;
            const ty = Math.abs(dy) > 1e-6 ? ((dy > 0 ? a.extent.y1 : a.extent.y0) - a.Y) / dy : Infinity;
            r0 = Math.min(Math.abs(tx), Math.abs(ty), Math.hypot(ex, ey)) + LABEL_GAP + 0.05 * FONT;
          }
          const r = r0 + ring * 0.18 * BL + Math.abs(dx) * w / 2 + Math.abs(dy) * h / 2;
          const cx = a.X + dx * r, cy = a.Y + dy * r;
          const box = boxOf(cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2);
          let cost = ring * 3 + angDiff(ang, open) / Math.PI * 2, hard = 0;
          for (const z of zones) if (boxesMeet(box, z)) { hard++; break; }
          /* it must be clearly nearer its own atom than any other, or it could be read as that atom's */
          const own = Math.hypot(cx - a.X, cy - a.Y);
          for (const o of atoms) if (o !== a && Math.hypot(cx - o.X, cy - o.Y) < own / 0.85) { hard++; break; }
          for (const sg of segs) if (segMeetsBox(sg[0], sg[1], sg[2], sg[3], grow(box, LW))) hard++;
          for (const bl of blocks) if (boxesMeet(box, bl)) hard++;
          for (const n of notes) if (boxesMeet(grow(box, NOTE_GAP), n.box)) hard++;
          cost += 100 * hard;
          if (!best || cost < best.cost) best = { cost, hard, cx, cy, box };
        }
        if (best && !best.hard) break;
      }
      if (best.hard) { info.overlapping++; return; }                 // nowhere free: left out, not drawn over things
      notes.push({ i, text, x: best.cx - w / 2, base: best.cy + 0.36 * NOTE, box: best.box });
    });
    info.placed = notes.length;
    return { notes, info };
  }

  /* ---- SVG elements of a depiction ---- */
  function lineEl(x1, y1, x2, y2, color, s) {
    const at = { x1: r2(x1), y1: r2(y1), x2: r2(x2), y2: r2(y2), stroke: color };
    if (s.width) at['stroke-width'] = s.width;
    if (s.butt) at['stroke-linecap'] = 'butt';
    if (s.dash) at['stroke-dasharray'] = `${r2(0.1 * BL)} ${r2(0.08 * BL)}`;
    return el('line', at);
  }
  function drawBonds(scene, style) {
    const mono = !!(style && style.mono);
    const g = el('g', { class: 'mv-2d-bonds', 'stroke-width': mono ? TEXTBOOK.lw : LW, 'stroke-linecap': 'round', fill: 'none' });
    for (const b of scene.bonds) {
      const bg = el('g', { 'data-bond': b.k });
      for (const seg of b.segs) {
        const s = mono && seg.width === HASH_LW ? Object.assign({}, seg, { width: TEXTBOOK.hashLw }) : seg;
        if (s.c1 === s.c2) { bg.appendChild(lineEl(s.x1, s.y1, s.x2, s.y2, s.c1, s)); continue; }
        const mx = (s.x1 + s.x2) / 2, my = (s.y1 + s.y2) / 2;
        bg.appendChild(lineEl(s.x1, s.y1, mx, my, s.c1, s));
        bg.appendChild(lineEl(mx, my, s.x2, s.y2, s.c2, s));
      }
      for (const p of b.polys) bg.appendChild(el('polygon', { points: p.pts.map(q => r2(q[0]) + ',' + r2(q[1])).join(' '), fill: p.color, stroke: p.color, 'stroke-width': 0.6, 'stroke-linejoin': 'round' }));
      for (const p of b.paths) bg.appendChild(el('path', { d: p.d, stroke: p.color, 'stroke-width': 1.6 }));
      g.appendChild(bg);
    }
    return g;
  }
  function drawLabels(scene, style) {
    const g = el('g', { class: 'mv-2d-labels', 'font-family': FAMILY, 'font-weight': style && style.mono ? TEXTBOOK.weight : WEIGHT });
    for (const a of scene.atoms) {
      if (!a.labelled) continue;
      const t = el('text', { 'data-atom': a.i, fill: a.color });
      /* in reading order (line, then left to right), so the text copies as "H2N", "NH3+" */
      const line = p => Math.round((p.base - a.Y - BASE) / HLINE);
      const order = a.pieces.slice().sort((p, q) => (line(p) - line(q)) || (p.x - q.x));
      for (const p of order) t.appendChild(el('tspan', { x: r2(p.x), y: r2(p.base), 'font-size': r2(p.size) }, p.text));
      g.appendChild(t);
      for (const d of a.dots) g.appendChild(el('circle', { cx: r2(d.x), cy: r2(d.y), r: r2(DOT_R), fill: a.color }));
    }
    return g;
  }
  function drawElectrons(marks, style) {
    const g = el('g', { class: 'mv-2d-electrons' });
    for (const m of marks) for (const d of m.dots)
      g.appendChild(el('circle', { cx: r2(d.x), cy: r2(d.y), r: r2(DOT_R), fill: style && style.mono ? TEXTBOOK.ink : ELECTRON_INK, 'data-atom': m.i }));
    return g;
  }
  function drawNotes(notes) {
    const g = el('g', { class: 'mv-2d-notes', 'font-family': FAMILY, 'font-weight': WEIGHT, 'font-style': 'italic', 'font-size': r2(NOTE), fill: NOTE_INK });
    for (const n of notes) g.appendChild(el('text', { x: r2(n.x), y: r2(n.base), 'data-atom': n.i }, n.text));
    return g;
  }
  function drawHits(scene) {
    const g = el('g', { class: 'mv-2d-hits' });
    for (const b of scene.bonds) {
      const h = b.hit;
      g.appendChild(el('line', { x1: r2(h.x1), y1: r2(h.y1), x2: r2(h.x2), y2: r2(h.y2), stroke: '#000', 'stroke-opacity': 0, 'stroke-width': HIT_BOND, 'stroke-linecap': 'butt', 'pointer-events': 'stroke', 'data-hit': 'b' + b.k }));
    }
    for (const a of scene.atoms) {
      const ag = el('g', { 'data-hit': 'a' + a.i });
      ag.appendChild(el('circle', { cx: r2(a.X), cy: r2(a.Y), r: HIT_ATOM, fill: '#000', 'fill-opacity': 0, 'pointer-events': 'all' }));
      if (a.extent) {
        const e = grow(a.extent, 0.1 * FONT);
        ag.appendChild(el('rect', { x: r2(e.x0), y: r2(e.y0), width: r2(e.x1 - e.x0), height: r2(e.y1 - e.y0), fill: '#000', 'fill-opacity': 0, 'pointer-events': 'all' }));
      }
      g.appendChild(ag);
    }
    return g;
  }
  function drawHalo(scene, kind, spec) {
    const style = HALO[kind] || HALO.select;
    const color = (spec && spec.color) || style.color;
    const g = el('g', { class: 'mv-2d-hl-' + kind, opacity: style.opacity, fill: color, stroke: color, 'pointer-events': 'none' });
    if (!spec || !scene) return g;
    const byK = new Map(scene.bonds.map(b => [b.k, b]));
    for (const k of spec.bonds || []) {
      const b = byK.get(k);
      if (!b) continue;
      const P = scene.atoms[b.a], Q = scene.atoms[b.b];
      g.appendChild(el('line', { x1: r2(P.X), y1: r2(P.Y), x2: r2(Q.X), y2: r2(Q.Y), 'stroke-width': HALO_BOND, 'stroke-linecap': 'round' }));
    }
    for (const i of spec.atoms || []) {
      const a = scene.atoms[i];
      if (!a) continue;
      g.appendChild(el('circle', { cx: r2(a.X), cy: r2(a.Y), r: HALO_ATOM, stroke: 'none' }));
      if (a.extent) {
        const e = grow(a.extent, 0.18 * FONT), rr = Math.min(0.45 * FONT, (e.y1 - e.y0) / 2);
        g.appendChild(el('rect', { x: r2(e.x0), y: r2(e.y0), width: r2(e.x1 - e.x0), height: r2(e.y1 - e.y0), rx: r2(rr), stroke: 'none' }));
      }
    }
    return g;
  }

  /* =============================== the view =============================== */
  function create(container, opts = {}) {
    if (!container) throw new Error('Mol2DView.create: no container');
    const maxBondPx = opts.maxBondPx || 56;
    const style = { mono: !!opts.textbook };        // textbook: black lines and labels (setStyle changes it)
    let lastMol = null;
    const autoHover = opts.autoHover !== false;
    const svg = el('svg', { class: 'mv-2d', tabindex: 0, role: 'img', 'aria-label': '2D structure' });
    svg.style.cssText = 'display:block;width:100%;height:100%;touch-action:none;user-select:none;-webkit-user-select:none;outline:none;cursor:grab';
    const bg = el('rect', { x: 0, y: 0, width: '100%', height: '100%', fill: '#000', 'fill-opacity': 0 });
    const view = el('g', { class: 'mv-2d-view' });
    const leftOut = el('text', { class: 'mv-2d-leftout', x: 8, y: '100%', dy: -6, 'font-family': FAMILY, 'font-size': 11, fill: '#555', 'pointer-events': 'none' });
    svg.append(bg, view, leftOut);
    container.appendChild(svg);

    let scene = null, notesText = null, electronCounts = null;
    const specs = { group: null, select: null, hover: null };
    const layers = {};
    let s = 1, tx = 0, ty = 0, autoFit = true, needFit = true, lastW = 0, lastH = 0;
    let hoverHit = null, press = null, disposed = false;
    const pointers = new Map();
    let pinch = null;

    const apply = () => view.setAttribute('transform', `matrix(${s} 0 0 ${s} ${tx} ${ty})`);
    const size = () => { const r = svg.getBoundingClientRect(); return { w: r.width, h: r.height, left: r.left, top: r.top }; };
    const clampS = v => Math.min(600 / BL, Math.max(6 / BL, v));

    function rebuildLayers() {
      view.replaceChildren();
      if (!scene) return;
      layers.group = drawHalo(scene, 'group', specs.group);
      layers.select = drawHalo(scene, 'select', specs.select);
      layers.hover = drawHalo(scene, 'hover', specs.hover);
      layers.bonds = drawBonds(scene, style);
      layers.labels = drawLabels(scene, style);
      layers.notes = drawNotes(scene.notes);
      layers.electrons = drawElectrons(scene.electrons || [], style);
      for (const k of ['bonds', 'labels', 'notes', 'electrons']) layers[k].setAttribute('pointer-events', 'none');
      layers.hits = drawHits(scene);
      view.append(layers.group, layers.select, layers.hover, layers.bonds, layers.labels, layers.electrons, layers.notes, layers.hits);
    }
    /* the drawing's extent, including annotations */
    function extent() {
      const b = Object.assign({}, scene.bbox);
      for (const n of scene.notes.concat(scene.electrons || [])) { b.x0 = Math.min(b.x0, n.box.x0); b.y0 = Math.min(b.y0, n.box.y0); b.x1 = Math.max(b.x1, n.box.x1); b.y1 = Math.max(b.y1, n.box.y1); }
      return b;
    }

    function fit() {
      const { w, h } = size();
      if (!(w > 0 && h > 0)) { needFit = true; return; }
      needFit = false; autoFit = true; lastW = w; lastH = h;
      if (!scene) { s = 1; tx = w / 2; ty = h / 2; apply(); return; }
      const b = extent(), m = 16, foot = leftOut.textContent ? 18 : 0;   // room for the 'left out' line
      const bw = Math.max(b.x1 - b.x0, 1), bh = Math.max(b.y1 - b.y0, 1);
      s = clampS(Math.min((w - 2 * m) / bw, (h - 2 * m - foot) / bh, maxBondPx / BL));
      tx = w / 2 - s * (b.x0 + b.x1) / 2; ty = (h - foot) / 2 - s * (b.y0 + b.y1) / 2;
      apply();
    }
    function resize() {
      const { w, h } = size();
      if (!(w > 0 && h > 0)) return;                                  // hidden: wait until shown
      if (autoFit || needFit) { fit(); return; }
      tx += (w - lastW) / 2; ty += (h - lastH) / 2;                    // keep the middle where it was
      lastW = w; lastH = h;
      apply();
    }
    const ro = typeof ResizeObserver !== 'undefined' ? new ResizeObserver(() => { if (!disposed) resize(); }) : null;
    if (ro) ro.observe(container);

    function zoomAt(f, cx, cy) {
      const ns = clampS(s * f);
      if (ns === s) return;
      tx = cx - (cx - tx) * ns / s; ty = cy - (cy - ty) * ns / s; s = ns;
      autoFit = false; apply();
    }

    /* ---- hit testing. The transparent hit layer (data-hit="a<i>" / "b<k>") says which atoms and bonds are
       under the pointer. Of those: an atom's label always picks the atom; along its own bond an atom reaches
       at most 35 % of the bond, so the middle of every bond picks the bond even when bonds are short; against
       any other bond the nearest wins (in a crowded layout an atom's circle can reach over another bond). ---- */
    function hitAt(cx, cy, target) {
      if (!scene) return null;
      let els = null;
      try { els = document.elementsFromPoint(cx, cy); } catch (err) { els = null; }
      if (!els) els = target ? [target] : [];
      const ids = new Set();
      for (const n of els) {
        const h = n && n.closest ? n.closest('[data-hit]') : null;
        if (h && view.contains(h)) ids.add(h.getAttribute('data-hit'));
      }
      if (!ids.size) return null;
      const r = size(), x = (cx - r.left - tx) / s, y = (cy - r.top - ty) / s;
      let atom = null, ad = Infinity;
      const near = [];                                                 // candidate bonds with their distances
      for (const id of ids) {
        if (id[0] === 'a') {
          const a = scene.atoms[+id.slice(1)];
          if (!a) continue;
          const e = a.extent && grow(a.extent, 0.1 * FONT);
          const d = e && x >= e.x0 && x <= e.x1 && y >= e.y0 && y <= e.y1 ? 0 : Math.hypot(x - a.X, y - a.Y);
          if (d < ad) { ad = d; atom = a; }
        } else {
          const b = scene.byK.get(+id.slice(1));
          if (!b) continue;
          const h = b.hit, dx = h.x2 - h.x1, dy = h.y2 - h.y1, L2 = dx * dx + dy * dy || 1;
          const t = Math.max(0, Math.min(1, ((x - h.x1) * dx + (y - h.y1) * dy) / L2));
          near.push({ b, d: Math.hypot(x - h.x1 - t * dx, y - h.y1 - t * dy) });
        }
      }
      /* the bonds that beat the atom (all of them when there is no atom); the nearest of those wins */
      const beats = near.filter(({ b, d }) => {
        if (!atom) return true;
        if (b.a === atom.i || b.b === atom.i) { const P = scene.atoms[b.a], Q = scene.atoms[b.b]; return ad > 0 && ad > 0.35 * Math.hypot(Q.X - P.X, Q.Y - P.Y); }
        return d < ad;
      });
      if (!beats.length) return atom ? { type: 'atom', index: atom.i } : null;
      const best = beats.reduce((p, q) => (q.d < p.d ? q : p));
      return { type: 'bond', index: best.b.k };
    }
    const sameHit = (p, q) => (!p && !q) || (p && q && p.type === q.type && p.index === q.index);
    function setHover(hit, cx, cy) {
      const changed = !sameHit(hit, hoverHit);
      if (changed) {
        hoverHit = hit;
        if (autoHover) highlight('hover', hit ? (hit.type === 'atom' ? { atoms: [hit.index] } : { bonds: [hit.index] }) : null);
        if (!press) svg.style.cursor = hit ? 'pointer' : 'grab';
      }
      if (opts.onHover && (hit || changed)) opts.onHover(hit, cx, cy);
    }

    const onPointerDown = e => {
      if (e.pointerType === 'mouse' && e.button !== 0) return;
      try { svg.focus({ preventScroll: true }); } catch (err) { /* old browsers */ }
      pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
      if (pointers.size === 1) press = { id: e.pointerId, x: e.clientX, y: e.clientY, tx, ty, moved: false, hit: hitAt(e.clientX, e.clientY, e.target) };
      else if (pointers.size === 2) {
        const [p, q] = [...pointers.values()];
        const r = size(), cx = (p.x + q.x) / 2 - r.left, cy = (p.y + q.y) / 2 - r.top;
        pinch = { d: Math.hypot(p.x - q.x, p.y - q.y) || 1, s, mx: (cx - tx) / s, my: (cy - ty) / s };   // the model point under the fingers
        if (press) press.moved = true;
      }
      try { svg.setPointerCapture(e.pointerId); } catch (err) { /* ignore */ }
    };
    const onPointerMove = e => {
      if (pointers.has(e.pointerId)) pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
      if (pinch && pointers.size >= 2) {
        const [p, q] = [...pointers.values()], r = size();
        const d = Math.hypot(p.x - q.x, p.y - q.y) || 1, cx = (p.x + q.x) / 2 - r.left, cy = (p.y + q.y) / 2 - r.top;
        s = clampS(pinch.s * d / pinch.d);
        tx = cx - s * pinch.mx; ty = cy - s * pinch.my;
        autoFit = false; apply();
        return;
      }
      if (press && e.pointerId === press.id) {
        const dx = e.clientX - press.x, dy = e.clientY - press.y;
        if (!press.moved && Math.hypot(dx, dy) > 4) {
          press.moved = true; svg.style.cursor = 'grabbing';
          if (hoverHit) setHover(null, e.clientX, e.clientY);
        }
        if (press.moved) { tx = press.tx + dx; ty = press.ty + dy; autoFit = false; apply(); }
        return;
      }
      if (!pointers.size) setHover(hitAt(e.clientX, e.clientY, e.target), e.clientX, e.clientY);
    };
    const onPointerUp = e => {
      const wasPress = press && e.pointerId === press.id;
      pointers.delete(e.pointerId);
      try { svg.releasePointerCapture(e.pointerId); } catch (err) { /* ignore */ }
      if (pointers.size < 2) pinch = null;
      if (wasPress) {
        const clicked = !press.moved && e.type === 'pointerup';
        const hit = press.hit;
        press = null;
        svg.style.cursor = 'grab';
        if (clicked && opts.onClick) opts.onClick(hit, e);
        if (e.pointerType === 'mouse') setHover(hitAt(e.clientX, e.clientY, null), e.clientX, e.clientY);
      }
      if (!pointers.size) press = null;
    };
    const onLeave = e => { if (!pointers.size) setHover(null, e.clientX, e.clientY); };
    const onEnter = () => { if (needFit) resize(); };                  // in case the ResizeObserver has not run yet
    const onWheel = e => {
      if (!scene) return;
      e.preventDefault();
      const r = size();
      const dy = e.deltaMode === 1 ? e.deltaY * 33 : e.deltaMode === 2 ? e.deltaY * 400 : e.deltaY;
      zoomAt(Math.exp(-dy * 0.0015), e.clientX - r.left, e.clientY - r.top);
    };
    const onDbl = e => { if (!hitAt(e.clientX, e.clientY, e.target)) fit(); };
    const onKey = e => {
      const { w, h } = size();
      const step = 40;
      switch (e.key) {
        case '+': case '=': zoomAt(1.25, w / 2, h / 2); break;
        case '-': case '_': zoomAt(0.8, w / 2, h / 2); break;
        case '0': fit(); break;
        case 'ArrowLeft': tx += step; autoFit = false; apply(); break;
        case 'ArrowRight': tx -= step; autoFit = false; apply(); break;
        case 'ArrowUp': ty += step; autoFit = false; apply(); break;
        case 'ArrowDown': ty -= step; autoFit = false; apply(); break;
        default: return;
      }
      e.preventDefault();
    };
    svg.addEventListener('pointerdown', onPointerDown);
    svg.addEventListener('pointermove', onPointerMove);
    svg.addEventListener('pointerup', onPointerUp);
    svg.addEventListener('pointercancel', onPointerUp);
    svg.addEventListener('pointerleave', onLeave);
    svg.addEventListener('pointerenter', onEnter);
    svg.addEventListener('wheel', onWheel, { passive: false });
    svg.addEventListener('dblclick', onDbl);
    svg.addEventListener('keydown', onKey);

    function setMolecule(mol) {
      specs.group = specs.select = specs.hover = null;
      notesText = null; electronCounts = null; hoverHit = null;
      lastMol = mol;
      scene = mol && Array.isArray(mol.atoms) && mol.atoms.length ? depict(mol, style) : null;
      rebuildLayers();
      showNotesLeftOut();
      fit();
      if (scene && scene.warnings.length && window.console) console.warn('Mol2DView:', scene.warnings.join('; '));
      return { warnings: scene ? scene.warnings.slice() : [] };
    }
    /* setStyle({textbook}): redraw in the other style, keeping the view (zoom, pan), highlights and annotations */
    function setStyle(o = {}) {
      if ('textbook' in o) style.mono = !!o.textbook;
      if (!lastMol || !scene) return;
      scene = depict(lastMol, style);
      applyElectrons();
      rebuildLayers();
      if (notesText) setAnnotations(notesText);
    }
    /* setElectrons(counts | null): counts[i] = how many LONE PAIRS to draw around atom i (the bonds are the bonding
       pairs). Returns {pairs, omitted}: how many were drawn, and how many had no free side. */
    function applyElectrons() {
      if (!scene) return { pairs: 0, omitted: 0 };
      const r = placeElectrons(scene, electronCounts);
      scene.electrons = r.marks; scene.electronsInfo = r.info;
      return r.info;
    }
    function setElectrons(counts) {
      electronCounts = Array.isArray(counts) ? counts.slice() : null;
      if (!scene) return { pairs: 0, omitted: 0 };
      const info = applyElectrons();
      const fresh = drawElectrons(scene.electrons, style);
      fresh.setAttribute('pointer-events', 'none');
      view.replaceChild(fresh, layers.electrons);
      layers.electrons = fresh;
      if (notesText) setAnnotations(notesText); else { showNotesLeftOut(); if (autoFit) fit(); }
      return Object.assign({}, info);
    }
    function setAnnotations(texts) {
      notesText = Array.isArray(texts) ? texts.slice() : null;
      if (!scene) return { placed: 0, overlapping: 0, hOmitted: 0 };
      const r = placeNotes(scene, notesText);
      scene.notes = r.notes; scene.notesInfo = r.info;
      const fresh = drawNotes(scene.notes);
      fresh.setAttribute('pointer-events', 'none');
      view.replaceChild(fresh, layers.notes);
      layers.notes = fresh;
      showNotesLeftOut();
      if (autoFit) fit();
      return Object.assign({}, r.info);
    }
    /* on screen only: say that some annotations are left out, so an atom without one is not taken for one
       without a value */
    function showNotesLeftOut() {
      const inf = scene && scene.notesInfo, parts = [];
      if (inf && inf.overlapping) parts.push(`${inf.overlapping} label${inf.overlapping === 1 ? '' : 's'} left out where there is no free space`);
      if (inf && inf.hOmitted) parts.push('H labels left out (more than 40 heavy atoms)');
      const e = scene && scene.electronsInfo;
      if (e && e.omitted) parts.push(`${e.omitted} lone pair${e.omitted === 1 ? '' : 's'} left out where there is no free space`);
      leftOut.textContent = parts.length ? parts.join(' · ') + ' — hover an atom for its values' : '';
    }
    function highlight(kind, spec) {
      if (!(kind in specs)) throw new Error('Mol2DView.highlight: unknown kind ' + kind);
      specs[kind] = spec && ((spec.atoms && spec.atoms.length) || (spec.bonds && spec.bonds.length)) ? { atoms: (spec.atoms || []).slice(), bonds: (spec.bonds || []).slice(), color: spec.color } : null;
      if (!scene) return;
      const fresh = drawHalo(scene, kind, specs[kind]);
      view.replaceChild(fresh, layers[kind]);
      layers[kind] = fresh;
    }

    /* ---- export: a standalone SVG of the drawing (no hit areas; halos only if asked) ---- */
    // the export's words for the notes the drawing leaves out (the screen's line, without "hover")
    function leftOutSentence(inf) {
      const parts = [];
      if (inf && inf.overlapping) parts.push(`${inf.overlapping} label${inf.overlapping === 1 ? '' : 's'} left out where there is no free space`);
      if (inf && inf.hOmitted) parts.push('H labels left out (more than 40 heavy atoms)');
      if (scene && scene.electronsInfo && scene.electronsInfo.omitted) parts.push(`${scene.electronsInfo.omitted} lone pair${scene.electronsInfo.omitted === 1 ? '' : 's'} left out where there is no free space`);
      return parts.length ? 'Not every atom is labelled: ' + parts.join('; ') + '.' : '';
    }
    function exportSVG(o = {}) { return svgExport(o).text; }
    function svgExport(o) {
      if (!scene) return { text: `<svg xmlns="${NS}" width="1" height="1"/>`, w: 1, h: 1 };
      const b = extent(), m = o.highlights ? 0.4 * BL : 0.25 * BL;
      const x0 = Math.floor(b.x0 - m), y0 = Math.floor(b.y0 - m);
      let w = Math.ceil(b.x1 + m) - x0, h = Math.ceil(b.y1 + m) - y0;
      /* a caption (e.g. the method behind the annotations) in lines under the drawing, never over it; notes left
         out (no free space; H notes with many heavy atoms) are said there too, as on screen: an atom without a note
         must not read as one without a value */
      const capLines = [], CAP = NOTE, CAP_LH = 1.35 * NOTE, CAP_X = 0.25 * BL;
      const paras = o.caption != null && o.caption !== '' ? [].concat(o.caption).map(String) : [];
      const left = leftOutSentence(scene.notesInfo);
      if (left) paras.push(left);
      if (paras.length) {
        const maxW = Math.max(w - 2 * CAP_X, 18 * CAP);
        for (const para of paras) {
          let line = '';
          for (const word of para.split(/\s+/).filter(Boolean)) {
            const t = line ? line + ' ' + word : word;
            if (line && textWidth(t, CAP) > maxW) { capLines.push(line); line = word; } else line = t;
          }
          if (line) capLines.push(line);
        }
      }
      const drawH = h;
      if (capLines.length) {
        w = Math.max(w, Math.ceil(Math.max(...capLines.map(t => textWidth(t, CAP))) + 2 * CAP_X));
        h += Math.ceil(capLines.length * CAP_LH + 0.3 * BL);
      }
      const root = el('svg', { width: w, height: h, viewBox: `${x0} ${y0} ${w} ${h}` });   // the serializer adds xmlns
      if (!o.transparent) root.appendChild(el('rect', { x: x0, y: y0, width: w, height: h, fill: '#fff' }));
      if (capLines.length) {
        const g = el('g', { class: 'mv-2d-caption', 'font-family': FAMILY, 'font-size': r2(CAP), fill: '#444' });
        capLines.forEach((t, k) => g.appendChild(el('text', { x: r2(x0 + CAP_X), y: r2(y0 + drawH + (k + 1) * CAP_LH) }, t)));
        root.appendChild(g);
      }
      if (o.highlights) {
        if (specs.group) root.appendChild(drawHalo(scene, 'group', specs.group));
        if (specs.select) root.appendChild(drawHalo(scene, 'select', specs.select));
      }
      root.appendChild(drawBonds(scene, style));
      root.appendChild(drawLabels(scene, style));
      if (scene.electrons && scene.electrons.length) root.appendChild(drawElectrons(scene.electrons, style));
      if (scene.notes.length) root.appendChild(drawNotes(scene.notes));
      return { text: new XMLSerializer().serializeToString(root), w, h };
    }
    function exportPNG(scale = 2, o = {}) {
      return new Promise((resolve, reject) => {
        if (!scene) { reject(new Error('Nothing to export.')); return; }
        const { text, w, h } = svgExport(o);
        const img = new Image();
        img.onload = () => {
          try {
            const cv = document.createElement('canvas');
            cv.width = Math.max(1, Math.round(w * scale)); cv.height = Math.max(1, Math.round(h * scale));
            const ctx = cv.getContext('2d');
            if (!o.transparent) { ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, cv.width, cv.height); }
            ctx.drawImage(img, 0, 0, cv.width, cv.height);
            cv.toBlob(blob => (blob ? resolve(blob) : reject(new Error('PNG encoding failed'))), 'image/png');
          } catch (err) { reject(err); }
        };
        img.onerror = () => reject(new Error('The drawing could not be rendered to PNG.'));
        img.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(text);
      });
    }

    function clientPosition(type, i) {
      if (!scene) return null;
      let x, y;
      if (type === 'atom') { const a = scene.atoms[i]; if (!a) return null; x = a.X; y = a.Y; }
      else { const b = scene.bonds.find(q => q.k === i); if (!b) return null; x = (b.hit.x1 + b.hit.x2) / 2; y = (b.hit.y1 + b.hit.y2) / 2; }
      const r = size();
      return { x: r.left + tx + s * x, y: r.top + ty + s * y };
    }

    function dispose() {
      if (disposed) return;
      disposed = true;
      if (ro) ro.disconnect();
      svg.removeEventListener('pointerdown', onPointerDown);
      svg.removeEventListener('pointermove', onPointerMove);
      svg.removeEventListener('pointerup', onPointerUp);
      svg.removeEventListener('pointercancel', onPointerUp);
      svg.removeEventListener('pointerleave', onLeave);
      svg.removeEventListener('pointerenter', onEnter);
      svg.removeEventListener('wheel', onWheel);
      svg.removeEventListener('dblclick', onDbl);
      svg.removeEventListener('keydown', onKey);
      svg.remove();
      scene = null;
    }

    apply();
    return { setMolecule, setStyle, setAnnotations, setElectrons, highlight, fit, resize, exportSVG, exportPNG, dispose, clientPosition, element: svg };
  }

  window.Mol2DView = { create };
})();
