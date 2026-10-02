/* MolDraw — draw a Chem graph as a small SVG picture (for result pages, not editing).
     MolDraw.svg(graph, width, height, opts) → <svg> element, the molecule laid out and fitted inside; it carries
       atomXY (atom id → {x, y, labelled}) and scale, so arrows and marks can be drawn on top (see svg() below).
     MolDraw.layoutKeepingStereo(graph) → Chem's layout without losing R/S or E/Z (shared with the sketcher).
   The graph is cloned and laid out, so the caller's coordinates are untouched. */
(() => {
  'use strict';
  const NS = 'http://www.w3.org/2000/svg';
  const el = (name, attrs = {}, text) => {
    const n = document.createElementNS(NS, name);
    for (const k in attrs) n.setAttribute(k, attrs[k]);
    if (text != null) n.textContent = text;
    return n;
  };


  /* Lay a molecule out again WITHOUT losing its stereochemistry.

     A wedge only means something next to the coordinates it was drawn on, so re-laying a
     molecule out and then putting the old wedges back can silently produce the mirror image
     (that is how L-alanine came back as D-alanine). Chem already solves this when it redraws
     an AGILES parse: read the configuration first (R/S per centre, E/Z per double bond),
     lay out, then draw that configuration again on the new coordinates and check it. If any
     centre or double bond would come out different, the old picture is kept instead.

     `graph` is modified in place; returns true if the new layout was kept. */
  function layoutKeepingStereo(graph) {
    const C = window.Chem;
    const rs = new Map(), ez = [];
    /* a graph with no real coordinates yet (everything on one spot) says nothing about
       configuration, so there is nothing to preserve: just lay it out */
    const xsAll = graph.atoms.map(a => a.x || 0), ysAll = graph.atoms.map(a => a.y || 0);
    const flat = Math.max(...xsAll) - Math.min(...xsAll) < 1e-6 && Math.max(...ysAll) - Math.min(...ysAll) < 1e-6;
    try {
      if (flat) throw new Error('no coordinates');
      graph.atoms.forEach(a => { if (C.isStereocenter(graph, a)) { const r = C.assignRS(graph, a.id); if (r && r.label) rs.set(a.id, r.label); } });
      graph.bonds.forEach(b => { if (b.order === 2) { const e = C.assignEZ(graph, b); if (e && e.label) ez.push({ a: b.a, b: b.b, label: e.label }); } });
    } catch (e) { rs.clear(); ez.length = 0; }

    const xy = graph.atoms.map(a => ({ x: a.x, y: a.y }));
    const marks = graph.bonds.map(b => ({ stereo: b.stereo, narrow: b.narrow }));
    const restore = () => {
      graph.atoms.forEach((a, i) => { a.x = xy[i].x; a.y = xy[i].y; });
      graph.bonds.forEach((b, i) => {
        if (marks[i].stereo) { b.stereo = marks[i].stereo; b.narrow = marks[i].narrow; }
        else { delete b.stereo; delete b.narrow; }
      });
    };

    try {
      C.layoutGraph(graph);
      if (!rs.size && !ez.length) return true;                    // nothing to preserve
      graph.bonds.forEach(b => { delete b.stereo; delete b.narrow; });
      for (const x of ez) {
        const bd = C.bondBetween(graph, x.a, x.b), e = bd && C.assignEZ(graph, bd);
        if (e && e.label !== x.label) C.flipAlkeneEnd(graph, x.a, x.b);
      }
      graph.atoms.forEach(a => { if (rs.has(a.id)) a.__rs = rs.get(a.id); });
      C.applyParsedStereo(graph);
      graph.atoms.forEach(a => { delete a.__rs; });
      for (const [id, label] of rs) { const r = C.assignRS(graph, id); if (!r || r.label !== label) { restore(); return false; } }
      for (const x of ez) {
        const bd = C.bondBetween(graph, x.a, x.b), e = bd && C.assignEZ(graph, bd);
        if (!e || e.label !== x.label) { restore(); return false; }
      }
      return true;
    } catch (e) { restore(); return false; }
  }

  /* svg(graph, width, height, opts)
       opts.layout false  keep the graph's own coordinates (a mechanism step parks its reagent deliberately, and
                          laying it out again would move it); the picture is still scaled and centred to fit
       opts.pad           extra margin in px, for arrows drawn outside the structure
     The returned <svg> carries `atomXY`, a Map of atom id → {x, y} in its own coordinates, so a caller can draw
     curly arrows or other marks on top, and `scale`, the px per graph unit it used. */
  function svg(graph, width = 180, height = 120, opts = {}) {
    const C = window.Chem;
    const root = el('svg', { viewBox: `0 0 ${width} ${height}`, width, height, class: 'moldraw' });
    root.atomXY = new Map();
    root.scale = 1;
    if (!graph || !graph.atoms || !graph.atoms.length) return root;
    const g = { atoms: graph.atoms.map(a => ({ ...a })), bonds: graph.bonds.map(b => ({ ...b })), nextId: graph.nextId };
    if (g.bonds.length && opts.layout !== false) layoutKeepingStereo(g);
    const atom = id => g.atoms.find(a => a.id === id);
    const bondsOf = id => g.bonds.filter(b => b.a === id || b.b === id);
    const neighbours = id => bondsOf(id).map(b => atom(b.a === id ? b.b : b.a));

    /* fit: scale so the molecule fills the box with a margin, bond length capped at 30px */
    const xs = g.atoms.map(a => a.x), ys = g.atoms.map(a => a.y);
    const w = Math.max(...xs) - Math.min(...xs), h = Math.max(...ys) - Math.min(...ys);
    const lens = g.bonds.map(b => Math.hypot(atom(b.a).x - atom(b.b).x, atom(b.a).y - atom(b.b).y)).filter(l => l > 0);
    const meanLen = lens.length ? lens.reduce((s, l) => s + l, 0) / lens.length : 30;
    const margin = 18 + (opts.pad || 0);
    const scale = Math.min(30 / meanLen, (width - 2 * margin) / (w || 1), (height - 2 * margin) / (h || 1));
    const cx = (Math.min(...xs) + Math.max(...xs)) / 2, cy = (Math.min(...ys) + Math.max(...ys)) / 2;
    g.atoms.forEach(a => { a.x = width / 2 + (a.x - cx) * scale; a.y = height / 2 + (a.y - cy) * scale; });
    root.scale = scale;
    g.atoms.forEach(a => root.atomXY.set(a.id, { x: a.x, y: a.y, labelled: false }));
    const fontPx = Math.max(9, Math.min(13, 13 * scale * meanLen / 30));

    const labelled = a => a.element !== 'C' || a.charge || a.radical || bondsOf(a.id).length === 0;
    const trim = fontPx * 0.55;
    for (const b of g.bonds) {
      const p = atom(b.a), q = atom(b.b);
      const dx = q.x - p.x, dy = q.y - p.y, L = Math.hypot(dx, dy) || 1;
      const ux = dx / L, uy = dy / L, nx = -uy, ny = ux;
      const tp = labelled(p) ? trim : 0, tq = labelled(q) ? trim : 0;
      const x1 = p.x + ux * tp, y1 = p.y + uy * tp, x2 = q.x - ux * tq, y2 = q.y - uy * tq;
      const grp = el('g', { class: 'bond' });
      if (b.stereo && b.order === 1) {
        const flip = b.narrow === q.id;
        const sx = flip ? x2 : x1, sy = flip ? y2 : y1, ex = flip ? x1 : x2, ey = flip ? y1 : y2;
        const wdt = 3 * scale * meanLen / 30 + 1;
        if (b.stereo === 'wedge') grp.appendChild(el('polygon', { class: 'wedge', points: `${sx},${sy} ${ex + nx * wdt},${ey + ny * wdt} ${ex - nx * wdt},${ey - ny * wdt}` }));
        else { const steps = 6; for (let k = 1; k <= steps; k++) { const t = k / steps, ww = wdt * t, px = sx + (ex - sx) * t, py = sy + (ey - sy) * t; grp.appendChild(el('line', { x1: px + nx * ww, y1: py + ny * ww, x2: px - nx * ww, y2: py - ny * ww })); } }
        root.appendChild(grp); continue;
      }
      const off = 3.6 * Math.min(1, scale * meanLen / 30);
      let offsets = { 1: [0], 2: [0, off * 1.7], 3: [-off * 1.7, 0, off * 1.7] }[b.order] || [0];
      if (b.order === 2) {
        const others = [...neighbours(p.id), ...neighbours(q.id)].filter(o => o !== p && o !== q);
        if (!others.length) offsets = [-off, off];
        else { const side = others.reduce((s, o) => s + Math.sign((o.x - p.x) * nx + (o.y - p.y) * ny), 0); offsets = [0, side < 0 ? -off * 1.7 : off * 1.7]; }
      }
      for (const o of offsets) {
        const shrink = b.order === 2 && o !== 0 ? 3 : 0;
        grp.appendChild(el('line', { x1: x1 + nx * o + ux * shrink, y1: y1 + ny * o + uy * shrink, x2: x2 + nx * o - ux * shrink, y2: y2 + ny * o - uy * shrink }));
      }
      root.appendChild(grp);
    }
    for (const a of g.atoms) {
      if (!labelled(a)) continue;
      root.atomXY.get(a.id).labelled = true;
      const hN = a.element === 'H' ? 0 : C.implicitH(g, a);
      const text = a.element + (hN ? 'H' : '') + (hN > 1 ? hN : '') + (a.charge ? (a.charge > 0 ? '+' : '−') : '');
      const bw = text.length * fontPx * 0.62 + 4;
      root.appendChild(el('rect', { class: 'label-bg', x: a.x - bw / 2, y: a.y - fontPx * 0.65, width: bw, height: fontPx * 1.3, rx: 2 }));
      const t = el('text', { class: 'label', x: a.x, y: a.y, style: `font-size:${fontPx}px`, 'data-atom': a.id }, a.element);
      if (hN) t.appendChild(el('tspan', { style: `font-size:${fontPx * 0.72}px` }, 'H' + (hN > 1 ? hN : '')));
      if (a.charge) t.appendChild(el('tspan', { style: `font-size:${fontPx * 0.7}px;baseline-shift:super` }, (Math.abs(a.charge) > 1 ? Math.abs(a.charge) : '') + (a.charge > 0 ? '+' : '−')));
      root.appendChild(t);
    }
    return root;
  }
  window.MolDraw = { svg, layoutKeepingStereo };
})();
