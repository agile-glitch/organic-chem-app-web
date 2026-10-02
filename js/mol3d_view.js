/* Mol3DView — the Molecule tab's 3D viewer. three.js r185 (vendor/three/three.js, window.THREE) is loaded by MolView
   before create() is called; nothing here touches THREE until then. It holds no chemistry: elements, radii, bond
   orders, rings and partial charges come from the caller in record order (the Molecule tab's one numbering: the
   molecule's atoms, then their hydrogens grouped by parent), and every index it reports back is that record index.

   Mol3DView.create(container, opts) → view
     opts: { onHover(hit, clientX, clientY), onClick(hit, event), onContextLost(), onContextRestored(),
             onSpinChange(on) (Space key), onSurface(result|null, error|null) (surface redone after setCoords),
             background ('#ffffff') }            hit = {type:'atom', index} | {type:'bond', index} | null
   view.setMolecule({atoms:[{el, z, rvdw, rcov, color?, isH?}], bonds:[{a, b, order:1|2|3, aromatic}], rings})
   view.setCoords(Float64Array(3N), {keepView = true})     the first call after setMolecule frames the molecule
   view.setStyle('ballstick'|'spacefill'|'sticks'|'wireframe')   view.setShowH(bool)   view.setLabels(string[N] | null)
   view.setElectrons(number[N] | null)   lone PAIRS per atom, drawn as pairs of dots around it (see lonePairDirs)
   view.setSpin(bool)   view.resetView()   view.setOrthographic(bool)
   view.setSurface(null | {mode:'plain'|'potential', charges?, opacity = 0.75, clamp?}) → Promise<{vmin, vmax, clamp, …} | null>
   view.highlight('hover'|'select'|'group', {atoms, bonds} | null)   view.setMeasurements([{atoms:[i,j(,k(,l))], text}] | null)
   view.exportPNG(scale = 2, {legend: {title, minLabel, maxLabel, midLabel?, note}, caption, transparent}) → Promise<Blob>
       (transparent: no background, the picture's own alpha only)
       the view with its labels; a legend (red–white–blue bar with the caller's labels and note) and/or a caption
       (string or lines) go into the bottom-left corner if it shows only background, else into a strip added
       below, never over the molecule. When the image leaves labels out (overlaps; H labels with more than 40
       heavy atoms), a sentence saying so is added to the caption, as the screen says it. Without any of these the
       image is exactly the view.
   view.pick(clientX, clientY) → hit   view.setPickMode('all'|'atoms')   view.resize()   view.dispose()
   Mol3DView.espColor(V, clamp) → '#rrggbb': the surface colour of potential V (for the legend).
   Behaviour the caller can rely on:
   - The hovered atom or bond gets the 'hover' halo by itself; highlight('hover', …) overrides it until the pointer moves.
     onHover fires once per frame while the pointer moves (and after a new molecule/conformer), with null on leave.
   - setMolecule clears labels, measurements, highlights and the surface (charges belong to one molecule).
   - "Show H" hides atom.isH when given, else each hydrogen bonded to exactly one non-hydrogen; hidden atoms and their
     bonds are neither drawn nor picked. The surface always includes every atom.
   - setSurface throws out (rejects) rather than guess: an atom without a vdW radius, or a missing/NaN charge in
     'potential' mode. Before any coordinates it resolves null and builds the surface when they arrive. New
     coordinates drop the old surface at once and rebuild it, reporting through onSurface. A superseded request
     resolves null. Result: {vmin, vmax, clamp, autoClamp, fixedClamp, p2, p98, units, vertices, triangles, gridStep, ms}.
   - Picking is analytic (ray–sphere, ray–cylinder on the pieces actually drawn), thin pieces widened to 4 px, plus a
     flat pick-only strip across the gaps between the pieces of a multiple bond; a bond wins only where it is in
     front of the atom; space-filling picks atoms only. setPickMode('atoms') (the Molecule tab's Measure mode) picks
     atoms only, in every style, for hover and click: the atom under the pointer even where a bond is drawn in front
     of it; on a bond, the end nearer the pointer; elsewhere, an atom whose edge is within 12 px.
   - Spin and the glide after a drag run only while the view can be seen: when its box is not displayed (another
     tab of the app) or the page is hidden, nothing is drawn and no frame is asked for; it carries on when shown.
   - Labels never overlap: measurement labels are placed first, then atom labels nearest first (the hovered atom,
     then heavy atoms, then H); one whose box would overlap a label already placed is left out, and the overlay
     says how many; so does it for labels of atoms hidden behind other atoms; with more than 40 heavy atoms, H labels
     show only on hover. stats().labels = {shown, overlapping, hOnHover, hidden}.
   - Controls: drag rotates (trackball), Shift+drag rolls, right/middle/Ctrl+drag pans, wheel zooms; touch: one finger
     rotates, two pinch-zoom, pan and roll. Keys while the canvas has focus: arrows rotate 10° (Shift 1°), + − zoom,
     0/Home reset, Space spin. An untouched view is framed again when the style, H or surface changes what must fit.

   Where the numbers come from:
   - Colours: Jmol's default element colours, read from window.MolData.element(el).color; JMOL below is a copy of the
     common ones for pages without MolData (tests/molview/view3d_test.html checks that the two agree).
   - Sizes: space-filling spheres are the caller's van der Waals radii. Ball-and-stick balls (0.23 × vdW) and the
     stick and bond radii are picture choices, not data.
   - Surface: the union of the van der Waals spheres, found as the zero level of the signed distance to that union
     by marching cubes on a 0.25–0.35 Å grid (at most 200³ points; coarser for bigger molecules). Each vertex is the
     exact zero of the distance on its grid edge, so it lies on the union; normals are the exact sphere normals.
     Crevices narrower than the grid are bridged and enclosed voids dropped, so the drawn area is 1–4 % below the
     exact vdW area (measured on tests/molview/view3d_test.html against a dense Shrake–Rupley).
   - Potential: V = 332.0637 Σ qᵢ/rᵢ with q in e and r in Å gives kcal/mol per e (332.0637 = e²/4πε₀ in
     kcal·Å/mol). It is the Coulomb potential of the caller's point charges, not a quantum-chemical ESP. Colours
     run white → red for negative and white → blue for positive, clamped symmetrically at the 98th percentile of |V|
     on the surface unless the caller fixes the clamp.
   Rendering choices, checked against this three.js build (r185): colour management is on (sRGB output; hex colours
   are read as sRGB and stored linear), lights use physical units (three.js's default since r155; the intensities
   below are set for that), and an InstancedMesh's bounding sphere is recomputed after its instances move (three.js
   computes it once, and frustum culling would otherwise use the old one).

   Marching-cubes tables (EDGE_TABLE, TRI_TABLE) are copied from three.js examples/jsm/objects/MarchingCubes.js
   (three@0.185.1), which took them from Paul Bourke (paulbourke.net/geometry/polygonise), who had them from Cory
   Gene Bloyd. three.js is MIT-licensed: Copyright © 2010-2026 three.js authors. Permission is hereby granted, free
   of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
   to deal in the Software without restriction, including without limitation the rights to use, copy, modify,
   merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
   Software is furnished to do so, subject to the following conditions: The above copyright notice and this
   permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS
   PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND (full text: vendor/three/LICENSE.txt). */
(() => {
  'use strict';

  const COULOMB = 332.0637;          // e²/4πε₀ in kcal·Å/mol: V in kcal/mol per e for q in e and r in Å
  const FOV = 30;                    // vertical field of view, degrees; narrow so shapes are not distorted much
  const MIN_PICK_PX = 4;             // anything thinner than this many CSS px is still hit within it
  const NEAR_ATOM_PX = 12;           // atoms-only picking: an atom whose edge is this close to the pointer is picked
  const PICK_TIE = 0.02;             // Å: a bond this close in depth to an atom counts as level with it
  const PICK_STRIDE = 11, RIBBON = 2; // picking shapes: see rebuildModel
  const BALL_SCALE = 0.23, BALL_MIN = 0.18, BALL_BOND = 0.14, STICK = 0.18, WIRE_ATOM = 0.25, WIRE_GAP = 0.2;
  const SPIN_RATE = 0.5;             // rad/s about the view's vertical axis
  const INERTIA_TAU = 0.25, INERTIA_STOP = 0.05;
  const LABEL_EAGER = 150000;        // label × atom tests per frame below which labels stay on while moving
  const LABEL_FONT = '600 12px system-ui, "Segoe UI", sans-serif', LABEL_H = 14, LABEL_PAD = 1.5;   // CSS px
  const H_LABELS_MAX_HEAVY = 40;     // above this many heavy atoms, H labels show only on hover
  const NOTE_W = 600, NOTE_H = 20;   // CSS px: the bottom-left corner kept for the 'labels left out' note
  const HALO = { hover: [0x1e90ff, 0.5], select: [0xff8c00, 0.55], group: [0x16a34a, 0.42] };
  const MEASURE_COLOR = '#c2185b';
  const ELECTRON_COLOR = '#3f6fd8';   // lone-pair dots
  const PLAIN_SURFACE = '#c9d3e6';
  const STYLES = ['ballstick', 'spacefill', 'sticks', 'wireframe'];

  // Jmol default colours (PAL.java argbsCpk; D and T from its isotope table), only for pages without window.MolData
  const JMOL = {
    H: '#FFFFFF', D: '#FFFFC0', T: '#FFFFA0', He: '#D9FFFF', Li: '#CC80FF', Be: '#C2FF00', B: '#FFB5B5', C: '#909090',
    N: '#3050F8', O: '#FF0D0D', F: '#90E050', Ne: '#B3E3F5', Na: '#AB5CF2', Mg: '#8AFF00', Al: '#BFA6A6',
    Si: '#F0C8A0', P: '#FF8000', S: '#FFFF30', Cl: '#1FF01F', Ar: '#80D1E3', K: '#8F40D4', Ca: '#3DFF00',
    Ti: '#BFC2C7', Cr: '#8A99C7', Mn: '#9C7AC7', Fe: '#E06633', Co: '#F090A0', Ni: '#50D050', Cu: '#C88033',
    Zn: '#7D80B0', Ga: '#C28F8F', Ge: '#668F8F', As: '#BD80E3', Se: '#FFA100', Br: '#A62929', Kr: '#5CB8D1',
    Rb: '#702EB0', Sr: '#00FF00', Pd: '#006985', Ag: '#C0C0C0', Sn: '#668080', Sb: '#9E63B5', Te: '#D47A00',
    I: '#940094', Xe: '#429EB0', Cs: '#57178F', Ba: '#00C900', Pt: '#D0D0E0', Au: '#FFD123', Hg: '#B8B8D0',
    Pb: '#575961', Bi: '#9E4FB5',
  };
  const UNKNOWN_COLOR = '#FF1493';

  function elementInfo(a) {
    const MD = window.MolData;
    if (!MD || !MD.element) return null;
    const sym = a.el === 'D' || a.el === 'T' ? 'H' : a.el;
    return MD.element(sym) || (a.z ? MD.element(a.z) : null);
  }
  function atomColor(a) {
    if (a.color != null && a.color !== '') return a.color;
    if (a.el === 'D' || a.el === 'T') return JMOL[a.el];
    const e = elementInfo(a);
    return (e && e.color) || JMOL[a.el] || UNKNOWN_COLOR;
  }
  function atomVdw(a) {
    if (Number.isFinite(a.rvdw) && a.rvdw > 0) return a.rvdw;
    const e = elementInfo(a);
    return e && Number.isFinite(e.rvdw) && e.rvdw > 0 ? e.rvdw : null;
  }

  /* ---------- pure helpers (no THREE; exported as Mol3DView._internals for the test page) ---------- */

  // let the page breathe during long work: a message round-trip (setTimeout(0) is clamped, and throttled to ~1 s in
  // background tabs)
  const tick = typeof MessageChannel === 'function'
    ? () => new Promise(r => { const ch = new MessageChannel(); ch.port1.onmessage = () => { ch.port1.close(); r(); }; ch.port2.postMessage(0); })
    : () => new Promise(r => setTimeout(r, 0));

  // eigen-decomposition of a symmetric 3×3 matrix (cyclic Jacobi) → [{value, vector}] sorted by value, largest first
  function eigenSym3(m) {
    const A = m.map(r => r.slice()), V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let sweep = 0; sweep < 60; sweep++) {
      if (A[0][1] ** 2 + A[0][2] ** 2 + A[1][2] ** 2 < 1e-22) break;
      for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
        if (Math.abs(A[p][q]) < 1e-300) continue;
        const th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const t = (th >= 0 ? 1 : -1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < 3; k++) { const x = A[k][p], y = A[k][q]; A[k][p] = c * x - s * y; A[k][q] = s * x + c * y; }
        for (let k = 0; k < 3; k++) { const x = A[p][k], y = A[q][k]; A[p][k] = c * x - s * y; A[q][k] = s * x + c * y; }
        for (let k = 0; k < 3; k++) { const x = V[k][p], y = V[k][q]; V[k][p] = c * x - s * y; V[k][q] = s * x + c * y; }
      }
    }
    return [0, 1, 2].map(j => ({ value: A[j][j], vector: [V[0][j], V[1][j], V[2][j]] })).sort((a, b) => b.value - a.value);
  }

  // principal axes of a point set: [e1, e2, e3], largest spread first, always right-handed (e3 = e1 × e2), so
  // turning the molecule onto them is a proper rotation and never a mirror image
  function principalAxes(X, n) {
    if (!n) return { center: [0, 0, 0], axes: [[1, 0, 0], [0, 1, 0], [0, 0, 1]], spread: [0, 0, 0] };
    let cx = 0, cy = 0, cz = 0;
    for (let i = 0; i < n; i++) { cx += X[3 * i]; cy += X[3 * i + 1]; cz += X[3 * i + 2]; }
    cx /= n; cy /= n; cz /= n;
    const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < n; i++) {
      const d = [X[3 * i] - cx, X[3 * i + 1] - cy, X[3 * i + 2] - cz];
      for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) C[r][c] += d[r] * d[c];
    }
    const e = eigenSym3(C), e1 = e[0].vector, e2 = e[1].vector;
    const e3 = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
    const l3 = Math.hypot(e3[0], e3[1], e3[2]) || 1;
    return { center: [cx, cy, cz], axes: [e1, e2, e3.map(v => v / l3)], spread: e.map(x => x.value / Math.max(n, 1)) };
  }

  /* Surface of the union of spheres (centres X, radii R) by marching cubes on the signed distance
     f(p) = min_k(|p − c_k| − r_k). f is exact at every grid point within r + 2h of an atom and +3 elsewhere, which
     is enough: a grid edge crossing zero has one end inside a sphere, so the other end is within h of it.
     Returns {positions, normals, index, owner (atom whose sphere each vertex lies on), h, dims} or null if cancelled. */
  async function buildSurface(X, R, o = {}) {
    const N = R.length, cancelled = o.cancelled || (() => false);
    if (!N) throw new Error('no atoms, so no surface');
    let h = o.h || (N <= 100 ? 0.25 : N <= 500 ? 0.3 : 0.35);   // finer where it costs little
    const maxPoints = o.maxPoints || 200 ** 3;
    let rmax = 0, x0 = Infinity, y0 = Infinity, z0 = Infinity, x1 = -Infinity, y1 = -Infinity, z1 = -Infinity;
    for (let a = 0; a < N; a++) {
      const x = X[3 * a], y = X[3 * a + 1], z = X[3 * a + 2];
      rmax = Math.max(rmax, R[a]);
      x0 = Math.min(x0, x); y0 = Math.min(y0, y); z0 = Math.min(z0, z);
      x1 = Math.max(x1, x); y1 = Math.max(y1, y); z1 = Math.max(z1, z);
    }
    let nx, ny, nz, ox, oy, oz;
    for (;;) {                                    // grid spacing grows until the grid fits the cap
      const pad = rmax + 3 * h;
      ox = x0 - pad; oy = y0 - pad; oz = z0 - pad;
      nx = Math.ceil((x1 + pad - ox) / h) + 1; ny = Math.ceil((y1 + pad - oy) / h) + 1; nz = Math.ceil((z1 + pad - oz) / h) + 1;
      if (nx * ny * nz <= maxPoints) break;
      h *= 1.05;
    }
    const nxy = nx * ny, total = nxy * nz, F = new Float32Array(total).fill(3), OWN = new Int32Array(total).fill(-1);
    const cand = new Int32Array(8);
    let last = performance.now();
    const pause = async () => {
      if (performance.now() - last < 25) return false;
      await tick(); last = performance.now();
      return cancelled();
    };
    for (let a = 0; a < N; a++) {                 // splat each sphere's signed distance into the cells near it
      const cx = X[3 * a], cy = X[3 * a + 1], cz = X[3 * a + 2], r = R[a], reach = r + 2 * h, reach2 = reach * reach;
      const i0 = Math.max(0, Math.floor((cx - reach - ox) / h)), i1 = Math.min(nx - 1, Math.ceil((cx + reach - ox) / h));
      const j0 = Math.max(0, Math.floor((cy - reach - oy) / h)), j1 = Math.min(ny - 1, Math.ceil((cy + reach - oy) / h));
      const k0 = Math.max(0, Math.floor((cz - reach - oz) / h)), k1 = Math.min(nz - 1, Math.ceil((cz + reach - oz) / h));
      for (let k = k0; k <= k1; k++) {
        const dz = oz + k * h - cz, dz2 = dz * dz;
        for (let j = j0; j <= j1; j++) {
          const dy = oy + j * h - cy, dyz2 = dy * dy + dz2;
          if (dyz2 > reach2) continue;
          let g = i0 + nx * (j + ny * k);
          for (let i = i0; i <= i1; i++, g++) {
            const dx = ox + i * h - cx, d2 = dx * dx + dyz2;
            if (d2 > reach2) continue;
            const f = Math.sqrt(d2) - r;
            if (f < F[g]) { F[g] = f; OWN[g] = a; }
          }
        }
      }
      if (await pause()) return null;
    }

    // growable output buffers
    let vCap = 1 << 15, vN = 0, P = new Float32Array(vCap * 3), NR = new Float32Array(vCap * 3), VO = new Int32Array(vCap);
    let tCap = 1 << 16, tN = 0, T = new Uint32Array(tCap * 3);
    const growV = () => {
      vCap *= 2;
      const p = new Float32Array(vCap * 3); p.set(P); P = p;
      const n = new Float32Array(vCap * 3); n.set(NR); NR = n;
      const w = new Int32Array(vCap); w.set(VO); VO = w;
    };
    const vertex = (g0, g1, f0, f1) => {          // the zero crossing on the grid edge g0 → g1 (g0 the lower index)
      if (vN === vCap) growV();
      const t = f0 / (f0 - f1);
      const i = g0 % nx, j = ((g0 - i) / nx) % ny, k = (g0 - i - nx * j) / nxy;
      const d = g1 - g0, ax = d === 1 ? 1 : 0, ay = d === nx ? 1 : 0, az = d === nxy ? 1 : 0;
      const x0 = ox + i * h, y0 = oy + j * h, z0 = oz + k * h;
      // Interpolating f linearly along the edge is off by up to ~0.06 Å where spheres meet (f has a kink there), so
      // the zero of the exact f is found on the edge instead (Illinois false position). The vertex stays on its edge,
      // so the table's triangles cannot fold. f is evaluated over candidate atoms: the nearest atoms of the edge's
      // end points and of their six grid neighbours (so a third sphere meeting there is seen too).
      let nc = 0;
      for (const g of [g0, g1]) {
        for (const o of [0, 1, -1, nx, -nx, nxy, -nxy]) {
          const gg = g + o;
          if (gg < 0 || gg >= total) continue;
          const a = OWN[gg];
          if (a < 0) continue;
          let seen = false;
          for (let c = 0; c < nc; c++) if (cand[c] === a) { seen = true; break; }
          if (!seen && nc < cand.length) cand[nc++] = a;
        }
      }
      let best = -1;
      const fAt = s => {                         // exact f at fraction s along the edge; sets best
        const px = x0 + s * ax * h, py = y0 + s * ay * h, pz = z0 + s * az * h;
        let m = Infinity;
        for (let c = 0; c < nc; c++) {
          const a = cand[c], ex = px - X[3 * a], ey = py - X[3 * a + 1], ez = pz - X[3 * a + 2];
          const f = Math.sqrt(ex * ex + ey * ey + ez * ez) - R[a];
          if (f < m) { m = f; best = a; }
        }
        return m;
      };
      let lo = 0, hi = 1, flo = f0, fhi = f1, s = t, side = 0;
      for (let it = 0; it < 30; it++) {
        const fs = fAt(s);
        if (Math.abs(fs) < 1e-6 || hi - lo < 1e-7) break;
        if ((fs < 0) === (flo < 0)) { lo = s; flo = fs; if (side === -1) fhi /= 2; side = -1; }
        else { hi = s; fhi = fs; if (side === 1) flo /= 2; side = 1; }
        s = lo + (hi - lo) * flo / (flo - fhi);
      }
      fAt(s);
      const px = x0 + s * ax * h, py = y0 + s * ay * h, pz = z0 + s * az * h;
      const q = 3 * vN;
      P[q] = px; P[q + 1] = py; P[q + 2] = pz;
      if (best >= 0) {
        const ex = px - X[3 * best], ey = py - X[3 * best + 1], ez = pz - X[3 * best + 2], l = Math.sqrt(ex * ex + ey * ey + ez * ez) || 1;
        NR[q] = ex / l; NR[q + 1] = ey / l; NR[q + 2] = ez / l;
      }
      VO[vN] = best;
      return vN++;
    };
    // Triangles keep the table's winding, which is the same for every cube, so neighbours stay consistent and the
    // mesh closed; which way that winding faces is decided once, at the end, from the enclosed volume's sign.
    const tri = (a, b, c) => {
      if (a === b || b === c || a === c) return;
      if (tN === tCap) { tCap *= 2; const t = new Uint32Array(tCap * 3); t.set(T); T = t; }
      const q = 3 * tN++;
      T[q] = a; T[q + 1] = b; T[q + 2] = c;
    };

    // Edge caches, so neighbouring cubes share vertices: X and Y edges on the cube layer's bottom (A) and top (B)
    // planes, Z edges between them; indexed by the edge's lower grid point within its plane.
    let exA = new Int32Array(nxy).fill(-1), eyA = new Int32Array(nxy).fill(-1);
    let exB = new Int32Array(nxy), eyB = new Int32Array(nxy);
    const ez = new Int32Array(nxy), ev = new Int32Array(12);
    for (let k = 0; k < nz - 1; k++) {
      exB.fill(-1); eyB.fill(-1); ez.fill(-1);
      for (let j = 0; j < ny - 1; j++) {
        for (let i = 0; i < nx - 1; i++) {
          const l = i + nx * j, g = l + nxy * k;
          const f0 = F[g], f1 = F[g + 1], f2 = F[g + 1 + nx], f3 = F[g + nx];
          const f4 = F[g + nxy], f5 = F[g + 1 + nxy], f6 = F[g + 1 + nx + nxy], f7 = F[g + nx + nxy];
          let ci = 0;                             // Bourke's corner order: bit set = inside
          if (f0 < 0) ci |= 1; if (f1 < 0) ci |= 2; if (f2 < 0) ci |= 4; if (f3 < 0) ci |= 8;
          if (f4 < 0) ci |= 16; if (f5 < 0) ci |= 32; if (f6 < 0) ci |= 64; if (f7 < 0) ci |= 128;
          if (ci === 0 || ci === 255) continue;
          const bits = EDGE_TABLE[ci];
          if (bits & 1) ev[0] = exA[l] >= 0 ? exA[l] : (exA[l] = vertex(g, g + 1, f0, f1));
          if (bits & 2) ev[1] = eyA[l + 1] >= 0 ? eyA[l + 1] : (eyA[l + 1] = vertex(g + 1, g + 1 + nx, f1, f2));
          if (bits & 4) ev[2] = exA[l + nx] >= 0 ? exA[l + nx] : (exA[l + nx] = vertex(g + nx, g + nx + 1, f3, f2));
          if (bits & 8) ev[3] = eyA[l] >= 0 ? eyA[l] : (eyA[l] = vertex(g, g + nx, f0, f3));
          if (bits & 16) ev[4] = exB[l] >= 0 ? exB[l] : (exB[l] = vertex(g + nxy, g + nxy + 1, f4, f5));
          if (bits & 32) ev[5] = eyB[l + 1] >= 0 ? eyB[l + 1] : (eyB[l + 1] = vertex(g + 1 + nxy, g + 1 + nx + nxy, f5, f6));
          if (bits & 64) ev[6] = exB[l + nx] >= 0 ? exB[l + nx] : (exB[l + nx] = vertex(g + nx + nxy, g + nx + 1 + nxy, f7, f6));
          if (bits & 128) ev[7] = eyB[l] >= 0 ? eyB[l] : (eyB[l] = vertex(g + nxy, g + nx + nxy, f4, f7));
          if (bits & 256) ev[8] = ez[l] >= 0 ? ez[l] : (ez[l] = vertex(g, g + nxy, f0, f4));
          if (bits & 512) ev[9] = ez[l + 1] >= 0 ? ez[l + 1] : (ez[l + 1] = vertex(g + 1, g + 1 + nxy, f1, f5));
          if (bits & 1024) ev[10] = ez[l + 1 + nx] >= 0 ? ez[l + 1 + nx] : (ez[l + 1 + nx] = vertex(g + 1 + nx, g + 1 + nx + nxy, f2, f6));
          if (bits & 2048) ev[11] = ez[l + nx] >= 0 ? ez[l + nx] : (ez[l + nx] = vertex(g + nx, g + nx + nxy, f3, f7));
          const base = ci << 4;
          for (let t = 0; TRI_TABLE[base + t] !== -1; t += 3) tri(ev[TRI_TABLE[base + t]], ev[TRI_TABLE[base + t + 1]], ev[TRI_TABLE[base + t + 2]]);
        }
      }
      let s = exA; exA = exB; exB = s;            // this layer's top plane is the next layer's bottom plane
      s = eyA; eyA = eyB; eyB = s;
      if (await pause()) return null;
    }
    // outward winding (counter-clockwise seen from outside) encloses a positive volume: flip everything if not
    let vol = 0;
    for (let t = 0; t < 3 * tN; t += 3) {
      const a = 3 * T[t], b = 3 * T[t + 1], c = 3 * T[t + 2];
      vol += P[a] * (P[b + 1] * P[c + 2] - P[b + 2] * P[c + 1]) - P[a + 1] * (P[b] * P[c + 2] - P[b + 2] * P[c]) + P[a + 2] * (P[b] * P[c + 1] - P[b + 1] * P[c]);
    }
    if (vol < 0) for (let t = 0; t < 3 * tN; t += 3) { const x = T[t + 1]; T[t + 1] = T[t + 2]; T[t + 2] = x; }
    return { positions: P.slice(0, 3 * vN), normals: NR.slice(0, 3 * vN), index: T.slice(0, 3 * tN), owner: VO.slice(0, vN), h, dims: [nx, ny, nz], volume: Math.abs(vol) / 6 };
  }

  // Coulomb potential at points P (Float32Array/Float64Array 3V) from charges Q at X: kcal/mol per e
  async function coulomb(P, X, Q, cancelled = () => false) {
    const V = new Float64Array(P.length / 3);
    const idx = [];                                      // only atoms that carry a charge
    for (let a = 0; a < Q.length; a++) if (Q[a] !== 0) idx.push(a);
    const n = idx.length, cx = new Float64Array(n), cy = new Float64Array(n), cz = new Float64Array(n), cq = new Float64Array(n);
    idx.forEach((a, j) => { cx[j] = X[3 * a]; cy[j] = X[3 * a + 1]; cz[j] = X[3 * a + 2]; cq[j] = Q[a]; });
    let last = performance.now();
    for (let v = 0; v < V.length; v++) {
      const px = P[3 * v], py = P[3 * v + 1], pz = P[3 * v + 2];
      let s = 0;
      for (let j = 0; j < n; j++) {
        const dx = px - cx[j], dy = py - cy[j], dz = pz - cz[j];
        s += cq[j] / Math.sqrt(dx * dx + dy * dy + dz * dz);
      }
      V[v] = COULOMB * s;
      if ((v & 1023) === 1023 && performance.now() - last > 25) {
        await tick(); last = performance.now();
        if (cancelled()) return null;
      }
    }
    return V;
  }

  // the value below which a fraction p of the (sorted copy of) values lie
  function percentile(values, p) {
    if (!values.length) return 0;
    const s = Float64Array.from(values).sort();
    return s[Math.min(s.length - 1, Math.max(0, Math.round(p * (s.length - 1))))];
  }

  // surface colour for the scaled potential s ∈ [−1, 1], as sRGB 0..1: negative white → red, positive white → blue
  function espRGB(s) {
    s = Math.max(-1, Math.min(1, s || 0));
    return s < 0 ? [1, 1 + s, 1 + s] : [1 - s, 1 - s, 1];
  }
  function espColor(V, clamp) {
    const c = espRGB(clamp > 0 ? V / clamp : 0);
    return '#' + c.map(x => Math.round(x * 255).toString(16).padStart(2, '0')).join('');
  }

  /* ---------- the view ---------- */

  let spareRenderer = null;          // one WebGL context per page: a disposed view hands its renderer to the next
  let styleInjected = false;
  function injectStyle() {
    if (styleInjected) return;
    styleInjected = true;
    const st = document.createElement('style');
    st.textContent = '.mv-3d-canvas{outline:none}.mv-3d-canvas:focus-visible{outline:2px solid #2f6fd8;outline-offset:-2px}';
    document.head.appendChild(st);
  }

  function create(container, opts = {}) {
    const THREE = window.THREE;
    if (!THREE || !THREE.WebGLRenderer) throw new Error('Mol3DView: three.js (window.THREE) is not loaded');
    if (!container || !container.appendChild) throw new Error('Mol3DView: no container element');
    injectStyle();
    const call = (name, ...args) => {
      const f = opts[name];
      if (typeof f !== 'function') return;
      try { f(...args); } catch (e) { console.error('Mol3DView ' + name + ' callback:', e); }
    };

    /* DOM: a wrapper holding the WebGL canvas and a 2D canvas for labels on top of it */
    const wrap = document.createElement('div');
    wrap.className = 'mv-3d';
    wrap.style.cssText = 'position:relative;width:100%;height:100%;overflow:hidden;';
    let renderer = spareRenderer;
    spareRenderer = null;
    if (renderer && renderer.getContext().isContextLost()) { renderer.dispose(); renderer = null; }
    if (!renderer) {
      // alpha: a PNG can be exported with a transparent background (on screen the scene background is opaque, as before)
      renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true, alpha: true });
      renderer.outputColorSpace = THREE.SRGBColorSpace;
      renderer.toneMapping = THREE.NoToneMapping;          // keep the element colours true
    }
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    renderer.setSize(1, 1, false);
    const canvas = renderer.domElement;
    canvas.className = 'mv-3d-canvas';
    canvas.tabIndex = 0;
    canvas.setAttribute('role', 'img');
    canvas.setAttribute('aria-label', '3D molecule. Drag to rotate, Shift+drag to roll, right-drag or Ctrl+drag to pan, ' +
      'wheel to zoom; keys: arrows rotate, + and − zoom, 0 resets, Space spins.');
    canvas.style.cssText = 'position:absolute;left:0;top:0;width:100%;height:100%;display:block;touch-action:none;';
    const overlay = document.createElement('canvas');
    overlay.className = 'mv-3d-labels';
    overlay.style.cssText = 'position:absolute;left:0;top:0;width:100%;height:100%;pointer-events:none;';
    wrap.append(canvas, overlay);
    container.appendChild(wrap);
    const octx = overlay.getContext('2d');

    /* scene: the camera never turns; the molecule group does (so the lights stay fixed relative to the viewer) */
    const bgColor = opts.background || '#ffffff';
    const scene = new THREE.Scene();
    scene.background = new THREE.Color(bgColor);
    const darkBg = (() => { const c = scene.background; return 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b < 0.2; })();
    const persp = new THREE.PerspectiveCamera(FOV, 1, 0.1, 1000);
    const ortho = new THREE.OrthographicCamera(-1, 1, 1, -1, 0.01, 1000);
    let camera = persp;
    scene.add(new THREE.AmbientLight(0xffffff, 0.5));
    const sun = new THREE.DirectionalLight(0xffffff, 2.0);
    sun.position.set(2, 3, 5);
    scene.add(sun);
    const molGroup = new THREE.Group();          // trackball rotation and pan
    const inner = new THREE.Group();             // position = −centroid, so rotation pivots on the centroid
    molGroup.add(inner);
    scene.add(molGroup);

    // shared geometries and materials (made again, never disposed, after a context loss: see onRestored)
    const geo = {}, mat = {};
    function makeResources() {
      Object.assign(geo, {
        sphereHi: new THREE.SphereGeometry(1, 32, 16), sphereLo: new THREE.SphereGeometry(1, 20, 12),
        cylOpen: new THREE.CylinderGeometry(1, 1, 1, 16, 1, true), cylClosed: new THREE.CylinderGeometry(1, 1, 1, 12, 1, false),
        haloSphere: new THREE.SphereGeometry(1, 32, 16), haloCyl: new THREE.CylinderGeometry(1, 1, 1, 16, 1, true),
      });
      Object.assign(mat, {
        solid: new THREE.MeshStandardMaterial({ roughness: 0.45, metalness: 0 }),
        line: new THREE.LineBasicMaterial({ vertexColors: true }),
        dash: new THREE.LineDashedMaterial({ vertexColors: true, dashSize: 0.15, gapSize: 0.1 }),
        measure: new THREE.LineDashedMaterial({ color: MEASURE_COLOR, dashSize: 0.2, gapSize: 0.12, depthTest: false, transparent: true }),
        halo: {},
        surfPre: new THREE.MeshBasicMaterial({ colorWrite: false, depthWrite: true, transparent: true }),
        electron: new THREE.MeshStandardMaterial({ color: ELECTRON_COLOR, roughness: 0.4, metalness: 0, emissive: ELECTRON_COLOR, emissiveIntensity: 0.35 }),
      });
      for (const k in HALO) {
        mat.halo[k] = new THREE.MeshBasicMaterial({ color: HALO[k][0], transparent: true, opacity: HALO[k][1], depthWrite: false, side: THREE.BackSide });
      }
    }
    makeResources();

    // scratch objects, reused so the hot paths do not allocate
    const _m = new THREE.Matrix4(), _q = new THREE.Quaternion(), _qi = new THREE.Quaternion(), _s = new THREE.Vector3();
    const _p = new THREE.Vector3(), _a = new THREE.Vector3(), _b = new THREE.Vector3(), _d = new THREE.Vector3();
    const _c = new THREE.Color(), _Y = new THREE.Vector3(0, 1, 0);
    const _ray = new THREE.Ray(), _inv = new THREE.Matrix4();
    const raycaster = new THREE.Raycaster(), _ndc = new THREE.Vector2();

    /* state */
    let mol = null;              // prepared molecule (prepare())
    let X = null;                // coordinates, Float64Array 3N
    let pca = null;              // principal axes of the current coordinates
    let style = 'ballstick', showH = true, labels = null, electrons = null, measures = [], spin = false, orthoOn = false, atomsOnly = false;
    let needFrame = false, framedWithoutSize = false, userMoved = false;
    let atomVis = null, bondVis = null, atomR = null, atomPickR = null, visDeg = null;
    let pickSegs = new Float64Array(0), pickSegBond = new Int32Array(0), nPickSegs = 0;
    let bondHaloR = null;
    const model = { atoms: null, atomOf: null, bonds: null, bondOf: null, extras: null, extraOf: null, lines: null, dashes: null, measure: null };
    const halos = {};            // kind → {sel, sphere: InstancedMesh, cyl: InstancedMesh}
    const view = { dist: 10, R: 5 };
    let cssW = 0, cssH = 0;
    let surf = null;             // {cfg, geom, mesh, pre, result, X}
    let surfCfg = null, surfGen = 0;
    let lost = false, disposed = false;
    let raf = 0, lastT = 0, inFrame = false, moving = false, settleTimer = 0;
    let unseen = false;          // a frame found the view not visible (frame()); cleared when it is shown again
    let inertia = null;          // THREE.Vector3 angular velocity (rad/s, world axes)
    const pointers = new Map();
    let drag = null;
    const hover = { x: 0, y: 0, inside: false, pending: false, key: '' };
    let labelVisCache = null;
    let labelStats = { shown: 0, overlapping: 0, hOnHover: 0, hidden: 0 };   // the last labels drawn (drawOverlay)
    let lastLabels = [];         // the labels on screen: [{index, text, box:[x0, y0, x1, y1] CSS px}] (for the tests)
    let lastMeasures = [];       // the measurement labels on screen: [{x, y, w, h, text}] (for the tests)
    const stats = { renderMs: 0, renderAvgMs: 0, frames: 0, pickMs: 0, pickAvgMs: 0, picks: 0, surfaceMs: 0, labelMs: 0 };

    /* ---------- molecule ---------- */

    function prepare(m) {
      if (!m || !Array.isArray(m.atoms) || !Array.isArray(m.bonds)) throw new Error('Mol3DView.setMolecule: expected {atoms, bonds, rings}');
      const N = m.atoms.length, M = m.bonds.length;
      const el = [], z = new Int32Array(N), rvdw = new Float64Array(N), hasVdw = new Uint8Array(N), colors = [];
      const adj = Array.from({ length: N }, () => []);
      m.atoms.forEach((a, i) => {
        if (!a || typeof a.el !== 'string') throw new Error('Mol3DView.setMolecule: atom ' + i + ' has no element');
        el.push(a.el);
        z[i] = a.z | 0;
        const r = atomVdw(a);
        hasVdw[i] = r ? 1 : 0;
        rvdw[i] = r || 2.0;                        // drawing fallback only; the surface refuses atoms without one
        colors.push(new THREE.Color(atomColor(a)));  // sRGB hex → linear, as colour management expects
      });
      const bonds = m.bonds.map((b, k) => {
        const a = b.a, c = b.b;
        if (!(a >= 0 && a < N && c >= 0 && c < N && a !== c)) throw new Error('Mol3DView.setMolecule: bond ' + k + ' has bad atoms');
        adj[a].push({ atom: c, bond: k }); adj[c].push({ atom: a, bond: k });
        return { a, b: c, order: b.order === 2 || b.order === 3 ? b.order : 1, aromatic: !!b.aromatic };
      });
      const rings = (m.rings || []).map(r => {
        if (!Array.isArray(r) || r.some(i => !(i >= 0 && i < N))) throw new Error('Mol3DView.setMolecule: bad ring');
        return r.slice();
      });
      // H that "Show H" hides: the caller's isH, else a hydrogen with one bond to a non-hydrogen
      const isHyd = i => z[i] === 1 || el[i] === 'H' || el[i] === 'D' || el[i] === 'T';
      const hideable = new Uint8Array(N);
      m.atoms.forEach((a, i) => {
        hideable[i] = typeof a.isH === 'boolean' ? (a.isH ? 1 : 0)
          : (isHyd(i) && adj[i].length === 1 && !isHyd(adj[i][0].atom) ? 1 : 0);
      });
      // the ring that places a bond's inner line: for an aromatic bond the smallest fully aromatic ring holding it,
      // otherwise the smallest ring holding it (fused rings: always the same one, the first smallest)
      const bondKey = (a, c) => a < c ? a * N + c : c * N + a;
      const bondIndex = new Map(bonds.map((b, k) => [bondKey(b.a, b.b), k]));
      const ringBonds = rings.map(r => {
        const ks = [];
        for (let i = 0; i < r.length; i++) {
          const k = bondIndex.get(bondKey(r[i], r[(i + 1) % r.length]));
          if (k !== undefined) ks.push(k);
        }
        return ks;
      });
      const ringAromatic = ringBonds.map((ks, ri) => ks.length === rings[ri].length && ks.every(k => bonds[k].aromatic));
      const bondRing = new Int32Array(M).fill(-1);
      bonds.forEach((b, k) => {
        let best = -1;
        rings.forEach((r, ri) => {
          if (!r.includes(b.a) || !r.includes(b.b)) return;
          const better = best < 0 ||
            (b.aromatic && ringAromatic[ri] && !ringAromatic[best]) ||
            ((!b.aromatic || ringAromatic[ri] === ringAromatic[best]) && r.length < rings[best].length);
          if (better) best = ri;
        });
        bondRing[k] = best;
      });
      return { N, M, el, z, rvdw, hasVdw, colors, adj, bonds, rings, hideable, bondRing };
    }

    function setMolecule(m) {
      const prepared = prepare(m);
      mol = prepared;
      X = null; pca = null;
      labels = null; measures = []; labelVisCache = null;
      for (const k in halos) halos[k].sel = null;
      surfCfg = null; surfGen++; removeSurface();
      clearModel(); updateAllHalos(); rebuildMeasures();
      needFrame = true;
      hover.key = '';
      hover.pending = hover.inside;                      // tell the caller again what the pointer is over
      requestRender();
    }

    function setCoords(coords, o = {}) {
      if (!mol) throw new Error('Mol3DView.setCoords: call setMolecule first');
      if (!coords || coords.length !== 3 * mol.N) throw new Error('Mol3DView.setCoords: expected ' + (3 * mol.N) + ' numbers, got ' + (coords ? coords.length : 0));
      for (let i = 0; i < coords.length; i++) if (!Number.isFinite(coords[i])) throw new Error('Mol3DView.setCoords: coordinate ' + i + ' is not a number');
      X = Float64Array.from(coords);
      pca = principalAxes(X, mol.N);
      inner.position.set(-pca.center[0], -pca.center[1], -pca.center[2]);
      view.R = boundRadius();
      rebuildModel();
      if (needFrame || o.keepView === false) { resetView(); needFrame = false; }
      labelVisCache = null;
      // a surface belongs to the coordinates it was built on: drop it now, rebuild it for the new ones
      if (surfCfg) {
        removeSurface();
        const gen = ++surfGen;
        makeSurface(surfCfg, gen).then(r => { if (gen === surfGen) call('onSurface', r, null); },
          e => { if (gen === surfGen) call('onSurface', null, e); });
      }
      updateCamera();
      hover.pending = hover.inside;
      requestRender();
    }

    function boundRadius() {                       // centroid → farthest vdW sphere surface, over all atoms
      let R = 0;
      const c = pca.center;
      for (let i = 0; i < mol.N; i++) R = Math.max(R, Math.hypot(X[3 * i] - c[0], X[3 * i + 1] - c[1], X[3 * i + 2] - c[2]) + mol.rvdw[i]);
      return Math.max(R, 1.5);
    }

    /* ---------- building the picture ---------- */

    function drawRadius(i) {
      const ball = Math.max(BALL_MIN, BALL_SCALE * mol.rvdw[i]);
      switch (style) {
        case 'spacefill': return mol.rvdw[i];
        case 'sticks': return visDeg[i] ? STICK : ball;
        case 'wireframe': return visDeg[i] ? 0 : ball;
        default: return ball;
      }
    }

    function disposeObj(o) {
      if (!o) return;
      if (o.parent) o.parent.remove(o);
      if (o.isInstancedMesh) o.dispose();                // frees the instance buffers
      else if (o.geometry && !Object.values(geo).includes(o.geometry)) o.geometry.dispose();
    }
    function clearModel() {
      for (const k of ['atoms', 'bonds', 'extras', 'lines', 'dashes', 'lonePairs']) { disposeObj(model[k]); model[k] = null; }
      model.atomOf = model.bondOf = model.extraOf = null;
      nPickSegs = 0;
      atomVis = bondVis = atomR = atomPickR = visDeg = bondHaloR = null;
    }

    // unit vector n perpendicular to bond k, in the plane its extra lines go: toward the ring centre for a ring bond,
    // else in the plane of a neighbouring bond (the σ framework, as a structure drawing puts a double bond's second line)
    function bondNormal(k, out) {
      const b = mol.bonds[k], A = 3 * b.a, B = 3 * b.b;
      _d.set(X[B] - X[A], X[B + 1] - X[A + 1], X[B + 2] - X[A + 2]).normalize();
      const perp = (vx, vy, vz) => {
        const t = vx * _d.x + vy * _d.y + vz * _d.z;
        out.set(vx - t * _d.x, vy - t * _d.y, vz - t * _d.z);
        return out.length();
      };
      const ri = mol.bondRing[k];
      if (ri >= 0) {                                   // toward the ring's centre
        const r = mol.rings[ri];
        let cx = 0, cy = 0, cz = 0;
        for (const i of r) { cx += X[3 * i]; cy += X[3 * i + 1]; cz += X[3 * i + 2]; }
        if (perp(cx / r.length - X[A], cy / r.length - X[A + 1], cz / r.length - X[A + 2]) > 1e-3) return out.normalize();
      }
      for (const [from, other] of [[b.a, b.b], [b.b, b.a]]) {   // in the plane of a neighbour (the σ framework)
        for (const nb of mol.adj[from]) {
          if (nb.atom === other) continue;
          const C = 3 * nb.atom, F = 3 * from;
          if (perp(X[C] - X[F], X[C + 1] - X[F + 1], X[C + 2] - X[F + 2]) > 0.1) return out.normalize();
        }
      }
      // nothing to go by (O=C=O, a linear chain): perpendicular to the flattest principal axis, so the lines show
      // side by side in the default view; else perpendicular to the coordinate axis least aligned with the bond
      const e3 = pca.axes[2];
      out.set(_d.y * e3[2] - _d.z * e3[1], _d.z * e3[0] - _d.x * e3[2], _d.x * e3[1] - _d.y * e3[0]);
      if (out.length() > 0.1) return out.normalize();
      const ax = Math.abs(_d.x), ay = Math.abs(_d.y), az = Math.abs(_d.z);
      const e = ax <= ay && ax <= az ? [1, 0, 0] : ay <= az ? [0, 1, 0] : [0, 0, 1];
      out.set(_d.y * e[2] - _d.z * e[1], _d.z * e[0] - _d.x * e[2], _d.x * e[1] - _d.y * e[0]);
      return out.normalize();
    }

    // the drawn pieces of bond k in the current style: [t0, t1, offset along n (Å), radius, closed?]
    function dashes(t0, t1, off, r, n = 3) {
      const len = (t1 - t0) / (n + 0.6 * (n - 1)), out = [];
      for (let i = 0; i < n; i++) { const s = t0 + i * 1.6 * len; out.push([s, s + len, off, r, true]); }
      return out;
    }
    function bondKind(k) {
      const b = mol.bonds[k];
      if (b.aromatic && mol.bondRing[k] >= 0) return 'ar';
      return b.order === 2 ? 'd' : b.order === 3 ? 't' : 's';
    }
    function bondParts(k) {
      const kind = bondKind(k);
      if (style === 'ballstick') {
        const r = BALL_BOND;
        if (kind === 'd') return [[0, 1, 0.8 * r, 0.5 * r, false], [0, 1, -0.8 * r, 0.5 * r, false]];
        if (kind === 't') return [[0, 1, 0, 0.4 * r, false], [0, 1, 1.1 * r, 0.4 * r, false], [0, 1, -1.1 * r, 0.4 * r, false]];
        if (kind === 'ar') return [[0, 1, 0, r, false], ...dashes(0.22, 0.78, 2.2 * r, 0.45 * r)];
        return [[0, 1, 0, r, false]];
      }
      if (style === 'sticks') {                        // a full stick, plus thin rails beside it for the extra bonds
        const r = STICK, off = 1.9 * r, rr = 0.38 * r;
        if (kind === 'd') return [[0, 1, 0, r, false], [0.2, 0.8, off, rr, true]];
        if (kind === 't') return [[0, 1, 0, r, false], [0.2, 0.8, off, rr, true], [0.2, 0.8, -off, rr, true]];
        if (kind === 'ar') return [[0, 1, 0, r, false], ...dashes(0.2, 0.8, off, rr)];
        return [[0, 1, 0, r, false]];
      }
      return [];
    }

    /* ---------- lone pairs ----------
       Each lone pair is two small spheres just outside the atom, on a direction found by pushing the pairs away from
       the atom's bonds and from each other on a sphere around it (the VSEPR idea: electron pairs get as far apart as
       they can). The bonds are held fixed and only the lone pairs move. Every real bond counts, drawn or not: with
       the H atoms hidden, water's O–H bonds still decide where its two pairs go (above and below the H–O–H plane),
       and ammonia's one pair still sits on the axis opposite the three N–H. Repulsion alone cannot say which way round
       a terminal atom's pairs lie (a carbonyl O's two pairs could turn freely about the C=O axis), so for a terminal
       atom with two lone pairs — three electron domains, trigonal — they are put in the plane of the neighbour's other
       bonds, 120° from the bond, as in an sp² carbonyl oxygen. With no bonds at all (a bare ion) the pairs spread over
       the whole sphere. */
    const bondDir = (i, j) => new THREE.Vector3(X[3 * j] - X[3 * i], X[3 * j + 1] - X[3 * i + 1], X[3 * j + 2] - X[3 * i + 2]);
    function lonePairDirs(i, k) {
      const fixed = [];
      for (const nb of mol.adj[i]) {
        const d = bondDir(i, nb.atom);
        if (d.lengthSq() > 1e-12) fixed.push(d.normalize());
      }
      if (fixed.length === 1 && k === 2) {               // trigonal terminal atom: in the neighbour's plane
        const j = mol.adj[i][0].atom, b = fixed[0];
        let p = null;
        for (const nb of mol.adj[j]) {
          if (nb.atom === i) continue;
          const v = bondDir(j, nb.atom), t = v.dot(b);
          v.addScaledVector(b, -t);                        // the part of the neighbour's other bond square to b
          if (v.lengthSq() > 1e-6) { p = v.normalize(); break; }
        }
        if (p) return [b.clone().multiplyScalar(-0.5).addScaledVector(p, 0.866), b.clone().multiplyScalar(-0.5).addScaledVector(p, -0.866)];
      }
      /* start from evenly spread points (a Fibonacci sphere), taking the k farthest from the bonds */
      const SEED = 64, cand = [];
      for (let s = 0; s < SEED; s++) {
        const y = 1 - 2 * (s + 0.5) / SEED, r = Math.sqrt(Math.max(0, 1 - y * y)), th = Math.PI * (1 + Math.sqrt(5)) * s;
        const v = new THREE.Vector3(Math.cos(th) * r, y, Math.sin(th) * r);
        let worst = -1;
        for (const f of fixed) worst = Math.max(worst, v.dot(f));
        cand.push({ v, worst });
      }
      cand.sort((a, b) => a.worst - b.worst);
      const pts = [];
      for (const c of cand) {
        if (pts.length >= k) break;
        if (pts.every(p => p.dot(c.v) < 0.5)) pts.push(c.v.clone());
      }
      while (pts.length < k) pts.push(cand[pts.length % cand.length].v.clone());
      /* relax: every other pair and every bond pushes, the bonds do not move */
      const force = new THREE.Vector3(), tmp = new THREE.Vector3();
      for (let it = 0; it < 120; it++) {
        const step = 0.25 * (1 - it / 120);
        for (let a = 0; a < pts.length; a++) {
          force.set(0, 0, 0);
          const push = o => { tmp.subVectors(pts[a], o); const d2 = Math.max(0.02, tmp.lengthSq()); force.addScaledVector(tmp.normalize(), 1 / d2); };
          for (let b = 0; b < pts.length; b++) if (b !== a) push(pts[b]);
          for (const f of fixed) push(f);
          if (force.lengthSq() < 1e-12) continue;
          pts[a].addScaledVector(force.normalize(), step).normalize();
        }
      }
      return pts;
    }
    function rebuildElectrons() {
      disposeObj(model.lonePairs); model.lonePairs = null;
      if (!electrons || !mol || !X) return;
      const list = [], up = new THREE.Vector3(), side = new THREE.Vector3();
      for (let i = 0; i < mol.N; i++) {
        const k = Math.round(+electrons[i] || 0);
        if (k <= 0 || !atomVis[i]) continue;
        /* big enough to read next to the atom, and clear of its surface (in space-filling the spheres are large, so
           the dots go a little further out); the two dots of a pair just touch */
        const r = Math.max(atomR[i], 0.12), rad = Math.min(0.13, Math.max(0.075, 0.16 * r)), gap = r + 2.1 * rad;
        const sep = 2.2 * rad;
        /* the local plane: through the atom's first two bonds, or for a terminal atom through its bond and the
           neighbour's next bond. The two dots of a pair lie side by side IN that plane (square to the pair's
           direction), so a flat molecule seen face-on shows both dots of every in-plane pair, as a Lewis drawing does */
        const nbs = mol.adj[i].map(nb => nb.atom);
        up.set(0, 0, 0);
        if (nbs.length >= 2) up.crossVectors(bondDir(i, nbs[0]), bondDir(i, nbs[1]));
        else if (nbs.length === 1) {
          const j = nbs[0], other = mol.adj[j].find(nb => nb.atom !== i);
          if (other) up.crossVectors(bondDir(i, j), bondDir(j, other.atom));
        }
        for (const dir of lonePairDirs(i, k)) {
          side.set(0, 0, 0);
          if (up.lengthSq() > 1e-9) side.crossVectors(dir, up);
          if (side.lengthSq() < 1e-9 && nbs.length) {     // a pair along the plane's normal (a p-type pair): dots along a bond's shadow
            side.copy(bondDir(i, nbs[0])).addScaledVector(dir, -bondDir(i, nbs[0]).dot(dir));
          }
          if (side.lengthSq() < 1e-9) side.crossVectors(dir, Math.abs(dir.y) < 0.9 ? _Y : new THREE.Vector3(1, 0, 0));
          side.normalize().multiplyScalar(sep / 2);
          for (const sgn of [1, -1]) list.push([
            X[3 * i] + dir.x * gap + sgn * side.x, X[3 * i + 1] + dir.y * gap + sgn * side.y, X[3 * i + 2] + dir.z * gap + sgn * side.z, rad]);
        }
      }
      if (!list.length) return;
      const mesh = new THREE.InstancedMesh(geo.sphereLo, mat.electron, list.length);
      list.forEach((d, j) => { _m.compose(_p.set(d[0], d[1], d[2]), _qi.identity(), _s.setScalar(d[3])); mesh.setMatrixAt(j, _m); });
      mesh.renderOrder = 1;
      finishInstances(mesh);
      model.lonePairs = mesh;
      inner.add(mesh);
    }
    function setElectrons(arr) {
      if (arr != null && (!mol || !Array.isArray(arr) || arr.length !== mol.N)) throw new Error('Mol3DView.setElectrons: expected ' + (mol ? mol.N : 0) + ' counts');
      electrons = arr ? arr.map(n => Math.max(0, Math.round(+n || 0))) : null;
      rebuildElectrons();
      requestRender();
    }

    function rebuildModel() {
      clearModel();
      if (!mol || !X) return;
      const N = mol.N, M = mol.M;
      atomVis = new Uint8Array(N); bondVis = new Uint8Array(M); visDeg = new Int32Array(N);
      for (let i = 0; i < N; i++) atomVis[i] = showH || !mol.hideable[i] ? 1 : 0;
      mol.bonds.forEach((b, k) => {
        bondVis[k] = atomVis[b.a] && atomVis[b.b] ? 1 : 0;
        if (bondVis[k]) { visDeg[b.a]++; visDeg[b.b]++; }
      });
      atomR = new Float64Array(N); atomPickR = new Float64Array(N);
      for (let i = 0; i < N; i++) {
        atomR[i] = drawRadius(i);
        atomPickR[i] = style === 'wireframe' && visDeg[i] ? WIRE_ATOM : atomR[i];
      }

      // atoms: one instanced sphere mesh
      const drawn = [];
      for (let i = 0; i < N; i++) if (atomVis[i] && atomR[i] > 0) drawn.push(i);
      if (drawn.length) {
        const mesh = new THREE.InstancedMesh(N <= 200 ? geo.sphereHi : geo.sphereLo, mat.solid, drawn.length);
        const of = new Int32Array(drawn.length);
        drawn.forEach((i, n) => {
          _m.compose(_p.set(X[3 * i], X[3 * i + 1], X[3 * i + 2]), _qi.identity(), _s.setScalar(atomR[i]));
          mesh.setMatrixAt(n, _m); mesh.setColorAt(n, mol.colors[i]); of[n] = i;
        });
        finishInstances(mesh);
        model.atoms = mesh; model.atomOf = of;
        inner.add(mesh);
      }

      // bonds: open cylinders for the main pieces, closed ones for dashes and rails; each piece split at the bond's
      // midpoint so each half takes its atom's colour
      // Pieces are [bond, t0, t1, offset, radius, colour atom, n] with n the bond's offset direction (bondNormal).
      // Every drawn cylinder is also a picking cylinder [bond, t0, t1, offset, radius, n, closed], so the picker
      // sees exactly what is drawn (thin pieces are widened to MIN_PICK_PX there).
      const mains = [], extras = [], segs = [];
      bondHaloR = new Float64Array(M);
      const n = new THREE.Vector3();
      if (style === 'ballstick' || style === 'sticks') {
        for (let k = 0; k < M; k++) {
          if (!bondVis[k]) continue;
          const parts = bondParts(k), b = mol.bonds[k];
          if (parts.some(p => p[2] !== 0)) bondNormal(k, n); else n.set(0, 0, 0);
          const nn = [n.x, n.y, n.z];
          for (const [t0, t1, off, r, closed] of parts) {
            segs.push([k, t0, t1, off, r, nn, closed]);
            const list = closed ? extras : mains;
            if (t0 < 0.5 && t1 > 0.5) { list.push([k, t0, 0.5, off, r, b.a, nn]); list.push([k, 0.5, t1, off, r, b.b, nn]); }
            else list.push([k, t0, t1, off, r, (t0 + t1) / 2 < 0.5 ? b.a : b.b, nn]);
          }
          pickEnvelope(k, parts, nn, segs);
          bondHaloR[k] = style === 'sticks' ? 2.5 * STICK : 2.2 * BALL_BOND;
        }
      } else if (style === 'wireframe') {
        buildWire(n, segs);
      }
      model.bonds = cylinderMesh(mains, geo.cylOpen, 'bondOf');
      model.extras = cylinderMesh(extras, geo.cylClosed, 'extraOf');

      // picking shapes in molecule coordinates: [ax, ay, az, bx, by, bz, r, kind, nx, ny, nz], kind 0 an open
      // cylinder, 1 a closed one, RIBBON a flat strip of half-width r along n
      nPickSegs = segs.length;
      pickSegs = new Float64Array(PICK_STRIDE * nPickSegs); pickSegBond = new Int32Array(nPickSegs);
      segs.forEach(([k, t0, t1, off, r, nn, kind], j) => {
        const b = mol.bonds[k], A = 3 * b.a, B = 3 * b.b, q = PICK_STRIDE * j;
        for (let c = 0; c < 3; c++) {
          const a = X[A + c], d = X[B + c] - a, o = off * nn[c];
          pickSegs[q + c] = a + t0 * d + o;
          pickSegs[q + 3 + c] = a + t1 * d + o;
          pickSegs[q + 8 + c] = nn[c];
        }
        pickSegs[q + 6] = r;
        pickSegs[q + 7] = kind === RIBBON ? RIBBON : kind ? 1 : 0;
        pickSegBond[j] = k;
      });
      updateAllHalos();
      rebuildMeasures();
      rebuildElectrons();
      labelVisCache = null;
    }

    // A double, triple or aromatic bond is drawn as side-by-side pieces with gaps between them, and the middle of a
    // drawn double bond is exactly such a gap: pointing there picked nothing. So each multiple bond also gets a
    // pick-only ribbon (never drawn): the flat strip between the axes of its outermost pieces, in their common plane,
    // over the stretch where they lie side by side. Flat, so it covers the gaps and nothing beyond the drawn pieces;
    // atoms still win wherever they are in front (pick()).
    function pickEnvelope(k, parts, nn, segs) {
      let t0 = Infinity, t1 = -Infinity, lo = Infinity, hi = -Infinity;
      for (const p of parts) if (p[2] !== 0) { t0 = Math.min(t0, p[0]); t1 = Math.max(t1, p[1]); }
      if (!(t1 > t0)) return;                            // a single bond: nothing side by side
      for (const [a, b, off] of parts) if (a < t1 && b > t0) { lo = Math.min(lo, off); hi = Math.max(hi, off); }
      if (hi > lo) segs.push([k, t0, t1, (lo + hi) / 2, (hi - lo) / 2, nn, RIBBON]);
    }

    function cylinderMesh(list, geometry, ofName) {
      if (!list.length) return null;
      const mesh = new THREE.InstancedMesh(geometry, mat.solid, list.length), of = new Int32Array(list.length);
      list.forEach(([k, t0, t1, off, r, colorAtom, n], j) => {
        const b = mol.bonds[k], A = 3 * b.a, B = 3 * b.b;
        const dx = X[B] - X[A], dy = X[B + 1] - X[A + 1], dz = X[B + 2] - X[A + 2];
        _a.set(X[A] + t0 * dx + off * n[0], X[A + 1] + t0 * dy + off * n[1], X[A + 2] + t0 * dz + off * n[2]);
        _b.set(X[A] + t1 * dx + off * n[0], X[A + 1] + t1 * dy + off * n[1], X[A + 2] + t1 * dz + off * n[2]);
        setCylinder(mesh, j, _a, _b, r);
        mesh.setColorAt(j, mol.colors[colorAtom]);
        of[j] = k;
      });
      finishInstances(mesh);
      model[ofName] = of;
      inner.add(mesh);
      return mesh;
    }
    function setCylinder(mesh, j, A, B, r) {
      _d.subVectors(B, A);
      const L = _d.length();
      if (L < 1e-9) { _m.makeScale(0, 0, 0); mesh.setMatrixAt(j, _m); return; }
      _q.setFromUnitVectors(_Y, _d.divideScalar(L));
      _m.compose(_p.addVectors(A, B).multiplyScalar(0.5), _q, _s.set(r, L, r));
      mesh.setMatrixAt(j, _m);
    }
    function finishInstances(mesh) {
      mesh.instanceMatrix.needsUpdate = true;
      if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;
      mesh.computeBoundingSphere();                      // three computes it once and never refreshes it
    }

    // wireframe: 1-px lines with half-bond vertex colours; extra bonds as shortened parallel lines, aromatic as dashes
    function buildWire(n, segs) {
      const pos = [], col = [], dpos = [], dcol = [];
      // Unlit 1-px lines have no shading to set them off the background: on a light background, colours brighter
      // than linear luminance 0.5 (white H, pale elements) are dimmed to it, so they stay visible (white → #BCBCBC).
      const lineCol = mol.colors.map(c => {
        const L = 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b;
        return !darkBg && L > 0.5 ? c.clone().multiplyScalar(0.5 / L) : c;
      });
      const push = (P, C, k, t0, t1, off) => {
        const b = mol.bonds[k], A = 3 * b.a, B = 3 * b.b, cuts = t0 < 0.5 && t1 > 0.5 ? [t0, 0.5, t1] : [t0, t1];
        for (let s = 0; s + 1 < cuts.length; s++) {
          const c = lineCol[(cuts[s] + cuts[s + 1]) / 2 < 0.5 ? b.a : b.b];
          for (const t of [cuts[s], cuts[s + 1]]) {
            for (let d = 0; d < 3; d++) P.push(X[A + d] + t * (X[B + d] - X[A + d]) + off * n.getComponent(d));
            C.push(c.r, c.g, c.b);                       // linear values, as vertex colours must be
          }
        }
      };
      for (let k = 0; k < mol.M; k++) {
        if (!bondVis[k]) continue;
        const kind = bondKind(k);
        if (kind !== 's') bondNormal(k, n); else n.set(0, 0, 0);
        const nn = [n.x, n.y, n.z];
        // lines: radius 0, so picked MIN_PICK_PX wide; [t0, t1, offset, radius] as bondParts gives them
        const parts = [[0, 1, 0, 0]];
        if (kind === 'd' || kind === 'ar') parts.push([0.15, 0.85, WIRE_GAP, 0]);
        if (kind === 't') parts.push([0.15, 0.85, WIRE_GAP, 0], [0.15, 0.85, -WIRE_GAP, 0]);
        parts.forEach(([t0, t1, off], j) => {
          push(kind === 'ar' && j ? dpos : pos, kind === 'ar' && j ? dcol : col, k, t0, t1, off);
          segs.push([k, t0, t1, off, 0, nn, true]);
        });
        pickEnvelope(k, parts, nn, segs);
        bondHaloR[k] = 0.15;
      }
      const lines = (P, C, material) => {
        if (!P.length) return null;
        const g = new THREE.BufferGeometry();
        g.setAttribute('position', new THREE.Float32BufferAttribute(P, 3));
        g.setAttribute('color', new THREE.Float32BufferAttribute(C, 3));
        const o = new THREE.LineSegments(g, material);
        if (material.isLineDashedMaterial) o.computeLineDistances();
        inner.add(o);
        return o;
      };
      model.lines = lines(pos, col, mat.line);
      model.dashes = lines(dpos, dcol, mat.dash);
    }

    /* ---------- halos (hover, selection, functional group) ---------- */

    function haloFor(kind) {
      if (!halos[kind]) halos[kind] = { sel: null, sphere: null, cyl: null, capS: 0, capC: 0 };
      return halos[kind];
    }
    function updateHalo(kind) {
      const H = haloFor(kind);
      const atoms = [], bonds = [];
      if (H.sel && mol && X && atomVis) {
        for (const i of H.sel.atoms || []) if (atomVis[i]) atoms.push(i);
        if (style !== 'spacefill') for (const k of H.sel.bonds || []) if (bondVis[k]) bonds.push(k);
      }
      const ensure = (key, capKey, geometry, need) => {
        if (need > H[capKey]) {
          disposeObj(H[key]);
          H[capKey] = Math.max(need, 16, 2 * H[capKey]);
          H[key] = new THREE.InstancedMesh(geometry, mat.halo[kind], H[capKey]);
          H[key].renderOrder = 1;
          H[key].frustumCulled = false;
          inner.add(H[key]);
        }
        if (H[key]) H[key].count = need;
        return H[key];
      };
      const sm = ensure('sphere', 'capS', geo.haloSphere, atoms.length);
      atoms.forEach((i, j) => {
        const r = style === 'spacefill' ? atomR[i] + 0.2 : style === 'wireframe' && visDeg[i] ? 0.3 : atomR[i] * 1.25 + 0.1;
        sm.setMatrixAt(j, _m.compose(_p.set(X[3 * i], X[3 * i + 1], X[3 * i + 2]), _qi.identity(), _s.setScalar(r)));
      });
      if (sm) sm.instanceMatrix.needsUpdate = true;
      const cm = ensure('cyl', 'capC', geo.haloCyl, bonds.length);
      bonds.forEach((k, j) => {
        const b = mol.bonds[k];
        _a.set(X[3 * b.a], X[3 * b.a + 1], X[3 * b.a + 2]); _b.set(X[3 * b.b], X[3 * b.b + 1], X[3 * b.b + 2]);
        setCylinder(cm, j, _a, _b, bondHaloR[k] || 0.3);
      });
      if (cm) cm.instanceMatrix.needsUpdate = true;
    }
    function updateAllHalos() { for (const k of ['hover', 'select', 'group']) updateHalo(k); }

    function highlight(kind, sel) {
      if (!HALO[kind]) throw new Error('Mol3DView.highlight: kind must be hover, select or group');
      if (sel) {
        const N = mol ? mol.N : 0, M = mol ? mol.M : 0;
        for (const i of sel.atoms || []) if (!(i >= 0 && i < N)) throw new Error('Mol3DView.highlight: atom index ' + i + ' out of range');
        for (const k of sel.bonds || []) if (!(k >= 0 && k < M)) throw new Error('Mol3DView.highlight: bond index ' + k + ' out of range');
      }
      haloFor(kind).sel = sel && ((sel.atoms && sel.atoms.length) || (sel.bonds && sel.bonds.length))
        ? { atoms: (sel.atoms || []).slice(), bonds: (sel.bonds || []).slice() } : null;
      updateHalo(kind);
      requestRender();
    }

    /* ---------- measurements ---------- */

    function setMeasurements(list) {
      const N = mol ? mol.N : 0;
      const out = [];
      for (const m of list || []) {
        if (!m || !Array.isArray(m.atoms) || m.atoms.length < 2 || m.atoms.length > 4 || m.atoms.some(i => !(i >= 0 && i < N))) {
          throw new Error('Mol3DView.setMeasurements: each needs 2 to 4 atom indices in range');
        }
        out.push({ atoms: m.atoms.slice(), text: String(m.text == null ? '' : m.text) });
      }
      measures = out;
      rebuildMeasures();
      requestRender();
    }
    function rebuildMeasures() {
      disposeObj(model.measure); model.measure = null;
      if (!measures.length || !X) return;
      const P = [];
      for (const m of measures) {
        for (let s = 0; s + 1 < m.atoms.length; s++) {
          const i = m.atoms[s], j = m.atoms[s + 1];
          P.push(X[3 * i], X[3 * i + 1], X[3 * i + 2], X[3 * j], X[3 * j + 1], X[3 * j + 2]);
        }
      }
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.Float32BufferAttribute(P, 3));
      const o = new THREE.LineSegments(g, mat.measure);
      o.computeLineDistances();
      o.renderOrder = 4;
      o.frustumCulled = false;
      inner.add(o);
      model.measure = o;
    }
    // where a measurement's label goes, in molecule coordinates
    function measureAnchor(m) {
      const p = i => [X[3 * i], X[3 * i + 1], X[3 * i + 2]], a = m.atoms;
      if (a.length === 2 || a.length === 4) {
        const u = p(a[a.length === 2 ? 0 : 1]), v = p(a[a.length === 2 ? 1 : 2]);
        return [(u[0] + v[0]) / 2, (u[1] + v[1]) / 2, (u[2] + v[2]) / 2];
      }
      const c = p(a[1]), u = p(a[0]), v = p(a[2]);       // an angle: inside it, near the vertex
      const du = Math.hypot(u[0] - c[0], u[1] - c[1], u[2] - c[2]) || 1, dv = Math.hypot(v[0] - c[0], v[1] - c[1], v[2] - c[2]) || 1;
      const w = [0, 1, 2].map(d => (u[d] - c[d]) / du + (v[d] - c[d]) / dv), lw = Math.hypot(w[0], w[1], w[2]) || 1;
      return [0, 1, 2].map(d => c[d] + 0.45 * w[d] / lw);
    }

    /* ---------- surface ---------- */

    function removeSurface() {
      if (!surf) return;
      if (surf.mesh) inner.remove(surf.mesh);
      if (surf.pre) inner.remove(surf.pre);
      if (surf.geom) surf.geom.dispose();
      if (surf.material) surf.material.dispose();
      surf = null;
    }

    function setSurface(cfg) {
      const gen = ++surfGen;
      if (!cfg) { surfCfg = null; removeSurface(); requestRender(); return Promise.resolve(null); }
      const mode = cfg.mode === 'potential' ? 'potential' : 'plain';
      const opacity = cfg.opacity == null ? 0.75 : Math.max(0.05, Math.min(1, +cfg.opacity));
      const next = { mode, opacity, charges: cfg.charges || null, clamp: Number.isFinite(cfg.clamp) && cfg.clamp > 0 ? +cfg.clamp : null };
      // no coordinates yet: kept, built when setCoords brings them, and reported through onSurface
      if (!mol || !X) { surfCfg = mol ? next : null; removeSurface(); return Promise.resolve(null); }
      // only the opacity changed: keep the geometry
      if (surf && surf.X === X && surfCfg && surfCfg.mode === next.mode && surfCfg.clamp === next.clamp && sameCharges(surfCfg.charges, next.charges)) {
        surfCfg = next;
        surf.opacity = opacity;
        styleSurface(surf, opacity);
        requestRender();
        return Promise.resolve(surf.result);
      }
      surfCfg = next;
      removeSurface();
      requestRender();
      return makeSurface(next, gen);
    }
    function sameCharges(a, b) {
      if (a === b) return true;
      if (!a || !b || a.length !== b.length) return false;
      for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) return false;
      return true;
    }

    async function makeSurface(cfg, gen) {
      if (!mol || !X) throw new Error('Mol3DView.setSurface: no molecule coordinates yet');
      const N = mol.N;
      for (let i = 0; i < N; i++) {
        if (!mol.hasVdw[i]) throw new Error('no van der Waals radius for ' + mol.el[i] + ' (atom ' + i + '), so no surface is drawn');
      }
      let Q = null;
      if (cfg.mode === 'potential') {
        if (!cfg.charges || cfg.charges.length !== N) throw new Error('the potential surface needs one partial charge per atom (' + N + ')');
        Q = Float64Array.from(cfg.charges);
        for (let i = 0; i < N; i++) if (!Number.isFinite(Q[i])) throw new Error('no partial charge for atom ' + i + ' (' + mol.el[i] + '), so the potential is not drawn');
      }
      const t0 = performance.now(), coords = X;
      const cancelled = () => gen !== surfGen || disposed || X !== coords;
      const S = await buildSurface(coords, mol.rvdw, { cancelled });
      if (!S || cancelled()) return null;
      let result = null, V = null;
      if (Q) {
        V = await coulomb(S.positions, coords, Q, cancelled);
        if (!V || cancelled()) return null;
        let vmin = Infinity, vmax = -Infinity;
        const absV = new Float64Array(V.length);
        for (let v = 0; v < V.length; v++) { vmin = Math.min(vmin, V[v]); vmax = Math.max(vmax, V[v]); absV[v] = Math.abs(V[v]); }
        const auto = percentile(absV, 0.98);
        result = { vmin, vmax, clamp: cfg.clamp || auto, autoClamp: auto, fixedClamp: !!cfg.clamp,
          p2: percentile(V, 0.02), p98: percentile(V, 0.98), units: 'kcal/mol per e' };
      }
      let colors = null;
      if (V) {
        colors = new Float32Array(3 * V.length);
        const cl = result.clamp;
        for (let v = 0; v < V.length; v++) {
          const rgb = espRGB(cl > 0 ? V[v] / cl : 0);
          _c.setRGB(rgb[0], rgb[1], rgb[2], THREE.SRGBColorSpace);   // stored linear, as vertex colours must be
          colors[3 * v] = _c.r; colors[3 * v + 1] = _c.g; colors[3 * v + 2] = _c.b;
        }
      }
      const ms = performance.now() - t0;
      stats.surfaceMs = ms;
      const info = { vertices: S.positions.length / 3, triangles: S.index.length / 3, gridStep: S.h, ms: Math.round(ms) };
      result = result ? Object.assign(result, info) : null;
      removeSurface();
      surf = { X: coords, result, V, colors, data: S, info, opacity: cfg.opacity, geom: null, material: null, mesh: null, pre: null };
      surfaceMeshes(surf);
      if (!userMoved) resetView();                       // untouched view: frame the surface too
      requestRender();
      return result;
    }
    // the GPU side of a surface, from its arrays (again after a context loss)
    function surfaceMeshes(s) {
      const g = new THREE.BufferGeometry(), S = s.data;
      g.setAttribute('position', new THREE.BufferAttribute(S.positions, 3));
      g.setAttribute('normal', new THREE.BufferAttribute(S.normals, 3));
      g.setIndex(new THREE.BufferAttribute(S.index, 1));
      if (s.colors) g.setAttribute('color', new THREE.BufferAttribute(s.colors, 3));
      g.computeBoundingSphere();
      s.geom = g;
      s.material = new THREE.MeshStandardMaterial({ roughness: 0.6, metalness: 0, vertexColors: !!s.colors, color: s.colors ? 0xffffff : PLAIN_SURFACE });
      s.mesh = new THREE.Mesh(g, s.material);
      s.pre = new THREE.Mesh(g, mat.surfPre);             // depth-only pass, so only the front layer is blended
      s.pre.renderOrder = 2; s.mesh.renderOrder = 3;
      styleSurface(s, s.opacity);
      inner.add(s.pre, s.mesh);
    }
    function styleSurface(s, opacity) {
      const m = s.material, see = opacity < 0.999;
      m.transparent = see; m.opacity = see ? opacity : 1; m.depthWrite = !see; m.depthFunc = THREE.LessEqualDepth;
      m.needsUpdate = true;
      s.pre.visible = see;
    }

    /* ---------- camera and view ---------- */

    function tanHalf() { return Math.tan(FOV * Math.PI / 360); }
    function updateCamera() {
      const w = cssW || 1, h = cssH || 1, aspect = w / h, R = view.R;
      if (orthoOn) {
        const hh = view.dist * tanHalf(), z = 3 * R + 10;
        Object.assign(ortho, { left: -hh * aspect, right: hh * aspect, top: hh, bottom: -hh, near: 0.01, far: z + 2 * R });
        ortho.position.set(0, 0, z);
        camera = ortho;
      } else {
        persp.aspect = aspect;
        persp.near = Math.max(0.05, view.dist - 1.6 * R);
        persp.far = view.dist + 1.6 * R;
        persp.position.set(0, 0, view.dist);
        camera = persp;
      }
      camera.updateProjectionMatrix();
      camera.updateMatrixWorld();
    }
    // Å per CSS pixel at distance t from the camera along a ray
    function aPerPx(t) {
      const h = cssH || 1;
      return orthoOn ? 2 * view.dist * tanHalf() / h : 2 * Math.max(t, 0) * tanHalf() / h;
    }

    function resetView() {
      inertia = null;
      molGroup.position.set(0, 0, 0);
      if (!mol || !X || !pca) { view.dist = 10; updateCamera(); requestRender(); return; }
      // principal axes onto the screen (largest spread across, flattest toward the viewer), tilted a little so
      // flat molecules still read as 3D and linear ones show their multiple bonds
      const [e1, e2, e3] = pca.axes;
      _m.set(e1[0], e1[1], e1[2], 0, e2[0], e2[1], e2[2], 0, e3[0], e3[1], e3[2], 0, 0, 0, 0, 1);
      molGroup.quaternion.setFromRotationMatrix(_m);
      molGroup.quaternion.premultiply(_q.setFromAxisAngle(_d.set(1, 0, 0), -0.35));
      molGroup.quaternion.premultiply(_q.setFromAxisAngle(_d.set(0, 1, 0), 0.25));
      molGroup.updateMatrixWorld(true);
      // distance at which every visible atom fits, with a margin
      const w = cssW || 1, h = cssH || 1, tv = tanHalf(), th = tv * w / h;
      let dist = 0;
      const vis = atomVis, c = pca.center;
      for (let i = 0; i < mol.N; i++) {
        if (vis && !vis[i]) continue;
        _p.set(X[3 * i] - c[0], X[3 * i + 1] - c[1], X[3 * i + 2] - c[2]).applyQuaternion(molGroup.quaternion);
        const r = Math.max((atomR ? atomR[i] : BALL_SCALE * mol.rvdw[i]) + 0.3, surf ? mol.rvdw[i] : 0);   // a surface must fit too
        dist = Math.max(dist, (Math.abs(_p.x) + r) * 1.06 / th + _p.z, (Math.abs(_p.y) + r) * 1.06 / tv + _p.z);
      }
      view.dist = Math.max(dist, 0.5 * view.R, 3);
      framedWithoutSize = !cssW || !cssH;
      userMoved = false;
      labelVisCache = null;
      updateCamera();
      requestRender();
    }

    function zoomBy(f) {
      view.dist = Math.max(0.5 * view.R, Math.min(20 * view.R, view.dist * f));
      userMoved = true;
      updateCamera(); labelVisCache = null; requestRender();
    }
    function rotateWorld(ax, ay, az, angle) {
      if (!angle) return;
      _d.set(ax, ay, az).normalize();
      molGroup.quaternion.premultiply(_q.setFromAxisAngle(_d, angle)).normalize();
      userMoved = true; labelVisCache = null;
    }
    function panBy(dxPx, dyPx) {
      const s = aPerPx(view.dist);
      molGroup.position.x += dxPx * s; molGroup.position.y -= dyPx * s;
      userMoved = true; labelVisCache = null;
    }
    function setOrthographic(on) { orthoOn = !!on; updateCamera(); labelVisCache = null; requestRender(); }
    function setSpin(on) { spin = !!on; if (spin) inertia = null; requestRender(); }

    /* ---------- picking ---------- */

    // the pointer's ray in molecule coordinates (the group is never scaled, so distances along it are in Å)
    function rayAt(clientX, clientY) {
      const rect = canvas.getBoundingClientRect();
      if (!rect.width || !rect.height) return null;
      _ndc.set(((clientX - rect.left) / rect.width) * 2 - 1, -((clientY - rect.top) / rect.height) * 2 + 1);
      raycaster.setFromCamera(_ndc, camera);
      molGroup.updateMatrixWorld(true);
      _inv.copy(inner.matrixWorld).invert();
      return _ray.copy(raycaster.ray).applyMatrix4(_inv);
    }
    // nearest atom along a ray: [atom, t] with radius per style (at least MIN_PICK_PX wide); -1 if none
    function hitAtoms(ray, radii, minPx, skip = -1) {
      const o = ray.origin, d = ray.direction;
      let best = Infinity, who = -1;
      for (let i = 0; i < mol.N; i++) {
        if (!atomVis[i] || i === skip) continue;
        const cx = X[3 * i] - o.x, cy = X[3 * i + 1] - o.y, cz = X[3 * i + 2] - o.z;
        const tca = cx * d.x + cy * d.y + cz * d.z;
        let r = radii[i];
        if (minPx) r = Math.max(r, MIN_PICK_PX * aPerPx(tca));
        if (r <= 0) continue;
        const d2 = cx * cx + cy * cy + cz * cz - tca * tca;
        if (d2 > r * r) continue;
        const thc = Math.sqrt(r * r - d2);
        let t = tca - thc;
        if (t < 0) t = tca + thc;
        if (t >= 0 && t < best) { best = t; who = i; }
      }
      return [who, best];
    }
    // where a ray first meets a bond cylinder (side; and the end discs if it is closed): t, or Infinity
    function hitCylinder(o, d, q) {
      const ax = pickSegs[q], ay = pickSegs[q + 1], az = pickSegs[q + 2];
      let ux = pickSegs[q + 3] - ax, uy = pickSegs[q + 4] - ay, uz = pickSegs[q + 5] - az;
      const L = Math.sqrt(ux * ux + uy * uy + uz * uz);
      if (L < 1e-9) return Infinity;
      ux /= L; uy /= L; uz /= L;
      const wx = o.x - ax, wy = o.y - ay, wz = o.z - az;
      const du = d.x * ux + d.y * uy + d.z * uz, wu = wx * ux + wy * uy + wz * uz;
      const px = d.x - du * ux, py = d.y - du * uy, pz = d.z - du * uz;      // parts across the axis
      const qx = wx - wu * ux, qy = wy - wu * uy, qz = wz - wu * uz;
      const a = px * px + py * py + pz * pz, bh = qx * px + qy * py + qz * pz, qq = qx * qx + qy * qy + qz * qz;
      const tc = a > 1e-12 ? -bh / a : Math.abs(wu);                          // closest approach, for the pixel size
      const r = Math.max(pickSegs[q + 6], MIN_PICK_PX * aPerPx(tc));
      let best = Infinity;
      if (a > 1e-12) {
        const disc = bh * bh - a * (qq - r * r);
        if (disc < 0) return Infinity;
        const t = (-bh - Math.sqrt(disc)) / a, s = wu + t * du;             // entry through the side
        if (t >= 0 && s >= 0 && s <= L) best = t;
      }
      if (pickSegs[q + 7] && Math.abs(du) > 1e-12) {                        // closed: the end discs
        for (const s of [0, L]) {
          const t = (s - wu) / du;
          if (t < 0 || t >= best) continue;
          const ex = qx + t * px, ey = qy + t * py, ez = qz + t * pz;
          if (ex * ex + ey * ey + ez * ez <= r * r) best = t;
        }
      }
      return best;
    }
    // where a ray meets a pick-only ribbon (pickEnvelope): t, or Infinity. Edge-on it is missed, and nothing is lost:
    // the pieces it joins then lie one behind the other on screen.
    function hitRibbon(o, d, q) {
      const ax = pickSegs[q], ay = pickSegs[q + 1], az = pickSegs[q + 2];
      let ux = pickSegs[q + 3] - ax, uy = pickSegs[q + 4] - ay, uz = pickSegs[q + 5] - az;
      const L = Math.sqrt(ux * ux + uy * uy + uz * uz);
      if (L < 1e-9) return Infinity;
      ux /= L; uy /= L; uz /= L;
      const nx = pickSegs[q + 8], ny = pickSegs[q + 9], nz = pickSegs[q + 10];
      const mx = uy * nz - uz * ny, my = uz * nx - ux * nz, mz = ux * ny - uy * nx;   // the ribbon's plane normal
      const den = d.x * mx + d.y * my + d.z * mz;
      if (Math.abs(den) < 1e-9) return Infinity;
      const t = ((ax - o.x) * mx + (ay - o.y) * my + (az - o.z) * mz) / den;
      if (t < 0) return Infinity;
      const wx = o.x + t * d.x - ax, wy = o.y + t * d.y - ay, wz = o.z + t * d.z - az;
      const s = wx * ux + wy * uy + wz * uz, v = wx * nx + wy * ny + wz * nz;
      return s >= 0 && s <= L && Math.abs(v) <= pickSegs[q + 6] ? t : Infinity;
    }
    // atoms-only picking, when the ray meets no atom: the atom whose edge (its pick radius) is nearest the ray in
    // screen pixels, within NEAR_ATOM_PX; on a tie the one nearer the viewer. → [atom, t] or [-1, Infinity]
    function nearestAtom(ray) {
      const o = ray.origin, d = ray.direction;
      let best = Infinity, who = -1, bt = Infinity;
      for (let i = 0; i < mol.N; i++) {
        if (!atomVis[i]) continue;
        const cx = X[3 * i] - o.x, cy = X[3 * i + 1] - o.y, cz = X[3 * i + 2] - o.z;
        const tca = cx * d.x + cy * d.y + cz * d.z;
        if (tca <= 0) continue;
        const px = aPerPx(tca), gap = (Math.sqrt(Math.max(0, cx * cx + cy * cy + cz * cz - tca * tca)) - atomPickR[i]) / px;
        if (gap > NEAR_ATOM_PX) continue;
        if (gap < best - 0.5 || (Math.abs(gap - best) <= 0.5 && tca < bt)) { best = gap; who = i; bt = tca; }
      }
      return [who, bt];
    }
    function setPickMode(mode) {
      if (mode !== 'all' && mode !== 'atoms') throw new Error("Mol3DView.setPickMode: 'all' or 'atoms'");
      atomsOnly = mode === 'atoms';
      hover.key = ''; hover.pending = hover.inside;       // the thing under the pointer may change
      requestRender();
    }
    let lastPick = null;
    function pick(clientX, clientY) {
      if (!mol || !X || !atomVis) return null;
      const t0 = performance.now();
      const ray = rayAt(clientX, clientY);
      if (!ray) return null;
      let [atom, ta] = hitAtoms(ray, atomPickR, true);
      let bond = -1, tb = Infinity;
      if (style !== 'spacefill') {                       // space-filling bonds are buried: atoms only
        for (let j = 0; j < nPickSegs; j++) {
          const q = PICK_STRIDE * j;
          const t = pickSegs[q + 7] === RIBBON ? hitRibbon(ray.origin, ray.direction, q) : hitCylinder(ray.origin, ray.direction, q);
          if (t < tb) { tb = t; bond = pickSegBond[j]; }
        }
      }
      let hit = null;
      if (atomsOnly) {
        // atoms only: the atom under the pointer, even behind a bond drawn in front of it (a C–H bond over a methyl
        // carbon); else, on a bond, the end nearer the pointer on screen; else the nearest atom within NEAR_ATOM_PX
        if (atom < 0 && bond >= 0) {
          const b = mol.bonds[bond], o = ray.origin, d = ray.direction;
          const off = i => { const cx = X[3 * i] - o.x, cy = X[3 * i + 1] - o.y, cz = X[3 * i + 2] - o.z, t = cx * d.x + cy * d.y + cz * d.z;
            return Math.sqrt(Math.max(0, cx * cx + cy * cy + cz * cz - t * t)) / aPerPx(Math.max(t, 1e-3)); };
          atom = off(b.a) <= off(b.b) ? b.a : b.b; ta = tb;
        }
        if (atom < 0) [atom, ta] = nearestAtom(ray);
        if (atom >= 0) hit = { type: 'atom', index: atom };
        bond = -1; tb = Infinity;
      }
      const ms = performance.now() - t0;
      stats.pickMs = ms; stats.picks++; stats.pickAvgMs += (ms - stats.pickAvgMs) / Math.min(stats.picks, 60);
      // a bond wins only if it is in front of the atom; where a bond's surface meets its atom's (a stick joint) the
      // two are level, and the atom wins
      if (atomsOnly) { /* decided above */ }
      else if (bond >= 0 && tb < ta - PICK_TIE) hit = { type: 'bond', index: bond };
      else if (atom >= 0) hit = { type: 'atom', index: atom };
      else if (bond >= 0) hit = { type: 'bond', index: bond };
      lastPick = hit ? { hit, t: hit.type === 'atom' ? ta : tb, tAtom: ta, tBond: tb } : null;
      return hit;
    }

    function doHover() {
      hover.pending = false;
      if (!hover.inside || drag || pointers.size) return;
      const hit = pick(hover.x, hover.y);
      const key = hit ? hit.type + hit.index : '';
      if (key !== hover.key) {
        hover.key = key;
        haloFor('hover').sel = hit ? (hit.type === 'atom' ? { atoms: [hit.index], bonds: [] } : { atoms: [], bonds: [hit.index] }) : null;
        updateHalo('hover');
        canvas.style.cursor = hit ? 'pointer' : '';
      }
      call('onHover', hit, hover.x, hover.y);
    }
    function clearHover() {
      if (hover.key) {
        hover.key = '';
        haloFor('hover').sel = null; updateHalo('hover');
        canvas.style.cursor = '';
        requestRender();
      }
      call('onHover', null, hover.x, hover.y);
    }

    /* ---------- labels (2D overlay) ---------- */

    function setLabels(arr) {
      if (arr != null) {
        if (!mol || !Array.isArray(arr) || arr.length !== mol.N) throw new Error('Mol3DView.setLabels: expected ' + (mol ? mol.N : 0) + ' labels');
        labels = arr.map(s => (s == null ? '' : String(s)));
      } else labels = null;
      labelVisCache = null;
      requestRender();
    }
    function labelCount() { let n = 0; if (labels) for (let i = 0; i < labels.length; i++) if (labels[i] && atomVis && atomVis[i]) n++; return n; }
    // which labels are not hidden behind another atom: cast from the camera to a point just inside the atom's
    // front surface; the label shows if no other atom is hit first
    function labelVisibility() {
      const t0 = performance.now(), N = mol.N, vis = new Uint8Array(N);
      molGroup.updateMatrixWorld(true);
      _inv.copy(inner.matrixWorld).invert();
      const camLocal = new THREE.Vector3().copy(camera.position).applyMatrix4(_inv);
      const back = new THREE.Vector3(0, 0, 1).transformDirection(_inv);   // toward the viewer, molecule axes
      const ray = new THREE.Ray(), toCam = new THREE.Vector3(), P = new THREE.Vector3();
      for (let i = 0; i < N; i++) {
        if (!labels[i] || !atomVis[i]) continue;
        if (style === 'wireframe' || atomR[i] <= 0) { vis[i] = 1; continue; }
        const c = _p.set(X[3 * i], X[3 * i + 1], X[3 * i + 2]);
        if (orthoOn) toCam.copy(back); else toCam.subVectors(camLocal, c).normalize();
        P.copy(c).addScaledVector(toCam, 0.9 * atomR[i]);
        if (orthoOn) ray.origin.copy(P).addScaledVector(back, 1000); else ray.origin.copy(camLocal);
        ray.direction.subVectors(P, ray.origin);
        const tP = ray.direction.length();
        ray.direction.divideScalar(tP);
        const [who, t] = hitAtoms(ray, atomR, false, i);
        vis[i] = who < 0 || t > tP ? 1 : 0;
      }
      stats.labelMs = performance.now() - t0;
      return vis;
    }
    function project(x, y, z, out) {
      _p.set(x, y, z).applyMatrix4(inner.matrixWorld).project(camera);
      out.x = (_p.x + 1) / 2 * cssW; out.y = (1 - _p.y) / 2 * cssH; out.z = _p.z;
      return out;
    }
    // labels onto a 2D context: the overlay canvas (cleared first), or the export image (final: drawn over the
    // 3D picture, and with every label's occlusion worked out now)
    function drawOverlay(ctx, scale, final) {
      ctx.setTransform(scale, 0, 0, scale, 0, 0);
      if (!final) ctx.clearRect(0, 0, cssW, cssH);
      labelStats = { shown: 0, overlapping: 0, hOnHover: 0, hidden: 0 };
      if (!final) { lastLabels = []; lastMeasures = []; }
      if (!mol || !X || !atomVis) return;
      molGroup.updateMatrixWorld(true);
      const pt = { x: 0, y: 0, z: 0 };
      ctx.textAlign = 'center'; ctx.textBaseline = 'middle'; ctx.lineJoin = 'round';
      // Boxes already taken on the overlay, in a grid of 48-px cells so the overlap tests stay linear. Nothing is
      // ever drawn over a label: a number half hidden under another could read as a different number.
      const taken = boxGrid(48);
      // measurement labels first: always shown, on top, with the caller's units
      const mBoxes = [];
      if (measures.length) {
        ctx.font = LABEL_FONT;
        for (const m of measures) {
          if (!m.text) continue;
          const a = measureAnchor(m);
          project(a[0], a[1], a[2], pt);
          if (Math.abs(pt.z) > 1) continue;
          const w = ctx.measureText(m.text).width + 10, h = 18;
          let x = pt.x - w / 2, y = pt.y - h / 2;
          // step down past measurement labels already placed, so nearby measurements stay readable
          for (let n = 0; n < 8 && mBoxes.some(r => x < r.x + r.w && r.x < x + w && y < r.y + r.h && r.y < y + h); n++) y += h + 2;
          mBoxes.push({ x, y, w, h, text: m.text });
          taken.add([x, y, x + w, y + h]);
        }
      }
      if (labels) {
        const heavy = labelCount() * mol.N > LABEL_EAGER;
        if (final || !heavy || !moving) {
          if (!labelVisCache || !heavy) labelVisCache = labelVisibility();
          const vis = labelVisCache;
          _inv.copy(inner.matrixWorld).invert();
          const camLocal = new THREE.Vector3().copy(camera.position).applyMatrix4(_inv);
          const back = new THREE.Vector3(0, 0, 1).transformDirection(_inv), toCam = new THREE.Vector3();
          ctx.font = LABEL_FONT;
          // with many heavy atoms the H labels crowd them out: then an H is labelled only while it is hovered
          let nHeavy = 0;
          for (let i = 0; i < mol.N; i++) if (!mol.hideable[i]) nHeavy++;
          const hOnlyOnHover = nHeavy > H_LABELS_MAX_HEAVY;
          const hovered = !final && hover.key.startsWith('atom') ? +hover.key.slice(4) : -1;
          const cand = [];
          for (let i = 0; i < mol.N; i++) {
            if (!labels[i]) continue;
            // behind another atom: left out, and counted (an H left for hover counts as that, not here)
            if (!vis[i]) { if (atomVis[i] && !(hOnlyOnHover && mol.hideable[i] && i !== hovered)) labelStats.hidden++; continue; }
            if (hOnlyOnHover && mol.hideable[i] && i !== hovered) { labelStats.hOnHover++; continue; }
            const x = X[3 * i], y = X[3 * i + 1], z = X[3 * i + 2];
            if (orthoOn) toCam.copy(back); else toCam.set(camLocal.x - x, camLocal.y - y, camLocal.z - z).normalize();
            const r = atomR[i];
            project(x + r * toCam.x, y + r * toCam.y, z + r * toCam.z, pt);   // the atom's point nearest the viewer
            if (Math.abs(pt.z) > 1) continue;
            cand.push({ i, x: pt.x, y: pt.y, z: pt.z, w: ctx.measureText(labels[i]).width, rank: i === hovered ? 0 : mol.hideable[i] ? 2 : 1 });
          }
          // on screen, the corner where the note below goes is kept clear
          const note = !final && cand.length ? [0, cssH - NOTE_H, Math.min(cssW, NOTE_W), cssH] : null;
          if (note) taken.add(note);
          // nearest first (the hovered atom, then heavy atoms, then H); a label whose box would overlap one
          // already placed is left out rather than printed over it
          cand.sort((p, q) => (p.rank - q.rank) || (p.z - q.z));
          const placed = [];
          for (const c of cand) {
            const b = [c.x - c.w / 2 - LABEL_PAD, c.y - LABEL_H / 2, c.x + c.w / 2 + LABEL_PAD, c.y + LABEL_H / 2];
            if (taken.hits(b)) { labelStats.overlapping++; continue; }
            taken.add(b);
            c.box = b;
            placed.push(c);
          }
          if (!final) lastLabels = placed.map(c => ({ index: c.i, text: labels[c.i], box: c.box }));
          ctx.lineWidth = 3;
          ctx.strokeStyle = darkBg ? 'rgba(0,0,0,0.85)' : 'rgba(255,255,255,0.9)';
          ctx.fillStyle = darkBg ? '#f2f2f2' : '#111';
          for (const c of placed) { ctx.strokeText(labels[c.i], c.x, c.y); ctx.fillText(labels[c.i], c.x, c.y); }
          labelStats.shown = placed.length;
          // on screen only (not in an export): say that some labels are left out, so an unlabelled atom is not
          // taken for one without a value
          if (note && (labelStats.overlapping || labelStats.hOnHover || labelStats.hidden)) {
            const s = labelStats.overlapping, hb = labelStats.hidden, left = [], parts = [];
            if (s) left.push(`${s} overlapping (zoom in or hover)`);
            if (hb) left.push(`${hb} on atoms hidden behind others (rotate)`);
            if (left.length) parts.push((s + hb === 1 ? 'Label left out: ' : 'Labels left out: ') + left.join(', '));
            if (labelStats.hOnHover) parts.push('H labels show on hover');
            ctx.save();
            ctx.font = '11px system-ui, "Segoe UI", sans-serif';
            ctx.textAlign = 'left'; ctx.textBaseline = 'bottom';
            ctx.fillStyle = darkBg ? '#c8c8c8' : '#555';
            const t = fitText(ctx, parts.join(' · '), note[2] - 12);
            ctx.strokeText(t, 8, cssH - 5);
            ctx.fillText(t, 8, cssH - 5);
            ctx.restore();
          }
        }
      }
      ctx.font = LABEL_FONT;
      if (!final) lastMeasures = mBoxes;
      for (const m of mBoxes) {
        ctx.fillStyle = darkBg ? 'rgba(20,20,20,0.85)' : 'rgba(255,255,255,0.9)';
        ctx.strokeStyle = MEASURE_COLOR; ctx.lineWidth = 1;
        roundRect(ctx, m.x, m.y, m.w, m.h, 4);
        ctx.fill(); ctx.stroke();
        ctx.fillStyle = darkBg ? '#ffb3d0' : MEASURE_COLOR;
        ctx.fillText(m.text, m.x + m.w / 2, m.y + m.h / 2 + 0.5);
      }
    }
    // screen boxes [x0, y0, x1, y1] bucketed by cell: add(b); hits(b) → it overlaps one added
    function boxGrid(cell) {
      const g = new Map();
      const each = (b, f) => {
        for (let gx = Math.floor(b[0] / cell); gx <= Math.floor(b[2] / cell); gx++) {
          for (let gy = Math.floor(b[1] / cell); gy <= Math.floor(b[3] / cell); gy++) if (f(gx + ',' + gy)) return true;
        }
        return false;
      };
      return {
        add: b => { each(b, k => { const l = g.get(k); if (l) l.push(b); else g.set(k, [b]); return false; }); },
        hits: b => each(b, k => { const l = g.get(k); return !!l && l.some(o => b[0] < o[2] && o[0] < b[2] && b[1] < o[3] && o[1] < b[3]); }),
      };
    }
    // the text, cut with an ellipsis to fit maxW (CSS px) in the context's font
    function fitText(ctx, text, maxW) {
      if (ctx.measureText(text).width <= maxW) return text;
      let s = text;
      while (s.length > 1 && ctx.measureText(s + '…').width > maxW) s = s.slice(0, -1);
      return s + '…';
    }
    function roundRect(ctx, x, y, w, h, r) {
      ctx.beginPath();
      ctx.moveTo(x + r, y); ctx.arcTo(x + w, y, x + w, y + h, r); ctx.arcTo(x + w, y + h, x, y + h, r);
      ctx.arcTo(x, y + h, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath();
    }

    /* ---------- render loop (on demand) ---------- */

    function requestRender() {
      if (inFrame || raf || lost || disposed) return;
      raf = requestAnimationFrame(frame);
    }
    // can the view be seen? Not while its box is not displayed (0×0: another tab of the app) or the page is hidden
    function visible() { return !document.hidden && wrap.clientWidth > 0 && wrap.clientHeight > 0; }
    function frame(now) {
      raf = 0;
      if (lost || disposed) return;
      // Not visible: draw nothing and ask for no next frame, so spin and glide stop costing CPU/GPU; resize() (the
      // box shown again) or visibilitychange (the page shown again) asks for a frame and they carry on from there.
      if (!visible()) { lastT = 0; unseen = true; return; }
      let animating = false, wasMoving = moving;
      inFrame = true;
      try {
        const dt = lastT ? Math.min(0.1, (now - lastT) / 1000) : 0;
        if (spin && !drag) { rotateWorld(0, 1, 0, SPIN_RATE * dt); animating = true; }
        else if (inertia && !drag) {
          const w = inertia.length();
          if (w < INERTIA_STOP) inertia = null;
          else { rotateWorld(inertia.x, inertia.y, inertia.z, w * dt); inertia.multiplyScalar(Math.exp(-dt / INERTIA_TAU)); animating = true; }
        }
        lastT = animating ? now : 0;
        moving = animating || !!drag;
        molGroup.updateMatrixWorld(true);
        if (hover.pending || (animating && hover.inside)) doHover();
        renderNow();
      } finally {
        inFrame = false;
      }
      if (animating) requestRender();
      if (wasMoving && !moving) settleLabels();
    }
    function renderNow() {
      if (!cssW || !cssH) return;
      const t0 = performance.now();
      renderer.render(scene, camera);
      drawOverlay(octx, overlay.width / cssW, false);
      const ms = performance.now() - t0;
      stats.renderMs = ms; stats.frames++;
      stats.renderAvgMs += (ms - stats.renderAvgMs) / Math.min(stats.frames, 60);
    }
    // heavy label sets are hidden while the molecule moves; redo their occlusion ~100 ms after it stops
    function settleLabels() {
      clearTimeout(settleTimer);
      settleTimer = setTimeout(() => { labelVisCache = null; requestRender(); }, 100);
    }

    /* ---------- input ---------- */

    const listeners = [];
    const on = (target, type, fn, o) => { target.addEventListener(type, fn, o); listeners.push([target, type, fn, o]); };

    on(canvas, 'pointerdown', e => {
      if (lost) return;
      canvas.focus({ preventScroll: true });
      try { canvas.setPointerCapture(e.pointerId); } catch (err) { /* pointer already gone */ }
      pointers.set(e.pointerId, { x: e.clientX, y: e.clientY, x0: e.clientX, y0: e.clientY, t0: performance.now() });
      inertia = null;
      if (pointers.size === 1) {
        const mode = e.button === 1 || e.button === 2 || e.ctrlKey || e.metaKey ? 'pan' : e.shiftKey ? 'roll' : 'rotate';
        drag = { mode, moved: 0, samples: [], button: e.button };
      } else if (pointers.size === 2) {
        drag = Object.assign(drag || { moved: 0, samples: [] }, { mode: 'pinch', moved: 99 });
        drag.pinch = pinchState();
      }
      hover.inside = false;
      clearHover();
      e.preventDefault();
    });
    function pinchState() {
      const [p, q] = [...pointers.values()];
      return { d: Math.hypot(q.x - p.x, q.y - p.y) || 1, a: Math.atan2(q.y - p.y, q.x - p.x), mx: (p.x + q.x) / 2, my: (p.y + q.y) / 2 };
    }
    on(canvas, 'pointermove', e => {
      const p = pointers.get(e.pointerId);
      if (!p) {                                          // no button down: hover
        if (e.pointerType === 'touch') return;
        hover.x = e.clientX; hover.y = e.clientY; hover.inside = true; hover.pending = true;
        requestRender();
        return;
      }
      const dx = e.clientX - p.x, dy = e.clientY - p.y;
      p.x = e.clientX; p.y = e.clientY;
      if (!drag) return;
      drag.moved += Math.abs(dx) + Math.abs(dy);
      const h = cssH || 1;
      if (drag.mode === 'pinch' && pointers.size >= 2) {
        const s = pinchState(), s0 = drag.pinch;
        zoomBy(s0.d / s.d);
        panBy(s.mx - s0.mx, s.my - s0.my);
        rotateWorld(0, 0, 1, -(s.a - s0.a));              // screen y points down, world y up
        drag.pinch = s;
      } else if (drag.mode === 'rotate') {
        if (!dx && !dy) return;
        const ang = Math.hypot(dx, dy) * Math.PI / h;
        rotateWorld(dy, dx, 0, ang);                    // the front surface follows the pointer
        const now = performance.now();
        drag.samples.push([now, dy * Math.PI / h, dx * Math.PI / h]);
        while (drag.samples.length && now - drag.samples[0][0] > 60) drag.samples.shift();
      } else if (drag.mode === 'roll') {
        const r = canvas.getBoundingClientRect(), cx = r.left + r.width / 2, cy = r.top + r.height / 2;
        const a1 = Math.atan2(-(e.clientY - cy), e.clientX - cx), a0 = Math.atan2(-(e.clientY - dy - cy), e.clientX - dx - cx);
        rotateWorld(0, 0, 1, a1 - a0);
      } else if (drag.mode === 'pan') {
        panBy(dx, dy);
      }
      requestRender();
    });
    function endPointer(e, cancelled) {
      const p = pointers.get(e.pointerId);
      if (!p) return;
      pointers.delete(e.pointerId);
      try { canvas.releasePointerCapture(e.pointerId); } catch (err) { /* already released */ }
      if (!drag) return;
      const wasClick = !cancelled && drag.mode !== 'pinch' && drag.moved < 5 && pointers.size === 0 && (drag.button === 0 || e.pointerType !== 'mouse');
      if (pointers.size === 1 && drag.mode === 'pinch') {   // one finger left: carry on rotating with it
        drag = { mode: 'rotate', moved: 99, samples: [] };
        return;
      }
      if (pointers.size) return;
      if (drag.mode === 'rotate' && drag.samples.length > 1 && !spin) {
        // angular velocity over the last ~50 ms of the drag (sample i covers the time since sample i − 1)
        const now = performance.now(), s = drag.samples, dt = (s[s.length - 1][0] - s[0][0]) / 1000;
        if (now - s[s.length - 1][0] < 50 && dt > 0.005) {
          let ax = 0, ay = 0;
          for (let i = 1; i < s.length; i++) { ax += s[i][1]; ay += s[i][2]; }
          inertia = new THREE.Vector3(ax / dt, ay / dt, 0);
          if (inertia.length() < INERTIA_STOP * 4) inertia = null;
        }
      }
      drag = null;
      moving = !!inertia || spin;
      if (!moving) settleLabels();
      if (wasClick) call('onClick', pick(e.clientX, e.clientY), e);
      if (e.pointerType !== 'touch' && !cancelled) { hover.x = e.clientX; hover.y = e.clientY; hover.inside = true; hover.pending = true; }
      requestRender();
    }
    on(canvas, 'pointerup', e => endPointer(e, false));
    on(canvas, 'pointercancel', e => endPointer(e, true));
    on(canvas, 'pointerleave', e => {
      if (pointers.has(e.pointerId)) return;              // captured drag: still ours
      hover.inside = false; hover.pending = false;
      clearHover();
    });
    on(canvas, 'contextmenu', e => e.preventDefault());   // right-drag pans
    on(canvas, 'wheel', e => {
      e.preventDefault();
      const dy = e.deltaY * (e.deltaMode === 1 ? 16 : e.deltaMode === 2 ? (cssH || 400) : 1);
      zoomBy(Math.exp(dy * 0.0015));
      settleLabels();
    }, { passive: false });
    on(canvas, 'keydown', e => {                          // only reaches us while the canvas has focus
      if (e.altKey || e.ctrlKey || e.metaKey) return;
      const step = (e.shiftKey ? 1 : 10) * Math.PI / 180;
      let handled = true;
      switch (e.key) {
        case 'ArrowLeft': rotateWorld(0, 1, 0, -step); break;
        case 'ArrowRight': rotateWorld(0, 1, 0, step); break;
        case 'ArrowUp': rotateWorld(1, 0, 0, -step); break;
        case 'ArrowDown': rotateWorld(1, 0, 0, step); break;
        case '+': case '=': zoomBy(1 / 1.15); break;
        case '-': case '_': zoomBy(1.15); break;
        case '0': case 'Home': resetView(); break;
        case ' ': case 'Spacebar': setSpin(!spin); call('onSpinChange', spin); break;
        default: handled = false;
      }
      if (!handled) return;
      e.preventDefault();
      settleLabels();
      requestRender();
    });

    /* ---------- size, context loss, export, dispose ---------- */

    function resize() {
      if (disposed) return;
      const w = wrap.clientWidth, h = wrap.clientHeight;
      if (!w || !h) { unseen = true; return; }           // hidden tab: keep the last size; frames stop (frame())
      // shown again: draw what changed meanwhile and let spin or glide carry on (before the same-size return below,
      // since a box that was hidden and shown again usually has its old size)
      if (unseen) { unseen = false; requestRender(); }
      const pr = Math.min(window.devicePixelRatio || 1, 2);
      if (w === cssW && h === cssH && pr === renderer.getPixelRatio()) return;
      cssW = w; cssH = h;
      renderer.setPixelRatio(pr);
      renderer.setSize(w, h, false);
      overlay.width = Math.round(w * pr); overlay.height = Math.round(h * pr);
      if (framedWithoutSize && !userMoved && X) resetView();   // it was framed while hidden: frame it again
      updateCamera();
      labelVisCache = null;
      renderNow();                                       // draw at once, so a resize never shows a blank frame
    }
    const ro = typeof ResizeObserver === 'function' ? new ResizeObserver(() => resize()) : null;
    if (ro) ro.observe(wrap); else on(window, 'resize', resize);
    on(document, 'visibilitychange', () => { if (!document.hidden && !disposed) { unseen = false; requestRender(); } });

    on(canvas, 'webglcontextlost', e => {
      e.preventDefault();                                // lets the browser restore it
      lost = true;
      cancelAnimationFrame(raf); raf = 0;
      call('onContextLost');
    });
    on(canvas, 'webglcontextrestored', () => {
      lost = false;
      // Rebuild the scene from our own data. Objects from before the loss are dropped, never disposed: disposing
      // them makes three delete buffers that belonged to the lost context, which WebGL reports as INVALID_OPERATION.
      for (const k of ['atoms', 'bonds', 'extras', 'lines', 'dashes', 'measure']) { if (model[k]) inner.remove(model[k]); model[k] = null; }
      for (const k in halos) {
        const H = halos[k];
        for (const o of [H.sphere, H.cyl]) if (o) inner.remove(o);
        H.sphere = H.cyl = null; H.capS = H.capC = 0;
      }
      if (surf) { inner.remove(surf.pre, surf.mesh); }
      makeResources();
      rebuildModel();                                    // also halos and measurement lines
      if (surf) surfaceMeshes(surf);
      requestRender();
      call('onContextRestored');
    });

    // exportPNG(scale, {legend: {title, minLabel, maxLabel, midLabel?, note}, caption}): the picture with its labels;
    // with a legend, a red–white–blue bar (Mol3DView.espColor at −1, 0, +1 of the clamp) with the caller's labels
    // under its ends and the note under it; a caption (string or lines) is text only. They go into the bottom-left
    // corner when that corner shows only background, else into a strip added below the picture, so they never
    // cover the molecule, a label or a measurement. Without them the image is exactly the view.
    function exportPNG(scale = 2, o = {}) {
      return new Promise((resolve, reject) => {
        if (disposed) return reject(new Error('the 3D view was closed'));
        if (lost) return reject(new Error('the WebGL context is lost; try again in a moment'));
        if (!cssW || !cssH) return reject(new Error('the 3D view is not visible'));
        const gl = renderer.getContext(), maxRB = gl.getParameter(gl.MAX_RENDERBUFFER_SIZE), dims = gl.getParameter(gl.MAX_VIEWPORT_DIMS);
        const k = Math.max(0.5, Math.min(+scale || 2, maxRB / cssW, maxRB / cssH, dims[0] / cssW, dims[1] / cssH, 8));
        const oldPR = renderer.getPixelRatio();
        let out;
        const transparent = !!(o && o.transparent), keepBg = scene.background;
        try {
          renderer.setPixelRatio(k);
          if (transparent) { scene.background = null; renderer.setClearColor(0x000000, 0); }
          renderer.render(scene, camera);
          out = document.createElement('canvas');
          out.width = canvas.width; out.height = canvas.height;
          const ctx = out.getContext('2d');
          ctx.drawImage(canvas, 0, 0);
          const saved = labelVisCache;
          drawOverlay(ctx, out.width / cssW, true);      // the labels, drawn again at the export scale
          labelVisCache = saved;
          exportLabels = Object.assign({}, labelStats);
          const legend = o && o.legend && typeof o.legend === 'object' ? o.legend : null;
          let caption = o && o.caption != null && o.caption !== '' ? [].concat(o.caption).map(String) : null;
          // labels this image leaves out (overlaps; H labels with many heavy atoms) are said in it, as on screen: an
          // unlabelled atom must not read as one without a value
          const left = leftOutSentence(exportLabels);
          if (left) caption = (caption || []).concat(left);
          if (legend || caption) out = withLegend(out, out.width / cssW, legend, caption, transparent);
        } catch (e) { reject(e); return; } finally {
          if (transparent) { scene.background = keepBg; renderer.setClearColor(0xffffff, 1); }
          renderer.setPixelRatio(oldPR);
          renderNow();
        }
        out.toBlob(b => (b ? resolve(b) : reject(new Error('PNG encoding failed'))), 'image/png');
      });
    }

    // the export's own counts (its labels are placed again at export time, without the on-screen note's corner)
    let exportLabels = null;
    function leftOutSentence(st) {
      const parts = [];
      if (st && st.overlapping) parts.push(`${st.overlapping} label${st.overlapping === 1 ? '' : 's'} left out to avoid overlaps`);
      if (st && st.hidden) parts.push(`${st.hidden} label${st.hidden === 1 ? '' : 's'} left out because ${st.hidden === 1 ? 'its atom is' : 'their atoms are'} hidden behind other atoms`);
      if (st && st.hOnHover) parts.push(`H labels left out (more than ${H_LABELS_MAX_HEAVY} heavy atoms)`);
      return parts.length ? 'Not every atom is labelled: ' + parts.join('; ') + '.' : '';
    }

    const LEGEND_TITLE = '600 12px system-ui, "Segoe UI", sans-serif', LEGEND_TEXT = '11px system-ui, "Segoe UI", sans-serif';
    // words of text in lines no wider than maxW (CSS px) in the context's font (a longer single word keeps its line)
    function wrapLines(ctx, text, maxW) {
      const lines = [];
      for (const para of String(text).split('\n')) {
        let line = '';
        for (const w of para.split(/\s+/).filter(Boolean)) {
          const t = line ? line + ' ' + w : w;
          if (line && ctx.measureText(t).width > maxW) { lines.push(line); line = w; } else line = t;
        }
        if (line) lines.push(line);
      }
      return lines;
    }
    // the legend block laid out in CSS px: {w, h, draw(ctx, x, y)}
    function legendBlock(ctx, L, caption, maxW) {
      const str = v => (v == null ? '' : String(v));
      const fg = darkBg ? '#f2f2f2' : '#111', muted = darkBg ? '#c8c8c8' : '#444';
      const rows = [];                                   // {kind:'text', font, color, text, x, y, align} | {kind:'bar', y, w}
      let y = 0, w = 0;
      const text = (font, color, lines) => {
        ctx.font = font;
        for (const t of lines) { rows.push({ kind: 'text', font, color, text: t, x: 0, y: y + 11, align: 'left' }); w = Math.max(w, ctx.measureText(t).width); y += 14; }
      };
      if (L) {
        ctx.font = LEGEND_TITLE;
        const title = str(L.title) ? wrapLines(ctx, str(L.title), maxW) : [];
        text(LEGEND_TITLE, fg, title);
        if (title.length) y += 3;
        ctx.font = LEGEND_TEXT;
        const minL = str(L.minLabel), maxL = str(L.maxLabel), midL = str(L.midLabel);
        const wMin = ctx.measureText(minL).width, wMax = ctx.measureText(maxL).width, wMid = midL ? ctx.measureText(midL).width : 0;
        const need = wMin + wMax + (midL ? wMid + 24 : 12);
        const barW = Math.min(maxW, Math.max(180, need));
        rows.push({ kind: 'bar', y, w: barW });
        y += 12 + 3;
        w = Math.max(w, barW);
        if (need <= barW) {                              // both ends' labels (and the middle one) on one line
          if (minL) rows.push({ kind: 'text', font: LEGEND_TEXT, color: fg, text: minL, x: 0, y: y + 11, align: 'left' });
          if (midL) rows.push({ kind: 'text', font: LEGEND_TEXT, color: fg, text: midL, x: barW / 2, y: y + 11, align: 'center' });
          if (maxL) rows.push({ kind: 'text', font: LEGEND_TEXT, color: fg, text: maxL, x: barW, y: y + 11, align: 'right' });
          if (minL || midL || maxL) y += 14;
        } else {                                         // too long for one line: the left end's, then the right end's
          if (minL) { rows.push({ kind: 'text', font: LEGEND_TEXT, color: fg, text: '← ' + minL, x: 0, y: y + 11, align: 'left' }); y += 14; }
          if (maxL) { rows.push({ kind: 'text', font: LEGEND_TEXT, color: fg, text: maxL + ' →', x: barW, y: y + 11, align: 'right' }); y += 14; }
          w = Math.max(w, wMin + 16, wMax + 16);
        }
        const note = str(L.note);
        if (note) { y += 3; ctx.font = LEGEND_TEXT; text(LEGEND_TEXT, muted, wrapLines(ctx, note, maxW)); }
      }
      if (caption && caption.length) {
        if (y) y += 4;
        ctx.font = LEGEND_TEXT;
        text(LEGEND_TEXT, muted, caption.flatMap(c => wrapLines(ctx, c, maxW)));
      }
      return {
        w: Math.ceil(w), h: Math.ceil(y),
        draw(g, x0, y0) {
          for (const r of rows) {
            if (r.kind === 'bar') {
              const grad = g.createLinearGradient(x0, 0, x0 + r.w, 0);
              [-1, -0.5, 0, 0.5, 1].forEach(v => grad.addColorStop((v + 1) / 2, espColor(v, 1)));
              g.fillStyle = grad;
              g.fillRect(x0, y0 + r.y, r.w, 12);
              g.strokeStyle = darkBg ? '#888' : '#999'; g.lineWidth = 1;
              g.strokeRect(x0 + 0.5, y0 + r.y + 0.5, r.w - 1, 11);
              continue;
            }
            g.font = r.font; g.fillStyle = r.color; g.textAlign = r.align; g.textBaseline = 'alphabetic';
            g.fillText(r.text, x0 + r.x, y0 + r.y);
          }
        },
      };
    }
    // the export image with the legend in its bottom-left corner, or in a strip added below it
    function withLegend(img, s, L, caption, transparent) {
      const PAD = 10, W = img.width / s, H = img.height / s;
      const ctx = img.getContext('2d');
      ctx.setTransform(1, 0, 0, 1, 0, 0);
      const B = legendBlock(ctx, L, caption, Math.max(120, Math.min(W - 2 * PAD, 440)));
      // is the corner (the block plus its margin) background only? then the legend covers nothing there
      const px0 = 0, py0 = Math.max(0, Math.floor((H - B.h - 2 * PAD) * s));
      const pw = Math.min(img.width, Math.ceil((B.w + 2 * PAD) * s)), ph = img.height - py0;
      let free = H - B.h - 2 * PAD >= 0 && pw > 0 && ph > 0;
      if (free) {
        const probe = document.createElement('canvas');
        probe.width = probe.height = 1;
        const pg = probe.getContext('2d');
        pg.fillStyle = bgColor; pg.fillRect(0, 0, 1, 1);
        const bg = pg.getImageData(0, 0, 1, 1).data, d = ctx.getImageData(px0, py0, pw, ph).data;
        for (let q = 0; q < d.length; q += 4) {
          if (transparent) { if (d[q + 3] > 8) { free = false; break; } continue; }   // transparent: nothing drawn there
          if (Math.abs(d[q] - bg[0]) > 8 || Math.abs(d[q + 1] - bg[1]) > 8 || Math.abs(d[q + 2] - bg[2]) > 8) { free = false; break; }
        }
      }
      if (free) {
        ctx.setTransform(s, 0, 0, s, 0, 0);
        B.draw(ctx, PAD, H - PAD - B.h);
        return img;
      }
      const ext = document.createElement('canvas');
      ext.width = img.width; ext.height = img.height + Math.round((B.h + 2 * PAD) * s);
      const g = ext.getContext('2d');
      if (!transparent) { g.fillStyle = bgColor; g.fillRect(0, 0, ext.width, ext.height); }
      g.drawImage(img, 0, 0);
      g.setTransform(s, 0, 0, s, 0, 0);
      B.draw(g, PAD, H + PAD);
      return ext;
    }

    function dispose() {
      if (disposed) return;
      disposed = true;
      surfGen++;
      cancelAnimationFrame(raf); raf = 0;
      clearTimeout(settleTimer);
      if (ro) ro.disconnect();
      for (const [t, type, fn, o] of listeners) t.removeEventListener(type, fn, o);
      clearModel();
      disposeObj(model.measure);
      removeSurface();
      for (const k in halos) { disposeObj(halos[k].sphere); disposeObj(halos[k].cyl); }
      Object.values(geo).forEach(g => g.dispose());
      [mat.solid, mat.line, mat.dash, mat.measure, mat.surfPre, mat.electron, ...Object.values(mat.halo)].forEach(m => m.dispose());
      renderer.renderLists.dispose();
      wrap.remove();
      // keep this renderer for the next view (browsers cap live WebGL contexts); release an older spare
      const release = r => { r.dispose(); try { r.forceContextLoss(); } catch (e) { /* already gone */ } };
      if (spareRenderer) release(spareRenderer);
      spareRenderer = null;
      if (lost) release(renderer); else spareRenderer = renderer;
    }

    // style and H change what has to fit: an untouched view is framed again
    function setStyle(s) {
      if (!STYLES.includes(s)) throw new Error('Mol3DView.setStyle: one of ' + STYLES.join(', '));
      style = s; rebuildModel(); if (X && !userMoved) resetView(); hover.pending = hover.inside; requestRender();
    }
    function setShowH(b) { showH = !!b; rebuildModel(); if (X && !userMoved) resetView(); hover.pending = hover.inside; requestRender(); }

    updateCamera();
    resize();

    // after dispose() every method is a harmless no-op (setSurface resolves null, exportPNG rejects, pick gives null)
    const live = f => (...a) => (disposed ? undefined : f(...a));
    return {
      setMolecule: live(setMolecule), setCoords: live(setCoords), setStyle: live(setStyle), setShowH: live(setShowH),
      setLabels: live(setLabels), setElectrons: live(setElectrons), setSpin: live(setSpin), resetView: live(resetView), setOrthographic: live(setOrthographic), setPickMode: live(setPickMode),
      setSurface: (...a) => (disposed ? Promise.resolve(null) : setSurface(...a)),
      highlight: live(highlight), setMeasurements: live(setMeasurements), exportPNG, resize: live(resize), dispose,
      pick: (...a) => (disposed ? null : pick(...a)),
      get canvas() { return canvas; },
      stats: () => Object.assign({ atoms: mol ? mol.N : 0, bonds: mol ? mol.M : 0, style, showH, spin, orthographic: orthoOn,
        instances: { atoms: model.atoms ? model.atoms.count : 0, bonds: model.bonds ? model.bonds.count : 0, extras: model.extras ? model.extras.count : 0 },
        drawCalls: renderer.info.render.calls, triangles: renderer.info.render.triangles,
        surface: surf ? surf.info : null, labels: Object.assign({}, labelStats) }, stats),
      // for tests/molview/view3d_test.html: projection of a molecule-coordinate point, the meshes, the surface data
      _debug: {
        project: (x, y, z) => {
          molGroup.updateMatrixWorld(true);
          const pt = project(x, y, z, {}), r = canvas.getBoundingClientRect();
          return { clientX: r.left + pt.x, clientY: r.top + pt.y, ndcZ: pt.z };
        },
        meshes: () => ({ atoms: model.atoms, atomOf: model.atomOf, bonds: model.bonds, bondOf: model.bondOf, extras: model.extras, extraOf: model.extraOf, camera, inner }),
        /* the lone-pair directions the electrons labels use, atom by atom (unit vectors), for tests */
        lonePairDirs: (i, k) => (mol && X ? lonePairDirs(i, k).map(v => [v.x, v.y, v.z]) : null),
        lonePairDots: () => (model.lonePairs ? model.lonePairs.count : 0),
        surface: () => (surf ? { positions: surf.data.positions, normals: surf.data.normals, index: surf.data.index, owner: surf.data.owner, V: surf.V, h: surf.data.h } : null),
        radii: () => ({ draw: atomR && Array.from(atomR), pick: atomPickR && Array.from(atomPickR), vis: atomVis && Array.from(atomVis) }),
        labelVisibility: () => (labels && mol && atomVis ? Array.from(labelVisibility()) : null),
        lastPick: () => lastPick,
        exportLabels: () => (exportLabels ? Object.assign({}, exportLabels) : null),   // the last export's label counts
        labelsDrawn: () => lastLabels.map(l => ({ index: l.index, text: l.text, box: l.box.slice() })),
        measuresDrawn: () => lastMeasures.map(m => [m.x, m.y, m.x + m.w, m.y + m.h]),
        // one frame drawn at once and waited for (a 1-pixel read makes the GPU finish): ms, for timing tests
        timeFrame: (turn = 0) => {
          if (turn) rotateWorld(0, 1, 0, turn);
          const gl = renderer.getContext(), px = new Uint8Array(4), t0 = performance.now();
          molGroup.updateMatrixWorld(true);
          renderNow();
          gl.readPixels(0, 0, 1, 1, gl.RGBA, gl.UNSIGNED_BYTE, px);
          return performance.now() - t0;
        },
        cameraLocal: () => {                             // camera position in molecule coordinates
          molGroup.updateMatrixWorld(true);
          return new THREE.Vector3().copy(camera.position).applyMatrix4(_inv.copy(inner.matrixWorld).invert());
        },
        pxPerA: (x, y, z) => {                           // CSS px per Å across the view at a point
          molGroup.updateMatrixWorld(true);
          const up = new THREE.Vector3(0, 1, 0).transformDirection(_inv.copy(inner.matrixWorld).invert());
          const p0 = project(x, y, z, {}), p1 = project(x + up.x, y + up.y, z + up.z, {});
          return Math.hypot(p1.x - p0.x, p1.y - p0.y);
        },
        renderer,
      },
    };
  }

  /* Marching-cubes tables: three@0.185.1 examples/jsm/objects/MarchingCubes.js (MIT; see the header), from Paul
     Bourke / Cory Gene Bloyd. Corner i of a cube and edge e use Bourke's numbering. */
  const EDGE_TABLE = new Int32Array([
    0x0,0x109,0x203,0x30a,0x406,0x50f,0x605,0x70c,0x80c,0x905,0xa0f,0xb06,0xc0a,0xd03,0xe09,0xf00,
    0x190,0x99,0x393,0x29a,0x596,0x49f,0x795,0x69c,0x99c,0x895,0xb9f,0xa96,0xd9a,0xc93,0xf99,0xe90,
    0x230,0x339,0x33,0x13a,0x636,0x73f,0x435,0x53c,0xa3c,0xb35,0x83f,0x936,0xe3a,0xf33,0xc39,0xd30,
    0x3a0,0x2a9,0x1a3,0xaa,0x7a6,0x6af,0x5a5,0x4ac,0xbac,0xaa5,0x9af,0x8a6,0xfaa,0xea3,0xda9,0xca0,
    0x460,0x569,0x663,0x76a,0x66,0x16f,0x265,0x36c,0xc6c,0xd65,0xe6f,0xf66,0x86a,0x963,0xa69,0xb60,
    0x5f0,0x4f9,0x7f3,0x6fa,0x1f6,0xff,0x3f5,0x2fc,0xdfc,0xcf5,0xfff,0xef6,0x9fa,0x8f3,0xbf9,0xaf0,
    0x650,0x759,0x453,0x55a,0x256,0x35f,0x55,0x15c,0xe5c,0xf55,0xc5f,0xd56,0xa5a,0xb53,0x859,0x950,
    0x7c0,0x6c9,0x5c3,0x4ca,0x3c6,0x2cf,0x1c5,0xcc,0xfcc,0xec5,0xdcf,0xcc6,0xbca,0xac3,0x9c9,0x8c0,
    0x8c0,0x9c9,0xac3,0xbca,0xcc6,0xdcf,0xec5,0xfcc,0xcc,0x1c5,0x2cf,0x3c6,0x4ca,0x5c3,0x6c9,0x7c0,
    0x950,0x859,0xb53,0xa5a,0xd56,0xc5f,0xf55,0xe5c,0x15c,0x55,0x35f,0x256,0x55a,0x453,0x759,0x650,
    0xaf0,0xbf9,0x8f3,0x9fa,0xef6,0xfff,0xcf5,0xdfc,0x2fc,0x3f5,0xff,0x1f6,0x6fa,0x7f3,0x4f9,0x5f0,
    0xb60,0xa69,0x963,0x86a,0xf66,0xe6f,0xd65,0xc6c,0x36c,0x265,0x16f,0x66,0x76a,0x663,0x569,0x460,
    0xca0,0xda9,0xea3,0xfaa,0x8a6,0x9af,0xaa5,0xbac,0x4ac,0x5a5,0x6af,0x7a6,0xaa,0x1a3,0x2a9,0x3a0,
    0xd30,0xc39,0xf33,0xe3a,0x936,0x83f,0xb35,0xa3c,0x53c,0x435,0x73f,0x636,0x13a,0x33,0x339,0x230,
    0xe90,0xf99,0xc93,0xd9a,0xa96,0xb9f,0x895,0x99c,0x69c,0x795,0x49f,0x596,0x29a,0x393,0x99,0x190,
    0xf00,0xe09,0xd03,0xc0a,0xb06,0xa0f,0x905,0x80c,0x70c,0x605,0x50f,0x406,0x30a,0x203,0x109,0x0,
  ]);
  const TRI_TABLE = new Int32Array([
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,8,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,1,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,8,3,9,8,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,8,3,1,2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    9,2,10,0,2,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    2,8,3,2,10,8,10,9,8,-1,-1,-1,-1,-1,-1,-1,
    3,11,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,11,2,8,11,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,9,0,2,3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,11,2,1,9,11,9,8,11,-1,-1,-1,-1,-1,-1,-1,
    3,10,1,11,10,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,10,1,0,8,10,8,11,10,-1,-1,-1,-1,-1,-1,-1,
    3,9,0,3,11,9,11,10,9,-1,-1,-1,-1,-1,-1,-1,
    9,8,10,10,8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    4,7,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    4,3,0,7,3,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,1,9,8,4,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    4,1,9,4,7,1,7,3,1,-1,-1,-1,-1,-1,-1,-1,
    1,2,10,8,4,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    3,4,7,3,0,4,1,2,10,-1,-1,-1,-1,-1,-1,-1,
    9,2,10,9,0,2,8,4,7,-1,-1,-1,-1,-1,-1,-1,
    2,10,9,2,9,7,2,7,3,7,9,4,-1,-1,-1,-1,
    8,4,7,3,11,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    11,4,7,11,2,4,2,0,4,-1,-1,-1,-1,-1,-1,-1,
    9,0,1,8,4,7,2,3,11,-1,-1,-1,-1,-1,-1,-1,
    4,7,11,9,4,11,9,11,2,9,2,1,-1,-1,-1,-1,
    3,10,1,3,11,10,7,8,4,-1,-1,-1,-1,-1,-1,-1,
    1,11,10,1,4,11,1,0,4,7,11,4,-1,-1,-1,-1,
    4,7,8,9,0,11,9,11,10,11,0,3,-1,-1,-1,-1,
    4,7,11,4,11,9,9,11,10,-1,-1,-1,-1,-1,-1,-1,
    9,5,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    9,5,4,0,8,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,5,4,1,5,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    8,5,4,8,3,5,3,1,5,-1,-1,-1,-1,-1,-1,-1,
    1,2,10,9,5,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    3,0,8,1,2,10,4,9,5,-1,-1,-1,-1,-1,-1,-1,
    5,2,10,5,4,2,4,0,2,-1,-1,-1,-1,-1,-1,-1,
    2,10,5,3,2,5,3,5,4,3,4,8,-1,-1,-1,-1,
    9,5,4,2,3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,11,2,0,8,11,4,9,5,-1,-1,-1,-1,-1,-1,-1,
    0,5,4,0,1,5,2,3,11,-1,-1,-1,-1,-1,-1,-1,
    2,1,5,2,5,8,2,8,11,4,8,5,-1,-1,-1,-1,
    10,3,11,10,1,3,9,5,4,-1,-1,-1,-1,-1,-1,-1,
    4,9,5,0,8,1,8,10,1,8,11,10,-1,-1,-1,-1,
    5,4,0,5,0,11,5,11,10,11,0,3,-1,-1,-1,-1,
    5,4,8,5,8,10,10,8,11,-1,-1,-1,-1,-1,-1,-1,
    9,7,8,5,7,9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    9,3,0,9,5,3,5,7,3,-1,-1,-1,-1,-1,-1,-1,
    0,7,8,0,1,7,1,5,7,-1,-1,-1,-1,-1,-1,-1,
    1,5,3,3,5,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    9,7,8,9,5,7,10,1,2,-1,-1,-1,-1,-1,-1,-1,
    10,1,2,9,5,0,5,3,0,5,7,3,-1,-1,-1,-1,
    8,0,2,8,2,5,8,5,7,10,5,2,-1,-1,-1,-1,
    2,10,5,2,5,3,3,5,7,-1,-1,-1,-1,-1,-1,-1,
    7,9,5,7,8,9,3,11,2,-1,-1,-1,-1,-1,-1,-1,
    9,5,7,9,7,2,9,2,0,2,7,11,-1,-1,-1,-1,
    2,3,11,0,1,8,1,7,8,1,5,7,-1,-1,-1,-1,
    11,2,1,11,1,7,7,1,5,-1,-1,-1,-1,-1,-1,-1,
    9,5,8,8,5,7,10,1,3,10,3,11,-1,-1,-1,-1,
    5,7,0,5,0,9,7,11,0,1,0,10,11,10,0,-1,
    11,10,0,11,0,3,10,5,0,8,0,7,5,7,0,-1,
    11,10,5,7,11,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    10,6,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,8,3,5,10,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    9,0,1,5,10,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,8,3,1,9,8,5,10,6,-1,-1,-1,-1,-1,-1,-1,
    1,6,5,2,6,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,6,5,1,2,6,3,0,8,-1,-1,-1,-1,-1,-1,-1,
    9,6,5,9,0,6,0,2,6,-1,-1,-1,-1,-1,-1,-1,
    5,9,8,5,8,2,5,2,6,3,2,8,-1,-1,-1,-1,
    2,3,11,10,6,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    11,0,8,11,2,0,10,6,5,-1,-1,-1,-1,-1,-1,-1,
    0,1,9,2,3,11,5,10,6,-1,-1,-1,-1,-1,-1,-1,
    5,10,6,1,9,2,9,11,2,9,8,11,-1,-1,-1,-1,
    6,3,11,6,5,3,5,1,3,-1,-1,-1,-1,-1,-1,-1,
    0,8,11,0,11,5,0,5,1,5,11,6,-1,-1,-1,-1,
    3,11,6,0,3,6,0,6,5,0,5,9,-1,-1,-1,-1,
    6,5,9,6,9,11,11,9,8,-1,-1,-1,-1,-1,-1,-1,
    5,10,6,4,7,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    4,3,0,4,7,3,6,5,10,-1,-1,-1,-1,-1,-1,-1,
    1,9,0,5,10,6,8,4,7,-1,-1,-1,-1,-1,-1,-1,
    10,6,5,1,9,7,1,7,3,7,9,4,-1,-1,-1,-1,
    6,1,2,6,5,1,4,7,8,-1,-1,-1,-1,-1,-1,-1,
    1,2,5,5,2,6,3,0,4,3,4,7,-1,-1,-1,-1,
    8,4,7,9,0,5,0,6,5,0,2,6,-1,-1,-1,-1,
    7,3,9,7,9,4,3,2,9,5,9,6,2,6,9,-1,
    3,11,2,7,8,4,10,6,5,-1,-1,-1,-1,-1,-1,-1,
    5,10,6,4,7,2,4,2,0,2,7,11,-1,-1,-1,-1,
    0,1,9,4,7,8,2,3,11,5,10,6,-1,-1,-1,-1,
    9,2,1,9,11,2,9,4,11,7,11,4,5,10,6,-1,
    8,4,7,3,11,5,3,5,1,5,11,6,-1,-1,-1,-1,
    5,1,11,5,11,6,1,0,11,7,11,4,0,4,11,-1,
    0,5,9,0,6,5,0,3,6,11,6,3,8,4,7,-1,
    6,5,9,6,9,11,4,7,9,7,11,9,-1,-1,-1,-1,
    10,4,9,6,4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    4,10,6,4,9,10,0,8,3,-1,-1,-1,-1,-1,-1,-1,
    10,0,1,10,6,0,6,4,0,-1,-1,-1,-1,-1,-1,-1,
    8,3,1,8,1,6,8,6,4,6,1,10,-1,-1,-1,-1,
    1,4,9,1,2,4,2,6,4,-1,-1,-1,-1,-1,-1,-1,
    3,0,8,1,2,9,2,4,9,2,6,4,-1,-1,-1,-1,
    0,2,4,4,2,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    8,3,2,8,2,4,4,2,6,-1,-1,-1,-1,-1,-1,-1,
    10,4,9,10,6,4,11,2,3,-1,-1,-1,-1,-1,-1,-1,
    0,8,2,2,8,11,4,9,10,4,10,6,-1,-1,-1,-1,
    3,11,2,0,1,6,0,6,4,6,1,10,-1,-1,-1,-1,
    6,4,1,6,1,10,4,8,1,2,1,11,8,11,1,-1,
    9,6,4,9,3,6,9,1,3,11,6,3,-1,-1,-1,-1,
    8,11,1,8,1,0,11,6,1,9,1,4,6,4,1,-1,
    3,11,6,3,6,0,0,6,4,-1,-1,-1,-1,-1,-1,-1,
    6,4,8,11,6,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    7,10,6,7,8,10,8,9,10,-1,-1,-1,-1,-1,-1,-1,
    0,7,3,0,10,7,0,9,10,6,7,10,-1,-1,-1,-1,
    10,6,7,1,10,7,1,7,8,1,8,0,-1,-1,-1,-1,
    10,6,7,10,7,1,1,7,3,-1,-1,-1,-1,-1,-1,-1,
    1,2,6,1,6,8,1,8,9,8,6,7,-1,-1,-1,-1,
    2,6,9,2,9,1,6,7,9,0,9,3,7,3,9,-1,
    7,8,0,7,0,6,6,0,2,-1,-1,-1,-1,-1,-1,-1,
    7,3,2,6,7,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    2,3,11,10,6,8,10,8,9,8,6,7,-1,-1,-1,-1,
    2,0,7,2,7,11,0,9,7,6,7,10,9,10,7,-1,
    1,8,0,1,7,8,1,10,7,6,7,10,2,3,11,-1,
    11,2,1,11,1,7,10,6,1,6,7,1,-1,-1,-1,-1,
    8,9,6,8,6,7,9,1,6,11,6,3,1,3,6,-1,
    0,9,1,11,6,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    7,8,0,7,0,6,3,11,0,11,6,0,-1,-1,-1,-1,
    7,11,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    7,6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    3,0,8,11,7,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,1,9,11,7,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    8,1,9,8,3,1,11,7,6,-1,-1,-1,-1,-1,-1,-1,
    10,1,2,6,11,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,2,10,3,0,8,6,11,7,-1,-1,-1,-1,-1,-1,-1,
    2,9,0,2,10,9,6,11,7,-1,-1,-1,-1,-1,-1,-1,
    6,11,7,2,10,3,10,8,3,10,9,8,-1,-1,-1,-1,
    7,2,3,6,2,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    7,0,8,7,6,0,6,2,0,-1,-1,-1,-1,-1,-1,-1,
    2,7,6,2,3,7,0,1,9,-1,-1,-1,-1,-1,-1,-1,
    1,6,2,1,8,6,1,9,8,8,7,6,-1,-1,-1,-1,
    10,7,6,10,1,7,1,3,7,-1,-1,-1,-1,-1,-1,-1,
    10,7,6,1,7,10,1,8,7,1,0,8,-1,-1,-1,-1,
    0,3,7,0,7,10,0,10,9,6,10,7,-1,-1,-1,-1,
    7,6,10,7,10,8,8,10,9,-1,-1,-1,-1,-1,-1,-1,
    6,8,4,11,8,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    3,6,11,3,0,6,0,4,6,-1,-1,-1,-1,-1,-1,-1,
    8,6,11,8,4,6,9,0,1,-1,-1,-1,-1,-1,-1,-1,
    9,4,6,9,6,3,9,3,1,11,3,6,-1,-1,-1,-1,
    6,8,4,6,11,8,2,10,1,-1,-1,-1,-1,-1,-1,-1,
    1,2,10,3,0,11,0,6,11,0,4,6,-1,-1,-1,-1,
    4,11,8,4,6,11,0,2,9,2,10,9,-1,-1,-1,-1,
    10,9,3,10,3,2,9,4,3,11,3,6,4,6,3,-1,
    8,2,3,8,4,2,4,6,2,-1,-1,-1,-1,-1,-1,-1,
    0,4,2,4,6,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,9,0,2,3,4,2,4,6,4,3,8,-1,-1,-1,-1,
    1,9,4,1,4,2,2,4,6,-1,-1,-1,-1,-1,-1,-1,
    8,1,3,8,6,1,8,4,6,6,10,1,-1,-1,-1,-1,
    10,1,0,10,0,6,6,0,4,-1,-1,-1,-1,-1,-1,-1,
    4,6,3,4,3,8,6,10,3,0,3,9,10,9,3,-1,
    10,9,4,6,10,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    4,9,5,7,6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,8,3,4,9,5,11,7,6,-1,-1,-1,-1,-1,-1,-1,
    5,0,1,5,4,0,7,6,11,-1,-1,-1,-1,-1,-1,-1,
    11,7,6,8,3,4,3,5,4,3,1,5,-1,-1,-1,-1,
    9,5,4,10,1,2,7,6,11,-1,-1,-1,-1,-1,-1,-1,
    6,11,7,1,2,10,0,8,3,4,9,5,-1,-1,-1,-1,
    7,6,11,5,4,10,4,2,10,4,0,2,-1,-1,-1,-1,
    3,4,8,3,5,4,3,2,5,10,5,2,11,7,6,-1,
    7,2,3,7,6,2,5,4,9,-1,-1,-1,-1,-1,-1,-1,
    9,5,4,0,8,6,0,6,2,6,8,7,-1,-1,-1,-1,
    3,6,2,3,7,6,1,5,0,5,4,0,-1,-1,-1,-1,
    6,2,8,6,8,7,2,1,8,4,8,5,1,5,8,-1,
    9,5,4,10,1,6,1,7,6,1,3,7,-1,-1,-1,-1,
    1,6,10,1,7,6,1,0,7,8,7,0,9,5,4,-1,
    4,0,10,4,10,5,0,3,10,6,10,7,3,7,10,-1,
    7,6,10,7,10,8,5,4,10,4,8,10,-1,-1,-1,-1,
    6,9,5,6,11,9,11,8,9,-1,-1,-1,-1,-1,-1,-1,
    3,6,11,0,6,3,0,5,6,0,9,5,-1,-1,-1,-1,
    0,11,8,0,5,11,0,1,5,5,6,11,-1,-1,-1,-1,
    6,11,3,6,3,5,5,3,1,-1,-1,-1,-1,-1,-1,-1,
    1,2,10,9,5,11,9,11,8,11,5,6,-1,-1,-1,-1,
    0,11,3,0,6,11,0,9,6,5,6,9,1,2,10,-1,
    11,8,5,11,5,6,8,0,5,10,5,2,0,2,5,-1,
    6,11,3,6,3,5,2,10,3,10,5,3,-1,-1,-1,-1,
    5,8,9,5,2,8,5,6,2,3,8,2,-1,-1,-1,-1,
    9,5,6,9,6,0,0,6,2,-1,-1,-1,-1,-1,-1,-1,
    1,5,8,1,8,0,5,6,8,3,8,2,6,2,8,-1,
    1,5,6,2,1,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,3,6,1,6,10,3,8,6,5,6,9,8,9,6,-1,
    10,1,0,10,0,6,9,5,0,5,6,0,-1,-1,-1,-1,
    0,3,8,5,6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    10,5,6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    11,5,10,7,5,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    11,5,10,11,7,5,8,3,0,-1,-1,-1,-1,-1,-1,-1,
    5,11,7,5,10,11,1,9,0,-1,-1,-1,-1,-1,-1,-1,
    10,7,5,10,11,7,9,8,1,8,3,1,-1,-1,-1,-1,
    11,1,2,11,7,1,7,5,1,-1,-1,-1,-1,-1,-1,-1,
    0,8,3,1,2,7,1,7,5,7,2,11,-1,-1,-1,-1,
    9,7,5,9,2,7,9,0,2,2,11,7,-1,-1,-1,-1,
    7,5,2,7,2,11,5,9,2,3,2,8,9,8,2,-1,
    2,5,10,2,3,5,3,7,5,-1,-1,-1,-1,-1,-1,-1,
    8,2,0,8,5,2,8,7,5,10,2,5,-1,-1,-1,-1,
    9,0,1,5,10,3,5,3,7,3,10,2,-1,-1,-1,-1,
    9,8,2,9,2,1,8,7,2,10,2,5,7,5,2,-1,
    1,3,5,3,7,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,8,7,0,7,1,1,7,5,-1,-1,-1,-1,-1,-1,-1,
    9,0,3,9,3,5,5,3,7,-1,-1,-1,-1,-1,-1,-1,
    9,8,7,5,9,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    5,8,4,5,10,8,10,11,8,-1,-1,-1,-1,-1,-1,-1,
    5,0,4,5,11,0,5,10,11,11,3,0,-1,-1,-1,-1,
    0,1,9,8,4,10,8,10,11,10,4,5,-1,-1,-1,-1,
    10,11,4,10,4,5,11,3,4,9,4,1,3,1,4,-1,
    2,5,1,2,8,5,2,11,8,4,5,8,-1,-1,-1,-1,
    0,4,11,0,11,3,4,5,11,2,11,1,5,1,11,-1,
    0,2,5,0,5,9,2,11,5,4,5,8,11,8,5,-1,
    9,4,5,2,11,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    2,5,10,3,5,2,3,4,5,3,8,4,-1,-1,-1,-1,
    5,10,2,5,2,4,4,2,0,-1,-1,-1,-1,-1,-1,-1,
    3,10,2,3,5,10,3,8,5,4,5,8,0,1,9,-1,
    5,10,2,5,2,4,1,9,2,9,4,2,-1,-1,-1,-1,
    8,4,5,8,5,3,3,5,1,-1,-1,-1,-1,-1,-1,-1,
    0,4,5,1,0,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    8,4,5,8,5,3,9,0,5,0,3,5,-1,-1,-1,-1,
    9,4,5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    4,11,7,4,9,11,9,10,11,-1,-1,-1,-1,-1,-1,-1,
    0,8,3,4,9,7,9,11,7,9,10,11,-1,-1,-1,-1,
    1,10,11,1,11,4,1,4,0,7,4,11,-1,-1,-1,-1,
    3,1,4,3,4,8,1,10,4,7,4,11,10,11,4,-1,
    4,11,7,9,11,4,9,2,11,9,1,2,-1,-1,-1,-1,
    9,7,4,9,11,7,9,1,11,2,11,1,0,8,3,-1,
    11,7,4,11,4,2,2,4,0,-1,-1,-1,-1,-1,-1,-1,
    11,7,4,11,4,2,8,3,4,3,2,4,-1,-1,-1,-1,
    2,9,10,2,7,9,2,3,7,7,4,9,-1,-1,-1,-1,
    9,10,7,9,7,4,10,2,7,8,7,0,2,0,7,-1,
    3,7,10,3,10,2,7,4,10,1,10,0,4,0,10,-1,
    1,10,2,8,7,4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    4,9,1,4,1,7,7,1,3,-1,-1,-1,-1,-1,-1,-1,
    4,9,1,4,1,7,0,8,1,8,7,1,-1,-1,-1,-1,
    4,0,3,7,4,3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    4,8,7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    9,10,8,10,11,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    3,0,9,3,9,11,11,9,10,-1,-1,-1,-1,-1,-1,-1,
    0,1,10,0,10,8,8,10,11,-1,-1,-1,-1,-1,-1,-1,
    3,1,10,11,3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,2,11,1,11,9,9,11,8,-1,-1,-1,-1,-1,-1,-1,
    3,0,9,3,9,11,1,2,9,2,11,9,-1,-1,-1,-1,
    0,2,11,8,0,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    3,2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    2,3,8,2,8,10,10,8,9,-1,-1,-1,-1,-1,-1,-1,
    9,10,2,0,9,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    2,3,8,2,8,10,0,1,8,1,10,8,-1,-1,-1,-1,
    1,10,2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    1,3,8,9,1,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,9,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0,3,8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
  ]);

  window.Mol3DView = {
    create, espColor,
    _internals: { buildSurface, coulomb, percentile, principalAxes, eigenSym3, espRGB, JMOL, COULOMB, EDGE_TABLE, TRI_TABLE },
  };
})();
