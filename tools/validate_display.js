// Label validation, the drawing: runs the Molecule tab itself (2D view) on a sample of molecules in every Labels mode and
// measures what a reader sees: labels left out ("no free space"), labels nearer another atom than their own (which atom
// is it for?), labels inside a ring, labels far from their atom. In 'angles' mode it also takes the 3D model's
// coordinates, so tools/validate_labels.py can compare the bond angles with RDKit's MMFF94 geometry.
//   python tools/serve.py   (another terminal)
//   node tools/validate_display.js [http://localhost:8000/app.html] [--n 120] [--angles 40] [--out tools/data/label_display.json]
const fs = require('fs'), path = require('path');
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const args = process.argv.slice(2), opt = (k, d) => { const i = args.indexOf(k); return i >= 0 ? args[i + 1] : d; };
const base = args.find(a => /^https?:/.test(a)) || 'http://localhost:8000/app.html';
const N = +opt('--n', 120), NA = +opt('--angles', 40), out = opt('--out', path.join(__dirname, 'data', 'label_display.json'));
const MODES = ['index', 'cip', 'hybridization', 'formal', 'partial', 'oxidation', 'lonepairs', 'pka', 'steric', 'electronic', 'taft'];

// the sample: the textbook list, then every k-th single-component molecule of data/pka.js with 6-30 heavy atoms
function sample() {
  const T = JSON.parse(fs.readFileSync(path.join(__dirname, 'label_textbook.json'), 'utf8'));
  const s = fs.readFileSync(path.join(__dirname, '..', 'data', 'pka.js'), 'utf8');
  const all = Object.keys(JSON.parse(s.slice(s.indexOf('{'), s.lastIndexOf('}') + 1)).molecules)
    .filter(x => !x.includes('.')).filter(x => { const h = (x.match(/\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnops]/g) || []).length; return h >= 6 && h <= 30; });
  const k = Math.max(1, Math.floor(all.length / N));
  const pick = all.filter((_, i) => i % k === 0).slice(0, N);
  return [...new Set(T.cases.map(c => c.smiles).filter(x => x.length > 3).concat(pick))];
}

(async () => {
  const smis = sample();
  const b = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const p = await (await b.newContext({ viewport: { width: 1400, height: 900 } })).newPage();
  const errs = []; p.on('pageerror', e => errs.push(e.message));
  await p.goto(base + '#mol', { waitUntil: 'load' });
  await p.waitForFunction(() => window.MolView, null, { timeout: 60000 });
  const res = {};
  let k = 0;
  for (const smi of smis) {
    k++;
    const ok = await p.evaluate(async smi => (await MolView.show(smi, { view: '2d' })).ok, smi);
    if (!ok) { res[smi] = { error: 'not shown' }; continue; }
    await p.waitForTimeout(50);
    const r = { modes: {} };
    for (const md of MODES) {
      r.modes[md] = await p.evaluate(md => {
        const S = MolView._debug.state, v2 = MolView._debug.views().v2, rec = S.cs.record, n = rec.n;
        const sel = document.getElementById('molLabels'); sel.value = md; sel.dispatchEvent(new Event('change'));
        const pos = []; for (let i = 0; i < n; i++) pos.push(v2.clientPosition('atom', i));
        const bl = []; for (const bd of rec.bonds) if (bd.a < n && bd.b < n) bl.push(Math.hypot(pos[bd.a].x - pos[bd.b].x, pos[bd.a].y - pos[bd.b].y));
        bl.sort((x, y) => x - y);
        const BL = bl.length ? bl[bl.length >> 1] : 40;
        const rings = (rec.rings || []).filter(rg => rg.every(i => i < n)).map(rg => rg.map(i => pos[i]));
        const inside = (pt, poly) => { let c = false; for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) { const a = poly[i], b = poly[j]; if ((a.y > pt.y) !== (b.y > pt.y) && pt.x < (b.x - a.x) * (pt.y - a.y) / (b.y - a.y) + a.x) c = !c; } return c; };
        const labelBox = {};
        for (const t of document.querySelectorAll('#mol2d .mv-2d-labels text[data-atom]')) labelBox[+t.getAttribute('data-atom')] = t.getBoundingClientRect();
        const notes = [...document.querySelectorAll('#mol2d .mv-2d-notes text')].map(t => {
          const R = t.getBoundingClientRect(), c = { x: R.left + R.width / 2, y: R.top + R.height / 2 }, own = +t.getAttribute('data-atom');
          let near = own, dn = Infinity;
          pos.forEach((q, i) => { const d = Math.hypot(q.x - c.x, q.y - c.y); if (d < dn) { dn = d; near = i; } });
          const d0 = Math.hypot(pos[own].x - c.x, pos[own].y - c.y);
          // the gap between an atom and the label's nearest edge: a long label is read from the end next to its atom
          // to an atom as drawn: its written label (OH, NH2, Cl) when it has one, else its point
          const edge = (q, i) => { const L = labelBox[i];
            return L ? Math.hypot(Math.max(L.left - R.right, 0, R.left - L.right), Math.max(L.top - R.bottom, 0, R.top - L.bottom))
              : Math.hypot(Math.max(R.left - q.x, 0, q.x - R.right), Math.max(R.top - q.y, 0, q.y - R.bottom)); };
          const e0 = edge(pos[own], own), edgeNear = pos.findIndex((q, i) => i !== own && edge(q, i) < e0);
          return { atom: own, text: t.textContent, nearest: near, ownDist: +(d0 / BL).toFixed(2), nearDist: +(dn / BL).toFixed(2), inRing: rings.some(pg => inside(c, pg)),
            edgeGap: +(e0 / BL).toFixed(2), edgeNearer: edgeNear };
        });
        const left = (document.querySelector('#mol2d .mv-2d-leftout') || {}).textContent || '';
        const m = /(\d+) labels? left out where there is no free space/.exec(left);
        return { drawn: notes.length, leftOut: m ? +m[1] : 0, notes };
      }, md);
    }
    if (k <= NA) {                                       // bond angles: wait for the 3D model, then take its coordinates
      await p.evaluate(() => { const sel = document.getElementById('molLabels'); sel.value = 'angles'; sel.dispatchEvent(new Event('change')); });
      try { await p.waitForFunction(() => { const B = MolView._debug.state.cs.b3d; return B && (B.done || B.error); }, null, { timeout: 90000 }); } catch (e) {}
      await p.evaluate(() => { const sel = document.getElementById('molLabels'); sel.dispatchEvent(new Event('change')); });
      r.angles = await p.evaluate(() => {
        const cs = MolView._debug.state.cs, B = cs.b3d, c = B && B.confs.find(x => x.key === B.cur);
        const notes = document.querySelectorAll('#mol2d .mv-2d-notes text').length;
        const left = (document.querySelector('#mol2d .mv-2d-leftout') || {}).textContent || '';
        return c ? { coords: Array.from(c.coords), optimised: !!c.optimised, energy: c.energy, drawn: notes, foot: left, n: cs.record.n } : { error: B && B.error ? String(B.error) : 'no conformer' };
      });
    }
    res[smi] = r;
    process.stdout.write(`\r${k}/${smis.length}`);
  }
  console.log();
  fs.mkdirSync(path.dirname(out), { recursive: true });
  fs.writeFileSync(out, JSON.stringify({ base, made: new Date().toISOString(), molecules: res }));
  console.log(`wrote ${out}` + (errs.length ? `; page errors: ${errs.slice(0, 3).join(' | ')}` : ''));
  await b.close();
})();
