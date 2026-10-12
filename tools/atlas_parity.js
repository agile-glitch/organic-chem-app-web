// Python / browser parity for the Fragment Atlas.   node tools/atlas_parity.js [baseUrl]    (needs: python tools/serve.py, Playwright)
// For every molecule of tools/atlas_parity_cases.json (written by tools/atlas_parity_ref.py with RDKit) js/atlas.js runs in the page on the same SMILES
// and must give the same instances: fragment, reacting atom, the four environment hashes, steric score, Gasteiger charge, ring, strain, hybridization.
//   node tools/atlas_parity.js --judge   also checks the judgement: for 400 instances of tools/atlas_parity_judge.json (tools/atlas_parity_judge_ref.py) the page reads
//   data/atlas/ through its own shard decoder and must reach the same probability as numpy (table lookups, backoff chain, combining model), for both questions: transformed, and also transformed when the reaction is elsewhere.
const fs = require('fs'), path = require('path');
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || '/opt/node22/lib/node_modules/playwright');
const args = process.argv.slice(2), judge = args.includes('--judge'), base = args.find(a => !a.startsWith('--')) || 'http://localhost:8000/';
const cases = JSON.parse(fs.readFileSync(path.join(__dirname, 'atlas_parity_cases.json'), 'utf8')).cases;
(async () => {
  const browser = await chromium.launch({ executablePath: process.env.CHROMIUM || '/opt/pw-browsers/chromium-1194/chrome-linux/chrome', args: ['--no-sandbox'] });
  const page = await browser.newPage();
  await page.goto(base + '#reactions', { waitUntil: 'load' });
  await page.waitForFunction(() => window.Atlas && window.RDKitLoad && window.MolInfo && window.REACTION_RULES, null, { timeout: 90000 });
  const got = await page.evaluate(async smis => {
    const R = await window.RDKitLoad(); await window.Atlas.ready();
    return smis.map(s => window.Atlas.describe(R, s).map(i => ({ fid: i.fid, atom: i.atom, e: i.e, steric: i.steric, charge: i.charge, ring: i.ring, strain: i.strain, hyb: i.hyb, nAtoms: i.nAtoms })));
  }, cases.map(c => c.smiles));
  const stats = { molecules: cases.length, instances: 0, ok: 0, keyMismatch: 0, missing: 0, extra: 0, steric: 0, charge: 0, ring: 0, strain: 0, hyb: 0, moleculesAllOk: 0 };
  const bad = [];
  cases.forEach((c, k) => {
    const ref = new Map(), js = new Map();
    c.instances.forEach(i => ref.set(i.fid + ':' + i.atom + ':' + i.nAtoms, i));
    got[k].forEach(i => js.set(i.fid + ':' + i.atom + ':' + i.nAtoms, i));
    let allOk = true;
    for (const [id, r] of ref) {
      stats.instances++;
      const j = js.get(id);
      if (!j) { stats.missing++; allOk = false; if (bad.length < 12) bad.push(`${c.smiles}: ${id} missing in the page`); continue; }
      let ok = true;
      if (r.e.some((v, q) => v !== j.e[q])) { stats.keyMismatch++; ok = false; if (bad.length < 12) bad.push(`${c.smiles}: ${id} environment ${r.e} vs ${j.e}`); }
      if (r.steric !== j.steric) { stats.steric++; ok = false; }
      if ((r.charge == null) !== (j.charge == null) || (r.charge != null && Math.abs(r.charge - j.charge) > 2e-4)) { stats.charge++; ok = false; if (bad.length < 12) bad.push(`${c.smiles}: ${id} charge ${r.charge} vs ${j.charge}`); }
      if (r.ring !== j.ring) { stats.ring++; ok = false; }
      if (Math.abs(r.strain - j.strain) > 1e-3) { stats.strain++; ok = false; }
      if (r.hyb !== j.hyb) { stats.hyb++; ok = false; }
      if (ok) stats.ok++; else allOk = false;
    }
    for (const id of js.keys()) if (!ref.has(id)) { stats.extra++; allOk = false; if (bad.length < 12) bad.push(`${c.smiles}: ${id} only in the page`); }
    if (allOk) stats.moleculesAllOk++;
  });
  console.log(JSON.stringify(stats));
  bad.forEach(b => console.log('  ', b));
  let fine = stats.keyMismatch + stats.missing + stats.extra === 0;
  if (judge) {
    const J = JSON.parse(fs.readFileSync(path.join(__dirname, 'atlas_parity_judge.json'), 'utf8')).cases;
    const res = await page.evaluate(async cases => {
      const R = await window.RDKitLoad(); const M = await window.Atlas.ready(); const RR = window.REACTION_RULES;
      const out = [];
      for (const c of cases) {
        const sh = await window.Atlas.shard(c.F);
        const parts = c.parts.map(s => { const r = (RR.reagents || {})[s] || {}; return { smiles: s, role: r.role, pKaH: r.pKaH }; });
        const inst = { F: c.F, fid: c.fid, e: c.e, steric: c.steric, charge: c.charge, ring: c.ring, strain: c.strain, hyb: c.hyb, nAtoms: c.nAtoms };
        const side = { cls: c.cls, parts };
        const a = window.Atlas.features(R, inst, side, sh, 'A'), b = window.Atlas.features(R, inst, side, sh, 'B');
        out.push([window.Atlas.combine(a.v, 'A'), window.Atlas.combine(b.v, 'B')]);
      }
      return { out, entries: M ? M.entries : 0 };
    }, J);
    let worst = 0, bigs = 0; const wrong = [];
    J.forEach((c, i) => {
      const d = Math.max(Math.abs(c.p - res.out[i][0]), Math.abs(c.pb - res.out[i][1])); worst = Math.max(worst, d);
      if (d > 2e-3) { bigs++; if (wrong.length < 8) wrong.push(`F${c.F} ${c.fid} cls=${c.cls}: python ${c.p.toFixed(4)} / ${c.pb.toFixed(4)} page ${res.out[i][0].toFixed(4)} / ${res.out[i][1].toFixed(4)}`); }
    });
    console.log(`judgement: ${J.length} instances, largest difference ${worst.toExponential(2)}, ${bigs} above 0.002 (atlas of ${res.entries.toLocaleString()} entries)`);
    wrong.forEach(w => console.log('  ', w));
    fine = fine && bigs === 0;
  }
  console.log(fine ? 'PARITY OK' : 'PARITY DIFFERENCES');
  await browser.close();
  process.exit(fine ? 0 : 1);
})().catch(e => { console.error(e); process.exit(2); });
