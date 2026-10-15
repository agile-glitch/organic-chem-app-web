// The AGILES table: every compound of the app's library with its AGILES code and identifiers (IUPAC name, SMILES, InChI,
// InChIKey, CAS), as the Molecule tab itself works them out (MolView._debug: no 3D model, so the codes hold no bond angles;
// the tab measures them when a code is shown). tools/agiles_table.py then writes the spreadsheet and checks every code.
//   python tools/serve.py   (another terminal)
//   node tools/agiles_table.js [http://localhost:8000/app.html] [--names "aspirin,caffeine"] [--limit N] [--out tools/data/agiles_table.json]
const fs = require('fs'), path = require('path');
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const args = process.argv.slice(2), opt = (k, d) => { const i = args.indexOf(k); return i >= 0 ? args[i + 1] : d; };
const base = args.find(a => /^https?:/.test(a)) || 'http://localhost:8000/app.html';
const limit = +opt('--limit', 0), out = opt('--out', path.join(__dirname, 'data', 'agiles_table.json')), names = opt('--names', '');

(async () => {
  const b = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const p = await (await b.newContext({ viewport: { width: 1200, height: 800 } })).newPage();
  const errs = []; p.on('pageerror', e => errs.push(e.message));
  await p.goto(base + '#mol', { waitUntil: 'load' });
  await p.waitForFunction(() => window.MolView && window.MolView.agilesFor && window.Chem && window.Chem.MOL_LIBRARY, null, { timeout: 60000 });
  // the library entries (name, SMILES, CAS), or the names asked for
  let list = names ? names.split(',').map(n => ({ name: n.trim(), input: n.trim() }))
    : await p.evaluate(() => window.Chem.MOL_LIBRARY.map(e => ({ name: e.name, input: e.smiles, smiles: e.smiles, libCas: e.cas || null })));
  if (limit) list = list.slice(0, limit);
  const rows = [], t0 = Date.now();
  for (let k = 0; k < list.length; k++) {
    const e = list[k];
    const r = await p.evaluate(async x => { try { return await MolView.agilesFor(x); } catch (err) { return { ok: false, error: String(err && err.message || err) }; } }, e.input);
    rows.push(Object.assign({ name: e.name, input: e.input, libCas: e.libCas || null }, r));
    if (k % 50 === 0 || k === list.length - 1) process.stdout.write(`\r${k + 1}/${list.length} molecules, ${((Date.now() - t0) / 1000).toFixed(0)} s`);
  }
  console.log();
  fs.mkdirSync(path.dirname(out), { recursive: true });
  fs.writeFileSync(out, JSON.stringify({ base, made: new Date().toISOString(), rows }));
  console.log(`wrote ${out}: ${rows.filter(r => r.ok).length} of ${rows.length} with a code` + (errs.length ? `; page errors: ${errs.slice(0, 3).join(' | ')}` : ''));
  await b.close();
})();
