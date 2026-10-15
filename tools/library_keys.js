// Step 2 of tools/library_add.py: the page's library key (Chem.libKeyHash(Chem.canonicalKey(...)), stereo-blind) for each
// candidate, and for unnamed ones (building blocks) the app's IUPAC name with whether it reads back as the same molecule.
//   python tools/serve.py   (another terminal)
// The name is kept even when the app's own parser cannot read it back (checked 'no'): tools/library_add.py then checks it
// with OPSIN. The namer takes about half a second a molecule, so --pages pages work in parallel.
//   node tools/library_keys.js [http://localhost:8000/app.html] [--in tools/data/library_candidates.json] [--out tools/data/library_keys.json] [--pages 3]
const fs = require('fs'), path = require('path');
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const args = process.argv.slice(2), opt = (k, d) => { const i = args.indexOf(k); return i >= 0 ? args[i + 1] : d; };
const base = args.find(a => /^https?:/.test(a)) || 'http://localhost:8000/app.html';
const inp = opt('--in', path.join(__dirname, 'data', 'library_candidates.json')), out = opt('--out', path.join(__dirname, 'data', 'library_keys.json'));

(async () => {
  const cand = JSON.parse(fs.readFileSync(inp, 'utf8'));
  const b = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const ctx = await b.newContext(), pages = [];
  for (let i = 0; i < +opt('--pages', 3); i++) {
    const p = await ctx.newPage();
    await p.goto(base + '#mol', { waitUntil: 'load' });
    await p.waitForFunction(() => window.Chem && window.MolInfo && window.RDKitLoad, null, { timeout: 60000 });
    await p.evaluate(async () => { window.__R = await window.RDKitLoad(); });
    pages.push(p);
  }
  const parts = [], t0 = Date.now(), N = 50;
  for (let k = 0; k < cand.length; k += N) parts.push(k);
  const done = new Array(parts.length); let next = 0, finished = 0;
  const work = async p => { for (;;) {
    const j = next++; if (j >= parts.length) return;
    const k = parts[j], part = cand.slice(k, k + N).map(c => ({ inchikey: c.inchikey, smiles: c.smiles, name: !c.name }));
    done[j] = await p.evaluate(list => list.map(c => {
      const r = { inchikey: c.inchikey };
      try { r.key = Chem.libKeyHash(Chem.canonicalKey(Chem.parseSmiles(c.smiles))); } catch (e) { r.error = String(e.message || e); return r; }
      if (c.name) {
        let M = null;
        try { M = __R.get_mol(c.smiles); const n = MolInfo.names(__R, M); r.iupac = n.iupac || null; r.checked = n.iupac ? (n.iupacChecked ? 'yes' : 'no') : null; }
        catch (e) { r.iupac = null; } finally { if (M) M.delete(); }
      }
      return r;
    }), part);
    finished += part.length;
    process.stdout.write(`\r${finished}/${cand.length}, ${((Date.now() - t0) / 1000).toFixed(0)} s`);
  } };
  await Promise.all(pages.map(work));
  const rows = [].concat(...done);
  console.log();
  fs.writeFileSync(out, JSON.stringify({ base, made: new Date().toISOString(), rows }));
  console.log(`wrote ${out}: ${rows.filter(r => r.key).length} keys, ${rows.filter(r => r.checked === 'yes').length} checked names`);
  await b.close();
})();
