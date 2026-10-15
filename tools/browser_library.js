// The compound library in the page (js/chem.js + data/library_more.js): every entry of data/library_more.js is found by
// name, its key is the page's own key of its SMILES, and the page writes back the same molecule, stereo included (RDKit
// canonical SMILES compared); a few named look-ups (drugs, controlled substances, a ChEMBL ID, the corrected stereo).
// Where the app's own reader loses a stereo mark (ring cis/trans in lindane, tranexamic acid ...), the Molecule tab must
// still show the entry's exact stereoisomer: it reads the entry's SMILES with RDKit. Those names are shown in the tab
// and its RDKit SMILES compared with the entry's.
//   python tools/serve.py &    node tools/browser_library.js [http://localhost:8000/app.html]
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const base = process.argv[2] || 'http://localhost:8000/app.html';
(async () => {
  const b = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const p = await b.newPage(); const errs = []; p.on('pageerror', e => errs.push(e.message));
  await p.goto(base + '#mol', { waitUntil: 'load' });
  await p.waitForFunction(() => window.Chem && window.Chem.MOL_LIBRARY && window.RDKitLoad, null, { timeout: 60000 });
  const r = await p.evaluate(async () => {
    const R = await RDKitLoad(), C = window.Chem, out = { fail: [], n: 0 };
    const canon = s => { let m = null; try { m = R.get_mol(s); return m ? m.get_smiles() : null; } finally { if (m) m.delete(); } };
    const more = (window.MOL_LIBRARY_MORE || '').split('\n').filter(Boolean).map(l => l.split('|'));
    out.more = more.length; out.total = C.MOL_LIBRARY.length;
    for (const [name, smi, key] of more) {
      out.n++;
      let g; try { g = C.parseSmiles(smi); } catch (e) { out.fail.push(`${name}: the page cannot read ${smi}`); continue; }
      if (C.libKeyHash(C.canonicalKey(g)) !== key) { out.fail.push(`${name}: key`); continue; }
      const e = C.MOL_BY_NAME.get(name.toLowerCase());
      if (!e || e.smiles !== smi) { out.fail.push(`${name}: not found by its name`); continue; }
      if (canon(C.toSmiles(g)) !== canon(smi)) (out.lossy = out.lossy || []).push([name, canon(smi)]);
    }
    out.tab = [];
    for (const [name, want] of (out.lossy || [])) {
      const r = await MolView.agilesFor(name);
      const got = r.ok && r.ids && r.ids.smiles_rdkit;
      // a salt or hydrate: the tab shows its main component, which must be one of the entry's components exactly
      if (got !== want && !(got && want.split('.').includes(got))) out.fail.push(`${name}: the Molecule tab shows ${got || r.error}, not ${want}`); else out.tab.push(name);
    }
    const look = {};
    for (const q of ['lsd', 'psilocybin', 'psilocin', 'mdma', 'fentanyl citrate', 'dronabinol', 'sodium oxybate', 'chembl2', 'methamphetamine', 'imatinib']) {
      try { const x = C.searchMolecule(q); look[q] = canon(C.toSmiles(x.graph)); } catch (e) { look[q] = 'NOT FOUND: ' + e.message.slice(0, 60); }
    }
    out.look = look;
    return out;
  });
  const want = { lsd: 'CCN(CC)C(=O)[C@@H]1C=C2c3cccc4[nH]cc(c34)C[C@H]2N(C)C1', methamphetamine: 'CN[C@@H](C)Cc1ccccc1' };
  let bad = r.fail.length;
  console.log(`${r.total} library compounds, ${r.more} from data/library_more.js; ${r.n - r.fail.length}/${r.n} pass`);
  r.fail.slice(0, 20).forEach(f => console.log('  FAIL', f));
  console.log(`  ${(r.lossy || []).length} entries lose a stereo mark in the app's own reader; the Molecule tab shows ${r.tab.length} of them exactly (RDKit reads the entry's SMILES)`);
  for (const [q, s] of Object.entries(r.look)) {
    const ok = want[q] ? s === want[q] : !/^NOT FOUND/.test(s);
    if (!ok) bad++;
    console.log(`  ${ok ? 'ok  ' : 'FAIL'} ${q}: ${s}`);
  }
  if (errs.length) { console.log('page errors:', errs.slice(0, 3)); bad++; }
  await b.close(); process.exit(bad ? 1 : 0);
})();
