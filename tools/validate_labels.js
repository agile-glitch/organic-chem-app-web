// Label validation, step 1 of 2: the values every Labels mode of the Molecule tab would show, for a set of molecules, as
// the page itself works them out (MolInfo.analyse, MoleculePka, MoleculeEffects). tools/validate_labels.py then checks
// them against independent references (RDKit Python, the IUPAC measured pKa values, a hand-checked textbook list).
//   python tools/serve.py   (another terminal)
//   node tools/validate_labels.js [http://localhost:8000/app.html] [--set pka|textbook|all] [--limit N] [--out tools/data/label_dump.json]
// The molecule set: 'pka' = every molecule of data/pka.js (IUPAC Digitized pKa Dataset, 8,894 molecules: acids, bases,
// 539 with stereocentres); 'textbook' = the SMILES of tools/label_textbook.json; 'all' = both.
const fs = require('fs'), path = require('path');
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const args = process.argv.slice(2), opt = (k, d) => { const i = args.indexOf(k); return i >= 0 ? args[i + 1] : d; };
const base = args.find(a => /^https?:/.test(a)) || 'http://localhost:8000/app.html';
const set = opt('--set', 'all'), limit = +opt('--limit', 0), out = opt('--out', path.join(__dirname, 'data', 'label_dump.json'));

function molecules() {
  const smis = [];
  if (set === 'textbook' || set === 'all') {
    const T = JSON.parse(fs.readFileSync(path.join(__dirname, 'label_textbook.json'), 'utf8'));
    for (const c of T.cases) smis.push(c.smiles);
    for (const sr of T.steric_series.series) smis.push(...sr.smiles);
  }
  if (set === 'pka' || set === 'all') {
    const s = fs.readFileSync(path.join(__dirname, '..', 'data', 'pka.js'), 'utf8');
    smis.push(...Object.keys(JSON.parse(s.slice(s.indexOf('{'), s.lastIndexOf('}') + 1)).molecules));
  }
  const u = [...new Set(smis)];
  return limit ? u.slice(0, limit) : u;
}

(async () => {
  const smis = molecules();
  const b = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const p = await b.newPage();
  const errs = []; p.on('pageerror', e => errs.push(e.message));
  await p.goto(base + '#mol', { waitUntil: 'load' });
  await p.waitForFunction(() => window.MolInfo && window.MoleculePka && window.MoleculeEffects && window.RDKitLoad, null, { timeout: 60000 });
  await p.evaluate(async () => { window.__R = await RDKitLoad(); });
  const res = {};
  const t0 = Date.now();
  for (let k = 0; k < smis.length; k += 200) {
    const part = await p.evaluate(batch => {
      const R = window.__R, MI = window.MolInfo, out = {};
      const MEAS = window.PKA_MEASURED;
      const sites = r => r ? ['acid', 'base'].reduce((o, kind) => (o[kind] = r[kind].map(x => ({ atoms: (x.atoms || [x.atom]).slice(), pKa: x.pKa, fragments: [...new Set(x.fragments)], source: x.source })), o), {}) : null;
      for (const smi of batch) {
        let m = null;
        try {
          m = R.get_mol(smi);
          if (!m || !m.is_valid()) { out[smi] = { error: 'RDKit.js cannot read it' }; continue; }
          const rec = MI.analyse(R, m), n = rec.n, A = rec.atoms.slice(0, n);
          const pka = sites(window.MoleculePka.forMol(R, m));
          window.PKA_MEASURED = { molecules: {} };                       // the rules alone: no measured value to fall back on
          let rules; try { rules = sites(window.MoleculePka.forMol(R, m)); } finally { window.PKA_MEASURED = MEAS; }
          const eff = (window.MoleculeEffects.forMol(R, m, rec) || []).map(e => ({ fragment: e.fragment, atom: e.atom, steric: e.steric.cls, score: e.steric.score,
            charge: e.electronic.charge, sigma: e.electronic.sigma, sigmaStar: e.taft.sigmaStar, EsMin: e.taft.EsMin,
            groups: e.taft.groups.map(g => [g.name, g.sigmaStar, g.Es]) }));
          out[smi] = {
            canonical: m.get_smiles(), n,
            labels: Object.fromEntries(['cip', 'hybridization', 'formal', 'partial', 'oxidation', 'lonepairs'].map(md => [md, MI.labelsFor(rec, md).slice(0, n)])),
            el: A.map(a => a.el), cip: A.map(a => a.cip), hyb: A.map(a => a.hybridization.label), rdHyb: A.map(a => a.rdkitHybridization),
            charge: A.map(a => a.charge), q: A.map(a => a.gasteiger), ox: A.map(a => a.oxidationState.value), oxFlag: A.map(a => a.oxidationState.flag || null),
            lp: A.map(a => a.lonePairs), ez: rec.bonds.filter(bd => bd.ez).map(bd => [bd.a, bd.b, bd.ez]),
            pka, rules, effects: eff,
          };
        } catch (e) { out[smi] = { error: String(e && e.message || e) }; }
        finally { if (m) m.delete(); }
      }
      return out;
    }, smis.slice(k, k + 200));
    Object.assign(res, part);
    process.stdout.write(`\r${Math.min(k + 200, smis.length)}/${smis.length} molecules, ${((Date.now() - t0) / 1000).toFixed(0)} s`);
  }
  console.log();
  fs.mkdirSync(path.dirname(out), { recursive: true });
  fs.writeFileSync(out, JSON.stringify({ base, set, made: new Date().toISOString(), molecules: res }));
  const bad = Object.values(res).filter(r => r.error).length;
  console.log(`wrote ${out}: ${smis.length} molecules, ${bad} not analysed` + (errs.length ? `; page errors: ${errs.slice(0, 3).join(' | ')}` : ''));
  await b.close();
})();
