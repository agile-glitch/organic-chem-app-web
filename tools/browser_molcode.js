// Molecule code checks, in a headless browser: export a molecule with all its labels as a code and as an SDF file, paste
// the code back, and see the same molecule with the same labels; a code whose labels differ from the app's own says so.
//   python tools/serve.py   (another terminal)        node tools/browser_molcode.js [http://localhost:8000/app.html]
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const base = process.argv[2] || 'http://localhost:8000/app.html'; let fails = 0;
const check = (n, ok, d) => { if (!ok) { fails++; console.log('FAIL', n, d === undefined ? '' : JSON.stringify(d).slice(0, 400)); } else console.log('ok  ', n); };
const MODES = ['cip', 'hybridization', 'formal', 'partial', 'oxidation', 'lonepairs', 'pka', 'steric', 'electronic', 'taft'];
(async () => {
  const b = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const ctx = await b.newContext({ viewport: { width: 1400, height: 900 }, acceptDownloads: true, permissions: ['clipboard-read', 'clipboard-write'] });
  const p = await ctx.newPage(); const errs = []; p.on('pageerror', e => errs.push(e.message)); p.on('console', m => { if (m.type() === 'error' && !/favicon/.test(m.location().url || '')) errs.push(m.text()); });
  await p.goto(base + '#mol', { waitUntil: 'load' });
  await p.waitForFunction(() => window.MolView && window.MolCode, null, { timeout: 30000 });
  const show = async text => { await p.locator('#molInput').fill(text); await p.locator('#molShow').click(); await p.waitForFunction(() => { const S = MolView._debug.state; return S.cs && S.cs.record && !document.querySelector('#molStatus:not([hidden]) .mv-spin'); }, null, { timeout: 60000 }); };
  // the labels drawn in the 2D view for each mode: {mode: [[atom, text], ...]}
  const drawn = () => p.evaluate(async modes => {
    const out = {}, sel = document.getElementById('molLabels');
    for (const m of modes.concat('angles')) {
      sel.value = m; sel.dispatchEvent(new Event('change'));
      out[m] = [...document.querySelectorAll('#mol2d .mv-2d-notes text')].map(t => [t.getAttribute('data-atom'), t.textContent]);
    }
    return out;
  }, MODES);
  const exportAs = async fmt => {
    await p.selectOption('#molExportFmt', fmt);
    if (fmt === 'code') { await p.locator('#molExport').click(); await p.waitForSelector('#molResolved .mv-code', { timeout: 90000 }); return p.inputValue('#molResolved .mv-code'); }
    const [dl] = await Promise.all([p.waitForEvent('download', { timeout: 90000 }), p.locator('#molExport').click()]);
    return require('fs').readFileSync(await dl.path(), 'utf8');
  };

  await show('aspirin');
  check('aspirin shown', await p.evaluate(() => MolView._debug.state.cs.record.n) === 13);
  const sdf = await exportAs('sdfcode');           // builds the 3D model the bond angles are measured on
  const before = await drawn();
  check('bond angles drawn before the code', before.angles.length > 10, before.angles);
  check('SDF has the MOL block and the format item', /M {2}END\n>  <ACM\.format>\nAgileCHEM molecule code 1/.test(sdf), sdf.slice(0, 300));
  check('SDF has every label mode', MODES.every(m => sdf.includes(`>  <ACM.labels.${m}>`)));
  check('SDF has the bond angles', /<ACM\.labels\.angles>\n\d+ \d+ \d+ [\d.]+°/.test(sdf));
  check('SDF atoms: 13 heavy, 8 H', sdf.includes('>  <ACM.atoms>\n13 8\n'));
  const code = await exportAs('code');
  check('code starts ACM1:', /^ACM1:[A-Za-z0-9_-]+$/.test(code), code.slice(0, 80));
  console.log('     aspirin: SDF', sdf.length, 'characters, code', code.length);

  // a different molecule, then the code back
  await show('ethanol');
  check('ethanol shown in between', await p.evaluate(() => MolView._debug.state.cs.record.n) === 3);
  await show(code);
  check('the code shows aspirin again', await p.evaluate(() => MolView._debug.state.cs.record.n) === 13);
  check('read from a molecule code', /molecule code/.test(await p.textContent('#molResolved')), await p.textContent('#molResolved'));
  const after = await drawn();
  for (const m of MODES.concat('angles')) check(`same ${m} labels after the code`, JSON.stringify(after[m]) === JSON.stringify(before[m]), { before: before[m], after: after[m] });
  check('banner: the labels are the code\'s and agree', /molecule code's own.*same as this version/.test(await p.textContent('#molBanners')), await p.textContent('#molBanners'));

  // pasting the code into the box shows it at once
  await show('ethanol');
  await p.evaluate(c => navigator.clipboard.writeText(c), code);
  await p.locator('#molInput').fill(''); await p.locator('#molInput').focus(); await p.keyboard.press('Control+V');
  await p.waitForFunction(() => MolView._debug.state.cs.record.n === 13, null, { timeout: 60000 });
  check('pasting the code shows the molecule', true);

  // the SDF file opened again keeps its labels too
  await show('ethanol');
  await p.setInputFiles('#molFile', { name: 'aspirin_labels.sdf', mimeType: 'chemical/x-mdl-sdfile', buffer: Buffer.from(sdf) });
  await p.waitForFunction(() => MolView._debug.state.cs.record.n === 13, null, { timeout: 60000 });
  check('the SDF file opens with the code\'s labels', !!(await p.evaluate(() => MolView._debug.state.cs.code)));

  // a code whose labels differ from this version's: the drawing shows the code's, the banner lists the difference
  check('the pKa labels hold "1 pKa 17"', /<ACM\.labels\.pka>\n(?:.+\n)*1 pKa 17\n/.test(sdf));
  const changed = sdf.replace(/(<ACM\.labels\.pka>\n(?:.+\n)*)1 pKa 17\n/, '$11 pKa 99\n');
  const code2 = await p.evaluate(s => MolCode.encode(s), changed);
  await show(code2);
  const pk = await p.evaluate(() => { const sel = document.getElementById('molLabels'); sel.value = 'pka'; sel.dispatchEvent(new Event('change'));
    return [...document.querySelectorAll('#mol2d .mv-2d-notes text')].map(t => t.textContent); });
  check('the code\'s changed label is drawn', pk.includes('pKa 99'), pk);
  check('banner lists the difference', /1 of them differ.*pKa on C1: code “pKa 99”, this version “pKa 17”/.test(await p.textContent('#molBanners')), await p.textContent('#molBanners'));

  // a damaged code is refused with a reason, and the molecule shown stays
  await p.locator('#molInput').fill(code.slice(0, 30)); await p.locator('#molShow').click();
  await p.waitForFunction(() => /damaged|cut short|not a molecule/.test(document.getElementById('molResolved').textContent), null, { timeout: 20000 });
  check('a damaged code is refused', true);

  // salts: a code is written for the component shown
  await show('CC(=O)[O-].[Na+]');
  const saltCode = await exportAs('code');
  await show(saltCode);
  check('a salt\'s code shows its component', await p.evaluate(() => MolView._debug.state.cs.record.n) === 4);

  check('no page errors', errs.length === 0, errs);
  await b.close();
  console.log(fails ? `${fails} FAILED` : 'all passed');
  process.exit(fails ? 1 : 0);
})();
