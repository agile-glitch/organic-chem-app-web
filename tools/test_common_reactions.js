/* Test the Common reactions panel of the Reactions page in a real headless browser.

     python tools/serve.py 8765 &            # serve the app
     NODE_PATH=<dir with playwright> node tools/test_common_reactions.js [http://localhost:8765/] [--all]

   Needs Playwright with Chromium (npm i playwright). It opens the Reactions tab, opens the list, checks the entries,
   the search and the keyboard use, loads reactions through the list and compares the page's products with the products
   recorded in each rule, and (with --all) loads every one of the rules. Exit code 1 if anything fails. */
const url = process.argv.find(a => a.startsWith('http')) || 'http://localhost:8765/';
const ALL = process.argv.includes('--all');
let chromium;
try { ({ chromium } = require('playwright')); } catch (e) { console.log('SKIP: Playwright is not installed (npm i playwright; set NODE_PATH)'); process.exit(0); }
const fails = [];
const check = (ok, what) => { console.log((ok ? 'PASS  ' : 'FAIL  ') + what); if (!ok) fails.push(what); };

(async () => {
  const browser = await chromium.launch();
  const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });
  const errors = [];
  page.on('pageerror', e => errors.push('pageerror: ' + e.message));
  page.on('console', m => { if (m.type() === 'error') errors.push('console: ' + m.text()); });
  await page.goto(url);
  await page.click('button.tab[data-page="reactions"]');
  await page.waitForSelector('#rxn-common-btn', { state: 'visible' });
  const rules = await page.evaluate(() => window.REACTION_RULES.rules.map(r => ({ id: r.id, name: r.name, code: r.code })));
  check(rules.length === 102, `the page has ${rules.length} workbook reactions`);

  // layout: the button sits left of the gear and the code box does not run under it
  const box = async sel => page.evaluate(s => { const r = document.querySelector(s).getBoundingClientRect(); return { l: r.left, r: r.right, t: r.top, b: r.bottom }; }, sel);
  const [bc, bb, bg] = [await box('#rxn-code'), await box('#rxn-common-btn'), await box('#rxn-gear')];
  check(bc.r <= bb.l + 0.5 && bb.r <= bg.l + 0.5 && bb.b - bb.t > 20, 'the code box, the Common reactions button and the gear sit side by side without overlap');

  // open the list
  await page.click('#rxn-common-btn');
  await page.waitForSelector('#rxn-common', { state: 'visible' });
  const n = await page.locator('.rxn-c-item').count();
  const heads = await page.locator('.rxn-c-h > span:first-child').allTextContents();
  check(n === rules.length, `the list shows all ${n} reactions`);
  console.log('      families: ' + (await page.evaluate(() => { const c = {}; for (const x of window.ReactionsCommon.list()) c[x.family] = (c[x.family] || 0) + 1; return Object.entries(c).map(([k, v]) => `${k} ${v}`).join('; '); })));
  const other = await page.evaluate(() => window.ReactionsCommon.list().filter(x => x.family === 'Other').map(x => x.id + ' ' + x.name));
  check(heads.length >= 6 && other.length <= 2, `reactions are grouped under ${heads.length} headings; ${other.length} left in "Other"` + (other.length ? ' (' + other.join('; ') + ': it consumes no fragment, so no family can be derived)' : ''));
  check((await page.textContent('#rxn-c-sum')).startsWith('102'), 'the count reads 102 reactions');
  check(await page.evaluate(() => document.activeElement.id) === 'rxn-c-q', 'the search box has focus on opening');
  const sub = await page.textContent('.rxn-c-item[data-id="R001"] .d');
  check(/H2SO4/i.test(sub) && /ethanol/.test(sub), `an entry shows its reagent and conditions ("${sub}")`);
  await page.screenshot({ path: process.env.SHOTS ? process.env.SHOTS + '/common_open.png' : '/tmp/common_open.png' });

  // search
  await page.fill('#rxn-c-q', 'grignard');
  const g = await page.locator('.rxn-c-item').count();
  check(g >= 3 && g < 20, `searching "grignard" leaves ${g} reactions`);
  await page.fill('#rxn-c-q', 'zzzzq');
  check(await page.locator('.rxn-c-item').count() === 0 && (await page.textContent('#rxn-c-list')).includes('No reaction matches'), 'a search with no match says so');
  await page.fill('#rxn-c-q', 'nabh4 ketone');
  check(await page.locator('.rxn-c-item').count() >= 1, 'two words must both match (reagent + fragment)');
  await page.fill('#rxn-c-q', '');
  check(await page.locator('.rxn-c-item').count() === rules.length, 'clearing the search brings every reaction back');

  // keyboard: Escape closes, ArrowDown moves into the list
  await page.keyboard.press('Escape');
  check(await page.isHidden('#rxn-common'), 'Escape closes the list');
  await page.click('#rxn-common-btn');
  await page.keyboard.press('ArrowDown');
  check(await page.evaluate(() => document.activeElement.classList.contains('rxn-c-item')), 'ArrowDown moves from the search box into the list');
  await page.click('#rxn-c-close');
  check(await page.isHidden('#rxn-common'), 'the × closes the list');

  // load reactions through the list and compare with the recorded products
  const canon = async list => page.evaluate(async smis => { const R = await (window.rdkitReady ? window.rdkitReady() : window.RDKit); return smis; }, list).catch(() => list);
  const productsOf = async () => page.evaluate(() => document.getElementById('rxn-prod').value.split('.').filter(Boolean).sort().join('.'));
  const runViaList = async id => {
    await page.evaluate(() => { document.getElementById('rxn-prod').value = ''; });
    await page.evaluate(i => window.ReactionsCommon.load(i), id);
    const t0 = Date.now();
    while (Date.now() - t0 < 60000) { if ((await productsOf())) break; await page.waitForTimeout(100); }
    await page.waitForTimeout(300);
    return { code: await page.inputValue('#rxn-code'), prods: await productsOf() };
  };
  const sample = ALL ? rules.map(r => r.id) : ['R001', 'R003', 'R017', 'R026', 'R033', 'R037', 'R040', 'R051'];
  let same = 0; const bad = [];
  for (const id of sample) {
    const r = rules.find(x => x.id === id), p = r.code.split('~');
    const out = await runViaList(id);
    const exp = p[3].split('.').filter(Boolean).sort().join('.');
    const codeOk = out.code === p.slice(0, 3).join('~') + '~' || out.code.startsWith(p[0].split('.')[0].slice(0, 3));
    if (out.prods === exp && codeOk) same++; else bad.push(`${id} ${r.name}: page gave "${out.prods}", recorded "${exp}"`);
  }
  check(!bad.length, `${same} of ${sample.length} reactions loaded from the list give the recorded products` + (bad.length ? '\n        ' + bad.slice(0, 6).join('\n        ') : ''));
  await page.screenshot({ path: process.env.SHOTS ? process.env.SHOTS + '/common_loaded.png' : '/tmp/common_loaded.png' });

  // the page's own regression tests still pass with the new code in place
  if (ALL) {
    await page.evaluate(() => window.ReactionsTests.run(false));
    const t0 = Date.now(); while (await page.evaluate(() => window.ReactionsTests.running()) && Date.now() - t0 < 1800000) await page.waitForTimeout(1000);
    const res = await page.evaluate(() => window.ReactionsTests.results().map(r => ({ ok: r.ok, n: r.t.name })));
    check(res.length > 100 && res.every(r => r.ok), `the page's own tests: ${res.filter(r => r.ok).length} of ${res.length} pass`);
  }
  check(errors.length === 0, 'no console or page errors' + (errors.length ? ': ' + errors.slice(0, 3).join(' | ') : ''));
  await browser.close();
  console.log('\n' + (fails.length ? fails.length + ' FAILED' : 'ALL TESTS PASSED'));
  process.exit(fails.length ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
