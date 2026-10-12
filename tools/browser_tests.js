// Run the Reactions page's own "Run all tests" in a headless browser and print the result.
//   python tools/serve.py            (another terminal; serves http://localhost:8000)
//   node tools/browser_tests.js [http://localhost:8000/app.html] [results.json]
// Needs Playwright (npm i playwright; npx playwright install chromium). PLAYWRIGHT_MODULE / CHROMIUM override where they are found.
// Exit status 1 when any case fails. The cases: every reaction of reaction_rules.xlsx, tools/reaction_tests.json, tools/learned_tests.json.
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const base = process.argv[2] || 'http://localhost:8000/app.html', out = process.argv[3] || '';
(async () => {
  const browser = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const page = await browser.newPage({ viewport: { width: 1500, height: 900 } });
  await page.goto(base + '#reactions', { waitUntil: 'load' });
  await page.waitForFunction(() => window.ReactionsTests, null, { timeout: 90000 });
  await page.evaluate(() => window.ReactionsTests.run(false));
  let last = '';
  for (;;) {
    await page.waitForTimeout(4000);
    const s = await page.evaluate(() => ({ running: window.ReactionsTests.running(), sum: document.getElementById('rxn-t-sum').textContent }));
    if (s.sum !== last) { console.log(new Date().toISOString().slice(11, 19), s.sum); last = s.sum; }
    if (!s.running && /passed/.test(s.sum)) break;
  }
  const res = await page.evaluate(() => window.ReactionsTests.results().map(x => ({ group: x.t.group, name: x.t.name, code: x.t.code, ok: x.ok, why: x.why })));
  if (out) require('fs').writeFileSync(out, JSON.stringify(res, null, 1));
  const groups = {};
  res.forEach(r => { const g = groups[r.group] = groups[r.group] || [0, 0]; g[1]++; if (r.ok) g[0]++; });
  console.log(JSON.stringify(groups));
  const bad = res.filter(r => !r.ok);
  console.log('TOTAL', res.length - bad.length, '/', res.length);
  bad.slice(0, 40).forEach(r => console.log('FAIL', r.group, '|', r.name, '|', r.why.join(' ; ').slice(0, 200)));
  await browser.close();
  process.exit(bad.length ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
