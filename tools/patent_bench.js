// The Reactions page against recorded reactions of US patents it never learned from.
//
//   python tools/serve.py &
//   node tools/patent_bench.js sample [n=1000] [seed=1]        # tools/patent_bench_cases.json: n reactions, one per held-out US patent
//   node tools/patent_bench.js run [http://localhost:8000/app.html] [cases.json] [results.json]
//   python tools/patent_bench_grade.py                          # tools/patent_bench_report.json / .md
//
// Held out: the 1 patent in 20 with md5(patent id) mod 20 == 0 (tools/rel_common.py, TEST). The relevance model never trained on them and
// nothing was tuned on them. The page's learned data (data/learned/rx) still holds their reactions for the exact look-up, so the run turns
// that look-up off: otherwise the page would read the answer instead of predicting it.
// Each case is the recorded code without its product (reactants>reagents> | cond=recorded solvents), typed into the code bar, then the arrow,
// exactly as a user would; the page's products and its reading are kept as they are and graded afterwards with RDKit.
// Needs Playwright (npm i playwright). PLAYWRIGHT_MODULE / CHROMIUM override where they are found.
const fs = require('fs'), path = require('path'), zlib = require('zlib'), crypto = require('crypto');
const RxnCode = require(path.join(__dirname, '..', 'js', 'rxn_code.js'));
const HERE = __dirname, APP = path.join(HERE, '..');
const LIMIT_S = 120;                                     // a reaction the page has not answered in this long counts as an error
const CASES = path.join(HERE, 'patent_bench_cases.json'), RESULTS = path.join(HERE, 'patent_bench_results.json');

const patentOf = source => String(source || '').trim().split(/\s+/)[0] || '?';
const heldOut = patent => BigInt('0x' + crypto.createHash('md5').update(patent).digest('hex')) % 20n === 0n;     // rel_common.split_of_patent == TEST
function rng(seed) { let s = seed >>> 0 || 1; return () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; }; }

function sample(n, seed) {
  const dir = path.join(APP, 'data', 'learned', 'rx'), byPatent = new Map();
  let total = 0;
  for (const f of fs.readdirSync(dir).sort()) {
    const s = fs.readFileSync(path.join(dir, f), 'utf8');
    const rows = JSON.parse(zlib.gunzipSync(Buffer.from(s.match(/"gz":"([^"]+)"/)[1], 'base64')));
    for (const r of rows) {
      total++;
      const pat = patentOf(r[7]);
      if (!/^US/.test(pat) || !heldOut(pat)) continue;
      (byPatent.get(pat) || byPatent.set(pat, []).get(pat)).push(r);
    }
  }
  const rand = rng(seed), patents = [...byPatent.keys()].sort();
  for (let i = patents.length - 1; i > 0; i--) { const j = Math.floor(rand() * (i + 1)); [patents[i], patents[j]] = [patents[j], patents[i]]; }
  const cases = patents.slice(0, n).map(pat => {
    const rows = byPatent.get(pat), r = rows[Math.floor(rand() * rows.length)], p = RxnCode.parse(r[1]);
    return { id: r[0], patent: pat, source: r[7], year: r[8], yield: r[9],
             code: RxnCode.build({ reactants: p.reactants.join('.'), agents: p.agents.join('.'), products: '', cond: p.cond }),
             reactants: p.reactants, reagent: p.agents.join('.'), solvents: p.cond, product: p.products.join('.') };
  });
  const heldRx = [...byPatent.values()].reduce((a, b) => a + b.length, 0);
  fs.writeFileSync(CASES, JSON.stringify({ about: `${cases.length} recorded reactions, one from each of ${cases.length} randomly chosen held-out US patents ` +
    `(seed ${seed}; ${byPatent.size.toLocaleString()} held-out US patents, ${heldRx.toLocaleString()} of their reactions among the ${total.toLocaleString()} learned). ` +
    'Made by tools/patent_bench.js sample.', seed, cases }, null, 1));
  console.log(`wrote ${CASES}: ${cases.length} cases from ${byPatent.size} held-out US patents (${heldRx} reactions of ${total})`);
}

async function run(base, casesPath, outPath) {
  const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
  const { cases } = JSON.parse(fs.readFileSync(casesPath, 'utf8'));
  const done = fs.existsSync(outPath) ? JSON.parse(fs.readFileSync(outPath, 'utf8')) : [];   // resumes where a stopped run left off
  const seen = new Set(done.map(x => x.id));
  const browser = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  let page = null;
  const open = async () => {
    if (page) await page.close().catch(() => {});
    page = await browser.newPage({ viewport: { width: 1500, height: 900 } });
    await page.goto(base + '#reactions', { waitUntil: 'load' });
    await page.waitForFunction(() => window.ReactionsTests && window.ReactionsTests.predict, null, { timeout: 120000 });
    await page.evaluate(() => window.ReactionsTests.lookup(false));
  };
  await open();
  const t0 = Date.now();
  let k = 0;
  for (const c of cases) {
    if (seen.has(c.id)) continue;
    const t1 = Date.now();
    let r;
    let timer;                                           // the page's own 90 s limit cannot fire while its main thread is busy: a hard limit here
    const hung = new Promise(res => { timer = setTimeout(() => res({ error: `hung (no answer in ${LIMIT_S} s)`, hung: true }), LIMIT_S * 1000); });
    try { r = await Promise.race([page.evaluate(code => window.ReactionsTests.predict(code), c.code), hung]); }
    catch (e) { r = { error: 'page: ' + String(e.message || e).slice(0, 200), hung: true }; }
    clearTimeout(timer);
    if (r.hung) { delete r.hung; await open(); }        // a crashed or hung page: a fresh one for the next case
    done.push({ id: c.id, ms: Date.now() - t1, ...r });
    if (++k % 10 === 0 || done.length === cases.length) {
      fs.writeFileSync(outPath, JSON.stringify(done));
      const rate = k / ((Date.now() - t0) / 1000);
      console.log(`${done.length} / ${cases.length}  (${rate.toFixed(2)} per s, about ${Math.round((cases.length - done.length) / rate / 60)} min left)`);
    }
  }
  fs.writeFileSync(outPath, JSON.stringify(done));
  await browser.close();
  console.log(`wrote ${outPath}: ${done.length} results`);
}

const [cmd, a, b, c] = process.argv.slice(2);
if (cmd === 'sample') sample(+(a || 1000), +(b || 1));
else if (cmd === 'run') run(a || 'http://localhost:8000/app.html', b || CASES, c || RESULTS).catch(e => { console.error(e); process.exit(1); });
else { console.log('node tools/patent_bench.js sample [n] [seed] | run [url] [cases.json] [results.json]'); process.exit(1); }
