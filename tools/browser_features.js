// Checks of the newer features, in a headless browser: the Atlas, the reagent suggestions, the learned retrosynthesis, the equation balancer,
// the Reactions page's one-line verdict and its "Open Mechanism in Draw" button, the cross-app links.
//   python tools/serve.py   (another terminal)        node tools/browser_features.js [http://localhost:8000/app.html]
// Each check states a chemical fact the program must get right; the exit status is 1 when one fails.
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const base = process.argv[2] || 'http://localhost:8000/app.html';
let fails = 0, n = 0;
const check = (name, ok, detail) => { n++; if (!ok) { fails++; console.log('FAIL', name, detail === undefined ? '' : '\n    ' + JSON.stringify(detail).slice(0, 400)); } else console.log('ok  ', name); };
(async () => {
  const browser = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const page = await browser.newPage({ viewport: { width: 1500, height: 1000 } });
  const errors = []; page.on('pageerror', e => errors.push(e.message));
  await page.goto(base + '#reactions', { waitUntil: 'load' });
  await page.waitForFunction(() => window.Advisor && window.RxnPage && window.Atlas && window.RetroLearned && window.Stoich, null, { timeout: 90000 });

  // ---- Atlas: the groups and the two questions
  const atlas = await page.evaluate(async () => {
    const R = await window.RDKitLoad(); await window.Atlas.ready();
    const J = await window.Atlas.judge(R, ['CCOC(=O)c1ccc([N+](=O)[O-])cc1'], '[BH4-].[Na+]');
    const by = {}; J.forEach(j => { by[j.fid] = j; });
    const pro = await window.Atlas.judge(R, ['CC(=O)Cl', 'Nc1ccccc1'], '');
    return { meta: !!window.Atlas.meta, entries: window.Atlas.meta.entries, ester: by.F004 && { p: by.F004.p, pb: by.F004.pb, n: by.F004.nReagentB }, nitro: Object.values(by).filter(j => /nitro/i.test(j.name)).map(j => ({ p: j.p, pb: j.pb })),
      acyl: pro.filter(j => /acyl chloride/i.test(j.name)).map(j => j.p) };
  });
  check('Atlas loads its 1.5 million entries', atlas.meta && atlas.entries > 1.4e6, atlas.entries);
  check('an acyl chloride is judged transformed (> 90 %)', atlas.acyl.length && Math.max(...atlas.acyl) > 0.9, atlas.acyl);
  check('NaBH4 spares an ester in a keto-ester (bystander < 40 %)', atlas.ester && atlas.ester.pb < 0.4, atlas.ester);

  // ---- the suggestions
  const adv = await page.evaluate(async () => {
    const out = {};
    for (const [k, smi] of [['nitroester', 'CCOC(=O)c1ccc([N+](=O)[O-])cc1'], ['ketoaldehyde', 'CC(=O)CCC=O']]) {
      const r = await window.Advisor.suggest({ subs: [smi], reagent: '', cond: '' });
      out[k] = r.textbook.map(t => ({ id: t.rule.id, plaus: t.plausible, rec: !!t.recorded, opts: t.options.map(o => ({ label: o.label, risk: o.collateral && o.collateral.risk, risky: o.collateral && o.collateral.risky.map(x => x.names[0]), unknown: o.collateral && (o.collateral.unknown || []).map(x => x.names[0]) })), recOpts: t.recorded ? t.recorded.options.map(o => o.label) : [] }));
    }
    return out;
  });
  const nitro = adv.nitroester, kald = adv.ketoaldehyde;
  check('nitro ester: nitro reduction is offered, backed by recorded reactions', nitro.some(t => t.id === 'R023' && t.rec && t.recOpts.some(l => /Pd/i.test(l))), nitro.map(t => t.id));
  check('nitro ester: saponification is offered', nitro.some(t => t.id === 'R047'));
  check('nitro ester: LiAlH4 is NOT called safe for the nitro group (no precedent or interference)', nitro.some(t => t.id === 'R038' && t.opts.some(o => (o.unknown || []).concat(o.risky || []).some(x => /nitro/i.test(x)))), nitro.filter(t => t.id === 'R038'));
  check('keto-aldehyde: NaBH4 is a top idea', kald.slice(0, 4).some(t => t.id === 'R037'), kald.map(t => t.id));
  check('keto-aldehyde: n-BuLi / hydroxylamine are flagged as also hitting the aldehyde', kald.filter(t => t.id === 'R076' || t.id === 'R094').every(t => t.opts.some(o => (o.risky || []).some(x => /aldehyde/i.test(x)))), kald.filter(t => t.id === 'R076' || t.id === 'R094'));

  // ---- the panel and the Use button
  // the equation block is test-only detail now (the page shows the verdict and the mechanism button): ask for it
  await page.evaluate(() => window.ReactionsTests.detail(true));
  await page.locator('#rxn-code').click(); await page.locator('#rxn-code').fill('CC(=O)CCC=O'); await page.locator('#rxn-code').press('Enter');
  await page.waitForTimeout(1200);
  await page.locator('#rxn-suggest').click();
  await page.waitForFunction(() => /ideas/.test((document.querySelector('.adv-status') || {}).textContent || ''), null, { timeout: 120000 });
  const nCards = await page.locator('.adv-card').count();
  check('the Suggest panel lists ideas', nCards >= 5, nCards);
  await page.locator('.adv-card', { hasText: 'Ketone reduction (NaBH4)' }).first().locator('.adv-use').first().click();
  await page.waitForFunction(() => { const v = document.querySelector('.rxn-verdict'); return v && v.textContent && !v.textContent.startsWith('Running'); }, null, { timeout: 90000 });
  await page.waitForTimeout(1500);
  const used = { bar: await page.locator('#rxn-code').inputValue(), prod: await page.locator('#rxn-prod').inputValue() };
  check('Use puts the reagent on the page and runs it (aldehyde reduced, ketone kept)', /BH4/.test(used.bar) && /CC\(=O\)CCCO/.test(used.prod), used);
  await page.waitForFunction(() => document.querySelector('.rxn-eqn'), null, { timeout: 30000 }).catch(() => {});
  const eq = await page.locator('.rxn-eqn').first().innerText().catch(() => '');
  check('the equation block balances the reduction with H2 and gives a table', /H₂/.test(eq) && /mmol/.test(eq), eq.slice(0, 300));

  // ---- equation balancer in the page (esterification needs water)
  const esterEq = await page.evaluate(() => { const S = window.Stoich; const sp = (n, f) => Object.assign({ name: n }, S.formula(f)); const r = S.balance([sp('a', 'C2H4O2'), sp('b', 'C2H6O')], [sp('c', 'C4H8O2')]); return r.ok ? r.added.map(a => a.name + ' ' + a.side) : 'fail'; });
  check('Stoich: esterification closes with water on the right', JSON.stringify(esterEq) === JSON.stringify(['H₂O right']), esterEq);

  // ---- the one-line verdict and "Open Mechanism in Draw" (Draw.showMechanism is watched; the real one, when loaded, still draws)
  await page.evaluate(() => window.ReactionsTests.detail(false));
  const mechRun = async code => {
    await page.evaluate(c => window.Reactions.load(c), code);
    await page.evaluate(() => document.getElementById('rxn-arrow').dispatchEvent(new MouseEvent('click', { bubbles: true })));
    await page.waitForFunction(() => { const v = document.querySelector('#rxn-rule .rxn-verdict'); return v && v.textContent && !v.textContent.startsWith('Running'); }, null, { timeout: 120000 });
    await page.waitForTimeout(500);
    return page.evaluate(() => { const n = document.getElementById('rxn-rule'), b = n.querySelector('.rxn-mech-open');
      return { kids: [...n.children].map(c => c.className), verdict: n.querySelector('.rxn-verdict').textContent, enabled: !!b && !b.disabled, title: b ? b.title : '', text: b ? b.textContent : '' }; });
  };
  await page.evaluate(() => { window.Draw = window.Draw || {}; const own = window.Draw.showMechanism;
    window.Draw.showMechanism = p => { window.__mechPayload = p; return typeof own === 'function' ? own.call(window.Draw, p) : true; }; });
  const m4 = await mechRun('C=CC>Br> | cond=light');
  check('Reactions shows only the verdict and the Open Mechanism in Draw button', JSON.stringify(m4.kids) === JSON.stringify(['rxn-verdict', 'rxn-mech-open']) && m4.text === 'Open Mechanism in Draw' && m4.enabled, m4);
  await page.locator('#rxn-rule .rxn-mech-open').click(); await page.waitForTimeout(300);
  const p4 = await page.evaluate(() => { const p = window.__mechPayload; let clone = false; try { clone = JSON.stringify(structuredClone(p)) === JSON.stringify(p); } catch (e) { clone = false; }
    return p && { title: p.title, ruleId: p.ruleId, onExample: p.onExample, steps: p.steps.map(s => [s.step, !!s.before, !!s.arrows, !!s.graph]), clone, tab: (document.querySelector('.tab.active') || {}).dataset.page }; });
  check('R004 mechanism goes to Draw: 4 workbook steps with arrows (initiation is two steps), on the Draw tab', p4 && p4.title === 'Mechanism for Radical hydrobromination (anti-Markovnikov)' && p4.ruleId === 'R004' && p4.onExample === true &&
    p4.steps.length === 4 && p4.steps.every(s => s[1] && s[2] && !s[3]) && p4.clone && p4.tab === 'draw', p4);
  await page.evaluate(() => { window.__mechPayload = null; });
  const mSN2 = await mechRun('CCC(C)Br>[C-]#N.[Na+]> | cond=DMSO');
  await page.locator('#rxn-rule .rxn-mech-open').click(); await page.waitForTimeout(300);
  const pSN2 = await page.evaluate(() => { const p = window.__mechPayload; let clone = false; try { clone = JSON.stringify(structuredClone(p)) === JSON.stringify(p); } catch (e) { clone = false; }
    return p && { title: p.title, onExample: p.onExample, graphs: p.steps.every(s => s.graph && s.graph.atoms.length && s.gArrows.length), clone, reaction: p.reaction }; });
  check('SN2 (engine) mechanism goes to Draw on the user\'s molecules: graph steps with arrows', /^SN2 /.test(mSN2.verdict) && pSN2 && pSN2.title === 'Mechanism for SN2 substitution' && pSN2.onExample === false && pSN2.graphs && pSN2.clone &&
    pSN2.reaction.reactants[0] === 'CCC(C)Br' && pSN2.reaction.conditions === 'DMSO', pSN2);
  const mL = await mechRun('NCCc1ccccc1[N+](=O)[O-]>[Pd]> | cond=ethanol');
  check('a learned answer has no mechanism: the button is disabled and says why', /L0000001/.test(mL.verdict) && !mL.enabled && /no mechanism/i.test(mL.title), mL);

  // ---- retrosynthesis from recorded reactions
  await page.locator('.tab[data-page="retro"]').click();
  await page.locator('#retroTarget').fill('CC(=O)Nc1ccc(O)cc1');
  await page.locator('#retroGo').click();
  await page.waitForFunction(() => { const s = document.querySelector('.rl-status'); return s && /disconnections, best first\. A proposal/.test(s.textContent); }, null, { timeout: 240000 });
  const cards = await page.locator('.rl-card').allInnerTexts();
  check('retro: paracetamol has an amide disconnection to 4-aminophenol and a buyable acid or acid chloride', cards.slice(0, 4).some(t => /4-aminophenol/i.test(t) && /(acetic acid|acetyl chloride|ethanoic anhydride|ethanoyl)/i.test(t)), cards.slice(0, 3).map(t => t.slice(0, 120)));
  check('retro: each disconnection shows recorded reagents', cards.slice(0, 3).every(t => /recorded/.test(t)));
  await page.locator('#retroTarget').fill('CCN(CC)CC(=O)Nc1c(C)cccc1C');
  await page.locator('#retroGo').click();
  await page.waitForFunction(() => { const s = document.querySelector('.rl-status'); return s && /disconnections, best first\. A proposal/.test(s.textContent); }, null, { timeout: 240000 });
  await page.locator('.rl-btns button', { hasText: 'Plan a whole route' }).click();
  await page.waitForFunction(() => document.querySelector('.rl-sec .card.verdict'), null, { timeout: 300000 });
  const route = await page.locator('.rl-sec .card.verdict').first().innerText();
  check('retro: lidocaine route has no piece without a route', /reactions?/.test(route) && !/no recorded route|no taught route/.test(route), route.slice(0, 200));

  // ---- cross-app links
  await page.locator('.tab[data-page="sketcher"]').click();
  await page.locator('#sketchImport').fill('c1ccccc1O');
  await page.locator('[data-act="import"]').click(); await page.waitForTimeout(600);
  await page.locator('[data-act="toReactions"]').click(); await page.waitForTimeout(2000);
  check('Sketcher → Reactions puts the drawing in a box', /Oc1ccccc1|c1ccc\(O\)cc1|c1ccccc1O/.test(await page.locator('#rxn-code').inputValue()), await page.locator('#rxn-code').inputValue());
  const ok = await page.evaluate(async () => window.Reactions.load('CC(=O)O.CCO>O=S(=O)(O)O>CCOC(C)=O | conc=98%w/w | cond=reflux, ethanol'));
  check('Reactions.load fills the page from a code', ok === true);
  check('no page error', errors.length === 0, errors.slice(0, 3));
  console.log(fails ? `${fails} of ${n} FAILED` : `all ${n} passed`);
  await browser.close();
  process.exit(fails ? 1 : 0);
})().catch(e => { console.error(e); process.exit(2); });
