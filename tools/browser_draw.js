// Draw tab checks, in a headless browser: molecule boxes, arrows, text, colour, resize, turn, captions, PNG export.
//   python tools/serve.py   (another terminal)        node tools/browser_draw.js [http://localhost:8000/app.html]
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const base = process.argv[2] || 'http://localhost:8000/app.html'; const out = require('os').tmpdir(); let fails = 0;
const check = (n, ok, d) => { if (!ok) { fails++; console.log('FAIL', n, d === undefined ? '' : JSON.stringify(d).slice(0, 300)); } else console.log('ok  ', n); };
(async () => {
  const b = await chromium.launch({ executablePath: process.env.CHROMIUM || undefined, args: ['--no-sandbox'] });
  const ctx = await b.newContext({ viewport: { width: 1400, height: 900 }, acceptDownloads: true });
  const p = await ctx.newPage(); const errs = []; p.on('pageerror', e => errs.push(e.message)); p.on('console', m => { if (m.type() === 'error') errs.push(m.text()); });
  await p.goto(base + '#draw', { waitUntil: 'load' });
  await p.waitForFunction(() => window.Draw, null, { timeout: 30000 });
  check('one tab named Draw', await p.locator('.tab', { hasText: /^Draw$/ }).count() === 1 && await p.locator('.tab', { hasText: /^Draw1$/ }).count() === 0);
  const svg = p.locator('.dr-scroll svg');
  const origin = await svg.boundingBox();
  const at = (x, y) => [origin.x + x, origin.y + y];
  const btn = t => p.locator('.dr-bar button', { hasText: t }).first();
  const drag = async (x1, y1, x2, y2, mods) => { await p.mouse.move(...at(x1, y1)); await p.mouse.down(); await p.mouse.move(...at((x1 + x2) / 2, (y1 + y2) / 2), { steps: 3 }); await p.mouse.move(...at(x2, y2), { steps: 3 }); await p.mouse.up(); };
  const importInto = async text => { await p.locator('.dr-dlg input[type=text]').fill(text); await p.locator('.dr-dlg [data-import]').click(); await p.waitForTimeout(150); };

  // molecule box 1: ethanol
  await btn('Molecule box').click(); await drag(60, 80, 300, 240);
  check('dialog opens after dragging a molecule box', await p.locator('.dr-dlg').count() === 1);
  await importInto('ethanol');
  check('ethanol imported into the box', await p.evaluate(() => Draw.items[0].smiles) === 'CCO', await p.evaluate(() => Draw.items[0].smiles));
  await p.locator('.dr-dlg [data-done]').click();
  // box 2: aspirin by name, box 3: benzene by SMILES
  await btn('Molecule box').click(); await drag(520, 80, 800, 260);
  await importInto('aspirin'); await p.locator('.dr-dlg [data-done]').click();
  await btn('Molecule box').click(); await drag(560, 340, 780, 500);
  await importInto('c1ccccc1'); await p.locator('.dr-dlg [data-done]').click();
  const info = await p.evaluate(() => Draw.items.map(i => ({ s: i.smiles, w: Math.round(i.w), h: Math.round(i.h), atoms: i.atomXY.length })));
  check('three molecule boxes with their molecules', info.length === 3 && info[1].atoms === 13 && info[2].atoms === 6, info);
  check('box auto-fits the molecule', info[2].w < 260 && info[2].h < 200, info[2]);
  // arrow
  await btn('→').click(); await drag(310, 160, 500, 168);
  const ar = await p.evaluate(() => { const a = Draw.items[3]; return { kind: a.kind, dy: Math.abs(a.p2.y - a.p1.y) }; });
  check('arrow snaps straight near horizontal', ar.kind === 'forward' && ar.dy < 0.5, ar);
  await btn('⇌').click(); await drag(310, 400, 500, 400);
  check('equilibrium arrow draws two lines', await p.evaluate(() => Draw.items[4].g.querySelectorAll('path').length) === 2);
  // electron + lone pair
  await btn('••').click(); await drag(120, 330, 140, 330);
  check('lone pair placed, turned across the drag', await p.evaluate(() => Draw.items[5].pair && Math.abs(Draw.items[5].angle - Math.PI / 2) < 0.3));
  await btn('•').click(); await drag(200, 330, 200, 330);
  // text with a subscript
  await btn(/^Text$/).click(); await p.mouse.click(...at(60, 420));
  await p.keyboard.type('H'); await btn('x₂').click(); await p.keyboard.type('2'); await btn('x₂').click(); await p.keyboard.type('O');
  await p.mouse.click(...at(900, 700));
  const tb = await p.evaluate(() => { const t = Draw.items.find(i => i instanceof Draw.TextBox); return t && t.runs; });
  check('text box keeps a subscript run', tb && tb.some(r => r.s === -1 && r.t === '2') && tb.map(r => r.t).join('') === 'H2O', tb);
  // select by rectangle, colour
  await p.mouse.move(...at(520, 60)); await p.mouse.down(); await p.mouse.move(...at(820, 280), { steps: 4 }); await p.mouse.up();
  check('rectangle selects the aspirin box', await p.evaluate(() => Draw.selection.length === 1 && Draw.selection[0].smiles.includes('C(=O)')), await p.evaluate(() => Draw.selection.length));
  await btn('Color').click(); await p.locator('.dr-color i[data-c="#e03131"]').first().click();
  check('colour applied to the molecule', await p.evaluate(() => Draw.selection[0].atomColors.size > 0 && [...Draw.selection[0].atomColors.values()][0] === '#e03131'));
  await p.mouse.click(...at(900, 650));
  // atom-level lasso: only the carboxyl part of aspirin
  const aspirin = await p.evaluate(() => { const i = Draw.items[1]; return { x: i.x, y: i.y, xy: i.atomXY.map(q => [q.x, q.y]) }; });
  // resize: select the benzene box, drag its right edge
  await p.mouse.click(...at(670, 420));
  const before = await p.evaluate(() => Draw.items[2].w);
  const bx = await p.evaluate(() => { const i = Draw.items[2]; return { r: i.x + i.w, y: i.y + i.h / 2 }; });
  await drag(bx.r, bx.y, bx.r + 60, bx.y);
  const after = await p.evaluate(() => Draw.items[2].w);
  check('selected box resizes from its border', after > before + 40, [before, after]);
  // rotate: drag in the band just outside the molecule
  const geo = await p.evaluate(() => { const i = Draw.items[2]; const xs = i.atomXY.map(q => q.x), ys = i.atomXY.map(q => q.y); return { x: i.x, y: i.y, cx: i.x + (Math.min(...xs) + Math.max(...xs)) / 2, cy: i.y + (Math.min(...ys) + Math.max(...ys)) / 2, r: Math.max(...xs) - Math.min(...xs) }; });
  await drag(geo.cx + geo.r / 2 + 22, geo.cy, geo.cx, geo.cy + geo.r / 2 + 40);
  const ang = await p.evaluate(() => Draw.items[2].angle);
  check('dragging beside the molecule turns it', Math.abs(ang) > 20, ang);
  // display options
  await p.mouse.dblclick(...at(geo.cx, geo.cy));
  await p.locator('.dr-dlg button', { hasText: 'IUPAC' }).click(); await p.locator('.dr-dlg button', { hasText: 'SMILES' }).click(); await p.locator('.dr-dlg button', { hasText: 'Show H' }).click();
  await p.locator('.dr-dlg [data-done]').click();
  const cap = await p.evaluate(() => Draw.items[2].captionLines());
  check('captions show IUPAC name and SMILES', cap.length === 2 && /benzene/i.test(cap[0]), cap);
  // delete + Esc
  await p.mouse.click(...at(160, 160)); await p.keyboard.press('Delete');
  check('Delete removes the selected molecule box', await p.evaluate(() => Draw.items.length) === 7, await p.evaluate(() => Draw.items.length));
  await p.screenshot({ path: out + '/d1.png' });
  const [dl] = await Promise.all([p.waitForEvent('download'), btn('Export PNG').click()]);
  const path = await dl.path(); const fs = require('fs'); const buf = fs.readFileSync(path);
  check('Export PNG saves a PNG', buf.slice(1, 4).toString() === 'PNG' && buf.length > 5000, buf.length);
  check('no page errors', errs.length === 0, errs);
  await b.close(); console.log(fails ? fails + ' FAILED' : 'all passed'); process.exit(fails ? 1 : 0);
})();
