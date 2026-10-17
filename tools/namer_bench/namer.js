// The app's IUPAC namer (js/chem.js), run in Node.js for tools/namer_bench/bench.py.
// One JSON request per line on standard input: {"id": ..., "molfile": ...} or {"id": ..., "smiles": ...}
// One JSON answer per line on standard output: {"id": ..., "name": ...} or {"id": ..., "error": ...}
const fs = require('fs'), vm = require('vm'), path = require('path'), readline = require('readline');
const ctx = { window: {}, console: { log() {}, warn() {}, error() {}, info() {} }, Math, JSON, Date, setTimeout, clearTimeout,
              performance: { now: () => Date.now() } };
ctx.self = ctx.globalThis = ctx;
vm.createContext(ctx);
const CHEM = process.env.CHEM_JS || path.join(__dirname, '..', '..', 'js', 'chem.js');
vm.runInContext(fs.readFileSync(CHEM, 'utf8'), ctx, { filename: 'chem.js' });
const C = ctx.window.Chem;

readline.createInterface({ input: process.stdin }).on('line', line => {
  let req = {}, out;
  try {
    req = JSON.parse(line);
    const g = req.molfile ? C.parseMolfile(req.molfile) : C.parseSmiles(req.smiles);
    const name = C.iupacName(g);
    out = (!name || C.nameFailed(name)) ? { id: req.id, error: name || 'no name' } : { id: req.id, name };
  } catch (e) { out = { id: req.id, error: 'crash: ' + String(e && e.message || e) }; }
  process.stdout.write(JSON.stringify(out) + '\n');
});
