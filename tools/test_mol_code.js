// Unit tests for js/mol_code.js (molecule code: SDF record with labels, and its compressed AGILES1: form).   node tools/test_mol_code.js
const path = require('path');
const MolCode = require(path.join(__dirname, '..', 'js', 'mol_code.js'));
let fails = 0, n = 0;
const eq = (name, got, want) => { n++; const a = JSON.stringify(got), b = JSON.stringify(want); if (a !== b) { fails++; console.log('FAIL', name, '\n   got ', a, '\n   want', b); } };

// ethanol: C1 C2 O3, then 6 H
const MB = '\n     RDKit          2D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n' +
  '    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n' +
  '    1.2990    0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n' +
  '    2.5981   -0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n' +
  '  1  2  1  0\n  2  3  1  0\nM  END\n';
const N = 9, blank = () => new Array(N).fill('');
const hyb = blank(); hyb[0] = 'sp³'; hyb[1] = 'sp³'; hyb[2] = 'sp³';
const pka = blank(); pka[2] = 'pKa 15.9';
const steric = blank(); steric[2] = 'open';
const partial = blank().map((_, i) => (i < 3 ? ['+0.03', '+0.08', '−0.39'][i] : '+0.05'));
const ids = { name: 'ethanol', iupac: 'ethanol', iupac_checked: 'yes', smiles: 'CCO', inchi: 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3', inchikey: 'LFQSCWFLJHTTHZ-UHFFFAOYSA-N', cas: '64-17-5' };
const sdf = MolCode.build({ molblock: MB, title: 'ethanol', app: '1.0.9 (2026-10-13)', nHeavy: 3, nH: 6, ids,
  labels: { hybridization: hyb, pka, steric, partial, formal: blank() }, angles: [{ a: 1, i: 0, j: 2, text: '109.4°' }], notes: ['a note\non two lines'] });

eq('title on line 1', sdf.split('\n')[0], 'ethanol');
eq('MOL block kept', sdf.includes('  1  2  1  0\n  2  3  1  0\nM  END\n>  <AGILES.format>\nAGILES 1\n'), true);
eq('ends with $$$$', sdf.endsWith('\n$$$$\n'), true);
eq('label lines', sdf.includes('>  <AGILES.labels.pka>\n3 pKa 15.9\n\n'), true);
eq('empty mode is an empty item', sdf.includes('>  <AGILES.labels.formal>\n\n'), true);
eq('label with a space', sdf.includes('>  <AGILES.labels.steric>\n3 open\n'), true);
eq('hasData', [MolCode.hasData(sdf), MolCode.hasData(MB)], [true, false]);

const p = MolCode.parse(sdf);
eq('parse ok', p.ok, true);
eq('atoms', [p.nHeavy, p.nH], [3, 6]);
eq('app', p.app, '1.0.9 (2026-10-13)');
eq('round trip hybridization', p.labels.hybridization, hyb);
eq('round trip pka', p.labels.pka, pka);
eq('round trip partial (H too)', p.labels.partial, partial);
eq('empty mode round trip', p.labels.formal, blank());
eq('mode not written stays absent', 'taft' in p.labels, false);
eq('identifiers round trip', p.ids, ids);
eq('identifier items', sdf.includes('>  <AGILES.inchikey>\nLFQSCWFLJHTTHZ-UHFFFAOYSA-N\n') && sdf.includes('>  <AGILES.cas>\n64-17-5\n'), true);
eq('unknown identifier left out', MolCode.build({ molblock: MB, nHeavy: 3, nH: 6, ids: { smiles: 'CCO', cas: null } }).includes('AGILES.cas'), false);
eq('angles', p.angles, [{ a: 1, i: 0, j: 2, text: '109.4°' }]);
eq('note on one line', p.notes, ['a note on two lines']);
eq('molblock', p.molblock.startsWith('ethanol\n') && p.molblock.endsWith('M  END\n'), true);
eq('CRLF read', MolCode.parse(sdf.replace(/\n/g, '\r\n')).labels.pka, pka);

// errors
eq('no M END', MolCode.parse('hello').ok, false);
eq('plain MOL, no data', MolCode.parse(MB).ok, false);
eq('newer version refused', /version 2/.test(MolCode.parse(sdf.replace('AGILES 1', 'AGILES 2')).error), true);
eq('atom out of range', MolCode.parse(sdf.replace('3 pKa 15.9', '12 pKa 15.9')).ok, false);
eq('build needs M END', (() => { try { MolCode.build({ molblock: 'x' }); return 'built'; } catch (e) { return 'threw'; } })(), 'threw');

(async () => {
  const code = await MolCode.encode(sdf);
  eq('code prefix', code.startsWith('AGILES1:'), true);
  eq('isCode', [MolCode.isCode(code), MolCode.isCode('  ' + code + ' '), MolCode.isCode('CCO'), MolCode.isCode('AGILES1:ab+c')], [true, true, false, false]);
  eq('code shorter than the record', code.length < sdf.length, true);
  eq('decode gives the record back', await MolCode.decode(code), sdf);
  // codes and records written before the rename (ACM1:, ACM.* items, "AgileCHEM molecule code 1") are still read
  const oldSdf = sdf.replace(/<AGILES\./g, '<ACM.').replace('AGILES 1', 'AgileCHEM molecule code 1');
  const oldCode = 'ACM1:' + (await MolCode.encode(oldSdf)).slice('AGILES1:'.length);
  eq('old ACM1: code read', MolCode.isCode(oldCode), true);
  eq('old record parsed', MolCode.parse(await MolCode.decode(oldCode)).labels.pka, pka);
  let err = '';
  try { await MolCode.decode(code.slice(0, 40)); } catch (e) { err = e.message; }
  eq('cut-short code refused', /damaged|cut short|not a molecule/.test(err), true);
  console.log(`${n - fails}/${n} passed` + (fails ? `, ${fails} FAILED` : ''));
  process.exit(fails ? 1 : 0);
})();
