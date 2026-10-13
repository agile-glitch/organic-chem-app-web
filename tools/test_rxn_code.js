// Unit tests for js/rxn_code.js (reaction code v2, the old '~' form it still reads, toV2) and js/reagent_conc.js.   node tools/test_rxn_code.js
const path = require('path');
const RxnCode = require(path.join(__dirname, '..', 'js', 'rxn_code.js'));
const RC = require(path.join(__dirname, '..', 'js', 'reagent_conc.js'));
let fails = 0, n = 0;
const eq = (name, got, want) => { n++; const a = JSON.stringify(got), b = JSON.stringify(want); if (a !== b) { fails++; console.log('FAIL', name, '\n   got ', a, '\n   want', b); } };

// ---- version 2
let p = RxnCode.parse('CC(=O)O.CCO>O=S(=O)(O)O>CCOC(C)=O | conc=98%w/w | cond=reflux, ethanol');
eq('v2 reactants', p.reactants, ['CC(=O)O', 'CCO']); eq('v2 agents', p.agents, ['O=S(=O)(O)O']); eq('v2 products', p.products, ['CCOC(C)=O']);
eq('v2 cond', p.cond, 'reflux, ethanol'); eq('v2 conc text', p.conc, '98%w/w'); eq('v2 version', p.version, 2);
eq('v2 build', RxnCode.build({ reactants: 'CC(=O)O.CCO', agents: 'O=S(=O)(O)O', products: 'CCOC(C)=O', cond: 'reflux, ethanol', conc: '98%w/w' }),
   'CC(=O)O.CCO>O=S(=O)(O)O>CCOC(C)=O | conc=98%w/w | cond=reflux, ethanol');
eq('v2 build, nothing', RxnCode.build({ reactants: '', agents: '', products: '' }), '');
eq('v2 build, no tail', RxnCode.build({ reactants: 'C=C', agents: '', products: '' }), 'C=C>>');
eq('v2 round trip', RxnCode.build({ reactants: p.reactants.join('.'), agents: p.agents.join('.'), products: p.products.join('.'), cond: p.cond, conc: p.conc }),
   'CC(=O)O.CCO>O=S(=O)(O)O>CCOC(C)=O | conc=98%w/w | cond=reflux, ethanol');
// ---- the first word is a reaction SMILES by itself
eq('reaction smiles', RxnCode.reactionSmiles('A.B>C>D | conc=5% | cond=x'), 'A.B>C>D');
// ---- empty fields, spaces
p = RxnCode.parse('CC=O>>');
eq('empty agents/products', [p.reactants, p.agents, p.products], [['CC=O'], [], []]);
p = RxnCode.parse('  C=C>Br>  |  cond=dichloromethane , 0 C  ');
eq('spaces', [p.reactants, p.agents, p.cond], [['C=C'], ['Br'], 'dichloromethane , 0 C']);
p = RxnCode.parse('C=C>Br>BrCCBr | note=hello | cond=rt');
eq('unknown key kept', p.extra, { note: 'hello' });
eq('bare molecules', RxnCode.parse('CCO.CC(=O)O').reactants, ['CCO', 'CC(=O)O']);
eq('too many >', RxnCode.parse('A>B>C>D').ok, false);
eq('empty', RxnCode.parse('   ').ok, false);
// ---- version 1 (read only)
p = RxnCode.parse('CC(=O)O.CCO~O=S(=O)(O)O~reflux in ethanol~');
eq('v1 version', p.version, 1); eq('v1 reactants', p.reactants, ['CC(=O)O', 'CCO']); eq('v1 agents', p.agents, ['O=S(=O)(O)O']); eq('v1 cond', p.cond, 'reflux in ethanol'); eq('v1 products', p.products, []);
p = RxnCode.parse('CC=C~Br~light');
eq('v1 three parts', [p.reactants, p.agents, p.cond, p.products], [['CC=C'], ['Br'], 'light', []]);
eq('v1 five parts', RxnCode.parse('a~b~c~d~e').ok, false);
p = RxnCode.parse('CC(C)=O~[BH4-].[Na+]~methanol, 0 C~CC(C)O');
eq('v1 with products', [p.reactants, p.agents, p.cond, p.products], [['CC(C)=O'], ['[BH4-]', '[Na+]'], 'methanol, 0 C', ['CC(C)O']]);
// ---- a dative bond '->' inside a molecule is not a field separator (RDKit reads it the same way)
p = RxnCode.parse('Cc1ccccc1Br>[NH2]->[Pd]<-[NH2]>Cc1ccccc1N | cond=DMSO');
eq('dative v2', [p.ok, p.reactants, p.agents, p.products, p.cond], [true, ['Cc1ccccc1Br'], ['[NH2]->[Pd]<-[NH2]'], ['Cc1ccccc1N'], 'DMSO']);
p = RxnCode.parse('Cc1ccccc1Br~[NH2]->[Pd]~DMSO~Cc1ccccc1N');
eq('dative v1', [p.version, p.agents, p.cond], [1, ['[NH2]->[Pd]'], 'DMSO']);
// ---- toV2: an old code rewritten with the same fields
eq('toV2', RxnCode.toV2('CC(C)=O~[BH4-].[Na+]~methanol, 0 C~CC(C)O'), 'CC(C)=O>[BH4-].[Na+]>CC(C)O | cond=methanol, 0 C');
eq('toV2 no products', RxnCode.toV2('CC=C~Br~light'), 'CC=C>Br> | cond=light');
eq('toV2 no conditions', RxnCode.toV2('CC=C~Br~~BrCC(C)Br'), 'CC=C>Br>BrCC(C)Br');
eq('toV2 dative', RxnCode.toV2('A~X->[Pd]~DMSO~P'), 'A>X->[Pd]>P | cond=DMSO');
eq('toV2 keeps a v2 code', RxnCode.toV2('C=C>Br>BrCCBr | conc=2M | cond=rt'), 'C=C>Br>BrCCBr | conc=2M | cond=rt');
eq('toV2 unreadable', RxnCode.toV2('a~b~c~d~e'), null);

// ---- concentrations
for (const [t, want] of [['37%', { v: 37, unit: 'ww' }], ['5%v/v', { v: 5, unit: 'vv' }], ['3 % w/v', { v: 3, unit: 'wv' }], ['2M', { v: 2, unit: 'M' }], ['1.5 M', { v: 1.5, unit: 'M' }],
                         ['conc', { cls: 'concentrated' }], ['Conc.', { cls: 'concentrated' }], ['dil', { cls: 'dilute' }], ['dilute', { cls: 'dilute' }],
                         ['0%', null], ['101%', null], ['abc', null], ['', null], ['98%w/w', { v: 98, unit: 'ww' }]])
  eq('parse ' + t, RC.parse(t), want);
for (const c of [{ v: 37, unit: 'ww' }, { v: 5, unit: 'vv' }, { v: 3, unit: 'wv' }, { v: 2, unit: 'M' }, { cls: 'concentrated' }, { cls: 'dilute' }, { v: 12.5, unit: 'wv' }])
  eq('code round trip ' + RC.code(c), RC.parse(RC.code(c)), c);
eq('format', [RC.format({ v: 37, unit: 'ww' }), RC.format({ v: 2, unit: 'M' }), RC.format({ cls: 'dilute' }), RC.short({ v: 37, unit: 'ww' }), RC.short({ cls: 'concentrated' })],
   ['37 % w/w', '2 M', 'dilute', '37%', 'conc.']);
const T = id => RC.TABLE.find(e => e.id === id);
const cls = (id, c) => RC.classify(T(id), c);
eq('H2SO4 98% conc', cls('h2so4', { v: 98, unit: 'ww' }), 'concentrated');
eq('H2SO4 50% neither', cls('h2so4', { v: 50, unit: 'ww' }), 'moderate');
eq('H2SO4 10% dilute', cls('h2so4', { v: 10, unit: 'ww' }), 'dilute');
eq('HCl 37% conc', cls('hcl', { v: 37, unit: 'ww' }), 'concentrated');
eq('HCl 10% dilute', cls('hcl', { v: 10, unit: 'ww' }), 'dilute');
eq('HCl 12 M conc', cls('hcl', { v: 12, unit: 'M' }), 'concentrated');
eq('HCl 1 M dilute', cls('hcl', { v: 1, unit: 'M' }), 'dilute');
eq('NaOH 50% conc', cls('naoh', { v: 50, unit: 'ww' }), 'concentrated');
eq('H2O2 30% conc', cls('h2o2', { v: 30, unit: 'ww' }), 'concentrated');
eq('H2O2 3% dilute', cls('h2o2', { v: 3, unit: 'ww' }), 'dilute');
eq('unknown reagent 100% conc', RC.classify(null, { v: 100, unit: 'ww' }), 'concentrated');
eq('unknown reagent 5% dilute', RC.classify(null, { v: 5, unit: 'ww' }), 'dilute');
eq('qualitative passes through', RC.classify(T('hcl'), { cls: 'dilute' }), 'dilute');
eq('nothing given', RC.classify(T('hcl'), null), null);
// the menu: own strengths first, generic ones after, no duplicates
const pr = RC.presets(T('h2so4'));
eq('first preset', pr[0].label, '98 % w/w'); eq('own flag', pr[0].own, true);
eq('no duplicate labels', new Set(pr.map(x => x.label)).size, pr.length);
eq('generic after own', pr.findIndex(x => !x.own) > pr.findLastIndex(x => x.own), true);
eq('every preset has a class or moderate', pr.every(x => ['concentrated', 'dilute', 'moderate'].includes(x.cls)), true);
// entry lookup through a canonicaliser
const canon = s => ({ '[OH-]': '[OH-]', 'O[Na]': 'O[Na]' }[s] || s);
eq('lookup HCl', (RC.entryFor(canon, 'Cl') || {}).id, 'hcl');
eq('lookup NaOH either order', [(RC.entryFor(canon, '[OH-].[Na+]') || {}).id, (RC.entryFor(canon, '[Na+].[OH-]') || {}).id], ['naoh', 'naoh']);
eq('lookup unknown', RC.entryFor(canon, 'CCO'), null);
eq('mixture Zn + HCl finds HCl', (RC.entryFor(canon, '[Zn].Cl') || {}).id, 'hcl');
eq('mixture HNO3 + H2SO4 finds the first', (RC.entryFor(canon, 'O=[N+]([O-])O.O=S(=O)(O)O') || {}).id, 'h2so4');
eq('NaOH is not found inside Na2CO3-like lines', RC.entryFor(canon, '[Na+].[Na+].O=C([O-])[O-]'), null);
eq('TABLE forms unique', (() => { const seen = new Set(); let dup = 0; RC.TABLE.forEach(e => e.forms.forEach(f => { if (seen.has(f)) dup++; seen.add(f); })); return dup; })(), 0);
console.log(fails ? `${fails} of ${n} FAILED` : `all ${n} passed`);
process.exit(fails ? 1 : 0);
