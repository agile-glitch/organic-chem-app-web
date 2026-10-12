// Unit tests for js/stoich.js (exact balancing, by-product inference, amounts).   node tools/test_stoich.js
const path = require('path');
const S = require(path.join(__dirname, '..', 'js', 'stoich.js'));
let fails = 0, n = 0;
const eq = (name, got, want) => { n++; const ser = x => JSON.stringify(x, (k, v) => (typeof v === 'bigint' ? v.toString() + 'n' : v)), a = ser(got), b = ser(want); if (a !== b) { fails++; console.log('FAIL', name, '\n   got ', a, '\n   want', b); } };
const sp = (name, f) => Object.assign({ name }, S.formula(f));
const bal = (L, R, o) => S.balance(L.map(x => sp(...x)), R.map(x => sp(...x)), o);
const coefs = r => r.ok ? [r.coefficients.left, r.coefficients.right, r.added.map(a => `${a.coef} ${a.name} ${a.side}`)] : 'FAIL';

// formulas
eq('formula', S.formula('Ca(OH)2'), { comp: { Ca: 1, O: 2, H: 2 }, charge: 0 });
eq('formula charge', S.formula('SO4 2-'), { comp: { S: 1, O: 4 }, charge: -2 });
eq('formula cation', S.formula('Fe3+'), { comp: { Fe: 1 }, charge: 3 });
eq('formula anion', S.formula('MnO4-'), { comp: { Mn: 1, O: 4 }, charge: -1 });
eq('formula nested', S.formula('Al2(SO4)3').comp, { Al: 2, S: 3, O: 12 });
// classics
eq('methane', coefs(bal([['CH4', 'CH4'], ['O2', 'O2']], [['CO2', 'CO2'], ['H2O', 'H2O']])), [[1, 2], [1, 2], []]);
eq('propane', coefs(bal([['C3H8', 'C3H8'], ['O2', 'O2']], [['CO2', 'CO2'], ['H2O', 'H2O']])), [[1, 5], [3, 4], []]);
eq('permanganate + HCl', coefs(bal([['KMnO4', 'KMnO4'], ['HCl', 'HCl']], [['KCl', 'KCl'], ['MnCl2', 'MnCl2'], ['Cl2', 'Cl2'], ['H2O', 'H2O']])), [[2, 16], [2, 2, 5, 8], []]);
eq('ionic redox (charge row)', coefs(bal([['MnO4-', 'MnO4-'], ['Fe2+', 'Fe2+'], ['H+', 'H+']], [['Mn2+', 'Mn2+'], ['Fe3+', 'Fe3+'], ['H2O', 'H2O']])), [[1, 5, 8], [1, 5, 4], []]);
eq('dichromate redox', coefs(bal([['Cr2O7 2-', 'Cr2O7 2-'], ['Fe2+', 'Fe2+'], ['H+', 'H+']], [['Cr3+', 'Cr3+'], ['Fe3+', 'Fe3+'], ['H2O', 'H2O']])), [[1, 6, 14], [2, 6, 7], []]);
eq('SN2', coefs(bal([['CH3Br', 'CH3Br'], ['OH-', 'OH-']], [['CH3OH', 'CH4O'], ['Br-', 'Br-']])), [[1, 1], [1, 1], []]);
eq('thermite', coefs(bal([['Fe2O3', 'Fe2O3'], ['Al', 'Al']], [['Al2O3', 'Al2O3'], ['Fe', 'Fe']])), [[1, 2], [1, 2], []]);
// organic: the missing small molecule is found and named
eq('esterification adds water', coefs(bal([['acetic acid', 'C2H4O2'], ['ethanol', 'C2H6O']], [['ethyl acetate', 'C4H8O2']])), [[1, 1], [1], ['1 H₂O right']]);
eq('amide coupling adds water', coefs(bal([['benzoic acid', 'C7H6O2'], ['methylamine', 'CH5N']], [['N-methylbenzamide', 'C8H9NO']])), [[1, 1], [1], ['1 H₂O right']]);
eq('ketone reduction needs H2', coefs(bal([['acetone', 'C3H6O']], [['isopropanol', 'C3H8O']])), [[1], [1], ['1 H₂ left']]);
eq('dehydration loses water', coefs(bal([['cyclohexanol', 'C6H12O']], [['cyclohexene', 'C6H10']])), [[1], [1], ['1 H₂O right']]);
eq('acyl chloride + alcohol loses HCl', coefs(bal([['acetyl chloride', 'C2H3ClO'], ['ethanol', 'C2H6O']], [['ethyl acetate', 'C4H8O2']])), [[1, 1], [1], ['1 HCl right']]);
eq('hydrolysis needs water', coefs(bal([['ethyl acetate', 'C4H8O2']], [['acetic acid', 'C2H4O2'], ['ethanol', 'C2H6O']])), [[1], [1, 1], ['1 H₂O left']]);
eq('alkene + Br2 needs nothing', coefs(bal([['propene', 'C3H6'], ['Br2', 'Br2']], [['1,2-dibromopropane', 'C3H6Br2']])), [[1, 1], [1], []]);
eq('decarboxylation loses CO2', coefs(bal([['malonic acid', 'C3H4O4']], [['acetic acid', 'C2H4O2']])), [[1], [1], ['1 CO₂ right']]);
eq('already balanced stays', S.balance([sp('A', 'C2H6O')], [sp('B', 'C2H6O')]).note, 'Balanced as written.');
// impossible and ambiguous
eq('impossible', S.balance([sp('carbon', 'C')], [sp('nitrogen', 'N2')]).ok, false);
eq('impossible message', /cannot be balanced/.test(S.balance([sp('carbon', 'C')], [sp('nitrogen', 'N2')]).note), true);
const amb = bal([['H2', 'H2'], ['O2', 'O2']], [['H2O', 'H2O'], ['H2O2', 'H2O2']]);
eq('independent reactions are reported', [amb.ok, amb.multiple], [true, true]);
eq('independent reactions: still whole, positive numbers', amb.ok && amb.coefficients.left.concat(amb.coefficients.right).every(x => Number.isInteger(x) && x > 0), true);
// the exact null space
eq('nullspace of a rank-deficient matrix', S.nullspace([[1n, 2n, 3n], [2n, 4n, 6n]], 3).length, 2);
eq('nullspace is exact (a big, awkward case)', S.nullspace([[1n, -2n, 0n], [0n, 3n, -7n]], 3), [[14n, 7n, 3n]]);
// amounts
const species = [{ name: 'acetic acid', mw: 60.05, side: 'left', coef: 1 }, { name: 'ethanol', mw: 46.07, side: 'left', coef: 1 }, { name: 'ethyl acetate', mw: 88.11, side: 'right', coef: 1 }, { name: 'water', mw: 18.015, side: 'right', coef: 1 }];
const t = S.table(species, { limiting: 0, mmol: 10 });
eq('table mmol', t.map(r => +r.mmol.toFixed(3)), [10, 10, 10, 10]);
eq('table mg of the product', +t[2].mg.toFixed(1), 881.1);
eq('table equiv', t.map(r => +r.equiv.toFixed(2)), [1, 1, 1, 1]);
eq('atom economy', +(S.atomEconomy(t, 'ethyl acetate') * 100).toFixed(1), 83.0);
eq('2:1 stoichiometry', S.table([{ name: 'a', mw: 10, side: 'left', coef: 2 }, { name: 'b', mw: 5, side: 'left', coef: 1 }, { name: 'c', mw: 20, side: 'right', coef: 1 }], { limiting: 0, mmol: 4 }).map(r => r.mmol), [4, 2, 2]);
console.log(fails ? `${fails} of ${n} FAILED` : `all ${n} passed`);
process.exit(fails ? 1 : 0);
