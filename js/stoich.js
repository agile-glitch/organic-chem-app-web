/* Stoich — balance a reaction exactly, and turn it into a table of amounts.

   Method (the standard algebraic one): every species is a column of element counts, with charge as one more row; the reactants count +, the products −;
   a balanced equation is a vector of positive whole numbers in the null space of that matrix. The null space is found by Gauss-Jordan elimination in exact
   rational arithmetic (BigInt fractions: no rounding, no floating-point surprises), so any equation that can be balanced is, ionic and redox equations
   included (charge is a row; electrons never need to be written).

   When the species given do not close (a hydrolysis that needs water, a reduction that needs H2, a condensation that loses HCl), the smallest set of small
   molecules that does close it is searched: none, then one, then two of  H2O, H2, O2, HCl, HBr, HI, HF, CO2, N2, NH3, H2O2, H+, OH-.  Each may land on either
   side. The answer names what was added, so nothing is assumed silently; and says when the equation cannot be balanced or has more than one solution.

   Stoich.formula(text)                  'C6H12O6' / 'Ca(OH)2' / 'SO4 2-' → {El: n}, charge
   Stoich.balance(left, right, opts)     left, right = [{name, comp:{El:n}, charge?}]  → {ok, coefficients:{left[], right[]}, added:[{name, side, coef}], multiple, note}
   Stoich.table(result, species, o)      amounts for a scale: o = {limiting: index in left, mmol}  → rows with mmol, mg, equivalents, theoretical yield
   Stoich.atomEconomy(rows)              mass of the wanted product / mass of everything used                                                          */
(function (root) {
  'use strict';
  const gcd = (a, b) => { a = a < 0n ? -a : a; b = b < 0n ? -b : b; while (b) { [a, b] = [b, a % b]; } return a; };
  const lcm = (a, b) => (a === 0n || b === 0n ? 0n : (a / gcd(a, b)) * b);

  class Frac {
    constructor(n, d = 1n) { if (d < 0n) { n = -n; d = -d; } const g = gcd(n, d) || 1n; this.n = n / g; this.d = d / g; }
    add(o) { return new Frac(this.n * o.d + o.n * this.d, this.d * o.d); }
    sub(o) { return new Frac(this.n * o.d - o.n * this.d, this.d * o.d); }
    mul(o) { return new Frac(this.n * o.n, this.d * o.d); }
    div(o) { return new Frac(this.n * o.d, this.d * o.n); }
    get zero() { return this.n === 0n; }
  }

  /* the null space of an integer matrix (rows x columns, BigInt entries): a list of integer vectors with gcd 1 */
  function nullspace(M, cols) {
    const A = M.map(r => r.map(x => new Frac(BigInt(x))));
    const pivots = [];
    let row = 0;
    for (let c = 0; c < cols && row < A.length; c++) {
      let p = -1;
      for (let r = row; r < A.length; r++) if (!A[r][c].zero) { p = r; break; }
      if (p < 0) continue;
      [A[row], A[p]] = [A[p], A[row]];
      const inv = A[row][c];
      A[row] = A[row].map(x => x.div(inv));
      for (let r = 0; r < A.length; r++) {
        if (r === row || A[r][c].zero) continue;
        const f = A[r][c];
        A[r] = A[r].map((x, j) => x.sub(f.mul(A[row][j])));
      }
      pivots.push(c); row++;
    }
    const free = [];
    for (let c = 0; c < cols; c++) if (!pivots.includes(c)) free.push(c);
    return free.map(fc => {
      const v = new Array(cols).fill(null).map(() => new Frac(0n));
      v[fc] = new Frac(1n);
      pivots.forEach((pc, r) => { v[pc] = new Frac(0n).sub(A[r][fc]); });
      let L = 1n;
      v.forEach(x => { L = lcm(L, x.d); });
      const ints = v.map(x => (x.n * (L / x.d)));
      let g = 0n; ints.forEach(x => { g = gcd(g, x); });
      return ints.map(x => x / (g || 1n));
    });
  }

  // ---------------------------------------------------------------- formulas
  function formula(text) {
    let s = String(text).trim(), charge = 0, m;
    // a charge: "Fe3+" (a single-element ion: the digits are the charge), "SO4 2-" or "SO4^2-" (a space or ^ before the number), "MnO4-" / "NH4+" (the sign alone: 1)
    if ((m = s.match(/^([A-Z][a-z]?)(\d*)([+-])$/))) { charge = (m[2] === '' ? 1 : +m[2]) * (m[3] === '+' ? 1 : -1); s = m[1]; }
    else if ((m = s.match(/(?:\s|\^)(\d*)([+-])$/))) { charge = (m[1] === '' ? 1 : +m[1]) * (m[2] === '+' ? 1 : -1); s = s.slice(0, m.index); }
    else if ((m = s.match(/([+-])$/))) { charge = m[1] === '+' ? 1 : -1; s = s.slice(0, m.index); }
    const stack = [{}];
    const re = /([A-Z][a-z]?)(\d*)|(\()|(\))(\d*)/g;
    let t;
    while ((t = re.exec(s))) {
      if (t[1]) { const top = stack[stack.length - 1]; top[t[1]] = (top[t[1]] || 0) + (t[2] === '' ? 1 : +t[2]); }
      else if (t[3]) stack.push({});
      else if (t[4]) { const inner = stack.pop(), k = t[5] === '' ? 1 : +t[5], top = stack[stack.length - 1]; for (const [e, n] of Object.entries(inner)) top[e] = (top[e] || 0) + n * k; }
    }
    return { comp: stack[0], charge };
  }
  const SMALL = [['H₂O', 'H2O'], ['H₂', 'H2'], ['O₂', 'O2'], ['HCl', 'HCl'], ['HBr', 'HBr'], ['HI', 'HI'], ['HF', 'HF'], ['CO₂', 'CO2'], ['N₂', 'N2'], ['NH₃', 'NH3'], ['H₂O₂', 'H2O2'], ['H⁺', 'H+'], ['OH⁻', 'OH-']]
    .map(([name, f]) => Object.assign({ name }, formula(f)));

  // ---------------------------------------------------------------- balancing
  /* try to balance exactly these species: fixed ones (left > 0, right > 0 required) and optional free ones (any sign, or 0). */
  function solve(left, right, free) {
    const elements = [...new Set([...left, ...right, ...free].flatMap(s => Object.keys(s.comp)))].sort();
    const charged = [...left, ...right, ...free].some(s => s.charge);
    const cols = [...left.map(s => ({ s, sign: 1 })), ...right.map(s => ({ s, sign: -1 })), ...free.map(s => ({ s, sign: 1 }))];
    const rows = elements.map(e => cols.map(c => BigInt(c.sign * (c.s.comp[e] || 0))));
    if (charged) rows.push(cols.map(c => BigInt(c.sign * (c.s.charge || 0))));
    const basis = nullspace(rows, cols.length);
    if (!basis.length) return null;
    const nL = left.length, nR = right.length;
    const fine = v => v.slice(0, nL + nR).every(x => x > 0n);
    const frees = v => v.slice(nL + nR);
    let sols = [];
    for (const b of basis) for (const sg of [1n, -1n]) { const v = b.map(x => x * sg); if (fine(v)) sols.push(v); }
    if (!sols.length && basis.length > 1) {                      // a combination of the basis vectors
      const range = [-3n, -2n, -1n, 0n, 1n, 2n, 3n];
      const k = Math.min(basis.length, 3);
      const rec = (i, acc) => {
        if (i === k) { if (acc && fine(acc)) sols.push(acc); return; }
        for (const a of range) rec(i + 1, acc ? acc.map((x, j) => x + a * basis[i][j]) : basis[i].map(x => x * a));
      };
      rec(0, null);
    }
    if (!sols.length) return null;
    // the smallest whole-number solution; free species must not be zero when they were asked for
    const cost = v => v.reduce((t, x) => t + (x < 0n ? -x : x), 0n);
    sols = sols.filter(v => v.some(x => x !== 0n));
    sols.sort((a, b) => (cost(a) < cost(b) ? -1 : cost(a) > cost(b) ? 1 : 0));
    const v = sols[0];
    let g = 0n; v.forEach(x => { g = gcd(g, x); });
    const w = v.map(x => x / (g || 1n));
    return { vector: w, multiple: basis.length > 1 };
  }
  const CANDIDATE_ORDER = ['H₂O', 'H₂', 'HCl', 'HBr', 'CO₂', 'O₂', 'N₂', 'NH₃', 'HI', 'HF', 'H₂O₂', 'H⁺', 'OH⁻'];
  function balance(left, right, opts = {}) {
    if (!left.length || !right.length) return { ok: false, note: 'Both sides need at least one species.' };
    const pool = (opts.pool || SMALL).slice().sort((a, b) => CANDIDATE_ORDER.indexOf(a.name) - CANDIDATE_ORDER.indexOf(b.name));
    const attempts = [[]];
    pool.forEach(c => attempts.push([c]));
    if (opts.maxAdded !== 1) for (let i = 0; i < pool.length; i++) for (let j = i + 1; j < pool.length; j++) attempts.push([pool[i], pool[j]]);
    for (const free of attempts) {
      const sol = solve(left, right, free);
      if (!sol) continue;
      const v = sol.vector, nL = left.length, nR = right.length;
      const added = free.map((c, i) => ({ name: c.name, comp: c.comp, charge: c.charge, coef: v[nL + nR + i] })).filter(a => a.coef !== 0n)
        .map(a => ({ name: a.name, comp: a.comp, charge: a.charge, side: a.coef > 0n ? 'left' : 'right', coef: Number(a.coef > 0n ? a.coef : -a.coef) }));
      if (free.length && !added.length) continue;               // the free species were not needed: an earlier attempt would have found it
      return {
        ok: true, coefficients: { left: v.slice(0, nL).map(Number), right: v.slice(nL, nL + nR).map(Number) }, added, multiple: sol.multiple,
        note: added.length ? `Balanced by adding ${added.map(a => `${a.coef > 1 ? a.coef + ' ' : ''}${a.name} (${a.side})`).join(' and ')}.` : 'Balanced as written.' + (sol.multiple ? ' More than one balanced form exists (independent reactions); the smallest is shown.' : ''),
      };
    }
    return { ok: false, note: 'These species cannot be balanced, with or without small molecules such as H₂O, H₂, HCl or CO₂: the atoms do not add up. Check the formulas.' };
  }

  // ---------------------------------------------------------------- amounts
  /* species: [{name, mw, side: 'left'|'right', coef, role}] in the equation order; o = {limiting: index among the left species, mmol: of that species} */
  function table(species, o) {
    const lim = species.filter(s => s.side === 'left')[o.limiting || 0];
    if (!lim) return [];
    const base = (o.mmol == null ? 1 : o.mmol) / lim.coef;
    return species.map(s => {
      const mmol = base * s.coef;
      return Object.assign({}, s, { mmol, mg: mmol * s.mw, equiv: mmol / (base * lim.coef) });
    });
  }
  function atomEconomy(rows, wanted) {
    const used = rows.filter(r => r.side === 'left').reduce((t, r) => t + r.coef * r.mw, 0);
    const w = rows.find(r => r.side === 'right' && (wanted == null || r.name === wanted)) || null;
    return w && used ? (w.coef * w.mw) / used : null;
  }

  const api = { formula, balance, table, atomEconomy, nullspace, SMALL };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.Stoich = api;
})(typeof window !== 'undefined' ? window : globalThis);
