/* RxnCode — the reaction code on the Reactions page, version 2.

     reactants>agents>products | conc=37%w/w | cond=reflux, ethanol, 2 equiv

   The first word is an ordinary reaction SMILES (reactants>agents>products, dot-separated molecules in the SMILES-line order), so any
   chemistry tool that reads the first word of a line (RDKit's ReactionFromSmarts(line.split()[0]), the USPTO / Reaxys / ORD exports) can
   use the code as it stands. What a reaction SMILES cannot say follows after a space, as  key=value  items joined by ' | ':

     conc=<strength of the reagent>   37%w/w, 5%v/v, 3%w/v, 2M, or the words conc / dil     (js/reagent_conc.js)
     cond=<the conditions box>        the comma list of the conditions box in its standard spelling ("reflux, ethanol, 2 equiv", "-78 C")
                                      (amounts are written here as "2 equiv", "excess", "catalytic")

   The agents field holds the reagent line (reagents and catalysts); solvents and temperatures stay in cond because the page reads them as
   words (reflux in ethanol), not as molecules. The products field is empty until the arrow has been clicked.

   Version 1 (still read, never written):   starting materials~reagent~conditions~products   (at most three '~').
   '~' never occurs in a SMILES, so neither format needs quoting. A '>' separates the three fields except in '->', a dative bond (a
   SMILES never ends in a bond, so '->' is never a separator; RDKit reads 'CC>[NH3]->[Pd]>CC' the same way). tools/rxn_code.py is the
   same definition in Python, for the tools that build the data.

   RxnCode.parse(text) -> { ok, version, reactants[], agents[], products[], cond, conc, extra{}, error }
   RxnCode.build({reactants, agents, products, cond, conc}) -> string   (the molecule fields are already-joined SMILES lines)
   RxnCode.toV2(text) -> the same reaction written as version 2 (null if it cannot be read)
   RxnCode.reactionSmiles(text) -> the first word, for copying into another program                                                       */
(function (root) {
  'use strict';
  const split = s => String(s || '').split('.').map(x => x.trim()).filter(Boolean);

  function parse(text) {
    const t = String(text || '').trim();
    if (!t) return { ok: false, error: 'empty' };
    // version 1: a '~' and no '>' in the first word
    const first = t.split(/\s+/)[0];
    if (!/(?<!-)>/.test(first) && t.includes('~')) {
      const parts = t.split('~').map(x => x.trim());
      if (parts.length > 4) return { ok: false, error: 'A reaction code is starting materials~reagent~conditions~products (at most three ~)' };
      while (parts.length < 4) parts.push('');
      return { ok: true, version: 1, reactants: split(parts[0]), agents: split(parts[1]), products: split(parts[3]), cond: parts[2], conc: null, extra: {} };
    }
    const rxn = first.split(/(?<!-)>/);                 // not the '>' of a dative bond '->'
    if (rxn.length === 1) {                              // a bare molecule list: starting materials only
      return { ok: true, version: 2, reactants: split(first), agents: [], products: [], cond: '', conc: null, extra: parseTail(t.slice(first.length)).extra, ...tailFields(t.slice(first.length)) };
    }
    if (rxn.length !== 3) return { ok: false, error: 'A reaction SMILES has exactly two ">" (reactants>agents>products)' };
    return { ok: true, version: 2, reactants: split(rxn[0]), agents: split(rxn[1]), products: split(rxn[2]), cond: '', conc: null, extra: {}, ...tailFields(t.slice(first.length)) };
  }

  function parseTail(tail) {
    const items = {};
    String(tail || '').split('|').map(x => x.trim()).filter(Boolean).forEach(seg => {
      const i = seg.indexOf('=');
      if (i > 0) items[seg.slice(0, i).trim().toLowerCase()] = seg.slice(i + 1).trim();
    });
    return { items, extra: {} };
  }
  function tailFields(tail) {
    const { items } = parseTail(tail), extra = {};
    let cond = '', conc = null;
    for (const [k, v] of Object.entries(items)) {
      if (k === 'cond') cond = v;
      else if (k === 'conc') conc = v;                  // parsed by ReagentConc where it is used (this file knows no units)
      else extra[k] = v;
    }
    return { cond, conc, extra };
  }

  function build(o) {
    const rxn = `${o.reactants || ''}>${o.agents || ''}>${o.products || ''}`;
    if (rxn === '>>' && !o.cond && !o.conc) return '';
    const tail = [];
    if (o.conc) tail.push('conc=' + o.conc);
    if (o.cond) tail.push('cond=' + o.cond);
    return rxn + (tail.length ? ' | ' + tail.join(' | ') : '');
  }

  const reactionSmiles = text => String(text || '').trim().split(/\s+/)[0] || '';

  /* a version 1 code (or any code) rewritten as version 2: the same molecules, conditions and strength */
  function toV2(text) {
    const p = parse(text);
    if (!p.ok) return null;
    return build({ reactants: p.reactants.join('.'), agents: p.agents.join('.'), products: p.products.join('.'), cond: p.cond, conc: p.conc });
  }

  const api = { parse, build, reactionSmiles, toV2 };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.RxnCode = api;
})(typeof window !== 'undefined' ? window : globalThis);
