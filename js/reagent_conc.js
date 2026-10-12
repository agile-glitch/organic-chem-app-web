/* ReagentConc — the strength of a reagent (37 % HCl, 98 % H2SO4, 3 % H2O2 ...), what it is called, and what it counts as.

   A concentration is {v, unit} (unit: 'ww' % w/w, 'vv' % v/v, 'wv' % w/v, 'M' mol/L) or just {cls} ('concentrated' | 'dilute') when the
   user only knows the textbook word. It belongs to the reagent as a whole (the line above the arrow), is written into the reaction code
   as  conc=37%w/w  /  conc=2M  /  conc=conc  /  conc=dil,  and reaches the rules as one of the conditions 'concentrated' / 'dilute'
   (sheet Conditions; tools/add_concentration.py), because what counts as concentrated depends on the reagent:
   37 % is conc. HCl, 98 % is conc. H2SO4, 50 % is conc. NaOH, and 30 % is already conc. H2O2.

   Strengths are the ones sold and written in procedures (sources: supplier specification sheets and the CRC Handbook's concentrated-reagent
   table: HCl 37 %, HNO3 68-70 %, H2SO4 95-98 %, H3PO4 85 %, HBr 48 %, HI 57 %, HF 48 %, NH3 28-30 %, H2O2 30-35 %, formalin 37 %, NaOH/KOH 45-50 %,
   glacial acetic acid 99.7 %). The thresholds are working definitions for the rules, not standards, and every one is listed here to be edited.

   ReagentConc.entryFor(canon, line)   -> the table entry for a reagent line (dot-joined SMILES), or null. canon(smiles) -> canonical SMILES (RDKit).
   ReagentConc.presets(entry)          -> [{label, v, unit, cls, note}]  the menu: this reagent's own strengths, then the generic ones
   ReagentConc.classify(entry, c)      -> 'concentrated' | 'dilute' | 'moderate' | null (null: no concentration given)
   ReagentConc.format(c) / short(c) / code(c) / parse(text)
   ReagentConc.EQUIV                   -> the amounts offered under "Amount" ([{label, text}])                                                       */
(function (root) {
  'use strict';

  /* forms: every way the reagent is typed. concMin / diluteMax in % (w/w, v/v and w/v are treated alike: the table is approximate by nature),
     concMinM / diluteMaxM in mol/L. Outside both: 'moderate'. */
  const TABLE = [
    { id: 'hcl', name: 'hydrochloric acid (HCl)', forms: ['Cl', '[H+].[Cl-]'], concMin: 30, diluteMax: 15, concMinM: 9, diluteMaxM: 4,
      presets: [{ v: 37, unit: 'ww', note: 'conc. HCl, about 12 M (fuming)' }, { v: 20, unit: 'ww', note: 'about 6 M' }, { v: 10, unit: 'ww', note: 'about 3 M' }, { v: 3.6, unit: 'ww', note: 'about 1 M' }],
      why: 'Conc. HCl protonates and chlorinates (alcohol to chloride, Clemmensen); dilute aqueous HCl hydrolyses and keeps water as the nucleophile.' },
    { id: 'h2so4', name: 'sulfuric acid (H2SO4)', forms: ['O=S(=O)(O)O', 'OS(=O)(=O)O'], concMin: 80, diluteMax: 40, concMinM: 14, diluteMaxM: 6,
      presets: [{ v: 98, unit: 'ww', note: 'conc. H2SO4, about 18 M' }, { v: 96, unit: 'ww', note: 'reagent grade' }, { v: 70, unit: 'ww' }, { v: 50, unit: 'ww' }, { v: 20, unit: 'ww', note: 'about 2.5 M' }, { v: 10, unit: 'ww', note: 'about 1 M' }],
      why: 'Conc. H2SO4 dehydrates and nitrates (it takes up water); dilute aqueous H2SO4 hydrates alkenes and hydrolyses esters, acetals and epoxides.' },
    { id: 'hno3', name: 'nitric acid (HNO3)', forms: ['O=[N+]([O-])O', 'O[N+](=O)[O-]'], concMin: 60, diluteMax: 25, concMinM: 12, diluteMaxM: 5,
      presets: [{ v: 70, unit: 'ww', note: 'conc. HNO3, about 16 M' }, { v: 65, unit: 'ww', note: 'about 14 M' }, { v: 90, unit: 'ww', note: 'fuming' }, { v: 10, unit: 'ww', note: 'about 1.6 M' }],
      why: 'Conc. HNO3 (with conc. H2SO4) nitrates arenes; dilute nitric acid is mainly an oxidising acid.' },
    { id: 'h3po4', name: 'phosphoric acid (H3PO4)', forms: ['O=P(O)(O)O', 'OP(=O)(O)O'], concMin: 75, diluteMax: 40, presets: [{ v: 85, unit: 'ww', note: 'conc. (syrupy)' }, { v: 50, unit: 'ww' }, { v: 10, unit: 'ww' }],
      why: 'Conc. H3PO4 dehydrates alcohols without charring.' },
    { id: 'acoh', name: 'acetic acid (AcOH)', forms: ['CC(=O)O'], concMin: 90, diluteMax: 30, concMinM: 15, diluteMaxM: 5,
      presets: [{ v: 100, unit: 'vv', note: 'glacial (99.7 %)' }, { v: 80, unit: 'vv' }, { v: 50, unit: 'vv' }, { v: 10, unit: 'vv' }, { v: 5, unit: 'vv', note: 'vinegar strength' }],
      why: 'Glacial acetic acid is an acidic, polar solvent for brominations and nitrations; aqueous acetic acid is a weak acid.' },
    { id: 'hbr', name: 'hydrobromic acid (HBr)', forms: ['Br', '[H+].[Br-]'], concMin: 40, diluteMax: 20, presets: [{ v: 48, unit: 'ww', note: 'conc. aqueous HBr' }, { v: 33, unit: 'ww', note: 'HBr in acetic acid' }, { v: 10, unit: 'ww' }],
      why: 'Conc. HBr cleaves ethers and makes alkyl bromides from alcohols.' },
    { id: 'hi', name: 'hydroiodic acid (HI)', forms: ['I', '[H+].[I-]'], concMin: 45, diluteMax: 20, presets: [{ v: 57, unit: 'ww', note: 'conc. aqueous HI' }, { v: 10, unit: 'ww' }],
      why: 'Conc. HI cleaves ethers and reduces many groups.' },
    { id: 'hf', name: 'hydrofluoric acid (HF)', forms: ['F'], concMin: 40, diluteMax: 15, presets: [{ v: 48, unit: 'ww', note: 'conc. aqueous HF' }, { v: 70, unit: 'ww', note: 'HF-pyridine' }, { v: 5, unit: 'ww' }],
      why: 'Toxic and corrosive; 48 % aqueous HF is the usual reagent.' },
    { id: 'naoh', name: 'sodium hydroxide (NaOH)', forms: ['[Na+].[OH-]', 'O[Na]', '[Na]O'], concMin: 30, diluteMax: 15, concMinM: 8, diluteMaxM: 4,
      presets: [{ v: 50, unit: 'ww', note: 'conc. NaOH, about 19 M' }, { v: 30, unit: 'ww' }, { v: 10, unit: 'ww', note: 'about 2.8 M' }, { v: 4, unit: 'ww', note: 'about 1 M' }, { v: 1, unit: 'ww' }],
      why: 'Conc. NaOH drives Cannizzaro and eliminations; dilute aqueous NaOH saponifies and runs aldol additions.' },
    { id: 'koh', name: 'potassium hydroxide (KOH)', forms: ['[K+].[OH-]', 'O[K]', '[K]O'], concMin: 30, diluteMax: 15, concMinM: 8, diluteMaxM: 4,
      presets: [{ v: 50, unit: 'ww', note: 'conc. KOH' }, { v: 45, unit: 'ww', note: 'commercial solution' }, { v: 10, unit: 'ww' }, { v: 5.6, unit: 'ww', note: 'about 1 M' }],
      why: 'As NaOH; also the strong base in the Wolff-Kishner reduction.' },
    { id: 'nh3', name: 'ammonia (NH3)', forms: ['N'], concMin: 20, diluteMax: 10, concMinM: 10, diluteMaxM: 4,
      presets: [{ v: 28, unit: 'ww', note: 'conc. aqueous ammonia, about 15 M (25-30 %)' }, { v: 10, unit: 'ww' }, { v: 3.4, unit: 'ww', note: 'about 2 M' }],
      why: 'Conc. aqueous ammonia is the nucleophile for aminations; dilute ammonia is a weak base (Tollens reagent).' },
    { id: 'h2o2', name: 'hydrogen peroxide (H2O2)', forms: ['OO'], concMin: 25, diluteMax: 12, presets: [{ v: 30, unit: 'ww', note: 'perhydrol; 35 % is also sold' }, { v: 35, unit: 'ww' }, { v: 50, unit: 'ww' }, { v: 3, unit: 'ww', note: 'pharmacy strength' }],
      why: '30 % H2O2 is the usual oxidant (hydroboration work-up, epoxidations); 3 % is a mild antiseptic.' },
    { id: 'naocl', name: 'sodium hypochlorite (bleach)', forms: ['[Na+].[O-]Cl', '[O-]Cl.[Na+]'], concMin: 10, diluteMax: 3, presets: [{ v: 12.5, unit: 'wv', note: 'industrial' }, { v: 6, unit: 'wv', note: 'household bleach (5-6 %)' }, { v: 1, unit: 'wv' }],
      why: 'Oxidises alcohols and sulfides; the strength decays on storage.' },
    { id: 'ch2o', name: 'formaldehyde (formalin)', forms: ['C=O'], concMin: 30, diluteMax: 10, presets: [{ v: 37, unit: 'ww', note: 'formalin (stabilised with 10-15 % methanol)' }, { v: 10, unit: 'ww' }],
      why: 'Aqueous formaldehyde is almost entirely the hydrate; the aldehyde is released on reaction.' },
    { id: 'tfa', name: 'trifluoroacetic acid (TFA)', forms: ['O=C(O)C(F)(F)F'], concMin: 90, diluteMax: 20, presets: [{ v: 100, unit: 'vv', note: 'neat' }, { v: 50, unit: 'vv', note: 'in CH2Cl2 (Boc removal)' }, { v: 25, unit: 'vv' }, { v: 10, unit: 'vv' }],
      why: 'Neat or 25-50 % TFA in CH2Cl2 removes Boc and tert-butyl groups.' },
    { id: 'n2h4', name: 'hydrazine', forms: ['NN'], concMin: 60, diluteMax: 20, presets: [{ v: 64, unit: 'ww', note: 'hydrazine hydrate (N2H4 content)' }, { v: 35, unit: 'ww', note: 'aqueous' }],
      why: 'Hydrazine hydrate is the usual form (Wolff-Kishner, hydrazone formation).' },
    { id: 'kmno4', name: 'potassium permanganate (KMnO4)', forms: ['[K+].[O-][Mn](=O)(=O)=O', '[K+].O=[Mn](=O)(=O)[O-]'], concMin: 4, diluteMax: 2,
      presets: [{ v: 1, unit: 'wv', cls: 'dilute', note: 'cold, dilute: syn-diol' }, { v: 5, unit: 'wv', note: 'near saturation (about 6 % at 20 °C): hot, concentrated cleaves' }],
      why: 'Cold dilute KMnO4 gives the syn-diol; hot concentrated KMnO4 cleaves the C=C and oxidises side chains.' },
  ];

  const GENERIC = [{ v: 100, unit: 'ww', note: 'neat / undiluted' }, { v: 50, unit: 'ww' }, { v: 25, unit: 'ww' }, { v: 10, unit: 'ww' }, { v: 5, unit: 'ww' }, { v: 1, unit: 'ww' }];
  const EQUIV = [
    { label: '0.05 equiv (5 mol %)', text: '0.05 equiv' }, { label: '0.1 equiv (10 mol %, catalytic)', text: '0.1 equiv' },
    { label: '1 equiv', text: '1 equiv' }, { label: '1.1 equiv', text: '1.1 equiv' }, { label: '1.5 equiv', text: '1.5 equiv' },
    { label: '2 equiv', text: '2 equiv' }, { label: '3 equiv', text: '3 equiv' }, { label: '5 equiv', text: '5 equiv' },
    { label: 'excess (reacts at every group it can)', text: 'excess' },
  ];
  const UNIT_LABEL = { ww: '% w/w', vv: '% v/v', wv: '% w/w', M: 'M' };
  UNIT_LABEL.wv = '% w/v';
  const UNIT_CODE = { ww: '%w/w', vv: '%v/v', wv: '%w/v', M: 'M' };

  let index = null, indexCanon = null, partsOf = null;
  /* The reagent line is matched to the table by its parts: exactly (NaOH as [Na+].[OH-]), or, for a mixture such as Zn + HCl (Clemmensen) or HNO3 + H2SO4
     (nitration), by the first table entry whose parts are all in the line. The strength of a mixture is read against that first reagent's thresholds. */
  function entryFor(canon, line) {
    if (!line) return null;
    if (!index || indexCanon !== canon) {
      index = new Map(); indexCanon = canon; partsOf = [];
      for (const e of TABLE) for (const f of e.forms) { index.set(key(canon, f), e); partsOf.push([e, new Set(key(canon, f).split('.'))]); }
    }
    const k = key(canon, line);
    if (index.has(k)) return index.get(k);
    const have = new Set(k.split('.'));
    for (const [e, ps] of partsOf) if ([...ps].every(x => have.has(x))) return e;
    return null;
  }
  function key(canon, line) {
    return [...new Set(String(line).split('.').filter(Boolean).map(p => { try { return canon(p) || p; } catch (e) { return p; } }))].sort().join('.');
  }

  const fmtNum = v => (Math.round(v * 100) / 100).toString();
  function format(c) {
    if (!c) return '';
    if (c.v == null) return c.cls === 'concentrated' ? 'concentrated' : c.cls === 'dilute' ? 'dilute' : '';
    return fmtNum(c.v) + (c.unit === 'M' ? ' M' : ' ' + (UNIT_LABEL[c.unit] || '% w/w'));
  }
  function short(c) {
    if (!c) return '';
    if (c.v == null) return c.cls === 'concentrated' ? 'conc.' : c.cls === 'dilute' ? 'dilute' : '';
    return fmtNum(c.v) + (c.unit === 'M' ? ' M' : '%');
  }
  function code(c) {
    if (!c) return '';
    if (c.v == null) return c.cls === 'concentrated' ? 'conc' : c.cls === 'dilute' ? 'dil' : '';
    return fmtNum(c.v) + (UNIT_CODE[c.unit] || '%w/w');
  }
  function parse(text) {
    const s = String(text || '').trim().toLowerCase().replace(/\s+/g, '');
    if (!s) return null;
    if (/^(conc|conc\.|concentrated)$/.test(s)) return { cls: 'concentrated' };
    if (/^(dil|dil\.|dilute|diluted)$/.test(s)) return { cls: 'dilute' };
    let m = s.match(/^(\d+(?:\.\d+)?)m(?:ol\/l|ol)?$/);
    if (m && +m[1] > 0) return { v: +m[1], unit: 'M' };
    m = s.match(/^(\d+(?:\.\d+)?)%(?:\(?(w\/w|v\/v|w\/v)\)?)?$/);
    if (m && +m[1] > 0 && +m[1] <= 100) return { v: +m[1], unit: m[2] === 'v/v' ? 'vv' : m[2] === 'w/v' ? 'wv' : 'ww' };
    m = s.match(/^(\d+(?:\.\d+)?)(w\/w|v\/v|w\/v)%$/);
    if (m && +m[1] > 0 && +m[1] <= 100) return { v: +m[1], unit: m[2] === 'v/v' ? 'vv' : m[2] === 'w/v' ? 'wv' : 'ww' };
    return null;
  }

  function classify(entry, c) {
    if (!c) return null;
    if (c.v == null) return c.cls || null;
    const e = entry || {};
    if (c.unit === 'M') {
      const hi = e.concMinM != null ? e.concMinM : 10, lo = e.diluteMaxM != null ? e.diluteMaxM : 2;
      return c.v >= hi ? 'concentrated' : c.v <= lo ? 'dilute' : 'moderate';
    }
    const hi = e.concMin != null ? e.concMin : 70, lo = e.diluteMax != null ? e.diluteMax : 15;
    return c.v >= hi ? 'concentrated' : c.v <= lo ? 'dilute' : 'moderate';
  }

  function presets(entry) {
    const out = [];
    const mk = (p, own) => {
      const c = { v: p.v, unit: p.unit };
      const cls = p.cls || classify(own ? entry : null, c);
      out.push({ label: format(c), v: p.v, unit: p.unit, cls, note: p.note || '', own: !!own });
    };
    if (entry) entry.presets.forEach(p => mk(p, true));
    GENERIC.filter(g => !entry || !entry.presets.some(p => p.v === g.v && p.unit === g.unit)).forEach(p => mk(p, false));
    return out;
  }

  const api = { TABLE, GENERIC, EQUIV, entryFor, presets, classify, format, short, code, parse, UNIT_LABEL };
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.ReagentConc = api;
})(typeof window !== 'undefined' ? window : globalThis);
