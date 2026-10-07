"""pKa of each fragment IN ITS MOLECULE: measured if the molecule has a measured value, else estimated from the
substituents (Hammett equation), else the fragment's typical value (Fragments sheet of molecule_data.xlsx).

    python tools/pka.py            # builds tools/data/pka_measured.json + data/pka.js and prints how well the estimates do

Measured values: IUPAC Digitized pKa Dataset v2.4 (tools/data/iupac_pka/), CC BY-NC 4.0: "Reproduced by permission
of International Union of Pure and Applied Chemistry" (Zheng & Lafontant-Joseph, doi:10.5281/zenodo.7236452).
Only water (no co-solvent) at 20-30 °C; the most reliable assessment available (Reliable > Approximate > not stated >
Uncertain); the median of those values. "pKa1, pKa2 ..." are the molecule's acidic sites; "pKaH1 ..." the pKa of the
protonated form of its basic sites.

Estimates (Hammett equation, pKa = pKa0 - rho * sum of sigma, for meta and para substituents on a benzene ring;
sigma, sigma- and rho from C. Hansch, A. Leo, R. W. Taft, Chem. Rev. 1991, 91, 165 and standard tables):
  phenol (F_EBFEFTNC):            pKa0 9.99, rho 2.23, para resonance-withdrawing groups use sigma- (4-nitrophenol: 9.99 - 2.23*1.27 = 7.16)
  benzoic acid (F_WH4T6PP5 on an aryl ring): pKa0 4.20, rho 1.00, sigma
  anilinium (F_UY2ANYI2, pKaH):   pKa0 4.60, rho 2.89, para resonance-withdrawing groups use sigma-
An ortho substituent, a fused ring (any ring atom shared with another ring), another heteroaromatic ring, or an unknown substituent: no estimate (the typical value is
used and the page says why). The page's fragmentPka (js/reactions.js) computes the same.
"""
import collections, csv, json, os, re, statistics, sys

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
SRC = os.path.join(HERE, 'data', 'iupac_pka', 'iupac_high-confidence_v2_4.csv')
CITE = ('IUPAC Digitized pKa Dataset v2.4 (Zheng & Lafontant-Joseph, doi:10.5281/zenodo.7236452), CC BY-NC 4.0: '
        'Reproduced by permission of International Union of Pure and Applied Chemistry')

# substituent: SMARTS from the ring carbon's neighbour outwards ([*:1] = the atom bonded to the ring), sigma_m, sigma_p, sigma_p-
SIGMA = [
    ('nitro',       '[N+:1](=O)[O-]',          0.71, 0.78, 1.27),
    ('cyano',       '[C:1]#N',                 0.56, 0.66, 1.00),
    ('trifluoromethyl', '[C:1](F)(F)F',        0.43, 0.54, 0.65),
    ('methylsulfonyl', '[S:1](=O)(=O)[#6]',    0.60, 0.72, 1.13),
    ('sulfonamide', '[S:1](=O)(=O)[NX3]',      0.53, 0.60, 0.94),
    ('formyl',      '[CH1:1]=O',               0.35, 0.42, 1.03),
    ('acetyl (ketone)', '[C:1](=O)[#6]',       0.38, 0.50, 0.84),
    ('carboxylic acid', '[C:1](=O)[OH]',       0.37, 0.45, 0.77),
    ('ester',       '[C:1](=O)O[#6]',          0.37, 0.45, 0.75),
    ('amide',       '[C:1](=O)[NX3]',          0.28, 0.36, 0.61),
    ('trifluoromethoxy', '[OX2;+0:1]C(F)(F)F',      0.38, 0.35, 0.35),
    ('fluoro',      '[F:1]',                   0.34, 0.06, -0.03),
    ('chloro',      '[Cl:1]',                  0.37, 0.23, 0.19),
    ('bromo',       '[Br:1]',                  0.39, 0.23, 0.25),
    ('iodo',        '[I:1]',                   0.35, 0.18, 0.27),
    ('trifluoromethylthio', '[SX2;+0:1]C(F)(F)F', 0.40, 0.50, 0.57),
    ('methylthio',  '[SX2;+0:1][CX4;!$(C~[!#6;!#1])]',              0.15, 0.00, 0.06),
    ('acetamido',   '[NX3;H1;+0:1]C(=O)',             0.21, 0.00, -0.46),
    ('dimethylamino', '[NX3;+0:1]([CX4])[CX4]',   -0.16, -0.83, -0.12),
    ('amino',       '[NX3;H2;+0:1]',                 -0.16, -0.66, -0.15),
    ('hydroxy',     '[OX2;H1;+0:1]',                  0.12, -0.37, -0.37),
    ('phenoxy',     '[OX2;+0:1]c',                  0.25, -0.03, -0.10),
    ('alkoxy',      '[OX2;+0:1][CX4]',              0.12, -0.27, -0.26),
    ('tert-butyl',  '[CX4;!$(C~[!#6;!#1]):1](C)(C)C',          -0.10, -0.20, -0.13),
    ('methyl',      '[CH3;!$(C~[!#6;!#1]):1]',                 -0.07, -0.17, -0.17),
    ('alkyl',       '[CX4;!$(C~[!#6;!#1]):1]',                 -0.07, -0.15, -0.15),
    ('phenyl',      '[c:1]',                   0.06, -0.01, 0.02),
    ('vinyl',       '[CX3:1]=[CX3]',           0.06, -0.04, -0.08),
]
# the acid/base sites estimated: fragment, which pKa, SMARTS of the site ([*:1] = the ring atom it sits on / the ring N)
HAMMETT = {
    'F_EBFEFTNC': dict(kind='acid', smarts='[c:1][OX2H1]', pKa0=9.99, rho=2.23, para='sigma-', name='phenol'),
    'F_WH4T6PP5': dict(kind='acid', smarts='[c:1][CX3](=O)[OX2H1]', pKa0=4.20, rho=1.00, para='sigma', name='benzoic acid'),
    'F_UY2ANYI2': dict(kind='base', smarts='[c:1][NX3;H2;+0]', pKa0=4.60, rho=2.89, para='sigma-', name='anilinium'),
}
# (pyridinium, pKa0 5.25 rho 5.90, was tested and dropped: within 0.5 of the measured value for only about half of 63
#  measured pyridines; pyridines use measured values or the fragment's typical value)
_Q = {}


def _q(sm):
    from rdkit import Chem
    if sm not in _Q: _Q[sm] = Chem.MolFromSmarts(sm)
    return _Q[sm]


def canon(smi):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m) if m else None


# ---------- measured ----------
RANK = {'Reliable': 0, 'Approximate': 1, 'Probably approximate': 1, '': 2, 'Unknown': 2, 'Uncertain': 3, 'Very uncertain': 4}


def measured_table():
    """{canonical SMILES: {'acid': [pKa, ...] most acidic first, 'base': [pKaH, ...] most basic first, 'n': measurements}}"""
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    best = collections.defaultdict(list)                 # (smiles, type) -> [(rank, value)]
    for r in csv.DictReader(open(SRC, encoding='utf-8')):
        if r['cosolvent'] or r['pka_type'] not in {f'pKa{i}' for i in range(1, 9)} | {f'pKaH{i}' for i in range(1, 9)}: continue
        try: t = float(r['T'])
        except ValueError: t = 25.0 if r['T'] in ('', 'not_stated') else None
        if t is None or not 20 <= t <= 30: continue
        try: v = float(r['pka_value'])
        except ValueError: continue
        best[(r['SMILES'], r['pka_type'])].append((RANK.get(r['assessment'], 3), v))
    out = {}
    cache = {}
    for (smi, typ), vals in best.items():
        if smi not in cache: cache[smi] = canon(smi)
        c = cache[smi]
        if not c: continue
        top = min(rk for rk, _ in vals)
        v = round(statistics.median(x for rk, x in vals if rk == top), 2)
        d = out.setdefault(c, {'acid': {}, 'base': {}, 'quality': {}})
        d['acid' if typ.startswith('pKa') and not typ.startswith('pKaH') else 'base'][typ] = v
        d['quality'][typ] = ['Reliable', 'Approximate', 'not stated', 'Uncertain', 'Very uncertain'][top]
    table = {}
    for c, d in out.items():
        acid = [d['acid'][k] for k in sorted(d['acid'], key=lambda k: int(k[3:]))]
        base = sorted(d['base'].values(), reverse=True)  # the most basic site first
        q = sorted(set(d['quality'].values()), key=lambda x: ['Reliable', 'Approximate', 'not stated', 'Uncertain', 'Very uncertain'].index(x))
        table[c] = {'acid': acid, 'base': base, 'quality': q[-1] if q else ''}
    return table


# ---------- estimated (Hammett) ----------
def hammett(m, fid):
    """[(estimate, explanation)] for each occurrence of the fragment's site in molecule m; None where it cannot be estimated"""
    h = HAMMETT.get(fid)
    if not h: return []
    out = []
    for match in m.GetSubstructMatches(_q(h['smarts'])):
        site = match[0]                                  # the ring atom carrying the group (or the ring N)
        a = m.GetAtomWithIdx(site)
        ri = m.GetRingInfo()
        rings = [r for r in ri.AtomRings() if site in r and len(r) == 6]
        if len(rings) != 1 or any(ri.NumAtomRings(i) != 1 for i in rings[0]):
            out.append((None, 'fused ring: no estimate')); continue
        ring = rings[0]
        if any(not m.GetAtomWithIdx(i).GetIsAromatic() for i in ring) or (fid != 'F_WR22PQX5' and any(m.GetAtomWithIdx(i).GetSymbol() != 'C' for i in ring)) \
                or (fid == 'F_WR22PQX5' and sum(1 for i in ring if m.GetAtomWithIdx(i).GetSymbol() != 'C') != 1):
            out.append((None, 'not a benzene (or simple pyridine) ring: no estimate')); continue
        pos = {i: min(abs(ring.index(i) - ring.index(site)), 6 - abs(ring.index(i) - ring.index(site))) for i in ring}
        own = set(match)
        total, parts, ok = 0.0, [], True
        for i in ring:
            if i == site: continue
            for nb in m.GetAtomWithIdx(i).GetNeighbors():
                j = nb.GetIdx()
                if j in ring or j in own: continue
                if pos[i] == 1: ok = False; parts.append('ortho substituent'); break
                kind = None
                for name, sm, sm_m, sm_p, sm_pm in SIGMA:
                    q = _q(sm)
                    if any(mt[[x.GetAtomMapNum() for x in q.GetAtoms()].index(1)] == j for mt in m.GetSubstructMatches(q)):
                        # sigma- only for a para group that withdraws by resonance (sigma- > sigma: NO2, CN, C=O, SO2):
                        # through-conjugation with the O-/NH2 site; donors (OMe, NHAc, NH2) take their ordinary sigma_p
                        para_minus = h['para'] == 'sigma-' and sm_pm > sm_p
                        kind = (name, sm_m if pos[i] == 2 else (sm_pm if para_minus else sm_p), pos[i] == 3 and para_minus); break
                if kind is None: ok = False; parts.append('a substituent with no sigma value'); break
                total += kind[1]
                parts.append(f'{"meta" if pos[i] == 2 else "para"}-{kind[0]} ({"σm" if pos[i] == 2 else ("σp⁻" if kind[2] else "σp")} {kind[1]:+.2f})')
            if not ok: break
        if not ok:
            out.append((None, f'{parts[-1]}: no estimate')); continue
        v = round(h['pKa0'] - h['rho'] * total, 2)
        out.append((v, f'Hammett: {h["name"]} {h["pKa0"]} − {h["rho"]} × ({" + ".join(parts) if parts else "no substituent"}) = {v}'))
    return out


def fragment_defaults():
    """{fragment ID: (pKa, pKaH)} typical values from the Fragments sheet of molecule_data.xlsx (first number of the cell)"""
    sys.path.insert(0, HERE)
    import learn_reactions as L
    num = lambda v: (lambda m: float(m.group(0).replace('−', '-')) if m else None)(re.search(r'[-−]?\d+(?:\.\d+)?', str(v or '')))
    return {str(f['Fragment ID']).strip(): (num(f.get('pKa (its most acidic H)')), num(f.get('pKa of its conjugate acid (how basic)')))
            for f in L.read_sheet(os.path.join(APP, 'molecule_data.xlsx'), 'Fragments') if str(f.get('Fragment ID', '')).startswith('F')}


def molecule_pkas(smiles, frags, defaults, measured):
    """Every acidic and basic SITE of a molecule with its pKa and where the number comes from.
    A site is one atom: the H-bearing O/N/S (else C) of an acidic fragment with a typical pKa <= 20, or the N of a
    basic fragment with a typical pKaH >= 0; several fragments matching the same atom are one site (typical value:
    the most acidic / most basic of them). Each site starts with its Hammett estimate or typical value; the molecule's
    measured values (IUPAC) then replace them in order: the most acidic measured value goes to the site expected to
    be most acidic, and so on (basic sites: the most basic first). Returns
    {'acid': [{atom, fragments, pKa, source}], 'base': [...]}."""
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    out = {'acid': {}, 'base': {}}
    for fid, q in frags:
        if q is None: continue
        dk, dh = defaults.get(fid, (None, None))
        for match in m.GetSubstructMatches(q, uniquify=True, maxMatches=1000):
            if dk is not None and dk <= 20:
                het = [i for i in match if m.GetAtomWithIdx(i).GetSymbol() in ('O', 'N', 'S') and m.GetAtomWithIdx(i).GetTotalNumHs() > 0]
                # a carbon acid's acidic H is on an sp3 carbon (the alpha C of an aldehyde, not its formyl C-H)
                car = [i for i in match if m.GetAtomWithIdx(i).GetSymbol() == 'C' and m.GetAtomWithIdx(i).GetTotalNumHs() > 0
                       and not m.GetAtomWithIdx(i).GetIsAromatic() and all(b.GetBondTypeAsDouble() == 1.0 for b in m.GetAtomWithIdx(i).GetBonds())]
                site = (het or car or [None])[0]
                if site is not None:
                    s_ = out['acid'].setdefault(site, {'atom': site, 'fragments': [], 'typical': dk})
                    s_['fragments'].append(fid); s_['typical'] = min(s_['typical'], dk)
            if dh is not None and dh >= 0:
                ns = [i for i in match if m.GetAtomWithIdx(i).GetSymbol() == 'N' and m.GetAtomWithIdx(i).GetFormalCharge() == 0]
                if ns:
                    s_ = out['base'].setdefault(ns[0], {'atom': ns[0], 'fragments': [], 'typical': dh})
                    s_['fragments'].append(fid); s_['typical'] = max(s_['typical'], dh)
    cls = symmetry_classes(m)
    for kind in ('acid', 'base'):                        # symmetry-equivalent sites (acetone's two CH3) are one site
        merged = {}
        for a in sorted(out[kind]):
            s_ = out[kind][a]
            key = cls[a]
            if key in merged:
                g = merged[key]; g['atoms'].append(a); g['fragments'] += s_['fragments']
                g['typical'] = min(g['typical'], s_['typical']) if kind == 'acid' else max(g['typical'], s_['typical'])
            else:
                merged[key] = dict(s_, atoms=[a])
        out[kind] = {g['atom']: g for g in merged.values()}
    for kind in ('acid', 'base'):
        for s_ in out[kind].values():
            s_['pKa'], s_['source'] = s_['typical'], 'typical value for ' + ', '.join(sorted(set(s_['fragments'])))
            for fid in s_['fragments']:
                if HAMMETT.get(fid, {}).get('kind') == kind:
                    for v, src in [e for e in hammett(m, fid)]:
                        if v is not None and _site_of(m, fid) in s_['atoms']:
                            s_['pKa'], s_['source'] = v, 'estimated: ' + src
                        elif v is None and _site_of(m, fid) in s_['atoms']:
                            s_['source'] += f' ({src})'
        out[kind] = sorted(out[kind].values(), key=lambda x: x['pKa'], reverse=(kind == 'base'))
    # measured values -> sites: the assignment closest to the expected values. Acid pKa and the pKaH of a basic site are
    # on one scale (a proton leaving), so they are matched together: an amino acid (recorded as its zwitterion, where
    # the dataset's 'acid' value is the NH3+) still gets 10.2 on the amine and 3.6 on the acid.
    # an aldehyde or ketone measured in water is partly hydrated (R2C=O + H2O <-> R2C(OH)2): with no O-H/N-H/S-H of its
    # own, a measured acid value below 15 is the hydrate's O-H (acetaldehyde 13.57, formaldehyde 13.3); it goes on the O
    meas_acid = list((measured or {}).get('acid') or [])
    hydrate = None
    co = m.GetSubstructMatches(_q('[CX3;$([CH1](=O)[#6]),$([CH2]=O),$(C(=O)([#6])[#6])]=[OX1]'))
    if co and meas_acid and min(meas_acid) < 15 and not any(m.GetAtomWithIdx(x['atom']).GetSymbol() in ('O', 'N', 'S') for x in out['acid']):
        v = min(meas_acid); meas_acid.remove(v)
        fid = 'F_ADOAMFFM' if m.GetAtomWithIdx(co[0][0]).GetTotalNumHs() > 0 else 'F_CFJ46WPI'
        hydrate = {'atom': co[0][1], 'atoms': [co[0][1]], 'fragments': [fid], 'typical': v, 'pKa': v,
                   'source': 'measured (IUPAC Digitized pKa Dataset) in water, where the C=O is partly hydrated: the O-H of the hydrate R2C(OH)2'}
    groups = out['acid'] + out['base']
    # one slot per atom: a group of equivalent sites (succinic acid's two COOH) takes as many measured values as it has
    # atoms (4.21 for the first proton, 5.67 for the second)
    slots = [g for g in groups for _ in g['atoms']]
    # a measured pKaH below 0 is the protonation of a very weak base (C=O, ether: acetone -7.2), not a basic site here
    meas = meas_acid + [v for v in ((measured or {}).get('base') or []) if v >= 0]
    if meas and len(meas) <= len(slots):
        import itertools
        best = None
        if len(slots) <= 8:
            for perm in itertools.permutations(range(len(slots)), len(meas)):
                cost = sum(abs(slots[j]['pKa'] - v) for j, v in zip(perm, meas))
                if best is None or cost < best[0]: best = (cost, perm)
            pairs = list(zip(best[1], meas))
        else:                                            # many sites: most acidic measured value to the most acidic site, ...
            order = sorted(range(len(slots)), key=lambda j: slots[j]['pKa'])
            pairs = list(zip(order, sorted(meas)))
        got = {}
        for j, v in pairs: got.setdefault(id(slots[j]), (slots[j], []))[1].append(v)
        for g, vals in got.values():
            vals.sort(reverse=g in out['base'])
            before = ('the estimate was ' if g['source'].startswith('estimated') else 'the typical value was ') + str(g['pKa'])
            g['source'] = 'measured (IUPAC Digitized pKa Dataset)' + (f'; {before}' if abs(vals[0] - g['pKa']) > 0.005 else '') + \
                ''.join(f'; then {v} for the next equivalent site' for v in vals[1:])
            g['pKa'] = vals[0]
    elif meas:
        for g in groups: g['source'] += f'; the molecule has {len(meas)} measured values for {len(slots)} site(s), not assigned'
    if hydrate: out['acid'].append(hydrate)
    return out


def symmetry_classes(m):
    """atoms that are equivalent by symmetry get the same class: start from element, heavy degree, H count, charge and
    aromaticity, then refine with the sorted classes of the neighbours until nothing changes (the page's symmetryClasses
    does the same, so both group the same sites)"""
    inv = [f'{a.GetSymbol()},{a.GetDegree()},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.GetIsAromatic())}' for a in m.GetAtoms()]
    cls = _renumber(inv)
    for _ in range(m.GetNumAtoms()):
        new = _renumber([f'{cls[i]}|' + ','.join(str(x) for x in sorted(cls[n.GetIdx()] for n in m.GetAtomWithIdx(i).GetNeighbors())) for i in range(m.GetNumAtoms())])
        if len(set(new)) == len(set(cls)): break
        cls = new
    return cls


def _renumber(labels):
    order = {lab: k for k, lab in enumerate(sorted(set(labels)))}
    return [order[lab] for lab in labels]


def _site_of(m, fid):
    """the atom a Hammett-estimated fragment's site is on (the OH oxygen, the COOH OH oxygen, the NH2 nitrogen)"""
    h = HAMMETT[fid]
    for match in m.GetSubstructMatches(_q(h['smarts'])):
        for i in match[1:]:
            a = m.GetAtomWithIdx(i)
            if (h['kind'] == 'acid' and a.GetSymbol() == 'O' and a.GetTotalNumHs() > 0) or (h['kind'] == 'base' and a.GetSymbol() == 'N'):
                return i
    return None


def main():
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
    table = measured_table()
    json.dump({'source': CITE, 'molecules': table}, open(os.path.join(HERE, 'data', 'pka_measured.json'), 'w', encoding='utf-8'), ensure_ascii=False)
    with open(os.path.join(APP, 'data', 'pka.js'), 'w', encoding='utf-8') as f:
        f.write('/* generated by tools/pka.py: measured pKa values. ' + CITE + ' */\n')
        f.write('window.PKA_MEASURED = ' + json.dumps({'source': CITE, 'molecules': table}, ensure_ascii=False, separators=(',', ':')) + ';\n')
    print(f'measured pKa for {len(table):,} molecules -> tools/data/pka_measured.json, data/pka.js')
    # how good are the estimates? molecules with exactly one site of the kind and one measured value of that kind
    for fid, h in HAMMETT.items():
        errs, worst = [], []
        for smi, d in table.items():
            vals = d['acid'] if h['kind'] == 'acid' else d['base']
            if len(vals) != 1: continue
            m = Chem.MolFromSmiles(smi)
            est = hammett(m, fid)
            if len(est) != 1 or est[0][0] is None: continue
            # only molecules with exactly one acidic (or basic) site, so the measured value belongs to this one
            sites = {'acid': '[$([OX2H1]),$([NX3;H1,H2]S(=O)=O),$([SX2H1]),$([nH]),$([NX4+]),$([NX3;+0;!$(N-C=O);!$(N-S(=O)=O);!$(N-a)])]',
                     'base': '[$([NX3;+0;!$(N-C=O);!$(N-S(=O)=O)]),$([nX2])]'}[h['kind']]
            if len(m.GetSubstructMatches(_q(sites))) != 1: continue
            if h['kind'] == 'base' and m.HasSubstructMatch(_q('[$([OX2H1]),$([NX3;H1,H2]S(=O)=O),$([SX2H1])]')): continue
            e = est[0][0] - vals[0]
            errs.append(abs(e)); worst.append((abs(e), smi, vals[0], est[0][0]))
        if errs:
            worst.sort(reverse=True)
            print(f'{h["name"]:13s}: {len(errs):4d} measured molecules; estimate within 0.5 of measured for {sum(1 for x in errs if x <= 0.5) / len(errs):.0%}, '
                  f'median error {statistics.median(errs):.2f}; worst: ' + '; '.join(f'{s} meas {mv} est {ev}' for _, s, mv, ev in worst[:3]))


if __name__ == '__main__':
    main()
