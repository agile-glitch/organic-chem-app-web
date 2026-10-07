"""Reactivity properties of every fragment, computed (never typed in) on an automatically chosen model compound.

    python tools/fragment_properties.py                 # all fragments -> tools/data/fragment_properties.json
    python tools/fragment_properties.py --only alcohol  # only fragments whose name contains this (testing)
    python tools/fragment_properties.py --workers 3     # fragments in parallel (one DFT thread each)
    python tools/fragment_properties.py --resume        # continue a stopped run / compute only fragments that have no finished row
    python tools/fragment_properties.py --refresh       # recompute only the table-based columns, keep the DFT values

The model compound of a fragment is picked by rule, not by hand, from the app's own molecule lists (the Molecules sheet of
molecule_data.xlsx, tools/data/reagents.json and the small molecules of the recorded reactions in data/learned/rx): among the
molecules that contain the fragment and have at least one heavy atom beyond its pattern, the one with, in this order, no charged
atom outside the pattern (unless the pattern is an ion), no atom of a 3-4-membered ring outside it, no ring inside the match
unless the pattern is a ring pattern, the least complexity outside the pattern (+1 per ring / aromatic / multiply bonded atom,
+2 per heteroatom), the fewest H on carbonyl-type pattern carbons (so an ester is an acetate, not a formate), the fewest heavy atoms, a single occurrence of the fragment, and the fewest other fragments (arene
C-H ignored). Recorded-reaction molecules count only if they occur at least 3 times (compounds the data really uses). The choice is cached in tools/data/fragment_models.json. Ionic salts: only the component that holds
the fragment. A property of the model compound is a property of the fragment in its simplest setting; the real value in a
real molecule moves with its substituents (tools/effects.py reads that, molecule by molecule).

Geometry: RDKit conformers, then GFN2-xTB (tblite) optimisation of the lowest few, lowest energy kept.
Electronic structure: B3LYP/def2-SVP (def2 ECP for Sn, I), ddCOSMO water continuum (lmax 4, Lebedev 11) so that anions behave, PySCF.
  * HOMO / LUMO energies are Kohn-Sham orbital energies: use them to RANK fragments, not as ionisation energies.
  * atomic shares of an orbital: Lowdin populations (orthogonalised AOs); they stand in for the Fukui functions
    f- (HOMO, where an electrophile attacks) and f+ (LUMO, where a nucleophile attacks), the frontier-orbital approximation.
  * hardness eta = (E_LUMO - E_HOMO)/2, softness S = 1/eta, local softness s = S * share, electrophilicity
    omega = mu^2 / (2 eta) with mu = (E_HOMO + E_LUMO)/2 (Parr, Pearson).
  * polarisability: static, isotropic, by finite electric field (forward difference, 0.002 a.u.) in the same continuum.
Tabulated values come from the app's own tables so the page and these numbers agree: electronegativity (Allred-Pauling),
average bond enthalpies (OpenStax) and typical bond lengths (Allen 1987) from js/mol_data.js, Hammett sigma from tools/pka.py
(Hansch, Leo, Taft 1991), steric / Taft values and the 1-5 steric class from tools/effects.py, pKa from the Fragments sheet.
"""
import argparse, collections, glob, json, math, os, re, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, 'data', 'fragment_properties.json')
ALLOWED = {'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I', 'Mg', 'Al', 'Sn', 'Li', 'Na', 'K'}
MAX_HEAVY = 20
EPS_WATER = 78.3553
BOHR = 0.529177210903
HARTREE_EV = 27.211386245988

# fragment (by name) -> the substituent of tools/pka.py SIGMA that the fragment is (as a ring substituent)
SIGMA_OF = {'nitroarene': 'nitro', 'nitrile': 'cyano', 'aldehyde': 'formyl', 'ketone': 'acetyl (ketone)',
            'carboxylic acid': 'carboxylic acid', 'ester': 'ester', 'amide': 'amide',
            'sulfonamide, primary (R–SO2NH2)': 'sulfonamide', 'aryl fluoride, SNAr-activated (o/p-EWG or ring N)': 'fluoro',
            'aryl chloride': 'chloro', 'aryl bromide': 'bromo', 'aryl iodide (Ar–I, incl. heteroaryl iodides)': 'iodo',
            'phenol': 'hydroxy', 'arylamine (aniline NH2)': 'amino', 'amine 3° (trialkylamine, N,N-dialkylaniline, enamine)': 'dimethylamino',
            'ether (dialkyl)': 'alkoxy', 'alkene (C=C)': 'vinyl', 'alkyl arene (aryl–C sp3)': 'alkyl', 'aryl methyl (benzylic C–H)': 'methyl'}


# ---------------------------------------------------------------- the app's own tables
def _mol_data():
    s = open(os.path.join(APP, 'js', 'mol_data.js'), encoding='utf-8').read()
    a = s.index('/*JSON-BEGIN*/') + len('/*JSON-BEGIN*/'); b = s.index('/*JSON-END*/')
    return json.loads(s[a:b])


_D = _mol_data()
EL = {e['sym']: e for e in _D['ELEMENTS']}
BE = {}
for _r in _D['BOND_ENTHALPY']:
    BE[(_r['a'], _r['b'], _r['order'])] = BE[(_r['b'], _r['a'], _r['order'])] = _r['kJ']


def en(sym):
    return (EL.get(sym) or {}).get('en')


def bond_enthalpy(a, b, order):
    return None if order == 1.5 else BE.get((a, b, order))


def bond_length(a, b, order, h1=None, h2=None):
    """typical length (Angstrom) from the Allen 1987 rows of mol_data.js: element pair, order, no env row; a hybridisation that
    contradicts the row is excluded, a matching one wins, then the row marked default"""
    best, score = None, -1
    for r in _D['BOND_LENGTH']:
        if r.get('env') or r['order'] != order: continue
        if {r['a'], r['b']} != {a, b}: continue
        hy = r.get('hyb') or [None, None]
        ends = [(r['a'], hy[0]), (r['b'], hy[1])]
        if r['a'] == a: want = [(a, h1), (b, h2)]
        else: want = [(b, h2), (a, h1)]
        if any(e[1] and w[1] and e[1] != w[1] for e, w in zip(ends, want)): continue
        sc = sum(1 for e, w in zip(ends, want) if e[1] and e[1] == w[1]) * 2 + (1 if r.get('def') else 0)
        if sc > score: best, score = r['A'], sc
    return best


# ---------------------------------------------------------------- fragments and model compounds
def sheet(path, name):
    import openpyxl
    ws = openpyxl.load_workbook(os.path.join(APP, path), read_only=True, data_only=True)[name]
    it = ws.iter_rows(values_only=True); head = [str(h or '').strip() for h in next(it)]
    return [dict(zip(head, r)) for r in it if r and r[0] not in (None, '')]


def fragments():
    return [f for f in sheet('molecule_data.xlsx', 'Fragments') if str(f['Fragment ID']).startswith('F')]


def base_pool():
    pool = []
    for r in sheet('molecule_data.xlsx', 'Molecules'):
        if r.get('SMILES'): pool.append((f"Molecules sheet: {r.get('Common name') or r['SMILES']}", str(r['SMILES'])))
    for r in json.load(open(os.path.join(HERE, 'data', 'reagents.json'), encoding='utf-8')):
        if r.get('smiles'): pool.append((f"reagents.json: {r.get('name') or r['smiles']}", r['smiles']))
    return pool


def recorded_pool(maxlen=26, minseen=3):
    """small molecules of the recorded reactions (starting materials and products of data/learned/rx) that occur at least
    minseen times: a model compound should be a compound the data really uses. The counts are cached in the temp dir."""
    import base64, gzip, tempfile
    cache = os.path.join(tempfile.gettempdir(), 'organic_chem_recorded_small_molecules.json')
    if os.path.exists(cache): count = json.load(open(cache, encoding='utf-8'))
    else:
        count = collections.Counter()
        for f in sorted(glob.glob(os.path.join(APP, 'data', 'learned', 'rx', '*.js'))):
            s = open(f, encoding='utf-8').read(); mm = re.search(r'"gz":"([^"]+)"', s)
            for row in json.loads(gzip.decompress(base64.b64decode(mm.group(1)))):
                parts = row[1].split('~')
                for smi in set(parts[0].split('.') + parts[-1].split('.')):
                    if smi and len(smi) <= maxlen: count[smi] += 1
        json.dump(count, open(cache, 'w', encoding='utf-8'))
    return [('recorded reactions (data/learned/rx)', smi) for smi, n in sorted(count.items()) if n >= minseen]


def components(smiles):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    if m is None: return []
    return [c for c in Chem.GetMolFrags(m, asMols=True, sanitizeFrags=True)]


MODELS = os.path.join(HERE, 'data', 'fragment_models.json')


def _ionic(smarts):
    """does the pattern itself ask for a charged atom ([O-], [N+], [AlH4-] ...)? A bond '-' is followed by an atom, a charge by ] ; , & :"""
    return bool(re.search(r'[+-](?:[1-9])?(?=[\];,&:)])', smarts))


def choose_models(frags, log=print, rechoose=False):
    """{fragment id: [source label, SMILES]} by the rule in the module docstring. The choice is cached in
    tools/data/fragment_models.json, keyed by fragment ID: the ID is a hash of the SMARTS, so an edited pattern gets
    a new ID and is chosen afresh; --rechoose forgets the cache."""
    from rdkit import Chem
    cache = {} if rechoose or not os.path.exists(MODELS) else json.load(open(MODELS, encoding='utf-8'))
    queries = {f['Fragment ID']: Chem.MolFromSmarts(f['SMARTS pattern']) for f in frags}
    name = {f['Fragment ID']: str(f['Fragment name']) for f in frags}
    ionic = {f['Fragment ID']: _ionic(str(f['SMARTS pattern'])) for f in frags}
    ignore = {i for i, n in name.items() if n.startswith('arene C–H')}
    ring_pattern = {}
    for i, q in queries.items():
        Chem.FastFindRings(q)
        ring_pattern[i] = q.GetRingInfo().NumRings() > 0 or bool(re.search(r'r\d', str(next(f['SMARTS pattern'] for f in frags if f['Fragment ID'] == i))))
    todo = {i for i in queries if i not in cache}
    if not todo: return {i: cache[i] for i in queries}
    log(f'choosing model compounds for {len(todo)} fragment(s) ...')

    def search(pool, fids):
        cand = collections.defaultdict(list)
        for label, smi in pool:
            for c in components(smi):
                if c.GetNumHeavyAtoms() > MAX_HEAVY or not {a.GetSymbol() for a in c.GetAtoms()} <= ALLOWED: continue
                hits = {i: c.GetSubstructMatch(queries[i]) for i in fids if c.HasSubstructMatch(queries[i])}
                if not hits: continue
                csmi = Chem.MolToSmiles(c); every = [i for i, q in queries.items() if c.HasSubstructMatch(q)]
                small = {k for r in c.GetRingInfo().AtomRings() if len(r) <= 4 for k in r}
                for i, mt in hits.items():
                    match = set(mt); atoms = c.GetAtoms()
                    charged = 0 if ionic[i] else sum(1 for a in atoms if a.GetFormalCharge())
                    strained = len(small - match)
                    ringy = {k for r in c.GetRingInfo().AtomRings() for k in r}
                    # complexity outside the pattern: +1 for a ring / aromatic / multiply bonded atom, +2 for a heteroatom
                    exotic = sum((1 if (a.GetIsAromatic() or a.GetIdx() in ringy or any(b.GetBondTypeAsDouble() > 1 for b in a.GetBonds())) else 0)
                                 + (2 if a.GetSymbol() != 'C' else 0) for a in atoms if a.GetIdx() not in match)
                    ring_in = 0 if ring_pattern[i] else len(match & ringy)       # a ring the pattern does not ask for
                    hsp2 = sum(c.GetAtomWithIdx(k).GetTotalNumHs() for k in match if c.GetAtomWithIdx(k).GetSymbol() == 'C'
                               and any(b.GetBondTypeAsDouble() == 2 and b.GetOtherAtom(c.GetAtomWithIdx(k)).GetSymbol() in ('O', 'N', 'S') for b in c.GetAtomWithIdx(k).GetBonds()))
                    others = len([j for j in every if j != i and j not in ignore])
                    cand[i].append((charged, strained, ring_in, exotic, hsp2, c.GetNumHeavyAtoms(), len(c.GetSubstructMatches(queries[i])), others,
                                    len(csmi), csmi, label, len(match)))
        out = {}
        for i, cs in cand.items():
            roomy = [x for x in cs if x[5] >= x[11] + 1] or cs                 # at least one heavy atom beyond the pattern
            best = min(roomy); out[i] = [best[10], best[9]]
        return out

    seen, base = set(), []
    for label, smi in base_pool() + recorded_pool():
        if smi not in seen: seen.add(smi); base.append((label, smi))
    cache.update(search(base, todo))
    missing = [i for i in todo if i not in cache]
    if missing:
        log(f'  {len(missing)} fragment(s) found no model among the common molecules; trying every recorded molecule: ' + ', '.join(name[i] for i in missing))
        cache.update(search([(l, s) for l, s in recorded_pool(minseen=1)], set(missing)))
    for i in todo:
        if i not in cache: log(f'  no model compound for {name[i]}')
    os.makedirs(os.path.dirname(MODELS), exist_ok=True)
    json.dump(cache, open(MODELS, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    return {i: cache[i] for i in queries if i in cache}


# ---------------------------------------------------------------- geometry
def _embed(mh, n=8):
    from rdkit.Chem import AllChem
    ps = AllChem.ETKDGv3(); ps.randomSeed = 0xF00D; ps.useRandomCoords = True; ps.pruneRmsThresh = 0.3
    cids = list(AllChem.EmbedMultipleConfs(mh, n, ps))
    if not cids:
        ps.useBasicKnowledge = False; ps.enforceChirality = False
        cids = list(AllChem.EmbedMultipleConfs(mh, n, ps))
    if not cids: raise RuntimeError('could not embed')
    from rdkit.Chem import rdForceFieldHelpers as FF
    res = []
    for c in cids:
        try:
            if FF.MMFFHasAllMoleculeParams(mh): AllChem.MMFFOptimizeMolecule(mh, confId=c, maxIters=500)
            else: AllChem.UFFOptimizeMolecule(mh, confId=c, maxIters=500)
        except Exception: pass
        res.append(c)
    return res


def _xtb(nums, xyz, charge):
    from tblite.interface import Calculator
    c = Calculator('GFN2-xTB', nums, xyz, charge=float(charge)); c.set('verbosity', 0)
    r = c.singlepoint()
    return float(r.get('energy')), r.get('gradient')


def optimise(mh, charge):
    """lowest-energy GFN2-xTB geometry (Bohr) among the first conformers; returns (numbers, xyz Bohr, energy Eh)"""
    import numpy as np
    from scipy.optimize import minimize
    nums = np.array([a.GetAtomicNum() for a in mh.GetAtoms()])
    cids = _embed(mh)
    pre = []
    for c in cids:
        xyz = mh.GetConformer(c).GetPositions() / BOHR
        try: pre.append((_xtb(nums, xyz, charge)[0], c))
        except Exception: pass
    if not pre: raise RuntimeError('xtb failed on every conformer')
    best = None
    for e0, c in sorted(pre)[:3]:
        x0 = mh.GetConformer(c).GetPositions() / BOHR
        fun = lambda x: (lambda eg: (eg[0], np.asarray(eg[1]).ravel()))(_xtb(nums, x.reshape(-1, 3), charge))
        r = minimize(fun, x0.ravel(), jac=True, method='L-BFGS-B', options=dict(maxiter=500, gtol=1e-4, ftol=1e-11))
        if best is None or r.fun < best[0]: best = (r.fun, r.x.reshape(-1, 3))
    return nums, best[1], best[0]


# ---------------------------------------------------------------- DFT
def dft(nums, xyz, charge, polarisability=True):
    import numpy as np
    from pyscf import gto, dft as pdft, lib, scf
    sym = {1: 'H', 3: 'Li', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl',
           19: 'K', 35: 'Br', 50: 'Sn', 53: 'I'}
    atoms = [(sym[int(z)], tuple(p)) for z, p in zip(nums, xyz)]
    ecp = {sym[int(z)]: 'def2-svp' for z in set(nums) if int(z) >= 37}
    mol = gto.M(atom=atoms, unit='Bohr', basis='def2-svp', ecp=ecp or None, charge=int(charge), spin=0, verbose=0)

    lib.num_threads(int(os.environ.get('FP_THREADS', '1')))

    def make(h=None, tol=1e-7):
        mf = pdft.RKS(mol); mf.xc = 'b3lyp'; mf.grids.level = 1; mf.conv_tol = tol; mf.max_cycle = 150
        mf = mf.DDCOSMO(); w = mf.with_solvent; w.eps = EPS_WATER
        w.lmax = 4; w.lebedev_order = 11; w.grids.atom_grid = (30, 110)      # coarser continuum grid: 2x faster, orbital energies move < 5 meV
        if h is not None: mf.get_hcore = lambda *a, h=h: h
        return mf

    def solve(mf, dm0=None):
        mf.kernel(dm0=dm0)
        if not mf.converged:
            mf.level_shift = 0.3; mf.kernel(dm0=dm0)
        if not mf.converged:
            mf = scf.newton(mf); mf.kernel(dm0=dm0)
        if not mf.converged: raise RuntimeError('SCF did not converge')
        return mf

    mf = solve(make())
    nocc = int(round(mol.nelectron / 2)); homo = nocc - 1
    e = mf.mo_energy * HARTREE_EV
    S = mol.intor('int1e_ovlp'); w, v = np.linalg.eigh(S); Sh = (v * np.sqrt(w)) @ v.T
    Co = Sh @ mf.mo_coeff                                   # Lowdin-orthogonalised coefficients
    sl = mol.aoslice_by_atom()
    share = lambda k: np.array([float((Co[p0:p1, k] ** 2).sum()) for (_, _, p0, p1) in sl])
    hs, ls = share(homo), share(homo + 1)
    dm = mf.make_rdm1()
    pop = np.array([float((Co[p0:p1, :nocc] ** 2).sum()) * 2 for (_, _, p0, p1) in sl])
    charges = np.array([mol.atom_charge(i) for i in range(mol.natm)]) - pop
    r1 = mol.intor('int1e_r', comp=3)
    mu = (mol.atom_charges()[:, None] * mol.atom_coords()).sum(0) - np.einsum('xij,ji->x', r1, dm)
    out = dict(E_homo=float(e[homo]), E_lumo=float(e[homo + 1]), homo_share=hs.tolist(), lumo_share=ls.tolist(), lowdin=charges.tolist(),
               converged=bool(mf.converged), E_homo_1=float(e[homo - 1]) if homo > 0 else None, E_lumo_1=float(e[homo + 2]) if len(e) > homo + 2 else None,
               dipole_D=float(np.linalg.norm(mu)) * 2.541746)
    if polarisability:                                      # forward difference in a field F along each axis, warm-started
        r = mol.intor('int1e_r', comp=3); h0 = mf.get_hcore(); F = 2e-3; alpha = []
        for a in range(3):
            m2 = solve(make(h0 + F * r[a], tol=1e-9), dm0=dm)
            alpha.append(-(float(np.einsum('ij,ji->', r[a], m2.make_rdm1())) - float(np.einsum('ij,ji->', r[a], dm))) / F)
        out['alpha_au'] = float(np.mean(alpha)); out['alpha_A3'] = float(np.mean(alpha)) * BOHR ** 3
    return out


# ---------------------------------------------------------------- descriptors
def oxidation_state(mh, x):
    a = mh.GetAtomWithIdx(x); e0 = en(a.GetSymbol()); v = float(a.GetFormalCharge())
    if e0 is None: return None
    for b in a.GetBonds():
        o = b.GetOtherAtom(a); e1 = en(o.GetSymbol())
        order = 1.5 if b.GetIsAromatic() else b.GetBondTypeAsDouble()
        if e1 is None: continue
        if e1 > e0: v += order
        elif e1 < e0: v -= order
    return v


def buried_volume(mh, conf_xyz_A, x, radius=3.5, scale=1.17, step=0.1):
    """percent of the sphere (radius A) round atom x that the other atoms' scaled van der Waals spheres fill (SambVca style)"""
    import numpy as np
    c = conf_xyz_A[x]
    g = np.arange(-radius, radius + 1e-9, step)
    P = np.stack(np.meshgrid(g, g, g, indexing='ij'), -1).reshape(-1, 3)
    P = P[(P ** 2).sum(1) <= radius ** 2] + c
    occ = np.zeros(len(P), bool)
    for i, a in enumerate(mh.GetAtoms()):
        if i == x: continue
        r = (EL[a.GetSymbol()].get('rvdw') or 1.7) * scale
        occ |= ((P - conf_xyz_A[i]) ** 2).sum(1) <= r * r
    return 100.0 * occ.mean()


def sigma_row(name):
    import pka as PK
    key = SIGMA_OF.get(name)
    if not key: return None
    return next((s for s in PK.SIGMA if s[0] == key), None)


def acid_base_partners():
    rows = sheet('reaction_rules.xlsx', 'Acids_Bases')
    out = []
    for r in rows:
        try: p = float(str(r['pKa (acid: its own; base: of its conjugate acid)']).replace('−', '-').split()[0].lstrip('≈~'))
        except Exception: continue
        out.append((str(r['Name']), str(r['Role']).strip().lower(), p))
    return out


def describe(fragment, model, label, defaults, partners, frags_q, quick=False):
    """one row of properties for one fragment"""
    import numpy as np
    from rdkit import Chem
    import effects as EF
    m = next((c for c in components(model)), None)
    q = frags_q[fragment['Fragment ID']]
    mh = Chem.AddHs(m)
    charge = sum(a.GetFormalCharge() for a in m.GetAtoms())
    ents = [e for e in EF.effects_of(m, [(i, qq) for i, qq in frags_q.items()]) if e['fragment'] == fragment['Fragment ID']]
    if not ents: raise RuntimeError(f'{fragment["Fragment name"]}: effects.py finds no reacting atom in {model}')
    ent = min(ents, key=lambda e: e['atom']); x = ent['atom']
    a = mh.GetAtomWithIdx(x); sym = a.GetSymbol()
    row = collections.OrderedDict()
    row['Fragment ID'] = fragment['Fragment ID']; row['Fragment name'] = fragment['Fragment name']
    row['model_smiles'] = Chem.MolToSmiles(m); row['model_source'] = label; row['reacting_atom'] = f'{sym}{x + 1}'
    # electronic
    hyb = {'SP': 'sp', 'SP2': 'sp2', 'SP3': 'sp3', 'SP3D': 'sp3d', 'SP3D2': 'sp3d2'}.get(str(a.GetHybridization()), str(a.GetHybridization()).lower())
    row['element'] = sym; row['hybridisation'] = hyb + (' (aromatic)' if a.GetIsAromatic() else '')
    row['EN'] = en(sym)
    nb = [(n.GetSymbol(), en(n.GetSymbol())) for n in a.GetNeighbors()]
    ens = [e for _, e in nb if e is not None]
    row['EN_neighbours'] = round(sum(ens) / len(ens), 2) if ens and row['EN'] is not None else None
    pol = max(((e - row['EN'], s) for s, e in nb if e is not None and row['EN'] is not None), key=lambda t: abs(t[0]), default=None)
    row['bond_polarity'] = None if pol is None else round(pol[0], 2); row['bond_polarity_partner'] = None if pol is None else pol[1]
    ch = ent['electronic']['charge']; row['charge_gasteiger'] = None if ch is None or ch != ch else round(float(ch), 3)
    row['oxidation_state'] = oxidation_state(mh, x)
    sg = sigma_row(fragment['Fragment name'])
    row['sigma_m'], row['sigma_p'], row['sigma_p_minus'] = (sg[2], sg[3], sg[4]) if sg else (None, None, None)
    if sg:
        ind = '−I' if sg[2] > 0.1 else '+I' if sg[2] < -0.05 else '≈0 I'
        d = sg[3] - sg[2]; res = '−R' if d >= 0.05 else '+R' if d <= -0.10 else '≈0 R'
        row['sigma_character'] = f'{ind} / {res}'
    else: row['sigma_character'] = None
    # steric
    st = ent['steric']; row['steric_class'] = st['cls']; row['steric_score'] = st['score']
    cn = [n for n in m.GetAtomWithIdx(x).GetNeighbors() if n.GetSymbol() == 'C']
    row['carbon_neighbours'] = len(cn)
    row['substitution'] = (['methyl', 'primary', 'secondary', 'tertiary', 'quaternary'][len(cn)] if sym == 'C' and hyb == 'sp3' and len(cn) <= 4 else None)
    row['taft_Es_min'] = ent['taft']['Es_min']; row['taft_sigma_star_sum'] = ent['taft']['sigma_star']
    ring_nb = [n for n in m.GetAtomWithIdx(x).GetNeighbors() if n.GetIsAromatic()]
    if m.GetAtomWithIdx(x).GetIsAromatic():
        ri = m.GetRingInfo(); ring = next((set(r) for r in ri.AtomRings() if x in r and all(m.GetAtomWithIdx(k).GetIsAromatic() for k in r)), set())
        ortho = [n.GetIdx() for n in m.GetAtomWithIdx(x).GetNeighbors() if n.GetIdx() in ring]
        row['ortho_substituents'] = sum(1 for o in ortho for n in m.GetAtomWithIdx(o).GetNeighbors() if n.GetIdx() not in ring)
    else: row['ortho_substituents'] = None
    # bonds at the reacting atom (average enthalpy, typical length)
    bonds = []
    for b in a.GetBonds():
        o = b.GetOtherAtom(a); order = 1.5 if b.GetIsAromatic() else b.GetBondTypeAsDouble()
        hy = lambda t: {'SP': 'sp', 'SP2': 'sp2', 'SP3': 'sp3'}.get(str(t.GetHybridization()))
        h1 = 'ar' if (a.GetIsAromatic() and sym == 'C' and a.IsInRingSize(6)) else hy(a)
        h2 = 'ar' if (o.GetIsAromatic() and o.GetSymbol() == 'C' and o.IsInRingSize(6)) else hy(o)
        kj = bond_enthalpy(sym, o.GetSymbol(), order); ln = bond_length(sym, o.GetSymbol(), order, h1, h2)
        mark = {1: '–', 2: '=', 3: '≡', 1.5: '∷'}[order]
        bonds.append((f'{sym}{mark}{o.GetSymbol()}', order, kj, ln))
    uniq = list(dict.fromkeys(bonds))
    row['bonds_text'] = '; '.join(f'{n} ' + (f'{kj:g} kJ/mol' if kj else 'aromatic' if o == 1.5 else 'no tabulated enthalpy') + (f', {ln:g} Å' if ln else '') for n, o, kj, ln in uniq)
    # the weakest bond is named only when every non-aromatic bond at the atom has a tabulated enthalpy (C-Mg, C-Sn, Al-H have none)
    weak = None if any(b[2] is None and b[1] != 1.5 for b in uniq) else min((b for b in uniq if b[2]), key=lambda b: b[2], default=None)
    row['weakest_bond'] = None if not weak else f'{weak[0]} {weak[2]:g} kJ/mol'; row['weakest_bond_kJ'] = None if not weak else weak[2]
    row['weakest_bond_length'] = None if not weak else weak[3]
    # acid / base partners from the app's own table
    pka, pkah = defaults.get(fragment['Fragment ID'], (None, None))
    row['pKa'] = pka; row['pKaH'] = pkah
    bases = [(p, n) for n, role, p in partners if role == 'base' and pka is not None and p >= pka]
    acids = [(p, n) for n, role, p in partners if role == 'acid' and pkah is not None and p <= pkah]
    row['weakest_base_that_deprotonates'] = (None if pka is None else (f'{min(bases)[1]} (ΔpKa {min(bases)[0] - pka:+.1f})' if bases else 'none in the table'))
    row['weakest_acid_that_protonates'] = (None if pkah is None else (f'{max(acids)[1]} (ΔpKa {pkah - max(acids)[0]:+.1f})' if acids else 'none in the table'))
    if quick: return row
    # geometry, steric volume, DFT
    nums, xyz, exb = optimise(mh, charge)
    row['buried_volume_pct'] = round(buried_volume(mh, xyz * BOHR, x), 1)
    d = dft(nums, xyz, charge)
    eh, el = d['E_homo'], d['E_lumo']; eta = (el - eh) / 2
    mu = (eh + el) / 2
    row['E_homo_eV'] = round(eh, 2); row['E_lumo_eV'] = round(el, 2); row['gap_eV'] = round(el - eh, 2)
    hs, ls = np.array(d['homo_share']), np.array(d['lumo_share'])
    row['homo_share_pct'] = round(100 * hs[x], 1); row['lumo_share_pct'] = round(100 * ls[x], 1)
    lab = lambda i: f'{mh.GetAtomWithIdx(int(i)).GetSymbol()}{int(i) + 1}'
    row['homo_top'] = f'{lab(hs.argmax())} {100 * hs.max():.0f}%'; row['lumo_top'] = f'{lab(ls.argmax())} {100 * ls.max():.0f}%'
    eh_, el_ = row['homo_share_pct'], row['lumo_share_pct']
    row['frontier_character'] = ('electrophile-like (LUMO)' if el_ >= 15 and el_ >= 1.5 * eh_ else 'nucleophile-like (HOMO)' if eh_ >= 15 and eh_ >= 1.5 * el_
                                 else 'mixed' if max(eh_, el_) >= 15 else 'frontier orbitals sit elsewhere')
    row['hardness_eV'] = round(eta, 2); row['softness_per_eV'] = round(1 / eta, 3)
    row['electrophilicity_eV'] = round(mu * mu / (2 * eta), 2)
    row['local_softness_minus'] = round(hs[x] / eta, 3); row['local_softness_plus'] = round(ls[x] / eta, 3)
    row['charge_lowdin'] = round(d['lowdin'][x], 3)
    row['polarisability_A3'] = round(d['alpha_A3'], 2)
    row['raw'] = dict(atoms=[at.GetSymbol() for at in mh.GetAtoms()], reacting_index=x, xyz_A=(xyz * BOHR).round(4).tolist(), E_xtb_Eh=exb,
                      E_homo_eV=eh, E_lumo_eV=el, E_homo_1_eV=d['E_homo_1'], E_lumo_1_eV=d['E_lumo_1'], dipole_D=round(d['dipole_D'], 3),
                      homo_share=np.round(hs, 4).tolist(), lumo_share=np.round(ls, 4).tolist(), lowdin=np.round(d['lowdin'], 4).tolist(),
                      alpha_au=d['alpha_au'], method='B3LYP/def2-SVP, ddCOSMO water (lmax 4, Lebedev 11), geometry GFN2-xTB')
    return row


def run_one(args):
    fragment, (label, model), defaults, partners, frags_smarts, quick = args
    os.environ['OMP_NUM_THREADS'] = '1'
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
    t0 = time.time()
    fq = {i: Chem.MolFromSmarts(s) for i, s in frags_smarts.items()}
    try:
        row = describe(fragment, model, label, defaults, partners, fq, quick)
        row['seconds'] = round(time.time() - t0, 1); return row
    except Exception as ex:
        return {'Fragment ID': fragment['Fragment ID'], 'Fragment name': fragment['Fragment name'], 'model_smiles': model, 'error': f'{type(ex).__name__}: {ex}'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', default=''); ap.add_argument('--workers', type=int, default=1)
    ap.add_argument('--quick', action='store_true', help='skip geometry and DFT (tables and rules only)')
    ap.add_argument('--rechoose', action='store_true', help='forget the cached model compounds and choose again')
    ap.add_argument('--resume', action='store_true', help='continue a run that stopped: reuse the finished rows of <out>.partial.jsonl')
    ap.add_argument('--reverse', action='store_true', help='work through the fragments from the last to the first (a second process can share a run)')
    ap.add_argument('--refresh', action='store_true', help='recompute the table-based columns of the existing output and keep its DFT values')
    ap.add_argument('--out', default=OUT)
    a = ap.parse_args()
    import pka as PK
    frags = fragments(); t0 = time.time()
    models = choose_models(frags, rechoose=a.rechoose)
    print(f'model compounds chosen for {len(models)} of {len(frags)} fragments ({time.time() - t0:.0f} s)', flush=True)
    todo = [f for f in frags if a.only.lower() in str(f['Fragment name']).lower() and f['Fragment ID'] in models]
    defaults, partners = PK.fragment_defaults(), acid_base_partners()
    smarts = {f['Fragment ID']: f['SMARTS pattern'] for f in frags}
    jobs = [(f, models[f['Fragment ID']], defaults, partners, smarts, a.quick) for f in todo]
    if a.reverse: jobs.reverse()
    if a.refresh:
        old = {r['Fragment ID']: r for r in json.load(open(a.out, encoding='utf-8'))}
        rows = []
        for j in jobs:
            fid = j[0]['Fragment ID']
            new = run_one(j[:5] + (True,))
            rows.append({**old.get(fid, {}), **new} if 'error' not in new else new)
        order = {f['Fragment ID']: k for k, f in enumerate(frags)}; rows.sort(key=lambda r: order[r['Fragment ID']])
        json.dump(rows, open(a.out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
        print(f'refreshed the table-based columns of {len(rows)} fragment(s) in {a.out}'); return
    part = a.out + '.partial.jsonl'
    rows = []
    if a.resume and os.path.exists(part):
        rows = [json.loads(l) for l in open(part, encoding='utf-8') if l.strip()]
        rows = [r for r in rows if not r.get('error')]
        done = {r['Fragment ID'] for r in rows}; jobs = [j for j in jobs if j[0]['Fragment ID'] not in done]
        print(f'resuming: {len(rows)} done, {len(jobs)} to do', flush=True)
    else:
        open(part, 'w').close()
    total = len(rows) + len(jobs)

    def keep(r):
        rows.append(r)
        open(part, 'a', encoding='utf-8').write(json.dumps(r, ensure_ascii=False) + '\n')
        print(f'  {len(rows)}/{total} {r["Fragment name"][:40]:40s} ' + (r.get('error') or f'{r.get("seconds")} s'), flush=True)
    if a.workers > 1:
        from multiprocessing import Pool
        with Pool(a.workers) as pool:
            for r in pool.imap_unordered(run_one, jobs): keep(r)
    else:
        for j in jobs: keep(run_one(j))
    order = {f['Fragment ID']: k for k, f in enumerate(frags)}
    rows.sort(key=lambda r: order[r['Fragment ID']])
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(rows, open(a.out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    bad = [r for r in rows if r.get('error')]
    print(f'wrote {a.out}: {len(rows)} fragment(s), {len(bad)} failed ({time.time() - t0:.0f} s)')
    for r in bad: print('  FAILED', r['Fragment name'], r['error'])


if __name__ == '__main__':
    main()
