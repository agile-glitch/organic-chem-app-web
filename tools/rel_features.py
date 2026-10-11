"""Features of a reaction for the template-relevance model: what the page knows before it predicts.

    python tools/rel_features.py --tag v1 --first 103                # the feature table of every learned reaction, original 103 fragments (about 10 min)
    python tools/rel_features.py --tag v2 --solvents                 # all fragments of the sheet now + one bit per recorded solvent (tools/rel_conditions.py first)
    python tools/rel_features.py --tag v2 --check                    # print the features of a few reactions, to compare with the page's

Input of the model = one bit vector per reaction (layout in LAYOUT, written next to the table so the trainer and the page read
the same one):
  sub_fp     Morgan fingerprint (radius 2, 2048 bits, bit = hash mod 2048) of every starting material, OR-ed together
  reag_fp    the same (512 bits) for every part of the reagent (ions and metals included)
  reag_id    one bit per distinct reagent part, chosen by FNV-1a of its InChIKey (1024 bits): K2CO3 and Cs2CO3 differ here
  roles      the Reagents sheet's role (and second role) of each part: 26 bits (base, oxidant, ligand ...)
  frag_sub   per fragment of the Fragments sheet (molecule_data.xlsx), 2 bits: found in the starting materials at least once /
             at least twice (the functional groups the molecules carry; the more fragments are learned the finer this is)
  frag_reag  per fragment, 1 bit: found in a part of the reagent
  meta       number of starting materials (1, 2, 3, 4+), a reagent is given, number of reagent parts (0, 1, 2, 3-4, 5+)
The page computes the same bits (js/relevance.js) from the molecules, the reagent and the Fragments sheet.
"""
import argparse, json, os, sys, tempfile, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C

FP_SUB, FP_REAG, ID_BITS = 2048, 512, 1024
META_BITS = 10
TMP = tempfile.gettempdir()


def table_paths(tag='v1'):
    return os.path.join(TMP, f'organic_chem_relfeat_{tag}.npy'), os.path.join(TMP, f'organic_chem_relfeat_{tag}.json')


def role_list():
    roles = json.load(open(os.path.join(HERE, 'data', 'reagents.json'), encoding='utf-8'))
    return sorted({x['role'] for x in roles if x.get('role')}), {x['smiles']: x for x in roles}


def fragment_list(first=0):
    """[(Fragment ID, name, SMARTS)] of the Fragments sheet, sorted by ID (the order of the model's inputs); first = only F001..F<first>
    (103 = the sheet as it was when the learned reactions were made, before tools/add_fragments.py)"""
    import learn_reactions as L
    return [f for f in sorted(L.fragments(), key=lambda f: f[0]) if not first or int(f[0][1:]) <= first]


def solvent_names():
    import learn_reactions as L
    return sorted(set(L.solvents().values()), key=str.lower)


def make_layout(n_frag, n_roles, n_solv=0):
    parts = [('sub_fp', FP_SUB), ('reag_fp', FP_REAG), ('reag_id', ID_BITS), ('roles', n_roles), ('frag_sub', 2 * n_frag), ('frag_reag', n_frag), ('meta', META_BITS)]
    if n_solv: parts.append(('solv', n_solv))                     # one bit per recorded solvent (the conditions box of the page)
    off, layout = 0, {}
    for name, n in parts:
        layout[name] = [off, n]; off += n
    layout['bits'] = off
    layout['bytes'] = (off + 7) // 8
    return layout


_W = {}


def _init(frag_smarts, roles_by_smiles, role_names, solvents=()):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator as FG
    RDLogger.DisableLog('rdApp.*')
    _W['frags'] = [Chem.MolFromSmarts(s) for s in frag_smarts]
    _W['roles'], _W['role_idx'] = roles_by_smiles, {r: i for i, r in enumerate(role_names)}
    _W['g_sub'], _W['g_reag'] = FG.GetMorganGenerator(radius=2, fpSize=FP_SUB), FG.GetMorganGenerator(radius=2, fpSize=FP_REAG)
    _W['reag_cache'] = {}
    _W['solv_idx'] = {n: i for i, n in enumerate(solvents)}


def mol_bits(gen, m):
    return np.asarray(gen.GetFingerprintAsNumPy(m), dtype=np.uint8)


def frag_counts(m, cap=2):
    return [len(m.GetSubstructMatches(q, uniquify=True, maxMatches=cap)) if q is not None else 0 for q in _W['frags']]


def reagent_block(reagent):
    """(reag_fp bits, reag_id bits, role bits, fragment bits) of a reagent string; cached"""
    hit = _W['reag_cache'].get(reagent)
    if hit is not None: return hit
    from rdkit import Chem
    fp = np.zeros(FP_REAG, np.uint8); ids = np.zeros(ID_BITS, np.uint8); rl = np.zeros(len(_W['role_idx']), np.uint8)
    fr = np.zeros(len(_W['frags']), np.uint8)
    parts = sorted(set(x for x in reagent.split('.') if x))
    for p in parts:
        m = Chem.MolFromSmiles(p)
        if m is None: continue
        fp |= mol_bits(_W['g_reag'], m)
        try: key = Chem.MolToInchiKey(m)
        except Exception: key = ''
        ids[C.fnv1a32(key or p) % ID_BITS] = 1
        info = _W['roles'].get(p)
        if info:
            for role in (info.get('role'), info.get('also')):
                if role in _W['role_idx']: rl[_W['role_idx'][role]] = 1
        for i, c in enumerate(frag_counts(m, 1)):
            if c: fr[i] = 1
    out = (fp, ids, rl, fr, len(parts))
    if len(_W['reag_cache']) < 60000: _W['reag_cache'][reagent] = out
    return out


def featurise(subs, reagent, layout, conds=''):
    """the bit vector (packed to bytes) of one reaction; conds = its recorded solvents, 'name, name' (only used when the layout has a solv block)"""
    from rdkit import Chem
    bits = np.zeros(layout['bytes'] * 8, np.uint8)
    nfrag = len(_W['frags'])
    sub_fp = np.zeros(FP_SUB, np.uint8); counts = np.zeros(nfrag, np.int32)
    for s in subs:
        m = Chem.MolFromSmiles(s)
        if m is None: continue
        sub_fp |= mol_bits(_W['g_sub'], m)
        counts += np.array(frag_counts(m), np.int32)
    fp, ids, rl, fr, nparts = reagent_block(reagent)
    o = layout
    bits[o['sub_fp'][0]:o['sub_fp'][0] + FP_SUB] = sub_fp
    bits[o['reag_fp'][0]:o['reag_fp'][0] + FP_REAG] = fp
    bits[o['reag_id'][0]:o['reag_id'][0] + ID_BITS] = ids
    bits[o['roles'][0]:o['roles'][0] + len(rl)] = rl
    f0 = o['frag_sub'][0]
    bits[f0:f0 + nfrag] = counts >= 1
    bits[f0 + nfrag:f0 + 2 * nfrag] = counts >= 2
    bits[o['frag_reag'][0]:o['frag_reag'][0] + nfrag] = fr
    m0 = o['meta'][0]
    ns = len(subs)
    bits[m0 + min(ns, 4) - 1] = 1                                # 1, 2, 3, 4+ starting materials
    bits[m0 + 4] = 1 if reagent else 0
    bits[m0 + 5 + (0 if nparts == 0 else 1 if nparts == 1 else 2 if nparts == 2 else 3 if nparts <= 4 else 4)] = 1
    if 'solv' in o:
        for name in [x for x in (conds or '').split(', ') if x]:
            if name in _W['solv_idx']: bits[o['solv'][0] + _W['solv_idx'][name]] = 1
    return np.packbits(bits)


def _chunk(args):
    items, layout = args
    return np.stack([featurise(subs, reagent, layout, conds) for subs, reagent, conds in items])


def build(tag='v1', workers=0, limit=0, first=0, solvents=False):
    """first: only the fragments F001..F<first> of the sheet (0 = all of them now)"""
    D = C.read_rows()
    rows = D['rows'][:limit] if limit else D['rows']
    fr = fragment_list(first)
    roles, by_smiles = role_list()
    solv = solvent_names() if solvents else []
    conds = None
    if solvents:
        import rel_conditions as RC
        conds = RC.load()
    layout = make_layout(len(fr), len(roles), len(solv))
    path, meta_path = table_paths(tag)
    workers = workers or max(1, (os.cpu_count() or 2) - 1)
    N, CH = len(rows), 4000
    mm = np.lib.format.open_memmap(path, mode='w+', dtype=np.uint8, shape=(N, layout['bytes']))
    t0 = time.time()
    chunks = ((([(r[1], r[2], (conds[i] if conds else '')) for i, r in enumerate(rows[a:a + CH], a)]), layout) for a in range(0, N, CH))
    with Pool(workers, initializer=_init, initargs=([f[2] for f in fr], by_smiles, roles, solv)) as pool:
        for k, arr in enumerate(pool.imap(_chunk, chunks)):
            mm[k * CH:k * CH + len(arr)] = arr
            if (k + 1) % 20 == 0: print(f'  {min(N, (k + 1) * CH):,} of {N:,} ({time.time() - t0:.0f} s)', flush=True)
    mm.flush()
    meta = {'rows': N, 'layout': layout, 'fragments': [f[0] for f in fr], 'roles': roles, 'solvents': solv, 'first': rows[0][0], 'last': rows[-1][0], 'built': time.strftime('%Y-%m-%d %H:%M')}
    json.dump(meta, open(meta_path, 'w'), indent=1)
    print(f'wrote {path} ({N:,} x {layout["bytes"]} bytes, {layout["bits"]} bits, {len(fr)} fragments, {len(solv)} solvents) in {time.time() - t0:.0f} s')


def describe(tag='v1', n=3):
    path, meta_path = table_paths(tag)
    meta = json.load(open(meta_path)); L = meta['layout']
    X = np.load(path, mmap_mode='r'); D = C.read_rows()
    for k in (0, 500000, 1000000)[:n]:
        v = np.unpackbits(X[k])[:L['bits']]
        print(D['rows'][k][0], D['rows'][k][1], '|', D['rows'][k][2])
        for name in ('sub_fp', 'reag_fp', 'reag_id', 'roles', 'frag_sub', 'frag_reag', 'meta'):
            a, w = L[name]
            print(f'   {name:9s} {int(v[a:a + w].sum()):4d} bits set of {w}')
        fs = meta['fragments']; f0 = L['frag_sub'][0]
        print('   fragments in the starting materials:', [fs[i] for i in range(len(fs)) if v[f0 + i]])


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', default='v1'); ap.add_argument('--workers', type=int, default=0)
    ap.add_argument('--limit', type=int, default=0, help='only the first N reactions (a quick test)')
    ap.add_argument('--solvents', action='store_true', help='add one bit per recorded solvent (needs tools/rel_conditions.py first)')
    ap.add_argument('--first', type=int, default=0, help='only the fragments F001..F<first> (103 = the v1 table, made before the 66 learned fragments were added)')
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    if a.check: describe(a.tag)
    else: build(a.tag, a.workers, a.limit, a.first, solvents=a.solvents)
