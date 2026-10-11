"""Which functional group reacts first? Learned from the recorded reactions, per reagent.

    python tools/rel_reactivity.py --counts --tag cur --first 103   # fragment counts (original 103 fragments) of every learned reaction's starting materials and product
    python tools/rel_reactivity.py --counts --tag new               # the same with every fragment of the sheet now
    python tools/rel_reactivity.py --table --tag new                # the reactivity table from the TRAINING patents (honest measurement)
    python tools/rel_reactivity.py --table --all --tag new          # ... from every reaction (what the page ships)
    python tools/rel_reactivity.py --export data/learned/m/prior.js --tag new   # the page's copy of the all-reactions table

For each recorded reaction the Fragments sheet (molecule_data.xlsx) says which fragments the starting materials have and which the
product has; a fragment whose count drops was CONSUMED (a fragment whose count rises was formed). Over many reactions that gives,
for a reagent, how often each fragment that is present is consumed: NaOH takes the methyl ester of a Cbz-protected amino ester, 1-2
equivalents of NaBH4 take the aldehyde before the ketone and leave the ester, HCl takes a MOM ether before an N-benzyl group ...
These are the same facts as the "Reactivity_Order" of reaction_rules.xlsx, found in the recorded reactions instead of written by hand,
so a better fragment list (more fragments, finer ones) gives a better table.

Rates are smoothed down a chain: the exact reagent -> reagents with the same roles (Reagents sheet) -> every reaction:
    p(exact) = (consumed + k p(roles)) / (present + k),   k = 5
A candidate product is then scored by what it consumed and what it left: sum of log-odds of the consumed fragments (present in the
starting materials) and sum of log(1 - rate) of the fragments it left untouched (tools/rel_fuse.py reads both).
"""
import argparse, json, os, pickle, sys, tempfile, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C

TMP = tempfile.gettempdir()
CAP = 8                                  # matches counted per molecule per fragment
K_SMOOTH = 5.0
SKIP_ROLES = {'counter-ion / spectator ion', 'solvent', 'salt (work-up or additive)', 'other', 'drying agent'}
_W = {}


def fragment_list(first=0):
    """the Fragments sheet sorted by ID; first = only F001..F<first> (103 = the sheet before tools/add_fragments.py)"""
    import learn_reactions as L
    return [f for f in sorted(L.fragments(), key=lambda f: f[0]) if not first or int(f[0][1:]) <= first]


def counts_path(tag): return os.path.join(TMP, f'organic_chem_relfragcounts_{tag}.npz')
def table_path(tag): return os.path.join(TMP, f'organic_chem_relprior_{tag}.pkl')


def _init(frag_smarts):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
    _W['q'] = [Chem.MolFromSmarts(s) for s in frag_smarts]


def mol_counts(m, cap=CAP):
    return np.array([min(cap, len(m.GetSubstructMatches(q, uniquify=True, maxMatches=cap))) if q is not None else 0 for q in _W['q']], np.int32)


def counts_of(smiles_list):
    from rdkit import Chem
    c = np.zeros(len(_W['q']), np.int32)
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if m is not None: c += mol_counts(m)
    return np.minimum(c, 255).astype(np.uint8)


def _chunk(items):
    return np.stack([counts_of(subs) for subs, _ in items]), np.stack([counts_of([prod]) for _, prod in items])


def build_counts(tag='cur', workers=0, first=0):
    D = C.read_rows(); rows = D['rows']
    fr = fragment_list(first)
    N, CH = len(rows), 4000
    workers = workers or max(1, (os.cpu_count() or 2) - 1)
    S = np.zeros((N, len(fr)), np.uint8); P = np.zeros((N, len(fr)), np.uint8)
    t0 = time.time()
    with Pool(workers, initializer=_init, initargs=([f[2] for f in fr],)) as pool:
        chunks = ([(rows[i][1], rows[i][3]) for i in range(a, min(N, a + CH))] for a in range(0, N, CH))
        for k, (s, p) in enumerate(pool.imap(_chunk, chunks)):
            S[k * CH:k * CH + len(s)] = s; P[k * CH:k * CH + len(p)] = p
            if (k + 1) % 20 == 0: print(f'  {min(N, (k + 1) * CH):,} of {N:,} ({time.time() - t0:.0f} s)', flush=True)
    np.savez(counts_path(tag), S=S, P=P, fragments=np.array([f[0] for f in fr]))
    print(f'wrote {counts_path(tag)} ({time.time() - t0:.0f} s)')


def role_signature(reagent, roles):
    """the jobs a reagent's parts do (Reagents sheet), spectators and solvents left out ('' = no reagent; None if a part is unknown)"""
    if not reagent: return ''
    sig = set()
    for p in reagent.split('.'):
        r = roles.get(p)
        if r is None: return None
        if r['role'] not in SKIP_ROLES: sig.add(r['role'])
    return '|'.join(sorted(sig)) if sig else None


def reagent_roles():
    return {x['smiles']: x for x in json.load(open(os.path.join(HERE, 'data', 'reagents.json'), encoding='utf-8'))}


def classes_of(rows, roles):
    """per reaction: (exact reagent key or '', role signature or None)"""
    return [(r[5] if r[2] and not r[5].startswith('?') else '', role_signature(r[2], roles)) for r in rows]


def build_table(tag='cur', use_all=False, min_exact=30):
    D = C.read_rows(); rows = D['rows']
    z = np.load(counts_path(tag)); S, P = z['S'], z['P']
    nf = S.shape[1]
    roles = reagent_roles()
    cls = classes_of(rows, roles)
    split = np.array(C.split_codes(rows))
    use = np.ones(len(rows), bool) if use_all else (split == C.TRAIN)
    present = (S >= 1)
    consumed = (S > P)
    formed = (P > S)
    tab = {'exact': {}, 'role': {}, 'global': None}
    g_np = present[use].sum(0).astype(np.float64); g_nc = (present & consumed)[use].sum(0).astype(np.float64); g_nf = formed[use].sum(0).astype(np.float64)
    tab['global'] = (g_np, g_nc, g_nf, int(use.sum()))
    ex_rows, rl_rows = {}, {}
    for i in np.flatnonzero(use):
        e, r = cls[i]
        if e: ex_rows.setdefault(e, []).append(i)
        if r is not None: rl_rows.setdefault(r, []).append(i)
    for name, groups, lim in (('exact', ex_rows, min_exact), ('role', rl_rows, 1)):
        for k, idx in groups.items():
            if len(idx) < lim: continue
            idx = np.array(idx)
            tab[name][k] = (present[idx].sum(0).astype(np.float32), (present & consumed)[idx].sum(0).astype(np.float32), formed[idx].sum(0).astype(np.float32), len(idx))
    meta = {'fragments': [str(x) for x in z['fragments']], 'all': use_all, 'reactions': int(use.sum()), 'built': time.strftime('%Y-%m-%d %H:%M')}
    pickle.dump({'table': tab, 'meta': meta}, open(table_path(tag + ('_all' if use_all else '')), 'wb'), protocol=4)
    print(f'reactivity table from {int(use.sum()):,} reactions: {len(tab["exact"]):,} exact reagents, {len(tab["role"]):,} role signatures -> {table_path(tag + ("_all" if use_all else ""))}')


class Prior:
    """consumption rates per fragment for a reagent, smoothed down the chain exact -> roles -> global.
    ids = use only these fragments of the table (e.g. the original 103) while the candidates carry counts for all of them"""
    def __init__(self, tag='cur', use_all=False, ids=None):
        d = pickle.load(open(table_path(tag + ('_all' if use_all else '')), 'rb'))
        self.t, self.meta = d['table'], d['meta']
        self.source_ids = self.meta['fragments']                      # the fragment list the candidates must have been counted with
        self.cols = np.array([self.source_ids.index(i) for i in ids]) if ids else np.arange(len(self.source_ids))
        self.frag_ids = [self.source_ids[i] for i in self.cols]
        g = self.t['global']; self.g_rate = (g[1][self.cols] + 1.0) / (g[0][self.cols] + 2.0)
        self.cache = {}

    def rates(self, exact_key, role_sig):
        ck = (exact_key, role_sig)
        hit = self.cache.get(ck)
        if hit is not None: return hit
        c = self.cols
        p = self.g_rate
        r = self.t['role'].get(role_sig) if role_sig is not None else None
        if r is not None: p = (r[1][c] + K_SMOOTH * p) / (r[0][c] + K_SMOOTH)
        e = self.t['exact'].get(exact_key) if exact_key else None
        if e is not None: p = (e[1][c] + K_SMOOTH * p) / (e[0][c] + K_SMOOTH)
        p = np.clip(p, 1e-3, 1 - 1e-3)
        if len(self.cache) < 20000: self.cache[ck] = p
        return p

    def features(self, s_counts, c_counts, exact_key, role_sig):
        """(sum log-odds of the consumed fragments, sum log(1-rate) of the fragments present and left, consumed, formed) for one candidate"""
        p = self.rates(exact_key, role_sig)
        s = np.asarray(s_counts)[self.cols]; c = np.asarray(c_counts)[self.cols]
        present = s >= 1
        cons = s > c
        formed = c > s
        a = float(np.log(p[cons] / (1 - p[cons])).sum()) if cons.any() else 0.0
        b = float(np.log(1 - p[present & ~cons]).sum()) if (present & ~cons).any() else 0.0
        return a, b, int(cons.sum()), int(formed.sum())


def export_prior(tag, out, use_all=True):
    """data/learned/m/prior.js: the table the page's Relevance.priorRates reads (built from every reaction: --table --all)"""
    from rel_fuse import write_json_part
    d = pickle.load(open(table_path(tag + ('_all' if use_all else '')), 'rb'))
    t, meta = d['table'], d['meta']
    pack = lambda v: {'np': [int(x) for x in v[0]], 'nc': [int(x) for x in v[1]]}
    g = t['global']
    obj = {'k': K_SMOOTH, 'fragments': meta['fragments'], 'global': pack(g), 'role': {k: pack(v) for k, v in t['role'].items()},
           'exact': {k: pack(v) for k, v in t['exact'].items()}, 'reactions': meta['reactions'], 'built': meta['built']}
    write_json_part(out, obj)
    print(f'wrote {out}: {len(obj["exact"]):,} exact reagents, {len(obj["role"]):,} role signatures ({os.path.getsize(out) / 1e6:.1f} MB)')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--counts', action='store_true'); ap.add_argument('--table', action='store_true'); ap.add_argument('--all', action='store_true')
    ap.add_argument('--export', default='', help='write the page\'s prior (data/learned/m/prior.js) from the --all table')
    ap.add_argument('--tag', default='cur'); ap.add_argument('--workers', type=int, default=0); ap.add_argument('--first', type=int, default=0, help='with --counts: only the fragments F001..F<first>'); ap.add_argument('--train-only', action='store_true', help='with --export: use the training-patents table instead of the all-reactions one')
    a = ap.parse_args()
    if a.counts: build_counts(a.tag, a.workers, a.first)
    if a.table: build_table(a.tag, a.all)
    if a.export: export_prior(a.tag, a.export, not a.train_only)
