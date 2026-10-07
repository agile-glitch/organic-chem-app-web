"""Compare ways of choosing the candidate templates for new molecules, on the same held-out reactions.

    python tools/selectivity_experiments.py [--test 4000]

The page tries the templates learned with the exact reagent (else with its parts). When none of them gives the
recorded product the prediction cannot be right, whatever the ranking. These variants try to raise that ceiling:
  exact      the page's list (exact reagent, else parts)
  cascade    exact reagent; only when none of its templates fits, the parts; then the same roles
  union      exact reagent + its parts, always
  roles      exact, else parts, else reagents with the same roles (Reagents sheet, tools/data/reagents.json:
             K2CO3 ~ Cs2CO3 as "base"; spectator ions, solvents and unknown parts ignored)
Each is ranked with the learned win rates (data/learned/w/winrates.js). Writes tools/selectivity_experiments.json.
"""
import argparse, base64, collections, gzip, json, os, sys, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import learn_selectivity as S

SKIP_ROLES = {'counter-ion / spectator ion', 'solvent', 'salt (work-up or additive)', 'other', 'drying agent'}


def role_sig(reagent, roles):
    sig = set()
    for p in reagent.split('.'):
        r = roles.get(p)
        if r is None: return None                      # an unknown part: no role signature
        if r['role'] in SKIP_ROLES: continue
        sig.add(r['role'])
    return '|'.join(sorted(sig)) if sig else None


_X = {}


def _init(smarts, idx, roles, role_index, solv_frags):
    S._init(smarts, idx, solv_frags)
    _X['roles'], _X['role_index'] = roles, role_index


def outcomes(r):
    import learn_reactions as L
    idx = S._G['idx']
    exact = S.exact_list(r, idx) if r[2] else None
    parts, _ = S.candidates([r[0], r[1], r[2], r[3], r[4], '#none#', r[6]], idx) if r[2] else ([], '')
    sig = role_sig(r[2], _X['roles']) if r[2] else None
    by_role = [x for x in (_X['role_index'].get(sig, []) if sig else []) if x[1] >= S.MIN_OBS]
    lists = {'exact': (exact or parts or [])[:S.TOP], 'union': ((exact or [])[:S.TOP] + parts[:S.TOP]),
             'roles': (exact or parts or by_role or [])[:S.TOP], '_parts': parts[:S.TOP], '_roles': by_role[:S.TOP]}
    allt = {}
    for name, lst in lists.items():
        for t, n in lst: allt.setdefault(t, n)
    extras = L.pool_extras(r[2], L._W['roles'])
    fps = L.pool_fps(r[1], extras)
    outs = {}
    for t in allt:
        o, _ = L.run_template(S._G['smarts'][t], r[1], extras, screen=fps)
        if o: outs[t] = ['.'.join(sorted(x)) for x in o]
    res = {name: [(t, n, outs[t]) for t, n in lst if t in outs] for name, lst in lists.items()}
    # cascade: the exact reagent's templates; only if none fits, its parts; only if none fits, the same roles
    res['cascade'] = res['exact'] or res['_parts'] or res['_roles']
    res['cascade-roles'] = res['exact'] or res['_roles'] or res['_parts']
    return r[0], r[3], res


def weak_cascade(res, tau, win):
    """the cascade, but a step whose best answer is weak (its template wins less than tau of the time it could react)
    does not stop the search: the next step's templates join and all are ranked together"""
    lst = []
    for stage in ('exact', '_parts', '_roles'):
        lst = lst + [x for x in res[stage] if x[0] not in {y[0] for y in lst}]
        if not lst: continue
        best = max(((win.get(t, (0, 0))[0] + 1) / (win.get(t, (0, 0))[1] + 2) for t, n, outs in lst), default=0)
        if best >= tau: break
    return lst


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--test', type=int, default=4000); args = ap.parse_args()
    import learn_reactions as L
    t0 = time.time()
    D = S.read_rows(); rows = D['rows']
    train = [r for r in rows if not S.is_test(r[0])]; test = [r for r in rows if S.is_test(r[0])]
    idx = S.build_index(train)
    roles = {x['smiles']: x for x in json.load(open(os.path.join(HERE, 'data', 'reagents.json'), encoding='utf-8'))}
    role_index = collections.defaultdict(collections.Counter)
    for r in train:
        if r[2]:
            sig = role_sig(r[2], roles)
            if sig: role_index[sig][r[4]] += 1
    role_index = {k: c.most_common() for k, c in role_index.items()}
    w = open(os.path.join(os.path.dirname(HERE), 'data', 'learned', 'w', 'winrates.js'), encoding='utf-8').read()
    W = json.loads(gzip.decompress(base64.b64decode(w.split('"gz":"')[1].split('"')[0])))
    win = {int(k): tuple(v) for k, v in W['win'].items()}
    sample = test[::max(1, len(test) // args.test)][:args.test]
    print(f'{len(sample):,} held-out reactions; role signatures: {len(role_index):,} ({time.time() - t0:.0f} s)', flush=True)
    tally = collections.defaultdict(lambda: collections.defaultdict(lambda: [0, 0, 0]))
    # 4 workers: each holds a copy of the (unlimited) index
    with Pool(4, initializer=_init, initargs=(D['smarts'], idx, roles, role_index, (L.solvents(), L.fragments()))) as pool:
        n = 0
        for lid, prod, res in pool.imap_unordered(outcomes, sample, chunksize=10):
            n += 1
            for tau in (0.1, 0.2, 0.3):
                res[f'cascade-weak{tau}'] = weak_cascade(res, tau, win)
            for name, lst in res.items():
                if name.startswith('_'): continue
                for m in S.METHODS:
                    rk = S.rank(lst, m, win)
                    t = tally[name][m]
                    if rk[:1] and prod in rk[0].split('.'): t[0] += 1
                    if any(prod in o.split('.') for o in rk[:3]): t[1] += 1
                    if any(prod in o.split('.') for o in rk): t[2] += 1
            if n % 500 == 0: print(f'  {n:,} ({time.time() - t0:.0f} s)', flush=True)
    out = {name: {m: {'top1': round(v[0] / n, 4), 'top3': round(v[1] / n, 4), 'ceiling': round(v[2] / n, 4)} for m, v in d.items()} for name, d in tally.items()}
    print(json.dumps({k: {m: v[m] for m in ('count', 'winrate-max', 'winrate-sum')} for k, v in out.items()}, indent=1))
    json.dump({'tested': n, 'variants': out}, open(os.path.join(HERE, 'selectivity_experiments.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
