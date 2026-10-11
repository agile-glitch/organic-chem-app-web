"""Regioselectivity contests from the recorded reactions: where a template could have reacted at several sites, which one did?

    python tools/int_contests.py --split train --n 300000 --out int_contests_train.pkl
    python tools/int_contests.py --split val --n 0 --out int_contests_val.pkl          # n 0 = every reaction of the split

For a recorded reaction the starting materials and the reagent are given to the reaction's own template with every unchanged leaf atom stripped
(tools/template_general.py level 3, the same relaxed form the page runs). Each distinct product set is an OUTCOME; the one that contains the recorded
product is the winner. A contest needs two or more outcomes and the recorded one among them. Every outcome carries the interference flags of its
reacting atoms (tools/interferences.py: the changed atoms are marked by the template, tools/int_marks.py). The learned data are positive-only
(a recorded reaction says what happened, never what failed), but a contest is a comparison inside one reaction: the competing sites were all
possible for the same template, and the recorded product shows which one the chemist got. tools/int_train.py fits the weights of the flags.
Reactions are de-duplicated on (starting materials, product): a patent often records the same step many times.
"""
import argparse, collections, os, pickle, sys, tempfile, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C
import template_general as TG

TMP = tempfile.gettempdir()
_E = {}


def _init(solv_frags, frag_smarts):
    import learn_reactions as L
    import rel_reactivity as RV
    import interferences as I
    L._init(*solv_frags, L.reagent_roles())
    RV._init(list(frag_smarts))
    _E['lib'] = I.Library()
    _E['gcache'] = {}


def run_task(task):
    """task = (row index, starting materials, reagent, recorded product, template SMARTS) -> (row index, family, consumed fragment IDs, [(hit, flags, product key)]) or None"""
    import learn_reactions as L
    import rel_reactivity as RV
    import int_run as R
    import interference_defs as D
    ri, subs, reagent, prod, smarts = task
    g = _E['gcache'].get(smarts)
    if g is None:
        g = _E['gcache'][smarts] = TG.generalise(smarts, 3) or smarts
        if len(_E['gcache']) > 20000: _E['gcache'].clear()
    extras = L.pool_extras(reagent, L._W['roles'])
    fps = L.pool_fps(subs, extras)
    outs = R.run_marked(g, subs, extras, fps, _E['lib'])
    R.clear_caches(L)
    if len(outs) < 2: return ri, None, None, len(outs)
    res = [(prod in cs, tuple(flags), '.'.join(sorted(cs)), '.'.join(m)) for cs, _u, m, flags in outs]       # (is the recorded one, flags, product key, the products with the changed atoms marked)
    if not any(r[0] for r in res): return ri, None, None, -1
    ids = [f[0] for f in RV.fragment_list()]
    s, p = RV.counts_of(subs), RV.counts_of([prod])
    consumed = [ids[i] for i in np.flatnonzero(s > p)]
    return ri, D.family_of(set(consumed)), consumed, res


def pick_rows(split_name, n, seed=11):
    D = C.read_rows(); rows = D['rows']
    split = np.array(C.split_codes(rows))
    code = {'train': C.TRAIN, 'val': C.VAL, 'test': C.TEST}[split_name]
    idx = np.flatnonzero(split == code)
    seen, keep = set(), []
    for i in np.random.default_rng(seed).permutation(idx):
        r = rows[i]
        if not r[2] and False: continue
        k = (tuple(sorted(r[1])), r[3])
        if k in seen: continue
        seen.add(k); keep.append(int(i))
        if n and len(keep) >= n: break
    return D, keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--split', choices=['train', 'val', 'test'], default='train'); ap.add_argument('--n', type=int, default=20000, help='0 = every (de-duplicated) reaction of the split')
    ap.add_argument('--out', required=True); ap.add_argument('--workers', type=int, default=0); ap.add_argument('--seed', type=int, default=11)
    a = ap.parse_args()
    t0 = time.time()
    D, keep = pick_rows(a.split, a.n, a.seed)
    rows, smarts = D['rows'], D['smarts']
    import rel_reactivity as RV
    import learn_reactions as L
    frs = RV.fragment_list()
    tasks = [(i, rows[i][1], rows[i][2], rows[i][3], smarts[rows[i][4]]) for i in keep]
    print(f'{len(tasks):,} {a.split} reactions (de-duplicated) ({time.time() - t0:.0f} s)', flush=True)
    workers = a.workers or max(1, (os.cpu_count() or 2) - 1)
    contests, stat = [], collections.Counter()
    with Pool(workers, initializer=_init, initargs=((L.solvents(), L.fragments()), [f[2] for f in frs])) as pool:
        for k, (ri, fam, consumed, res) in enumerate(pool.imap_unordered(run_task, tasks, chunksize=8), 1):
            if fam is None:
                stat['one outcome' if res == 1 else 'no outcome' if res == 0 else 'recorded product not among them' if res == -1 else 'other'] += 1
            else:
                contests.append({'ri': ri, 'family': fam, 'consumed': consumed, 'outs': res}); stat['contest'] += 1
            if k % 5000 == 0: print(f'  {k:,} of {len(tasks):,}: {stat["contest"]:,} contests ({time.time() - t0:.0f} s)', flush=True)
    pickle.dump({'contests': contests, 'stat': dict(stat), 'split': a.split, 'when': time.strftime('%Y-%m-%d %H:%M')}, open(os.path.join(TMP, a.out), 'wb'), protocol=4)
    print(dict(stat), f'-> {a.out} ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
