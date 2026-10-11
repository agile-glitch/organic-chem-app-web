"""Recompute the interference flags of a candidate dump (tools/eval_top1.py --marks) with the CURRENT pattern library.

    python tools/int_reflag_dump.py --in val_full2m2.pkl --out val_full2m2.pkl

The dump keeps, for every candidate, its products with the changed atoms marked ('mk'); after a pattern of tools/interference_defs.py is corrected or added the
flags are read off those again (a few minutes on 7 workers) instead of running every template once more (15 minutes per split).
"""
import argparse, os, pickle, sys, tempfile, time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

TMP = tempfile.gettempdir()
_E = {}


def _init():
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    import interferences as I
    _E['lib'] = I.Library()


def _flags(chunk):
    return [tuple(_E['lib'].flags_of_smiles(s)) for s in chunk]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--in', dest='inp', required=True); ap.add_argument('--out', required=True); ap.add_argument('--workers', type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    d = pickle.load(open(os.path.join(TMP, a.inp), 'rb'))
    uniq = {}
    n_c = 0
    for cands in d['res'].values():
        for c in cands.values():
            if 'flags' not in c: continue
            if 'mk' not in c: raise SystemExit('this dump does not keep the marked products: make it again with tools/eval_top1.py --marks')
            uniq.setdefault(c['mk'], None); n_c += 1
    keys = list(uniq)
    print(f'{n_c:,} candidates, {len(keys):,} different marked products ({time.time() - t0:.0f} s)', flush=True)
    workers = a.workers or max(1, (os.cpu_count() or 2) - 1)
    chunks = [keys[i:i + 2000] for i in range(0, len(keys), 2000)]
    got = []
    with Pool(workers, initializer=_init) as pool:
        for k, r in enumerate(pool.imap(_flags, chunks), 1):
            got.extend(r)
            if k % 25 == 0: print(f'  {k * 2000:,} of {len(keys):,} ({time.time() - t0:.0f} s)', flush=True)
    fl = dict(zip(keys, got))
    import interference_defs as D
    site_keys = [s['key'] for s in D.SITE]
    kidx = {k: i for i, k in enumerate(site_keys)}
    for cands in d['res'].values():
        for c in cands.values():
            if 'flags' in c: c['flags'] = sorted(kidx[k] for k in fl[c['mk']] if k in kidx)
    d['site_keys'] = site_keys
    d['reflagged'] = time.strftime('%Y-%m-%d %H:%M')
    pickle.dump(d, open(os.path.join(TMP, a.out), 'wb'), protocol=4)
    print(f'wrote {a.out} ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
