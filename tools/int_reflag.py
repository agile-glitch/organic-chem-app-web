"""Recompute the interference flags of stored contests (tools/int_contests.py keeps each outcome's marked products) with the CURRENT pattern library.

    python tools/int_reflag.py --in int_contests_train.pkl --out int_contests_train.pkl

After a pattern of tools/interference_defs.py is corrected or added, the contests need no new template runs: the flags are read off the stored marked
products again (a minute or two on 7 workers), and tools/int_train.py is run again.
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
    for c in d['contests']:
        for o in c['outs']:
            if len(o) < 4: raise SystemExit('these contests do not keep the marked products: make them again with tools/int_contests.py')
            uniq.setdefault(o[3], None)
    keys = list(uniq)
    print(f'{len(d["contests"]):,} contests, {len(keys):,} different marked products ({time.time() - t0:.0f} s)', flush=True)
    workers = a.workers or max(1, (os.cpu_count() or 2) - 1)
    chunks = [keys[i:i + 2000] for i in range(0, len(keys), 2000)]
    with Pool(workers, initializer=_init) as pool:
        got = []
        for k, r in enumerate(pool.imap(_flags, chunks), 1):
            got.extend(r)
            if k % 20 == 0: print(f'  {k * 2000:,} of {len(keys):,} ({time.time() - t0:.0f} s)', flush=True)
    fl = dict(zip(keys, got))
    for c in d['contests']:
        c['outs'] = [(o[0], fl[o[3]], o[2], o[3]) for o in c['outs']]
    d['reflagged'] = time.strftime('%Y-%m-%d %H:%M')
    pickle.dump(d, open(os.path.join(TMP, a.out), 'wb'), protocol=4)
    print(f'wrote {a.out} ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
