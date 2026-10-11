"""Give an evaluation dump (tools/eval_top1.py --dump) the fragment counts of a LONGER Fragments sheet, without running the templates again.

    python tools/rel_recount.py --dump val_full2.pkl --out val_full2_f205.pkl

A dump keeps, per reaction, the fragment counts of the starting materials and of every candidate product (tools/rel_fuse.py needs them for the
reactivity prior). When fragments are added to the sheet only the counts of the NEW fragments are missing: they are computed from the molecules
(the candidates are kept as SMILES) and appended, and the dump's fragment list is extended. The dump's old fragments must be the first ones of
the new list (new IDs always come last).
"""
import argparse, os, pickle, sys, tempfile, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C
import rel_reactivity as RV

TMP = tempfile.gettempdir()


def work(chunk):
    out = []
    for ri, subs, cands in chunk:
        s = RV.counts_of(subs).tolist()
        out.append((ri, s, [(k, RV.counts_of(k.split('.')).tolist()) for k in cands]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dump', required=True); ap.add_argument('--out', required=True); ap.add_argument('--workers', type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    d = pickle.load(open(os.path.join(TMP, a.dump), 'rb'))
    rows = C.read_rows()['rows']
    frs = RV.fragment_list()
    old = d['fragments']
    ids = [f[0] for f in frs]
    assert ids[:len(old)] == old, 'the dump\'s fragments are not the first ones of the sheet'
    new = frs[len(old):]
    print(f'{len(old)} fragments in the dump, {len(new)} new ones to count ({time.time() - t0:.0f} s)', flush=True)
    chunks, cur = [], []
    for ri, cands in d['res'].items():
        with_counts = [k for k, c in cands.items() if 'counts' in c]
        cur.append((ri, rows[ri][1], with_counts))
        if len(cur) >= 400: chunks.append(cur); cur = []
    if cur: chunks.append(cur)
    with Pool(a.workers or max(1, (os.cpu_count() or 2) - 1), initializer=RV._init, initargs=([f[2] for f in new],)) as pool:
        for k, res in enumerate(pool.imap_unordered(work, chunks), 1):
            for ri, s_new, cs in res:
                info = d['info'][ri]
                if info.get('s_counts'): info['s_counts'] = list(info['s_counts']) + s_new
                for key, c_new in cs: d['res'][ri][key]['counts'] = list(d['res'][ri][key]['counts']) + c_new
            if k % 10 == 0: print(f'  {k * 400:,} of {len(d["res"]):,} reactions ({time.time() - t0:.0f} s)', flush=True)
    d['fragments'] = ids
    pickle.dump(d, open(os.path.join(TMP, a.out), 'wb'), protocol=4)
    print(f'wrote {a.out} with {len(ids)} fragments ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
