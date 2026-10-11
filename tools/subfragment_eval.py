"""Measure the sub-fragments of tools/new_subfragments.py against the recorded reactions.

    python tools/subfragment_eval.py                         # about 2 minutes; writes tools/subfragment_scored.json

For each sub-fragment:
  checks      the SMARTS compiles, matches every pos molecule and none of the neg ones
  support     reactions (of the sample) whose starting materials have the parent AND match it
  rates       how often the parent group is consumed (its count drops from the starting materials to the product) with / without the sub-fragment
  gain        log-likelihood (nats) the split adds to "is the parent consumed?" given the reagent role signature, rates learned on the TRAIN reactions and
              scored on the VALIDATION reactions (reactions the rates never saw)
  overlap     the other fragment of the sheet that matches the most of the same molecules (Jaccard); a sub-fragment that IS another fragment is not needed
"""
import argparse, collections, json, os, sys, tempfile, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C
import rel_reactivity as RV
import subfragment_scan as SS
import subfragment_select as SEL
import new_subfragments as NS

TMP = tempfile.gettempdir()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=200000); ap.add_argument('--held', type=int, default=60000); ap.add_argument('--seed', type=int, default=13)
    ap.add_argument('--workers', type=int, default=0); ap.add_argument('--out', default=os.path.join(HERE, 'subfragment_scored.json'))
    a = ap.parse_args()
    t0 = time.time()
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
    defs = NS.SUBFRAGMENTS
    D = C.read_rows(); rows = D['rows']
    z = np.load(os.path.join(TMP, 'organic_chem_relfragcounts_new.npz'))
    S, P = z['S'], z['P']
    frs = RV.fragment_list(); fidx = {f[0]: i for i, f in enumerate(frs)}
    roles = RV.reagent_roles(); cls = RV.classes_of(rows, roles)
    sigs = np.array([('none' if r[1] == '' else ('?' if r[1] is None else r[1])) for r in cls], dtype=object)
    split = np.array(C.split_codes(rows))
    freq = collections.Counter(sigs[split == C.TRAIN]); keep = {s: i for i, (s, _) in enumerate(freq.most_common(SS.NSIG))}
    sigb = np.array([keep.get(s, SS.NSIG) for s in sigs], np.int64)
    rng = np.random.default_rng(a.seed)
    tr = np.sort(rng.permutation(np.flatnonzero(split == C.TRAIN))[:a.n]); va = np.sort(rng.permutation(np.flatnonzero(split == C.VAL))[:a.held])
    allr = np.concatenate([tr, va]); is_tr = np.r_[np.ones(len(tr), bool), np.zeros(len(va), bool)]
    pos_in = {int(r): i for i, r in enumerate(allr)}
    out = []
    for d in defs:
        q = Chem.MolFromSmarts(d['smarts'])
        row = {'slug': d['slug'], 'name': d['name'], 'parent': d['parent'], 'smarts': d['smarts']}
        if q is None: row['error'] = 'SMARTS does not compile'; out.append(row); print('!!', d['slug'], row['error']); continue
        bad = []
        for s in d.get('pos', []):
            m = Chem.MolFromSmiles(s)
            if m is None or not m.HasSubstructMatch(q): bad.append('misses ' + s)
        for s in d.get('neg', []):
            m = Chem.MolFromSmiles(s)
            if m is not None and m.HasSubstructMatch(q): bad.append('wrongly matches ' + s)
        row['check_problems'] = bad
        out.append(row)
    ok = [i for i, r in enumerate(out) if 'error' not in r]
    parents = {i: [fidx[p] for p in defs[i]['parent']] for i in ok}
    by_parent = collections.defaultdict(list)
    for i in ok:
        for p in parents[i]: by_parent[p].append(i)
    items = []
    for r in allr:
        pr = np.flatnonzero(S[r] >= 1)
        cl = sorted({c for f in pr for c in by_parent.get(int(f), ())})
        if cl: items.append((int(r), rows[r][1], cl))
    chunks = [items[k:k + 3000] for k in range(0, len(items), 3000)]
    present = np.zeros((len(allr), len(defs)), bool)
    print(f'{len(ok)} sub-fragments; matching in {len(items):,} reactions ({time.time() - t0:.0f} s)', flush=True)
    with Pool(a.workers or max(1, (os.cpu_count() or 2) - 1), initializer=SEL._init, initargs=([d['smarts'] for d in defs],)) as pool:
        for k, res in enumerate(pool.imap_unordered(SEL.work, chunks), 1):
            for ri, c in res: present[pos_in[ri], c] = True
    pres_frag = S[allr] >= 1
    Sa, Pa = S[allr], P[allr]
    for i in ok:
        pcols = parents[i]
        parent_present = pres_frag[:, pcols].any(1)
        y_all = (Sa[:, pcols] > Pa[:, pcols]).any(1).astype(np.float64)
        rows_f = np.flatnonzero(parent_present)
        z_ = present[rows_f, i].astype(np.int64)
        y = y_all[rows_f]; sg = sigb[allr[rows_f]]; trm = is_tr[rows_f]
        n1 = int(z_[trm].sum()); n0 = int((1 - z_)[trm].sum())
        r = out[i]
        r['support_train'] = n1; r['without_train'] = n0
        r['support_held'] = int(z_[~trm].sum())
        r['rate_with'] = round(float(y[trm][z_[trm] == 1].mean()), 3) if n1 else None
        r['rate_without'] = round(float(y[trm][z_[trm] == 0].mean()), 3) if n0 else None
        if n1 >= 100 and n0 >= 100:
            g_va, g_tr = SEL.cell_gain(np.zeros(len(rows_f), np.int64), z_, sg, y, trm)
            r['gain_held_nats'] = round(g_va, 1); r['gain_train_nats'] = round(g_tr, 1)
        col = present[:, i]; best = (0.0, None)
        for g in range(pres_frag.shape[1]):
            if g in pcols: continue
            inter = int((col & pres_frag[:, g]).sum())
            if inter:
                j = inter / (int(col.sum()) + int(pres_frag[:, g].sum()) - inter)
                if j > best[0]: best = (j, g)
        r['overlap'] = {'fragment': frs[best[1]][0], 'name': frs[best[1]][1][:40], 'jaccard': round(best[0], 2)} if best[1] is not None else None
        # how much of the parent's molecules it covers
        r['share_of_parent'] = round(float(z_.mean()), 3)
    json.dump(out, open(a.out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print(f'\n{"slug":24s} {"parent":8s} {"with":>7s} {"rate":>5s} {"without":>8s} {"rate":>5s} {"held nats":>9s} {"share":>6s}  overlap')
    for r in sorted([r for r in out if 'error' not in r], key=lambda r: -(r.get('gain_held_nats') or -1e9)):
        o = r.get('overlap') or {}
        print(f"{r['slug']:24s} {','.join(r['parent'])[:8]:8s} {r.get('support_train', 0):7,d} {r.get('rate_with') if r.get('rate_with') is not None else float('nan'):5.2f} {r.get('without_train', 0):8,d} "
              f"{r.get('rate_without') if r.get('rate_without') is not None else float('nan'):5.2f} {r.get('gain_held_nats', float('nan')):9.1f} {r.get('share_of_parent', 0):6.1%}  "
              f"{o.get('fragment', '')} {o.get('name', '')} J{o.get('jaccard', '')}" + ('  CHECK: ' + '; '.join(r['check_problems']) if r['check_problems'] else ''))
    print(f'({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
