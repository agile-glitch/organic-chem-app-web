"""Turn the sub-fragment scan (tools/subfragment_scan.py) into a short list of sub-fragments worth naming.

    python tools/subfragment_select.py --top 500            # reads tools/subfragment_candidates.json, writes tools/subfragment_selected.json

1. Each candidate (a parent fragment + a feature of one atom of its match) becomes a SMARTS: the parent with the atom constraint added, or with one more
   atom bonded to that atom (an atom outside the match). The SMARTS is run on the reactions the scan used and must find the same support.
2. A candidate that matches (almost) the same molecules as ANOTHER fragment of the sheet is dropped: that group already has a name.
3. Greedy choice per parent fragment: the candidate that adds the most held-out log-likelihood over what is already chosen (consumption rates learned
   on the TRAIN reactions, scored on the VALIDATION reactions, per reagent role signature), until the next one adds less than --min-gain nats.
"""
import argparse, collections, json, os, sys, tempfile, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C
import rel_reactivity as RV
import subfragment_scan as SS

TMP = tempfile.gettempdir()
_W = {}


def make_smarts(parent, positions, kind, feature):
    """parent SMARTS + one atom's constraint (kind 'own') or one extra bonded atom (kind 'nbr'); the atom is the first of the symmetric positions"""
    from rdkit import Chem
    q = Chem.MolFromSmarts(parent)
    k = positions[0]
    rw = Chem.RWMol(q)
    if kind == 'own':
        prim = dict(SS.OWN)[feature]
        inner = rw.GetAtomWithIdx(k).GetSmarts()
        body = inner[1:-1] if inner.startswith('[') else inner
        rw.ReplaceAtom(k, Chem.AtomFromSmarts('[' + body + ';' + prim + ']'))
    else:
        new = rw.AddAtom(Chem.AtomFromSmarts(dict(SS.NBR)[feature]))
        rw.AddBond(k, new, Chem.BondType.SINGLE)
        rw.ReplaceBond(rw.GetBondBetweenAtoms(k, new).GetIdx(), Chem.BondFromSmarts('~'))
    s = Chem.MolToSmarts(rw)
    return s if Chem.MolFromSmarts(s) is not None else None


def _init(smarts):
    from rdkit import RDLogger, Chem
    RDLogger.DisableLog('rdApp.*')
    _W['q'] = [Chem.MolFromSmarts(s) for s in smarts]


def work(items):
    """items = [(reaction index, substrate smiles, [candidate numbers to test])] -> [(reaction index, candidate number)] that match"""
    from rdkit import Chem
    out = []
    for ri, subs, cands in items:
        mols = [Chem.MolFromSmiles(s) for s in subs]
        for c in cands:
            q = _W['q'][c]
            if any(m is not None and m.HasSubstructMatch(q) for m in mols): out.append((ri, c))
    return out


def ll(c, n, p):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return c * np.log(p) + (n - c) * np.log(1 - p)


def cell_gain(code, extra, sg, y, tr):
    """held-out (and train) log-likelihood a further binary feature adds to the grouping code (x reagent signature): rates from TRAIN, shrunk to the coarser cell"""
    fine = (code * 2 + extra) * 4096 + sg
    coarse = code * 4096 + sg
    uf, fi = np.unique(fine, return_inverse=True)
    uc, ci = np.unique(coarse, return_inverse=True)
    n_f = np.bincount(fi[tr], minlength=len(uf)).astype(np.float64); c_f = np.bincount(fi[tr], weights=y[tr], minlength=len(uf))
    n_c = np.bincount(ci[tr], minlength=len(uc)).astype(np.float64); c_c = np.bincount(ci[tr], weights=y[tr], minlength=len(uc))
    overall = (y[tr].sum() + 1) / (tr.sum() + 2)
    p_c = (c_c + 4 * overall) / (n_c + 4)
    # the coarse rate belonging to each fine cell
    fine_to_coarse = np.zeros(len(uf), np.int64); fine_to_coarse[fi] = ci
    p_f = (c_f + 6 * p_c[fine_to_coarse]) / (n_f + 6)
    va = ~tr
    n_fv = np.bincount(fi[va], minlength=len(uf)).astype(np.float64); c_fv = np.bincount(fi[va], weights=y[va], minlength=len(uf))
    n_cv = np.bincount(ci[va], minlength=len(uc)).astype(np.float64); c_cv = np.bincount(ci[va], weights=y[va], minlength=len(uc))
    g_va = float((ll(c_fv, n_fv, p_f) - ll(c_fv, n_fv, p_c[fine_to_coarse])).sum())
    g_tr = float((ll(c_f, n_f, p_f) - ll(c_f, n_f, p_c[fine_to_coarse])).sum())
    return g_va, g_tr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--top', type=int, default=500); ap.add_argument('--n', type=int, default=200000); ap.add_argument('--held', type=int, default=60000); ap.add_argument('--seed', type=int, default=13)
    ap.add_argument('--min-gain', type=float, default=30.0, help='held-out nats a further sub-fragment of the same parent must add')
    ap.add_argument('--min-train-gain', type=float, default=80.0); ap.add_argument('--max-per-parent', type=int, default=4); ap.add_argument('--dup', type=float, default=0.85)
    ap.add_argument('--candidates', default=os.path.join(HERE, 'subfragment_candidates.json')); ap.add_argument('--out', default=os.path.join(HERE, 'subfragment_selected.json'))
    ap.add_argument('--workers', type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    J = json.load(open(a.candidates, encoding='utf-8'))
    cands = [r for r in J['candidates'] if r['gain_held_nats'] > 0 and r['gain_train_nats'] >= a.min_train_gain][:a.top]
    for r in cands: r['sub_smarts'] = make_smarts(r['smarts'], r['positions'], r['kind'], r['feature'])
    cands = [r for r in cands if r['sub_smarts']]
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
    by_parent = collections.defaultdict(list)
    for ci, r in enumerate(cands): by_parent[fidx[r['fragment']]].append(ci)
    items = []
    for i in allr:
        pr = np.flatnonzero(S[i] >= 1)
        cl = [c for f in pr for c in by_parent.get(int(f), ())]
        if cl: items.append((int(i), rows[i][1], cl))
    chunks = [items[k:k + 3000] for k in range(0, len(items), 3000)]
    present = np.zeros((len(allr), len(cands)), bool)
    print(f'{len(cands)} candidates; matching their SMARTS in {len(items):,} reactions ({time.time() - t0:.0f} s)', flush=True)
    with Pool(a.workers or max(1, (os.cpu_count() or 2) - 1), initializer=_init, initargs=([r['sub_smarts'] for r in cands],)) as pool:
        for k, res in enumerate(pool.imap_unordered(work, chunks), 1):
            for ri, c in res: present[pos_in[ri], c] = True
            if k % 10 == 0: print(f'  {min(len(items), k * 3000):,} of {len(items):,} ({time.time() - t0:.0f} s)', flush=True)
    for ci, r in enumerate(cands): r['support_check'] = int(present[is_tr, ci].sum())
    good = [ci for ci, r in enumerate(cands) if abs(r['support_check'] - r['support']) <= 0.03 * r['support'] + 20]
    print(f'the SMARTS find the support the scan counted (within 3 %) for {len(good)} of {len(cands)} candidates', flush=True)
    cands = [cands[ci] for ci in good]; present = present[:, good]
    pres_frag = S[allr] >= 1
    # candidates that are the same molecules as another fragment of the sheet
    for ci, r in enumerate(cands):
        col = present[:, ci]; best = (0.0, None)
        for g in range(pres_frag.shape[1]):
            if g == fidx[r['fragment']]: continue
            inter = int((col & pres_frag[:, g]).sum())
            if inter:
                j = inter / (int(col.sum()) + int(pres_frag[:, g].sum()) - inter)
                if j > best[0]: best = (j, g)
        r['same_as'] = {'fragment': frs[best[1]][0], 'name': frs[best[1]][1], 'jaccard': round(best[0], 2)} if best[1] is not None and best[0] >= 0.5 else None
        r['duplicate'] = bool(best[0] >= a.dup)
    # greedy choice per parent: the candidate that adds the most held-out nats over those already chosen
    selected = []
    parents = collections.defaultdict(list)
    for ci, r in enumerate(cands): parents[r['fragment']].append(ci)
    for fid, cis in sorted(parents.items()):
        f = fidx[fid]
        rows_f = np.flatnonzero(pres_frag[:, f])
        if len(rows_f) < 800: continue
        y = (S[allr[rows_f], f] > P[allr[rows_f], f]).astype(np.float64); sg = sigb[allr[rows_f]]; trm = is_tr[rows_f]
        code = np.zeros(len(rows_f), np.int64); picked = []
        while len(picked) < a.max_per_parent:
            best = None
            for ci in cis:
                if ci in picked or cands[ci]['duplicate']: continue
                g_va, g_tr = cell_gain(code, present[rows_f, ci].astype(np.int64), sg, y, trm)
                if best is None or g_va > best[0]: best = (g_va, g_tr, ci)
            if best is None or best[0] < a.min_gain or best[1] < a.min_train_gain: break
            picked.append(best[2])
            cands[best[2]].update({'added_held_nats': round(best[0], 1), 'added_train_nats': round(best[1], 1), 'rank_in_parent': len(picked)})
            selected.append(cands[best[2]])
            code = code * 2 + present[rows_f, best[2]].astype(np.int64)
    # a few example molecules per selected sub-fragment (from the validation reactions)
    from rdkit import Chem
    for r in selected:
        q = Chem.MolFromSmarts(r['sub_smarts']); ex = []
        ci = cands.index(r)
        for pi in np.flatnonzero(present[:, ci])[:40]:
            for s in rows[int(allr[pi])][1]:
                m = Chem.MolFromSmiles(s)
                if m is not None and m.HasSubstructMatch(q) and len(s) < 70: ex.append(s); break
            if len(ex) >= 3: break
        r['examples'] = ex
    selected.sort(key=lambda r: -r['added_held_nats'])
    json.dump({'when': time.strftime('%Y-%m-%d'), 'selected': selected, 'all': cands}, open(a.out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print(f'\n{len(selected)} sub-fragments chosen over {len({r["fragment"] for r in selected})} parents ({time.time() - t0:.0f} s)')
    for r in selected[:80]:
        print(f"  +{r['added_held_nats']:7.1f} nats  {r['fragment']} {r['name'][:30]:30s} {r['kind']} {r['feature']:22s} support {r['support']:7,d}  rate {r['rate_with']:.2f} vs {r['rate_without']:.2f}"
              + (f"  [same as {r['same_as']['fragment']} {r['same_as']['name'][:20]}, J {r['same_as']['jaccard']}]" if r['same_as'] else ''))


if __name__ == '__main__':
    main()
