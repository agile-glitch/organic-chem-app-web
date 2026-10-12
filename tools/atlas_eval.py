"""Fragment Atlas, step 2: does knowing the SURROUNDINGS of a group (and its sterics, charge, ring angle and pKa) help to say
whether a reaction will transform it? Measured on reactions held out by scaffold.

    python tools/atlas_eval.py <instances.npz> [--out tools/atlas_eval.json]

The task: for a functional-group instance in a molecule that reacts, given the reagents, will the reaction transform that group
(y = 1) or leave it alone (y = 0)? That is the question behind every interference ("the ester survives NaBH4", "LiAlH4 takes the
nitro group too") and behind choosing a reagent that reaches one group and spares another.

Counts come from TRAIN reactions; the small combining model is fitted on VALIDATION; every figure is on TEST (scaffold split,
tools/atlas_mine.py). Models, each adding to the one before it:

  M0  the group alone                       P(y | group)                                  what a bare fragment list knows
  M1  + the kind of reagent                 P(y | group, reagent class)                    the earlier reactivity table
  M2  + the group's surroundings            P(y | group, class, environment radius 1..3)   contextual fragments
  M3  + the named reagents                  P(y | group, reagent part [, environment])
  M4  + steric / charge / ring-angle / pKa  the continuous descriptors of the reacting atom
  GB  gradient boosting on the same inputs  an upper look at what a non-linear combiner could add (not shipped)

Each count is smoothed down a chain (coarse -> fine): p_k = (hits_k + m p_{k-1}) / (n_k + m), so a rare environment borrows from
the coarser one and a well-measured environment speaks for itself. The combiner is a logistic regression over the logits of those
chain estimates (+ the amount of evidence behind each) and the descriptors.
"""
import argparse, json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import atlas_keys as K


class Table:
    """(sorted unique key, n, hits) from the training instances, keeping the keys met at least min_n times; a key never kept gives (0, 0)"""
    def __init__(self, keys, y, min_n=1):
        self.k, inv = np.unique(keys, return_inverse=True)
        self.n = np.bincount(inv).astype(np.float64)
        self.h = np.bincount(inv, weights=y).astype(np.float64)
        keep = self.n >= min_n
        self.k, self.n, self.h = self.k[keep], self.n[keep], self.h[keep]

    def get(self, keys):
        i = np.searchsorted(self.k, keys)
        i[i >= len(self.k)] = 0
        ok = self.k[i] == keys
        return np.where(ok, self.n[i], 0.0), np.where(ok, self.h[i], 0.0)


def shrink(prev, n, h, m):
    return (h + m * prev) / (n + m)


def logit(p):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


def logloss(y, p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return float(-np.mean(y * np.log(p) + (1 - y) * np.log(1 - p)))


def auc(y, s):
    order = np.argsort(s, kind='mergesort'); r = np.empty(len(s)); r[order] = np.arange(1, len(s) + 1)
    # average ranks for ties
    _, inv, cnt = np.unique(s, return_inverse=True, return_counts=True)
    if len(cnt) < len(s):
        sums = np.bincount(inv, weights=r); r = (sums / cnt)[inv]
    n1 = y.sum(); n0 = len(y) - n1
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)) if n1 and n0 else float('nan')


def fit_logistic(X, y, l2=1.0, iters=30):
    """L2 logistic regression by Newton / IRLS. X has no intercept column (one is added). Returns coef (intercept first), mean, std."""
    mu, sd = X.mean(0), X.std(0) + 1e-9
    Z = np.hstack([np.ones((len(X), 1)), (X - mu) / sd])
    w = np.zeros(Z.shape[1])
    reg = np.eye(Z.shape[1]) * l2; reg[0, 0] = 0
    for _ in range(iters):
        p = 1 / (1 + np.exp(-np.clip(Z @ w, -30, 30)))
        g = Z.T @ (p - y) + reg @ w
        H = (Z * (p * (1 - p))[:, None]).T @ Z + reg + 1e-9 * np.eye(Z.shape[1])
        step = np.linalg.solve(H, g)
        w -= step
        if np.abs(step).max() < 1e-7: break
    return w, mu, sd


def predict_logistic(model, X):
    w, mu, sd = model
    z = w[0] + ((X - mu) / sd) @ w[1:]
    return 1 / (1 + np.exp(-np.clip(z, -30, 30)))


def within_fragment_auc(y, s, frag, min_n=300):
    """AUC inside each fragment (so it measures what the context adds beyond 'which group it is'), averaged over fragments weighted by size"""
    tot, wsum, per = 0.0, 0.0, {}
    for f in np.unique(frag):
        m = frag == f
        if m.sum() < min_n: continue
        a = auc(y[m], s[m])
        if a == a:
            per[int(f)] = (a, int(m.sum()))
            tot += a * m.sum(); wsum += m.sum()
    return (tot / wsum if wsum else float('nan')), per


def acid_pka_by_smiles():
    import learn_reactions as L
    from rdkit import Chem
    out = {}
    for r in L.read_sheet(os.path.join(L.APP, 'reaction_rules.xlsx'), 'Acids_Bases'):
        if str(r.get('Role')).strip() != 'acid': continue
        try: pk = float(r['pKa (acid: its own; base: of its conjugate acid)'])
        except Exception:
            ks = [k for k in r if str(k).startswith('pKa')]
            try: pk = float(r[ks[0]])
            except Exception: continue
        m = Chem.MolFromSmiles(str(r['SMILES (as typed in the reagent box)'])) if 'SMILES (as typed in the reagent box)' in r else None
        if m is None:
            ks = [k for k in r if str(k).startswith('SMILES')]
            m = Chem.MolFromSmiles(str(r[ks[0]]))
        if m is not None: out[Chem.MolToSmiles(m)] = pk
    return out


def build_features(d, P, train_mask, rows, m_chain=K.M_CHAIN, min_n=K.MIN_N):
    """chain estimates + evidence + descriptors for the `rows` (boolean mask), counts from train_mask. Returns (names, groups, X, row indices)"""
    return build_features_with_base(d, P, train_mask, rows, float(d['y'].astype(np.float64)[train_mask].mean()), m_chain, min_n)


def build_features_with_base(d, P, train_mask, rows, base, m_chain=K.M_CHAIN, min_n=K.MIN_N):
    y = d['y'].astype(np.float64)
    tr_idx = np.where(train_mask)[0]
    ktr = K.family_keys(P, tr_idx)
    ytr = y[tr_idx]
    tabs = {}
    for name in ('F', 'R0', 'R1', 'R2', 'R3', 'U1', 'U2', 'U3'):
        tabs[name] = Table(ktr[name], ytr, min_n[name])
    for lvl in ('P0', 'P1'):                                  # a part is a part: the three slots share one table
        keys = np.concatenate([ktr[(lvl, j)][P['has'][j][tr_idx]] for j in range(3)])
        yy = np.concatenate([ytr[P['has'][j][tr_idx]] for j in range(3)])
        tabs[lvl] = Table(keys, yy, min_n[lvl])
    ridx = np.where(rows)[0]
    kr = K.family_keys(P, ridx)
    g = base
    cols, names = [], []

    def add(nm, v):
        cols.append(np.asarray(v, dtype=np.float64)); names.append(nm)

    nF, hF = tabs['F'].get(kr['F']); pF = shrink(g, nF, hF, 2.0)
    n, h = tabs['R0'].get(kr['R0']); pR0 = shrink(pF, n, h, m_chain); nR0 = n
    chain = {'F': (pF, nF), 'R0': (pR0, nR0)}
    prev = pR0
    for lv in (1, 2, 3):
        n, h = tabs[f'R{lv}'].get(kr[f'R{lv}']); prev = shrink(prev, n, h, m_chain); chain[f'R{lv}'] = (prev, n)
    prev = pF
    for lv in (1, 2, 3):
        n, h = tabs[f'U{lv}'].get(kr[f'U{lv}']); prev = shrink(prev, n, h, m_chain); chain[f'U{lv}'] = (prev, n)
    # named reagents: the part with the most evidence
    best_n = np.full(len(ridx), -1.0); pP0 = pR0.copy(); pP1 = pR0.copy(); nP0 = np.zeros(len(ridx)); nP1 = np.zeros(len(ridx))
    for j in range(3):
        has = P['has'][j][ridx]
        n0, h0 = tabs['P0'].get(kr[('P0', j)]); n1, h1 = tabs['P1'].get(kr[('P1', j)])
        p0 = shrink(pF, n0, h0, m_chain); p1 = shrink(p0, n1, h1, m_chain)
        better = has & (n0 > best_n)
        pP0 = np.where(better, p0, pP0); pP1 = np.where(better, p1, pP1); nP0 = np.where(better, n0, nP0); nP1 = np.where(better, n1, nP1)
        best_n = np.where(better, n0, best_n)
    chain['P0'] = (pP0, nP0); chain['P1'] = (pP1, nP1)
    groups = {'M0': ['F'], 'M1': ['R0'], 'M2': ['R1', 'R2', 'R3', 'U1', 'U2', 'U3'], 'M3': ['P0', 'P1']}
    colgroup = []
    for gname, members in groups.items():
        for nm in members:
            p, nn = chain[nm]
            add('logit_' + nm, logit(p)); colgroup.append(gname)
            if nm not in ('F',):
                add('logn_' + nm, np.log1p(nn)); colgroup.append(gname)
    # descriptors
    st = d['steric'][ridx].astype(np.float64)
    q = d['charge'][ridx].astype(np.float64)
    ring = d['ring'][ridx]; strain = d['strain'][ridx].astype(np.float64); hyb = d['hyb'][ridx]
    nat = d['natoms'][ridx].astype(np.float64)
    desc = [('steric', np.clip(st, 0, 9)), ('steric_sq', np.clip(st, 0, 9) ** 2 / 9),
            ('charge', np.nan_to_num(q, nan=0.0)), ('charge_missing', np.isnan(q).astype(float)),
            ('abs_charge', np.abs(np.nan_to_num(q, nan=0.0))),
            ('ring3', (ring == 3).astype(float)), ('ring4', (ring == 4).astype(float)), ('ring5', (ring == 5).astype(float)),
            ('ring6', (ring == 6).astype(float)), ('ring7p', (ring >= 7).astype(float)), ('in_ring', (ring > 0).astype(float)),
            ('strain', strain / 30.0), ('sp', (hyb == 1).astype(float)), ('sp2', (hyb == 2).astype(float)), ('match_atoms', np.clip(nat, 0, 12) / 12)]
    # pKa: the reagent's strength against the group's acidity / basicity
    import atlas_core as A
    fp = A.fragment_pka(); frag_ids = list(d['frag_ids'])
    f_acid = np.array([np.nan if fp.get(i, (None, None))[0] is None else fp[i][0] for i in frag_ids], np.float64)
    f_base = np.array([np.nan if fp.get(i, (None, None))[1] is None else fp[i][1] for i in frag_ids], np.float64)
    acid_tab = acid_pka_by_smiles()
    parts = list(d['parts']); roles = list(d['part_role']); pkas = d['part_pka'].astype(np.float64)
    part_base = np.full(len(parts) + 1, np.nan); part_acid = np.full(len(parts) + 1, np.nan)
    for i, (s, r, pk) in enumerate(zip(parts, roles, pkas), start=1):
        if r == 'base' and pk == pk: part_base[i] = pk
        if r == 'acid' and s in acid_tab: part_acid[i] = acid_tab[s]
    pid = [d['p1'], d['p2'], d['p3']]                                     # part ids (0 = none) into the vocabulary
    pb = np.full(len(ridx), np.nan); pa = np.full(len(ridx), np.nan)
    for j in range(3):
        v = part_base[pid[j][ridx]]; pb = np.fmax(pb, v)                 # the strongest base among the parts
        v = part_acid[pid[j][ridx]]; pa = np.fmin(pa, v) if j else v      # the strongest acid
    fi = d['frag'][ridx]
    dpb = pb - f_acid[fi]; dpa = f_base[fi] - pa
    desc += [('has_base', (~np.isnan(pb)).astype(float)), ('dpka_base', np.clip(np.nan_to_num(dpb, nan=0.0), -25, 35) / 10),
             ('has_acid', (~np.isnan(pa)).astype(float)), ('dpka_acid', np.clip(np.nan_to_num(dpa, nan=0.0), -25, 35) / 10)]
    for nm, v in desc:
        add(nm, v); colgroup.append('M4')
    X = np.column_stack(cols)
    return names, colgroup, X, ridx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('npz'); ap.add_argument('--out', default=os.path.join(HERE, 'atlas_eval.json'))
    ap.add_argument('--m', type=float, default=5.0)
    ap.add_argument('--target', choices=('any', 'bystander'), default='any',
                    help="any: is the group transformed (all groups of reacting molecules); bystander: the reaction happens elsewhere in the molecule - is this group ALSO transformed "
                         "(only groups that touch no site, or sit on one of two or more; tools/atlas_mine.py)")
    ap.add_argument('--split', choices=('scaffold', 'molecule'), default='scaffold', help='scaffold: related molecules stay on one side; molecule: a held-out molecule never occurs in a training reaction')
    a = ap.parse_args()
    t0 = time.time()
    d = dict(np.load(a.npz, allow_pickle=True))
    split = d['split_mol'] if a.split == 'molecule' else d['split']
    d['split'] = split
    el = np.ones(len(split), bool)
    if a.target == 'bystander':
        el = d['elig'] == 1
        d['y'] = d['yb']
    y = d['y'].astype(np.float64)
    print(f'{len(y):,} instances ({a.target}: {int(el.sum()):,} used); train {int(((split == 2) & el).sum()):,}, validation {int(((split == 1) & el).sum()):,}, test {int(((split == 0) & el).sum()):,}', flush=True)
    train = (split == 2) & el
    vmask, tmask = (split == 1) & el, (split == 0) & el
    res = {'target': a.target, 'split': a.split, 'instances': int(len(y)), 'train': int(train.sum()), 'validation': int(vmask.sum()), 'test': int(tmask.sum()), 'm': a.m}
    P = K.prepare(d)
    names, grp, Xv, iv = build_features(d, P, train, vmask, a.m)
    _, _, Xt, it = build_features(d, P, train, tmask, a.m)
    yv, yt = y[iv], y[it]
    ft = d['frag'][it]
    print(f'features built in {time.time() - t0:.0f}s: {len(names)} columns', flush=True)
    order = ['M0', 'M1', 'M2', 'M3', 'M4']
    out = {}
    preds = {}
    base_rate = float(y[train].mean())
    out['constant'] = {'logloss': logloss(yt, np.full(len(yt), base_rate)), 'auc': 0.5}
    for k, upto in enumerate(order):
        cols = [i for i, g in enumerate(grp) if g in order[:k + 1]]
        mdl = fit_logistic(Xv[:, cols], yv, l2=2.0)
        p = predict_logistic(mdl, Xt[:, cols]); preds[upto] = p
        wf, _ = within_fragment_auc(yt, p, ft)
        out[upto] = {'logloss': logloss(yt, p), 'auc': auc(yt, p), 'within_fragment_auc': wf, 'columns': len(cols),
                     'brier': float(np.mean((p - yt) ** 2))}
        print(f'{upto}: logloss {out[upto]["logloss"]:.4f}  AUC {out[upto]["auc"]:.4f}  within-fragment AUC {wf:.4f}  ({len(cols)} columns)', flush=True)
    # the descriptors alone on top of M0 (do they help without the environment keys?)
    cols = [i for i, g in enumerate(grp) if g in ('M0', 'M1', 'M4')]
    mdl = fit_logistic(Xv[:, cols], yv, l2=2.0); p = predict_logistic(mdl, Xt[:, cols]); wf, _ = within_fragment_auc(yt, p, ft)
    out['M1+M4 (no environment keys)'] = {'logloss': logloss(yt, p), 'auc': auc(yt, p), 'within_fragment_auc': wf}
    print(f'M1+M4: logloss {out["M1+M4 (no environment keys)"]["logloss"]:.4f}  AUC {out["M1+M4 (no environment keys)"]["auc"]:.4f}  within {wf:.4f}', flush=True)
    # a non-linear combiner on everything (headroom)
    try:
        from sklearn.ensemble import HistGradientBoostingClassifier
        gb = HistGradientBoostingClassifier(max_iter=250, learning_rate=0.08, max_leaf_nodes=31, l2_regularization=1.0, random_state=0)
        gb.fit(Xv, yv); p = gb.predict_proba(Xt)[:, 1]; preds['GB'] = p; wf, _ = within_fragment_auc(yt, p, ft)
        out['GB'] = {'logloss': logloss(yt, p), 'auc': auc(yt, p), 'within_fragment_auc': wf}
        print(f'GB: logloss {out["GB"]["logloss"]:.4f}  AUC {out["GB"]["auc"]:.4f}  within {wf:.4f}', flush=True)
    except Exception as ex:
        print('GB skipped:', ex)
    res['results'] = out
    # how each model does where the surroundings were NEVER seen in training (the 'unproven reaction' case): the descriptors have to carry it there
    ix = {nm: i for i, nm in enumerate(names)}
    seen = lambda lv: (Xt[:, ix[f'logn_R{lv}']] > 0) | (Xt[:, ix[f'logn_U{lv}']] > 0)
    subsets = {'all': np.ones(len(yt), bool), 'seen at radius 3': seen(3), 'seen at radius 2, not 3': seen(2) & ~seen(3),
               'seen at radius 1, not 2': seen(1) & ~seen(2), 'radius 1 never seen': ~seen(1)}
    sub = {}
    for nm, mk in subsets.items():
        if mk.sum() < 500: continue
        row = {'share': float(mk.mean()), 'n': int(mk.sum()), 'hit_rate': float(yt[mk].mean())}
        for mname, p in preds.items():
            row[mname] = {'logloss': logloss(yt[mk], p[mk]), 'auc': auc(yt[mk], p[mk])}
        sub[nm] = row
        print(f'  {nm:26s} {mk.mean() * 100:5.1f} % of test | ' + ' | '.join(f'{k} {v["logloss"]:.3f}/{v["auc"]:.3f}' for k, v in row.items() if isinstance(v, dict)), flush=True)
    res['by_novelty'] = sub
    # SITE CHOICE: a molecule with several instances of the same group (two esters, four aromatic C-H ...) of which exactly one was transformed.
    # Which one does each model pick? (This is where steric hindrance, partial charge and ring angle are expected to matter: the group, the
    # reagent and the molecule are the same, only the surroundings differ.)
    if 'rx' in d and 'mol' in d and a.target == 'any':
        g = (d['rx'][it].astype(np.uint64) * np.uint64(1000003)) ^ (d['mol'][it].astype(np.uint64) * np.uint64(7919)) ^ (ft.astype(np.uint64) * np.uint64(104729))
        _, gid = np.unique(g, return_inverse=True)
        size = np.bincount(gid); hits = np.bincount(gid, weights=yt)
        ok = (size[gid] >= 2) & (hits[gid] == 1)
        sel = np.where(ok)[0]
        ng = len(np.unique(gid[sel]))
        site = {'groups': int(ng), 'random_top1': float(np.mean(1.0 / size[np.unique(gid[sel])]))}
        # score of each model; ties broken at random, so a model with no information scores like chance
        rng = np.random.default_rng(0)
        order = sel[np.argsort(gid[sel], kind='stable')]
        starts = np.r_[0, np.flatnonzero(np.diff(gid[order])) + 1]
        ends = np.r_[starts[1:], len(order)]
        for mname, p in preds.items():
            jitter = p + rng.random(len(p)) * 1e-9
            right = 0
            for a0, b0 in zip(starts, ends):
                idx = order[a0:b0]
                right += int(yt[idx[np.argmax(jitter[idx])]] == 1)
            site[mname] = right / len(starts)
        # the descriptors on their own, with the group (what they know without any counted context)
        dcols = [i for i, gname in enumerate(grp) if gname == 'M4']
        mD = fit_logistic(np.column_stack([Xv[:, [ix['logit_F']]], Xv[:, dcols]]), yv, l2=2.0)
        pD = predict_logistic(mD, np.column_stack([Xt[:, [ix['logit_F']]], Xt[:, dcols]]))
        jitter = pD + rng.random(len(pD)) * 1e-9
        site['descriptors_only'] = sum(int(yt[order[a0:b0][np.argmax(jitter[order[a0:b0]])]] == 1) for a0, b0 in zip(starts, ends)) / len(starts)
        res['site_choice'] = site
        print('  site choice (%d molecules with several of the same group, one transformed): ' % ng + ' | '.join(f'{k} {v:.3f}' for k, v in site.items() if k != 'groups'), flush=True)
    # the shipped combiner (M4, all columns) fitted on validation
    cols = list(range(len(names)))
    mdl = fit_logistic(Xv, yv, l2=2.0)
    res['model'] = {'names': names, 'coef': [float(v) for v in mdl[0]], 'mean': [float(v) for v in mdl[1]], 'std': [float(v) for v in mdl[2]]}
    # per-fragment gain M4 vs M0 on test
    p4 = predict_logistic(mdl, Xt)
    c0 = [i for i, g in enumerate(grp) if g == 'M0']
    m0 = fit_logistic(Xv[:, c0], yv, l2=2.0); p0 = predict_logistic(m0, Xt[:, c0])
    ids = list(d['frag_ids'])
    rows = []
    for f in np.unique(ft):
        mk = ft == f
        if mk.sum() < 400: continue
        rows.append({'fragment': ids[int(f)], 'n': int(mk.sum()), 'auc_M0': auc(yt[mk], p0[mk]), 'auc_M4': auc(yt[mk], p4[mk]), 'hit_rate': float(yt[mk].mean())})
    res['per_fragment'] = sorted(rows, key=lambda r: -r['n'])
    # calibration of M4 on test (deciles)
    qs = np.quantile(p4, np.linspace(0, 1, 11)); cal = []
    for lo, hi in zip(qs[:-1], qs[1:]):
        mk = (p4 >= lo) & (p4 <= hi)
        if mk.sum(): cal.append({'stated': float(p4[mk].mean()), 'measured': float(yt[mk].mean()), 'n': int(mk.sum())})
    res['calibration'] = cal
    json.dump(res, open(a.out, 'w'), indent=1)
    print('wrote', a.out, f'({time.time() - t0:.0f}s)')


if __name__ == '__main__':
    main()
