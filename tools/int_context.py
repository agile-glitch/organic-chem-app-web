"""Chemoselectivity from the molecule's other groups: does a competitor change whether a group reacts?

    python tools/int_context.py --tag new --out ctx_new            # fit on the training patents, report on the validation patents

tools/rel_reactivity.py learns, per reagent, how often each functional group is consumed when it is present (NaOH takes the ester, NaBH4 the aldehyde
before the ketone ...). That rate does not look at what ELSE is in the molecule. Here each fragment F gets a logistic model of "F is consumed" whose
offset is that reagent rate and whose extra terms are the other fragments present (G): beta[F, G] < 0 means that G being there protects F (G is the
site the reagent prefers, or G quenches/consumes the reagent), beta > 0 that G helps (G activates F, or the reagent is used in excess when G is present).
This is the "competing site" interference of the Interferences sheet, learned instead of written: the pairs that the Mantel–Haenszel scan
(tools/interference_scan.py) found by stratifying on the reagent are what these coefficients put into the prediction.

The offsets are CROSS-FITTED (the reagent rate of a training reaction comes from the other four folds of the training patents), so the coefficients are
not fitted against rates that already contain the reaction they are scored on; the validation/test offsets come from the full training table, as on the page.
"""
import argparse, collections, hashlib, json, os, pickle, sys, tempfile, time

import numpy as np
from scipy import sparse
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C
import rel_reactivity as RV

TMP = tempfile.gettempdir()
K = RV.K_SMOOTH
FOLDS = 5


def fold_of(patent):
    return int(hashlib.md5(('ctx' + str(patent)).encode()).hexdigest(), 16) % FOLDS


def meta_path(): return os.path.join(TMP, 'organic_chem_ctx_meta.pkl')


def load_meta():
    """per reaction: split code, fold, exact reagent id (-1 none), role signature id (-1 none/unknown); small arrays, cached so the 1.5 GB of rows need not be loaded again"""
    if os.path.exists(meta_path()) and os.path.getmtime(meta_path()) > os.path.getmtime(C.CACHE):
        m = pickle.load(open(meta_path(), 'rb'))
        if 'patent' in m and 'parts_ids' in m: return m
    D = C.read_rows(); rows = D['rows']
    roles = RV.reagent_roles()
    split = np.array(C.split_codes(rows), np.int8)
    cache = {}
    fold = np.array([cache.setdefault(r[8], fold_of(r[8])) for r in rows], np.int8)
    cls = RV.classes_of(rows, roles)
    ex_ids, rl_ids = {}, {}
    ex = np.array([ex_ids.setdefault(e, len(ex_ids)) if e else -1 for e, _ in cls], np.int32)
    rl = np.array([rl_ids.setdefault(r, len(rl_ids)) if r is not None else -1 for _, r in cls], np.int32)
    pid = {}
    patent = np.array([pid.setdefault(r[8], len(pid)) for r in rows], np.int32)
    vocab, flat, offs = {}, [], [0]
    for r in rows:
        flat.extend(sorted({vocab.setdefault(x, len(vocab)) for x in str(r[2] or '').split('.') if x}))
        offs.append(len(flat))
    out = {'split': split, 'fold': fold, 'exact': ex, 'role': rl, 'exact_keys': list(ex_ids), 'role_keys': list(rl_ids), 'patent': patent,
           'parts_vocab': list(vocab), 'parts_ids': np.array(flat, np.int32), 'parts_offsets': np.array(offs, np.int64)}
    pickle.dump(out, open(meta_path(), 'wb'), protocol=4)
    del D, rows
    return out


def group_sums(M, gid, n_groups):
    """sum the rows of M (N x F, bool/uint8) by group id (rows with gid < 0 ignored): (n_groups x F) float32"""
    use = gid >= 0
    order = np.argsort(gid[use], kind='stable')
    idx = np.flatnonzero(use)[order]
    g = gid[idx]
    out = np.zeros((n_groups, M.shape[1]), np.float32)
    if not len(idx): return out
    starts = np.r_[0, np.flatnonzero(np.diff(g)) + 1]
    out[g[starts]] = np.add.reduceat(M[idx], starts, axis=0, dtype=np.float32)
    return out


def rate_tables(present, pc, meta, mask, min_exact=30):
    """smoothed consumption-rate tables from the reactions in mask (the chain of tools/rel_reactivity.py Prior.rates): global (F,), role (R x F), and the exact
    reagents that have at least min_exact records in mask (compact ids: tab['emap'][exact id] -> row of ep/ec, -1 = not enough)"""
    gp = present[mask].sum(0).astype(np.float64); gc = pc[mask].sum(0).astype(np.float64)
    g_rate = (gc + 1.0) / (gp + 2.0)
    R = len(meta['role_keys'])
    rid = np.where(mask, meta['role'], -1)
    rp, rc = group_sums(present, rid, R), group_sums(pc, rid, R)
    rn = np.bincount(rid[rid >= 0], minlength=R)
    role_rate = np.where(rn[:, None] > 0, (rc + K * g_rate[None, :]) / (rp + K), g_rate[None, :])
    eid = np.where(mask, meta['exact'], -1)
    en = np.bincount(eid[eid >= 0], minlength=len(meta['exact_keys']))
    keep = np.flatnonzero(en >= min_exact)
    emap = np.full(len(meta['exact_keys']), -1, np.int32); emap[keep] = np.arange(len(keep))
    cid = np.where(eid >= 0, emap[np.maximum(eid, 0)], -1)
    ep, ec = group_sums(present, cid, len(keep)), group_sums(pc, cid, len(keep))
    return {'g': g_rate, 'role': role_rate.astype(np.float64), 'role_has': rn > 0, 'ep': ep, 'ec': ec, 'emap': emap}


def rates_for(tab, meta, idx):
    """the smoothed rate vector (len(idx) x F) of the reactions idx: exact if the reagent has enough records, else its roles, else global"""
    n = len(idx)
    base = np.tile(tab['g'], (n, 1))
    rid = meta['role'][idx]
    has_r = (rid >= 0)
    has_r[has_r] = tab['role_has'][rid[has_r]]
    base[has_r] = tab['role'][rid[has_r]]
    eid = meta['exact'][idx]
    cid = np.where(eid >= 0, tab['emap'][np.maximum(eid, 0)], -1)
    has_e = cid >= 0
    ep, ec = tab['ep'][cid[has_e]], tab['ec'][cid[has_e]]
    base[has_e] = (ec + K * base[has_e]) / (ep + K)
    return np.clip(base, 1e-3, 1 - 1e-3)


def fit_one(Xs, y, off, lam, maxiter=200):
    """logistic regression with an offset: p = sigmoid(off + b + X beta); returns (b, beta)"""
    n, d = Xs.shape
    def f(w):
        z = off + w[0] + Xs @ w[1:]
        p = 1 / (1 + np.exp(-z))
        loss = -(y * np.log(p + 1e-12) + (1 - y) * np.log(1 - p + 1e-12)).sum() + 0.5 * lam * (w[1:] ** 2).sum()
        r = p - y
        g = np.r_[r.sum(), Xs.T @ r + lam * w[1:]]
        return loss, g
    res = minimize(f, np.zeros(d + 1), jac=True, method='L-BFGS-B', options={'maxiter': maxiter})
    return res.x[0], res.x[1:]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tag', default='new'); ap.add_argument('--out', default='ctx_new'); ap.add_argument('--lam', type=float, default=5.0)
    ap.add_argument('--grid', action='store_true', help='choose the penalty per group from 2, 10, 50 (on one fold of the training patents)'); ap.add_argument('--maxiter', type=int, default=150)
    ap.add_argument('--first', type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    meta = load_meta()
    z = np.load(RV.counts_path(a.tag)); S, P = z['S'], z['P']
    frag_ids = [str(x) for x in z['fragments']]
    nF = S.shape[1]
    present = S >= 1
    cons = S > P
    pc = present & cons
    train = meta['split'] == C.TRAIN; val = meta['split'] == C.VAL
    print(f'{S.shape[0]:,} reactions x {nF} fragments; train {int(train.sum()):,} validation {int(val.sum()):,} ({time.time() - t0:.0f} s)', flush=True)

    def entries(tab, idx):
        """(reaction, fragment, offset) for every fragment present in the reactions idx, the offset being the logit of the reagent rate (cross-fitted tables: see the module text)"""
        er, ef, eo = [], [], []
        for a0 in range(0, len(idx), 20000):
            sub = idx[a0:a0 + 20000]
            r = rates_for(tab, meta, sub)
            i, f = np.nonzero(present[sub])
            er.append(sub[i].astype(np.int64)); ef.append(f.astype(np.int16)); eo.append(np.log(r[i, f] / (1 - r[i, f])).astype(np.float32))
        return np.concatenate(er), np.concatenate(ef), np.concatenate(eo)
    tr_parts = []
    for k in range(FOLDS):
        tab = rate_tables(present, pc, meta, train & (meta['fold'] != k))
        tr_parts.append(entries(tab, np.flatnonzero(train & (meta['fold'] == k))))
        print(f'  fold {k}: offsets ({time.time() - t0:.0f} s)', flush=True)
        del tab
    tr_r, tr_f, tr_o = [np.concatenate([p[j] for p in tr_parts]) for j in range(3)]
    del tr_parts
    tab_full = rate_tables(present, pc, meta, train)
    va_r, va_f, va_o = entries(tab_full, np.flatnonzero(val))
    del tab_full
    o_tr = np.argsort(tr_f, kind='stable'); tr_r, tr_f, tr_o = tr_r[o_tr], tr_f[o_tr], tr_o[o_tr]
    o_va = np.argsort(va_f, kind='stable'); va_r, va_f, va_o = va_r[o_va], va_f[o_va], va_o[o_va]
    tr_start = np.searchsorted(tr_f, np.arange(nF + 1)); va_start = np.searchsorted(va_f, np.arange(nF + 1))
    print(f'offsets for {len(tr_r):,} training and {len(va_r):,} validation (reaction, group) pairs ({time.time() - t0:.0f} s)', flush=True)
    pr_i, pr_j = np.nonzero(present)
    Pres = sparse.csr_matrix((np.ones(len(pr_i), np.float32), (pr_i, pr_j)), shape=present.shape)
    del pr_i, pr_j
    mult2 = (S >= 2)
    b = np.zeros(nF); beta = np.zeros((nF, nF)); mult = np.zeros(nF)
    lams = [a.lam] if not a.grid else [2.0, 10.0, 50.0]
    stats = []
    ll = lambda p, y: float(-(y * np.log(p + 1e-12) + (1 - y) * np.log(1 - p + 1e-12)).mean())
    sig = lambda z: 1 / (1 + np.exp(-z))
    for f in range(nF):
        if a.first and f >= a.first: break
        a0, a1 = tr_start[f], tr_start[f + 1]; c0, c1 = va_start[f], va_start[f + 1]
        rows_t, rows_v = tr_r[a0:a1], va_r[c0:c1]
        if len(rows_t) < 200 or len(rows_v) < 20: continue
        keep = np.delete(np.arange(nF), f)
        def design(rows):
            return sparse.hstack([Pres[rows][:, keep], sparse.csr_matrix(mult2[rows, f].astype(np.float32)[:, None])]).tocsr()
        Xt, Xv = design(rows_t), design(rows_v)
        yt = cons[rows_t, f].astype(np.float64); yv = cons[rows_v, f].astype(np.float64)
        ot = tr_o[a0:a1].astype(np.float64); ov = va_o[c0:c1].astype(np.float64)
        # the penalty and whether the context helps at all are chosen on one fold of the TRAINING patents (fold 0), so the validation patents stay untouched
        hold = meta['fold'][rows_t] == 0
        l_hold0 = ll(sig(ot[hold]), yt[hold]) if hold.any() else None
        chosen, best_l = None, l_hold0
        if hold.sum() >= 50 and (~hold).sum() >= 100:
            for lam in lams:
                bb, be = fit_one(Xt[~hold], yt[~hold], ot[~hold], lam, a.maxiter)
                l = ll(sig(ot[hold] + bb + Xt[hold] @ be), yt[hold])
                if l < best_l - 1e-5: chosen, best_l = lam, l
        pv0 = sig(ov); l0 = ll(pv0, yv)
        if chosen is None:                                     # the other groups do not help this one on held-out training patents: leave the reagent rate alone
            stats.append((frag_ids[f], len(rows_v), l0, l0, 0)); continue
        bb, be = fit_one(Xt, yt, ot, chosen, a.maxiter)
        l1 = ll(sig(ov + bb + Xv @ be), yv)
        b[f] = bb; beta[f, keep] = be[:-1]; mult[f] = be[-1]
        stats.append((frag_ids[f], len(rows_v), l0, l1, chosen))
        if f % 20 == 0: print(f'  fragment {f + 1}/{nF} ({time.time() - t0:.0f} s)', flush=True)
    tot0 = sum(s[1] * s[2] for s in stats); tot1 = sum(s[1] * s[3] for s in stats); n = sum(s[1] for s in stats)
    print(f'validation log-loss of "is this group consumed" over {n:,} (reaction, fragment) pairs: reagent rate alone {tot0 / n:.4f}; with the other groups {tot1 / n:.4f}')
    for s in sorted(stats, key=lambda s: -(s[2] - s[3]) * s[1])[:25]:
        print(f'   {s[0]} n={s[1]:6d}  {s[2]:.4f} -> {s[3]:.4f}  (lambda {s[4]})')
    # strongest coefficients
    print('strongest coefficients (F consumed given G present):')
    flat = [(abs(beta[f, g]), f, g) for f in range(nF) for g in range(nF) if beta[f, g] != 0]
    for _, f, g in sorted(flat, reverse=True)[:30]:
        print(f'   {frag_ids[f]} given {frag_ids[g]}: {beta[f, g]:+.2f}')
    out = {'fragments': frag_ids, 'b': b, 'beta': beta.astype(np.float32), 'mult2': mult, 'stats': stats, 'val_logloss': {'rate_only': tot0 / n, 'with_context': tot1 / n, 'pairs': n}, 'built': time.strftime('%Y-%m-%d %H:%M'), 'tag': a.tag}
    pickle.dump(out, open(os.path.join(TMP, f'organic_chem_{a.out}.pkl'), 'wb'), protocol=4)
    print(f'wrote organic_chem_{a.out}.pkl ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
