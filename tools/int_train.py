"""Learn how much each interference matters, from the regioselectivity contests (tools/int_contests.py).

    python tools/int_train.py --train int_contests_train.pkl --val int_contests_val.pkl --tag v1
    python tools/int_train.py ... --grid                   # choose the two penalties on the validation contests

Model: a conditional logit over the outcomes of one contest (the softmax of  score = sum of weights of the flags the outcome has  among the outcomes;
the recorded one is the target). Each flag has a GLOBAL weight shared by every kind of reaction and a weight per reaction FAMILY (tools/interference_defs.py
FAMILIES: SNAr and coupling, electrophilic substitution, SN2/SN1, acyl substitution ...) that is pulled toward zero by a heavier penalty, so a family
deviates from the global value only where its contests say so. The ring-substituent flags also feed three Hammett terms (the sum of the substituent
constants σ at the ortho, meta and para positions) with a global and a per-family coefficient ρ: the reaction constant of physical organic chemistry,
learned from the recorded reactions instead of tabulated.

Everything is fitted on contests of the TRAINING patents; the penalties are chosen on the validation patents; the figures printed are on validation
contests (the shipped figures on the test patents come from tools/rel_fuse.py with the model applied to the whole candidate list).
"""
import argparse, collections, json, os, pickle, sys, tempfile, time

import numpy as np
from scipy import sparse
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import interference_defs as D

TMP = tempfile.gettempdir()
KEYS = [s['key'] for s in D.SITE]
KIDX = {k: i for i, k in enumerate(KEYS)}
FAMS = D.FAMILY_NAMES
FIDX = {f: i for i, f in enumerate(FAMS)}
NF, NFAM = len(KEYS), len(FAMS)
POS = ['o', 'm', 'p']
SIG_OF = {}                        # flag index -> (position index, sigma)
for s in D.SITE:
    if s['sigma'] is not None:
        rel = s['key'].split('-')[1]
        SIG_OF[KIDX[s['key']]] = ({'ortho': 0, 'meta': 1, 'para': 2}[rel], s['sigma'])
NSIG = 3


def sigma_terms(flag_idx):
    z = [0.0, 0.0, 0.0]
    for f in flag_idx:
        t = SIG_OF.get(f)
        if t is not None: z[t[0]] += t[1]
    return z


def load_contests(name, min_varying=True):
    d = pickle.load(open(os.path.join(TMP, name), 'rb'))
    out = []
    for c in d['contests']:
        outs = [(bool(o[0]), tuple(sorted(KIDX[k] for k in o[1] if k in KIDX))) for o in c['outs']]
        if min_varying and len({fl for _, fl in outs}) < 2: continue
        out.append({'ri': c['ri'], 'fam': FIDX[D.family_of(set(c['consumed']))], 'outs': outs})        # the family from the consumed fragments with the CURRENT family list
    return out, d['stat']


class Design:
    """sparse design for a list of contests: columns = [global flags, global ρ(o,m,p) | per family: flags, ρ(o,m,p)]"""
    def __init__(self, contests):
        self.contests = contests
        rows_i, cols_g, cols_f, y, first, fam_of_row = [], [], [], [], [], []
        z = []
        n = 0
        for c in contests:
            first.append(n)
            for h, fl in c['outs']:
                for f in fl:
                    rows_i.append(n); cols_g.append(f); cols_f.append(NF + NSIG + c['fam'] * (NF + NSIG) + f)
                zt = sigma_terms(fl)
                z.append(zt)
                y.append(1.0 if h else 0.0); fam_of_row.append(c['fam'])
                n += 1
        self.n = n; self.first = np.array(first, np.int64)
        self.gid = np.repeat(np.arange(len(contests)), [len(c['outs']) for c in contests])
        self.y = np.array(y); self.fam = np.array(fam_of_row)
        ncol = (NF + NSIG) * (1 + NFAM)
        data = np.ones(len(rows_i))
        Xg = sparse.csr_matrix((data, (rows_i, cols_g)), shape=(n, ncol))
        Xf = sparse.csr_matrix((data, (rows_i, cols_f)), shape=(n, ncol))
        Z = np.array(z) if z else np.zeros((0, 3))
        zi, zc_g, zc_f, zv = [], [], [], []
        for r in range(n):
            for k in range(NSIG):
                if Z[r, k] != 0:
                    zi.append(r); zc_g.append(NF + k); zc_f.append(NF + NSIG + self.fam[r] * (NF + NSIG) + NF + k); zv.append(Z[r, k])
        Zg = sparse.csr_matrix((zv, (zi, zc_g)), shape=(n, ncol)); Zf = sparse.csr_matrix((zv, (zi, zc_f)), shape=(n, ncol))
        self.A = (Xg + Xf + Zg + Zf).tocsr()
        self.Ag = (Xg + Zg).tocsr()                                  # the global part alone (for the "no families" baseline)

    def penalty(self, lam_g, lam_d, lam_rho_g, lam_rho_d):
        lam = np.zeros(self.A.shape[1])
        lam[:NF] = lam_g; lam[NF:NF + NSIG] = lam_rho_g
        for f in range(NFAM):
            a = NF + NSIG + f * (NF + NSIG)
            lam[a:a + NF] = lam_d; lam[a + NF:a + NF + NSIG] = lam_rho_d
        return lam

    def scores(self, w, A=None):
        return (self.A if A is None else A) @ w

    def top1(self, w, A=None):
        s = self.scores(w, A)
        mx = np.maximum.reduceat(s, self.first)[self.gid]
        best = (s >= mx - 1e-12)
        share = np.add.reduceat(best.astype(float), self.first)[self.gid]
        right = np.add.reduceat(best * self.y / share, self.first)                      # ties split evenly
        return right

    def loglik(self, w, A=None):
        s = self.scores(w, A)
        mx = np.maximum.reduceat(s, self.first)[self.gid]
        e = np.exp(s - mx)
        z = np.add.reduceat(e, self.first)
        zh = np.add.reduceat(e * self.y, self.first)
        return np.log(zh + 1e-300) - np.log(z)


def fit(des, lam, maxiter=300, w0=None):
    A, y, first, gid = des.A, des.y, des.first, des.gid
    n_c = len(first)

    def f(w):
        s = A @ w
        mx = np.maximum.reduceat(s, first)[gid]
        e = np.exp(s - mx)
        z = np.add.reduceat(e, first)[gid]
        p = e / z
        eh = e * y
        zh = np.add.reduceat(eh, first)[gid]
        ph = eh / zh
        loss = -(np.log(np.add.reduceat(eh, first)) - np.log(np.add.reduceat(e, first))).sum() / n_c + 0.5 * (lam * w * w).sum() / n_c
        grad = A.T @ (p - ph) / n_c + lam * w / n_c
        return loss, grad
    r = minimize(f, np.zeros(A.shape[1]) if w0 is None else w0, jac=True, method='L-BFGS-B', options={'maxiter': maxiter})
    return r.x


def effective(w):
    """per family: (weights of the flags, rho(o,m,p)) with the global part added"""
    g = w[:NF + NSIG]
    out = {}
    for f, name in enumerate(FAMS):
        a = NF + NSIG + f * (NF + NSIG)
        out[name] = g + w[a:a + NF + NSIG]
    return g, out


def evidence(contests, patent_of=None):
    """per flag and family: contests where the flag is present in some but not all outcomes, how often the recorded outcome had it, and how often a random pick would;
    with patent_of (reaction index -> patent), also the number of different patents those contests come from (a patent repeats the same step many times)"""
    stat = collections.defaultdict(lambda: [0, 0.0, 0.0])
    pats = collections.defaultdict(set)
    for c in contests:
        outs = c['outs']; m = len(outs)
        cnt = collections.Counter()
        for _, fl in outs:
            for f in fl: cnt[f] += 1
        win = [fl for h, fl in outs if h]
        for f, k in cnt.items():
            if k == m: continue
            wf = np.mean([f in fl for fl in win])
            for key in ((f, 'all'), (f, FAMS[c['fam']])):
                s = stat[key]; s[0] += 1; s[1] += wf; s[2] += k / m
                if patent_of is not None: pats[key].add(int(patent_of[c['ri']]))
    return {k: dict({'n': v[0], 'win': round(v[1] / v[0], 4), 'expected': round(v[2] / v[0], 4)}, **({'patents': len(pats[k])} if patent_of is not None else {})) for k, v in stat.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--train', required=True); ap.add_argument('--val', required=True); ap.add_argument('--tag', default='v1')
    ap.add_argument('--lam-g', type=float, default=3.0); ap.add_argument('--lam-d', type=float, default=30.0)
    ap.add_argument('--grid', action='store_true'); ap.add_argument('--maxiter', type=int, default=300)
    a = ap.parse_args()
    t0 = time.time()
    tr, st = load_contests(a.train); va, sv = load_contests(a.val)
    import int_context as IC
    patent_of = IC.load_meta()['patent']
    print(f'contests with flags that differ: train {len(tr):,} (of {st.get("contest", 0):,}), validation {len(va):,} (of {sv.get("contest", 0):,}) ({time.time() - t0:.0f} s)', flush=True)
    dt, dv = Design(tr), Design(va)
    print(f'design built ({time.time() - t0:.0f} s): train {dt.n:,} outcomes x {dt.A.shape[1]} columns', flush=True)
    rand = float(np.mean([np.mean([h for h, _ in c['outs']]) for c in va]))
    print(f'validation: a random pick among the outcomes is right {rand:.3f}')
    best = (a.lam_g, a.lam_d)
    if a.grid:
        res = []
        for lg in (0.3, 1.0, 3.0, 10.0):
            for ld in (3.0, 10.0, 30.0, 100.0):
                w = fit(dt, dt.penalty(lg, ld, 0.3, ld / 10), a.maxiter)
                ll = float(dv.loglik(w).mean()); t1 = float(dv.top1(w).mean())
                res.append((ll, t1, lg, ld)); print(f'  lam_g {lg:5.1f} lam_d {ld:6.1f}: validation log-lik {ll:.4f} top-1 {t1:.4f}', flush=True)
        best = max(res)[2:]
        print('chosen penalties', best)
    lam = dt.penalty(best[0], best[1], 0.3, best[1] / 10)
    w = fit(dt, lam, a.maxiter)
    # baselines on validation: the global part alone
    lam_g_only = dt.penalty(best[0], 1e9, 0.3, 1e9)
    wg = fit(dt, lam_g_only, a.maxiter)
    t_full, t_glob = float(dv.top1(w).mean()), float(dv.top1(wg).mean())
    print(f'validation contests: random {rand:.3f}; global weights only {t_glob:.3f}; with families {t_full:.3f}; mean log-lik full {float(dv.loglik(w).mean()):.4f} global {float(dv.loglik(wg).mean()):.4f} random {float(np.mean([np.log(np.mean([h for h, _ in c["outs"]])) for c in va])):.4f}')
    by_fam = collections.defaultdict(lambda: [0, 0.0, 0.0, 0.0])
    rt, rg = dv.top1(w), dv.top1(wg)
    for c, t1, tg in zip(va, rt, rg):
        r = by_fam[FAMS[c['fam']]]
        r[0] += 1; r[1] += t1; r[2] += tg; r[3] += np.mean([h for h, _ in c['outs']])
    print('by family (validation contests): n, random, global only, with families')
    for k, r in sorted(by_fam.items(), key=lambda kv: -kv[1][0]):
        print(f'   {k:20s} {r[0]:6d}  {r[3] / r[0]:.3f}  {r[2] / r[0]:.3f}  {r[1] / r[0]:.3f}')
    g, eff = effective(w)
    print('global Hammett coefficients rho (ortho, meta, para):', np.round(g[NF:], 3))
    for fam in FAMS:
        print(f'   rho {fam:20s}', np.round(eff[fam][NF:], 3))
    order = np.argsort(-np.abs(g[:NF]))[:30]
    print('largest global weights:')
    for i in order: print(f'   {KEYS[i]:34s} {g[i]:+.3f}')
    model = {'w': w, 'keys': KEYS, 'families': FAMS, 'penalties': best, 'val': {'random': rand, 'global': t_glob, 'families': t_full, 'n': len(va)}, 'built': time.strftime('%Y-%m-%d %H:%M'),
             'evidence': evidence(tr, patent_of), 'train_contests': len(tr)}
    pickle.dump(model, open(os.path.join(TMP, f'organic_chem_intmodel_{a.tag}.pkl'), 'wb'), protocol=4)
    print(f'wrote organic_chem_intmodel_{a.tag}.pkl ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
