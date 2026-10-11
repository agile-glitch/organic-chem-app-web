"""Rank the candidate products of a reaction and say how far to trust the first one.

    python tools/rel_fuse.py --val val_cands.pkl --test test_cands.pkl --prior new --sources mlp,g3        # fit on validation patents, score on test patents
    python tools/rel_fuse.py ... --out data/learned/m/fusion.js --report tools/eval_top1.json            # also write the page's weights and the measured figures

Input: the candidate dumps of tools/eval_top1.py (--dump): per reaction the products the relevance model's templates give (and the
same templates with every unchanged leaf atom stripped, source g3), with the probability mass behind each, how often the reaction's
reagent used the template in the training patents (source rc: a feature, never a candidate) and its Morgan difference fingerprint
(environments formed / destroyed).

Ranking = a conditional logit over the candidates of one reaction (softmax of  w . features  among its candidates; the recorded
product is the target; stereo variants of one structure are merged):
  dense    log mass behind the product (sum, max over templates), log number of templates, log share of the reaction's mass,
           top-of-list flag, size of the product relative to the starting materials, number of molecules, the mass x
           (reagent given / model unsure) interactions, what the product consumed / left of the fragments present judged by the
           reagent's recorded reactivity (tools/rel_reactivity.py), the relaxed-template mass, how often this reagent used the template
  sparse   one weight per formed / destroyed Morgan environment (2 x 2048): which new bonds and groups a reaction may make
The ranker is a bag of K models, each fitted on K-1 folds of the validation patents (folds by patent); the page gets their average.
Confidence = a logistic model of "the first product is right" on the first product's softmax share, margin, the model's own
probability of "none of the common templates", whether a reagent was given, whether that reagent was ever recorded with the template ...
followed by an isotonic map (the logistic alone was too sure of itself in the middle). Everything is fitted on the VALIDATION patents
only (the confidence on out-of-fold scores) and measured on the TEST patents; nothing learned has seen a test patent.
"""
import argparse, collections, hashlib, json, os, pickle, sys, tempfile, time

import numpy as np
from scipy import sparse
from scipy.optimize import minimize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C

TMP = tempfile.gettempdir()
DIFF_BITS = 2048
DENSE = ['log_sum', 'log_max', 'log_n', 'log_share', 'is_first', 'log_size', 'extra_mols', 'sum_x_reagent', 'sum_x_other', 'share_x_reagent',
         'prior_consumed', 'prior_left', 'prior_both', 'n_consumed', 'n_formed', 'log_sum_g3', 'only_gen', 'frac_used', 'rc_log', 'rc_seen']
DENSE_INT = ['int_raw', 'int_rel', 'int_lshare']                         # the learned interference score of the site (tools/int_model.py IntModel) and its gap to the best site of the same relaxed template
DENSE_CTX = ['ctx_consumed', 'ctx_left', 'ctx_both', 'ctx_ll']       # what the candidate consumed / left, judged with the competing groups present (tools/int_model.py CtxModel)
EXTRA = []                                                 # which of them this run uses (main() sets it from --int-model / --ctx-model); the page reads the names from the weights file
EPS = 1e-5
# What the Reactions page shows (js/reactions.js reads these from the figures stored with the model, so it does what was measured): an answer only when
# its stated chance reaches the floor of its kind, and never for starting materials this small (the learned data has almost none: 5 of 12,000
# held-out reactions have fewer than 8 heavy atoms in all, so no figure can be measured there and the network only extrapolates).
MIN_HEAVY = 8
FLOORS = {'reagent': 0.5, 'no_reagent': 0.7, 'plain_reagent': 0.85}        # plain_reagent = a reagent without a reactive role in the Reagents sheet (spectator ions, parts it does not list)


def reagent_kind(info):
    return 'no_reagent' if not info['has_reagent'] else ('plain_reagent' if info.get('role_sig') is None else 'reagent')


def page_shows(info, p):
    """would the page show this answer? (the guard on size, then the floor of the reagent's kind)"""
    return (info.get('sub_heavy') or 0) >= MIN_HEAVY and p >= FLOORS[reagent_kind(info)]
# what the page predicted before the relevance model (templates learned with the exact reagent, then with its parts, then with reagents of the same roles),
# measured by the earlier harness (tools/eval_levers.py of the branch that wrote it; not part of this change) on 2,000 reactions of these same held-out
# patents, with its win rates learned from 30,000 reactions. Not re-measured by tools/eval_top1.py, which is why the figures are stated, not computed.
BASELINE = {'top1': 0.452, 'with_reagent_top1': 0.5885, 'n': 2000, 'note': 'the page\'s earlier prediction (templates learned with the exact reagent, then its parts, then reagents with the same roles), '
            'measured on 2,000 reactions of the same held-out patents; it answers nothing when no reagent is given'}


def load_dump(name):
    d = pickle.load(open(os.path.join(TMP, name), 'rb'))
    return d


def write_json_part(path, obj):
    """data/learned/m/<name>.js in the form the page's loader reads: window.LEARNED_PARTS["m/<name>"] = {gz: base64(gzip(json))}
    (a path ending in .json is written as plain JSON instead)"""
    if path.endswith('.json'):
        json.dump(obj, open(path, 'w'), separators=(',', ':')); return
    import base64, gzip
    name = os.path.basename(path)[:-3]
    gz = base64.b64encode(gzip.compress(json.dumps(obj, separators=(',', ':')).encode('utf-8'), 9)).decode('ascii')
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'(window.LEARNED_PARTS=window.LEARNED_PARTS||{{}})["m/{name}"]={{"gz":"{gz}"}};\n')


def build_groups(dump, rows, prior=None, sources=None, intm=None, ctxm=None):
    """one group per reaction: the flat-merged candidates with their features. rows = the learned rows (for the patent id and the recorded product).
    prior = a rel_reactivity.Prior: adds what each candidate consumed (and left) judged by the reagent's recorded reactivity order"""
    groups = []
    if prior is not None and dump.get('fragments') != prior.source_ids:
        raise SystemExit('the candidates and the reactivity table were made with different fragment lists: rebuild one of them')
    for ri, cands in dump['res'].items():
        info = dump['info'][ri]
        merged = collections.OrderedDict()
        for key, c in cands.items():
            if sources is not None: c = dict(c, src={k: v for k, v in c['src'].items() if k in sources or k == 'rc'})      # 'rc' is a feature, never a candidate source
            if not [k for k in c['src'] if k != 'rc']: continue
            m = c['src'].get('mlp')
            w = (m[0] if m else 0.0) + 1e-3 * sum(v[0] for k, v in c['src'].items() if k != 'mlp')       # which stereo variant speaks for the group
            g = merged.get(c['flat'])
            if g is None:
                g = merged[c['flat']] = {'key': key, 'flat': c['flat'], 'best': -1.0, 'sum': 0.0, 'max': 0.0, 'n': 0, 'hit': False, 'hit_exact': False, 'formed': None, 'destroyed': None, 'heavy': None,
                                         'n_mols': None, 'counts': None, 'srcs': {}, 'used': 0, 'flags': None, 'grp': 0}
            for name, v in c['src'].items():
                a = g['srcs'].setdefault(name, [0.0, 0.0, 0]); a[0] += v[0]; a[1] = max(a[1], v[1]); a[2] += v[2]
            if m: g['sum'] += m[0]; g['max'] = max(g['max'], m[1]); g['n'] += m[2]
            g['hit'] = g['hit'] or c['hit_flat']; g['hit_exact'] = g['hit_exact'] or c['hit']
            if w >= g['best'] or g['formed'] is None:
                if w >= g['best']: g['best'] = w; g['key'] = key; g['used'] = c.get('used', 0); g['flags'] = c.get('flags'); g['grp'] = c.get('grp', 0)
                if 'formed' in c: g['formed'], g['destroyed'], g['heavy'], g['n_mols'], g['counts'] = c['formed'], c['destroyed'], c.get('heavy'), c.get('n_mols'), c.get('counts')
        cs = list(merged.values())
        if prior is not None and info.get('s_counts'):
            for c in cs:
                c['prior'] = prior.features(info['s_counts'], c['counts'], info['exact_key'], info['role_sig']) if c['counts'] is not None else (0.0, 0.0, 0, 0)
        if ctxm is not None and prior is not None and info.get('s_counts'):
            p_ctx = ctxm.rates(prior.rates(info['exact_key'], info['role_sig']), info['s_counts'])
            for c in cs:
                c['ctx'] = ctxm.features(p_ctx, info['s_counts'], c['counts']) if c['counts'] is not None else (0.0, 0.0, 0, 0)
        if intm is not None:
            import interference_defs as ID
            site_keys = dump.get('site_keys')
            if site_keys is None: raise SystemExit('the candidate dump has no interference flags: make it with tools/eval_top1.py --marks')
            remap = [intm.kidx.get(k, -1) for k in site_keys]
            ids = dump['fragments']
            s_arr = np.asarray(info['s_counts']) if info.get('s_counts') else None
            for c in cs:
                fam = 'other'
                if c['counts'] is not None and s_arr is not None:
                    fam = ID.family_of({ids[i] for i in np.flatnonzero(s_arr > np.asarray(c['counts']))})
                idx = [remap[i] for i in (c['flags'] or []) if remap[i] >= 0]
                c['fam'] = fam
                c['has_int'] = c['flags'] is not None and c['counts'] is not None       # (the page works the family out from the fragments the product consumed, which it only counts for the leading 200 candidates)
                c['int_raw'] = intm.score_idx(idx, fam) if c['has_int'] else 0.0
            best = {}
            for c in cs:
                if c['has_int']: best[c['grp']] = max(best.get(c['grp'], -1e9), c['int_raw'])
            lse = {}
            for c in cs:
                if c['has_int']: lse.setdefault(c['grp'], []).append(c['int_raw'])
            lse = {g: float(np.log(np.exp(np.array(v) - max(v)).sum()) + max(v)) for g, v in lse.items()}
            for c in cs:
                c['int_rel'] = (c['int_raw'] - best[c['grp']]) if c['has_int'] else 0.0
                c['int_lshare'] = (c['int_raw'] - lse[c['grp']]) if c['has_int'] else 0.0           # log of the site's share among the sites of its relaxed template (what the contest model fits)
        groups.append({'ri': ri, 'cands': cs, 'info': info, 'patent': rows[ri][8]})
    return groups


def extra_value(c, name):
    if name == 'int_raw': return c.get('int_raw', 0.0)
    if name == 'int_rel': return c.get('int_rel', 0.0)
    if name == 'int_lshare': return c.get('int_lshare', 0.0)
    ctx = c.get('ctx', (0.0, 0.0, 0, 0, 0.0))
    if name == 'ctx_consumed': return ctx[0]
    if name == 'ctx_left': return ctx[1]
    if name == 'ctx_both': return ctx[0] + ctx[1]
    if name == 'ctx_ll': return ctx[4] if len(ctx) > 4 else 0.0
    raise KeyError(name)


def design(groups):
    """(dense matrix, sparse formed/destroyed matrix, group index per candidate, group start offsets, label per candidate)"""
    dense, rows_i, cols_i, gid, starts, label, exact = [], [], [], [], [], [], []
    n = 0
    for gi, g in enumerate(groups):
        cs = g['cands']
        starts.append(n)
        if not cs: continue
        tot = sum(c['sum'] for c in cs) or 1.0
        order = sorted(range(len(cs)), key=lambda i: -cs[i]['sum'])
        first = order[0]
        inf = g['info']
        has_r = 1.0 if inf['has_reagent'] else 0.0
        sub_heavy = max(1, inf.get('sub_heavy') or 1)
        for i, c in enumerate(cs):
            share = c['sum'] / tot
            size = np.log(max(1, c['heavy'] or sub_heavy) / sub_heavy)
            pa, pb, nc, nf = c.get('prior', (0.0, 0.0, 0, 0))
            dense.append([np.log(c['sum'] + EPS), np.log(c['max'] + EPS), np.log1p(c['n']), np.log(share + EPS), 1.0 if i == first else 0.0, float(np.clip(size, -3, 3)),
                          float((c['n_mols'] or 1) - 1), np.log(c['sum'] + EPS) * has_r, np.log(c['sum'] + EPS) * inf['p_other'], np.log(share + EPS) * has_r,
                          pa, pb, pa + pb, float(nc), float(nf),
                          np.log(c['srcs'].get('g3', [0.0])[0] + EPS), 1.0 if c['sum'] == 0 else 0.0,
                          (c.get('used', 0) / max(1, inf['n_subs'])) if c.get('used') is not None else 1.0,
                          np.log1p(c['srcs'].get('rc', [0.0, 0.0])[1]), 1.0 if (has_r and c['srcs'].get('rc', [0.0, 0.0])[1] > 0) else 0.0]
                         + [extra_value(c, name) for name in EXTRA])
            for b in (c['formed'] or []): rows_i.append(n); cols_i.append(b)
            for b in (c['destroyed'] or []): rows_i.append(n); cols_i.append(DIFF_BITS + b)
            gid.append(gi); label.append(1.0 if c['hit'] else 0.0); exact.append(1.0 if c['hit_exact'] else 0.0)
            n += 1
    sp = sparse.csr_matrix((np.ones(len(rows_i), np.float32), (rows_i, cols_i)), shape=(n, 2 * DIFF_BITS))
    return (np.array(dense, np.float64), sp, np.array(gid, np.int64), np.array(starts, np.int64), np.array(label, np.float64), np.array(exact, np.float64))


class Fusion:
    def __init__(self, lam_dense=1e-3, lam_sparse=3.0): self.lam_d, self.lam_s = lam_dense, lam_sparse; self.w = None; self.mu = None; self.sd = None

    def _stack(self, dense, sp, fit=False):
        if fit: self.mu = dense.mean(0); self.sd = np.where(dense.std(0) < 1e-6, 1.0, dense.std(0))
        d = (dense - self.mu) / self.sd
        return d, sp

    def _scores(self, d, sp, w):
        return d @ w[:d.shape[1]] + sp @ w[d.shape[1]:]

    def fit(self, groups, maxiter=200, stats=None):
        """stats = (mean, sd) of the dense features, shared by the models of a bag so that the average of their weights is their average score"""
        dense, sp, gid, starts, label, _ = design(groups)
        if stats is not None: self.mu, self.sd = stats
        d, sp = self._stack(dense, sp, fit=stats is None)
        has_hit = np.zeros(len(groups), bool); np.logical_or.at(has_hit, gid, label > 0)
        use = has_hit[gid]                                                   # groups whose recorded product was generated at all
        d, sp, gid, label = d[use], sp[use], gid[use], label[use]
        gids, inv = np.unique(gid, return_inverse=True)
        first = np.r_[0, np.flatnonzero(np.diff(inv)) + 1]
        nd = d.shape[1]; w0 = np.zeros(nd + sp.shape[1])
        lam = np.r_[np.full(nd, self.lam_d), np.full(sp.shape[1], self.lam_s)]
        n_groups = len(gids)

        def f(w):
            s = self._scores(d, sp, w)
            m = np.maximum.reduceat(s, first)[inv]
            e = np.exp(s - m)
            z = np.add.reduceat(e, first)[inv]
            p = e / z
            eh = e * label
            zh = np.add.reduceat(eh, first)[inv]
            ph = eh / zh
            loss = -(np.log(np.add.reduceat(eh, first)) - np.log(np.add.reduceat(e, first))).sum() / n_groups + 0.5 * (lam * w * w).sum() / n_groups
            r = p - ph
            grad = np.r_[d.T @ r, sp.T @ r] / n_groups + lam * w / n_groups
            return loss, grad

        res = minimize(f, w0, jac=True, method='L-BFGS-B', options={'maxiter': maxiter})
        self.w = res.x
        return res

    def score_groups(self, groups):
        """per group: list of (candidate, score) best first; the softmax share of each"""
        dense, sp, gid, starts, label, exact = design(groups)
        d, sp = self._stack(dense, sp)
        s = self._scores(d, sp, self.w)
        out, k = [], 0
        for g in groups:
            n = len(g['cands'])
            if not n: out.append(None); continue
            sc = s[k:k + n]; k += n
            e = np.exp(sc - sc.max()); p = e / e.sum()
            out.append((sc, p))
        return out


def top1_stats(groups, scored):
    n = len(groups); hit_flat = hit_exact = answered = generated = 0
    for g, sc in zip(groups, scored):
        if sc is None: continue
        answered += 1
        cs = g['cands']
        generated += any(c['hit'] for c in cs)
        i = int(np.argmax(sc[0]))
        hit_flat += cs[i]['hit']; hit_exact += cs[i]['hit_exact']
    return {'n': n, 'answered': round(answered / n, 4), 'generated': round(generated / n, 4), 'top1_flat': round(hit_flat / n, 4), 'top1_exact': round(hit_exact / n, 4)}


_FORMULA = {}


def formula_of(flat):
    """molecular formula of a candidate (its products as one string)"""
    f = _FORMULA.get(flat)
    if f is None:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors
        m = Chem.MolFromSmiles(flat)
        f = _FORMULA[flat] = rdMolDescriptors.CalcMolFormula(m) if m is not None else '?'
    return f


def error_buckets(groups, scored):
    """what the wrong first answers are: the recorded product is not among the candidates at all (no ranking can help); an ISOMER of the recorded product was put first (the
    same atoms bonded another way: the wrong site, the regioselectivity); a product of another formula was put first (another reaction, the chemoselectivity)"""
    n = len(groups); c = {'right': 0, 'not generated': 0, 'wrong isomer (site)': 0, 'wrong reaction (formula differs)': 0, 'no answer': 0}
    for g, sc in zip(groups, scored):
        if sc is None: c['no answer'] += 1; continue
        cs = g['cands']
        i = int(np.argmax(sc[0]))
        if cs[i]['hit']: c['right'] += 1; continue
        hits = [x for x in cs if x['hit']]
        if not hits: c['not generated'] += 1; continue
        if formula_of(cs[i]['flat']) in {formula_of(x['flat']) for x in hits}: c['wrong isomer (site)'] += 1
        else: c['wrong reaction (formula differs)'] += 1
    return {k: v for k, v in c.items()}, {k: round(v / n, 4) for k, v in c.items()}


def bootstrap_ci(groups, scored, n_boot=1000, seed=7):
    """95 % intervals of the two top-1 figures, resampling whole PATENTS (reactions of one patent are alike, so resampling reactions would be too optimistic)"""
    by, order = collections.defaultdict(lambda: [0, 0, 0]), []
    for g, sc in zip(groups, scored):
        r = by[g['patent']]
        r[2] += 1
        if sc is not None:
            i = int(np.argmax(sc[0])); r[0] += g['cands'][i]['hit']; r[1] += g['cands'][i]['hit_exact']
    arr = np.array(list(by.values()), np.float64)
    rng = np.random.default_rng(seed)
    flat, exact = [], []
    for _ in range(n_boot):
        pick = arr[rng.integers(0, len(arr), len(arr))]
        flat.append(pick[:, 0].sum() / pick[:, 2].sum()); exact.append(pick[:, 1].sum() / pick[:, 2].sum())
    q = lambda v: [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]
    return {'top1_flat': q(flat), 'top1_exact': q(exact), 'patents': len(arr), 'resamples': n_boot}


def confidence_features(g, sc):
    """features of the first answer for the 'is it right' model"""
    scores, p = sc
    order = np.argsort(-scores)
    i = order[0]
    p1 = p[i]; p2 = p[order[1]] if len(order) > 1 else 0.0
    inf = g['info']
    c = g['cands'][i]
    unseen = 1.0 if (inf['has_reagent'] and c['srcs'].get('rc', [0.0, 0.0])[1] <= 0) else 0.0
    return np.array([np.log(p1 + 1e-4), np.log(p1 - p2 + 1e-4), np.log(c['sum'] + 1e-5), inf['p_other'], 1.0 if inf['has_reagent'] else 0.0, np.log1p(len(order)), inf['top_p'], unseen,
                     np.log(max(1, inf.get('sub_heavy') or 1))])


CONF_NAMES = ['log_share', 'log_margin', 'log_mass', 'p_other', 'has_reagent', 'log_n_candidates', 'top_template_p', 'reagent_never_with_template', 'log_sub_heavy']


def fit_logistic(F, y, lam=1e-2):
    mu, sd = F.mean(0), F.std(0) + 1e-9
    Z = np.c_[(F - mu) / sd, np.ones(len(F))]

    def f(w):
        z = Z @ w; p = 1 / (1 + np.exp(-z))
        loss = -(y * np.log(p + 1e-12) + (1 - y) * np.log(1 - p + 1e-12)).mean() + 0.5 * lam * (w[:-1] ** 2).sum() / len(y)
        g = Z.T @ (p - y) / len(y) + np.r_[lam * w[:-1], 0] / len(y)
        return loss, g
    r = minimize(f, np.zeros(Z.shape[1]), jac=True, method='L-BFGS-B')
    return {'w': r.x.tolist(), 'mu': mu.tolist(), 'sd': sd.tolist()}


def isotonic_knots(p, y, bins=40):
    """a monotone map from a score p to the measured chance: equal-count bins, then pool-adjacent-violators on the bin means"""
    order = np.argsort(p)
    ps, ys = p[order], y[order]
    edges = np.linspace(0, len(p), bins + 1).astype(int)
    xs = [float(ps[a:b].mean()) for a, b in zip(edges[:-1], edges[1:]) if b > a]
    vs = [float(ys[a:b].mean()) for a, b in zip(edges[:-1], edges[1:]) if b > a]
    ws = [float(b - a) for a, b in zip(edges[:-1], edges[1:]) if b > a]
    blocks = []                                              # [value, weight, x sum (weighted), count]
    for x, v, w in zip(xs, vs, ws):
        blocks.append([v, w, x * w])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            b2, b1 = blocks.pop(), blocks.pop()
            blocks.append([(b1[0] * b1[1] + b2[0] * b2[1]) / (b1[1] + b2[1]), b1[1] + b2[1], b1[2] + b2[2]])
    return [round(b[2] / b[1], 5) for b in blocks], [round(b[0], 5) for b in blocks]


def apply_calibrated(model, F):
    p = apply_logistic(model, F)
    iso = model.get('iso')
    return np.interp(p, iso['x'], iso['y']) if iso else p


def apply_logistic(model, F):
    Z = np.c_[(F - np.array(model['mu'])) / np.array(model['sd']), np.ones(len(F))]
    return 1 / (1 + np.exp(-(Z @ np.array(model['w']))))


def fold_of(patent, k=4):
    return int(hashlib.md5(('fold' + str(patent)).encode()).hexdigest(), 16) % k


def reliability(p, y, edges=(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.01)):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (p >= lo) & (p < hi)
        if m.any(): out.append({'stated': f'{lo:.2f}-{min(hi, 1):.2f}', 'n': int(m.sum()), 'mean_stated': round(float(p[m].mean()), 3), 'measured': round(float(y[m].mean()), 3)})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--val', required=True); ap.add_argument('--test', required=True); ap.add_argument('--out', default='')
    ap.add_argument('--lam-sparse', type=float, default=3.0); ap.add_argument('--maxiter', type=int, default=150); ap.add_argument('--folds', type=int, default=4)
    ap.add_argument('--prior', default='', help='tag of the reactivity table (tools/rel_reactivity.py --table), e.g. cur')
    ap.add_argument('--parity', default='', help='write what this ranker answers for the first 40 test reactions, to compare with the browser (js/relevance.js)')
    ap.add_argument('--report', default='', help='write the measured figures (test patents) here, for tools/rel_export.py --report and the Validation tab')
    ap.add_argument('--prior-ids', default='', help='use only these fragments in the prior: first103 (the original sheet) or a comma list of IDs')
    ap.add_argument('--sources', default='', help='candidate sources to use, e.g. mlp or mlp,g3 (default: every source in the dump)')
    ap.add_argument('--int-model', default='', help='tag of the learned interference weights (tools/int_train.py): adds the features int_raw, int_rel (the dumps must come from eval_top1.py --marks)')
    ap.add_argument('--analysis', default='', help='write the ranked candidates of every test reaction (top 6 with their features) to this pickle in the temp folder, to look at the errors')
    ap.add_argument('--results-json', default='', help='write the measured figures of this run (validation out of fold, test, error buckets) here, for an ablation table')
    ap.add_argument('--int-feats', default='raw,rel', help='which interference features to use, any of: raw (the learned score of the site), rel (its gap to the best site of the same relaxed template), lshare (log of its share among those sites)')
    ap.add_argument('--ctx-feats', default='consumed,left,both', help='which of ctx_consumed, ctx_left, ctx_both, ctx_ll to use')
    ap.add_argument('--ctx-model', default='', help='tag of the context-aware reactivity prior (tools/int_context.py): adds ctx_consumed, ctx_left, ctx_both')
    a = ap.parse_args()
    t0 = time.time()
    val, test = load_dump(a.val), load_dump(a.test)
    rows = C.rows_subset(set(val['res']) | set(test['res']))
    import rel_reactivity as RV
    pids = None
    if a.prior_ids == 'first103': pids = ['F%03d' % i for i in range(1, 104)]
    elif a.prior_ids: pids = a.prior_ids.split(',')
    prior = RV.Prior(a.prior, ids=pids) if a.prior else None
    srcs = set(a.sources.split(',')) if a.sources else None
    import int_model as IM
    intm = IM.IntModel(a.int_model) if a.int_model else None
    ctxm = IM.CtxModel(a.ctx_model) if a.ctx_model else None
    if ctxm is not None and prior is None: raise SystemExit('--ctx-model needs --prior (the reagent rates are its offset)')
    global EXTRA
    EXTRA = (['int_' + x for x in a.int_feats.split(',') if x] if intm is not None else []) + (['ctx_' + x for x in a.ctx_feats.split(',') if x] if ctxm is not None else [])
    gv, gt = build_groups(val, rows, prior, srcs, intm, ctxm), build_groups(test, rows, prior, srcs, intm, ctxm)
    print(f'validation {len(gv):,} reactions, test {len(gt):,} ({time.time() - t0:.0f} s)', flush=True)
    # the model's own ranking (sum of template mass), for comparison
    base = [(np.array([c['sum'] for c in g['cands']]), None) if g['cands'] else None for g in gt]
    mass_sum = top1_stats(gt, base)
    print('mass-sum ranking on test:', json.dumps(mass_sum))
    # a bag of K models, each fitted on K-1 folds of the validation patents (folds by patent); the confidence is fitted on their
    # out-of-fold scores, and the model the page ships is the average of the K weight vectors (= the average of their scores, since the
    # standardisation is shared), so the scores the confidence was fitted on and the scores it is applied to are of the same kind
    K = a.folds
    dense_all = np.vstack([design([g])[0] for g in gv if g['cands']])
    stats = (dense_all.mean(0), np.where(dense_all.std(0) < 1e-6, 1.0, dense_all.std(0)))
    folds = [[g for g in gv if fold_of(g['patent'], K) == k] for k in range(K)]
    oof, bag = {}, []
    for k in range(K):
        fu = Fusion(lam_sparse=a.lam_sparse); fu.fit([g for j in range(K) if j != k for g in folds[j]], a.maxiter, stats)
        bag.append(fu.w)
        sc_k = fu.score_groups(folds[k])
        for g, s_ in zip(folds[k], sc_k): oof[g['ri']] = s_
        print(f'  fold {k}: top-1 flat {top1_stats(folds[k], sc_k)["top1_flat"]:.4f} (fitted on the other {K - 1})', flush=True)
    oof_all = top1_stats(gv, [oof[g['ri']] for g in gv])
    print('validation, out of fold:', json.dumps(oof_all), flush=True)         # the figure to choose features by (the test patents are looked at once, at the end)
    F, y = [], []
    for g in gv:
        s_ = oof[g['ri']]
        if s_ is None: continue
        F.append(confidence_features(g, s_)); y.append(1.0 if g['cands'][int(np.argmax(s_[0]))]['hit'] else 0.0)
    F, y = np.array(F), np.array(y)
    cal = fit_logistic(F, y)
    xs_, ys_ = isotonic_knots(apply_logistic(cal, F), y)
    cal['iso'] = {'x': xs_, 'y': ys_}
    final = Fusion(lam_sparse=a.lam_sparse); final.mu, final.sd = stats; final.w = np.mean(bag, axis=0)
    print(f'bagged ranker = average of {K} fold models ({time.time() - t0:.0f} s)', flush=True)
    sc = final.score_groups(gt)
    out = {'validation_reactions': len(gv), 'test_reactions': len(gt), 'test': top1_stats(gt, sc), 'mass_sum_ranking': mass_sum, 'ci95': bootstrap_ci(gt, sc)}
    out['errors_test'] = error_buckets(gt, sc)[1]; print('test, what the answers are:', json.dumps(error_buckets(gt, sc)[0]))
    sc_v = [oof[g['ri']] for g in gv]; print('validation (out of fold), what the answers are:', json.dumps(error_buckets(gv, sc_v)[0]))
    Ft, yt, ptest = [], [], []
    for g, s in zip(gt, sc):
        if s is None: continue
        Ft.append(confidence_features(g, s)); yt.append(1.0 if g['cands'][int(np.argmax(s[0]))]['hit'] else 0.0)
    pt = apply_calibrated(cal, np.array(Ft)); yt = np.array(yt)
    out['reliability_test'] = reliability(pt, yt)
    N = len(gt)
    out['by_confidence'] = {f'>= {th}': {'answers': int((pt >= th).sum()), 'share_of_all': round(float((pt >= th).sum()) / N, 4), 'right': round(float(yt[pt >= th].mean()), 4) if (pt >= th).any() else None}
                            for th in (0.5, 0.7, 0.8, 0.9, 0.95)}
    out['ece_test'] = round(float(np.abs(pt - yt).mean()), 4)
    infos = [g['info'] for g, s_ in zip(gt, sc) if s_ is not None]
    shown = np.array([page_shows(inf, p) for inf, p in zip(infos, pt)])
    out['page_policy'] = {'floor': FLOORS['reagent'], 'floor_no_reagent': FLOORS['no_reagent'], 'floor_plain_reagent': FLOORS['plain_reagent'], 'min_heavy': MIN_HEAVY,
                          'shown_share_of_all': round(float(shown.sum()) / N, 4),
                          'shown_right': round(float(yt[shown].mean()), 4) if shown.any() else None, 'not_shown': int((~shown).sum())}
    heavy = np.array([inf.get('sub_heavy') or 0 for inf in infos]); kinds = np.array([reagent_kind(inf) for inf in infos])
    top = [g['cands'][int(np.argmax(s_[0]))] for g, s_ in zip(gt, sc) if s_ is not None]
    ex = np.array([1.0 if c['hit_exact'] else 0.0 for c in top])

    def slice_of(mask):
        if not mask.any(): return None
        sh = shown & mask
        return {'n': int(mask.sum()), 'share': round(float(mask.sum()) / N, 4), 'top1_exact': round(float(ex[mask].mean()), 4), 'stated': round(float(pt[mask].mean()), 4), 'right': round(float(yt[mask].mean()), 4),
                'shown_share': round(float(sh.sum()) / float(mask.sum()), 4), 'shown_right': round(float(yt[sh].mean()), 4) if sh.any() else None}
    out['by_size'] = {label: slice_of((heavy >= lo) & (heavy < hi)) for label, lo, hi in (('< 8', 0, 8), ('8-11', 8, 12), ('12-19', 12, 20), ('20-29', 20, 30), ('30 or more', 30, 10 ** 6))}
    out['by_reagent'] = {label: slice_of(kinds == k) for label, k in (('a reagent with a reactive role', 'reagent'), ('a reagent without one (spectator ions, parts the Reagents sheet lacks)', 'plain_reagent'), ('no reagent', 'no_reagent'))}
    grad = {g['ri']: C.gradable(rows[g['ri']]) for g in gt}
    for name, keep in (('with a reagent', lambda g: g['info']['has_reagent']), ('no reagent', lambda g: not g['info']['has_reagent']),
                       ('gradable records (Validation tab data hygiene)', lambda g: grad[g['ri']])):
        idx = [i for i, g in enumerate(gt) if keep(g)]
        out['test / ' + name] = top1_stats([gt[i] for i in idx], [sc[i] for i in idx])
    print(json.dumps(out, indent=1))
    if a.analysis:
        ana = []
        for g, s_ in zip(gt, sc):
            if s_ is None: continue
            order = np.argsort(-s_[0])[:6]
            r = rows[g['ri']]
            ana.append({'ri': g['ri'], 'patent': g['patent'], 'subs': r[1], 'reagent': r[2], 'product': r[3], 'info': {k: g['info'].get(k) for k in ('has_reagent', 'p_other', 'top_p', 'sub_heavy', 'role_sig')},
                        'cands': [{'flat': g['cands'][i]['flat'], 'score': float(s_[0][i]), 'share': float(s_[1][i]), 'hit': g['cands'][i]['hit'], 'hit_exact': g['cands'][i].get('hit_exact'), 'sum': g['cands'][i]['sum'], 'n': g['cands'][i]['n'],
                                   'prior': g['cands'][i].get('prior'), 'ctx': g['cands'][i].get('ctx'), 'flags': g['cands'][i].get('flags'), 'grp': g['cands'][i].get('grp'), 'fam': g['cands'][i].get('fam'),
                                   'int_raw': g['cands'][i].get('int_raw'), 'int_rel': g['cands'][i].get('int_rel'), 'used': g['cands'][i].get('used')} for i in order]})
        pickle.dump(ana, open(os.path.join(TMP, a.analysis), 'wb'), protocol=4)
    if a.results_json:
        json.dump({'validation_oof': oof_all, 'validation_errors': error_buckets(gv, sc_v)[0], 'test': out['test'], 'test_errors': error_buckets(gt, sc)[0], 'ci95': out['ci95'], 'features': DENSE + EXTRA,
                   'ece_test': out['ece_test'], 'shown': out['page_policy']}, open(a.results_json, 'w'), indent=1)
    if a.parity:
        par = []
        for g, s_ in zip(gt, sc):
            if s_ is None or len(par) >= 40: continue
            i = int(np.argmax(s_[0])); r = rows[g['ri']]
            order = np.argsort(-s_[0])
            par.append({'id': r[0], 'subs': r[1], 'reagent': r[2], 'flat': g['cands'][i]['flat'],
                        'key': g['cands'][i]['key'], 'share': float(s_[1][i]), 'score': float(s_[0][i]), 'sum': g['cands'][i]['sum'], 'n': len(g['cands']),
                        'conf': float(apply_calibrated(cal, confidence_features(g, s_)[None, :])[0]), 'second': g['cands'][int(order[1])]['key'] if len(order) > 1 else None})
        json.dump(par, open(a.parity, 'w'))
        print('wrote', a.parity)
    if a.report:
        gk = [k for k in out if k.startswith('test / gradable')]
        rep = {'n': out['test_reactions'], 'patents': len({g['patent'] for g in gt}), 'when': time.strftime('%Y-%m-%d'), 'sources': a.sources or 'all in the dump',
               'top1_exact': out['test']['top1_exact'], 'top1_flat': out['test']['top1_flat'], 'generated': out['test']['generated'], 'answered': out['test']['answered'],
               'with_reagent': out['test / with a reagent'], 'no_reagent': out['test / no reagent'], 'gradable': out[gk[0]] if gk else None,
               'mass_sum_ranking': out['mass_sum_ranking'], 'ci95': out['ci95'], 'by_confidence': out['by_confidence'], 'reliability': out['reliability_test'], 'floor': out['page_policy']['floor'], 'floor_no_reagent': out['page_policy']['floor_no_reagent'], 'floor_plain_reagent': out['page_policy']['floor_plain_reagent'], 'min_heavy': MIN_HEAVY,
               'shown': out['page_policy'], 'by_size': out['by_size'], 'by_reagent': out['by_reagent'],
               'baseline': BASELINE}
        json.dump(rep, open(a.report, 'w', encoding='utf-8'), indent=1)
        print('wrote', a.report)
    if a.out:
        sd = final.sd.tolist(); mu = final.mu.tolist()
        nd = len(DENSE) + len(EXTRA)
        weights = {'dense': DENSE + EXTRA, 'mu': mu, 'sd': sd, 'w_dense': final.w[:nd].tolist(), 'w_formed': [round(float(x), 5) for x in final.w[nd:nd + DIFF_BITS]],
                   'w_destroyed': [round(float(x), 5) for x in final.w[nd + DIFF_BITS:]], 'confidence': cal, 'confidence_features': CONF_NAMES,
                   'test_measured': out['test'], 'built': time.strftime('%Y-%m-%d %H:%M')}
        write_json_part(a.out, weights)
        print('wrote', a.out)


if __name__ == '__main__':
    main()
