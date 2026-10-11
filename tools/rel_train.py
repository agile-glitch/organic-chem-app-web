"""Train the template-relevance model: which recorded reaction template does a reaction (starting materials + reagent) use?

    python tools/rel_train.py --tag exp1 --subset 250000 --epochs 3          # a quick experiment, measured on held-out patents
    python tools/rel_train.py --tag final --all --epochs 6 --export           # learn from EVERY learned reaction, write the page's data

The model (plain numpy, no framework): bits (tools/rel_features.py) -> hidden -> hidden -> softmax over the templates that
at least --minc training reactions used, plus one "other" class (a reaction whose template is rarer: the model then says
"none of the common ones", which is what keeps its confidence honest).

Measured, never claimed: with the default split 1 patent in 20 is held out and NOTHING learned sees it (the templates' class set,
the weights, the calibration). Reported per run: template top-1 / top-5 / top-20 / top-100 over the held-out reactions, and over
those whose template is one of the classes (the ceiling is how many held-out reactions use a template the training data has).
--drop removes feature blocks (e.g. --drop frag_sub,frag_reag) to measure what a block is worth.
"""
import argparse, json, os, sys, tempfile, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C
import rel_features as F

TMP = tempfile.gettempdir()


LABELS = os.path.join(TMP, 'organic_chem_rel_labels.npz')


def load_labels(meta):
    """template number and split of every learned reaction (cached, so a training run does not need the 1.3 GB of row objects)"""
    if os.path.exists(LABELS):
        z = np.load(LABELS)
        if z['y'].shape[0] == meta['rows'] and str(z['first']) == meta['first'] and str(z['last']) == meta['last']:
            return z['y'], z['split'], int(z['n_templates'])
    D = C.read_rows()
    rows = D['rows']
    assert len(rows) == meta['rows'] and rows[0][0] == meta['first'] and rows[-1][0] == meta['last'], 'feature table does not match the rows'
    y = np.array([r[4] for r in rows], dtype=np.int32)
    split = np.array(C.split_codes(rows), dtype=np.int8)               # 0 train, 1 validation, 2 test (by patent)
    np.savez(LABELS, y=y, split=split, n_templates=len(D['smarts']), first=rows[0][0], last=rows[-1][0])
    return y, split, len(D['smarts'])


def load(tag_feat='v1'):
    """(feature table (memmap), template numbers, split codes, feature meta, number of templates)"""
    global NBITS
    path, meta_path = F.table_paths(tag_feat)
    meta = json.load(open(meta_path))
    NBITS = meta['layout']['bits']
    X = np.load(path, mmap_mode='r')
    assert X.shape[0] == meta['rows']
    y, split, n_templates = load_labels(meta)
    return X, y, split, meta, n_templates


def column_index(meta, drop):
    L = meta['layout']
    keep = np.ones(L['bits'], bool)
    for name in [d for d in drop.split(',') if d]:
        a, w = L[name]; keep[a:a + w] = False
    return np.flatnonzero(keep)


def class_set(y_train, minc, n_templates):
    cnt = np.bincount(y_train, minlength=n_templates)
    ids = np.flatnonzero(cnt >= minc)
    ids = ids[np.argsort(-cnt[ids], kind='stable')]
    cls_of = np.full(cnt.shape[0], -1, np.int32); cls_of[ids] = np.arange(len(ids))
    return ids, cls_of


class Net:
    def __init__(self, d, hidden, n_cls, prior, seed=0):
        rng = np.random.default_rng(seed)
        self.sizes = [d] + list(hidden) + [n_cls]
        self.W, self.b = [], []
        for i in range(len(self.sizes) - 1):
            fan = self.sizes[i]
            std = np.sqrt(2.0 / fan) if i < len(self.sizes) - 2 else 0.02
            self.W.append((rng.standard_normal((self.sizes[i], self.sizes[i + 1]), dtype=np.float32) * np.float32(std)))
            self.b.append(np.zeros(self.sizes[i + 1], np.float32))
        self.b[-1] = np.log(prior + 1e-9).astype(np.float32)
        self.m = [np.zeros_like(a) for a in self.W + self.b]; self.v = [np.zeros_like(a) for a in self.W + self.b]; self.t = 0

    def forward(self, x, drop=0.0, rng=None, keep_act=False):
        acts, masks, a = [x], [], x
        for i in range(len(self.W) - 1):
            z = a @ self.W[i] + self.b[i]
            np.maximum(z, 0, out=z)
            if drop > 0:
                mk = (rng.random(z.shape, dtype=np.float32) >= drop).astype(np.float32) / np.float32(1 - drop)
                z *= mk; masks.append(mk)
            acts.append(z); a = z
        logits = a @ self.W[-1] + self.b[-1]
        return (logits, acts, masks) if keep_act else logits

    def step(self, x, y, lr, drop, rng, wd=0.0, smooth=0.0):
        logits, acts, masks = self.forward(x, drop, rng, keep_act=True)
        B = x.shape[0]
        logits -= logits.max(axis=1, keepdims=True)
        np.exp(logits, out=logits)
        s = logits.sum(axis=1, keepdims=True)
        logits /= s                                           # now the probabilities
        loss = float(-np.log(logits[np.arange(B), y] + 1e-12).mean())
        if smooth > 0: logits -= np.float32(smooth / logits.shape[1])            # gradient of the smoothed target: p - eps/C - (1-eps) onehot
        logits[np.arange(B), y] -= np.float32(1 - smooth)
        g = logits / np.float32(B)
        grads_W, grads_b = [None] * len(self.W), [None] * len(self.b)
        for i in range(len(self.W) - 1, -1, -1):
            grads_W[i] = acts[i].T @ g; grads_b[i] = g.sum(axis=0)
            if i > 0:
                g = g @ self.W[i].T
                if masks: g *= masks[i - 1]
                g *= (acts[i] > 0)
        self.adam(grads_W + grads_b, lr, wd)
        return loss

    def adam(self, grads, lr, wd, b1=0.9, b2=0.999, eps=1e-8):
        self.t += 1
        c1, c2 = 1 - b1 ** self.t, 1 - b2 ** self.t
        params = self.W + self.b
        for k, (p, gr) in enumerate(zip(params, grads)):
            m, v = self.m[k], self.v[k]
            m *= b1; m += (1 - b1) * gr
            v *= b2; v += (1 - b2) * gr * gr
            if wd and k < len(self.W): p *= (1 - lr * wd)
            p -= (lr / c1) * m / (np.sqrt(v / c2) + eps)

    def probs(self, x):
        z = self.forward(x)
        z -= z.max(axis=1, keepdims=True); np.exp(z, out=z); z /= z.sum(axis=1, keepdims=True)
        return z


NBITS = None                                            # set by load(): the number of real bits of a row (the last byte is padded)


def batch(X, idx, cols, dtype=np.float32, mask=None):
    """the unpacked bits of rows idx (columns cols), after mask(bits) if given (tools/eval_top1.py clears the solvent bits with it)"""
    xb = np.unpackbits(X[idx], axis=1, count=NBITS)
    if mask is not None: xb = mask(xb)
    return xb[:, cols].astype(dtype) if cols is not None else xb.astype(dtype)


def evaluate(net, X, cols, idx, y_cls, ks=(1, 3, 5, 10, 20, 50, 100), bs=2048):
    """template-level top-k on reactions idx (labels y_cls: class number, or the 'other' class)"""
    n_cls = net.sizes[-1]; other = n_cls - 1
    hit = {k: 0 for k in ks}; seen = 0; in_set = 0; loss = 0.0
    for a in range(0, len(idx), bs):
        j = idx[a:a + bs]
        p = net.probs(batch(X, j, cols))
        yy = y_cls[j]
        p_real = p[:, :other]
        top = np.argpartition(-p_real, max(ks), axis=1)[:, :max(ks)]
        order = np.take_along_axis(top, np.argsort(-np.take_along_axis(p_real, top, 1), axis=1), 1)
        for k in ks:
            hit[k] += int((order[:, :k] == yy[:, None]).any(axis=1).sum())
        in_set += int((yy != other).sum()); seen += len(j)
        loss += float(-np.log(p[np.arange(len(j)), yy] + 1e-12).sum())
    return {'n': seen, 'template_in_class_set': round(in_set / seen, 4), 'loss': round(loss / seen, 4),
            **{f'top{k}': round(hit[k] / seen, 4) for k in ks}, **{f'top{k}_of_covered': round(hit[k] / max(in_set, 1), 4) for k in ks[:4]}}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--feat', default='v1'); ap.add_argument('--tag', default='exp')
    ap.add_argument('--subset', type=int, default=0, help='train on this many training reactions (0 = all of them)')
    ap.add_argument('--all', action='store_true', help='no held-out patents: learn from every reaction, validation and test patents too (the model the page ships)')
    ap.add_argument('--minc', type=int, default=5); ap.add_argument('--hidden', default='1024,512')
    ap.add_argument('--epochs', type=int, default=3); ap.add_argument('--batch', type=int, default=1024)
    ap.add_argument('--lr', type=float, default=1e-3); ap.add_argument('--drop', type=float, default=0.2); ap.add_argument('--wd', type=float, default=0.0)
    ap.add_argument('--smooth', type=float, default=0.0)
    ap.add_argument('--dropblocks', default='', help='feature blocks to leave out, e.g. frag_sub,frag_reag')
    ap.add_argument('--test', type=int, default=20000); ap.add_argument('--seed', type=int, default=0)
    a = ap.parse_args()
    t0 = time.time()
    X, y, split, meta, n_templates = load(a.feat)
    rng = np.random.default_rng(a.seed)
    tr_all = np.flatnonzero(split == C.TRAIN) if not a.all else np.arange(len(y))
    te_all = np.flatnonzero(split == C.TEST)
    tr = rng.permutation(tr_all)[:a.subset] if a.subset else tr_all
    te = np.sort(rng.permutation(te_all)[:a.test])
    ids, cls_of = class_set(y[tr], a.minc, n_templates)
    n_cls = len(ids) + 1
    y_cls = np.where(cls_of[y] >= 0, cls_of[y], n_cls - 1).astype(np.int32)
    cols = column_index(meta, a.dropblocks) if a.dropblocks else None
    d = len(cols) if cols is not None else meta['layout']['bits']
    cnt = np.bincount(y_cls[tr], minlength=n_cls).astype(np.float64)
    print(f'train {len(tr):,} reactions, test {len(te):,} (held-out patents: {not a.all}); classes {n_cls - 1:,} templates seen >= {a.minc}x + other; '
          f'inputs {d}; "other" share of training {cnt[-1] / cnt.sum():.1%}; test reactions with a class template '
          f'{(y_cls[te] != n_cls - 1).mean():.1%} ({time.time() - t0:.0f} s)', flush=True)
    hidden = [int(h) for h in a.hidden.split(',') if h]
    net = Net(d, hidden, n_cls, cnt / cnt.sum(), a.seed)
    steps_per_epoch = len(tr) // a.batch
    total = steps_per_epoch * a.epochs
    done = 0
    for ep in range(a.epochs):
        order = rng.permutation(tr)
        run = 0.0
        for s in range(steps_per_epoch):
            j = order[s * a.batch:(s + 1) * a.batch]
            lr = a.lr * (0.5 * (1 + np.cos(np.pi * done / total)) * 0.95 + 0.05)          # cosine decay to 5 %
            run += net.step(batch(X, j, cols), y_cls[j], lr, a.drop, rng, a.wd, a.smooth)
            done += 1
            if (s + 1) % 100 == 0:
                print(f'  epoch {ep + 1}/{a.epochs} step {s + 1}/{steps_per_epoch} loss {run / 100:.3f} ({time.time() - t0:.0f} s)', flush=True); run = 0.0
        if not a.all:
            r = evaluate(net, X, cols, te, y_cls)
            print(f'== epoch {ep + 1}: held-out patents: ' + json.dumps(r), flush=True)
        np.savez(os.path.join(TMP, f'organic_chem_relmodel_{a.tag}.npz'), **{f'W{i}': w for i, w in enumerate(net.W)}, **{f'b{i}': b for i, b in enumerate(net.b)},
                 classes=ids, cols=(cols if cols is not None else np.arange(d)), minc=a.minc, feat=a.feat)
    print(f'done ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
