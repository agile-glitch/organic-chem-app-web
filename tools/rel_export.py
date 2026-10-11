"""Write the relevance model for the page: data/learned/m/*.js  (the page loads them on first use, as the other data/learned/ files).

    python tools/rel_export.py --model full1 --out data/learned/m
    python tools/rel_export.py --model expB1 --check            # only measure what int8 rounding does to the predictions

Files (each is a script that fills window.LEARNED_PARTS["m/..."]; binary parts are base64, gzipped when that helps):
  meta        feature layout, fragment ids, role names, sizes, the held-out figures measured for this model
  w1 w2       first two layers: int8 weights (rows = inputs) + one float scale per output unit + float bias
  w3_0 ...    the output layer, cut into column chunks (a chunk = some template classes): int8 + scale + bias
  classes     the learned-template number of each class (the page's data/learned/t/ numbering)
  gen3        each class template with every unchanged leaf atom stripped (the 'relaxed' candidates)
  (fusion, prior: written by tools/rel_fuse.py --out and tools/rel_reactivity.py --export)
Weights are rounded to int8 with one scale per output unit; --check measures that the first choice of the model hardly ever changes.
"""
import argparse, base64, gzip, json, os, pickle, sys, tempfile, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import rel_common as C
import rel_features as F
import rel_train as T

TMP = tempfile.gettempdir()
CHUNK_CLASSES = 10000


def quantise(W):
    """int8 weights with one scale per output unit (column)"""
    s = np.abs(W).max(axis=0) / 127.0
    s[s == 0] = 1.0
    q = np.clip(np.rint(W / s), -127, 127).astype(np.int8)
    return q, s.astype(np.float32)


def dequantise(q, s):
    return q.astype(np.float32) * s


def b64(raw, z=True):
    data = gzip.compress(raw, 9) if z else raw
    return base64.b64encode(data).decode('ascii'), z


def write_part(out_dir, name, obj):
    """obj: {'meta': json-able, 'arrays': {key: ndarray}} -> data/learned/m/<name>.js"""
    parts, offset = {}, 0
    blob = b''
    for k, a in obj['arrays'].items():
        raw = np.ascontiguousarray(a).tobytes()
        pad = (-offset) % 8                                # every array starts on an 8-byte boundary (typed-array views need it)
        blob += b'\0' * pad; offset += pad
        parts[k] = {'dtype': str(a.dtype), 'shape': list(a.shape), 'offset': offset, 'bytes': len(raw)}
        blob += raw; offset += len(raw)
    s, z = b64(blob)
    text = json.dumps({'gz': z, 'arrays': parts, **obj.get('meta', {}), 'b64': s}, separators=(',', ':'))
    path = os.path.join(out_dir, name.replace('/', '_') + '.js')
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'(window.LEARNED_PARTS=window.LEARNED_PARTS||{{}})["m/{name}"]={text};\n')
    return os.path.getsize(path)


def write_gen3(out_dir, classes):
    """m/gen3.js: for each class template, the same template with every unchanged leaf atom stripped (tools/template_general.py level 3), as
    {u: unique relaxed SMARTS, i: index into u per class, -1 = no relaxed form}. The page runs it next to the template as learned."""
    import template_general as TG
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    smarts = C.read_rows()['smarts']
    cache_path = os.path.join(TMP, 'organic_chem_gen3_cache.pkl')
    cache = pickle.load(open(cache_path, 'rb')) if os.path.exists(cache_path) else {}
    uniq, where, ids = [], {}, []
    for t in classes:
        t = int(t)
        if t not in cache: cache[t] = TG.generalise(smarts[t], 3) or ''
        g = cache[t]
        if not g: ids.append(-1); continue
        if g not in where: where[g] = len(uniq); uniq.append(g)
        ids.append(where[g])
    pickle.dump(cache, open(cache_path, 'wb'))
    from rel_fuse import write_json_part
    path = os.path.join(out_dir, 'gen3.js')
    write_json_part(path, {'u': uniq, 'i': ids})
    print(f'  relaxed templates: {sum(1 for i in ids if i >= 0):,} of {len(ids):,} classes have one ({len(uniq):,} different)')
    return os.path.getsize(path)


def load_model(tag):
    Z = np.load(os.path.join(TMP, f'organic_chem_relmodel_{tag}.npz'))
    nl = len([k for k in Z.files if k.startswith('W')])
    W = [Z[f'W{i}'].astype(np.float32) for i in range(nl)]
    b = [Z[f'b{i}'].astype(np.float32) for i in range(nl)]
    return W, b, Z['classes'], Z['cols']


def check(tag, n=6000):
    """does int8 change the model's answers on held-out reactions?"""
    W, b, classes, cols = load_model(tag)
    import eval_top1 as E
    X, y, split, meta, n_templates = T.load(E.model_feat(tag))
    test = np.flatnonzero(split == C.TEST)
    rng = np.random.default_rng(11)
    idx = np.sort(rng.permutation(test)[:n])
    cols_arg = None if (len(cols) == meta['layout']['bits'] and (cols == np.arange(len(cols))).all()) else cols
    Wq = []
    for w in W:
        q, s = quantise(w); Wq.append(dequantise(q, s))

    def run(Ws):
        net = T.Net.__new__(T.Net); net.W = Ws; net.b = b; net.sizes = [Ws[0].shape[0]] + [w.shape[1] for w in Ws]
        out = []
        for a in range(0, len(idx), 1024):
            out.append(net.probs(T.batch(X, idx[a:a + 1024], cols_arg)))
        return np.concatenate(out)
    P0, P1 = run(W), run(Wq)
    t0, t1 = P0[:, :-1].argmax(1), P1[:, :-1].argmax(1)
    top5_0 = np.argsort(-P0[:, :-1], axis=1)[:, :5]; top5_1 = np.argsort(-P1[:, :-1], axis=1)[:, :5]
    print(f'{n} held-out reactions: first-choice template unchanged by int8 for {(t0 == t1).mean():.2%}; same top-5 set for '
          f'{np.mean([set(a) == set(c) for a, c in zip(top5_0, top5_1)]):.2%}; mean |dp| of the first choice {np.abs(P0.max(1) - P1.max(1)).mean():.4f}')


def export(tag, out_dir, report=''):
    os.makedirs(out_dir, exist_ok=True)
    measured = json.load(open(report, encoding='utf-8')) if report else {}
    W, b, classes, cols = load_model(tag)
    import eval_top1 as E
    feat_path, meta_path = F.table_paths(E.model_feat(tag))
    fmeta = json.load(open(meta_path))
    sizes, total = {}, 0
    q1, s1 = quantise(W[0]); q2, s2 = quantise(W[1])
    sizes['w1'] = write_part(out_dir, 'w1', {'arrays': {'q': q1, 's': s1, 'b': b[0]}})
    sizes['w2'] = write_part(out_dir, 'w2', {'arrays': {'q': q2, 's': s2, 'b': b[1]}})
    n_cls = W[-1].shape[1]
    chunks = []
    for k, a in enumerate(range(0, n_cls, CHUNK_CLASSES)):
        sl = slice(a, min(n_cls, a + CHUNK_CLASSES))
        q3, s3 = quantise(W[-1][:, sl])
        sizes[f'w3_{k}'] = write_part(out_dir, f'w3_{k}', {'arrays': {'q': q3, 's': s3, 'b': b[-1][sl]}, 'meta': {'start': a, 'n': sl.stop - a}})
        chunks.append([a, sl.stop - a])
    sizes['classes'] = write_part(out_dir, 'classes', {'arrays': {'t': classes.astype(np.uint32)}})
    sizes['gen3'] = write_gen3(out_dir, classes)
    meta = {'layout': fmeta['layout'], 'fragments': fmeta['fragments'], 'roles': fmeta['roles'], 'hidden': [W[0].shape[1], W[1].shape[1]], 'classes': int(len(classes)),
            'chunks': chunks, 'input_cols': (None if (len(cols) == fmeta['layout']['bits'] and (cols == np.arange(len(cols))).all()) else cols.tolist()),
            'fp_bits': {'sub': F.FP_SUB, 'reag': F.FP_REAG, 'id': F.ID_BITS}, 'solvents': fmeta.get('solvents', []), 'model': tag, 'measured': measured or {}, 'built': time.strftime('%Y-%m-%d %H:%M')}
    text = json.dumps(meta, separators=(',', ':'))
    path = os.path.join(out_dir, 'meta.js')
    with open(path, 'w', encoding='utf-8') as f: f.write('window.RELEVANCE_META=' + text + ';\n')
    sizes['meta'] = os.path.getsize(path)
    print('wrote', out_dir, {k: f'{v / 1e6:.2f} MB' for k, v in sizes.items()}, f'total {sum(sizes.values()) / 1e6:.1f} MB')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True); ap.add_argument('--out', default=os.path.join(APP, 'data', 'learned', 'm'))
    ap.add_argument('--check', action='store_true')
    ap.add_argument('--report', default='', help='the measured figures written by tools/rel_fuse.py --report; stored in meta.js for the Validation tab')
    a = ap.parse_args()
    if a.check: check(a.model)
    else: export(a.model, a.out, a.report)
