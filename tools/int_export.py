"""Write the interference data for the page: data/learned/i/*.js  (loaded on first use by js/interference.js)

    python tools/int_export.py --model full2 --int-model v1 --out data/learned/i            # maps + weights
    python tools/int_export.py --maps-only ... / --weights-only ...

Files (each fills window.LEARNED_PARTS["i/<name>"] with gzipped base64 JSON, as data/learned/m/):
  maps      per class template of the relevance model ('c') and per relaxed form ('g'), the atom-map numbers of the atoms the template changes
            ("1,3"), which the page needs to mark the product side of the template (tools/int_marks.py changed_maps; the page cannot work it out: it has
            no reaction-template parser)
  weights   the learned weight of each site pattern (global + per family), the Hammett coefficients, the family definitions, the figures measured on
            validation contests (tools/int_train.py)
The patterns themselves, and the advisories, are the sheet "Interferences" of reaction_rules.xlsx (tools/rules_to_js.py puts them in data/reaction_rules.js).
"""
import argparse, base64, gzip, json, os, pickle, sys, tempfile, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
sys.path.insert(0, HERE)

TMP = tempfile.gettempdir()


def write_part(path, name, obj):
    gz = base64.b64encode(gzip.compress(json.dumps(obj, separators=(',', ':')).encode('utf-8'), 9)).decode('ascii')
    with open(path, 'w', encoding='utf-8') as f:
        f.write(f'(window.LEARNED_PARTS=window.LEARNED_PARTS||{{}})["i/{name}"]={{"gz":"{gz}"}};\n')
    return os.path.getsize(path)


def templates_smarts():
    """the SMARTS of every learned template by template number, read from the page's own shards (data/learned/t/*.js): 0.1 GB instead of the 1.5 GB of the learned rows"""
    import re
    d = os.path.join(APP, 'data', 'learned')
    meta = open(os.path.join(d, 'learned_meta.js'), encoding='utf-8').read()
    chunk = int(re.search(r'"tmpl_chunk":\s*(\d+)', meta).group(1))
    n = int(re.search(r'"n_templates":\s*(\d+)', meta).group(1))
    out = []
    for c in range((n + chunk - 1) // chunk):
        txt = open(os.path.join(d, 't', f'{c:03d}.js'), encoding='utf-8').read()
        j = json.loads(gzip.decompress(base64.b64decode(re.search(r'"gz":"([^"]+)"', txt).group(1))))
        out.extend(j['smarts'])
    return out


def relaxed_forms(classes):
    """the unique relaxed templates in the order data/learned/m/gen3.js lists them (rel_export.write_gen3) and each class's index into them (-1 = none)"""
    smarts = templates_smarts()
    cache = pickle.load(open(os.path.join(TMP, 'organic_chem_gen3_cache.pkl'), 'rb'))
    uniq, where, ids = [], {}, []
    import template_general as TG
    for t in classes:
        t = int(t)
        if t not in cache: cache[t] = TG.generalise(smarts[t], 3) or ''
        g = cache[t]
        if not g: ids.append(-1); continue
        if g not in where: where[g] = len(uniq); uniq.append(g)
        ids.append(where[g])
    return smarts, uniq, ids


def export_maps(model, out_dir):
    import int_marks as M
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    Z = np.load(os.path.join(TMP, f'organic_chem_relmodel_{model}.npz'))
    classes = Z['classes']
    smarts, uniq, ids = relaxed_forms(classes)
    t0 = time.time()
    c = []
    for k, t in enumerate(classes):
        c.append(','.join(str(x) for x in sorted(M.changed_maps(smarts[int(t)]))))
        if (k + 1) % 10000 == 0: print(f'  class templates {k + 1:,} of {len(classes):,} ({time.time() - t0:.0f} s)', flush=True)
    g = []
    for k, u in enumerate(uniq):
        g.append(','.join(str(x) for x in sorted(M.changed_maps(u))))
        if (k + 1) % 10000 == 0: print(f'  relaxed forms {k + 1:,} of {len(uniq):,} ({time.time() - t0:.0f} s)', flush=True)
    n_empty = sum(1 for x in c if not x), sum(1 for x in g if not x)
    print(f'maps: {len(c):,} class templates ({n_empty[0]:,} without a changed atom RDKit could find), {len(g):,} relaxed forms ({n_empty[1]:,} without)')
    size = write_part(os.path.join(out_dir, 'maps.js'), 'maps', {'c': c, 'g': g, 'built': time.strftime('%Y-%m-%d %H:%M')})
    print(f'wrote {out_dir}/maps.js ({size / 1e6:.2f} MB)')


def export_weights(tag, out_dir):
    import interference_defs as D
    import int_model as IM
    m = IM.IntModel(tag)
    nf = m.nf
    obj = {'keys': m.keys, 'families': m.fams, 'family_fragments': [[n, ids] for n, _, ids in D.FAMILIES],
           'g': [round(float(x), 4) for x in m.g[:nf]], 'fam': {n: [round(float(x), 4) for x in m.d[n][:nf]] for n in m.fams},
           'rho_g': [round(float(x), 4) for x in m.g[nf:]], 'rho_f': {n: [round(float(x), 4) for x in m.d[n][nf:]] for n in m.fams},
           'sigma': {s['key']: [{'ortho': 0, 'meta': 1, 'para': 2}[s['key'].split('-')[1]], s['sigma']] for s in D.SITE if s['sigma'] is not None},
           'val': m.m['val'], 'penalties': m.m['penalties'], 'train_contests': m.m['train_contests'], 'built': m.m['built']}
    size = write_part(os.path.join(out_dir, 'weights.js'), 'weights', obj)
    print(f'wrote {out_dir}/weights.js ({size / 1e3:.0f} KB): {nf} patterns x {len(m.fams)} families')


def export_ctx(tag, out_dir):
    import int_model as IM
    m = IM.CtxModel(tag)
    nf = len(m.fragments)
    beta = [[[int(g), round(float(m.beta[f, g]), 4)] for g in np.flatnonzero(m.beta[f])] for f in range(nf)]
    obj = {'fragments': m.fragments, 'b': [round(float(x), 4) for x in m.b], 'mult': [round(float(x), 4) for x in m.mult], 'beta': beta,
           'val_logloss': m.m['val_logloss'], 'built': m.m['built']}
    size = write_part(os.path.join(out_dir, 'ctx.js'), 'ctx', obj)
    print(f'wrote {out_dir}/ctx.js ({size / 1e3:.0f} KB): {sum(len(r) for r in beta):,} coefficients over {nf} groups')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='full2'); ap.add_argument('--int-model', default='v4'); ap.add_argument('--out', default=os.path.join(APP, 'data', 'learned', 'i'))
    ap.add_argument('--ctx-model', default='new2'); ap.add_argument('--maps-only', action='store_true'); ap.add_argument('--weights-only', action='store_true')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    if not a.weights_only: export_maps(a.model, a.out)
    if not a.maps_only: export_weights(a.int_model, a.out); export_ctx(a.ctx_model, a.out)
