"""What the PYTHON computes for some held-out reactions, to compare with the page (js/relevance.js): the network's input bits and its answer.

    python tools/rel_check_js.py --model full2 --n 40 --out tools/_parity_bits.json
    python tools/rel_check_js.py --model full2 --n 40 --conds "toluene, water" --out tools/_parity_bits_solv.json     # the conditions box says this

Then in the browser console (Reactions tab); the page must give the same input bits and, up to int8 rounding of the weights, the same most probable templates:

    const P = await (await fetch('tools/_parity_bits.json')).json(); let bitsOk = 0, topOk = 0;
    for (const p of P.reactions) {
      const o = await ReactionsRelevance.inspect(p.subs, p.reagent, 5, P.conds);
      if (JSON.stringify(o.bits) === JSON.stringify(p.bits)) bitsOk++;
      if (o.top[0][0] === p.top[0][0]) topOk++;
    }
    console.log(bitsOk, 'of', P.reactions.length, 'input vectors identical;', topOk, 'first-choice templates equal');

(delete the json afterwards; whole predictions are compared with rel_fuse.py --parity, see tools/README_relevance.md)
"""
import argparse, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C, rel_features as F, rel_train as T, rel_export as X, eval_top1 as E


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True); ap.add_argument('--n', type=int, default=40); ap.add_argument('--seed', type=int, default=21)
    ap.add_argument('--conds', default='', help="the conditions box text; its solvent names set the solvent bits ('toluene, water')")
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    D = C.read_rows(); rows = D['rows']
    feat = E.model_feat(a.model)
    meta = json.load(open(F.table_paths(feat)[1])); layout = meta['layout']
    fr = [f for f in F.fragment_list() if f[0] in meta['fragments']]
    assert [f[0] for f in fr] == meta['fragments'], 'the Fragments sheet no longer has the fragments this model was trained with'
    roles, by_smiles = F.role_list()
    F._init([f[2] for f in fr], by_smiles, roles, meta.get('solvents', []))
    W, b, classes, cols = X.load_model(a.model)
    net = T.Net.__new__(T.Net); net.W = W; net.b = b; net.sizes = [W[0].shape[0]] + [w.shape[1] for w in W]
    cols = None if (len(cols) == layout['bits'] and (cols == np.arange(len(cols))).all()) else cols
    test = [i for i in range(len(rows)) if C.split_of_patent(rows[i][8]) == C.TEST]
    pick = sorted(np.random.default_rng(a.seed).permutation(test)[:a.n].tolist())
    conds = ', '.join(x.strip() for x in a.conds.split(',') if x.strip())      # the solvent names as the page reads them (Solvents sheet names, ', ' between)
    res = []
    for ri in pick:
        r = rows[ri]
        packed = F.featurise(r[1], r[2], layout, conds)
        bits = np.flatnonzero(np.unpackbits(packed)[:layout['bits']]).tolist()
        x = np.zeros((1, layout['bits']), np.float32); x[0, bits] = 1
        p = net.probs(x[:, cols] if cols is not None else x)[0]
        top = np.argsort(-p[:-1])[:5]
        res.append({'id': r[0], 'subs': r[1], 'reagent': r[2], 'bits': bits, 'top': [[int(classes[c]), float(p[c])] for c in top], 'p_other': float(p[-1]), 'template': int(r[4])})
    json.dump({'model': a.model, 'conds': conds, 'reactions': res}, open(a.out, 'w'))
    print(len(res), 'reactions;', sum(len(x['bits']) for x in res) // len(res), 'input bits each on average ->', a.out)


if __name__ == '__main__':
    main()
