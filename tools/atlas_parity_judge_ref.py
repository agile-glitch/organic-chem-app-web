"""Fragment Atlas: the reference side of the judgement parity test (table lookups, backoff chain, combining model).

    python tools/atlas_parity_judge_ref.py <instances.npz> <atlas_eval.json>     writes tools/atlas_parity_judge.json   (then: node tools/atlas_parity.js --judge)

Picks random fragment instances, computes with numpy exactly what the page must compute from the shipped tables (counts from ALL instances, pruned at the
minimum of each family) and writes the instance, its reagent side and the expected probability. tools/atlas_parity.js feeds the same instance to
Atlas.features / Atlas.combine in the browser, which reads data/atlas/ through its own shard decoder and its own hashing.
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atlas_keys as K
import atlas_eval as E


def main():
    d = dict(np.load(sys.argv[1], allow_pickle=True))
    ev = json.load(open(sys.argv[2], encoding='utf-8'))
    n = len(d['y'])
    split = d['split_mol'] if ev.get('split') == 'molecule' else d['split']
    P = K.prepare(d)
    rng = np.random.default_rng(11)
    rows = np.zeros(n, bool); rows[rng.choice(n, 400, replace=False)] = True
    # the page's tables come from every instance, with the training base rate as the prior (as tools/atlas_build.py writes them)
    d2 = dict(d); d2['y'] = d['y']
    names, grp, X, ridx = E.build_features(d, P, np.ones(n, bool), rows, K.M_CHAIN, K.MIN_N)
    train_base_all = float(d['y'].mean()); train_base = float(d['y'][split == 2].mean())
    # build_features used the all-instance mean as the prior; the page uses the training mean: recompute the expected p with the page's prior
    names, grp, X, ridx = E.build_features_with_base(d, P, np.ones(n, bool), rows, train_base, K.M_CHAIN, K.MIN_N)
    M = ev['model']
    assert names == M['names'], 'the model was fitted with other columns'
    w, mu, sd = np.array(M['coef']), np.array(M['mean']), np.array(M['std'])
    p = 1 / (1 + np.exp(-np.clip(w[0] + ((X - mu) / sd) @ w[1:], -30, 30)))
    cls_list, parts, roles = list(d['classes']), list(d['parts']), list(d['part_role'])
    cases = []
    for i, r in enumerate(ridx.tolist()):
        side_parts = [parts[int(d[k][r]) - 1] for k in ('p1', 'p2', 'p3') if int(d[k][r]) > 0]
        q = float(d['charge'][r])
        cases.append({'F': int(d['frag'][r]), 'e': [int(d['e0'][r]), int(d['e1'][r]), int(d['e2'][r]), int(d['e3'][r])], 'steric': int(d['steric'][r]),
                      'charge': None if q != q else q, 'ring': int(d['ring'][r]), 'strain': float(d['strain'][r]), 'hyb': int(d['hyb'][r]), 'nAtoms': int(d['natoms'][r]),
                      'cls': cls_list[int(d['cls'][r])], 'parts': side_parts, 'p': float(p[i]), 'fid': str(d['frag_ids'][int(d['frag'][r])])})
    out = os.path.join(HERE, 'atlas_parity_judge.json')
    json.dump({'about': 'expected probabilities computed by tools/atlas_eval.py build_features + the fitted model; checked by tools/atlas_parity.js --judge', 'cases': cases}, open(out, 'w'), separators=(',', ':'))
    print(len(cases), 'instances ->', out, f'(mean p {np.mean(p):.3f})')


if __name__ == '__main__':
    main()
