"""Held-out top-1 of the product predictor: how often is the recorded product the first one the page would show?

    python tools/eval_top1.py --model exp1 --split test --n 2000        # the relevance model's templates, run by RDKit
    python tools/eval_top1.py --model exp1 --split val --n 20000 --dump val_cands.pkl     # candidates kept, to fit the fusion weights

For every sampled reaction of the held-out patents (rel_common.py: 1 in 20 test, 1 in 20 validation; nothing learned has seen them)
the starting materials and the reagent are given, exactly what the Reactions page knows. The model proposes templates, each is run
on the molecules by RDKit the way the page runs a learned template (learn_reactions.run_template: every order, every site, the
carbon-free reagent pieces taking part), and the products are ranked. The prediction counts as right when the recorded product is
one of the molecules of the first-ranked outcome. A model with solvent bits (feature table v2) is scored with those bits CLEARED by default
(--solvents none): the page's conditions box is usually empty, so that is what a user gets; --solvents recorded gives it the recorded
solvent of each reaction (worth about 0.2 points of template top-1). Two readings, both reported:
  exact   the canonical SMILES with stereo must match
  flat    the structure without stereo must match (what the Validation tab and the older precedent figures use; RDKit's
          reaction runner does not carry R/S or E/Z through a template reliably)
'generated' is the ceiling: the recorded product is among the candidates at all (no ranking can fix a product that is not there).

Each candidate (a set of product molecules) keeps what the fusion step (tools/rel_fuse.py) needs: the model's probability mass
behind it, and its Morgan "difference fingerprint": the environments the reaction made (in the products, not in the starting
materials) and destroyed. The page can compute the same from the products it gets.
"""
import argparse, collections, json, os, pickle, sys, tempfile, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C
import rel_features as F
import rel_train as T
import template_general as TG

TMP = tempfile.gettempdir()
GEN_LEVELS = ()                                 # generalisation levels tried as extra candidate sources ((3,) = every unchanged leaf stripped); set by --gen
DIFF_BITS = 2048
MAX_DIFF = 200                                  # candidates (best by the model's mass) that get a difference fingerprint (js/relevance.js does the same)
_E = {}


def _init(solv_frags, frag_smarts=(), gen_levels=(), marks=False):
    global GEN_LEVELS
    GEN_LEVELS = tuple(gen_levels)                # (a worker process imports this module afresh, so the setting travels as an argument)
    import learn_reactions as L
    import rel_reactivity as RV
    L._init(*solv_frags, L.reagent_roles())
    if marks:                                     # interference flags (tools/interferences.py): the templates run with their changed atoms marked
        import interferences as I
        _E['lib'] = I.Library()
        _E['site_keys'] = _E['lib'].keys
    from rdkit.Chem import rdFingerprintGenerator as FG
    _E['gen'] = FG.GetMorganGenerator(radius=2, fpSize=DIFF_BITS)
    RV._init(list(frag_smarts))                 # the fragments, to see which one a candidate consumes (rel_reactivity.counts_of)


def clear_template_caches(L, limit=2500):
    """learn_reactions keeps every compiled template (unbounded) and, for the quick screen, the query fingerprints keyed by id(rx).
    Both must be emptied together: a freed rx object's id is reused, and a stale fingerprint would wrongly reject another template."""
    if len(L._W['rxn']) > limit:
        L._W['rxn'].clear()
        L._W.get('qfp', {}).clear()


def run_used(smarts, subs, extras, screen):
    """learn_reactions.run_template, but each distinct outcome also says how many of the STARTING MATERIALS (not the reagent pieces) its
    first producing order used: [(product SMILES list, number of starting materials used)]. js/reactions.js runLearned reports the same.
    With the interference library loaded (--marks) the template runs with its changed atoms marked and each outcome is (products, used, flags)."""
    if _E.get('lib') is not None:
        import int_run as R
        return [(cs, used, flags, '.'.join(mk)) for cs, used, mk, flags in R.run_marked(smarts, subs, extras, screen, _E['lib'])]
    import itertools
    import learn_reactions as L
    from rdkit import Chem
    from rdkit.Chem import AllChem
    rx = L._W['rxn'].get(smarts)
    if rx is None:
        try: rx = AllChem.ReactionFromSmarts(smarts); rx.Initialize()
        except Exception: rx = False
        L._W['rxn'][smarts] = rx
    if rx is False: return []
    n = rx.GetNumReactantTemplates()
    pool = list(subs) + [x for x in extras if x not in subs]
    if n > len(pool): return []
    mols = [Chem.MolFromSmiles(s) for s in pool]
    if screen is not None and not L.could_fit(rx, screen): return []
    outs, seen = [], set()
    for order in itertools.permutations(range(len(pool)), n):
        if min(order) >= len(subs): continue
        try: res = rx.RunReactants(tuple(mols[i] for i in order), 50)
        except Exception: continue
        used = len({i for i in order if i < len(subs)})
        for prods in res:
            cs = []
            for pm in prods:
                try: cs.append(L._canon(Chem.MolToSmiles(pm)))
                except Exception: cs.append(None)
            if not cs or any(c is None for c in cs): continue
            key = '.'.join(sorted(cs))
            if key not in seen: seen.add(key); outs.append((cs, used))
    return outs


def flat_smiles(smi):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    return Chem.MolToSmiles(m, isomericSmiles=False) if m is not None else smi


def union_bits(smiles_list):
    from rdkit import Chem
    v = np.zeros(DIFF_BITS, np.uint8)
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if m is not None: v |= np.asarray(_E['gen'].GetFingerprintAsNumPy(m), dtype=np.uint8)
    return v


def run_task(task):
    """task = (row index, starting materials, reagent, recorded product, [(template smarts, {source: weight})])
    returns (row index, {product set: {'flat': flat key, 'src': {source: [sum, max, n templates]}, 'hit': ..., 'hit_flat': ...,
    'formed': [bit], 'destroyed': [bit]}})"""
    import learn_reactions as L
    ri, subs, reagent, prod, tl = task
    extras = L.pool_extras(reagent, L._W['roles'])
    fps = L.pool_fps(subs, extras)
    prod_flat = flat_smiles(prod)
    cands = {}
    flat_cache = {}

    site_idx = {k: i for i, k in enumerate(_E['site_keys'])} if _E.get('lib') is not None else None

    def add(outs, ws, grp=0):
        for o in outs:
            cs, used = o[0], o[1]
            key = '.'.join(sorted(cs))
            c = cands.get(key)
            if c is None:
                fl = []
                for x in cs:
                    if x not in flat_cache: flat_cache[x] = flat_smiles(x)
                    fl.append(flat_cache[x])
                c = cands[key] = {'flat': '.'.join(sorted(fl)), 'src': {}, 'hit': prod in cs, 'hit_flat': prod_flat in fl, 'used': used}
                if site_idx is not None:                          # the interference flags of the first (most probable) template that gives this product, and which relaxed template that was
                    c['flags'] = sorted(site_idx[k] for k in o[2]); c['grp'] = grp; c['mk'] = o[3]        # mk: its products with the changed atoms marked (tools/int_reflag_dump.py re-reads the flags from it)
            c['used'] = max(c['used'], used)
            for s, w in ws.items():
                a = c['src'].setdefault(s, [0.0, 0.0, 0])
                a[0] += w / max(1, len(outs)); a[1] = max(a[1], w); a[2] += 1

    gcache = _E.setdefault('gcache', {})
    import zlib
    for smarts, ws in tl:
        grp = 0
        for level in GEN_LEVELS[:1] if site_idx is not None else ():       # the relaxed form identifies the group of candidates that are different sites of one template
            g0 = gcache.get((smarts, level), 0)
            if g0 == 0:
                g0 = gcache[(smarts, level)] = TG.generalise(smarts, level)
            grp = zlib.crc32((g0 or smarts).encode()) & 0x7FFFFFFF
        outs = run_used(smarts, subs, extras, fps)
        add(outs, ws, grp)
        for level in GEN_LEVELS:                              # the same template with less of its surroundings required (tools/template_general.py)
            g = gcache.get((smarts, level), 0)
            if g == 0:
                g = gcache[(smarts, level)] = TG.generalise(smarts, level)
            if g:
                gouts = run_used(g, subs, extras, fps)
                add(gouts, {f'g{level}': ws['mlp']}, grp)
    if len(gcache) > 30000: gcache.clear()
    from rdkit import Chem
    import rel_reactivity as RV
    extra = {'sub_heavy': sum(m.GetNumHeavyAtoms() for m in (Chem.MolFromSmiles(s) for s in subs) if m is not None), 'n_cands': len(cands),
             's_counts': RV.counts_of(subs).tolist() if RV._W.get('q') else []}
    if cands:
        sub_bits = union_bits(subs)
        best = sorted(cands, key=lambda k: -max((v[0] for v in cands[k]['src'].values()), default=0))[:MAX_DIFF]
        for k in best:
            parts = k.split('.')
            pb = union_bits(parts)
            cands[k]['formed'] = np.flatnonzero(pb & (1 - sub_bits)).astype(np.int16).tolist()
            cands[k]['destroyed'] = np.flatnonzero(sub_bits & (1 - pb)).astype(np.int16).tolist()
            cands[k]['n_mols'] = len(parts)
            cands[k]['heavy'] = sum(m.GetNumHeavyAtoms() for m in (Chem.MolFromSmiles(x) for x in parts) if m is not None)
            if RV._W.get('q'): cands[k]['counts'] = RV.counts_of(parts).tolist()
    clear_template_caches(L)
    return ri, cands, extra


def rank_by(cands, src, how='sum'):
    """candidates best first by one source's score; ties by the number of templates, then the text (deterministic)"""
    i = 0 if how == 'sum' else 1
    sc = [(c['src'][src][i], c['src'][src][2], k) for k, c in cands.items() if src in c['src']]
    sc.sort(key=lambda x: (-x[0], -x[1], x[2]))
    return [k for _, _, k in sc]


def flat_rank(cands, order):
    """the same ranking with stereo variants of one structure merged (the best one speaks for the group): flat keys best first"""
    seen, out = set(), []
    for k in order:
        f = cands[k]['flat']
        if f not in seen: seen.add(f); out.append(k)
    return out


class Tally:
    def __init__(self): self.n = 0; self.a = collections.Counter()
    def add(self, cands, order):
        self.n += 1
        if not order: return
        hit = [cands[k]['hit'] for k in order]
        fo = flat_rank(cands, order)
        hff = [cands[k]['hit_flat'] for k in fo]
        self.a['answered'] += 1
        for name, h in (('exact', hit), ('flat', hff)):
            if h[:1] == [True]: self.a[name + '_top1'] += 1
            if any(h[:3]): self.a[name + '_top3'] += 1
            if any(h[:5]): self.a[name + '_top5'] += 1
        if any(hit): self.a['exact_generated'] += 1
        if any(cands[k]['hit_flat'] for k in order): self.a['flat_generated'] += 1
    def out(self):
        d = max(self.n, 1)
        return {'n': self.n, **{k: round(v / d, 4) for k, v in sorted(self.a.items())}}


def sample_rows(split_code, n, seed, feat='v1'):
    X, y, split, meta, n_templates = T.load(feat)
    D = C.read_rows()
    idx = np.flatnonzero(split == split_code)
    rng = np.random.default_rng(seed)
    return np.sort(rng.permutation(idx)[:n]), X, y, split, meta, D


def solvent_mask(meta, solvents):
    """a function that clears the solvent bits of a batch (unpacked, all bits) unless solvents == 'recorded'"""
    L = meta['layout']
    if solvents == 'recorded' or 'solv' not in L: return lambda xb: xb
    a, w = L['solv']
    def clear(xb):
        xb[:, a:a + w] = 0
        return xb
    return clear


def model_feat(tag):
    """the feature table (tag) a model was trained on; models saved before this was recorded used 'v1'"""
    Z = np.load(os.path.join(TMP, f'organic_chem_relmodel_{tag}.npz'))
    return str(Z['feat']) if 'feat' in Z.files else 'v1'


def load_net(tag, meta):
    Z = np.load(os.path.join(TMP, f'organic_chem_relmodel_{tag}.npz'))
    nl = len([k for k in Z.files if k.startswith('W')])
    net = T.Net.__new__(T.Net); net.W = [Z[f'W{i}'] for i in range(nl)]; net.b = [Z[f'b{i}'] for i in range(nl)]
    net.sizes = [net.W[0].shape[0]] + [w.shape[1] for w in net.W]
    cols = Z['cols']
    cols = None if len(cols) == meta['layout']['bits'] and (cols == np.arange(len(cols))).all() else cols
    return net, Z['classes'], cols


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', required=True); ap.add_argument('--feat', default='')
    ap.add_argument('--split', choices=['test', 'val'], default='test'); ap.add_argument('--n', type=int, default=2000)
    ap.add_argument('--topn', type=int, default=150, help='templates tried per reaction (the model\'s most probable)')
    ap.add_argument('--workers', type=int, default=0); ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--dump', default='', help='keep every reaction\'s candidates in this file (temp folder), to fit/score later')
    ap.add_argument('--gen', default='', help='also try these generalisation levels of each template, e.g. 3 or 2,3 (extra sources g2, g3)')
    ap.add_argument('--solvents', choices=['none', 'recorded'], default='none', help='for a model with solvent bits: clear them (what a user with an empty conditions box gets) or give the recorded solvents')
    ap.add_argument('--marks', action='store_true', help='run the templates with their changed atoms marked and keep the interference flags of each candidate (tools/interferences.py)')
    ap.add_argument('--note', default='')
    a = ap.parse_args()
    t0 = time.time()
    global GEN_LEVELS
    GEN_LEVELS = tuple(int(x) for x in a.gen.split(',') if x)
    code = C.TEST if a.split == 'test' else C.VAL
    a.feat = a.feat or model_feat(a.model)
    sel, X, y, split, meta, D = sample_rows(code, a.n, a.seed, a.feat)
    rows, smarts = D['rows'], D['smarts']
    net, classes, cols = load_net(a.model, meta)
    mask = solvent_mask(meta, a.solvents)
    import rel_reactivity as RV
    roles_by_smiles = RV.reagent_roles()
    frs = RV.fragment_list()
    # how many TRAINING reactions used this exact reagent with each template (the page's data/learned/g/ shards hold the same counts)
    by_reagent = collections.defaultdict(dict)
    for i, r in enumerate(rows):
        if split[i] == C.TRAIN and r[2] and not r[5].startswith('?'):
            d = by_reagent[r[5]]; d[r[4]] = d.get(r[4], 0) + 1
    tasks, in_topn, info = [], 0, {}
    for s in range(0, len(sel), 1024):
        j = sel[s:s + 1024]
        Pall = net.probs(T.batch(X, j, cols, mask=mask))
        P = Pall[:, :-1]
        top = np.argsort(-P, axis=1)[:, :a.topn]
        for q, ri in enumerate(j):
            r = rows[ri]
            rk = r[5] if r[2] and not r[5].startswith('?') else ''
            cnts = by_reagent.get(rk, {})
            tl = [(smarts[classes[c]], dict({'mlp': float(P[q, c])}, **({'rc': float(cnts[classes[c]])} if cnts.get(classes[c]) else {}))) for c in top[q]]
            tasks.append((int(ri), r[1], r[2], r[3], tl))
            in_topn += int(any(classes[c] == y[ri] for c in top[q]))
            info[int(ri)] = {'p_other': float(Pall[q, -1]), 'mass': float(P[q, top[q]].sum()), 'top_p': float(P[q, top[q][0]]), 'has_reagent': bool(r[2]), 'n_subs': len(r[1]),
                             'template_in_topn': bool(any(classes[c] == y[ri] for c in top[q])),
                             'exact_key': rk, 'rc_known': bool(cnts), 'role_sig': RV.role_signature(r[2], roles_by_smiles)}
    print(f'{len(tasks):,} {a.split} reactions; the recorded template is among the model\'s top {a.topn} for {in_topn / len(tasks):.1%} ({time.time() - t0:.0f} s)', flush=True)
    import learn_reactions as L
    workers = a.workers or max(1, (os.cpu_count() or 2) - 1)
    res = {}
    with Pool(workers, initializer=_init, initargs=((L.solvents(), L.fragments()), [f[2] for f in frs], GEN_LEVELS, a.marks)) as pool:
        for k, (ri, cands, extra) in enumerate(pool.imap_unordered(run_task, tasks, chunksize=4), 1):
            res[ri] = cands; info[ri].update(extra)
            if k % 250 == 0: print(f'  {k:,} of {len(tasks):,} ({time.time() - t0:.0f} s)', flush=True)
    if a.dump:
        site_keys = None
        if a.marks:
            import interference_defs as ID
            site_keys = [x['key'] for x in ID.SITE]
        pickle.dump({'res': res, 'info': info, 'sel': sel.tolist(), 'model': a.model, 'topn': a.topn, 'fragments': [f[0] for f in frs], 'site_keys': site_keys}, open(os.path.join(TMP, a.dump), 'wb'), protocol=4)
    out = {'model': a.model, 'split': a.split, 'n': len(tasks), 'topn': a.topn, 'solvents': a.solvents, 'recorded_template_in_topn': round(in_topn / len(tasks), 4), 'note': a.note, 'when': time.strftime('%Y-%m-%d %H:%M')}
    groups = {'all': lambda r: True, 'with a reagent': lambda r: bool(r[2]), 'no reagent': lambda r: not r[2]}
    for how in ('sum', 'max'):
        for gname, keep in groups.items():
            t = Tally()
            for ri, cands in res.items():
                if keep(rows[ri]): t.add(cands, rank_by(cands, 'mlp', how))
            out[f'{gname} / {how}'] = t.out()
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(TMP, f'organic_chem_eval_{a.model}_{a.split}.json'), 'w'), indent=1)
    print(f'({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
