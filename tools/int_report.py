"""Put the interference layer's measured figures into tools/eval_top1.json (key "interference"), which tools/rel_export.py stores in data/learned/m/meta.js and the
Validation tab shows.

    python tools/int_report.py --int-model v4 --ctx-model new2 --paired ana_base.pkl,ana_both.pkl --runs base=relwork/abl_base.json,site=relwork/abl_site.json,context=relwork/abl_ctx.json,both=relwork/abl_both.json

Run it after tools/rel_fuse.py --report has written tools/eval_top1.json. Nothing here is computed on the test patents except what the --runs files hold
(tools/rel_fuse.py --results-json: each ranker was fitted on the validation patents and scored on the test patents once).
"""
import argparse, json, os, pickle, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import interference_defs as D

TMP = tempfile.gettempdir()


def paired_difference(path_a, path_b, n_boot=4000, seed=11):
    """first-answer accuracy of ranker B minus ranker A on the SAME test reactions, with a 95 % interval from resampling whole patents (the reactions of one patent are alike).
    gained / lost = reactions B gets right and A wrong / the other way round"""
    import numpy as np
    A = {x['ri']: x for x in pickle.load(open(os.path.join(TMP, path_a), 'rb'))}
    B = {x['ri']: x for x in pickle.load(open(os.path.join(TMP, path_b), 'rb'))}
    ris = sorted(set(A) & set(B))
    pat = np.array([A[r]['patent'] for r in ris])
    uniq, inv = np.unique(pat, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uniq)).astype(float)
    rng = np.random.default_rng(seed)
    picks = rng.integers(0, len(uniq), size=(n_boot, len(uniq)))
    out = {'n': len(ris), 'patents': int(len(uniq))}
    for field, label in (('hit_exact', 'top1_exact'), ('hit', 'top1_flat')):
        ya = np.array([1.0 if A[r]['cands'] and A[r]['cands'][0].get(field) else 0.0 for r in ris])
        yb = np.array([1.0 if B[r]['cands'] and B[r]['cands'][0].get(field) else 0.0 for r in ris])
        d = yb - ya
        sums = np.bincount(inv, weights=d, minlength=len(uniq))
        boot = sums[picks].sum(1) / cnt[picks].sum(1)
        lo, hi = np.percentile(boot, [2.5, 97.5])
        out[label] = {'diff': round(float(d.mean()), 4), 'ci95': [round(float(lo), 4), round(float(hi), 4)], 'gained': int((d > 0).sum()), 'lost': int((d < 0).sum()),
                      'share_of_resamples_above_zero': round(float((boot > 0).mean()), 3)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--int-model', default='v4'); ap.add_argument('--ctx-model', default='new2')
    ap.add_argument('--runs', required=True, help='name=path of rel_fuse --results-json files, comma separated; names: base, site, context, both')
    ap.add_argument('--evidence', default=os.path.join(HERE, 'interference_evidence.json'))
    ap.add_argument('--report', default=os.path.join(HERE, 'eval_top1.json'))
    ap.add_argument('--note', default='')
    ap.add_argument('--paired', default='', help='base_analysis.pkl,other_analysis.pkl (rel_fuse --analysis of two runs on the same test reactions): the paired difference with a patent-level interval')
    a = ap.parse_args()
    import int_model as IM
    im = IM.IntModel(a.int_model)
    ctx = pickle.load(open(os.path.join(TMP, f'organic_chem_ctx_{a.ctx_model}.pkl'), 'rb'))
    rep = json.load(open(a.report, encoding='utf-8'))
    m = im.m
    out = {'built': time.strftime('%Y-%m-%d'), 'patterns': len(D.SITE), 'advisories': len(D.ADVISORY), 'families': len(D.FAMILIES),
           'contests': {'train': m['train_contests']},
           'site_choice': {'random': round(m['val']['random'], 4), 'global': round(m['val']['global'], 4), 'learned': round(m['val']['families'], 4), 'n': m['val']['n']},
           'chemoselectivity': {'rate_only': round(ctx['val_logloss']['rate_only'], 4), 'with_context': round(ctx['val_logloss']['with_context'], 4), 'pairs': ctx['val_logloss']['pairs']},
           'note': a.note}
    # Hammett reaction constants learned per family
    import numpy as np
    out['rho'] = {'global': [round(float(x), 2) for x in im.g[im.nf:]], **{f: [round(float(x), 2) for x in im.eff[f][im.nf:]] for f in im.fams if f in ('aryl_halide', 'arene_ch', 'acyl_substitution', 'alkyl_substitution')}}
    top = {}
    errors = {}
    for item in a.runs.split(','):
        name, _, path = item.partition('=')
        r = json.load(open(path, encoding='utf-8'))
        top[name] = {'top1_exact': r['test']['top1_exact'], 'top1_flat': r['test']['top1_flat'], 'validation_oof_exact': r['validation_oof']['top1_exact'], 'features': r['features'][20:]}
        errors[name] = r['test_errors']
    out['top1'], out['errors'] = top, errors
    if a.paired:
        pa, pb = a.paired.split(',')
        out['paired'] = paired_difference(pa, pb)
        print('paired difference (second minus first):', json.dumps(out['paired']))
    if os.path.exists(a.evidence):
        ev = json.load(open(a.evidence, encoding='utf-8'))
        out['advisory_evidence'] = {k: {x: v[x] for x in ('reactions_with_reagent_and_a_group', 'consumed_with', 'consumed_other', 'would_be_shown') if x in v} for k, v in ev.items() if v.get('measured')}
    rep['interference'] = out
    json.dump(rep, open(a.report, 'w', encoding='utf-8'), indent=1)
    print('wrote the interference figures to', a.report)


if __name__ == '__main__':
    main()
