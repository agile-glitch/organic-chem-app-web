"""What the PYTHON gives for the interference layer, to compare with the page (js/interference.js).

    python tools/int_check_js.py --model full2 --n 1500 --dump val_full2m2.pkl --n-real 3000 --out tools/_parity_int.json

Then in the browser console of the Reactions tab (Interference must be loaded; data/reaction_rules.js must hold the Interferences sheet):

    const P = await (await fetch('tools/_parity_int.json')).json(); const R = await ReactionsRelevance.rdkit(); await Interference.load();
    let mk = 0, fl = 0; const bad = [];
    for (const t of P.marking) if (Interference.markTemplate(t.smarts, t.maps) === t.marked) mk++; else bad.push(['mark', t.smarts]);
    for (const c of P.flags) if (JSON.stringify(Interference.flags(R, [c.smiles])) === JSON.stringify(c.flags)) fl++; else bad.push(['flags', c.smiles, c.flags, Interference.flags(R, [c.smiles])]);
    let rl = 0; for (const c of P.real || []) if (JSON.stringify(Interference.flags(R, [c.smiles])) === JSON.stringify(c.flags)) rl++; else bad.push(['real', c.smiles, c.flags, Interference.flags(R, [c.smiles])]);
    console.log(mk, 'of', P.marking.length, 'marked templates identical;', fl, 'of', P.flags.length, 'flag lists identical;', rl, 'of', (P.real || []).length, 'real candidate products identical', bad.slice(0, 5));

(delete the json afterwards)
"""
import argparse, json, os, sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import int_marks as M
import interferences as I
import interference_tests as T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='full2'); ap.add_argument('--n', type=int, default=1500); ap.add_argument('--seed', type=int, default=5)
    ap.add_argument('--out', required=True)
    ap.add_argument('--dump', default='', help='a candidate dump of tools/eval_top1.py --marks (in the temp folder): also take real marked products from it')
    ap.add_argument('--n-real', type=int, default=3000)
    a = ap.parse_args()
    import tempfile
    Z = np.load(os.path.join(tempfile.gettempdir(), f'organic_chem_relmodel_{a.model}.npz'))
    import int_export as X
    smarts = X.templates_smarts()
    rng = np.random.default_rng(a.seed)
    pick = rng.choice(Z['classes'], size=a.n, replace=False)
    marking = []
    for t in pick:
        s = smarts[int(t)]
        maps = sorted(M.changed_maps(s))
        marking.append({'smarts': s, 'maps': maps, 'marked': M.mark_template(s, only=set(maps))})
    lib = I.Library()
    flags = [{'smiles': smi, 'flags': lib.flags_of_smiles(smi)} for smi, _y, _n in T.all_cases()]
    real = []
    if a.dump:
        import pickle
        d = pickle.load(open(os.path.join(tempfile.gettempdir(), a.dump), 'rb'))
        mks = sorted({c['mk'] for r in d['res'].values() for c in r.values() if c.get('mk')})
        rng2 = np.random.default_rng(a.seed + 1)
        for i in rng2.choice(len(mks), size=min(a.n_real, len(mks)), replace=False):
            real.append({'smiles': mks[int(i)], 'flags': lib.flags_of_smiles(mks[int(i)])})
    json.dump({'marking': marking, 'flags': flags, 'real': real}, open(a.out, 'w'))
    print(len(marking), 'templates,', len(flags), 'marked molecules,', len(real), 'real candidate products ->', a.out)


if __name__ == '__main__':
    main()
