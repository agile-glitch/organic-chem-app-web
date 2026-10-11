"""Find the reacting groups the Fragments sheet does not describe yet, from the recorded reactions.

    python tools/discover_fragments.py --n 60000                 # a sample of training reactions; writes tools/fragment_candidates.json

For each sampled reaction (training patents only) the reaction's own template is run on its starting materials (as the page runs
a learned template); the atoms of each starting material that MADE, BROKE or CHANGED a bond are the reaction centre. The Fragments
sheet (molecule_data.xlsx) is then asked which fragments cover that centre (tools/learn_reactions.reacted_fragments). A centre is
"described" when a SPECIFIC fragment covers it; the generic ones (an alkene, an arene C-H, a pyridine N, an enolisable C-H, an
alkyl arene, an azole N-H ...) match so many molecules that they say little about what reacted.
The per-substrate records (reaction, reagent, SMILES, changed atoms, fragments that reacted) are kept in the temp folder
(organic_chem_fragment_records.pkl): tools/fragment_eval.py scores candidate fragments against them.
Every centre that is not described is written down as its environment (the changed atoms and their neighbours, as a SMILES with the
changed atoms marked), and the environments are counted over all sampled reactions. The commonest are the groups worth adding.
"""
import argparse, collections, json, os, sys, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C

GENERIC = {'F005', 'F009', 'F038', 'F089', 'F091', 'F101', 'F035', 'F037', 'F102'}      # alkene, arene C-H, alkyl arene, azole N-H, pyridine N, enolisable C-H, ...


def _init(solv_frags):
    import learn_reactions as L
    L._init(*solv_frags, L.reagent_roles())


def centre_key(mol, atoms):
    """SMILES of the changed atoms and their neighbours, the changed ones marked with isotope 9; hydrogens counted on the changed atoms"""
    from rdkit import Chem
    m = Chem.RWMol(mol)
    env = set(atoms)
    for a in atoms:
        for nb in mol.GetAtomWithIdx(a).GetNeighbors(): env.add(nb.GetIdx())
    for a in atoms:
        at = m.GetAtomWithIdx(a); at.SetIsotope(9)
    try:
        return Chem.MolFragmentToSmiles(m, atomsToUse=sorted(env), canonical=True, isomericSmiles=True, allHsExplicit=False)
    except Exception:
        return ''


def work(task):
    import learn_reactions as L
    from rdkit import Chem
    ri, subs, reagent, prod, smarts = task
    extras = L.pool_extras(reagent, L._W['roles'])
    try:
        outs, leaving, changed = L.run_template(smarts, subs, extras, want=prod, detail=True)
    except Exception:
        return ri, []
    res = []
    for j, s in enumerate(subs):
        atoms = changed.get(j)
        m = Chem.MolFromSmiles(s)
        if m is None: continue
        if not atoms:
            res.append(([], '', s, []))                      # a starting material that did not react: kept so a candidate fragment's support can be counted
            continue
        fr = L.reacted_fragments(s, atoms)
        res.append((fr, centre_key(m, atoms), s, sorted(atoms)))
    if len(L._W['rxn']) > 2500:                      # both caches together (could_fit keys its fingerprints by id(rx))
        L._W['rxn'].clear(); L._W.get('qfp', {}).clear()
    return ri, res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=60000); ap.add_argument('--seed', type=int, default=5)
    ap.add_argument('--workers', type=int, default=0); ap.add_argument('--out', default=os.path.join(HERE, 'fragment_candidates.json'))
    a = ap.parse_args()
    t0 = time.time()
    import learn_reactions as L
    D = C.read_rows(); rows, smarts = D['rows'], D['smarts']
    rng = np.random.default_rng(a.seed)
    train = [i for i in range(len(rows)) if C.split_of_patent(rows[i][8]) == C.TRAIN]
    pick = sorted(rng.permutation(train)[:a.n].tolist())
    tasks = [(i, rows[i][1], rows[i][2], rows[i][3], smarts[rows[i][4]]) for i in pick]
    ids = {i: rows[i][0] for i in pick}
    reagents = {i: rows[i][2] for i in pick}
    del D, rows
    cover = collections.Counter(); sets = collections.Counter(); keys = collections.Counter(); ex = {}; n_centres = n_desc = n_reactions = 0
    by_frag = collections.Counter(); records = []
    with Pool(a.workers or max(1, (os.cpu_count() or 2) - 1), initializer=_init, initargs=((L.solvents(), L.fragments()),)) as pool:
        for k, (ri, res) in enumerate(pool.imap_unordered(work, tasks, chunksize=50), 1):
            n_reactions += 1
            for fr, key, s, atoms in res:
                records.append((ids[ri], reagents[ri], s, atoms, fr))
                if not atoms: continue
                n_centres += 1
                specific = [f for f in fr if f not in GENERIC]
                sets[' + '.join(fr) or '(none)'] += 1
                for f in fr: by_frag[f] += 1
                if specific: n_desc += 1; continue
                keys[key] += 1
                if key not in ex: ex[key] = (ids[ri], s, atoms, fr)
            if k % 10000 == 0: print(f'  {k:,} of {len(tasks):,} ({time.time() - t0:.0f} s)', flush=True)
    out = {'reactions': n_reactions, 'centres': n_centres, 'described_by_a_specific_fragment': round(n_desc / max(n_centres, 1), 4),
           'undescribed_centres': n_centres - n_desc, 'distinct_undescribed_environments': len(keys),
           'reacted_fragment_sets_top': sets.most_common(60), 'fragments_reacting': by_frag.most_common(),
           'top_environments': [{'env': k, 'count': c, 'example': ex[k]} for k, c in keys.most_common(400)]}
    json.dump(out, open(a.out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    import pickle, tempfile
    pickle.dump(records, open(os.path.join(tempfile.gettempdir(), 'organic_chem_fragment_records.pkl'), 'wb'), protocol=4)
    print(f'{n_centres:,} reaction centres in {n_reactions:,} reactions; {out["described_by_a_specific_fragment"]:.1%} described by a specific fragment; '
          f'{len(keys):,} distinct undescribed environments ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
