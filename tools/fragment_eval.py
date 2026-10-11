"""Score candidate fragments before they are added to the Fragments sheet.

    python tools/fragment_eval.py                      # tools/new_fragments.py CANDIDATES against the recorded reactions

For each candidate (name + SMARTS, tools/new_fragments.py):
  checks      the SMARTS compiles, matches every positive example and none of the negative ones
  support     starting materials of the sampled reactions (tools/discover_fragments.py) it matches
  reacted     ... of which an atom of the match changed in the recorded reaction (it took part)
  new         ... of which no SPECIFIC fragment of the current sheet described the reaction centre (the information it adds)
  reagent     learned reactions whose reagent has a part it matches (for reagent-type fragments)
and prints the most common reagents of the reactions where it reacted.
"""
import collections, json, os, pickle, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C
import discover_fragments as DF
import new_fragments as NF
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')


def main():
    t0 = time.time()
    recs = pickle.load(open(os.path.join(tempfile.gettempdir(), 'organic_chem_fragment_records.pkl'), 'rb'))
    mols = {}
    def mol(s):
        if s not in mols: mols[s] = Chem.MolFromSmiles(s)
        return mols[s]
    D = C.read_rows(); rows = D['rows']
    parts = collections.Counter()
    for r in rows:
        for p in set(x for x in r[2].split('.') if x): parts[p] += 1
    print(f'{len(recs):,} sampled starting materials; {len(parts):,} distinct reagent parts ({time.time() - t0:.0f} s)', flush=True)
    out = []
    for c in NF.CANDIDATES:
        q = Chem.MolFromSmarts(c['smarts'])
        row = {'slug': c['slug'], 'name': c['name'], 'smarts': c['smarts'], 'kind': c['kind']}
        if q is None:
            row['error'] = 'SMARTS does not compile'; out.append(row); print('!!', c['slug'], row['error']); continue
        bad = []
        for s in c.get('pos', []):
            m = Chem.MolFromSmiles(s)
            if m is None or not m.HasSubstructMatch(q): bad.append('misses ' + s)
        for s in c.get('neg', []):
            m = Chem.MolFromSmiles(s)
            if m is not None and m.HasSubstructMatch(q): bad.append('wrongly matches ' + s)
        row['check_problems'] = bad
        sup = rea = new = 0
        reag = collections.Counter(); ex = []
        for rid, reagent, s, atoms, fr in recs:
            m = mol(s)
            if m is None: continue
            ms = m.GetSubstructMatches(q, uniquify=True, maxMatches=200)
            if not ms: continue
            sup += 1
            if atoms and any(set(atoms) & set(mm) for mm in ms):
                rea += 1
                reag[reagent or '(none)'] += 1
                if not [f for f in fr if f not in DF.GENERIC]: new += 1
                if len(ex) < 2: ex.append(s)
        n_r = 0
        for p, n in parts.items():
            m = mol(p)
            if m is not None and m.HasSubstructMatch(q): n_r += n
        row.update({'support': sup, 'reacted': rea, 'new': new, 'reagent_reactions': n_r, 'top_reagents': reag.most_common(3), 'examples': ex})
        out.append(row)
        flag = ' PROBLEMS: ' + '; '.join(bad) if bad else ''
        print(f"{c['slug']:22s} support {sup:5d} reacted {rea:5d} new {new:4d} as reagent in {n_r:7,d} reactions{flag}", flush=True)
    json.dump(out, open(os.path.join(HERE, 'fragment_candidates_scored.json'), 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print(f'({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
