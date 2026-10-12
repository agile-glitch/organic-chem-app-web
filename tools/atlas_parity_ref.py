"""Fragment Atlas: the reference side of the Python / browser parity test.

    python tools/atlas_parity_ref.py [--n 300]        writes tools/atlas_parity_cases.json   (then: node tools/atlas_parity.js)

For a few hundred molecules (the molecule workbook plus reactants of recorded reactions) this computes, with RDKit, exactly what the page must
compute with RDKit.js: every fragment instance of the Fragments sheet, its reacting atom, the four environment hashes, the steric score, the
Gasteiger charge, the smallest ring, the angle strain and the hybridization (tools/atlas_core.py). tools/atlas_parity.js then runs js/atlas.js in a
browser on the same SMILES and compares. A key that differs would look up the wrong counts without any error, so this is the test that matters.
"""
import argparse, json, os, random, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atlas_core as A


def instances_of(smiles):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smiles)
    if m is None: return None
    wl = A.wl_hashes(m)
    q = A.gasteiger(m)
    out = []
    for fid, atoms, x in A.instances(m):
        out.append({'fid': fid, 'atom': x, 'atoms': list(atoms), 'e': [int(wl[k][x]) for k in range(4)], 'steric': A.steric_score(m, x),
                    'charge': None if q[x] is None else round(q[x], 5), 'ring': A.smallest_ring(m, x), 'strain': round(A.angle_strain(m, x), 4),
                    'hyb': A.hybridization(m.GetAtomWithIdx(x)), 'nAtoms': len(atoms)})
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--n', type=int, default=300)
    ap.add_argument('--src', default=os.environ.get('ATLAS_SRC', ''), help='folder with the USPTO-MIT test.txt (reactants are sampled from it)')
    a = ap.parse_args()
    import learn_reactions as L
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
    smiles = []
    for r in L.read_sheet(os.path.join(L.APP, 'molecule_data.xlsx'), 'Molecules'):
        s = str(r.get('SMILES') or '').strip()
        m = Chem.MolFromSmiles(s) if s else None
        if m is not None: smiles.append(Chem.MolToSmiles(m))
    if a.src and os.path.exists(os.path.join(a.src, 'test.txt')):
        rnd = random.Random(7)
        lines = open(os.path.join(a.src, 'test.txt'), encoding='utf-8').read().split('\n')
        rnd.shuffle(lines)
        for ln in lines:
            if len(smiles) >= a.n + 187: break
            try:
                for part in ln.split(' ')[0].split('>>')[0].split('.'):
                    m = Chem.MolFromSmiles(part)
                    if m is None or m.GetNumHeavyAtoms() < 6: continue
                    for at in m.GetAtoms(): at.SetAtomMapNum(0)
                    smiles.append(Chem.MolToSmiles(m)); break
            except Exception:
                pass
    seen, cases = set(), []
    for s in smiles:
        if s in seen: continue
        seen.add(s)
        ins = instances_of(s)
        if ins is not None: cases.append({'smiles': s, 'instances': ins})
    out = os.path.join(HERE, 'atlas_parity_cases.json')
    json.dump({'about': 'reference instances computed by tools/atlas_core.py; checked against js/atlas.js by tools/atlas_parity.js', 'cases': cases}, open(out, 'w'), separators=(',', ':'))
    print(len(cases), 'molecules,', sum(len(c['instances']) for c in cases), 'instances ->', out)


if __name__ == '__main__':
    main()
