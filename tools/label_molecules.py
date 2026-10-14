"""The molecule set for the label validation: every molecule the app's own data holds, plus a hand-made list of unusual
ones, each tagged with the classes it belongs to, so every check can be read class by class.

  sources: data/pka.js (IUPAC measured pKa, small acids and bases), data/reactions.js (patent reactions: drug-like
  molecules, reagents, catalysts), data/iupac_names.js (the namer's systematic names), tools/data/reagents.json
  (reagents, salts, metals), tools/label_special.json (ions, radicals, isotopes, metals, hypervalent, ylides, stereo of
  every kind, macrocycles, polycyclic and heteroaromatics, natural products).

  python tools/label_molecules.py   -> tools/data/label_molecules.json  [{smiles, name, sources, classes}]
"""
import json, os
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog('rdApp.*')
HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
ORGANIC = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 34, 35, 53}


def js_obj(path):
    s = open(path, encoding='utf-8').read()
    return json.loads(s[s.index('=') + 1:].strip().rstrip(';'))


def classes(m):
    c = []
    heavy = m.GetNumHeavyAtoms()
    frags = Chem.GetMolFrags(m)
    if len(frags) > 1: c.append('salt or mixture')
    charges = [a.GetFormalCharge() for a in m.GetAtoms()]
    if any(charges):
        c.append('zwitterion or dipole' if any(x > 0 for x in charges) and any(x < 0 for x in charges) and sum(charges) == 0 and len(frags) == 1 else 'ion')
    if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()): c.append('radical or carbene')
    if any(a.GetIsotope() for a in m.GetAtoms()): c.append('isotope label')
    els = {a.GetAtomicNum() for a in m.GetAtoms()}
    if els - ORGANIC: c.append('metal or other element')
    if any((a.GetSymbol() in ('S', 'P', 'Cl', 'Br', 'I', 'Se', 'Xe') and a.GetTotalValence() > {'S': 2, 'P': 3, 'Se': 2}.get(a.GetSymbol(), 1)) for a in m.GetAtoms()):
        c.append('hypervalent S / P / halogen')
    centres = Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False)
    if any(lbl != '?' for _, lbl in centres): c.append('stereocentres')
    if any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in m.GetBonds()): c.append('E/Z double bond')
    ri = m.GetRingInfo()
    if any(len(r) >= 8 for r in ri.AtomRings()): c.append('macrocycle (ring of 8+)')
    if any(len(r) <= 4 for r in ri.AtomRings()): c.append('small ring (3-4)')
    narom = rdMolDescriptors.CalcNumAromaticRings(m)
    if narom >= 3: c.append('3+ aromatic rings')
    if any(a.GetIsAromatic() and a.GetAtomicNum() != 6 for a in m.GetAtoms()): c.append('heteroaromatic')
    if any(a.GetAtomicNum() == 6 and sum(1 for b in a.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE) >= 2 for a in m.GetAtoms()): c.append('cumulene / allene / ketene')
    if heavy >= 40: c.append('large (40+ heavy atoms)')
    elif heavy <= 6: c.append('tiny (6 or fewer heavy atoms)')
    if not c: c.append('ordinary organic')
    return c


def main():
    got = {}
    def add(smi, src, name=''):
        if not smi or len(smi) > 400: return
        m = Chem.MolFromSmiles(smi)
        if m is None or m.GetNumAtoms() == 0: return
        e = got.setdefault(smi, {'smiles': smi, 'name': name, 'sources': [], 'classes': classes(m)})
        if src not in e['sources']: e['sources'].append(src)
        if name and not e['name']: e['name'] = name
    for smi in js_obj(os.path.join(APP, 'data', 'pka.js'))['molecules']: add(smi, 'IUPAC pKa set')
    for x in js_obj(os.path.join(APP, 'data', 'reactions.js'))['molecules']: add(x['s'], 'patent reactions')
    for name, smi in js_obj(os.path.join(APP, 'data', 'iupac_names.js')): add(smi, 'named compounds', name)
    for x in json.load(open(os.path.join(HERE, 'data', 'reagents.json'))): add(x['smiles'], 'reagents', x.get('name', ''))
    for name, smi in json.load(open(os.path.join(HERE, 'label_special.json')))['molecules']: add(smi, 'unusual (hand-made)', name)
    out = list(got.values())
    json.dump(out, open(os.path.join(HERE, 'data', 'label_molecules.json'), 'w'), ensure_ascii=False)
    from collections import Counter
    cc = Counter(c for e in out for c in e['classes'])
    sc = Counter(s for e in out for s in e['sources'])
    print(f'{len(out):,} molecules ->tools/data/label_molecules.json')
    print('sources:', dict(sc))
    for k, v in cc.most_common(): print(f'  {k:30s} {v:6,}')


if __name__ == '__main__':
    main()
