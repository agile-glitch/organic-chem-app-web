"""Grow the compound library: approved drugs (all of them, controlled or not), research compounds, and building blocks.
The new entries go to data/library_more.js (js/chem.js reads it after its own list); js/chem.js itself is not edited.

  python tools/library_add.py candidates      -> tools/data/library_candidates.json   (step 1)
  node   tools/library_keys.js                -> tools/data/library_keys.json         (step 2: the page's keys and names)
  python tools/library_add.py write           -> data/library_more.js                 (step 3)

Sources:
  drugs           every approved small-molecule drug in ChEMBL (datamol's copy, data/chembl_approved_drugs.parquet:
                  2,628 drugs with their ChEMBL preferred names, approved salts as their own entries), scheduled
                  (controlled) substances included. ChEMBL, EMBL-EBI, CC BY-SA 3.0.
  research        compounds studied as medicines that are not (or not yet) approved, many of them controlled
                  (RESEARCH below; each formula checked, each stereocentre set from its CIP label).
  building blocks the starting materials the app's learned patent reactions (data/learned/rx, US patents 1976-2016)
                  use most: every single-component neutral starting material used in at least MIN_USES reactions. They
                  are named by the app's own IUPAC namer, and kept only when the name reads back as the same molecule
                  (full InChIKey) with the app's own name parser or with OPSIN (py2opsin; needs Java).
A candidate already in the library (the same full InChIKey, so the same stereoisomer) is left out; a stereoisomer of a
library compound is a different molecule and is added.
"""
import base64, collections, glob, gzip, json, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
sys.path.insert(0, HERE)
CAND = os.path.join(HERE, 'data', 'library_candidates.json')
KEYS = os.path.join(HERE, 'data', 'library_keys.json')
COUNTS = os.path.join(HERE, 'data', 'starting_material_counts.json')
OUT = os.path.join(APP, 'data', 'library_more.js')
MIN_USES = 10
MAX_HEAVY_BB = 40

# name, SMILES (stereo written or set below from CIP), formula, {atom map number: 'R'|'S'}, alternative names
RESEARCH = [
    ('psilocin', 'CN(C)CCc1c[nH]c2cccc(O)c12', 'C12H16N2O', {}, ['4-hydroxy-N,N-dimethyltryptamine', '4-HO-DMT']),
    ('N,N-dimethyltryptamine', 'CN(C)CCc1c[nH]c2ccccc12', 'C12H16N2', {}, ['DMT']),
    ('5-methoxy-N,N-dimethyltryptamine', 'COc1ccc2[nH]cc(CCN(C)C)c2c1', 'C13H18N2O', {}, ['5-MeO-DMT']),
    ('bufotenine', 'CN(C)CCc1c[nH]c2ccc(O)cc12', 'C12H16N2O', {}, ['5-HO-DMT']),
    ('3,4-methylenedioxyamphetamine', 'CC(N)Cc1ccc2c(c1)OCO2', 'C10H13NO2', {}, ['MDA']),
    ('2C-B', 'COc1cc(CCN)c(OC)cc1Br', 'C10H14BrNO2', {}, ['4-bromo-2,5-dimethoxyphenethylamine']),
    ('cathinone', 'C[CH:1](N)C(=O)c1ccccc1', 'C9H11NO', {1: 'S'}, ['(S)-cathinone']),
    ('methcathinone', 'CN[CH:1](C)C(=O)c1ccccc1', 'C10H13NO', {1: 'S'}, ['ephedrone']),
    ('harmine', 'COc1ccc2c(c1)[nH]c1c(C)nccc12', 'C13H12N2O', {}, []),
    ('harmaline', 'COc1ccc2c(c1)[nH]c1c(C)=NCCc12', 'C13H14N2O', {}, []),
    ('muscimol', 'NCc1cc(=O)[nH]o1', 'C4H6N2O2', {}, []),
    ('ibotenic acid', 'N[CH:1](C(=O)O)c1cc(=O)[nH]o1', 'C5H6N2O4', {1: 'S'}, []),
    ('cannabinol', 'CCCCCc1cc(O)c2c(c1)OC(C)(C)c1ccc(C)cc1-2', 'C21H26O2', {}, ['CBN']),
    ('cannabigerol', 'CCCCCc1cc(O)c(C/C=C(\\C)CCC=C(C)C)c(O)c1', 'C21H32O2', {}, ['CBG']),
    ('tetrahydrocannabivarin', 'CCCc1cc(O)c2c(c1)OC(C)(C)[CH:1]1CCC(C)=C[CH:2]21', 'C19H26O2', {1: 'R', 2: 'R'}, ['THCV']),
    ('lysergic acid', 'CN1C[CH:1](C(=O)O)C=C2c3cccc4[nH]cc(c34)C[CH:2]21', 'C16H16N2O2', {1: 'R', 2: 'R'}, ['D-lysergic acid']),
    ('lysergic acid amide', 'CN1C[CH:1](C(N)=O)C=C2c3cccc4[nH]cc(c34)C[CH:2]21', 'C16H17N3O', {1: 'R', 2: 'R'}, ['ergine', 'LSA']),
    ('3,4-methylenedioxy-N-ethylamphetamine', 'CCNC(C)Cc1ccc2c(c1)OCO2', 'C12H17NO2', {}, ['MDEA', 'MDE']),
    ('mephedrone', 'CN[CH:1](C)C(=O)c1ccc(C)cc1', 'C11H15NO', {}, ['4-methylmethcathinone', '4-MMC']),
    ('norketamine', 'N[C:1]1(c2ccccc2Cl)CCCCC1=O', 'C12H14ClNO', {}, []),
    ('γ-butyrolactone', 'O=C1CCCO1', 'C4H6O2', {}, ['GBL']),
    ('1,4-butanediol', 'OCCCCO', 'C4H10O2', {}, []),
    ('tryptamine', 'NCCc1c[nH]c2ccccc12', 'C10H12N2', {}, []),
    ('phenethylamine', 'NCCc1ccccc1', 'C8H11N', {}, ['2-phenylethylamine']),
    ('arecoline', 'COC(=O)C1=CCCN(C)C1', 'C8H13NO2', {}, []),
    ('cytisine', 'O=c1cccc2n1C[CH:1]1CNC[CH:2]2C1', 'C11H14N2O', {1: 'S', 2: 'R'}, ['cytisinicline']),
]

TAGS = [("alcohols", "[CX4][OX2H]"), ("acids", "[CX3](=O)[OX2H1,OX1-]"), ("amines", "[NX3;!$(NC=O);!$(N-[S,P]=O);!$(N=*);!$(N-N);!a][#6]"),
        ("aldehydes", "[CX3H1](=O)[#6]"), ("ketones", "[#6][CX3](=O)[#6]"), ("ethers", "[OD2]([#6])[#6]"), ("esters", "[#6][CX3](=O)[OX2][#6]"),
        ("amides", "[CX3](=O)[NX3]"), ("halides", "[F,Cl,Br,I]"), ("alkenes", "[CX3]=[CX3]"), ("alkynes", "[CX2]#[CX2]"), ("nitriles", "[NX1]#[CX2]"),
        ("thiols", "[SX2H]"), ("aromatics", "a")]
METALS = re.compile(r'\[(Li|Na|K|Rb|Cs|Mg|Ca|Sr|Ba|Al|Zn|Cu|Fe|Co|Ni|Mn|Cr|Pd|Pt|Rh|Ru|Ir|Ag|Au|Hg|Sn|Pb|Ti|Zr|Ce|La|Bi|Sb|Tl|Cd|Os|Re|W|Mo|V|Ga|In|Ge|Sc|Y)')


def rd():
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
    return Chem


def library():
    """(name, SMILES) of every line of js/chem.js MOL_LIBRARY_RAW"""
    js = open(os.path.join(APP, 'js', 'chem.js'), encoding='utf-8').read()
    raw = js[js.index('const MOL_LIBRARY_RAW = `') + 25:]
    raw = raw[:raw.index('`;')]
    return [l.split('|') for l in raw.split('\n') if l.strip()]


def tags_of(m, extra):
    Chem = rd()
    out = [t for t, sm in TAGS if m.HasSubstructMatch(Chem.MolFromSmarts(sm))]
    if m.GetRingInfo().NumRings(): out.append('rings')
    if Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False): out.append('chiral')
    if any(a.GetFormalCharge() for a in m.GetAtoms()): out.append('ions')
    return ','.join(out + extra)


def kekule(m):
    Chem = rd()
    m = Chem.Mol(m)
    try: Chem.Kekulize(m, clearAromaticFlags=True)
    except Exception: return Chem.MolToSmiles(m)
    return Chem.MolToSmiles(m, kekuleSmiles=True)


def _count_chunk(f):
    import rxn_code
    t = open(f, encoding='utf-8').read()
    j = json.loads(t[t.index('=', t.index(']')) + 1:].rstrip().rstrip(';'))
    c = collections.Counter()
    for row in json.loads(gzip.decompress(base64.b64decode(j['gz']))):
        p = rxn_code.parse(row[1])
        if p['ok']:
            for s in set(p['reactants']): c[s] += 1
    return c


def starting_material_counts():
    """{SMILES: number of learned reactions using it as a starting material}, from data/learned/rx (cached)"""
    if os.path.exists(COUNTS): return dict(json.load(open(COUNTS)))
    from multiprocessing import Pool
    C = collections.Counter()
    with Pool() as pool:
        for c in pool.imap_unordered(_count_chunk, sorted(glob.glob(os.path.join(APP, 'data', 'learned', 'rx', '*.js')))): C.update(c)
    out = [(s, n) for s, n in C.most_common() if n >= 3]
    json.dump(out, open(COUNTS, 'w'))
    return dict(out)


def set_cip(smi, want):
    """the stereoisomer of smi whose mapped atoms have the CIP labels asked for"""
    Chem = rd()
    from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
    m = Chem.MolFromSmiles(smi)
    if not want:
        for a in m.GetAtoms(): a.SetAtomMapNum(0)
        return m
    from rdkit.Chem import rdCIPLabeler
    for iso in EnumerateStereoisomers(m, options=StereoEnumerationOptions(onlyUnassigned=True, unique=True)):
        rdCIPLabeler.AssignCIPLabels(iso)            # the new (exact) CIP labeller, not the legacy approximation
        lab = {a.GetAtomMapNum(): a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None for a in iso.GetAtoms() if a.GetAtomMapNum()}
        if all(lab.get(k) == v for k, v in want.items()):
            for a in iso.GetAtoms(): a.SetAtomMapNum(0)
            # only the centres asked for keep their stereo
            return Chem.MolFromSmiles(Chem.MolToSmiles(iso))
    raise ValueError(f'no stereoisomer of {smi} has {want}')


def candidates():
    Chem = rd()
    from rdkit.Chem import inchi, rdMolDescriptors
    import pandas as pd
    have, names = set(), set()
    for f in library():
        names.add(f[0].lower())
        for a in (f[5].split(';') if len(f) > 5 and f[5] else []): names.add(a.lower())
        m = Chem.MolFromSmiles(f[1])
        if m is not None: have.add(inchi.MolToInchiKey(m))
    out, why = [], collections.Counter()

    def add(name, m, src, tags, alt=(), uses=None):
        if m is None: why[f'{src}: unreadable SMILES'] += 1; return
        key = inchi.MolToInchiKey(m)
        if not key: why[f'{src}: no InChIKey'] += 1; return
        if key in have: why[f'{src}: already in the library'] += 1; return
        if name and name.lower() in names: why[f'{src}: name already used by another compound'] += 1; return
        have.add(key)
        if name: names.add(name.lower())
        out.append({'name': name, 'smiles': kekule(m), 'inchikey': key, 'formula': rdMolDescriptors.CalcMolFormula(m),
                    'tags': tags_of(m, tags), 'alt': [a for a in alt if a.lower() not in names], 'src': src, **({'uses': uses} if uses else {})})
        for a in alt: names.add(a.lower())

    d = pd.read_parquet(os.path.join(os.path.dirname(__import__('datamol').__file__), 'data', 'chembl_approved_drugs.parquet'))
    for _, r in d.sort_values('first_approval').iterrows():
        add(r.pref_name.lower(), Chem.MolFromSmiles(r.smiles), 'drug', ['drugs'], [r.molecule_chembl_id])
    for name, smi, formula, want, alt in RESEARCH:
        m = set_cip(smi, want)
        f = rdMolDescriptors.CalcMolFormula(m)
        if f != formula: raise SystemExit(f'{name}: formula {f}, expected {formula}')
        add(name, m, 'research', ['drugs'], alt)
    counts = starting_material_counts()
    for smi, n in sorted(counts.items(), key=lambda x: -x[1]):
        if n < MIN_USES: break
        if '.' in smi or METALS.search(smi): why['building block: salt, mixture or metal'] += 1; continue
        m = Chem.MolFromSmiles(smi)
        if m is None: why['building block: unreadable SMILES'] += 1; continue
        if any(a.GetFormalCharge() for a in m.GetAtoms()) and Chem.GetFormalCharge(m): why['building block: an ion'] += 1; continue
        if any(a.GetIsotope() or a.GetNumRadicalElectrons() for a in m.GetAtoms()):
            why['building block: isotope or radical'] += 1; continue
        if m.GetNumHeavyAtoms() > MAX_HEAVY_BB or m.GetNumHeavyAtoms() < 2: why['building block: size'] += 1; continue
        add(None, m, 'building block', ['building blocks'], uses=n)
    json.dump(out, open(CAND, 'w'), ensure_ascii=False)
    c = collections.Counter(x['src'] for x in out)
    print(f'{len(out)} candidates -> {os.path.relpath(CAND, APP)}:', dict(c))
    for k, v in sorted(why.items()): print(f'  left out, {k}: {v}')


def opsin_check(cand, keys):
    """{InChIKey: True} for each unnamed candidate whose app-made IUPAC name OPSIN reads as the same molecule (full InChIKey)"""
    import warnings
    from py2opsin import py2opsin
    Chem = rd()
    from rdkit.Chem import inchi
    todo = [(c['inchikey'], keys[c['inchikey']]['iupac']) for c in cand
            if not c['name'] and c['inchikey'] in keys and keys[c['inchikey']].get('iupac') and keys[c['inchikey']].get('checked') != 'yes']
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        for i in range(0, len(todo), 2000):
            part = todo[i:i + 2000]
            for (ik, _), smi in zip(part, py2opsin([n for _, n in part])):
                m = Chem.MolFromSmiles(smi) if smi else None
                out[ik] = bool(m is not None and inchi.MolToInchiKey(m) == ik)
    return out


def write():
    cand = json.load(open(CAND))
    keys = {r['inchikey']: r for r in json.load(open(KEYS))['rows']}
    opsin = opsin_check(cand, keys)
    lines, why, names = [], collections.Counter(), {f[0].lower() for f in library()}
    for c in cand:
        k = keys.get(c['inchikey'])
        if not k or not k.get('key'): why[f"{c['src']}: the page cannot read it"] += 1; continue
        name = c['name']
        if not name:
            if not k.get('iupac'): why['building block: the app\'s namer gives no name'] += 1; continue
            if k.get('checked') == 'yes': why['(kept) building block name read back by the app\'s parser'] += 1
            elif opsin.get(c['inchikey']): why['(kept) building block name read back by OPSIN'] += 1
            else: why['building block: its name does not read back as the same molecule'] += 1; continue
            name = k['iupac']
        if name.lower() in names or '|' in name or '`' in name: why[f"{c['src']}: name taken"] += 1; continue
        names.add(name.lower())
        lines.append('|'.join([name, c['smiles'], k['key'], c['formula'], c['tags'], ';'.join(c['alt']), '']))
    head = ('/* generated by tools/library_add.py: do not edit by hand. More compounds for the library (js/chem.js reads them after its\n'
            '   own list, same line format: name|SMILES|key|formula|tags|other names|CAS). Approved drugs from ChEMBL (EMBL-EBI,\n'
            '   CC BY-SA 3.0; the ChEMBL ID is an other name), research compounds, and the building blocks the learned patent\n'
            '   reactions use most, named by the app\'s IUPAC namer. */\n')
    open(OUT, 'w', encoding='utf-8').write(head + 'window.MOL_LIBRARY_MORE = `' + '\n'.join(lines) + '`;\n')
    c = collections.Counter('drugs' if 'building blocks' not in l.split('|')[4] else 'building blocks' for l in lines)
    print(f'{len(lines)} compounds -> {os.path.relpath(OUT, APP)} ({os.path.getsize(OUT) / 1e6:.2f} MB):', dict(c))
    for k, v in sorted(why.items()): print(f'  {k}: {v}' if k.startswith('(kept)') else f'  left out, {k}: {v}')


if __name__ == '__main__':
    {'candidates': candidates, 'write': write}[sys.argv[1] if len(sys.argv) > 1 else 'candidates']()
