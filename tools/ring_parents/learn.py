#!/usr/bin/env python3
"""Learn ring parent hydrides for the IUPAC namer (js/chem.js) from OPSIN, the open-source IUPAC name reader.

The namer names a ring system by matching it to a parent hydride and taking that parent's numbering. Writing those
numberings by hand does not scale (and is easy to get wrong), so they are learned: candidate names are written by the
rules of IUPAC nomenclature, OPSIN reads each one (its CML output gives every atom's locant, ring-fusion atoms such as 4a
included), and every name OPSIN accepts becomes a template: its atoms, bonds, which atoms are sp3 in the parent (the
indicated hydrogen), and the locant of every atom. Kinds:

  mono    heteromonocycles: the retained names (pyridine, imidazole, morpholine ...) and Hantzsch-Widman names
          (1,3-oxazole, 1,2,4-oxadiazole, 1,4-oxazepane ...), each with its saturated name (piperidine, 1,3-dioxolane)
  fused   ortho-fused systems: retained names (quinoline, purine, carbazole ...), benzo names (1,3-benzoxazole,
          2H-1-benzopyran) and fusion names (thieno[3,2-b]pyridine, imidazo[1,2-a]pyridine ...), the base component
          chosen by the seniority rules of IUPAC P-25.3.1.3, the fusion descriptor the earliest one that gives the system
  vb      von Baeyer skeletons (bicyclo[2.2.1]heptane ...; adamantane): heteroatoms are named by replacement prefixes
          (7-oxabicyclo[2.2.1]heptane), so a skeleton template is all carbon
  spiro   monospiro skeletons (spiro[4.5]decane), heteroatoms by replacement prefixes (2-oxa-6-azaspiro[3.3]heptane)

    python tools/ring_parents/learn.py            writes the generated table into js/chem.js (between its RING_PARENTS
                                                  markers), and the full table to tools/ring_parents/out/ring_parents.json
Needs Java and OPSIN (as tools/namer_bench/bench.py) and RDKit (to tell structures apart).
"""
import itertools, json, os, re, subprocess, sys
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(REPO, 'tools', 'namer_bench'))
from bench import java, opsin_jar                                   # noqa: E402

CML = '{http://www.xml-cml.org/schema}'
SENIOR = ['O', 'S', 'Se', 'Te', 'N', 'P']                         # replacement / Hantzsch-Widman order
HW_PREFIX = {'O': 'oxa', 'S': 'thia', 'Se': 'selena', 'Te': 'tellura', 'N': 'aza', 'P': 'phospha'}
MULT = {1: '', 2: 'di', 3: 'tri', 4: 'tetra', 5: 'penta', 6: 'hexa'}
ALKANE = {3: 'propane', 4: 'butane', 5: 'pentane', 6: 'hexane', 7: 'heptane', 8: 'octane', 9: 'nonane', 10: 'decane',
          11: 'undecane', 12: 'dodecane', 13: 'tridecane', 14: 'tetradecane', 15: 'pentadecane', 16: 'hexadecane'}


# ---------------------------------------------------------------------------------------------------- OPSIN
def opsin_parse(names, chunk=4000):
    """for each name: None, or {'atoms': [(element, locant, hydrogens)], 'bonds': [(i, j, order)]} (heavy atoms only)"""
    out = []
    for k in range(0, len(names), chunk):
        part = names[k:k + chunk]
        p = subprocess.run([java(), '-jar', opsin_jar(), '-ocml'], input='\n'.join(part) + '\n', capture_output=True,
                           text=True, encoding='utf-8')
        root = ET.fromstring(p.stdout)
        mols = root.findall(CML + 'molecule')
        if len(mols) != len(part): raise RuntimeError('OPSIN gave %d molecules for %d names' % (len(mols), len(part)))
        for m in mols:
            arr = m.find(CML + 'atomArray')
            if arr is None: out.append(None); continue
            idx, atoms = {}, []
            for a in arr:
                if a.get('elementType') == 'H': continue
                locs = [l.get('value') for l in a.findall(CML + 'label') if l.get('dictRef') == 'cmlDict:locant']
                num = [l for l in locs if re.match(r'^\d+[a-z]?\'*$', l)]
                idx[a.get('id')] = len(atoms)
                atoms.append([a.get('elementType'), num[0] if num else None, int(a.get('hydrogenCount') or -1), a.get('formalCharge')])
            bonds = []
            for b in (m.find(CML + 'bondArray') if m.find(CML + 'bondArray') is not None else []):
                x, y = b.get('atomRefs2').split()
                if x in idx and y in idx: bonds.append((idx[x], idx[y], {'S': 1, 'D': 2, 'T': 3, 'A': 4}[b.get('order')]))
            hcount = [0] * len(atoms)
            for b in (m.find(CML + 'bondArray') if m.find(CML + 'bondArray') is not None else []):
                x, y = b.get('atomRefs2').split()
                if x in idx and y not in idx: hcount[idx[x]] += 1
                if y in idx and x not in idx: hcount[idx[y]] += 1
            for i, a in enumerate(atoms): a[2] = hcount[i]
            out.append({'atoms': atoms, 'bonds': bonds})
    return out


def ring_atoms(mol):
    """the atoms in rings (a parent hydride is all ring), by repeatedly pruning atoms of degree 1"""
    n = len(mol['atoms'])
    adj = [set() for _ in range(n)]
    for i, j, _ in mol['bonds']: adj[i].add(j); adj[j].add(i)
    alive = set(range(n))
    changed = True
    while changed:
        changed = False
        for i in list(alive):
            if len(adj[i] & alive) <= 1: alive.discard(i); changed = True
    return alive


def template(name, mol, kind, extra=None):
    """the template of a parsed parent (None when it is not a single ring system of ring atoms only)"""
    atoms, bonds = mol['atoms'], mol['bonds']
    if any(a[3] for a in atoms): return None                        # charged
    if len(ring_atoms(mol)) != len(atoms): return None              # something not in a ring
    if any(a[1] is None for a in atoms): return None
    unsat = [False] * len(atoms)
    for i, j, o in bonds:
        if o >= 2: unsat[i] = unsat[j] = True
    # an atom with no double bond in the parent is either an indicated-hydrogen position (a C or N carrying H: the
    # 1H of 1H-indole, the CH2 of 2H-1,3-benzodioxole) or always sp3 (O, S, an N with no H such as indolizine's)
    ih = [not unsat[i] and a[2] > 0 and a[0] in ('C', 'N') for i, a in enumerate(atoms)]
    fixed = [not unsat[i] and not ih[i] for i in range(len(atoms))]
    t = {'nm': re.sub(r'^(\d+[a-z]?H,?)+-', '', name), 'kind': kind, 'el': [a[0] for a in atoms], 'b': [[i, j] for i, j, _ in bonds],
         'sat': [1 if f else 0 for f in fixed], 'ih': sum(ih), 'loc': [a[1] for a in atoms]}
    if extra: t.update(extra)
    return t


def numbering_of(t):
    """the template's numbering, independent of atom order: (locant, element, neighbour locants) of every atom"""
    nb = {i: [] for i in range(len(t['el']))}
    for i, j in t['b']: nb[i].append(t['loc'][j]); nb[j].append(t['loc'][i])
    return sorted((t['loc'][i], t['el'][i], tuple(sorted(nb[i]))) for i in range(len(t['el'])))


def key_of(t, elements=True):
    """a canonical key of a template's graph (RDKit), elements included or not"""
    from rdkit import Chem
    m = Chem.RWMol()
    for e in t['el']:
        a = Chem.Atom(e if elements else 'C'); a.SetNoImplicit(True); m.AddAtom(a)
    for i, j in t['b']: m.AddBond(i, j, Chem.BondType.SINGLE)
    return Chem.MolToSmiles(m.GetMol(), canonical=True)


# ---------------------------------------------------------------------------------------------------- monocycles
HW_STEM = {  # ring size: (mancude, saturated) for a ring whose last-cited heteroatom is N / is not N (6: by class)
    3: {'N': ('irine', 'iridine'), '*': ('irene', 'irane')}, 4: {'N': ('ete', 'etidine'), '*': ('ete', 'etane')},
    5: {'N': ('ole', 'olidine'), '*': ('ole', 'olane')}, 6: {'N': ('ine', 'inane'), 'P': ('inine', 'inane'), '*': ('ine', 'ane')},
    7: {'*': ('epine', 'epane')}, 8: {'*': ('ocine', 'ocane')}, 9: {'*': ('onine', 'onane')}, 10: {'*': ('ecine', 'ecane')}}
RETAINED_MONO = [  # (mancude name, saturated name or None): IUPAC 2013 retained names
    ('furan', None), ('thiophene', None), ('1H-pyrrole', 'pyrrolidine'), ('1H-imidazole', 'imidazolidine'),
    ('1H-pyrazole', 'pyrazolidine'), ('pyridine', 'piperidine'), ('pyrazine', 'piperazine'), ('pyrimidine', None),
    ('pyridazine', None), ('2H-pyran', None), ('2H-1,4-oxazine', 'morpholine'), ('2H-1,4-thiazine', 'thiomorpholine')]


def elide(prefix, stem):
    return prefix[:-1] + stem if prefix.endswith('a') and stem[0] in 'aeiou' else prefix + stem


def hw_name(n, het, saturated):
    """the Hantzsch-Widman name of an n-membered ring with heteroatoms het {locant: element}"""
    els = sorted(set(het.values()), key=SENIOR.index)
    last = els[-1]
    row = HW_STEM[n]
    if n == 6: col = 'N' if last == 'N' else 'P' if last == 'P' else '*'          # 6A / 6B / 6C: the last-cited heteroatom
    else: col = 'N' if 'N' in els and 'N' in row else '*'                          # 3-5: whether the ring holds N
    stem = row[col][1 if saturated else 0]
    locs, pre = [], ''
    for e in els:
        ls = sorted(l for l, x in het.items() if x == e)
        locs += ls
        m = MULT[len(ls)]
        p = HW_PREFIX[e]
        pre = elide(pre, m + p) if pre else m + p
    word = elide(pre, stem)
    word = word.replace('tetraaz', 'tetraz').replace('triaaz', 'triaz')
    return (','.join(map(str, locs)) + '-' if len(het) > 1 else '') + word


def monocycle_names():
    names = []
    for n in range(3, 9):
        for k in range(1, min(4, n - 1) + 1):
            if n >= 7 and k > 3: continue
            for pos in itertools.combinations(range(1, n + 1), k):
                if pos[0] != 1: continue
                for els in itertools.product(['O', 'S', 'N'], repeat=k):
                    het = dict(zip(pos, els))
                    if els[0] != sorted(els, key=SENIOR.index)[0]: continue      # the most senior heteroatom at 1
                    names.append((hw_name(n, het, False), hw_name(n, het, True)))
    return names


# ---------------------------------------------------------------------------------------------------- fused
RETAINED_FUSED = """pentalene 1H-indene naphthalene azulene heptalene biphenylene as-indacene s-indacene acenaphthylene 9H-fluorene
1H-phenalene phenanthrene anthracene fluoranthene acephenanthrylene aceanthrylene triphenylene pyrene chrysene tetracene
pleiadene picene perylene pentaphene pentacene 1H-indole 2H-isoindole indolizine 1H-indazole 7H-purine 9H-purine 1H-purine
3H-purine 1H-pyrrolizine quinoline isoquinoline 4H-quinolizine phthalazine 1,8-naphthyridine 1,5-naphthyridine
1,6-naphthyridine 1,7-naphthyridine 2,6-naphthyridine 2,7-naphthyridine quinoxaline quinazoline cinnoline pteridine
9H-carbazole acridine phenazine phenanthridine 1,10-phenanthroline 1,7-phenanthroline 4,7-phenanthroline perimidine
10H-phenothiazine 10H-phenoxazine phenoxathiine thianthrene 9H-xanthene 9H-thioxanthene 1H-benzimidazole
dibenzo[b,d]furan dibenzo[b,d]thiophene 9H-pyrido[3,4-b]indole 9H-pyrido[2,3-b]indole 5H-pyrido[4,3-b]indole
5H-dibenzo[b,f]azepine 5H-dibenzo[a,d][7]annulene 6H-benzo[c]chromene
dibenzo[b,e][1,4]dioxine 5H-benzo[7]annulene benzo[a]anthracene
benzo[c]phenanthrene benzo[a]pyrene""".split()
PREFIXES = [  # (prefix, its own locant count (0: benzo/naphtho by letters only), kind) — monocyclic first
    ('benzo', 0), ('naphtho', 8), ('cyclopenta', 0), ('cyclohepta', 0), ('cyclobuta', 0), ('pyrido', 6), ('pyrimido', 6),
    ('pyrazino', 6), ('pyridazino', 6), ('pyrrolo', 5), ('imidazo', 5), ('pyrazolo', 5), ('furo', 5), ('thieno', 5),
    ('[1,3]oxazolo', 5), ('[1,3]thiazolo', 5), ('[1,2]oxazolo', 5), ('[1,2]thiazolo', 5), ('[1,2,4]triazolo', 5),
    ('[1,2,3]triazolo', 5), ('tetrazolo', 5), ('[1,2,4]triazino', 6), ('[1,3]dioxolo', 5), ('[1,4]dioxino', 6),
    ('[1,4]oxazino', 6), ('[1,4]thiazino', 6), ('pyrano', 6), ('thiopyrano', 6), ('azepino', 7), ('oxepino', 7),
    ('[1,2,5]oxadiazolo', 5), ('[1,2,5]thiadiazolo', 5), ('[1,3,4]thiadiazolo', 5), ('[1,3,4]oxadiazolo', 5),
    ('indolo', 9), ('quinolino', 10), ('isoquinolino', 10)]
BASES = ['pyridine', 'pyrimidine', 'pyrazine', 'pyridazine', '1,2,4-triazine', '1,3,5-triazine', 'pyrrole', 'imidazole',
         'pyrazole', 'furan', 'thiophene', '1,3-oxazole', '1,3-thiazole', '1,2-oxazole', '1,2-thiazole', '1,2,4-triazole',
         '1,2,3-triazole', '1,2,5-oxadiazole', '1,2,5-thiadiazole', '1,3,4-thiadiazole', '1,3,4-oxadiazole', 'tetrazole',
         'pyran', 'thiopyran', 'azepine', 'oxepine', '1,4-oxazine', '1,4-thiazine', '1,4-dioxine', '1,3-dioxole',
         '1,4-diazepine', '1,4-oxazepine', '1,3-diazepine', 'indole', 'quinoline', 'isoquinoline', 'quinazoline',
         'quinoxaline', 'benzimidazole', 'indazole', 'purine', 'carbazole', 'phthalazine', 'pteridine', 'naphthyridine']


def component_info(mol):
    """seniority facts of a parsed component: (has N, heteroatoms in seniority order, rings, ring sizes desc, ...)"""
    els = [a[0] for a in mol['atoms']]
    het = sorted([e for e in els if e != 'C'], key=SENIOR.index)
    nrings = len(mol['bonds']) - len(els) + 1
    return {'N': 'N' in het, 'het': het, 'rings': nrings, 'size': len(els)}


def senior(a, b):
    """True when component a is senior to b as a fusion base component (IUPAC P-25.3.1.3, its first criteria)"""
    if a['N'] != b['N']: return a['N']
    if not a['N']:
        ha = [SENIOR.index(e) for e in a['het']] or [99]; hb = [SENIOR.index(e) for e in b['het']] or [99]
        if min(ha) != min(hb): return min(ha) < min(hb)
    if a['rings'] != b['rings']: return a['rings'] > b['rings']
    if a['size'] != b['size']: return a['size'] > b['size']
    if len(a['het']) != len(b['het']): return len(a['het']) > len(b['het'])
    if len(set(a['het'])) != len(set(b['het'])): return len(set(a['het'])) > len(set(b['het']))
    for x, y in zip(sorted(SENIOR.index(e) for e in a['het']), sorted(SENIOR.index(e) for e in b['het'])):
        if x != y: return x < y
    return False


PREFIX_PARENT = {'benzo': 'benzene', 'naphtho': 'naphthalene', 'cyclopenta': 'cyclopentadiene', 'cyclohepta': 'cycloheptatriene',
                 'cyclobuta': 'cyclobutadiene', 'pyrido': 'pyridine', 'pyrimido': 'pyrimidine', 'pyrazino': 'pyrazine',
                 'pyridazino': 'pyridazine', 'pyrrolo': 'pyrrole', 'imidazo': 'imidazole', 'pyrazolo': 'pyrazole',
                 'furo': 'furan', 'thieno': 'thiophene', '[1,3]oxazolo': '1,3-oxazole', '[1,3]thiazolo': '1,3-thiazole',
                 '[1,2]oxazolo': '1,2-oxazole', '[1,2]thiazolo': '1,2-thiazole', '[1,2,4]triazolo': '1,2,4-triazole',
                 '[1,2,3]triazolo': '1,2,3-triazole', 'tetrazolo': 'tetrazole', '[1,2,4]triazino': '1,2,4-triazine',
                 '[1,3]dioxolo': '1,3-dioxole', '[1,4]dioxino': '1,4-dioxine', '[1,4]oxazino': '1,4-oxazine',
                 '[1,4]thiazino': '1,4-thiazine', 'pyrano': 'pyran', 'thiopyrano': 'thiopyran', 'azepino': 'azepine',
                 'oxepino': 'oxepine', '[1,2,5]oxadiazolo': '1,2,5-oxadiazole', '[1,2,5]thiadiazolo': '1,2,5-thiadiazole',
                 '[1,3,4]thiadiazolo': '1,3,4-thiadiazole', '[1,3,4]oxadiazolo': '1,3,4-oxadiazole', 'indolo': 'indole',
                 'quinolino': 'quinoline', 'isoquinolino': 'isoquinoline'}


def fusion_names(info):
    """candidate fusion names prefix[locants-letter]base with the base senior to the prefix; in order of preference
    (earliest letter, then lowest prefix locants), so the first name OPSIN accepts for a structure is kept"""
    out = []
    letters = 'abcdefghijklmnop'
    for pre, nloc in PREFIXES:
        pi = info.get(PREFIX_PARENT[pre])
        for base in BASES:
            bi = info.get(base)
            if not pi or not bi or not senior(bi, pi): continue
            nb_ = bi['size'] + bi['rings'] - 1                       # peripheral bonds (enough letters)
            bname = re.sub(r"^([0-9,]+)-", lambda m: "[" + m.group(1) + "]", base)   # 1,3-oxazole is cited as [1,3]oxazole
            for L in letters[:min(nb_, 14)]:
                if nloc == 0:
                    out.append(((L, ()), f'{pre}[{L}]{bname}'))
                    continue
                for x in range(1, nloc + 1):
                    for y in (x - 1, x + 1):
                        if 1 <= y <= nloc: out.append(((L, (min(x, y), x, y)), f'{pre}[{x},{y}-{L}]{bname}'))
    out.sort(key=lambda t: t[0])
    return [n for _, n in out]


def benzo_name(t):
    """the benzo name of benzo + a heteromonocycle (1,3-benzoxazole, 2H-1-benzopyran), from the fused template"""
    m = re.match(r'^(\d+H-)?benzo\[[a-z]\](?:\[([\d,]+)\])?(.+)$', t['nm'])
    if not m or len(t['b']) - len(t['el']) + 1 != 2: return None      # benzo + one ring only
    ring = m.group(3)
    het = [(SENIOR.index(e), int(l)) for e, l in zip(t['el'], t['loc']) if e != 'C' and l.isdigit()]
    locs = ','.join(str(l) for _, l in sorted(het))
    word = ('benz' + ring) if ring[0] in 'aeiou' else ('benzo' + ring)
    return (m.group(1) or '') + (locs + '-' if locs else '') + word


# ---------------------------------------------------------------------------------------------------- von Baeyer, spiro
def vb_names():
    out = []
    for n in range(4, 15):
        for a in range(1, n):
            for b in range(1, a + 1):
                c = n - 2 - a - b
                if c < 0 or c > b: continue
                if c == 0 and min(a, b) + 2 >= 5: continue          # two rings of 5 or more: a fusion name instead
                out.append(f'bicyclo[{a}.{b}.{c}]{ALKANE[n]}')
    out += ['adamantane', 'tricyclo[3.3.1.1~3,7~]decane', 'tricyclo[2.2.1.0~2,6~]heptane', 'tricyclo[4.4.0.0~3,8~]decane',
            'tricyclo[5.2.1.0~2,6~]decane', 'tricyclo[3.2.1.0~2,4~]octane', 'tricyclo[3.3.0.0~2,8~]octane',
            'tricyclo[4.3.1.1~3,8~]undecane', 'pentacyclo[4.2.0.0~2,5~.0~3,8~.0~4,7~]octane', 'tricyclo[1.1.1.0~1,3~]pentane',
            'tricyclo[3.1.1.0~3,6~]heptane', 'tricyclo[4.1.0.0~2,7~]heptane', 'bicyclo[2.2.2]octane']
    return out


def spiro_names():
    return [f'spiro[{a}.{b}]{ALKANE[a + b + 1]}' for a in range(2, 8) for b in range(a, 8)]


# ---------------------------------------------------------------------------------------------------- learning
def learn():
    tpl = {}

    def keep(t, rank):
        if not t: return
        k = (t['kind'] in ('vb', 'spiro'), key_of(t, elements=t['kind'] not in ('vb', 'spiro')))
        old = tpl.get(k)
        if old is None or rank < old[0]: tpl[k] = (rank, t)

    # monocycles: retained first (rank 0), then Hantzsch-Widman (rank 1)
    hw = monocycle_names()
    mono_m = opsin_parse([m for m, _ in hw]); mono_s = opsin_parse([s for _, s in hw])
    pairs = [(m, s, pm, ps) for (m, s), pm, ps in zip(hw, mono_m, mono_s)]
    ret_m = opsin_parse([m for m, _ in RETAINED_MONO]); ret_s = opsin_parse([s or 'x' for _, s in RETAINED_MONO])
    for (m, s), pm, ps in zip(RETAINED_MONO, ret_m, ret_s):
        t = template(m, pm, 'mono') if pm else None
        if t:
            sat_name = s if (s and ps and key_of(template(s, ps, 'mono') or t) == key_of(t)) else None
            if sat_name is None:                                     # its Hantzsch-Widman saturated name
                for mm, ss, pmm, pss in pairs:
                    if pss and pmm and key_of(template(ss, pss, 'mono') or t) == key_of(t): sat_name = ss; break
            if sat_name: t['satName'] = sat_name
            keep(t, 0)
    seen_hw = 0
    for m, s, pm, ps in pairs:
        if not ps: continue
        ts = template(s, ps, 'mono')
        if not ts: continue
        t = template(m, pm, 'mono') if pm else None
        if t and key_of(t) == key_of(ts): t['satName'] = s
        else:                                                         # no mancude form: the saturated ring only
            t = dict(ts); t['satOnly'] = True
        keep(t, 1); seen_hw += 1
    print(f'monocycles: {len(hw)} Hantzsch-Widman name pairs tried, {seen_hw} read', flush=True)

    # fused: retained (0), benzo names (1), fusion names (2)
    parsed = opsin_parse(RETAINED_FUSED)
    for nm, pm in zip(RETAINED_FUSED, parsed):
        if pm: keep(template(nm, pm, 'fused'), 0)
        else: print('  OPSIN does not read', nm)
    comps = sorted(set(BASES) | set(PREFIX_PARENT.values()))
    info = {c: component_info(p) for c, p in zip(comps, opsin_parse(comps)) if p}
    cands = fusion_names(info)
    print(f'fused: {len(cands)} fusion names tried', flush=True)
    for nm, pm in zip(cands, opsin_parse(cands)):
        if not pm: continue
        t = template(nm, pm, 'fused')
        if not t: continue
        bz = benzo_name(t) if nm.startswith('benzo[') else None
        if bz:
            t['fusionName'] = t['nm']; t['nm'] = re.sub(r'^(\d+[a-z]?H,?)+-', '', bz)
            keep(t, 1)
        else:
            keep(t, 2)
    # the benzo names, read back by OPSIN: a name it reads as another system is dropped
    bz = [(k, v) for k, v in tpl.items() if v[1].get('fusionName')]
    for (k, (rank, t)), pm in zip(bz, opsin_parse([t['nm'] for _, (_, t) in bz])):
        tt = template(t['nm'], pm, 'fused') if pm else None
        if not tt or numbering_of(tt) != numbering_of(t):
            t['nm'] = t.pop('fusionName')
        else:
            t.pop('fusionName')

    # von Baeyer and spiro skeletons (all carbon)
    for kind, names in (('vb', vb_names()), ('spiro', spiro_names())):
        ok = 0
        for nm, pm in zip(names, opsin_parse(names)):
            t = template(nm, pm, kind)
            if t: keep(t, 0); ok += 1
        print(f'{kind}: {ok} of {len(names)} skeletons read', flush=True)
    out = [t for _, t in sorted(tpl.values(), key=lambda v: (v[1]['kind'], len(v[1]['el']), v[1]['nm']))]
    return out


def to_js(tpls):
    """the compact table written into js/chem.js: one string per template
       kind|name (no indicated hydrogen)|elements(,)|bonds i-j(,)|always-sp3 atoms(,)|indicated hydrogens|locants(,)|
       saturated name (monocycles; '!' + name: the saturated ring only)"""
    rows = []
    for t in tpls:
        rows.append('|'.join([t['kind'], t['nm'], ','.join(t['el']), ','.join(f'{i}-{j}' for i, j in t['b']),
                              ','.join(str(i) for i, s in enumerate(t['sat']) if s), str(t['ih']), ','.join(t['loc']),
                              t.get('satName', '') if not t.get('satOnly') else '!' + t['nm']]))
    return rows


def main():
    tpls = learn()
    counts = {}
    for t in tpls: counts[t['kind']] = counts.get(t['kind'], 0) + 1
    print('templates:', counts, flush=True)
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)                 # (the full table, for inspection; not kept in git)
    with open(os.path.join(HERE, 'out', 'ring_parents.json'), 'w', encoding='utf-8') as f:
        json.dump(tpls, f, ensure_ascii=False, indent=0)
    rows = to_js(tpls)
    js = os.path.join(REPO, 'js', 'chem.js')
    s = open(js, encoding='utf-8').read()
    a, b = '/* RING_PARENTS:BEGIN', '/* RING_PARENTS:END */'
    if a not in s: raise SystemExit('js/chem.js has no RING_PARENTS markers')
    i, j = s.index(a), s.index(b)
    block = (a + ' (generated by tools/ring_parents/learn.py from OPSIN: do not edit by hand) */\n'
             '    const RING_PARENTS_RAW = ' + json.dumps(rows, ensure_ascii=False, separators=(',', ':')).replace('","', '",\n      "') + ';\n    ')
    s = s[:i] + block + s[j:]
    open(js, 'w', encoding='utf-8', newline='\n').write(s)
    print(f'js/chem.js: {len(rows)} ring parents, {len(block) / 1e3:.0f} kB')


if __name__ == '__main__':
    if hasattr(sys.stdout, 'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
    main()
