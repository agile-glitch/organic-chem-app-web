"""Read an AGILES code (AGILES1:... or the older ACM1:...) and print everything in it: the identifiers, every atom with its
element, name (C1, O2 ...), 2D position and its label in every Labels mode, the bonds, the bond angles and the notes.
Standard library only (base64, zlib); docs/AGILES1_format.md describes the format.

  python tools/agiles_read.py "AGILES1:..."        (or pipe codes in, one per line; --json for machine-readable output)
"""
import base64, json, sys, zlib

MODES = ['cip', 'hybridization', 'formal', 'partial', 'oxidation', 'lonepairs', 'pka', 'steric', 'electronic', 'taft']
IDS = ['name', 'iupac', 'iupac_checked', 'smiles', 'smiles_rdkit', 'inchi', 'inchikey', 'cas', 'formula', 'isotopic_formula', 'mw', 'monoisotopic']


def unpack(code):
    code = code.strip()
    prefix, _, body = code.partition(':')
    if prefix not in ('AGILES1', 'ACM1') or not body: raise ValueError('not an AGILES code: it must start with AGILES1: (or ACM1:)')
    return zlib.decompress(base64.urlsafe_b64decode(body + '=' * (-len(body) % 4)), -15).decode('utf-8')   # raw DEFLATE


def parse(sdf):
    sdf = sdf.replace('\r\n', '\n')
    molblock, rest = sdf.split('M  END', 1)
    lines = molblock.split('\n')
    title = lines[0]
    counts = lines[3]
    if 'V3000' in counts: raise ValueError('V3000 MOL block (over 999 atoms): read it with RDKit')
    na, nb = int(counts[0:3]), int(counts[3:6])
    atoms = [{'x': float(l[0:10]), 'y': float(l[10:20]), 'el': l[31:34].strip(), 'chg_field': int(l[36:39] or 0)} for l in lines[4:4 + na]]
    bonds = [{'a': int(l[0:3]), 'b': int(l[3:6]), 'order': int(l[6:9]), 'stereo': int(l[9:12] or 0)} for l in lines[4 + na:4 + na + nb]]
    charges = {}
    for l in lines[4 + na + nb:]:                          # "M  CHG" lines: the formal charges
        if l.startswith('M  CHG'):
            f = l.split()[3:]
            for i in range(0, len(f), 2): charges[int(f[i])] = int(f[i + 1])
    items, cur = {}, None
    for l in rest.split('\n'):
        if l.startswith('$$$$'): break
        if l.startswith('>') and '<' in l:
            cur = l[l.index('<') + 1:l.index('>', l.index('<'))]
            cur = cur.replace('ACM.', 'AGILES.', 1)        # before the rename
            items.setdefault(cur, [])
        elif cur and l.strip(): items[cur].append(l)
        else: cur = None
    one = lambda k: ' '.join(items.get('AGILES.' + k, [])).strip() or None
    nh, nH = map(int, one('atoms').split())
    N = nh + nH
    # the atom list: the MOL block's atoms (1..nh), then one H per implicit H, grouped by parent in ascending parent order
    # (the parents are stored in AGILES.h_parents; for older codes they follow from each heavy atom's valence, as below)
    # each heavy atom's H count from its valence (the smallest allowed valence at or above its bond-order sum, shifted by
    # its charge as RDKit does); the H are then numbered parent by parent. Checked: the counts must add up to nH.
    VAL = {'H': [1], 'B': [3], 'C': [4], 'N': [3, 5], 'O': [2], 'F': [1], 'Si': [4], 'P': [3, 5], 'S': [2, 4, 6], 'Cl': [1, 3, 5, 7],
           'Se': [2, 4, 6], 'Br': [1, 3, 5], 'I': [1, 3, 5, 7]}
    ISO = {1: {'N': 'C', 'O': 'N', 'S': 'P', 'C': 'B', 'B': 'Be'}, -1: {'C': 'N', 'N': 'O', 'O': 'F', 'S': 'Cl', 'B': 'C'}}   # isoelectronic element
    bo = [0] * nh
    for b_ in bonds:
        o = {1: 1, 2: 2, 3: 3, 4: 1.5}.get(b_['order'], 0)
        if b_['order'] == 9: o = 0                         # dative
        bo[b_['a'] - 1] += o; bo[b_['b'] - 1] += o
    rad = {}
    for l in lines[4 + na + nb:]:
        if l.startswith('M  RAD'):
            f = l.split()[3:]
            for i in range(0, len(f), 2): rad[int(f[i])] = {1: 2, 2: 1, 3: 2}.get(int(f[i + 1]), 0)
    hcount = []
    for i, a_ in enumerate(atoms):
        c, el = charges.get(i + 1, 0), a_['el']
        if c and abs(c) == 1 and el in ISO[c]: el = ISO[c][el]
        vals = VAL.get(el)
        if not vals: hcount.append(0); continue
        used = bo[i] + rad.get(i + 1, 0)
        v = next((x for x in vals if x >= used), None)
        hcount.append(max(0, int(round(v - used))) if v is not None else 0)
    parents = [i + 1 for i, k in enumerate(hcount) for _ in range(k)] if sum(hcount) == nH else None
    if one('h_parents'): parents = [int(x) for x in one('h_parents').split()]   # stored since 2026-10-15: no need to work it out
    names, count = [], {}
    els = [a['el'] for a in atoms] + ['H'] * nH
    for e in els: count[e] = count.get(e, 0) + 1; names.append(f'{e}{count[e]}')
    labels = {}
    for m in MODES:
        if 'AGILES.labels.' + m not in items: continue
        arr = [''] * N
        for l in items['AGILES.labels.' + m]:
            k, _, txt = l.strip().partition(' ')
            arr[int(k) - 1] = txt
        labels[m] = arr
    angles = [dict(zip(('centre', 'i', 'j'), map(int, l.split()[:3])), text=l.split(None, 3)[3]) for l in items.get('AGILES.labels.angles', [])]
    return {'title': title, 'format': one('format'), 'app': one('app'), 'ids': {k: one(k) for k in IDS if one(k)},
            'heavy_atoms': nh, 'implicit_H': nH,
            'atoms': [{'n': i + 1, 'name': names[i], 'el': els[i], **({'x': atoms[i]['x'], 'y': atoms[i]['y'], 'charge': charges.get(i + 1, 0)} if i < nh
                                                                    else {'parent': parents[i - nh] if parents else None}),
                       'labels': {m: labels[m][i] for m in labels if labels[m][i]}} for i in range(N)],
            'bonds': bonds, 'angles': angles, 'notes': [' '.join(items[k]) for k in items if k == 'AGILES.note'] or [], 'modes_present': list(labels)}


def show(d):
    print(f"{d['title']}   ({d['format']}, written by {d['app']})")
    for k, v in d['ids'].items(): print(f'  {k:14s} {v}')
    print(f"  {d['heavy_atoms']} heavy atoms (with 2D coordinates), then {d['implicit_H']} H")
    cols = d['modes_present']
    print('  ' + ' | '.join(['atom', 'name', 'x', 'y'] + cols))
    for a in d['atoms']:
        xy = [f"{a['x']:.3f}", f"{a['y']:.3f}"] if 'x' in a else ['on ' + (d['atoms'][a['parent'] - 1]['name'] if a.get('parent') else '?'), '']
        print('  ' + ' | '.join([str(a['n']), a['name']] + xy + [a['labels'].get(m, '') for m in cols]))
    print('  bonds (atom-atom order; order 1 single, 2 double, 3 triple, 9 dative; stereo 1 wedge, 6 hash, 4 either (single bond), 3 either (double bond: E/Z not given), narrow end at the first atom):')
    print('   ' + ', '.join(f"{b['a']}-{b['b']} {b['order']}" + (f" s{b['stereo']}" if b['stereo'] else '') for b in d['bonds']))
    if d['angles']: print('  bond angles (centre: i-centre-j):', ', '.join(f"{x['centre']}: {x['i']}-{x['centre']}-{x['j']} {x['text']}" for x in d['angles']))
    for n in d['notes']: print('  note:', n)


if __name__ == '__main__':
    codes = [a for a in sys.argv[1:] if not a.startswith('--')] or [l for l in sys.stdin.read().split() if ':' in l]
    for c in codes:
        d = parse(unpack(c))
        print(json.dumps(d, ensure_ascii=False, indent=1)) if '--json' in sys.argv else show(d)
        print()
