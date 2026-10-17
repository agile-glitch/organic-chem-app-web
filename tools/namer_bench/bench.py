#!/usr/bin/env python3
"""How well the app's IUPAC namer does, checked by an independent name reader: every molecule is named by the namer
(js/chem.js, in Node.js), each name is read back by OPSIN (the open-source IUPAC name reader), and RDKit compares what
OPSIN read with the molecule (full InChIKey, so stereo counts). Each molecule ends up as one of:

    ok          the name reads back to the same molecule, stereo included
    stereo      the right skeleton, but the stereo read back differs (or is missing)
    wrong       OPSIN reads the name as a DIFFERENT molecule (the worst outcome: a wrong name)
    unreadable  OPSIN cannot read the name
    unnamed     the namer gives no name (its reason is counted)

    python tools/namer_bench/bench.py [library] [--smiles FILE] [--codes FILE] [--sample N] [--seed S] [--show KIND]

library: the AGILES spreadsheet's molecules (AGILES_codes.xlsx, named from each code's drawing). --smiles FILE: one
SMILES per line (a tab, then anything). --codes FILE: one AGILES1 code per line. --show KIND: list the molecules of a kind
(ok, stereo, wrong, unreadable, unnamed). A report of every molecule goes to tools/namer_bench/out/<set>[.TAG].json
(--tag TAG); CHEM_JS=path names with another copy of chem.js (e.g. the one before a change, to compare).
Needs Node.js, Java (JAVA, else java on PATH), OPSIN (OPSIN_JAR: the opsin-cli jar-with-dependencies) and RDKit.
"""
import base64, collections, glob, json, os, random, re, subprocess, sys, zipfile, zlib
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
OUT = os.path.join(HERE, 'out')


def java():
    if os.environ.get('JAVA'): return os.environ['JAVA']
    for pat in (r'C:\Program Files\Eclipse Adoptium\jre-*\bin\java.exe', r'C:\Program Files\Eclipse Adoptium\jdk-*\bin\java.exe'):
        hits = sorted(glob.glob(pat))
        if hits: return hits[-1]
    return 'java'


def opsin_jar():
    if os.environ.get('OPSIN_JAR'): return os.environ['OPSIN_JAR']
    hits = sorted(glob.glob(os.path.join(os.path.expanduser('~'), 'tools', 'opsin', 'opsin-cli-*-jar-with-dependencies.jar')))
    if not hits: sys.exit('OPSIN not found: set OPSIN_JAR to the opsin-cli jar-with-dependencies.')
    return hits[-1]


def unpack(code):
    body = code.split(':', 1)[1]
    return zlib.decompress(base64.urlsafe_b64decode(body + '=' * (-len(body) % 4)), -15).decode('utf-8')


def molblock(code):
    return unpack(code).split('M  END', 1)[0] + 'M  END\n'


def library():
    """(id, molfile) of every molecule in the AGILES spreadsheet"""
    ns = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'
    out = []
    with zipfile.ZipFile(os.path.join(REPO, 'AGILES_codes.xlsx')) as z, z.open('xl/worksheets/sheet1.xml') as f:
        for _, el in ET.iterparse(f):
            if el.tag != ns + 'row': continue
            row = {}
            for c in el:
                col = ''.join(ch for ch in c.get('r', '') if ch.isalpha())
                t = c.find(f'{ns}is/{ns}t')
                v = t.text if t is not None else c.findtext(ns + 'v')
                if v: row[col] = v.strip()
            el.clear()
            if row.get('B', '').startswith('AGILES1:'): out.append((row.get('A') or row['B'][:20], {'molfile': molblock(row['B'])}))
    return out


def name_all(items):
    """{id: (name or None, error or None)} from the namer"""
    p = subprocess.run(['node', os.path.join(HERE, 'namer.js')], input='\n'.join(json.dumps({'id': i, **req}) for i, (_, req) in enumerate(items)) + '\n',
                       capture_output=True, text=True, encoding='utf-8')
    out = {}
    for line in p.stdout.splitlines():
        d = json.loads(line); out[d['id']] = (d.get('name'), (d.get('why') and ('namer: ' + d['why'])) or d.get('error'))
    return out


def opsin(names):
    """OPSIN's SMILES of each name ('' where it cannot read it), in order"""
    if not names: return []
    p = subprocess.run([java(), '-jar', opsin_jar(), '-osmi'], input='\n'.join(names) + '\n', capture_output=True, text=True, encoding='utf-8')
    lines = p.stdout.replace('\r', '').split('\n')
    return (lines + [''] * len(names))[:len(names)]


def main(argv):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
    sample, seed, show = None, 1, None
    if '--sample' in argv: sample = int(argv[argv.index('--sample') + 1])
    if '--seed' in argv: seed = int(argv[argv.index('--seed') + 1])
    if '--show' in argv: show = argv[argv.index('--show') + 1]
    tag = argv[argv.index('--tag') + 1] if '--tag' in argv else ''
    if '--smiles' in argv:
        path = argv[argv.index('--smiles') + 1]; label = os.path.splitext(os.path.basename(path))[0]
        items = [(l.split('\t')[0].strip(), {'smiles': l.split('\t')[0].strip()}) for l in open(path, encoding='utf-8') if l.strip() and not l.startswith('#')]
    elif '--codes' in argv:
        path = argv[argv.index('--codes') + 1]; label = os.path.splitext(os.path.basename(path))[0]
        items = [(l.strip()[:30], {'molfile': molblock(l.strip())}) for l in open(path, encoding='utf-8') if l.strip().startswith('AGILES1:')]
    else:
        label, items = 'library', library()
    if sample and sample < len(items):
        random.Random(seed).shuffle(items); items = items[:sample]
    names = name_all(items)
    ids = [i for i in range(len(items)) if names.get(i, (None,))[0]]
    smis = dict(zip(ids, opsin([names[i][0] for i in ids])))
    report, count, reasons = [], collections.Counter(), collections.Counter()
    for i, (key, req) in enumerate(items):
        nm, err = names.get(i, (None, 'no answer'))
        ref = Chem.MolFromMolBlock(req['molfile']) if 'molfile' in req else Chem.MolFromSmiles(req['smiles'])
        refkey = Chem.MolToInchiKey(ref) if ref is not None else ''
        r = {'key': key, 'name': nm, 'ref_smiles': Chem.MolToSmiles(ref) if ref is not None else ''}
        if not nm:
            kind = 'unnamed'; r['reason'] = err; reasons[re.sub(r'\d+', 'N', err or '')[:100]] += 1
        else:
            smi = smis.get(i, '')
            back = Chem.MolFromSmiles(smi) if smi else None
            if back is None: kind = 'unreadable'
            else:
                k = Chem.MolToInchiKey(back); r['opsin_smiles'] = smi
                kind = 'ok' if k == refkey else 'stereo' if k[:14] == refkey[:14] else 'wrong'
        r['kind'] = kind; count[kind] += 1; report.append(r)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, label + (('.' + tag) if tag else '') + '.json'), 'w', encoding='utf-8') as f: json.dump(report, f, ensure_ascii=False, indent=0)
    n = len(items)
    print(f'{label}: {n} molecules')
    for k in ('ok', 'stereo', 'wrong', 'unreadable', 'unnamed'):
        print(f'  {k:11s} {count[k]:6d}  {100 * count[k] / max(n, 1):5.1f}%')
    if reasons:
        print('  why unnamed:'); [print(f'    {v:6d}  {k}') for k, v in reasons.most_common(12)]
    if show:
        for r in report:
            if r['kind'] == show: print(r['ref_smiles'], '|', r.get('name') or r.get('reason'), '|', r.get('opsin_smiles', ''))


if __name__ == '__main__':
    if hasattr(sys.stdout, 'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
    main(sys.argv[1:])
