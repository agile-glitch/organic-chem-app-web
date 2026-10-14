"""Label validation, step 2 of 2: checks the values the Molecule tab's Labels modes show (tools/validate_labels.js dumps
them from the page) against references the app does not use:

  formal charge   RDKit Python's formal charge on the same atom
  lone pairs      (valence electrons - formal charge - bonds (Kekulé orders, H included) - unpaired electrons) / 2
  oxidation       RDKit Python CalcOxidationNumbers (its own implementation; it reads one Kekulé structure, so atoms
                  whose value depends on which Kekulé structure is drawn are counted apart)
  partial charge  RDKit Python Gasteiger-Marsili charges (12 iterations), the model the app ports
  hybridization   a textbook rule written here from scratch: steric number (σ bonds + lone pairs) 4/3/2 = sp3/sp2/sp,
                  sp2 for an atom whose lone pair is conjugated with a π bond on a neighbour or is part of an aromatic
                  ring; also the app's port of RDKit's own hybridization against RDKit Python
  R/S, E/Z        RDKit Python's new CIP labeller (rdCIPLabeler); unspecified centres from FindMolChiralCenters
  pKa             the app's rules alone (no measured value to fall back on) against the IUPAC Digitized pKa Dataset
                  (water), for molecules with one site of the kind and one measured value of the kind
  Taft, steric    tools/label_textbook.json (Taft table values, crowding series)
  every mode      tools/label_textbook.json: hand-checked textbook values

  python tools/validate_labels.py [tools/data/label_dump.json]   -> tools/validate_labels_report.md, tools/data/label_validation.json
"""
import json, os, statistics, sys
from collections import Counter, defaultdict
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdCIPLabeler, rdMolDescriptors

RDLogger.DisableLog('rdApp.*')
HERE = os.path.dirname(os.path.abspath(__file__))
SUP = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹', '0123456789')
PT = Chem.GetPeriodicTable()
MAX_EXAMPLES = 8


def norm_hyb(h):
    return (h or '').translate(SUP).replace(' ', '') or None


def lone_pairs(m, i):
    a = m.GetAtomWithIdx(i)
    if a.GetAtomicNum() == 1: return 0
    bonds = sum(b.GetBondTypeAsDouble() for b in a.GetBonds()) + a.GetTotalNumHs()
    v = PT.GetNOuterElecs(a.GetAtomicNum()) - a.GetFormalCharge() - bonds - a.GetNumRadicalElectrons()
    return int(v // 2) if v >= 0 and v == int(v) and v % 2 == 0 else None


PI_ATOMS = {6, 7, 8}


def p_pi(b):
    """a pπ-pπ bond: a double, triple or aromatic bond between C, N and O (an S=O or P=O is a polar single bond in
    reality, Reed & Schleyer 1990, so a lone pair next to one is not conjugated)"""
    return (b.GetBondTypeAsDouble() > 1 or b.GetIsAromatic()) and b.GetBeginAtom().GetAtomicNum() in PI_ATOMS and b.GetEndAtom().GetAtomicNum() in PI_ATOMS


def textbook_hyb(m, i, lp):
    """steric number rule; a lone pair next to a C/N/O π bond (or in an aromatic ring) makes the atom sp2; halogens
    and other terminal atoms with lone pairs are left out (their label is a convention: there is no angle to test)"""
    a = m.GetAtomWithIdx(i)
    if a.GetAtomicNum() in (1, 9, 17, 35, 53) or lp is None: return None
    sigma = a.GetDegree() + a.GetTotalNumHs()
    sn = sigma + lp
    if sn > 4: return None                               # hypervalent: no simple label
    base = {4: 'sp3', 3: 'sp2', 2: 'sp', 1: None}.get(sn)
    if base == 'sp3' and lp > 0:
        if a.GetIsAromatic(): return 'sp2'
        own_pi = any(b.GetBondTypeAsDouble() > 1 for b in a.GetBonds())
        nb_pi = any(any(p_pi(b) for b in n.GetBonds() if b.GetOtherAtomIdx(n.GetIdx()) != i) for n in a.GetNeighbors())
        if not own_pi and nb_pi: return 'sp2'
    return base


def environment(m, i):
    a = m.GetAtomWithIdx(i)
    nb_pi = any(any(p_pi(b) for b in n.GetBonds() if b.GetOtherAtomIdx(n.GetIdx()) != i) for n in a.GetNeighbors())
    return f'{a.GetSymbol()}{"(arom)" if a.GetIsAromatic() else ""} {a.GetDegree() + a.GetTotalNumHs()} σ-bonds' + \
        (f' charge {a.GetFormalCharge():+d}' if a.GetFormalCharge() else '') + (', next to a π bond' if nb_pi else '')


def kekule_dependent(m, i):
    """an atom with an aromatic bond to an atom of another element: its oxidation state depends on the Kekulé structure"""
    a = m.GetAtomWithIdx(i)
    return any(b.GetIsAromatic() and b.GetOtherAtom(a).GetAtomicNum() != a.GetAtomicNum() for b in a.GetBonds())


class Tally:
    def __init__(self, name, what):
        self.name, self.what, self.n, self.ok, self.classes = name, what, 0, 0, defaultdict(list)

    def add(self, ok, cls=None, example=None):
        self.n += 1
        if ok: self.ok += 1
        else: self.classes[cls].append(example)

    def summary(self):
        return {'name': self.name, 'what': self.what, 'checked': self.n, 'agree': self.ok,
                'disagreements': {k: {'count': len(v), 'examples': v[:MAX_EXAMPLES]} for k, v in sorted(self.classes.items(), key=lambda kv: -len(kv[1]))}}


def main():
    dump_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'data', 'label_dump.json')
    D = json.load(open(dump_path))['molecules']
    T = json.load(open(os.path.join(HERE, 'label_textbook.json')))
    s = open(os.path.join(HERE, '..', 'data', 'pka.js')).read()
    MEAS = json.loads(s[s.index('{'):s.rindex('}') + 1])['molecules']

    t = {k: Tally(k, w) for k, w in [
        ('formal', 'formal charge vs RDKit Python'), ('lonepairs', 'lone pairs vs the electron count'),
        ('oxidation', 'oxidation state vs RDKit CalcOxidationNumbers (atoms whose value does not depend on the Kekulé structure)'),
        ('oxidation_kekule', 'oxidation state, atoms whose value depends on the Kekulé structure (the app averages; RDKit takes one structure)'),
        ('partial', 'Gasteiger-Marsili charge vs RDKit Python (difference below 0.005; the label shows 2 decimals)'),
        ('hybridization', 'displayed hybridization vs the textbook steric-number rule'),
        ('rdkit_hybridization', "the app's port of RDKit's hybridization vs RDKit Python"),
        ('cip', 'R/S vs RDKit Python CIP labeller (atoms either one labels)'),
        ('ez', 'E/Z vs RDKit Python CIP labeller (bonds either one labels)')]}
    analysed = failed = 0
    fails = []
    for smi, r in D.items():
        if r.get('error'):
            failed += 1; fails.append((smi, r['error'])); continue
        m = Chem.MolFromSmiles(smi)
        if m is None or m.GetNumAtoms() != r['n']:
            failed += 1; fails.append((smi, 'RDKit Python reads a different atom count')); continue
        analysed += 1
        mk = Chem.Mol(m)
        try: Chem.Kekulize(mk, clearAromaticFlags=False)
        except Exception: mk = m
        try: rdMolDescriptors.CalcOxidationNumbers(m); have_ox = True
        except Exception: have_ox = False
        mq = Chem.Mol(m)
        try: AllChem.ComputeGasteigerCharges(mq, nIter=12); have_q = True
        except Exception: have_q = False
        rdCIPLabeler.AssignCIPLabels(m)
        for i in range(r['n']):
            a = m.GetAtomWithIdx(i)
            ex = f'{smi} atom {i} ({a.GetSymbol()})'
            t['formal'].add(r['charge'][i] == a.GetFormalCharge(), 'different', f'{ex}: app {r["charge"][i]}, RDKit {a.GetFormalCharge()}')
            lp = lone_pairs(mk, i)
            if lp is not None or r['lp'][i] is not None:
                t['lonepairs'].add(lp == r['lp'][i], 'unusual electron count' if lp is None or r['lp'][i] is None else environment(m, i), f'{ex}: app {r["lp"][i]}, count {lp}')
            if have_ox and r['ox'][i] is not None:
                ref = int(a.GetProp('OxidationNumber'))
                key = 'oxidation_kekule' if kekule_dependent(m, i) else 'oxidation'
                t[key].add(abs(r['ox'][i] - ref) < 1e-6, environment(m, i), f'{ex}: app {r["ox"][i]}, RDKit {ref}')
            if have_q and r['q'][i] is not None:
                ref = float(mq.GetAtomWithIdx(i).GetProp('_GasteigerCharge'))
                if ref == ref:
                    t['partial'].add(abs(r['q'][i] - ref) < 0.005, 'different', f'{ex}: app {r["q"][i]:.4f}, RDKit {ref:.4f}')
            th = textbook_hyb(mk, i, lp)
            ah = norm_hyb(r['hyb'][i])
            if th is not None and ah is not None:
                t['hybridization'].add(th == ah, f'{environment(m, i)}: app {ah}, rule {th}', ex)
            rh = str(a.GetHybridization()).replace('SP3D2', 'sp3d2').replace('SP3D', 'sp3d').replace('SP3', 'sp3').replace('SP2', 'sp2').replace('SP', 'sp')
            if r['rdHyb'][i] is not None:
                t['rdkit_hybridization'].add(str(r['rdHyb'][i]).lower() == rh.lower(), f'app {r["rdHyb"][i]}, RDKit {rh}', ex)
            ref_cip = a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None
            app_cip = r['cip'][i] if r['cip'][i] not in ('?', None) else None
            if ref_cip or app_cip:
                t['cip'].add(ref_cip == app_cip, f'app {app_cip}, RDKit {ref_cip}', ex)
        app_ez = {tuple(sorted(x[:2])): x[2] for x in r['ez'] if x[2] in ('E', 'Z')}
        ref_ez = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))): b.GetProp('_CIPCode') for b in m.GetBonds() if b.HasProp('_CIPCode') and b.GetProp('_CIPCode') in ('E', 'Z')}
        for k in set(app_ez) | set(ref_ez):
            t['ez'].add(app_ez.get(k) == ref_ez.get(k), f'app {app_ez.get(k)}, RDKit {ref_ez.get(k)}', f'{smi} bond {k}')

    # ---- pKa: the rules alone against measured values
    pk = {'acid': [], 'base': []}
    missed = {'acid': Counter(), 'base': Counter()}
    shown_err = []
    for smi, r in D.items():
        if r.get('error') or smi not in MEAS or not r.get('rules'): continue
        M = MEAS[smi]
        n_vals = len(M.get('acid') or []) + len([v for v in M.get('base') or [] if v >= 0])
        n_sites = len(r['rules']['acid']) + len(r['rules']['base'])
        for kind in ('acid', 'base'):
            vals = [v for v in M.get(kind) or [] if kind == 'acid' or v >= 0]
            sites = r['rules'][kind]
            if vals and not sites: missed[kind][M.get('quality', '?')] += 1
            # one ionisable group in all (one site, one measured value, the same kind): the value belongs to that site
            if len(vals) == 1 and len(sites) == 1 and n_vals == 1 and n_sites == 1:
                s_ = sites[0]
                method = 'Hammett estimate' if s_['source'].startswith('estimated') else 'typical value'
                pk[kind].append({'smiles': smi, 'measured': vals[0], 'rules': s_['pKa'], 'err': s_['pKa'] - vals[0], 'method': method,
                                 'fragments': '+'.join(sorted(set(s_['fragments']))), 'quality': M.get('quality')})
        for kind in ('acid', 'base'):
            for s_ in r['pka'][kind]:
                if s_['source'].startswith('measured'): shown_err.append(0)

    def stats(rows):
        if not rows: return None
        e = [abs(x['err']) for x in rows]
        return {'n': len(rows), 'median_abs_err': round(statistics.median(e), 2), 'within_1': round(sum(x <= 1 for x in e) / len(e), 3),
                'within_2': round(sum(x <= 2 for x in e) / len(e), 3), 'mean_signed_err': round(statistics.mean(x['err'] for x in rows), 2)}
    pka_summary = {}
    for kind in ('acid', 'base'):
        rows = pk[kind]
        by_frag = defaultdict(list)
        for x in rows: by_frag[x['fragments']].append(x)
        pka_summary[kind] = {
            'all': stats(rows),
            'by_method': {mth: stats([x for x in rows if x['method'] == mth]) for mth in ('Hammett estimate', 'typical value')},
            'by_fragment': {f: dict(stats(v), name=f) for f, v in sorted(by_frag.items(), key=lambda kv: -len(kv[1])) if len(v) >= 10},
            'worst': sorted(rows, key=lambda x: -abs(x['err']))[:15],
            'molecules_with_a_measured_value_and_no_site': dict(missed[kind]),
        }

    # ---- textbook list
    tb_rows = []
    for c in T['cases']:
        r = D.get(c['smiles'])
        for e in c['expect']:
            row = {'name': c['name'], 'smiles': c['smiles'], 'mode': e['mode'], 'expected': None, 'app': None, 'ok': False, 'why': e.get('why', '')}
            if not r or r.get('error'):
                row['app'] = 'not analysed: ' + (r or {}).get('error', 'missing from the dump'); tb_rows.append(row); continue
            i = e.get('atom')
            md = e['mode']
            if md == 'hybridization': row['expected'], row['app'] = e['value'], norm_hyb(r['hyb'][i])
            elif md == 'formal': row['expected'], row['app'] = e['value'], r['charge'][i]
            elif md == 'lonepairs': row['expected'], row['app'] = e['value'], r['lp'][i]
            elif md == 'oxidation': row['expected'], row['app'] = e['value'], r['ox'][i]
            elif md == 'cip': row['expected'], row['app'] = e['value'], r['cip'][i]
            elif md == 'ez':
                row['expected'] = e['value']
                row['app'] = next((x[2] for x in r['ez'] if sorted(x[:2]) == sorted(e['bond'])), None)
            elif md == 'pka':
                kind = e.get('kind', 'acid')
                site = next((s_ for s_ in r['pka'][kind] if i in s_['atoms']), None)
                row['app'] = None if site is None else f'{site["pKa"]} ({site["source"][:60]})'
                if e.get('none'):
                    row['expected'] = 'no label'; row['ok'] = site is None; tb_rows.append(row); continue
                row['expected'] = f'{e["range"][0]}–{e["range"][1]}'
                row['ok'] = site is not None and e['range'][0] <= site['pKa'] <= e['range'][1]; tb_rows.append(row); continue
            elif md == 'taft':
                eff = next((x for x in r['effects'] if x['atom'] == i), None)
                row['expected'] = f'σ* {e["sigmaStar"]:+.2f} / Es {e["Es"]:+.2f}'
                row['app'] = None if not eff else f'σ* {eff["sigmaStar"]} / Es {eff["EsMin"]}'
                row['ok'] = bool(eff) and eff['sigmaStar'] is not None and eff['EsMin'] is not None and abs(eff['sigmaStar'] - e['sigmaStar']) < 0.006 and abs(eff['EsMin'] - e['Es']) < 0.006
                tb_rows.append(row); continue
            row['ok'] = row['app'] == row['expected'] if not isinstance(row['expected'], (int, float)) or row['app'] is None else abs(row['app'] - row['expected']) < 1e-6
            tb_rows.append(row)
    series_rows = []
    for sr in T['steric_series']['series']:
        scores = []
        for smi in sr['smiles']:
            r = D.get(smi) or {}
            eff = [x for x in r.get('effects') or [] if x['fragment'] != 'F009']
            # the reacting atom of the molecule's main group: the C=O carbon, else the C–Br carbon, else the N
            pick = next((x for x in eff if x['atom'] in [i for i, el in enumerate(r.get('el') or []) if el == 'C']), None) if 'N' not in smi else next((x for x in eff if (r.get('el') or [])[x['atom']] == 'N'), None)
            scores.append((smi, pick['score'] if pick else None, pick['steric'] if pick else None))
        ok = all(s1 is not None for _, s1, _ in scores) and all(scores[k][1] <= scores[k + 1][1] for k in range(len(scores) - 1))
        series_rows.append({'name': sr['name'], 'ok': ok, 'molecules': scores})

    display = check_display(os.path.join(HERE, 'data', 'label_display.json'))
    out = {'dump': os.path.relpath(dump_path, HERE), 'display': display, 'analysed': analysed, 'not_analysed': failed, 'not_analysed_examples': fails[:15],
           'checks': {k: v.summary() for k, v in t.items()}, 'pka': pka_summary, 'textbook': tb_rows, 'steric_series': series_rows}
    os.makedirs(os.path.join(HERE, 'data'), exist_ok=True)
    json.dump(out, open(os.path.join(HERE, 'data', 'label_validation.json'), 'w'), indent=1, ensure_ascii=False)
    write_report(out)


def angle(X, i, j, k):
    import math
    a = [X[3 * i + t] - X[3 * j + t] for t in range(3)]; b = [X[3 * k + t] - X[3 * j + t] for t in range(3)]
    d = sum(x * y for x, y in zip(a, b)) / (math.sqrt(sum(x * x for x in a)) * math.sqrt(sum(x * x for x in b)))
    return math.degrees(math.acos(max(-1, min(1, d))))


def check_display(path):
    """tools/validate_display.js: labels left out, nearer another atom than their own, inside a ring, far from their
    atom (distances in bond lengths); bond angles of the app's 3D model against RDKit ETKDG + MMFF94"""
    if not os.path.exists(path): return None
    D = json.load(open(path))['molecules']
    modes = defaultdict(lambda: Counter())
    examples = defaultdict(list)
    for smi, r in D.items():
        if r.get('error'): continue
        for md, x in r['modes'].items():
            c = modes[md]
            c['molecules'] += 1; c['drawn'] += x['drawn']; c['left_out'] += x['leftOut']
            if x['leftOut']: c['molecules_with_left_out'] += 1
            for nt in x['notes']:
                amb = nt['nearest'] != nt['atom']
                c['ambiguous'] += amb; c['in_ring'] += nt['inRing']; c['far'] += nt['ownDist'] > 1.0
                if amb and len(examples[md]) < 6: examples[md].append(f"{smi}: '{nt['text']}' for atom {nt['atom']} sits nearer atom {nt['nearest']} ({nt['ownDist']} vs {nt['nearDist']} bond lengths)")
    # bond angles
    from rdkit.Chem import AllChem as AC
    diffs, worst, n_mol, notes_left = [], [], 0, Counter()
    for smi, r in D.items():
        A = r.get('angles')
        if not A or A.get('error') or not A.get('optimised'): continue
        m = Chem.AddHs(Chem.MolFromSmiles(smi))
        if AC.EmbedMolecule(m, randomSeed=0xF00D) != 0: continue
        AC.MMFFOptimizeMolecule(m, maxIters=5000)
        X = [v for i in range(m.GetNumAtoms()) for v in m.GetConformer().GetAtomPosition(i)]
        n_mol += 1
        heavy = A['n']
        if '(' in A.get('foot', '') or 'left out' in A.get('foot', ''): notes_left[A['foot'][:80]] += 1
        for j in range(heavy):
            nb = [x.GetIdx() for x in m.GetAtomWithIdx(j).GetNeighbors() if x.GetIdx() < heavy]
            for p_ in range(len(nb)):
                for q in range(p_ + 1, len(nb)):
                    a1, a2 = angle(A['coords'], nb[p_], j, nb[q]), angle(X, nb[p_], j, nb[q])
                    diffs.append(abs(a1 - a2)); worst.append((abs(a1 - a2), smi, (nb[p_], j, nb[q]), round(a1, 1), round(a2, 1)))
    worst.sort(reverse=True)
    ang = None
    if diffs:
        ds = sorted(diffs)
        ang = {'molecules': n_mol, 'angles': len(ds), 'median': round(statistics.median(ds), 2), 'p95': round(ds[int(0.95 * (len(ds) - 1))], 2),
               'max': round(ds[-1], 2), 'over_3': sum(d > 3 for d in ds), 'over_5': sum(d > 5 for d in ds),
               'worst': [{'smiles': w[1], 'atoms': w[2], 'app': w[3], 'rdkit': w[4]} for w in worst[:8]]}
    return {'modes': {k: dict(v) for k, v in modes.items()}, 'ambiguous_examples': dict(examples), 'angles': ang}


def pct(a, n): return f'{100 * a / n:.2f} %' if n else '—'


def write_report(o):
    L = ['# Label validation (Molecule tab)', '',
         f'Made by `tools/validate_labels.js` + `tools/validate_labels.py` from `{o["dump"]}`: {o["analysed"]:,} molecules analysed, {o["not_analysed"]} not.', '']
    if os.path.exists(os.path.join(HERE, 'validate_labels_findings.md')):
        L += [open(os.path.join(HERE, 'validate_labels_findings.md')).read().rstrip(), '']
    L += [
         '## Agreement with independent references', '', '| Label | Reference | Checked | Agree | |', '|---|---|---:|---:|---:|']
    for k, c in o['checks'].items():
        L.append(f'| {k} | {c["what"]} | {c["checked"]:,} | {c["agree"]:,} | {pct(c["agree"], c["checked"])} |')
    for k, c in o['checks'].items():
        if not c['disagreements']: continue
        L += ['', f'### {k}: disagreements by kind', '']
        for cls, d in list(c['disagreements'].items())[:12]:
            L.append(f'- **{cls}**: {d["count"]:,} — e.g. ' + '; '.join(f'`{x}`' for x in d['examples'][:3]))
    L += ['', '## pKa: the rules alone against measured values (IUPAC, water)', '',
          'Molecules with one ionisable group in all: the app finds one site (acid or base) and the dataset has one measured value, of the same kind, so the value belongs to that site. The app works its value out with no measured value to fall back on (on screen a measured value, when there is one, is shown instead).', '']
    for kind in ('acid', 'base'):
        P = o['pka'][kind]
        if not P['all']: continue
        a = P['all']
        L += [f'**{kind}**: {a["n"]:,} molecules; median error {a["median_abs_err"]}, within 1 unit {a["within_1"]:.0%}, within 2 units {a["within_2"]:.0%}, mean signed error {a["mean_signed_err"]:+}.', '']
        L += ['| | n | median error | within 1 | within 2 | mean signed |', '|---|---:|---:|---:|---:|---:|']
        for mth, st in P['by_method'].items():
            if st: L.append(f'| {mth} | {st["n"]:,} | {st["median_abs_err"]} | {st["within_1"]:.0%} | {st["within_2"]:.0%} | {st["mean_signed_err"]:+} |')
        for f, st in list(P['by_fragment'].items())[:15]:
            L.append(f'| {f} | {st["n"]:,} | {st["median_abs_err"]} | {st["within_1"]:.0%} | {st["within_2"]:.0%} | {st["mean_signed_err"]:+} |')
        L += ['', 'Largest errors: ' + '; '.join(f'`{x["smiles"]}` measured {x["measured"]}, rules {x["rules"]} ({x["fragments"]})' for x in P['worst'][:6]) + '.',
              f'Molecules with a measured {kind} value but no {kind} site found: ' + ', '.join(f'{v} ({q})' for q, v in P['molecules_with_a_measured_value_and_no_site'].items()) + '.', '']
    Dp = o.get('display')
    if Dp:
        L += ['## The drawing (`tools/validate_display.js`)', '',
              'Each label mode on the 2D drawing, as the tab shows it (Show H off). *Left out*: labels with a value that were not drawn ("no free space"). *Nearer another atom*: the label\'s centre is closer to another atom than to its own, so a reader can take it for that atom. *Inside a ring*: the label sits inside a ring. *Far*: more than one bond length from its atom.', '',
              '| Mode | Molecules | Labels drawn | Left out | Molecules with labels left out | Nearer another atom | Inside a ring | Far |', '|---|---:|---:|---:|---:|---:|---:|---:|']
        for md, c in Dp['modes'].items():
            dr = c.get('drawn', 0)
            L.append(f"| {md} | {c.get('molecules', 0)} | {dr:,} | {c.get('left_out', 0):,} ({pct(c.get('left_out', 0), dr + c.get('left_out', 0))}) | {c.get('molecules_with_left_out', 0)} | "
                     f"{c.get('ambiguous', 0):,} ({pct(c.get('ambiguous', 0), dr)}) | {c.get('in_ring', 0):,} ({pct(c.get('in_ring', 0), dr)}) | {c.get('far', 0):,} |")
        L += ['', 'Examples of labels nearer another atom: ' + '; '.join(f'`{x}`' for md in ('partial', 'index') for x in Dp['ambiguous_examples'].get(md, [])[:3]), '']
        a = Dp.get('angles')
        if a:
            L += ['### Bond angles: the app\'s 3D model against RDKit (ETKDG + MMFF94)', '',
                  f"{a['angles']:,} angles between heavy atoms in {a['molecules']} molecules: median difference {a['median']}°, 95 % within {a['p95']}°, largest {a['max']}°; {a['over_3']} differ by more than 3°, {a['over_5']} by more than 5°.", '',
                  'Largest: ' + '; '.join(f"`{w['smiles']}` atoms {w['atoms']}: app {w['app']}°, RDKit {w['rdkit']}°" for w in a['worst'][:5]) + '.', '']
    n_ok = sum(r['ok'] for r in o['textbook'])
    L += ['## Textbook list (`tools/label_textbook.json`)', '', f'{n_ok} of {len(o["textbook"])} hand-checked values agree.', '']
    bad = [r for r in o['textbook'] if not r['ok']]
    if bad:
        L += ['| Molecule | Label | Expected | App | Note |', '|---|---|---|---|---|']
        for r in bad: L.append(f'| {r["name"]} `{r["smiles"]}` | {r["mode"]} | {r["expected"]} | {r["app"]} | {r["why"]} |')
        L.append('')
    L += ['## Steric class along crowding series', '']
    for s_ in o['steric_series']:
        L.append(f'- {s_["name"]}: {"ok" if s_["ok"] else "**out of order**"} — ' + ', '.join(f'`{a}` {c} ({b})' for a, b, c in s_['molecules']))
    open(os.path.join(HERE, 'validate_labels_report.md'), 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L[:20]))


if __name__ == '__main__':
    main()
