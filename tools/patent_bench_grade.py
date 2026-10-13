"""Grade a patent benchmark run (tools/patent_bench.js) with RDKit and write tools/patent_bench_report.json and tools/patent_bench_report.md.

    python tools/patent_bench_grade.py [cases.json] [results.json]

For each recorded reaction, every molecule of the recorded product must be among the products the page showed. Three readings, strictest first:
  exact     canonical SMILES with stereo
  flat      the same without stereo (R/S, E/Z)
  compound  the same compound: the InChIKey's first block (the skeleton; a salt or protonated form counts, e.g. R-NH3+ Br- for R-NH2)
What the page did with each one:
  answered  it showed products (graded right or wrong at each reading)
  declined  it showed none: "Unknown" (nothing it knows covers these molecules, or the model was not sure enough) or "No reaction"
  error     the page could not take the code or did not finish
Which part of the page answered is read from the first line of its reading: a workbook rule (R001..., SN1/SN2/E1/E2), or the relevance model
(learned from the patents of the training split). The model's stated chance is read from the same line, to check it against how often it was right.
"""
import collections, json, os, re, sys

from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = os.path.dirname(os.path.abspath(__file__))
CASES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'patent_bench_cases.json')
RESULTS = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'patent_bench_results.json')
OUT = os.path.join(HERE, 'patent_bench_report')

_cache = {}


def forms(smi):
    """(exact, flat, compound) keys of one molecule, None where RDKit cannot read it"""
    if smi in _cache: return _cache[smi]
    m = Chem.MolFromSmiles(smi)
    if m is None: v = (None, None, None)
    else:
        ik = Chem.MolToInchiKey(m) or ''
        v = (Chem.MolToSmiles(m), Chem.MolToSmiles(m, isomericSmiles=False), ik[:14] or None)
    _cache[smi] = v
    return v


def heavy(smiles_list):
    n = 0
    for s in smiles_list:
        m = Chem.MolFromSmiles(s)
        if m: n += m.GetNumHeavyAtoms()
    return n


def tier(note):
    if re.match(r'(R\d{3}|SN1|SN2|E1|E2)\b', note): return 'workbook rule'
    if note.startswith('Predicted by the relevance model'): return 'relevance model (learned)'
    if re.match(r'L\d{6,7}\b', note): return 'recorded reaction'              # the exact look-up (off in a benchmark run, so never expected)
    if note.startswith('No reaction'): return 'no reaction (workbook)'
    if note.startswith('Unknown'): return 'unknown'
    return 'other'


def grade(case, res):
    if res.get('error'): return {'outcome': 'error', 'tier': 'error', 'why': res['error']}
    got = res.get('products') or []
    note = res.get('note') or ''
    t = tier(note)
    want = [x for x in case['product'].split('.') if x]

    def match(products):
        gf = [forms(x) for x in products]
        r = {}
        for k, name in enumerate(('exact', 'flat', 'compound')):
            have = {f[k] for f in gf if f[k]}
            r[name] = bool(want) and all(forms(w)[k] and forms(w)[k] in have for w in want)
        return r
    if not got:
        out = {'outcome': 'declined', 'tier': t}
        g = re.search(r"best guess was (\S+) \(template \S+\), but its estimated chance of being right is (\d+)%", res.get('reading') or '')
        if g: out['guess'] = {**match(g.group(1).split('.')), 'stated': int(g.group(2)) / 100}       # the model's answer the page held back
        return out
    out = {'outcome': 'answered', 'tier': t, **match(got)}
    m = re.search(r'(\d+)% estimated chance', note)
    if m: out['stated'] = int(m.group(1)) / 100
    return out


def pct(a, b): return f'{100 * a / b:.1f} %' if b else '-'


def summarise(items):
    n = len(items)
    ans = [g for _, g in items if g['outcome'] == 'answered']
    d = {'n': n, 'answered': len(ans), 'declined': sum(g['outcome'] == 'declined' for _, g in items), 'error': sum(g['outcome'] == 'error' for _, g in items)}
    for k in ('exact', 'flat', 'compound'):
        d[k] = sum(g[k] for g in ans)
    return d


def main():
    cases = {c['id']: c for c in json.load(open(CASES, encoding='utf-8'))['cases']}
    results = json.load(open(RESULTS, encoding='utf-8'))
    graded = [(cases[r['id']], grade(cases[r['id']], r)) for r in results if r['id'] in cases]
    for (c, g), r in zip(graded, [r for r in results if r['id'] in cases]):
        g['ms'] = r.get('ms'); g['note'] = (r.get('note') or '')[:300]; g['products'] = r.get('products') or []
    total = summarise(graded)
    by_tier = {t: summarise([x for x in graded if x[1]['tier'] == t]) for t in sorted({g['tier'] for _, g in graded})}

    def size_bin(c):
        h = heavy(c['reactants'])
        return '< 8' if h < 8 else '8-11' if h < 12 else '12-19' if h < 20 else '20-29' if h < 30 else '30 or more'
    by_size = collections.OrderedDict((b, summarise([x for x in graded if size_bin(x[0]) == b])) for b in ('< 8', '8-11', '12-19', '20-29', '30 or more'))
    by_reagent = {'reagent given': summarise([x for x in graded if x[0]['reagent']]), 'no reagent': summarise([x for x in graded if not x[0]['reagent']])}
    cal = []                                           # the model's stated chance against how often its first answer was right (exact)
    for lo, hi in ((0.5, 0.7), (0.7, 0.8), (0.8, 0.9), (0.9, 0.95), (0.95, 1.01)):
        xs = [g for _, g in graded if g['outcome'] == 'answered' and g.get('stated') is not None and lo <= g['stated'] < hi]
        if xs: cal.append({'stated': f'{lo:.2f}-{min(hi, 1):.2f}', 'n': len(xs), 'mean_stated': round(sum(g['stated'] for g in xs) / len(xs), 3),
                           'right_exact': round(sum(g['exact'] for g in xs) / len(xs), 3), 'right_compound': round(sum(g['compound'] for g in xs) / len(xs), 3)})
    held = [g['guess'] for _, g in graded if g.get('guess')]
    withheld = {'declined_with_a_guess': len(held), 'guess_right_exact': sum(x['exact'] for x in held), 'guess_right_compound': sum(x['compound'] for x in held),
                'by_stated': [{'stated': f'{lo:.2f}-{hi:.2f}', 'n': len(xs), 'right_exact': sum(x['exact'] for x in xs)}
                              for lo, hi in ((0, 0.3), (0.3, 0.5), (0.5, 0.7), (0.7, 0.85), (0.85, 1.01))
                              for xs in [[x for x in held if lo <= x['stated'] < hi]] if xs]}
    wrong = [{'id': c['id'], 'patent': c['patent'], 'code': c['code'], 'recorded': c['product'], 'page': '.'.join(g['products']), 'note': g['note'][:160], 'tier': g['tier']}
             for c, g in graded if g['outcome'] == 'answered' and not g['compound']]
    rep = {'about': json.load(open(CASES, encoding='utf-8'))['about'], 'total': total, 'by_tier': by_tier, 'by_size': by_size, 'by_reagent': by_reagent,
           'calibration': cal, 'withheld': withheld, 'median_ms': sorted(g['ms'] or 0 for _, g in graded)[len(graded) // 2] if graded else None,
           'wrong_examples': wrong[:40],
           'declined_reasons': collections.Counter(re.sub(r'\s+', ' ', g['note'])[:70] for _, g in graded if g['outcome'] == 'declined').most_common(12),
           'errors': collections.Counter(g.get('why', '') for _, g in graded if g['outcome'] == 'error').most_common(5)}
    json.dump(rep, open(OUT + '.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

    T = total
    L = ['# The Reactions page against 1,000 held-out US patent reactions', '', rep['about'], '',
         f"Each reaction's starting materials, reagent and recorded solvents were typed into the code bar and the arrow clicked, with the exact look-up of "
         f"recorded reactions turned off. Graded with RDKit (tools/patent_bench_grade.py). Median time per reaction: {rep['median_ms']} ms.", '',
         '## Overall', '', '| | reactions | share of all |', '|---|---|---|',
         f"| answered (showed products) | {T['answered']} | {pct(T['answered'], T['n'])} |",
         f"| declined (Unknown / No reaction) | {T['declined']} | {pct(T['declined'], T['n'])} |",
         f"| error | {T['error']} | {pct(T['error'], T['n'])} |", '',
         '| recorded product among the page\'s products | of all | of those answered |', '|---|---|---|']
    for k, label in (('exact', 'exact (with stereo)'), ('flat', 'ignoring stereo'), ('compound', 'same compound (salt / protonated form counts)')):
        L.append(f"| {label} | {T[k]} ({pct(T[k], T['n'])}) | {pct(T[k], T['answered'])} |")
    def table(title, groups):
        L.extend(['', f'## {title}', '', '| | reactions | answered | right, exact | right, same compound | right when answered (exact) |', '|---|---|---|---|---|---|'])
        for name, d in groups.items():
            L.append(f"| {name} | {d['n']} | {d['answered']} | {d['exact']} ({pct(d['exact'], d['n'])}) | {d['compound']} ({pct(d['compound'], d['n'])}) | {pct(d['exact'], d['answered'])} |")
    table('By the part of the page that answered', by_tier)
    table('By size of the starting materials (heavy atoms)', by_size)
    table('With and without a reagent', by_reagent)
    L.extend(['', '## Is the model\'s stated chance honest?', '', '| stated chance | answers | mean stated | right (exact) | right (same compound) |', '|---|---|---|---|---|'])
    for c in cal: L.append(f"| {c['stated']} | {c['n']} | {100 * c['mean_stated']:.0f} % | {100 * c['right_exact']:.0f} % | {100 * c['right_compound']:.0f} % |")
    W = withheld
    L.extend(['', '## The answers it held back', '',
              f"When the page declines with \"Unknown\", its reading still names the relevance model's best guess and its stated chance (below the floor the page "
              f"needs). {W['declined_with_a_guess']} declined reactions had one; {W['guess_right_exact']} of those guesses were right exactly "
              f"({pct(W['guess_right_exact'], W['declined_with_a_guess'])}), {W['guess_right_compound']} the same compound.", '',
              '| stated chance | held-back guesses | right (exact) |', '|---|---|---|'] +
             [f"| {b['stated']} | {b['n']} | {b['right_exact']} ({pct(b['right_exact'], b['n'])}) |" for b in W['by_stated']])
    L.extend(['', '## Why it declined (most common first lines)', ''] + [f'* {n} × {t}' for t, n in rep['declined_reasons']])
    L.extend(['', '## Wrong answers (first 40: the page\'s product is not even the same compound)', '', '| patent | starting materials > reagent | recorded | page | answered by |', '|---|---|---|---|---|'])
    for w in wrong[:40]: L.append(f"| {w['patent']} | `{w['code']}` | `{w['recorded']}` | `{w['page']}` | {w['tier']} |")
    open(OUT + '.md', 'w', encoding='utf-8').write('\n'.join(L) + '\n')
    print('\n'.join(L[:24]))
    print(f'wrote {OUT}.json and {OUT}.md')


if __name__ == '__main__':
    main()
