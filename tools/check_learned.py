"""Re-check the learned reactions in learned_reactions.xlsx, independently of how they were made.

    python tools/check_learned.py              # a spread sample of 20,000 rows
    python tools/check_learned.py --all        # every row (about as long as learning them)

For each row:
  1. the code's four parts agree with the Starting materials / Reagent / Conditions / Product columns;
  2. the template (Templates sheet), run by RDKit on the starting materials plus the carbon-free reagent pieces,
     exactly as the page runs it, gives the product;
  3. the fragments consumed / formed, recounted with molecule_data.xlsx, are the ones written in the row;
  4. the look-up keys: InChIKey of every distinct molecule (starting materials + reagent parts) and of the reagent's
     distinct parts.
Prints how many rows pass each check and lists the first failures.
"""
import argparse, os, sys, time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SRC = os.path.join(os.path.dirname(HERE), 'learned_reactions.xlsx')

_S = {}


def _init(smarts, solv_frags):
    import learn_reactions as L
    L._init(*solv_frags, L.reagent_roles())
    _S['smarts'] = smarts


def check(row):
    import learn_reactions as L
    lid, code, subs, reagent, conds, prod, t, consumed, formed, skey, rkey = row
    bad = []
    if code != '~'.join([subs, reagent, conds, prod]): bad.append('code differs from its columns')
    sl = subs.split('.')
    extras = L.pool_extras(reagent, L._W['roles'])
    outs, _ = L.run_template(_S['smarts'][t], sl, extras)
    if not any(prod in o for o in outs): bad.append('template does not give the product')
    a0, b0 = L.frag_counts(sl), L.frag_counts([prod])
    if ' + '.join(sorted(f for f in a0 if b0[f] < a0[f])) != consumed: bad.append('fragments consumed differ')
    if ' + '.join(sorted(f for f in b0 if b0[f] > a0[f])) != formed: bad.append('fragments formed differ')
    if L.lookup_key(sl, reagent) != skey or L.reagent_key(reagent) != rkey: bad.append('look-up keys differ')
    return lid, bad


def main():
    import openpyxl, learn_reactions as L
    ap = argparse.ArgumentParser()
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--sample', type=int, default=20000)
    args = ap.parse_args()
    t0 = time.time()
    wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
    trows = [r for r in list(wb['Templates'].iter_rows(values_only=True))[1:] if r and r[0]]
    tindex = {r[0]: k for k, r in enumerate(trows)}
    rows = []
    for name in wb.sheetnames:
        if not name.startswith('Learned'): continue
        head = None
        for r in wb[name].iter_rows(values_only=True):
            if head is None: head = {h: k for k, h in enumerate(r)}; continue
            if not r or not r[0]: continue
            g = lambda col: r[head[col]] or ''
            rows.append((g('Learned ID'), g('Reaction code (starting materials~reagent~conditions~product)'), g('Starting materials'), g('Reagent'),
                         g('Conditions (recorded solvents)'), g('Product'), tindex[r[head['Template ID']]], g('Fragments consumed'), g('Fragments formed'),
                         g('Look-up key (InChIKey of starting materials + reagent parts)'), g('Reagent key (InChIKey of its distinct parts)')))
    pick = rows if args.all else rows[::max(1, len(rows) // args.sample)][:args.sample]
    print(f'{len(rows):,} learned reactions; checking {len(pick):,} ({time.time() - t0:.0f} s)', flush=True)
    fails, kinds, n = [], {}, 0
    with Pool(max(1, (os.cpu_count() or 2) - 1), initializer=_init, initargs=([r[1] for r in trows], (L.solvents(), L.fragments()))) as pool:
        for lid, bad in pool.imap_unordered(check, pick, chunksize=200):
            n += 1
            for b in bad: kinds[b] = kinds.get(b, 0) + 1
            if bad: fails.append((lid, bad))
            if n % 50000 == 0: print(f'  {n:,} checked ({time.time() - t0:.0f} s)', flush=True)
    print(f'\n{n - len(fails):,} of {n:,} learned reactions pass every check ({time.time() - t0:.0f} s)')
    for k, v in kinds.items(): print(f'  {k}: {v:,}')
    for lid, bad in sorted(fails)[:20]: print(' ', lid, '; '.join(bad))


if __name__ == '__main__':
    main()
