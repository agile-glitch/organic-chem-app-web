"""Which groups, present in the same starting materials, change how readily another group is consumed by the SAME kind of reagent?
(competing sites, quenching / incompatible groups, bystanders that hinder or help)   -- from the fragment count tables, no RDKit needed.

    python tools/rel_reactivity.py --counts --tag new      # the count tables it reads (about 17 minutes)
    python tools/interference_scan.py                      # about 4 minutes; writes interference_pairs.json in the temp folder and prints the strongest pairs

For a target fragment F and another fragment G (both in the starting materials), per reagent role signature (stratum):
   consumed(F) = count of F drops from starting materials to product.
   G state: absent / present and left alone (bystander) / present and consumed (competitor or co-reactant).
Mantel-Haenszel log odds ratios of "F consumed" for each G state against G absent, over the strata. Discovered on TRAIN patents, confirmed on the
held-out (validation + test) patents: the pair must have the same sign there and |z| >= 4.
"""
import sys, os, time, json, tempfile, collections
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, 'tools'))
import rel_common as C
import rel_reactivity as RV

T = tempfile.gettempdir()
t0 = time.time()
z = np.load(os.path.join(T, 'organic_chem_relfragcounts_new.npz'))
S, P, ids = z['S'], z['P'], [str(x) for x in z['fragments']]
nf = len(ids)
D = C.read_rows(); rows = D['rows']
roles = RV.reagent_roles()
cls = RV.classes_of(rows, roles)
sigs = [('none' if r[1] == '' else ('?' if r[1] is None else r[1])) for r in cls]
sig_index = {s: i for i, s in enumerate(sorted(set(sigs)))}
sig = np.array([sig_index[s] for s in sigs], np.int32)
split = np.array(C.split_codes(rows))
print(f'loaded {len(rows):,} reactions, {len(sig_index)} reagent role signatures ({time.time() - t0:.0f} s)', flush=True)

names = {}
import openpyxl
wb = openpyxl.load_workbook(os.path.join(REPO, 'molecule_data.xlsx'), read_only=True, data_only=True)
for r in list(wb['Fragments'].iter_rows(values_only=True))[1:]:
    if r[0]: names[str(r[0]).strip()] = str(r[1])

present = S >= 1
consumed = S > P
NS = len(sig_index)


def mh(n, c, g):
    """Mantel-Haenszel log OR (state g vs state 0) over strata, and its standard error (Robins-Breslow-Greenland); n, c: (NS, 3) arrays"""
    a = c[:, g]; b = n[:, g] - c[:, g]; c0 = c[:, 0]; d = n[:, 0] - c[:, 0]
    N = n[:, g] + n[:, 0]
    ok = (n[:, g] > 0) & (n[:, 0] > 0)
    if not ok.any(): return None
    a, b, c0, d, N = a[ok], b[ok], c0[ok], d[ok], N[ok]
    R = a * d / N; Sx = b * c0 / N
    if R.sum() <= 0 or Sx.sum() <= 0: return None
    lor = np.log(R.sum() / Sx.sum())
    P_ = (a + d) / N; Q_ = (b + c0) / N
    var = (np.sum(P_ * R) / (2 * R.sum() ** 2)) + (np.sum(P_ * Sx + Q_ * R) / (2 * R.sum() * Sx.sum())) + (np.sum(Q_ * Sx) / (2 * Sx.sum() ** 2))
    return lor, float(np.sqrt(var)), int(n[ok][:, g].sum())


def stats_for(mask_rows, F):
    """per G: {state: (lor, se, n)} for the reactions in mask_rows where F is present"""
    idx = np.flatnonzero(mask_rows & present[:, F])
    if len(idx) < 2000: return None
    y = consumed[idx, F].astype(np.float64)
    sg = sig[idx] * 3
    gstate = np.zeros((len(idx), nf), np.uint8)
    pr = present[idx]; co = consumed[idx]
    gstate[pr] = 1
    gstate[pr & co] = 2
    out = {}
    for G in range(nf):
        if G == F: continue
        comb = sg + gstate[:, G]
        n = np.bincount(comb, minlength=NS * 3).reshape(NS, 3).astype(np.float64)
        c = np.bincount(comb, weights=y, minlength=NS * 3).reshape(NS, 3)
        res = {}
        for g in (1, 2):
            if n[:, g].sum() >= 300 and n[:, 0].sum() >= 300:
                m = mh(n, c, g)
                if m: res[g] = m
        if res: out[G] = res
    return out


train = split == C.TRAIN
held = split != C.TRAIN
found = []
for F in range(nf):
    a = stats_for(train, F)
    if not a: continue
    b = stats_for(held, F) or {}
    for G, res in a.items():
        # how much of the time is G a part of what F is (a sub-fragment relation rather than a bystander)?
        both = (present[:, F] & present[:, G]).sum(); g_only = present[:, G].sum()
        implies = both / max(1, g_only)
        for g, (lor, se, n) in res.items():
            hb = b.get(G, {}).get(g)
            if not hb: continue
            lor2, se2, n2 = hb
            zt, zh = lor / max(se, 1e-9), lor2 / max(se2, 1e-9)
            if abs(zt) >= 6 and abs(zh) >= 4 and np.sign(lor) == np.sign(lor2):
                found.append({'F': ids[F], 'Fname': names.get(ids[F], ''), 'G': ids[G], 'Gname': names.get(ids[G], ''), 'state': 'bystander' if g == 1 else 'consumed',
                              'logOR_train': round(float(lor), 3), 'logOR_held': round(float(lor2), 3), 'z_held': round(float(zh), 1), 'n_pairs': n, 'G_implies_F': round(float(implies), 2)})
    if F % 20 == 19: print(f'  {F + 1}/{nf} targets ({time.time() - t0:.0f} s), {len(found)} pairs so far', flush=True)

json.dump(found, open(os.path.join(T, 'interference_pairs.json'), 'w'), indent=1)
print(f'{len(found)} confirmed pairs ({time.time() - t0:.0f} s)')


def show(title, lst, k=25):
    print('\n==', title)
    for r in lst[:k]:
        print(f"  {r['F']} {r['Fname'][:30]:30s} | {r['G']} {r['Gname'][:30]:30s} {r['state']:9s} OR {np.exp(r['logOR_held']):5.2f} (train {np.exp(r['logOR_train']):5.2f}) z {r['z_held']:6.1f} n {r['n_pairs']:7,d}  G=>F {r['G_implies_F']}")

ok = [r for r in found if r['G_implies_F'] < 0.9]            # drop sub-fragment relations (G always comes with F)
by = sorted(ok, key=lambda r: -abs(r['z_held']))
show('G present and LEFT ALONE makes F LESS likely to react (hinders)', [r for r in by if r['state'] == 'bystander' and r['logOR_held'] < 0])
show('G present and LEFT ALONE makes F MORE likely to react (helps)', [r for r in by if r['state'] == 'bystander' and r['logOR_held'] > 0])
show('G present and CONSUMED, F less likely to react (a competitor)', [r for r in by if r['state'] == 'consumed' and r['logOR_held'] < 0])
show('G present and CONSUMED, F more likely to react (co-reacting)', [r for r in by if r['state'] == 'consumed' and r['logOR_held'] > 0], 12)
