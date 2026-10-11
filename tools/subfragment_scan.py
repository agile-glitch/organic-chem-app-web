"""Which fragments lump together groups that react differently?  Data-driven sub-fragment discovery.

    python tools/subfragment_scan.py --n 200000 --held 60000        # about 20 minutes; writes tools/subfragment_candidates.json

A fragment of the Fragments sheet matches many different surroundings: a ketone can be an aryl ketone, an enone, a ring ketone, a methyl ketone ...
For every fragment F and every way its match can be told apart,
    own atom:   the matched atom is aromatic / aliphatic / in a ring (of 3..7 atoms) / has 0-3 hydrogens
    neighbour:  an atom outside the match, attached to the matched atom, is an aryl C, CH3, CH2, CH, quaternary C, vinyl C, carbonyl C, alkynyl C,
                O, N, charged N, S, a halogen, Si/B/P or a metal
the scan asks: when the starting materials have F, does the recorded reaction CONSUME F (its count drops from the starting materials to the product)
differently with this feature than without it, given the same kind of reagent (the reagent's role signature, Reagents sheet)? A feature that changes
how F reacts is a sub-fragment worth naming: its SMARTS (the parent with the atom constraint or the extra neighbour atom added) is written out.

Honest measurement: the rates are learned on TRAIN patents (--n reactions) and the gain in log-likelihood (nats) is scored on reactions of the
VALIDATION patents (--held) that the rates never saw. Candidates are ranked by that held-out gain.
"""
import argparse, collections, json, os, sys, tempfile, time
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rel_common as C
import rel_reactivity as RV

TMP = tempfile.gettempdir()
NSIG = 160                                    # reagent role signatures kept (the most frequent), the rest share one bucket
OWN = [('aromatic', 'a'), ('aliphatic', 'A'), ('in a ring', 'R'), ('not in a ring', '!R'), ('in a 3-ring', 'r3'), ('in a 4-ring', 'r4'), ('in a 5-ring', 'r5'),
       ('in a 6-ring', 'r6'), ('in a 7-ring', 'r7'), ('no H', 'H0'), ('1 H', 'H1'), ('2 H', 'H2'), ('3 H', 'H3')]
NBR = [('aryl C', 'c'), ('heteroaryl atom', '[n,o,s;+0]'), ('CH3', '[CX4;H3]'), ('CH2', '[CX4;H2]'), ('CH (branched)', '[CX4;H1]'), ('quaternary C', '[CX4;H0]'),
       ('vinyl C', '[CX3;!$(C=[O,N,S])]'), ('carbonyl / imine C', '[CX3;$(C=[O,N,S])]'), ('alkynyl / nitrile C', '[CX2]'), ('O', '[OX2]'),
       ('N', '[NX3,nX3]'), ('charged N', '[N+,n+]'), ('S', '[#16]'), ('F', 'F'), ('Cl', 'Cl'), ('Br', 'Br'), ('I', 'I'), ('Si / B / P', '[Si,B,P]'),
       ('metal', '[Li,Na,K,Mg,Zn,Cu,Sn,Pd,Al,Ca]')]
_W = {}


def own_labels(mol, atom):
    out = ['aromatic' if atom.GetIsAromatic() else 'aliphatic']
    ri = mol.GetRingInfo()
    if atom.IsInRing():
        out.append('in a ring')
        n = ri.MinAtomRingSize(atom.GetIdx())
        if 3 <= n <= 7: out.append(f'in a {n}-ring')
    else:
        out.append('not in a ring')
    h = atom.GetTotalNumHs()
    if h <= 3: out.append('no H' if h == 0 else f'{h} H')
    return out


def nbr_labels(mol, nb):
    s = nb.GetSymbol()
    if s == 'C':
        if nb.GetIsAromatic(): return ['aryl C']
        dbl = [b.GetOtherAtom(nb).GetSymbol() for b in nb.GetBonds() if b.GetBondTypeAsDouble() == 2.0]
        if any(x in ('O', 'N', 'S') for x in dbl): return ['carbonyl / imine C']
        if 'C' in dbl: return ['vinyl C']
        if any(b.GetBondTypeAsDouble() == 3.0 for b in nb.GetBonds()): return ['alkynyl / nitrile C']
        h = nb.GetTotalNumHs()
        return ['CH3' if h == 3 else 'CH2' if h == 2 else 'CH (branched)' if h == 1 else 'quaternary C']
    if s in ('N', 'O', 'S') and nb.GetIsAromatic():
        return ['charged N' if (s == 'N' and nb.GetFormalCharge() > 0) else 'heteroaryl atom']
    if s == 'O': return ['O']
    if s == 'N': return ['charged N' if nb.GetFormalCharge() > 0 else 'N']
    if s == 'S': return ['S']
    if s in ('F', 'Cl', 'Br', 'I'): return [s]
    if s in ('Si', 'B', 'P'): return ['Si / B / P']
    if s in ('Li', 'Na', 'K', 'Mg', 'Zn', 'Cu', 'Sn', 'Pd', 'Al', 'Ca'): return ['metal']
    return []


def vocabulary(frags):
    """[(fragment index, position class, kind 'own'/'nbr', label)] and, per fragment, the position classes: [[positions of class 0], ...]"""
    from rdkit import Chem
    vocab, classes = [], []
    for fi, (fid, name, smarts) in enumerate(frags):
        q = Chem.MolFromSmarts(smarts)
        if q is None or q.GetNumAtoms() == 0: classes.append([]); continue
        q.UpdatePropertyCache(False); Chem.FastFindRings(q)
        ranks = list(Chem.CanonicalRankAtoms(q, breakTies=False))
        groups = collections.OrderedDict()
        for k, r in enumerate(ranks): groups.setdefault(r, []).append(k)
        classes.append(list(groups.values()))
        for ci in range(len(groups)):
            for lab, _ in OWN: vocab.append((fi, ci, 'own', lab))
            for lab, _ in NBR: vocab.append((fi, ci, 'nbr', lab))
    return vocab, classes


def _init(frag_smarts, classes, vocab):
    from rdkit import RDLogger, Chem
    RDLogger.DisableLog('rdApp.*')
    _W['q'] = [Chem.MolFromSmarts(s) for s in frag_smarts]
    _W['classes'] = classes
    _W['tok'] = {v: i for i, v in enumerate(vocab)}
    _W['nv'] = len(vocab)


def work(task):
    """task = (list of (row index, substrate smiles list, set of fragment indices present, {fragment: consumed}, signature bucket)) -> counts"""
    from rdkit import Chem
    items = task
    nv = _W['nv']
    idx1, y1 = [], []                                # token x signature cells of 'F present with this feature' and whether F was consumed
    idxF, yF = [], []                                # fragment x signature cells of 'F present'
    nF = len(_W['q'])
    for ri, subs, frs, cons, sb in items:
        mols = [Chem.MolFromSmiles(s) for s in subs]
        for f in frs:
            q = _W['q'][f]
            if q is None: continue
            toks = set()
            for m in mols:
                if m is None: continue
                for match in m.GetSubstructMatches(q, uniquify=True, maxMatches=30):
                    ms = set(match)
                    for ci, positions in enumerate(_W['classes'][f]):
                        for k in positions:
                            at = m.GetAtomWithIdx(match[k])
                            for lab in own_labels(m, at): toks.add(_W['tok'].get((f, ci, 'own', lab), -1))
                            for nb in at.GetNeighbors():
                                if nb.GetIdx() in ms: continue
                                for lab in nbr_labels(m, nb): toks.add(_W['tok'].get((f, ci, 'nbr', lab), -1))
            y = 1 if cons[f] else 0
            idxF.append(f * NSIG1 + sb); yF.append(y)
            for t in toks:
                if t >= 0: idx1.append(t * NSIG1 + sb); y1.append(y)
    def sparse(idx, y, size):
        if not idx: return np.zeros(0, np.int64), np.zeros(0, np.int32), np.zeros(0, np.int32)
        idx = np.array(idx, np.int64)
        n = np.bincount(idx, minlength=size); c = np.bincount(idx, weights=np.array(y, np.float64), minlength=size)
        nz = np.flatnonzero(n)
        return nz, n[nz].astype(np.int32), c[nz].astype(np.int32)
    return sparse(idx1, y1, nv * NSIG1), sparse(idxF, yF, nF * NSIG1)


NSIG1 = NSIG + 1


def collect(rows_idx, rows, S, P, sigb, frags, vocab, classes, workers, label):
    t0 = time.time()
    items = []
    for i in rows_idx:
        present = np.flatnonzero(S[i] >= 1)
        if not len(present): continue
        items.append((int(i), rows[i][1], set(int(f) for f in present), {int(f): bool(S[i, f] > P[i, f]) for f in present}, int(sigb[i])))
    chunks = [items[a:a + 4000] for a in range(0, len(items), 4000)]
    nv, nF = len(vocab), len(frags)
    N1 = np.zeros(nv * NSIG1, np.int64); C1 = np.zeros(nv * NSIG1, np.int64); NF = np.zeros(nF * NSIG1, np.int64); CF = np.zeros(nF * NSIG1, np.int64)
    with Pool(workers, initializer=_init, initargs=([f[2] for f in frags], classes, vocab)) as pool:
        for k, ((i1, n1, c1), (iF, nf, cf)) in enumerate(pool.imap_unordered(work, chunks), 1):
            N1[i1] += n1; C1[i1] += c1; NF[iF] += nf; CF[iF] += cf
            if k % 5 == 0: print(f'  {label}: {min(len(items), k * 4000):,} of {len(items):,} reactions ({time.time() - t0:.0f} s)', flush=True)
    return N1.reshape(nv, NSIG1), C1.reshape(nv, NSIG1), NF.reshape(nF, NSIG1), CF.reshape(nF, NSIG1)


def ll(c, n, p):
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return c * np.log(p) + (n - c) * np.log(1 - p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=200000); ap.add_argument('--held', type=int, default=60000); ap.add_argument('--seed', type=int, default=13)
    ap.add_argument('--workers', type=int, default=0); ap.add_argument('--out', default=os.path.join(HERE, 'subfragment_candidates.json'))
    ap.add_argument('--min-support', type=int, default=300)
    a = ap.parse_args()
    t0 = time.time()
    D = C.read_rows(); rows = D['rows']
    z = np.load(os.path.join(TMP, 'organic_chem_relfragcounts_new.npz'))
    S, P, ids = z['S'], z['P'], [str(x) for x in z['fragments']]
    frs = RV.fragment_list()
    assert [f[0] for f in frs] == ids, 'the fragment count table was made with a different fragment list: python tools/rel_reactivity.py --counts --tag new'
    roles = RV.reagent_roles()
    cls = RV.classes_of(rows, roles)
    sigs = np.array([('none' if r[1] == '' else ('?' if r[1] is None else r[1])) for r in cls], dtype=object)
    split = np.array(C.split_codes(rows))
    freq = collections.Counter(sigs[split == C.TRAIN])
    keep = {s: i for i, (s, _) in enumerate(freq.most_common(NSIG))}
    sigb = np.array([keep.get(s, NSIG) for s in sigs], np.int32)
    rng = np.random.default_rng(a.seed)
    tr = np.sort(rng.permutation(np.flatnonzero(split == C.TRAIN))[:a.n]); va = np.sort(rng.permutation(np.flatnonzero(split == C.VAL))[:a.held])
    vocab, classes = vocabulary(frs)
    print(f'{len(frs)} fragments, {len(vocab):,} features; learning rates on {len(tr):,} TRAIN reactions, scoring on {len(va):,} VALIDATION ones ({time.time() - t0:.0f} s)', flush=True)
    workers = a.workers or max(1, (os.cpu_count() or 2) - 1)
    A = collect(tr, rows, S, P, sigb, frs, vocab, classes, workers, 'train')
    B = collect(va, rows, S, P, sigb, frs, vocab, classes, workers, 'held-out')
    N1a, C1a, NFa, CFa = [x.astype(np.float64) for x in A]
    N1b, C1b, NFb, CFb = [x.astype(np.float64) for x in B]
    fi = np.array([v[0] for v in vocab])
    N0a, C0a = NFa[fi] - N1a, CFa[fi] - C1a
    N0b, C0b = NFb[fi] - N1b, CFb[fi] - C1b
    # rates per stratum from the TRAIN counts, shrunk toward the stratum rate of the fragment (and that toward its overall rate)
    allrate = (CFa.sum(1) + 1) / (NFa.sum(1) + 2)
    ps = (CFa + 4 * allrate[:, None]) / (NFa + 4)
    PS = ps[fi]
    p1 = (C1a + 6 * PS) / (N1a + 6); p0 = (C0a + 6 * PS) / (N0a + 6)
    gain_train = (ll(C1a, N1a, p1) + ll(C0a, N0a, p0) - ll(C1a, N1a, PS) - ll(C0a, N0a, PS)).sum(1)
    gain_held = (ll(C1b, N1b, p1) + ll(C0b, N0b, p0) - ll(C1b, N1b, PS) - ll(C0b, N0b, PS)).sum(1)
    n1, n0 = N1a.sum(1), N0a.sum(1)
    c1, c0 = C1a.sum(1), C0a.sum(1)
    out = []
    for t, v in enumerate(vocab):
        if n1[t] < a.min_support or n0[t] < a.min_support: continue
        # the plain effect: consumption rate with / without the feature, averaged over reagent strata with weights of the fragment's own strata
        w = N1a[t] + N0a[t]
        eff = float(((p1[t] - p0[t]) * w).sum() / max(w.sum(), 1))
        out.append({'fragment': ids[v[0]], 'name': frs[v[0]][1], 'position_class': v[1], 'positions': classes[v[0]][v[1]], 'kind': v[2], 'feature': v[3],
                    'support': int(n1[t]), 'complement': int(n0[t]), 'rate_with': round(float(c1[t] / n1[t]), 3), 'rate_without': round(float(c0[t] / n0[t]), 3),
                    'effect_within_reagent': round(eff, 3), 'gain_train_nats': round(float(gain_train[t]), 1), 'gain_held_nats': round(float(gain_held[t]), 1),
                    'smarts': frs[v[0]][2]})
    out.sort(key=lambda r: -r['gain_held_nats'])
    json.dump({'when': time.strftime('%Y-%m-%d'), 'train_reactions': len(tr), 'held_out_reactions': len(va), 'candidates': out}, open(a.out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print(f'{len(out):,} candidates with support; best held-out gains:')
    for r in out[:40]:
        print(f"  {r['gain_held_nats']:8.1f} (train {r['gain_train_nats']:8.1f})  {r['fragment']} {r['name'][:34]:34s} pos{r['positions']} {r['kind']:3s} {r['feature']:22s} with {r['support']:7,d} rate {r['rate_with']:.2f} vs {r['rate_without']:.2f}")
    print(f'({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
