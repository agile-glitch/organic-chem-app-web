"""Fragment Atlas, step 1: read atom-mapped recorded reactions and write down, for every functional-group instance in every
reacting molecule, what the reaction did to it.

    python tools/atlas_mine.py --src <folder with train.txt valid.txt test.txt> --out <instances.npz> [--n 20000] [--workers 4]

Input: the USPTO-MIT set of Jin, Coley et al. (NIPS 2017; ~480,000 US-patent reactions text-mined by D. M. Lowe, CC0, every atom
mapped), one reaction per line:  mapped reactants>>mapped product  <space>  changed-bond pairs.
From it the reaction's own atom maps say exactly which atoms changed, so no template and no guess is needed:

  * a reactant molecule none of whose atoms is in the product is a SPECTATOR (reagent, catalyst, base, solvent);
  * an atom CHANGED if its element charge, H count, bonds (to which mapped neighbours, of which kind) differ between the reactant
    and the product side, or if it is not in the product at all;
  * a molecule that has atoms in the product and changed atoms is a PARTICIPANT; its fragment instances (the 205 named groups of
    the Fragments sheet, every match) are what is counted. An instance is HIT (y = 1) if any atom of its match changed.
  * the reagent side of an instance = the spectators and the other participants that the Reagents sheet knows (role, pKa).

One row per instance is written (the group, the four environment hashes, the continuous descriptors, up to three reagent parts,
the reagent class, hit or not, and which scaffold split the reaction falls in), plus the part vocabulary. tools/atlas_eval.py
fits and measures on it; tools/atlas_build.py turns the counts into the tables the page loads.

Split (scaffold, so close relatives never straddle the line): hash of the Bemis-Murcko scaffold of the largest participant;
0 = TEST (5 %), 1 = VALIDATION (5 %), the rest TRAIN. Reaction-level splits leak badly here: one patent records dozens of
near-identical reactions.
"""
import argparse, os, sys, time
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atlas_core as A

_W = {}


def _init():
    from rdkit import RDLogger
    RDLogger.DisableLog('rdApp.*')
    _W['frags'] = A.fragments()
    _W['fidx'] = {f[0]: k for k, f in enumerate(_W['frags'])}
    _W['roles'] = A.reagent_roles()


def split_of(key):
    h = A.fnv32(key) % 20
    return 0 if h == 0 else 1 if h == 1 else 2          # 0 test, 1 validation, 2 train


def _bond_cls(b):
    return 'A' if b.GetIsAromatic() else str(int(b.GetBondTypeAsDouble()))


def _sigs(mols):
    sig = {}
    for m in mols:
        for a in m.GetAtoms():
            mp = a.GetAtomMapNum()
            if not mp: continue
            nb = tuple(sorted((n.GetAtomMapNum(), _bond_cls(m.GetBondBetweenAtoms(a.GetIdx(), n.GetIdx()))) for n in a.GetNeighbors() if n.GetAtomMapNum()))
            sig[mp] = (a.GetAtomicNum(), a.GetFormalCharge(), a.GetTotalNumHs(), nb)
    return sig


def process(line):
    """-> None, or (split, class string, [instance tuples]) for one reaction line"""
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    try:
        rxn = line.split(' ')[0]
        lhs, rhs = rxn.split('>>')
        rmols = [Chem.MolFromSmiles(s) for s in lhs.split('.')]
        pmols = [Chem.MolFromSmiles(s) for s in rhs.split('.')]
    except Exception:
        return None
    if any(m is None for m in rmols) or any(m is None for m in pmols): return None
    rs, ps = _sigs(rmols), _sigs(pmols)
    if not rs or not ps: return None
    changed = {mp for mp, s in rs.items() if mp not in ps or ps[mp] != s}
    # molecules: spectator / participant
    info = []
    for m in rmols:
        maps = [a.GetAtomMapNum() for a in m.GetAtoms() if a.GetAtomMapNum()]
        inprod = [mp for mp in maps if mp in ps]
        ch = [mp for mp in maps if mp in changed]
        info.append((m, maps, bool(inprod), bool(ch)))
    parts = []                                           # (smiles, role, pKa, owner index or -1 for spectators)
    roles = _W['roles']
    for k, (m, maps, inprod, ch) in enumerate(info):
        mm = Chem.Mol(m)
        for a in mm.GetAtoms(): a.SetAtomMapNum(0)
        smi = Chem.MolToSmiles(mm)
        r = roles.get(smi)
        info[k] = (m, maps, inprod, ch, mm, smi)
        if r is not None and r[0] not in A.SKIP_ROLES:
            parts.append((smi, r[0], r[1], -1 if not inprod else k))
    participants = [k for k, t in enumerate(info) if t[2] and t[3]]
    if not participants: return None
    big = max(participants, key=lambda k: info[k][4].GetNumHeavyAtoms())
    try:
        sc = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(info[big][4]))
    except Exception:
        sc = ''
    split = split_of(sc if sc else info[big][5])
    # molecule-disjoint split: a reaction is TEST if ANY participant molecule hashes into the test bucket, else VALIDATION if any hashes into
    # that one, else TRAIN. A held-out molecule then never occurs in a training reaction at all (as substrate or as partner).
    mol_h = {k: A.fnv32(info[k][5]) % 20 for k in participants}
    split_mol = 0 if any(h == 0 for h in mol_h.values()) else 1 if any(h == 1 for h in mol_h.values()) else 2
    out = []
    frags, fidx = _W['frags'], _W['fidx']
    for k in participants:
        m, maps, _ip, _ch, mm, smi = info[k]
        # a molecule already seen would only repeat itself: salts and metals have no fragments anyway
        try:
            Chem.SanitizeMol(mm)
        except Exception:
            continue
        mapof = {a.GetIdx(): a.GetAtomMapNum() for a in m.GetAtoms()}
        wl = A.wl_hashes(mm)
        q = A.gasteiger(mm)
        other = [p for p in parts if p[3] != k]
        cls = A.reagent_class([p[1] for p in other])
        chosen = sorted(other, key=lambda p: (p[1] in ('solvent',), p[0]))[:3]
        for fid, atoms, x in A.instances(mm, frags):
            y = int(any(mapof.get(i, 0) in changed for i in atoms))
            out.append((fidx[fid], x, wl[0][x], wl[1][x], wl[2][x], wl[3][x], y, A.steric_score(mm, x), q[x] if q[x] is not None else np.nan,
                        A.smallest_ring(mm, x), A.angle_strain(mm, x), A.hybridization(mm.GetAtomWithIdx(x)),
                        cls, tuple((p[0], p[1], p[2]) for p in chosen), len(atoms), A.fnv32(smi), split_mol))
    if not out: return None
    return split, out


def work(args):
    lines, = args
    res = []
    for ln in lines:
        try:
            r = process(ln)
        except Exception:
            r = None
        res.append(r)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--n', type=int, default=0, help='only the first N lines of each file (a quick test)')
    ap.add_argument('--workers', type=int, default=4)
    a = ap.parse_args()
    lines = []
    for f in ('train.txt', 'valid.txt', 'test.txt'):
        with open(os.path.join(a.src, f), encoding='utf-8') as fh:
            for i, ln in enumerate(fh):
                if a.n and i >= a.n: break
                lines.append(ln.strip())
    print(f'{len(lines):,} reactions', flush=True)
    t0 = time.time()
    chunk = 400
    jobs = [(lines[i:i + chunk],) for i in range(0, len(lines), chunk)]
    cols = {k: [] for k in ('rx', 'frag', 'x', 'e0', 'e1', 'e2', 'e3', 'y', 'steric', 'charge', 'ring', 'strain', 'hyb', 'cls', 'p1', 'p2', 'p3', 'split', 'natoms', 'mol', 'split_mol')}
    cls_voc, part_voc, part_info = {}, {}, {}
    rid = 0
    used = 0
    with Pool(a.workers, initializer=_init) as pool:
        for ci, res in enumerate(pool.imap(work, jobs, chunksize=1)):
            for r in res:
                rid += 1
                if r is None: continue
                used += 1
                split, inst = r
                for (fi, x, e0, e1, e2, e3, y, st, q, ring, strain, hyb, cls, pts, na, mh, sm) in inst:
                    cols['rx'].append(rid); cols['frag'].append(fi); cols['x'].append(x)
                    cols['e0'].append(e0); cols['e1'].append(e1); cols['e2'].append(e2); cols['e3'].append(e3)
                    cols['y'].append(y); cols['steric'].append(st); cols['charge'].append(q); cols['ring'].append(ring)
                    cols['strain'].append(strain); cols['hyb'].append(hyb); cols['natoms'].append(na)
                    cols['cls'].append(cls_voc.setdefault(cls, len(cls_voc)))
                    ids = [0, 0, 0]
                    for j, (s, role, pka) in enumerate(pts):
                        if s not in part_voc: part_voc[s] = len(part_voc) + 1; part_info[s] = (role, pka)
                        ids[j] = part_voc[s]
                    cols['p1'].append(ids[0]); cols['p2'].append(ids[1]); cols['p3'].append(ids[2])
                    cols['split'].append(split); cols['mol'].append(mh); cols['split_mol'].append(sm)
            if ci % 50 == 0:
                print(f'  {rid:,} reactions read, {used:,} usable, {len(cols["y"]):,} instances, {time.time() - t0:.0f}s', flush=True)
    dt = {'rx': np.uint32, 'frag': np.uint16, 'x': np.uint16, 'e0': np.uint32, 'e1': np.uint32, 'e2': np.uint32, 'e3': np.uint32, 'y': np.uint8,
          'steric': np.uint8, 'charge': np.float32, 'ring': np.uint8, 'strain': np.float32, 'hyb': np.uint8, 'cls': np.uint16,
          'p1': np.uint32, 'p2': np.uint32, 'p3': np.uint32, 'split': np.uint8, 'natoms': np.uint8, 'mol': np.uint32, 'split_mol': np.uint8}
    arrays = {k: np.asarray(v, dtype=dt[k]) for k, v in cols.items()}
    inv_parts = sorted(part_voc, key=part_voc.get)
    inv_cls = sorted(cls_voc, key=cls_voc.get)
    np.savez_compressed(a.out, **arrays, parts=np.array(inv_parts, dtype=object), part_role=np.array([part_info[s][0] for s in inv_parts], dtype=object),
                        part_pka=np.array([np.nan if part_info[s][1] is None else part_info[s][1] for s in inv_parts], dtype=np.float32),
                        classes=np.array(inv_cls, dtype=object), frag_ids=np.array([f[0] for f in A.fragments()], dtype=object), allow_pickle=True)
    print(f'done: {rid:,} reactions, {used:,} usable, {len(arrays["y"]):,} instances, {len(inv_parts)} reagent parts, {len(inv_cls)} classes, '
          f'{time.time() - t0:.0f}s -> {a.out}', flush=True)


if __name__ == '__main__':
    main()
