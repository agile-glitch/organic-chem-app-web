"""What the recorded reactions say about each advisory of the Interferences sheet: does the group really react when that reagent is used?

    python tools/int_evidence.py --out tools/interference_evidence.json

For an advisory (a reagent side and a list of groups, tools/interference_defs.py ADVISORY) the recorded reactions whose reagent matches it (a role of the
Reagents sheet, a fragment found in a reagent part, or a reagent part as typed) are split into those that have one of the
groups in their starting materials and the rest. For each group the share of those reactions in which it was consumed (fragment count dropped,
tools/rel_reactivity.py) is compared with the share in reactions with other reagents:
    'also reacts' / 'quenches' advisories: the group is consumed much more often with the reagent (LiAlH4 takes the ester, the Grignard the O–H)
    'spared' advisories: much less often (NaBH4 leaves the ester)
Descriptive figures over every recorded reaction (nothing is predicted from them); an alternative of the reagent side that names a solvent type is not counted here (the recorded conditions are not in the feature table).
"""
import argparse, json, os, sys, tempfile, time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import interference_defs as D
import rel_features as F
import rel_reactivity as RV

TMP = tempfile.gettempdir()


def unpack_columns(mm, a, w, chunk=200000):
    """bits a .. a+w-1 of every row of the packed feature table, as a bool matrix (N x w)"""
    b0, b1 = a // 8, (a + w + 7) // 8
    out = np.zeros((mm.shape[0], w), bool)
    for s in range(0, mm.shape[0], chunk):
        bits = np.unpackbits(mm[s:s + chunk, b0:b1], axis=1)
        out[s:s + chunk] = bits[:, a - b0 * 8:a - b0 * 8 + w].astype(bool)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True); ap.add_argument('--counts', default='new'); ap.add_argument('--feat', default='v2')
    a = ap.parse_args()
    t0 = time.time()
    path, meta_path = F.table_paths(a.feat)
    meta = json.load(open(meta_path)); L = meta['layout']
    mm = np.load(path, mmap_mode='r')
    roles = meta['roles']; frs = meta['fragments']
    R = unpack_columns(mm, L['roles'][0], L['roles'][1])
    RF = unpack_columns(mm, L['frag_reag'][0], L['frag_reag'][1])
    import int_context as IC
    from scipy import sparse
    meta2 = IC.load_meta()
    vocab = {p: i for i, p in enumerate(meta2['parts_vocab'])}
    n_rows = len(meta2['parts_offsets']) - 1
    PM = sparse.csr_matrix((np.ones(len(meta2['parts_ids']), bool), meta2['parts_ids'], meta2['parts_offsets']), shape=(n_rows, len(vocab))).tocsc()
    import re as _re
    def elem_mask(sym):
        rx = _re.compile(r'\[(?:\d+)?' + _re.escape(sym) + r'(?![a-z])')
        idx = [j for p, j in vocab.items() if rx.search(p)]
        if not idx: return np.zeros(n_rows, bool)
        return np.asarray(PM[:, idx].sum(1)).ravel() > 0
    def part_mask(smiles):
        j = vocab.get(smiles)
        if j is None: return np.zeros(n_rows, bool)
        return np.asarray(PM[:, j].todense()).ravel().astype(bool)
    print(f'reagent roles, fragments and parts read ({time.time() - t0:.0f} s)', flush=True)
    z = np.load(RV.counts_path(a.counts)); S, P = z['S'], z['P']
    cf = [str(x) for x in z['fragments']]
    if cf != frs[:len(cf)] and cf != frs: print('warning: the count table and the feature table list different fragments')
    present = S >= 1
    cons = S > P
    cidx = {f: i for i, f in enumerate(cf)}
    ridx = {r: i for i, r in enumerate(roles)}
    fidx = {f: i for i, f in enumerate(frs)}
    out = {}
    for adv in D.ADVISORY:
        if not adv['reagent'] or not adv['groups']: continue
        mask = np.zeros(S.shape[0], bool); used = []
        if not adv['reagent']: mask[:] = True; used.append('any reagent')
        for alt in adv['reagent']:
            if any(t.startswith('solvent:') for t in alt): continue
            m = np.ones(S.shape[0], bool); ok = True
            for t in alt:
                neg = t.startswith('not:')
                kind, _, val = (t[4:] if neg else t).partition(':')
                if kind == 'role' and val in ridx: has = R[:, ridx[val]]
                elif kind == 'frag' and val in fidx: has = RF[:, fidx[val]]
                elif kind == 'smiles': has = part_mask(val)
                elif kind == 'elem': has = elem_mask(val)
                else: ok = False; continue
                m &= ~has if neg else has
            if ok: mask |= m; used.append(' + '.join(alt))
        if not used: out[adv['key']] = {'measured': False, 'why': 'its reagent side names a solvent type: the recorded conditions are not in the feature table'}; continue
        groups = [g for g in adv['groups'] if g in cidx]
        if not groups: continue
        gi = np.array([cidx[g] for g in groups])
        pg = present[:, gi]; cg = cons[:, gi]
        anyp = pg.any(1)
        fires = mask & (pg & ~cg).any(1)                                            # the advisory would be shown: the reagent matches and a group is present that the product did not consume
        n_with = int((mask & anyp).sum()); n_reag = int(mask.sum())
        per = {}
        for j, g in enumerate(groups):
            pw, po = pg[:, j] & mask, pg[:, j] & ~mask
            nw, no = int(pw.sum()), int(po.sum())
            if nw >= 30: per[g] = {'n': nw, 'with': round(float((cg[:, j] & pw).sum()) / nw, 4), 'other': round(float((cg[:, j] & po).sum()) / max(no, 1), 4), 'n_other': no}
        pw_all = int((pg & mask[:, None]).sum()); po_all = int((pg & ~mask[:, None]).sum())
        out[adv['key']] = {'measured': True, 'reagent_alternatives_used': used, 'reactions_with_reagent': n_reag, 'reactions_with_reagent_and_a_group': n_with, 'would_be_shown': round(float(fires.sum()) / S.shape[0], 4),
                           'pairs_with': pw_all, 'consumed_with': round(float((cg & mask[:, None]).sum()) / max(pw_all, 1), 4),
                           'pairs_other': po_all, 'consumed_other': round(float((cg & ~mask[:, None]).sum()) / max(po_all, 1), 4), 'per_group': per}
    json.dump(out, open(a.out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print(f'wrote {a.out}: {sum(1 for v in out.values() if v.get("measured"))} advisories measured ({time.time() - t0:.0f} s)')
    for k, v in out.items():
        if v.get('measured'): print(f'   {k:36s} n={v["reactions_with_reagent_and_a_group"]:7d}  group consumed {v["consumed_with"]:.2f} with the reagent vs {v["consumed_other"]:.2f} with other reagents; shown in {100 * v["would_be_shown"]:.1f} % of reactions')


if __name__ == '__main__':
    main()
