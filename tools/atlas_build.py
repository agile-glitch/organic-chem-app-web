"""Fragment Atlas, step 3: turn the mined instances into the tables the page loads (data/atlas/).

    python tools/atlas_build.py <instances.npz> <atlas_eval.json> [--out data/atlas]

Every instance (a functional group in a molecule that reacted, and whether the reaction transformed it) adds one to the counts of
ten table families (tools/atlas_keys.py): the group alone, with the kind of reagent, with its surroundings at radius 1, 2 and 3,
with and without the reagent, and with a named reagent. An entry is (key, met, transformed). Entries met fewer times than the
minimum of their family are dropped (they would only add noise and size); the page treats a missing entry as "never met" and backs
off to the coarser one, so a context nobody has recorded is still judged, from what its parts are like.

Layout, so a page only loads what a molecule needs:
  data/atlas/meta.js     counts, thresholds, the combining model (fitted by tools/atlas_eval.py), the measured figures
  data/atlas/f/NNN.js    every entry of fragment number NNN (the order of the Fragments sheet):
                         gzip( u32 count, then per entry, sorted by key: varint(key - previous key), varint(met), varint(transformed) ), in base64
"""
import argparse, base64, gzip, json, os, struct, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atlas_keys as K


def varint(v, out):
    v = int(v)
    while v >= 0x80:
        out.append((v & 0x7f) | 0x80); v >>= 7
    out.append(v)


def pack_entries(keys, n, h):
    order = np.argsort(keys, kind='stable')
    keys, n, h = keys[order], n[order], h[order]
    out = bytearray(struct.pack('<I', len(keys)))
    prev = 0
    for k, a, b in zip(keys.tolist(), n.tolist(), h.tolist()):
        varint(k - prev, out); varint(a, out); varint(b, out); prev = k
    return bytes(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('npz'); ap.add_argument('eval_json')
    ap.add_argument('--out', default=os.path.join(os.path.dirname(HERE), 'data', 'atlas'))
    a = ap.parse_args()
    t0 = time.time()
    d = dict(np.load(a.npz, allow_pickle=True))
    ev = json.load(open(a.eval_json, encoding='utf-8'))
    y = d['y'].astype(np.float64)
    split = d['split_mol'] if ev.get('split') == 'molecule' else d['split']
    train_base = float(y[split == 2].mean())             # the prior the combining model was fitted with (tools/atlas_eval.py)
    P = K.prepare(d)
    keys = K.family_keys(P)
    nfrag = len(d['frag_ids'])
    fam_arrays = {}                                     # family -> (keys, F of each key, met, transformed)
    for fam in ('F', 'R0', 'R1', 'R2', 'R3', 'U1', 'U2', 'U3'):
        kk = keys[fam]
        u, first, inv = np.unique(kk, return_index=True, return_inverse=True)
        n = np.bincount(inv); h = np.bincount(inv, weights=y).astype(np.int64)
        keep = n >= K.MIN_N[fam]
        fam_arrays[fam] = (u[keep], P['F'][first][keep], n[keep], h[keep])
    for fam in ('P0', 'P1'):
        kk = np.concatenate([keys[(fam, j)][P['has'][j]] for j in range(3)])
        yy = np.concatenate([y[P['has'][j]] for j in range(3)])
        ff = np.concatenate([P['F'][P['has'][j]] for j in range(3)])
        u, first, inv = np.unique(kk, return_index=True, return_inverse=True)
        n = np.bincount(inv); h = np.bincount(inv, weights=yy).astype(np.int64)
        keep = n >= K.MIN_N[fam]
        fam_arrays[fam] = (u[keep], ff[first][keep], n[keep], h[keep])
    os.makedirs(os.path.join(a.out, 'f'), exist_ok=True)
    for old in os.listdir(os.path.join(a.out, 'f')):
        os.remove(os.path.join(a.out, 'f', old))
    shard_entries, shard_bytes, total = {}, {}, 0
    for f in range(nfrag):
        ks, ns, hs = [], [], []
        for fam, (u, ff, n, h) in fam_arrays.items():
            m = ff == f
            ks.append(u[m]); ns.append(n[m]); hs.append(h[m])
        k = np.concatenate(ks); n = np.concatenate(ns); h = np.concatenate(hs)
        if not len(k): continue
        # two families sharing a 32-bit key by chance: add their counts (a 1-in-a-billion event per pair)
        u, inv = np.unique(k, return_inverse=True)
        if len(u) != len(k):
            n = np.bincount(inv, weights=n).astype(np.int64); h = np.bincount(inv, weights=h).astype(np.int64); k = u
        raw = pack_entries(k, n, h)
        b64 = base64.b64encode(gzip.compress(raw, 9)).decode('ascii')
        name = f'{f:03d}'
        with open(os.path.join(a.out, 'f', name + '.js'), 'w', encoding='utf-8') as fh:
            fh.write(f'(window.ATLAS_PARTS=window.ATLAS_PARTS||{{}})["f/{name}"]={{"gz":"{b64}"}};\n')
        shard_entries[f] = int(len(k)); shard_bytes[f] = len(b64); total += len(k)
    meta = {
        'version': 1, 'built': time.strftime('%Y-%m-%d'), 'm': K.M_CHAIN, 'minN': K.MIN_N, 'families': K.FAM,
        'instances': int(len(y)), 'base': train_base, 'allBase': float(y.mean()), 'entries': total, 'fragments': [str(i) for i in d['frag_ids']],
        'shardEntries': shard_entries, 'shardBytes': shard_bytes,
        'source': 'USPTO-MIT (Jin, Coley et al. 2017; US patent reactions 1976-2016 text-mined by D. M. Lowe, CC0), every atom mapped',
        'model': ev['model'], 'eval': {k: ev[k] for k in ('split', 'instances', 'train', 'validation', 'test', 'results', 'by_novelty', 'site_choice') if k in ev},
    }
    with open(os.path.join(a.out, 'meta.js'), 'w', encoding='utf-8') as fh:
        fh.write('window.ATLAS_META=' + json.dumps(meta, separators=(',', ':'), ensure_ascii=False) + ';\n')
    mb = sum(shard_bytes.values()) / 1e6
    print(f'{total:,} entries in {len(shard_entries)} shards, {mb:.1f} MB (base64), {time.time() - t0:.0f}s -> {a.out}')
    big = sorted(shard_bytes.items(), key=lambda t: -t[1])[:5]
    print('largest shards:', ', '.join(f'F{f}: {b / 1e6:.2f} MB' for f, b in big))


if __name__ == '__main__':
    main()
