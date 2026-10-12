"""Fragment Atlas: the table keys, one definition for the evaluation, the builder and (through tools/atlas_parity.js) the page.

A table entry counts, for one situation, how many fragment instances were met (n) and how many of them the reaction transformed (h).
The situation is the fragment F (0..204, the order of the Fragments sheet) with a context, one of ten FAMILIES:

  F   the fragment alone
  R0  F + the kind of reagent                          (reagent class: the Reagents-sheet roles, sorted, e.g. "base|electrophile / alkylating agent")
  R1  F + reagent class + environment radius 1         (the reacting atom and its neighbours)
  R2  ...                  radius 2
  R3  ...                  radius 3
  U1  F + environment radius 1, any reagent
  U2  F + environment radius 2, any reagent
  U3  F + environment radius 3, any reagent
  P0  F + one named reagent part                       (its InChIKey)
  P1  F + one named reagent part + environment radius 1

All keys are 32-bit: mix32 chains of small integers, so numpy here and Math.imul in js/atlas.js give the same number.
    key = mix32(mix32(mix32(FAMILY, F), a), b)          (a, b = 0 when the family does not use them)
The environment hashes (e1..e3) come from tools/atlas_core.py wl_hashes; class and part hashes are fnv32 of the class text and of the
part's InChIKey. Everything shipped is derived from these definitions; tools/atlas_parity.js checks the page against them.
"""
import numpy as np

FAM = {'F': 1, 'R0': 2, 'R1': 3, 'R2': 4, 'R3': 5, 'U1': 6, 'U2': 7, 'U3': 8, 'P0': 9, 'P1': 10}
C1, C2, C3 = np.uint32(0x9E3779B9), np.uint32(0x85EBCA6B), np.uint32(0xC2B2AE35)
# the page keeps every entry that has been met at least this often (rare deep contexts would only add noise and size)
MIN_N = {'F': 1, 'R0': 1, 'R1': 1, 'U1': 1, 'P0': 1, 'R2': 3, 'R3': 3, 'U2': 3, 'U3': 3, 'P1': 3}
M_CHAIN = 5.0                                           # smoothing: p_k = (h_k + M p_(k-1)) / (n_k + M)


def mix32(a, b):
    """murmur3-style mix of two uint32 values (arrays or scalars); wraps like Math.imul on the page"""
    a = np.atleast_1d(np.asarray(a)).astype(np.uint32)
    b = np.atleast_1d(np.asarray(b)).astype(np.uint32)
    x = a ^ ((b ^ C1) * C2)
    x ^= x >> np.uint32(16)
    x = x * C2
    x ^= x >> np.uint32(13)
    x = x * C3
    x ^= x >> np.uint32(16)
    return x


def key(fam, F, a=0, b=0):
    return mix32(mix32(mix32(np.uint32(FAM[fam]), F), a), b)


def fnv32_text(s):
    h = 0x811c9dc5
    for c in s.encode('utf-8'):
        h = ((h ^ c) * 0x01000193) & 0xffffffff
    return h


def part_hashes(parts):
    """uint32 hash of each reagent part (index 0 = no part): fnv32 of its InChIKey (what the page's keyOf gives)"""
    from rdkit import Chem
    out = np.zeros(len(parts) + 1, np.uint32)
    for i, s in enumerate(parts, start=1):
        m = Chem.MolFromSmiles(s)
        try:
            k = Chem.MolToInchiKey(m) if m is not None else ''
        except Exception:
            k = ''
        out[i] = fnv32_text(k or s)
    return out


def class_hashes(classes):
    return np.array([fnv32_text(c) for c in classes], np.uint32)


def prepare(d):
    """per-instance uint32 arrays used by every family: F, class hash, parts' hashes, environments"""
    cl = class_hashes(list(d['classes']))
    ph = part_hashes(list(d['parts']))
    return {
        'F': d['frag'].astype(np.uint32), 'cls': cl[d['cls']],
        'e': [d['e0'].astype(np.uint32), d['e1'].astype(np.uint32), d['e2'].astype(np.uint32), d['e3'].astype(np.uint32)],
        'p': [ph[d['p1']], ph[d['p2']], ph[d['p3']]], 'has': [d['p1'] > 0, d['p2'] > 0, d['p3'] > 0],
    }


def family_keys(P, idx=None):
    """{family: uint32 keys} for the instances idx (all when None). The part families have one key per part slot: ('P0', j) / ('P1', j)."""
    sel = (lambda v: v) if idx is None else (lambda v: v[idx])
    F, cls, e = sel(P['F']), sel(P['cls']), [sel(x) for x in P['e']]
    k = {'F': key('F', F), 'R0': key('R0', F, cls)}
    for lv in (1, 2, 3):
        k[f'R{lv}'] = key(f'R{lv}', F, cls, e[lv])
        k[f'U{lv}'] = key(f'U{lv}', F, e[lv])
    for j in range(3):
        pj = sel(P['p'][j])
        k[('P0', j)] = key('P0', F, pj)
        k[('P1', j)] = key('P1', F, pj, e[1])
    return k
