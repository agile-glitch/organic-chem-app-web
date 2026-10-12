"""Fragment Atlas, shared core: how a functional group IN ITS MOLECULE is described, identically here and in js/atlas.js.

A fragment instance = a functional group (one of the 205 named fragments of the Fragments sheet) at a particular place in a
particular molecule. Its identity is the group plus the surroundings of its reacting atom, written at four radii:

    e0  the reacting atom alone               (element, aromatic, heavy-atom degree, H count, charge, in a ring)
    e1  + its neighbours and the bonds to them
    e2  + their neighbours                    (Weisfeiler-Lehman refinement, 32-bit FNV-1a, see wl_hashes)
    e3  + one more shell

The same ester in two different molecules shares e0 and usually e1, and differs at e2 or e3, so the number of distinct
fragments the program can tell apart is not a list somebody wrote: it is every environment that occurs (millions to billions
of possible ones). What is known about a fragment comes from the most specific radius that has enough recorded reactions
behind it, backing off to the coarser ones (tools/atlas_mine.py counts, js/atlas.js looks up).

Alongside the key each instance carries the continuous facts a chemist reasons with: how crowded the reacting atom is,
its Gasteiger partial charge, the ring it sits in and how far the ring forces its bond angle from the ideal, and the pKa of
the group measured against the reagent's. tools/atlas_eval.py measures what each of them adds on reactions held out by scaffold.

Parity: every function here that has a JavaScript twin is checked by tools/atlas_parity.js against this file.
"""
import os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

SKIP_ROLES = {'counter-ion / spectator ion', 'solvent', 'salt (work-up or additive)', 'other', 'drying agent'}
RADII = 4                                      # e0..e3


# ------------------------------------------------------------------ hashing (the JS twin uses Math.imul)
def fnv32(s):
    h = 0x811c9dc5
    for b in s.encode('utf-8'):
        h = ((h ^ b) * 0x01000193) & 0xffffffff
    return h


def mix(*parts):
    """one 32-bit key from several strings / ints (the page builds table keys the same way)"""
    return fnv32('|'.join(str(p) for p in parts))


# ------------------------------------------------------------------ the sheets
_FRAG = None
_ROLES = None


def fragments():
    """[(id, name, smarts, compiled query)] of the Fragments sheet, sorted by id"""
    global _FRAG
    if _FRAG is None:
        from rdkit import Chem
        import learn_reactions as L
        _FRAG = []
        for fid, name, sm in sorted(L.fragments(), key=lambda f: f[0]):
            q = Chem.MolFromSmarts(sm)
            if q is not None: _FRAG.append((fid, name, sm, q))
    return _FRAG


def _first_number(text):
    m = re.search(r'-?\d+(?:\.\d+)?', str(text or '').replace('−', '-'))
    return float(m.group(0)) if m else None


_FRAG_PKA = None


def fragment_pka():
    """{fragment id: pKa of its most acidic H} (first number of the sheet cell, the sheet's own 'first number wins' rule)"""
    global _FRAG_PKA
    if _FRAG_PKA is None:
        import learn_reactions as L
        _FRAG_PKA = {}
        for f in L.read_sheet(os.path.join(L.APP, 'molecule_data.xlsx'), 'Fragments'):
            _FRAG_PKA[str(f['Fragment ID']).strip()] = (_first_number(f.get('pKa (its most acidic H)')), _first_number(f.get('pKa of its conjugate acid (how basic)')))
    return _FRAG_PKA


def reagent_roles():
    """{canonical SMILES of one reagent part: (role, pKa of the conjugate acid or None, name)}"""
    global _ROLES
    if _ROLES is None:
        from rdkit import Chem
        import learn_reactions as L
        _ROLES = {}
        for r in L.read_sheet(os.path.join(L.APP, 'reaction_rules.xlsx'), 'Reagents'):
            s = str(r.get('SMILES (one reagent part)') or '').strip()
            m = Chem.MolFromSmiles(s) if s else None
            if m is None: continue
            _ROLES[Chem.MolToSmiles(m)] = (str(r.get('Role') or 'other').strip(), _first_number(r.get('pKa of conjugate acid (bases)')),
                                           str(r.get('Name') or '').strip())
    return _ROLES


# ------------------------------------------------------------------ environments
def _atom_inv(a):
    return f'{a.GetAtomicNum()},{int(a.GetIsAromatic())},{sum(1 for n in a.GetNeighbors() if n.GetAtomicNum() > 1)},{a.GetTotalNumHs()},{a.GetFormalCharge()},{int(a.IsInRing())}'


def _bond_cls(b):
    return 'A' if b.GetIsAromatic() else str(int(b.GetBondTypeAsDouble()))


def wl_hashes(m, radii=RADII):
    """hashes[k][i] = environment hash of atom i at radius k (k = 0 .. radii-1). Heavy atoms only count as neighbours."""
    n = m.GetNumAtoms()
    cur = [fnv32(_atom_inv(m.GetAtomWithIdx(i))) for i in range(n)]
    out = [cur]
    for _ in range(radii - 1):
        nxt = []
        for i in range(n):
            a = m.GetAtomWithIdx(i)
            nb = sorted(f'{_bond_cls(b)}:{cur[b.GetOtherAtomIdx(i)]}' for b in a.GetBonds() if b.GetOtherAtom(a).GetAtomicNum() > 1)
            nxt.append(fnv32(f'{cur[i]}[' + ','.join(nb) + ']'))
        cur = nxt
        out.append(cur)
    return out


# ------------------------------------------------------------------ the reacting atom of a fragment match (effects.py pick)
def _order(b):
    return 1.5 if b.GetIsAromatic() else b.GetBondTypeAsDouble()


def reactive_atom(m, atoms):
    A = m.GetAtomWithIdx
    heavy = lambda i: [n.GetIdx() for n in A(i).GetNeighbors() if n.GetAtomicNum() > 1]
    bo = lambda i, j: _order(m.GetBondBetweenAtoms(i, j))
    isC = lambda i: A(i).GetSymbol() == 'C'
    tests = (lambda i: isC(i) and any(A(j).GetSymbol() in ('O', 'N', 'S') and bo(i, j) >= 2 for j in heavy(i)),
             lambda i: isC(i) and any(A(j).GetSymbol() in ('F', 'Cl', 'Br', 'I') for j in heavy(i)),
             lambda i: A(i).GetSymbol() in ('N', 'O', 'S') and A(i).GetTotalNumHs() > 0,
             lambda i: A(i).GetSymbol() in ('N', 'O', 'S'),
             lambda i: isC(i) and any(A(j).GetSymbol() == 'C' and bo(i, j) >= 2 and bo(i, j) != 1.5 for j in heavy(i)))
    for t in tests:
        for i in atoms:
            if t(i): return i
    return atoms[0] if atoms else None


# ------------------------------------------------------------------ continuous descriptors of the reacting atom
def steric_score(m, x):
    """effects.py steric(): heavy neighbours (not a doubly bonded O/S) + the most further heavy neighbours any of them carries
    (an aromatic neighbour: only the substituents on the ring atoms next to it)"""
    A = m.GetAtomWithIdx
    heavy = lambda i: [n.GetIdx() for n in A(i).GetNeighbors() if n.GetAtomicNum() > 1]
    ri = m.GetRingInfo()
    partner = [j for j in heavy(x) if not (A(j).GetSymbol() in ('O', 'S') and _order(m.GetBondBetweenAtoms(x, j)) == 2)]
    k = 0
    for j in partner:
        if A(j).GetIsAromatic():
            ring = next((list(r) for r in ri.AtomRings() if len(r) == 6 and j in r and all(A(q).GetIsAromatic() for q in r)), None)
            if ring:
                p = ring.index(j)
                nb = {ring[(p + 1) % 6], ring[(p - 1) % 6]}
                kk = sum(len([q for q in heavy(o) if q not in ring]) for o in nb)
            else:
                kk = 0
        else:
            kk = len([q for q in heavy(j) if q != x])
        k = max(k, kk)
    return len(partner) + k


def smallest_ring(m, x):
    sizes = [len(r) for r in m.GetRingInfo().AtomRings() if x in r]
    return min(sizes) if sizes else 0


def hybridization(a):
    """1 = sp, 2 = sp2 / aromatic, 3 = sp3 (heavy-atom view: from the bonds, not RDKit's lone-pair rules)"""
    if a.GetIsAromatic(): return 2
    orders = [b.GetBondTypeAsDouble() for b in a.GetBonds()]
    if any(o == 3 for o in orders) or sum(1 for o in orders if o == 2) >= 2: return 1
    if any(o == 2 for o in orders): return 2
    return 3


def angle_strain(m, x):
    """how far the ring forces the bond angle at x from the ideal for its hybridization, in degrees (0 for a chain atom):
    regular-polygon interior angle (n-2)*180/n against 109.5 (sp3), 120 (sp2) or 180 (sp)"""
    n = smallest_ring(m, x)
    if not n: return 0.0
    ideal = {3: 109.5, 2: 120.0, 1: 180.0}[hybridization(m.GetAtomWithIdx(x))]
    return abs((n - 2) * 180.0 / n - ideal)


def gasteiger(m):
    """list of partial charges (None where RDKit has none); the app's molinfo.js reproduces these to 5e-11"""
    from rdkit.Chem import AllChem
    try:
        AllChem.ComputeGasteigerCharges(m)
    except Exception:
        return [None] * m.GetNumAtoms()
    out = []
    for a in m.GetAtoms():
        q = a.GetDoubleProp('_GasteigerCharge') if a.HasProp('_GasteigerCharge') else None
        out.append(None if q is None or q != q or abs(q) == float('inf') else q)
    return out


# ------------------------------------------------------------------ fragment instances of a molecule
def instances(m, frag_list=None, cap=8):
    """[(fragment id, matched atoms, reacting atom)] for every match of every fragment (at most `cap` matches per fragment)"""
    out = []
    for fid, _name, _sm, q in (frag_list or fragments()):
        try:
            ms = m.GetSubstructMatches(q, uniquify=True, maxMatches=cap)
        except Exception:
            continue
        for mt in ms:
            x = reactive_atom(m, list(mt))
            if x is not None: out.append((fid, tuple(mt), x))
    return out


def reagent_class(role_list):
    """the jobs a reagent's active parts do, sorted and joined (the Reagents-sheet roles; tools/learned_to_js.role_signature)"""
    s = sorted({r for r in role_list if r not in SKIP_ROLES})
    return '|'.join(s[:3]) if s else '-'
