"""The interference engine: which interferences does a candidate product (or a reagent and its molecules) have?

    import interferences as I
    lib = I.Library()
    lib.flags_of_smiles('Nc1ccc([900Cl])cc1')               # site patterns true at a marked (reacting) atom of a product with marks
    I.advisories(ctx, present, consumed)                    # reagent + group advisories

The definitions are in tools/interference_defs.py (the page reads the same rows from the sheet Interferences); the isotope marks are made by tools/int_marks.py
(the page does both the same way: js/interference.js). A site pattern is written for its reacting atom (the first atom); the engine puts the mark on it.
"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import interference_defs as D
import int_marks as M

ANCHOR_PRE = {'aromatic C': '[900c]', 'aromatic n': '[900n]', 'sp3 C': '[900C;X4]', 'carbonyl C': '[900C;X3;$(C=[O,N,S])]', 'alkene C': '[900C;X3;$(C=C)]',
              'nitrogen': '[900N]', 'oxygen': '[900O]', 'sulfur': '[900S]'}


def mark_first_atom(smarts):
    """[c;$(c:c-X)] -> [900c;$(c:c-X)]  (the isotope primitive in front of the first atom's specification)"""
    if not smarts.startswith('['): raise ValueError(smarts)
    return '[' + str(M.MARK) + smarts[1:]


class Library:
    def __init__(self, site=None):
        from rdkit import Chem
        self.Chem = Chem
        self.site = list(site if site is not None else D.SITE)
        self.q = []
        self.by_anchor = {}
        for s in self.site:
            q = Chem.MolFromSmarts(mark_first_atom(s['smarts']))
            if q is None: raise ValueError('bad pattern ' + s['key'] + ': ' + s['smarts'])
            self.q.append(q)
            self.by_anchor.setdefault(s['anchor'], []).append(len(self.q) - 1)
        self.pre = {k: Chem.MolFromSmarts(v) for k, v in ANCHOR_PRE.items()}
        self.keys = [s['key'] for s in self.site]
        self.cache = {}
        self.sub_pre = {}                                    # pattern index -> (prefilter SMARTS, compiled): shared by the three positions of one substituent
        cache = {}
        for i, s in enumerate(self.site):
            pf = s.get('prefilter')
            if pf:
                if pf not in cache:
                    cache[pf] = Chem.MolFromSmarts(pf)
                    if cache[pf] is None: raise ValueError('bad prefilter ' + s['key'] + ': ' + pf)
                self.sub_pre[i] = (pf, cache[pf])

    def flags(self, mol):
        """the keys of the site patterns true at some marked atom of the molecule (an RDKit molecule with the isotope marks), sorted"""
        out = []
        has_sub = {}
        for anchor, idx in self.by_anchor.items():
            pre = self.pre.get(anchor)
            if pre is not None and not mol.HasSubstructMatch(pre): continue
            for i in idx:
                sp = self.sub_pre.get(i)
                if sp is not None:
                    if sp[0] not in has_sub: has_sub[sp[0]] = mol.HasSubstructMatch(sp[1])
                    if not has_sub[sp[0]]: continue
                if mol.HasSubstructMatch(self.q[i]): out.append(self.keys[i])
        return sorted(out)

    def flags_of_smiles(self, smiles):
        """the same for a product SMILES that carries the marks ([900C]...); several molecules joined by '.' give the union of their flags"""
        hit = self.cache.get(smiles)
        if hit is not None: return hit
        Chem = self.Chem
        got = set()
        for part in smiles.split('.'):
            if '900' not in part: continue
            part = _HISO.sub('[H]', part)                 # deuterium and tritium are atoms to RDKit and would count in D / X: read them as hydrogen
            m = Chem.MolFromSmiles(part)
            if m is None: continue
            got.update(self.flags(m))
        out = sorted(got)
        if len(self.cache) < 200000: self.cache[smiles] = out
        return out


# ---------------------------------------------------------------------------------------------------- advisories
_ELEM = {}


def has_element(parts, sym):
    """some reagent part, as typed, has this element as a bracket atom ([Pd], [Pd+2], [Pd-4] ...); meant for metals"""
    if sym not in _ELEM: _ELEM[sym] = re.compile(r'\[(?:\d+)?' + re.escape(sym) + r'(?![a-z])')
    return any(_ELEM[sym].search(p) for p in parts)


_HISO = re.compile(r'\[[23]H\]')


def token_holds(token, ctx):
    """ctx = {'roles': set of Reagents-sheet roles of the reagent parts (role and 'also'), 'frags': set of fragment IDs found in the reagent parts,
              'parts': set of the reagent parts as typed, 'solvent_types': set of Solvents-sheet types of the solvents given}"""
    if token.startswith('not:'): return not token_holds(token[4:], ctx)                 # "not:smiles:[Ni]": the reagent must NOT have this part
    kind, _, val = token.partition(':')
    if kind == 'role': return val in ctx['roles']
    if kind == 'frag': return val in ctx['frags']
    if kind == 'smiles': return val in ctx['parts']
    if kind == 'elem': return has_element(ctx['parts'], val)
    if kind == 'solvent': return val in ctx.get('solvent_types', ())
    return False


def reagent_matches(adv, ctx):
    """an advisory with no reagent side applies whatever the reagent is"""
    return not adv['reagent'] or any(all(token_holds(t, ctx) for t in alt) for alt in adv['reagent'])


def advisories(ctx, present, consumed=(), only_keys=None):
    """the advisories a reagent and the groups present trigger. present = fragment IDs found in the starting materials; consumed = the ones the leading product
    consumed. Returns [{key, groups: [IDs a bystander], target: [IDs the product consumed]}], only those the advisory wants to show (its 'show' field)"""
    out = []
    for a in D.ADVISORY:
        if not a['groups'] or a['needs'] in ('size', 'reagent', 'conditions') or not a.get('page', True): continue
        if only_keys is not None and a['key'] not in only_keys: continue
        if a['needs'] == 'solvent' and not ctx.get('solvent_types'): continue
        if not reagent_matches(a, ctx): continue
        have = set(present) | (set(ctx['frags']) if a['group_where'] in ('reagent', 'both') else set())
        if a['group_where'] == 'reagent': have = set(ctx['frags'])
        hit = [g for g in a['groups'] if g in have]
        if len(hit) < a['min_distinct']: continue
        by, tg = [g for g in hit if g not in consumed], [g for g in hit if g in consumed]
        shown = by if a['show'] == 'bystander' else tg if a['show'] == 'target' else hit
        if not shown: continue
        out.append({'key': a['key'], 'groups': by, 'target': tg})
    return out
