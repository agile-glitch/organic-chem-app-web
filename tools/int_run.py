"""Run a learned template the way the page does, with the atoms it changes marked, and read the interference flags off each product.

    outs = run_marked(smarts, subs, extras, screen, lib)
    -> [(product SMILES (canonical, as eval_top1.run_used gives them), number of starting materials used, marked product SMILES, flags)]

This is tools/eval_top1.run_used with one difference: the template that RDKit runs has an isotope on every product atom the template changes
(tools/int_marks.py), so each product tells where the reaction happened; the canonical SMILES are taken with the marks removed, so the products
are the ones the unmarked run gives (tools/int_marks.py was checked for that: 4000 of 4000). `lib` is a tools/interferences.Library.
"""
import itertools
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import int_marks as M

_MAPS = {}                         # template SMARTS -> the marked template (cached; the unmarked one when RDKit cannot read the template)
_LIMIT = 60000


def marked_template(smarts):
    m = _MAPS.get(smarts)
    if m is None:
        try: m = M.mark_template(smarts, only=M.changed_maps(smarts))
        except Exception: m = smarts
        if len(_MAPS) > _LIMIT: _MAPS.clear()
        _MAPS[smarts] = m
    return m


def run_marked(smarts, subs, extras, screen, lib, want=None):
    import learn_reactions as L
    from rdkit import Chem
    from rdkit.Chem import AllChem
    msmarts = marked_template(smarts)
    rx = L._W['rxn'].get(msmarts)
    if rx is None:
        try: rx = AllChem.ReactionFromSmarts(msmarts); rx.Initialize()
        except Exception: rx = False
        L._W['rxn'][msmarts] = rx
    if rx is False: return []
    n = rx.GetNumReactantTemplates()
    pool = list(subs) + [x for x in extras if x not in subs]
    if n > len(pool): return []
    mols = [Chem.MolFromSmiles(s) for s in pool]
    if any(m is None for m in mols[:len(subs)]): return []
    if screen is not None and not L.could_fit(rx, screen): return []
    outs, seen = [], set()
    for order in itertools.permutations(range(len(pool)), n):
        if min(order) >= len(subs): continue
        try: res = rx.RunReactants(tuple(mols[i] for i in order), 50)
        except Exception: continue
        used = len({i for i in order if i < len(subs)})
        for prods in res:
            marked, cs = [], []
            for pm in prods:
                try:
                    ms = Chem.MolToSmiles(pm)
                    cs.append(L._canon(M.strip_marks(ms)))
                    marked.append(ms)
                except Exception:
                    cs.append(None)
            if not cs or any(c is None for c in cs): continue
            key = '.'.join(sorted(cs))
            if key in seen: continue
            seen.add(key)
            flags = lib.flags_of_smiles('.'.join(marked)) if lib is not None else []
            outs.append((cs, used, marked, flags))
    return outs


def clear_caches(L, limit=2500):
    """learn_reactions keeps every compiled template (unbounded) and the quick-screen fingerprints keyed by id(rx): empty both together (eval_top1.clear_template_caches)"""
    if len(L._W['rxn']) > limit:
        L._W['rxn'].clear()
        L._W.get('qfp', {}).clear()
