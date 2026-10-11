"""Generalise a learned template: the same reaction centre with less of its surroundings required.

(The idea and level 1-2 come from the earlier lever experiment, tools/template_general.py in the xylose worktree; level 3 is new.)

A learned template (rdchiral style) is  reactants>>products  where
  * the atoms that CHANGE carry a full description        [C;H0;D3;+0:1]   (H count, degree, charge)
  * unchanged neighbours carry just the element           [C:2]  [c:3]  [#7;a:4]
  * "special group" atoms (unchanged, but kept for context) carry degree / H:   [O;D1;H0:5]  [Cl;D1;H0:6]
  * atoms that LEAVE are unmapped                         C-O-
so a template learned from one crowded molecule demands that molecule's whole neighbourhood, and 28 % of the reactions of
unseen patents use a template no training reaction used (only 2 points of that are won back by stripping plain neighbours).

  level 0  the template as learned
  level 1  strip PLAIN unchanged leaf atoms, outside in
  level 2  level 1, and the changed atoms on the reactant side described by element and aromaticity only
  level 3  strip EVERY unchanged leaf atom (plain or special group), outside in, as long as it is the same atom with the same
           bond on both sides; what is left is the changed atoms, the leaving atoms and whatever connects them.
           RDKit carries every atom the template no longer mentions through to the product, so the product is unchanged.
  level 4  level 3 + level 2's loosening of the changed atoms

generalise(smarts, level) returns a reaction SMARTS (None when it cannot be built or nothing changed).
class_key(smarts, level) is the template with its atom-map numbers removed: templates with the same key are the same chemistry.
"""
import re
from rdkit import Chem
from rdkit.Chem import AllChem

PLAIN = re.compile(r'^\[(?:#\d+|[A-Za-z][a-z]?)(?:;a)?:\d+\]$')      # [C:2]  [c:3]  [#7:4]  [#7;a:5] : an unchanged neighbour
FULL = re.compile(r';[HD]\d|\+0|;[+-]')                              # a described atom
MAPNUM = re.compile(r':\d+\]')


def _is_plain(atom):
    return atom.GetAtomMapNum() > 0 and bool(PLAIN.match(atom.GetSmarts()))


def _same_atom(a, b):
    """the same atom described the same way on both sides (map numbers aside)"""
    return MAPNUM.sub(']', a.GetSmarts()) == MAPNUM.sub(']', b.GetSmarts())


def _strip_leaves(left, right, plain_only):
    """remove unchanged mapped atoms that are leaves (degree 1) on both sides, with the same bond to the same neighbour, repeatedly"""
    changed = True
    while changed:
        changed = False
        rmap = {a.GetAtomMapNum(): a.GetIdx() for a in right.GetAtoms() if a.GetAtomMapNum()}
        for a in list(left.GetAtoms()):
            m = a.GetAtomMapNum()
            if not m or a.GetDegree() != 1 or m not in rmap: continue
            b = right.GetAtomWithIdx(rmap[m])
            if b.GetDegree() != 1: continue
            if plain_only:
                if not (_is_plain(a) and _is_plain(b)): continue
            elif not _same_atom(a, b): continue
            ln, rn = a.GetNeighbors()[0], b.GetNeighbors()[0]
            if ln.GetAtomMapNum() != rn.GetAtomMapNum() or ln.GetAtomMapNum() == 0: continue
            lb, rb = left.GetBondBetweenAtoms(a.GetIdx(), ln.GetIdx()), right.GetBondBetweenAtoms(b.GetIdx(), rn.GetIdx())
            if lb.GetSmarts() != rb.GetSmarts(): continue
            left.RemoveAtom(a.GetIdx()); right.RemoveAtom(b.GetIdx()); changed = True
            break                                                    # indices moved: rescan
    return left, right


def _simplify_reactant_centres(left):
    """[C;H0;D3;+0:1] -> [C:1]; [c;H0;D3;+0:2] -> [c:2]; plain atoms are kept as they are"""
    for a in left.GetAtoms():
        s = a.GetSmarts()
        if a.GetAtomMapNum() and FULL.search(s):
            m = re.match(r'^\[(#?\w+?)(?:[;&@+\-].*?)?:(\d+)\]$', s)
            if not m: continue
            a.SetQuery(Chem.AtomFromSmarts('[' + m.group(1) + ']'))
            a.SetAtomMapNum(int(m.group(2)))
    return left


def generalise(smarts, level=1):
    try:
        L, R = smarts.split('>>')
        lm, rm = Chem.MolFromSmarts(L), Chem.MolFromSmarts(R)
        if lm is None or rm is None: return None
        lm, rm = Chem.RWMol(lm), Chem.RWMol(rm)
        before = lm.GetNumAtoms()
        if level >= 1: lm, rm = _strip_leaves(lm, rm, plain_only=level in (1, 2))
        if level in (2, 4): lm = _simplify_reactant_centres(lm)
        if level in (1, 3) and lm.GetNumAtoms() == before: return None            # nothing to strip
        if lm.GetNumBonds() == 0 and lm.GetNumAtoms() < 2: return None
        out = Chem.MolToSmarts(lm) + '>>' + Chem.MolToSmarts(rm)
        if out == smarts: return None
        rx = AllChem.ReactionFromSmarts(out)
        if rx is None or rx.Validate()[0] != 0: return None
        return out
    except Exception:
        return None


def class_key(smarts, level):
    g = generalise(smarts, level) if level else smarts
    return MAPNUM.sub(']', g or smarts)


if __name__ == '__main__':
    tests = [
        ('C-C-[O;H0;D2;+0:1]-[C:2]=[O;D1;H0:3]>>[O;D1;H0:3]=[C:2]-[OH;D1;+0:1]', 'CCOC(=O)c1ccccc1'),
        ('O=[N+;H0;D3:1](-[O-])-[c:2]>>[NH2;D1;+0:1]-[c:2]', 'O=[N+]([O-])c1ccc(Cl)cc1'),
        ('Br-[CH2;D2;+0:1]-[c:2].[#7;a:3]:[c:4](-[NH;D2;+0:5]-[C:6]=[O;D1;H0:7]):[#8;a:8]>>[#7;a:3]:[c:4](-[N;H0;D3;+0:5](-[C:6]=[O;D1;H0:7])-[CH2;D2;+0:1]-[c:2]):[#8;a:8]',
         'Cc1coc(NC(=O)C(C)C)n1.Clc1cccc(CBr)c1'),
    ]
    for t, m in tests:
        print('L0', t)
        for lv in (1, 2, 3, 4):
            g = generalise(t, lv); print(f'L{lv}', g)
            if g:
                rx = AllChem.ReactionFromSmarts(g)
                subs = [Chem.MolFromSmiles(x) for x in m.split('.')]
                n = rx.GetNumReactantTemplates()
                import itertools
                ps = set()
                for order in itertools.permutations(range(len(subs)), n):
                    for p in rx.RunReactants(tuple(subs[i] for i in order)):
                        try: ps.add(Chem.MolToSmiles(p[0]))
                        except Exception: pass
                print('   on', m, '->', ps)
