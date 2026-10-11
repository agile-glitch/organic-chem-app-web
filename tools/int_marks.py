"""Mark the atoms a template changes, so that the product says WHERE the reaction happened.

A learned template (tools/learn_reactions.py, written by rdchiral) has a convention: the atoms that CHANGE carry a strict specification
([N;H0;D3;+0:4]: hydrogens, degree, charge), the atoms around them that stay the same are plain ([c:2], [C:3]); atoms the reaction takes away or
adds have no map number. RDKit copies an isotope written on a PRODUCT-side template atom onto the product atom, so giving every changed product atom
an isotope (900) makes the reaction centre visible in the product molecule (RDKit.js does the same: no atom-mapping API is needed):

    mark_template('[C:1](=[O:2])[OH].[N;H2:3][c:4]>>[C;H0;D3;+0:1](=[O:2])[N;H1;D2;+0:3][c:4]')
        -> '...>>[900C;H0;D3;+0:1](=[O:2])[900N;H1;D2;+0:3][c:4]'

js/interference.js has the same function; tools/int_check_js.py compares them. strip_marks() removes the isotope again from a product SMILES.
"""
import re

MARK = 900
_STRICT = re.compile(r';(?:D\d|H\d|[+-]\d?)')           # H, D or a charge in the specification = a changed atom
_ATOM = re.compile(r'\[([^\[\]]*)\]')


def _split_sides(smarts):
    r, _, p = smarts.partition('>>')
    return r, p


def changed_maps(smarts):
    """the atom-map numbers of the atoms the template really changes: an atom whose specification differs between the reactant and the product side, or
    that gains or loses a neighbour or changes a bond (so the carbonyl O of a special group that merely sits there is not one). {} if RDKit cannot read it"""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    try: rxn = AllChem.ReactionFromSmarts(smarts)
    except Exception: return set()
    def side(mols):
        info = {}
        for t in mols:
            for a in t.GetAtoms():
                mp = a.GetAtomMapNum()
                if not mp: continue
                spec = re.sub(r':\d+\]$', ']', a.GetSmarts())
                nb = []
                for b in a.GetBonds():
                    o = b.GetOtherAtom(a)
                    omp = o.GetAtomMapNum()
                    nb.append((omp if omp else 'u:' + re.sub(r':\d+\]$', ']', o.GetSmarts()), b.GetSmarts()))
                info[mp] = (spec, sorted(nb, key=str))
        return info
    R = side([rxn.GetReactantTemplate(i) for i in range(rxn.GetNumReactantTemplates())])
    P = side([rxn.GetProductTemplate(i) for i in range(rxn.GetNumProductTemplates())])
    return {mp for mp in P if mp in R and R[mp] != P[mp]}


def mark_product_side(p, mark_new=True, only=None):
    """the product side with an isotope on every changed atom (and, with mark_new, on every atom the template creates: bracket atoms without a map
    number and unbracketed ones such as Cl). only = the map numbers that count as changed (changed_maps); None: every atom with a strict
    specification (rdchiral also writes those for the atoms of a special group that does not change)"""
    out, i, n = [], 0, len(p)
    while i < n:
        ch = p[i]
        if ch == '[':
            depth, j = 1, i + 1
            while j < n and depth:
                if p[j] == '[': depth += 1
                elif p[j] == ']': depth -= 1
                j += 1
            body = p[i + 1:j - 1]
            mm = re.search(r':(\d+)$', body)
            changed = ((int(mm.group(1)) in only) if only is not None else bool(_STRICT.search(body))) if mm else mark_new
            if changed and not re.match(r'\d', body):
                body = f'{MARK}{body}'
            out.append(f'[{body}]'); i = j
        elif p[i:i + 2] in ('Cl', 'Br') or (ch in 'BCNOPSFI' and not (i + 1 < n and p[i + 1] == 'l' and ch == 'C')) or ch in 'bcnops*':
            two = p[i:i + 2] in ('Cl', 'Br')
            tok = p[i:i + 2] if two else ch
            out.append(f'[{MARK}{tok}]' if mark_new else tok); i += len(tok)
        else:
            out.append(ch); i += 1
    return ''.join(out)


def mark_template(smarts, mark_new=True, only=None):
    r, p = _split_sides(smarts)
    return r + '>>' + mark_product_side(p, mark_new, only)


def strip_marks(smiles):
    """a product SMILES without the isotope the marks left ([900C] -> [C], [900CH2] -> [CH2]); canonicalise afterwards"""
    return re.sub(rf'\[{MARK}', '[', smiles)
