"""Steric and electronic effects of each fragment IN ITS MOLECULE, read at the fragment's reacting atom.
The Python twin of effectsOfMol in js/reactions.js (same rules), used by tools/learned_molecules.py.

Reacting atom of a fragment match: the C of a C=O / C=N / C#N; else a C carrying a halogen; else an N, O or S (one with
an H first); else a C of a C=C / C#C; else the first atom.
Steric score = the reacting atom's heavy neighbours (not counting a doubly bonded O/S partner) + the most of those
neighbours' own further heavy neighbours (an aromatic neighbour: only its ortho substituents count, the ring is flat):
<= 2 open, 3 moderately hindered, 4 hindered, >= 5 very hindered.
Taft: each group R on the reacting atom (not its own ring, not a =O/=S partner) gets sigma* (polar) and Es (steric) from TAFT.
Electronic: Gasteiger-Marsili partial charge (RDKit); for an atom on, or bonded to, a benzene ring the meta/para
substituents' Hammett sigma (tools/pka.py SIGMA) and their sum; conjugation; a carbon's oxidation state.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# Taft parameters of a group R on the reacting atom (CH3 = 0): sigma* (polar: > 0 pulls electrons, < 0 pushes) and Es (steric:
# more negative = bulkier). Taft, in Steric Effects in Organic Chemistry (1956); Hansch, Leo, Taft, Chem. Rev. 91, 165 (1991);
# Hansch & Leo, Substituent Constants for Correlation Analysis (1979). None = no reliable value.
# [name, SMARTS (first atom = the one bonded to the reacting atom), sigma*, Es]; the first match wins, so specific groups come first.
TAFT = [
    ('methyl', '[CH3X4]', 0.00, 0.00),
    ('ethyl', '[CH2X4;D2][CH3]', -0.10, -0.07),
    ('n-propyl', '[CH2X4;D2][CH2X4;D2][CH3]', -0.115, -0.36),
    ('n-butyl', '[CH2X4;D2][CH2X4;D2][CH2X4;D2][CH3]', -0.13, -0.39),
    ('isopropyl', '[CH1X4;D3]([CH3])[CH3]', -0.19, -0.47),
    ('isobutyl', '[CH2X4;D2][CH1X4;D3]([CH3])[CH3]', -0.125, -0.93),
    ('sec-butyl', '[CH1X4;D3]([CH3])[CH2X4;D2][CH3]', -0.21, -1.13),
    ('tert-butyl', '[CX4;D4]([CH3])([CH3])[CH3]', -0.30, -1.54),
    ('neopentyl', '[CH2X4;D2][CX4;D4]([CH3])([CH3])[CH3]', -0.165, -1.74),
    ('cyclohexyl', '[CH1X4;D3;R1]1[CH2;D2][CH2;D2][CH2;D2][CH2;D2][CH2;D2]1', -0.15, -0.79),
    ('benzyl', '[CH2X4;D2][c;D3]1[cH][cH][cH][cH][cH]1', 0.215, -0.38),
    ('trifluoromethyl', '[CX4;D4](F)(F)F', 2.61, -1.16),
    ('trichloromethyl', '[CX4;D4](Cl)(Cl)Cl', 2.65, -2.06),
    ('dichloromethyl', '[CH1X4;D3](Cl)Cl', 1.94, -1.54),
    ('fluoromethyl', '[CH2X4;D2]F', 1.10, -0.24),
    ('chloromethyl', '[CH2X4;D2]Cl', 1.05, -0.24),
    ('bromomethyl', '[CH2X4;D2]Br', 1.00, -0.27),
    ('iodomethyl', '[CH2X4;D2]I', 0.85, -0.37),
    ('hydroxymethyl', '[CH2X4;D2][OX2H1]', 0.555, -0.03),
    ('methoxymethyl', '[CH2X4;D2][OX2;D2][CH3]', 0.52, -0.19),
    ('phenyl', '[c;D3]1[cH][cH][cH][cH][cH]1', 0.60, -1.01),
    ('aryl (taken as phenyl; its ring substituents are in the Hammett reading)', '[c;D3]', 0.60, -1.01),
    ('fluoro', '[F]', None, 0.78),
    ('chloro', '[Cl]', None, 0.27),
    ('bromo', '[Br]', None, 0.08),
    ('iodo', '[I]', None, -0.16),
    ('hydroxy', '[OX2H1]', None, 0.69),
    ('methoxy', '[OX2;D2][CH3]', None, 0.69),
]
# a carbon group not in the table: Es estimated from its branching (sigma* left out: it depends on what the group carries)
TAFT_BRANCH = {2: ('a CH2R group (Es taken as n-propyl)', -0.36), 3: ('a CHR2 group (Es taken as isopropyl)', -0.47),
               4: ('a CR3 group (Es taken as tert-butyl)', -1.54)}

EN = {'H': 2.20, 'B': 2.04, 'C': 2.55, 'N': 3.04, 'O': 3.44, 'F': 3.98, 'Si': 1.90, 'P': 2.19, 'S': 2.58, 'Cl': 3.16, 'Br': 2.96, 'I': 2.66}


def _order(b):
    return 1.5 if b.GetIsAromatic() else b.GetBondTypeAsDouble()


def effects_of(m, frags):
    """[{fragment, atom, steric: {score, cls, text}, electronic: {charge, sigma, conjugated, text}}]"""
    from rdkit.Chem import AllChem
    import pka as PK
    AllChem.ComputeGasteigerCharges(m)
    A = m.GetAtomWithIdx
    heavy = lambda i: [n.GetIdx() for n in A(i).GetNeighbors() if n.GetAtomicNum() > 1]
    bo = lambda i, j: _order(m.GetBondBetweenAtoms(i, j))
    rings = [r for r in m.GetRingInfo().AtomRings() if len(r) == 6 and all(A(k).GetIsAromatic() for k in r)]
    ring_of = lambda i: next((list(r) for r in rings if i in r), None)
    dist = lambda ring, a, b: min(abs(ring.index(a) - ring.index(b)), 6 - abs(ring.index(a) - ring.index(b)))

    def pick(atoms):
        isC = lambda i: A(i).GetSymbol() == 'C'
        for test in (lambda i: isC(i) and any(A(j).GetSymbol() in ('O', 'N', 'S') and bo(i, j) >= 2 for j in heavy(i)),
                     lambda i: isC(i) and any(A(j).GetSymbol() in ('F', 'Cl', 'Br', 'I') for j in heavy(i)),
                     lambda i: A(i).GetSymbol() in ('N', 'O', 'S') and A(i).GetTotalNumHs() > 0,
                     lambda i: A(i).GetSymbol() in ('N', 'O', 'S'),
                     lambda i: isC(i) and any(A(j).GetSymbol() == 'C' and bo(i, j) >= 2 and bo(i, j) != 1.5 for j in heavy(i))):
            for i in atoms:
                if test(i): return i
        return atoms[0] if atoms else None

    def steric(x):
        partner = [j for j in heavy(x) if not (A(j).GetSymbol() in ('O', 'S') and bo(x, j) == 2)]
        k = 0
        for j in partner:
            if A(j).GetIsAromatic():
                ring = ring_of(j)
                kk = sum(len([q for q in heavy(o) if q not in ring]) for o in ring if o != j and dist(ring, o, j) == 1) if ring else 0
            else:
                kk = len([q for q in heavy(j) if q != x])
            k = max(k, kk)
        score = len(partner) + k
        cls = 'open' if score <= 2 else 'moderately hindered' if score == 3 else 'hindered' if score == 4 else 'very hindered'
        return {'score': score, 'cls': cls, 'text': f'{cls} (score {score}: {len(partner)} heavy neighbour{"" if len(partner) == 1 else "s"}; the most crowded one carries {k} more)'}

    def oxidation(x):                                   # carbon: +order to each more electronegative neighbour, -1 per H, +charge
        a = A(x)
        v = a.GetFormalCharge() - a.GetTotalNumHs()
        for b in a.GetBonds():
            o = b.GetOtherAtom(a); e1, e2 = EN.get(a.GetSymbol()), EN.get(o.GetSymbol())
            if e1 is None or e2 is None or o.GetSymbol() == 'C': continue
            v += _order(b) if e2 > e1 else -_order(b)
        return v

    def electronic(x):
        a = A(x); parts = []
        q = a.GetDoubleProp('_GasteigerCharge') if a.HasProp('_GasteigerCharge') else None
        if q is not None and q == q: parts.append(f'partial charge {q:+.2f}')
        else: q = None
        site = x if a.GetIsAromatic() else next((j for j in heavy(x) if A(j).GetIsAromatic()), None)
        sigma = None
        ring = ring_of(site) if site is not None else None
        if ring and all(A(k).GetSymbol() == 'C' for k in ring):
            t, named = 0.0, []
            for r in ring:
                if r == site: continue
                d = dist(ring, r, site)
                for j in heavy(r):
                    if j in ring or j == x: continue
                    if d == 1: named.append('ortho substituent'); continue
                    hit = None
                    for name, sm, sm_m, sm_p, _ in PK.SIGMA:
                        if any(mt[0] == j for mt in m.GetSubstructMatches(PK._q(sm))): hit = (name, sm_m, sm_p); break
                    if not hit: named.append(f'{"meta" if d == 2 else "para"} group with no σ value'); continue
                    v = hit[1] if d == 2 else hit[2]; t += v
                    named.append(f'{"meta" if d == 2 else "para"}-{hit[0]} (σ {v:+.2f})')
            sigma = round(t, 2)
            parts.append(f'ring: {", ".join(named)}; Σσ {sigma:+.2f} → ' + ('electron-poor ring' if sigma > 0.1 else 'electron-rich ring' if sigma < -0.1 else 'about neutral ring')
                         if named else 'ring: no other substituent')
        conj = any(b.GetIsConjugated() for b in a.GetBonds())
        if conj: parts.append('lone pair conjugated (delocalised: less available)' if a.GetSymbol() in ('N', 'O', 'S') else 'conjugated')
        if a.GetSymbol() == 'C': parts.append(f'oxidation state {oxidation(x):+g}'.replace('+0', '0'))
        return {'charge': q, 'sigma': sigma, 'conjugated': conj, 'text': '; '.join(parts)}

    def taft(x):                                        # the groups on the reacting atom, not counting its own ring or a =O/=S partner
        xr = [set(r) for r in m.GetRingInfo().AtomRings() if x in r]
        groups = []
        for j in heavy(x):
            if A(j).GetSymbol() in ('O', 'S') and bo(x, j) == 2: continue
            if any(j in r for r in xr): continue
            hit = None
            for name, sm, ss, es in TAFT:
                if any(mt[0] == j for mt in m.GetSubstructMatches(PK._q(sm), uniquify=False)): hit = (name, ss, es, False); break
            if not hit and A(j).GetSymbol() == 'C' and not A(j).GetIsAromatic() and all(_order(bd) == 1 for bd in A(j).GetBonds()) and len(heavy(j)) in TAFT_BRANCH:
                nm, es = TAFT_BRANCH[len(heavy(j))]; hit = (nm, None, es, True)
            if hit: groups.append({'atom': j, 'name': hit[0], 'sigma_star': hit[1], 'Es': hit[2], 'estimated': hit[3]})
        ss = [g['sigma_star'] for g in groups if g['sigma_star'] is not None]
        es = [g['Es'] for g in groups if g['Es'] is not None]
        s_sum = round(sum(ss), 3) if ss else None
        e_min = min(es) if es else None
        f = lambda v: f'{v:+.2f}'.replace('+0.00', '0.00')
        st = ('Taft Es: ' + ', '.join(f'{g["name"]} {f(g["Es"])}' for g in groups if g['Es'] is not None) + f'; bulkiest {f(e_min)}') if es else ''
        el = ('Taft σ*: ' + ', '.join(f'{g["name"]} {f(g["sigma_star"])}' for g in groups if g['sigma_star'] is not None) + f'; Σσ* {f(s_sum)} → '
              + ('its groups pull electrons away' if s_sum >= 0.5 else 'its groups push electrons in' if s_sum <= -0.1 else 'about like methyl groups')) if ss else ''
        return {'groups': groups, 'sigma_star': s_sum, 'Es_min': e_min, 'steric_text': st, 'electronic_text': el}

    out = {}
    for fid, q in frags:
        if q is None: continue
        for mt in m.GetSubstructMatches(q, uniquify=True, maxMatches=1000):
            x = pick(list(mt))
            if x is None or (fid, x) in out: continue
            t = taft(x); st_, el_ = steric(x), electronic(x)
            if t['steric_text']: st_['text'] += '; ' + t['steric_text']
            if t['electronic_text']: el_['text'] = (el_['text'] + '; ' if el_['text'] else '') + t['electronic_text']
            out[(fid, x)] = {'fragment': fid, 'atom': x, 'steric': st_, 'electronic': el_, 'taft': t}
    return sorted(out.values(), key=lambda e: (e['atom'], e['fragment']))
