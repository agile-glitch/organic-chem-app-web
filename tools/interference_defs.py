"""Interferences: everything that hinders (or helps) a template, written down as patterns the program can test.

An interference is a pattern that, when it sits in a particular relation to the atoms a template changes, makes that reaction faster, slower,
blocked, redirected or quenched. Two kinds are tested on every candidate product of a learned template (tools/int_marks.py gives every atom the
template changes the isotope 900 in the product, so a pattern can say "the reacting atom" as [900*]):

  SITE patterns  A SMARTS whose FIRST atom is the reacting atom (tools/interferences.py adds the mark). True when the surroundings of the reacting
                 atom look like that: an ortho nitro group on the ring carbon that reacts, a branch on the carbon that is attacked, an aryl group on
                 the nitrogen that is acylated ... [!900*] means "an atom the reaction does not change", i.e. a substituent the reacting atom had
                 before. The patterns are tested on the PRODUCT, so they describe what the reaction leaves in place around the atoms it changed (a
                 group the reaction consumed or created is not seen as a substituent). The weights (which pattern matters for which kind of reaction,
                 how much, and in which direction) are LEARNED from the recorded reactions in which a template could have reacted at several sites
                 (tools/int_contests.py, tools/int_train.py): the recorded product shows which site won.
  ADVISORIES     A reagent (a fragment, a role or a part as typed) and a group that is present: a Grignard and an alcohol, LiAlH4 and an ester, TFA
                 and a Boc group, Pd/H2 and a thioether ... They do not point at one atom: the reagent also does something to a group that is not
                 the target, or the group stops the reagent working. How often the group really reacts is MEASURED on the recorded reactions
                 (tools/int_evidence.py), and so is how often the advisory would be shown.

tools/add_interferences.py writes them, with the evidence, to the sheet "Interferences" of reaction_rules.xlsx; js/interference.js tests the same
patterns in the browser (tools/int_check_js.py checks that the two give the same answers).
"""
import re

U = '!900*'                      # an atom the reaction does not change (a substituent the reacting atom had before)
C_ = f'[#6;{U}]'                 # an unchanged carbon neighbour

# ---------------------------------------------------------------------------------------------------- reaction families
# The learned weight of a pattern is kept per family (and as a global value the families share). The family of a candidate product is read from the
# fragments it consumed (the Fragments sheet): the first family, in this order, that has a fragment among the consumed ones.
FAMILIES = [
    ('aryl_halide', 'substitution or coupling at an aryl / vinyl halide, sulfonate, diazonium or organometallic partner (SNAr, Sandmeyer, Pd / Cu / Ni coupling)',
     ['F013', 'F047', 'F084', 'F085', 'F086', 'F087', 'F065', 'F027', 'F170', 'F171', 'F172', 'F173', 'F175', 'F205', 'F081', 'F082', 'F083']),
    ('alkyl_substitution', 'substitution at an sp3 carbon that carries a halide, sulfonate or similar leaving group, or opening of an epoxide (SN2, SN1)',
     ['F010', 'F011', 'F012', 'F043', 'F044', 'F045', 'F046', 'F036', 'F174', 'F066', 'F125', 'F032', 'F033', 'F023']),
    ('deprotection', 'removal of a protecting group (Boc, Cbz, benzyl, PMB, allyl, trityl, Fmoc, silyl, acetal, phosphonate or sulfonyl) or cleavage of an ether',
     ['F077', 'F078', 'F130', 'F131', 'F132', 'F133', 'F134', 'F135', 'F052', 'F053', 'F039', 'F024', 'F050', 'F051', 'F103', 'F069']),
    ('acyl_substitution', 'acyl and sulfonyl substitution (acids, esters, amides, acyl chlorides, anhydrides, activated esters, sulfonyl chlorides, trifluoroacetyl)',
     ['F001', 'F004', 'F020', 'F021', 'F041', 'F058', 'F118', 'F119', 'F120', 'F121', 'F122', 'F123', 'F124', 'F136', 'F137', 'F138', 'F139', 'F184', 'F185', 'F186', 'F187', 'F188', 'F202',
      'F079', 'F074', 'F100', 'F080', 'F064', 'F088']),
    ('carbonyl_addition', 'addition to, or reduction of, an aldehyde, ketone, imine or nitrile; enolate chemistry (reductive amination, Grignard, Wittig, aldol)',
     ['F018', 'F019', 'F116', 'F117', 'F031', 'F049', 'F196', 'F204', 'F158', 'F022', 'F191', 'F101', 'F102', 'F192', 'F193', 'F194', 'F195']),
    ('conjugate_addition', 'addition to an α,β-unsaturated carbonyl (Michael, conjugate reduction)', ['F040', 'F189', 'F190']),
    ('unsaturated', 'addition to, or oxidation of, an alkene, alkyne or diene (also click chemistry: azide + alkyne)', ['F005', 'F006', 'F007', 'F008', 'F030', 'F037', 'F126', 'F127', 'F203']),
    ('redox_other', 'reduction or oxidation of a nitro, azide, nitroso, sulfur or similar group', ['F025', 'F056', 'F055', 'F057', 'F092', 'F072', 'F073', 'F071', 'F070']),
    ('n_nucleophile', 'a nitrogen nucleophile reacts with an electrophile the sheet has no fragment for (also N-oxidation and quaternisation of ring nitrogens)',
     ['F002', 'F034', 'F026', 'F093', 'F094', 'F095', 'F096', 'F097', 'F098', 'F099', 'F115', 'F089', 'F090', 'F109', 'F110', 'F111', 'F067', 'F068', 'F176', 'F177', 'F178', 'F179', 'F180', 'F181', 'F182', 'F183',
      'F091', 'F112', 'F113']),
    ('o_nucleophile', 'an oxygen nucleophile reacts, or an alcohol is oxidised or substituted',
     ['F003', 'F014', 'F015', 'F016', 'F017', 'F028', 'F129', 'F197', 'F198', 'F199']),
    ('arene_ch', 'reaction at an aromatic C–H (electrophilic substitution, metalation, C–H functionalisation); last in the list because a lost aryl group (benzyl, Cbz, trityl) also counts as a consumed C–H',
     ['F009', 'F104', 'F105', 'F106', 'F107', 'F108']),
]
FAMILY_NAMES = [f[0] for f in FAMILIES] + ['other']
FAMILY_TEXT = {f[0]: f[1] for f in FAMILIES}
FAMILY_TEXT['other'] = 'a reaction none of the families above describes'


def family_of(consumed):
    """the family of a candidate that consumed these fragment IDs (strings such as 'F013'): the first family in the list that has one of them"""
    for name, _, ids in FAMILIES:
        if any(i in consumed for i in ids): return name
    return 'other'


# ---------------------------------------------------------------------------------------------------- SITE patterns
# substituent on the ring: (label, key, SMARTS of the substituent with its first atom bonded to the ring, sigma_meta, sigma_para, class, why, example SMILES)
# Hammett constants: Hansch, Leo, Taft, Chem. Rev. 91, 165 (1991)
SUBSTITUENTS = [
    ('nitro', 'nitro', '[N+](=O)[O-]', 0.71, 0.78, 'EWG', 'the strongest common electron-withdrawing group (induction and resonance)', '[N+](=O)[O-]'),
    ('cyano', 'cyano', 'C#N', 0.56, 0.66, 'EWG', 'withdraws by induction and by resonance', 'C#N'),
    ('trifluoromethyl', 'cf3', 'C(F)(F)F', 0.43, 0.54, 'EWG', 'withdraws by induction (three C–F dipoles); no resonance', 'C(F)(F)F'),
    ('sulfonyl (sulfone, sulfonamide)', 'sulfonyl', '[SX4](=O)=O', 0.60, 0.72, 'EWG', 'withdraws by induction and by resonance', 'S(C)(=O)=O'),
    ('formyl', 'formyl', '[CH1]=O', 0.35, 0.42, 'EWG', 'the C=O pulls electron density out of the ring, most at ortho and para', 'C=O'),
    ('acyl (ketone)', 'acyl', '[CX3](=O)[#6]', 0.38, 0.50, 'EWG', 'the C=O pulls electron density out of the ring, most at ortho and para', 'C(C)=O'),
    ('carboxy (acid, ester)', 'carboxy', '[CX3](=O)[OX2]', 0.37, 0.45, 'EWG', 'an acid or ester C=O withdraws electrons from the ring', 'C(=O)OC'),
    ('carboxamide', 'carboxamide', '[CX3](=O)[NX3]', 0.28, 0.36, 'EWG', 'the amide C=O withdraws, less than an ester does', 'C(=O)NC'),
    ('trifluoromethoxy', 'ocf3', '[OX2]C(F)(F)F', 0.38, 0.35, 'EWG', 'withdraws by induction; the O lone pair is held back by the CF3', 'OC(F)(F)F'),
    ('sulfonyloxy (OMs, OTs)', 'osulfonyl', '[OX2][SX4](=O)=O', 0.39, 0.36, 'EWG', 'withdraws by induction; the O lone pair is held back by the sulfonyl group (triflate is stronger: σp 0.53)', 'OS(C)(=O)=O'),
    ('fluoro', 'fluoro', 'F', 0.34, 0.06, 'halogen', 'strongly withdrawing by induction but a lone-pair donor: directs ortho/para in electrophilic substitution', 'F'),
    ('chloro', 'chloro', 'Cl', 0.37, 0.23, 'halogen', 'withdrawing by induction, a weak lone-pair donor: deactivates but directs ortho/para', 'Cl'),
    ('bromo', 'bromo', 'Br', 0.39, 0.23, 'halogen', 'withdrawing by induction, a weak lone-pair donor: deactivates but directs ortho/para', 'Br'),
    ('iodo', 'iodo', 'I', 0.35, 0.18, 'halogen', 'withdrawing by induction, a weak lone-pair donor: deactivates but directs ortho/para', 'I'),
    ('thioether', 'sr', '[SX2;!$(S-C(F)(F)F);!$(S-C#N)][#6]', 0.15, 0.00, 'EDG', 'a weak donor by resonance, a weak withdrawer by induction', 'SC'),
    ('acylamino (amide N–H)', 'nhcor', '[NX3;H1]C(=O)', 0.21, 0.00, 'EDG', 'the N lone pair is shared with the C=O: a moderate donor', 'NC(C)=O'),
    ('sulfonylamino (N–H sulfonamide)', 'nhso2', '[NX3;H1][SX4](=O)=O', 0.20, 0.03, 'EDG', 'the N lone pair is shared with the sulfonyl group: about neutral, a weak ortho/para director', 'NS(C)(=O)=O'),
    ('dialkylamino', 'nr2', '[NX3;H0;+0;!$(N-C=O);!$(N-S=O);!$(N(-c)-c)]([#6])[#6]', -0.16, -0.83, 'EDG', 'the N lone pair donates strongly into the ring (resonance)', 'N(C)C'),
    ('amino (NH2, NHR)', 'nh', '[NX3;H2,H1;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-[#7,#8])]', -0.16, -0.66, 'EDG', 'the N lone pair donates strongly into the ring (resonance)', '[NH2]'),
    ('hydroxy', 'oh', '[OX2H1]', 0.12, -0.37, 'EDG', 'the O lone pair donates into the ring (resonance)', '[OH]'),
    ('alkoxy', 'or', '[OX2][CX4;!$(C(F)(F)F)]', 0.12, -0.27, 'EDG', 'the O lone pair donates into the ring (resonance)', 'OC'),
    ('aryloxy', 'oar', '[OX2]c', 0.25, -0.03, 'EDG', 'the O lone pair is shared with two rings: a weak donor', 'Oc1ccccc1'),
    ('acyloxy', 'ocor', '[OX2]C(=O)', 0.39, 0.31, 'mixed', 'the O lone pair is shared with the C=O: it withdraws by induction (σp is positive) yet donates to a cation by resonance (σp⁺ −0.19)', 'OC(C)=O'),
    ('methyl', 'methyl', '[CH3]', -0.07, -0.17, 'EDG', 'a weak donor (induction, hyperconjugation)', 'C'),
    ('alkyl (not methyl)', 'alkyl', '[CX4;!H3;!$(C-[F,Cl,Br,I,$(C#N)])]', -0.07, -0.15, 'EDG', 'a weak donor (induction, hyperconjugation); bulkier than methyl', 'CC'),
    ('aryl', 'aryl', 'c', 0.06, -0.01, 'EDG', 'a second ring: a weak donor by resonance', 'c1ccccc1'),
    ('vinyl', 'vinyl', '[CX3]=[CX3]', 0.06, -0.04, 'EDG', 'conjugated C=C: a weak donor', 'C=C'),
    ('ethynyl', 'ethynyl', 'C#C', 0.21, 0.23, 'EWG', 'an sp carbon: weakly electron-withdrawing', 'C#C'),
    ('silyl', 'silyl', '[SiX4]', -0.04, -0.07, 'EDG', 'a weak donor (σ-donation, β-effect)', '[Si](C)(C)C'),
    ('ammonium', 'nr3plus', '[NX4+]', 0.88, 0.82, 'EWG', 'a positive charge next to the ring withdraws strongly', '[N+](C)(C)C'),
    ('phenoxide (O⁻)', 'ominus', '[OX1-]', -0.47, -0.81, 'EDG', 'an anionic O donates more strongly than OH', '[O-]'),
    ('carboxylate', 'co2minus', '[CX3](=O)[OX1-]', -0.10, 0.00, 'neutral', 'an anionic carboxylate: neutral to weakly donating', 'C(=O)[O-]'),
]
RELATIONS = [('o', 'ortho'), ('m', 'meta'), ('p', 'para')]
SUBST_KEYS = [s[1] for s in SUBSTITUENTS]
SIGMA = {s[1]: {'m': s[3], 'p': s[4], 'o': s[4]} for s in SUBSTITUENTS}       # σ used for ortho = σp (ortho effects are electronic like para, plus steric, which has its own patterns)

SITE = []                        # dicts: key, name, group, kind, anchor, smarts, relation, effect, why, source, sigma (substituent flags), prefilter


def site(key, name, group, kind, anchor, smarts, relation, effect, why, source='', sigma=None, prefilter=None):
    SITE.append(dict(key=key, name=name, group=group, kind=kind, anchor=anchor, smarts=smarts, relation=relation, effect=effect, why=why, source=source, sigma=sigma, prefilter=prefilter))


def guard(sub):
    """the substituent SMARTS with '!900*' added to its first atom (an atom the reaction did not change). It is written INLINE in the path patterns: inside a
    recursive $() the substituent would be matched on its own, and its first neighbour could then be the ring atom itself (a carboxylate would pass as a ketone)"""
    if sub.startswith('['):
        depth = 0
        for i, ch in enumerate(sub):
            if ch == '[': depth += 1
            elif ch == ']':
                depth -= 1
                if depth == 0: return sub[:i] + ';!900*' + sub[i:]
    tok = re.match(r'Cl|Br', sub)
    tok = tok.group(0) if tok else sub[0]
    return '[' + tok + ';!900*]' + sub[len(tok):]


def ring_patterns(sub):
    """the ortho / meta / para patterns of a ring substituent (SMARTS of the group, first atom bonded to the ring): the substituent must be an atom the reaction did not
    change; meta and para are measured inside ONE six-membered ring (a peri or ring-fused position is not meta or para), the intermediate atoms may be ring heteroatoms
    (the 6-position of a 2-pyridyl carbon is meta, through the N), and five-membered rings have their own paths (1,3 = meta-like, 1,4 = para-like around the ring)"""
    g = guard(sub)
    return {'o': f'[c;$(c:c-{g})]',
            'm': f'[c;$(c1:a:c(-{g}):a:a:a:1),$([c;r5]:[c;r5]:[c;r5]-{g})]',
            'p': f'[c;$(c1:a:a:c(-{g}):a:a:1),$([c;r5]:[c;r5]:[c;r5]:[c;r5]-{g})]'}


_EFFECT = {'EWG': 'withdraws electrons: speeds nucleophilic attack at the ring carbon (SNAr, oxidative addition, metalation next to it), slows electrophilic substitution',
           'EDG': 'donates electrons: speeds electrophilic substitution at the ring carbon, slows nucleophilic attack (SNAr, oxidative addition)',
           'halogen': 'withdraws by induction but donates by resonance: slows electrophilic substitution yet directs it ortho/para; speeds SNAr a little',
           'neutral': 'has almost no net electronic effect on the ring carbon (σ close to 0); a carboxylate is a weak meta director in electrophilic substitution',
           'mixed': 'withdraws by induction (σp > 0: speeds nucleophilic attack a little) yet donates to a cation by resonance (σp⁺ < 0): a moderate activator of electrophilic substitution that directs ortho/para'}
_POS_NOTE = {'o': ' At ortho the electronic effect is like para, and the group also crowds the carbon.',
             'm': ' At meta mostly the inductive part acts; resonance does not reach this position.',
             'p': ' At para induction and resonance both act.'}
_HLT = 'Hansch, Leo, Taft, Chem. Rev. 91, 165 (1991); Clayden, Organic Chemistry (aromatic substitution)'
for label, key, sub, sm, sp, kind, why, _ex in SUBSTITUENTS:
    pats = ring_patterns(sub)
    for tag, rel in RELATIONS:
        sig = sm if tag == 'm' else sp
        where = {'o': 'on the ring atom next to the reacting carbon (ortho)', 'm': 'on the ring atom two bonds from the reacting carbon, inside the same six-membered ring (meta)',
                 'p': 'on the ring atom opposite the reacting carbon, inside the same six-membered ring (para)'}[tag]
        site(f'ring-{rel}-{key}', f'{label}, {rel} to the reacting aromatic carbon', 'ring substituent (position)', 'electronic', 'aromatic C', pats[tag], where,
             _EFFECT[kind] + '.' + _POS_NOTE[tag] + (' A ring at ortho also twists out of the plane; with three or more ortho groups the biaryl axis is locked (atropisomers).' if key == 'aryl' and tag == 'o' else ''),
             f'{label} ({"σm" if tag == "m" else "σp"} {sig:+.2f}): {why}', _HLT, sigma=sig, prefilter=f'[c]-{sub}')

_JM = 'Joule & Mills, Heterocyclic Chemistry'
site('ring-alpha-N', 'reacting carbon next to a pyridine-type ring nitrogen (α)', 'ring type', 'electronic', 'aromatic C', '[c;$(c:[nX2])]', 'the reacting carbon is bonded to a pyridine-type ring N (2-position of pyridine, C2 of azoles)',
     'electron-poor carbon: nucleophiles add here (SNAr, Chichibabin, Minisci) and the C–X bond is activated for oxidative addition', 'the ring N pulls electron density from the neighbouring carbon and takes the negative charge of the intermediate', _JM)
site('ring-beta-N', 'reacting carbon one atom from a pyridine-type ring nitrogen (β)', 'ring type', 'electronic', 'aromatic C', '[c;$(c1:c:[nX2]:a:a:a:1),$([c;r5]:[c;r5]:[nX2;r5])]', 'the reacting carbon is two bonds from a pyridine-type ring N, in the same ring (3-position of pyridine)',
     'the least electron-poor position of an azine: electrophilic substitution goes here (when it goes at all)', 'no resonance structure puts the charge on N, so the β position behaves like a deactivated benzene carbon', _JM)
site('ring-gamma-N', 'reacting carbon opposite a pyridine-type ring nitrogen (γ)', 'ring type', 'electronic', 'aromatic C', '[c;$(c1:a:a:[nX2]:a:a:1)]', 'the reacting carbon is para to a pyridine-type ring N in the same six-membered ring (4-position of pyridine)',
     'electron-poor carbon: nucleophiles add here (SNAr, Minisci) and oxidative addition is easy', 'the ring N takes the negative charge of the intermediate by resonance', _JM)
site('ring-between-2N', 'reacting carbon between two ring nitrogens', 'ring type', 'electronic', 'aromatic C', '[c;$(c(:n):n)]', 'C2 of pyrimidine, imidazole, triazole: a ring N on each side',
     'the most electron-poor carbon of an azole (C2 of imidazole, thiazole, triazole): its C–H is the acidic one and lithiates or couples first; in pyrimidines, quinazolines and purines SNAr and oxidative addition go to C4 / C6 before C2 (2,4-dichloropyrimidine reacts first at C4)', 'two ring N atoms withdraw electrons and stabilise the anion', _JM)
site('ring-alpha-azinium', 'reacting carbon next to a ring N⁺ (pyridinium, N-oxide)', 'ring type', 'electronic', 'aromatic C', '[c;$(c:[n+])]', 'the reacting carbon is bonded to a positively charged ring N (C2 of an N-oxide or N-alkylpyridinium)',
     'very electron-poor: nucleophiles add here (POCl3 / Ac2O chlorination and acyloxylation of N-oxides, Minisci, Reissert), metalation and Pd C–H activation are directed here', 'the positive charge withdraws strongly and takes the incoming nucleophile', _JM)
site('ring-gamma-azinium', 'reacting carbon opposite a ring N⁺ (γ of pyridinium, N-oxide)', 'ring type', 'electronic', 'aromatic C', '[c;$(c1:a:a:[n+]:a:a:1)]', 'the reacting carbon is para to a positively charged ring N',
     'electron-poor: nucleophiles add here (4-substitution of pyridine N-oxides and pyridiniums), nitration of N-oxides goes here', 'the positive charge withdraws and resonance carries it to this position', _JM)
site('ring-alpha-S', 'reacting carbon next to a ring sulfur (α)', 'ring type', 'electronic', 'aromatic C', '[c;$(c:s)]', 'a carbon bonded to a ring S: C2/C5 of thiophene, C2/C5 of thiazole',
     'thiophene: electrophilic substitution and lithiation go to the α carbon first; thiazole: C2 (between S and N) is the acidic, lithiated and nucleophile-attacked carbon, C5 is the electrophilic-substitution site', 'the S lone pair stabilises the cation formed by attack at the α carbon better than at β', _JM)
site('ring-alpha-O', 'reacting carbon next to a ring oxygen (α)', 'ring type', 'electronic', 'aromatic C', '[c;$(c:o)]', 'a carbon bonded to a ring O: C2/C5 of furan, C2 of oxazole',
     'furan: electrophilic substitution goes to the α carbon first (furan is more reactive than thiophene); oxazole: C2 is the acidic, lithiated carbon, C5 the electrophilic-substitution site', 'the O lone pair stabilises the cation formed by attack at the α carbon', _JM)
site('ring-beta-SO', 'reacting carbon one atom from a ring sulfur or oxygen (β, five-membered ring)', 'ring type', 'electronic', 'aromatic C', '[c;r5;$(c:c:[s,o])]', 'C3 of thiophene or furan (a five-membered ring carbon two bonds from S or O)',
     'the less reactive position of the ring toward electrophiles; reached when both α positions are taken', 'the cation formed by β attack has fewer resonance forms', _JM)
site('ring-alpha-NR', 'reacting carbon next to a pyrrole-type ring nitrogen (N–H or N–R)', 'ring type', 'electronic', 'aromatic C', '[c;r5;$(c:[nX3;r5;+0])]', 'a five-membered-ring carbon bonded to a neutral pyrrole-type N: C2/C5 of pyrrole, C2 of indole, C2/C5 of N-substituted imidazole',
     'pyrrole: electrophilic substitution goes to C2 first; indole: C2 is the carbon that lithiates and couples (C–H activation) while electrophiles go to C3; the N substituent directs metalation to C2', 'the N lone pair is part of the aromatic sextet and donates strongly to the α carbon', _JM)
site('ring-beta-NR', 'reacting carbon one atom from a pyrrole-type ring nitrogen (β)', 'ring type', 'electronic', 'aromatic C', '[c;r5;$([c;r5]:[c;r5]:[nX3;r5;+0])]', 'C3 of indole or pyrrole',
     'indole C3 is the most nucleophilic carbon (Vilsmeier, Mannich, halogenation, Friedel–Crafts); pyrrole C3 reacts when C2 is blocked', 'attack at C3 of indole keeps the benzene ring aromatic in the cation', _JM)
site('ring-5-membered', 'reacting carbon in an electron-rich five-membered ring (furan, thiophene, pyrrole and their benzo-fused forms)', 'ring type', 'electronic', 'aromatic C',
     '[c;r5;$(c1:c:c:c:[o,s,nX3;+0]:1),$(c1:c:c:[o,s,nX3;+0]:c:1)]', 'the reacting carbon is in a five-membered ring made of one O, S or pyrrole-type N and four carbons',
     'electron-rich ring (six electrons on five atoms): electrophilic substitution is far faster than for benzene (pyrrole ≫ furan > thiophene)', 'π-excessive heteroaromatics', _JM)
site('ring-fused-alpha', 'reacting carbon next to a ring-fusion atom (α of naphthalene, quinoline …)', 'ring type', 'electronic', 'aromatic C', '[c;r6;$(c:[c;R2;$(c(:a)(:a):a)])]', 'a six-membered-ring carbon bonded to a ring-fusion carbon whose other neighbours are aromatic too (C1 of naphthalene, C5/C8 of quinoline, C4/C7 of indole; not the benzo carbons of benzodioxole, tetralin or indoline)',
     'naphthalene: the α position is the kinetic site of electrophilic substitution (more hindered, peri, for large reagents); in indole, benzimidazole and benzofuran the benzo carbons C4/C7 are less reactive than the hetero ring', 'more resonance forms of the naphthalene intermediate keep the other ring aromatic', 'Clayden, Organic Chemistry')
site('ring-peri', 'reacting carbon with a substituent in the peri position', 'ring sterics', 'steric', 'aromatic C', '[c;$(c:[c;R2]:[c;$(c-!@[!#1;!900*])])]', 'a substituent sits on the ring atom of the neighbouring ring that is next to the fusion carbon (C8 for C1 of naphthalene, C4 for C3 of indole)',
     'peri strain: attack and coupling at this carbon are slowed (1,8-disubstituted naphthalenes, 3,4-disubstituted indoles), large reagents go elsewhere', 'the peri substituent lies in the plane of the ring, about 2.5 Å from the reacting carbon', 'Clayden, Organic Chemistry')

# sterics around the reacting aromatic carbon
site('ring-ortho-1', 'at least one ortho substituent', 'ring sterics', 'steric', 'aromatic C', '[c;$(c:[c;$(c-[!#1;!900*])])]', 'a ring neighbour of the reacting carbon carries a substituent (also a ring-fused alkyl group)',
     'crowding next to the reacting carbon: slows attack and coupling, favours the less hindered position', 'the new bond has to form in the plane of the ring beside a substituent', 'Clayden, Organic Chemistry')
site('ring-ortho-2', 'two ortho substituents (2,6-disubstituted)', 'ring sterics', 'steric', 'aromatic C', '[c;$(c(:[c;$(c-[!#1;!900*])]):[c;$(c-[!#1;!900*])])]', 'both ring neighbours of the reacting carbon carry a substituent',
     'strongly hindered: couplings and substitutions need forcing conditions or special ligands; other positions win when they exist', 'the reacting carbon is flanked on both sides', 'Clayden, Organic Chemistry')
site('ring-ortho-branched', 'a branched ortho substituent', 'ring sterics', 'steric', 'aromatic C', '[c;$(c:c-[CX4;D3,D4;!900*])]', 'an ortho substituent that is itself branched (isopropyl, tert-butyl, CHR2)',
     'much bulkier than an ortho methyl: blocks the position for most reagents', 'the branched group sits over the reacting carbon', 'Taft steric constants Es: i-Pr −0.47, t-Bu −1.54')

# an aromatic nitrogen is the reacting atom (N-alkylation, N-arylation, N-acylation of azoles and azines)
site('azole-n-hindered', 'azole nitrogen with a substituent on a neighbouring carbon', 'ring nitrogen', 'steric', 'aromatic n', '[n;r5;$(n:c-[!#1;!900*])]', 'the reacting five-membered-ring N is next to a ring carbon that carries a substituent',
     'slower alkylation / arylation: the less hindered ring N of a pyrazole, imidazole or triazole reacts', 'the substituent crowds the lone pair or the new N–R bond', _JM)
site('azole-n-hindered-2', 'azole nitrogen with substituents on both neighbouring carbons', 'ring nitrogen', 'steric', 'aromatic n', '[n;r5;$(n(:c-[!#1;!900*]):c-[!#1;!900*])]', 'both ring carbons next to the reacting N carry substituents',
     'strongly hindered: this ring N is the last to react', 'both sides of the lone pair are crowded', _JM)
site('azole-n-next-to-N', 'azole nitrogen next to another ring nitrogen', 'ring nitrogen', 'electronic', 'aromatic n', '[n;r5;$(n:n)]', 'the reacting N is bonded to another ring N (pyrazole, indazole, triazole, tetrazole)',
     'N-alkylation gives N1 / N2 regioisomer mixtures: the product ratio follows sterics and the tautomer that reacts, so the less hindered N usually wins', 'the anion is delocalised over both N atoms', _JM)
site('azole-n-1,3', 'azole nitrogen separated from another ring nitrogen by one carbon', 'ring nitrogen', 'electronic', 'aromatic n', '[n;r5;$(n:c:n)]', 'imidazole, benzimidazole, purine: N–C–N',
     'the basic pyridine-type N is the nucleophile of neutral imidazoles; the anion alkylates at either N', 'the amidine-like N–C=N unit', _JM)
site('azole-n-fused', 'azole nitrogen next to a ring-fusion carbon', 'ring nitrogen', 'steric', 'aromatic n', '[n;r5;$(n:[c;R2])]', 'the reacting N is bonded to a carbon shared with a second ring (indole N1, benzimidazole, indazole N1)',
     'more hindered by the peri C–H of the fused ring, but conjugation keeps the N–R product stable', 'the fused ring crowds one side', _JM)
site('azine-n', 'six-membered ring nitrogen', 'ring nitrogen', 'electronic', 'aromatic n', '[n;r6]', 'the reacting N is the N of pyridine, pyrimidine or pyridone',
     'pyridine N: quaternises with alkyl halides, forms N-oxides, binds metals; pyridone N–H vs O: N-alkylation with soft electrophiles, O-alkylation with hard ones or Ag salts', 'the lone pair lies in the ring plane, not in the π system', _JM)
site('azine-n-hindered', 'six-membered ring nitrogen with a substituent on a neighbouring carbon', 'ring nitrogen', 'steric', 'aromatic n', '[n;r6;$(n:c-[!#1;!900*])]', 'the reacting ring N is next to a carbon that carries a substituent (2-substituted pyridine)',
     'slower quaternisation, N-oxidation and metal binding (2-picoline is less nucleophilic at N than pyridine)', 'the substituent crowds the lone pair', _JM)
site('azine-n-hindered-2', 'six-membered ring nitrogen with substituents on both neighbouring carbons', 'ring nitrogen', 'steric', 'aromatic n', '[n;r6;$(n(:c-[!#1;!900*]):c-[!#1;!900*])]', 'both carbons next to the reacting ring N carry substituents (2,6-disubstituted pyridine)',
     'a hindered base that does not alkylate or bind metals (2,6-lutidine, 2,6-di-tert-butylpyridine)', 'both sides of the lone pair are crowded', _JM)

# ---------------------------------------------------------------------------------------------------- sp3 carbon is the reacting atom
_CL = 'Clayden, Organic Chemistry'
site('sp3-methyl', 'methyl-type carbon (no carbon substituent)', 'sp3 carbon: substitution degree', 'steric', 'sp3 C', f'[C;X4;!$(C-{C_})]', 'the attacked carbon had no carbon substituent (CH3–X); only unchanged carbons are counted, the sp3-<degree>-at-<atom> rows below count all carbon neighbours',
     'the fastest SN2 and the least hindered position; no cation possible', 'no substituent blocks the backside; no alkyl group stabilises a cation', _CL + ' (SN2)')
site('sp3-primary', 'primary carbon (one carbon substituent)', 'sp3 carbon: substitution degree', 'steric', 'sp3 C', f'[C;X4;$(C-{C_});!$(C(-{C_})-{C_})]', 'the attacked carbon had one carbon substituent (R–CH2–X)',
     'open: SN2 fast, E2 with bulky bases', 'one substituent only: the backside is open', _CL + ' (SN2)')
site('sp3-secondary', 'secondary carbon (two carbon substituents)', 'sp3 carbon: substitution degree', 'steric', 'sp3 C', f'[C;X4;$(C(-{C_})-{C_});!$(C(-{C_})(-{C_})-{C_})]', 'the attacked carbon had two carbon substituents (R2CH–X)',
     'SN2 slower; SN1, E1 and E2 compete', 'two substituents crowd the backside and stabilise a cation', _CL + ' (SN2 / SN1)')
site('sp3-tertiary', 'tertiary carbon (three carbon substituents)', 'sp3 carbon: substitution degree', 'steric', 'sp3 C', f'[C;X4;$(C(-{C_})(-{C_})-{C_})]', 'the attacked carbon had three carbon substituents (R3C–X)',
     'no SN2: SN1 / E1 through the cation, or E2', 'the backside is shielded; three alkyl groups stabilise the cation', _CL + ' (SN1)')
site('sp3-benzylic', 'benzylic carbon', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[c;!900*])]', 'the reacting carbon is attached to an aromatic ring',
     'SN1 and SN2 both fast; radical and oxidation chemistry at this C–H; easily inserted by metals', 'the ring stabilises the cation, the radical and the SN2 transition state', _CL)
site('sp3-allylic', 'allylic carbon', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[C;X3;!$(C=O);!900*]=C)]', 'the reacting carbon is attached to a C=C',
     'SN1 and SN2 fast; SN2′ and allylic transposition possible; radical reactions at this C–H', 'the π bond stabilises the cation, the radical and the transition state', _CL)
site('sp3-propargylic', 'propargylic carbon', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[C;X2;!900*]#[#6])]', 'the reacting carbon is attached to a C≡C',
     'SN2 fast (like allylic); SN2′ gives allenes', 'the π system activates the transition state', _CL)
site('sp3-alpha-carbonyl', 'carbon α to a carbonyl group', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[C;X3;!900*]=[O,N,S])]', 'the reacting carbon is next to a C=O or C=N',
     'SN2 much faster (α-halo ketones are among the most reactive alkylating agents); the C–H is acidic: enolate chemistry', 'the π* of the C=O lowers the SN2 transition state; the C–H is acidic', _CL)
site('sp3-alpha-nitrile', 'carbon α to a nitrile', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[C;X2;!900*]#N)]', 'the reacting carbon is next to a C≡N',
     'SN2 fast; the C–H is acidic (acetonitrile pKa 31 in DMSO)', 'the nitrile lowers the SN2 transition state and acidifies the C–H', _CL)
site('sp3-alpha-CF3', 'carbon next to a CF3 group', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[C;!900*](F)(F)F)]', 'the reacting carbon carries a CF3 group it already had (a CF3 the reaction put there does not count)',
     'SN2 and SN1 slowed: the CF3 destabilises positive charge and the transition state', 'the CF3 group withdraws strongly by induction (Taft σ* +2.61)', 'Hansch, Leo, Taft; Taft σ* values')
site('sp3-alpha-halogen', 'a halogen stays on the reacting carbon', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[F,Cl,Br,I;!900*])]', 'the reacting carbon keeps a halogen (CHX2, CX3)',
     'SN2 slowed (inductive, steric); SN1 favoured by lone-pair donation (Cl, Br)', 'a halogen both withdraws and donates a lone pair', _CL)
site('sp3-alpha-oxygen', 'an oxygen stays on the reacting carbon', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[OX2;!900*])]', 'the reacting carbon keeps an OR group (acetal carbon, α-alkoxy halide)',
     'SN1 very fast (oxocarbenium ion): α-alkoxy halides such as MOMCl are strongly activated; SN2 at an acetal carbon is unusual', 'the O lone pair stabilises the cation', _CL)
site('sp3-alpha-nitrogen', 'a nitrogen stays on the reacting carbon', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[NX3;+0;!900*])]', 'the reacting carbon keeps an amino group (aminal, α-amino)',
     'SN1 very fast (iminium ion); the α-amino C–H is easily oxidised', 'the N lone pair stabilises the cation', _CL)
site('sp3-neopentyl', 'neopentyl-type carbon (β-quaternary)', 'sp3 carbon: neighbour effect', 'steric', 'sp3 C', '[C;X4;$(C-[C;X4;D4;!900*;!$(C(F)(F)F)])]', 'the carbon next to the reacting one is quaternary (CH2–CMe3)',
     'SN2 essentially blocked (about 10^5 slower than ethyl); SN1 gives rearrangement', 'the tert-butyl group shields the backside of the CH2', _CL + ' (SN2); Streitwieser relative rates; Taft Es −1.74')
site('sp3-beta-branched', 'β-branched carbon (isobutyl-type)', 'sp3 carbon: neighbour effect', 'steric', 'sp3 C', '[C;X4;$(C-[C;X4;H1;!900*]([#6])[#6])]', 'the carbon next to the reacting one carries two further carbons (CH2–CHR2)',
     'SN2 about 30-fold slower than for an unbranched chain; E2 gets competitive', 'the branch crowds the approach of the nucleophile', _CL + ' (SN2); Streitwieser relative rates; Taft Es −0.93')
site('sp3-beta-opened', 'carbon next to a carbon that now bears the O, N or S the reaction released from it (ring opening)', 'sp3 carbon: neighbour effect', 'electronic', 'sp3 C', '[C;X4;$(C-[CX4;!900*]-[900O,900N,900S]);!$(C1-[CX4;!900*]-[900O,900N,900S]-1)]', 'the neighbouring carbon keeps a heteroatom whose bond to the reacting carbon was broken (epoxide, aziridine, cyclic sulfate opened at this carbon; a three-membered ring closed onto this carbon does not count)',
     'regioselectivity of ring opening: nucleophiles attack the less hindered carbon of an epoxide under basic or neutral conditions (SN2), the more substituted carbon under acid; the β-hydroxy / β-amino group left behind can chelate', 'the strained ring is opened at the carbon the pattern marks', _CL + ' (epoxide opening)')
for n, label, extra in ((3, 'three-membered ring (epoxide, aziridine, cyclopropane)', 'a strained ring was formed or kept at this carbon: three-membered rings form by intramolecular SN2 (halohydrin to epoxide), carbene or peracid addition, and open easily with nucleophiles or acid'),
                        (4, 'four-membered ring (oxetane, azetidine, cyclobutane)', 'a strained ring was formed or kept at this carbon; four-membered rings form slowly and open easily under acid'),
                        (5, 'five-membered ring', 'a five-membered ring carbon: five-membered rings close easily (5-exo); SN2 at cyclopentyl is not slowed by ring strain the way it is at cyclohexyl'),
                        (6, 'six-membered ring', 'a six-membered ring carbon: SN2 at cyclohexyl is much slower than at cyclopentyl or at an open-chain secondary carbon; axial or equatorial position decides E2 (anti-periplanar H) and the face of attack')):
    site(f'sp3-ring{n}', f'carbon in a {label}', 'sp3 carbon: ring', 'geometry', 'sp3 C', f'[C;X4;r{n}]', f'the reacting carbon is in a ring of {n} atoms in the product', extra, 'ring size sets strain and the geometry of the attack', _CL + ' (ring strain)')

# the carbon's degree together with the kind of atom that took its place (the changed neighbour): tells apart the two ends of an addition across a C=C or C=O.
# With a heteroatom or main-group partner ALL carbon neighbours count (the true degree of the carbon in the product); with a carbon partner only the unchanged ones
# (the partner carbon is itself one of the carbon neighbours)
_PART = {'O': '[900O]', 'N': '[900N]', 'S': '[900S]', 'X': '[900F,900Cl,900Br,900I]', 'M': '[900B,900Si,900P]', 'C': '[900C,900c]'}
_PART_TEXT = {'O': 'an oxygen', 'N': 'a nitrogen', 'S': 'a sulfur', 'X': 'a halogen', 'M': 'a boron, silicon or phosphorus', 'C': 'a carbon'}
_DEG_TEXT = {'methyl': 'no carbon substituent', 'primary': 'one carbon substituent', 'secondary': 'two carbon substituents', 'tertiary': 'three carbon substituents'}


def _deg_head(deg, nb):
    return {'methyl': f'[C;X4;!$(C-{nb})', 'primary': f'[C;X4;$(C-{nb});!$(C(-{nb})-{nb})',
            'secondary': f'[C;X4;$(C(-{nb})-{nb});!$(C(-{nb})(-{nb})-{nb})', 'tertiary': f'[C;X4;$(C(-{nb})(-{nb})-{nb})'}[deg]


for deg in ('methyl', 'primary', 'secondary', 'tertiary'):
    for p, ps in _PART.items():
        nb = C_ if p == 'C' else '[#6]'
        site(f'sp3-{deg}-at-{p}', f'{deg} carbon that now carries {_PART_TEXT[p]}', 'sp3 carbon: where the new group went', 'steric', 'sp3 C', f'{_deg_head(deg, nb)};$(C-{ps})]',
             f'the reacting carbon ({_DEG_TEXT[deg]}{" besides the new bond" if p == "C" else ""}) is bonded to {_PART_TEXT[p]} the reaction put there',
             'regioselectivity of additions and ring openings: the new group goes to the more substituted carbon in cationic / SN1-like reactions (Markovnikov: hydration, hydrohalogenation, halohydrin, acid-catalysed epoxide opening), '
             'to the less substituted carbon in SN2-like, radical and hydroboration reactions (anti-Markovnikov: base-catalysed epoxide opening, hydroboration, radical HBr)',
             'the two ends of the old double bond (or ring) differ in how well they carry positive charge or how open they are to attack', _CL + ' (alkene additions, epoxide opening)')

# ---------------------------------------------------------------------------------------------------- carbonyl carbon is the reacting atom (the C=O is still there in the product)
site('co-aryl', 'aryl on the carbonyl carbon', 'carbonyl carbon', 'electronic', 'carbonyl C', '[C;X3;$(C=[O,N,S]);$(C-[c;!900*])]', 'the retained carbonyl carbon (acyl substitution) is bonded to an aromatic ring',
     'less electrophilic than an alkyl carbonyl (conjugation) and more hindered; ortho substituents twist it out of conjugation', 'the ring donates into the C=O by resonance', _CL)
site('co-vinyl', 'C=C on the carbonyl carbon (enone, acrylate)', 'carbonyl carbon', 'electronic', 'carbonyl C', '[C;X3;$(C=[O,N,S]);$(C-[C;X3;!$(C=O);!900*]=C)]', 'the retained carbonyl carbon is bonded to a C=C',
     'hard nucleophiles add 1,2; soft nucleophiles add 1,4 at the β carbon instead', 'conjugation spreads the electrophilicity to the β carbon', _CL)
site('co-alpha-branched', 'α-branched carbonyl (CHR2 next to C=O)', 'carbonyl carbon', 'steric', 'carbonyl C', '[C;X3;$(C=[O,N,S]);$(C-[CX4;H1;!900*](-[#6])-[#6])]', 'the retained carbonyl carbon is bonded to a CH that carries two more carbons',
     'additions and acyl substitutions slower; the α C–H is more easily epimerised', 'the branch crowds the approach to the carbonyl', _CL)
site('co-alpha-quaternary', 'α-quaternary carbonyl (CR3 next to C=O)', 'carbonyl carbon', 'steric', 'carbonyl C', '[C;X3;$(C=[O,N,S]);$(C-[CX4;H0;!900*](-[#6])-[#6])]', 'the retained carbonyl carbon is bonded to a carbon with no H and at least two further carbons (pivaloyl, Aib, α-hydroxyisobutyryl)',
     'strongly hindered: additions and hydrolysis are slow (pivalates survive conditions that cleave acetates); no α C–H, so no enolisation', 'the tert-alkyl group shields the carbonyl carbon', _CL + '; Taft Es t-Bu −1.54')
site('co-alpha-hetero', 'α-heteroatom (halogen, O, N) on the carbon next to the carbonyl', 'carbonyl carbon', 'electronic', 'carbonyl C', '[C;X3;$(C=[O,N,S]);$(C-[CX4;$(C-[F,Cl,Br,I,O,N]);!$(C(F)(F)F);!900*])]', 'the carbon next to the retained carbonyl carbon carries a halogen, O or N (a CF3 has its own row)',
     'α-halo carbonyls are more electrophilic at C=O and far more reactive in SN2; α-alkoxy / α-amino groups chelate Lewis acids (chelation-controlled addition)', 'induction withdraws electron density; the heteroatom can bind the metal', _CL)
site('co-alpha-CF3', 'CF3 on the carbonyl carbon', 'carbonyl carbon', 'electronic', 'carbonyl C', '[C;X3;$(C=[O,N,S]);$(C-[C;!900*](F)(F)F)]', 'the retained carbonyl carbon is bonded to a CF3 it already had (trifluoromethyl ketone, trifluoroacetamide)',
     'much more electrophilic: hydrates readily and adds nucleophiles without catalysts; trifluoroacetamides are cleaved easily', 'the CF3 withdraws strongly (Taft σ* +2.61)', 'Hansch, Leo, Taft; Taft σ* values')
site('co-two-hetero', 'two heteroatoms on the carbonyl carbon (carbamate, urea, carbonate)', 'carbonyl carbon', 'electronic', 'carbonyl C', '[C;X3;$(C=O);$(C(-[O,N])-[O,N])]', 'the retained carbonyl carbon carries two O/N substituents',
     'poorly electrophilic: carbamates and ureas survive conditions that cleave amides and esters', 'two lone pairs donate into the C=O', _CL)
site('co-heteroaryl', 'azinyl on the carbonyl carbon', 'carbonyl carbon', 'electronic', 'carbonyl C', '[C;X3;$(C=[O,N,S]);$(C-[c;$(c:[nX2]);!900*])]', 'the retained carbonyl carbon is bonded to a ring carbon next to a pyridine-type ring N (picolinoyl-type; pyrrole-type N–H / N–R do not count)',
     'more electrophilic than a phenyl carbonyl; the ring N chelates metals and directs', 'the azine is electron-poor and can bind', _JM)
for n in (4, 5, 6):
    site(f'co-ring{n}', f'carbonyl in a {n}-membered ring', 'carbonyl carbon', 'geometry', 'carbonyl C', f'[C;X3;$(C=[O,N,S]);r{n}]', f'the retained carbonyl carbon is in a ring of {n} atoms (lactam, lactone, cyclic ketone)',
         {4: 'strained ring carbonyl: β-lactams and cyclobutanones are very reactive to nucleophiles',
          5: 'γ-lactones and γ-lactams are the least strained and the most stable lactones and lactams; cyclopentanone adds nucleophiles more slowly than cyclohexanone (adding to the C=O increases eclipsing in the five-membered ring)',
          6: 'δ-valerolactone is strained and opens faster than γ-butyrolactone; cyclohexanone adds nucleophiles faster than cyclopentanone and acyclic ketones (the sp2 → sp3 change relieves torsional strain); axial or equatorial attack decides the diastereomer'}[n],
         'ring size sets strain and the face that is attacked', _CL)
# ---------------------------------------------------------------------------------------------------- N, O or S is the reacting atom
_ARYL_REL = {'N': 'aniline, N-aryl amide', 'O': 'phenol, aryl ether'}
_ARYL_EFFECT = {
    'N': 'a weaker nucleophile and base than an alkylamine (the lone pair is delocalised into the ring: aniline pKaH 4.6 against 10.6), but a stronger acid; in acylation it reacts after the alkylamine and before an alcohol',
    'O': 'a weaker nucleophile than an alcohol but a much stronger acid (phenol pKa 10 against 16 for an alcohol): easily deprotonated, so alkylation with a mild base goes at a phenol before an alcohol; in acylation it reacts after an alcohol'}
_HETARYL_REL = {'N': '2-aminopyridine, 2-aminopyrimidine, 2-aminothiazole', 'O': '2-alkoxypyridine, 2-aryloxypyridine'}
_HETARYL_EFFECT = {
    'N': 'a much weaker nucleophile at the exocyclic N (amidine-like delocalisation); alkylation tends to go on the ring N instead',
    'O': 'a poorer donor than an alkoxy on benzene (the lone pair is shared with the C=N of the ring); the C–O bond is cleaved by acid, giving the pyridone'}
_RING_REL = {'N': 'piperidine, pyrrolidine, morpholine, lactam', 'O': 'THF, THP, dioxane, lactone, cyclic ether'}
_RING_EFFECT = {
    'N': 'ring amines are more nucleophilic than acyclic ones of the same basicity (the alkyl groups are tied back)',
    'O': 'ring ethers are better donors than acyclic ones (the alkyl groups are tied back) and open under acid'}
for e, nm in (('N', 'nitrogen'), ('O', 'oxygen')):
    site(f'{e.lower()}-aryl', f'{nm} on an aromatic ring', nm, 'electronic', nm, f'[{e};$({e}-[c;!900*])]', f'the reacting {nm} was bonded to an aromatic carbon ({_ARYL_REL[e]})',
         _ARYL_EFFECT[e], 'conjugation of the lone pair with the ring',
         _CL + '; Reactivity_Order: amine > aniline > alcohol > phenol for acylation')
    site(f'{e.lower()}-heteroaryl', f'{nm} on an azine carbon', nm, 'electronic', nm, f'[{e};$({e}-[c;$(c:[nX2]);!900*])]', f'the reacting {nm} was bonded to a ring carbon next to a pyridine-type ring N ({_HETARYL_REL[e]})',
         _HETARYL_EFFECT[e], 'the lone pair is shared with the C=N of the ring', _JM)
    site(f'{e.lower()}-alpha-branched', f'{nm} on a branched carbon', nm, 'steric', nm, f'[{e};$({e}-[CX4;D3,D4;!900*])]', f'the reacting {nm} was bonded to a CH or quaternary carbon (isopropyl, tert-butyl)',
         'slower acylation, alkylation and coupling; tert-alkyl amines and alcohols are poor nucleophiles; the less hindered nucleophile of the molecule reacts first', 'the branch crowds the lone pair', 'Taft steric constants Es: i-Pr −0.47, t-Bu −1.54')
    site(f'{e.lower()}-benzylic', f'{nm} on a benzylic carbon', nm, 'electronic', nm, f'[{e};$({e}-[CX4;$(C-c);!900*])]', f'the reacting {nm} was bonded to a CH2 or CH attached to an aromatic ring',
         'an ordinary nucleophile; the benzylic C–' + e + ' bond is cleaved by hydrogenolysis (protecting group) and by acid for electron-rich rings', 'the ring stabilises the benzylic cation and radical', 'Greene, Protective Groups in Organic Synthesis')
    site(f'{e.lower()}-alpha-ewg', f'{nm} next to an electron-withdrawing carbon', nm, 'electronic', nm, f'[{e};$({e}-[CX4;$(C-[$(C(F)(F)F),$(C#N),$(C=O)]);!900*])]', f'the reacting {nm} was bonded to a carbon that carries CF3, CN or C=O',
         'weaker nucleophile and base (2,2,2-trifluoroethylamine pKaH 5.7 against 10.7 for ethylamine) and stronger acid (CF3CH2OH pKa 12.4 against 15.9 for ethanol)', 'the withdrawing group pulls on the lone pair', 'Taft σ*; Evans pKa table')
    site(f'{e.lower()}-in-ring', f'{nm} in a ring', nm, 'steric', nm, f'[{e};R]', f'the reacting {nm} is part of a ring ({_RING_REL[e]})',
         _RING_EFFECT[e], 'the substituents are held away from the lone pair', _CL)
    site(f'{e.lower()}-ring-small', f'{nm} in a three- or four-membered ring', nm, 'geometry', nm, f'[{e};r3,r4]', f'the reacting {nm} is in an epoxide / aziridine / oxetane / azetidine ring in the product',
         'a strained ring was formed or kept: forming it needs a reactive pair (intramolecular SN2, peracid, nitrene); the ring opens easily under acid or with nucleophiles', 'ring strain of about 26 kcal/mol', _CL + ' (ring strain)')
    site(f'{e.lower()}-ring-medium', f'{nm} in a seven- to nine-membered ring', nm, 'geometry', nm, f'[{e};r7,r8,r9]', f'the reacting {nm} is in a ring of seven to nine atoms in the product',
         'medium rings close slowly (entropy and transannular strain): seven- to nine-membered lactams and ethers need high dilution or templates', 'transannular strain and a low effective molarity', _CL)
    site(f'{e.lower()}-ring-macro', f'{nm} in a macrocycle', nm, 'geometry', nm, f'[{e};R;!r3;!r4;!r5;!r6;!r7;!r8;!r9]', f'the reacting {nm} is in a ring of ten or more atoms in the product (macrolactamisation, macrolactonisation, ring-closing)',
         'macrocycle closures need high dilution (0.001–0.01 M), slow addition or a template; competing oligomerisation is the usual failure', 'entropy of ring closure against intermolecular reaction', _CL)
site('n-acyl-origin', 'nitrogen that was an amide / carbamate / sulfonamide N', 'nitrogen', 'electronic', 'nitrogen', '[N;X3;$(N-[C,S;X3,X4;!900*]=[O,S])]', 'the reacting N was bonded to a C=O, C=S or S=O (amide, carbamate, sulfonamide, imide)',
     'a poor nucleophile (lone pair delocalised) but an acid: reacts after deprotonation (NaH, Cs2CO3) or under Mitsunobu / Buchwald / Chan–Lam conditions', 'the lone pair is shared with the carbonyl or sulfonyl group', _CL)
site('n-adjacent-hetero', 'nitrogen next to another heteroatom (hydrazine, hydroxylamine)', 'nitrogen', 'electronic', 'nitrogen', '[N;X3;$(N-[N,O,n;!900*])]', 'the reacting N is bonded to another N or O (α-effect nucleophile)',
     'more nucleophilic than its pKaH suggests (the α-effect): hydrazines and hydroxylamines react at the terminal N', 'repulsion between adjacent lone pairs raises the HOMO', _CL)
site('n-primary-origin', 'nitrogen with one carbon substituent (primary amine before)', 'nitrogen', 'steric', 'nitrogen', f'[N;X3;$(N-{C_});!$(N(-{C_})-{C_})]', 'the reacting N had one carbon substituent (primary amine)',
     'the least hindered amine: acylation, alkylation and reductive amination are fastest; over-alkylation is the usual problem', 'two N–H bonds: little steric hindrance', _CL)
site('n-secondary-origin', 'nitrogen with two carbon substituents (secondary amine before)', 'nitrogen', 'steric', 'nitrogen', f'[N;X3;$(N(-{C_})-{C_});!$(N(-{C_})(-{C_})-{C_})]', 'the reacting N had two carbon substituents (secondary amine)',
     'more nucleophilic than primary (more electron-rich) but more hindered: acylates well, alkylates cleanly once', 'two alkyl groups push electrons in but add bulk', _CL)
site('n-tertiary-origin', 'nitrogen with three carbon substituents (tertiary amine before)', 'nitrogen', 'steric', 'nitrogen', f'[N;X3,X4;$(N(-{C_})(-{C_})-{C_})]', 'the reacting N had three carbon substituents (tertiary amine)',
     'a base and nucleophilic catalyst; quaternises with alkyl halides; hindered tertiary amines (DIPEA) only deprotonate', 'three alkyl groups: strongest base, most hindered', _CL)
site('o-primary-origin', 'oxygen on a primary carbinol carbon (R–CH2–O or CH3–O)', 'oxygen', 'steric', 'oxygen', '[O;X2;$(O-[CX4;H2,H3;!900*]);!$(O(-[#6;!900*])-[#6;!900*])]', 'the reacting O is bonded to a CH2 or CH3 and to nothing but the new bond (a primary alcohol or methanol before; ethers are excluded)',
     'esterifies, oxidises, substitutes and alkylates faster at a primary carbinol than at a secondary or tertiary one', 'primary alcohols are the least hindered', 'Reactivity_Order: esterification, oxidation, substitution')
site('s-aryl', 'sulfur on an aromatic ring', 'sulfur', 'electronic', 'sulfur', '[S;$(S-[c;!900*])]', 'the reacting S was bonded to an aromatic carbon (thiophenol, aryl sulfide)',
     'arylthiolates (thiophenol pKa 6.6) are weaker bases than alkylthiolates (pKa about 10.5) and still soft nucleophiles: they do SNAr and Pd / Cu couplings; aryl sulfides oxidise more slowly', 'conjugation of the S lone pair with the ring', _CL)
site('s-alkyl', 'sulfur on an alkyl carbon', 'sulfur', 'electronic', 'sulfur', '[S;$(S-[CX4;!900*])]', 'the reacting S was bonded to an sp3 carbon (thiol, dialkyl sulfide)',
     'soft and strongly nucleophilic: alkylates, adds to Michael acceptors, oxidises to sulfoxide and sulfone', 'a large polarisable lone pair', _CL)
for e, nm in (('N', 'nitrogen'), ('O', 'oxygen'), ('S', 'sulfur')):
    site(f'{e.lower()}-aryl-ortho2', f'{nm} on an aromatic ring with two ortho substituents', nm, 'steric', nm, f'[{e};$({e}-c(:c-[!#1;!900*]):c-[!#1;!900*])]', f'the reacting {nm} was bonded to a ring carbon whose both neighbours carry substituents (2,6-disubstituted aniline, phenol or thiophenol)',
         'strongly hindered: acylation, alkylation and coupling are slow (2,6-dimethylaniline acylates far more slowly than aniline; Buchwald needs bulky ligands)', 'the two ortho groups flank the lone pair', _CL)

# an aryl ring on the reacting atom: the Hammett effect of its substituents on the nucleophilicity / electrophilicity of that atom
EWG_SUBS = [x[2] for x in SUBSTITUENTS if x[4] >= 0.30]          # the substituents with sigma-para of at least +0.30 ...
EDG_SUBS = [x[2] for x in SUBSTITUENTS if x[4] <= -0.15]          # ... and of at most -0.15 (the Hammett effect of the ring on the atom it is bonded to)
_ARYL_ANCHORS = [('n', 'N', 'nitrogen', 'nitrogen', 'N'), ('o', 'O', 'oxygen', 'oxygen', 'O'), ('s', 'S', 'sulfur', 'sulfur', 'S'),
                 ('co', 'C;X3;$(C=[O,N,S])', 'carbonyl C', 'the retained carbonyl carbon', 'C'), ('bz', 'C;X4', 'sp3 C', 'the reacting sp3 carbon', 'C')]
_ARYL_TEXT = {
    ('n', 'ewg'): 'weaker nucleophile and base (p-nitroaniline pKaH 1.0 against 4.6 for aniline) and stronger acid: reacts after the more electron-rich amine of the molecule; needs a stronger base or a more reactive electrophile',
    ('n', 'edg'): 'stronger nucleophile and base (p-anisidine pKaH 5.3, p-toluidine 5.1 against 4.6 for aniline): reacts before less electron-rich anilines',
    ('o', 'ewg'): 'a weaker nucleophile and a stronger acid (p-nitrophenol pKa 7.2 against 10.0 for phenol): easily deprotonated (K2CO3 is enough) but a poorer nucleophile; a better leaving group',
    ('o', 'edg'): 'a more nucleophilic, less acidic phenol (p-methoxyphenol pKa 10.2)', ('s', 'ewg'): 'a weaker, less basic thiolate nucleophile; stronger acid',
    ('s', 'edg'): 'a more basic, more nucleophilic thiolate',
    ('co', 'ewg'): 'more electrophilic carbonyl: acyl substitution and addition are faster (p-nitrobenzoyl chloride reacts much faster than benzoyl chloride)',
    ('co', 'edg'): 'less electrophilic carbonyl (the ring donates into the C=O): slower acyl substitution and addition',
    ('bz', 'ewg'): 'SN2 slightly faster, SN1 much slower (the benzylic cation is destabilised)',
    ('bz', 'edg'): 'SN1 far faster (p-methoxybenzyl chloride solvolyses orders of magnitude faster than benzyl chloride: the PMB group is acid-labile); SN2 about as for benzyl'}
# an ortho group acts through the ring too, but it also crowds the atom and can hydrogen-bond to it, so the para text does not apply to the nitrogen, oxygen and sulfur rows
_ARYL_TEXT_ORTHO = {
    ('n', 'ewg'): 'a far weaker nucleophile and base (o-nitroaniline pKaH −0.3 against 4.6 for aniline; the nitro group also hydrogen-bonds to the N–H and crowds it): needs a stronger base or a more reactive electrophile; reacts after other amines',
    ('n', 'edg'): 'about as basic as aniline (o-toluidine pKaH 4.4, o-anisidine 4.5 against 4.6) and more hindered: the ortho group crowds the N and twists it out of conjugation, so acylation and coupling are slower than for the para isomer',
    ('o', 'ewg'): 'a stronger acid and a weaker nucleophile (o-nitrophenol pKa 7.2 against 10.0 for phenol; the nitro group hydrogen-bonds to the O–H and crowds the O)',
    ('o', 'edg'): 'about as acidic as phenol (o-cresol pKa 10.3, guaiacol 10.0) and more hindered; an ortho alkoxy group can chelate a metal or a Lewis acid',
    ('s', 'ewg'): 'a weaker, less basic thiolate nucleophile and a stronger acid; the ortho group also crowds the S',
    ('s', 'edg'): 'a thiolate about as nucleophilic as thiophenolate, more hindered by the ortho group'}
for ak, head, anchor, what, el in _ARYL_ANCHORS:
    for cls, subs_, clsname in (('ewg', EWG_SUBS, 'electron-withdrawing'), ('edg', EDG_SUBS, 'electron-donating')):
        for pos, ringfmt, posname in (('p', '[c;!900*]1:a:a:c(-{g}):a:a:1', 'para'), ('o', '[c;!900*]:c-{g}', 'ortho')):
            alts = ','.join(f'$({el}-{ringfmt.format(g=guard(x))})' for x in subs_)
            site(f'{ak}-aryl-{pos}-{cls}', f'{what} on an aryl ring with a {posname} {clsname} group', 'aryl ring on the reacting atom', 'electronic', anchor,
                 f'[{head};{alts}]', f'the reacting atom is bonded to an aromatic carbon (unchanged) whose {posname} position carries {"nitro, cyano, CF3, a carbonyl, a sulfonyl, an ammonium, an acyloxy, OCF3 or a sulfonyloxy group" if cls == "ewg" else "an alkoxy, hydroxy, phenoxide, amino or alkyl group"}',
                 (_ARYL_TEXT_ORTHO.get((ak, cls)) if pos == 'o' and (ak, cls) in _ARYL_TEXT_ORTHO else _ARYL_TEXT[(ak, cls)] + (' An ortho group also crowds the atom.' if pos == 'o' else '')),
                 'Hammett effect through the ring (σ values of the substituent)', _HLT)

# ---------------------------------------------------------------------------------------------------- an alkene carbon that is still an alkene is the reacting atom
site('alkene-terminal', 'carbon of a vinyl group (CH=CH2)', 'alkene carbon', 'steric', 'alkene C', '[C;X3;H1;$(C=[C;H2])]', 'the reacting carbon is the CH of a vinyl group',
     'the least hindered alkene: Heck, metathesis and hydroboration are fast; electrophiles end up on the CH2 (Markovnikov)', 'one substituent only', _CL)
site('alkene-conjugated-carbonyl', 'alkene conjugated with an electron-withdrawing group (enone, acrylate, acrylonitrile, nitroalkene, vinyl sulfone)', 'alkene carbon', 'electronic', 'alkene C',
     '[C;X3;$(C=C-[$([CX3]=[O,N]),$([CX2]#N),$([N+](=O)[O-]),$([SX4](=O)=O)]);!$(C=O)]', 'the reacting carbon is part of a C=C next to a C=O, C≡N, NO2 or SO2',
     'electron-poor alkene: nucleophiles add at the β carbon, electrophiles react slowly; a good dienophile', 'the acceptor withdraws electrons from the C=C', _CL)
site('alkene-styrenyl', 'alkene conjugated with an aromatic ring', 'alkene carbon', 'electronic', 'alkene C', '[C;X3;$(C=C-c);!$(C=O)]', 'the reacting carbon is part of a C=C next to an aromatic ring (styrene)',
     'electrophiles add to give the benzylic cation; hydrogenation and radical additions are fast', 'the ring stabilises the cation or radical', _CL)
site('alkene-electron-rich', 'alkene with an oxygen or nitrogen donor (enol ether, enamine)', 'alkene carbon', 'electronic', 'alkene C', '[C;X3;$(C=C-[O,N;X2,X3;+0;!$(*-[C,S,P]=[O,S,N])])]', 'the reacting carbon is part of an enol ether, enamine or a simple vinyl ether / vinyl amine (not an enamide, enol ester or enol triflate)',
     'very electron-rich: reacts with electrophiles at the β carbon (Vilsmeier, halogenation, Mannich) and hydrolyses in mild acid', 'the donor lone pair raises the HOMO of the C=C', _CL)
site('alkene-two-carbons', 'alkene carbon with two carbon substituents (no H)', 'alkene carbon', 'steric', 'alkene C', f'[C;X3;$(C(=C)(-{C_})-{C_})]', 'the reacting carbon has two carbon substituents besides the double bond (1,1-disubstituted, trisubstituted or tetrasubstituted alkene)',
     'more electron-rich (faster with electrophiles: epoxidation, halogenation) but more hindered for catalysts (hydrogenation, hydroboration slower)', 'alkyl groups donate and crowd', _CL)
site('alkene-new-bond-to-C', 'retained alkene carbon that took a new carbon', 'alkene carbon', 'steric', 'alkene C', '[C;X3;$(C=C);$(C-[900C,900c])]', 'the alkene carbon is bonded to a carbon the reaction put there (Heck, Suzuki or Negishi of a vinyl partner)',
     'Heck: the aryl goes to the less substituted carbon of the alkene (neutral Pd) or the more substituted one (cationic Pd, electron-rich alkene)', 'sterics against polarisation in the insertion step', _CL)

# The patterns that really are hindrance: the page's "Hindrance" line shows only these, and only when the recorded reactions penalise them. The other 'steric' rows describe the degree of the
# reacting atom (methyl, primary, secondary), which end of an addition it is, or the ring it sits in: they are not hindrance by themselves.
HINDRANCE_KEYS = {'ring-peri', 'ring-ortho-1', 'ring-ortho-2', 'ring-ortho-branched', 'azole-n-hindered', 'azole-n-hindered-2', 'azine-n-hindered', 'azine-n-hindered-2', 'sp3-neopentyl', 'sp3-beta-branched',
                  'co-alpha-branched', 'co-alpha-quaternary', 'n-alpha-branched', 'o-alpha-branched', 'n-aryl-ortho2', 'o-aryl-ortho2', 's-aryl-ortho2'}
for _s in SITE:
    if _s['key'] in HINDRANCE_KEYS: _s['kind'] = 'hindrance'
SITE_BY_KEY = {s['key']: s for s in SITE}
assert len(SITE_BY_KEY) == len(SITE), 'duplicate pattern keys'

# ---------------------------------------------------------------------------------------------------- ADVISORIES
# reagent side = alternatives (any one); each alternative is a list of tokens that must ALL hold for the reagent as typed:
#   frag:<fragment ID found in a reagent part>   smiles:<a reagent part exactly as the Reagents sheet types it>   role:<role (or second role) of a part in the Reagents sheet>
#   solvent:<Solvents-sheet type of a solvent in the conditions box: protic, polar aprotic ...>
# [] = any reagent (or none). group side = fragment IDs (any one present in the starting materials; group_where 'reagent' or 'both' also looks in the reagent parts);
# show = 'bystander' (a group the leading product did NOT consume: what else the reagent does), 'target' (a group it consumed) or 'both'; min_distinct = at least that many
# different groups present; page = False: knowledge only, not tested on the page
T_HYDROXIDE = [['frag:F060', 'not:smiles:CCOC([O-])[O-]', 'not:smiles:COC([O-])[O-]'], ['smiles:[Li][OH]'], ['smiles:CC(C)(C)[O][K]'], ['smiles:CC(C)(C)[O][Na]'], ['smiles:C[O][Na]'], ['smiles:C[Si](C)(C)[O-]'], ['smiles:C[Si](C)(C)[O][K]']]            # hydroxide, alkoxides, silanolate (F060 finds HO⁻, MeO⁻, EtO⁻, t-BuO⁻ as ions; two mis-drawn orthoformate parts of the Reagents sheet look like alkoxides and are left out)
T_ORGANOMETALLIC = [['frag:F029'], ['frag:F140'], ['frag:F141'], ['smiles:[CH3][Mg+]'], ['smiles:C[CH2][Mg+]'], ['smiles:[CH2-]CCC']]                  # Grignard reagents and organolithiums (also the cation / anion forms the Reagents sheet types)
T_EXCHANGE = [['smiles:[Li][CH2]CCC'], ['smiles:[Li][C](C)(C)C'], ['smiles:[Li][CH](C)CC'], ['smiles:C[CH](C)[Mg][Cl]'], ['smiles:C[CH](C)[Mg][Br]']]                 # n-BuLi, t-BuLi, s-BuLi, iPrMgCl, iPrMgBr: the reagents that swap a halogen for a metal (MeMgX, EtMgX, PhMgX, MeLi add or deprotonate instead)
T_STRONG_BASE = T_HYDROXIDE + [['frag:F142'], ['frag:F143'], ['frag:F140'], ['frag:F029'], ['frag:F141'], ['smiles:[NH2-]'], ['smiles:[Li][N]([Si](C)(C)C)[Si](C)(C)C']]
T_STRONG_ACID = [['smiles:O=C(O)C(F)(F)F'], ['smiles:Cl'], ['smiles:Br'], ['smiles:I'], ['smiles:F'], ['smiles:O=S(=O)(O)O'], ['smiles:Cc1ccc(S(=O)(=O)O)cc1'], ['smiles:CS(=O)(=O)O'], ['smiles:O=S(=O)(O)C(F)(F)F'], ['smiles:CC1(C)C2CCC1(CS(=O)(=O)O)C(=O)C2'], ['smiles:[O-][Cl+3]([O-])([O-])O'], ['smiles:O=CO'], ['smiles:O=P(O)(O)O'], ['smiles:[H+]']]                  # not pyridinium (PPTS) or bisulfate: too weak to remove a Boc group
T_LEWIS = [['smiles:FB(F)F'], ['smiles:BrB(Br)Br'], ['smiles:ClB(Cl)Cl'], ['smiles:[Cl][Ti]([Cl])([Cl])[Cl]'], ['smiles:[Cl][Sn]([Cl])([Cl])[Cl]'], ['smiles:[Cl][Fe]([Cl])[Cl]'], ['smiles:[Al+3]', 'smiles:[Cl-]'], ['smiles:[Zn+2]', 'smiles:[Cl-]'], ['smiles:C[Si](C)(C)OS(=O)(=O)C(F)(F)F']]
T_ETHER_CLEAVING = [['smiles:BrB(Br)Br'], ['smiles:ClB(Cl)Cl'], ['smiles:Br'], ['smiles:I'], ['smiles:C[Si](C)(C)I'], ['smiles:C[Si](C)(C)Br']]
T_ALH = [['frag:F042'], ['frag:F143', 'smiles:[Al+3]'], ['smiles:COCC[O][AlH2-][O]CCOC']]                    # LiAlH4 (typed [AlH4-], or as hydride + Al3+), Red-Al
T_HYDRIDE = T_ALH + [['frag:F143'], ['frag:F145'], ['smiles:CCC(C)[BH-](C(C)CC)C(C)CC']]
T_NABH4 = [['smiles:[BH4-]', 'smiles:[Na+]', 'not:smiles:[Li+]', 'not:role:Lewis acid', 'not:elem:Co', 'not:elem:Ni', 'not:elem:Cu', 'not:elem:Pd', 'not:elem:Pt', 'not:elem:Ca', 'not:smiles:II']]                  # NaBH4 alone: with LiCl, CaCl2, CoCl2, NiCl2, I2, a Lewis acid or a metal catalyst it is another reagent (LiBH4, borane, a metal boride)
T_LIBH4 = [['smiles:[BH4-]', 'smiles:[Li+]']]
T_H2_PGM = [['frag:F164', 'smiles:[H][H]']]
T_H2_NOT_NI = [['frag:F164', 'smiles:[H][H]', 'not:elem:Ni']]
T_PGM = [['frag:F164']]
T_PD = [['elem:Pd']]                                                                       # any palladium part (Pd/C, Pd(PPh3)4, Pd(OAc)2 ...): only Pd(0) ionises allyl groups
T_PEROXIDE = [['smiles:O=C(OO)c1cccc(Cl)c1'], ['smiles:CC(=O)OO'], ['smiles:OO'], ['smiles:O=S(=O)([O-])O[O-]'], ['smiles:O=S([O-])OO']]              # mCPBA, peracetic acid, H2O2, Oxone
T_METAL_OXO = [['smiles:CC(=O)OI1(OC(C)=O)(OC(C)=O)OC(=O)c2ccccc21'], ['smiles:O=C(O)c1ccccc1I(=O)=O'], ['smiles:CC(=O)OI(OC(C)=O)c1ccccc1'], ['smiles:O=C(OI(OC(=O)C(F)(F)F)c1ccccc1)C(F)(F)F'],
               ['smiles:[O]=[Cr](=[O])([O-])[Cl]'], ['smiles:[O]=[Cr](=[O])=[O]'], ['smiles:[O]=[Cr](=[O])([O-])[O][Cr](=[O])(=[O])[O-]'], ['smiles:[O]=[Cr](=[O])([OH])[OH]'],
               ['smiles:[O]=[Mn](=[O])(=[O])[O-]'], ['smiles:[O]=[Ru](=[O])(=[O])[O-]'], ['smiles:[O]=[Ru]=[O]']]
# Dess–Martin, IBX, PIDA, PIFA, PCC, CrO3, dichromate, chromic acid, permanganate, perruthenate, RuO2: the oxidants of alcohols. Not periodate / iodate (they cleave diols), MnO2 (allylic and benzylic alcohols only), OsO4, SeO2
T_ACYLATING = [['frag:F020'], ['frag:F041'], ['frag:F064'], ['role:coupling reagent (activates an acid)'], ['smiles:CC(=O)OC(C)=O']]
T_ALKYLATING = [['role:electrophile / alkylating agent', 'frag:' + f] for f in ('F010', 'F011', 'F012', 'F036', 'F174', 'F066', 'F125')]            # alkyl halides, sulfonates and α-halo carbonyls the Reagents sheet calls alkylating agents (not CCl4, CBr4, dichloroethane, Selectfluor, triphosgene ... which match the same fragments)
T_HALOGENATING = [['smiles:BrBr'], ['smiles:O=C1CCC(=O)N1Br'], ['smiles:CC1(C)C(=O)N(Br)C(=O)N1Br'], ['smiles:ClBr']]                  # Br2, NBS, DBDMH, BrCl: the brominating agents (not I2, ICl, NCS, NIS, Cl2)
T_DEHYDRATING = [['smiles:O=S(Cl)Cl'], ['smiles:O=P(Cl)(Cl)Cl'], ['smiles:ClP(Cl)(Cl)(Cl)Cl'], ['smiles:O=C(Cl)C(=O)Cl'], ['smiles:O=S(=O)(OS(=O)(=O)C(F)(F)F)C(F)(F)F'], ['smiles:O=S(Br)Br'], ['smiles:BrP(Br)Br'], ['smiles:ClP(Cl)Cl'], ['smiles:O=P(Br)(Br)Br']]
T_STRONG_BASE_ONLY = [['frag:F142'], ['frag:F143'], ['frag:F140'], ['frag:F029'], ['frag:F141'], ['smiles:[NH2-]'], ['smiles:[Li][N]([Si](C)(C)C)[Si](C)(C)C'], ['smiles:CC(C)(C)[O-]'], ['smiles:CCC(C)(C)[O-]'],
                      ['smiles:CC(C)(C)[O][K]'], ['smiles:CC(C)(C)[O][Na]'], ['smiles:C[Si](C)(C)[O][K]']]                 # strong enough to take a proton from an alcohol, an amide N–H or a ketone: not hydroxide, methoxide or ethoxide
ACIDIC_XH = ['F001', 'F118', 'F003', 'F014', 'F015', 'F016', 'F017', 'F129', 'F197', 'F198', 'F199', 'F028', 'F072', 'F002', 'F034', 'F026', 'F094', 'F176', 'F177', 'F178', 'F181', 'F115', 'F182', 'F183',
             'F096', 'F097', 'F095', 'F120', 'F185', 'F186', 'F187', 'F067', 'F068', 'F089', 'F090', 'F109', 'F110', 'F111', 'F114', 'F008', 'F048', 'F169', 'F102', 'F049', 'F099', 'F056', 'F150', 'F098']
# the same groups without the ones that are a subset of another in the list (an alcohol is F014 and, if primary, F003 too): one entry per kind of acidic site, for advisories that count sites
ACIDIC_LEAF = ['F001', 'F014', 'F028', 'F072', 'F176', 'F026', 'F034', 'F094', 'F096', 'F097', 'F095', 'F067', 'F068', 'F089', 'F090', 'F008', 'F048', 'F169', 'F102', 'F192', 'F193', 'F194', 'F195', 'F049', 'F099', 'F056', 'F098', 'F150']
CARBONYL_LIKE = ['F004', 'F139', 'F138', 'F021', 'F120', 'F184', 'F185', 'F202', 'F022', 'F191', 'F001', 'F118', 'F020', 'F041', 'F018', 'F117', 'F019', 'F116', 'F031', 'F040', 'F124', 'F121', 'F136', 'F074', 'F079', 'F100', 'F122']
NUCLEOPHILES = ['F002', 'F034', 'F026', 'F094', 'F176', 'F177', 'F178', 'F181', 'F115', 'F003', 'F014', 'F015', 'F028', 'F072', 'F090', 'F129', 'F017', 'F197', 'F198']
LEAVING_GROUPS = ['F010', 'F011', 'F012', 'F013', 'F047', 'F084', 'F065', 'F066', 'F087', 'F085']
ADVISORY = []


def advice(key, name, kind, effect_kind, reagent, groups, effect, why, source, needs='', group_where='substrate', show='bystander', min_distinct=1, page=True):
    ADVISORY.append(dict(key=key, name=name, kind=kind, effect_kind=effect_kind, reagent=reagent, groups=list(groups), effect=effect, why=why, source=source, needs=needs,
                         group_where=group_where, show=show, min_distinct=min_distinct, page=page))


advice('quench-organometallic-acidic-XH', 'Grignard or organolithium + a group with an acidic X–H', 'quencher', 'quenches', T_ORGANOMETALLIC, ACIDIC_XH,
       'the reagent is destroyed first: each O–H, N–H, S–H, CO2H or terminal C–H protonates one equivalent of R–M (R–H forms); protect the group or use an extra equivalent per acidic H',
       'the carbon bound to the metal is a very strong base (conjugate acid pKa 45–50); an O–H, N–H or S–H is deprotonated at once (C–H acids such as ketone α-H compete with addition to the carbonyl, and addition usually wins for aldehydes and unhindered ketones)',
       'Clayden, Organic Chemistry (organometallics); Fragment_Conditions F029')
advice('quench-strong-base-acidic-XH', 'a strong base + several different acidic groups', 'competing site', 'competes', T_STRONG_BASE_ONLY, ACIDIC_LEAF,
       'the base removes the most acidic proton first and a second equivalent the next: with several sites the product follows the acidity order (carboxylic acid 4–5, phenol 10, alcohol 16, amide N–H 17, ketone α-C–H 19–20, ester α-C–H 25, terminal alkyne 25)',
       'the equilibrium constant of a proton transfer is 10^ΔpKa', 'Evans pKa table; Clayden, Organic Chemistry', min_distinct=2)
advice('incompatible-hydride-acidic-XH', 'a metal hydride (LiAlH4, DIBAL-H, NaH) + O–H, N–H or CO2H', 'quencher', 'quenches', T_HYDRIDE, ACIDIC_XH,
       'H2 is evolved and one hydride used up per acidic proton; the alkoxide / amide formed may react more slowly than the free group (carboxylic acids are reduced via the carboxylate only with excess LiAlH4)',
       'hydride is a strong base (conjugate acid H2, pKa about 35)', 'Clayden, Organic Chemistry (reduction); Fragment_Conditions F042')
advice('also-reduced-LiAlH4', 'LiAlH4 + other reducible groups', 'incompatible group', 'also_reacts', T_ALH, CARBONYL_LIKE + ['F023', 'F025', 'F056', 'F055', 'F010', 'F011', 'F012', 'F044', 'F045', 'F066', 'F077', 'F049', 'F125'],
       'also reduced: esters, amides, nitriles, acids, acyl chlorides, epoxides, azides, carbamates (to N-methyl) and most alkyl halides and sulfonates; aliphatic nitro compounds go to amines and nitroarenes to azo compounds; unless the amount is limited the molecule is reduced at every one of them',
       'a very strong hydride donor with little selectivity between polar functional groups', 'Clayden, Organic Chemistry (reduction)')
advice('spared-by-NaBH4', 'NaBH4 + an amide, carbamate or nitro group', 'selectivity', 'spared', T_NABH4, ['F021', 'F077', 'F120', 'F025'],
       'usually survive: NaBH4 alone reduces aldehydes, ketones, acyl chlorides and (with acid) imines at 0–25 °C but leaves amides, carbamates and nitro groups; esters, nitriles and acids also survive without an additive, but in the recorded reactions they are reduced about as often as with other reagents (they are the target, or LiCl / CaCl2 / CoCl2 / I2 / heat is added)',
       'the N lone pair makes an amide or carbamate C=O far less electrophilic than a ketone, and NaBH4 is a mild hydride donor', 'Clayden, Organic Chemistry; Known_No_Reaction sheet (NaBH4)')
advice('also-reduced-LiBH4', 'LiBH4 + an ester, lactone, epoxide or acid derivative', 'incompatible group', 'also_reacts', T_LIBH4, ['F004', 'F139', 'F138', 'F023', 'F124', 'F020', 'F041', 'F121'],
       'also reduced: unlike NaBH4, LiBH4 reduces esters, lactones, epoxides and acyl derivatives (the Li⁺ activates the C=O); it still leaves acids, amides, nitriles and nitro groups largely alone',
       'the lithium ion is a better Lewis acid than sodium', 'Clayden, Organic Chemistry (reduction)')
advice('also-oxidised-by-peroxide', 'a peroxide or peracid (H2O2, mCPBA, peracetic acid, Oxone) + a sulfide, tertiary amine, alkene or electron-rich ring', 'incompatible group', 'also_reacts', T_PEROXIDE,
       ['F073', 'F072', 'F093', 'F179', 'F180', 'F005', 'F006', 'F126', 'F127', 'F203', 'F106', 'F105', 'F104', 'F028', 'F026', 'F019', 'F091'],
       'also oxidised: sulfides to sulfoxide and sulfone, tertiary amines and pyridines to N-oxides, alkenes to epoxides, enol ethers and enamines, electron-rich heteroarenes, phenols and anilines (to quinones and tars), ketones to esters (Baeyer–Villiger)',
       'a peracid is an electrophilic oxygen donor: it takes electrons from whichever group is the most electron-rich first', 'Clayden, Organic Chemistry (oxidation)')
advice('also-oxidised-by-metal-oxo', 'a chromium(VI), permanganate, ruthenium or hypervalent-iodine oxidant + other oxidisable groups', 'incompatible group', 'also_reacts', T_METAL_OXO,
       ['F003', 'F015', 'F199', 'F073', 'F072', 'F130'],
       'also oxidised: primary and secondary alcohols (to aldehydes / ketones; tertiary alcohols are not oxidised), sulfides and thiols, and benzyl ethers (chromic acid, permanganate, ruthenium tetroxide); permanganate and chromic acid carry aldehydes on to acids, while Dess–Martin periodinane and IBX stop at the aldehyde or ketone; alkenes and benzylic C–H are attacked only by the strongest of these (permanganate, RuO4), which the recorded reactions do not show often enough to list them',
       'high-valent metal oxo and iodine(V) reagents take hydrides or electrons from whatever is most oxidisable', 'Clayden, Organic Chemistry (oxidation)')
advice('cleaved-by-acid', 'a strong acid (TFA, HCl, HBr, H2SO4, TsOH, MsOH, TfOH, formic acid) + an acid-labile protecting group', 'protecting group', 'cleaved', T_STRONG_ACID, ['F078', 'F039', 'F052', 'F134', 'F132', 'F126', 'F200', 'F201'],
       'removed or opened: Boc and tert-butyl esters / ethers (TFA, HCl), acetals, THP and MOM ethers, trityl, PMB, enol ethers and silyl ethers are cleaved by acid; benzyl carbamates (Cbz) survive TFA and HCl but go with HBr / AcOH; choose a base-mediated step or a milder acid when the group must survive',
       'the cation formed on protonation is stabilised (tert-butyl, oxocarbenium, trityl, benzylic)', 'Greene, Protective Groups in Organic Synthesis', show='both')
advice('cleaved-by-boron-halide', 'BBr3 / BCl3 / TMSI / TMSBr / HBr + an ether, benzyl group or acid-labile group', 'protecting group', 'cleaved', T_ETHER_CLEAVING, ['F050', 'F051', 'F024', 'F130', 'F131', 'F078', 'F039', 'F052', 'F103'],
       'cleaved: methyl and benzyl ethers (BBr3 and BCl3 at −78 °C to rt), aryl alkyl ethers, benzyl carbamates and esters, tert-butyl groups, acetals, silyl ethers and phosphonate esters (TMSBr)',
       'the boron halide binds the ether O, and bromide displaces the alkyl group', 'Greene, Protective Groups in Organic Synthesis', show='both')
advice('cleaved-by-base-nucleophile', 'hydroxide or an alkoxide + a base-labile group', 'protecting group', 'cleaved', T_HYDROXIDE,
       ['F004', 'F139', 'F138', 'F136', 'F135', 'F124', 'F079', 'F074', 'F125', 'F023', 'F041', 'F020', 'F121', 'F080', 'F069', 'F065', 'F066', 'F064'],
       'attacked: esters hydrolyse or transesterify (with an alkoxide the alkoxy group of the ester exchanges), Fmoc and trifluoroacetyl come off, carbonates, thioesters and active esters are cleaved, α-halo carbonyls and epoxides are opened, sulfonates are displaced',
       'hydroxide and alkoxides are strong nucleophiles toward carbonyl carbons; Fmoc is removed by E1cB (fluorenyl C–H pKa about 23)', 'Greene, Protective Groups in Organic Synthesis')
advice('cleaved-by-Pd0', 'a palladium catalyst + an allyl ester, carbamate or ether', 'protecting group', 'cleaved', T_PD, ['F133'],
       'cleaved: Pd(0) ionises allyl carbonates, carbamates (Alloc), esters and ethers into a π-allyl complex that a nucleophile (amine, morpholine, dimedone, silane) takes; allyl groups are otherwise stable to acid and base',
       'oxidative addition of Pd(0) into the allylic C–O bond', 'Greene, Protective Groups in Organic Synthesis')
advice('also-hydrogenated', 'H2 with Pd / Pt / Rh / Ni + other reducible groups', 'incompatible group', 'also_reacts', T_H2_PGM,
       ['F005', 'F006', 'F007', 'F008', 'F030', 'F130', 'F131', 'F132', 'F025', 'F056', 'F055', 'F013', 'F084', 'F047', 'F086', 'F170', 'F175', 'F036', 'F174', 'F022', 'F191', 'F031', 'F049', 'F196', 'F204', 'F018', 'F117', 'F057'],
       'also reduced or cleaved: alkenes and alkynes, benzyl ethers / benzyl esters / Cbz / PMB (hydrogenolysis), nitro groups and azides, aryl iodides, bromides and (slowly) chlorides, benzylic halides, nitriles, imines and oximes',
       'a metal surface binds π systems and activates H2; weak benzylic C–O / C–N bonds and weak C–X bonds are cleaved', 'Greene, Protective Groups in Organic Synthesis; Clayden, Organic Chemistry')
advice('desulfurised-by-Raney-Ni', 'Raney nickel + sulfur', 'incompatible group', 'also_reacts', [['smiles:[Ni]']], ['F072', 'F073', 'F076', 'F104', 'F074'],
       'desulfurised: Raney Ni replaces C–S bonds by C–H (thioethers, thiols, thiophenes, thioacetals), whatever else is being hydrogenated',
       'nickel binds sulfur strongly and hydrogenolyses the C–S bond', 'Clayden, Organic Chemistry')
advice('poisoned-by-sulfur', 'hydrogenation over Pd / Pt / Rh / Ru / Ir + sulfur', 'catalyst poison', 'poisons', T_H2_NOT_NI, ['F072', 'F073', 'F076', 'F074', 'F075', 'F104'],
       'the catalyst is poisoned or inhibited: thiols, sulfides, thioureas and thiophenes bind Pd and Pt irreversibly (use much more catalyst, a sulfur-tolerant catalyst, or do the step before the sulfur is introduced); cross-couplings of thiophenes and aryl sulfides usually work with a suitable ligand',
       'soft sulfur lone pairs bind the soft metal more strongly than the substrate', 'Clayden, Organic Chemistry; Hartwig, Organotransition Metal Chemistry')
advice('inhibited-by-chelating-N', 'a palladium or copper catalyst + a chelating nitrogen group', 'catalyst poison', 'poisons', T_PGM + [['frag:F165']], ['F096', 'F097', 'F098', 'F099', 'F114', 'F115', 'F182', 'F183'],
       'hydrazines, hydroxylamines, amidines, tetrazoles and 2-aminoazines chelate the metal and slow couplings (use a bidentate or bulky ligand, protect the nitrogen, or add a Lewis acid to tie it up)',
       'a bidentate N,N group competes with the substrate and the ligand for the metal', 'Hartwig, Organotransition Metal Chemistry')
advice('consumed-lewis-acid', 'a Lewis acid (AlCl3, TiCl4, BF3, ZnCl2, BBr3) + a basic nitrogen or oxygen', 'catalyst poison', 'poisons', T_LEWIS,
       ['F002', 'F034', 'F093', 'F026', 'F094', 'F091', 'F112', 'F113', 'F115', 'F021', 'F184', 'F095', 'F004', 'F019', 'F018', 'F024', 'F050', 'F028', 'F003', 'F014', 'F015', 'F001', 'F118'],
       'the Lewis acid is bound and used up by the basic site (an amine, azine or amide takes one equivalent), so the group it was meant to activate stays untouched: Friedel–Crafts fails on anilines and pyridines; use one extra equivalent per basic site',
       'a free lone pair coordinates the Lewis acid; the complexed group is then a strong electron-withdrawing group', 'Clayden, Organic Chemistry (Friedel–Crafts); Directing_Groups: blocks Friedel–Crafts', min_distinct=1)
advice('competing-nucleophile-acylation', 'an acylating, sulfonylating or alkylating agent + several nucleophilic sites', 'competing site', 'competes', T_ACYLATING + T_ALKYLATING, NUCLEOPHILES,
       'which site reacts decides the product: amines are acylated before alcohols before phenols; with a base, alkylation follows acidity (phenol O before alcohol O) and softness (S and N before O with soft electrophiles; hard electrophiles go to O); ambident anions (pyridone, amide, azole) give N- and O-products',
       'nucleophilicity falls with delocalisation of the lone pair and rises with basicity and polarisability; with a base the most acidic site is the one deprotonated', 'Clayden, Organic Chemistry; Reactivity_Order sheet', min_distinct=1)
advice('competing-radical-or-ring-bromination', 'NBS or Br2 + an electron-rich ring, alkene or benzylic C–H', 'competing site', 'competes', T_HALOGENATING, ['F035', 'F037', 'F038', 'F005', 'F006', 'F028', 'F026', 'F050', 'F106', 'F105', 'F104'],
       'radical bromination goes at the weakest C–H (allylic, benzylic) but the same reagent adds to alkenes and brominates electron-rich rings by electrophilic substitution; light or a radical initiator in CCl4 / PhCF3 favours the radical path, polar solvents and acid the ring',
       'conditions decide between a radical chain and electrophilic bromine', 'Clayden, Organic Chemistry; Condition_Suggestions sheet')
advice('competing-diazonium-coupling', 'a diazonium salt (or nitrite making one) + a phenol, aniline or active methylene', 'competing site', 'competes', [['frag:F027'], ['smiles:O=N[O-]'], ['smiles:CC(C)(C)ON=O'], ['smiles:CC(C)CCON=O']], ['F028', 'F026', 'F002', 'F034', 'F102', 'F106', 'F105'],
       'azo coupling: the diazonium ion attacks phenols (pH 8–10, ortho / para to O⁻), anilines (pH 4–6) and active methylene compounds (Japp–Klingemann) before it can be used for Sandmeyer or reduction; keep it at 0–5 °C and add the nucleophile on purpose',
       'the diazonium ion is a weak electrophile that picks the most electron-rich partner', 'Clayden, Organic Chemistry (diazonium salts); Fragment_Conditions F027')
advice('activates-every-OH', 'SOCl2, POCl3, PCl5, oxalyl chloride or Tf2O + several O–H / N–H groups', 'competing site', 'also_reacts', T_DEHYDRATING,
       ['F001', 'F118', 'F003', 'F014', 'F015', 'F016', 'F017', 'F129', 'F197', 'F198', 'F199', 'F028', 'F120', 'F089', 'F049'],
       'every hydroxyl group, acid, primary amide (to nitrile), azole N–H and oxime (Beckmann) can react: acids become acyl chlorides, alcohols chlorides or sulfonates, primary amides dehydrate; amines are protonated or sulfonylated',
       'these reagents turn any O–H into a leaving group', 'Clayden, Organic Chemistry')
advice('cleaved-by-fluoride', 'fluoride (TBAF, CsF, HF) + silyl ethers', 'protecting group', 'cleaved', [['frag:F154'], ['smiles:F']], ['F052', 'F053'],
       'silyl ethers (TMS, TES, TBS, TIPS, TBDPS) and C-silyl groups are cleaved, fastest for the least hindered (TMS ≫ TES > TBS > TIPS ≈ TBDPS); the alkoxide formed can act as a base and move acyl or silyl groups',
       'the Si–F bond (about 140 kcal/mol) is far stronger than Si–O', 'Greene, Protective Groups in Organic Synthesis', show='both')
advice('also-attacked-by-organometallic', 'a Grignard or organolithium + an ester, nitrile, amide, acyl chloride, epoxide or carbonyl elsewhere in the molecule', 'incompatible group', 'also_reacts', T_ORGANOMETALLIC,
       ['F004', 'F022', 'F191', 'F021', 'F184', 'F123', 'F020', 'F041', 'F023', 'F018', 'F117', 'F019', 'F116', 'F040', 'F031', 'F124', 'F136', 'F121'],
       'also attacked: esters add twice (to the tertiary alcohol), acyl chlorides over-add, nitriles give ketones after hydrolysis, aldehydes and ketones add, epoxides open, enones take 1,2- or 1,4-addition; protect the group, use a more selective organometallic (organozinc, organocuprate), or use an acyl electrophile that stops at the ketone (Weinreb amide). In the recorded reactions an ester is consumed about as often with a Grignard as with other reagents: chemists choose conditions that spare it, or it is the target',
       'carbon nucleophiles of this strength attack every electrophilic carbonyl-type carbon', 'Clayden, Organic Chemistry (organometallics)')
advice('halogen-metal-exchange', 'n-BuLi, s-BuLi, t-BuLi or iPrMgCl / iPrMgBr + an aryl bromide / iodide or a vinyl halide', 'incompatible group', 'also_reacts', T_EXCHANGE, ['F084', 'F013', 'F087', 'F012', 'F171'],
       'halogen–metal exchange: n-BuLi, t-BuLi and iPrMgCl swap with aryl iodides and bromides (and vinyl halides) at −78 °C faster than they add to an ester or nitrile or deprotonate; aldehydes and ketones are still attacked first; the new organometallic then reacts with whatever electrophile is present (aryl chlorides do not exchange)',
       'the exchange is fast and favours the more stabilised carbanion (sp2 over sp3)', 'Clayden, Organic Chemistry (organolithiums)')
advice('epimerised-by-base', 'a strong base + an α-amino acid derivative', 'stereochemistry', 'epimerises', T_STRONG_BASE, ['F119'],
       'epimerisation: a base that can remove the α C–H racemises or epimerises the stereocentre through the enolate (amino acid esters and amides first; α-aryl ketones and esters too); use a hindered or weaker base, a lower temperature or a shorter time. Knowledge only: the page does not test it, because the fragments cannot tell whether the α carbon is a stereocentre (F119 finds the free acid, which a base turns into a carboxylate that does not enolise)',
       'enolisation is fast when the α C–H pKa is below the pKaH of the base', _CL, page=False)
advice('protic-solvent-slows-SN2', 'a protic solvent (water, alcohols) with an anionic nucleophile', 'medium', 'slows', [['solvent:protic', 'frag:F153'], ['solvent:protic', 'frag:F155'], ['solvent:protic', 'frag:F061'], ['solvent:protic', 'frag:F060']], ['F010', 'F011', 'F012', 'F043', 'F044', 'F045'],
       'SN2 with anionic nucleophiles (halide, azide, cyanide, hydroxide, alkoxide) is slower in protic solvents (hydrogen-bonded nucleophile: 10^3–10^5 slower than in DMF, DMSO or acetone); SN1 and solvolysis are favoured',
       'protic solvents solvate small anions strongly', 'Clayden, Organic Chemistry (SN1 / SN2)', needs='solvent', show='both')
advice('water-destroys-reagent', 'water or a protic solvent with an acyl chloride, anhydride, organometallic or hydride', 'medium', 'quenches', [['solvent:protic']], ['F029', 'F140', 'F141', 'F042', 'F143', 'F145', 'F020', 'F041', 'F124'],
       'hydrolysed or quenched: acyl chlorides and anhydrides are hydrolysed, Grignards, organolithiums and hydrides destroyed; use an anhydrous aprotic solvent (THF, ether, DCM) under inert gas',
       'these reagents react faster with water than with most organic substrates', 'Fragment_Conditions (F020, F029, F042)', needs='solvent', group_where='both', show='both')
advice('leaving-group-order', 'several leaving groups in one molecule', 'leaving group', 'competes', [], LEAVING_GROUPS,
       'which one goes first: Pd / Ni / Cu couplings follow I > Br ≈ OTf > Cl ≫ F (oxidative addition); SN2 follows OTf ≫ I > Br ≈ OTs > Cl ≫ F; SNAr, when addition is rate-limiting, follows F > Cl ≈ Br > I and is helped by ortho / para electron-withdrawing groups and ring nitrogens. With several leaving groups in one molecule the best one for that mechanism reacts first',
       'bond strength and the stability of the leaving anion (pKa of HX: HI −10, HBr −9, HCl −7, HF 3.2)', 'Clayden, Organic Chemistry', show='both', min_distinct=2)
advice('weak-base-for-the-acid', 'a base weaker than the proton it should remove', 'reagent strength', 'slows', [['role:base']], ['F101', 'F102', 'F003', 'F014', 'F028', 'F001', 'F008', 'F095', 'F192', 'F193'],
       'little or no reaction when the pKaH of the base is more than 3–5 units below the pKa of the proton: K2CO3 (pKaH 10.3) deprotonates phenols and acids but not alcohols or ketones (use NaH, KOt-Bu or LDA); Et3N (10.8) does not deprotonate amides',
       'ΔpKa sets the equilibrium (K = 10^ΔpKa); a catalytic amount of the anion is enough when the next step is fast', 'Acids_Bases sheet; Evans pKa table', page=False)
advice('temperature-selectivity', 'temperature, time and equivalents', 'medium', 'redirects', [], [],
       'kinetic against thermodynamic control: low temperature and short times give the faster product (the less substituted enolate, the 1,2-adduct, ortho substitution), heat and long times the more stable one (the more substituted enolate, the 1,4-adduct, para, rearranged); equivalents decide mono- against bis-reaction',
       'the recorded reactions list the solvent only; temperature, time and amounts are not in the learned data, so this interference is advisory', 'Clayden, Organic Chemistry', needs='conditions', page=False)
advice('data-support-size', 'starting materials with fewer than 8 heavy atoms', 'data support', 'unreliable', [], [],
       'the learned predictor was built from patent reactions and has seen almost none this small, so no figure can say how far to trust it there and the page predicts nothing. {SIZE_FIGURE}',
       'a model can only be as good as the molecules it saw', 'tools/eval_top1.json by_size', needs='size', page=False)
advice('data-support-reagent-role', 'a reagent with no reactive role in the Reagents sheet', 'data support', 'unreliable', [], [],
       'spectator ions, solvents and parts the sheet does not list tell the model nothing about what the reagent does, so the page asks for a higher stated chance of being right before it shows an answer. {REAGENT_FIGURE}',
       'measured on held-out patents', 'tools/eval_top1.json by_reagent', needs='reagent', page=False)
ADVISORY_BY_KEY = {a['key']: a for a in ADVISORY}
assert len(ADVISORY_BY_KEY) == len(ADVISORY), 'duplicate advisory keys'
