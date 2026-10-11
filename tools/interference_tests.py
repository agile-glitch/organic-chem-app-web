"""Positive and negative examples for every site pattern of tools/interference_defs.py, and for the advisories.

    python tools/interference_tests.py            # run them (exit status 1 on a failure)
    python tools/interference_tests.py --json out.json     # also write the cases and what the Python gives, for js/interference.js to run the same ones (tools/int_check_js.py)

Each site case is (SMILES with the isotope-900 marks on the reacting atoms, flags that must be present, flags that must be absent).
The ring-substituent cases are generated: every substituent in each of ortho, meta and para position of a marked aromatic carbon (present), the other two
positions absent, and an unsubstituted ring with none. Many hand-written cases come from an independent review of the patterns (fused rings, ring heteroatoms,
charged atoms, groups the reaction itself created ...).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import interference_defs as D
import interferences as I


def ring_cases():
    out = []
    for label, key, _sm, _ss, _sp, kind, _why, ex in D.SUBSTITUENTS:
        for tag, rel in D.RELATIONS:
            pos = {'o': 1, 'm': 2, 'p': 3}[tag]
            atoms = ['[900cH]', 'c', 'c', 'c', 'c', 'c']
            atoms[pos] = f'c({ex})'
            smi = atoms[0] + '9' + ''.join(atoms[1:]) + '9'
            others = [f'ring-{r2}-{key}' for t2, r2 in D.RELATIONS if t2 != tag]
            out.append((smi, [f'ring-{rel}-{key}'], others))
    plain = '[900cH]9ccccc9'
    out.append((plain, [], [f'ring-{r}-{k}' for k in D.SUBST_KEYS for _, r in D.RELATIONS]))
    return out


def ring_hetero_cases():
    out = []
    def ring(k_n, element='n', anchor=0):
        atoms = ['c'] * 6
        atoms[anchor] = '[900cH]'
        for k in ([k_n] if isinstance(k_n, int) else k_n): atoms[k] = element
        return atoms[0] + '9' + ''.join(atoms[1:]) + '9'
    out.append((ring(1), ['ring-alpha-N'], ['ring-beta-N', 'ring-gamma-N', 'ring-between-2N']))
    out.append((ring(2), ['ring-beta-N'], ['ring-alpha-N', 'ring-gamma-N', 'ring-between-2N']))
    out.append((ring(3), ['ring-gamma-N'], ['ring-alpha-N', 'ring-beta-N', 'ring-between-2N']))
    out.append((ring([1, 5]), ['ring-alpha-N', 'ring-between-2N'], ['ring-gamma-N']))
    out.append(('[900cH]1ccc[s]1', ['ring-alpha-S', 'ring-5-membered'], ['ring-beta-SO']))
    out.append(('[900cH]1cc[s]c1', ['ring-beta-SO', 'ring-5-membered'], ['ring-alpha-S']))
    out.append(('[900cH]1ccc[o]1', ['ring-alpha-O', 'ring-5-membered'], ['ring-beta-SO']))
    out.append(('[900cH]1ccc[nH]1', ['ring-alpha-NR', 'ring-5-membered'], ['ring-alpha-N', 'ring-beta-NR']))
    out.append(('Cn1[900cH]ccc1', ['ring-alpha-NR', 'ring-5-membered'], ['ring-alpha-N']))
    out.append(('Cn1[900cH]nc2ccccc12', ['ring-alpha-NR', 'ring-between-2N'], ['ring-5-membered']))
    out.append(('[900cH]1c[nH]c2ccccc12', ['ring-beta-NR', 'ring-5-membered'], ['ring-alpha-NR', 'ring-peri']))
    out.append(('[900cH]1c[nH]c2cccc([N+](=O)[O-])c12', ['ring-beta-NR', 'ring-peri'], ['ring-meta-nitro', 'ring-para-nitro', 'ring-ortho-nitro']))
    out.append(('[900cH]1cccc2cccc([N+](=O)[O-])c12', ['ring-peri', 'ring-fused-alpha'], ['ring-meta-nitro', 'ring-para-nitro', 'ring-ortho-nitro']))
    out.append(('[900cH]1cccc2ccc([N+](=O)[O-])cc12', ['ring-fused-alpha'], ['ring-peri', 'ring-meta-nitro', 'ring-para-nitro']))
    out.append(('[900cH]1cc([N+](=O)[O-])cc2ccccc12', ['ring-fused-alpha', 'ring-meta-nitro'], ['ring-peri', 'ring-para-nitro']))
    out.append(('[900cH]1cccc2ncccc12', ['ring-fused-alpha'], ['ring-gamma-N', 'ring-beta-N', 'ring-alpha-N']))
    out.append(('[900cH]1cccc(C#N)n1', ['ring-alpha-N', 'ring-meta-cyano'], ['ring-ortho-cyano', 'ring-para-cyano']))
    out.append(('[900cH]1c(C#N)cccn1', ['ring-alpha-N', 'ring-ortho-cyano'], ['ring-meta-cyano']))
    out.append(('[O-][n+]1[900cH]cccc1', ['ring-alpha-azinium'], ['ring-alpha-NR', 'ring-alpha-N']))
    out.append(('C[n+]1cc[900cH]cc1', ['ring-gamma-azinium'], ['ring-alpha-azinium']))
    out.append(('[900cH]1cncs1', ['ring-alpha-S'], ['ring-5-membered']))
    out.append(('[900cH]1coc2ccccc12', ['ring-5-membered', 'ring-beta-SO'], ['ring-alpha-O']))
    out.append(('[900cH]1ccco1', ['ring-alpha-O', 'ring-5-membered'], []))
    out.append(('c1ccc2[900cH]cccc2c1', ['ring-fused-alpha'], []))
    out.append(('[900cH]1ccc2ccccc2c1', [], ['ring-fused-alpha']))
    out.append(('[900cH]1ccc2ccccc2n1', ['ring-alpha-N'], ['ring-fused-alpha']))
    return out


def hand_cases():
    S = []
    add = lambda smi, yes, no=(): S.append((smi, list(yes), list(no)))
    # ring substituents: groups the reaction itself created are not substituents; the groups' own definitions
    add('[900cH]9c([900NH2])cccc9', [], ['ring-ortho-nh'])
    add('[900cH]9c(N(c1ccccc1)C)cccc9', ['ring-ortho-1'], ['ring-ortho-nr2', 'ring-ortho-nh'])
    add('[900cH]9c(NC(=S)N)cccc9', [], ['ring-ortho-nh'])
    add('[900cH]9c(OC(F)(F)F)cccc9', ['ring-ortho-ocf3'], ['ring-ortho-or'])
    add('[900cH]9c(SC(F)(F)F)cccc9', [], ['ring-ortho-sr'])
    add('[900cH]9c(CCl)cccc9', ['ring-ortho-1'], ['ring-ortho-alkyl', 'ring-ortho-methyl'])
    add('[900cH]9c(C(F)(F)F)cccc9', ['ring-ortho-cf3'], ['ring-ortho-alkyl'])
    add('[900cH]9c(OS(C)(=O)=O)cccc9', ['ring-ortho-osulfonyl'], ['ring-ortho-sulfonyl', 'ring-ortho-or'])
    add('[900cH]9c(NS(C)(=O)=O)cccc9', ['ring-ortho-nhso2'], ['ring-ortho-nh', 'ring-ortho-sulfonyl'])
    add('[900cH]9c(OC(C)=O)cccc9', ['ring-ortho-ocor'], ['ring-ortho-or', 'ring-ortho-acyl'])
    add('[900cH]1c(C)cccc1C', ['ring-ortho-1', 'ring-ortho-2'])
    add('[900cH]1cccc2c1CCCC2', ['ring-ortho-1'], ['ring-ortho-2'])
    add('[900cH]9ccccc9', [], ['ring-ortho-1', 'ring-ortho-2', 'ring-ortho-branched'])
    add('[900cH]1c(C)cccc1', ['ring-ortho-1', 'ring-ortho-methyl'], ['ring-ortho-2', 'ring-ortho-branched'])
    add('[900cH]1c(C(C)C)cccc1', ['ring-ortho-1', 'ring-ortho-branched'])
    add('[900cH]1c(-c2ccccc2)cccc1', ['ring-ortho-aryl', 'ring-ortho-1'])
    # azole / azine nitrogen
    add('C[900n]1ccnc1', ['azole-n-1,3'], ['azole-n-next-to-N', 'azole-n-hindered'])
    add('C[900n]1nccc1', ['azole-n-next-to-N'], ['azole-n-1,3', 'azole-n-hindered'])
    add('C[900n]1nc(C)cc1C', ['azole-n-next-to-N', 'azole-n-hindered'], ['azole-n-hindered-2'])
    add('C[900n]1c(C)ccc1C', ['azole-n-hindered', 'azole-n-hindered-2'])
    add('C[900n]1ccc2ccccc12', ['azole-n-fused'])
    add('C[900n+]1ccccc1', ['azine-n'], ['azole-n-hindered', 'azine-n-hindered'])
    add('O=c1cccc[900nH]1', ['azine-n'])
    add('[O-][900n+]1ccccc1', ['azine-n'])
    add('Cc1cccc(C)[900n+]1[O-]', ['azine-n-hindered', 'azine-n-hindered-2'])
    add('Cc1cccc[900n+]1[O-]', ['azine-n-hindered'], ['azine-n-hindered-2'])
    # sp3 carbon
    add('CC[900CH2][900NH]C(=O)c1ccccc1', ['sp3-primary', 'sp3-primary-at-N'], ['sp3-methyl', 'sp3-secondary', 'sp3-tertiary', 'sp3-methyl-at-N'])
    add('C[900CH](C)[900NH]C(=O)C', ['sp3-secondary', 'sp3-secondary-at-N'], ['sp3-primary', 'sp3-tertiary'])
    add('C[900C](C)(C)[900Br]', ['sp3-tertiary', 'sp3-tertiary-at-X'], ['sp3-secondary'])
    add('[900CH3][900NH]CC', ['sp3-methyl', 'sp3-methyl-at-N'], ['sp3-primary'])
    add('c1ccccc1[900CH2][900Br]', ['sp3-benzylic', 'sp3-primary', 'sp3-primary-at-X'], ['sp3-allylic', 'sp3-alpha-carbonyl'])
    add('C=C[900CH2][900Br]', ['sp3-allylic'], ['sp3-benzylic'])
    add('C#C[900CH2][900Br]', ['sp3-propargylic'], ['sp3-allylic'])
    add('N#C[900CH2][900Br]', ['sp3-alpha-nitrile'], ['sp3-alpha-carbonyl', 'sp3-propargylic'])
    add('CC(=O)[900CH2][900Br]', ['sp3-alpha-carbonyl'], ['sp3-allylic'])
    add('FC(F)(F)[900CH2][900NH]C', ['sp3-alpha-CF3'], ['sp3-neopentyl'])
    add('Cl[900CH](C)[900NH]C', ['sp3-alpha-halogen'])
    add('CO[900CH](C)[900Br]', ['sp3-alpha-oxygen'], ['sp3-alpha-nitrogen'])
    add('CN(C)[900CH](C)[900Br]', ['sp3-alpha-nitrogen'], ['sp3-alpha-oxygen'])
    add('[O-][N+](=O)[900CH2][900Br]', [], ['sp3-alpha-nitrogen'])
    add('CC(C)(C)[900CH2][900Br]', ['sp3-neopentyl', 'sp3-primary'], ['sp3-beta-branched'])
    add('CC(C)[900CH2][900Br]', ['sp3-beta-branched', 'sp3-primary'], ['sp3-neopentyl'])
    add('C[CH]([900OH])[900CH2][900NH]C', ['sp3-beta-opened'], [])
    add('C1[900O][900CH]1C', ['sp3-ring3'], ['sp3-ring4', 'sp3-ring5'])
    add('C1CO[900CH]1C', ['sp3-ring4'], ['sp3-ring3', 'sp3-ring5'])
    add('C1CCO[900CH]1C', ['sp3-ring5'], ['sp3-ring3', 'sp3-ring4'])
    add('C1CC[900CH](C1)[900O]C', ['sp3-ring5', 'sp3-secondary-at-O'], ['sp3-ring6'])
    add('C1CCC[900CH](C1)[900O]C', ['sp3-ring6'], ['sp3-ring5'])
    add('CC[900CH]([900OH])C', ['sp3-secondary-at-O'], ['sp3-primary-at-O'])
    add('CC[900CH2][900OH]', ['sp3-primary-at-O'], ['sp3-secondary-at-O', 'sp3-methyl-at-O'])
    add('c1ccccc1[900CH]([900Br])C', ['sp3-secondary-at-X', 'sp3-benzylic'])
    add('C[900CH2][900c]1ccccc1', ['sp3-primary-at-C'])
    add('C[900CH2][900B](O)O', ['sp3-primary-at-M'])
    # carbonyl carbon
    add('CC(C)[900C](=O)[900NH]C', ['co-alpha-branched'], ['co-aryl', 'co-alpha-quaternary'])
    add('CC(C)(C)[900C](=O)[900NH]C', ['co-alpha-quaternary'], ['co-alpha-branched'])
    add('CC(C)(N)[900C](=O)[900NH]C', ['co-alpha-quaternary'])
    add('CC(C)(O)[900C](=O)[900NH]C', ['co-alpha-quaternary'])
    add('FC(F)(F)[900C](=O)[900NH]C', ['co-alpha-CF3'], ['co-alpha-quaternary', 'co-alpha-branched'])
    add('ClC(Cl)(Cl)[900C](=O)[900NH]C', [], ['co-alpha-quaternary'])
    add('FC(F)[900C](=O)[900NH]C', [], ['co-alpha-branched'])
    add('O=[900C](c1ccccc1)[900NH]C', ['co-aryl'], ['co-vinyl'])
    add('C=C[900C](=O)[900NH]C', ['co-vinyl'], ['co-aryl'])
    add('ClC[900C](=O)[900NH]C', ['co-alpha-hetero'], ['co-alpha-CF3'])
    add('O=[900C](OC)[900NH]C', ['co-two-hetero'])
    add('O=[900C](c1ccccn1)[900NH]C', ['co-heteroaryl', 'co-aryl'])
    add('O=[900C]1CCC[900NH]1', ['co-ring5'], ['co-ring4', 'co-ring6'])
    add('O=[900C]1CCCC[900NH]1', ['co-ring6'], ['co-ring5'])
    add('O=[900C]1CC[900N]1C', ['co-ring4'], ['co-ring5'])
    # N, O, S
    add('c1ccccc1[900NH][900C](C)=O', ['n-aryl', 'n-primary-origin'], ['n-secondary-origin', 'n-heteroaryl', 'n-benzylic'])
    add('c1cccnc1[900NH][900C](C)=O', ['n-heteroaryl'])
    add('CC(C)[900NH][900C](C)=O', ['n-alpha-branched', 'n-primary-origin'], ['n-aryl'])
    add('c1ccccc1C[900NH][900C](C)=O', ['n-benzylic', 'n-primary-origin'], ['n-aryl'])
    add('FC(F)(F)C[900NH][900C](C)=O', ['n-alpha-ewg'])
    add('C1CCC[900N](C1)[900C](C)=O', ['n-in-ring', 'n-secondary-origin'])
    add('C[900N](CC)[900C](C)=O', ['n-secondary-origin'], ['n-primary-origin'])
    add('C[900N+](C)(C)[900CH2]c1ccccc1', ['n-tertiary-origin'])
    add('CC(=O)[900N](C(C)=O)[900C]', ['n-acyl-origin'])
    add('CN[900NH][900C](C)=O', ['n-adjacent-hetero'])
    add('C1C[900N]1[900C](C)=O', ['n-ring-small', 'n-in-ring'], ['n-ring-medium', 'n-ring-macro'])
    add('C1CCCCCCCCCC[900N]1[900C](C)=O', ['n-ring-macro', 'n-in-ring'], ['n-ring-medium', 'n-ring-small'])
    add('C1CCCCCCCCCC[900O]1', ['o-ring-macro', 'o-in-ring'], ['o-ring-medium'])
    add('C1CCCCC[900N]1[900C](C)=O', ['n-ring-medium', 'n-in-ring'], ['n-ring-small'])
    add('C1CC[900O]1', ['o-ring-small', 'o-in-ring'])
    add('c1ccccc1[900O][900C](C)=O', ['o-aryl'], ['o-alpha-branched', 'o-primary-origin'])
    add('CC[900O][900C](C)=O', ['o-primary-origin'], ['o-aryl'])
    add('C[900O][900C](C)=O', ['o-primary-origin'])
    add('CC(C)[900O][900C](C)=O', ['o-alpha-branched'], ['o-aryl', 'o-primary-origin'])
    add('CC(C)(C)[900O][900C](C)=O', ['o-alpha-branched'], ['o-primary-origin'])
    add('c1ccccc1C[900O][900C](C)=O', ['o-benzylic', 'o-primary-origin'])
    add('FC(F)(F)C[900O][900C](C)=O', ['o-alpha-ewg', 'o-primary-origin'])
    add('C1CCC[900O]1', ['o-in-ring'], ['o-primary-origin'])
    add('c1ccccc1[900O]C', ['o-aryl'], ['o-primary-origin'])
    add('C[900O]c1ccccn1', ['o-heteroaryl', 'o-aryl'])
    add('c1ccccc1[900S][900C](C)=O', ['s-aryl'], ['s-alkyl'])
    add('CC[900S]c1ccccc1', ['s-aryl', 's-alkyl'])
    add('CC[900S][900C](C)=O', ['s-alkyl'], ['s-aryl'])
    add('Cc1cccc(C)c1[900NH][900C](C)=O', ['n-aryl', 'n-aryl-ortho2'])
    add('Cc1cccc(C)c1[900O][900C](C)=O', ['o-aryl', 'o-aryl-ortho2'])
    add('Cc1ccccc1[900NH][900C](C)=O', ['n-aryl'], ['n-aryl-ortho2'])
    # an aryl ring on the reacting atom: Hammett
    add('O=[N+]([O-])c1ccc([900NH][900C](C)=O)cc1', ['n-aryl-p-ewg'], ['n-aryl-p-edg', 'n-aryl-o-ewg'])
    add('CC[900NH][900c]1ccc([N+](=O)[O-])cc1', [], ['n-aryl-p-ewg', 'n-aryl'])
    add('[O-]C(=O)c1ccc([900NH][900C](C)=O)cc1', ['n-aryl'], ['n-aryl-p-ewg', 'n-aryl-p-edg'])
    add('FC(F)(F)Oc1ccc([900NH][900C](C)=O)cc1', ['n-aryl-p-ewg'], ['n-aryl-p-edg'])
    add('Clc1ccc([900NH][900C](C)=O)cc1', ['n-aryl'], ['n-aryl-p-ewg', 'n-aryl-p-edg'])
    add('COc1ccc([900NH][900C](C)=O)cc1', ['n-aryl-p-edg'], ['n-aryl-p-ewg'])
    add('O=[N+]([O-])c1ccccc1[900NH][900C](C)=O', ['n-aryl-o-ewg'], ['n-aryl-p-ewg'])
    add('Cc1ccccc1[900NH][900C](C)=O', ['n-aryl-o-edg'], ['n-aryl-p-edg', 'n-aryl-o-ewg'])
    add('O=[N+]([O-])c1ccc([900O][900C](C)=O)cc1', ['o-aryl-p-ewg'])
    add('COc1ccc([900O][900C](C)=O)cc1', ['o-aryl-p-edg'])
    add('FC(F)(F)c1ccc([900S][900CH3])cc1', ['s-aryl-p-ewg'])
    add('O=[900C](c1ccc([N+](=O)[O-])cc1)[900NH]C', ['co-aryl-p-ewg', 'co-aryl'])
    add('O=[900C](c1ccc(OC)cc1)[900NH]C', ['co-aryl-p-edg', 'co-aryl'])
    add('COc1ccc([900CH2][900Br])cc1', ['bz-aryl-p-edg', 'sp3-benzylic'], ['bz-aryl-p-ewg'])
    add('O=[N+]([O-])c1ccc([900CH2][900Br])cc1', ['bz-aryl-p-ewg', 'sp3-benzylic'], ['bz-aryl-p-edg'])
    add('C1CCCCC[900O]1', ['o-ring-medium', 'o-in-ring'], ['o-ring-small'])
    add('Cc1cccc(C)c1[900S]C', ['s-aryl-ortho2', 's-aryl'])
    add('O=[N+]([O-])c1ccccc1[900O][900C](C)=O', ['o-aryl-o-ewg'], ['o-aryl-p-ewg'])
    add('Cc1ccccc1[900O][900C](C)=O', ['o-aryl-o-edg'], ['o-aryl-p-edg'])
    add('O=[N+]([O-])c1ccccc1[900S]C', ['s-aryl-o-ewg'])
    add('COc1ccc([900S]C)cc1', ['s-aryl-p-edg'], ['s-aryl-p-ewg'])
    add('COc1ccccc1[900S]C', ['s-aryl-o-edg'])
    add('O=[900C](c1ccccc1[N+](=O)[O-])[900NH]C', ['co-aryl-o-ewg'])
    add('O=[900C](c1ccccc1OC)[900NH]C', ['co-aryl-o-edg'])
    add('O=[N+]([O-])c1ccccc1[900CH2][900Br]', ['bz-aryl-o-ewg'])
    add('COc1ccccc1[900CH2][900Br]', ['bz-aryl-o-edg'])
    # alkene
    add('[900CH2]=[900CH]C', ['alkene-terminal'])
    add('C[900CH]=[900CH]C(=O)C', ['alkene-conjugated-carbonyl'], ['alkene-terminal'])
    add('N#CC=[900CH]c1ccccc1', ['alkene-conjugated-carbonyl'], ['alkene-styrenyl'])
    add('C[900CH]=[900CH]c1ccccc1', ['alkene-styrenyl'], ['alkene-conjugated-carbonyl'])
    add('C[900CH]=[900CH]OC', ['alkene-electron-rich'])
    add('[900CH](C)=C[N+](=O)[O-]', ['alkene-conjugated-carbonyl'], ['alkene-electron-rich'])
    add('C[900CH]=[900CH]OC(C)=O', [], ['alkene-electron-rich'])
    add('C[900C](C)=[900CH]C', ['alkene-two-carbons'])
    add('C(=O)(OC)[900CH]=[900CH][900c]1ccccc1', ['alkene-new-bond-to-C', 'alkene-styrenyl', 'alkene-conjugated-carbonyl'])
    add('c1ccccc1[900CH]=[900CH]c1ccccc1', ['alkene-styrenyl'])
    # review round 3: fused rings, N-heteroaryl, terminal alkene, CF3 the reaction created, ring closure, deuterium
    add('[900cH]1cccc2OCOc12', [], ['ring-fused-alpha'])
    add('[900cH]1cccc2CCCCc12', [], ['ring-fused-alpha'])
    add('[900cH]1cccc2CCCc12', [], ['ring-fused-alpha'])
    add('[900cH]1cccc2NCCc12', [], ['ring-fused-alpha'])
    add('[900cH]1cccc2C(=O)c3ccccc3-c12', [], ['ring-fused-alpha'])
    add('[900cH]1cccc2ccccc12', ['ring-fused-alpha'])
    add('[900cH]1cccc2[nH]ccc12', ['ring-fused-alpha'])
    add('O=[900C](c1cc2ccccc2[nH]1)[900NH]C', [], ['co-heteroaryl'])
    add('O=[900C](c1ccc[nH]1)[900NH]C', [], ['co-heteroaryl'])
    add('O=[900C](c1ccccn1)[900NH]C', ['co-heteroaryl'])
    add('CC(=O)[900NH]c1ccc[nH]1', [], ['n-heteroaryl'])
    add('CC(=O)[900NH]c1ccccn1', ['n-heteroaryl'])
    add('C[900O]c1cc2ccccc2[nH]1', [], ['o-heteroaryl'])
    add('C[900O]c1ccccn1', ['o-heteroaryl'])
    add('C=[900C](C)C', ['alkene-two-carbons'], ['alkene-terminal'])
    add('C[900CH]=C', ['alkene-terminal'])
    add('C[900C]([900OH])(c1ccccc1)[900C]([900F])([900F])[900F]', [], ['sp3-alpha-CF3'])
    add('FC(F)(F)[900CH2][900Br]', ['sp3-alpha-CF3'])
    add('FC(F)(F)C[900CH2][900Br]', [], ['sp3-alpha-CF3'])
    add('O=[900C](C(F)(F)F)[900NH]C', ['co-alpha-CF3'], ['co-alpha-hetero'])
    add('O=[900C]([900C]([900F])([900F])[900F])c1ccccc1', [], ['co-alpha-CF3'])
    add('O=[900C](C(Cl)C)[900NH]C', ['co-alpha-hetero'])
    add('CC1[900CH2][900O]1', [], ['sp3-beta-opened'])
    add('CC1[900CH2][900NH]1', [], ['sp3-beta-opened'])
    add('C[CH]([900OH])[900CH2][900NH]C', ['sp3-beta-opened'])
    add('[2H]C([2H])([2H])[900NH][900C](C)=O', [], ['n-alpha-branched'])
    add('[2H]C([2H])([2H])[900O][900C](C)=O', [], ['o-alpha-branched'])
    add('[900cH]1c(C([2H])([2H])[2H])cccc1', ['ring-ortho-methyl'], ['ring-ortho-branched'])
    return S


def sp3_partner_cases():
    heads = {'methyl': '[900CH3]{p}', 'primary': 'C[900CH2]{p}', 'secondary': 'C[900CH](C){p}', 'tertiary': 'C[900C](C)(C){p}'}
    parts = {'O': '[900OH]', 'N': '[900NH2]', 'S': '[900SH]', 'X': '[900Br]', 'C': '[900CH3]', 'M': '[900B](O)O'}
    out = []
    for deg, h in heads.items():
        for pk, ps in parts.items():
            no = [f'sp3-{deg}-at-{q}' for q in parts if q != pk] + [f'sp3-{d2}-at-{pk}' for d2 in heads if d2 != deg and not (pk == 'C' and d2 == 'methyl')]
            out.append((h.format(p=ps), [f'sp3-{deg}-at-{pk}'], no))
    return out


def injectivity_cases():
    """A substituent's definition names a neighbour atom (the [#6] of a ketone, the c of an aryl ether, the [#6] of a thioether). That neighbour has to be a different
    atom from the ring carbon the substituent hangs on: inside a recursive $() the match is made alone, and a carboxylate passed as a ketone, a phenol as an aryl
    ether and a thiol as a thioether. Each probe gives its own ring flag and none of the flags of the group it could be mistaken for."""
    out = []
    probes = [  # (substituent SMILES, ring key it must give or None, ring keys it must not give, Hammett class of the aryl ring on an amide N: 'ewg', 'edg' or None)
        ('C(=O)[O-]', 'co2minus', ['acyl', 'carboxy', 'carboxamide', 'formyl'], None),
        ('C(=O)O', 'carboxy', ['acyl', 'co2minus', 'carboxamide'], 'ewg'),
        ('C(=O)OC', 'carboxy', ['acyl', 'carboxamide'], 'ewg'),
        ('C(=O)N', 'carboxamide', ['acyl', 'carboxy'], 'ewg'),
        ('C(=O)C', 'acyl', ['carboxy', 'carboxamide', 'formyl'], 'ewg'),
        ('C(=O)c8ccccc8', 'acyl', ['aryl', 'carboxy'], 'ewg'),
        ('C=O', 'formyl', ['acyl', 'vinyl'], 'ewg'),
        ('O', 'oh', ['oar', 'or', 'ocor', 'ocf3', 'osulfonyl'], 'edg'),
        ('OC', 'or', ['oar', 'oh'], 'edg'),
        ('Oc8ccccc8', 'oar', ['or', 'oh'], None),
        ('OC(=O)C', 'ocor', ['or', 'acyl'], 'ewg'),
        ('S', None, ['sr', 'oar', 'oh', 'nh'], None),
        ('SC', 'sr', ['oar', 'or'], None),
        ('Sc8ccccc8', 'sr', ['oar', 'aryl'], None),
        ('N', 'nh', ['nr2', 'nhcor', 'nhso2'], 'edg'),
        ('NC', 'nh', ['nr2'], 'edg'),
        ('NC(=O)C', 'nhcor', ['nh', 'acyl', 'carboxamide'], None),
        ('NS(C)(=O)=O', 'nhso2', ['nh', 'sulfonyl'], None),
        ('C#N', 'cyano', ['acyl', 'ethynyl'], 'ewg'),
        ('c8ccccc8', 'aryl', ['vinyl'], None),
        ('C=C', 'vinyl', ['aryl'], None),
    ]
    for x, yes, no, cls in probes:
        for tag, rel, tpl in (('p', 'para', '[900cH]9ccc({x})cc9'), ('m', 'meta', '[900cH]9cc({x})ccc9'), ('o', 'ortho', '[900cH]9c({x})cccc9')):
            out.append((tpl.format(x=x), [f'ring-{rel}-{yes}'] if yes else [], [f'ring-{rel}-{k}' for k in no]))
        nv = f'CC(=O)[900NH]c9ccc({x})cc9'
        out.append((nv, ['n-aryl'] + ([f'n-aryl-p-{cls}'] if cls else []), [f'n-aryl-p-{c}' for c in ('ewg', 'edg') if c != cls]))
    return out


def all_cases():
    return ring_cases() + ring_hetero_cases() + hand_cases() + sp3_partner_cases() + injectivity_cases()


def advisory_cases():
    """(name, ctx, present, consumed, keys that must be shown, keys that must not)"""
    def ctx(frags=(), parts=(), roles=(), solv=()): return {'frags': set(frags), 'parts': set(parts), 'roles': set(roles), 'solvent_types': set(solv)}
    return [
        ('Grignard + alcohol, aldehyde is the target', ctx(['F029']), {'F003', 'F014', 'F018'}, {'F018'}, ['quench-organometallic-acidic-XH'], ['also-attacked-by-organometallic']),
        ('Grignard typed as a cation + ester bystander', ctx(parts=['[CH3][Mg+]']), {'F004', 'F018'}, {'F018'}, ['also-attacked-by-organometallic'], []),
        ('LiAlH4 + nitro bystander', ctx(['F042']), {'F004', 'F025'}, {'F004'}, ['also-reduced-LiAlH4'], []),
        ('LiAlH4 typed as hydride + Al3+', ctx(['F143'], ['[Al+3]', '[Li+]']), {'F004', 'F025'}, {'F004'}, ['also-reduced-LiAlH4'], []),
        ('NaBH4 leaves the amide', ctx(parts=['[BH4-]', '[Na+]']), {'F021', 'F018'}, {'F018'}, ['spared-by-NaBH4'], ['also-reduced-LiBH4']),
        ('NaBH4 with LiCl is LiBH4 in situ', ctx(parts=['[BH4-]', '[Na+]', '[Li+]', '[Cl-]']), {'F004', 'F021', 'F018'}, {'F018'}, ['also-reduced-LiBH4'], ['spared-by-NaBH4']),
        ('LiBH4 reduces the ester too', ctx(parts=['[BH4-]', '[Li+]']), {'F004', 'F018'}, {'F018'}, ['also-reduced-LiBH4'], ['spared-by-NaBH4']),
        ('TFA and Boc', ctx(parts=['O=C(O)C(F)(F)F']), {'F078'}, {'F078'}, ['cleaved-by-acid'], []),
        ('TFA leaves Cbz alone', ctx(parts=['O=C(O)C(F)(F)F']), {'F131', 'F078'}, {'F078'}, ['cleaved-by-acid'], ['cleaved-by-boron-halide']),
        ('HBr cleaves Cbz', ctx(parts=['Br']), {'F131'}, set(), ['cleaved-by-boron-halide'], []),
        ('triethylamine does not cleave esters', ctx(parts=['CCN(CC)CC'], roles=['base']), {'F004'}, set(), [], ['cleaved-by-base-nucleophile', 'cleaved-by-acid', 'epimerised-by-base']),
        ('lithium hydroxide, an ester is the target and a trifluoroacetamide is also there', ctx(parts=['[Li][OH]']), {'F004', 'F136'}, {'F004'}, ['cleaved-by-base-nucleophile'], []),
        ('lithium hydroxide hydrolyses the only ester: nothing else to say', ctx(parts=['[Li][OH]']), {'F004'}, {'F004'}, [], ['cleaved-by-base-nucleophile']),
        ('epimerisation is knowledge only (the fragments cannot see a stereocentre)', ctx(['F142']), {'F119'}, set(), [], ['epimerised-by-base']),
        ('LDA and ethyl acetate: no stereocentre', ctx(['F142']), {'F101', 'F193'}, set(), [], ['epimerised-by-base']),
        ('Pd/H2 and a benzyl ether, an alkene is also there', ctx(['F164'], ['[H][H]']), {'F130', 'F005'}, {'F130'}, ['also-hydrogenated'], []),
        ('Pd coupling without H2', ctx(['F164']), {'F130', 'F005'}, set(), [], ['also-hydrogenated']),
        ('Pd/H2 and a thioether', ctx(['F164'], ['[H][H]']), {'F073'}, set(), ['poisoned-by-sulfur'], []),
        ('Pd coupling of a thiophene is not hydrogenation', ctx(['F164']), {'F104'}, set(), [], ['poisoned-by-sulfur']),
        ('Raney nickel, H2 and a thioether: desulfurised, not poisoned', ctx(['F164'], ['[H][H]', '[Ni]']), {'F073'}, set(), ['desulfurised-by-Raney-Ni'], ['poisoned-by-sulfur']),
        ('Pd(0) and an allyl carbamate', ctx(['F164'], ['[Pd]']), {'F133'}, set(), ['cleaved-by-Pd0'], []),
        ('Pd(PPh3)4 as typed in the Reagents sheet and an allyl ether', ctx(['F164'], ['c1ccc([P](c2ccccc2)(c2ccccc2)[Pd]([P](c2ccccc2)(c2ccccc2)c2ccccc2)([P](c2ccccc2)(c2ccccc2)c2ccccc2)[P](c2ccccc2)(c2ccccc2)c2ccccc2)cc1']), {'F133'}, set(), ['cleaved-by-Pd0'], []),
        ('platinum and H2 do not ionise an allyl group', ctx(['F164'], ['[Pt]', '[H][H]']), {'F133', 'F005'}, set(), ['also-hydrogenated'], ['cleaved-by-Pd0']),
        ('nickel does not ionise an allyl group', ctx(['F164'], ['[Ni]']), {'F133'}, set(), [], ['cleaved-by-Pd0']),
        ('two leaving groups', ctx(), {'F013', 'F047'}, {'F013'}, ['leaving-group-order'], []),
        ('one leaving group', ctx(), {'F013'}, {'F013'}, [], ['leaving-group-order']),
        ('benzyl bromide alone is one leaving group', ctx(), {'F010', 'F036'}, {'F010'}, [], ['leaving-group-order']),
        ('2-chloropyridine alone is one leaving group', ctx(), {'F047', 'F086'}, {'F047'}, [], ['leaving-group-order']),
        ('one alcohol and a strong base', ctx(parts=['CC(C)(C)[O][K]']), {'F014', 'F003'}, set(), [], ['quench-strong-base-acidic-XH']),
        ('alcohol + acid and a strong base', ctx(parts=['CC(C)(C)[O][K]']), {'F014', 'F001'}, set(), ['quench-strong-base-acidic-XH'], []),
        ('protic solvent, a halide anion and an alkyl bromide', ctx(['F153'], solv=['protic']), {'F010'}, {'F010'}, ['protic-solvent-slows-SN2'], []),
        ('protic solvent without an anionic nucleophile', ctx(solv=['protic']), {'F010'}, {'F010'}, [], ['protic-solvent-slows-SN2']),
        ('no solvent given', ctx(['F153']), {'F010'}, {'F010'}, [], ['protic-solvent-slows-SN2']),
        ('water and a Grignard typed as the reagent', ctx(['F029'], solv=['protic']), set(), set(), ['water-destroys-reagent'], []),
        ('n-butyllithium + aryl bromide', ctx(['F140'], ['[Li][CH2]CCC']), {'F084'}, set(), ['halogen-metal-exchange'], []),
        ('isopropylmagnesium chloride + aryl bromide', ctx(['F141'], ['C[CH](C)[Mg][Cl]']), {'F084'}, set(), ['halogen-metal-exchange'], []),
        ('methylmagnesium bromide adds, it does not exchange', ctx(['F029'], ['[CH3][Mg][Br]']), {'F084', 'F004', 'F018'}, {'F018'}, ['also-attacked-by-organometallic'], ['halogen-metal-exchange']),
        ('methyllithium does not exchange with an aryl bromide', ctx(['F140'], ['[Li][CH3]']), {'F084'}, set(), [], ['halogen-metal-exchange']),
        ('organolithium + aryl chloride', ctx(['F140'], ['[Li][CH2]CCC']), {'F047', 'F170'}, set(), [], ['halogen-metal-exchange']),
        ('acid chloride + amino alcohol, amine is the target', ctx(['F020']), {'F002', 'F003'}, {'F002'}, ['competing-nucleophile-acylation'], []),
        ('nitrite is not an alkylating agent', ctx(parts=['O=N[O-]'], roles=['electrophile / alkylating agent']), {'F002', 'F003'}, {'F002'}, [], ['competing-nucleophile-acylation']),
        ('thionyl chloride + alcohols', ctx(parts=['O=S(Cl)Cl']), {'F001', 'F003'}, {'F001'}, ['activates-every-OH'], []),
        ('PBr3 + alcohols', ctx(parts=['BrP(Br)Br']), {'F001', 'F003'}, {'F001'}, ['activates-every-OH'], []),
        ('fluoride and a silyl ether', ctx(['F154']), {'F052'}, {'F052'}, ['cleaved-by-fluoride'], []),
        ('Raney nickel and a thioether', ctx(parts=['[Ni]']), {'F073'}, set(), ['desulfurised-by-Raney-Ni'], []),
        ('aluminium chloride typed as ions + an amine', ctx(parts=['[Al+3]', '[Cl-]']), {'F002'}, set(), ['consumed-lewis-acid'], []),
        ('BF3 + an ether', ctx(parts=['FB(F)F']), {'F024'}, set(), ['consumed-lewis-acid'], []),
        ('copper iodide is not a Lewis acid here', ctx(parts=['[Cu][I]']), {'F002'}, set(), [], ['consumed-lewis-acid']),
        ('tetrafluoroborate is not a Lewis acid here', ctx(parts=['F[B-](F)(F)F']), {'F002'}, set(), [], ['consumed-lewis-acid']),
        ('PCC and a tertiary alcohol: nothing to oxidise', ctx(parts=['[O]=[Cr](=[O])([O-])[Cl]']), {'F014', 'F016'}, set(), [], ['also-oxidised-by-metal-oxo']),
        ('PCC and a secondary alcohol, another group reacts', ctx(parts=['[O]=[Cr](=[O])([O-])[Cl]']), {'F015', 'F014', 'F005'}, {'F005'}, ['also-oxidised-by-metal-oxo'], []),
        ('Dess-Martin and a primary alcohol', ctx(parts=['CC(=O)OI1(OC(C)=O)(OC(C)=O)OC(=O)c2ccccc21']), {'F003', 'F014', 'F073'}, {'F073'}, ['also-oxidised-by-metal-oxo'], []),
        ('periodate cleaves diols, it does not oxidise alcohols', ctx(['F162'], ['[O-][I+3]([O-])([O-])[O-]']), {'F003', 'F014', 'F073'}, {'F073'}, [], ['also-oxidised-by-metal-oxo']),
        ('manganese dioxide is not a general alcohol oxidant', ctx(['F163'], ['[O]=[Mn]=[O]']), {'F003', 'F014', 'F073'}, {'F073'}, [], ['also-oxidised-by-metal-oxo']),
        ('NaBH4 with NiCl2 reduces nitro groups', ctx(parts=['[BH4-]', '[Na+]', '[Ni+2]', '[Cl-]']), {'F025', 'F019'}, {'F019'}, [], ['spared-by-NaBH4']),
        ('NaBH4 with a Lewis acid is borane-like', ctx(parts=['[BH4-]', '[Na+]', 'FB(F)F'], roles=['Lewis acid']), {'F021', 'F019'}, {'F019'}, [], ['spared-by-NaBH4']),
        ('NaBH4 with iodine is borane', ctx(parts=['[BH4-]', '[Na+]', 'II']), {'F021', 'F019'}, {'F019'}, [], ['spared-by-NaBH4']),
        ('NaBH4 alone leaves a carbamate', ctx(parts=['[BH4-]', '[Na+]']), {'F077', 'F019'}, {'F019'}, ['spared-by-NaBH4'], []),
        ('carbon tetrabromide is a brominating agent, not an alkylating agent', ctx(['F010'], ['BrC(Br)(Br)Br', 'c1ccc(P(c2ccccc2)c2ccccc2)cc1'], ['halogenating agent']), {'F003', 'F002'}, {'F003'}, [], ['competing-nucleophile-acylation']),
        ('iodomethane alkylates amines and alcohols', ctx(['F012'], ['CI'], ['electrophile / alkylating agent']), {'F003', 'F002'}, {'F002'}, ['competing-nucleophile-acylation'], []),
        ('Selectfluor matches the alkyl chloride fragment but is not an alkylating agent', ctx(['F011'], ['F[N+]12CC[N+](CCl)(CC1)CC2'], ['halogenating agent']), {'F003', 'F002'}, {'F002'}, [], ['competing-nucleophile-acylation']),
        ('iodine is not a bromination or radical halogenation reagent', ctx(['F152'], ['II'], ['halogenating agent']), {'F005', 'F028'}, {'F005'}, [], ['competing-radical-or-ring-bromination']),
        ('NBS and an alkene', ctx(['F151'], ['O=C1CCC(=O)N1Br'], ['halogenating agent']), {'F005', 'F028'}, {'F028'}, ['competing-radical-or-ring-bromination'], []),
        ('sodium hydroxide is not strong enough to rank the acidic sites', ctx(['F060'], ['[OH-]']), {'F001', 'F014', 'F028'}, {'F001'}, [], ['quench-strong-base-acidic-XH']),
        ('sodium hydride and two acidic groups', ctx(['F143'], ['[H-]']), {'F014', 'F028'}, set(), ['quench-strong-base-acidic-XH'], []),
        ('a mis-drawn orthoformate is not an alkoxide base', ctx(['F060'], ['CCOC([O-])[O-]']), {'F001', 'F014', 'F004'}, {'F004'}, [], ['cleaved-by-base-nucleophile', 'quench-strong-base-acidic-XH']),
        ('PPTS does not remove a Boc group', ctx(parts=['c1cc[nH+]cc1']), {'F078', 'F039'}, set(), [], ['cleaved-by-acid']),
        ('nitrite and an aniline + a phenol (diazotisation with a coupling partner)', ctx(parts=['O=N[O-]']), {'F026', 'F028'}, {'F026'}, ['competing-diazonium-coupling'], []),
    ]


def main():
    cases = all_cases()
    lib = I.Library()
    bad = 0
    for smi, yes, no in cases:
        got = set(lib.flags_of_smiles(smi))
        miss = [k for k in yes if k not in got]
        extra = [k for k in no if k in got]
        for k in yes + no:
            if k not in D.SITE_BY_KEY: print('UNKNOWN KEY in a test:', k); bad += 1
        if miss or extra:
            bad += 1
            print(f'FAIL {smi}: missing {miss}; wrongly present {extra}; got {sorted(got)}')
    covered = {k for _, y, _n in cases for k in y}
    uncovered = [k for k in lib.keys if k not in covered]
    na = 0
    for name, ctx, present, consumed, yes, no in advisory_cases():
        got = {a['key'] for a in I.advisories(ctx, present, consumed)}
        miss, extra = [k for k in yes if k not in got], [k for k in no if k in got]
        na += 1
        if miss or extra:
            bad += 1
            print(f'FAIL advisory "{name}": missing {miss}; wrongly shown {extra}; got {sorted(got)}')
    print(f'{len(cases)} site cases and {na} advisory cases, {bad} failures; {len(covered & set(lib.keys))} of {len(lib.keys)} site patterns have a positive example')
    if uncovered: print('no positive example for:', ', '.join(uncovered))
    if '--json' in sys.argv:
        out = sys.argv[sys.argv.index('--json') + 1]
        adv = [{'name': name, 'ctx': {k: sorted(v) for k, v in ctx.items()}, 'present': sorted(present), 'consumed': sorted(consumed), 'yes': yes, 'no': no,
                'got': sorted(a['key'] for a in I.advisories(ctx, present, consumed))} for name, ctx, present, consumed, yes, no in advisory_cases()]
        json.dump({'flags': [{'smiles': s, 'yes': y, 'no': n, 'flags': lib.flags_of_smiles(s)} for s, y, n in cases], 'advisories': adv}, open(out, 'w'))
        print('wrote', out)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
