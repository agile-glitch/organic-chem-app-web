"""Sub-fragments: finer descriptions of groups that tools/subfragment_scan.py showed react differently in the recorded reactions.

Each entry splits a PARENT fragment of the Fragments sheet by something its match can be told apart by (the atom it sits on, what is next to it, what the
rest of the ring carries). The SMARTS is written out in full (not as the parent + a constraint) so it reads as a fragment of its own; pos / neg are molecules
it must / must not match. tools/subfragment_eval.py measures each against the recorded reactions: how many starting materials carry it, how often the parent
group is consumed with and without it for the same kind of reagent, and how much that tells on reactions the rates never saw.
"""

EWG = '[$([N+](=[OX1])[O-]),$([CX2]#[NX1]),$([CX4](F)(F)F),$([CX3]=[OX1]),$([SX4](=[OX1])=[OX1])]'          # nitro, cyano, CF3, carbonyl, sulfonyl
EDG = '[$([OX2H1]),$([OX2][#6]),$([NX3;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-[#7,#8])])]'                          # OH, OR, amino (not acylated / sulfonylated)


def ortho_para(group):
    """a ring atom with `group` on a neighbouring (ortho) or opposite (para) ring atom"""
    return f'$(c:c-{group}),$(c:c:c:c-{group})'


SUBFRAGMENTS = [
    # ---- halides: what the ring carries decides whether they react ----
    dict(slug='azine_azole_chloride', parent=['F047'], name='aryl chloride next to a ring N (2-chloropyridine, 2-chlorothiazole, 2-chloropyrimidine)',
         smarts='[Cl;X1]-c:[nX2;+0]', pos=['Clc1ccccn1', 'Clc1nccs1', 'Clc1ncccn1'], neg=['Clc1ccccc1', 'Clc1ccc(Cl)cc1']),
    dict(slug='azine_azole_bromide', parent=['F013'], name='aryl bromide next to a ring N (2-bromopyridine, 2-bromothiazole)',
         smarts='[Br;X1]-c:[nX2;+0]', pos=['Brc1ccccn1', 'Brc1nccs1'], neg=['Brc1ccccc1']),
    dict(slug='ewg_aryl_chloride', parent=['F047'], name='aryl chloride with an ortho/para electron-withdrawing group (SNAr-activated)',
         smarts=f'[Cl;X1]-[c;{ortho_para(EWG)}]', pos=['Clc1ccc([N+](=O)[O-])cc1', 'Clc1ccccc1C#N', 'Clc1ccc(C(C)=O)cc1'], neg=['Clc1ccccc1', 'Clc1ccc(OC)cc1']),
    dict(slug='edg_aryl_chloride', parent=['F047'], name='aryl chloride with an ortho/para electron-donating group (deactivated)',
         smarts=f'[Cl;X1]-[c;{ortho_para(EDG)}]', pos=['Clc1ccc(OC)cc1', 'Clc1ccc(N(C)C)cc1', 'Clc1ccccc1O'], neg=['Clc1ccccc1', 'Clc1ccc([N+](=O)[O-])cc1']),
    dict(slug='benzylic_chloride', parent=['F011'], name='benzylic chloride (Ar–CH2–Cl)', smarts='c-[CX4;H2]-[Cl;X1]', pos=['ClCc1ccccc1', 'ClCc1ccccn1'], neg=['CCCl']),
    dict(slug='azole_halide_5ring', parent=['F047', 'F013', 'F084'], name='2-halo-azole (Cl/Br/I on a 5-ring carbon next to a ring N: thiazole, imidazole, oxazole)',
         smarts='[Cl,Br,I;X1]-[c;r5]:[nX2;r5;+0]', pos=['Clc1nccs1', 'Brc1ncc[nH]1', 'Brc1nc2ccccc2s1'], neg=['Clc1ccccn1', 'Clc1ccsc1']),

    # ---- arene C–H: ring electronics decide electrophilic substitution ----
    # ---- amines: nucleophilicity and basicity depend on what the N sits on ----
    dict(slug='alkylamine_1', parent=['F002'], name='primary alkylamine (R–NH2 on sp3 C; not an aniline)',
         smarts='[NX3;H2;+0;!$(N~[!#6;!#1]);!$(N-[#6]=[!#6]);!$(N-[#6]#[!#6])]-[CX4]', pos=['CCN', 'NCc1ccccc1', 'NC1CCCCC1'], neg=['Nc1ccccc1', 'NC(C)=O']),
    dict(slug='cyclic_amine_2', parent=['F034'], name='cyclic secondary amine (piperidine, pyrrolidine, morpholine, piperazine N–H)',
         smarts='[NX3;H1;R;!$(NC=O);!$(N-c)]([CX4])[CX4]', pos=['C1CCNCC1', 'C1COCCN1', 'C1CCNC1'], neg=['CCNCC', 'CNc1ccccc1']),
    dict(slug='acyclic_amine_2', parent=['F034'], name='acyclic secondary amine (R2NH: Et2NH, N-methylbenzylamine)',
         smarts='[NX3;H1;!R;!$(NC=O);!$(N-c)]([CX4])[CX4]', pos=['CCNCC', 'CNCc1ccccc1'], neg=['C1CCNCC1']),
    dict(slug='trialkylamine_acyclic', parent=['F093'], name='trialkylamine, all three substituents acyclic (Et3N, DIPEA, Me2NR): a base',
         smarts='[NX3;H0;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N);!$(N-[#7,#8,#16])]([CX4;!R])([CX4;!R])[CX4;!R]', pos=['CCN(CC)CC', 'CC(C)N(CC)C(C)C', 'CN(C)CCO'], neg=['CN1CCCCC1', 'CN(C)c1ccccc1']),
    dict(slug='dialkylaniline', parent=['F093'], name='N,N-dialkylaniline (Ar–NR2: dimethylaniline, N-arylpiperazine/piperidine)',
         smarts='[NX3;H0;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-C#N);!$(N-[#7,#8,#16])](-c)([CX4])[CX4]', pos=['CN(C)c1ccccc1', 'C1CCN(CC1)c1ccccc1', 'c1ccccc1N1CCNCC1'], neg=['CCN(CC)CC', 'CN1CCCCC1']),
    dict(slug='azine_secondary_amine', parent=['F094'], name='secondary amine on a ring N-heteroarene (2-alkylamino-pyridine / pyrimidine / thiazole: Ar(het)–NH–R)',
         smarts='[NX3;H1;+0;!$(N-[C,S,P]=[O,S,N]);!$(N-[#7,#8,#16])](-c:[nX2;+0])-[#6]', pos=['CNc1ccccn1', 'CCNc1ncccn1', 'CNc1nccs1'], neg=['CNc1ccccc1']),
    dict(slug='aminoazole_nh2', parent=['F115'], name='2-aminoazole NH2 (2-aminothiazole, 2-aminoimidazole, 2-aminooxazole: NH2 on a 5-ring carbon next to a ring N)',
         smarts='[nX2;r5]:c-[NX3;H2;+0]', pos=['Nc1nccs1', 'Nc1ncc[nH]1', 'Nc1nccn1C'], neg=['Nc1ccccn1', 'Nc1ncccn1']),
    dict(slug='aminoazine_nh2', parent=['F115'], name='aminopyridine / aminopyrimidine NH2 (NH2 on a 6-ring carbon next to a ring N)',
         smarts='[nX2;r6]:c-[NX3;H2;+0]', pos=['Nc1ccccn1', 'Nc1ncccn1'], neg=['Nc1nccs1', 'Nc1ccccc1']),

    # ---- amides, imides, carbamates ----
    dict(slug='tertiary_amide', parent=['F021'], name='tertiary amide (R–C(=O)–NR2: DMF, DMA, Weinreb amide, N-acylpiperidine, N-alkyl lactam)',
         smarts='[CX3;!$(C(-[OX2,SX2,NX3])-[OX2,SX2,NX3])](=[OX1])[NX3;H0;+0]', pos=['CN(C)C=O', 'CON(C)C(C)=O', 'O=C(C)N1CCCCC1', 'CN1CCCC1=O'], neg=['CNC(C)=O', 'NC(C)=O', 'CN(C)C(=O)OC', 'CN(C)C(=O)N(C)C']),
    dict(slug='secondary_amide', parent=['F021'], name='secondary amide (R–C(=O)–NH–R)',
         smarts='[#6][CX3](=[OX1])[NX3;H1;+0]-[#6]', pos=['CNC(C)=O', 'O=C(NC)c1ccccc1', 'O=C1CCCN1'], neg=['CN(C)C=O', 'NC(C)=O']),
    dict(slug='ring5_nh_amide', parent=['F095'], name='5-membered ring N–H amide / carbamate / urea (γ-lactam, oxazolidinone, hydantoin, thiazolidinedione N–H)',
         smarts='[NX3;H1;+0;r5;!a]-[CX3;!a]=[OX1]', pos=['O=C1NCCO1', 'O=C1CSC(=O)N1', 'O=C1CNC(=O)N1', 'O=C1CCCN1'], neg=['O=C1CCCCN1', 'CNC(C)=O']),
    dict(slug='imide_nh', parent=['F121'], name='imide N–H (succinimide, glutarimide, thiazolidinedione, maleimide N–H): an acid',
         smarts='[CX3](=[OX1])[NX3;+0;H1][CX3]=[OX1]', pos=['O=C1CCC(=O)N1', 'O=C1CSC(=O)N1', 'O=C1C=CC(=O)N1'], neg=['O=C1c2ccccc2C(=O)N1C', 'O=C1CCC(=O)N1C']),
    dict(slug='imide_nsub', parent=['F121'], name='N-substituted imide (N-alkyl phthalimide, N-alkyl succinimide: masked amine)',
         smarts='[CX3](=[OX1])[NX3;+0;H0]([#6,#7,#8])[CX3]=[OX1]', pos=['O=C1c2ccccc2C(=O)N1CC', 'O=C1CCC(=O)N1C', 'O=C1c2ccccc2C(=O)N1O'], neg=['O=C1CCC(=O)N1', 'O=C1CCC(=O)N1Br']),

    # ---- esters, acids, carbonyls ----
    dict(slug='enone_acyclic', parent=['F040'], name='acyclic enone (chalcone, methyl vinyl ketone: a strong Michael acceptor)',
         smarts='[CX3;!R]=[CX3][CX3](=[OX1])[#6]', pos=['C=CC(C)=O', 'O=C(C=Cc1ccccc1)c1ccccc1', 'CC(=O)C=CC'], neg=['O=C1C=CCCC1', 'CCC(C)=O']),
    dict(slug='enone_cyclic', parent=['F040'], name='cyclic enone (cyclohexenone, cyclopentenone, chromone: ring C=C–C=O)',
         smarts='[CX3;R]=[CX3;R][CX3;R](=[OX1])[#6]', pos=['O=C1C=CCCC1', 'O=C1C=CCC1'], neg=['C=CC(C)=O']),
    dict(slug='aryl_nitrile', parent=['F022'], name='aryl nitrile (Ar–C≡N: benzonitrile, cyanopyridine)', smarts='c-[CX2]#[NX1]', pos=['N#Cc1ccccc1', 'N#Cc1ccccn1'], neg=['CC#N', 'N#CCc1ccccc1']),
    dict(slug='ketone_alpha_ch', parent=['F101'], name='α-C–H of a ketone (R–CO–CH: pKa ≈ 19–20; aldol, enolate alkylation, α-halogenation)',
         smarts='[CX4;!H0;!$(C-[F,Cl,Br,I]);!$(C(-[CX3]=[OX1])-[CX3]=[OX1]);!$(C(-[CX3]=[OX1])-C#N);$(C-[CX3](=[OX1])-[#6])]', pos=['CC(C)=O', 'O=C1CCCCC1', 'CCC(=O)c1ccccc1'], neg=['CC(=O)OC', 'CC#N']),
    dict(slug='ester_alpha_ch', parent=['F101'], name='α-C–H of an ester or lactone (pKa ≈ 25; needs LDA/NaH-strength base)',
         smarts='[CX4;!H0;!$(C-[F,Cl,Br,I]);!$(C(-[CX3]=[OX1])-[CX3]=[OX1]);!$(C(-[CX3]=[OX1])-C#N);$(C-[CX3](=[OX1])-[OX2;H0;+0][#6])]', pos=['CC(=O)OC', 'CCC(=O)OCC', 'O=C1CCCO1'], neg=['CC(C)=O']),
    dict(slug='nitrile_alpha_ch', parent=['F101'], name='α-C–H of a nitrile (pKa ≈ 25; cyanomethyl anion)',
         smarts='[CX4;!H0;!$(C-[F,Cl,Br,I]);!$(C(-[CX3]=[OX1])-C#N);!$(C(-C#N)-C#N);$(C-[CX2]#[NX1])]', pos=['CC#N', 'N#CCc1ccccc1'], neg=['CC(C)=O']),
    dict(slug='aldehyde_alpha_ch', parent=['F101'], name='α-C–H of an aldehyde (pKa ≈ 17; enolises most easily)',
         smarts='[CX4;!H0;!$(C-[F,Cl,Br,I]);$(C-[CX3H1]=[OX1])]', pos=['CC=O', 'O=CC(C)C'], neg=['CC(C)=O']),
    dict(slug='imine_aldimine', parent=['F031'], name='aldimine / aldoxime (H–C=N: one substituent on carbon)', smarts='[CX3;H1]=[NX2]', pos=['CC=NC', 'CC=NO', 'c1ccccc1C=Nc1ccccc1'], neg=['CC(C)=NC']),
    # ---- alcohols, ethers, acetals: small ones behave as reagents ----
    dict(slug='methanol', parent=['F014'], name='methanol (CH3–OH): reagent and solvent', smarts='[CH3]-[OX2H1]', pos=['CO'], neg=['CCO']),
    dict(slug='ethanol', parent=['F003'], name='ethanol (CH3CH2–OH): reagent and solvent', smarts='[CH3]-[CH2]-[OX2H1]', pos=['CCO'], neg=['CO', 'CCCO']),
    dict(slug='benzyl_alcohol', parent=['F130'], name='benzylic alcohol (Ar–CH2–OH)', smarts='c-[CX4;H2]-[OX2H1]', pos=['OCc1ccccc1', 'OCc1ccncc1'], neg=['COCc1ccccc1']),
    dict(slug='orthoester', parent=['F039'], name='orthoester (R–C(OR)3: triethyl orthoformate, trimethyl orthoacetate)', smarts='[CX4]([OX2][#6])([OX2][#6])[OX2][#6]', pos=['CCOC(OCC)OCC', 'COC(C)(OC)OC'], neg=['COC(C)OC', 'C1OCCO1']),
    dict(slug='amide_acetal', parent=['F039'], name='amide acetal / aminal ether (DMF-DMA: Me2N–CH(OMe)2)', smarts='[CX4]([OX2][#6])([OX2][#6])[NX3]', pos=['COC(OC)N(C)C'], neg=['CCOC(OCC)OCC']),

    dict(slug='lactam', parent=['F021'], name='lactam, 5- to 7-membered (pyrrolidinone, piperidinone, caprolactam, oxindole)',
         smarts='[#6][CX3;r5,r6,r7](=[OX1])[NX3;r5,r6,r7;+0]', pos=['O=C1CCCN1', 'O=C1CCCCN1', 'CN1CCCC1=O', 'O=C1CCCCCN1'], neg=['CNC(C)=O', 'CN(C)C=O', 'O=C1CCN1']),
    dict(slug='enamine_acyclic', parent=['F127'], name='acyclic enamine / enaminone (N–CH=C: Me2N–CH=CH–C(=O)R from DMF-DMA)',
         smarts='[CX3]=[CX3;!R]-[NX3;+0;!$(N-C=O);!$(N-S(=O)=O);!$(N-c)]', pos=['CN(C)C=CC(C)=O', 'CN(C)C=C'], neg=['C1=CCCN1C', 'CC(=O)NC=C']),
    dict(slug='imine_acyclic', parent=['F031'], name='acyclic imine / oxime / hydrazone (C=N not in a ring)', smarts='[CX3]=[NX2;!R]', pos=['CC=NC', 'CC(C)=NO', 'CC=NNC'], neg=['C1=NCCC1']),

    # ---- sulfur ----
    dict(slug='azine_sulfone', parent=['F071'], name='2-sulfonyl / 2-sulfinyl azine (SO2R on a ring carbon next to a ring N: an SNAr leaving group)',
         smarts='[#6]-[S;$([SX3]=[OX1]),$([SX4](=[OX1])=[OX1])]-c:[nX2;+0]', pos=['CS(=O)(=O)c1ncccn1', 'CS(=O)c1ncccn1', 'CS(=O)(=O)c1nccs1'], neg=['CS(=O)(=O)c1ccccc1']),
]
