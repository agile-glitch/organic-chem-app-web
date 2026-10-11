"""The descriptions of the candidate fragments in tools/new_fragments.py, in the columns of the Fragments sheet of molecule_data.xlsx:
reactive (Reactive atom(s)), role (Role), acid (pKa of its most acidic H), base (pKa of its conjugate acid), lg (Leaving group),
steric (Steric hindrance at the reactive atom), notes (Notes).

The pKa cells follow the sheet's convention (tools/rules_to_js.py first_num, tools/pka.py): the FIRST number in the cell is the typical
value the page uses (acid <= 20 makes the fragment an acidic site, base >= 0 a basic site on an N), so a cell with a group that is
not an acid / base starts with words and no number, or with a number outside those ranges. Typical values in water unless stated
(Evans pKa table; Bordwell DMSO scale where marked; Clayden, Organic Chemistry): approximate, the rest of the molecule shifts them.
Evidence (how often the group is found, how often it reacts) is added by tools/add_fragments.py from tools/fragment_candidates_scored.json.
"""

NONE = 'None'
NOACID = 'no acidic H'

TEXT = {
 # ---------------- ring C-H sites
 'thiophene_ch': dict(
    reactive='the ring C–H (C2/C5 first, then C3/C4)',
    role='nucleophile: electrophilic aromatic substitution, much faster than on benzene and at C2/C5 (NBS or Br2 bromination, Vilsmeier–Haack formylation POCl3/DMF, Friedel–Crafts acylation with SnCl4 or TiCl4, nitration with Ac2O/HNO3). C2–H is lithiated by n-BuLi (THF, −78 °C → 0 °C) and trapped with DMF, I2, R3SiCl or B(OR)3; direct Pd-catalysed C–H arylation (Pd(OAc)2, PivOH, K2CO3)',
    acid='ring C–H ≈ 33 (C2–H, DMSO estimate); not a useful acid without a strong base', base=NONE, lg=NONE,
    steric='C2/C5 are open; a 2-substituent sends electrophiles to C5, and C3 reacts only when both α positions are blocked or by directed lithiation',
    notes='Covers thiophene, benzothiophene and thienopyridines. The ring C–H a generic arene C–H (F009) cannot tell apart from a benzene C–H'),
 'furan_ch': dict(
    reactive='the ring C–H (C2/C5 first)',
    role='nucleophile: electrophilic substitution (Vilsmeier–Haack formylation, mild halogenation NBS/DMF, Friedel–Crafts acylation with mild Lewis acids BF3·OEt2) at C2/C5; lithiated at C2 by n-BuLi. The ring is a diene (Diels–Alder with maleimides) and opens in strong aqueous acid',
    acid='ring C–H ≈ 35 (C2–H, estimate)', base=NONE, lg=NONE,
    steric='C2/C5 are open; avoid strong aqueous acid, which opens the ring to a 1,4-dicarbonyl',
    notes='Covers furan and benzofuran. More acid-sensitive than thiophene'),
 'pyrrole_ch': dict(
    reactive='the ring C–H (C2 of pyrrole, C3 of indole)',
    role='nucleophile: the most reactive common arene toward electrophiles: Vilsmeier–Haack formylation (POCl3/DMF) , Mannich (CH2O + R2NH, AcOH), halogenation (NBS or NIS, 0 °C, DMF), Friedel–Crafts acylation with Et2AlCl or acid chlorides, azo coupling, Michael addition to enones (indole C3, Lewis acid). Strong acid polymerises pyrrole; an N-sulfonyl or N-Boc group tames it',
    acid='ring C–H ≈ 39 (estimate)', base='protonates on carbon: ≈ −3.8 (pyrrole) / −3.6 (indole)', lg=NONE,
    steric='indole C3 is open; a 3-substituent sends electrophiles to C2, 2,3-disubstituted indoles to C6 or N',
    notes='Covers N–H and N-substituted pyrroles and indoles. Reactions at the ring N–H are in F089 and the pyrrole N–H fragments'),
 'azole_c2h': dict(
    reactive='the ring C–H between the two heteroatoms (C2)',
    role='the most acidic ring C–H of an azole: lithiated at C2 by n-BuLi or LDA (THF, −78 °C) and trapped (DMF, I2, R3SiCl, ZnCl2 then Negishi); direct Pd/Cu C–H arylation (Pd(OAc)2, CuI, Cs2CO3, PivOH); N-heterocyclic carbene precursor after N-alkylation at N3 (imidazolium C2–H ≈ 21 in DMSO); electrophilic substitution is slow here (it goes to C4/C5)',
    acid='≈ 30 (C2–H of 1-methylimidazole and thiazole, DMSO estimates); imidazolium C2–H ≈ 21 (DMSO)', base=NONE, lg=NONE,
    steric='open unless C4/C5 carry large groups', notes='Covers imidazole, benzimidazole, thiazole, benzothiazole, oxazole and benzoxazole C2–H, and 1,2,4-triazole C–H'),
 'azine_ch': dict(
    reactive='the ring C–H (C2/C4/C6 attacked by nucleophiles and radicals, C3/C5 by electrophiles)',
    role='an electron-poor ring: electrophilic substitution needs forcing conditions and goes to C3/C5 (nitration, halogenation at 150–300 °C), while radicals attack C2/C4 (Minisci: R• from RCO2H/AgNO3/(NH4)2S2O8 or alkyl halides, in acid) and organolithiums add at C2 (Chichibabin: NaNH2 gives the 2-amine). Deprotonated ortho to a directing group with LDA, LiTMP or TMPMgCl·LiCl; the N-oxide opens C2 to POCl3 (2-chloropyridine) and Pd C–H arylation',
    acid='ring C–H ≈ 40 (estimate)', base=NONE, lg=NONE,
    steric='positions next to the ring N or a substituent are crowded', notes='Covers pyridine, pyrimidine, pyrazine, pyridazine, quinoline and isoquinoline ring C–H (the ring that carries the N). Benzo-ring C–H of quinoline is an ordinary arene C–H (F009)'),
 # ---------------- azole N-H and N
 'pyrazole_nh': dict(
    reactive='the ring N–H (after deprotonation, N1 or N2)',
    role='N–H acid, then N-nucleophile: N-alkylation (R–X or R–OMs with K2CO3/Cs2CO3 in DMF or MeCN, or NaH in THF/DMF), Mitsunobu, Michael addition to acrylates, N-arylation (Chan–Lam Cu(OAc)2/ArB(OH)2/air; Ullmann CuI/diamine/K3PO4; SNAr on activated ArF), N-protection (SEM-Cl, THP by DHP/TsOH, Boc2O, trityl, PMB-Cl). Unsymmetrical pyrazoles and indazoles give N1/N2 mixtures',
    acid='pyrazole 14.2 (19.8 in DMSO); indazole 13.9; 1,2,3-triazole 9.3; 4-nitropyrazole 9.6', base='pyrazole 2.5; indazole 1.3 (the ring N lone pair)', lg='none here (the azolide is a good leaving group in N-acyl and N-sulfonyl azoles)',
    steric='alkylation goes mostly to the N next to the smaller neighbour (the less hindered N); a 3-substituent on pyrazole favours N1 alkylation, indazole N1 (thermodynamic, NaH/THF) vs N2 (kinetic, Mitsunobu, acidic alkylation)',
    notes='N–H next to another ring N. A narrower split of F089'),
 'imidazole_nh': dict(
    reactive='the ring N–H (after deprotonation, either N)',
    role='N–H acid, then N-nucleophile: N-alkylation (R–X, K2CO3 or Cs2CO3 in DMF/MeCN or NaH), N-arylation (Ullmann CuI/ligand, Buchwald Pd, SNAr), N-acylation and N-sulfonylation (gives the activated azolides used as acyl transfer agents: CDI, N-acylbenzotriazoles), N-protection (trityl, SEM, Boc). Unsymmetrical imidazoles and benzimidazoles give N1/N3 mixtures',
    acid='imidazole 14.5; benzimidazole 12.8; purine 8.9; 7-azaindole ≈ 15', base='imidazole 7.0; benzimidazole 5.5; 7-azaindole 4.6 (the pyridine-type N)', lg='none here (the azolide is a good leaving group in N-acyl azoles)',
    steric='alkylation of an unsymmetrical imidazole goes mostly to the N far from the bulkier substituent (and the more basic N under neutral conditions)',
    notes='N–H next to a C=N: imidazole, benzimidazole, purine, azaindole, 1,2,4-triazole. A narrower split of F089'),
 'pyrrole_nh': dict(
    reactive='the ring N–H (after deprotonation)',
    role='weaker N–H acid than the azoles, so it needs a stronger base: NaH or KOtBu or KOH/DMSO or Cs2CO3 (heated) for N-alkylation, N-sulfonylation (ArSO2Cl/NaH or Bu4NHSO4/NaOH), N-acylation, N-Boc (Boc2O/DMAP), N-arylation (Ullmann, Buchwald, Chan–Lam). Under neutral or acidic conditions the C-nucleophilic positions react instead (indole C3)',
    acid='pyrrole 17.5 (23 in DMSO); indole 16.2 (21 in DMSO); carbazole ≈ 20 (DMSO)', base='protonates on carbon: ≈ −3.8 (pyrrole) / −3.6 (indole), not a base', lg=NONE,
    steric='groups at C2 (pyrrole) or C2/C7 (indole) slow N-alkylation and N-arylation', notes='N–H with no second ring N next to it: pyrrole, indole, carbazole. A narrower split of F089'),
 'imidazole_n': dict(
    reactive='the ring N lone pair (sp2, in the ring plane)',
    role='base and nucleophile: N-alkylation to imidazolium salts (R–X in MeCN, heated), N-acylation (the intermediate of CDI and DMAP-type catalysis), metal ligand (Zn, Cu, Fe, Pd), and the proton acceptor in general base catalysis; the less hindered N of an unsymmetrical imidazole is alkylated',
    acid=NOACID, base='imidazole 7.0; benzimidazole 5.5; purine N 2.4; 1-methylimidazole 7.1', lg=NONE,
    steric='a substituent at C2 or C4/C5 next to the N slows alkylation', notes='sp2 N of imidazole, benzimidazole and purine (next to a C that carries the substituted ring N). A narrower view of what F091 does for 6-ring azines'),
 'azole_n': dict(
    reactive='the ring N lone pair (sp2, in the ring plane)',
    role='weak base and nucleophile, metal ligand: N-alkylation (R–X, Me3O⁺ BF4⁻), N-oxidation is rare, coordination to Pd/Cu directs C–H activation next to it; protonation at this N activates C2 of thiazole/oxazole to nucleophiles; the ring N of a 2-substituted thiazole or oxazole is not alkylated easily',
    acid=NOACID, base='pyrazole 2.5; thiazole 2.5; 1,2,4-triazole 2.2; 1,2,3-triazole 1.2; oxazole 0.8; tetrazole ≈ −3', lg=NONE,
    steric='substituents at the carbons next to the N slow alkylation and metal binding', notes='sp2 N of pyrazole, thiazole, oxazole, triazole and tetrazole rings'),
 'tetrazole': dict(
    reactive='the ring N–H (acid) and the ring N atoms (alkylation)',
    role='carboxylic-acid isostere: the N–H is as acidic as a carboxylic acid. Made from a nitrile and NaN3 (NH4Cl or ZnBr2 or Et3N·HCl, DMF or water, 100–130 °C) by [3+2]; alkylation gives N1/N2 mixtures (N2 favoured for bulky electrophiles), protected as N-trityl, N-PMB or N-cyanoethyl; stable to most reagents (not to strong reducing conditions)',
    acid='tetrazole N–H 4.9 (5-aryltetrazoles ≈ 4.5)', base='≈ −3 (very weak base)', lg=NONE,
    steric='alkylation at N1 is hindered by a 5-substituent, so N2 products dominate for 5-aryl tetrazoles', notes='The aromatic 1,2,3,4-tetrazole ring'),
 'hetaryl_nh2': dict(
    reactive='the ring N (protonation, alkylation) and the NH2 N',
    role='a much weaker nucleophile than an aniline or alkylamine, because the NH2 lone pair is delocalised into the amidine-like N–C=N: amide coupling needs acid chlorides, HATU or T3P with heat or a stronger base (NaH, LiHMDS, pyridine/DMAP); Buchwald–Hartwig N-arylation (Pd2dba3, Xantphos or BINAP, Cs2CO3); diazotisation (NaNO2/HCl or t-BuONO/CuBr) to the halide; alkyl halides react on the ring N first',
    acid='N–H ≈ 28 (2-aminopyridine, DMSO)', base='about 6.9 for 2-aminopyridine; 5.4 for 2-aminothiazole; 3.5 for 2-aminopyrimidine; 3.1 for 2-aminopyrazine (the ring N is the basic site)', lg=NONE,
    steric='a substituent next to the NH2 (position 3 of a 2-aminopyridine) slows acylation and arylation', notes='NH2 on a carbon next to a ring N (2-aminopyridine, aminopyrimidine, aminothiazole, aminopyrazine). The aniline NH2 (F026) reacts quite differently'),
 # ---------------- carbonyl and related groups
 'aryl_ketone': dict(
    reactive='the carbonyl C (electrophile) and the α C–H',
    role='electrophile less reactive than a dialkyl ketone (conjugated): NaBH4 reduction (MeOH, 0 °C), Grignard/RLi addition, reductive amination (Ti(OiPr)4 then NaBH4), Wolff–Kishner (N2H4, KOH) or Et3SiH/TFA to the CH2, Baeyer–Villiger (the aryl group migrates less than a tertiary alkyl), α-bromination (CuBr2, Br2/AcOH, NBS/TsOH), Claisen–Schmidt aldol with ArCHO (NaOH/EtOH), Willgerodt; the enolate (LDA, LiHMDS) alkylates',
    acid='α C–H 18.3 (acetophenone; 24.7 in DMSO)', base='≈ −6 (O-protonated ketone)', lg=NONE,
    steric='an ortho substituent twists the carbonyl out of conjugation and slows additions', notes='Ar–C(=O)–C: acetophenones, aryl ketones, benzophenones (diaryl). A narrower split of F019'),
 'aryl_aldehyde': dict(
    reactive='the carbonyl C (electrophile)',
    role='electrophile with no α C–H (cannot enolise): reductive amination (NaBH(OAc)3 or NaBH3CN, AcOH, DCE/MeOH), Wittig and Horner–Wadsworth–Emmons olefination, Knoevenagel (malononitrile, piperidine), Claisen–Schmidt aldol with ketones, Henry (MeNO2), imine, oxime and hydrazone formation, reduction (NaBH4, 0 °C), oxidation to the acid (NaClO2/NaH2PO4 Pinnick; KMnO4), Cannizzaro in conc. NaOH, acetal formation (HC(OMe)3/TsOH)',
    acid=NOACID + ' (no α C–H)', base='≈ −7 (O-protonated aldehyde)', lg=NONE,
    steric='ortho substituents slow additions; electron-poor aldehydes (4-nitro, pyridine-carbaldehyde) are more reactive and hydrate', notes='Ar–CHO: benzaldehydes and heteroaryl aldehydes. A narrower split of F018'),
 'aryl_acid': dict(
    reactive='the carboxylic acid C (after activation) and the O–H',
    role='acylating agent after activation: amide coupling (HATU/DIPEA, EDC/HOBt, T3P, SOCl2 or (COCl)2/DMF to the acid chloride), esterification (MeOH/H2SO4, MeI/K2CO3, TMSCHN2), Curtius (DPPA, Et3N then t-BuOH or amine), reduction (BH3·THF, LiAlH4), Weinreb amide; salts (ArCO2⁻) direct ortho-lithiation (s-BuLi/TMEDA); decarboxylation needs ortho EWG or heat/Cu (quinoline, 200 °C)',
    acid='benzoic acid 4.2; 4-nitro 3.4; 4-methoxy 4.5; nicotinic acid 4.8', base='≈ −7 (O-protonated acid)', lg='OH is made a leaving group by activation (EDC, HATU, SOCl2); the carboxylate itself is not',
    steric='ortho substituents slow coupling and esterification', notes='Ar–CO2H and Ar–CO2⁻. A narrower split of F001'),
 'amino_acid': dict(
    reactive='the amine N (nucleophile), the carboxylic acid C (after activation) and the α C–H (racemisation)',
    role='peptide building block: N-protection (Boc2O/NaOH or Et3N, Fmoc-OSu or Fmoc-Cl/NaHCO3, CbzCl/Na2CO3), amide coupling of the CO2H (EDC/HOBt, HATU/DIPEA, PyBOP, T3P), esterification (SOCl2/MeOH, TMSCl/MeOH), reduction to the amino alcohol (LiAlH4, BH3·SMe2). The α C–H racemises through the oxazolone or enolate with excess base or over-activation (use HOAt/Oxyma, collidine, 0 °C)',
    acid='α-carboxylic acid ≈ 3.6 (N-acyl and N-carbamate amino acids; the free zwitterion 2.3)', base=NONE, lg=NONE,
    steric='α,α-disubstituted and N-methyl amino acids couple slowly (use HATU/HOAt, PyBrOP, acid fluorides)', notes='N–CH–CO2H, free or N-protected (Boc, Cbz, Fmoc, acetyl); includes proline. The amine of a free amino acid is also F002/F034 (its pKaH ≈ 9.6)'),
 'primary_amide': dict(
    reactive='the carbonyl C, the N–H and the N',
    role='dehydration to the nitrile (POCl3, SOCl2, TFAA/pyridine, Burgess reagent, cyanuric chloride), Hofmann rearrangement to the amine or carbamate (NaOBr, PIDA or PIFA, Br2/NaOMe), hydrolysis to the acid (6 M HCl heat, NaOH/H2O2), reduction to the amine (LiAlH4, BH3), N-arylation (Goldberg CuI/DMEDA, Buchwald), conversion to the N-acylsulfonamide',
    acid='N–H ≈ 17 (25.5 in DMSO)', base='≈ −0.5 (O-protonated amide)', lg=NONE,
    steric='open', notes='C(=O)NH2 on carbon. A narrower split of F021 and F095'),
 'imide': dict(
    reactive='the carbonyl C (ring opening), the N–H (acid) and the N',
    role='N–H acid: the potassium salt of phthalimide alkylates with R–Br (K2CO3 or KOH, DMF) in the Gabriel synthesis, and the N-alkylphthalimide releases the primary amine with hydrazine hydrate (EtOH, reflux) or MeNH2 (Ing–Manske). Also Mitsunobu on the N–H, ring opening by NaOH/MeOH, NaBH4 reduction of one C=O to the hydroxylactam, N-bromosuccinimide and N-hydroxysuccinimide chemistry',
    acid='phthalimide 8.3; succinimide 9.6; glutarimide 11.4', base='≈ −2 (O-protonated imide)', lg='the imide N⁻ (pKa 8–10) leaves from acyl imides',
    steric='open at the N–H; N-substituents slow the ring opening', notes='C(=O)–N–C(=O): phthalimide, succinimide, glutarimide, maleimide, hydantoin-type imides. A narrower split of F021 and F095'),
 'urea': dict(
    reactive='the N–H and the carbonyl C',
    role='formed from isocyanates and amines, or CDI/triphosgene/phosgene-equivalent and two amines; very stable to hydrolysis; N-alkylation of the N–H (NaH, R–X), cyclic ureas (imidazolidin-2-ones) from 1,2-diamines and CDI, intramolecular cyclisations to hydantoins and dihydrouracils; ureas are good H-bond donors (organocatalysis)',
    acid='N–H ≈ 27 (urea, DMSO)', base='protonates on O, not on N', lg=NONE,
    steric='N,N-disubstituted ureas are twisted and unreactive', notes='N–C(=O)–N: ureas, cyclic ureas, carbamoyl ureas. Overlaps F021 and F095'),
 'hydroxamic': dict(
    reactive='the N–OH (acid) and the carbonyl C',
    role='metal chelator (HDAC and metalloprotease inhibitors). Made from an ester and NH2OH·HCl/KOH (MeOH), or an acid chloride/activated acid and O-protected hydroxylamine (THPONH2, BnONH2, TBSONH2) then deprotected (H2/Pd–C, TFA, TBAF). Lossen rearrangement of an O-activated hydroxamate gives the isocyanate; the Weinreb amide (N-methoxy-N-methyl) adds one RMgX/RLi to give a ketone, or DIBAL/LiAlH4 to give an aldehyde, without over-addition (chelated tetrahedral intermediate)',
    acid='acetohydroxamic acid 9.4; benzohydroxamic acid 8.8 (the O–H)', base=NONE, lg='the N–O⁻ after O-activation (Lossen); MeO(Me)N–H leaves from the Weinreb tetrahedral intermediate on work-up',
    steric='N,N-disubstituted (Weinreb) amides are open at C=O', notes='C(=O)–N–OH and C(=O)–N–OR, including Weinreb amides and O-benzyl hydroxamates'),
 'active_ester': dict(
    reactive='the carbonyl C',
    role='selective acylating agent: reacts with amines (and hydrazines) at rt in DMF, DMSO, CH2Cl2 or aqueous buffer (pH 7–9) to give amides, much more selectively than acid chlorides and without racemisation; stable enough to isolate (NHS, pentafluorophenyl, p-nitrophenyl esters). HOBt/HOAt esters are the intermediates inside EDC/HATU couplings. Made from the acid with NHS/EDC or DCC, HOBt/EDC, PFP-OH/EDC, or 4-nitrophenyl chloroformate',
    acid=NOACID, base='≈ −7 (O-protonated ester)', lg='N-hydroxysuccinimide anion (HOSu pKa 6.0); HOBt 4.6; HOAt 3.3; pentafluorophenoxide 5.5; 4-nitrophenoxide 7.1',
    steric='hindered amines (secondary, α,α-disubstituted) react slowly: heat or use HOAt/HATU', notes='C(=O)–O–N or C(=O)–O–Ar(EWG): N-hydroxysuccinimide, HOBt/HOAt, pentafluorophenyl and p-nitrophenyl esters'),
 'alpha_halocarbonyl': dict(
    reactive='the halogen-bearing α C (SN2) and the carbonyl C',
    role='a very reactive SN2 electrophile (the C=O stabilises the transition state): amines (K2CO3 or DIPEA, MeCN or DMF, 0 °C → rt), thiolates, azide, carboxylates (Cs2CO3/DMF) and phenoxides (K2CO3, acetone) displace the halide; Hantzsch thiazole synthesis with thioamides and thioureas (EtOH, reflux), imidazo[1,2-a]pyridines with 2-aminopyridines, Darzens (with aldehydes, base), Favorskii (α-halo ketone + alkoxide), Reformatsky (Zn); lachrymatory and toxic',
    acid='the α C–H next to the halogen ≈ 16–17 (estimate)', base=NONE, lg='I⁻ > Br⁻ > Cl⁻ (HI −10, HBr −9, HCl −7)',
    steric='tertiary α-halo carbonyls (α-bromoisobutyrates) are hindered for SN2 but ionise and eliminate; secondary react more slowly', notes='X–C–C(=O): α-halo ketones (phenacyl bromide), esters (bromoacetate), amides (chloroacetamide)'),
 'enol_ether': dict(
    reactive='the β C (nucleophile) and the acetal-forming α C',
    role='an electron-rich alkene and acid-labile: aqueous acid (HCl, TFA, PPTS/MeOH) hydrolyses it to the carbonyl and the alcohol; dihydropyran (DHP) protects alcohols as THP ethers (PPTS or TsOH cat.); also reacts with electrophiles at the β carbon (NBS/H2O to the α-bromo ketone, Vilsmeier, Mukaiyama aldol of silyl enol ethers), Diels–Alder (inverse demand), [2+2], Heck arylation with Pd, ozonolysis. Not covering enol esters (vinyl acetate) or silyl enol ethers',
    acid=NOACID, base='protonates on carbon: ≈ −2 (alkyl vinyl ether, estimate)', lg='the alcohol R–OH leaves on hydrolysis',
    steric='open', notes='C=C–O–R (not C=C–O–C=O): vinyl ethers, dihydropyran, dihydrofuran, glycals'),
 'enamine': dict(
    reactive='the β C (nucleophile) and the N',
    role='nucleophilic at carbon (Stork enamine alkylation, acylation and Michael addition) and hydrolysed by aqueous acid back to the carbonyl and amine; enaminones (β-aminoenones, from ketones with DMF-DMA or Bredereck reagent) are 1,3-bielectrophiles that close pyrazoles (hydrazines), isoxazoles (NH2OH), pyrimidines (amidines, guanidine) and pyridones; Vilsmeier–Haack reaction; reduced by NaBH3CN/AcOH to the amine',
    acid=NOACID, base='iminium after C-protonation ≈ 9–11 (alkyl enamines); enaminones ≈ 3', lg='the amine R2NH leaves in condensations and hydrolysis',
    steric='open; 2,2-disubstituted enamines are slower', notes='C=C–N where N is not an amide or sulfonamide or aniline: enamines, enaminones, DMF-DMA adducts'),
 'diazo': dict(
    reactive='the diazo C',
    role='carbene or 1,3-dipole precursor: methyl esters of carboxylic acids and phenols (CH2N2 or TMSCHN2 in MeOH/benzene), Arndt–Eistert homologation (ArCOCl + CH2N2 → diazoketone → Wolff rearrangement, Ag2O/H2O), Rh2(OAc)4 or Cu-catalysed cyclopropanation and X–H insertion, 1,3-dipolar cycloaddition to alkenes and alkynes (pyrazolines, pyrazoles), Bamford–Stevens; toxic and explosive (diazomethane): TMS-diazomethane is the safer reagent',
    acid=NOACID, base='protonates on carbon: ≈ 10 (CH2N2 → CH3N2⁺, estimate)', lg='N2 leaves (carbene or diazonium)',
    steric='open', notes='R2C=N⁺=N⁻: diazomethane, TMS-diazomethane, ethyl diazoacetate, diazoketones'),
 'amino_alcohol': dict(
    reactive='the amine N (nucleophile) and the O–H',
    role='N,O-bifunctional: the amine acylates, Boc-protects and sulfonylates before the alcohol (Boc2O, Et3N, 0 °C; CbzCl, NaHCO3); with CDI, triphosgene or diethyl carbonate/K2CO3 it closes the oxazolidin-2-one; with aldehydes/ketones the oxazolidine; with SOCl2 or Mitsunobu the aziridine; with DAST the fluoride; the alcohol oxidises (TEMPO/BAIB, Swern, DMP) once the N is protected. Made by reduction of amino acids (LiAlH4, BH3·SMe2, NaBH4/I2) and epoxide aminolysis',
    acid='O–H ≈ 15.5 (alcohol)', base='ethanolamine 9.5 (the ammonium)', lg='OH is made a leaving group by TsCl, MsCl, Mitsunobu',
    steric='the amine is the more nucleophilic site, so selective N-acylation works even with an unprotected alcohol', notes='free N–C–C–OH (1,2-amino alcohol; also with N-alkyl, not an amide or carbamate N)'),
 # ---------------- protecting groups
 'benzyl_ons': dict(
    reactive='the benzylic C–O / C–N / C–S bond (the CH2)',
    role='protecting group removed by hydrogenolysis (H2, 10% Pd/C or Pd(OH)2/C, MeOH/EtOH/EtOAc, 1 atm, rt; with AcOH or HCl for benzylamines; transfer: ammonium formate or cyclohexadiene), BCl3 or BBr3 (CH2Cl2, −78 °C, with pentamethylbenzene), Na/NH3, DDQ (slow for plain Bn), or ACE-Cl for N-benzyl amines. Put on with BnBr/NaH, BnBr/Ag2O (selective), or benzyl trichloroacetimidate/TfOH (base-sensitive)',
    acid=NOACID, base=NONE, lg='the alkoxide, amide or thiolate leaves on hydrogenolysis (as the alcohol, amine or thiol)',
    steric='hindered benzyl ethers hydrogenolyse more slowly; sulfur poisons Pd catalysts', notes='Ar–CH2–O, Ar–CH2–N, Ar–CH2–S (not benzyl esters, Cbz or benzyl halides)'),
 'benzyl_ester': dict(
    reactive='the benzylic C–O bond',
    role='Cbz carbamate and benzyl ester/carbonate: removed by hydrogenolysis (H2, Pd/C, rt, neutral, gives toluene + CO2 + the amine or acid), HBr/AcOH (33%), TMSI, BCl3, or Na/NH3; stable to TFA and to mild bases; the benzyl ester also saponifies (LiOH). Cbz is put on with CbzCl and Na2CO3 or NaHCO3 (Schotten–Baumann) or Cbz-OSu',
    acid=NOACID, base='≈ −6.5 (O-protonated ester)', lg='the carboxylate or carbamate anion (leaves as CO2 + amine)',
    steric='open', notes='Ar–CH2–O–C(=O): benzyl esters, benzyl carbamates (Cbz, Z) and benzyl carbonates'),
 'pmb': dict(
    reactive='the benzylic C–O / C–N bond and the electron-rich ring',
    role='protecting group removed oxidatively (DDQ in CH2Cl2/H2O, 0 °C → rt, or CAN in MeCN/H2O), by TFA (with anisole or triethylsilane as cation scavenger) or by hydrogenolysis; stable to bases; put on with PMB-Cl/NaH or PMB trichloroacetimidate/TfOH. Orthogonal to benzyl (DDQ takes PMB and leaves Bn)',
    acid=NOACID, base=NONE, lg='the alkoxide or amide leaves from the PMB cation',
    steric='open', notes='4-MeO–C6H4–CH2–O, –N or –S'),
 'allyl_ons': dict(
    reactive='the allylic C–O / C–N bond and the C=C',
    role='allyl ethers, allyl esters and allyl carbamates (Alloc): removed by Pd(0) (Pd(PPh3)4 with a scavenger: morpholine, dimedone, PhSiH3, Bu3SnH, N,N′-dimethylbarbituric acid) in THF/CH2Cl2 at rt, or isomerised (KOtBu/DMSO, Rh or Ir) to the enol ether and hydrolysed with mild acid; Tsuji–Trost allylic substitution; aryl allyl ethers undergo the Claisen rearrangement (heat, 180–220 °C); the C=C is cleaved by ozonolysis or cross-metathesis',
    acid=NOACID, base=NONE, lg='the alkoxide, carboxylate or carbamate leaves as the π-allyl palladium complex forms',
    steric='open', notes='CH2=CH–CH2–O, –N, –S'),
 'trityl': dict(
    reactive='the trityl C–O / C–N bond',
    role='bulky protecting group for primary OH, NH, SH and the N of imidazoles, pyrazoles and tetrazoles; removed by mild acid (1–5% TFA/CH2Cl2 with Et3SiH, 80% AcOH, HCl/MeOH, formic acid) through the stable, orange trityl cation; put on with TrCl, pyridine or Et3N/DMAP; stable to bases and nucleophiles',
    acid=NOACID, base=NONE, lg='the alcohol, amine or azole leaves from the trityl cation (pKR+ of Ph3C⁺ ≈ −6.6)',
    steric='very bulky: selective for primary over secondary alcohols', notes='Ph3C–N, Ph3C–O, Ph3C–S'),
 'fmoc': dict(
    reactive='the fluorenyl C9–H (acidic) and the carbamate C',
    role='base-labile amine protecting group: 20% piperidine in DMF (rt, 5–20 min), DBU, morpholine or Et2NH eliminate through the fluorenyl anion to dibenzofulvene (trapped by piperidine) + CO2 + the free amine; stable to TFA and, mostly, to H2/Pd; put on with Fmoc-Cl (NaHCO3, dioxane/H2O) or Fmoc-OSu; the basis of solid-phase peptide synthesis',
    acid='fluorenyl C–H ≈ 22.6 (fluorene, DMSO)', base=NONE, lg='dibenzofulvene + CO2 + the amine',
    steric='open', notes='Fmoc carbamate: fluorenyl-CH2–O–C(=O)–N'),
 'tfa_acyl': dict(
    reactive='the acid O–H (trifluoroacetic acid) or the carbonyl C (amide, ester)',
    role='trifluoroacetic acid (TFA): strong organic acid used neat or 20–50% in CH2Cl2 (rt, 0.5–2 h) to remove Boc, tert-butyl esters and ethers and trityl (with scavengers: Et3SiH, anisole, thioanisole, water); salts of amines are trifluoroacetates. The trifluoroacetamide protects amines (TFAA/Et3N or ethyl trifluoroacetate), is removed by K2CO3/MeOH–H2O or NaBH4, and the trifluoroacetate ester is labile to methanol',
    acid='trifluoroacetic acid 0.5 (the O–H)', base='≈ −7 (O-protonated acid)', lg='CF3CO2⁻ (pKa of the acid 0.5) is a good leaving group; the trifluoroacetamide N⁻ is a poorer one',
    steric='open', notes='CF3–C(=O)–O/N: trifluoroacetic acid, trifluoroacetamides, trifluoroacetate esters'),
 'acetamide': dict(
    reactive='the carbonyl C and the N–H',
    role='N-acetyl protecting group or the product of acetylation (Ac2O or AcCl with Et3N/pyridine): robust; removed by 6 M HCl or aqueous NaOH/EtOH at reflux, hydrazine, Boc2O/DMAP then base, or Meerwein salt then hydrolysis; on an aniline it is a moderate ortho/para director (acetanilide) and a weaker activator than NH2',
    acid='N–H ≈ 17 (secondary amide; acetanilide 21.5 in DMSO)', base='≈ −0.5 (O-protonated amide)', lg=NONE,
    steric='N-aryl and N,N-disubstituted acetamides hydrolyse more slowly', notes='CH3–C(=O)–N: acetamides and acetanilides. A narrower split of F021'),
 'acetate_ester': dict(
    reactive='the carbonyl C',
    role='acetyl protecting group on an alcohol or phenol (Ac2O, pyridine, DMAP): removed by K2CO3/MeOH (rt), NaOMe/MeOH (Zemplén), NH3/MeOH, LiOH/THF–H2O or lipases; more labile than benzoate and pivalate (so it can be removed selectively), less than formate; acetyl groups migrate between neighbouring hydroxyls under base; phenyl acetates are cleaved even by mild base or amines',
    acid='α C–H ≈ 25 (ester enolate, DMSO)', base='≈ −6.5 (O-protonated ester)', lg='the alkoxide or phenoxide (pKa of the alcohol 15–16, phenol 10)',
    steric='primary acetates cleave faster than secondary and tertiary', notes='CH3–C(=O)–O–C: alkyl and aryl acetates. A narrower split of F004'),
 'methyl_ethyl_ester': dict(
    reactive='the carbonyl C and the α C–H',
    role='the common ester protecting group and synthetic intermediate: saponified by LiOH (THF/MeOH/H2O, rt), NaOH, KOH or K2CO3/MeOH–H2O to the acid; aminolysis with AlMe3 (Weinreb) or by heating with amines; reduction to the alcohol (LiAlH4, LiBH4, NaBH4/CaCl2) or to the aldehyde (DIBAL at −78 °C); Grignard double addition to a tertiary alcohol; Claisen and Dieckmann condensations (NaOEt, NaH, KOtBu); α-alkylation (LDA, −78 °C); cleaved with LiI/pyridine, BBr3, Me3SnOH or TMSOK in the presence of other esters. Hydroxide takes Me/Et esters in the presence of tert-butyl and benzyl esters',
    acid='α C–H ≈ 25 (ester enolate, DMSO)', base='≈ −6.5 (O-protonated ester)', lg='the alkoxide (MeO⁻, EtO⁻: pKa of the alcohol ≈ 15.5)',
    steric='methyl esters hydrolyse about 2–3 times faster than ethyl esters; α-branching and ortho substituents slow it', notes='C–C(=O)–O–CH3 and C–C(=O)–O–CH2CH3. A narrower split of F004'),
 # ---------------- reagents
 'organolithium': dict(
    reactive='the C–Li carbon',
    role='very strong base and carbon nucleophile (n-BuLi in hexanes/THF, s-BuLi, t-BuLi, MeLi, PhLi): halogen–lithium exchange of ArBr/ArI (−78 °C, THF or Et2O; 2 equiv t-BuLi), directed ortho-metalation (s-BuLi/TMEDA; directing groups OMe, CONR2, CO2⁻, oxazoline, F, Cl), deprotonation of terminal alkynes and weak C–H acids, addition to aldehydes, ketones, nitriles, Weinreb amides, CO2, DMF and B(OR)3; n-BuLi slowly cleaves THF above about −20 °C; t-BuLi is pyrophoric',
    acid='conjugate acid, the C–H ≈ 50 (butane)', base=NONE, lg=NONE,
    steric='t-BuLi and s-BuLi are bulkier and more basic; n-BuLi adds to hindered carbonyls less well than it deprotonates', notes='[Li]–C. Destroyed by water, alcohols, N–H and acidic C–H'),
 'grignard_ci': dict(
    reactive='the C–Mg carbon',
    role='carbon nucleophile and base: RMgCl and RMgI (MeMgCl, EtMgBr is F029, iPrMgCl, PhMgCl, vinylMgBr): addition to aldehydes, ketones (tertiary alcohols), esters (double addition), Weinreb amides (ketones), nitriles (ketones after hydrolysis), epoxides, CO2 (acids); iPrMgCl·LiCl (Knochel turbo-Grignard) does halogen–magnesium exchange on ArBr/ArI at −20 to 25 °C and tolerates esters and nitriles; Kumada coupling (Ni or Pd, Fe); Cu(I) catalysis gives 1,4-addition',
    acid='conjugate acid, the C–H ≈ 45–50', base=NONE, lg=NONE,
    steric='bulky Grignards (tBuMgCl) reduce rather than add to hindered ketones', notes='C–Mg–Cl and C–Mg–I (the bromides are F029). Destroyed by water, alcohols and N–H'),
 'amide_base': dict(
    reactive='the anionic N (a base, not a nucleophile)',
    role='strong, hindered, non-nucleophilic bases: LDA, LiHMDS, NaHMDS, KHMDS and LiTMP (THF, −78 °C): kinetic enolates of ketones, esters and amides, α-alkylation, aldol, Claisen and Dieckmann cyclisations, nitrile and sulfone α-deprotonation, ortho-lithiation with LiTMP; LiHMDS with Pd catalysts aminates ArX as an ammonia equivalent; KHMDS/NaHMDS give the more reactive, more thermodynamic (E/Z) enolates',
    acid='conjugate acid, diisopropylamine ≈ 36 (THF); 2,2,6,6-tetramethylpiperidine ≈ 37; hexamethyldisilazane ≈ 26', base=NONE, lg=NONE,
    steric='very hindered, so it deprotonates rather than adds', notes='[N⁻] with two silyl or two alkyl groups (lithium, sodium or potassium amides)'),
 'hydride': dict(
    reactive='the H⁻ (a base; rarely a hydride donor)',
    role='metal hydride as base: NaH (60% in mineral oil; DMF, THF, 0 °C → rt, H2 evolves) deprotonates alcohols, phenols, amides, azole N–H, thiols and 1,3-dicarbonyl compounds before alkylation or acylation (Williamson, N-alkylation, Claisen); KH is more soluble and basic; CaH2 dries solvents; reacts violently with water',
    acid='conjugate acid, dihydrogen ≈ 35', base=NONE, lg=NONE,
    steric='a heterogeneous base; reaction starts at the solid surface', notes='H⁻: sodium hydride, potassium hydride, calcium hydride'),
 'borohydride': dict(
    reactive='the B–H hydride',
    role='mild hydride donors: NaBH4 reduces aldehydes and ketones (MeOH or EtOH, 0 °C; THF/H2O), acid chlorides and imines/iminiums, but not esters, amides, nitriles or acids unless activated (CaCl2, LiCl, I2, hot MeOH/THF, CeCl3 for Luche 1,2-reduction of enones); NaBH3CN reductively aminates at pH 4–6 (MeOH, AcOH), selective for the iminium over the carbonyl; NaBH(OAc)3 (STAB) reductive amination in DCE/THF with AcOH, tolerating acetals, esters, nitro; LiBH4 reduces esters (THF)',
    acid=NOACID, base=NONE, lg=NONE,
    steric='NaBH(OAc)3 is bulky and selective; NaBH4 reduces unhindered carbonyls fastest', notes='BH4⁻, BH3CN⁻, BH(OAc)3⁻ (sodium, lithium, potassium salts and tetraalkylammonium borohydrides)'),
 'alkyl_alh': dict(
    reactive='the Al–H hydride',
    role='DIBAL-H (1 M in toluene, hexane or CH2Cl2): ester → aldehyde (1 equiv, −78 °C, stops at the hemiacetal), nitrile → aldehyde (after hydrolysis), lactone → lactol, Weinreb amide → aldehyde; with excess or warmer: alcohols; α,β-unsaturated esters → allylic alcohols (2+ equiv); it is a Lewis acid and also opens epoxides and acetals. Red-Al (sodium bis(2-methoxyethoxy)aluminium hydride) reduces like LiAlH4 and dissolves in toluene. Work-up with Rochelle salt or the Fieser method',
    acid=NOACID, base=NONE, lg=NONE,
    steric='bulky and Lewis-acidic: coordinates the carbonyl O first', notes='R2Al–H: DIBAL-H and relatives (LiAlH4 is F042)'),
 'phosphine': dict(
    reactive='the P lone pair',
    role='nucleophile and reductant: triphenylphosphine in the Mitsunobu reaction (with DEAD or DIAD: alcohol → ester, azide, ether, imide with inversion), Appel (CBr4 or CCl4: alcohol → halide), Staudinger (azide + PPh3 then H2O → amine), Wittig salt formation (+ R–X → phosphonium), deoxygenation of epoxides, sulfoxides and N-oxides; as a ligand on Pd/Rh/Ru/Ir (PPh3, PCy3, P(tBu)3, XPhos, SPhos, dppf, BINAP, Xantphos) for Suzuki, Heck, Sonogashira, Buchwald–Hartwig; oxidised to Ph3P=O (remove by crystallisation, ZnCl2 or silica)',
    acid=NOACID, base='conjugate acids: triphenylphosphine 2.7; tributylphosphine 8.4; tricyclohexylphosphine 9.7; tri-tert-butylphosphine 11.4', lg=NONE,
    steric='cone angle sets the ligand character: PPh3 145°, PCy3 170°, P(tBu)3 182°', notes='trialkyl and triaryl phosphines (P with three carbons), not phosphites or phosphine oxides (F103)'),
 'phosphonium': dict(
    reactive='the ylide C (nucleophile) and the P⁺',
    role='Wittig reagents: a phosphonium salt R3P⁺–CH2R′ X⁻ with a base (n-BuLi, NaHMDS, KOtBu, NaH; K2CO3 for stabilised salts) gives the ylide, which condenses with an aldehyde or ketone to the alkene + Ph3P=O (non-stabilised ylides: mostly Z; ester- or ketone-stabilised: E; Schlosser modification: E); also the Horner–Wadsworth–Emmons relatives (F103); phosphonium salts are phase-transfer catalysts and Mitsunobu intermediates',
    acid='methyltriphenylphosphonium ≈ 22 (DMSO); the ester-stabilised ylide precursor ≈ 8.5 (DMSO)', base=NONE, lg='Ph3P=O leaves (the driving force)',
    steric='hindered ketones react poorly with stabilised ylides', notes='R4P⁺ and R3P=CR2: phosphonium salts and phosphorus ylides'),
 'carbodiimide': dict(
    reactive='the central C',
    role='coupling and dehydrating agent: EDC·HCl (water-soluble urea, with HOBt/HOAt/Oxyma, DMAP cat.), DCC (insoluble DCU), DIC (soluble urea; solid-phase) activate a carboxylic acid to the O-acylisourea, which gives amides with amines, esters with alcohols (Steglich, DMAP), anhydrides, or rearranges to the unreactive N-acylurea (suppressed by HOBt); also Pfitzner–Moffatt oxidation (DCC/DMSO/H3PO4) and dehydration of primary amides and aldoximes to nitriles; sensitiser (DCC)',
    acid=NOACID, base=NONE, lg='the urea (DCU, EDU) is the by-product',
    steric='DCC and DIC are hindered nucleophiles but good dehydrating agents', notes='N=C=N: EDC, DCC, DIC and the other carbodiimides'),
 'uronium_coupling': dict(
    reactive='the uronium C / P⁺ (acyl transfer after the carboxylate attacks it)',
    role='amide coupling reagents: HATU (HOAt-based, fastest, least racemisation), HBTU, TBTU, HCTU (uronium/guanidinium), BOP, PyBOP, PyAOP (phosphonium), with DIPEA, NMM or Et3N in DMF, NMP or CH2Cl2 (0 °C → rt, 1–16 h); the carboxylate attacks the reagent to give the OBt/OAt active ester, which acylates the amine; excess reagent guanidinylates the amine (use ≤ 1 equiv or pre-activate the acid 1–5 min); BOP gives carcinogenic HMPA (use PyBOP); by-product tetramethylurea',
    acid=NOACID, base=NONE, lg='HOAt (pKa 3.3) or HOBt (4.6) leaves from the active ester; tetramethylurea or tris(pyrrolidino)phosphine oxide is the by-product',
    steric='works for hindered acids and amines better than carbodiimides; use HATU/HOAt for N-methyl and α,α-disubstituted amino acids', notes='C(=N⁺)(N)O–N / P⁺(N)(N)(N)O–N: HATU, HBTU, TBTU, HCTU, PyBOP, BOP'),
 'hobt': dict(
    reactive='the N–O⁻ / O–H oxygen',
    role='additive that makes carbodiimide couplings clean: traps the O-acylisourea as the HOBt/HOAt active ester, which acylates amines with little racemisation and without N-acylurea; HOAt is the more reactive (neighbouring-group effect of the pyridine N). Anhydrous HOBt is explosive (sold as the hydrate); OxymaPure and 6-Cl-HOBt are alternatives',
    acid='HOBt 4.6; HOAt 3.3 (the O–H)', base=NONE, lg='HOBt/HOAt leave from their active esters',
    steric='open', notes='N-hydroxybenzotriazole and N-hydroxy-7-azabenzotriazole'),
 'n_halo_imide': dict(
    reactive='the N–X halogen (X⁺ source)',
    role='electrophilic halogen: NBS (MeCN, DMF or CH2Cl2, 0 °C → rt) brominates activated arenes (anisoles, phenols, anilines, indoles, pyrroles, thiophenes, pyrazole C4), ketone α-positions (acid or TsOH cat.) and, with a radical initiator (AIBN, BPO or light; CCl4, PhCF3 or MeCN, reflux), allylic and benzylic C–H (Wohl–Ziegler); bromohydrins and bromoethers from alkenes (NBS/H2O–DMSO or ROH); NCS chlorinates, NIS iodinates (arenes with TFA/TfOH or alkynes with Ag⁺); use 1.0–1.1 equiv to avoid polyhalogenation; succinimide by-product',
    acid=NOACID, base=NONE, lg='the imide anion (succinimide N–H pKa 9.6) leaves; Br⁺ is transferred',
    steric='open; benzylic bromination selects the most stabilised radical position', notes='N(C=O)–Cl, –Br, –I: N-bromosuccinimide, N-chlorosuccinimide, N-iodosuccinimide, DBDMH, trichloroisocyanuric acid, chloramine-T'),
 'dihalogen': dict(
    reactive='the X–X bond',
    role='Br2 (CH2Cl2 or AcOH, 0 °C → rt): anti addition to alkenes and alkynes (via the bromonium ion), α-bromination of ketones (acid-catalysed enol; AcOH/HBr cat.), aromatic bromination (FeBr3 or Fe for benzene, none needed for phenols and anilines), Hell–Volhard–Zelinsky of carboxylic acids (PBr3 cat., Br2); Cl2 similar (radical chlorination of alkanes with light); I2 (with HIO3, H5IO6, CAN or Ag⁺ for arenes; base for phenols and imidazoles) for iodination, iodolactonisation (I2/NaHCO3), the haloform and iodoform reactions of methyl ketones (X2/NaOH); ICl and IBr are more electrophilic',
    acid=NOACID, base=NONE, lg='X⁻ (HBr −9, HCl −7, HI −10)',
    steric='open', notes='Br–Br, Cl–Cl, I–I, I–Cl, I–Br'),
 'halide_anion': dict(
    reactive='the anion (nucleophile, counter-ion)',
    role='NaI in acetone (Finkelstein: R–Cl or R–Br → R–I, driven by the precipitation of NaCl/NaBr); LiCl or LiBr in DMF (SN2 on tosylates and mesylates; Krapcho decarboxylation with LiCl, DMSO, H2O, 150 °C; dealkylation of methyl esters with LiI); halide opens epoxides and aziridines; Cu-catalysed aromatic Finkelstein (CuI/diamine, NaI); with an oxidant (Oxone, NaOCl, H2O2) Br⁻ and I⁻ become Br⁺ and I⁺ equivalents; Bu4NBr/Bu4NI as phase-transfer catalysts; as the counter-ion of hydrochlorides and of ammonium and phosphonium salts it is a spectator',
    acid='conjugate acids HI −10, HBr −9, HCl −7', base=NONE, lg=NONE,
    steric='iodide is the most nucleophilic in protic solvents; chloride the most in DMF', notes='Cl⁻, Br⁻, I⁻'),
 'fluoride_anion': dict(
    reactive='the anion',
    role='F⁻ sources (TBAF in THF; CsF; KF/18-crown-6; Et3N·3HF; HF·pyridine): removes silyl groups (TBS, TIPS, TBDPS, TMS; TBAF/THF, 0 °C → rt; buffer with AcOH for base-sensitive substrates) and C–SiMe3 groups (alkyne and aryl desilylation), SNAr fluorination of activated aryl chlorides and nitroarenes (halex: CsF or KF, DMSO or sulfolane, 120–180 °C), SN2 on alkyl sulfonates (TBAF, tert-amyl alcohol), activates boronic acids and silanes in Suzuki and Hiyama couplings; naked F⁻ in aprotic solvents is strongly basic (E2 elimination)',
    acid='HF 3.2', base=NONE, lg=NONE,
    steric='TBAF is hygroscopic and basic; anhydrous F⁻ is a poor nucleophile in protic solvents (strong H-bonding)', notes='F⁻'),
 'azide_anion': dict(
    reactive='the terminal N (nucleophile)',
    role='NaN3 (DMF, DMSO or water/acetone, 50–80 °C; TMSN3): SN2 on alkyl halides and sulfonates (with inversion) to alkyl azides, which are reduced to amines (H2/Pd, PPh3/H2O Staudinger, Zn/NH4Cl) or clicked with alkynes (CuAAC, RuAAC); epoxide opening (NaN3, NH4Cl, MeOH/H2O); SNAr on activated haloarenes and halo-azines; acyl azides (Curtius) from acid chlorides; tetrazoles from nitriles (NaN3, NH4Cl or ZnBr2, 100–130 °C); toxic; HN3 is explosive (do not acidify; avoid CH2Cl2 and CHCl3, which form diazidomethane)',
    acid='hydrazoic acid 4.7', base=NONE, lg=NONE,
    steric='a small, strong nucleophile: SN2 even at hindered secondary centres', notes='N3⁻ (sodium, potassium azide; azide from TMSN3 and fluoride)'),
 'carbonate_anion': dict(
    reactive='the O⁻ (base), seldom the carbonyl C',
    role='the common mild inorganic base: K2CO3, Cs2CO3, Na2CO3, NaHCO3 in DMF, MeCN, acetone, THF/H2O or dioxane/H2O: O-, N- and S-alkylation of phenols, azoles, amides and thiols (Cs2CO3 is more soluble and gives a more naked anion), Suzuki, Buchwald–Hartwig, Heck and Ullmann couplings, SNAr, ester hydrolysis in MeOH/H2O (K2CO3), and neutralisation or work-up (saturated NaHCO3); saturated NaHCO3 buffers at pH ≈ 8.3',
    acid='bicarbonate → carbonate 10.3 (water); carbonic acid 6.4', base=NONE, lg=NONE,
    steric='K2CO3 is poorly soluble in organic solvents (solid-liquid); Cs2CO3 and K3PO4 dissolve better', notes='CO3²⁻ and HCO3⁻. Recorded reactions often list carbonate as a reactant because the atom mapping assigned one of its O atoms to the product: such templates are artefacts'),
 'carbon_dioxide': dict(
    reactive='the carbon (electrophile)',
    role='electrophile for carbanions: ArLi, RLi and Grignard reagents add to CO2 (dry ice or gas; THF or Et2O, −78 °C → rt) to give, after acid, the carboxylic acid; Kolbe–Schmitt (phenoxide + CO2, 125 °C, 5 atm) gives salicylic acids; amines add reversibly to give carbamate salts (R2N–CO2⁻), which Cs2CO3/R–X/Bu4NI alkylate to carbamates; it is also what leaves in decarboxylations (Boc removal, malonic acid, Krapcho) and the C1 source of cyclic carbonates (epoxide + CO2)',
    acid=NOACID, base=NONE, lg=NONE,
    steric='open', notes='O=C=O'),
 'formaldehyde': dict(
    reactive='the carbonyl C (electrophile)',
    role='the most reactive aldehyde, with no α C–H: formalin (37% aqueous with methanol), paraformaldehyde (depolymerises with heat or base), trioxane. Hydroxymethylation of enolates and active methylenes (aldol, NaOH or K2CO3), Mannich (R2NH·HCl + CH2O + ketone), Eschweiler–Clarke N-methylation (CH2O/HCO2H), reductive amination to N-methyl (CH2O + NaBH3CN or NaBH(OAc)3), Pictet–Spengler (tryptamines, acid), Prins (alkenes, acid), Blanc chloromethylation (CH2O/HCl/ZnCl2), Cannizzaro (conc. NaOH), MOM and methylene acetals',
    acid=NOACID, base='≈ −7 (O-protonated aldehyde)', lg=NONE,
    steric='open', notes='CH2=O: formaldehyde in all its forms'),
 'cs2_cscl2': dict(
    reactive='the thiocarbonyl C (electrophile)',
    role='CS2 with an amine and base (Et3N, NaOH) gives the dithiocarbamate (then S-alkylation to R2N–C(=S)–SR, or with EDC/TsCl/Boc2O the isothiocyanate); with an alkoxide the xanthate (Chugaev elimination; Barton–McCombie deoxygenation: ROC(=S)SMe then Bu3SnH/AIBN); with o-aminophenols, hydrazides and anthranilic acids the benzoxazole-, oxadiazole- and quinazoline-thiones; thiophosgene (CSCl2; CHCl3 or CH2Cl2/water, NaHCO3) turns R–NH2 into R–N=C=S and 1,2-diamines into cyclic thioureas; toxic; CS2 flash point −30 °C',
    acid=NOACID, base=NONE, lg='S²⁻/HS⁻ (from CS2) or Cl⁻ (from thiophosgene)',
    steric='open', notes='S=C=S and S=C(Cl)Cl'),
 'borate_ester': dict(
    reactive='the B (electrophile)',
    role='trimethyl and triisopropyl borate trap aryl- and alkyllithiums or Grignards (−78 °C, THF) to give the boronate, hydrolysed with aqueous HCl or NH4Cl to the boronic acid for Suzuki coupling (in situ quench: LDA or LiTMP + ArX in the presence of B(OiPr)3 tolerates halides); pinacol esterification; oxidation to phenols (H2O2/NaOH, NaBO3); B(OMe)3 also traps enolates',
    acid=NOACID, base=NONE, lg='the alkoxide (MeO⁻, iPrO⁻) leaves',
    steric='B(OiPr)3 is bulkier and gives less over-addition (borinic acids) than B(OMe)3', notes='B(OR)3 with three alkyl groups (not boronic acids F081 or pinacol esters F082)'),
 'dast': dict(
    reactive='the S–F fluorine (F⁻ donor) and the S',
    role='deoxofluorination: DAST (Et2NSF3), Deoxo-Fluor ((MeOCH2CH2)2NSF3), XtalFluor-E/M, PyFluor: R–OH → R–F by SN2 with inversion (CH2Cl2, −78 °C → rt), aldehydes and ketones → gem-difluorides (rt to reflux, slower), carboxylic acids → acyl fluorides; elimination and rearrangement compete for β-branched, benzylic and tertiary alcohols; DAST decomposes explosively above ≈ 90 °C (Deoxo-Fluor and XtalFluor are safer); generates HF: glass, quench with cold saturated NaHCO3',
    acid=NOACID, base=NONE, lg='R2N–S(O)F after transfer of F⁻',
    steric='hindered alcohols eliminate instead of fluorinating', notes='R2N–SF3'),
 'hypervalent_iodine': dict(
    reactive='the I(III) / I(V) centre (an oxidant and group-transfer agent)',
    role='Dess–Martin periodinane (DMP: CH2Cl2, rt, optionally buffered with NaHCO3): primary alcohols → aldehydes, secondary → ketones, no over-oxidation, tolerates most functional groups; IBX (DMSO; or EtOAc/MeCN suspension at 80 °C): alcohols → carbonyls, ketones → enones with excess; PIDA/PIFA: Hofmann rearrangement of primary amides, oxidative dearomatisation of phenols (spirocycles), C–H amination and oxidative couplings, TEMPO/PIDA alcohol oxidation; iodonium salts transfer aryl groups',
    acid=NOACID, base=NONE, lg='the iodoarene (ArI) leaves',
    steric='DMP is mild; IBX is insoluble in most solvents', notes='I with three or more neighbours: DMP, IBX, PIDA, PIFA, diaryliodonium salts'),
 'metal_oxo_oxidant': dict(
    reactive='the M=O oxygen (oxidant)',
    role='Jones reagent (CrO3, H2SO4, acetone: secondary alcohols → ketones, primary → acids), PCC/PDC (CH2Cl2 ± Celite or NaOAc: primary → aldehydes), KMnO4 (cold alkaline: alkenes → diols; hot or acidic: cleavage; alkylarenes → ArCO2H; aldehydes → acids), OsO4 (cat. with NMO: alkene → syn-diol; Sharpless asymmetric dihydroxylation with (DHQ)2PHAL/K3Fe(CN)6), RuO4 (RuCl3/NaIO4: alkene cleavage, ether → ester, arene → acid), SeO2 (allylic oxidation; α-oxidation of ketones to 1,2-diketones, dioxane/H2O reflux); toxic: Cr(VI), Os, Se',
    acid=NOACID, base=NONE, lg=NONE,
    steric='OsO4 adds to the less hindered face and the more electron-rich alkene', notes='Cr, Mn, Os, Ru, Se with an oxo or oxide ligand'),
 'pgm_catalyst': dict(
    reactive='the metal (catalyst)',
    role='Pd/C or Pd(OH)2/C with H2 (1 atm–50 psi; EtOH, MeOH, EtOAc): alkene and alkyne → alkane (Lindlar: cis-alkene), nitro, azide, nitrile and C=N → amine, hydrogenolysis of benzyl, Cbz and benzylic C–X, dehalogenation of aryl halides; Pd(PPh3)4, PdCl2(dppf), Pd(OAc)2, Pd2(dba)3 with phosphines: cross couplings (Suzuki with aqueous K2CO3 or K3PO4 in dioxane or toluene/EtOH, 80–110 °C; Buchwald–Hartwig C–N/C–O with XPhos, BINAP or Xantphos and NaOtBu or Cs2CO3; Heck; Sonogashira with CuI; Stille; Negishi; carbonylation); Pt/PtO2 (Adams) and Rh/Al2O3 hydrogenate arenes; Raney Ni desulfurises and reduces; Ru: Noyori transfer hydrogenation and metathesis; Ir: C–H borylation and asymmetric hydrogenation',
    acid=NOACID, base=NONE, lg=NONE,
    steric='ligand bulk sets which halide couples (ArCl needs bulky electron-rich phosphines)', notes='Pd, Pt, Rh, Ru, Ir, Ni in any form (metal on carbon, salts, complexes)'),
 'copper': dict(
    reactive='the metal (catalyst, reagent)',
    role='CuI, CuBr, CuCl, CuCN, Cu(OAc)2, Cu2O: Ullmann and Goldberg C–N, C–O and C–S couplings (CuI + 1,10-phenanthroline, DMEDA or L-proline, K2CO3/K3PO4/Cs2CO3, DMSO or DMF, 90–130 °C), Chan–Lam (Cu(OAc)2, ArB(OH)2 + N–H/O–H, air, pyridine or Et3N, rt), Sandmeyer (ArN2⁺ + CuX/CuCN), Rosenmund–von Braun (ArBr + CuCN, DMF/NMP, 150–200 °C), Sonogashira co-catalyst, click (CuAAC: CuSO4/ascorbate or CuI/DIPEA), organocuprate conjugate additions (R2CuLi; RMgX + cat. CuI), aromatic Finkelstein, trifluoromethylation',
    acid=NOACID, base=NONE, lg=NONE,
    steric='ligands make Cu(I) soluble and lower the temperature', notes='Cu in any oxidation state'),
 'reducing_metal': dict(
    reactive='the metal (electron donor)',
    role='Zn dust (activated with HCl, TMSCl or 1,2-dibromoethane): nitro → amine (Zn/NH4Cl, AcOH), Clemmensen (Zn(Hg)/HCl), Reformatsky, Negishi organozinc insertion, Troc removal; Fe powder with NH4Cl, AcOH or HCl (Béchamp): nitro → aniline, tolerating halogens, esters and nitriles; Sn or SnCl2·2H2O (EtOH reflux): nitro → amine; Mg turnings (activate with I2, 1,2-dibromoethane): Grignard formation in THF or Et2O, Mg/MeOH reductions; Li, Na, K (Birch: Na or Li/NH3 gives 1,4-dienes and trans-alkenes; Bouveault–Blanc; alkoxide formation); Al, Mn, In, Sm (SmI2: ketyl couplings)',
    acid=NOACID, base=NONE, lg=NONE,
    steric='surface area and activation decide the rate', notes='zero-valent Zn, Fe, Sn, Mg, Na, K, Li, Al, Sm, In, Mn, Ca (dust, powder, turnings, wire)'),
 'lewis_acid_halide': dict(
    reactive='the metal or B (electron-pair acceptor)',
    role='BBr3 and BCl3 (CH2Cl2, −78 °C → rt) cleave aryl methyl ethers, benzyl ethers and tert-butyl ethers (1 equiv per Lewis-basic site); AlCl3 drives Friedel–Crafts alkylation and acylation (stoichiometric for acylation), the Fries rearrangement and ArOMe demethylation next to C=O; BF3·OEt2 activates epoxides, acetals and ketals, silyl enol ethers (Mukaiyama aldol) and glycosylations; TiCl4 for Mukaiyama and Evans aldol, TiCl3/Zn or LiAlH4 for McMurry; ZnCl2 in the Fischer indole synthesis and as the Negishi transmetalation partner; SnCl4, FeCl3 for Friedel–Crafts and oxidative couplings; MgBr2·OEt2 chelation control; CeCl3 for Luche reduction',
    acid=NOACID, base=NONE, lg=NONE,
    steric='bulky Lewis acids (MAD, BBr3 at low temperature) are more selective', notes='B, Al, Ti, Zn, Sn, Fe, Mg, Zr, In, Cu, Sc, Yb, Ce, Ga, Mn, Co, Cd, Sb bonded to F, Cl, Br or I and to no carbon'),
 'quat_ammonium': dict(
    reactive='the N⁺ (cation) and the counter-anion',
    role='tetraalkylammonium salts (Bu4N⁺ Br⁻, I⁻, HSO4⁻, F⁻, OH⁻, N3⁻; Aliquat 336; benzyltriethylammonium chloride) are phase-transfer catalysts (1–10 mol%): they carry OH⁻, CN⁻, N3⁻, F⁻, phenoxides and enolates into organic solvents for alkylation, Darzens, epoxidation, hypochlorite oxidations and biphasic NaOH(aq)/toluene or CH2Cl2 reactions; Bu4NI catalyses SN2 by iodide exchange; Hofmann elimination of R4N⁺ OH⁻ on heating gives the least-substituted alkene',
    acid=NOACID, base=NONE, lg='the tertiary amine leaves in Hofmann elimination and SN2 at the α C',
    steric='long chains make the ion pair soluble in organic solvents', notes='R4N⁺ with four carbons on N (not a protonated amine, see ammonium salt)'),
 'ammonium_salt': dict(
    reactive='the N⁺–H (acid)',
    role='amine hydrochlorides and other ammonium salts are the storage form of amines: they need a base (Et3N, DIPEA, NaHCO3, K2CO3) to free the amine before acylation, alkylation or reductive amination; NH4Cl is the mild acid for quenching, with Zn or Fe for nitro reduction, and the ammonia source for primary amides (NH4Cl + acid chloride, ester or anhydride + Et3N) and Strecker reactions; NH4OAc is the imine/enamine source in reductive amination to primary amines (NaBH3CN), Hantzsch and Debus–Radziszewski imidazole syntheses',
    acid='ammonium ion 9.25; triethylammonium 10.75; anilinium 4.6; pyridinium 5.2; piperidinium 11.1', base=NONE, lg=NONE,
    steric='open', notes='R3NH⁺, R2NH2⁺, RNH3⁺, NH4⁺ (a protonated amine or ammonia)'),
}
