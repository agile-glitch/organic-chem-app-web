"""The descriptions of the sub-fragments in tools/new_subfragments.py, in the columns of the Fragments sheet of molecule_data.xlsx:
reactive (Reactive atom(s)), role (Role), acid (pKa of its most acidic H), base (pKa of its conjugate acid), lg (Leaving group),
steric (Steric hindrance at the reactive atom), ring (Ring strain; default 'none: not a ring fragment'), notes (Notes).

The pKa cells follow the sheet's convention (tools/rules_to_js.py first_num): the FIRST number in the cell is the typical value the page uses,
so a group that is not an acid / base starts with words and no number. Typical values in water unless stated (Evans pKa table; Bordwell DMSO
scale where marked; Clayden, Organic Chemistry): approximate, the rest of the molecule shifts them.
The evidence (how many reactions, how often the parent group reacts with and without this feature) is added by tools/add_fragments.py --sub from
tools/subfragment_scored.json.
"""

NONE = 'None'
NOACID = 'no acidic H'

TEXT = {
 # ---------------- halides: what the ring carries decides whether they react
 'azine_azole_chloride': dict(
    reactive='the ring C–Cl carbon (attacked at the ipso position)',
    role='SNAr electrophile: the ring N takes the negative charge of the Meisenheimer complex, so amines (DIPEA, n-BuOH or NMP, 80–150 °C; 2-chloropyrimidines already at rt), alkoxides and thiolates displace Cl⁻ where an ordinary aryl chloride does not react. Suzuki, Buchwald–Hartwig and Negishi couplings work with ordinary Pd catalysts (Pd(PPh3)4, Pd(dppf)Cl2; XPhos or RuPhos for hindered partners). Hydrolysis (aq. HCl or NaOH, heat) to the pyridone-type lactam is the usual side reaction; H2/Pd-C or Zn/AcOH remove the Cl',
    acid=NONE, base='≈ 0.7 (2-chloropyridinium; 2-chloropyrimidine is far less basic)', lg='Cl⁻ leaves after addition; the rate follows how well the ring stabilises the anion (pyrimidine, N-oxide and nitro-pyridine > pyridine)',
    steric='a substituent next to the C–Cl slows attack; a 3-substituent has little effect on a coupling',
    notes='Cl–c:n. 2-Chloropyridine, 2-chloropyrimidine, 2-chloropyrazine, 2-chlorothiazole, 2-chlorobenzothiazole, 1-chloroisoquinoline. A narrower split of F047 (aryl chloride), which it separates from the ordinary chloroarene that needs a Pd catalyst'),
 'azine_azole_bromide': dict(
    reactive='the ring C–Br carbon',
    role='the easiest ring bromide to couple: Suzuki–Miyaura, Stille, Negishi, Sonogashira and Buchwald–Hartwig reactions (Pd(PPh3)4 or Pd(dppf)Cl2, 80–100 °C) and halogen–metal exchange (n-BuLi, −78 °C, or i-PrMgCl·LiCl) at the C–Br; SNAr with amines and alkoxides is slower than for the chloride or fluoride at the same position but works in hot DMSO or NMP',
    acid=NONE, base='≈ 0.9 (2-bromopyridinium)', lg='Br⁻ leaves from the Meisenheimer complex or by oxidative addition to Pd(0)',
    steric='a substituent next to the C–Br slows coupling; the ring N can bind Pd, so bidentate ligands (dppf, Xantphos) help',
    notes='Br–c:n. 2-Bromopyridine, 2-bromopyrimidine, 2-bromothiazole, 2-bromobenzothiazole. A narrower split of F013 (aryl bromide)'),
 'ewg_aryl_chloride': dict(
    reactive='the ring C–Cl carbon (ipso)',
    role='SNAr electrophile: an ortho or para electron-withdrawing group (NO2 ≫ CN, C(=O)R, SO2R, CF3) stabilises the Meisenheimer anion, so amines (K2CO3 or DIPEA, DMSO/DMF/NMP, 60–120 °C), alkoxides and thiolates displace Cl⁻. Pd-catalysed coupling is easier than for an electron-rich chloride (Suzuki: Pd(dppf)Cl2). A nitro group is reduced first (Fe/NH4Cl, SnCl2, H2/Pd-C leaves the Cl at rt)',
    acid=NONE, base=NONE, lg='Cl⁻ after addition (rate-limiting step is the addition)',
    steric='ortho substituents beyond the EWG slow SNAr and coupling; a second EWG (2,4-dinitro) makes it very reactive',
    notes='Cl–c with a nitro, cyano, carbonyl, CF3 or sulfonyl group ortho or para on the same ring. A narrower split of F047 (aryl chloride); the fluoride is F085'),
 'edg_aryl_chloride': dict(
    reactive='the ring C–Cl carbon',
    role='a poor SNAr electrophile (the donor destabilises the Meisenheimer anion) and a slow oxidative-addition partner: couplings need an electron-rich ligand (Pd2dba3 with XPhos, SPhos, RuPhos or t-Bu3P; Ni(cod)2 with a phosphine) at 100–110 °C; Grignard formation needs Mg/THF with an initiator; t-BuLi may deprotonate ortho to an OMe instead of exchanging the Cl; Pd/C–H2 does not remove it at rt',
    acid=NONE, base=NONE, lg='Cl⁻ only by special routes',
    steric='an ortho OMe or NR2 also coordinates Pd and can help or hinder',
    notes='Cl–c with an alkoxy, hydroxy or amino group ortho or para on the same ring (chloroanisole, chloroaniline). A narrower split of F047 (aryl chloride)'),
 'benzylic_chloride': dict(
    reactive='the CH2–Cl carbon',
    role='SN2 and SN1 electrophile, faster than a plain primary alkyl chloride (electron-rich aryl groups favour SN1 and solvolysis): amines (K2CO3 or Et3N, MeCN or DMF, rt–60 °C), alkoxides and phenoxides (Williamson), cyanide (NaCN, DMSO), azide (NaN3, DMF), thiols, malonate (NaH, DMF); catalytic NaI (acetone) speeds it up (Finkelstein). Arbuzov (P(OEt)3, heat) gives the benzylphosphonate; Pd(0) couples it with ArB(OH)2 to diarylmethanes. Lachrymatory, hydrolyses slowly in hot water',
    acid=NONE, base=NONE, lg='Cl⁻ (pKa of HCl ≈ −7): good',
    steric='primary: open; an ortho substituent slows SN2 slightly',
    notes='c–CH2–Cl. A narrower split of F011 (alkyl chloride) and F044 (1° alkyl halide); the bromide is F036'),
 'azole_halide_5ring': dict(
    reactive='the ring C–X carbon (C2 of the azole)',
    role='SNAr electrophile at C2 (2-chlorobenzothiazole and 2-chlorobenzoxazole react with amines, alkoxides and thiolates at rt–80 °C; 2-chloroimidazoles need protection of the N–H and heat); Pd coupling (Suzuki, Sonogashira, Buchwald) with standard catalysts; halogen–metal exchange (n-BuLi, −78 °C) at C2 gives the 2-metallo azole, which also forms by deprotonating the C2–H of the unsubstituted azole',
    acid=NONE, base=NONE, lg='Cl⁻, Br⁻ or I⁻ after addition',
    steric='a 4-substituent next to the ring N slows attack; N–H azoles need SEM, Boc or THP protection before metalation',
    notes='Cl, Br or I on a five-membered ring carbon next to a ring N (2-chlorothiazole, 2-bromoimidazole, 2-chlorobenzothiazole). A narrower split of F047, F013 and F084 (aryl halides); the six-membered rings are F086'),

 # ---------------- amines: nucleophilicity and basicity depend on what the N sits on
 'alkylamine_1': dict(
    reactive='the N (lone pair)',
    role='the strongest common nucleophile among primary amines: amide coupling (acid chlorides with Et3N or DIPEA at 0 °C; HATU, EDC/HOBt or T3P with DIPEA at rt), reductive amination (aldehyde, NaBH(OAc)3, AcOH, DCE; or NaBH3CN, MeOH), SNAr with halo-azines (DIPEA, n-BuOH or NMP, 80–150 °C), sulfonylation (RSO2Cl, pyridine or Et3N), carbamate and urea formation (Boc2O; isocyanates), SN2 alkylation (over-alkylation is common: use excess amine or Cs2CO3, DMF), aza-Michael addition, imine formation with aldehydes and ketones',
    acid='≈ 35 (N–H; only very strong bases remove it)', base='≈ 10.6 (RNH3⁺; ethylammonium 10.7, benzylammonium 9.3, 2,2,2-trifluoroethylammonium 5.7)', lg=NONE,
    steric='low: NH2 is the least hindered amine; branching at the α carbon slows acylation and SN2, and a tert-alkyl amine is a poor nucleophile',
    notes='NH2 on an sp3 carbon (not an aniline, amide or heteroaryl amine). A narrower split of F002 (amine 1°), whose SMARTS also matches anilines, which have their own fragments (F026, F115)'),
 'cyclic_amine_2': dict(
    reactive='the N–H nitrogen',
    role='a better nucleophile than an acyclic dialkylamine of similar basicity (the alkyl groups are tied back): SNAr with halo-arenes and halo-azines (K2CO3 or DIPEA, DMF/DMSO/NMP, 80–140 °C), Buchwald–Hartwig N-arylation (Pd2dba3 with RuPhos, NaOt-Bu), acylation and sulfonylation (rt), reductive amination (NaBH(OAc)3), alkylation (K2CO3, MeCN), carbamate protection (Boc2O, Cbz-Cl). Morpholine is less basic and less nucleophilic; piperazine reacts at both N unless one is protected',
    acid='≈ 36 (N–H)', base='≈ 11.1 (piperidinium 11.1; pyrrolidinium 11.3; morpholinium 8.4; piperazinium 9.7)', lg=NONE,
    steric='low to moderate; 2-substituents on the ring slow acylation and arylation',
    ring='low: 5- and 6-membered rings',
    notes='N–H in a ring with two sp3 carbons: piperidine, pyrrolidine, morpholine, piperazine, azepane. A narrower split of F034 (amine 2°)'),
 'acyclic_amine_2': dict(
    reactive='the N–H nitrogen',
    role='nucleophile and base: acylation, sulfonylation, SNAr, reductive amination, Mannich reactions and enamine formation with ketones and aldehydes. More hindered than a ring amine of the same size, so diisopropylamine is the precursor of the non-nucleophilic base LDA (n-BuLi, THF, −78 °C), and diethylamine is used as a volatile nucleophile',
    acid='≈ 36 (N–H)', base='≈ 10.9 (Et2NH2⁺ 10.9; i-Pr2NH2⁺ 11.1; N-methylbenzylammonium 9.6)', lg=NONE,
    steric='moderate; two α-branches (diisopropylamine) make it non-nucleophilic',
    notes='R2N–H with both substituents acyclic sp3 carbon. A narrower split of F034 (amine 2°); the ring amines are a separate fragment'),
 'trialkylamine_acyclic': dict(
    reactive='the N lone pair',
    role='an amine base and acid scavenger for acylations, sulfonylations, SN2 and SNAr (Et3N pKaH 10.75; DIPEA ≈ 11, hindered and non-nucleophilic); a nucleophilic catalyst when unhindered (DABCO, Et3N; DMAP is stronger). With alkyl halides it quaternises (Menshutkin) and the ammonium salt can undergo Hofmann elimination; oxidised to the N-oxide (mCPBA, H2O2); N-dealkylated by chloroformates (ACE-Cl, then MeOH)',
    acid='no N–H (α C–H ≈ 40)', base='≈ 10.75 (Et3NH⁺; DIPEA-H⁺ ≈ 11; N,N-dimethylbenzylamine 9.0)', lg=NONE,
    steric='DIPEA is shielded (non-nucleophilic); Et3N moderately; Me2NR is open',
    notes='N(CX4)(CX4)CX4 with no ring bond on the three carbons: Et3N, DIPEA, Me2N–R, Bu3N, DMAE. Recorded reactions list these as starting materials mainly when they act as a base, so the group counts as consumed when it is protonated. A narrower split of F093 (amine 3°)'),
 'dialkylaniline': dict(
    reactive='the para carbon of the ring and the N lone pair',
    role='an electron-rich arene: electrophilic substitution para to the amine (Vilsmeier–Haack formylation POCl3/DMF, NBS or Br2 bromination, azo coupling with diazonium salts, nitration only with care); a weak base and nucleophile (pKaH 5.1 for N,N-dimethylaniline); oxidised by mCPBA to the N-oxide or by one-electron oxidants to the radical cation; N-dealkylated by chloroformates (von Braun); directed ortho-lithiation (n-BuLi, TMEDA). In arylpiperazines and arylpiperidines the N is usually a spectator',
    acid='no N–H (α C–H ≈ 40)', base='≈ 5.1 (PhNHMe2⁺ 5.1; N-phenylpiperidinium 5.2; 4-nitro analogue ≈ 0.6)', lg=NONE,
    steric='an ortho substituent twists the NR2 out of conjugation and raises the basicity and nucleophilicity of N',
    notes='Ar–N(C)C with N bonded to an aromatic carbon and two sp3 carbons (N,N-dimethylaniline, N-arylpiperazine, N-arylpiperidine, N-aryl morpholine). A narrower split of F093 (amine 3°)'),
 'azine_secondary_amine': dict(
    reactive='the N–H nitrogen and the ring N',
    role='a weak N-nucleophile, because the lone pair is delocalised into the ring (amidine-like): acylation needs acid chlorides, or HATU with heat or a stronger base (NaH, LiHMDS); Buchwald–Hartwig N-arylation (Pd2dba3, Xantphos, Cs2CO3); alkyl halides react at the ring N unless the N–H is deprotonated first (NaH, DMF). Made by SNAr of a 2-halo-azine with a primary amine (DIPEA, 100–140 °C) or by Buchwald coupling',
    acid='≈ 26 (N–H; 2-(methylamino)pyridine, estimate in DMSO)', base='≈ 7 (2-(alkylamino)pyridinium 7.0; the ring N is the basic site)', lg=NONE,
    steric='a 3-substituent on the ring slows acylation and arylation',
    notes='R–NH–c:n. 2-(Alkylamino)pyridine, 2-(alkylamino)pyrimidine, 2-(alkylamino)thiazole, aminopurines. A narrower split of F094 (secondary aryl amine); the NH2 versions are F115'),
 'aminoazole_nh2': dict(
    reactive='the NH2 nitrogen and the ring N',
    role='an amidine-like, weakly nucleophilic NH2: acylation with acid chlorides (pyridine or Et3N, DCM, 0 °C → rt) or HATU, urea formation with isocyanates, N-arylation (Buchwald: Pd2dba3/Xantphos), diazotisation to the 2-halo-azole (t-BuONO, CuBr2), Hantzsch-type condensation with α-halo ketones to fused imidazoles; the ring N is the more nucleophilic site for alkyl halides; Boc protection needs Boc2O with DMAP',
    acid='≈ 26 (N–H, estimate in DMSO)', base='≈ 5.4 (2-aminothiazolium 5.4; 2-aminobenzimidazolium ≈ 7.5; 2-aminoimidazolium ≈ 8.5)', lg=NONE,
    steric='a 4-substituent next to the ring N slows acylation',
    ring='none: not a ring fragment',
    notes='NH2 on a five-membered ring carbon next to a ring N (2-aminothiazole, 2-aminoimidazole, 2-aminobenzimidazole, 2-aminooxazole). A narrower split of F115'),
 'aminoazine_nh2': dict(
    reactive='the NH2 nitrogen and the ring N',
    role='a weak nucleophile (amidine-like): acylation and sulfonylation need acid or sulfonyl chlorides with pyridine, or HATU/T3P with DIPEA and heat (bis-acylation is a side reaction); Buchwald–Hartwig coupling (Pd2dba3, Xantphos, Cs2CO3, dioxane, 100 °C); as the nucleophile in SNAr it needs the anion (NaH or LiHMDS); diazotisation to the 2-halo or 2-hydroxy azine; condensation with chloroacetaldehyde gives imidazo[1,2-a]pyridine through the ring N; electrophilic halogenation at C5 (NBS, MeCN)',
    acid='≈ 28 (N–H; 2-aminopyridine, DMSO)', base='≈ 6.9 (2-aminopyridinium 6.9; 2-aminopyrimidinium 3.5; 2-aminopyrazinium 3.1)', lg=NONE,
    steric='a 3-substituent next to the NH2 slows acylation and arylation',
    notes='NH2 on a six-membered ring carbon next to a ring N (2-aminopyridine, 2-aminopyrimidine, 2-aminopyrazine). A narrower split of F115'),

 # ---------------- amides, imides, lactams
 'secondary_amide': dict(
    reactive='the carbonyl C (hydrolysis, reduction), the N–H and the N',
    role='the most robust carbonyl group: stable to base, mild acid and most nucleophiles; hydrolysed only by hot 6 M HCl or NaOH (reflux, hours); reduced by LiAlH4 (THF, reflux) or BH3·THF to the amine; N–H is a weak acid, deprotonated by NaH or LDA for N-alkylation (NaH, RX, DMF) or N-arylation (Goldberg: CuI, DMEDA, K3PO4; Pd/Xantphos); O-alkylation with Meerwein salt (Me3O⁺ BF4⁻) gives the imidate; Tf2O/2-chloropyridine activates it for nucleophiles (Movassaghi); POCl3 gives the imidoyl chloride (Bischler–Napieralski with an N-phenethyl group); a directing group for ortho-lithiation',
    acid='≈ 17 (N–H; 25.9 in DMSO for N-methylacetamide)', base='≈ −0.5 (O-protonated amide)', lg='NR⁻: very poor',
    steric='N-alkyl groups, especially tert-alkyl, slow N-functionalisation; an α-quaternary centre protects the carbonyl',
    notes='R–C(=O)–NH–R with carbon on the carbonyl. A narrower split of F021 (amide) and F095 (amide N–H)'),
 'tertiary_amide': dict(
    reactive='the carbonyl C (hydrolysis, reduction, addition)',
    role='no N–H. The Weinreb amide (N-methoxy-N-methyl) gives a ketone with RMgX or RLi (one addition, a stable chelate) and an aldehyde with LiAlH4 or DIBAL-H. Other tertiary amides are reduced by LiAlH4 (THF, reflux) or BH3 to the amine, or by Schwartz reagent to the aldehyde; hydrolysed only under forcing conditions (6 M HCl or NaOH, reflux). DMF with POCl3 is the Vilsmeier reagent and with ArLi gives the aldehyde (Bouveault); DMF, DMA and NMP are polar solvents. Tf2O/2-fluoropyridine activation allows addition of nucleophiles; N,N-diethylbenzamides undergo ortho-lithiation (s-BuLi, TMEDA)',
    acid='no N–H (α C–H ≈ 30)', base='≈ −0.3 (O-protonated DMF; N,N-dimethylacetamide −0.3)', lg='NR2⁻: very poor',
    steric='N,N-dialkyl groups shield the carbonyl; an α-quaternary centre makes hydrolysis very slow',
    notes='R–C(=O)–NR2 (carbon or H on the carbonyl, two non-H substituents on N; not a carbamate or urea). A narrower split of F021 (amide)'),
 'lactam': dict(
    reactive='the carbonyl C (ring opening, reduction) and the N',
    role='ring amides are stable to nucleophiles; ring-opening hydrolysis needs hot strong acid or base; N–H lactams are alkylated (NaH, RX, DMF) or acylated (Boc2O, DMAP: the N-Boc lactam opens with LiOH or amines); N-substituted lactams are reduced by LiAlH4 or BH3 to the cyclic amine; Lawesson reagent gives the thiolactam; Meerwein salt gives the lactim ether; the α C–H forms the lactam enolate (LDA, THF, −78 °C; the N is protected)',
    acid='≈ 17 (N–H; 2-pyrrolidinone 24.2 in DMSO)', base='≈ −0.2 (O-protonated 2-pyrrolidinone)', lg='NR⁻: very poor',
    steric='moderate; ring substituents next to the N or the carbonyl slow acylation and hydrolysis',
    ring='low for 5- to 7-membered lactams (a β-lactam is strained and reacts very differently: it is not included)',
    notes='Cyclic amide in a 5-, 6- or 7-membered ring: pyrrolidinone, piperidinone, caprolactam, oxindole, N-aryl lactams. A narrower split of F021 (amide)'),
 'ring5_nh_amide': dict(
    reactive='the N–H (an acid) and the carbonyl C',
    role='an acidic N–H: deprotonated by K2CO3, Cs2CO3 or NaH for N-alkylation (RX, DMF), N-arylation (Chan–Lam: Cu(OAc)2, ArB(OH)2; Buchwald or CuI/DMEDA) and N-acylation (the Evans oxazolidinone: n-BuLi, or Et3N/DMAP, then RCOCl); thiazolidine-2,4-dione condenses at C5 with aldehydes (Knoevenagel: piperidine, AcOH); hydantoins open with aqueous base to the amino acid',
    acid='≈ 17 (N–H of a γ-lactam; oxazolidinone ≈ 20 in DMSO, hydantoin 9.1, thiazolidine-2,4-dione 6.5)', base='≈ −0.5 (O-protonated amide)', lg='NR⁻: poor',
    steric='open at the N–H; a 4-substituent (Evans auxiliary) controls the face of the N-acyl enolate',
    ring='low: 5-membered ring',
    notes='N–H and C=O in the same five-membered ring: 2-pyrrolidinone, 2-oxazolidinone, hydantoin, thiazolidinedione, imidazolidinone. A narrower split of F095 (amide N–H)'),
 'imide_nh': dict(
    reactive='the N–H (acid) and the carbonyl C',
    role='an acidic N–H (pKa 6–11): the anion alkylates cleanly (K2CO3, RX, DMF; or Mitsunobu with an alcohol: PPh3, DIAD), the Gabriel route to primary amines; N-halogenation gives NBS and NCS (Br2 or Cl2 with NaOH); N-hydroxylation gives NHS; NaOH opens the ring to the amic acid; one C=O is reduced selectively (NaBH4 gives the hydroxylactam, Et3SiH/TFA the lactam)',
    acid='8.3 (phthalimide; succinimide 9.6; glutarimide 11.4; thiazolidinedione 6.5)', base='≈ −2 (O-protonated imide)', lg='the imide N⁻ (pKa 8–10) leaves from N-acyl and N-sulfonyl imides',
    steric='open at the N–H',
    notes='C(=O)–NH–C(=O): phthalimide, succinimide, glutarimide, maleimide, thiazolidinedione with a free N–H. A narrower split of F121 (imide)'),
 'imide_nsub': dict(
    reactive='the carbonyl C (ring opening) and the C–N bond (deprotection)',
    role='a masked amine: the N-alkylphthalimide releases the primary amine with hydrazine hydrate (EtOH, reflux; Ing–Manske), MeNH2 (EtOH, rt) or ethylenediamine; NaBH4 then AcOH (Osby) is a mild alternative; NaOH or HCl hydrolysis is slow. N-Hydroxyphthalimide esters are redox-active esters (Ni, Zn or photoredox decarboxylative coupling); N-acyl imides are active acyl donors (LiOH/H2O2 cleaves an Evans auxiliary)',
    acid='no N–H (α C–H of a succinimide ≈ 25)', base='≈ −2 (O-protonated imide)', lg='the phthalimide or succinimide N⁻ leaves from the N-acyl form',
    steric='the N-substituent slows the ring opening; tert-alkyl groups are very slow to deprotect',
    notes='C(=O)–N(C, N or O)–C(=O): N-alkyl or N-aryl phthalimide, succinimide, N-hydroxyphthalimide, N-acyl imides (not N-halo: F151). A narrower split of F121 (imide)'),

 # ---------------- unsaturated carbonyls and nitriles
 'enone_acyclic': dict(
    reactive='the β carbon (conjugate addition) and the carbonyl C (1,2-addition)',
    role='a strong Michael acceptor: soft nucleophiles add at the β carbon (thiols with Et3N; amines, aza-Michael, often uncatalysed; malonate and nitroalkanes with DBU or K2CO3; organocuprates R2CuLi, THF, −78 °C); a dienophile (Lewis-acid-catalysed Diels–Alder); hard nucleophiles (RLi, LiAlH4, NaBH4 with CeCl3: Luche) add 1,2 to the carbonyl; nucleophilic epoxidation (H2O2, NaOH; Weitz–Scheffer); selective C=C reduction (H2, Pd/C; NaBH4/NiCl2; L-Selectride; Stryker reagent); Nazarov cyclisation and Robinson annulation. It polymerises (methyl vinyl ketone) and adds water slowly',
    acid=NONE, base='≈ −3 (O-protonated enone)', lg=NONE,
    steric='α- and β-substituents slow conjugate addition; a β,β-disubstituted enone is far less reactive',
    notes='C=C–C(=O)–C with the C=C outside a ring: chalcone, methyl vinyl ketone, benzylideneacetone, mesityl oxide. A narrower split of F040 (enone)'),
 'enone_cyclic': dict(
    reactive='the β carbon (conjugate addition) and the carbonyl C',
    role='a Michael acceptor, less reactive than an acyclic enone but readily attacked by cuprates (R2CuLi with TMSCl: the classic 1,4-addition to cyclohexenone), thiols and amines under mild catalysis; Luche reduction (NaBH4, CeCl3, MeOH) gives the allylic alcohol; hydrogenation (H2, Pd/C) reduces the C=C first; nucleophilic epoxidation (H2O2, NaOH); Diels–Alder dienophile; photochemical [2+2]; the product of a Robinson annulation. Ring substituents at C4 and C5 direct the face of attack',
    acid=NONE, base='≈ −3 (O-protonated enone)', lg=NONE,
    steric='ring conformation and 4-/5-substituents direct the face of attack; a 3-substituent (β) slows conjugate addition',
    ring='low: 5- and 6-membered rings',
    notes='C=C–C(=O) all in one ring: cyclohexenone, cyclopentenone, chromone, steroid 4-en-3-ones. A narrower split of F040 (enone)'),
 'aryl_nitrile': dict(
    reactive='the nitrile C (electrophile)',
    role='an electrophile at carbon, deactivated by conjugation: hydrolysed to the amide (H2O2/K2CO3, DMSO; or conc. H2SO4, heat) and the acid (6 M HCl or NaOH, reflux); reduced to the benzylamine (LiAlH4, THF; H2/Raney Ni with NH3; BH3·THF; CoCl2/NaBH4) or, with DIBAL-H (toluene, −78 °C), to the aldehyde; Grignard or organolithium addition gives an aryl ketone after hydrolysis; [3+2] with NaN3 (NH4Cl, DMF, 120 °C, or ZnBr2) gives the tetrazole; Pinner reaction (HCl, ROH) gives the imidate; hydroxylamine gives the amidoxime. As a strong EWG it activates ortho/para halides for SNAr and survives most Pd couplings',
    acid=NONE, base='≈ −10 (N-protonated nitrile)', lg=NONE,
    steric='ortho substituents slow additions',
    notes='Ar–C≡N: benzonitrile, cyanopyridine, cyanoarenes. A narrower split of F022 (nitrile)'),
 'ketone_alpha_ch': dict(
    reactive='the α carbon and its C–H',
    role='forms the enol (acid) or enolate (base): aldol (NaOH/EtOH, or LDA, THF, −78 °C), α-alkylation (LDA kinetic at the less substituted side; NaH or KOt-Bu thermodynamic at the more substituted side; then R–X), α-halogenation (Br2/AcOH, NBS/TsOH, CuBr2), Mannich (HCHO, R2NH), Michael addition, α-arylation (Pd, Cs2CO3), Claisen–Schmidt condensation with ArCHO, Stork enamine alkylation, silyl enol ether formation (TMSCl with Et3N, or LDA then TMSCl)',
    acid='≈ 19 (acetone 19.3; acetophenone 18.3; 26.5 and 24.7 in DMSO)', base='≈ −7 (O-protonated ketone)', lg=NONE,
    steric='a bulky base (LDA) takes the less substituted α C–H first; the more substituted gives the more stable enol or enolate',
    notes='C–H on the carbon next to a ketone C=O (not an α-halo carbon, not between two carbonyls: F102). A narrower split of F101 (enolisable α-C–H)'),
 'ester_alpha_ch': dict(
    reactive='the α carbon and its C–H',
    role='needs a strong base (LDA or LiHMDS, THF, −78 °C; NaH or KOt-Bu in DMF) for the enolate: alkylation (R–I or R–Br), aldol, Claisen condensation (NaOEt in EtOH, matching alkoxide), Michael addition, α-hydroxylation (Davis oxaziridine), α-halogenation (LDA then NBS or I2), α-arylation (Pd, LiHMDS), Ireland–Claisen (silyl ketene acetal). Hydroxide or a non-matching alkoxide hydrolyses or transesterifies the ester instead',
    acid='≈ 25 (ethyl acetate 25.6; 29.5 in DMSO)', base='≈ −6.5 (O-protonated ester)', lg=NONE,
    steric='a bulky base takes the less substituted α C–H; an α,α-disubstituted ester forms the enolate slowly',
    notes='C–H on the carbon next to an ester or lactone C=O (not between two carbonyls: F102). A narrower split of F101 (enolisable α-C–H)'),
 'aldehyde_alpha_ch': dict(
    reactive='the α carbon and its C–H',
    role='the most acidic carbonyl C–H: enolises under mild conditions, so aldol condensations run with catalytic base (NaOH, K2CO3) or amine catalysis (proline or pyrrolidine through the enamine: Mannich, α-functionalisation, α-halogenation with NCS or NBS); self-condensation and polymerisation are the side reactions, so enolisable aldehydes are added slowly; α-stereocentres epimerise easily',
    acid='≈ 17 (acetaldehyde 16.7)', base='≈ −7 (O-protonated aldehyde)', lg=NONE,
    steric='open; α-branching slows enolisation',
    notes='C–H on the carbon next to an aldehyde C=O. A narrower split of F101 (enolisable α-C–H)'),
 'nitrile_alpha_ch': dict(
    reactive='the α carbon and its C–H',
    role='needs LDA, LiHMDS, KOt-Bu or NaH for the cyanomethyl anion: alkylation (LDA, THF, −78 °C, R–X), Thorpe self-condensation, conjugate addition, α-arylation (Pd/ZnF2, or SNAr with aryl fluorides and KHMDS). Benzylic nitriles (ArCH2CN, pKa 21.9 in DMSO) alkylate under phase-transfer conditions (50 % NaOH, Bu4N⁺ Br⁻, R–X)',
    acid='≈ 25 (acetonitrile; 31.3 in DMSO; phenylacetonitrile 21.9 in DMSO)', base='≈ −10 (N-protonated nitrile)', lg=NONE,
    steric='open; α-branching slows alkylation',
    notes='C–H on the carbon next to a nitrile C≡N (not between two carbonyls or nitriles: F102). A narrower split of F101 (enolisable α-C–H)'),

 # ---------------- imines, enamines
 'imine_acyclic': dict(
    reactive='the C of the C=N (electrophile) and the N',
    role='an electrophile at carbon, reduced to the amine (NaBH4, NaBH3CN, NaBH(OAc)3 in reductive amination; H2/Pd; LiAlH4); Strecker (TMSCN), Mannich and Grignard/organolithium addition; hydrolysed back to the carbonyl in aqueous acid. Oximes undergo Beckmann rearrangement (TsCl, PCl5, H2SO4), reduction to the amine (LiAlH4, H2/Raney Ni, Zn/AcOH) and, as aldoximes, dehydration to the nitrile; hydrazones give Wolff–Kishner and Shapiro reactions; Ellman N-sulfinyl imines add RMgX diastereoselectively; the α C–H forms the aza-enolate (LDA)',
    acid='≈ 11 (oxime O–H: acetone oxime 12.4; benzaldoxime 10.7; an N-alkyl imine has no acidic H)', base='≈ 7 (N-alkyl imine; oxime ≈ 1; hydrazone ≈ 3)', lg=NONE,
    steric='ketimines are more hindered than aldimines; N-aryl imines are less basic and more hydrolytically labile',
    notes='C=N with the nitrogen outside a ring: Schiff bases, oximes, hydrazones, N-sulfinyl imines. A narrower split of F031 (imine)'),
 'imine_aldimine': dict(
    reactive='the C of the C=N (electrophile) and the N',
    role='more electrophilic and less hindered than a ketimine: reduced and attacked under mild conditions (NaBH3CN at pH 6; Strecker; Mannich; allyl- and crotylboronates); formed reversibly from aldehydes and amines (MgSO4 or 4 Å sieves, or Dean–Stark) and hydrolysed easily; aldoximes dehydrate (Ac2O, Burgess reagent, SOCl2) to nitriles and are oxidised (NCS/pyridine, NaOCl) to nitrile oxides for [3+2] cycloadditions; with an α C–H it tautomerises to the enamine',
    acid='≈ 11 (aldoxime O–H 10.7–11.3)', base='≈ 7 (aldimine; aldoxime ≈ 1)', lg=NONE,
    steric='open',
    notes='H–C=N: aldimines, aldoximes, aldehyde hydrazones (formaldehyde-derived CH2=N is not included). A narrower split of F031 (imine)'),
 'enamine_acyclic': dict(
    reactive='the β carbon (nucleophile) and, in an enaminone, the N-bearing carbon (electrophile)',
    role='nucleophilic at the β carbon (Stork enamine alkylation and acylation, then H3O⁺ hydrolysis to the carbonyl); the enaminone from a methyl ketone and DMF-DMA (toluene or neat, 80–110 °C) is the standard precursor of heterocycles: hydrazines give pyrazoles, amidines and guanidines give pyrimidines, hydroxylamine gives isoxazoles, amines exchange the Me2N group by addition–elimination; hydrolysed by aqueous acid; reduced by NaBH3CN',
    acid=NONE, base='≈ 9 (simple enamine; an enaminone ≈ 3 on O)', lg='Me2NH after addition–elimination',
    steric='open; the (E)-isomer dominates',
    notes='N–CH=C with the carbon bonded to N outside a ring: enamines of acyclic ketones, DMF-DMA enaminones (Me2N–CH=CH–C(=O)R), push–pull alkenes. A narrower split of F127 (enamine / enaminone)'),

 # ---------------- alcohols, orthoesters, amide acetals
 'methanol': dict(
    reactive='the O–H (nucleophile) and, activated, the CH3–O bond',
    role='solvent and reagent: Fischer esterification (H2SO4 or HCl, reflux; SOCl2 or TMSCl in MeOH; DCC/DMAP or EDC for acids), transesterification, methanolysis of acid chlorides and anhydrides, methyl ether formation (NaH then MeI; Mitsunobu with phenols), dimethyl acetal and ketal formation (TsOH or CSA, with HC(OMe)3), hydrogenation and NaBH4 reductions run in it. In base it is NaOMe; it quenches Grignard and organolithium reagents',
    acid='15.5 (methanol; 29.0 in DMSO)', base='≈ −2.2 (CH3OH2⁺)', lg='MeOH2⁺ leaves from the protonated alcohol (SN2 at CH3)',
    steric='the smallest alcohol: the most nucleophilic towards hindered electrophiles',
    notes='CH3–OH. Recorded reactions list it as a starting material mainly when it is the reagent. A narrower split of F014 (alcohol on sp3 C)'),
 'ethanol': dict(
    reactive='the O–H (nucleophile) and, activated, the CH3CH2–O bond',
    role='solvent and reagent: Fischer esterification (H2SO4 or HCl, reflux), transesterification, ethanolysis of acid chlorides and anhydrides, ethyl ether formation (NaH then EtI), diethyl acetal and ketal formation (TsOH, HC(OEt)3), NaOEt for Claisen and malonic ester syntheses (matching the ester avoids transesterification), a protic solvent for NaBH4 reductions and hydrogenations. It quenches Grignard and organolithium reagents',
    acid='15.9 (ethanol; 29.8 in DMSO)', base='≈ −2.4 (CH3CH2OH2⁺)', lg='EtOH2⁺ leaves from the protonated alcohol',
    steric='small: slightly less nucleophilic than methanol',
    notes='CH3–CH2–OH. Recorded reactions list it as a starting material mainly when it is the reagent. A narrower split of F003 (alcohol 1°)'),
 'benzyl_alcohol': dict(
    reactive='the O–H (nucleophile), the benzylic C–H (oxidation) and the C–O bond (activation)',
    role='oxidised to the aldehyde (MnO2, selective for benzylic and allylic; Dess–Martin periodinane; Swern; TEMPO/PhI(OAc)2; Cu/TEMPO/air) and to the acid (Jones, KMnO4); converted to the benzylic chloride or bromide readily (SOCl2, PBr3, MsCl then LiBr, CBr4/PPh3; fast and SN1-like for electron-rich aryl groups); Mitsunobu; ether and ester formation (NaH, RX; RCOCl, Et3N); hydrogenolysis of the C–O bond (H2, Pd/C) for neutral and electron-rich aryl groups; acid-catalysed Friedel–Crafts alkylation through the benzylic cation',
    acid='15.4 (benzyl alcohol; 27.0 in DMSO)', base='≈ −2 (ArCH2OH2⁺)', lg='H2O after protonation or activation (OMs, OTs, OTf); benzylic cation stabilised by the ring',
    steric='open (primary)',
    notes='Ar–CH2–OH. A narrower split of F130 (benzyl group on O, N or S) and of F003 (alcohol 1°)'),
 'orthoester': dict(
    reactive='the central carbon (electrophile) and the C–O bonds',
    role='a reagent: a water scavenger and a source of an ester, acetal or ethoxymethylene group. Acid-catalysed (TsOH, CSA) formation of acetals and ketals from carbonyls (HC(OMe)3, MeOH); esters from acids (heat); the Johnson–Claisen rearrangement with allylic alcohols (CH3C(OEt)3, cat. propionic acid, 140 °C); ethoxymethylene derivatives of active methylenes (HC(OEt)3, Ac2O); amidines and heterocycles with amines (with NaN3 gives tetrazoles; with 2-aminobenzamides gives quinazolinones). Hydrolyses at once in aqueous acid to the ester and alcohol',
    acid=NONE, base=NONE, lg='ROH after protonation (a stabilised dialkoxycarbenium ion)',
    steric='open',
    notes='C(OR)3 on carbon: triethyl orthoformate, trimethyl orthoacetate, trimethyl orthobenzoate. A reagent: it is consumed in nearly every recorded reaction that lists it. A narrower split of F039 (acetal / ketal)'),
 'amide_acetal': dict(
    reactive='the central carbon (electrophile) and the OMe groups (methylation)',
    role='a formylating and methylating reagent: with active methylene groups and methyl ketones it gives the enaminone (toluene or neat, 80–110 °C), the usual route to pyrazoles, pyrimidines, isoxazoles and indoles (Leimgruber–Batcho: o-nitrotoluene, DMF-DMA, pyrrolidine, then reduction); it methylates carboxylic acids, phenols, thiols and N–H heterocycles on heating; with primary amides it gives the N-acyl formamidine; transaminates with amines; hydrolyses to DMF and methanol in water. Bredereck reagent (t-BuO–CH(NMe2)2) is the stronger, more hindered relative',
    acid=NONE, base=NONE, lg='MeOH or Me2NH',
    steric='open',
    notes='C(OR)(OR)N: N,N-dimethylformamide dimethyl acetal (DMF-DMA), DMF diethyl acetal, Bredereck reagent. A reagent: it is consumed in nearly every recorded reaction that lists it. A narrower split of F039 (acetal / ketal)'),

 # ---------------- sulfur
 'azine_sulfone': dict(
    reactive='the ring carbon bearing the sulfonyl or sulfinyl group (attacked at the ipso position)',
    role='an SNAr electrophile with an excellent leaving group (sulfinate), more reactive than the chloro analogue: amines (DIPEA, THF or MeCN, rt–60 °C), alkoxides, thiolates, cyanide and carbanions displace it at rt; 2-sulfonylpyrimidines arylate cysteine selectively at pH 7. Made by oxidising the methylthio-azine (mCPBA, DCM; Oxone, H2O/MeOH; H2O2 with Na2WO4). The other use of a heteroaryl sulfone is the Julia–Kocienski olefination (benzothiazolyl and phenyltetrazolyl sulfones, KHMDS, DME, −55 °C, with an aldehyde)',
    acid='no acidic H on the ring; the α C–H of the sulfone ≈ 29 (DMSO)', base='≈ −2 (the ring N; weakly basic)', lg='RSO2⁻ (sulfinate; pKa of the sulfinic acid ≈ 2): excellent',
    steric='a substituent next to the C–S slows attack',
    notes='R–S(=O)(=O)–c:n and R–S(=O)–c:n: 2-(methylsulfonyl)pyrimidine, 2-sulfonylpyridine, 2-(methylsulfonyl)benzothiazole. A narrower split of F071 (sulfone / sulfoxide)'),
}
