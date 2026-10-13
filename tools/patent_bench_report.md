# The Reactions page against 1,000 held-out US patent reactions

1000 recorded reactions, one from each of 1000 randomly chosen held-out US patents (seed 1; 4,956 held-out US patents, 49,913 of their reactions among the 1,090,843 learned). Made by tools/patent_bench.js sample.

Each reaction's starting materials, reagent and recorded solvents were typed into the code bar and the arrow clicked, with the exact look-up of recorded reactions turned off. Graded with RDKit (tools/patent_bench_grade.py). Median time per reaction: 729 ms.

## Overall

| | reactions | share of all |
|---|---|---|
| answered (showed products) | 573 | 57.3 % |
| declined (Unknown / No reaction) | 422 | 42.2 % |
| error | 5 | 0.5 % |

| recorded product among the page's products | of all | of those answered |
|---|---|---|
| exact (with stereo) | 443 (44.3 %) | 77.3 % |
| ignoring stereo | 453 (45.3 %) | 79.1 % |
| same compound (salt / protonated form counts) | 460 (46.0 %) | 80.3 % |

## By the part of the page that answered

| | reactions | answered | right, exact | right, same compound | right when answered (exact) |
|---|---|---|---|---|---|
| error | 5 | 0 | 0 (0.0 %) | 0 (0.0 %) | - |
| no reaction (workbook) | 42 | 0 | 0 (0.0 %) | 0 (0.0 %) | - |
| relevance model (learned) | 546 | 546 | 430 (78.8 %) | 445 (81.5 %) | 78.8 % |
| unknown | 380 | 0 | 0 (0.0 %) | 0 (0.0 %) | - |
| workbook rule | 27 | 27 | 13 (48.1 %) | 15 (55.6 %) | 48.1 % |

## By size of the starting materials (heavy atoms)

| | reactions | answered | right, exact | right, same compound | right when answered (exact) |
|---|---|---|---|---|---|
| < 8 | 4 | 0 | 0 (0.0 %) | 0 (0.0 %) | - |
| 8-11 | 47 | 24 | 15 (31.9 %) | 16 (34.0 %) | 62.5 % |
| 12-19 | 277 | 155 | 110 (39.7 %) | 122 (44.0 %) | 71.0 % |
| 20-29 | 369 | 216 | 168 (45.5 %) | 171 (46.3 %) | 77.8 % |
| 30 or more | 303 | 178 | 150 (49.5 %) | 151 (49.8 %) | 84.3 % |

## With and without a reagent

| | reactions | answered | right, exact | right, same compound | right when answered (exact) |
|---|---|---|---|---|---|
| reagent given | 813 | 492 | 374 (46.0 %) | 391 (48.1 %) | 76.0 % |
| no reagent | 187 | 81 | 69 (36.9 %) | 69 (36.9 %) | 85.2 % |

## Is the model's stated chance honest?

| stated chance | answers | mean stated | right (exact) | right (same compound) |
|---|---|---|---|---|
| 0.50-0.70 | 55 | 60 % | 42 % | 53 % |
| 0.70-0.80 | 49 | 75 % | 57 % | 61 % |
| 0.80-0.90 | 114 | 86 % | 76 % | 78 % |
| 0.90-0.95 | 100 | 93 % | 83 % | 85 % |
| 0.95-1.00 | 228 | 97 % | 92 % | 93 % |

## The answers it held back

When the page declines with "Unknown", its reading still names the relevance model's best guess and its stated chance (below the floor the page needs). 371 declined reactions had one; 78 of those guesses were right exactly (21.0 %), 83 the same compound.

| stated chance | held-back guesses | right (exact) |
|---|---|---|
| 0.00-0.30 | 195 | 18 (9.2 %) |
| 0.30-0.50 | 134 | 39 (29.1 %) |
| 0.50-0.70 | 30 | 14 (46.7 %) |
| 0.70-0.85 | 12 | 7 (58.3 %) |

## Why it declined (most common first lines)

* 37 × Unknown: nothing in reaction_rules.xlsx covers F001 carboxylic acid, F
* 15 × Unknown: nothing in reaction_rules.xlsx covers F002 amine (1°), F026 a
* 14 × Unknown: nothing in reaction_rules.xlsx covers F004 ester, F005 alkene
* 14 × Unknown: nothing in reaction_rules.xlsx covers F014 alcohol (on sp3 C)
* 14 × Unknown: nothing in reaction_rules.xlsx covers F035 aryl methyl (benzy
* 10 × No reaction under these conditions: R035 · Acid-catalysed dehydration 
* 9 × Unknown: nothing in reaction_rules.xlsx covers F003 alcohol (1°), F014
* 9 × Unknown: nothing in reaction_rules.xlsx covers F038 alkyl arene (aryl–
* 8 × Unknown: nothing in reaction_rules.xlsx covers F011 alkyl chloride, F0
* 8 × Unknown: nothing in reaction_rules.xlsx covers F005 alkene (C=C), F019
* 7 × Unknown: nothing in reaction_rules.xlsx covers F021 amide, F038 alkyl 
* 7 × Unknown: nothing in reaction_rules.xlsx covers F019 ketone, F038 alkyl

## Wrong answers (first 40: the page's product is not even the same compound)

| patent | starting materials > reagent | recorded | page | answered by |
|---|---|---|---|---|
| US06982266B2 | `COc1cc2c(N3CCN(C(=O)OC(C)(C)C)CC3)ncnc2cc1O.Cc1ccc(S(=O)(=O)OC(C)Cl)cc1>O=C([O-])[O-].[Cs+]> | cond=DMF` | `COc1cc2c(N3CCN(C(=O)OC(C)(C)C)CC3)ncnc2cc1OCCCl` | `COc1cc2c(N3CCN(C(=O)OC(C)(C)C)CC3)ncnc2cc1OC(C)Cl.Cc1ccc(S(=O)(=O)O)cc1` | relevance model (learned) |
| US06951867B2 | `CC(C)(C)OC(=O)N1CCC(CN)CC1.Cc1nc2c(C(=O)O)cc(Cl)c(N)n2c1C>>` | `Cc1cn2c(N)c(Cl)cc(C(=O)NCC3CCN(C(=O)OC(C)(C)C)CC3)c2n1` | `Cc1nc2c(C(=O)NCC3CCN(C(=O)OC(C)(C)C)CC3)cc(Cl)c(N)n2c1C` | relevance model (learned) |
| US08653269B2 | `Cc1ccccc1.Nc1cccc(CCl)n1.O=C(OC(Cl)(Cl)Cl)OC(Cl)(Cl)Cl>[Na+].[OH-]> | cond=tert-butanol, water` | `CC(C)(C)OC(=O)Nc1cccc(CCl)n1` | `Nc1cccc(CO)n1.[Cl-]` | workbook rule |
| US20140194430A1 | `Clc1ccc2ncc(-c3cc4ccccc4o3)n2n1.CNCC(O)c1cccnc1>Cl.[H-].[Na+]> | cond=DMF` | `CN(CC(O)c1cccnc1)c1ccc2ncc(-c3cc4ccccc4o3)n2n1` | `CNCC(Oc1ccc2ncc(-c3cc4ccccc4o3)n2n1)c1cccnc1.Cl` | relevance model (learned) |
| US04092474 | `CC(=O)OC(C1=CS[C@@H]2C(NC(=O)Cc3cccs3)C(=O)N2C1C(=O)OC(c1ccccc1)c1ccccc1)c1ccccc1>O=C(OO)c1cccc(Cl)c1> | cond=dichloromethane, isopropanol` | `CC(=O)OC(C1=C(C(=O)OC(c2ccccc2)c2ccccc2)N2C(=O)C(NC(=O)Cc3cccs3)[C@H]2SC1)c1ccccc1` | `CC(=O)OC(c1ccccc1)C12OC1S[C@@H]1C(NC(=O)Cc3cccs3)C(=O)N1C2C(=O)OC(c1ccccc1)c1ccccc1.O=C(O)c1cccc(Cl)c1` | workbook rule |
| US07157408B2 | `OCCc1c(Br)cc(Br)c2cccnc12.O=C(OC(Cl)(Cl)Cl)OC(Cl)(Cl)Cl>> | cond=1,4-dioxane` | `COC(=O)OCCc1c(Br)cc(Br)c2cccnc12` | `O=C(Cl)OCCc1c(Br)cc(Br)c2cccnc12.O=C=O.Cl` | relevance model (learned) |
| US08313923B2 | `CC(C)=CCCC(C)CC=O>OC[C@H]1O[C@@H]2O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@H]4[C@H](O)[C@@H](O)[C@@H](O[C@H]5[C@H](O)[C@@H](O)[C@@H](O[C@H]6[C@H](O)[C@@H](O)[C@@H](O[C@H]7[C@H](O)[C@@H](O)[C@@H](O[C@H]8[C@H](O)[C@@H](O)[C@@H](O[C@H]1[C@H](O)[C@H]2O)O[C@@H]8CO)O[C@@H]7CO)O[C@@H]6CO)O[C@@H]5CO)O[C@@H]4CO)O[C@@H]3CO>` | `CC(C)=CCCC(C)=CC=O` | `CC(C)=CCCC(C)CCO` | relevance model (learned) |
| US05312840 | `Cc1ccccc1NCCNc1ccccc1C.N#CBr>[Na+].[OH-]> | cond=ethanol` | `Cc1ccccc1N1CCN(c2ccccc2C)C1=N` | `Cc1ccccc1NCCN(C#N)c1ccccc1C.Br` | relevance model (learned) |
| US07585876B2 | `CC(C)(C)OC(=O)CC#N.O=[N+]([O-])c1cccnc1Cl>O=C([O-])[O-].[K+]> | cond=tetrahydrofuran` | `N#CCc1ncccc1[N+](=O)[O-]` | `CC(C)(C)OC(=O)C(C#N)c1ncccc1[N+](=O)[O-].Cl` | relevance model (learned) |
| US06930124B2 | `O=Cc1cccn1-c1cc(Cl)cc(Cl)c1.CC(=O)Nc1ccc2c(c1)CC(=O)N2>C1CCNCC1> | cond=ethanol` | `CC(=O)Nc1ccc2c(c1)C(c1cccn1-c1cc(Cl)cc(Cl)c1)C(=O)N2` | `CC(=O)Nc1ccc2c(c1)C(=Cc1cccn1-c1cc(Cl)cc(Cl)c1)C(=O)N2` | relevance model (learned) |
| US07071198B1 | `CC(=O)Cc1ccncc1.COc1ccc(C=O)cc1F>C1CCNCC1> | cond=toluene` | `C=CC(=O)Cc1cccnc1` | `COc1ccc(C=C(C(C)=O)c2ccncc2)cc1F` | relevance model (learned) |
| US20090232907A1 | `CCOC(=O)NC(=O)Nc1cc2c3c(cccc3c1)C(=O)N(CCN(C)C)C2=O>NC(N)=O.[Na+].[OH-]>` | `CN(C)CCN1C(=O)c2cccc3cc(NC(N)=O)cc(c23)C1=O` | `CN(C)CCN1C(=O)c2cccc3cc(NC(=O)NC(=O)O)cc(c23)C1=O.CCO` | relevance model (learned) |
| US07863305B2 | `CS(=O)(=O)OC(CCCl)c1ccccc1.N#Cc1cccc2[nH]ncc12>[H-].[Na+]> | cond=NMP` | `N#Cc1cccc2c1cnn2C(CCCl)c1ccccc1` | `CS(=O)(=O)OC(CCn1ncc2c(C#N)cccc21)c1ccccc1.Cl` | relevance model (learned) |
| US04806625 | `Cc1cc(C)c(O)c(C)c1>Br.BrBr.ClC(Cl)(Cl)Cl.[Br-]>` | `Cc1c(O)c(C)c(Br)c(C)c1Br` | `Cc1cc(C)c(Br)c(C)c1O.Br` | relevance model (learned) |
| US06159971 | `COc1ccc([N+](=O)[O-])c2c1CCC(C(=O)O)C2.Nc1ccc(C(=O)N2CCOCC2)cc1.CCN(CC)CC>O=S(Cl)Cl> | cond=dichloromethane, DMF, toluene, water` | `COc1ccc([N+](=O)[O-])c2c1CCC(C(=O)Nc1ccc(C(=O)N3CCOCC3)cc1)C2` | `COc1ccc([N+](=O)[O-])c2c1CCC(C(=O)Cl)C2.O=S=O.Cl` | workbook rule |
| US08691782B2 | `COc1ccc(C(Cl)(c2ccccc2)c2ccc(OC)cc2)cc1.O=C(OCC1c2ccccc2-c2ccccc21)N1CCC[C@H]1CO>c1ccncc1> | cond=chloroform, methanol` | `COc1ccc(C(C(=O)OCC2c3ccccc3-c3ccccc32)(c2ccccc2)c2ccc(OC)cc2)cc1` | `COc1ccc(C(OC[C@@H]2CCCN2C(=O)OCC2c3ccccc3-c3ccccc32)(c2ccccc2)c2ccc(OC)cc2)cc1.Cl` | relevance model (learned) |
| US05599833 | `CCCCS(=O)(=O)c1ccc(-c2sc3cc(S(=O)(=O)CCCC)ccc3c2C(=O)c2ccc(OCCN3CCCCC3)cc2)cc1>Cl> | cond=diethyl ether, ethyl acetate` | `CCCCSc1ccc2c(C(=O)c3ccc(OCCN4CCCCC4)cc3)c(-c3ccc(S(=O)(=O)CCCC)cc3)sc2c1` | `CCCCS(=O)(=O)c1ccc(-c2sc3cc(S(=O)(=O)CCCC)ccc3c2Cc2ccc(OCCN3CCCCC3)cc2)cc1` | relevance model (learned) |
| US04247550 | `NCC=CCc1cc(N)c(Oc2ccccc2)c(S(N)(=O)=O)c1.O=Cc1ccccc1>[BH4-].[Na+]> | cond=acetic acid, methanol, water` | `NCC=CCc1cc(NCc2ccccc2)c(Oc2ccccc2)c(S(N)(=O)=O)c1` | `OCc1ccccc1` | workbook rule |
| US08288552B2 | `COC(=O)c1c(Cl)ccc(N)c1Br.O=C([O-])[O-]>CB1OB(C)OB(C)O1.[Cl][Pd][Cl].[Fe+2].[K+].c1ccc(P(c2ccccc2)[c-]2cccc2)cc1> | cond=1,4-dioxane, water` | `COC(=O)c1cc(N)c(C)cc1Cl` | `COC(=O)c1c(Cl)ccc(N)c1C` | relevance model (learned) |
| US20030108892A1 | `ClCCc1cc2ccccc2c2ccccc12.O=C1NC(=O)c2ccccc21>CCCCCCCCCCCCCCCC[P+](CCCC)(CCCC)CCCC.[Br-].[K]> | cond=toluene, water` | `O=C1NC(=O)c2c(CCc3cc4ccccc4c4ccccc34)cccc21` | `O=C1c2ccccc2C(=O)N1CCc1cc2ccccc2c2ccccc12.Cl` | relevance model (learned) |
| US04954600 | `Oc1ccc(-c2ccc(O)cc2)cc1.C=CCOCCCCCCCCBr>[K+].[OH-]> | cond=methanol` | `C=CCOCCCCCCCCOc1ccc(-c2ccc(O)cc2)cc1` | `C=CCOCCCCCCCCO.[Br-]` | workbook rule |
| US07074789B2 | `Nc1nc2ccccc2[nH]1.Clc1ncnc(Cl)n1>CCN(C(C)C)C(C)C> | cond=DMF, water` | `Nc1nc2ccccc2n1-c1ncnc(Cl)n1` | `Clc1ncnc(Nc2nc3ccccc3[nH]2)n1.Cl` | relevance model (learned) |
| US06465444B2 | `O=Cc1ccc2ccc(CO)cc2c1.CS(=O)(=O)Cl>CCN(CC)CC.[Cl-].[NH4+]> | cond=dichloromethane` | `CS(=O)(=O)OCc1ccc2ccc(C=O)cc2c1` | `O=Cc1ccc2ccc(CCl)cc2c1.CS(=O)(=O)O` | relevance model (learned) |
| US05583149 | `CNS(=O)(=O)Cc1ccc(N)c(I)c1.C[Si](C)(C)C#C[Si](C)(C)C>O=C([O-])[O-].[Na+]> | cond=1,2-dimethoxyethane` | `CNS(=O)(=O)Cc1ccc2[nH]c([Si](C)(C)C)cc2c1` | `CNS(=O)(=O)Cc1ccc2[nH]c([Si](C)(C)C)c([Si](C)(C)C)c2c1.I` | relevance model (learned) |
| US09181514B2 | `C=C(C)CC1(C)CCCC(C)C1=O>COCC[O][AlH2-][O]CCOC.[Na+]> | cond=toluene` | `CC1CCCC2(C)CC(C)(C)OC12` | `C=C(C)CC1(C)CCCC(C)[C@@H]1O` | relevance model (learned) |
| US06251927B1 | `CC(=O)N1C(C(=O)O)CSC1c1cc(O)ccc1O.C1CSCN1>>` | `COC(=O)C1CSC(c2cc(O)ccc2O)N1C(C)=O` | `CC(=O)N1C(C(=O)N2CCSC2)CSC1c1cc(O)ccc1O` | relevance model (learned) |
| US05103036 | `Cc1cccc[c]1[Mg][Br].Brc1ccc(I)cc1>c1ccc([P](c2ccccc2)(c2ccccc2)[Pd]([P](c2ccccc2)(c2ccccc2)c2ccccc2)([P](c2ccccc2)(c2ccccc2)c2ccccc2)[P](c2ccccc2)(c2ccccc2)c2ccccc2)cc1> | cond=tetrahydrofuran` | `Cc1ccccc1-c1ccc(Br)cc1` | `Cc1ccccc1-c1ccc(I)cc1.[Br-].[Mg+2]` | relevance model (learned) |
| US20150166505A1 | `COc1cc([N+](=O)[O-])ncc1Br.OC1CC1>[H-].[Na+]> | cond=diethyl ether, NMP` | `COc1cc(OC2CC2)ncc1Br` | `COc1cc([N+](=O)[O-])ncc1OC1CC1.Br` | relevance model (learned) |
| US05977131 | `O=Cc1c[nH]c2ncccc12.COC(=O)CCC[N+](=O)[O-]>CC(=O)[O-].[NH4+]> | cond=tetrahydrofuran` | `COC(=O)CCC(=Cc1c[nH]c2ncccc12)[N+](=O)[O-]` | `O=C1CC[C@@H]([N+](=O)[O-])[C@H](c2c[nH]c3ncccc23)N1` | relevance model (learned) |
| US04061746 | `N[C@@H]1CCC[C@@H]1c1ccccc1.N[C@@H]1CCC[C@H]1c1ccccc1>Cl.Cl>` | `c1ccc([C@@H]2CCC[C@H]2N=C2CCCCCN2)cc1` | `N[C@@H]1CCC[C@@H]1C1CCCCC1` | relevance model (learned) |
| US08143252B2 | `O=C1OC(=O)c2c(F)c(F)c(F)c(F)c21.CN1C(=O)NC(=O)C(C)(N)C1=O>CCN(CC)CC.Cl> | cond=DMSO, water` | `CN1C(=O)NC(=O)C(C)(NC(=O)c2cc(F)c(F)c(F)c2F)C1=O` | `CN1C(=O)NC(=O)C(C)(N2C(=O)c3c(F)c(F)c(F)c(F)c3C2=O)C1=O` | relevance model (learned) |
| US07850742B2 | `CCCn1nc(Br)c([N+](=O)[O-])c1NCc1ccccc1.Oc1ccccc1>O=C([O-])[O-].[Cu][I]> | cond=NMP` | `Cn1nc(Oc2ccccc2)c([N+](=O)[O-])c1NCc1ccccc1` | `CCCn1nc(Oc2ccccc2)c([N+](=O)[O-])c1NCc1ccccc1.Br` | relevance model (learned) |
| US07968538B2 | `COc1ccc(Cl)cc1C(C)=O.NC=O>O=CO> | cond=water` | `COc1ccc(Cl)cc1C(C)N` | `COc1ccc(Cl)cc1C(C)NC=O` | relevance model (learned) |
| US04239702 | `O=C([O-])OCCC(=O)c1ccc(Cl)cc1>[Pd]> | cond=acetic acid` | `OCc1ccc(Cl)cc1` | `O=C([O-])OCCCc1ccc(Cl)cc1` | relevance model (learned) |
| US04704387 | `O=c1c2ccccc2c(Cc2ccc(Cl)cc2)nn1C1CCCNCC1.COCCCl>Br.O=C([O-])[O-].[K+]> | cond=DMA` | `COC(C)N1CCCC(n2nc(Cc3ccc(Cl)cc3)c3ccccc3c2=O)CC1` | `COCCN1CCCC(n2nc(Cc3ccc(Cl)cc3)c3ccccc3c2=O)CC1.Cl` | relevance model (learned) |
| US20140349997A1 | `CC(C)N(C[C@@H](C(=O)O)c1ccc(Cl)cc1)C(=O)OC(C)(C)C.C[C@@H]1C[C@@H](O)c2ncnc(N3CCNCC3)c21>CCN(C(C)C)C(C)C.CN(C)C(On1nnc2ccccc21)=[N+](C)C.Cl.F[P-](F)(F)(F)(F)F> | cond=dichloromethane` | `CC(C)NC(=O)OCCC=O` | `CC(C)N(C[C@@H](C(=O)N1CCN(c2ncnc3c2[C@H](C)C[C@H]3O)CC1)c1ccc(Cl)cc1)C(=O)OC(C)(C)C` | relevance model (learned) |
| US07560455B2 | `CC(C)(C)OC(=O)N[C@@H](CC(=O)O)Cc1cc(F)c(F)cc1F.CC1NCCCn2c1nnc2C(F)(F)F.ClCCCl>CCN(CC)CC.Cl.O=C([O-])O.On1nnc2ccccc21.[Na+]> | cond=DMF` | `CC1c2nnc(C(F)(F)F)n2CCCN1C(=O)C[C@H](N)Cc1cc(F)c(F)cc1F` | `CC1c2nnc(C(F)(F)F)n2CCCN1C(=O)C[C@@H](Cc1cc(F)c(F)cc1F)NC(=O)OC(C)(C)C` | relevance model (learned) |
| US07956216B2 | `CCCNC(=O)c1ccc(Br)cc1.C#Cc1ccccc1>C1CCNCC1.[Cl][Pd]([Cl])([P](c1ccccc1)(c1ccccc1)c1ccccc1)[P](c1ccccc1)(c1ccccc1)c1ccccc1> | cond=dichloromethane, water` | `CCCNC(=O)c1cccc(C#Cc2ccccc2)c1` | `CCCNC(=O)c1ccc(C#Cc2ccccc2)cc1.Br` | relevance model (learned) |
| US05198556 | `FC(F)(F)C1(C(F)(F)F)OC(F)(CCI)C(F)(CCI)O1>[K+].[OH-]>` | `C=CC1(F)OC(C(F)(F)F)(C(F)(F)F)OC1(F)C=C` | `OCCC1(F)OC(C(F)(F)F)(C(F)(F)F)OC1(F)CCI.[I-]` | workbook rule |
| US08778970B2 | `COCCOc1cc(CC2CCNCC2)ccc1Br>BrB(Br)Br.Cl> | cond=dichloromethane, methanol` | `OCCOc1cc(CC2CCNCC2)ccc1Br` | `Oc1cc(CC2CCNCC2)ccc1Br.CO` | relevance model (learned) |
