"""Common molecules to add to the molecule library (js/chem.js MOL_LIBRARY_RAW), so the Sketcher / Molecule search box suggests them.

    python tools/add_molecules.py --check         # every SMILES must give its stated molecular formula (RDKit); lists the problems
    python tools/add_molecules.py --lines KEYS.json   # print the library lines (keys: name -> lookup key, computed in the browser by Chem.libKeyHash(Chem.canonicalKey(...)))

Each entry: (name, [alternative names], SMILES, molecular formula). Stereo is written only where it is known; a flat structure
means the configuration was not specified. No CAS numbers are given here (none is better than a wrong one).
Left out on purpose, too big or too uncertain to write by hand: insulin, heparin, cellulose, starch, PEG, nylon, kevlar,
chlorophyll, heme, cobalamin, vancomycin, gentamicin, azithromycin, erythromycin, ivermectin, digoxin, vincristine, paclitaxel,
NAD+, acetyl-CoA, phylloquinone, stevioside, sucralose, strychnine, cisplatin, potassium permanganate, thiamine, maltose,
doxycycline (stereo), charged species such as paraquat and fumarate.
"""
import json, sys

M = [
 # --- amino acids (L) and relatives
 ("arginine", ["l-arginine", "arg"], "N[C@@H](CCCNC(N)=N)C(=O)O", "C6H14N4O2"),
 ("asparagine", ["l-asparagine", "asn"], "N[C@@H](CC(N)=O)C(=O)O", "C4H8N2O3"),
 ("glutamine", ["l-glutamine", "gln"], "N[C@@H](CCC(N)=O)C(=O)O", "C5H10N2O3"),
 ("histidine", ["l-histidine", "his"], "N[C@@H](Cc1cnc[nH]1)C(=O)O", "C6H9N3O2"),
 ("isoleucine", ["l-isoleucine", "ile"], "CC[C@H](C)[C@H](N)C(=O)O", "C6H13NO2"),
 ("lysine", ["l-lysine", "lys"], "NCCCC[C@H](N)C(=O)O", "C6H14N2O2"),
 ("methionine", ["l-methionine", "met"], "CSCC[C@H](N)C(=O)O", "C5H11NO2S"),
 ("proline", ["l-proline", "pro"], "OC(=O)[C@@H]1CCCN1", "C5H9NO2"),
 ("threonine", ["l-threonine", "thr"], "C[C@@H](O)[C@H](N)C(=O)O", "C4H9NO3"),
 ("tryptophan", ["l-tryptophan", "trp"], "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", "C11H12N2O2"),
 ("tyrosine", ["l-tyrosine", "tyr"], "N[C@@H](Cc1ccc(O)cc1)C(=O)O", "C9H11NO3"),
 ("gamma-aminobutyric acid", ["GABA", "4-aminobutanoic acid"], "NCCCC(=O)O", "C4H9NO2"),
 ("glutathione", ["gsh"], "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O", "C10H17N3O6S"),
 ("betaine", ["trimethylglycine"], "C[N+](C)(C)CC(=O)[O-]", "C5H11NO2"),
 ("aspartame", [], "COC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CC(=O)O", "C14H18N2O5"),
 ("lipoic acid", ["alpha-lipoic acid", "thioctic acid"], "OC(=O)CCCC[C@@H]1CCSS1", "C8H14O2S2"),
 # --- nucleosides and nucleotides
 ("adenosine", [], "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O", "C10H13N5O4"),
 ("guanosine", [], "Nc1nc2c(ncn2[C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]1", "C10H13N5O5"),
 ("cytidine", [], "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1", "C9H13N3O5"),
 ("uridine", [], "O=c1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)[nH]1", "C9H12N2O6"),
 ("thymidine", [], "Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O", "C10H14N2O5"),
 ("adenosine monophosphate", ["AMP", "adenosine 5'-monophosphate"], "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O", "C10H14N5O7P"),
 ("adenosine diphosphate", ["ADP", "adenosine 5'-diphosphate"], "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O", "C10H15N5O10P2"),
 ("zidovudine", ["AZT"], "Cc1cn([C@H]2C[C@H](N=[N+]=[N-])[C@@H](CO)O2)c(=O)[nH]c1=O", "C10H13N5O4"),
 ("acyclovir", ["aciclovir"], "Nc1nc2n(COCCO)cnc2c(=O)[nH]1", "C8H11N5O3"),
 ("5-fluorouracil", ["fluorouracil", "5-FU"], "O=C1NC=C(F)C(=O)N1", "C4H3FN2O2"),
 ("capecitabine", [], "CCCCCOC(=O)Nc1nc(=O)n([C@@H]2O[C@H](C)[C@@H](O)[C@H]2O)cc1F", "C15H22FN3O6"),
 ("methotrexate", [], "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1", "C20H22N8O5"),
 # --- analgesics, anti-inflammatories, anaesthetics
 ("diclofenac", [], "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", "C14H11Cl2NO2"),
 ("ketoprofen", [], "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1", "C16H14O3"),
 ("indomethacin", [], "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1", "C19H16ClNO4"),
 ("celecoxib", [], "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1", "C17H14F3N3O2S"),
 ("oxycodone", [], "COc1ccc2C[C@H]3N(C)CC[C@@]45[C@@H](Oc1c24)C(=O)CC[C@@]35O", "C18H21NO4"),
 ("heroin", ["diacetylmorphine", "diamorphine"], "CC(=O)O[C@H]1C=C[C@H]2[C@H]3Cc4ccc(OC(C)=O)c5O[C@@H]1[C@]2(CCN3C)c45", "C21H23NO5"),
 ("fentanyl", [], "CCC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1", "C22H28N2O"),
 ("methadone", [], "CCC(=O)C(CC(C)N(C)C)(c1ccccc1)c1ccccc1", "C21H27NO"),
 ("tramadol", [], "COc1cccc([C@]2(O)CCCC[C@H]2CN(C)C)c1", "C16H25NO2"),
 ("ketamine", [], "CNC1(c2ccccc2Cl)CCCCC1=O", "C13H16ClNO"),
 ("cocaine", [], "COC(=O)[C@H]1[C@@H]2CC[C@H](C[C@@H]1OC(=O)c1ccccc1)N2C", "C17H21NO4"),
 ("procaine", ["novocaine"], "CCN(CC)CCOC(=O)c1ccc(N)cc1", "C13H20N2O2"),
 ("bupivacaine", [], "CCCCN1CCCCC1C(=O)Nc1c(C)cccc1C", "C18H28N2O"),
 # --- psychiatric and neurological drugs
 ("diazepam", ["valium"], "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", "C16H13ClN2O"),
 ("lorazepam", [], "OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O", "C15H10Cl2N2O2"),
 ("alprazolam", ["xanax"], "Cc1nnc2CN=C(c3ccccc3)c3cc(Cl)ccc3-n12", "C17H13ClN4"),
 ("clonazepam", [], "O=C1CN=C(c2ccccc2Cl)c2cc([N+](=O)[O-])ccc2N1", "C15H10ClN3O3"),
 ("fluoxetine", ["prozac"], "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", "C17H18F3NO"),
 ("sertraline", ["zoloft"], "CN[C@H]1CC[C@@H](c2ccc(Cl)c(Cl)c2)c2ccccc12", "C17H17Cl2N"),
 ("paroxetine", [], "Fc1ccc([C@@H]2CCNC[C@H]2COc2ccc3OCOc3c2)cc1", "C19H20FNO3"),
 ("citalopram", [], "CN(C)CCCC1(c2ccc(F)cc2)OCc2cc(C#N)ccc21", "C20H21FN2O"),
 ("escitalopram", [], "CN(C)CCC[C@@]1(c2ccc(F)cc2)OCc2cc(C#N)ccc21", "C20H21FN2O"),
 ("venlafaxine", [], "COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1", "C17H27NO2"),
 ("amitriptyline", [], "CN(C)CCC=C1c2ccccc2CCc2ccccc12", "C20H23N"),
 ("imipramine", [], "CN(C)CCCN1c2ccccc2CCc2ccccc21", "C19H24N2"),
 ("haloperidol", [], "O=C(CCCN1CCC(O)(c2ccc(Cl)cc2)CC1)c1ccc(F)cc1", "C21H23ClFNO2"),
 ("chlorpromazine", [], "CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc21", "C17H19ClN2S"),
 ("risperidone", [], "Cc1nc2n(c(=O)c1CCN1CCC(c3noc4cc(F)ccc34)CC1)CCCC2", "C23H27FN4O2"),
 ("olanzapine", [], "Cc1cc2c(s1)Nc1ccccc1N=C2N1CCN(C)CC1", "C17H20N4S"),
 ("quetiapine", [], "OCCOCCN1CCN(C2=Nc3ccccc3Sc3ccccc32)CC1", "C21H25N3O2S"),
 ("aripiprazole", [], "O=C1CCc2ccc(OCCCCN3CCN(c4cccc(Cl)c4Cl)CC3)cc2N1", "C23H27Cl2N3O2"),
 ("carbamazepine", [], "NC(=O)N1c2ccccc2C=Cc2ccccc21", "C15H12N2O"),
 ("phenytoin", [], "O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1", "C15H12N2O2"),
 ("phenobarbital", [], "CCC1(c2ccccc2)C(=O)NC(=O)NC1=O", "C12H12N2O3"),
 ("valproic acid", ["valproate", "2-propylpentanoic acid"], "CCCC(CCC)C(=O)O", "C8H16O2"),
 ("gabapentin", [], "NCC1(CC(=O)O)CCCCC1", "C9H17NO2"),
 ("pregabalin", [], "CC(C)C[C@H](CN)CC(=O)O", "C8H17NO2"),
 ("levetiracetam", [], "CC[C@H](C(N)=O)N1CCCC1=O", "C8H14N2O2"),
 ("lamotrigine", [], "Nc1nnc(-c2cccc(Cl)c2Cl)c(N)n1", "C9H7Cl2N5"),
 ("lithium carbonate", [], "O=C([O-])[O-].[Li+].[Li+]", "CLi2O3"),
 ("melatonin", [], "COc1ccc2[nH]cc(CCNC(C)=O)c2c1", "C13H16N2O2"),
 ("norepinephrine", ["noradrenaline"], "NC[C@H](O)c1ccc(O)c(O)c1", "C8H11NO3"),
 ("psilocybin", [], "CN(C)CCc1c[nH]c2cccc(OP(=O)(O)O)c12", "C12H17N2O4P"),
 ("mescaline", [], "COc1cc(CCN)cc(OC)c1OC", "C11H17NO3"),
 ("lsd", ["lysergic acid diethylamide"], "CCN(CC)C(=O)[C@H]1C=C2c3cccc4[nH]cc(c34)C[C@@H]2N(C)C1", "C20H25N3O"),
 ("mdma", ["ecstasy", "3,4-methylenedioxymethamphetamine"], "CC(NC)Cc1ccc2OCOc2c1", "C11H15NO2"),
 ("methamphetamine", [], "C[C@H](Cc1ccccc1)NC", "C10H15N"),
 ("ephedrine", [], "C[C@@H](NC)[C@H](O)c1ccccc1", "C10H15NO"),
 ("pseudoephedrine", [], "C[C@H](NC)[C@H](O)c1ccccc1", "C10H15NO"),
 ("atropine", [], "CN1C2CCC1CC(C2)OC(=O)C(CO)c1ccccc1", "C17H23NO3"),
 ("scopolamine", ["hyoscine"], "CN1C2CC(CC1C1OC21)OC(=O)C(CO)c1ccccc1", "C17H21NO4"),
 ("delta-9-tetrahydrocannabinol", ["THC", "dronabinol"], "CCCCCc1cc(O)c2c(c1)OC(C)(C)[C@@H]1CCC(C)=C[C@H]21", "C21H30O2"),
 ("cannabidiol", ["CBD"], "C=C(C)[C@@H]1CC=C(C)C[C@H]1c1c(O)cc(CCCCC)cc1O", "C21H30O2"),
 # --- antibiotics, antivirals, antifungals, antimalarials
 ("amoxicillin", [], "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccc(O)cc3)C(=O)N2[C@H]1C(=O)O", "C16H19N3O5S"),
 ("ampicillin", [], "CC1(C)S[C@@H]2[C@H](NC(=O)[C@H](N)c3ccccc3)C(=O)N2[C@H]1C(=O)O", "C16H19N3O4S"),
 ("penicillin g", ["benzylpenicillin"], "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O", "C16H18N2O4S"),
 ("penicillin v", ["phenoxymethylpenicillin"], "CC1(C)S[C@@H]2[C@H](NC(=O)COc3ccccc3)C(=O)N2[C@H]1C(=O)O", "C16H18N2O5S"),
 ("cephalexin", ["cefalexin"], "CC1=C(C(=O)O)N2C(=O)[C@@H](NC(=O)[C@H](N)c3ccccc3)[C@H]2SC1", "C16H17N3O4S"),
 ("ciprofloxacin", [], "OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O", "C17H18FN3O3"),
 ("levofloxacin", [], "C[C@H]1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(=O)O)cn1c23", "C18H20FN3O4"),
 ("tetracycline", [], "CC1(O)c2cccc(O)c2C(=O)C2=C(O)C3(O)C(=O)C(C(N)=O)=C(O)C(N(C)C)C3CC21", "C22H24N2O8"),
 ("metronidazole", [], "Cc1ncc([N+](=O)[O-])n1CCO", "C6H9N3O3"),
 ("sulfamethoxazole", [], "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1", "C10H11N3O3S"),
 ("trimethoprim", [], "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC", "C14H18N4O3"),
 ("isoniazid", [], "NNC(=O)c1ccncc1", "C6H7N3O"),
 ("fluconazole", [], "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F", "C13H12F2N6O"),
 ("oseltamivir", ["tamiflu"], "CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H](N)C1", "C16H28N2O4"),
 ("chloroquine", [], "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12", "C18H26ClN3"),
 ("hydroxychloroquine", [], "CCN(CCO)CCCC(C)Nc1ccnc2cc(Cl)ccc12", "C18H26ClN3O"),
 ("quinine", [], "C=C[C@H]1CN2CC[C@H]1C[C@H]2[C@H](O)c1ccnc2ccc(OC)cc12", "C20H24N2O2"),
 ("artemisinin", [], "C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]32OO4", "C15H22O5"),
 # --- cardiovascular, metabolic, GI, respiratory, allergy
 ("metformin", [], "CN(C)C(=N)NC(N)=N", "C4H11N5"),
 ("glipizide", [], "Cc1cnc(C(=O)NCCc2ccc(S(=O)(=O)NC(=O)NC3CCCCC3)cc2)cn1", "C21H27N5O4S"),
 ("atorvastatin", [], "CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O", "C33H35FN2O5"),
 ("simvastatin", [], "CCC(C)(C)C(=O)O[C@H]1C[C@@H](C)C=C2C=C[C@H](C)[C@H](CC[C@@H]3C[C@@H](O)CC(=O)O3)[C@@H]12", "C25H38O5"),
 ("rosuvastatin", [], "CC(C)c1nc(N(C)S(C)(=O)=O)nc(-c2ccc(F)cc2)c1/C=C/[C@@H](O)C[C@@H](O)CC(=O)O", "C22H28FN3O6S"),
 ("lisinopril", [], "NCCCC[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1CCC[C@H]1C(=O)O", "C21H31N3O5"),
 ("enalapril", [], "CCOC(=O)[C@H](CCc1ccccc1)N[C@@H](C)C(=O)N1CCC[C@H]1C(=O)O", "C20H28N2O5"),
 ("captopril", [], "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O", "C9H15NO3S"),
 ("losartan", [], "CCCCc1nc(Cl)c(CO)n1Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1", "C22H23ClN6O"),
 ("valsartan", [], "CCCCC(=O)N(Cc1ccc(-c2ccccc2-c2nn[nH]n2)cc1)[C@H](C(=O)O)C(C)C", "C24H29N5O3"),
 ("amlodipine", [], "CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl", "C20H25ClN2O5"),
 ("nifedipine", [], "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]", "C17H18N2O6"),
 ("verapamil", [], "COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC", "C27H38N2O4"),
 ("diltiazem", [], "COc1ccc([C@@H]2Sc3ccccc3N(CCN(C)C)C(=O)[C@@H]2OC(C)=O)cc1", "C22H26N2O4S"),
 ("metoprolol", [], "COCCc1ccc(OCC(O)CNC(C)C)cc1", "C15H25NO3"),
 ("atenolol", [], "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1", "C14H22N2O3"),
 ("propranolol", [], "CC(C)NCC(O)COc1cccc2ccccc12", "C16H21NO2"),
 ("carvedilol", [], "COc1ccccc1OCCNCC(O)COc1cccc2[nH]c3ccccc3c12", "C24H26N2O4"),
 ("furosemide", [], "NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl", "C12H11ClN2O5S"),
 ("hydrochlorothiazide", [], "NS(=O)(=O)c1cc2c(cc1Cl)NCNS2(=O)=O", "C7H8ClN3O4S2"),
 ("spironolactone", [], "CC(=O)S[C@@H]1CC2=CC(=O)CC[C@]2(C)[C@H]2CC[C@@]3(C)[C@@H](CC[C@@]34CCC(=O)O4)[C@@H]12", "C24H32O4S"),
 ("warfarin", [], "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O", "C19H16O4"),
 ("dicoumarol", [], "O=C1Oc2ccccc2C(O)=C1Cc1c(O)c2ccccc2oc1=O", "C19H12O6"),
 ("clopidogrel", [], "COC(=O)[C@H](c1ccccc1Cl)N1CCc2sccc2C1", "C16H16ClNO2S"),
 ("nitroglycerin", ["glyceryl trinitrate"], "[O-][N+](=O)OCC(CO[N+](=O)[O-])O[N+](=O)[O-]", "C3H5N3O9"),
 ("sildenafil", ["viagra"], "CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12", "C22H30N6O4S"),
 ("tadalafil", [], "CN1CC(=O)N2[C@H](Cc3c([nH]c4ccccc34)[C@@H]2c2ccc3OCOc3c2)C1=O", "C22H19N3O4"),
 ("omeprazole", [], "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1", "C17H19N3O3S"),
 ("esomeprazole", [], "COc1ccc2[nH]c([S@](=O)Cc3ncc(C)c(OC)c3C)nc2c1", "C17H19N3O3S"),
 ("ranitidine", [], "CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1", "C13H22N4O3S"),
 ("famotidine", [], "NC(N)=NC1=NC(CSCCC(N)=NS(N)(=O)=O)=CS1", "C8H15N7O2S3"),
 ("loratadine", [], "CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc32)CC1", "C22H23ClN2O2"),
 ("cetirizine", [], "OC(=O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1", "C21H25ClN2O3"),
 ("fexofenadine", [], "CC(C)(C(=O)O)c1ccc(C(O)CCCN2CCC(C(O)(c3ccccc3)c3ccccc3)CC2)cc1", "C32H39NO4"),
 ("diphenhydramine", ["benadryl"], "CN(C)CCOC(c1ccccc1)c1ccccc1", "C17H21NO"),
 ("montelukast", [], "CC(C)(O)c1ccccc1CC[C@H](SCC1(CC(=O)O)CC1)c1cccc(/C=C/c2ccc3ccc(Cl)cc3n2)c1", "C35H36ClNO3S"),
 ("salbutamol", ["albuterol"], "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1", "C13H21NO3"),
 ("theophylline", [], "Cn1c(=O)c2[nH]cnc2n(C)c1=O", "C7H8N4O2"),
 # --- steroids and hormones
 ("testosterone", [], "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2O", "C19H28O2"),
 ("estradiol", ["oestradiol"], "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CC[C@@H]2O", "C18H24O2"),
 ("estrone", [], "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CCC2=O", "C18H22O2"),
 ("estriol", [], "C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1C[C@H](O)[C@@H]2O", "C18H24O3"),
 ("progesterone", [], "CC(=O)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C", "C21H30O2"),
 ("hydrocortisone", ["cortisol"], "C[C@]12CCC(=O)C=C1CC[C@@H]1[C@@H]2[C@@H](O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO", "C21H30O5"),
 ("cortisone", [], "C[C@]12CCC(=O)C=C1CC[C@@H]1[C@@H]2C(=O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO", "C21H28O5"),
 ("prednisolone", [], "C[C@]12C=CC(=O)C=C1CC[C@@H]1[C@@H]2[C@@H](O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO", "C21H28O5"),
 ("prednisone", [], "C[C@]12C=CC(=O)C=C1CC[C@@H]1[C@@H]2C(=O)C[C@@]2(C)[C@H]1CC[C@]2(O)C(=O)CO", "C21H26O5"),
 ("dexamethasone", [], "C[C@@H]1C[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]2(C)[C@@]1(O)C(=O)CO", "C22H29FO5"),
 ("levothyroxine", ["thyroxine", "T4"], "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(=O)O", "C15H11I4NO4"),
 # --- cancer drugs
 ("cyclophosphamide", [], "ClCCN(CCCl)P1(=O)NCCCO1", "C7H15Cl2N2O2P"),
 ("tamoxifen", [], "CC/C(=C(\\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1", "C26H29NO"),
 ("imatinib", ["gleevec"], "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1", "C29H31N7O"),
 ("gefitinib", [], "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1", "C22H24ClFN4O3"),
 ("erlotinib", [], "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC", "C22H23N3O4"),
 ("doxorubicin", ["adriamycin"], "COc1cccc2c1C(=O)c1c(O)c3c(c(O)c1C2=O)C[C@@](O)(C(=O)CO)C[C@@H]3O[C@H]1C[C@H](N)[C@H](O)[C@H](C)O1", "C27H29NO11"),
 ("camptothecin", [], "CC[C@@]1(O)C(=O)OCc2c1cc1-c3nc4ccccc4cc3Cn1c2=O", "C20H16N2O4"),
 ("podophyllotoxin", [], "COc1cc([C@H]2c3cc4OCOc4cc3[C@H](O)[C@@H]3COC(=O)[C@H]23)cc(OC)c1OC", "C22H22O8"),
 # --- vitamins and lipids
 ("niacin", [], "OC(=O)c1cccnc1", "C6H5NO2"),
 ("riboflavin", ["vitamin b2"], "Cc1cc2nc3c(=O)[nH]c(=O)nc-3n(C[C@H](O)[C@H](O)[C@H](O)CO)c2cc1C", "C17H20N4O6"),
 ("retinol", ["vitamin a"], "C/C(=C\\CO)/C=C/C=C(\\C)/C=C/C1=C(C)CCCC1(C)C", "C20H30O"),
 ("cholecalciferol", ["vitamin d3"], "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2/C(=C/C=C3/C[C@@H](O)CCC3=C)CCC[C@]12C", "C27H44O"),
 ("alpha-tocopherol", ["tocopherol", "vitamin e"], "Cc1c(C)c2c(c(C)c1O)CC[C@@](C)(CCC[C@H](C)CCC[C@H](C)CCCC(C)C)O2", "C29H50O2"),
 ("menadione", ["vitamin k3"], "CC1=CC(=O)c2ccccc2C1=O", "C11H8O2"),
 ("coenzyme q10", ["ubiquinone", "ubidecarenone"], "COC1=C(OC)C(=O)C(C/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CC/C=C(\\C)CCC=C(C)C)=C(C)C1=O", "C59H90O4"),
 ("beta-carotene", [], "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C2=C(C)CCCC2(C)C)C(C)(C)CCC1", "C40H56"),
 ("lycopene", [], "CC(C)=CCC/C(C)=C/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C=C(C)/CCC=C(C)C", "C40H56"),
 ("palmitic acid", ["hexadecanoic acid"], "CCCCCCCCCCCCCCCC(=O)O", "C16H32O2"),
 ("stearic acid", ["octadecanoic acid"], "CCCCCCCCCCCCCCCCCC(=O)O", "C18H36O2"),
 ("oleic acid", [], "CCCCCCCC/C=C\\CCCCCCCC(=O)O", "C18H34O2"),
 ("linoleic acid", [], "CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O", "C18H32O2"),
 ("arachidonic acid", [], "CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCC(=O)O", "C20H32O2"),
 ("sebacic acid", ["decanedioic acid"], "OC(=O)CCCCCCCCC(=O)O", "C10H18O4"),
 # --- natural products and flavours
 ("capsaicin", [], "COc1cc(CNC(=O)CCCC/C=C/C(C)C)ccc1O", "C18H27NO3"),
 ("curcumin", [], "COc1cc(/C=C/C(=O)CC(=O)/C=C/c2ccc(O)c(OC)c2)ccc1O", "C21H20O6"),
 ("resveratrol", [], "Oc1ccc(/C=C/c2cc(O)cc(O)c2)cc1", "C14H12O3"),
 ("quercetin", [], "O=c1c(O)c(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12", "C15H10O7"),
 ("catechin", [], "Oc1cc(O)c2c(c1)O[C@H](c1ccc(O)c(O)c1)[C@@H](O)C2", "C15H14O6"),
 ("epigallocatechin gallate", ["EGCG"], "O=C(O[C@H]1Cc2c(O)cc(O)cc2O[C@@H]1c1cc(O)c(O)c(O)c1)c1cc(O)c(O)c(O)c1", "C22H18O11"),
 ("piperine", [], "O=C(/C=C/C=C/c1ccc2OCOc2c1)N1CCCCC1", "C17H19NO3"),
 ("sulforaphane", [], "CS(=O)CCCCN=C=S", "C6H11NOS2"),
 ("salicin", [], "OC[C@H]1O[C@@H](Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O", "C13H18O7"),
 ("isatin", [], "O=C1Nc2ccccc2C1=O", "C8H5NO2"),
 # --- pesticides and explosives
 ("glyphosate", [], "OC(=O)CNCP(=O)(O)O", "C3H8NO5P"),
 ("atrazine", [], "CCNc1nc(Cl)nc(NC(C)C)n1", "C8H14ClN5"),
 ("malathion", [], "CCOC(=O)CC(SP(=S)(OC)OC)C(=O)OCC", "C10H19O6PS2"),
 ("parathion", [], "CCOP(=S)(OCC)Oc1ccc([N+](=O)[O-])cc1", "C10H14NO5PS"),
 ("permethrin", [], "CC1(C)C(C=C(Cl)Cl)C1C(=O)OCc1cccc(Oc2ccccc2)c1", "C21H20Cl2O3"),
 ("ddt", ["dichlorodiphenyltrichloroethane"], "ClC(Cl)(Cl)C(c1ccc(Cl)cc1)c1ccc(Cl)cc1", "C14H9Cl5"),
 ("tnt", ["trinitrotoluene", "2,4,6-trinitrotoluene"], "Cc1c([N+](=O)[O-])cc([N+](=O)[O-])cc1[N+](=O)[O-]", "C7H5N3O6"),
 ("rdx", ["cyclonite", "hexogen"], "[O-][N+](=O)N1CN([N+](=O)[O-])CN([N+](=O)[O-])C1", "C3H6N6O6"),
 ("petn", ["pentaerythritol tetranitrate"], "[O-][N+](=O)OCC(CO[N+](=O)[O-])(CO[N+](=O)[O-])CO[N+](=O)[O-]", "C5H8N4O12"),
 # --- solvents, reagents, simple compounds
 ("n,n-dimethylformamide", ["DMF", "dimethylformamide"], "CN(C)C=O", "C3H7NO"),
 ("dimethyl sulfoxide", ["DMSO"], "CS(C)=O", "C2H6OS"),
 ("1,4-dioxane", ["dioxane"], "C1COCCO1", "C4H8O2"),
 ("propylene glycol", ["propane-1,2-diol", "1,2-propanediol"], "CC(O)CO", "C3H8O2"),
 ("o-xylene", ["1,2-dimethylbenzene", "ortho-xylene"], "Cc1ccccc1C", "C8H10"),
 ("m-xylene", ["1,3-dimethylbenzene", "meta-xylene"], "Cc1cccc(C)c1", "C8H10"),
 ("p-xylene", ["1,4-dimethylbenzene", "para-xylene"], "Cc1ccc(C)cc1", "C8H10"),
 ("pyrene", [], "c1cc2ccc3cccc4ccc(c1)c2c34", "C16H10"),
 ("benzothiazole", [], "c1ccc2scnc2c1", "C7H5NS"),
 ("1,2,4-triazole", ["triazole"], "c1nc[nH]n1", "C2H3N3"),
 ("tetrazole", [], "c1nnn[nH]1", "CH2N4"),
 ("melamine", [], "Nc1nc(N)nc(N)n1", "C3H6N6"),
 ("hydrazine", [], "NN", "H4N2"),
 ("hydrogen peroxide", [], "OO", "H2O2"),
 ("phosphoric acid", [], "OP(=O)(O)O", "H3O4P"),
 ("sodium chloride", ["salt", "table salt"], "[Na+].[Cl-]", "ClNa"),
 ("sodium bicarbonate", ["sodium hydrogen carbonate", "baking soda"], "[Na+].OC([O-])=O", "CHNaO3"),
 ("butylated hydroxytoluene", ["BHT"], "Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1", "C15H24O"),
 ("dabco", ["1,4-diazabicyclo[2.2.2]octane"], "C1CN2CCN1CC2", "C6H12N2"),
 ("dbu", ["1,8-diazabicyclo[5.4.0]undec-7-ene"], "C1CCC2=NCCCN2CC1", "C9H16N2"),
 ("dmap", ["4-dimethylaminopyridine"], "CN(C)c1ccncc1", "C7H10N2"),
 ("edta", ["ethylenediaminetetraacetic acid"], "OC(=O)CN(CCN(CC(=O)O)CC(=O)O)CC(=O)O", "C10H16N2O8"),
]


def check():
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    bad = 0
    names = set()
    for name, alts, smi, formula in M:
        m = Chem.MolFromSmiles(smi)
        if m is None: print("NOT PARSED :", name, smi); bad += 1; continue
        f = rdMolDescriptors.CalcMolFormula(m)
        if f != formula: print(f"FORMULA    : {name}: {f} (stated {formula})"); bad += 1
        for n in [name] + alts:
            if n.lower() in names: print("DUPLICATE  :", n); bad += 1
            names.add(n.lower())
    print(len(M), "entries,", bad, "problem(s)")
    return bad


TAGS = [("alcohols", "[CX4][OX2H]"), ("acids", "[CX3](=O)[OX2H1,OX1-]"), ("amines", "[NX3;!$(NC=O);!$(N-[S,P]=O);!$(N=*);!$(N-N);!a][#6]"),
        ("aldehydes", "[CX3H1](=O)[#6]"), ("ketones", "[#6][CX3](=O)[#6]"), ("ethers", "[OD2]([#6])[#6]"), ("esters", "[#6][CX3](=O)[OX2][#6]"),
        ("amides", "[CX3](=O)[NX3]"), ("halides", "[F,Cl,Br,I]"), ("alkenes", "[CX3]=[CX3]"), ("alkynes", "[CX2]#[CX2]"), ("nitriles", "[NX1]#[CX2]"),
        ("thiols", "[SX2H]"), ("aromatics", "a")]


def tags_of(smi):
    from rdkit import Chem
    m = Chem.MolFromSmiles(smi)
    out = [t for t, sm in TAGS if m.HasSubstructMatch(Chem.MolFromSmarts(sm))]
    if m.GetRingInfo().NumRings(): out.append("rings")
    if Chem.FindMolChiralCenters(m, includeUnassigned=True, useLegacyImplementation=False): out.append("chiral")
    return ",".join(out)


def lines(keys):
    from rdkit import Chem
    out = []
    for name, alts, smi, formula in M:
        if name not in keys: continue                  # already in the library under another name: only gets extra search names
        m = Chem.MolFromSmiles(smi)
        Chem.Kekulize(m, clearAromaticFlags=True)
        ksmi = Chem.MolToSmiles(m, kekuleSmiles=True)
        out.append("|".join([name, ksmi, keys[name], formula, tags_of(smi), ";".join(alts), ""]))
    return out


if __name__ == "__main__":
    if "--check" in sys.argv: sys.exit(1 if check() else 0)
    if "--smiles" in sys.argv:
        print(json.dumps([[n, s] for n, _, s, _ in M]))
    elif "--lines" in sys.argv:
        keys = json.load(open(sys.argv[sys.argv.index("--lines") + 1]))
        print("\n".join(lines(keys)))
