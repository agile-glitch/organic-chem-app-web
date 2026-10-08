/* More names for the Sketcher's search box.

   1. EXTRA_COMMON: common / trivial names the molecule library lacks, added to Chem.COMMON_NAMES (name → SMILES).
   2. IupacGen: systematic names. Families of simple molecules (chains, branched alkanes, haloalkanes, alcohols,
      ketones, acids and their derivatives, ethers, esters, amines, cycloalkanes, benzenes ...) are enumerated as SMILES
      and named by the app's own IUPAC namer (Chem.iupacName), so every suggestion is a name the namer writes and the
      parser reads back. Built lazily in small slices after the first use, so the page never stalls. */
(function () {
  'use strict';
  const C = window.Chem; if (!C) return;

  const EXTRA_COMMON = {
    'acenaphthene': 'C1Cc2cccc3cccc1c23', 'acenaphthylene': 'C1=Cc2cccc3cccc1c23', 'acetal': 'CC(OCC)OCC', 'acetoin': 'CC(O)C(C)=O',
    'acetoacetic acid': 'CC(=O)CC(=O)O', 'acetoxime': 'CC(C)=NO', 'acetyl bromide': 'CC(=O)Br', 'acetyl fluoride': 'CC(=O)F', 'acetyl iodide': 'CC(=O)I',
    'acetylglycine': 'CC(=O)NCC(=O)O', 'acetylcysteine': 'CC(=O)NC(CS)C(=O)O', 'acetazolamide': 'CC(=O)Nc1nnc(s1)S(N)(=O)=O',
    'aceclofenac': 'OC(=O)COC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl', 'acesulfame': 'CC1=CC(=O)NS(=O)(=O)O1', 'acetophenone oxime': 'CC(=NO)c1ccccc1',
    'acetonylacetone': 'CC(=O)CCC(C)=O', 'acetoxyacetic acid': 'CC(=O)OCC(=O)O', 'acetylsalicylic acid': 'CC(=O)Oc1ccccc1C(=O)O',
    'acetanisole': 'COc1ccc(cc1)C(C)=O', 'acetaldoxime': 'CC=NO', 'acetone cyanohydrin': 'CC(C)(O)C#N', 'acetone azine': 'CC(C)=NN=C(C)C',
    'acetylene dicarboxylic acid': 'OC(=O)C#CC(=O)O', 'acetylacetonate': 'CC(=O)[CH-]C(C)=O', 'acetyl phosphate': 'CC(=O)OP(=O)(O)O',
    'acetonitrile oxide': 'CC#[N+][O-]', 'aceturic acid': 'CC(=O)NCC(=O)O', 'acetamidine': 'CC(N)=N', 'acetohydrazide': 'CC(=O)NN',
    'acetyl isocyanate': 'CC(=O)N=C=O', 'acetylpyridine': 'CC(=O)c1ccccn1', 'acetylthiophene': 'CC(=O)c1cccs1', 'acetylferrocene': 'CC(=O)[c-]1cccc1',
    'acetylpyrrole': 'CC(=O)c1ccc[nH]1', 'acetylfuran': 'CC(=O)c1ccco1', 'acetylcyclohexane': 'CC(=O)C1CCCCC1', 'acetylacetic ester': 'CCOC(=O)CC(C)=O',
    'acetoacetanilide': 'CC(=O)CC(=O)Nc1ccccc1', 'acetonitrile-d3': '[2H]C([2H])([2H])C#N',
    'dimethylformamide': 'CN(C)C=O', 'diethylformamide': 'CCN(CC)C=O', 'dimethylacetamide': 'CC(=O)N(C)C', 'hexamethylphosphoramide': 'CN(C)P(=O)(N(C)C)N(C)C',
    'tetrahydrofuran': 'C1CCOC1', 'tetrahydrothiophene': 'C1CCSC1', 'diglyme': 'COCCOCCOC', 'triglyme': 'COCCOCCOCCOC', 'glyme': 'COCCOC',
    'diethylene glycol': 'OCCOCCO', 'triethylene glycol': 'OCCOCCOCCO', 'propylene oxide': 'CC1CO1', 'epichlorohydrin': 'ClCC1CO1', 'styrene oxide': 'C1OC1c1ccccc1',
    'cyclohexene oxide': 'C1CCC2OC2C1', 'cyclopentene oxide': 'C1CC2OC2C1', 'isobutylene oxide': 'CC1(C)CO1', 'glycidol': 'OCC1CO1', 'butylene oxide': 'CCC1CO1',
    'methyl ethyl ketone': 'CCC(C)=O', 'methyl isobutyl ketone': 'CC(C)CC(C)=O', 'methyl propyl ketone': 'CCCC(C)=O', 'diacetyl': 'CC(=O)C(C)=O',
    'glyoxal': 'O=CC=O', 'methylglyoxal': 'CC(=O)C=O', 'glutaraldehyde': 'O=CCCCC=O', 'succinaldehyde': 'O=CCCC=O', 'crotonaldehyde': 'CC=CC=O',
    'cinnamaldehyde': 'O=CC=Cc1ccccc1', 'furfural': 'O=Cc1ccco1', 'citral': 'CC(C)=CCCC(C)=CC=O', 'chloral': 'ClC(Cl)(Cl)C=O', 'paraformaldehyde unit': 'C=O',
    'mesityl oxide': 'CC(C)=CC(C)=O', 'diacetone alcohol': 'CC(=O)CC(C)(C)O', 'phorone': 'CC(C)=CC(=O)C=C(C)C', 'isophorone': 'CC1=CC(=O)CC(C)(C)C1',
    'dibenzylideneacetone': 'O=C(C=Cc1ccccc1)C=Cc1ccccc1', 'benzil': 'O=C(C(=O)c1ccccc1)c1ccccc1', 'benzoin': 'O=C(C(O)c1ccccc1)c1ccccc1', 'chalcone': 'O=C(C=Cc1ccccc1)c1ccccc1',
    'propiophenone': 'CCC(=O)c1ccccc1', 'butyrophenone': 'CCCC(=O)c1ccccc1', 'benzylacetone': 'CC(=O)CCc1ccccc1', 'phenylacetone': 'CC(=O)Cc1ccccc1',
    'phenylacetic acid': 'OC(=O)Cc1ccccc1', 'phenylacetaldehyde': 'O=CCc1ccccc1', 'mandelic acid': 'OC(C(=O)O)c1ccccc1', 'cinnamic acid': 'OC(=O)C=Cc1ccccc1',
    'hydrocinnamic acid': 'OC(=O)CCc1ccccc1', 'phenylpropanoic acid': 'OC(=O)CCc1ccccc1', 'tropic acid': 'OCC(C(=O)O)c1ccccc1', 'atropic acid': 'C=C(C(=O)O)c1ccccc1',
    'o-toluic acid': 'Cc1ccccc1C(=O)O', 'm-toluic acid': 'Cc1cccc(c1)C(=O)O', 'p-toluic acid': 'Cc1ccc(cc1)C(=O)O', 'anthranilic acid': 'Nc1ccccc1C(=O)O',
    'p-aminobenzoic acid': 'Nc1ccc(cc1)C(=O)O', 'gallic acid': 'OC(=O)c1cc(O)c(O)c(O)c1', 'vanillic acid': 'COc1cc(ccc1O)C(=O)O', 'veratric acid': 'COc1ccc(cc1OC)C(=O)O',
    'protocatechuic acid': 'OC(=O)c1ccc(O)c(O)c1', 'p-hydroxybenzoic acid': 'OC(=O)c1ccc(O)cc1', 'caffeic acid': 'OC(=O)C=Cc1ccc(O)c(O)c1', 'ferulic acid': 'COc1cc(C=CC(=O)O)ccc1O',
    'coumaric acid': 'OC(=O)C=Cc1ccc(O)cc1', 'sinapic acid': 'COc1cc(C=CC(=O)O)cc(OC)c1O', 'isophthalic acid': 'OC(=O)c1cccc(c1)C(=O)O', 'trimellitic acid': 'OC(=O)c1ccc(C(=O)O)c(c1)C(=O)O',
    'pyromellitic acid': 'OC(=O)c1cc(C(=O)O)c(cc1C(=O)O)C(=O)O', 'naphthoic acid': 'OC(=O)c1cccc2ccccc12', 'picolinic acid': 'OC(=O)c1ccccn1', 'isonicotinic acid': 'OC(=O)c1ccncc1',
    'quinolinic acid': 'OC(=O)c1cccnc1C(=O)O', 'furoic acid': 'OC(=O)c1ccco1', 'thiophene-2-carboxylic acid': 'OC(=O)c1cccs1', 'pyrrole-2-carboxylic acid': 'OC(=O)c1ccc[nH]1',
    'glycolic acid': 'OCC(=O)O', 'glyoxylic acid': 'O=CC(=O)O', 'tartaric acid': 'OC(C(O)C(=O)O)C(=O)O', 'malic acid': 'OC(CC(=O)O)C(=O)O', 'fumaric acid': 'OC(=O)/C=C/C(=O)O',
    'itaconic acid': 'C=C(CC(=O)O)C(=O)O', 'glutaric acid': 'OC(=O)CCCC(=O)O', 'adipic acid': 'OC(=O)CCCCC(=O)O', 'pimelic acid': 'OC(=O)CCCCCC(=O)O', 'suberic acid': 'OC(=O)CCCCCCC(=O)O',
    'azelaic acid': 'OC(=O)CCCCCCCC(=O)O', 'sebacic acid': 'OC(=O)CCCCCCCCC(=O)O', 'crotonic acid': 'C/C=C/C(=O)O', 'methacrylic acid': 'CC(=C)C(=O)O', 'sorbic acid': 'C/C=C/C=C/C(=O)O',
    'oxaloacetic acid': 'OC(=O)CC(=O)C(=O)O', 'ketoglutaric acid': 'OC(=O)CCC(=O)C(=O)O', 'levulinic acid': 'CC(=O)CCC(=O)O', 'isobutyric acid': 'CC(C)C(=O)O', 'pivalic acid': 'CC(C)(C)C(=O)O',
    'valeric acid': 'CCCCC(=O)O', 'caproic acid': 'CCCCCC(=O)O', 'enanthic acid': 'CCCCCCC(=O)O', 'caprylic acid': 'CCCCCCCC(=O)O', 'pelargonic acid': 'CCCCCCCCC(=O)O',
    'capric acid': 'CCCCCCCCCC(=O)O', 'lauric acid': 'CCCCCCCCCCCC(=O)O', 'myristic acid': 'CCCCCCCCCCCCCC(=O)O', 'palmitic acid': 'CCCCCCCCCCCCCCCC(=O)O', 'stearic acid': 'CCCCCCCCCCCCCCCCCC(=O)O',
    'oleic acid': 'CCCCCCCC/C=C\\CCCCCCCC(=O)O', 'linoleic acid': 'CCCCC/C=C\\C/C=C\\CCCCCCCC(=O)O', 'arachidic acid': 'CCCCCCCCCCCCCCCCCCCC(=O)O', 'chloroacetic acid': 'OC(=O)CCl',
    'dichloroacetic acid': 'OC(=O)C(Cl)Cl', 'trichloroacetic acid': 'OC(=O)C(Cl)(Cl)Cl', 'bromoacetic acid': 'OC(=O)CBr', 'iodoacetic acid': 'OC(=O)CI', 'fluoroacetic acid': 'OC(=O)CF',
    'thioglycolic acid': 'OC(=O)CS', 'thioacetic acid': 'CC(=O)S', 'peroxyacetic acid': 'CC(=O)OO', 'performic acid': 'O=COO', 'carbonic acid': 'OC(O)=O', 'carbamic acid': 'NC(=O)O',
    'carbamide': 'NC(N)=O', 'thiourea': 'NC(N)=S', 'semicarbazide': 'NNC(N)=O', 'biuret': 'NC(=O)NC(N)=O', 'allantoin': 'NC(=O)NC1NC(=O)NC1=O', 'creatine': 'CN(CC(=O)O)C(N)=N',
    'creatinine': 'CN1CC(=O)NC1=N', 'urethane': 'CCOC(N)=O', 'ethyl carbamate': 'CCOC(N)=O', 'methyl carbamate': 'COC(N)=O', 'ethyl formate': 'CCOC=O', 'methyl formate': 'COC=O',
    'propyl acetate': 'CCCOC(C)=O', 'isopropyl acetate': 'CC(C)OC(C)=O', 'butyl acetate': 'CCCCOC(C)=O', 'isobutyl acetate': 'CC(C)COC(C)=O', 'tert-butyl acetate': 'CC(=O)OC(C)(C)C',
    'amyl acetate': 'CCCCCOC(C)=O', 'vinyl acetate': 'C=COC(C)=O', 'allyl acetate': 'C=CCOC(C)=O', 'benzyl acetate': 'CC(=O)OCc1ccccc1', 'phenyl acetate': 'CC(=O)Oc1ccccc1',
    'ethyl propionate': 'CCOC(=O)CC', 'ethyl butyrate': 'CCCC(=O)OCC', 'methyl butyrate': 'CCCC(=O)OC', 'ethyl benzoate': 'CCOC(=O)c1ccccc1', 'methyl salicylate': 'COC(=O)c1ccccc1O',
    'ethyl acrylate': 'C=CC(=O)OCC', 'methyl acrylate': 'C=CC(=O)OC', 'butyl acrylate': 'C=CC(=O)OCCCC', 'ethyl cinnamate': 'CCOC(=O)C=Cc1ccccc1', 'ethyl lactate': 'CCOC(=O)C(C)O',
    'dimethyl oxalate': 'COC(=O)C(=O)OC', 'diethyl oxalate': 'CCOC(=O)C(=O)OCC', 'dimethyl malonate': 'COC(=O)CC(=O)OC', 'dimethyl succinate': 'COC(=O)CCC(=O)OC', 'diethyl succinate': 'CCOC(=O)CCC(=O)OCC',
    'dimethyl phthalate': 'COC(=O)c1ccccc1C(=O)OC', 'diethyl phthalate': 'CCOC(=O)c1ccccc1C(=O)OCC', 'dibutyl phthalate': 'CCCCOC(=O)c1ccccc1C(=O)OCCCC', 'dimethyl terephthalate': 'COC(=O)c1ccc(cc1)C(=O)OC',
    'diethyl carbonate': 'CCOC(=O)OCC', 'ethylene carbonate': 'O=C1OCCO1', 'propylene carbonate': 'CC1COC(=O)O1', 'vinylene carbonate': 'O=C1OC=CO1', 'diphenyl carbonate': 'O=C(Oc1ccccc1)Oc1ccccc1',
    'triacetin': 'CC(=O)OCC(COC(C)=O)OC(C)=O', 'tributyrin': 'CCCC(=O)OCC(COC(=O)CCC)OC(=O)CCC', 'tripalmitin': 'CCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC',
    'acetyl chloride': 'CC(=O)Cl', 'propionyl chloride': 'CCC(=O)Cl', 'butyryl chloride': 'CCCC(=O)Cl', 'isobutyryl chloride': 'CC(C)C(=O)Cl', 'pivaloyl chloride': 'CC(C)(C)C(=O)Cl',
    'oxalyl chloride': 'ClC(=O)C(Cl)=O', 'thionyl chloride': 'ClS(Cl)=O', 'sulfuryl chloride': 'ClS(Cl)(=O)=O', 'phosphorus oxychloride': 'ClP(Cl)(Cl)=O', 'phosphorus trichloride': 'ClP(Cl)Cl',
    'phosphorus pentachloride': 'ClP(Cl)(Cl)(Cl)Cl', 'phosphoric acid': 'OP(O)(O)=O', 'phosphorous acid': 'OP(O)O', 'triphenyl phosphate': 'O=P(Oc1ccccc1)(Oc1ccccc1)Oc1ccccc1',
    'triethyl phosphate': 'CCOP(=O)(OCC)OCC', 'triphenylphosphine oxide': 'O=P(c1ccccc1)(c1ccccc1)c1ccccc1', 'tributylphosphine': 'CCCCP(CCCC)CCCC', 'trimethylphosphine': 'CP(C)C',
    'dimethyl sulfoxide': 'CS(C)=O', 'dimethyl sulfone': 'CS(C)(=O)=O', 'diethyl sulfide': 'CCSCC', 'diphenyl sulfide': 'c1ccc(cc1)Sc1ccccc1', 'diphenyl sulfoxide': 'O=S(c1ccccc1)c1ccccc1',
    'diphenyl sulfone': 'O=S(=O)(c1ccccc1)c1ccccc1', 'thioanisole': 'CSc1ccccc1', 'thiophenol': 'Sc1ccccc1', 'benzyl mercaptan': 'SCc1ccccc1', 'dimethyl disulfide': 'CSSC',
    'diethyl disulfide': 'CCSSCC', 'carbon disulfide': 'S=C=S', 'carbonyl sulfide': 'O=C=S', 'thioacetamide': 'CC(N)=S', 'thiobenzophenone': 'S=C(c1ccccc1)c1ccccc1',
    'lawesson reagent': 'COc1ccc(cc1)P1(=S)SP(=S)(S1)c1ccc(OC)cc1', 'methanesulfonamide': 'CS(N)(=O)=O', 'sulfanilamide': 'Nc1ccc(cc1)S(N)(=O)=O', 'saccharin': 'O=C1NS(=O)(=O)c2ccccc12',
    'taurine': 'NCCS(=O)(=O)O', 'cysteamine': 'NCCS', 'homocysteine': 'NC(CCS)C(=O)O', 'cystine': 'NC(CSSCC(N)C(=O)O)C(=O)O', 'ethionine': 'CCSCCC(N)C(=O)O',
    'methylamine hydrochloride': 'C[NH3+].[Cl-]', 'ammonium chloride': '[NH4+].[Cl-]', 'tetramethylammonium chloride': 'C[N+](C)(C)C.[Cl-]', 'tetrabutylammonium bromide': 'CCCC[N+](CCCC)(CCCC)CCCC.[Br-]',
    'choline': 'C[N+](C)(C)CCO', 'acetylcholine': 'CC(=O)OCC[N+](C)(C)C', 'betaine': 'C[N+](C)(C)CC(=O)[O-]', 'carnitine': 'C[N+](C)(C)CC(O)CC(=O)[O-]', 'muscarine': 'CC1OC(C[N+](C)(C)C)CC1O',
    'propylamine': 'CCCN', 'isopropylamine': 'CC(C)N', 'butylamine': 'CCCCN', 'sec-butylamine': 'CCC(C)N', 'tert-butylamine': 'CC(C)(C)N', 'isobutylamine': 'CC(C)CN', 'pentylamine': 'CCCCCN',
    'hexylamine': 'CCCCCCN', 'cyclohexylamine': 'NC1CCCCC1', 'cyclopentylamine': 'NC1CCCC1', 'benzylamine': 'NCc1ccccc1', 'phenethylamine': 'NCCc1ccccc1', 'amphetamine': 'CC(N)Cc1ccccc1',
    'methamphetamine': 'CNC(C)Cc1ccccc1', 'allylamine': 'C=CCN', 'propargylamine': 'C#CCN', 'ethanolamine': 'NCCO', 'diethanolamine': 'OCCNCCO', 'triethanolamine': 'OCCN(CCO)CCO',
    'dipropylamine': 'CCCNCCC', 'dibutylamine': 'CCCCNCCCC', 'tributylamine': 'CCCCN(CCCC)CCCC', 'tripropylamine': 'CCCN(CCC)CCC', 'diisopropylethylamine': 'CCN(C(C)C)C(C)C',
    'dicyclohexylamine': 'C1CCC(CC1)NC1CCCCC1', 'diphenylamine': 'c1ccc(cc1)Nc1ccccc1', 'triphenylamine': 'c1ccc(cc1)N(c1ccccc1)c1ccccc1', 'n-methylaniline': 'CNc1ccccc1',
    'n,n-dimethylaniline': 'CN(C)c1ccccc1', 'n,n-diethylaniline': 'CCN(CC)c1ccccc1', 'toluidine': 'Cc1ccccc1N', 'o-toluidine': 'Cc1ccccc1N', 'p-toluidine': 'Cc1ccc(N)cc1',
    'm-toluidine': 'Cc1cccc(N)c1', 'o-phenylenediamine': 'Nc1ccccc1N', 'm-phenylenediamine': 'Nc1cccc(N)c1', 'p-phenylenediamine': 'Nc1ccc(N)cc1', 'benzidine': 'Nc1ccc(cc1)-c1ccc(N)cc1',
    '1-naphthylamine': 'Nc1cccc2ccccc12', '2-naphthylamine': 'Nc1ccc2ccccc2c1', '1-naphthol': 'Oc1cccc2ccccc12', '2-naphthol': 'Oc1ccc2ccccc2c1', 'p-anisidine': 'COc1ccc(N)cc1',
    'phenetidine': 'CCOc1ccc(N)cc1', 'phenacetin': 'CCOc1ccc(NC(C)=O)cc1', 'sulfanilic acid': 'Nc1ccc(cc1)S(O)(=O)=O', 'dapsone': 'Nc1ccc(cc1)S(=O)(=O)c1ccc(N)cc1',
    'nitroethane': 'CC[N+](=O)[O-]', 'nitropropane': 'CCC[N+](=O)[O-]', '2-nitropropane': 'CC(C)[N+](=O)[O-]', 'nitroglycerin': '[O-][N+](=O)OCC(CO[N+]([O-])=O)O[N+]([O-])=O',
    'trinitrotoluene': 'Cc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]', 'picric acid': 'Oc1c(cc(cc1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-]', '2,4-dinitrophenol': 'Oc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]',
    '2,4-dinitrophenylhydrazine': 'NNc1ccc(cc1[N+](=O)[O-])[N+](=O)[O-]', 'o-nitrophenol': 'Oc1ccccc1[N+](=O)[O-]', 'm-nitroaniline': 'Nc1cccc(c1)[N+](=O)[O-]', 'o-nitroaniline': 'Nc1ccccc1[N+](=O)[O-]',
    'm-dinitrobenzene': '[O-][N+](=O)c1cccc(c1)[N+](=O)[O-]', 'nitrotoluene': 'Cc1ccccc1[N+](=O)[O-]', 'p-nitrotoluene': 'Cc1ccc(cc1)[N+](=O)[O-]', 'nitronaphthalene': '[O-][N+](=O)c1cccc2ccccc12',
    'azobisisobutyronitrile': 'CC(C)(C#N)N=NC(C)(C)C#N', 'azoxybenzene': '[O-][N+](=Nc1ccccc1)c1ccccc1', 'hydrazobenzene': 'c1ccc(NNc2ccccc2)cc1', 'benzaldehyde oxime': 'ON=Cc1ccccc1',
    'cyclohexanone oxime': 'ON=C1CCCCC1', 'dimethylglyoxime': 'CC(=NO)C(C)=NO', 'semicarbazone': 'NC(=O)NN=C', 'phenylhydrazone': 'C=NNc1ccccc1', 'hydroxyurea': 'NC(=O)NO',
    'urea nitrate': 'NC(N)=O.O[N+]([O-])=O', 'cyanamide': 'NC#N', 'dicyandiamide': 'NC(N)=NC#N', 'melamine': 'Nc1nc(N)nc(N)n1', 'cyanuric chloride': 'Clc1nc(Cl)nc(Cl)n1', 'cyanuric acid': 'O=c1[nH]c(=O)[nH]c(=O)[nH]1',
    'hydrogen cyanide': 'C#N', 'cyanogen': 'N#CC#N', 'malononitrile': 'N#CCC#N', 'succinonitrile': 'N#CCCC#N', 'adiponitrile': 'N#CCCCCC#N', 'propionitrile': 'CCC#N', 'butyronitrile': 'CCCC#N',
    'benzyl cyanide': 'N#CCc1ccccc1', 'isocyanic acid': 'N=C=O', 'methyl isocyanate': 'CN=C=O', 'methyl isothiocyanate': 'CN=C=S', 'phenyl isothiocyanate': 'S=C=Nc1ccccc1', 'toluene diisocyanate': 'Cc1ccc(N=C=O)cc1N=C=O',
    'hexamethylene diisocyanate': 'O=C=NCCCCCCN=C=O', 'methylene diphenyl diisocyanate': 'O=C=Nc1ccc(Cc2ccc(cc2)N=C=O)cc1', 'azidomethane': 'CN=[N+]=[N-]', 'phenyl azide': '[N-]=[N+]=Nc1ccccc1', 'sodium azide': '[N-]=[N+]=[N-].[Na+]',
    'ethyl diazoacetate': 'CCOC(=O)C=[N+]=[N-]', 'trimethylsilyldiazomethane': 'C[Si](C)(C)C=[N+]=[N-]', 'benzenediazonium chloride': '[N+](#Nc1ccccc1).[Cl-]', 'diazald': 'Cc1ccc(cc1)S(=O)(=O)N(C)N=O',
    'chloromethane': 'CCl', 'bromomethane': 'CBr', 'fluoromethane': 'CF', 'dibromomethane': 'BrCBr', 'diiodomethane': 'ICI', 'bromoform': 'BrC(Br)Br', 'iodoform': 'IC(I)I', 'tetrachloroethylene': 'ClC(Cl)=C(Cl)Cl',
    'trichloroethylene': 'ClC=C(Cl)Cl', '1,1-dichloroethylene': 'C=C(Cl)Cl', 'vinyl bromide': 'C=CBr', 'allyl chloride': 'C=CCCl', 'chloroprene': 'C=CC(Cl)=C', 'tetrafluoroethylene': 'FC(F)=C(F)F',
    'hexafluoropropylene': 'FC(F)=C(F)C(F)(F)F', 'difluoromethane': 'FCF', 'freon-12': 'FC(F)(Cl)Cl', 'freon-11': 'FC(Cl)(Cl)Cl', 'halothane': 'FC(F)(F)C(Cl)Br', 'sevoflurane': 'FCOC(C(F)(F)F)C(F)(F)F',
    'isoflurane': 'FC(F)OC(Cl)C(F)(F)F', 'desflurane': 'FC(F)OC(F)C(F)(F)F', 'chlorobutane': 'CCCCCl', '1-chlorobutane': 'CCCCCl', '1-bromobutane': 'CCCCBr', '1-bromopropane': 'CCCBr', '1-iodobutane': 'CCCCI',
    '1-chloropropane': 'CCCCl', '2-chloropropane': 'CC(C)Cl', 'isopropyl bromide': 'CC(C)Br', 'isopropyl iodide': 'CC(C)I', 'sec-butyl bromide': 'CCC(C)Br', 'tert-butyl bromide': 'CC(C)(C)Br', 'tert-butyl iodide': 'CC(C)(C)I',
    'neopentyl bromide': 'CC(C)(C)CBr', 'benzyl chloride': 'ClCc1ccccc1', 'benzal chloride': 'ClC(Cl)c1ccccc1', 'benzotrichloride': 'ClC(Cl)(Cl)c1ccccc1', 'benzyl iodide': 'ICc1ccccc1', 'trityl chloride': 'ClC(c1ccccc1)(c1ccccc1)c1ccccc1',
    'chlorodiphenylmethane': 'ClC(c1ccccc1)c1ccccc1', 'allyl iodide': 'C=CCI', 'propargyl bromide': 'C#CCBr', 'bromoacetone': 'CC(=O)CBr', 'chloroacetone': 'CC(=O)CCl', 'chloroacetyl chloride': 'ClCC(Cl)=O',
    'bromoacetyl bromide': 'BrCC(Br)=O', 'ethyl bromoacetate': 'CCOC(=O)CBr', 'ethyl chloroacetate': 'CCOC(=O)CCl', 'methyl chloroformate': 'COC(Cl)=O', 'benzyl chloroformate': 'ClC(=O)OCc1ccccc1',
    'di-tert-butyl dicarbonate': 'CC(C)(C)OC(=O)OC(=O)OC(C)(C)C', 'boc anhydride': 'CC(C)(C)OC(=O)OC(=O)OC(C)(C)C', 'fmoc chloride': 'ClC(=O)OCC1c2ccccc2-c2ccccc12', 'trifluoroacetic anhydride': 'FC(F)(F)C(=O)OC(=O)C(F)(F)F',
    'trifluoromethanesulfonic acid': 'OS(=O)(=O)C(F)(F)F', 'triflic anhydride': 'FC(F)(F)S(=O)(=O)OS(=O)(=O)C(F)(F)F', 'trifluoroethanol': 'OCC(F)(F)F', 'hexafluoroisopropanol': 'OC(C(F)(F)F)C(F)(F)F',
    'perfluorooctanoic acid': 'OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F', 'trifluoroacetone': 'CC(=O)C(F)(F)F', 'hexafluoroacetone': 'O=C(C(F)(F)F)C(F)(F)F', 'trifluorotoluene': 'FC(F)(F)c1ccccc1',
    'pentafluorophenol': 'Oc1c(F)c(F)c(F)c(F)c1F', 'hexafluorobenzene': 'Fc1c(F)c(F)c(F)c(F)c1F', 'hexachlorobenzene': 'Clc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl', 'dichlorobenzene': 'Clc1ccccc1Cl', 'p-dichlorobenzene': 'Clc1ccc(Cl)cc1',
    'trichlorobenzene': 'Clc1cccc(Cl)c1Cl', 'chlorotoluene': 'Cc1ccccc1Cl', 'p-chlorotoluene': 'Cc1ccc(Cl)cc1', 'bromotoluene': 'Cc1ccccc1Br', 'p-bromotoluene': 'Cc1ccc(Br)cc1', 'p-bromoaniline': 'Nc1ccc(Br)cc1',
    'p-chloroaniline': 'Nc1ccc(Cl)cc1', 'p-chlorophenol': 'Oc1ccc(Cl)cc1', 'p-bromophenol': 'Oc1ccc(Br)cc1', 'pentachlorophenol': 'Oc1c(Cl)c(Cl)c(Cl)c(Cl)c1Cl', 'triclosan': 'Oc1cc(Cl)ccc1Oc1ccc(Cl)cc1Cl',
    'chlorobenzaldehyde': 'O=Cc1ccccc1Cl', 'p-anisaldehyde': 'COc1ccc(C=O)cc1', 'salicylaldehyde': 'O=Cc1ccccc1O', 'p-tolualdehyde': 'Cc1ccc(C=O)cc1', 'p-nitrobenzaldehyde': 'O=Cc1ccc(cc1)[N+](=O)[O-]',
    'p-hydroxybenzaldehyde': 'O=Cc1ccc(O)cc1', 'piperonal': 'O=Cc1ccc2OCOc2c1', 'veratraldehyde': 'COc1ccc(C=O)cc1OC', 'isovanillin': 'COc1ccc(C=O)cc1O', 'syringaldehyde': 'COc1cc(C=O)cc(OC)c1O',
    'ethyl vanillin': 'CCOc1cc(C=O)ccc1O', 'terephthalaldehyde': 'O=Cc1ccc(C=O)cc1', 'phthalaldehyde': 'O=Cc1ccccc1C=O', 'isophthalaldehyde': 'O=Cc1cccc(C=O)c1', '1-naphthaldehyde': 'O=Cc1cccc2ccccc12',
    'benzocaine': 'CCOC(=O)c1ccc(N)cc1', 'lidocaine': 'CCN(CC)CC(=O)Nc1c(C)cccc1C', 'tetracaine': 'CCCCNc1ccc(cc1)C(=O)OCCN(C)C', 'benzyl benzoate': 'O=C(OCc1ccccc1)c1ccccc1', 'phenyl benzoate': 'O=C(Oc1ccccc1)c1ccccc1',
    'benzoic anhydride': 'O=C(OC(=O)c1ccccc1)c1ccccc1', 'benzoyl peroxide': 'O=C(OOC(=O)c1ccccc1)c1ccccc1', 'dicumyl peroxide': 'CC(C)(OOC(C)(C)c1ccccc1)c1ccccc1', 'tert-butyl hydroperoxide': 'CC(C)(C)OO',
    'cumene hydroperoxide': 'CC(C)(OO)c1ccccc1', 'di-tert-butyl peroxide': 'CC(C)(C)OOC(C)(C)C', 'ozone': 'O=[O+][O-]', 'tert-butyl peroxybenzoate': 'CC(C)(C)OOC(=O)c1ccccc1',
    'diethyl azodicarboxylate': 'CCOC(=O)N=NC(=O)OCC', 'diisopropyl azodicarboxylate': 'CC(C)OC(=O)N=NC(=O)OC(C)C', 'dicyclohexylcarbodiimide': 'C1CCC(CC1)N=C=NC1CCCCC1', 'diisopropylcarbodiimide': 'CC(C)N=C=NC(C)C',
    'edc': 'CCN=C=NCCCN(C)C', 'carbonyldiimidazole': 'O=C(n1ccnc1)n1ccnc1', 'hobt': 'On1nnc2ccccc12', 'nhs': 'ON1C(=O)CCC1=O', 'n-hydroxysuccinimide': 'ON1C(=O)CCC1=O', 'pyridine hydrochloride': 'c1cc[nH+]cc1.[Cl-]',
    'imidazole': 'c1c[nH]cn1', 'n-methylimidazole': 'Cn1ccnc1', 'n-methylmorpholine': 'CN1CCOCC1', 'n-methylpiperidine': 'CN1CCCCC1', 'n-methylpyrrolidine': 'CN1CCCC1', 'n-methylpyrrole': 'Cn1cccc1',
    'quinuclidine': 'C1CN2CCC1CC2', 'tropane': 'CN1C2CCC1CCC2', 'proton sponge': 'CN(C)c1cccc2cccc(N(C)C)c12', 'lutidine': 'Cc1cccc(C)n1', '2,6-lutidine': 'Cc1cccc(C)n1', 'collidine': 'Cc1cc(C)nc(C)c1',
    '2,2-bipyridine': 'c1ccc(nc1)-c1ccccn1', 'bipyridine': 'c1ccc(nc1)-c1ccccn1', '1,10-phenanthroline': 'c1cnc2c(c1)ccc1cccnc12', 'terpyridine': 'c1ccc(nc1)-c1cccc(n1)-c1ccccn1', 'pyridoxine': 'Cc1ncc(CO)c(CO)c1O',
    'thiamine': 'Cc1ncc(C[n+]2csc(CCO)c2C)c(N)n1', 'riboflavin': 'Cc1cc2nc3c(nc(=O)[nH]c3=O)n(CC(O)C(O)C(O)CO)c2cc1C', 'nicotinamide': 'NC(=O)c1cccnc1', 'niacin': 'OC(=O)c1cccnc1',
    'biotin': 'OC(=O)CCCCC1SCC2NC(=O)NC12', 'folic acid': 'Nc1nc2NCC(CNc3ccc(cc3)C(=O)NC(CCC(=O)O)C(=O)O)Nc2c(=O)[nH]1', 'ascorbic acid': 'OCC(O)C1OC(=O)C(O)=C1O', 'vitamin c': 'OCC(O)C1OC(=O)C(O)=C1O',
    'retinol': 'CC1=C(C=CC(C)=CC=CC(C)=CCO)C(C)(C)CCC1', 'beta-carotene': 'CC1=C(C(CCC1)(C)C)C=CC(=CC=CC(=CC=CC=C(C=CC=C(C=CC2=C(CCCC2(C)C)C)C)C)C)C', 'tocopherol': 'CC1=C(C)C2=C(CCC(C)(CCCC(C)CCCC(C)CCCC(C)C)O2)C(C)=C1O',
    'caffeine': 'Cn1cnc2c1c(=O)n(C)c(=O)n2C', 'theobromine': 'Cn1cnc2c1c(=O)[nH]c(=O)n2C', 'theophylline': 'Cn1c2nc[nH]c2c(=O)n(C)c1=O', 'xanthine': 'O=c1[nH]c(=O)c2[nH]cnc2[nH]1', 'uric acid': 'O=c1[nH]c2[nH]c(=O)[nH]c2c(=O)[nH]1',
    'hypoxanthine': 'O=c1[nH]cnc2[nH]cnc12', 'orotic acid': 'OC(=O)c1cc(=O)[nH]c(=O)[nH]1', 'allopurinol': 'O=c1[nH]cnc2[nH]ncc12', 'barbital': 'CCC1(CC)C(=O)NC(=O)NC1=O', 'phenobarbital': 'CCC1(C(=O)NC(=O)NC1=O)c1ccccc1',
    'phenytoin': 'O=C1NC(=O)C(N1)(c1ccccc1)c1ccccc1', 'carbamazepine': 'NC(=O)N1c2ccccc2C=Cc2ccccc12', 'valproic acid': 'CCCC(CCC)C(=O)O', 'gabapentin': 'NCC1(CC(=O)O)CCCCC1', 'levetiracetam': 'CCC(N1CCCC1=O)C(N)=O',
    'diazepam': 'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12', 'lorazepam': 'OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2NC1=O', 'alprazolam': 'Cc1nnc2CN=C(c3ccccc3)c3cc(Cl)ccc3-n12', 'haloperidol': 'OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c1ccc(Cl)cc1',
    'fluoxetine': 'CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1', 'sertraline': 'CNC1CCC(c2ccc(Cl)c(Cl)c2)c2ccccc12', 'citalopram': 'CN(C)CCCC1(OCc2cc(ccc12)C#N)c1ccc(F)cc1', 'amitriptyline': 'CN(C)CCC=C1c2ccccc2CCc2ccccc12',
    'imipramine': 'CN(C)CCCN1c2ccccc2CCc2ccccc12', 'chlorpromazine': 'CN(C)CCCN1c2ccccc2Sc2ccc(Cl)cc12', 'clozapine': 'CN1CCN(CC1)C1=Nc2cc(Cl)ccc2Nc2ccccc12', 'bupropion': 'CC(NC(C)(C)C)C(=O)c1cccc(Cl)c1',
    'morphine': 'CN1CCC23C4Oc5c(O)ccc(CC1C2C=CC4O)c35', 'codeine': 'COc1ccc2CC3N(C)CCC45C(Oc1c24)C(O)C=CC35', 'atropine': 'CN1C2CCC1CC(C2)OC(=O)C(CO)c1ccccc1', 'scopolamine': 'CN1C2CC(OC(=O)C(CO)c3ccccc3)CC1C1OC21',
    'quinine': 'COc1ccc2nccc(C(O)C3CC4CCN3CC4C=C)c2c1', 'capsaicin': 'COc1cc(CNC(=O)CCCCC=CC(C)C)ccc1O', 'menthone': 'CC(C)C1CCC(C)CC1=O', 'carvone': 'CC1=CCC(CC1=O)C(C)=C', 'thymol': 'Cc1ccc(C(C)C)c(O)c1',
    'carvacrol': 'Cc1ccc(C(C)C)cc1O', 'eugenol': 'COc1cc(CC=C)ccc1O', 'geraniol': 'CC(C)=CCCC(C)=CCO', 'linalool': 'CC(C)=CCCC(C)(O)C=C', 'farnesol': 'CC(C)=CCCC(C)=CCCC(C)=CCO', 'pinene': 'CC1=CCC2CC1C2(C)C',
    'alpha-pinene': 'CC1=CCC2CC1C2(C)C', 'camphene': 'CC1(C)C2CCC(C2)C1=C', 'borneol': 'CC1(C)C2CCC1(C)C(O)C2', 'cineole': 'CC1(C)OC2CCC1(C)CC2', 'eucalyptol': 'CC1(C)OC2CCC1(C)CC2', 'squalene': 'CC(C)=CCCC(C)=CCCC(C)=CCCC=C(C)CCC=C(C)CCC=C(C)C',
    'testosterone': 'CC12CCC3C(CCC4=CC(=O)CCC34C)C1CCC2O', 'estradiol': 'CC12CCC3C(CCc4cc(O)ccc34)C1CCC2O', 'progesterone': 'CC(=O)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C', 'cortisol': 'CC12CCC(=O)C=C1CCC1C2C(O)CC2(C)C1CCC2(O)C(=O)CO',
    'cholesterol': 'CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C', 'ergosterol': 'CC(C=CC(C)C(C)C)C1CCC2C3=CC=C4CC(O)CCC4(C)C3CCC12C', 'cholic acid': 'CC(CCC(=O)O)C1CCC2C3C(O)CC4CC(O)CCC4(C)C3CC(O)C12C',
    'aspirin': 'CC(=O)Oc1ccccc1C(=O)O', 'paracetamol': 'CC(=O)Nc1ccc(O)cc1', 'ibuprofen': 'CC(C)Cc1ccc(cc1)C(C)C(=O)O', 'naproxen': 'COc1ccc2cc(ccc2c1)C(C)C(=O)O', 'metformin': 'CN(C)C(=N)NC(N)=N',
    'warfarin': 'CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O', 'penicillin g': 'CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O', 'amoxicillin': 'CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O',
    'ciprofloxacin': 'OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O', 'sulfamethoxazole': 'Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1', 'chloramphenicol': 'OCC(NC(=O)C(Cl)Cl)C(O)c1ccc(cc1)[N+](=O)[O-]', 'isoniazid': 'NNC(=O)c1ccncc1',
    'nitrofurantoin': 'O=C1CN(N=Cc2ccc(o2)[N+](=O)[O-])C(=O)N1', 'metronidazole': 'Cc1ncc([N+](=O)[O-])n1CCO', 'fluconazole': 'OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F', 'omeprazole': 'COc1ccc2[nH]c(nc2c1)S(=O)Cc1ncc(C)c(OC)c1C',
    'ranitidine': 'CNC(NCCSCc1ccc(CN(C)C)o1)=C[N+](=O)[O-]', 'propranolol': 'CC(C)NCC(O)COc1cccc2ccccc12', 'atenolol': 'CC(C)NCC(O)COc1ccc(CC(N)=O)cc1', 'salbutamol': 'CC(C)(C)NCC(O)c1ccc(O)c(CO)c1',
    'albuterol': 'CC(C)(C)NCC(O)c1ccc(O)c(CO)c1', 'amlodipine': 'CCOC(=O)C1=C(COCCN)NC(C)=C(C1c1ccccc1Cl)C(=O)OC', 'nifedipine': 'COC(=O)C1=C(C)NC(C)=C(C1c1ccccc1[N+](=O)[O-])C(=O)OC', 'captopril': 'CC(CS)C(=O)N1CCCC1C(=O)O',
    'furosemide': 'NS(=O)(=O)c1cc(C(=O)O)c(NCc2ccco2)cc1Cl', 'hydrochlorothiazide': 'NS(=O)(=O)c1cc2c(NCNS2(=O)=O)cc1Cl', 'sildenafil': 'CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1', 'diphenhydramine': 'CN(C)CCOC(c1ccccc1)c1ccccc1',
    'loratadine': 'CCOC(=O)N1CCC(=C2c3ccc(Cl)cc3CCc3cccnc23)CC1', 'cetirizine': 'OC(=O)COCCN1CCN(CC1)C(c1ccccc1)c1ccc(Cl)cc1', 'lysergic acid diethylamide': 'CCN(CC)C(=O)C1C=C2c3cccc4[nH]cc(CC2N(C)C1)c34', 'lsd': 'CCN(CC)C(=O)C1C=C2c3cccc4[nH]cc(CC2N(C)C1)c34',
    'mdma': 'CC(NC)Cc1ccc2OCOc2c1', 'thc': 'CCCCCc1cc(O)c2C3C=C(C)CCC3C(C)(C)Oc2c1', 'cannabidiol': 'CCCCCc1cc(O)c(C2C=C(C)CCC2C(C)=C)c(O)c1', 'psilocin': 'CN(C)CCc1c[nH]c2cccc(O)c12',
    'melatonin': 'COc1ccc2[nH]cc(CCNC(C)=O)c2c1', 'tryptamine': 'NCCc1c[nH]c2ccccc12', 'histamine': 'NCCc1c[nH]cn1', 'tyramine': 'NCCc1ccc(O)cc1', 'noradrenaline': 'NCC(O)c1ccc(O)c(O)c1', 'norepinephrine': 'NCC(O)c1ccc(O)c(O)c1',
    'epinephrine': 'CNCC(O)c1ccc(O)c(O)c1', 'levodopa': 'NC(Cc1ccc(O)c(O)c1)C(=O)O', 'l-dopa': 'NC(Cc1ccc(O)c(O)c1)C(=O)O', 'gaba': 'NCCCC(=O)O', 'glycine betaine': 'C[N+](C)(C)CC(=O)[O-]',
    'ethylene': 'C=C', 'propylene': 'CC=C', 'butylene': 'CCC=C', 'isobutylene': 'CC(C)=C', 'cyclohexene': 'C1=CCCCC1', 'cyclopentene': 'C1=CCCC1', 'cyclooctene': 'C1=CCCCCCC1', 'norbornene': 'C1CC2C=CC1C2',
    'norbornadiene': 'C1C2C=CC1C=C2', 'cyclopentadiene': 'C1=CCC=C1', 'dicyclopentadiene': 'C1=CC2C3CC(C=C3)C2C1', '1,5-cyclooctadiene': 'C1CC=CCCC=C1', 'cycloheptatriene': 'C1=CC=CCC=C1', 'styrene': 'C=Cc1ccccc1',
    'alpha-methylstyrene': 'CC(=C)c1ccccc1', 'stilbene': 'C(=Cc1ccccc1)c1ccccc1', 'phenylacetylene': 'C#Cc1ccccc1', 'diphenylacetylene': 'c1ccc(cc1)C#Cc1ccccc1', 'tolane': 'c1ccc(cc1)C#Cc1ccccc1', 'divinylbenzene': 'C=Cc1ccccc1C=C',
    'naphthalene': 'c1ccc2ccccc2c1', 'anthracene': 'c1ccc2cc3ccccc3cc2c1', 'phenanthrene': 'c1ccc2c(c1)ccc1ccccc12', 'tetracene': 'c1ccc2cc3cc4ccccc4cc3cc2c1', 'chrysene': 'c1ccc2c(c1)ccc1c2ccc2ccccc12',
    'triphenylene': 'c1ccc2c(c1)c1ccccc1c1ccccc21', 'perylene': 'c1cc2cccc3c4cccc5cccc(c(c1)c23)c54', 'coronene': 'c1cc2ccc3ccc4ccc5ccc6ccc1c7c2c3c4c5c67', 'fullerene': 'c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c%11c6c6c7c7c2c2c8c8c3c3c9c9c4c4c%10c5c5c1c1c%11c6c6c7c2c2c8c3c3c9c4c4c5c1c6c2c34',
    'bicyclohexyl': 'C1CCC(CC1)C1CCCCC1', 'biphenylene': 'c1ccc2c(c1)c1ccccc12', 'terphenyl': 'c1ccc(cc1)-c1ccc(cc1)-c1ccccc1', 'durene': 'Cc1cc(C)c(C)cc1C', 'hexamethylbenzene': 'Cc1c(C)c(C)c(C)c(C)c1C',
    'pseudocumene': 'Cc1ccc(C)c(C)c1', 'cymene': 'Cc1ccc(C(C)C)cc1', 'p-cymene': 'Cc1ccc(C(C)C)cc1', 'tert-butylbenzene': 'CC(C)(C)c1ccccc1', 'butylbenzene': 'CCCCc1ccccc1', 'propylbenzene': 'CCCc1ccccc1',
    'bibenzyl': 'c1ccc(CCc2ccccc2)cc1', 'dibenzyl ether': 'c1ccc(COCc2ccccc2)cc1', 'benzyl ether': 'c1ccc(COCc2ccccc2)cc1', 'phenetole': 'CCOc1ccccc1', 'veratrole': 'COc1ccccc1OC', 'guaiacol': 'COc1ccccc1O',
    'dimethoxybenzene': 'COc1ccc(OC)cc1', '1,4-dimethoxybenzene': 'COc1ccc(OC)cc1', '1,3-dimethoxybenzene': 'COc1cccc(OC)c1', 'pyrogallol': 'Oc1cccc(O)c1O', 'phloroglucinol': 'Oc1cc(O)cc(O)c1', 'hydroxyquinol': 'Oc1cc(O)c(O)cc1',
    'bisphenol a': 'CC(C)(c1ccc(O)cc1)c1ccc(O)cc1', 'bht': 'Cc1cc(c(O)c(c1)C(C)(C)C)C(C)(C)C', 'butylated hydroxytoluene': 'Cc1cc(c(O)c(c1)C(C)(C)C)C(C)(C)C', 'butylated hydroxyanisole': 'COc1ccc(O)c(c1)C(C)(C)C',
    '4-tert-butylphenol': 'CC(C)(C)c1ccc(O)cc1', 'nonylphenol': 'CCCCCCCCCc1ccc(O)cc1', 'o-cresol': 'Cc1ccccc1O', 'm-cresol': 'Cc1cccc(O)c1', 'xylenol': 'Cc1cccc(C)c1O', '2,6-dimethylphenol': 'Cc1cccc(C)c1O',
    'ethyl alcohol': 'CCO', 'methyl alcohol': 'CO', 'propyl alcohol': 'CCCO', 'butyl alcohol': 'CCCCO', 'isobutyl alcohol': 'CC(C)CO', 'sec-butyl alcohol': 'CCC(C)O', 'tert-butyl alcohol': 'CC(C)(C)O',
    'amyl alcohol': 'CCCCCO', 'isoamyl alcohol': 'CC(C)CCO', 'hexanol': 'CCCCCCO', 'heptanol': 'CCCCCCCO', 'octanol': 'CCCCCCCCO', '2-ethylhexanol': 'CCCCC(CC)CO', 'decanol': 'CCCCCCCCCCO', 'lauryl alcohol': 'CCCCCCCCCCCCO',
    'cetyl alcohol': 'CCCCCCCCCCCCCCCCO', 'stearyl alcohol': 'CCCCCCCCCCCCCCCCCCO', 'benzhydrol': 'OC(c1ccccc1)c1ccccc1', 'triphenylmethanol': 'OC(c1ccccc1)(c1ccccc1)c1ccccc1', 'cyclopentanol': 'OC1CCCC1',
    'cycloheptanol': 'OC1CCCCCC1', 'menthol': 'CC(C)C1CCC(C)CC1O', 'furfuryl alcohol': 'OCc1ccco1', 'tetrahydrofurfuryl alcohol': 'OCC1CCCO1', 'cinnamyl alcohol': 'OCC=Cc1ccccc1', 'propargyl alcohol': 'C#CCO',
    'butynediol': 'OCC#CCO', '1,3-propanediol': 'OCCCO', '1,4-butanediol': 'OCCCCO', '1,5-pentanediol': 'OCCCCCO', '1,6-hexanediol': 'OCCCCCCO', '1,2-propanediol': 'CC(O)CO', '2,3-butanediol': 'CC(O)C(C)O',
    'pinacol': 'CC(C)(O)C(C)(C)O', 'neopentyl glycol': 'CC(C)(CO)CO', 'pentaerythritol': 'OCC(CO)(CO)CO', 'trimethylolpropane': 'CCC(CO)(CO)CO', 'erythritol': 'OCC(O)C(O)CO', 'xylitol': 'OCC(O)C(O)C(O)CO',
    'sorbitol': 'OCC(O)C(O)C(O)C(O)CO', 'mannitol': 'OCC(O)C(O)C(O)C(O)CO', 'inositol': 'OC1C(O)C(O)C(O)C(O)C1O', 'dulcitol': 'OCC(O)C(O)C(O)C(O)CO', 'sucrose': 'OCC1OC(OC2(CO)OC(CO)C(O)C2O)C(O)C(O)C1O',
    'lactose': 'OCC1OC(OC2C(O)C(O)C(O)OC2CO)C(O)C(O)C1O', 'maltose': 'OCC1OC(OC2C(O)C(O)C(O)OC2CO)C(O)C(O)C1O', 'cellobiose': 'OCC1OC(OC2C(O)C(O)C(O)OC2CO)C(O)C(O)C1O', 'trehalose': 'OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O',
    'glucose': 'OCC1OC(O)C(O)C(O)C1O', 'dextrose': 'OCC1OC(O)C(O)C(O)C1O', 'galactose': 'OCC1OC(O)C(O)C(O)C1O', 'mannose': 'OCC1OC(O)C(O)C(O)C1O', 'fructose': 'OCC1(O)OCC(O)C(O)C1O', 'ribose': 'OCC1OC(O)C(O)C1O',
    'deoxyribose': 'OCC1OC(O)CC1O', 'xylose': 'OC1COC(O)C(O)C1O', 'arabinose': 'OC1COC(O)C(O)C1O', 'sorbose': 'OCC1(O)OCC(O)C(O)C1O', 'fucose': 'CC1OC(O)C(O)C(O)C1O', 'rhamnose': 'CC1OC(O)C(O)C(O)C1O',
    'glucosamine': 'NC1C(O)OC(CO)C(O)C1O', 'n-acetylglucosamine': 'CC(=O)NC1C(O)OC(CO)C(O)C1O', 'gluconic acid': 'OCC(O)C(O)C(O)C(O)C(=O)O', 'glucuronic acid': 'OC1C(O)C(O)OC(C(=O)O)C1O', 'sialic acid': 'CC(=O)NC1C(O)CC(O)(OC1C(O)C(O)CO)C(=O)O',
    'glycerol': 'OCC(O)CO', 'glycerin': 'OCC(O)CO', 'glycerine': 'OCC(O)CO', 'ethylene glycol': 'OCCO', 'antifreeze': 'OCCO', 'vinegar': 'CC(=O)O', 'dry ice': 'O=C=O', 'laughing gas': '[N-]=[N+]=O', 'nitrous oxide': '[N-]=[N+]=O',
    'carbon monoxide': '[C-]#[O+]', 'sulfur dioxide': 'O=S=O', 'sulfur trioxide': 'O=S(=O)=O', 'nitric oxide': '[N]=O', 'nitrogen dioxide': 'O=[N]=O', 'ammonia': 'N', 'hydrogen sulfide': 'S', 'hydrogen chloride': 'Cl',
    'hydrogen bromide': 'Br', 'hydrogen iodide': 'I', 'hydrogen fluoride': 'F', 'sulfuric acid': 'OS(O)(=O)=O', 'nitric acid': 'O[N+](=O)[O-]', 'perchloric acid': 'OCl(=O)(=O)=O', 'hypochlorous acid': 'OCl',
    'sodium hydroxide': '[Na+].[OH-]', 'potassium hydroxide': '[K+].[OH-]', 'sodium bicarbonate': '[Na+].OC([O-])=O', 'sodium carbonate': '[Na+].[Na+].[O-]C([O-])=O', 'sodium acetate': '[Na+].CC([O-])=O',
    'sodium chloride': '[Na+].[Cl-]', 'potassium carbonate': '[K+].[K+].[O-]C([O-])=O', 'sodium ethoxide': '[Na+].CC[O-]', 'sodium methoxide': '[Na+].C[O-]', 'potassium tert-butoxide': '[K+].CC(C)(C)[O-]',
    'sodium hydride': '[Na+].[H-]', 'lithium aluminum hydride': '[Li+].[AlH4-]', 'sodium borohydride': '[Na+].[BH4-]', 'lithium diisopropylamide': '[Li+].CC(C)[N-]C(C)C', 'butyllithium': '[Li]CCCC', 'n-butyllithium': '[Li]CCCC',
    'methyllithium': '[Li]C', 'phenyllithium': '[Li]c1ccccc1', 'methylmagnesium bromide': 'C[Mg]Br', 'phenylmagnesium bromide': 'Br[Mg]c1ccccc1', 'ethylmagnesium bromide': 'CC[Mg]Br', 'vinylmagnesium bromide': 'C=C[Mg]Br',
    'diethylzinc': 'CC[Zn]CC', 'trimethylaluminum': 'C[Al](C)C', 'triethylaluminum': 'CC[Al](CC)CC', 'tributyltin hydride': 'CCCC[SnH](CCCC)CCCC', 'tetramethylsilane': 'C[Si](C)(C)C', 'trimethylsilyl chloride': 'C[Si](C)(C)Cl',
    'hexamethyldisiloxane': 'C[Si](C)(C)O[Si](C)(C)C', 'hexamethyldisilazane': 'C[Si](C)(C)N[Si](C)(C)C', 'tetraethyl orthosilicate': 'CCO[Si](OCC)(OCC)OCC', 'tetraethoxysilane': 'CCO[Si](OCC)(OCC)OCC',
    'trimethyl borate': 'COB(OC)OC', 'triethylborane': 'CCB(CC)CC', 'borane': 'B', 'diborane': 'BB', 'catecholborane': 'B1Oc2ccccc2O1', 'pinacolborane': 'CC1(C)OBOC1(C)C', 'bis(pinacolato)diboron': 'CC1(C)OB(OC1(C)C)B1OC(C)(C)C(C)(C)O1',
    'ferrocene': '[Fe+2].c1cc[cH-]c1.c1cc[cH-]c1', 'cisplatin': 'N.N.Cl[Pt]Cl', 'tetrakis(triphenylphosphine)palladium': 'c1ccc(cc1)P(c1ccccc1)c1ccccc1', 'dimethyl sulfate': 'COS(=O)(=O)OC', 'diethyl sulfate': 'CCOS(=O)(=O)OCC',
    'methyl triflate': 'COS(=O)(=O)C(F)(F)F', 'methyl iodide': 'CI', 'ethyl iodide': 'CCI', 'trimethyloxonium tetrafluoroborate': 'C[O+](C)C.F[B-](F)(F)F', 'trityl cation': 'c1ccc(cc1)[C+](c1ccccc1)c1ccccc1',
    'benzyne': 'C1=CC#CC=C1', 'carbene': '[CH2]', 'methyl radical': '[CH3]', 'tert-butyl cation': 'C[C+](C)C', 'methyl cation': '[CH3+]', 'hydronium': '[OH3+]', 'hydroxide': '[OH-]', 'ammonium': '[NH4+]', 'proton': '[H+]',
    'tms': 'C[Si](C)(C)C', 'dioxane': 'C1COCCO1', '1,3-dioxane': 'C1COCOC1', 'dioxolane': 'C1COCO1', 'trioxane': 'C1OCOCO1', 'crown ether 12-crown-4': 'C1COCCOCCOCCO1', '15-crown-5': 'C1COCCOCCOCCOCCO1', '18-crown-6': 'C1COCCOCCOCCOCCOCCO1',
    'cryptand 222': 'C1COCCOCCN2CCOCCOCCN1CCOCCOCC2', 'cyclodextrin unit': 'OCC1OC(O)C(O)C(O)C1O', 'tetrahydropyran': 'C1CCOCC1', 'dihydropyran': 'C1=COCCC1', '3,4-dihydro-2h-pyran': 'C1=COCCC1', '2,3-dihydrofuran': 'C1=COCC1',
    'furan': 'c1ccoc1', 'thiophene': 'c1ccsc1', 'pyrrole': 'c1cc[nH]c1', 'pyridine': 'c1ccncc1', 'pyrimidine': 'c1cncnc1', 'pyrazine': 'c1cnccn1', 'pyridazine': 'c1ccnnc1', 'triazine': 'c1ncncn1', 'tetrazine': 'c1nncnn1',
    'oxazole': 'c1cocn1', 'isoxazole': 'c1cnoc1', 'thiazole': 'c1cscn1', 'isothiazole': 'c1cnsc1', 'pyrazole': 'c1cn[nH]c1', 'triazole': 'c1nc[nH]n1', 'tetrazole': 'c1nnn[nH]1', 'oxadiazole': 'c1nnco1', 'thiadiazole': 'c1nncs1',
    'indole': 'c1ccc2[nH]ccc2c1', 'isoindole': 'c1ccc2c[nH]cc2c1', 'indolizine': 'c1ccn2cccc2c1', 'benzofuran': 'c1ccc2occc2c1', 'isobenzofuran': 'c1ccc2cocc2c1', 'benzothiophene': 'c1ccc2sccc2c1', 'quinoline': 'c1ccc2ncccc2c1',
    'isoquinoline': 'c1ccc2cnccc2c1', 'acridine': 'c1ccc2nc3ccccc3cc2c1', 'phenazine': 'c1ccc2nc3ccccc3nc2c1', 'phenothiazine': 'c1ccc2Sc3ccccc3Nc2c1', 'phenoxazine': 'c1ccc2Oc3ccccc3Nc2c1', 'pteridine': 'c1cnc2ncncc2n1',
    'purine': 'c1ncc2[nH]cnc2n1', 'cinnoline': 'c1ccc2nnccc2c1', 'phthalazine': 'c1ccc2cnncc2c1', 'carboline': 'c1ccc2c(c1)[nH]c1cnccc12', 'porphine': 'c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3',
  };

  /* ---- generated IUPAC names ---- */
  const GEN = { names: [], map: new Map(), ready: false, started: false, onProgress: null };

  const ALKYL = { 1: 'C', 2: 'CC', 3: 'CCC', 4: 'CCCC', 5: 'CCCCC', 6: 'CCCCCC' };
  // chain of n carbons, with optional branch strings per position and bond orders after each position
  function chain(n, br, bo) {
    let s = '';
    for (let i = 1; i <= n; i++) {
      s += 'C' + (br && br[i] ? br[i].map(b => '(' + b + ')').join('') : '');
      if (i < n) s += bo && bo[i] === 2 ? '=' : bo && bo[i] === 3 ? '#' : '';
    }
    return s;
  }
  const add = (smi) => { const k = smi; if (!GEN.seen.has(k)) { GEN.seen.add(k); GEN.todo.push(smi); } };

  function enumerate() {
    GEN.seen = new Set(); GEN.todo = [];
    const SUB = ['F', 'Cl', 'Br', 'I', 'O', 'S', 'N', '[N+](=O)[O-]', 'C#N', 'OC', 'C(=O)O', 'C(=O)N', 'C=O', 'C(=O)C', 'OC(C)=O'];
    // 1. unbranched chains: alkane / alkene / alkyne at each position, and one substituent at each position
    for (let n = 1; n <= 12; n++) {
      add(chain(n));
      for (let i = 1; i < n; i++) { add(chain(n, null, { [i]: 2 })); add(chain(n, null, { [i]: 3 })); }
      for (let i = 1; i < n && n <= 8; i++) for (let j = i + 2; j < n; j++) add(chain(n, null, { [i]: 2, [j]: 2 }));
      if (n <= 10) for (let i = 1; i <= n; i++) for (const x of SUB) { const br = { [i]: [x] }; add(chain(n, br)); }
    }
    // 2. terminal groups: acids, aldehydes, nitriles, amides, acid halides, esters, anhydrides, ketones
    for (let n = 1; n <= 14; n++) {
      const R = ALKYL[n] || 'C'.repeat(n);
      for (const g of ['C(=O)O', 'C=O', 'C#N', 'C(=O)N', 'C(=O)Cl', 'C(=O)Br', 'C(=O)NC', 'C(=O)N(C)C', 'C(=O)OC', 'C(=O)OCC', 'C(=O)OC(C)C', 'C(=O)OC(C)(C)C', 'C(=O)OC=C', 'C(=O)OOC', 'C(=S)N', 'C(=O)S', 'C(N)=N', 'S(=O)(=O)O', 'S(=O)(=O)Cl', 'S(=O)(=O)N', 'OS(=O)(=O)C', 'OC(=O)C', 'OC(=O)CC', 'OC(=O)c1ccccc1', 'C(=O)c1ccccc1', 'c1ccccc1', 'Cc1ccccc1'])
        add(R + g);
      if (n <= 8) { add('O=C(O)' + R + 'C(=O)O'); add('N#C' + R + 'C#N'); add('O=C' + R + 'C=O'); }
    }
    // 3. ketones, ethers, esters, amines, thioethers, anhydrides over small R groups
    for (let n = 1; n <= 7; n++) for (let m = n; m <= 7; m++) {
      const A = 'C'.repeat(n), B = 'C'.repeat(m);
      add(A + 'C(=O)' + B); add(A + 'O' + B); add(A + 'S' + B); add(A + 'NC'.slice(0, 1) + B);
      add(A + 'C(=O)O' + B); add(B + 'C(=O)O' + A); add(A + 'C(=O)OC(=O)' + B); add(A + 'NC' + B);
      add(A + 'N(C)' + B); add(A + 'C(=O)N' + B); add(A + 'C(=O)NC' + B);
    }
    for (let n = 1; n <= 6; n++) for (let m = 1; m <= 6; m++) { add('C'.repeat(n) + 'N(' + 'C'.repeat(m) + ')C'); add('C'.repeat(n) + 'N(CC)' + 'C'.repeat(m)); }
    // 4. branched alkanes and their substituted forms: methyl / ethyl branches on chains up to C9
    for (let n = 3; n <= 9; n++) {
      for (let i = 2; i < n; i++) {
        for (const b of ['C', 'CC', 'C(C)C', 'CCC']) {
          add(chain(n, { [i]: [b] }));
          for (let j = i; j < n; j++) { add(chain(n, { [i]: [b], [j]: ['C'] })); if (j > i || b === 'C') add(chain(n, { [i]: ['C', 'C'] })); }
          if (n <= 7) for (const x of ['O', 'Cl', 'Br', 'N', 'S']) for (let k = 1; k <= n; k++) add(chain(n, { [i]: [b], [k]: [x] }));
          if (n <= 7) for (let k = 1; k < n; k++) { add(chain(n, { [i]: [b] }, { [k]: 2 })); add(chain(n, { [i]: [b] }, { [k]: 3 })); }
        }
        for (let j = i; j < n; j++) for (let k = j; k < n; k++) if (n <= 8) add(chain(n, { [i]: ['C'], [j]: ['C'], [k]: ['C'] }));
      }
    }
    // 5. di- and poly-substituted chains: diols, dihalides, diamines, hydroxy acids, amino acids, haloacids, keto acids
    for (let n = 2; n <= 8; n++) {
      for (const x of ['O', 'Cl', 'Br', 'F', 'N', 'S', 'I']) for (let i = 1; i <= n; i++) for (let j = i; j <= n; j++) {
        if (i === j && x !== 'Cl' && x !== 'Br' && x !== 'F' && x !== 'I') continue;
        const br = {}; br[i] = [x]; br[j] = (br[j] || []).concat(i === j ? [x] : [x]);
        add(chain(n, br));
      }
      for (let i = 1; i <= n; i++) for (const x of ['O', 'N', 'Cl', 'Br', 'F', 'C']) { const body = chain(n - 1, { [i]: i <= n - 1 ? [x] : null }); if (i <= n - 1) add(body + 'C(=O)O'); }
      for (let i = 1; i < n; i++) { add(chain(n - 1, { [i]: ['=O'] }) + 'C(=O)O'); add(chain(n - 1, { [i]: ['=O'] }) + 'C=O'); }
      for (let i = 1; i <= n; i++) for (let j = i + 1; j <= n; j++) { add(chain(n, { [i]: ['=O'], [j]: ['=O'] })); add(chain(n, { [i]: ['O'], [j]: ['=O'] })); add(chain(n, { [i]: ['N'], [j]: ['=O'] })); }
    }
    // 6. cycloalkanes, cycloalkenes: parent, one substituent, two methyls, ketones, alcohols, acids
    for (let n = 3; n <= 8; n++) {
      const ring = (subs, dbl) => {
        let s = '';
        for (let i = 0; i < n; i++) { s += 'C' + (i === 0 ? '1' : '') + (subs[i] ? subs[i].map(b => '(' + b + ')').join('') : '') + (i === n - 1 ? '1' : (dbl === i ? '=' : '')); }
        return s;
      };
      const ring1 = (subs, dbl) => {                // ring closure after the last atom
        let s = '';
        for (let i = 0; i < n; i++) s += 'C' + (i === 0 ? '1' : '') + (subs[i] ? subs[i].map(b => '(' + b + ')').join('') : '') + (i === n - 1 ? '1' : '') + (dbl === i && i < n - 1 ? '=' : '');
        return s;
      };
      void ring;
      add(ring1([])); add(ring1([], 0)); if (n >= 5) add(ring1([], 1)); if (n >= 6) { add(ring1([], 2)); }
      for (const x of SUB.concat(['C', 'CC', 'C(C)C', 'C(C)(C)C', 'C=C', 'c1ccccc1'])) { add(ring1([[x]])); if (n >= 4) add(ring1([[x]], 0)); }
      for (let j = 1; j < n; j++) for (const [a, b] of [['C', 'C'], ['O', 'O'], ['Cl', 'Cl'], ['Br', 'Br'], ['O', 'C'], ['N', 'C'], ['Cl', 'C'], ['O', '=O'], ['C', '=O'], ['C(=O)O', 'O'], ['N', 'N'], ['O', 'N'], ['Cl', 'O'], ['Br', 'C'], ['F', 'F'], ['C(=O)O', 'C(=O)O']]) { const s = []; s[0] = [a]; s[j] = [b]; add(ring1(s)); }
      for (let j = 1; j < n; j++) for (let k = j + 1; k < n; k++) { const s = []; s[0] = ['C']; s[j] = ['C']; s[k] = ['C']; add(ring1(s)); }
    }
    // 7. benzene: one, two (o / m / p) and three substituents from the common list
    const B = ['C', 'CC', 'C(C)C', 'C(C)(C)C', 'F', 'Cl', 'Br', 'I', 'O', 'OC', 'N', 'N(C)C', '[N+](=O)[O-]', 'C#N', 'C=O', 'C(=O)C', 'C(=O)O', 'C(=O)N', 'C(=O)OC', 'S', 'S(=O)(=O)O', 'C(F)(F)F', 'C=C', 'C#C', 'OC(C)=O', 'NC(C)=O', 'CO', 'CCl', 'CBr', 'Cc1ccccc1', 'c1ccccc1', 'Oc1ccccc1'];
    const benz = (subs) => { const a = ['c', 'c', 'c', 'c', 'c', 'c']; let s = ''; for (let i = 0; i < 6; i++) s += a[i] + (i === 0 ? '1' : '') + (subs[i] ? '(' + subs[i] + ')' : '') + (i === 5 ? '1' : ''); return s; };
    add(benz([]));
    for (const x of B) add(benz([x]));
    for (let a = 0; a < B.length; a++) for (let b = a; b < B.length; b++) for (const pos of [1, 2, 3]) { const s = []; s[0] = B[a]; s[pos] = B[b]; add(benz(s)); }
    for (const x of B.slice(0, 14)) for (const y of B.slice(0, 14)) for (const z of ['C', 'Cl', 'Br', '[N+](=O)[O-]', 'O', 'N']) { const s = []; s[0] = x; s[1] = y; s[2] = z; if (x !== y || x === z) add(benz(s)); const t = []; t[0] = x; t[2] = y; t[4] = z; add(benz(t)); }
    // 8. aromatic heterocycles with one substituent
    for (const ring of ['c1ccncc1', 'c1ccoc1', 'c1ccsc1', 'c1cc[nH]c1']) for (const x of ['C', 'CC', 'Cl', 'Br', 'F', 'O', 'N', 'C=O', 'C(=O)O', 'C(=O)C', '[N+](=O)[O-]', 'C#N', 'OC', 'S']) {
      add(ring.replace('c1', 'c1(' + x + ')')); add(ring.replace('cc', 'cc(' + x + ')')); add(ring.replace('cc', 'c(' + x + ')c'));
    }
    // 9. polyfunctional small molecules: haloforms, polyhalides, glycols, hydroxy-ketones, etc.
    for (let n = 1; n <= 6; n++) for (const x of ['F', 'Cl', 'Br', 'I']) { const R = 'C'.repeat(n); add(x + R + x); add(x + R.slice(0, -1) + 'C(' + x + ')' + x); add(R.slice(0, -1) + 'C(' + x + ')(' + x + ')' + x); }
    return GEN.todo;
  }

  function nameOf(smi) {
    try {
      const g = C.parseSmiles(smi);
      if (!g.atoms.length) return null;
      const n = C.iupacName(g);
      if (!n || typeof n !== 'string' || /[?]/.test(n)) return null;
      /* a name the parser reads as some OTHER molecule is wrong: drop it. One it cannot read at all is kept (the
         suggestion carries its structure) */
      try { if (C.canonicalKey(C.parseIupacName(n)) !== C.canonicalKey(g)) return null; } catch (e) { /* unreadable: keep */ }
      return n;
    } catch (e) { return null; }
  }

  function start() {
    if (GEN.started) return; GEN.started = true;
    let todo; try { todo = enumerate(); } catch (e) { return; }
    let i = 0;
    const seen = new Set();
    const slice = () => {
      const t0 = performance.now();
      while (i < todo.length && performance.now() - t0 < 12) {
        const smi = todo[i++];
        const nm = nameOf(smi);
        if (nm) { const k = nm.toLowerCase(); if (!seen.has(k)) { seen.add(k); GEN.names.push({ label: nm, smiles: smi }); GEN.map.set(k, smi); } }
      }
      if (i < todo.length) setTimeout(slice, 8); else { GEN.ready = true; todo.length = 0; }
    };
    setTimeout(slice, 300);
  }

  // add the extra common names (never overriding one the app already has)
  if (C.COMMON_NAMES) for (const k of Object.keys(EXTRA_COMMON)) if (!(k in C.COMMON_NAMES)) C.COMMON_NAMES[k] = EXTRA_COMMON[k];

  window.IupacGen = { start, get names() { return GEN.names; }, smilesFor: k => GEN.map.get(k.toLowerCase()), get ready() { return GEN.ready; } };
  // begin once the page is idle
  (window.requestIdleCallback || (f => setTimeout(f, 1500)))(start);
})();
