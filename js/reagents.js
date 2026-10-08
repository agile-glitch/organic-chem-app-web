/* Reagents: reading a reagent the way a chemist writes it, and mapping it to an RXN rule.
   Lifted from organic-chem-tool/src/shell.html (RULE_REAGENTS, REAGENT_PHRASE, resolvePredReagent,
   moleculeRule, ruleApplies, ACRONYMS, PRED_COND). One change: the drawn molecule the old app used to
   disambiguate a reagent (window.ChemApp.sketch) is now window.REAGENT_SUBSTRATE, set by the page.
   Needs the Chem and RXN globals. */
(() => {
const RULE_REAGENTS = {
  /* --- the modern catalog. Names are matched with spaces, dots, hyphens and brackets stripped,
         so 'Dess-Martin' and 'dess martin' both arrive here as 'dessmartin'. --- */
  dmp:      { smiles: null, names: ['dmp', 'dessmartin', 'dessmartinperiodinane', 'periodinane'], formulas: [] },
  swern:    { smiles: null, names: ['swern', 'swernoxidation', 'dmso/oxalylchloride', 'oxalylchloride/dmso', 'cocl22/dmso'], formulas: [] },
  tempo:    { smiles: null, names: ['tempo', 'tempo/naocl', 'tempo/bleach', 'anelli', 'anellioxidation'], formulas: [] },
  jones:    { smiles: null, names: ['jones', 'jonesreagent', 'jonesoxidation', 'cro3/h2so4', 'chromicacid'], formulas: [] },
  soCl2:           { smiles: null, names: ['socl2', 'thionylchloride', 'oxalylchloride', 'cocl2cocl', 'c2o2cl2'], formulas: ['socl2'] },
  pcl5:            { smiles: null, names: ['pcl5', 'pcl3', 'phosphoruspentachloride', 'phosphorustrichloride'], formulas: ['pcl5', 'pcl3'] },
  ch2n2:           { smiles: null, names: ['ch2n2', 'diazomethane', 'tmsdiazomethane', 'tmschn2'], formulas: ['ch2n2'] },
  radpoly:         { smiles: null, names: ['radpoly', 'peroxideinitiator', 'radicalinitiator', 'aibn', 'bpo', 'radicalpolymerisation', 'radicalpolymerization', 'roor'], formulas: [] },
  cationic:        { smiles: null, names: ['cationicpolymerisation', 'cationicpolymerization', 'bf3h2o'], formulas: [] },
  anionic:         { smiles: null, names: ['anionicpolymerisation', 'anionicpolymerization', 'livingpolymerisation'], formulas: [] },
  ziegler:         { smiles: null, names: ['ziegler', 'zieglernatta', 'ticl4et2alcl', 'zieglernattacatalyst'], formulas: [] },
  nabh3cn:         { smiles: null, names: ['nabh3cn', 'sodiumcyanoborohydride', 'nabh(oac)3', 'sodiumtriacetoxyborohydride', 'reductiveamination', 'nacnbh3'], formulas: ['nabh3cn'] },
  tautomer:        { smiles: null, names: ['tautomer', 'ketoenol', 'tautomerisation', 'tautomerization', 'hplusoroh', 'traceacid', 'tracebase', 'enolisation', 'enolization'], formulas: [] },
  lda:             { smiles: null, names: ['lda', 'lithiumdiisopropylamide', 'kineticenolate', 'ldathf'], formulas: [] },
  enolate_base:    { smiles: null, names: ['thermodynamicenolate', 'naoetenolate', 'enolate'], formulas: [] },
  thermal:         { smiles: null, names: ['thermal', 'decarboxylation', 'heatonly', 'justheat', 'heatalone'], formulas: [] },
  ar_cl2:          { smiles: null, names: ['cl2fecl3', 'cl2alcl3', 'aromaticchlorination', 'chlorinationaromatic'], formulas: [] },
  ar_i2:           { smiles: null, names: ['i2hno3', 'aromaticiodination'], formulas: [] },
  epox_h2o:        { smiles: null, names: ['epoxidehydrolysis', 'h2oh2so4epoxide', 'antidiol'], formulas: [] },
  epox_hx:         { smiles: null, names: ['halohydrin', 'hxepoxide', 'epoxidehx'], formulas: [] },
  epox_roh:        { smiles: null, names: ['ch3ohh2so4', 'methanolacid', 'epoxidemethanol'], formulas: [] },
  clemmensen:      { smiles: null, names: ['clemmensen', 'znhghcl', 'zinchgamalgam', 'clemmensenreduction', 'zinamalgam'], formulas: [] },
  wolffkishner:    { smiles: null, names: ['wolffkishner', 'n2h4koh', 'nh2nh2koh', 'hydrazinekoh', 'wolffkishnerreduction', 'hydrazine'], formulas: [] },
  ozonolysis:      { smiles: null, names: ['o3', 'ozone', 'ozonolysis', 'o3me2s', 'o3zn', 'ozonolysisreductive'], formulas: [] },
  ozonolysis_ox:   { smiles: null, names: ['o3h2o2', 'ozonolysisoxidative'], formulas: [] },
  alkynehyd:       { smiles: null, names: ['hgso4', 'hgso4h2so4', 'alkynehydration', 'mercurichydration'], formulas: [] },
  alkynehydb:      { smiles: null, names: ['sia2bh', 'disiamylborane', '9bbn', 'alkyneantimarkovnikov'], formulas: [] },
  allylbr:         { smiles: null, names: ['nbs', 'nbromosuccinimide', 'allylicbromination', 'benzylicbromination', 'br2hv', 'br2light'], formulas: [] },
  aqacid:          { smiles: null, names: ['h3o', 'h3oplus', 'aqueousacid', 'diluteacid', 'hplus', 'h2oh', 'workupacid'], formulas: [] },
  ag2o:            { smiles: null, names: ['ag2o', 'silveroxide', 'tollens', 'tollensreagent'], formulas: [] },
  alkoxyhg:        { smiles: null, names: ['hgoac2', 'mercuricacetate', 'alkoxymercuration', 'oxymercuration'], formulas: [] },
  nano2:        { smiles: null, names: ['nano2', 'sodiumnitrite', 'nitrousacid', 'hno2', 'diazotisation', 'diazotization', 'nano2hcl'], formulas: ['nano2', 'hno2'] },
  sandmeyer_cl: { smiles: null, names: ['cucl', 'coppericchloride', 'copperichloride', 'sandmeyercl', 'sandmeyerchloride'], formulas: ['cucl'] },
  sandmeyer_br: { smiles: null, names: ['cubr', 'copperibromide', 'sandmeyerbr', 'sandmeyerbromide'], formulas: ['cubr'] },
  sandmeyer_cn: { smiles: null, names: ['cucn', 'coppericyanide', 'sandmeyercn', 'sandmeyercyanide'], formulas: ['cucn'] },
  diazo_i:      { smiles: null, names: ['ki', 'potassiumiodide', 'aryliodide'], formulas: ['ki'] },
  diazo_f:      { smiles: null, names: ['hbf4', 'tetrafluoroboricacid', 'balzschiemann', 'balzschiemannreaction'], formulas: ['hbf4'] },
  diazo_oh:     { smiles: null, names: ['warmwater', 'h2owarm', 'diazoniumhydrolysis'], formulas: [] },
  diazo_h:      { smiles: null, names: ['h3po2', 'hypophosphorousacid', 'phosphinicacid', 'reductivedeamination'], formulas: ['h3po2'] },
  mno2:     { smiles: null, names: ['mno2', 'manganesedioxide', 'activatedmno2', 'manganeseivoxide'], formulas: ['mno2'] },
  dibal:    { smiles: null, names: ['dibal', 'dibalh', 'diisobutylaluminiumhydride', 'diisobutylaluminumhydride'], formulas: [] },
  lindlar:  { smiles: null, names: ['lindlar', 'lindlarcatalyst', 'h2/lindlar', 'pd/caco3', 'poisonedpalladium'], formulas: [] },
  na_nh3:   { smiles: null, names: ['na/nh3', 'nanh3', 'sodiuminammonia', 'li/nh3', 'dissolvingmetal', 'dissolvingmetalreduction', 'birch'], formulas: [] },
  nitrored: { smiles: null, names: ['fe/hcl', 'zn/hcl', 'sn/hcl', 'sncl2', 'nitroreduction', 'reducethenitro', 'h2/pdc', 'h2/pd', 'h2pdc', 'hydrogenolysis'], formulas: [] },
  boc_on:   { smiles: null, names: ['boc2o', 'bocanhydride', 'ditertbutyldicarbonate', 'bocprotection', 'protectwithboc'], formulas: [] },
  boc_off:  { smiles: null, names: ['tfa/ch2cl2', 'bocremoval', 'bocdeprotection', 'removeboc', 'hcl/dioxane', 'tfa'], formulas: [] },
  cbz_on:   { smiles: null, names: ['cbzcl', 'benzylchloroformate', 'cbzprotection', 'protectwithcbz'], formulas: [] },
  cbz_off:  { smiles: null, names: ['cbzremoval', 'cbzdeprotection', 'removecbz', 'hydrogenolysis', 'h2/pdc', 'h2/pd', 'h2pdc'], formulas: [] },
  fmoc_on:  { smiles: null, names: ['fmoccl', 'fmocosu', 'fmocprotection', 'protectwithfmoc'], formulas: [] },
  fmoc_off: { smiles: null, names: ['piperidine/dmf', 'fmocremoval', 'fmocdeprotection', 'removefmoc', '20piperidine'], formulas: [] },
  tbs_on:   { smiles: null, names: ['tbscl', 'tbdmscl', 'tertbutyldimethylsilylchloride', 'silylprotection', 'tbsotf'], formulas: [] },
  tbs_off:  { smiles: null, names: ['tbaf', 'tetrabutylammoniumfluoride', 'hf/pyridine', 'desilylation', 'removetbs', 'removethesilyl'], formulas: [] },
  bn_on:    { smiles: null, names: ['bnbr', 'benzylbromide', 'benzylprotection', 'bncl', 'benzylchloride'], formulas: [] },
  bn_off:   { smiles: null, names: ['debenzylation', 'removethebenzyl', 'removebn', 'h2/pdc', 'h2/pd', 'h2pdc'], formulas: [] },
  thp_on:   { smiles: null, names: ['dhp', 'dihydropyran', '34dihydro2hpyran', 'thpprotection'], formulas: [] },
  ac_on:    { smiles: null, names: ['ac2o/pyridine', 'aceticanhydride/pyridine', 'acetylation', 'acetylate'], formulas: [] },
  ms_on:    { smiles: null, names: ['mscl', 'methanesulfonylchloride', 'mesylchloride', 'mesylation'], formulas: [] },
  ts_on:    { smiles: null, names: ['tscl', 'tosylchloride', 'ptoluenesulfonylchloride', 'tosylation'], formulas: [] },
  acetal_on:{ smiles: null, names: ['ethyleneglycol', 'ethyleneglycol/tsoh', 'acetalprotection', 'protectthecarbonyl', 'dioxolane'], formulas: [] },
  acetal_off:{ smiles: null, names: ['acetalhydrolysis', 'removetheacetal', 'aqueousacid', 'diluteacid'], formulas: [] },
  cyclise:  { smiles: null, names: ['cyclise', 'cyclize', 'cyclisation', 'cyclization', 'closethering', 'intramolecular', 'lactamisation', 'lactonisation'], formulas: [] },
  hbr:   { smiles: 'Br',  names: ['hbr', 'hydrogenbromide', 'hydrobromicacid'], formulas: ['hbr', 'brh'] },
  hcl:   { smiles: 'Cl',  names: ['hcl', 'hydrogenchloride', 'hydrochloricacid'], formulas: ['hcl', 'clh'] },
  hi:    { smiles: 'I',   names: ['hi', 'hydrogeniodide', 'hydroiodicacid'], formulas: ['hi', 'ih'] },
  cl2:   { smiles: 'ClCl', names: ['cl2', 'chlorine', 'dichlorine'], formulas: ['cl2'] },
  br2:   { smiles: 'BrBr', names: ['br2', 'bromine', 'dibromine'], formulas: ['br2'] },
  sulfonation: { smiles: 'O=S(=O)=O', names: ['so3', 'sulfurtrioxide', 'sulfonation', 'fumingsulfuricacid', 'oleum', 'h2so4/so3', 'so3/h2so4'], formulas: ['so3', 'o3s'] },
  nitration: { smiles: 'O[N+](=O)[O-]', names: ['hno3', 'nitricacid', 'nitration', 'hno3/h2so4', 'nitratingmixture'], formulas: ['hno3', 'hno3/h2so4'] },
  naocl: { smiles: '[Na+].[O-]Cl', names: ['naocl', 'bleach', 'sodiumhypochlorite', 'hypochlorite', 'naocl/acoh', 'bleach/aceticacid', 'naocl/aceticacid', 'clorox', 'householdbleach', 'hocl', 'hypochlorousacid', 'bleachoxidation', 'haloform', 'haloformreaction', 'naocl/naoh'], formulas: ['naocl', 'naclo', 'clnao', 'hocl', 'hclo'] },
  etard: { smiles: 'O=[Cr](=O)(Cl)Cl', names: ['cro2cl2', 'chromylchloride', 'etard', 'étard', 'etardreagent', 'étardreagent', 'etardreaction', 'étardreaction', 'etardoxidation', 'chromiumdioxidedichloride'], formulas: ['cro2cl2', 'cl2cro2', 'cro2cl2/cs2'] },
  ar_br2: { smiles: 'BrBr', names: ['br2/febr3', 'febr3', 'br2febr3', 'ringbromination', 'bromination'], formulas: ['br2/febr3', 'febr3'] },
  hydration: { smiles: 'O', names: ['water', 'hydration', 'h2o/h2so4', 'acidwater'], formulas: ['h2o', 'oh2'] },
  h2:    { names: ['h2', 'hydrogen', 'h2/pd', 'h2pd', 'h2/pdc', 'h2pdc', 'hydrogenation'], formulas: ['h2'] },
  naoh:  { smiles: '[OH-]', names: ['naoh', 'sodiumhydroxide', 'hydroxide', 'koh', 'potassiumhydroxide'], formulas: ['naoh', 'hnao', 'koh', 'hko'] },
  naoet: { smiles: 'CC[O-]', names: ['naoet', 'sodiumethoxide', 'ethoxide', 'etona'], formulas: ['c2h5nao', 'c2h5ona', 'nac2h5o'] },
  tbuok: { smiles: 'CC(C)(C)[O-]', names: ['tbuok', 'kotbu', 'potassiumtertbutoxide', 'potassiumtert-butoxide', 'tertbutoxide', 'tbutoxide'], formulas: ['c4h9ko', 'koc4h9'] },
  nacn:  { smiles: '[C-]#N', names: ['nacn', 'sodiumcyanide', 'cyanide', 'kcn'], formulas: ['nacn', 'cnna', 'kcn', 'cnk'] },
  pbr3:  { smiles: 'P(Br)(Br)Br', names: ['pbr3', 'phosphorustribromide'], formulas: ['pbr3', 'br3p'] },
  mcpba: { smiles: 'OOC(=O)c1cccc(Cl)c1', names: ['mcpba', 'm-cpba', 'metachloroperoxybenzoicacid', 'metachloroperbenzoicacid', 'peracid'], formulas: ['c7h5clo3'] },
  oso4:  { names: ['oso4', 'osmiumtetroxide', 'osmiumtetraoxide'], formulas: ['oso4', 'o4os'] },
  kmno4: { names: ['kmno4', 'potassiumpermanganate', 'permanganate', 'cro3', 'jonesreagent', 'chromicacid', 'h2cro4'], formulas: ['kmno4', 'mno4k', 'cro3'] },
  nabh4: { names: ['nabh4', 'sodiumborohydride', 'borohydride'], formulas: ['nabh4', 'bh4na'] },
  lialh4: { names: ['lialh4', 'lithiumaluminiumhydride', 'lithiumaluminumhydride', 'lah'], formulas: ['lialh4', 'alh4li'] },
  h2so4_dehydrate: { smiles: 'OS(=O)(=O)O', names: ['h2so4', 'sulfuricacid', 'sulphuricacid', 'dehydration'], formulas: ['h2so4', 'h2o4s'] },
  hydroboration: { names: ['bh3', 'borane', 'hydroboration', 'bh3thf', '9bbn', 'bh3thenh2o2/naoh', 'bh3thenh2o2naoh', 'bh3/thfthenh2o2/naoh', 'bh3h2o2', 'bh3/h2o2', 'bh3thf/h2o2', 'boranethenperoxide', 'hydroborationoxidation', '1bh32h2o2naoh', '1bh3thf2h2o2naoh', 'bh3thenh2o2', 'b2h6thenh2o2/naoh'], formulas: ['bh3', 'b2h6'] },
  pcc:   { names: ['pcc', 'pyridiniumchlorochromate', 'collinsreagent'], formulas: [] },
  hbr_peroxide: { names: ['hbrperoxide', 'hbrroor', 'hbr+roor', 'hbr+peroxide', 'hbr+peroxides', 'hbrperoxides'], formulas: [] },
};
let _ruleKeyCache = null;
/* =====================================================================
   READING A REAGENT AS A CHEMIST WRITES IT

   A reagent over an arrow is not an entry in a name table. It is a phrase:
   "3 H2 / Pd-C", "mCPBA (1 equiv)", "Zn(Hg) + HCl", "H2O2 / CH3COOH",
   "CH3COCl / AlCl3", "(1) O3; (2) Me2S". Looking those up by name fails, and
   failing silently is worse than being wrong, because the student is told
   "not a reagent I know" about the most standard reagent in the chapter.

   Two things fix it. First this PHRASE table, which matches the way pairs are
   written and is tried BEFORE the name table — because "H2O2 / AcOH" means one
   specific reaction and neither half of it means that reaction on its own, and
   "CH3COCl / AlCl3" is a Friedel-Crafts acylation, not an acyl chloride.
   Second `reagentForms` below, which peels a phrase back to the plain name:
   drops a trailing parenthetical, a leading stoichiometric coefficient, a
   named solvent, and finally tries each half of an X / Y pair.
   ===================================================================== */
const REAGENT_PHRASE = [
  /* hydrogenation: any metal, any stoichiometry */
  { re: /^\s*\d*\s*h2\s*[\/,+]\s*(?:pd|pt|ni|raney|rh|ru|lindlar|c\b|pd-?c|pd\/c)/i, id: 'h2' },
  { re: /^\s*h2\s*[\/,+]\s*pd\s*[\/,-]?\s*c?\s*,?\s*(?:quinoline|lindlar)/i, id: 'lindlar' },
  /* aromatic halogenation needs the Lewis acid; the alkene version does not */
  { re: /^\s*br2\s*[\/,+]\s*(?:febr3|fecl3|alcl3|fe\b)/i, id: 'ar_br2' },
  { re: /^\s*cl2\s*[\/,+]\s*(?:fecl3|alcl3|fe\b)/i, id: 'ar_cl2' },
  { re: /^\s*br2\s*[\/,+]\s*(?:ccl4|ch2cl2|dcm|cs2|h2o)/i, id: 'br2' },
  { re: /^\s*(?:i2|br2)\s*[\/,+]\s*(?:h2so4|hno3)/i, id: 'ar_br2' },
  /* nitration and sulfonation */
  { re: /^\s*hno3\s*[\/,+]\s*h2so4/i, id: 'nitration' },
  { re: /^\s*(?:so3|h2so4)\s*[\/,+]\s*(?:h2so4|so3)/i, id: 'sulfonation' },
  /* Friedel-Crafts and any other Lewis-acid- or acid-promoted reaction: the catalyst is what
     makes it one, and the OTHER half is the molecule that reacts. Read that half as the
     reagent molecule — whatever the "use common reagents" setting says, because the phrase
     is explicit about what is happening. */
  { re: /^(.+?)\s*[\/,+]\s*(?:alcl3|fecl3|febr3|sncl4|ticl4|bf3|bf3\s*\u00b7?\s*(?:et2o|oet2)|hf|zncl2|h3po4)\s*(?:\([^)]*\))?\s*$/i, moleculeFrom: 1 },
  { re: /^\s*alcl3\s*$/i, id: null },
  /* peracid from hydrogen peroxide and acetic acid — and the plain peracids */
  { re: /^\s*h2o2\s*[\/,+]\s*(?:ch3co2h|ch3cooh|acoh|aceticacid|h2o|naoh)/i, id: 'mcpba' },
  { re: /^\s*(?:mcpba|m-?cpba|peraceticacid|peroxyaceticacid|ch3co-?o-?oh|ch3coooh|mmpp|magnesiummonoperoxyphthalate|peroxybenzoicacid|rcooh\b)/i, id: 'mcpba' },
  /* the C=O to CH2 reductions */
  { re: /^\s*(?:nh2nh2|n2h4|hydrazine)\s*[\/,+;]?\s*(?:koh|naoh|kot-?bu|base|heat|Δ|,)?/i, id: 'wolffkishner' },
  { re: /^\s*zn\s*\(?\s*hg\s*\)?\s*[\/,+]?\s*(?:hcl|h2so4|conc)?/i, id: 'clemmensen' },
  /* Swern, however it is written */
  { re: /(?:dmso|cs\(c\)=o)\s*[\/,+;]\s*\(?\s*cocl\)?2|oxalylchloride/i, id: 'swern' },
  /* ozonolysis: reductive workup gives the carbonyls, oxidative gives the acid */
  { re: /^\s*o3\b(?!.*(?:h2o2|kmno4))/i, id: 'ozonolysis' },
  { re: /^\s*o3\b.*h2o2/i, id: 'ozonolysis_ox' },
  /* alkyne hydration, both regiochemistries */
  { re: /hgso4|hg\(oac\)2\s*[\/,+]\s*h2so4/i, id: 'alkynehyd' },
  { re: /sia2bh|disiamylborane|\(sia\)2bh|9-?bbn/i, id: 'alkynehydb' },
  /* aqueous acid and aqueous base as the general workup-or-hydrolysis reagent */
  { re: /^\s*(?:h3o\+?|h\+\s*[\/,+]\s*h2o|h2o\s*[\/,+]\s*h(?:\+|2so4|cl)|aq\.?\s*acid|dilute\s*acid|h2o\s*,?\s*h\+)/i, id: 'aqacid' },
  { re: /^\s*h2o\s*(?:\(\s*(?:large\s*)?excess\s*\))?\s*,?\s*h\+/i, id: 'aqacid' },
  { re: /^\s*(?:-?oh|ho-|hydroxide|oh-|aq\.?\s*base|naoh\s*[\/,+]\s*h2o)/i, id: 'naoh' },
  /* the metals that make an organometallic, and the bases written as formulas */
  { re: /^\s*(?:mg|magnesium)\s*(?:[\/,+]\s*(?:et2o|ether|thf))?\s*$/i, id: null, molecule: '[Mg]' },
  { re: /^\s*(?:li|lithium)\s*(?:metal)?\s*(?:[\/,+]\s*(?:et2o|ether|pentane))?\s*$/i, id: null, molecule: '[Li]' },
  { re: /^\s*(?:naoch3|naome|sodiummethoxide|ch3ona|ch3o-?)\s*$/i, id: null, molecule: 'C[O-].[Na+]' },
  { re: /^\s*(?:naoch2ch3|naoet|sodiumethoxide)\s*$/i, id: 'naoet' },
  { re: /^\s*(?:nasch3|sodiumthiomethoxide|ch3s-?)\s*$/i, id: null, molecule: 'C[S-].[Na+]' },
  { re: /^\s*(?:nh3|ammonia)\s*$/i, id: null, molecule: 'N' },
  { re: /^\s*(?:ch3i|mei|methyliodide|iodomethane)\s*$/i, id: null, molecule: 'CI' },
  { re: /^\s*(?:ch3br|mebr|bromomethane|methylbromide)\s*$/i, id: null, molecule: 'CBr' },
  { re: /^\s*(?:buli|n-?buli|butyllithium|n-?butyllithium)\s*$/i, id: null, molecule: 'CCCC[Li]' },
  { re: /^\s*h\+?\s*(?:or|\/|,)\s*-?\s*oh\b/i, id: 'tautomer' },
  { re: /^\s*(?:na\+?\s*)?-?\s*o(?:h|-)\s*(?:\(|,|$)/i, id: 'naoh' },
  { re: /^\s*(?:na\+?\s*)?-?\s*och2ch3|^\s*(?:na\+?\s*)?-?\s*oc2h5|^\s*(?:na\+?\s*)?-?\s*oet\b/i, id: 'naoet' },
  { re: /^\s*(?:na\+?\s*)?-?\s*hco3|nahco3|sodiumbicarbonate|na2co3|k2co3\s*$/i, id: 'naoh' },
  { re: /^\s*(?:hcn|-?\s*cn)\b/i, id: 'nacn' },
  { re: /^\s*(?:cl2|br2|i2)\s*,?\s*(?:excess\s*)?[\/,+]?\s*-?\s*oh\b/i, id: 'br2' },
  { re: /^\s*socl2|thionylchloride|oxalylchloride/i, id: 'soCl2' },
  { re: /^\s*pcl[35]/i, id: 'pcl5' },
  { re: /^\s*ch2n2|diazomethane/i, id: 'ch2n2' },
  { re: /tempo/i, id: 'tempo' },
  { re: /nabh3cn|nacnbh3|na\s*bh\s*\(\s*oac\s*\)\s*3|cyanoborohydride|triacetoxyborohydride/i, id: 'nabh3cn' },
  { re: /^\s*(?:dcc|edc|dic|hatu|hbtu|pybop|t3p|cdi)\b/i, id: null, coupling: true },
  { re: /^\s*\(?\s*(?:ch3|me)\s*\)?2\s*culi/i, id: null, molecule: 'C[Cu]C.[Li]' },
  { re: /^\s*(?:ph|c6h5)2\s*culi/i, id: null, molecule: '[Li].c1ccccc1[Cu]c1ccccc1' },
  { re: /^\s*\(?\s*(?:bu|ch3ch2ch2ch2)\s*\)?2\s*culi/i, id: null, molecule: 'CCCC[Cu]CCCC.[Li]' },
  { re: /^\s*\(?\s*(?:ch2=ch|vinyl)\s*\)?2\s*culi/i, id: null, molecule: 'C=C[Cu]C=C.[Li]' },
  { re: /^\s*(?:ag2o|silveroxide|silver\(i\)oxide)/i, id: 'ag2o' },
  { re: /^\s*(?:tbaf|bu4nf|\(ch3ch2ch2ch2\)4n\+?f-?|tetrabutylammoniumfluoride)/i, id: 'tbs_off' },
  { re: /^\s*(?:tbdms-?cl|tbscl|tbdmscl|\(ch3\)3csi\(ch3\)2cl|t-?bu\s*me2\s*sicl)/i, id: 'tbs_on' },
  { re: /^\s*(?:hg\(oac\)2|mercuric\s*acetate|mercury\(ii\)\s*acetate)/i, id: 'alkoxyhg' },
  { re: /^\s*(?:lialh\(o-?t-?bu\)3|lithiumtri-?tert-?butoxyaluminiumhydride)/i, id: 'dibal' },
  { re: /^\s*(?:nbs|n-?bromosuccinimide)\s*[\/,+]?\s*(?:hv|hν|light|roor|peroxide|bpo|aibn|heat|ccl4)?/i, id: 'allylbr' },
  { re: /^\s*br2\s*[\/,+]\s*(?:hv|hν|light|heat|Δ)/i, id: 'allylbr' },
  { re: /^\s*(?:kmno4|potassiumpermanganate|na2cr2o7|k2cr2o7|cro3)\s*(?:[\/,+]\s*(?:h2so4|h3o\+|h2o|naoh|koh|heat))?\s*$/i, id: 'kmno4' },
  { re: /^\s*(?:tsoh|p-?tsoh|h2so4|hcl|camphorsulfonicacid|csa)\s*\(?\s*cat/i, id: 'acidcat' },
  { re: /^\s*(?:ticl4|tiCl4)\s*[\/,+]/i, id: 'ziegler' },
  { re: /^\s*(?:aibn|bpo|roor|benzoylperoxide|radicalinitiator)\s*(?:,?\s*heat)?\s*$/i, id: 'radpoly' },
];

/* every progressively plainer reading of one written reagent */
function reagentForms(t0) {
  const out = [];
  const push = x => { x = String(x || '').trim().replace(/[,;]+$/, ''); if (x && out.indexOf(x) < 0) out.push(x); };
  push(t0);
  /* a trailing parenthetical is a gloss or a stoichiometry, not part of the name:
     "mCPBA (m-chloroperoxybenzoic acid)", "HBr (1 equiv)", "NaOH (or any base)" */
  push(String(t0).replace(/\s*\([^()]*\)\s*$/, ''));
  push(String(t0).replace(/\s*\([^()]*\)/g, ' '));
  /* a leading coefficient: "3 H2 / Pd-C", "2 CH3MgBr", "2 R'Li" */
  push(String(t0).replace(/^\s*\d+(?:\.\d+)?\s*(?:equiv\.?|eq\.?|x)?\s+/i, ''));
  /* a named solvent or temperature tacked on with a comma is handled by the splitter,
     but "mCPBA, CH2Cl2, 0 C" reaching here as one piece should still work */
  push(String(t0).replace(/\s*,\s*(?:ch2cl2|dcm|thf|et2o|ether|hexane|toluene|dmf|dmso|meoh|etoh|h2o|water|ccl4|dioxane|benzene)\b.*$/i, ''));
  /* each half of a pair, longest first — "CH3COCl / AlCl3" already matched a phrase above,
     so what reaches here is a pair where one half IS the reagent ("KI / acetone") */
  const parts = String(t0).split(/\s*[\/+]\s*/).map(x => x.trim()).filter(Boolean);
  if (parts.length > 1) parts.slice().sort((a, b) => b.length - a.length).forEach(push);
  return out;
}

/* the phrase table, applied to one written reagent */
function reagentPhrase(text) {
  const raw = String(text || '').trim();
  if (!raw) return null;
  const flat = raw.replace(/\s+/g, '');
  for (const p of REAGENT_PHRASE) {
    if (!(p.re.test(raw) || p.re.test(flat))) continue;
    if (p.id) { if (Chem.RULES.some(r => r.id === p.id)) return { id: p.id, via: 'phrase' }; continue; }
    if (p.molecule) {
      try { const rule = moleculeRule(Chem.parseSmiles(p.molecule), raw); if (rule) return { id: rule.id, via: 'phrase — read as the molecule ' + p.molecule }; } catch (e) {}
    }
    if (p.coupling) return { cond: true, via: 'a coupling agent — it activates the acid; the amine or alcohol is the other reactant' };
    if (p.moleculeFrom) {
      const m = p.re.exec(raw) || p.re.exec(flat);
      const half = m && m[p.moleculeFrom] ? m[p.moleculeFrom].trim() : '';
      if (half) {
        /* the half may itself be a taught reagent name ("CH3COCl") or a structure */
        try { const r = resolvePredReagent(half, true); if (r && r.id) return { id: r.id, via: 'phrase — the reacting half, \u201c' + half + '\u201d' }; } catch (e) {}
        try { const g = Chem.searchMolecule(half).graph; const rule = moleculeRule(g, half); if (rule) return { id: rule.id, via: 'phrase — \u201c' + half + '\u201d as the molecule, with the Lewis acid as the catalyst' }; } catch (e) {}
      }
    }
  }
  return null;
}
function resolvePredReagent(text, _noForms) {
  const raw = String(text || '').trim();
  if (!raw) return null;
  { const ct = casReagentText(raw);
    if (ct) { const r = resolvePredReagent(ct); if (r && r.id) return { id: r.id, via: 'CAS ' + raw }; } }
  /* a written PHRASE first: a pair of reagents usually means one named reaction */
  { const ph = reagentPhrase(raw); if (ph && ph.id) return ph; if (ph && ph.cond) return null; }
  const norm = raw.toLowerCase().replace(/[\s\u00b7.()\[\]\-]/g, '');
  /* One written reagent can mean several different reactions now: "H₂ / Pd–C" reduces an alkene,
     takes a Cbz off, takes a benzyl ether off, or reduces a nitro group, depending on what is in
     the flask. So collect every rule the name matches and prefer one that actually has a site on
     the drawn molecule, rather than taking whichever happens to be first in the table. */
  {
    const hits = [];
    for (const id in RULE_REAGENTS) {
      const r = RULE_REAGENTS[id];
      if (r.names.includes(norm)) { hits.push(id); continue; }
      const rule = Chem.RULES.find(x => x.id === id);
      if (rule && rule.label.toLowerCase().replace(/[\s\u00b7.()\[\]\-]/g, '') === norm) hits.push(id);
    }
    if (hits.length === 1) return { id: hits[0], via: 'name' };
    if (hits.length > 1) {
      const sk = window.REAGENT_SUBSTRATE;
      if (sk && sk.atoms && sk.atoms.length) {
        const fits = hits.find(id => { const rl = Chem.RULES.find(x => x.id === id); try { return rl && ruleApplies(rl, sk); } catch (e) { return false; } });
        if (fits) return { id: fits, via: 'name — the one that fits what is drawn' };
      }
      /* nothing drawn fits any of them: prefer the plainest reading of the name */
      const plain = ['h2', 'nabh4', 'lialh4'].find(x => hits.indexOf(x) >= 0);
      return { id: plain || hits[0], via: 'name' };
    }
  }
  /* organic-chemistry acronyms (Reich 2002): mapped to a taught reagent
     where that is chemically honest, otherwise a definition instead of a guess */
  {
    const ac = acronymLookup(raw);
    if (ac) {
      if (ac.rule) {
        /* with “common reagents” off, an acronym that has a real structure (LDA, NaH …) whose taught
           stand-in does not fit the drawn molecule is used as the molecule itself instead */
        if (window._prefs && window._prefs.commonReagents === false && ac.s) {
          try {
            const rule0 = Chem.RULES.find(x => x.id === ac.rule); const sk = window.REAGENT_SUBSTRATE;
            if (rule0 && sk && sk.atoms && sk.atoms.length && !ruleApplies(rule0, sk)) { const rule = moleculeRule(Chem.parseSmiles(ac.s), raw); if (rule) return { id: rule.id, via: 'molecule' }; }
          } catch (e) {}
        }
        return { id: ac.rule, via: 'acronym', acro: ac };
      }
      return { error: raw.toUpperCase() + ' = ' + ac.f + ' — ' + ac.r +
        '. It is not one of the taught reagents, so no product is predicted for it.' };
    }
  }
  for (const id in RULE_REAGENTS) {
    if (RULE_REAGENTS[id].formulas.includes(norm)) return { id, via: 'formula' };
  }
  /* structure: SMILES or IUPAC name — spectator Na+/K+/Li+ stripped,
     then matched by canonical structure against the taught reagents */
  let g = null;
  try { g = Chem.searchMolecule(raw).graph; } catch (e) { return { error: e.message }; }
  const strip = gr => {
    const keep = gr.atoms.filter(a => !['Na', 'K', 'Li'].includes(a.element));
    if (keep.length === gr.atoms.length || !keep.length) return gr;
    const ids = new Set(keep.map(a => a.id));
    return { atoms: keep, bonds: gr.bonds.filter(b => ids.has(b.a) && ids.has(b.b)), nextId: gr.nextId };
  };
  if (!_ruleKeyCache) {
    _ruleKeyCache = {};
    for (const id in RULE_REAGENTS) {
      const s = RULE_REAGENTS[id].smiles;
      if (!s) continue;
      try { _ruleKeyCache[id] = Chem.canonicalKey(strip(Chem.parseSmiles(s))); } catch (e) {}
    }
  }
  let key = null;
  try { key = Chem.canonicalKey(strip(g)); } catch (e) { return { error: 'Could not read that as a reagent.' }; }
  for (const id in _ruleKeyCache) if (_ruleKeyCache[id] === key) return { id, via: 'structure' };
  /* any other molecule: with "common reagents" switched off it becomes a reagent of its own */
  if (window._prefs && window._prefs.commonReagents === false) {
    const rule = moleculeRule(g, raw);
    if (rule) return { id: rule.id, via: 'molecule' };
  }
  /* Before giving up: try the plainer readings of the same phrase — without its trailing
     parenthetical, without a leading coefficient, one half of a pair. This is what makes
     "HBr (1 equiv)" and "3 H2 / Pd-C" work without a table entry for every variant. */
  if (!_noForms) {
    for (const form of reagentForms(raw)) {
      if (form === raw) continue;
      let r = null;
      try { r = resolvePredReagent(form, true); } catch (e) { r = null; }
      if (r && r.id) return { id: r.id, via: 'read as “' + form + '”' };
    }
  }
  if (!_noForms) { const mr = moleculeRule(g, raw); if (mr) return { id: mr.id, via: 'molecule' }; }
  return { error: 'Read the structure fine, but it is not one of the taught reagents — the buttons below list what the predictor knows.' };
}
/* A whole molecule as a reagent: a dynamic rule the engine's molecule+molecule
   family handles. Registered once per structure, so it shows up as a chip. */
function moleculeRule(g, typed) {
  let key = null; try { key = Chem.libKeyHash(Chem.canonicalKey(g)); } catch (e) { return null; }
  const id = 'mol:' + key;
  let rule = Chem.RULES.find(r => r.id === id);
  if (rule) return rule;
  let nm = ''; try { nm = Chem.displayNameFor(g) || ''; } catch (e) {}
  if (!nm) { try { nm = Chem.formula(g); } catch (e) { nm = typed || 'molecule'; } }
  rule = { id, mechClass: 'molecule', label: nm, needs: 'any', molecule: g, typed: typed || nm,
    apply() { return { error: 'This molecule needs a partner with a matching functional group.' }; } };
  Chem.RULES.push(rule);
  return rule;
}

const CAS_REAGENT_TEXT = {
  '10035-10-6': 'HBr', '7647-01-0': 'HCl', '10034-85-2': 'HI', '7782-50-5': 'Cl2', '7726-95-6': 'Br2',
  '937-14-4': 'mCPBA', '20816-12-0': 'OsO4', '7664-93-9': 'H2SO4', '7789-60-8': 'PBr3', '143-33-9': 'NaCN',
  '16853-85-3': 'LiAlH4', '14044-65-6': 'BH3', '7722-84-1': 'H2O2', '1333-74-0': 'H2', '7440-05-3': 'H2',
  '1310-73-2': 'NaOH', '141-52-6': 'NaOEt', '865-47-4': 't-BuOK', '26299-14-9': 'PCC',
  '7722-64-7': 'KMnO4', '16940-66-2': 'NaBH4', '7732-18-5': 'H2O',
};
function casReagentText(raw) {
  const t = String(raw || '').trim();
  try { if (Chem.looksLikeCAS(t) && CAS_REAGENT_TEXT[t]) return CAS_REAGENT_TEXT[t]; } catch (e) {}
  return null;
}

function hasBenzylicCH(g) {
  try {
    return g.atoms.some(a => a.element === 'C' && !Chem.isAromaticCarbon(g, a.id) && Chem.implicitH(g, a) >= 2 &&
      Chem.neighbors(g, a.id).every(n => n.bond.order === 1 && n.atom.element === 'C') && Chem.neighbors(g, a.id).some(n => Chem.isAromaticCarbon(g, n.atom.id)));
  } catch (e) { return false; }
}
function engineFamily(id) { try { return RXN.familyOf[id] || null; } catch (e) { return null; } }
function engineHasSite(id, g) {
  try {
    const fam = RXN.familyOf[id]; if (!fam) return false;
    const env = RXN.resolveEnv(id, {});
    const live = f => { try { return (f.sites(g, id, env) || []).some(x => !x.inert); } catch (e) { return false; } };
    if (live(fam)) return true;
    /* a reagent's SECOND job counts too: BH3 reduces an amide when there is no alkene,
       aqueous HBr hydrolyses an ester when there is no C=C */
    return [].concat(fam.fallback || []).some(k => RXN.families[k] && live(RXN.families[k]));
  } catch (e) { return false; }
}
function ruleApplies(rule, g) {
  const gs = Chem.findGroups(g);
  const have = new Set(gs.map(x => x.type));
  /* the structural families (protecting groups, cross-coupling, nitro reduction) know their own
     substrates far better than a functional-group census does — ask them */
  const famS = engineFamily(rule.id);
  if (famS && famS.structural) return engineHasSite(rule.id, g);
  if (famS && famS.fallback && engineHasSite(rule.id, g)) return true;
  if (rule.needs === 'carbonyl-or-nitrile') return have.has('carbonyl') || have.has('nitrile') || engineHasSite(rule.id, g);
  if (rule.needs === 'benzylic') return hasBenzylicCH(g);
  if (rule.needs === 'alkene-or-alkyne') return have.has('alkene') || have.has('alkyne');
  if (rule.needs === 'arene') return gs.some(x => x.type === 'alkene' && x.aromatic);
  return have.has(rule.needs);
}
const NEED_NAMES = { 'arylhalide': 'an aryl or vinyl halide (C–X on an sp² carbon), or an aryl triflate',
                     'nitro': 'a nitro group (–NO₂)',
                     'protectable': 'a group this reagent can protect — a free N–H amine, an O–H alcohol, or a ketone / aldehyde',
                     'protected': 'the protecting group this reagent removes',
                     'arene': 'an aromatic ring with a free C–H', 'alkene': 'an alkene (C=C)', 'alkyne': 'an alkyne', 'alkene-or-alkyne': 'an alkene or alkyne',
                     'halide': 'an alkyl halide (C–Br/Cl/I/F)', 'alcohol': 'an alcohol (C–OH)', 'carbonyl': 'a carbonyl (C=O)',
                     'carbonyl-or-nitrile': 'a carbonyl (C=O) or a nitrile (C≡N)', 'amide': 'an amide (C=O with the nitrogen on it)',
                     'amine': 'an amine (C–NH₂, C–NHR or C–NR₂)', 'diazonium': 'an aryl diazonium salt (Ar–N₂⁺) — make one first with NaNO₂ / HCl at 0–5 °C',
                     'epoxide': 'an epoxide (a three-membered ring with an oxygen in it)',
                     'ether': 'an ether (C–O–C, not an ester)',
                     'alkene-or-arene': 'an allylic or benzylic C–H — an sp³ carbon next to a C=C or an aromatic ring',
                     'any': 'nothing in particular; this reagent looks at whatever is there', 'benzylic': 'a benzylic CH₃ / CH₂ (alkyl carbon on an aromatic ring)' };

const ACRONYMS = {
  'ac': { f: 'acetyl (CH3C=O)', r: 'substituent/protecting group — acetyl on O or N' },
  'acoh': { f: 'acetic acid', r: 'weak acid, solvent', s: 'CC(=O)O' },
  'ac2o': { f: 'acetic anhydride', r: 'acylating agent — makes esters and amides', s: 'CC(=O)OC(C)=O' },
  'aibn': { f: 'azobis(isobutyronitrile)', r: 'radical initiator — splits on heating to start radical chains' },
  '9-bbn': { f: '9-borabicyclo[3.3.1]nonane', r: 'hindered hydroboration reagent — very regioselective anti-Markovnikov', rule: 'hydroboration' },
  'bht': { f: 'butylated hydroxytoluene (2,6-di-tert-butyl-4-methylphenol)', r: 'radical inhibitor / antioxidant', s: 'Cc1cc(C(C)(C)C)c(O)c(C(C)(C)C)c1' },
  'binap': { f: '2,2\'-bis(diphenylphosphino)-1,1\'-binaphthyl', r: 'chiral ligand for asymmetric catalysis' },
  'bms': { f: 'borane dimethyl sulfide (BH3·SMe2)', r: 'stable borane source for hydroboration', rule: 'hydroboration' },
  'bn': { f: 'benzyl (PhCH2-)', r: 'substituent/protecting group' },
  'boc': { f: 'tert-butyloxycarbonyl', r: 'amine protecting group (removed with acid)' },
  'bz': { f: 'benzoyl (PhC=O) — caution: sometimes used for benzyl', r: 'substituent/protecting group' },
  'cbz': { f: 'carbobenzyloxy (BnO-C=O)', r: 'amine protecting group (removed by hydrogenolysis)' },
  'can': { f: 'ceric ammonium nitrate', r: 'one-electron oxidant' },
  'csa': { f: 'camphorsulfonic acid', r: 'chiral organic acid catalyst' },
  'dabco': { f: '1,4-diazabicyclo[2.2.2]octane', r: 'small, unhindered amine base/nucleophilic catalyst', s: 'C1CN2CCN1CC2' },
  'dast': { f: '(diethylamino)sulfur trifluoride', r: 'converts OH and C=O to fluorides' },
  'dbn': { f: '1,5-diazabicyclo[4.3.0]non-5-ene', r: 'strong, hindered amidine base — elimination without substitution', s: 'C1CN=C2CCCN2C1', rule: 'tbuok' },
  'dbu': { f: '1,8-diazabicyclo[5.4.0]undec-7-ene', r: 'strong, hindered amidine base — the classic E2-without-SN2 base', s: 'C1CCC2=NCCCN2CC1', rule: 'tbuok' },
  'dcc': { f: 'N,N-dicyclohexylcarbodiimide', r: 'coupling agent — activates acids for amide/ester formation', s: 'C(=NC1CCCCC1)=NC1CCCCC1' },
  'ddq': { f: '2,3-dichloro-5,6-dicyano-1,4-benzoquinone', r: 'strong organic oxidant', s: 'N#CC1=C(C#N)C(=O)C(Cl)=C(Cl)C1=O' },
  'dead': { f: 'diethyl azodicarboxylate', r: 'Mitsunobu partner with PPh3', s: 'CCOC(=O)N=NC(=O)OCC' },
  'dhp': { f: '3,4-dihydro-2H-pyran', r: 'makes THP ethers (alcohol protection)', s: 'C1=COCCC1' },
  'dibal': { f: 'diisobutylaluminium hydride (also DIBAL-H, DIBAH)', r: 'bulky partial reducing agent — ester to aldehyde at low T; not in the taught set (use LiAlH4/NaBH4 rules)' },
  'diglyme': { f: 'diethylene glycol dimethyl ether', r: 'high-boiling ether solvent', s: 'COCCOCCOC' },
  'dmac': { f: 'N,N-dimethylacetamide', r: 'polar aprotic solvent', s: 'CC(=O)N(C)C' },
  'dmad': { f: 'dimethyl acetylenedicarboxylate', r: 'reactive dienophile/Michael acceptor', s: 'COC(=O)C#CC(=O)OC' },
  'dmap': { f: '4-(dimethylamino)pyridine', r: 'acylation catalyst — supercharged pyridine', s: 'CN(C)c1ccncc1' },
  'dme': { f: '1,2-dimethoxyethane (glyme)', r: 'ether solvent', s: 'COCCOC' },
  'dmf': { f: 'N,N-dimethylformamide', r: 'polar aprotic solvent — great for SN2', s: 'CN(C)C=O' },
  'dmpu': { f: 'N,N\'-dimethylpropyleneurea', r: 'polar aprotic cosolvent (HMPA substitute)' },
  'dmso': { f: 'dimethyl sulfoxide', r: 'polar aprotic solvent — great for SN2; Swern oxidant', s: 'CS(C)=O' },
  'da': { f: 'Diels-Alder reaction', r: 'named reaction — diene + dienophile cycloaddition' },
  'edta': { f: 'ethylenediaminetetraacetic acid', r: 'metal chelator' },
  'ee': { f: 'enantiomeric excess', r: 'measure of optical purity: %major − %minor' },
  'de': { f: 'diastereomeric excess', r: 'same idea as ee, for diastereomers' },
  'et3n': { f: 'triethylamine (TEA)', r: 'common mild amine base', s: 'CCN(CC)CC' },
  'fc': { f: 'Friedel-Crafts reaction', r: 'named reaction — aromatic alkylation/acylation with AlCl3' },
  'fmoc': { f: '9-fluorenylmethoxycarbonyl', r: 'amine protecting group (removed with base)' },
  'glyme': { f: '1,2-dimethoxyethane (= DME)', r: 'ether solvent', s: 'COCCOC' },
  'hmpa': { f: 'hexamethylphosphoric triamide', r: 'strongly cation-solvating cosolvent (carcinogen)' },
  'hsab': { f: 'hard-soft acid-base principle', r: 'hard prefers hard, soft prefers soft' },
  'hvz': { f: 'Hell-Volhard-Zelinsky reaction', r: 'named reaction — alpha-bromination of carboxylic acids' },
  'kda': { f: 'potassium diisopropylamide', r: 'very strong hindered base', rule: 'tbuok' },
  'lah': { f: 'lithium aluminium hydride (LiAlH4)', r: 'strong hydride reducing agent', rule: 'lialh4' },
  'lda': { f: 'lithium diisopropylamide', r: 'THE strong, hindered, non-nucleophilic base — deprotonates, eliminates, never substitutes', rule: 'tbuok', s: 'CC(C)[N-]C(C)C.[Li+]' },
  'lhmds': { f: 'lithium hexamethyldisilazide (LiN(SiMe3)2)', r: 'strong hindered silazide base', rule: 'tbuok' },
  'litmp': { f: 'lithium tetramethylpiperidide', r: 'very hindered strong base', rule: 'tbuok' },
  'lta': { f: 'lead tetraacetate', r: 'oxidant — oxidative cleavages' },
  'ma': { f: 'maleic anhydride', r: 'classic dienophile', s: 'O=C1C=CC(=O)O1' },
  'mcpba': { f: 'meta-chloroperoxybenzoic acid (MCPBA)', r: 'peracid — epoxidation, Baeyer-Villiger', s: 'O=C(OO)c1cccc(Cl)c1', rule: 'mcpba' },
  'mem': { f: '2-methoxyethoxymethyl', r: 'alcohol protecting group' },
  'mes': { f: 'mesityl (2,4,6-trimethylphenyl)', r: 'bulky aryl substituent' },
  'mom': { f: 'methoxymethyl (CH3OCH2-)', r: 'alcohol protecting group' },
  'ms': { f: 'methanesulfonyl (mesyl, CH3SO2-)', r: 'makes OH into a mesylate leaving group' },
  'mscl': { f: 'methanesulfonyl chloride', r: 'reagent that installs the mesylate leaving group', s: 'CS(=O)(=O)Cl' },
  'mvk': { f: 'methyl vinyl ketone', r: 'Michael acceptor', s: 'C=CC(C)=O' },
  'nbs': { f: 'N-bromosuccinimide', r: 'low-concentration Br source — allylic/benzylic radical bromination (not the taught Br2 addition)', s: 'O=C1CCC(=O)N1Br' },
  'ncs': { f: 'N-chlorosuccinimide', r: 'low-concentration Cl source', s: 'O=C1CCC(=O)N1Cl' },
  'nis': { f: 'N-iodosuccinimide', r: 'electrophilic I source', s: 'O=C1CCC(=O)N1I' },
  'nmo': { f: 'N-methylmorpholine N-oxide', r: 'co-oxidant that recycles OsO4 in dihydroxylation', rule: 'oso4' },
  'nmp': { f: 'N-methylpyrrolidone', r: 'polar aprotic solvent', s: 'CN1CCCC1=O' },
  'pcc': { f: 'pyridinium chlorochromate', r: 'mild oxidant — alcohol to aldehyde/ketone, stops there', rule: 'pcc' },
  'cro2cl2': { f: 'chromyl chloride (Étard reagent)', r: 'oxidises a benzylic CH₃ to the aldehyde via the Étard complex — no over-oxidation', rule: 'etard' },
  'pdc': { f: 'pyridinium dichromate', r: 'chromium oxidant, slightly stronger than PCC', rule: 'pcc' },
  'pmb': { f: 'para-methoxybenzyl', r: 'alcohol protecting group' },
  'ppa': { f: 'polyphosphoric acid', r: 'strong dehydrating acid medium' },
  'ppts': { f: 'pyridinium p-toluenesulfonate', r: 'mild acid catalyst' },
  'ptc': { f: 'phase-transfer catalyst', r: 'shuttles ions between water and organic phases' },
  'py': { f: 'pyridine', r: 'base/solvent/catalyst', s: 'c1ccncc1' },
  'rt': { f: 'room temperature', r: '~20-25 °C, 293-298 K' },
  'tbaf': { f: 'tetra-n-butylammonium fluoride', r: 'fluoride source — removes silyl protecting groups' },
  'tbhp': { f: 'tert-butyl hydroperoxide', r: 'oxidant (Sharpless epoxidation partner)', s: 'CC(C)(C)OO' },
  'tbs': { f: 'tert-butyldimethylsilyl (= TBDMS)', r: 'robust alcohol protecting group' },
  'tbdps': { f: 'tert-butyldiphenylsilyl', r: 'bulkier silyl protecting group' },
  'tea': { f: 'triethylamine (= Et3N)', r: 'common mild amine base', s: 'CCN(CC)CC' },
  'tes': { f: 'triethylsilyl', r: 'silyl protecting group' },
  'tf': { f: 'triflyl (CF3SO2-)', r: 'triflate — among the best leaving groups known' },
  'tfa': { f: 'trifluoroacetic acid', r: 'strong organic acid — removes Boc', s: 'O=C(O)C(F)(F)F' },
  'thf': { f: 'tetrahydrofuran', r: 'the everyday ether solvent', s: 'C1CCOC1' },
  'thp': { f: 'tetrahydropyranyl (from DHP)', r: 'alcohol protecting group' },
  'tips': { f: 'triisopropylsilyl', r: 'bulky silyl protecting group' },
  'tmeda': { f: 'N,N,N\',N\'-tetramethylethylenediamine', r: 'chelating amine — activates organolithiums', s: 'CN(C)CCN(C)C' },
  'tms': { f: 'trimethylsilyl (also tetramethylsilane, the NMR zero)', r: 'smallest silyl protecting group' },
  'tsoh': { f: 'p-toluenesulfonic acid (tosic acid)', r: 'strong organic acid catalyst — dehydrations', s: 'Cc1ccc(S(=O)(=O)O)cc1', rule: 'h2so4_dehydrate' },
  'tscl': { f: 'p-toluenesulfonyl chloride', r: 'installs the tosylate leaving group on alcohols', s: 'Cc1ccc(S(=O)(=O)Cl)cc1' },
  'ts': { f: 'tosyl (p-CH3C6H4SO2-)', r: 'tosylate — turns OH into a great leaving group' },
  'trityl': { f: 'triphenylmethyl (Tr)', r: 'very bulky protecting group' },
  'wk': { f: 'Wolff-Kishner reduction', r: 'named reaction — C=O to CH2 with hydrazine/base' },
};
function acronymLookup(text) {
  const k = String(text || '').trim().toLowerCase().replace(/[\s.]/g, '');
  return ACRONYMS[k] ? { key: k, ...ACRONYMS[k] } : null;
}

const PRED_COND = { cyclise: 'dilute, base or Δ', dmp: 'CH₂Cl₂, rt', swern: '−78 °C; then Et₃N', tempo: 'NaOCl, KBr, 0 °C', jones: 'acetone/H₂O',
  soCl2: 'reflux, neat or in toluene', pcl5: 'rt', ch2n2: 'ether, 0 °C',
  radpoly: '60–150 °C, pressure', cationic: 'low temperature', anionic: 'dry, aprotic', ziegler: 'heterogeneous, mild',
  nabh3cn: 'MeOH, pH 6–7', tautomer: 'a trace of acid or base', lda: 'THF, −78 °C', enolate_base: 'EtOH, rt', thermal: '100–150 °C',
  ar_cl2: 'FeCl₃ (cat.), 0–25 °C', ar_i2: 'HNO₃ as the oxidant', epox_h2o: 'dilute H₂SO₄, H₂O',
  epox_hx: 'cold, inert solvent', epox_roh: 'the alcohol as solvent, H₂SO₄ cat.',
  clemmensen: 'Zn(Hg), conc. HCl, reflux', wolffkishner: 'then KOH, diethylene glycol, 180 °C',
  ozonolysis: 'CH₂Cl₂, −78 °C, then Me₂S', ozonolysis_ox: 'CH₂Cl₂, −78 °C, then H₂O₂',
  alkynehyd: 'H₂SO₄ / HgSO₄, H₂O, 60 °C', alkynehydb: 'then H₂O₂ / NaOH',
  allylbr: 'AIBN or hν, CCl₄, reflux', aqacid: 'dilute, rt (heat for a nitrile)',
  ag2o: 'THF / H₂O, rt', alkoxyhg: 'then NaBH₄',
  nano2: 'aq. HCl, 0–5 °C', sandmeyer_cl: 'CuCl, aq. HCl, 0 °C → rt', sandmeyer_br: 'CuBr, aq. HBr, 0 °C → rt',
  sandmeyer_cn: 'CuCN / KCN, pH 7', diazo_i: 'KI (aq.), 0 °C → rt', diazo_f: 'HBF₄, then heat dry 100–140 °C',
  diazo_oh: 'H₂O / H₂SO₄, 50–100 °C', diazo_h: 'H₃PO₂, 0 °C → rt',
  mno2: 'CH₂Cl₂, rt', dibal: 'toluene, −78 °C', lindlar: 'Pb(OAc)₂ / quinoline', na_nh3: 'NH₃ (l), −78 °C',
  nitrored: 'MeOH, 1 atm', boc_on: 'Et₃N, CH₂Cl₂', boc_off: '25 % in CH₂Cl₂', cbz_on: 'NaHCO₃ (aq.)', cbz_off: 'MeOH, 1 atm',
  fmoc_on: 'NaHCO₃ (aq.)', fmoc_off: '20 % in DMF', tbs_on: 'imidazole, DMF', tbs_off: 'THF, rt', bn_on: 'NaH, DMF',
  bn_off: 'MeOH, 1 atm', thp_on: 'cat. TsOH', ac_on: 'pyridine', ms_on: 'Et₃N, 0 °C', ts_on: 'pyridine, 0 °C',
  acetal_on: 'cat. TsOH, Dean–Stark', acetal_off: 'acetone/H₂O',
  suzuki_ph: 'K₂CO₃, 80 °C', buchwald_morph: 'Cs₂CO₃, 100 °C', sonogashira_ph: 'Et₃N, 60 °C',
  naocl: 'AcOH, 0\u201325 \u00b0C', etard: 'CS\u2082, 0 \u00b0C; then H\u2082O', sulfonation: 'fuming H\u2082SO\u2084', nitration: 'H\u2082SO\u2084, < 50 \u00b0C', ar_br2: 'FeBr\u2083 (cat.)', h2so4_dehydrate: '\u0394 (heat)', naoet: 'EtOH, \u0394', tbuok: 'tBuOH, \u0394',
  hbr_peroxide: 'ROOR, \u0394 or h\u03bd', h2: 'Pd/C', hydroboration: 'then H\u2082O\u2082, NaOH',
  kmno4: 'H\u2083O\u207a workup', pcc: 'CH\u2082Cl\u2082', mcpba: 'CH\u2082Cl\u2082', oso4: 'NMO',
  nacn: 'DMSO', lialh4: 'Et\u2082O; then H\u2083O\u207a', nabh4: 'MeOH', cl2: 'CH\u2082Cl\u2082',
  br2: 'CH\u2082Cl\u2082', pbr3: 'Et\u2082O, 0 \u00b0C', naoh: 'H\u2082O', hydration: 'H\u2082SO\u2084 (cat.)' };

window.Reagents = { resolve: resolvePredReagent, moleculeRule, ruleApplies, engineHasSite, reagentForms,
                    RULE_REAGENTS, ACRONYMS, PRED_COND, NEED_NAMES };
})();
