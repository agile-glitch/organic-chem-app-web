/* MolInfo — what the Molecule tab says about a molecule, about each atom and about each bond. Pure computation, no
   DOM. Every label, note and caveat the views show is written here, so the 2D view, the 3D view, the hover card and
   the side panel all say the same thing.

   Inputs: window.MolData (element and bond tables, js/mol_data.js), an RDKit.js module R and an RDKit molecule M
   (hydrogens suppressed; isotopic H stay as atoms), window.Chem (optional, for names).

   Atom numbering (shared by every view): atoms 0..n-1 are M's atoms; then ONE explicit H per implicit H, appended
   grouped by parent in ascending parent index (the order RDKit add_hs uses). Bonds 0..m-1 are M's bonds, then one X–H
   bond per appended H in the same order. Display numbers are 1-based ("C1" is atom index 0).

   Where the facts come from (all checked; see tools/test_molinfo.js):
   - Conjugation, hybridization (RDKit's), Gasteiger–Marsili charges and the strict rotatable-bond rule are ports of
     RDKit 2026.03.6 (ConjugHybrid.cpp, GasteigerCharges.cpp/GasteigerParams.cpp, Lipinski.cpp strict SMARTS). The
     port agreed with RDKit Python on 21,869 molecules (hybridization, conjugation and rotatable flags: 0 mismatches;
     Gasteiger: largest difference 5e-11; tools/test_molinfo.js --bulk repeats the comparison).
   - Lone pairs, octet, steric number and VSEPR geometry (Gillespie's VSEPR rules), the displayed hybridization (from
     the steric number, sp² for a lone pair conjugated with a π system), oxidation state (IUPAC 2016 rule: bonding
     electrons to the more electronegative atom, here with Pauling values; aromatic bonds averaged over the Kekulé
     structures), σ/π counts, orbital overlap, polarity class (ΔEN cut-offs, TEXT.polarityCutoffs), E/Z candidates,
     dihedral names (IUPAC Gold Book, Klyne & Prelog 1960): 489 hand-checked textbook values in the fixtures.
   - Molecule-level statistics (RDKit descriptors; CIAAW weights), the chirality rule (section 10) and the
     functional-group pattern table (section 10): checked against RDKit Python on 18 reference molecules
     (612 fields), 41 chirality cases and 521+ functional-group tests.
   - Atomic weights (CIAAW 2024), electronegativity (Pauling, Allen), radii, bond enthalpies and typical bond lengths
     come from MolData with their own citations.
   - Dipole: no size or direction (a Gasteiger estimate was measured and failed, see DIPOLE); the verdict "zero by
     symmetry / zero in this conformation only / polar" comes from the conformer's symmetry (section 12b).
   - After the review of the Molecule tab: names only when the app's parser reads them back (section 15); tautomeric
     P/S oxo centres are not stereocentres; typical lengths never borrowed from another group's row; resonance bonds get
     no single average enthalpy; textbook rotation separate from RDKit's rotatable-bond count.
   - After the second review round: six-membered aromatic length rows never used for a bond in a five-membered ring;
     phosphine / phosphine-oxide / phosphate-triester rows from P's own neighbours; acetal, hemiacetal and epoxide C–O
     shown with the ether entry and carbonate / carbamate alkyl–O with the ester's (each saying so); ring fusion next
     to a biaryl axis counts as an ortho group; ring cis/trans at fusion atoms (decalin) and the one-configuration
     bridgeheads of small bridged rings (norbornene); flexible rings never "zero by symmetry"; axial chirality through
     spiro atoms and exocyclic C=C; oxidation states of delocalised ions marked "as drawn".
   - After the third review round: molecule-level 3D values withheld when the force-field check flags a heavy-atom
     bond, and angles / dihedrals next to a flagged bond marked; "linear" / "planar" before the NPR shape class;
     phenyl rows for aromatic bonds shared with a non-aromatic ring; primary / tertiary amide, ortho-substituted
     biaryl rows; no vinyl row for an acyl halide C–X or a vinyl ester's O–C; stand-in notes name the entry shown;
     "as drawn" oxidation states for charge-separated resonance structures (1,3-dipoles).

   API, window.MolInfo:
     analyse(R, M, opts?) → record            (throws on an empty molecule; M is not deleted — the caller owns it)
     geometry(record, coords, opts?) → geo    coords: Float64Array(3N) in record order. opts: optimised (false → no
         measured numbers at all, except bonds in flaggedBonds [[a, b, text]], which keep their length and the text),
         radii (per-atom override, for tests), mcPerAtom, areaPoints. geo.molecule.vdwVolume / vdwArea are
         computed on first read; the whole result is cached for the same coords array and options.
         geo.molecule.dipole (computed on first read) = {verdict, text, note, symmetry, torsions?, ring?}: verdict 'zero' |
         'zero-conformer' | 'polar' | 'not assessed' | 'ion'; text is the value to show, note its caveat (TEXT.dipole);
         symmetry lists the operations found ('i', 'σ', 'C2' …); ring: true when a flexible ring (not a torsion) makes
         the zero conformer-only. No size and no direction.
         geo.molecule.form = 'linear' | 'planar' | '3d' (every atom within 0.05 Å of a line / 0.1 Å of a plane; formDev
         = {line, plane} largest distances, Å); geo.molecule.shape stays the NPR class ('rod' | 'disc' | 'sphere').
         Flagged bonds (optimised model): geo.flagged (Map bond → text), geo.flaggedAtoms (Set of their atoms);
         angles[c][t].unreliable / dihedral[k].unreliable mark values that involve such an atom. When a bond between two
         heavy atoms is flagged, geo.molecule.withheld = {bonds, names, text, note} and every molecule-level value
         (rg, pmi, npr1, npr2, shape, form, extent, vdwVolume, vdwArea, centre) is null and the dipole verdict is 'not
         assessed' (a flagged model's shape can be wrong: the app's MMFF94 model of SO₃ is pyramidal, O–S–O 98°).
     record.stereo = {atoms: [{index, label, kind: 'centre'|'ring', ring?: 'substituent'|'fusion'|'bridgehead'}], bonds,
         ignored: [{index, reason}]}: ignored atoms are centres RDKit proposes that are not stereocentres (two or more H
         on P; tautomeric P/S oxo; the bridgeheads of a small bridged ring such as norbornene, which have one possible
         configuration); a view listing open stereo from another source (the 3D engine) should drop them too.
     record.bonds[k].rotation = {free, text} (textbook), .rotatable = {yes, reason, code} (RDKit strict descriptor),
         .conjugated (textbook: carbon radicals count) and .rdkitConjugated, .resonance (no single average enthalpy).
     atomRows(record, i, geo|null), bondRows(record, k, geo|null) → [{key, label, value, kind, note?, main?}]
         main: true marks the ~10 rows for a compact hover card.
     moleculeSections(record, geo|null, extra?) → [{title, rows}]
         extra = {names: names(), conformer: {energy, relEnergy, source}, component: {index, count, formula}, combined}
     labelsFor(record, mode) → string[N]  mode 'element'|'index'|'cip'|'hybridization'|'formal'|'partial'|'oxidation'|'lonepairs'
     bondLabelsFor(record, 'cip'|'order') → string[bonds]   (E / Z / ? on double bonds)
     names(R, M) → {iupac, iupacNote, iupacChecked, common, cas}   via window.Chem, each checked to describe exactly
         this structure (the systematic name: see NAMES_RULE and namerRefusal)
     combined(R, fullMol) → {formula, formulaText, componentsText, mw, charge, components:[{formula, formulaText, charge,
         heavyAtoms, atoms}]}   for salts and mixtures (components in atom order, as get_frags gives them)
     TEXT (shared caveat wording), KINDS, DIPOLE (measured figures), FG_TABLE, util {angle, dihedral, klynePrelog, …}
   Row kinds: 'exact' (follows from the structure), 'model' (a model's estimate), 'table' (a tabulated average or
   constant), '3d' (measured on the current 3D model), 'convention' (a definition or bookkeeping rule).
   Every RDKit.js object created here is deleted before returning.                                                   */
(() => {
  'use strict';

  /* ================================================================ 0. constants and element access */
  // RDKit 2026.03.6 periodic table (atomic_data.cpp): outer-shell electrons and default valence (-1 = any), index Z.
  // Needed so the conjugation / hybridization / Gasteiger ports reproduce RDKit exactly (checked in the test).
  const RD_NOUTER = [0, 1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2,
    3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
    13, 14, 15, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
  const RD_DV = [0, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, 3, 4, 3, 2, 1, 0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 3, 4, 3, 2, 1, 0, 1, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 3, 2, 3, 2, 1, 0, 1, 2, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 2, 3, 2, 1, 0, 1, 2, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1];
  const NOUT = z => RD_NOUTER[z] || 0;
  const DV = z => (RD_DV[z] === undefined ? -1 : RD_DV[z]);

  const MD = () => {
    if (!window.MolData) throw new Error('MolInfo needs js/mol_data.js (window.MolData)');
    return window.MolData;
  };
  const elCache = new Map();
  function EL(z) {                                      // element facts, with a stub for the dummy atom '*'
    if (elCache.has(z)) return elCache.get(z);
    const e = z > 0 ? MD().element(z) : null;
    const v = e || { Z: z, sym: z === 0 ? '*' : '?', name: z === 0 ? 'unspecified atom' : 'unknown', en: null, aw: null,
      rcov: null, rvdw: null, src: {} };
    elCache.set(z, v);
    return v;
  }
  const SYM = z => EL(z).sym;
  const EN = z => (z > 0 ? EL(z).en : null);            // Pauling (Allred 1961), MolData
  const ALLEN = z => (z > 0 && EL(z).enAllen != null ? EL(z).enAllen : null);

  const KINDS = ['exact', 'model', 'table', '3d', 'convention'];
  const MINUS = '−';
  const signed = (v, dp) => {                            // +0.210 / −0.397 / 0 (typographic minus)
    const r = Number(v.toFixed(dp));
    if (r === 0) return (0).toFixed(dp);
    return (r > 0 ? '+' : MINUS) + Math.abs(r).toFixed(dp);
  };
  const num = (v, dp) => { const s = Math.abs(v).toFixed(dp); return v < 0 && Number(s) !== 0 ? MINUS + s : s; };   // never "−0.0"
  const trimNum = (v, dp) => { const s = String(Number(v.toFixed(dp))); return s.startsWith('-') ? MINUS + s.slice(1) : s; };
  const intSigned = v => (v > 0 ? '+' + v : v < 0 ? MINUS + (-v) : '0');
  const HYB_TEXT = { sp: 'sp', sp2: 'sp²', sp3: 'sp³', sp3d: 'sp³d', sp3d2: 'sp³d²' };
  const hybText = h => (h ? HYB_TEXT[h] || h : '');
  const SUB = '₀₁₂₃₄₅₆₇₈₉', SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹';
  const subDigits = s => s.replace(/\d/g, d => SUB[+d]);
  const supDigits = s => s.replace(/\d/g, d => SUP[+d]);
  const chargeSup = c => (c === 0 ? '' : (Math.abs(c) > 1 ? supDigits(String(Math.abs(c))) : '') + (c > 0 ? '⁺' : '⁻'));
  const chargeLabel = c => (c === 0 ? '' : (Math.abs(c) > 1 ? Math.abs(c) : '') + (c > 0 ? '+' : MINUS));
  const atomName = (record, i) => record.atoms[i].el + (i + 1);
  const plural = (n, one, many) => n + ' ' + (n === 1 ? one : many || one + 's');
  const WORDS = ['no', 'one', 'two', 'three', 'four', 'five', 'six'];
  const countWord = n => WORDS[n] || String(n);

  /* Wording shared by every view (accuracy caveats). */
  const TEXT = {
    asDrawn: 'for the Lewis structure as drawn (an S=O and an S⁺–O⁻ drawing of the same molecule give different values)',
    osAsDrawn: 'for the Lewis structure as drawn: the formal charge is delocalised by resonance, and a drawing with the ' +
      'charge on another atom can give a different value here',
    gasteiger: 'Gasteiger–Marsili partial charge (empirical model)',
    gasteigerShort: 'Gasteiger–Marsili (empirical model)',
    gasteigerNote: 'An empirical charge model (PEOE, RDKit parameters, 12 iterations), not a measured or quantum-chemical ' +
      'charge. The formal charge is not the partial charge: N in NH₄⁺ comes out −0.37.',
    gasteigerSO: 'Unreliable for S=O / P=O compounds: the value depends on how the structure is drawn (DMSO sulfur is ' +
      '+0.015 drawn as S=O and −0.060 drawn as S⁺–O⁻).',
    geometry3d: 'force-field (MMFF94) model geometry of this conformer, not an experimental structure',
    notOptimised: '— (geometry not force-field optimised)',
    notOptimisedBanner: 'approximate geometry (not force-field optimised): no measured lengths, angles, volume or polarity are shown',
    // a force-field model with a bond the MMFF94 parameter-gap check flagged (CO₂, SO₂, SO₃, N₂O)
    unreliableHere: '(model geometry not reliable here)',
    unreliableNote: 'a value that involves an atom of a bond the force-field check flagged (MMFF94 has no specific ' +
      'parameters for that bond) is not reliable',
    flaggedValues: names => `— (model geometry not reliable: flagged bond${names.length === 1 ? '' : 's'} ${names.join(', ')})`,
    flaggedValuesNote: 'The force-field check flagged a bond of this model as not reliable (MMFF94 has no specific ' +
      'parameters for it), so the whole shape can be wrong (the MMFF94 model of SO₃ comes out pyramidal instead of ' +
      'planar). No molecule-level 3D values are shown: longest extent, radius of gyration, principal moments, shape, ' +
      'van der Waals volume and area, and the dipole verdict. Bond angles and dihedrals that involve an atom of the ' +
      'flagged bond are marked "model geometry not reliable here".',
    mmffEnergy: 'MMFF94 steric energy (kcal/mol) — only comparable between conformers of this molecule; not a heat of ' +
      'formation; no entropy or solvent',
    esp: 'Coulomb potential of Gasteiger point charges on the van der Waals surface — a model, not a quantum-chemical ESP',
    namer: "name from the app's namer, checked: the app's name parser reads it back to exactly this structure (stereo " +
      'included). It may differ from the preferred IUPAC name in style',
    namerUnchecked: "from the app's namer (not independently checked)",
    aromaticModel: "RDKit's aromaticity model (it also counts 2-pyridone, uracil and caffeine's six-membered ring as aromatic)",
    sssr: 'rings of the smallest set of smallest rings (SSSR)',
    cip: "CIP rules (RDKit's new CIP labeller)",
    rotation: 'textbook description: rotation about a single bond is free at room temperature (a barrier of a few kJ/mol, ' +
      'e.g. about 12 kJ/mol in ethane); rotation is restricted about double and triple bonds, ring bonds, the partial ' +
      'double bonds of amides and esters, and the axis of a hindered biaryl (atropisomers)',
    rotatable: "RDKit's strict rotatable-bond rule, a drug-likeness descriptor that counts bonds whose rotation changes the " +
      'heavy-atom shape (it does not say whether a bond can turn): bonds to terminal atoms, CH₃/CF₃/CCl₃/CBr₃/tert-butyl ' +
      'ends, bonds next to a triple bond, ring bonds and amide/ester-type C(=O)–N/O/S bonds are not counted; unlike RDKit, ' +
      'formyl C(=O)–N/O bonds (formamide, formates) are not counted either',
    conjugated: 'conjugated = part of an unbroken chain of p orbitals: a π bond next to another π bond, a lone pair ' +
      '(N, O, C⁻), an empty p orbital (C⁺) or the half-filled p orbital of a carbon radical (benzyl, allyl). Otherwise ' +
      "RDKit's definition: halogens and S or P with two or more neighbours never count, boron's empty p orbital is not " +
      'counted, and a CH₃ is never conjugated (no hyperconjugation). RDKit itself does not count radical centres.',
    polarityCutoffs: 'ΔEN (Pauling) < 0.4 nonpolar covalent (OpenStax Chemistry 2e §7.2); 0.4 to < 2.0 polar covalent; ' +
      '≥ 2.0 largely ionic (OpenStax Organic Chemistry §2.1); a bond from Na, K, Rb, Cs, Ca, Sr or Ba to a non-metal is ' +
      'largely ionic',
    enthalpy: 'average over many compounds; the real value depends on the molecule',
    length: 'mean over organic crystal structures (Allen et al. 1987); one bond in one molecule can differ by a few ' +
      'hundredths of an ångström',
    dihedral: 'for this force-field conformer; other conformers have other values',
    functionalGroups: null,      // set after the pattern table (section 10), which gives the count
    lipinski: 'Lipinski counts: donors = H on N or O; acceptors = N + O atoms',
    tpsa: "TPSA (Ertl 2000), counting N and O only; PubChem's TPSA also counts S and P, so it is larger for sulfur and " +
      'phosphorus compounds',
    logp: 'Crippen logP (estimate; fit σ = 0.68 log units on 9,920 training molecules, Wildman & Crippen 1999)',
    mr: 'Crippen molar refractivity (estimate)',
    ro5: 'Lipinski rule of five: pass with at most one violation (MW > 500, logP > 5, H donors > 5, H acceptors > 10); ' +
      'logP by Crippen (Lipinski used CLogP)',
    veber: 'Veber: rotatable bonds ≤ 10 and TPSA ≤ 140 Å²; rotatable bonds by RDKit strict, not Veber\'s own count',
    volume: 'union of atomic spheres with the van der Waals radii listed per atom (H 1.10 Å); seeded Monte Carlo, ' +
      'statistical error < 1%',
    area: 'van der Waals surface (probe radius 0), golden-spiral points on each sphere; within 1% of a 20,000-point reference',
    shape: 'linear: every atom within 0.05 Å of one straight line; planar: every atom within 0.1 Å of one plane; ' +
      'otherwise the nearest corner of the normalised principal-moment ratios (NPR, Sauer & Schwarz 2003): rod (0,1), ' +
      'disc (0.5,0.5), sphere (1,1). A planar molecule always lies on the rod–disc edge of that triangle (NPR1 + NPR2 = 1), ' +
      'so rod / disc / sphere is used only for molecules that are not flat',
  };

  /* DIPOLE. No size and no direction is shown. A dipole estimated from Gasteiger charges on a force-field conformer
     was measured against the 19 experimental gas-phase dipoles of the reference table (NIST CCCBDB experimental list;
     RDKit ETKDGv3 seed 42, lowest of 30 MMFF94 conformers, RDKit Gasteiger charges = this file's): mean absolute error
     0.59 D, median signed error −13 %, Spearman rank correlation 0.78, and about 4× too small for esters and pyrroles
     (methyl acetate 0.29 D against 1.72 D). So the size is not shown, and a size-based class ("nonpolar / weakly polar /
     polar") is not either: it called esters and pyrrole "weakly polar". What is shown instead is what symmetry decides
     exactly (section 12b): whether the conformer's symmetry operations leave any direction unchanged. If none, the
     dipole is zero; otherwise symmetry allows one (its size can still be tiny: propane 0.08 D). Validation of the
     verdicts: DIPOLE.CHECK (tests/molview fixtures, RDKit MMFF94 and the app's OpenChemLib conformers). */
  const DIPOLE = { MAE: 0.59, SPEARMAN: 0.78, MEDIAN_PCT: -13, N_EXP: 19, TOL: 0.15,
    CHECK: { molecules: 92, conformers: 229, agree: 229 } };      // tests/molview fixtures; tools/test_molinfo.js recounts
  TEXT.dipoleVerdict = {
    zero: 'zero by symmetry (nonpolar)',
    zeroConformer: 'zero in this conformation only (other conformations can be polar)',
    zeroRing: 'zero in this conformation; other ring conformations may not cancel',
    polar: "polar: the molecule's symmetry allows a dipole moment (size not estimated)",
    polarConformer: 'polar in this conformation: its symmetry allows a dipole moment (size not estimated)',
    notAssessed: 'not assessed: too many atoms for the symmetry search',
    hypervalent: 'not assessed: MMFF94 has no parameters for an atom with five or six electron domains, so the 3D shape ' +
      'around it (and with it the symmetry) is not reliable',
    ion: 'not defined for an ion (it depends on the choice of origin)',
    flagged: names => `not assessed: the model geometry is not reliable (flagged bond${names.length === 1 ? '' : 's'} ${names.join(', ')})`,
  };
  TEXT.dipole = 'Decided from the symmetry of this 3D conformer: a molecule has no dipole moment exactly when its symmetry ' +
    'operations leave no direction unchanged (for example a centre of inversion, or two different rotation axes). Atoms ' +
    'are matched within ' + DIPOLE.TOL + ' Å, and CH₃-type end groups count as freely rotating. Symmetry decides only whether ' +
    'a dipole can exist, not its size: alkanes have tiny dipoles (propane 0.08 D). "Zero in this conformation only" means ' +
    'that turning a single bond can remove the symmetry (hydrogen peroxide, hydrazine, butane); "zero in this ' +
    'conformation" for a molecule with a flexible ring means that another ring conformation can lack it (the chair of ' +
    'cyclohexane-1,4-dione has no dipole, but less symmetric ring conformations can have one), so zero by symmetry is ' +
    'never claimed for it. No size or direction is shown: a dipole estimated from Gasteiger charges was off by ' +
    DIPOLE.MAE + ' D on average against ' + DIPOLE.N_EXP +
    ' experimental dipoles (NIST CCCBDB). Checked on ' + DIPOLE.CHECK.molecules + ' molecules with known dipoles (' +
    DIPOLE.CHECK.conformers + ' force-field conformers, RDKit and OpenChemLib): the verdict agreed with experiment for ' +
    (DIPOLE.CHECK.agree === DIPOLE.CHECK.conformers ? 'all of them' : DIPOLE.CHECK.agree) + '; no molecule with a dipole was called zero.';

  /* ================================================================ 1. graph from RDKit.js get_json */
  function fromRDKitJson(J) {
    const d = J.defaults, m = J.molecules[0];
    const ext = (m.extensions || []).find(e => e.name === 'rdkitRepresentation') || {};
    const aA = new Set(ext.aromaticAtoms || []), aB = new Set(ext.aromaticBonds || []);
    const pick = (o, k, dflt) => (o[k] === undefined ? dflt[k] : o[k]);
    const atoms = m.atoms.map((a, i) => ({
      z: pick(a, 'z', d.atom), hs: pick(a, 'impHs', d.atom), chg: pick(a, 'chg', d.atom), rad: pick(a, 'nRad', d.atom),
      iso: pick(a, 'isotope', d.atom), stereo: pick(a, 'stereo', d.atom), arom: aA.has(i) }));
    const bonds = (m.bonds || []).map((b, k) => ({ a: b.atoms[0], b: b.atoms[1], order: pick(b, 'bo', d.bond), arom: aB.has(k) }));
    const g = prepare({ atoms, bonds, rings: ext.atomRings || [] });
    g.cipRanks = ext.cipRanks || null;
    return g;
  }
  function prepare(g) {
    g.adj = g.atoms.map(() => []);
    g.bonds.forEach((b, k) => { g.adj[b.a].push({ nb: b.b, bond: k }); g.adj[b.b].push({ nb: b.a, bond: k }); });
    return g;
  }
  const degree = (g, i) => g.adj[i].length;                          // RDKit getDegree
  const totalDegree = (g, i) => g.adj[i].length + g.atoms[i].hs;     // RDKit getTotalDegree
  const nH = (g, i) => g.atoms[i].hs + g.adj[i].filter(e => g.atoms[e.nb].z === 1).length;
  const sumBO = (g, i) => g.adj[i].reduce((s, e) => s + g.bonds[e.bond].order, 0) + g.atoms[i].hs;   // Kekulé valence
  const period = z => (z <= 2 ? 1 : z <= 10 ? 2 : z <= 18 ? 3 : z <= 36 ? 4 : z <= 54 ? 5 : z <= 86 ? 6 : 7);
  function group(z) {                                   // IUPAC group 1-18, null for the f block
    if (z === 1) return 1; if (z === 2) return 18;
    const p = period(z), first = [0, 1, 3, 11, 19, 37, 55, 87][p];
    let k = z - first;
    if (p <= 3) return k < 2 ? k + 1 : k + 11;
    if (p >= 6) { if (k >= 3 && k < 17) return null; if (k >= 17) k -= 14; }
    return k + 1;
  }

  /* ================================================================ 2. RDKit conjugation (== GetIsConjugated) */
  function countAtomElec(g, i) {
    const a = g.atoms[i], dv = DV(a.z);
    if (dv <= 1) return -1;                              // univalent or "any valence" elements
    const deg = degree(g, i) + a.hs;
    if (deg > 3) return -1;
    const nlp = Math.max(NOUT(a.z) - dv - a.chg, 0);
    return (dv - deg) + nlp - a.rad;                     // only the sign (> 0) matters
  }
  function isConjCand(g, i) {
    const a = g.atoms[i], dv = DV(a.z);
    if (!a.chg && dv >= 0 && sumBO(g, i) > dv) return false;     // neutral hypervalent (sulfone S, phosphate P)
    const no = NOUT(a.z);
    return (a.z <= 10 || (no !== 5 && no !== 6) || (no === 6 && totalDegree(g, i) < 2)) && countAtomElec(g, i) > 0;
  }
  const piLike = b => b.arom || b.order >= 2;
  // extraCand(i) (optional): more atoms that count as conjugation candidates. Without it this is exactly RDKit's rule
  // (used for hybridization and charges); the displayed flag adds carbon radical centres (textbookConjugation).
  function rdkitConjugation(g, extraCand) {
    const n = g.atoms.length, conj = g.bonds.map(b => !!b.arom), cand = [], nsub = [];
    for (let i = 0; i < n; i++) { cand[i] = isConjCand(g, i) || !!(extraCand && extraCand(i)); nsub[i] = cand[i] ? degree(g, i) + g.atoms[i].hs : 0; }
    for (let i = 0; i < n; i++) {
      if (!cand[i] || nsub[i] < 2 || nsub[i] > 3) continue;
      for (const e1 of g.adj[i]) {
        if (!piLike(g.bonds[e1.bond]) || !cand[e1.nb]) continue;
        for (const e2 of g.adj[i]) {
          if (e2.bond === e1.bond || nsub[e2.nb] > 3) continue;
          if (cand[e2.nb]) { conj[e1.bond] = true; conj[e2.bond] = true; }
        }
      }
    }
    return conj;
  }
  // a carbon radical with three σ bonds and no lone pair (planar, its half-filled p orbital can overlap a π bond)
  const carbonRadicalCentre = (g, i) => { const a = g.atoms[i];
    return a.z === 6 && a.rad === 1 && totalDegree(g, i) === 3 && sumBO(g, i) === 3 && a.chg === 0; };
  const textbookConjugation = g => rdkitConjugation(g, i => carbonRadicalCentre(g, i));

  /* ================================================================ 3. RDKit hybridization (== GetHybridization) */
  function rdkitHybridization(g, conj) {
    return g.atoms.map((a, i) => {
      if (a.z === 0) return 'UNSPECIFIED';
      const td = totalDegree(g, i);
      if ((a.stereo === 'cw' || a.stereo === 'ccw') && td === 4) return 'SP3';
      let norbs;
      if (a.z < 89) {
        if (a.z <= 1) norbs = td;
        else {
          const tv = sumBO(g, i), nout = NOUT(a.z), free = nout - (tv + a.chg);
          if (tv + nout - a.chg < 8) norbs = td + Math.trunc((free - a.rad) / 2) + a.rad;
          else norbs = td + Math.trunc(free / 2);
        }
      } else norbs = td;
      switch (norbs) {
        case 0: case 1: return 'S';
        case 2: return 'SP';
        case 3: return 'SP2';
        case 4: return (td > 3 || !g.adj[i].some(e => conj[e.bond])) ? 'SP3' : 'SP2';
        case 5: return 'SP3D';
        case 6: return 'SP3D2';
        default: return 'UNSPECIFIED';
      }
    });
  }

  /* ================================================================ 4. Gasteiger–Marsili (== RDKit ComputeGasteigerCharges, 12 iterations) */
  // [a, b, c]: chi(q) = a + b q + c q²; chi+ = a + b + c (H: fixed 20.02). RDKit GasteigerParams.cpp.
  const GPARAM = {
    'H *': [7.17, 6.24, -0.56],
    'C sp3': [7.98, 9.18, 1.88], 'C sp2': [8.79, 9.32, 1.51], 'C sp': [10.39, 9.45, 0.73],
    'N sp3': [11.54, 10.82, 1.36], 'N sp2': [12.87, 11.15, 0.85], 'N sp': [15.68, 11.7, -0.27],
    'O sp3': [14.18, 12.92, 1.39], 'O sp2': [17.07, 13.79, 0.47],
    'F sp3': [14.66, 13.85, 2.31], 'Cl sp3': [11.00, 9.69, 1.35], 'Br sp3': [10.08, 8.47, 1.16], 'I sp3': [9.9, 7.96, 0.96],
    'S sp3': [10.14, 9.13, 1.38], 'S so': [10.14, 9.13, 1.38], 'S so2': [12.00, 10.81, 1.20], 'S sp2': [10.88, 9.49, 1.33],
    'P sp3': [8.90, 8.24, 0.96], 'P sp2': [9.665, 8.530, 0.735],
    'Si sp3': [7.300, 6.567, 0.657], 'Si sp2': [7.905, 6.748, 0.443], 'Si sp': [9.065, 7.027, -0.002],
    'B sp3': [5.980, 6.820, 1.605], 'B sp2': [6.420, 6.807, 1.322],
    'Be sp3': [3.845, 6.755, 3.165], 'Be sp2': [4.005, 6.725, 3.035],
    'Mg sp2': [3.565, 5.572, 2.197], 'Mg sp3': [3.300, 5.587, 2.447], 'Mg sp': [4.040, 5.472, 1.823],
    'Al sp3': [5.375, 4.953, 0.867], 'Al sp2': [5.795, 5.020, 0.695],
  };
  const IONXH = 20.02, DAMP = 0.5, DAMP_SCALE = 0.5, EPS = 2.220446049250313e-16;

  function gasteiger(g, conj, hyb, nIter) {
    nIter = nIter || 12;
    const n = g.atoms.length, q = new Array(n).fill(0), hq = new Array(n).fill(0), missing = [];
    for (let i = 0; i < n; i++) {                       // formal charges spread over same-element atoms 2 conjugated bonds away
      let formal = g.atoms[i].chg;
      if (Math.abs(formal) > EPS && Math.abs(q[i]) < EPS) {
        const marker = [i];
        for (const e1 of g.adj[i]) {
          if (!conj[e1.bond]) continue;
          for (const e2 of g.adj[e1.nb]) {
            if (e2.bond === e1.bond || !conj[e2.bond]) continue;
            if (g.atoms[e2.nb].z === g.atoms[i].z) { formal += g.atoms[e2.nb].chg; marker.push(e2.nb); }
          }
        }
        for (const k of marker) q[k] = formal / marker.length;
      }
    }
    const P = [], ionX = [];
    for (let i = 0; i < n; i++) {
      const a = g.atoms[i], sym = a.z ? SYM(a.z) : '*';
      let mode = { SP3: 'sp3', SP2: 'sp2', SP: 'sp' }[hyb[i]];
      if (!mode) {
        if (a.z === 1) mode = '*';
        else if (a.z === 16) {
          const no = g.adj[i].filter(e => g.atoms[e.nb].z === 8).length;
          mode = no === 2 ? 'so2' : no === 1 ? 'so' : 'sp3';
        } else mode = '';
      }
      P[i] = GPARAM[sym + ' ' + mode];
      if (!P[i]) {                                      // RDKit silently uses zeros: NaN/Inf wherever they meet a bond.
        P[i] = [0, 0, 0];                               // An isolated ion keeps its formal charge (as RDKit gives).
        if (g.adj[i].length || a.hs) missing.push(sym);
      }
      ionX[i] = a.z === 1 ? IONXH : P[i][0] + P[i][1] + P[i][2];
    }
    const hP = GPARAM['H *'], energ = new Array(n).fill(0);
    let damp = DAMP;
    for (let it = 0; it < nIter; it++) {
      for (let i = 0; i < n; i++) energ[i] = P[i][0] + q[i] * (P[i][1] + P[i][2] * q[i]);
      for (let i = 0; i < n; i++) {
        let dq = 0;
        for (const e of g.adj[i]) {
          const dx = energ[e.nb] - energ[i], sgn = dx < 0 ? 0 : 1;
          dq += dx / (sgn * (ionX[i] - ionX[e.nb]) + ionX[e.nb]);
        }
        const niHs = g.atoms[i].hs;
        if (niHs > 0) {
          const qHs = hq[i] / niHs, enr = hP[0] + qHs * (hP[1] + hP[2] * qHs), dx = enr - energ[i], sgn = dx < 0 ? 0 : 1;
          const dqH = dx / (sgn * (ionX[i] - IONXH) + IONXH);
          dq += niHs * dqH;
          hq[i] -= niHs * dqH * damp;
        }
        q[i] += damp * dq;
      }
      damp *= DAMP_SCALE;
    }
    const ok = !missing.length && q.every(Number.isFinite) && hq.every(Number.isFinite);
    return { q, hq, ok, missing: [...new Set(missing)] };
  }

  /* ================================================================ 5. textbook atom values (Lewis, VSEPR, hybridization) */
  const LEWIS = new Set([1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 32, 33, 34, 35, 36, 51, 52, 53, 54]);
  function valenceElectrons(z) {
    const gr = group(z);
    if (gr === 1 || gr === 2) return gr;
    if (gr >= 13) return z === 2 ? 2 : gr - 10;
    return null;
  }
  function lewis(g, i) {
    const a = g.atoms[i], V = valenceElectrons(a.z);
    if (V === null) return null;
    const B = sumBO(g, i), nonbond = V - a.chg - B, lp2 = nonbond - a.rad;
    const lp = lp2 >= 0 && lp2 % 2 === 0 ? lp2 / 2 : null;
    return { V, bondOrderSum: B, nonbonding: nonbond, lp, radicals: a.rad, shell: 2 * B + nonbond };
  }
  const VSEPR = {
    '2,0': ['linear', 'linear', '180°'],
    '3,0': ['trigonal planar', 'trigonal planar', '120°'], '3,1': ['trigonal planar', 'bent', 'slightly less than 120°'],
    '4,0': ['tetrahedral', 'tetrahedral', '109.5°'], '4,1': ['tetrahedral', 'trigonal pyramidal', 'slightly less than 109.5°'],
    '4,2': ['tetrahedral', 'bent', 'less than 109.5°'],
    '5,0': ['trigonal bipyramidal', 'trigonal bipyramidal', '90°, 120° and 180°'], '5,1': ['trigonal bipyramidal', 'seesaw', 'about 90° and 120°'],
    '5,2': ['trigonal bipyramidal', 'T-shaped', 'about 90°'], '5,3': ['trigonal bipyramidal', 'linear', '180°'],
    '6,0': ['octahedral', 'octahedral', '90° and 180°'], '6,1': ['octahedral', 'square pyramidal', 'about 90°'],
    '6,2': ['octahedral', 'square planar', '90°'],
  };
  const EDG_ONLY = { 2: 'linear', 3: 'trigonal planar', 4: 'tetrahedral', 5: 'trigonal bipyramidal', 6: 'octahedral' };
  function isCarbocation(g, j) {
    const a = g.atoms[j]; if (a.z !== 6 || a.chg !== 1 || a.rad) return false;
    return totalDegree(g, j) === 3 && sumBO(g, j) === 3;
  }
  function isAmideN(g, i) {                              // N single-bonded to a C that is double-bonded to O or S
    if (g.atoms[i].z !== 7) return false;
    return g.adj[i].some(e => { const b = g.bonds[e.bond]; if (b.order !== 1 || b.arom || g.atoms[e.nb].z !== 6) return false;
      return g.adj[e.nb].some(f => f.nb !== i && g.bonds[f.bond].order === 2 && !g.bonds[f.bond].arom && [8, 16].includes(g.atoms[f.nb].z)); });
  }
  // elements at the other end of the non-aromatic double/triple bonds of atom i (one entry per bond)
  const multipleBondPartners = (g, i) => g.adj[i].filter(e => g.bonds[e.bond].order >= 2 && !g.bonds[e.bond].arom).map(e => SYM(g.atoms[e.nb].z));
  function atomTextbook(g, i, ctx) {
    const a = g.atoms[i];
    if (a.z === 0) return { lewis: null, vsepr: null, hybridization: { label: null, note: 'unspecified atom (*): not assigned' } };
    const out = { lewis: LEWIS.has(a.z) ? lewis(g, i) : null, vsepr: null };
    const nSigma = totalDegree(g, i), L = out.lewis;
    if (!L || L.lp === null) {
      out.hybridization = { label: null, note: LEWIS.has(a.z) ? 'unusual electron count: not assigned' : 'metal: not assigned' };
      return out;
    }
    if (a.z === 1) {
      out.hybridization = { label: null, orbital: '1s', note: 'hydrogen uses its 1s orbital: not hybridized' };
      out.octet = nSigma === 0 ? (L.shell === 2 ? 'hydride ion: 2 electrons' : null) : (L.shell === 2 ? 'duet' : 'not 2');
      return out;
    }
    if (nSigma === 0) { out.hybridization = { label: null, note: 'single atom / ion: not assigned' }; return out; }
    if (a.rad) {
      const planarC = a.z === 6 && nSigma === 3 && L.lp === 0;
      out.hybridization = planarC ? { label: 'sp2', conventional: true, note: 'carbon radical: (nearly) planar, conventionally sp²' }
                                  : { label: null, note: 'radical: simple VSEPR does not apply' };
      return out;
    }
    const SN = nSigma + L.lp;
    out.stericNumber = SN;
    let promoted = null;
    if (SN === 4 && L.lp >= 1 && nSigma <= 3) {
      if (a.arom) promoted = 'aromatic';
      else if (g.adj[i].some(e => ctx.conj[e.bond])) promoted = isAmideN(g, i) ? 'amide' : 'conjugated';
      else if (a.z <= 10 && g.adj[i].some(e => isCarbocation(g, e.nb))) promoted = 'cation';
    }
    const effSN = promoted ? 3 : SN, effLP = promoted ? L.lp - 1 : L.lp;
    const hybLabel = { 2: 'sp', 3: 'sp2', 4: 'sp3', 5: 'sp3d', 6: 'sp3d2' }[effSN] || null;
    const hypervalent = L.shell > 8 && a.z > 10;
    const terminal = nSigma === 1 && a.z !== 1;
    out.hybridization = { label: hybLabel, promoted,
      conventional: terminal || promoted === 'conjugated' || promoted === 'cation' || hypervalent || effSN >= 5,
      labelOnly: effSN >= 5, terminal, hypervalent, partners: hypervalent ? multipleBondPartners(g, i) : [] };
    out.effSN = effSN; out.effLP = effLP;
    const key = effSN + ',' + effLP;
    if (terminal) out.vsepr = { electronDomain: EDG_ONLY[effSN] || null, molecular: 'terminal atom (one bond: no bond angle)', ideal: null };
    else if (VSEPR[key]) out.vsepr = { electronDomain: VSEPR[key][0], molecular: VSEPR[key][1], ideal: VSEPR[key][2] };
    else out.vsepr = { electronDomain: EDG_ONLY[effSN] || null, molecular: null, ideal: null };
    if (!terminal && out.vsepr.ideal) {
      if (promoted === 'conjugated' || promoted === 'cation') {
        out.vsepr.ideal = null;
        out.vsepr.idealNote = 'between the sp³ (109.5°) and sp² (120°) values: the lone pair is only partly delocalised; see the measured angles';
      } else if (period(a.z) >= 3 && effSN === 4 && effLP >= 1) {
        out.vsepr.ideal = 'well below 109.5° (about 90–107° for P, S and heavier atoms)';
      }
      const mb = multipleBondPartners(g, i);
      if (hypervalent && effSN === 4 && effLP === 0 && mb.length) out.vsepr.idealNote = `VSEPR counts ${mb.length === 1 ? 'the drawn' : 'each drawn'} ` +
        `${[...new Set(mb)].map(x => SYM(a.z) + '=' + x).join(' / ')} bond as one domain; real angles can differ by several degrees (see the measured angles)`;
      const smallRings = g.rings.filter(r => r.includes(i) && r.length <= 5);
      if (smallRings.length) {
        const k = Math.min(...smallRings.map(r => r.length)), poly = { 3: 60, 4: 90, 5: 108 }[k];
        // the regular-polygon angle holds roughly only for rings of C, N and O: bonds to S, P and heavier atoms are long,
        // which closes the angle at that atom (thiophene C–S–C 92.2°, microwave structure) and opens the others
        // (thiophene's carbons about 112°), so no polygon value is given for any atom of such a ring
        const heavy = [...new Set(smallRings.filter(r => r.length === k).flat().filter(x => period(g.atoms[x].z) >= 3).map(x => SYM(g.atoms[x].z)))];
        out.vsepr.ringNote = period(a.z) >= 3
          ? `inside the ${k}-membered ring the angle at ${SYM(a.z)} is forced well below the regular-polygon ${poly}°, because bonds to ` +
            `${SYM(a.z)} are long` + (k === 5 ? ' (C–S–C is about 92° in thiophene)' : '') + '; see the measured angle'
          : heavy.length ? `inside a ${k}-membered ring that contains ${heavy.join(' and ')}: its long bonds to ${heavy.length === 1 ? 'that atom' : 'those atoms'} ` +
            'make the ring angles unequal, so the regular-polygon value does not apply; see the measured angle'
          : `inside the ${k}-membered ring the angle is forced to about ${poly}°`;
      }
    }
    out.octet = L.shell === 8 ? 'octet' : L.shell < 8 ? 'incomplete octet' : 'expanded octet (as drawn)';
    return out;
  }

  /* ================================================================ 6. Kekulé structures, oxidation state, aromatic π count */
  function kekuleStats(g, cap) {
    cap = cap || 200000;
    const frac = g.bonds.map(b => (b.arom ? null : b.order)), dblFrac = g.bonds.map(() => null), totalAt = g.bonds.map(() => null);
    const V = new Set();
    g.bonds.forEach(b => { if (b.arom && b.order === 2) { V.add(b.a); V.add(b.b); } });
    const seen = new Set(), comps = [];
    for (const s of V) {
      if (seen.has(s)) continue;
      const comp = [], st = [s]; seen.add(s);
      while (st.length) { const u = st.pop(); comp.push(u);
        for (const e of g.adj[u]) if (g.bonds[e.bond].arom && V.has(e.nb) && !seen.has(e.nb)) { seen.add(e.nb); st.push(e.nb); } }
      comps.push(comp.sort((x, y) => x - y));
    }
    const compCount = new Map(); let capped = false;
    for (const comp of comps) {
      const inC = new Set(comp), matched = new Map(), stackB = [], dbl = new Map(); let total = 0;
      const rec = () => {
        if (capped) return;
        let u = -1; for (const x of comp) if (!matched.has(x)) { u = x; break; }
        if (u < 0) { total++; for (const k of stackB) dbl.set(k, (dbl.get(k) || 0) + 1); if (total > cap) capped = true; return; }
        for (const e of g.adj[u]) {
          const b = g.bonds[e.bond];
          if (!b.arom || !inC.has(e.nb) || matched.has(e.nb)) continue;
          matched.set(u, e.nb); matched.set(e.nb, u); stackB.push(e.bond);
          rec();
          stackB.pop(); matched.delete(u); matched.delete(e.nb);
        }
      };
      rec();
      for (const x of comp) compCount.set(x, total);
      g.bonds.forEach((b, k) => { if (b.arom && inC.has(b.a) && inC.has(b.b)) {
        frac[k] = 1 + (dbl.get(k) || 0) / total; dblFrac[k] = (dbl.get(k) || 0); totalAt[k] = total; } });
    }
    g.bonds.forEach((b, k) => { if (b.arom && frac[k] === null) frac[k] = 1; });   // pyrrole-type aromatic bonds: always single
    return { capped, avgOrder: frac, doubleIn: dblFrac, structuresOf: totalAt, structuresAt: i => compCount.get(i) || 1 };
  }

  // Oxidation state of an atom (element z, formal charge chg) from its neighbours [{z, order}] (order may be a Kekulé
  // average); scale 'pauling' (shown) or 'allen' (only to state what the IUPAC 2016 scale would give).
  function oxidationFrom(z, chg, nbrs, scale) {
    const en = scale === 'allen' ? ALLEN : EN;
    if (en(z) == null) return null;
    let os = chg; const borderline = new Map();
    for (const { z: zn, order } of nbrs) {
      if (en(zn) == null) return null;
      if (zn === z) continue;
      const d = en(zn) - en(z);
      if (scale !== 'allen') {
        const al = ALLEN(zn) != null && ALLEN(z) != null ? ALLEN(zn) - ALLEN(z) : null;
        const scalesDisagree = al !== null && Math.sign(al) !== Math.sign(d);
        if (Math.abs(d) < 0.05 || scalesDisagree)
          borderline.set(zn, { with: SYM(zn), z: zn, dPauling: Math.round(Math.abs(d) * 100) / 100, allenAgrees: al === null ? null : !scalesDisagree });
      }
      os += d > 0 ? order : d < 0 ? -order : 0;         // equal EN, different element: split (0)
    }
    return { value: Math.round(os * 1e6) / 1e6, borderline: [...borderline.values()] };
  }
  function neighbourOrders(g, i, kek) {
    const out = g.adj[i].map(e => ({ z: g.atoms[e.nb].z, order: kek.avgOrder[e.bond] }));
    for (let t = 0; t < g.atoms[i].hs; t++) out.push({ z: 1, order: 1 });
    return out;
  }

  function aromaticSystems(g) {
    const n = g.atoms.length, seen = new Array(n).fill(false), out = [];
    for (let s = 0; s < n; s++) {
      if (!g.atoms[s].arom || seen[s]) continue;
      const comp = [], st = [s]; seen[s] = true;
      while (st.length) { const u = st.pop(); comp.push(u);
        for (const e of g.adj[u]) if (g.bonds[e.bond].arom && !seen[e.nb]) { seen[e.nb] = true; st.push(e.nb); } }
      let pi = 0;
      for (const u of comp) {
        const inRingDouble = g.adj[u].some(e => g.bonds[e.bond].arom && g.bonds[e.bond].order === 2);
        const exoDouble = g.adj[u].some(e => !g.bonds[e.bond].arom && g.bonds[e.bond].order >= 2);
        const L = lewis(g, u);
        pi += inRingDouble ? 1 : exoDouble ? 0 : (L && L.lp >= 1 ? 2 : 0);
      }
      out.push({ atoms: comp.sort((x, y) => x - y), piElectrons: pi, huckelN: (pi - 2) % 4 === 0 && pi >= 2 ? (pi - 2) / 4 : null });
    }
    return out;
  }

  /* ================================================================ 7. bond values */
  const METAL_IONIC = new Set([11, 19, 37, 55, 87, 20, 38, 56, 88]);   // Na K Rb Cs Fr Ca Sr Ba Ra
  function polarity(za, zb) {
    const ea = EN(za), eb = EN(zb);
    if (ea == null || eb == null) return null;
    const d = Math.round(Math.abs(ea - eb) * 100) / 100;
    let cls = d < 0.4 ? 'nonpolar' : d < 2.0 ? 'polar covalent' : 'largely ionic';
    if ((METAL_IONIC.has(za) || METAL_IONIC.has(zb)) && !(METAL_IONIC.has(za) && METAL_IONIC.has(zb))) cls = 'largely ionic';
    return { dEN: d, cls, negativeEnd: ea > eb ? 'a' : eb > ea ? 'b' : null };
  }
  function ringBondSet(g) {
    if (!g._ringBonds) {
      g._ringBonds = new Set();
      for (const r of g.rings) for (let t = 0; t < r.length; t++) {
        const u = r[t], v = r[(t + 1) % r.length];
        const e = g.adj[u].find(x => x.nb === v); if (e) g._ringBonds.add(e.bond);
      }
    }
    return g._ringBonds;
  }
  const inRingBond = (g, k) => ringBondSet(g).has(k);
  /* RDKit's strict rotatable-bond descriptor per bond. The reasons are the descriptor's rules, not physics: ethane's
     C–C is "not counted" although it rotates freely (see rotationInfo for the textbook description). */
  const ROT_REASON = {
    arom: 'aromatic ring bond', multiple: 'double or triple bond', ring: 'ring bond',
    terminal: 'bond to a terminal heavy atom: turning it moves only hydrogen atoms',
    triple: 'next to a triple bond: the atoms along it lie on one straight line',
    symmetric: 'symmetric end group (CH₃, CF₃, CCl₃, CBr₃, tert-butyl): turning it gives an equivalent position',
    amide: 'amide / ester / thioester / amidine-type bond (partial double bond)',
    formyl: "formamide / formate-type C(=O)–N/O bond (partial double bond); RDKit's strict count includes it, the app does not",
    hydrogen: 'bond to a hydrogen atom',
  };
  function strictRotatable(g, k) {    // == RDKit NumRotatableBonds (Strict) per bond, plus the formyl deviation
    const b = g.bonds[k];
    if (!(b.order === 1 || b.arom) || (b.arom && inRingBond(g, k))) return { rotatable: false, rdkitStrict: false, code: b.arom ? 'arom' : 'multiple', reason: ROT_REASON[b.arom ? 'arom' : 'multiple'] };
    if (inRingBond(g, k)) return { rotatable: false, rdkitStrict: false, code: 'ring', reason: ROT_REASON.ring };
    const D = i => degree(g, i);
    const aliph = i => !g.atoms[i].arom;
    const isC = i => g.atoms[i].z === 6 && aliph(i);
    const inTriple = i => g.adj[i].some(e => g.bonds[e.bond].order === 3 && !g.bonds[e.bond].arom);
    const cx3 = i => isC(i) && [9, 17, 35].some(x => g.adj[i].filter(e => g.atoms[e.nb].z === x && (g.bonds[e.bond].order === 1 || g.bonds[e.bond].arom)).length >= 3);
    const isCH3 = i => isC(i) && nH(g, i) === 3;
    const tBu = i => isC(i) && g.adj[i].filter(e => isCH3(e.nb) && (g.bonds[e.bond].order === 1 || g.bonds[e.bond].arom)).length >= 3;
    const common = i => D(i) !== 1 && !inTriple(i) && !cx3(i) && !tBu(i) && !isCH3(i);
    const single = e => g.bonds[e.bond].order === 1 && !g.bonds[e.bond].arom;
    const dbl = e => g.bonds[e.bond].order === 2 && !g.bonds[e.bond].arom;
    const nonRingSingle = e => single(e) && !inRingBond(g, e.bond);
    const isXhet = j => g.atoms[j].z === 7 || (aliph(j) && g.atoms[j].z === 8) || (aliph(j) && g.atoms[j].z === 16 && D(j) !== 1);
    const acylC = (c, plus) => isC(c) && D(c) === 3 && g.adj[c].some(e => dbl(e) && aliph(e.nb) &&
      (plus ? (g.atoms[e.nb].z === 7 && g.atoms[e.nb].chg === 1) : [7, 8, 16].includes(g.atoms[e.nb].z)));
    const amideLike = i =>
      (acylC(i, false) && g.adj[i].some(e => nonRingSingle(e) && isXhet(e.nb))) ||
      (isXhet(i) && g.adj[i].some(e => nonRingSingle(e) && acylC(e.nb, false))) ||
      (acylC(i, true) && g.adj[i].some(e => nonRingSingle(e) && g.atoms[e.nb].z === 7 && D(e.nb) !== 1)) ||
      (g.atoms[i].z === 7 && D(i) !== 1 && g.adj[i].some(e => nonRingSingle(e) && acylC(e.nb, true)));
    const first = i => common(i) && !amideLike(i);
    const ok = (first(b.a) && common(b.b)) || (first(b.b) && common(b.a));
    if (ok) {
      const formylC = i => isC(i) && D(i) === 2 && nH(g, i) === 1 &&
        g.adj[i].some(e => dbl(e) && aliph(e.nb) && [7, 8, 16].includes(g.atoms[e.nb].z));
      if ((formylC(b.a) && isXhet(b.b)) || (formylC(b.b) && isXhet(b.a)))
        return { rotatable: false, rdkitStrict: true, code: 'formyl', reason: ROT_REASON.formyl };
      return { rotatable: true, rdkitStrict: true, code: null, reason: null };
    }
    const code = [b.a, b.b].some(i => isCH3(i) || cx3(i) || tBu(i)) ? 'symmetric' : D(b.a) === 1 || D(b.b) === 1 ? 'terminal'
      : inTriple(b.a) || inTriple(b.b) ? 'triple' : 'amide';
    return { rotatable: false, rdkitStrict: false, code, reason: ROT_REASON[code] };
  }
  /* Textbook description of rotation about a bond (the owner's "free rotation" question), separate from the
     descriptor above: single bonds rotate freely (ethane), except ring bonds, the partial double bonds of amides,
     esters, acids, thioesters and amidines, and the axis of a hindered biaryl (atropisomers). */
  function partialDoubleBond(g, k) {
    const b = g.bonds[k];
    if (b.order !== 1 || b.arom || inRingBond(g, k)) return false;
    const acyl = c => g.atoms[c].z === 6 && !g.atoms[c].arom && g.adj[c].some(e => { const x = g.bonds[e.bond];
      return x.order === 2 && !x.arom && [7, 8, 16].includes(g.atoms[e.nb].z); });
    const donor = x => { const a = g.atoms[x];
      return !a.arom && a.chg === 0 && (a.z === 7 || a.z === 8 || (a.z === 16 && degree(g, x) <= 2)) && sumBO(g, x) === DV(a.z); };
    return (acyl(b.a) && donor(b.b)) || (acyl(b.b) && donor(b.a));
  }
  const TERMINAL_ROTATION = 'not applicable: bond to a terminal atom (turning it changes nothing)';
  function rotationInfo(g, k, ctx, hvDrawn) {         // hvDrawn: a "double" bond of a hypervalent atom (S=O, P=O)
    const b = g.bonds[k], terminal = i => totalDegree(g, i) === 1;
    if (b.arom) return { free: false, text: 'restricted: aromatic ring bond' };
    if (terminal(b.a) || terminal(b.b)) return { free: null, text: TERMINAL_ROTATION };
    if (b.order === 3) return { free: null, text: 'no conformations: the atoms along a triple bond lie on one straight line' };
    if (b.order === 2) return { free: false, text: hvDrawn ? 'restricted (drawn as a double bond)' : 'restricted: the π bond prevents rotation' };
    if (inRingBond(g, k)) return { free: false, text: 'restricted: ring bond (only limited ring flexing)' };
    if (partialDoubleBond(g, k)) return { free: false, text: 'restricted: partial double bond (amide / ester resonance)' };
    if (ctx.biaryl && ctx.biaryl.has(k)) return { free: false, text: 'restricted: hindered biaryl axis (atropisomers possible)' };
    return { free: true, text: 'free' };
  }
  /* Non-ring aryl–aryl bonds with at least three ortho substituents: rotation about them can be slow enough for
     atropisomers (BINOL). Used for the rotation row and by the chirality rule (axialPlanarFlags). An ortho atom with a
     third heavy neighbour counts once; a ring-fusion ortho atom of two six-membered rings, next to an axis atom in a
     six-membered ring, counts twice: the fused ring's peri C–H points at the axis like a second ortho group, which is
     why 1,1′-binaphthyl (two such fusions, no other ortho group) gives atropisomers. */
  function orthoCount(g, x, other) {
    let n = 0;
    const six = r => r.length === 6;
    for (const e of g.adj[x]) {
      const o = e.nb;
      if (o === other || !g.atoms[o].arom || degree(g, o) < 3) continue;
      n++;
      const ringsO = g.rings.filter(r => r.includes(o));
      if (ringsO.length >= 2 && ringsO.every(six) && g.rings.some(r => six(r) && r.includes(o) && r.includes(x))) n++;
    }
    return n;
  }
  function hinderedBiarylBonds(g) {
    const out = new Set();
    g.bonds.forEach((b, k) => {
      if (b.order !== 1 || b.arom || inRingBond(g, k) || !g.atoms[b.a].arom || !g.atoms[b.b].arom) return;
      if (orthoCount(g, b.a, b.b) + orthoCount(g, b.b, b.a) >= 3) out.add(k);
    });
    return out;
  }
  /* Bonds that resonance makes equivalent, so that no single average bond enthalpy fits the Kekulé order drawn:
     X(=E)(–E⁻) (carboxylate, nitro, nitrate, sulfonate, phosphate anions), azide and diazo N⁺ (C=N⁺=N⁻, N=N⁺=N⁻,
     C⁻–N⁺≡N), amidinium / guanidinium C(=N⁺)–N. */
  function resonanceBonds(g) {
    const out = new Set(), terminal = x => degree(g, x) === 1;
    g.atoms.forEach((a, i) => {
      const nbs = g.adj[i].map(e => ({ nb: e.nb, bond: e.bond, b: g.bonds[e.bond], x: g.atoms[e.nb] })).filter(o => !o.b.arom);
      for (const E of [8, 16]) {
        const dbl = nbs.filter(o => o.x.z === E && terminal(o.nb) && o.b.order === 2 && o.x.chg === 0);
        const neg = nbs.filter(o => o.x.z === E && terminal(o.nb) && o.b.order === 1 && o.x.chg === -1);
        if (dbl.length && neg.length) dbl.concat(neg).forEach(o => out.add(o.bond));
      }
      if (a.z === 7 && a.chg === 1 && nbs.length === 2 && degree(g, i) === 2 && a.hs === 0 && nbs[0].b.order + nbs[1].b.order === 4) {
        const t = nbs.find(o => o.x.z === 7 && terminal(o.nb)), other = nbs.find(o => o !== t);
        if (t && other && (t.x.chg === -1 || other.x.chg === -1 || other.b.order === 2)) nbs.forEach(o => out.add(o.bond));
      }
      if (a.z === 6 && a.chg === 0 && !a.arom) {
        const plus = nbs.filter(o => o.x.z === 7 && o.x.chg === 1 && o.b.order === 2);
        const lp = nbs.filter(o => o.x.z === 7 && o.x.chg === 0 && o.b.order === 1 && !o.x.arom && sumBO(g, o.nb) === 3);
        if (plus.length && lp.length) plus.concat(lp).forEach(o => out.add(o.bond));
      }
    });
    return out;
  }
  function bondTextbook(g, k, ctx) {
    const b = g.bonds[k], A = g.atoms[b.a], B = g.atoms[b.b], out = {};
    const hyb = i => {
      const t = ctx.atoms[i].hybridization;
      if (g.atoms[i].z === 1) return '1s';
      return t && t.label && !t.labelOnly ? hybText(t.label) : null;
    };
    const p = z => period(z) + 'p';
    const aSym = SYM(A.z), bSym = SYM(B.z);
    const hv = [b.a, b.b].some(i => ctx.atoms[i].hybridization && ctx.atoms[i].hybridization.hypervalent);
    if (b.arom) {
      const f = ctx.kek.avgOrder[k];
      out.type = 'aromatic'; out.kekule = b.order === 2 ? 'double' : 'single';
      out.bondOrder = !ctx.kek.capped && Math.abs(f - 1.5) < 1e-9 ? 1.5 : null;
      out.sigma = 1; out.pi = null;
    } else {
      out.type = { 1: 'single', 2: 'double', 3: 'triple' }[b.order] || 'other';
      out.kekule = null; out.bondOrder = b.order; out.sigma = 1; out.pi = b.order - 1;
      if (hv && b.order >= 2) out.piNote = 'drawn as a double bond; better described as a polar single bond ' +
        (ctx.atoms[b.a].hybridization && ctx.atoms[b.a].hybridization.hypervalent ? aSym + '⁺–' + bSym + '⁻' : bSym + '⁺–' + aSym + '⁻') +
        ' (no d-orbital π bond)';
    }
    const ha = hyb(b.a), hb = hyb(b.b);
    out.sigmaOverlap = ha && hb ? `σ: ${aSym}(${ha})–${bSym}(${hb})` : null;
    out.piOverlap = !b.arom && b.order >= 2 && !out.piNote ? `π: ${aSym}(${p(A.z)})–${bSym}(${p(B.z)})` + (b.order === 3 ? ' ×2 (perpendicular)' : '') : null;
    out.polarity = polarity(A.z, B.z);
    out.conjugated = ctx.conjShown[k];
    out.rdkitConjugated = ctx.conj[k];
    out.rotatable = strictRotatable(g, k);
    out.rotation = rotationInfo(g, k, ctx, !!out.piNote);
    out.amideCN = b.order === 1 && !b.arom && ((A.z === 6 && B.z === 7 && isAmideN(g, b.b) && g.adj[b.a].some(e => g.bonds[e.bond].order === 2 && [8, 16].includes(g.atoms[e.nb].z))) ||
                                               (B.z === 6 && A.z === 7 && isAmideN(g, b.a) && g.adj[b.b].some(e => g.bonds[e.bond].order === 2 && [8, 16].includes(g.atoms[e.nb].z))));
    return out;
  }

  /* ================================================================ 8. geometry primitives */
  function angle(p0, p1, p2) {                          // degrees, angle at p1
    const ux = p0[0] - p1[0], uy = p0[1] - p1[1], uz = p0[2] - p1[2], vx = p2[0] - p1[0], vy = p2[1] - p1[1], vz = p2[2] - p1[2];
    const c = (ux * vx + uy * vy + uz * vz) / Math.sqrt((ux * ux + uy * uy + uz * uz) * (vx * vx + vy * vy + vz * vz));
    return Math.acos(Math.max(-1, Math.min(1, c))) * 180 / Math.PI;
  }
  function dihedral(p0, p1, p2, p3) {                   // IUPAC sign: + = clockwise looking from p1 to p2 (== RDKit)
    const sub = (p, q) => [p[0] - q[0], p[1] - q[1], p[2] - q[2]], dot = (p, q) => p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
    const b0 = sub(p0, p1), b1 = sub(p2, p1), b2 = sub(p3, p2);
    const n1 = Math.sqrt(dot(b1, b1)), u = [b1[0] / n1, b1[1] / n1, b1[2] / n1];
    const v = sub(b0, u.map(x => x * dot(b0, u))), w = sub(b2, u.map(x => x * dot(b2, u)));
    const cx = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]];
    return Math.atan2(dot(cx, w), dot(v, w)) * 180 / Math.PI;
  }
  // IUPAC Gold Book "torsion angle" (Klyne & Prelog 1960). `common` is the everyday word; gauche and eclipsed are
  // words for sp³–sp³ bonds (butane), so they are given only when sp3sp3 is true (syn / anti are used for any bond).
  function klynePrelog(t, sp3sp3) {
    const a = Math.abs(t);
    if (a < 30) return { code: 'sp', name: 'synperiplanar', common: sp3sp3 ? 'syn; eclipsed if ≈0°' : 'syn' };
    if (a < 90) return { code: 'sc', name: 'synclinal', common: sp3sp3 ? 'gauche' : null };
    if (a < 150) return { code: 'ac', name: 'anticlinal', common: sp3sp3 ? 'eclipsed if ≈120°' : null };
    return { code: 'ap', name: 'antiperiplanar', common: 'anti' };
  }
  function refNeighbours(g, from, to) {                 // reference atoms on `from`, looking away from `to`, best first
    const cand = g.adj[from].filter(e => e.nb !== to && g.atoms[e.nb].z > 1).map(e => e.nb);
    const key = x => [g.atoms[x].z].concat(
      g.adj[x].filter(e => e.nb !== from).map(e => g.atoms[e.nb].z).concat(new Array(g.atoms[x].hs).fill(1)).sort((p, q) => q - p));
    const cmp = (x, y) => { const kx = key(x), ky = key(y);
      for (let t = 0; t < Math.max(kx.length, ky.length); t++) { const d = (ky[t] || 0) - (kx[t] || 0); if (d) return d; }
      return x - y; };
    return cand.sort(cmp);
  }
  function torsionAtoms(g, k) {
    const b = g.bonds[k];
    if (b.order !== 1 || b.arom) return null;
    const c1 = refNeighbours(g, b.a, b.b), c2 = refNeighbours(g, b.b, b.a);
    if (!c1.length || !c2.length) return null;
    // the two end atoms must differ (in a 3-membered ring both ends see the same third atom first: C2–C1–C3–C2 is not
    // a dihedral); take the next candidate on either side, else there is no dihedral to show
    for (const r1 of c1) for (const r2 of c2) if (r1 !== r2) return [r1, b.a, b.b, r2];
    return null;
  }

  /* ================================================================ 9. stereo labels (get_stereo_tags = RDKit's new CIP labeller) */
  /* A P, As, S, Se or Te atom with two terminal O (or two terminal S) of the same isotope that exchange by moving a
     proton or a charge (=O / –OH / –O⁻) is not a stereocentre: phosphate and phosphodiester P (ATP's Pα, Pβ), phosphonic
     and sulfinic acids. RDKit's FindPotentialStereo does not know this; standard InChI drops such centres from its /t
     layer. One O and one S (phosphorothioates, ATPαS) stay real stereocentres. */
  const OXO_CENTRE = new Set([15, 33, 16, 34, 52]);
  function tautomericOxoCentre(g, i) {
    if (!OXO_CENTRE.has(g.atoms[i].z)) return null;
    const byKey = new Map();
    for (const e of g.adj[i]) {
      const x = g.atoms[e.nb], b = g.bonds[e.bond];
      if ((x.z !== 8 && x.z !== 16) || degree(g, e.nb) !== 1 || b.arom) continue;          // terminal O / S only
      const ok = (b.order === 2 && x.chg === 0) || (b.order === 1 && (x.chg === -1 || (x.chg === 0 && x.hs === 1)));
      if (!ok) continue;
      const key = x.z + ':' + (x.iso || 0);
      byKey.set(key, (byKey.get(key) || 0) + 1);
    }
    for (const [key, c] of byKey) if (c >= 2) {
      const E = key.startsWith('8:') ? 'O' : 'S';
      return `not a stereocentre: its two terminal ${E} atoms exchange by moving a proton or a charge (=${E} / –${E}H / –${E}⁻)`;
    }
    return null;
  }
  /* A 2D copy of M in which each atom of `marks` can be given a wedge (1) or a hash (6): it gets one bond of its own
     that carries no stereo mark yet, turned so that the atom is the bond's first atom (specified stereo keeps its own
     marks). → {read(flags, fn), close()} or null (flags: Map atom → 1 | 6; fn gets the RDKit molecule read back) */
  function wedgeProbe(R, M, g, marks) {
    let c = null;
    try {
      c = M.copy(); c.set_new_coords();
      const lines = c.get_molblock().split('\n');
      const na = +lines[3].slice(0, 3), nb = +lines[3].slice(3, 6);
      if (na !== g.atoms.length || nb !== g.bonds.length) return null;
      const bl = lines.slice(4 + na, 4 + na + nb), lineOf = new Map(), used = new Set();
      for (const i of marks) {
        for (let k = 0; k < bl.length; k++) {
          if (used.has(k) || +bl[k].slice(9, 12)) continue;
          const x = +bl[k].slice(0, 3) - 1, y = +bl[k].slice(3, 6) - 1;
          if (x !== i && y !== i) continue;
          bl[k] = String(i + 1).padStart(3) + String((x === i ? y : x) + 1).padStart(3) + bl[k].slice(6);
          lineOf.set(i, k); used.add(k); break;
        }
      }
      const atomOfLine = new Map([...lineOf].map(([i, k]) => [k, i]));
      const withFlags = flags => lines.slice(0, 4 + na).concat(bl.map((line, k) => {
        const i = atomOfLine.get(k);
        return i !== undefined && flags.has(i) ? line.slice(0, 9) + String(flags.get(i)).padStart(3) + line.slice(12) : line;
      }), lines.slice(4 + na + nb)).join('\n');
      return {
        marked: i => lineOf.has(i),
        read: (flags, fn) => { let m2 = null; try { m2 = R.get_mol(withFlags(flags)); return m2 ? fn(m2) : null; } catch (e) { return null; } finally { if (m2) m2.delete(); } },
      };
    } catch (e) { return null; }
    finally { if (c) c.delete(); }
  }
  /* Are a and b the two bridgeheads of a small bridged ring system? Three separate a–b paths (bridges) of at least one
     atom each, and both smallest rings through a and b have ≤ 7 members (norbornane, bicyclo[2.2.2]octane, [3.2.1],
     [3.3.1], [4.4.1]): then only one relative configuration of the two fits (no trans bridge, no in/out isomers).
     The same rule as js/mol3d_engine.js smallBridgeheads (greedy shortest paths: a miss only keeps a centre open). */
  function smallBridgeheads(g, a, b) {
    if (g.adj[a].some(e => e.nb === b)) return false;
    const used = new Set(), lens = [];
    for (let t = 0; t < 3; t++) {
      const prev = new Map([[a, -1]]), queue = [a];
      let found = false;
      for (let h = 0; h < queue.length && !found; h++) {
        for (const { nb: y } of g.adj[queue[h]]) {
          if (prev.has(y) || used.has(y)) continue;
          prev.set(y, queue[h]);
          if (y === b) { found = true; break; }
          queue.push(y);
        }
      }
      if (!found) return false;
      let n = 0;
      for (let x = prev.get(b); x !== a; x = prev.get(x)) { used.add(x); n++; }
      lens.push(n);
    }
    lens.sort((x, y) => x - y);
    return lens[0] >= 1 && lens[0] + lens[2] + 2 <= 7;
  }
  const BRIDGEHEAD_REASON = 'not a stereocentre choice: a bridgehead of a small bridged ring system, where only one ' +
    'configuration of the two bridgeheads fits (inverting both gives the same compound)';
  // sp³ atoms that could be stereocentres (four σ bonds, at most one H): the ones a probe holds fixed
  const stereoCapable = (g, i) => { const a = g.atoms[i];
    return a.z > 1 && !a.arom && totalDegree(g, i) === 4 && nH(g, i) <= 1 && g.adj[i].every(e => g.bonds[e.bond].order === 1 && !g.bonds[e.bond].arom); };
  /* RDKit proposes the two bridgeheads of norbornene as open stereocentres ('?'), but in a small bridged ring only one
     configuration of the pair fits. A pair is not a choice when, with every other possible stereocentre held fixed,
     inverting both bridgeheads gives the same compound (canonical SMILES and InChI) for one of the pair's two relative
     configurations (the other one is the impossible trans bridge). Camphor's bridgeheads stay open: inverting both
     gives its enantiomer. → Set of atom indices (M's numbering) */
  function fixedBridgeheads(R, M, g, tags) {
    const out = new Set();
    const open = tags.CIP_atoms.filter(([, lab]) => lab === '(?)').map(t => t[0]);
    const pairs = [];
    for (let x = 0; x < open.length; x++) for (let y = x + 1; y < open.length; y++)
      if (smallBridgeheads(g, open[x], open[y])) pairs.push([open[x], open[y]]);
    if (!pairs.length) return out;
    const specified = new Set(tags.CIP_atoms.filter(([, lab]) => lab !== '(?)').map(t => t[0]));
    const others = g.atoms.map((_, i) => i).filter(i => !specified.has(i) && stereoCapable(g, i));
    const P = wedgeProbe(R, M, g, others);
    if (!P) return out;
    const id = flags => P.read(flags, m2 => ({ smi: m2.get_smiles(), inchi: m2.get_inchi() }));
    const flip = f => (f === 1 ? 6 : 1);
    for (const [a, b] of pairs) {
      if (!P.marked(a) || !P.marked(b)) continue;
      const base = new Map(others.filter(i => i !== a && i !== b && P.marked(i)).map(i => [i, 1]));
      for (const [fa, fb] of [[1, 1], [1, 6]]) {
        const s1 = id(new Map([...base, [a, fa], [b, fb]])), s2 = id(new Map([...base, [a, flip(fa)], [b, flip(fb)]]));
        if (s1 && s2 && s1.smi === s2.smi && s1.inchi === s2.inchi) { out.add(a); out.add(b); break; }
      }
    }
    return out;
  }
  /* Ring cis/trans stereo that RDKit's CIP labeller does not list when it is unspecified: 1,4-disubstituted
     cyclohexanes (centres that become r/s once specified) and ring fusions (decalin, hydrindane). Ring atoms with
     two or three ring bonds, ≤ 1 H and no multiple bond are given a wedge in a 2D copy; the ones RDKit keeps as
     stereocentres are the open cis/trans atoms. Left out: the bridgeheads of small bridged rings (norbornane: one
     configuration only) and fusions where a trans ring junction is impossible (both rings together < 9 atoms:
     bicyclo[2.2.0]hexane, [3.1.0]hexane). With up to four such atoms every configuration of them is written out as
     SMILES (`isomers`), so that the chirality rule can test each ring isomer instead of guessing.
     → {atoms: [{index, label, ring: 'substituent' | 'fusion' | 'bridgehead'}], isomers: [smiles] | null} */
  const RING_ISOMER_MAX = 4;
  function openRingStereo(R, M, g, known) {
    const none = { atoms: [], isomers: null };
    const ringBonds = i => g.adj[i].filter(e => inRingBond(g, e.bond));
    const bondRingSizes = k => { const b = g.bonds[k];
      return g.rings.filter(r => { const ia = r.indexOf(b.a), ib = r.indexOf(b.b); return ia >= 0 && ib >= 0 && (Math.abs(ia - ib) === 1 || Math.abs(ia - ib) === r.length - 1); }).map(r => r.length); };
    let cand = [];
    for (let i = 0; i < g.atoms.length; i++) {
      if (known.has(i) || !stereoCapable(g, i)) continue;
      const nr = ringBonds(i).length;
      if (nr === 2 || nr === 3) cand.push(i);
    }
    const drop = new Set();
    for (const x of cand) for (const y of cand) if (x < y && smallBridgeheads(g, x, y)) { drop.add(x); drop.add(y); }
    const fusionPartner = i => ringBonds(i).find(e => cand.includes(e.nb) && ringBonds(e.nb).length === 3 && ringBonds(i).length === 3);
    for (const i of cand) {
      const f = fusionPartner(i);
      if (f) { const s = bondRingSizes(f.bond).sort((p, q) => p - q); if (s.length >= 2 && s[0] + s[1] < 9) drop.add(i); }
    }
    cand = cand.filter(i => !drop.has(i));
    if (cand.length < 2) return none;
    try {
      const P = wedgeProbe(R, M, g, cand);
      if (!P) return none;
      const t2 = P.read(new Map(cand.filter(P.marked).map(i => [i, 1])), m2 => JSON.parse(m2.get_stereo_tags()));
      if (!t2) return none;
      const atoms = t2.CIP_atoms.filter(([i, lab]) => cand.includes(i) && /\((R|S|r|s)\)/.test(lab))
        .map(([i, lab]) => ({ index: i, label: lab.replace(/[()]/g, ''),
          ring: ringBonds(i).length === 2 ? 'substituent' : fusionPartner(i) ? 'fusion' : 'bridgehead' }))
        .sort((x, y) => x.index - y.index);
      if (!atoms.length || atoms.length > RING_ISOMER_MAX) return { atoms, isomers: null };
      const isomers = new Set();
      for (let mask = 0; mask < (1 << atoms.length); mask++) {
        const smi = P.read(new Map(atoms.map((x, t) => [x.index, (mask >> t) & 1 ? 6 : 1])), m2 => m2.get_smiles());
        if (!smi) return { atoms, isomers: null };
        isomers.add(smi);
      }
      return { atoms, isomers: [...isomers] };
    } catch (e) { return none; }
  }
  const RING_TEXT = { substituent: 'cis/trans ring stereo: not specified', fusion: 'cis/trans ring fusion: not specified',
    bridgehead: 'bridgehead configuration (in/out): not specified' };
  function stereoInfo(g, tags, ringOpen, fixedBridge) {
    const ranks = g.cipRanks;
    const atoms = new Array(g.atoms.length).fill(null), bonds = new Array(g.bonds.length).fill(null), ignored = [];
    for (const [i, lab] of tags.CIP_atoms) {
      const L = lab.replace(/[()]/g, '');
      const why = tautomericOxoCentre(g, i);
      if (why) { ignored.push({ index: i, reason: why }); continue; }
      if (L === '?') {
        if (nH(g, i) >= 2) { ignored.push({ index: i, reason: `not a stereocentre: ${countWord(nH(g, i))} of its substituents are hydrogen` }); continue; }
        if (fixedBridge && fixedBridge.has(i)) { ignored.push({ index: i, reason: BRIDGEHEAD_REASON }); continue; }
        atoms[i] = { label: '?', kind: 'centre', text: 'possible stereocentre: configuration not specified' }; }
      else atoms[i] = { label: L, kind: 'centre', text: L === 'r' || L === 's' ? 'pseudoasymmetric centre (' + L + ')' : L };
    }
    for (const { index: i, ring } of ringOpen ? ringOpen.atoms : []) if (!atoms[i]) atoms[i] = { label: '?', kind: 'ring', ring, text: RING_TEXT[ring] };
    for (const [a, b, lab] of tags.CIP_bonds) {
      const e = g.adj[a].find(x => x.nb === b); if (!e) continue;
      const L = lab.replace(/[()]/g, '');
      bonds[e.bond] = L === 'E' || L === 'Z' ? { label: L, text: L } : { label: 'unspecified', text: 'E/Z possible: not specified in the input' };
    }
    const ringSize = k => { let best = 0; const b = g.bonds[k];
      for (const r of g.rings) { const ia = r.indexOf(b.a), ib = r.indexOf(b.b);
        if (ia >= 0 && ib >= 0 && (Math.abs(ia - ib) === 1 || Math.abs(ia - ib) === r.length - 1)) best = best ? Math.min(best, r.length) : r.length; }
      return best; };
    const endOK = (u, v) => {
      if (g.adj[u].filter(e => g.bonds[e.bond].order === 2 && !g.bonds[e.bond].arom).length > 1) return false;   // cumulated
      const subs = g.adj[u].filter(e => e.nb !== v).map(e => e.nb), h = g.atoms[u].hs + subs.filter(x => g.atoms[x].z === 1).length;
      const heavy = subs.filter(x => g.atoms[x].z > 1);
      if (g.atoms[u].z === 7 && g.atoms[u].chg === 0 && heavy.length + h === 1) return true;   // imine / oxime / azo N
      if (heavy.length + h !== 2 || h === 2) return false;
      if (h === 1) return true;
      return !!ranks && ranks[heavy[0]] !== ranks[heavy[1]];
    };
    g.bonds.forEach((b, k) => {
      if (bonds[k] || b.order !== 2 || b.arom) return;
      const rs = ringSize(k); if (rs && rs < 8) return;
      if (endOK(b.a, b.b) && endOK(b.b, b.a)) bonds[k] = { label: 'unspecified', text: 'E/Z possible: not specified in the input' };
    });
    return { atoms, bonds, ignored };
  }

  /* ================================================================ 10. molecule level: formula, functional groups, chirality */
  function hillOrder(counts) {                            // carbon first, then H, then alphabetical (no C: all alphabetical)
    const els = Object.keys(counts).filter(e => counts[e] > 0).sort();
    return counts.C ? ['C', ...(counts.H ? ['H'] : []), ...els.filter(e => e !== 'C' && e !== 'H')] : els;
  }
  function hill(counts, charge) {                         // "C2H3O2-", "H4N+", "Ca2+" (the reference format)
    let f = hillOrder(counts).map(e => e + (counts[e] > 1 ? counts[e] : '')).join('');
    if (charge) f += (Math.abs(charge) > 1 ? Math.abs(charge) : '') + (charge > 0 ? '+' : '-');
    return f;
  }
  const hillText = (counts, charge) =>                    // "C₂H₃O₂⁻", "Ca²⁺"
    hillOrder(counts).map(e => e + (counts[e] > 1 ? subDigits(String(counts[e])) : '')).join('') + chargeSup(charge);
  function isotopicFormula(atoms) {                      // CHD3O-style: isotopes listed after their element (D, T for H)
    const c = new Map();
    const add = (sym, iso, k) => { const key = sym + '|' + (iso || 0); c.set(key, (c.get(key) || 0) + k); };
    for (const a of atoms) { if (a.z > 0) add(SYM(a.z), a.iso, 1); if (a.hs) add('H', 0, a.hs); }
    const syms = [...new Set([...c.keys()].map(k => k.split('|')[0]))].sort();
    const order = syms.includes('C') ? ['C', ...(syms.includes('H') ? ['H'] : []), ...syms.filter(e => e !== 'C' && e !== 'H')] : syms;
    let f = '';
    for (const s of order) {
      const isos = [...c.keys()].filter(k => k.split('|')[0] === s).map(k => +k.split('|')[1]).sort((x, y) => x - y);
      for (const iso of isos) {
        const k = c.get(s + '|' + iso);
        const label = !iso ? s : s === 'H' && iso === 2 ? 'D' : s === 'H' && iso === 3 ? 'T' : supDigits(String(iso)) + s;
        f += label + (k > 1 ? subDigits(String(k)) : '');
      }
    }
    return f;
  }

  /* The functional-group table (order = precedence; SMARTS written for this app and checked against RDKit Python's
     perception of the same patterns and 521+ positive / negative test molecules in the fixtures). kind 'overlay'
     groups never claim atoms; 'claim' groups are dropped when one of their `claims` atoms (query indices; whole match
     if absent) was claimed by an earlier group (unless that group is listed in `sharesWith`). Highlight
     `atomsToHighlight`, else the whole match. Review changes: 'diazo', a neutral imine N, acetal O sharing. */
  const FG = [
    { name: 'arene (benzene ring)', smarts: 'c1ccccc1', kind: 'overlay' },
    { name: 'aromatic heterocycle (6-membered)', smarts: '[a;!#6]1aaaaa1', kind: 'overlay' },
    { name: 'aromatic heterocycle (5-membered)', smarts: '[a;!#6]1aaaa1', kind: 'overlay' },
    { name: 'α,β-unsaturated carbonyl', smarts: '[CX3]=[CX3]-[#6X3]=[OX1]', kind: 'overlay' },
    { name: 'acyl halide', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[F,Cl,Br,I]', kind: 'claim' },
    { name: 'anhydride', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[#8X2][#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])]=[OX1]', kind: 'claim' },
    { name: 'carbamate', smarts: '[#7X3][#6X3](=[OX1])[#8X2][#6]', kind: 'claim', claims: [1, 2, 3] },
    { name: 'urea', smarts: '[#7X3][#6X3](=[OX1])[#7X3]', kind: 'claim', claims: [1, 2] },
    { name: 'carbonate ester', smarts: '[#6][#8X2][#6X3](=[OX1])[#8X2][#6]', kind: 'claim', claims: [1, 2, 3, 4] },
    { name: 'imide', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[#7X3][#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])]=[OX1]', kind: 'claim', claims: [0, 1, 3, 4] },
    { name: 'carboxylic acid', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[OX2H1]', kind: 'claim' },
    { name: 'carboxylate', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[OX1-]', kind: 'claim' },
    { name: 'thioester', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[SX2][#6]', kind: 'claim', claims: [0, 1, 2] },
    { name: 'lactone', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])@[#8X2]', kind: 'claim' },
    { name: 'ester', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[#8X2][#6]', kind: 'claim', claims: [0, 1, 2] },
    { name: 'lactam', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])@[#7X3]', kind: 'claim', claims: [0, 1] },
    { name: 'amide (1°)', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[NX3;H2]', kind: 'claim', claims: [0, 1] },
    { name: 'amide (2°)', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[#7X3;H1]', kind: 'claim', claims: [0, 1] },
    { name: 'amide (3°)', smarts: '[#6X3;$([#6](=[OX1])~[#6]),$([#6;H1]=[OX1])](=[OX1])[#7X3;H0]', kind: 'claim', claims: [0, 1] },
    { name: 'sulfonic acid', smarts: '[#6][SX4](=[OX1])(=[OX1])[OX2H1]', kind: 'claim', claims: [1, 2, 3, 4] },
    { name: 'sulfonate (anion)', smarts: '[#6][SX4](=[OX1])(=[OX1])[OX1-]', kind: 'claim', claims: [1, 2, 3, 4] },
    { name: 'sulfonate ester', smarts: '[#6][SX4](=[OX1])(=[OX1])[OX2][#6]', kind: 'claim', claims: [1, 2, 3, 4] },
    { name: 'sulfonamide', smarts: '[#6][SX4](=[OX1])(=[OX1])[#7X3]', kind: 'claim', claims: [1, 2, 3, 4] },
    { name: 'sulfone', smarts: '[#6][SX4](=[OX1])(=[OX1])[#6]', kind: 'claim', claims: [1, 2, 3] },
    { name: 'sulfoxide', smarts: '[#6][SX3;$([SX3]=[OX1]),$([SX3+]-[OX1-])](~[OX1])[#6]', kind: 'claim', claims: [1, 2] },
    { name: 'phosphate ester', smarts: '[PX4](=[OX1])([#8X2][#6])([#8])[#8]', kind: 'claim', claims: [0, 1, 2, 4, 5] },
    { name: 'nitro', smarts: '[#6][NX3+](=[OX1])[OX1-]', kind: 'claim', claims: [1, 2, 3] },
    { name: 'nitrate ester', smarts: '[#6]-[#8X2]-[NX3+](=[OX1])[OX1-]', kind: 'claim', claims: [1, 2, 3, 4] },
    { name: 'azide', smarts: '[#6]-[#7X2]~[#7X2+]~[#7X1]', kind: 'claim', claims: [1, 2, 3] },
    { name: 'azo', smarts: '[#6]-[NX2]=[NX2]-[#6]', kind: 'claim', claims: [1, 2] },
    { name: 'isocyanate', smarts: '[#6]-[NX2]=[CX2]=[OX1]', kind: 'claim', claims: [1, 2, 3] },
    // C=N⁺=N⁻ or C⁻–N⁺≡N (diazomethane, diazo esters); not a diazonium Ar–N⁺≡N (its C is not charged)
    { name: 'diazo', smarts: '[#6X3;$([#6]=[#7X2+]=[#7X1-]),$([#6-]-[#7X2+]#[#7X1])]~[#7X2+]~[#7X1]', kind: 'claim' },
    { name: 'oxime', smarts: '[CX3]=[NX2]-[OX2H1]', kind: 'claim' },
    { name: 'hydrazone', smarts: '[CX3]=[NX2]-[NX3]', kind: 'claim' },
    { name: 'hydrazine', smarts: '[NX3;!$(N-[#6]=[!#6]);!$(N-[#6]#[!#6]);!$(N-[SX4])]-[NX3;!$(N-[#6]=[!#6]);!$(N-[#6]#[!#6]);!$(N-[SX4])]', kind: 'claim' },
    { name: 'guanidine', smarts: '[#7;!$([#7]~[!#6;!#1])]~[CX3;!R](~[#7;!$([#7]~[!#6;!#1])])~[#7;!$([#7]~[!#6;!#1])]', kind: 'claim', claims: [0, 1, 2, 3] },
    { name: 'amidine', smarts: '[CX3;$(C-[#6]),$([CH1]);!$(C(~[#7])(~[#7])~[#7])](=[NX2;!$(N-[!#6;!#1])])-[NX3;!$(N-[!#6;!#1])]', kind: 'claim' },
    { name: 'imine', smarts: '[CX3;!$(C-[!#6;!#1])]=[NX2;+0;!$(N-[!#6;!#1])]', kind: 'claim' },
    { name: 'nitrile', smarts: '[#6]-[CX2]#[NX1]', kind: 'claim', claims: [1, 2] },
    { name: 'quaternary ammonium', smarts: '[NX4+;H0](-[#6])(-[#6])(-[#6])-[#6]', kind: 'claim', claims: [0] },
    { name: 'ammonium (protonated amine)', smarts: '[NX4+;!H0;!$([#7]~[!#6;!#1])]', kind: 'claim' },
    { name: 'enamine', smarts: '[NX3;+0;!$(N~[!#6;!#1]);!$(N-[#6]=[!#6]);!$(N-[#6]#[!#6])]-[CX3]=[CX3]', kind: 'claim' },
    { name: 'aniline (aryl amine)', smarts: '[NX3;+0;!$(N~[!#6;!#1]);!$(N-[#6]=[!#6]);!$(N-[#6]#[!#6])]-c', kind: 'claim', claims: [0] },
    { name: 'amine (3°)', smarts: '[NX3;+0;!$(N~[!#6;!#1]);!$(N-[#6]=[!#6]);!$(N-[#6]#[!#6])](-[#6])(-[#6])-[#6]', kind: 'claim', claims: [0] },
    { name: 'amine (2°)', smarts: '[NX3;H1;+0;!$(N~[!#6;!#1]);!$(N-[#6]=[!#6]);!$(N-[#6]#[!#6])](-[#6])-[#6]', kind: 'claim', claims: [0] },
    { name: 'amine (1°)', smarts: '[NX3;H2;+0;!$(N~[!#6;!#1]);!$(N-[#6]=[!#6]);!$(N-[#6]#[!#6])]-[#6]', kind: 'claim', claims: [0] },
    { name: 'epoxide', smarts: '[#6]1-[#8X2]-[#6]1', kind: 'claim' },
    // sharesWith: an O already claimed by another acetal centre may be claimed again, so a glycosidic O shared by two
    // acetal centres (sucrose) does not hide the second one; both O stay claimed, so 'ether' cannot take them, and an O
    // claimed by an ester still blocks an acylal
    { name: 'acetal / ketal', smarts: '[#6]-[#8X2]-[CX4;!$(C(-[#8X2])(-[#8X2])-[!#6;!#1])]-[#8X2]-[#6]', kind: 'claim', claims: [1, 2, 3], sharesWith: ['acetal / ketal'] },
    { name: 'hemiacetal', smarts: '[#6]-[#8X2]-[CX4;!$(C(-[#8X2])(-[#8X2])-[!#6;!#1])]-[OX2H1]', kind: 'claim', claims: [1, 2, 3] },
    { name: 'enol', smarts: '[OX2H1]-[CX3]=[CX3]', kind: 'claim' },
    { name: 'phenol', smarts: '[OX2H1]-c', kind: 'claim', claims: [0] },
    { name: 'alcohol (3°)', smarts: '[OX2H1]-[CX4;!$(C(-[OX2H1])-[!#6;!#1])](-[#6])(-[#6])-[#6]', kind: 'claim', claims: [0] },
    { name: 'alcohol (2°)', smarts: '[OX2H1]-[CX4;!$(C(-[OX2H1])-[!#6;!#1]);!$(C(-[#6])(-[#6])-[#6])](-[#6])-[#6]', kind: 'claim', claims: [0] },
    { name: 'alcohol (1°)', smarts: '[OX2H1]-[CX4;!$(C(-[OX2H1])-[!#6;!#1]);!$(C(-[#6])-[#6])]-[#6]', kind: 'claim', claims: [0] },
    { name: 'alcohol (methanol)', smarts: '[OX2H1]-[CX4;!$(C(-[OX2H1])-[!#1])]', kind: 'claim', claims: [0] },
    { name: 'ether', smarts: '[#6]-[OX2;!$(O-[#6]=[!#6])]-[#6]', kind: 'claim', claims: [1] },
    { name: 'aldehyde', smarts: '[CX3;$([CH1]-[#6]),$([CH2])]=[OX1]', kind: 'claim' },
    { name: 'ketone', smarts: '[#6][#6X3](=[OX1])[#6]', kind: 'claim', claims: [1, 2] },
    { name: 'thiol', smarts: '[#6;!$([#6]=[!#6])]-[SX2H1]', kind: 'claim', claims: [1] },
    { name: 'disulfide', smarts: '[#6]-[SX2]-[SX2]-[#6]', kind: 'claim', claims: [1, 2] },
    { name: 'sulfide', smarts: '[#6]-[SX2;!$(S-[#6]=[!#6])]-[#6]', kind: 'claim', claims: [1] },
    { name: 'aryl halide', smarts: 'c-[F,Cl,Br,I]', kind: 'claim', claims: [1] },
    { name: 'vinyl halide', smarts: '[CX3](=[CX3])-[F,Cl,Br,I]', kind: 'claim', claims: [2], atomsToHighlight: [0, 2] },
    { name: 'alkyl halide', smarts: '[CX4]-[F,Cl,Br,I]', kind: 'claim', claims: [1] },
    { name: 'alkene', smarts: '[CX3]=[CX3]', kind: 'claim' },
    { name: 'alkyne', smarts: '[CX2]#[CX2]', kind: 'claim' },
  ];
  /* functional group → MolData.typicalLength environment (only where Allen 1987 has matching rows). A function entry
     picks the environment per bond from the match (query indices) and may add a stand-in note, which names the table
     entry actually shown: the C–O bonds of an acetal, hemiacetal or epoxide are C–O–C ether bonds (the ether entry
     for that carbon: dialkyl or aryl alkyl), the alkyl–O or aryl–O bond of a carbonate or carbamate is that bond of an
     ester (the ester entry for that carbon); their other bonds have no entry of their own. Amides: the C(=O)–N bond by
     the number of H on N (Allen 1987, acyclic amides: NH₂–C=O 1.325, C*–NH–C=O 1.334, (C*)₂N–C=O 1.346 Å); the
     N–C* entry is written for C*–NH–C=O, so a tertiary amide's N–C* bonds get none. Lactams: Allen's lactam values
     differ from the acyclic ones (β-lactam C(=O)–N 1.385 Å) and are not in the table, so a lactam's ring C(=O)–N and
     its N's other bonds get no amide entry, and a β-lactam's C=O none either. Phosphorus environments (phosphine,
     phosphine oxide, phosphate triester) come from P's neighbours (analyse: pEnv), ortho-substituted biaryl bonds from
     analyse: biarylEnv. */
  const STAND_IN = {       // note(t: the table row shown, aryl: the bond's carbon is an aromatic-ring carbon)
    ether: kind => t => `this ${kind} C–O bond is a C–O–C ether bond: the ether entry "${t.sub}" is shown (the table has no ${kind} entry)`,
    ester: kind => (t, aryl) => { const x = aryl ? 'aryl' : 'alkyl';
      return `the ${x}–O bond of this ${kind} is the ${x}–O–C(=O) bond of an ester: the ester entry "${t.sub}" is shown (the table has no ${kind} entry)`; },
  };
  const pick = (match, a, b, pairs) => {                // pairs: [[query i, query j, env, note]]
    for (const [i, j, env, note] of pairs) if ((match[i] === a && match[j] === b) || (match[i] === b && match[j] === a)) return env ? { env, note } : null;
    return null;
  };
  const FG_ENV = { 'carboxylic acid': 'acid', carboxylate: 'carboxylate', ester: 'ester', lactone: 'ester',
    'amide (2°)': 'amide', 'alcohol (1°)': 'alcohol',
    'alcohol (2°)': 'alcohol', 'alcohol (3°)': 'alcohol', 'alcohol (methanol)': 'alcohol', phenol: 'phenol', enol: 'enol',
    aldehyde: 'aldehyde', ketone: 'ketone', ether: 'ether', 'amine (1°)': 'amine', 'amine (2°)': 'amine', 'amine (3°)': 'amine',
    'aniline (aryl amine)': 'amine', 'ammonium (protonated amine)': 'ammonium', 'quaternary ammonium': 'ammonium',
    imine: 'imine', oxime: 'oxime', nitrile: 'nitrile', nitro: 'nitro', sulfide: 'sulfide',
    thiol: 'thiol', sulfone: 'sulfone', sulfoxide: 'sulfoxide', sulfonamide: 'sulfonamide',
    // [C](=[O])[NH2]: the C(=O)–N bond has its own entry; the N–H bonds (a === b) keep the general N–H entry
    'amide (1°)': (m, a, b) => (a === b ? { env: 'amide', note: null } : pick(m, a, b, [[0, 2, 'amide-primary', null], [0, 1, 'amide', null]])),
    // [C](=[O])[N](C)C: no entry for the N–C* bonds (a === b)
    'amide (3°)': (m, a, b) => (a === b ? null : pick(m, a, b, [[0, 2, 'amide-tertiary', null], [0, 1, 'amide', null]])),
    // [C](=[O])@[N]: only the C=O, and only outside a β-lactam (a four-membered ring through the C–N bond)
    lactam: (m, a, b, G) => (a !== b && ((m[0] === a && m[1] === b) || (m[0] === b && m[1] === a)) && !G.ring4 ? { env: 'amide', note: null } : null),
    // [C]-[O]-[C]-[O]-[C]
    'acetal / ketal': (m, a, b) => pick(m, a, b, [[0, 1], [1, 2], [2, 3], [3, 4]].map(([i, j]) => [i, j, 'ether', STAND_IN.ether('acetal')])),
    // [C]-[O]-[C]-[OH]: the C–OH bond (and its O–H) is an alcohol's
    hemiacetal: (m, a, b) => (a === b ? (a === m[3] ? { env: 'alcohol', note: null } : null)
      : pick(m, a, b, [[2, 3, 'alcohol', null], [0, 1, 'ether', STAND_IN.ether('hemiacetal')], [1, 2, 'ether', STAND_IN.ether('hemiacetal')]])),
    // [C]1-[O]-[C]1
    epoxide: (m, a, b) => pick(m, a, b, [[0, 1, 'ether', STAND_IN.ether('epoxide')], [1, 2, 'ether', STAND_IN.ether('epoxide')]]),
    // [C][O][C](=[O])[O][C]
    'carbonate ester': (m, a, b) => pick(m, a, b, [[0, 1, 'ester', STAND_IN.ester('carbonate')], [4, 5, 'ester', STAND_IN.ester('carbonate')]]),
    // [N][C](=[O])[O][C]
    carbamate: (m, a, b) => pick(m, a, b, [[3, 4, 'ester', STAND_IN.ester('carbamate')]]) };
  const envFor = (G, a, b) => { const e = FG_ENV[G.name];
    if (!e) return null;
    return typeof e === 'function' ? e(G.match, a, b, G) : { env: e, note: null }; };
  const ENV_CC = new Set(['nitrile']);                  // environments with C–C rows (C*–C≡N, Car–C≡N)
  // groups whose bonds have no row in Allen 1987 and must not borrow one written for another group: the azo N=N row
  // does not fit an azide or a diazo compound (terminal N=N 1.13 Å, not 1.24 Å). Sulfonic acids and sulfonates need
  // no entry here: they have no FG_ENV, so MolData uses only env-free rows or its covalent-radius estimate for them.
  const NO_LENGTH_ROW_NONE = new Set(['azide', 'diazo']);        // resonance N–N / C–N bonds: no row and no estimate
  TEXT.functionalGroups = `from the app's ${FG.length}-group pattern table (checked by tools/test_molinfo.js); groups it does not cover are not reported`;

  function substructMatches(M, q) {
    const r = JSON.parse(M.get_substruct_matches(q, JSON.stringify({ maxMatches: 1000 })));
    return Array.isArray(r) ? r.map(x => x.atoms) : [];   // "{}" when nothing matches
  }
  function perceiveGroups(R, M) {
    const out = [], claimed = new Map();                // atom → the group that claimed it
    for (const f of FG) {
      const q = R.get_qmol(f.smarts);
      if (!q) continue;
      let matches;
      try { matches = substructMatches(M, q); } finally { q.delete(); }
      const seen = new Set();
      for (const a of matches) {
        const key = [...a].sort((x, y) => x - y).join(','); if (seen.has(key)) continue; seen.add(key);
        const hl = (f.atomsToHighlight || a.map((_, i) => i)).map(i => a[i]);
        if (f.kind === 'overlay') { out.push({ name: f.name, atoms: hl, match: a.slice(), kind: 'overlay' }); continue; }
        const cl = (f.claims || a.map((_, i) => i)).map(i => a[i]);
        if (cl.some(x => claimed.has(x) && !(f.sharesWith && f.sharesWith.includes(claimed.get(x))))) continue;
        cl.forEach(x => { if (!claimed.has(x)) claimed.set(x, f.name); });
        out.push({ name: f.name, atoms: hl, match: a.slice(), kind: 'claim' });
      }
    }
    return out;
  }

  /* Chirality (41 reference cases in the fixtures): refuse for '*' / non-tetrahedral stereo; flag and refuse allenes,
     hindered biaryls and E double bonds in 8–11 rings; one unspecified centre → chiral; a molecule is chiral when its
     mirror image (every @ ↔ @@) has a different canonical SMILES, cross-checked with standard InChI; meso (IUPAC Gold
     Book: an achiral member of a set of stereoisomers that also has chiral members) if inverting a subset of the
     centres gives a chiral isomer (≤ 10 centres). */
  const mirrorSmiles = s => s.replace(/@@|@/g, t => (t === '@' ? '@@' : '@'));
  function canonPair(R, s) {
    const m = R.get_mol(s); if (!m) return null;
    try { return { smi: m.get_smiles(), inchi: m.get_inchi() }; } finally { m.delete(); }
  }
  function isChiralFully(R, smi) {
    const a = canonPair(R, smi), b = canonPair(R, mirrorSmiles(smi));
    if (!a || !b) return { chiral: false, agree: false };
    const bySmiles = a.smi !== b.smi, byInchi = a.inchi !== b.inchi;
    return { chiral: bySmiles, agree: bySmiles === byInchi };
  }
  /* Axial chirality through a ring, which has no tetrahedral stereocentre (so neither RDKit's labels nor the mirror
     test by SMILES and InChI can see it): a spiro atom joining two rings (2,6-disubstituted spiro[3.3]heptanes such as
     Fecht's acid, which has been resolved), or a ring atom with an exocyclic C=C to a carbon with two different
     substituents (4-methylcyclohexylideneacetic acid, resolved in 1909). Each end of the axis must be desymmetrised:
     a monocyclic even-membered ring whose two branches leaving the axis atom are alike (else that atom is an ordinary
     stereocentre) and whose opposite atom carries two different substituents (or an exocyclic C=C to such a carbon). */
  function ringAxisFlags(g) {
    const ranks = g.cipRanks, out = [];
    if (!ranks) return out;
    const ringCount = i => g.rings.filter(r => r.includes(i)).length;
    const twoDifferent = (x, notThese) => {          // two substituents outside `notThese` (implicit H counts), unlike
      const ex = g.adj[x].filter(e => !notThese.includes(e.nb)).map(e => e.nb), h = g.atoms[x].hs;
      if (ex.length + h !== 2 || h === 2) return false;
      return h === 1 ? true : ranks[ex[0]] !== ranks[ex[1]];
    };
    const desymmetrised = (r, s) => {
      if (r.length % 2 || r.some(x => x !== s && ringCount(x) !== 1)) return false;
      const k = r.indexOf(s), n1 = r[(k + 1) % r.length], n2 = r[(k + r.length - 1) % r.length], o = r[(k + r.length / 2) % r.length];
      if (ranks[n1] !== ranks[n2]) return false;
      const a = g.atoms[o];
      if (totalDegree(g, o) === 4 && sumBO(g, o) === 4 && (a.z === 6 || a.z === 14 || (a.z === 7 && a.chg === 1))) return twoDifferent(o, r);
      const ex = g.adj[o].find(e => !r.includes(e.nb) && g.bonds[e.bond].order === 2 && !g.bonds[e.bond].arom && g.atoms[e.nb].z === 6);
      return !!ex && totalDegree(g, ex.nb) === 3 && twoDifferent(ex.nb, [o]);
    };
    for (const s of spiroAtoms(g)) {
      const rs = g.rings.filter(r => r.includes(s));
      if (rs.length === 2 && desymmetrised(rs[0], s) && desymmetrised(rs[1], s)) { out.push('spiro compound: axial chirality possible, not assessed'); break; }
    }
    const alkylidene = g.bonds.some(b => b.order === 2 && !b.arom && [[b.a, b.b], [b.b, b.a]].some(([r0, x]) => {
      const rs = g.rings.filter(r => r.includes(r0));
      return rs.length === 1 && !rs[0].includes(x) && g.atoms[x].z === 6 && g.atoms[r0].z === 6 && totalDegree(g, x) === 3 &&
        twoDifferent(x, [r0]) && desymmetrised(rs[0], r0);
    }));
    if (alkylidene) out.push('ring with an exocyclic double bond (alkylidene ring): axial chirality possible, not assessed');
    return out;
  }
  function axialPlanarFlags(R, m, g) {
    const flags = [];
    const has = sm => { const q = R.get_qmol(sm); if (!q) return []; try { return substructMatches(m, q); } finally { q.delete(); } };
    if (has('[CX3;!H2]=[CX2]=[CX3;!H2]').length) flags.push('allene: axial chirality possible, not assessed');
    if (hinderedBiarylBonds(g).size) flags.push('hindered biaryl: atropisomerism (axial chirality) possible, not assessed');
    flags.push(...ringAxisFlags(g));
    const tags = JSON.parse(m.get_stereo_tags());
    for (const [a, b, lab] of tags.CIP_bonds) if (lab === '(E)' && g.rings.some(r => r.length >= 8 && r.length <= 11 && r.includes(a) && r.includes(b)))
      flags.push('E double bond in a medium ring: planar chirality possible, not assessed');
    return flags;
  }
  // SMILES atom tokens, in the order RDKit numbers the atoms of a molecule read from that SMILES
  const SMILES_ATOM = /\[[^\]]+\]|Br|Cl|[BCNOPSFIbcnops*]/g;
  // the canonical SMILES without the stereo tags of tautomeric P/S oxo centres (a phosphate P tag is not a configuration)
  function stripOxoTags(R, input) {
    const m = R.get_mol(input); if (!m) return input;
    let g, tags;
    try { g = fromRDKitJson(JSON.parse(m.get_json())); tags = JSON.parse(m.get_stereo_tags()); } finally { m.delete(); }
    const drop = new Set(tags.CIP_atoms.map(t => t[0]).filter(i => tautomericOxoCentre(g, i)));
    if (!drop.size) return input;
    let k = -1;
    const s2 = input.replace(SMILES_ATOM, tok => { k++; return drop.has(k) ? tok.replace(/@@?(TH\d|AL\d|SP\d|TB\d+|OH\d+)?/, '') : tok; });
    const c = canonPair(R, s2);
    return c ? c.smi : input;
  }
  /* opts.ringOpen: {atoms, isomers} from openRingStereo (unspecified cis/trans ring stereo). The open centres are
     counted from the same filtered list the record shows (no centre with two H: PH3, CH3PH2; no tautomeric P/S oxo
     centre: phosphates), and the tags of such P/S atoms are removed before the mirror tests. */
  function chirality(R, input, opts) {
    opts = opts || {};
    const m = R.get_mol(input); if (!m) return { verdict: 'not determined', reason: 'could not read the structure' };
    let tags, flags, g, fixedBridge;
    try {
      tags = JSON.parse(m.get_stereo_tags()); g = fromRDKitJson(JSON.parse(m.get_json())); flags = axialPlanarFlags(R, m, g);
      fixedBridge = fixedBridgeheads(R, m, g, tags);
    } finally { m.delete(); }
    // not stereocentres: tautomeric P/S oxo centres (their tags are removed below) and the bridgeheads of a small
    // bridged ring whose configuration is the only one possible (unspecified, so they have no tag to remove)
    const drop = new Set(tags.CIP_atoms.map(t => t[0]).filter(i => tautomericOxoCentre(g, i)));
    const skip = new Set([...drop, ...fixedBridge]);
    const smi = drop.size ? stripOxoTags(R, input) : input;
    if (/@(TH|AL|SP|TB|OH)/.test(smi) || /\*/.test(smi)) return { verdict: 'not determined', reason: 'non-tetrahedral stereo or an unspecified atom (*)' };
    const nUn = tags.CIP_atoms.filter(([i, lab]) => lab === '(?)' && !skip.has(i) && nH(g, i) < 2).length;
    const nTok = (smi.match(/@/g) || []).length - (smi.match(/@@/g) || []).length;
    const upper = tags.CIP_atoms.filter(t => /\((R|S)\)/.test(t[1]) && !drop.has(t[0])).length;
    const ring = opts.ringOpen || { atoms: [], isomers: null };
    if (flags.length) return { verdict: 'not determined', reason: flags.join('; ') };
    if (nUn > 0) {
      if (nUn === 1 && nTok === 0) return { verdict: 'chiral', reason: 'one stereocentre (configuration not given): chiral whichever it is' };
      return { verdict: 'not determined', reason: `${nUn} stereocentre${nUn > 1 ? 's' : ''} without a configuration` };
    }
    if (ring.atoms.length) {
      // open cis/trans ring stereo: the verdict must hold for every ring isomer, each tested by the mirror rule
      if (!ring.isomers) return { verdict: 'not determined', reason: `cis/trans ring stereo at ${ring.atoms.length} atoms not specified` };
      const res = ring.isomers.map(s => isChiralFully(R, stripOxoTags(R, s)));
      if (res.some(x => !x.agree)) return { verdict: 'not determined', reason: 'SMILES and InChI mirror tests disagree' };
      // (decalin's fusion atoms are chirality centres, so "no chirality centres" would be wrong there)
      if (res.every(x => !x.chiral)) return { verdict: 'achiral', reason: 'every cis/trans ring isomer (the input does not say which) is superimposable on its mirror image' };
      if (res.every(x => x.chiral)) return { verdict: 'chiral', reason: 'chiral whichever cis/trans ring isomer it is (the input does not say)' };
      return { verdict: 'not determined', reason: 'some cis/trans ring isomers are chiral and some are not; the input does not say which' };
    }
    if (nTok === 0) return { verdict: 'achiral', reason: 'no stereocentres' };
    const t = isChiralFully(R, smi);
    if (!t.agree) return { verdict: 'not determined', reason: 'SMILES and InChI mirror tests disagree' };
    if (t.chiral) return { verdict: 'chiral', reason: 'not superimposable on its mirror image' };
    if (nTok > 10) return { verdict: 'achiral', reason: 'superimposable on its mirror image (meso test skipped: more than 10 centres)' };
    const toks = []; smi.replace(/@@|@/g, (x, off) => { toks.push([off, x.length]); return x; });
    // inverting only one bridgehead of a small bridged ring gives an impossible trans bridge, not a stereoisomer
    // (specified norbornene would otherwise be called meso): such subsets are skipped
    const tied = [];
    { const g2 = (() => { const m2 = R.get_mol(smi); if (!m2) return null; try { return fromRDKitJson(JSON.parse(m2.get_json())); } finally { m2.delete(); } })();
      if (g2) {
        const atomAt = []; let k = -1;
        smi.replace(SMILES_ATOM, (tok, off) => { k++; if (/@/.test(tok)) atomAt.push([off, k]); return tok; });
        const tokAtom = toks.map(([off]) => { const hit = atomAt.find(([o]) => o < off && smi.indexOf(']', o) > off); return hit ? hit[1] : -1; });
        for (let x = 0; x < toks.length; x++) for (let y = x + 1; y < toks.length; y++)
          if (tokAtom[x] >= 0 && tokAtom[y] >= 0 && smallBridgeheads(g2, tokAtom[x], tokAtom[y])) tied.push([x, y]);
      } }
    for (let mask = 1; mask < (1 << toks.length); mask++) {
      if (tied.some(([x, y]) => ((mask >> x) & 1) !== ((mask >> y) & 1))) continue;
      let s = '', last = 0;
      toks.forEach(([off, len], k) => { s += smi.slice(last, off); const tok = smi.substr(off, len); s += (mask >> k) & 1 ? (tok === '@' ? '@@' : '@') : tok; last = off + len; });
      s += smi.slice(last);
      const c = canonPair(R, s); if (!c) continue;
      if (isChiralFully(R, c.smi).chiral) return { verdict: 'meso', reason: `superimposable on its mirror image although it has stereocentres (${upper} R/S)` };
    }
    return { verdict: 'achiral', reason: 'superimposable on its mirror image; no chiral stereoisomer exists' };
  }

  /* ================================================================ 11. analyse → record */
  const VAL_DOU = { C: 4, Si: 4, N: 3, P: 3, O: 2, S: 2, H: 1, F: 1, Cl: 1, Br: 1, I: 1 };
  const RO5_ELEMENTS = new Set(['H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Se', 'Br', 'I']);

  function countComponents(n, bonds) {
    const p = [...Array(n).keys()], f = x => (p[x] === x ? x : (p[x] = f(p[x])));
    bonds.forEach(b => { p[f(b.a)] = f(b.b); });
    const roots = new Map();
    for (let i = 0; i < n; i++) { const r = f(i); if (!roots.has(r)) roots.set(r, []); roots.get(r).push(i); }
    return [...roots.values()];
  }
  function atomMass(z, iso) {                           // for mass-weighted 3D values
    if (iso) return iso;                                 // isotope: its mass number (within 0.02 u of the true mass)
    const e = EL(z);
    if (e.aw != null) return e.aw;
    if (e.massNumber && e.massNumber.length) return e.massNumber[0];
    return null;
  }
  function weightStats(g, desc) {                         // CIAAW formula weight, else RDKit amw (labelled)
    const counts = {}; let isotopes = false, missing = null;
    for (const a of g.atoms) {
      const s = a.z ? SYM(a.z) : '*'; counts[s] = (counts[s] || 0) + 1;
      if (a.hs) counts.H = (counts.H || 0) + a.hs;
      if (a.iso) isotopes = true;
      if (!a.z || EL(a.z).aw == null) missing = s;
    }
    let mw, method, note = null, composition = null;
    let rounded;
    if (!isotopes && !missing) {
      // Exact decimal sum (integer units of 1e-9 u), rounded to 2 dp with ties to the even digit (ISO 80000-1):
      // C6H10O = 98.145 → 98.14, C13H18O2 = 206.285 → 206.28, C8H9NO2 = 151.165 → 151.16, as PubChem shows them.
      const nano = Object.entries(counts).reduce((s, [e, k]) => s + Math.round(MD().element(e).aw * 1e9) * k, 0);
      const cents = Math.floor(nano / 1e7), rem = nano - cents * 1e7;
      mw = nano / 1e9;
      rounded = (rem > 5e6 || (rem === 5e6 && cents % 2 === 1) ? cents + 1 : cents) / 100;
      method = 'CIAAW';
      composition = hillOrder(counts).map(e => ({ el: e, pct: 100 * counts[e] * MD().element(e).aw / mw }));
    } else {
      mw = desc.amw; method = 'RDKit'; rounded = Math.round(mw * 100) / 100;
      note = isotopes ? "RDKit's average molecular weight: isotope-labelled atoms use their isotope's mass"
                      : "RDKit's average molecular weight: " + (missing === '*' ? 'the structure has an unspecified atom (*)' : missing + ' has no standard atomic weight');
    }
    return { counts, mw, rounded, method, note, composition, isotopes };
  }

  function analyse(R, M, opts) {
    opts = opts || {};
    if (!R || !M) throw new Error('MolInfo.analyse needs the RDKit module and a molecule');
    const J = JSON.parse(M.get_json());
    const g = fromRDKitJson(J);
    const n = g.atoms.length;
    if (!n) throw new Error('empty molecule');
    const tags = JSON.parse(M.get_stereo_tags());
    const desc = JSON.parse(M.get_descriptors());

    // per-atom and per-bond chemistry on the H-suppressed graph (the form every port was validated on)
    const conj = rdkitConjugation(g), rdHyb = rdkitHybridization(g, conj), kek = kekuleStats(g, opts.kekuleCap || 200000);
    const ctx = { conj, conjShown: textbookConjugation(g), rdHyb, kek, atoms: [], biaryl: hinderedBiarylBonds(g) };
    for (let i = 0; i < n; i++) ctx.atoms[i] = atomTextbook(g, i, ctx);
    const gast = gasteiger(g, conj, rdHyb, 12);
    const tb = g.bonds.map((b, k) => bondTextbook(g, k, ctx));
    const ringOpen = openRingStereo(R, M, g, new Set(tags.CIP_atoms.map(t => t[0])));
    const stereo = stereoInfo(g, tags, ringOpen, fixedBridgeheads(R, M, g, tags));
    const groups = perceiveGroups(R, M);
    const resonant = resonanceBonds(g);

    // appended H: one per implicit H, grouped by parent in ascending parent index (RDKit add_hs order)
    const hParent = [];
    for (let i = 0; i < n; i++) for (let t = 0; t < g.atoms[i].hs; t++) hParent.push(i);
    const N = n + hParent.length, m0 = g.bonds.length;

    const ringsOf = i => [...new Set(g.rings.filter(r => r.includes(i)).map(r => r.length))].sort((x, y) => x - y);
    const bondRings = k => { const b = g.bonds[k]; const s = new Set();
      for (const r of g.rings) { const ia = r.indexOf(b.a), ib = r.indexOf(b.b);
        if (ia >= 0 && ib >= 0 && (Math.abs(ia - ib) === 1 || Math.abs(ia - ib) === r.length - 1)) s.add(r.length); }
      return [...s].sort((x, y) => x - y); };
    const groupsOfAtom = i => [...new Set(groups.filter(G => G.atoms.includes(i)).map(G => G.name))];
    // a β-lactam: the lactam's C(=O)–N bond lies in a four-membered ring (its C=O gets no amide entry)
    for (const G of groups) if (G.name === 'lactam') {
      const e = g.adj[G.match[0]].find(x => x.nb === G.match[2]);
      G.ring4 = !!e && bondRings(e.bond).includes(4);
    }
    // isotope labels: the isotope's mass from RDKit's isotope table (NIST); a label RDKit has no mass for keeps null
    const isoMass = new Map();
    for (const a of g.atoms) if (a.iso && a.z > 0 && !isoMass.has(a.z + ':' + a.iso)) {
      let m1 = null, v = null;
      try { m1 = R.get_mol(`[${a.iso}${SYM(a.z)}]`); v = m1 ? JSON.parse(m1.get_descriptors()).exactmw : null; } catch (e) { v = null; } finally { if (m1) m1.delete(); }
      // RDKit returns the bare mass number for an isotope it does not know ([99C] → 99); 12C is exactly 12 by definition
      isoMass.set(a.z + ':' + a.iso, v && Math.abs(v - a.iso) < 0.5 && (v !== a.iso || (a.z === 6 && a.iso === 12)) ? v : null);
    }

    /* Oxidation states that depend on where the drawing puts a delocalised charge (tropylium, Cp⁻, C₃H₃⁺, allyl and
       benzyl ions, enolates, pyridinium): atoms of a conjugated system that carries a formal charge, and that have a
       conjugated bond to an atom of the same element. The other atoms keep their value in every such drawing (an
       atom more electronegative than its partners always gets the bonding electrons: enolate O, pyridinium N). A charge
       paired with an opposite charge on a bonded atom (nitro, N-oxide, azide) does not count here: charge-separated
       systems are tested below by building their other drawings. */
    const osAsDrawn = new Set();
    {
      const parent = [...Array(n).keys()], f = x => (parent[x] === x ? x : (parent[x] = f(parent[x])));
      g.bonds.forEach((b, k) => { if (ctx.conjShown[k]) parent[f(b.a)] = f(b.b); });
      const lone = i => g.atoms[i].chg && !g.adj[i].some(e => g.atoms[e.nb].chg * g.atoms[i].chg < 0);
      const charged = new Set();
      for (let i = 0; i < n; i++) if (lone(i) && g.adj[i].some(e => ctx.conjShown[e.bond])) charged.add(f(i));
      for (let i = 0; i < n; i++) if (charged.has(f(i)) && g.adj[i].some(e => ctx.conjShown[e.bond] && g.atoms[e.nb].z === g.atoms[i].z)) osAsDrawn.add(i);

      /* Charge-separated resonance structures (1,3-dipoles: diazomethane C=N⁺=N⁻ ↔ C⁻–N⁺≡N, ozone, N₂O, azides,
         nitrile oxides, nitrones): a negative atom with a lone pair, bonded to a positive atom or in one conjugated
         system with one. Its lone pair is pushed into the π system along an alternating path (Z⁻–Y=X ↔ Z=Y–X⁻, and
         on along –W=V), which moves the negative charge without making a new one; each such drawing is built and the
         oxidation states of the atoms it changes are recomputed. An atom gets the caveat when its value changes; one
         whose value is the same in every drawing keeps none (nitro N and O, N₂O's O, azide's middle N, ozone's middle
         O: bonding electrons between atoms of one element are split, so moving them changes nothing there). Paths
         follow the drawn Kekulé orders, through aromatic rings too (pyridine N-oxide, pyridinium N-imides, sydnones). */
      if (!kek.capped) {
        const posComp = new Set();
        for (let i = 0; i < n; i++) if (g.atoms[i].chg > 0) posComp.add(f(i));
        const osNow = i => { const o = oxidationFrom(g.atoms[i].z, g.atoms[i].chg, neighbourOrders(g, i, kek), 'pauling'); return o ? o.value : null; };
        const osWith = (i, dq, dOrd) => {
          const nb = g.adj[i].map(e => ({ z: g.atoms[e.nb].z, order: kek.avgOrder[e.bond] + (dOrd.get(e.bond) || 0) }));
          for (let t = 0; t < g.atoms[i].hs; t++) nb.push({ z: 1, order: 1 });
          const o = oxidationFrom(g.atoms[i].z, g.atoms[i].chg + dq, nb, 'pauling');
          return o ? o.value : null;
        };
        const order = (k, dOrd) => g.bonds[k].order + (dOrd.get(k) || 0);
        const compare = (Z, end, dOrd) => {                // Z gave its lone pair (+1), end took one (−1)
          const touched = new Set([Z, end]);
          dOrd.forEach((d, k) => { touched.add(g.bonds[k].a); touched.add(g.bonds[k].b); });
          for (const i of touched) {
            const before = osNow(i), after = osWith(i, i === Z ? 1 : i === end ? -1 : 0, dOrd);
            if (before !== null && after !== null && Math.abs(before - after) > 1e-9) osAsDrawn.add(i);
          }
        };
        const push = (Z, u, dOrd, seen, depth) => {        // u holds the lone pair; form u=v, break v=w
          for (const e1 of g.adj[u]) {
            const v = e1.nb;
            if (seen.has(v) || order(e1.bond, dOrd) >= 3) continue;
            for (const e2 of g.adj[v]) {
              const w = e2.nb;
              if (w === u || seen.has(w) || order(e2.bond, dOrd) < 2) continue;
              const d2 = new Map(dOrd);
              d2.set(e1.bond, (d2.get(e1.bond) || 0) + 1); d2.set(e2.bond, (d2.get(e2.bond) || 0) - 1);
              compare(Z, w, d2);
              if (depth < 8) push(Z, w, d2, new Set([...seen, v, w]), depth + 1);
            }
          }
        };
        for (let i = 0; i < n; i++) {
          const a = g.atoms[i], L = ctx.atoms[i].lewis;
          if (a.chg >= 0 || !L || !(L.lp >= 1)) continue;
          if (!g.adj[i].some(e => g.atoms[e.nb].chg > 0) && !posComp.has(f(i))) continue;
          push(i, i, new Map(), new Set([i]), 0);
        }
      }
    }

    /* ---- atoms */
    const atoms = [];
    for (let i = 0; i < n; i++) {
      const a = g.atoms[i], T = ctx.atoms[i], L = T.lewis, e = EL(a.z), H = T.hybridization || {};
      // too many Kekulé structures to enumerate: no averaged value for atoms on aromatic bonds (never a partial average)
      const kekBlocked = kek.capped && g.adj[i].some(e => g.bonds[e.bond].arom);
      const os = kekBlocked ? null : oxidationFrom(a.z, a.chg, neighbourOrders(g, i, kek), 'pauling');
      atoms.push({
        index: i, el: e.sym, z: a.z, isH: a.z === 1, parent: -1, appended: false, charge: a.chg, isotope: a.iso || 0,
        radical: a.rad, hCount: nH(g, i), degree: totalDegree(g, i), neighbours: [], aromatic: a.arom, ringSizes: ringsOf(i),
        ringList: g.rings.filter(r => r.includes(i)).map(r => r.length).sort((x, y) => x - y),
        isotopeMass: a.iso ? isoMass.get(a.z + ':' + a.iso) : null,
        cip: stereo.atoms[i] ? stereo.atoms[i].label : null, cipText: stereo.atoms[i] ? stereo.atoms[i].text : null,
        groups: groupsOfAtom(i),
        valenceElectrons: L ? L.V : null, lonePairs: L ? L.lp : null, electronsAround: L && L.lp !== null ? L.shell : null,
        octetText: T.octet || (L && L.lp !== null ? (a.z <= 2 ? (L.shell === 2 ? 'duet' : 'not 2') : L.shell === 8 ? 'octet' : L.shell < 8 ? 'incomplete octet' : 'expanded octet (as drawn)') : null),
        stericNumber: T.stericNumber || null, effectiveStericNumber: T.effSN || null,
        domainGeometry: T.vsepr ? T.vsepr.electronDomain : null, molecularGeometry: T.vsepr ? T.vsepr.molecular : null,
        idealAngleText: T.vsepr ? T.vsepr.ideal : null, idealAngleNote: T.vsepr && T.vsepr.idealNote || null,
        ringAngleNote: T.vsepr && T.vsepr.ringNote || null,
        hybridization: hybridizationInfo(a, H), rdkitHybridization: rdHyb[i],
        oxidationState: kekBlocked ? { value: null, text: '—', flag: null, note: 'not computed: more than 200,000 Kekulé structures to average over', structures: null, borderline: [] }
          : withAsDrawn(oxidationInfo(a.z, a.chg, os, kek.structuresAt(i), () => oxidationFrom(a.z, a.chg, neighbourOrders(g, i, kek), 'allen'), nH(g, i)), osAsDrawn.has(i)),
        gasteiger: gast.ok ? gast.q[i] : null,
        en: e.en != null ? e.en : null, enAllen: e.enAllen != null ? e.enAllen : null, rcov: e.rcov != null ? e.rcov : null,
        rvdw: e.rvdw != null ? e.rvdw : null, mass: a.iso && isoMass.get(a.z + ':' + a.iso) ? isoMass.get(a.z + ':' + a.iso) : atomMass(a.z, a.iso), name: e.name,
      });
    }
    for (let t = 0; t < hParent.length; t++) {
      const p = hParent[t], zp = g.atoms[p].z, H = EL(1);
      const os = oxidationFrom(1, 0, [{ z: zp, order: 1 }], 'pauling');
      atoms.push({
        index: n + t, el: 'H', z: 1, isH: true, parent: p, appended: true, charge: 0, isotope: 0, radical: 0, hCount: 0, degree: 1,
        neighbours: [], aromatic: false, ringSizes: [], ringList: [], isotopeMass: null, cip: null, cipText: null, groups: [],
        valenceElectrons: 1, lonePairs: 0, electronsAround: 2, octetText: 'duet', stericNumber: null, effectiveStericNumber: null,
        domainGeometry: null, molecularGeometry: null, idealAngleText: null, idealAngleNote: null, ringAngleNote: null,
        hybridization: { label: null, qualifier: null, promoted: null, basis: 'hydrogen uses its 1s orbital: not hybridized',
          note: 'not hybridized: hydrogen bonds with its 1s orbital', hypervalent: false, terminal: false },
        rdkitHybridization: 'S',
        oxidationState: oxidationInfo(1, 0, os, 1, () => oxidationFrom(1, 0, [{ z: zp, order: 1 }], 'allen')),
        gasteiger: gast.ok ? gast.hq[p] / g.atoms[p].hs : null,
        en: H.en, enAllen: H.enAllen != null ? H.enAllen : null, rcov: H.rcov, rvdw: H.rvdw, mass: H.aw, name: H.name,
      });
    }

    /* ---- bonds */
    const bonds = [];
    const inGroup = (names, a, b) => groups.some(G => names.has(G.name) && G.match.includes(a) && G.match.includes(b));
    // phosphorus environments of the Allen tables, from P's own neighbours (no functional group names them):
    // phosphine R3P (C2–P–C*, C2–P–Car), phosphine oxide R3P=O (C2–P(=O)–C*, C3–P=O), phosphate triester (RO)3P=O
    const pEnv = p => {
      if (p >= n || g.atoms[p].z !== 15 || g.atoms[p].chg || g.atoms[p].rad || g.atoms[p].hs) return null;
      const single = g.adj[p].filter(e => g.bonds[e.bond].order === 1 && !g.bonds[e.bond].arom);
      const oxo = g.adj[p].filter(e => g.bonds[e.bond].order === 2 && !g.bonds[e.bond].arom && g.atoms[e.nb].z === 8 && degree(g, e.nb) === 1 && !g.atoms[e.nb].chg);
      const isC = e => g.atoms[e.nb].z === 6;
      const esterO = e => g.atoms[e.nb].z === 8 && !g.atoms[e.nb].chg && degree(g, e.nb) === 2 && g.adj[e.nb].some(f => f.nb !== p && g.atoms[f.nb].z === 6);
      if (g.adj[p].length === 3 && single.length === 3 && single.every(isC)) return 'phosphine';
      if (g.adj[p].length === 4 && oxo.length === 1 && single.length === 3 && single.every(isC)) return 'phosphine-oxide';
      if (g.adj[p].length === 4 && oxo.length === 1 && single.length === 3 && single.every(esterO)) return 'phosphate';
      return null;
    };
    // Car–Car bond between two rings with a non-H group on a ring atom next to the bond (a ring fusion counts): Allen
    // 1987 gives 1.490 Å (n = 212) for these and 1.487 Å (n = 30) with only H at the ortho positions
    const biarylEnv = (a, b) => {
      const e = g.adj[a].find(x => x.nb === b), bd = e && g.bonds[e.bond];
      if (!bd || bd.arom || bd.order !== 1 || inRingBond(g, e.bond) || lengthHyb(a, e.bond) !== 'ar' || lengthHyb(b, e.bond) !== 'ar') return null;
      const orthoGroup = (x, other) => g.adj[x].some(o => o.nb !== other && g.atoms[o.nb].arom && g.adj[o.nb].filter(f => g.atoms[f.nb].z !== 1).length >= 3);
      return orthoGroup(a, b) || orthoGroup(b, a) ? { env: 'biaryl-ortho', note: null } : null;
    };
    const envOf = (a, b) => {                            // environment tag for MolData.typicalLength
      const het = x => x < n && ![1, 6].includes(g.atoms[x].z);
      const pa = pEnv(a), pb = pEnv(b);
      if ((pa || pb) && a !== b) {
        // only the bonds those rows describe: C–P (phosphines, phosphine oxides) and P=O (phosphine oxides, phosphate
        // triesters); a phosphate's P–O and O–C bonds have no row
        const other = pa ? b : a, zo = g.atoms[other].z, bond = g.adj[a].find(e => e.nb === b);
        if (zo === 6 || (zo === 8 && bond && g.bonds[bond.bond].order === 2)) return { env: pa || pb, note: null };
      }
      const cc = !het(a) && !het(b);
      if (cc && a === b) return null;                    // C–H: no environment rows
      if (cc && a < n && b < n) { const bi = biarylEnv(a, b); if (bi) return bi; }
      for (const G of groups) {
        if (!G.match.includes(a) || !G.match.includes(b)) continue;
        const e = envFor(G, a, b);
        if (e && (!cc || ENV_CC.has(e.env))) return e;
      }
      if (cc) return null;
      // a bond from a group's heteroatom to an atom outside the match (the N–CH₃ of an N-methylamide)
      for (const G of groups) for (const h of [a, b]) if (het(h) && G.match.includes(h)) { const e = envFor(G, h, h); if (e) return e; }
      return null;
    };
    // carbon type for the Allen tables: 'ar' only for a benzenoid six-membered ring (Allen's Car). A ring whose atoms
    // carry an exocyclic double bond (2-pyridone, uracil, cytosine) is aromatic to RDKit but not benzene-like: its
    // bonds alternate (C5=C6 1.34 Å, C4–C5 1.43 Å), so it gets no 'ar' row. An aromatic bond gets the 'ar' rows (phenyl
    // C≃C, pyridine C≃N) only when every AROMATIC ring holding the bond is benzenoid: the five-membered-ring C–N of
    // indole, benzimidazole or purine at the ring fusion is 1.36–1.38 Å, not pyridine's 1.337. A non-aromatic ring
    // fused on (tetralin, chromane, indene, 9,10-dihydroanthracene) does not change the benzene ring's own bonds.
    const benzenoid = r => r.length === 6 && r.every(x => g.atoms[x].arom &&
      !g.adj[x].some(e => !g.bonds[e.bond].arom && g.bonds[e.bond].order >= 2));
    const aromaticRing = r => r.every((x, t) => { const e = g.adj[x].find(q => q.nb === r[(t + 1) % r.length]); return !!e && g.bonds[e.bond].arom; });
    const bondRingList = k => { const b = g.bonds[k];
      return g.rings.filter(r => { const ia = r.indexOf(b.a), ib = r.indexOf(b.b); return ia >= 0 && ib >= 0 && (Math.abs(ia - ib) === 1 || Math.abs(ia - ib) === r.length - 1); }); };
    // A ring bond joining two aryl rings (fluorene's C4a–C4b, in a five-membered ring) is not the free biaryl bond of
    // Allen's biphenyl row: 'none' is no carbon type, so MolData gives no value for it.
    const lengthHyb = (i, k) => {
      if (i >= n || g.atoms[i].z !== 6) return null;
      if (g.atoms[i].arom && g.rings.some(r => r.includes(i) && benzenoid(r))) {
        if (k === undefined || k >= g.bonds.length) return 'ar';
        const b = g.bonds[k];
        if (!b.arom) return inRingBond(g, k) && g.atoms[b.a].arom && g.atoms[b.b].arom ? 'none' : 'ar';
        const rl = bondRingList(k).filter(aromaticRing);
        return rl.length && rl.every(benzenoid) ? 'ar' : 'sp2';
      }
      if (g.atoms[i].arom) return 'sp2';
      const h = atoms[i].hybridization.label;
      return ['sp3', 'sp2', 'sp'].includes(h) ? h : null;
    };
    /* typical length: never a row written for another group. Bonds of groups without rows (azide, diazo) get none;
       a bond to a charged atom or a radical gets a value only from a row of its own group (carboxylate, nitro,
       ammonium), since phenoxide C–O⁻ or a nitroxide N–O• are far from the neutral-bond rows. A conjugated bond gets
       no covalent-radius estimate: resonance shortens or lengthens it by up to 0.12 Å (urea C–N 1.35 Å against a sum
       of 1.46 Å; the C5=C6 of uracil 1.34 Å against 1.42 Å).
       The table's sp² carbon rows for C–X single bonds (C=C–F, C=C–Cl, Br–Csp², and an ester's acyl C(=O)–O) describe
       a C=C or acyl carbon, so two bonds get none: the C–halogen bond of an acyl or imidoyl halide (C=O, C=N, C=S
       carbon; Allen 1987 has no entry for it), and the O–C bond from an ester oxygen to a C=C carbon (vinyl esters,
       enol carbonates: no entry either, and conjugated, so no estimate). */
    const HALOGEN = new Set([9, 17, 35, 53]);
    const heteroDouble = c => g.atoms[c].z === 6 && g.adj[c].some(e => !g.bonds[e.bond].arom && g.bonds[e.bond].order === 2 && [7, 8, 16].includes(g.atoms[e.nb].z));
    const typicalOf = (a, b, order, isH, conjugated, k) => {
      const A = atoms[a], B = atoms[b];
      if (!(A.z > 0 && B.z > 0)) return null;
      if (!isH && inGroup(NO_LENGTH_ROW_NONE, a, b)) return null;
      if (!isH && order === 1 && [[a, b], [b, a]].some(([c, x]) => HALOGEN.has(atoms[x].z) && heteroDouble(c))) return null;
      const e = isH ? envOf(a, a) : envOf(a, b), env = e ? e.env : null;
      // the carbon of a C–O bond (for the ester and stand-in rules)
      const cO = isH ? -1 : A.z === 6 && B.z === 8 ? a : B.z === 6 && A.z === 8 ? b : -1;
      if (env === 'ester' && order === 1 && cO >= 0 && lengthHyb(cO, k) === 'sp2' && !heteroDouble(cO)) return null;
      const t = MD().typicalLength(A.el, B.el, order, lengthHyb(a, k), isH ? null : lengthHyb(b, k), env);
      const special = A.charge || B.charge || A.radical || B.radical;
      if (special && !(t && !t.estimate && env && t.env === env)) return null;
      if (t && t.estimate && conjugated) return null;
      // a stand-in entry (acetal C–O shown with an ether entry) says so, naming the entry shown
      if (t && e && e.note && t.env === env) {
        const note = typeof e.note === 'function' ? e.note(t, cO >= 0 && lengthHyb(cO, k) === 'ar') : e.note;
        return Object.assign({}, t, { note: [note, t.note].filter(Boolean).join('; ') });
      }
      return t;
    };
    for (let k = 0; k < m0; k++) {
      const b = g.bonds[k], T = tb[k], A = atoms[b.a], B = atoms[b.b];
      const pol = T.polarity;
      const lengthOrder = b.arom ? 1.5 : b.order;
      bonds.push(bondRecord({
        index: k, a: b.a, b: b.b, order: b.order, aromatic: b.arom, T, pol, A, B,
        doubleIn: kek.capped ? null : kek.doubleIn[k], structures: kek.capped ? null : kek.structuresOf[k],
        enthalpy: b.arom || resonant.has(k) ? null : MD().bondEnthalpy(A.el, B.el, b.order), resonance: resonant.has(k),
        typical: typicalOf(b.a, b.b, lengthOrder, false, b.arom || T.rdkitConjugated, k),
        rot: T.rotatable, rotation: T.rotation, conjugated: T.conjugated, rdkitConjugated: T.rdkitConjugated, ringSizes: bondRings(k),
        ez: stereo.bonds[k] ? stereo.bonds[k].label : null, torsion: torsionAtoms(g, k),
        groups: [...new Set(groups.filter(G => G.atoms.includes(b.a) && G.atoms.includes(b.b)).map(G => G.name))],
      }));
    }
    for (let t = 0; t < hParent.length; t++) {
      const p = hParent[t], h = n + t, A = atoms[p], B = atoms[h], k = m0 + t;
      const ha = A.hybridization.label && !['sp3d', 'sp3d2'].includes(A.hybridization.label) ? hybText(A.hybridization.label) : null;
      const T = { type: 'single', kekule: null, bondOrder: 1, sigma: 1, pi: 0, piNote: null,
        sigmaOverlap: ha ? `σ: ${A.el}(${ha})–H(1s)` : null, piOverlap: null, polarity: polarity(A.z, 1), amideCN: false };
      bonds.push(bondRecord({
        index: k, a: p, b: h, order: 1, aromatic: false, T, pol: T.polarity, A, B, doubleIn: null, structures: null,
        enthalpy: MD().bondEnthalpy(A.el, 'H', 1), resonance: false,
        typical: typicalOf(p, h, 1, true),
        rot: { rotatable: false, rdkitStrict: false, code: 'hydrogen', reason: ROT_REASON.hydrogen },
        rotation: { free: null, text: TERMINAL_ROTATION },
        conjugated: false, rdkitConjugated: false, ringSizes: [], ez: null, torsion: null, groups: [], isHBond: true,
      }));
    }
    // neighbour lists in record numbering
    bonds.forEach(b => {
      atoms[b.a].neighbours.push({ atom: b.b, bond: b.index, order: b.order, aromatic: b.aromatic });
      atoms[b.b].neighbours.push({ atom: b.a, bond: b.index, order: b.order, aromatic: b.aromatic });
    });

    /* ---- molecule level */
    const W = weightStats(g, desc);
    let charge = 0, radicals = 0, chargedAtoms = 0, dummy = false;
    for (const a of g.atoms) { charge += a.chg; radicals += a.rad; if (a.chg) chargedAtoms++; if (a.z === 0) dummy = true; }
    const comps = countComponents(n, g.bonds);
    const formula = dummy ? null : hill(W.counts, charge);
    let dou = null, douNote = null;
    if (dummy) douNote = 'not defined: the structure has an unspecified atom (*)';
    else if (comps.length > 1) douNote = 'not defined for more than one molecule or ion';
    else if (charge !== 0) douNote = 'not defined for a net-charged species';
    else if (radicals) douNote = 'not defined for radicals';
    else if (Object.keys(W.counts).some(e => VAL_DOU[e] === undefined)) douNote = 'not defined: only for C, H, N, O, S, P, Si and halogens';
    else {
      dou = (2 + Object.entries(W.counts).reduce((s, [e, k]) => s + k * (VAL_DOU[e] - 2), 0)) / 2;
      // the formula assumes S and O divalent, P and N trivalent, halogens univalent: with an expanded octet (S=O in
      // DMSO, P=O, SF4, PCl5) its value is not the number of rings + π bonds, and can even be negative (SF6 −2)
      const expanded = ctx.atoms.some((T, i) => g.atoms[i].z > 10 && T.lewis && T.lewis.lp !== null && T.lewis.shell > 8);
      if (dou < 0 || expanded) { dou = null; douNote = 'not defined: an atom has more bonds than the formula assumes (hypervalent: drawn with more than 8 valence electrons)'; }
    }
    const douStructure = g.bonds.length - n + comps.length + g.bonds.reduce((s, b) => s + (b.order - 1), 0);
    if (dou !== null && douStructure !== dou) douNote = `the drawn structure has ${douStructure} ring${douStructure === 1 ? '' : 's'} + π bond${douStructure === 1 ? '' : 's'} in total`;
    const isomeric = M.get_smiles();
    // without stereo but with the isotope labels (RDKit's doIsomericSmiles:false would drop both): the stereo marks
    // are removed from the isomeric SMILES and the result is canonicalised again
    let noStereo = isomeric;
    { const s0 = isomeric.replace(/\[[^\]]+\]/g, t => t.replace(/@@?(TH\d|AL\d|SP\d|TB\d+|OH\d+)?/, '')).replace(/[/\\]/g, '');
      if (s0 !== isomeric) { const c = canonPair(R, s0); noStereo = c ? c.smi : null; } }
    const smiles = { canonical: M.get_smiles(JSON.stringify({ doIsomericSmiles: false })), isomeric, noStereo };
    let inchi = null, inchikey = null;
    try { inchi = M.get_inchi(); inchikey = inchi ? R.get_inchikey_for_inchi(inchi) : null; } catch (e) { inchi = null; }
    const rotCount = bonds.reduce((s, b) => s + (b.rotatable.yes ? 1 : 0), 0);
    const elementsSet = Object.keys(W.counts);
    const organic = comps.length === 1 && W.counts.C > 0 && elementsSet.every(e => RO5_ELEMENTS.has(e));
    const ro5v = [];
    if (W.mw > 500) ro5v.push('MW > 500'); if (desc.CrippenClogP > 5) ro5v.push('logP > 5');
    if (desc.lipinskiHBD > 5) ro5v.push('H donors > 5'); if (desc.lipinskiHBA > 10) ro5v.push('H acceptors > 10');
    const stereoAtoms = [], stereoBonds = [];
    stereo.atoms.forEach((s, i) => { if (s) stereoAtoms.push(Object.assign({ index: i, label: s.label, kind: s.kind }, s.ring ? { ring: s.ring } : {})); });
    stereo.bonds.forEach((s, k) => { if (s) stereoBonds.push({ index: k, label: s.label }); });
    const chir = chirality(R, smiles.isomeric, { ringOpen });

    const record = {
      n, N, atoms, bonds, rings: g.rings.map(r => r.slice()), aromaticSystems: aromaticSystems(g),
      groups: groups.map(G => ({ name: G.name, atoms: G.atoms, match: G.match, kind: G.kind })),
      // stereo.ignored: centres RDKit proposes that are not stereocentres (two H; tautomeric P/S oxo). A view that
      // lists open stereo from another source (the 3D engine) should leave these atoms out too.
      stereo: { atoms: stereoAtoms, bonds: stereoBonds, ignored: stereo.ignored }, chirality: chir,
      smiles, inchi, inchikey, formula, formulaText: formula ? hillText(W.counts, charge) : null,
      isotopicFormula: W.isotopes ? isotopicFormula(g.atoms) + chargeSup(charge) : null,
      charge, chargedAtoms, radicals, components: comps.length,
      mw: { value: W.mw, rounded: W.rounded, method: W.method, note: W.note }, composition: W.composition,
      monoisotopic: desc.exactmw, heavyAtoms: desc.NumHeavyAtoms, totalAtoms: desc.NumAtoms,
      dou: { value: dou, note: douNote, structure: douStructure },
      ringCounts: { total: desc.NumRings, aromatic: desc.NumAromaticRings, aliphatic: desc.NumAliphaticRings,
        saturated: desc.NumSaturatedRings, heterocycles: desc.NumHeterocycles, aromaticHeterocycles: desc.NumAromaticHeterocycles,
        spiro: spiroAtoms(g).length, spiroRDKit: desc.NumSpiroAtoms, bridgehead: desc.NumBridgeheadAtoms },
      rotatableBonds: { count: rotCount, rdkit: desc.NumRotatableBonds },
      hbd: desc.lipinskiHBD, hba: desc.lipinskiHBA, tpsa: desc.tpsa, logp: desc.CrippenClogP, mr: desc.CrippenMR,
      fsp3: W.counts.C ? desc.FractionCSP3 : null,
      ro5: organic ? { violations: ro5v, pass: ro5v.length <= 1 } : null,
      veber: organic ? { pass: desc.NumRotatableBonds <= 10 && desc.tpsa <= 140 } : null,
      gasteigerNote: gast.ok ? null : 'not available: no Gasteiger parameters for ' + (gast.missing.join(', ') || 'this structure'),
      // S=O / P=O present in either drawing (double bond, or S⁺–O⁻ / P⁺–O⁻): Gasteiger charges depend on the drawing
      hasHypervalentSP: bonds.some(b => { const x = atoms[b.a], y = atoms[b.b], sp = t => t.el === 'S' || t.el === 'P';
        const [c, o] = sp(x) && y.el === 'O' ? [x, y] : sp(y) && x.el === 'O' ? [y, x] : [null, null];
        return !!c && (b.order >= 2 || (c.charge > 0 && o.charge < 0)); }),
      descriptors: desc, kekuleCapped: kek.capped,
    };
    return record;
  }

  /* Spiro atoms by the IUPAC definition: the only atom common to two rings, i.e. removing it splits its ring system
     into parts with no ring bond between them (spiro[4.5]decane C5). RDKit's NumSpiroAtoms is SSSR-based and also
     counts atoms where two SSSR rings of one fused system touch (morphine C13). */
  function spiroAtoms(g) {
    const rb = ringBondSet(g), out = [];
    for (let i = 0; i < g.atoms.length; i++) {
      const rn = g.adj[i].filter(e => rb.has(e.bond)).map(e => e.nb);
      if (rn.length < 4) continue;
      const seen = new Set([i]); let parts = 0;
      for (const s0 of rn) {
        if (seen.has(s0)) continue;
        parts++; const st = [s0]; seen.add(s0);
        while (st.length) { const u = st.pop(); for (const e of g.adj[u]) if (rb.has(e.bond) && !seen.has(e.nb)) { seen.add(e.nb); st.push(e.nb); } }
      }
      if (parts >= 2) out.push(i);
    }
    return out;
  }
  function hybridizationInfo(a, H) {
    const label = H.label || null, el = a.z ? SYM(a.z) : '*';
    let qualifier = H.conventional ? 'conventionally' : null, note = H.note || null;
    if (label) {
      if (H.labelOnly) note = hybText(label) + ' is a textbook label only: d orbitals contribute little to the bonding';
      else if (H.hypervalent) {
        // name the bonds this atom really has (S=O, P=O, P=C of an ylide, P=S …), singular or plural by count
        const P = H.partners || [], kinds = [...new Set(P)];
        const REF = ', with no d orbitals involved (Reed & Schleyer, J. Am. Chem. Soc. 1990, 112, 1434)';
        if (!P.length) note = 'drawn with more than 8 electrons: the bonding is better described by polar bonds' + REF;
        else {
          const one = P.length === 1, drawn = kinds.map(x => `${el}=${x}`).join(' and '), polar = kinds.map(x => `${el}⁺–${x}⁻`).join(' and ');
          note = `drawn with more than 8 electrons: ${one ? 'the' : 'the ' + countWord(P.length)} ${drawn} bond${one ? ' is' : 's are'} better ` +
            `described as ${one ? 'a polar' : 'polar'} ${polar} single bond${one ? '' : 's'}` + (kinds.includes('C') && el === 'P' ? ' (an ylide)' : '') + REF;
        }
      }
      else if (H.promoted === 'aromatic') note = 'the lone pair is part of the aromatic π system, so it sits in a p orbital';
      else if (H.promoted === 'amide') note = 'amide N: planar; its lone pair is delocalised into the C=O group';
      else if (H.promoted === 'conjugated') note = 'the lone pair is conjugated with the neighbouring π system; the real atom lies ' +
        'between sp³ and sp² (slightly pyramidal)';
      else if (H.promoted === 'cation') note = 'the lone pair is shared with the neighbouring empty p orbital (C⁺)';
      else if (H.terminal) note = 'terminal atom: the label is a convention (there is no bond angle to test it)';
      else if (!note) note = 'from the steric number (σ bonds + lone pairs)';
    }
    return { label, qualifier, promoted: H.promoted || null, basis: H.promoted || H.note || null, note,
      hypervalent: !!H.hypervalent, terminal: !!H.terminal };
  }
  function oxidationInfo(z, chg, os, structures, allenFn, nHbonds) {   // nHbonds: bonds from this atom to H
    if (!os) return { value: null, text: '—', flag: null, note: 'no electronegativity for this element or a neighbour', structures, borderline: [] };
    const v = os.value;
    const text = Number.isInteger(v) ? intSigned(v) : (v > 0 ? '+' : MINUS) + trimNum(Math.abs(v), 2);
    const notes = [];
    if (!Number.isInteger(v)) notes.push(`average over ${structures} Kekulé structures`);
    let flag = null;
    if (os.borderline.length) {
      const al = allenFn();
      const parts = os.borderline.map(bd => {
        const A = SYM(z), X = bd.with, eA = EN(z), eX = EN(bd.z);
        if (bd.allenAgrees === false) return `the electronegativity scales disagree for ${A}–${X}: Pauling ${A} ${eA.toFixed(2)}, ${X} ${eX.toFixed(2)}; ` +
          `Allen ${A} ${ALLEN(z).toFixed(3)}, ${X} ${ALLEN(bd.z).toFixed(3)}`;
        let s = `${A} and ${X} have almost the same Pauling electronegativity (${eA.toFixed(2)} vs ${eX.toFixed(2)})`;
        if (bd.allenAgrees) s += '; the Allen scale orders them the same way, so both scales give this value';
        if ((z === 1 || bd.z === 1) && [5, 14, 15, 32, 33, 51, 52].includes(z === 1 ? bd.z : z) && EN(z === 1 ? bd.z : z) < EN(1)) {
          // the common textbook rule counts H as +1 with any non-metal: each such H bond moves 2 units (P in PH3: +3 → −3)
          const alt = z === 1 ? 1 : Number.isInteger(v) && nHbonds ? v - 2 * nHbonds : null;
          s += "; most textbooks instead treat H as +1 when it is bonded to a non-metal" + (alt !== null ? `, which gives ${intSigned(alt)} here` : '');
        }
        return s;
      });
      flag = 'borderline: ' + parts.join('; ');
      if (al && al.value !== v) flag += `. IUPAC's 2016 definition uses Allen electronegativities, which give ${Number.isInteger(al.value) ? intSigned(al.value) : trimNum(al.value, 2)}`;
    }
    return { value: v, text, flag, note: notes.join('; ') || null, structures,
      borderline: os.borderline.map(bd => ({ with: bd.with, dPauling: bd.dPauling, allenAgrees: bd.allenAgrees })) };
  }
  const withAsDrawn = (os, yes) => (yes && os.value !== null ? Object.assign(os, { note: [os.note, TEXT.osAsDrawn].filter(Boolean).join('; ') }) : os);
  function bondRecord(o) {
    const { T, pol, A, B } = o;
    const b = {
      index: o.index, a: o.a, b: o.b, order: o.order, aromatic: o.aromatic,
      typeText: T.type, kekule: T.kekule, bondOrder: T.bondOrder,
      orderText: o.aromatic ? (T.bondOrder === 1.5 ? '1.5' : 'between 1 and 2') : String(o.order),
      sigma: T.sigma, pi: T.pi,
      sigmaPiText: o.aromatic ? "1 σ + a share of the ring's delocalised π system" : ['1 σ', '1 σ + 1 π', '1 σ + 2 π'][o.order - 1] || '—',
      overlap: [T.sigmaOverlap, T.piOverlap].filter(Boolean), piNote: T.piNote || null,
      kekuleDouble: o.aromatic && o.structures ? { times: o.doubleIn, of: o.structures } : null,
      dEN: pol ? pol.dEN : null,
      polarity: pol ? polarityInfo(pol, o.a, o.b, A, B) : null,
      enthalpy: o.enthalpy ? { kJ: o.enthalpy.kJ, source: o.enthalpy.source, note: o.enthalpy.note } : null,
      resonance: !!o.resonance,      // equivalent by resonance to another bond of its group: no single average enthalpy
      typical: o.typical ? { A: o.typical.A, sd: o.typical.sd, n: o.typical.n, source: o.typical.source, sub: o.typical.sub,
        note: o.typical.note, estimate: !!o.typical.estimate, env: o.typical.env || null } : null,
      // rotatable = RDKit's strict descriptor (counted or not, with the rule's reason); rotation = the textbook answer
      rotatable: { yes: o.rot.rotatable, reason: o.rot.reason, rdkitStrict: o.rot.rdkitStrict, code: o.rot.code || null },
      rotation: o.rotation ? { free: o.rotation.free, text: o.rotation.text } : null,
      conjugated: !!o.conjugated, rdkitConjugated: !!o.rdkitConjugated, inRing: o.ringSizes.length > 0, ringSizes: o.ringSizes, ez: o.ez,
      amideNote: T.amideCN ? 'partial double-bond character (resonance): rotation is restricted, slow at room temperature' : null,
      torsionAtoms: o.torsion, groups: o.groups, isHBond: !!o.isHBond,
    };
    return b;
  }
  function polarityInfo(pol, ia, ib, A, B) {
    const neg = pol.negativeEnd === 'a' ? ia : pol.negativeEnd === 'b' ? ib : null;
    const pos = neg === null ? null : neg === ia ? ib : ia;
    const nm = i => (i === ia ? A : B).el + (i + 1);
    let text;
    if (pol.cls === 'largely ionic') text = `largely ionic (ΔEN ${pol.dEN.toFixed(2)})`;
    else if (pol.cls === 'polar covalent') text = `polar covalent (ΔEN ${pol.dEN.toFixed(2)}): δ+ on ${nm(pos)}, δ− on ${nm(neg)}`;
    else if (pol.dEN === 0) text = neg === null ? 'nonpolar covalent (same electronegativity)' : 'nonpolar covalent';
    else {
      const ch = [A.el, B.el].sort().join('') === 'CH';
      text = `nonpolar covalent${ch ? ' (conventionally)' : ''}: ΔEN ${pol.dEN.toFixed(2)} is very small` +
        (neg !== null ? ` (${nm(neg)} slightly δ−)` : '');
    }
    return { cls: pol.cls, text, positiveEnd: pos, negativeEnd: neg };
  }

  /* ================================================================ 12. geometry from 3D coordinates */
  function mulberry32(seed) {                             // small seeded PRNG (int32 state keeps it fast)
    let s = seed | 0;
    const imul = Math.imul;                               // local: global lookups are slow in some hosts (Node vm)
    return () => {
      s = (s + 0x6D2B79F5) | 0;
      let t = imul(s ^ (s >>> 15), 1 | s);
      t = (t + imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }
  function jacobiEigen3(A) {                              // eigenvalues of a symmetric 3×3 matrix, ascending
    const a = [A[0].slice(), A[1].slice(), A[2].slice()];
    for (let sweep = 0; sweep < 50; sweep++) {
      const off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
      if (off < 1e-22) break;
      for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < 3; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
        for (let k = 0; k < 3; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
      }
    }
    return [a[0][0], a[1][1], a[2][2]].sort((x, y) => x - y);
  }
  /* van der Waals volume by Monte Carlo inside each sphere: a random point of sphere i counts when no sphere j < i
     contains it, so V = Σ Vi·(fraction counted). Sampling only inside the spheres keeps the error far below that of a
     bounding-box estimate for long molecules; the variance is summed per atom for the standard error. Seeded. */
  function mcVolume(X, r, perAtom, seed) {
    const n = r.length, rnd = mulberry32(seed);
    let V = 0, varV = 0;
    for (let i = 0; i < n; i++) {
      const nb = [];
      for (let j = 0; j < i; j++) {
        const dx = X[3 * i] - X[3 * j], dy = X[3 * i + 1] - X[3 * j + 1], dz = X[3 * i + 2] - X[3 * j + 2];
        if (dx * dx + dy * dy + dz * dz < (r[i] + r[j]) * (r[i] + r[j])) nb.push(j);
      }
      const Vi = 4 / 3 * Math.PI * r[i] * r[i] * r[i];
      if (!nb.length) { V += Vi; continue; }
      let kept = 0;
      for (let s = 0; s < perAtom; s++) {
        let ux, uy, uz, q;
        do { ux = 2 * rnd() - 1; uy = 2 * rnd() - 1; uz = 2 * rnd() - 1; q = ux * ux + uy * uy + uz * uz; } while (q > 1);
        const px = X[3 * i] + r[i] * ux, py = X[3 * i + 1] + r[i] * uy, pz = X[3 * i + 2] + r[i] * uz;
        let covered = false;
        for (let t = 0; t < nb.length; t++) {
          const j = nb[t], dx = px - X[3 * j], dy = py - X[3 * j + 1], dz = pz - X[3 * j + 2];
          if (dx * dx + dy * dy + dz * dz <= r[j] * r[j]) { covered = true; break; }
        }
        if (!covered) kept++;
      }
      const f = kept / perAtom;
      V += Vi * f; varV += Vi * Vi * f * (1 - f) / perAtom;
    }
    return { value: V, stdErr: Math.sqrt(varV), points: perAtom * n };
  }
  const spiralCache = new Map();
  function goldenSpiral(n) {                             // golden-section spiral: n near-uniform points on a sphere
    if (spiralCache.has(n)) return spiralCache.get(n);
    const P = new Float64Array(3 * n), ga = Math.PI * (1 + Math.sqrt(5));
    for (let k = 0; k < n; k++) {
      const i = k + 0.5, phi = Math.acos(1 - 2 * i / n), th = ga * i;
      P[3 * k] = Math.cos(th) * Math.sin(phi); P[3 * k + 1] = Math.sin(th) * Math.sin(phi); P[3 * k + 2] = Math.cos(phi);
    }
    spiralCache.set(n, P);
    return P;
  }
  function vdwArea(X, r, nPts) {                          // Shrake–Rupley, probe 0
    const n = r.length, P = goldenSpiral(nPts);
    let tot = 0;
    for (let i = 0; i < n; i++) {
      const nb = [];
      for (let j = 0; j < n; j++) {
        if (j === i) continue;
        const dx = X[3 * i] - X[3 * j], dy = X[3 * i + 1] - X[3 * j + 1], dz = X[3 * i + 2] - X[3 * j + 2];
        if (dx * dx + dy * dy + dz * dz < (r[i] + r[j]) * (r[i] + r[j])) nb.push(j);
      }
      let free = 0, last = 0;
      for (let k = 0; k < nPts; k++) {
        const px = X[3 * i] + r[i] * P[3 * k], py = X[3 * i + 1] + r[i] * P[3 * k + 1], pz = X[3 * i + 2] + r[i] * P[3 * k + 2];
        let buried = false;
        for (let t = 0; t < nb.length; t++) {
          const j = nb[(t + last) % nb.length], dx = px - X[3 * j], dy = py - X[3 * j + 1], dz = pz - X[3 * j + 2];
          if (dx * dx + dy * dy + dz * dz <= r[j] * r[j]) { buried = true; last = (t + last) % nb.length; break; }
        }
        if (!buried) free++;
      }
      tot += 4 * Math.PI * r[i] * r[i] * free / nPts;
    }
    return tot;
  }
  const AREA_POINTS = 5000, MC_PER_ATOM = 4000, MC_MIN_POINTS = 200000, MC_SEED = 20260928;
  const geoCache = new WeakMap();

  function geometry(record, coords, opts) {
    opts = opts || {};
    const N = record.N;
    if (!coords || coords.length < 3 * N) throw new Error('MolInfo.geometry: expected ' + 3 * N + ' coordinates in record order');
    const optimised = opts.optimised !== false;
    let sum = 0; for (let i = 0; i < 3 * N; i++) sum += coords[i] * (i % 7 + 1);
    const cached = geoCache.get(record);
    const key = [coords, sum, optimised, opts.flaggedBonds, opts.radii, opts.mcPerAtom, opts.areaPoints];
    if (cached && cached.key.every((v, t) => v === key[t])) return cached.geo;
    const flagged = new Map();
    for (const f of opts.flaggedBonds || []) {
      const k = record.bonds.findIndex(b => (b.a === f[0] && b.b === f[1]) || (b.a === f[1] && b.b === f[0]));
      if (k >= 0) flagged.set(k, f[2] || 'flagged by the force-field check');
    }
    // atoms of flagged bonds: an angle or dihedral that involves one is marked; a flagged bond between two heavy atoms
    // (the engine checks only those) withholds every molecule-level value, since the whole model can be wrong there
    const flaggedAtoms = new Set();
    flagged.forEach((t, k) => { flaggedAtoms.add(record.bonds[k].a); flaggedAtoms.add(record.bonds[k].b); });
    const heavy = [...flagged.keys()].filter(k => !record.atoms[record.bonds[k].a].isH && !record.atoms[record.bonds[k].b].isH);
    const P = i => [coords[3 * i], coords[3 * i + 1], coords[3 * i + 2]];
    const dist = (i, j) => Math.hypot(coords[3 * i] - coords[3 * j], coords[3 * i + 1] - coords[3 * j + 1], coords[3 * i + 2] - coords[3 * j + 2]);
    let geo;
    if (!optimised) {
      // no measured numbers, except lengths of bonds the engine explicitly flagged (their text says why)
      const bl = new Float64Array(record.bonds.length).fill(NaN);
      flagged.forEach((t, k) => { bl[k] = dist(record.bonds[k].a, record.bonds[k].b); });
      geo = { optimised: false, bondLength: bl, flagged, flaggedAtoms, angles: null, dihedral: null, molecule: null };
    } else {
      const bl = new Float64Array(record.bonds.length);
      record.bonds.forEach((b, k) => { bl[k] = dist(b.a, b.b); });
      const angles = record.atoms.map((a, c) => {
        const out = [], nb = a.neighbours;
        for (let x = 0; x < nb.length; x++) for (let y = x + 1; y < nb.length; y++) {
          const t = { i: nb[x].atom, j: nb[y].atom, deg: angle(P(nb[x].atom), P(c), P(nb[y].atom)) };
          if ([t.i, c, t.j].some(q => flaggedAtoms.has(q))) t.unreliable = true;
          out.push(t);
        }
        return out;
      });
      const dih = record.bonds.map(b => {
        if (!b.torsionAtoms) return null;
        const [p, q, r, s] = b.torsionAtoms, deg = dihedral(P(p), P(q), P(r), P(s));
        const kp = klynePrelog(deg, record.atoms[q].hybridization.label === 'sp3' && record.atoms[r].hybridization.label === 'sp3');
        const d = { deg, name: kp.name, code: kp.code, common: kp.common, atoms: b.torsionAtoms.slice() };
        if (b.torsionAtoms.some(x => flaggedAtoms.has(x))) d.unreliable = true;
        return d;
      });
      const withheld = heavy.length ? (() => {
        const names = heavy.map(k => `${atomName(record, record.bonds[k].a)}–${atomName(record, record.bonds[k].b)}`);
        return { bonds: heavy.slice(), names, text: TEXT.flaggedValues(names),
          note: TEXT.flaggedValuesNote + ' ' + heavy.map((k, t) => `${names[t]}: ${flagged.get(k)}`).join('. ') + '.' };
      })() : null;
      geo = { optimised: true, bondLength: bl, flagged, flaggedAtoms, angles, dihedral: dih, molecule: moleculeGeometry(record, coords, opts, withheld) };
    }
    geoCache.set(record, { key, geo });
    return geo;
  }

  /* linear / planar: largest distance of an atom from the least-squares line and plane through all atoms (H included) */
  const LINEAR_TOL = 0.05, PLANAR_TOL = 0.1;
  function flatness(X, N) {
    if (N < 2) return { form: null, dev: null };
    const c = [0, 0, 0];
    for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) c[d] += X[3 * i + d] / N;
    const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let i = 0; i < N; i++) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) C[a][b] += (X[3 * i + a] - c[a]) * (X[3 * i + b] - c[b]);
    const E = eigenSym3(C), n = E[0].vector, u = E[2].vector;
    let plane = 0, line = 0;
    for (let i = 0; i < N; i++) {
      const v = [0, 1, 2].map(d => X[3 * i + d] - c[d]), along = v[0] * u[0] + v[1] * u[1] + v[2] * u[2];
      plane = Math.max(plane, Math.abs(v[0] * n[0] + v[1] * n[1] + v[2] * n[2]));
      line = Math.max(line, Math.hypot(v[0] - along * u[0], v[1] - along * u[1], v[2] - along * u[2]));
    }
    return { form: line <= LINEAR_TOL ? 'linear' : plane <= PLANAR_TOL ? 'planar' : '3d', dev: { line, plane } };
  }

  function moleculeGeometry(record, X, opts, withheld) {
    const N = record.N, masses = record.atoms.map(a => a.mass);
    const out = {};
    if (withheld) {
      // a flagged heavy-atom bond: no size, shape or symmetry value from this model (same keys, all null)
      Object.assign(out, { withheld, centre: null, rg: null, pmi: null, npr1: null, npr2: null, shape: null, shapeScores: null,
        form: null, formDev: null, extent: null, vdwVolume: null, vdwVolumeStdErr: null, vdwArea: null, radiiNote: null,
        dipole: { verdict: 'not assessed', text: TEXT.dipoleVerdict.flagged(withheld.names), note: withheld.note, symmetry: null } });
      return out;
    }
    out.withheld = null;
    if (masses.every(m => m != null && m > 0)) {
      const W = masses.reduce((s, m) => s + m, 0), com = [0, 0, 0];
      for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) com[d] += masses[i] * X[3 * i + d] / W;
      let rg = 0; const I = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
      for (let i = 0; i < N; i++) {
        const x = X[3 * i] - com[0], y = X[3 * i + 1] - com[1], z = X[3 * i + 2] - com[2], m = masses[i];
        rg += m * (x * x + y * y + z * z);
        I[0][0] += m * (y * y + z * z); I[1][1] += m * (x * x + z * z); I[2][2] += m * (x * x + y * y);
        I[0][1] -= m * x * y; I[0][2] -= m * x * z; I[1][2] -= m * y * z;
      }
      I[1][0] = I[0][1]; I[2][0] = I[0][2]; I[2][1] = I[1][2];
      out.centre = com;
      out.rg = Math.sqrt(rg / W);
      const ev = jacobiEigen3(I);
      out.pmi = ev.map(v => Math.max(0, v));
      if (out.pmi[2] > 1e-6) {
        out.npr1 = out.pmi[0] / out.pmi[2]; out.npr2 = out.pmi[1] / out.pmi[2];
        const sc = { rod: out.npr2 - out.npr1, disc: 2 - 2 * out.npr2, sphere: out.npr1 + out.npr2 - 1 };
        out.shapeScores = sc;
        out.shape = Object.keys(sc).reduce((best, k) => (sc[k] > sc[best] ? k : best), 'rod');
      } else { out.npr1 = null; out.npr2 = null; out.shape = null; out.shapeScores = null; }
    } else { out.rg = null; out.pmi = null; out.npr1 = null; out.npr2 = null; out.shape = null; out.centre = null; }
    const fl = flatness(X, N);
    out.form = fl.form; out.formDev = fl.dev;
    const radii = opts.radii || record.atoms.map(a => a.rvdw);
    const haveR = radii.every(v => v != null && v > 0);
    let dMax = 0, eMax = 0;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) {
      const d = Math.hypot(X[3 * i] - X[3 * j], X[3 * i + 1] - X[3 * j + 1], X[3 * i + 2] - X[3 * j + 2]);
      if (d > dMax) dMax = d;
      if (haveR && d + radii[i] + radii[j] > eMax) eMax = d + radii[i] + radii[j];
    }
    if (N === 1 && haveR) eMax = 2 * radii[0];
    out.extent = { nuclei: dMax, vdw: haveR ? eMax : null };
    // van der Waals volume and area: computed on first use, then cached
    let vol, area;
    Object.defineProperty(out, 'vdwVolume', { enumerable: true, get() {
      if (vol === undefined) vol = haveR ? mcVolume(X, radii, opts.mcPerAtom || Math.max(MC_PER_ATOM, Math.ceil(MC_MIN_POINTS / N)), MC_SEED) : null;
      return vol ? vol.value : null; } });
    Object.defineProperty(out, 'vdwVolumeStdErr', { enumerable: true, get() { if (vol === undefined) void out.vdwVolume; return vol ? vol.stdErr : null; } });
    Object.defineProperty(out, 'vdwArea', { enumerable: true, get() {
      if (area === undefined) area = haveR ? vdwArea(X, radii, opts.areaPoints || AREA_POINTS) : null;
      return area; } });
    out.radiiNote = haveR ? null : 'no van der Waals radius for ' + record.atoms.filter((a, i) => !(radii[i] > 0)).map(a => a.el)[0];
    // dipole verdict from the conformer's symmetry: computed on first read, then cached
    let dip;
    Object.defineProperty(out, 'dipole', { enumerable: true, get() { if (dip === undefined) dip = dipoleInfo(record, X); return dip; } });
    return out;
  }

  /* ================================================================ 12b. dipole: what the conformer's symmetry allows
     A dipole moment μ is unchanged by every symmetry operation R of the rigid conformer: Rμ = μ. If no nonzero
     vector is left unchanged by all the operations found, μ = 0 exactly. A missed operation can only leave a vector
     unchanged, never remove one, so the search errs towards "a dipole is allowed", never towards a false zero.
     Operations are searched about the centroid: inversion, mirror planes, and rotations C2–C6 and S4 (C5 about
     principal axes). Candidate axes and plane normals: the principal axes, centroid→atom directions, and for every
     pair of atoms of the same element and isotope at the same distance from the centroid, the direction of their
     midpoint (a C2 axis that swaps them) and of their difference (the normal of a mirror that swaps them). An operation
     is accepted when every atom lands within DIPOLE.TOL of an atom of the same element and isotope; operations of one
     kind whose axes lie within 15° count once (real point groups up to n = 6 put them ≥ 20° apart), so a tilted copy
     of one axis cannot pass for a second axis. A CH₃-type end group (three identical terminal atoms on an sp³ atom
     bonded to a non-terminal atom) turns freely, and the dipoles of its three identical bonds cancel across its axis,
     so it enters as one point at the centre of its three atoms: the frozen orientation of a methyl group does not
     decide the answer (p-xylene, (E)-but-2-ene). */
  const SYM_MAX_ATOMS = 400, SYM_WORK = 4e6, SYM_SAME_AXIS = Math.cos(15 * Math.PI / 180), SYM_INVARIANT = 0.1;
  function eigenSym3(A) {                                // Jacobi: eigenvalues ascending with unit eigenvectors
    const a = [A[0].slice(), A[1].slice(), A[2].slice()], V = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    for (let sweep = 0; sweep < 50; sweep++) {
      if (a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2] < 1e-24) break;
      for (const [p, q] of [[0, 1], [0, 2], [1, 2]]) {
        if (Math.abs(a[p][q]) < 1e-300) continue;
        const th = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1)), c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < 3; k++) { const akp = a[k][p], akq = a[k][q]; a[k][p] = c * akp - s * akq; a[k][q] = s * akp + c * akq; }
        for (let k = 0; k < 3; k++) { const apk = a[p][k], aqk = a[q][k]; a[p][k] = c * apk - s * aqk; a[q][k] = s * apk + c * aqk; }
        for (let k = 0; k < 3; k++) { const vkp = V[k][p], vkq = V[k][q]; V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq; }
      }
    }
    return [0, 1, 2].map(k => ({ value: a[k][k], vector: [V[0][k], V[1][k], V[2][k]] })).sort((x, y) => x.value - y.value);
  }
  function symOp(u, ang, improper) {                    // rotation by ang about unit u, then (improper) the mirror ⊥ u
    const [x, y, z] = u, c = Math.cos(ang), s = Math.sin(ang), t = 1 - c;
    const Rm = [[t * x * x + c, t * x * y - s * z, t * x * z + s * y], [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
      [t * x * z - s * y, t * y * z + s * x, t * z * z + c]];
    if (!improper) return Rm;
    const Sm = [[1 - 2 * x * x, -2 * x * y, -2 * x * z], [-2 * x * y, 1 - 2 * y * y, -2 * y * z], [-2 * x * z, -2 * y * z, 1 - 2 * z * z]];
    return [0, 1, 2].map(i => [0, 1, 2].map(j => Sm[i][0] * Rm[0][j] + Sm[i][1] * Rm[1][j] + Sm[i][2] * Rm[2][j]));
  }
  // smallest eigenvalue of Σ (M − I)ᵀ(M − I) is ≈ 0 exactly when some direction is left unchanged by every operation
  function invariantLeft(ops) {
    const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const { M } of ops) for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++)
      for (let k = 0; k < 3; k++) A[i][j] += (M[k][i] - (k === i ? 1 : 0)) * (M[k][j] - (k === j ? 1 : 0));
    return jacobiEigen3(A)[0] < SYM_INVARIANT;
  }
  /* types: a string per point (element + isotope); X: Float64Array(3n). → {zero, complete, ops: [kinds found]}
     complete false: the work cap was reached before a zero was found (the answer is then unknown). */
  function symmetrySearch(types, X, tol) {
    const n = types.length, c = [0, 0, 0];
    for (let i = 0; i < n; i++) for (let d = 0; d < 3; d++) c[d] += X[3 * i + d] / n;
    const P = new Float64Array(3 * n), rad = new Float64Array(n);
    for (let i = 0; i < n; i++) { for (let d = 0; d < 3; d++) P[3 * i + d] = X[3 * i + d] - c[d]; rad[i] = Math.hypot(P[3 * i], P[3 * i + 1], P[3 * i + 2]); }
    const byType = new Map();
    for (let i = 0; i < n; i++) { if (!byType.has(types[i])) byType.set(types[i], []); byType.get(types[i]).push(i); }
    for (const [t, list] of byType) { list.sort((p, q) => rad[p] - rad[q]); byType.set(t, { idx: list, rs: list.map(i => rad[i]) }); }
    const tol2 = tol * tol;
    let work = 0;
    const find = (ty, qx, qy, qz) => {                  // an atom of type ty within tol of (qx, qy, qz)?
      const T = byType.get(ty), rq = Math.hypot(qx, qy, qz);
      let lo = 0, hi = T.rs.length;
      while (lo < hi) { const mid = (lo + hi) >> 1; if (T.rs[mid] < rq - tol) lo = mid + 1; else hi = mid; }
      for (let k = lo; k < T.rs.length && T.rs[k] <= rq + tol; k++) {
        work++;
        const j = T.idx[k], dx = P[3 * j] - qx, dy = P[3 * j + 1] - qy, dz = P[3 * j + 2] - qz;
        if (dx * dx + dy * dy + dz * dz <= tol2) return true;
      }
      return false;
    };
    const order = [...Array(n).keys()].sort((p, q) => rad[q] - rad[p]);       // outer atoms first: they fail fastest
    const test = M => {
      for (const i of order) {
        const x = P[3 * i], y = P[3 * i + 1], z = P[3 * i + 2];
        if (!find(types[i], M[0][0] * x + M[0][1] * y + M[0][2] * z, M[1][0] * x + M[1][1] * y + M[1][2] * z,
          M[2][0] * x + M[2][1] * y + M[2][2] * z)) return false;
      }
      return true;
    };
    const ops = [];
    const accept = (kind, u, M) => {                    // → true when no direction is left unchanged any more
      if (u && ops.some(o => o.kind === kind && Math.abs(o.u[0] * u[0] + o.u[1] * u[1] + o.u[2] * u[2]) > SYM_SAME_AXIS)) return false;
      ops.push({ kind, u, M });
      return !invariantLeft(ops);
    };
    const result = zero => ({ zero, complete: zero || work <= SYM_WORK, ops: ops.map(o => o.kind) });
    const inv = symOp([0, 0, 1], Math.PI, true);
    if (test(inv) && accept('i', null, inv)) return result(true);
    const seen = new Set();
    const tryAxis = (x, y, z, principal) => {           // → true when the dipole is shown to be zero
      const L = Math.hypot(x, y, z);
      if (L < 1e-6) return false;
      x /= L; y /= L; z /= L;
      if (x < -1e-9 || (Math.abs(x) <= 1e-9 && (y < -1e-9 || (Math.abs(y) <= 1e-9 && z < 0)))) { x = -x; y = -y; z = -z; }
      const key = Math.round(x * 60) + ',' + Math.round(y * 60) + ',' + Math.round(z * 60);
      if (seen.has(key)) return false;
      seen.add(key);
      const u = [x, y, z];
      const T = (kind, ang, imp) => { const M = symOp(u, ang, imp); return test(M) ? (accept(kind, u, M) ? 'zero' : 'found') : null; };
      if (T('σ', 0, true) === 'zero') return true;
      const c2 = T('C2', Math.PI, false); if (c2 === 'zero') return true;
      const c3 = T('C3', 2 * Math.PI / 3, false); if (c3 === 'zero') return true;
      if (c2 && (T('C4', Math.PI / 2, false) === 'zero' || T('S4', Math.PI / 2, true) === 'zero')) return true;
      if (c2 && c3 && T('C6', Math.PI / 3, false) === 'zero') return true;
      if (principal && T('C5', 2 * Math.PI / 5, false) === 'zero') return true;
      return false;
    };
    const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];         // principal axes of the point cloud (unit weights)
    for (let i = 0; i < n; i++) for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) C[a][b] += P[3 * i + a] * P[3 * i + b];
    for (const { vector: v } of eigenSym3(C)) if (tryAxis(v[0], v[1], v[2], true)) return result(true);
    for (const i of order) {
      if (work > SYM_WORK) return result(false);
      if (rad[i] > 0.3 && tryAxis(P[3 * i], P[3 * i + 1], P[3 * i + 2], false)) return result(true);
    }
    for (const { idx, rs } of byType.values()) {
      for (let p = 0; p < idx.length; p++) for (let q = p + 1; q < idx.length && rs[q] - rs[p] <= 2 * tol; q++) {
        if (work > SYM_WORK) return result(false);
        const i = idx[p], j = idx[q];
        if (tryAxis(P[3 * i] + P[3 * j], P[3 * i + 1] + P[3 * j + 1], P[3 * i + 2] + P[3 * j + 2], false)) return result(true);
        if (tryAxis(P[3 * i] - P[3 * j], P[3 * i + 1] - P[3 * j + 1], P[3 * i + 2] - P[3 * j + 2], false)) return result(true);
      }
    }
    return result(false);
  }
  /* Topology for the dipole verdict (cached per record):
     tops — CH₃-type groups: an sp³ atom with three identical terminal atoms and one non-terminal neighbour;
     torsions — "dipole-changing" bonds: non-ring single bonds whose two ends each have another neighbour, where
       neither end is linear (sp) nor a symmetric top (three identical substituents: CH₃, CF₃, CCl₃, NH₃⁺, and
       tert-butyl-like groups made of such tops). Turning such a bond can change the molecule's symmetry. */
  const dipTopoCache = new WeakMap();
  function dipoleTopology(record) {
    if (dipTopoCache.has(record)) return dipTopoCache.get(record);
    const A = record.atoms, key = i => A[i].el + ':' + (A[i].isotope || 0);
    const single4 = u => A[u].neighbours.length === 4 && A[u].neighbours.every(x => x.order === 1 && !x.aromatic);
    const sig = (x, from, depth) => {                    // signature of the branch at x seen from `from`, or null
      if (A[x].neighbours.length === 1) return key(x);
      if (depth <= 0 || !single4(x)) return null;
      const s = A[x].neighbours.filter(y => y.atom !== from).map(y => sig(y.atom, x, depth - 1));
      return s.length === 3 && s[0] !== null && s.every(t => t === s[0]) ? key(x) + '(' + s[0] + ')' : null;
    };
    const topEnd = (u, v) => A[u].neighbours.length > 1 && sig(u, v, 2) !== null;
    const linear = u => A[u].neighbours.some(x => x.order === 3 && !x.aromatic) ||
      A[u].neighbours.filter(x => x.order === 2 && !x.aromatic).length >= 2 || A[u].hybridization.label === 'sp';
    const tops = [];
    A.forEach((a, u) => {
      if (!single4(u)) return;
      const term = a.neighbours.filter(x => A[x.atom].neighbours.length === 1).map(x => x.atom);
      const rest = a.neighbours.filter(x => A[x.atom].neighbours.length > 1).map(x => x.atom);
      if (rest.length !== 1 || term.length !== 3 || !term.every(t => key(t) === key(term[0]))) return;
      tops.push({ centre: u, members: term, partner: rest[0], type: 'top:' + key(term[0]) });
    });
    const torsions = [];
    record.bonds.forEach((b, k) => {
      if (b.order !== 1 || b.aromatic || b.inRing || b.isHBond) return;
      if (A[b.a].neighbours.length < 2 || A[b.b].neighbours.length < 2) return;
      if (linear(b.a) || linear(b.b) || topEnd(b.a, b.b) || topEnd(b.b, b.a)) return;
      torsions.push(k);
    });
    // a flexible ring also changes shape (cyclohexane-1,4-dione: its C2h chair has no dipole, but less symmetric ring
    // conformations can have one):
    // four or more members, not aromatic, three or more sp³ atoms, and a ring bond between two atoms that are in no
    // other ring. The last condition leaves out rigid cages (adamantane, cubane), where every such bond has an end
    // shared with another ring; it keeps fused and bridged rings that can flip (decalin, norbornane's C2–C3).
    const ringCount = new Map();
    record.rings.forEach(r => r.forEach(i => ringCount.set(i, (ringCount.get(i) || 0) + 1)));
    const flexibleRing = record.rings.some(r => r.length >= 4 && !r.every(i => A[i].aromatic) &&
      r.filter(i => A[i].hybridization.label === 'sp3').length >= 3 &&
      r.some((i, t) => ringCount.get(i) === 1 && ringCount.get(r[(t + 1) % r.length]) === 1));
    // five or six electron domains (SF4, ClF3, PCl5, SF6): MMFF94 has no types for them, and a geometry from elsewhere
    // cannot be checked by the force field
    const hypervalent = A.some(a => !a.isH && a.stericNumber >= 5);
    const topo = { tops, torsions, flexibleRing, hypervalent };
    dipTopoCache.set(record, topo);
    return topo;
  }
  function symmetryOf(record, topo, X) {                 // the points the search sees: each top folded into one point
    const removed = new Set(), types = [], pts = [];
    for (const t of topo.tops) t.members.forEach(i => removed.add(i));
    record.atoms.forEach((a, i) => { if (!removed.has(i)) { types.push(a.el + ':' + (a.isotope || 0)); pts.push(X[3 * i], X[3 * i + 1], X[3 * i + 2]); } });
    for (const t of topo.tops) {
      const p = [0, 1, 2].map(d => t.members.reduce((s, i) => s + X[3 * i + d], 0) / 3);
      types.push(t.type); pts.push(p[0], p[1], p[2]);
    }
    return symmetrySearch(types, Float64Array.from(pts), DIPOLE.TOL);
  }
  // Does the dipole stay zero at any angle of torsion bond k? Tested by turning one side to three unrelated angles
  // (biphenyl stays D2 at every twist; hydrogen peroxide and butane become polar).
  function torsionKeepsZero(record, topo, X, k) {
    const b = record.bonds[k], side = new Set([b.b]), st = [b.b];
    while (st.length) {
      const u = st.pop();
      for (const x of record.atoms[u].neighbours) if (x.atom !== b.a && !side.has(x.atom)) { side.add(x.atom); st.push(x.atom); }
    }
    if (side.has(b.a)) return false;
    const ax = [0, 1, 2].map(d => X[3 * b.b + d] - X[3 * b.a + d]), L = Math.hypot(ax[0], ax[1], ax[2]), u = ax.map(v => v / L);
    for (const deg of [37, 83, 131]) {
      const M = symOp(u, deg * Math.PI / 180, false), Y = Float64Array.from(X);
      for (const i of side) {
        const v = [0, 1, 2].map(d => X[3 * i + d] - X[3 * b.b + d]);
        for (let d = 0; d < 3; d++) Y[3 * i + d] = X[3 * b.b + d] + M[d][0] * v[0] + M[d][1] * v[1] + M[d][2] * v[2];
      }
      if (!symmetryOf(record, topo, Y).zero) return false;
    }
    return true;
  }
  function dipoleInfo(record, X) {
    const V = TEXT.dipoleVerdict;
    if (record.charge !== 0) return { verdict: 'ion', text: V.ion, note: null, symmetry: null };
    if (record.N > SYM_MAX_ATOMS) return { verdict: 'not assessed', text: V.notAssessed, note: TEXT.dipole, symmetry: null };
    const topo = dipoleTopology(record);
    if (topo.hypervalent) return { verdict: 'not assessed', text: V.hypervalent, note: TEXT.dipole, symmetry: null };
    const s = symmetryOf(record, topo, X);
    // a polar conformer of a molecule that can change shape (azobenzene twisted by MMFF94, cyclopentane's envelope)
    // says so: the verdict is about this conformation
    if (!s.zero) return s.complete ? { verdict: 'polar', text: topo.torsions.length || topo.flexibleRing ? V.polarConformer : V.polar,
      note: TEXT.dipole, symmetry: s.ops } : { verdict: 'not assessed', text: V.notAssessed, note: TEXT.dipole, symmetry: null };
    let open = topo.torsions;
    if (open.length === 1 && torsionKeepsZero(record, topo, X, open[0])) open = [];
    if (open.length) return { verdict: 'zero-conformer', text: V.zeroConformer, note: TEXT.dipole, symmetry: s.ops, torsions: open.slice() };
    // a ring that can change conformation: the symmetry found holds for this conformer only
    if (topo.flexibleRing) return { verdict: 'zero-conformer', text: V.zeroRing, note: TEXT.dipole, symmetry: s.ops, ring: true };
    return { verdict: 'zero', text: V.zero, note: TEXT.dipole, symmetry: s.ops };
  }
  /* The Gasteiger-charge dipole (debye) the DIPOLE figures were measured with. Internal (tests): never shown. */
  function gasteigerDipole(record, X, charges) {
    const q = charges || record.atoms.map(a => a.gasteiger);
    if (q.some(v => v == null || !Number.isFinite(v))) return null;
    const p = [0, 0, 0];
    for (let i = 0; i < record.N; i++) for (let d = 0; d < 3; d++) p[d] += q[i] * X[3 * i + d];
    return Math.hypot(p[0], p[1], p[2]) * 4.80320;        // e·Å → debye
  }

  /* ================================================================ 13. rows (hover card, pinned details, tables) */
  const row = (key, label, value, kind, note, main) => {
    const r = { key, label, value: value === null || value === undefined ? '—' : String(value), kind };
    if (note) r.note = note;
    if (main) r.main = true;
    return r;
  };
  const listNames = (record, idx) => idx.map(i => atomName(record, i)).join(', ');
  // SSSR rings through an atom, e.g. "in a 6-membered ring", "in two rings (5- and 6-membered)", "(both 6-membered)"
  const ringText = list => {
    if (!list.length) return 'not in a ring';
    if (list.length === 1) return `in a ${list[0]}-membered ring`;
    const same = list.every(x => x === list[0]);
    return `in ${countWord(list.length)} rings (` + (same ? `${list.length === 2 ? 'both' : 'all'} ${list[0]}-membered`
      : list.slice(0, -1).join('-, ') + '- and ' + list[list.length - 1] + '-membered') + ')';
  };

  function atomRows(record, i, geo) {
    const a = record.atoms[i], rows = [];
    if (!a) return rows;
    const e = EL(a.z);
    rows.push(row('atom', 'Atom', `${a.el}${i + 1} · ${a.name.toLowerCase()}` + (a.appended ? ` (on ${atomName(record, a.parent)})` : ''), 'exact', null, true));
    rows.push(row('z', 'Atomic number', a.z, 'exact'));
    if (!a.appended) rows.push(row('charge', 'Formal charge', intSigned(a.charge), 'convention', 'bookkeeping ' + TEXT.asDrawn +
      ', not the real charge distribution (see the partial charge)', a.charge !== 0));
    if (a.isotope) rows.push(row('isotope', 'Isotope', supDigits(String(a.isotope)) + a.el + (a.el === 'H' && a.isotope === 2 ? ' (deuterium)' : a.el === 'H' && a.isotope === 3 ? ' (tritium)' : ''), 'exact', null, true));
    if (a.radical) rows.push(row('radical', 'Unpaired electrons', plural(a.radical, 'unpaired electron') + ' (radical)', 'exact', null, true));
    const nbText = a.neighbours.length ? `${plural(a.neighbours.length, 'atom')}: ${listNames(record, a.neighbours.map(x => x.atom))}` +
      (a.hCount && a.neighbours.length > 1 ? ` (${a.hCount} of them H)` : '') : 'nothing (a single atom or ion)';
    rows.push(row('neighbours', 'Bonded to', nbText, 'exact', null, true));
    if (!a.isH) rows.push(row('hcount', 'Hydrogens attached', a.hCount, 'exact'));
    if (!a.isH) {
      rows.push(row('aromatic', 'Aromatic', a.aromatic ? 'yes' : 'no', 'convention', TEXT.aromaticModel, a.aromatic));
      rows.push(row('ring', 'Ring', ringText(a.ringList), 'exact', TEXT.sssr));
    }
    if (a.cip) rows.push(row('cip', 'Stereocentre', a.cip === '?' ? a.cipText :
      a.cip === 'r' || a.cip === 's' ? `pseudoasymmetric (${a.cip})` : `(${a.cip})`, a.cip === '?' ? 'exact' : 'convention',
      a.cip === 'r' || a.cip === 's' ? TEXT.cip + '; lowercase r/s marks a pseudoasymmetric centre, whose two ligands differ only in configuration' : TEXT.cip, true));
    if (a.valenceElectrons !== null) rows.push(row('valence', 'Valence electrons', a.valenceElectrons, 'exact'));
    if (a.valenceElectrons !== null && a.isH) rows.push(row('octet', 'Electrons around the atom', `${a.electronsAround} (${a.octetText})`, 'convention'));
    else if (a.valenceElectrons !== null) {
      if (a.lonePairs === null) rows.push(row('lonepairs', 'Lone pairs', 'unusual electron count', 'convention', TEXT.asDrawn, true));
      else {
        rows.push(row('lonepairs', 'Lone pairs', a.lonePairs, 'convention', 'lone pairs ' + TEXT.asDrawn, true));
        rows.push(row('octet', 'Electrons around the atom', a.octetText.includes('(') ? `${a.electronsAround}, ${a.octetText}` : `${a.electronsAround} (${a.octetText})`, 'convention',
          'bonding pairs counted fully plus non-bonding electrons, ' + TEXT.asDrawn));
      }
    }
    if (a.stericNumber !== null) {
      const sigma = a.degree, lp = a.lonePairs;
      let v = `${a.stericNumber} (${plural(sigma, 'σ bond')} + ${plural(lp, 'lone pair')})`;
      if (a.effectiveStericNumber !== a.stericNumber) v += `; effectively ${a.effectiveStericNumber}: one lone pair is delocalised`;
      rows.push(row('sn', 'Steric number', v, 'model', 'VSEPR model'));
    }
    if (a.domainGeometry || a.molecularGeometry) {
      const mg = a.molecularGeometry;
      rows.push(row('geometry', 'Shape around the atom (VSEPR)', mg && mg.startsWith('terminal') ? `${a.domainGeometry} electron domains; terminal atom, no bond angle` :
        `${a.domainGeometry || '—'} electron domains` + (mg ? `; ${mg}` : ''), 'model', 'VSEPR model', true));
    }
    if (a.idealAngleText || a.idealAngleNote || a.ringAngleNote) {
      const notes = [a.idealAngleNote, a.ringAngleNote].filter(Boolean).join('; ');
      rows.push(row('ideal', 'Ideal bond angle (VSEPR)', a.idealAngleText || 'no single ideal value', 'model', notes || 'VSEPR model'));
    }
    const hy = a.hybridization;
    if (a.isH) rows.push(row('hyb', 'Hybridization', 'not hybridized (1s orbital)', 'convention', null, true));
    else rows.push(row('hyb', 'Hybridization', hy.label ? hybText(hy.label) + (hy.qualifier ? ' (conventionally)' : '') : 'not assigned', 'model',
      hy.label ? hy.note : hy.note, true));
    const os = a.oxidationState;
    rows.push(row('os', 'Oxidation state', os.text, 'convention',
      [os.value !== null ? "bonding electrons go to the more electronegative atom (Pauling scale here; IUPAC's 2016 definition " +
        'uses Allen electronegativities, which differ only where flagged)' : os.note,
        os.value !== null ? os.note : null, os.flag].filter(Boolean)
        .map((s, i) => (i ? s.charAt(0).toUpperCase() + s.slice(1) : s)).join('. '), true));   // each joined part starts a sentence
    if (a.gasteiger !== null) {
      const notes = [TEXT.gasteigerNote];
      if (record.hasHypervalentSP) notes.push(TEXT.gasteigerSO);
      rows.push(row('gasteiger', 'Partial charge', signed(a.gasteiger, 3), 'model', TEXT.gasteiger + '. ' + notes.join(' '), true));
    } else rows.push(row('gasteiger', 'Partial charge', record.gasteigerNote || 'not available', 'model', TEXT.gasteiger, true));
    if (a.en !== null || a.enAllen !== null) {
      const v = (a.en !== null ? a.en.toFixed(2) + ' (Pauling)' : 'no Pauling value') + (a.enAllen !== null ? `; Allen ${a.enAllen.toFixed(3)}` : '');
      rows.push(row('en', 'Electronegativity', v, 'table', 'Pauling scale: Allred 1961 values' +
        (a.enAllen !== null ? '; Allen scale (average valence-electron energy): Allen 1989 / Mann, Meek & Allen 2000' : '') + (e.enNote ? '. ' + e.enNote : '')));
    }
    if (a.isotope) {
      const lab = supDigits(String(a.isotope)) + a.el;
      if (a.isotopeMass) rows.push(row('mass', 'Atomic mass', `${a.isotopeMass.toFixed(4)} u (${lab})`, 'table', "mass of this isotope: RDKit's isotope table (NIST atomic masses)"));
      else rows.push(row('mass', 'Mass number', `${a.isotope} (${lab})`, 'exact', 'the isotope label; no tabulated mass for this isotope'));
    }
    else if (e.aw != null) rows.push(row('mass', 'Atomic mass', `${e.aw} (standard atomic weight)`, 'table', 'CIAAW 2024' + (e.awKind === 'conventional' ? ' (conventional value)' : '')));
    else if (e.massNumber) rows.push(row('mass', 'Atomic mass', `[${e.massNumber.join(' or ')}] (mass number of the longest-lived isotope)`, 'table', 'no standard atomic weight: radioactive element (NUBASE2020)'));
    if (a.rcov !== null) {
      let v = `${a.rcov.toFixed(2)} Å`;
      if (e.rcovHyb) v += ` (sp³ ${e.rcovHyb.sp3.toFixed(2)}, sp² ${e.rcovHyb.sp2.toFixed(2)}, sp ${e.rcovHyb.sp.toFixed(2)})`;
      rows.push(row('rcov', 'Covalent radius', v, 'table', 'Cordero et al. 2008'));
    }
    if (a.rvdw !== null) rows.push(row('rvdw', 'van der Waals radius', `${a.rvdw.toFixed(2)} Å`, 'table',
      ({ Bondi1964: 'Bondi 1964', RowlandTaylor1996: 'Rowland & Taylor 1996', Mantina2009: 'Mantina et al. 2009', Alvarez2013: 'Alvarez 2013' }[e.src && e.src.rvdw] || '') +
      (e.vdwNote ? '. ' + e.vdwNote : '')));
    if (a.groups.length) rows.push(row('groups', 'Functional group', a.groups.join('; '), 'convention', TEXT.functionalGroups, true));
    if (a.neighbours.length >= 2) {
      if (!geo) { /* no 3D model shown */ }
      else if (!geo.optimised) rows.push(row('angles', 'Bond angles (this conformer)', TEXT.notOptimised, '3d', TEXT.notOptimisedBanner, true));
      else {
        const fmt = x => `${atomName(record, x.i)}–${a.el}${i + 1}–${atomName(record, x.j)} ${x.deg.toFixed(1)}°`;
        const L = geo.angles[i], bad = L.filter(x => x.unreliable), good = L.filter(x => !x.unreliable);
        // angles that involve an atom of a flagged bond say so (all of them: once, at the end; some: listed apart)
        const v = !bad.length ? L.map(fmt).join(', ') : !good.length ? L.map(fmt).join(', ') + ' ' + TEXT.unreliableHere
          : good.map(fmt).join(', ') + '; ' + TEXT.unreliableHere.slice(1, -1) + ': ' + bad.map(fmt).join(', ');
        rows.push(row('angles', 'Bond angles (this conformer)', v, '3d', bad.length ? TEXT.geometry3d + '; ' + TEXT.unreliableNote : TEXT.geometry3d, true));
      }
    }
    return rows;
  }

  function bondRows(record, k, geo) {
    const b = record.bonds[k], rows = [];
    if (!b) return rows;
    const A = record.atoms[b.a], B = record.atoms[b.b], nm = `${atomName(record, b.a)}–${atomName(record, b.b)}`;
    rows.push(row('bond', 'Bond', nm, 'exact', null, true));
    if (b.aromatic) {
      // the 3D view draws every aromatic bond alike, so the wording names the 2D drawing explicitly
      rows.push(row('type', 'Bond type', `aromatic (drawn ${b.kekule} in the 2D Kekulé structure; π electrons delocalised)`, 'convention', TEXT.aromaticModel, true));
      const kd = b.kekuleDouble;
      rows.push(row('order', 'Bond order', b.bondOrder === 1.5 ? '1.5' : 'between single and double', 'convention',
        kd ? `double in ${kd.times} of the ${kd.of} Kekulé structures` : null));
    } else {
      rows.push(row('type', 'Bond type', b.typeText, 'exact', null, true));
      rows.push(row('order', 'Bond order', b.orderText, 'exact'));
    }
    rows.push(row('sigmapi', 'σ and π bonds', b.sigmaPiText, 'convention', null, true));
    if (b.overlap.length || b.piNote) rows.push(row('overlap', 'Orbital overlap', b.overlap.join('; ') || '—', 'model', b.piNote ||
      'valence-bond picture using the hybridization shown for each atom'));
    if (b.polarity) rows.push(row('polarity', 'Polarity', b.polarity.text, 'convention', TEXT.polarityCutoffs, true));
    if (b.aromatic) rows.push(row('enthalpy', 'Bond enthalpy (average)', 'no single average (resonance)', 'table', MD().NOTES.aromaticEnthalpy));
    else if (b.resonance) rows.push(row('enthalpy', 'Bond enthalpy (average)', 'no single average (resonance-delocalised)', 'table',
      'this bond and its partners in the group (nitro, carboxylate, sulfonate, azide, diazo, amidinium …) are equivalent by ' +
      'resonance, so neither the single- nor the double-bond average applies'));
    else if (b.enthalpy) rows.push(row('enthalpy', 'Bond enthalpy (average)', `${b.enthalpy.kJ} kJ/mol`, 'table',
      [TEXT.enthalpy, b.enthalpy.source, b.enthalpy.note].filter(Boolean).join('. ')));
    if (b.typical) {
      const t = b.typical;
      const v = `${t.A.toFixed(3)} Å` + (t.estimate ? ' (rough estimate)' : t.sd != null && t.n != null ? ` (s.d. ${t.sd}, n = ${t.n})` : ' (no spread given)');
      rows.push(row('typical', 'Typical length', v, 'table', t.estimate ? t.note : [TEXT.length, t.sub ? 'entry: ' + t.sub : null, t.note].filter(Boolean).join('. ')));
    }
    if (geo) {
      if (geo.flagged && geo.flagged.has(k) && Number.isFinite(geo.bondLength[k]))
        rows.push(row('length', 'Length (this conformer)', `${geo.bondLength[k].toFixed(3)} Å (flagged: not reliable)`, '3d', geo.flagged.get(k), true));
      else if (!geo.optimised) rows.push(row('length', 'Length (this conformer)', TEXT.notOptimised, '3d', TEXT.notOptimisedBanner, true));
      else rows.push(row('length', 'Length (this conformer)', `${geo.bondLength[k].toFixed(3)} Å`, '3d', TEXT.geometry3d, true));
      if (b.torsionAtoms) {
        if (!geo.optimised) rows.push(row('dihedral', 'Dihedral (this conformer)', TEXT.notOptimised, '3d', TEXT.notOptimisedBanner));
        else {
          const d = geo.dihedral[k];
          rows.push(row('dihedral', 'Dihedral (this conformer)', `${num(d.deg, 1)}° ${d.atoms.map(i => atomName(record, i)).join('–')}: ${d.name}` + (d.common ? `; ${d.common}` : '') +
            (d.unreliable ? ' ' + TEXT.unreliableHere : ''), '3d', 'IUPAC sign convention; names after Klyne & Prelog; ' + TEXT.dihedral + (d.unreliable ? '; ' + TEXT.unreliableNote : ''), true));
        }
      }
    }
    if (b.rotation) rows.push(row('rotation', 'Rotation about this bond', b.rotation.text, 'exact', TEXT.rotation, true));
    rows.push(row('rotatable', 'Counted as rotatable (RDKit strict)', b.rotatable.yes ? 'yes' : 'no: ' + b.rotatable.reason, 'convention', TEXT.rotatable));
    if (b.amideNote) rows.push(row('amide', 'Amide C–N', b.amideNote, 'convention', null, true));
    rows.push(row('conjugated', 'Conjugated', b.conjugated ? 'yes' : 'no', 'convention', TEXT.conjugated));
    if (!b.isHBond) rows.push(row('ring', 'In ring', b.inRing ? 'yes (' + b.ringSizes.join('- and ') + '-membered)' : 'no', 'exact', TEXT.sssr));
    if (b.ez) rows.push(row('ez', 'E/Z', b.ez === 'unspecified' ? 'E/Z possible: not specified in the input' : `(${b.ez})`, b.ez === 'unspecified' ? 'exact' : 'convention', TEXT.cip, true));
    if (b.groups.length) rows.push(row('groups', 'Functional group', b.groups.join('; '), 'convention', TEXT.functionalGroups));
    void A; void B;
    return rows;
  }

  function moleculeSections(record, geo, extra) {
    extra = extra || {};
    const S = [], r = record;
    // identity
    const id = [];
    if (extra.component) id.push(row('component', 'Component', `${extra.component.index + 1} of ${extra.component.count}` +
      (extra.component.formula ? ` (${extra.component.formula})` : ''), 'exact', 'a salt or mixture: each component is analysed on its own', true));
    const nm = extra.names;
    if (nm) {
      if (nm.common) id.push(row('common', 'Common name', nm.common, 'exact', "from the app's compound library (identity checked by SMILES)", true));
      id.push(row('iupac', 'Systematic name', nm.iupac || '—', 'convention', nm.iupac ? TEXT.namer : nm.iupacNote, true));
      if (nm.cas) id.push(row('cas', 'CAS number', nm.cas, 'exact', "from the app's compound library (identity checked by SMILES)"));
    }
    id.push(row('formula', 'Molecular formula', r.formulaText || 'not defined (unspecified atom *)', 'exact', 'Hill order', true));
    if (r.isotopicFormula) id.push(row('isoformula', 'Isotopic formula', r.isotopicFormula, 'exact'));
    id.push(row('smiles', 'SMILES (RDKit canonical)', r.smiles.isomeric, 'exact', 'canonical within RDKit only; other programs write other SMILES'));
    if (r.smiles.noStereo && r.smiles.noStereo !== r.smiles.isomeric) id.push(row('smiles0', 'SMILES without stereo', r.smiles.noStereo, 'exact',
      'the same structure with its stereo marks removed' + (r.isotopicFormula ? '; isotope labels kept' : '')));
    if (r.inchi) id.push(row('inchi', 'InChI', r.inchi, 'exact', r.chargedAtoms ? 'standard InChI normalises protons and charges (a zwitterion gets the same InChI as the neutral form)' : 'standard InChI'));
    if (r.inchikey) id.push(row('inchikey', 'InChIKey', r.inchikey, 'exact'));
    if (extra.combined && extra.combined.components.length > 1) id.push(row('combined', 'Whole input', `${extra.combined.formulaText}, ${extra.combined.mw.rounded.toFixed(2)} g/mol`, 'exact',
      extra.combined.mw.method === 'CIAAW' ? 'all components together; CIAAW 2024 atomic weights' : extra.combined.mw.note));
    S.push({ title: 'Identity', rows: id });

    // mass and composition
    const ms = [];
    ms.push(row('mw', 'Molecular weight', `${r.mw.rounded.toFixed(2)} g/mol`, r.mw.method === 'CIAAW' ? 'exact' : 'table',
      r.mw.method === 'CIAAW' ? 'sum of CIAAW 2024 standard atomic weights' : r.mw.note, true));
    // RDKit's exact mass already subtracts the electron mass (0.000549 Da) per positive charge and adds it per negative
    // charge: the value is the mass of the ion itself (checked against Descriptors.ExactMolWt)
    ms.push(row('mono', r.isotopicFormula ? 'Exact mass (this isotopologue)' : 'Monoisotopic mass', `${r.monoisotopic.toFixed(4)} Da`, 'exact',
      (r.isotopicFormula ? "labelled atoms use their isotope's mass, the others the most abundant isotope of their element (RDKit)"
        : 'most abundant isotope of each element (RDKit)') +
      (r.charge ? '; the mass of the ion: includes the electrons removed or added (0.000549 Da each)' : '')));
    if (r.composition) ms.push(row('composition', 'Composition (by mass)', r.composition.map(c => `${c.el} ${c.pct.toFixed(2)}%`).join(', '), 'exact', 'isotopes merged'));
    ms.push(row('atoms', 'Atoms', `${r.heavyAtoms} heavy, ${r.totalAtoms} in total (with H)`, 'exact', 'D and T count as hydrogen'));
    ms.push(row('charge', 'Net charge', intSigned(r.charge) + (r.chargedAtoms ? ` (${plural(r.chargedAtoms, 'charged atom')})` : ''), 'exact', null, r.charge !== 0 || r.chargedAtoms > 0));
    if (r.radicals) ms.push(row('radicals', 'Unpaired electrons', r.radicals, 'exact'));
    S.push({ title: 'Mass and composition', rows: ms });

    // structure
    const st = [];
    st.push(row('dou', 'Degree of unsaturation', r.dou.value !== null ? trimNum(r.dou.value, 1) : '—', 'exact',
      r.dou.value !== null ? 'rings + π bonds from the formula: C + Si + 1 + (N + P)/2 − (H + X)/2' + (r.dou.note ? '; ' + r.dou.note : '') : r.dou.note, true));
    const rc = r.ringCounts;
    st.push(row('rings', 'Rings', rc.total ? `${rc.total} (${rc.saturated} saturated; ${rc.heterocycles} heterocyclic)` : '0', 'exact', TEXT.sssr, true));
    if (rc.total) st.push(row('aromrings', 'Aromatic rings', `${rc.aromatic} aromatic, ${rc.aliphatic} not aromatic`, 'convention', TEXT.aromaticModel));
    if (rc.spiro || rc.bridgehead) st.push(row('spiro', 'Spiro / bridgehead atoms', `${rc.spiro} / ${rc.bridgehead}`, 'convention',
      'spiro atom: the only atom shared by two ring systems (IUPAC)' + (rc.spiroRDKit !== rc.spiro ? `; RDKit's SSSR-based NumSpiroAtoms gives ${rc.spiroRDKit}` : '') +
      "; bridgehead atoms: RDKit's NumBridgeheadAtoms"));
    for (const s of r.aromaticSystems) st.push(row('aromsys', 'Aromatic π system', `${listNames(r, s.atoms)}: ${s.piElectrons} π electrons` +
      (s.huckelN !== null ? ` (4n+2, n = ${s.huckelN})` : ''), 'convention', 'counting 1 per atom with an in-ring double bond, 2 per lone-pair atom (pyrrole N, furan O), 0 for C⁺ or an exocyclic C=O'));
    const rb = r.rotatableBonds;
    st.push(row('rotatable', 'Rotatable bonds', rb.count + (rb.count !== rb.rdkit ? ` (RDKit's descriptor gives ${rb.rdkit})` : ''), 'convention', TEXT.rotatable, true));
    if (r.fsp3 !== null) st.push(row('fsp3', 'Fraction sp³ carbon', r.fsp3.toFixed(2), 'convention', 'RDKit FractionCSP3'));
    S.push({ title: 'Structure', rows: st });

    // stereochemistry
    const sc = [];
    sc.push(row('chirality', 'Chirality', r.chirality.verdict + ': ' + r.chirality.reason, 'exact',
      r.chirality.verdict === 'meso' ? 'meso (IUPAC): an achiral member of a set of stereoisomers that also has chiral members' : null, true));
    const atomsS = r.stereo.atoms.filter(s => s.label !== '?'), unspec = r.stereo.atoms.filter(s => s.label === '?');
    const RING_WORD = { substituent: 'cis/trans', fusion: 'cis/trans ring fusion', bridgehead: 'in/out bridgehead' };
    sc.push(row('centres', 'Stereocentres', r.stereo.atoms.length ? r.stereo.atoms.map(s => `${atomName(r, s.index)} (${s.label === '?' ?
      'not specified' + (s.kind === 'ring' ? ', ' + (RING_WORD[s.ring] || 'cis/trans') : '') : s.label})`).join(', ') : 'none',
      'convention', TEXT.cip + (atomsS.some(s => s.label === 'r' || s.label === 's') ? '; lowercase r/s = pseudoasymmetric' : '') +
      (r.stereo.atoms.some(s => s.kind === 'ring') ? '; cis/trans: which face of the ring a substituent or a ring junction is on (RDKit lists these only once they are specified)' : '') +
      (r.stereo.ignored && r.stereo.ignored.length ? '. Not counted: ' + r.stereo.ignored.map(x => `${atomName(r, x.index)} (${x.reason})`).join('; ') : '')));
    const ez = r.stereo.bonds.filter(s => s.label !== 'unspecified'), ezU = r.stereo.bonds.filter(s => s.label === 'unspecified');
    if (ez.length) sc.push(row('ez', 'Double-bond geometry', ez.map(s => `${atomName(r, r.bonds[s.index].a)}=${atomName(r, r.bonds[s.index].b)} (${s.label})`).join(', '), 'convention', TEXT.cip));
    const unC = unspec.filter(s => s.kind !== 'ring'), unRing = kind => unspec.filter(s => s.kind === 'ring' && s.ring === kind);
    if (unspec.length || ezU.length) sc.push(row('unspecified', 'Unspecified stereo', [unC.length ? plural(unC.length, 'stereocentre') + ': ' + unC.map(s => atomName(r, s.index)).join(', ') : null,
      ...['substituent', 'fusion', 'bridgehead'].map(k => unRing(k).length ? `${{ substituent: 'cis/trans ring stereo', fusion: 'cis/trans ring fusion', bridgehead: 'in/out bridgehead' }[k]} at ` +
        unRing(k).map(s => atomName(r, s.index)).join(', ') : null),
      ezU.length ? plural(ezU.length, 'double bond') + ': ' + ezU.map(s => `${atomName(r, r.bonds[s.index].a)}=${atomName(r, r.bonds[s.index].b)}`).join(', ') : null].filter(Boolean).join('; '),
      'exact', 'the input does not say which stereoisomer; the 3D model picks one', true));
    S.push({ title: 'Stereochemistry', rows: sc });

    // functional groups
    S.push({ title: 'Functional groups', rows: r.groups.length ? r.groups.map((G, t) => row('group' + t, G.name, listNames(r, G.atoms), 'convention', TEXT.functionalGroups))
      : [row('group', 'Functional groups', 'none recognised', 'convention', TEXT.functionalGroups)] });

    // estimates
    const pr = [];
    pr.push(row('hbd', 'H-bond donors (Lipinski)', r.hbd, 'convention', TEXT.lipinski));
    pr.push(row('hba', 'H-bond acceptors (Lipinski)', r.hba, 'convention', TEXT.lipinski));
    pr.push(row('tpsa', 'Polar surface area (TPSA)', `${r.tpsa.toFixed(1)} Å²`, 'model', TEXT.tpsa));
    pr.push(row('logp', 'logP (octanol/water)', num(r.logp, 1), 'model', TEXT.logp));
    pr.push(row('mr', 'Molar refractivity', `${r.mr.toFixed(1)} cm³/mol`, 'model', TEXT.mr));
    if (r.ro5) pr.push(row('ro5', 'Rule of five', r.ro5.pass ? 'pass' + (r.ro5.violations.length ? ' (1 violation: ' + r.ro5.violations[0] + ')' : '') : 'fail: ' + r.ro5.violations.join(', '), 'convention', TEXT.ro5));
    if (r.veber) pr.push(row('veber', 'Veber', r.veber.pass ? 'pass' : 'fail', 'convention', TEXT.veber));
    S.push({ title: 'Properties (estimates and rules)', rows: pr });

    // 3D
    if (geo || extra.conformer) {
      const g3 = [];
      const cf = extra.conformer;
      if (cf) {
        if (cf.energy != null && Number.isFinite(cf.energy)) g3.push(row('energy', 'MMFF94 energy', `${num(cf.energy, 1)} kcal/mol` +
          (cf.relEnergy != null && Number.isFinite(cf.relEnergy) ? ` (ΔE ${num(cf.relEnergy, 1)} vs the lowest found)` : ''), '3d', TEXT.mmffEnergy, true));
        if (cf.source) g3.push(row('source', '3D model from', cf.source, '3d', null));
      }
      if (geo && !geo.optimised) g3.push(row('geometry', '3D values', TEXT.notOptimised, '3d', TEXT.notOptimisedBanner, true));
      else if (geo && geo.molecule && geo.molecule.withheld) {
        // a flagged heavy-atom bond: one row naming it instead of the size, shape and dipole values
        g3.push(row('geometry', '3D values', geo.molecule.withheld.text, '3d', geo.molecule.withheld.note, true));
      } else if (geo && geo.molecule) {
        const m = geo.molecule;
        const d = m.dipole;
        // symmetry decides whether a dipole can exist: an exact consequence of the conformer's geometry ('3d')
        if (d) g3.push(row('polarity', 'Dipole moment', d.text, '3d', d.note || null, true));
        if (m.rg !== null) g3.push(row('rg', 'Radius of gyration', `${m.rg.toFixed(2)} Å`, '3d', 'mass-weighted, all atoms; ' + TEXT.geometry3d));
        if (m.pmi) g3.push(row('pmi', 'Principal moments of inertia', m.pmi.map(v => v.toFixed(1)).join(', ') + ' amu·Å²', '3d', 'about the centre of mass; ' + TEXT.geometry3d));
        // linear / planar first: every planar molecule lies on the rod–disc edge of the NPR triangle
        const flat2 = m.form === 'linear' || m.form === 'planar';
        if (m.shape || flat2) g3.push(row('shape', 'Shape', (flat2 ? m.form : `${m.shape}-like`) +
          (m.npr1 !== null ? ` (NPR1 ${m.npr1.toFixed(2)}, NPR2 ${m.npr2.toFixed(2)})` : ''), '3d', TEXT.shape + '; ' + TEXT.geometry3d, true));
        g3.push(row('extent', 'Longest extent', `${m.extent.nuclei.toFixed(1)} Å between nuclei` + (m.extent.vdw !== null ? `; ${m.extent.vdw.toFixed(1)} Å including van der Waals radii` : ''), '3d', TEXT.geometry3d));
        const vol = m.vdwVolume, area = m.vdwArea;
        g3.push(row('volume', 'van der Waals volume', vol !== null ? `${Math.round(vol)} Å³` : m.radiiNote, '3d', TEXT.volume + '; ' + TEXT.geometry3d, true));
        g3.push(row('area', 'van der Waals surface area', area !== null ? `${Math.round(area)} Å²` : m.radiiNote, '3d', TEXT.area + '; ' + TEXT.geometry3d));
      }
      S.push({ title: '3D model (this conformer)', rows: g3 });
    }
    return S;
  }

  /* ================================================================ 14. labels for the viewers */
  function labelsFor(record, mode) {
    const A = record.atoms;
    switch (mode) {
      case 'element': return A.map(a => a.el);
      case 'index': return A.map((a, i) => a.el + (i + 1));
      case 'cip': return A.map(a => (a.cip ? a.cip : ''));
      case 'hybridization': return A.map(a => (a.isH || !a.hybridization.label ? '' : hybText(a.hybridization.label)));
      case 'formal': return A.map(a => chargeLabel(a.charge));
      case 'partial': return A.map(a => (a.gasteiger === null ? '' : signed(a.gasteiger, 2)));
      case 'oxidation': return A.map(a => (a.oxidationState.value === null ? '' : a.oxidationState.text));
      case 'lonepairs': return A.map(a => (a.isH || a.lonePairs === null ? '' : String(a.lonePairs)));
      default: return A.map(() => '');
    }
  }
  function bondLabelsFor(record, mode) {
    return record.bonds.map(b => {
      if (mode === 'cip' || mode === 'ez') return b.ez === 'E' || b.ez === 'Z' ? b.ez : b.ez === 'unspecified' ? '?' : '';
      if (mode === 'order') return b.aromatic ? (b.bondOrder === 1.5 ? '1.5' : 'ar') : String(b.order);
      return '';
    });
  }

  /* ================================================================ 15. names (window.Chem, round-trip checked) and whole-input totals */
  /* The namer (js/chem.js iupacName) writes neutral closed-shell names: it drops charges and radicals, turns C=N and
     N=N into single bonds, misreads cumulated C=C=X and N–O, and can give a wrong R/S. Its name is therefore shown only
     when (1) the structure has none of those features, (2) the name's R/S/E/Z descriptors are, as a multiset, RDKit's
     CIP labels, and (3) the app's name parser reads the name back to exactly this structure (canonical isomeric SMILES
     compared by RDKit; a structure without stereo marks is compared without stereo).
     How the rule was chosen (a measurement made for this app, 2026-09-28): the 10,929 neutral closed-shell products of
     data/reactions.json plus 32 review cases. The namer named 5,915; 26 lost the exact structure; the guards refused
     443 (C=N/N=N/N=O 224, R/S/E/Z unlike RDKit's 117, N–O 93, charge 6, cumulated 3); of the other 5,446 the parser
     read back 103: 100 the same structure, 3 not (2 wrong names that leave out a C=C, 1 right name to which the parser
     adds an E); 5,343 could not be read back. 2 wrong of 103 read back = 1.9 % > 0.5 %, so a name that cannot be read
     back is not shown either (showUnchecked false): only read-back names are shown. */
  const NAMES_RULE = { showUnchecked: false, measured: { structures: 10961, readBack: 103, wrong: 2, unreadable: 5343 } };
  function namerRefusal(g, tags, name) {                // → reason the namer's name cannot be trusted, or null
    if (g.atoms.some(a => a.chg)) return 'the structure has a formal charge, which the namer does not write';
    if (g.atoms.some(a => a.rad)) return 'the structure is a radical, which the namer does not write';
    if (g.atoms.some(a => a.iso)) return 'the structure has isotope labels, which the namer does not write';
    for (const b of g.bonds) {
      const za = g.atoms[b.a].z, zb = g.atoms[b.b].z, pair = [za, zb].sort((x, y) => x - y).join('-');
      if (!b.arom && b.order === 2 && (pair === '6-7' || pair === '7-7' || pair === '7-8')) return 'the structure has a C=N, N=N or N=O bond, which the namer does not name correctly';
      if (pair === '7-8' && (b.order === 1 || b.arom)) return 'the structure has an N–O bond, which the namer does not name correctly';
    }
    if (g.atoms.some((a, i) => (a.z === 6 || a.z === 7) && g.adj[i].filter(e => g.bonds[e.bond].order === 2 && !g.bonds[e.bond].arom).length >= 2))
      return 'the structure has cumulated double bonds (C=C=X), which the namer does not name correctly';
    const want = new Map(), got = new Map(), add = (m, k) => m.set(k, (m.get(k) || 0) + 1);
    tags.CIP_atoms.forEach(([, lab]) => { const L = lab.replace(/[()]/g, ''); if (L !== '?') add(want, L); });
    tags.CIP_bonds.forEach(([, , lab]) => { const L = lab.replace(/[()]/g, ''); if (L === 'E' || L === 'Z') add(want, L); });
    for (const m of name.matchAll(/\(([^()]*)\)/g)) {
      const toks = m[1].split(',').map(t => t.trim());
      if (!toks.every(t => /^(\d+[a-z]?'*)?[RSEZrs]\*?$/.test(t))) continue;
      toks.forEach(t => add(got, t.replace(/\*$/, '').slice(-1)));
    }
    const keys = new Set([...want.keys(), ...got.keys()]);
    if ([...keys].some(k => (want.get(k) || 0) !== (got.get(k) || 0))) return "its stereodescriptors do not match the structure's CIP labels (RDKit)";
    return null;
  }
  function names(R, M, opts) {
    const out = { iupac: null, iupacNote: null, iupacChecked: null, common: null, cas: null };
    const C = window.Chem;
    if (!C || !C.parseSmiles) { out.iupacNote = "the app's namer is not loaded"; return out; }
    const canon = s => { let m = null; try { m = R.get_mol(s); return m && m.get_num_atoms() ? m.get_smiles() : null; } catch (e) { return null; } finally { if (m) m.delete(); } };
    const smi = M.get_smiles();
    let g;
    try { g = C.parseSmiles(smi); } catch (e) { out.iupacNote = "the app's namer cannot read this structure"; return out; }
    let back = null;
    try { back = canon(C.toSmiles(g)); } catch (e) { back = null; }
    if (back !== smi) { out.iupacNote = "no name: the app's namer does not keep this exact structure (isotopes, radicals, charges or stereo would be lost)"; return out; }
    let nm = null;
    if (smi.includes('.')) out.iupacNote = 'systematic names are given for each component separately';   // the namer writes "C2H3O2 + Na"
    else {
      try { nm = C.iupacName(g); } catch (e) { nm = null; }
      if (nm && !C.nameFailed(nm)) {
        let why = null, parsed = 'unparsed';
        { const J = JSON.parse(M.get_json()), g0 = fromRDKitJson(J); why = namerRefusal(g0, JSON.parse(M.get_stereo_tags()), nm); }
        if (!why) {
          // read the name back with the app's parser: 'same' structure, a 'different' one, or 'unparsed'
          // A structure without any stereo marks is compared without stereo: the parser gives a double bond that the
          // name leaves open an E geometry of its own (the name has no descriptor, and guard 2 found none expected).
          const flat = x => x && x.replace(/\[[^\]]+\]/g, t => t.replace(/@@?/, '')).replace(/[/\\]/g, '');
          try {
            const pg = C.parseIupacName(nm);
            if (pg) {
              let back = canon(C.toSmiles(pg));
              if (back && !/[@/\\]/.test(smi) && /[@/\\]/.test(back)) back = canon(flat(back));
              parsed = back === smi ? 'same' : 'different';
            }
          } catch (e) { parsed = 'unparsed'; }
        }
        if (opts && opts.diagnostics) out.diagnostics = { name: nm, refusal: why, parsed };
        if (why) out.iupacNote = "no name: the app's namer was not used because " + why;
        else if (parsed === 'same') { out.iupac = nm; out.iupacChecked = true; out.iupacNote = TEXT.namer; }
        else if (parsed === 'different') out.iupacNote = "no name: the app's namer's name reads back as a different structure";
        else if (NAMES_RULE.showUnchecked) { out.iupac = nm; out.iupacChecked = false; out.iupacNote = TEXT.namerUnchecked; }
        else out.iupacNote = "no name: the app's namer's name could not be checked (the app's name parser cannot read it back)";
      }
      else out.iupacNote = "the app's namer cannot name this structure" + (nm ? ': ' + nm : '');
    }
    // common name: Chem's everyday name, kept only if that name resolves back to exactly this structure (stereo
    // included: "alanine" is refused for L-alanine); otherwise the library entry's name when the entry's SMILES is
    // this structure and the name is not just a numbered systematic name ("ethan-1-ol"). CAS only from such an entry.
    let common = null;
    try {
      const cn = C.commonNameFor ? C.commonNameFor(g) : null;
      if (cn) { const r = C.searchMolecule(cn); if (r && r.graph && !r.warn && canon(C.toSmiles(r.graph)) === smi) common = cn; }
    } catch (e) { common = null; }
    try {
      const e = C.libEntryFor(g);
      if (e && e.smiles && canon(e.smiles) === smi) {
        if (!common && e.name && !/\d/.test(e.name)) common = e.name;
        out.cas = e.cas || null;
      }
    } catch (e) { /* library lookup is optional */ }
    if (common && out.iupac && common.toLowerCase() === out.iupac.toLowerCase()) common = null;
    out.common = common;
    return out;
  }

  function combined(R, full) {
    const J = JSON.parse(full.get_json()), g = fromRDKitJson(J);
    const desc = JSON.parse(full.get_descriptors());
    const W = weightStats(g, desc);
    let charge = 0, dummy = false;
    g.atoms.forEach(a => { charge += a.chg; if (a.z === 0) dummy = true; });
    const comps = countComponents(g.atoms.length, g.bonds).map(atomsIdx => {
      const counts = {}; let ch = 0, heavy = 0;
      for (const i of atomsIdx) { const a = g.atoms[i], s = a.z ? SYM(a.z) : '*'; counts[s] = (counts[s] || 0) + 1; if (a.hs) counts.H = (counts.H || 0) + a.hs; ch += a.chg; if (a.z > 1) heavy++; }
      return { formula: hill(counts, ch), formulaText: hillText(counts, ch), charge: ch, heavyAtoms: heavy, atoms: atomsIdx };
    });
    const formula = dummy ? null : hill(W.counts, charge);
    return { formula, formulaText: formula ? hillText(W.counts, charge) : null,
      componentsText: comps.map(c => c.formulaText).join(' · '),
      mw: { value: W.mw, rounded: W.rounded, method: W.method, note: W.note }, charge, components: comps };
  }

  window.MolInfo = {
    analyse, geometry, atomRows, bondRows, moleculeSections, labelsFor, bondLabelsFor, names, combined,
    TEXT, KINDS, DIPOLE, FG_TABLE: FG, RD_NOUTER, RD_DV,
    util: { angle, dihedral, klynePrelog, hybText, signed, num, hill, hillText },
    // internals, exposed for tools/test_molinfo.js (bulk comparison with RDKit Python); not for page code
    _internal: { fromRDKitJson, rdkitConjugation, rdkitHybridization, gasteiger, strictRotatable, perceiveGroups, chirality,
      gasteigerDipole, symmetrySearch, dipoleTopology },
  };
})();
