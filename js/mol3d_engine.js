/* Mol3DEngine: 3D models for the Molecule tab. It uses OpenChemLib (vendor/openchemlib/, BSD-3-Clause): the
   conformer generator (torsions and bond lengths from COD crystal statistics), then MMFF94 minimisation.

     Mol3DEngine.build(molblock, {k, timeBudgetMs, seed, onProgress, previous}) → {promise, cancel()}
     Mol3DEngine.energy(molblock3DwithH) → Promise<{ok, energy, error?, singleAtom?, note?}>   single-point MMFF94 energy, kcal/mol
     Mol3DEngine.molfile(record, coords, title, {chiral?}) / xyz(record, coords, title) / sdf(record, conformers, title, {chiral?})
     Mol3DEngine.chiralFlag(record) → the molfile chiral flag molfile() writes by default (see there)
     Mol3DEngine.ready() → Promise<'worker'|'main'>;  Mol3DEngine.status();  Mol3DEngine.options (set before first use)

   Atom order of every result:
   - first, the molblock's atoms in file order, including any H the file lists;
   - then one H per implicit H, grouped by parent in ascending parent index.
   This is the order RDKit's add_hs uses, so the record numbering from js/molinfo.js holds for every index.

   Result = { ok, error?, source, table, seed, N, n, elements, parents (-1 for file atoms), bonds [[a, b, order]] (file
   bonds in file order, then one X–H per added H), charges, isotopes, radicals (radical electrons),
   conformers [{id, coords, energy, relEnergy, optimised, rc, tries, flaggedBonds [[a, b, text]], warnings, ms}],
   stereo {checked, wrong, unspecified}, warnings, atomFlags, generated, potential, singleAtom,
   search {tried, target, stoppedBy, stereoDropped, otherConfiguration, notConverged, approximate, notes},
   merged [[droppedId, keptId]], partial, cancelled?, ms, msFirst, worker }.
   - conformer ids stay the same within a job; merged says which id replaced a duplicate;
   - the promise never rejects;
   - cancel() resolves it at once with what was found so far (cancelled: true);
   - singleAtom: true for a lone atom or ion ([Na+], [Cl-]). Its one conformer is `optimised` (MMFF94 has a type
     for it, and a point has no geometry to change) with energy 0 by definition, and carries
     note: 'single atom: no geometry to optimise', which a page shows instead of "MMFF94-minimised, 0.00 kcal/mol".
     No conformer search is made (search.stoppedBy 'single-atom'). A lone atom MMFF94 cannot type ([He]) is
     "not optimised" as usual, and is also singleAtom;
   - opts.previous (an earlier Result for the same molblock) keeps those conformers and searches on with another seed;
     generated and the search counts then cover both searches;
   - atomFlags [{atom, kind, text, source?}]: every atom whose geometry is approximate (see below); warnings has one line
     per cause; search.stoppedBy 'approximate': MMFF94 cannot optimise this structure (below), so no search is made.

   OpenChemLib was checked against RDKit on 131 inputs before this file was written (tools/test_mol3d_engine.js keeps
   the checks). This file guards what that showed:
   - OCL's MMFF94 energies equal RDKit's to 0.001 kcal/mol, and its atom types matched in 54 of 54 molecules.
   - ForceFieldMMFF94.minimise writes coordinates back only when it returns 0 (converged). Exactly linear start
     geometries (CO2, HCN, acetylene, azide) return 2 until the atoms are nudged.
   - Molecule.toMolfile rescales coordinates when the mean bond length is < 1 Å or > 3 Å, or there are no bonds.
     Coordinates are therefore always read atom by atom.
   - OCL's internal y and z axes are the molfile's, negated. Negating both is a rotation, not a mirror.
   - The generator silently builds (E)-cyclooctene as Z. Each conformer's stereo is read back from its coordinates
     and compared with the input's, and a model that fails is never returned.
   - Adding H to 3D heavy-atom coordinates inverted stereocentres in 4 of 6 tests. Freezing the parities into a 2D
     drawing and rebuilding worked in 6 of 6, so every input goes through that path.
   - Where MMFF94 lacks real parameters (CO2, SO2, SO3, N2O, azide, nitrite, cyanide), minimisation moves a bond
     > 0.17 Å from the generator's crystal-statistics length; nothing else moved more than 0.153 Å. Such bonds are
     flagged. MMFF94 has no carbocation, carbanion or radical types, so those atoms are flagged too.
   - Multi-component inputs are refused. OCL would place the components 2–3 Å apart and include the electrostatics
     between them.
   Found in review (tools/test_mol3d_engine.js has a regression input for each):
   - MMFF94 can converge on a saddle point from a symmetric start: H2O2 at H–O–O–H exactly 180°, 15.30 kcal/mol (the
     gauche minimum is 14.96). Every converged model is nudged and minimised again (minimise()).
   - MMFF94 has no torsion term across an allene's sp carbon: penta-2,3-diene flattens from 90° to 21°, 1,3-dichloro-
     allene to 180° (RDKit's MMFF94 does the same). Such models are set back and marked approximate (optimise()).
   - SF4-type atoms (a lone pair and 4 or more σ bonds) are typed like a sulfone S and built as a tetrahedron, and
     SiF6²⁻ is "minimised" to a distorted shape: every atom with a lone pair and steric number ≥ 5, or with ≥ 5 σ
     bonds, is marked approximate without running MMFF94 (hypervalent()). A lone centre with identical terminal
     ligands (I3−, XeF2, SF4, SF6) is given its ideal VSEPR shape (idealiseVSEPR()).
   - I3− was refused as an "atom-mapping" failure: OpenChemLib's molfile round trip adds an H to a charged atom with an
     abnormal valence (keepValence()).
   - A lone atom has no geometry to optimise: Result.singleAtom, conformer.note, no search.
   - A charged or radical carbon in a conjugated system spreads over the whole system (tropylium), so the whole
     system is flagged, not one atom. A lone atom MMFF94 cannot type (He, H+) is "not optimised", not "energy 0".
   - OCL lists norbornane's bridgeheads as open stereocentres; atoms that are not a real choice are dropped from the
     "not specified" warning (realOpen()).
   Found while building this file (tools/test_mol3d_engine.js):
   - A molfile with chiral flag 0 (RDKit's) makes OCL read every stereocentre as racemic. Its fragment cache then
     returned the first-built enantiomer for the second, so stereocentres are set to absolute before building.
   - inventCoordinates drops H atoms that a file lists; they are re-added and mapped back to their file positions.
   - Where the input leaves stereo open, the first model's configuration is kept for every conformer. Otherwise the
     "conformers" of CC=CC would include both E and Z.

   The work runs in a Web Worker made from a Blob: OCL_FACTORY.toString() (the whole library as text) plus the
   engine core below. The page stays responsive during multi-second builds. If no Worker starts, the same core runs
   on the main thread and yields between conformers. */
window.Mol3DEngine = (() => {
  'use strict';

  /* ================= engine core ================= */
  // Self-contained: the Worker rebuilds it from this function's source text, so it may use only OCL and built-ins.
  function engineCore(OCL) {
    'use strict';
    const Mol = OCL.Molecule, CG = OCL.ConformerGenerator, FF = OCL.ForceFieldMMFF94;
    const TABLE = 'MMFF94', MAX_ITS = 4000, TRIES = 3;
    const GAP = 0.17;                 // Å, the MMFF94 parameter-gap threshold (see header)
    const E_TOL = 0.01, RMSD_CUT = 0.3;  // duplicate conformers: energy within 0.01 kcal/mol, or heavy-atom RMSD < 0.3 Å
    // saddle escape: up to 3 re-minimisations from a copy nudged by up to ±0.05 Å per coordinate; a result counts as
    // lower when it drops > 0.001 kcal/mol (below that it is the minimiser's own tolerance), and the escape goes on
    // while it drops > 0.01
    const SADDLE_TRIES = 3, SADDLE_MIN = 0.001, SADDLE_DROP = 0.01, NUDGE = 0.1;
    const CUMULENE_TOL = 10;          // degrees: end planes of an allene within 90 ± 10°, of a butatriene within 0 ± 10° (or 180)
    const MAPFAIL = 'internal atom-mapping check failed';
    const SINGLE_NOTE = 'single atom: no geometry to optimise';
    const RAD_E = { 16: 2, 32: 1, 48: 2 };  // OCL radical state (singlet, doublet, triplet) → radical electrons (RDKit's count)
    const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());
    const errText = e => String((e && e.message) || e).replace(/^Class\$S\d+: /, '');
    let registered = false;

    // valence electrons of a main-group element (groups 1, 2, 13–18); null for H, He and the d and f blocks
    function valenceElectrons(z) {
      if (z >= 3 && z <= 10) return z - 2;
      if (z >= 11 && z <= 18) return z - 10;
      if (z === 19 || z === 37 || z === 55) return 1;
      if (z === 20 || z === 38 || z === 56) return 2;
      if (z >= 31 && z <= 36) return z - 28;
      if (z >= 49 && z <= 54) return z - 46;
      if (z >= 81 && z <= 86) return z - 78;
      return null;
    }

    function register(resources) {
      if (registered) return;
      if (!resources) throw new Error('OpenChemLib data tables missing (openchemlib-resources.js)');
      OCL.Resources.register(resources);
      registered = true;
    }

    // mulberry32: the nudges come from a seeded generator, so the same seed repeats the same result
    function rng(seed) {
      let a = seed >>> 0;
      return () => {
        a = (a + 0x6D2B79F5) | 0;
        let t = Math.imul(a ^ (a >>> 15), 1 | a);
        t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
      };
    }

    const dist = (m, i, j) => Math.hypot(m.getAtomX(i) - m.getAtomX(j), m.getAtomY(i) - m.getAtomY(j), m.getAtomZ(i) - m.getAtomZ(j));
    function bondLengths(m) {
      const r = new Float64Array(m.getAllBonds());
      for (let b = 0; b < r.length; b++) r[b] = dist(m, m.getBondAtom(0, b), m.getBondAtom(1, b));
      return r;
    }
    function xyzOf(m) {
      const n = m.getAllAtoms(), a = new Float64Array(3 * n);
      for (let i = 0; i < n; i++) { a[3 * i] = m.getAtomX(i); a[3 * i + 1] = m.getAtomY(i); a[3 * i + 2] = m.getAtomZ(i); }
      return a;
    }
    function setXYZ(m, a) {
      for (let i = 0; i < m.getAllAtoms(); i++) { m.setAtomX(i, a[3 * i]); m.setAtomY(i, a[3 * i + 1]); m.setAtomZ(i, a[3 * i + 2]); }
    }
    function nudge(m, rand) {
      for (let i = 0; i < m.getAllAtoms(); i++) {
        m.setAtomX(i, m.getAtomX(i) + (rand() - 0.5) * NUDGE);
        m.setAtomY(i, m.getAtomY(i) + (rand() - 0.5) * NUDGE);
        m.setAtomZ(i, m.getAtomZ(i) + (rand() - 0.5) * NUDGE);
      }
    }
    // signed dihedral i–j–k–l in degrees (IUPAC sign, as RDKit's GetDihedralDeg), from flat coordinates X
    function dihedral(X, i, j, k, l) {
      const v = (p, q) => [X[3 * q] - X[3 * p], X[3 * q + 1] - X[3 * p + 1], X[3 * q + 2] - X[3 * p + 2]];
      const cr = (u, w) => [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
      const dot = (u, w) => u[0] * w[0] + u[1] * w[1] + u[2] * w[2];
      const b1 = v(i, j), b2 = v(j, k), b3 = v(k, l), n1 = cr(b1, b2), n2 = cr(b2, b3), len = Math.hypot(b2[0], b2[1], b2[2]) || 1;
      return Math.atan2(dot(cr(n1, n2), b2) / len, dot(n1, n2)) * 180 / Math.PI;
    }

    /* OpenChemLib's molfile writer and reader disagree about the valence field of a charged atom with an abnormal
       valence: the central I of I3− (I[I-]I) is written as valence 3 and read back as 4, which adds a hydrogen
       (I[IH-]I); ICl2−, IF4−, ICl4− and SiF6²⁻ do the same. After every molfile round trip inside the engine (same atom
       indices on both sides) the valence is copied back from the molecule that was written. */
    function keepValence(from, to) {
      const n = Math.min(from.getAllAtoms(), to.getAllAtoms());
      for (let j = 0; j < n; j++) {
        const v = from.getAtomAbnormalValence(j);
        if (to.getAtomAbnormalValence(j) !== v) to.setAtomAbnormalValence(j, v);
      }
      return to;
    }

    // OCL's molfile with the real coordinates written back (toMolfile may have rescaled them). Only used to read
    // stereo parities back from a model; the page's coordinates are read atom by atom.
    function molfile3D(m) {
      const lines = m.toMolfile().split('\n'), f = v => v.toFixed(4).padStart(10);
      if (!/V2000/.test(lines[3] || '')) throw new Error('unexpected molfile layout from OpenChemLib');
      for (let i = 0; i < m.getAllAtoms(); i++) lines[4 + i] = f(m.getAtomX(i)) + f(-m.getAtomY(i)) + f(-m.getAtomZ(i)) + lines[4 + i].slice(30);
      return lines.join('\n');
    }

    // The molfile's own bond list, in file order (the record's bond order). null if the text is not laid out as expected.
    function fileBondList(text, nAtoms) {
      const lines = String(text).split(/\r\n|\r|\n/);
      let c = -1;
      for (let i = 0; i < Math.min(lines.length, 6); i++) if (/V[23]000\s*$/.test(lines[i])) { c = i; break; }
      if (c < 0) return null;
      if (/V3000/.test(lines[c])) {
        const v = []; let buf = '';
        for (let i = c + 1; i < lines.length; i++) {
          let l = lines[i];
          if (/^M  END/.test(l)) break;
          if (!l.startsWith('M  V30 ')) continue;
          l = buf + l.slice(7); buf = '';
          if (/-\s*$/.test(l)) { buf = l.replace(/-\s*$/, ''); continue; }   // continuation line
          v.push(l.trim());
        }
        const out = []; let inBonds = false;
        for (const l of v) {
          if (/^BEGIN BOND/.test(l)) inBonds = true;
          else if (/^END BOND/.test(l)) break;
          else if (inBonds) { const f = l.split(/\s+/); out.push([+f[2] - 1, +f[3] - 1]); }
        }
        return out;
      }
      const na = parseInt(lines[c].slice(0, 3), 10), nb = parseInt(lines[c].slice(3, 6), 10);
      if (na !== nAtoms || !(nb >= 0) || lines.length < c + 1 + na + nb) return null;
      const out = [];
      for (let k = 0; k < nb; k++) {
        const l = lines[c + 1 + na + k];
        out.push([parseInt(l.slice(0, 3), 10) - 1, parseInt(l.slice(3, 6), 10) - 1]);
      }
      return out;
    }

    function isPermutation(map, n) {
      if (map.length !== n) return false;
      const seen = new Uint8Array(n);
      for (const v of map) { if (!(v >= 0 && v < n) || seen[v]) return false; seen[v] = 1; }
      return true;
    }

    /* Everything that does not depend on coordinates: parse, the frozen-parity 2D rebuild, the H-complete molecule
       every conformer starts from, and the output numbering. Model index space: the parsed molecule's (m0) atoms other
       than plain hydrogens keep their index; every plain H (implicit, or listed in the file) is re-added after them.
       Each step's effect on the atom list is checked, never assumed. */
    function prepare(text) {
      let parsed;
      try { parsed = Mol.fromMolfileWithAtomMap(String(text)); } catch (e) { return { error: 'could not read the structure: ' + errText(e) }; }
      const m0 = parsed.molecule, map0 = Array.from(parsed.map || []), n0 = m0.getAllAtoms();
      if (!n0) return { error: 'empty structure: no atoms' };
      if (!isPermutation(map0, n0)) return { error: MAPFAIL + ' (the molfile reader changed the atom list)' };
      if (m0.getFragmentNumbers([], false, false) > 1) return { error: 'several separate components: a 3D model is built for one component at a time' };
      // The configuration drawn is the one to build. A molfile with chiral flag 0 (RDKit always writes 0) makes OCL mark
      // every stereocentre "AND" (racemic), and the generator's fragment cache then treats R and S fragments as one:
      // after (R)-butan-2-ol, (S)-butan-2-ol came out R. Absolute ESR types avoid that.
      m0.ensureHelperArrays(Mol.cHelperParities);
      for (let j = 0; j < n0; j++) if (m0.getAtomESRType(j) !== Mol.cESRTypeAbs) m0.setAtomESR(j, Mol.cESRTypeAbs, 0);
      for (let b = 0; b < m0.getAllBonds(); b++) if (m0.getBondESRType(b) !== Mol.cESRTypeAbs) m0.setBondESR(b, Mol.cESRTypeAbs, 0);
      m0.ensureHelperArrays(Mol.cHelperParities);         // parities from the file (wedges in 2D, coordinates in 3D)
      const inv0 = new Int32Array(n0);
      map0.forEach((j, i) => { inv0[j] = i; });
      let missingH = 0;
      for (let j = 0; j < n0; j++) missingH += m0.getImplicitHydrogens(j);
      const fromFile = m0.is3D() && missingH === 0;       // 3D with every H: its own coordinates are the first conformer

      const fb = fileBondList(text, n0);
      let fileBonds = fb;
      if (fb) {
        if (fb.length !== m0.getAllBonds()) return { error: MAPFAIL + ' (bond count)' };
        for (const [a, b] of fb) if (!(a >= 0 && a < n0 && b >= 0 && b < n0) || m0.getBond(map0[a], map0[b]) < 0) return { error: MAPFAIL + ' (bond list)' };
      } else {
        fileBonds = [];
        for (let b = 0; b < m0.getAllBonds(); b++) fileBonds.push([inv0[m0.getBondAtom(0, b)], inv0[m0.getBondAtom(1, b)]]);
      }

      // freeze the parities into a fresh 2D drawing, then re-read it (parities now come from its wedges);
      // inventCoordinates drops plain H atoms listed in the file: they come back with the implicit ones below
      const flat = m0.getCompactCopy();
      flat.ensureHelperArrays(Mol.cHelperParities);
      flat.inventCoordinates(); flat.setStereoBondsFromParity(); flat.setParitiesValid(0);
      let r1;
      try { r1 = Mol.fromMolfileWithAtomMap(flat.toMolfile()); } catch (e) { return { error: 'could not prepare the structure: ' + errText(e) }; }
      const ref = r1.molecule, nRef = ref.getAllAtoms();
      if (nRef > n0 || Array.from(r1.map).some((v, i) => v !== i)) return { error: MAPFAIL + ' (the 2D rebuild changed the atom order)' };
      keepValence(flat, ref);                             // I3−: see keepValence()
      for (let j = 0; j < nRef; j++) {
        if (ref.getAtomicNo(j) !== m0.getAtomicNo(j) || ref.getAtomMass(j) !== m0.getAtomMass(j) || ref.getAtomCharge(j) !== m0.getAtomCharge(j) ||
            ref.getAtomRadical(j) !== m0.getAtomRadical(j)) return { error: MAPFAIL + ' (the 2D rebuild changed an atom)' };
      }
      m0.ensureHelperArrays(Mol.cHelperRings);
      const fileH = new Int32Array(nRef);                  // plain H listed in the file, per parent
      for (let j = nRef; j < n0; j++) {
        if (m0.getAtomicNo(j) !== 1 || m0.getAllConnAtoms(j) !== 1 || m0.getConnAtom(j, 0) >= nRef) return { error: MAPFAIL + ' (the 2D rebuild dropped an atom that is not a plain hydrogen)' };
        fileH[m0.getConnAtom(j, 0)]++;
      }
      ref.ensureHelperArrays(Mol.cHelperParities);

      const withH = ref.getCompactCopy();
      withH.addImplicitHydrogens();
      withH.ensureHelperArrays(Mol.cHelperRings);          // getConnAtom on hydrogens needs the ring helpers
      const N = withH.getAllAtoms();
      if (N !== n0 + missingH) return { error: MAPFAIL + ' (hydrogen count)' };
      for (let j = 0; j < nRef; j++) if (withH.getAtomicNo(j) !== m0.getAtomicNo(j) || withH.getAtomMass(j) !== m0.getAtomMass(j)) return { error: MAPFAIL + ' (adding hydrogens moved an atom)' };
      const addedOn = Array.from({ length: nRef }, () => []);
      for (let h = nRef; h < N; h++) {
        if (withH.getAtomicNo(h) !== 1 || withH.getAllConnAtoms(h) !== 1) return { error: MAPFAIL + ' (added hydrogen)' };
        const p = withH.getConnAtom(h, 0);
        if (p >= nRef) return { error: MAPFAIL + ' (added hydrogen parent)' };
        addedOn[p].push(h);
      }
      for (let p = 0; p < nRef; p++) if (addedOn[p].length !== m0.getImplicitHydrogens(p) + fileH[p]) return { error: MAPFAIL + ' (hydrogens per atom)' };

      // output order: file atoms in file order (a listed H takes one of the H re-added on its parent), then the other
      // added H grouped by parent, ascending parent index in output order
      const out2model = new Int32Array(N), model2out = new Int32Array(N), parents = new Int32Array(N).fill(-1);
      for (let i = 0; i < n0; i++) {
        const j = map0[i];
        out2model[i] = j < nRef ? j : addedOn[m0.getConnAtom(j, 0)].shift();
      }
      const rest = [];
      for (let p = 0; p < nRef; p++) for (const h of addedOn[p]) rest.push([inv0[p], h]);
      rest.sort((x, y) => x[0] - y[0] || x[1] - y[1]);
      rest.forEach(([po, h], k) => { out2model[n0 + k] = h; parents[n0 + k] = po; });
      model2out.fill(-1);
      for (let o = 0; o < N; o++) model2out[out2model[o]] = o;
      if (model2out.some(v => v < 0)) return { error: MAPFAIL + ' (output numbering)' };
      const elements = [], charges = new Int32Array(N), isotopes = new Int32Array(N), radicals = new Int32Array(N), heavy = [];
      for (let o = 0; o < N; o++) {
        const j = out2model[o];
        elements.push(withH.getAtomLabel(j));
        if (o < n0) { const j0 = map0[o]; charges[o] = m0.getAtomCharge(j0); isotopes[o] = m0.getAtomMass(j0); radicals[o] = RAD_E[m0.getAtomRadical(j0)] || 0; }
        if (withH.getAtomicNo(j) > 1) heavy.push(o);
      }
      const bonds = fileBonds.map(([a, b]) => [a, b, m0.getBondOrder(m0.getBond(map0[a], map0[b]))]);
      rest.forEach(([po], k) => bonds.push([po, n0 + k, 1]));
      for (const [a, b] of bonds) if (withH.getBond(out2model[a], out2model[b]) < 0) return { error: MAPFAIL + ' (bond missing in the model)' };
      if (bonds.length !== withH.getAllBonds()) return { error: MAPFAIL + ' (model bond count)' };
      const lab = o => elements[o] + (o + 1);             // "C3": 1-based, like a molfile

      // stereo the input leaves open (the model then shows one arbitrary configuration)
      const unspecified = [], open = { atoms: [], bonds: [] };      // output indices / model indices
      for (let j = 0; j < m0.getAtoms(); j++) {
        const p = m0.getAtomParity(j);
        if (p === 3 || (p === 0 && m0.isAtomStereoCenter(j))) { unspecified.push({ atom: inv0[j] }); open.atoms.push(j); }
      }
      for (let b = 0; b < m0.getBonds(); b++) {
        if (m0.getBondParity(b) !== 3) continue;
        const a0 = m0.getBondAtom(0, b), a1 = m0.getBondAtom(1, b);
        unspecified.push({ bond: [inv0[a0], inv0[a1]] }); open.bonds.push([a0, a1]);
      }

      // Atoms MMFF94 has no types for (radical, carbocation, carbanion): its geometry there is approximate. In a
      // conjugated or aromatic system the charge or odd electron is spread over the whole system (tropylium: MMFF94 gives
      // every ring C–C 1.34 Å, not ~1.40), so every atom of that system is flagged, not only the one drawn charged.
      const sources = [];
      for (let j = 0; j < n0; j++) {
        if (m0.getAtomRadical(j)) sources.push({ j, kind: 'radical', why: 'MMFF94 has no radical atom types' });
        else if (m0.getAtomicNo(j) === 6 && m0.getAtomCharge(j) !== 0 && m0.getAtomPi(j) === 0)
          sources.push({ j, kind: m0.getAtomCharge(j) > 0 ? 'carbocation' : 'carbanion', why: 'MMFF94 has no parameters for it' });
      }
      const isSource = new Set(sources.map(s => s.j)), piAtom = j => m0.getAtomPi(j) > 0 || m0.isAromaticAtom(j) || isSource.has(j);
      const atomFlags = [], flagWarnings = [], flagged = new Set();
      const labList = os => (os.length > 10 ? os.slice(0, 10).map(lab).join(', ') + ` and ${os.length - 10} more` : os.map(lab).join(', '));
      for (const s of sources) {
        const sys = [s.j], seen = new Set(sys);          // the conjugated system: out through atoms with π bonds
        for (let q = 0; q < sys.length; q++) {
          for (let k = 0; k < m0.getConnAtoms(sys[q]); k++) {
            const y = m0.getConnAtom(sys[q], k);
            if (!seen.has(y) && piAtom(y)) { seen.add(y); sys.push(y); }
          }
        }
        const o = inv0[s.j], others = sys.slice(1).map(j => inv0[j]).sort((x, y) => x - y);
        const text = others.length
          ? `${lab(o)} is a ${s.kind} in a conjugated system (${labList([o].concat(others).sort((x, y) => x - y))}): ${s.why}, so the geometry of the whole system is approximate`
          : `${lab(o)} is a ${s.kind}: ${s.why}, so the geometry around it is approximate`;
        flagWarnings.push(text);
        if (!flagged.has(o)) { flagged.add(o); atomFlags.push({ atom: o, kind: s.kind, text }); }
        for (const x of others) {
          if (flagged.has(x)) continue;
          flagged.add(x);
          atomFlags.push({ atom: x, kind: 'conjugated', source: o, text: `${lab(x)} is in the conjugated system of the ${s.kind} ${lab(o)}: ${s.why}, so the geometry here is approximate` });
        }
      }
      atomFlags.sort((x, y) => x.atom - y.atom);

      // model index = m0 index for every atom that is not a plain H
      const hyperInfo = hypervalent(m0).map(h => ({ out: inv0[h.j], model: h.j, sigma: h.sigma, lp: h.lp })).sort((x, y) => x.out - y.out);
      const hyper = hyperInfo.map(h => h.out);

      // Cumulated double bonds at sp carbons (allenes, ketenimines, carbodiimides, butatrienes), in model indices.
      // MMFF94 has no torsion term across an sp centre, so minimisation turns the end planes freely (penta-2,3-diene
      // 90° → 21°, 1,3-dichloroallene → 180°); they are checked after every minimisation (optimise()).
      const cumulenes = [], inChain = new Set();
      const isCentre = j => j < nRef && withH.getAtomicNo(j) === 6 && withH.getAtomCharge(j) === 0 && withH.getAllConnAtoms(j) === 2 &&
        withH.getConnBondOrder(j, 0) === 2 && withH.getConnBondOrder(j, 1) === 2;
      for (let c = 0; c < nRef; c++) {
        if (!isCentre(c) || inChain.has(c)) continue;
        const chain = [c];
        inChain.add(c);
        for (const dir of [0, 1]) {                        // walk out to both ends of the chain
          let prev = c, cur = withH.getConnAtom(c, dir);
          while (isCentre(cur) && !inChain.has(cur)) {
            inChain.add(cur);
            if (dir) chain.push(cur); else chain.unshift(cur);
            const nx = withH.getConnAtom(cur, 0) === prev ? withH.getConnAtom(cur, 1) : withH.getConnAtom(cur, 0);
            prev = cur; cur = nx;
          }
          if (dir) chain.push(cur); else chain.unshift(cur);
        }
        const a = chain[0], b = chain[chain.length - 1];
        if (a === b || isCentre(a) || isCentre(b)) continue;
        const subs = (e, nb) => { const s = []; for (let k = 0; k < withH.getAllConnAtoms(e); k++) if (withH.getConnAtom(e, k) !== nb) s.push(withH.getConnAtom(e, k)); return s; };
        const sa = subs(a, chain[1]), sb = subs(b, chain[chain.length - 2]);
        if (!sa.length || !sb.length) continue;            // ketene, isocyanate, CO2: nothing to turn
        // the atoms on b's side, turned as one piece to set the ends back; null when a ring joins the two ends
        let side = [];
        const seen = new Set(chain.concat(sb)), stack = sb.slice();
        while (stack.length && side) {
          const x = stack.pop();
          side.push(x);
          for (let k = 0; k < withH.getAllConnAtoms(x); k++) {
            const y = withH.getConnAtom(x, k);
            if (y === b && sb.includes(x)) continue;
            if (chain.includes(y)) { side = null; break; }
            if (!seen.has(y)) { seen.add(y); stack.push(y); }
          }
        }
        const nDouble = chain.length - 1, allC = chain.every(j => withH.getAtomicNo(j) === 6);
        // does the configuration matter (axial chirality, or E/Z)? Only when neither end carries two equal substituents
        if (!cumulenes.length) withH.ensureHelperArrays(Mol.cHelperSymmetrySimple);
        const unequal = s2 => s2.length === 1 || withH.getSymmetryRank(s2[0]) !== withH.getSymmetryRank(s2[1]);
        // an even number of cumulated double bonds (allene: 2) has perpendicular end planes, an odd one (butatriene: 3) one plane
        cumulenes.push({ chain, a, b, sa: sa[0], sb: sb[0], side, perpendicular: nDouble % 2 === 0, stereo: unequal(sa) && unequal(sb),
                         name: nDouble === 2 ? (allC ? 'allene' : 'heteroallene') : 'cumulene',
                         label: chain.map(j => lab(model2out[j])).join('=') });
      }
      // OpenChemLib keeps an allene's configuration (a parity on its central atom) but not a butatriene's E/Z or a
      // pentatetraene's axis: the generator picks one, so say so (and search() keeps the first model's choice)
      const tracked = q => q.chain.some((j, i) => (i > 0 && i < q.chain.length - 1 && (m0.isAtomStereoCenter(j) || m0.getAtomParity(j) !== 0)) ||
        (i > 0 && m0.getBond(q.chain[i - 1], j) >= 0 && m0.getBondParity(m0.getBond(q.chain[i - 1], j)) !== 0));
      const cumuleneWarnings = cumulenes.filter(q => q.stereo && !tracked(q))
        .map(q => `${q.perpendicular ? 'the axial configuration of' : 'E/Z at'} the ${q.name} ${q.label} is not taken from the input: the model shows the configuration the builder chose`);
      return { m0, map0, inv0, ref, withH, n0, N, fromFile, out2model, model2out, parents, elements, charges, isotopes, radicals, bonds, heavy, lab,
               unspecified, open, atomFlags, flagWarnings, hyper, hyperInfo, cumulenes, cumuleneWarnings };
    }

    /* Hypervalent atoms: with a lone pair and VSEPR steric number ≥ 5 (SF4, ClF3, XeF2, I3−, PhI(OAc)2), or with 5 or
       more σ bonds (PCl5, SF6, PF6−, SiF6²⁻). MMFF94 has no type for either: OpenChemLib (like RDKit) types SF4's S as a
       sulfone-like S and builds a regular tetrahedron, not the seesaw, and it types SiF6²⁻'s Si and "minimises" it to
       a shape with 73°, 88° and 133° F–Si–F angles (RDKit's MMFF94 gives the same energy, 227.5 kcal/mol), not the
       octahedron. Steric number = σ neighbours + lone pairs, the lone pairs from valence electrons, bond orders and
       charge (as js/molinfo.js counts them). Returns [{j, sigma, lp}]. */
    function hypervalent(m) {
      const out = [];
      m.ensureHelperArrays(Mol.cHelperNeighbours);
      for (let j = 0; j < m.getAllAtoms(); j++) {
        const V = valenceElectrons(m.getAtomicNo(j));
        if (V == null) continue;
        let sigma = m.getImplicitHydrogens(j), B = sigma;
        for (let k = 0; k < m.getAllConnAtoms(j); k++) { sigma++; B += m.getConnBondOrder(j, k); }
        const free = V - m.getAtomCharge(j) - B - (RAD_E[m.getAtomRadical(j)] || 0);
        const lp = free >= 2 && free % 2 === 0 ? free / 2 : 0;
        if ((lp && sigma + lp >= 5) || sigma >= 5) out.push({ j, sigma, lp });
      }
      return out;
    }
    // "S2, a hypervalent atom with a lone pair" / "P2, a 6-coordinate atom" (h: [{name, sigma, lp}])
    function hyperText(h) {
      const kind = x => (x.lp ? 'a hypervalent atom with a lone pair' : `a ${x.sigma}-coordinate atom`);
      if (h.length === 1) return `${h[0].name}, ${kind(h[0])}`;
      if (h.every(x => x.lp)) return `${h.map(x => x.name).join(' and ')}, hypervalent atoms with a lone pair`;
      return h.map(x => `${x.name} (${kind(x)})`).join(' and ');
    }

    /* The ideal VSEPR shape for the simplest hypervalent species: one centre whose σ neighbours are all terminal atoms
       of one element, charge and bond order (I3−, ICl2−, XeF2, SF4, ClF3, XeF4, IF4−, BrF5, SF6, PF6−, SiF6²⁻). Lone
       pairs take the equatorial sites of a trigonal bipyramid and trans sites of an octahedron; the generator's bond
       lengths are kept. m: the model (every atom); c: the centre (model index). Returns the shape's name, or null
       when the species is not of that kind (the generator's geometry is then shown, as approximate). */
    function idealiseVSEPR(m, c, sigma, lp) {
      const n = m.getAllAtoms(), L = [];
      if (m.getAllConnAtoms(c) !== sigma || n !== sigma + 1) return null;
      for (let k = 0; k < sigma; k++) L.push(m.getConnAtom(c, k));
      const same = j => m.getAtomicNo(j) === m.getAtomicNo(L[0]) && m.getAtomCharge(j) === m.getAtomCharge(L[0]) &&
        m.getAtomMass(j) === m.getAtomMass(L[0]) && m.getAllConnAtoms(j) === 1 && m.getBondOrder(m.getBond(c, j)) === m.getBondOrder(m.getBond(c, L[0]));
      if (!L.every(same)) return null;
      const ax = [[0, 0, 1], [0, 0, -1]], eq3 = [0, 1, 2].map(k => [Math.cos(2 * Math.PI * k / 3), Math.sin(2 * Math.PI * k / 3), 0]);
      const eq5 = [0, 1, 2, 3, 4].map(k => [Math.cos(2 * Math.PI * k / 5), Math.sin(2 * Math.PI * k / 5), 0]);
      const oct = { px: [1, 0, 0], mx: [-1, 0, 0], py: [0, 1, 0], my: [0, -1, 0], pz: [0, 0, 1], mz: [0, 0, -1] };
      const SHAPES = {                                    // steric number: {ligands: [name, directions]}
        5: { 5: ['trigonal bipyramidal', ax.concat(eq3)], 4: ['seesaw', ax.concat(eq3.slice(0, 2))], 3: ['T-shaped', ax.concat(eq3.slice(0, 1))], 2: ['linear', ax] },
        6: { 6: ['octahedral', Object.values(oct)], 5: ['square pyramidal', [oct.pz, oct.px, oct.mx, oct.py, oct.my]], 4: ['square planar', [oct.px, oct.mx, oct.py, oct.my]],
             3: ['T-shaped', [oct.px, oct.mx, oct.py]], 2: ['linear', [oct.px, oct.mx]] },
        7: { 7: ['pentagonal bipyramidal', ax.concat(eq5)], 5: ['pentagonal planar', eq5] },
      };
      const s = (SHAPES[sigma + lp] || {})[sigma];
      if (!s) return null;
      const cx = m.getAtomX(c), cy = m.getAtomY(c), cz = m.getAtomZ(c);
      L.forEach((j, k) => {
        const r = dist(m, c, j) || 1.5, d = s[1][k];
        m.setAtomX(j, cx + r * d[0]); m.setAtomY(j, cy + r * d[1]); m.setAtomZ(j, cz + r * d[2]);
      });
      return s[0];
    }

    // MMFF94 minimisation. Coordinates change only when rc = 0; exactly linear start geometries give rc 2 until nudged.
    // Then the saddle escape: a symmetric start can converge on a stationary point that is not a minimum (H2O2 from
    // the generator: H–O–O–H exactly 180°, 15.30 kcal/mol; the minimum is gauche, 14.96). So after convergence a nudged
    // copy is minimised again and kept when its energy is lower (see SADDLE_* above).
    function minimise(m, rand) {
      if (m.getAllAtoms() < 2) {
        // a lone atom: nothing to optimise (minimise() itself returns 2 for Na+), but MMFF94 must have a type for it:
        // He, U and H+ have none, and are then "not optimised" like any untypable atom
        try { new FF(m, TABLE); } catch (e) { return { ok: false, rc: null, energy: NaN, error: errText(e) }; }
        return { ok: true, rc: 0, energy: 0, tries: 0, single: true };
      }
      let ff;
      try { ff = new FF(m, TABLE); } catch (e) { return { ok: false, rc: null, energy: NaN, error: errText(e) }; }
      let rc = ff.minimise({ maxIts: MAX_ITS }), t = 1;
      for (; rc !== 0 && t < TRIES; t++) { nudge(m, rand); ff = new FF(m, TABLE); rc = ff.minimise({ maxIts: MAX_ITS }); }
      if (rc !== 0) return { ok: false, rc, energy: NaN, tries: t };
      let e = ff.getTotalEnergy(), escapes = 0;           // (equals a fresh field's energy on these coordinates)
      for (let s = 0; s < SADDLE_TRIES; s++) {
        const keep = xyzOf(m);
        nudge(m, rand);
        const f2 = new FF(m, TABLE), rc2 = f2.minimise({ maxIts: MAX_ITS }), e2 = rc2 === 0 ? f2.getTotalEnergy() : NaN;
        if (!(e2 < e - SADDLE_MIN)) { setXYZ(m, keep); break; }   // not lower (or not converged): back to the kept minimum
        escapes++;
        const drop = e - e2;
        e = e2;
        if (drop <= SADDLE_DROP) break;
      }
      const energy = new FF(m, TABLE).getTotalEnergy();   // kcal/mol, fresh field on the final coordinates
      return { ok: isFinite(energy), rc: 0, energy, tries: t, escapes };
    }

    // Angle between the end planes of a cumulene (0–90°), each plane through the axis and the end's first substituent.
    function endPlanes(X, q) { const d = Math.abs(dihedral(X, q.sa, q.a, q.b, q.sb)); return Math.min(d, 180 - d); }
    const cumuleneOK = (X, q) => Math.abs(endPlanes(X, q) - (q.perpendicular ? 90 : 0)) <= CUMULENE_TOL;
    // a model's configuration at each cumulene whose ends are unequal: the sense of twist (allene) or cis/trans (butatriene)
    const cumuleneConfig = (P, X) => P.cumulenes.map(q => { if (!q.stereo) return 0; const d = dihedral(X, q.sa, q.a, q.b, q.sb); return q.perpendicular ? Math.sign(d) : (Math.abs(d) < 90 ? 1 : -1); }).join(',');
    // Turn the atoms on q's b side about the a→b axis so that the dihedral sa–a–b–sb becomes `target` degrees.
    function turnTo(m, q, target) {
      const rot = deg => {
        const X = xyzOf(m), ax = [X[3 * q.b] - X[3 * q.a], X[3 * q.b + 1] - X[3 * q.a + 1], X[3 * q.b + 2] - X[3 * q.a + 2]];
        const L = Math.hypot(ax[0], ax[1], ax[2]), k = ax.map(x => x / L), t = deg * Math.PI / 180, c = Math.cos(t), s = Math.sin(t);
        for (const i of q.side) {
          const v = [X[3 * i] - X[3 * q.b], X[3 * i + 1] - X[3 * q.b + 1], X[3 * i + 2] - X[3 * q.b + 2]];
          const kv = k[0] * v[0] + k[1] * v[1] + k[2] * v[2], kx = [k[1] * v[2] - k[2] * v[1], k[2] * v[0] - k[0] * v[2], k[0] * v[1] - k[1] * v[0]];
          m.setAtomX(i, X[3 * q.b] + v[0] * c + kx[0] * s + k[0] * kv * (1 - c));
          m.setAtomY(i, X[3 * q.b + 1] + v[1] * c + kx[1] * s + k[1] * kv * (1 - c));
          m.setAtomZ(i, X[3 * q.b + 2] + v[2] * c + kx[2] * s + k[2] * kv * (1 - c));
        }
      };
      const off = () => { const d = dihedral(xyzOf(m), q.sa, q.a, q.b, q.sb) - target; return ((d % 360) + 540) % 360 - 180; };
      const delta = -off();
      rot(delta);
      if (Math.abs(off()) > 1) rot(-2 * delta);            // the other sense of rotation
    }

    /* minimise(), then the cumulene check. MMFF94 has no torsion term across an sp carbon, so an allene's end planes
       (90° apart in reality) turn freely, towards whatever the other terms prefer: RDKit's MMFF94 gives the same. When
       they end up more than 10° off, the b end is turned back (keeping the sense of twist the generator built, i.e. the
       configuration of an axially chiral allene) and minimised again. If they stay (tetramethylallene: MMFF94's twist
       surface is flat there, 0.03 kcal/mol from 21° to 90°), the model is a converged MMFF94 structure with the right
       shape; if they drift off again (penta-2,3-diene, 1,3-dichloroallene, carbodiimides), the turned-back model is kept
       and marked approximate, not force-field optimised. A ring through both ends cannot be turned: the generator's own
       geometry is kept, also marked approximate. */
    function optimise(m, rand, P, pre) {
      const r = minimise(m, rand);
      if (!r.ok || !P.cumulenes.length) return r;
      let X = xyzOf(m);
      const bad = P.cumulenes.filter(q => !cumuleneOK(X, q));
      if (!bad.length) return r;
      const was = bad.map(q => endPlanes(X, q)), turnable = bad.every(q => q.side);
      if (turnable) {
        for (const q of bad) {
          const d = dihedral(pre, q.sa, q.a, q.b, q.sb);
          turnTo(m, q, q.perpendicular ? (d >= 0 ? 90 : -90) : (Math.abs(d) <= 90 ? 0 : 180));
        }
        // minimised again with the saddle escape too, so that a turned-back model sitting on a stationary point by
        // symmetry is not taken for a minimum
        const turned = xyzOf(m), r2 = minimise(m, rand);
        X = xyzOf(m);
        if (r2.ok && P.cumulenes.every(q => cumuleneOK(X, q))) return Object.assign({}, r2, { tries: r.tries });
        setXYZ(m, turned);
      } else setXYZ(m, pre);
      const what = bad.map((q, i) => `the ${q.name} ${q.label} (${q.perpendicular ? 'perpendicular' : 'planar'} ends; minimisation gave ${Math.round(was[i])}° between the end planes)`).join(' and ');
      const why = `MMFF94 has no torsion term across an sp carbon, so it cannot hold ${what}; ` +
        (turnable ? 'the ends were set back by hand, so this model is not a force-field minimum' : "a ring joins the ends, so the conformer generator's geometry is shown");
      return { ok: false, rc: 0, energy: NaN, tries: r.tries, why, stop: 'approximate' };
    }

    // The model's own stereo parities, perceived from its coordinates (same atom indices as the model).
    function readBack(model) {
      const back = keepValence(model, Mol.fromMolfile(molfile3D(model)));
      back.ensureHelperArrays(Mol.cHelperParities);
      return back;
    }
    // Stereo elements the reference specifies (parity 1/2) must come out of the model's coordinates unchanged.
    // Same OCL index space in ref and model (the model only appends hydrogens).
    function stereoCheck(ref, back) {
      ref.ensureHelperArrays(Mol.cHelperParities);
      const res = { checked: 0, atoms: [], bonds: [] };
      for (let i = 0; i < ref.getAtoms(); i++) {
        const p = ref.getAtomParity(i);
        if (p === 1 || p === 2) { res.checked++; if (back.getAtomParity(i) !== p) res.atoms.push(i); }
      }
      for (let b = 0; b < ref.getBonds(); b++) {
        const p = ref.getBondParity(b);
        if (p !== 1 && p !== 2) continue;
        res.checked++;
        const a0 = ref.getBondAtom(0, b), a1 = ref.getBondAtom(1, b), bb = back.getBond(a0, a1);
        if (bb < 0 || back.getBondParity(bb) !== p) res.bonds.push([a0, a1]);
      }
      return res;
    }

    // Where the input leaves stereo open, the first model picks a configuration; the conformer list keeps that one
    // (otherwise "conformers" of CC=CC would include both E and Z, i.e. two different compounds).
    function pinOf(open, back) {
      return { atoms: open.atoms.map(j => [j, back.getAtomParity(j)]),
               bonds: open.bonds.map(([a0, a1]) => { const bb = back.getBond(a0, a1); return [a0, a1, bb < 0 ? -1 : back.getBondParity(bb)]; }) };
    }
    function samePin(pin, back) {
      for (const [j, p] of pin.atoms) if (back.getAtomParity(j) !== p) return false;
      for (const [a0, a1, p] of pin.bonds) { const bb = back.getBond(a0, a1); if ((bb < 0 ? -1 : back.getBondParity(bb)) !== p) return false; }
      return true;
    }

    // OpenChemLib's canonical code of a model with the parities of `flips` inverted (via the same frozen-parity 2D
    // drawing as prepare(); every stereo element absolute, so enantiomers get different codes).
    function stereoId(model, flips) {
      const c = model.getCompactCopy();
      c.ensureHelperArrays(Mol.cHelperParities);
      for (const j of flips) { const p = c.getAtomParity(j); if (p === 1 || p === 2) c.setAtomParity(j, 3 - p, c.isAtomParityPseudo(j)); }
      c.inventCoordinates(); c.setStereoBondsFromParity(); c.setParitiesValid(0);
      const r = keepValence(c, Mol.fromMolfile(c.toMolfile()));
      for (let j = 0; j < r.getAllAtoms(); j++) if (r.getAtomESRType(j) !== Mol.cESRTypeAbs) r.setAtomESR(j, Mol.cESRTypeAbs, 0);
      for (let b = 0; b < r.getAllBonds(); b++) if (r.getBondESRType(b) !== Mol.cESRTypeAbs) r.setBondESR(b, Mol.cESRTypeAbs, 0);
      return r.getIDCode();
    }
    // Are a and b the two bridgeheads of a small bridged ring system? Three separate a–b paths (bridges), each with at
    // least one atom, and both smallest rings through them ≤ 7 members (norbornane, [2.2.2], [3.2.1], [3.3.1]): then only
    // one relative configuration fits, and in/out isomers are impossible. (Greedy shortest paths: a miss only keeps a warning.)
    function smallBridgeheads(m, a, b) {
      if (m.getBond(a, b) >= 0) return false;
      const used = new Set(), lens = [];
      for (let t = 0; t < 3; t++) {
        const prev = new Map([[a, -1]]), queue = [a];
        let found = false;
        for (let h = 0; h < queue.length && !found; h++) {
          const x = queue[h];
          for (let k = 0; k < m.getConnAtoms(x); k++) {
            const y = m.getConnAtom(x, k);
            if (prev.has(y) || used.has(y)) continue;
            prev.set(y, x);
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
      return lens[0] + lens[2] + 2 <= 7;
    }
    /* Is the compound with the configuration at `flips` inverted the same compound as the model `back`? OpenChemLib's
       ID code says so for 7-methylnorbornane's C7 (true: a rotation turns one into the other), but also for C2 of
       2,6-dimethylspiro[3.3]heptane (false: that gives the enantiomer, an axially chiral compound; the ID code cannot
       express chirality through a spiro atom, and RDKit's canonical SMILES cannot either). So realOpen() takes the ID
       code's "same" only when this independent test agrees. From the model's own coordinates it reads the configuration
       of every tetrahedral centre (4 neighbours with at most one H, including spiro atoms; 3-coordinate P, S, As, Se, Sb,
       Te; bridgehead N with three ring bonds) and the cis/trans of every double bond outside rings of < 8 atoms, then
       searches for a symmetry of the graph (elements, charges, isotopes, H counts, bond orders, aromaticity) that maps
       the model's configuration onto the inverted one. A 3-coordinate N that is not a bridgehead inverts by itself, so
       its "configuration" is ignored. Returns false when there is no such symmetry, a flipped atom is not a
       configuration this test reads (an allene axis), or the search passes 200 000 steps: the warning then stays. */
    function sameWhenInverted(back, flips) {
      const n = back.getAllAtoms();
      back.ensureHelperArrays(Mol.cHelperRings);
      const isH = j => back.getAtomicNo(j) === 1 && !back.getAtomMass(j) && back.getAllConnAtoms(j) === 1;
      const heavy = [], hi = new Int32Array(n).fill(-1);
      for (let j = 0; j < n; j++) if (!isH(j)) { hi[j] = heavy.length; heavy.push(j); }
      const H = heavy.length, X = xyzOf(back);
      const nbr = heavy.map(() => []), hOn = new Int32Array(H), hAt = new Int32Array(H).fill(-1);
      for (let b = 0; b < back.getAllBonds(); b++) {
        const a0 = back.getBondAtom(0, b), a1 = back.getBondAtom(1, b);
        if (hi[a0] < 0 || hi[a1] < 0) { const p = hi[a0] < 0 ? hi[a1] : hi[a0]; if (p >= 0) { hOn[p]++; hAt[p] = hi[a0] < 0 ? a0 : a1; } continue; }
        const o = back.isAromaticBond(b) ? 4 : back.getBondOrder(b);
        nbr[hi[a0]].push([hi[a1], o]); nbr[hi[a1]].push([hi[a0], o]);
      }
      const adj = heavy.map(() => new Map());
      nbr.forEach((l, v) => { l.sort((x, y) => x[0] - y[0]); for (const [w, o] of l) adj[v].set(w, o); });
      // colours: element, charge, isotope, H count, then refined by neighbour colours (automorphisms keep them)
      let col = heavy.map((j, v) => [back.getAtomicNo(j), back.getAtomCharge(j), back.getAtomMass(j), back.getAtomRadical(j), hOn[v]].join(','));
      for (let it = 0; it < H; it++) {
        const next = col.map((c, v) => c + '|' + nbr[v].map(([w, o]) => o + ':' + col[w]).sort().join(';'));
        const ids = new Map(); next.forEach(c => { if (!ids.has(c)) ids.set(c, ids.size); });
        const renamed = next.map(c => String(ids.get(c))), before = new Set(col).size;
        col = renamed;
        if (ids.size === before) break;
      }
      // configurations read from the coordinates
      const P3 = v => { const j = heavy[v]; return [X[3 * j], X[3 * j + 1], X[3 * j + 2]]; };
      const sub = (p, q) => [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
      const det = (a, b, c) => a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);
      const PYR = new Set([15, 16, 33, 34, 51, 52]);
      const tet = new Int8Array(H);                       // ±1 per tetrahedral centre (0: none), neighbours in index order
      for (let v = 0; v < H; v++) {
        const j = heavy[v], deg = nbr[v].length + hOn[v];
        if (hOn[v] > 1) continue;
        const z = back.getAtomicNo(j), ringBonds = nbr[v].filter(([w]) => back.getBond(j, heavy[w]) >= 0 && back.isRingBond(back.getBond(j, heavy[w]))).length;
        if (!(deg === 4 || (deg === 3 && !back.isAromaticAtom(j) && (PYR.has(z) || (z === 7 && ringBonds === 3))))) continue;
        const pts = nbr[v].map(([w]) => P3(w));
        if (hOn[v]) pts.push([X[3 * hAt[v]], X[3 * hAt[v] + 1], X[3 * hAt[v] + 2]]);   // the one H goes last (it maps to the H)
        const d = pts.length === 4 ? det(sub(pts[1], pts[0]), sub(pts[2], pts[0]), sub(pts[3], pts[0])) : det(sub(pts[0], P3(v)), sub(pts[1], P3(v)), sub(pts[2], P3(v)));
        if (Math.abs(d) > 1e-3) tet[v] = d > 0 ? 1 : -1;
      }
      const dbl = new Map();                              // 'v,w' (v < w) → ±1 cis/trans of the lowest-index neighbours
      const refOf = (v, w) => { let r = -1; for (const [x] of nbr[v]) if (x !== w && (r < 0 || x < r)) r = x; return r; };
      for (let b = 0; b < back.getAllBonds(); b++) {
        const a0 = back.getBondAtom(0, b), a1 = back.getBondAtom(1, b), v = hi[a0], w = hi[a1];
        if (v < 0 || w < 0 || back.isAromaticBond(b) || back.getBondOrder(b) !== 2 || (back.isRingBond(b) && back.getBondRingSize(b) < 8)) continue;
        if (nbr[v].length + hOn[v] !== 3 || nbr[w].length + hOn[w] !== 3 || hOn[v] > 1 || hOn[w] > 1) continue;
        const rv = refOf(v, w), rw = refOf(w, v);
        if (rv < 0 || rw < 0) continue;
        const dh = dihedral(X, heavy[rv], a0, a1, heavy[rw]);
        dbl.set(Math.min(v, w) + ',' + Math.max(v, w), Math.abs(dh) < 90 ? 1 : -1);
      }
      // the target: the model with `flips` inverted
      const want = Int8Array.from(tet);
      for (const j of flips) {
        const v = hi[j];
        if (v < 0) return false;
        if (tet[v]) want[v] = -tet[v];
        else if (!(back.getAtomicNo(j) === 7 && nbr[v].length + hOn[v] === 3)) return false;   // a configuration not read here
      }
      // search: map heavy atoms in breadth-first order; check each configuration as soon as all its atoms are mapped
      const order = [], seen = new Uint8Array(H);
      for (let s = 0; s < H; s++) {
        if (seen[s]) continue;
        seen[s] = 1; order.push(s);
        for (let q = order.length - 1; q < order.length; q++) for (const [w] of nbr[order[q]]) if (!seen[w]) { seen[w] = 1; order.push(w); }
      }
      const pos = new Int32Array(H);
      order.forEach((v, k) => { pos[v] = k; });
      const checksAt = order.map(() => []);             // configurations completed when order[k] is mapped
      for (let v = 0; v < H; v++) if (tet[v]) checksAt[Math.max(pos[v], ...nbr[v].map(([w]) => pos[w]))].push(['t', v]);
      for (const key of dbl.keys()) {
        const [v, w] = key.split(',').map(Number);
        checksAt[Math.max(pos[v], pos[w], ...nbr[v].map(([x]) => pos[x]), ...nbr[w].map(([x]) => pos[x]))].push(['d', v, w]);
      }
      const sigma = new Int32Array(H).fill(-1), used = new Uint8Array(H);
      const parity = list => { let p = 1; for (let i = 0; i < list.length; i++) for (let k = i + 1; k < list.length; k++) if (list[i] > list[k]) p = -p; return p; };
      const okCheck = c => {
        if (c[0] === 't') {
          const v = c[1], u = sigma[v];
          if (!tet[u]) return false;
          return tet[v] * parity(nbr[v].map(([w]) => sigma[w])) === want[u];
        }
        const [, v, w] = c, u = sigma[v], x = sigma[w], key = Math.min(u, x) + ',' + Math.max(u, x);
        if (!dbl.has(key)) return false;
        const f = (a, b, ia, ib) => (sigma[refOf(a, b)] === refOf(ia, ib) ? 1 : -1);
        return dbl.get(Math.min(v, w) + ',' + Math.max(v, w)) * f(v, w, u, x) * f(w, v, x, u) === dbl.get(key);
      };
      let steps = 0;
      const go = k => {
        if (k === H) return true;
        if (++steps > 200000) throw new Error('symmetry search too large');
        const v = order[k];
        for (let u = 0; u < H; u++) {
          if (used[u] || col[u] !== col[v]) continue;
          let fit = true, mappedV = 0, mappedU = 0;
          for (const [w, o] of nbr[v]) if (sigma[w] >= 0) { mappedV++; if (adj[u].get(sigma[w]) !== o) { fit = false; break; } }
          if (!fit) continue;
          for (const [x] of nbr[u]) if (used[x]) mappedU++;
          if (mappedU !== mappedV) continue;
          sigma[v] = u; used[u] = 1;
          if (checksAt[k].every(okCheck) && go(k + 1)) return true;
          sigma[v] = -1; used[u] = 0;
        }
        return false;
      };
      try { return go(0); } catch (e) { return false; }
    }

    /* The stereo the input leaves open, as far as it is a real choice. OpenChemLib also lists atoms whose configuration
       is not one: the bridgeheads of a small bridged ring (norbornane: inverting both gives the same compound, inverting
       one is impossible), and atoms that give the same compound either way (7-methylnorbornane's C7). Those are left out
       of the "not specified" warning; the model's configuration there is the only one. Everything else stays, including
       ring cis/trans that RDKit's CIP labels miss (1,4-dimethylcyclohexane), camphor's bridgeheads (two enantiomers) and
       the ring atoms of 2,6-disubstituted spiro[3.3]heptanes (axially chiral; see sameWhenInverted()).
       `back`: a model read back from its coordinates (model indices; the file's atoms keep theirs). */
    function realOpen(P, back) {
      const open = P.open.atoms, keepAll = P.unspecified.slice();
      if (!open.length || open.length > 16) return keepAll;
      const drop = new Set();
      try {
        const base = stereoId(back, []);
        for (let x = 0; x < open.length; x++) for (let y = x + 1; y < open.length; y++) {
          if (smallBridgeheads(P.m0, open[x], open[y]) && stereoId(back, [open[x], open[y]]) === base && sameWhenInverted(back, [open[x], open[y]])) { drop.add(open[x]); drop.add(open[y]); }
        }
        for (const j of open) if (!drop.has(j) && stereoId(back, [j]) === base && sameWhenInverted(back, [j])) drop.add(j);
      } catch (e) { return keepAll; }
      const out = new Set([...drop].map(j => P.inv0[j]));
      return keepAll.filter(u => u.atom == null || !out.has(u.atom));
    }

    // Horn's quaternion method: RMSD after optimal superposition. A and B are flat [x,y,z,...] of the same atoms.
    function rmsd(A, B) {
      const n = A.length / 3;
      if (n < 1) return 0;
      const ca = [0, 0, 0], cb = [0, 0, 0];
      for (let i = 0; i < n; i++) for (let k = 0; k < 3; k++) { ca[k] += A[3 * i + k] / n; cb[k] += B[3 * i + k] / n; }
      let ga = 0, gb = 0;
      const S = [0, 0, 0, 0, 0, 0, 0, 0, 0];
      for (let i = 0; i < n; i++) {
        const a = [A[3 * i] - ca[0], A[3 * i + 1] - ca[1], A[3 * i + 2] - ca[2]], b = [B[3 * i] - cb[0], B[3 * i + 1] - cb[1], B[3 * i + 2] - cb[2]];
        ga += a[0] * a[0] + a[1] * a[1] + a[2] * a[2]; gb += b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
        for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) S[3 * r + c] += a[r] * b[c];
      }
      const [xx, xy, xz, yx, yy, yz, zx, zy, zz] = S;
      const K = [[xx + yy + zz, yz - zy, zx - xz, xy - yx], [yz - zy, xx - yy - zz, xy + yx, zx + xz],
                 [zx - xz, xy + yx, -xx + yy - zz, yz + zy], [xy - yx, zx + xz, yz + zy, -xx - yy + zz]];
      const shift = (ga + gb) / 2 + 1e-9;      // power iteration on K + shift·I finds K's largest eigenvalue
      let v = [1, 0.1, 0.01, 0.001], lam = 0;
      for (let it = 0; it < 1000; it++) {
        const w = [0, 1, 2, 3].map(r => K[r][0] * v[0] + K[r][1] * v[1] + K[r][2] * v[2] + K[r][3] * v[3] + shift * v[r]);
        const norm = Math.hypot(w[0], w[1], w[2], w[3]), nl = norm - shift;
        v = w.map(x => x / norm);
        if (Math.abs(nl - lam) < 1e-12 * Math.max(1, Math.abs(nl))) { lam = nl; break; }
        lam = nl;
      }
      return Math.sqrt(Math.max(0, (ga + gb - 2 * lam) / n));
    }

    function failure(e) {
      return { ok: false, error: '3D engine error: ' + errText(e), source: 'OpenChemLib', conformers: [], warnings: [], worker: false };
    }

    /* One build = parse + first conformer, then (step by step) the conformer search. The driver (Worker or main
       thread) calls step() until it returns false, posting snapshot(true) whenever takeDirty() says the list changed. */
    function newBuild(text, opts, limits) {
      opts = opts || {}; limits = limits || {};
      const t0 = now();
      const prev = opts.previous && opts.previous.ok && opts.previous.conformers && opts.previous.conformers.length ? opts.previous : null;
      const k = Math.max(0, Math.floor(opts.k == null ? 10 : +opts.k) || 0);
      const budget = opts.timeBudgetMs == null ? 20000 : Math.max(0, +opts.timeBudgetMs || 0);
      const seed = (opts.seed == null ? (prev ? (prev.seed >>> 0) + 1 : 12345) : +opts.seed) >>> 0;
      const rand = rng(seed ^ 0x9E3779B9);
      const S = { phase: 'first', P: null, error: null, confs: [], nextId: 0, tried: 0, stereoDropped: 0, notConverged: 0, approximate: 0, cg: null,
                  target: k, potential: null, stop: null, lastStepMs: 0, dirty: false, msFirst: null, stereo: null, pin: null, otherConfig: 0, notes: [] };

      function fail(msg) { S.error = msg; S.phase = 'done'; S.dirty = true; }

      function conformerFrom(m, id, raw, r, ms) {
        const P = S.P, coords = new Float64Array(3 * P.N);
        for (let o = 0; o < P.N; o++) {
          const j = P.out2model[o];
          coords[3 * o] = m.getAtomX(j); coords[3 * o + 1] = -m.getAtomY(j); coords[3 * o + 2] = -m.getAtomZ(j);
        }
        const flaggedBonds = [], warnings = [];
        if (!r.ok) {
          let why = r.why || `MMFF94 did not converge (code ${r.rc})`;
          if (!r.why && r.error) {
            const x = /atom (\d+) \(([A-Za-z]+)\)/.exec(r.error), o = x ? P.model2out[+x[1]] : -1;
            why = x && o >= 0 && P.elements[o] === x[2] ? `MMFF94 could not assign an atom type to ${P.lab(o)}` : 'MMFF94: ' + r.error;
          }
          warnings.push('approximate geometry (not force-field optimised): ' + why);
        } else if (raw) {
          for (let b = 0; b < m.getAllBonds(); b++) {
            const i = m.getBondAtom(0, b), j = m.getBondAtom(1, b);
            if (m.getAtomicNo(i) === 1 || m.getAtomicNo(j) === 1) continue;
            const d = dist(m, i, j) - raw[b];
            if (Math.abs(d) > GAP) {
              const oi = P.model2out[i], oj = P.model2out[j];
              flaggedBonds.push([Math.min(oi, oj), Math.max(oi, oj), `MMFF94 has no specific parameters for this bond: minimisation made it ${Math.abs(d).toFixed(2)} Å ${d > 0 ? 'longer' : 'shorter'} than the crystal-structure length, so its length is unreliable`]);
            }
          }
          flaggedBonds.sort((x, y) => x[0] - y[0] || x[1] - y[1]);
        }
        const heavyXYZ = new Float64Array(3 * P.heavy.length);
        P.heavy.forEach((o, q) => { heavyXYZ[3 * q] = coords[3 * o]; heavyXYZ[3 * q + 1] = coords[3 * o + 1]; heavyXYZ[3 * q + 2] = coords[3 * o + 2]; });
        const c = { id, coords, energy: r.ok ? r.energy : NaN, relEnergy: NaN, optimised: r.ok, rc: r.rc, tries: r.tries || 0, flaggedBonds, warnings, ms: Math.round(ms), heavyXYZ };
        if (P.N === 1) c.note = SINGLE_NOTE;
        return c;
      }

      // stereo of model m against the frozen reference and against the file's own parities
      function checkStereo(m, back) {
        const P = S.P, a = stereoCheck(P.ref, back), b = stereoCheck(P.m0, back);
        const atoms = [...new Set(a.atoms.concat(b.atoms))].map(j => P.model2out[j]).sort((x, y) => x - y);
        const seen = new Set(), bonds = [];
        for (const [i, j] of a.bonds.concat(b.bonds)) {
          const oi = P.model2out[i], oj = P.model2out[j], key = Math.min(oi, oj) + ',' + Math.max(oi, oj);
          if (!seen.has(key)) { seen.add(key); bonds.push([Math.min(oi, oj), Math.max(oi, oj)]); }
        }
        return { checked: Math.max(a.checked, b.checked), wrong: atoms.map(o => ({ atom: o })).concat(bonds.map(x => ({ bond: x }))) };
      }
      // Names of stereo elements: C3 (a stereocentre), C2=C3 (E/Z), "the allene axis C2=C3=C4" (OCL puts an allene's
      // parity on its central atom), "the C45–C50 axis" (a single bond OCL treats as stereo: an atropisomeric biaryl axis).
      function stereoName(w) {
        const P = S.P;
        if (w.atom != null) {
          const j = P.out2model[w.atom], q = P.cumulenes.find(c => c.chain.indexOf(j) > 0 && c.chain.indexOf(j) < c.chain.length - 1);
          return q ? `the ${q.name} axis ${q.label}` : P.lab(w.atom);
        }
        const [a, b] = w.bond, bd = P.bonds.find(x => (x[0] === a && x[1] === b) || (x[0] === b && x[1] === a));
        return bd && bd[2] === 1 ? `the ${P.lab(a)}–${P.lab(b)} axis` : P.lab(a) + '=' + P.lab(b);
      }
      const isAxis = w => w.bond != null && /axis$/.test(stereoName(w));
      const stereoText = list => list.map(stereoName).join(', ');
      // "stereochemistry at C3, C5=C6 is …" / "the axial (atropisomeric) configuration about C45–C50 is …"
      function openStereoWarning(list) {
        const plain = list.filter(w => !isAxis(w)), axes = list.filter(isAxis);
        const parts = [];
        if (plain.length) parts.push('stereochemistry at ' + stereoText(plain));
        if (axes.length) parts.push(`the axial (atropisomeric) configuration${axes.length > 1 ? 's' : ''} about ` + axes.map(w => S.P.lab(w.bond[0]) + '–' + S.P.lab(w.bond[1])).join(', '));
        return `${parts.join(' and ')} ${parts.length > 1 || axes.length > 1 ? 'are' : 'is'} not specified in the input: the model shows the configuration the builder chose`;
      }

      function first() {
        const P = prepare(text);
        if (P.error) return fail(P.error);
        S.P = P;
        if (limits.maxHeavyAtoms && P.heavy.length > limits.maxHeavyAtoms)
          return fail(`${P.heavy.length} heavy atoms is too large to build without a background worker (limit ${limits.maxHeavyAtoms})`);
        if (prev) return restore();
        const ts = now();
        const m = P.withH.getCompactCopy();
        let raw = null;
        if (P.fromFile) {                                  // the file's own coordinates (every atom is in the file)
          for (let i = 0; i < P.n0; i++) { const j = P.out2model[i], j0 = P.map0[i]; m.setAtomX(j, P.m0.getAtomX(j0)); m.setAtomY(j, P.m0.getAtomY(j0)); m.setAtomZ(j, P.m0.getAtomZ(j0)); }
          S.notes.push("first conformer: the file's own 3D coordinates, then MMFF94");
        } else {
          let c;
          try { c = new CG(seed).getOneConformerAsMolecule(m); } catch (e) { return fail('the conformer generator failed: ' + errText(e)); }
          if (!c) return fail('OpenChemLib could not generate a 3D structure for this molecule (unusual valence or structure)');
          if (m.getAllAtoms() !== P.N) return fail(MAPFAIL + ' (the conformer generator changed the atom list)');
          raw = bondLengths(m);
        }
        // a hypervalent atom: MMFF94 would build the wrong shape, so it is not run at all (like an atom it cannot type).
        // A lone centre with identical terminal ligands (I3−, SF4, SF6) gets its ideal VSEPR shape; otherwise the
        // generator's (or the file's) geometry is shown. Either way the model is marked approximate.
        const shape = P.hyper.length === 1 && !P.fromFile ? idealiseVSEPR(m, P.hyperInfo[0].model, P.hyperInfo[0].sigma, P.hyperInfo[0].lp) : null;
        const pre = xyzOf(m);
        const r = P.hyper.length ? { ok: false, rc: null, energy: NaN, tries: 0, stop: 'approximate',
                                     why: `MMFF94 has no parameters for ${hyperText(P.hyperInfo.map(h => ({ name: P.lab(h.out), sigma: h.sigma, lp: h.lp })))}` +
                                          (shape ? `; the model shows the ideal VSEPR shape (${shape}) with the conformer generator's bond lengths` : ' (VSEPR shape not modelled)') }
                                 : optimise(m, rand, P, pre);
        const back = readBack(m), st = checkStereo(m, back);
        if (st.wrong.length) {
          S.stereo = { checked: st.checked, wrong: st.wrong, unspecified: P.unspecified.slice() };
          return fail(`OpenChemLib built the wrong stereoisomer (at ${stereoText(st.wrong)}): no 3D model is shown, because it would be a different compound`);
        }
        S.stereo = { checked: st.checked, wrong: [], unspecified: realOpen(P, back) };
        S.pin = pinOf(P.open, back);
        S.cumPin = cumuleneConfig(P, xyzOf(m));
        const conf = conformerFrom(m, S.nextId++, raw, r, now() - ts);
        S.confs.push(conf);
        S.msFirst = Math.round(now() - t0);
        S.dirty = true;
        if (r.single) {                                   // [Na+]: MMFF94 types it; there is nothing to optimise or search
          S.notes.push(SINGLE_NOTE + ' (MMFF94 has a type for it; its energy is 0 by definition), so no conformer search');
          S.stop = 'single-atom'; S.phase = 'done';
          return;
        }
        if (!r.ok && (r.rc === null || r.stop)) {          // no energies: nothing to rank, and every conformer would be the same case
          S.stop = r.stop || 'untypable'; S.phase = 'done';
          if (r.stop) S.notes.push('no conformer search: MMFF94 cannot optimise this structure (see the conformer warning)');
          return;
        }
        S.phase = k > 0 ? 'init' : 'done';
        if (!k) S.stop = 'k';
      }

      // "Search more": keep the earlier result's conformers (same topology checked) and search with another seed
      function restore() {
        const P = S.P;
        const same = prev.N === P.N && prev.elements && prev.elements.every((e, i) => e === P.elements[i]) &&
          prev.parents && Array.from(prev.parents).every((p, i) => p === P.parents[i]) &&
          prev.bonds && prev.bonds.length === P.bonds.length && prev.bonds.every((b, i) => b[0] === P.bonds[i][0] && b[1] === P.bonds[i][1]);
        if (!same) return fail(MAPFAIL + ' (the earlier result belongs to a different structure)');
        for (const c of prev.conformers) {
          const coords = Float64Array.from(c.coords), heavyXYZ = new Float64Array(3 * P.heavy.length);
          P.heavy.forEach((o, q) => { heavyXYZ[3 * q] = coords[3 * o]; heavyXYZ[3 * q + 1] = coords[3 * o + 1]; heavyXYZ[3 * q + 2] = coords[3 * o + 2]; });
          S.confs.push({ id: c.id, coords, energy: c.energy, relEnergy: NaN, optimised: c.optimised, rc: c.rc, tries: c.tries || 0,
                         flaggedBonds: (c.flaggedBonds || []).map(x => x.slice()), warnings: (c.warnings || []).slice(), ms: c.ms, heavyXYZ });
          if (P.N === 1) S.confs[S.confs.length - 1].note = SINGLE_NOTE;
          S.nextId = Math.max(S.nextId, c.id + 1);
        }
        const m = P.withH.getCompactCopy(), c0 = S.confs[0].coords;       // the earlier model's configuration, read back
        for (let o = 0; o < P.N; o++) { const j = P.out2model[o]; m.setAtomX(j, c0[3 * o]); m.setAtomY(j, -c0[3 * o + 1]); m.setAtomZ(j, -c0[3 * o + 2]); }
        const back = readBack(m);
        S.stereo = { checked: prev.stereo ? prev.stereo.checked : 0, wrong: [], unspecified: realOpen(P, back) };
        S.pin = pinOf(P.open, back);
        S.cumPin = cumuleneConfig(P, xyzOf(m));
        S.msFirst = 0; S.dirty = true;
        if (P.N === 1 && S.confs[0].optimised) { S.stop = 'single-atom'; S.phase = 'done'; return; }
        if (!S.confs.some(c => c.optimised)) { S.stop = (prev.search && prev.search.stoppedBy === 'approximate') ? 'approximate' : 'untypable'; S.phase = 'done'; return; }
        S.phase = k > 0 ? 'init' : 'done';
        if (!k) S.stop = 'k';
      }

      function init() {
        // setting up the search costs up to ~0.65 × the first conformer's time (cholesterol 4.0 s → 2.6 s, erythromycin
        // 5.1 s → 3.0 s): do not start it when it would clearly overrun the budget
        if (now() - t0 + 0.75 * (S.msFirst || 0) > budget) { S.stop = 'time'; S.phase = 'done'; return; }
        const base = S.P.withH.getCompactCopy();
        const cg = new CG(seed);
        let ok;
        try { ok = cg.initializeConformers(base, { strategy: CG.STRATEGY_LIKELY_RANDOM, maxTorsionSets: 100000, use60degreeSteps: false }); } catch (e) { ok = false; }
        if (!ok) { S.stop = 'no-search'; S.notes.push('the conformer search could not start for this structure'); S.phase = 'done'; return; }
        S.cg = cg;
        S.potential = cg.getPotentialConformerCount();
        if (S.potential > 0 && S.potential < S.target) S.target = S.potential;   // e.g. benzene: 1 torsion set
        S.phase = 'search';
      }

      function search() {
        const P = S.P;
        if (S.tried >= S.target) { S.stop = S.target < k ? 'exhausted' : 'k'; S.phase = 'done'; return; }
        // do not start a conformer that would probably overrun the budget
        if (now() - t0 + S.lastStepMs > budget) { S.stop = 'time'; S.phase = 'done'; return; }
        const ts = now();
        const c = S.cg.getNextConformerAsMolecule(null);   // a fresh copy each time, same atom order; null = nothing new
        if (!c) { S.stop = 'exhausted'; S.phase = 'done'; return; }
        S.tried++;
        if (c.getAllAtoms() !== P.N) { S.notes.push(MAPFAIL + ' for a search conformer: skipped'); S.lastStepMs = now() - ts; return; }
        const raw = bondLengths(c), r = optimise(c, rand, P, xyzOf(c));
        if (!r.ok) { if (r.stop) S.approximate++; else S.notConverged++; S.lastStepMs = now() - ts; return; }
        const back = readBack(c);
        if (checkStereo(c, back).wrong.length) S.stereoDropped++;
        else if (!samePin(S.pin, back) || cumuleneConfig(P, xyzOf(c)) !== S.cumPin) S.otherConfig++;
        else { S.confs.push(conformerFrom(c, S.nextId++, raw, r, now() - ts)); S.dirty = true; }
        S.lastStepMs = now() - ts;
      }

      function step() {
        if (S.phase === 'first') first();
        else if (S.phase === 'init') init();
        else if (S.phase === 'search') search();
        return S.phase !== 'done';
      }

      // ranked, duplicates removed; the lowest-energy member of each duplicate set is kept
      function ranked() {
        const opt = S.confs.filter(c => c.optimised && isFinite(c.energy)).sort((a, b) => a.energy - b.energy || a.id - b.id);
        const kept = [], merged = [];
        for (const f of opt) {
          const dup = kept.find(q => Math.abs(q.energy - f.energy) <= E_TOL || rmsd(f.heavyXYZ, q.heavyXYZ) < RMSD_CUT);
          if (dup) merged.push([f.id, dup.id]); else kept.push(f);
        }
        if (!kept.length) return { kept: S.confs.slice(0, 1), merged };   // nothing optimised: the first (approximate) model only
        return { kept, merged };
      }

      function snapshot(partial) {
        const P = S.P;
        const res = { ok: !S.error, source: 'OpenChemLib', table: TABLE, seed, partial: !!partial, worker: false, ms: Math.round(now() - t0), msFirst: S.msFirst };
        if (P && !P.error) {
          Object.assign(res, { N: P.N, n: P.n0, singleAtom: P.N === 1, elements: P.elements.slice(), parents: Int32Array.from(P.parents), bonds: P.bonds.map(b => b.slice()),
                               charges: Int32Array.from(P.charges), isotopes: Int32Array.from(P.isotopes), radicals: Int32Array.from(P.radicals) });
        }
        res.stereo = S.stereo ? { checked: S.stereo.checked, wrong: S.stereo.wrong.map(x => Object.assign({}, x)), unspecified: S.stereo.unspecified.map(x => Object.assign({}, x)) } : null;
        const warnings = [];
        if (P && !P.error) {
          // OCL's list less the atoms that are not a real choice (realOpen: norbornane's bridgeheads); it keeps ring
          // stereo that RDKit's CIP labels miss (cis/trans-1,4-dimethylcyclohexane)
          if (S.stereo && S.stereo.unspecified.length) warnings.push(openStereoWarning(S.stereo.unspecified));
          for (const w of P.cumuleneWarnings) warnings.push(w);
          for (const w of P.flagWarnings) warnings.push(w);
          res.atomFlags = P.atomFlags.map(f => Object.assign({}, f));
        }
        res.warnings = warnings;
        // counts cover the whole list: after "search more" (prev), the earlier search's counts are carried forward
        const ps = prev && prev.search ? prev.search : {}, add = key => (ps[key] | 0);
        res.generated = (prev ? (prev.generated | 0) : (S.confs.length ? 1 : 0)) + S.tried;
        res.potential = S.potential;
        res.search = { tried: add('tried') + S.tried, target: add('target') + S.target, stoppedBy: partial ? null : S.stop, stereoDropped: add('stereoDropped') + S.stereoDropped,
                       otherConfiguration: add('otherConfiguration') + S.otherConfig, notConverged: add('notConverged') + S.notConverged,
                       approximate: add('approximate') + S.approximate, notes: S.notes.slice() };
        if (S.error) { res.error = S.error; res.conformers = []; return res; }
        const { kept, merged } = ranked();
        const e0 = kept.length && kept[0].optimised ? kept[0].energy : NaN;
        res.conformers = kept.map(c => Object.assign({ id: c.id, coords: Float64Array.from(c.coords), energy: c.energy, relEnergy: c.optimised ? c.energy - e0 : NaN,
                                          optimised: c.optimised, rc: c.rc, tries: c.tries, flaggedBonds: c.flaggedBonds.map(x => x.slice()), warnings: c.warnings.slice(), ms: c.ms },
                                          c.note ? { note: c.note } : {}));
        res.merged = merged;
        return res;
      }

      return { step, snapshot, takeDirty() { const d = S.dirty; S.dirty = false; return d; } };
    }

    // Single-point MMFF94 energy of a 3D structure that already has every hydrogen (e.g. PubChem's coordinates).
    function energy(text) {
      let parsed;
      try { parsed = Mol.fromMolfileWithAtomMap(String(text)); } catch (e) { return { ok: false, energy: NaN, error: 'could not read the structure: ' + errText(e) }; }
      const m = parsed.molecule;
      if (!m.getAllAtoms()) return { ok: false, energy: NaN, error: 'empty structure' };
      const header = String(text).split(/\r\n|\r|\n/)[1] || '';
      if (header.slice(20, 22) === '2D' || (!m.is3D() && header.slice(20, 22) !== '3D')) return { ok: false, energy: NaN, error: 'needs 3D coordinates' };
      m.ensureHelperArrays(Mol.cHelperNeighbours);
      let missing = 0;
      for (let i = 0; i < m.getAllAtoms(); i++) missing += m.getImplicitHydrogens(i);
      if (missing) return { ok: false, energy: NaN, error: `${missing} hydrogen atom${missing > 1 ? 's are' : ' is'} missing (MMFF94 needs every hydrogen)` };
      if (m.getFragmentNumbers([], false, false) > 1) return { ok: false, energy: NaN, error: 'several separate components' };
      const hv = hypervalent(m), fileNo = j => Array.from(parsed.map || []).indexOf(j) + 1;   // m's index → 1-based file number
      if (hv.length) return { ok: false, energy: NaN, error: `MMFF94 has no parameters for ${hyperText([Object.assign({ name: m.getAtomLabel(hv[0].j) + (fileNo(hv[0].j) || '') }, hv[0])])} (VSEPR shape not modelled)` };
      try {
        if (m.getAllAtoms() < 2) { new FF(m, TABLE); return { ok: true, energy: 0, table: TABLE, singleAtom: true, note: SINGLE_NOTE }; }   // typing test only (He, U, H+ fail)
        const e = new FF(m, TABLE).getTotalEnergy();
        return isFinite(e) ? { ok: true, energy: e, table: TABLE } : { ok: false, energy: NaN, error: 'MMFF94 gave no finite energy' };
      } catch (e) { return { ok: false, energy: NaN, error: 'MMFF94: ' + errText(e) }; }
    }

    return { register, newBuild, energy, failure, rmsd, version: OCL.version };
  }

  /* ================= worker side ================= */
  // Also rebuilt from source text inside the Worker; `core` is engineCore(OCL) made there.
  function workerMain(self, core) {
    'use strict';
    const msg = e => String((e && e.message) || e);
    self.onmessage = ev => {
      const d = ev.data || {};
      if (d.type === 'init') {
        try { core.register(d.resources); self.postMessage({ type: 'ready', version: core.version }); }
        catch (e) { self.postMessage({ type: 'initError', error: msg(e) }); }
      } else if (d.type === 'build') {
        let b;
        try { b = core.newBuild(d.molblock, d.opts, null); } catch (e) { self.postMessage({ type: 'done', id: d.id, result: core.failure(e) }); return; }
        const tick = () => {
          let more;
          try { more = b.step(); } catch (e) { self.postMessage({ type: 'done', id: d.id, result: core.failure(e) }); return; }
          if (b.takeDirty()) { const r = b.snapshot(true); if (r.ok) self.postMessage({ type: 'progress', id: d.id, result: r }); }
          if (more) setTimeout(tick, 0);
          else self.postMessage({ type: 'done', id: d.id, result: b.snapshot(false) });
        };
        tick();
      }
    };
  }

  /* ================= page side: loading, the Worker, the job queue ================= */
  const options = { base: 'vendor/openchemlib/', worker: true, mainThreadMaxHeavyAtoms: 120, workerStartTimeoutMs: 30000 };
  const S = { loading: null, mode: null, worker: null, workerURL: null, workerStarting: null, workerError: null, mainCore: null, queue: [], active: null, seq: 0 };
  const message = e => String((e && e.message) || e);

  function loadLibrary() {
    if (window.OCL_FACTORY && window.OCL && window.OCL_RESOURCES) return Promise.resolve();
    const L = window.RDKitLoad;
    if (!L || typeof L.script !== 'function') return Promise.reject(new Error('the script loader (js/rdkit_load.js) is not on this page'));
    return Promise.all([
      window.OCL_FACTORY ? null : L.script(options.base + 'openchemlib.js'),
      window.OCL_RESOURCES ? null : L.script(options.base + 'openchemlib-resources.js'),
    ]).then(() => { if (!window.OCL_FACTORY || !window.OCL || !window.OCL_RESOURCES) throw new Error('OpenChemLib did not load'); });
  }
  function mainCore() {
    if (!S.mainCore) { const c = engineCore(window.OCL); c.register(window.OCL_RESOURCES); S.mainCore = c; }
    return S.mainCore;
  }
  function workerSource() {
    return '/* Mol3DEngine worker: OpenChemLib and the engine core, rebuilt from their source text */\n' +
      'var OCL = (' + window.OCL_FACTORY.toString() + ')();\n' +
      'var CORE = (' + engineCore.toString() + ')(OCL);\n' +
      '(' + workerMain.toString() + ')(self, CORE);\n';
  }
  function startWorker() {
    return new Promise((resolve, reject) => {
      let w, timer = null;
      const fail = err => { clearTimeout(timer); try { w && w.terminate(); } catch (e) {} reject(err); };
      try {
        if (!S.workerURL) S.workerURL = URL.createObjectURL(new Blob([workerSource()], { type: 'text/javascript' }));
        w = new Worker(S.workerURL);
      } catch (e) { fail(e); return; }
      timer = setTimeout(() => fail(new Error(`the 3D worker did not start within ${options.workerStartTimeoutMs / 1000} s`)), options.workerStartTimeoutMs);
      w.onerror = ev => { if (ev && ev.preventDefault) ev.preventDefault(); fail(new Error('the 3D worker could not start: ' + ((ev && ev.message) || 'error'))); };
      w.onmessage = ev => {
        const d = ev.data || {};
        if (d.type === 'ready') { clearTimeout(timer); w.onmessage = onWorkerMessage; w.onerror = onWorkerError; resolve(w); }
        else if (d.type === 'initError') fail(new Error(d.error));
      };
      w.postMessage({ type: 'init', resources: window.OCL_RESOURCES });
    });
  }
  function ensureWorker() {
    if (S.worker) return Promise.resolve(S.worker);
    if (!S.workerStarting) {
      S.workerStarting = startWorker().then(w => { S.workerStarting = null; S.worker = w; return w; },
                                            e => { S.workerStarting = null; throw e; });
    }
    return S.workerStarting;
  }
  function useMainThread(reason) { S.workerError = S.workerError || reason; S.mode = 'main'; }

  // Loads OpenChemLib and starts the Worker (or settles for the main thread). Resolves to 'worker' or 'main'.
  function ready() {
    if (!S.loading) {
      S.loading = (async () => {
        await loadLibrary();
        const canWorker = options.worker && typeof Worker !== 'undefined' && typeof Blob !== 'undefined' && typeof URL !== 'undefined' && URL.createObjectURL;
        if (canWorker) {
          try { await ensureWorker(); S.mode = 'worker'; return S.mode; } catch (e) { useMainThread(message(e)); }
        } else useMainThread(options.worker ? 'Web Workers are not available here' : 'worker disabled by options');
        mainCore();
        return S.mode;
      })();
      S.loading.catch(() => { S.loading = null; });
    }
    return S.loading;
  }

  function errorResult(error, extra) {
    return Object.assign({ ok: false, error, source: 'OpenChemLib', conformers: [], warnings: [], worker: S.mode === 'worker' }, extra || {});
  }
  function finishJob(job, result) {
    if (job.done) return;
    job.done = true;
    if (S.active === job) S.active = null;
    job.resolve(result);
    pump();
  }
  function progress(job, r) {
    if (job.done) return;
    job.last = r;
    if (job.onProgress) { try { job.onProgress(r); } catch (e) { console.error(e); } }
  }

  async function pump() {
    if (S.active || !S.queue.length) return;
    const job = S.queue.shift();
    S.active = job;
    try { await ready(); } catch (e) { finishJob(job, errorResult('the 3D engine could not load: ' + message(e))); return; }
    if (job.done) return;                                         // cancelled while loading
    if (S.mode === 'worker') {
      try { await ensureWorker(); } catch (e) { useMainThread('the 3D worker could not restart: ' + message(e)); }
      if (job.done) return;
    }
    if (S.mode === 'worker') { job.worker = S.worker; S.worker.postMessage({ type: 'build', id: job.id, molblock: job.molblock, opts: job.opts }); }
    else runOnMainThread(job);
  }
  function onWorkerMessage(ev) {
    const d = ev.data || {}, job = S.active;
    if (!job || job.id !== d.id || job.done) return;
    const r = d.result || errorResult('empty reply from the 3D worker');
    r.worker = true;
    if (d.type === 'progress') progress(job, r);
    else if (d.type === 'done') finishJob(job, r);
  }
  function onWorkerError(ev) {
    if (ev && ev.preventDefault) ev.preventDefault();
    const w = S.worker;
    S.worker = null;
    try { w && w.terminate(); } catch (e) {}
    const job = S.active;
    if (job && !job.done) finishJob(job, errorResult('the 3D worker stopped: ' + ((ev && ev.message) || 'unknown error'), { worker: true }));
  }
  async function runOnMainThread(job) {
    let b;
    try {
      b = mainCore().newBuild(job.molblock, job.opts, { maxHeavyAtoms: options.mainThreadMaxHeavyAtoms });
      for (;;) {
        if (job.done) return;                                     // cancelled between steps
        const more = b.step();
        if (job.done) return;
        if (b.takeDirty()) { const r = b.snapshot(true); if (r.ok) progress(job, r); }
        if (!more) break;
        await new Promise(res => setTimeout(res, 0));             // let the page breathe between conformers
      }
      finishJob(job, b.snapshot(false));
    } catch (e) { finishJob(job, Object.assign(mainCore().failure(e), { worker: false })); }
  }

  function cancelJob(job) {
    if (job.done) return false;
    const qi = S.queue.indexOf(job);
    if (qi >= 0) S.queue.splice(qi, 1);
    else if (S.active === job && job.worker && S.worker === job.worker) {
      S.worker = null;                                            // a Worker cannot be interrupted mid-step: end it
      try { job.worker.terminate(); } catch (e) {}                // (a new one is made for the next job)
    }
    const r = job.last ? Object.assign({}, job.last, { partial: false, cancelled: true, search: Object.assign({}, job.last.search, { stoppedBy: 'cancelled' }) })
                       : errorResult('cancelled', { cancelled: true });
    finishJob(job, r);
    return true;
  }

  /* Build a 3D model: first conformer (reported through onProgress as soon as it exists), then up to k more within
     timeBudgetMs, each MMFF94-minimised and stereo-checked, duplicates removed, sorted by energy. The promise always
     resolves (never rejects); after cancel() it resolves with the conformers found so far and cancelled: true. */
  function build(molblock, opts) {
    opts = opts || {};
    const job = {
      id: ++S.seq, molblock: String(molblock || ''), done: false, last: null, worker: null,
      onProgress: typeof opts.onProgress === 'function' ? opts.onProgress : null,
      opts: { k: opts.k == null ? 10 : opts.k, timeBudgetMs: opts.timeBudgetMs == null ? 20000 : opts.timeBudgetMs,
              seed: opts.seed == null ? null : opts.seed, previous: opts.previous || null },
    };
    job.promise = new Promise(res => { job.resolve = res; });
    if (!job.molblock.trim()) { job.done = true; job.resolve(errorResult('no structure given')); }
    else { S.queue.push(job); pump(); }
    return { promise: job.promise, cancel: () => cancelJob(job) };
  }

  // Single-point MMFF94 energy (kcal/mol) of 3D coordinates with every H. Runs on the page's own copy of OCL:
  // one force-field set-up, not a search, so it does not wait behind a running build.
  function energy(molblock) {
    return loadLibrary().then(() => mainCore().energy(String(molblock || '')))
      .catch(e => ({ ok: false, energy: NaN, error: message(e) }));
  }

  /* ================= writers (export; also how tests read the models back) ================= */
  // record: MolInfo's record ({atoms:[{el, charge, isotope, radical}], bonds:[{a, b, order}]}) or a build Result
  // ({elements, charges, isotopes, radicals, bonds:[[a, b, order]]}). radical = radical electrons (RDKit's count).
  function topology(rec) {
    if (rec && Array.isArray(rec.atoms)) {
      return { N: rec.atoms.length, atoms: rec.atoms.map(a => ({ el: a.el, charge: a.charge | 0, isotope: a.isotope | 0, radical: a.radical | 0 })),
               bonds: (rec.bonds || []).map(b => [b.a, b.b, b.order, !!b.aromatic]) };
    }
    if (rec && rec.elements) {
      const g = (arr, i) => (arr ? arr[i] | 0 : 0);
      return { N: rec.elements.length, atoms: rec.elements.map((el, i) => ({ el, charge: g(rec.charges, i), isotope: g(rec.isotopes, i), radical: g(rec.radicals, i) })),
               bonds: (rec.bonds || []).map(b => [b[0], b[1], b[2], false]) };
    }
    throw new Error('molfile writer: no atoms in the record');
  }
  const pad = (v, w) => String(v).padStart(w);
  const oneLine = s => String(s == null ? '' : s).replace(/[\r\n]+/g, ' ').slice(0, 80);
  function propLines(tag, list) {                           // M  CHG / M  ISO / M  RAD, at most 8 entries per line
    const out = [];
    for (let i = 0; i < list.length; i += 8) {
      const part = list.slice(i, i + 8);
      out.push(`M  ${tag}${pad(part.length, 3)}` + part.map(([a, v]) => ` ${pad(a + 1, 3)} ${pad(v, 3)}`).join(''));
    }
    return out;
  }
  /* The counts line's chiral flag. 0 makes OpenChemLib (and ChemDraw, Marvin) read every stereocentre as relative or
     racemic ("both enantiomers"), so it is 1 when the record's stereocentres are all specified: a MolInfo record with
     at least one R/S/r/s centre and no '?', or a build Result with checked stereo and none left open. It stays 0 when
     the input left a centre open (the model's configuration there was chosen by the builder). */
  function chiralFlag(rec) {
    const st = rec && rec.stereo;
    if (!st) return false;
    if (Array.isArray(st.atoms)) return st.atoms.some(s => /^[RSrs]$/.test(s.label)) && !st.atoms.some(s => s.label === '?');
    return st.checked > 0 && Array.isArray(st.unspecified) && !st.unspecified.some(u => u.atom != null);
  }
  function molfile(record, coords, title, opt) {
    const T = topology(record), N = T.N;
    const chiral = opt && opt.chiral != null ? !!opt.chiral : chiralFlag(record);
    if (N > 999 || T.bonds.length > 999) throw new Error('a V2000 molfile holds at most 999 atoms and 999 bonds');
    if (!coords || coords.length < 3 * N) throw new Error('molfile writer: coordinates missing');
    const d = new Date(), two = v => String(v).padStart(2, '0');
    const stamp = two(d.getMonth() + 1) + two(d.getDate()) + two(d.getFullYear() % 100) + two(d.getHours()) + two(d.getMinutes());
    const f = v => (+v).toFixed(4).padStart(10);
    const L = [oneLine(title), '  MolView ' + stamp + '3D', oneLine(opt && opt.comment),
               `${pad(N, 3)}${pad(T.bonds.length, 3)}  0  0${pad(chiral ? 1 : 0, 3)}  0  0  0  0  0999 V2000`];
    const chg = [], iso = [], rad = [];
    T.atoms.forEach((a, i) => {
      L.push(`${f(coords[3 * i])}${f(coords[3 * i + 1])}${f(coords[3 * i + 2])} ${String(a.el).padEnd(3)} 0  0  0  0  0  0  0  0  0  0  0  0`);
      if (a.charge) chg.push([i, a.charge]);
      if (a.isotope) iso.push([i, a.isotope]);
      if (a.radical) rad.push([i, a.radical === 1 ? 2 : 3]);   // 1 electron = doublet (2); 2 = triplet (3), as RDKit writes
    });
    for (const [a, b, order, arom] of T.bonds) L.push(`${pad(a + 1, 3)}${pad(b + 1, 3)}${pad(order >= 1 && order <= 3 ? order : (arom ? 4 : 1), 3)}  0`);
    L.push(...propLines('CHG', chg), ...propLines('ISO', iso), ...propLines('RAD', rad), 'M  END');
    return L.join('\n') + '\n';
  }
  function xyz(record, coords, title) {
    const T = topology(record), g = v => (+v).toFixed(5).padStart(11);
    const lines = [String(T.N), oneLine(title)];
    T.atoms.forEach((a, i) => lines.push(`${String(a.el).padEnd(2)} ${g(coords[3 * i])} ${g(coords[3 * i + 1])} ${g(coords[3 * i + 2])}`));
    return lines.join('\n') + '\n';
  }
  function sdf(record, conformers, title, opt) {
    const single = topology(record).N === 1;             // a lone atom: no MMFF94 energy to report (see singleAtom)
    return (conformers || []).map((c, k) => {
      const items = [];
      if (single) items.push(['NOTE', 'single atom: no geometry to optimise']);
      else if (c.optimised && isFinite(c.energy)) {
        items.push(['MMFF94_ENERGY_KCAL', (+c.energy).toFixed(4)]);
        if (isFinite(c.relEnergy)) items.push(['MMFF94_REL_ENERGY_KCAL', (+c.relEnergy).toFixed(4)]);
      } else items.push(['MMFF94_OPTIMISED', 'no (approximate geometry)']);
      if (c.id != null) items.push(['CONFORMER_ID', String(c.id)]);
      const comment = single ? 'single atom: no geometry to optimise' : c.optimised ? 'OpenChemLib conformer, MMFF94-minimised' : 'OpenChemLib conformer, not force-field optimised';
      return molfile(record, c.coords, `${oneLine(title)}${title ? ' ' : ''}conformer ${k + 1}`, Object.assign({}, opt, { comment })) +
        items.map(([k2, v]) => `> <${k2}>\n${v}\n`).join('\n') + '\n$$$$\n';
    }).join('');
  }

  function status() { return { mode: S.mode, workerError: S.workerError, busy: !!S.active, queued: S.queue.length }; }

  return { build, energy, molfile, xyz, sdf, chiralFlag, ready, status, options,
           _engineCore: engineCore, _workerMain: workerMain, _workerSource: workerSource };
})();
