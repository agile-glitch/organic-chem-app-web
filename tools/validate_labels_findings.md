## Findings (2026-10-14)

**Values that follow from the structure are right.** Formal charge, lone pairs, oxidation state (every atom whose
value does not depend on which Kekulé structure is drawn), Gasteiger–Marsili partial charge, the RDKit hybridization
port, R/S and E/Z agree with RDKit Python on every one of 116,300 atoms / 1,357 stereocentres / 87 double bonds of
8,936 molecules. The displayed hybridization agrees with a textbook steric-number rule on 99.5 % of atoms; the rest
are conventions (below). 132 of 132 hand-checked textbook values agree (130 before fix 4 below). Bond angles of the app's 3D model agree with RDKit's MMFF94 geometry within 5° (218 angles, 40 molecules).

**Problems found, most important first:**

1. **Steric, electronic and Taft labels are dropped although there is room** — 26 %, 20 % and 36 % of the labels that
   have a value were not drawn ("left out where there is no free space"). The placement kept a label's centre
   clearly nearer its own atom than any other (`placeNotes`, mol2d_view.js), and these texts are long; on a carbonyl
   carbon with three neighbours no position passed. (Aspirin, Taft mode: both labels left out on an empty canvas.)
   **Fixed (2026-10-14):** a long label may also pass by the edge it is read from (that edge within half a bond of its
   atom and clearly nearer it than any other atom), may start or end at its atom and run outward, and is looked for up
   to 0.7 bond further out. Left out now 12 % / 10 % / 17 % (the rest is real crowding: CCl4's carbon between four Cl,
   several long labels on neighbouring atoms); labels whose nearest edge is nearer another atom went down
   (steric 39 → 31, Taft 21 → 10); short labels (hybridization, charges) are placed as before.
2. **pKaH of anilines the Hammett estimate cannot reach: 10.7 instead of about 4.6.** The primary-amine fragment
   (F002, 10.7) also matches an aryl NH2, and a basic site takes the most basic of its fragments, so it beats the
   aniline fragment (F026, 4.6) whenever there is an ortho substituent, a fused ring or a heteroaryl ring
   (2,6-dibromoaniline: measured 0.17, rules 10.7). 228 measured molecules: median error 6.4.
3. **pKaH of N,N-dialkylanilines: 10 instead of about 5.** Same pattern: the trialkylamine fragment (F093, 10.0)
   beats F180 (5.1). 76 measured molecules: median error 6.1, none within 1 unit.
4. **Measured pKa table includes high-pressure measurements.** tools/pka.py took the median of every 20–30 °C value,
   including ones at 1,000–3,000 bar: acetic acid showed 4.58 (4.76 at 1 bar), benzoic acid 4.02 (4.21). 11 molecules.
   **Fixed (2026-10-14):** measurements above 1 bar are left out; the 11 molecules now have their 1-bar values.
5. **Gaps in the pKa rules** (a site the rules miss or value they cannot know; shown only when no measured value
   exists): a CH between two C=O (active methine, about 11–13; phenylbutazone measured 2.7, rules 29); an N–H between
   C=O and SO2 (acyl sulfonamide; saccharin measured 1.3, rules 17); xanthic / dithiocarboxylic acid S–H (about 1.5,
   rules 13); N-substituted amidines and guanidines (pKaH 12–14, rules 7 as an imine); gem-dinitro / polynitro C–H.
   Overall the rules alone are within 1 unit for 42 % of acids and 35 % of bases (median error 1.3 / 1.7); the
   Hammett estimates are within 1 unit for 100 % / 93 %.
6. **Labels inside rings**: 5–8 % of per-atom labels (14 % of R/S labels) sit inside a ring — still nearest their own
   atom, but harder to read (your aspirin screenshots: C8's label inside the benzene ring).

**Conventions, not errors** (worth an ⓘ note, or a decision):
- Oxidation state of an atom with an aromatic bond to another element (pyridine C2, pyrimidine N) is the average over
  the Kekulé structures, so it can be fractional (+0.5); RDKit and most textbooks read one Kekulé structure (72 % of
  these atoms happen to agree).
- Hybridization: an S or Se on an aromatic ring (thiophenol, thioanisole) stays sp3; a lone pair on N or O next to
  C=S (thioamide, xanthate) is sp2 like an amide. Halogens are labelled by convention.
- Atoms with no simple electron count (Na⁺, Zn²⁺) get no lone-pair label.
- Steric classes saturate: isobutyraldehyde and pivaldehyde are both "hindered"; propanal is already "moderately
  hindered".
