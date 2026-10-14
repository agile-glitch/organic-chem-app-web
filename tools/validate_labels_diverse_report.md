# Label validation (Molecule tab)

Made by `tools/validate_labels.js` + `tools/validate_labels.py` from `data/label_dump_diverse.json`: 36,708 molecules analysed, 0 not.

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
   **Fixed (2026-10-14):** a catch-all amine fragment gives no value when an aryl-amine fragment is on the same N
   (SITE_RULES.baseYield in js/reactions.js, base_yield in tools/pka.py). Median error 6.4 → 0.7, within 1 unit 35 % → 57 %.
3. **pKaH of N,N-dialkylanilines: 10 instead of about 5.** Same pattern: the trialkylamine fragment (F093, 10.0)
   beats F180 (5.1). 76 measured molecules: median error 6.1, none within 1 unit.
   **Fixed (2026-10-14)** with 2: median error 6.1 → 1.6, within 1 unit 0 % → 32 %.
4. **Measured pKa table includes high-pressure measurements.** tools/pka.py took the median of every 20–30 °C value,
   including ones at 1,000–3,000 bar: acetic acid showed 4.58 (4.76 at 1 bar), benzoic acid 4.02 (4.21). 11 molecules.
   **Fixed (2026-10-14):** measurements above 1 bar are left out; the 11 molecules now have their 1-bar values.
5. **Gaps in the pKa rules** (a site the rules miss or value they cannot know; shown only when no measured value
   exists): a CH between two C=O (active methine, about 11–13; phenylbutazone measured 2.7, rules 29); an N–H between
   C=O and SO2 (acyl sulfonamide; saccharin measured 1.3, rules 17); xanthic / dithiocarboxylic acid S–H (about 1.5,
   rules 13); N-substituted amidines and guanidines (pKaH 12–14, rules 7 as an imine); gem-dinitro / polynitro C–H.
   **Fixed (2026-10-14):** site rules for these groups (SITE_RULES.extra in js/reactions.js, extra in tools/pka.py),
   values from the measured data where the group is the only ionisable site, else textbook values: acyclic C–H
   between two C=O / C≡N 11, the same C–H in a 5-/6-ring 5, acyl sulfonamide N–H 4.5, xanthic S–H 2.0,
   dithiocarbamic 4.0, dithiocarboxylic 2.5, thioic S-acid 3.5, gem-dinitro C–H 3.8, HC(NO2)3 0.1, guanidine
   pKaH 13.6, N-aryl guanidine 10.8, amidine 11.5. Rules alone against measured: xanthic acids median error 0.5
   (was about 11), guanidines 0.2 (about 6.5), amidines 0.8 (about 4.5), gem-dinitro 1.15 (about 6), active
   methine 1.5 (7–25). Found while retesting and fixed with them: the imine fragments F196 / F204 (pKa 11, an
   oxime's O–H) put a "pKa 11" on the =NH of guanidines and amidines; they no longer mark a site.
   Overall the rules alone are now within 1 unit for 46 % of acids and 41 % of bases (median error 1.14 / 1.32;
   42 % / 35 % and 1.3 / 1.7 before fixes 2, 3 and 5); the Hammett estimates are within 1 unit for 100 % / 93 %.
   17 textbook values for these groups added to tools/label_textbook.json: 149 of 149 agree.
6. **Labels inside rings, and labels whose edge sits by another atom**: 5–8 % of per-atom labels (14 % of R/S labels)
   sat inside a ring (aspirin: C8's label inside the benzene ring, beside C3's), and up to 6 % of pKa / steric labels
   had their nearest edge closer to another atom than to their own. The cause of the first: for an atom whose bonds
   balance out (a ring carbon with a substituent) the "open side" defaulted to straight down, often into the ring.
   **Fixed (2026-10-14):** such an atom's open side is away from its rings; a place inside a ring costs more; and every
   label's nearest edge must be clearly nearer its own atom (its written label: OH, NH2) than any other atom. Inside a
   ring now 0.3–3.5 % (only where nothing outside is free: fused-ring atoms); edge nearer another atom 0–0.6 % (was up
   to 6 %); labels left out about the same (pKa 7 → 2).

**Conventions, not errors** (worth an ⓘ note, or a decision):
- Oxidation state of an atom with an aromatic bond to another element (pyridine C2, pyrimidine N) is the average over
  the Kekulé structures, so it can be fractional (+0.5); RDKit and most textbooks read one Kekulé structure (72 % of
  these atoms happen to agree).
- Hybridization: an S or Se on an aromatic ring (thiophenol, thioanisole) stays sp3; a lone pair on N or O next to
  C=S (thioamide, xanthate) is sp2 like an amide. Halogens are labelled by convention.
- Atoms with no simple electron count (Na⁺, Zn²⁺) get no lone-pair label.
- Steric classes saturate: isobutyraldehyde and pivaldehyde are both "hindered"; propanal is already "moderately
  hindered".

## Agreement with independent references

| Label | Reference | Checked | Agree | |
|---|---|---:|---:|---:|
| formal | formal charge vs RDKit Python | 593,502 | 593,502 | 100.00 % |
| lonepairs | lone pairs vs the electron count | 593,462 | 593,080 | 99.94 % |
| oxidation | oxidation state vs RDKit CalcOxidationNumbers (atoms whose value does not depend on the Kekulé structure) | 535,149 | 535,132 | 100.00 % |
| oxidation_kekule | oxidation state, atoms whose value depends on the Kekulé structure (the app averages; RDKit takes one structure) | 58,344 | 44,915 | 76.98 % |
| partial | Gasteiger-Marsili charge vs RDKit Python (difference below 0.005; the label shows 2 decimals) | 588,681 | 588,681 | 100.00 % |
| hybridization | displayed hybridization vs the textbook steric-number rule | 573,482 | 571,484 | 99.65 % |
| rdkit_hybridization | the app's port of RDKit's hybridization vs RDKit Python | 593,502 | 593,498 | 100.00 % |
| cip | R/S vs RDKit Python CIP labeller (atoms either one labels) | 6,697 | 6,697 | 100.00 % |
| ez | E/Z vs RDKit Python CIP labeller (bonds either one labels) | 609 | 609 | 100.00 % |
| oxidation_sum | self-consistency: the oxidation states of all atoms (H included) add up to the net charge | 36,705 | 36,538 | 99.55 % |
| partial_sum | self-consistency: the partial charges of all atoms (H included) add up to the net charge (within 0.01) | 36,333 | 36,332 | 100.00 % |
| pka_acid_site | self-consistency: an acid pKa sits on an atom that has an H | 32,124 | 32,121 | 99.99 % |
| pka_base_site | self-consistency: a pKaH sits on a neutral N, O, P or S with a lone pair | 33,070 | 33,070 | 100.00 % |
| pka_range | self-consistency: every pKa / pKaH is a number between -12 and 52 | 62,636 | 62,636 | 100.00 % |
| electronic_partial | self-consistency: the electronic label's δ is the partial charge of the same atom | 363,098 | 363,098 | 100.00 % |
| steric_class | self-consistency: the steric word matches its crowding score (≤2 open, 3 moderately, 4 hindered, ≥5 very) | 365,170 | 365,170 | 100.00 % |

### Agreement by kind of molecule

Percent of checked atoms (or sites, bonds, molecules) that agree; — = nothing of that kind to check. Molecules: ordinary organic 15,541, heteroaromatic 11,736, 3+ aromatic rings 4,459, stereocentres 2,774, hypervalent S / P / halogen 2,649, zwitterion or dipole 2,532, tiny (6 or fewer heavy atoms) 1,804, small ring (3-4) 1,285, ion 1,003, large (40+ heavy atoms) 619, E/Z double bond 501, metal or other element 449, macrocycle (ring of 8+) 349, cumulene / allene / ketene 214, salt or mixture 149, radical or carbene 82, isotope label 12.

| Kind of molecule | formal | lonepairs | oxidation | oxidation_sum | partial | partial_sum | hybridization | cip | ez | pka_acid_site | pka_base_site | electronic_partial |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heteroaromatic | 100 | 99.99 | 100.00 | 98.58 | 100 | 99.99 | 99.61 | 100 | 100 | 100 | 100 | 100 |
| ordinary organic | 100 | 100 | 100 | 100 | 100 | 100 | 99.63 | — | — | 100 | 100 | 100 |
| 3+ aromatic rings | 100 | 99.98 | 100.00 | 98.79 | 100 | 100 | 99.82 | 100 | 100 | 100 | 100 | 100 |
| stereocentres | 100 | 100 | 100.00 | 99.75 | 100 | 100 | 99.73 | 100 | 100 | 100 | 100 | 100 |
| hypervalent S / P / halogen | 100 | 99.95 | 100.00 | 99.96 | 100 | 100 | 99.81 | 100 | 100 | 100 | 100 | 100 |
| zwitterion or dipole | 100 | 100.00 | 100 | 99.92 | 100 | 99.96 | 99.78 | 100 | 100 | 100 | 100 | 100 |
| large (40+ heavy atoms) | 100 | 99.93 | 100.00 | 100 | 100 | 100 | 99.76 | 100 | 100 | 100 | 100 | 100 |
| small ring (3-4) | 100 | 99.99 | 100 | 99.92 | 100 | 100 | 99.47 | 100 | 100 | 100 | 100 | 100 |
| ion | 100 | 99.30 | 99.99 | 100 | 100 | 100 | 99.59 | 100 | 100 | 100 | 100 | 100 |
| E/Z double bond | 100 | 99.94 | 100 | 100 | 100 | 100 | 99.77 | 100 | 100 | 100 | 100 | 100 |
| tiny (6 or fewer heavy atoms) | 100 | 97.37 | 99.80 | 100 | 100 | 100 | 99.11 | 100 | 100 | 99.62 | 100 | 100 |
| macrocycle (ring of 8+) | 100 | 100 | 100 | 100 | 100 | 100 | 99.92 | 100 | 100 | 100 | 100 | 100 |
| metal or other element | 100 | 91.78 | 99.98 | 100 | 100 | 100 | 99.87 | 100 | 100 | 100 | 100 | 100 |
| cumulene / allene / ketene | 100 | 100 | 100 | 100 | 100 | 100 | 99.84 | 100 | — | 100 | 100 | 100 |
| salt or mixture | 100 | 99.45 | 100 | 100 | 100 | 100 | 99.79 | 100 | 100 | 100 | 100 | 100 |
| radical or carbene | 100 | 93.44 | 100 | 100 | 100 | 100 | 99.85 | — | — | 100 | 100 | 100 |
| isotope label | 100 | 100 | 74.60 | 100 | 100 | 100 | 100 | 100 | — | 57.14 | — | 100 |

### lonepairs: disagreements by kind

- **unusual electron count**: 381 — e.g. `CCN(CC)c1ccc2nc3cc(C)c(N)cc3[o+]c2c1.CCN(CC)c1ccc2nc3cc(C)c(N)cc3[o+]c2c1.[Cl-].[Cl-].[Cl-].[Cl-].[Zn+2] atom 46 (Zn): app None, count 0`; `CN(C)c1ccc(/N=N/c2ccc(S(=O)(=O)[O-])cc2)cc1.[Na+] atom 21 (Na): app None, count 0`; `O=S(=O)(O)c1ccc(N=Nc2ccc(Nc3ccccc3)cc2)cc1.[Na] atom 25 (Na): app None, count 0`
- **H 0 σ-bonds charge -1**: 1 — e.g. `[H-] atom 0 (H): app 1, count 0`

### oxidation: disagreements by kind

- **C(arom) 3 σ-bonds, next to a π bond**: 6 — e.g. `[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H] atom 1 (C): app -1, RDKit 0`; `[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H] atom 2 (C): app -1, RDKit 0`; `[2H]c1c([2H])c([2H])c([2H])c([2H])c1[2H] atom 4 (C): app -1, RDKit 0`
- **C 4 σ-bonds**: 5 — e.g. `[2H]C([2H])([2H])O atom 1 (C): app -2, RDKit 1`; `[2H][C@@H](C)O atom 1 (C): app -1, RDKit 0`; `[2H]C(Cl)(Cl)Cl atom 1 (C): app 2, RDKit 3`
- **O 2 σ-bonds**: 2 — e.g. `[2H]O[2H] atom 1 (O): app -2, RDKit 0`; `[3H]O atom 1 (O): app -2, RDKit -1`
- **C 4 σ-bonds, next to a π bond**: 2 — e.g. `[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H] atom 1 (C): app -3, RDKit 0`; `[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H] atom 6 (C): app -3, RDKit 0`
- **Ir 6 σ-bonds charge +1, next to a π bond**: 1 — e.g. `CC(C)(C)c1cc[n]2->[Ir+]34([c]5cc(F)cc(F)c5-c5ccc(C(F)(F)F)c[n]->35)([c]3cc(F)cc(F)c3-c3ccc(C(F)(F)F)c[n]->43)<-[n]3ccc(C(C)(C)C)cc3-c2c1 atom 8 (Ir): app 71, RDKit 3`
- **O 2 σ-bonds, next to a π bond**: 1 — e.g. `CC(=O)O[2H] atom 3 (O): app -2, RDKit -1`

### oxidation_kekule: disagreements by kind

- **C(arom) 3 σ-bonds, next to a π bond**: 13,091 — e.g. `Nc1ccc(S(=O)(=O)Nc2ccccn2)cc1 atom 9 (C): app 2.5, RDKit 2`; `Nc1ccc(S(=O)(=O)Nc2ccccn2)cc1 atom 13 (C): app 0.5, RDKit 1`; `c1cncc([C@@H]2CCCCN2)c1 atom 1 (C): app 0.5, RDKit 0`
- **N(arom) 2 σ-bonds, next to a π bond**: 333 — e.g. `c1ccnnc1 atom 3 (N): app -1.5, RDKit -2`; `c1ccnnc1 atom 4 (N): app -1.5, RDKit -2`; `Nc1cccnn1 atom 5 (N): app -1.5, RDKit -1`
- **N(arom) 3 σ-bonds, next to a π bond**: 4 — e.g. `CC(C)(C)c1cc[n]2->[Ir+]34([c]5cc(F)cc(F)c5-c5ccc(C(F)(F)F)c[n]->35)([c]3cc(F)cc(F)c3-c3ccc(C(F)(F)F)c[n]->43)<-[n]3ccc(C(C)(C)C)cc3-c2c1 atom 7 (N): app -20, RDKit -3`; `CC(C)(C)c1cc[n]2->[Ir+]34([c]5cc(F)cc(F)c5-c5ccc(C(F)(F)F)c[n]->35)([c]3cc(F)cc(F)c3-c3ccc(C(F)(F)F)c[n]->43)<-[n]3ccc(C(C)(C)C)cc3-c2c1 atom 26 (N): app -20, RDKit -3`; `CC(C)(C)c1cc[n]2->[Ir+]34([c]5cc(F)cc(F)c5-c5ccc(C(F)(F)F)c[n]->35)([c]3cc(F)cc(F)c3-c3ccc(C(F)(F)F)c[n]->43)<-[n]3ccc(C(C)(C)C)cc3-c2c1 atom 44 (N): app -20, RDKit -3`
- **N(arom) 3 σ-bonds charge +1, next to a π bond**: 1 — e.g. `[O-][n+]1nc(Nc2ccc(Cl)cc2)c2ccccc2c1Cc1ccncc1 atom 1 (N): app 0.333333, RDKit 0`

### hybridization: disagreements by kind

- **S 2 σ-bonds, next to a π bond: app sp3, rule sp2**: 1,551 — e.g. `Nc1ccccc1S atom 7 (S)`; `CSc1cccc(N)c1 atom 1 (S)`; `CSc1ccc(N)cc1 atom 1 (S)`
- **N 3 σ-bonds: app sp2, rule sp3**: 280 — e.g. `CC(N)=S atom 2 (N)`; `NC(N)=S atom 0 (N)`; `NC(N)=S atom 2 (N)`
- **P 3 σ-bonds, next to a π bond: app sp3, rule sp2**: 77 — e.g. `OP(O)c1cccc(Cl)c1 atom 1 (P)`; `OP(O)c1ccc(Cl)cc1 atom 1 (P)`; `OP(O)c1cccc(Br)c1 atom 1 (P)`
- **O 2 σ-bonds: app sp2, rule sp3**: 40 — e.g. `COC(=S)S atom 1 (O)`; `C=CCOC(=S)S atom 3 (O)`; `CC(C)OC(=S)S atom 3 (O)`
- **Se 2 σ-bonds, next to a π bond: app sp3, rule sp2**: 23 — e.g. `C[Se]c1cccc2cccnc12 atom 1 (Se)`; `[SeH]c1cccc2cccnc12 atom 0 (Se)`; `O=C(O)C[Se]c1ccccc1 atom 4 (Se)`
- **S 3 σ-bonds charge +1, next to a π bond: app sp3, rule sp2**: 10 — e.g. `C[S+](C)c1cccc(O)c1 atom 1 (S)`; `C[S+](C)c1ccc(O)cc1 atom 1 (S)`; `Oc1cccc([SH2+])c1.[I-] atom 6 (S)`
- **S 1 σ-bonds charge -1: app sp2, rule sp3**: 8 — e.g. `CC(=S)[S-] atom 3 (S)`; `CCOC(=S)[S-] atom 5 (S)`; `O=S(=O)([SH-]C(=S)[S-])C(Cl)(Cl)Cl atom 6 (S)`
- **O 1 σ-bonds charge -1: app sp2, rule sp3**: 5 — e.g. `O=[Si]([O-])[O-] atom 2 (O)`; `O=[Si]([O-])[O-] atom 3 (O)`; `[O-]C(=S)c1ccccc1 atom 0 (O)`
- **As 3 σ-bonds, next to a π bond: app sp3, rule sp2**: 4 — e.g. `O=C(O)C[As](CC(=O)O)c1ccccc1 atom 4 (As)`; `CC(=O)Nc1cc([As]2SCC(CO)S2)ccc1O atom 7 (As)`; `O=C(O)CC[As](CCC(=O)O)c1ccc(Cl)cc1 atom 5 (As)`

### rdkit_hybridization: disagreements by kind

- **app UNSPECIFIED, RDKit sp2**: 4 — e.g. `CC(C)(C)c1cc[n]2->[Ir+]34([c]5cc(F)cc(F)c5-c5ccc(C(F)(F)F)c[n]->35)([c]3cc(F)cc(F)c3-c3ccc(C(F)(F)F)c[n]->43)<-[n]3ccc(C(C)(C)C)cc3-c2c1 atom 7 (N)`; `CC(C)(C)c1cc[n]2->[Ir+]34([c]5cc(F)cc(F)c5-c5ccc(C(F)(F)F)c[n]->35)([c]3cc(F)cc(F)c3-c3ccc(C(F)(F)F)c[n]->43)<-[n]3ccc(C(C)(C)C)cc3-c2c1 atom 26 (N)`; `CC(C)(C)c1cc[n]2->[Ir+]34([c]5cc(F)cc(F)c5-c5ccc(C(F)(F)F)c[n]->35)([c]3cc(F)cc(F)c3-c3ccc(C(F)(F)F)c[n]->43)<-[n]3ccc(C(C)(C)C)cc3-c2c1 atom 44 (N)`

### oxidation_sum: disagreements by kind

- **sum differs**: 167 — e.g. `Oc1nc2cccnc2nc1O: Σ 1.0000000010279564e-06, charge 0`; `Oc1ccnc2nccnc12: Σ 1.0000000010279564e-06, charge 0`; `c1cnc2ncncc2n1: Σ 1.000000000139778e-06, charge 0`

### partial_sum: disagreements by kind

- **sum differs**: 1 — e.g. `O=[N+]([O-])c1c[nH]c([N+](=O)[O-])n1: Σ -0.500, charge 0`

### pka_acid_site: disagreements by kind

- **no H on the atom**: 3 — e.g. `[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H] atom 1 (C): pKa 19 from F019+F101+F192`; `[2H]C([2H])([2H])C(=O)C([2H])([2H])[2H] atom 6 (C): pKa 19 from F019+F101+F192`; `CC(=O)O[2H] atom 3 (O): pKa 4.8 from F001`

## pKa: the rules alone against measured values (IUPAC, water)

Molecules with one ionisable group in all: the app finds one site (acid or base) and the dataset has one measured value, of the same kind, so the value belongs to that site. The app works its value out with no measured value to fall back on (on screen a measured value, when there is one, is shown instead).

**acid**: 1,171 molecules; median error 1.14, within 1 unit 46%, within 2 units 70%, mean signed error +1.58.

| | n | median error | within 1 | within 2 | mean signed |
|---|---:|---:|---:|---:|---:|
| Hammett estimate | 96 | 0.08 | 100% | 100% | +0.05 |
| typical value | 1,075 | 1.28 | 41% | 68% | +1.72 |
| F028 | 357 | 1.72 | 38% | 55% | +2.11 |
| F001 | 308 | 1.06 | 46% | 83% | +1.01 |
| F001+F118 | 139 | 0.5 | 70% | 88% | +0.66 |
| F056+gem-dinitro C–H (C(NO2)2) | 45 | 1.12 | 47% | 87% | -0.21 |
| F021+F095+F121+F185+F187+F202 | 33 | 0.39 | 91% | 94% | +0.21 |
| F021+F095+F123 | 27 | 0.72 | 70% | 82% | +0.99 |
| F072 | 24 | 0.64 | 71% | 96% | +0.57 |
| F123 | 23 | 1.32 | 39% | 87% | +1.42 |
| F003+F014 | 15 | 1.69 | 40% | 60% | +2.27 |
| F063 | 14 | 3.45 | 21% | 36% | +2.36 |
| F076+xanthic acid S–H (RO–C(=S)–SH) | 13 | 0.48 | 54% | 100% | -0.52 |
| F068 | 11 | 0.52 | 73% | 91% | +0.61 |
| F056 | 10 | 1.94 | 30% | 50% | +1.24 |
| F076+dithiocarbamic acid S–H (R2N–C(=S)–SH) | 10 | 1.21 | 40% | 80% | -1.28 |
| F021+F095+F121+F185+F186+F187+F202 | 10 | 1.36 | 10% | 70% | +0.41 |

Largest errors: `COC(=O)C[SeH]` measured 4.7, rules 25 (F101+F139+F193); `Cc1ccc2c(c1)S(=O)(=O)CS2(=O)=O` measured 12.63, rules 29 (F071+F101); `CC1=C(O)C(C)C1=O` measured 2.8, rules 19 (F019+F101+F192); `O=C1CS(=O)(=O)c2ccccc21` measured 2.88, rules 18.3 (F019+F071+F101+F116+F192); `N#CC1C=CC=C1` measured 9.78, rules 25 (F101+F194); `CCOC(=O)C1=C(O)CNC1=O` measured 2.34, rules 17 (F021+F095+F185+F186+F202).
Molecules with a measured acid value but no acid site found: 4 (Reliable), 206 (Uncertain), 118 (Approximate), 15 (Very uncertain), 7 (not stated).

**base**: 1,393 molecules; median error 1.32, within 1 unit 41%, within 2 units 66%, mean signed error +0.82.

| | n | median error | within 1 | within 2 | mean signed |
|---|---:|---:|---:|---:|---:|
| Hammett estimate | 81 | 0.22 | 93% | 100% | +0.07 |
| typical value | 1,312 | 1.44 | 38% | 64% | +0.86 |
| F091 | 303 | 2.1 | 30% | 48% | -0.09 |
| F002+F026 | 228 | 0.69 | 57% | 79% | +0.77 |
| F002+F176 | 164 | 1.11 | 46% | 84% | +1.22 |
| F093+F179 | 121 | 2.48 | 17% | 43% | +2.87 |
| F093 | 106 | 0.86 | 52% | 73% | +1.03 |
| F093+F180 | 76 | 1.61 | 32% | 67% | +0.74 |
| F113 | 60 | 0.9 | 53% | 88% | +0.12 |
| F034+F177 | 56 | 0.75 | 68% | 80% | +0.99 |
| F034+F178 | 55 | 0.6 | 60% | 82% | +1.03 |
| F094 | 53 | 0.92 | 55% | 74% | +0.68 |
| F112 | 40 | 1.67 | 25% | 68% | +1.67 |
| F031 | 34 | 1.83 | 9% | 56% | +1.06 |
| F031+F204 | 27 | 2.15 | 26% | 37% | -1.53 |
| F031+F196 | 17 | 2.38 | 12% | 12% | -0.02 |

Largest errors: `NC1=NN(c2ccccc2)C(=O)C1` measured 1.2, rules 12.4 (F031+F098); `NC1=CC(=O)C(=O)c2ccccc21` measured 1.57, rules 10.7 (F002+F127); `C#CCN(CC#C)CC#C` measured 3.09, rules 10.75 (F093+F179); `ClCCN(CCCl)Cc1ccc2ccccc2c1` measured 3.32, rules 10.75 (F093+F179); `NC1(N)C=CC(S(=O)(=O)O)=CC1` measured 3.29, rules 10.7 (F002+F176); `Fc1cccc(CN(CCCl)CCCl)c1` measured 3.35, rules 10.75 (F093+F179).
Molecules with a measured base value but no base site found: 60 (Approximate), 23 (not stated), 81 (Uncertain), 2 (Reliable), 4 (Very uncertain).

## The drawing (`tools/validate_display.js`)

Each label mode on the 2D drawing, as the tab shows it (Show H off). *Left out*: labels with a value that were not drawn ("no free space"). *Centre / edge nearer another atom*: the label's centre (or its nearest edge) is closer to another atom than to its own, so a reader can take it for that atom. *Inside a ring*: the label sits inside a ring. *Far*: more than one bond length from its atom.

| Mode | Molecules | Labels drawn | Left out | Molecules with labels left out | Centre nearer another atom | Edge nearer another atom | Inside a ring | Far |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| index | 181 | 1,908 | 6 (0.31 %) | 5 | 0 (0.00 %) | 5 (0.26 %) | 6 (0.31 %) | 0 |
| cip | 181 | 43 | 0 (0.00 %) | 0 | 0 (0.00 %) | 0 (0.00 %) | 0 (0.00 %) | 0 |
| hybridization | 181 | 1,906 | 8 (0.42 %) | 7 | 0 (0.00 %) | 11 (0.58 %) | 66 (3.46 %) | 0 |
| formal | 181 | 53 | 0 (0.00 %) | 0 | 0 (0.00 %) | 0 (0.00 %) | 0 (0.00 %) | 0 |
| partial | 181 | 1,888 | 11 (0.58 %) | 8 | 2 (0.11 %) | 6 (0.32 %) | 46 (2.44 %) | 0 |
| oxidation | 181 | 1,906 | 8 (0.42 %) | 7 | 0 (0.00 %) | 3 (0.16 %) | 7 (0.37 %) | 0 |
| lonepairs | 181 | 1,910 | 4 (0.21 %) | 4 | 1 (0.05 %) | 4 (0.21 %) | 2 (0.10 %) | 0 |
| pka | 181 | 327 | 2 (0.61 %) | 2 | 10 (3.06 %) | 0 (0.00 %) | 0 (0.00 %) | 29 |
| steric | 181 | 498 | 83 (14.29 %) | 65 | 81 (16.27 %) | 1 (0.20 %) | 15 (3.01 %) | 158 |
| electronic | 181 | 524 | 57 (9.81 %) | 39 | 43 (8.21 %) | 1 (0.19 %) | 11 (2.10 %) | 68 |
| taft | 181 | 277 | 49 (15.03 %) | 42 | 47 (16.97 %) | 1 (0.36 %) | 13 (4.69 %) | 38 |

Examples of labels nearer another atom: `O=C1NS(=O)(=O)c2ccccc21: '+0.26' for atom 3 sits nearer atom 6 (0.8 vs 0.67 bond lengths)`; `CCCCC[PH](=O)O: '+0.19' for atom 5 sits nearer atom 4 (0.92 vs 0.89 bond lengths)`

### Bond angles: the app's 3D model against RDKit (ETKDG + MMFF94)

218 angles between heavy atoms in 40 molecules: median difference 0.0°, 95 % within 2.35°, largest 4.35°; 8 differ by more than 3°, 0 by more than 5°.

Largest: `CC(=O)Oc1ccccc1C(=O)O` atoms (10, 9, 4): app 120.8°, RDKit 125.1°; `CC(=O)Oc1ccccc1C(=O)O` atoms (8, 9, 10): app 120.8°, RDKit 116.6°; `CC(=O)c1ccccc1` atoms (1, 3, 4): app 118.4°, RDKit 122.3°; `CC(=O)c1ccccc1` atoms (1, 3, 8): app 122.3°, RDKit 118.4°; `Oc1ccc([N+](=O)[O-])cc1` atoms (0, 1, 9): app 117.5°, RDKit 121.4°.

## Textbook list (`tools/label_textbook.json`)

149 of 149 hand-checked values agree.

## Steric class along crowding series

- aldehydes: ok — `CC=O` open (2), `CCC=O` moderately hindered (3), `CC(C)C=O` hindered (4), `CC(C)(C)C=O` hindered (4)
- alkyl bromides (SN2 order): ok — `CBr` open (1), `CCBr` open (2), `CC(C)Br` moderately hindered (3), `CC(C)(C)Br` hindered (4)
- amines: ok — `CCCCN` open (2), `CCC(C)N` moderately hindered (3), `CC(C)(C)N` hindered (4)
- ketones: ok — `CC(C)=O` moderately hindered (3), `CCC(C)=O` hindered (4), `CC(C)C(C)=O` very hindered (5), `CC(C)(C)C(C)=O` very hindered (5)
