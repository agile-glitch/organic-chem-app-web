# Label validation (Molecule tab)

Made by `tools/validate_labels.js` + `tools/validate_labels.py` from `data/label_dump_diverse.json`: 36,708 molecules analysed, 0 not.

## Findings on every kind of molecule (2026-10-14)

36,708 molecules: the IUPAC pKa set (8,894), the patent-reaction molecules (21,868: drug-like molecules, reagents,
catalysts), the named compounds (6,350), the reagent list (1,200) and 173 hand-picked unusual ones
(`tools/label_special.json`), in 17 kinds (table below). Every one was analysed by the tab; none failed.

**Agreement with independent references: complete wherever the answer is not a convention.** Formal charge,
partial charge, oxidation state, the RDKit hybridization port, R/S (6,697 centres) and E/Z (609 bonds) agree with RDKit
on every atom of every kind of molecule; oxidation states add up to the net charge in all 36,705 molecules, partial
charges to RDKit's own total in all 36,333; every pKa sits on an atom with an H (acid) or a neutral atom with a lone
pair (base) and lies in range; the electronic label's δ is always the atom's partial charge and the steric word always
its score; 149 of 149 hand-checked textbook values agree.

**Fixed while testing:** a dative bond (L→M, as in metal complexes written with arrows) was read as bond order 17, so
a metal amine's N showed oxidation state −20 and the metal +36 (an iridium photocatalyst +71). Now cisplatin shows
Pt +2 / N −3 and the N's lone pair counts as the bond's.

**Differences that are conventions, not errors:**
- Oxidation state next to an aromatic bond to another element is averaged over the Kekulé structures (pyridine C2
  +0.5); RDKit takes one structure (77 % of such atoms happen to agree).
- Metal atoms (Mg, Li, Cu, Pd, Zn, Na, Sn …, 385 atoms) get no lone-pair label: a lone pair is not a useful count for
  a metal.
- Hybridization (99.65 % agree with a plain steric-number rule): S, Se and P(III) on an aromatic ring stay sp³;
  N and O next to C=S (thioamides, thioureas, xanthates) are sp² like amides; O⁻ / S⁻ next to C=S or Si=O are sp².
- One molecule (2,4-dinitroimidazole) whose partial charges add to −0.5 instead of 0: RDKit's model does the same; the
  app ports it exactly.

**Bond angles** (55 molecules of every kind, 1,125 angles, the app's 3D model against RDKit MMFF94): at each centre
the angles agree to a median of 0.1°, 95 % within 2.9°; 12 differ by more than 5°, all in strained fused rings
(a cephalosporin-type C–S–C 105° vs 96.5°) or a charged dye — the two force-field programs, not the labels.

**The drawing** (250 molecules, every kind, every label mode): per-atom labels are left out for 0.3–1 % (charges,
oxidation states, lone pairs, hybridization), pKa 0.3 %; steric, electronic and Taft labels — long words — for
12–17 %, where the drawing is crowded (several on neighbouring atoms, metal centres). Labels whose edge is nearer
another atom: 0–1 %.

## Agreement with independent references

| Label | Reference | Checked | Agree | |
|---|---|---:|---:|---:|
| formal | formal charge vs RDKit Python | 593,502 | 593,502 | 100.00 % |
| lonepairs | lone pairs vs the electron count | 593,466 | 593,081 | 99.94 % |
| oxidation | oxidation state vs RDKit CalcOxidationNumbers (atoms whose value does not depend on the Kekulé structure) | 535,122 | 535,122 | 100.00 % |
| oxidation_kekule | oxidation state, atoms whose value depends on the Kekulé structure (the app averages; RDKit takes one structure) | 58,344 | 44,919 | 76.99 % |
| partial | Gasteiger-Marsili charge vs RDKit Python (difference below 0.005; the label shows 2 decimals) | 588,681 | 588,681 | 100.00 % |
| hybridization | displayed hybridization vs the textbook steric-number rule | 573,482 | 571,484 | 99.65 % |
| rdkit_hybridization | the app's port of RDKit's hybridization vs RDKit Python | 593,502 | 593,502 | 100.00 % |
| cip | R/S vs RDKit Python CIP labeller (atoms either one labels) | 6,697 | 6,697 | 100.00 % |
| ez | E/Z vs RDKit Python CIP labeller (bonds either one labels) | 609 | 609 | 100.00 % |
| oxidation_sum | self-consistency: the oxidation states of all atoms (H included) add up to the net charge | 36,705 | 36,705 | 100.00 % |
| partial_sum | the partial charges of all atoms (H included) add up to the same total as RDKit's model (within 0.01) | 36,333 | 36,333 | 100.00 % |
| pka_acid_site | self-consistency: an acid pKa sits on an atom that has an H | 32,124 | 32,124 | 100.00 % |
| pka_base_site | self-consistency: a pKaH sits on a neutral N, O, P or S with a lone pair | 33,070 | 33,070 | 100.00 % |
| pka_range | self-consistency: every pKa / pKaH is a number between -12 and 52 | 62,636 | 62,636 | 100.00 % |
| electronic_partial | self-consistency: the electronic label's δ is the partial charge of the same atom | 363,098 | 363,098 | 100.00 % |
| steric_class | self-consistency: the steric word matches its crowding score (≤2 open, 3 moderately, 4 hindered, ≥5 very) | 365,170 | 365,170 | 100.00 % |

### Agreement by kind of molecule

Percent of checked atoms (or sites, bonds, molecules) that agree; — = nothing of that kind to check. Molecules: ordinary organic 15,541, heteroaromatic 11,736, 3+ aromatic rings 4,459, stereocentres 2,774, hypervalent S / P / halogen 2,649, zwitterion or dipole 2,532, tiny (6 or fewer heavy atoms) 1,804, small ring (3-4) 1,285, ion 1,003, large (40+ heavy atoms) 619, E/Z double bond 501, metal or other element 449, macrocycle (ring of 8+) 349, cumulene / allene / ketene 214, salt or mixture 149, radical or carbene 82, isotope label 12.

| Kind of molecule | formal | lonepairs | oxidation | oxidation_sum | partial | partial_sum | hybridization | cip | ez | pka_acid_site | pka_base_site | electronic_partial |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| heteroaromatic | 100 | 99.99 | 100 | 100 | 100 | 100 | 99.61 | 100 | 100 | 100 | 100 | 100 |
| ordinary organic | 100 | 100 | 100 | 100 | 100 | 100 | 99.63 | — | — | 100 | 100 | 100 |
| 3+ aromatic rings | 100 | 99.97 | 100 | 100 | 100 | 100 | 99.82 | 100 | 100 | 100 | 100 | 100 |
| stereocentres | 100 | 100 | 100 | 100 | 100 | 100 | 99.73 | 100 | 100 | 100 | 100 | 100 |
| hypervalent S / P / halogen | 100 | 99.95 | 100 | 100 | 100 | 100 | 99.81 | 100 | 100 | 100 | 100 | 100 |
| zwitterion or dipole | 100 | 100.00 | 100 | 100 | 100 | 100 | 99.78 | 100 | 100 | 100 | 100 | 100 |
| large (40+ heavy atoms) | 100 | 99.91 | 100 | 100 | 100 | 100 | 99.76 | 100 | 100 | 100 | 100 | 100 |
| small ring (3-4) | 100 | 99.99 | 100 | 100 | 100 | 100 | 99.47 | 100 | 100 | 100 | 100 | 100 |
| ion | 100 | 99.28 | 100 | 100 | 100 | 100 | 99.59 | 100 | 100 | 100 | 100 | 100 |
| E/Z double bond | 100 | 99.94 | 100 | 100 | 100 | 100 | 99.77 | 100 | 100 | 100 | 100 | 100 |
| tiny (6 or fewer heavy atoms) | 100 | 97.38 | 100 | 100 | 100 | 100 | 99.11 | 100 | 100 | 100 | 100 | 100 |
| macrocycle (ring of 8+) | 100 | 100 | 100 | 100 | 100 | 100 | 99.92 | 100 | 100 | 100 | 100 | 100 |
| metal or other element | 100 | 91.70 | 100 | 100 | 100 | 100 | 99.87 | 100 | 100 | 100 | 100 | 100 |
| cumulene / allene / ketene | 100 | 100 | 100 | 100 | 100 | 100 | 99.84 | 100 | — | 100 | 100 | 100 |
| salt or mixture | 100 | 99.45 | 100 | 100 | 100 | 100 | 99.79 | 100 | 100 | 100 | 100 | 100 |
| radical or carbene | 100 | 93.44 | 100 | 100 | 100 | 100 | 99.85 | — | — | 100 | 100 | 100 |
| isotope label | 100 | 100 | 100 | 100 | 100 | 100 | 100 | 100 | — | 100 | — | 100 |

### lonepairs: disagreements by kind

- **unusual electron count**: 385 — e.g. `CCN(CC)c1ccc2nc3cc(C)c(N)cc3[o+]c2c1.CCN(CC)c1ccc2nc3cc(C)c(N)cc3[o+]c2c1.[Cl-].[Cl-].[Cl-].[Cl-].[Zn+2] atom 46 (Zn): app None, count 0`; `CN(C)c1ccc(/N=N/c2ccc(S(=O)(=O)[O-])cc2)cc1.[Na+] atom 21 (Na): app None, count 0`; `O=S(=O)(O)c1ccc(N=Nc2ccc(Nc3ccccc3)cc2)cc1.[Na] atom 25 (Na): app None, count 0`

### oxidation_kekule: disagreements by kind

- **C(arom) 3 σ-bonds, next to a π bond**: 13,091 — e.g. `Nc1ccc(S(=O)(=O)Nc2ccccn2)cc1 atom 9 (C): app 2.5, RDKit 2`; `Nc1ccc(S(=O)(=O)Nc2ccccn2)cc1 atom 13 (C): app 0.5, RDKit 1`; `c1cncc([C@@H]2CCCCN2)c1 atom 1 (C): app 0.5, RDKit 0`
- **N(arom) 2 σ-bonds, next to a π bond**: 333 — e.g. `c1ccnnc1 atom 3 (N): app -1.5, RDKit -2`; `c1ccnnc1 atom 4 (N): app -1.5, RDKit -2`; `Nc1cccnn1 atom 5 (N): app -1.5, RDKit -1`
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
| index | 250 | 4,314 | 27 (0.62 %) | 20 | 0 (0.00 %) | 14 (0.32 %) | 58 (1.34 %) | 0 |
| cip | 250 | 147 | 0 (0.00 %) | 0 | 0 (0.00 %) | 0 (0.00 %) | 3 (2.04 %) | 0 |
| hybridization | 250 | 4,248 | 28 (0.65 %) | 24 | 0 (0.00 %) | 18 (0.42 %) | 135 (3.18 %) | 0 |
| formal | 250 | 89 | 0 (0.00 %) | 0 | 0 (0.00 %) | 0 (0.00 %) | 2 (2.25 %) | 0 |
| partial | 250 | 3,958 | 42 (1.05 %) | 35 | 14 (0.35 %) | 23 (0.58 %) | 103 (2.60 %) | 0 |
| oxidation | 250 | 4,311 | 30 (0.69 %) | 26 | 0 (0.00 %) | 7 (0.16 %) | 39 (0.90 %) | 0 |
| lonepairs | 250 | 4,270 | 15 (0.35 %) | 15 | 5 (0.12 %) | 5 (0.12 %) | 35 (0.82 %) | 0 |
| pka | 250 | 339 | 1 (0.29 %) | 1 | 4 (1.18 %) | 0 (0.00 %) | 5 (1.47 %) | 30 |
| steric | 250 | 825 | 167 (16.83 %) | 105 | 142 (17.21 %) | 4 (0.48 %) | 31 (3.76 %) | 194 |
| electronic | 250 | 857 | 115 (11.83 %) | 81 | 94 (10.97 %) | 8 (0.93 %) | 23 (2.68 %) | 134 |
| taft | 250 | 516 | 80 (13.42 %) | 69 | 81 (15.70 %) | 4 (0.78 %) | 19 (3.68 %) | 43 |

By kind of molecule (all label modes together):

| Kind of molecule | Labels drawn | Left out | Edge nearer another atom | Inside a ring |
|---|---:|---:|---:|---:|
| 3+ aromatic rings | 10,267 | 220 (2.10 %) | 50 (0.49 %) | 245 (2.39 %) |
| heteroaromatic | 9,113 | 153 (1.65 %) | 36 (0.40 %) | 132 (1.45 %) |
| stereocentres | 6,024 | 91 (1.49 %) | 41 (0.68 %) | 80 (1.33 %) |
| large (40+ heavy atoms) | 5,474 | 113 (2.02 %) | 43 (0.79 %) | 166 (3.03 %) |
| hypervalent S / P / halogen | 4,322 | 192 (4.25 %) | 21 (0.49 %) | 139 (3.22 %) |
| small ring (3-4) | 3,304 | 69 (2.05 %) | 16 (0.48 %) | 75 (2.27 %) |
| ion | 3,292 | 71 (2.11 %) | 4 (0.12 %) | 110 (3.34 %) |
| E/Z double bond | 2,130 | 28 (1.30 %) | 6 (0.28 %) | 18 (0.85 %) |
| zwitterion or dipole | 1,765 | 43 (2.38 %) | 5 (0.28 %) | 56 (3.17 %) |
| metal or other element | 1,679 | 73 (4.17 %) | 10 (0.60 %) | 94 (5.60 %) |
| macrocycle (ring of 8+) | 1,550 | 11 (0.70 %) | 1 (0.06 %) | 93 (6.00 %) |
| salt or mixture | 1,500 | 45 (2.91 %) | 2 (0.13 %) | 40 (2.67 %) |
| cumulene / allene / ketene | 1,147 | 26 (2.22 %) | 1 (0.09 %) | 18 (1.57 %) |
| ordinary organic | 1,067 | 17 (1.57 %) | 0 (0.00 %) | 16 (1.50 %) |
| tiny (6 or fewer heavy atoms) | 1,015 | 37 (3.52 %) | 3 (0.30 %) | 5 (0.49 %) |
| radical or carbene | 962 | 49 (4.85 %) | 5 (0.52 %) | 61 (6.34 %) |
| isotope label | 273 | 24 (8.08 %) | 3 (1.10 %) | 2 (0.73 %) |

Examples of labels nearer another atom: `C[N+](C)(C)c1cccc(O)c1.[Cl-]: '−0.30' for atom 1 sits nearer atom 3 (0.81 vs 0.71 bond lengths)`; `CCCCCCCCN(CC[N+]1(C)CCCCC1)c1ccccc1: '−0.32' for atom 11 sits nearer atom 13 (0.86 vs 0.5 bond lengths)`; `C[N+](C)(C)CCc1cccnc1.[Br-]: '−0.33' for atom 1 sits nearer atom 3 (0.81 vs 0.71 bond lengths)`

### Bond angles: the app's 3D model against RDKit (ETKDG + MMFF94)

1,125 angles between heavy atoms in 55 molecules: median difference 0.19°, 95 % within 4.13°, largest 13.74°; 103 differ by more than 3°, 30 by more than 5°.

Largest: `COc1cc(CC2=C(N)S[C@@H]3C(C(=O)Cc4ccccc4)C(=O)N3C2C(=O)OC(c2ccccc2)c2ccccc2)ccc1O` atoms (1, 2, 43): app 128.8°, RDKit 115.1°; `COc1cc(CC2=C(N)S[C@@H]3[C@H](C(=O)Cc4ccccc4)C(=O)N3C2C(=O)O)ccc1O` atoms (1, 2, 30): app 128.7°, RDKit 115.0°; `COc1cc(CC2=C(N)S[C@@H]3C(C(=O)Cc4ccccc4)C(=O)N3C2C(=O)OC(c2ccccc2)c2ccccc2)ccc1O` atoms (1, 2, 3): app 112.7°, RDKit 126.2°; `COc1cc(CC2=C(N)S[C@@H]3[C@H](C(=O)Cc4ccccc4)C(=O)N3C2C(=O)O)ccc1O` atoms (1, 2, 3): app 112.8°, RDKit 126.2°; `COc1ccc(N=Nc2c(S(=O)(=O)O)cc3cc(S(=O)(=O)O)c(N=Nc4ccc(OC)cc4)c(O)c3c2O)cc1` atoms (27, 28, 29): app 114.7°, RDKit 126.1°.

The same angles compared at each centre in order of size (a group pointing the other way in one model swaps two angles without either being wrong): median 0.11°, 95 % within 2.88°, largest 9.5°; 54 differ by more than 3°, 12 by more than 5°. Largest: `CC(C)CCN1C=C/C(=C\c2cc[n+](CCC(C)C)c3ccccc23)c2ccccc21.[I-]` atom 14: app 109.4°, RDKit 118.9°; `COc1cc(CC2=C(N)S[C@@H]3C(C(=O)Cc4ccccc4)C(=O)N3C2C(=O)OC(c2ccccc2)c2ccccc2)ccc1O` atom 9: app 105.2°, RDKit 96.5°; `COc1cc(CC2=C(N)S[C@@H]3C(C(=O)Cc4ccccc4)C(=O)N3C2C(=O)OC(c2ccccc2)c2ccccc2)ccc1O` atom 10: app 110.8°, RDKit 103.0°; `CC(C)CCN1C=C/C(=C\c2cc[n+](CCC(C)C)c3ccccc23)c2ccccc21.[I-]` atom 16: app 109.7°, RDKit 103.2°.

## Textbook list (`tools/label_textbook.json`)

149 of 149 hand-checked values agree.

## Steric class along crowding series

- aldehydes: ok — `CC=O` open (2), `CCC=O` moderately hindered (3), `CC(C)C=O` hindered (4), `CC(C)(C)C=O` hindered (4)
- alkyl bromides (SN2 order): ok — `CBr` open (1), `CCBr` open (2), `CC(C)Br` moderately hindered (3), `CC(C)(C)Br` hindered (4)
- amines: ok — `CCCCN` open (2), `CCC(C)N` moderately hindered (3), `CC(C)(C)N` hindered (4)
- ketones: ok — `CC(C)=O` moderately hindered (3), `CCC(C)=O` hindered (4), `CC(C)C(C)=O` very hindered (5), `CC(C)(C)C(C)=O` very hindered (5)
