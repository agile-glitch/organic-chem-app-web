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
