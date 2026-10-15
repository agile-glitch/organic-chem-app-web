# AGILES1 — format specification (for reading AGILES codes without the app)

Audience: a reader (a person, a program, or an AI model in a new conversation) given this document and one or more AGILES
codes, who must report everything the code says — the molecule, its identifiers, and every atom label exactly as the
AgileCHEM Molecule tab displays it. Follow this document literally; where it says a value is copied as text, report the
text verbatim (including its Unicode minus sign U+2212 "−", superscripts "²" "³", and abbreviations).

Format version described: **AGILES 1** (code prefix `AGILES1:`), as written by AgileCHEM 1.1.0 and later.
Older codes with prefix `ACM1:` are the same format under its earlier name (see §10).

---

## 1. Recognising a code

An AGILES code is one line of text: the prefix `AGILES1:` (or `ACM1:`) followed by characters from `A–Z a–z 0–9 - _`
(base64url, no padding). Typical length 700–2,500 characters. Whitespace around it is not part of it. A code cut short
cannot be read (decompression fails): ask for the whole line.

## 2. Unpacking (needs code execution)

The part after the first `:` is **base64url** (RFC 4648 §5, `-` and `_`, padding `=` removed) of the record compressed
with **raw DEFLATE** (RFC 1951, no zlib or gzip header). Unpacked, it is UTF-8 text: one SDF record (§3).

You cannot unpack DEFLATE by reading; run code. In Python:

```python
import base64, zlib
body = code.split(':', 1)[1]
sdf = zlib.decompress(base64.urlsafe_b64decode(body + '=' * (-len(body) % 4)), -15).decode('utf-8')
```

In JavaScript (browser or Node 18+): `new Response(new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate-raw'))).text()`.
The complete reader in §12 does the unpacking and all the parsing below, and prints a table; prefer running it.

## 3. The record

An ordinary SDF record (any SDF reader, e.g. RDKit `Chem.MolFromMolBlock`, reads the molecule and keeps the data items):

1. a **MOL block** (V2000; V3000 only for molecules over 999 atoms) ending with the line `M  END`;
2. **data items**, each: a header line `>  <NAME>`, then zero or more value lines, then one blank line;
3. the line `$$$$`.

Item names all start with `AGILES.` (older codes: `ACM.`, same meaning). A reader must ignore items it does not know
(later versions of format 1 may add items).

### 3.1 The MOL block (V2000)

| Line | Content |
|---|---|
| 1 | title: the molecule's formula as shown in the tab (e.g. `C₉H₈O₄`), or its name |
| 2 | program line (ignore) |
| 3 | comment (ignore) |
| 4 | counts line: columns 1–3 = number of atoms *n*, 4–6 = number of bonds *b* |
| 5 … 4+n | atoms: x (cols 1–10), y (11–20), z (21–30, always 0), element symbol (32–34) |
| next b lines | bonds: first atom (cols 1–3), second atom (4–6), order (7–9: 1 single, 2 double, 3 triple, 9 dative; never 4, aromatic rings are written Kekulé), stereo (10–12: 1 wedge, 6 hash, 4 either/wavy single bond; on a double bond 3 = "either": its E/Z is not given, whatever the drawn geometry; 0 none). A wedge or hash has its **narrow end at the first atom** (the stereocentre), as in every MOL file. A double bond without flag 3 has the E/Z its drawn geometry shows |
| then | property lines: `M  CHG` (count, then pairs: atom number, charge), `M  RAD` (radicals: 2 = one unpaired electron), `M  ISO` (isotope mass), until `M  END` |

- Coordinates are the **2D drawing exactly as the tab lays it out**, in RDKit units (a bond is 1.5 units), **y pointing
  up** (screen y is −y). Aromatic rings are written with alternating single/double bonds (Kekulé form); the tab draws
  them that way.
- Atom numbers 1…n are the **heavy atoms** (every atom written in the block; an explicit isotopic H such as D is one of
  them). Implicit H are not in the block: see §4.
- The structure in the block is authoritative; the identifiers (§5) describe it.

### 3.2 Data items, in the order written

| Item | Value | Present |
|---|---|---|
| `AGILES.format` | `AGILES 1` (older: `AgileCHEM molecule code 1`). A number above 1 means a newer format: stop and say so | always |
| `AGILES.app` | app version and date that wrote the code, e.g. `1.1.0 (2026-10-15)` | usually |
| `AGILES.atoms` | two integers: *nHeavy nH* (heavy atoms = the MOL block's atoms; implicit H added after them) | always |
| `AGILES.h_parents` | *nH* integers in one line: the parent atom number of H atom nHeavy+1, nHeavy+2, … | when nH > 0 (since 2026-10-15) |
| `AGILES.name` | common name | when known |
| `AGILES.iupac` | IUPAC name written by the app's namer | when the namer can name it |
| `AGILES.iupac_checked` | `yes` = the name was read back to exactly this structure; `no` = it could not be checked (the app shows it tagged "unchecked") | with `AGILES.iupac` |
| `AGILES.smiles` | the SMILES as the molecule was given (typed, or the library's for a named compound; else the app's own writer) | when known |
| `AGILES.smiles_rdkit` | RDKit canonical isomeric SMILES: how the tab reads the molecule | when known |
| `AGILES.inchi` | standard InChI | when computable |
| `AGILES.inchikey` | standard InChIKey | when computable |
| `AGILES.cas` | CAS Registry Number from the app's compound library | only for library compounds |
| `AGILES.formula` | molecular formula, Hill order, plain digits, net charge at the end (`C9H8O4`, `C2H3O2-`) | when computable |
| `AGILES.isotopic_formula` | the formula with its isotope labels, as the tab writes it, with subscript digits (`CHD₃O`: D, T for ²H, ³H; other isotopes after their element) | only for an isotope-labelled molecule |
| `AGILES.mw` | molecular weight as the tab shows it: `180.16 g/mol` (CIAAW standard atomic weights, 2 decimals) | when computable |
| `AGILES.monoisotopic` | monoisotopic mass as the tab shows it: `180.0423 Da` (most abundant isotopes, 4 decimals; for an isotope-labelled molecule, the exact mass of that isotopologue) | when computable |
| `AGILES.labels.MODE` | the labels of one Labels mode (§6), one line per labelled atom | per mode |
| `AGILES.labels.angles` | bond angles (§6.12) | when a 3D model was built |
| `AGILES.note` | free text (e.g. why there are no bond angles); may occur several times | optional |

**Salts and mixtures.** The tab shows one component at a time (the largest first); a code describes the component that
was shown, so sodium acetate's code is the acetate ion (`formula` `C2H3O2-`).

**Missing values.** An identifier that is unknown is simply absent (no item): e.g. no IUPAC name when the namer cannot
name the molecule, no CAS for a compound outside the library. Absence is not an error; report "not given".

## 4. Atom numbering and names

- Atoms **1 … nHeavy** = the MOL block's atoms, in its order.
- Atoms **nHeavy+1 … nHeavy+nH** = the implicit hydrogens, one per H, **grouped by parent in ascending parent order**
  (all H of atom 1 first, then those of atom 2, …), as RDKit `AddHs` numbers them. Their parents: `AGILES.h_parents`.
  (Codes without that item: give each heavy atom H = its usual valence − its bond-order sum, adjusted for charge and
  radicals; the reader in §12 does this and checks the total equals nH.)
- **Atom names** (as the tab shows them in "atom number" mode and in its panels): element symbol + a counter per
  element, in atom-number order over all atoms: `C1, C2, O1, …`; the H continue the H count (`H1, H2, …` after any H in
  the block). Use these names when reporting.

## 5. Identifiers

Report them as given. `smiles` and `smiles_rdkit` are the same molecule written two ways (the first as it was entered,
the second as RDKit canonicalises it). The InChIKey is a checksum of the structure; the tab warns when a pasted code's
InChIKey is not its structure's.

## 6. The labels (Labels menu of the Molecule tab)

Each `AGILES.labels.MODE` item lists one line per labelled atom: `<atom number> <label text>` (the first space separates
them; the label text may contain spaces). An atom with no line has **no label** in that mode (blank in the tab). An item
that is present with no lines means the mode labels nothing on this molecule (e.g. no stereocentres). An item that is
absent means the code does not store that mode (the tab then works it out itself). Labels of H atoms (numbers above
nHeavy) are shown only when "Show H" is on.

Number formats used in label text: a negative number uses **U+2212 "−"**; signed values carry "+" when positive;
zero is written without sign.

### 6.1 `cip` — "R/S and E/Z" mode
`R`, `S` (CIP descriptors, RDKit's new CIP labeller), `r`, `s` (pseudoasymmetric centre), `?` (a possible stereocentre,
or a ring cis/trans centre, whose configuration the input leaves open). Displayed as the bare letter (`S`, not `(S)`). Only atoms are labelled here; E/Z of double bonds
is not a stored label — read it from `smiles_rdkit` (`/` `\`) or the MOL block geometry if needed.

### 6.2 `hybridization`
`sp`, `sp²`, `sp³`, `sp³d`, `sp³d²` (superscript digits). From the steric number (σ bonds + lone pairs; VSEPR); an atom
whose lone pair is conjugated with a π system or in an aromatic ring is `sp²` (amide N, phenol O, ester O, pyrrole N);
S, Se and P(III) on an aromatic ring stay `sp³`. H atoms and atoms with no simple count (metals) have no label.

### 6.3 `formal` — formal charge
`+`, `−`, `2+`, `2−`, … ; neutral atoms have no label.

### 6.4 `partial` — partial charge
Gasteiger–Marsili charge (empirical model: RDKit parameters, 12 iterations), 2 decimals, signed: `+0.31`, `−0.43`,
`0.00`. H atoms included. A model estimate, not a measurement.

### 6.5 `oxidation` — oxidation state
IUPAC-style: bonding electrons go to the more electronegative atom (Pauling values); bonds between like atoms split.
Integers `+3`, `−2`, `0`. An atom with an aromatic bond to another element gets the **average over the Kekulé
structures**, which can be fractional: `+0.5`, `−1.5` (up to 2 decimals). A dative bond (L→M) adds nothing (the ligand
keeps its pair). No label = not computable (an element without an electronegativity value).

### 6.6 `lonepairs`
An integer: lone pairs for the Lewis structure as drawn (`0`, `1`, `2`, `3`). A donor atom's pair used in a dative bond
counts as bonding. No label for H, metals and atoms with an unusual electron count. The "electrons" mode of the tab draws
these numbers as pairs of dots around the atom; it has no item of its own.

### 6.7 `pka`
On each acidic site's atom (the H-bearing O, N or S; for a C–H acid, the carbon): `pKa 4.76`. On each basic site's
atom (a neutral N, O, P or S with a lone pair): `pKaH 10.7`. Both on one atom: `pKa 9.1 / pKaH 2.3`. Water, 25 °C.
The number is a measured value (IUPAC Digitized pKa Dataset) when the molecule has one, else a Hammett estimate, else a
typical value for the group (the tab's hover card says which; the code stores only the number). Acid sites above 20 are
shown only for C–H next to a carbonyl, nitrile or sulfone.

### 6.8 `steric`
At each functional group's reacting atom: `open` (crowding score ≤ 2), `mod. hindered` (3, "moderately hindered"),
`hindered` (4), `v. hindered` (≥ 5, "very hindered"). Score = the atom's heavy neighbours (a doubly bonded O/S not
counted) + the most further heavy neighbours any one of them carries.

### 6.9 `electronic`
At each functional group's reacting atom: `δ +0.31` (its partial charge, as in §6.4), optionally followed by
` Σσ +0.12` when the atom is on or next to a benzene ring: the sum of the meta/para substituents' Hammett σ
(> 0 electron-poor ring, < 0 electron-rich).

### 6.10 `taft`
At each functional group's reacting atom: `σ*/Es`, e.g. `+0.60/−1.01`: Σσ* of the groups on that atom (Taft polar
constant; CH3 = 0; > 0 pulls electrons) / Es of the bulkiest group (Taft steric constant; CH3 = 0; more negative =
bulkier). Either part is `–` (U+2013) when no group has a value: `–/+0.69`. "The groups on that atom" are its heavy
neighbours, not counting a doubly bonded O or S or an atom of the reacting atom's own ring. A group in Taft's table
(methyl, ethyl, n-/iso-/tert-butyl, phenyl, CF3, CH2Cl, OH, OMe, halogens …) gives both values; any other carbon group
gives only an Es estimated from its branching (CH2R as n-propyl −0.36, CHR2 as isopropyl −0.47, CR3 as tert-butyl −1.54)
and no σ*, which is why an NH2 on a CH carbon reads `–/−0.47`.

### 6.11 `element`, `index` (no items)
"element" shows every atom's symbol as a label, carbons included; "atom number" shows the names of §4 (`C1`, `O2` …).
Both follow from the MOL block.

### 6.12 `labels.angles` — bond angles
One line per angle: `<centre> <i> <j> <text>` = the angle i–centre–j, e.g. `11 10 12 120.3°` (degrees, 1 decimal),
measured on the tab's 3D model (lowest-energy conformer found, MMFF94 force field). The tab draws the number inside the
wedge between the bonds centre–i and centre–j. Absent when the code was written without a 3D model (an `AGILES.note`
says so); the tab then measures the angles itself.

## 7. How the tab displays all this

The 2D drawing is the MOL block (§3.1): carbon atoms are drawn as line ends and corners without a symbol (unless
charged, isotopic or radical), other atoms by their symbol with their H written beside them (OH, NH2); implicit H are not
drawn unless "Show H" is on (it is off by default), and then the H labels of §6 are shown too. The selected Labels mode puts each atom's label as small italic text beside the
atom, on the side away from its bonds, never over bonds, atom symbols or other labels; a label with no free place is left
out and the drawing says "N labels left out where there is no free space" (hovering the atom shows the value). Where it
is placed is the app's choice and is not stored. When a pasted code's labels differ from what the current app version
works out, the tab still shows the code's labels and a banner lists the differences.

## 8. Telling molecules apart (isomers, stereoisomers)

A code describes exactly one molecule, itself: the drawing carries the stereo (wedges and hashes for tetrahedral centres,
the drawn geometry of each double bond for E/Z), and the identifiers say the same thing in other notations. Changing one
bond, or one centre from R to S, makes a different molecule with a different InChIKey. To decide whether two codes (or a
code and a structure) are the same molecule, compare in this order:

1. **Full InChIKey** (`AGILES.inchikey`, 27 characters `XXXXXXXXXXXXXX-YYYYYYYYSA-Z`). Equal keys: the same molecule,
   stereo included. Different keys: different molecules (or the same molecule in another protonation state: the last
   letter `N` = neutral as drawn, other letters = protons added or removed).
2. **First block** (14 letters before the first `-`): the constitution (which atoms, which bonds). Different first
   blocks with the same formula = **constitutional isomers** (butan-1-ol and diethyl ether; ortho/meta/para).
3. **Same first block, different second block** (8 letters): same constitution, different stereo or isotopes —
   **stereoisomers**. Then compare the stereo in `smiles_rdkit` (`@`/`@@` per centre, `/` `\` per double bond) or the
   `cip` labels: every centre inverted (R↔S at all of them) = **enantiomers**; some but not all inverted, or a double bond
   E↔Z = **diastereomers** (cis/trans isomers included); a molecule identical to its mirror image despite stereocentres =
   **meso** (shown by `r`/`s` or by equal keys for the "two" forms). The second block `UHFFFAOYSA` means no stereo or
   isotope layer at all.
4. **Specified vs unspecified.** A double bond with bond stereo 3 has no E/Z (it is drawn one way only because a drawing
   must draw it somehow). A centre labelled `?` is drawn without configuration: the code is "configuration not
   given", a different InChIKey from either the R or the S form, and must not be read as either.

Limits of the InChIKey (say so when it matters): standard InChI treats mobile-H tautomers as one molecule (2-pyridone and
2-hydroxypyridine get one key), and does not encode axial chirality (atropisomers such as BINOL, chiral allenes) —
compare `smiles_rdkit` and the drawing for those. Conformers (rotations about single bonds, chair flips) are not
different molecules and are not in the code.

The drawing, not the layout, is what counts: the same molecule drawn rotated, flipped (with its wedges flipped too) or
laid out differently has different coordinates and so a different code text, but the same InChIKey. Each code's drawing
is checked to read back as its own molecule and stereoisomer (`tools/agiles_table.py`: the MOL block alone, read by
RDKit, gives the stored InChIKey).

## 9. Reporting a code

Give: the identifiers (§3.2, "not given" for absent ones); a table of atoms (number, name, element, and each stored
mode's label, blank where none); the bonds if asked; the angles if present; the notes. Answer questions by atom name
(e.g. "C2: oxidation state +3, partial charge +0.31").

## 10. Worked example (aspirin)

Code (one line):

<!-- example:begin -->
```
AGILES1:pVRNbtw2FN7zFG_hxQxqcfhI_VBFLcAZxx0jiQeNmyJeGbKGiQhrJEGS4ThLtwtnmyvkBm1P0AOkd_BJCooa1bJVoBlriHnf43v89H0kZ-Z3N58Wdze3y7ub3wiY5_XBC91A__ADQgCFGcBGRxiG8AtnjJn1DqcilMKsFJQLJg1ilDHGYA7_RXF_tCxIWRjy9v3UlVuyAKMYhnJUy3IbLYyGckuWe1rY9vsCSN1QylEt38DCaRBK9lQtPYuDlPubHdrakcMpw21Z-t19ipZ_T_opWu6xCOoNWb7h1vW_o8eO_jcLAvD22yQcQLTfXeL2FRfA6xMPwO_bfICgrwQAsq9IgLCvhICsr3htYivIALFPEJBv2kwiusorgOfHByQC-GH_x6OXz0_ou6Jax01EbApI7hfjsowIUqQMJpxx30HmoDcd9jTFuo4ICpCD-fSsjCuVN6YG5uNDANIYEIO-PF6riMR1qSudDyr6soyTiHBnopo0zovrrPhwPT1X-cdCJxAnevW4_SxJVXKhVhHJi0G1XutM1RGZzyd7y-lyjnvzeTvQTow0n1WrC930SxJMzDPWr_Mk1RE5yufp0R6ezObhQi7dWYKOP0E2ReFIx3OEwx3XCSZyGk4Qp8hnKXe8xS4uxO5ksYu4i3z6mPdCXUfk2cnp8evTV29-fvvs7clPzpvF4eHh_vL0ZN85HqxI4joiHnMC6XDy8JwvszgiVtugtr6KCEpG0Yf3s3WRDYtFXui6aIpSJ7aNuVzAQTzoyuJzldU00WU0Np9en1d6pT_GjS5ycyHq8q8_CDfhdyJscG3wbPBtCGyQNoQ2IOsidrEjwo5pREB7zbNRbWVcNdrUEL5jlAnCTRRIBNzdfmaUe8S1yBXEMzV0iW9nGCdBh3wiexR2yDVSzQJmtBpW14jd0GL3BlcSdLuXo7cB_gYELQgIyhb4BMMOcLYBlp2N-Ss-6FW_73e3n1uDwrrj1ho3vtCaQusIrR20XtAYsSZE54B38nmrHVvh2KrGVjK2erEVi61SbGXimMasyFUZ66r9s2CEAyMCjDijjBEfGAmAEQmMhEZHJ4YZLa0QPnqyF7HhK1_EYHfbIEGFHGuuG1W1NxzWxYpCqvOVqtSK8Ac5PmwY4VKZSpqqyC3f1z_7i2WhQMPSYRe-fvn717Zj9Pia-F0T2TvJ2My22cvks5k5GsqGW5oXjYrIcQHnRb6COH-fqfp7uKp006gcrnSTFpcNxCAOjBGVwaRJFTTxOaxVXF9WqoYmVWu4SlVuECTFSoGuoU6Lq3xKCdnZ2dkh_wA
```

What `python agiles_read.py` prints for it:

```
C₉H₈O₄   (AGILES 1, written by 1.1.0 (2026-10-15))
  name           aspirin
  iupac          2-(ethanoyloxy)benzoic acid
  iupac_checked  no
  smiles         CC(=O)OC1=CC=CC=C1C(=O)O
  smiles_rdkit   CC(=O)Oc1ccccc1C(=O)O
  inchi          InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)
  inchikey       BSYNRYMUTXBXSQ-UHFFFAOYSA-N
  cas            50-78-2
  formula        C9H8O4
  mw             180.16 g/mol
  monoisotopic   180.0423 Da
  13 heavy atoms (with 2D coordinates), then 8 H
  atom | name | x | y | cip | hybridization | formal | partial | oxidation | lonepairs | pka | steric | electronic | taft
  1 | C1 | -2.398 | 3.231 |  | sp³ |  | +0.03 | −3 | 0 | pKa 25 | mod. hindered | δ +0.03 | 
  2 | C2 | -1.099 | 2.481 |  | sp² |  | +0.31 | +3 | 0 |  | mod. hindered | δ +0.31 | +0.00/+0.00
  3 | O1 | 0.200 | 3.231 |  | sp² |  | −0.25 | −2 | 2 |  |  |  | 
  4 | O2 | -1.099 | 0.981 |  | sp² |  | −0.43 | −2 | 2 |  |  |  | 
  5 | C3 | 0.200 | 0.231 |  | sp² |  | +0.14 | +1 | 0 |  |  |  | 
  6 | C4 | 1.499 | 0.981 |  | sp² |  | −0.02 | −1 | 0 |  |  |  | 
  7 | C5 | 2.798 | 0.231 |  | sp² |  | −0.06 | −1 | 0 |  |  |  | 
  8 | C6 | 2.798 | -1.269 |  | sp² |  | −0.06 | −1 | 0 |  |  |  | 
  9 | C7 | 1.499 | -2.019 |  | sp² |  | −0.04 | −1 | 0 |  |  |  | 
  10 | C8 | 0.200 | -1.269 |  | sp² |  | +0.10 | 0 | 0 |  |  |  | 
  11 | C9 | -1.099 | -2.019 |  | sp² |  | +0.34 | +3 | 0 |  | mod. hindered | δ +0.34 Σσ +0.00 | +0.60/−1.01
  12 | O3 | -1.099 | -3.519 |  | sp² |  | −0.25 | −2 | 2 |  |  |  | 
  13 | O4 | -2.398 | -1.269 |  | sp² |  | −0.48 | −2 | 2 | pKa 3.38 |  |  | 
  14 | H1 | on C1 |  |  |  |  | +0.03 | +1 |  |  |  |  | 
  15 | H2 | on C1 |  |  |  |  | +0.03 | +1 |  |  |  |  | 
  16 | H3 | on C1 |  |  |  |  | +0.03 | +1 |  |  |  |  | 
  17 | H4 | on C4 |  |  |  |  | +0.07 | +1 |  |  |  |  | 
  18 | H5 | on C5 |  |  |  |  | +0.06 | +1 |  |  |  |  | 
  19 | H6 | on C6 |  |  |  |  | +0.06 | +1 |  |  |  |  | 
  20 | H7 | on C7 |  |  |  |  | +0.06 | +1 |  |  |  |  | 
  21 | H8 | on O4 |  |  |  |  | +0.30 | +1 |  |  |  |  | 
  bonds (atom-atom order; order 1 single, 2 double, 3 triple, 9 dative; stereo 1 wedge, 6 hash, 4 either (single bond), 3 either (double bond: E/Z not given), narrow end at the first atom):
   1-2 1, 2-3 2, 2-4 1, 4-5 1, 5-6 2, 6-7 1, 7-8 2, 8-9 1, 9-10 2, 5-10 1, 10-11 1, 11-12 2, 11-13 1
  note: No bond angles: written without a 3D model (the tab measures them when the code is shown).
```
<!-- example:end -->



## 11. Versions and compatibility

- The format number is in the prefix (`AGILES1:`) and in `AGILES.format`. Within format 1, items may be **added**; a
  reader ignores items it does not know and treats absent ones as "not given". A change that alters the meaning of an
  existing item or the packing gets a new number (`AGILES2:`), and this document gets a new section.
- History:
  - 2026-10-13 — first version, named "AgileCHEM molecule code", prefix `ACM1:`, items `ACM.*`, format line
    `AgileCHEM molecule code 1`: MOL block, `atoms`, `app`, `labels.*`, `labels.angles`, `note`.
  - 2026-10-15 — renamed **AGILES** (prefix `AGILES1:`, items `AGILES.*`, format line `AGILES 1`; `ACM1:` codes still
    read). Added the identifiers `name`, `iupac`, `iupac_checked`, `smiles`, `smiles_rdkit`, `inchi`, `inchikey`, `cas`,
    `formula`, `isotopic_formula`, `mw`, `monoisotopic` and `h_parents`. Double bonds whose E/Z is not given carry bond
    stereo 3 ("either") in the MOL block, so the drawing reads back as exactly its molecule. For a salt, `smiles` is the
    part of the given SMILES that is the component shown. Steric labels written `mod. hindered` / `v. hindered` (were `moderately hindered` / `very
    hindered`).
- Label values are those of the app version that wrote the code (`AGILES.app`); a later version may work some out
  differently (e.g. 2026-10-15: amino-acid pKa moved from the NH2 to the CO2H oxygen). The tab, given an older code,
  shows the code's labels and a banner listing what it now works out differently.

## 12. Reference reader (Python 3, standard library only)

Run it with a code as the argument (or codes on standard input); `--json` prints the parsed structure. It is the file
`tools/agiles_read.py` of the app's repository; the copy below is kept identical to it by `tools/agiles_doc.py`.

<!-- reader:begin -->
```python
"""Read an AGILES code (AGILES1:... or the older ACM1:...) and print everything in it: the identifiers, every atom with its
element, name (C1, O2 ...), 2D position and its label in every Labels mode, the bonds, the bond angles and the notes.
Standard library only (base64, zlib); docs/AGILES1_format.md describes the format.

  python tools/agiles_read.py "AGILES1:..."        (or pipe codes in, one per line; --json for machine-readable output)
"""
import base64, json, sys, zlib

MODES = ['cip', 'hybridization', 'formal', 'partial', 'oxidation', 'lonepairs', 'pka', 'steric', 'electronic', 'taft']
IDS = ['name', 'iupac', 'iupac_checked', 'smiles', 'smiles_rdkit', 'inchi', 'inchikey', 'cas', 'formula', 'isotopic_formula', 'mw', 'monoisotopic']


def unpack(code):
    code = code.strip()
    prefix, _, body = code.partition(':')
    if prefix not in ('AGILES1', 'ACM1') or not body: raise ValueError('not an AGILES code: it must start with AGILES1: (or ACM1:)')
    return zlib.decompress(base64.urlsafe_b64decode(body + '=' * (-len(body) % 4)), -15).decode('utf-8')   # raw DEFLATE


def parse(sdf):
    sdf = sdf.replace('\r\n', '\n')
    molblock, rest = sdf.split('M  END', 1)
    lines = molblock.split('\n')
    title = lines[0]
    counts = lines[3]
    if 'V3000' in counts: raise ValueError('V3000 MOL block (over 999 atoms): read it with RDKit')
    na, nb = int(counts[0:3]), int(counts[3:6])
    atoms = [{'x': float(l[0:10]), 'y': float(l[10:20]), 'el': l[31:34].strip(), 'chg_field': int(l[36:39] or 0)} for l in lines[4:4 + na]]
    bonds = [{'a': int(l[0:3]), 'b': int(l[3:6]), 'order': int(l[6:9]), 'stereo': int(l[9:12] or 0)} for l in lines[4 + na:4 + na + nb]]
    charges = {}
    for l in lines[4 + na + nb:]:                          # "M  CHG" lines: the formal charges
        if l.startswith('M  CHG'):
            f = l.split()[3:]
            for i in range(0, len(f), 2): charges[int(f[i])] = int(f[i + 1])
    items, cur = {}, None
    for l in rest.split('\n'):
        if l.startswith('$$$$'): break
        if l.startswith('>') and '<' in l:
            cur = l[l.index('<') + 1:l.index('>', l.index('<'))]
            cur = cur.replace('ACM.', 'AGILES.', 1)        # before the rename
            items.setdefault(cur, [])
        elif cur and l.strip(): items[cur].append(l)
        else: cur = None
    one = lambda k: ' '.join(items.get('AGILES.' + k, [])).strip() or None
    nh, nH = map(int, one('atoms').split())
    N = nh + nH
    # the atom list: the MOL block's atoms (1..nh), then one H per implicit H, grouped by parent in ascending parent order
    # (the parents are stored in AGILES.h_parents; for older codes they follow from each heavy atom's valence, as below)
    # each heavy atom's H count from its valence (the smallest allowed valence at or above its bond-order sum, shifted by
    # its charge as RDKit does); the H are then numbered parent by parent. Checked: the counts must add up to nH.
    VAL = {'H': [1], 'B': [3], 'C': [4], 'N': [3, 5], 'O': [2], 'F': [1], 'Si': [4], 'P': [3, 5], 'S': [2, 4, 6], 'Cl': [1, 3, 5, 7],
           'Se': [2, 4, 6], 'Br': [1, 3, 5], 'I': [1, 3, 5, 7]}
    ISO = {1: {'N': 'C', 'O': 'N', 'S': 'P', 'C': 'B', 'B': 'Be'}, -1: {'C': 'N', 'N': 'O', 'O': 'F', 'S': 'Cl', 'B': 'C'}}   # isoelectronic element
    bo = [0] * nh
    for b_ in bonds:
        o = {1: 1, 2: 2, 3: 3, 4: 1.5}.get(b_['order'], 0)
        if b_['order'] == 9: o = 0                         # dative
        bo[b_['a'] - 1] += o; bo[b_['b'] - 1] += o
    rad = {}
    for l in lines[4 + na + nb:]:
        if l.startswith('M  RAD'):
            f = l.split()[3:]
            for i in range(0, len(f), 2): rad[int(f[i])] = {1: 2, 2: 1, 3: 2}.get(int(f[i + 1]), 0)
    hcount = []
    for i, a_ in enumerate(atoms):
        c, el = charges.get(i + 1, 0), a_['el']
        if c and abs(c) == 1 and el in ISO[c]: el = ISO[c][el]
        vals = VAL.get(el)
        if not vals: hcount.append(0); continue
        used = bo[i] + rad.get(i + 1, 0)
        v = next((x for x in vals if x >= used), None)
        hcount.append(max(0, int(round(v - used))) if v is not None else 0)
    parents = [i + 1 for i, k in enumerate(hcount) for _ in range(k)] if sum(hcount) == nH else None
    if one('h_parents'): parents = [int(x) for x in one('h_parents').split()]   # stored since 2026-10-15: no need to work it out
    names, count = [], {}
    els = [a['el'] for a in atoms] + ['H'] * nH
    for e in els: count[e] = count.get(e, 0) + 1; names.append(f'{e}{count[e]}')
    labels = {}
    for m in MODES:
        if 'AGILES.labels.' + m not in items: continue
        arr = [''] * N
        for l in items['AGILES.labels.' + m]:
            k, _, txt = l.strip().partition(' ')
            arr[int(k) - 1] = txt
        labels[m] = arr
    angles = [dict(zip(('centre', 'i', 'j'), map(int, l.split()[:3])), text=l.split(None, 3)[3]) for l in items.get('AGILES.labels.angles', [])]
    return {'title': title, 'format': one('format'), 'app': one('app'), 'ids': {k: one(k) for k in IDS if one(k)},
            'heavy_atoms': nh, 'implicit_H': nH,
            'atoms': [{'n': i + 1, 'name': names[i], 'el': els[i], **({'x': atoms[i]['x'], 'y': atoms[i]['y'], 'charge': charges.get(i + 1, 0)} if i < nh
                                                                    else {'parent': parents[i - nh] if parents else None}),
                       'labels': {m: labels[m][i] for m in labels if labels[m][i]}} for i in range(N)],
            'bonds': bonds, 'angles': angles, 'notes': [' '.join(items[k]) for k in items if k == 'AGILES.note'] or [], 'modes_present': list(labels)}


def show(d):
    print(f"{d['title']}   ({d['format']}, written by {d['app']})")
    for k, v in d['ids'].items(): print(f'  {k:14s} {v}')
    print(f"  {d['heavy_atoms']} heavy atoms (with 2D coordinates), then {d['implicit_H']} H")
    cols = d['modes_present']
    print('  ' + ' | '.join(['atom', 'name', 'x', 'y'] + cols))
    for a in d['atoms']:
        xy = [f"{a['x']:.3f}", f"{a['y']:.3f}"] if 'x' in a else ['on ' + (d['atoms'][a['parent'] - 1]['name'] if a.get('parent') else '?'), '']
        print('  ' + ' | '.join([str(a['n']), a['name']] + xy + [a['labels'].get(m, '') for m in cols]))
    print('  bonds (atom-atom order; order 1 single, 2 double, 3 triple, 9 dative; stereo 1 wedge, 6 hash, 4 either (single bond), 3 either (double bond: E/Z not given), narrow end at the first atom):')
    print('   ' + ', '.join(f"{b['a']}-{b['b']} {b['order']}" + (f" s{b['stereo']}" if b['stereo'] else '') for b in d['bonds']))
    if d['angles']: print('  bond angles (centre: i-centre-j):', ', '.join(f"{x['centre']}: {x['i']}-{x['centre']}-{x['j']} {x['text']}" for x in d['angles']))
    for n in d['notes']: print('  note:', n)


if __name__ == '__main__':
    codes = [a for a in sys.argv[1:] if not a.startswith('--')] or [l for l in sys.stdin.read().split() if ':' in l]
    for c in codes:
        d = parse(unpack(c))
        print(json.dumps(d, ensure_ascii=False, indent=1)) if '--json' in sys.argv else show(d)
        print()
```
<!-- reader:end -->
