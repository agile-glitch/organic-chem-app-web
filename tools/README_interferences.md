# Interferences: what hinders a template, and how the Reactions page uses it

A learned template says "molecules that look like this gave a product that looks like that". It does not say how well the reaction goes at one particular atom,
which of two competing groups reacts, or what else in the molecule gets in the way. Those facts are **interferences**: the electron-withdrawing group next to the
reacting carbon, the branch that crowds it, the alcohol that quenches a Grignard reagent, the sulfur that poisons a palladium catalyst. This layer writes them down as
patterns the program can test, measures what the recorded reactions say about each one, and cross-references them with the templates when the page decides **where**
a reaction happens (regioselectivity) and **which group** reacts (chemoselectivity). Everything is measured on patents that were kept out of everything learned.

## What is in it

The sheet **Interferences** of `reaction_rules.xlsx` (263 rows; `tools/rules_to_js.py` puts them in `data/reaction_rules.js`, which the page reads):

| kind | rows | one row is |
|---|---|---|
| site patterns | 231 | a SMARTS whose first atom is the atom that reacts: "a nitro group para to the reacting ring carbon", "the attacked carbon is neopentyl", "the nitrogen was an amide N", "a peri substituent" … |
| advisories | 32 | a reagent and a group that is present in the molecules: "a Grignard + an O–H" (quenched), "LiAlH4 + a nitro group" (also reduced), "TFA + a Boc group" (cleaved), "Pd + a thioether" (poisoned) … |

The families of site patterns: ring substituents in the ortho / meta / para position with their Hammett constants (32 groups × 3), the kind of ring (next to a
pyridine-type N, between two N, azinium, thiophene / furan / pyrrole α and β, fused α, peri), ortho sterics, aromatic nitrogens (azole N1/N2 sterics), the carbon that is attacked
(methyl / primary / secondary / tertiary, benzylic, allylic, α to a carbonyl, neopentyl, β-branched, ring size, which end of an addition the new group went to),
the carbonyl carbon (aryl, vinyl, α-branched, α-quaternary, α-CF3, two heteroatoms, ring size), N / O / S (aryl, heteroaryl, branched, benzylic, in a ring, next to an
electron-withdrawing group, with an aryl ring that carries an EWG or EDG, 2,6-disubstituted), and retained alkenes.
The advisories: quenchers, incompatible groups, protecting groups, catalyst poisons, competing sites, leaving-group order, medium effects, and two "data support" notes.
Every row says what it means (`Where / when`), what it does to the reaction (`Effect`), why (`Why`) and where that comes from (`Source`).

## How a template is read: marking the atoms it changes

A learned template (rdchiral) writes the atoms it changes with a strict specification (`[N;H0;D3;+0:4]`) and the atoms around them that stay the same plainly (`[c:2]`).
`tools/int_marks.py` finds the changed atoms (those whose specification or neighbours differ between the reactant and the product side) and gives each of them
**the isotope 900 on the product side of the template**. RDKit (and RDKit.js) copy that isotope onto the product atom, so every product says where the reaction
happened: `CC[900CH2][900NH]C(=O)c1ccccc1`. The atoms are marked the same way on the page: `js/interference.js markTemplate` is the same function, and the changed-atom lists
of the 39,110 templates and 21,794 relaxed forms are shipped in `data/learned/i/maps.js` (the page has no reaction-template parser). A site pattern is tested at the marked
atoms; `[!900*]` in a pattern means "an atom the reaction did not change", i.e. a substituent the reacting atom had before. **The patterns are tested on the product**, so they
describe what the reaction leaves in place around the atoms it changed: a group the reaction created or consumed is not seen as a substituent.
Checked in the browser (RDKit.js) against the Python on the final library: the marked templates are identical for 1,500 of 1,500 templates; the flag lists are identical for the 394 hand-written test molecules (`tools/interference_tests.py`: a positive example for each of the 231 patterns, negative examples for the near misses) and for 3,000 real candidate products of the validation patents (`tools/int_check_js.py`); the advisories give the same answers in 68 of 68 cases; whole predictions (`tools/rel_fuse.py --parity`) give the same first answer for 40 of 40 held-out reactions. When the marks were written, 4,000 of 4,000 products came out identical with and without them.

## What is learned, and from what

The recorded reactions are **positive-only**: a record says what happened, never what failed to happen, so "this group hinders the reaction" cannot be read off directly.
Two things can be learned honestly:

1. **Which site won when several were possible** (regioselectivity). For each recorded reaction the starting materials and the reagent are given to the reaction's own template
   with its unchanged leaf atoms stripped (`tools/template_general.py`, the relaxed form the page also runs). Every distinct product is an outcome; the recorded product is the winner.
   154,494 training reactions (about 16 % of all of them) are such contests, and 137,098 of those separate their outcomes by at least one pattern; only those carry information about the patterns. A conditional logit over the outcomes (`tools/int_train.py`) gives each pattern a **global weight** and a **weight per kind of reaction** (11 + "other" families, read from the fragments
   the product consumed), the latter pulled toward the global value, plus three Hammett terms per family: ρ × (sum of the σ constants of the substituents at the ortho, meta and para position).
   The weights are fitted on the training patents, the two penalties are chosen on the validation patents, and the figures below are on validation contests of patents never seen.
2. **Which group reacts when several are present** (chemoselectivity). `tools/int_context.py` gives every fragment F a logistic model of "F is consumed" whose offset is the reagent's recorded
   rate for F (`tools/rel_reactivity.py`, the Reactivity_Order of the sheets found in the data) and whose other terms are the fragments G present in the molecule: β[F,G] < 0 means that G
   being there protects F. The offsets of training reactions are cross-fitted (the rate of a reaction comes from the other four folds of patents), the penalty and whether the context helps at all are chosen
   on one fold of the training patents, and the figures are on the validation patents.

## Results (patents kept out of everything learned)

**1. Which site wins (regioselectivity).** On 8,976 validation contests (those in which the competing sites differ in at least one pattern; 137,098 training contests set the weights):

| | the recorded site picked |
|---|---|
| a random pick among the outcomes | 42.5 % |
| the global weights alone (one weight per pattern, the same in every kind of reaction) | 79.4 % |
| global + per-family weights + Hammett ρ per family | **84.7 %** |

**2. Which group reacts (chemoselectivity).** On 590,119 (reaction, group) cases of validation patents, the log-loss of "this group is consumed" is
0.3674 with the reagent's recorded rate alone and **0.2766** with the other groups present (lower is better).

**3. What it is worth for the first answer.** The features were chosen on the validation patents (out of fold); the test patents were scored once.

| ranking of the candidate products | first answer right (exact) | without stereo | validation, out of fold | wrong site (an isomer) | wrong reaction |
|---|---|---|---|---|---|
| the relevance model alone | 71.51 % | 71.88 % | 73.03 % | 348 | 772 |
| + the interference weights of the site | 71.89 % | 72.28 % | 73.39 % | 313 | 759 |
| + the other groups present | 71.81 % | 72.18 % | 73.29 % | 336 | 749 |
| + both (installed) | 72.00 % | 72.39 % | 73.61 % | 304 | 755 |

The same 11,942 test reactions (2,795 patents) ranked with and without the features: **+0.49 points** of first answers right
(95 % interval +0.22 to +0.76, resampling whole patents; 146 reactions gained, 87 lost; without stereo +0.51).

**4. The advisories.** For 26 of the 32 advisories the recorded reactions say how often the group is really consumed when that reagent is used, against all other
reagents (descriptive figures over all 1,090,000 recorded reactions; nothing is predicted from them). The largest differences:

| advisory | reactions with the reagent and a group | group consumed with the reagent | with other reagents | page would show it in |
|---|---|---|---|---|
| `cleaved-by-fluoride` | 5,118 | 98 % | 29 % | 0.0 % of reactions |
| `cleaved-by-Pd0` | 628 | 90 % | 24 % | 0.0 % of reactions |
| `also-oxidised-by-metal-oxo` | 4,068 | 91 % | 37 % | 0.1 % of reactions |
| `cleaved-by-acid` | 62,021 | 75 % | 23 % | 1.4 % of reactions |
| `cleaved-by-boron-halide` | 9,165 | 70 % | 20 % | 0.3 % of reactions |
| `also-reduced-LiBH4` | 1,387 | 98 % | 50 % | 0.0 % of reactions |
| `also-reduced-LiAlH4` | 9,169 | 91 % | 49 % | 0.1 % of reactions |
| `cleaved-by-base-nucleophile` | 63,418 | 89 % | 48 % | 0.7 % of reactions |
| `also-hydrogenated` | 8,320 | 79 % | 44 % | 0.2 % of reactions |
| `desulfurised-by-Raney-Ni` | 352 | 64 % | 31 % | 0.0 % of reactions |
| `activates-every-OH` | 21,500 | 92 % | 62 % | 0.3 % of reactions |
| `halogen-metal-exchange` | 7,453 | 89 % | 59 % | 0.1 % of reactions |

The rest are properties the records cannot show (conditions, order of addition, the amount of reagent): the row says so.


## What the page does with it

* **Ranking.** Each candidate product of a learned template carries its pattern flags, its family, the learned log-score of its site (`int_lshare`: the log of its share among the sites of the same relaxed
  template) and what it consumed or left judged with the competing groups present (`ctx_consumed`, `ctx_left`, `ctx_both`, `ctx_ll`). These are features of the ranker of `js/relevance.js`
  (fitted on the validation patents, `tools/rel_fuse.py`), next to the template probabilities, the formed / destroyed environments and the reagent's reactivity order.
* **Explanations** (the reading under the arrow, "Interference ·" lines): why this site was preferred over another site of the same template (the patterns that differ, with their learned weights),
  a hindered carbon or ring-nitrogen site (the patterns of kind `hindrance` that the recorded reactions penalise), and the advisories that apply (what the reagent also does to a group that the leading product did not change).
* **Validation tab**: the figures above, stored with the model.

### Examples (what the Reactions page says, reading under the arrow)

Each is a line starting "Interference ·" in the reading; the reaction codes are in `tools/reaction_tests.json` (feature "interferences (site choice, groups present)"), so *Run all tests* checks them.

* **2,4-dichloropyrimidine + morpholine (K2CO3, DMSO)** → the C4 product, 90 % estimated chance of being right.
  *Site choice: another place the same template could react gives Clc1ccnc(N2CCOCC2)n1 (12 % of the ranked weight): the leading product has reacting carbon opposite a pyridine-type ring nitrogen (γ) (+1.1); the other has reacting carbon between two ring nitrogens (−0.7) — weights learned from recorded reactions of the kind "aryl halide".*
* **2,4-difluoronitrobenzene + morpholine (DIPEA, DMSO)** → substitution of the fluorine ortho to the nitro group, 92 %.
  *Site choice: … the leading product has nitro, ortho to the reacting aromatic carbon (+2.9); the other has nitro, para to the reacting aromatic carbon (+1.2).*
* **Boc-glycine benzyl ester + TFA** → the free amine, benzyl ester kept, 88 %.
  *Reagent and groups present: a strong acid (TFA, HCl, HBr, H2SO4, TsOH, MsOH, TfOH, formic acid) + an acid-labile protecting group: F078 tert-butyl protecting group … — removed or opened: Boc and tert-butyl esters / ethers (TFA, HCl), acetals, THP and MOM ethers … benzyl carbamates (Cbz) survive TFA and HCl but go with HBr / AcOH.*
* **4-(benzyloxy)but-1-ene + H2, Pd** → the alkene is reduced, 94 %.
  *Reagent and groups present: H2 with Pd / Pt / Rh / Ni + other reducible groups: F130 benzyl group on O, N or S … — also reduced or cleaved: alkenes and alkynes, benzyl ethers / benzyl esters / Cbz / PMB (hydrogenolysis) …*
* **benzoic acid + 2-aminoethanol (EDC, HOBt)** → the amide, 95 %.
  *Reagent and groups present: an acylating, sulfonylating or alkylating agent + several nucleophilic sites: F003 alcohol (1°), F014 alcohol (on sp3 C) — which site reacts decides the product: amines are acylated before alcohols before phenols …*
* **4-(2-hydroxyethyl)phenol + EtBr (K2CO3, DMF)** → the aryl ethyl ether, 89 %.
  *Site choice: another place the same template could react gives CCOCCc1ccc(O)cc1 (0 % of the ranked weight): the leading product has oxygen on an aromatic ring (+2.6); the other has oxygen on a primary carbinol carbon (+0.6).*
* **Plain Suzuki coupling, N-benzylation of acetanilide, 2-chloro-5-nitropyridine + ethylamine**: no line (nothing the data or the sheet can add).

## Using and extending

* To add or correct an interference, edit the sheet (a new `Key` for a new row; the pattern columns are what is tested) and run `python tools/rules_to_js.py`. Site patterns need weights:
  `python tools/int_reflag.py` (re-reads the stored contests with the new pattern list, a minute or two) and `python tools/int_train.py ...`, then `python tools/int_export.py`.
  `tools/interference_defs.py` is the source the sheet is written from (`python tools/add_interferences.py`); the tests (`python tools/interference_tests.py`) hold a positive and a negative
  example for every pattern.
* The whole chain, in order: `tools/int_contests.py` (contests, 35 min) → `tools/int_train.py` → `tools/int_context.py` → `tools/eval_top1.py --marks` (candidate dumps) →
  `tools/rel_fuse.py --int-model … --ctx-model …` → `tools/int_evidence.py` → `tools/int_export.py` → `tools/add_interferences.py` → `tools/rules_to_js.py` → `tools/int_report.py`.

## Limits

* The data are positive-only. A weight is learned only where a template could have reacted at several sites and the recorded product shows which one won; a pattern that never separates two
  outcomes has no learned weight (its row says so) and stays textbook knowledge.
* The patterns look at the product. What a reaction consumed (an alkene that was hydrogenated, a leaving group) can only be seen through what it left behind; a ring closed onto
  the reacting carbon is told from a ring opened there only for three-membered rings (`sp3-beta-opened`).
* An advisory compares sets of fragments: a group counts as consumed when its count drops, so a second ester of the same molecule is not reported as a bystander when the first one reacts;
  and `epimerised-by-base` is knowledge only (the fragments cannot tell whether an α carbon is a stereocentre, so the page does not test it).
* Reagent tokens are exact strings of the Reagents sheet (`smiles:`), roles, fragments, or an element (`elem:Pd`: some part has a bracket atom of that element); `not:` in front of any of them
  excludes it. A reagent typed in a form the sheet does not list matches by fragment or element only, so an advisory may stay silent for it.
* Temperature, time, amounts and most conditions are not in the learned data, so kinetic against thermodynamic control, equivalents and solvent effects are advisories only.
* Patent records contain errors and multi-step sequences recorded as one step; a share of the "wrong" answers below are of that kind.
* The ceiling is the candidates: 18.8 % of the held-out reactions have no candidate that contains the recorded product, and the weights here can only reorder candidates that exist. Of the reactions that do have one, the wrong answers are 304 isomers of the right product (the site) and 755 different reactions (of 12,000); this layer addresses the first kind and only partly the second. 100 % is not reachable on this data.

## Files

| file | what it is |
|---|---|
| `reaction_rules.xlsx` sheet `Interferences` | the knowledge: patterns, advisories, text, sources, the learned weights and evidence |
| `tools/interference_defs.py` | the same rows as Python (the sheet is written from it): site patterns, advisories, families |
| `tools/int_marks.py`, `tools/int_run.py`, `tools/interferences.py` | marking the changed atoms, running a template with the marks, testing the patterns and advisories |
| `tools/interference_tests.py`, `tools/int_check_js.py` | positive / negative cases; Python against the browser |
| `tools/int_contests.py`, `tools/int_reflag.py`, `tools/int_train.py` | the contests of the recorded reactions, re-flagging them, the weights |
| `tools/int_context.py`, `tools/int_model.py` | the context-aware reactivity table and the loaders of both models |
| `tools/int_evidence.py`, `tools/int_report.py`, `tools/add_interferences.py` | what the recorded reactions say about each advisory; the figures stored with the model; the sheet writer |
| `tools/int_export.py` | `data/learned/i/maps.js`, `weights.js`, `ctx.js` |
| `js/interference.js` | the page side: marking, flags, scores, explanations, advisories |
