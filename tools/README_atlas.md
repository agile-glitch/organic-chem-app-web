# The Fragment Atlas, reagent suggestions and the reaction code

What a chemist reasons with is not "an ester" but *this* ester: how crowded it is, what is next to it, its charge, the ring it sits in. The Atlas makes that the unit of knowledge, learns from
recorded reactions what happens to such groups, and uses it for three things the Reactions page does: **judge interferences**, **suggest reagents and conditions**, and **plan routes**.
Everything below is measured on reactions held out of the learning; the figures are those stored in `data/atlas/meta.js` and shown here.

## 1. Contextual fragments: how the program understands a group

A *fragment instance* is one of the 205 named groups of the Fragments sheet at one place in one molecule. Its identity is the group plus the surroundings of its reacting atom at four radii
(`tools/atlas_core.py`, `js/atlas.js`): a Weisfeiler-Lehman refinement of element, aromaticity, branching, H count, charge, ring membership and bond kinds. The same ester in two molecules
shares radius 1 and differs at radius 2 or 3, so the number of distinguishable fragments is not a list somebody wrote. Counted in the mined reactions there are
**1,565,749 contexts that were met often enough to keep** (and about 4.7 million distinct ones before pruning); every other context is still judged, by backing off to the coarser ones
(`p_k = (hits_k + 5 p_(k-1)) / (n_k + 5)`), and says how much evidence it stands on. Alongside the key every instance carries what a chemist looks at: how crowded the reacting atom is
(Taft-style neighbour count), its Gasteiger partial charge, the smallest ring and the angle strain it forces (regular-polygon angle against 109.5 / 120 / 180 degrees), and the pKa of the group
against the reagent's.

## 2. Where the knowledge comes from

`tools/atlas_mine.py` reads the **479,035 atom-mapped US-patent reactions of USPTO-MIT** (W. Jin, C. W. Coley, R. Barzilay, T. Jaakkola, NeurIPS 2017; derived from D. M. Lowe's patent data, CC0).
The atom maps say exactly which atoms changed, so no template and no guess is needed. A molecule none of whose atoms reaches the product is a spectator (reagent, catalyst, base, solvent);
the others are participants, and every fragment instance of a participant is recorded with what the reaction did to it, the reagents on the other side, and four environment hashes:
**8,061,104 instances**. Splits are by molecule scaffold or, the strictest, by molecule (a reaction is test if *any* of its molecules is a test molecule, so a held-out molecule never occurs in training).

Two questions are counted for every context (`tools/atlas_keys.py`, ten table families: the group alone, with the kind of reagent, with the surroundings at radius 1/2/3, with and without the reagent, with a named reagent):

* **Transformed**: of the groups met in reacting molecules, how many did the reaction transform (this includes the groups the reaction is *for*).
* **Also transformed (the interference question)**: the reaction takes place somewhere else in the molecule; is *this* group also transformed? The changed atoms of a molecule form *sites*
  (connected pieces); a group touching no site is spared, a group on one of two or more sites was transformed too, a group on or next to the only site is the reaction itself and is not counted.
  5,587,744 instances qualify; the base rate is about 1 %. (A first version counted only "transformed" and reported a phenol as "also reacting" in 95 % of amide couplings: it mixed in the cases where the phenol is the reaction site,
  for example alkylations. The two questions are separate for that reason.)

## 3. What was measured (held out by molecule)

**Is the group transformed?** (6,680,737 training / 626,917 validation / 753,450 test instances; counts from training, the small combining model fitted on validation, figures on test.)

| model | log-loss | AUC | AUC within one group type |
|---|---|---|---|
| the group alone | 0.3369 | 0.881 | 0.500 |
| + the kind of reagent | 0.2979 | 0.908 | 0.677 |
| + the group's surroundings (radius 1, 2, 3) | 0.1965 | 0.962 | 0.884 |
| + the named reagents | 0.1847 | 0.966 | 0.897 |
| + steric, charge, ring-angle and pKa descriptors | 0.1838 | 0.966 | 0.897 |

Which of several same-type groups in one molecule reacts (12,513 molecules, one transformed): chance 41 %, the group and reagent alone 41 %, with the surroundings **70 %**,
the descriptors alone 47 %.

**Is the group also transformed when the reaction is elsewhere?** (4,637,234 / 432,368 / 518,142)

| model | log-loss | AUC | AUC within one group type |
|---|---|---|---|
| the group alone | 0.0395 | 0.903 | 0.500 |
| + the kind of reagent | 0.0372 | 0.923 | 0.662 |
| + the group's surroundings (radius 1, 2, 3) | 0.0242 | 0.973 | 0.885 |
| + the named reagents | 0.0230 | 0.976 | 0.898 |
| + steric, charge, ring-angle and pKa descriptors | 0.0229 | 0.976 | 0.895 |

How much the surroundings that were *never seen* matter (log-loss on the part of the test set whose radius-3 / 2 / 1 context never occurred in training):

| the test instance's context | share | reagent kind only | + surroundings and reagents | + descriptors |
|---|---|---|---|---|
| all | 100.0 % | 0.298 | 0.185 | 0.184 |
| seen at radius 3 | 82.4 % | 0.289 | 0.169 | 0.168 |
| seen at radius 2, not 3 | 13.1 % | 0.322 | 0.238 | 0.239 |
| seen at radius 1, not 2 | 3.6 % | 0.428 | 0.352 | 0.352 |
| radius 1 never seen | 0.8 % | 0.240 | 0.169 | 0.139 |

The scaffold split (related molecules stay on one side; `tools/atlas_eval_scaffold.json`, measured on an earlier mining of the same reactions) gives the same picture for "transformed":
log-loss 0.317 (reagent kind) to 0.207 (surroundings) to 0.194 (named reagents), site choice 71 % against 44 % by chance.

Reading the numbers honestly:

* The surroundings do most of the work (log-loss 0.298 to 0.185 for "transformed", 0.037 to 0.023 for "also transformed"). That is the point: a group is understood through where it sits.
* The steric, charge, angle and pKa descriptors add almost nothing where the surroundings were seen, because a radius-2 or 3 environment already determines most of them. They help where the
  surroundings were never seen (radius 1 never seen: log-loss 0.169 to 0.139), which is the case you care about: a fragment that has not been recorded. That is 0.8 % of the test instances.
  A non-linear combiner (gradient boosting, not shipped) reaches 0.100 there, so there is more to gain from the physics in the unfamiliar corner.
* The records are positive-only: they show what was tried and worked, never what failed. "Spared in 40 of 40 recorded cases" is evidence; "no recorded case" is silence, and the page says
  **"no precedent: unknown, not safe"**, never "spared". (LiAlH4 with a nitro group elsewhere reads that way, which is right.)
* There is no temperature, time or equivalents in the patent records, so the rates are over the conditions people used.

## 4. What the pages do with it

| where | what |
|---|---|
| Reactions, under the arrow | **Groups in their surroundings**: each group as target, "may / likely also react", "recorded as spared" or "no precedent", with the estimate, the number of recorded cases and a confidence, and a table of crowding, charge, ring and angle strain, pKa. |
| Reactions, "Suggest reagents" | Textbook rules that fit, and what was recorded for the kinds of reaction the relevance model proposes (reagent lines, solvents, median yield, from `data/learned/c/`, built by `tools/cond_build.py`); each idea screened by the Atlas and ranked by how likely the reagent is to transform the target at all, then by how well it spares the rest. "Use" puts it on the page and runs it. |
| Reactions, right-click the reagent | Strength (37 % HCl, 98 % H2SO4 ...: `js/reagent_conc.js`) and amount (equiv, excess, catalytic). The strength reaches the rules as `concentrated` / `dilute` (workbook sheet Conditions): 10 % H2SO4 at 160 C no longer dehydrates an alcohol, and the page says why. A rule is assumed to have its textbook strength when none is given. |
| Reactions, equation | **Equation and amounts** (`js/stoich.js`): exact integer null-space balancing, charge as a row, the smallest set of small molecules that closes an equation named; scalable table of mmol, mg, equivalents, atom economy. |
| Retrosynthesis | **Disconnections from recorded reactions** (`js/retro_learned.js`): the 15,000 commonest recorded reaction types run backwards on the target, ranked by how often each was recorded, how much simpler the precursors are and how many can be bought; recorded reagents, solvents and Atlas screening for each; "Open this step in Reactions"; "Plan a whole route". |
| Molecule | **Groups in context**, and the reagents recorded most often transforming a group like it. |
| Sketcher, Possibilities | Send a drawing to Reactions, Retrosynthesis or Molecule; "Open in Reactions" now works. |

## 5. The reaction code (version 2)

```
reactants>reagents>products | conc=37%w/w | cond=reflux, ethanol, 2 equiv
```

The first word is an ordinary reaction SMILES (RDKit's `ReactionFromSmarts(line.split()[0])`, USPTO and ORD exports, an electronic notebook can read it). What a reaction SMILES cannot say follows after a
space as `key=value` items: `conc` (the reagent's strength: `37%w/w`, `5%v/v`, `3%w/v`, `2M`, `conc`, `dil`) and `cond` (the conditions box in its standard spelling). The old
`starting materials~reagent~conditions~products` is still read (typed or pasted into the code bar, or passed to `Reactions.load`) and is no longer written or stored: every stored code (the workbook's
102 reactions, the test cases and the 1.09 million recorded reactions) was rewritten as version 2 by `tools/convert_codes_v2.py`, which checks that each one reads back with the same molecules
and conditions. A `>` inside a molecule's `->` (a dative bond, as in some palladium precatalysts) is not a field separator, as in RDKit. `js/rxn_code.js` is the whole definition;
`tools/rxn_code.py` is the same in Python for the tools that build the data.

## 6. Research this follows, and where it differs

* Recommending reagents, solvents and catalysts from recorded reactions: H. Gao et al., *ACS Cent. Sci.* 2018, 4, 1465 (a neural network over about 10 million Reaxys reactions, top-10 match 69.6 %, temperature within 20 C in 60-70 %).
  Here the same idea is done by counting what was recorded, per reaction type, which is explainable and runs in the page; it has no temperature because the patent data have none.
  Wang et al., *Research* 2023 (Parrot): a Transformer condition predictor with open condition datasets and attention to the reaction centre.
* Functional-group compatibility from reaction databases: ChemiRise (arXiv 2108.04682) scores compatibility with a naive-Bayes product over groups, assuming independence. The Atlas conditions on each group's
  surroundings and the named reagent, and separates "also transformed" from "transformed".
* Reactivity from structure: Mayr's N/E/sN parameters (log k = sN (N + E); over 1,300 nucleophiles and 350 electrophiles) are the physical-organic version of this question; they are not in the program
  (the database was not reachable and values are not reproduced from memory). The Atlas learns the same kind of dependence from outcomes instead.
* Equation balancing by the null space of the element (and charge) matrix with exact arithmetic is the standard algebraic method; `js/stoich.js` adds the search for the missing small molecules.
* Template-based retrosynthesis from reversed, popularity-ranked reaction templates is the approach of ASKCOS-style planners; here the templates are the 39,000 learned ones and the ranking adds buyability
  and a simplification term. A proposal is a candidate, not proof: check it in the Reactions tab.

## 7. Rebuilding and checking

```
python tools/atlas_mine.py --src <USPTO-MIT folder with train/valid/test.txt> --out instances.npz          # 14 min, 4 cores
python tools/atlas_eval.py instances.npz --split molecule --target any       --out tools/atlas_eval_molecule.json
python tools/atlas_eval.py instances.npz --split molecule --target bystander --out tools/atlas_eval_bystander_molecule.json
python tools/atlas_build.py instances.npz tools/atlas_eval_molecule.json tools/atlas_eval_bystander_molecule.json   # data/atlas/
python tools/cond_build.py                                                                                  # data/learned/c/
python tools/atlas_parity_ref.py && python tools/atlas_parity_judge_ref.py instances.npz tools/atlas_eval_molecule.json tools/atlas_eval_bystander_molecule.json
python tools/serve.py &                                                                                     # then:
node tools/atlas_parity.js http://localhost:8000/ --judge       # keys and probabilities, Python against the page (3,606 instances; 400 probabilities, max difference 2e-10)
node tools/browser_tests.js http://localhost:8000/              # the Reactions page's own suite (231 cases)
node tools/browser_features.js http://localhost:8000/           # the Atlas, suggestions, retrosynthesis, balancer, links (18 checks)
node tools/test_rxn_code.js && node tools/test_stoich.js        # reaction code, reagent strengths, equation balancing
```

## 8. Limits

* Records are positive-only and carry no temperature, time or equivalents (see section 3). Textbook conditions come from the workbook rules.
* "Transformed" is judged on the atoms of the group's pattern, so a group whose pattern includes ring atoms can read as transformed when the ring is substituted nearby; overlapping group definitions
  (ester / methyl-ethyl ester) are merged per atom on the page.
* An aromatic or Kekule bond-label change counts as a change; it is 0.2 % of changed atoms in a sample of 6,000 reactions.
* The learned retrosynthesis ranks by popularity, so routine steps (protections, hydrolyses) are favoured over clever ones; precursors that are not in the catalogue are called "make", not "buy".
* USPTO-MIT is a subset (479,000) of the 1.09 million reactions behind the relevance model. The Atlas was not used to re-rank the relevance model's first answer; whether it would add to its 72.0 % top-1 is not measured.
