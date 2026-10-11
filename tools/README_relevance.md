# The relevance model: how the Reactions page predicts a learned reaction

When no workbook rule (reaction_rules.xlsx) applies, the Reactions page predicts a product from the 1.09 million recorded reactions
in `learned_reactions.xlsx` (US patents + the Open Reaction Database). Until now it ran the templates recorded with the same reagent and
ranked their products by win rate. It now asks a small neural network which of the ~39,000 common recorded templates this reaction
uses, runs the most probable ones, ranks the products, and says how far to trust the first one. It works with or without a reagent.

Everything here is **measured on patents the model never saw** (see "What top-1 means"), and the figure shown on the Validation tab is the
figure stored with the installed model (`data/learned/m/meta.js`), so the number on screen is the number that was measured.

## What "top-1" means here

*Top-1* = how often the first product the page would show is the product recorded for that reaction. Rules of the measurement:

* **Split by patent, never by reaction.** A patent records many near-identical reactions, so a reaction-level split leaks. `rel_common.split_of_patent`:
  md5(patent id) mod 20 -> 0 = TEST, 1 = VALIDATION, otherwise TRAIN. The network learns from TRAIN only. The ranking and the confidence
  are fitted on VALIDATION only. All figures below are on TEST: `12,000` reactions sampled from `2,813` test patents.
* **Input = what the page knows**: the starting materials and the reagent (when there is one). The network also has one input bit per recorded
  solvent (read from the conditions box); the measurement leaves those bits **empty** (`eval_top1.py --solvents none`, the default), which is what a user with an
  empty conditions box gets. Giving the recorded solvent is worth about 0.1 point of template top-1 (67.4 % against 67.3 % of the covered reactions).
* **exact** counts the canonical SMILES with stereo, **flat** the structure without stereo (RDKit's reaction runner does not carry R/S or E/Z
  through a template reliably). Both are reported; the Validation tab leads with *exact*.
* **ceiling** ("generated") = the recorded product is among the candidates at all. No ranking can fix a product that is not there.
* **gradable** = the same figure after leaving out the records the Validation tab calls data errors (the recorded product is a solvent or a
  fragment of a much larger reactant). Only about 1 reaction in 1,000 is dropped by that rule.

## Measured (snapshot; `tools/eval_top1.json` holds the current figures)

| | exact | flat |
|---|---|---|
| first answer right, all reactions | **72.0 %** | **72.4 %** |
| with a reagent given (9,352 reactions) | 73.9 % | 74.3 % |
| no reagent given (2,648 reactions) | 65.4 % | 65.8 % |
| ceiling: recorded product among the candidates | 81.2 % | |
| the page's earlier prediction, same held-out patents | 45.2 % (58.9 % with a reagent, nothing without) | |

How far to trust an answer (the stated chance, calibrated on VALIDATION, checked on TEST):

| the model says | share of all reactions | first answer right |
|---|---|---|
| ≥ 0.95 sure | 32.8 % | 96.7 % |
| ≥ 0.9 sure | 46.6 % | 95.8 % |
| ≥ 0.8 sure | 62.0 % | 93.3 % |
| ≥ 0.7 sure | 68.4 % | 91.1 % |
| ≥ 0.5 sure | 75.7 % | 87.5 % |

The 95 % interval of the exact figure, resampling whole patents (reactions of one patent are alike), is 70.4 to 73.3 %.

The Reactions page shows an answer only when the stated chance reaches the floor of its kind: 50 % with a reagent that has a reactive role in the
Reagents sheet, 70 % with no reagent, 85 % with a reagent that has none (spectator ions, parts the sheet does not list), and never when the starting
materials have fewer than 8 heavy atoms in all. It then answers 71.8 % of reactions and is right 89.3 % of the time; the rest are "Unknown".
The page reads the floors and the size limit from the figures stored with the model (`data/learned/m/meta.js`), so it does what was measured
(`REL_POLICY` in `js/reactions.js` holds only fallbacks).

### Where the figure applies

The learned data are patent reactions. Starting materials of 20 or more heavy atoms are the norm; textbook-size molecules are almost absent
(5 of the 12,000 held-out reactions have fewer than 8 heavy atoms in all). Such inputs are outside what the network learned: it extrapolates and is
confident anyway (acetone + NaCl gave isopropanol at 88 % before the guard; benzene + HBr gives bromobenzene, because a few patents record that). So the
confidence model takes the size of the starting materials as an input (it fixed the mis-calibration of small molecules, test log-loss 0.3621 -> 0.3594) and the page
refuses below the limit. The workbook rules, which come first, answer the textbook reactions they cover.

Held-out figures by size and by kind of reagent (the stated chance is the model's own estimate; "right" is the measured first-answer accuracy):

| heavy atoms in the starting materials | reactions | first answer right | stated | page shows | right when shown |
|---|---|---|---|---|---|
| < 8 | 5 | 40.0 % | 59.2 % | 0.0 % | - |
| 8-11 | 241 | 67.2 % | 68.5 % | 68.5 % | 80.6 % |
| 12-19 | 2,362 | 69.0 % | 70.5 % | 69.4 % | 87.1 % |
| 20-29 | 4,257 | 72.1 % | 73.2 % | 71.6 % | 89.4 % |
| 30 or more | 5,077 | 74.4 % | 75.7 % | 74.1 % | 90.6 % |

| reagent | reactions | first answer right | stated | page shows | right when shown |
|---|---|---|---|---|---|
| a reagent with a reactive role | 8,540 | 76.8 % | 77.8 % | 80.5 % | 89.2 % |
| a reagent without one (spectator ions, parts the Reagents sheet lacks) | 768 | 44.9 % | 48.5 % | 21.2 % | 91.4 % |
| no reagent | 2,634 | 65.8 % | 67.6 % | 60.0 % | 89.5 % |

### What each part is worth (the same held-out reactions each time)

| | ceiling | exact | flat |
|---|---|---|---|
| network only: its own ranking of the products (sum of template probability) | 76.8 % | 67.2 % | 67.6 % |
| + a learned ranking of the products (formed / destroyed environments, template counts, reagent memory) | 76.8 % | 69.5 % | 69.9 % |
| + every template also run relaxed (unchanged leaf atoms stripped) | 81.2 % | 71.2 % | 71.5 % |
| + the reactivity prior (original 103 fragments) | 81.2 % | 71.6 % | 71.9 % |
| + the 66 new fragments in the prior | 81.2 % | 71.6 % | 72.0 % |
| + the 36 sub-fragments in the prior (205 fragments; the baseline of the rows below) | 81.2 % | 71.5 % | 71.9 % |
| + interference weights for the site (regioselectivity) | 81.2 % | 71.9 % | 72.3 % |
| + the other groups present (chemoselectivity) | 81.2 % | 71.8 % | 72.2 % |
| + both = the installed ranker | 81.2 % | 72.0 % | 72.4 % |

### Why 100 % is not reachable on this data

* About 19 % of the reactions have a recorded product that no candidate contains. Looking at them: records that join two transformations
  (a hydrolysis and an alkylation written as one reaction), reactions whose template is not in the learned data (rare chemistry), records with
  an unrecorded reagent or an odd structure. The model cannot generate what the learned templates cannot.
* Of the 12,000 held-out reactions, 2,196 have no candidate with the recorded product, 304 are answered with an isomer of it (the wrong site: the para
  rather than the ortho product, N- rather than O-alkylation) and 755 with a different reaction altogether (a different group reacted, or the record joins
  several steps; reactions that depend on conditions the page is not given: temperature, amounts, order of addition). Several of those records are themselves doubtful.
  The interference layer (below) is aimed at the second and third kind; it moved the first answer by +0.49 points (95 % interval +0.22 to +0.76 points), not by tens.
* The honest way to get close to 100 % is to answer only when sure: the answers rated 90 % or better are right 96 % of the time.

## How it works

1. **Features** (`rel_features.py`, `js/relevance.js`): Morgan bits of every starting material (2048), of the reagent parts (512), one bit
   per distinct reagent part (FNV-1a of its InChIKey, 1024), the Reagents-sheet roles, the Fragments-sheet groups found in the starting
   materials (once / at least twice) and in the reagent, the recorded solvents of the conditions box, and a few counts.
2. **Network** (`rel_train.py`, numpy only): bits -> 1024 -> 512 -> probability for each template seen at least 3 times in TRAIN
   (39,110 classes) plus "none of them". Trained on 982,384 reactions.
3. **Candidates**: the `100` most probable templates are run on the molecules (`runLearned`, the same code the older cascade used: every
   order, every site, the carbon-free and group-transfer parts of the reagent taking part). Each template is also run **relaxed**: with every
   unchanged leaf atom stripped (`template_general.py`, level 3), which lets a reaction use a template that was recorded with different
   neighbours. `data/learned/m/gen3.js` holds the relaxed forms. The probability behind a product is the sum over the templates that give it.
4. **Ranking** (`rel_fuse.py`): a conditional logit over the candidates of one reaction (a bag of four models; the page averages them).
   Features: the mass behind the product, how many templates give it, its size, the environments the reaction made and destroyed (Morgan
   difference, 2 x 2048 sparse weights), how often this exact reagent was recorded with the template (`rc`, from `data/learned/g/`), and which
   functional groups the product consumed against how readily the reagent usually takes them (`rel_reactivity.py`: the recorded reactivity
   order of each reagent, from the Fragments sheet).
   The ranker also gets the **interference features** (`tools/README_interferences.md`): the learned score of the site the template reacted at, relative to the other
   sites of the same template, and what the product consumed or left judged with the other groups present.
5. **Confidence**: a logistic model on the first product's share, margin and the network's own "none of the common templates" probability,
   followed by an isotonic map. It was too sure of itself in the middle range before the map; the reliability table on the Validation tab shows
   stated against measured.

## The 66 new fragments (F104-F169)

`discover_fragments.py` looked for reaction centres the Fragments sheet did not describe (on a sample of TRAIN reactions: 90.6 % of 70,446
centres were already covered by a specific fragment), `new_fragments.py` holds 66 candidates with positive and negative examples,
`fragment_eval.py` scored them against the recorded reactions, `add_fragments.py` appended them to `molecule_data.xlsx` (idempotent; manifest
`tools/data/fragments_new_ids.json`; then `tools/rules_to_js.py` regenerates `data/reaction_rules.js`). Their pKa cells follow the sheet's
"first number wins" rule. The new fragments enter the model through the reactivity prior and the fragment bits of the v2 feature table.
Measured: the v2 network (all 169 fragments + solvent bits) puts the recorded template first for 67.3 % of the covered held-out reactions against 66.3 % for
the first network (103 fragments), and 0.0 points of product top-1 more on the validation patents (73.49 % against 73.46 % flat, cross-validated): the ranker's
reactivity prior already carried what the new fragments tell. The page ships the v2 network.
`learned_reactions.xlsx` was built with the first 103 fragments and was not rebuilt, so its fragment columns do not mention F104+
(`tools/check_learned.py` accounts for that).

## The 36 sub-fragments (F170-F205)

Where one fragment lumps together groups that react differently (an aryl chloride is either a Pd-coupling partner or, next to a ring N, an
SNAr electrophile), `subfragment_scan.py` looks for the split in the recorded reactions. For every fragment and every way its match can be
told apart (the matched atom is aromatic / in a ring of 3-7 atoms / has 0-3 H; a neighbouring atom outside the match is aryl C, CH3, CH2, CH,
quaternary C, vinyl C, carbonyl C, alkynyl C, O, N, charged N, S, a halogen, Si/B/P or a metal) it asks whether the recorded reaction
consumes the group differently with and without that feature **for the same kind of reagent** (the role signature of the Reagents sheet). Rates are
learned on 200,000 TRAIN reactions and the gain in log-likelihood (nats) is scored on 58,000 VALIDATION reactions (13,920 features tried, 1,324 with
enough support; the held-out gains track the training gains with correlation 0.94). `subfragment_select.py` turns the best into SMARTS, checks
that they find the molecules the scan counted, drops those that are the same molecules as a fragment already on the sheet and picks the best few per
parent. From that list 36 sub-fragments were written as chemical classes (`new_subfragments.py`), measured again as written (`subfragment_eval.py`:
every SMARTS matches its positive examples and none of its negative ones), described in the columns of the sheet (`new_subfragments_text.py`: reactive
atom, role with typical reagents, pKa, leaving group, sterics) and added with `add_fragments.py --sub`, each with the evidence in its Notes cell.

Examples (the parent group is consumed in the recorded reaction this often, with the feature / without it, same kind of reagent):
aryl chloride next to a ring N 73 % / 13 %; secondary amine on an azine 11 % / 39 %; primary alkylamine 94 % / 74 %; α-C–H of an aldehyde 93 % / 51 %;
acyclic enone 76 % / 40 %; acyclic trialkylamine (a base) 38 % / 9 %; imide N–H 13 % / 63 %. Summed held-out gain of the 36: 7,764 nats
(0.13 nats per validation reaction).

Measured effect on product top-1: **none**. The ranking with the reactivity table of all 205 fragments scores 71.5 % exact on the test reactions against
71.6 % with 169 (the first 66 new fragments moved it by 0.05 points as well, and the v2 network, which has them as inputs, gained a point of template
top-1 but nothing at product level): the ranker already gets this information from the formed / destroyed environments. The page ships the 205-fragment table (the same within noise; the interference features were fitted with it) and
the 169-fragment network that was measured. What the sub-fragments change is what the app can say: the Molecule page and "How it was read" on the
Reactions page now tell them apart, each with its pKa, role and sterics. `rel_reactivity.py --counts --tag new2` and `rel_recount.py` rebuild the
tables and the evaluation dumps for the longer sheet if a later change makes them worth shipping.

## Interferences

Built: see **`tools/README_interferences.md`** (the sheet `Interferences` of `reaction_rules.xlsx`, the learned weights, the context table, the page lines and the figures).
What a template only sees through its own atoms, and what else decides whether and how it fires (the checklist the layer was built from):

* **electronic**: electron-donating / withdrawing groups (the position matters: ortho, para, meta; α, β), electron-poor rings (azines) against rich ones (pyrrole, furan), conjugation, charge state (ammonium, carboxylate, salts);
* **steric**: branching at the reacting atom, ortho substituents, bulky protecting groups, axial / equatorial attack, bridgeheads, bulky reagents;
* **competing sites**: other groups the same reagent could attack (aldehyde > ketone > ester for hydrides, I > Br > Cl in couplings), intramolecular alternatives;
* **quenching / incompatible groups**: acidic H (OH, NH, CO2H, terminal alkyne) against organometallics and strong bases; oxidants and thioethers; hydrides and esters; acid and Boc / acetals / silyl ethers; H2/Pd and benzyl, nitro, aryl halides;
* **protecting groups**: they block a site, and the conditions can also remove them;
* **coordination**: catalyst poisons (S, free amines, pyridines), Lewis-basic groups that quench Lewis acids, directing and chelating groups;
* **geometry**: ring strain, Baldwin and Bredt limits, the anti-periplanar requirement for E2, facial selectivity from existing stereocentres (the exact-vs-flat gap);
* **conditions**: solvent, water, air, temperature, time, equivalents, order of addition, pH. The learned data records only the solvent;
* **leaving group and reagent strength**: poor leaving groups need activation, a base can be too weak for the pKa;
* **record artefacts** that look like interference: unrecorded reagents, two steps in one record, molecules outside the data (the size guard), template-extraction errors.

The recorded reactions contain only reactions that happened, so an interference can only be learned from **which site or pathway won when several were possible**, or from curated
knowledge; the layer does both, and says in each row which. Two first measurements preceded it: `interference_scan.py` (which groups in the same starting materials change how readily
another group is consumed by the same kind of reagent: competition and activation, confirmed on held-out patents; e.g. when a 1° or 2° amine is consumed an alcohol in the same molecule
is spared, odds ratio 0.12; when an aryl bromide is consumed the aryl chloride is spared, 0.13; an aryl chloride is 8 times as likely to react when a ring N is present and 3.5 times when
a nitro group is) and the sub-fragments above (the interference that sits on the fragment's own atoms).

## Rebuilding (in this order; nothing here changes the committed data until the last step)

```
python tools/rel_conditions.py                                   # recorded solvents of every learned reaction (about 10 min)
python tools/rel_features.py --tag v2 --solvents                 # the model's input bits for every learned reaction (about 10 min)
python tools/rel_reactivity.py --counts --tag new2               # fragment counts of every reaction (all 205 fragments)
python tools/rel_reactivity.py --table --tag new2                # reactivity table from TRAIN patents  (honest measurement)
python tools/rel_reactivity.py --table --all --tag new2          # ... from every reaction              (what the page ships)
python tools/rel_train.py --feat v2 --tag full2 --minc 3 --hidden 1024,512 --epochs 4 --test 10000    # about an hour
python tools/eval_top1.py --model full2 --split val  --n 15000 --topn 100 --gen 3 --marks --dump val_full2m2.pkl     # --marks: the interference flags of every candidate
python tools/eval_top1.py --model full2 --split test --n 12000 --topn 100 --gen 3 --marks --dump test_full2m2.pkl
# the interference weights and the context table (tools/README_interferences.md): int_contests.py (35 min), int_train.py --grid, int_context.py --grid, then
python tools/rel_fuse.py --val val_full2m2.pkl --test test_full2m2.pkl --prior new2 --sources mlp,g3 --int-model v4 --int-feats lshare --ctx-model new2 --ctx-feats consumed,left,both,ll --out data/learned/m/fusion.js --report tools/eval_top1.json
python tools/int_evidence.py --out tools/interference_evidence.json --counts new2
python tools/int_report.py --int-model v4 --ctx-model new2 --runs base=...,site=...,context=...,both=... --paired ana_base.pkl,ana_both.pkl      # the figures of the Validation tab
python tools/rel_reactivity.py --export data/learned/m/prior.js --tag new2
python tools/rel_export.py --model full2 --out data/learned/m --report tools/eval_top1.json
python tools/int_export.py --model full2 --int-model v4 --ctx-model new2
python tools/add_interferences.py --int-model v4 ; python tools/rules_to_js.py
```

`learned_reactions.xlsx` (git-ignored, rebuilt by `python tools/learn_reactions.py`) is the source of every step; its parsed rows are cached in the temp folder
(`organic_chem_*` files are caches; every step says what it wrote). `python tools/rel_features.py --tag v1 --first 103` rebuilds the
original 103-fragment table (the first model, `full1`, used it).

## Files

| | |
|---|---|
| `js/relevance.js` | the browser side: features, network, candidates, ranking, confidence |
| `js/reactions.js` | `relevanceReaction` (the hook in `learnedReaction`, after the exact look-up), the floors, `window.ReactionsRelevance` |
| `js/validation.js` | the Validation tab card (reads `data/learned/m/meta.js`) |
| `data/learned/m/` | `meta` (layout, measured figures), `w1 w2 w3_*` (int8 weights, one scale per output unit), `classes`, `gen3`, `fusion`, `prior` (about 29 MB) |
| `tools/rel_*.py`, `eval_top1.py`, `template_general.py` | features, training, evaluation, ranking, reactivity table, export |
| `tools/subfragment_scan.py`, `subfragment_select.py`, `subfragment_eval.py`, `new_subfragments*.py` | sub-fragment discovery, selection, measurement and descriptions |
| `tools/interference_scan.py`, `rel_recount.py` | which groups hinder or help another group; fragment counts of a longer sheet for an evaluation dump |
| `js/interference.js`, `data/learned/i/` | the interference layer on the page: marking, flags, scores, advisories; the changed-atom lists, learned weights and context table |
| `tools/interference_defs.py`, `interferences.py`, `int_*.py`, `add_interferences.py` | the interference patterns, engine, learning, export and sheet writer (see `tools/README_interferences.md`) |
| `tools/eval_top1.json` | the measured figures |
| `tools/reaction_tests.json` | "relevance model" regression cases for the page's *Run all tests* |

## Checking the browser port

The page recomputes the model's input bits and the whole prediction in JavaScript. Two comparisons with the Python, on held-out reactions:

* **Input vectors and the network** (`tools/rel_check_js.py`, the console snippet is in its docstring): the page must give the same input bits for every
  reaction (also with solvents in the conditions box, `--conds "toluene, water"`) and, up to int8 rounding of the weights, the same most probable templates.
  Measured on 40 reactions: 40 of 40 input vectors identical, first-choice template equal for 40 of 40.
* **Whole predictions** (`python tools/rel_fuse.py ... --parity tools/_parity.json`, then in the console with the Reactions tab open):

```js
const P = await (await fetch('tools/_parity.json')).json(); let same = 0;
for (const p of P) { const o = await ReactionsRelevance.predict(p.subs, p.reagent, ''); if (o && o.candidates[0] && o.candidates[0].flat === p.flat) same++; }
console.log(same, 'of', P.length, 'first answers equal');
```

  Measured: 40 of 40 equal. The page's stated confidence for these (test) reactions is a few points higher on average than the Python's, because the page's
  reagent counts (`data/learned/g/`) and reactivity table (`prior.js`) are built from every reaction and so include them; for a new reaction there is no such gap.
  Delete the json files afterwards.

`ReactionsRelevance.inspect(subs, reagent, n, conditionsText)` returns the input bits and the most probable templates.

## Limits

* `29 MB` of model data are loaded on the first learned prediction (plus the template chunks it needs, up to 15 MB); one prediction then takes about half a second.
* The network learned from 90 % of the patents (TRAIN), the ranking and the confidence from another 5 % (VALIDATION); the last 5 % (TEST) was used for nothing
  but the figures above. The reactivity table and the "recorded with this reagent" counts shipped with the page are built from every reaction, which can only
  help; the figures above use the TRAIN-only versions, as does the shipped network.
* The conditions box contributes only its solvents; temperature, amounts and time are not model inputs.
* Rebuilding the learned data (`tools/learn_reactions.py`) would use all 205 fragments; the model's fragment bits then need `rel_features.py` and a retrain.
