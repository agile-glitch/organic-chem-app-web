# Organic chemistry app (website)

The pages AgileCHEM shows inside its window, and that run in any browser (no server): **Sketcher**, **Draw**, **Molecule**, **Retrosynthesis**, **Reactions**, **Simulate**, **Possibilities**, **Validation**.
Open `index.html` (or `python tools/serve.py` and visit http://localhost:8000): the home screen. Each app icon opens the app in a tab of its own,
so several can be open at once (Draw1, Draw2, Molecule1 …); a tab is `app.html?embed=<page>` in a frame (`js/shell.js`, `js/app.js`).
`app.html` on its own is the whole app with its original tab bar, and `index.html#reactions` opens that app on arrival.

What is where:

* `tools/README_atlas.md` — the Fragment Atlas (groups understood in their surroundings, learned from 479,000 atom-mapped patent reactions), the interference judgement, reagent and condition
  suggestions, the learned retrosynthesis, the reaction code format, what was measured and the limits. Start here.
* `tools/README_relevance.md`, `tools/README_interferences.md` — the learned reaction predictor (72.0 % top-1 on held-out patents) and the hand-written interference sheet it builds on.
* `docs/AGILES1_format.md` — how to read an AGILES code without the app (give it to a new Claude conversation with a code); `tools/agiles_read.py` is its reader, `AGILES_codes.xlsx` lists example molecules (the library's 6,987 compounds) with their codes and identifiers; every molecule the Molecule tab can show has a code.
* `js/mol_code.js` — AGILES, the Molecule tab's molecule code: Export → "AGILES file (SDF with all labels)" (the 2D MOL block plus every label mode's labels per atom and the bond angles, as SDF data items any SDF reader keeps) or "AGILES code (copy)" (that record DEFLATE-compressed into one `AGILES1:` line; older `ACM1:` codes still read). Pasting the code into the molecule box, or opening the SDF, shows the molecule with the stored labels; a banner lists any label this version works out differently.
* `reaction_rules.xlsx`, `molecule_data.xlsx` — the workbooks: 102 textbook reactions with mechanisms, 205 fragments, 1,200 reagent parts, conditions (now including `concentrated` / `dilute`).
* `data/` — everything the pages load: `learned/` (recorded reactions, templates, the relevance model, recipes per reaction type), `atlas/` (the Fragment Atlas tables), `rxnkb/`, `sources/`.

Checks (all run the real page in a headless browser or the real code in Node):

```
python tools/serve.py &
node tools/browser_tests.js    http://localhost:8000/app.html  # the Reactions page's own suite: every workbook reaction + feature cases (231)
node tools/browser_features.js http://localhost:8000/app.html  # Atlas, reagent suggestions, learned retrosynthesis, equation balancer, links between tabs (18)
node tools/atlas_parity.js     http://localhost:8000/app.html --judge   # the page against Python: keys of 3,606 fragment instances, 400 probabilities
node tools/test_rxn_code.js && node tools/test_stoich.js   # reaction code, reagent strengths; exact equation balancing
node tools/test_mol_code.js && node tools/browser_molcode.js http://localhost:8000/app.html   # AGILES molecule code: SDF with every label, AGILES1: paste code, round trip
python tools/agiles_doc.py                                 # docs/AGILES1_format.md (how to read an AGILES code) still matches the code
python tools/run_all_checks.py                             # the workbook reactions through RDKit (needs rdkit, openpyxl)
node tools/patent_bench.js run && python tools/patent_bench_grade.py   # 1,000 held-out US patent reactions through the page (tools/patent_bench_report.md)
```
