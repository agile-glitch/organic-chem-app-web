# Notes for Claude sessions working on this repository

## AGILES (the molecule code)

`js/mol_code.js` writes and reads AGILES codes (`AGILES1:…`, older `ACM1:…`); `js/molview.js` (`codeRecord`, `codeIds`)
decides what goes in them. `docs/AGILES1_format.md` is the specification a reader (or a new Claude conversation) uses to
read a code without the app, and `tools/agiles_read.py` is its reference reader.

Whenever the AGILES code changes — an item added, renamed or reformatted, a label's wording or number format changed
(the label texts come from `MolInfo.labelsFor` in `js/molinfo.js` and `pkaLabels` / `effectLabels` in `js/molview.js`),
the packing or the prefix changed — in the same change:
1. update `docs/AGILES1_format.md` (the item table, the label section, and a line in §10 History);
2. update `tools/agiles_read.py` if reading changes;
3. run `python tools/agiles_doc.py --write` (refreshes the document's reader copy and worked example), then
   `python tools/agiles_doc.py` (must say it is in step), `node tools/test_mol_code.js` and
   `node tools/browser_molcode.js http://localhost:8000/app.html`.
Within format 1 only add items; a change of meaning or packing is a new format number (`AGILES2:`) with its own section.
