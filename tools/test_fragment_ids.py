"""Test the hash fragment IDs end to end:   python tools/test_fragment_ids.py

1. every ID equals the hash of its SMARTS, and no ID or SMARTS repeats
2. a real reaction (Fischer esterification, R001): the fragments it consumes/forms are found by their ID's SMARTS
3. editing a SMARTS in a throw-away copy of the repo, then '--fix': the old ID disappears everywhere and the new one
   shows up in the reaction rules, the workbook and the learned-reaction data
"""
import json, os, re, shutil, sys, tempfile, gzip, base64, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fragment_id as FI
from rdkit import Chem
fails = []
def check(ok, what):
    print(("PASS  " if ok else "FAIL  ") + what)
    if not ok: fails.append(what)

# 1 ---------------------------------------------------------------------------------------------------------------
rows = FI.sheet_rows()
ids = [str(r["Fragment ID"]) for r in rows]; sm = [str(r["SMARTS pattern"]).strip() for r in rows]
check(all(i == FI.fragment_id(s) for i, s in zip(ids, sm)), f"all {len(ids)} IDs equal the hash of their SMARTS")
check(len(set(ids)) == len(ids) and len(set(sm)) == len(sm), "no ID and no SMARTS repeats")
check(FI.fragment_id(" [CX3]=O ") == FI.fragment_id("[CX3]=O"), "same pattern -> same ID (whitespace ignored)")
check(FI.fragment_id("[CX3]=O") != FI.fragment_id("[CX3]=S"), "different pattern -> different ID")

# 2 ---------------------------------------------------------------------------------------------------------------
by_id = {str(r["Fragment ID"]): Chem.MolFromSmarts(str(r["SMARTS pattern"])) for r in rows}
by_name = {str(r["Fragment name"]): str(r["Fragment ID"]) for r in rows}
acid, ester = by_name["carboxylic acid"], by_name["ester"]
t = next(x for x in json.load(open(os.path.join(FI.ROOT, "tools", "reaction_tests.json"), encoding="utf-8"))["reactions"] if x["id"] == "R001")
print("      R001:", ".".join(t["substrates"]), ">>", ".".join(t["products"]))
check(any(Chem.MolFromSmiles(s).HasSubstructMatch(by_id[acid]) for s in t["substrates"]), f"{acid} (carboxylic acid) found in a substrate")
check(any(Chem.MolFromSmiles(s).HasSubstructMatch(by_id[ester]) for s in t["products"]), f"{ester} (ester) found in the product")
rules = open(os.path.join(FI.ROOT, "data", "reaction_rules.js"), encoding="utf-8").read()
check(acid in rules and ester in rules, "the generated rules JS refers to those IDs")

# 3 ---------------------------------------------------------------------------------------------------------------
tmp = tempfile.mkdtemp(); real = FI.ROOT
try:
    for f in FI.TEXT_FILES + FI.BOOKS:
        os.makedirs(os.path.dirname(os.path.join(tmp, f)), exist_ok=True); shutil.copy(os.path.join(real, f), os.path.join(tmp, f))
    shutil.copytree(os.path.join(real, "data", "learned", "rx"), os.path.join(tmp, "data", "learned", "rx"))
    os.makedirs(os.path.join(tmp, "tools", "data"), exist_ok=True)
    FI.ROOT = tmp
    # edit the ketone SMARTS inside the copy of the workbook (inline-string cell)
    import zipfile
    old_sm = str(next(r for r in rows if str(r["Fragment ID"]) == by_name["ketone"])["SMARTS pattern"])
    new_sm = old_sm.replace("[CX3](=O)", "[CX3;!$(C(=O)[OX2])](=O)")
    p = os.path.join(tmp, "molecule_data.xlsx"); zin = zipfile.ZipFile(p); out = p + ".new"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zo:
        for it in zin.infolist():
            d = zin.read(it.filename)
            if it.filename.startswith("xl/worksheets/"):
                esc = lambda s: s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
                d = d.decode().replace(">" + esc(old_sm) + "<", ">" + esc(new_sm) + "<").encode()
            zo.writestr(it, d)
    zin.close(); os.replace(out, p)
    ketone_old = by_name["ketone"]; ketone_new = FI.fragment_id(new_sm)
    stale = FI.stale()
    check(stale == {ketone_old: ketone_new}, f"editing one SMARTS makes exactly that ID stale ({ketone_old} -> {ketone_new})")
    FI.apply_map(stale)
    check(not FI.stale(), "after the fix every ID matches its SMARTS again")
    check(ketone_old not in open(os.path.join(tmp, "data", "reaction_rules.js"), encoding="utf-8").read()
          and ketone_new in open(os.path.join(tmp, "data", "reaction_rules.js"), encoding="utf-8").read(), "rules JS: old ID gone, new ID present")
    left = new = 0
    for f in glob.glob(os.path.join(tmp, "data", "learned", "rx", "*.js")):
        rws = json.loads(gzip.decompress(base64.b64decode(re.search(r'"gz":"([^"]+)"', open(f, encoding="utf-8").read()).group(1))))
        for r in rws:
            left += ketone_old in r[3] + " " + r[4]; new += ketone_new in r[3] + " " + r[4]
    check(left == 0 and new > 0, f"learned reactions: old ID gone, new ID in {new} rows")
finally:
    FI.ROOT = real; shutil.rmtree(tmp, ignore_errors=True)

print("\n" + ("ALL TESTS PASSED" if not fails else f"{len(fails)} FAILED"))
sys.exit(1 if fails else 0)
