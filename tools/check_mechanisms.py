"""Replay check for the curved arrows in reaction_rules.xlsx (sheet Mechanism).

For every step, the arrows are APPLIED to 'Structure the arrows act on' (atom-mapped SMILES) with electron
bookkeeping, and every species the arrows touched must come out as one of the species in
'Intermediate / product after the step'. That proves the arrows really turn the start of the step into its result.

Arrow notation (same as the sheet):  a>b  a>b-c  a-b>c  a-b>c-d   (~ in front = one electron, fishhook)
  a>b       lone pair on a forms a bond a–b                      a: +1 charge, b: -1
  a>a-b     the same, written as a new bond                       (the lone pair must be on a or b)
  a-b>a     bond a–b breaks, both electrons go to a               a: -1, b: +1
  a-b>b-c   the bond's electrons shift to a new bond b–c          a (left behind): +1, c (gained): -1
  ~a-b>a    one electron of a–b goes to a (homolysis half)        bond -1/2, a gains a radical
  ~a>a-b    a's unpaired electron goes into a new bond a–b        bond +1/2, a loses a radical
  ~a-b>b-c  one electron shifts from a–b to b–c                   bond orders ±1/2
  ~a>b      one electron moves from a to b (electron transfer)    a: +1, b: -1, radical moves
Hydrogens move only as numbered [H:n] atoms; every other H count stays as written.

Usage:
    python tools/check_mechanisms.py            report
    python tools/check_mechanisms.py --write    report, and write the 'Mechanism status' column into reaction_rules.xlsx
"""
import os, re, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOK = re.compile(r"^(~)?(\d+)(?:-(\d+))?>(\d+)(?:-(\d+))?$")

def _rdkit():
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        return Chem
    except ImportError:
        return None

def canon_set(Chem, smiles):
    """canonical SMILES of each species, explicit H removed, stereo ignored"""
    out = set()
    for part in str(smiles or "").split("."):
        if not part.strip():
            continue
        m = Chem.MolFromSmiles(part)
        if m is None:
            out.add("UNPARSEABLE:" + part)
            continue
        Chem.RemoveStereochemistry(m)
        out.add(Chem.MolToSmiles(Chem.RemoveHs(m)))
    return out

def replay(Chem, before, arrows, after):
    """apply the arrows; return (ok, message)"""
    ps = Chem.SmilesParserParams(); ps.removeHs = False
    raw = Chem.SmilesParserParams(); raw.removeHs = False; raw.sanitize = False
    mol = Chem.MolFromSmiles(before, raw)
    if mol is None:
        return False, "the 'structure the arrows act on' does not parse"
    try:   # keep the single/double bonds exactly as written (no re-assigned aromatic bonds)
        Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    except Exception as e:
        return False, f"the 'structure the arrows act on' is not a valid structure ({str(e).splitlines()[0][:80]})"
    Chem.Kekulize(mol, clearAromaticFlags=True)            # arrows act on the drawn single/double bonds of a ring
    hcount = {a.GetIdx(): a.GetTotalNumHs() for a in mol.GetAtoms()}   # H counts stay as written
    idx = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    order = defaultdict(float)                     # bond orders (by atom-index pair), as numbers
    for b in mol.GetBonds():
        i, j = sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))
        order[(i, j)] = b.GetBondTypeAsDouble()
    charge = {a.GetIdx(): a.GetFormalCharge() for a in mol.GetAtoms()}
    rad = {a.GetIdx(): a.GetNumRadicalElectrons() for a in mol.GetAtoms()}
    touched = set()
    def bond(i, j, d):
        k = tuple(sorted((i, j))); order[k] += d
    # one electron arriving pairs with an unpaired one already there; one leaving a pair leaves the other unpaired
    def add_e(i): rad[i] = rad[i] - 1 if rad[i] > 0 else rad[i] + 1
    def take_e(i): rad[i] = rad[i] - 1 if rad[i] > 0 else rad[i] + 1
    for t in [x.strip() for x in str(arrows or "").split(";") if x.strip()]:
        m = TOK.match(t)
        if not m:
            return False, f"arrow '{t}' is not in the notation"
        fish = bool(m.group(1))
        try:
            a = idx[int(m.group(2))]; b = idx[int(m.group(3))] if m.group(3) else None
            c = idx[int(m.group(4))]; d = idx[int(m.group(5))] if m.group(5) else None
        except KeyError as e:
            return False, f"arrow '{t}' uses atom {e.args[0]}, which is not numbered in the structure"
        touched.update(x for x in (a, b, c, d) if x is not None)
        h = 0.5 if fish else 1.0
        if b is None:                               # from a lone pair (or a single electron) on a
            if d is None:                           # a>c
                if fish and c != a:                 # ~a>c: electron transfer
                    charge[a] += 1; charge[c] -= 1; take_e(a); add_e(c)
                else:
                    bond(a, c, h); charge[a] += (0 if fish else 1); charge[c] -= (0 if fish else 1)
                    if fish: take_e(a)
            else:                                   # a>c-d: the lone pair must be on c or d
                if a not in (c, d):
                    return False, f"arrow '{t}': a lone pair can only form a bond to its own atom (a>a-b)"
                o = d if a == c else c
                bond(a, o, h)
                if fish: take_e(a)
                else: charge[a] += 1; charge[o] -= 1
        else:                                       # from the a–b bond
            bond(a, b, -h)
            if d is None:                           # to an atom
                if c not in (a, b):
                    return False, f"arrow '{t}' is ambiguous: bond electrons going to a third atom; write which new bond forms (a-b>b-c)"
                o = b if c == a else a
                if fish: add_e(c)
                else: charge[c] -= 1; charge[o] += 1
            else:                                   # to a new bond c–d
                shared = {a, b} & {c, d}
                if len(shared) != 1:
                    return False, f"arrow '{t}': the new bond must share one atom with the old one"
                s = shared.pop(); left = ({a, b} - {s}).pop(); gained = ({c, d} - {s}).pop()
                bond(s, gained, h)
                if not fish: charge[left] += 1; charge[gained] -= 1
    # rebuild the molecule
    rw = Chem.RWMol(mol)
    for bd in list(rw.GetBonds()):
        rw.RemoveBond(bd.GetBeginAtomIdx(), bd.GetEndAtomIdx())
    types = {1.0: Chem.BondType.SINGLE, 2.0: Chem.BondType.DOUBLE, 3.0: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC}
    for (i, j), o in order.items():
        if abs(o) < 1e-9:
            continue
        if o < 0 or (o not in types):
            return False, f"the arrows leave a bond with order {o:g} between atoms {mol.GetAtomWithIdx(i).GetAtomMapNum() or i} and {mol.GetAtomWithIdx(j).GetAtomMapNum() or j}"
        rw.AddBond(i, j, types[o])
    for a in rw.GetAtoms():
        a.SetFormalCharge(int(charge[a.GetIdx()]))
        if rad[a.GetIdx()] < 0:
            return False, f"atom {a.GetAtomMapNum() or a.GetIdx()} gives away an electron it does not have (radical count below 0)"
        a.SetNumRadicalElectrons(int(rad[a.GetIdx()]))
        a.SetNoImplicit(True)                       # H counts stay as written (H moves only as [H:n] atoms)
        a.SetNumExplicitHs(hcount[a.GetIdx()])
        a.SetIsAromatic(False)
    res = rw.GetMol()
    try:
        Chem.SanitizeMol(res)
    except Exception as e:
        return False, f"the result is not a valid structure ({str(e).splitlines()[0][:90]})"
    frags = Chem.GetMolFrags(res, asMols=False)
    got_touched = set()
    for f in frags:
        if not touched & set(f):
            continue
        sub = Chem.RWMol(Chem.PathToSubmol(res, [b.GetIdx() for b in res.GetBonds() if b.GetBeginAtomIdx() in f]) if len(f) > 1 else None) if False else None
        fm = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(res, atomsToUse=list(f), canonical=True), ps)
        if fm is None:
            return False, "a species made by the arrows could not be read back"
        for a in fm.GetAtoms():
            a.SetAtomMapNum(0)
        Chem.RemoveStereochemistry(fm)
        got_touched.add(Chem.MolToSmiles(Chem.RemoveHs(fm)))
    expected = canon_set(Chem, after)
    # species of the start that no arrow touches stay as they are (spectators)
    unchanged = set()
    for f in Chem.GetMolFrags(mol, asMols=False):
        if touched & set(f):
            continue
        fm = Chem.MolFromSmiles(Chem.MolFragmentToSmiles(mol, atomsToUse=list(f), canonical=True), ps)
        if fm is not None:
            for a in fm.GetAtoms(): a.SetAtomMapNum(0)
            Chem.RemoveStereochemistry(fm)
            unchanged.add(Chem.MolToSmiles(Chem.RemoveHs(fm)))
    # the test: everything the step claims as its result must be made by the arrows (or be an untouched spectator)
    missing = sorted(expected - got_touched - unchanged)
    if missing:
        return False, f"the step's result has {'.'.join(missing)}, but the arrows make {'.'.join(sorted(got_touched)) or 'nothing new'}"
    extra = sorted(got_touched - expected)
    return True, ("ok" if not extra else "ok; not listed in the step's result: " + ".".join(extra))

def check_rows(rows, Chem):
    """rows: dicts with id, step, before, arrows, after. Returns {rid: {"status", "problems", "steps"}}"""
    by = defaultdict(list)
    for r in rows:
        by[r["id"]].append(r)
    out = {}
    for rid, steps in by.items():
        probs, notes, noarrow = [], [], 0
        for st in sorted(steps, key=lambda x: x["step"]):
            if not st["before"] or not st["arrows"]:
                noarrow += 1
                continue
            ok, msg = replay(Chem, st["before"], st["arrows"], st["after"])
            if not ok:
                probs.append(f"step {st['step']}: {msg}")
            elif msg != "ok":
                notes.append(f"step {st['step']}: {msg[4:]}")
        if probs:
            status = "drafted (replay failed)"
        elif noarrow == len(steps):
            status = "drafted (no arrows yet)"
        elif noarrow:
            status = f"verified (replay passed; {noarrow} step{'s' if noarrow > 1 else ''} without arrows)"
        else:
            status = "verified (replay passed)"
        out[rid] = {"status": status, "problems": probs, "notes": notes, "steps": len(steps)}
    return out

def read_workbook(path=None):
    from openpyxl import load_workbook
    path = path or os.path.join(ROOT, "reaction_rules.xlsx")
    wb = load_workbook(path, read_only=True)
    ws = wb["Mechanism"]
    it = ws.iter_rows(values_only=True)
    head = [str(h) for h in next(it)]
    col = lambda start: next(i for i, h in enumerate(head) if h.startswith(start))
    ci, cs, ca, cb, cr = col("Reaction ID"), col("Step"), col("Intermediate"), col("Structure the arrows"), col("Arrows (")
    rows = []
    for r in it:
        if not r or not isinstance(r[cs], (int, float)) or not r[ci]:
            continue
        rows.append({"id": str(r[ci]), "step": int(r[cs]), "after": r[ca] or "", "before": r[cb] or "", "arrows": r[cr] or ""})
    rx = wb["Reactions"]
    ids = [str(r[0]) for r in rx.iter_rows(min_row=2, values_only=True) if r and r[0]]
    return rows, ids

def report(results, ids):
    ver = [i for i in ids if results.get(i, {}).get("status", "").startswith("verified")]
    bad = [i for i in ids if results.get(i, {}).get("status", "") == "drafted (replay failed)"]
    noar = [i for i in ids if results.get(i, {}).get("status", "") == "drafted (no arrows yet)"]
    none = [i for i in ids if i not in results]
    lines = [f"mechanisms: {len(ver)} of {len(ids)} verified by replay, {len(bad)} failing, {len(noar)} without arrows, {len(none)} with no mechanism yet"]
    for i in bad:
        lines += [f"  {i}: {p}" for p in results[i]["problems"]]
    if noar:
        lines.append("  no arrows yet: " + ", ".join(noar))
    if none:
        lines.append("  no mechanism yet: " + ", ".join(none))
    return "\n".join(lines)

def statuses(path=None):
    Chem = _rdkit()
    if Chem is None:
        return None, "mechanism replay skipped (RDKit is not installed)"
    rows, ids = read_workbook(path)
    res = check_rows(rows, Chem)
    full = {i: (res[i] if i in res else {"status": "none", "problems": [], "steps": 0}) for i in ids}
    return full, report(res, ids)

def write_status(full, path=None):
    from openpyxl import load_workbook
    path = path or os.path.join(ROOT, "reaction_rules.xlsx")
    wb = load_workbook(path)
    ws = wb["Reactions"]
    head = [c.value for c in ws[1]]
    if "Mechanism status" in head:
        c = head.index("Mechanism status") + 1
    else:
        c = len(head) + 1
        src = ws.cell(1, c - 1)
        cell = ws.cell(1, c, "Mechanism status")
        from copy import copy
        cell.font, cell.fill, cell.border, cell.alignment = copy(src.font), copy(src.fill), copy(src.border), copy(src.alignment)
        ws.column_dimensions[cell.column_letter].width = 30
    for r in range(2, ws.max_row + 1):
        rid = ws.cell(r, 1).value
        if rid in full:
            ws.cell(r, c, full[rid]["status"])
    wb.save(path)

if __name__ == "__main__":
    full, text = statuses()
    print(text)
    if full is not None and "--write" in sys.argv:
        write_status(full)
        print("wrote the 'Mechanism status' column in reaction_rules.xlsx")
