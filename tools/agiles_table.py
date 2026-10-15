"""The AGILES spreadsheet: every compound of the app's library with its AGILES code and identifiers, from
tools/agiles_table.js (the Molecule tab's own values), each code checked before it is written:
  - the code unpacks (raw DEFLATE, base64url) to an AGILES record whose items are the row's identifiers;
  - RDKit Python reads the row's SMILES (as given) and RDKit SMILES as the same molecule, whose InChIKey is the row's;
  - the drawing represents itself: the code's MOL block alone (atoms, bonds, wedges, 2D coordinates, nothing else) is
    read by RDKit Python as a molecule with the same full InChIKey, so its stereo (R/S from the wedges, E/Z from the
    double bonds' drawn geometry) is that molecule's and no other's.

  python tools/agiles_table.py [tools/data/agiles_table.json]   -> AGILES_codes.xlsx, AGILES_codes.csv
"""
import base64, csv, json, os, sys, zlib
from rdkit import Chem, RDLogger
from rdkit.Chem import inchi as rdinchi

RDLogger.DisableLog('rdApp.*')
HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
COLS = [('Name', 'name'), ('AGILES code', 'code'), ('IUPAC name', 'iupac'), ('IUPAC name checked', 'iupac_checked'),
        ('SMILES (as given)', 'smiles'), ('SMILES (RDKit canonical)', 'smiles_rdkit'), ('InChI', 'inchi'), ('InChIKey', 'inchikey'), ('CAS number', 'cas'),
        ('Formula (with charge)', 'formula'), ('Isotopic formula', 'isotopic_formula'), ('Molecular weight', 'mw'), ('Monoisotopic mass', 'monoisotopic')]


NOTE = ('These are EXAMPLE molecules: the {n:,} compounds of the app\'s library. AGILES1 codes are not limited to them. '
        'Every molecule the Molecule tab can show (typed, named, opened from a file, or drawn in the Sketcher) has an AGILES1 '
        'code: show it and choose Export → AGILES code (copy).')


def unpack(code):
    b = code.split(':', 1)[1]
    return zlib.decompress(base64.urlsafe_b64decode(b + '=' * (-len(b) % 4)), -15).decode('utf-8')


def items(sdf):
    out, cur = {}, None
    for line in sdf.split('M  END', 1)[1].split('\n'):
        if line.startswith('>'):
            cur = line[line.index('<') + 1:line.index('>', line.index('<'))]; out[cur] = []
        elif cur is not None and line.strip() and not line.startswith('$$$$'): out[cur].append(line)
        else: cur = None
    return {k: ' '.join(v) for k, v in out.items()}


def macrocycle_only(md):
    """the drawing's InChIKey when the E/Z of double bonds in rings of 8 or more atoms is not read from the drawing:
    InChI reads one for such a ring bond (a porphyrin's), where the ring allows only the geometry drawn, while the same
    molecule from SMILES has none"""
    m = Chem.Mol(md)
    ri = m.GetRingInfo()
    for b in m.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE and any(len(r) >= 8 and b.GetIdx() in r for r in ri.BondRings()):
            b.SetStereo(Chem.BondStereo.STEREOANY)
    m.RemoveAllConformers()
    return rdinchi.MolToInchiKey(m)


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, 'data', 'agiles_table.json')
    rows = json.load(open(src))['rows']
    good, bad, macro = [], [], []
    for r in rows:
        if not r.get('ok'):
            bad.append((r['name'], r.get('error') or 'not shown')); continue
        ids = r['ids']
        row = {'name': r['name'], 'code': r['code'], **{k: ids.get(k) or '' for _, k in COLS if k not in ('name', 'code')}}
        why = None
        try:
            sdf = unpack(r['code'])
            it = items(sdf)
            for _, k in COLS[2:]:
                if (it.get('AGILES.' + k) or '') != row[k]: why = f'the code\'s {k} is not the row\'s'; break
        except Exception as e:
            why = f'the code does not unpack ({e})'
        if not why:
            m1, m2 = Chem.MolFromSmiles(row['smiles']), Chem.MolFromSmiles(row['smiles_rdkit'])
            if m1 is None or m2 is None: why = 'RDKit Python cannot read a SMILES'
            elif Chem.MolToSmiles(m1) != Chem.MolToSmiles(m2): why = 'the two SMILES are not the same molecule'
            elif row['inchikey'] and rdinchi.MolToInchiKey(m2) != row['inchikey']: why = 'the InChIKey is not the SMILES\' molecule'
            else:
                md = Chem.MolFromMolBlock(sdf.split('M  END')[0] + 'M  END\n')
                if md is None: why = 'RDKit Python cannot read the drawing (MOL block)'
                elif row['inchikey'] and rdinchi.MolToInchiKey(md) != row['inchikey']:
                    if macrocycle_only(md) == row['inchikey']: macro.append(row['name'])
                    else: why = f'the drawing reads as another molecule or stereoisomer ({rdinchi.MolToInchiKey(md)})'
        (bad.append((r['name'], why)) if why else good.append(row))
    good.sort(key=lambda x: x['name'].lower())
    xlsx, csvp = os.path.join(APP, 'AGILES_codes.xlsx'), os.path.join(APP, 'AGILES_codes.csv')
    with open(csvp, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f); w.writerow([c for c, _ in COLS])
        for g in good: w.writerow([g[k] for _, k in COLS])
    global NOTE
    NOTE = NOTE.format(n=len(good))
    from openpyxl import Workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    wb = Workbook(); ws = wb.active; ws.title = 'AGILES codes'
    ws.append([NOTE]); ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(COLS))
    ws['A1'].font = Font(bold=True, color='7A5B00'); ws['A1'].fill = PatternFill('solid', fgColor='FFF7E0')
    ws['A1'].alignment = Alignment(wrap_text=True, vertical='center'); ws.row_dimensions[1].height = 48
    ws.append([c for c, _ in COLS])
    for c in ws[2]: c.font = Font(bold=True, color='FFFFFF'); c.fill = PatternFill('solid', fgColor='2F5D50')
    for g in good: ws.append([g[k] for _, k in COLS])
    for col, wdt in zip('ABCDEFGHIJKLM', (28, 60, 40, 10, 40, 40, 50, 30, 14, 16, 14, 14, 16)): ws.column_dimensions[col].width = wdt
    ws.freeze_panes = 'C3'; ws.auto_filter.ref = f'A2:{chr(64 + len(COLS))}{ws.max_row}'
    for row in ws.iter_rows(min_row=3, min_col=2, max_col=2):
        for c in row: c.alignment = Alignment(wrap_text=False)
    about = wb.create_sheet('About')
    for line in [
        NOTE,
        'AGILES codes for the compounds of the app\'s library, made by the Molecule tab itself (tools/agiles_table.js, tools/agiles_table.py).',
        'Paste an AGILES code (it starts with AGILES1:) into the Molecule tab\'s box: the molecule comes back drawn as it was, with every label and these identifiers.',
        'The codes hold no bond angles (no 3D model was built for them); the tab measures the angles when a code is shown.',
        'IUPAC name: the app\'s namer; "no" in the checked column means its name could not be read back to confirm it (shown "unchecked" in the app). A blank cell: no name.',
        'SMILES (as given): the library\'s own SMILES. SMILES (RDKit canonical): how the tab reads the molecule and numbers its atoms in the drawing.',
        'CAS number: from the app\'s library; a blank cell: the library has none for this compound.',
        f'Every code was unpacked and checked: its identifiers are the row\'s, both SMILES are one molecule whose InChIKey (RDKit) is the row\'s, and the drawing alone (the MOL block: atoms, bonds, wedges, coordinates) reads back as that same molecule and stereoisomer. {len(good):,} rows; {len(bad)} compounds left out (no code or a check failed).',
        f'{len(macro)} of them have a double bond in a ring of 8 or more atoms (a porphyrin\'s): InChI reads an E/Z for it from the drawing that the ring allows only one way, so their drawing\'s InChIKey has that extra layer; without it, it is the row\'s. ' + ', '.join(macro[:20]),
    ]: about.append([line])
    about.column_dimensions['A'].width = 160
    wb.save(xlsx)
    print(f'{len(good):,} rows -> {os.path.relpath(xlsx, APP)}, {os.path.relpath(csvp, APP)}; {len(bad)} left out')
    for n, w in bad[:20]: print('  left out:', n, '-', w)


if __name__ == '__main__':
    main()
