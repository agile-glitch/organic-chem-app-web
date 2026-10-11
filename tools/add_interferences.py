"""Write the sheet "Interferences" of reaction_rules.xlsx.

    python tools/add_interferences.py --int-model v2 --evidence tools/interference_evidence.json
    python tools/rules_to_js.py                     # then regenerate data/reaction_rules.js, which the page reads

One row per interference, from tools/interference_defs.py (the patterns and advisories written down from textbooks and the literature), with what the
recorded reactions say about each, measured by tools/int_train.py (site patterns: how often the recorded winner had the pattern when the outcomes
differed, the learned weight, the weight per reaction family) and tools/int_evidence.py (advisories: how often the group is consumed when that reagent is used,
and in how many reactions the advisory would be shown).
The pattern, reagent side and group columns are what the program tests; edit them here and run tools/rules_to_js.py. The learned columns are rewritten by this
script; the sheet is the knowledge, js/interference.js reads it, and the weights (data/learned/i/weights.js) are keyed by the Key column.
IDs are kept in tools/data/interference_ids.json so a row keeps its ID when the list changes.
"""
import argparse, json, os, sys, tempfile, time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import interference_defs as D

SHEET = 'Interferences'
HEAD = ['ID', 'Key', 'Name', 'Group', 'Kind', 'Applies to (the atom that reacts, or: reagent + group)',
        'Pattern (SMARTS; first atom = the atom that reacts; [900*] = an atom the reaction changes, [!900*] = one it does not)',
        'Prefilter (SMARTS the product must contain somewhere; it only makes the test faster)',
        'Reagent side (any alternative; + = all of them; frag:F… = a fragment found in a reagent part, smiles:… = a part as typed, role:… = a Reagents-sheet role, solvent:… = a Solvents-sheet type, elem:Pd = some part has that element; not:… = the reagent must NOT have it; empty = any reagent)',
        'Groups (fragment IDs; any one present)', 'Groups are looked for in (substrate, reagent, both)', 'Shown for (bystander = a group the product did not consume; target; both)',
        'Distinct groups needed', 'Needs (solvent, conditions, size, reagent)', 'Used by the page (yes/no)', 'Effect kind', 'Where / when', 'Effect on the reaction', 'Why', 'Hammett σ',
        'Learned weight (log-odds of this site being chosen; + = favoured; for a ring substituent the pattern plus ρ × σ)', 'Learned by reaction family (only where the data say more than the global weight)',
        'Contests where this pattern told the outcomes apart', 'Patents those contests come from', 'Recorded winner had it', 'A random pick would have it',
        'Recorded reactions with the reagent and a group', 'Advisory would be shown in this share of all recorded reactions', 'Group consumed with this reagent', 'Group consumed with other reagents',
        'Evidence note', 'Source']
WIDTH = [8, 30, 44, 26, 18, 18, 52, 30, 44, 34, 12, 14, 9, 12, 9, 12, 50, 70, 50, 9, 14, 40, 14, 12, 12, 12, 14, 14, 14, 14, 50, 44]
PCT = {'Recorded winner had it', 'A random pick would have it', 'Advisory would be shown in this share of all recorded reactions', 'Group consumed with this reagent', 'Group consumed with other reagents'}


def ids_for(keys):
    path = os.path.join(HERE, 'data', 'interference_ids.json')
    ids = json.load(open(path, encoding='utf-8')) if os.path.exists(path) else {}
    n = max([int(v[1:]) for v in ids.values()] + [0])
    for k in keys:
        if k not in ids:
            n += 1
            ids[k] = f'I{n:03d}'
    json.dump(ids, open(path, 'w', encoding='utf-8'), indent=0, sort_keys=True)
    return ids


def frag_names():
    import openpyxl
    wb = openpyxl.load_workbook(os.path.join(ROOT, 'molecule_data.xlsx'), read_only=True)
    out = {}
    it = wb['Fragments'].iter_rows(values_only=True)
    next(it)
    for r in it:
        if r and r[0]: out[str(r[0]).strip()] = str(r[1])
    return out


def reagent_text(alts):
    return ' | '.join(' + '.join(t) for t in alts)


def figure_texts():
    """the measured figures the two 'data support' advisories quote, from tools/eval_top1.json"""
    p = os.path.join(HERE, 'eval_top1.json')
    if not os.path.exists(p): return '', ''
    j = json.load(open(p, encoding='utf-8'))
    pc = lambda x: f'{100 * x:.0f} %' if x is not None else '—'
    s = (j.get('by_size') or {}).get('8-11') or {}
    small = (j.get('by_size') or {}).get('< 8') or {}
    size = (f'Measured on {j.get("n", 0):,} held-out reactions: {small.get("n", 0)} had fewer than 8 heavy atoms (too few to measure anything), and starting materials of 8–11 heavy atoms got a right first answer '
            f'{pc(s.get("top1_exact"))} of the time against a stated chance of {pc(s.get("stated"))}.') if s else ''
    r = (j.get('by_reagent') or {})
    plain = next((v for k, v in r.items() if k.startswith('a reagent without')), None)
    reag = (f'Measured: reactions whose reagent has no reactive role in the Reagents sheet got a right first answer {pc(plain.get("top1_exact"))} of the time overall; the page shows an answer for them '
            f'only at a stated chance of 85 %, and the answers it shows were right {pc(plain.get("shown_right"))} of the time.') if plain else ''
    return size, reag


def build_rows(intm, ev_path):
    names = frag_names()
    ev_adv = json.load(open(ev_path, encoding='utf-8')) if ev_path and os.path.exists(ev_path) else {}
    size_text, reag_text = figure_texts()
    keys = [s['key'] for s in D.SITE] + [a['key'] for a in D.ADVISORY]
    ids = ids_for(keys)
    rows = []
    fams = intm.fams if intm else []
    evid = intm.m.get('evidence', {}) if intm else {}
    kidx = intm.kidx if intm else {}
    for s in D.SITE:
        key = s['key']
        g = n_all = pat = winner = rand = None
        fam_text = ''
        if intm is not None and key in kidx:
            i = kidx[key]
            g = round(intm.total(key), 2)                                  # the pattern's weight, plus rho x sigma for a ring substituent
            e = evid.get((i, 'all'))
            if e: n_all, pat, winner, rand = e['n'], e.get('patents'), e['win'], e['expected']
            parts = []
            for fam in fams:
                ef = evid.get((i, fam))
                if not ef or ef['n'] < 25: continue
                eff = intm.total(key, fam)
                if abs(eff - intm.total(key)) >= 0.25 or abs(eff) >= 0.5:
                    parts.append((abs(eff), f'{fam.replace("_", " ")} {eff:+.2f} (n={ef["n"]}; winner had it {100 * ef["win"]:.0f} % vs {100 * ef["expected"]:.0f} % at random)'))
            fam_text = ' · '.join(t for _, t in sorted(parts, reverse=True)[:6])
        note = ('learned from the training patents (tools/int_train.py)' if n_all and n_all >= 30 else
                'too few contests in the recorded reactions to learn a weight (the weight is shrunk to about zero): the textbook effect is advisory' if intm is not None else '')
        rows.append([ids[key], key, s['name'], s['group'], s['kind'], s['anchor'], s['smarts'], s.get('prefilter') or '', '', '', '', '', None, '', 'yes', '', s['relation'], s['effect'], s['why'], s['sigma'],
                     g, fam_text, n_all, pat, winner, rand, None, None, None, None, note, s['source']])
    for a in D.ADVISORY:
        key = a['key']
        e = ev_adv.get(key, {})
        n = shown = c_with = c_other = None
        note = ''
        if e.get('measured'):
            n, c_with, c_other, shown = e['reactions_with_reagent_and_a_group'], e['consumed_with'], e['consumed_other'], e.get('would_be_shown')
            note = f'{e["pairs_with"]:,} (reaction, group) pairs with this reagent; alternatives counted: {"; ".join(e["reagent_alternatives_used"])}'
            if a['effect_kind'] in ('poisons', 'slows', 'epimerises', 'redirects', 'unreliable', 'cleaved', 'spared', 'competes', 'quenches'):
                note += '. A group being consumed measures this effect only indirectly: the figures say how often the group reacts at all with this reagent'
        elif e:
            note = e.get('why', '')
        elif not a['reagent'] and not a['groups']:
            note = 'not a reagent-and-group pair: a property of the data or of the conditions the records lack'
        effect = a['effect'].replace('{SIZE_FIGURE}', size_text).replace('{REAGENT_FIGURE}', reag_text).strip()
        rows.append([ids[key], key, a['name'], 'reagent + group', a['kind'], 'reagent + group', '', '', reagent_text(a['reagent']), ' | '.join(f'{g} {names.get(g, "")}'.strip() for g in a['groups']),
                     a['group_where'], a['show'], a['min_distinct'], a['needs'], 'yes' if a.get('page', True) and a['needs'] not in ('size', 'reagent', 'conditions') else 'no', a['effect_kind'],
                     '', effect, a['why'], None, None, '', None, None, None, None, n, shown, c_with, c_other, note, a['source']])
    return rows


def write_sheet(rows, path):
    from openpyxl import load_workbook
    from openpyxl.styles import Alignment, Font, PatternFill
    from openpyxl.utils import get_column_letter
    wb = load_workbook(path)
    if SHEET in wb.sheetnames: del wb[SHEET]
    ws = wb.create_sheet(SHEET, index=wb.sheetnames.index('Acids_Bases') + 1)
    ws.append(HEAD)
    for r in rows: ws.append(r)
    head_font = Font(name='Arial', size=10, bold=True, color='FFFFFFFF')
    fill = PatternFill('solid', fgColor='FF1F4E78')
    for c in ws[1]:
        c.font = head_font; c.fill = fill; c.alignment = Alignment(wrap_text=True, horizontal='center', vertical='center')
    body = Font(name='Arial', size=10)
    for row in ws.iter_rows(min_row=2):
        for c in row:
            c.font = body; c.alignment = Alignment(wrap_text=True, vertical='top')
    for j, w in enumerate(WIDTH, 1): ws.column_dimensions[get_column_letter(j)].width = w
    ws.row_dimensions[1].height = 92
    ws.freeze_panes = 'C2'
    ws.auto_filter.ref = f'A1:{get_column_letter(len(HEAD))}{len(rows) + 1}'
    for j, h in enumerate(HEAD, 1):
        letter = get_column_letter(j)
        if h.startswith('Learned weight'):
            for c in ws[letter][1:]:
                if isinstance(c.value, float): c.number_format = '+0.00;-0.00;0.00'
        if h in PCT:
            for c in ws[letter][1:]:
                if isinstance(c.value, float): c.number_format = '0.0%'
    # the explanation in "How to use"
    hw = wb['How to use']
    text = ('Interferences: things that make a template go faster, slower, elsewhere or not at all, to be cross-referenced with the templates when the program decides WHERE a reaction happens (regioselectivity) and WHICH group reacts '
            '(chemoselectivity). A SITE row has a Pattern: a SMARTS whose first atom is the atom the template changes (an ortho nitro group on the ring carbon that reacts, a branch on the carbon that is attacked). '
            'The page runs a learned template with its changed atoms marked (isotope 900 on the product side) and tests the pattern at the marked atoms; [!900*] in a pattern is an atom the reaction leaves alone. '
            'The patterns are tested on the PRODUCT, so they describe what the reaction leaves in place around the atoms it changed. '
            'An ADVISORY row has a Reagent side (a fragment found in a reagent part, a reagent as typed in the Reagents sheet, a role, or a solvent type; "+" = all of them, "|" = any alternative) and Groups (fragment IDs): '
            'when the reagent is used and one of the groups is present in the starting materials but is not what the product consumed, the page says what the reagent also does to it (a Grignard is quenched by an O–H, '
            'LiAlH4 also reduces an ester, Pd/H2 also cleaves a benzyl ether). The learned columns are measured on the recorded reactions (tools/int_train.py, tools/int_evidence.py) and rewritten by tools/add_interferences.py; '
            'the effect, why and source columns are written from textbooks. The data are positive-only (a record shows what happened, never what failed), so a site weight is learned only where a template could have reacted at '
            'several sites and the recorded product shows which one won. To add an interference: add a row (new Key), run tools/rules_to_js.py, and re-learn the weights with tools/int_train.py if it is a site pattern.')
    row = next((r for r in range(1, hw.max_row + 1) if str(hw.cell(row=r, column=1).value or '').startswith('Interferences:')), None)
    if row is None:
        hw.append([text]); row = hw.max_row
    else:
        hw.cell(row=row, column=1).value = text
    c = hw.cell(row=row, column=1)
    c.font = Font(name='Arial', size=10); c.alignment = Alignment(wrap_text=True, vertical='top')
    wb.save(path)
    return len(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--int-model', default='v4'); ap.add_argument('--evidence', default=os.path.join(HERE, 'interference_evidence.json'))
    ap.add_argument('--book', default=os.path.join(ROOT, 'reaction_rules.xlsx'))
    ap.add_argument('--no-model', action='store_true', help='write the sheet without learned columns (before the weights exist)')
    a = ap.parse_args()
    intm = None
    if not a.no_model:
        import int_model as IM
        intm = IM.IntModel(a.int_model)
    rows = build_rows(intm, a.evidence)
    n = write_sheet(rows, a.book)
    print(f'wrote sheet {SHEET} with {n} rows to {a.book} ({time.strftime("%H:%M")})')


if __name__ == '__main__':
    main()
