"""Build Reaction_prediction_and_functional_groups.docx (run tools/fg_census.py first)."""
import json, os, tempfile
import numpy as np
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.section import WD_ORIENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from openpyxl import load_workbook

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TMP = tempfile.gettempdir()
C = json.load(open(os.path.join(TMP, 'fg_census.json')))
J = json.load(open(os.path.join(REPO, 'tools', 'eval_top1.json'), encoding='utf-8'))
I = J['interference']
z = np.load(os.path.join(TMP, 'organic_chem_relfragcounts_new2.npz'), allow_pickle=True)
S = z['S']; fr = [str(x) for x in z['fragments']]
names = {r[0]: r[1] for r in load_workbook(os.path.join(REPO, 'molecule_data.xlsx'), read_only=True)['Fragments'].iter_rows(min_row=2, values_only=True) if r[0]}
mean_all = S.mean(0); share = (S >= 1).mean(0); mean_pres = np.where(share > 0, mean_all / np.maximum(share, 1e-12), 0)
N = S.shape[0]
pct = lambda x, d=1: f'{100 * x:.{d}f} %'

doc = Document()
st = doc.styles['Normal']; st.font.name = 'Calibri'; st.font.size = Pt(10.5)
for s in doc.sections: s.left_margin = s.right_margin = Cm(2); s.top_margin = s.bottom_margin = Cm(2)


def shade(cell, color):
    tcPr = cell._tc.get_or_add_tcPr(); sh = OxmlElement('w:shd'); sh.set(qn('w:val'), 'clear'); sh.set(qn('w:color'), 'auto'); sh.set(qn('w:fill'), color); tcPr.append(sh)


def table(rows, header, widths=None, size=9, align_right_from=1):
    t = doc.add_table(rows=1, cols=len(header)); t.style = 'Table Grid'; t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, h in enumerate(header):
        c = t.rows[0].cells[i]; c.text = ''; p = c.paragraphs[0]; r = p.add_run(h); r.bold = True; r.font.size = Pt(size); r.font.color.rgb = RGBColor(255, 255, 255); shade(c, '1F4E78')
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ''; p = cells[i].paragraphs[0]; r = p.add_run(str(v)); r.font.size = Pt(size)
            if i >= align_right_from: p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    if widths:
        for row in t.rows:
            for i, w in enumerate(widths): row.cells[i].width = Cm(w)
    # repeat header row
    trPr = t.rows[0]._tr.get_or_add_trPr(); h = OxmlElement('w:tblHeader'); h.set(qn('w:val'), 'true'); trPr.append(h)
    return t


def bullets(items):
    for x in items:
        p = doc.add_paragraph(style='List Bullet')
        if isinstance(x, tuple): r = p.add_run(x[0]); r.bold = True; p.add_run(x[1])
        else: p.add_run(x)


def numbered(items):
    for x in items:
        p = doc.add_paragraph(style='List Number')
        r = p.add_run(x[0]); r.bold = True; p.add_run(x[1])


doc.add_heading('How the app predicts reactions, and what is in the fragment, interference and template lists', 0)
doc.add_paragraph('Organic chemistry app, Reactions page. All accuracy figures are measured on patents the model never learned from. '
                  'Counts of functional groups are counts of definitions (fragments, interference rows, templates), explained in section 4.')

doc.add_heading('1. How a reaction is predicted', 1)
doc.add_paragraph('The page tries four tiers in order and uses the first one that applies.')
numbered([
    ('Hand-written workbook rules. ', 'The 102 textbook reactions in reaction_rules.xlsx (NaBH4 ketone reduction, esterification and so on). If the reagent and conditions match one, it answers and shows a mechanism where one exists.'),
    ('Exact look-up in the patent data. ', 'If the exact molecules and reagent appear in a recorded patent reaction, the recorded product is returned.'),
    ('The relevance model (the main learned predictor). ', 'A small neural network trained on about 982,000 recorded US patent reactions guesses which of about 39,000 common patent templates the reaction uses, from the molecules, the reagent and the fragments they contain. '
     'The page runs the 100 most probable templates (each also in a looser form), a ranking model scores the products, and a confidence model says how likely the first answer is to be right.'),
    ('Older fallback. ', 'The earlier reagent-keyed templates answer when the relevance model is not installed or fails.'),
])
doc.add_paragraph('A relevance-model answer is shown only if its confidence passes a floor: 50 % with a reagent that has a reactive role, 70 % with no reagent, 85 % with a reagent the sheet does not classify. '
                  'It never predicts when the starting materials have fewer than 8 heavy atoms in all, because the patent data has almost nothing that small. Below the floor the page says "Unknown".')
doc.add_heading('What the ranking model uses', 2)
bullets([
    ('Template probabilities ', 'from the network, summed over the templates that give the same product.'),
    ('Fragments. ', 'Each fragment is a named functional-group pattern (205 in the sheet). They are inputs to the network, they feed a reactivity prior (which groups a given reagent usually consumes), and they supply the groups for the interference advisories.'),
    ('Interferences. ', f'{I["patterns"]} site patterns (electron-donating or withdrawing groups at ortho/meta/para, crowding, ring type, carbon and carbonyl substitution) with weights learned from reactions where a template could have reacted at several sites, plus {I["advisories"]} reagent-and-group advisories (Grignard with an alcohol, LiAlH4 with a nitro group, TFA with a Boc group...). A second table learns how much the other groups present change whether a group reacts.'),
])

doc.add_heading('2. How well it works', 1)
doc.add_paragraph('Measured on 12,000 reactions from patents held out of everything learned (split by patent, never by reaction):')
table([
    ['First answer right (exact, with stereochemistry)', pct(J['top1_exact'])],
    ['First answer right (stereochemistry ignored)', pct(J['top1_flat'])],
    ['With a reagent given / with no reagent', f'{pct(J["with_reagent"]["top1_exact"])} / {pct(J["no_reagent"]["top1_exact"])}'],
    ['Ceiling: recorded product is among the candidates at all', pct(J['generated'])],
    ['Without the interference layer', pct(I['top1']['base']['top1_exact'])],
    ['Gain from the interference layer (paired, 95 % interval)', f'+{100 * I["paired"]["top1_exact"]["diff"]:.2f} points ({100 * I["paired"]["top1_exact"]["ci95"][0]:+.2f} to {100 * I["paired"]["top1_exact"]["ci95"][1]:+.2f})'],
    ['Answers rated 90 % sure or better: share of reactions / right', f'{pct(J["by_confidence"][">= 0.9"]["share_of_all"])} / {pct(J["by_confidence"][">= 0.9"]["right"])}'],
    ['Page shows an answer for / right when shown', f'{pct(J["shown"]["shown_share_of_all"])} / {pct(J["shown"]["shown_right"])}'],
    ['Separate spot check: 1,000 fresh random held-out reactions, exact match', '68.9 % (689 of 1,000); 698 when stereochemistry is ignored'],
    ['Site choice when a template could react at several places', f'{pct(I["site_choice"]["learned"])} right (random pick {pct(I["site_choice"]["random"])}), {I["site_choice"]["n"]:,} cases'],
], ['Measure', 'Value'], widths=[11, 6], size=9.5, align_right_from=1)
doc.add_paragraph()
doc.add_paragraph('In the 1,000-reaction spot check the misses were: 200 where no template produced the recorded product, 74 where the model predicted a different reaction, 22 where it picked the wrong isomer (wrong site), and 6 with no answer. '
                  'Reaching 100 % is not possible on this data: about 19 % of recorded products cannot be produced by any learned template (records that combine two steps, unusual chemistry, data errors).')

doc.add_heading('3. Fragments in the patent reactions', 1)
doc.add_paragraph(f'Counting every fragment in the starting materials of all {N:,} recorded reactions (reagents not counted): on average a reaction has {S.sum(1).mean():.1f} fragment hits spread over {(S >= 1).sum(1).mean():.1f} different fragments. '
                  f'Averaged over all 205 fragments, a fragment occurs {mean_all.mean():.2f} times per reaction, and {mean_pres[share > 0].mean():.2f} times when it is present. The 25 most common:')
order = np.argsort(-mean_all)
table([[fr[i], names.get(fr[i], '')[:60], f'{mean_all[i]:.2f}', pct(share[i]), f'{mean_pres[i]:.2f}'] for i in order[:25]],
      ['ID', 'Fragment', 'Average per reaction', 'Reactions containing it', 'Average when present'], widths=[1.5, 7.5, 2.8, 2.9, 2.6], size=9, align_right_from=2)
doc.add_paragraph('The complete list of all 205 fragments is in the appendix. Fragments that score zero (LiAlH4, DIBAL-H, NaH, Pd, halide, azide and fluoride anions) are reagent-type groups and occur only in reagents.')

# ------------------------------------------------------------------------------ landscape section for the census
sec = doc.add_section(); sec.orientation = WD_ORIENT.LANDSCAPE; sec.page_width, sec.page_height = sec.page_height, sec.page_width
sec.left_margin = sec.right_margin = Cm(1.5)
doc.add_heading('4. How often each functional group appears', 1)
doc.add_heading('Two words used in the table: entries and occurrences', 2)
doc.add_paragraph('The table counts each functional group in two ways. Both are simple counts.')
p = doc.add_paragraph(style='List Bullet'); r = p.add_run('Entry. '); r.bold = True
p.add_run('An entry is one item in a list: one fragment, one interference row, or one template. "Entries" answers the question: in how many items does this group appear at all?')
p = doc.add_paragraph(style='List Bullet'); r = p.add_run('Occurrence. '); r.bold = True
p.add_run('An occurrence is one single appearance of the group. "Occurrences" answers the question: how many times does this group appear in total, counting repeats inside the same item?')
doc.add_paragraph('Example. Suppose three fragments contain an alcohol (an OH on a carbon), and one of those three has two OH groups. Then the alcohol has 3 entries (three fragments contain it) and 4 occurrences '
                  '(one in each of two fragments, and two in the third). Occurrences is always equal to or higher than entries.')
doc.add_heading('Why it helps to read both', 3)
bullets([
    ('Entries show how widely a group is covered. ', 'A group with many entries is described in many places, so the program has many separate chances to recognise it. A group with few entries is only described in a few places, however often it repeats there.'),
    ('Occurrences show how heavily a group is represented. ', 'A big gap between occurrences and entries means the group usually comes in repeats inside one item (for example benzene rings: many rings per template).'),
    ('Together they show where coverage is thin. ', 'A group with many template entries but few fragment or interference entries is common in the patent chemistry but weakly described by our fragments and interferences. That is where a new fragment or interference would add the most.'),
    ('One number alone can mislead. ', 'A single large template, or a single row that lists many substituents, can push occurrences up without making the group any better covered. Entries are not affected by that, so reading the two side by side keeps the picture honest.'),
])
doc.add_paragraph('Both numbers count definitions (fragments, interference rows, templates), not molecules. How often a group occurs in the real patent molecules is the separate "average per reaction" table in section 3.')
nf, ns, na, nt = C['n_fragments'], C['n_site'], C['n_adv'], C['n_templates']
doc.add_paragraph(f'The table gives two numbers for each functional group and list: the number of entries that contain the group, and the total number of times the group occurs (an entry with two alcohols counts once as an entry and twice as occurrences). The lists are: {nf} fragments (the named group patterns), {ns + na} interference rows ({ns} site patterns + {na} advisories) and {nt:,} templates (the reaction templates the network chooses between).')
bullets([
    ('Fragments: ', 'the fragment contains the group in its pattern.'),
    ('Interferences: ', 'a site pattern contains the group in its pattern or in the surroundings it tests; an advisory contains it through the fragments it names as groups. The site-pattern and advisory counts are added.'),
    ('Fragments + interferences: ', f'the two lists added (entries out of {nf + ns + na}).'),
    ('Occurrences: ', 'distinct matches of the group inside the skeleton of each entry, added over the list. For a site pattern every group named in its surroundings is counted; for an advisory, through the fragments it names. For a template, occurrences on the two sides are added in the "either side" column.'),
    ('Templates: ', 'the group is in the reacting part of the template, on the reactant side (what the reaction consumes), the product side (what it makes), or either side.'),
])
rows = sorted(C['rows'], key=lambda r: -r['tpl_any'])
f = lambda v: f'{v:,}'
table([[r['group'], f(r['fragments']), f(r['o_fragments']), f(r['interferences']), f(r['o_interferences']), f(r['both']), f(r['o_both']),
        f(r['tpl_reactant']), f(r['o_tpl_reactant']), f(r['tpl_product']), f(r['o_tpl_product']), f(r['tpl_any']), f(r['o_tpl_any'])] for r in rows],
      ['Functional group', f'Fragments: entries (of {nf})', 'Fragments: occurrences', f'Interferences: entries (of {ns + na})', 'Interferences: occurrences',
       f'Fragments + interferences: entries (of {nf + ns + na})', 'Fragments + interferences: occurrences',
       f'Templates, reactant side: entries (of {nt:,})', 'Templates, reactant side: occurrences', 'Templates, product side: entries', 'Templates, product side: occurrences',
       'Templates, either side: entries', 'Templates, either side: occurrences'],
      widths=[5.0, 1.8, 1.9, 1.9, 1.9, 2.2, 2.2, 2.0, 2.1, 1.9, 2.1, 2.0, 2.1], size=7.5, align_right_from=1)
doc.add_paragraph()
doc.add_heading('Reading the table', 2)
bullets([
    f'Templates are by far the largest list, so their counts are in the thousands (an aromatic carbon is in {rows[0]["tpl_any"]:,} of {nt:,}); fragments and interference rows are hundreds at most.',
    'Groups that are common in templates but rare in the fragments or interferences are where a new fragment or interference would add the most coverage (for example furan / thiophene rings, ureas, thiols, isocyanates, boron and phosphorus groups, which have no or one or two fragments).',
    'Alcohol: ' + f'{[r for r in rows if r["group"].startswith("Alcohol")][0]["fragments"]} fragments ({[r for r in rows if r["group"].startswith("Alcohol")][0]["o_fragments"]} occurrences), {[r for r in rows if r["group"].startswith("Alcohol")][0]["interferences"]} interference rows ({[r for r in rows if r["group"].startswith("Alcohol")][0]["o_interferences"]} occurrences), {[r for r in rows if r["group"].startswith("Alcohol")][0]["tpl_any"]:,} templates ({[r for r in rows if r["group"].startswith("Alcohol")][0]["o_tpl_any"]:,} occurrences). '
    'Epoxide: ' + f'{[r for r in rows if r["group"].startswith("Epoxide")][0]["fragments"]} fragment, {[r for r in rows if r["group"].startswith("Epoxide")][0]["interferences"]} interference rows, {[r for r in rows if r["group"].startswith("Epoxide")][0]["tpl_any"]:,} templates.',
])

doc.add_heading('5. How the counts were made, and their limits', 1)
bullets([
    'Each fragment, site pattern and template is written as a SMARTS pattern, not as a molecule. Each pattern was turned into a concrete skeleton molecule (elements, charges, aromaticity, bond orders and the hydrogen counts the pattern states) and the 60 functional-group definitions were searched in it.',
    'This is a count of definitions, not of molecules. It does not say how often a group occurs in the patent data; for that use the fragment averages in section 3.',
    'A template contains only the atoms it changes and their neighbours, so its groups are the ones that take part or sit next to the reacting centre; a group elsewhere in a molecule is not in the template.',
    'Recursive conditions inside site patterns ("a nitro group para to the reacting carbon") were read as positive groups; conditions that exclude a group ("not an amide") were ignored. The 20 aryl electron-withdrawing / donating rows each list all the groups of their substituent set, so groups such as nitro, nitrile and sulfonyl are counted in them.',
    'The groups overlap by design (a primary amine is also an aniline-type N when on a ring; an ester is also an ether-like C–O). A skeleton built from a pattern that leaves an atom unspecified assumes carbon.',
    'The functional-group list is a chosen list of 60 common groups; groups outside it are not counted.',
])

doc.add_page_break()
doc.add_heading('Appendix: all 205 fragments', 1)
doc.add_paragraph(f'Average count per recorded reaction (starting materials of {N:,} reactions), share of reactions containing the fragment, and average count when present.')
table([[fr[i], names.get(fr[i], '')[:95], f'{mean_all[i]:.3f}', pct(share[i], 2), f'{mean_pres[i]:.2f}'] for i in range(len(fr))],
      ['ID', 'Fragment', 'Average per reaction', 'Reactions containing it', 'Average when present'], widths=[1.8, 15.5, 3.2, 3.6, 3.2], size=8, align_right_from=2)

out = os.path.join(REPO, 'Reaction_prediction_and_functional_groups.docx')
doc.save(out)
print('saved', out)
