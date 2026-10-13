"""What did chemists use for this kind of reaction? Tally the recorded reagents, solvents and yields of every learned reaction type.

    python tools/cond_build.py            writes data/learned/c/NNN.js  (and prints what it kept)

data/learned/rx/ holds the 1.09 million recorded reactions the Reactions page learned from (US patents 1976-2016 text-mined by D. M. Lowe, CC0, and the
Open Reaction Database, CC-BY-SA 4.0), each as  [id, "starting materials>reagent>product | cond=solvents", template #, consumed, formed, leaving, times recorded,
source, year, yield %]. A template is one KIND of reaction (the changed atoms and their neighbours); the page's relevance model proposes templates for
a set of starting materials. This reads every record once and keeps, per template with at least MIN_N records:

    n          records (each counted as many times as it was recorded)
    yield      the median recorded yield (%), where recorded
    recipes    the reagent lines used, most common first: [reagent SMILES line, count, a patent that used it, its yield]
    solvents   the recorded solvent sets, most common first: [solvents as written in the conditions box, count]
    consumed, formed   the groups the reaction type most often consumes and forms ("F004 + F095", Fragments sheet ids)

so that "NaBH4 and methanol" can be offered for a ketone reduction because that is what 4,000 recorded reductions used, and shown with that number.
Patent records have no temperature or time (the text-mining did not keep them), so none is offered from here. data/learned/c/NNN.js holds templates
NNN*2000 .. NNN*2000+1999 (gzip inside base64, like the other learned parts).
"""
import base64, collections, gzip, json, os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rxn_code  # noqa: E402
import statistics

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
RX = os.path.join(APP, 'data', 'learned', 'rx')
OUT = os.path.join(APP, 'data', 'learned', 'c')
PER_SHARD = 2000
MIN_N = 3
TOP_RECIPES, TOP_SOLVENTS = 6, 5
MAX_YIELDS = 300


def main():
    files = sorted(f for f in os.listdir(RX) if f.endswith('.js'))
    recipes = collections.defaultdict(collections.Counter)
    example = {}                                        # (template, recipe) -> (yield, source) of the best-yielding record
    solv = collections.defaultdict(collections.Counter)
    yields = collections.defaultdict(list)
    cons, form = collections.defaultdict(collections.Counter), collections.defaultdict(collections.Counter)
    count = collections.Counter()
    total = 0
    for fn in files:
        s = open(os.path.join(RX, fn), encoding='utf-8').read()
        body = json.loads(gzip.decompress(base64.b64decode(re.search(r'\{"gz":"([^"]+)"', s).group(1))))
        for r in body:
            code, t, times, source, yl = r[1], r[2], int(r[6] or 1), r[7], r[9]
            p = rxn_code.parse(code)                       # version 2 (an old ~ code reads the same)
            if not p['ok']: continue
            reagent, cond = '.'.join(p['agents']), p['cond']
            cons[t][r[3] or ''] += times; form[t][r[4] or ''] += times
            total += times
            count[t] += times
            recipes[t][reagent] += times
            solv[t][cond] += times
            try: y = float(yl)
            except (TypeError, ValueError): y = None
            if y is not None and 0 < y <= 100:
                if len(yields[t]) < MAX_YIELDS: yields[t].append(y)
                k = (t, reagent)
                if k not in example or y > example[k][0]: example[k] = (y, source)
    os.makedirs(OUT, exist_ok=True)
    for old in os.listdir(OUT): os.remove(os.path.join(OUT, old))
    shards = collections.defaultdict(dict)
    kept = 0
    for t, n in count.items():
        if n < MIN_N: continue
        rec = []
        for reagent, c in recipes[t].most_common(TOP_RECIPES):
            ex = example.get((t, reagent))
            rec.append([reagent, c, ex[1] if ex else '', ex[0] if ex else None])
        sol = [[k, c] for k, c in solv[t].most_common(TOP_SOLVENTS)]          # [solvents as written, count]
        ys = yields[t]
        shards[t // PER_SHARD][str(t)] = [n, round(statistics.median(ys), 1) if ys else None, rec, sol, cons[t].most_common(1)[0][0], form[t].most_common(1)[0][0]]
        kept += 1
    for k, d in shards.items():
        raw = json.dumps(d, separators=(',', ':'), ensure_ascii=False).encode('utf-8')
        b64 = base64.b64encode(gzip.compress(raw, 9)).decode('ascii')
        with open(os.path.join(OUT, f'{k:03d}.js'), 'w', encoding='utf-8') as fh:
            fh.write(f'(window.LEARNED_PARTS=window.LEARNED_PARTS||{{}})["c/{k:03d}"]={{"gz":"{b64}"}};\n')
    size = sum(os.path.getsize(os.path.join(OUT, f)) for f in os.listdir(OUT))
    print(f'{total:,} recorded reactions, {len(count):,} templates, {kept:,} kept (>= {MIN_N} records), {len(shards)} files, {size / 1e6:.1f} MB')


if __name__ == '__main__':
    main()
