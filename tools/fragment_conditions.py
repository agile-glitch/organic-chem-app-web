"""What the recorded reactions say about each fragment: its conditions and the sites it competes with.

    python tools/fragment_conditions.py                       # -> tools/data/fragment_conditions.json
    python tools/fragment_conditions.py --sample 2500         # reactions sampled per fragment for the competition count

Source: the learned reactions in data/learned/rx (US patents 1976-2016, text-mined by D. M. Lowe, CC0; Open Reaction Database,
CC-BY-SA 4.0): one row per distinct reaction, the starting materials, the reagent, the recorded solvents, the product, the
fragments CONSUMED (their count drops from the starting materials to the product) and FORMED, and the recorded yield.
What is NOT in that data, and therefore not reported: temperature, concentration, time, equivalents.

Per fragment, over the reactions that consumed it:
  * how many, the median recorded yield;
  * the commonest solvents and the share of protic / polar aprotic / nonpolar solvent systems (Solvents sheet of reaction_rules.xlsx;
    a mixture counts as protic if any solvent is protic, else polar aprotic if any is);
  * the commonest bases, acids and Lewis acids, metal catalysts and other reagents (roles from tools/data/reagents.json).

Competing sites, in the same molecule: for fragments f and g that sit in one starting material of a recorded reaction,
  s(f->g) = share of the reactions that consumed f in which g, present in the same molecule, was NOT consumed
          = how often g survives while f reacts.
  f is said to react in preference to g when g survives in >= 60% of the (>= 15) reactions that consumed f and f survives in at
  least 30 points fewer of the reactions that consumed g; when g never gets a turn (fewer than 15 such reactions) it takes >= 90% of >= 30.
It is a count of what was recorded: the chemist chose reagents that favour one site, so it measures chemoselectivity under the
conditions people used, not an intrinsic rate. Reactions are sampled (--sample per fragment) and arene C-H is left out of the lists.
"""
import argparse, base64, collections, glob, gzip, json, os, random, re, statistics, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
APP = os.path.dirname(HERE)
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, 'data', 'fragment_conditions.json')
SKIP_ROLES = {'counter-ion / spectator ion', 'solvent', 'gas', 'salt (work-up or additive)'}
GROUPS = [('bases', {'base'}), ('acids', {'acid', 'Lewis acid'}), ('catalysts', {'catalyst (metal)', 'oxidizing catalyst', 'phase-transfer catalyst'}),
          ('other', None)]


def load_rows():
    rows = []
    for f in sorted(glob.glob(os.path.join(APP, 'data', 'learned', 'rx', '*.js'))):
        s = open(f, encoding='utf-8').read(); mm = re.search(r'"gz":"([^"]+)"', s)
        for r in json.loads(gzip.decompress(base64.b64decode(mm.group(1)))):
            parts = r[1].split('~')
            rows.append((parts[0], parts[1], parts[2], [x for x in r[3].split(' + ') if x], [x for x in r[4].split(' + ') if x], r[9]))
    return rows


def pct(n, d):
    return f'{100 * n / d:.0f}%' if d else '-'


def top(counter, n, total, names=None):
    return '; '.join(f'{(names or {}).get(k, k)} {pct(v, total)}' for k, v in counter.most_common(n))


def conditions(rows, names):
    import fragment_properties as FP
    solv = {str(r['Solvent']): str(r['Type (protic / polar aprotic / nonpolar)']).strip() for r in FP.sheet('reaction_rules.xlsx', 'Solvents')}
    reag = {}
    for r in json.load(open(os.path.join(HERE, 'data', 'reagents.json'), encoding='utf-8')):
        reag[r['smiles']] = (r.get('abbreviation') or r.get('name') or r['smiles'], r['role'])
    acc = collections.defaultdict(lambda: dict(n=0, formed=0, yields=[], solvent=collections.Counter(), with_solvent=0, types=collections.Counter(),
                                               reag={g: collections.Counter() for g, _ in GROUPS}))
    formed = collections.Counter()
    for sm, reagent, cond, consumed, made, y in rows:
        for fid in made: formed[fid] += 1
        if not consumed: continue
        sv = [x for x in cond.split(', ') if x]
        kinds = {solv.get(x) for x in sv if solv.get(x)}
        typ = 'protic' if 'protic' in kinds else 'polar aprotic' if 'polar aprotic' in kinds else 'nonpolar' if kinds else None
        agents = {}
        for c in reagent.split('.'):
            if c in reag and reag[c][1] not in SKIP_ROLES: agents[c] = reag[c]
        for fid in consumed:
            a = acc[fid]; a['n'] += 1
            if isinstance(y, (int, float)): a['yields'].append(y)       # a missing yield is ''
            if sv:
                a['with_solvent'] += 1
                for x in set(sv): a['solvent'][x] += 1
                if typ: a['types'][typ] += 1
            for c, (nm, role) in agents.items():
                g = next((g for g, rs in GROUPS if rs and role in rs), 'other')
                a['reag'][g][nm] += 1
    out = {}
    for fid in names:
        a = acc.get(fid)
        if not a: out[fid] = dict(n_consumed=0, n_formed=formed.get(fid, 0)); continue
        n, ws = a['n'], a['with_solvent']
        out[fid] = dict(n_consumed=n, n_formed=formed.get(fid, 0),
                        median_yield=round(statistics.median(a['yields']), 1) if a['yields'] else None,
                        solvents=top(a['solvent'], 5, ws), solvent_recorded=pct(ws, n),
                        solvent_type=', '.join(f'{k} {pct(v, sum(a["types"].values()))}' for k, v in a['types'].most_common()) or None,
                        **{g: top(a['reag'][g], 5 if g == 'other' else 4, n) or None for g, _ in GROUPS})
    return out


def competition(rows, frags, sample, log=print):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog('rdApp.*')
    ids = [f['Fragment ID'] for f in frags]
    q = {f['Fragment ID']: Chem.MolFromSmarts(f['SMARTS pattern']) for f in frags}
    arene = next(f['Fragment ID'] for f in frags if str(f['Fragment name']).startswith('arene C–H'))
    memo = {}

    def present(smi):
        v = memo.get(smi)
        if v is None:
            m = Chem.MolFromSmiles(smi)
            v = frozenset(i for i in ids if m is not None and m.HasSubstructMatch(q[i]))
            memo[smi] = v
        return v
    by = collections.defaultdict(list)
    for k, r in enumerate(rows):
        for fid in r[3]: by[fid].append(k)
    wins = collections.defaultdict(int); both = collections.defaultdict(int)       # (g, f): g consumed; f survived / f consumed too
    t0 = time.time()
    for n, g in enumerate(ids):
        ks = by.get(g, [])
        ks = random.Random(g).sample(ks, sample) if len(ks) > sample else ks
        for k in ks:
            sm, _, _, consumed, _, _ = rows[k]; used = set(consumed)
            for smi in set(sm.split('.')):
                p = present(smi)
                if g not in p: continue
                for f in p:
                    if f == g: continue
                    (both if f in used else wins)[(g, f)] += 1
        if n % 10 == 9: log(f'  competition: {n + 1}/{len(ids)} fragments, {len(memo):,} molecules read ({time.time() - t0:.0f} s)')
    name = {f['Fragment ID']: str(f['Fragment name']).split(' (')[0] for f in frags}
    res = {}
    for f in ids:
        beats, loses = [], []
        for g in ids:
            if g == f or g == arene: continue
            nfg = wins[(f, g)] + both[(f, g)]; ngf = wins[(g, f)] + both[(g, f)]
            sfg = wins[(f, g)] / nfg if nfg >= 15 else None            # g survives while f reacts
            sgf = wins[(g, f)] / ngf if ngf >= 15 else None            # f survives while g reacts
            if sfg is not None and sgf is not None:
                if sfg >= 0.6 and sfg - sgf >= 0.3: beats.append((sfg - sgf, nfg, g, sfg, nfg, sgf, ngf))
                elif sgf >= 0.6 and sgf - sfg >= 0.3: loses.append((sgf - sfg, ngf, g, sgf, ngf, sfg, nfg))
            elif sfg is not None and sfg >= 0.9 and nfg >= 30 and ngf < 15: beats.append((sfg, nfg, g, sfg, nfg, None, ngf))
            elif sgf is not None and sgf >= 0.9 and ngf >= 30 and nfg < 15: loses.append((sgf, ngf, g, sgf, ngf, None, nfg))
        fmt = lambda t: f'{name[t[2]]} [{100 * t[3]:.0f}% of {t[4]}' + (f' vs {100 * t[5]:.0f}% of {t[6]}]' if t[5] is not None else ']')
        beats.sort(reverse=True); loses.sort(reverse=True)
        res[f] = dict(reacts_before=('; '.join(fmt(t) for t in beats[:4]) or None), gives_way_to=('; '.join(fmt(t) for t in loses[:4]) or None),
                      n_beats=len(beats), n_loses=len(loses))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sample', type=int, default=2500); ap.add_argument('--out', default=OUT)
    a = ap.parse_args(); t0 = time.time()
    import fragment_properties as FP
    frags = FP.fragments()
    names = {f['Fragment ID']: str(f['Fragment name']) for f in frags}
    rows = load_rows(); print(f'{len(rows):,} recorded reactions read ({time.time() - t0:.0f} s)', flush=True)
    cond = conditions(rows, names); print(f'conditions counted ({time.time() - t0:.0f} s)', flush=True)
    comp = competition(rows, frags, a.sample, log=lambda s: print(s, flush=True))
    out = {fid: {**cond[fid], **comp[fid]} for fid in names}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(out, open(a.out, 'w', encoding='utf-8'), indent=1, ensure_ascii=False)
    print(f'wrote {a.out} ({time.time() - t0:.0f} s)')


if __name__ == '__main__':
    main()
