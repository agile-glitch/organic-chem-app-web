#!/usr/bin/env python3
"""What changed between two bench reports: python tools/namer_bench/compare.py out/A.json out/B.json"""
import json, sys
if hasattr(sys.stdout, 'reconfigure'): sys.stdout.reconfigure(encoding='utf-8')
a = {r['ref_smiles']: r for r in json.load(open(sys.argv[1], encoding='utf-8'))}
b = {r['ref_smiles']: r for r in json.load(open(sys.argv[2], encoding='utf-8'))}
GOOD = {'ok': 4, 'stereo': 3, 'unnamed': 2, 'unreadable': 1, 'wrong': 0}
up = down = 0
for k, rb in b.items():
    ra = a.get(k)
    if not ra or ra['kind'] == rb['kind']: continue
    better = GOOD[rb['kind']] > GOOD[ra['kind']]
    up += better; down += not better
    print('BETTER' if better else 'WORSE ', ra['kind'], '->', rb['kind'], '|', k, '|', ra.get('name') or ra.get('reason'), '=>', rb.get('name') or rb.get('reason'))
print(f'{up} better, {down} worse')
