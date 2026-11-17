"""Sanity: my P0s rule at full depth reproduces the frozen scorer (score.score) row by row."""
import csv, json, sys, collections
sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent))
import lib as L
import score as S0
C = L.C
T, rows, avail = S0.score(['D-real', 'T-real', 'T-sim'], ['A', 'B'], 'x')
frozen = {(r['readset'], r['donor'], r['gene'], r['haplotype'], r['arm']): r for r in rows}
mine = {}
with open(L.TS / 'results/rows.tsv') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['rule'] == 'P0s' and r['depth'] == 'full' and r['arm'] in ('A', 'B'):
            mine[(r['readset'], r['donor'], r['gene'], r['haplotype'], r['arm'])] = r
diff = collections.Counter()
ex = []
for k, r in frozen.items():
    m = mine.get(k)
    if m is None:
        diff['missing_in_mine'] += 1
        continue
    for a, b in (('named4', 'ok4'), ('named3_decoded', 'ok3'), ('named2_decoded', 'ok2')):
        if int(bool(r[a])) != int(m[b]):
            diff[a] += 1
            if len(ex) < 20:
                ex.append((k, a, r['call'], m['chosen']))
print(json.dumps(dict(n_frozen=len(frozen), n_mine=len(mine), diff=dict(diff), examples=ex), default=str))
C.save_json(L.TS / 'results/scorer_reproduction.json', dict(n_frozen=len(frozen), n_mine=len(mine), diff=dict(diff), examples=ex))
