"""Make manifests (local, small): python s3_manifest.py pilot|full OUT.json   -- reads results/cells.tsv, donors.tsv"""
import csv, hashlib, json, sys
from pathlib import Path
R = Path(__file__).resolve().parent.parent / 'results'
donors = list(csv.DictReader(open(R / 'donors.tsv'), delimiter='\t'))
cells = {}
for r in csv.DictReader(open(R / 'cells.tsv'), delimiter='\t'):
    cells.setdefault(r['donor'], []).append(r['cell'])
mode, out = sys.argv[1], sys.argv[2]
drop = set(sys.argv[3].split(',')) if len(sys.argv) > 3 else set()
h = lambda d: hashlib.sha256((d + '|pilot').encode()).hexdigest()
if mode == 'pilot':
    sel = []
    for g in ('ME', 'EAS', 'SAS'):
        sel += sorted([d['donor'] for d in donors if d['group'] == g], key=h)[:2]
else:
    sel = [d['donor'] for d in donors]
tasks = [dict(donor=d, cells=['REPRO_A', 'REPRO_B'] * (mode == 'pilot') + [c for c in cells[d] if c not in drop]) for d in sel]
json.dump(dict(mode=mode, tasks=tasks), open(out, 'w'), indent=1)
print([(t['donor'], len(t['cells'])) for t in tasks])
