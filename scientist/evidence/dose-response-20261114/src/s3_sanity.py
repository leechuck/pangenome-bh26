"""Pilot sanity checks and timing. (i) REPRO_A reproduces the existing arm-A table; (ii) REPRO_B reproduces
the existing arm-B table. If not byte-identical, quantify call differences. Timing per cell from COMPLETE.json.
  python s3_sanity.py MANIFEST TAG
"""
import collections, csv, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import s3_lib as L
C = L.C
man = C.load_json(sys.argv[1]); tag = sys.argv[2]
runs = L.donor_runs()
rows, timing = [], []
def cmp_tables(a, b):
    ra, rb = C.read_t1k_raw(a), C.read_t1k_raw(b)
    ca, cb = C.read_t1k(a), C.read_t1k(b)
    diff_genes = [g for g in sorted(set(ca) | set(cb)) if sorted(map(sorted, ca.get(g, []))) != sorted(map(sorted, cb.get(g, [])))]
    ab = [g for g in sorted(set(ra) | set(rb)) if ra.get(g) != rb.get(g)]
    return diff_genes, ab
for t in man['tasks']:
    d = t['donor']
    for cell in t['cells']:
        w = L.STAGE / 'work' / d / cell
        if (w / 'COMPLETE.json').exists():
            rec = C.load_json(w / 'COMPLETE.json')
            timing.append(dict(donor=d, readset=runs[d]['readset'], cell=cell, t1k_seconds=rec['t1k_seconds'],
                               n_added=rec['n_added'], n_drawn=rec.get('n_drawn', ''), status='complete'))
            if cell in ('REPRO_A', 'REPRO_B'):
                orig = (Path(runs[d]['armA_run']) if cell == 'REPRO_A' else Path(L.armB_spec(runs[d]['readset'], d)['run'])) / 'calls_genotype.tsv'
                same = rec['identical']
                dg, ab = ([], []) if same else cmp_tables(orig, w / 'calls_genotype.tsv')
                rows.append(dict(donor=d, check=cell, byte_identical=same, genes_with_different_calls=len(dg),
                                 genes_with_different_rows_incl_abundance=len(ab), different_genes=','.join(dg)))
        else:
            timing.append(dict(donor=d, readset=runs[d]['readset'], cell=cell, t1k_seconds='', n_added='', n_drawn='',
                               status='FAILED' if (w / 'FAILURE.json').exists() else 'missing'))
R = L.STAGE / 'results'
for name, data in (('sanity', rows), ('timing', timing)):
    if data:
        c = list(data[0])
        with open(R / f'pilot_{name}_{tag}.tsv', 'w') as f:
            f.write('\t'.join(c) + '\n')
            for r in data: f.write('\t'.join(str(r[x]) for x in c) + '\n')
print(json.dumps(rows, indent=1)); print(sum(1 for t in timing if t['status']!='complete'), 'not complete of', len(timing))
