"""Step 1: donors, read files, pair counts, arm availability, original depth statistics."""
import gzip, json, sys
sys.path.insert(0, '/encrypted/e3008/hohndor-dogohla-20261020/thin-data-20261114/src')
import lib as L
C = L.C

folds = C.load_json(C.S / 'inputs/FOLDS.json')
rows = []
def npairs(p):
    n = 0
    with gzip.open(p, 'rt') as f:
        for _ in f:
            n += 1
    return n // 4
for rs, donors in L.readsets().items():
    for d in donors:
        info = {a: L.arm_info(rs, d, a) for a in L.ARMS_ALL}
        a = info['A']
        row = dict(readset=rs, donor=d, population=folds['donors'][d]['population'],
                   resource=folds['donors'][d]['resource'], r1=str(a['r1']), r2=str(a['r2']),
                   threads_A=a['threads'], seconds_A=round(a['seconds']))
        for k in ('B', 'N', 'M'):
            row[k + '_ok'] = info[k] is not None
            if info[k]:
                row[k + '_same_reads'] = (info[k]['r1'] == a['r1'] and info[k]['r2'] == a['r2'])
                row[k + '_seconds'] = round(info[k]['seconds'])
        row['pairs_r1'] = npairs(a['r1']); row['pairs_r2'] = npairs(a['r2'])
        row.update({k: v for k, v in L.log_stats(a['run'] / 'run.log').items()})
        rows.append(row)
        print(json.dumps(row), flush=True)
out = L.TS / 'results/donor_table.json'
C.save_json(out, rows)
