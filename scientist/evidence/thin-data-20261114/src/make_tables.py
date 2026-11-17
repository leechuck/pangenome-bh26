"""Local: markdown tables from results/figure_table.tsv and contrasts.tsv (small summary tables only)."""
import csv, sys, collections
R = '/mnt/data1/dogohla/pangenome/scientist/evidence/thin-data-20261114/results/'
D = ['full', '1/2', '1/4', '1/8', '1/16']
fig = {}
for r in csv.DictReader(open(R + 'figure_table.tsv'), delimiter='\t'):
    fig[(r['read_type'], r['group'], r['depth'], r['arm'], r['rule'], r['metric'])] = r
def cell(k, pct=True):
    r = fig.get(k)
    if not r: return 'NA'
    return '%.1f [%.1f-%.1f]' % (float(r['value_pct']), float(r['lo95']), float(r['hi95']))
out = []
for rt, group, arms in (('real', 'REAL', 'ABN'), ('real', 'EAS', 'ABN'), ('real', 'SAS', 'ABN'), ('sim', 'ME', 'ABMN')):
    for metric in ('acc2', 'acc3', 'acc4'):
        out.append(f'\n**{rt} reads, group {group}, {metric} (percent, 95% interval)**\n')
        out.append('| arm / rule | ' + ' | '.join(D) + ' |')
        out.append('|---|' + '---|' * len(D))
        for arm in arms:
            for rule in ('P0f', 'P0s', 'PM', 'PO'):
                row = [cell((rt, group, d, arm, rule, metric)) for d in D]
                if all(x == 'NA' for x in row): continue
                out.append(f'| {arm} / {rule} | ' + ' | '.join(row) + ' |')
open(R + 'tables_accuracy.md', 'w').write('\n'.join(out) + '\n')
print('\n'.join(out)[:3000])
