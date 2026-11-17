"""Steps 1 and 2: donor table, panel draws for every (donor, cell), availability, truth-only supply.
No typing, no reads. Writes tables (no sequences) to results/ and DRAWS.json to work/.
"""
import collections
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s3_lib as L  # noqa: E402

C = L.C
R = L.STAGE / 'results'
R.mkdir(exist_ok=True)


def tsv(path, rows, cols):
    with open(path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join('' if r.get(c) is None else str(r.get(c)) for c in cols) + '\n')


def main():
    folds, pop, grp, comp_of = L.load_cohort()
    runs = L.donor_runs()
    records = C.load_records()
    by_donor = collections.defaultdict(list)
    for r in records:
        by_donor[r['donor']].append(r)
    missing = [d for d in pop if d not in by_donor]
    assert not missing, missing
    comps = L.components_by_group(grp, comp_of)
    comp_size = collections.Counter(comp_of[d] for d in pop)
    nonasian = L.load_nonasian()
    assert not (set(nonasian) & set(pop)), 'non-Asian donors overlap the benchmark donors'
    na_comps = {d: [d] for d in nonasian}
    na_complete = sum(1 for d, v in nonasian.items() if len(v) == 16)
    ipd = L.ipd_sequence_set()

    # ---- step 1: donor table
    drows = []
    for d in sorted(pop):
        absent = collections.Counter(r['gene'] for r in by_donor[d] if not r['genomic_labels'])
        x = runs[d]
        drow = dict(donor=d, population=pop[d], group=grp[d], readtype='sim' if x['readset'] == 'T-sim' else 'real',
                    readset=x['readset'], component=comp_of[d], component_size=comp_size[comp_of[d]],
                    n_copies=len(by_donor[d]), n_ipd_absent=sum(absent.values()),
                    r1=x.get('r1'), r2=x.get('r2'), armA_ok=x.get('armA_run') is not None)
        for g in L.GENES:
            drow['absent_' + g] = absent.get(g, 0)
        try:
            b = L.armB_spec(x['readset'], d)
            drow['armB_ok'] = True
        except Exception as e:  # noqa: BLE001
            drow['armB_ok'] = False
        drows.append(dict(drow))
    dcols = ['donor', 'population', 'group', 'readtype', 'readset', 'component', 'component_size', 'n_copies',
             'n_ipd_absent'] + ['absent_' + g for g in L.GENES] + ['armA_ok', 'armB_ok', 'r1', 'r2']
    tsv(R / 'donors.tsv', drows, dcols)

    # ---- draws and panels
    draws = {}
    supply, cellrows = [], []
    panel_cache = {}
    for dr in drows:
        d, g = dr['donor'], dr['group']
        draws[d] = {}
        for cell, src, n in L.cell_names_for(g, with_nonasian=True):
            if src == L.NONASIAN:
                drawn, used = L.draw_panel(d, cell, None, na_comps, n)
                recs_of = lambda x: nonasian[x]  # noqa: E731
            else:
                drawn, used = L.draw_panel(d, cell, comp_of[d], comps[src], n)
                recs_of = lambda x: by_donor[x]  # noqa: E731
            assert d not in drawn
            assert all(comp_of.get(x, x) != comp_of[d] for x in drawn if src != L.NONASIAN)
            panel_recs = [r for x in drawn for r in recs_of(x)]
            adds, st = L.build_additions(panel_recs, ipd)
            sset = collections.defaultdict(set)
            bset = collections.defaultdict(set)
            for a in adds:
                sset[a['gene']].add(a['sha'])
                bset[a['gene']].add(a['body_sha'])
            draws[d][cell] = dict(source=src, n_requested=n, n_drawn=len(drawn), drawn=drawn,
                                  components=used, n_added=st['added'], stats=st)
            cellrows.append(dict(donor=d, group=g, readtype=dr['readtype'], readset=dr['readset'], cell=cell,
                                 source=src, n_requested=n, n_drawn=len(drawn), n_components=len(used),
                                 n_added=st['added'], skipped_ipd=st['skipped_ipd'], skipped_N=st['skipped_N'],
                                 short=int(len(drawn) < n)))
            for r in by_donor[d]:
                seq = r['core_sequence'].upper()
                ex = r['core_exons']
                absent = not r['genomic_labels']
                supply.append(dict(donor=d, group=g, readtype=dr['readtype'], readset=dr['readset'], cell=cell,
                                   source=src, n_requested=n, n_drawn=len(drawn), gene=r['gene'],
                                   copy=r['haplotype'], ipd_absent=int(absent),
                                   supplied_byte=int(absent and L.sha(seq) in sset[r['gene']]),
                                   supplied_body=int(absent and L.sha(L.gene_body(seq, ex)) in bset[r['gene']])))
    tsv(R / 'cells.tsv', cellrows, ['donor', 'group', 'readtype', 'readset', 'cell', 'source', 'n_requested',
                                    'n_drawn', 'n_components', 'n_added', 'skipped_ipd', 'skipped_N', 'short'])
    tsv(R / 'supply_rows.tsv', supply, ['donor', 'group', 'readtype', 'readset', 'cell', 'source', 'n_requested',
                                        'n_drawn', 'gene', 'copy', 'ipd_absent', 'supplied_byte', 'supplied_body'])
    C.save_json(L.STAGE / 'work/DRAWS.json', draws)

    # ---- summaries
    agg = collections.defaultdict(lambda: collections.Counter())
    for s in supply:
        a = agg[(s['group'], s['readtype'], s['cell'])]
        a['copies'] += 1
        a['absent'] += s['ipd_absent']
        a['sup_byte'] += s['supplied_byte']
        a['sup_body'] += s['supplied_body']
    donors_n = collections.Counter()
    drawn_sum = collections.Counter()
    short = collections.Counter()
    for c in cellrows:
        k = (c['group'], c['readtype'], c['cell'])
        donors_n[k] += 1
        drawn_sum[k] += c['n_drawn']
        short[k] += c['short']
    srows = []
    for k in sorted(agg):
        a = agg[k]
        srows.append(dict(group=k[0], readtype=k[1], cell=k[2], test_donors=donors_n[k],
                          mean_n_drawn=round(drawn_sum[k] / donors_n[k], 1), donors_short=short[k],
                          ipd_absent_copies=a['absent'], supplied_byte=a['sup_byte'],
                          supplied_byte_share=round(a['sup_byte'] / a['absent'], 4) if a['absent'] else None,
                          supplied_body=a['sup_body'],
                          supplied_body_share=round(a['sup_body'] / a['absent'], 4) if a['absent'] else None))
    tsv(R / 'truth_ceiling_summary.tsv', srows, list(srows[0]))
    info = dict(donors=len(pop), by_group=dict(collections.Counter(grp.values())),
                by_readset=dict(collections.Counter(x['readset'] for x in runs.values())),
                armA_missing=[d for d, x in runs.items() if x.get('armA_run') is None],
                armB_missing=[r['donor'] for r in drows if not r['armB_ok']],
                nonasian_donors=len(nonasian), nonasian_donors_with_16_copies=na_complete,
                n_cell_rows=len(cellrows), n_cells_short=sum(c['short'] for c in cellrows),
                components_by_group={g: len(v) for g, v in comps.items()})
    C.save_json(R / 'inventory_info.json', info)
    print(json.dumps(info, indent=1))


if __name__ == '__main__':
    main()
