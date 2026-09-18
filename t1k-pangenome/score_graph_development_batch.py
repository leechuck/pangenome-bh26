"""Score all predefined balanced development donors, preserving pending denominators."""
import json
from score_graph_recruitment import HERE,run,write_tsv


if __name__=='__main__':
    plan=json.loads((HERE/'graph-development-batch.json').read_text())
    donors=[plan['existing_pilot'],*plan['donors']]
    output=HERE/'development/graph-development-batch'
    conditions=('all','hprc/sampled','hprc/unsampled','hprc_asian/sampled','hprc_asian/unsampled')
    rows=[]
    for donor in donors:
        scores,_=run(HERE/'work/graph-recruitment-snapshot',output/'donors'/donor['donor'],
                     donor['donor'],conditions=conditions,quiet=True)
        for row in scores:row.update(stratum=donor['stratum'],family=donor['family'])
        rows.extend(scores)
    summary=[]
    for method in ('T1K','ipd_genome',*conditions):
        for stratum in ('ALL',*sorted({r['stratum'] for r in donors})):
            for fields in (2,4):
                subset=[r for r in rows if r['method']==method and r['fields']==fields and
                        (stratum=='ALL' or r['stratum']==stratum)]
                summary.append(dict(method=method,stratum=stratum,fields=fields,
                                    planned_donors=len({r['donor'] for r in subset}),
                                    completed_donors=len({r['donor'] for r in subset if r['state']=='complete'}),
                                    **{key:sum(r[key] for r in subset) for key in ('eligible','called','correct','allele_matches')}))
    write_tsv(output/'gene_scores.tsv',rows);write_tsv(output/'summary.tsv',summary)
    print(json.dumps([r for r in summary if r['stratum']=='ALL'],indent=2))
