"""Score the full exposed-cohort ablation; never claim independent validation."""
import csv
import json
from pathlib import Path
from score_graph_pair_refinement import load
from score_linear_controls import HERE,REPO,GENES,truth_a,genomic_truth,truth_slots,Nomenclature,prediction,compare,sha,write_tsv


def run():
    plan_path=HERE/'alignment-only-ablation-plan.json'
    plan=json.loads(plan_path.read_text())
    out=HERE/'development/completed-alignment-only-ablation-v1'
    if out.exists():raise FileExistsError(out)
    snapshot=HERE/'work/alignment-only-ablation-v1'
    batches={}
    p=json.loads((HERE/'graph-development-batch.json').read_text())
    batches['development']=[p['existing_pilot'],*p['donors']]
    batches['exposed_validation']=[dict(r,stratum=r['superpopulation']) for r in
        csv.DictReader((HERE/'validation/cohort.tsv').open(),delimiter='\t')]
    expected={(group,r['donor'],panel) for group,donors in batches.items()
              for r in donors for panel in ('hprc','hprc_asian')}
    assert {(r['cohort'],r['donor'],r['panel']) for r in plan['cases']}==expected
    assert len(plan['cases'])==len(expected)
    calls_by_case={};provenance={}
    for case in plan['cases']:
        group,d,panel=(case[k] for k in ('cohort','donor','panel'))
        baseline=(HERE/'work/linear-snapshot/ipd_genome'/d if group=='development'
                  else HERE/'work/validation-snapshot/genomic-ipd-v1'/d)
        folder=snapshot/group/d/panel/'unsampled'
        state,calls,record=load(folder,baseline)
        assert state=='complete'
        assert record['driver_sha256']==plan['code_sha256']['run_graph_pair_alignment_only.py']
        assert record['model_sha256']==plan['code_sha256']['graph_pair_refinement_alignment_only.py']
        assert record['base_model_sha256']==plan['code_sha256']['graph_pair_refinement.py']
        assert record['parameters']['alignment_only'] and record['parameters']['internal_pairs'] and record['parameters']['pair_specific']
        calls_by_case[group,d,panel]=calls
        provenance['/'.join((group,d,panel))]=sha(folder/'COMPLETE.json')
    truth=truth_a(REPO/'hla-analysis/results/sequence_catalogue.tsv',
                  {d for _,d,_ in expected})
    nom=Nomenclature();rows=[]
    sources={'development':HERE/'development/completed-reference-controls-v2/gene_scores.tsv',
             'exposed_validation':HERE/'validation/completed-frozen-v1/gene_scores.tsv'}
    for group,donors in batches.items():
        ids={d['donor'] for d in donors}
        for r in csv.DictReader(sources[group].open(),delimiter='\t'):
            if r['donor'] in ids and r['method'] in ('T1K','ipd_genome'):
                rows.append(dict({k:int(v) if k in ('fields','eligible','called','correct','allele_matches') else v
                                  for k,v in r.items()},cohort=group))
        for donor in donors:
            d=donor['donor']
            for panel in ('hprc','hprc_asian'):
                calls=calls_by_case[group,d,panel]
                for gene in GENES:
                    records=truth.get((d,gene))
                    for fields in (2,4):
                        slots=(genomic_truth(records) if fields==4 else truth_slots(records,'two_field',nom,gene)) if records else None
                        called,correct,matches=compare(prediction(calls,gene,fields),slots) if slots is not None else (0,0,0)
                        rows.append(dict(cohort=group,donor=d,family=donor['family'],stratum=donor['stratum'],
                            method='alignment_only_'+panel,gene=gene,fields=fields,state='complete',
                            eligible=int(slots is not None),called=called,correct=correct,allele_matches=matches))
    summary=[]
    for group,donors in batches.items():
        for method in ('T1K','ipd_genome','alignment_only_hprc','alignment_only_hprc_asian'):
            for fields in (2,4):
                subset=[r for r in rows if (r['cohort'],r['method'],r['fields'])==(group,method,fields)]
                assert len(subset)==len(donors)*len(GENES)
                summary.append(dict(cohort=group,method=method,fields=fields,donors=len(donors),
                    **{k:sum(r[k] for r in subset) for k in ('eligible','called','correct','allele_matches')}))
    out.mkdir(parents=True)
    write_tsv(out/'gene_scores.tsv',rows);write_tsv(out/'summary.tsv',summary)
    (out/'provenance.json').write_text(json.dumps(dict(scope='Post-unblinding exploratory ablation; not independent validation',
        plan_sha256=sha(plan_path),scorer_sha256=sha(Path(__file__)),runs=provenance),indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':run()
