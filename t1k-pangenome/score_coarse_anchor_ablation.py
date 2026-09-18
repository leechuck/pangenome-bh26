"""Evaluate truth-free coarse anchoring on exposed donors only."""
import csv
import json
from pathlib import Path
from anchor_t1k_coarse import anchor
from score_graph_pair_refinement import load
from score_linear_controls import HERE,REPO,GENES,truth_a,genomic_truth,truth_slots,Nomenclature,prediction,compare,sha,write_tsv,decode


def run():
    out=HERE/'development/completed-coarse-anchor-ablation-v1'
    if out.exists():raise FileExistsError(out)
    p=json.loads((HERE/'graph-development-batch.json').read_text())
    cohorts={'development':[p['existing_pilot'],*p['donors']],
             'exposed_validation':[dict(r,stratum=r['superpopulation']) for r in
                  csv.DictReader((HERE/'validation/cohort.tsv').open(),delimiter='\t')]}
    plan=json.loads((HERE/'native-locus-ablation-plan.json').read_text())
    calls={};provenance={};decisions={}
    for group,donors in cohorts.items():
        for d in (r['donor'] for r in donors):
            if group=='development':
                original=REPO/'hla-spechla-pg/series20260918/snapshot/t1k4'/d
                marker=original/'COMPLETE';hashes=json.loads(marker.read_text())['outputs']
                baseline=HERE/'work/linear-snapshot/ipd_genome'/d
            else:
                original=HERE/'work/validation-snapshot/frozen-t1k-v1'/d
                marker=original/'COMPLETE.json';hashes=json.loads(marker.read_text())['output_sha256']
                baseline=HERE/'work/validation-snapshot/genomic-ipd-v1'/d
            table=original/(d+'_genotype.tsv')
            assert sha(table)==hashes[table.name]
            assert json.loads((original/'manifest.json').read_text())['status']=='complete'
            native=decode(table.read_text(),{})
            calls[group,d,'T1K']=native
            bm=json.loads((baseline/'COMPLETE.json').read_text())
            assert sha(baseline/'calls.json')==bm['output_sha256']['calls.json']
            candidates={'anchored_ipd_genome':json.loads((baseline/'calls.json').read_text())}
            provenance[group+'/'+d]=dict(original_sha256=sha(marker),genomic_sha256=sha(baseline/'COMPLETE.json'))
            for panel in ('hprc','hprc_asian'):
                folder=HERE/'work/native-locus-ablation-v1'/group/d/panel/'unsampled'
                state,candidate,r=load(folder,baseline)
                assert state=='complete'
                assert r['model_sha256']==plan['code_sha256']['graph_pair_refinement_native_locus.py']
                assert r['driver_sha256']==plan['code_sha256']['run_graph_pair_native_locus.py']
                candidates['anchored_'+panel]=candidate
                provenance[group+'/'+d][panel+'_sha256']=sha(folder/'COMPLETE.json')
            for method,candidate in candidates.items():
                updated,reasons=anchor(native,candidate)
                calls[group,d,method]=updated;decisions[group+'/'+d+'/'+method]=reasons
                for gene in GENES:
                    assert prediction(updated,gene,2)==prediction(native,gene,2)
    truth=truth_a(REPO/'hla-analysis/results/sequence_catalogue.tsv',{d for _,d,_ in calls})
    nom=Nomenclature();rows=[]
    for group,donors in cohorts.items():
        for donor in donors:
            d=donor['donor']
            for method in ('T1K','anchored_ipd_genome','anchored_hprc','anchored_hprc_asian'):
                for gene in GENES:
                    records=truth.get((d,gene))
                    for fields in (2,4):
                        slots=(genomic_truth(records) if fields==4 else truth_slots(records,'two_field',nom,gene)) if records else None
                        called,correct,matches=compare(prediction(calls[group,d,method],gene,fields),slots) if slots is not None else (0,0,0)
                        rows.append(dict(cohort=group,donor=d,family=donor['family'],stratum=donor['stratum'],
                            method=method,gene=gene,fields=fields,state='complete',eligible=int(slots is not None),
                            called=called,correct=correct,allele_matches=matches))
    summary=[]
    for group,donors in cohorts.items():
        for method in ('T1K','anchored_ipd_genome','anchored_hprc','anchored_hprc_asian'):
            for fields in (2,4):
                subset=[r for r in rows if (r['cohort'],r['method'],r['fields'])==(group,method,fields)]
                assert len(subset)==len(donors)*len(GENES)
                summary.append(dict(cohort=group,method=method,fields=fields,donors=len(donors),
                    **{k:sum(r[k] for r in subset) for k in ('eligible','called','correct','allele_matches')}))
    out.mkdir(parents=True)
    write_tsv(out/'gene_scores.tsv',rows);write_tsv(out/'summary.tsv',summary)
    (out/'decisions.json').write_text(json.dumps(decisions,indent=2)+'\n')
    (out/'provenance.json').write_text(json.dumps(dict(scope='Exploratory coarse-anchor composition on exposed donors; not independent validation',
        anchor_sha256=sha(HERE/'anchor_t1k_coarse.py'),scorer_sha256=sha(Path(__file__)),runs=provenance),indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':run()
