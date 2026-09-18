"""Score graph-refined calls on the predefined, exposed development panel only."""
import csv
import json
from pathlib import Path
from score_linear_controls import HERE,REPO,GENES,truth_a,genomic_truth,truth_slots,Nomenclature,prediction,compare,sha,write_tsv


CONDITIONS=('hprc/sampled','hprc/unsampled','hprc_asian/sampled','hprc_asian/unsampled')


def load(folder,baseline):
    if not (folder/'manifest.json').exists():return 'not_started',{},None
    report=json.loads((folder/'manifest.json').read_text())
    bm=json.loads((baseline/'COMPLETE.json').read_text())
    if report['donor']!=bm['donor'] or report['reads_sha256']!=bm['reads_sha256'] or report['baseline_sha256']!=sha(baseline/'COMPLETE.json'):
        raise ValueError('Refinement baseline mismatch')
    if report['status']!='complete':return report['status'],{},report
    if json.loads((folder/'COMPLETE.json').read_text())!=report:raise ValueError('Completion mismatch')
    for name in ('calls.json','decisions.json'):
        if sha(folder/name)!=report['output_sha256'][name]:raise ValueError('Refined output changed')
    if sha(baseline/'calls.json')!=bm['output_sha256']['calls.json']:raise ValueError('Baseline calls changed')
    native=json.loads((baseline/'calls.json').read_text());calls=json.loads((folder/'calls.json').read_text())
    decisions=json.loads((folder/'decisions.json').read_text())
    if set(native)!=set(calls):raise ValueError('Refinement gene set changed')
    if set(decisions)!=set(GENES):raise ValueError('Missing or unexpected refinement decisions')
    for gene in native:
        if native[gene]['resolutions']['2']!=calls[gene]['resolutions']['2']:
            raise ValueError('Refinement changed two-field calls')
        if not decisions.get(gene,{'changed':False})['changed'] and native[gene]!=calls[gene]:
            raise ValueError('Unrecorded call change')
    return report['status'],calls,report


def run(snapshot=None,output=None):
    snapshot=Path(snapshot) if snapshot is not None else HERE/'work/graph-pair-refinement-snapshot'
    plan=json.loads((HERE/'graph-development-batch.json').read_text())
    donors=[plan['existing_pilot'],*plan['donors']]
    exposed={r['donor'] for r in csv.DictReader((HERE/'execution-cohort.tsv').open(),delimiter='\t')}
    if not {r['donor'] for r in donors}<=exposed:raise ValueError('Only exposed development donors permitted')
    truth=truth_a(REPO/'hla-analysis/results/sequence_catalogue.tsv',{r['donor'] for r in donors})
    nom=Nomenclature();rows=[];provenance={}
    original=list(csv.DictReader((HERE/'development/completed-reference-controls-v2/gene_scores.tsv').open(),delimiter='\t'))
    for row in original:
        if row['donor'] in {r['donor'] for r in donors} and row['method'] in ('T1K','ipd_genome'):
            rows.append({k:int(v) if k in ('fields','eligible','called','correct','allele_matches') else v for k,v in row.items()})
    for donor in donors:
        d=donor['donor'];baseline=HERE/'work/linear-snapshot/ipd_genome'/d
        for condition in CONDITIONS:
            state,calls,report=load(snapshot/d/condition,baseline)
            if report:provenance[d+'/'+condition]=report
            for gene in GENES:
                records=truth.get((d,gene))
                for fields in (2,4):
                    slots=(genomic_truth(records) if fields==4 else truth_slots(records,'two_field',nom,gene)) if records else None
                    called,correct,matches=compare(prediction(calls,gene,fields),slots) if slots is not None else (0,0,0)
                    rows.append(dict(donor=d,stratum=donor['stratum'],family=donor['family'],method=condition,gene=gene,
                                     fields=fields,state=state,eligible=int(slots is not None),called=called,correct=correct,allele_matches=matches))
    summary=[]
    for method in ('T1K','ipd_genome',*CONDITIONS):
        for fields in (2,4):
            subset=[r for r in rows if r['method']==method and r['fields']==fields]
            summary.append(dict(method=method,fields=fields,planned_donors=len(donors),
                                completed_donors=len({r['donor'] for r in subset if r['state']=='complete'}),
                                **{k:sum(r[k] for r in subset) for k in ('eligible','called','correct','allele_matches')}))
    out=Path(output) if output is not None else HERE/'development/graph-pair-refinement';out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/'gene_scores.tsv',rows);write_tsv(out/'summary.tsv',summary)
    (out/'provenance.json').write_text(json.dumps(dict(scope='Exposed development only; partial comparisons are not accuracy claims',runs=provenance,scorer_sha256=sha(Path(__file__))),indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':run()
