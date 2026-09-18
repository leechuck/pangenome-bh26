"""Score graph recruitment development pilots; never an independent validation."""
import argparse
import csv
import json
from pathlib import Path
from score_linear_controls import (HERE,SERIES,REPO,GENES,Nomenclature,truth_a,
    truth_slots,genomic_truth,compare,prediction,decode,sha,write_tsv)


CONDITIONS=('native','all','hprc/sampled','hprc/unsampled','hprc_asian/sampled','hprc_asian/unsampled')


def load_control(folder,donor,expected_reads):
    manifest=folder/'manifest.json'
    if not manifest.exists():return 'not_started',{},None
    report=json.loads(manifest.read_text())
    if report['donor']!=donor or report['reads_sha256']!=expected_reads:
        raise ValueError('Different donor or read inputs')
    state=report['status']
    if state!='complete':return state,{},report
    if json.loads((folder/'COMPLETE.json').read_text())!=report:
        raise ValueError('Completion record mismatch')
    # Only small outputs are transferred. FASTQ hashes remain in run provenance.
    for name in ('calls.json',donor+'_genotype.tsv'):
        if sha(folder/name)!=report['output_sha256'][name]:raise ValueError('Changed genotype output')
    calls=json.loads((folder/'calls.json').read_text())
    if calls!=decode((folder/(donor+'_genotype.tsv')).read_text(),{}):
        raise ValueError('Decoded calls disagree with genotype table')
    return state,calls,report


def run(snapshot,output,donor,conditions=CONDITIONS,quiet=False):
    donors={r['donor']:r for r in csv.DictReader((HERE/'execution-cohort.tsv').open(),delimiter='\t')}
    if donor not in donors:raise ValueError('Only exposed development donors allowed')
    base=SERIES/'snapshot/t1k4'/donor
    bm=json.loads((base/'manifest.json').read_text())
    if bm['status']!='complete':raise ValueError('Baseline incomplete')
    reads=bm['configuration']['reads'];table=base/(donor+'_genotype.tsv')
    completion=json.loads((base/'COMPLETE').read_text())
    if sha(table)!=completion['outputs'][table.name]:raise ValueError('Baseline genotype changed')
    controls={'T1K':('complete',decode(table.read_text(),{}),None)}
    controls['ipd_genome']=load_control(HERE/'work/linear-snapshot/ipd_genome'/donor,donor,reads)
    for condition in conditions:
        controls[condition]=load_control(snapshot/donor/condition,donor,reads)
    catalogue=REPO/'hla-analysis/results/sequence_catalogue.tsv'
    truth=truth_a(catalogue,{donor});nom=Nomenclature();rows=[];summaries=[]
    for method,(state,calls,report) in controls.items():
        for fields in (2,4):
            subset=[]
            for gene in GENES:
                records=truth.get((donor,gene))
                slots=(genomic_truth(records) if fields==4 else truth_slots(records,'two_field',nom,gene)) if records else None
                called,correct,matches=compare(prediction(calls,gene,fields),slots) if slots is not None else (0,0,0)
                subset.append(dict(donor=donor,method=method,state=state,gene=gene,fields=fields,
                                   eligible=int(slots is not None),called=called,correct=correct,allele_matches=matches))
            rows.extend(subset)
            summaries.append(dict(donor=donor,method=method,state=state,fields=fields,
                                  **{key:sum(r[key] for r in subset) for key in ('eligible','called','correct','allele_matches')}))
    output.mkdir(parents=True,exist_ok=True)
    write_tsv(output/'gene_scores.tsv',rows);write_tsv(output/'summary.tsv',summaries)
    provenance=dict(scope='One exposed development donor; no independent accuracy claim',donor=donor,
                    baseline_genotype_sha256=sha(table),catalogue_sha256=sha(catalogue),
                    scorer_sha256=sha(Path(__file__)),
                    methods={method:report for method,(_,_,report) in controls.items() if report})
    (output/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    if not quiet:print(json.dumps(summaries,indent=2))
    return rows,summaries


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot',type=Path,default=HERE/'work/graph-recruitment-snapshot')
    p.add_argument('--output',type=Path,default=HERE/'development/graph-recruitment-pilot')
    p.add_argument('--donor',default='HG00658')
    a=p.parse_args();run(a.snapshot,a.output,a.donor)
