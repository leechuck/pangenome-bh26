"""Describe residual genomic-IPD errors using exposed development truth only."""
import collections
import csv
import json
from pathlib import Path
from score_linear_controls import HERE,REPO,truth_a,genomic_truth,write_tsv,sha
from score_comparators import four_field


def run():
    cohort=list(csv.DictReader((HERE/'execution-cohort.tsv').open(),delimiter='\t'))
    source=HERE/'development/completed-reference-controls-v2/gene_scores.tsv'
    scores=list(csv.DictReader(source.open(),delimiter='\t'))
    original={(r['donor'],r['gene'],r['fields']):r for r in scores if r['method']=='ipd_genome'}
    truth=truth_a(REPO/'hla-analysis/results/sequence_catalogue.tsv',{r['donor'] for r in cohort})
    available={};hashes={}
    for panel in ('hprc','hprc_asian'):
        for path in (HERE/'references/observed-v2'/panel).glob('HLA-*.json'):
            metadata=json.loads(path.read_text())
            if not isinstance(metadata,list):continue
            gene=path.stem.removeprefix('HLA-')
            available[panel,gene]={four_field(label) for row in metadata for label in row['genomic_labels'] if four_field(label)}
            hashes[str(path.relative_to(HERE))]=sha(path)
    errors=[];changes=[]
    all_scores=list(csv.DictReader((HERE/'development/all-read-control/gene_scores.tsv').open(),delimiter='\t'))
    assert all(r['state']=='complete' for r in all_scores)
    all_index={(r['donor'],r['gene'],r['fields']):r for r in all_scores if r['method']=='all'}
    for row in scores:
        if row['method']!='ipd_genome':continue
        donor,gene,fields=row['donor'],row['gene'],row['fields']
        other=all_index[donor,gene,fields]
        if row['correct']!=other['correct']:
            changes.append(dict(donor=donor,gene=gene,fields=fields,genomic_ipd_correct=row['correct'],all_read_correct=other['correct']))
        if fields!='4' or row['eligible']!='1' or row['correct']=='1':continue
        slots=genomic_truth(truth[donor,gene]);assert slots is not None
        calls=json.loads((HERE/'work/linear-snapshot/ipd_genome'/donor/'calls.json').read_text())[gene]
        predicted=calls['resolutions']['4']
        item=dict(donor=donor,gene=gene,two_field_correct=original[donor,gene,'2']['correct'],
                  called=row['called'],allele_matches=row['allele_matches'],copy_count=calls['copy_count'],
                  truth1=','.join(sorted(slots[0])),truth2=','.join(sorted(slots[1])),
                  predicted1=predicted[0]['raw'],predicted2=predicted[1]['raw'],
                  disjoint_truth_four_field_labels=int(slots[0].isdisjoint(slots[1])),
                  all_read_correct=other['correct'])
        for panel in ('hprc','hprc_asian'):
            item[panel+'_truth_haplotypes_represented']=sum(bool(slot & available[panel,gene]) for slot in slots)
        errors.append(item)
    output=HERE/'development/residual-error-audit';output.mkdir(exist_ok=True)
    write_tsv(output/'errors.tsv',errors)
    if changes:write_tsv(output/'all-read-changes.tsv',changes)
    summary=dict(scope='Exposed development errors only; graph label availability is not proof of read-level identifiability or rescue',
                 errors=len(errors),by_gene=dict(collections.Counter(r['gene'] for r in errors)),
                 two_field_correct=sum(r['two_field_correct']=='1' for r in errors),
                 disjoint_truth_but_homozygous_call=sum(r['disjoint_truth_four_field_labels'] and r['copy_count']==1 for r in errors),
                 graph_representation={p:dict(collections.Counter(str(r[p+'_truth_haplotypes_represented']) for r in errors)) for p in ('hprc','hprc_asian')},
                 source_sha256=sha(source),metadata_sha256=hashes,driver_sha256=sha(Path(__file__)))
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='metadata_sha256'},indent=2))


if __name__=='__main__':run()
