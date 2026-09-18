"""Run the genomic-IPD control on reserved reads without evaluation truth."""
import argparse
import csv
import json
import os
from pathlib import Path
from evidence_io import sha
from run_linear_control import run as linear


def run(index,assets,output):
    if int(os.environ.get('SLURM_CPUS_PER_TASK','0'))<4:raise ValueError('Requires four Slurm CPUs')
    plan=json.loads((assets/'validation-genomic-plan.json').read_text())
    for filename,key in (('cohort.tsv','cohort_sha256'),('READ_PREPARATION.json','preparation_sha256'),('RESERVATION.json','reservation_sha256')):
        if sha(assets/filename)!=plan[key]:raise ValueError('Reserved metadata changed')
    cohort=list(csv.DictReader((assets/'cohort.tsv').open(),delimiter='\t'))
    if not 0<=index<len(cohort):raise ValueError('Invalid donor index')
    donor=cohort[index]['donor']
    if not donor or Path(donor).name!=donor or donor in ('.','..'):raise ValueError('Unsafe donor')
    prep=json.loads((assets/'READ_PREPARATION.json').read_text())
    if prep['status']!='read_preparation_complete' or prep['cohort_sha256']!=plan['cohort_sha256']:
        raise ValueError('Preparation mismatch')
    records=[r for r in prep['donors'] if r['donor']==donor]
    if len(records)!=1:raise ValueError('Missing or duplicate read record')
    root=Path('/home/leechuck/hla/t1k-pangenome')
    reads=root/'hgsvc-validation-v1/validation-reads'/donor
    if {name:sha(reads/name) for name in ('r1.fq.gz','r2.fq.gz')}!=records[0]['reads_sha256']:
        raise ValueError('Reserved reads changed')
    reference=root/'t1k-references/genome-v2'
    for name,digest in plan['reference_sha256'].items():
        if sha(reference/name)!=digest:raise ValueError('Genomic reference changed')
    if sha(Path('/home/leechuck/hla/mm/envs/typing/bin/run-t1k'))!=plan['tool_sha256']:
        raise ValueError('T1K wrapper changed')
    linear(assets/'cohort.tsv',index,reference,root/'hgsvc-validation-v1/validation-reads',output,4,
           scope='Reserved genomic-IPD control; no truth loaded or outcomes scored')
    folder=output/donor
    report=json.loads((folder/'COMPLETE.json').read_text())
    report.update(validation_plan_sha256=sha(assets/'validation-genomic-plan.json'),
                  validation_driver_sha256=sha(Path(__file__)))
    for name in ('manifest.json','COMPLETE.json'):
        (folder/name).write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index',type=int,required=True);p.add_argument('--assets',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.index,a.assets,a.output)
