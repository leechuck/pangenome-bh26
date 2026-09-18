"""Execute frozen T1K on reserved reads without loading truth or scoring calls."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time
from evidence_io import sha


def run(index,assets,output_root,threads=4):
    if not 1<=threads<=int(os.environ.get('SLURM_CPUS_PER_TASK','0')):
        raise ValueError('Requires Slurm allocation')
    plan=json.loads((assets/'validation-baseline-plan.json').read_text())
    for filename,key in (('cohort.tsv','cohort_sha256'),('READ_PREPARATION.json','preparation_sha256'),('RESERVATION.json','reservation_sha256')):
        if sha(assets/filename)!=plan[key]:raise ValueError('Reserved input metadata changed: '+filename)
    cohort=list(csv.DictReader((assets/'cohort.tsv').open(),delimiter='\t'))
    if not 0<=index<len(cohort):raise ValueError('Invalid reserved donor index')
    donor=cohort[index]['donor']
    if not donor or Path(donor).name!=donor or donor in ('.','..'):raise ValueError('Unsafe donor')
    prepared=json.loads((assets/'READ_PREPARATION.json').read_text())
    if prepared['status']!='read_preparation_complete' or prepared['cohort_sha256']!=plan['cohort_sha256']:
        raise ValueError('Read preparation does not match reserved cohort')
    read_record=[r for r in prepared['donors'] if r['donor']==donor]
    if len(read_record)!=1:raise ValueError('Missing or repeated preparation donor')
    root=Path('/home/leechuck/hla')
    reads=root/'t1k-pangenome/hgsvc-validation-v1/validation-reads'/donor
    if not (reads/'COMPLETE').exists():raise ValueError('Read preparation not complete')
    actual={name:sha(reads/name) for name in ('r1.fq.gz','r2.fq.gz')}
    if actual!=read_record[0]['reads_sha256']:raise ValueError('Reserved reads changed')
    reference=root/'hla-typer/t1kdb/hla365_dna_seq.fa'
    exe=root/'mm/envs/typing/bin/run-t1k'
    recipe=plan['recipe']
    if sha(reference)!=recipe['database_sha256'] or sha(exe)!=recipe['tool_sha256']:
        raise ValueError('Frozen T1K reference or wrapper changed')
    if recipe['alleleDigitUnits']!=4 or recipe['alleleDelimiter']!=':' or recipe['preset']!='hla-wgs':
        raise ValueError('Unexpected baseline recipe')
    output=output_root/donor
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    command=[str(exe),'-1',str(reads/'r1.fq.gz'),'-2',str(reads/'r2.fq.gz'),
             '--preset','hla-wgs','-f',str(reference),'-t',str(threads),
             '--alleleDigitUnits','4','--alleleDelimiter',':','-o',str(output/donor)]
    record=dict(status='running',donor=donor,started=time.time(),plan_sha256=sha(assets/'validation-baseline-plan.json'),
                driver_sha256=sha(Path(__file__)),reads_sha256=actual,command=command,
                recipe=recipe,scope='Reserved baseline calls only; no truth loaded or outcomes scored')
    def save():(output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        env=dict(os.environ,PATH=str(exe.parent)+':'+os.environ['PATH'])
        with (output/'run.log').open('w') as log:
            subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        table=output/(donor+'_genotype.tsv')
        if not table.is_file() or not table.stat().st_size:raise ValueError('Missing genotype output')
        record.update(status='complete',finished=time.time(),output_sha256={table.name:sha(table)})
        save();(output/'COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
    except BaseException as error:
        record.update(status='failed',finished=time.time(),error=repr(error));save();raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--index',type=int,required=True)
    p.add_argument('--assets',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.index,a.assets,a.output)
