"""Expand the all-read control across the fixed development cohort."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit
from run_all_read_control import COHORT_SHA


if __name__=='__main__':
    ledger=HERE/'all-read-control-launch.json'
    if ledger.exists():raise FileExistsError('Inspect existing jobs first')
    if remote('sha256sum '+B+'/hla/dogohla-series20260918/cohort.tsv').split()[0]!=COHORT_SHA:
        raise ValueError('Development cohort path or contents changed')
    previous=json.loads((HERE/'graph-recruitment-launch.json').read_text())
    directory=B+'/hla/t1k-pangenome/all-read-control-v1-code'
    files=(*previous['code_sha256'],'run_all_read_control.py')
    hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    if any(hashes[f]!=digest for f,digest in previous['code_sha256'].items()):
        raise ValueError('Pilot inference code changed')
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(HERE/f) for f in files],HOST+':'+directory+'/'],check=True)
    actual=remote('sha256sum '+' '.join(directory+'/'+f for f in files))
    if {line.split()[1].rsplit('/',1)[-1]:line.split()[0] for line in actual.splitlines()}!=hashes:
        raise ValueError('Code transfer mismatch')
    root='/home/leechuck/hla/t1k-pangenome'
    pilot=previous['jobs']['native']
    completed=remote('test -f '+B+'/hla/t1k-pangenome/development/graph-recruitment-v1/HG00658/native/COMPLETE.json && echo yes || echo no')=='yes'
    extra=['--array=1-63']
    if not completed:
        if not remote('squeue -h -j '+pilot+' -o %T'):raise ValueError('Native parity not complete or live')
        extra+=['--dependency=afterok:'+pilot]
    job=submit('t1k-all-read-development',4,'12G','02:00:00',
               ['python3',root+'/all-read-control-v1-code/run_all_read_control.py','--index','$SLURM_ARRAY_TASK_ID',
                '--cohort','/home/leechuck/hla/dogohla-series20260918/cohort.tsv'],extra)
    ledger.write_text(json.dumps(dict(job=job,pilot_job=previous['jobs']['all'],code_sha256=hashes,
                      scope='All 64 fixed development donors; index zero uses existing pilot; no validation donors'),indent=2)+'\n')
    print(ledger.read_text())
