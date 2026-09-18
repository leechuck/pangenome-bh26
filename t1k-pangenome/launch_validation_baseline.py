"""Launch outcome-blinded frozen baseline typing for all reserved donors."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit


if __name__=='__main__':
    ledger=HERE/'validation-baseline-launch.json'
    if ledger.exists():raise FileExistsError('Inspect existing jobs before resubmitting')
    files=[HERE/name for name in ('run_validation_baseline.py','evidence_io.py','validation-baseline-plan.json')]
    files += [HERE/'validation'/name for name in ('cohort.tsv','READ_PREPARATION.json','RESERVATION.json')]
    directory=B+'/hla/t1k-pangenome/validation-baseline-v1-code'
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(p) for p in files],HOST+':'+directory+'/'],check=True)
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    actual=remote('sha256sum '+' '.join(directory+'/'+p.name for p in files))
    if {line.split()[1].rsplit('/',1)[-1]:line.split()[0] for line in actual.splitlines()}!=hashes:
        raise ValueError('Transferred baseline assets differ')
    root='/home/leechuck/hla/t1k-pangenome';code=root+'/validation-baseline-v1-code'
    command=['python3',code+'/run_validation_baseline.py','--assets',code,
             '--output',root+'/validation-runs/frozen-t1k-v1']
    record=dict(code_sha256=hashes,scope='84 reserved baseline predictions only; no truth or scoring; candidate not yet frozen')
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    pilot=submit('t1k-validation-baseline-pilot',4,'12G','02:00:00',command+['--index','0'])
    record['pilot_job']=pilot;save()
    record['array_job']=submit('t1k-validation-baseline',4,'12G','02:00:00',
                               command+['--index','$SLURM_ARRAY_TASK_ID'],
                               ['--array=1-83','--dependency=afterok:'+pilot]);save()
    print(json.dumps(record,indent=2))
