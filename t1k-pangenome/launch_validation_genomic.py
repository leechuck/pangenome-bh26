"""Submit outcome-blinded genomic-IPD control with a pilot gate."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit


if __name__=='__main__':
    ledger=HERE/'validation-genomic-launch.json'
    if ledger.exists():raise FileExistsError('Inspect existing jobs before resubmission')
    files=[HERE/f for f in ('run_validation_genomic.py','run_linear_control.py','evidence_io.py',
           'build_graph.py','build_panel.py','decode_t1k.py','t1k_reference.py','validation-genomic-plan.json')]
    files += [HERE/'validation'/f for f in ('cohort.tsv','READ_PREPARATION.json','RESERVATION.json')]
    directory=B+'/hla/t1k-pangenome/validation-genomic-v1-code'
    remote('mkdir '+shlex.quote(directory))
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    subprocess.run(['scp',*[str(p) for p in files],HOST+':'+directory+'/'],check=True)
    actual=remote('sha256sum '+' '.join(directory+'/'+p.name for p in files))
    if {r.split()[1].rsplit('/',1)[-1]:r.split()[0] for r in actual.splitlines()}!=hashes:
        raise ValueError('Transferred assets differ')
    root='/home/leechuck/hla/t1k-pangenome';code=root+'/validation-genomic-v1-code'
    command=['python3',code+'/run_validation_genomic.py','--assets',code,
             '--output',root+'/validation-runs/genomic-ipd-v1']
    record=dict(code_sha256=hashes,scope='84 reserved genomic-IPD predictions; no outcome fetch or scoring')
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    pilot=submit('t1k-validation-genomic-pilot',4,'12G','02:00:00',command+['--index','0'],['--partition=batch,debug'])
    record['pilot_job']=pilot;save()
    record['array_job']=submit('t1k-validation-genomic',4,'12G','02:00:00',
        command+['--index','$SLURM_ARRAY_TASK_ID'],['--array=1-83','--dependency=afterok:'+pilot,'--partition=batch,debug'])
    save();print(json.dumps(record,indent=2))
