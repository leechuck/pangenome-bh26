"""Queue the unchanged refinement candidate after matched graph evidence joins."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit


if __name__=='__main__':
    ledger=HERE/'graph-pair-batch-launch.json'
    if ledger.exists():raise FileExistsError('Inspect existing jobs before resubmission')
    pilot=json.loads((HERE/'graph-pair-refinement-launch.json').read_text())
    files=(*pilot['code_sha256'],'run_graph_pair_batch.py','graph-development-batch.json')
    hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    if any(hashes[f]!=digest for f,digest in pilot['code_sha256'].items()):
        raise ValueError('Inference code differs from completed pilot')
    for condition in pilot['jobs']:
        result=json.loads(remote('cat '+B+'/hla/t1k-pangenome/development/graph-pair-refinement-v1/HG00658/'+condition+'/COMPLETE.json'))
        if result['status']!='complete':raise ValueError('Pilot not complete')
    directory=B+'/hla/t1k-pangenome/graph-pair-batch-v1-code'
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(HERE/f) for f in files],HOST+':'+directory+'/'],check=True)
    actual=remote('sha256sum '+' '.join(directory+'/'+f for f in files))
    if {line.split()[1].rsplit('/',1)[-1]:line.split()[0] for line in actual.splitlines()}!=hashes:
        raise ValueError('Code transfer mismatch')
    plan=json.loads((HERE/'graph-development-batch.json').read_text())
    conditions=('hprc/sampled','hprc/unsampled','hprc_asian/sampled','hprc_asian/unsampled')
    paths=[B+'/hla/t1k-pangenome/development/join-v3/'+d['donor']+'/'+c+'/COMPLETE.json' for d in plan['donors'] for c in conditions]
    script='import json,os; paths='+repr(paths)+'; print(all(os.path.isfile(p) and json.load(open(p))["status"]=="complete" for p in paths))'
    complete=remote('python3 -c '+shlex.quote(script))=='True'
    extra=['--array=0-27']
    dependency=json.loads((HERE/'graph-development-batch-launch.json').read_text())['jobs']['join']
    if not complete:
        if not remote('squeue -h -j '+dependency+' -o %T'):raise ValueError('Incomplete join without a live job')
        extra+=['--dependency=aftercorr:'+dependency]
    code='/home/leechuck/hla/t1k-pangenome/graph-pair-batch-v1-code'
    job=submit('t1k-pair-refine-batch',1,'12G','02:00:00',
               ['python3',code+'/run_graph_pair_batch.py','--config',code+'/graph-development-batch.json',
                '--index','$SLURM_ARRAY_TASK_ID'],extra)
    ledger.write_text(json.dumps(dict(job=job,code_sha256=hashes,scope='Fixed seven-donor development expansion; same refinement as HG00658 pilot'),indent=2)+'\n')
    print(job)
