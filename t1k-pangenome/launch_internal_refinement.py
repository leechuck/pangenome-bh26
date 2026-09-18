"""Compare gene-internal paired evidence on the unchanged eight-donor panel."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit


if __name__=='__main__':
    ledger=HERE/'graph-internal-refinement-launch.json'
    if ledger.exists():raise FileExistsError('Inspect existing jobs before resubmission')
    plan=json.loads((HERE/'graph-development-batch.json').read_text())
    donors=[plan['existing_pilot'],*plan['donors']]
    conditions=('hprc/sampled','hprc/unsampled','hprc_asian/sampled','hprc_asian/unsampled')
    paths=[B+'/hla/t1k-pangenome/development/join-v3/'+d['donor']+'/'+c+'/COMPLETE.json'
           for d in donors for c in conditions]
    script='import json; paths='+repr(paths)+'; print(all(json.load(open(p))["status"]=="complete" for p in paths))'
    if remote('python3 -c '+shlex.quote(script))!='True':raise ValueError('Incomplete evidence joins')
    files=(*json.loads((HERE/'graph-pair-refinement-launch.json').read_text())['code_sha256'],
           'run_graph_pair_batch.py','graph-development-batch.json')
    directory=B+'/hla/t1k-pangenome/graph-internal-refinement-v1-code'
    remote('mkdir '+shlex.quote(directory))
    hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    subprocess.run(['scp',*[str(HERE/f) for f in files],HOST+':'+directory+'/'],check=True)
    actual=remote('sha256sum '+' '.join(directory+'/'+f for f in files))
    if {r.split()[1].rsplit('/',1)[-1]:r.split()[0] for r in actual.splitlines()}!=hashes:
        raise ValueError('Transferred code mismatch')
    code='/home/leechuck/hla/t1k-pangenome/graph-internal-refinement-v1-code'
    remote('mkdir -p '+B+'/hla/t1k-pangenome/development/graph-pair-refinement-v2-internal')
    job=submit('t1k-internal-refinement',1,'12G','00:30:00',
               ['python3',code+'/run_graph_pair_batch.py','--config',code+'/graph-development-batch.json',
                '--index','$SLURM_ARRAY_TASK_ID','--internal-pairs','--include-pilot'],
               ['--array=0-31','--partition=batch,debug'])
    ledger.write_text(json.dumps(dict(job=job,code_sha256=hashes,scope='Fixed eight development donors; require both mates inside exon envelope at every placement of the selected locus; no validation outcomes'),indent=2)+'\n')
    print(job)
