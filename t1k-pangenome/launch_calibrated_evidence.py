"""Project calibrated mappings and join native evidence for all graph conditions."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit


def main():
    ledger=HERE/'development-evidence-v3-launch.json'
    if ledger.exists():raise FileExistsError('Inspect existing jobs before resubmitting')
    root='/home/leechuck/hla/t1k-pangenome'
    files=('run_development_evidence.py','run_development_mapping.py','map_personalized.py',
           'build_graph.py','evidence_io.py','fragment_support.py','path_support.py',
           'join_evidence.py','t1k_evidence.py','projection_graph.py')
    directory=B+'/hla/t1k-pangenome/development-evidence-v3-code'
    exists=remote('test -d '+shlex.quote(directory)+' && echo yes || echo no')
    if exists=='no':
        remote('mkdir '+shlex.quote(directory))
        subprocess.run(['scp',*[str(HERE/f) for f in files],HOST+':'+directory+'/'],check=True)
    hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    actual=remote('sha256sum '+' '.join(directory+'/'+f for f in files))
    if {line.split()[1].rsplit('/',1)[-1]:line.split()[0] for line in actual.splitlines()}!=hashes:
        raise ValueError('Code transfer mismatch')
    # Inputs are complete already; do not depend on aged-out Slurm IDs.
    checks=[root+'/development/mapping-v3/HG00658/'+panel+'/'+selection+'/'+gene+'/COMPLETE.json'
            for panel in ('hprc','hprc_asian') for selection in ('sampled','unsampled')
            for gene in ('A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1')]
    native=root+'/development/native-evidence-v1/HG00658'
    paths=[path.replace('/home/leechuck/hla',B+'/hla',1) for path in [native+'/COMPLETE.json',*checks]]
    script='import json,os; paths='+repr(paths)+'; print(json.dumps([json.load(open(p)) if os.path.exists(p) else None for p in paths]))'
    inputs=json.loads(remote('python3 -c '+shlex.quote(script)))
    if len(inputs)!=len(paths) or inputs[0]['status']!='complete':
        raise ValueError('Incomplete native input')
    mapping_jobs=json.loads((HERE/'development-map-calibrated-launch.json').read_text())['jobs']
    dependencies={}
    for i,key in enumerate(('hprc/sampled','hprc/unsampled','hprc_asian/sampled','hprc_asian/unsampled')):
        group=inputs[1+8*i:1+8*(i+1)]
        dependencies[key]=[]
        if any(item is None or item['status']!='complete' for item in group):
            job=mapping_jobs[key]
            if not remote('squeue -h -j '+job+' -o %T'):
                raise ValueError('Incomplete mapping without live job: '+key)
            dependencies[key]=['--dependency=afterok:'+job]
    record=dict(code_sha256=hashes,jobs={},scope='Calibrated graph evidence; no genotype truth')
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    save();code=root+'/development-evidence-v3-code'
    for panel in ('hprc','hprc_asian'):
        for selection in ('sampled','unsampled'):
            key=panel+'/'+selection
            job=submit('t1k-evidence-v3-'+panel+'-'+selection,1,'8G','02:00:00',
                       ['python3',code+'/run_development_evidence.py','--version','v3',
                        '--panel',panel,'--selection',selection,'--index','$SLURM_ARRAY_TASK_ID'],
                       ['--array=0-7','--partition=debug',*dependencies[key]])
            record['jobs'][key]=dict(fragments=job);save()
            prefix=root+'/development/fragments-v3/HG00658/'+key
            genes=('A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1')
            command=['python3',code+'/join_evidence.py','--native',native,
                     '--native-assignments',native+'/HG00658_assign.tsv',
                     '--output',root+'/development/join-v3/HG00658/'+key,
                     '--graph-fragments',*[prefix+'/'+gene+'/fragments' for gene in genes],
                     '--graph-support',*[prefix+'/'+gene+'/support' for gene in genes],
                     '--graph-mapping',*[root+'/development/mapping-v3/HG00658/'+key+'/'+gene for gene in genes]]
            record['jobs'][key]['join']=submit('t1k-join-v3-'+panel+'-'+selection,1,'8G','02:00:00',command,
                                             ['--partition=debug','--dependency=afterok:'+job]);save()
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
