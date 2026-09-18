"""Retry development evidence with GBZ-coordinate paths; preserve old attempts."""
import hashlib
import json
import shlex
import subprocess
import sys
from launch_development_evidence import HERE,B,HOST,remote,submit


def main():
    ledger=HERE/'development-evidence-v2-launch.json'
    resume='--resume' in sys.argv
    if ledger.exists() and not resume:raise FileExistsError('Inspect existing jobs before resubmitting')
    if resume and not ledger.exists():raise FileNotFoundError(ledger)
    files=('run_development_evidence.py','run_development_mapping.py','map_personalized.py',
           'build_graph.py','evidence_io.py','fragment_support.py','path_support.py',
           'join_evidence.py','t1k_evidence.py','projection_graph.py')
    directory=B+'/hla/t1k-pangenome/development-evidence-v2-code'
    if not resume:
        remote('mkdir '+shlex.quote(directory))
        subprocess.run(['scp',*[str(HERE/name) for name in files],HOST+':'+directory+'/'],check=True)
    hashes={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files}
    actual=remote('sha256sum '+' '.join(shlex.quote(directory+'/'+name) for name in files))
    if {line.split()[1].rsplit('/',1)[-1]:line.split()[0] for line in actual.splitlines()}!=hashes:
        raise ValueError('Code transfer mismatch')
    root='/home/leechuck/hla/t1k-pangenome';code=root+'/development-evidence-v2-code'
    record=json.loads(ledger.read_text()) if resume else dict(scope='Repair GBZ/GFA node-coordinate mismatch; all observed paths sequence-verified; no genotype calls',code_sha256=hashes,jobs={})
    if record['code_sha256']!=hashes:raise ValueError('Cannot resume with changed inference code')
    native=json.loads(remote('cat '+B+'/hla/t1k-pangenome/development/native-evidence-v1/HG00658/COMPLETE.json'))
    if native!=json.loads((HERE/'native-evidence-development-result.json').read_text()):
        raise ValueError('Completed native evidence changed')
    record['native_dependency']='Verified complete manifest; completed Slurm job has aged out of live dependency records'
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    base=['python3',code+'/run_development_evidence.py']
    pilot=record.get('pilot') or submit('t1k-projection-repair-pilot',1,'8G','02:00:00',
                 base+['--panel','hprc','--selection','unsampled','--index','7'],['--partition=debug'])
    record['pilot']=pilot;save()
    genes=('A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1')
    for panel in ('hprc','hprc_asian'):
        for selection in ('sampled','unsampled'):
            key=panel+'/'+selection
            array='0-6' if key=='hprc/unsampled' else '0-7'
            job=record['jobs'].get(key,{}).get('fragments') or submit('t1k-projection-v2-'+panel+'-'+selection,1,'8G','02:00:00',
                       base+['--panel',panel,'--selection',selection,'--index','$SLURM_ARRAY_TASK_ID'],
                       ['--array='+array,'--dependency=afterok:'+pilot,'--partition=debug'])
            record['jobs'].setdefault(key,{})['fragments']=job;save()
            if 'join' in record['jobs'][key]:continue
            prefix=root+'/development/fragments-v2/HG00658/'+key
            command=['python3',code+'/join_evidence.py','--native',root+'/development/native-evidence-v1/HG00658',
                     '--native-assignments',root+'/development/native-evidence-v1/HG00658/HG00658_assign.tsv',
                     '--output',root+'/development/join-v2/HG00658/'+key,
                     '--graph-fragments',*[prefix+'/'+gene+'/fragments' for gene in genes],
                     '--graph-support',*[prefix+'/'+gene+'/support' for gene in genes],
                     '--graph-mapping',*[root+'/development/mapping-v2/HG00658/'+key+'/'+gene for gene in genes]]
            joined=submit('t1k-join-v2-'+panel+'-'+selection,1,'8G','02:00:00',command,
                          ['--dependency=afterok:'+job+':'+pilot,'--partition=debug'])
            record['jobs'][key]['join']=joined;save()
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
