"""Submit balanced development graph comparisons with explicit stage dependencies."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit


def main():
    ledger=HERE/'graph-development-batch-launch.json'
    if ledger.exists():raise FileExistsError('Inspect existing stage jobs before resubmission')
    config=json.loads((HERE/'graph-development-batch.json').read_text())
    if remote('sha256sum '+B+'/hla/dogohla-series20260918/cohort.tsv').split()[0]!=config['cohort_sha256']:
        raise ValueError('Cohort missing or changed')
    if len({r['family'] for r in [config['existing_pilot'],*config['donors']]})!=8:
        raise ValueError('Development families overlap')
    files=('run_graph_development_batch.py','graph-development-batch.json','run_development_mapping.py',
           'map_personalized.py','run_development_evidence.py','kmer_coverage.py',
           'library_calibration.py','run_native_evidence.py','join_evidence.py',
           'graph_recruitment.py','evidence_io.py','fragment_support.py','path_support.py',
           'projection_graph.py','build_graph.py','build_panel.py','t1k_evidence.py',
           'decode_t1k.py','t1k_reference.py')
    directory=B+'/hla/t1k-pangenome/graph-development-batch-v1-code'
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(HERE/f) for f in files],HOST+':'+directory+'/'],check=True)
    hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    actual=remote('sha256sum '+' '.join(directory+'/'+f for f in files))
    if {line.split()[1].rsplit('/',1)[-1]:line.split()[0] for line in actual.splitlines()}!=hashes:
        raise ValueError('Transferred source differs')
    root='/home/leechuck/hla/t1k-pangenome/graph-development-batch-v1-code'
    base=['python3',root+'/run_graph_development_batch.py','--config',root+'/graph-development-batch.json',
          '--index','$SLURM_ARRAY_TASK_ID']
    record=dict(code_sha256=hashes,jobs={},donors=config['donors'],scope=config['scope'])
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    def stage(name,count,cpus,mem,dependencies=()):
        extra=['--array=0-'+str(count-1),'--partition=debug']
        if dependencies:extra+=['--dependency='+','.join(kind+':'+record['jobs'][parent] for kind,parent in dependencies)]
        job=submit('t1k-batch-'+name,cpus,mem,'02:00:00',base+['--stage',name],extra)
        record['jobs'][name]=job;save()
    save()
    n=len(config['donors'])
    stage('coverage',n,1,'8G')
    stage('native',n,4,'12G')
    stage('bootstrap',n*8,4,'12G',[('afterok','coverage')])
    stage('bootstrap_evidence',n*8,1,'8G',[('aftercorr','bootstrap')])
    stage('calibration',n,8,'16G',[('afterok','bootstrap_evidence')])
    stage('mapping',n*32,4,'12G',[('afterok','calibration')])
    stage('evidence',n*32,1,'8G',[('aftercorr','mapping')])
    stage('join',n*4,1,'8G',[('afterok','evidence'),('afterok','native')])
    stage('recruit',n*4,4,'12G',[('aftercorr','join')])
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
