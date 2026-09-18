"""Run native parity, all-read, and four graph-recruitment development controls."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit


def main():
    ledger=HERE/'graph-recruitment-launch.json'
    if ledger.exists():raise FileExistsError('Inspect existing jobs before resubmission')
    root='/home/leechuck/hla/t1k-pangenome'
    directory=B+'/hla/t1k-pangenome/graph-recruitment-v1-code'
    files=('graph_recruitment.py','evidence_io.py','join_evidence.py','t1k_evidence.py',
           'build_graph.py','decode_t1k.py','t1k_reference.py','build_panel.py')
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(HERE/f) for f in files],HOST+':'+directory+'/'],check=True)
    hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    actual=remote('sha256sum '+' '.join(directory+'/'+f for f in files))
    if {line.split()[1].rsplit('/',1)[-1]:line.split()[0] for line in actual.splitlines()}!=hashes:
        raise ValueError('Code transfer mismatch')
    record=dict(code_sha256=hashes,jobs={},scope='HG00658 development: graph read selection with native T1K inference')
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    base=['python3',root+'/graph-recruitment-v1-code/graph_recruitment.py',
          '--baseline',root+'/development/linear-v1/ipd_genome/HG00658',
          '--reads','/home/leechuck/hla/hla-typer/reads/HG00658',
          '--reference',root+'/t1k-references/genome-v2','--tool',root+'/tools/t1k-evidence-v1','--threads','4']
    output=root+'/development/graph-recruitment-v1/HG00658/'
    native=submit('t1k-recruit-native',4,'12G','02:00:00',base+['--mode','native','--output',output+'native'],['--partition=debug'])
    record['jobs']['native']=native;save()
    record['jobs']['all']=submit('t1k-recruit-all',4,'12G','02:00:00',base+['--mode','all','--output',output+'all'],
                                ['--partition=debug','--dependency=afterok:'+native]);save()
    joins=json.loads((HERE/'development-evidence-v3-launch.json').read_text())['jobs']
    for key,jobs in joins.items():
        graph=root+'/development/join-v3/HG00658/'+key
        hostgraph=graph.replace('/home/leechuck/hla',B+'/hla',1)
        complete=remote('test -f '+shlex.quote(hostgraph+'/COMPLETE.json')+' && echo yes || echo no')=='yes'
        dependency=native
        if not complete:
            job=jobs['join']
            if not remote('squeue -h -j '+job+' -o %T'):raise ValueError('Missing join without live job: '+key)
            dependency+=':'+job
        record['jobs'][key]=submit('t1k-recruit-'+key.replace('/','-'),4,'12G','02:00:00',
                                  base+['--mode','graph','--graph',graph,'--output',output+key],
                                  ['--partition=debug','--dependency=afterok:'+dependency]);save()
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
