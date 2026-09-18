"""Queue fragment evidence and all-locus native/graph joins after mapping."""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'hla-spechla-pg/series20260918'))
from launch_ibex import B, HOST, remote, submit


def main():
    ledger = HERE/'development-evidence-launch.json'
    if ledger.exists():
        raise FileExistsError('Inspect existing jobs before resubmitting')
    mappings = json.loads((HERE/'development-map-all-loci-launch.json').read_text())['jobs']
    native = json.loads((HERE/'native-evidence-development-launch.json').read_text())['job']
    files = ('run_development_evidence.py','run_development_mapping.py','map_personalized.py',
             'build_graph.py','evidence_io.py','fragment_support.py','path_support.py',
             'join_evidence.py','t1k_evidence.py')
    directory = B+'/hla/t1k-pangenome/development-evidence-v1-code'
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(HERE/name) for name in files],HOST+':'+directory+'/'],check=True)
    hashes = {name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files}
    # Verify all small code files in one remote metadata operation.
    actual = remote('sha256sum '+' '.join(shlex.quote(directory+'/'+name) for name in files))
    observed = {Path(line.split()[1]).name:line.split()[0] for line in actual.splitlines()}
    if observed != hashes:
        raise ValueError('Transferred code hash mismatch')
    root = '/home/leechuck/hla/t1k-pangenome'
    code = root+'/development-evidence-v1-code'
    record = dict(donor='HG00658', scope='Development fragment evidence only, without genotype calls',
                  code_sha256=hashes, jobs={})
    def save():
        ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    smoke = submit('t1k-fragment-disk-smoke',1,'4G','00:20:00',
                   ['python3',code+'/run_development_evidence.py','--smoke'])
    record['smoke'] = smoke;save()
    genes = ('A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1')
    for key, mapping in mappings.items():
        panel, selection = key.split('/')
        dependency = 'aftercorr:'+mapping['job']+',afterok:'+smoke
        job = submit('t1k-fragments-'+panel+'-'+selection,1,'8G','02:00:00',
                     ['python3',code+'/run_development_evidence.py','--panel',panel,
                      '--selection',selection,'--index','$SLURM_ARRAY_TASK_ID'],
                     ['--array=0-7','--dependency='+dependency])
        record['jobs'][key] = dict(fragments=job,dependency=dependency);save()
        base = root+'/development/fragments-v1/HG00658/'+key
        command = ['python3',code+'/join_evidence.py','--native',root+'/development/native-evidence-v1/HG00658',
                   '--native-assignments',root+'/development/native-evidence-v1/HG00658/HG00658_assign.tsv',
                   '--output',root+'/development/join-v1/HG00658/'+key,
                   '--graph-fragments',*[base+'/'+gene+'/fragments' for gene in genes],
                   '--graph-support',*[base+'/'+gene+'/support' for gene in genes],
                   '--graph-mapping',*[root+'/development/mapping-v2/HG00658/'+key+'/'+gene for gene in genes]]
        joined = submit('t1k-join-'+panel+'-'+selection,1,'8G','02:00:00',command,
                        ['--dependency=afterok:'+job+':'+native])
        record['jobs'][key]['join'] = joined;save()
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()
