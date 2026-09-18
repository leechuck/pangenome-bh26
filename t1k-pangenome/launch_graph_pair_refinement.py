"""Gate four development refinement conditions on one real-data pilot."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit


if __name__=='__main__':
    ledger=HERE/'graph-pair-refinement-launch.json'
    if ledger.exists():raise FileExistsError('Inspect recorded jobs before another launch')
    files=('run_graph_pair_refinement.py','graph_pair_refinement.py','evidence_io.py','build_graph.py',
           'build_panel.py','t1k_reference.py','run_development_mapping.py','map_personalized.py')
    directory=B+'/hla/t1k-pangenome/graph-pair-refinement-v1-code'
    remote('mkdir -p '+B+'/hla/t1k-pangenome/development/graph-pair-refinement-v1')
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(HERE/f) for f in files],HOST+':'+directory+'/'],check=True)
    hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    actual=remote('sha256sum '+' '.join(directory+'/'+f for f in files))
    if {line.split()[1].rsplit('/',1)[-1]:line.split()[0] for line in actual.splitlines()}!=hashes:
        raise ValueError('Transferred code mismatch')
    root='/home/leechuck/hla/t1k-pangenome'
    base=['python3',root+'/graph-pair-refinement-v1-code/run_graph_pair_refinement.py',
          '--baseline',root+'/development/linear-v1/ipd_genome/HG00658',
          '--calibration',root+'/development/library-calibration-v1/HG00658/COMPLETE.json']
    record=dict(code_sha256=hashes,jobs={},scope='HG00658 development graph four-field refinement; no validation outcomes')
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    def command(panel,selection):
        return base+['--joined',root+'/development/join-v3/HG00658/'+panel+'/'+selection,
                     '--reference',root+'/references/observed-v2/'+panel,
                     '--output',root+'/development/graph-pair-refinement-v1/HG00658/'+panel+'/'+selection]
    save()
    pilot=submit('t1k-pair-refine-pilot',1,'12G','02:00:00',command('hprc_asian','unsampled'),['--partition=batch,debug'])
    record['jobs']['hprc_asian/unsampled']=pilot;save()
    for panel,selection in (('hprc','sampled'),('hprc','unsampled'),('hprc_asian','sampled')):
        record['jobs'][panel+'/'+selection]=submit('t1k-pair-refine-'+panel+'-'+selection,1,'12G','02:00:00',
                                                 command(panel,selection),['--dependency=afterok:'+pilot,'--partition=batch,debug']);save()
    print(json.dumps(record,indent=2))
