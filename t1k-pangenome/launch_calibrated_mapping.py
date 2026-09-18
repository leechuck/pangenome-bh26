"""Remap development reads using one shared, truth-free insert calibration."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE, B, HOST, remote, submit


def main():
    ledger = HERE/'development-map-calibrated-launch.json'
    if ledger.exists():
        raise FileExistsError('Inspect recorded jobs before another submission')
    root = '/home/leechuck/hla/t1k-pangenome'
    calibration = root+'/development/library-calibration-v1/HG00658/COMPLETE.json'
    observed = json.loads(remote('cat '+B+'/hla/t1k-pangenome/development/library-calibration-v1/HG00658/COMPLETE.json'))
    if observed != json.loads((HERE/'library-calibration-result.json').read_text()):
        raise ValueError('Calibration changed')
    if observed['status'] != 'complete' or not observed['estimate']['usable']:
        raise ValueError('Calibration unusable')
    files = ('run_development_mapping.py','map_personalized.py','build_graph.py')
    directory = B+'/hla/t1k-pangenome/development-map-v3-code'
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(HERE/f) for f in files],HOST+':'+directory+'/'],check=True)
    hashes = {f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in files}
    actual = remote('sha256sum '+' '.join(directory+'/'+f for f in files))
    if {line.split()[1].rsplit('/',1)[-1]:line.split()[0] for line in actual.splitlines()} != hashes:
        raise ValueError('Code transfer mismatch')
    command = ['python3',root+'/development-map-v3-code/run_development_mapping.py',
               '--calibration',calibration]
    record = dict(code_sha256=hashes,jobs={},scope='Matched graph mapping with shared insert parameters; no genotype truth')
    def save():
        ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    pilot = submit('t1k-calibrated-map-pilot',4,'12G','02:00:00',
                   command+['--panel','hprc','--selection','unsampled','--index','0'],
                   ['--partition=debug'])
    record['pilot']=pilot;save()
    for panel in ('hprc','hprc_asian'):
        for selection in ('sampled','unsampled'):
            array = '1-7' if (panel,selection)==('hprc','unsampled') else '0-7'
            job = submit('t1k-calibrated-'+panel+'-'+selection,4,'12G','02:00:00',
                         command+['--panel',panel,'--selection',selection,'--index','$SLURM_ARRAY_TASK_ID'],
                         ['--array='+array,'--partition=debug','--dependency=afterok:'+pilot])
            record['jobs'][panel+'/'+selection]=job;save()
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()
