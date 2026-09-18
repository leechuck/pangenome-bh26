"""Queue all-locus graph evidence after the legacy benchmark and mapping pilots."""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'hla-spechla-pg/series20260918'))
from launch_ibex import B, HOST, remote, submit


def main():
    ledger = HERE/'development-map-all-loci-launch.json'
    if ledger.exists():
        raise FileExistsError('Inspect existing jobs before resubmitting')
    pilots = json.loads((HERE/'development-map-launch.json').read_text())['jobs']
    files = ('run_development_mapping.py', 'map_personalized.py', 'build_graph.py')
    directory = B+'/hla/t1k-pangenome/development-map-v2-code'
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp', *[str(HERE/name) for name in files],
                    HOST+':'+directory+'/'], check=True)
    hashes = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files}
    for name, digest in hashes.items():
        actual = remote('sha256sum '+shlex.quote(directory+'/'+name)).split()[0]
        if actual != digest:
            raise ValueError('Transferred code hash mismatch: '+name)
    record = dict(donor='HG00658', genes=['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1'],
                  scope='Development alignment evidence only; no validation truth or genotype claims',
                  coverage='coverage/development-genome-screen-v1/HG00658.json',
                  code_sha256=hashes, jobs={},
                  legacy_dependency='52039862', coverage_dependency='52045735')
    def save():
        ledger.write_text(json.dumps(record, indent=2)+'\n')
    save()
    for panel in ('hprc', 'hprc_asian'):
        for selection in ('sampled', 'unsampled'):
            key = panel+'/'+selection
            dependency = 'afterany:52039862,afterok:52045735:'+pilots[key]
            command = ['python3', '/home/leechuck/hla/t1k-pangenome/development-map-v2-code/run_development_mapping.py',
                       '--panel', panel, '--selection', selection, '--index', '$SLURM_ARRAY_TASK_ID']
            job = submit('t1k-map8-'+panel+'-'+selection, 4, '12G', '02:00:00', command,
                         ['--array=0-7', '--dependency='+dependency])
            record['jobs'][key] = dict(job=job, dependency=dependency)
            save()
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
