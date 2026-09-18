"""Queue matched-64 development controls after the existing benchmark jobs."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
SERIES = HERE.parent/'hla-spechla-pg/series20260918'
sys.path.insert(0,str(SERIES))
from launch_ibex import B, HOST, R, submit


def main():
    ledger = HERE/'linear-control-launch.json'
    if ledger.exists():
        raise FileExistsError('Inspect existing job IDs before resubmitting')
    previous = json.loads((SERIES/'ibex-launch.json').read_text())['jobs']
    dependencies = [previous[key] for key in (
        'gourraud_packed_full','gourraud_packed_native',
        'gourraud_packed_asian_matched','gourraud_asian_NA20790_arpack_retry')]
    files = ['run_linear_control.py','decode_t1k.py','t1k_reference.py','build_panel.py','build_graph.py']
    target = '/home/leechuck/hla/t1k-pangenome'
    subprocess.run(['scp',*[str(HERE/name) for name in files],
                    HOST+':'+B+'/hla/t1k-pangenome/'],check=True)
    record = dict(scope='Matched-64 development only; native T1K with observed sequence additions',
                  dependencies=dependencies, dependency_policy='afterany: retain failed benchmark samples',
                  code_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
                  jobs={})
    def save():
        ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    for panel in ('hprc','hprc_asian'):
        base = ['python3',target+'/run_linear_control.py','--cohort',R+'/cohort.tsv',
                '--reference',target+'/t1k-references/v2/'+panel,
                '--reads-root','/home/leechuck/hla/hla-typer/reads',
                '--output-root',target+'/development/linear-v1/'+panel,'--threads','4']
        pilot = submit('t1k-linear-'+panel+'-pilot',4,'8G','01:00:00',
                       base+['--index','0'],['--dependency=afterany:'+':'.join(dependencies)])
        record['jobs'][panel] = dict(pilot=pilot)
        save()
        array = submit('t1k-linear-'+panel,4,'8G','01:00:00',
                       base+['--index','$SLURM_ARRAY_TASK_ID'],
                       ['--array=1-63','--dependency=afterok:'+pilot])
        record['jobs'][panel]['array'] = array
        save()
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()
