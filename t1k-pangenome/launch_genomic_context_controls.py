"""Separate observed-context effects from restoration of genomic IPD sequence."""
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'hla-spechla-pg/series20260918'))
from launch_ibex import B, HOST, R, remote, submit


def main():
    ledger = HERE/'genomic-context-control-launch.json'
    if ledger.exists():
        raise FileExistsError('Inspect existing job IDs before resubmitting')
    files = ('t1k_reference.py','build_graph.py','build_panel.py','run_linear_control.py','decode_t1k.py')
    directory = B+'/hla/t1k-pangenome/genomic-context-v1-code'
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(HERE/name) for name in files],HOST+':'+directory+'/'],check=True)
    hashes = {name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files}
    actual = remote('sha256sum '+' '.join(shlex.quote(directory+'/'+name) for name in files))
    if {Path(line.split()[1]).name:line.split()[0] for line in actual.splitlines()} != hashes:
        raise ValueError('Code transfer mismatch')
    root = '/home/leechuck/hla/t1k-pangenome'
    code = root+'/genomic-context-v1-code'
    record = dict(scope='Matched-64 development controls: full-genomic IPD plus observed contexts; no graph alignment',
                  code_sha256=hashes,genomic_reference_sha256='9202c3d37a89045a8f01f08dfd3ba21c465835b72ac2e9e7a731ef4dfc2c1d6a',jobs={})
    def save():
        ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    for panel in ('hprc','hprc_asian'):
        reference = root+'/t1k-references/genome-context-v1/'+panel
        prep = submit('t1k-genomic-context-ref-'+panel,2,'8G','00:30:00',
                      ['python3',code+'/t1k_reference.py','--reference',root+'/references/observed-v2',
                       '--panel',panel,'--ipd',root+'/t1k-references/genome-v2/reference.fa','--output',reference],
                      ['--partition=debug'])
        record['jobs'][panel] = dict(reference=prep);save()
        base = ['python3',code+'/run_linear_control.py','--cohort',R+'/cohort.tsv',
                '--reference',reference,'--reads-root','/home/leechuck/hla/hla-typer/reads',
                '--output-root',root+'/development/linear-v1/genome_'+panel,'--threads','4']
        pilot = submit('t1k-genomic-context-pilot-'+panel,4,'12G','02:00:00',base+['--index','0'],
                       ['--dependency=afterok:'+prep,'--partition=debug'])
        record['jobs'][panel]['pilot'] = pilot;save()
        array = submit('t1k-genomic-context-'+panel,4,'12G','02:00:00',
                       base+['--index','$SLURM_ARRAY_TASK_ID'],['--array=1-63','--dependency=afterok:'+pilot])
        record['jobs'][panel]['array'] = array;save()
    print(json.dumps(record,indent=2))


if __name__ == '__main__':
    main()
