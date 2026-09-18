"""Run a development donor through native T1K with observed context additions."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import time
from build_graph import sha
from decode_t1k import decode


def run(cohort, index, reference, reads_root, output_root, threads,
        scope='Development linear-reference control; no graph alignment'):
    if not os.environ.get('SLURM_CPUS_PER_TASK'):
        raise RuntimeError('Run under Slurm')
    if threads < 1 or threads > int(os.environ['SLURM_CPUS_PER_TASK']):
        raise ValueError('Invalid thread allocation')
    rows = list(csv.DictReader(cohort.open(), delimiter='\t'))
    if index < 0 or index >= len(rows):
        raise ValueError('Cohort index out of range')
    donor = rows[index]['donor']
    if not donor or Path(donor).name != donor or donor in ('.','..'):
        raise ValueError('Unsafe donor identifier')
    output = output_root/donor
    if output.exists():
        raise FileExistsError(output)
    manifest = json.loads((reference/'COMPLETE.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('Reference preparation incomplete')
    for name in ('reference.fa','aliases.json'):
        if sha(reference/name) != manifest['output_sha256'][name]:
            raise ValueError('Reference changed: '+name)
    reads = reads_root/donor
    if not (reads/'COMPLETE').is_file():
        raise ValueError('Development reads not complete')
    exe = Path('/home/leechuck/hla/mm/envs/typing/bin/run-t1k')
    record = dict(status='running', donor=donor, started=time.time(),
                  scope=scope,
                  reference_manifest_sha256=sha(reference/'COMPLETE.json'),
                  cohort_sha256=sha(cohort), tool_sha256=sha(exe),
                  driver_sha256=sha(Path(__file__)),
                  decoder_sha256=sha(Path(__file__).with_name('decode_t1k.py')),
                  reads_sha256={str(i):sha(reads/f'r{i}.fq.gz') for i in (1,2)})
    command = [str(exe),'-1',str(reads/'r1.fq.gz'),'-2',str(reads/'r2.fq.gz'),
               '--preset','hla-wgs','-f',str(reference/'reference.fa'),'-t',str(threads),
               '--alleleDigitUnits','4','--alleleDelimiter',':','-o',str(output/donor)]
    record['command'] = command
    output.mkdir(parents=True)
    def save():
        (output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        env = dict(os.environ, PATH=str(exe.parent)+':'+os.environ['PATH'])
        with (output/'run.log').open('w') as log:
            subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        table = output/(donor+'_genotype.tsv')
        calls = decode(table.read_text(),json.loads((reference/'aliases.json').read_text()))
        if not calls:
            raise ValueError('Empty genotype table')
        (output/'calls.json').write_text(json.dumps(calls,indent=2)+'\n')
        record.update(status='complete',finished=time.time(),
                      output_sha256={name:sha(output/name) for name in (table.name,'calls.json')})
        save()
        (output/'COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
    except BaseException as error:
        record.update(status='failed',finished=time.time(),error=repr(error))
        save()
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort',type=Path,required=True)
    p.add_argument('--index',type=int,required=True)
    p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--reads-root',type=Path,required=True)
    p.add_argument('--output-root',type=Path,required=True)
    p.add_argument('--threads',type=int,default=4)
    a = p.parse_args()
    run(a.cohort,a.index,a.reference,a.reads_root,a.output_root,a.threads)
