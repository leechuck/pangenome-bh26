"""Retry sampling-index construction with an explicit indexing backbone."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from build_graph import sha


def repair(source, output, threads):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Requires Slurm allocation')
    if output.exists():
        raise FileExistsError(output)
    source = source.resolve()
    previous = json.loads((source/'manifest.json').read_text())
    if previous['status'] != 'failed' or 'top-level chain 0 is a loop' not in (source/'build.log').read_text():
        raise ValueError('Not the diagnosed top-level loop failure')
    if not previous.get('exported_path_names'):
        raise ValueError('No successful input-path verification')
    output.mkdir(parents=True)
    output = output.resolve()
    record = dict(previous)
    record.update(status='running', started=time.time(), job=os.environ['SLURM_JOB_ID'],
                  repair_driver_sha256=sha(Path(__file__)),
                  previous_manifest_sha256=sha(source/'manifest.json'), previous_attempt=str(source),
                  repair='Explicit path prefix for distance-index backbone', commands=[])
    record.pop('error', None)
    record.pop('finished', None)
    os.environ['OMP_NUM_THREADS'] = str(threads)
    vg = '/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin/vg'
    if sha(Path(vg)) != previous['tools']['vg']['sha256']:
        raise ValueError('vg binary changed')

    def save():
        (output/'manifest.json').write_text(json.dumps(record, indent=2)+'\n')

    save()
    try:
        for name in ('graph.gbz', 'graph.ri'):
            shutil.copy2(source/name, output/name)
            if sha(source/name) != sha(output/name):
                raise ValueError('Copy mismatch: ' + name)
        commands = [[vg,'index','-j','graph.dist','--no-nested-distance','-t',str(threads),
                     '-P',previous['indexing_backbone'],'graph.gbz'],
                    [vg,'haplotypes','-v','2','-t',str(threads),'-H','graph.hapl','graph.gbz']]
        for command in commands:
            record['commands'].append(command);save()
            with (output/'build.log').open('ab') as log:
                subprocess.run(command, cwd=output, stdout=log, stderr=log, check=True)
        products = ('graph.gbz','graph.ri','graph.dist','graph.hapl')
        if any(not (output/p).stat().st_size for p in products):
            raise ValueError('Empty index')
        record.update(status='complete', finished=time.time(),
                      verified_paths=len(previous['exported_path_names']),
                      output_sha256={p:sha(output/p) for p in products})
        save()
        (output/'COMPLETE.json').write_text(json.dumps(record, indent=2)+'\n')
    except Exception as error:
        record.update(status='failed', finished=time.time(), error=repr(error));save()
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--threads', type=int, default=8)
    a = p.parse_args()
    repair(a.source, a.output, a.threads)
