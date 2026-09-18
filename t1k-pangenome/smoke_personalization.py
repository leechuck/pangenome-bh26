"""Exercise personalized graph mapping on synthetic reads from training paths.

This checks software integration only, never benchmark accuracy or generalization.
"""
import argparse
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import time
from build_graph import sequences, sha


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def smoke(graph, reference, output, threads):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Smoke test requires a Slurm allocation')
    if output.exists():
        raise FileExistsError(output)
    graph, reference = graph.resolve(), reference.resolve()
    os.environ['OMP_NUM_THREADS'] = str(threads)
    os.environ['PATH'] = '/home/leechuck/hla/t1k-pangenome/tools/kmc-3.2.4:' + os.environ['PATH']
    kmc = shutil.which('kmc')
    if not kmc:
        raise FileNotFoundError('KMC is required by vg integrated k-mer counting')
    manifest = json.loads((graph/'COMPLETE.json').read_text())
    if manifest['status'] != 'complete' or manifest['source_sha256'] != sha(reference):
        raise ValueError('Graph/reference mismatch')
    for name, digest in manifest['output_sha256'].items():
        if sha(graph/name) != digest:
            raise ValueError('Graph product changed: ' + name)
    paths = sequences(reference)
    chosen = sorted(name for name, seq in paths.items() if 'N' not in seq)[:2]
    if len(chosen) != 2:
        raise ValueError('Need two unambiguous training sequences')
    output.mkdir(parents=True)
    output = output.resolve()
    rng = random.Random(20260918)
    pair_count = 0
    with (output/'r1.fq').open('w') as r1, (output/'r2.fq').open('w') as r2:
        for name in chosen:
            seq = paths[name]
            for _ in range(40*len(seq)//300):
                start = rng.randrange(len(seq)-350+1)
                fragment = seq[start:start+350]
                for mate, read, stream in ((1, fragment[:150], r1),
                                           (2, reverse_complement(fragment[-150:]), r2)):
                    stream.write(f'@synthetic{pair_count}/{mate}\n{read}\n+\n'+('I'*150)+'\n')
                pair_count += 1
    # Isolate automatically generated indexes from the verified source graph.
    for name in ('graph.gbz', 'graph.hapl'):
        shutil.copy2(graph/name, output/name)
    vg = '/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin/vg'
    command = [vg, 'giraffe', '--haplotype-sampling', '--no-diploid-sampling',
               '--num-haplotypes', '8', '-Z', 'graph.gbz', '--haplotype-name', 'graph.hapl',
               '-f', 'r1.fq', '-f', 'r2.fq', '-t', str(threads), '-N', 'synthetic']
    record = dict(status='running', started=time.time(), command=command,
                  scope='Synthetic training-path integration test; not typing accuracy',
                  graph_manifest_sha256=sha(graph/'COMPLETE.json'),
                  driver_sha256=sha(Path(__file__)), vg_sha256=sha(Path(vg)),
                  kmc_sha256=sha(Path(kmc)),
                  selected_training_paths=chosen, pairs=pair_count, seed=20260918,
                  reads_sha256={p:sha(output/p) for p in ('r1.fq','r2.fq')})

    def save():
        (output/'manifest.json').write_text(json.dumps(record, indent=2)+'\n')

    save()
    try:
        with (output/'mapping.gam').open('wb') as gam, (output/'mapping.log').open('wb') as log:
            subprocess.run(command, cwd=output, stdout=gam, stderr=log, check=True)
        with (output/'alignments.jsonl').open('wb') as stream:
            subprocess.run([vg,'view','-aj',str(output/'mapping.gam')], stdout=stream, check=True)
        count = mapped = 0
        with (output/'alignments.jsonl').open() as stream:
            for line in stream:
                alignment = json.loads(line)
                count += 1
                mapped += bool(alignment.get('path', {}).get('mapping'))
        record.update(alignments=count, mapped=mapped)
        if count != 2*pair_count or mapped/count < 0.95:
            raise ValueError('Synthetic mapping completeness gate failed')
        record.update(status='complete', finished=time.time(), gam_sha256=sha(output/'mapping.gam'))
        save()
        (output/'COMPLETE.json').write_text(json.dumps(record, indent=2)+'\n')
    except Exception as error:
        record.update(status='failed', finished=time.time(), error=repr(error))
        save()
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--graph', type=Path, required=True)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--threads', type=int, default=8)
    a = p.parse_args()
    smoke(a.graph, a.reference, a.output, a.threads)
