"""Build one observed-panel graph and require vg haplotype-sampling preprocessing.

Reference preparation only: no sample reads, allele calls or evaluation truth.
Run inside the existing IBEX container under a Slurm allocation.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sequences(path):
    result, name, parts = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            if name is not None:
                if name in result:
                    raise ValueError('Duplicate path: ' + name)
                result[name] = ''.join(parts).upper()
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    if name is not None:
        if name in result:
            raise ValueError('Duplicate path: ' + name)
        result[name] = ''.join(parts).upper()
    if not result or any(not s for s in result.values()):
        raise ValueError('Empty graph input/path')
    return result


def verify_paths(expected, actual):
    # GBZ exports haplotype paths with an additional phase-block field (#0).
    # Accept only that explicit representation change, with a bijective map.
    aliases = {}
    normalized = {}
    for name, sequence in actual.items():
        target = name if name in expected else name[:-2] if name.endswith('#0') else name
        if target not in expected or target in normalized:
            raise ValueError('Unknown or duplicate exported path: ' + name)
        aliases[name] = target
        normalized[target] = sequence
    if expected != normalized:
        raise ValueError('Graph path names/sequences differ from input: '
                         f'{len(expected)} expected, {len(actual)} observed')
    return aliases


def build(reference, panel, gene, output, threads):
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Graph construction requires a Slurm allocation')
    if output.exists():
        raise FileExistsError(output)
    reference = reference.resolve()
    manifest = json.loads((reference / 'COMPLETE.json').read_text())
    relative = f'{panel}/HLA-{gene}.fa'
    source = reference / relative
    if manifest['status'] != 'complete' or sha(source) != manifest['output_sha256'][relative]:
        raise ValueError('Reference manifest/input mismatch')
    expected = sequences(source)
    if len(expected) < 2 or min(map(len, expected.values())) < 2000:
        raise ValueError('Pilot alignment parameters require >=2 paths of >=2000 bp')
    output.mkdir(parents=True)
    output = output.resolve()
    os.environ['PATH'] = '/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin:' \
        '/home/leechuck/hla/mm/envs/pggb053/bin:' + os.environ['PATH']
    tools = ('bgzip', 'samtools', 'wfmash', 'seqwish', 'smoothxg', 'gfaffix', 'odgi', 'vg')
    record = dict(status='running', started=time.time(), panel=panel, gene=gene,
                  job=os.environ['SLURM_JOB_ID'], reference_manifest_sha256=sha(reference/'COMPLETE.json'),
                  source_sha256=sha(source), builder_sha256=sha(Path(__file__)),
                  parameters=dict(segment_length=500, minimum_block_length=1000, identity=90),
                  commands=[], tools={})

    def save():
        (output/'manifest.json').write_text(json.dumps(record, indent=2)+'\n')

    def run(args, stdout=None):
        args = list(map(str, args))
        record['commands'].append(args)
        save()
        with (output/'build.log').open('ab') as log:
            if stdout:
                with (output/stdout).open('wb') as out:
                    subprocess.run(args, cwd=output, stdout=out, stderr=log, check=True)
            else:
                subprocess.run(args, cwd=output, stdout=log, stderr=log, check=True)

    save()
    try:
        for tool in tools:
            path = shutil.which(tool)
            if not path:
                raise FileNotFoundError(tool)
            record['tools'][tool] = dict(path=path, sha256=sha(Path(path)))
        n = len(expected)
        run(['bgzip', '-c', '-@', threads, source], 'input.fa.gz')
        run(['samtools', 'faidx', 'input.fa.gz'])
        common = ['wfmash', '-s', 500, '-l', 1000, '-p', 90, '-n', n-1,
                  '-k', 19, '-H', 0.001, '-X', '-t', threads, '--tmp-base', output]
        run(common + ['input.fa.gz', '--approx-map'], 'mappings.paf')
        run(common + ['input.fa.gz', '-i', 'mappings.paf', '--invert-filtering'], 'alignments.paf')
        if not (output/'alignments.paf').stat().st_size:
            raise ValueError('No cross-path alignments')
        run(['seqwish', '-s', 'input.fa.gz', '-p', 'alignments.paf', '-k', 19, '-f', 0,
             '-g', 'raw.gfa', '-B', 10000000, '-t', threads, '--temp-dir', output, '-P'])
        run(['smoothxg', '-t', threads, '-T', threads, '-g', 'raw.gfa', '-r', n,
             '--base', output, '--chop-to', 100, '-I', 0.9, '-R', 0, '-j', 0, '-e', 0,
             '-l', '700,900,1100', '-p', '1,19,39,3,81,1', '-O', 0.001,
             '-Y', 100*n, '-d', 0, '-D', 0, '-V', '-o', 'smooth.gfa'])
        run(['gfaffix', 'smooth.gfa', '-o', 'fixed.gfa'], 'affixes.tsv')
        run(['odgi', 'build', '-t', threads, '-g', 'fixed.gfa', '-o', 'built.og', '-O'])
        run(['odgi', 'unchop', '-t', threads, '-i', 'built.og', '-o', 'unchopped.og'])
        run(['odgi', 'sort', '-p', 'Ygs', '--temp-dir', output, '-t', threads,
             '-i', 'unchopped.og', '-o', 'sorted.og'])
        run(['odgi', 'view', '-i', 'sorted.og', '-g'], 'graph.gfa')
        run(['vg', 'gbwt', '-G', 'graph.gfa', '--gbz-format', '-g', 'raw.gbz', '--num-threads', threads])
        # Deterministic indexing backbone only; this is not a sample genotype.
        backbone = max(expected, key=lambda name:(len(expected[name]), name)).split('#')[0]
        record['indexing_backbone'] = backbone
        run(['vg', 'gbwt', '-Z', 'raw.gbz', '--set-reference', backbone,
             '--gbz-format', '-g', 'graph.gbz'])
        run(['vg', 'paths', '-x', 'graph.gbz', '-F'], 'paths.fa')
        record['exported_path_names'] = verify_paths(expected, sequences(output/'paths.fa'))
        run(['vg', 'index', '-j', 'graph.dist', '--no-nested-distance', '-t', threads, 'graph.gbz'])
        run(['vg', 'gbwt', '-p', '--num-threads', threads, '-r', 'graph.ri', '-Z', 'graph.gbz'])
        run(['vg', 'haplotypes', '-v', 2, '-t', threads, '-H', 'graph.hapl', 'graph.gbz'])
        products = ('graph.gbz', 'graph.dist', 'graph.ri', 'graph.hapl')
        if any(not (output/p).stat().st_size for p in products):
            raise ValueError('Empty sampling index')
        record.update(status='complete', finished=time.time(), verified_paths=n,
                      output_sha256={p:sha(output/p) for p in products})
        save()
        (output/'COMPLETE.json').write_text(json.dumps(record, indent=2)+'\n')
    except Exception as error:
        record.update(status='failed', finished=time.time(), error=repr(error))
        save()
        raise


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--panel', choices=('hprc', 'hprc_asian'), required=True)
    p.add_argument('--gene', choices=('A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1'), required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--threads', type=int, default=8)
    a = p.parse_args()
    if a.threads < 1:
        p.error('--threads must be positive')
    build(a.reference, a.panel, a.gene, a.output, a.threads)
