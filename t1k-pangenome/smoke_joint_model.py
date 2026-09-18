"""Fit paired synthetic training-path evidence; not an accuracy benchmark."""
import argparse
import hashlib
import json
import os
from pathlib import Path
from build_graph import sequences, sha
from joint_model import emission_matrix, fit


def run(source, reference, output):
    if not os.environ.get('SLURM_CPUS_PER_TASK'):
        raise RuntimeError('Run under Slurm')
    if output.exists():
        raise FileExistsError(output)
    manifest = json.loads((source/'COMPLETE.json').read_text())
    evidence = source/'fragments.jsonl'
    if manifest['status'] != 'complete' or sha(evidence) != manifest['output_sha256']:
        raise ValueError('Fragment evidence changed')
    seqs = sequences(reference)
    metadata = json.loads(reference.with_suffix('.json').read_text())
    labels = {record['path']: record for record in metadata}
    if set(labels) != set(seqs):
        raise ValueError('Reference metadata/path mismatch')
    for path, sequence in seqs.items():
        if hashlib.sha256(sequence.encode()).hexdigest() != labels[path]['sequence_sha256']:
            raise ValueError('Reference sequence changed')
    fragments = [json.loads(line) for line in evidence.read_text().splitlines()]
    paths = sorted(seqs)
    # Known synthetic library insert length; this is not a real-library estimator.
    lengths = {path:len(seqs[path])-350+1 for path in paths}
    matrix = emission_matrix(fragments, paths, lengths)
    result = fit(matrix, ['A']*len(paths))
    for gene_result in result['genes'].values():
        gene_result['path_pairs'] = [[paths[i] for i in pair] for pair in gene_result.pop('pairs')]
        gene_result['genomic_label_pairs'] = [
            [labels[path]['genomic_labels'] for path in pair]
            for pair in gene_result['path_pairs']]
    report = dict(status='complete', result=result, fragments=len(fragments),
                  candidate_paths=len(paths), effective_insert=350, temperature=10.0,
                  source_manifest_sha256=sha(source/'COMPLETE.json'),
                  reference_sha256=sha(reference), metadata_sha256=sha(reference.with_suffix('.json')),
                  model_sha256=sha(Path(__file__).with_name('joint_model.py')),
                  driver_sha256=sha(Path(__file__)),
                  scope='Synthetic training-path integration only; no all-IPD fallback or cross-locus decoys')
    if not result['converged']:
        raise RuntimeError('Synthetic joint fit did not converge')
    output.mkdir(parents=True)
    (output/'COMPLETE.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.source, a.reference, a.output)
