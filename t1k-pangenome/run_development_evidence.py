"""Project development graph alignments and retain paired-fragment evidence."""
import argparse
import json
import os
from pathlib import Path
from evidence_io import sha
from fragment_support import build
from path_support import project
from run_development_mapping import GENES, ROOT, graph_path


def smoke():
    output = ROOT/'smoke/fragment-support-disk-v1'
    build(ROOT/'smoke/path-support-v1', output, 1000)
    previous = ROOT/'smoke/fragment-support-v2/fragments.jsonl'
    if sha(output/'fragments.jsonl') != sha(previous):
        raise ValueError('Disk fragment grouping changed synthetic evidence')
    report = dict(status='complete', previous_sha256=sha(previous),
                  current_sha256=sha(output/'fragments.jsonl'),
                  scope='Byte-identical fragment evidence on the synthetic integration library')
    (output/'PARITY.json').write_text(json.dumps(report, indent=2)+'\n')


def run(panel, selection, index):
    gene = GENES[index]
    mapping = ROOT/'development/mapping-v2/HG00658'/panel/selection/gene
    graph = graph_path(panel, gene)
    mm = json.loads((mapping/'COMPLETE.json').read_text())
    gm = json.loads((graph/'COMPLETE.json').read_text())
    if mm['status'] != 'complete' or gm['status'] != 'complete':
        raise ValueError('Incomplete mapping or graph')
    if mm['graph_manifest_sha256'] != sha(graph/'COMPLETE.json'):
        raise ValueError('Mapping graph provenance changed')
    source = Path(gm['previous_attempt']) if 'previous_attempt' in gm else graph
    if source != graph and sha(source/'manifest.json') != gm['previous_manifest_sha256']:
        raise ValueError('Original graph construction manifest changed')
    if sha(source/'graph.gbz') != gm['output_sha256']['graph.gbz']:
        raise ValueError('Path source and mapping graph differ')
    if sha(mapping/'alignments.jsonl') != mm['output_sha256']['alignments.jsonl']:
        raise ValueError('Alignment evidence changed')
    output = ROOT/'development/fragments-v1/HG00658'/panel/selection/gene
    project(source/'graph.gfa', mapping/'alignments.jsonl', output/'support')
    build(output/'support', output/'fragments', 1000)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--panel', choices=('hprc','hprc_asian'))
    parser.add_argument('--selection', choices=('sampled','unsampled'))
    parser.add_argument('--index', type=int, choices=range(len(GENES)))
    args = parser.parse_args()
    if not os.environ.get('SLURM_CPUS_PER_TASK'):
        parser.error('Evidence processing requires Slurm')
    if args.smoke:
        smoke()
    elif args.panel is None or args.selection is None or args.index is None:
        parser.error('Panel, selection and index required for development evidence')
    else:
        run(args.panel,args.selection,args.index)
