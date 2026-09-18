"""All-locus development mapping using full-genomic-IPD-screened coverage."""
import argparse
from pathlib import Path
from map_personalized import run

GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')
ROOT = Path('/home/leechuck/hla/t1k-pangenome')


def graph_path(panel, gene):
    if panel == 'hprc' and gene == 'A':
        version = 'pilot-v2'
    elif gene in ('DQA1', 'DRB1'):
        version = 'build-v2'
    else:
        version = 'build-v1'
    return ROOT/'graphs'/version/panel/gene


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--panel', choices=('hprc', 'hprc_asian'), required=True)
    parser.add_argument('--selection', choices=('sampled', 'unsampled'), required=True)
    parser.add_argument('--index', type=int, choices=range(len(GENES)), required=True)
    parser.add_argument('--calibration', type=Path)
    args = parser.parse_args()
    gene = GENES[args.index]
    reads = Path('/home/leechuck/hla/hla-typer/reads/HG00658')
    run(graph_path(args.panel, gene),
        ROOT/'coverage/development-genome-screen-v1/HG00658.json', gene,
        [reads/'r1.fq.gz', reads/'r2.fq.gz'],
        ROOT/'development'/('mapping-v3' if args.calibration else 'mapping-v2')/'HG00658'/args.panel/args.selection/gene,
        threads=4, disable_selection=args.selection == 'unsampled',
        calibration_file=args.calibration)
