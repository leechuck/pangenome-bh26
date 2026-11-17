"""Summarise SCAN.json files (counts, contig names only).  python scan_table.py SCAN_DIR"""
import json
import sys
from pathlib import Path
print('| donor | copy | IPD-absent | k-mer col range | own reads (all, primary) | own reads in hole | pairs in hole | in old recruited fastq | pulled by old rule | contig class | MAPQ | proper pair |')
print('|---|---|---|---|---:|---:|---:|---:|---:|---|---|---:|')
for p in sorted(Path(sys.argv[1]).glob('*/SCAN.json')):
    a = json.loads(p.read_text())
    for c in a['copies']:
        print(f"| {a['donor']} | {c['hap']} | {c['ipd_absent']} | {a['col_range'][0]}-{a['col_range'][1]} | {c['own_reads_all']} | "
              f"{c['own_reads_in_hole']} | {c['own_pairs_in_hole']} | {c['in_hole_in_recruited_fastq']} | {c['in_hole_pulled_by_rule']} | "
              f"{'; '.join(f'{k} {v}' for k, v in c['in_hole_contig_class'].items())} | {'; '.join(f'{k}:{v}' for k, v in sorted(c['in_hole_mapq'].items()))} | {c['in_hole_proper_pair']} |")
