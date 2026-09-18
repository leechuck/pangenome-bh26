"""Compose frozen HGSVC follow-up calls without accessing evaluation truth."""
import argparse
import csv
import json
import os
from pathlib import Path
from evidence_io import sha
from decode_t1k import decode
from anchor_t1k_coarse import anchor

ROOT = Path('/home/leechuck/hla/t1k-pangenome/hgsvc-validation-v1')


def run(index, assets):
    if not os.environ.get('SLURM_CPUS_PER_TASK'):raise ValueError('Requires Slurm')
    cohort = list(csv.DictReader((assets/'cohort.tsv').open(),delimiter='\t'))
    if not 0 <= index < len(cohort)*3:raise ValueError('Invalid index')
    i,m = divmod(index,3)
    method = ('ipd_genome','graph_hprc','graph_hprc_asian')[m]
    donor = cohort[i]['donor']
    original = ROOT/'validation-runs/frozen-t1k-v1'/donor
    baseline = ROOT/'validation-runs/genomic-ipd-v1'/donor
    candidate = baseline if m==0 else ROOT/'validation-runs/graph-v1/calls'/donor/('hprc','hprc_asian')[m-1]/'unsampled'
    om = json.loads((original/'COMPLETE.json').read_text())
    cm = json.loads((candidate/'COMPLETE.json').read_text())
    if om['status']!='complete' or cm['status']!='complete':raise ValueError('Incomplete inputs')
    table = original/(donor+'_genotype.tsv')
    if sha(table)!=om['output_sha256'][table.name]:raise ValueError('Original calls changed')
    if sha(candidate/'calls.json')!=cm['output_sha256']['calls.json']:raise ValueError('Candidate calls changed')
    if sorted(om['reads_sha256'].values())!=sorted(cm['reads_sha256'].values()):raise ValueError('Different reads')
    native = decode(table.read_text(),{})
    calls,decisions = anchor(native,json.loads((candidate/'calls.json').read_text()))
    output = ROOT/'validation-runs/anchored-v1'/method/donor
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    (output/'calls.json').write_text(json.dumps(calls,indent=2)+'\n')
    (output/'decisions.json').write_text(json.dumps(decisions,indent=2)+'\n')
    record = dict(status='complete',donor=donor,method=method,
        reads_sha256=om['reads_sha256'],original_sha256=sha(original/'COMPLETE.json'),
        candidate_sha256=sha(candidate/'COMPLETE.json'),driver_sha256=sha(Path(__file__)),
        anchor_sha256=sha(Path(__file__).with_name('anchor_t1k_coarse.py')),
        output_sha256={n:sha(output/n) for n in ('calls.json','decisions.json')})
    (output/'COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index',type=int,required=True)
    parser.add_argument('--assets',type=Path,required=True)
    args=parser.parse_args();run(args.index,args.assets)
