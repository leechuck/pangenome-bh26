"""Estimate locus k-mer depth using conserved markers screened against IPD.

Specificity is relative to supplied HLA/decoy references, not the whole genome.
Coverage is a diagnostic for personalization, never an allele call.
"""
import argparse
import collections
import gzip
import json
import math
import os
from pathlib import Path
import statistics
from build_graph import sequences, sha
from build_panel import GENES


def kmers(sequence, k):
    if k < 1:
        raise ValueError('Positive k required')
    mask, shift = (1 << (2*k))-1, 2*(k-1)
    forward = reverse = length = 0
    bases = {'A':0,'C':1,'G':2,'T':3}
    for base in sequence.upper():
        value = bases.get(base)
        if value is None:
            forward = reverse = length = 0
            continue
        forward = ((forward << 2)|value)&mask
        reverse = (reverse >> 2)|((3-value) << shift)
        length += 1
        if length >= k:
            yield min(forward,reverse)


def conserved_markers(paths, k, fraction):
    if not paths or not 0 < fraction <= 1:
        raise ValueError('Nonempty paths and conservation fraction in (0,1] required')
    occurrences, repeated = collections.Counter(), set()
    for sequence in paths:
        counts = collections.Counter(kmers(sequence,k))
        occurrences.update(counts.keys())
        repeated.update(key for key,count in counts.items() if count > 1)
    threshold = math.ceil(fraction*len(paths))
    return {key for key,count in occurrences.items() if count >= threshold and key not in repeated}


def build(reference, ipd, output, k=29, fraction=.9):
    if output.exists():
        raise FileExistsError(output)
    manifest = json.loads((reference/'COMPLETE.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('Incomplete panel')
    owners = collections.defaultdict(set)
    before = {}
    for gene in GENES:
        name = 'hprc_asian/HLA-'+gene+'.fa'
        if sha(reference/name) != manifest['output_sha256'][name]:
            raise ValueError('Reference changed')
        markers = conserved_markers(list(sequences(reference/name).values()),k,fraction)
        before[gene] = len(markers)
        for key in markers:
            owners[key].add(gene)
    # Check every panel path, including markers not conserved in the other locus.
    for gene in GENES:
        for sequence in sequences(reference/'hprc_asian'/('HLA-'+gene+'.fa')).values():
            for key in kmers(sequence,k):
                if key in owners:
                    owners[key].add(gene)
    for name,sequence in sequences(ipd).items():
        gene = name.split('*')[0].removeprefix('HLA-')
        for key in kmers(sequence,k):
            if key in owners:
                owners[key].add(gene)
    by_gene = {gene:sorted(format(key,'x') for key,genes in owners.items() if genes=={gene}) for gene in GENES}
    report = dict(status='complete',k=k,conservation_fraction=fraction,
                  markers=by_gene,counts_before_specificity=before,
                  counts={gene:len(values) for gene,values in by_gene.items()},
                  reference_manifest_sha256=sha(reference/'COMPLETE.json'),ipd_sha256=sha(ipd),
                  driver_sha256=sha(Path(__file__)),
                  specificity_scope='Observed HLA panels and supplied all-IPD reference only; not genome-wide')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({key:value for key,value in report.items() if key!='markers'},indent=2))


def estimate(counts, minimum_markers=100, minimum_depth=20):
    depth = statistics.median(counts) if counts else 0
    zero_fraction = counts.count(0)/len(counts) if counts else 1
    usable = len(counts)>=minimum_markers and depth>=minimum_depth and zero_fraction<=.25
    return dict(markers=len(counts),median_kmer_depth=depth,zero_fraction=zero_fraction,
                usable=usable,personalization_coverage=int(round(depth)) if usable else None,
                fallback_reason=None if usable else 'Insufficient or inconsistent conserved-marker coverage')


def measure(markers, reads, output):
    if output.exists():
        raise FileExistsError(output)
    reference = json.loads(markers.read_text())
    if reference['status'] != 'complete':
        raise ValueError('Marker preparation incomplete')
    target = {int(value,16):0 for values in reference['markers'].values() for value in values}
    read_count = 0
    for path in reads:
        opener = gzip.open if path.suffix=='.gz' else open
        with opener(path,'rt') as stream:
            while True:
                name = stream.readline()
                if not name:
                    break
                sequence,plus,quality = [stream.readline().rstrip('\n') for _ in range(3)]
                if not name.startswith('@') or not plus.startswith('+') or not sequence or len(sequence)!=len(quality):
                    raise ValueError('Malformed FASTQ')
                read_count += 1
                for key in kmers(sequence,reference['k']):
                    if key in target:
                        target[key] += 1
    result = dict(status='complete',marker_sha256=sha(markers),driver_sha256=sha(Path(__file__)),
                  reads_sha256={str(p):sha(p) for p in reads},reads=read_count,
                  genes={gene:estimate([target[int(value,16)] for value in values])
                         for gene,values in reference['markers'].items()},
                  specificity_scope=reference['specificity_scope'],
                  scope='Coverage diagnostic only; low support must retain unsampled graph/all-IPD fallback')
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='mode',required=True)
    b = sub.add_parser('build')
    for name in ('reference','ipd','output'):
        b.add_argument('--'+name,type=Path,required=True)
    m = sub.add_parser('measure')
    for name in ('markers','output'):
        m.add_argument('--'+name,type=Path,required=True)
    m.add_argument('--reads',type=Path,nargs='+',required=True)
    a = p.parse_args()
    if not os.environ.get('SLURM_CPUS_PER_TASK'):
        p.error('Reference screening and read counting require Slurm')
    if a.mode=='build':
        build(a.reference,a.ipd,a.output)
    else:
        measure(a.markers,a.reads,a.output)
