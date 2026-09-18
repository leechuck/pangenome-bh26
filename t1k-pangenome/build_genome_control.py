"""Replace truncated T1K entries with available full genomic IPD sequences.

Retain every original candidate, including partial-reference fallbacks. This is
an IPD-information control, not a graph method or Asian population prior.
"""
import argparse
import collections
import json
import os
from pathlib import Path
import subprocess
from build_graph import sha


def records(text, validate_coordinates=True):
    result = {}
    name,header,parts = None,None,[]
    def save():
        if name is None:
            return
        if name in result:
            raise ValueError('Duplicate allele: '+name)
        sequence = ''.join(parts).upper()
        fields = header.split()
        count = int(fields[1])
        coordinates = list(map(int,fields[2:]))
        if not sequence or count<1 or len(coordinates)!=2*count:
            raise ValueError('Incomplete reference annotation: '+name)
        previous = -1
        for start,end in zip(coordinates[::2],coordinates[1::2]):
            if validate_coordinates and not previous < start <= end < len(sequence):
                raise ValueError('Invalid reference coordinates: '+name)
            previous = end
        result[name] = (header,sequence)
    for line in text.splitlines():
        if line.startswith('>'):
            save()
            header = line[1:]
            name,parts = header.split()[0],[]
        else:
            parts.append(line.strip())
    save()
    if not result:
        raise ValueError('Empty reference')
    return result


def invalid_annotations(data):
    invalid = {}
    for name,(header,sequence) in data.items():
        try:
            records('>'+header+'\n'+sequence+'\n')
        except ValueError as error:
            invalid[name] = str(error)
    return invalid


def combine(original, genomic):
    result = dict(original)
    replaced = []
    for name,record in genomic.items():
        if name in result:
            result[name] = record
            replaced.append(name)
    return result,sorted(replaced)


def build(ipd, baseline, parser, output):
    if not os.environ.get('SLURM_CPUS_PER_TASK'):
        raise RuntimeError('Reference parsing requires Slurm')
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    command = ['perl',str(parser),str(ipd),'--mode','genome','--gene','HLA']
    with (output/'genomic.fa').open('w') as target, (output/'parse.log').open('w') as log:
        subprocess.run(command,stdout=target,stderr=log,check=True)
    # Preserve legacy entries exactly, including inherited annotation defects.
    # Never introduce an invalid genomic replacement to "repair" them silently.
    original = records(baseline.read_text(),validate_coordinates=False)
    genomic = records((output/'genomic.fa').read_text(),validate_coordinates=False)
    original_invalid = invalid_annotations(original)
    genomic_invalid = invalid_annotations(genomic)
    genomic = {name:record for name,record in genomic.items() if name not in genomic_invalid}
    combined,replaced = combine(original,genomic)
    if not replaced or set(combined)!=set(original):
        raise ValueError('Candidate-preserving genome replacement failed')
    with (output/'reference.fa').open('w') as target:
        for header,sequence in combined.values():
            target.write('>'+header+'\n'+sequence+'\n')
    (output/'aliases.json').write_text('{}\n')
    report = dict(status='complete',panel='ipd_genome',ipd_sha256=sha(ipd),
                  baseline_sha256=sha(baseline),parser_sha256=sha(parser),
                  driver_sha256=sha(Path(__file__)),command=command,
                  candidates=len(original),genomic_replacements=len(replaced),
                  retained_original_candidates=len(original)-len(replaced),
                  unused_genomic_candidates=len(set(genomic)-set(original)),
                  original_annotation_exceptions=original_invalid,
                  rejected_genomic_annotations=genomic_invalid,
                  replacements_by_gene=dict(collections.Counter(name.split('*')[0] for name in replaced)),
                  output_sha256={name:sha(output/name) for name in ('reference.fa','aliases.json','genomic.fa')},
                  scope='Full genomic IPD information control; all original candidate IDs retained; no graph or new imputation')
    (output/'COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('ipd','baseline','parser','output'):
        p.add_argument('--'+name,type=Path,required=True)
    a = p.parse_args()
    build(a.ipd,a.baseline,a.parser,a.output)
