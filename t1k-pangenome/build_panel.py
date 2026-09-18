"""Build observed, family-excluded HLA path references for T1K graph experiments.

No read calling or truth scoring. An exact CDS match is not a four-field label.
Identical sequences are stored once, with all allowed source paths in a sidecar.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

GENES = ('A', 'B', 'C', 'DPA1', 'DPB1', 'DQA1', 'DQB1', 'DRB1')
FLANK = 2000


def digest(data):
    return hashlib.sha256(data).hexdigest()


def rows(path):
    with path.open() as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def fasta(path, allowed):
    name, seq = None, []
    with path.open() as stream:
        for line in stream:
            if line.startswith('>'):
                if name in allowed:
                    yield name, ''.join(seq).upper()
                name = line[1:].split()[0]
                seq = []
            elif name in allowed:
                seq.append(line.strip())
    if name in allowed:
        yield name, ''.join(seq).upper()


def validate_sequence(name, sequence, annotation, exons):
    if len(sequence) <= 2 * FLANK or set(sequence) - set('ACGTN'):
        raise ValueError('Invalid flanked sequence: ' + name)
    body = sequence[FLANK:-FLANK]
    if digest(body.encode()) != annotation['gene_sha256']:
        raise ValueError('Gene sequence does not match annotation: ' + name)
    intervals = sorted(exons)
    for i, (start, end) in enumerate(intervals):
        if not FLANK <= start < end <= len(sequence) - FLANK:
            raise ValueError('Coding interval outside gene: ' + name)
        if i and start < intervals[i-1][1]:
            raise ValueError('Overlapping coding intervals: ' + name)
    if intervals:
        cds = ''.join(sequence[a:b] for a, b in intervals)
        if digest(cds.encode()) != annotation['cds_sha256']:
            raise ValueError('Coding sequence does not match annotation: ' + name)
    return intervals


def group_paths(records):
    grouped = {}
    for record in records:
        sequence = record['sequence']
        key = digest(sequence.encode())
        if key not in grouped:
            grouped[key] = dict(sequence=sequence, sequence_sha256=key,
                                sources=[], genomic_labels=set(), cds_labels=set(), exons=set())
        group = grouped[key]
        group['sources'].append(record['name'])
        group['genomic_labels'].update(record['genomic_labels'])
        group['cds_labels'].update(record['cds_labels'])
        if record['exons']:
            group['exons'].add(tuple(record['exons']))
    for group in grouped.values():
        group['sources'].sort()
        if len(group['exons']) > 1:
            raise ValueError('Identical sequences have inconsistent coding annotations')
        group['exons'] = list(next(iter(group['exons']), ()))
        group['genomic_labels'] = sorted(group['genomic_labels'])
        group['cds_labels'] = sorted(group['cds_labels'])
    return sorted(grouped.values(), key=lambda g:g['sequence_sha256'])


def build(source, metadata_path, catalogue_path, exons_path, reservation, output):
    if output.exists():
        raise FileExistsError('Reference output already exists: ' + str(output))
    lock = json.loads((reservation/'RESERVATION.json').read_text())
    for file, field in [('cohort.tsv', 'cohort_sha256'), ('excluded_paths.tsv', 'excluded_paths_sha256')]:
        if digest((reservation/file).read_bytes()) != lock[field]:
            raise ValueError('Reservation changed: ' + file)
    if digest(metadata_path.read_bytes()) != lock['source_sha256']['metadata']:
        raise ValueError('Donor metadata changed since reservation')
    excluded = {r['hap_id'] for r in rows(reservation/'excluded_paths.tsv')}
    metadata = {r['hap_id']:r for r in rows(metadata_path)}
    allowed_haps = set(metadata) - excluded
    annotations = {}
    for row in rows(catalogue_path):
        if row['hap_id'] not in allowed_haps:
            continue
        annotations[row['name']] = row
    exons = collections.defaultdict(list)
    for row in rows(exons_path):
        if row['name'] in annotations:
            exons[row['name']].append((int(row['start0']), int(row['end0'])))
    output.mkdir(parents=True)
    manifest = dict(status='building', exclusion_count=len(excluded),
                    builder_sha256=digest(Path(__file__).read_bytes()), flank_bases=FLANK,
                    observation_policy='Assembly sequences only; no imputed bases; N positions remain unknown.',
                    labels='Exact genomic and exact CDS labels kept separate; no new allele names assigned.',
                    input_sha256={str(p):digest(p.read_bytes()) for p in
                                  (metadata_path,catalogue_path,exons_path,reservation/'RESERVATION.json')},
                    genes={}, output_sha256={})
    for gene in GENES:
        source_fasta = source/f'HLA-{gene}.fa'
        manifest['input_sha256'][str(source_fasta)] = digest(source_fasta.read_bytes())
        records = []
        for name, seq in fasta(source_fasta, annotations):
            ann = annotations[name]
            if ann['gene'] != 'HLA-' + gene:
                raise ValueError('Wrong locus annotation: ' + name)
            hap = '#'.join(name.split('#')[:2])
            if hap not in allowed_haps or hap != ann['hap_id']:
                raise ValueError('Excluded or inconsistent path: ' + name)
            iv = validate_sequence(name, seq, ann, exons[name])
            records.append(dict(name=name, sequence=seq, exons=iv, cohort=metadata[hap]['cohort'],
                                genomic_labels=[n for n in ann['exact_genomic_alleles'].split(';') if n],
                                cds_labels=[n for n in ann['exact_cds_alleles'].split(';') if n]))
        if not records:
            raise ValueError('Empty panel: ' + gene)
        for panel in ('hprc', 'hprc_asian'):
            selected = [r for r in records if panel == 'hprc_asian'
                        or r['cohort'].startswith('HPRC') or r['cohort'] == 'REF']
            groups = group_paths(selected)
            if not groups:
                raise ValueError('Empty panel: ' + panel + '/' + gene)
            directory = output/panel
            directory.mkdir(exist_ok=True)
            fa, sidecar = directory/f'HLA-{gene}.fa', directory/f'HLA-{gene}.json'
            info = []
            with fa.open('w') as stream:
                for group in groups:
                    # Sequence-addressed identity is stable across panel variants.
                    path = 'PG' + group['sequence_sha256'] + '#0#HLA_' + gene
                    stream.write('>' + path + '\n' + group['sequence'] + '\n')
                    info.append(dict(path=path, **{k:v for k,v in group.items() if k != 'sequence'},
                                     n_unknown=group['sequence'].count('N')))
            sidecar.write_text(json.dumps(info, indent=2)+'\n')
            for p in (fa,sidecar):
                manifest['output_sha256'][str(p.relative_to(output))] = digest(p.read_bytes())
            manifest['genes'][panel+'/'+gene] = dict(source_paths=len(selected), unique_sequences=len(groups),
                source_cohorts=dict(collections.Counter(r['cohort'] for r in selected)))
    manifest['status'] = 'complete'
    (output/'COMPLETE.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('source','metadata','catalogue','exons','reservation','output'):
        parser.add_argument('--'+key, type=Path, required=True)
    args = parser.parse_args()
    result = build(args.source,args.metadata,args.catalogue,args.exons,args.reservation,args.output)
    print(json.dumps(result['genes'], indent=2))
