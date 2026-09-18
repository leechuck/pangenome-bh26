"""Add observed genomic contexts to a complete T1K IPD reference.

This is the linear-reference control, not graph alignment. PG identifiers are
internal sequence candidates, never official HLA names. Retain the alias map.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
from build_graph import sequences, sha
from build_panel import GENES


def header(alias, intervals, length):
    previous = -1
    positions = []
    if not intervals:
        raise ValueError('Missing coding annotation')
    for start, end in intervals:
        if not 0 <= start < end <= length or start < previous:
            raise ValueError('Invalid coding interval')
        positions.extend((start, end-1))
        previous = end
    return '>'+alias+' '+str(len(intervals))+' '+' '.join(map(str, positions))+'\n'


def gene_namespace(gene, names):
    matches = {name.split('*')[0] for name in names
               if '*' in name and name.split('*')[0].removeprefix('HLA-') == gene}
    if len(matches) != 1:
        raise ValueError('Missing or ambiguous IPD gene namespace: '+gene)
    return matches.pop()


def resolve_alias(alias, aliases, fields):
    """Return numeric-field alternatives and an explicit unresolved flag.

    Full CDS identity supports two fields only. Unknown candidates must remain
    unresolved even if another candidate in a genotype has a known label.
    """
    if fields not in (2, 4):
        raise ValueError('Supported resolutions: 2 or 4 fields')
    record = aliases[alias]
    labels = record['genomic_labels'] if fields == 4 else record['cds_labels']
    calls = set()
    unresolved = not labels
    for label in labels:
        match = re.fullmatch(r'([A-Z0-9]+)\*([0-9:]+)([A-Z]*)', label)
        if not match or len(match[2].split(':')) < fields:
            unresolved = True
        else:
            calls.add(match[1]+'*'+':'.join(match[2].split(':')[:fields]))
    return dict(alleles=sorted(calls), unresolved=unresolved,
                original_labels=labels, resolution=fields)


def build(reference, panel, ipd, output):
    if output.exists():
        raise FileExistsError(output)
    manifest = json.loads((reference/'COMPLETE.json').read_text())
    if manifest['status'] != 'complete':
        raise ValueError('Incomplete observed reference')
    original = ipd.read_bytes()
    if not original.startswith(b'>') or not original.endswith(b'\n'):
        raise ValueError('IPD reference must be FASTA ending in a newline')
    # Existing IPD records (including paralogs) stay byte-for-byte unchanged.
    ipd_sequences = sequences(ipd)
    existing_sequences = set(ipd_sequences.values())
    aliases, skipped, additions = {}, [], []
    for gene in GENES:
        namespace = gene_namespace(gene, ipd_sequences)
        fa, sidecar = f'{panel}/HLA-{gene}.fa', f'{panel}/HLA-{gene}.json'
        for name in (fa, sidecar):
            if sha(reference/name) != manifest['output_sha256'][name]:
                raise ValueError('Observed reference changed: '+name)
        seqs = sequences(reference/fa)
        records = json.loads((reference/sidecar).read_text())
        if len(records) != len(seqs) or {r['path'] for r in records} != set(seqs):
            raise ValueError('Sidecar/path mismatch')
        for record in records:
            seq = seqs[record['path']]
            digest = hashlib.sha256(seq.encode()).hexdigest()
            if digest != record['sequence_sha256']:
                raise ValueError('Sequence/annotation mismatch')
            reason = ('missing_coding_annotation' if not record['exons'] else
                      'unknown_bases' if 'N' in seq else
                      'sequence_already_in_ipd' if seq in existing_sequences else None)
            if reason:
                skipped.append(dict(path=record['path'], reason=reason));continue
            alias = namespace+'*PG'+digest[:12]
            if alias in aliases or alias in ipd_sequences:
                raise ValueError('Internal candidate identifier collision: '+alias)
            aliases[alias] = dict(record, gene=gene, internal_identifier=True)
            additions.append(header(alias, record['exons'], len(seq))+seq+'\n')
    output.mkdir(parents=True)
    (output/'reference.fa').write_bytes(original+''.join(additions).encode())
    (output/'aliases.json').write_text(json.dumps(aliases, indent=2)+'\n')
    report = dict(status='complete', panel=panel, driver_sha256=sha(Path(__file__)),
                  source_manifest_sha256=sha(reference/'COMPLETE.json'), ipd_sha256=sha(ipd),
                  ipd_records=len(ipd_sequences), added_contexts=len(aliases), skipped=skipped,
                  scope='Linear sequence-information control; no graph alignment or ancestry prior',
                  limitation='T1K treats contexts as separate sequence candidates. Its weighting '
                             'and tie-breaking may depend on reference multiplicity; this is not '
                             'a calibrated allele-level graph likelihood.',
                  output_sha256={p:sha(output/p) for p in ('reference.fa','aliases.json')})
    (output/'COMPLETE.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='skipped'}, indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path, required=True)
    p.add_argument('--panel', choices=('hprc','hprc_asian'), required=True)
    p.add_argument('--ipd', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    build(a.reference, a.panel, a.ipd, a.output)
