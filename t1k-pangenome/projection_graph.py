"""Retain GBZ node coordinates while verifying and normalizing exported paths."""
from pathlib import Path
from build_graph import sequences, verify_paths
from evidence_io import sha


def normalize(source, output, reference):
    if output.exists():
        raise FileExistsError(output)
    segments, paths = {}, {}
    lines = source.read_text().splitlines()
    complement = str.maketrans('ACGTRYMKSWBDHVN','TGCAYRKMSWVHDBN')
    for line in lines:
        fields = line.split('\t')
        if fields[0] == 'S':
            if fields[1] in segments or fields[2] == '*':
                raise ValueError('Missing/duplicate graph segment sequence')
            segments[fields[1]] = fields[2].upper()
        elif fields[0] == 'P':
            if fields[1] in paths:
                raise ValueError('Duplicate exported path')
            paths[fields[1]] = fields[2].split(',')
    actual = {}
    for name, walk in paths.items():
        parts = []
        for node in walk:
            if node[-1:] not in ('+','-') or node[:-1] not in segments:
                raise ValueError('Invalid exported path walk')
            sequence = segments[node[:-1]]
            parts.append(sequence if node[-1]=='+' else sequence.translate(complement)[::-1])
        actual[name] = ''.join(parts)
    aliases = verify_paths(sequences(reference),actual)
    with output.open('w') as target:
        for line in lines:
            fields = line.split('\t')
            if fields[0] == 'P':
                fields[1] = aliases[fields[1]]
            target.write('\t'.join(fields)+'\n')
    return dict(export_sha256=sha(source),reference_sha256=sha(reference),
                gfa_sha256=sha(output),verified_paths=len(aliases),path_aliases=aliases,
                driver_sha256=sha(Path(__file__)))
