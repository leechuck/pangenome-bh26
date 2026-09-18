"""Translate internal T1K context candidates without inventing allele resolution."""
import argparse
import json
from pathlib import Path
from build_graph import sha
from t1k_reference import resolve_alias


def decode_token(token, aliases, fields):
    alternatives, unresolved = set(), False
    for candidate in token.split(','):
        if candidate not in aliases:
            candidate = candidate.removeprefix('HLA-')
        source = aliases if candidate in aliases else {
            candidate:dict(genomic_labels=[candidate], cds_labels=[candidate])}
        call = resolve_alias(candidate, source, fields)
        alternatives.update(call['alleles'])
        unresolved |= call['unresolved']
    return dict(raw=token, alleles=sorted(alternatives), unresolved=unresolved)


def decode(table, aliases):
    calls = {}
    for line in table.splitlines():
        if not line or line.startswith('#'):
            continue
        values = line.split('\t')
        gene = values[0].removeprefix('HLA-')
        if gene in calls:
            raise ValueError('Duplicate gene in T1K output: '+gene)
        copies = int(values[1])
        tokens = []
        if copies in (1,2):
            tokens.append(values[2] if float(values[4]) > 0 else '.')
            tokens.append(tokens[0] if copies == 1 else values[5] if float(values[7]) > 0 else '.')
        else:
            tokens = ['.','.']
        resolutions = {}
        for fields in (2,4):
            pair = [decode_token(t, aliases, fields) for t in tokens]
            for call in pair:
                if any(not allele.startswith(gene+'*') for allele in call['alleles']):
                    raise ValueError('Candidate assigned to wrong gene')
            resolutions[str(fields)] = pair
        calls[gene] = dict(copy_count=copies, resolutions=resolutions)
    return calls


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--genotypes', type=Path, required=True)
    p.add_argument('--aliases', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    aliases = json.loads(a.aliases.read_text()) if a.aliases else {}
    result = dict(genotype_sha256=sha(a.genotypes), aliases_sha256=sha(a.aliases) if a.aliases else None,
                  decoder_sha256=sha(Path(__file__)),
                  calls=decode(a.genotypes.read_text(), aliases))
    a.output.write_text(json.dumps(result, indent=2)+'\n')
