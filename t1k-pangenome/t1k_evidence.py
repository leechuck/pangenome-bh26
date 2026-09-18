"""Parse instrumented native T1K assignments without dropping competing loci."""
import math
import argparse
import collections
import json
from pathlib import Path
from build_graph import sha


def parse(lines):
    """Yield each native fragment's assignments, preserving uncalibrated weights.

    Export is contiguous by fragment. Repeated fragment IDs in separated blocks
    are rejected because silently combining them could double-count evidence.
    """
    seen, current, candidates = set(), None, {}
    for line in lines:
        values = line.rstrip('\n').split('\t')
        if len(values) != 7:
            raise ValueError('Expected seven instrumented assignment columns')
        fragment, allele = values[:2]
        if not fragment or '*' not in allele:
            raise ValueError('Missing fragment or allele identity')
        start,end = map(int,values[2:4])
        weight,quality,adjusted = map(float,values[4:])
        if start < 0 or end < start or any(not math.isfinite(x) or x < 0 for x in (weight,quality,adjusted)):
            raise ValueError('Invalid native assignment coordinates or weights')
        if fragment != current:
            if fragment in seen:
                raise ValueError('Noncontiguous repeated fragment ID')
            if current is not None:
                yield dict(fragment=current,candidates=candidates)
            current,candidates = fragment,{}
            seen.add(fragment)
        assignment = dict(start=start,end=end,weight=weight,quality=quality,adjusted_weight=adjusted)
        bucket = candidates.setdefault(allele,[])
        if assignment not in bucket:
            bucket.append(assignment)
    if current is not None:
        yield dict(fragment=current,candidates=candidates)


def audit(source, output):
    if output.exists():
        raise FileExistsError(output)
    fragments = assignments = cross_gene = 0
    genes = collections.Counter()
    with source.open() as stream:
        for record in parse(stream):
            fragments += 1
            assignments += sum(map(len,record['candidates'].values()))
            competing = {name.split('*')[0] for name in record['candidates']}
            cross_gene += len(competing)>1
            genes.update(competing)
    result = dict(status='complete',source_sha256=sha(source),parser_sha256=sha(Path(__file__)),
                  fragments=fragments,assignments=assignments,cross_gene_fragments=cross_gene,
                  fragments_by_gene=dict(genes),
                  scope='Native assignments before coalescing; weights are not posterior probabilities')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a = p.parse_args()
    audit(a.source,a.output)
