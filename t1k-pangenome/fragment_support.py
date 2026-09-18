"""Join graph placements into paired-fragment evidence without double counting."""
import argparse
import collections
import json
from pathlib import Path
from build_graph import sha


def pair_key(record):
    previous = (record.get('fragment_prev') or {}).get('name')
    following = (record.get('fragment_next') or {}).get('name')
    if bool(previous) == bool(following):
        raise ValueError('Expected exactly one explicit mate link for '+record['name'])
    return tuple(sorted((record['name'],previous or following))), 2 if previous else 1


def pair_placements(first, second, maximum_insert):
    candidates = collections.defaultdict(dict)
    for left in first:
        for right in second:
            by_path = collections.defaultdict(list)
            for placement in right['placements']:
                by_path[placement['path']].append(placement)
            for a in left['placements']:
                for b in by_path.get(a['path'],[]):
                    if a['reverse'] == b['reverse']:
                        continue
                    forward, reverse = (b,a) if a['reverse'] else (a,b)
                    span = reverse['end']-forward['start']
                    if not (forward['start'] <= reverse['start'] and
                            forward['end'] <= reverse['end'] and 0 < span <= maximum_insert):
                        continue
                    placement = dict(first_start=a['start'], first_end=a['end'],
                                     second_start=b['start'], second_end=b['end'],
                                     first_reverse=a['reverse'], insert_size=span,
                                     alignment_score=left['score']+right['score'])
                    # Duplicate alignments/placements do not create extra evidence.
                    key = (a['start'],a['end'],b['start'],b['end'],a['reverse'])
                    previous = candidates[a['path']].get(key)
                    if previous is None or placement['alignment_score'] > previous['alignment_score']:
                        candidates[a['path']][key] = placement
    return {path:list(placements.values()) for path,placements in sorted(candidates.items())}


def collect(records, maximum_insert):
    fragments = collections.defaultdict(lambda:collections.defaultdict(list))
    for record in records:
        key, slot = pair_key(record)
        fragments[key][slot].append(record)
    for names, mates in sorted(fragments.items()):
        candidates = pair_placements(mates[1],mates[2],maximum_insert) if mates[1] and mates[2] else {}
        status = 'supported' if candidates else 'missing_mate' if not mates[1] or not mates[2] else 'no_concordant_path'
        yield dict(read_names=names, status=status, candidates=candidates,
                   alignment_records={str(i):len(mates[i]) for i in (1,2)})


def build(source, output, maximum_insert):
    if output.exists():
        raise FileExistsError(output)
    if maximum_insert < 1:
        raise ValueError('Maximum insert must be positive')
    manifest = json.loads((source/'COMPLETE.json').read_text())
    support = source/'support.jsonl'
    if manifest['status']!='complete' or sha(support)!=manifest['output_sha256']:
        raise ValueError('Read-level evidence changed')
    output.mkdir(parents=True)
    counts = collections.Counter()
    with support.open() as stream, (output/'fragments.jsonl').open('w') as target:
        for fragment in collect((json.loads(line) for line in stream), maximum_insert):
            counts[fragment['status']] += 1
            target.write(json.dumps(fragment)+'\n')
    report = dict(status='complete', source_manifest_sha256=sha(source/'COMPLETE.json'),
                  driver_sha256=sha(Path(__file__)), maximum_insert=maximum_insert,
                  outcomes=dict(counts), output_sha256=sha(output/'fragments.jsonl'),
                  scope='Concordant paired placements; no calibrated likelihoods or genotype calls')
    (output/'COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--maximum-insert', type=int, default=1000)
    a = p.parse_args()
    build(a.source,a.output,a.maximum_insert)
