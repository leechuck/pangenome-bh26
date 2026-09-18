"""Project graph alignments onto compatible observed paths, retaining ambiguity.

This is an evidence adapter, not a genotype caller. All compatible placements
are retained; unsupported graph walks are explicit rather than forced to a path.
"""
import argparse
import collections
import json
from pathlib import Path
from evidence_io import sha


class PathIndex:
    def __init__(self, lengths, paths):
        self.lengths = lengths
        self.walks = {}
        self.starts = collections.defaultdict(list)
        for name, walk in paths.items():
            for reverse in (False, True):
                oriented = tuple(-n for n in reversed(walk)) if reverse else tuple(walk)
                prefix = [0]
                for node in oriented:
                    prefix.append(prefix[-1]+lengths[abs(node)])
                self.walks[name, reverse] = (oriented, prefix)
                for i, node in enumerate(oriented):
                    self.starts[node].append((name, reverse, i))

    @classmethod
    def from_gfa(cls, path):
        lengths, paths = {}, {}
        with path.open() as stream:
            for line in stream:
                fields = line.rstrip('\n').split('\t')
                if fields[0] == 'S':
                    if fields[2] == '*':
                        raise ValueError('GFA requires explicit segment sequences')
                    node = int(fields[1])
                    if node <= 0 or node in lengths:
                        raise ValueError('Invalid/duplicate GFA node')
                    lengths[node] = len(fields[2])
                elif fields[0] == 'P':
                    if fields[1] in paths:
                        raise ValueError('Duplicate GFA path')
                    paths[fields[1]] = tuple(int(s[:-1])*(1 if s[-1]=='+' else -1)
                                             for s in fields[2].split(','))
        if not paths:
            raise ValueError('No named P paths in GFA')
        return cls(lengths, paths)

    def placements(self, alignment):
        mappings = alignment.get('path', {}).get('mapping', [])
        if not mappings:
            return [], 'unmapped'
        walk, bounds = [], []
        for i, mapping in enumerate(mappings):
            position = mapping['position']
            node = int(position['node_id'])
            offset = int(position.get('offset', 0))
            length = sum(int(e.get('from_length',0)) for e in mapping.get('edit',[]))
            if node not in self.lengths:
                raise ValueError('Alignment node absent from indexed graph')
            if not 0 <= offset < offset+length <= self.lengths[node]:
                return [], 'unsupported_mapping_bounds'
            if i and offset != 0 or i < len(mappings)-1 and offset+length != self.lengths[node]:
                return [], 'noncontiguous_mapping_bounds'
            walk.append(-node if position.get('is_reverse',False) else node)
            bounds.append((offset,length))
        result = []
        for name, reverse, start in self.starts.get(walk[0],[]):
            oriented, prefix = self.walks[name,reverse]
            if oriented[start:start+len(walk)] != tuple(walk):
                continue
            begin = prefix[start]+bounds[0][0]
            end = prefix[start+len(walk)-1]+sum(bounds[-1])
            if reverse:
                begin,end = prefix[-1]-end,prefix[-1]-begin
            result.append(dict(path=name, reverse=reverse, start=begin, end=end))
        return result, 'supported' if result else 'no_observed_path'


def project(gfa, alignments, output):
    if output.exists():
        raise FileExistsError(output)
    index = PathIndex.from_gfa(gfa)
    output.mkdir(parents=True)
    counts = collections.Counter()
    with alignments.open() as source, (output/'support.jsonl').open('w') as target:
        for line in source:
            alignment = json.loads(line)
            placements, status = index.placements(alignment)
            counts[status] += 1
            target.write(json.dumps(dict(name=alignment['name'], status=status,
                read_length=len(alignment.get('sequence','')),
                fragment_prev=alignment.get('fragment_prev'), fragment_next=alignment.get('fragment_next'),
                secondary=alignment.get('is_secondary',False),
                score=alignment.get('score'), identity=alignment.get('identity'),
                mapping_quality=alignment.get('mapping_quality'), placements=placements))+'\n')
    report = dict(status='complete', gfa_sha256=sha(gfa), alignments_sha256=sha(alignments),
                  driver_sha256=sha(Path(__file__)), outcomes=dict(counts),
                  output_sha256=sha(output/'support.jsonl'),
                  scope='Path-compatible placements only; no genotype inference or independent accuracy claim')
    (output/'COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--gfa', type=Path, required=True)
    p.add_argument('--alignments', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    project(a.gfa,a.alignments,a.output)
