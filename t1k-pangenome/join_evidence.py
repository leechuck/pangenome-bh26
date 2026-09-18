"""Join native and graph fragment evidence without treating them as independent reads."""
import argparse
import collections
from contextlib import closing
import json
from pathlib import Path
import sqlite3
from build_graph import sha
from t1k_evidence import parse


def fragment_id(names):
    if len(names)!=2 or any(not name or any(c.isspace() for c in name) for name in names):
        raise ValueError('Expected two explicit FASTQ identifiers')
    normalized = [name[:-2] if name.endswith(('/1','/2')) else name for name in names]
    if normalized[0]!=normalized[1] or not normalized[0]:
        raise ValueError('Graph mates do not share a native T1K fragment identity')
    return normalized[0]


def merge(connection, native, graphs):
    connection.execute('CREATE TABLE evidence (fragment TEXT, source TEXT, record TEXT, PRIMARY KEY(fragment,source))')
    def insert(identifier,source,record):
        try:
            connection.execute('INSERT INTO evidence VALUES (?,?,?)',
                               (identifier,source,json.dumps(record,separators=(',',':'))))
        except sqlite3.IntegrityError as error:
            raise ValueError('Repeated fragment in one evidence source: '+identifier) from error
    for record in native:
        insert(record['fragment'],'native',record['candidates'])
    for number,records in enumerate(graphs):
        for record in records:
            insert(fragment_id(record['read_names']),'graph:'+str(number),record)
    connection.commit()
    current,record = None,None
    for identifier,source,payload in connection.execute('SELECT fragment,source,record FROM evidence ORDER BY fragment,source'):
        if identifier!=current:
            if record is not None:
                yield record
            current = identifier
            record = dict(fragment=identifier,native_candidates={},graph_sources={})
        if source=='native':
            record['native_candidates'] = json.loads(payload)
        else:
            record['graph_sources'][source] = json.loads(payload)
    if record is not None:
        yield record


def build(native, native_file, fragment_dirs, support_dirs, mapping_dirs, output):
    if output.exists():
        raise FileExistsError(output)
    if not fragment_dirs or len(fragment_dirs)!=len(support_dirs) or len(fragment_dirs)!=len(mapping_dirs):
        raise ValueError('One complete provenance chain per graph required')
    native_record = json.loads((native/'COMPLETE.json').read_text())
    native_key = str(native_file.resolve().relative_to(native.resolve()))
    if native_record['status']!='complete' or sha(native_file)!=native_record['output_sha256'][native_key]:
        raise ValueError('Native evidence changed or incomplete')
    native_reads = sorted(digest for name,digest in native_record['input_sha256'].items()
                          if name.endswith(('.fq','.fq.gz','.fastq','.fastq.gz')))
    if len(native_reads)!=2:
        raise ValueError('Expected native input paired-read provenance')
    chains = []
    for fragments,support,mapping in zip(fragment_dirs,support_dirs,mapping_dirs):
        fm,sm,mm = [json.loads((folder/'COMPLETE.json').read_text()) for folder in (fragments,support,mapping)]
        if any(record['status']!='complete' for record in (fm,sm,mm)):
            raise ValueError('Incomplete graph evidence')
        if sorted(mm['reads_sha256'].values())!=native_reads:
            raise ValueError('Native and graph evidence came from different reads')
        if fm['source_manifest_sha256']!=sha(support/'COMPLETE.json') or fm['output_sha256']!=sha(fragments/'fragments.jsonl'):
            raise ValueError('Fragment evidence chain changed')
        if sm['alignments_sha256']!=sha(mapping/'alignments.jsonl') or sm['output_sha256']!=sha(support/'support.jsonl'):
            raise ValueError('Graph alignment evidence chain changed')
        chains.append({str(folder):sha(folder/'COMPLETE.json') for folder in (fragments,support,mapping)})
    output.mkdir(parents=True)
    counts = collections.Counter()
    handles = []
    try:
        native_stream = native_file.open();handles.append(native_stream)
        streams = []
        for folder in fragment_dirs:
            stream = (folder/'fragments.jsonl').open();handles.append(stream)
            streams.append(map(json.loads,stream))
        with closing(sqlite3.connect(output/'join.sqlite')) as connection, (output/'evidence.jsonl').open('w') as target:
            for record in merge(connection,parse(native_stream),streams):
                has_native = bool(record['native_candidates'])
                has_graph = any(source['candidates'] for source in record['graph_sources'].values())
                category = 'both' if has_native and has_graph else 'native_only' if has_native else 'graph_only' if has_graph else 'neither'
                counts[category] += 1
                target.write(json.dumps(record,separators=(',',':'))+'\n')
        report = dict(status='complete',native_manifest_sha256=sha(native/'COMPLETE.json'),
                      graph_chains=chains,driver_sha256=sha(Path(__file__)),
                      fragments=sum(counts.values()),support_categories=dict(counts),
                      output_sha256=sha(output/'evidence.jsonl'),
                      scope='One row per fragment; native weights and graph scores remain distinct, uncalibrated evidence')
        (output/'COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2))
    finally:
        for handle in handles:
            handle.close()


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for name in ('native','native-assignments','output'):
        p.add_argument('--'+name,type=Path,required=True)
    for name in ('graph-fragments','graph-support','graph-mapping'):
        p.add_argument('--'+name,type=Path,nargs='+',required=True)
    a = p.parse_args()
    build(a.native,a.native_assignments,a.graph_fragments,a.graph_support,a.graph_mapping,a.output)
