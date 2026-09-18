"""Supplement native T1K candidate pairs with graph-supported pairs.

Graph scores select reads only. Native T1K performs all allele competition and
four-field genotyping against the unchanged full genomic IPD reference.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import os
from pathlib import Path
import subprocess
import time
from evidence_io import sha
from join_evidence import fragment_id
from decode_t1k import decode


def fastq(path):
    opener = gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt') as stream:
        while True:
            head=stream.readline()
            if not head:return
            seq,plus,quality=(stream.readline() for _ in range(3))
            if not head.startswith('@') or not plus.startswith('+') or not quality:
                raise ValueError('Malformed FASTQ')
            if len(seq.rstrip())!=len(quality.rstrip()):raise ValueError('FASTQ length mismatch')
            yield head,seq,plus,quality


def paired(first,second):
    for a,b in itertools.zip_longest(fastq(first),fastq(second)):
        if a is None or b is None:raise ValueError('Unequal paired FASTQ lengths')
        names=[r[0][1:].split()[0] for r in (a,b)]
        yield fragment_id(names),a,b


def content_digest(a,b):
    return hashlib.sha256(''.join((a[1],a[3],b[1],b[3])).encode()).hexdigest()


def select_reads(raw,candidates,output,graph_ids=(),all_reads=False):
    native={}
    for identifier,a,b in paired(*candidates):
        if identifier in native:raise ValueError('Duplicate native fragment')
        native[identifier]=content_digest(a,b)
    wanted=set(native)|set(graph_ids)
    counts=dict(native_pairs=len(native),graph_supported_pairs=len(set(graph_ids)),
                input_pairs=0,output_pairs=0,added_pairs=0)
    seen=set()
    with (output/'r1.fq').open('w') as first,(output/'r2.fq').open('w') as second:
        for identifier,a,b in paired(*raw):
            if identifier in seen:raise ValueError('Duplicate raw fragment')
            seen.add(identifier);counts['input_pairs']+=1
            if identifier in native and content_digest(a,b)!=native[identifier]:
                raise ValueError('Native candidate sequence/quality differs from raw reads')
            if all_reads or identifier in wanted:
                first.writelines(a);second.writelines(b)
                counts['output_pairs']+=1;counts['added_pairs']+=identifier not in native
    if not wanted<=seen:raise ValueError('Candidate fragment absent from raw reads')
    return counts


def graph_candidates(folder,raw_hashes):
    report=json.loads((folder/'COMPLETE.json').read_text())
    if report['status']!='complete':raise ValueError('Incomplete graph join')
    for chain in report['graph_chains']:
        mapping=[Path(p) for p in chain if '/mapping-' in p]
        if len(mapping)!=1:raise ValueError('Missing mapping provenance')
        mm=mapping[0]/'COMPLETE.json'
        if sha(mm)!=chain[str(mapping[0])]:raise ValueError('Mapping manifest changed')
        record=json.loads(mm.read_text())
        if record['status']!='complete' or sorted(record['reads_sha256'].values())!=sorted(raw_hashes.values()):
            raise ValueError('Graph and native reads differ')
    evidence=folder/'evidence.jsonl'
    if sha(evidence)!=report['output_sha256']:raise ValueError('Graph evidence changed')
    identifiers=set()
    with evidence.open() as stream:
        for line in stream:
            row=json.loads(line)
            if any(source['candidates'] for source in row['graph_sources'].values()):
                identifiers.add(row['fragment'])
    return identifiers


def run(baseline,reads,reference,tool,output,mode,graph,threads):
    if not 1<=threads<=int(os.environ.get('SLURM_CPUS_PER_TASK','0')):
        raise ValueError('Requires Slurm allocation')
    if output.exists():raise FileExistsError(output)
    if (mode=='graph')!=(graph is not None):raise ValueError('Graph input required only for graph mode')
    original=json.loads((baseline/'COMPLETE.json').read_text())
    if original['status']!='complete':raise ValueError('Incomplete native control')
    donor=original['donor']
    raw=[reads/'r1.fq.gz',reads/'r2.fq.gz']
    hashes={str(i):sha(p) for i,p in enumerate(raw,1)}
    rm=json.loads((reference/'COMPLETE.json').read_text())
    if (hashes!=original['reads_sha256'] or rm['status']!='complete' or
        sha(reference/'COMPLETE.json')!=original['reference_manifest_sha256'] or
        sha(reference/'reference.fa')!=rm['output_sha256']['reference.fa']):
        raise ValueError('Baseline inputs changed')
    for name,digest in original['output_sha256'].items():
        if sha(baseline/name)!=digest:raise ValueError('Baseline output changed')
    build=json.loads((tool/'COMPLETE.json').read_text());exe=tool/'original/genotyper'
    if build['status']!='complete' or sha(exe)!=build['output_sha256']['original/genotyper']:
        raise ValueError('Native executable changed')
    ids=graph_candidates(graph,hashes) if graph else set()
    output.mkdir(parents=True)
    record=dict(status='running',started=time.time(),donor=donor,mode=mode,
                baseline_sha256=sha(baseline/'COMPLETE.json'),reads_sha256=hashes,
                reference_sha256=sha(reference/'reference.fa'),tool_sha256=sha(exe),
                driver_sha256=sha(Path(__file__)),graph_manifest_sha256=sha(graph/'COMPLETE.json') if graph else None,
                scope='Graph-assisted recruitment with unchanged native T1K allele inference; no genotype truth')
    def save():(output/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    save()
    try:
        candidates=[baseline/(donor+'_candidate_'+str(i)+'.fq') for i in (1,2)]
        record['selection']=select_reads(raw,candidates,output,ids,mode=='all')
        command=[str(exe),'-s','0.97','--alleleDigitUnits','4','--alleleDelimiter',':',
                 '-f',str(reference/'reference.fa'),'-1',str(output/'r1.fq'),'-2',str(output/'r2.fq'),
                 '-t',str(threads),'-o',str(output/donor)]
        record['command']=command;save()
        with (output/'run.log').open('w') as log:
            subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
        table=donor+'_genotype.tsv'
        if mode=='native' and (output/table).read_bytes()!=(baseline/table).read_bytes():
            raise ValueError('Native recruitment parity failed')
        calls=decode((output/table).read_text(),{})
        if not calls:raise ValueError('Empty native genotype calls')
        (output/'calls.json').write_text(json.dumps(calls,indent=2)+'\n')
        record.update(status='complete',finished=time.time(),output_sha256={
            name:sha(output/name) for name in (table,'calls.json','r1.fq','r2.fq')})
        save();(output/'COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
    except BaseException as error:
        record.update(status='failed',finished=time.time(),error=repr(error));save();raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('baseline','reads','reference','tool','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--mode',choices=('native','all','graph'),required=True)
    p.add_argument('--graph',type=Path)
    p.add_argument('--threads',type=int,default=4)
    a=p.parse_args();run(a.baseline,a.reads,a.reference,a.tool,a.output,a.mode,a.graph,a.threads)
