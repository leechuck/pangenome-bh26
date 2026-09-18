"""Estimate insert size from concordant graph placements across competing loci.

Development calibration only. A fragment contributes once, regardless of how
many graph paths support it. No allele labels or evaluation truth are read.
"""
import argparse
import collections
from concurrent.futures import ProcessPoolExecutor
import json
import math
import os
from pathlib import Path
import statistics
import time
from evidence_io import sha
from fragment_support import collect_disk
from join_evidence import fragment_id
from build_panel import GENES


def candidate(fragment, gene):
    placements=[p for values in fragment['candidates'].values() for p in values]
    if not placements:return None
    best=max(p['alignment_score'] for p in placements)
    spans=[p['insert_size'] for p in placements if p['alignment_score']==best]
    span=statistics.median(spans) if max(spans)-min(spans)<=5 else None
    return dict(fragment=fragment_id(fragment['read_names']),gene=gene,score=best,span=span)


def select(records, margin=10):
    by_fragment=collections.defaultdict(dict)
    for row in records:
        group=by_fragment[row['fragment']]
        if row['gene'] in group:raise ValueError('Repeated fragment in one locus')
        group[row['gene']]=row
    selected=[];counts=collections.Counter()
    for group in by_fragment.values():
        ordered=sorted(group.values(),key=lambda row:row['score'],reverse=True)
        top=ordered[0]
        if len(ordered)>1 and top['score']-ordered[1]['score']<margin:
            counts['competing_locus']+=1;continue
        if top['span'] is None:
            counts['ambiguous_insert']+=1;continue
        selected.append(top);counts['selected']+=1
    return selected,dict(counts)


def estimate(spans):
    if len(spans)<200:return dict(usable=False,reason='Fewer than 200 independent fragments')
    center=statistics.median(spans)
    scale=1.4826*statistics.median(abs(value-center) for value in spans)
    if scale<=0:return dict(usable=False,reason='Degenerate insert distribution')
    core=[value for value in spans if abs(value-center)<=4*scale]
    mean,stdev=statistics.mean(core),statistics.pstdev(core)
    return dict(usable=len(core)>=200 and mean>0 and stdev>0,
                fragment_mean=mean,fragment_stdev=stdev,median=center,robust_scale=scale,
                input_fragments=len(spans),retained_fragments=len(core),trimmed_fragments=len(spans)-len(core),
                observed_min=min(spans),observed_max=max(spans))


def one_gene(arguments):
    gene,source,mapping,output=arguments
    projection=json.loads((source/'COMPLETE.json').read_text())
    support=source/'support'
    sm=json.loads((support/'COMPLETE.json').read_text())
    mm=json.loads((mapping/'COMPLETE.json').read_text())
    if any(m['status']!='complete' for m in (projection,sm,mm)):
        raise ValueError('Incomplete graph evidence')
    if (projection['mapping_manifest_sha256']!=sha(mapping/'COMPLETE.json') or
        sm['alignments_sha256']!=mm['output_sha256']['alignments.jsonl'] or
        sm['gfa_sha256']!=projection['projection_graph']['gfa_sha256'] or
        sha(support/'support.jsonl')!=sm['output_sha256']):
        raise ValueError('Graph evidence chain changed')
    counts=collections.Counter();rows=[]
    def filtered(stream):
        for line in stream:
            row=json.loads(line)
            identity=float(row.get('identity') or 0);quality=float(row.get('mapping_quality') or 0)
            if not math.isfinite(identity) or not math.isfinite(quality):raise ValueError('Nonfinite alignment quality')
            counts['alignment_records']+=1
            if row['placements'] and identity>=.98 and quality>=30:
                counts['retained_alignment_records']+=1
                yield row
    with (support/'support.jsonl').open() as stream:
        for fragment in collect_disk(filtered(stream),1000,output):
            row=candidate(fragment,gene)
            if row is not None:rows.append(row)
    return dict(gene=gene,rows=rows,counts=dict(counts),reads_sha256=mm['reads_sha256'],
                support_manifest_sha256=sha(support/'COMPLETE.json'),
                projection_manifest_sha256=sha(source/'COMPLETE.json'))


def run(source,mapping,output,threads):
    if output.exists():raise FileExistsError(output)
    if not 1<=threads<=int(os.environ.get('SLURM_CPUS_PER_TASK','0')):
        raise ValueError('Requires sufficient Slurm allocation')
    output.mkdir(parents=True)
    started=time.time()
    tasks=[(gene,source/gene,mapping/gene,output) for gene in GENES]
    with ProcessPoolExecutor(max_workers=threads) as pool:
        results=list(pool.map(one_gene,tasks))
    reads=results[0]['reads_sha256']
    if any(result['reads_sha256']!=reads for result in results):raise ValueError('Loci use different input reads')
    chosen,counts=select(row for result in results for row in result['rows'])
    fit=estimate([row['span'] for row in chosen])
    genes=collections.Counter(row['gene'] for row in chosen)
    if len(genes)<2:fit.update(usable=False,reason='Support from fewer than two loci')
    report=dict(status='complete' if fit['usable'] else 'failed',started=started,finished=time.time(),
                reads_sha256=reads,estimate=fit,selection_counts=counts,selected_by_gene=dict(genes),
                parameters=dict(minimum_identity=.98,minimum_mapq=30,locus_score_margin=10,
                                maximum_span_disagreement=5,maximum_insert=1000,minimum_fragments=200),
                sources=[{k:v for k,v in result.items() if k!='rows'} for result in results],
                driver_sha256=sha(Path(__file__)),
                scope='Insert calibration only; one fragment once; no genotype truth or ancestry prior',
                limitations='Conditional on alignable conserved graph context and inserts <=1000; selected-locus bias remains possible')
    (output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    if not fit['usable']:raise RuntimeError('Insufficient library calibration')
    (output/'COMPLETE.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('source','mapping','output'):parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--threads',type=int,default=8)
    args=parser.parse_args();run(args.source,args.mapping,args.output,args.threads)
