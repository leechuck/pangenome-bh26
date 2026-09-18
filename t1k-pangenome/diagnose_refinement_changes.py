"""Partition changed-pair score evidence; diagnostic only, never changes calls."""
import argparse
import collections
import json
import math
import os
from pathlib import Path
from build_graph import sequences
from build_panel import GENES
from evidence_io import sha
from graph_pair_refinement import assign_fragment
from run_graph_pair_refinement import candidate_metadata, native_pair


def run(root, donor, condition, output):
    if not os.environ.get('SLURM_CPUS_PER_TASK'):
        raise RuntimeError('Run under Slurm')
    folder=root/'development/graph-pair-refinement-v1'/donor/condition
    report=json.loads((folder/'COMPLETE.json').read_text())
    decisions=json.loads((folder/'decisions.json').read_text())
    if sha(folder/'decisions.json')!=report['output_sha256']['decisions.json']:
        raise ValueError('Decisions changed')
    base=root/'development/linear-v1/ipd_genome'/donor
    if sha(base/'COMPLETE.json')!=report['baseline_sha256']:
        raise ValueError('Baseline changed')
    bm=json.loads((base/'COMPLETE.json').read_text())
    if sha(base/'calls.json')!=bm['output_sha256']['calls.json']:
        raise ValueError('Baseline calls changed')
    calls=json.loads((base/'calls.json').read_text())
    calibration=root/'development/library-calibration-v1'/donor/'COMPLETE.json'
    if sha(calibration)!=report['calibration_sha256']:
        raise ValueError('Calibration changed')
    insert=json.loads(calibration.read_text())['estimate']['fragment_mean']
    panel=condition.split('/')[0]
    reference=root/'references/observed-v2'
    rm=json.loads((reference/'COMPLETE.json').read_text())
    if sha(reference/'COMPLETE.json')!=report['reference_manifest_sha256']:
        raise ValueError('Reference changed')
    paths={};spans={}
    for gene in GENES:
        fasta=reference/panel/f'HLA-{gene}.fa';metadata=fasta.with_suffix('.json')
        for file in (fasta,metadata):
            if sha(file)!=rm['output_sha256'][str(file.relative_to(reference))]:
                raise ValueError('Reference component changed')
        records=json.loads(metadata.read_text())
        pp,_=candidate_metadata(gene,records,sequences(fasta),native_pair(calls.get(gene,{})),insert)
        paths.update(pp)
        for r in records:
            spans[r['path']]=(min(a for a,b in r['exons']),max(b for a,b in r['exons'])) if r['exons'] else None
    joined=root/'development/join-v3'/donor/condition
    jm=json.loads((joined/'COMPLETE.json').read_text())
    source=joined/'evidence.jsonl'
    if sha(joined/'COMPLETE.json')!=report['joined_sha256'] or sha(source)!=jm['output_sha256']:
        raise ValueError('Joined evidence changed')
    changed={g:d for g,d in decisions.items() if d['changed']}
    totals={g:collections.defaultdict(lambda:dict(fragments=0,score_gain=0.0)) for g in changed}
    with source.open() as stream:
        for line in stream:
            item=json.loads(line);assigned=assign_fragment(item,paths)
            if assigned is None or assigned[0] not in changed:continue
            gene,values=assigned;old=native_pair(calls[gene]);new=changed[gene]['pair']
            emission=lambda pair:.01+.99*.5*sum(values.get(a,0) for a in pair)
            gain=math.log(emission(new))-math.log(emission(old))
            overlap=False
            for graph in item['graph_sources'].values():
                for path,placements in graph['candidates'].items():
                    if paths[path]['gene']!=gene or not set(paths[path]['candidates'])&set(old+new):continue
                    span=spans[path]
                    if span is None:continue
                    for p in placements:
                        if any(p[m+'_start']<span[1] and p[m+'_end']>span[0] for m in ('first','second')):
                            overlap=True
            for category in ('all', 'overlaps_gene_span' if overlap else 'outside_gene_span',
                             'both_native_alleles_supported' if all(values.get(a,0)>0 for a in old)
                             else 'some_native_allele_unsupported'):
                totals[gene][category]['fragments']+=1
                totals[gene][category]['score_gain']+=gain
    for gene in changed:
        if not math.isclose(totals[gene]['all']['score_gain'],changed[gene]['score_gain_over_native'],abs_tol=1e-6):
            raise ValueError('Diagnostic failed to reproduce recorded pair score')
    result=dict(donor=donor,condition=condition,decisions=changed,score_partitions=totals,
                scope='Development diagnostics; gene span means envelope of annotated exons, including introns but excluding terminal UTRs',
                manifest_sha256=sha(folder/'COMPLETE.json'),driver_sha256=sha(Path(__file__)))
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--donor',required=True);p.add_argument('--condition',required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.root,a.donor,a.condition,a.output)
