"""Post-unblinding diagnostics of fixed graph decisions; never emit new calls."""
import argparse
import collections
import json
import math
import os
from pathlib import Path

from build_graph import sequences
from build_panel import GENES
from evidence_io import sha
from graph_pair_refinement import assign_fragment, internal_fragment
from run_graph_pair_refinement import candidate_metadata, native_pair


def pair_log(values, pair, noise=.01):
    # Reproduce the frozen implementation exactly, including noise placement.
    return math.log(noise + .5*(1-noise)*sum(values.get(a, 0) for a in pair))


def run(root, donor, panel, output):
    if not os.environ.get('SLURM_CPUS_PER_TASK'):
        raise RuntimeError('Run evidence processing under Slurm')
    if output.exists():
        raise FileExistsError(output)
    base=root/'validation-runs/genomic-ipd-v1'/donor
    graph=root/'validation-runs/graph-v1'
    folder=graph/'calls'/donor/panel/'unsampled'
    record=json.loads((folder/'COMPLETE.json').read_text())
    assert record['status']=='complete'
    assert record['model_sha256']==sha(Path(__file__).with_name('graph_pair_refinement.py'))
    assert record['driver_sha256']==sha(Path(__file__).with_name('run_graph_pair_refinement.py'))
    assert record['parameters']['internal_pairs'] and record['parameters']['pair_specific']
    assert record['baseline_sha256']==sha(base/'COMPLETE.json')
    baseline=json.loads((base/'COMPLETE.json').read_text())
    assert sha(base/'calls.json')==baseline['output_sha256']['calls.json']
    calls=json.loads((base/'calls.json').read_text())
    assert sha(folder/'decisions.json')==record['output_sha256']['decisions.json']
    decisions=json.loads((folder/'decisions.json').read_text())
    changed={g:d for g,d in decisions.items() if d['changed']}
    calibration=graph/'library-calibration-v1'/donor/'COMPLETE.json'
    assert sha(calibration)==record['calibration_sha256']
    insert=json.loads(calibration.read_text())['estimate']['fragment_mean']
    reference=root/'references/observed-v2'
    assert sha(reference/'COMPLETE.json')==record['reference_manifest_sha256']
    rm=json.loads((reference/'COMPLETE.json').read_text())
    paths={}
    for gene in GENES:
        fasta=reference/panel/f'HLA-{gene}.fa';metadata=fasta.with_suffix('.json')
        for p in (fasta,metadata):
            assert sha(p)==rm['output_sha256'][str(p.relative_to(reference))]
        pp,_=candidate_metadata(gene,json.loads(metadata.read_text()),sequences(fasta),
                                native_pair(calls.get(gene,{})),insert)
        paths.update(pp)
    unit_paths={p:dict(m,effective_length=1.) for p,m in paths.items()}
    joined=graph/'join-v3'/donor/panel/'unsampled'
    assert sha(joined/'COMPLETE.json')==record['joined_sha256']
    jm=json.loads((joined/'COMPLETE.json').read_text())
    source=joined/'evidence.jsonl'
    assert sha(source)==jm['output_sha256']
    totals={g:collections.defaultdict(lambda:dict(fragments=0,frozen_gain=0.,equal_length_gain=0.))
            for g in changed}
    with source.open() as stream:
        for line in stream:
            item=json.loads(line)
            assigned=assign_fragment(item,paths)
            if assigned is None or assigned[0] not in changed:
                continue
            gene,values=assigned
            if not internal_fragment(item,paths,gene):
                continue
            equal=assign_fragment(item,unit_paths)
            assert equal is not None and equal[0]==gene
            old=native_pair(calls[gene]);new=changed[gene]['pair']
            gain=pair_log(values,new)-pair_log(values,old)
            equal_gain=pair_log(equal[1],new)-pair_log(equal[1],old)
            categories=['all', 'length_only_pair_difference' if abs(equal_gain)<=1e-8
                        else 'alignment_pair_difference',
                        'native_allele_absent' if any(values.get(a,0)==0 for a in old)
                        else 'both_native_alleles_represented']
            for category in categories:
                row=totals[gene][category]
                row['fragments']+=1;row['frozen_gain']+=gain;row['equal_length_gain']+=equal_gain
    for gene,d in changed.items():
        assert totals[gene]['all']['fragments']==d['fragments']
        assert math.isclose(totals[gene]['all']['frozen_gain'],d['score_gain_over_native'],abs_tol=1e-6)
    result=dict(donor=donor,panel=panel,decisions=changed,score_partitions=totals,
                scope='Exploratory post-unblinding fixed-pair diagnostics; no new calls or validation claim',
                manifest_sha256=sha(folder/'COMPLETE.json'),joined_sha256=sha(joined/'COMPLETE.json'),
                driver_sha256=sha(Path(__file__)))
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--donor',required=True)
    p.add_argument('--panel',choices=('hprc','hprc_asian'),required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.root,a.donor,a.panel,a.output)
