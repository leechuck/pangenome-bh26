"""Unblind and score a complete, provenance-verified frozen validation snapshot."""
import argparse
import csv
import datetime
import json
from pathlib import Path
from evidence_io import sha
from validation_handoff import HERE,METHOD_DIRS,verify_grid
from score_linear_controls import decode,truth_a,Nomenclature,write_tsv
from score_graph_pair_refinement import load as load_graph
from evaluate_reserved import score


def run(snapshot,output):
    if output.exists():raise FileExistsError('Preserve the existing evaluation; do not overwrite it')
    ep_path=HERE/'validation/EVALUATION_FREEZE.json';ep=json.loads(ep_path.read_text())
    hp_path=HERE/'validation/HANDOFF_FREEZE.json';hp=json.loads(hp_path.read_text())
    if hp['evaluation_freeze_sha256']!=sha(ep_path):raise ValueError('Handoff/evaluation freeze mismatch')
    for filename,digest in {**ep['code_sha256'],**ep['input_sha256'],**hp['code_sha256']}.items():
        if sha(HERE.parent/filename)!=digest:raise ValueError('Frozen evaluation asset changed')
    frozen=HERE/'validation/FROZEN_VALIDATION.json';plan=json.loads(frozen.read_text())
    if sha(frozen)!=ep['inference_freeze_sha256']:raise ValueError('Inference freeze mismatch')
    ready_file=snapshot/'EXECUTION_READY.json';ready=json.loads(ready_file.read_text())
    if (not ready['ready'] or ready['inference_freeze_sha256']!=sha(frozen) or
        ready['evaluation_freeze_sha256']!=sha(ep_path)):
        raise ValueError('Snapshot is not ready for the frozen evaluation')
    cohort_path=HERE/'validation/cohort.tsv'
    if sha(cohort_path)!=plan['metadata_sha256']['cohort.tsv']:raise ValueError('Cohort changed')
    cohort=[dict(r,stratum=r['superpopulation']) for r in csv.DictReader(cohort_path.open(),delimiter='\t')]
    indexed,complete=verify_grid(ready['runs'],[r['donor'] for r in cohort])
    if not complete:raise ValueError('Incomplete snapshot; failure audit required before scoring')
    # Validate every result before reading any genotype table or evaluation truth.
    for item in indexed.values():
        relative=METHOD_DIRS[item['method']].format(donor=item['donor'])
        if item['relative']!=relative:raise ValueError('Unexpected snapshot path')
        folder=snapshot/relative
        if sha(folder/'COMPLETE.json')!=item['manifest_sha256']:raise ValueError('Completion metadata changed')
        r=json.loads((folder/'COMPLETE.json').read_text())
        if json.loads((folder/'manifest.json').read_text())!=r:raise ValueError('Manifest/completion mismatch')
        for name,digest in r['output_sha256'].items():
            if Path(name).name!=name or name in ('.','..'):raise ValueError('Unsafe output name')
            if sha(folder/name)!=digest:raise ValueError('Prediction output changed after handoff')
    unblinding=HERE/'validation/UNBLINDED.json'
    marker=dict(started=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                execution_ready_sha256=sha(ready_file),evaluation_freeze_sha256=sha(ep_path),
                handoff_freeze_sha256=sha(hp_path),scope='Reserved predictions and truth first evaluated after complete verified handoff')
    if unblinding.exists():
        prior=json.loads(unblinding.read_text())
        for key in ('execution_ready_sha256','evaluation_freeze_sha256','handoff_freeze_sha256'):
            if prior[key]!=marker[key]:raise ValueError('Previously unblinded under different inputs')
    else:unblinding.write_text(json.dumps(marker,indent=2)+'\n')
    results={}
    for (donor,method),item in indexed.items():
        folder=snapshot/item['relative']
        if method=='T1K':calls=decode((folder/(donor+'_genotype.tsv')).read_text(),{})
        elif method=='ipd_genome':
            calls=json.loads((folder/'calls.json').read_text())
            if decode((folder/(donor+'_genotype.tsv')).read_text(),{})!=calls:
                raise ValueError('Genomic table/decoded calls differ')
        else:
            state,calls,_=load_graph(folder,snapshot/METHOD_DIRS['ipd_genome'].format(donor=donor))
            if state!='complete':raise ValueError('Graph call incomplete')
        results[donor,method]=dict(state='complete',calls=calls)
    truth=truth_a(HERE.parent/'hla-analysis/results/sequence_catalogue.tsv',{r['donor'] for r in cohort})
    rows,report=score(cohort,results,truth,Nomenclature(),ep['bootstrap_replicates'],ep['bootstrap_seed'])
    report.update(inference_freeze_sha256=sha(frozen),evaluation_freeze_sha256=sha(ep_path),
                  handoff_freeze_sha256=sha(hp_path),execution_ready_sha256=sha(ready_file))
    output.mkdir(parents=True)
    write_tsv(output/'gene_scores.tsv',rows);write_tsv(output/'summary.tsv',report['summary'])
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(summary=report['summary'],primary_checks=report['contrasts']['T1K']['protocol_checks']),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.snapshot,a.output)
