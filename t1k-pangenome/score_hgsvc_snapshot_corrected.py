"""Score all reserved HGSVC donors after complete, frozen provenance review."""
import argparse
import csv
import datetime
import gzip
import json
from pathlib import Path
from evidence_io import sha
from hgsvc_handoff import HERE,OUT,LOCATIONS,check_freezes,verify_reports
from hgsvc_truth import GENES,fasta,reference_index,derive
from score_linear_controls import decode,prediction,write_tsv
from anchor_t1k_coarse import anchor
from evaluate_hgsvc_corrected import score,METHODS


def run(snapshot,output):
    if output.exists():raise FileExistsError('Preserve existing evaluation')
    plan=check_freezes()
    correction=OUT/'SCORING_CORRECTION.json'
    amendment=json.loads(correction.read_text())
    if amendment['original_evaluation_freeze_sha256']!=sha(OUT/'EVALUATION_FREEZE.json'):
        raise ValueError('Correction refers to different frozen evaluation')
    for name,digest in amendment['corrected_code_sha256'].items():
        if sha(HERE/name)!=digest:raise ValueError('Corrected evaluation code changed')
    ep=json.loads((OUT/'EVALUATION_FREEZE.json').read_text())
    ready_file=snapshot/'EXECUTION_READY.json'
    ready=json.loads(ready_file.read_text())
    for key,name in [('inference_freeze_sha256','FROZEN_VALIDATION.json'),('evaluation_freeze_sha256','EVALUATION_FREEZE.json'),('handoff_freeze_sha256','HANDOFF_FREEZE.json')]:
        if ready[key]!=sha(OUT/name):raise ValueError('Snapshot freeze mismatch: '+key)
    if not ready['ready']:raise ValueError('Snapshot not ready')
    cohort=[dict(r,stratum=r['superpopulation']) for r in csv.DictReader((OUT/'cohort.tsv').open(),delimiter='\t')]
    prepared={r['donor']:r['reads_sha256'] for r in json.loads((OUT/'READ_PREPARATION.json').read_text())['donors']}
    bp=OUT/'validation-baseline-plan.json';gp=OUT/'validation-genomic-plan.json'
    if not verify_reports(ready['runs'],[r['donor'] for r in cohort],plan,prepared,json.loads(bp.read_text()),json.loads(gp.read_text()),sha(bp),sha(gp)):
        raise ValueError('Incomplete snapshot requires failure audit')
    # Validate every manifest and output digest before loading any calls or truth.
    for item in ready['runs']:
        relative=LOCATIONS[item['method']].format(donor=item['donor'])
        if item['relative']!=relative:raise ValueError('Unexpected output path')
        folder=snapshot/relative
        if sha(folder/'COMPLETE.json')!=item['manifest_sha256']:raise ValueError('Manifest changed')
        record=json.loads((folder/'COMPLETE.json').read_text())
        if record!=item['record']:raise ValueError('Snapshot record mismatch')
        if (folder/'manifest.json').exists() and json.loads((folder/'manifest.json').read_text())!=record:
            raise ValueError('Manifest/completion mismatch')
        for name,digest in record['output_sha256'].items():
            if Path(name).name!=name or name in ('.','..'):raise ValueError('Unsafe output filename')
            if sha(folder/name)!=digest:raise ValueError('Output changed')
    marker=dict(started=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                execution_ready_sha256=sha(ready_file),evaluation_freeze_sha256=sha(OUT/'EVALUATION_FREEZE.json'),
                handoff_freeze_sha256=sha(OUT/'HANDOFF_FREEZE.json'),scope='First evaluation of reserved predictions and genomic truth')
    unblinded=OUT/'UNBLINDED.json'
    if unblinded.exists():
        previous=json.loads(unblinded.read_text())
        if any(previous[k]!=marker[k] for k in ('execution_ready_sha256','evaluation_freeze_sha256','handoff_freeze_sha256')):
            raise ValueError('Previous evaluation used different inputs')
    else:unblinded.write_text(json.dumps(marker,indent=2)+'\n')
    results={};anchor_counts={}
    for row in cohort:
        donor=row['donor']
        original=snapshot/LOCATIONS['T1K'].format(donor=donor)
        native=decode((original/(donor+'_genotype.tsv')).read_text(),{})
        results[donor,'T1K']=dict(state='complete',calls=native)
        for method in METHODS[1:]:
            folder=snapshot/LOCATIONS[method].format(donor=donor)
            calls=json.loads((folder/'calls.json').read_text())
            raw='raw_genomic' if method=='ipd_genome' else 'raw_'+method.removeprefix('graph_')
            candidate=json.loads((snapshot/LOCATIONS[raw].format(donor=donor)/'calls.json').read_text())
            recomputed,decisions=anchor(native,candidate)
            if calls!=recomputed or decisions!=json.loads((folder/'decisions.json').read_text()):
                raise ValueError('Anchor composition differs from frozen algorithm')
            for gene in GENES:
                if prediction(calls,gene,2)!=prediction(native,gene,2):raise ValueError('Coarse prediction changed')
            for reason in decisions.values():
                key=method+'/'+reason;anchor_counts[key]=anchor_counts.get(key,0)+1
            results[donor,method]=dict(state='complete',calls=calls)
    references={}
    source=HERE.parent/'hla-analysis/source'
    for gene in GENES:
        with (source/'imgt/fasta'/(gene+'_gen.fasta')).open() as stream:
            references[gene]=reference_index(fasta(stream),gene)
    with gzip.open(source/'HGSVCv3_HLA_alleles.fasta.gz','rt') as stream:
        truth=derive(fasta(stream),[r['donor'] for r in cohort],references)
    rows,report=score(cohort,results,truth,ep['bootstrap_replicates'],ep['bootstrap_seed'])
    report.update(scoring_correction_sha256=sha(correction),anchor_decisions=anchor_counts,coarse_prediction_identity_verified=True,
                  inference_freeze_sha256=sha(OUT/'FROZEN_VALIDATION.json'),evaluation_freeze_sha256=sha(OUT/'EVALUATION_FREEZE.json'),
                  handoff_freeze_sha256=sha(OUT/'HANDOFF_FREEZE.json'),execution_ready_sha256=sha(ready_file),
                  truth_scope='Exact whole genomic HGSVC gene cores against IPD 3.65; two-field truth also genomic-derived')
    truth_rows=[dict(donor=d,gene=g,haplotype=r['haplotype'],status=r['status'],core_length=r['core_length'],exact_genomic_alleles=';'.join(r['exact_genomic_alleles'])) for (d,g),pair in truth.items() for r in pair]
    output.mkdir(parents=True)
    write_tsv(output/'gene_scores.tsv',rows);write_tsv(output/'summary.tsv',report['summary']);write_tsv(output/'truth_qc.tsv',truth_rows)
    (output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Independent HGSVC28 validation — corrected scoring','',
           'A documented post-unblinding correction removes the gene prefix from two-field truth to match the existing prediction parser. Predictions, eligibility, four-field scoring and statistical gates are unchanged. The original erroneous report is preserved; see ../SCORING_CORRECTION.json.','',
           'All 28 reserved donors; four methods; eight loci. Four-field truth requires exact whole-genomic matches. Two-field truth is also genomic-derived. Original two-field prediction identity was verified for each anchored method.','',
           '| Method | Fields | Correct | Eligible | Failed donors |','|---|---:|---:|---:|---:|']
    for r in report['summary']:
        lines.append('| {method} | {fields} | {correct} | {eligible} | {failed_donors} |'.format(**r))
    lines+=['','Primary combined protocol pass: **'+str(report['contrasts']['T1K']['protocol_pass'])+'**.','',
            'The original-T1K contrast measures the complete pipeline change. Incremental genomic-IPD and HPRC contrasts are exploratory; reference-only gains do not prove an Asian graph advantage. Per-ancestry intervals, failure accounting and paired contrasts are in report.json. The earlier 84-donor failed combined gate remains a separate result.','']
    (output/'REPORT.md').write_text('\n'.join(lines))
    print(json.dumps(dict(summary=report['summary'],primary_checks=report['contrasts']['T1K']['protocol_checks']),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--snapshot',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();run(a.snapshot,a.output)
