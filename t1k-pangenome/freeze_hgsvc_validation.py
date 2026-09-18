"""Bind the selected HGSVC candidate to verified reads and executable code."""
import csv
import datetime
import json
from pathlib import Path
import shlex
from build_graph import sha
from launch_development_evidence import B,remote

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
OUT = HERE/'validation-next'


def save(path, value):
    path.write_text(json.dumps(value,indent=2)+'\n')


def main():
    frozen = OUT/'FROZEN_VALIDATION.json'
    if frozen.exists():raise FileExistsError(frozen)
    selected = json.loads((OUT/'PRIMARY_CANDIDATE.json').read_text())
    for name,digest in selected['code_sha256'].items():
        if sha(HERE/name)!=digest:raise ValueError('Selected code changed: '+name)
    for name,digest in selected['metadata_sha256'].items():
        if sha(OUT/name)!=digest:raise ValueError('Selected metadata changed: '+name)
    for name,digest in selected['truth_source_sha256'].items():
        if sha(REPO/name)!=digest:raise ValueError('Selected truth source changed: '+name)
    cohort = list(csv.DictReader((OUT/'cohort.tsv').open(),delimiter='\t'))
    root = B+'/hla/t1k-pangenome/hgsvc-validation-v1/validation-reads'
    # Metadata only on login: full FASTQ verification was performed inside Slurm.
    code = '''import json,pathlib,hashlib
root=pathlib.Path(ROOT_LITERAL)
rows=[]
for donor in DONORS_LITERAL:
 p=root/donor/'VERIFIED.json'
 if not p.exists():
  rows.append(dict(donor=donor,status='not_verified'));continue
 raw=p.read_bytes();record=json.loads(raw)
 if (root/donor/'manifest.json').read_bytes()!=raw:raise ValueError('Manifest differs')
 rows.append(dict(donor=donor,status=record['status'],record=record,verified_sha256=hashlib.sha256(raw).hexdigest()))
print(json.dumps(rows))'''.replace('ROOT_LITERAL',repr(root)).replace('DONORS_LITERAL',repr([r['donor'] for r in cohort]))
    rows = json.loads(remote('python3 -c '+shlex.quote(code)))
    pending = [r['donor'] for r in rows if r['status']!='complete']
    if pending:
        print(json.dumps(dict(status='waiting_for_verified_reads',complete=len(rows)-len(pending),pending=pending)))
        return
    launch = json.loads((OUT/'READ_PREPARATION_LAUNCH.json').read_text())
    for row in rows:
        record=row['record']
        if (record['donor']!=row['donor'] or record['sample_ids']!=[row['donor']] or
            record['pairs']<=0 or record['cohort_sha256']!=sha(OUT/'cohort.tsv') or
            record['driver_sha256']!=launch['code_and_input_sha256']['prepare_validation_reads.py'] or
            record['script_sha256']!=launch['code_and_input_sha256']['recruit.sh']):
            raise ValueError('Read provenance differs: '+row['donor'])
    prep=dict(status='read_preparation_complete',cohort_sha256=sha(OUT/'cohort.tsv'),
              donors=[dict(donor=r['donor'],pairs=r['record']['pairs'],reads_sha256=r['record']['output_sha256'],verified_sha256=r['verified_sha256']) for r in rows])
    save(OUT/'READ_PREPARATION.json',prep)
    for kind in ('baseline','genomic'):
        name='validation-'+kind+'-plan.json'
        plan=json.loads((HERE/name).read_text())
        plan.update(cohort_sha256=sha(OUT/'cohort.tsv'),preparation_sha256=sha(OUT/'READ_PREPARATION.json'),reservation_sha256=sha(OUT/'RESERVATION.json'),scope='Independent HGSVC28 execution; predictions remain unscored')
        save(OUT/name,plan)
    old=json.loads((HERE/'validation/FROZEN_VALIDATION.json').read_text())
    names=set(old['code_sha256'])-{'run_validation_graph.py'}
    names.update(['run_hgsvc_stage.py','run_hgsvc_baseline.py','run_hgsvc_genomic.py','run_hgsvc_graph.py','run_hgsvc_anchor.py','run_linear_control.py','decode_t1k.py','run_graph_pair_native_locus.py','graph_pair_refinement_native_locus.py','graph_pair_refinement_alignment_only.py','anchor_t1k_coarse.py'])
    metadata=['cohort.tsv','RESERVATION.json','PRIMARY_CANDIDATE.json','READ_PREPARATION.json','validation-baseline-plan.json','validation-genomic-plan.json','PROTOCOL.md','GRAPH_SOURCE_AUDIT.json']
    # Reuse audited training references, never HGSVC held-out assembly sequences.
    assets=old['asset_sha256']
    remote_root=B+'/hla/t1k-pangenome'
    actual=remote('sha256sum '+' '.join(shlex.quote(remote_root+'/'+n) for n in assets))
    got={line.split()[1][len(remote_root)+1:]:line.split()[0] for line in actual.splitlines()}
    if got!=assets:raise ValueError('Remote reference manifest changed')
    now=datetime.datetime.now(datetime.timezone.utc).isoformat()
    plan=dict(status='frozen_before_outcome_inspection',frozen_at=now,donors=len(cohort),
              primary=selected['primary'],parameters=old['parameters'],
              code_sha256={n:sha(HERE/n) for n in sorted(names)},
              metadata_sha256={n:sha(OUT/n) for n in metadata},asset_sha256=assets,
              output_root='/home/leechuck/hla/t1k-pangenome/hgsvc-validation-v1',
              limitations='Known-family exclusion only; exact genomic assembly truth eligibility; see PROTOCOL.md')
    save(frozen,plan)
    original_eval=json.loads((HERE/'validation/EVALUATION_FREEZE.json').read_text())
    eval_names=set(original_eval['code_sha256'])-{'t1k-pangenome/evaluate_reserved.py'}
    eval_names.update(['t1k-pangenome/evaluate_hgsvc.py','t1k-pangenome/hgsvc_truth.py'])
    evaluation=dict(status='endpoint_code_frozen_before_outcome_inspection',frozen_at=now,
        inference_freeze_sha256=sha(frozen),code_sha256={n:sha(REPO/n) for n in sorted(eval_names)},
        input_sha256=selected['truth_source_sha256'],primary='graph_hprc_asian minus T1K',
        secondary_exploratory=['graph_hprc_asian minus ipd_genome','graph_hprc_asian minus graph_hprc'],
        bootstrap_replicates=10000,bootstrap_seed=20260918,
        failure_handling=original_eval['failure_handling'])
    save(OUT/'EVALUATION_FREEZE.json',evaluation)
    save(OUT/'HANDOFF_FREEZE.json',dict(status='handoff_frozen_before_outcome_inspection',frozen_at=now,
        evaluation_freeze_sha256=sha(OUT/'EVALUATION_FREEZE.json'),
        code_sha256={n:sha(HERE/n) for n in ('hgsvc_handoff.py','score_hgsvc_snapshot.py')}))
    print(json.dumps(dict(status='frozen',donors=len(cohort),freeze_sha256=sha(frozen))))


if __name__=='__main__':main()
