"""Verify execution metadata before fetching any reserved prediction contents."""
import argparse
import collections
import csv
import json
import shlex
import subprocess
from pathlib import Path
from evidence_io import sha
from launch_development_evidence import HERE,B,HOST,remote

METHOD_DIRS={'T1K':'frozen-t1k-v1/{donor}', 'ipd_genome':'genomic-ipd-v1/{donor}',
             'graph_hprc':'graph-v1/calls/{donor}/hprc/unsampled',
             'graph_hprc_asian':'graph-v1/calls/{donor}/hprc_asian/unsampled'}


def verify_grid(reports,donors):
    expected={(d,m) for d in donors for m in METHOD_DIRS}
    if len(reports)!=len(expected):raise ValueError('Incomplete execution metadata grid')
    indexed={(r['donor'],r['method']):r for r in reports}
    if set(indexed)!=expected:raise ValueError('Duplicate or unexpected execution metadata')
    # Missing/pending/failed results must never be silently treated as complete.
    return indexed,all(r['status']=='complete' for r in reports)


def collect(fetch=False):
    frozen=HERE/'validation/FROZEN_VALIDATION.json';plan=json.loads(frozen.read_text())
    evaluation=HERE/'validation/EVALUATION_FREEZE.json';ep=json.loads(evaluation.read_text())
    handoff=HERE/'validation/HANDOFF_FREEZE.json';hp=json.loads(handoff.read_text())
    if hp['evaluation_freeze_sha256']!=sha(evaluation):raise ValueError('Handoff freeze mismatch')
    for filename,digest in hp['code_sha256'].items():
        if sha(HERE.parent/filename)!=digest:raise ValueError('Frozen handoff code changed')
    if ep['inference_freeze_sha256']!=sha(frozen):raise ValueError('Inference/evaluation freeze mismatch')
    for filename,digest in {**ep['code_sha256'],**ep['input_sha256']}.items():
        if sha(HERE.parent/filename)!=digest:raise ValueError('Frozen endpoint asset changed')
    cohort=HERE/'validation/cohort.tsv'
    if sha(cohort)!=plan['metadata_sha256']['cohort.tsv']:raise ValueError('Reserved cohort changed')
    donors=[r['donor'] for r in csv.DictReader(cohort.open(),delimiter='\t')]
    if len(donors)!=plan['donors'] or len(set(donors))!=len(donors):raise ValueError('Reserved cohort count changed')
    for d in donors:
        if not d or Path(d).name!=d or d in ('.','..'):raise ValueError('Unsafe donor')
    # Metadata and small output digests only. No prediction text is returned.
    script='''import pathlib,json,hashlib
root=pathlib.Path(ROOT)
donors=DONORS
locations=LOCATIONS
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
reports=[]
keep=('donor','reads_sha256','driver_sha256','model_sha256','reference_manifest_sha256',
      'baseline_sha256','parameters','panel','plan_sha256','recipe','cohort_sha256',
      'validation_plan_sha256','validation_driver_sha256','tool_sha256','output_sha256')
for donor in donors:
 for method,template in locations.items():
  relative=template.format(donor=donor);folder=root/relative
  completed=folder/'COMPLETE.json';manifest=folder/'manifest.json'
  item=dict(donor=donor,method=method,relative=relative,status='not_started')
  if manifest.exists():
   record=json.load(manifest.open());item['status']=record['status']
   if record['status']=='complete':
    assert completed.exists() and json.load(completed.open())==record
    for name,value in record['output_sha256'].items():
     assert pathlib.Path(name).name==name and name not in ('.','..')
     assert digest(folder/name)==value
    item.update(manifest_sha256=digest(completed),record={k:record[k] for k in keep if k in record})
  reports.append(item)
print(json.dumps(reports))
'''.replace('ROOT',repr(B+'/hla/t1k-pangenome/validation-runs')).replace('DONORS',repr(donors)).replace('LOCATIONS',repr(METHOD_DIRS))
    reports=json.loads(remote('python3 -c '+shlex.quote(script)))
    indexed,ready=verify_grid(reports,donors)
    preparation=HERE/'validation/READ_PREPARATION.json'
    if sha(preparation)!=plan['metadata_sha256']['READ_PREPARATION.json']:raise ValueError('Read metadata changed')
    prepared={r['donor']:r['reads_sha256'] for r in json.loads(preparation.read_text())['donors']}
    baseline_plan=HERE/'validation-baseline-plan.json';genomic_plan=HERE/'validation-genomic-plan.json'
    bp=json.loads(baseline_plan.read_text());gp=json.loads(genomic_plan.read_text())
    bl=json.loads((HERE/'validation-baseline-launch.json').read_text())['code_sha256']
    gl=json.loads((HERE/'validation-genomic-launch.json').read_text())['code_sha256']
    for item in reports:
        if item['status']!='complete':continue
        r=item['record'];d=item['donor'];m=item['method']
        reads={str(i):r['reads_sha256'].get(str(i),r['reads_sha256'].get(f'r{i}.fq.gz')) for i in (1,2)}
        if r['donor']!=d or reads!={str(i):prepared[d][f'r{i}.fq.gz'] for i in (1,2)}:
            raise ValueError('Completed run donor/read mismatch')
        if m=='T1K':
            if r['plan_sha256']!=sha(baseline_plan) or r['recipe']!=bp['recipe'] or r['driver_sha256']!=bl['run_validation_baseline.py']:
                raise ValueError('Original T1K execution differs from frozen baseline')
        elif m=='ipd_genome':
            if (r['validation_plan_sha256']!=sha(genomic_plan) or r['tool_sha256']!=gp['tool_sha256'] or
                r['driver_sha256']!=gl['run_linear_control.py'] or r['validation_driver_sha256']!=gl['run_validation_genomic.py'] or
                r['reference_manifest_sha256']!=plan['asset_sha256']['t1k-references/genome-v2/COMPLETE.json'] or
                r['cohort_sha256']!=plan['metadata_sha256']['cohort.tsv']):
                raise ValueError('Genomic control differs from frozen execution')
        else:
            control=indexed[d,'ipd_genome']
            if control['status']!='complete' or r['baseline_sha256']!=control['manifest_sha256']:
                raise ValueError('Graph/native provenance mismatch')
            expected={k:plan['parameters'][k] for k in ('temperature','locus_margin','minimum_fragments','minimum_gap','noise','internal_pairs','pair_specific')}
            if (r['parameters']!=expected or r['panel']!=m.removeprefix('graph_') or
                r['driver_sha256']!=plan['code_sha256']['run_graph_pair_refinement.py'] or
                r['model_sha256']!=plan['code_sha256']['graph_pair_refinement.py'] or
                r['reference_manifest_sha256']!=plan['asset_sha256']['references/observed-v2/COMPLETE.json']):
                raise ValueError('Graph inference differs from frozen candidate')
    report=dict(ready=ready,states=dict(collections.Counter(r['status'] for r in reports)),
                inference_freeze_sha256=sha(frozen),evaluation_freeze_sha256=sha(evaluation),runs=reports,
                scope='Execution metadata and output hashes only; no prediction contents fetched unless explicitly requested after readiness')
    (HERE/'validation/execution-handoff.json').write_text(json.dumps(report,indent=2)+'\n')
    if fetch:
        if not ready:raise ValueError('Not ready: pending or failed runs require completion/failure audit before unblinding')
        snapshot=HERE/'work/validation-snapshot';snapshot.mkdir(parents=True,exist_ok=True)
        # Fixed subdirectories only; never fetch intermediate alignments or reads.
        for name in ('frozen-t1k-v1','genomic-ipd-v1','graph-v1/calls'):
            dest=snapshot/name;dest.mkdir(parents=True,exist_ok=True)
            subprocess.run(['rsync','-az','--include=*/','--include=manifest.json','--include=COMPLETE.json',
                '--include=*_genotype.tsv','--include=calls.json','--include=decisions.json','--exclude=*',
                HOST+':'+B+'/hla/t1k-pangenome/validation-runs/'+name+'/',str(dest)+'/'],check=True)
        (snapshot/'EXECUTION_READY.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(ready=ready,states=report['states'],predictions_fetched=bool(fetch)),indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fetch',action='store_true')
    collect(p.parse_args().fetch)
