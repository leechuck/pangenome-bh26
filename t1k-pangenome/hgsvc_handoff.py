"""Audit the complete HGSVC output chain before fetching predictions."""
import argparse
import collections
import csv
import json
from pathlib import Path
import shlex
import subprocess
from evidence_io import sha
from launch_development_evidence import B,HOST,remote

HERE=Path(__file__).resolve().parent
OUT=HERE/'validation-next'
LOCATIONS={'T1K':'frozen-t1k-v1/{donor}', 'raw_genomic':'genomic-ipd-v1/{donor}',
 'raw_hprc':'graph-v1/calls/{donor}/hprc/unsampled',
 'raw_hprc_asian':'graph-v1/calls/{donor}/hprc_asian/unsampled',
 **{m:'anchored-v1/'+m+'/{donor}' for m in ('ipd_genome','graph_hprc','graph_hprc_asian')}}


def check_freezes():
    freeze=OUT/'FROZEN_VALIDATION.json'
    evaluation=OUT/'EVALUATION_FREEZE.json'
    ep=json.loads(evaluation.read_text())
    if ep['inference_freeze_sha256']!=sha(freeze):raise ValueError('Freeze mismatch')
    for name,digest in {**ep['code_sha256'],**ep['input_sha256']}.items():
        if sha(HERE.parent/name)!=digest:raise ValueError('Frozen evaluator/input changed: '+name)
    plan=json.loads(freeze.read_text())
    for name,digest in plan['code_sha256'].items():
        if sha(HERE/name)!=digest:raise ValueError('Frozen inference code changed: '+name)
    for name,digest in plan['metadata_sha256'].items():
        if sha(OUT/name)!=digest:raise ValueError('Frozen metadata changed: '+name)
    hp=json.loads((OUT/'HANDOFF_FREEZE.json').read_text())
    for name,digest in hp['code_sha256'].items():
        if sha(HERE/name)!=digest:raise ValueError('Handoff code changed: '+name)
    if hp['evaluation_freeze_sha256']!=sha(evaluation):raise ValueError('Handoff/evaluation freeze mismatch')
    return plan


def verify_reports(reports,donors,plan,prepared,bp,gp,bp_sha,gp_sha):
    expected={(d,m) for d in donors for m in LOCATIONS}
    indexed={(r['donor'],r['method']):r for r in reports}
    if len(reports)!=len(expected) or set(indexed)!=expected:raise ValueError('Incomplete or duplicated provenance grid')
    code=plan['code_sha256'];assets=plan['asset_sha256']
    def require(condition,reason):
        if not condition:raise ValueError(reason)
    for item in reports:
        if item['status']!='complete':continue
        d,m=item['donor'],item['method'];r=item['record']
        reads={str(i):r['reads_sha256'].get(str(i),r['reads_sha256'].get('r'+str(i)+'.fq.gz')) for i in (1,2)}
        require(r['donor']==d and reads=={str(i):prepared[d]['r'+str(i)+'.fq.gz'] for i in (1,2)},'Donor/read mismatch')
        if m=='T1K':
            require(r['plan_sha256']==bp_sha and r['recipe']==bp['recipe'] and r['driver_sha256']==code['run_hgsvc_baseline.py'],'Original T1K differs')
        elif m=='raw_genomic':
            require(r['validation_plan_sha256']==gp_sha and r['tool_sha256']==gp['tool_sha256'] and r['driver_sha256']==code['run_linear_control.py'] and r['validation_driver_sha256']==code['run_hgsvc_genomic.py'] and r['reference_manifest_sha256']==assets['t1k-references/genome-v2/COMPLETE.json'] and r['cohort_sha256']==plan['metadata_sha256']['cohort.tsv'],'Genomic control differs')
        elif m.startswith('raw_'):
            base=indexed[d,'raw_genomic']
            require(base['status']=='complete' and r['baseline_sha256']==base['manifest_sha256'],'Graph/native link mismatch')
            parameters={k:plan['parameters'][k] for k in ('temperature','locus_margin','minimum_fragments','minimum_gap','noise','internal_pairs','pair_specific')}
            parameters.update(alignment_only=True,native_locus_guard=True)
            require(r['parameters']==parameters and r['panel']==m.removeprefix('raw_'),'Graph parameters differ')
            for key,filename in [('driver_sha256','run_graph_pair_native_locus.py'),('model_sha256','graph_pair_refinement_native_locus.py'),('alignment_model_sha256','graph_pair_refinement_alignment_only.py'),('base_model_sha256','graph_pair_refinement.py')]:
                require(r[key]==code[filename],'Graph code differs: '+key)
            require(r['reference_manifest_sha256']==assets['references/observed-v2/COMPLETE.json'],'Graph reference differs')
        else:
            raw='raw_genomic' if m=='ipd_genome' else 'raw_'+m.removeprefix('graph_')
            original=indexed[d,'T1K'];candidate=indexed[d,raw]
            require(original['status']==candidate['status']=='complete','Anchor parents incomplete')
            require(r['original_sha256']==original['manifest_sha256'] and r['candidate_sha256']==candidate['manifest_sha256'],'Anchor parent link mismatch')
            require(r['driver_sha256']==code['run_hgsvc_anchor.py'] and r['anchor_sha256']==code['anchor_t1k_coarse.py'] and r['method']==m,'Anchor code/method differs')
    return all(r['status']=='complete' for r in reports)


def collect(fetch=False):
    plan=check_freezes()
    donors=[r['donor'] for r in csv.DictReader((OUT/'cohort.tsv').open(),delimiter='\t')]
    if len(donors)!=plan['donors'] or len(set(donors))!=len(donors):raise ValueError('Invalid cohort')
    root=B+'/hla/t1k-pangenome/hgsvc-validation-v1/validation-runs'
    script='''import json,pathlib,hashlib
root=pathlib.Path(ROOT_LITERAL)
locations=LOCATIONS_LITERAL
reports=[]
for donor in DONORS_LITERAL:
 for method,template in locations.items():
  relative=template.format(donor=donor);folder=root/relative
  complete=folder/'COMPLETE.json';manifest=folder/'manifest.json'
  item=dict(donor=donor,method=method,relative=relative,status='not_started')
  if manifest.exists():item['status']=json.load(manifest.open())['status']
  if complete.exists():
   raw=complete.read_bytes();record=json.loads(raw)
   assert record['status']=='complete'
   if manifest.exists():assert json.load(manifest.open())==record
   for name,digest in record['output_sha256'].items():
    assert pathlib.Path(name).name==name and name not in ('.','..')
    assert hashlib.sha256((folder/name).read_bytes()).hexdigest()==digest
   item.update(status='complete',record=record,manifest_sha256=hashlib.sha256(raw).hexdigest())
  reports.append(item)
print(json.dumps(reports))'''.replace('ROOT_LITERAL',repr(root)).replace('LOCATIONS_LITERAL',repr(LOCATIONS)).replace('DONORS_LITERAL',repr(donors))
    reports=json.loads(remote('python3 -c '+shlex.quote(script)))
    prepared={r['donor']:r['reads_sha256'] for r in json.loads((OUT/'READ_PREPARATION.json').read_text())['donors']}
    bp=OUT/'validation-baseline-plan.json';gp=OUT/'validation-genomic-plan.json'
    ready=verify_reports(reports,donors,plan,prepared,json.loads(bp.read_text()),json.loads(gp.read_text()),sha(bp),sha(gp))
    report=dict(ready=ready,states=dict(collections.Counter(r['status'] for r in reports)),runs=reports,
        inference_freeze_sha256=sha(OUT/'FROZEN_VALIDATION.json'),evaluation_freeze_sha256=sha(OUT/'EVALUATION_FREEZE.json'),
        handoff_freeze_sha256=sha(OUT/'HANDOFF_FREEZE.json'),scope='Metadata and output hashes; no prediction text inspected')
    (OUT/'execution-handoff.json').write_text(json.dumps(report,indent=2)+'\n')
    if fetch:
        if not ready:raise ValueError('Incomplete grid: failures require explicit audit before scoring')
        snapshot=HERE/'work/hgsvc-validation-snapshot'
        if snapshot.exists():raise FileExistsError('Inspect existing snapshot before fetching again')
        snapshot.mkdir(parents=True)
        for name in ('frozen-t1k-v1','genomic-ipd-v1','graph-v1/calls','anchored-v1'):
            dest=snapshot/name;dest.mkdir(parents=True,exist_ok=True)
            subprocess.run(['rsync','-az','--include=*/','--include=manifest.json','--include=COMPLETE.json','--include=*_genotype.tsv','--include=calls.json','--include=decisions.json','--exclude=*',HOST+':'+root+'/'+name+'/',str(dest)+'/'],check=True)
        (snapshot/'EXECUTION_READY.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(ready=ready,states=report['states'],fetched=fetch)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--fetch',action='store_true')
    collect(p.parse_args().fetch)
