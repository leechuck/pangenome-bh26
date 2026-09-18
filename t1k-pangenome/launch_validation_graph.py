"""Queue frozen graph validation; complete pilot gates the remaining cohort."""
import hashlib
import json
import shlex
import subprocess
from launch_development_evidence import HERE,B,HOST,remote,submit


if __name__=='__main__':
    ledger=HERE/'validation-graph-launch.json'
    if ledger.exists():raise FileExistsError('Inspect existing validation jobs before resubmission')
    frozen=HERE/'validation/FROZEN_VALIDATION.json';plan=json.loads(frozen.read_text())
    digest=hashlib.sha256(frozen.read_bytes()).hexdigest()
    files=[HERE/f for f in plan['code_sha256']]
    files += [HERE/'validation'/f if (HERE/'validation'/f).exists() else HERE/f for f in plan['metadata_sha256']]
    files += [frozen]
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    if any(hashes[f]!=h for f,h in {**plan['code_sha256'],**plan['metadata_sha256']}.items()):
        raise ValueError('Local files differ from freeze')
    directory=B+'/hla/t1k-pangenome/validation-graph-v1-code'
    remote('mkdir '+shlex.quote(directory))
    subprocess.run(['scp',*[str(p) for p in files],HOST+':'+directory+'/'],check=True)
    actual=remote('sha256sum '+' '.join(directory+'/'+p.name for p in files))
    if {r.split()[1].rsplit('/',1)[-1]:r.split()[0] for r in actual.splitlines()}!=hashes:
        raise ValueError('Transferred validation assets differ')
    root='/home/leechuck/hla/t1k-pangenome';code=root+'/validation-graph-v1-code'
    base=['python3',code+'/run_validation_graph.py','--assets',code,'--freeze-sha256',digest,
          '--index','$SLURM_ARRAY_TASK_ID']
    genomic=json.loads((HERE/'validation-genomic-launch.json').read_text())
    # Finished jobs may age out of Slurm; use verified completion metadata then.
    script=('import csv,json,pathlib; root=pathlib.Path('+repr(B+'/hla/t1k-pangenome')+'); '
            'donors=[r["donor"] for r in csv.DictReader((root/"validation-genomic-v1-code/cohort.tsv").open(),delimiter="\\t")]; '
            'paths=[root/"validation-runs/genomic-ipd-v1"/d/"COMPLETE.json" for d in donors]; '
            'print(json.dumps([p.exists() and json.load(p.open())["status"]=="complete" for p in paths]))')
    completed=json.loads(remote('python3 -c '+shlex.quote(script)))
    if len(completed)!=84:raise ValueError('Genomic cohort count mismatch')
    record=dict(freeze_sha256=digest,transferred_sha256=hashes,jobs={},
                scope='84 reserved donors; two unsampled panels; no prediction fetch or truth scoring')
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    def stage(group,name,width,cpus,mem,parents=(),external=()):
        start,end=(0,width-1) if group=='pilot' else (width,84*width-1)
        extra=['--array='+str(start)+'-'+str(end)+'%300','--partition=batch,debug']
        deps=[kind+':'+record['jobs'][group][p] for kind,p in parents]+list(external)
        if group=='cohort':deps+=['afterok:'+record['jobs']['pilot']['refine']]
        if deps:extra+=['--dependency='+','.join(deps)]
        job=submit('t1k-val-'+group+'-'+name,cpus,mem,'02:00:00',base+['--stage',name],extra)
        record['jobs'][group][name]=job;save()
    for group in ('pilot','cohort'):
        record['jobs'][group]={};save()
        stage(group,'coverage',1,1,'8G')
        # The genomic control pilot is a scalar job; the expansion is indexed 1-83.
        genomic_dep=([] if completed[0] else ['afterok:'+genomic['pilot_job']]) if group=='pilot' else (
            [] if all(completed[1:]) else ['aftercorr:'+genomic['array_job']])
        stage(group,'native',1,4,'12G',external=genomic_dep)
        stage(group,'bootstrap',8,4,'12G',[('afterok','coverage')])
        stage(group,'bootstrap_evidence',8,1,'8G',[('aftercorr','bootstrap')])
        stage(group,'calibration',1,8,'16G',[('afterok','bootstrap_evidence')])
        if group=='cohort':
            record['dispatch_policy']='Remaining stages submitted by advance_validation_graph.py after predecessor completion to respect queued-job cap'
            save()
            continue
        stage(group,'mapping',16,4,'12G',[('afterok','calibration')])
        stage(group,'evidence',16,1,'8G',[('aftercorr','mapping')])
        stage(group,'join',2,1,'8G',[('afterok','evidence'),('afterok','native')])
        stage(group,'refine',2,1,'12G',[('aftercorr','join')])
    print(json.dumps(record['jobs'],indent=2))
