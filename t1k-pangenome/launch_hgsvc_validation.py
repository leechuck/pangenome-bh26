"""Stage the HGSVC freeze and submit baseline, calibration and bootstrap jobs."""
import json
import shlex
import subprocess
from pathlib import Path
from build_graph import sha
from launch_development_evidence import B,HOST,remote,submit

HERE=Path(__file__).resolve().parent
OUT=HERE/'validation-next'
CODE='/home/leechuck/hla/t1k-pangenome/hgsvc-validation-v1/code'


def main():
    ledger=OUT/'EXECUTION_LAUNCH.json'
    if ledger.exists():raise FileExistsError('Inspect existing jobs before any resubmission')
    frozen=OUT/'FROZEN_VALIDATION.json'
    plan=json.loads(frozen.read_text())
    evaluation=json.loads((OUT/'EVALUATION_FREEZE.json').read_text())
    if evaluation['inference_freeze_sha256']!=sha(frozen):raise ValueError('Evaluation freeze mismatch')
    files=[HERE/n for n in plan['code_sha256']]+[OUT/n for n in plan['metadata_sha256']]+[frozen]
    hashes={p.name:sha(p) for p in files}
    if len(hashes)!=len(files):raise ValueError('Duplicate staged basename')
    for n,h in {**plan['code_sha256'],**plan['metadata_sha256']}.items():
        if hashes[n]!=h:raise ValueError('Local freeze changed: '+n)
    destination=B+'/hla/t1k-pangenome/hgsvc-validation-v1/code'
    remote('mkdir '+shlex.quote(destination))
    subprocess.run(['scp',*map(str,files),HOST+':'+destination+'/'],check=True)
    actual=remote('sha256sum '+' '.join(shlex.quote(destination+'/'+n) for n in hashes))
    if {Path(line.split()[1]).name:line.split()[0] for line in actual.splitlines()}!=hashes:
        raise ValueError('Transfer mismatch')
    record=dict(freeze_sha256=sha(frozen),transferred_sha256=hashes,jobs={},pending_submission=None,
                scope='All 28 reserved donors; no truth or typing output inspection')
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    save()
    stages=[('baseline',1,4,'8G',[]),('genomic',1,4,'8G',[]),('coverage',1,1,'8G',[]),
            ('native',1,4,'12G',[('aftercorr','genomic')]),
            ('bootstrap',8,4,'12G',[('afterok','coverage')]),
            ('bootstrap_evidence',8,1,'8G',[('aftercorr','bootstrap')]),
            ('calibration',1,8,'16G',[('afterok','bootstrap_evidence')])]
    for stage,width,cpus,mem,parents in stages:
        name='t1k-hgsvc28-'+stage
        extra=['--array=0-'+str(plan['donors']*width-1)+'%200','--partition=batch,debug']
        if parents:extra+=['--dependency='+','.join(kind+':'+record['jobs'][parent] for kind,parent in parents)]
        command=['python3',CODE+'/run_hgsvc_stage.py','--assets',CODE,'--freeze-sha256',sha(frozen),
                 '--index','$SLURM_ARRAY_TASK_ID','--stage',stage]
        record['pending_submission']=dict(stage=stage,name=name,command=command,extra=extra);save()
        job=submit(name,cpus,mem,'02:00:00',command,extra)
        record['jobs'][stage]=job;record['pending_submission']=None;save()
        print(stage,job,flush=True)


if __name__=='__main__':main()
