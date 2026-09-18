"""Advance frozen HGSVC jobs from Slurm state; never read predictions or truth."""
import json
from pathlib import Path
from launch_development_evidence import remote,submit
from launch_hgsvc_validation import OUT,CODE


def main():
    ledger=OUT/'EXECUTION_LAUNCH.json'
    if not ledger.exists():raise FileNotFoundError(ledger)
    record=json.loads(ledger.read_text())
    if record.get('pending_submission'):raise RuntimeError('Reconcile pending submission before advancing')
    jobs=record['jobs']
    raw=remote('sacct -X -n -P -j '+','.join(jobs.values())+' --format=JobID,State,ExitCode')
    states={}
    for line in raw.splitlines():
        fields=line.split('|')
        if len(fields)>=3:states[fields[0]]=(fields[1],fields[2])
    widths={'baseline':1,'genomic':1,'coverage':1,'native':1,'bootstrap':8,'bootstrap_evidence':8,'calibration':1,'mapping':16,'evidence':16,'join':2,'refine':2,'anchor':3}
    def done(stage):
        return stage in jobs and all(states.get(jobs[stage]+'_'+str(i))==('COMPLETED','0:0') for i in range(28*widths[stage]))
    def save():ledger.write_text(json.dumps(record,indent=2)+'\n')
    next_stage=None
    if 'mapping' not in jobs and done('calibration'):next_stage=('mapping',4,'12G',[])
    elif 'mapping' in jobs and 'evidence' not in jobs:next_stage=('evidence',1,'8G',['--dependency=aftercorr:'+jobs['mapping']])
    elif 'join' not in jobs and done('evidence') and done('native'):next_stage=('join',1,'8G',[])
    elif 'refine' not in jobs and done('join'):next_stage=('refine',1,'8G',[])
    elif 'anchor' not in jobs and done('refine') and done('baseline'):next_stage=('anchor',1,'4G',[])
    if next_stage:
        stage,cpus,mem,extra=next_stage
        extra+=['--array=0-'+str(28*widths[stage]-1)+'%200','--partition=batch,debug']
        name='t1k-hgsvc28-'+stage
        command=['python3',CODE+'/run_hgsvc_stage.py','--assets',CODE,'--freeze-sha256',record['freeze_sha256'],'--index','$SLURM_ARRAY_TASK_ID','--stage',stage]
        record['pending_submission']=dict(stage=stage,name=name,command=command,extra=extra);save()
        job=submit(name,cpus,mem,'02:00:00',command,extra)
        jobs[stage]=job;record['pending_submission']=None;save()
        print('SUBMITTED',stage,job,flush=True)
    failures={job:state for job,state in states.items() if state[0].split()[0] in ('FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','NODE_FAIL')}
    print(json.dumps(dict(completed_stages=[s for s in jobs if done(s)],failures=failures,terminal=done('anchor'))),flush=True)
    return done('anchor'),bool(failures)


if __name__=='__main__':main()
