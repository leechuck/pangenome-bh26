#!/usr/bin/env python3
"""Submit packed Gourraud arrays after each method passes its execution check."""
import json,shlex,time
from launch_ibex import B,HERE,R,remote,submit

LEDGER=HERE/'ibex-launch.json'
ARMS=('native','full','hprc','asian_matched')

def ready_reads():
    script="""import csv,json,subprocess
from pathlib import Path
b=Path(BPATH)/'hla'
rows=list(csv.DictReader((b/'dogohla-series20260918/gourraud/cohort.tsv').open(),delimiter='\\t'))
ready=[]
for i,r in enumerate(rows):
 p=b/'hla-typer/reads'/r['donor']
 if (p/'COMPLETE').exists() and all((p/f'r{h}.fq.gz').is_file() and (p/f'r{h}.fq.gz').stat().st_size>0 for h in (1,2)):ready.append(i)
print(json.dumps(ready))
""".replace('BPATH',repr(B))
    return set(json.loads(remote('python3 -c '+shlex.quote(script))))

def state(job):
    lines=remote('sacct -X -n -P -j '+job+' --format=JobID,State').splitlines()
    return next((s.split('|')[1] for s in lines if s.split('|')[0]==job),'UNKNOWN')

def main():
    from pack_pending import submit_cohort
    while True:
        jobs=json.loads(LEDGER.read_text())['jobs']
        ready=ready_reads()
        print('Recruitment ready',len(ready),'/946',flush=True)
        if len(ready)!=946:
            time.sleep(60);continue
        for arm in ARMS:
            key='gourraud_packed_'+arm
            if key in jobs:continue
            gate=jobs['native_smoke'] if arm=='native' else jobs[arm+'_smoke']
            current=state(gate)
            if current!='COMPLETED':
                print('Waiting for',arm,'execution gate',gate,current,flush=True);continue
            queued_count=int(remote('squeue -r -h -u hohndor | wc -l'))
            if queued_count+119>1900:
                print('Waiting for queue capacity',queued_count,flush=True);continue
            panel='full' if arm=='native' else arm
            submit_cohort(key,'native' if arm=='native' else 'dogohla',R+'/gourraud/arms/'+panel,panel)
            jobs=json.loads(LEDGER.read_text())['jobs']
        if all('gourraud_packed_'+arm in jobs for arm in ARMS):
            print('All Gourraud typing arrays submitted',flush=True);return
        time.sleep(60)

if __name__=='__main__':main()
