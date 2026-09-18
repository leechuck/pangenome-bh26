#!/usr/bin/env python3
"""Submit checked preparation and matched 64-donor experiments on IBEX."""
import json,shlex,subprocess
from pathlib import Path
HERE=Path(__file__).resolve().parent
B='/ibex/scratch/projects/c2014/rob/dogohla-benchmark'
R='/home/leechuck/hla/dogohla-series20260918'
HOST='hohndor@ilogin.ibex.kaust.edu.sa'

def remote(cmd):
 r=subprocess.run(['ssh','-o','BatchMode=yes',HOST,cmd],capture_output=True,text=True,check=True)
 return r.stdout.strip()

def submit(name,cpus,mem,wall,command,extra=(),wrap=False):
 args=['sbatch','--parsable','--partition=batch','--account=pi-hohndor','--chdir='+B,'--job-name='+name,'--cpus-per-task='+str(cpus),'--mem='+mem,'--time='+wall,'--output='+B+'/'+name+'-%A_%a.log',*extra]
 if not wrap and command[0]=='python3':command=['/home/leechuck/hla/mm/envs/asian50-spechla/bin/python3',*command[1:]]
 if not wrap and '$SLURM_ARRAY_TASK_ID' in command:
  command=['bash','-c',shlex.join(command).replace("'$SLURM_ARRAY_TASK_ID'", '"$SLURM_ARRAY_TASK_ID"')]
 args += ['--wrap='+command] if wrap else [B+'/run_ibex.sbatch',*command]
 return remote(shlex.join(args)).split(';')[0]

def main():
 ledger=HERE/'ibex-launch.json'
 if ledger.exists():raise FileExistsError('Already submitted; inspect ledger before changing jobs')
 m=json.loads((HERE/'custom-transfer-complete.json').read_text())
 digest=m['sha256']
 extract=f"cat {B}/custom-parts/part[0-9][0-9][0-9][0-9] > {B}/custom.tar.gz && echo '{digest}  {B}/custom.tar.gz' | sha256sum -c - && tar xzf {B}/custom.tar.gz -C {B}/hla && bash {B}/run_ibex.sbatch bash {R}/series/prepare_ibex.sh"
 jobs={}
 def save():ledger.write_text(json.dumps(dict(workdir=B,jobs=jobs,scope='Matched 32 Asian / 32 non-Asian; graph-panel arms and four-field T1K. Native-long and extended Gourraud cohort pending separate launch.',custom_sha256=digest),indent=2)+'\n')
 jobs['prepare']=submit('dogo-prepare365',8,'16G','00:45:00',extract,['--dependency=afterok:52032352'],True);save()
 jobs['graphs']=submit('dogo-graphs',8,'24G','02:00:00',['bash',R+'/series/prepare_arm.sh'],['--array=0-2','--dependency=afterok:'+jobs['prepare']]);save()
 recruitment=(HERE/'ibex-recruit-job.txt').read_text().strip()
 deps=['--array=0-63','--dependency=afterok:'+jobs['prepare']+',aftercorr:'+recruitment]
 jobs['native']=submit('dogo-native365',4,'8G','01:30:00',['python3',R+'/code/dogohla.py','--root',R+'/arms/full','--native-only','--threads','4','--cohort-index','$SLURM_ARRAY_TASK_ID'],deps)
 save()
 jobs['t1k4']=submit('dogo-t1k4',4,'8G','01:00:00',['python3',R+'/series/run_t1k4.py',R,'$SLURM_ARRAY_TASK_ID'],deps);save()
 for i,arm in enumerate(('full','hprc','asian_matched')):
  jobs[arm]=submit('dogo-'+arm,4,'8G','02:00:00',['python3',R+'/code/dogohla.py','--root',R+'/arms/'+arm,'--threads','4','--cohort-index','$SLURM_ARRAY_TASK_ID'],['--array=0-63','--export=ALL,ARM='+arm,'--dependency=afterok:'+jobs['graphs']+'_'+str(i)+',aftercorr:'+recruitment]);save()
 print(ledger.read_text())

if __name__=='__main__':main()
