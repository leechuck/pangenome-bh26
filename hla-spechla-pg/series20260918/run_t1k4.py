#!/usr/bin/env python3
"""Run T1K with explicit four-field reporting and preserve the real genotype table."""
import argparse,csv,json,os,subprocess,time
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'code'))
from experiment import sha,start,finish,failed
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('index',type=int);a=p.parse_args()
row=list(csv.DictReader((a.root/'cohort.tsv').open(),delimiter='\t'))[a.index];d=row['donor']
reads=Path('/home/leechuck/hla/hla-typer/reads')/d
assert (reads/'COMPLETE').exists()
out=a.root/'t1k4'/d;db=Path('/home/leechuck/hla/hla-typer/t1kdb/hla365_dna_seq.fa')
exe=Path('/home/leechuck/hla/mm/envs/typing/bin/run-t1k')
config=dict(donor=d,alleleDigitUnits=4,alleleDelimiter=':',preset='hla-wgs',ipd='3.65.0',database_sha256=sha(db),tool_sha256=sha(exe),reads={str(i):sha(reads/f'r{i}.fq.gz') for i in (1,2)})
if not start(out,config):sys.exit(0)
try:
 env=dict(os.environ);env['PATH']=str(exe.parent)+':'+env['PATH']
 args=[str(exe),'-1',str(reads/'r1.fq.gz'),'-2',str(reads/'r2.fq.gz'),'--preset','hla-wgs','-f',str(db),'-t','4','--alleleDigitUnits','4','--alleleDelimiter',':','-o',str(out/d)]
 with (out/'run.log').open('w') as log:subprocess.run(args,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 result=out/(d+'_genotype.tsv')
 if not result.is_file() or not result.stat().st_size:raise ValueError('Missing genotype table')
 finish(out,config,{result.name:sha(result)})
except BaseException as exc:failed(out,exc);raise
