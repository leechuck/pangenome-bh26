#!/usr/bin/env python3
"""Fetch compact, completed result files and retain manifests; no alignments or reads."""
import io,subprocess,tarfile,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
script=r'''
import csv,io,json,sys,tarfile
from pathlib import Path
root=Path('/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/dogohla-series20260918')
COHORT_ROOT
donors=list(csv.DictReader((root/'cohort.tsv').open(),delimiter='\t'))
with tarfile.open(fileobj=sys.stdout.buffer,mode='w|gz') as tar:
 for name in ('cohort.tsv','PREPARED.json','code/MANIFEST.json'):
  if (root/name).exists():tar.add(root/name,arcname=name)
 for r in donors:
  d=r['donor'];paths=[root/'t1k4'/d,root/'arms/full/runs'/d/'native',root/'arms/full/runs'/d/'native-long']
  paths += [root/'arms'/a/'final/runs'/d/m for a in ('full','hprc','asian_matched') for m in ('DogoHLA','DogoHLA-no-graph')]
  for p in paths:
   for n in ('manifest.json','COMPLETE','TERMINAL_FAILURE.json'):
    if (p/n).exists():tar.add(p/n,arcname=str((p/n).relative_to(root)))
   if (p/'COMPLETE').exists():
    m=json.loads((p/'COMPLETE').read_text())
    for n in m['outputs']:
     if Path(n).is_absolute() or '..' in Path(n).parts:raise ValueError('Unsafe output path')
     tar.add(p/n,arcname=str((p/n).relative_to(root)))
'''
extended='--gourraud' in sys.argv
script=script.replace('COHORT_ROOT',"root=root/'gourraud'" if extended else '')
r=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','hohndor@ilogin.ibex.kaust.edu.sa','python3 -'],input=script.encode(),capture_output=True,check=True,timeout=90)
out=HERE/'snapshot'/('gourraud' if extended else '');out.mkdir(parents=True,exist_ok=True)
with tarfile.open(fileobj=io.BytesIO(r.stdout),mode='r:gz') as tar:tar.extractall(out,filter='data')
print('Fetched',len(r.stdout),'bytes')
