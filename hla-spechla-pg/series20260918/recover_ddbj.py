#!/usr/bin/env python3
"""Recover execution failures without changing truth, thresholds, or successful outputs."""
import argparse,csv,json,shutil,sys
from pathlib import Path
from experiment import completed
from native_long_control import run as native_run
from structural_overlay import run as overlay
from graph_diagnostic import diagnose
p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('mode',choices=['native-long','overlay']);p.add_argument('index',type=int);a=p.parse_args()
root=a.root.resolve();r=list(csv.DictReader((root/'cohort.tsv').open(),delimiter='\t'))[a.index];d=r['donor'];fold=int(r['fold'])
def ready(out):
 if (out/'COMPLETE').exists():
  m=json.loads((out/'manifest.json').read_text())
  if completed(out,m['configuration']):return True
 return False

def archive(out):
 if not out.exists():return
 dest=root/'attempts/20260918-execution-fixes'/out.relative_to(root)
 dest.parent.mkdir(parents=True,exist_ok=True)
 if dest.exists():raise FileExistsError('Previous retry exists: '+str(dest))
 shutil.move(out,dest)

if a.mode=='native-long':
 out=root/'runs'/d/'native-long'
 if ready(out):print('Previously verified complete',d);sys.exit(0)
 archive(out);archive(root/'runs'/d/'native-long_work')
 native_run(root,d,fold,4)
else:
 source=root/'phase/runs'/d/'AE-ipd365-pairing'
 if not ready(source):raise ValueError('Phasing not complete')
 for name,genes in [('DogoHLA-no-graph',[]),('DogoHLA',['DRB1'])]:
  out=root/'final/runs'/d/name
  if ready(out):continue
  archive(out)
  if genes:diagnose(root,root/'diagnostics'/d/'DRB1',fold,'DRB1',d,'pggb',1000000,4,300)
  overlay(source,root/'diagnostics'/d,d,out,genes,4)
