#!/usr/bin/env python3
"""Reuse immutable, donor-free panel references with separate cohort outputs."""
import csv,json
from pathlib import Path
BASE=Path('/home/leechuck/hla/dogohla-series20260918')
ROOT=BASE/'gourraud';ROOT.mkdir(exist_ok=True)
source=ROOT/'source';source.mkdir(exist_ok=True)
for p in (BASE/'source').iterdir():
 if p.name=='donors_folds.tsv':continue
 t=source/p.name
 if not t.exists():t.symlink_to(p)
rows=list(csv.DictReader((ROOT/'cohort.tsv').open(),delimiter='\t'))
with (source/'donors_folds.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=['donor','fold','cram'],delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
for arm in ('full','hprc','asian_matched'):
 root=ROOT/'arms'/arm;root.mkdir(parents=True,exist_ok=True);prepared=BASE/'arms'/arm
 links={'code':BASE/'code','scripts':BASE/'code','source':source,'cohort.tsv':ROOT/'cohort.tsv','reads':Path('/home/leechuck/hla/hla-typer/reads'),'db':prepared/'db','graphs':prepared/'graphs'}
 for name,p in links.items():
  t=root/name
  if not t.exists():t.symlink_to(p)
 phase=root/'phase';phase.mkdir(exist_ok=True)
 for name in ('code','vendor','db','PREPARED.json'):
  t=phase/name
  if not t.exists():t.symlink_to(prepared/'phase'/name)
print('GOURRAUD_OUTPUT_ROOTS_READY')
