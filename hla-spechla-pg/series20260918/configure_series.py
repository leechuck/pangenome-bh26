#!/usr/bin/env python3
"""Initialize isolated roots after the checked asset bundle and software are ready."""
from pathlib import Path
import shutil
ROOT=Path('/home/leechuck/hla/dogohla-series20260918')
OLD=Path('/home/leechuck/hla/spechla-pg')
for p in (OLD/'source').iterdir():
 target=ROOT/'source'/p.name
 if not target.exists():target.symlink_to(p)
for arm in ('full','hprc','asian_matched'):
 root=ROOT/'arms'/arm;root.mkdir(parents=True,exist_ok=True)
 for name,source in [('code',ROOT/'code'),('scripts',ROOT/'code'),('source',ROOT/'source'),('cohort.tsv',ROOT/'cohort.tsv'),('reads',Path('/home/leechuck/hla/hla-typer/reads'))]:
  t=root/name
  if not t.exists():t.symlink_to(source)
 (root/'phase').mkdir(exist_ok=True)
 if not (root/'phase/code').exists():(root/'phase/code').symlink_to(ROOT/'code')
(ROOT/'logs').mkdir(exist_ok=True)
Path('/home/leechuck/hla/hla-typer/reads').mkdir(exist_ok=True)
