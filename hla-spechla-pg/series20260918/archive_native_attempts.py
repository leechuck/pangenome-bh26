from pathlib import Path
import csv,json,shutil
root=Path('/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/dogohla-series20260918')
for r in csv.DictReader((root/'cohort.tsv').open(),delimiter='\t'):
 for n in ('native','native_work'):
  p=root/'arms/full/runs'/r['donor']/n
  if not p.exists():continue
  if (p/'COMPLETE').exists():raise ValueError('Do not archive successful output: '+str(p))
  dest=root/'attempts/empty-ipd-records'/p.relative_to(root)
  dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists():raise FileExistsError(dest)
  shutil.move(p,dest)
print('FAILED_NATIVE_ATTEMPTS_PRESERVED')
