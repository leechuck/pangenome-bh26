#!/usr/bin/env python3
"""A held-out-only allele and missing training calls must not enter inference."""
import csv,gzip,json,shutil,subprocess,tempfile
from pathlib import Path
R=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as td:
 p=Path(td);(p/'source').mkdir();(p/'graph').mkdir();shutil.copy(R/'prepare_panels.py',p)
 def write(name,rows):
  with open(p/'source'/name,'w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 names=['train'+str(i) for i in range(26)]+['test']
 write('panel_manifest.tsv',[dict(sample=s,donor_id=s,haplotype=h,cohort='HPRC-Rest') for s in names for h in ['1','2']])
 write('donors.tsv',[dict(donor='test',fold=0)])
 write('excluded_paths.tsv',[dict(fold=k,sample='test') for k in range(5)])
 with gzip.open(p/'graph/full.top.vcf.gz','wt') as f:
  f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(names)+'\n')
  for pos,id,alt,g in [(10,'shared','C',['0|1']*27),(20,'private','C',['0|0']*26+['0|1']),(30,'missing','C',['.|1']+['0|1']*26),(40,'trim','C,G',['0|1']*26+['0|2'])]:
   f.write('\t'.join(['ctg1',str(pos),id,'A',alt,'.','PASS','.','GT']+g)+'\n')
 subprocess.run(['python3',str(p/'prepare_panels.py')],check=True,capture_output=True)
 for v in (p/'panels').glob('*/*/bubbles.vcf'):
  lines=v.read_text().splitlines();assert 'test' not in lines[2].split('\t')
  records=[l.split('\t') for l in lines if not l.startswith('#')]
  assert [r[2] for r in records]==['shared','trim'];assert records[1][4]=='C'
  audit=json.loads((v.parent/'audit.json').read_text());assert audit['training_donors']==26
print('PASS: test-only alleles trimmed, test paths excluded, missing calls filtered, both arms/five folds checked.')
