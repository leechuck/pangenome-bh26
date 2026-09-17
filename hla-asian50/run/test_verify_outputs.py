#!/usr/bin/env python3
"""Accept PanGenie's missing GT, but reject out-of-range allele indices."""
import csv,json,shutil,subprocess,tempfile
from pathlib import Path
R=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as tmp:
 p=Path(tmp);(p/'source').mkdir();shutil.copy(R/'verify_outputs.py',p)
 with open(p/'source/donors.tsv','w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['donor','stratum'])
  for i in range(40):w.writerow(['s'+str(i),'EAS' if i<20 else 'SAS'])
 for i in range(40):
  d='s'+str(i);t=p/'t1k'/d;t.mkdir(parents=True);(t/(d+'_genotype.tsv')).write_text('HLA-A\t2\tA*01:01\t1\t1\tA*02:01\t1\t1\n')
  t=p/'spechla'/d;t.mkdir(parents=True);(t/'hla.result.txt').write_text('\t'.join(['Sample']+['allele'+str(j) for j in range(16)])+'\n'+'\t'.join([d]+['A*01:01']*16)+'\n')
  for arm in ['full','hprc']:
   t=p/'pangenie'/d/arm;t.mkdir(parents=True);(t/'calls_genotyping.vcf').write_text('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+d+'\nctg1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t.\n')
 subprocess.run(['python3',str(p/'verify_outputs.py')],check=True,capture_output=True)
 rr=list(csv.DictReader(open(p/'source/output_integrity.tsv'),delimiter='\t'))
 assert len(rr)==40 and all(r['full_missing_GT']=='1' and r['hprc_missing_GT']=='1' for r in rr)
 bad=p/'pangenie/s0/full/calls_genotyping.vcf';bad.write_text(bad.read_text().replace('\tGT\t.\n','\tGT\t2/2\n'))
 r=subprocess.run(['python3',str(p/'verify_outputs.py')],capture_output=True,text=True)
 assert r.returncode!=0 and 'allele range' in r.stderr
print('PASS: missing GT preserved/counts as no-call; out-of-range GT rejected.')
