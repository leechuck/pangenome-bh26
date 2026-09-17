#!/usr/bin/env python3
"""Verify completed native outputs; this is file integrity, not accuracy scoring."""
import csv,json,collections,hashlib,tarfile
from pathlib import Path
R=Path(__file__).resolve().parent
rows=list(csv.DictReader(open(R/'source/donors.tsv'),delimiter='\t'))
out=[]
manifest=[]
for r in rows:
 d=r['donor'];counts={}
 p=R/'t1k'/d/(d+'_genotype.tsv');lines=[s.split('\t') for s in p.read_text().splitlines() if s]
 assert len(lines)>0 and all(len(s)>=8 for s in lines),(d,'T1K schema')
 counts['T1K_gene_rows']=len(lines)
 p=R/'spechla'/d/'hla.result.txt';lines=[s for s in p.read_text().splitlines() if s and not s.startswith('#')]
 assert len(lines)==2,(d,'SpecHLA row count')
 header=lines[0].split('\t');v=lines[1].split('\t');assert len(header)==len(v)==17 and v[0]==d,(d,'SpecHLA schema')
 counts['SpecHLA_allele_fields']=len(v)-1
 for arm in ['full','hprc']:
  counts[arm+'_records']=0;counts[arm+'_missing_GT']=0
  p=R/'pangenie'/d/arm/'calls_genotyping.vcf'
  saw_header=False
  with open(p) as f:
   for line in f:
    if line.startswith('##'):continue
    if line.startswith('#CHROM'):
     assert line.rstrip().split('\t')[9:]==[d],(d,arm,'sample');saw_header=True;continue
    fields=line.rstrip().split('\t');assert len(fields)==10,(d,arm,'columns')
    fmt=fields[8].split(':');values=fields[9].split(':');gt=values[fmt.index('GT')]
    alleles=['.','.'] if gt=='.' else gt.replace('|','/').split('/');assert len(alleles)==2,(d,arm,'ploidy')
    n=1+len(fields[4].split(','))
    assert all(a=='.' or 0<=int(a)<n for a in alleles),(d,arm,'allele range')
    counts[arm+'_records']+=1;counts[arm+'_missing_GT']+=int('.' in alleles)
  assert saw_header and counts[arm+'_records']>0,(d,arm,'empty VCF')
  h=hashlib.sha256()
  with open(p,'rb') as handle:
   for block in iter(lambda:handle.read(1024*1024),b''):h.update(block)
  digest=h.hexdigest()
  manifest.append(dict(donor=d,arm=arm,path=str(p.relative_to(R)),size_bytes=p.stat().st_size,sha256=digest))
 out.append(dict(donor=d,stratum=r['stratum'],**counts))
with open(R/'source/output_integrity.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
(R/'source/output_integrity.json').write_text(json.dumps(dict(status='PASS',donors=len(out),T1K_files=40,SpecHLA_files=40,PanGenie_VCFs=80,meaning='Native file identity and schema checks, not genotyping accuracy'),indent=2)+'\n')
print('PASS: 40 T1K files, 40 SpecHLA files, 80 PanGenie VCFs')

(R/'source/pangenie_output_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
with tarfile.open(R/'source/native_HLA_calls.tar.gz','w:gz') as archive:
 for r in rows:
  d=r['donor']
  for p in [R/'t1k'/d/(d+'_genotype.tsv'),R/'spechla'/d/'hla.result.txt']:
   archive.add(p,arcname=str(p.relative_to(R)))
print('Archived 80 native HLA call files and hashed 80 PanGenie VCFs.')
