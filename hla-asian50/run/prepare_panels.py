#!/usr/bin/env python3
"""Create conservative, complete-genotype panels from full-graph top-level bubbles.

Global graph construction is transductive; test family paths/alleles never enter
these reference panels. Missing calls are filtered, never changed to reference.
"""
import csv,gzip,json,collections,re
from pathlib import Path
R=Path(__file__).resolve().parent
rows=list(csv.DictReader(open(R/'source/panel_manifest.tsv'),delimiter='\t'))
donors=list(csv.DictReader(open(R/'source/donors.tsv'),delimiter='\t'))
exclude=list(csv.DictReader(open(R/'source/excluded_paths.tsv'),delimiter='\t'))
by=collections.defaultdict(list)
for r in rows:by[r['sample']].append(r)
# Canonical donors only, with two input haplotypes; do not count alias assemblies twice.
valid={s for s,rs in by.items() if len(rs)==2 and s==rs[0]['donor_id'] and {r['haplotype'] for r in rs}=={'1','2'}}
contexts=[]
for k in range(5):
 ex={r['sample'] for r in exclude if int(r['fold'])==k}
 for arm in ['full','hprc']:
  ss={s for s in valid-ex if arm=='full' or by[s][0]['cohort'].startswith('HPRC')}
  path=R/'panels'/f'fold{k}'/arm;path.mkdir(parents=True,exist_ok=True)
  contexts.append(dict(fold=k,arm=arm,samples=ss,path=path,stats=collections.Counter(),last_end=-1))
with gzip.open(R/'graph/full.top.vcf.gz','rt') as f:
 for line in f:
  if line.startswith('##'):continue
  if line.startswith('#CHROM'):
   columns=line.rstrip().split('\t');samples=columns[9:]
   for c in contexts:
    c['idx']=[(i,s) for i,s in enumerate(samples) if s in c['samples']]
    assert len(c['idx'])>20
    (c['path']/'samples.txt').write_text('\n'.join(s for _,s in c['idx'])+'\n')
    c['file']=open(c['path']/'bubbles.vcf','w')
    c['map']=gzip.open(c['path']/'allele_map.tsv.gz','wt')
    c['map'].write('ID\toriginal_allele_indices\n')
    c['file'].write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype">\n'+'\t'.join(columns[:9]+[s for _,s in c['idx']])+'\n')
   continue
  v=line.rstrip().split('\t');alleles=[v[3]]+v[4].split(',');pos=int(v[1]);fmt=v[8].split(':');gi=fmt.index('GT')
  for c in contexts:
   st=c['stats'];st['input_sites']+=1
   gts=[v[i+9].split(':')[gi] for i,_ in c['idx']]
   if not all(re.fullmatch(r'\d+\|\d+',g) for g in gts):st['missing_or_nonphased']+=1;continue
   used={int(a) for g in gts for a in g.split('|')};used.add(0)
   if len(used)==1:st['no_training_alt']+=1;continue
   if any(re.search('[^ACGT]',alleles[a]) for a in used):st['non_ACGT']+=1;continue
   if pos<=c['last_end']:st['overlapping']+=1;continue
   order=sorted(used);remap={a:i for i,a in enumerate(order)}
   # Do not normalize/change coordinates across folds; keep an explicit allele map.
   genos=['|'.join(str(remap[int(a)]) for a in g.split('|')) for g in gts]
   c['file'].write('\t'.join(v[:4]+[','.join(alleles[a] for a in order[1:]),'.','PASS','.','GT']+genos)+'\n')
   c['map'].write(v[2]+'\t'+','.join(map(str,order))+'\n')
   c['last_end']=pos+len(v[3])-1;st['retained_sites']+=1
   if any(abs(len(alleles[a])-len(v[3]))>=50 for a in order[1:]):st['retained_SV_length_change_sites']+=1
   elif all(len(alleles[a])==1 for a in order):st['retained_SNV_sites']+=1
for c in contexts:
 c['file'].close();c['map'].close();assert c['stats']['retained_sites']>0
 summary=dict(fold=c['fold'],arm=c['arm'],training_donors=len(c['idx']),counts=dict(c['stats']),caveat='Global graph topology; complete phased training genotypes only; no independent SNV/SV truth established yet')
 (c['path']/'audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary)
