#!/usr/bin/env python3
"""Reuse frozen training membership/labels; retain true region boundaries and shared sequence."""
import csv,gzip,json,hashlib,collections
from pathlib import Path
from prepare_hla import GENES,fasta,labels
R=Path(__file__).resolve().parent
catalog=collections.defaultdict(list)
for row in csv.DictReader(open(R/'hla_source/sequence_catalogue.tsv'),delimiter='\t'):catalog[row['name']].append(row)
reference=next(iter(fasta(R/'graph/reference.fa').values()))
seqs={g:fasta(R/f'hla_source/HLA-{g}.fa') for g in GENES}
regions={g:{r['name']:r for r in csv.DictReader(open(R/f'hla_source/HLA-{g}.regions.tsv'),delimiter='\t')} for g in GENES}
audit=[]
for fold in range(5):
 for arm in ['full','hprc']:
  samples=(R/f'panels/fold{fold}/{arm}/samples.txt').read_text().splitlines();path=R/f'locityper/panels/fold{fold}/{arm}';path.mkdir(parents=True,exist_ok=True)
  with open(path/'targets.bed','w') as bed,gzip.open(path/'labels.tsv.gz','wt') as out:
   w=csv.DictWriter(out,fieldnames=['gene','haplotype','labels','source_haplotypes','sha256'],delimiter='\t',lineterminator='\n');w.writeheader()
   for g in GENES:
    rr=regions[g][f'GRCh38#0#HLA-{g}'];lo,hi=map(int,rr['region'].rsplit(':',1)[1].split('-'));group=collections.defaultdict(list)
    for s,h in [('GRCh38','0')]+[(s,h) for s in samples for h in ['1','2']]:
     name=f'{s}#{h}#HLA-{g}';seq=seqs[g].get(name)
     if seq and set(seq)<=set('ACGT'):
      if rr['strand']=='-':seq=seq.translate(str.maketrans('ACGT','TGCA'))[::-1]
      group[seq].append(name)
    assert reference[lo-1:hi] in group,('GRCh38 orientation/extraction mismatch',g)
    fa=path/f'HLA-{g}.fa'
    with open(fa,'w') as f:
     for i,(seq,names) in enumerate(sorted(group.items())):
      hap=f'P{i:04d}';f.write(f'>{hap}\n{seq}\n');ls=set()
      for name in names:
       rs=catalog[name];ls.update(labels(rs[0]) if len(rs)==1 else set())
      w.writerow(dict(gene=g,haplotype=hap,labels=';'.join(sorted(ls)),source_haplotypes=';'.join(names),sha256=hashlib.sha256(seq.encode()).hexdigest()))
    bed.write(f'ctg1\t{lo-1}\t{hi}\tHLA-{g}\t{fa}\n')
    audit.append(dict(fold=fold,arm=arm,gene=g,haplotypes=len(group),training_donors=len(samples)))
(R/'locityper/panel_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
