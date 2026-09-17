#!/usr/bin/env python3
"""Score named HLA calls with frozen truth/label rules; retain failed PanGenie experiment."""
import csv,gzip,json,collections,argparse
from pathlib import Path
from prepare_hla import GENES
from score_hla import norm,compare,write,table

def main():
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.run;a.out.mkdir(parents=True,exist_ok=True)
 donors=table(r/'source/donors.tsv');mapping={};panelsets={}
 for fold in range(5):
  for arm in ['full','hprc']:
   with gzip.open(r/f'locityper/panels/fold{fold}/{arm}/labels.tsv.gz','rt') as f:rs=list(csv.DictReader(f,delimiter='\t'))
   for x in rs:mapping[fold,arm,x['gene'],x['haplotype']]={norm(v) for v in x['labels'].split(';') if v}
   for g in GENES:panelsets[fold,arm,g]=set().union(*(mapping[fold,arm,x['gene'],x['haplotype']] for x in rs if x['gene']==g))
 predictions={};details=[]
 for d in donors:
  s=d['donor'];fold=int(d['fold'])
  for arm in ['full','hprc']:
   for g in GENES:
    path=r/f'locityper/calls/{s}/{arm}/loci/HLA-{g}/res.json.gz'
    with gzip.open(path,'rt') as f:res=json.load(f)
    raw=res.get('genotype') or '';haplotypes=raw.split(',') if isinstance(raw,str) else raw
    pair=[None,None]
    if len(haplotypes)==2:
     for i,hap in enumerate(haplotypes):
      ls=mapping[fold,arm,g,hap];pair[i]=next(iter(ls)) if len(ls)==1 and None not in ls else None
    method='Locityper_'+arm;predictions[s,g,method]=tuple(pair)
    details.append(dict(donor=s,stratum=d['stratum'],gene=g,method=method,haplotypes=raw,allele1=pair[0] or '',allele2=pair[1] or '',called=int(all(pair)),quality=res.get('quality',''),warnings=json.dumps(res.get('warnings',[])),total_reads=res.get('total_reads',''),unexpl_reads=res.get('unexpl_reads','')))
 write(a.out/'locityper_predictions.tsv',details)
 original=table(r/'results/hla_scores.tsv');scores=[]
 for x in original:
  y=dict(x);y['method']=y['method'].replace('graph_HLA_','PanGenie_whole_locus_');scores.append(y)
 truth=[x for x in original if x['method']=='T1K'];ds={d['donor']:d for d in donors}
 for x in truth:
  s=x['donor'];g=x['gene'];ts=[set(v.split('/')) for v in x['truth'].split(';')]
  for arm in ['full','hprc']:
   y=dict(x);method='Locityper_'+arm;pair=predictions[s,g,method];called,correct,matches=compare(pair,ts)
   y.update(method=method,prediction=';'.join(v or 'NO_CALL' for v in pair),called=called,correct=correct,allele_matches=matches,truth_labels_in_panel=int(all(t&panelsets[int(ds[s]['fold']),arm,g] for t in ts)))
   scores.append(y)
 write(a.out/'named_HLA_scores.tsv',scores)
 groups=collections.defaultdict(list)
 for x in scores:
  for g in [x['gene'],'ALL']:groups[x['stratum'],x['endpoint'],x['method'],g].append(x)
 summary=[]
 for (s,e,m,g),rs in sorted(groups.items()):
  c=sum(int(x['correct']) for x in rs);ac=sum(int(x['allele_matches']) for x in rs)
  summary.append(dict(stratum=s,endpoint=e,method=m,gene=g,donors=len({x['donor'] for x in rs}),genotypes=len(rs),called=sum(int(x['called']) for x in rs),correct=c,genotype_accuracy_pct=100*c/len(rs),alleles=2*len(rs),allele_matches=ac,allele_accuracy_pct=50*ac/len(rs)))
 write(a.out/'named_HLA_summary.tsv',summary);print(json.dumps([x for x in summary if x['gene']=='ALL'],indent=2))
if __name__=='__main__':main()
