#!/usr/bin/env python3
"""Paired donor bootstrap for frozen variant and HLA comparisons."""
import csv,collections,hashlib,json
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parent/'results'
def read(name):
 with open(R/name) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows):
 with open(R/name,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def bootstrap(arms,baseline,key):
 ids=sorted(arms['full']);assert ids==sorted(arms[baseline]);n=len(ids)
 aa=np.array([[*arms['full'][s],*arms[baseline][s]] for s in ids],dtype=float)
 assert np.all(aa[:,1]==aa[:,3])
 if np.any(aa[:,1]==0):return None
 seed=int.from_bytes(hashlib.sha256(('20260917|'+str(key)+'|'+baseline).encode()).digest()[:8],'little')
 rng=np.random.default_rng(seed);bb=aa[rng.integers(0,n,(10000,n))].sum(axis=1)
 delta=100*(bb[:,0]/bb[:,1]-bb[:,2]/bb[:,3]);lo,hi=np.quantile(delta,[.025,.975]);total=aa.sum(axis=0)
 pd=aa[:,0]/aa[:,1]-aa[:,2]/aa[:,3]
 return dict(baseline=baseline,donors=n,full_minus_baseline_pp=100*(total[0]/total[1]-total[2]/total[3]),ci95_low_pp=lo,ci95_high_pp=hi,donors_better=int((pd>0).sum()),donors_tied=int((pd==0).sum()),donors_worse=int((pd<0).sum()))
def main():
 if (R/'variant_per_donor.tsv').exists():
  groups=collections.defaultdict(list)
  for r in read('variant_per_donor.tsv'):groups[r['stratum'],r['universe'],r['variant_class']].append(r)
  sums=[];pairs=[]
  for key,rs in sorted(groups.items()):
   s,u,k=key
   for arm in sorted({r['arm'] for r in rs}):
    rr=[r for r in rs if r['arm']==arm];c={col:sum(int(r[col]) for r in rr) for col in ['n','site_present','called','correct','nonref_n','nonref_correct']}
    sums.append(dict(stratum=s,universe=u,variant_class=k,arm=arm,**c,accuracy_pct=100*c['correct']/c['n'],nonref_accuracy_pct=100*c['nonref_correct']/c['nonref_n'] if c['nonref_n'] else ''))
   for endpoint,num,den in [('all','correct','n'),('nonref','nonref_correct','nonref_n')]:
    arms=collections.defaultdict(dict)
    for r in rs:arms[r['arm']][r['donor']]=(int(r[num]),int(r[den]))
    for baseline in sorted(set(arms)-{'full'}):
     b=bootstrap(arms,baseline,(*key,endpoint))
     if b:pairs.append(dict(stratum=s,universe=u,variant_class=k,endpoint=endpoint,**b))
  write('variant_summary.tsv',sums);write('variant_paired.tsv',pairs)
 if (R/'hla_scores.tsv').exists():
  groups=collections.defaultdict(list)
  named=(R/'named_HLA_scores.tsv').exists()
  for r in read('named_HLA_scores.tsv' if named else 'hla_scores.tsv'):groups[r['stratum'],r['endpoint']].append(r)
  pairs=[]
  for (s,e),rs in sorted(groups.items()):
   arms=collections.defaultdict(lambda:collections.defaultdict(lambda:[0,0]))
   for r in rs:
    method='full' if r['method']==('Locityper_full' if named else 'graph_HLA_full') else r['method'];arms[method][r['donor']][0]+=int(r['correct']);arms[method][r['donor']][1]+=1
   for baseline in sorted(set(arms)-{'full'}):
    b=bootstrap(arms,baseline,(s,e,'HLA'))
    pairs.append(dict(stratum=s,endpoint=e,**b))
  write('hla_paired.tsv',pairs)
 (R/'bootstrap.json').write_text(json.dumps(dict(replicates=10000,seed='SHA256 of 20260917 + endpoint key + baseline',unit='paired donor resampling within ancestry',caveat='Conditional on fixed test donors/sites/folds and overlapping training panels; no multiple-testing correction'),indent=2)+'\n')
if __name__=='__main__':main()
