#!/usr/bin/env python3
"""Aggregate counts and paired donor bootstrap; no variant independence assumption."""
import csv,collections,json
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parent/'results'
def load(name):
    rs=list(csv.DictReader(open(R/name),delimiter='\t'))
    for r in rs:
        for k in ['n','correct','called','site_present','represented','nonref_n','nonref_correct']:
            if k in r:r[k]=int(r[k])
    return rs

def write(name,rows):
    with open(R/name,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def main():
    rs=load('panel_per_donor.tsv')
    if (R/'linear_per_donor.tsv').exists():rs+=load('linear_per_donor.tsv')
    groups=collections.defaultdict(list)
    for r in rs:groups[r['stratum'],r['universe'],r['variant_class']].append(r)
    sums=[];pairs=[]
    rng=np.random.default_rng(20260917)
    for (stratum,universe,kind),data in sorted(groups.items()):
        arms=collections.defaultdict(dict)
        for r in data:arms[r['arm']][r['donor']]=r
        for arm,dr in sorted(arms.items()):
            c={k:(sum(r[k] for r in dr.values()) if all(k in r for r in dr.values()) else '') for k in ['n','correct','called','site_present','represented','nonref_n','nonref_correct']}
            sums.append(dict(stratum=stratum,universe=universe,variant_class=kind,arm=arm,donors=len(dr),**c,accuracy_pct=100*c['correct']/c['n'],nonref_accuracy_pct=100*c['nonref_correct']/c['nonref_n'] if c['nonref_n'] else ''))
        for baseline in ['hprc','linear']:
            if baseline not in arms:continue
            ids=sorted(arms['full']);assert ids==sorted(arms[baseline]);n=len(ids)
            for endpoint,num,den in [('all','correct','n'),('nonref','nonref_correct','nonref_n')]:
                aa=np.array([[arms['full'][s][num],arms['full'][s][den],arms[baseline][s][num],arms[baseline][s][den]] for s in ids],dtype=float)
                # Own all-truth denominators are identical; all shared denominators must be too.
                assert np.all(aa[:,1]==aa[:,3])
                if np.any(aa[:,1]==0):continue
                delta=100*(aa[:,0]/aa[:,1]-aa[:,2]/aa[:,3])
                boot=aa[rng.integers(0,n,(10000,n))].sum(axis=1)
                bd=100*(boot[:,0]/boot[:,1]-boot[:,2]/boot[:,3]);lo,hi=np.quantile(bd,[.025,.975])
                totals=aa.sum(axis=0);point=100*(totals[0]/totals[1]-totals[2]/totals[3])
                pairs.append(dict(stratum=stratum,universe=universe,variant_class=kind,baseline=baseline,endpoint=endpoint,donors=n,full_minus_baseline_pp=point,ci95_low_pp=lo,ci95_high_pp=hi,donors_better=int((delta>0).sum()),donors_tied=int((delta==0).sum()),donors_worse=int((delta<0).sum())))
    write('summary.tsv',sums);write('paired_bootstrap.tsv',pairs)
    (R/'bootstrap.json').write_text(json.dumps(dict(seed=20260917,replicates=10000,unit='donor, within ancestry; paired arms; percentile CI; conditional on fixed sites/folds; no multiplicity correction'),indent=2)+'\n')
if __name__=='__main__':main()
