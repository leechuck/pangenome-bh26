#!/usr/bin/env python3
"""Standalone figures from frozen metrics; shared-site and recovery endpoints separate."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R=Path(__file__).resolve().parent/'results'
rs=list(csv.DictReader(open(R/'summary.tsv'),delimiter='\t'))
def value(s,u,k,a):return float(next(r['accuracy_pct'] for r in rs if (r['stratum'],r['universe'],r['variant_class'],r['arm'])==(s,u,k,a)))
fig,axs=plt.subplots(1,3,figsize=(12.4,4.6))
for ax,(universe,kind,title,ylim) in zip(axs,[('shared_sites','SNV','SNV concordance\nShared graph sites',(98.8,100)),('shared_sites','truth_SV_length','SV-bearing genotype concordance\nShared graph sites',(65,90)),('all_truth','truth_SV_length','Correct SV-bearing genotypes recovered\nAll eligible truth sites',(55,90))]):
    x=np.arange(2)
    for j,(a,label,col) in enumerate([('hprc','HPRC-only','#697582'),('full','Full panel','#087f8c')]):
        vals=[value(s,universe,kind,a) for s in ['EAS','SAS']]
        bars=ax.bar(x+(j-.5)*.34,vals,width=.32,label=label,color=col)
        for b,v in zip(bars,vals):ax.text(b.get_x()+b.get_width()/2,v+(ylim[1]-ylim[0])*.015,f'{v:.2f}',ha='center',fontsize=9)
    ax.set_xticks(x,['East Asian','South Asian']);ax.set_ylim(*ylim);ax.set_title(title,fontsize=11);ax.set_ylabel('Percent');ax.spines[['top','right']].set_visible(False)
axs[0].legend(loc='lower left',frameon=False)
fig.suptitle('Added haplotypes help at shared sites; current filtering loses coverage',fontsize=13)
fig.text(.5,.015,'20 donors per ancestry · held-out assembly labels on shared full-graph topology · exact unphased diploid sequence matches\nSV-bearing: ≥50 bp truth-allele length change; not SV event precision/recall. Note truncated y-axes.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.09,1,.94));fig.savefig(R/'panel_comparison.png',dpi=180);fig.savefig(R/'panel_comparison.pdf');plt.close(fig)
