#!/usr/bin/env python3
"""Standalone figure from scored TSVs; distinct truth sources stay separate."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R=Path(__file__).resolve().parent/'results'
def read(name):
 with open(R/name) as f:return list(csv.DictReader(f,delimiter='\t'))
v=read('variant_summary.tsv');h=read('named_HLA_summary.tsv')
fig,axes=plt.subplots(1,3,figsize=(14,5),layout='constrained')
colors=['#939393','#4477AA','#228877','#CC6677']
sv=[('old_full','Old full'),('hprc','Repaired HPRC'),('full','Repaired full')]
methods=[('Locityper_hprc','HPRC + Locityper'),('Locityper_full','Full + Locityper'),('T1K','T1K'),('SpecHLA','SpecHLA')]
for i,(method,label) in enumerate(sv):
 vals=[float(next(r['accuracy_pct'] for r in v if r['stratum']==s and r['universe']=='all_truth' and r['variant_class']=='truth_SV_length' and r['arm']==method)) for s in ['EAS','SAS']]
 bars=axes[0].bar(np.arange(2)+(i-1)*.25,vals,width=.24,label=label,color=colors[i]);axes[0].bar_label(bars,fmt='%.1f',fontsize=8)
axes[0].set_xticks([0,1],['EAS (1,465)','SAS (1,577)'])
axes[0].set_title('SV-bearing genotype recovery\nGraph-derived assembly truth',fontsize=11)
for ai,endpoint,strata in [(1,'assembly_exact_CDS_twofield',['EAS','SAS']),(2,'experimental_twofield',['EAS'])]:
 ax=axes[ai]
 for i,(method,label) in enumerate(methods):
  vals=[float(next(r['genotype_accuracy_pct'] for r in h if r['stratum']==s and r['endpoint']==endpoint and r['method']==method and r['gene']=='ALL')) for s in strata]
  bars=ax.bar(np.arange(len(strata))+(i-1.5)*.19,vals,width=.18,label=label,color=[colors[1],colors[2],'#DD9933',colors[3]][i]);ax.bar_label(bars,fmt='%.1f',fontsize=8)
 ax.set_xticks(np.arange(len(strata)),['EAS (159)','SAS (160)'] if ai==1 else ['EAS (39 genotypes, 8 donors)'])
 ax.set_title('Two-field HLA allele-pair accuracy\n'+('Assembly-derived labels' if ai==1 else 'Historical experimental labels'),fontsize=11)
for ax in axes:
 ax.set_ylim(0,105);ax.set_yticks([0,20,40,60,80,100]);ax.set_ylabel('Correct genotypes (%)')
 ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
 ax.legend(loc='upper center',bbox_to_anchor=(.5,-.13),frameon=False,fontsize=8,ncol=2)
fig.suptitle('Callable-site repair and named HLA typing — 40 held-out 1000 Genomes donors',fontsize=13)
fig.savefig(R/'comparison.svg');fig.savefig(R/'comparison.pdf');fig.savefig(R/'comparison.png',dpi=180)

svg=R/'comparison.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
