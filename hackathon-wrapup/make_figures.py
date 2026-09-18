#!/usr/bin/env python3
"""Figures for the BioHackathon 2026 wrap-up slides. Every number is read from result tables in this
repository; nothing is typed in by hand except GRCh38 gene coordinates for the C4 annotation.

Outputs (figures/):
  fig_cohorts.png          754 haplotypes by cohort (slide 1)
  fig_graph_density.png    node density and haplotype sharing along the reference path (slide 2)
  fig_variant_classes.png  PanGenie genotype concordance per variant class, full vs HPRC-only, and SNVs vs linear (slide 2)
  fig_hla_dev.png          Stage B per-gene two-field accuracy vs T1K on 40 dev donors; SpecHLA step-A sequence accuracy (slide 3)
"""
import csv,gzip,collections,sys,statistics
from pathlib import Path
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

R=Path(__file__).resolve().parent.parent          # repository root
OUT=Path(__file__).resolve().parent/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
FULL,HPRC,ASIAN,LIN,T1K,OTHER='#00838f','#6b7280','#e08a1e','#b03a48','#3b6fb6','#8e6bb3'

def tsv(p,**kw):return list(csv.DictReader(open(p),delimiter='\t',**kw))

# ---------------------------------------------------------------- slide 1: cohorts
def fig_cohorts():
    c=collections.Counter(r['cohort'] for r in tsv(R/'hla-typer/source/haplotype_donors.tsv'))
    order=[('APR','Arab (APR)'),('JaSaPaGe-Saudi','Saudi (JaSaPaGe)'),('JaSaPaGe-Japanese','Japanese (JaSaPaGe)'),
           ('HPRC-Japanese','Japanese (HPRC r2)'),('HPRC-EastAsian','other East Asian (HPRC r2)'),('KPanRef-Korean','Korean (K-PanRef graph)'),
           ('CPC-Chinese','Chinese (CPC graph)'),('HPRC-Jewish','Ashkenazi (HPRC r2)'),('HPRC-Rest','rest of HPRC r2'),('REF','GRCh38, CHM13')]
    vals=[c[k] for k,_ in order];labels=[l for _,l in order]
    asian={'APR','JaSaPaGe-Saudi','JaSaPaGe-Japanese','HPRC-Japanese','HPRC-EastAsian','KPanRef-Korean','CPC-Chinese'}
    cols=[FULL if k in asian else HPRC for k,_ in order]
    fig,ax=plt.subplots(figsize=(5.2,4.2))
    y=np.arange(len(order))[::-1]
    ax.barh(y,vals,color=cols)
    for yi,v in zip(y,vals):ax.text(v+4,yi,str(v),va='center',fontsize=10)
    ax.set_yticks(y);ax.set_yticklabels(labels,fontsize=10);ax.set_xlabel('MHC haplotypes')
    ax.set_xlim(0,max(vals)*1.15)
    n=sum(vals);na=sum(c[k] for k in asian)
    ax.set_title(f'{n} haplotypes; {na} Asian or Arab (teal)',fontsize=11,loc='left')
    fig.tight_layout();fig.savefig(OUT/'fig_cohorts.png',dpi=200);plt.close(fig)
    print('cohorts',n,na)

# ---------------------------------------------------------------- slide 2: graph along the reference
OFF=28_410_700   # ctg1 position + OFF = chr6 (GRCh38)
def fig_graph_density():
    rows=[]
    for r in csv.DictReader(gzip.open(R/'hla-typer/results/graph_scan/ref_nodes.tsv.gz','rt'),delimiter='\t'):
        rows.append((int(r['ref_start0']),int(r['ref_end0']),int(r['length']),int(r['haplotypes'])))
    L=rows[-1][1];W=20_000;nb=L//W+1
    nodes=np.zeros(nb);share=np.zeros(nb);bases=np.zeros(nb)
    for s,e,l,h in rows:
        b=s//W;nodes[b]+=1;share[b]+=h*l;bases[b]+=l
    share=np.where(bases>0,share/np.maximum(bases,1)/754*100,np.nan)
    x=(np.arange(nb)*W+W/2+OFF)/1e6
    fig,(a1,a2)=plt.subplots(2,1,figsize=(9,4.0),sharex=True,gridspec_kw={'hspace':0.12})
    a1.fill_between(x,nodes/(W/1000),color=FULL,alpha=0.85,lw=0);a1.set_ylabel('reference nodes\nper kb')
    a2.fill_between(x,share,100,color=LIN,alpha=0.75,lw=0);a2.set_ylabel('haplotypes sharing\nreference nodes (%)');a2.set_ylim(40,100)
    a2.set_xlabel('GRCh38 chr6 position (Mb)')
    spec=tsv(R/'hla-typer/source/loci_spec.tsv')
    genes=[(r['locus'].replace('HLA-',''),int(r['chr6_start1']),int(r['chr6_end1'])) for r in spec]
    genes.append(('C4A/C4B',31_982_057,32_035_418))   # GRCh38 RefSeq C4A..C4B span (RCCX module)
    genes.sort(key=lambda g:g[1]);top=nodes.max()/(W/1000)
    for i,(name,s,e) in enumerate(genes):
        for a in (a1,a2):a.axvspan(s/1e6,e/1e6,color='k',alpha=0.08,lw=0)
        dx={'DRB345':-0.09,'DQ':0.09,'DRB1':0.0}.get(name,0.0);lvl={'DRB345':1.08,'DRB1':1.30,'DQ':1.08}.get(name,1.08 if i%2 else 1.22)
        a1.text((s+e)/2e6+dx,top*lvl,name,ha='center',va='bottom',fontsize=8.5)
    a1.set_ylim(0,top*1.45)
    a1.set_title(f'Minigraph-Cactus graph, 754 haplotypes: {len(rows):,} nodes on the {L/1e6:.2f} Mb GRCh38 path, 20 kb windows',fontsize=10.5,loc='left',pad=14)
    fig.tight_layout();fig.savefig(OUT/'fig_graph_density.png',dpi=200);plt.close(fig)
    print('graph windows',nb,'median share',np.nanmedian(share))

# ---------------------------------------------------------------- slide 2: variant classes
def fig_variant_classes():
    rows=tsv(R/'hla-asian50/refined/results/variant_summary.tsv')
    def get(stratum,universe,cls,arm,col='accuracy_pct'):
        for r in rows:
            if (r['stratum'],r['universe'],r['variant_class'],r['arm'])==(stratum,universe,cls,arm):return float(r[col]),int(r['n'])
        raise KeyError((stratum,universe,cls,arm))
    classes=[('SNV','SNV'),('small_complex','small\ncomplex'),('SV_length_site','SV site'),('truth_SV_length','SV-bearing\ngenotype')]
    fig,axes=plt.subplots(1,3,figsize=(12.5,3.1),gridspec_kw={'width_ratios':[2.2,2.2,1.3],'wspace':0.3})
    for ax,strat in zip(axes[:2],['EAS','SAS']):
        xs=np.arange(len(classes));w=0.38
        for i,(arm,col,lab) in enumerate([('hprc',HPRC,'HPRC-only panel'),('full',FULL,'full panel')]):
            v=[get(strat,'frozen_HPRC_sites',c,arm)[0] for c,_ in classes]
            b=ax.bar(xs+(i-0.5)*w,v,w,color=col,label=lab)
            for xi,vi in zip(xs+(i-0.5)*w,v):ax.text(xi,vi+0.6,f'{vi:.1f}',ha='center',fontsize=8.2)
        ax.set_xticks(xs);ax.set_xticklabels([l for _,l in classes],fontsize=9);ax.set_ylim(60,102)
        n=[get(strat,'frozen_HPRC_sites',c,'full')[1] for c,_ in classes]
        ax.set_title(f'{"East" if strat=="EAS" else "South"} Asian, 20 donors\nn = {n[0]:,} / {n[1]:,} / {n[2]:,} / {n[3]:,} genotypes',fontsize=9.5)
        if strat=='EAS':ax.set_ylabel('exact diploid genotype concordance (%)');ax.legend(fontsize=8.5,loc='lower center',frameon=False,bbox_to_anchor=(0.5,-0.32),ncol=2)
    ax=axes[2];xs=np.arange(2);w=0.26
    for i,(arm,col,lab) in enumerate([('linear',LIN,'published linear\n(NYGC GATK)'),('hprc',HPRC,'HPRC-only'),('full',FULL,'full panel')]):
        v=[get(s,'frozen_threeway_SNV','SNV',arm)[0] for s in ('EAS','SAS')]
        ax.bar(xs+(i-1)*w,v,w,color=col,label=lab)
        for xi,vi in zip(xs+(i-1)*w,v):ax.text(xi,vi+0.005,f'{vi:.2f}',ha='center',fontsize=7.5,rotation=90,va='bottom')
    ax.set_xticks(xs);ax.set_xticklabels(['East Asian','South Asian'],fontsize=9);ax.set_ylim(99.5,100.05)
    n=get('EAS','frozen_threeway_SNV','SNV','full')[1]
    ax.set_title(f'SNVs shared by all three callsets\nn = {n:,} genotypes per stratum',fontsize=9.5);ax.legend(fontsize=7.5,loc='upper left',frameon=False,bbox_to_anchor=(0,1.0))
    fig.savefig(OUT/'fig_variant_classes.png',dpi=200,bbox_inches='tight');plt.close(fig)

# ---------------------------------------------------------------- slide 3: HLA typer development numbers
def fig_hla_dev():
    sys.path.insert(0,str(R/'hla-bench'));sys.path.insert(0,str(R/'hla-typer/stageB'))
    from names import Nomenclature,pred_slot,compare
    import score_calls as sc
    rows=tsv(R/'hla-typer/results/stageB/scores.tsv')
    genes=['A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1','DRB3','DRB4','DRB5']
    acc=collections.defaultdict(lambda:[0,0])
    for r in rows:
        if r['level']=='two_field':a=acc[r['arm'],r['gene']];a[0]+=int(r['correct']);a[1]+=1
    # T1K on the same reads, pinned IPD 3.65.0; truth from the donor's own assemblies via hla-bench names
    nom=Nomenclature();by_hap,hd=sc.load_labels(R/'hla-typer/source/haplotype_labels.tsv.gz')
    dh=collections.defaultdict(list)
    for h,(d,c) in hd.items():dh[d].append((c,h))
    donors=sorted({r['donor'] for r in rows})
    for d in donors:
        p=R/f'hla-typer/results/comparators/t1k/{d}/{d}_genotype.tsv'
        if not p.exists():continue
        calls={}
        for line in open(p):
            v=line.rstrip('\n').split('\t');calls[v[0].removeprefix('HLA-')]=v
        own=(sorted(h for c,h in dh[d] if not c.startswith('JaSaPaGe')) or sorted(h for c,h in dh[d]))[:2]
        for g in genes:
            truth=[sc.slot(by_hap[h].get(g),'two_field',nom,g,True) for h in own]
            if len(truth)!=2 or None in truth:continue
            v=calls.get(g);n=int(v[1]) if v else 0
            a1=pred_slot(v[2],'two_field',nom,g) if v and n>=1 and float(v[4])>0 else None
            a2=a1 if (v and n==1) else (pred_slot(v[5],'two_field',nom,g) if v and n==2 and float(v[7])>0 else None)
            if n==0 and g in sc.PRESENCE_GENES:a1=a2=[sc.ABSENT]
            ok=compare([a1,a2],truth)[1];a=acc['T1K',g];a[0]+=int(bool(ok));a[1]+=1
    fig,(ax,ax2)=plt.subplots(1,2,figsize=(12.5,3.3),gridspec_kw={'width_ratios':[3,1.6],'wspace':0.28})
    xs=np.arange(len(genes));w=0.2
    arms=[('hprc',HPRC,'HPRC-only panel (466)'),('asian_matched',ASIAN,'Asian size-matched panel (466)'),('full',FULL,'full panel (754)'),('T1K',T1K,'T1K, IPD-IMGT/HLA 3.65')]
    for i,(arm,col,lab) in enumerate(arms):
        v=[100*acc[arm,g][0]/acc[arm,g][1] if acc[arm,g][1] else np.nan for g in genes]
        ax.bar(xs+(i-1.5)*w,v,w,color=col,label=lab)
    ax.set_xticks(xs);ax.set_xticklabels(genes,fontsize=9);ax.set_ylim(0,128);ax.set_yticks([0,20,40,60,80,100]);ax.set_ylabel('two-field genotypes correct (%)')
    ax.axvline(7.5,color='k',lw=0.6,ls=':');ax.text(9.5,104,'presence + type',ha='center',fontsize=8)
    ax.legend(fontsize=8,ncol=2,loc='upper left',frameon=False)
    ax.set_title('Locityper, leave-one-donor-out, 40 development donors, truth = own assemblies',fontsize=9.5,loc='left')
    # pooled numbers into the caption file
    pooled={}
    for arm,_,_ in arms:
        c8=[acc[arm,g] for g in genes[:8]];d345=[acc[arm,g] for g in genes[8:]]
        pooled[arm]=(sum(a for a,_ in c8),sum(b for _,b in c8),sum(a for a,_ in d345),sum(b for _,b in d345))
    # SpecHLA step A: reconstructed gene sequences vs own assembly
    sq=tsv(R/'hla-spechla-pg/results/sequence_scores.tsv')
    g8=['A','B','C','DQA1','DQB1','DPA1','DPB1','DRB1']
    ex=collections.defaultdict(lambda:[0,0]);ed=collections.defaultdict(list)
    for r in sq:
        k=(r['mode'],r['gene']);ex[k][0]+=int(r['exact_haplotypes']);ex[k][1]+=2;ed[k]+=[int(r['edits_hap1']),int(r['edits_hap2'])]
    xs=np.arange(len(g8));w=0.38
    for i,(mode,col,lab) in enumerate([('native',LIN,'SpecHLA 1.0.12, IMGT-only\nread binning'),('A','FULL','+ pangenome panel in\nread binning (step A)')]):
        col=FULL if col=='FULL' else col
        v=[ex[mode,g][0] for g in g8]
        ax2.bar(xs+(i-0.5)*w,v,w,color=col,label=lab)
    ax2.set_xticks(xs);ax2.set_xticklabels(g8,fontsize=8,rotation=45);ax2.set_ylim(0,20);ax2.set_yticks([0,4,8,12,16]);ax2.set_ylabel('haplotypes reconstructed exactly (of 16)')
    ax2.legend(fontsize=7.5,loc='upper left',frameon=False)
    ax2.set_title('SpecHLA sequences vs own assembly, 8 donors',fontsize=9.5,loc='left')
    fig.savefig(OUT/'fig_hla_dev.png',dpi=200,bbox_inches='tight');plt.close(fig)
    with open(OUT/'hla_dev_numbers.txt','w') as f:
        for arm,(a,b,c,d) in pooled.items():f.write(f'{arm}\tclassical8 {a}/{b} = {100*a/b:.1f}%\tDRB345 {c}/{d} = {100*c/d:.1f}%\n')
        for mode in ('native','A'):
            tot=sum(ex[mode,g][0] for g in g8);n=sum(ex[mode,g][1] for g in g8);m=statistics.mean(e for g in g8 for e in ed[mode,g])
            f.write(f'spechla {mode}\texact {tot}/{n}\tmean edits {m:.1f}\n')
    print(open(OUT/'hla_dev_numbers.txt').read())

if __name__=='__main__':
    fig_cohorts();fig_graph_density();fig_variant_classes();fig_hla_dev()
