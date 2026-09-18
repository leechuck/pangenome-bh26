#!/usr/bin/env python3
"""Score Stage B (panel haplotype-pair) calls against the donor's own assembly labels.

Prediction for a gene = the catalogue labels of the gene copy carried by each selected panel
haplotype. Truth = the donor's own two haplotypes (primary assembly, not the JaSaPaGe re-assembly).
Levels: two_field and g_group (hla-bench/names.py, IPD 3.65.0), cds_exact and gene_exact
(sequence identity via catalogue sha256). DRB3/4/5: copy number per haplotype (0/1) is part of
the genotype; a haplotype without the gene contributes the value ABSENT. No-calls are errors.
"""
import argparse,csv,gzip,json,sys,collections
from pathlib import Path
H=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(H.parent/'hla-bench'))
from names import Nomenclature,truth_slot,pred_slot,compare
csv.field_size_limit(10**9)
LOCUS_GENES={'HLA-A':['A'],'HLA-B':['B'],'HLA-C':['C'],'HLA-DRB1':['DRB1'],'HLA-DQ':['DQA1','DQB1'],
             'HLA-DP':['DPA1','DPB1'],'HLA-DRB345':['DRB3','DRB4','DRB5']}
PRESENCE_GENES={'DRB3','DRB4','DRB5'}
ABSENT=frozenset(['ABSENT'])

def load_labels(path):
    by_hap=collections.defaultdict(dict);hap_donor={}
    with gzip.open(path,'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            by_hap[r['haplotype']][r['gene'].removeprefix('HLA-')]=r;hap_donor[r['haplotype']]=(r['donor_id'],r['cohort'])
    return by_hap,hap_donor

def slot(row,level,nom,gene,truth):
    """One haplotype's value for a gene at a level; None if not determinable."""
    if row is None:return ABSENT if gene in PRESENCE_GENES else None
    if level in ('cds_exact','gene_exact'):
        v=row['cds_sha256' if level=='cds_exact' else 'gene_sha256']
        ok=row['cds_complete']=='1' or level=='gene_exact'
        return (frozenset([v]) if truth else [frozenset([v])]) if v and ok else None
    names=row['exact_cds_alleles'] if row['cds_complete']=='1' and row['exact_cds_alleles'] else row['exact_protein_alleles']
    if not names:return None
    if truth:return truth_slot('/'.join(n.split('*',1)[1] for n in names.split(';')),level,nom,gene)
    return pred_slot(','.join(names.split(';')),level,nom,gene)

def candidate_rows(labels_path,by_hap,ipd_dir=H/'source/ipd_cds'):
    """Stage C1 candidates -> {candidate: {gene: catalogue-like row}}. Grafted genes take the IPD names of their CDS."""
    ipd={}
    for f in ipd_dir.glob('*.tsv'):
        for r in csv.DictReader(open(f),delimiter='\t'):ipd[r['cds_sha256']]=r['alleles']
    out=collections.defaultdict(dict)
    for r in csv.DictReader(open(labels_path),delimiter='\t'):
        base=by_hap.get(r['backbone'],{})
        out[r['candidate']].update({g:v for g,v in base.items() if g not in out[r['candidate']]})
        if r['source']!='backbone':
            names=r.get('alleles') or ipd[r['cds_sha256']]
            out[r['candidate']][r['gene']]=dict(cds_complete='1',exact_cds_alleles=';'.join(names.split(';')),exact_protein_alleles='',
                                                  cds_sha256=r['cds_sha256'],gene_sha256='',source=r['source'])
        elif r['allele']=='ABSENT':out[r['candidate']].pop(r['gene'],None)
    return out

def pick_sources(picks,lookup,gene):
    """'backbone', 't1k', 'ipd_near', 'ipd_partial', 't1k_joint' per picked haplotype for this gene."""
    return ','.join(lookup.get(x,{}).get(gene,{}).get('source','backbone') for x in picks)

def main():
    p=argparse.ArgumentParser();p.add_argument('--calls',type=Path,default=H/'results/calls')
    p.add_argument('--c1',action='store_true',help='CALLS is a c1 root: <arm-variant>/<donor>/calls/loci/<locus>/res.json.gz')
    p.add_argument('--labels',type=Path,default=H/'source/haplotype_labels.tsv.gz');p.add_argument('--out',type=Path,default=H/'results/stageB')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    nom=Nomenclature();by_hap,hap_donor=load_labels(a.labels)
    donor_haps=collections.defaultdict(list)
    for h,(d,c) in hap_donor.items():donor_haps[d].append((c,h))
    rows=[]
    pattern='*/*/calls/loci/*/res.json.gz' if a.c1 else '*/*/loci/*/res.json.gz'
    cand_cache={}
    for res in sorted(a.calls.glob(pattern)):
        locus=res.parent.name
        if a.c1:
            donor=res.parents[3].name;arm=res.parents[4].name
            key=res.parents[3]/'candidates/labels.tsv'
            if key not in cand_cache:cand_cache[key]=candidate_rows(key,by_hap)
            lookup=cand_cache[key]
        else:
            donor=res.parent.parent.parent.name;arm=res.parent.parent.parent.parent.name;lookup=by_hap
        r=json.load(gzip.open(res,'rt'));gt=r.get('genotype') or ''
        picks=gt.split(',') if isinstance(gt,str) and gt else []
        own=sorted(donor_haps[donor],key=lambda x:x[0].startswith('JaSaPaGe'))
        own=[h for c,h in own if not c.startswith('JaSaPaGe')] or [h for c,h in own]
        own=sorted(own)[:2]
        for gene in LOCUS_GENES[locus]:
            for level in ['two_field','g_group','cds_exact','gene_exact']:
                truth=[slot(by_hap[h].get(gene),level,nom,gene,True) for h in own]
                if len(truth)!=2 or None in truth:continue
                pair=[slot(lookup.get(x,{}).get(gene),level,nom,gene,False) for x in picks] if len(picks)==2 else [None,None]
                pair=[[p] if isinstance(p,frozenset) else p for p in pair]
                called,correct,matches=compare(pair,truth)
                rows.append(dict(donor=donor,arm=arm,locus=locus,gene=gene,level=level,picks=gt,quality=round(r.get('quality',0),2),
                                 truth_haplotypes=','.join(own),called=called,correct=correct,allele_matches=matches,
                                 sources=pick_sources(picks,lookup,gene) if a.c1 else 'backbone'))
    if not rows:sys.exit('no calls found')
    with open(a.out/'scores.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    agg=collections.defaultdict(lambda:[0,0,0])
    for r in rows:
        for g in [r['gene'],'ALL']:
            x=agg[r['arm'],r['level'],g];x[0]+=1;x[1]+=r['correct'];x[2]+=r['called']
    with open(a.out/'summary.tsv','w') as f:
        f.write('arm\tlevel\tgene\tgenotypes\tcorrect\tcalled\taccuracy_pct\n')
        for (arm,level,g),(n,c,k) in sorted(agg.items()):f.write(f'{arm}\t{level}\t{g}\t{n}\t{c}\t{k}\t{100*c/n:.2f}\n')
    for (arm,level,g),(n,c,k) in sorted(agg.items()):
        if g=='ALL':print(arm,level,f'{c}/{n}',f'{100*c/n:.1f}%')
if __name__=='__main__':main()
