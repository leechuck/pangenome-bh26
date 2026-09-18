#!/usr/bin/env python3
"""Genotype priors for a Stage C1 candidate set (Locityper --priors: <locus> <hap1,hap2> <log10 prior>).

Prior of a diploid pair = f(a) f(b) (x2 if a != b), where f(h) is the product over the locus genes of the
two-field allele frequency in the leave-one-out panel (Laplace pseudo-count so grafted alleles absent from
the panel keep a floor). This counters the likelihood's preference for heterozygous pairs of near-identical
haplotypes over a homozygous common allele.
"""
import argparse,csv,gzip,math,collections,sys
from pathlib import Path
H=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(H/'stageB'))
from score_calls import load_labels,candidate_rows,LOCUS_GENES

def two_field(row):
    if not row:return 'ABSENT'
    names=row['exact_cds_alleles'] if row['cds_complete']=='1' and row['exact_cds_alleles'] else row['exact_protein_alleles']
    if not names:return None
    x=names.split(';')[0].split('*')[-1].split(':')
    return ':'.join(x[:2])

def panel_frequencies(by_hap,hap_donor,donor,genes,pseudo=0.5):
    """{gene: {two_field: freq}} over leave-one-out panel haplotypes (ABSENT counted for paralogs)."""
    out={}
    for g in genes:
        c=collections.Counter()
        for h,rows in by_hap.items():
            if hap_donor[h][0]==donor or hap_donor[h][1]=='REF':continue
            t=two_field(rows.get(g))
            if t:c[t]+=1
        n=sum(c.values());k=len(c)+1
        out[g]=(lambda c=c,n=n,k=k:{'_n':n,'_k':k,**{a:(v+pseudo)/(n+pseudo*k) for a,v in c.items()}})()
        out[g]['_floor']=pseudo/(n+pseudo*k)
    return out

def hap_prior(cand_genes,freqs,genes):
    p=1.0
    for g in genes:
        t=two_field(cand_genes.get(g)) if cand_genes.get(g) else 'ABSENT'
        f=freqs[g].get(t,freqs[g]['_floor']) if t else freqs[g]['_floor']
        p*=f
    return p

def main():
    p=argparse.ArgumentParser();p.add_argument('--donor',required=True);p.add_argument('--candidates',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--pseudo',type=float,default=0.5)
    p.add_argument('--db',type=Path,required=True,help='locityper database built from these candidates (identical sequences collapsed)');a=p.parse_args()
    by_hap,hap_donor=load_labels(H/'source/haplotype_labels.tsv.gz')
    cand=candidate_rows(a.candidates/'labels.tsv',by_hap)
    loci=collections.defaultdict(list)
    for r in csv.DictReader(open(a.candidates/'labels.tsv'),delimiter='\t'):
        if r['candidate'] not in loci[r['locus']]:loci[r['locus']].append(r['candidate'])
    n=0
    with open(a.out,'w') as f:
        for L,names in loci.items():
            genes=LOCUS_GENES[L];freqs=panel_frequencies(by_hap,hap_donor,a.donor,genes,a.pseudo)
            # identical candidate sequences are one haplotype in the database: sum their priors under the representative
            rep={}
            for line in open(a.db/'loci'/L/'discarded_haplotypes.txt'):
                r,others=line.split('=');r=r.strip()
                for o in others.split(','):
                    if o.strip():rep[o.strip()]=r
            hp=collections.defaultdict(float)
            for h in names:hp[rep.get(h,h)]+=hap_prior(cand.get(h,{}),freqs,genes)
            names=sorted(hp)
            for i,h1 in enumerate(names):
                for h2 in names[i:]:
                    pr=hp[h1]*hp[h2]*(1 if h1==h2 else 2)
                    f.write(f'{L}\t{h1},{h2}\t{math.log10(pr):.4f}\n');n+=1
    print('priors written:',n)
if __name__=='__main__':main()
