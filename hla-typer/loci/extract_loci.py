#!/usr/bin/env python3
"""Cut every haplotype walk of MHC.full.gfa.gz between the boundary nodes of each locus.

For each locus and haplotype: status ok when both boundary nodes occur exactly once on the
same contig walk in the right order; the sequence runs from the start of the left node to
the end of the right node. Each annotated gene copy (hla-analysis/source/HLA-<gene>.fa minus
its 2 kb flanks, in genomic orientation) is searched in the locus sequence as a check that
the locus contains the gene annotated on that haplotype.
"""
import argparse,csv,gzip,re,collections
from pathlib import Path
STEP=re.compile(r'([><])([^><]+)')
COMP=str.maketrans('ACGTNacgtn','TGCANtgcan')
FLANK=2000

def fasta(path):
    out={};name=None
    for line in open(path):
        if line[0]=='>':name=line[1:].split()[0];out[name]=[]
        else:out[name].append(line.strip().upper())
    return {k:''.join(v) for k,v in out.items()}

def main():
    p=argparse.ArgumentParser();p.add_argument('--gfa',default='/home/asianhla/data/upload/HLA/mhc_graph/MHC.full.gfa.gz')
    p.add_argument('--spec',type=Path,default='source/loci_spec.tsv');p.add_argument('--genes',type=Path,default='source/genes')
    p.add_argument('--out',type=Path,default='results/loci');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    spec=list(csv.DictReader(open(a.spec),delimiter='\t'))
    seq={};haps=set()
    nL=collections.Counter();nR=collections.Counter();cut={};same_line=set()
    with gzip.open(a.gfa,'rt') as f:
        for line in f:
            if line[0]=='S':
                x=line.rstrip('\n').split('\t');seq[x[1]]=x[2].upper()
            elif line[0]=='W':
                x=line.rstrip('\n').split('\t');hap='#'.join(x[1:3]);haps.add(hap)
                steps=STEP.findall(x[6]);pos=collections.defaultdict(list)
                for i,(o,n) in enumerate(steps):pos[n].append(i)
                for s in spec:
                    k=(s['locus'],hap);li=pos.get(s['left_node'],[]);ri=pos.get(s['right_node'],[])
                    nL[k]+=len(li);nR[k]+=len(ri)
                    if len(li)==1 and len(ri)==1:
                        same_line.add(k);i,j=li[0],ri[0]
                        if i<=j and steps[i][0]=='>' and steps[j][0]=='>':
                            cut[k]=''.join(seq[n] if o=='>' else seq[n].translate(COMP)[::-1] for o,n in steps[i:j+1])
    rows=[];seqs={}
    for s in spec:
        with gzip.open(a.out/f"{s['locus']}.fa.gz",'wt') as fa:
            for hap in sorted(haps):
                k=(s['locus'],hap);out=''
                if nL[k]==0 and nR[k]==0:status='no_boundary'
                elif nL[k]!=1 or nR[k]!=1:status='boundary_missing_or_repeated'
                elif k not in same_line:status='boundaries_on_different_contigs'
                elif k not in cut:status='boundary_order_or_orientation'
                else:status='ok';out=cut[k];seqs[k]=out;fa.write(f'>{hap}\n{out}\n')
                rows.append(dict(locus=s['locus'],hap_id=hap,status=status,length=len(out),N_bases=out.count('N'),non_ACGT=sum(c not in 'ACGT' for c in out)))
    # gene containment check
    genes={g for s in spec for g in s['genes'].split(',')}
    found_genes=collections.defaultdict(list)
    for g in sorted(genes):
        recs=fasta(a.genes/f'HLA-{g}.fa')
        regions={r['name']:r for r in csv.DictReader(open(a.genes/f'HLA-{g}.regions.tsv'),delimiter='\t')}
        locus=next(s['locus'] for s in spec if g in s['genes'].split(','))
        for name,body in recs.items():
            hap='#'.join(name.split('#')[:2])
            if len(body)<=2*FLANK:continue
            body=body[FLANK:-FLANK]
            if regions[name]['strand']=='-':body=body.translate(COMP)[::-1]
            ls=seqs.get((locus,hap),'')
            found_genes[locus,hap].append(f"{name.split('#')[2]}:{int(bool(ls) and body in ls)}")
    for r in rows:r['genes_contained']=';'.join(found_genes.get((r['locus'],r['hap_id']),[]))
    with open(a.out/'locus_haplotypes.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    summary=collections.Counter((r['locus'],r['status']) for r in rows)
    for k,v in sorted(summary.items()):print(*k,v)
    for s in spec:
        rs=[r for r in rows if r['locus']==s['locus'] and r['status']=='ok']
        tot=sum(r['genes_contained'].count(':') for r in rs);ok=sum(r['genes_contained'].count(':1') for r in rs)
        L=sorted(r['length'] for r in rs)
        print(s['locus'],'ok',len(rs),'length min/median/max',L[0],L[len(L)//2],L[-1],'acgt_only',sum(r['non_ACGT']==0 for r in rs),'gene copies contained',ok,'/',tot)
if __name__=='__main__':main()
