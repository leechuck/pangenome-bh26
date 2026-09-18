#!/usr/bin/env python3
"""Locityper --reg-weights (per-haplotype subregion weights) and CDS intervals inside each locus sequence.

Exons of catalogued HLA genes weigh EXON, the rest of each gene span INTRON, everything else
INTERGENIC (defaults = Locityper paper weighted mode: 1, 0.1, 0.005). Exon coordinates come
from the Immuannot GTF of each haplotype (MHC-segment frame); the gene body (catalogue gene
FASTA minus 2 kb flanks, genomic orientation) is located inside the locus sequence to shift them.
A gene whose body is not found exactly once (e.g. DRB1 3' tail clipped by the locus boundary)
is placed by its longest found prefix/suffix and flagged in the audit.
"""
import argparse,csv,gzip,re,collections,json
from pathlib import Path
COMP=str.maketrans('ACGT','TGCA');FLANK=2000
LOCUS_GENES={'HLA-A':['A'],'HLA-B':['B'],'HLA-C':['C'],'HLA-DRB1':['DRB1'],'HLA-DQ':['DQA1','DQB1'],
             'HLA-DP':['DPA1','DPB1'],'HLA-DRB345':['DRB3','DRB4','DRB5']}

def fasta(path):
    out={};name=None
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt') as f:
        for line in f:
            if line[0]=='>':name=line[1:].split()[0];out[name]=[]
            else:out[name].append(line.strip().upper())
    return {k:''.join(v) for k,v in out.items()}

def exons(gtf,gene,start,end,features=('exon',)):
    """Feature intervals (1-based closed, segment frame) of HLA-<gene> within start..end."""
    out=[]
    with gzip.open(gtf,'rt') as f:
        for line in f:
            x=line.split('\t')
            if len(x)<9 or x[2] not in features or f'gene_name "HLA-{gene}"' not in x[8]:continue
            a,b=int(x[3]),int(x[4])
            if a>=start and b<=end:out.append((a,b))
    return sorted(set(out))

def merge(iv):
    out=[]
    for a,b in sorted(iv):
        if out and a<=out[-1][1]+1:out[-1][1]=max(out[-1][1],b)
        else:out.append([a,b])
    return [tuple(x) for x in out]

def tile(length,genes,w):
    """genes: list of (gene_start0, gene_end0, [(exon_start0, exon_end0)]) in locus coords -> BED rows."""
    marks=[(0,length,w['intergenic'])]
    for gs,ge,ex in genes:marks.append((gs,ge,w['intron']))
    for gs,ge,ex in genes:
        for a,b in ex:marks.append((a,b,w['exon']))
    # paint in order (later = higher priority), then run-length encode
    cuts=sorted({0,length}|{max(0,min(length,p)) for m in marks for p in m[:2]})
    rows=[]
    for a,b in zip(cuts,cuts[1:]):
        v=w['intergenic']
        for s,e,val in marks:
            if s<=a and b<=e:v=val
        if rows and rows[-1][2]==v and rows[-1][1]==a:rows[-1][1]=b
        else:rows.append([a,b,v])
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--loci',type=Path,default='results/loci');p.add_argument('--genes',type=Path,default='source/genes')
    p.add_argument('--gtf',type=Path,default='/home/leechuck/hla/gtf_all');p.add_argument('--donors',type=Path,default='source/haplotype_donors.tsv')
    p.add_argument('--out',type=Path,default='weights');p.add_argument('--exon',type=float,default=1);p.add_argument('--intron',type=float,default=0.1)
    p.add_argument('--intergenic',type=float,default=0.005);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    w=dict(exon=a.exon,intron=a.intron,intergenic=a.intergenic)
    name={r['hap_id']:r['haplotype'] for r in csv.DictReader(open(a.donors),delimiter='\t')}
    audit=collections.Counter();index=[];cds_rows=[]
    for locus,genes in LOCUS_GENES.items():
        seqs=fasta(a.loci/f'{locus}.fa.gz');placed=collections.defaultdict(list)
        for g in genes:
            recs=fasta(a.genes/f'HLA-{g}.fa')
            regions={r['name']:r for r in csv.DictReader(open(a.genes/f'HLA-{g}.regions.tsv'),delimiter='\t')}
            for rec,s in recs.items():
                hap='#'.join(rec.split('#')[:2]);L=seqs.get(hap)
                if not L or len(s)<=2*FLANK:audit[locus,g,'no_locus_or_short']+=1;continue
                body=s[FLANK:-FLANK];rg=regions[rec]
                if rg['strand']=='-':body=body.translate(COMP)[::-1]
                seg_start=int(rg['region'].rsplit(':',1)[1].split('-')[0])+FLANK   # 1-based gene start in segment
                off=L.find(body);clip=0
                if off<0 or L.find(body,off+1)>=0:
                    # boundary clipped the gene: place by the longest unique prefix or suffix
                    for k in range(1,200):
                        if (i:=L.find(body[k:]))>=0 and L.find(body[k:],i+1)<0:off=i-k;clip=k;break
                        if (i:=L.find(body[:-k]))>=0 and L.find(body[:-k],i+1)<0:off=i;clip=-k;break
                    else:audit[locus,g,'not_placed']+=1;continue
                    audit[locus,g,'clipped_at_boundary']+=1
                gtf=a.gtf/f"{hap.replace('#','_')}.gtf.gz";end1=seg_start+len(body)-1
                ex=[(off+x-seg_start,off+y-seg_start+1) for x,y in exons(gtf,g,seg_start,end1)]
                if not ex:audit[locus,g,'no_exons']+=1;continue
                cds=[(off+x-seg_start,off+y-seg_start+1) for x,y in merge(exons(gtf,g,seg_start,end1,('CDS','stop_codon')))]
                cds_rows.append(dict(locus=locus,haplotype=name[hap],hap_id=hap,gene=g,strand=rg['strand'],clipped=clip,
                                     gene_start0=off,gene_end0=off+len(body),cds=','.join(f'{x}-{y}' for x,y in cds)))
                audit[locus,g,'ok']+=1
                placed[hap].append((max(0,off),min(len(L),off+len(body)),[(max(0,x),min(len(L),y)) for x,y in ex]))
        bed=a.out/f'{locus}.bed'
        with open(bed,'w') as f:
            for hap,L in sorted(seqs.items()):
                for s,e,v in tile(len(L),placed.get(hap,[]),w):f.write(f'{name[hap]}\t{s}\t{e}\t{v}\n')
        index.append((locus,bed.resolve()))
    with open(a.out/'cds_intervals.tsv','w') as f:
        wr=csv.DictWriter(f,fieldnames=list(cds_rows[0]),delimiter='\t',lineterminator='\n');wr.writeheader();wr.writerows(cds_rows)
    with open(a.out/'reg_weights.tsv','w') as f:
        for locus,bed in index:f.write(f'{locus}\t{bed}\n')
    (a.out/'audit.json').write_text(json.dumps(dict(weights=w,counts={'|'.join(k):v for k,v in sorted(audit.items())}),indent=1)+'\n')
    for k,v in sorted(audit.items()):print(*k,v)
if __name__=='__main__':main()
