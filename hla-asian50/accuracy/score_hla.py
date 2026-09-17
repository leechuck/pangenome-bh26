#!/usr/bin/env python3
"""Secondary numeric two-field HLA endpoint against independent experimental labels.
No graph-derived HLA calls exist in this run, so this does not rank graph panels.
"""
import csv,tarfile,re,collections
from pathlib import Path
from score_panels import write
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent/'results'
GENES=['A','B','C','DRB1','DQB1']
def norm(x):
    x=x.split('*')[-1];m=re.match(r'^(\d+):(\d+)',x)
    return ':'.join(m.groups()) if m and int(m[2]) else None

def singleton(x):
    vals={norm(a) for a in re.split('[,;/]',x)}
    return next(iter(vals)) if len(vals)==1 and None not in vals else None

def main():
    ds={r['donor']:r for r in csv.DictReader(open(ROOT/'hla-asian50/run/source/donors.tsv'),delimiter='\t')}
    truth=list(csv.DictReader(open(ROOT/'hla-pilot/source/1000G_HLA.txt'),delimiter=' '));rows=[]
    with tarfile.open(ROOT/'hla-asian50/run/source/native_HLA_calls.tar.gz') as tar:
        for t in truth:
            s=t['id']
            if s not in ds:continue
            t1k={}
            for line in tar.extractfile(f't1k/{s}/{s}_genotype.tsv').read().decode().splitlines():
                v=line.split('\t');g=v[0].removeprefix('HLA-');num=int(v[1]);p1=singleton(v[2]) if float(v[4])>0 else None
                # For the five classical diploid loci, one reported allele is interpreted as homozygous.
                # This does not prove CN=2 or distinguish biological allele dropout.
                p2=p1 if num==1 else singleton(v[5]) if num==2 and float(v[7])>0 else None
                t1k[g]=(p1,p2) if num else (None,None)
            lines=[l for l in tar.extractfile(f'spechla/{s}/hla.result.txt').read().decode().splitlines() if not l.startswith('#')]
            sp=dict(zip(lines[0].split('\t'),lines[1].split('\t')))
            for g in GENES:
                allowed=[{norm(a) for a in t[k].split('/')} for k in [g,g+'.1']]
                if any(None in vals for vals in allowed):continue
                for method,pair in [('T1K',t1k.get(g,(None,None))),('SpecHLA',tuple(singleton(sp[f'HLA_{g}_{h}']) for h in [1,2]))]:
                    called=all(pair);correct=bool(called and ((pair[0] in allowed[0] and pair[1] in allowed[1]) or (pair[0] in allowed[1] and pair[1] in allowed[0])))
                    allele_matches=max(sum(x is not None and x in vals for x,vals in zip(ps,allowed)) for ps in [pair,pair[::-1]])
                    rows.append(dict(donor=s,stratum=ds[s]['stratum'],gene=g,method=method,truth=t[g]+';'+t[g+'.1'],prediction=';'.join(x or 'NO_CALL' for x in pair),called=int(called),correct=int(correct),allele_matches=allele_matches))
    write(OUT/'hla_experimental.tsv',rows);sums=[]
    for method in ['T1K','SpecHLA']:
        for gene in GENES+['ALL']:
            rs=[r for r in rows if r['method']==method and (gene=='ALL' or r['gene']==gene)]
            sums.append(dict(method=method,gene=gene,donors=len({r['donor'] for r in rs}),n=len(rs),called=sum(r['called'] for r in rs),correct=sum(r['correct'] for r in rs),allele_n=2*len(rs),allele_matches=sum(r['allele_matches'] for r in rs)))
    write(OUT/'hla_summary.tsv',sums)
if __name__=='__main__':main()
