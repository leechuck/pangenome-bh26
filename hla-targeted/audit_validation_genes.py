#!/usr/bin/env python3
"""Post-freeze audit of validation assembly labels, never used for marker discovery."""
import argparse,csv,gzip,json,re
from pathlib import Path
import numpy as np
from Bio import SeqIO
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rs):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
def prepare():
    assert (ROOT/'source/validation_freeze.json').exists(),'Freeze inference first'
    donors={r['donor'] for r in read(ROOT/'source/validation_donors.tsv')}
    panel=[r for r in read(ROOT.parent/'hla-structural/source/catalogue.tsv') if r['locus']=='RCCX' and r['donor'] in donors]
    assert len(panel)==2*len(donors)
    with gzip.open(ROOT.parent/'hla-structural/source/locus_sequences.fa.gz','rt') as f:seq={r.id:str(r.seq) for r in SeqIO.parse(f,'fasta') if r.id in {x['sequence_id'] for x in panel}}
    genes=[]
    with open(ROOT/'source/validation_C4_genes.fa','w') as out:
        for r in panel:
            n=0;entry=r['hap_id'].replace('#','_');offset=int(r['start0'])
            with gzip.open(ROOT.parent/f'hla-analysis/source/gtf_all/{entry}.gtf.gz','rt') as f:
                for line in f:
                    if line.startswith('#'):continue
                    x=line.rstrip().split('\t')
                    if len(x)<9 or x[2]!='gene':continue
                    m=re.search(r'gene_name "(C4[AB][LS])"',x[8])
                    if not m:continue
                    assert x[0]==r['contig'];a=int(x[3])-1;b=int(x[4]);s=seq[r['sequence_id']][a-offset:b-offset];assert len(s)==b-a
                    if x[6]=='-':s=s.translate(str.maketrans('ACGT','TGCA'))[::-1]
                    n+=1;name=entry+'__C4_'+str(n);out.write('>'+name+'\n'+s+'\n')
                    genes.append(dict(gene_id=name,hap_id=r['hap_id'],donor=r['donor'],annotated_form=m[1]))
    write(ROOT/'source/validation_C4_genes.tsv',genes)
    print('Prepared',len(genes),'validation C4 sequences for label audit')
def audit():
    assert (ROOT/'source/validation_freeze.json').exists()
    meta=json.loads((ROOT/'source/target_definitions.json').read_text());genes={r.id:str(r.seq).upper() for r in SeqIO.parse(ROOT/'source/validation_C4_genes.fa','fasta')};align={}
    for line in open(ROOT/'source/validation_C4_genes.paf'):
        x=line.rstrip().split('\t')
        if x[4]!='+' or (int(x[3])-int(x[2]))/int(x[1])<.95:continue
        if x[0] not in align or int(x[10])>int(align[x[0]][10]):align[x[0]]=x
    result=[]
    for r in read(ROOT/'source/validation_C4_genes.tsv'):
        name=r['gene_id'];s=genes[name]
        if name not in align:result.append(dict(**r,motif='',inferred_form='',status='alignment_incomplete'));continue
        x=align[name];tags={v.split(':',2)[0]:v.split(':',2)[2] for v in x[12:]};q=int(x[2]);t=int(x[7]);mapping=np.full(21058,-1,int)
        for length,op in re.findall(r'(\d+)([MIDNSHP=X])',tags['cg']):
            n=int(length)
            if op in 'M=X':mapping[t:t+n]=np.arange(q,q+n);q+=n;t+=n
            elif op in 'DN':t+=n
            elif op in 'IS':q+=n
        mapped=[int(mapping[i]) for i in meta['diagnostic_sites0']];motif=''.join(s[i] for i in mapped) if min(mapped)>=0 else ''
        ab='A' if motif==meta['A_motif'] else 'B' if motif==meta['B_motif'] else '?'
        left,right=int(mapping[meta['HERV_start0']-1]),int(mapping[meta['HERV_end0']]);span=right-left-1 if min(left,right)>=0 else -1
        ls='L' if span>=6000 else 'S' if 0<=span<100 else '?';form='C4'+ab+ls
        result.append(dict(**r,motif=motif,inferred_form=form,status='match' if form==r['annotated_form'] else 'diagnostic_annotation_disagreement'))
    write(ROOT/'results/validation_C4_diagnostic_audit.tsv',result)
    print('Audited',len(result),'genes;',sum(r['status']!='match' for r in result),'nonmatching/incomplete')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('step',choices=['prepare','audit']);args=p.parse_args();prepare() if args.step=='prepare' else audit()
