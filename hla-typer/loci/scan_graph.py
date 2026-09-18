#!/usr/bin/env python3
"""One pass over MHC.full.gfa.gz: per-node haplotype support along the GRCh38 reference walk.

For every node on GRCh38#0#ctg1 report how many haplotypes traverse it, the total
number of traversals and how many are reverse. Candidate locus boundaries are nodes
traversed exactly once, forward, by (nearly) every haplotype.
"""
import argparse,gzip,json,re
from pathlib import Path
STEP=re.compile(r'([><])([^><]+)')

def main():
    p=argparse.ArgumentParser();p.add_argument('--gfa',default='/home/asianhla/data/upload/HLA/mhc_graph/MHC.full.gfa.gz')
    p.add_argument('--ref',default='GRCh38#0#ctg1');p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    index={};length=[]
    haps=[];last=None;seen_haps=set();nonconsecutive=[]
    ref=None;trav=None
    with gzip.open(a.gfa,'rt') as f:
        for line in f:
            if line[0]=='S':
                assert trav is None,'S line after path lines'
                x=line.split('\t',3);index[x[1]]=len(length);length.append(len(x[2].rstrip('\n')))
            elif line[0] in 'WP':
                if trav is None:
                    n=len(length);trav=[0]*n;rev=[0]*n;hapc=[0]*n;lasthap=[-1]*n
                x=line.rstrip('\n').split('\t')
                if line[0]=='W':
                    hap='#'.join(x[1:3]);name='#'.join(x[1:4]);steps=STEP.findall(x[6])
                else:
                    name=x[1];hap='#'.join(name.split('#')[:2]);steps=[('>' if s[-1]=='+' else '<',s[:-1]) for s in x[2].split(',')]
                if hap!=last:
                    if hap in seen_haps:nonconsecutive.append(hap)
                    seen_haps.add(hap);haps.append(hap);last=hap
                h=len(haps)-1
                for o,node in steps:
                    i=index[node];trav[i]+=1
                    if o=='<':rev[i]+=1
                    if lasthap[i]!=h:lasthap[i]=h;hapc[i]+=1
                if name==a.ref:ref=[(o,node) for o,node in steps]
    assert ref is not None,'reference walk not found'
    ids={v:k for k,v in index.items()}
    with gzip.open(a.out/'ref_nodes.tsv.gz','wt') as out:
        out.write('ref_start0\tref_end0\tnode\tlength\tref_orientation\thaplotypes\ttraversals\treverse_traversals\n')
        pos=0
        for o,node in ref:
            i=index[node];L=length[i]
            out.write(f'{pos}\t{pos+L}\t{node}\t{L}\t{o}\t{hapc[i]}\t{trav[i]}\t{rev[i]}\n');pos+=L
    summary=dict(gfa=a.gfa,ref=a.ref,nodes=len(length),haplotypes=len(set(haps)),reference_length=pos,
                 nonconsecutive_haplotype_lines=sorted(set(nonconsecutive)))
    (a.out/'scan_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (a.out/'haplotypes.txt').write_text('\n'.join(sorted(set(haps)))+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
