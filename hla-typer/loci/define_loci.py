#!/usr/bin/env python3
"""Choose locus boundaries as graph nodes shared by (nearly) all haplotypes.

A boundary node lies on GRCh38#0#ctg1, is >= MIN_LEN bp, forward in every traversal,
and traversed exactly once by each haplotype that has it (traversals == haplotypes).
In a search window, nodes carried by >= MIN_FRACTION of haplotypes are preferred, then
the node closest to the target flank distance. Loci are defined in ctg1 coordinates
(0-based, half-open), with the boundary nodes included. HLA-DRB345 and HLA-DRB1 share
the DRB1 3' boundary node (one node of overlap).
"""
import argparse,csv,gzip,json
from pathlib import Path
R=Path(__file__).resolve().parent.parent
MIN_LEN=32
MIN_FRACTION=0.99
# Gene bodies on ctg1 (0-based half-open) = hla-analysis/source/HLA-*.regions.tsv GRCh38 rows without the 2 kb flanks.
BODY={'A':(1531553,1535055),'C':(2857553,2861871),'B':(2942661,2946742),'DRB1':(4168068,4179148),
      'DQA1':(4226016,4232500),'DQB1':(4249330,4256432),'DPA1':(4653618,4663393),'DPB1':(4664975,4676501)}
DRA_END=4075000   # universal anchors end here; the DR structural block follows (see gap analysis)
FLANK=10000
# locus: (left search window, left target, right search window, right target, genes)
SPEC={
 'HLA-A':   ((BODY['A'][0]-25000,BODY['A'][0]-3000),BODY['A'][0]-FLANK,(BODY['A'][1]+3000,BODY['A'][1]+25000),BODY['A'][1]+FLANK,['A']),
 'HLA-C':   ((BODY['C'][0]-25000,BODY['C'][0]-3000),BODY['C'][0]-FLANK,(BODY['C'][1]+3000,BODY['C'][1]+25000),BODY['C'][1]+FLANK,['C']),
 'HLA-B':   ((BODY['B'][0]-25000,BODY['B'][0]-3000),BODY['B'][0]-FLANK,(BODY['B'][1]+3000,BODY['B'][1]+25000),BODY['B'][1]+FLANK,['B']),
 'HLA-DRB345':((DRA_END-8000,DRA_END+1000),DRA_END,None,None,['DRB3','DRB4','DRB5']),
 'HLA-DRB1':((BODY['DRB1'][0]-1000,BODY['DRB1'][0]+1000),BODY['DRB1'][0],(BODY['DRB1'][1]+3000,BODY['DRB1'][1]+25000),BODY['DRB1'][1]+FLANK,['DRB1']),
 'HLA-DQ':  ((BODY['DQA1'][0]-25000,BODY['DQA1'][0]-3000),BODY['DQA1'][0]-FLANK,(BODY['DQB1'][1]+3000,BODY['DQB1'][1]+25000),BODY['DQB1'][1]+FLANK,['DQA1','DQB1']),
 'HLA-DP':  ((BODY['DPA1'][0]-25000,BODY['DPA1'][0]-3000),BODY['DPA1'][0]-FLANK,(BODY['DPB1'][1]+3000,BODY['DPB1'][1]+25000),BODY['DPB1'][1]+FLANK,['DPA1','DPB1']),
}

def pick(nodes,window,target,side):
    lo,hi=window
    c=[n for n in nodes if lo<=n['ref_start0'] and n['ref_end0']<=hi]
    if not c:return None
    edge=(lambda n:n['ref_end0']) if side=='left' else (lambda n:n['ref_start0'])
    # Near-universal nodes (>= MIN_FRACTION) are equivalent; then prefer the target flank distance.
    return max(c,key=lambda n:(n['haplotypes']>=MIN_FRACTION*754,-abs(edge(n)-target),n['haplotypes'],n['length']))

def main():
    p=argparse.ArgumentParser();p.add_argument('--scan',type=Path,default=R/'results/graph_scan/ref_nodes.tsv.gz')
    p.add_argument('--out',type=Path,default=R/'source/loci_spec.tsv');a=p.parse_args()
    nodes=[]
    for r in csv.DictReader(gzip.open(a.scan,'rt'),delimiter='\t'):
        n={k:(int(v) if k not in ('node','ref_orientation') else v) for k,v in r.items()}
        if n['length']>=MIN_LEN and n['ref_orientation']=='>' and n['reverse_traversals']==0 and n['traversals']==n['haplotypes']:nodes.append(n)
    total=754;out=[];chosen={}
    for locus,(lw,lt,rw,rt,genes) in SPEC.items():
        left=pick(nodes,lw,lt,'left');chosen[locus]=[left,None]
    for locus,(lw,lt,rw,rt,genes) in SPEC.items():
        # DRB345 ends with the DRB1 left boundary node, so the two loci tile the DR block.
        right=chosen['HLA-DRB1'][0] if rw is None else pick(nodes,rw,rt,'right')
        left=chosen[locus][0]
        assert left and right,locus
        start,end=left['ref_start0'],right['ref_end0']
        out.append(dict(locus=locus,genes=','.join(genes),left_node=left['node'],left_ref_start0=left['ref_start0'],left_length=left['length'],left_haplotype_fraction=round(left['haplotypes']/total,4),
                        right_node=right['node'],right_ref_end0=right['ref_end0'],right_length=right['length'],right_haplotype_fraction=round(right['haplotypes']/total,4),
                        ref_start0=start,ref_end0=end,ref_length=end-start,chr6_start1=start+28410701,chr6_end1=end+28410700))
    with open(a.out,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
    for r in out:print(r['locus'],r['ref_start0'],r['ref_end0'],r['ref_length'],r['left_haplotype_fraction'],r['right_haplotype_fraction'],r['left_length'],r['right_length'])
if __name__=='__main__':main()
