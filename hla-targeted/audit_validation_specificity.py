#!/usr/bin/env python3
"""Post-evaluation full-MHC check of held-out diagnostic probe multiplicity."""
import csv,json
from pathlib import Path
from collections import defaultdict
import numpy as np
from Bio import SeqIO
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def main():
    definitions=[tuple(map(int,s.split())) for s in (ROOT/'source/marker_groups.tsv').read_text().splitlines()]
    keys={k:i for i,(k,m) in enumerate(definitions) if m&95};groups=json.loads((ROOT/'source/marker_design.json').read_text())['groups'];genes={r['gene_id']:r for r in read(ROOT/'source/validation_C4_genes.tsv')}
    counts=defaultdict(lambda:np.zeros(len(definitions),dtype=np.uint32));bases={'A':0,'C':1,'G':2,'T':3}
    for gene in SeqIO.parse(ROOT/'source/validation_C4_genes.fa','fasta'):
        target=counts[genes[gene.id]['hap_id']];f=rev=n=0
        for c in str(gene.seq).upper():
            if c not in bases:f=rev=n=0;continue
            b=bases[c];f=((f<<2)|b)&((1<<62)-1);rev=(rev>>2)|((3-b)<<60);n+=1
            if n>=31:
                j=keys.get(min(f,rev))
                if j is not None:target[j]+=1
    results=[];extra=[]
    for hap,expected in sorted(counts.items()):
        actual=np.fromfile(ROOT/'source/validation_background'/f'{hap.replace("#","_")}.bin',dtype='<u4');assert len(actual)==len(definitions)
        for group,bit in groups.items():
            if group=='BOUNDARY':continue
            ids=[j for j,(_,mask) in enumerate(definitions) if mask&bit]
            delta=actual[ids].astype(np.int64)-expected[ids].astype(np.int64);assert np.all(delta>=0),(hap,group)
            results.append(dict(hap_id=hap,group=group,probes=len(ids),probes_with_extra_matches=int((delta>0).sum()),extra_matches=int(delta.sum())))
            if np.any(delta):extra.append((hap,group,int(delta.sum())))
    with open(ROOT/'results/validation_MHC_specificity.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
    summary=dict(haplotypes=len(counts),C4_genes=len(genes),diagnostic_probes=len(keys),haplotype_group_combinations_with_extra_matches=len(extra),extra_matches=extra,scope='Whole held-out MHC versus annotated C4 genes; excludes module-start contexts and does not establish whole-genome specificity',analysis_stage='post-evaluation diagnostic; frozen caller unchanged')
    (ROOT/'results/validation_MHC_specificity.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
