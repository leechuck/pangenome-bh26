#!/usr/bin/env python3
"""Dense diagnostic probes from training C4 genes; no held-out marker discovery."""
import csv,gzip,json,re
from collections import defaultdict,Counter
from pathlib import Path
from Bio import SeqIO
import numpy as np
ROOT=Path(__file__).resolve().parent
GROUPS={'A':1,'B':2,'LONG_LEFT':4,'LONG_RIGHT':8,'SHORT':16,'BOUNDARY':32,'NONCANONICAL':64}
COMP=str.maketrans('ACGT','TGCA')
def key(s):
    s=min(s,s.translate(COMP)[::-1]);v=0
    for c in s:v=(v<<2)|'ACGT'.index(c)
    return v
def write(p,rows):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def main():
    meta=json.loads((ROOT/'source/target_definitions.json').read_text());sites=meta['diagnostic_sites0'];left=meta['HERV_start0'];right=meta['HERV_end0']
    genes={r.id:str(r.seq).upper() for r in SeqIO.parse(ROOT/'source/training_C4_genes.fa','fasta')}
    labels={r['gene_id']:r for r in csv.DictReader(open(ROOT/'source/training_C4_genes.tsv'),delimiter='\t')}
    align={}
    for line in open(ROOT/'source/training_C4_genes.paf'):
        x=line.rstrip().split('\t')
        if x[4]!='+' or (int(x[3])-int(x[2]))/int(x[1])<.95:continue
        if x[0] not in align or int(x[10])>int(align[x[0]][10]):align[x[0]]=x
    probes=defaultdict(set);records=[]
    def add(s,starts,group):
        for start in starts:
            if start<0 or start+31>len(s):continue
            k=s[start:start+31]
            if set(k)<=set('ACGT'):probes[key(k)].add(group)
    for name,s in genes.items():
        if name not in align:
            records.append(dict(gene_id=name,annotated_form=labels[name]['annotated_form'],motif='',inferred_form='',status='alignment_incomplete'));continue
        x=align[name];tags={v.split(':',2)[0]:v.split(':',2)[2] for v in x[12:]};q=int(x[2]);t=int(x[7]);mapping=np.full(21058,-1,dtype=int)
        for length,op in re.findall(r'(\d+)([MIDNSHP=X])',tags['cg']):
            n=int(length)
            if op in 'M=X':mapping[t:t+n]=np.arange(q,q+n);q+=n;t+=n
            elif op in 'DN':t+=n
            elif op in 'IS':q+=n
        mapped=[int(mapping[i]) for i in sites]
        motif=''.join(s[i] for i in mapped) if min(mapped)>=0 else ''
        ab='A' if motif==meta['A_motif'] else 'B' if motif==meta['B_motif'] else '?'
        qleft,qright=int(mapping[left-1]),int(mapping[right]);span=qright-qleft-1 if qleft>=0 and qright>=0 else -1
        ls='L' if span>=6000 else 'S' if 0<=span<100 else '?'
        form='C4'+ab+ls;status='match' if form==labels[name]['annotated_form'] else 'diagnostic_annotation_disagreement'
        records.append(dict(gene_id=name,annotated_form=labels[name]['annotated_form'],motif=motif,inferred_form=form,status=status))
        if ab!='?' and mapped[-1]-mapped[0]<=30:add(s,range(mapped[-1]-30,mapped[0]+1),ab)
        elif motif and mapped[-1]-mapped[0]<=30:add(s,range(mapped[-1]-30,mapped[0]+1),'NONCANONICAL')
        if ls=='L':
            add(s,range(qleft+1-23,qleft+1-7),'LONG_LEFT');add(s,range(qright-23,qright-7),'LONG_RIGHT')
        elif ls=='S':add(s,range(qleft+1-23,qleft+1-7),'SHORT')
    for r in SeqIO.parse(ROOT/'source/module_boundary_contexts.fa','fasta'):add(str(r.seq),range(0,len(r.seq)-30,8),'BOUNDARY')
    # A probe cannot be evidence for mutually exclusive biological states.
    conflicts={k for k,g in probes.items() if {'A','B'}<=g or ('SHORT' in g and bool(g&{'LONG_LEFT','LONG_RIGHT'}))}
    for k in conflicts:del probes[k]
    keys=set(int(x) for x in (ROOT.parent/'hla-structural/source/markers.txt').read_text().splitlines())|set(probes)
    keys=sorted(keys)
    (ROOT/'source/markers.txt').write_text(''.join(str(k)+'\n' for k in keys))
    (ROOT/'source/marker_groups.tsv').write_text(''.join(str(k)+'\t'+str(sum(GROUPS[g] for g in probes.get(k,set())))+'\n' for k in keys))
    write(ROOT/'results/training_C4_diagnostic_audit.tsv',records)
    summary=dict(total_markers=len(keys),targeted_markers=len(probes),group_memberships=dict(Counter(g for groups in probes.values() for g in groups)),conflicting_probes_removed=len(conflicts),training_gene_audit=dict(Counter(r['status'] for r in records)),groups=GROUPS)
    (ROOT/'source/marker_design.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
