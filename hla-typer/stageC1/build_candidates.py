#!/usr/bin/env python3
"""Stage C1: open-set candidate haplotypes for one donor.

Backbones = every haplotype in the Stage B options with probability >= --min-prob (at least the
top --top-k genotypes), from the leave-one-out panel, so no donor sequence enters. Candidates per
backbone and gene: IPD-IMGT/HLA complete CDS within --max-edits of the backbone CDS, plus the
alleles T1K called for the donor (ensemble arm; --no-t1k for the ablation). Each candidate CDS is
grafted onto the backbone (graft.py). For two-gene loci (DQ, DP) the T1K alleles of both genes are
also grafted jointly. Writes a Locityper target set plus labels and region weights.
"""
import argparse,csv,gzip,json,hashlib,collections,sys
from pathlib import Path
import edlib
from graft import graft,cds_of

def tile(length,marks):
    """Paint (start,end,weight) marks in order onto [0,length) and run-length encode."""
    cuts=sorted({0,length}|{max(0,min(length,x)) for m in marks for x in m[:2]});rows=[]
    for u,v in zip(cuts,cuts[1:]):
        w=[m[2] for m in marks if m[0]<=u and v<=m[1]][-1]
        if rows and rows[-1][2]==w:rows[-1][1]=v
        else:rows.append([u,v,w])
    return rows
H=Path(__file__).resolve().parent.parent
def safe(allele):return allele.replace('*','_').replace(':','-')
LOCUS_GENES={'HLA-A':['A'],'HLA-B':['B'],'HLA-C':['C'],'HLA-DRB1':['DRB1'],'HLA-DQ':['DQA1','DQB1'],
             'HLA-DP':['DPA1','DPB1'],'HLA-DRB345':['DRB3','DRB4','DRB5']}
W=dict(exon=1,intron=0.1,intergenic=0.005)

def fasta_gz(path):
    out={};name=None
    with gzip.open(path,'rt') as f:
        for line in f:
            if line[0]=='>':name=line[1:].split()[0];out[name]=[]
            else:out[name].append(line.strip())
    return {k:''.join(v) for k,v in out.items()}

def load_ipd(d):
    idx={}
    for g in [g for gs in LOCUS_GENES.values() for g in gs]:
        rows=list(csv.DictReader(open(d/f'{g}.tsv'),delimiter='\t'))
        by_name={a:r for r in rows for a in r['alleles'].split(';')}
        idx[g]=(rows,by_name)
    return idx

def partial_near(cds,rows,k):
    """Partial IPD CDS (usually exon 2 only) -> full CDS with that span replaced, if it fits the
    backbone within k edits. Returns {new_cds: (representative allele, all alleles)}."""
    out={}
    for r in rows:
        res=edlib.align(r['sequence'],cds,mode='HW',task='locations',k=k)
        if res['editDistance']<0 or not res['locations']:continue
        i,j=res['locations'][0]
        new=cds[:i]+r['sequence']+cds[j+1:]
        if new!=cds:out[new]=(r['alleles'].split(';')[0],res['editDistance'],r['alleles'])
    return out

def t1k_alleles(path):
    """{gene: [allele names]} from a T1K genotype file (quality > 0, all listed alternatives)."""
    out=collections.defaultdict(list)
    if not path or not Path(path).exists():return out
    for line in open(path):
        v=line.rstrip('\n').split('\t');g=v[0].removeprefix('HLA-')
        for names,q in [(v[2],v[4]),(v[5],v[7])]:
            if names!='.' and float(q)>0:out[g]+=[n.removeprefix('HLA-') for n in names.split(',')]
    return out

def t1k_backbones(t1k_names,gene,labels,donor,hap,limit=3):
    """Leave-one-out panel haplotypes whose catalogue label matches a T1K two-field call."""
    want={norm2(n) for n in t1k_names if norm2(n)}
    out=[]
    for h,rows in labels.items():
        r=rows.get(gene)
        if not r or r['donor_id']==donor or h not in hap:continue
        names=r['exact_cds_alleles'] if r['cds_complete']=='1' and r['exact_cds_alleles'] else r['exact_protein_alleles']
        if names and {norm2(x) for x in names.split(';')}&want:out.append(h)
    # distinct sequences first, deterministic order
    seen=set();res=[]
    for h in sorted(out):
        if hap[h] in seen:continue
        seen.add(hap[h]);res.append(h)
    return res[:limit*max(1,len(want))]

def norm2(name):
    x=name.split('*')[-1].split(':')
    return ':'.join(x[:2]) if len(x)>=2 else None

def load_labels(path):
    out=collections.defaultdict(dict)
    with gzip.open(path,'rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):out[r['haplotype']][r['gene'].removeprefix('HLA-')]=r
    return out

def t1k_cds(names,by_name):
    """T1K may report 2- or 3-field names: take every complete IPD CDS whose allele name extends them."""
    seqs={}
    for n in names:
        for a,r in by_name.items():
            if a==n or a.startswith(n+':'):seqs[r['sequence']]=r['alleles'].split(';')[0]
    return seqs

def near(cds,rows,k):
    """{new_cds: (representative allele, edit distance)} for complete IPD CDS within k edits."""
    out={}
    for r in rows:
        if abs(len(r['sequence'])-len(cds))>k:continue
        d=edlib.align(r['sequence'],cds,mode='NW',task='distance',k=k)['editDistance']
        if 0<d<=k:out[r['sequence']]=(r['alleles'].split(';')[0],d)
    return out

def cap_per_type(cands,per_type):
    """cands: {seq: (allele, edits[, ...])}; keep the per_type nearest sequences of each two-field type."""
    by=collections.defaultdict(list)
    for seq,v in cands.items():by[norm2(v[0])].append((v[1],seq))
    keep={}
    for t,lst in by.items():
        for _,seq in sorted(lst)[:per_type]:keep[seq]=cands[seq]
    return keep

def main():
    p=argparse.ArgumentParser();p.add_argument('--donor',required=True);p.add_argument('--arm',default='full')
    p.add_argument('--top-k',type=int,default=5);p.add_argument('--min-prob',type=float,default=1e-3);p.add_argument('--max-edits',type=int,default=3)
    p.add_argument('--t1k',type=Path);p.add_argument('--no-t1k',action='store_true');p.add_argument('--partial',action='store_true',help='also graft partial IPD CDS (exon-level, e.g. DRB1)')
    p.add_argument('--per-type',type=int,default=2,help='IPD grafts kept per backbone and two-field type');p.add_argument('--max-candidates',type=int,default=500);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--loci',nargs='+',default=list(LOCUS_GENES));a=p.parse_args()
    ipd=load_ipd(H/'source/ipd_cds');partial_idx=load_ipd(H/'source/ipd_partial')
    cds_iv=collections.defaultdict(dict)
    for r in csv.DictReader(open(H/'weights/cds_intervals.tsv'),delimiter='\t'):
        cds_iv[r['haplotype']][r['gene']]=(r['strand'],[tuple(map(int,x.split('-'))) for x in r['cds'].split(',')])
    t1k={} if a.no_t1k else t1k_alleles(a.t1k)
    panel_labels=load_labels(H/'source/haplotype_labels.tsv.gz') if t1k else {}
    a.out.mkdir(parents=True,exist_ok=True);audit={}
    with open(a.out/'targets.bed','w') as bed,open(a.out/'labels.tsv','w') as lab,open(a.out/'reg_weights.tsv','w') as rw:
        lab.write('locus\tcandidate\tbackbone\tgene\tsource\tallele\tcds_sha256\talleles\n')
        for L in a.loci:
            res=json.load(gzip.open(H/f'calls/{a.arm}/{a.donor}/loci/{L}/res.json.gz','rt'))
            opts=sorted(res.get('options',[]),key=lambda o:-o['prob'])
            keep=[o for i,o in enumerate(opts) if i<a.top_k or o['prob']>=a.min_prob]
            backbones=sorted({h for o in keep for h in o['genotype'].split(',')})
            hap=fasta_gz(H/f'db/{a.arm}/loci/{L}/haplotypes.fa.gz')
            # locityper collapses identical haplotypes; --leave-out can name a collapsed alias
            for line in open(H/f'db/{a.arm}/loci/{L}/discarded_haplotypes.txt'):
                rep,others=line.split('=');rep=rep.strip()
                for o in others.split(','):
                    o=o.strip()
                    if o:hap[o]=hap[rep]
            if 'GRCh38.0' not in backbones:backbones.append('GRCh38.0')   # reference haplotype kept for locityper target
            extra_bb=[]
            for g in LOCUS_GENES[L]:
                for h in t1k_backbones(t1k.get(g,[]),g,panel_labels,a.donor,hap):
                    if h not in backbones and h not in extra_bb:extra_bb.append(h)
            backbones+=extra_bb
            seqs={};labels=[];ivs={}
            for b in backbones:
                ivs[b]={g:cds_iv[b][g][1] for g in LOCUS_GENES[L] if g in cds_iv[b]}
                seqs[b]=hap[b]
                for g in LOCUS_GENES[L]:
                    if g in cds_iv[b]:
                        st,iv=cds_iv[b][g];labels.append((b,b,g,'backbone','',hashlib.sha256(cds_of(hap[b],iv,st).encode()).hexdigest(),''))
                    else:labels.append((b,b,g,'backbone','ABSENT','',''))
            joint=collections.defaultdict(dict);pending=[]
            for b in backbones:
                for g in LOCUS_GENES[L]:
                    if g not in cds_iv[b]:continue
                    st,iv=cds_iv[b][g];old=cds_of(hap[b],iv,st)
                    cands=cap_per_type(near(old,ipd[g][0],a.max_edits),a.per_type)
                    tk=t1k_cds(t1k.get(g,[]),ipd[g][1]);joint[b][g]=tk
                    part=cap_per_type(partial_near(old,partial_idx[g][0],a.max_edits),a.per_type) if a.partial else {}
                    # priority: T1K grafts (0), complete IPD grafts by edits (1+d), partial grafts (10+d)
                    todo=[(0,0,s,al,'t1k','') for s,al in tk.items()]
                    todo+=[(1,d,s,al,'ipd_near','') for s,(al,d) in cands.items()]
                    todo+=[(2,d,s,al,'ipd_partial',alist) for s,(al,d,alist) in part.items()]
                    for pri,d,s,allele,src,alist in todo:
                        if s==old:continue
                        name=f"{b}~{'p' if src=='ipd_partial' else ''}{safe(allele)}"
                        if name in seqs:continue
                        pending.append((pri,d,name,b,g,s,src,allele,alist,iv,st))
                # joint grafts of T1K alleles for two-gene loci
                if len(LOCUS_GENES[L])==2 and all(joint[b].get(g) for g in LOCUS_GENES[L]):
                    g1,g2=LOCUS_GENES[L]
                    for s1,a1 in joint[b][g1].items():
                        for s2,a2 in joint[b][g2].items():
                            st,iv=cds_iv[b][g1];x,niv1=graft(hap[b],iv,st,s1)
                            # g2 intervals shift only if g1 lies before g2 and changed length
                            shift=len(x)-len(hap[b]);st2,iv2=cds_iv[b][g2]
                            if iv2[0][0]>=iv[-1][1]:iv2=[(u+shift,v+shift) for u,v in iv2]
                            name=f'{b}~{safe(a1)}~{safe(a2)}'
                            if name in seqs:continue
                            y,niv2=graft(x,iv2,st2,s2)
                            if niv1[0][0]>=iv2[-1][1]:niv1=[(u+len(y)-len(x),v+len(y)-len(x)) for u,v in niv1]
                            seqs[name]=y;ivs[name]={g1:niv1,g2:niv2}
                            labels.append((name,b,g1,'t1k_joint',a1,hashlib.sha256(s1.encode()).hexdigest(),''))
                            labels.append((name,b,g2,'t1k_joint',a2,hashlib.sha256(s2.encode()).hexdigest(),''))
            budget=max(0,a.max_candidates-len(seqs))
            for pri,d,name,b,g,s,src,allele,alist,iv,st in sorted(pending,key=lambda x:(x[0],x[1],x[2]))[:budget]:
                seqs[name],niv=graft(hap[b],iv,st,s);shift=len(seqs[name])-len(hap[b])
                ivs[name]={h:(niv if h==g else [(u+shift,v+shift) if u>=iv[-1][1] else (u,v) for u,v in x]) for h,x in ivs[b].items()}
                labels.append((name,b,g,src,allele,hashlib.sha256(s.encode()).hexdigest(),alist))
            strand={g:cds_iv[b][g][0] for b in backbones for g in cds_iv[b]}
            for name,_,g,src,allele,sha,_al in labels:   # every labelled CDS must be recoverable from its candidate
                if sha:assert hashlib.sha256(cds_of(seqs[name],ivs[name][g],strand[g]).encode()).hexdigest()==sha,(name,g)
            fa=a.out/f'{L}.fa'
            with open(fa,'w') as f:
                for n,s in seqs.items():f.write(f'>{n}\n{s}\n')
            for row in labels:lab.write(L+'\t'+'\t'.join(row)+'\n')
            wbed=a.out/f'{L}.weights.bed'
            with open(wbed,'w') as f:
                for n,sq in seqs.items():
                    marks=[(0,len(sq),W['intergenic'])]
                    for x in ivs.get(n,{}).values():
                        marks.append((x[0][0],x[-1][1],W['intron']));marks+=[(u,v,W['exon']) for u,v in x]
                    for u,v,w in tile(len(sq),marks):f.write(f'{n}\t{u}\t{v}\t{w}\n')
            rw.write(f'{L}\t{wbed.resolve()}\n')
            region=[x for x in open(H/'source/loci_spec.tsv') if x.startswith(L+'\t')][0].split('\t')
            spec=dict(zip(open(H/'source/loci_spec.tsv').readline().rstrip('\n').split('\t'),region))
            bed.write(f"ctg1\t{spec['ref_start0']}\t{spec['ref_end0'].strip()}\t{L}\t{fa.resolve()}\n")
            audit[L]=dict(backbones=len(backbones),t1k_backbones=len(extra_bb),candidates=len(seqs),proposed=len(pending)+len(seqs),options_kept=len(keep))
    (a.out/'audit.json').write_text(json.dumps(audit,indent=1)+'\n');print(json.dumps(audit))
if __name__=='__main__':main()
