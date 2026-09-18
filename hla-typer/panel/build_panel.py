#!/usr/bin/env python3
"""Build Locityper input panels (one per arm) from anchor-bounded locus sequences.

Arms (reference haplotypes GRCh38/CHM13 always included; donors sampled as whole donors):
  full            every panel donor
  hprc            HPRC r2 donors only
  asian_matched   all non-HPRC donors + random HPRC donors, same haplotype count as hprc
  random_matched  random donors, same haplotype count as hprc
Test donors are never removed here: leave-one-out happens at genotyping time (locityper --leave-out).
A haplotype enters a locus panel only if its locus cut is status ok and ACGT-only.
"""
import argparse,csv,gzip,json,random,hashlib,collections
from pathlib import Path
SEED='hla-typer-panels-20260917'

def fasta_gz(path):
    out={};name=None
    with gzip.open(path,'rt') as f:
        for line in f:
            if line[0]=='>':name=line[1:].split()[0];out[name]=[]
            else:out[name].append(line.strip())
    return {k:''.join(v) for k,v in out.items()}

def arms(donors):
    """donors: {unit: (cohort, [hap_id,...])} -> {arm: set(hap_id)}"""
    ref={h for d,(c,hs) in donors.items() if c=='REF' for h in hs}
    real={d:v for d,v in donors.items() if v[0]!='REF'}
    hprc=sorted(d for d,(c,_) in real.items() if c.startswith('HPRC'))
    other=sorted(d for d in real if d not in set(hprc))
    nh=sum(len(real[d][1]) for d in hprc)
    def fill(fixed,pool,tag):
        rng=random.Random(SEED+tag);pool=sorted(pool);rng.shuffle(pool);chosen=list(fixed)
        n=sum(len(real[d][1]) for d in chosen)
        for d in pool:
            if n>=nh:break
            chosen.append(d);n+=len(real[d][1])
        return chosen
    sel={'full':sorted(real),'hprc':hprc,'asian_matched':fill(other,hprc,'asian'),'random_matched':fill([],sorted(real),'random')}
    return {a:ref|{h for d in ds for h in real[d][1]} for a,ds in sel.items()}

def main():
    p=argparse.ArgumentParser();p.add_argument('--loci',type=Path,default='results/loci');p.add_argument('--spec',type=Path,default='source/loci_spec.tsv')
    p.add_argument('--donors',type=Path,default='source/haplotype_donors.tsv');p.add_argument('--reference',type=Path,required=True)
    p.add_argument('--out',type=Path,default='panels');a=p.parse_args()
    hd=list(csv.DictReader(open(a.donors),delimiter='\t'));name={r['hap_id']:r['haplotype'] for r in hd}
    # Sampling unit = donor x assembly source, so the JaSaPaGe re-assemblies of HPRC JPT donors never enter the
    # HPRC arm; leave-one-out at genotyping time removes every unit of a donor.
    donors={}
    for r in hd:
        unit=(r['donor_id'],'HPRC' if r['cohort'].startswith('HPRC') else r['cohort'])
        donors.setdefault(unit,(r['cohort'],[]))[1].append(r['hap_id'])
    A=arms(donors)
    status={(r['locus'],r['hap_id']):r for r in csv.DictReader(open(a.loci/'locus_haplotypes.tsv'),delimiter='\t')}
    ref=''.join(l.strip() for l in open(a.reference) if l[0]!='>')
    spec=list(csv.DictReader(open(a.spec),delimiter='\t'));audit=[]
    seqs={s['locus']:fasta_gz(a.loci/f"{s['locus']}.fa.gz") for s in spec}
    for s in spec:
        refseq=ref[int(s['ref_start0']):int(s['ref_end0'])]
        assert seqs[s['locus']]['GRCh38#0']==refseq,('reference cut differs from ctg1 interval',s['locus'])
    for arm,haps in A.items():
        d=a.out/arm;d.mkdir(parents=True,exist_ok=True)
        with open(d/'targets.bed','w') as bed,open(d/'members.tsv','w') as mem:
            mem.write('locus\thap_id\thaplotype\tincluded\treason\tsha256\n')
            for s in spec:
                L=s['locus'];fa=d/f'{L}.fa';n=0;reasons=collections.Counter()
                with open(fa,'w') as out:
                    for h in sorted(haps):
                        st=status.get((L,h));seq=seqs[L].get(h,'')
                        reason='' if st and st['status']=='ok' and seq and set(seq)<=set('ACGT') else (st['status'] if st and st['status']!='ok' else 'non_ACGT')
                        sha=hashlib.sha256(seq.encode()).hexdigest() if seq else ''
                        mem.write(f"{L}\t{h}\t{name[h]}\t{int(not reason)}\t{reason}\t{sha}\n")
                        if reason:reasons[reason]+=1;continue
                        out.write(f'>{name[h]}\n{seq}\n');n+=1
                bed.write(f"ctg1\t{s['ref_start0']}\t{s['ref_end0']}\t{L}\t{fa.resolve()}\n")
                audit.append(dict(arm=arm,locus=L,haplotypes_in_arm=len(haps),included=n,excluded=dict(reasons)))
    (a.out/'panel_audit.json').write_text(json.dumps(dict(seed=SEED,arm_sizes={k:len(v) for k,v in A.items()},loci=audit),indent=1)+'\n')
    for x in audit:print(x['arm'],x['locus'],x['included'],x['excluded'])
if __name__=='__main__':main()
