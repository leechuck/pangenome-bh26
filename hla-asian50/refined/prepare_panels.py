#!/usr/bin/env python3
"""Retain callable sites with PanGenie's native phased-missing GT support.
No imputation; known haplotype alleles survive another haplotype's missingness.
All training/test membership and coordinates remain frozen from the first run.
"""
import csv,gzip,json,collections,re,argparse
from pathlib import Path

def clean_gt(gt,alleles,conflict=False):
    if conflict or not re.fullmatch(r'(?:\d+|\.)\|(?:\d+|\.)',gt):return (None,None)
    out=[]
    for x in gt.split('|'):
        if x=='.':out.append(None);continue
        i=int(x);assert i<len(alleles)
        out.append(i if re.fullmatch('[ACGT]+',alleles[i]) else None)
    return tuple(out)

def remap_site(alleles,gts):
    used={0}|{a for g in gts for a in g if a is not None}
    if len(used)==1:return None
    order=sorted(used);remap={a:str(i) for i,a in enumerate(order)}
    return order,['|'.join('.' if a is None else remap[a] for a in g) for g in gts]

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    rows=list(csv.DictReader(open(a.run/'source/panel_manifest.tsv'),delimiter='\t'));exclude=list(csv.DictReader(open(a.run/'source/excluded_paths.tsv'),delimiter='\t'))
    by=collections.defaultdict(list)
    for r in rows:by[r['sample']].append(r)
    valid={s for s,rs in by.items() if len(rs)==2 and s==rs[0]['donor_id'] and {r['haplotype'] for r in rs}=={'1','2'}}
    contexts=[]
    for k in range(5):
        ex={r['sample'] for r in exclude if int(r['fold'])==k}
        for arm in ['full','hprc']:
            ss={s for s in valid-ex if arm=='full' or by[s][0]['cohort'].startswith('HPRC')}
            assert not ss&ex
            path=a.out/'panels'/f'fold{k}'/arm;path.mkdir(parents=True,exist_ok=True)
            contexts.append(dict(fold=k,arm=arm,samples=ss,path=path,stats=collections.Counter(),last_end=-1,keys=set()))
    with gzip.open(a.run/'graph/full.top.vcf.gz','rt') as f:
        for line in f:
            if line.startswith('##'):continue
            v=line.rstrip().split('\t')
            if line.startswith('#CHROM'):
                for c in contexts:
                    c['idx']=[(i+9,s) for i,s in enumerate(v[9:]) if s in c['samples']]
                    assert len(c['idx'])==len(c['samples'])
                    (c['path']/'samples.txt').write_text('\n'.join(s for _,s in c['idx'])+'\n')
                    c['file']=open(c['path']/'bubbles.vcf','w');c['map']=gzip.open(c['path']/'allele_map.tsv.gz','wt')
                    c['map'].write('ID\toriginal_allele_indices\n')
                    c['file'].write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased genotype; missing haplotype alleles retained">\n'+'\t'.join(v[:9]+[s for _,s in c['idx']])+'\n')
                continue
            aa=[v[3]]+v[4].split(',');pos=int(v[1]);gi=v[8].split(':').index('GT');conflicts=set()
            for tag in v[7].split(';'):
                if tag.startswith('CONFLICT='):conflicts.update(tag.split('=',1)[1].split(','))
            for c in contexts:
                st=c['stats'];st['input_sites']+=1
                if not re.fullmatch('[ACGT]+',aa[0]):st['invalid_reference']+=1;continue
                gts=[clean_gt(v[i].split(':')[gi],aa,s in conflicts) for i,s in c['idx']]
                result=remap_site(aa,gts)
                if result is None:st['no_known_training_alt']+=1;continue
                # Full top-level VCF should not overlap; fail instead of silently losing baseline sites.
                assert pos>c['last_end'],('overlapping top-level sites',c['fold'],c['arm'],pos)
                order,genos=result
                c['file'].write('\t'.join(v[:4]+[','.join(aa[x] for x in order[1:]),'.','PASS','.','GT']+genos)+'\n')
                c['map'].write(v[2]+'\t'+','.join(map(str,order))+'\n')
                c['last_end']=pos+len(aa[0])-1;c['keys'].add((pos,aa[0]));st['retained_sites']+=1
                missing=sum(x is None for g in gts for x in g);st['missing_haplotype_alleles_retained']+=missing;st['sites_with_missing_retained']+=bool(missing)
    for c in contexts:
        c['file'].close();c['map'].close()
        old=set()
        with open(a.run/'panels'/f"fold{c['fold']}"/c['arm']/'bubbles.vcf') as f:
            for line in f:
                if not line.startswith('#'):
                    v=line.split('\t');old.add((int(v[1]),v[3]))
        assert old<=c['keys'],('lost original callable sites',c['fold'],c['arm'],len(old-c['keys']))
        full=next(x for x in contexts if x['fold']==c['fold'] and x['arm']=='full')
        assert c['keys']<=full['keys']
        audit=dict(fold=c['fold'],arm=c['arm'],training_donors=len(c['idx']),counts=dict(c['stats']),old_sites_preserved=len(old),new_sites=len(c['keys']-old),old_sites_lost=0,full_contains_HPRC_sites=True,missing_policy='native phased missing alleles; no reference filling or imputation')
        (c['path']/'audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(audit,flush=True)
if __name__=='__main__':main()
