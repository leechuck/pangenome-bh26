#!/usr/bin/env python3
"""Whole-locus PanGenie representation for named HLA typing from graph-input paths.
Eight primary loci, real gene +/-2kb regions. First/last64 bp are reference context;
the inner region is one sequence-resolved allele. Unrepresented genotypes stay missing.
Labels come only from each training panel plus GRCh38, never held-out labels.
"""
import csv,json,gzip,re,hashlib,collections,argparse
from pathlib import Path
GENES=['A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1']
def fasta(path):
    out={};name=None
    with open(path) as f:
        for l in f:
            if l.startswith('>'):name=l[1:].split()[0];assert name not in out;out[name]=[]
            else:out[name].append(l.strip().upper())
    return {k:''.join(v) for k,v in out.items()}
def labels(row):
    if not row or row['cds_complete']!='1':return set()
    result=set()
    for allele in row['exact_cds_alleles'].split(';'):
        m=re.fullmatch(r'(?:HLA-)?([A-Z0-9]+)\*(\d+):(\d+)(?::\d+)*[A-Z]*',allele)
        if m:result.add(m[1]+'*'+m[2]+':'+m[3])
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    catalog=collections.defaultdict(list)
    for r in csv.DictReader(open(a.source/'sequence_catalogue.tsv'),delimiter='\t'):catalog[r['name']].append(r)
    seqs={g:fasta(a.source/f'HLA-{g}.fa') for g in GENES}
    regions={g:{r['name']:r for r in csv.DictReader(open(a.source/f'HLA-{g}.regions.tsv'),delimiter='\t')} for g in GENES}
    reference=next(iter(fasta(a.run/'graph/reference.fa').values()));background=list(reference)
    refgenes={};checks=[]
    for g in GENES:
        name=f'GRCh38#0#HLA-{g}';rr=regions[g][name];lo,hi=map(int,rr['region'].rsplit(':',1)[1].split('-'));s=reference[lo-1:hi]
        if rr['strand']=='-':s=s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
        assert s==seqs[g][name],('reference extraction mismatch',g)
        assert hi-lo+1==len(s);background[lo-1:hi]='N'*(hi-lo+1);refgenes[g]=s
        checks.append(dict(gene=g,lo=lo,hi=hi,strand=rr['strand'],length=len(s),sha256=hashlib.sha256(s.encode()).hexdigest()))
    a.out.mkdir(parents=True,exist_ok=True)
    with open(a.out/'reference.fa','w') as f:
        f.write('>MHC_background\n'+''.join(background)+'\n')
        for g in GENES:f.write(f'>HLA-{g}\n{refgenes[g]}\n')
    (a.out/'reference_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    stats=[]
    for fold in range(5):
        for arm in ['full','hprc']:
            samples=(a.run/'panels'/f'fold{fold}'/arm/'samples.txt').read_text().splitlines()
            path=a.out/f'fold{fold}'/arm;path.mkdir(parents=True,exist_ok=True)
            (path/'samples.txt').write_text('\n'.join(samples)+'\n')
            with open(path/'bubbles.vcf','w') as out,gzip.open(path/'allele_labels.tsv.gz','wt') as lf:
                out.write('##fileformat=VCFv4.2\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Phased whole-locus haplotype">\n')
                out.write('\t'.join(['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']+samples)+'\n')
                writer=csv.DictWriter(lf,fieldnames=['gene','allele_index','sequence_sha256','twofield_labels','source_haplotypes','unlabelled_sources'],delimiter='\t',lineterminator='\n');writer.writeheader()
                for g in GENES:
                    ref=refgenes[g][64:-64];alleles=[ref];idx={ref:0};sources=collections.defaultdict(list);sources[0].append(f'GRCh38#0#HLA-{g}')
                    gt=[];missing=0
                    for sample in samples:
                        hapg=[]
                        for hap in ['1','2']:
                            name=f'{sample}#{hap}#HLA-{g}';seq=seqs[g].get(name,'');seq=seq[64:-64]
                            if not seq or not re.fullmatch('[ACGT]+',seq):hapg.append('.');missing+=1;continue
                            if seq not in idx:idx[seq]=len(alleles);alleles.append(seq)
                            i=idx[seq];sources[i].append(name);hapg.append(str(i))
                        gt.append('|'.join(hapg))
                    assert len(alleles)>1
                    out.write('\t'.join([f'HLA-{g}','65',f'HLA-{g}',ref,','.join(alleles[1:]),'.','PASS','.','GT']+gt)+'\n')
                    unknown=0
                    for i,seq in enumerate(alleles):
                        ls=set();bad=0
                        for name in sources[i]:
                            rr=catalog.get(name,[]);ll=labels(rr[0]) if len(rr)==1 else set()
                            ls.update(ll);bad+=not bool(ll)
                        # An unlabelled copy cannot invalidate an identical sequence with a known label.
                        unknown+=len(ls)!=1
                        writer.writerow(dict(gene=g,allele_index=i,sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),twofield_labels=';'.join(sorted(ls)),source_haplotypes=';'.join(sources[i]),unlabelled_sources=bad))
                    stats.append(dict(fold=fold,arm=arm,gene=g,training_donors=len(samples),alleles=len(alleles),missing_haplotypes=missing,ambiguous_or_unknown_label_alleles=unknown))
    (a.out/'audit.json').write_text(json.dumps(stats,indent=2)+'\n');print('Built',len(stats),'fold/arm/locus panels',flush=True)
if __name__=='__main__':main()
