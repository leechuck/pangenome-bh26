#!/usr/bin/env python3
"""Named HLA comparison: graph-derived whole-locus PanGenie vs T1K/SpecHLA.
Experimental labels and assembly-sequence labels are distinct endpoints.
All eligible donor/loci remain in denominators, including prediction no-calls.
"""
import csv,gzip,tarfile,json,hashlib,re,collections,argparse,itertools
from pathlib import Path
from prepare_hla import GENES,labels

def table(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def norm(x):
    x=x.split('*')[-1];m=re.match(r'^(\d+):(\d+)',x)
    return ':'.join(m.groups()) if m and int(m[2]) else None

def unique(x):
    vals={norm(a) for a in re.split('[,;/]',x)}
    return next(iter(vals)) if len(vals)==1 and None not in vals else None

def compare(pair,truth):
    called=all(pair)
    correct=bool(called and any(all(p in t for p,t in zip(ps,truth)) for ps in [pair,pair[::-1]]))
    matches=max(sum(p is not None and p in t for p,t in zip(ps,truth)) for ps in [pair,pair[::-1]])
    return int(called),int(correct),matches

def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();r=a.run;a.out.mkdir(parents=True,exist_ok=True)
    donors=table(r/'source/donors.tsv');catalog=collections.defaultdict(list)
    for row in table(r/'hla_source/sequence_catalogue.tsv'):catalog[row['name']].append(row)
    preds={};rows=[];mapping={};panel_labels={}
    for fold in range(5):
        for arm in ['full','hprc']:
            with gzip.open(r/f'hla_panels/fold{fold}/{arm}/allele_labels.tsv.gz','rt') as f:
                rr=list(csv.DictReader(f,delimiter='\t'))
            for row in rr:
                vals=row['twofield_labels'].split(';') if row['twofield_labels'] else []
                mapping[fold,arm,row['gene'],row['sequence_sha256']]=[norm(x) for x in vals]
            for g in GENES:panel_labels[fold,arm,g]={norm(x) for row in rr if row['gene']==g for x in row['twofield_labels'].split(';') if x}
    with tarfile.open(r/'native_HLA_calls.tar.gz') as tar:
        for d in donors:
            s=d['donor'];fold=int(d['fold'])
            for line in tar.extractfile(f't1k/{s}/{s}_genotype.tsv').read().decode().splitlines():
                v=line.split('\t');g=v[0].removeprefix('HLA-')
                if g not in GENES:continue
                n=int(v[1]);p1=unique(v[2]) if float(v[4])>0 else None
                p2=p1 if n==1 else unique(v[5]) if n==2 and float(v[7])>0 else None
                preds[s,g,'T1K']=(p1,p2) if n else (None,None)
            lines=[x for x in tar.extractfile(f'spechla/{s}/hla.result.txt').read().decode().splitlines() if not x.startswith('#')]
            sp=dict(zip(lines[0].split('\t'),lines[1].split('\t')))
            for g in GENES:preds[s,g,'SpecHLA']=tuple(unique(sp[f'HLA_{g}_{h}']) for h in [1,2])
            for arm in ['full','hprc']:
                method='graph_HLA_'+arm;seen=set()
                with open(r/f'hla_calls/{s}/{arm}/calls_genotyping.vcf') as f:
                    for line in f:
                        if line.startswith('#'):continue
                        v=line.rstrip().split('\t');g=v[0].removeprefix('HLA-');assert g in GENES and g not in seen;seen.add(g)
                        aa=[v[3]]+v[4].split(',');gt=v[9].split(':')[v[8].split(':').index('GT')];parts=re.split(r'[/|]',gt);pair=[]
                        for x in parts if len(parts)==2 else ['.','.']:
                            if not x.isdigit():pair.append(None);continue
                            seq=aa[int(x)];ls=mapping[fold,arm,g,hashlib.sha256(seq.encode()).hexdigest()]
                            pair.append(ls[0] if len(ls)==1 else None)
                        preds[s,g,method]=tuple(pair)
                assert seen==set(GENES),(s,arm,seen)
            for g in GENES:
                for method in ['graph_HLA_full','graph_HLA_hprc','T1K','SpecHLA']:
                    pair=preds.get((s,g,method),(None,None))
                    rows.append(dict(donor=s,stratum=d['stratum'],population=d['population'],fold=fold,gene=g,method=method,allele1=pair[0] or '',allele2=pair[1] or '',called=int(all(pair))))
    write(a.out/'hla_predictions.tsv',rows)
    with open(r/'experimental_HLA.txt') as f:exper={x['id']:x for x in csv.DictReader(f,delimiter=' ')}
    scores=[];elig=[]
    for d in donors:
        s=d['donor'];fold=int(d['fold'])
        for g in GENES:
            truthsets={}
            if s in exper and g in exper[s]:
                ts=[{norm(x) for x in exper[s][k].split('/')} for k in [g,g+'.1']]
                if all(t and None not in t for t in ts):truthsets['experimental_twofield']=ts
            ts=[]
            for hap in ['1','2']:
                rr=catalog.get(f'{s}#{hap}#HLA-{g}',[]);ls=labels(rr[0]) if len(rr)==1 else set();ts.append({norm(x) for x in ls})
            if all(len(t)==1 and None not in t for t in ts):truthsets['assembly_exact_CDS_twofield']=ts
            for endpoint in ['experimental_twofield','assembly_exact_CDS_twofield']:
                elig.append(dict(donor=s,stratum=d['stratum'],gene=g,endpoint=endpoint,eligible=int(endpoint in truthsets)))
            for endpoint,truth in truthsets.items():
                for method in ['graph_HLA_full','graph_HLA_hprc','T1K','SpecHLA']:
                    pair=preds.get((s,g,method),(None,None));called,correct,matches=compare(pair,truth)
                    covered=''
                    if method.startswith('graph_HLA_'):
                        arm=method.removeprefix('graph_HLA_');covered=int(all(t & panel_labels[fold,arm,g] for t in truth))
                    scores.append(dict(donor=s,stratum=d['stratum'],population=d['population'],gene=g,method=method,endpoint=endpoint,truth=';'.join('/'.join(sorted(t)) for t in truth),prediction=';'.join(x or 'NO_CALL' for x in pair),called=called,correct=correct,allele_matches=matches,truth_labels_in_panel=covered))
    write(a.out/'hla_scores.tsv',scores);write(a.out/'hla_truth_eligibility.tsv',elig)
    groups=collections.defaultdict(list)
    for x in scores:
        for g in [x['gene'],'ALL']:groups[x['stratum'],x['endpoint'],x['method'],g].append(x)
    summary=[]
    for (s,e,m,g),rr in sorted(groups.items()):
        summary.append(dict(stratum=s,endpoint=e,method=m,gene=g,donors=len({x['donor'] for x in rr}),genotypes=len(rr),called=sum(x['called'] for x in rr),correct=sum(x['correct'] for x in rr),genotype_accuracy_pct=100*sum(x['correct'] for x in rr)/len(rr),alleles=2*len(rr),allele_matches=sum(x['allele_matches'] for x in rr),allele_accuracy_pct=50*sum(x['allele_matches'] for x in rr)/len(rr)))
    write(a.out/'hla_summary.tsv',summary)
    print(json.dumps([x for x in summary if x['gene']=='ALL'],indent=2))
if __name__=='__main__':main()
