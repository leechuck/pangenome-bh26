import csv,gzip,json,re,collections,sys,edlib
csv.field_size_limit(10**9)
R='/home/leechuck/Public/software/pangenome/hla-asian50/refined/results/'
def fasta(p):
    o={};n=None
    for l in open(p):
        if l[0]=='>':n=l[1:].split()[0];o[n]=[]
        else:o[n].append(l.strip().upper())
    return {k:''.join(v) for k,v in o.items()}
GENES=['A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1']
seqs={g:fasta(f'hla_source/HLA-{g}.fa') for g in GENES}
cat={}
for r in csv.DictReader(open('hla_source/sequence_catalogue.tsv'),delimiter='\t'):cat.setdefault(r['name'],r)
def lab(name):
    r=cat.get(name)
    if not r:return 'ABSENT'
    s=set()
    for a in r['exact_cds_alleles'].split(';'):
        m=re.fullmatch(r'(?:HLA-)?([A-Z0-9]+)\*(\d+):(\d+)(?::\d+)*[A-Z]*',a)
        if m:s.add(m[2]+':'+m[3])
    return '/'.join(sorted(s)) or f"NOLABEL(cds_complete={r['cds_complete']},old={r['old_consensus']},{r['assessment']})"
donors={r['donor']:r for r in csv.DictReader(open('/home/leechuck/Public/software/pangenome/hla-asian50/run/source/donors.tsv'),delimiter='\t')}
panel={}
for f in range(5):
    for r in csv.DictReader(gzip.open(f'panels/fold{f}/full/labels.tsv.gz','rt'),delimiter='\t'):panel[f,r['gene'],r['haplotype']]=r
def d(a,b):return edlib.align(a,b,task='distance')['editDistance']
rows=[r for r in csv.DictReader(open(R+'named_HLA_scores.tsv'),delimiter='\t') if r['method']=='Locityper_full' and r['endpoint'].startswith('assembly') and r['correct']=='0']
for r in rows:
    s,g=r['donor'],r['gene'];f=int(donors[s]['fold'])
    res=json.load(gzip.open(f'loc/locityper/calls/{s}/full/loci/HLA-{g}/res.json.gz','rt'))
    picks=res['genotype'].split(',')
    print(f"\n== {s} {g} truth={r['truth']} pred={r['prediction']} q={res['quality']:.0f} reads={res['total_reads']}")
    P={h:seqs[g][x['source_haplotypes'].split(';')[0]] for (ff,gg,h),x in panel.items() if ff==f and gg==g}
    for hap in ['1','2']:
        n=f'{s}#{hap}#HLA-{g}';t=seqs[g].get(n)
        if not t:print('  truth hap',hap,'missing');continue
        ds=sorted((d(t,q),h) for h,q in P.items())
        best=ds[0];tl=lab(n)
        print(f"  truth h{hap} {tl} len={len(t)}; nearest panel: "+', '.join(f"{h}[{panel[f,g,h]['labels'] or 'NOLABEL'}]d={x}" for x,h in ds[:3])+" | picked: "+', '.join(f"{h}[{panel[f,g,h]['labels'] or 'NOLABEL'}]d={d(t,P[h])}" for h in picks))
    for h in picks:
        x=panel[f,g,h]
        if not x['labels'] or ';' in x['labels']:
            print('  unlabelled pick',h,'sources:',[ (n,lab(n)) for n in x['source_haplotypes'].split(';')][:4])
