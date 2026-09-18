import csv,gzip,re,collections,edlib
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
donors={r['donor']:r for r in csv.DictReader(open('/home/leechuck/Public/software/pangenome/hla-asian50/run/source/donors.tsv'),delimiter='\t')}
panel=collections.defaultdict(list)
for f in range(5):
    for r in csv.DictReader(gzip.open(f'panels/fold{f}/full/labels.tsv.gz','rt'),delimiter='\t'):panel[f,r['gene']].append(r)
rows=[r for r in csv.DictReader(open(R+'named_HLA_scores.tsv'),delimiter='\t') if r['method']=='Locityper_full' and r['endpoint'].startswith('assembly')]
tot=collections.Counter();out=[]
for r in rows:
    s,g=r['donor'],r['gene'];f=int(donors[s]['fold']);truth=r['truth'].split(';');pred=[]
    dists=[]
    for hap in ['1','2']:
        t=seqs[g].get(f'{s}#{hap}#HLA-{g}')
        if not t:pred.append(None);continue
        best=(10**9,None)
        for x in panel[f,g]:
            q=seqs[g][x['source_haplotypes'].split(';')[0]]
            d=edlib.align(t,q,task='distance',k=best[0])['editDistance']
            if d>=0 and d<best[0]:best=(d,x)
        ls={re.sub(r'^.*\*','',v) for v in best[1]['labels'].split(';') if v};pred.append(next(iter(ls)) if len(ls)==1 else None);dists.append(best[0])
    ok=sorted(map(str,pred))==sorted(truth)
    tot[g,'n']+=1;tot[g,'oracle']+=ok;tot[g,'loc']+=int(r['correct']);tot[g,'exact']+=all(d==0 for d in dists) and len(dists)==2;tot[g,'near']+=all(d<=5 for d in dists) and len(dists)==2
    out.append((s,g,r['truth'],pred,dists,ok,r['correct']))
for g in GENES:print(g,{k:tot[g,k] for k in ['n','loc','oracle','exact','near']})
print('ALL',{k:sum(tot[g,k] for g in GENES) for k in ['n','loc','oracle','exact','near']})
for o in out:
    if not o[5]:print('oracle-miss',*o)
