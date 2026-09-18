import csv,collections,sys
sys.path.insert(0,'../hla-bench');sys.path.insert(0,'stageB')
from names import Nomenclature,pred_slot,compare
import score_calls as sc
arm=sys.argv[1] if len(sys.argv)>1 else 'full-ens'
rows=list(csv.DictReader(open('results/c1_scores/scores.tsv'),delimiter='\t'))
ens=[r for r in rows if r['arm']==arm]
loci=collections.defaultdict(set)
for r in ens:loci[r['donor']].add(r['locus'])
done={d for d in loci if len(loci[d])==7}
C8={'A','B','C','DRB1','DQA1','DQB1','DPA1','DPB1'}
sb={(r['donor'],r['gene'],r['level']):int(r['correct']) for r in csv.DictReader(open('results/stageB/scores.tsv'),delimiter='\t') if r['arm']=='full'}
nom=Nomenclature();by_hap,hd=sc.load_labels('source/haplotype_labels.tsv.gz')
dh=collections.defaultdict(list)
for h,(d,c) in hd.items():dh[d].append((c,h))
def t1k_ok(d,g,level):
    if level not in ('two_field','g_group'):return None
    calls={}
    for line in open(f'results/comparators/t1k/{d}/{d}_genotype.tsv'):
        v=line.rstrip('\n').split('\t');calls[v[0].removeprefix('HLA-')]=v
    own=sorted(h for c,h in dh[d] if not c.startswith('JaSaPaGe'))[:2]
    truth=[sc.slot(by_hap[h].get(g),level,nom,g,True) for h in own]
    if len(truth)!=2 or None in truth:return None
    v=calls.get(g);n=int(v[1]) if v else 0
    a1=pred_slot(v[2],level,nom,g) if v and n>=1 and float(v[4])>0 else None
    a2=a1 if (v and n==1) else (pred_slot(v[5],level,nom,g) if v and n==2 and float(v[7])>0 else None)
    if n==0 and g in sc.PRESENCE_GENES:a1=a2=[sc.ABSENT]
    return compare([a1,a2],truth)[1]
agg=collections.defaultdict(lambda:[0,0,0,0]);src=collections.Counter();flips=[]
for r in ens:
    if r['donor'] not in done:continue
    k=(r['donor'],r['gene'],r['level']);scope='classical8' if r['gene'] in C8 else 'DRB345'
    t=t1k_ok(r['donor'],r['gene'],r['level'])
    a=agg[scope,r['level']];a[0]+=int(r['correct']);a[1]+=sb.get(k,0);a[2]+=(t or 0);a[3]+=1
    if r['level']=='two_field':
        src[r['sources']]+=1
        if int(r['correct'])!=sb.get(k,0):flips.append((r['donor'],r['gene'],'ens_right' if r['correct']=='1' else 'stageB_right',r['picks'],r['sources'],'T1K_ok' if t else 'T1K_wrong'))
print(arm,'donors complete:',len(done))
for (scope,level),(e,b,t,n) in sorted(agg.items()):print(f'{scope:10s} {level:10s} n={n:3d} {arm} {e} ({100*e/n:.1f}%)  stageB {b} ({100*b/n:.1f}%)  T1K {t if level in ("two_field","g_group") else "-"}')
print('pick sources (two-field rows):',dict(src))
print('flips vs stage B:');[print(' ',*f) for f in flips]
