#!/usr/bin/env python3
"""Predict without reading held-out assembly labels; dosage and imputation kept separate."""
import argparse,csv,json,os,re
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def pair_scores(X,y):
    X=X.astype(float);gram=X@X.T;z=np.diag(gram)-2*(X@y)
    s=z[:,None]+z[None,:]+2*gram+float(y@y)
    s[np.tril_indices(len(X),-1)]=np.inf
    return s

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--cohort',choices=['development','validation'],required=True);args=parser.parse_args()
    model=json.loads((ROOT/'source/dosage_model.json').read_text())
    rows=json.loads((ROOT/'source/training_rows.json').read_text());X=np.load(ROOT/'source/training_profiles.npy')
    B=np.stack([np.fromfile(ROOT/'source/background'/f'{r["hap_id"].replace("#","_")}.bin',dtype='<u4') for r in rows])
    keys=np.array([int(x) for x in (ROOT/'source/markers.txt').read_text().splitlines()],dtype=np.uint64)
    original=set(map(int,(ROOT.parent/'hla-structural/source/markers.txt').read_text().splitlines()))
    old=np.array([int(k) in original for k in keys])
    masks=np.array([int(s.split()[1]) for s in (ROOT/'source/marker_groups.tsv').read_text().splitlines()])
    # Sequence-derived training categories include the noncanonical diagnostic motif.
    genes={r['gene_id']:r for r in read(ROOT/'source/training_C4_genes.tsv')}
    dosage={r['hap_id']:np.zeros(6,dtype=int) for r in rows}
    for g in read(ROOT/'results/training_C4_diagnostic_audit.tsv'):
        f=g['inferred_form'];d=dosage[genes[g['gene_id']]['hap_id']];d[0]+=1;d[1 if f[2]=='A' else 2 if f[2]=='B' else 3]+=1;d[4 if f[-1]=='L' else 5]+=1
    D=np.stack([dosage[r['hap_id']] for r in rows])
    if args.cohort=='validation': queries=read(ROOT/'source/validation_donors.tsv')
    else:
        families={r['donor']:r['family'] for r in rows}
        queries=[dict(donor=d,family=families[d]) for d in (ROOT/'source/all_development_donors.txt').read_text().splitlines()]
    results=[];pairresults=[]
    for query in queries:
        donor=query['donor'];family=query['family'];ids=[i for i,r in enumerate(rows) if r['family']!=family]
        T=X[ids];ds=D[ids];types=np.array([rows[i]['structure_id'] for i in ids]);cn=ds[:,0]
        support=np.zeros(X.shape[1],int)
        for fam in {rows[i]['family'] for i in ids}:support+=(T[[j for j,i in enumerate(ids) if rows[i]['family']==fam]]>0).any(axis=0)
        keep=(support>=3)&(T==B[ids]).all(axis=0)&(T.max(axis=0)<=8)
        anchors=keep&old&((T==1).mean(axis=0)>=.98)
        conserved=keep&old&((T==cn[:,None]).mean(axis=0)>=.98)
        raw=np.fromfile(ROOT/f'source/read_counts/{args.cohort}/{donor}.bin',dtype='<u4').astype(float)
        assert len(raw)==len(keys)
        scale=float(np.median(raw[anchors])/2)
        estimate=float(np.median(raw[conserved])/scale/model['total_C4_factor']) if scale else float('nan')
        n=int(np.floor(estimate+.5)) if np.isfinite(estimate) else -1
        g={r['metric']:int(r['fragments']) for r in read(ROOT/f'source/read_counts/{args.cohort}/{donor}.groups.tsv')}
        ab=np.array([g['A_only']/model['A_relative_B']['factor'],g['B_only'],g['OTHER_only']/model['OTHER_relative_B']])
        ls=np.array([(g['LONG_LEFT']+g['LONG_RIGHT'])/2/model['L_relative_S']['factor'],g['SHORT']])
        imbalance=abs(g['LONG_LEFT']-g['LONG_RIGHT'])/max(g['LONG_LEFT']+g['LONG_RIGHT'],1)
        abest=n*ab/ab.sum() if ab.sum() else np.full(3,np.nan);lsest=n*ls/ls.sum() if ls.sum() else np.full(2,np.nan)
        # Sum-preserving nearest integer allocation, not independent rounding.
        def allocation(v,total):
            if not np.isfinite(v).all() or total<0:return np.full(len(v),-1)
            a=np.floor(v).astype(int)
            for j in np.argsort(-(v-a))[:total-int(a.sum())]:a[j]+=1
            return a
        abc=allocation(abest,n);lsc=allocation(lsest,n)
        cnok=scale>=3 and anchors.sum()>=30 and conserved.sum()>=10 and 1<=n<=8 and abs(estimate-n)<=model['maximum_rounding_distance']
        abok=cnok and (g['A_only']+g['B_only']+g['OTHER_only'])>=model['minimum_diagnostic_fragments'] and np.max(abs(abest-abc))<=model['maximum_rounding_distance'] and g['AB_ambiguous']==0
        lsok=cnok and (g['LONG_LEFT']+g['LONG_RIGHT']+g['SHORT'])>=model['minimum_diagnostic_fragments'] and np.max(abs(lsest-lsc))<=model['maximum_rounding_distance'] and imbalance<=model['maximum_long_junction_imbalance'] and g['LS_ambiguous']==0
        rec=dict(donor=donor,cohort=args.cohort,haploid_marker_depth=scale,C4_estimate=estimate,C4_call=n if cnok else '',AB_status='call' if abok else 'ambiguous',LS_status='call' if lsok else 'ambiguous',long_junction_imbalance=imbalance,BOUNDARY_fragments=g['BOUNDARY'],NONCANONICAL_fragments=g['NONCANONICAL'])
        for j,label in enumerate(['A','B','OTHER']):rec[label+'_estimate']=abest[j];rec[label+'_call']=int(abc[j]) if abok else ''
        for j,label in enumerate(['L','S']):rec[label+'_estimate']=lsest[j];rec[label+'_call']=int(lsc[j]) if lsok else ''
        results.append(rec)
        # The full compatible set expresses dosage ambiguity. A ranked path is imputation.
        allowed=np.triu(np.ones((len(ids),len(ids)),bool))
        if cnok:allowed &= cn[:,None]+cn[None,:]==n
        else:allowed[:]=False
        cn_allowed=allowed.copy()
        if abok:
            for k in range(3):allowed &= ds[:,k+1,None]+ds[None,:,k+1]==abc[k]
        if lsok:allowed &= ds[:,4,None]+ds[None,:,4]==lsc[0]
        def pairs(a):return sorted({tuple(sorted((types[i],types[j]))) for i,j in np.argwhere(a)})
        compatible=pairs(allowed)
        for name,features,candidates in [('old_sketch_CN',keep&old,cn_allowed),('targeted_dosage_paths',keep,allowed),('targeted_without_module_context',keep&((masks&32)==0),allowed)]:
            scores=pair_scores(T[:,features],raw[features]/max(scale,1e-9));scores[~candidates]=np.inf
            pair='';margin=''
            if np.isfinite(scores).any():
                ij=np.unravel_index(np.argmin(scores),scores.shape);best=scores[ij];bestpairs=pairs(np.isclose(scores,best,rtol=1e-9,atol=1e-7))
                if len(bestpairs)==1:
                    pair=';'.join(bestpairs[0]);a,b=bestpairs[0];same=((types[:,None]==a)&(types[None,:]==b))|((types[:,None]==b)&(types[None,:]==a));margin=float((np.min(np.where(same,np.inf,scores))-best)/max(best,1))
            method_compatible=pairs(candidates)
            pairresults.append(dict(donor=donor,method=name,ranked_pair=pair,status='imputed_not_physically_phased' if pair else 'no_unique_ranked_pair',relative_margin=margin,dosage_compatible_pair_count=len(method_compatible),dosage_compatible_pairs=json.dumps(method_compatible,separators=(',',':')),fully_called_marginals=int(cnok and abok and lsok)))
        print(donor,rec['C4_call'],rec['AB_status'],rec['LS_status'],len(compatible),flush=True)
    write(ROOT/f'results/{args.cohort}_dosage_predictions.tsv',results);write(ROOT/f'results/{args.cohort}_path_predictions.tsv',pairresults)
if __name__=='__main__':main()
