#!/usr/bin/env python3
"""Compare frozen predictions with assembly annotations; never called by inference."""
import argparse,csv,gzip,json,re
from collections import defaultdict
from pathlib import Path
import numpy as np
from design_dosage_model import labels
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows):
    with open(p,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cohort',choices=['development','validation'],required=True);args=ap.parse_args();cohort=args.cohort
    rows=[r for r in read(ROOT.parent/'hla-structural/source/catalogue.tsv') if r['locus']=='RCCX'];bydonor=defaultdict(list)
    for r in rows:bydonor[r['donor']].append(r)
    truth={}
    for donor,rs in bydonor.items():
        fs=[g for r in rs for g in re.findall(r'C4([AB][LS])[+-]',r['structural_signature'])]
        truth[donor]={'total':len(fs),'A':sum(f[0]=='A' for f in fs),'B':sum(f[0]=='B' for f in fs),'OTHER':0,'L':sum(f[1]=='L' for f in fs),'S':sum(f[1]=='S' for f in fs)}
    # Development diagnostic truth distinguishes the three noncanonical genes.
    truth.update(labels())
    pred=read(ROOT/f'results/{cohort}_dosage_predictions.tsv');metrics=[];errors=[]
    def metric(method,feature,items):
        metrics.append(dict(cohort=cohort,method=method,feature=feature,n=len(items),called=sum(p is not None for _,p,t in items),correct=sum(p is not None and p==t for _,p,t in items)))
        for donor,p,t in items:
            if p is None or p!=t:errors.append(dict(donor=donor,method=method,feature=feature,prediction=p,truth=t))
    def integer(s):return int(float(s)) if s not in ('',None) else None
    for feature in ['total','A','B','OTHER','L','S']:
        col=('C4' if feature=='total' else feature)+'_call'
        metric('targeted_markers',feature,[(r['donor'],integer(r[col]),truth[r['donor']][feature]) for r in pred])
    # Dedicated comparator: count all intended samples, including failed/missing outputs.
    base={}
    for r in pred:
        donor=r['donor'];p=ROOT/f'source/baseline/{cohort}/{donor}/results/C4Investigator_c4_summary.csv'
        base[donor]=list(csv.DictReader(open(p)))[0] if p.exists() else {}
    baseline_donors=set((ROOT/'source/development_donors.txt').read_text().splitlines()) if cohort=='development' else {r['donor'] for r in pred}
    for feature in ['total','A','B','L','S']:
        col='C4'+('' if feature=='total' else feature)+'_copy'
        metric('C4Investigator',feature,[(r['donor'],integer(base[r['donor']].get(col)),truth[r['donor']][feature]) for r in pred if r['donor'] in baseline_donors])
    for method in {r['method'] for r in read(ROOT/f'results/{cohort}_path_predictions.tsv')}:
        ps=[r for r in read(ROOT/f'results/{cohort}_path_predictions.tsv') if r['method']==method]
        metric(method,'full_signature_pair',[(r['donor'],r['ranked_pair'] or None,';'.join(sorted(x['structure_id'] for x in bydonor[r['donor']]))) for r in ps])
        metric(method,'truth_in_dosage_compatible_set',[(r['donor'],True,sorted(x['structure_id'] for x in bydonor[r['donor']]) in json.loads(r['dosage_compatible_pairs'])) for r in ps])
    # Fair reference-depth baseline, using factors already fitted in earlier pilot.
    factors={r['control']:float(r['calibration_factor']) for r in read(ROOT.parent/'hla-structural/results/calibrated_depth_baseline.tsv')}
    positions=defaultdict(set)
    for line in open(ROOT.parent/'hla-analysis/source/rccx_depth_targets.bed'):
        _,a,b,name=line.rstrip().split('\t');positions[name].update(range(int(a)+1,int(b)+1))
    depth=[]
    for r in pred:
        donor=r['donor'];p=ROOT/f'source/read_counts/{cohort}/{donor}.depth.tsv.gz'
        with gzip.open(p,'rt') as f:d={int(x[1]):int(x[2]) for x in (line.split('\t') for line in f)}
        assert all(s<=d.keys() for s in positions.values()),donor
        means={k:float(np.mean([d[i] for i in s])) for k,s in positions.items()}
        for control,factor in factors.items():
            value=2*(means['C4A']+means['C4B'])/means[control]/factor if means[control] else None
            depth.append(dict(donor=donor,control=control,estimate=value,call=int(np.floor(value+.5)) if value is not None else None,truth=truth[donor]['total']))
    for control in factors:metric('reference_depth_'+control,'total',[(r['donor'],r['call'],r['truth']) for r in depth if r['control']==control])
    write(ROOT/f'results/{cohort}_metrics.tsv',metrics)
    if errors:write(ROOT/f'results/{cohort}_errors.tsv',errors)
    write(ROOT/f'results/{cohort}_depth_baseline.tsv',depth)
    print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
