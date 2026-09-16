#!/usr/bin/env python3
"""Development-only calibrations for fragment-group dosage; no validation read access."""
import csv,json,re
from collections import defaultdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def labels():
    meta={r['gene_id']:r for r in read(ROOT/'source/training_C4_genes.tsv')}
    truth=defaultdict(lambda:dict(total=0,A=0,B=0,OTHER=0,L=0,S=0))
    for r in read(ROOT/'results/training_C4_diagnostic_audit.tsv'):
        d=truth[meta[r['gene_id']]['donor']]; f=r['inferred_form'];d['total']+=1
        assert len(f)==4 and f[-1] in 'LS',r
        d[f[2] if f[2] in 'AB' else 'OTHER']+=1;d[f[-1]]+=1
    return dict(truth)
def main():
    truth=labels();donors=(ROOT/'source/all_development_donors.txt').read_text().splitlines();records=[]
    for donor in donors:
        groups={r['metric']:int(r['fragments']) for r in read(ROOT/f'source/read_counts/development/{donor}.groups.tsv')}
        t=truth[donor]
        records.append(dict(donor=donor,**t,**groups))
    # Relative detection efficiency estimated on mixed-copy development donors only.
    # Ratios are per-fragment, avoiding overlap-probe double counting.
    def efficiency(g1,g2,t1,t2):
        ratios=[r[g1]/r[g2]*r[t2]/r[t1] for r in records if r[t1]>0 and r[t2]>0 and r[g1]>=10 and r[g2]>=10]
        assert len(ratios)>=10
        return dict(factor=float(np.median(ratios)),development_donors=len(ratios))
    for r in records:r['LONG_MEAN']=(r['LONG_LEFT']+r['LONG_RIGHT'])/2
    model=dict(A_relative_B=efficiency('A_only','B_only','A','B'),L_relative_S=efficiency('LONG_MEAN','SHORT','L','S'),OTHER_relative_B=1.0,minimum_diagnostic_fragments=20,maximum_rounding_distance=0.4,maximum_long_junction_imbalance=0.5,notes='Exploratory development fit. Thresholds are heuristic, not calibrated error probabilities. OTHER efficiency assumed, not fitted. Total-C4 factor inherited from the earlier 18-donor pilot.',total_C4_factor=json.loads((ROOT.parent/'hla-structural/source/dosage_calibration.json').read_text())['factor'])
    (ROOT/'source/dosage_model.json').write_text(json.dumps(model,indent=2)+'\n')
    with open(ROOT/'results/development_calibration_inputs.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
    print(json.dumps(model,indent=2))
if __name__=='__main__':main()
