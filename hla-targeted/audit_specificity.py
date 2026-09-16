#!/usr/bin/env python3
"""Compare every probe's full-MHC count with its RCCX interval count."""
import csv,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def main():
    rows=json.loads((ROOT/'source/training_rows.json').read_text())
    X=np.load(ROOT/'source/training_profiles.npy')
    B=np.stack([np.fromfile(ROOT/'source/background'/f'{r["hap_id"].replace("#","_")}.bin',dtype='<u4') for r in rows])
    assert B.shape==X.shape
    assert np.all(B>=X),'Interval counts exceed whole-assembly counts; check sources/order'
    probes=[list(map(int,s.split())) for s in (ROOT/'source/marker_groups.tsv').read_text().splitlines()]
    groups=json.loads((ROOT/'source/marker_design.json').read_text())['groups']
    records=[]
    for j,(key,mask) in enumerate(probes):
        if not mask:continue
        records.append(dict(marker=key,groups=';'.join(g for g,b in groups.items() if mask&b),off_locus_haplotypes=int((B[:,j]>X[:,j]).sum()),maximum_extra_copies=int((B[:,j]-X[:,j]).max()),present_haplotypes=int((X[:,j]>0).sum())))
    with open(ROOT/'results/probe_specificity.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
    summary=dict(training_haplotypes=len(rows),targeted_probes=len(records),off_locus_probes=sum(r['off_locus_haplotypes']>0 for r in records),scope='MHC sequences only; does not prove whole-genome specificity',groups={g:dict(probes=sum(bool(m&b) for _,m in probes),off_locus_probes=sum(bool(m&b) and bool((B[:,j]>X[:,j]).any()) for j,(_,m) in enumerate(probes))) for g,b in groups.items()})
    (ROOT/'results/probe_specificity.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
