#!/usr/bin/env python3
"""Secondary fairness diagnostic: same pilot-only scalar recipe as other CN baselines.

The native upstream calls remain the primary published-tool comparison. This
supplement scales the continuous C4Investigator total-C4 estimate only. Its
recipe is recorded before inspecting any held-out read predictions.
"""
import csv,datetime,hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def main():
    plan=json.loads((ROOT/'source/comparator_calibration_plan.json').read_text())
    assert hashlib.sha256(Path(__file__).read_bytes()).hexdigest()==plan['script_sha256']
    catalogue=[r for r in read(ROOT.parent/'hla-structural/source/catalogue.tsv') if r['locus']=='RCCX']
    records=[]
    for donor in (ROOT/'source/development_donors.txt').read_text().splitlines():
        p=ROOT/f'source/baseline/development/{donor}/results/C4Investigator_c4_detailed.csv'
        r=list(csv.DictReader(open(p)));assert len(r)==1
        estimate=float(r[0]['wgs_c4_copy']);truth=sum(int(r['copy_number']) for r in catalogue if r['donor']==donor)
        assert np.isfinite(estimate) and estimate>0 and truth>0
        records.append(dict(donor=donor,estimate=estimate,truth=truth))
    factor=float(np.median([r['estimate']/r['truth'] for r in records]));out=ROOT/'source/comparator_total_calibration.json'
    calibration=dict(factor=factor,pilot_records=records,recipe='median(continuous C4Investigator total estimate / assembly total) on the original 18 pilot donors only',plan_utc=plan['utc'])
    if out.exists():assert json.loads(out.read_text())==calibration
    else:out.write_text(json.dumps(calibration,indent=2)+'\n')
    validation=[]
    for r in read(ROOT/'source/validation_donors.tsv'):
        donor=r['donor'];p=ROOT/f'source/baseline/validation/{donor}/results/C4Investigator_c4_detailed.csv'
        rs=list(csv.DictReader(open(p)));assert len(rs)==1
        estimate=float(rs[0]['wgs_c4_copy'])/factor;truth=sum(int(x['copy_number']) for x in catalogue if x['donor']==donor)
        call=int(np.floor(estimate+.5)) if np.isfinite(estimate) else None
        validation.append(dict(donor=donor,calibrated_estimate=estimate,call=call,truth=truth,correct=int(call==truth)))
    with open(ROOT/'results/validation_C4Investigator_calibrated_total.tsv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(validation[0]),delimiter='\t');w.writeheader();w.writerows(validation)
    n=len(validation);correct=sum(r['correct'] for r in validation)
    text=f'''# Secondary C4Investigator total-dosage calibration

The unmodified C4Investigator calls remain in the primary report. This fairness diagnostic gives its continuous total-C4 estimate the same pilot-only scalar-calibration recipe used for the depth and marker baselines.

The recipe was fixed at {plan['utc']}, before inspecting held-out read predictions. The fitted factor is {factor:.6g}, using the original 18 pilot donors only. No held-out outcomes enter the fit. This supplement does not recalibrate A/B or long/short component calls or alter the frozen targeted assay.

Calibrated total-copy accuracy: **{correct}/{n}** attempted held-out donors. The per-donor estimates and calls are in `results/validation_C4Investigator_calibrated_total.tsv`.

This comparison is conditional on MHC-recruited reads and the recorded references and software environment. Differences from the native caller do not establish its performance on unrestricted whole-WGS inputs.
'''
    (ROOT/'COMPARATOR_CALIBRATION.md').write_text(text);print(correct,'/',n,'factor',factor)
if __name__=='__main__':main()
