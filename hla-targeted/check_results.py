#!/usr/bin/env python3
"""Completion invariants, separate from accuracy interpretation."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def main():
    freeze=json.loads((ROOT/'source/validation_freeze.json').read_text())
    for name,expected in freeze['files'].items():
        p=ROOT/name;assert p.exists(),name
        assert hashlib.sha256(p.read_bytes()).hexdigest()==expected,('changed after freeze',name)
    validation=read(ROOT/'source/validation_donors.tsv');donors={r['donor'] for r in validation};families={r['family'] for r in validation}
    development=set((ROOT/'source/all_development_donors.txt').read_text().splitlines())
    rows=json.loads((ROOT/'source/training_rows.json').read_text());assert len(rows)==543
    assert len(donors)==24 and len(development)==106 and not donors&development
    assert not families & {r['family'] for r in rows}
    priorfamilies=set(json.loads((ROOT/'source/validation_selection.json').read_text())['excluded_previous_families']);assert not families&priorfamilies
    assert not families & {r['family'] for r in read(ROOT/'source/training_C4_genes.tsv')}
    keys=(ROOT/'source/markers.txt').read_text().splitlines();assert len(keys)==len(set(keys))==12214
    assert [r.split()[0] for r in (ROOT/'source/marker_groups.tsv').read_text().splitlines()]==keys
    predictions=read(ROOT/'results/validation_dosage_predictions.tsv');assert len(predictions)==24 and {r['donor'] for r in predictions}==donors
    for r in predictions:
        d=r['donor'];raw=np.fromfile(ROOT/f'source/read_counts/validation/{d}.bin',dtype='<u4');assert len(raw)==len(keys) and raw.sum()>0
        groups=read(ROOT/f'source/read_counts/validation/{d}.groups.tsv');assert int(next(x['fragments'] for x in groups if x['metric']=='all_fragments'))>0
        baseline=ROOT/f'source/baseline/validation/{d}/results/C4Investigator_c4_summary.csv';assert baseline.exists(),d
        b=list(csv.DictReader(open(baseline)));assert len(b)==1 and b[0]['rn'].startswith(d+'_')
        if r['AB_status']=='call':assert sum(int(r[g+'_call']) for g in ['A','B','OTHER'])==int(r['C4_call'])
        if r['LS_status']=='call':assert sum(int(r[g+'_call']) for g in ['L','S'])==int(r['C4_call'])
    paths=read(ROOT/'results/validation_path_predictions.tsv');assert len(paths)==72
    for method in {r['method'] for r in paths}:
        ps=[r for r in paths if r['method']==method];assert len(ps)==24 and {r['donor'] for r in ps}==donors
        for r in ps:
            candidates=json.loads(r['dosage_compatible_pairs']);assert len(candidates)==int(r['dosage_compatible_pair_count'])
            if r['ranked_pair']:assert r['ranked_pair'].split(';') in candidates and r['status']=='imputed_not_physically_phased'
    metrics=read(ROOT/'results/validation_metrics.tsv');assert all(int(r['n'])==24 for r in metrics)
    assert all(0<=int(r['correct'])<=int(r['called'])<=int(r['n']) for r in metrics)
    audit=read(ROOT/'results/validation_C4_diagnostic_audit.tsv');assert {r['donor'] for r in audit}==donors
    assert all(r['status']!='alignment_incomplete' for r in audit)
    summary=dict(status='passed',validation_donors=24,development_donors=106,training_haplotypes=543,checks=['frozen_assay_hashes','family_exclusion','sample_split','complete_measurements','complete_dedicated_comparator','sum_preserving_dosage','ranked_pairs_in_compatible_set','all_attempted_denominators','validation_assembly_diagnostic_audit'],limitations='These checks establish provenance and accounting, not clinical truth, external-cohort validity, phase or calibrated uncertainty.')
    (ROOT/'results/completion_checks.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
