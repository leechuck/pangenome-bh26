#!/usr/bin/env python3
"""Freeze assay after development and before inspecting held-out predictions."""
import csv,datetime,hashlib,json,subprocess,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def main():
    out=ROOT/'source/validation_freeze.json'
    assert not out.exists(),'Freeze already exists; preserve it and record subsequent changes separately'
    validation=read(ROOT/'source/validation_donors.tsv');families={r['family'] for r in validation}
    rows=json.loads((ROOT/'source/training_rows.json').read_text())
    assert not families & {r['family'] for r in rows}
    development=set((ROOT/'source/all_development_donors.txt').read_text().splitlines())
    assert not development & {r['donor'] for r in validation}
    assert len(read(ROOT/'results/development_dosage_predictions.tsv'))==len(development)==106
    assert not (ROOT/'results/validation_dosage_predictions.tsv').exists()
    assert json.loads((ROOT/'results/probe_specificity.json').read_text())['off_locus_probes']==0
    subprocess.run(['python3',str(ROOT/'test_count_groups.py')],check=True)
    shutil.copyfile(ROOT/'results/development_metrics.tsv',ROOT/'source/frozen_development_metrics.tsv')
    paths=list(ROOT.glob('*.py'))+list(ROOT.glob('*.cpp'))+list(ROOT.glob('*.sbatch'))+list(ROOT.glob('*.sh'))
    paths += [ROOT/'source'/p for p in ['markers.txt','marker_groups.tsv','marker_design.json','target_definitions.json','training_rows.json','training_profiles.npy','dosage_model.json','tool_version.json','validation_donors.tsv','validation_selection.json','training_C4_genes.tsv','frozen_development_metrics.tsv']]
    paths += [ROOT/'results'/p for p in ['training_C4_diagnostic_audit.tsv','development_calibration_inputs.tsv','development_dosage_predictions.tsv','probe_specificity.json']]
    paths += list((ROOT/'source/background').glob('*.bin'))
    assert len(list((ROOT/'source/background').glob('*.bin')))==len(rows)
    paths += [ROOT.parent/'hla-structural/source/markers.txt',ROOT.parent/'hla-structural/source/dosage_calibration.json',ROOT.parent/'hla-structural/results/calibrated_depth_baseline.tsv',ROOT.parent/'hla-analysis/source/rccx_depth_targets.bed']
    hashes={str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else '../'+str(p.relative_to(ROOT.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}
    record=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),validation_donors=len(validation),development_donors=len(development),excluded_families=sorted(families),training_haplotypes=len(rows),scope='All selected families excluded from assay discovery and candidate reference paths. Validation samples remain part of the original graph construction; this is new-donor read validation within the existing panel, not external-cohort validation.',files=hashes)
    out.write_text(json.dumps(record,indent=2)+'\n');print('Frozen',len(hashes),'files at',record['utc'])
if __name__=='__main__':main()
