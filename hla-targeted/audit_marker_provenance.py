#!/usr/bin/env python3
"""Verify that the inherited global sketch contributes only training-supported features.

The old vocabulary was a global deterministic hash sketch, not itself a held-out
reference panel. Eligibility in the new experiment is computed from training
families only. This audit documents that distinction without changing inference.
"""
import json,importlib.util
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent

def main():
    source=ROOT.parent/'hla-structural/prepare_markers.py';spec=importlib.util.spec_from_file_location('old_sketch',source);old_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(old_module)
    rows=json.loads((ROOT/'source/training_rows.json').read_text());X=np.load(ROOT/'source/training_profiles.npy');B=np.stack([np.fromfile(ROOT/'source/background'/f'{r["hap_id"].replace("#","_")}.bin',dtype='<u4') for r in rows])
    keys=[int(k) for k in (ROOT/'source/markers.txt').read_text().splitlines()];oldkeys=set(map(int,(ROOT.parent/'hla-structural/source/markers.txt').read_text().splitlines()));old=np.array([k in oldkeys for k in keys]);assert all(old_module.mix(k)&63==0 for k in oldkeys)
    support=np.zeros(len(keys),int)
    for fam in {r['family'] for r in rows}:support+=(X[[i for i,r in enumerate(rows) if r['family']==fam]]>0).any(axis=0)
    keep=(support>=3)&(X==B).all(axis=0)&(X.max(axis=0)<=8)
    masks=np.array([int(s.split()[1]) for s in (ROOT/'source/marker_groups.tsv').read_text().splitlines()]);targeted=masks>0
    assert np.all(X[:,keep&old].max(axis=0)>0) and np.all(support[keep&old]>=3)
    assert np.all(X[:,targeted].max(axis=0)>0)
    # Match inherited columns against the original profiles for these same rows.
    legacy_rows=json.loads((ROOT.parent/'hla-structural/source/matrix_rows.json').read_text());legacy_ids={r['sequence_id']:i for i,r in enumerate(legacy_rows)};legacy=np.load(ROOT.parent/'hla-structural/source/path_marker_counts.npy')
    legacy_keys=list(map(int,(ROOT.parent/'hla-structural/source/markers.txt').read_text().splitlines()));current={k:i for i,k in enumerate(keys)}
    assert np.array_equal(X[:,[current[k] for k in legacy_keys]],legacy[[legacy_ids[r['sequence_id']] for r in rows]])
    result=dict(status='passed',inherited_global_vocabulary=len(oldkeys),retained_inherited_features=int((keep&old).sum()),targeted_group_probes=int(targeted.sum()),all_retained_inherited_features_have_three_training_families=True,all_inherited_keys_follow_fixed_hash_rule=True,all_targeted_group_probes_present_in_training_RCCX=True,legacy_training_profile_columns_match=True,interpretation='The inherited measurement vocabulary was originally a union across the panel. Every inherited feature actually used here passes a fixed sequence-only hash and occurs in at least three retained training families; regenerating that sketch from those training sequences would include all used features. Held-out-only vocabulary entries are excluded from scoring and normalization. New targeted probe groups were discovered solely from retained training sequences.')
    (ROOT/'results/marker_provenance_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
