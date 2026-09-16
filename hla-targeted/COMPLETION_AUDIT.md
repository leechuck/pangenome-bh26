# Completion audit

The following requirements were checked against the recorded artifacts. This audit establishes completion of this experiment, not clinical accuracy or a global novelty claim.

| Requirement | Evidence | Disposition |
|---|---|---|
| A/B, long/short and module-context typing | Frozen probe definitions, dosage model and validation predictions | Complete |
| Development-only use of prior 106 donors | Selection manifest, calibration inputs and family-exclusion checks | Complete |
| New read validation | 24 new donors; selected families absent from 543 reference paths; 24 complete measurement sets | Complete within the existing panel; not external-cohort validation |
| Preserve ambiguity and distinguish phase | Complete compatible-pair sets, imputation labels, structural error counterexamples | Complete; no unique dosage-only structure in this set |
| Dedicated comparator | Native outputs for 18 development and 24 validation donors; independently specified total-dosage calibration | Complete |
| Module-context contribution | Frozen ablation evaluated on all 24 donors | Complete; no accuracy gain |
| Sequence-label and probe controls | 91-gene diagnostic audit, complete probe coverage, 48-haplotype whole-MHC specificity audit | Complete within MHC; whole-genome specificity unproven |
| Reproducibility | Frozen file hashes, training-supported inherited-feature audit, software manifests, 934-file input archive, byte-identical replay of 24 dosage and 72 path rows | Complete for derived-input inference; raw-read and sequence reconstruction require original inputs |
| Reporting | Main report, calibration/configuration supplements, error/candidate records and standalone figures | Complete |

`results/completion_checks.json` verifies the primary invariants, while `results/replay_check.json` and `results/validation_MHC_specificity.json` record the additional controls. The named reports explain the limits of each check. UK Biobank phenotype validation remains a later study, as agreed at the outset; no disease association or validated novel named allele is claimed here.
