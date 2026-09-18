# Benchmark and T1K-improvement completion audit

Audited against the active objective: complete the started benchmarks, execute the
planned T1K improvement experiments, and establish a provable improvement over
original T1K. The supported claim is **four-field improvement under the specified
reference/inference configuration**, with preserved two-field predictions. This
is not evidence that Asian graph paths add accuracy beyond genomic-IPD T1K.

| Requirement | Authoritative evidence | Outcome |
|---|---|---|
| Matched Asian/non-Asian cohort and baselines | `../hla-spechla-pg/series20260918/cohort.tsv`, `analysis/REPORT.md`, ancestry-stratified `analysis/summary.tsv` | 16 EAS + 16 SAS versus 16 EUR + 16 AFR; all nine method/ablation entries complete 64/64 |
| Expanded experimental-truth benchmark | `../hla-spechla-pg/series20260918/analysis/gourraud/REPORT.md` and accompanying tables/audits | 946 donors, 660 families; all methods terminal; investigated final failures retained in denominator |
| T1K four-field output and matched IPD database | `../hla-analysis/source/imgt/release_version.txt`, baseline/reference plans, result tables | IPD 3.65.0 across comparisons; explicit four-field T1K output; no padding of short names |
| HPRC versus Asian/additive panels | Matched/Gourraud panel ablations; `development/completed-reference-controls-v2/`; HGSVC frozen execution | Replacement and additive/reference controls completed; HPRC and HPRC+Asian graph inference executed independently |
| T1K graph/refinement experiments | Development sampling, recruitment, linear-reference and pair-refinement artifacts; original 84-donor report | Experiments completed; original 84-donor combined gate failed and is preserved; exposed-sample changes are exploratory |
| New independent cohort and graph exclusions | `validation-next/RESERVATION.json`, exposure audits, `GRAPH_SOURCE_AUDIT.json` | All 28 reserved donors retained; 3,515 actual source entries checked; zero donor/known-family overlap; missing pedigree limitations stated |
| Fixed follow-up candidate and endpoint | `validation-next/PRIMARY_CANDIDATE.json`, `FROZEN_VALIDATION.json`, `EVALUATION_FREEZE.json`, `HANDOFF_FREEZE.json` | Candidate, inference, eligibility and endpoint rules fixed before outcome evaluation |
| Complete execution and provenance | `validation-next/EXECUTION_COMPLETE.json`, `execution-handoff.json` | All 1,680 inference tasks completed 0:0; all 196 original/intermediate/final output records verified; no failed donor omitted |
| Honest scoring correction | `validation-next/SCORING_CORRECTION.json`, `SCORING_CORRECTION_VERIFIED.json`, retained original report | Post-unblinding two-field name-format bug corrected; predictions, eligibility and every four-field row unchanged; correction explicitly disclosed |
| Improvement over original T1K | `validation-next/completed-scoring-correction-v2/report.json`, `summary.tsv`, `gene_scores.tsv` | Four-field 120/138 versus 71/138; +35.51 points, paired family-bootstrap 95% CI +28.46 to +42.64; 51 gains, 2 losses |
| Two-field preservation | Same report; recomputed anchoring of verified snapshot | 142/144 for all methods; exact per-gene two-field prediction identity; all three predeclared acceptance checks pass |
| Attribution to Asian graph | Same report and target-locus call comparison | No incremental benefit: anchored genomic-IPD, HPRC and HPRC+Asian calls are identical across all eight loci on this cohort |
| Reviewable outputs | Corrected `REPORT.md`, `calls.tsv`, `truth_qc.tsv`, `discordant_genotypes.tsv`; executable scorer and correction tests | Two- and four-field calls exported; full eligibility and discordance records retained; 13 HGSVC tests passed |
| HLA-HD | User explicitly postponed indefinitely | Excluded from completion requirements; no invented results |

The graph comparison is complete but negative for incremental benefit. The positive
result belongs to the full-genomic-IPD/coarse-anchor enhancement over the original
T1K DNA-reference configuration. The selected graph candidate shares that gain.
Do not claim an independent graph-specific improvement, clinical validation, or
uniform superiority by ancestry/locus. Four-field truth is assembly-derived and
eligible for 138/224 genotypes; ineligible truth receives no correctness credit.
The original 84-donor failure remains part of the experimental history.

The scheduler monitor exits after successful anchoring. Its absence as a transient
systemd unit is not the completion evidence; the task-level Slurm audit and verified
output chains are. A late scheduling-helper error was reconciled using its remote
incremental record and all 896 completed mapping/projection tasks, without retries.

All result paths in this audit are relative to `t1k-pangenome/` unless explicitly
prefixed otherwise. The matched benchmark's `analysis/` paths are relative to its
`series20260918/` directory. Large read/alignment files stay on IBEX; the local
verified prediction snapshot is under `work/hgsvc-validation-snapshot/`.
