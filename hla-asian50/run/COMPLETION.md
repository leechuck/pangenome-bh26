# Submitted calling run complete

The current 1000 Genomes run has finished: **20 East Asian and 20 South Asian donors**, with all 120 calling jobs completed successfully. The first all-complete monitoring snapshot was recorded at **2026-09-17 02:35:43 UTC**. The subsequent output-integrity job also completed successfully. Saudi/Arabian validation remains deferred.

| Output | Completed |
|---|---:|
| T1K donor jobs / native genotype files | 40/40 |
| SpecHLA donor jobs / native genotype files | 40/40 |
| PanGenie donor jobs, each running full and HPRC-only panels | 40/40 |
| PanGenie VCFs | 80/80 |
| Family-excluded panel indices | 10/10 |

[Slurm completion records](source/job_completion.tsv) retain job IDs, exit codes, start/end times and allocated CPUs. All 120 calling jobs exited `0:0`. SpecHLA concurrency was increased as resources became available (1 → 4 → 5 → 6 jobs), without changing any caller settings.

## Verification and saved outputs

The [integrity check](source/output_integrity.json) passed: expected donor identities and file schemas, diploid called-genotype allele indices, and nonempty VCFs were checked. All 80 VCFs were hashed. Missing genotypes remain missing; no prediction files were edited.

- [Native T1K and SpecHLA call archive](source/native_HLA_calls.tar.gz): 80 files, with [SHA256](source/native_HLA_calls.sha256).
- [PanGenie output manifest](source/pangenie_output_manifest.json): remote relative paths, sizes and SHA256 hashes for all 80 VCFs. Large VCFs remain under `/home/leechuck/hla/codex-asian50/run/` on DDBJ.
- [Per-donor output checks](source/output_integrity.tsv): record counts and missing-genotype counts, retaining East/South Asian strata.
- [Panel preparation audits](source/panel_audits): per-fold sample counts and site-filtering denominators.

Across the 40 donors, full-panel VCFs contain 3,250,704 records, including 63 missing genotypes; HPRC-only VCFs contain 3,504,632 records, including 120 missing genotypes. These are output counts, **not accuracy scores**. The panels retain different site sets because of training-panel missingness and allele filtering; these totals cannot rank methods or establish SV improvement. No-calls must remain in the appropriate attempted-site denominators during evaluation.

Two historical failed attempts are retained in the logs: the original graph reader could not read GBZ v2 and was replaced by the compatible build; the initial integrity checker rejected PanGenie's valid missing genotype `.` because it assumed two explicit fields. The checker was corrected and regression-tested to count `.` as a no-call while continuing to reject out-of-range alleles. Successful verification job: **20634221**. Neither fix modified predictions or required rerunning a completed caller.

## Accuracy evaluation follow-up

The call-generation milestone above preceded accuracy scoring. See the subsequent [accuracy comparison](../accuracy/REPORT.md) for held-out assembly concordance, paired donor comparisons, coverage losses, the published linear SNV baseline and experimental specialist HLA scores. Independent SNV/SV truth and callable masks, a harmonized linear SV baseline, and classical-HLA translation of graph predictions remain outstanding. Native HLA databases differ (SpecHLA outputs identify IPD-IMGT/HLA 3.38.0); database provenance must be considered when interpreting comparisons. Initial reads were recruited through reference alignments, so this run does not establish whole-WGS recruitment improvements.
