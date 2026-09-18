# Independent HGSVC28 validation — corrected scoring

A documented post-unblinding correction removes the gene prefix from two-field truth to match the existing prediction parser. Predictions, eligibility, four-field scoring and statistical gates are unchanged. The original erroneous report is preserved; see ../SCORING_CORRECTION.json.

All 28 reserved donors; four methods; eight loci. Four-field truth requires exact whole-genomic matches. Two-field truth is also genomic-derived. Original two-field prediction identity was verified for each anchored method.

| Method | Fields | Correct | Eligible | Failed donors |
|---|---:|---:|---:|---:|
| T1K | 2 | 142 | 144 | 0 |
| T1K | 4 | 71 | 138 | 0 |
| ipd_genome | 2 | 142 | 144 | 0 |
| ipd_genome | 4 | 120 | 138 | 0 |
| graph_hprc | 2 | 142 | 144 | 0 |
| graph_hprc | 4 | 120 | 138 | 0 |
| graph_hprc_asian | 2 | 142 | 144 | 0 |
| graph_hprc_asian | 4 | 120 | 138 | 0 |

Primary combined protocol pass: **True**.

The original-T1K contrast measures the complete pipeline change. Incremental genomic-IPD and HPRC contrasts are exploratory; reference-only gains do not prove an Asian graph advantage. Per-ancestry intervals, failure accounting and paired contrasts are in report.json. The earlier 84-donor failed combined gate remains a separate result.

## Effect and attribution

The four-field gain over original T1K is **35.51 percentage points** (paired
family-bootstrap 95% interval **28.46 to 42.64**), with **51 gains and 2 losses**.
Two-field accuracy is **142/144 for every method**, and the two-field prediction
for every target gene was verified identical to original T1K. All 28 donors and
all four methods completed successfully.

Anchored genomic-IPD, HPRC and HPRC+Asian have **identical calls across all eight
target loci**, including ineligible loci—not merely equal accuracy totals. This
cohort therefore demonstrates the full-genomic-reference/coarse-anchor improvement
over original T1K, but no incremental benefit from graph refinement or Asian paths.
The selected graph mode is unsampled; this result is not evidence that k-mer path
selection improves HLA typing.

## Scope and limitations

There are 28 donors in 28 known families: 5 EAS, 8 SAS, 11 AFR, 2 EUR and 2 AMR.
All eight loci are reported. Four-field truth is eligible for 138/224 genotypes;
86 remain ineligible rather than receiving correctness credit. Two-field truth is
eligible for 144/224 and is genomic-derived in this follow-up. Across 448 truth
haplotypes, 339 have exact whole-gene IPD matches and 109 do not; unmatched does
not establish novelty. Truth is assembly-derived, not clinical four-field typing.

The graph audit found no reserved-donor or known-family overlaps in 3,515 source
entries. Pedigrees are unavailable for some graph donors, so unknown relationships
cannot be ruled out. The small, ancestry-imbalanced cohort does not establish an
advantage for every population or locus.

The two-field name-format correction was made after unblinding and is explicitly
recorded. It changes no predictions, eligibility, four-field scores, bootstrap
settings or acceptance gates. The original erroneous evaluation is preserved in
`../completed-frozen-v1/`; use this corrected report for interpretation. Exact
row-by-row invariance checks are in `../SCORING_CORRECTION_VERIFIED.json`.
