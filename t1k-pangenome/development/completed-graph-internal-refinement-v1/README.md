# Completed gene-internal refinement comparison

All 32 donor/condition jobs completed. Archived provenance was checked against
the submitted driver/model hashes and the internal-pairs parameter.

| Method | Four-field correct /47 |
|---|---:|
| Original T1K | 21 |
| Genomic IPD T1K | 41 |
| HPRC sampled / unsampled | 41 / 41 |
| HPRC + Asian sampled / unsampled | 40 / 43 |

The unsampled additive graph corrects DQA1 in HG03516 and NA19185 without losing
a correct genotype. It restores three correct alleles across these two genotypes.
HPRC unsampled restores one allele in NA19185 but leaves its diploid genotype
incorrect. The earlier erroneous unsampled HLA-A change no longer occurs. The
sampled additive DPB1 error remains.

This is a development result, after a revision motivated by development errors.
It is not independent evidence of superiority. All graph conditions preserve the
genomic-IPD two-field calls, 62/63 versus original T1K 63/63 on these donors.
Neither the primary candidate nor its validation outcome is frozen/scored yet.
The global informative-fragment count remains a known limitation: variation
across all candidates is not necessarily variation between the proposed/native
pair. Pair-specific discrimination needs checking before validation.
