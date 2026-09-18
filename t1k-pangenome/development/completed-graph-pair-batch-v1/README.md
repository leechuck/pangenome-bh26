# Completed eight-donor graph pair refinement

The fixed development panel contains two donors each from EAS, SAS, EUR and AFR.
All four graph conditions completed for all eight donors. Four-field correct
genotypes out of 47 eligible comparisons:

| Method | Correct |
|---|---:|
| Original T1K | 21 |
| Genomic IPD T1K | 41 |
| HPRC sampled / unsampled | 41 / 41 |
| HPRC + Asian sampled / unsampled | 40 / 40 |

HPRC makes no call changes. Each additive condition changes one previously
correct call incorrectly, at different donor/locus pairs: sampled HG02132 DPB1
and unsampled HG01530 A. Neither rescues any baseline error. These small
development differences do not establish population-level superiority or
inferiority; they provide no evidence of additional graph benefit. The large
improvement over original T1K belongs to genomic IPD, not graph refinement.

All refinements preserve the genomic-IPD two-field calls (62/63 versus original
T1K 63/63). Native alleles missing from the graph, unresolved native calls, ties
and unknown winning paths trigger fallback. The decision audit records these
reasons separately from accuracy. Independent validation outcomes remain
unexamined. This version is not frozen as the primary validation candidate.
