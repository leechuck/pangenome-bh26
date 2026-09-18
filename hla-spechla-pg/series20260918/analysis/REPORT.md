# Matched 64-donor benchmark status

All methods use IPD 3.65. Four-field truth requires exact genomic identity for both haplotypes; two-field truth requires exact CDS identity. These are assembly-derived labels. Independent Gourraud concordance is reported separately.

| Method | Complete | Two-field genotype pairs | Four-field genotype pairs |
|---|---:|---:|---:|
| SpecHLA-IPD365-noFreq | 64/64 | 377/509 | 103/323 |
| SpecHLA-IPD365-long-noFreq | 64/64 | 379/509 | 102/323 |
| T1K-four-field | 64/64 | 495/509 | 163/323 |
| full:DogoHLA | 64/64 | 404/509 | 108/323 |
| full:DogoHLA-no-graph | 64/64 | 402/509 | 108/323 |
| hprc:DogoHLA | 64/64 | 407/509 | 109/323 |
| hprc:DogoHLA-no-graph | 64/64 | 402/509 | 109/323 |
| asian_matched:DogoHLA | 64/64 | 400/509 | 108/323 |
| asian_matched:DogoHLA-no-graph | 64/64 | 399/509 | 108/323 |

`summary.tsv` also separates EAS, SAS, EUR and AFR. Pending methods retain all planned truth denominators in TSV output; headline accuracy is withheld until the method finishes all 64 donors. Calls on eligible truth that remain unresolved are failures, not removed from the denominator.
