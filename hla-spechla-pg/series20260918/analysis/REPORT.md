# Matched 64-donor benchmark status

All methods use IPD 3.65. Four-field truth requires exact genomic identity for both haplotypes; two-field truth requires exact CDS identity. These are assembly-derived labels. Independent Gourraud concordance is reported separately.

| Method | Complete | Two-field genotype pairs | Four-field genotype pairs |
|---|---:|---:|---:|
| SpecHLA-IPD365-noFreq | 0/64 | pending | pending |
| SpecHLA-IPD365-long-noFreq | 0/64 | pending | pending |
| T1K-four-field | 64/64 | 495/509 | 163/323 |
| full:DogoHLA | 0/64 | pending | pending |
| full:DogoHLA-no-graph | 0/64 | pending | pending |
| hprc:DogoHLA | 0/64 | pending | pending |
| hprc:DogoHLA-no-graph | 0/64 | pending | pending |
| asian_matched:DogoHLA | 0/64 | pending | pending |
| asian_matched:DogoHLA-no-graph | 0/64 | pending | pending |

`summary.tsv` also separates EAS, SAS, EUR and AFR. Pending methods retain all planned truth denominators in TSV output; headline accuracy is withheld until the method finishes all 64 donors. Calls on eligible truth that remain unresolved are failures, not removed from the denominator.
