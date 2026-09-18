# Matched T1K comparison and four-field performance

Supplemental analysis requested after the first six outcomes were available. Primary frozen scoring is unchanged. Figures are diploid genotype accuracy, not individual alleles.

## first_six: 6 donors

| Method | Completed | Two-field correct | Four-field correct | Four-field resolved |
|---|---:|---:|---:|---:|
| native | 6/6 | 42/48 | 12/33 | 27/33 |
| native-long | 6/6 | 43/48 | 12/33 | 27/33 |
| DogoHLA-no-graph | 6/6 | 46/48 | 14/33 | 27/33 |
| DogoHLA | 6/6 | 46/48 | 14/33 | 27/33 |
| T1K | 6/6 | 47/48 | Not resolved | 0/33 |

## completed_matched: 32 donors

| Method | Completed | Two-field correct | Four-field correct | Four-field resolved |
|---|---:|---:|---:|---:|
| native | 32/32 | 232/256 | 74/165 | 146/165 |
| native-long | 32/32 | 233/256 | 75/165 | 146/165 |
| DogoHLA-no-graph | 32/32 | 242/256 | 76/165 | 147/165 |
| DogoHLA | 32/32 | 242/256 | 76/165 | 147/165 |
| T1K | 32/32 | 252/256 | Not resolved | 0/165 |

Four-field truth requires an exact genomic match for both haplotypes. Unmatched/novel truth sequences are excluded for every method and are counted in the separate whole-gene endpoint. Shorter allele names are not silently expanded. T1K's archived output may stop at three fields; that is an output-resolution limit, not proof its reconstructed sequence is wrong.

T1K uses IPD 3.65 whereas native/DōgoHLA allele naming uses IPD 3.38. These are comparisons of the deployed configurations, not an isolated algorithm comparison. No matched OptiType, HLA*LA or other established short-read caller outputs were found in the local project.
