# Gourraud experimental-truth benchmark

946 donors; 660 family groups; graph-donor and known-family overlap zero. Four-field calls are retained, but this truth supports two-field scoring only.

| Method | Successful | Final failures | Two-field genotype pairs |
|---|---:|---:|---:|
| SpecHLA-IPD365-noFreq | 946/946 | 0 | 3093/4723 |
| T1K-four-field | 945/946 | 1 | 4437/4723 |
| full:DogoHLA | 945/946 | 1 | 3255/4723 |
| full:DogoHLA-no-graph | 945/946 | 1 | 3260/4723 |
| hprc:DogoHLA | 944/946 | 2 | 3246/4723 |
| hprc:DogoHLA-no-graph | 944/946 | 2 | 3245/4723 |
| asian_matched:DogoHLA | 945/946 | 1 | 3235/4723 |
| asian_matched:DogoHLA-no-graph | 945/946 | 1 | 3240/4723 |

Pending/failed calls retain planned eligible denominators. A final failure requires a hash-bound audit record after investigation/retry; it receives zero correct credit. Summary TSV includes ancestry strata. Intervals appear only after all planned samples have successful or final-failure outcomes in both compared methods.
