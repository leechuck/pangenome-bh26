# Gourraud experimental-truth benchmark

946 donors; 660 family groups; graph-donor and known-family overlap zero. Four-field calls are retained, but this truth supports two-field scoring only.

| Method | Successful | Final failures | Two-field genotype pairs |
|---|---:|---:|---:|
| SpecHLA-IPD365-noFreq | 0/946 | 0 | pending |
| T1K-four-field | 945/946 | 1 | 4437/4723 |
| full:DogoHLA | 0/946 | 0 | pending |
| full:DogoHLA-no-graph | 0/946 | 0 | pending |
| hprc:DogoHLA | 0/946 | 0 | pending |
| hprc:DogoHLA-no-graph | 0/946 | 0 | pending |
| asian_matched:DogoHLA | 0/946 | 0 | pending |
| asian_matched:DogoHLA-no-graph | 0/946 | 0 | pending |

Pending/failed calls retain planned eligible denominators. A final failure requires a hash-bound audit record after investigation/retry; it receives zero correct credit. Summary TSV includes ancestry strata. Intervals appear only after all planned samples have successful or final-failure outcomes in both compared methods.
