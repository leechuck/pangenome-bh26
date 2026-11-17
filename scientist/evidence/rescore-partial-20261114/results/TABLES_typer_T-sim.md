# Mosaic-typer test, read set T-sim (all 8 genes), partial-credit re-scoring

975 true copies per arm.

**all copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 975 | 2 | 757 (77.6%) | 924 (94.8%) | 868 (89.0%) | 934 (95.8%) | 963 (98.8%) | 960 (98.5%) | 966 (99.1%) | 969 (99.4%) |
| B | 975 | 2 | 784 (80.4%) | 923 (94.7%) | 884 (90.7%) | 939 (96.3%) | 963 (98.8%) | 961 (98.6%) | 966 (99.1%) | 969 (99.4%) |
| C | 975 | 2 | 762 (78.2%) | 924 (94.8%) | 869 (89.1%) | 934 (95.8%) | 963 (98.8%) | 960 (98.5%) | 966 (99.1%) | 969 (99.4%) |
| D | 975 | 2 | 790 (81.0%) | 924 (94.8%) | 885 (90.8%) | 939 (96.3%) | 963 (98.8%) | 961 (98.6%) | 966 (99.1%) | 969 (99.4%) |
| E | 975 | 2 | 849 (87.1%) | 940 (96.4%) | 933 (95.7%) | 959 (98.4%) | 973 (99.8%) | 967 (99.2%) | 972 (99.7%) | 973 (99.8%) |

**all copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 973 | 5.3 / 0 | 0.3 / 0 | 5.1 / 0 | 0 |
| B | 973 | 5.2 / 0 | 0.3 / 0 | 5.0 / 0 | 0 |
| C | 973 | 5.3 / 0 | 0.3 / 0 | 5.0 / 0 | 0 |
| D | 973 | 5.2 / 0 | 0.3 / 0 | 5.0 / 0 | 0 |
| E | 973 | 0.3 / 0 | 0.1 / 0 | 0.2 / 0 | 0 |

**IPD-absent copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 162 | 1 | 34 (21.0%) | 142 (87.7%) | 77 (47.5%) | 127 (78.4%) | 152 (93.8%) | 152 (93.8%) | 155 (95.7%) | 158 (97.5%) |
| B | 162 | 1 | 62 (38.3%) | 142 (87.7%) | 94 (58.0%) | 131 (80.9%) | 152 (93.8%) | 152 (93.8%) | 155 (95.7%) | 158 (97.5%) |
| C | 162 | 1 | 40 (24.7%) | 143 (88.3%) | 78 (48.1%) | 127 (78.4%) | 152 (93.8%) | 152 (93.8%) | 155 (95.7%) | 158 (97.5%) |
| D | 162 | 1 | 67 (41.4%) | 143 (88.3%) | 94 (58.0%) | 131 (80.9%) | 152 (93.8%) | 152 (93.8%) | 155 (95.7%) | 158 (97.5%) |
| E | 162 | 0 | 129 (79.6%) | 158 (97.5%) | 143 (88.3%) | 153 (94.4%) | 162 (100.0%) | 160 (98.8%) | 162 (100.0%) | 162 (100.0%) |

**IPD-absent copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 161 | 30.7 / 2 | 1.1 / 0 | 29.6 / 2 | 0 |
| B | 161 | 30.2 / 1 | 1.1 / 0 | 29.2 / 1 | 0 |
| C | 161 | 30.6 / 2 | 1.1 / 0 | 29.5 / 2 | 0 |
| D | 161 | 30.2 / 1 | 1.1 / 0 | 29.1 / 1 | 0 |
| E | 162 | 0.8 / 0 | 0.1 / 0 | 0.7 / 0 | 0 |

**IPD-present copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 813 | 1 | 723 (88.9%) | 782 (96.2%) | 791 (97.3%) | 807 (99.3%) | 811 (99.8%) | 808 (99.4%) | 811 (99.8%) | 811 (99.8%) |
| B | 813 | 1 | 722 (88.8%) | 781 (96.1%) | 790 (97.2%) | 808 (99.4%) | 811 (99.8%) | 809 (99.5%) | 811 (99.8%) | 811 (99.8%) |
| C | 813 | 1 | 722 (88.8%) | 781 (96.1%) | 791 (97.3%) | 807 (99.3%) | 811 (99.8%) | 808 (99.4%) | 811 (99.8%) | 811 (99.8%) |
| D | 813 | 1 | 723 (88.9%) | 781 (96.1%) | 791 (97.3%) | 808 (99.4%) | 811 (99.8%) | 809 (99.5%) | 811 (99.8%) | 811 (99.8%) |
| E | 813 | 2 | 720 (88.6%) | 782 (96.2%) | 790 (97.2%) | 806 (99.1%) | 811 (99.8%) | 807 (99.3%) | 810 (99.6%) | 811 (99.8%) |

**IPD-present copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 812 | 0.3 / 0 | 0.1 / 0 | 0.2 / 0 | 0 |
| B | 812 | 0.3 / 0 | 0.1 / 0 | 0.2 / 0 | 0 |
| C | 812 | 0.3 / 0 | 0.1 / 0 | 0.2 / 0 | 0 |
| D | 812 | 0.3 / 0 | 0.1 / 0 | 0.2 / 0 | 0 |
| E | 811 | 0.2 / 0 | 0.1 / 0 | 0.1 / 0 | 0 |

**By gene**

| gene | A: copies; gb exact; mean gb / mean cod | B: copies; gb exact; mean gb / mean cod | C: copies; gb exact; mean gb / mean cod | D: copies; gb exact; mean gb / mean cod | E: copies; gb exact; mean gb / mean cod |
|---|---|---|---|---|---|
| A | 122; 113; 0.8 / 0.4 | 122; 113; 0.8 / 0.4 | 122; 113; 0.8 / 0.4 | 122; 113; 0.8 / 0.4 | 122; 113; 0.1 / 0.0 |
| B | 122; 117; 0.0 / 0.0 | 122; 117; 0.0 / 0.0 | 122; 117; 0.0 / 0.0 | 122; 117; 0.0 / 0.0 | 122; 119; 0.0 / 0.0 |
| C | 122; 112; 1.0 / 0.4 | 122; 113; 1.0 / 0.4 | 122; 113; 1.0 / 0.4 | 122; 114; 1.0 / 0.4 | 122; 113; 0.2 / 0.1 |
| DPA1 | 122; 89; 0.8 / 0.3 | 122; 90; 0.8 / 0.3 | 122; 90; 0.8 / 0.3 | 122; 91; 0.8 / 0.3 | 122; 92; 0.3 / 0.1 |
| DPB1 | 122; 102; 0.2 / 0.0 | 122; 99; 0.5 / 0.0 | 122; 102; 0.2 / 0.0 | 122; 99; 0.5 / 0.0 | 122; 116; 0.1 / 0.1 |
| DQA1 | 122; 85; 0.7 / 0.0 | 122; 97; 0.4 / 0.0 | 122; 88; 0.6 / 0.0 | 122; 101; 0.4 / 0.0 | 122; 91; 0.5 / 0.0 |
| DQB1 | 121; 103; 7.7 / 0.5 | 121; 103; 7.7 / 0.5 | 121; 103; 7.7 / 0.5 | 121; 103; 7.7 / 0.5 | 121; 107; 0.1 / 0.0 |
| DRB1 | 122; 36; 31.1 / 0.3 | 122; 52; 30.6 / 0.3 | 122; 36; 31.1 / 0.3 | 122; 52; 30.6 / 0.3 | 122; 98; 1.0 / 0.0 |

**Sensitivity to the tie rule**

| arm | copies with tied alleles | gb exact (worst tied allele; primary) | gb exact (best tied allele) | mean gb wrong bases (worst) | mean gb wrong bases (best) |
|---|---|---|---|---|---|
| A | 11 | 757 (77.6%) | 757 (77.6%) | 5.31 | 5.31 |
| B | 12 | 784 (80.4%) | 784 (80.4%) | 5.24 | 5.21 |
| C | 11 | 762 (78.2%) | 762 (78.2%) | 5.30 | 5.30 |
| D | 11 | 790 (81.0%) | 790 (81.0%) | 5.24 | 5.21 |
| E | 9 | 849 (87.1%) | 849 (87.1%) | 0.29 | 0.29 |

**Paired contrasts**

| contrast | subset | metric | n paired | values (a vs b) | difference | 95% CI (family bootstrap) | sign-flip p |
|---|---|---|---|---|---|---|---|
| B - A | all copies | share gene body exact (pp) | 975 | 784 vs 757 (+34/-7) | +2.77 pp | [+1.49, +3.90] | 0.00035 |
| B - A | all copies | share coding exact (pp) | 975 | 923 vs 924 (+5/-6) | -0.10 pp | [-0.94, +0.59] | 1 |
| B - A | all copies | mean wrong bases, gene body | 973 | 5.24 vs 5.31 (closer 54/farther 18) | -0.06 | [-0.13, +0.03] | 0.182 |
| B - A | all copies | mean wrong bases, coding | 973 | 0.25 vs 0.25 (closer 5/farther 6) | +0.00 | [-0.01, +0.01] | 1 |
| B - A | IPD-absent | share gene body exact (pp) | 162 | 62 vs 34 (+30/-2) | +17.28 pp | [+11.46, +22.93] | 5e-05 |
| B - A | IPD-absent | share coding exact (pp) | 162 | 142 vs 142 (+2/-2) | +0.00 pp | [-2.44, +2.52] | 1 |
| B - A | IPD-absent | mean wrong bases, gene body | 161 | 30.24 vs 30.66 (closer 45/farther 11) | -0.42 | [-0.74, -0.01] | 0.054 |
| B - A | IPD-absent | mean wrong bases, coding | 161 | 1.06 vs 1.06 (closer 2/farther 2) | +0.00 | [-0.03, +0.02] | 1 |
| B - A | IPD-present | share gene body exact (pp) | 813 | 722 vs 723 (+4/-5) | -0.12 pp | [-0.87, +0.63] | 1 |
| B - A | IPD-present | share coding exact (pp) | 813 | 781 vs 782 (+3/-4) | -0.12 pp | [-0.87, +0.45] | 1 |
| B - A | IPD-present | mean wrong bases, gene body | 812 | 0.29 vs 0.28 (closer 9/farther 7) | +0.01 | [-0.03, +0.05] | 0.891 |
| B - A | IPD-present | mean wrong bases, coding | 812 | 0.09 vs 0.09 (closer 3/farther 4) | +0.00 | [-0.01, +0.01] | 1 |
| D - B | all copies | share gene body exact (pp) | 975 | 790 vs 784 (+6/-0) | +0.62 pp | [+0.13, +1.03] | 0.0625 |
| D - B | all copies | mean wrong bases, gene body | 973 | 5.24 vs 5.24 (closer 6/farther 0) | -0.01 | [-0.01, -0.00] | 0.0625 |
| D - B | IPD-absent | share gene body exact (pp) | 162 | 67 vs 62 (+5/-0) | +3.09 pp | [+0.85, +5.33] | 0.0625 |
| D - B | IPD-absent | mean wrong bases, gene body | 161 | 30.20 vs 30.24 (closer 5/farther 0) | -0.03 | [-0.05, -0.01] | 0.0625 |
| D - B | IPD-present | share gene body exact (pp) | 813 | 723 vs 722 (+1/-0) | +0.12 pp | [+0.00, +0.30] | 1 |
| D - B | IPD-present | mean wrong bases, gene body | 812 | 0.29 vs 0.29 (closer 1/farther 0) | -0.00 | [-0.01, +0.00] | 1 |
| E - A | all copies | share gene body exact (pp) | 975 | 849 vs 757 (+98/-6) | +9.44 pp | [+7.64, +11.19] | 5e-05 |
| E - A | all copies | mean wrong bases, gene body | 972 | 0.29 vs 5.23 (closer 101/farther 11) | -4.93 | [-13.29, -0.53] | 5e-05 |
| E - A | IPD-absent | share gene body exact (pp) | 162 | 129 vs 34 (+95/-0) | +58.64 pp | [+51.01, +65.85] | 5e-05 |
| E - A | IPD-absent | mean wrong bases, gene body | 161 | 0.80 vs 30.66 (closer 96/farther 3) | -29.86 | [-78.52, -3.29] | 5e-05 |
| E - A | IPD-present | share gene body exact (pp) | 813 | 720 vs 723 (+3/-6) | -0.37 pp | [-1.48, +0.46] | 0.672 |
| E - A | IPD-present | mean wrong bases, gene body | 811 | 0.19 vs 0.18 (closer 5/farther 8) | +0.01 | [-0.01, +0.04] | 0.359 |

## Validation against the frozen scorer

| arm | IPD-absent copies | frozen recovered_final | frozen recovered (any slot) | mine whole-seq exact (worst tied) | mine whole-seq exact (any tied) | mine gene-body exact | IPD-present copies with one 4-field truth | frozen named4 correct | mine whole-seq exact | mine gene-body exact | all IPD-present copies | mine gene-body exact (all present) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 162 | 0 | 0 | 0 | 0 | 34 | 804 | 707 | 621 | 718 | 813 | 723 |
| B | 162 | 19 | 20 | 19 | 19 | 62 | 804 | 694 | 623 | 717 | 813 | 722 |
| C | 162 | 2 | 2 | 2 | 2 | 40 | 804 | 706 | 620 | 717 | 813 | 722 |
| D | 162 | 20 | 22 | 20 | 20 | 67 | 804 | 694 | 623 | 718 | 813 | 723 |
| E | 162 | 104 | 106 | 104 | 104 | 129 | 804 | 695 | 630 | 714 | 813 | 720 |

Per-copy disagreements (frozen four-field correct vs gene-body identical on IPD-present copies with one 4-field truth; frozen recovered_final vs whole-sequence identical on IPD-absent copies): count by arm and kind

- ('A', 'gb_identical_frozen_wrong(tied)'): 7
- ('A', 'gb_identical_other_4field_name_or_ends'): 4
- ('B', 'gb_identical_frozen_wrong(tied)'): 5
- ('B', 'gb_identical_other_4field_name_or_ends'): 18
- ('C', 'gb_identical_frozen_wrong(tied)'): 7
- ('C', 'gb_identical_other_4field_name_or_ends'): 4
- ('D', 'gb_identical_frozen_wrong(tied)'): 5
- ('D', 'gb_identical_other_4field_name_or_ends'): 19
- ('E', 'gb_identical_frozen_wrong(tied)'): 4
- ('E', 'gb_identical_other_4field_name_or_ends'): 15

First disagreements (arm, donor, gene, copy, kind, frozen, mine, n_tied): ('A', 'APR-F', 'B', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('A', 'APR-F', 'C', 1, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('A', 'APR-M', 'DPA1', 2, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('A', 'APR020', 'DPA1', 1, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('A', 'APR026', 'DPA1', 1, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('A', 'APR026', 'DPA1', 2, 'gb_identical_frozen_wrong(tied)', 0, 1, 5); ('A', 'APR035', 'DPA1', 1, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('A', 'APR044', 'DPA1', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('A', 'APR044', 'DPA1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('A', 'APR048', 'DPA1', 2, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('A', 'KSA001', 'DPB1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'APR-F', 'C', 1, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('B', 'APR-M', 'DPA1', 2, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('B', 'APR004', 'DPA1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'APR007', 'DRB1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'APR008', 'B', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'APR012', 'B', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'APR020', 'DPA1', 1, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('B', 'APR026', 'DPA1', 1, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('B', 'APR026', 'DPA1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1)

Rows with unresolved call labels: 0
