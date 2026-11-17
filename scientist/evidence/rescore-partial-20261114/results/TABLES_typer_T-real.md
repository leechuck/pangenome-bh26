# Mosaic-typer test, read set T-real (all 8 genes), partial-credit re-scoring

1056 true copies per arm.

**all copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 1056 | 0 | 902 (85.4%) | 1038 (98.3%) | 954 (90.3%) | 1015 (96.1%) | 1053 (99.7%) | 1054 (99.8%) | 1055 (99.9%) | 1056 (100.0%) |
| B | 1056 | 0 | 923 (87.4%) | 1039 (98.4%) | 972 (92.0%) | 1025 (97.1%) | 1054 (99.8%) | 1054 (99.8%) | 1055 (99.9%) | 1056 (100.0%) |
| C | 1056 | 0 | 905 (85.7%) | 1040 (98.5%) | 955 (90.4%) | 1015 (96.1%) | 1053 (99.7%) | 1054 (99.8%) | 1055 (99.9%) | 1056 (100.0%) |
| D | 1056 | 0 | 925 (87.6%) | 1040 (98.5%) | 973 (92.1%) | 1025 (97.1%) | 1054 (99.8%) | 1054 (99.8%) | 1055 (99.9%) | 1056 (100.0%) |
| E | 1056 | 0 | 968 (91.7%) | 1045 (99.0%) | 1004 (95.1%) | 1028 (97.3%) | 1053 (99.7%) | 1054 (99.8%) | 1055 (99.9%) | 1056 (100.0%) |

**all copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 1056 | 0.7 / 0 | 0.0 / 0 | 0.7 / 0 | 0 |
| B | 1056 | 0.5 / 0 | 0.0 / 0 | 0.5 / 0 | 0 |
| C | 1056 | 0.7 / 0 | 0.0 / 0 | 0.7 / 0 | 0 |
| D | 1056 | 0.5 / 0 | 0.0 / 0 | 0.5 / 0 | 0 |
| E | 1056 | 0.5 / 0 | 0.0 / 0 | 0.4 / 0 | 0 |

**IPD-absent copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 134 | 0 | 24 (17.9%) | 124 (92.5%) | 51 (38.1%) | 99 (73.9%) | 131 (97.8%) | 132 (98.5%) | 133 (99.3%) | 134 (100.0%) |
| B | 134 | 0 | 47 (35.1%) | 125 (93.3%) | 71 (53.0%) | 107 (79.9%) | 132 (98.5%) | 132 (98.5%) | 133 (99.3%) | 134 (100.0%) |
| C | 134 | 0 | 27 (20.1%) | 126 (94.0%) | 52 (38.8%) | 99 (73.9%) | 131 (97.8%) | 132 (98.5%) | 133 (99.3%) | 134 (100.0%) |
| D | 134 | 0 | 49 (36.6%) | 126 (94.0%) | 72 (53.7%) | 107 (79.9%) | 132 (98.5%) | 132 (98.5%) | 133 (99.3%) | 134 (100.0%) |
| E | 134 | 0 | 90 (67.2%) | 131 (97.8%) | 101 (75.4%) | 112 (83.6%) | 131 (97.8%) | 132 (98.5%) | 133 (99.3%) | 134 (100.0%) |

**IPD-absent copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 134 | 4.6 / 2 | 0.1 / 0 | 4.4 / 2 | 0 |
| B | 134 | 3.5 / 1 | 0.1 / 0 | 3.4 / 1 | 0 |
| C | 134 | 4.5 / 2 | 0.1 / 0 | 4.4 / 2 | 0 |
| D | 134 | 3.5 / 1 | 0.1 / 0 | 3.3 / 1 | 0 |
| E | 134 | 2.6 / 0 | 0.1 / 0 | 2.6 / 0 | 0 |

**IPD-present copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 922 | 0 | 878 (95.2%) | 914 (99.1%) | 903 (97.9%) | 916 (99.3%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) |
| B | 922 | 0 | 876 (95.0%) | 914 (99.1%) | 901 (97.7%) | 918 (99.6%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) |
| C | 922 | 0 | 878 (95.2%) | 914 (99.1%) | 903 (97.9%) | 916 (99.3%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) |
| D | 922 | 0 | 876 (95.0%) | 914 (99.1%) | 901 (97.7%) | 918 (99.6%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) |
| E | 922 | 0 | 878 (95.2%) | 914 (99.1%) | 903 (97.9%) | 916 (99.3%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) | 922 (100.0%) |

**IPD-present copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 922 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |
| B | 922 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |
| C | 922 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |
| D | 922 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |
| E | 922 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |

**By gene**

| gene | A: copies; gb exact; mean gb / mean cod | B: copies; gb exact; mean gb / mean cod | C: copies; gb exact; mean gb / mean cod | D: copies; gb exact; mean gb / mean cod | E: copies; gb exact; mean gb / mean cod |
|---|---|---|---|---|---|
| A | 132; 129; 0.0 / 0.0 | 132; 129; 0.0 / 0.0 | 132; 129; 0.0 / 0.0 | 132; 129; 0.0 / 0.0 | 132; 130; 0.0 / 0.0 |
| B | 132; 125; 0.1 / 0.0 | 132; 125; 0.1 / 0.0 | 132; 125; 0.1 / 0.0 | 132; 125; 0.1 / 0.0 | 132; 126; 0.1 / 0.0 |
| C | 132; 131; 0.0 / 0.0 | 132; 131; 0.0 / 0.0 | 132; 131; 0.0 / 0.0 | 132; 131; 0.0 / 0.0 | 132; 132; 0.0 / 0.0 |
| DPA1 | 132; 123; 0.2 / 0.0 | 132; 123; 0.2 / 0.0 | 132; 124; 0.2 / 0.0 | 132; 124; 0.2 / 0.0 | 132; 126; 0.2 / 0.0 |
| DPB1 | 132; 119; 0.3 / 0.0 | 132; 120; 0.2 / 0.0 | 132; 119; 0.3 / 0.0 | 132; 120; 0.2 / 0.0 | 132; 129; 0.0 / 0.0 |
| DQA1 | 132; 112; 0.4 / 0.0 | 132; 119; 0.2 / 0.0 | 132; 114; 0.3 / 0.0 | 132; 120; 0.2 / 0.0 | 132; 120; 0.3 / 0.0 |
| DQB1 | 132; 130; 0.0 / 0.0 | 132; 130; 0.0 / 0.0 | 132; 130; 0.0 / 0.0 | 132; 130; 0.0 / 0.0 | 132; 130; 0.0 / 0.0 |
| DRB1 | 132; 33; 4.6 / 0.1 | 132; 46; 3.6 / 0.1 | 132; 33; 4.6 / 0.1 | 132; 46; 3.6 / 0.1 | 132; 75; 3.0 / 0.1 |

**Sensitivity to the tie rule**

| arm | copies with tied alleles | gb exact (worst tied allele; primary) | gb exact (best tied allele) | mean gb wrong bases (worst) | mean gb wrong bases (best) |
|---|---|---|---|---|---|
| A | 7 | 902 (85.4%) | 902 (85.4%) | 0.69 | 0.69 |
| B | 9 | 923 (87.4%) | 923 (87.4%) | 0.55 | 0.54 |
| C | 7 | 905 (85.7%) | 905 (85.7%) | 0.69 | 0.68 |
| D | 9 | 925 (87.6%) | 925 (87.6%) | 0.54 | 0.53 |
| E | 7 | 968 (91.7%) | 968 (91.7%) | 0.45 | 0.45 |

**Paired contrasts**

| contrast | subset | metric | n paired | values (a vs b) | difference | 95% CI (family bootstrap) | sign-flip p |
|---|---|---|---|---|---|---|---|
| B - A | all copies | share gene body exact (pp) | 1056 | 923 vs 902 (+24/-3) | +1.99 pp | [+1.04, +3.03] | 0.0003 |
| B - A | all copies | share coding exact (pp) | 1056 | 1039 vs 1038 (+3/-2) | +0.09 pp | [-0.19, +0.38] | 1 |
| B - A | all copies | mean wrong bases, gene body | 1056 | 0.55 vs 0.69 (closer 43/farther 5) | -0.15 | [-0.21, -0.09] | 5e-05 |
| B - A | all copies | mean wrong bases, coding | 1056 | 0.03 vs 0.03 (closer 3/farther 2) | -0.00 | [-0.00, +0.00] | 1 |
| B - A | IPD-absent | share gene body exact (pp) | 134 | 47 vs 24 (+23/-0) | +17.16 pp | [+10.07, +24.62] | 1.53e-05 |
| B - A | IPD-absent | share coding exact (pp) | 134 | 125 vs 124 (+2/-1) | +0.75 pp | [+0.00, +2.33] | 1 |
| B - A | IPD-absent | mean wrong bases, gene body | 134 | 3.49 vs 4.56 (closer 40/farther 2) | -1.07 | [-1.52, -0.63] | 5e-05 |
| B - A | IPD-absent | mean wrong bases, coding | 134 | 0.14 vs 0.15 (closer 2/farther 1) | -0.01 | [-0.02, +0.00] | 1 |
| B - A | IPD-present | share gene body exact (pp) | 922 | 876 vs 878 (+1/-3) | -0.22 pp | [-0.54, +0.00] | 0.5 |
| B - A | IPD-present | share coding exact (pp) | 922 | 914 vs 914 (+1/-1) | +0.00 pp | [-0.33, +0.32] | 1 |
| B - A | IPD-present | mean wrong bases, gene body | 922 | 0.12 vs 0.13 (closer 3/farther 3) | -0.01 | [-0.04, +0.01] | 0.562 |
| B - A | IPD-present | mean wrong bases, coding | 922 | 0.01 vs 0.01 (closer 1/farther 1) | +0.00 | [-0.00, +0.00] | 1 |
| D - B | all copies | share gene body exact (pp) | 1056 | 925 vs 923 (+2/-0) | +0.19 pp | [+0.00, +0.47] | 0.5 |
| D - B | all copies | mean wrong bases, gene body | 1056 | 0.54 vs 0.55 (closer 2/farther 0) | -0.00 | [-0.01, +0.00] | 0.5 |
| D - B | IPD-absent | share gene body exact (pp) | 134 | 49 vs 47 (+2/-0) | +1.49 pp | [+0.00, +3.82] | 0.5 |
| D - B | IPD-absent | mean wrong bases, gene body | 134 | 3.46 vs 3.49 (closer 2/farther 0) | -0.04 | [-0.11, +0.00] | 0.5 |
| D - B | IPD-present | share gene body exact (pp) | 922 | 876 vs 876 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| D - B | IPD-present | mean wrong bases, gene body | 922 | 0.12 vs 0.12 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| E - A | all copies | share gene body exact (pp) | 1056 | 968 vs 902 (+66/-0) | +6.25 pp | [+4.83, +7.77] | 5e-05 |
| E - A | all copies | mean wrong bases, gene body | 1056 | 0.45 vs 0.69 (closer 69/farther 1) | -0.24 | [-0.32, -0.17] | 5e-05 |
| E - A | IPD-absent | share gene body exact (pp) | 134 | 90 vs 24 (+66/-0) | +49.25 pp | [+40.00, +58.04] | 5e-05 |
| E - A | IPD-absent | mean wrong bases, gene body | 134 | 2.65 vs 4.56 (closer 69/farther 1) | -1.91 | [-2.52, -1.38] | 5e-05 |
| E - A | IPD-present | share gene body exact (pp) | 922 | 878 vs 878 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| E - A | IPD-present | mean wrong bases, gene body | 922 | 0.13 vs 0.13 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |

## Validation against the frozen scorer

| arm | IPD-absent copies | frozen recovered_final | frozen recovered (any slot) | mine whole-seq exact (worst tied) | mine whole-seq exact (any tied) | mine gene-body exact | IPD-present copies with one 4-field truth | frozen named4 correct | mine whole-seq exact | mine gene-body exact | all IPD-present copies | mine gene-body exact (all present) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 134 | 0 | 0 | 0 | 0 | 24 | 914 | 861 | 818 | 871 | 922 | 878 |
| B | 134 | 28 | 29 | 28 | 28 | 47 | 914 | 844 | 812 | 869 | 922 | 876 |
| C | 134 | 2 | 3 | 2 | 2 | 27 | 914 | 861 | 818 | 871 | 922 | 878 |
| D | 134 | 29 | 30 | 29 | 29 | 49 | 914 | 844 | 812 | 869 | 922 | 876 |
| E | 134 | 85 | 89 | 84 | 84 | 90 | 914 | 857 | 822 | 871 | 922 | 878 |

Per-copy disagreements (frozen four-field correct vs gene-body identical on IPD-present copies with one 4-field truth; frozen recovered_final vs whole-sequence identical on IPD-absent copies): count by arm and kind

- ('A', 'frozen_right_gb_not'): 1
- ('A', 'gb_identical_frozen_wrong(tied)'): 5
- ('A', 'gb_identical_other_4field_name_or_ends'): 6
- ('B', 'frozen_right_gb_not'): 1
- ('B', 'gb_identical_frozen_wrong(tied)'): 5
- ('B', 'gb_identical_other_4field_name_or_ends'): 21
- ('C', 'frozen_right_gb_not'): 1
- ('C', 'gb_identical_frozen_wrong(tied)'): 5
- ('C', 'gb_identical_other_4field_name_or_ends'): 6
- ('D', 'frozen_right_gb_not'): 1
- ('D', 'gb_identical_frozen_wrong(tied)'): 5
- ('D', 'gb_identical_other_4field_name_or_ends'): 21
- ('E', 'absent'): 1
- ('E', 'frozen_right_gb_not'): 1
- ('E', 'gb_identical_frozen_wrong(tied)'): 5
- ('E', 'gb_identical_other_4field_name_or_ends'): 10

First disagreements (arm, donor, gene, copy, kind, frozen, mine, n_tied): ('A', 'HG02132', 'B', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('A', 'HG02135', 'B', 2, 'frozen_right_gb_not', 1, 0, 1); ('A', 'HG02135', 'C', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('A', 'HG03239', 'B', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('A', 'HG03688', 'B', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('A', 'NA18608', 'DQB1', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('A', 'NA18940', 'C', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('A', 'NA18959', 'DPA1', 1, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('A', 'NA18959', 'DPA1', 2, 'gb_identical_frozen_wrong(tied)', 0, 1, 4); ('A', 'NA18960', 'DPA1', 2, 'gb_identical_frozen_wrong(tied)', 0, 1, 5); ('A', 'NA18982', 'DPA1', 1, 'gb_identical_frozen_wrong(tied)', 0, 1, 2); ('A', 'NA18982', 'DPA1', 2, 'gb_identical_frozen_wrong(tied)', 0, 1, 4); ('B', 'HG00408', 'B', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'HG00597', 'B', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'HG02027', 'B', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'HG02071', 'DQB1', 1, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'HG02083', 'DPA1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'HG02132', 'B', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'HG02132', 'DPA1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('B', 'HG02135', 'B', 2, 'frozen_right_gb_not', 1, 0, 1)

Rows with unresolved call labels: 0
