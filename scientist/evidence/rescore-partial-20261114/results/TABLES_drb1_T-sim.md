# DRB1 patch-builder, read set T-sim (test), partial-credit re-scoring

122 true DRB1 copies per arm.

**all copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 122 | 0 | 36 (29.5%) | 111 (91.0%) | 46 (37.7%) | 92 (75.4%) | 117 (95.9%) | 118 (96.7%) | 120 (98.4%) | 121 (99.2%) |
| B | 122 | 0 | 52 (42.6%) | 110 (90.2%) | 61 (50.0%) | 98 (80.3%) | 117 (95.9%) | 118 (96.7%) | 120 (98.4%) | 121 (99.2%) |
| C4 | 122 | 0 | 55 (45.1%) | 113 (92.6%) | 66 (54.1%) | 100 (82.0%) | 117 (95.9%) | 118 (96.7%) | 120 (98.4%) | 121 (99.2%) |
| D4 | 122 | 2 | 64 (52.5%) | 111 (91.0%) | 73 (59.8%) | 102 (83.6%) | 115 (94.3%) | 116 (95.1%) | 118 (96.7%) | 119 (97.5%) |

**all copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 122 | 31.1 / 2 | 0.3 / 0 | 30.8 / 2 | 0 |
| B | 122 | 30.6 / 1.5 | 0.3 / 0 | 30.2 / 1 | 0 |
| C4 | 122 | 30.3 / 1 | 0.3 / 0 | 30.0 / 1 | 0 |
| D4 | 120 | 30.3 / 0 | 0.3 / 0 | 30.0 / 0 | 0 |

**IPD-absent copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 76 | 0 | 0 (0.0%) | 66 (86.8%) | 10 (13.2%) | 49 (64.5%) | 71 (93.4%) | 72 (94.7%) | 74 (97.4%) | 75 (98.7%) |
| B | 76 | 0 | 16 (21.1%) | 65 (85.5%) | 24 (31.6%) | 53 (69.7%) | 71 (93.4%) | 72 (94.7%) | 74 (97.4%) | 75 (98.7%) |
| C4 | 76 | 0 | 14 (18.4%) | 67 (88.2%) | 24 (31.6%) | 57 (75.0%) | 71 (93.4%) | 72 (94.7%) | 74 (97.4%) | 75 (98.7%) |
| D4 | 76 | 1 | 24 (31.6%) | 66 (86.8%) | 32 (42.1%) | 58 (76.3%) | 70 (92.1%) | 71 (93.4%) | 73 (96.1%) | 74 (97.4%) |

**IPD-absent copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 76 | 49.4 / 4 | 0.5 / 0 | 48.9 / 4 | 0 |
| B | 76 | 48.7 / 2 | 0.5 / 0 | 48.1 / 2 | 0 |
| C4 | 76 | 48.2 / 3 | 0.5 / 0 | 47.7 / 2 | 0 |
| D4 | 75 | 48.3 / 2 | 0.5 / 0 | 47.8 / 2 | 0 |

**IPD-present copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 46 | 0 | 36 (78.3%) | 45 (97.8%) | 36 (78.3%) | 43 (93.5%) | 46 (100.0%) | 46 (100.0%) | 46 (100.0%) | 46 (100.0%) |
| B | 46 | 0 | 36 (78.3%) | 45 (97.8%) | 37 (80.4%) | 45 (97.8%) | 46 (100.0%) | 46 (100.0%) | 46 (100.0%) | 46 (100.0%) |
| C4 | 46 | 0 | 41 (89.1%) | 46 (100.0%) | 42 (91.3%) | 43 (93.5%) | 46 (100.0%) | 46 (100.0%) | 46 (100.0%) | 46 (100.0%) |
| D4 | 46 | 1 | 40 (87.0%) | 45 (97.8%) | 41 (89.1%) | 44 (95.7%) | 45 (97.8%) | 45 (97.8%) | 45 (97.8%) | 45 (97.8%) |

**IPD-present copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 46 | 0.9 / 0 | 0.0 / 0 | 0.9 / 0 | 0 |
| B | 46 | 0.7 / 0 | 0.0 / 0 | 0.7 / 0 | 0 |
| C4 | 46 | 0.6 / 0 | 0.0 / 0 | 0.6 / 0 | 0 |
| D4 | 45 | 0.4 / 0 | 0.0 / 0 | 0.4 / 0 | 0 |

**Sensitivity to the tie rule**

| arm | copies with tied alleles | gb exact (worst tied allele; primary) | gb exact (best tied allele) | mean gb wrong bases (worst) | mean gb wrong bases (best) |
|---|---|---|---|---|---|
| A | 1 | 36 (29.5%) | 36 (29.5%) | 31.11 | 31.10 |
| B | 4 | 52 (42.6%) | 52 (42.6%) | 30.58 | 30.34 |
| C4 | 0 | 55 (45.1%) | 55 (45.1%) | 30.26 | 30.26 |
| D4 | 3 | 64 (52.5%) | 64 (52.5%) | 30.33 | 30.17 |

**Paired contrasts**

| contrast | subset | metric | n paired | values (a vs b) | difference | 95% CI (family bootstrap) | sign-flip p |
|---|---|---|---|---|---|---|---|
| C4 - A | all copies | share gene body exact (pp) | 122 | 55 vs 36 (+19/-0) | +15.57 pp | [+9.33, +23.15] | 1.53e-05 |
| C4 - A | all copies | share coding exact (pp) | 122 | 113 vs 111 (+2/-0) | +1.64 pp | [+0.00, +4.46] | 0.5 |
| C4 - A | all copies | mean wrong bases, gene body | 122 | 30.26 vs 31.11 (closer 35/farther 0) | -0.85 | [-1.13, -0.54] | 5e-05 |
| C4 - A | all copies | mean wrong bases, coding | 122 | 0.31 vs 0.33 (closer 2/farther 0) | -0.02 | [-0.04, +0.00] | 0.5 |
| C4 - A | IPD-absent | share gene body exact (pp) | 76 | 14 vs 0 (+14/-0) | +18.42 pp | [+10.26, +28.38] | 0.000488 |
| C4 - A | IPD-absent | share coding exact (pp) | 76 | 67 vs 66 (+1/-0) | +1.32 pp | [+0.00, +4.84] | 1 |
| C4 - A | IPD-absent | mean wrong bases, gene body | 76 | 48.24 vs 49.39 (closer 29/farther 0) | -1.16 | [-1.59, -0.67] | 5e-05 |
| C4 - A | IPD-absent | mean wrong bases, coding | 76 | 0.50 vs 0.51 (closer 1/farther 0) | -0.01 | [-0.05, +0.00] | 1 |
| C4 - A | IPD-present | share gene body exact (pp) | 46 | 41 vs 36 (+5/-0) | +10.87 pp | [+2.33, +20.00] | 0.0625 |
| C4 - A | IPD-present | share coding exact (pp) | 46 | 46 vs 45 (+1/-0) | +2.17 pp | [+0.00, +6.98] | 1 |
| C4 - A | IPD-present | mean wrong bases, gene body | 46 | 0.57 vs 0.91 (closer 6/farther 0) | -0.35 | [-0.66, -0.08] | 0.0625 |
| C4 - A | IPD-present | mean wrong bases, coding | 46 | 0.00 vs 0.02 (closer 1/farther 0) | -0.02 | [-0.07, +0.00] | 1 |
| D4 - B | all copies | share gene body exact (pp) | 122 | 64 vs 52 (+14/-2) | +9.84 pp | [+4.17, +16.36] | 0.00653 |
| D4 - B | all copies | share coding exact (pp) | 122 | 111 vs 110 (+3/-2) | +0.82 pp | [-4.00, +5.00] | 1 |
| D4 - B | all copies | mean wrong bases, gene body | 120 | 30.33 vs 31.03 (closer 26/farther 2) | -0.70 | [-1.08, -0.39] | 5e-05 |
| D4 - B | all copies | mean wrong bases, coding | 120 | 0.32 vs 0.34 (closer 3/farther 0) | -0.03 | [-0.06, +0.00] | 0.25 |
| D4 - B | IPD-absent | share gene body exact (pp) | 76 | 24 vs 16 (+9/-1) | +10.53 pp | [+3.08, +18.31] | 0.0352 |
| D4 - B | IPD-absent | share coding exact (pp) | 76 | 66 vs 65 (+2/-1) | +1.32 pp | [-3.23, +6.35] | 1 |
| D4 - B | IPD-absent | mean wrong bases, gene body | 75 | 48.29 vs 49.23 (closer 20/farther 2) | -0.93 | [-1.55, -0.46] | 0.000793 |
| D4 - B | IPD-absent | mean wrong bases, coding | 75 | 0.51 vs 0.53 (closer 2/farther 0) | -0.03 | [-0.07, +0.00] | 0.5 |
| D4 - B | IPD-present | share gene body exact (pp) | 46 | 40 vs 36 (+5/-1) | +8.70 pp | [+0.00, +18.60] | 0.219 |
| D4 - B | IPD-present | share coding exact (pp) | 46 | 45 vs 45 (+1/-1) | +0.00 pp | [-6.82, +6.38] | 1 |
| D4 - B | IPD-present | mean wrong bases, gene body | 45 | 0.40 vs 0.71 (closer 6/farther 0) | -0.31 | [-0.62, -0.06] | 0.0625 |
| D4 - B | IPD-present | mean wrong bases, coding | 45 | 0.00 vs 0.02 (closer 1/farther 0) | -0.02 | [-0.07, +0.00] | 1 |
| B - A | all copies | share gene body exact (pp) | 122 | 52 vs 36 (+17/-1) | +13.11 pp | [+7.00, +18.10] | 0.000244 |
| B - A | all copies | mean wrong bases, gene body | 122 | 30.58 vs 31.11 (closer 33/farther 9) | -0.53 | [-0.83, -0.16] | 0.0124 |
| B - A | IPD-absent | share gene body exact (pp) | 76 | 16 vs 0 (+16/-0) | +21.05 pp | [+12.68, +28.75] | 0.000244 |
| B - A | IPD-absent | mean wrong bases, gene body | 76 | 48.67 vs 49.39 (closer 29/farther 8) | -0.72 | [-1.15, -0.11] | 0.0372 |
| B - A | IPD-present | share gene body exact (pp) | 46 | 36 vs 36 (+1/-1) | +0.00 pp | [+0.00, +0.00] | 1 |
| B - A | IPD-present | mean wrong bases, gene body | 46 | 0.70 vs 0.91 (closer 4/farther 1) | -0.22 | [-0.48, +0.00] | 0.25 |

**C4 final novel (patched) calls vs the A call for the same copy (gene-body distance to truth)**

| read set | copy type | copies with novel final call | comparable | not comparable (arm A call missing) | closer than A | equal | farther |
|---|---|---|---|---|---|---|---|
| T-sim | IPD-absent | 28 | 28 | 0 | 28 | 0 | 0 |
| T-sim | IPD-present | 1 | 1 | 0 | 1 | 0 | 0 |

Change in gene-body wrong bases (C4 minus A); negative = patched call is closer:

| read set | copy type | n | min | median | max | <=-20 | -19..-5 | -4..-1 | 0 | +1..+4 | +5..+19 | >=+20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T-sim | IPD-absent | 28 | -11 | -2.0 | -1 | 0 | 6 | 22 | 0 | 0 | 0 | 0 |
| T-sim | IPD-present | 1 | -1 | -1 | -1 | 0 | 0 | 1 | 0 | 0 | 0 | 0 |

**D4 final novel (patched) calls vs the B call for the same copy (gene-body distance to truth)**

| read set | copy type | copies with novel final call | comparable | not comparable (arm B call missing) | closer than B | equal | farther |
|---|---|---|---|---|---|---|---|
| T-sim | IPD-absent | 45 | 45 | 0 | 17 | 27 | 1 |
| T-sim | IPD-present | 4 | 4 | 0 | 1 | 3 | 0 |

Change in gene-body wrong bases (D4 minus B); negative = patched call is closer:

| read set | copy type | n | min | median | max | <=-20 | -19..-5 | -4..-1 | 0 | +1..+4 | +5..+19 | >=+20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T-sim | IPD-absent | 45 | -14 | 0 | 1 | 0 | 4 | 13 | 27 | 1 | 0 | 0 |
| T-sim | IPD-present | 4 | -1 | 0.0 | 0 | 0 | 0 | 1 | 3 | 0 | 0 | 0 |

## Validation against the frozen scorer (DRB1 only)

| arm | IPD-absent copies | frozen recovered_final | frozen recovered (any slot) | mine whole-seq exact (worst tied) | mine whole-seq exact (any tied) | mine gene-body exact | IPD-present copies with one 4-field truth | frozen named4 correct | mine whole-seq exact | mine gene-body exact | all IPD-present copies | mine gene-body exact (all present) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 76 | 0 | 0 | 0 | 0 | 0 | 45 | 36 | 29 | 36 | 46 | 36 |
| B | 76 | 16 | 17 | 16 | 16 | 16 | 45 | 35 | 30 | 36 | 46 | 36 |
| C4 | 76 | 10 | 10 | 10 | 10 | 14 | 45 | 40 | 33 | 40 | 46 | 41 |
| D4 | 76 | 22 | 23 | 21 | 21 | 24 | 45 | 38 | 34 | 39 | 46 | 40 |

Per-copy disagreements (frozen four-field correct vs gene-body identical on IPD-present copies with one 4-field truth; frozen recovered_final vs whole-sequence identical on IPD-absent copies): count by arm and kind

- ('B', 'gb_identical_other_4field_name_or_ends'): 1
- ('D4', 'absent'): 1
- ('D4', 'gb_identical_other_4field_name_or_ends'): 1

First disagreements (arm, donor, gene, copy, kind, frozen, mine, n_tied): ('B', 'APR007', 'DRB1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('D4', 'APR038', 'DRB1', 2, 'absent', 1, 0, 1); ('D4', 'APR007', 'DRB1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1)

Gene-body recovery on IPD-absent DRB1 copies (frozen "gene body exact" column vs mine):

| arm | absent copies | frozen gene body exact (candidate in ANY slot) | same, restricted to FINAL calls (any slot) | mine gene-body exact (assigned slot, worst tied) | mine (any tied) |
|---|---|---|---|---|---|
| A | 76 | 0 | 0 | 0 | 0 |
| B | 76 | 17 | 17 | 16 | 16 |
| C4 | 76 | 16 | 14 | 14 | 14 |
| D4 | 76 | 26 | 26 | 24 | 24 |
