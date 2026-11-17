# DRB1 patch-builder, read set T-real (test), partial-credit re-scoring

132 true DRB1 copies per arm.

**all copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 132 | 0 | 33 (25.0%) | 125 (94.7%) | 48 (36.4%) | 95 (72.0%) | 129 (97.7%) | 130 (98.5%) | 131 (99.2%) | 132 (100.0%) |
| B | 132 | 0 | 46 (34.8%) | 125 (94.7%) | 64 (48.5%) | 104 (78.8%) | 130 (98.5%) | 130 (98.5%) | 131 (99.2%) | 132 (100.0%) |
| C4 | 132 | 0 | 44 (33.3%) | 127 (96.2%) | 63 (47.7%) | 105 (79.5%) | 130 (98.5%) | 131 (99.2%) | 131 (99.2%) | 132 (100.0%) |
| D4 | 132 | 0 | 56 (42.4%) | 127 (96.2%) | 78 (59.1%) | 110 (83.3%) | 130 (98.5%) | 131 (99.2%) | 131 (99.2%) | 132 (100.0%) |

**all copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 132 | 4.6 / 2 | 0.1 / 0 | 4.4 / 2 | 0 |
| B | 132 | 3.6 / 2 | 0.1 / 0 | 3.5 / 2 | 0 |
| C4 | 132 | 3.5 / 2 | 0.1 / 0 | 3.4 / 2 | 0 |
| D4 | 132 | 2.9 / 1 | 0.1 / 0 | 2.8 / 1 | 0 |

**IPD-absent copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 85 | 0 | 2 (2.4%) | 79 (92.9%) | 15 (17.6%) | 51 (60.0%) | 82 (96.5%) | 83 (97.6%) | 84 (98.8%) | 85 (100.0%) |
| B | 85 | 0 | 16 (18.8%) | 79 (92.9%) | 32 (37.6%) | 59 (69.4%) | 83 (97.6%) | 83 (97.6%) | 84 (98.8%) | 85 (100.0%) |
| C4 | 85 | 0 | 6 (7.1%) | 81 (95.3%) | 24 (28.2%) | 61 (71.8%) | 83 (97.6%) | 84 (98.8%) | 84 (98.8%) | 85 (100.0%) |
| D4 | 85 | 0 | 17 (20.0%) | 81 (95.3%) | 38 (44.7%) | 65 (76.5%) | 83 (97.6%) | 84 (98.8%) | 84 (98.8%) | 85 (100.0%) |

**IPD-absent copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 85 | 6.5 / 4 | 0.2 / 0 | 6.3 / 4 | 0 |
| B | 85 | 5.0 / 2 | 0.2 / 0 | 4.8 / 2 | 0 |
| C4 | 85 | 5.0 / 3 | 0.1 / 0 | 4.9 / 3 | 0 |
| D4 | 85 | 4.1 / 2 | 0.1 / 0 | 4.0 / 2 | 0 |

**IPD-present copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 47 | 0 | 31 (66.0%) | 46 (97.9%) | 33 (70.2%) | 44 (93.6%) | 47 (100.0%) | 47 (100.0%) | 47 (100.0%) | 47 (100.0%) |
| B | 47 | 0 | 30 (63.8%) | 46 (97.9%) | 32 (68.1%) | 45 (95.7%) | 47 (100.0%) | 47 (100.0%) | 47 (100.0%) | 47 (100.0%) |
| C4 | 47 | 0 | 38 (80.9%) | 46 (97.9%) | 39 (83.0%) | 44 (93.6%) | 47 (100.0%) | 47 (100.0%) | 47 (100.0%) | 47 (100.0%) |
| D4 | 47 | 0 | 39 (83.0%) | 46 (97.9%) | 40 (85.1%) | 45 (95.7%) | 47 (100.0%) | 47 (100.0%) | 47 (100.0%) | 47 (100.0%) |

**IPD-present copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 47 | 1.1 / 0 | 0.0 / 0 | 1.1 / 0 | 0 |
| B | 47 | 1.1 / 0 | 0.0 / 0 | 1.0 / 0 | 0 |
| C4 | 47 | 0.7 / 0 | 0.0 / 0 | 0.7 / 0 | 0 |
| D4 | 47 | 0.6 / 0 | 0.0 / 0 | 0.6 / 0 | 0 |

**Sensitivity to the tie rule**

| arm | copies with tied alleles | gb exact (worst tied allele; primary) | gb exact (best tied allele) | mean gb wrong bases (worst) | mean gb wrong bases (best) |
|---|---|---|---|---|---|
| A | 1 | 33 (25.0%) | 33 (25.0%) | 4.56 | 4.55 |
| B | 3 | 46 (34.8%) | 46 (34.8%) | 3.60 | 3.56 |
| C4 | 0 | 44 (33.3%) | 44 (33.3%) | 3.49 | 3.49 |
| D4 | 1 | 56 (42.4%) | 56 (42.4%) | 2.88 | 2.87 |

**Paired contrasts**

| contrast | subset | metric | n paired | values (a vs b) | difference | 95% CI (family bootstrap) | sign-flip p |
|---|---|---|---|---|---|---|---|
| C4 - A | all copies | share gene body exact (pp) | 132 | 44 vs 33 (+11/-0) | +8.33 pp | [+3.79, +12.88] | 0.000977 |
| C4 - A | all copies | share coding exact (pp) | 132 | 127 vs 125 (+2/-0) | +1.52 pp | [+0.00, +3.79] | 0.5 |
| C4 - A | all copies | mean wrong bases, gene body | 132 | 3.49 vs 4.56 (closer 42/farther 1) | -1.07 | [-1.55, -0.62] | 5e-05 |
| C4 - A | all copies | mean wrong bases, coding | 132 | 0.09 vs 0.13 (closer 3/farther 0) | -0.04 | [-0.08, +0.00] | 0.25 |
| C4 - A | IPD-absent | share gene body exact (pp) | 85 | 6 vs 2 (+4/-0) | +4.71 pp | [+1.12, +9.52] | 0.125 |
| C4 - A | IPD-absent | share coding exact (pp) | 85 | 81 vs 79 (+2/-0) | +2.35 pp | [+0.00, +5.88] | 0.5 |
| C4 - A | IPD-absent | mean wrong bases, gene body | 85 | 5.01 vs 6.47 (closer 35/farther 1) | -1.46 | [-2.14, -0.82] | 5e-05 |
| C4 - A | IPD-absent | mean wrong bases, coding | 85 | 0.13 vs 0.19 (closer 3/farther 0) | -0.06 | [-0.13, +0.00] | 0.25 |
| C4 - A | IPD-present | share gene body exact (pp) | 47 | 38 vs 31 (+7/-0) | +14.89 pp | [+6.00, +26.09] | 0.0156 |
| C4 - A | IPD-present | share coding exact (pp) | 47 | 46 vs 46 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| C4 - A | IPD-present | mean wrong bases, gene body | 47 | 0.74 vs 1.11 (closer 7/farther 0) | -0.36 | [-0.67, -0.12] | 0.0156 |
| C4 - A | IPD-present | mean wrong bases, coding | 47 | 0.02 vs 0.02 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| D4 - B | all copies | share gene body exact (pp) | 132 | 56 vs 46 (+10/-0) | +7.58 pp | [+3.79, +12.12] | 0.00195 |
| D4 - B | all copies | share coding exact (pp) | 132 | 127 vs 125 (+2/-0) | +1.52 pp | [+0.00, +3.79] | 0.5 |
| D4 - B | all copies | mean wrong bases, gene body | 132 | 2.88 vs 3.60 (closer 33/farther 1) | -0.72 | [-1.14, -0.36] | 0.0001 |
| D4 - B | all copies | mean wrong bases, coding | 132 | 0.09 vs 0.13 (closer 3/farther 0) | -0.04 | [-0.08, +0.00] | 0.25 |
| D4 - B | IPD-absent | share gene body exact (pp) | 85 | 17 vs 16 (+1/-0) | +1.18 pp | [+0.00, +3.75] | 1 |
| D4 - B | IPD-absent | share coding exact (pp) | 85 | 81 vs 79 (+2/-0) | +2.35 pp | [+0.00, +5.88] | 0.5 |
| D4 - B | IPD-absent | mean wrong bases, gene body | 85 | 4.13 vs 5.00 (closer 24/farther 1) | -0.87 | [-1.49, -0.33] | 0.00315 |
| D4 - B | IPD-absent | mean wrong bases, coding | 85 | 0.13 vs 0.19 (closer 3/farther 0) | -0.06 | [-0.13, +0.00] | 0.25 |
| D4 - B | IPD-present | share gene body exact (pp) | 47 | 39 vs 30 (+9/-0) | +19.15 pp | [+8.33, +32.00] | 0.00391 |
| D4 - B | IPD-present | share coding exact (pp) | 47 | 46 vs 46 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| D4 - B | IPD-present | mean wrong bases, gene body | 47 | 0.62 vs 1.06 (closer 9/farther 0) | -0.45 | [-0.77, -0.19] | 0.00391 |
| D4 - B | IPD-present | mean wrong bases, coding | 47 | 0.02 vs 0.02 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| B - A | all copies | share gene body exact (pp) | 132 | 46 vs 33 (+14/-1) | +9.85 pp | [+4.55, +15.91] | 0.00293 |
| B - A | all copies | mean wrong bases, gene body | 132 | 3.60 vs 4.56 (closer 32/farther 3) | -0.96 | [-1.44, -0.55] | 0.0001 |
| B - A | IPD-absent | share gene body exact (pp) | 85 | 16 vs 2 (+14/-0) | +16.47 pp | [+8.33, +25.68] | 0.000488 |
| B - A | IPD-absent | mean wrong bases, gene body | 85 | 5.00 vs 6.47 (closer 31/farther 2) | -1.47 | [-2.10, -0.86] | 5e-05 |
| B - A | IPD-present | share gene body exact (pp) | 47 | 30 vs 31 (+0/-1) | -2.13 pp | [-6.82, +0.00] | 1 |
| B - A | IPD-present | mean wrong bases, gene body | 47 | 1.06 vs 1.11 (closer 1/farther 1) | -0.04 | [-0.25, +0.12] | 1 |

**C4 final novel (patched) calls vs the A call for the same copy (gene-body distance to truth)**

| read set | copy type | copies with novel final call | comparable | not comparable (arm A call missing) | closer than A | equal | farther |
|---|---|---|---|---|---|---|---|
| T-real | IPD-absent | 36 | 36 | 0 | 35 | 0 | 1 |
| T-real | IPD-present | 0 | 0 | 0 | 0 | 0 | 0 |

Change in gene-body wrong bases (C4 minus A); negative = patched call is closer:

| read set | copy type | n | min | median | max | <=-20 | -19..-5 | -4..-1 | 0 | +1..+4 | +5..+19 | >=+20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T-real | IPD-absent | 36 | -14 | -2.0 | 12 | 0 | 10 | 25 | 0 | 0 | 1 | 0 |

**D4 final novel (patched) calls vs the B call for the same copy (gene-body distance to truth)**

| read set | copy type | copies with novel final call | comparable | not comparable (arm B call missing) | closer than B | equal | farther |
|---|---|---|---|---|---|---|---|
| T-real | IPD-absent | 53 | 53 | 0 | 24 | 28 | 1 |
| T-real | IPD-present | 2 | 2 | 0 | 0 | 2 | 0 |

Change in gene-body wrong bases (D4 minus B); negative = patched call is closer:

| read set | copy type | n | min | median | max | <=-20 | -19..-5 | -4..-1 | 0 | +1..+4 | +5..+19 | >=+20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T-real | IPD-absent | 53 | -14 | 0 | 12 | 0 | 5 | 19 | 28 | 0 | 1 | 0 |
| T-real | IPD-present | 2 | 0 | 0.0 | 0 | 0 | 0 | 0 | 2 | 0 | 0 | 0 |

## Validation against the frozen scorer (DRB1 only)

| arm | IPD-absent copies | frozen recovered_final | frozen recovered (any slot) | mine whole-seq exact (worst tied) | mine whole-seq exact (any tied) | mine gene-body exact | IPD-present copies with one 4-field truth | frozen named4 correct | mine whole-seq exact | mine gene-body exact | all IPD-present copies | mine gene-body exact (all present) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 85 | 0 | 0 | 0 | 0 | 2 | 47 | 31 | 28 | 31 | 47 | 31 |
| B | 85 | 14 | 15 | 14 | 14 | 16 | 47 | 29 | 28 | 30 | 47 | 30 |
| C4 | 85 | 2 | 3 | 2 | 2 | 6 | 47 | 38 | 35 | 38 | 47 | 38 |
| D4 | 85 | 15 | 17 | 15 | 15 | 17 | 47 | 38 | 37 | 39 | 47 | 39 |

Per-copy disagreements (frozen four-field correct vs gene-body identical on IPD-present copies with one 4-field truth; frozen recovered_final vs whole-sequence identical on IPD-absent copies): count by arm and kind

- ('B', 'gb_identical_other_4field_name_or_ends'): 1
- ('D4', 'gb_identical_other_4field_name_or_ends'): 1

First disagreements (arm, donor, gene, copy, kind, frozen, mine, n_tied): ('B', 'HG03816', 'DRB1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1); ('D4', 'HG03816', 'DRB1', 2, 'gb_identical_other_4field_name_or_ends', 0, 1, 1)

Gene-body recovery on IPD-absent DRB1 copies (frozen "gene body exact" column vs mine):

| arm | absent copies | frozen gene body exact (candidate in ANY slot) | same, restricted to FINAL calls (any slot) | mine gene-body exact (assigned slot, worst tied) | mine (any tied) |
|---|---|---|---|---|---|
| A | 85 | 0 | 0 | 2 | 2 |
| B | 85 | 15 | 14 | 16 | 16 |
| C4 | 85 | 5 | 4 | 6 | 6 |
| D4 | 85 | 17 | 15 | 17 | 17 |
