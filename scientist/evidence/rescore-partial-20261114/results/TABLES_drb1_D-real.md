# DRB1 patch-builder, read set D-real (DEVELOPMENT), partial-credit re-scoring

40 true DRB1 copies per arm.

**all copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 40 | 0 | 8 (20.0%) | 40 (100.0%) | 17 (42.5%) | 33 (82.5%) | 40 (100.0%) | 40 (100.0%) | 40 (100.0%) | 40 (100.0%) |
| B | 40 | 0 | 12 (30.0%) | 40 (100.0%) | 21 (52.5%) | 35 (87.5%) | 40 (100.0%) | 40 (100.0%) | 40 (100.0%) | 40 (100.0%) |
| C4 | 40 | 0 | 15 (37.5%) | 40 (100.0%) | 22 (55.0%) | 34 (85.0%) | 40 (100.0%) | 40 (100.0%) | 40 (100.0%) | 40 (100.0%) |
| D4 | 40 | 0 | 19 (47.5%) | 40 (100.0%) | 26 (65.0%) | 36 (90.0%) | 40 (100.0%) | 40 (100.0%) | 40 (100.0%) | 40 (100.0%) |

**all copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 40 | 3.2 / 2 | 0.0 / 0 | 3.2 / 2 | 0 |
| B | 40 | 2.8 / 1 | 0.0 / 0 | 2.8 / 1 | 0 |
| C4 | 40 | 2.5 / 1 | 0.0 / 0 | 2.5 / 1 | 0 |
| D4 | 40 | 2.1 / 1 | 0.0 / 0 | 2.1 / 1 | 0 |

**IPD-absent copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 25 | 0 | 0 (0.0%) | 25 (100.0%) | 9 (36.0%) | 20 (80.0%) | 25 (100.0%) | 25 (100.0%) | 25 (100.0%) | 25 (100.0%) |
| B | 25 | 0 | 4 (16.0%) | 25 (100.0%) | 13 (52.0%) | 22 (88.0%) | 25 (100.0%) | 25 (100.0%) | 25 (100.0%) | 25 (100.0%) |
| C4 | 25 | 0 | 3 (12.0%) | 25 (100.0%) | 10 (40.0%) | 20 (80.0%) | 25 (100.0%) | 25 (100.0%) | 25 (100.0%) | 25 (100.0%) |
| D4 | 25 | 0 | 6 (24.0%) | 25 (100.0%) | 13 (52.0%) | 22 (88.0%) | 25 (100.0%) | 25 (100.0%) | 25 (100.0%) | 25 (100.0%) |

**IPD-absent copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 25 | 4.2 / 3 | 0.0 / 0 | 4.2 / 3 | 0 |
| B | 25 | 3.5 / 1 | 0.0 / 0 | 3.5 / 1 | 0 |
| C4 | 25 | 3.6 / 2 | 0.0 / 0 | 3.6 / 2 | 0 |
| D4 | 25 | 3.1 / 1 | 0.0 / 0 | 3.1 / 1 | 0 |

**IPD-present copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 15 | 0 | 8 (53.3%) | 15 (100.0%) | 8 (53.3%) | 13 (86.7%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) |
| B | 15 | 0 | 8 (53.3%) | 15 (100.0%) | 8 (53.3%) | 13 (86.7%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) |
| C4 | 15 | 0 | 12 (80.0%) | 15 (100.0%) | 12 (80.0%) | 14 (93.3%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) |
| D4 | 15 | 0 | 13 (86.7%) | 15 (100.0%) | 13 (86.7%) | 14 (93.3%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) |

**IPD-present copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 15 | 1.5 / 0 | 0.0 / 0 | 1.5 / 0 | 0 |
| B | 15 | 1.6 / 0 | 0.0 / 0 | 1.6 / 0 | 0 |
| C4 | 15 | 0.7 / 0 | 0.0 / 0 | 0.7 / 0 | 0 |
| D4 | 15 | 0.5 / 0 | 0.0 / 0 | 0.5 / 0 | 0 |

**Sensitivity to the tie rule**

| arm | copies with tied alleles | gb exact (worst tied allele; primary) | gb exact (best tied allele) | mean gb wrong bases (worst) | mean gb wrong bases (best) |
|---|---|---|---|---|---|
| A | 1 | 8 (20.0%) | 8 (20.0%) | 3.17 | 3.12 |
| B | 1 | 12 (30.0%) | 12 (30.0%) | 2.77 | 2.73 |
| C4 | 1 | 15 (37.5%) | 15 (37.5%) | 2.50 | 2.45 |
| D4 | 1 | 19 (47.5%) | 19 (47.5%) | 2.12 | 2.08 |

**Paired contrasts**

| contrast | subset | metric | n paired | values (a vs b) | difference | 95% CI (family bootstrap) | sign-flip p |
|---|---|---|---|---|---|---|---|
| C4 - A | all copies | share gene body exact (pp) | 40 | 15 vs 8 (+7/-0) | +17.50 pp | [+7.50, +27.50] | 0.0156 |
| C4 - A | all copies | share coding exact (pp) | 40 | 40 vs 40 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| C4 - A | all copies | mean wrong bases, gene body | 40 | 2.50 vs 3.17 (closer 10/farther 0) | -0.68 | [-1.15, -0.30] | 0.00391 |
| C4 - A | all copies | mean wrong bases, coding | 40 | 0.00 vs 0.00 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| C4 - A | IPD-absent | share gene body exact (pp) | 25 | 3 vs 0 (+3/-0) | +12.00 pp | [+0.00, +25.00] | 0.25 |
| C4 - A | IPD-absent | share coding exact (pp) | 25 | 25 vs 25 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| C4 - A | IPD-absent | mean wrong bases, gene body | 25 | 3.60 vs 4.20 (closer 6/farther 0) | -0.60 | [-1.08, -0.17] | 0.0312 |
| C4 - A | IPD-absent | mean wrong bases, coding | 25 | 0.00 vs 0.00 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| C4 - A | IPD-present | share gene body exact (pp) | 15 | 12 vs 8 (+4/-0) | +26.67 pp | [+7.14, +47.06] | 0.125 |
| C4 - A | IPD-present | share coding exact (pp) | 15 | 15 vs 15 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| C4 - A | IPD-present | mean wrong bases, gene body | 15 | 0.67 vs 1.47 (closer 4/farther 0) | -0.80 | [-1.71, -0.14] | 0.125 |
| C4 - A | IPD-present | mean wrong bases, coding | 15 | 0.00 vs 0.00 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| D4 - B | all copies | share gene body exact (pp) | 40 | 19 vs 12 (+7/-0) | +17.50 pp | [+7.50, +27.50] | 0.0156 |
| D4 - B | all copies | share coding exact (pp) | 40 | 40 vs 40 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| D4 - B | all copies | mean wrong bases, gene body | 40 | 2.12 vs 2.77 (closer 9/farther 0) | -0.65 | [-1.15, -0.25] | 0.00781 |
| D4 - B | all copies | mean wrong bases, coding | 40 | 0.00 vs 0.00 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| D4 - B | IPD-absent | share gene body exact (pp) | 25 | 6 vs 4 (+2/-0) | +8.00 pp | [+0.00, +20.00] | 0.5 |
| D4 - B | IPD-absent | share coding exact (pp) | 25 | 25 vs 25 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| D4 - B | IPD-absent | mean wrong bases, gene body | 25 | 3.08 vs 3.48 (closer 4/farther 0) | -0.40 | [-0.88, -0.04] | 0.125 |
| D4 - B | IPD-absent | mean wrong bases, coding | 25 | 0.00 vs 0.00 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| D4 - B | IPD-present | share gene body exact (pp) | 15 | 13 vs 8 (+5/-0) | +33.33 pp | [+13.33, +57.14] | 0.0625 |
| D4 - B | IPD-present | share coding exact (pp) | 15 | 15 vs 15 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| D4 - B | IPD-present | mean wrong bases, gene body | 15 | 0.53 vs 1.60 (closer 5/farther 0) | -1.07 | [-2.00, -0.29] | 0.0625 |
| D4 - B | IPD-present | mean wrong bases, coding | 15 | 0.00 vs 0.00 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| B - A | all copies | share gene body exact (pp) | 40 | 12 vs 8 (+4/-0) | +10.00 pp | [+2.50, +20.00] | 0.125 |
| B - A | all copies | mean wrong bases, gene body | 40 | 2.77 vs 3.17 (closer 7/farther 2) | -0.40 | [-0.93, +0.23] | 0.219 |
| B - A | IPD-absent | share gene body exact (pp) | 25 | 4 vs 0 (+4/-0) | +16.00 pp | [+3.85, +30.43] | 0.125 |
| B - A | IPD-absent | mean wrong bases, gene body | 25 | 3.48 vs 4.20 (closer 7/farther 1) | -0.72 | [-1.59, +0.21] | 0.156 |
| B - A | IPD-present | share gene body exact (pp) | 15 | 8 vs 8 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| B - A | IPD-present | mean wrong bases, gene body | 15 | 1.60 vs 1.47 (closer 0/farther 1) | +0.13 | [+0.00, +0.43] | 1 |

**C4 final novel (patched) calls vs the A call for the same copy (gene-body distance to truth)**

| read set | copy type | copies with novel final call | comparable | not comparable (arm A call missing) | closer than A | equal | farther |
|---|---|---|---|---|---|---|---|
| D-real | IPD-absent | 6 | 6 | 0 | 6 | 0 | 0 |
| D-real | IPD-present | 0 | 0 | 0 | 0 | 0 | 0 |

Change in gene-body wrong bases (C4 minus A); negative = patched call is closer:

| read set | copy type | n | min | median | max | <=-20 | -19..-5 | -4..-1 | 0 | +1..+4 | +5..+19 | >=+20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D-real | IPD-absent | 6 | -5 | -2.0 | -1 | 0 | 1 | 5 | 0 | 0 | 0 | 0 |

**D4 final novel (patched) calls vs the B call for the same copy (gene-body distance to truth)**

| read set | copy type | copies with novel final call | comparable | not comparable (arm B call missing) | closer than B | equal | farther |
|---|---|---|---|---|---|---|---|
| D-real | IPD-absent | 12 | 12 | 0 | 4 | 8 | 0 |
| D-real | IPD-present | 0 | 0 | 0 | 0 | 0 | 0 |

Change in gene-body wrong bases (D4 minus B); negative = patched call is closer:

| read set | copy type | n | min | median | max | <=-20 | -19..-5 | -4..-1 | 0 | +1..+4 | +5..+19 | >=+20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D-real | IPD-absent | 12 | -5 | 0.0 | 0 | 0 | 1 | 3 | 8 | 0 | 0 | 0 |

## Validation against the frozen scorer (DRB1 only)

| arm | IPD-absent copies | frozen recovered_final | frozen recovered (any slot) | mine whole-seq exact (worst tied) | mine whole-seq exact (any tied) | mine gene-body exact | IPD-present copies with one 4-field truth | frozen named4 correct | mine whole-seq exact | mine gene-body exact | all IPD-present copies | mine gene-body exact (all present) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 25 | 0 | 0 | 0 | 0 | 0 | 15 | 8 | 8 | 8 | 15 | 8 |
| B | 25 | 4 | 4 | 4 | 4 | 4 | 15 | 8 | 8 | 8 | 15 | 8 |
| C4 | 25 | 2 | 3 | 2 | 2 | 3 | 15 | 12 | 12 | 12 | 15 | 12 |
| D4 | 25 | 6 | 6 | 6 | 6 | 6 | 15 | 13 | 13 | 13 | 15 | 13 |

Per-copy disagreements (frozen four-field correct vs gene-body identical on IPD-present copies with one 4-field truth; frozen recovered_final vs whole-sequence identical on IPD-absent copies): count by arm and kind

- none

First disagreements (arm, donor, gene, copy, kind, frozen, mine, n_tied): 

Gene-body recovery on IPD-absent DRB1 copies (frozen "gene body exact" column vs mine):

| arm | absent copies | frozen gene body exact (candidate in ANY slot) | same, restricted to FINAL calls (any slot) | mine gene-body exact (assigned slot, worst tied) | mine (any tied) |
|---|---|---|---|---|---|
| A | 25 | 0 | 0 | 0 | 0 |
| B | 25 | 4 | 4 | 4 | 4 |
| C4 | 25 | 4 | 3 | 3 | 3 |
| D4 | 25 | 6 | 6 | 6 | 6 |
