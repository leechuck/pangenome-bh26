# DRB1 patch-builder, read set D-sim (DEVELOPMENT), partial-credit re-scoring

40 true DRB1 copies per arm.

**all copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 40 | 0 | 10 (25.0%) | 39 (97.5%) | 19 (47.5%) | 35 (87.5%) | 40 (100.0%) | 39 (97.5%) | 40 (100.0%) | 40 (100.0%) |
| C4 | 40 | 0 | 27 (67.5%) | 39 (97.5%) | 30 (75.0%) | 37 (92.5%) | 40 (100.0%) | 39 (97.5%) | 40 (100.0%) | 40 (100.0%) |

**all copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 40 | 2.5 / 2 | 0.1 / 0 | 2.5 / 1.5 | 0 |
| C4 | 40 | 1.2 / 0 | 0.1 / 0 | 1.2 / 0 | 0 |

**IPD-absent copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 25 | 0 | 0 (0.0%) | 24 (96.0%) | 9 (36.0%) | 21 (84.0%) | 25 (100.0%) | 24 (96.0%) | 25 (100.0%) | 25 (100.0%) |
| C4 | 25 | 0 | 14 (56.0%) | 24 (96.0%) | 17 (68.0%) | 23 (92.0%) | 25 (100.0%) | 24 (96.0%) | 25 (100.0%) | 25 (100.0%) |

**IPD-absent copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 25 | 3.4 / 3 | 0.1 / 0 | 3.3 / 2 | 0 |
| C4 | 25 | 1.6 / 0 | 0.1 / 0 | 1.5 / 0 | 0 |

**IPD-present copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| A | 15 | 0 | 10 (66.7%) | 15 (100.0%) | 10 (66.7%) | 14 (93.3%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) |
| C4 | 15 | 0 | 13 (86.7%) | 15 (100.0%) | 13 (86.7%) | 14 (93.3%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) | 15 (100.0%) |

**IPD-present copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| A | 15 | 1.1 / 0 | 0.0 / 0 | 1.1 / 0 | 0 |
| C4 | 15 | 0.7 / 0 | 0.0 / 0 | 0.7 / 0 | 0 |

**Sensitivity to the tie rule**

| arm | copies with tied alleles | gb exact (worst tied allele; primary) | gb exact (best tied allele) | mean gb wrong bases (worst) | mean gb wrong bases (best) |
|---|---|---|---|---|---|
| A | 1 | 10 (25.0%) | 10 (25.0%) | 2.50 | 2.45 |
| C4 | 1 | 27 (67.5%) | 27 (67.5%) | 1.25 | 1.20 |

**Paired contrasts**

| contrast | subset | metric | n paired | values (a vs b) | difference | 95% CI (family bootstrap) | sign-flip p |
|---|---|---|---|---|---|---|---|
| C4 - A | all copies | share gene body exact (pp) | 40 | 27 vs 10 (+17/-0) | +42.50 pp | [+35.00, +50.00] | 1.53e-05 |
| C4 - A | all copies | share coding exact (pp) | 40 | 39 vs 39 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| C4 - A | all copies | mean wrong bases, gene body | 40 | 1.25 vs 2.50 (closer 20/farther 0) | -1.25 | [-1.77, -0.80] | 1.53e-05 |
| C4 - A | all copies | mean wrong bases, coding | 40 | 0.05 vs 0.05 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| C4 - A | IPD-absent | share gene body exact (pp) | 25 | 14 vs 0 (+14/-0) | +56.00 pp | [+40.00, +72.73] | 0.000122 |
| C4 - A | IPD-absent | share coding exact (pp) | 25 | 24 vs 24 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| C4 - A | IPD-absent | mean wrong bases, gene body | 25 | 1.60 vs 3.36 (closer 17/farther 0) | -1.76 | [-2.68, -1.04] | 6.1e-05 |
| C4 - A | IPD-absent | mean wrong bases, coding | 25 | 0.08 vs 0.08 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| C4 - A | IPD-present | share gene body exact (pp) | 15 | 13 vs 10 (+3/-0) | +20.00 pp | [+0.00, +40.00] | 0.25 |
| C4 - A | IPD-present | share coding exact (pp) | 15 | 15 vs 15 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| C4 - A | IPD-present | mean wrong bases, gene body | 15 | 0.67 vs 1.07 (closer 3/farther 0) | -0.40 | [-0.80, +0.00] | 0.25 |
| C4 - A | IPD-present | mean wrong bases, coding | 15 | 0.00 vs 0.00 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |

**C4 final novel (patched) calls vs the A call for the same copy (gene-body distance to truth)**

| read set | copy type | copies with novel final call | comparable | not comparable (arm A call missing) | closer than A | equal | farther |
|---|---|---|---|---|---|---|---|
| D-sim | IPD-absent | 17 | 17 | 0 | 17 | 0 | 0 |
| D-sim | IPD-present | 0 | 0 | 0 | 0 | 0 | 0 |

Change in gene-body wrong bases (C4 minus A); negative = patched call is closer:

| read set | copy type | n | min | median | max | <=-20 | -19..-5 | -4..-1 | 0 | +1..+4 | +5..+19 | >=+20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| D-sim | IPD-absent | 17 | -8 | -2 | -1 | 0 | 3 | 14 | 0 | 0 | 0 | 0 |

**D4 final novel (patched) calls vs the B call for the same copy (gene-body distance to truth)**

| read set | copy type | copies with novel final call | comparable | not comparable (arm B call missing) | closer than B | equal | farther |
|---|---|---|---|---|---|---|---|

Change in gene-body wrong bases (D4 minus B); negative = patched call is closer:

| read set | copy type | n | min | median | max | <=-20 | -19..-5 | -4..-1 | 0 | +1..+4 | +5..+19 | >=+20 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|

## Validation against the frozen scorer (DRB1 only)

| arm | IPD-absent copies | frozen recovered_final | frozen recovered (any slot) | mine whole-seq exact (worst tied) | mine whole-seq exact (any tied) | mine gene-body exact | IPD-present copies with one 4-field truth | frozen named4 correct | mine whole-seq exact | mine gene-body exact | all IPD-present copies | mine gene-body exact (all present) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A | 25 | 0 | 0 | 0 | 0 | 0 | 15 | 10 | 10 | 10 | 15 | 10 |
| C4 | 25 | 6 | 7 | 6 | 6 | 14 | 15 | 13 | 13 | 13 | 15 | 13 |

Per-copy disagreements (frozen four-field correct vs gene-body identical on IPD-present copies with one 4-field truth; frozen recovered_final vs whole-sequence identical on IPD-absent copies): count by arm and kind

- none

First disagreements (arm, donor, gene, copy, kind, frozen, mine, n_tied): 

Gene-body recovery on IPD-absent DRB1 copies (frozen "gene body exact" column vs mine):

| arm | absent copies | frozen gene body exact (candidate in ANY slot) | same, restricted to FINAL calls (any slot) | mine gene-body exact (assigned slot, worst tied) | mine (any tied) |
|---|---|---|---|---|---|
| A | 25 | 0 | 0 | 0 | 0 |
| C4 | 25 | 15 | 14 | 14 | 14 |
