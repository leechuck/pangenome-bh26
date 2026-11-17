# Confirmatory test, real reads (partial-credit re-scoring)

## Set P: 16 donors, 254 true copies per arm

**all copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| ipd | 254 | 0 | NA | 243 (95.7%) | NA | NA | NA | 245 (96.5%) | 246 (96.9%) | 246 (96.9%) |
| genomic | 254 | 0 | 219 (86.2%) | 251 (98.8%) | 231 (90.9%) | 249 (98.0%) | 254 (100.0%) | 254 (100.0%) | 254 (100.0%) | 254 (100.0%) |
| genomic_core_hprc_asian | 254 | 0 | 230 (90.6%) | 253 (99.6%) | 238 (93.7%) | 249 (98.0%) | 254 (100.0%) | 254 (100.0%) | 254 (100.0%) | 254 (100.0%) |
| genomic_core_nonasian | 254 | 0 | 213 (83.9%) | 249 (98.0%) | 229 (90.2%) | 250 (98.4%) | 254 (100.0%) | 253 (99.6%) | 254 (100.0%) | 254 (100.0%) |

**all copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| ipd | 254 | NA | 2.2 / 0 | NA | 0 |
| genomic | 254 | 0.5 / 0 | 0.0 / 0 | 0.5 / 0 | 0 |
| genomic_core_hprc_asian | 254 | 0.5 / 0 | 0.0 / 0 | 0.5 / 0 | 0 |
| genomic_core_nonasian | 254 | 0.5 / 0 | 0.0 / 0 | 0.5 / 0 | 0 |

**IPD-absent copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| ipd | 35 | 0 | NA | 33 (94.3%) | NA | NA | NA | 34 (97.1%) | 35 (100.0%) | 35 (100.0%) |
| genomic | 35 | 0 | 8 (22.9%) | 33 (94.3%) | 17 (48.6%) | 30 (85.7%) | 35 (100.0%) | 35 (100.0%) | 35 (100.0%) | 35 (100.0%) |
| genomic_core_hprc_asian | 35 | 0 | 19 (54.3%) | 35 (100.0%) | 24 (68.6%) | 30 (85.7%) | 35 (100.0%) | 35 (100.0%) | 35 (100.0%) | 35 (100.0%) |
| genomic_core_nonasian | 35 | 0 | 8 (22.9%) | 32 (91.4%) | 15 (42.9%) | 31 (88.6%) | 35 (100.0%) | 34 (97.1%) | 35 (100.0%) | 35 (100.0%) |

**IPD-absent copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| ipd | 35 | NA | 0.1 / 0 | NA | 0 |
| genomic | 35 | 3.3 / 2 | 0.1 / 0 | 3.3 / 2 | 0 |
| genomic_core_hprc_asian | 35 | 2.8 / 0 | 0.0 / 0 | 2.8 / 0 | 0 |
| genomic_core_nonasian | 35 | 3.1 / 2 | 0.2 / 0 | 3.0 / 2 | 0 |

**IPD-present copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| ipd | 219 | 0 | NA | 210 (95.9%) | NA | NA | NA | 211 (96.3%) | 211 (96.3%) | 211 (96.3%) |
| genomic | 219 | 0 | 211 (96.3%) | 218 (99.5%) | 214 (97.7%) | 219 (100.0%) | 219 (100.0%) | 219 (100.0%) | 219 (100.0%) | 219 (100.0%) |
| genomic_core_hprc_asian | 219 | 0 | 211 (96.3%) | 218 (99.5%) | 214 (97.7%) | 219 (100.0%) | 219 (100.0%) | 219 (100.0%) | 219 (100.0%) | 219 (100.0%) |
| genomic_core_nonasian | 219 | 0 | 205 (93.6%) | 217 (99.1%) | 214 (97.7%) | 219 (100.0%) | 219 (100.0%) | 219 (100.0%) | 219 (100.0%) | 219 (100.0%) |

**IPD-present copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| ipd | 219 | NA | 2.6 / 0 | NA | 0 |
| genomic | 219 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |
| genomic_core_hprc_asian | 219 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |
| genomic_core_nonasian | 219 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |

**By gene**

| gene | ipd: copies; gb exact; mean gb / mean cod | genomic: copies; gb exact; mean gb / mean cod | genomic_core_hprc_asian: copies; gb exact; mean gb / mean cod | genomic_core_nonasian: copies; gb exact; mean gb / mean cod |
|---|---|---|---|---|
| A | 32; NA; NA / 0.0 | 32; 32; 0.0 / 0.0 | 32; 32; 0.0 / 0.0 | 32; 31; 0.0 / 0.0 |
| B | 32; NA; NA / 0.1 | 32; 31; 0.0 / 0.0 | 32; 31; 0.0 / 0.0 | 32; 31; 0.0 / 0.0 |
| C | 32; NA; NA / 0.0 | 32; 32; 0.0 / 0.0 | 32; 32; 0.0 / 0.0 | 32; 32; 0.0 / 0.0 |
| DPA1 | 32; NA; NA / 0.0 | 32; 32; 0.0 / 0.0 | 32; 32; 0.0 / 0.0 | 32; 32; 0.0 / 0.0 |
| DPB1 | 32; NA; NA / 0.0 | 32; 28; 0.4 / 0.0 | 32; 28; 0.4 / 0.0 | 32; 28; 0.4 / 0.0 |
| DQA1 | 32; NA; NA / 0.0 | 32; 23; 0.6 / 0.1 | 32; 27; 0.3 / 0.0 | 32; 20; 0.7 / 0.1 |
| DQB1 | 32; NA; NA / 17.6 | 32; 31; 0.0 / 0.0 | 32; 31; 0.0 / 0.0 | 32; 31; 0.0 / 0.0 |
| DRB1 | 30; NA; NA / 0.1 | 30; 10; 3.4 / 0.0 | 30; 17; 3.1 / 0.0 | 30; 8; 3.2 / 0.1 |

**Sensitivity to the tie rule (gene body, arms with a full gene body)**

| arm | copies with tied alleles | gb exact (worst tied allele; primary) | gb exact (best tied allele) | mean gb wrong bases (worst) | mean gb wrong bases (best) |
|---|---|---|---|---|---|
| genomic | 3 | 219 (86.2%) | 219 (86.2%) | 0.53 | 0.52 |
| genomic_core_hprc_asian | 9 | 230 (90.6%) | 230 (90.6%) | 0.46 | 0.41 |
| genomic_core_nonasian | 16 | 213 (83.9%) | 218 (85.8%) | 0.53 | 0.48 |

**Paired contrasts (set P); for exactness a - b in percentage points, for wrong bases a - b (negative = a is closer)**

| contrast | subset | metric | n paired | values (a vs b) | difference | 95% CI (family bootstrap) | sign-flip p |
|---|---|---|---|---|---|---|---|
| genomic - ipd | all copies | share coding exact (pp) | 254 | 251 vs 243 (+11/-3) | +3.15 pp | [+0.00, +6.30] | 0.117 |
| genomic - ipd | all copies | mean wrong bases, coding | 254 | 0.01 vs 2.24 (closer 11/farther 3) | -2.23 | [-4.42, -0.56] | 0.0547 |
| genomic - ipd | IPD-absent | share coding exact (pp) | 35 | 33 vs 33 (+2/-2) | +0.00 pp | [-8.57, +7.32] | 1 |
| genomic - ipd | IPD-absent | mean wrong bases, coding | 35 | 0.06 vs 0.11 (closer 2/farther 2) | -0.06 | [-0.23, +0.09] | 1 |
| genomic - ipd | IPD-present | share coding exact (pp) | 219 | 218 vs 210 (+9/-1) | +3.65 pp | [+0.45, +7.31] | 0.0938 |
| genomic - ipd | IPD-present | mean wrong bases, coding | 219 | 0.00 vs 2.58 (closer 9/farther 1) | -2.58 | [-5.17, -0.63] | 0.0938 |
| genomic_core_hprc_asian - genomic | all copies | share gene body exact (pp) | 254 | 230 vs 219 (+11/-0) | +4.33 pp | [+1.57, +7.87] | 0.0156 |
| genomic_core_hprc_asian - genomic | all copies | share coding exact (pp) | 254 | 253 vs 251 (+2/-0) | +0.79 pp | [+0.00, +1.97] | 0.5 |
| genomic_core_hprc_asian - genomic | all copies | mean wrong bases, gene body | 254 | 0.46 vs 0.53 (closer 13/farther 6) | -0.07 | [-0.31, +0.19] | 0.615 |
| genomic_core_hprc_asian - genomic | all copies | mean wrong bases, coding | 254 | 0.00 vs 0.01 (closer 2/farther 0) | -0.01 | [-0.02, +0.00] | 0.5 |
| genomic_core_hprc_asian - genomic | IPD-absent | share gene body exact (pp) | 35 | 19 vs 8 (+11/-0) | +31.43 pp | [+13.33, +48.65] | 0.0156 |
| genomic_core_hprc_asian - genomic | IPD-absent | share coding exact (pp) | 35 | 35 vs 33 (+2/-0) | +5.71 pp | [+0.00, +12.50] | 0.5 |
| genomic_core_hprc_asian - genomic | IPD-absent | mean wrong bases, gene body | 35 | 2.80 vs 3.34 (closer 13/farther 6) | -0.54 | [-2.00, +1.55] | 0.615 |
| genomic_core_hprc_asian - genomic | IPD-absent | mean wrong bases, coding | 35 | 0.00 vs 0.06 (closer 2/farther 0) | -0.06 | [-0.12, +0.00] | 0.5 |
| genomic_core_hprc_asian - genomic | IPD-present | share gene body exact (pp) | 219 | 211 vs 211 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| genomic_core_hprc_asian - genomic | IPD-present | share coding exact (pp) | 219 | 218 vs 218 (+0/-0) | +0.00 pp | [+0.00, +0.00] | 1 |
| genomic_core_hprc_asian - genomic | IPD-present | mean wrong bases, gene body | 219 | 0.08 vs 0.08 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| genomic_core_hprc_asian - genomic | IPD-present | mean wrong bases, coding | 219 | 0.00 vs 0.00 (closer 0/farther 0) | +0.00 | [+0.00, +0.00] | 1 |
| genomic_core_hprc_asian - genomic_core_nonasian | all copies | share gene body exact (pp) | 254 | 230 vs 213 (+17/-0) | +6.69 pp | [+3.57, +10.24] | 0.000977 |
| genomic_core_hprc_asian - genomic_core_nonasian | all copies | share coding exact (pp) | 254 | 253 vs 249 (+4/-0) | +1.57 pp | [+0.00, +3.94] | 0.5 |
| genomic_core_hprc_asian - genomic_core_nonasian | all copies | mean wrong bases, gene body | 254 | 0.46 vs 0.53 (closer 19/farther 5) | -0.07 | [-0.30, +0.18] | 0.619 |
| genomic_core_hprc_asian - genomic_core_nonasian | all copies | mean wrong bases, coding | 254 | 0.00 vs 0.03 (closer 4/farther 0) | -0.03 | [-0.07, +0.00] | 0.5 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-absent | share gene body exact (pp) | 35 | 19 vs 8 (+11/-0) | +31.43 pp | [+13.33, +48.65] | 0.0156 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-absent | share coding exact (pp) | 35 | 35 vs 32 (+3/-0) | +8.57 pp | [+0.00, +18.75] | 0.5 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-absent | mean wrong bases, gene body | 35 | 2.80 vs 3.14 (closer 13/farther 5) | -0.34 | [-1.73, +1.68] | 0.745 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-absent | mean wrong bases, coding | 35 | 0.00 vs 0.17 (closer 3/farther 0) | -0.17 | [-0.43, +0.00] | 0.5 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-present | share gene body exact (pp) | 219 | 211 vs 205 (+6/-0) | +2.74 pp | [+0.94, +4.52] | 0.0312 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-present | share coding exact (pp) | 219 | 218 vs 217 (+1/-0) | +0.46 pp | [+0.00, +1.41] | 1 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-present | mean wrong bases, gene body | 219 | 0.08 vs 0.11 (closer 6/farther 0) | -0.03 | [-0.05, -0.01] | 0.0312 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-present | mean wrong bases, coding | 219 | 0.00 vs 0.01 (closer 1/farther 0) | -0.00 | [-0.01, +0.00] | 1 |

## Set A: 33 donors, 526 true copies per arm

**all copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| ipd | 526 | 0 | NA | 495 (94.1%) | NA | NA | NA | 503 (95.6%) | 504 (95.8%) | 504 (95.8%) |
| genomic | 526 | 0 | 443 (84.2%) | 517 (98.3%) | 472 (89.7%) | 507 (96.4%) | 524 (99.6%) | 525 (99.8%) | 526 (100.0%) | 526 (100.0%) |
| genomic_core_hprc_asian | 526 | 0 | 457 (86.9%) | 519 (98.7%) | 484 (92.0%) | 510 (97.0%) | 524 (99.6%) | 524 (99.6%) | 526 (100.0%) | 526 (100.0%) |
| genomic_core_nonasian | 526 | 0 | 431 (81.9%) | 513 (97.5%) | 470 (89.4%) | 507 (96.4%) | 524 (99.6%) | 523 (99.4%) | 526 (100.0%) | 526 (100.0%) |

**all copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| ipd | 526 | NA | 2.7 / 0 | NA | 0 |
| genomic | 526 | 0.9 / 0 | 0.0 / 0 | 0.8 / 0 | 0 |
| genomic_core_hprc_asian | 526 | 0.8 / 0 | 0.0 / 0 | 0.8 / 0 | 0 |
| genomic_core_nonasian | 526 | 0.9 / 0 | 0.0 / 0 | 0.8 / 0 | 0 |

**IPD-absent copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| ipd | 81 | 0 | NA | 78 (96.3%) | NA | NA | NA | 79 (97.5%) | 80 (98.8%) | 80 (98.8%) |
| genomic | 81 | 0 | 15 (18.5%) | 76 (93.8%) | 34 (42.0%) | 62 (76.5%) | 79 (97.5%) | 81 (100.0%) | 81 (100.0%) | 81 (100.0%) |
| genomic_core_hprc_asian | 81 | 0 | 33 (40.7%) | 79 (97.5%) | 47 (58.0%) | 65 (80.2%) | 79 (97.5%) | 80 (98.8%) | 81 (100.0%) | 81 (100.0%) |
| genomic_core_nonasian | 81 | 0 | 15 (18.5%) | 74 (91.4%) | 32 (39.5%) | 62 (76.5%) | 79 (97.5%) | 79 (97.5%) | 81 (100.0%) | 81 (100.0%) |

**IPD-absent copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| ipd | 81 | NA | 0.4 / 0 | NA | 0 |
| genomic | 81 | 5.2 / 2 | 0.1 / 0 | 5.2 / 2 | 0 |
| genomic_core_hprc_asian | 81 | 4.7 / 1 | 0.0 / 0 | 4.7 / 1 | 0 |
| genomic_core_nonasian | 81 | 5.0 / 2 | 0.1 / 0 | 4.9 / 2 | 0 |

**IPD-present copies: shares (denominator = all copies; a missing call counts as not correct)**

| arm | copies | missing calls | gene body exact | coding exact | gb <=1 | gb <=5 | gb <=20 | cod <=1 | cod <=5 | cod <=20 |
|---|---|---|---|---|---|---|---|---|---|---|
| ipd | 445 | 0 | NA | 417 (93.7%) | NA | NA | NA | 424 (95.3%) | 424 (95.3%) | 424 (95.3%) |
| genomic | 445 | 0 | 428 (96.2%) | 441 (99.1%) | 438 (98.4%) | 445 (100.0%) | 445 (100.0%) | 444 (99.8%) | 445 (100.0%) | 445 (100.0%) |
| genomic_core_hprc_asian | 445 | 0 | 424 (95.3%) | 440 (98.9%) | 437 (98.2%) | 445 (100.0%) | 445 (100.0%) | 444 (99.8%) | 445 (100.0%) | 445 (100.0%) |
| genomic_core_nonasian | 445 | 0 | 416 (93.5%) | 439 (98.7%) | 438 (98.4%) | 445 (100.0%) | 445 (100.0%) | 444 (99.8%) | 445 (100.0%) | 445 (100.0%) |

**IPD-present copies: wrong bases per called copy, mean / median (missing calls excluded)**

| arm | called copies | gene body | coding | non-coding | nc negative (clipped to 0) |
|---|---|---|---|---|---|
| ipd | 445 | NA | 3.2 / 0 | NA | 0 |
| genomic | 445 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |
| genomic_core_hprc_asian | 445 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |
| genomic_core_nonasian | 445 | 0.1 / 0 | 0.0 / 0 | 0.1 / 0 | 0 |

**By gene**

| gene | ipd: copies; gb exact; mean gb / mean cod | genomic: copies; gb exact; mean gb / mean cod | genomic_core_hprc_asian: copies; gb exact; mean gb / mean cod | genomic_core_nonasian: copies; gb exact; mean gb / mean cod |
|---|---|---|---|---|
| A | 66; NA; NA / 0.0 | 66; 65; 0.0 / 0.0 | 66; 65; 0.0 / 0.0 | 66; 64; 0.0 / 0.0 |
| B | 66; NA; NA / 0.1 | 66; 64; 0.1 / 0.1 | 66; 64; 0.1 / 0.1 | 66; 64; 0.1 / 0.1 |
| C | 66; NA; NA / 0.0 | 66; 65; 0.0 / 0.0 | 66; 65; 0.0 / 0.0 | 66; 65; 0.0 / 0.0 |
| DPA1 | 66; NA; NA / 0.0 | 66; 65; 1.6 / 0.0 | 66; 65; 1.6 / 0.0 | 66; 65; 1.6 / 0.0 |
| DPB1 | 66; NA; NA / 0.0 | 66; 57; 0.3 / 0.0 | 66; 57; 0.3 / 0.0 | 66; 56; 0.3 / 0.0 |
| DQA1 | 66; NA; NA / 0.0 | 66; 48; 0.6 / 0.1 | 66; 52; 0.4 / 0.0 | 66; 42; 0.7 / 0.1 |
| DQB1 | 66; NA; NA / 21.7 | 66; 61; 0.1 / 0.0 | 66; 61; 0.1 / 0.0 | 66; 61; 0.1 / 0.0 |
| DRB1 | 64; NA; NA / 0.0 | 64; 18; 4.4 / 0.0 | 64; 28; 4.1 / 0.0 | 64; 14; 4.2 / 0.1 |

**Sensitivity to the tie rule (gene body, arms with a full gene body)**

| arm | copies with tied alleles | gb exact (worst tied allele; primary) | gb exact (best tied allele) | mean gb wrong bases (worst) | mean gb wrong bases (best) |
|---|---|---|---|---|---|
| genomic | 3 | 443 (84.2%) | 443 (84.2%) | 0.87 | 0.86 |
| genomic_core_hprc_asian | 23 | 457 (86.9%) | 459 (87.3%) | 0.80 | 0.69 |
| genomic_core_nonasian | 27 | 431 (81.9%) | 441 (83.8%) | 0.86 | 0.83 |

**Paired contrasts (set A); for exactness a - b in percentage points, for wrong bases a - b (negative = a is closer)**

| contrast | subset | metric | n paired | values (a vs b) | difference | 95% CI (family bootstrap) | sign-flip p |
|---|---|---|---|---|---|---|---|
| genomic - ipd | all copies | share coding exact (pp) | 526 | 517 vs 495 (+31/-9) | +4.18 pp | [+1.79, +6.28] | 0.00701 |
| genomic - ipd | all copies | mean wrong bases, coding | 526 | 0.02 vs 2.74 (closer 30/farther 9) | -2.71 | [-4.14, -1.52] | 0.00125 |
| genomic - ipd | IPD-absent | share coding exact (pp) | 81 | 76 vs 78 (+3/-5) | -2.47 pp | [-8.64, +3.17] | 0.688 |
| genomic - ipd | IPD-absent | mean wrong bases, coding | 81 | 0.06 vs 0.42 (closer 3/farther 5) | -0.36 | [-1.19, +0.08] | 0.656 |
| genomic - ipd | IPD-present | share coding exact (pp) | 445 | 441 vs 417 (+28/-4) | +5.39 pp | [+2.67, +7.96] | 0.00281 |
| genomic - ipd | IPD-present | mean wrong bases, coding | 445 | 0.02 vs 3.16 (closer 27/farther 4) | -3.14 | [-4.88, -1.70] | 0.00128 |
| genomic_core_hprc_asian - genomic | all copies | share gene body exact (pp) | 526 | 457 vs 443 (+18/-4) | +2.66 pp | [+0.76, +5.32] | 0.0142 |
| genomic_core_hprc_asian - genomic | all copies | share coding exact (pp) | 526 | 519 vs 517 (+4/-2) | +0.38 pp | [-0.26, +1.34] | 0.625 |
| genomic_core_hprc_asian - genomic | all copies | mean wrong bases, gene body | 526 | 0.80 vs 0.87 (closer 25/farther 15) | -0.06 | [-0.21, +0.08] | 0.43 |
| genomic_core_hprc_asian - genomic | all copies | mean wrong bases, coding | 526 | 0.02 vs 0.02 (closer 4/farther 2) | -0.00 | [-0.01, +0.01] | 1 |
| genomic_core_hprc_asian - genomic | IPD-absent | share gene body exact (pp) | 81 | 33 vs 15 (+18/-0) | +22.22 pp | [+12.09, +34.29] | 0.000244 |
| genomic_core_hprc_asian - genomic | IPD-absent | share coding exact (pp) | 81 | 79 vs 76 (+4/-1) | +3.70 pp | [+0.00, +8.06] | 0.25 |
| genomic_core_hprc_asian - genomic | IPD-absent | mean wrong bases, gene body | 81 | 4.74 vs 5.23 (closer 25/farther 11) | -0.49 | [-1.30, +0.48] | 0.312 |
| genomic_core_hprc_asian - genomic | IPD-absent | mean wrong bases, coding | 81 | 0.04 vs 0.06 (closer 4/farther 1) | -0.02 | [-0.08, +0.02] | 0.625 |
| genomic_core_hprc_asian - genomic | IPD-present | share gene body exact (pp) | 445 | 424 vs 428 (+0/-4) | -0.90 pp | [-1.64, +0.00] | 0.25 |
| genomic_core_hprc_asian - genomic | IPD-present | share coding exact (pp) | 445 | 440 vs 441 (+0/-1) | -0.22 pp | [-0.51, +0.00] | 1 |
| genomic_core_hprc_asian - genomic | IPD-present | mean wrong bases, gene body | 445 | 0.09 vs 0.07 (closer 0/farther 4) | +0.02 | [+0.00, +0.04] | 0.25 |
| genomic_core_hprc_asian - genomic | IPD-present | mean wrong bases, coding | 445 | 0.02 vs 0.02 (closer 0/farther 1) | +0.00 | [+0.00, +0.01] | 1 |
| genomic_core_hprc_asian - genomic_core_nonasian | all copies | share gene body exact (pp) | 526 | 457 vs 431 (+27/-1) | +4.94 pp | [+3.14, +7.14] | 3.05e-05 |
| genomic_core_hprc_asian - genomic_core_nonasian | all copies | share coding exact (pp) | 526 | 519 vs 513 (+7/-1) | +1.14 pp | [+0.19, +2.48] | 0.125 |
| genomic_core_hprc_asian - genomic_core_nonasian | all copies | mean wrong bases, gene body | 526 | 0.80 vs 0.86 (closer 34/farther 13) | -0.05 | [-0.20, +0.10] | 0.494 |
| genomic_core_hprc_asian - genomic_core_nonasian | all copies | mean wrong bases, coding | 526 | 0.02 vs 0.04 (closer 7/farther 1) | -0.02 | [-0.05, +0.00] | 0.188 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-absent | share gene body exact (pp) | 81 | 33 vs 15 (+18/-0) | +22.22 pp | [+12.09, +34.29] | 0.000244 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-absent | share coding exact (pp) | 81 | 79 vs 74 (+6/-1) | +6.17 pp | [+1.20, +12.50] | 0.125 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-absent | mean wrong bases, gene body | 81 | 4.74 vs 5.04 (closer 25/farther 11) | -0.30 | [-1.15, +0.70] | 0.561 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-absent | mean wrong bases, coding | 81 | 0.04 vs 0.15 (closer 6/farther 1) | -0.11 | [-0.27, +0.01] | 0.25 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-present | share gene body exact (pp) | 445 | 424 vs 416 (+9/-1) | +1.80 pp | [+0.28, +2.91] | 0.0625 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-present | share coding exact (pp) | 445 | 440 vs 439 (+1/-0) | +0.22 pp | [+0.00, +0.81] | 1 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-present | mean wrong bases, gene body | 445 | 0.09 vs 0.10 (closer 9/farther 2) | -0.01 | [-0.02, +0.00] | 0.289 |
| genomic_core_hprc_asian - genomic_core_nonasian | IPD-present | mean wrong bases, coding | 445 | 0.02 vs 0.02 (closer 1/farther 0) | -0.00 | [-0.01, +0.00] | 1 |

## Validation against the frozen scorer

**Set P: frozen scorer vs this scorer (distance = 0)**

| arm | frozen IPD-absent recovered / copies | mine: whole-sequence identical | mine: gene body identical | frozen four-field correct / eligible donor-genes | mine: both copies whole-sequence identical (core differs from IPD allele ends in many copies) | mine: both copies gene body identical | mine: both copies coding identical |
|---|---|---|---|---|---|---|---|
| genomic | 0/35 | 0 | 8 | 93/99 | 84 | 93 | 98 |
| genomic_core_hprc_asian | 14/35 | 14 | 19 | 93/99 | 81 | 93 | 98 |
| genomic_core_nonasian | 2/35 | 2 | 8 | 88/99 | 79 | 90 | 97 |
| ipd | 0/35 | 0 | 0 | 53/99 | 0 | NA | 94 |

ipd arm: frozen four-field endpoint counts a donor-gene as called only if both slots carry one four-field name: 92 of the eligible donor-genes were called that way but far fewer match the truth name; the ipd reference has no or little intron/UTR content, so four-field (non-coding) differences cannot be resolved there. Coding-identical is a different, weaker test and is not expected to match the frozen count.

Gene-body-identical vs frozen four-field correct, mismatches at donor-gene level, counts by arm and kind: genomic_core_nonasian/frozen_not_called(tie/non-4-field): 2

First mismatches (arm, donor, gene, kind, frozen, mine): ('genomic_core_nonasian', 'HG03792', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2); ('genomic_core_nonasian', 'HG03919', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2)

**Set A: frozen scorer vs this scorer (distance = 0)**

| arm | frozen IPD-absent recovered / copies | mine: whole-sequence identical | mine: gene body identical | frozen four-field correct / eligible donor-genes | mine: both copies whole-sequence identical (core differs from IPD allele ends in many copies) | mine: both copies gene body identical | mine: both copies coding identical |
|---|---|---|---|---|---|---|---|
| genomic | 0/81 | 0 | 15 | 185/197 | 163 | 185 | 194 |
| genomic_core_hprc_asian | 21/81 | 21 | 33 | 182/197 | 156 | 184 | 194 |
| genomic_core_nonasian | 2/81 | 2 | 15 | 173/197 | 152 | 179 | 193 |
| ipd | 0/81 | 0 | 0 | 110/197 | 0 | NA | 181 |

ipd arm: frozen four-field endpoint counts a donor-gene as called only if both slots carry one four-field name: 181 of the eligible donor-genes were called that way but far fewer match the truth name; the ipd reference has no or little intron/UTR content, so four-field (non-coding) differences cannot be resolved there. Coding-identical is a different, weaker test and is not expected to match the frozen count.

Gene-body-identical vs frozen four-field correct, mismatches at donor-gene level, counts by arm and kind: genomic_core_hprc_asian/frozen_not_called(tie/non-4-field): 2, genomic_core_nonasian/frozen_not_called(tie/non-4-field): 6

First mismatches (arm, donor, gene, kind, frozen, mine): ('genomic_core_hprc_asian', 'HG03882', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2); ('genomic_core_hprc_asian', 'HG04153', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2); ('genomic_core_nonasian', 'HG03792', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2); ('genomic_core_nonasian', 'HG03868', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2); ('genomic_core_nonasian', 'HG03919', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2); ('genomic_core_nonasian', 'HG04144', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2); ('genomic_core_nonasian', 'HG04216', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2); ('genomic_core_nonasian', 'NA20858', 'C', 'frozen_not_called(tie/non-4-field)', '0', 1, 2)

Truth copies without Immuannot coding exons (coding/gene-body undefined): 0
Rows with unresolved call labels: 0
Gene-body defined share per arm: ipd 0/526, genomic 526/526, genomic_core_hprc_asian 526/526, genomic_core_nonasian 526/526
