# Endpoints, tag test

**f. fallbacks:** 0 of 254 tasks (PASS; limit 2). Gene-level non-complete by kind: none. Detail: none


## T-sim

**a. novel final calls on IPD-present copies (all genes):** {'C4': 1, 'C4R4n': 1, 'C4R1n': 1, 'D4': 1, 'D4R4n': 1, 'D4R1n': 1} -> PASS (primary arms <= 1)

**b. four-field named accuracy:** PASS

- C4R4n vs A: total 707 -> 723 of 804 (+16 of 804; +2.0 pp [+1.1, +3.0], p=0.000488, 50 clusters); per-gene change {'A': 1, 'B': 1, 'C': 0, 'DPA1': 2, 'DPB1': 0, 'DQA1': 7, 'DQB1': 1, 'DRB1': 4}
- D4R4n vs B: total 694 -> 711 of 804 (+17 of 804; +2.1 pp [+1.2, +3.2], p=0.000244, 50 clusters); per-gene change {'A': 1, 'B': 1, 'C': 0, 'DPA1': 2, 'DPB1': 0, 'DQA1': 8, 'DQB1': 1, 'DRB1': 4}

**c. gene-body recovery of IPD-absent copies, C4R4n minus A, all genes:** 0 -> 22 of 162; +22 of 162; +13.6 pp [+8.8, +18.7], p=1.53e-05, 49 clusters -> PASS

**c. gene-body recovery of IPD-absent copies, C4R4n minus A, non-DRB1 genes (secondary):** 0 -> 9 of 86; +9 of 86; +10.5 pp [+5.0, +16.1], p=0.00391, 41 clusters

**d. gene-body recovery of IPD-absent copies, D4R4n minus B, all genes:** 39 -> 56 of 162; +17 of 162; +10.5 pp [+6.1, +15.1], p=0.000244, 49 clusters

**d. gene-body recovery of IPD-absent copies, D4R4n minus B, non-DRB1 genes (secondary):** 22 -> 31 of 86; +9 of 86; +10.5 pp [+5.2, +15.8], p=0.00391, 41 clusters

**e. novel final calls on IPD-absent copies: right gene body / not**
- C4: 18/34 = 53% right; C4R4n: 15/23 = 65% right; right calls kept 15/18 -> PASS
- D4: 12/22 = 55% right; D4R4n: 10/14 = 71% right; right calls kept 10/12 -> PASS

**Secondary: per gene (T-sim)**

| gene | arm | IPD-absent | body | byte-exact final | 4-field named | novel on present | novel absent right/wrong | mean wrong bases per copy | closer/equal/farther |
|---|---|---|---|---|---|---|---|---|---|
| A | A | 2 | 0 | 0 | 112/119 | 0 | 0/0 | 0.80 (2 missing) |  |
| A | B | 2 | 0 | 0 | 112/119 | 0 | 0/0 | 0.80 (2 missing) |  |
| A | C4 | 2 | 0 | 0 | 113/119 | 0 | 0/1 | 0.72 (2 missing) | 0/0/0 |
| A | C4R4n | 2 | 0 | 0 | 113/119 | 0 | 0/0 | 0.79 (2 missing) | 0/0/0 |
| A | C4R1n | 2 | 0 | 0 | 113/119 | 0 | 0/1 | 0.72 (2 missing) | 0/0/0 |
| A | D4 | 2 | 0 | 0 | 113/119 | 0 | 0/1 | 0.73 (2 missing) | 0/0/0 |
| A | D4R4n | 2 | 0 | 0 | 113/119 | 0 | 0/0 | 0.79 (2 missing) | 0/0/0 |
| A | D4R1n | 2 | 0 | 0 | 113/119 | 0 | 0/0 | 0.79 (2 missing) | 0/0/0 |
| B | A | 12 | 0 | 0 | 105/109 | 0 | 0/0 | 0.05 |  |
| B | B | 12 | 0 | 0 | 102/109 | 0 | 0/0 | 0.05 |  |
| B | C4 | 12 | 0 | 0 | 106/109 | 0 | 0/0 | 0.03 | 0/0/0 |
| B | C4R4n | 12 | 0 | 0 | 106/109 | 0 | 0/0 | 0.03 | 0/0/0 |
| B | C4R1n | 12 | 0 | 0 | 106/109 | 0 | 0/0 | 0.03 | 0/0/0 |
| B | D4 | 12 | 0 | 0 | 103/109 | 0 | 0/0 | 0.03 | 0/0/0 |
| B | D4R4n | 12 | 0 | 0 | 103/109 | 0 | 0/0 | 0.03 | 0/0/0 |
| B | D4R1n | 12 | 0 | 0 | 103/109 | 0 | 0/0 | 0.03 | 0/0/0 |
| C | A | 6 | 0 | 0 | 107/115 | 0 | 0/0 | 1.05 |  |
| C | B | 6 | 0 | 0 | 104/115 | 0 | 0/0 | 1.03 |  |
| C | C4 | 6 | 1 | 0 | 107/115 | 0 | 0/0 | 1.05 | 0/0/0 |
| C | C4R4n | 6 | 1 | 0 | 107/115 | 0 | 0/0 | 1.05 | 0/0/0 |
| C | C4R1n | 6 | 1 | 0 | 107/115 | 0 | 0/0 | 1.05 | 0/0/0 |
| C | D4 | 6 | 1 | 0 | 104/115 | 0 | 0/0 | 1.03 | 0/0/0 |
| C | D4R4n | 6 | 1 | 0 | 104/115 | 0 | 0/0 | 1.03 | 0/0/0 |
| C | D4R1n | 6 | 1 | 0 | 104/115 | 0 | 0/0 | 1.03 | 0/0/0 |
| DPA1 | A | 16 | 0 | 0 | 71/104 | 0 | 0/0 | 0.81 |  |
| DPA1 | B | 16 | 2 | 0 | 68/104 | 0 | 0/0 | 0.83 |  |
| DPA1 | C4 | 16 | 1 | 0 | 73/104 | 0 | 1/0 | 0.70 | 1/0/0 |
| DPA1 | C4R4n | 16 | 1 | 0 | 73/104 | 0 | 1/0 | 0.70 | 1/0/0 |
| DPA1 | C4R1n | 16 | 1 | 0 | 73/104 | 0 | 1/0 | 0.70 | 1/0/0 |
| DPA1 | D4 | 16 | 3 | 0 | 70/104 | 0 | 1/1 | 0.69 | 2/0/0 |
| DPA1 | D4R4n | 16 | 3 | 0 | 70/104 | 0 | 1/1 | 0.69 | 2/0/0 |
| DPA1 | D4R1n | 16 | 3 | 0 | 70/104 | 0 | 1/1 | 0.69 | 2/0/0 |
| DPB1 | A | 15 | 0 | 0 | 97/104 | 0 | 0/0 | 0.24 |  |
| DPB1 | B | 15 | 0 | 0 | 95/104 | 0 | 0/0 | 0.48 |  |
| DPB1 | C4 | 15 | 3 | 0 | 97/104 | 0 | 0/1 | 0.22 | 1/0/0 |
| DPB1 | C4R4n | 15 | 3 | 0 | 97/104 | 0 | 0/1 | 0.22 | 1/0/0 |
| DPB1 | C4R1n | 15 | 3 | 0 | 97/104 | 0 | 0/1 | 0.22 | 1/0/0 |
| DPB1 | D4 | 15 | 3 | 0 | 95/104 | 0 | 0/1 | 0.47 | 1/0/0 |
| DPB1 | D4R4n | 15 | 3 | 0 | 95/104 | 0 | 0/1 | 0.47 | 1/0/0 |
| DPB1 | D4R1n | 15 | 3 | 0 | 95/104 | 0 | 0/1 | 0.47 | 1/0/0 |
| DQA1 | A | 32 | 0 | 0 | 76/90 | 0 | 0/0 | 0.66 |  |
| DQA1 | B | 32 | 20 | 3 | 75/90 | 0 | 0/0 | 0.43 |  |
| DQA1 | C4 | 32 | 4 | 2 | 83/90 | 0 | 3/0 | 0.52 | 3/0/0 |
| DQA1 | C4R4n | 32 | 4 | 2 | 83/90 | 0 | 3/0 | 0.52 | 3/0/0 |
| DQA1 | C4R1n | 32 | 4 | 2 | 83/90 | 0 | 3/0 | 0.52 | 3/0/0 |
| DQA1 | D4 | 32 | 24 | 4 | 83/90 | 0 | 3/0 | 0.32 | 3/0/0 |
| DQA1 | D4R4n | 32 | 24 | 4 | 83/90 | 0 | 3/0 | 0.32 | 3/0/0 |
| DQA1 | D4R1n | 32 | 24 | 4 | 83/90 | 0 | 3/0 | 0.32 | 3/0/0 |
| DQB1 | A | 3 | 0 | 0 | 103/118 | 0 | 0/0 | 7.69 |  |
| DQB1 | B | 3 | 0 | 0 | 103/118 | 0 | 0/0 | 7.69 |  |
| DQB1 | C4 | 3 | 0 | 0 | 104/118 | 0 | 0/0 | 7.63 | 0/0/0 |
| DQB1 | C4R4n | 3 | 0 | 0 | 104/118 | 0 | 0/0 | 7.63 | 0/0/0 |
| DQB1 | C4R1n | 3 | 0 | 0 | 104/118 | 0 | 0/0 | 7.63 | 0/0/0 |
| DQB1 | D4 | 3 | 0 | 0 | 104/118 | 0 | 0/0 | 7.63 | 0/0/0 |
| DQB1 | D4R4n | 3 | 0 | 0 | 104/118 | 0 | 0/0 | 7.63 | 0/0/0 |
| DQB1 | D4R1n | 3 | 0 | 0 | 104/118 | 0 | 0/0 | 7.63 | 0/0/0 |
| DRB1 | A | 76 | 0 | 0 | 36/45 | 0 | 0/0 | 31.11 |  |
| DRB1 | B | 76 | 17 | 16 | 35/45 | 0 | 0/0 | 30.42 |  |
| DRB1 | C4 | 76 | 16 | 10 | 40/45 | 1 | 14/14 | 30.26 | 29/0/0 |
| DRB1 | C4R4n | 76 | 13 | 7 | 40/45 | 1 | 11/7 | 30.41 | 19/0/0 |
| DRB1 | C4R1n | 76 | 15 | 9 | 40/45 | 1 | 13/11 | 30.32 | 25/0/0 |
| DRB1 | D4 | 76 | 27 | 23 | 39/45 | 1 | 8/7 | 29.71 | 16/0/0 |
| DRB1 | D4R4n | 76 | 25 | 21 | 39/45 | 1 | 6/2 | 29.84 | 9/0/0 |
| DRB1 | D4R1n | 76 | 26 | 22 | 39/45 | 1 | 7/4 | 29.80 | 12/0/0 |
| all | A | 162 | 0 | 0 | 707/804 | 0 | 0/0 | 5.31 (2 missing) |  |
| all | B | 162 | 39 | 19 | 694/804 | 0 | 0/0 | 5.22 (2 missing) |  |
| all | C4 | 162 | 25 | 12 | 723/804 | 1 | 18/16 | 5.15 (2 missing) | 34/0/0 |
| all | C4R4n | 162 | 22 | 9 | 723/804 | 1 | 15/8 | 5.18 (2 missing) | 24/0/0 |
| all | C4R1n | 162 | 24 | 11 | 723/804 | 1 | 17/13 | 5.16 (2 missing) | 30/0/0 |
| all | D4 | 162 | 58 | 27 | 711/804 | 1 | 12/10 | 5.08 (2 missing) | 22/0/0 |
| all | D4R4n | 162 | 56 | 25 | 711/804 | 1 | 10/4 | 5.11 (2 missing) | 15/0/0 |
| all | D4R1n | 162 | 57 | 26 | 711/804 | 1 | 11/6 | 5.10 (2 missing) | 18/0/0 |
| non-DRB1 | A | 86 | 0 | 0 | 671/759 | 0 | 0/0 | 1.61 (2 missing) |  |
| non-DRB1 | B | 86 | 22 | 3 | 659/759 | 0 | 0/0 | 1.61 (2 missing) |  |
| non-DRB1 | C4 | 86 | 9 | 2 | 683/759 | 0 | 4/2 | 1.55 (2 missing) | 5/0/0 |
| non-DRB1 | C4R4n | 86 | 9 | 2 | 683/759 | 0 | 4/1 | 1.56 (2 missing) | 5/0/0 |
| non-DRB1 | C4R1n | 86 | 9 | 2 | 683/759 | 0 | 4/2 | 1.55 (2 missing) | 5/0/0 |
| non-DRB1 | D4 | 86 | 31 | 4 | 672/759 | 0 | 4/3 | 1.55 (2 missing) | 6/0/0 |
| non-DRB1 | D4R4n | 86 | 31 | 4 | 672/759 | 0 | 4/2 | 1.56 (2 missing) | 6/0/0 |
| non-DRB1 | D4R1n | 86 | 31 | 4 | 672/759 | 0 | 4/2 | 1.56 (2 missing) | 6/0/0 |

**Secondary: mean wrong bases per copy, arm minus base, paired, family-clustered (T-sim)**

- C4 minus A, all genes, all copies: -156 of 973; -0.16 [-0.21, -0.11], p=5e-05, 50 clusters bases per copy
- C4 minus A, all genes, IPD-absent copies: -99 of 161; -0.61 [-0.85, -0.35], p=5e-05, 49 clusters bases per copy
- C4 minus A, non-DRB1, all copies: -52 of 851; -0.06 [-0.10, -0.04], p=5e-05, 50 clusters bases per copy
- C4 minus A, non-DRB1, IPD-absent copies: -11 of 85; -0.13 [-0.22, -0.06], p=0.00391, 41 clusters bases per copy
- C4R4n minus A, all genes, all copies: -129 of 973; -0.13 [-0.17, -0.10], p=5e-05, 50 clusters bases per copy
- C4R4n minus A, all genes, IPD-absent copies: -81 of 161; -0.50 [-0.72, -0.30], p=5e-05, 49 clusters bases per copy
- C4R4n minus A, non-DRB1, all copies: -43 of 851; -0.05 [-0.07, -0.04], p=5e-05, 50 clusters bases per copy
- C4R4n minus A, non-DRB1, IPD-absent copies: -11 of 85; -0.13 [-0.22, -0.06], p=0.00391, 41 clusters bases per copy
- C4R1n minus A, all genes, all copies: -149 of 973; -0.15 [-0.20, -0.11], p=5e-05, 50 clusters bases per copy
- C4R1n minus A, all genes, IPD-absent copies: -92 of 161; -0.57 [-0.81, -0.32], p=5e-05, 49 clusters bases per copy
- C4R1n minus A, non-DRB1, all copies: -52 of 851; -0.06 [-0.10, -0.04], p=5e-05, 50 clusters bases per copy
- C4R1n minus A, non-DRB1, IPD-absent copies: -11 of 85; -0.13 [-0.22, -0.06], p=0.00391, 41 clusters bases per copy
- D4 minus B, all genes, all copies: -136 of 973; -0.14 [-0.19, -0.10], p=5e-05, 50 clusters bases per copy
- D4 minus B, all genes, IPD-absent copies: -82 of 161; -0.51 [-0.76, -0.31], p=5e-05, 49 clusters bases per copy
- D4 minus B, non-DRB1, all copies: -50 of 851; -0.06 [-0.09, -0.04], p=5e-05, 50 clusters bases per copy
- D4 minus B, non-DRB1, IPD-absent copies: -10 of 85; -0.12 [-0.24, -0.03], p=0.0312, 41 clusters bases per copy
- D4R4n minus B, all genes, all copies: -113 of 973; -0.12 [-0.16, -0.08], p=5e-05, 50 clusters bases per copy
- D4R4n minus B, all genes, IPD-absent copies: -66 of 161; -0.41 [-0.67, -0.21], p=6.1e-05, 49 clusters bases per copy
- D4R4n minus B, non-DRB1, all copies: -43 of 851; -0.05 [-0.07, -0.03], p=5e-05, 50 clusters bases per copy
- D4R4n minus B, non-DRB1, IPD-absent copies: -10 of 85; -0.12 [-0.24, -0.03], p=0.0312, 41 clusters bases per copy
- D4R1n minus B, all genes, all copies: -118 of 973; -0.12 [-0.16, -0.09], p=5e-05, 50 clusters bases per copy
- D4R1n minus B, all genes, IPD-absent copies: -71 of 161; -0.44 [-0.68, -0.24], p=6.1e-05, 49 clusters bases per copy
- D4R1n minus B, non-DRB1, all copies: -43 of 851; -0.05 [-0.07, -0.03], p=5e-05, 50 clusters bases per copy
- D4R1n minus B, non-DRB1, IPD-absent copies: -10 of 85; -0.12 [-0.24, -0.03], p=0.0312, 41 clusters bases per copy

**Secondary: byte-exact recovery (final call) of IPD-absent copies, all genes (T-sim)**

- C4 minus A: +12 of 162; +7.4 pp [+3.6, +11.4], p=0.00195, 49 clusters
- C4R4n minus A: +9 of 162; +5.6 pp [+2.8, +8.6], p=0.00391, 49 clusters
- C4R1n minus A: +11 of 162; +6.8 pp [+3.1, +10.8], p=0.00391, 49 clusters
- D4 minus B: +8 of 162; +4.9 pp [+1.8, +8.1], p=0.0312, 49 clusters
- D4R4n minus B: +6 of 162; +3.7 pp [+0.7, +6.7], p=0.125, 49 clusters
- D4R1n minus B: +7 of 162; +4.3 pp [+1.3, +7.4], p=0.0625, 49 clusters

## T-real

**a. novel final calls on IPD-present copies (all genes):** {'C4': 0, 'C4R4n': 0, 'C4R1n': 0, 'D4': 0, 'D4R4n': 0, 'D4R1n': 0} -> PASS (primary arms <= 1)

**b. four-field named accuracy:** PASS

- C4R4n vs A: total 861 -> 873 of 914 (+12 of 914; +1.3 pp [+0.7, +2.0], p=0.000488, 66 clusters); per-gene change {'A': 1, 'B': 1, 'C': 0, 'DPA1': 1, 'DPB1': 1, 'DQA1': 1, 'DQB1': 0, 'DRB1': 7}
- D4R4n vs B: total 844 -> 858 of 914 (+14 of 914; +1.5 pp [+0.9, +2.3], p=0.000122, 66 clusters); per-gene change {'A': 1, 'B': 1, 'C': 0, 'DPA1': 1, 'DPB1': 1, 'DQA1': 1, 'DQB1': 0, 'DRB1': 9}

**c. gene-body recovery of IPD-absent copies, C4R4n minus A, all genes:** 0 -> 12 of 134; +12 of 134; +9.0 pp [+4.8, +13.4], p=0.000488, 60 clusters (reported)

**c. gene-body recovery of IPD-absent copies, C4R4n minus A, non-DRB1 genes (secondary):** 0 -> 7 of 49; +7 of 49; +14.3 pp [+5.7, +25.0], p=0.0156, 35 clusters

**d. gene-body recovery of IPD-absent copies, D4R4n minus B, all genes:** 31 -> 38 of 134; +7 of 134; +5.2 pp [+2.1, +9.1], p=0.0156, 60 clusters

**d. gene-body recovery of IPD-absent copies, D4R4n minus B, non-DRB1 genes (secondary):** 16 -> 21 of 49; +5 of 49; +10.2 pp [+2.2, +19.6], p=0.0625, 35 clusters

**e. novel final calls on IPD-absent copies: right gene body / not**
- C4: 10/42 = 24% right; C4R4n: 8/22 = 36% right; right calls kept 8/10 -> PASS
- D4: 3/25 = 12% right; D4R4n: 2/12 = 17% right; right calls kept 2/3 -> FAIL

**Secondary: per gene (T-real)**

| gene | arm | IPD-absent | body | byte-exact final | 4-field named | novel on present | novel absent right/wrong | mean wrong bases per copy | closer/equal/farther |
|---|---|---|---|---|---|---|---|---|---|
| A | A | 2 | 0 | 0 | 128/130 | 0 | 0/0 | 0.02 |  |
| A | B | 2 | 0 | 0 | 128/130 | 0 | 0/0 | 0.02 |  |
| A | C4 | 2 | 0 | 0 | 129/130 | 0 | 0/0 | 0.02 | 0/0/0 |
| A | C4R4n | 2 | 0 | 0 | 129/130 | 0 | 0/0 | 0.02 | 0/0/0 |
| A | C4R1n | 2 | 0 | 0 | 129/130 | 0 | 0/0 | 0.02 | 0/0/0 |
| A | D4 | 2 | 0 | 0 | 129/130 | 0 | 0/0 | 0.02 | 0/0/0 |
| A | D4R4n | 2 | 0 | 0 | 129/130 | 0 | 0/0 | 0.02 | 0/0/0 |
| A | D4R1n | 2 | 0 | 0 | 129/130 | 0 | 0/0 | 0.02 | 0/0/0 |
| B | A | 8 | 0 | 0 | 114/121 | 0 | 0/0 | 0.06 |  |
| B | B | 8 | 1 | 0 | 110/121 | 0 | 0/0 | 0.06 |  |
| B | C4 | 8 | 0 | 0 | 115/121 | 0 | 0/0 | 0.02 | 0/0/0 |
| B | C4R4n | 8 | 0 | 0 | 115/121 | 0 | 0/0 | 0.02 | 0/0/0 |
| B | C4R1n | 8 | 0 | 0 | 115/121 | 0 | 0/0 | 0.02 | 0/0/0 |
| B | D4 | 8 | 1 | 0 | 111/121 | 0 | 0/0 | 0.02 | 0/0/0 |
| B | D4R4n | 8 | 1 | 0 | 111/121 | 0 | 0/0 | 0.02 | 0/0/0 |
| B | D4R1n | 8 | 1 | 0 | 111/121 | 0 | 0/0 | 0.02 | 0/0/0 |
| C | A | 2 | 0 | 0 | 126/128 | 0 | 0/0 | 0.01 |  |
| C | B | 2 | 1 | 1 | 126/128 | 0 | 0/0 | 0.01 |  |
| C | C4 | 2 | 0 | 0 | 126/128 | 0 | 0/0 | 0.01 | 0/0/0 |
| C | C4R4n | 2 | 0 | 0 | 126/128 | 0 | 0/0 | 0.01 | 0/0/0 |
| C | C4R1n | 2 | 0 | 0 | 126/128 | 0 | 0/0 | 0.01 | 0/0/0 |
| C | D4 | 2 | 1 | 1 | 126/128 | 0 | 0/0 | 0.01 | 0/0/0 |
| C | D4R4n | 2 | 1 | 1 | 126/128 | 0 | 0/0 | 0.01 | 0/0/0 |
| C | D4R1n | 2 | 1 | 1 | 126/128 | 0 | 0/0 | 0.01 | 0/0/0 |
| DPA1 | A | 10 | 0 | 0 | 111/122 | 0 | 0/0 | 0.21 |  |
| DPA1 | B | 10 | 1 | 0 | 102/122 | 0 | 0/0 | 0.21 |  |
| DPA1 | C4 | 10 | 1 | 0 | 112/122 | 0 | 0/0 | 0.20 | 0/0/0 |
| DPA1 | C4R4n | 10 | 1 | 0 | 112/122 | 0 | 0/0 | 0.20 | 0/0/0 |
| DPA1 | C4R1n | 10 | 1 | 0 | 112/122 | 0 | 0/0 | 0.20 | 0/0/0 |
| DPA1 | D4 | 10 | 2 | 0 | 103/122 | 0 | 0/0 | 0.20 | 0/0/0 |
| DPA1 | D4R4n | 10 | 2 | 0 | 103/122 | 0 | 0/0 | 0.20 | 0/0/0 |
| DPA1 | D4R1n | 10 | 2 | 0 | 103/122 | 0 | 0/0 | 0.20 | 0/0/0 |
| DPB1 | A | 12 | 0 | 0 | 117/118 | 0 | 0/0 | 0.28 |  |
| DPB1 | B | 12 | 1 | 1 | 117/118 | 0 | 0/0 | 0.25 |  |
| DPB1 | C4 | 12 | 5 | 3 | 118/118 | 0 | 3/0 | 0.20 | 3/0/0 |
| DPB1 | C4R4n | 12 | 3 | 1 | 118/118 | 0 | 1/0 | 0.27 | 1/0/0 |
| DPB1 | C4R1n | 12 | 5 | 3 | 118/118 | 0 | 3/0 | 0.20 | 3/0/0 |
| DPB1 | D4 | 12 | 5 | 3 | 118/118 | 0 | 2/0 | 0.20 | 2/0/0 |
| DPB1 | D4R4n | 12 | 4 | 2 | 118/118 | 0 | 1/0 | 0.23 | 1/0/0 |
| DPB1 | D4R1n | 12 | 5 | 3 | 118/118 | 0 | 2/0 | 0.20 | 2/0/0 |
| DQA1 | A | 15 | 0 | 0 | 106/117 | 0 | 0/0 | 0.36 |  |
| DQA1 | B | 15 | 12 | 12 | 105/117 | 0 | 0/0 | 0.18 |  |
| DQA1 | C4 | 15 | 3 | 3 | 107/117 | 0 | 3/0 | 0.30 | 3/0/0 |
| DQA1 | C4R4n | 15 | 3 | 3 | 107/117 | 0 | 3/0 | 0.30 | 3/0/0 |
| DQA1 | C4R1n | 15 | 3 | 3 | 107/117 | 0 | 3/0 | 0.30 | 3/0/0 |
| DQA1 | D4 | 15 | 13 | 13 | 106/117 | 0 | 1/0 | 0.14 | 1/0/0 |
| DQA1 | D4R4n | 15 | 13 | 13 | 106/117 | 0 | 1/0 | 0.14 | 1/0/0 |
| DQA1 | D4R1n | 15 | 13 | 13 | 106/117 | 0 | 1/0 | 0.14 | 1/0/0 |
| DQB1 | A | 0 | 0 | 0 | 128/131 | 0 | 0/0 | 0.02 |  |
| DQB1 | B | 0 | 0 | 0 | 127/131 | 0 | 0/0 | 0.02 |  |
| DQB1 | C4 | 0 | 0 | 0 | 128/131 | 0 | 0/0 | 0.01 | 0/0/0 |
| DQB1 | C4R4n | 0 | 0 | 0 | 128/131 | 0 | 0/0 | 0.01 | 0/0/0 |
| DQB1 | C4R1n | 0 | 0 | 0 | 128/131 | 0 | 0/0 | 0.01 | 0/0/0 |
| DQB1 | D4 | 0 | 0 | 0 | 127/131 | 0 | 0/0 | 0.01 | 0/0/0 |
| DQB1 | D4R4n | 0 | 0 | 0 | 127/131 | 0 | 0/0 | 0.01 | 0/0/0 |
| DQB1 | D4R1n | 0 | 0 | 0 | 127/131 | 0 | 0/0 | 0.01 | 0/0/0 |
| DRB1 | A | 85 | 0 | 0 | 31/47 | 0 | 0/0 | 4.56 |  |
| DRB1 | B | 85 | 15 | 14 | 29/47 | 0 | 0/0 | 3.60 |  |
| DRB1 | C4 | 85 | 5 | 2 | 38/47 | 0 | 4/32 | 3.49 | 35/0/1 |
| DRB1 | C4R4n | 85 | 5 | 2 | 38/47 | 0 | 4/14 | 3.81 | 17/0/1 |
| DRB1 | C4R1n | 85 | 4 | 2 | 38/47 | 0 | 3/25 | 3.58 | 27/0/1 |
| DRB1 | D4 | 85 | 17 | 15 | 38/47 | 0 | 0/22 | 2.88 | 21/0/1 |
| DRB1 | D4R4n | 85 | 17 | 15 | 38/47 | 0 | 0/10 | 3.03 | 9/0/1 |
| DRB1 | D4R1n | 85 | 17 | 15 | 38/47 | 0 | 0/18 | 2.92 | 17/0/1 |
| all | A | 134 | 0 | 0 | 861/914 | 0 | 0/0 | 0.69 |  |
| all | B | 134 | 31 | 28 | 844/914 | 0 | 0/0 | 0.54 |  |
| all | C4 | 134 | 14 | 8 | 873/914 | 0 | 10/32 | 0.53 | 41/0/1 |
| all | C4R4n | 134 | 12 | 6 | 873/914 | 0 | 8/14 | 0.58 | 21/0/1 |
| all | C4R1n | 134 | 13 | 8 | 873/914 | 0 | 9/25 | 0.54 | 33/0/1 |
| all | D4 | 134 | 39 | 32 | 858/914 | 0 | 3/22 | 0.44 | 24/0/1 |
| all | D4R4n | 134 | 38 | 31 | 858/914 | 0 | 2/10 | 0.46 | 11/0/1 |
| all | D4R1n | 134 | 39 | 32 | 858/914 | 0 | 3/18 | 0.44 | 20/0/1 |
| non-DRB1 | A | 49 | 0 | 0 | 830/867 | 0 | 0/0 | 0.14 |  |
| non-DRB1 | B | 49 | 16 | 14 | 815/867 | 0 | 0/0 | 0.11 |  |
| non-DRB1 | C4 | 49 | 9 | 6 | 835/867 | 0 | 6/0 | 0.11 | 6/0/0 |
| non-DRB1 | C4R4n | 49 | 7 | 4 | 835/867 | 0 | 4/0 | 0.12 | 4/0/0 |
| non-DRB1 | C4R1n | 49 | 9 | 6 | 835/867 | 0 | 6/0 | 0.11 | 6/0/0 |
| non-DRB1 | D4 | 49 | 22 | 17 | 820/867 | 0 | 3/0 | 0.09 | 3/0/0 |
| non-DRB1 | D4R4n | 49 | 21 | 16 | 820/867 | 0 | 2/0 | 0.09 | 2/0/0 |
| non-DRB1 | D4R1n | 49 | 22 | 17 | 820/867 | 0 | 3/0 | 0.09 | 3/0/0 |

**Secondary: mean wrong bases per copy, arm minus base, paired, family-clustered (T-real)**

- C4 minus A, all genes, all copies: -168 of 1056; -0.16 [-0.22, -0.10], p=5e-05, 66 clusters bases per copy
- C4 minus A, all genes, IPD-absent copies: -140 of 134; -1.04 [-1.54, -0.59], p=5e-05, 60 clusters bases per copy
- C4 minus A, non-DRB1, all copies: -27 of 924; -0.03 [-0.05, -0.01], p=0.000122, 66 clusters bases per copy
- C4 minus A, non-DRB1, IPD-absent copies: -16 of 49; -0.33 [-0.62, -0.10], p=0.0156, 35 clusters bases per copy
- C4R4n minus A, all genes, all copies: -118 of 1056; -0.11 [-0.17, -0.06], p=5e-05, 66 clusters bases per copy
- C4R4n minus A, all genes, IPD-absent copies: -90 of 134; -0.67 [-1.11, -0.28], p=0.00145, 60 clusters bases per copy
- C4R4n minus A, non-DRB1, all copies: -19 of 924; -0.02 [-0.03, -0.01], p=0.000488, 66 clusters bases per copy
- C4R4n minus A, non-DRB1, IPD-absent copies: -8 of 49; -0.16 [-0.38, -0.02], p=0.0625, 35 clusters bases per copy
- C4R1n minus A, all genes, all copies: -156 of 1056; -0.15 [-0.21, -0.09], p=5e-05, 66 clusters bases per copy
- C4R1n minus A, all genes, IPD-absent copies: -128 of 134; -0.96 [-1.43, -0.51], p=5e-05, 60 clusters bases per copy
- C4R1n minus A, non-DRB1, all copies: -27 of 924; -0.03 [-0.05, -0.01], p=0.000122, 66 clusters bases per copy
- C4R1n minus A, non-DRB1, IPD-absent copies: -16 of 49; -0.33 [-0.62, -0.10], p=0.0156, 35 clusters bases per copy
- D4 minus B, all genes, all copies: -114 of 1056; -0.11 [-0.16, -0.06], p=5e-05, 66 clusters bases per copy
- D4 minus B, all genes, IPD-absent copies: -83 of 134; -0.62 [-1.06, -0.25], p=0.001, 60 clusters bases per copy
- D4 minus B, non-DRB1, all copies: -19 of 924; -0.02 [-0.04, -0.01], p=0.000977, 66 clusters bases per copy
- D4 minus B, non-DRB1, IPD-absent copies: -9 of 49; -0.18 [-0.43, +0.00], p=0.25, 35 clusters bases per copy
- D4R4n minus B, all genes, all copies: -90 of 1056; -0.09 [-0.14, -0.04], p=0.00025, 66 clusters bases per copy
- D4R4n minus B, all genes, IPD-absent copies: -59 of 134; -0.44 [-0.88, -0.08], p=0.0283, 60 clusters bases per copy
- D4R4n minus B, non-DRB1, all copies: -15 of 924; -0.02 [-0.03, -0.01], p=0.00195, 66 clusters bases per copy
- D4R4n minus B, non-DRB1, IPD-absent copies: -5 of 49; -0.10 [-0.31, +0.00], p=0.5, 35 clusters bases per copy
- D4R1n minus B, all genes, all copies: -109 of 1056; -0.10 [-0.16, -0.05], p=5e-05, 66 clusters bases per copy
- D4R1n minus B, all genes, IPD-absent copies: -78 of 134; -0.58 [-1.02, -0.21], p=0.00275, 60 clusters bases per copy
- D4R1n minus B, non-DRB1, all copies: -19 of 924; -0.02 [-0.04, -0.01], p=0.000977, 66 clusters bases per copy
- D4R1n minus B, non-DRB1, IPD-absent copies: -9 of 49; -0.18 [-0.43, +0.00], p=0.25, 35 clusters bases per copy

**Secondary: byte-exact recovery (final call) of IPD-absent copies, all genes (T-real)**

- C4 minus A: +8 of 134; +6.0 pp [+2.3, +9.7], p=0.00781, 60 clusters
- C4R4n minus A: +6 of 134; +4.5 pp [+1.5, +7.9], p=0.0312, 60 clusters
- C4R1n minus A: +8 of 134; +6.0 pp [+2.3, +9.7], p=0.00781, 60 clusters
- D4 minus B: +4 of 134; +3.0 pp [+0.7, +6.0], p=0.125, 60 clusters
- D4R4n minus B: +3 of 134; +2.2 pp [+0.0, +5.0], p=0.25, 60 clusters
- D4R1n minus B: +4 of 134; +3.0 pp [+0.7, +6.0], p=0.125, 60 clusters
