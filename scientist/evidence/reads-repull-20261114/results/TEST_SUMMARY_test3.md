# Endpoints, tag test3

**f. fallbacks:** 4 of 136 tasks (FAIL; limit 2). Gene-level non-complete by kind: {'None': 3}. Detail: {'T-real3/HG04184/C_attempt1': ['task_failed'], 'T-real3/HG04199/D_attempt1': ['task_failed'], 'T-real3/NA18971/C_attempt1': ['task_failed'], 'T-real3/NA18971/D_attempt1': ['task_failed']}


## T-real3

**a. novel final calls on IPD-present copies (all genes):** {'C4': 1, 'C4R4n': 1, 'C4R1n': 0, 'D4': 0, 'D4R4n': 0, 'D4R1n': 0} -> PASS (primary arms <= 1)

**b. four-field named accuracy:** PASS

- C4R4n vs A: total 868 -> 879 of 914 (+11 of 914; +1.2 pp [+0.6, +1.9], p=0.000977, 66 clusters); per-gene change {'A': 1, 'B': 1, 'C': 0, 'DPA1': 1, 'DPB1': 0, 'DQA1': 1, 'DQB1': 0, 'DRB1': 7}
- D4R4n vs B: total 852 -> 865 of 914 (+13 of 914; +1.4 pp [+0.8, +2.2], p=0.000244, 66 clusters); per-gene change {'A': 1, 'B': 1, 'C': 0, 'DPA1': 1, 'DPB1': 0, 'DQA1': 1, 'DQB1': 0, 'DRB1': 9}

**c. gene-body recovery of IPD-absent copies, C4R4n minus A, all genes:** 0 -> 18 of 134; +18 of 134; +13.4 pp [+8.1, +19.4], p=3.05e-05, 60 clusters (reported)

**c. gene-body recovery of IPD-absent copies, C4R4n minus A, non-DRB1 genes (secondary):** 0 -> 7 of 49; +7 of 49; +14.3 pp [+5.7, +25.0], p=0.0156, 35 clusters

**d. gene-body recovery of IPD-absent copies, D4R4n minus B, all genes:** 31 -> 45 of 134; +14 of 134; +10.4 pp [+5.1, +16.3], p=0.000488, 60 clusters

**d. gene-body recovery of IPD-absent copies, D4R4n minus B, non-DRB1 genes (secondary):** 16 -> 21 of 49; +5 of 49; +10.2 pp [+2.2, +19.6], p=0.0625, 35 clusters

**e. novel final calls on IPD-absent copies: right gene body / not**
- C4: 18/56 = 32% right; C4R4n: 14/28 = 50% right; right calls kept 14/18 -> PASS
- D4: 7/29 = 24% right; D4R4n: 2/11 = 18% right; right calls kept 2/7 -> FAIL

**Secondary: per gene (T-real3)**

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
| DPA1 | A | 10 | 0 | 0 | 118/122 | 0 | 0/0 | 0.19 |  |
| DPA1 | B | 10 | 1 | 0 | 109/122 | 0 | 0/0 | 0.20 |  |
| DPA1 | C4 | 10 | 1 | 0 | 119/122 | 0 | 0/0 | 0.18 | 0/0/0 |
| DPA1 | C4R4n | 10 | 1 | 0 | 119/122 | 0 | 0/0 | 0.18 | 0/0/0 |
| DPA1 | C4R1n | 10 | 1 | 0 | 119/122 | 0 | 0/0 | 0.18 | 0/0/0 |
| DPA1 | D4 | 10 | 2 | 0 | 110/122 | 0 | 0/0 | 0.19 | 0/0/0 |
| DPA1 | D4R4n | 10 | 2 | 0 | 110/122 | 0 | 0/0 | 0.19 | 0/0/0 |
| DPA1 | D4R1n | 10 | 2 | 0 | 110/122 | 0 | 0/0 | 0.19 | 0/0/0 |
| DPB1 | A | 12 | 0 | 0 | 117/118 | 0 | 0/0 | 0.28 |  |
| DPB1 | B | 12 | 1 | 1 | 118/118 | 0 | 0/0 | 0.24 |  |
| DPB1 | C4 | 12 | 5 | 3 | 117/118 | 0 | 3/0 | 0.21 | 3/0/0 |
| DPB1 | C4R4n | 12 | 3 | 1 | 117/118 | 0 | 1/0 | 0.27 | 1/0/0 |
| DPB1 | C4R1n | 12 | 5 | 3 | 117/118 | 0 | 3/0 | 0.21 | 3/0/0 |
| DPB1 | D4 | 12 | 5 | 3 | 118/118 | 0 | 2/0 | 0.20 | 2/0/0 |
| DPB1 | D4R4n | 12 | 4 | 2 | 118/118 | 0 | 1/0 | 0.23 | 1/0/0 |
| DPB1 | D4R1n | 12 | 5 | 3 | 118/118 | 0 | 2/0 | 0.20 | 2/0/0 |
| DQA1 | A | 15 | 0 | 0 | 106/117 | 0 | 0/0 | 0.32 |  |
| DQA1 | B | 15 | 12 | 12 | 105/117 | 0 | 0/0 | 0.18 |  |
| DQA1 | C4 | 15 | 3 | 3 | 107/117 | 0 | 3/0 | 0.27 | 3/0/0 |
| DQA1 | C4R4n | 15 | 3 | 3 | 107/117 | 0 | 3/0 | 0.27 | 3/0/0 |
| DQA1 | C4R1n | 15 | 3 | 3 | 107/117 | 0 | 3/0 | 0.27 | 3/0/0 |
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
| DRB1 | A | 85 | 0 | 0 | 31/47 | 0 | 0/0 | 4.33 |  |
| DRB1 | B | 85 | 15 | 14 | 29/47 | 0 | 0/0 | 3.36 |  |
| DRB1 | C4 | 85 | 13 | 4 | 38/47 | 1 | 12/38 | 3.13 | 50/0/1 |
| DRB1 | C4R4n | 85 | 11 | 2 | 38/47 | 1 | 10/14 | 3.58 | 24/0/1 |
| DRB1 | C4R1n | 85 | 6 | 4 | 38/47 | 0 | 5/25 | 3.39 | 29/0/1 |
| DRB1 | D4 | 85 | 28 | 26 | 38/47 | 0 | 4/22 | 2.47 | 25/0/1 |
| DRB1 | D4R4n | 85 | 24 | 22 | 38/47 | 0 | 0/9 | 2.70 | 8/0/1 |
| DRB1 | D4R1n | 85 | 28 | 26 | 38/47 | 0 | 4/17 | 2.53 | 20/0/1 |
| all | A | 134 | 0 | 0 | 868/914 | 0 | 0/0 | 0.65 |  |
| all | B | 134 | 31 | 28 | 852/914 | 0 | 0/0 | 0.51 |  |
| all | C4 | 134 | 22 | 10 | 879/914 | 1 | 18/38 | 0.48 | 56/0/1 |
| all | C4R4n | 134 | 18 | 6 | 879/914 | 1 | 14/14 | 0.54 | 28/0/1 |
| all | C4R1n | 134 | 15 | 10 | 879/914 | 0 | 11/25 | 0.51 | 35/0/1 |
| all | D4 | 134 | 50 | 43 | 865/914 | 0 | 7/22 | 0.38 | 28/0/1 |
| all | D4R4n | 134 | 45 | 38 | 865/914 | 0 | 2/9 | 0.41 | 10/0/1 |
| all | D4R1n | 134 | 50 | 43 | 865/914 | 0 | 7/17 | 0.39 | 23/0/1 |
| non-DRB1 | A | 49 | 0 | 0 | 837/867 | 0 | 0/0 | 0.13 |  |
| non-DRB1 | B | 49 | 16 | 14 | 823/867 | 0 | 0/0 | 0.10 |  |
| non-DRB1 | C4 | 49 | 9 | 6 | 841/867 | 0 | 6/0 | 0.10 | 6/0/0 |
| non-DRB1 | C4R4n | 49 | 7 | 4 | 841/867 | 0 | 4/0 | 0.11 | 4/0/0 |
| non-DRB1 | C4R1n | 49 | 9 | 6 | 841/867 | 0 | 6/0 | 0.10 | 6/0/0 |
| non-DRB1 | D4 | 49 | 22 | 17 | 827/867 | 0 | 3/0 | 0.08 | 3/0/0 |
| non-DRB1 | D4R4n | 49 | 21 | 16 | 827/867 | 0 | 2/0 | 0.09 | 2/0/0 |
| non-DRB1 | D4R1n | 49 | 22 | 17 | 827/867 | 0 | 3/0 | 0.08 | 3/0/0 |

**Secondary: mean wrong bases per copy, arm minus base, paired, family-clustered (T-real3)**

- C4 minus A, all genes, all copies: -182 of 1056; -0.17 [-0.23, -0.12], p=5e-05, 66 clusters bases per copy
- C4 minus A, all genes, IPD-absent copies: -155 of 134; -1.16 [-1.62, -0.75], p=5e-05, 60 clusters bases per copy
- C4 minus A, non-DRB1, all copies: -24 of 924; -0.03 [-0.04, -0.01], p=0.000488, 66 clusters bases per copy
- C4 minus A, non-DRB1, IPD-absent copies: -15 of 49; -0.31 [-0.60, -0.07], p=0.0312, 35 clusters bases per copy
- C4R4n minus A, all genes, all copies: -115 of 1056; -0.11 [-0.16, -0.06], p=5e-05, 66 clusters bases per copy
- C4R4n minus A, all genes, IPD-absent copies: -88 of 134; -0.66 [-1.07, -0.31], p=0.0003, 60 clusters bases per copy
- C4R4n minus A, non-DRB1, all copies: -16 of 924; -0.02 [-0.03, -0.01], p=0.00195, 66 clusters bases per copy
- C4R4n minus A, non-DRB1, IPD-absent copies: -7 of 49; -0.14 [-0.35, -0.02], p=0.125, 35 clusters bases per copy
- C4R1n minus A, all genes, all copies: -148 of 1056; -0.14 [-0.20, -0.09], p=5e-05, 66 clusters bases per copy
- C4R1n minus A, all genes, IPD-absent copies: -122 of 134; -0.91 [-1.35, -0.53], p=5e-05, 60 clusters bases per copy
- C4R1n minus A, non-DRB1, all copies: -24 of 924; -0.03 [-0.04, -0.01], p=0.000488, 66 clusters bases per copy
- C4R1n minus A, non-DRB1, IPD-absent copies: -15 of 49; -0.31 [-0.60, -0.07], p=0.0312, 35 clusters bases per copy
- D4 minus B, all genes, all copies: -136 of 1056; -0.13 [-0.18, -0.08], p=5e-05, 66 clusters bases per copy
- D4 minus B, all genes, IPD-absent copies: -105 of 134; -0.78 [-1.22, -0.42], p=5e-05, 60 clusters bases per copy
- D4 minus B, non-DRB1, all copies: -18 of 924; -0.02 [-0.03, -0.01], p=0.00195, 66 clusters bases per copy
- D4 minus B, non-DRB1, IPD-absent copies: -9 of 49; -0.18 [-0.43, +0.00], p=0.25, 35 clusters bases per copy
- D4R4n minus B, all genes, all copies: -102 of 1056; -0.10 [-0.15, -0.05], p=5e-05, 66 clusters bases per copy
- D4R4n minus B, all genes, IPD-absent copies: -71 of 134; -0.53 [-0.95, -0.18], p=0.00315, 60 clusters bases per copy
- D4R4n minus B, non-DRB1, all copies: -14 of 924; -0.02 [-0.03, -0.01], p=0.00391, 66 clusters bases per copy
- D4R4n minus B, non-DRB1, IPD-absent copies: -5 of 49; -0.10 [-0.31, +0.00], p=0.5, 35 clusters bases per copy
- D4R1n minus B, all genes, all copies: -128 of 1056; -0.12 [-0.18, -0.07], p=5e-05, 66 clusters bases per copy
- D4R1n minus B, all genes, IPD-absent copies: -97 of 134; -0.72 [-1.14, -0.36], p=5e-05, 60 clusters bases per copy
- D4R1n minus B, non-DRB1, all copies: -18 of 924; -0.02 [-0.03, -0.01], p=0.00195, 66 clusters bases per copy
- D4R1n minus B, non-DRB1, IPD-absent copies: -9 of 49; -0.18 [-0.43, +0.00], p=0.25, 35 clusters bases per copy

**Secondary: byte-exact recovery (final call) of IPD-absent copies, all genes (T-real3)**

- C4 minus A: +10 of 134; +7.5 pp [+3.7, +11.5], p=0.00195, 60 clusters
- C4R4n minus A: +6 of 134; +4.5 pp [+1.5, +7.9], p=0.0312, 60 clusters
- C4R1n minus A: +10 of 134; +7.5 pp [+3.7, +11.5], p=0.00195, 60 clusters
- D4 minus B: +15 of 134; +11.2 pp [+5.4, +17.4], p=0.000488, 60 clusters
- D4R4n minus B: +10 of 134; +7.5 pp [+3.0, +12.5], p=0.00391, 60 clusters
- D4R1n minus B: +15 of 134; +11.2 pp [+5.4, +17.4], p=0.000488, 60 clusters
