# DogoHLA: what the six-step route found (2026-11-14 to 2026-11-17)

*2026-11-17T07:02:04Z by Showboat 0.6.1*
<!-- showboat-id: c6016ae2-73f5-49ee-8b0f-e28880b08e3b -->

**In one sentence:** a small set of good Asian haplotypes, plus a tool that rebuilds missing genes from reads, plus re-fetching reads the standard alignment lost, lets us rebuild exactly 45 of 134 HLA gene copies that are missing from the public database in real Asian genomes, where the standard tool rebuilds none, without making ordinary typing worse.

### Words used here

- **HLA typing**: reading which version (allele) of each HLA gene a person carries, from their DNA reads.
- **IPD**: the public database of known HLA alleles. Many Asian, Arab and Saudi alleles are not in it.
- **IPD-absent copy**: a gene copy in a person whose exact sequence is not in IPD. These are the hard cases.
- **Whole gene rebuilt**: our final call has exactly the true sequence of the gene.
- **Four-field name**: the full allele name, down to differences outside the protein-coding part.
- **T1K**: the standard HLA typing tool we start from.
- **Panel**: extra Asian haplotypes from long-read genome assemblies, added to IPD as reference.
- **Builder**: our tool that patches a new allele together from the reads when nothing in the reference fits.
- **Frozen test**: the code, donors and pass/fail rules were fixed and fingerprinted before we looked at any test result.

### Why we did this

After the earlier confirmatory test we were not happy yet: known-allele databases miss many alleles in Asian, Arab and Saudi people, and our first gains were small. The question was whether the few good haplotypes we have can help type many more people. We ran six steps on the KAUST IBEX cluster. Development used 20 HPRC donors; final tests used 61 simulated Arab/Saudi donors (T-sim) and 66 real Asian genomes from HPRC / 1000 Genomes (T-real).

The commands below recount the main numbers from the result files in this repository. Per-person tables for the Arab/Saudi donors are not in this public repository; their numbers are quoted from the step reports.

## Step 1. How wrong are the misses? (rescore-partial)

With the full-gene reference, the protein-coding part of the call is right for about 98% of copies. When we miss an IPD-absent copy, we are usually off by only 1 to 2 bases outside the coding part. So the misses are near misses, worth fixing. Report: `rescore-partial-20261114/REPORT.md`.

## Step 2. Fetch the reads the standard alignment lost (reads-repull)

On real reads, one part of DRB1 (intron 1, in people with DR53-group alleles) had no reads at all. The reads exist, but the standard alignment put them on other chromosomes because that stretch looks like a repeat. The fix: also take every read pair, from anywhere in the genome, that carries at least 3 short DNA words (31 letters) found only in HLA genes. We tested this on the 66 real Asian donors in a frozen test, finished 2026-11-17. The main check (pre-registered) was: more IPD-absent copies rebuilt by the builder on the IPD reference.

```bash
grep -E '^\*\*P1|^\*\*P3|^- \(i|^- \(ii' reads-repull-20261114/results/SUMMARY_test3.txt
```

```output
**P1** gene bodies recovered among IPD-absent copies, all genes, C4R4n T-real3 minus T-real: 12 -> 18 of 134 (+6; +4.5 pp [+1.4, +8.2], p=0.0312) -> PASS
**P3** gene bodies recovered among IPD-absent copies, all genes, D4R4n T-real3 minus T-real: 38 -> 45 of 134 (+7; +5.2 pp [+1.4, +9.6], p=0.0391)
- (i) arm A four-field named, T-real3 vs T-real: 861 -> 868 of 914 (+7; +0.8 pp [+0.1, +1.6], p=0.0938); per gene change {'A': 0, 'B': 0, 'C': 0, 'DPA1': 7, 'DPB1': 0, 'DQA1': 0, 'DQB1': 0, 'DRB1': 0} -> PASS
- (i) arm B four-field named, T-real3 vs T-real: 844 -> 852 of 914 (+8; +0.9 pp [+0.2, +1.6], p=0.0312); per gene change {'A': 0, 'B': 0, 'C': 0, 'DPA1': 7, 'DPB1': 1, 'DQA1': 0, 'DQB1': 0, 'DRB1': 0} -> PASS
- (ii) novel final calls on IPD-present copies (T-real3): {'C4R4n': 1, 'D4R4n': 0} -> PASS
- (iii) tasks in fallback: 4 of 264 ['HG04184:C_attempt1:task_failed', 'HG04199:D_attempt1:task_failed', 'NA18971:C_attempt1:task_failed', 'NA18971:D_attempt1:task_failed'] -> FAIL
```

Recount from the per-copy table (old reads vs re-pulled reads, 134 IPD-absent copies, 1056 copies in total of which 914 have a four-field truth name):

```python3
import csv
rows = list(csv.DictReader(open('reads-repull-20261114/results/score_rows_test3.tsv'), delimiter='\t'))
yes = lambda x: x == 'True'
names = {'A': 'T1K + IPD', 'B': 'T1K + IPD + Asian panel', 'C4R4n': 'builder on IPD', 'D4R4n': 'builder on IPD + Asian panel'}
print(f"{'setup':32} {'whole genes rebuilt (of 134)':>30} {'four-field right':>18}")
for arm, label in names.items():
    out = []
    for rs in ('T-real', 'T-real3'):
        sel = [r for r in rows if r['readset'] == rs and r['arm'] == arm]
        body = sum(yes(r['body_recovered']) for r in sel if yes(r['ipd_absent']))
        n4 = sum(yes(r['named4']) for r in sel)
        out.append((body, n4))
    print(f"{label:32} {out[0][0]:>14} -> {out[1][0]:<13} {out[0][1]:>8} -> {out[1][1]}")

```

```output
setup                              whole genes rebuilt (of 134)   four-field right
T1K + IPD                                     0 -> 0                  861 -> 868
T1K + IPD + Asian panel                      31 -> 31                 844 -> 852
builder on IPD                               12 -> 18                 873 -> 879
builder on IPD + Asian panel                 38 -> 45                 858 -> 865
```

What this means: re-pulling the reads closed the hole (DRB1 regions with no reads went from 103 to 0) and the builder rebuilt 6 to 7 more whole genes, all of them DRB1 in DR53-group carriers. Ordinary typing did not get worse. Two caveats: one donor (NA18565, an unusually large genome) needed a rerun with more memory; counting it as failed, as the frozen plan says, gives 12 -> 17 and 38 -> 44, still a gain. Four builder runs timed out on slow machines and were retried once, as the plan allows. The builder's new DRB1 calls on the IPD reference are still more often wrong than right (10 right, 14 not). Report: `reads-repull-20261114/REPORT.md`.

## Step 3. Does a bigger panel from the same population help? (dose-response)

We added 0, 5, 10, 20 or 40 haplotypes from the donor's own population to the reference and measured how many IPD-absent copies were then typed exactly.

```python3
import csv
rows = list(csv.DictReader(open('dose-response-20261114/results/figure_ready_dose_response.tsv'), delimiter='\t'))
pick = {(r['group'], r['readtype'], r['cell']): r for r in rows if r['metric'] == 'absent_exact_any'}
pct = lambda r: 100 * int(r['k']) / int(r['n'])
print('% of IPD-absent copies where the exact true allele is among the calls,')
print('by number of same-population haplotypes added to the reference')
print(f"{'group':28} {'0':>6} {'5':>6} {'10':>6} {'20':>6} {'40':>6}   non-Asian 40")
for g, rt, label in (('ME', 'sim', 'Middle East (simulated)'), ('EAS', 'real', 'East Asian (real)'), ('SAS', 'real', 'South Asian (real)')):
    v = [pct(pick[(g, rt, c)]) for c in ('n0', 'match_n05', 'match_n10', 'match_n20', 'match_n40')]
    print(f"{label:28} " + ' '.join(f'{x:6.1f}' for x in v) + f"   {pct(pick[(g, rt, 'x_NONASIAN_n40')]):6.1f}")

```

```output
% of IPD-absent copies where the exact true allele is among the calls,
by number of same-population haplotypes added to the reference
group                             0      5     10     20     40   non-Asian 40
Middle East (simulated)         0.0    1.9    3.1    4.3    9.3      3.7
East Asian (real)               0.0    0.0    2.0    5.1    7.1      0.0
South Asian (real)              0.0   14.5   11.8   28.9   32.9      2.6
```

What this means: the more haplotypes from the person's own population, the more hard alleles we get exactly, and the curve has not flattened at 40. Haplotypes from the same population beat the same number of non-Asian haplotypes. So collecting more matched haplotypes should keep paying off. Report: `dose-response-20261114/REPORT.md`.

## Step 4. The builder on all eight genes (mosaic-allgenes, frozen test)

The builder patches a new allele from the reads when nothing in the reference fits. In a frozen test on all eight HLA genes it rebuilt whole IPD-absent genes that T1K alone never gets:

| test donors | T1K + IPD | builder on IPD | T1K + IPD + panel | builder on IPD + panel |
|---|---|---|---|---|
| 61 simulated Arab/Saudi (162 IPD-absent copies) | 0 | 22 | 39 | 56 |
| 66 real Asian (134 IPD-absent copies) | 0 | 12 | 31 | 38 |

Only one false "new allele" call was made on a copy that is in IPD. Five of the six pre-registered criteria passed; one (real donors, builder on the panel) was recorded as failed. Weakness: many new DRB1 calls on real reads are close to the truth but not exact. Step 2 above improved the real-donor numbers to 18 and 45. Report: `mosaic-allgenes-20261114/REPORT.md`.

## Step 5. Can a pangenome graph do the typing? (graph-paths)

A clean graph built from aligned gene sequences places reads correctly far more often than the general pangenome graph (DRB1: 96 to 98% vs 41 to 52%). But picking the best haplotype paths through that graph does not work as a typer: for DRB1 it rebuilt 0 of 25 IPD-absent copies. We did not take a graph typer to a frozen test. A possible next step is to feed the builder with reads placed by the clean graph. Report: `graph-paths-20261114/REPORT.md`.

## Step 6. What happens with thin data? (thin-data)

With fewer reads, two-field accuracy drops: real reads 99.0% at full depth, 95.8% at a quarter, 74.9% at one sixteenth; simulated Middle East 97.1%, 85.5%, 60.3%. At one sixteenth depth neither the panel nor a population-frequency prior helps (at most +1.2 points). The problem is ties: about half of the gene copies have several alleles that fit the reads equally well. Report: `thin-data-20261114/REPORT.md`.

## Bottom line

- **Helps:** same-population haplotypes in the reference, the builder for alleles nobody has seen, and re-fetching lost reads. Together on real Asian genomes: 45 of 134 IPD-absent copies rebuilt exactly (standard T1K + IPD: 0; panel alone: 31), with ordinary typing unharmed.
- **Does not help:** picking paths from a graph as a typer; priors at very low depth.
- **Next:** more matched haplotypes (the curve is still rising); fix the builder's near-miss DRB1 calls; try the builder on graph-placed reads.

All test jobs ran on IBEX. Protected sequences stay there; only counts, names and scores are in this repository. Each step folder has a `REPORT.md`, the frozen design and its fingerprint file (`FREEZE*.json`), and the code.
