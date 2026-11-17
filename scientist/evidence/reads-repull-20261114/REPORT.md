# Step 2: why real read sets have a hole in DRB1 intron 1, and a re-pull rule that closes it

Stage `reads-repull-20261114` (IBEX `/encrypted/e3008/hohndor-dogohla-20261020/reads-repull-20261114`).
Development donors only (20 HPRC donors). T-real was not touched. Written by the supervisor from the
implementing agent's hand-back (the agent could not write this file); headline counts recounted by the
supervisor from `results/score_rows_repull_rule*.tsv`.

Naming: on IBEX the adopted read set is `D-real3` (rule 2). `D-real2` is the rejected rule 1.

## Outcome

- **Cause.** The donor's own intron-1 reads exist in the original CRAM but were never pulled. The NYGC
  alignment (bwa 0.7.15 -Y, GATK 3.5) placed them on other chromosomes, mostly as proper pairs with
  MAPQ 20-59. They are not unmapped and not on decoy/alt/HLA contigs, so no contig list catches them.
  About 2 kb of intron 1 in the DR53-group copies (*09, *07, one *12-like copy) is close to a high-copy
  interspersed repeat. T1K and the builder lose nothing on top of this.
- **Rule 1, rejected: give T1K all whole-genome reads** (5.2-6.7 M pairs per donor). Closes the hole but
  wrecks DRB1 typing: repeat reads pile onto DR53 alleles (DRB1 four-field 8 -> 0 of 15).
- **Rule 2, adopted: old regional pull plus whole-genome pairs carrying MHC-specific IPD 31-mers.** Adds
  86-594 pairs per donor, closes the hole, raises C4 gene-body recovery from 4 to 9 of 25, and improves
  standard typing (arm A all genes 253 -> 258 of 276) with no gene losing a copy.

## Part A: diagnosis

**A1. How D-real was made.** `samtools view --fetch-pairs --region-file recruitment.bed -F 0x900`, collate,
fastq. Regions: chr6:28,500,000-33,500,000 plus every `chr6_*` and `HLA-*` contig (542 bed lines). No MAPQ
or truth filter. Mates outside the region are kept; fully unmapped pairs, decoys, chrUn/random and other
chromosomes are not pulled. 17 donors use NYGC CRAMs; NA18944/45/70 use this project's own bwa realignment.

**A2. Own-copy coverage per 100-base window** (`src/cov.py`, job 55582799; `results/cov_D-real_vs_D-sim.md`).
A read counts for a true copy if it fits it over its full length with <= 1 edit and fits the donor's other
true copy worse. 10 donors have a hole (real <= 1 where simulated >= 5) at alignment columns ~4,275-5,800:
*09 HG00558, HG00609, HG00621, HG02040, HG02523; *07 HG02080, HG02602, HG03710, HG03874; *12-like HG03784.
Median coverage inside the hole 1-3 (real) vs 18-27 (simulated). DRB1*04 copies show no hole; neither do
the self-realigned Japanese donors. T1K's candidate reads give the same numbers, so reads are missing
before T1K.

**A3. Searching the whole CRAM** (jobs 55584289, 55584752). 31-mers from the hole stretch of each true copy,
absent from the other copy, minus k-mers seen on chr1:1-60 Mb. Reads in the hole fitting the copy with <= 2 edits:

| donor | reads in hole | in old pull | where the rest sit |
|---|---:|---:|---|
| HG00558 | 9,586 | 14 | all chromosomes; MAPQ 20-59: 6,991 |
| HG00621 | 9,054 | 10 | all chromosomes |
| HG02080 | 22,936 | 43 | all chromosomes |
| HG03784 | 29 | 9 | scattered |
| HG00438 (control) | 44 + 20 | 11 + 7 | scattered |

These counts are 30-70 times expected depth (hence: a repeat). None of HG00558's 9,586 reads fits the copy
better than the pass-1 IPD allele. Diagnostic k-mers (copy vs pass-1 alleles; job 55584921): HG00558 5 reads
in the old pull, 8 elsewhere; HG00621 1 and 14; they sit at columns 5,500-6,000, the site the builder patches
on simulated reads. HG02080 and HG03784 have no diagnostic k-mer (FAILURE.json by design).

**A4. Correction to the earlier diagnosis.** The blind-region guard blocked 0 patches on old D-real; patches
failed earlier for lack of own-version reads.

## Part B: fix and measure

**B5. Rule 2** (`src/bait.py`, `src/select3.py`; bait job 55594839). No truth used; thresholds fixed beforehand.
1. Keep the old regional pull unchanged.
2. Add every whole-genome pair in which either mate has >= 3 bait 31-mers among the 31-mers starting at
   every 4th position (canonical).
3. Bait: the 931,070 canonical 31-mers of all 29,429 IPD genomic alleles minus the 43,021 that occur in
   GRCh38 outside the MHC (chr6:28.5-33.5 Mb masked, `*_alt`/`HLA-*` left out; decoys, unplaced, EBV
   included; jellyfish). 888,049 remain.
4. Whole-genome pairs are scanned through T1K's existing full-WGS candidate set as a prefilter.

**B6. Read pairs.** Old 625,859-1,047,319 per donor; rule 2 adds 86-594 (7,384 total); rule 1 5.20-6.68 M.

**B7. Scores** (`src/score2.py` = frozen scoring logic; `src/t1k_a.py`; `src/run_b.py` imports `drb1_run.run`).

| 20 dev donors | old D-real | rule 2 | rule 1 | D-sim |
|---|---|---|---|---|
| IPD-absent DRB1 copies | 25 | 25 | 25 | 25 |
| C4 gene body recovered | 4 | 9 (C4 - A +36 pp [19.2, 56.5]) | 0 | 15 |
| C4 byte-exact (any slot / final) | 3 / 2 | 3 / 2 | 0 | 7 / 6 |
| DRB1 four-field named, A -> C4 (of 15) | 8 -> 12 | 8 -> 12 | 0 -> 0 | 10 -> 13 |
| Novel final calls on IPD-present DRB1 copies | 0 | 0 | 1 | 0 |
| All genes four-field, arm A (of 276) | 253 | 258 | 250 | 244 |
| All genes four-field, C4 (of 276) | 257 | 262 | 250 | 247 |

Arm A per gene old -> rule 2: DPA1 28 -> 32, DPB1 35 -> 36, others unchanged. Of C4's novel calls on absent
copies with rule 2: 6 right gene body, 4 not.

**B8. Hole and guard.** Every hole donor now has own coverage 12-27 in columns 3,300-5,600 (simulated 18-27).
Blind regions 32 -> 0 (rule 1: 120); patches kept 15 -> 23 (D-sim 34).

## Deviations and open points

- Own-copy rule: "fits the other true copy worse" (the brief's stricter version removes nearly all reads
  where the copy equals the IPD allele locally).
- Japanese donors: first rule-2 selection duplicated reads (`/1` suffix mismatch); fixed and rerun.
- The T1K-candidate prefilter could miss baited reads; not measured.
- Repeat family not identified. Rule 1 coverage: 5 tasks timed out, not rerun (rule rejected anyway).
- Rule 2 gains gene bodies but no extra byte-exact copies (3 -> 3).

## Cost of re-pulling the 66 T-real donors

40 donors have full-WGS T1K candidate files (selection 1-2 min each). 26 (`remaining26`) need the T1K
extractor streamed over the CRAM (~3 h, 8 CPUs, 16 GB each; ~8 h wall at 12 at a time). Then typing and
builder for all 66, about 20 h wall.

## Step 2b: frozen test of the re-pull on the 66 real Asian donors (T-real3), 2026-11-17

Design and freeze: `DESIGN_TEST.txt`, `FREEZE_TEST.json` (frozen 2026-11-15 13:12, before any test job).
Result tables: `results/SUMMARY_test3.txt`, `results/endpoints_test3.json`, `results/TEST_SUMMARY_test3.md`,
`results/score_rows_test3.tsv` (16,896 rows), `results/TEST_DEVIATIONS.txt`. Scoring job 55701014.
The scorer reproduced the old T-real rows exactly (8,448 of 8,448 identical).

### Outcome

| Endpoint | Result | Verdict |
|---|---|---|
| P1 (primary): whole genes rebuilt for copies missing from IPD, builder on IPD reference (C4R4n), new reads vs old | 12 -> 18 of 134 (+4.5 pp, 95% interval +1.4 to +8.2, p=0.031) | PASS |
| P2 i: T1K four-field names, IPD reference | 861 -> 868 of 914; no gene down (DPA1 +7) | PASS |
| P2 i: T1K four-field names, IPD + Asian panel | 844 -> 852 of 914; no gene down (DPA1 +7, DPB1 +1) | PASS |
| P2 ii: false new-allele calls on copies that are in IPD | C4R4n 1, D4R4n 0 (limit 1 each) | PASS |
| P2 iii: tasks in fallback (limit 2 of 264) | see "Deviations" | FAIL as printed; see below |
| P3: same as P1 with the Asian panel (D4R4n) | 38 -> 45 of 134 (+5.2 pp, +1.4 to +9.6, p=0.039) | reported |
| P4: builder vs T1K alone on new reads | IPD: 0 -> 18 of 134 (+13.4 pp, +8.1 to +19.4); panel: 31 -> 45 (+10.4 pp, +5.2 to +16.3); four-field 868 -> 879 and 852 -> 865 of 914 | reported |
| P5: new-allele calls on IPD-absent DRB1 copies, right / wrong | C4R4n 4/14 -> 10/14; D4R4n 0/10 -> 0/9; all genes C4R4n 8/14 -> 14/14 | reported |
| P6: DRB1 | all gains are in DRB1; DR53-group carriers: C4R4n 2 -> 8 of 50, D4R4n 6 -> 13 of 50; other DRB1 carriers unchanged (3 -> 3, 11 -> 11 of 35); blind regions DRB1 103 -> 0, DPA1 6 -> 0; pairs added per donor median 365 (78-782) | reported |

Arm A two- and three-field accuracy is unchanged except DPA1 +2 copies and DPB1 -1 copy (3-field 131 -> 130, 2-field 132 -> 131 of 132).

### Deviations (all logged in `results/TEST_DEVIATIONS.txt` before any score existed)

- D-T1 / D-T2: the read extraction for donor NA18565 ran out of memory at 32 GB (it has about 1.17 billion reads).
  It was rerun once with 64 GB and the same code and completed (10 h 07 min). The frozen design says a failed
  pre-processing task makes that donor's four tasks failures, so both readings are given:
  - with the rerun: P1 12 -> 18, P3 38 -> 45 (table above);
  - strict (NA18565 counted as failed): P1 12 -> 17, P3 38 -> 44. The frozen interval was not recomputed for this
    reading; a plain donor-level bootstrap gives +0.7 to +7.5 pp (P1) and +0.7 to +9.0 pp (P3), still above 0.
    Under this reading P2 iii fails (4 failed tasks, limit 2).
- Four builder tasks (HG04184 C, NA18971 C and D, HG04199 D) hit the 1800 s timeout while loading inputs on slow
  nodes and were retried once on another node, as the design allows; all four then completed in under 2 minutes.
  The endpoint script prints "4 of 264 in fallback -> FAIL" because it counts the saved first attempts
  (`*_attempt1`). No final result comes from a fallback: all 264 final tasks are complete. So P2 iii is
  0 of 264 on the design's own definition (with the rerun), and 4 of 264 on the strict reading.

### Supervisor recount (from `results/score_rows_test3.tsv`)

- body_recovered among IPD-absent copies: A 0 -> 0, B 31 -> 31, C4R4n 12 -> 18, D4R4n 38 -> 45 of 134. Confirmed.
- named4 totals: A 861 -> 868, B 844 -> 852, C4R4n 873 -> 879, D4R4n 858 -> 865. Confirmed.
- All six C4R4n gains are DRB1, in six different donors; no donor lost a copy. D4R4n: 8 donors up, 1 down.
- The one false new-allele call: HG02135 DRB1 (true allele DRB1*09:01:02:08), arm C4R4n.

### Reading

Re-pulling the reads fixes the DR53-group hole on untouched real donors: it removes every blind region and lets the
builder rebuild 6-7 more whole genes that IPD does not contain, without hurting ordinary typing. The gain is
modest in absolute terms (18 of 134 with IPD alone, 45 of 134 with the Asian panel plus builder) and is confined
to DRB1. The builder still makes as many wrong as right new DRB1 calls on the IPD reference (10 right, 14 not).
