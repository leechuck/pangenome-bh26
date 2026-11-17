# Step 4: patch builder for all eight genes (stage `mosaic-allgenes-20261114`)

IBEX: `/encrypted/e3008/hohndor-dogohla-20261020/mosaic-allgenes-20261114` (`src/` identical to local `src/`,
checked by sha256). Development donors only in this part (20 donors of `mosaic-drb1-20261113/work/dev_tasks.txt`,
read sets D-sim and D-real, never pooled). No test-donor truth opened; the only test-donor runs in the development
part are the item-1 hang probes (frozen builder, unscored).

## Summary (development)

* Hang: not in the code. The frozen builder finishes T-sim APR051 arm D in 7-18 s in 6 repeats (3 hash seeds, the
  same node cn605-04-l) with the same output hash. The hung job stalled on Lustre file reads; I reproduced the same
  stall live (one task sat 16 min in `cl_sync_io_wait` reading `genomic.fa`, then finished). Fix: wall-clock limits
  per stage and per gene plus deterministic work caps, each with an explicit FAILURE record and a pass-1 fallback.
  Results unchanged (regression check).
* Regression (DRB1, abstain off, frozen confounded list): identical to the frozen stage: 120/120 haplotype records,
  60/60 OVERRIDES and 60/60 additions files; all `SCORES_dev3.md` DRB1 counts reproduced exactly (A, B, C4, D4,
  both read sets). Repeated with all three code versions (regress, regress2, regress3).
* All genes, same PARAMS: 480/480 gene runs, 60/60 tasks complete, no failures, no fallbacks. No novel final call on
  an IPD-present copy in any gene, arm or read set. Gains come from DRB1, DQA1 and DPB1; class I has almost no
  IPD-absent copies among development donors (A 0, B 0, C 1), so class I is evidence of safety only.
* Recommended abstain rule R4n (post hoc, see item 3): removes 5 of 6 novel calls that are not the true gene body
  (C4, both read sets) and 2 of 2 (D4), at the cost of 1 right call (DPB1, D-sim); never touches renames.

## Item 1: the hang

Truth-free, unscored runs (`src/hang_probe.py`, imports the frozen DRB1 modules read-only, writes only `work/hang/`,
times each step, dumps Python stacks every 120 s):

| job | what | result |
|---|---|---|
| 55582726_0 | T-sim APR051 arm D, frozen builder | 12.3 s, proposals sha 0cf27e01... |
| 55582807 / 08 / 09 | same, PYTHONHASHSEED 0 / 1 / 777 | 7.9 / 6.7 / 7.2 s, same sha |
| 55582810 | same, pinned to cn605-04-l (node of the hung task) | 18.0 s, same sha |
| 55582811 | T-sim APR051 arm C (control) | 6.7 s; proposals sha 730ab960... = frozen C3 proposals.json |
| 55582726_1 | D-real HG00438 arm D | 11.7 s |

Every builder step takes <= 1 s apart from mapping (~6 s) and loading (~2-3 s). The hung task (55562791_48, job
55563540): 3 h wall, 28 min CPU; its `D3/` folder was left empty with mtime 21:06 (start 20:10), i.e. the mapping
folder was removed after 56 min and no proposals.json was written in the next 2 h, although the builder needs ~5 s of
CPU. The job was mostly waiting. My `learn` run reproduced the same stall (task 55583287_5, D-sim HG02015): 954 s
reading `genomic.fa` on /ibex/project (Lustre); process state `cl_sync_io_wait`, `read_bytes` 32 GB for ~0.6 GB of
reads (the client re-fetched the file repeatedly); then it finished normally. The same donor took 15 s two minutes
earlier on another node. Cause: transient Lustre client stalls on input reads, not a loop in the code. Not proved
for the original job (its stack was not recorded).

Fix (`ag_run.py`, `drb1_patch2.py`):
- `signal.alarm` of 1800 s on loading+mapping and 1800 s per gene. Per-gene expiry writes `<gene>/FAILURE.json`
  (`failure_kind=timeout`, `fallback='pass-1 call kept'`, no OVERRIDES); task-level expiry writes `FAILURE.json`
  with `fallback_calls_table` (pass-1 table), which the scorer uses.
- Work caps per gene (family rows 20,000; events 20,000; alignments 20 M) -> `BoundExceeded`, `failure_kind=bound`,
  same fallback. Development maxima 5.4 k / 0.7 k / 52 k, so they never fire. A first cap of 3,000 family rows did
  fire on HLA-A/B/C (run 55583287); that run was discarded (`work/learn_cap3000_discard`) and repeated.
- Stack dump to the log every 600 s.

## Item 2: changes relative to the frozen DRB1 stage

PARAMS byte-identical. Unchanged: `drb1_lib.py`, `learn_confounded.py` and the other copied files.

| file | function | change | why |
|---|---|---|---|
| drb1_diag2.py | `S2` | this stage | write here only |
| drb1_diag2.py | new `load_shared()` | gene-independent part of `setup()`: inputs, pass-1 table, aliases, slots of all genes, one bwa mapping to the pass-1 alleles of all genes (target never depended on the gene) | map once per donor x arm |
| drb1_diag2.py | `setup(..., gene, shared)` | gene as parameter (MSA, slots, reads, mate-on-other-gene filter, leakage check); returns `rows_msa` | generalisation |
| drb1_run.py | `run(..., gene, shared, confounded_path)` | gene parameter; per-gene confounded list; timings, work counters, `pass1_msa_row` in proposals.json | generalisation, timing, evaluation |
| drb1_run.py | `unresolved_status()`, `ABSTAIN_RULES` | per-haplotype counts of guard rejections, ambiguous versions, still-supported own version, alternatives kept, blind regions, untestable events; rules R1-R4 | item 3 |
| drb1_run.py | end of `run()` | `OVERRIDES.json` as before; plus `OVERRIDES_Rk.json` (drops novel and rename overrides of unresolved haplotypes) and `OVERRIDES_Rkn.json` (drops novel final calls only) | rules evaluated without re-running |
| drb1_patch2.py | `WORK`, `WORK_CAP`, `BoundExceeded`, `_work()` | deterministic work caps | item 1 |
| ag_run.py (new) | | one task = read set x donor x arm, all 8 genes, per-gene COMPLETE/FAILURE, time limits, fallback records, no T1K pass 2 | items 1, 2, 4 |
| ag_learn_confounded.py (new) | | `learn_confounded.py` logic per gene | item 2 |
| ag_score.py (new) | | `drb1_score.py` for every gene; arms C4/D4 with rule suffix; gene-body columns; edlib NW edit distance of final call gene body vs true gene body, own min-sum copy-to-slot assignment | item 6 |
| ag_eval_patches.py (new) | | `drb1_eval_patches.py` for every gene and both panels | item 6 |
| ag_regress.py, ag_proposals_table.py, hang_probe.py, summarise.py, s4.sbatch, post.sbatch (new) | | regression, truth-free haplotype table, probe, local summary, job templates | |

Deviations: (a) C3/D3 dropped (never better than C4/D4 in the DRB1 test; ~10 min T1K per task). (b) No structural
change needed for class I (HLA-C families up to 5.4 k rows over both haplotypes, still <= 24 s per gene); paralog
read filter is generic. (c) Relearned confounded windows (`inputs/confounded_v1/`, D-sim arm C, job 55584444):
A 0, B 0, C 2, DPA1 3, DPB1 0, DQA1 0, DQB1 0, DRB1 13; DRB1 list identical to the frozen one. (d) Rule R4 added
after R1-R3 had been evaluated (post hoc).

## Item 5: regression check

`results/regression_regress{,2,3}.tsv`: generalised code, DRB1 only, frozen confounded list, abstain off; D-sim C,
D-real C, D-real D (60 tasks): 120/120 haplotype records equal (new fields ignored); 60/60 OVERRIDES.json and 60/60
additions.fa byte-identical; pass1, same_family, reads_assigned, reads_used, expected_depth all equal. Scores
recounted from `score_rows_regress.tsv` equal `SCORES_dev3.md`: D-sim C4 exact any 7 / final 6 / body 15, named
13/15, novel absent 6/8/3, all genes 247/276; D-real C4 3/2/4, 12/15, 2/1/3, 257/276; D-real D4 6/6/6, 13/15,
1/0/2, 255/276; D-real B 4/4/4, 250/276. No difference.

## Item 6: development results (run main2, jobs 55585119 + post 55585621)

Per gene (`results/compact_main2.md`); arm B/D exists only for D-real.

| read set | gene | IPD-absent | body A / C4 / C4-R4n | byte-exact final A / C4 | 4-field A / C4 / C4-R4n | novel on present C4 | novel absent right/wrong C4 -> C4-R4n | body B / D4 / D4-R4n | 4-field B / D4 | right/wrong D4 -> D4-R4n |
|---|---|---|---|---|---|---|---|---|---|---|
| D-sim | A | 0 | 0/0/0 | 0/0 | 39/39/39 of 40 | 0 | 0/0 -> 0/0 | | | |
| D-sim | B | 0 | 0/0/0 | 0/0 | 40/40/40 of 40 | 0 | 0/0 -> 0/0 | | | |
| D-sim | C | 1 | 0/0/0 | 0/0 | 36/37/37 of 39 | 0 | 0/0 -> 0/0 | | | |
| D-sim | DPA1 | 5 | 0/0/0 | 0/0 | 15/15/15 of 34 | 0 | 0/0 -> 0/0 | | | |
| D-sim | DPB1 | 3 | 0/1/0 | 0/1 | 36/36/36 of 36 | 0 | 1/0 -> 0/0 | | | |
| D-sim | DQA1 | 6 | 0/3/3 | 0/2 | 30/32/32 of 34 | 0 | 2/0 -> 2/0 | | | |
| D-sim | DQB1 | 1 | 0/0/0 | 0/0 | 38/38/38 of 38 | 0 | 0/0 -> 0/0 | | | |
| D-sim | DRB1 | 25 | 0/15/15 | 0/6 | 10/13/13 of 15 | 0 | 14/3 -> 14/1 | | | |
| D-real | A | 0 | 0/0/0 | 0/0 | 40/40/40 of 40 | 0 | 0/0 -> 0/0 | 0/0/0 | 40/40 | 0/0 -> 0/0 |
| D-real | B | 0 | 0/0/0 | 0/0 | 39/39/39 of 40 | 0 | 0/0 -> 0/0 | 0/0/0 | 38/38 | 0/0 -> 0/0 |
| D-real | C | 1 | 0/0/0 | 0/0 | 37/38/38 of 39 | 0 | 0/0 -> 0/0 | 0/0/0 | 37/38 | 0/0 -> 0/0 |
| D-real | DPA1 | 5 | 0/0/0 | 0/0 | 28/29/29 of 34 | 0 | 0/0 -> 0/0 | 1/1/1 | 25/26 | 0/0 -> 0/0 |
| D-real | DPB1 | 3 | 0/1/1 | 0/0 | 35/35/35 of 36 | 0 | 0/0 -> 0/0 | 1/2/2 | 35/35 | 0/0 -> 0/0 |
| D-real | DQA1 | 6 | 0/4/4 | 0/3 | 28/30/30 of 34 | 0 | 3/0 -> 3/0 | 5/6/6 | 29/31 | 1/0 -> 1/0 |
| D-real | DQB1 | 1 | 0/0/0 | 0/0 | 38/38/38 of 38 | 0 | 0/0 -> 0/0 | 0/0/0 | 38/38 | 0/0 -> 0/0 |
| D-real | DRB1 | 25 | 0/4/4 | 0/2 | 8/12/12 of 15 | 0 | 3/3 -> 3/0 | 4/6/6 | 8/13 | 1/2 -> 1/0 |

All genes (`results/tables_main2_full.md`):

| read set | arm | body | exact any/final | 4-field | novel on present | novel absent right/wrong | unresolved haplotypes (with patched call) | mean wrong bases per copy all / absent | novel closer/farther/same |
|---|---|---|---|---|---|---|---|---|---|
| D-sim | A | 0 | 0/0 | 244/276 | 0 | 0/0 | | 5.4 / 2.6 | |
| D-sim | C4 | 19 | 11/9 | 250/276 | 0 | 17/3 | | 5.2 / 1.3 | 20/0/0 |
| D-sim | C4R1 | 11 | 9/7 | 250 | 0 | 9/2 | 12 (10) | 5.2 / 1.7 | 11/0/0 |
| D-sim | C4R2 | 10 | 8/6 | 249 | 0 | 8/2 | 39 (21) | 5.2 / 1.9 | 10/0/0 |
| D-sim | C4R3 | 5 | 4/2 | 248 | 0 | 3/0 | 221 (30) | 5.3 / 2.5 | 3/0/0 |
| D-sim | C4R4n | 18 | 10/8 | 250 | 0 | 16/1 | 113 (5) | 5.2 / 1.5 | 17/0/0 |
| D-real | A | 0 | 0/0 | 253 | 0 | 0/0 | | 5.4 / 3.0 | |
| D-real | C4 | 9 | 8/5 | 261 | 0 | 6/3 | | 5.2 / 2.5 | 9/0/0 |
| D-real | C4R1 | 8 | 8/5 | 261 | 0 | 5/2 | 12 (2) | 5.2 / 2.6 | 7/0/0 |
| D-real | C4R2 | 8 | 8/5 | 259 | 0 | 5/2 | 32 (10) | 5.3 / 2.6 | 7/0/0 |
| D-real | C4R3 | 6 | 6/3 | 257 | 0 | 3/0 | 203 (16) | 5.3 / 2.9 | 3/0/0 |
| D-real | C4R4n | 9 | 8/5 | 261 | 0 | 6/0 | 118 (6) | 5.3 / 2.7 | 6/0/0 |
| D-real | B | 11 | 11/11 | 250 | 0 | 0/0 | | 6.2 / 2.3 | |
| D-real | D4 | 15 | 15/14 | 259 | 0 | 2/2 | | 6.1 / 2.0 | 4/0/0 |
| D-real | D4R4n | 15 | 15/14 | 259 | 0 | 2/0 | 117 (5) | 6.1 / 2.2 | 2/0/0 |

Every novel final call, wrong ones included, is closer to the truth than the A/B call it replaces. DQB1's mean of 39
wrong bases comes from one copy far off in every arm.

Patch verdicts (`results/verdicts_main2.md`): every applied patch of the final-call kind (`own_version_unsupported`)
is right at patch level (C arm: C 2+3, DPA1 5+3, DPB1 0+1, DQA1 13+12, DRB1 13+26 for D-real+D-sim; D arm: C 2,
DPA1 7, DQA1 6, DRB1 13). Guard rejections: DRB1 D-real 15 (C) and 17 (D) right; C D-real 1 right, 1 neither;
wrong rejections C D-sim 1, DPA1 D-sim 2. `shared_with_partner` patches whose truth is the own version: DPB1 D4 44,
DPA1 D4 11, DPA1 D-sim C 8, DRB1 D-sim C 6 (alternatives only; no final call changes). A and B: no events at all.

## Item 3: abstain rule

Statuses (from the builder only): guard (confounded / near-blind rejection), ambig (two supported versions /
pile-up), still (own version still supported), alt (alternative kept), blind (>= 100 bases without own reads),
untestable (a family version > 120 bases differs, untestable with reads).

| rule | unresolved when | C4 novel absent right/wrong (D-sim; D-real) | D4 | cost |
|---|---|---|---|---|
| none | | 17/3; 6/3 | 2/2 | |
| R1 (planned) | guard or ambig | 9/2; 5/2 | 2/2 | loses 9 right to drop 2 wrong (+T DRB1*09 copies always carry ambig/still notes) |
| R2 (planned) | R1 or still or alt | 8/2; 5/2 | 2/2 | as R1 plus rename losses (261 -> 259) |
| R3 (planned) | R2 or blind or untestable >= 1 | 3/0; 3/0 | 1/0 | flags 221/297 haplotypes |
| R4 (post hoc) | guard >= 1 or untestable >= 2 | 16/1; 6/0 | 2/0 | loses 1 right (DPB1 D-sim) |

`n` variants spare renames; dropping renames costs four-field accuracy with no gain. Recommendation: R4n. Caveats:
fitted to 6 wrong cases; flags almost every class I haplotype (A 118/120), so it effectively turns off novel class I
calls; the remaining wrong call (D-sim HG03710 DRB1) has only ambig/still notes.

## Item 7: run time

Per task (8 genes) median 35 s, max 95 s. Per gene median/max s: A 5.2/14.6, B 5.5/16.3, C 7.7/24.4, DPA1 2.6/6.8,
DPB1 6.1/18.4, DQA1 1.7/5.2, DQB1 2.6/6.5, DRB1 4.1/16.0. Loading+mapping median 2.4 s (954 s during one stall).
Test projection: 254 tasks x ~40 s = ~11 core-hours; 1-2 h wall with %12.

## Open problems (development)

- Hang cause shown by a live reproduction, not by a stack of the original job.
- Class I recovery untested on development (0/0/1 IPD-absent copies for A/B/C).
- `shared_with_partner` fires often against the truth for DPB1/DPA1 on the cohort panel (alternatives only).
- Edit distance uses its own min-sum copy-to-slot assignment and the first label of an ambiguous slot.
- No intervals for the per-gene development tables.
- I ran two tiny `python3 -c` commands on the IBEX login node (a no-op and one JSON read) before switching to grep.

## Development run log

| job | what |
|---|---|
| 55582726, 55582807-55582811 | hang probes |
| 55582956 | smoke test |
| 55583286 / 55584344 | regress + post |
| 55583287 | learn with cap 3,000 (15 gene failures, discarded) |
| 55584444 / 55584584 | learn + confounded lists |
| 55584750 / 55584753 / 55585075 | main / regress2 / post |
| 55585119 / 55585120 / 55585621 | main2 (final code) / regress3 / post |

All tasks wrote COMPLETE.json; gene-level failures only in the discarded run.

# Test (frozen code, FREEZE.json 2026-11-14T19:16:03+03:00)

Design and pass criteria: `DESIGN.md` (written and hashed into `FREEZE.json` before any test-donor task of this
stage; local and IBEX copies identical). Test array 55585788 (254 tasks: T-sim 61 + T-real 66, arms C and D, %12,
1 h limit); retry check `resubmit_timeouts.sh` found no timeout (0 retries); scoring job 55590531
(`ag_proposals_table.py`, `ag_score.py`, `ag_test_summary.py`, all frozen). 254/254 tasks COMPLETE, 2,032/2,032
gene runs complete, 0 FAILURE records, 1,016/1,016 arm tables scored without fallback. All numbers below are from
`results/TEST_SUMMARY_test.md`, recomputed from `score_rows_test.tsv` (16,248 rows) and `tasks_test.tsv`.
The truth of the test donors was opened only by the frozen scorer after the run. Patch-level verdicts were not
computed on test (the frozen `ag_eval_patches.py` refuses test read sets).

## Pass / fail per criterion

| criterion | T-sim | T-real |
|---|---|---|
| a. novel final calls on IPD-present copies, all genes (<= 1) | C4-R4n 1, D4-R4n 1 (both DRB1*16:02:01:03, 1 base off) -> PASS | 0 / 0 -> PASS |
| b. four-field named, total not lower, no gene -1 or worse | C4-R4n vs A 707 -> 723 of 804 (+2.0 pp [+1.1, +3.0]); D4-R4n vs B 694 -> 711 (+2.1 pp [+1.2, +3.2]); no gene lower -> PASS | C4-R4n 861 -> 873 of 914 (+1.3 pp [+0.7, +2.0]); D4-R4n 844 -> 858 (+1.5 pp [+0.9, +2.3]); no gene lower -> PASS |
| c. primary: gene body of IPD-absent copies, all genes, C4-R4n minus A | 0 -> 22 of 162, +13.6 pp [+8.8, +18.7], p = 1.5e-05 -> PASS | 0 -> 12 of 134, +9.0 pp [+4.8, +13.4], p = 0.0005 (reported) |
| d. complementarity, D4-R4n minus B | 39 -> 56 of 162, +10.5 pp [+6.1, +15.1] | 31 -> 38 of 134, +5.2 pp [+2.1, +9.1] |
| e. R4n keeps >= 75% of right novel calls | C4 18/34 right (53%) -> R4n 15/23 (65%), kept 15/18 = 83% PASS; D4 12/22 (55%) -> 10/14 (71%), kept 10/12 = 83% PASS | C4 10/42 (24%) -> 8/22 (36%), kept 8/10 = 80% PASS; **D4 3/25 (12%) -> 2/12 (17%), kept 2/3 = 67% FAIL** |
| f. fallbacks <= 2 of 254 | 0 of 254 -> PASS | |

Criterion e fails in one of four cells (T-real D4: one of three right calls lost; small numbers). All other criteria
pass. Intervals: family-clustered bootstrap 95%, p = family sign-flip (frozen `score.uncertainty`).

## Per gene (test)

| read set | gene | IPD-absent | body A / C4 / C4-R4n | exact final A / C4-R4n | 4-field A / C4-R4n (of n) | novel present C4-R4n | right/wrong C4 -> C4-R4n | body B / D4 / D4-R4n | exact final B / D4-R4n | 4-field B / D4-R4n | novel present D4-R4n | right/wrong D4 -> D4-R4n |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| T-sim | A | 2 | 0 / 0 / 0 | 0 / 0 | 112 / 113 (of 119) | 0 | 0/1 -> 0/0 | 0 / 0 / 0 | 0 / 0 | 112 / 113 | 0 | 0/1 -> 0/0 |
| T-sim | B | 12 | 0 / 0 / 0 | 0 / 0 | 105 / 106 (of 109) | 0 | 0/0 -> 0/0 | 0 / 0 / 0 | 0 / 0 | 102 / 103 | 0 | 0/0 -> 0/0 |
| T-sim | C | 6 | 0 / 1 / 1 | 0 / 0 | 107 / 107 (of 115) | 0 | 0/0 -> 0/0 | 0 / 1 / 1 | 0 / 0 | 104 / 104 | 0 | 0/0 -> 0/0 |
| T-sim | DPA1 | 16 | 0 / 1 / 1 | 0 / 0 | 71 / 73 (of 104) | 0 | 1/0 -> 1/0 | 2 / 3 / 3 | 0 / 0 | 68 / 70 | 0 | 1/1 -> 1/1 |
| T-sim | DPB1 | 15 | 0 / 3 / 3 | 0 / 0 | 97 / 97 (of 104) | 0 | 0/1 -> 0/1 | 0 / 3 / 3 | 0 / 0 | 95 / 95 | 0 | 0/1 -> 0/1 |
| T-sim | DQA1 | 32 | 0 / 4 / 4 | 0 / 2 | 76 / 83 (of 90) | 0 | 3/0 -> 3/0 | 20 / 24 / 24 | 3 / 4 | 75 / 83 | 0 | 3/0 -> 3/0 |
| T-sim | DQB1 | 3 | 0 / 0 / 0 | 0 / 0 | 103 / 104 (of 118) | 0 | 0/0 -> 0/0 | 0 / 0 / 0 | 0 / 0 | 103 / 104 | 0 | 0/0 -> 0/0 |
| T-sim | DRB1 | 76 | 0 / 16 / 13 | 0 / 7 | 36 / 40 (of 45) | 1 | 14/14 -> 11/7 | 17 / 27 / 25 | 16 / 21 | 35 / 39 | 1 | 8/7 -> 6/2 |
| T-real | A | 2 | 0 / 0 / 0 | 0 / 0 | 128 / 129 (of 130) | 0 | 0/0 -> 0/0 | 0 / 0 / 0 | 0 / 0 | 128 / 129 | 0 | 0/0 -> 0/0 |
| T-real | B | 8 | 0 / 0 / 0 | 0 / 0 | 114 / 115 (of 121) | 0 | 0/0 -> 0/0 | 1 / 1 / 1 | 0 / 0 | 110 / 111 | 0 | 0/0 -> 0/0 |
| T-real | C | 2 | 0 / 0 / 0 | 0 / 0 | 126 / 126 (of 128) | 0 | 0/0 -> 0/0 | 1 / 1 / 1 | 1 / 1 | 126 / 126 | 0 | 0/0 -> 0/0 |
| T-real | DPA1 | 10 | 0 / 1 / 1 | 0 / 0 | 111 / 112 (of 122) | 0 | 0/0 -> 0/0 | 1 / 2 / 2 | 0 / 0 | 102 / 103 | 0 | 0/0 -> 0/0 |
| T-real | DPB1 | 12 | 0 / 5 / 3 | 0 / 1 | 117 / 118 (of 118) | 0 | 3/0 -> 1/0 | 1 / 5 / 4 | 1 / 2 | 117 / 118 | 0 | 2/0 -> 1/0 |
| T-real | DQA1 | 15 | 0 / 3 / 3 | 0 / 3 | 106 / 107 (of 117) | 0 | 3/0 -> 3/0 | 12 / 13 / 13 | 12 / 13 | 105 / 106 | 0 | 1/0 -> 1/0 |
| T-real | DQB1 | 0 | 0 / 0 / 0 | 0 / 0 | 128 / 128 (of 131) | 0 | 0/0 -> 0/0 | 0 / 0 / 0 | 0 / 0 | 127 / 127 | 0 | 0/0 -> 0/0 |
| T-real | DRB1 | 85 | 0 / 5 / 5 | 0 / 2 | 31 / 38 (of 47) | 0 | 4/32 -> 4/14 | 15 / 17 / 17 | 14 / 15 | 29 / 38 | 0 | 0/22 -> 0/10 |

## Secondary (all genes and seven non-DRB1 genes)

All genes / non-DRB1 (from `TEST_SUMMARY_test.md`):

| read set | arm | body (all; non-DRB1) | exact final (all; non-DRB1) | novel absent right/wrong (all; non-DRB1) | novel closer/equal/farther (all) |
|---|---|---|---|---|---|
| T-sim | C4 | 25; 9 | 12; 2 | 18/16; 4/2 | 34/0/0 |
| T-sim | C4-R4n | 22; 9 | 9; 2 | 15/8; 4/1 | 24/0/0 |
| T-sim | C4-R1n | 24; 9 | 11; 2 | 17/13; 4/2 | 30/0/0 |
| T-sim | D4 | 58; 31 | 27; 4 | 12/10; 4/3 | 22/0/0 |
| T-sim | D4-R4n | 56; 31 | 25; 4 | 10/4; 4/2 | 15/0/0 |
| T-sim | D4-R1n | 57; 31 | 26; 4 | 11/6; 4/2 | 18/0/0 |
| T-real | C4 | 14; 9 | 8; 6 | 10/32; 6/0 | 41/0/1 |
| T-real | C4-R4n | 12; 7 | 6; 4 | 8/14; 4/0 | 21/0/1 |
| T-real | C4-R1n | 13; 9 | 8; 6 | 9/25; 6/0 | 33/0/1 |
| T-real | D4 | 39; 22 | 32; 17 | 3/22; 3/0 | 24/0/1 |
| T-real | D4-R4n | 38; 21 | 31; 16 | 2/10; 2/0 | 11/0/1 |
| T-real | D4-R1n | 39; 22 | 32; 17 | 3/18; 3/0 | 20/0/1 |

Non-DRB1 gene-body recovery (secondary, first look): C4-R4n minus A T-sim 0 -> 9 of 86, +10.5 pp [+5.0, +16.1];
T-real 0 -> 7 of 49, +14.3 pp [+5.7, +25.0]. D4-R4n minus B T-sim 22 -> 31, +10.5 pp [+5.2, +15.8]; T-real 16 -> 21,
+10.2 pp [+2.2, +19.6]. HLA-A, HLA-B and DQB1 gain no gene bodies; gains are in DQA1, DPB1, DPA1, C.

Mean wrong bases per copy (gene body, paired, family-clustered), arm minus base:
- T-sim C4-R4n - A: all copies -0.13 [-0.17, -0.10]; IPD-absent -0.50 [-0.72, -0.30]; non-DRB1 all -0.05 [-0.07, -0.04].
- T-sim D4-R4n - B: all -0.12 [-0.16, -0.08]; IPD-absent -0.41 [-0.67, -0.21].
- T-real C4-R4n - A: all -0.11 [-0.17, -0.06]; IPD-absent -0.67 [-1.11, -0.28]; non-DRB1 all -0.02 [-0.03, -0.01].
- T-real D4-R4n - B: all -0.09 [-0.14, -0.04]; IPD-absent -0.44 [-0.88, -0.08].
Byte-exact (final) recovery, C4-R4n - A: T-sim +9 of 162 (+5.6 pp [+2.8, +8.6]); T-real +6 of 134 (+4.5 pp
[+1.5, +7.9]). D4-R4n - B: T-sim +6 (+3.7 pp [+0.7, +6.7]); T-real +3 (+2.2 pp [+0.0, +5.0]).

## Notes on the test

- Real reads, DRB1: most novel calls on IPD-absent DRB1 copies are still not the true gene body (C4 4/32 right).
  R4n halves the wrong ones (32 -> 14, 22 -> 10) but they remain the majority. Outside DRB1 every T-real novel call
  on an IPD-absent copy is right (C4-R4n 4/0, D4-R4n 2/0).
- The single novel call on an IPD-present copy (T-sim, DRB1*16:02:01:03, 1 base off) is present in every arm,
  including R1n and R4n (R4 does not flag it).
- One T-real novel call (all C/D variants) is farther from the truth than the arm-A/B call; all others are closer.
- Two T-sim copies have no edit distance (no final call matched; listed as missing).
- Abstaining from renames is not tested here (all primary arms keep renames).
- Run time: T-sim tasks mean 46 s (max 108 s), T-real mean 60 s (max 146 s).

## Supervisor check and decision (2026-11-14)

- Recounted from `results/score_rows_test.tsv` (column `body_recovered`, `named4`, `false_mos`, `final_mos`):
  all headline counts in the Test section reproduce (T-sim 0 -> 22 and 39 -> 56 of 162; T-real 0 -> 12 and
  31 -> 38 of 134; four-field 707 -> 723, 694 -> 711, 861 -> 873, 844 -> 858; 254 of 254 tasks complete).
- Freeze order checked on IBEX: `FREEZE.json` written 19:14:07, test array 55585788 submitted 19:14:21;
  no `src/` file is newer than 19:09:44.
- Criterion e stays recorded as **failed** for T-real arm D4 (2 of 3 right calls kept). No rule is changed.
  The rule R4 was chosen after seeing development outcomes; the unfiltered arms C4/D4 are reported beside it.
- Open weakness: on real reads most novel DRB1 calls are still not the true gene body (C4-R4n 4 right, 14 not).
  A novel call must be read as "this copy differs from the catalogue", not as a finished sequence.
