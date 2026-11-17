# All-gene patch builder: design and pre-registered test (stage mosaic-allgenes-20261114)

Follows `mosaic-drb1-20261113` (frozen DRB1 builder v3.1). Development: 20 donors (D-sim, D-real); see REPORT.md.
Test: T-sim (61 simulated Arab/Saudi donors) and T-real (66 real HPRC Asian donors), never pooled. DRB1 on these
test donors is a third look (frozen arms, then the DRB1 builder test); the seven other genes are scored with the
builder for the first time.

## 1. What changed versus the frozen DRB1 stage

- Gene is a parameter everywhere (reads, families, events, confounded list, structural filter, scoring). PARAMS
  unchanged: max_edit 2, sup_share 0.25, min_sup 3, cur_ratio 0.25, het_share 0.3, pile 2.5, max_indel_cols 150,
  hom_ratio 1.5, blind_min_run 100, blind_near 150, blind_max_nm 2, n_max 4, family_theta 0.01.
- Reads are mapped once per donor x arm to the pass-1 alleles of all genes (the target never depended on the gene).
- Truth-free confounded lists relearned per gene from D-sim arm C (`inputs/confounded_v1/`; DRB1 list identical to
  the frozen one).
- Robustness: 1800 s wall-clock limit on loading+mapping and per gene; deterministic work caps per gene. Each yields
  a FAILURE record and the pass-1 call as fallback. (The DRB1 test hang was a Lustre read stall, not the code.)
- No T1K second pass (C3/D3 dropped). Arms: C4 = arm A pass-1 + builder (IPD panel); D4 = arm B pass-1 + builder
  (IPD + cohort panel of the donor's fold).
- With the abstain rule off, the DRB1 output reproduces the frozen DRB1 stage exactly on development (REPORT.md).

## 2. Abstain rules (truth-free, from builder statuses)

Per haplotype: guard (patch rejected as confounded or near a blind region), ambig (two supported versions /
pile-up), still (own version still supported), alt (alternative kept), blind (>= 100 bases without own reads),
untestable (a family version longer than 120 bases differs). A flagged haplotype keeps its pass-1 named call;
suffix `n` = only novel final calls are withheld, renames to a known allele stay.

- R1 (planned before development outcomes): guard >= 1 or ambig >= 1.
- R4 (chosen AFTER seeing development outcomes; post hoc): guard >= 1 or untestable >= 2.

Primary arms: C4-R4n and D4-R4n. Secondary, reported beside them: C4 / D4 (no abstaining) and C4-R1n / D4-R1n.

## 3. Execution rules fixed before the test

- 254 tasks: (T-sim 61 + T-real 66) x arms C and D, `%12`, `--time=01:00:00`, frozen `src/` only.
- A task that times out (task- or gene-level timeout, or killed without a record) is resubmitted once on a different
  node (`src/resubmit_timeouts.sh`, first node excluded). A second timeout is not retried: its fallback (pass-1
  calls) is scored and counted as a failure.
- Failure counting (criterion f): a task ends in a fallback if the task did not complete, or any gene ended with
  failure_kind timeout / bound / error, or the scorer had to use the pass-1 table. Genes without a pass-1 call
  (`no_pass1_call`) are listed separately, not counted as fallbacks.
- Scoring: `ag_score.py` (frozen scorer `score.py` underneath), then `ag_proposals_table.py` and
  `ag_test_summary.py`; family-clustered bootstrap and sign-flip test from the frozen `score.uncertainty`.

## 4. Pass criteria (per read set, T-sim and T-real separately)

a. Novel final calls on IPD-present copies, all genes: at most 1 per read set, for C4-R4n and for D4-R4n.
b. Four-field named accuracy, C4-R4n vs A and D4-R4n vs B: all-gene total not lower; no gene loses more than 1 copy.
c. Primary: gene-body recovery of IPD-absent copies over all eight genes, C4-R4n minus A: 95% interval above 0 on
   T-sim; reported with interval on T-real.
d. Complementarity: D4-R4n minus B, same endpoint, reported with interval for both read sets.
e. Of novel final calls on IPD-absent copies, the share with the true gene body, with and without R4n; R4n must keep
   at least 75% of the right calls (C4 and D4).
f. At most 2 of 254 tasks end in a fallback.

Secondary (labelled as such): the same endpoints for the seven non-DRB1 genes only; per-gene counts; byte-exact
recovery; mean wrong bases per copy (edlib global edit distance of the final call's gene body to the true gene body,
min-sum copy-to-slot assignment) with paired family-clustered intervals; closer/equal/farther counts for novel
final calls versus the base arm.
