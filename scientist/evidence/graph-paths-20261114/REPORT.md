# Step 5: graph path picking (development donors only)

Stage `graph-paths-20261114` (IBEX `/encrypted/e3008/hohndor-dogohla-20261020/graph-paths-20261114`).
20 development donors, read sets D-sim and D-real; test donors not touched. Written by the supervisor from
the implementing agent's hand-back; placement and sampling counts recounted from `results/placement_rows.tsv`
and `results/sampling_rows.tsv`. Full tables: `results/TABLES.txt`.

## Outcome

1. **A clean acyclic graph fixes read placement.** Graphs built straight from the per-gene alignment (IPD,
   or IPD + the fold's cohort copies; no loops). Giraffe places reads that carry the true local version
   exactly: DRB1 98% simulated / 96% real (old PGGB graph 41% / 52%); intron-2 repeat 99% / 98% (45% / 54%);
   reads needing a new combination of stored pieces 184/184 and 113/116 (17% / 15%). DQA1 86-91%; B 18/18.
2. **Picking two paths with `vg haplotypes` (v1.76.1) does not give a typer.** DRB1: 0 of 25 IPD-absent
   copies recovered at any setting, and only 1-4 of 15 IPD-present copies returned exactly. Shorter blocks
   make every gene worse (B: 38-39 of 40 exact with whole paths, 4-6 of 40 at 300 bases). DQA1 with whole
   paths on the IPD + cohort graph recovers 5 of 6 absent copies, which only equals arm B (5 of 6).
3. **Fallback (edge coverage + two-route search) was worse than arm A** everywhere (DRB1 0/40 exact).
4. **No frozen test of a graph typer.** The patch builder stays the DRB1 method.

## Graphs
18 graphs; each stored sequence is one path, spelled exactly, acyclic. DRB1: 912 IPD paths + 113-123 cohort,
~8,900 nodes. Intron-1 homopolymer: 7-8 parallel branches. Intron-2 repeat: 174-206 versions as 221-269 routes.

## Not explained
Why whole-path DRB1 sampling replaces 14 of 15 IPD-present copies with another stored allele (vg's own
coverage estimate and scoring options were not tried). 90 of 480 fallback rows were Giraffe 900 s timeouts.

## Problems during the run
The first donor arrays hit the time limit (Giraffe single-threaded on full candidate read sets). Fixed with
hard timeouts on every vg call, per-setting checkpoints, fragment-length settings, and giving Giraffe only
read pairs that contain a tested local version. Final array 55600566: 40/40, no failed setting.

## Possible follow-up (not part of this step)
Use clean-graph placement as the read-evidence source for the patch builder (reads are placed correctly at
exactly the repeat sites where the builder is weakest on real reads).
