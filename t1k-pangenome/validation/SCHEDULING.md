# Validation scheduling audit

The frozen inference and evaluation code remain unchanged. This note records
scheduler-only changes made on 18 September 2026 to avoid unnecessary barriers.

- The pilot refinement completed successfully before donor-specific dependencies
  were introduced for the remaining cohort.
- Of 747 attempted dependency updates, 107 succeeded: 24 bootstrap mapping tasks
  and all 83 donor calibration tasks. Each calibration task requires its own
  eight bootstrap projection tasks to finish successfully.
- The other 640 mapping dependency updates were rejected while the jobs were
  starting, completing, or disappearing from the controller cache. A subsequent
  `sacct -X` reconciliation verified **every one of those 640 tasks as
  COMPLETED with exit code 0:0**. None requires resubmission.
- All 664 cohort bootstrap mappings subsequently completed successfully.
- Projection task 52051990_200 was still pending on its `aftercorr` dependency
  even though matching mapping task 52051978_200 had completed successfully.
  Pending projection dependencies were cleared only after their exact matching
  mapping task is verified COMPLETED, exit code 0:0, in Slurm accounting.
  Running jobs are not restarted and no inference parameter is changed.

Detailed scheduler responses are retained on IBEX under
`/ibex/scratch/projects/c2014/rob/dogohla-benchmark/` in
`validation-donor-dependency-update.json` and
`validation-projection-dependency-release.json`.
The local donor-update audit additionally contains the completed-task
reconciliation. These are execution records, not genotype results.

The automatic dispatcher checks predecessor completion; projection can also use
matching-task Slurm dependencies, as documented below.
Prediction retrieval and scoring remain gated on the complete comparison grid.

## Verified outcome

All 664 bootstrap mappings, all 664 bootstrap projections and all 83 donor
calibrations completed with exit code 0:0. The projection release attempted
216 updates: 32 succeeded; 184 were rejected because tasks had already started
or finished. Every affected task was subsequently verified complete, including
all rejected updates. The detailed completion grid is retained in
`preparation-execution-complete.json`.

The dispatcher then submitted calibrated mapping array **52055293**, indices
16–1343 (1,328 tasks: 83 donors × eight loci × two panels), capped at 300
concurrent four-CPU tasks on `batch,debug`. The next stages remain automatically
gated on successful predecessor completion. No validation outcomes have been
inspected.

## Overlapping calibrated mapping and projection

Projection array **52056935** uses indices 16–1343 and
`aftercorr:52055293`, so each projection waits for successful completion of its
matching donor/locus/panel mapping. This removes the whole-cohort mapping
barrier while retaining the same prerequisites. Both arrays use identical
indices; projection is capped at 300 concurrent one-CPU tasks. The join stage
still requires all 1,328 projections and all native exports to complete.

The dispatcher regression checks passed (2 tests), including a check that
projection submission supplies the matching-task dependency and that an
incomplete projection cohort cannot release the join stage. Frozen inference,
read inputs, and endpoint scoring are unchanged.

The controller initially retained the `aftercorr` dependency even for matching
mappings already recorded complete. An accounting-verified release was therefore
started: only pending projection task IDs whose exact mapping ID was COMPLETED
with exit code 0:0 were passed to `scontrol update ... Dependency=`. This produced
actual overlap (47 projection tasks running while 10 mappings were still running).
The release SSH connection ended with exit 255 before its audit file was
written. No task was restarted. Subsequent independent accounting verified all
1,328 mappings and all 1,328 projections COMPLETED with exit code 0:0; see
`mapping-projection-execution-complete.json`. No dependency update or retry is
needed for these completed arrays.

## Join submission reconciliation

Join array **52058649**, indices 2–167, was accepted by Slurm but the monitor's
55-second dispatcher timeout expired before its job ID reached the local ledger.
Slurm accounting's `SubmitLine` verified the exact stage, index range and frozen
inference digest. The existing job was adopted into the ledger; no duplicate
was submitted. The monitor now allows the dispatcher 300 seconds to finish
submission and record its result. Its unknown-submission guard is retained.
