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

The automatic dispatcher submits the remaining calibrated mapping, projection,
join, and refinement stages only after checking every required predecessor.
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
