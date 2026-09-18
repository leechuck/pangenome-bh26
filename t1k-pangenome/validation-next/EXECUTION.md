# HGSVC follow-up execution

`monitor_hgsvc_validation.py` runs as the local user service
`dogohla-hgsvc28-monitor.service`. It polls every 60 seconds. The service was
confirmed active on 18 September 2026 with 19/28 read preparations verified;
the remaining nine were still running. No typing results have been inspected.

When all donor reads have VERIFIED.json records, `freeze_hgsvc_validation.py`
checks the selected candidate, truth-source hashes and read provenance, then
writes separate inference and evaluation freezes. It reuses the previous graph
asset manifests after verifying their remote hashes; the source-path audit binds
them to the same observed panels. Every execution stage checks the freeze.

`launch_hgsvc_validation.py` stages and hashes the code and submits original T1K,
genomic-IPD, coverage, native evidence, bootstrap mapping/projection and library
calibration. `advance_hgsvc_validation.py` subsequently submits calibrated mapping,
projection, locus joins, guarded refinement and coarse anchoring. All 28 donors
are retained. Mapping and projection each cover 448 donor/panel/locus tasks.
Mapping arrays allow 200 simultaneous tasks, four CPUs each. This is a concurrency
limit, not a promise of instantaneous allocation. Jobs use pi-hohndor and the
batch/debug partitions, with two-hour walltimes inherited from the validated
pipeline. Read preparation retains its separate four-hour limit.

The launcher records intent before each submission and the job ID immediately
after acknowledgement. An uncertain submission must be reconciled using Slurm;
it must not be retried blindly. The dispatcher advances from complete task-level
Slurm records, not elapsed time, and stops on detected terminal failures. A failed
monitor does not imply failed Slurm jobs: inspect the ledger, live queue and
accounting before recovery. The monitor never evaluates calls.

New execution paths are under `t1k-pangenome/hgsvc-validation-v1` on IBEX. The old
frozen execution drivers and output namespaces remain unchanged. New graph-driver
differences are limited to this namespace and the selected native-locus refinement
adapter; baseline and genomic drivers change only read locations. A shared stage
entry point verifies code, metadata and asset-manifest hashes. Coarse anchoring
also verifies original/candidate output hashes and identical input-read hashes.

Nine HGSVC tests passed: exact whole-core truth matching, rejected missing,
duplicate and ambiguous truth; strict name resolution; complete failure accounting;
tamper detection; full two-panel mapping coverage; prerequisite completeness; and
protection against duplicate uncertain submissions. Existing coarse-anchor and
native-locus tests also passed. Live inference success and independent improvement
remain to be established. Full provenance review and the complete terminal method
grid are required before scoring.
