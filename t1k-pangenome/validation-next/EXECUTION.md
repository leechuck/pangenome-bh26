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

## Provenance handoff and scoring

`hgsvc_handoff.py` checks all seven outputs per donor (original T1K, three raw
candidates and three anchored candidates) before fetching any prediction text.
It verifies the complete grid, read identities, frozen model/driver hashes,
parameters, reference manifests and parent-output links. The remote check hashes
small result files but returns only metadata. Failed or pending outputs block
fetching until their disposition is explicitly audited; they are never omitted.

`score_hgsvc_snapshot.py` rechecks downloaded output hashes, recomputes every
anchor from its parents, verifies per-gene two-field identity, and only then
constructs HGSVC truth and invokes the fixed evaluator. It writes an unblinding
record, complete gene-level results, summary, truth QC and paired contrasts.
These handoff/scoring scripts receive their own code freeze when the execution
freeze is written. They are not invoked automatically by the monitor.

Eleven HGSVC tests now pass, including rejection of altered parent links,
swapped reads, missing native-locus guards and incomplete provenance grids.
The scorer has not yet been run on held-out predictions or truth.

## Frozen launch, 18 September 2026

All 28 read jobs completed successfully, yielding **20,772,142 verified pairs**.
The execution freeze is
`4afbe16398deefe3e1a26e7a8b974c7d63219eda1220ccd3bf6919f277d272be`.
Local freeze consistency checks passed before recording this launch.

The first seven arrays are submitted (588 tasks total):

| Stage | Slurm job |
|---|---|
| Original T1K | 52062787 |
| Genomic-IPD T1K | 52062822 |
| Coverage | 52062857 |
| Native evidence | 52062887 |
| Bootstrap mapping | 52062901 |
| Bootstrap projection | 52062902 |
| Library calibration | 52062917 |

All 28 original-T1K, 28 genomic-IPD and 28 coverage tasks were confirmed RUNNING
at the initial check. Dependent arrays are queued behind their prerequisites.
The monitor remains active and will release calibrated mapping and subsequent
stages. This is execution status, not evidence of accuracy; no predictions or
held-out truth have been scored.
