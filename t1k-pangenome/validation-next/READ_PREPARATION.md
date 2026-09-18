# HGSVC28 read preparation

Submitted as IBEX array **52061497**, indices 0–27, all 28 tasks observed RUNNING
on 18 September 2026. Each receives 4 CPUs and 8 GB RAM, with a four-hour walltime.
The walltime is a resource limit, not an estimated completion time.

The fixed 28-donor reservation is unchanged. Public CRAM URLs come from the
previously recorded availability audit; the staging manifest additionally binds
that audit to the reserved cohort hash. Transferred code and inputs were checked
against local SHA-256 hashes before submission. The launcher records submission
intent before calling Slurm, so an interrupted acknowledgement must be reconciled
rather than resubmitted blindly.

Output namespace inside the existing container:
`/home/leechuck/hla/t1k-pangenome/hgsvc-validation-v1/validation-reads/`.
Host storage remains under `/ibex/scratch/projects/c2014/rob/dogohla-benchmark/hla/`.
Slurm uses the working `pi-hohndor` allocation.

Recruitment retains the existing benchmark's MHC, chr6 alternate/HLA contig and
unmapped-pair selection. Each task writes its own region BED, avoiding a shared
file race. The existing validation driver checks the reserved donor/source list,
CRAM read-group sample identity, FASTQ syntax, nonempty identical mate-name lists,
and final mate-file hashes. A COMPLETE recruitment marker alone is insufficient:
`VERIFIED.json` is required for downstream processing.

This is read preparation only. No HGSVC truth labels or candidate typing outputs
are consumed. The HGSVC endpoint implementation, complete inference freeze and
actual graph-path exclusion audit remain prerequisites for the next evaluation.
Original frozen validation and post-hoc development results remain separate.

Validation: four existing read-verification unit tests passed; Python compilation
and shell syntax checks passed. Live job state confirmed all 28 tasks RUNNING;
read completion has not yet been established.
