# Implementation checkpoint — 2026-09-18

The workflow review is in `review-2026-09-18/REVIEW.md`.

In progress:

- `phase_linkage.py` corrects maximum shared-allele scoring and best-HSP selection.
- `experiment.py` creates isolated original/corrected SpecHLA script copies,
  same-release IPD-only/panel phase databases, and four phase experiment arms.
- Completion requires validated output hashes and matching configuration;
  failed/partial directories are never silently reused.
- `jobs/phase_experiment.sbatch` schedules the paired development runs.

The five focused local regression/integrity tests pass. These new experiment
runners have not yet been validated on DDBJ at this checkpoint. Graph diagnostics,
structural controls, result scoring and end-to-end verification remain in progress.

The repository keeps code, small metadata, result tables, reviews and slides.
Downloaded IPD sequence databases, extracted loci and native run snapshots remain
local/on DDBJ and are excluded through `.gitignore`; existing fetch/build scripts
and source version manifests record the inputs.
