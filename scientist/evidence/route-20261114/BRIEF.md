# Shared brief for the six-step route (2026-11-14). Read fully before doing anything.

Project: DogoHLA, HLA typing from short reads for Asian populations (PI Robert Hoehndorf, KAUST).
You are an implementer. The supervisor (main session) designs and validates. Follow your step brief
exactly; where it is silent, choose the simplest correct option and record the choice in your REPORT.md.

## Hard rules
- All compute runs on IBEX through `sbatch` (`ssh ibex`, user hohndor). On the login node only
  ls/cat/head/squeue/sacct/sbatch/scp/small edits. No python/samtools/etc. on the login node. No compute on this workstation
  beyond editing files and parsing small summary tables.
- Slurm: `--account=pi-hohndor` only. Give every job the name prefix of your step (`s1-`, `s2-`, ...).
  Never cancel or modify a job that you did not submit yourself; never `scancel -u`. Other sessions run
  jobs under the same user (`yld-*`, `pan-*`, `eqh_*`, and the other steps). Limit arrays to `%12`.
  The `debug` partition has a low per-user CPU cap; use `batch` for arrays. Set realistic `--time`.
- Never write to `/tmp` (here or on IBEX). Use `$TMPDIR` (~/tmp) locally and `<stage>/tmp` on IBEX
  (`export TMPDIR=<stage>/tmp`). Never hardcode `/tmp/...` anywhere.
- Never read or print `/mnt/data1/agents/secrets` or the `remote-connect` skill's `resources/` directory.
- Sequences of donors (truth alleles, reads) stay on IBEX. Only summary tables (counts, positions,
  lengths, allele names, hashes, scores) may be copied to this workstation.
- `confirmatory-test-20261113` is read-only and may be used only where your step brief says so.
- Do not modify any existing stage (code, results, FREEZE files). Import their modules read-only or copy
  code into your own stage. Do not git commit. Do not spawn sub-agents.
- ssh output starts with a banner of `#` lines (`grep -v '^#'`). Local `grep` is ugrep. Foreground
  `sleep N; cmd` is blocked: wait with `until ...; do sleep 30; done` loops.
- No silent failure: every array task writes COMPLETE.json or FAILURE.json; count them; report failures.
- Hygiene: development donors = the 20 in `mosaic-drb1-20261113/work/dev_tasks.txt` (read sets D-sim,
  D-real). Test donors = T-sim (61 simulated Arab/Saudi) and T-real (66 real HPRC Asian). Never look at
  truth of test donors while developing or tuning; only frozen code may be scored on them, and only when
  your brief says so. Simulated and real read sets are never pooled.
- Relatedness: `FOLDS.json` `groups` are family components. A donor's own component must never
  contribute sequences to a panel/graph used to type that donor.

## Places
- IBEX protected root `E=/encrypted/e3008/hohndor-dogohla-20261020` (group ibex-e3008, umask 0002).
  Your stage: `$E/<stage-name>` with `src/ work/ results/ logs/ tmp/`.
- Public scratch `B=/ibex/scratch/projects/c2014/rob/dogohla-benchmark`, `ZEN=$B/hprc-zenodo-20261004`.
- Local mirror of your stage: `/mnt/data1/dogohla/pangenome/scientist/evidence/<stage-name>/`
  (keep `src/` identical to IBEX; put summary tables and REPORT.md there).
- Python on IBEX: `$E/mosaic-read-mapping-20261104/env/bin/python` (numpy, edlib? check; pysam? check).
  Job template to copy: `$E/mosaic-drb1-20261113/src/diag.sbatch`.
- Container for T1K/vg/bwa etc.: `$B/containers/ubuntu24.sif`; see how
  `$E/mosaic-typer-20261113/src/run_donor.py` (`BIND`, `SIF`, `RUN_T1K`, `t1k`) and
  `$E/mosaic-drb1-20261113/src/drb1_run.py` (`t1k`) call it. vg: `/home/leechuck/hla/cactus/cactus-bin-v3.3.0/bin/vg`
  inside the container with `$B/hla` bound at `/home/leechuck/hla` (see `mosaic-drb1-20261113/src/drb1_graphtest.py`).
- Shared library (read-only): `$E/mosaic-typer-20261113/src/common.py` (paths, `load_records` = truth
  records per donor/gene with `core_sequence`, `core_exons`, `genomic_labels`; `load_msas` = per-gene
  alignment of all IPD genomic alleles and cohort copies; `load_classes`; `readset_runs`; T1K table
  readers), `run_donor.py` (`pass1_spec(readset, donor, arm)` gives existing T1K runs and read files),
  `panel_select.py`, `score.py` (frozen scorer with family-clustered bootstrap).
- Read sets and existing runs: see `common.readset_runs()` and `run_donor.pass1_spec`. Arms so far:
  A = T1K with IPD full-gene reference (`genomic`), B = A plus cohort panel of the donor's fold
  (`genomic_core_hprc_asian`), built per fold in `$ZEN/pangenome-reference-cv-20261014` (FOLDS.json there).
- Earlier stages, with docs worth reading (local, under `scientist/evidence/`):
  `mosaic-drb1-20261113/DESIGN.md`, `results/DRB1_TEST_RESULT.md`, `graphtest/GRAPH_TEST.md`;
  `mosaic-typer-20261113/DESIGN.md`, `results/YIELD.md`; `confirmatory-test-20261113/PREREGISTRATION.md`,
  `CONFIRMATORY_RESULT.md`, `code/`.
- Populations (from FOLDS.json `donors`): Middle East = Arab (APR*, 52) + Saudi (KSA*, 9), simulated reads only;
  East Asian = CDX, CHB, CHS, KHV, JPT (50), South Asian = BEB, GIH, ITU, PJL, STU (36), real reads.

## Deliverable of every step
`REPORT.md` in your local stage dir: what was run (job ids), counts of completed/failed tasks, the result
tables, deviations from the brief, and open problems. Numbers must be recounted from the row-level output
files, not from memory. State plainly what you did not verify. Your final message to the supervisor: the
key numbers, the path of REPORT.md, and anything that needs a decision. Keep it factual.
