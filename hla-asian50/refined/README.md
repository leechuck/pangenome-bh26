# Callable-site repair and named HLA rerun

See [REPORT.md](REPORT.md) for the completed results and [DESIGN.md](DESIGN.md) for frozen comparison rules and limitations. Results
are separated into variant concordance, experimental HLA accuracy, and agreement
with assembly-derived HLA labels.

Remote working directory: `/home/leechuck/hla/codex-asian50/run/refined` on DDBJ.
The parent directory holds the original run and frozen donor/read inputs.
All compute runs through Slurm; the supplied batch scripts encode the allocation.

## Reproduce

1. Link `source`, `graph`, `reads`, and `accuracy` to the corresponding original
   run directories. Populate `hla_source` from the eight archived gene-region
   FASTAs, region tables, and exact-CDS sequence catalogue. Preserve
   `experimental_HLA.txt` and `native_HLA_calls.tar.gz` from the original run.
   The eight gene FASTAs and region tables are under `hla-analysis/source/`; the
   sequence catalogue is `hla-analysis/results/sequence_catalogue.tsv`.
2. Submit `build.sbatch` to build all ten repaired variant panels and indices.
   Alternatively, run the panel builder in a batch allocation, then submit
   `index.sbatch` as array `0-9` after the builder completes.
3. Submit `pangenie.sbatch` as array `0-39` after indexing, then `score.sbatch`.
   The scorer asserts preservation of old sites and exact reproduction of the
   old metrics using an independent implementation.
4. The initial whole-locus PanGenie HLA attempt is retained for reproducibility:
   `hla_build.sbatch`, `hla_genotype.sbatch` (array `0-39`), `score_hla.sbatch`.
   It is a failed representation experiment, not the main named-HLA method.
5. Submit `locityper_setup.sbatch`, then `locityper_genotype.sbatch` as array
   `0-39`, then `score_locityper.sbatch`. The latter uses the frozen eligible
   truth rows in `results/hla_scores.tsv`; it does not select test labels to
   resolve prediction ambiguity.
6. Submit `archive.sbatch` after variant/initial-HLA scoring and
   `archive_locityper.sbatch` after named-HLA scoring. Download `results` and
   run `python3 summarize.py` locally (requires NumPy).

Use `afterok` dependencies for each dependent step and concurrency limits that
respect the shared partition. `source/jobs.json` records the actual run IDs,
including superseded jobs. Local tests: `python3 -m unittest discover -s .`.

The repair uses native PanGenie phased-missing support; it does not fill missing
haplotypes with reference alleles or change GBZ topology. Locityper types intact
source-assembly locus haplotypes and maps them to frozen numeric two-field HLA
labels. It does not infer HLA names directly from unphased variant calls.
