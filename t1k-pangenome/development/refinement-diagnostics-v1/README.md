# Changed-call score diagnostics

Both diagnostic jobs completed and reproduced the recorded score changes to
floating-point tolerance. They use only exposed development results.

| Case | Total score favouring wrong pair | Overlapping gene span | Outside gene span |
|---|---:|---:|---:|
| HG01530 A, additive unsampled | +78.585 | -44.008 | +122.592 |
| HG02132 DPB1, additive sampled | +22.666 | +22.666 | 0.000 |

For HLA-A, surrounding context overwhelms gene-overlapping evidence favouring
the correct native pair. Preventing context-only support from deciding allele
labels is therefore a concrete next revision. The exon envelope includes introns
and excludes terminal UTRs; this diagnostic is not a variant-specific likelihood.
It classifies a fragment as overlapping if any relevant placement of either mate
overlaps that envelope. A revised caller will need a stricter rule that excludes
flanking contributions from the full paired score, not just this diagnostic flag.

For DPB1, the error is internal and remains unexplained. Sampled versus unsampled
mapping and support at the distinguishing variants need inspection. Lack of
support for both native alleles in one fragment is not itself evidence of bias:
a fragment normally originates from a single haplotype.

Initial diagnostic jobs failed before reading data because the container does
not forward SLURM_JOB_ID. The corrected jobs use its forwarded
SLURM_CPUS_PER_TASK allocation marker, matching the inference adapter. Both
attempts and code hashes are retained in graph-refinement-diagnostic-launch.json.
