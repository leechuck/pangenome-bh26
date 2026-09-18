# Exploratory diagnosis after validation unblinding

This analysis explains score changes in the completed frozen validation. It does not replace that result or establish a new validated improvement. All seven call-changing decisions across both graph panels were examined.

Every diagnostic reproduced the frozen native-to-proposed pair score gain (absolute tolerance 1e-6) and fragment count before writing results. The comparison then sets path effective lengths equal while preserving the same mapped placements, selected fragments, proposed pair and native pair. This is a fixed-pair ablation, not a new genotype search.

| Donor | Panel | Gene | Frozen gain | Equal-length gain | Accuracy change |
|---|---|---|---:|---:|---|
| HG01891 | hprc | DQA1 | 37.828 | 38.031 | rescue |
| HG02080 | hprc_asian | C | 15.202 | -27.454 | loss |
| HG02257 | hprc_asian | A | 19.711 | 19.610 | loss |
| HG02258 | hprc_asian | DQA1 | 16.615 | 16.574 | rescue |
| HG02615 | hprc | C | 15.770 | -19.610 | ineligible truth |
| HG03270 | hprc_asian | DQA1 | 28.190 | 28.188 | ineligible truth |
| HG03704 | hprc_asian | DQA1 | 10.137 | -15.688 | loss |

For HG02080 HLA-C, fragments with equal pair alignment emissions contribute +42.143 through path-length effects, overwhelming −26.941 from alignment-discriminating fragments. For HG03704 DQA1, the corresponding contributions are +25.825 and −15.688. Thus path normalization reverses the evidence in two of the three accuracy losses. The HLA-A loss (HG02257) remains and needs a different explanation.

The frozen pair-specific fragment counter also counts path-length-induced differences, so its minimum-20 gate is not equivalent to 20 sequence-discriminating fragments. HG02080 has only 13 fragments with differing alignment pair emissions; HG03704 has eight. HG01891 has 13 as well, so removing this bias may also discard a correct HPRC-only refinement. That trade-off must be measured across the full cohort.

In the JSON outputs, `length_only_pair_difference` includes all fragments whose equal-length pair scores are identical, including fragments with zero contribution under both models. `native_allele_absent` means at least one native allele lacks a placement for that fragment; it does not mean both native alleles are absent from the graph.

The first diagnostic submission (52059427) stopped at the exact-replay assertion because the diagnostic, not inference, transcribed the noise term incorrectly. The corrected expression was tested against the frozen scorer, including zero emissions; all seven tasks in 52059510 then passed. Original validation calls and scores remain unchanged.

An isolated alignment-only candidate is implemented in `graph_pair_refinement_alignment_only.py` and `run_graph_pair_alignment_only.py`. The original frozen files are untouched. It retains coarse genotypes, fallback rules, locus competition and thresholds, while removing full-path length normalization from the selected internal-read emissions. Its manifest records both the new model and the frozen base-model hashes.

The planned ablation covers all eight development donors and all 84 exposed validation donors, with both graph panels (184 cases). It is explicitly exploratory. A gain on these exposed donors would still require a fresh, predeclared independent evaluation.
