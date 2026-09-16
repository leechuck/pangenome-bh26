# Dedicated comparator configuration and calibration check

C4Investigator runs from upstream commit `a3bd9eb860bf9e7bb712554c61d8775e24eeeabb` with its supplied references and unmodified algorithm. The upstream script explicitly sets `run.mode <- 'WGS'`. Our invocation uses the documented FASTQ directory, output directory, filename-pattern and thread arguments; it does not accidentally select targeted-sequencing mode.

The native summary reports `C4_copy` as the sum of its inferred C4A and C4B counts. The detailed output also provides a continuous WGS total estimate (`wgs_c4_copy`). Both are retained. These are distinct outputs; the supplementary calibrated result is derived from the continuous estimate, not a replacement for the native algorithm.

Across the original 18 pilot donors, native total-copy calls match 3/18 assembly labels. The predeclared scalar calibration has factor **1.1596045197740117** and matches 18/18 pilot totals after correction. This is a development fit, not held-out performance. The recipe was recorded before inspecting held-out read predictions; it uses only the original pilot samples and no validation outcomes.

The native component calls can inherit total-dosage bias. Therefore differences between native C4Investigator component calls and the calibrated targeted assay must not be interpreted as general superiority of the targeted method. The final report includes the separately calibrated total-copy comparison as a fairness diagnostic. Inputs are MHC-recruited WGS reads under the recorded references and environment; this experiment does not establish performance on unrestricted whole-WGS inputs.

See `results/comparator_development_review.tsv`, `source/comparator_calibration_plan.json`, the original upstream source, and `c4investigator.sbatch` for evidence and commands.
