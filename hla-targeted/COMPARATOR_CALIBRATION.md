# Secondary C4Investigator total-dosage calibration

The unmodified C4Investigator calls remain in the primary report. This fairness diagnostic gives its continuous total-C4 estimate the same pilot-only scalar-calibration recipe used for the depth and marker baselines.

The recipe was fixed at 2026-09-16T04:24:01.358699+00:00, before inspecting held-out read predictions. The fitted factor is 1.1596, using the original 18 pilot donors only. No held-out outcomes enter the fit. This supplement does not recalibrate A/B or long/short component calls or alter the frozen targeted assay.

Calibrated total-copy accuracy: **24/24** attempted held-out donors. The per-donor estimates and calls are in `results/validation_C4Investigator_calibrated_total.tsv`.

This comparison is conditional on MHC-recruited reads and the recorded references and software environment. Differences from the native caller do not establish its performance on unrestricted whole-WGS inputs.
