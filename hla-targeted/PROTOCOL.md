# Targeted RCCX/C4 experiment

Status: development in progress; fresh-donor evaluation has not started.

## Question and endpoints

Can targeted C4 diagnostic probes improve marginal copy-number typing and structural-pair imputation beyond the prior sketch-based assay? Does module-start context add information? Neither graph storage nor a uniquely ranked reference pair establishes physical phase.

Report exact total C4, A, B, noncanonical diagnostic-motif, long and short dosage; full RCCX annotation-signature pair accuracy; coverage of the true signature pair by the dosage-compatible candidate set; and the size of that set. Report all-attempted denominators and abstentions. Joint AL/AS/BL/BS dosage, order and chromosome assignment must not be inferred directly from unphased A/B and long/short marginal counts.

## Data split

The previously examined 106 donors are development data, including all earlier pilot and validation samples. The 24 new donors were selected from public-WGS availability before inspecting their read-based predictions. All recorded families of these 24 are excluded from assay discovery and candidate paths. Their assemblies were nevertheless present in the existing graph construction, so this is new-donor read validation within the existing panel, not external-cohort validation. Family labels do not prove absence of all cryptic relatedness.

Development results guide the exploratory assay; they are not estimates of prospective performance. Freeze the code, probes, reference profiles, calibration and sample selection before inspecting fresh-donor predictions. Preserve the original freeze and identify any subsequent corrective analyses separately.

## Evidence and limitations

A/B probes span the five published diagnostic positions. Noncanonical motifs are retained as a separate category and do not establish novel named HLA or functional C4 alleles. Long/short probes span HERV insertion/deletion junctions. Module-start contexts are centred on WHR1-family annotations and are not proven SV breakpoints. All targeted probes were audited against the full training MHC sequences; whole-genome specificity has not been established.

Reads are public WGS recruited to the MHC using existing linear-reference alignments. Count exact canonical Q20 31-mers once per fragment and count each diagnostic group once per fragment. This is conditional on MHC recruitment and excludes reads absent from those intervals. It does not test unbiased whole-WGS recruitment.

Fit relative A/B and long/short fragment detection efficiencies using mixed-copy development donors. Retain the earlier pilot-derived total-C4 dosage correction. Evidence thresholds and nearest-integer calls are exploratory heuristics, not calibrated confidence probabilities. Explicitly report noncalls and discrepancies between the two long-junction measurements.

## Comparators and structural inference

Run published C4Investigator at the recorded upstream commit, with its supplied references and algorithm. Give it the same available MHC-recruited read source; retain its native calls. Compare calibrated ordinary reference depth using the earlier pilot-derived factors.

Compare the old marker sketch constrained by total copy number with the targeted marker model constrained by supported marginal dosage. For an ablation, omit module-context probes from the targeted score. Retain all candidate structural pairs compatible with accepted dosage calls; score-based ranking within that set is imputation. Missing reference structures and inaccurate marginal dosage can exclude the true pair, so compatibility is not a validated confidence set or unknown-structure detector.

Assembly annotations are the comparison labels, not independent clinical truth. Independently audit diagnostic motifs in any new-donor disagreements; report label corrections transparently and preserve original comparisons. UK Biobank phenotype and disease validation remains a post-hackathon study.
