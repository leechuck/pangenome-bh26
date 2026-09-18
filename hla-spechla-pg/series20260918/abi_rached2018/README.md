# Abi-Rached / Paganini 2018: secondary concordance reference

Source: [IGSR HLA collection](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HLA_types/),
[paper](https://doi.org/10.1371/journal.pone.0206512), and
[original S2 table](https://doi.org/10.1371/journal.pone.0206512.s010).
Original downloads and checksums are retained. Rebuild the audit with
`python3 hla-spechla-pg/series20260918/audit_abi_rached.py` from the repository root.

These are PolyPheMe v1.2 exome calls using IMGT 3.28 for five loci (A, B, C,
DQB1, DRB1). Selected discrepancies were investigated, including new laboratory
typing for 21 individuals. The entire set is not independent laboratory truth.
Use it for secondary two-field concordance; keep Gourraud experimental accuracy
and assembly/four-field validation separate. Report overlapping donors once per
analysis and never add the denominators from different truth sources together.

The downloaded table contains 2,693 unique donors; 996 overlap our Gourraud file
and 40 overlap the matched 64-donor experiment. The independently eligible set has
2,191 donors: AFR 563, AMR 268, EAS 435, EUR 473, SAS 452. Exclusions: 109 graph
donors, 268 other members of graph families, 125 without verified pedigree records.
These counts describe metadata eligibility, not completed or submitted runs.
Further sampling must use metadata, never prediction correctness.

Calls retain raw source values and expand slash shorthand and space-separated
allele groups. Both allele copies must be scored as an unordered diploid pair;
prediction ambiguity must not gain credit merely by containing the expected call.
Missing reference calls are excluded from truth eligibility; missing predictions
on eligible truth count as failures. Report ambiguity-compatible and strictly
resolved comparisons separately. This table supplies no four-field truth.

There are 620 missing allele slots, 132 ambiguous slots, 110 asterisk-marked slots,
and two malformed numeric values. Asterisk calls remain unresolved until the
source annotation is established. Both HG00592 DQB1 values are decimal fractions
in the FTP file and numeric cells in S2_Table.xlsx; do not guess intended alleles.

Graph-family exclusion uses canonical donor metadata, including alternate assembly
aliases, and the 1000G pedigree. The matched 64 donors require their separate
held-out graphs, even when reference labels exist here. Keep all labels out of
reference construction and parameter tuning.
