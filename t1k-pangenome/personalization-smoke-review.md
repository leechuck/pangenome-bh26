# Synthetic personalization integration result

IBEX job 52043166 completed. Two training HLA-A paths generated 2,004 paired
fragments (4,008 reads); all 4,008 produced mapped alignment records. vg used
3,415 informative k-mers and sampled eight haplotypes in one chain, retaining
one reference path. This demonstrates execution of k-mer counting, graph
reduction, indexing and paired-read mapping, not correct HLA typing.

The log warns that its k-mer coverage estimate is unreliable. The graph contains
only one short locus; automatic coverage estimation needs further investigation
before using selection on real samples. High synthetic read coverage alone does
not establish a reliable estimate. This result must not be described as evidence
of selection accuracy, calibrated confidence or improvement over T1K.

The previous attempt (52043073, personalization-v1) failed because KMC was absent.
Its failed manifest and logs are preserved on IBEX. The monitor now tracks the
successful v2 attempt; its command, tool hashes and read hashes are recorded in
personalization-smoke-result.json.
