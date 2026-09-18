# Native-locus audit of graph-discriminating fragments

Every frozen changed decision was audited after unblinding. All seven diagnostic
tasks (52059934) completed and reproduced their original pair scores. Fragment
identities, graph placements and native T1K assignments are retained remotely and
in `work/postvalidation-read-audit-v1`; `summary.json` records their hashes and
positive-native-assignment locus counts.

For **HG02257 HLA-A**, 16 fragments favour the incorrect proposed pair and 11
favour the correct native pair. All 16 pro-proposal fragments have positive
native T1K assignments exclusively to **HLA-E**. Ten of the opposing fragments
have native assignments to HLA-A; the remaining one has no native assignment.
The eight-locus graph competition omits HLA-E, allowing paralogous evidence to
be counted as HLA-A-specific support. The fixed-pair score therefore mistakes
missing-locus competition for strong allele evidence.

This is not unique to that sample: in HG02080 HLA-C, the three fragments whose
alignment emissions favour the wrong proposal have exclusive native assignments
to HLA-Y, HLA-F and HLA-J, respectively. The ten opposing fragments have native
HLA-C assignments. Conversely, 29 of the 30 fragments favouring the HG02258 DQA1
rescue have native DQA1 assignments, with one graph-only fragment.

The new isolated candidate vetoes a winning graph locus only when positive
native assignments support exclusively other loci. It preserves graph-only
fragments and does not claim to resolve mixed native-locus evidence. This is a
conservative integration of the already-exported native evidence, not calibrated
joint likelihood, and is being evaluated across all 92 exposed donors and both
panels. Original validation outputs and frozen sources are unchanged.
