# Callable-site repair and named HLA comparison

This rerun keeps the original 40 donors (20 EAS, 20 SAS), five test-family folds,
read recruitment, caller versions and native specialist predictions frozen.
Saudi/Arabian test samples remain deferred. The global graph is still transductive.

## Variant-panel repair

The original complete-case filter was unnecessarily strict: PanGenie 4.2.1
explicitly supports `.|.`, `0|.` and other phased missing genotypes. Preserve
known alleles and missing haplotypes instead of dropping the site. Invalid or
unphased training genotypes, conflicting paths, and non-ACGT haplotype alleles
are conservatively masked as missing; they are never converted to reference.
Only alleles observed in eligible training haplotypes, plus the reference allele,
enter the panel. Held-out-only alternate alleles remain excluded.

Both HPRC-only and full panels are rebuilt with the same rule. The builder asserts
that all old sites survive and that full-panel sites contain HPRC-panel sites.
The scorer checks these properties again on actual output VCFs. This is a repair
of the graph's genotyping-panel representation, not a change to the GBZ topology.

Evaluation uses unphased allele-sequence identity, masking missing/conflicting
assembly truth exactly as before. Fixed universes:

- All eligible original assembly-truth sites, including omissions as failures.
- Original HPRC-callable sites, fixed independently of new predictions.
- Original three-way shared PASS SNVs, retaining the same published linear calls.
- Shared new-panel sites, explicitly a different denominator from the first run.

Report EAS and SAS separately, SNVs and length-changing SV-bearing genotypes
separately, and paired donor bootstrap intervals. The new independent SHA256-based
scorer must reproduce all old all-truth counts exactly. SHA256 represents allele
sequence identity; numeric ALT indices are never compared across VCFs.

## Named HLA endpoint

Create eight whole-locus multiallelic PanGenie loci from the phased assembly
paths used to construct the graph: HLA-A, B, C, DRB1, DQA1, DQB1, DPA1, DPB1.
Use the existing strand-oriented gene +/-2 kb region FASTAs. The first and last
64 bp provide fixed GRCh38 context; the inner region is a whole-locus allele.
This is a new locus-level representation, not a claim that unphased per-variant
calls have already been translated into phased classical HLA types.

The same training donors and test-family exclusions apply. Each arm retains only
its training haplotypes plus the fixed GRCh38 reference. Missing gene sequences
remain missing. Identical allele sequences are deduplicated. Retain sequences
with unresolved labels rather than artificially restricting the prediction space.

For each training allele, attach the exact current CDS-match labels from the
existing frozen sequence catalogue. Translate a predicted allele to a numeric
two-field HLA type only when these labels collapse to a single type. Unknown or
ambiguous labels become no-calls; do not select a label using the test donor.
The GRCh38 reference label is allowed in both arms as a fixed external reference.
Include the remaining MHC reference as background with the eight extracted
regions masked, to retain off-locus reference sequence in k-mer counting. Verify
all eight reference extractions against the MHC reference exactly.

Compare named HLA predictions with the existing T1K 1.0.6 and SpecHLA calls:

- Primary independent endpoint: historical experimental labels, eight EAS donors,
  39 eligible diploid genotypes across five loci. Retain historical ambiguities.
- Secondary endpoint: complete exact-CDS assembly labels, 159 EAS and 160 SAS
  diploid genotypes across eight loci. These are assembly agreement, not independent
  experimental accuracy. All methods use identical eligible donor/locus pairs.

Report exact unordered two-field genotype matches and allele matches. Keep all
no-calls in accuracy denominators. Report label representation coverage separately.
Numeric resolution does not evaluate expression suffixes or higher-field typing.
T1K uses its native quality >0 recommendation, and one reported allele at these
classical loci is interpreted as homozygous (not a copy-number inference).
SpecHLA retains its native IPD-IMGT/HLA 3.38.0 database. Database differences,
historical assay resolution, limited experimental sample size, and lack of an
unseen-graph test prohibit a general claim of superiority over specialist typers.
No tuning against the observed test labels is permitted.

## Whole-locus PanGenie failure and replacement (before Locityper scoring)

The initial whole-locus PanGenie test yielded almost entirely missing genotypes.
The native caller reported UK=0 at seven of eight loci for the inspected donor.
Source inspection identifies the mechanism: `select_kmers` in
`stepwiseuniquekmercomputer.cpp` rejects k-mers shared by multiple alleles.
Most HLA sequence markers are shared by multiple whole-gene haplotypes; the
few whole-haplotype-private or artificial boundary markers do not transfer to
test donors. Preserve this failed test; do not interpret its no-calls as HLA
accuracy of the graph itself or tune the marker length against truth.

Use Locityper 1.7.4 on the same full and HPRC-only training memberships as the
replacement locus genotyper. It uses intact, strand-oriented region sequences,
without the artificial boundary joins in the PanGenie whole-locus VCF.
The existing exact-CDS label rules, truth sets and no-call denominators stay fixed.
Deduplicate identical whole-region sequences and retain all source-haplotype
identities in the allele-label map. No external IPD sequences are added to the
prediction panel; GRCh38 is the same fixed external reference in both arms.

Count canonical 25-mers genome-wide on GRCh38 primary chromosomes (1–22, X, Y, M),
excluding alternate representations/decoys so they are not mistaken for extra
physical copies. Follow the official repeat-count settings (`--lower-count 2`,
`--out-counter-len 2`). Use the identical MHC FASTQs for preprocessing and genotyping.
Because these FASTQs contain recruited MHC reads only, calibrate per-sample depth,
insert sizes and errors within the covered 3 Mb MHC interval ctg1:100000–3100000,
instead of Locityper's default chr17 background. The MHC reference is used for
regional mapping, while the genome-wide repeat counts mask duplicated k-mers.
This regional calibration is a pilot limitation, not a standard whole-WGS
Locityper benchmark. Default genotyping settings and seed 20260917 are fixed.

Official instructions: https://locityper.vercel.app/target,
https://locityper.vercel.app/preproc, https://locityper.vercel.app/genotype.

Before evaluating Locityper predictions, reverse-complement all gene-strand
FASTA sequences at reference-minus-strand loci to genomic-plus orientation.
Assert that the exact GRCh38 interval sequence is present in each target panel.
The initial target-build warnings revealed this orientation mismatch; that build
was cancelled before genotyping. Annotation-defined flank endpoints and genuine
polymorphism can still differ from the reference boundaries. Retain target-build
warnings and report this limitation, along with loci shorter than the suggested
10 kb target length. No boundaries are adjusted using test results.
