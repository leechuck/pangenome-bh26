Last modified: November 18th, 2022.
This directory contains the high-coverage 1000 Genomes Project (1kGP) phased WGS panel
consisting of high-quality SNV and INDEL calls across 3,202 1kGP samples.
Note: this is not the final phased panel. The phased panel included here was later expanded to
include structural variant calls. The final and most comprehensive version of the phased panel
based on the 3,202 high-coverage 1kGP samples is included in the
20220422_3202_phased_SNV_INDEL_SV/ directory.
Filtering criteria applied to the small variant callset prior to phasing:
● FILTER=PASS.
● Genotype missingness < 5%.
● Pass HWE test (i.e. HWE p-value > 1e-10 in at least one of the 5 super-populations).
● Mendelian error rate ≤ 5%.
● Minor allele count (MAC) ≥ 2.
For haplotype phasing of high-quality autosomal SNVs and INDELs we used statistical phasing
with pedigree-based correction, as implemented in the SHAPEIT2-duohmm software (Delaneau
et al., 2011; O'Connell et al., 2014). This approach leverages the inclusion of 602 trios in the
dataset to increase the accuracy of phasing.
SHAPEIT2 does not handle multiallelic variant phasing. To phase both biallelic and multiallelic
variants, we first split the multiallelics into separate rows, while left-aligning and normalizing
INDELs, using bcftools norm tool (Li, 2011). Then we shifted the position of multiallelic variants
(2nd, 3rd, etc ALT alleles) by 1 or more bp (depending on how many ALT alleles there are at a
given position) to ensure a unique start position for all variants, which is required for SHAPEIT2.
We shifted the positions back to the original ones after phasing.
Phasing with SHAPEIT2-duohmm was performed per chromosome using the default settings,
except for the window size "-W" which was increased from 2Mb (default) to 5Mb to account for
increased amounts of shared IBD due to pedigrees being present in the dataset (as
recommended in the SHAPEIT2 manual).
SHAPEIT2-duohmm supports phasing of autosomal variants only. Therefore, to phase variants
on chrX we used statistical phasing as implemented in the Eagle2 software (Loh et al., 2016).
Phasing with Eagle2 was performed using default parameters. No shifting of positions for
multiallelics was needed as Eagle2 supports phasing of variants with the same start site. Calls
on chrX and chrY were generated using sex-dependent ploidy settings.

