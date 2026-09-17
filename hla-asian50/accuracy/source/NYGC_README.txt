1000 Genomes Processing README
This README contains information relating to data associated with the 1000 Genomes
resequencing done at New York Genome Center.

Alignment, post-processing and variant calling
Alignment and post-processing are performed exactly as outlined by the Center for
Common Disease Genomics project: https://github.com/CCDG/PipelineStandardization/blob/master/PipelineStandard.md .
Programs and reference data
The data was aligned to the reference genome using the following programs and reference
datasets:
1.
2.
3.
4.
5.

BWA-MEM
Samtools-1.3.1
Picard-2.4.1
GATK-3.5-0
Resource files
– All the resource files used in the analysis can be obtained here:
https://console.cloud.google.com/storage/browser/genomics-publicdata/resources/broad/hg38/v0/ .

Reference genome: GRCh38 with alternative sequences, plus decoys and HLA
The reference genome that the data was aligned to can be obtained here:
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome
/GRCh38_full_analysis_set_plus_decoy_hla.fa

Command lines
1.

Alignment at lane level
bwa mem
-Y
-K 100000000
-t 16
-R $rg_string
$reference_fasta
$fastq_file(1)
$fastq_file(2) | samtools view -Shb -o $bam_aligned -

2.

Fix mate information in the BAM
java $jvm_args -jar picard.jar
FixMateInformation
MAX_RECORDS_IN_RAM=2000000

VALIDATION_STRINGENCY=SILENT
ADD_MATE_CIGAR=True
ASSUME_SORTED=true
I=$bam_aligned
O=$bam_fixmate

3.

Merging lane level bam files to Sample level bam files
java $jvm_args -jar picard.jar
MergeSamFiles
USE_THREADING=true
MAX_RECORDS_IN_RAM=2000000
VALIDATION_STRINGENCY=SILENT
SORT_ORDER=queryname
INPUT=$bam_fixmate_1
INPUT=$bam_fixmate_2
OUTPUT=$bam_merged

4.

Mark duplicates and coordinate sort BAM
java $jvm_args -jar picard.jar
MarkDuplicates
MAX_RECORDS_IN_RAM=2000000
VALIDATION_STRINGENCY=SILENT
M=$dedup_metrics
I=$bam_merged
O=$bam_dedup
java $jvm_args -jar picard.jar
SortSam
MAX_RECORDS_IN_RAM=2000000
VALIDATION_STRINGENCY=SILENT
SORT_ORDER=coordinate
CREATE_INDEX=true
I=$bam_dedup
O=$bam_sorted

5.

Recalibrate base quality scores using known SNPs
java $jvm_args -jar GenomeAnalysisTK.jar
-T BaseRecalibrator
-downsample_to_fraction 0.1
-nct 4
--preserve_qscores_less_than 6
-L $autosomes
-R $reference_fasta
-o $recal_data.table
-I $bam_sorted
-knownSites $known_snps_from_dbSNP138
-knownSites $known_indels
-knownSites $known_indels_from_mills_1000genomes

java $jvm_args -jar GenomeAnalysisTK.jar
-T PrintReads
-nct 4
--disable_indel_quals
--preserve_qscores_less_than 6
-SQQ 10
-SQQ 20
-SQQ 30
-rf BadCigar
-R $reference_fasta
-o $bam_recalibrated
-I $bam_sorted
-BQSR $recal_data.table

6.

Creating CRAM files
samtools view
-C
-T $reference_fasta
-o $cram
$bam_recalibrated
samtools index $cram

7.

Raw variant calls using HaplotypeCaller on single sample

For variant discovery we used HaplotypeCaller in GVCF mode (Poplin et al., 2017) with sexdependent ploidy settings on chromosome X and Y. Specifically, variant discovery on chrX
was performed using diploid settings in females, diploid settings on PAR regions in males,
and haploid settings on non-PAR regions in males. Variant discovery on chrY was
performed with haploid settings in males and skipped in females.
java $jvm_args -jar GenomeAnalysisTK.jar
-T HaplotypeCaller
--genotyping_mode DISCOVERY
-A AlleleBalanceBySample
-A DepthPerAlleleBySample
-A DepthPerSampleHC
-A InbreedingCoeff
-A MappingQualityZeroBySample
-A StrandBiasBySample
-A Coverage
-A FisherStrand
-A HaplotypeScore
-A MappingQualityRankSumTest
-A MappingQualityZero
-A QualByDepth
-A RMSMappingQuality
-A ReadPosRankSumTest
-A VariantType
-l INFO

--emitRefConfidence GVCF
-rf BadCigar
--variant_index_parameter 128000
--variant_index_type LINEAR
-R $reference_fasta
-nct 1
-L $interval
-ploidy $ploidy
-I $bam_recalibrated
-o $gvcf

8.

Jointly recalibrate Genotype Quality score of all samples
java $jvm_args -jar GenomeAnalysisTK.jar
-T GenotypeGVCFs
-R $reference_fasta
-nt 5
--disable_auto_index_creation_and_locking_when_reading_rods
--variant $gvcf
-o $vcf_genotyped

9.

Variant Quality Score Recalibration (VQSR) to assign FILTER status
java $jvm_args -jar GenomeAnalysisTK.jar
-T VariantRecalibrator
-R $reference_fasta
-nt 5
-input $vcf_genotyped
-mode SNP
-recalFile $vqsr_snp.recal
-tranchesFile $vqsr_snp.tranches
-rscriptFile $vqsr_snp_plots.R
-resource:hapmap,known=false,training=true,truth=true,prior=15.0

$hapmap
-resource:omni,known=false,training=true,truth=true,prior=12.0
$kg_omni
-resource:1000G,known=false,training=true,truth=false,prior=10.0
$kg_snps
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0
$dbsnp
-an QD
-an MQ
-an FS
-an MQRankSum
-an ReadPosRankSum
-an SOR
-an DP
-tranche 100.0
-tranche 99.8
-tranche 99.6
-tranche 99.4
-tranche 99.2

-tranche 99.0
-tranche 95.0
-tranche 90.0
java $jvm_args -jar GenomeAnalysisTK.jar
-T VariantRecalibrator
-R $reference_fasta
-nt 5
-input $vcf_genotyped
-mode INDEL
-recalFile $recalibrate_indel.recal
-tranchesFile $recalibrate_indel.tranches
-rscriptFile $recalibrate_indel_plots.R
-resource:mills,known=true,training=true,truth=true,prior=12.0
$kg_mills
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0
$dbsnp
-an QD
-an FS
-an ReadPosRankSum
-an MQRankSum
-an SOR
-an DP
-tranche 100.0
-tranche 99.0
-tranche 95.0
-tranche 92.0
-tranche 90.0
--maxGaussians 4
java $jvm_args -jar GenomeAnalysisTK.jar
-T ApplyRecalibration
-R $reference_fasta
-nt 5
-input $vcf_genotyped
-mode SNP
--ts_filter_level 99.80
-recalFile $recalibrate_SNP.recal
-tranchesFile $recalibrate_SNP.tranches
-o $vcf_recalibrated_snp
java $jvm_args -jar GenomeAnalysisTK.jar
-T ApplyRecalibration
-R $reference_fasta
-nt 5
-input $vcf_recalibrated_snp
-mode INDEL
--ts_filter_level 99.0
-recalFile $recalibrate_INDEL.recal

-tranchesFile $recalibrate_INDEL.tranches
-o $vcf_recalibrated_snp_indel

Definitions of delivered files
1.
2.

3.
4.

5.

6.

recalibrated_variants.vcf.gz[.tbi]
– All variants in Variant Call Format (VCF) file along with index.
recalibrated_variants.annotated.vcf.gz[.tbi]
– Normalized VCF stripped of genotype calls and annotated using snpEff and
BCFTools.
recalibrated_variants.annotated.txt
– Variant annotations in a tab-delimited file.
recalibrated_variants.annotated.coding.txt
– All annotated variants with HIGH/MODERATE impact in a tab-delimited file.
– High impact - The variant is assumed to have high (disruptive) impact in the
protein, probably causing protein truncation, loss of function or triggering
nonsense mediated decay. e.g. stop_gained, frameshift_variant.
– Moderate impact - A non-disruptive variant that might change protein
effectiveness. e.g. missense_variant, inframe_deletion
recalibrated_variants.annotated.coding_rare.txt
– All HIGH/MODERATE annotated variants with less than 5% allele frequency in
1000genomes and ExAC in a tab-delimited file.
recalibrated_variants.annotated.clinical.txt
– All low frequency HIGH/MODERATE annotated variants with possible clinical
impact from ClinVar in a tab-delimited file.

Acknowledgements
The following cell lines/DNA samples were obtained from the NIGMS Human Genetic Cell
Repository at the Coriell Institute for Medical Research: [NA06984, NA06985, NA06986,
NA06989, NA06991, NA06993, NA06994, NA06995, NA06997, NA07000, NA07014,
NA07019, NA07022, NA07029, NA07031, NA07034, NA07037, NA07045, NA07048,
NA07051, NA07055, NA07056, NA07340, NA07345, NA07346, NA07347, NA07348,
NA07349, NA07357, NA07435, NA10830, NA10831, NA10835, NA10836, NA10837,
NA10838, NA10839, NA10840, NA10842, NA10843, NA10845, NA10846, NA10847,
NA10850, NA10851, NA10852, NA10853, NA10854, NA10855, NA10856, NA10857,
NA10859, NA10860, NA10861, NA10863, NA10864, NA10865, NA11829, NA11830,
NA11831, NA11832, NA11839, NA11840, NA11843, NA11881, NA11882, NA11891,
NA11892, NA11893, NA11894, NA11917, NA11918, NA11919, NA11920, NA11930,
NA11931, NA11932, NA11933, NA11992, NA11993, NA11994, NA11995, NA12003,
NA12004, NA12005, NA12006, NA12043, NA12044, NA12045, NA12046, NA12056,
NA12057, NA12058, NA12144, NA12145, NA12146, NA12154, NA12155, NA12156,
NA12234, NA12239, NA12248, NA12249, NA12264, NA12272, NA12273, NA12274,
NA12275, NA12282, NA12283, NA12286, NA12287, NA12329, NA12335, NA12336,
NA12340, NA12341, NA12342, NA12343, NA12344, NA12347, NA12348, NA12375,

NA12376, NA12383, NA12386, NA12399, NA12400, NA12413, NA12414, NA12485,
NA12489, NA12546, NA12707, NA12708, NA12716, NA12717, NA12718, NA12739,
NA12740, NA12748, NA12749, NA12750, NA12751, NA12752, NA12753, NA12760,
NA12761, NA12762, NA12763, NA12766, NA12767, NA12775, NA12776, NA12777,
NA12778, NA12801, NA12802, NA12812, NA12813, NA12814, NA12815, NA12817,
NA12818, NA12827, NA12828, NA12829, NA12830, NA12832, NA12842, NA12843,
NA12864, NA12865, NA12872, NA12873, NA12874, NA12875, NA12877, NA12878,
NA12889, NA12890, NA12891, NA12892]. The data were generated at the New York
Genome Center with funds provided by NHGRI Grants 3UM1HG008901-03S1 and
3UM1HG008901-04S2.
If you have any questions please email service@nygenome.org .

