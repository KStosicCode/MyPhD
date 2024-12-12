#!/bin/bash

#directories
ref="/home/ulb/gastrexp/kstosic/hg38/Homo_sapiens_assembly38.fasta" 
known_sites="/home/ulb/gastrexp/kstosic/known_sites" 
aligned_reads="/home/ulb/gastrexp/kstosic/Aligned_Reads_cfDNA"
results="/home/ulb/gastrexp/kstosic/results"
germlineresource_and_pon="/home/ulb/gastrexp/kstosic/germlineresource_and_pon"



# Base quality score recalibration
echo "Base quality score recalibration"

# 1. Build the model on each library
gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO1.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO1_recal_report.table \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites ${known-sites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites ${known-sites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz

gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO2.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO2_recal_report.table \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites ${known-sites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites ${known-sites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz

gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO3.bam \
  -R ${ref} 
  -O ${aligned_reads}/SKO3_recal_report.table \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites ${known-sites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites ${known-sites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz

gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO4.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO4_recal_report.table \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites ${known-sites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites ${known-sites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz

gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO5.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO5_recal_report.table \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites ${known-sites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites ${known-sites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz

gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO6.bam \
  -R ${ref} 
  -O ${aligned_reads}/SKO6_recal_report.table \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/Homo_sapiens_assembly38.known_indels.vcf.gz \
  --known-sites ${known-sites}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
  --known-sites ${known-sites}/1000G_phase1.snps.high_confidence.hg38.vcf.gz



# 2. Apply the model to each sample to adjust the base quality scores
gatk ApplyBQSR \
  -I ${aligned_reads}/SKO1.bam \
  -R ${ref} \
  --bqsr-recal-file ${aligned_reads}/SKO1_recal_report.table \
  -O ${aligned_reads}/SKO1_bqsr.bam

gatk ApplyBQSR \
  -I ${aligned_reads}/SKO2.bam \
  -R ${ref} \
  --bqsr-recal-file ${aligned_reads}/SKO2_recal_report.table \
  -O ${aligned_reads}/SKO2_bqsr.bam

gatk ApplyBQSR \
  -I ${aligned_reads}/SKO3.bam \
  -R ${ref} \
  --bqsr-recal-file ${aligned_reads}/SKO3_recal_report.table \
  -O ${aligned_reads}/SKO3_bqsr.bam

gatk ApplyBQSR \
  -I ${aligned_reads}/SKO4.bam \
  -R ${ref} \
  --bqsr-recal-file ${aligned_reads}/SKO4_recal_report.table \
  -O ${aligned_reads}/SKO4_bqsr.bam

gatk ApplyBQSR \
  -I ${aligned_reads}/SKO5.bam \
  -R ${ref} \
  --bqsr-recal-file ${aligned_reads}/SKO5_recal_report.table \
  -O ${aligned_reads}/SKO5_bqsr.bam

gatk ApplyBQSR \
  -I ${aligned_reads}/SKO6.bam \
  -R ${ref} \
  --bqsr-recal-file ${aligned_reads}/SKO6_recal_report.table \
  -O ${aligned_reads}/SKO6_bqsr.bam



# Call variants - gatk Mutect2
echo "Call variants - gatk Mutect2"

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO1.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO1_f1r2.tar.gz \
  --stats ${results}/SKO1.raw_variants.vcf.stat \
  -O ${results}/SKO1_raw_variants.vcf 

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO2.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO2_f1r2.tar.gz \
  --stats ${results}/SKO2.raw_variants.vcf.stats \
  -O ${results}/SKO2_raw_variants.vcf

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO3.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO3_f1r2.tar.gz \
  --stats ${results}/SKO3.raw_variants.vcf.stats \
  -O ${results}/SKO3_raw_variants.vcf 
 

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO4.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO4_f1r2.tar.gz \
  -O ${results}/SKO4_raw_variants.vcf \
  --stats ${results}/SKO4.raw_variants.vcf.stats

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO5.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO5_f1r2.tar.gz \
  --stats ${results}/SKO5.raw_variants.vcf.stats \
  -O ${results}/SKO5_raw_variants.vcf 

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO6.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO6_f1r2.tar.gz \
  --stats ${results}/SKO6.raw_variants.vcf.stats \
  -O ${results}/SKO6_raw_variants.vcf 



# Merge the stats files 
echo "Merge the stats"

gatk MergeMutectStats \
  --stats ${results}/SKO1.raw_variants.vcf.stats
  --stats ${results}/SKO2.raw_variants.vcf.stats
  --stats ${results}/SKO3.raw_variants.vcf.stats
  --stats ${results}/SKO4.raw_variants.vcf.stats
  --stats ${results}/SKO5.raw_variants.vcf.stats
  --stats ${results}/SKO6.raw_variants.vcf.stats
  -O ${results}/merged.stats



# Passing the raw data to LearnReadOrientationModel
echo "Passing the raw data to LearnReadOrientationModel"

gatk LearnReadOrientationModel \
  -I ${results}/SKO1_f1r2.tar.gz \
  -I ${results}/SKO2_f1r2.tar.gz \
  -I ${results}/SKO3_f1r2.tar.gz \
  -I ${results}/SKO4_f1r2.tar.gz \
  -I ${results}/SKO5_f1r2.tar.gz \
  -I ${results}/SKO6_f1r2.tar.gz \
  -O ${results}/reads-orientation-model.tar.gz



# Generating the Intervals List from the VCF file of common germline variants
echo "Generating the Intervals List from the VCF file of common germline variants"

# 1. Selecting Variants from gnomad VCF file

gatk SelectVariants \
  -V ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --select-type-to-include SNP \
  -O ${results}/common_variants_selected.vcf.gz



# 2. Applying selected variants to generate the intervals list

gatk ScatterIntervalsByNs \
  -R ${ref} \
  -V ${results}/common_variants_selected.vcf.gz \
  -O ${results}/common_variants_intervals.list



# Estimating contamination levels
echo "Estimating contamination levels"

gatk GetPileupSummaries \
  -I ${aligned_reads}/SKO1.bam \
  -V ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  -L ${results}/common_variants_intervals.list \
  -O ${results}/SKO1_pileups.table

gatk GetPileupSummaries \
  -I ${aligned_reads}/SKO2.bam \
  -V ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  -L ${results}/common_variants_intervals.list \
  -O ${results}/SKO2_pileups.table

gatk GetPileupSummaries \
  -I ${aligned_reads}/SKO3.bam \
  -V ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  -L ${results}/common_variants_intervals.list \
  -O ${results}/SKO3_pileups.table

gatk GetPileupSummaries \
  -I ${aligned_reads}/SKO4.bam \
  -V ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  -L ${results}/common_variants_intervals.list \
  -O ${results}/SKO4_pileups.table

gatk GetPileupSummaries \
  -I ${aligned_reads}/SKO5.bam \
  -V ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  -L ${results}/common_variants_intervals.list \
  -O ${results}/SKO5_pileups.table

gatk GetPileupSummaries \
  -I ${aligned_reads}/SKO6.bam \
  -V ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  -L ${results}/common_variants_intervals.list \
  -O ${results}/SKO6_pileups.table



# Quantifying contamination
echo "Quantifying contamination"

gatk CalculateContamination \
  -I ${results}/SKO1_pileups.table \
  -O ${results}/SKO1_calculatedcontamination.table

gatk CalculateContamination \
  -I ${results}/SKO2_pileups.table \
  -O ${results}/SKO2_calculatedcontamination.table

gatk CalculateContamination \
  -I ${results}/SKO3_pileups.table \
  -O ${results}/SKO3_calculatedcontamination.table

gatk CalculateContamination \
  -I ${results}/SKO4_pileups.table \
  -O ${results}/SKO4_calculatedcontamination.table

gatk CalculateContamination \
  -I ${results}/SKO5_pileups.table \
  -O ${results}/SKO5_calculatedcontamination.table

gatk CalculateContamination \
  -I ${results}/SKO6_pileups.table \
  -O ${results}/SKO6_calculatedcontamination.table



# Filtering out artifacts using the generated models and stats
echo "Filtering out artifacts using the generated models and stats"

gatk FilterMutectCalls \
  -V ${results}/SKO1_raw_variants.vcf \
  -R ${ref} \
  --contamination-table ${results}/SKO1_calculatedcontamination.table \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO1_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO2_raw_variants.vcf \
  -R ${ref} \
  --contamination-table ${results}/SKO2_calculatedcontamination.table \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO2_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO3_raw_variants.vcf \
  -R ${ref} \
  --contamination-table ${results}/SKO3_calculatedcontamination.table \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO3_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO4_raw_variants.vcf \
  -R ${ref} \
  --contamination-table ${results}/SKO4_calculatedcontamination.table \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO4_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO5_raw_variants.vcf \
  -R ${ref} \
  --contamination-table ${results}/SKO5_calculatedcontamination.table \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO5_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO6_raw_variants.vcf \
  -R ${ref} \
  --contamination-table ${results}/SKO6_calculatedcontamination.table \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO6_filtered_variants.vcf



# Downloading the Funcotator

gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --hg38 --extract-after-download



# Annotate Variants - GATK4 Funcotator
echo "Annotate Variants - GATK4 Funcotator"

gatk Funcotator \
  --variant ${results}/SKO1_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.7.20200521g \
  --output ${results}/SKO1_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO2_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.7$
  --output ${results}/SKO2_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO3_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.7$
  --output ${results}/SKO3_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO4_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.7$
  --output ${results}/SKO4_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO5_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.7$
  --output ${results}/SKO5_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO6_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.7$
  --output ${results}/SKO6_filtered_variants_funcotated.maf \
  --output-file-format MAF

  
