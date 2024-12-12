#!/bin/bash

# directories
ref="/home/ulb/gastrexp/kstosic/hg38/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"
known_sites="/home/ulb/gastrexp/kstosic/known_sites"
aligned_reads="/home/ulb/gastrexp/kstosic/Aligned_Reads_cfDNA"
results="/home/ulb/gastrexp/kstosic/results"



# Base quality score recalibration
echo "Base quality score recalibration"

# 1. Build the model on each library
gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO1.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO1_recal_report.table \ 
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf

gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO2.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO2_recal_report.table \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf

gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO3.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO3_recal_report.table \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf
    
gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO4.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO4_recal_report.table \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf
 
gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO5.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO5_recal_report.table \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf

gatk BaseRecalibrator \
  -I ${aligned_reads}/SKO6.bam \
  -R ${ref} \
  -O ${aligned_reads}/SKO6_recal_report.table \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
  --known-sites ${known_sites}/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf \
  --known-sites ${known-sites}/resources_broad_hg38_v0_1000G_phase1.snps.high_confidence.hg38.vcf



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



# Call Variants - gatk haplotype caller
echo "Call Variants - gatk haplotype caller"

gatk HaplotypeCaller \
  -R ${ref} \
  -I ${aligned_reads}/SKO1_bqsr.bam \ 
  -O ${results}/SKO1_raw_variants.vcf

gatk HaplotypeCaller \
  -R ${ref} \
  -I ${aligned_reads}/SKO2_bqsr.bam \
  -O ${results}/SKO2_raw_variants.vcf

gatk HaplotypeCaller \
  -R ${ref} \
  -I ${aligned_reads}/SKO3_bqsr.bam \
  -O ${results}/SKO3_raw_variants.vcf

gatk HaplotypeCaller \
  -R ${ref} \
  -I ${aligned_reads}/SKO4_bqsr.bam \
  -O ${results}/SKO4_raw_variants.vcf

gatk HaplotypeCaller \
  -R ${ref} \
  -I ${aligned_reads}/SKO5_bqsr.bam \
  -O ${results}/SKO5_raw_variants.vcf

gatk HaplotypeCaller \
  -R ${ref} \
  -I ${aligned_reads}/SKO6_bqsr.bam \
  -O ${results}/SKO6_raw_variants.vcf



# Extract and split separately SNPs & INDELs
echo "Extract and split separately SNPs & INDELs"

gatk SelectVariants -R ${ref} -V ${results}/SKO1_raw_variants.vcf --select-type SNP -O ${results}/SKO1_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/SKO1_raw_variants.vcf --select-type INDEL -O ${results}/SKO1_raw_indels.vcf

gatk SelectVariants -R ${ref} -V ${results}/SKO2_raw_variants.vcf --select-type SNP -O ${results}/SKO2_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/SKO2_raw_variants.vcf --select-type INDEL -O ${results}/SKO2_raw_indels.vcf

gatk SelectVariants -R ${ref} -V ${results}/SKO3_raw_variants.vcf --select-type SNP -O ${results}/SKO3_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/SKO3_raw_variants.vcf --select-type INDEL -O ${results}/SKO3_raw_indels.vcf

gatk SelectVariants -R ${ref} -V ${results}/SKO4_raw_variants.vcf --select-type SNP -O ${results}/SKO4_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/SKO4_raw_variants.vcf --select-type INDEL -O ${results}/SKO4_raw_indels.vcf

gatk SelectVariants -R ${ref} -V ${results}/SKO5_raw_variants.vcf --select-type SNP -O ${results}/SKO5_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/SKO5_raw_variants.vcf --select-type INDEL -O ${results}/SKO5_raw_indels.vcf

gatk SelectVariants -R ${ref} -V ${results}/SKO6_raw_variants.vcf --select-type SNP -O ${results}/SKO6_raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/SKO6_raw_variants.vcf --select-type INDEL -O ${results}/SKO6_raw_indels.vcf



# Filter Variants - GATK4
echo "Filter Variants - GATK4"

# Filter SNPs
gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO1_raw_snps.vcf \
         -O ${results}/SKO1_filtered_snps.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 60.0" \
         -filter-name "MQ_filter" -filter  "MQ < 40.0" \
         -filter-name "SOR_filter" -filter  "SOR > 4.0" \
         -filter-name "MQRankSum_filter" -filter  "MQRankSum < -12.5" \
         -filter-name "ReadPosRankSum_filter" -filter  "ReadPosRankSum < -8.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO2_raw_snps.vcf \
         -O ${results}/SKO2_filtered_snps.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 60.0" \
         -filter-name "MQ_filter" -filter  "MQ < 40.0" \
         -filter-name "SOR_filter" -filter  "SOR > 4.0" \
         -filter-name "MQRankSum_filter" -filter  "MQRankSum < -12.5" \
         -filter-name "ReadPosRankSum_filter" -filter  "ReadPosRankSum < -8.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO3_raw_snps.vcf \
         -O ${results}/SKO3_filtered_snps.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 60.0" \
         -filter-name "MQ_filter" -filter  "MQ < 40.0" \
         -filter-name "SOR_filter" -filter  "SOR > 4.0" \
         -filter-name "MQRankSum_filter" -filter  "MQRankSum < -12.5" \
         -filter-name "ReadPosRankSum_filter" -filter  "ReadPosRankSum < -8.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO4_raw_snps.vcf \
         -O ${results}/SKO4_filtered_snps.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 60.0" \
         -filter-name "MQ_filter" -filter  "MQ < 40.0" \
         -filter-name "SOR_filter" -filter  "SOR > 4.0" \
         -filter-name "MQRankSum_filter" -filter  "MQRankSum < -12.5" \
         -filter-name "ReadPosRankSum_filter" -filter  "ReadPosRankSum < -8.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO5_raw_snps.vcf \
         -O ${results}/SKO5_filtered_snps.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 60.0" \
         -filter-name "MQ_filter" -filter  "MQ < 40.0" \
         -filter-name "SOR_filter" -filter  "SOR > 4.0" \
         -filter-name "MQRankSum_filter" -filter  "MQRankSum < -12.5" \
         -filter-name "ReadPosRankSum_filter" -filter  "ReadPosRankSum < -8.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO6_raw_snps.vcf \
         -O ${results}/SKO6_filtered_snps.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 60.0" \
         -filter-name "MQ_filter" -filter  "MQ < 40.0" \
         -filter-name "SOR_filter" -filter  "SOR > 4.0" \
         -filter-name "MQRankSum_filter" -filter  "MQRankSum < -12.5" \
         -filter-name "ReadPosRankSum_filter" -filter  "ReadPosRankSum < -8.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

# Filter INDELs
gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO1_raw_indels.vcf \
         -O ${results}/SKO1_filtered_indels.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 200.0" 
         -filter-name "SOR_filter" -filter  "SOR > 10.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO2_raw_indels.vcf \
         -O ${results}/SKO2_filtered_indels.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 200.0"
         -filter-name "SOR_filter" -filter  "SOR > 10.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO3_raw_indels.vcf \
         -O ${results}/SKO3_filtered_indels.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 200.0"
         -filter-name "SOR_filter" -filter  "SOR > 10.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO4_raw_indels.vcf \
         -O ${results}/SKO4_filtered_indels.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 200.0"
         -filter-name "SOR_filter" -filter  "SOR > 10.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO5_raw_indels.vcf \
         -O ${results}/SKO5_filtered_indels.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 200.0"
         -filter-name "SOR_filter" -filter  "SOR > 10.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"

gatk VariantFiltration \
         -R ${ref} \
         -V ${results}/SKO6_raw_indels.vcf \
         -O ${results}/SKO6_filtered_indels.vcf \
         -filter-name "QD_filter" -filter  "QD < 2.0" \
         -filter-name "FS_filter" -filter  "FS > 200.0"
         -filter-name "SOR_filter" -filter  "SOR > 10.0" \
         -genotype-filter-expression "DP < 10" \
         -genotype-filter-name "DP_filter" \
         -genotype-filter-expression "GQ < 10" \
         -genotype-filter-name "GQ_filter"



# Select Variants that PASSED filters
# SNPs
gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO1_filtered_snps.vcf \
        -O ${results}/SKO1_analysis-ready-snps.vcf

gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO2_filtered_snps.vcf \
        -O ${results}/SKO2_analysis-ready-snps.vcf

gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO3_filtered_snps.vcf \
        -O ${results}/SKO3_analysis-ready-snps.vcf

gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO4_filtered_snps.vcf \
        -O ${results}/SKO4_analysis-ready-snps.vcf

gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO5_filtered_snps.vcf \
        -O ${results}/SKO5_analysis-ready-snps.vcf


gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO6_filtered_snps.vcf \
        -O ${results}/SKO6_analysis-ready-snps.vcf

# INDELs
gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO1_filtered_indels.vcf \
        -O ${results}/SKO1_analysis-ready-indels.vcf

gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO2_filtered_indels.vcf \
        -O ${results}/SKO2_analysis-ready-indels.vcf2

gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO3_filtered_indels.vcf \
        -O ${results}/SKO3_analysis-ready-indels.vcf

gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO4_filtered_indels.vcf \
        -O ${results}/SKO4_analysis-ready-indels.vcf

gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO5_filtered_indels.vcf \
        -O ${results}/SKO5_analysis-ready-indels.vcf

gatk SelectVariants \
        --exclude-filtered \
        -V ${results}/SKO6_filtered_indels.vcf \
        -O ${results}/SKO6_analysis-ready-indels.vcf



# To exclude variants that failed genotype filters
cat SKO1_analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > SKO1_analysis-ready-snps-filteredGT.vcf
cat SKO1_analysis-ready-indels.vcf|grep -v -E "DP_filter|GQ_filter" > SKO1_analysis-ready-indels-filteredGT.vcf

cat SKO2_analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > SKO2_analysis-ready-snps-filteredGT.vcf
cat SKO2_analysis-ready-indels.vcf|grep -v -E "DP_filter|GQ_filter" > SKO2_analysis-ready-indels-filteredGT.vcf

cat SKO3_analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > SKO3_analysis-ready-snps-filteredGT.vcf
cat SKO3_analysis-ready-indels.vcf|grep -v -E "DP_filter|GQ_filter" > SKO3_analysis-ready-indels-filteredGT.vcf

cat SKO4_analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > SKO4_analysis-ready-snps-filteredGT.vcf
cat SKO4_analysis-ready-indels.vcf|grep -v -E "DP_filter|GQ_filter" > SKO4_analysis-ready-indels-filteredGT.vcf

cat SKO5_analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > SKO5_analysis-ready-snps-filteredGT.vcf
cat SKO5_analysis-ready-indels.vcf|grep -v -E "DP_filter|GQ_filter" > SKO5_analysis-ready-indels-filteredGT.vcf

cat SKO6_analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > SKO6_analysis-ready-snps-filteredGT.vcf
cat SKO6_analysis-ready-indels.vcf|grep -v -E "DP_filter|GQ_filter" > SKO6_analysis-ready-indels-filteredGT.vcf



# Downloading the Funcotator

./gatk FuncotatorDataSourceDownloader --germline --validate-integrity --extract-after-download


# Annotate Variants - GATK4 Funcotator
echo "Annotate Variants - GATK4 Funcotator"

# SNPs
gatk Funcotator \
        --variant ${results}/SKO1_analysis-ready-snps-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSources.v1.7.20200521g \
        --output ${results}/SKO1_analysis-ready-snps-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO2_analysis-ready-snps-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO2_analysis-ready-snps-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO3_analysis-ready-snps-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO3_analysis-ready-snps-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO4_analysis-ready-snps-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO4_analysis-ready-snps-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO5_analysis-ready-snps-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO5_analysis-ready-snps-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO6_analysis-ready-snps-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO6_analysis-ready-snps-filteredGT-funcotated.vcf \
        --output-file-format VCF

# INDELs
gatk Funcotator \
        --variant ${results}/SKO1_analysis-ready-indels-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO1_analysis-ready-indels-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO2_analysis-ready-indels-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO2_analysis-ready-indels-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO3_analysis-ready-indels-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO3_analysis-ready-indels-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO4_analysis-ready-indels-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO4_analysis-ready-indels-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO5_analysis-ready-indels-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO5_analysis-ready-indels-filteredGT-funcotated.vcf \
        --output-file-format VCF

gatk Funcotator \
        --variant ${results}/SKO6_analysis-ready-indels-filteredGT.vcf \
        --reference ${ref} \
        --ref-version hg38 \
        --data-sources-path /home/ulb/gastrexp/kstosic/funcotator/functotator_prepackaged_sources/funcotator/hg38/funcotator_dataSource$
        --output ${results}/SKO6_analysis-ready-indels-filteredGT-funcotated.vcf \
        --output-file-format VCF



# Extracting fields from VCF files to create tab-delimited tables
# SNPs
gatk VariantsToTable \
        -V ${results}/SKO1_analysis-ready-snps-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO1_output_snps.table

gatk VariantsToTable \
        -V ${results}/SKO2_analysis-ready-snps-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO2_output_snps.table

gatk VariantsToTable \
        -V ${results}/SKO3_analysis-ready-snps-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO3_output_snps.table

gatk VariantsToTable \
        -V ${results}/SKO4_analysis-ready-snps-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO4_output_snps.table

gatk VariantsToTable \
        -V ${results}/SKO5_analysis-ready-snps-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO5_output_snps.table

gatk VariantsToTable \
        -V ${results}/SKO6_analysis-ready-snps-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO6_output_snps.table


# INDELs
gatk VariantsToTable \
        -V ${results}/SKO1_analysis-ready-indels-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO1_output_indels.table

gatk VariantsToTable \
        -V ${results}/SKO2_analysis-ready-indels-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO2_output_indels.table

gatk VariantsToTable \
        -V ${results}/SKO3_analysis-ready-indels-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO3_output_indels.table

gatk VariantsToTable \
        -V ${results}/SKO4_analysis-ready-indels-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO4_output_indels.table

gatk VariantsToTable \
        -V ${results}/SKO5_analysis-ready-indels-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO5_output_indels.table

gatk VariantsToTable \
        -V ${results}/SKO6_analysis-ready-indels-filteredGT-funcotated.vcf -F AC -F AN -F DP -F AF -F FUNCOTATION \
        -O ${results}/SKO6_output_indels.table


