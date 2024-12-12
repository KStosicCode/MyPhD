#!/bin/bash

#directories
ref="/scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/hg38/Homo_sapiens_assembly38.fasta" 
known_sites="/scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/known_sites" 
aligned_reads="/scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/Aligned_Reads_cfDNA"
reads="/scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/reads"
results="/scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/results"
germlineresource_and_pon="/scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/germlineresource_and_pon"



# Create a fasta index
samtools faidx ${ref}



# Create a sequence dictionary
gatk CreateSequenceDictionary \
  -R ${ref} \
  -O ${ref%.fasta}.dict



# Quality control - FastQC
echo "Quality control - FastQC"

fastqc ${reads}/ULB-CB-ka-hs90-SKO1-cfDNA-pancreas_S1_R1_001.fastq.gz \
       ${reads}/ULB-CB-ka-hs90-SKO1-cfDNA-pancreas_S1_R2_001.fastq.gz \
       -o ${reads}/

fastqc ${reads}/ULB-CB-ka-hs90-SKO2-cfDNA-pancreas_S2_R1_001.fastq.gz \
       ${reads}/ULB-CB-ka-hs90-SKO2-cfDNA-pancreas_S2_R2_001.fastq.gz \
       -o ${reads}/

fastqc ${reads}/ULB-CB-ka-hs90-SKO3-cfDNA-pancreas_S3_R1_001.fastq.gz \
       ${reads}/ULB-CB-ka-hs90-SKO3-cfDNA-pancreas_S3_R2_001.fastq.gz \
       -o ${reads}/

fastqc ${reads}/ULB-CB-ka-hs90-SKO4-cfDNA-pancreas_S4_R1_001.fastq.gz \
       ${reads}/ULB-CB-ka-hs90-SKO4-cfDNA-pancreas_S4_R2_001.fastq.gz \
       -o ${reads}/

fastqc ${reads}/ULB-CB-ka-hs90-SKO5-cfDNA-pancreas_S5_R1_001.fastq.gz \
       ${reads}/ULB-CB-ka-hs90-SKO5-cfDNA-pancreas_S5_R2_001.fastq.gz \
       -o ${reads}/

fastqc ${reads}/ULB-CB-ka-hs90-SKO6-cfDNA-pancreas_S6_R1_001.fastq.gz \
       ${reads}/ULB-CB-ka-hs90-SKO6-cfDNA-pancreas_S6_R2_001.fastq.gz \
       -o ${reads}/



# Trimming - Cutadapt
echo "Trimming - Cutadapt"

cutadapt -q 20 -o ${reads}/SKO1_R1_trimmed.fastq.gz \
         -p ${reads}/SKO1_R2_trimmed.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO1-cfDNA-pancreas_S1_R1_001.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO1-cfDNA-pancreas_S1_R2_001.fastq.gz

cutadapt -q 20 -o ${reads}/SKO2_R1_trimmed.fastq.gz \
         -p ${reads}/SKO2_R2_trimmed.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO2-cfDNA-pancreas_S2_R1_001.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO2-cfDNA-pancreas_S2_R2_001.fastq.gz

cutadapt -q 20 -o ${reads}/SKO3_R1_trimmed.fastq.gz \
         -p ${reads}/SKO3_R2_trimmed.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO3-cfDNA-pancreas_S3_R1_001.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO3-cfDNA-pancreas_S3_R2_001.fastq.gz

cutadapt -q 20 -o ${reads}/SKO4_R1_trimmed.fastq.gz \
         -p ${reads}/SKO4_R2_trimmed.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO4-cfDNA-pancreas_S4_R1_001.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO4-cfDNA-pancreas_S4_R2_001.fastq.gz

cutadapt -q 20 -o ${reads}/SKO5_R1_trimmed.fastq.gz \
         -p ${reads}/SKO5_R2_trimmed.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO5-cfDNA-pancreas_S5_R1_001.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO5-cfDNA-pancreas_S5_R2_001.fastq.gz

cutadapt -q 20 -o ${reads}/SKO6_R1_trimmed.fastq.gz \
         -p ${reads}/SKO6_R2_trimmed.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO6-cfDNA-pancreas_S6_R1_001.fastq.gz \
         ${reads}/ULB-CB-ka-hs90-SKO6-cfDNA-pancreas_S6_R2_001.fastq.gz



# Map to the reference - BWA-MEM
echo "Map to the reference - BWA-MEM"

# 1. BWA index reference 

bwa index ${ref}



# 2. BWA alignment 

bwa mem ${ref} \
    ${reads}/SKO1_R1_trimmed.fastq.gz \
    ${reads}/SKO1_R2_trimmed.fastq.gz > ${aligned_reads}/SKO1.sam

bwa mem ${ref} \
    ${reads}/SKO2_R1_trimmed.fastq.gz \
    ${reads}/SKO2_R2_trimmed.fastq.gz > ${aligned_reads}/SKO2.sam

bwa mem ${ref} \
    ${reads}/SKO3_R1_trimmed.fastq.gz \
    ${reads}/SKO3_R2_trimmed.fastq.gz > ${aligned_reads}/SKO3.sam

bwa mem ${ref} \
    ${reads}/SKO4_R1_trimmed.fastq.gz \
    ${reads}/SKO4_R2_trimmed.fastq.gz > ${aligned_reads}/SKO4.sam

bwa mem ${ref} \
    ${reads}/SKO5_R1_trimmed.fastq.gz \
    ${reads}/SKO5_R2_trimmed.fastq.gz > ${aligned_reads}/SKO5.sam

bwa mem ${ref} \
    ${reads}/SKO6_R1_trimmed.fastq.gz \
    ${reads}/SKO6_R2_trimmed.fastq.gz > ${aligned_reads}/SKO6.sam



# Converting SAM to BAM, sorting BAM, indexing BAM and marking duplicates - gatk MarkDuplicatesSpark
echo "Converting SAM to BAM, sorting BAM, indexing BAM and marking duplicates - gatk MarkDuplicatesSpark"

gatk MarkDuplicatesSpark \
    -I ${aligned_reads}/SKO1.sam \
    -O ${aligned_reads}/SKO1.bam \
    -M ${aligned_reads}/SKO1_duplicates_metrics.txt \
    --CREATE_INDEX true

gatk MarkDuplicatesSpark \
    -I ${aligned_reads}/SKO2.sam \
    -O ${aligned_reads}/SKO2.bam \
    -M ${aligned_reads}/SKO2_duplicates_metrics.txt \
    --CREATE_INDEX true

gatk MarkDuplicatesSpark \
    -I ${aligned_reads}/SKO3.sam \
    -O ${aligned_reads}/SKO3.bam \
    -M ${aligned_reads}/SKO3_duplicates_metrics.txt \
    --CREATE_INDEX true

gatk MarkDuplicatesSpark \
    -I ${aligned_reads}/SKO4.sam \
    -O ${aligned_reads}/SKO4.bam \
    -M ${aligned_reads}/SKO4_duplicates_metrics.txt \
    --CREATE_INDEX true

gatk MarkDuplicatesSpark \
    -I ${aligned_reads}/SKO5.sam \
    -O ${aligned_reads}/SKO5.bam \
    -M ${aligned_reads}/SKO5_duplicates_metrics.txt \
    --CREATE_INDEX true

gatk MarkDuplicatesSpark \
    -I ${aligned_reads}/SKO6.sam \
    -O ${aligned_reads}/SKO6.bam \
    -M ${aligned_reads}/SKO6_duplicates_metrics.txt \
    --CREATE_INDEX true



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
  -O ${results}/SKO1_raw_variants.vcf 

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO2.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO2_f1r2.tar.gz \
  -O ${results}/SKO2_raw_variants.vcf

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO3.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO3_f1r2.tar.gz \
  -O ${results}/SKO3_raw_variants.vcf 
 
gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO4.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO4_f1r2.tar.gz \
  -O ${results}/SKO4_raw_variants.vcf 

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO5.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO5_f1r2.tar.gz \
  -O ${results}/SKO5_raw_variants.vcf 

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/SKO6.bam \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.hg38.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/1000g_pon.hg38.vcf.gz \
  --f1r2-tar-gz ${results}/SKO6_f1r2.tar.gz \
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



# Filtering out artifacts using the generated models and stats
echo "Filtering out artifacts using the generated models and stats"

gatk FilterMutectCalls \
  -V ${results}/SKO1_raw_variants.vcf \
  -R ${ref} \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO1_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO2_raw_variants.vcf \
  -R ${ref} \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO2_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO3_raw_variants.vcf \
  -R ${ref} \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO3_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO4_raw_variants.vcf \
  -R ${ref} \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO4_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO5_raw_variants.vcf \
  -R ${ref} \
  --ob-priors ${results}/reads-orientation-model.tar.gz \
  --stats ${results}/merged.stats \
  -O ${results}/SKO5_filtered_variants.vcf

gatk FilterMutectCalls \
  -V ${results}/SKO6_raw_variants.vcf \
  -R ${ref} \
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
  --data-sources-path /scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/funcotator/funcotator_dataSources.v1.8.hg38.20230908s/ \
  --output ${results}/SKO1_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO2_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/funcotator/funcotator_dataSources.v1.8.hg38.20230908s/ \
  --output ${results}/SKO2_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO3_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/funcotator/funcotator_dataSources.v1.8.hg38.20230908s/ \
  --output ${results}/SKO3_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO4_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/funcotator/funcotator_dataSources.v1.8.hg38.20230908s/ \
  --output ${results}/SKO4_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO5_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/funcotator/funcotator_dataSources.v1.8.hg38.20230908s/ \
  --output ${results}/SKO5_filtered_variants_funcotated.maf \
  --output-file-format MAF

gatk Funcotator \
  --variant ${results}/SKO6_filtered_variants.vcf \
  --reference ${ref} \
  --ref-version hg38 \
  --data-sources-path /scratch/ulb/gastrexp/kstosic/DNA_Variant_Calling/funcotator/funcotator_dataSources.v1.8.hg38.20230908s/ \
  --output ${results}/SKO6_filtered_variants_funcotated.maf \
  --output-file-format MAF

