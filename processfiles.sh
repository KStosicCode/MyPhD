

directory=/home/kobina/exp2/bams

for file in `ls $directory/*.bam`
do


samplename=$(basename $file .bam)

gatk Mutect2 \
  -R ${ref} \
  -I ${aligned_reads}/${file} \
  --germline-resource ${germlineresource_and_pon}/af-only-gnomad.vcf.gz \
  --panel-of-normals ${germlineresource_and_pon}/pon.vcf.gz \
  --f1r2-tar-gz ${results}/${samplename}_f1r2.tar.gz \
  --stats ${results}/${samplename}.raw_variants.vcf.stat \
  -O ${results}/${samplename}_raw_variants.vcf


done 
