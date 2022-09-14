#! /bin/sh

if [ $# -ne 5 ]; then
  echo "Error: 5 arguments must be contained."
  exit 1;
fi

REF=$1
INPUT=$2
OUTDIR=$3
NAME=$4
VCF=$5

if [ ! -e ${REF} ] || [ ! -e ${INPUT} ]; then
  echo "Error: Input file does not exist."
  exit 1;
fi

if [ ! -d ${OUTDIR} ]; then
  echo "Error: output directory does not exist."
fi

gatk IndexFeatureFile -I ${VCF}
gatk HaplotypeCaller \
       -R ${REF} \
       -I ${INPUT} \
       -O ${OUTDIR}/${NAME}_raw_variants_recal.vcf \
       --alleles ${VCF}

gatk SelectVariants \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_variants_recal.vcf \
       -select-type SNP \
       -O ${OUTDIR}/${NAME}_raw_snps_recal.vcf

gatk SelectVariants \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_variants_recal.vcf \
       -select-type INDEL \
       -O ${OUTDIR}/${NAME}_raw_indels_recal.vcf

gatk VariantFiltration \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_snps_recal.vcf \
       -O ${OUTDIR}/${NAME}_filtered_snps_recal.vcf \
       -filter-name "QD_filter" -filter "QD < 20.0" \
       -filter-name "FS_filter" -filter "FS > 60.0" \
       -filter-name "MQ_filter" -filter "MQ < 40.0" \
       -filter-name "SOR_filter" -filter "SOR > 4.0" \
       -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
       -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \

gatk VariantFiltration \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_indels_recal.vcf \
       -O ${OUTDIR}/${NAME}_filtered_indels_recal.vcf \
       -filter-name "QD_filter" -filter "QD < 20.0" \
       -filter-name "FS_filter" -filter "FS > 200.0" \
       -filter-name "SOR_filter" -filter "SOR > 10.0"

#
bgzip -c ${OUTDIR}/${NAME}_filtered_snps_recal.vcf \
      > ${OUTDIR}/${NAME}_filtered_snps_recal.vcf.gz
bcftools index ${OUTDIR}/${NAME}_filtered_snps_recal.vcf.gz

bgzip -c ${OUTDIR}/${NAME}_filtered_indels_recal.vcf \
      > ${OUTDIR}/${NAME}_filtered_indels_recal.vcf.gz
bcftools index ${OUTDIR}/${NAME}_filtered_indels_recal.vcf.gz

bcftools concat -o ${OUTDIR}/${NAME}_filtered_variants_recal.vcf.gz \
                -a -O z \
                ${OUTDIR}/${NAME}_filtered_snps_recal.vcf.gz \
                ${OUTDIR}/${NAME}_filtered_indels_recal.vcf.gz
bcftools index ${OUTDIR}/${NAME}_filtered_variants_recal.vcf.gz
