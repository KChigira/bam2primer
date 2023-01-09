#! /bin/sh

if [ $# -ne 6 ]; then
  echo "Error: 6 arguments must be contained."
  exit 1;
fi

REF=$1
INPUT=$2
OUTDIR=$3
NAME=$4
LOG=$5
VCF=$6

if [ ! -e ${REF} ] || [ ! -e ${INPUT} ]; then
  echo "Error: Input file does not exist."
  exit 1;
fi

if [ ! -d ${OUTDIR} ]; then
  echo "Error: output directory does not exist."
fi

touch ${LOG}

samtools index ${INPUT}

gatk SelectVariants \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_variants.vcf \
       -select-type SNP \
       -O ${OUTDIR}/${NAME}_raw_snps.vcf \
       >> ${LOG}

gatk SelectVariants \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_variants.vcf \
       -select-type INDEL \
       -O ${OUTDIR}/${NAME}_raw_indels.vcf \
       >> ${LOG}

gatk VariantFiltration \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_snps.vcf \
       -O ${OUTDIR}/${NAME}_filtered_snps.vcf \
       -filter-name "QD_filter" -filter "QD < 20.0" \
       -filter-name "FS_filter" -filter "FS > 60.0" \
       -filter-name "MQ_filter" -filter "MQ < 40.0" \
       -filter-name "SOR_filter" -filter "SOR > 4.0" \
       -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
       -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
       >> ${LOG}

gatk VariantFiltration \
       -R ${REF} \
       -V ${OUTDIR}/${NAME}_raw_indels.vcf \
       -O ${OUTDIR}/${NAME}_filtered_indels.vcf \
       -filter-name "QD_filter" -filter "QD < 20.0" \
       -filter-name "FS_filter" -filter "FS > 200.0" \
       -filter-name "SOR_filter" -filter "SOR > 10.0" \
       >> ${LOG}

#
bgzip -c ${OUTDIR}/${NAME}_filtered_snps.vcf \
      > ${OUTDIR}/${NAME}_filtered_snps.vcf.gz
bcftools index ${OUTDIR}/${NAME}_filtered_snps.vcf.gz

bgzip -c ${OUTDIR}/${NAME}_filtered_indels.vcf \
      > ${OUTDIR}/${NAME}_filtered_indels.vcf.gz
bcftools index ${OUTDIR}/${NAME}_filtered_indels.vcf.gz

bcftools concat -o ${OUTDIR}/${NAME}_filtered_variants.vcf.gz \
                -a -O z \
                ${OUTDIR}/${NAME}_filtered_snps.vcf.gz \
                ${OUTDIR}/${NAME}_filtered_indels.vcf.gz
bcftools index ${OUTDIR}/${NAME}_filtered_variants.vcf.gz
