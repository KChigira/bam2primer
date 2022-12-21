#! /bin/bash
function usage {
  cat <<EOM
Usage: $(basename "$0") -R referance_genome.fasta [OPTION]... Input_BAM_1.bam Input_BAM_2.bam Name_of_BAM_1 Name_of_BAM_2
  -h Display help
  -R [Referance Genome].fasta
  -n analysis_identify_name                default:"sample"
  -d minimum_sequence_depth     default:10
  -p maximum_sequence_depth     default:60
  -t [primer3 template file]    default:"template/template_indel.txt"
EOM
  exit 2
}

REF="f"
BAM_1="f"
BAM_2="f"
NAME_BAM_1="A"
NAME_BAM_2="B"
NAME="sample"
MIN_DEP=10
MAX_DEP=60
TEMP="template/template_indel.txt"

while getopts ":R:n:d:p:h" optKey; do
  case "$optKey" in
    R) REF=${OPTARG};;
    n) NAME=${OPTARG};;
    d) MIN_DEP=${OPTARG};;
    p) MAX_DEP=${OPTARG};;
    t) TEMP=${OPTARG};;
    '-h'|'--help'|* ) usage;;
  esac
done

#${OPTIND} is number of formatted arguments.
#so, shift arguments ${OPTIND} - 1
shift `expr "${OPTIND}" - 1`
if [ $# -ge 4 ]; then
  BAM_1=${1}
  BAM_2=${2}
  NAME_BAM_1=${3}
  NAME_BAM_2=${4}
else
  echo "Not enough arguments."
  usage
fi

#check
[ ! -e ${REF} ] && usage
[ ! -e ${BAM_1} ] && usage
[ ! -e ${BAM_2} ] && usage
[ ${NAME} == "" ] && usage
expr $MIN_DEP + 1 >&/dev/null
[ $? -ge 2 ] && usage
expr $MAX_DEP + 1 >&/dev/null
[ $? -ge 2 ] && usage

#get absolute path
REF=$(cd $(dirname "$REF") && pwd)/$(basename "$REF")
BAM_1=$(cd $(dirname "$BAM_1") && pwd)/$(basename "$BAM_1")
BAM_2=$(cd $(dirname "$BAM_2") && pwd)/$(basename "$BAM_2")

#move to program file directory
CURRENT=$(pwd)
cd $(dirname $0)
[ ! -e ${TEMP} ] && usage

#make log file
LOG=${CURRENT}/${NAME}_log
mkdir ${LOG}

#make a index file of reference genome
samtools faidx ${REF}
picard CreateSequenceDictionary -R ${REF}

#Make VCF from BAM
if [ -d ${CURRENT}/${NAME}_vcf_1st ]; then
  echo "The result of same sample name exists."
  usage
fi
mkdir ${CURRENT}/${NAME}_vcf_1st
bash mkvcf_indel.sh ${REF} \
                    ${BAM_1} \
                    ${CURRENT}/${NAME}_vcf_1st \
                    ${NAME_BAM_1} \
                    ${LOG}/mkvcf_indel_1.log \
                    none \
                    &
bash mkvcf_indel.sh ${REF} \
                    ${BAM_2} \
                    ${CURRENT}/${NAME}_vcf_1st \
                    ${NAME_BAM_2} \
                    ${LOG}/mkvcf_indel_2.log \
                    none \
                    &
wait
#
if test $? -ne 0 ; then
  echo "Haplotype calling was failed."
  exit 1
fi
#
bcftools merge -o ${CURRENT}/${NAME}_vcf_1st/${NAME}_merged_filtered_variants.vcf.gz \
               -O z \
               ${CURRENT}/${NAME}_vcf_1st/${NAME_BAM_1}_filtered_variants.vcf.gz \
               ${CURRENT}/${NAME}_vcf_1st/${NAME_BAM_2}_filtered_variants.vcf.gz
bcftools index ${CURRENT}/${NAME}_vcf_1st/${NAME}_merged_filtered_variants.vcf.gz
#force calling the position where alleles exists.
#To make genotype not "./." but "0/0"
mkdir ${CURRENT}/${NAME}_vcf_2nd
bash mkvcf_indel.sh ${REF} \
                    ${BAM_1} \
                    ${CURRENT}/${NAME}_vcf_2nd \
                    ${NAME_BAM_1} \
                    ${LOG}/mkvcf_indel_3.log \
                    ${CURRENT}/${NAME}_vcf_1st/${NAME}_merged_filtered_variants.vcf.gz \
                    &
bash mkvcf_indel.sh ${REF} \
                    ${BAM_2} \
                    ${CURRENT}/${NAME}_vcf_2nd \
                    ${NAME_BAM_2} \
                    ${LOG}/mkvcf_indel_4.log \
                    ${CURRENT}/${NAME}_vcf_1st/${NAME}_merged_filtered_variants.vcf.gz \
                    &
wait

#
if test $? -ne 0 ; then
  echo "2nd Haplotype calling was failed."
  exit 1
fi

bcftools merge -o ${CURRENT}/${NAME}_merged_final_variants.vcf \
               -O v \
               ${CURRENT}/${NAME}_vcf_2nd/${NAME_BAM_1}_filtered_variants.vcf.gz \
               ${CURRENT}/${NAME}_vcf_2nd/${NAME_BAM_2}_filtered_variants.vcf.gz

#
echo "Making merged vcf was done."
