#! /bin/bash
function usage {
  cat <<EOM
Usage: $(basename "$0") -R referance_genome.fasta [OPTION]... Input_BAM_1.bam Input_BAM_2.bam Name_of_BAM_1 Name_of_BAM_2
  -h Display help
  -R [Referance Genome].fasta
  -n analysis_identify_name                default:"sample"
EOM
  exit 2
}

REF="f"
BAM_1="f"
BAM_2="f"
NAME_BAM_1="A"
NAME_BAM_2="B"
NAME="sample"

while getopts ":R:n:h" optKey; do
  case "$optKey" in
    R) REF=${OPTARG};;
    n) NAME=${OPTARG};;
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

#get absolute path
REF=$(cd $(dirname "$REF") && pwd)/$(basename "$REF")
BAM_1=$(cd $(dirname "$BAM_1") && pwd)/$(basename "$BAM_1")
BAM_2=$(cd $(dirname "$BAM_2") && pwd)/$(basename "$BAM_2")

#move to program file directory
CURRENT=$(pwd)
cd $(dirname $0)

#make log file
LOG=${CURRENT}/${NAME}_log

bash mkvcf_indel_221222only.sh ${REF} \
                    ${BAM_1} \
                    ${CURRENT}/${NAME}_vcf_2nd \
                    ${NAME_BAM_1} \
                    ${LOG}/mkvcf_indel_3.log \
                    ${CURRENT}/${NAME}_vcf_1st/${NAME}_merged_filtered_variants.vcf.gz \
                    &
bash mkvcf_indel_221222only.sh ${REF} \
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
