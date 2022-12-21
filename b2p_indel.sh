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

#select candidates of snp following set values
perl perl/select_snp.pl ${CURRENT}/${NAME}_merged_final_variants.vcf \
                        ${CURRENT}/${NAME}_selected.vcf \
                        ${MIN_DEP} ${MAX_DEP} ${BETWEEN}
if test $? -ne 0 ; then
  echo "Selecting variants was failed."
  exit 1
fi

#To search primer using primer3, extracting target sequences is neseccary.
#To extract sequences, make the list of target position
perl perl/make_samtools_data.pl ${CURRENT}/${NAME}_selected.vcf \
                                ${CURRENT}/${NAME}_list_samtools.txt \
                                ${SCOPE}
if test $? -ne 0 ; then
  echo "Making data for extracting sequence was failed."
  exit 1
fi

#make fasta file with target sequences.
touch ${CURRENT}/${NAME}_sequences.fasta
while read line
do
  samtools faidx ${REF} ${line} \
            >> ${CURRENT}/${NAME}_sequences.fasta
  if test $? -ne 0 ; then
    echo "Extracting sequence from refernce was failed."
    exit 1
  fi
done < ${CURRENT}/${NAME}_list_samtools.txt
echo 'Extracting sequence from refernce has done.'

#make format for primer3. a lot of text file will be generated in the folder.
mkdir ${CURRENT}/${NAME}_primer3
perl perl/make_format.pl ${CURRENT}/${NAME}_sequences.fasta \
                         ${TEMP} ${SCOPE} ${MARGIN} \
                         ${CURRENT}/${NAME}_primer3
echo 'Making format for primer3 has done.'

#making database for blastn.
if [ ! -e "${REF}.nhr" ]; then
  makeblastdb -in ${REF} \
              -parse_seqids \
              -dbtype nucl
fi

#search primers by primer3 and search possibility of non-specific PCR Products by blastn.　
mkdir ${CURRENT}/${NAME}_blastn
CNT_ALL=0
CNT_FIL=0
while read line
do
    CNT_ALL=`expr $CNT_ALL + 1`
    primer3_core --output ${CURRENT}/${NAME}_primer3/result${line}.txt \
                 ${CURRENT}/${NAME}_primer3/format${line}.txt
    if test $? -ne 0 ; then
      echo "Make primer was failed."
      exit 1
    fi
    rm ${CURRENT}/${NAME}_primer3/format${line}.txt

    perl perl/get_seq_blastn.pl ${CURRENT}/${NAME}_primer3/result${line}.txt
    if test $? -ne 0 ; then
      echo "Extract primer sequence was failed."
      exit 1
    fi

    #If primers were found, the file named "result00000000_query.fasta" will be generated.
    if [ -e ${CURRENT}/${NAME}_primer3/result${line}_query.fasta ]; then
      CNT_FIL=`expr $CNT_FIL + 1`
      blastn -db ${REF} \
             -query ${CURRENT}/${NAME}_primer3/result${line}_query.fasta \
             -out ${CURRENT}/${NAME}_blastn/result${line}_specif.txt \
             -outfmt "6 std qseq sseq sstrand" \
             -evalue 30000 \
             -word_size 7 \
             -num_alignments 50000 \
             -penalty -1 \
             -reward 1 \
             -ungapped

      perl perl/filter_blastn_result.pl \
             ${CURRENT}/${NAME}_blastn/result${line}_specif.txt \
             ${CURRENT}/${NAME}_primer3/result${line}_query.fasta \
             5 \
             1 \
             4000
      if test $? -ne 0 ; then
        echo "Searching possibility of non-specific PCR Products by blastn was failed."
        exit 1
      fi

      rm ${CURRENT}/${NAME}_blastn/result${line}_specif.txt
    fi

    if test $(( ${CNT_ALL} % 50 )) -eq 0 ; then
      echo "${CNT_ALL} variants are processed."
    fi
done < ${CURRENT}/${NAME}_primer3/index.txt
echo "Searching primer by primer3 was finished."
echo "${CNT_FIL} variants with primer in ${CNT_ALL} variants."


perl perl/add_primer_to_vcf.pl ${CURRENT}/${NAME}_primer3 \
                               ${CURRENT}/${NAME}_blastn \
                               ${CURRENT}/${NAME}_selected.vcf \
                               ${CURRENT}/${NAME}_selected_with_primer.vcf

#[ -e "${CURRENT}/format_for_primer3.txt" ] && rm ${CURRENT}/format_for_primer3.txt
#[ -e "${CURRENT}/primer3_result.txt" ] && rm ${CURRENT}/primer3_result.txt

echo 'All done.'
