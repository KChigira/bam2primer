#! /bin/bash
function usage {
  cat <<EOM
Usage: $(basename "$0") -R referance_genome.fasta [OPTION]... Input_VCF.vcf
  -h Display help
  -R [Referance Genome].fasta
  -n analysis_identify_name                default:"sample"
  -d minimum_sequence_depth     default:10
  -p maximum_sequence_depth     default:60
  -b minimum_length_between_2_variants    default:150
  -s scope where primer designed          default:140
  -m margin of target sequence where primers don't locate    default:10
  -t [primer3 template file]    default:"template/template_snp.txt"
EOM
  exit 2
}

REF="f"
VCF="f"
NAME="sample"
MIN_DEP=10
MAX_DEP=60
BETWEEN=150
SCOPE=140
MARGIN=10
TEMP="template/template_snp.txt"

while getopts ":R:n:d:p:b:s:m:h" optKey; do
  case "$optKey" in
    R) REF=${OPTARG};;
    n) NAME=${OPTARG};;
    d) MIN_DEP=${OPTARG};;
    p) MAX_DEP=${OPTARG};;
    b) BETWEEN=${OPTARG};;
    s) SCOPE=${OPTARG};;
    m) MARGIN=${OPTARG};;
    t) TEMP=${OPTARG};;
    '-h'|'--help'|* ) usage;;
  esac
done

#${OPTIND} is number of formatted arguments.
#so, shift arguments ${OPTIND} - 1
shift `expr "${OPTIND}" - 1`
if [ $# -ge 1 ]; then
  VCF=${1}
else
  echo "Not enough arguments."
  usage
fi

#check
[ ! -e ${REF} ] && usage
[ ! -e ${VCF} ] && usage
[ ${NAME} == "" ] && usage
expr $MIN_DEP + 1 >&/dev/null
[ $? -ge 2 ] && usage
expr $MAX_DEP + 1 >&/dev/null
[ $? -ge 2 ] && usage
expr $BETWEEN + 1 >&/dev/null
[ $? -ge 2 ] && usage
expr $SCOPE + 1 >&/dev/null
[ $? -ge 2 ] && usage
expr $MARGIN + 1 >&/dev/null
[ $? -ge 2 ] && usage

#get absolute path
REF=$(cd $(dirname "$REF") && pwd)/$(basename "$REF")
VCF=$(cd $(dirname "$VCF") && pwd)/$(basename "$VCF")

#move to program file directory
CURRENT=$(pwd)
cd $(dirname $0)
[ ! -e ${TEMP} ] && usage


#select candidates of snp following set values
perl perl/select_snp_3vs1.pl ${VCF} \
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
