#! /bin/bash
function usage {
  cat <<EOM
Usage: $(basename "$0") -R referance_genome.fasta [OPTION]... Input_VCF.vcf(Made by b2p)
  -h Display help
  -R [Referance Genome].fasta
EOM
  exit 2
}

REF="f"
VCF="f"

while getopts ":R:h" optKey; do
  case "$optKey" in
    R) REF=${OPTARG};;
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

#get absolute path
REF=$(cd $(dirname "$REF") && pwd)/$(basename "$REF")
VCF=$(cd $(dirname "$VCF") && pwd)/$(basename "$VCF")

#move to program file directory
CURRENT=$(pwd)
cd $(dirname $0)

#make a index file of reference genome if absent.
if [ ! -e ${REF}.fai ]; then
  samtools faidx ${REF}
fi

perl perl/select_markers_indel.pl \
       ${VCF}

#
if test $? -ne 0 ; then
  echo "Select markers was failed."
  exit 1
fi

VCF_STEM=`echo $VCF | sed -e "s/\.[^.]*$//"` #delete extention
VCF_NEW=${VCF_STEM}_available.vcf
PNG=${VCF_STEM}__available.png
Rscript R/visualize_SNP.R ${VCF_NEW} \
                          ${REF}.fai \
                          ${PNG}

echo 'All done.'
