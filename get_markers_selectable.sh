#! /bin/bash
function usage {
  cat <<EOM
Usage: $(basename "$0") -R referance_genome.fasta [OPTION]... Input_VCF.vcf(Made by b2p)
  -h Display help
  -R [Referance Genome].fasta
  -n number of markers required (Integer)      default:384
  -i interval of 2 markers (bp)                default:1000000
  -x select column number of vcf to use as genotype A   default:9
  -y select column number of vcf to use as genotype B   default:10
EOM
  exit 2
}

REF="f"
VCF="f"
NUM_MARKER=384
INTERVAL=1000000
COL_X=9
COL_Y=10

while getopts ":R:n:i:x:y:h" optKey; do
  case "$optKey" in
    R) REF=${OPTARG};;
    n) NUM_MARKER=${OPTARG};;
    i) INTERVAL=${OPTARG};;
    x) COL_X=${OPTARG};;
    y) COL_Y=${OPTARG};;
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
expr $NUM_MARKER + 1 >&/dev/null
[ $? -ge 2 ] && usage
expr $INTERVAL + 1 >&/dev/null
[ $? -ge 2 ] && usage

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

perl perl/select_markers_to_use_selectable.pl \
       ${VCF} \
       ${REF}.fai \
       ${NUM_MARKER} \
       ${INTERVAL} ${COL_X} ${COL_Y}

#
if test $? -ne 0 ; then
  echo "Select markers was failed."
  exit 1
fi

VCF_STEM=`echo $VCF | sed -e "s/\.[^.]*$//"` #delete extention
VCF_NEW=${VCF_STEM}_select_${NUM_MARKER}_mk_or_${INTERVAL}_bp.vcf
PNG=${VCF_STEM}_select_${NUM_MARKER}_mk_or_${INTERVAL}_bp.png
Rscript R/visualize_SNP.R ${VCF_NEW} \
                          ${REF}.fai \
                          ${PNG}

echo 'All done.'
