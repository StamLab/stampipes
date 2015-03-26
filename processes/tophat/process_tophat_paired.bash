source $MODULELOAD
module load bowtie/1.0.0
module load java/jdk1.7.0_05
module load bedops/2.4.2
module load tophat/2.0.13
module load picard/1.118
module load cufflinks/2.2.1
module load samtools/0.1.19
module load bedtools/2.16.2
module load gcc/4.7.2

export SCRIPT_DIR="$STAMPIPES/scripts/tophat"
export REF_DIR=$(dirname "$BWAINDEX")

filesize=$( du --total *fastq.gz | tail -n1 | cut -f1)

# For small files, prioritize overall throughput
# For big files, we want them to finish someday.
if [ "$filesize" -gt 2000000 ] ; then
    THREADS=4
    SLOTS=2
    JOBS=4
else
    THREADS=1
    SLOTS=1
    JOBS=2
fi

if [ ! -e ${SAMPLE_NAME}_R1_fastqc -o ! -e ${SAMPLE_NAME}_R2_fastqc ]; then
qsub -N ".fq${SAMPLE_NAME}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail
  echo "Hostname: " `hostname`
  
  make -f $STAMPIPES/makefiles/fastqc.mk
__SCRIPT__
fi

qsub -cwd -V -q all.q -N .th-$SAMPLE_NAME -now no -pe threads $SLOTS -S /bin/bash <<__MAKE__
  set -e

  if [[ ( -n "$ADAPTER_P7" ) && ( -n "ADAPTER_P5" ) ]] ; then
    echo -e "P7\t$ADAPTER_P7\nP5\t$ADAPTER_P5" > ${SAMPLE_NAME}.adapters.txt
  fi

    make --keep-going -f \$STAMPIPES/makefiles/tophat/tophat_and_cufflinks.mk -j "$JOBS"
__MAKE__
