source $PYTHON3_ACTIVATE

source $MODULELOAD
module load bedops/2.4.19
module load jdk/1.8.0_92
module load gcc/4.7.2
module load R/3.2.5
module load picard/1.120
module load samtools/1.3
module load git/2.3.3
module load coreutils/8.25

module load python/2.7.11
module load cufflinks/2.2.1

cd $AGGREGATION_FOLDER

export LIBRARY_NAME=LN${LIBRARY}
export REF_DIR=$(dirname "$GENOME_INDEX")
export SCRIPT_DIR="$STAMPIPES/scripts/tophat"

BAM_COUNT=`ls $BAM_FILES | wc -l`

JOB_BASENAME=".AGG#${AGGREGATION_ID}"

export FINAL_BAM=${LIBRARY_NAME}.all.bam

if [ -n "$REDO_AGGREGATION" ]; then
  # TODO: Make
  bash $STAMPIPES/scripts/tophat/aggregate/basic/reset.bash
fi

MERGE_DUP_JOBNAME=${JOB_BASENAME}_merge_dup
COUNT_JOBNAME=${JOB_BASENAME}_count
MAKE_JOBNAME=${JOB_BASENAME}_make

# Check out files match first
python3 $STAMPIPES/scripts/utility/md5check.py bamfiles.txt || exit 1

PROCESSING=""

if [ ! -e ${FINAL_BAM} ]; then

  PROCESSING="$PROCESSING,${MERGE_DUP_JOBNAME}"

  qsub ${SUBMIT_SLOTS} -N "${MERGE_DUP_JOBNAME}" -V -cwd -S /bin/bash > /dev/stderr << __MERGE__
    set -x -e -o pipefail

    TMPBAM="\$TMPDIR/$LIBRARY_NAME.bam"

    echo "Hostname: "
    hostname

    echo "START: "
    date

    $STAMPIPES/scripts/tophat/merge_or_copy_bam.sh \$TMPBAM  $BAM_FILES

    java -Xmx24g -jar $PICARDPATH/MarkDuplicates.jar \
      INPUT=\$TMPBAM \
      METRICS_FILE=$LIBRARY_NAME.dups.txt \
      OUTPUT=$FINAL_BAM \
      REMOVE_DUPLICATES=false \
      ASSUME_SORTED=true \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'


    echo "FINISH: "
    date

__MERGE__

fi

# Currently single-threaded
make_strand(){
  name=$1
  targets="$name.bam $name.starch $name.bw $name.bam.bai"
  makefile="$STAMPIPES/makefiles/tophat/tophat_and_cufflinks.mk"

  if ! make -q -f $makefile $targets ; then
    PROCESSING="$PROCESSING,${MAKE_JOBNAME}-$name"
    qsub -N "${MAKE_JOBNAME}-$name" -hold_jid ${MERGE_DUP_JOBNAME} -V -cwd -S /bin/bash > /dev/stderr << __STRAND__
      set -x -e -o pipefail

      echo "Hostname: "
      hostname

      echo "START: "
      date

      make \
        -f $STAMPIPES/makefiles/tophat/tophat_and_cufflinks.mk \
        $targets

      echo "FINISH: "
      date

__STRAND__
  fi
}

do_cufflinks(){
  targets="${LIBRARY_NAME}_cufflinks/finished.txt"
  makefile="$STAMPIPES/makefiles/tophat/tophat_and_cufflinks.mk"

  export SAMPLE_NAME=$LIBRARY_NAME
  if ! make -q -f $makefile $targets ; then

    jobname="${MAKE_JOBNAME}-cufflinks" 
    PROCESSING="$PROCESSING,$jobname"
    tmpbamlink="$LIBRARY_NAME.all.$GENOME.bam"

    qsub -N $jobname -hold_jid ${MERGE_DUP_JOBNAME} -V -cwd -S /bin/bash > /dev/stderr -pe threads 4 << __CUFFLINKS__
    set -x -e -o pipefail

    echo "Hostname: "
    hostname

    echo "START: "
    date

    rm -f $tmpbamlink
    ln -s $FINAL_BAM  $tmpbamlink

    make \
      -f $STAMPIPES/makefiles/tophat/tophat_and_cufflinks.mk \
      $targets

    rm $tmpbamlink

    echo "FINISH: "
    date

__CUFFLINKS__
  fi
}

## Main ##

make_strand $LIBRARY_NAME.all
make_strand $LIBRARY_NAME.pos
make_strand $LIBRARY_NAME.neg

do_cufflinks


if [ -n ${PROCESSING} ]; then

  qsub -N "${JOB_BASENAME}_complete" -hold_jid "${PROCESSING}" -V -cwd -S /bin/bash > /dev/stderr << __SCRIPT__
  set -x -e -o pipefail

  echo "Hostname: "
  hostname

  echo "START: "
  date

  bash $STAMPIPES/scripts/tophat/aggregate/basic/checkcomplete.bash
  bash $STAMPIPES/scripts/tophat/aggregate/basic/attachfiles.bash

  echo "FINISH: "
  date

__SCRIPT__

fi
