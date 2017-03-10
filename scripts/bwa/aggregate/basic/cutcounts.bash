set -x -e -o pipefail

export LIBRARY_NAME=LN${LIBRARY}

cd $AGGREGATION_FOLDER

SAMPLE_NAME=$LIBRARY_NAME


GENOME_MAPPABILITY_FILE=/net/fileserv0/vol7/annotations/data/${GENOME}/${GENOME}.K${READ_LENGTH}.mappable_only.bed

export CUTS_BED=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.cuts.sorted.bed.starch
export CUTCOUNTS=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.cutcounts.sorted.bed.starch
export FRAGMENTS=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.fragments.sorted.bed.starch
export CUTCOUNTS_BIGWIG=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.cutcounts.$READ_LENGTH.bw

# temp files
COUNTFILETMP=$TMPDIR/base-count.$LIBRARY_NAME.perBase.$READ_LENGTH.$GENOME.bed
ALLBASE=$TMPDIR/all-perBase.$LIBRARY_NAME.$READ_LENGTH.$GENOME.temp
BEDTMP=$TMPDIR/$LIBRARY_NAME.count.uniques.$GENOME.bed
WIGTMP=$TMPDIR/$LIBRARY_NAME.count.uniques.$GENOME.wig
CHROMTMP=$TMPDIR/chroms.bed
CUTSTMP=$TMPDIR/cuts.bed
FRAGMENTSTMP=$TMPDIR/fragments.bed

date

# Create cut counts and fragments if they don't exist
if [ ! -s $CUTS_BED ]; then

# Create the uniques BAM if it doesn't currently exist
if [ ! -o $AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.uniques.sorted.bam ]; then
  export UNIQUES_BAMFILE=$TMPDIR/$LIBRARY_NAME.${GENOME}.uniques.sorted.bam
  bash $STAMPIPES/scripts/bwa/aggregate/basic/bamfilter.bash $AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.sorted.bam $UNIQUES_BAMFILE
else
  export UNIQUES_BAMFILE=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.uniques.sorted.bam
fi

# Convert a BAM file into fragments and cut counts
time bam2bed --do-not-sort < $UNIQUES_BAMFILE \
  | awk -v cutfile=$CUTSTMP -v fragmentfile=$FRAGMENTSTMP -f ${STAMPIPES}/scripts/bwa/aggregate/basic/cutfragments.awk

if [ -s "$FRAGMENTSTMP" ] ; then
  sort-bed --max-mem 16G $FRAGMENTSTMP | starch - > $FRAGMENTS
fi
sort-bed --max-mem 16G $CUTSTMP | starch - > $CUTS_BED 

fi

if [ ! -s $CUTCOUNTS ]; then

time unstarch $CUTS_BED \
  | cut -f1-3 \
  | bedops -m - \
  | awk '{ for(i = $2; i < $3; i += 1) { print $1"\t"i"\t"i + 1 }}' \
  > $ALLBASE

time unstarch $CUTS_BED | \
  bedmap --echo --count --delim "\t" ${ALLBASE} - \
  | awk '{print $1"\t"$2"\t"$3"\tid-"NR"\t"$4}' \
  > ${COUNTFILETMP}

time starch $COUNTFILETMP > $CUTCOUNTS

fi

if [ ! -s $CUTCOUNTS_BIGWIG ]; then

# removing mitochonrial, very large chrM datasets can cause problems with the browser.
time awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$5}' $COUNTFILETMP | grep -v chrM > $WIGTMP

time wigToBigWig -clip $WIGTMP ${GENOME_INDEX}.fai $CUTCOUNTS_BIGWIG

fi

date
