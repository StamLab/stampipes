source $MODULELOAD
module load samtools/1.3
module load gcc/4.7.2     # R dependent on this
module load R/3.2.5
module load STAR/2.4.2a   # Just for densities
module load bedops/2.4.19
module load subread/1.5.1 # for featureCounts
module load cufflinks/2.2.1 # for cuffLinks

export REFDIR="$(dirname $GENOME_INDEX)"
export STARrefDir="$REFDIR/${STAR_DIR}"
export TARGET_BAM=Aligned.toTranscriptome.out.bam
export GENOME_BAM=Aligned.toGenome.out.bam

export QUEUE=queue0

numbam=$(wc -w <<< $BAM_FILES)
# Temporary
if [ ! -s "$TARGET_BAM" ] ; then
  if [ $numbam -eq 1 ] ; then
    cp "$BAM_FILES" "$TARGET_BAM"
  else
    samtools merge -n -f "$TARGET_BAM" $BAM_FILES
  fi
fi

if [ ! -s "$GENOME_BAM" ] ; then
  GENOME_BAM_FILES=$(sed 's/toTranscriptome/sortedByCoord/g' <<< "$BAM_FILES")
  $STAMPIPES/scripts/tophat/merge_or_copy_bam.sh "$GENOME_BAM" $GENOME_BAM_FILES
  samtools index "$GENOME_BAM"
fi

# density information
if [ ! -s "Signal.UniqueMultiple.str+.starch" ] ; then
  qsub -cwd -V -N ".AGG${AGGREGATION_ID}.star_den" -q $QUEUE -S /bin/bash <<'__DEN__'

    # Write starch and bigwig to .tmp files
    function convertBedGraph(){
      in="$1"
      base="$2"
      chrom="$in.onlyChr.bg"
      grep '^chr' "$in" | sort -k1,1 -k2,2n > $chrom
      bedGraphToBigWig "$chrom" chrNL.txt "$base.bw.tmp"
      starch "$chrom" > "$base.starch.tmp"
    }

    set -x -e -o pipefail

    mkdir -p $TMPDIR/Signal

    echo STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir $TMPDIR/STAR
    STAR --runMode inputAlignmentsFromBAM --inputBAMfile $GENOME_BAM --outWigType bedGraph --outWigStrand Stranded --outFileNamePrefix $TMPDIR/Signal/ --outWigReferencesPrefix chr --outTmpDir $TMPDIR/STAR

    grep '^chr' $STARrefDir/chrNameLength.txt > chrNL.txt

    convertBedGraph $TMPDIR/Signal/Signal.Unique.str1.out.bg         Signal.Unique.str-
    convertBedGraph $TMPDIR/Signal/Signal.Unique.str2.out.bg         Signal.Unique.str+
    convertBedGraph $TMPDIR/Signal/Signal.UniqueMultiple.str1.out.bg Signal.UniqueMultiple.str-
    convertBedGraph $TMPDIR/Signal/Signal.UniqueMultiple.str2.out.bg Signal.UniqueMultiple.str+

    for i in Signal*.tmp ; do
      mv $i ${i/.tmp/}
    done

__DEN__
fi

# cufflinks
if [ ! -s "genes.fpkm_tracking" ] ; then
qsub -cwd -V -N ".AGG${AGGREGATION_ID}.star_cuff" -q $QUEUE -S /bin/bash <<'__CUFF__'
    CUFF=cufflinks
    CUFF_COMMON="--no-update-check --library-type fr-firststrand"
    $CUFF $CUFF_COMMON --GTF $ANNOTATION $GENOME_BAM
__CUFF__
fi

# featureCounts
if [ ! -s "feature_counts.txt" ] ; then
qsub -cwd -V -N ".AGG${AGGREGATION_ID}.star_fcounts" -q $QUEUE -S /bin/bash <<'__FCOUNTS__'
    FCOUNTS=featureCounts
    FCOUNTS_COMMON="--primary -B -C -p -P --fracOverlap .5 -s 2"
    $FCOUNTS $FCOUNTS_COMMON -t 'exon' -g 'gene_id' -a $ANNOTATION -o feature_counts.txt $GENOME_BAM
__FCOUNTS__
fi

qsub -N ".AGG#${AGGREGATION_ID}.complete" -hold_jid ".star_*" -V -cwd -q $QUEUE -S /bin/bash > /dev/stderr <<__SCRIPT__
    bash $STAMPIPES/scripts/rna-star/aggregate/checkcomplete.sh
    bash $STAMPIPES/scripts/rna-star/aggregate/attachfiles.sh
__SCRIPT__
