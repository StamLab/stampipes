#!/bin/bash

version=2.4.2
export NXF_VER=18.10.1  # The version of nextflow to run. 18.10.1 includes conda

cd "$(dirname "$0")"

source "$MODULELOAD"
module purge
module load jdk
module load nextflow
module load python/3.5.1
module load anaconda/2.1.0-dev

if [[ $(wc -w <<< "$BAM_FILES") -gt 1 ]] ; then
  bamfiles="$(sed 's/\s\+/,/g' <<< "$BAM_FILES")"
else
  bamfiles=$BAM_FILES
fi
outdir="output_${version}"
workdir="work"

WIN=75
BINI=20


export MAPPABLE_REGIONS=${MAPPABLE_REGIONS:-$GENOME_INDEX.K${READ_LENGTH}.mappable_only.bed}
export CHROM_SIZES=${CHROM_SIZES:-$GENOME_INDEX.chrom_sizes.bed}
export CENTER_SITES=${CENTER_SITES:-$GENOME_INDEX.K${READ_LENGTH}.center_sites.n100.nuclear.starch}
export NUCLEAR_CHR=${NUCLEAR_CHR:-$GENOME_INDEX.nuclear.txt}
export CHROM_BUCKET=$STAMPIPES_DATA/densities/chrom-buckets.$GENOME.${WIN}_${BINI}.bed.starch

HOTSPOT_INDEX=${HOTSPOT_INDEX:-.}

# Remove old stuff if necessary
if [[ -n "$REDO_AGGREGATION" ]] ; then
  rm -rf "$outdir"
  rm -rf "$workdir"
  echo \
  python3 "$STAMPIPES/scripts/lims/upload_data.py" \
    --clear_aggregation_stats \
    --aggregation_id "$AGGREGATION_ID"
fi

# Tell LIMS we're starting alignment
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  -a "$LIMS_API_URL" \
  -t "$LIMS_API_TOKEN" \
  --aggregation_id "$AGGREGATION_ID" \
  --start_aggregation

# Default peak caller for back-compat
if [[ -z "$PEAK_CALLER" ]] ; then
  # TODO: Match on ASSAY_CATEGORY instead
  if [[ "$ASSAY" =~ ChIP ]] || [[ "$ASSAY" == "Crosslinked" ]] ; then
    PEAK_CALLER=macs2
  else
    PEAK_CALLER=hotspot2
  fi
fi

# Run the whole process
nextflow run \
  "$STAMPIPES/processes/bwa/aggregate/basic.nf" \
  -c "$STAMPIPES/nextflow.config" \
  -w "$workdir" \
  --bams "$bamfiles" \
  --genome "$GENOME_INDEX" \
  --mappable "$MAPPABLE_REGIONS" \
  --chrom_sizes "$CHROM_SIZES" \
  --peakcaller "$PEAK_CALLER" \
  --centers "$CENTER_SITES" \
  --chrom_bucket "$CHROM_BUCKET" \
  --hotspot_index "$HOTSPOT_INDEX" \
  --hotspot_id "AG$AGGREGATION_ID" \
  --bias "$STAMPIPES_DATA/footprints/vierstra_et_al.txt" \
  --UMI="$UMI" \
  --outdir "$outdir" \
  --threads 3 \
  -profile cluster,modules \
  -with-report nextflow.report.html \
  -with-dag nextflow.flowchart.html \
  -with-timeline nextflow.timeline.html \
  -resume

export PEAK_CALLER
( cd "$outdir" \
    && bash "$STAMPIPES/scripts/bwa/aggregate/basic/attachfiles_nextflow.bash" \
    && bash "$STAMPIPES/scripts/bwa/aggregate/basic/uploadcounts_nextflow.bash"
)
# Upload results
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  -a "$LIMS_API_URL" \
  -t "$LIMS_API_TOKEN" \
  --aggregation_id "$AGGREGATION_ID" \
  --complete_aggregation

# TODO: remove below this line
####
exit
####


source $MODULELOAD
module load kentutil/302
module load bedops/2.4.35-typical
module load jdk/1.8.0_92
module load gcc/4.7.2
module load R/3.2.5
module load picard/2.8.1
module load samtools/1.3
module load git/2.3.3
module load coreutils/8.25
module load modwt/1.0
module load bedtools/2.25.0
module load python/3.5.1
module load pysam/0.9.0
module load htslib/1.6.0
module load hotspot2/2.1.1

module load numpy/1.11.0
module load atlas-lapack/3.10.2
module load scipy/1.0.0
module load scikit-learn/0.18.1

WIN=75
BINI=20

export FINAL_BAM=${LIBRARY_NAME}.${GENOME}.sorted.bam
export FINAL_BAM_MARKED=${LIBRARY_NAME}.${GENOME}.sorted.marked.bam
export FINAL_UNIQUES_BAM=${LIBRARY_NAME}.${GENOME}.uniques.sorted.bam
export TAGCOUNTS_FILE=${LIBRARY_NAME}.tagcounts.txt
export ADAPTER_COUNT_FILE=${LIBRARY_NAME}.adaptercounts.txt
export DENSITY_STARCH=${LIBRARY_NAME}.${WIN}_${BINI}.${GENOME}.uniques-density.bed.starch
export DENSITY_BIGWIG=${LIBRARY_NAME}.${WIN}_${BINI}.${GENOME}.bw
export NORM_DENSITY_STARCH=${LIBRARY_NAME}.${WIN}_${BINI}.normalized.${GENOME}.uniques-density.bed.starch
export NORM_DENSITY_BIGWIG=${LIBRARY_NAME}.${WIN}_${BINI}.normalized.${GENOME}.bw
export CUTCOUNTS_STARCH=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.cutcounts.sorted.bed.starch
export CUTCOUNTS_BIGWIG=$AGGREGATION_FOLDER/$LIBRARY_NAME.${GENOME}.cutcounts.$READ_LENGTH.bw
export INSERT_FILE=${LIBRARY_NAME}.CollectInsertSizeMetrics.picard
export DUPS_FILE=${LIBRARY_NAME}.MarkDuplicates.picard
export PRESEQ_HIST=${LIBRARY_NAME}.uniques.duphist.txt
export PRESEQ_RES=${LIBRARY_NAME}.uniques.preseq.txt
export PRESEQ_TRGT=${LIBRARY_NAME}.uniques.preseq.targets.txt

export HOTSPOT2_DIR=peaks_v2_1_1
HOTSPOT_PREFIX=$(basename "$FINAL_UNIQUES_BAM" .bam)
export HOTSPOT_CALLS=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.hotspots.fdr0.05.starch
export HOTSPOT_CALLS_01=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.hotspots.fdr0.01.starch
export HOTSPOT_CALLS_001=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.hotspots.fdr0.001.starch
export HOTSPOT_DENSITY=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.density.bw
export HOTSPOT_ALLCALLS=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.allcalls.starch
export HOTSPOT_CUTCOUNTS=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.cutcounts.starch
export HOTSPOT_CLEAVAGES=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.cleavage.total
export HOTSPOT_ISPOT=$HOTSPOT_PREFIX.ispot.info
export HOTSPOT_PROXDIST=$HOTSPOT_PREFIX.proxdist.info
export HOTSPOT_SCRIPT="hotspot2.sh"
export MAPPABLE_REGIONS=${MAPPABLE_REGIONS:-$GENOME_INDEX.K${READ_LENGTH}.mappable_only.bed}
export CHROM_SIZES=${CHROM_SIZES:-$GENOME_INDEX.chrom_sizes.bed}
export CENTER_SITES=${CENTER_SITES:-$GENOME_INDEX.K${READ_LENGTH}.center_sites.n100.nuclear.starch}
export NUCLEAR_CHR=${NUCLEAR_CHR:-$GENOME_INDEX.nuclear.txt}
export HOTSPOT_PEAKS_DEF=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.peaks.starch
export HOTSPOT_NPEAKS_DEF=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.narrowpeaks.starch
export HOTSPOT_NPEAKS_05=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.narrowpeaks.fdr0.05.starch
export HOTSPOT_PEAKS_05=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.peaks.fdr0.05.starch
export HOTSPOT_PEAKS_01=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.peaks.fdr0.01.starch
export HOTSPOT_PEAKS_001=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.peaks.fdr0.001.starch

# hard-coded until we retroactively change aggregation templates to their genome-specific template
export ALTIUS_MASTERLIST="/net/seq/data/genomes/human/GRCh38/noalts/ref/masterlist_DHSs_WM20180313_all_indexIDs.665samples.txt"
export FIMO_TRANSFAC_1E4="/net/seq/data/genomes/human/GRCh38/noalts/ref/fimo.combined.1e-4.parsed.starch"
export FIMO_NAMES="/net/seq/data/genomes/human/GRCh38/noalts/ref/fimo.transfac.names.txt"
export TSS="/net/seq/data/genomes/human/GRCh38/noalts/ref/refGene.CombinedTxStarts.bed"

JOB_BASENAME=".AGG#${AGGREGATION_ID}"
MERGE_JOBNAME=${JOB_BASENAME}_merge
DUP_JOBNAME=${JOB_BASENAME}_dup
PROCESS_BAM_JOBNAME=${JOB_BASENAME}_pb
HOTSPOT_JOBNAME=${JOB_BASENAME}_hotspot
SPOTSCORE_JOBNAME=${JOB_BASENAME}_spotscore
COUNT_JOBNAME=${JOB_BASENAME}_count
ADAPTERCOUNT_JOBNAME=${JOB_BASENAME}_adaptercount
PRESEQ_JOBNAME=${JOB_BASENAME}_preseq
DENSITY_JOBNAME=${JOB_BASENAME}_density
CUTCOUNTS_JOBNAME=${JOB_BASENAME}_cutcounts
INSERT_JOBNAME=${JOB_BASENAME}_insert

cd $AGGREGATION_FOLDER
BAM_COUNT=`ls $BAM_FILES | wc -l`

# record version
cp $STAMPIPES/version.json .

if [ -n "$REDO_AGGREGATION" ]; then
        bash $STAMPIPES/scripts/bwa/aggregate/basic/reset.bash
fi

# Check out files match first
#python3 $STAMPIPES/scripts/utility/md5check.py bamfiles.txt || exit 1

PROCESSING=""

# merge bams
if [[ ! -s "$FINAL_BAM.bai" ]]; then
	merge_jobid=$(sbatch --export=ALL -J "$MERGE_JOBNAME" -o "$MERGE_JOBNAME.o%A" -e "$MERGE_JOBNAME.e%A" --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: merging BAMs"
date

if [[ $BAM_COUNT -eq 1 ]]; then
	rsync ${BAM_FILES} ${FINAL_BAM}
else
	samtools merge ${FINAL_BAM} ${BAM_FILES}
fi
samtools index ${FINAL_BAM}

echo "FINISH: merging BAMs"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$merge_jobid"
fi

# set final bam dependencies
if [[ -n $merge_jobid ]]; then
   dependencies_merge=$(echo $merge_jobid | sed -e 's/^/--dependency=afterok:/g')
fi
	
# process BAM file
if [[ ! -s "$FINAL_UNIQUES_BAM.bai" ]]; then	
	pb_jobid=$(sbatch --export=ALL -J "$DUP_JOBNAME" -o "$DUP_JOBNAME.o%A" -e "$DUP_JOBNAME.e%A" $dependencies_merge --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: processing BAM"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

if [[ "$UMI" == "True" && -n "$PAIRED" ]]; then
	make -f "$STAMPIPES/makefiles/picard/dups_cigarumi.mk" SAMPLE_NAME="${LIBRARY_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM="${FINAL_BAM_MARKED}"
	mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
	samtools view -b -F 1536 ${FINAL_BAM} > ${FINAL_UNIQUES_BAM}
        samtools view -F 512 -u ${FINAL_BAM} | python3 /home/solexa/stampipes-hpc/scripts/bam/mark_dups.py -o /dev/null --hist "$PRESEQ_HIST"
elif [[ -n "$PAIRED" ]]; then
	make -f "$STAMPIPES/makefiles/picard/dups_cigar.mk" SAMPLE_NAME="${LIBRARY_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM="${FINAL_BAM_MARKED}"
	mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
	samtools view -b -F 512 ${FINAL_BAM} > ${FINAL_UNIQUES_BAM}
        samtools index ${FINAL_UNIQUES_BAM}
        cat ${NUCLEAR_CHR} | xargs samtools view -b ${FINAL_UNIQUES_BAM} > \${TMPDIR}/${FINAL_UNIQUES_BAM}.nuclear.bam
        python3 $STAMPIPES/scripts/bam/mark_dups.py -i \${TMPDIR}/${FINAL_UNIQUES_BAM}.nuclear.bam -o /dev/null --hist "$PRESEQ_HIST"
else
        make -f $STAMPIPES/makefiles/picard/dups.mk SAMPLE_NAME="${LIBRARY_NAME}" BAMFILE="${FINAL_BAM}" OUTBAM=${FINAL_BAM_MARKED}
        mv ${FINAL_BAM_MARKED} ${FINAL_BAM}
        samtools view -b -F 512 ${FINAL_BAM} > ${FINAL_UNIQUES_BAM}
        samtools index ${FINAL_UNIQUES_BAM}
        cat ${NUCLEAR_CHR} | xargs samtools view -b ${FINAL_UNIQUES_BAM} > \${TMPDIR}/${FINAL_UNIQUES_BAM}.nuclear.bam
        python3 $STAMPIPES/scripts/bam/mark_dups.py -i \${TMPDIR}/${FINAL_UNIQUES_BAM}.nuclear.bam -o /dev/null --hist "$PRESEQ_HIST"
fi

samtools index ${FINAL_BAM}
samtools index ${FINAL_UNIQUES_BAM}

rm -rf "\$TMPDIR"

echo "FINISH: processing BAM"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$pb_jobid"
fi

# set final bam dependencies
if [[ -n $pb_jobid ]]; then
   dependencies_pb=$(echo $pb_jobid | sed -e 's/^/--dependency=afterok:/g')
fi

# Run Hotspot2
if [[ ! -s "$HOTSPOT_CALLS" || ! -s "$HOTSPOT_CALLS_001" ]] ; then
	HOTSPOT_SPOT=$HOTSPOT2_DIR/$HOTSPOT_PREFIX.SPOT.txt
	jobid=$(sbatch --export=ALL -J "$HOTSPOT_JOBNAME" -o "$HOTSPOT_JOBNAME.o%A" -e "$HOTSPOT_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

echo "Hostname: "
hostname

echo "START: hotspot2"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

NUCLEAR_MAPPABLE_REGIONS=\$(mktemp)
NUCLEAR_CHROM_SIZES=\$(mktemp)

join <(sort "$NUCLEAR_CHR") "$MAPPABLE_REGIONS" | sed 's/\\s\\+/\\t/g' > "\$NUCLEAR_MAPPABLE_REGIONS"
join <(sort "$NUCLEAR_CHR") "$CHROM_SIZES" | sed 's/\\s\\+/\\t/g' > "\$NUCLEAR_CHROM_SIZES"

"$HOTSPOT_SCRIPT"  -F 0.5 -p varWidth_20_${LIBRARY_NAME} -M "\$NUCLEAR_MAPPABLE_REGIONS" -c "\$NUCLEAR_CHROM_SIZES" -C "$CENTER_SITES" "$FINAL_UNIQUES_BAM"  "$HOTSPOT2_DIR"
"$STAMPIPES/scripts/SPOT/info.sh" "$HOTSPOT_CALLS" hotspot2 \$(cat $HOTSPOT_SPOT) > "$HOTSPOT_PREFIX.hotspot2.info"
hsmerge.sh -f 0.01 $HOTSPOT_ALLCALLS $HOTSPOT_CALLS_01
hsmerge.sh -f 0.001 $HOTSPOT_ALLCALLS $HOTSPOT_CALLS_001

# peaks at FDRs
cleave_count=\$(cat $HOTSPOT_CLEAVAGES)
mv $HOTSPOT_PEAKS_DEF $HOTSPOT_PEAKS_05
mv $HOTSPOT_NPEAKS_DEF $HOTSPOT_NPEAKS_05
density-peaks.bash \$TMPDIR varWidth_20_${LIBRARY_NAME} $HOTSPOT_CUTCOUNTS $HOTSPOT_CALLS_01 $CHROM_SIZES $HOTSPOT_DENSITY $HOTSPOT_PEAKS_01 \$cleave_count
density-peaks.bash \$TMPDIR varWidth_20_${LIBRARY_NAME} $HOTSPOT_CUTCOUNTS $HOTSPOT_CALLS_001 $CHROM_SIZES $HOTSPOT_DENSITY $HOTSPOT_PEAKS_001 \$cleave_count
# create iSPOT
totalcuts=\$(cat ${HOTSPOT_CLEAVAGES})
if [[ -n "$ALTIUS_MASTERLIST" ]]; then
    echo -ne "ispot\t" > $HOTSPOT_ISPOT
    bedops -e 1 ${HOTSPOT_CUTCOUNTS} ${ALTIUS_MASTERLIST} | awk -v total=\$totalcuts '{sum += \$5} END {print sum/total}' >> $HOTSPOT_ISPOT
fi

# create prox/distal estimates
closest-features --dist --delim '\t' --closest $HOTSPOT_CALLS $TSS > \$TMPDIR/closest.txt
cat \$TMPDIR/closest.txt | grep -v "NA$" | awk -F"\t" '{print \$NF}' | sed -e 's/-//g' > \$TMPDIR/closest.clean.txt
echo -ne "percent-proximal-0bp\npercent-proximal-1000bp\npercent-proximal-2500bp\npercent-proximal-5000bp\npercent-proximal-10000bp\n" > \$TMPDIR/row_one.txt
cat \$TMPDIR/closest.clean.txt | awk '{ if (\$1 > 0) sum+= 1 } END {print sum/NR}' > \$TMPDIR/row_two.txt
cat \$TMPDIR/closest.clean.txt | awk '{ if (\$1 > 1000) sum+= 1 } END {print sum/NR}' >> \$TMPDIR/row_two.txt
cat \$TMPDIR/closest.clean.txt | awk '{ if (\$1 > 2500) sum+= 1 } END {print sum/NR}' >> \$TMPDIR/row_two.txt
cat \$TMPDIR/closest.clean.txt | awk '{ if (\$1 > 5000) sum+= 1 } END {print sum/NR}' >> \$TMPDIR/row_two.txt
cat \$TMPDIR/closest.clean.txt | awk '{ if (\$1 > 10000) sum+= 1 } END {print sum/NR}' >> \$TMPDIR/row_two.txt
paste \$TMPDIR/row_one.txt \$TMPDIR/row_two.txt > $HOTSPOT_PROXDIST

# create sparse motifs
bedmap --echo --echo-map-id --fraction-map 1 --delim '\t' $HOTSPOT_CALLS $FIMO_TRANSFAC_1E4 > \$TMPDIR/temp.bedmap.txt
python $STAMPIPES/scripts/bwa/aggregate/basic/sparse_motifs.py $FIMO_NAMES \$TMPDIR/temp.bedmap.txt
# temporary name change here
mv hs_motifs_svmlight.txt $HOTSPOT_PREFIX.hs_motifs_svmlight.txt
mv hs_motifs_svmlight.rows.txt $HOTSPOT_PREFIX.hs_motifs_svmlight.rows.txt
mv hs_motifs_svmlight.cols.txt $HOTSPOT_PREFIX.hs_motifs_svmlight.cols.txt

echo "FINISH: hotspot2"
date

rm -rf "\$TMPDIR"

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# SPOT score
if [[ -n "$PAIRED" && ! -e "$LIBRARY_NAME.$GENOME.R1.rand.uniques.sorted.spot.info" ]] || [[ ! -n "$PAIRED" && ! -e "$LIBRARY_NAME.$GENOME.rand.uniques.sorted.spot.info" ]]; then
	jobid=$(sbatch --export=ALL -J "$SPOTSCORE_JOBNAME" -o "$SPOTSCORE_JOBNAME.o%A" -e "$SPOTSCORE_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: HOTSPOT1 SPOT SCORE"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

if [[ -n "$PAIRED" ]]; then
	make -f $STAMPIPES/makefiles/SPOT/spot-R1-paired.mk BWAINDEX=$GENOME_INDEX ASSAY=$ASSAY GENOME=$GENOME READLENGTH=$READ_LENGTH SAMPLE_NAME="$LIBRARY_NAME.$GENOME"
else
	make -f $STAMPIPES/makefiles/SPOT/spot-single.mk BWAINDEX=$GENOME_INDEX ASSAY=$ASSAY GENOME=$GENOME READLENGTH=$READ_LENGTH SAMPLE_NAME="$LIBRARY_NAME.$GENOME"
fi

rm -rf "\$TMPDIR"

echo "FINISH: HOTSPOT1 SPOT SCORE"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# bam counts
if [ ! -e $TAGCOUNTS_FILE ]; then
	jobid=$(sbatch --export=ALL -J "$COUNT_JOBNAME" -o "$COUNT_JOBNAME.o%A" -e "$COUNT_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=4000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash

set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: tag counts"
date

python3 $STAMPIPES/scripts/bwa/bamcounts.py $FINAL_BAM $TAGCOUNTS_FILE

echo "FINISH: tag counts"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# adapter counts
if [[ ! -s "$ADAPTER_COUNT_FILE" ]]; then
	jobid=$(sbatch --export=ALL -J "$ADAPTERCOUNT_JOBNAME" -o "$ADAPTERCOUNT_JOBNAME.o%A" -e "$ADAPTERCOUNT_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=10000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: adapter counts"
date

adaptercount=\$(bash "$STAMPIPES/scripts/bam/count_adapters.sh" "$FINAL_BAM")
if [ -n \$adapter_count ]; then
	echo -e "adapter\t\$adaptercount" > "$ADAPTER_COUNT_FILE"
fi

echo "FINISH: adapter counts"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# preseq
if [[ ! -s "$PRESEQ_TRGT" ]]; then
        jobid=$(sbatch --export=ALL -J "$PRESEQ_JOBNAME" -o "$PRESEQ_JOBNAME.o%A" -e "$PRESEQ_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=4000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: preseq"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

# get preseq metric
preseq lc_extrap -hist "$PRESEQ_HIST" -extrap 1.001e9 -s 1e6 -v > "$PRESEQ_RES" || echo "NA" > "$PRESEQ_RES"

# write out preseq targets
bash "$STAMPIPES/scripts/utility/preseq_targets.sh" "${PRESEQ_RES}" "${PRESEQ_TRGT}"

rm -rf "\$TMPDIR"

echo "FINISH: preseq"
date

__SCRIPT__
)
        PROCESSING="$PROCESSING,$jobid"
fi

# density tracks
if [[ ! -s "$DENSITY_BIGWIG" || ! -s "$NORM_DENSITY_BIGWIG" ]]; then
	jobid=$(sbatch --export=ALL -J "$DENSITY_JOBNAME" -o "$DENSITY_JOBNAME.o%A" -e "$DENSITY_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: density"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

make -f $STAMPIPES/makefiles/densities/density.mk FAI=${GENOME_INDEX}.fai SAMPLE_NAME=${LIBRARY_NAME} GENOME=${GENOME} \
	BAMFILE=${FINAL_UNIQUES_BAM} STARCH_OUT=${DENSITY_STARCH} BIGWIG_OUT=${DENSITY_BIGWIG}
make -f "$STAMPIPES/makefiles/densities/normalize-density.mk" BAMFILE=${FINAL_UNIQUES_BAM} SAMPLE_NAME=${LIBRARY_NAME} FAI=${GENOME_INDEX}.fai

# tabix
unstarch $NORM_DENSITY_STARCH | bgzip > $NORM_DENSITY_STARCH.bgz
tabix -p bed $NORM_DENSITY_STARCH.bgz
unstarch $DENSITY_STARCH | bgzip > $DENSITY_STARCH.bgz
tabix -p bed $DENSITY_STARCH.bgz

rm -rf "\$TMPDIR"

echo "FINISH: density"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# cutcounts
if [ ! -e "$CUTCOUNTS_BIGWIG" ]; then
	jobid=$(sbatch --export=ALL -J "$CUTCOUNTS_JOBNAME" -o "$CUTCOUNTS_JOBNAME.o%A" -e "$CUTCOUNTS_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=32000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: cutcounts"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

bash $STAMPIPES/scripts/bwa/aggregate/basic/cutcounts.bash

# tabix
unstarch $CUTCOUNTS_STARCH | bgzip > $CUTCOUNTS_STARCH.bgz
tabix -p bed $CUTCOUNTS_STARCH.bgz

rm -rf "\$TMPDIR"

echo "FINISH: cutcounts"
date

__SCRIPT__
)
	PROCESSING="$PROCESSING,$jobid"
fi

# calculate insert sizes
if [[ ! -s "${INSERT_FILE}.info" && -n "$PAIRED" ]]; then
        jobid=$(sbatch --export=ALL -J "$INSERT_JOBNAME" -o "$INSERT_JOBNAME.o%A" -e "$INSERT_JOBNAME.e%A" $dependencies_pb --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=16000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: insert file"
date

export TMPDIR=/tmp/slurm.\$SLURM_JOB_ID
mkdir -p \$TMPDIR

cat ${NUCLEAR_CHR} | xargs samtools view -b ${FINAL_UNIQUES_BAM} > \${TMPDIR}/${FINAL_UNIQUES_BAM}.nuclear.bam
picard CollectInsertSizeMetrics INPUT=\${TMPDIR}/${FINAL_UNIQUES_BAM}.nuclear.bam OUTPUT=${LIBRARY_NAME}.CollectInsertSizeMetrics.picard \
    HISTOGRAM_FILE=${LIBRARY_NAME}.CollectInsertSizeMetrics.picard.pdf \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true && echo Picard stats >&2

cat ${LIBRARY_NAME}.CollectInsertSizeMetrics.picard | awk '/## HISTOGRAM/{x=1;next}x' | sed 1d > \$TMPDIR/temp_hist.txt
python3 $STAMPIPES/scripts/utility/picard_inserts_process.py \$TMPDIR/temp_hist.txt > ${LIBRARY_NAME}.CollectInsertSizeMetrics.picard.info

rm -rf "\$TMPDIR"

echo "FINISH: insert file"
date

__SCRIPT__
)
        PROCESSING="$PROCESSING,$jobid"
fi

# get complete dependencies and run complete even with some failures
complete_dependencies=$(echo $PROCESSING | sed -e 's/,/,afterany:/g' | sed -e 's/^,afterany/--dependency=afterany/g')

# check completion and upload data
if [[ -n "${PROCESSING}" ]]; then
	UPLOAD_SCRIPT=$STAMPIPES/scripts/lims/upload_data.py
	jobid=$(sbatch --export=ALL -J "${JOB_BASENAME}_complete" -o "${JOB_BASENAME}_complete.o%A" -e "${JOB_BASENAME}_complete.e%A" $complete_dependencies --partition=$QUEUE --cpus-per-task=1 --ntasks=1 --mem-per-cpu=1000 --parsable --oversubscribe <<__SCRIPT__
#!/bin/bash
set -x -e -o pipefail

echo "Hostname: "
hostname

echo "START: "
date

bash "$STAMPIPES/scripts/bwa/aggregate/basic/checkcomplete.bash"
bash "$STAMPIPES/scripts/bwa/aggregate/basic/attachfiles.bash"
bash "$STAMPIPES/scripts/bwa/aggregate/basic/uploadcounts.bash"

echo "FINISH: "
date

__SCRIPT__
)
fi
