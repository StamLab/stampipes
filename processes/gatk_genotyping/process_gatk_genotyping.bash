# Dependencies
source "$MODULELOAD"
module purge
module load jdk/11.0.16
module load nextflow
module load python/3.5.1
module load bwa/0.7.17 
module load samtools/1.7

# these may be  need to be defined to work with production runs
#source "$PYTHON3_ACTIVATE"
#source "$STAMPIPES/scripts/sentry/sentry-lib.bash"

FLOWCELL=$1
SAMPLE=$2
LIBRARY_NUMBER=$3
FLOWCELL_PATH=($(echo "/net/seq/data2/flowcells/FC${FLOWCELL}_*"))
SEQ_FILE_DIR=($(echo "$FLOWCELL_PATH/Project_Lab/Sample_${SAMPLE}"))

WORKDIR='./work/'
outDir='./results'


# Run the pipeline
NXF_VER=21.10.6 nextflow \
  #-c $STAMPIPES/nextflow.config \
  run "$STAMPIPES"/processes/genotyping/bwa_gatk.nf \
  -with-trace \
  -ansi-log false \
  -profile docker,cluster \
  -resume \
  -work-dir "$WORKDIR" \
  --inputDir "$SEQ_FILE_DIR"
  --outDir "$outdir"