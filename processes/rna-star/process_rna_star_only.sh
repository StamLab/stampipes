export VERSION=1.2-alpha
export OUT_DIR=output_$VERSION
mkdir -p "$OUT_DIR"

# Check for special UMTs
if [[ "$LIBRARY_KIT" == "Smarter Stranded Total v3 Pico RNASeq with RNA Isolation" ]] ; then
  UMI=True
  UMI_METHOD=takara-umt
fi

ADAPTER_FLAG=--use_fastp
if [[ "$LIBRARY_KIT" == "Agilent Sure Select XT HS2" ]] ; then
  ADAPTER_FLAG=--use_agent
fi


source "$MODULELOAD"
module load openssl-dev jdk nextflow
source "$PYTHON3_ACTIVATE"

cd "$(dirname "$0")"
VERSION_FILE="${OUT_DIR}/${SAMPLE_NAME}.versions.txt"

# $STAR_DIR is set by the process template, and are relative to the reference directory
REFDIR=$(dirname "$BWAINDEX")
export STARdir="$REFDIR/$STAR_DIR"
bash "$STAMPIPES/scripts/versions.bash" &> "$VERSION_FILE"

# Let LIMS know the alignment is starting
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  -a "$LIMS_API_URL"             \
  -t "$LIMS_API_TOKEN"           \
  --alignment_id "$ALIGNMENT_ID" \
  --start_alignment_progress     \
  --version_file "$VERSION_FILE"

# Run the code
NXF_VER=21.04.1 nextflow run \
  "$STAMPIPES/processes/rna-star/star_alignment.nf" \
  --r1 "$R1_FASTQ" --r2 "$R2_FASTQ" \
  --starIndexDir "$STARdir" \
  --readlength "$READLENGTH" \
  --umimethod "$UMI_METHOD" \
  --outdir "$OUT_DIR" \
  --adapter_p5 "$ADAPTER_P5" \
  --adapter_p7 "$ADAPTER_P7" \
  "$ADAPTER_FLAG" \
  -profile modules,cluster \
  -resume \
  "$@"

# Check for completeness and upload files.
set -e
bash "$STAMPIPES/scripts/rna-star/checkcomplete.bash"
bash "$STAMPIPES/scripts/rna-star/attachfiles.sh"

# Signal all-complete
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  -a "$LIMS_API_URL"             \
  -t "$LIMS_API_TOKEN"           \
  --alignment_id "$ALIGNMENT_ID" \
  --finish_alignment
