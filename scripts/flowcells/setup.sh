#!/bin/bash

set -o errexit
set -o pipefail

# Dependencies
source $MODULELOAD
source $PYTHON3_ACTIVATE

#########
# Options
#########

NODES=40

usage(){
cat << EOF
usage: $0 [-f flowcell_label]

Setup analysis for a flowcell

OPTIONS:
-h        Show this message
-v        Verbose
-f        Flowcell Label
-d        Requires by-hand demuxing
EOF
}

verbose=
demux=
while getopts ":hvdf:" opt ; do
    case $opt in
    h)
        usage
        exit 0
        ;;
    v)
        verbose=true
        ;;
    d)
        demux=true
        ;;
    f)
        flowcell="$OPTARG"
        ;;
    :)
        echo "Option -$OPTARG requires an argument" >&2
        usage >&2
        exit 1
        ;;
    \?)
        echo "Invalid option: -$OPTARG" >&2
        usage >&2
        exit 1
        ;;
    esac
done

if [ -z "$flowcell" ] ; then
    echo "No flowcell label specified"
    flowcell=$(basename $(pwd) | cut -f4 -d_ | cut -c2-6)
    echo "Guessing $flowcell..."
    #usage >&2
    #exit 1
fi

#######################
# Samplesheet functions
#######################

make_hiseq_samplesheet(){
  echo "FCID,Lane,SampleID,SampleRef,Index,Description,Control,Recipe,Operator,SampleProject"

  if [ -z "$demux" ] ; then
    # ( X | tostring) syntax is for non-string fields
    # Default values (if field is false or null) come after //
    jq -r --arg flowcell "$flowcell" '
    .libraries as $l
    | $l
    | map( select(.failed == false) )
    | map( .lane as $num
    | .barcode_index =
    if (  $l | map(select( $num  == .lane )) | length == 1 ) then
      "NoIndex"
    else
      .barcode_index
      end )
      | .[] | [
      "FC" + $flowcell,
      (.lane | tostring),
      .samplesheet_name,
      .alignments[0].genome_index // "contam",
      .barcode_index              // "NoIndex",
      .cell_type                  // "None"  ,
      "N",
      .assay                      // "N/A"   ,
      "orders",
      .project
      ] | join(",") ' "$json"
    else
      for i in $(seq 8) ; do
        echo "FC$flowcell,$i,none,none,GGGGGGGG-GGGGGGGG,none,N,none,none,none"
      done
    fi

  }

make_nextseq_samplesheet(){
  name=Stamlab
  date=$(date '+%m/%d/%Y')
  cat <<__SHEET__
[Header]
Investigator Name,$name
Project Name,$name
Experiment Name,$name
Date,$date
Workflow,GenerateFASTQ

[Settings]

[Data]
SampleID,SampleName,index,index2
none,none,GGGGGGGG,GGGGGGGG
__SHEET__

if [ -z "$demux" ] ; then
  # This bit of cryptic magic generates the samplesheet part.
  jq -r '.libraries[] | select(.failed == false) | [.samplesheet_name,.samplesheet_name,.barcode_index,""] | join(",") ' "$json" \
    | sed 's/\([ACTG]\+\)-\([ACTG]\+\),$/\1,\2/'  # Changes dual-index barcodes to proper format
fi

}

########
# Main #
########

json="processing.json"
illumina_dir=$(pwd)

link_command="#no linking to do"

source "$STAMPIPES/scripts/lims/api_functions.sh"
lims_put_by_url "$(lims_get_all "flowcell_run/?label=$flowcell" | jq -r .url)prepare_for_processing/"

# Get and read the processing script
python3 "$STAMPIPES/scripts/lims/get_processing.py" -f "$flowcell" -o "$json"
run_type=$(     jq -r '.flowcell.run_type'          "$json" )
analysis_dir=$( jq -r '.alignment_group.directory'  "$json" )
mask=$(         jq -r '.alignment_group.bases_mask' "$json" )
run_type=$(     jq -r '.flowcell.run_type'          "$json" )
has_umi=$(      jq -r '.libraries | map(.barcode1.umi) | any' "$json")

if [ -z "$demux" ] ; then
  bcl_mask=$mask
  mismatches=$(python3 $STAMPIPES/scripts/flowcells/max_mismatch.py --ignore_failed_lanes)
  if [ "$has_umi" == "true" ] ; then
    echo "---WARNING---"
    echo "Flowcell contains UMI samples, but -d param was not specified"
    echo "You probably want to re-run with -d"
    echo "-------------"
  fi

else # Set some options for manual demultiplexing
  bcl_mask=$(tr Nn Ii <<< $mask)
  mismatches="0,0"
  dmx_mismatches=$(python3 $STAMPIPES/scripts/flowcells/max_mismatch.py --ignore_failed_lanes | cut -c1 )
fi

case $run_type in
"NextSeq 500")
    echo "Regular NextSeq 500 run detected"
    parallel_env="-pe threads 6"
    link_command="python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o ."
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--nextseq"
    queue="queue2"
    make_nextseq_samplesheet > SampleSheet.csv
    bcl_tasks=1

    # The quadruple-backslash syntax on this is messy and gross.
    # It works, though, and the output is readable.
    # read -d '' always exits with status 1, so we ignore error

    # The NSLOTS lines are for scaling the various threads (2 per slot).
    # WARNING: Does not work for threads < 4
    # Table:
    # NSLOTS  l w d p   total
    # 4       1 1 2 4 = 8
    # 5       1 1 2 5 = 9
    # 6       2 2 3 6 = 13
    # 7       2 2 3 7 = 14
    # 8       2 2 4 8 = 16
    set +e
    read -d '' unaligned_command  << _U_
    bcl2fastq \\\\
      --input-dir "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
      --use-bases-mask "$bcl_mask" \\\\
      --output-dir "$fastq_dir" \\\\
      --barcode-mismatches "$mismatches" \\\\
      --loading-threads        \\\$(( NSLOTS / 4 )) \\\\
      --writing-threads        \\\$(( NSLOTS / 4 )) \\\\
      --demultiplexing-threads \\\$(( NSLOTS / 2 )) \\\\
      --processing-threads     \\\$(( NSLOTS ))
_U_
    set -e
    ;;
"HiSeq 4000")
    echo "Hiseq 4000 run detected"
    parallel_env="-pe threads 6"
    link_command="python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o ."
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--hiseq4k"
    queue="queue2"
    make_nextseq_samplesheet > SampleSheet.csv
    bcl_tasks=1-8

    # The quadruple-backslash syntax on this is messy and gross.
    # It works, though, and the output is readable.
    # read -d '' always exits with status 1, so we ignore error

    # The NSLOTS lines are for scaling the various threads (2 per slot).
    # WARNING: Does not work for threads < 4
    # Table:
    # NSLOTS  l w d p   total
    # 4       1 1 2 4 = 8
    # 5       1 1 2 5 = 9
    # 6       2 2 3 6 = 13
    # 7       2 2 3 7 = 14
    # 8       2 2 4 8 = 16
    set +e
    read -d '' unaligned_command  << _U_
    bcl2fastq \\\\
      --input-dir "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
      --use-bases-mask "$bcl_mask" \\\\
      --output-dir "$fastq_dir.L00\\\$SGE_TASK_ID" \\\\
      --barcode-mismatches "$mismatches" \\\\
      --loading-threads        \\\$(( NSLOTS / 4 )) \\\\
      --writing-threads        \\\$(( NSLOTS / 4 )) \\\\
      --demultiplexing-threads \\\$(( NSLOTS / 2 )) \\\\
      --processing-threads     \\\$(( NSLOTS ))     \\\\
      --tiles s_\\\$SGE_TASK_ID
_U_
  ;;
"MiniSeq High Output Kit")
    # Identical to nextseq processing
    echo "High-output MiniSeq run detected"
    parallel_env="-pe threads 6"
    link_command="python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o ."
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--miniseq"
    queue="queue2"
    make_nextseq_samplesheet > SampleSheet.csv
    bcl_tasks=1
    set +e
    read -d '' unaligned_command  << _U_
    bcl2fastq \\\\
      --input-dir "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
      --use-bases-mask "$bcl_mask" \\\\
      --output-dir "$fastq_dir" \\\\
      --barcode-mismatches "$mismatches" \\\\
      --loading-threads        \\\$(( NSLOTS / 4 )) \\\\
      --writing-threads        \\\$(( NSLOTS / 4 )) \\\\
      --demultiplexing-threads \\\$(( NSLOTS / 2 )) \\\\
      --processing-threads     \\\$(( NSLOTS ))
_U_
    set -e
    ;;
    #TODO: Add HISEQ V3 on hiseq 2500 (rapid run mode)
"HISEQ V4")
    echo "Regular HiSeq 2500 run detected"
    echo "HiSeq 2500 processing not supported on the new cluster! (Does not have old version of bcl2fastq)"
    exit 2
    #parallel_env=""
    #link_command='#no linking to do'
    #samplesheet=$(pwd)/Data/Intensities/BaseCalls/SampleSheet.csv
    #mkdir -p $(dirname "$samplesheet")
    #make_hiseq_samplesheet > "$samplesheet"
    #fastq_dir="$illumina_dir/Unaligned/"  # Trailing slash is important for rsync!
    #bc_flag="--hiseq"
    #bcl_tasks=1

    #set +e
    #read -d '' unaligned_command <<_U_
    #if [ ! -e "$fastq_dir" ] ; then
    #        configureBclToFastq.pl \\\\
    #          --mismatches "$mismatches" \\\\
    #          --output-dir "$fastq_dir" \\\\
    #          --fastq-cluster-count 16000000 \\\\
    #          --with-failed-reads --sample-sheet $samplesheet \\\\
    #          --use-bases-mask "$bcl_mask"  \\\\
    #          --input-dir "$illumina_dir/Data/Intensities/BaseCalls"
    #fi

    #cd "$fastq_dir"
    #qmake -now no -cwd -q all.q -V -- -j "$NODES"
#_U_
    #set -e
    ;;
*)
    echo "Unrecognized run type '$run_type'"
    exit 1
    ;;
esac

copy_from_dir="$fastq_dir"
if [ -n "$demux" ] ; then
  copy_from_dir="$(pwd)/Demultiplexed/"
  demux_cmd="$STAMPIPES/scripts/flowcells/demux_flowcell.sh -i $fastq_dir -o $copy_from_dir -p $json -q $queue -m $dmx_mismatches"
  link_command="#Demuxing happened, no linking to do"
fi


flowcell_id=$( curl \
  "$LIMS_API_URL"/flowcell_run/?label=$flowcell \
  -H "Authorization: Token $LIMS_API_TOKEN" \
  -k \
  2>/dev/null \
  | jq '.results[] | .id'
)

# The final script is below:
cat > run_bcl2fastq.sh <<__BCL2FASTQ__
#!/bin/bash

source $MODULELOAD
module load bcl2fastq2/2.17.1.14
source $PYTHON3_ACTIVATE

source $STAMPIPES/scripts/lims/api_functions.sh

# Register the file directory
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  --attach_directory "$analysis_dir" \
  --attach_file_contenttype sequencingdata.flowcellrun \
  --attach_file_purpose flowcell-directory \
  --attach_file_objectid $flowcell_id

# register as "Sequencing" in LIMS
lims_patch "flowcell_run/$flowcell_id/" "status=https://lims.stamlab.org/api/flowcell_run_status/2/"

while [ ! -e "$illumina_dir/RTAComplete.txt" ] ; do sleep 60 ; done

# register as "Processing" in LIMS
lims_patch "flowcell_run/$flowcell_id/" "status=https://lims.stamlab.org/api/flowcell_run_status/3/"
lims_patch "flowcell_run/$flowcell_id/" "folder_name=${PWD##*/}"

# Submit a barcode job for each mask
for bcmask in $(python $STAMPIPES/scripts/flowcells/barcode_masks.py | xargs) ; do
    export bcmask
    qsub -cwd -N "bc-$flowcell" -q $queue $parallel_env -V -S /bin/bash <<'__BARCODES__'
    GOMAXPROCS=\$(( NSLOTS * 2 )) bcl_barcode_count --mask=\$bcmask $bc_flag > barcodes.\$bcmask.json

    python3 $STAMPIPES/scripts/lims/upload_data.py --barcode_report barcodes.\$bcmask.json
__BARCODES__
done

qsub -cwd -N "u-$flowcell" -q $queue $parallel_env -V -t $bcl_tasks -S /bin/bash  <<'__FASTQ__'

set -x -e -o pipefail

cd "$illumina_dir"

$unaligned_command

__FASTQ__

qsub -cwd -N "dmx-$flowcell" -q $queue -hold_jid "u-$flowcell" -V -S /bin/bash <<__DMX__
  $demux_cmd
  # Wait for jobs to finish
  while ( squeue -o "%j" | grep -q '^.dmx') ; do
    sleep 60
  done
__DMX__

qsub -cwd -N "c-$flowcell" -q $queue -hold_jid "dmx-$flowcell" -V -S /bin/bash <<__COPY__
mkdir -p "$analysis_dir"

rsync -avP "$copy_from_dir" "$analysis_dir/"
rsync -avP "$illumina_dir/InterOp" "$analysis_dir/"
rsync -avP "$illumina_dir/RunInfo.xml" "$analysis_dir/"
rsync -avP "$samplesheet" "$analysis_dir"

cd "$analysis_dir"

$link_command

# Remove existing scripts if they exist (to avoid appending)
rm -f fastqc.bash collate.bash run.bash

# Create fastqc scripts
python3 "$STAMPIPES/scripts/laneprocess.py" \
  --script_template "$STAMPIPES/processes/fastq/fastqc.bash" \
  --qsub-prefix .fq \
  --queue $queue \
  --sample-script-basename fastqc.bash \
  --flowcell_label "$flowcell" \
  --outfile fastqc.bash

# Create collation scripts
python3 "$STAMPIPES/scripts/laneprocess.py" \
  --script_template "$STAMPIPES/processes/fastq/collate_fastq.bash" \
  --qsub-prefix .collatefq \
  --queue $queue \
  --sample-script-basename "collate.bash" \
  --flowcell_label "$flowcell" \
  --outfile collate.bash

bash collate.bash

qsub -N .run$flowcell -q $queue -hold_jid '.collatefq*$flowcell*' -cwd -V -S /bin/bash <<__ALIGN__
  bash fastqc.bash

  # Create alignment scripts
  python3 $STAMPIPES/scripts/alignprocess.py \
    --flowcell $flowcell \
    --outfile run.bash

  bash run.bash
__ALIGN__

__COPY__

__BCL2FASTQ__

if [ -e "RTAComplete.txt" ] ; then
    echo -e "Setup complete. To kick everything off, type:\n\nbash run_bcl2fastq.sh"
else
    echo -e "Setup complete, sequencing still in progress. To queue everything up, type:\n\nnohup bash run_bcl2fastq.sh &"
fi
