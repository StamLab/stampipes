#!/bin/bash
# shellcheck disable=SC1090
# shellcheck disable=SC2162

set -o errexit
set -o pipefail

# Dependencies
source "$MODULELOAD"
source "$PYTHON3_ACTIVATE"

source "$STAMPIPES/scripts/sentry/sentry-lib.bash"

#########
# Options
#########

usage(){
cat << EOF
usage: $0 [-f flowcell_label]

Setup analysis for a flowcell

OPTIONS:
-h        Show this message
-v        Verbose
-f        Flowcell Label
-d        Requires by-hand demuxing
-x        No sleep (running by hand)
EOF
}

verbose=
demux=
nosleep=
while getopts ":hvdxf:" opt ; do
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
    x)
        nosleep=true
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
    flowcell=$(basename "$PWD" | cut -f4 -d_ | cut -c2-10)
    echo "Guessing $flowcell..."
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

# placeholder
make_miniseq_samplesheet(){
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

# placeholder
make_miniseq_guideseq_samplesheet(){
sleep 10
}

########
# Main #
########

json="processing.json"
illumina_dir=$(pwd)

link_command="#no linking to do"

source "$STAMPIPES/scripts/lims/api_functions.sh"
(
  set -e
  url=$(lims_get_all "flowcell_run/?label=$flowcell" | jq -r .url)
  lims_put_by_url "${url}prepare_for_processing/"
) || (
  # Make sure that "Prepare for Processing" has completed.
  if [[ -z "$nosleep" ]] ; then
    echo "Prepare for processing is taking a while to complete..."
    echo "sleeping for 5 minutes to wait for LIMS to set up... (skip with -x if you're sure it's ready)"
    sleep 300
  fi
)

# Get and read the processing script
python3 "$STAMPIPES/scripts/lims/get_processing.py" -f "$flowcell" -o "$json"
run_type=$(     jq -r '.flowcell.run_type'          "$json" )
analysis_dir=$( jq -r '.alignment_group.directory'  "$json" )
mask=$(         jq -r '.alignment_group.bases_mask' "$json" )
run_type=$(     jq -r '.flowcell.run_type'          "$json" )
has_umi=$(      jq -r '.libraries | map(.barcode1.umi) | any' "$json")

# Novaseq runs always use native bcl2fastq demuxing
if [[ $run_type =~ Novaseq ]] || [[ $run_type =~ NovaSeq ]] ; then
  unset demux
fi

# Check if read1length=0 -> that means altseq
# Handle specially
# TODO: Check this from processing.json
flowcell_data=$(lims_get_all "flowcell_run/?label=$flowcell")
read1length=$(echo $flowcell_data | jq -r .read1_length | head -n1)
if [[  "$read1length" = "0" ]] ; then
  echo "Alt-seq run detected"
  mkdir -p "$analysis_dir"
  cp processing.json "$analysis_dir/"
  runscript="$analysis_dir/run.bash"
  (
    echo "#!/bin/bash"
    echo "export FLOWCELL=$flowcell"
    echo "export STAMPIPES=$STAMPIPES"
    cat "$STAMPIPES"/processes/altseq/process_altseq.bash
  ) > "$runscript"

  # Create wrapper for cronjob to call
  cat > run_bcl2fastq.sh <<__BCL2FASTQ__
#!/bin/bash
sbatch --cpus 1 \
  --mem '4G'  \
  --partition hpcz-2 \
  --job-name "altseq-$flowcell-supervisor" <<EOF
#!/bin/bash
cd "$analysis_dir"
bash "$runscript"
EOF
__BCL2FASTQ__
  echo "Run $runscript to start analysis!"

  exit 0
fi

if [ -z "$demux" ] ; then
  bcl_mask=$mask
  mismatches=$(python3 $STAMPIPES/scripts/flowcells/max_mismatch.py --ignore_failed_lanes --allow_collisions)
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

# Long command definitions
# The quadruple-backslash syntax on this is messy and gross.
# It works, though, and the output is readable.
# read -d '' always exits with status 1, so we ignore error
# We split threads equally between processing and loading+writing.
set +e
read -d '' regular_bcl_command  << _REG_BCL_CMD_
    PATH=/home/nelsonjs/src/bcl2fastq2/bin/:\$PATH
    bcl2fastq \\\\
      --input-dir "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
      --use-bases-mask "$bcl_mask" \\\\
      --output-dir "$fastq_dir" \\\\
      --barcode-mismatches "$mismatches" \\\\
      --writing-threads        0                       \\\\
      --loading-threads        \\\$SLURM_CPUS_PER_TASK \\\\
      --processing-threads     \\\$SLURM_CPUS_PER_TASK
_REG_BCL_CMD_

read -d '' novaseq_bcl_command  << _NOVA_BCL_CMD_
    PATH=/home/nelsonjs/src/bcl2fastq2/bin/:\$PATH
    for samplesheet in SampleSheet.withmask*csv ; do
      bcl_mask=\$(sed 's/.*withmask\\.//;s/\\.csv//' <<< \$samplesheet)
      fastq_dir=\$(sed 's/,/-/g' <<< "fastq-withmask-\$bcl_mask")
      bcl2fastq \\\\
        --input-dir          "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
        --output-dir         "${illumina_dir}/\$fastq_dir"                \\\\
        --use-bases-mask     "\$bcl_mask"                                 \\\\
        --barcode-mismatches "$mismatches"                                \\\\
        --sample-sheet       "${illumina_dir}/\$samplesheet"              \\\\
        --writing-threads    0                                            \\\\
        --loading-threads    \\\$SLURM_CPUS_PER_TASK                      \\\\
        --processing-threads \\\$SLURM_CPUS_PER_TASK
    done
_NOVA_BCL_CMD_

# TODO: Remove hardcoded queue here!
# The issue that that 'queue' isn't set until later in the script, but is needed for NOVA_SUBMIT_CMD
queue=hpcz-2

# This is a variant where we submit one job for each lane
read -d '' novaseq_submit_command <<_NOVA_SUBMIT_CMD_
# Run bcl2fastq in parallel, for each samplesheet and lane
PROCESSING=
for samplesheet in SampleSheet.withmask*csv ; do
  for lane in {1..8} ; do
    # Skip submission if lane not in this samplesheet
    if ! (cut -d, -f1 \$samplesheet | sort -u | grep -q \$lane) ; then
      echo "Lane \$lane not in samplesheet \$samplesheet, skipping"
      continue
    fi
    bcl_mask=\$(sed 's/.*withmask\\.//;s/\\.csv//' <<< \$samplesheet)
    fastq_dir=\$(sed 's/,/-/g' <<< "fastq-withmask-\$bcl_mask-lane-00\$lane")
    jobname=u-$flowcell-\$bcl_mask-L00\$lane
    bcl_jobid=\$(sbatch --export=ALL -J "\$jobname" -o "\$jobname.o%A" -e "\$jobname.e%A" --partition=$queue --ntasks=1 --cpus-per-task=40 --mem-per-cpu=8000 --parsable --oversubscribe <<__FASTQ__
#!/bin/bash
      set -x -e -o pipefail
      cd "${illumina_dir}"
      PATH=/home/nelsonjs/src/bcl2fastq2/bin/:\$PATH
      bcl2fastq \\\\
        --input-dir          "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
        --output-dir         "${illumina_dir}/\\\$fastq_dir"              \\\\
        --use-bases-mask     "\\\$bcl_mask"                               \\\\
        --tiles              "s_\\\$lane"                                 \\\\
        --barcode-mismatches "$mismatches"                                \\\\
        --sample-sheet       "${illumina_dir}/\$samplesheet"              \\\\
        --writing-threads    0                                            \\\\
        --loading-threads    \\\\\$SLURM_CPUS_PER_TASK                    \\\\
        --processing-threads \\\\\$SLURM_CPUS_PER_TASK
__FASTQ__
)
    PROCESSING="\$PROCESSING,\$bcl_jobid"
  done
done
if [[ -n "\$PROCESSING" ]]; then
  bcl_dependency=\$(echo \$PROCESSING | sed -e 's/,/,afterok:/g' | sed -e 's/^,afterok/--dependency=afterok/g')
fi

_NOVA_SUBMIT_CMD_

read -d '' novaseq_link_command  <<'_NOVA_LINK_CMD_'
for fq_dir in fastq-withmask-* ; do
  [[ -d $fq_dir ]] || continue
  python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i "$fq_dir" -o Demultiplexed -p processing.json
done
_NOVA_LINK_CMD_
set -e

if [ -z "$flowcell" ] ; then
    echo "No flowcell label specified"
    flowcell=$(basename "$PWD" | cut -f4 -d_ | cut -c2-10)
    echo "Guessing $flowcell..."
fi


case $run_type in

"Novaseq 6000 S1")
    echo "Novaseq 6000: S1 (non-pooled)"
    unset demux
    parallel_env="-pe threads 6"
    link_command=$novaseq_link_command
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--novaseq"
    queue="hpcz-2"
    python "$STAMPIPES/scripts/flowcells/make_samplesheets.py" --reverse_barcode1 -p processing.json
    bcl_tasks=1
    #unaligned_command=$novaseq_bcl_command
    submit_bcl2fastq_cmd=$novaseq_submit_command
;;

"Novaseq 6000 S2")
    echo "Novaseq 6000: S2 (non-pooled)"
    unset demux
    parallel_env="-pe threads 6"
    link_command=$novaseq_link_command
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--novaseq"
    queue="hpcz-2"
    python "$STAMPIPES/scripts/flowcells/make_samplesheets.py" --reverse_barcode1 -p processing.json
    bcl_tasks=1
    #unaligned_command=$novaseq_bcl_command
    submit_bcl2fastq_cmd=$novaseq_submit_command

;;
"Novaseq 6000 S4")
    echo "Novaseq 6000: S4 (non-pooled)"
    unset demux
    parallel_env="-pe threads 6"
    link_command=$novaseq_link_command
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--novaseq"
    queue="hpcz-2"
    python "$STAMPIPES/scripts/flowcells/make_samplesheets.py" --reverse_barcode1 -p processing.json
    bcl_tasks=1
    #unaligned_command=$novaseq_bcl_command
    submit_bcl2fastq_cmd=$novaseq_submit_command

;;
"NovaSeq X 1.5B")
    echo "NovaSeq X: 1.5B"
    unset demux
    parallel_env="-pe threads 6"
    link_command=$novaseq_link_command
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--novaseq"
    queue="hpcz-2"
    python "$STAMPIPES/scripts/flowcells/make_samplesheets.py" --reverse_barcode1 -p processing.json
    bcl_tasks=1
    #unaligned_command=$novaseq_bcl_command
    submit_bcl2fastq_cmd=$novaseq_submit_command

;;
"NovaSeq X 10B")
    echo "NovaSeq X: 10B"
    unset demux
    parallel_env="-pe threads 6"
    link_command=$novaseq_link_command
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--novaseq"
    queue="hpcz-2"
    python "$STAMPIPES/scripts/flowcells/make_samplesheets.py" --reverse_barcode1 -p processing.json
    bcl_tasks=1
    #unaligned_command=$novaseq_bcl_command
    submit_bcl2fastq_cmd=$novaseq_submit_command

;;
"NovaSeq X 25B")
    echo "NovaSeq X: 25B"
    unset demux
    parallel_env="-pe threads 6"
    link_command=$novaseq_link_command
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--novaseq"
    queue="hpcz-2"
    python "$STAMPIPES/scripts/flowcells/make_samplesheets.py" --reverse_barcode1 -p processing.json
    bcl_tasks=1
    #unaligned_command=$novaseq_bcl_command
    submit_bcl2fastq_cmd=$novaseq_submit_command

;;

"Novaseq 6000 SP")
    echo "Novaseq 6000: SP (non-pooled)"
    unset demux
    parallel_env="-pe threads 6"
    link_command=$novaseq_link_command
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--novaseq"
    queue="hpcz-2"
    python "$STAMPIPES/scripts/flowcells/make_samplesheets.py" --reverse_barcode1 -p processing.json
    bcl_tasks=1
    unaligned_command=$novaseq_bcl_command

;;
"NextSeq 500")

    echo "Regular NextSeq 500 run detected"
    parallel_env="-pe threads 6"
    link_command="python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o . --merge-across-lanes"
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--nextseq"
    queue="queue0"
    make_nextseq_samplesheet > SampleSheet.csv
    bcl_tasks=1
    unaligned_command=$regular_bcl_command
    ;;
"HiSeq 4000")
    echo "Hiseq 4000 run detected"
    parallel_env="-pe threads 6"
    link_command="python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o ."
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--hiseq4k"
    queue="queue0"
    make_nextseq_samplesheet > SampleSheet.csv
    bcl_tasks=1-8
    unaligned_command=$regular_bcl_command
  ;;
"MiniSeq High Output Kit DNase")
    # Identical to nextseq processing
    echo "High-output MiniSeq run detected for DNase"
    parallel_env="-pe threads 6"
    link_command="python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o . --merge-across-lanes"
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--miniseq"
    queue="queue0"
    make_nextseq_samplesheet > SampleSheet.csv
    bcl_tasks=1
    unaligned_command=$regular_bcl_command
    ;;
"MiniSeq Mid Output Kit GUIDEseq")
    # Identical to nextseq processing
    echo "Mid-output MiniSeq run detected for GUIDEseq"
    parallel_env="-pe threads 6"
    link_command="python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o . --merge-across-lanes"
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--miniseq"
    queue="queue0"
    minidemux="True"
    # placeholder
    cp /home/dchee7/projects/guide-seq/data/samplesheets/SampleSheet.csv SampleSheet.csv
    bcl_tasks=1
    set +e
    read -d '' unaligned_command  << _U_
    bcl2fastq \\\\
      --input-dir "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
      --output-dir "$fastq_dir" \\\\
      --create-fastq-for-index-reads
_U_
    set -e
    ;;
"MiniSeq Mid Output Kit")
    # Identical to nextseq processing
    echo "Mid-output MiniSeq run detected"
    parallel_env="-pe threads 6"
    link_command="python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o . --merge-across-lanes"
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--miniseq"
    queue="queue0"
    minidemux="True"
    make_miniseq_samplesheet > SampleSheet.csv
    bcl_tasks=1
    set +e
    read -d '' unaligned_command  << _U_
    bcl2fastq \\\\
      --input-dir "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
      --output-dir "$fastq_dir" \\\\
      --no-lane-splitting
_U_
    set -e
    ;;
"MiniSeq High Output Kit")
    # Identical to nextseq processing
    echo "High-output MiniSeq run detected"
    parallel_env="-pe threads 6"
    link_command="python3 $STAMPIPES/scripts/flowcells/link_nextseq.py -i fastq -o . --merge-across-lanes"
    samplesheet="SampleSheet.csv"
    fastq_dir="$illumina_dir/fastq"  # Lack of trailing slash is important for rsync!
    bc_flag="--miniseq"
    queue="queue0"
    minidemux="True"
    # placeholder
    cat /net/fileserv0/projects/vol2/dchee7/datastore/talens/sample_sheets/SampleSheet.csv > SampleSheet.csv
    #make_nextseq_samplesheet > SampleSheet.csv
    bcl_tasks=1
    set +e
    read -d '' unaligned_command  << _U_
    bcl2fastq \\\\
      --input-dir "${illumina_dir}/Data/Intensities/BaseCalls" \\\\
      --output-dir "$fastq_dir" \\\\
      --no-lane-splitting
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
  # obsolete now?
  demux_cmd="$STAMPIPES/scripts/flowcells/demux_flowcell.sh -i $fastq_dir -o $copy_from_dir -p $json -q $queue -m $dmx_mismatches"
  link_command="#Demuxing happened, no linking to do"
elif [[ "$bc_flag" == "--novaseq" ]] ; then
  copy_from_dir="$(pwd)/Demultiplexed/"
fi

flowcell_id=$( curl \
  "$LIMS_API_URL"/flowcell_run/?label=$flowcell \
  -H "Authorization: Token $LIMS_API_TOKEN" \
  -k \
  2>/dev/null \
  | jq '.results[] | .id'
)

# The final script is below:
if [[ -n "$minidemux" ]]; then
  # If miniseq, demux flowcell with fixed/known barcodes.
cat > run_bcl2fastq.sh <<__BCL2FASTQ__
#!/bin/bash
source "$STAMPIPES/scripts/sentry/sentry-lib.bash"

source $MODULELOAD
module load bcl2fastq2/2.17.1.14
source $PYTHON3_ACTIVATE
source $STAMPIPES/scripts/lims/api_functions.sh

# Register the file directory
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  --attach_directory "$analysis_dir" \
  --attach_file_contenttype SequencingData.flowcellrun \
  --attach_file_purpose flowcell-directory \
  --attach_file_objectid $flowcell_id

# Register as "Sequencing" in LIMS
lims_patch "flowcell_run/$flowcell_id/" "status=https://lims.stamlab.org/api/flowcell_run_status/2/"

# Wait for CopyComplete
while [ ! -e "$illumina_dir/CopyComplete.txt" ] ; do sleep 60 ; done

# Register as "Processing" in LIMS
lims_patch "flowcell_run/$flowcell_id/" "status=https://lims.stamlab.org/api/flowcell_run_status/3/"
lims_patch "flowcell_run/$flowcell_id/" "folder_name=${PWD##*/}"

# bcl2fastq
bcl_jobid=\$(sbatch --export=ALL -J "u-$flowcell" -o "u-$flowcell.o%A" -e "u-$flowcell.e%A" \$dependencies_barcodes --partition=$queue --ntasks=1 --cpus-per-task=20 --mem-per-cpu=8000 --parsable --oversubscribe <<'__FASTQ__'
#!/bin/bash

set -x -e -o pipefail
cd "$illumina_dir"

$unaligned_command

# if the run is for GUIDEseq, swap the indexes
if cat processing.json | grep -q "MiniSeq Mid Output Kit GUIDEseq"; then
    zcat fastq/Undetermined_S0_L001_I2_001.fastq.gz | awk '{if(NR % 4 == 2) {x=(substr(\$0,9,16)); y=(substr(\$0,0,8)); print x y; } else print; }' > fastq/Undetermined_S0_L001_I2_001.rev.fastq
    gzip fastq/Undetermined_S0_L001_I2_001.rev.fastq
fi

__FASTQ__
)

__BCL2FASTQ__

else # If not miniseq

# Default (slow) bcl2fastq cmd
if [[ -z "$submit_bcl2fastq_cmd" ]] ; then
  # If we haven't created the submit command yet, wrap up the unaligned_command
  # This is the "old" way of submitting one job that does the whole flowcell
  submit_bcl2fastq_cmd=<<__SUBMIT_BCL2FASTQ_CMD__
# bcl2fastq
bcl_jobid=\$(sbatch --export=ALL -J "u-$flowcell" -o "u-$flowcell.o%A" -e "u-$flowcell.e%A"  --partition=$queue --ntasks=1 --cpus-per-task=20 --mem-per-cpu=8000 --parsable --oversubscribe <<'__FASTQ__'
#!/bin/bash
set -x -e -o pipefail
cd "$illumina_dir"

$unaligned_command
__FASTQ__
)
# Wait for bcl2fastq to complete
if [[ -n \$bcl_jobid ]]; then
   bcl_dependency=\$(echo \$bcl_jobid | sed -e 's/^/--dependency=afterok:/g')
fi
__SUBMIT_BCL2FASTQ_CMD__
fi

# Not miniseq
cat > run_bcl2fastq.sh <<__BCL2FASTQ__
#!/bin/bash

source $MODULELOAD
module load bcl2fastq2/2.20.0.422
source $PYTHON3_ACTIVATE
source $STAMPIPES/scripts/lims/api_functions.sh

# Register the file directory
python3 "$STAMPIPES/scripts/lims/upload_data.py" \
  --attach_directory "$analysis_dir" \
  --attach_file_contenttype SequencingData.flowcellrun \
  --attach_file_purpose flowcell-directory \
  --attach_file_objectid $flowcell_id

# Register as "Sequencing" in LIMS
lims_patch "flowcell_run/$flowcell_id/" "status=https://lims.stamlab.org/api/flowcell_run_status/2/"

# Wait for CopyComplete
while [ ! -e "$illumina_dir/CopyComplete.txt" ] ; do sleep 60 ; done

# Register as "Processing" in LIMS
lims_patch "flowcell_run/$flowcell_id/" "status=https://lims.stamlab.org/api/flowcell_run_status/3/"
lims_patch "flowcell_run/$flowcell_id/" "folder_name=${PWD##*/}"

# Submit a barcode job for each mask
for bcmask in $(python $STAMPIPES/scripts/flowcells/barcode_masks.py | xargs) ; do
    export bcmask
    bcjobid=\$(sbatch --export=ALL -J "bc-$flowcell" -o "bc-$flowcell.o%A" -e "bc-$flowcell.e%A" --partition=$queue --cpus-per-task=10 --ntasks=1 --mem-per-cpu=6400 --parsable --oversubscribe --mail-type=FAIL --mail-user=sequencing@altius.org <<'__BARCODES__'
#!/bin/bash
bcl_barcode_count --mask=\$bcmask $bc_flag > barcodes.\$bcmask.json
python3 $STAMPIPES/scripts/lims/upload_data.py --barcode_report barcodes.\$bcmask.json
bctest=\$(python $STAMPIPES/scripts/flowcells/barcode_check.py --barcodes barcodes.\$bcmask.json --processing processing.json --bcmask \$bcmask)
if [ \$bctest = "FALSE" ];
then
    exit 1
fi

__BARCODES__
)
    PROCESSING="\$PROCESSING,\$bcjobid"
done

$submit_bcl2fastq_cmd

sbatch --export=ALL -J queuedemux-$flowcell -o "queuedemux-$flowcell.o%A" -e "queuedemux-$flowcell.e%A" \$bcl_dependency --partition $queue --ntasks=1 --cpus-per-task=1 --mem-per-cpu=1000 --parsable --oversubscribe <<__PART2__
#!/bin/bash
bash run_bcl2fastq_2.sh
__PART2__

__BCL2FASTQ__

cat > run_bcl2fastq_2.sh <<__BCL2FASTQ2__
# !/bin/bash
source "$STAMPIPES/scripts/sentry/sentry-lib.bash"
source "$MODULELOAD"
source "$PYTHON3_ACTIVATE"
source "$STAMPIPES/scripts/lims/api_functions.sh"

if [[ -n "$demux" ]] ; then
# demultiplex
if [ -d "$fastq_dir.L001" ] ; then
  inputfiles=(\$(find $fastq_dir.L00[1-9] -name "*Undetermined_*fastq.gz" -size +0 ))
else
  inputfiles=(\$(find $fastq_dir          -name "*Undetermined_*fastq.gz"))
fi
run_type=\$(jq -r .flowcell.run_type < "$json")

for i in "\${inputfiles[@]}" ; do
   lane=\$(sed 's/.*_L\(00[0-9]\)_.*/\1/' <(basename "\$i" ))
   if [[ "\$run_type" == "NextSeq 500" ]] ; then
      suffix="--autosuffix"
   else
      suffix="--suffix \$(sed 's/.*_L00[0-9]\(_R[12]_.*\).fastq.gz/\1/' <(basename "\$i" ))"
   fi

   jobid=\$(sbatch --export=ALL -J dmx\$(basename "\$i") -o .dmx\$(basename "\$i").o%A -e .dmx\$(basename "\$i").e%A --partition $queue --ntasks=1 --cpus-per-task=1 --mem-per-cpu=4000 --parsable --oversubscribe <<__DEMUX__
#!/bin/bash
    source "$STAMPIPES/scripts/sentry/sentry-lib.bash"
    python3 $STAMPIPES/scripts/flowcells/demux_fastq.py   \
      \$suffix                     \
      --processing "$json"             \
      --outdir "$copy_from_dir"        \
      --mismatches "$dmx_mismatches"   \
      --lane "\$lane"                  \
      --ignore_failed_lanes            \
      "\$i"
__DEMUX__
   )
   DEMUX_JOBIDS="\$DEMUX_JOBIDS,\$jobid"
done
fi

if [[ -n \$DEMUX_JOBIDS ]]; then
   dmx_dependency=\$(echo \$DEMUX_JOBIDS | sed -e 's/,/,afterok:/g' | sed -e 's/^,afterok/--dependency=afterok/g')
fi

# copy files and prep collation/fastqc
copy_jobid=\$(sbatch --export=ALL -J "c-$flowcell" \$dmx_dependency -o "c-$flowcell.o%A" -e "c-$flowcell.e%A" --partition=$queue --cpus-per-task=1 --ntasks=1 --mem-per-cpu=1000 --parsable --oversubscribe <<'__COPY__'
#!/bin/bash
source "$STAMPIPES/scripts/sentry/sentry-lib.bash"
$link_command

# copy files
mkdir -p "$analysis_dir"
rsync -avP "$illumina_dir/InterOp" "$analysis_dir/"
rsync -avP "$illumina_dir/RunInfo.xml" "$analysis_dir/"
rsync -avP "$illumina_dir"/SampleSheet*.csv "$analysis_dir/"


# Copy each sample by itself, checking to see if we have a project_share_directory set
# This is very important to keep customer data separate from internal data.
(
    cd "$copy_from_dir"
    for dir in Project*/Sample* ; do
        [[ -d \$dir ]] || continue
        samp_number=\$(sed 's/.*DS\([0-9]*\).*/\1/' <<< "\$dir")
        [[ -n "\$samp_number" ]]
        destination=\$(jq -c -r ".libraries[] | select(.sample == \$samp_number) | .project_share_directory" ../processing.json)
        if [[ -z "\$destination" ]] ; then
            destination=$analysis_dir
        elif [[ ! -d "\$destination" ]] ; then
            echo "Destination \$destination does not exist! Please create it." >&2
            exit 1
        else
            destination=\$destination/fastq
        fi
        destination=\$destination/\$dir
        mkdir -p "\$destination"
        rsync -aL "\$dir/" "\$destination/"
    done
    for dir in Project*/LibraryPool* ; do
        [[ -d \$dir ]] || continue
        destination=$analysis_dir
        destination=\$destination/\$dir
        mkdir -p "\$destination"
        rsync -aL "\$dir/" "\$destination/"
    done
)


# create fastqc and collation scripts
cd "$analysis_dir"
# Remove existing scripts if they exist (to avoid appending)
rm -f fastqc.bash collate.bash run.bash

# Create fastqc scripts
python3 "$STAMPIPES/scripts/apilaneprocess.py" \
  --script_template "$STAMPIPES/processes/fastq/fastqc.bash" \
  --qsub-prefix .fq \
  --queue $queue \
  --sample-script-basename fastqc.bash \
  --flowcell_label "$flowcell" \
  --outfile fastqc.bash

# Create collation scripts
python3 "$STAMPIPES/scripts/apilaneprocess.py" \
  --script_template "$STAMPIPES/processes/fastq/collate_fastq.bash" \
  --qsub-prefix .collatefq \
  --queue $queue \
  --sample-script-basename "collate.bash" \
  --flowcell_label "$flowcell" \
  --outfile collate.bash

__COPY__
)

if [[ -n \$copy_jobid ]]; then
   copy_dependency=\$(echo \$copy_jobid | sed -e 's/^/--dependency=afterok:/g')
fi

# Collate
sbatch --export=ALL -J "collate-$flowcell" \$copy_dependency -o "collate-$flowcell.o%A" -e "collate-$flowcell.e%A" --partition=$queue --cpus-per-task=1 --ntasks=1 --mem-per-cpu=1000 --parsable --oversubscribe <<'__COLLATE__'
#!/bin/bash
source "$STAMPIPES/scripts/sentry/sentry-lib.bash"

cd "$analysis_dir"
# Remove existing scripts if they exist (to avoid appending)
rm -f fastqc.bash collate.bash run_alignments.bash run_aggregations.bash

# Create fastqc scripts
python3 "$STAMPIPES/scripts/apilaneprocess.py" \
  --script_template "$STAMPIPES/processes/fastq/fastqc.bash" \
  --qsub-prefix .fq \
  --queue $queue \
  --sample-script-basename fastqc.bash \
  --flowcell_label "$flowcell" \
  --outfile fastqc.bash

# Create collation scripts
python3 "$STAMPIPES/scripts/apilaneprocess.py" \
  --script_template "$STAMPIPES/processes/fastq/collate_fastq.bash" \
  --qsub-prefix .collatefq \
  --queue $queue \
  --sample-script-basename "collate.bash" \
  --flowcell_label "$flowcell" \
  --outfile collate.bash

bash collate.bash

# Wait for collation jobs to finish
while ( squeue -o "%j" | grep -q '^.collatefq*$flowcell') ; do
   sleep 60
done

# Run fastQC
bash fastqc.bash

# Set up of flowcell alignments
python3 "$STAMPIPES/scripts/alignprocess.py" \
  --flowcell "$flowcell"                     \
  --auto_aggregate                           \
  --qsub-queue queue0                        \
  --outfile run_alignments.bash

python3 "$STAMPIPES/scripts/poolprocess.py" --flowcell "$flowcell" --outfile run_pools.bash

# Set up of flowcell aggregations
curl -X POST "$LIMS_API_URL/flowcell_run/$flowcell_id/autoaggregate/" -H "Authorization: Token \$LIMS_API_TOKEN"

# Run alignments
bash run_alignments.bash
bash run_pools.bash

__COLLATE__

__BCL2FASTQ2__

fi

if [ -e "CopyComplete.txt" ] ; then
    echo -e "Setup complete. To kick everything off, type:\n\nbash run_bcl2fastq.sh"
else
    echo -e "Setup complete, sequencing still in progress. To queue everything up, type:\n\nnohup bash run_bcl2fastq.sh &"
fi
