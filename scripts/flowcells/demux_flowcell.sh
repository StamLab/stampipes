#!/bin/bash
set -e

source $MODULELOAD
module unload anaconda
module load python

usage(){
cat << EOF
usage: $0 -i <input_dir> -o <output_dir> -p <processing.json>

Setup analysis for a flowcell

OPTIONS:
-h        Show this message
-i [dir]  Unaligned input directory
-o [dir]  Unaligned output directory
-p [json] processing.json
-l [int]  lane; integer (default: all lanes)
-n        Dry-run - do not actually demux
EOF
}

check_folder(){
  if [ ! -d "$1" ] ; then
    echo "Directory $1 does not exist, aborting." >&2
    exit 1
  fi
}

demux_script="$STAMPIPES/scripts/flowcells/demux_fastq.py"
indir=
outdir=
processing=
LANE=
dryrun=
mismatches=0
while getopts ":hi:o:p:m:l:n" opt ; do
  case $opt in
    h)
      usage
      exit 0
      ;;
    i)
      indir=$(readlink -f "$OPTARG")
      ;;
    o)
      outdir="$OPTARG"
      ;;
    p)
      processing=$(readlink -f "$OPTARG")
      ;;
    m)
      mismatches="$OPTARG"
      ;;
    l)
      LANE="L00$OPTARG"
      ;;
    n)
      dryrun="--dry-run"
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

if [ -z "$processing" ] || [ -z "$indir" ] || [ -z "$outdir" ] ; then
  usage >&2
  exit 1
fi

if [ ! -s "$processing" ] ; then
  echo "$processing is not a non-empty file" >&2
  exit 1
fi

inputfiles=($(find "$indir" -name "*Undetermined_*$LANE*fastq.gz"))

for i in "${inputfiles[@]}" ; do
  lane=$( sed 's/.*_L\(00.\)_.*/\1/' <(basename "$i" ))

  # If dryrun; we want the output on stdout, not qsubbed.
  if [ -z "$dryrun" ] ; then
    submitjob="qsub -cwd -V -q all.q -N .dmx$(basename $i)"
  else
    submitjob="bash"
  fi

  cat <<__DEMUX__ | $submitjob
    python "$demux_script"        \
      --autosuffix                \
      --processing "$processing"  \
      --outdir "$outdir"          \
      --mismatches "$mismatches"  \
      --lane "$lane"              \
      $dryrun                     \
      $i
__DEMUX__

done
