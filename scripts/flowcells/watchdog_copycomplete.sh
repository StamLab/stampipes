#!/bin/bash
#SBATCH --job-name=watchdog-copycomplete
#SBATCH --partition=hpcz-test,hpcy-test
#SBATCH --time=00:55:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G

# Watchdog script to ensure flowcells are set up and bcl-convert is launched
# when CopyComplete.txt arrives.
# Intended to be run hourly via scrontab; can also be submitted directly with sbatch.

# Default configuration
SEQUENCER_MOUNT=${SEQUENCER_MOUNT:-/net/seq/data2/sequencers}
DAYS_BACK=${DAYS_BACK:-10}

set -e -u -o pipefail

if [ -z "$SEQUENCER_MOUNT" ] ; then
  echo "Sequencer mount not set; set \$SEQUENCER_MOUNT or supply an argument" >&2
  exit 1
fi

if [ ! -d "$SEQUENCER_MOUNT" ] ; then
  echo "Sequencer mount '$SEQUENCER_MOUNT' does not exist" >&2
  exit 1
fi

# Takes an integer and gets the datecode for that many days ago
# Returns both 4-digit and 2-digit year formats
fmt_days_ago() {
  local day=$1
  local y4
  local y2
  y4=$(date --date "$day days ago" +"%Y%m%d")
  y2=$(date --date "$day days ago" +"%y%m%d")
  echo "$y4 $y2"
}

# Takes flowcell_root and days_back
get_recent_flowcells() {
  local flowcell_root=$1
  local days_back=$2
  for days_ago in $(seq "$days_back" -1 0); do
    # Get both 4-digit and 2-digit year formats
    read -r y4 y2 <<< $(fmt_days_ago "$days_ago")
    ls -d "$flowcell_root/${y4}_"* "$flowcell_root/${y2}_"* 2>/dev/null || true
  done
}

# Takes a directory as argument
flowcell_is_setup() {
  [ -s "$1/run_bcl2fastq.sh" ]
}

ts() { date +"%Y-%m-%d %H:%M:%S"; }

log() {
  local fc=$1
  shift
  echo "$(ts)  $(basename "$fc")  $*"
}

logmsg() {
  echo "$(ts)  $*"
}

# Takes a directory as argument
setup_flowcell() {
  local fc_dir="$1"
  # Run in a subshell to not affect current path
  (
    set -e
    cd "$fc_dir"
    logmsg "Running setup.sh in $fc_dir..."
    bash "$STAMPIPES/scripts/flowcells/setup.sh" || exit 1
    logmsg "Launching run_bcl2fastq.sh in $fc_dir..."
    bash run_bcl2fastq.sh
  )
}

logmsg "Starting watchdog check"
logmsg "Scanning $SEQUENCER_MOUNT for flowcells from the last $DAYS_BACK days"

for fc in $(get_recent_flowcells "$SEQUENCER_MOUNT" "$DAYS_BACK"); do
  [ -d "$fc" ] || continue
  # Only act if CopyComplete.txt exists
  if [ ! -e "$fc/CopyComplete.txt" ] ; then
    continue
  fi

  if flowcell_is_setup "$fc"; then
    log "$fc" "Already set up"
    continue
  fi

  log "$fc" "CopyComplete.txt found but not set up. Setting up..."
  if setup_flowcell "$fc"; then
    log "$fc" "Setup successful"
  else
    log "$fc" "Setup failed"
  fi
done

logmsg "Watchdog check complete"
