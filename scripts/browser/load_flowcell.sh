#!/bin/bash

set -e -o pipefail

flowcell_name=${1:(-5)}
priority=$2

flowcell_dir=$(ls $FLOWCELLS/FC${flowcell_name}_*tag -d)

die(){
  echo $@ >&2
  exit 1
}

if [[ ! -s $HOME/.hg.conf ]] ; then
  die "Must be on the same server as the browser"
fi

if [[ -z "$priority" ]] ; then
  die "Usage: $0 flowcell_label priority"
fi

if [[ -z "$flowcell_dir" ]] ; then
  die "Couldn't find flowcell dir for $flowcell"
fi

source "$PYTHON3_ACTIVATE"

cd "$flowcell_dir"

python "$STAMPIPES/scripts/lims/get_processing.py" -f "$flowcell_name" --quiet

python "$STAMPIPES/scripts/browser/make_trackhubs.py" -p "$priority" --quiet