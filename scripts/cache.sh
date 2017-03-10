#!/bin/bash
#
# This script will cache a file or directory (recursively) and output the name
# of the cached version.
#
# Destination can be controlled by setting the $CACHE_DIR envvar
#
# If the cache cannot be created, the original filename will be returned, and
# exit status will still be zero

usage="Usage: $0 input_file"

CACHE_DIR="${CACHE_DIR:-/node-data/cache/}"

error(){
  echo "ERROR:" "$1" 1>&2
  echo "$src"
  exit 0
}

# Check arguments
if [ $# -lt 1 ] ; then
  echo $usage
  exit 1
fi

# Check that src argument exists
src=$(readlink -f "$1")
if [ ! -e "$src" ] ; then
  error "Target file does not exist: $1"
fi

# Check that cache dir exists
mkdir -p "$CACHE_DIR"
if [ ! -d "$CACHE_DIR" ] ; then
  error "Cache directory could not be created: $CACHE_DIR"
fi

# Name destination & lockfile
dest="$(readlink -m $CACHE_DIR/$src)"
mkdir -p "$(dirname $dest)"
lockfile="$CACHE_DIR.lock"

# Add trailing slash for directories
if [ -d "$src" ] ; then
  src="$src/"
fi

# Lock & rsync
(
flock -n "$lockfile" rsync --quiet --recursive --archive "$src" "$dest"
) || error "Could not lock/transfer"

# Return destination
echo "$dest"
