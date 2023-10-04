#!/bin/bash
# This script finds and gzips all StarSOLO matrix/barcodes/features files.
root_directory=${1:-.}
threads=${2:-10}

#shellcheck disable=SC2037
name_query=( '(' -name '*.mtx' -o -name barcodes.tsv -o -name features.tsv ')' )

# Gzip regular files
find "$root_directory" -type f \
  "${name_query[@]}" \
  -print0 \
  | xargs --no-run-if-empty -0 -n 1 -P "$threads" gzip

# Gzip any targets pointed to by a symlink
find "$root_directory" -type l \
  "${name_query[@]}" \
  -print0 \
  | xargs -0 -n1 readlink -f \
  | sort -u \
  | xargs --no-run-if-empty -n 1 -P "$threads" gzip

# Create new symlinks that point to the new target
find "$root_directory" -type l \
  "${name_query[@]}" \
  | while read -r symlink ; do
    target=$(readlink "$symlink")
    realtarget=$(readlink -f "$symlink")
    if [[ -e "$realtarget.gz" ]] ; then
      ln -s "$target.gz" "$symlink.gz"
      rm "$symlink"
    fi
done
