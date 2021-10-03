#!/bin/bash

for dir in "$@" ; do
  [[ -d "$dir/work" ]] || continue
  dir=$(readlink -f "$dir")
  find "$dir/work" -type f \
    | grep -v -F -f <( \
      find "$dir"/output* -type l -print0 \
      | xargs --no-run-if-empty -0 readlink -f \
    ) \
    | xargs --no-run-if-empty du -scb | tail -n1  | sed "s!total!$dir!"
done
