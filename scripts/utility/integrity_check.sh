#!/bin/bash
# Performs quick-checking of file integrity, based on file extension
# Exit status is the number of failing files (0 if all are good)

if [[ "$#" -eq 0 ]]; then
  echo "Usage: $0 file [files...]" >&2
  exit 0
fi

exitstatus=0

# Checks the exit status of the last command
report_status() {
  if [[ $? -ne 0 ]]; then
    echo "$arg: FAIL - invalid file" >&2
    exitstatus=$((exitstatus + 1))
  fi
}

# Loop through each argument and check file integrity
for arg in "$@"; do
  # Make sure the file exists
  if [[ ! -f "$arg" ]]; then
    echo "$arg: FAIL - is not a file" >&2
    continue
  fi

  ext=${arg##*.}

  case "$ext" in
    "bam")
      samtools quickcheck "$arg"
      ;;
    "starch")
      [[ "$(unstarch --is-starch "$arg")" == "1" ]]
      ;;
    *)
      echo "$arg: WARN - Checking filetype $ext is not supported" >&2
      ;;
  esac

  report_status
done

exit $exitstatus
