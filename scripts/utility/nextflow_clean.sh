#!/bin/bash


# 1   visit directory containing work/ dir
# 2   find all symlinks within top-level dir (and subdirs), excluding those in work dir.
# 3   for each symlink found:
# 3A    get the target of the symlink, resolving any symlink chains (e.g: the symlink may point to another symlink, etc)
# 3B    if the symlink does not point into a work/ directory, skip it.
# 3C    make a copy of the link for backup purposes in case we get interrupted (e.g: mysymlink.bak)
# 3D    replace the symlink with a hard link to the real target
# 3E    remove the symlink backup
# 4   repeat step 2, exiting with error if any symlinks remain.
# 5   remove work directory
# 6   set LIMS status


set -euo pipefail

AGG_DIR_REGEX=/net/seq/data/aggregations/LN[0-9]+/aggregation-[0-9]+

#########################

function usage(){
  echo "Usage: $0 [-f] dir [dirs...]"  >&2
  exit 2
}

# Parse arguments
FORCE=
while getopts ":fh" arg; do
  case $arg in
    h)
      usage
      ;;
    f)
      FORCE=TRUE
      ;;
    *)
      usage
      ;;

  esac
done
shift $((OPTIND -1))
[[ -z "$*" ]] && usage

# Print arguments to stderr (printf style format string) and exit
DIE_F(){
  printf "ERROR: "
  # shellcheck disable=SC2059
  printf "$@"
  printf "\n"
  exit 1
}
INFO_F(){
  printf "INFO: "
  # shellcheck disable=SC2059
  printf "$@"
  printf "\n"
}

RUN_IF_FORCE(){
  if [[ "$FORCE" == TRUE ]] ; then
    printf "INFO: Running:' "
    printf "%q " "$@"
    printf "'\n"
    "$@"
  else
    printf "INFO: Dry-run, would have run: '"
    printf "%q " "$@"
    printf "'\n"
  fi
}

lims_patch_by_url() {
  /usr/bin/curl \
  --request PATCH \
  --data "$2" \
  "$1" \
  -H "Authorization: Token $LIMS_API_TOKEN" \
  -k \
  2>/dev/null
}

# $1: url, $2: field=val
lims_patch () {
  lims_patch_by_url "$LIMS_API_URL/$1" "$2"
}

lock_aggregation() {
  local aggdir
  local aggid
  aggdir=$(readlink -f $1)
  if  [[ ! "$aggdir" =~ $AGG_DIR_REGEX ]] ; then
    INFO_F "Not an aggregation directory, not locking"
    return
  fi

  # Get aggregation ID from directory name
  aggid=${aggdir/*aggregation-/}
  if [[ ! $aggid  =~ ^[0-9]+$ ]] ; then
    DIE_F "Couldn't parse aggregation directory '%s' ('%s') into aggid (got '%s')" "$1" "$aggdir" "$aggid"
  fi

  lims_patch "aggregation/$aggid/" 'locked=True'
}

# Find output symlinks in a directory
# Expects them to be in subdirectories starting with 'output'
find_output_symlinks(){
  # find all symlinks within top-level dir (and subdirs), excluding those in work dir.
  find "$1" -maxdepth 1 -name 'output*' -type d -print0 \
    | xargs --no-run-if-empty -0 -I '{}' find '{}' -type l
}

has_broken_symlink(){
  local broken_symlink
  broken_symlink=$(find "$1" -maxdepth 1 -name 'output*' -type d -print0 \
    | xargs --no-run-if-empty -0 -I '{}' find '{}' -type l -xtype l -print -quit)
  [[ -n "$broken_symlink" ]]
}

# Following symlinks, see if the two arguments point to the same region of the disk.
is_same_inode(){
  [[ "$1" -ef "$2" ]]
}

process_dir() {
  local dir=$1
  # Start subshell to isolate our changes/work
  # Also redirecting logs
  (
    INFO_F "Processing dir '%s'" "$dir"
    # Make sure dir exists
    [[ -d "$dir" ]] || DIE_F "Directory does not exist: '%s'" "$dir"
    # It should also contain a workdir
    [[ -d "$dir/work" ]] || DIE_F "Work directory does not exist: '%s/work'" "$dir"

    set -euo pipefail  # This is inherited, but in case somebody changes the top declaration.

    # 1) visit directory containing work/ dir
    cd "$dir"

    # 1.5) Check for broken symlinks
    has_broken_symlink . &&
      DIE_F "Broken symlink detected. dir=%s" "$dir"

    # 2) find symlinks
    local symlinks=()
    readarray -t symlinks < <( find_output_symlinks . )

    # Check 
    # https://stackoverflow.com/questions/7577052/bash-empty-array-expansion-with-set-u
    # shellcheck disable=SC2199
    if [[ -z ${symlinks[@]+"${symlinks[@]}"} ]] ; then
      INFO_F "No input symlinks found. dir=%s" "$dir"
      symlinks=()
    else

      # 3)
      for symlink in "${symlinks[@]}" ; do
        # 3A) Get symlink target
        local target
        target=$(readlink -f "$symlink")
        # 3B) Die if not in work directory
        [[ "$target" =~ /work/ ]] ||
          DIE_F "Target not in workdir. dir='%s';symlink='%s';target='%s'" "$dir" "$symlink" "$target"

        # 3?) Make sure symlink isn't broken
        [[ -e "$target" ]] ||
          DIE_F "Target not found. dir='%s';symlink='%s';target='%s'" "$dir" "$symlink" "$target"

        is_same_inode "$target" "$symlink" ||
          DIE_F "Somehow, target & symlink inodes differ. dir='%s';symlink='%s';target='%s'" "$dir" "$symlink" "$target"

        # 3C) make a copy of the link for backup purposes in case we get interrupted (e.g: mysymlink.bak)
        local symlink_bak=$symlink.bak
        RUN_IF_FORCE cp --no-dereference "$symlink" "$symlink_bak"

        # 3D) Create hard link to real target, replacing soft link
        local hardlink=$symlink
        RUN_IF_FORCE ln -f --logical "$target" "$hardlink"

        # Triple-check, make sure the hard link and backup symlink point to the same place
        if [[ "$FORCE" == TRUE ]] ; then
          is_same_inode "$symlink_bak" "$hardlink" ||
            DIE_F "Somehow, the symlink and hard link point to different files. hardlink='%s', target='%s'" "$hardlink" "$target"
        fi

        # 3E) remove backup symlink
        RUN_IF_FORCE rm "$symlink_bak"
      done
    fi

    # Only run check for symlinks if we deleted
    if [[ $FORCE == TRUE ]] ; then
      # 4)
      local remaining_symlinks
      remaining_symlinks=$(find_output_symlinks .)
      [[ -z "$remaining_symlinks" ]] ||
        DIE_F "Some symlinks remain, not removing workdir. dir='%s';remaining='%q'" "$dir" "${remaining_symlinks[*]}"
    else
      INFO_F "Dry-run, skipping the check for remaining symlinks"
    fi

    RUN_IF_FORCE rm -r --one-file-system --preserve-root "work"

    RUN_IF_FORCE lock_aggregation .

    INFO_F "Done. dir='%s'" "$dir"

  ) | tee "$dir"/freezing.log
}

for directory in "$@" ; do
  if [[ -d "$directory/work" ]]  ; then
    process_dir "$(readlink -f $directory)"
  else
    INFO_F "Dir '%s' does not contain a workdir" "$directory"
  fi
done
