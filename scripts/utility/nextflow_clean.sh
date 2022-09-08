#!/bin/bash

# It seems clear that we're using a lot more space in Nextflow work/ dirs than I understood, and that the benefits of keeping it around are outweighed by the costs.
# I'll have to put together a script today to safely clean out Nextflow work directories for the pipeline.
# brainstorming pseudo-steps:

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

#      very cool.  I think we add in an accounting step which will set the LIMS representation to "locked" for aggregations (and alignments?) were we do this?
#      The term I'm thinking that describes the nextflow work dir/symlink/file move effort is "fossilize", cause I like a good metaphor.


set -euo pipefail

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

# Find output symlinks in a directory
# Expects them to be in subdirectories starting with 'output'
find_output_symlinks(){
  # find all symlinks within top-level dir (and subdirs), excluding those in work dir.
  find "$1" -maxdepth 1 -name 'output*' -type d -print0 \
    | xargs -0 -I '{}' find '{}' -type l
}

# Gets the device ID and inode of a file, following symlinks
# If two files have the same return value, they point to the same bytes-on-disk.
get_inode(){
  stat -L -c '%d:%i' "$1"
}
is_same_inode(){
  [[ $(get_inode "$1") == $(get_inode "$2") ]]
}

process_dir() {
  local dir=$1
  INFO_F "Processing dir '%s'" "$dir"
  # Make sure dir exists
  [[ -d "$dir" ]] || DIE_F "Directory does not exist: '%s'" "$dir"
  # It should also contain a workdir
  [[ -d "$dir/work" ]] || DIE_F "Work directory does not exist: '%s/work'" "$dir"

  # Start subshell to isolate our changes/work
  (
    set -euo pipefail  # This is inherited, but in case somebody changes the top declaration.

    # 1) visit directory containing work/ dir
    cd "$dir"

    # 2) find symlinks
    local symlinks=()
    readarray -t symlinks < <( find_output_symlinks . )

    # Check 
    # https://stackoverflow.com/questions/7577052/bash-empty-array-expansion-with-set-u
    # shellcheck disable=SC2199
    [[ -z ${symlinks[@]+"${symlinks[@]}"} ]] &&
      INFO_F "No input symlinks found. dir=%s" "$dir"

    # 3)
    for symlink in ${symlinks[@]+"${symlinks[@]}"} ; do
      # 3A) Get symlink target
      local target
      target=$(readlink -f "$symlink" || true)
      [[ -n "$target" ]] ||
        DIE_F "Target empty, likely broken symlink. dir='%s';symlink='%s';target='%s'" "$dir" "$symlink" "$target"

      # 3B) Die if not in work directory
      [[ "$target" =~ /work/ ]] ||
        DIE_F "Target not in workdir. dir='%s';symlink='%s';target='%s'" "$dir" "$symlink" "$target"

      # TODO: redundant with -n check?
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
      #if [[ -d "$target" ]] ; then
      #  # TODO: This is probably slow, is there a faster way?
      #  # e.g: remove link, create dir, then hardlink the contents?
      #  RUN_IF_FORCE rm "$hardlink"
      #  RUN_IF_FORCE cp -r "$target" "$hardlink"
      #else
        RUN_IF_FORCE ln -f --logical "$target" "$hardlink"
      #fi

      # Triple-check, make sure the hard link and backup symlink point to the same place
      if [[ "$FORCE" == TRUE ]] ; then
        is_same_inode "$symlink_bak" "$hardlink" ||
          DIE_F "Somehow, the symlink and hard link point to different files. hardlink='%s', target='%s'" "$hardlink" "$target"
      fi

      # 3E) remove backup symlink
      RUN_IF_FORCE rm "$symlink_bak"
    done

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

    INFO_F "Done. dir='%s'" "$dir"

  ) | tee "$dir"/freezing.log
}

for directory in "$@" ; do
  process_dir "$directory"
done
