BEGIN {
  readsep = "~";
  unmatched_reads = 0;
}
{
  if( NR == 1 ) {
    previous = $0;
    seekname = $1;
    next;
  }

  if( seekname ) {
    if ( seekname != $1 ) {
      # print "UNMATCHED: " seekname > "/dev/stderr";
      previous = $0;
      seekname = $1;
      unmatched_reads = unmatched_reads + 1;
      next;
    }

    print $0 readsep previous;
    seekname = "";
    next;
  }
  
  previous = $0;
  seekname = $1;
}
END {
  if( unmatched_reads ) {
    print "UNMATCHED READS: " unmatched_reads > "/dev/stderr";
  }
}
