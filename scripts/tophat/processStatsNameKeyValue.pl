#!/usr/bin/perl -w

my %unames = ();
my %db = ();

# Read stuff produced by other scripts
while (<STDIN>) {
    chomp;
    my ($name,$key,$value) = split /\t/;
    my $prevValue = $db{"$name-$key"}; 
    unless (defined($prevValue) and ($prevValue > 0)) { 
        $db{"$name-$key"} = $value;
    }
    $unames{$name} = 1;
}

#CORRECT_STRAND_READS INCORRECT_STRAND_READS PCT_CORRECT_STRAND_READS

#my $header = "lane\tinput_reads\tmapped\t\%rRNA\t\%duplicates\tlib_size_est\texon:intron\t\%intergenic\t\%chrM\n";
my $header = "sample\tinput_reads\tmapped\t\%rRNA\t\%duplicates\texon:intron\t\%intergenic\t\%chrM\t\%correct_strand\n";
print $header;
foreach my $name (sort keys %unames) {
    # where stuff aligned
    my $exonIntronRatio = "TBD";
    my $percentIntergenic = "TBD";
    my $alignedBases = $db{"$name-PF_ALIGNED_BASES"};
    unless (defined($alignedBases)) {
        # different for single-end?
        $alignedBases = $db{"$name-PF_BASES"};
    }
    if (defined($alignedBases)) {
        # probably have the rest of these
        my $coding = decomma( $db{"$name-CODING_BASES"} );
        my $UTR = decomma( $db{"$name-UTR_BASES"} );
        my $intronic = decomma( $db{"$name-INTRONIC_BASES"} );
        my $intergenic = decomma( $db{"$name-INTERGENIC_BASES"} );
        $exonIntronRatio = round2( ($coding + $UTR) / (1 + $intronic) );
        $percentIntergenic = $intergenic / decomma( $alignedBases );
    }
    my $inputReads = get( $name, "total_reads" );
    my $mappedReads = get( $name, "mapped_reads" );
    # duplicates
    my $pctDup = "TBD";
    my $libSizeEst;
    if (defined($db{"$name-PERCENT_DUPLICATION"} ) ) {
        $pctDup = pct( get( $name, "PERCENT_DUPLICATION" ) );
        $libSizeEst = get( $name, "ESTIMATED_LIBRARY_SIZE" );
    } else {
        my $duplicates = get( $name, "duplicates" );
        if (defined($mappedReads) and ($mappedReads>0)) {
            $pctDup = pct( $duplicates / $mappedReads );
        }
    }
    #
    my $rRNA = "TBD";
    my $rRNAcount = get($name, "rRNA" );
    if (defined($rRNAcount) and defined($mappedReads) and ($mappedReads>0)) {
        $rRNA = pct( $rRNAcount / ($rRNAcount + $mappedReads ) );
    }
    #
    my $chrM = "TBD";
    if (defined($mappedReads) and ($mappedReads>0)) {
        $chrM = $db{"$name-chrM"};
        unless (defined($chrM)) {
            $chrM = pct( get( $name, "mapped_chrM" ) / $mappedReads );
        }
    }
    my $pctCorrectStrand = get( $name, "PCT_CORRECT_STRAND_READS" );
    if (defined( $pctCorrectStrand )) {
        $pctCorrectStrand = pct( $pctCorrectStrand );
    } else {
        $pctCorrectStrand = "N/A"; # ?
    }
    #
    #my @fields = ($inputReads,$mappedReads,$rRNA, $pctDup, $libSizeEst,$exonIntronRatio, pct($percentIntergenic),$chrM );
    my @fields = ($inputReads,$mappedReads,$rRNA, $pctDup, $exonIntronRatio, pct($percentIntergenic),$chrM, $pctCorrectStrand );
    my $line = "$name";
    foreach my $field (@fields) {
        $line .= "\t" . (defined($field) ? $field : "TBD" );
    }
    print "$line\n";
}

# The value of this function is that it complains
sub get {
    my ($name,$key) = @_;
    my $result = $db{"$name-$key"};
    unless (defined($result)) {
        warn "No $key value for $name\n";
    }
    return $result;
}

    
sub round2 {
    my $value = shift;
    return int( 0.5 + (100 * $value) ) / 100;
}

sub pct {
    my $value = shift;
    unless (defined($value) and $value ne "TBD") {
        return "undef";
    }
    my $pct = 100 * $value;
    my $rounded = int( 0.5 + (100 * $pct) ) / 100;
    return $rounded;
}

sub decomma {
    my $str = shift;
    return 0 + join( "", split( /,/, $str ) );
}

#  ::::::::::::::
#  DS18745-FCC04N5-3.dir/xp.picard.DS18745-FCC04N5-3.RnaSeq.txt
#  ::::::::::::::
#  PF_ALIGNED_BASES	30,489,721,211
#  CODING_BASES	10,385,438,735
#  UTR_BASES	7,252,288,405
#  INTRONIC_BASES	7,048,608,265
#  INTERGENIC_BASES	5,803,385,806
#  PCT_CODING_BASES	0.340621
#  PCT_UTR_BASES	0.23786
#  PCT_INTRONIC_BASES	0.23118
#  PCT_INTERGENIC_BASES	0.190339
#  PCT_MRNA_BASES	0.578481
#  ::::::::::::::
#  DS18745-FCC04N5-3.dir/xp.picard.DS18745-FCC04N5-3.InsertSize.txt
#  ::::::::::::::
#  MEDIAN_INSERT_SIZE	234	insert_size	73	74	75	76
#  MEDIAN_ABSOLUTE_DEVIATION	76	fr_count	5	9	49	5,111
#  MIN_INSERT_SIZE	73
#  MAX_INSERT_SIZE	248,574,133
#  MEAN_INSERT_SIZE	248.078493
#  STANDARD_DEVIATION	145.455678
#  READ_PAIRS	181,634,767
#  PAIR_ORIENTATION	FR
#  ::::::::::::::
#  DS18745-FCC04N5-3.dir/xp.picard.DS18745-FCC04N5-3.AlignmentSummary.txt
#  ::::::::::::::
#  CATEGORY	FIRST_OF_PAIR	SECOND_OF_PAIR	PAIR
#  TOTAL_READS	201,871,579	199,346,589	401,218,168
#  PF_READS	201,871,579	199,346,589	401,218,168
#  MEAN_READ_LENGTH	76	76	76
#  
#  duplicates: 
#  UNPAIRED_READS_EXAMINED	22,308,423	VALUE	1
#  READ_PAIRS_EXAMINED	179,259,621
#  UNPAIRED_READ_DUPLICATES	17,767,993
#  READ_PAIR_DUPLICATES	53,915,693
#  PERCENT_DUPLICATION	0.329806
#  ESTIMATED_LIBRARY_SIZE	234,649,023
#  
#  DS18745-FCC04N5-3	reads	521,172,344
#  DS18745-FCC04N5-3	rRNA	1,703,488
#  DS18745-FCC04N5-3	mapped	428,920,747
#  
#  
#  IV. RNA-seq Sequence Experiment QC Metrics.
#  The following QC metrics should be determined and monitored to ensure that REMC
#  RNA-seq libraries are of high quality:
#  
#  1. Exon:Intron ratio;
#  -assessed to detect potential genomic contamination (the mRNA-seq protocol
#  should include a standard DNAse treatment step).
#  AND
#  6. Fraction of reads mapping to intergenic regions; and
#  -assessed to detect potential genomic contamination (the RNA-seq protocol
#  should include a standard DNAse treatment step).
#  
#      EDH:  Picard file:  DS18745-FCC04N5-3.dir/xp.picard.DS18745-FCC04N5-3.RnaSeq.txt
#  
#  
#  2. Overall transcript coverage;
#  -assessed to detect mRNA degradation and degree to which full-length mRNA
#  was randomly sampled. Uniform coverage should be obtained for 1Kb transcripts
#  
#      TODO  (how is this represented?)  use the pdf?
#  
#  
#  3. Total number of duplicate reads (identical forward and reverse read starts);
#  -a measure of library diversity. Tissue and cell type dependent but outliers should
#  be investigated as potential library failures.
#      
#      EDH:  Picard MarkDuplicates 
#  
#  4. Fraction of reads mapping to mitochondrial transcripts;
#  -assessed to detect RNA degradation. Tissue and cell type dependent but outliers
#  should be investigated as potential library failures.
#      
#      EDH:  Counting...  *.dir/chromosome_counts.txt
#  
#  
#  5. Fraction of reads mapping to ribosomal transcripts;
#  -assessed to detect degree of polyA enrichment. Tissue and cell type dependent
#  but outliers should be investigated as potential library failures.
#  
#      EDH:  Separate Bowtie alignment.  eventually (NEXT TIME) by segmented reads...
#      
#  
#  7. Percentage of reads inappropriately aligned to anti-sense strands should be <=1%.
#  This can be determined by enumerating the fraction of exon-exon junction reads aligned
#  to the anti-sense strand. Alternatively known poly-adenylated foreign RNA standards
#  can be spiked into the RNA prior to library construction and used to determine the level
#  of artifact anti-sense transcription.
#  
#      EDH:  Is the protocol stranded?  TODO:  Browser view...
#  
#  
#  

