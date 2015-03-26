#!/usr/bin/perl -w

#read a bunch of files with bits of info, make a big table

if (@ARGV < 2) {
    die "Usage:  $0 [base directory for find] [sample name]\n";
}
my ($searchDir,$name) = @ARGV;
unless (-d $searchDir) {
    die  "$0:  Failed to find directory ($searchDir)\n";
}

# For one reason or another, these fail rather than producing a real result
my @ignoreList = qw(  BAD_CYCLES CATEGORY IGNORED_READS PCT_ADAPTER PCT_CHIMERAS PCT_PF_READS_ALIGNED PCT_READS_ALIGNED_IN_PAIRS PF_ALIGNED_BASES PF_HQ_ALIGNED_BASES PF_HQ_ALIGNED_Q20_BASES PF_HQ_ALIGNED_READS PF_HQ_ERROR_RATE PF_HQ_MEDIAN_MISMATCHES PF_INDEL_RATE PF_MISMATCH_RATE PF_NOISE_READS PF_READS_ALIGNED READS_ALIGNED_IN_PAIRS STRAND_BALANCE );
my %ignore = ();
foreach my $key (@ignoreList) {
    $ignore{$key} = 1; 
}

# These are from Picard programs:
open INFILES, "find $searchDir -name 'xp.*.txt' |" or die;
while (<INFILES>) { 
    chomp:
    my $xpfile = $_;
    #warn "$xpfile\n";
    open IN, $xpfile or die "Failed to read $xpfile\n";
    while (<IN>) {
        chomp;
        my ($key,$value,@ignore) = split /\t/;
        unless (defined($ignore{$key})) {
            if (defined($value)) {
                if ($value =~ /^[\-\+0-9E\.,]+$/i) {
                    $value = decomma( $value );
                    print "$name\t$key\t$value\n";
                } else {
                    warn "Skipping non-numeric(?) key-value pair ($key,$value)\n";
                }
            }
        }
    }
    close IN;
}
close INFILES;

sub decomma {
    my $str = shift;
    if ($str =~ /Library/i) {
        # that's something else
        return $str;
    }
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

