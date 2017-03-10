###################
# This makefile takes one set of FastQ files and aligns them using BWA, using a trimmer.
# It then provides a sorted BAM file.  The TMPDIR is best passed in as an empty directory
# to prevent collisions.
###################
# These variables must be passed in or set for the makefile to work.  If the genome's
# FAI file is not at $(BWAINDEX).fai, then it must also be specified under FAI.
###################
# BWAINDEX=/path/to/genome
# FASTQ_FILE=/path/to/R.fq.gz
# OUTBAM=/path/to/out.sorted.bam
###################
# REQUIRED MODULES
###################
# module load bwa/0.7.12
# module load samtools
# module load python
###################

SHELL = bash

# Use the default PATH version of each of these programs unless overridden
SAMTOOLS ?= samtools
BWA ?= bwa

FAI ?= $(BWAINDEX).fai

READ_LENGTH ?= 36

# Ideally this will be set in the environment or on the command line to an actual temp directory
TMPDIR ?= $(shell pwd)

# Discard reads failing chastity filter, use 32 bases as a seed
# and set the mismatch level appropriate to the read length
# NOTE: THESE ARE NOT WORKING FOR CHANGING IN QMAKE SCRIPT
BWA_ALN_OPTIONS ?= -Y -l 32 -n 0.04

# Include up to 10 matches for each read and
# set the inset size to a generous 750
BWA_SAMSE_OPTIONS ?= -n 10

# List the intermediate files that can be deleted when we finish up
.INTERMEDIATE : 

all : info align 

info : 
	@echo "------"
	@echo "PROGRAM VERSIONS"
	@echo "------"
	@echo "bwa:"
	@$(BWA) 2>&1 | grep "Version"
	@echo "samtools:"
	@$(SAMTOOLS) 2>&1 | grep "Version"
	@echo "------"
	@echo "METADATA"
	@echo "------"
	@echo "SAMPLE_NAME: " $(SAMPLE_NAME)
	@echo "BWAINDEX: " $(BWAINDEX)
	@echo "FAI: " $(FAI)
	@echo "READ LENGTH: " $(READ_LENGTH)
	@echo "------"
	@echo "PROGRAM OPTIONS"
	@echo "------"
	@echo "BWA_ALN_OPTIONS: " $(BWA_ALN_OPTIONS)
	@echo "BWA_SAMSE_OPTIONS: " $(BWA_SAMSE_OPTIONS)
	@echo "------"

align : $(OUTBAM)

# Index the resulting BAM
$(OUTBAM).bai : $(OUTBAM)
	$(SAMTOOLS) index $^

# Copy the final sorted bam to its finished place
$(OUTBAM) : $(TMPDIR)/align.sorted.bam
	rsync $^ $@

# Create sorted BAM files
$(TMPDIR)/align.sorted.bam : $(TMPDIR)/align.filtered.bam
	time $(SAMTOOLS) sort $^ > $@ && echo made $@ >&2

# Create unsorted filtered BAM files
$(TMPDIR)/align.filtered.bam : $(TMPDIR)/align.unsorted.bam
	time $(PYTHON) $(STAMPIPES)/scripts/bwa/filter_reads.py $^ $@

# Create unsorted raw BAM files
$(TMPDIR)/align.unsorted.bam : $(TMPDIR)/align.sam
	time $(SAMTOOLS) view -bT $(FAI) $^ > $@ && echo made $(TMPDIR)/align.unsorted.bam >&2

# Create the SAM files from each pair of SAI and FASTQ files
$(TMPDIR)/align.sam : $(TMPDIR)/align.sai $(FASTQ_FILE)
	time $(BWA) samse $(BWA_SAMSE_OPTIONS) $(BWAINDEX) $^ > $@ && echo made $(TMPDIR)/align.sam >&2

# Make SAI alignments for each FASTQ file
$(TMPDIR)/%.sai : $(FASTQ_FILE)
	time $(BWA) aln $(BWA_ALN_OPTIONS) $(BWAINDEX) $(FASTQ_FILE) > $(TMPDIR)/$*.sai && echo made $(TMPDIR)/$*.sai >&2
