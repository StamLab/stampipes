###################
# These variables must be passed in or set for the makefile to work.  If the genome's
# FAI file is not at $(BWAINDEX).fai, then it must also be specified under FAI.
###################
# SAMPLE_NAME := Sample_NoIndex_L007
# BWAINDEX := /path/to/genome
###################
# REQUIRED MODULES
###################
# module load samtools
# module load python
###################

SHELL = bash

# Use the default PATH version of each of these programs unless overridden
SAMTOOLS ?= samtools

FAI ?= $(BWAINDEX).fai

READ_LENGTH ?= 36

MIN_MAP_QUALITY ?= 10
MAX_MISMATCHES ?= 2

EXCLUDE_FLAG ?= 512

# Filter to a minimum mapping quality, where both pairs
# are properly paired in the mapping
SAMTOOL_OPTIONS ?= -F $(EXCLUDE_FLAG)

# Ideally this will be set to something else set in the environment or on the command line
TMPDIR ?= $(shell pwd)

# where our results files go
OUTDIR ?= $(shell pwd)

INBAM ?= $(OUTDIR)/$(SAMPLE_NAME).sorted.bam
OUTBAM ?= $(OUTDIR)/$(SAMPLE_NAME).uniques.sorted.bam

all : info metrics uniques $(INSERTMETRICS)

info : 
	@echo "------"
	@echo "PROGRAM VERSIONS"
	@echo "------"
	@echo "samtools:"
	@$(SAMTOOLS) 2>&1 | grep "Version"
	@echo "------"
	@echo "METADATA"
	@echo "------"
	@echo "SAMPLE_NAME: " $(SAMPLE_NAME)
	@echo "BWAINDEX: " $(BWAINDEX)
	@echo "FAI: " $(FAI)
	@echo "READ LENGTH: " $(READ_LENGTH)

metrics : 

uniques : $(INBAM).bai $(OUTBAM).bai

# Index uniquely mapping reads 
$(OUTBAM).bai : $(OUTBAM)
	time $(SAMTOOLS) index $^
                
# Sorted uniquely mapping reads BAM
$(OUTBAM) : $(INBAM)
	time $(SAMTOOLS) view $(SAMTOOL_OPTIONS) -b $(INBAM) > $@

# Index sorted BAM file
$(INBAM).bai : $(INBAM)
	time $(SAMTOOLS) index $^	
