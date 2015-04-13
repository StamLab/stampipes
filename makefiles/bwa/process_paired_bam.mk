###################
# These variables must be passed in or set for the makefile to work.  If the genome's
# FAI file is not at $(BWAINDEX).fai, then it must also be specified under FAI.
###################
# SAMPLE_NAME := Sample_NoIndex_L007
# BWAINDEX := /path/to/genome
###################
# REQUIRED MODULES
###################
# module load java
# module load samtools
# module load python
###################

SHELL = bash

# Use the default PATH version of each of these programs unless overridden
SAMTOOLS ?= samtools

FAI ?= $(BWAINDEX).fai

READ_LENGTH ?= 36

MIN_MAP_QUALITY ?= 30

# Filter to a minimum mapping quality, where both pairs
# are properly paired in the mapping
SAMTOOL_OPTIONS ?= -F 512

# Ideally this will be set to something else set in the environment or on the command line
TMPDIR ?= $(shell pwd)

# where our results files go
OUTDIR ?= $(shell pwd)

INBAM ?= $(OUTDIR)/$(SAMPLE_NAME).sorted.bam
OUTBAM ?= $(OUTDIR)/$(SAMPLE_NAME).uniques.sorted.bam
INSERTMETRICS ?= $(OUTDIR)/$(SAMPLE_NAME).CollectInsertSizeMetrics.picard

all : info metrics uniques

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

metrics : $(INSERTMETRICS)

uniques : $(INBAM).bai $(OUTBAM).bai

# Sometimes this will report errors about a read not mapping that should have a mapq of 0
# See this for more info: http://seqanswers.com/forums/showthread.php?t=4246
$(INSERTMETRICS) : $(INBAM) 
	time java -Xmx1000m -jar `which CollectInsertSizeMetrics.jar` INPUT=$^ OUTPUT=$@ \
                HISTOGRAM_FILE=$(INSERTMETRICS).pdf \
                VALIDATION_STRINGENCY=LENIENT \
                ASSUME_SORTED=true && echo Picard stats >&2

# Index uniquely mapping reads 
$(OUTBAM).bai : $(OUTBAM)
	time $(SAMTOOLS) index $^
                
# Sorted uniquely mapping reads BAM
$(OUTBAM) : $(INBAM)
	time $(SAMTOOLS) view $(SAMTOOL_OPTIONS) -b $(INBAM) > $@

# Index sorted BAM file
$(INBAM).bai : $(INBAM)
	time $(SAMTOOLS) index $^	
