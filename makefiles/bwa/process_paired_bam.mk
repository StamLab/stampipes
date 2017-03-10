###################
# These variables must be passed in or set for the makefile to work.  If the genome's
# FAI file is not at $(BWAINDEX).fai, then it must also be specified under FAI.
###################
# SAMPLE_NAME := Sample_NoIndex_L007
# BWAINDEX := /path/to/genome
###################
# REQUIRED MODULES
###################
# module load jdk
# module load picard
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
INSERTMETRICS ?= $(OUTDIR)/$(SAMPLE_NAME).CollectInsertSizeMetrics.picard

all : info metrics uniques indices

single_ended : info uniques indices

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

indices : $(INBAM).bai $(OUTBAM).bai

uniques : $(INBAM) $(OUTBAM)

NODE_RAM_GB ?= 8
NSLOTS ?= 1
JAVA_HEAP = $(shell echo "($(NODE_RAM_GB)*$(NSLOTS) - 2) * 1024" | bc)m
# Sometimes this will report errors about a read not mapping that should have a mapq of 0
# See this for more info: http://seqanswers.com/forums/showthread.php?t=4246
$(INSERTMETRICS) : $(OUTBAM) 
	time java -Xmx$(JAVA_HEAP) -jar $(PICARDPATH)/CollectInsertSizeMetrics.jar INPUT=$^ OUTPUT=$@ \
                HISTOGRAM_FILE=$(INSERTMETRICS).pdf \
                VALIDATION_STRINGENCY=LENIENT \
                ASSUME_SORTED=true && echo Picard stats >&2

# Index uniquely mapping reads 
%.bam.bai : %.bam
	time $(SAMTOOLS) index $^
                
# Sorted uniquely mapping reads BAM
$(OUTBAM) : $(INBAM)
	time $(SAMTOOLS) view $(SAMTOOL_OPTIONS) -b $(INBAM) > $@
