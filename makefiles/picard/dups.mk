###################
# Creates a duplication score from a sorted BAM file.
###################
# REQUIRED VARIABLES
###################
# SAMPLE_NAME=Example_NoIndex_L007
###################
# OPTIONAL VARIABLES
# Note: SAMPLE_NAME is not necessary if BAMFILE and DUP_OUT are specified.
###################
# BAMFILE : The input BAM files (default $(SAMPLE_NAME).sorted.bam)
# DUP_OUT : The output file (default $(SAMPLE_NAME).sorted.bam.picard.dup)
# OUTBAM : The output bam (default: /dev/null)
###################
# REQUIRED MODULES
###################
# module load java
###################

BAMFILE ?= $(SAMPLE_NAME).sorted.bam
DUP_OUT ?= $(SAMPLE_NAME).MarkDuplicates.picard
OUTBAM ?= /dev/null

all : $(DUP_OUT)

# Calculate the duplication score of the random sample
$(DUP_OUT) : $(BAMFILE)
	java -jar `which MarkDuplicates.jar` INPUT=$(BAMFILE) OUTPUT=$(OUTBAM) \
	  METRICS_FILE=$(DUP_OUT) ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
