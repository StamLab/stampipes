###################
# These variables must be passed in or set for the makefile to work.
###################
# SAMPLE_NAME=Example_NoIndex_L007
###################
# If the data is paired, then indicate this
###################
# PAIRED=True
###################
# If uploading data to the LIMS, then indicate the LIMS_API_TOKEN to use and FLOWCELL uploading name
# and the ALIGNMENT_ID for the counts
###################
# LIMS_API_URL=http://lims.com
# LIMS_API_TOKEN=aethua87ao997i9a9u
# ALIGNMENT_ID=1435
# FLOWCELL=CU8TE
###################

# Use the default PATH version of each of these programs unless overridden
FASTQC ?= fastqc

# Use four threads by default; each thread uses 250MB of memory, so we are staying
# well within qsub memory limits
FASTQC_OPTIONS ?= -t 4

INDIR ?= $(shell pwd)

OUTDIR ?= $(shell pwd)

STAMPIPES ?= ~/stampipes

ifdef PAIRED
 TARGETS ?= $(OUTDIR)/$(SAMPLE_NAME)_R1_fastqc $(OUTDIR)/$(SAMPLE_NAME)_R2_fastqc
else
 TARGETS ?= $(OUTDIR)/$(SAMPLE_NAME)_R1_fastqc 
endif 

NUM_FILES = $(shell ls $(INDIR)/$(SAMPLE_NAME)_R1_???.fastq.gz | wc -l)

.PHONY : upload fastqc info

all : info fastqc upload

fastqc : $(TARGETS) 
	@echo "Targets: " $(TARGETS)
	
info :
	@echo "------"
	@echo "PROGRAM VERSIONS"
	@echo "------"
	@echo "FastQC:"
	@$(FASTQC) -v
	@echo "------"
	@echo "METADATA"
	@echo "------"
	@echo "SAMPLE_NAME: " $(SAMPLE_NAME)
	@echo "------"
	@echo "PROGRAM OPTIONS"
	@echo "------"
	@echo "FASTQC_OPTIONS: " $(FASTQC_OPTIONS)
	@echo "------"

# Make fastqc results
# Also move the folder if we only have one file pair so we don't have the extra 001
$(OUTDIR)/$(SAMPLE_NAME)_%_fastqc :
	@echo FastQ files: $(NUM_FILES)
	time $(FASTQC) $(FASTQC_OPTIONS) $(INDIR)/$(SAMPLE_NAME)_$*_???.fastq.gz --casava && echo FastQC stats >&2
ifeq ($(NUM_FILES), 1)
	mv $(OUTDIR)/$(SAMPLE_NAME)_$*_001_fastqc $(OUTDIR)/$(SAMPLE_NAME)_$*_fastqc
	mv $(OUTDIR)/$(SAMPLE_NAME)_$*_001_fastqc.zip $(OUTDIR)/$(SAMPLE_NAME)_$*_fastqc.zip
endif

ifdef LIMS_API_TOKEN
upload : $(TARGETS)
	@echo "Uploading FastQC data"
ifdef PAIRED
	@echo "Paired end upload"
	python $(STAMPIPES)/scripts/lims/upload_data.py -f $(FLOWCELL) --alignment_id $(ALIGNMENT_ID) --flowcell_lane_id=$(FLOWCELL_LANE_ID) \
	  --fastqcfile $(OUTDIR)/$(SAMPLE_NAME)_R1_fastqc/fastqc_data.txt --fastqcfile $(OUTDIR)/$(SAMPLE_NAME)_R2_fastqc/fastqc_data.txt \
	  --fastqc_counts
else
	@echo "Single end upload"
	python $(STAMPIPES)/scripts/lims/upload_data.py -f $(FLOWCELL) --alignment_id $(ALIGNMENT_ID) --flowcell_lane_id=$(FLOWCELL_LANE_ID) \
	  --fastqcfile $(OUTDIR)/$(SAMPLE_NAME)_R1_fastqc/fastqc_data.txt \
	  --fastqc_counts
endif
endif
