###################
# These variables must be passed in or set for the makefile to work.  If the genome's
# FAI file is not at $(BWAINDEX).fai, then it must also be specified under FAI.
###################
# SAMPLE_NAME=Example_NoIndex_L007
# BWAINDEX=/path/to/genome/hg19/hg19
# GENOME=hg19
# READLENGTH=36
# ASSAY=DNaseI
###################
# REQUIRED MODULES
###################
# module load jdk
# module load picard
# module load samtools
# module load python
# module load bedops
# module load bedtools
###################

FAI ?= $(BWAINDEX).fai
SAMPLE_SIZE ?= 5000000
BAMFILE ?= $(SAMPLE_NAME).uniques.sorted.bam
STAMPIPES ?= ~/stampipes
HOTSPOT_DIR ?= ~/hotspot/hotspot-distr

TMPDIR ?= $(shell pwd)
OUTDIR ?= $(shell pwd)

SPOTDIR ?= $(TMPDIR)/$(SAMPLE_NAME)_spot_R1

all : calcdup calcspot

SPOT_OUT ?= $(OUTDIR)/$(SAMPLE_NAME).rand.uniques.sorted.spot.out
DUP_OUT ?= $(OUTDIR)/$(SAMPLE_NAME).rand.uniques.sorted.spotdups.txt

NUCLEAR_SAMPLE_BAM ?= $(TMPDIR)/$(SAMPLE_NAME).nuclear.uniques.sorted.bam
RANDOM_SAMPLE_BAM ?= $(TMPDIR)/$(SAMPLE_NAME).rand.uniques.sorted.bam

calcspot : $(SPOT_OUT)
calcdup : $(DUP_OUT)

# exclude chrM*, chrC and random
# Note: awk needs to have $$ to escape make's interpretation
$(NUCLEAR_SAMPLE_BAM) : $(BAMFILE)
	samtools view $< \
		| awk '{if( ! index($$3, "chrM") && $$3 != "chrC" && $$3 != "random"){print}}' \
		| samtools view -uS -t $(FAI) - \
		> $@

# Make a random sample from the filtered BAM
$(RANDOM_SAMPLE_BAM) : $(NUCLEAR_SAMPLE_BAM)
	bash -e $(STAMPIPES)/scripts/SPOT/randomsample.bash $(SAMPLE_SIZE) $(FAI) $^ $@

$(SPOT_OUT) : $(SPOTDIR)/$(SAMPLE_NAME).rand.uniques.sorted.spot.out
	cp $(SPOTDIR)/$(SAMPLE_NAME).rand.uniques.sorted.spot.out $(SPOT_OUT)

# run the SPOT program
$(SPOTDIR)/$(SAMPLE_NAME).rand.uniques.sorted.spot.out : $(RANDOM_SAMPLE_BAM)
	bash -e $(STAMPIPES)/scripts/SPOT/runhotspot.bash $(HOTSPOT_DIR) $(SPOTDIR) $(RANDOM_SAMPLE_BAM) $(GENOME) $(READLENGTH) $(ASSAY)

# Calculate the duplication score of the random sample
$(DUP_OUT) : $(RANDOM_SAMPLE_BAM)
	java -jar $(PICARDPATH)/MarkDuplicates.jar INPUT=$(RANDOM_SAMPLE_BAM) OUTPUT=/dev/null \
	  METRICS_FILE=$(DUP_OUT) ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
