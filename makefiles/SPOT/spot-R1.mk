###################
# This only samples from R1.
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
# module load java
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

SPOT_OUT ?= $(OUTDIR)/$(SAMPLE_NAME).R1.rand.uniques.sorted.spot.out
DUP_OUT ?= $(OUTDIR)/$(SAMPLE_NAME).R1.rand.uniques.sorted.spotdups.txt

RANDOM_SAMPLE_BAM ?= $(TMPDIR)/$(SAMPLE_NAME).R1.rand.uniques.sorted.bam

calcspot : $(SPOT_OUT)
calcdup : $(DUP_OUT)

# Only use Read 1 from uniquely mapping reads; exclude chrM and random 
# Note: awk needs to have $$ to escape make's interpretation
$(TMPDIR)/$(SAMPLE_NAME).R1.uniques.sorted.bam : $(BAMFILE)
	samtools view -f 0x0040 $(SAMPLE_NAME).uniques.sorted.bam | \
		awk '{if($$3 != "chrM" && $$3 != "random"){print}}' | \
		samtools view -bS -t $(FAI) - > \
		$(TMPDIR)/$(SAMPLE_NAME).R1.uniques.sorted.bam

# Make a random sample from the filtered BAM 
$(RANDOM_SAMPLE_BAM) : $(TMPDIR)/$(SAMPLE_NAME).R1.uniques.sorted.bam
	bash -e $(STAMPIPES)/scripts/SPOT/randomsample.bash $(SAMPLE_SIZE) $(FAI) $^ $@

$(SPOT_OUT) : $(SPOTDIR)/$(SAMPLE_NAME).R1.rand.uniques.sorted.spot.out
	cp $(SPOTDIR)/$(SAMPLE_NAME).R1.rand.uniques.sorted.spot.out $(SPOT_OUT)

# run the SPOT program
$(SPOTDIR)/$(SAMPLE_NAME).R1.rand.uniques.sorted.spot.out : $(RANDOM_SAMPLE_BAM)
	bash -e $(STAMPIPES)/scripts/SPOT/runhotspot.bash $(HOTSPOT_DIR) $(SPOTDIR) $(RANDOM_SAMPLE_BAM) $(GENOME) $(READLENGTH) $(ASSAY)

# picard has trouble calculating dups when we are just working with a sampling
# from R1 on paired end data; converting back and forth with BED will
# strip that information
$(RANDOM_SAMPLE_BAM).bed : $(RANDOM_SAMPLE_BAM)
	bamToBed -i $^ > $@

$(RANDOM_SAMPLE_BAM).bed.sorted.bam : $(RANDOM_SAMPLE_BAM).bed
	bedToBam -i $^ -g $(FAI) > $(RANDOM_SAMPLE_BAM).bed.bam
	samtools sort $(RANDOM_SAMPLE_BAM).bed.bam $(RANDOM_SAMPLE_BAM).bed.sorted

# Calculate the duplication score of the random sample
$(DUP_OUT) : $(RANDOM_SAMPLE_BAM).bed.sorted.bam
	java -jar `which MarkDuplicates.jar` INPUT=$(RANDOM_SAMPLE_BAM).bed.sorted.bam OUTPUT=$(TMPDIR)/$(SAMPLE_NAME).R1.rand.uniques.dup \
	  METRICS_FILE=$(OUTDIR)/$(SAMPLE_NAME).R1.rand.uniques.sorted.spotdups.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT
