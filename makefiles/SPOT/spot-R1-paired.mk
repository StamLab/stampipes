###################
# This only samples from R1.  Duplicates are calculated from paired data.
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
SAMPLE_SIZE ?= 5000000 # Number of fragments (read pairs)

INDIR ?= $(shell pwd)

BAMFILE ?= $(INDIR)/$(SAMPLE_NAME).uniques.sorted.bam
STAMPIPES ?= ~/stampipes
HOTSPOT_DIR ?= ~/hotspot/hotspot-distr

OUTDIR ?= $(shell pwd)
TMPDIR ?= $(OUTDIR)

SPOTDIR ?= $(TMPDIR)/$(SAMPLE_NAME)_spot_R1

all : calcdup calcspot

# Files produced by hotspot have this prefix
SPOTPREFIX=$(SAMPLE_NAME).R1.rand.uniques.sorted

SPOT_OUT ?= $(OUTDIR)/$(SPOTPREFIX).spot.out
DUP_OUT ?= $(OUTDIR)/$(SAMPLE_NAME).rand.uniques.sorted.spotdups.txt

SPOT_INFO ?= $(OUTDIR)/$(SPOTPREFIX).spot.info

PROPERLY_PAIRED_BAM ?= $(TMPDIR)/$(SAMPLE_NAME).properlypaired.sorted.bam
RANDOM_SAMPLE_BAM ?= $(TMPDIR)/$(SAMPLE_NAME).rand.uniques.sorted.bam
RANDOM_SAMPLE_BAM_R1 ?= $(TMPDIR)/$(SPOTPREFIX).bam


# Files produced by hotspot
HOTSPOT_SPOT = $(SPOTDIR)/$(SPOTPREFIX).spot.out
HOTSPOT_WIG = $(SPOTDIR)/$(SPOTPREFIX)-both-passes/$(SPOTPREFIX).hotspot.twopass.zscore.wig
HOTSPOT_STARCH = $(SPOTDIR)/$(SPOTPREFIX).hotspots.starch

calcspot : $(SPOT_OUT) $(SPOT_INFO)
calcdup : $(DUP_OUT)

$(RANDOM_SAMPLE_BAM) : $(BAMFILE)
	samtools view -h -F 12 -f 3 $^ \
		| awk '{if( ! index($$3, "chrM") && $$3 != "chrC" && $$3 != "random"){print}}' \
		| samtools view -uS - \
		> $(PROPERLY_PAIRED_BAM)
	bash $(STAMPIPES)/scripts/bam/random_sample.sh $(PROPERLY_PAIRED_BAM) $@ $(SAMPLE_SIZE)

# Only use Read 1 from our sample for SPOT score
$(RANDOM_SAMPLE_BAM_R1) : $(RANDOM_SAMPLE_BAM)
	samtools view -f 0x0040 $^ | samtools view -bSt $(FAI) - > $@

$(SPOT_OUT) : $(SPOTDIR)/$(SAMPLE_NAME).R1.rand.uniques.sorted.spot.out
	cp $(SPOTDIR)/$(SAMPLE_NAME).R1.rand.uniques.sorted.spot.out $(SPOT_OUT)

# run the SPOT program
$(HOTSPOT_SPOT) : $(RANDOM_SAMPLE_BAM_R1)
	bash -e $(STAMPIPES)/scripts/SPOT/runhotspot.bash $(HOTSPOT_DIR) $(SPOTDIR) $(RANDOM_SAMPLE_BAM_R1) $(GENOME) $(READLENGTH) $(ASSAY)

$(SPOT_INFO) : $(HOTSPOT_STARCH) $(HOTSPOT_SPOT)
	$(STAMPIPES)/scripts/SPOT/info.sh $(HOTSPOT_STARCH) hotspot1 $(HOTSPOT_SPOT) > $@

$(HOTSPOT_STARCH) : $(HOTSPOT_WIG)
	starch --header $(HOTSPOT_WIG) > "$@"

# Dummy rule
$(HOTSPOT_WIG) : $(HOTSPOT_SPOT)
	@

# Calculate the duplication score of the random sample
$(DUP_OUT) : $(RANDOM_SAMPLE_BAM)
	java -jar $(PICARDPATH)/MarkDuplicates.jar INPUT=$(RANDOM_SAMPLE_BAM) OUTPUT=$(TMPDIR)/$(SAMPLE_NAME).R1.rand.uniques.dup \
		METRICS_FILE=$(OUTDIR)/$(SPOTPREFIX).spotdups.txt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
		READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
