###################
# This creates BigWig files from a BAM file for loading into a UCSC browser.
# 
###################
# SAMPLE_NAME=Example_NoIndex_L007
# BWAINDEX=/path/to/genome/hg19/hg19
# GENOME=hg19
# READLENGTH=36
###################
# Optional arguments
###################
# WIN=75
# BIN=20
# STAMPIPES=~/stampipes # Location of stampipes installation
###################
# REQUIRED MODULES
###################
# module load samtools
# module load python
# module load bedops
# module load bedtools
###################

FAI ?= $(BWAINDEX).fai
BAMFILE ?= $(SAMPLE_NAME).uniques.sorted.bam
STAMPIPES ?= ~/stampipes

TMPDIR ?= $(shell pwd)

# Window
WIN=75
# Bin interval
BINI=20
TAGEXP=$(FLOWCELL)-$(SAMPLE_NAME)-$(FORG)

BUCKETS_DIR ?= $(STAMPIPES)/data/densities/

CHROM_BUCKET ?= $(BUCKETS_DIR)/chrom-buckets.$(GENOME).$(WIN).$(BINI).bed.starch

all : $(SAMPLE_NAME).$(WIN).$(BINI).uniques-density.$(READLENGTH).$(GENOME).bed.starch $(SAMPLE_NAME).$(WIN)_$(BINI).$(GENOME).bw

# clip tags to single 5' end base
$(TMPDIR)/$(SAMPLE_NAME).m$(READLENGTH).uniques.$(BINI).bed : $(BAMFILE)
	bamToBed -i $^ \
    | awk '{if($$6=="+"){s=$$2; e=$$2+1}else{s=$$3-1; e=$$3}print $$1"\t"s"\t"e}' \
    | sort-bed - \
    > $@

$(SAMPLE_NAME).$(WIN).$(BINI).uniques-density.$(READLENGTH).$(GENOME).bed.starch : $(CHROM_BUCKET) $(SAMPLE_NAME).m$(READLENGTH).uniques.$(BINI).bed   
	unstarch $(CHROM_BUCKET) \
    | bedmap --faster --range $(WIN) --echo --count --delim "\t" - $(SAMPLE_NAME).m$(READLENGTH).uniques.$(BINI).bed \
    | starch - \
    > $@

$(TMPDIR)/$(SAMPLE_NAME).$(WIN)_$(BINI).$(GENOME).wig : $(SAMPLE_NAME).$(WIN).$(BINI).uniques-density.$(READLENGTH).$(GENOME).bed.starch
	unstarch $^ | awk -v binI=$(BINI) -f $(STAMPIPES)/awk/bedToWig.awk > $@
      
$(SAMPLE_NAME).$(WIN)_$(BINI).$(GENOME).bw : $(TMPDIR)/$(SAMPLE_NAME).$(WIN)_$(BINI).$(GENOME).wig
	wigToBigWig -clip $^ $(FAI) $@

