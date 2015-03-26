###################
# This creates the chromosome buckets to use for making densities to make bigWig files
# for the UCSC browser.
# Needs to be given either BWAINDEX with a .fai or a FAI file.
###################
# BWAINDEX=/path/to/genome/hg19/hg19
# GENOME=hg19
###################
# REQUIRED MODULES
###################
# module load bedops
###################

FAI ?= $(BWAINDEX).fai
# WINDOW
WIN ?= 75
# BIN INTERVAL
BINI ?= 20
BUCKETS_DIR ?= $(STAMPIPES)/data/densities
BUCKETS_FILE ?= $(BUCKETS_DIR)/chrom-buckets.$(GENOME).$(WIN)_$(BINI).bed.starch

TMPDIR ?= $(shell pwd)

all : $(BUCKETS_FILE)

$(BUCKETS_FILE) :
	awk -v w=$(WIN) '{print $$1"\t0\t"$$2-w}' $(FAI) | sort-bed - \
  | awk -v binI=$(BINI) -v win=$(WIN) '{ for(i = $$2 + win; i < $$3; i += binI) { print $$1"\t"i - win"\t"i + win }}' \
  | starch - > $@
