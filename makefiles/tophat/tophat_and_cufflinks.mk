SHELL=/bin/bash

BOWTIE ?= bowtie

THREADS ?= 1

TOPHAT ?= tophat
TOPHAT_REF ?= $(GENOME)
REF_SEQ ?= $(TOPHAT_REF)
LIBTYPE ?= "fr-firststrand" #This is used for alignment
STRAND_SPEC ?= "SECOND_READ_TRANSCRIPTION_STRAND" # Used for metrics

ANNOT_GTF ?= $(REF_DIR)/genes.gtf
ANNOT_GENEPRED ?= $(REF_DIR)/refFlat.txt
CHROM_SIZES ?= $(REF_DIR)/chrom_sizes.txt

TMPDIR ?= $(shell pwd)

SAMTOOLS ?= samtools
JAVA ?= java

# This might not work in qmake??
MARK_DUPS ?= $(PICARDPATH)/MarkDuplicates.jar

SCRIPT_DIR ?= $(STAMPIPES)/scripts/tophat

RIBOSOMAL_REF ?= $(REF_DIR)/contamination/hg_rRNA
CONTROL_REF   ?= $(REF_DIR)/spikeInControlRNA/ERCC92


ADAPTER_FILE ?= $(SAMPLE_NAME).adapters.txt

TRIM_DIR ?= $(TMPDIR)/trimmed
R1_trimmed_fastq ?= $(addprefix $(TRIM_DIR)/,$(notdir $(R1_FASTQ)))
R2_trimmed_fastq ?= $(addprefix $(TRIM_DIR)/,$(notdir $(R2_FASTQ)))

ribosomal_file ?= $(addprefix $(TMPDIR)/, $(addsuffix .ribosomalRNA.txt, $(notdir $(R1_FASTQ))))
control_file   ?= $(addprefix $(TMPDIR)/, $(addsuffix .spikeInControlRNA.bowtie.txt, $(notdir $(R1_FASTQ))))
tophat_file    ?= $(TMPDIR)/tophat/accepted_hits.bam

ribo_txt = $(SAMPLE_NAME).rRNAcounts.txt
control_txt = $(SAMPLE_NAME).spikeInControlCounts.txt
readcount_txt = $(SAMPLE_NAME).read_counts.txt
bamcount_txt = $(SAMPLE_NAME).bam_counts.txt

summary_txt = $(SAMPLE_NAME).sample_summary.txt
upload_txt = $(SAMPLE_NAME).rna_metrics.txt

marked_bam = $(SAMPLE_NAME).all.$(GENOME).bam

cufflinks_finished = $(SAMPLE_NAME)_cufflinks/finished.txt
coverage_types = all pos neg
strand_prefix = $(addsuffix .$(GENOME), $(addprefix $(SAMPLE_NAME)., $(coverage_types)))
strand_bam = $(addsuffix .bam, $(strand_prefix))
bai = $(addsuffix .bam.bai, $(strand_prefix))
bigwig = $(addsuffix .bw, $(strand_prefix))
starch = $(addsuffix .starch, $(strand_prefix))

.DELETE_ON_ERROR:

.SUFFIXES:

.INTERMEDIATE: $(tophat_file) $(R1_trimmed_fastq) $(R2_trimmed_fastq) $(ribosomal_file) $(control_file)

.PHONY: default all summary cufflinks coverage ribosomal alignment upload

#### Targets

default: all

all: upload cufflinks coverage alignment

summary: $(summary_txt)

cufflinks: $(cufflinks_finished)

coverage: $(bigwig) $(starch)

ribosomal: $(ribosomal_files) $(ribo_txt)

alignment: $(marked_bam) $(strand_bam) $(bai)

upload: $(upload_txt)
	python3 $(STAMPIPES)/scripts/lims/upload_data.py --rnafile $^ --alignment_id $(ALIGNMENT_ID)


#### Coverage files

%.bw : $(TMPDIR)/%.bed
	bedGraphToBigWig $^ $(CHROM_SIZES) $@

%.starch : $(TMPDIR)/%.bed
	starch $^ > $@

$(TMPDIR)/%.bed : %.bam
	bam2bed --split < $^ | $(SCRIPT_DIR)/bed_coverage.pl > $@

%.bam.bai : %.bam
	samtools index $^

#### Strand-specific BAM

%.neg.$(GENOME).bam : %.all.$(GENOME).bam
	samtools merge -f $@ \
	<( samtools view -h  -u -f 0x90 -F 0x40 $^) \
	<( samtools view     -u -F 0x90 -f 0x40 $^) \

%.pos.$(GENOME).bam : %.all.$(GENOME).bam
	samtools merge -f $@ \
	<( samtools view -h  -u -f 0x50 -F 0x80 $^) \
	<( samtools view     -u -F 0x50 -f 0x80 $^) \

# Alternate rules
%.neg.bam : %.all.bam
	samtools merge -f $@ \
	<( samtools view -h  -u -f 0x90 -F 0x40 $^) \
	<( samtools view     -u -F 0x90 -f 0x40 $^) \

%.pos.bam : %.all.bam
	samtools merge -f $@ \
	<( samtools view -h  -u -f 0x50 -F 0x80 $^) \
	<( samtools view     -u -F 0x50 -f 0x80 $^) \


$(cufflinks_finished) : $(marked_bam)
	$(SCRIPT_DIR)/cufflinks.sh $(marked_bam) $(REF_SEQ) $(LIBTYPE) $(ANNOT_GTF) $(dir $@) $(THREADS) GTF && \
	touch $@

#### Counts

$(upload_txt) : $(summary_txt)
	$(SCRIPT_DIR)/transposeTable.pl $^ > $@

$(summary_txt) : $(bamcount_txt) $(ribo_txt) $(readcount_txt)
	cat $^ | $(SCRIPT_DIR)/processStatsNameKeyValue.pl > $@

# BAM counts
$(bamcount_txt) : $(marked_bam)
	$(SCRIPT_DIR)/summarizeBAM.sh $(SAMPLE_NAME) $< $(STRAND_SPEC) $(ANNOT_GENEPRED) ; \
	$(SCRIPT_DIR)/mappedBamStats.pl $(marked_bam) $(SAMPLE_NAME) >> $@ ; \
	$(SCRIPT_DIR)/collectStatsFromTransposedPicardTables.pl . $(SAMPLE_NAME) >> $@

# Read counts:
$(readcount_txt) : $(R1_FASTQ) $(R2_FASTQ)
	$(SCRIPT_DIR)/countReads.pl $(SAMPLE_NAME) total_reads $^ > $@

# Ribosomal counts:
$(ribo_txt) : $(ribosomal_file)
	$(SCRIPT_DIR)/countLines.pl $(SAMPLE_NAME) rRNA $^ > $@

# Spike-in counts:
$(control_txt) : $(control_file)
	cat $^ | cut -f3 | $(SCRIPT_DIR)/countFew.pl > $@

#### Alignment

# Final (all) BAM
$(marked_bam) $(SAMPLE_NAME).spotdups.txt : $(tophat_file)
	$(JAVA) -Xmx24g -jar $(MARK_DUPS) \
		INPUT=$< \
		METRICS_FILE=$(SAMPLE_NAME).spotdups.txt \
		OUTPUT=$(marked_bam) \
		REMOVE_DUPLICATES=false \
		ASSUME_SORTED=true \
		READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'

$(tophat_file) : $(R1_trimmed_fastq) $(R2_trimmed_fastq)
	$(SCRIPT_DIR)/tophatPE.sh $^ $(TOPHAT_REF) $(LIBTYPE) $(ANNOT_GTF) $(TMPDIR)/tophat \
		$(THREADS)

$(ribosomal_file) : $(R1_trimmed_fastq)
	zcat -f $^ | $(BOWTIE) --threads $(THREADS) -n 3 -e 140 \
		$(RIBOSOMAL_REF) - $@

$(control_file) : $(R1_trimmed_fastq)
	zcat -f $^ | $(BOWTIE) --threads $(THREADS) -n 3 -e 140 \
		$(CONTROL_REF) - $@

$(R1_trimmed_fastq) $(R2_trimmed_fastq) : $(R1_FASTQ) $(R2_FASTQ)
	mkdir $(TRIM_DIR) -p ; \
	trim-adapters-illumina \
		-f $(ADAPTER_FILE) \
		-1 P5 -2 P7 \
		$+ \
		$(R1_trimmed_fastq) $(R2_trimmed_fastq)

