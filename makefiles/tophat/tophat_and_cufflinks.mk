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

SAMTOOLS ?= samtools
JAVA ?= java

# This might not work in qmake??
MARK_DUPS ?= $(shell which MarkDuplicates.jar)

SCRIPT_DIR ?= $(STAMPIPES)/scripts/tophat

RIBOSOMAL_REF ?= $(REF_DIR)/contamination/hg_rRNA
CONTROL_REF   ?= $(REF_DIR)/spikeInControlRNA/ERCC92

R1_FASTQ_FILES ?= $(wildcard $(SAMPLE_NAME)_R1_*.fastq.gz)
R2_FASTQ_FILES ?= $(wildcard $(SAMPLE_NAME)_R2_*.fastq.gz)

ADAPTERFILE ?= $(SAMPLE_NAME).adapters.txt

TRIM_DIR ?= trimmed
R1_trimmed_fastq_files ?= $(prefix $(TRIM_DIR)/$(R1_FASTQ_FILES))
R2_trimmed_fastq_files ?= $(prefix $(TRIM_DIR)/$(R2_FASTQ_FILES))

ribosomal_files ?= $(patsubst $(SAMPLE_NAME)_R1_%.fastq.gz,ribosomalRNA.$(SAMPLE_NAME)_%.bowtie.txt,$(R1_FASTQ_FILES))
control_files     ?= $(patsubst $(SAMPLE_NAME)_R1_%.fastq.gz,spikeInControlRNA.$(SAMPLE_NAME)_%.bowtie.txt,$(R1_FASTQ_FILES))
tophat_files    ?= $(patsubst $(SAMPLE_NAME)_R1_%.fastq.gz,$(SAMPLE_NAME)_tophat_%/accepted_hits.bam,$(R1_FASTQ_FILES))

ribo_txt = $(SAMPLE_NAME).rRNAcounts.txt
control_txt = $(SAMPLE_NAME).spikeInControlCounts.txt
readcount_txt = $(SAMPLE_NAME).read_counts.txt
bamcount_txt = $(SAMPLE_NAME).bam_counts.txt

summary_txt = $(SAMPLE_NAME).sample_summary.txt
upload_txt = $(SAMPLE_NAME).rna_metrics.txt

marked_bam = $(SAMPLE_NAME).all.$(GENOME).bam
cufflinks_finished = $(SAMPLE_NAME)_cufflinks/finished.txt
coverage_types = all pos neg
bigwig = $(addsuffix .$(GENOME).bw, $(addprefix $(SAMPLE_NAME)., $(coverage_types)))
starch = $(addsuffix .$(GENOME).starch, $(addprefix $(SAMPLE_NAME)., $(coverage_types)))

.DELETE_ON_ERROR:

.SECONDARY:

.SUFFIXES:

.INTERMEDIATE: %.bed


default: all

all: upload $(cufflinks_finished) $(bigwig)

upload: $(upload_txt)
	python $(STAMPIPES)/scripts/lims/upload_data.py --rnafile $^ --alignment_id $(ALIGNMENT_ID)

$(upload_txt) : $(summary_txt)
	$(SCRIPT_DIR)/transposeTable.pl $^ > $@

$(summary_txt) : $(bamcount_txt) $(ribo_txt) $(readcount_txt)
	cat $^ | $(SCRIPT_DIR)/processStatsNameKeyValue.pl > $@

%.bw : %.bed
	bedGraphToBigWig $^ $(CHROM_SIZES) $@

%.bed : %.starch
	bedops --ec -u $^ | $(SCRIPT_DIR)/singleBedFileBaseCoverage.sh | $(SCRIPT_DIR)/compressBed4.pl > $@


%.$(GENOME).starch : %.$(GENOME).bam
	samtools view -u $^ | bedtools bamtobed -split -i stdin | cut -f1-3 | sort-bed - | starch - > $@


%.pos.$(GENOME).bam : %.all.$(GENOME).bam
	samtools merge $@ \
	<( samtools view -h  -u -f 0x90 -F 0x40 $^) \
	<( samtools view     -u -F 0x90 -f 0x40 $^) \

%.neg.$(GENOME).bam : %.all.$(GENOME).bam
	samtools merge $@ \
	<( samtools view -h  -u -f 0x50 -F 0x80 $^) \
	<( samtools view     -u -F 0x50 -f 0x80 $^) \



$(cufflinks_finished) : $(marked_bam)
	 $(SCRIPT_DIR)/cufflinks.sh $(marked_bam) $(REF_SEQ) $(LIBTYPE) $(ANNOT_GTF) $(dir $@) $(THREADS) GTF && \
	touch $@

# Mark dups
$(marked_bam) $(SAMPLE_NAME).spotdups.txt : $(SAMPLE_NAME).tophat_merged.bam
	$(JAVA) -Xmx24g -jar $(MARK_DUPS) \
		INPUT=$< \
		METRICS_FILE=$(SAMPLE_NAME).spotdups.txt \
		OUTPUT=$(marked_bam) \
		REMOVE_DUPLICATES=false \
		ASSUME_SORTED=true

# Merge alignments together
$(SAMPLE_NAME).tophat_merged.bam : $(tophat_files)
	$(SCRIPT_DIR)/merge_or_copy_bam.sh $@ $^

# Read counts:
$(readcount_txt) : $(R1_FASTQ_FILES) $(R2_FASTQ_FILES)
	$(SCRIPT_DIR)/countReads.pl $(SAMPLE_NAME) total_reads . > $@

# Ribosomal counts:
$(ribo_txt) : $(ribosomal_files)
	$(SCRIPT_DIR)/countLines.pl $(SAMPLE_NAME) rRNA $^ > $@

# Spike-in counts:
$(control_txt) : $(control_files)
	cat $^ | cut -f3 | $(SCRIPT_DIR)/countFew.pl > $@

$(bamcount_txt) : $(marked_bam)
	$(SCRIPT_DIR)/summarizeBAM.sh $(SAMPLE_NAME) $< $(STRAND_SPEC) $(ANNOT_GENEPRED) ; \
	$(SCRIPT_DIR)/mappedBamStats.pl $(marked_bam) $(SAMPLE_NAME) >> $@ ; \
	$(SCRIPT_DIR)/collectStatsFromTransposedPicardTables.pl . $(SAMPLE_NAME) >> $@

# File-specific processing
# TODO: These rules are hecka ugly. There's probably a nice make-way to write
# the filetype transformation, but my make-fu is weak.

$(SAMPLE_NAME)_tophat_%/accepted_hits.bam : $(TRIM_DIR)/$(SAMPLE_NAME)_R1_%.fastq.gz $(TRIM_DIR)/$(SAMPLE_NAME)_R2_%.fastq.gz
	 $(SCRIPT_DIR)/tophatPE.sh $^ $(TOPHAT_REF) $(LIBTYPE) $(ANNOT_GTF) $(SAMPLE_NAME)_tophat_$* \
		$(THREADS)

ribosomalRNA.$(SAMPLE_NAME)_%.bowtie.txt : $(TRIM_DIR)/$(SAMPLE_NAME)_R1_%.fastq.gz
	 zcat -f $^ | $(BOWTIE) --threads $(THREADS) -n 3 -e 140 \
		$(RIBOSOMAL_REF) - $@

spikeInControlRNA.$(SAMPLE_NAME)_%.bowtie.txt : $(TRIM_DIR)/$(SAMPLE_NAME)_R1_%.fastq.gz
	 zcat -f $^ | $(BOWTIE) --threads $(THREADS) -n 3 -e 140 \
		$(CONTROL_REF) - $@

$(TRIM_DIR)/$(SAMPLE_NAME)_R1_%.fastq.gz $(TRIM_DIR)/$(SAMPLE_NAME)_R2_%.fastq.gz : $(SAMPLE_NAME)_R1_%.fastq.gz $(SAMPLE_NAME)_R2_%.fastq.gz
	mkdir $(TRIM_DIR) -p ; \
	trim-adapters-illumina \
	  -f $(ADAPTERFILE) \
	  -1 P5 -2 P7 \
	  $+ \
	  $(TRIM_DIR)/$(SAMPLE_NAME)_R1_$*.fastq.gz $(TRIM_DIR)/$(SAMPLE_NAME)_R2_$*.fastq.gz

