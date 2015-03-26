TALLY_SCRIPT ?= $(STAMPIPES)/scripts/flowcells/tallybarcodes.sh
REPORT_SCRIPT ?= $(STAMPIPES)/scripts/flowcells/barcode_report.py

PROCESSING ?= $(FLOWCELL_DIR)/processing.json

SAMPLE_FASTQ = $(shell ls $(FLOWCELL_DIR)/Project_*/Sample*/*_R1_???.fastq.gz)
UNDET_FASTQ = $(shell ls $(FLOWCELL_DIR)/Undetermined_indices/Sample*/*_R1_???.fastq.gz)
FASTQ ?= $(SAMPLE_FASTQ) $(UNDET_FASTQ)
barcodes = $(patsubst %.fastq.gz,%.barcodes.txt,$(FASTQ))
report = $(FLOWCELL_DIR)/barcode_report.txt

all : info barcodes report

info :
	@echo "------" && \
	 echo "PROGRAM VERSIONS" && \
	 echo "------" && \
	 echo "METADATA" && \
	 echo "------" && \
	 echo "FLOWCELL_NAME: " $(FLOWCELL) && \
	 echo "------" && \
	 echo "FASTQ_COUNT: " $(words $(FASTQ)) && \
	 echo "------"

report: $(report)
	
$(report): $(barcodes)
	SGE_RREQ=" -N .tbr$(FLOWCELL)-report " python "$(REPORT_SCRIPT)" -p "$(PROCESSING)" -b "$(FLOWCELL_DIR)" > $@

barcodes : $(barcodes)
	

%.barcodes.txt : %.fastq.gz
	SGE_RREQ=" -N .tb$(FLOWCELL)-fasta " bash $(TALLY_SCRIPT) $^ > $@
