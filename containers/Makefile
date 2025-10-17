# Makefile for building BWA process containers
#
# Usage:
#   make all           - Build all containers
#   make align         - Build just the align container
#   make clean         - Remove all .sif files
#   make list          - List all available targets

# Container definitions and their corresponding .sif files
CONTAINERS = align samtools-1.12 pysam-samtools mark_duplicates insert_size spot_score density_files macs2 bwa-aggregate-base

# Directory where containers are stored
CONTAINER_DIR = bwa

# Default target
.PHONY: all
all: $(addprefix $(CONTAINER_DIR)/, $(addsuffix .sif, $(CONTAINERS)))

# Pattern rule for building .sif files from .def files
$(CONTAINER_DIR)/%.sif: $(CONTAINER_DIR)/%.def
	@echo "Building container: $@"
	apptainer build --force "$@" "$<"
	@echo "Successfully built: $@"

# Individual container targets for convenience
.PHONY: align samtools-1.12 pysam-samtools mark_duplicates insert_size spot_score density_files macs2 bwa-aggregate-base
align: $(CONTAINER_DIR)/align.sif
samtools-1.12: $(CONTAINER_DIR)/samtools-1.12.sif
pysam-samtools: $(CONTAINER_DIR)/pysam-samtools.sif
mark_duplicates: $(CONTAINER_DIR)/mark_duplicates.sif
insert_size: $(CONTAINER_DIR)/insert_size.sif
density_files: $(CONTAINER_DIR)/density_files.sif
spot_score: $(CONTAINER_DIR)/spot_score.sif
macs2: $(CONTAINER_DIR)/macs2.sif
bwa-aggregate-base: $(CONTAINER_DIR)/bwa-aggregate-base.sif

# Clean target to remove all .sif files
.PHONY: clean
clean:
	@echo "Removing all .sif files..."
	rm -f $(CONTAINER_DIR)/*.sif
	@echo "Clean complete."

# List available targets
.PHONY: list
list:
	@echo "Available container targets:"
	@echo "  all              - Build all containers"
	@echo "  align            - Build align container"
	@echo "  samtools-1.12    - Build samtools-1.12 container"
	@echo "  pysam-samtools   - Build pysam-samtools container"
	@echo "  mark_duplicates  - Build mark_duplicates container"
	@echo "  insert_size      - Build insert_size container"
	@echo "  density_files    - Build density_files container"
	@echo "  spot_score       - Build spot_score container"
	@echo "  macs2            - Build macs2 container"
	@echo "  bwa-aggregate-base - Build bwa-aggregate-base container"
	@echo "  clean            - Remove all .sif files"
	@echo "  list             - Show this help"

# Help target
.PHONY: help
help: list
