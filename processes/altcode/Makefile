.PHONY: all
all: containers

.PHONY: containers
containers: scanpy.sif

%.sif: %.apptainer.def %.environment.yaml
	apptainer build --force "$@" "$<"
