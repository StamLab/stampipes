This is a focused collection of makefiles and scripts to for analysis on Illumina sample data. While this probably can't be directly used by someone not in our lab, it might be a useful reference for others.

Set up
-------------

In order to set up these scripts, checkout the repository.  There are some environment variables to set:

* STAMPIPES -- full path to the installation directory
* LIMS_API_URL -- URL to the LIMS API
* LIMS_API_TOKEN -- Auth to the LIMS API
* MODULELOAD -- Full path to the GNU modules init script
* HOTSPOT_DIR -- Path to HOTSPOT, defaults to `~/hotspot/hotspot-distr` - download at [https://github.com/rthurman/hotspot]
* PYTHON3_ACTIVATE -- A virtualenv activation script for python3.

A python version that includes pip needs to be installed for proper python package version recording.

Modules
-------------

The current modules to load for working the pipeline are:

* `bcl2fastq/1.8.4`
* `bcl2fastq2/2.15.0.4`
* `bedops/2.4.19`
* `bedtools/2.25.0`
* `bowtie/1.0.0`
* `bwa/0.7.12`
* `coreutils/8.25`
* `cufflinks/2.2.1`
* `fastqc/0.11.5`
* `gcc/4.7.2`
* `git/2.3.3`
* `jdk/1.8.0_92`
* `picard/1.120`
* `pigz/2.3.3`
* `pysam/0.9.0`
* `python/2.7.11`
* `python/3.5.1`
* `R/3.2.5`
* `RSEM/1.2.30`
* `samtools/1.3`
* `STAR/2.4.2a`
* `tophat/2.0.13`

If $MODULELOAD is set, then the processing scripts should load the required modules automatically.
