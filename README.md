This is a focused collection of makefiles and scripts to for analysis on Illumina sample data. While this probably can't be directly used by someone not in our lab, it might be a useful reference for others.

Set up
-------------

In order to set up these scripts, checkout the repository.  There are some environment variables to set:

* STAMPIPES -- full path to the installation directory
* LIMS_API_URL -- URL to the LIMS API
* LIMS_API_TOKEN -- Auth to the LIMS API
* MODULELOAD -- Full path to the GNU modules init script
* HOTSPOT_DIR -- Path to HOTSPOT, defaults to `~/hotspot/hotspot-distr` - download at [https://github.com/rthurman/hotspot]

A python version that includes pip needs to be installed for proper python package version recording.

Modules
-------------

The current modules to load for working the pipeline are:

	module load bwa/0.7.0
	module load samtools/0.1.19
	module load bedops/2.4.2
	module load bedtools/2.16.2
	module load python/2.7.3
	module load java/jdk1.7.0_05
	module load picard/1.118

If $MODULELOAD is set, then the processing scripts should load the required modules automatically.
