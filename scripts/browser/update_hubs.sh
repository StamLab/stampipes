# !/bin/bash

# keep master hubs locations updated

FLOWCELL=$1
OUTPUTDIR=/net/seq/data/flowcells/
OUTPUTDIR_MASTERS=/net/seq/data/flowcells/trackhubs/masters/
OUTPUTDIR_FC=/net/seq/data/flowcells/trackhubs/flowcells/

# this is the list of names to create the master track from
buildfile="master_flowcell_list.txt"

# create flowcell link and add to master hub list
if [ ! -L $OUTPUTDIR_FC$FLOWCELL ];
then
    LOC=$(pwd)
    ln -s $LOC/flowcell_trackhub $OUTPUTDIR_FC$FLOWCELL
    echo "$FLOWCELL" >> $OUTPUTDIR_MASTERS$buildfile
fi

# make master hubs.txt if doesn't exist
if [ ! -f ${OUTPUTDIR_MASTERS}hub.txt ];
then
    echo "hub master_hub" > ${OUTPUTDIR_MASTERS}hub.txt
    echo "shortLabel MasterHub" >> ${OUTPUTDIR_MASTERS}hub.txt
    echo "longLabel all tracks" >> ${OUTPUTDIR_MASTERS}hub.txt
    echo "genomesFile genomes.txt" >> ${OUTPUTDIR_MASTERS}hub.txt
    echo "email placeholder@altiusinstitute.org" >> ${OUTPUTDIR_MASTERS}hub.txt
fi

# make master genomes.txt if it doesn't exist
if [ ! -f ${OUTPUTDIR_MASTERS}genomes.txt ];
then 
   echo "" > ${OUTPUTDIR_MASTERS}genomes.txt
fi

# for each flowcell / genome, make sure each genome is represented and remove any existing tracks
while IFS= read line
do
    cd $OUTPUTDIR_FC$line
    for GENOME in */
    do
        if [ ! -d $OUTPUTDIR_MASTERS$GENOME ]
        then
	    GENOME=$(echo "$GENOME" | sed -e 's/\///g')
	    mkdir $OUTPUTDIR_MASTERS$GENOME
	    echo "genome $GENOME" >> ${OUTPUTDIR_MASTERS}genomes.txt
	    echo "trackDb $GENOME/super_track.txt" >> ${OUTPUTDIR_MASTERS}genomes.txt
	    echo "" >> ${OUTPUTDIR_MASTERS}genomes.txt
        else
	    if [ -f $OUTPUTDIR_MASTERS$GENOME/super_track.txt ];
	    then
		echo "" > $OUTPUTDIR_MASTERS$GENOME/super_track.txt
	    fi
	fi
    done
    cd ..
done < $OUTPUTDIR_MASTERS$buildfile

# now for each flowcell / genome, go and add the super tracks togethers
while IFS= read line
do
    cd $OUTPUTDIR_FC$line
    for GENOME in */
    do
        GENOME=$(echo "$GENOME" | sed -e 's/\///g')
        cat $OUTPUTDIR_FC$line/$GENOME/super_track.txt >> $OUTPUTDIR_MASTERS$GENOME/super_track.txt
    done
    cd ..
done < $OUTPUTDIR_MASTERS$buildfile