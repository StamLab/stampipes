# !/bin/bash

# create hub with given name from a text file listing all flowcells desired for that hub
# ex: update_hubs.sh june_2016 list_of_fc_from_june2016.txt

HUBNAME=$1
HUBFILE=$2

OUTPUTDIR_FC=/net/seq/data/flowcells/trackhubs/flowcells/
OUTPUTDIR_MASTERHUBS=/net/seq/data/flowcells/trackhubs/master_hubs/

# check if hub name exists
if [ -d $OUTPUTDIR_MASTERHUBS$HUBNAME ];
then
    echo "a hub with that name already exists, please remove existing directory or use a different name"
    exit;
else
    mkdir $OUTPUTDIR_MASTERHUBS$HUBNAME
fi

# create hub.txt file
echo "hub master_hub_${HUBNAME}" > ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/hub.txt
echo "shortLabel ${HUBNAME}" >> ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/hub.txt
echo "longLabel all tracks" >> ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/hub.txt
echo "genomesFile genomes.txt" >> ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/hub.txt
echo "email placeholder@altiusinstitute.org" >> ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/hub.txt

# create genome.txt file
echo "" > ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/genomes.txt

# for each flowcell / genome, make sure each genome is represented
while IFS= read line
do
    cd $OUTPUTDIR_FC$line
    for GENOME in */
    do
	if [ ! -d ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/${GENOME} ];
	then
	   GENOME=$(echo "$GENOME" | sed -e 's/\///g')
	   mkdir ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/$GENOME
	   echo "genome $GENOME" >> ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/genomes.txt
	   echo "trackDb $GENOME/super_track.txt" >> ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/genomes.txt
	   echo "" >> ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/genomes.txt
	fi
	cat ${OUTPUTDIR_FC}${line}/$GENOME/super_track.txt >> ${OUTPUTDIR_MASTERHUBS}${HUBNAME}/$GENOME/super_track.txt
    done
done < $HUBFILE