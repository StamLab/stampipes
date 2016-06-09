# !/bin/bash
# generate trackhub at the level of the flowcell

FLOW=$1
OUTPUTDIR_FC=/net/seq/data/flowcells/trackhubs/flowcells/
LOC=$(pwd)

# make directory for the concatenated flowcell hub
if [ ! -d flowcell_trackhub ];
then
    mkdir flowcell_trackhub
fi
PROJECTS=(browser-load*/)
for PROJECT in ${PROJECTS[*]}
do
    cd $PROJECT
    for GENOME in */
    do
        if [ ! -d ../flowcell_trackhub/$GENOME ];
	then
	    mkdir ../flowcell_trackhub/$GENOME
        fi
	echo "" > ../flowcell_trackhub/$GENOME/master_track.txt
    done
    cd ..
done

# for each projects/genome folder, add tracks to the master flowcell tracks
PROJECTS=(browser-load*/)
for PROJECT in ${PROJECTS[*]}
do
    cd $PROJECT
    for GENOME in */
    do
        # cat to master flowcell track
        GROUP=$(cat $GENOME\trackDb*.txt | grep "group" | head -1 | sed -e 's/group //g')
        cat $GENOME\trackDb*.txt | sed -e "0,/shortLabel/ s/shortLabel.*/shortLabel ${FLOW}_$GROUP/" >> ../flowcell_trackhub/$GENOME/master_track.txt
    done
    cd ..
done

cd flowcell_trackhub

# make new hub.txt file
echo "hub ${FLOW}_hub" > hub.txt
echo "shortLabel ${FLOW}_hub" >> hub.txt
echo "longLabel hub for all flowcell data" >> hub.txt
echo "genomesFile genomes.txt" >> hub.txt
echo "email placeholder@altiusinstitute.org" >> hub.txt

# make new genomes.txt from genomes seen
echo "" > genomes.txt
for GENOME in */
do

    GENOMENAME=$(echo "$GENOME" | sed 's/.$//')
    echo "genome $GENOMENAME" >> genomes.txt
    echo "trackDb ${GENOME}master_track.txt" >> genomes.txt
    echo "" >> genomes.txt

    # and create supertrack
    echo "track "$FLOW"_super" > $GENOME/super_track.txt
    echo "superTrack on" >> $GENOME/super_track.txt
    echo "group master" >> $GENOME/super_track.txt
    echo "shortLabel $FLOW" >> $GENOME/super_track.txt
    echo "longLabel all project tracks for $FLOW" >> $GENOME/super_track.txt
    echo "" >> $GENOME/super_track.txt
       
    cat $GENOME\master_track.txt | sed -e "s/compositeTrack on/parent ${FLOW}_super\ncompositeTrack on/g" >> $GENOME/super_track.txt

done

# create flowcell link in web-accessible location
if [ ! -L $OUTPUTDIR_FC$FLOWCELL ];
then
    ln -s $LOC/flowcell_trackhub $OUTPUTDIR_FC$FLOW
fi
