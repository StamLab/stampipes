#!/bin/bash  

if [ "$#" -lt "4" ] ; then
    echo "Usage:  $0  name file.bam strand-specificity refFlat_genePred_annotations" 1>&2
    exit
fi
name=$1
bam=$2
strandSpecificity=$3
REF_FLAT=$4


pushd `dirname $bam`
bam=`basename $bam`
picardfolder=$PICARDPATH

opspecs=( "InsertSize?HISTOGRAM_FILE=histo_insert_size_$name.pdf"  \
         "AlignmentSummary?" )
if [ -e "$REF_FLAT" ] ; then
    opspecs+=( "RnaSeq?REF_FLAT=$REF_FLAT STRAND_SPECIFICITY=$strandSpecificity CHART_OUTPUT=coverage_$name.pdf" )

fi
    #not used in this function, but should be:  "AlignmentSummary?INCLUDE_SECONDARY_ALIGNMENTS=false"
for opspec in "${opspecs[@]}"
do 
    op=`echo $opspec | cut -f1 -d'?'`
    params=`echo $opspec | cut -f2- -d'?'`
    output=picard.$name.$op.txt
    if ! [ -f  "$output" ] ; then
        operationjar=$picardfolder/Collect${op}Metrics.jar
        java -Xmx1000m -jar $operationjar INPUT=$bam OUTPUT=$output $params VALIDATION_STRINGENCY=SILENT
    fi
    if [ -s "$output" ] ; then
        cat $output | grep -v '^#' | grep -A4 _ | $SCRIPT_DIR/transposeTable.pl | $SCRIPT_DIR/encomma.pl  > xp.$output
    else
        echo "Failed to create $output" 1>&2
    fi
done


popd
exit

labels="TOTAL_READS READ_PAIRS MEAN_INSERT_SIZE STANDARD_DEVIATION MEDIAN_INSERT_SIZE PCT_CODING_BASES PCT_UTR_BASES PCT_MRNA_BASES PCT_INTRONIC_BASES PCT_INTERGENIC_BASES"
#header=`echo -e "dataset $labels" | tr ' ' '	'`
#echo -e "$header" > summarytable.txt
#for bam in `ls *.bam | grep 100.200` ; do
#    name=${bam%%.bam}
    line="$name"
    for label in $labels ; do
        value=`cat xp.picard.$name* | grep $label | awk '{print $2}'`
        echo "$label ($value)" 1>&2
        if [ "$value" == "" ] ; then
            value="missing"
        fi
        line="$line\t$value"
    done
    echo -e "$line" >> summarytable.txt
#done


#  
#  ERROR: Option 'INPUT' is required.
#  
#  Reads a SAM or BAM file and writes a file containing summary alignment metrics.
#  
#  
#  Options:
#  
#  --help
#  -h                            Displays options specific to this tool.
#  
#  --stdhelp
#  -H                            Displays options specific to this tool AND options common to all Picard command line 
#                                tools.
#  
#  MAX_INSERT_SIZE=Integer       Paired end reads above this insert size will be considered chimeric along with 
#                                inter-chromosomal pairs.  Default value: 100000. This option can be set to 'null' to 
#                                clear the default value. 
#  
#  ADAPTER_SEQUENCE=String       This option may be specified 0 or more times. This option can be set to 'null' to clear 
#                                the default list. 
#  
#  IS_BISULFITE_SEQUENCED=Boolean
#  BS=Boolean                    Whether the SAM or BAM file consists of bisulfite sequenced reads.    Default value: 
#                                false. This option can be set to 'null' to clear the default value. Possible values: 
#                                {true, false} 
#  
#  INPUT=File
#  I=File                        Input SAM or BAM file.  Required. 
#  
#  OUTPUT=File
#  O=File                        File to write the output to.  Required. 
#  
#  REFERENCE_SEQUENCE=File
#  R=File                        Reference sequence fasta  Default value: null. 
#  
#  ASSUME_SORTED=Boolean
#  AS=Boolean                    If true (default), then the sort order in the header file will be ignored.  Default 
#                                value: true. This option can be set to 'null' to clear the default value. Possible 
#                                values: {true, false} 
#  
#  STOP_AFTER=Integer            Stop after processing N reads, mainly for debugging.  Default value: 0. This option can 
#                                be set to 'null' to clear the default value. 
#  

#INCLUDE_SECONDARY_ALIGNMENTS=Boolean    If false, do not write secondary alignments to output. Default value: true. This option can be set to 'null' to clear the default value. Possible values: {true, false} 


