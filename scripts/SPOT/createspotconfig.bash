OUTDIR=$1
TAGS=$2
GENOME=$3
K=$4

CONFIGOUT=$OUTDIR/runall.tokens.txt

if [ -e $CONFIGOUT ];
    rm $CONFIGOUT
fi

echo "_TAGS_ = $TAGS" >> $CONFIGOUT
