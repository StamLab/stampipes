# Dependencies
source $MODULELOAD
module load python/2.7.3

MAKEFILE="$STAMPIPES/makefiles/barcodes.mk"
bash "$STAMPIPES/scripts/versions.bash" &> barcode.versions.txt

make -q -f "$MAKEFILE" report \
    || qmake -cwd -N .tb$FLOWCELL -V -q all.q -now no -- -f "$MAKEFILE" -j 40
