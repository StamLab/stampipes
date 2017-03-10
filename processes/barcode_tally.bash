# Dependencies
source $MODULELOAD
module load python/2.7.11

MAKEFILE="$STAMPIPES/makefiles/barcodes.mk"
bash "$STAMPIPES/scripts/versions.bash" &> barcode.versions.txt

make -q -f "$MAKEFILE" report \
    || qmake -cwd -N .tb$FLOWCELL -V -now no -- -f "$MAKEFILE" -j 40
