echo "---"
echo "STAMPIPES:" 
echo "VERSION: " `jq '.version' $STAMPIPES/version.json`
echo "SHA: " `git --git-dir $STAMPIPES/.git rev-parse HEAD`
echo "---"

if [[ `command -v bedops` ]]; then
echo "BEDOPS:"
bedops --help | grep "version" | head -n1
echo "---"
fi

if [[ `command -v bwa` ]]; then
echo "BWA:"
bwa  2>&1 | grep "Version"
echo "---"
fi

if [[ `command -v samtools` ]]; then
echo "SAMTOOLS:"
samtools  2>&1 | grep "Version"
echo "---"
fi

if [[ `command -v fastqc` ]]; then
echo "FASTQC:"
fastqc --version
echo "---"
fi

CollectInsertSizeMetrics=`command -v CollectInsertSizeMetrics.jar`
MarkDuplicates=`command -v MarkDuplicates.jar`

if [[ $CollectInsertSizeMetrics || $MarkDuplicates ]]; then
echo "PICARD APPS: "
fi

if [[ $CollectInsertSizeMetrics ]]; then
echo "CollectInsertSizeMetrics"
java -jar $CollectInsertSizeMetrics  --version
fi

if [[ $MarkDuplicates ]]; then
echo "MarkDuplicates"
java -jar $MarkDuplicates --version
fi

if [[ $CollectInsertSizeMetrics || $MarkDuplicates ]]; then
echo "---"
fi

if [[ `command -v awk` ]]; then
echo "AWK:"
awk --version | head -n 1
echo "---"
fi

if [[ `command -v jq` ]]; then
echo "JQ:"
jq --version
echo "---"
fi

if [[ `command -v python` ]]; then
echo "PYTHON:" 
python -V 2>&1
python $STAMPIPES/scripts/versions.py
echo "---"
fi

if [[ `command -v java` ]]; then
echo "JAVA:"
java -version 2>&1
echo "---"
fi

if [[ `command -v bash` ]]; then
echo "BASH:"
bash --version
echo "---"
fi

if [[ `command -v R` ]]; then
echo "R:"
R --version | grep "R version"
fi
