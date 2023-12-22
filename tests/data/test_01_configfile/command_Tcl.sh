#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


# - Check that users can select only a subset of the annotation columns provided by AnnotSV by using the $ANNOTSV/etc/AnnotSV/configfile
#   => decrease the number of output annotation column
# - Check that no column is shifted by using the configfile



# Run of AnnotSV, asking for the minimum of annotation columns (Using a "configfile" located in the same directory as the input file)
cp $ANNOTSV/etc/AnnotSV/configfile.minAnnotation ./input/configfile
rm -f "./output/test_minAnn_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 5  -SVinputInfo 1 -outputFile "./output/test_minAnn_tcl.annotated.tsv" -genomeBuild GRCh37
rm ./input/configfile

minAnn=`head -1 ./output/test_minAnn_tcl.annotated.tsv | wc -w`
nLinesForMinAnn=`wc -l ./output/test_minAnn_tcl.annotated.tsv | awk '{print $1}'`


# Run of AnnotSV, with all annotation columns asked
rm -f "./output/test_allAnn_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 5  -SVinputInfo 1 -outputFile "./output/test_allAnn_tcl.annotated.tsv" -genomeBuild GRCh37

allAnn=`head -1 "./output/test_allAnn_tcl.annotated.tsv" | wc -w`
nLinesForAllAnn=`wc -l "./output/test_allAnn_tcl.annotated.tsv" | awk '{print $1}'`


if [ $nLinesForMinAnn != $nLinesForAllAnn ] 
then
	echo "error 1: not the same number of lines in both .tsv files"
	exit 1
fi

if [[ $minAnn > 71 ]]
then
        echo "error 1: too much annotation columns ($minAnn)"
        exit 1
fi

if [ $minAnn >= $allAnn ]
then
        echo "error 1: minimum of annotation columns not respected"
        exit 1
fi



# Check that no column is shifted
for v in `$cut "./output/test_minAnn_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "error1: some columns semms to be shifted"
                exit 1
        fi
done


echo "ok - Finished"


