#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"


# Check with a bed containing 2 SV on MT chromosome
# input/test.bed
# M       649     1603
# M       987     1611

rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -genomeBuild GRCh37

test=`$cut "./output/test_tcl.annotated.tsv" "AnnotSV type"`

if [ `echo $test | grep -c "No column to display"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: AnnotSV should not annotate Sv on MT chromosome"
        exit 1
fi


echo "ok - Finished"

