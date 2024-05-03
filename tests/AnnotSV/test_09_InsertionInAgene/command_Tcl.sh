#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"

# Check the annotation of an insertion in the BBS1 gene:
########################################################
# 11      66278204        66278205        104252245       M        <INS>  M       M       M
 

rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtBEDcol 6 -genomeBuild GRCh37


if [ `$cut "./output/test_tcl.annotated.tsv" "Gene_name" | grep -c BBS1` != 2 ]
then
        echo "error1: Insertion is not annotated in the BBS1 gene."
        exit 1
fi

echo "ok - Finished"

