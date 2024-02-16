#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


# Test with a large SV (105 Mb)
###############################
# (overlapping too many TAD boundaries => previously, caused a bug in excel display (split the line with a "\n")
# chr1:2806107-107058351


rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtBEDcol 6 -genomeBuild GRCh37


# Check that no column is shifted
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "error1: Error with Large SV. Columns shifted."
                exit 1
        fi
done


echo "ok - Finished"

