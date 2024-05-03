#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}



# Test of the BED file example available at the AnnotSV website (=/home/geoffroy/www/AnnotSV/web/Documentation/test.bed)
########################################################################################################################

rm -f "./input/test.bed" "./output/test_tcl.annotated.tsv"
cp "/home/geoffroy/www/AnnotSV/web/Documentation/test.bed" "./input"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37

# Look at the last column of the file (SV ranking) and check that there is no shift:

for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 3 4 5 ACMG_class full=1 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with the BED file example available at the AnnotSV website: shift in the ranking"
                echo "error 1"
                exit 1
        fi
done



echo "ok - Finished"

