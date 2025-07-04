#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


mkdir -p ./output


# Check that no column is shifted
#################################

# RefSeq
rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 4 -genomeBuild CHM13 -outputFile "./output/test_tcl.annotated.tsv" -Tx RefSeq

for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "error1: (test39)"
                exit 1
        fi
done


# ENSEMBL
rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 4 -genomeBuild CHM13 -outputFile "./output/test_tcl.annotated.tsv" -Tx ENSEMBL

for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "error2: (test39)"
                exit 1
        fi
done


echo "ok - Finished"


