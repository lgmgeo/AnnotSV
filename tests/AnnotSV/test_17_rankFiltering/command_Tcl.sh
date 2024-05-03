#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}



#########################################
# Check the "-rankFiltering" option
#########################################

# -> 1  3 4 5 (pas de 2)
rm -f "./output/test.1_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test.1_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37
for v in `$cut "./output/test.1_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 3 4 5 ACMG_class full=1 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with ./input/test.bed"
                echo "error 1"
                exit 1
        fi
done


# -> 1
rm -f "./output/test.2_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test.2_tcl.annotated.tsv" -svtBEDcol 4 -rankFiltering "1-2" -genomeBuild GRCh37
for v in `$cut "./output/test.2_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 ACMG_class full=1 full=2" " " "$v"
        then
                echo "Error with ./input/test.bed"
                echo "error 1"
                exit 1
        fi
done


# -> 1 3
rm -f "./output/test.3_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test.3_tcl.annotated.tsv" -svtBEDcol 4 -rankFiltering "1-3" -genomeBuild GRCh37
for v in `$cut "./output/test.3_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 3 ACMG_class full=1 full=3" " " "$v"
        then
                echo "Error with ./input/test.bed"
                echo "error 1"
                exit 1
        fi
done


# -> 1 3
rm -f "./output/test.4_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test.4_tcl.annotated.tsv" -svtBEDcol 4 -rankFiltering "1,2,3" -genomeBuild GRCh37
for v in `$cut "./output/test.4_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 3 ACMG_class full=1 full=3" " " "$v"
        then
                echo "Error with ./input/test.bed"
                echo "error 1"
                exit 1
        fi
done


# -> 1 3 4
rm -f "./output/test.5_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test.5_tcl.annotated.tsv" -svtBEDcol 4 -rankFiltering "1,2-4" -genomeBuild GRCh37
for v in `$cut "./output/test.5_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 3 4 ACMG_class full=1 full=3 full=4" " " "$v"
        then
                echo "Error with ./input/test.bed"
                echo "error 1"
                exit 1
        fi
done



echo "ok - Finished"

