#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


# # BED file séparé par des espaces (au lieu de tabulation) : space-separated-input.bed
# #####################################################################################
# 
# $ANNOTSV/bin/AnnotSV -SVinputFile space-separated-input.bed -SVinputInfo 1 -outputFile ./space-separated-input_tcl.annotated.tsv -svtBEDcol -1 -genomeBuild GRCh37
# => plante :
# -- refGeneAnnotation --
# bedtools intersect -sorted -a ./space-separated-input.formatted.sorted.bed -b /home/geoffroy/analysisOfSV/AnnotSV/AnnotSV_dev/share/AnnotSV/Annotations_Human/RefGene/GRCh37/refGene.sorted.bed -wa -wb > ./space-separ
# ated-input_tcl.annotated.tsv.tmp.tmp
# Error: unable to open file or unable to determine types for file ./space-separated-input.formatted.sorted.bed
# Exit with error



# TSV: Check that no column is shifted
######################################
rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/tab-separated.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37 

for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with ./input/test.bed"
                echo "error 1"
                exit 1
        fi
done

# Space separated value: Error message
######################################
# The above "set -e" option instructs bash to immediately exit if any command has a non-zero exit status.
# So here, as the "AnnotSV" command fails, the "echo" command will run anyway thanks to the "||"
rm -f "./output/space-separated_tcl.annotated.tsv"
checks=`$ANNOTSV/bin/AnnotSV -SVinputFile "./input/space-separated.bed" -SVinputInfo 1 -outputFile "./output/space-separated_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37 || echo "AnnotSV error"`

if [ `echo $checks | grep -c "has non positional records, which are only valid for the groupBy tool"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: \"non positional records\" error message NOT found"
        exit 1
fi

rm -f space-separated.NA.bed

echo "ok - Finished"


