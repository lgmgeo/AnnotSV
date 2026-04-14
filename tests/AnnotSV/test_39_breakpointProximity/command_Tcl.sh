#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


# Check that user can place an annotation BED file in ".../Users/*/BNDproximity/" and use this BNDproximity mode
################################################################################################################

mkdir -p "./output"
rm -f "./output/*.annotated.tsv"

# Run of AnnotSV with a file in ".../Users/GRCh38/BNDproximity/"
cp ./input/B.bed $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh38/BNDproximity/
$ANNOTSV/bin/AnnotSV -SVinputFile ./input/A.bed -outputFile ./output/A.annotated.tsv -svtBEDcol 4 -genomeBuild GRCh38
rm $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh38/BNDproximity/B.formatted.sorted.bed $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh38/BNDproximity/B.header.tsv


result=`$cut "./output/A.annotated.tsv" "AnnotSV_ID;B_annotations_BNDprox" | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2 }'`
if [ "$result" != "AnnotSV_ID B_annotations_BNDprox 1_101_200_DEL_1 B_ann2 1_151_170_DUP_1 B_ann1;B_ann2 " ]
then
        echo "Error with BNDproximity mode"
        exit 1
fi

# Run of AnnotSV with a file in ".../Users/GRCh38/BNDproximity/", and with the "-breakpointProximity" option
cp ./input/B.bed $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh38/BNDproximity/
$ANNOTSV/bin/AnnotSV -SVinputFile ./input/A.bed -outputFile ./output/A.annotated.tsv -svtBEDcol 4 -genomeBuild GRCh38 -breakpointProximity 60
rm $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh38/BNDproximity/B.formatted.sorted.bed $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh38/BNDproximity/B.header.tsv


result=`$cut "./output/A.annotated.tsv" "AnnotSV_ID;B_annotations_BNDprox" | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2 }'`
if [ "$result" != "AnnotSV_ID B_annotations_BNDprox 1_101_200_DEL_1  1_151_170_DUP_1 B_ann2 " ]
then
        echo "Error with the "-breakpointProximity" option"
        exit 1
fi



# Check that no column is shifted
for v in `$cut "./output/A.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "error1: some columns semms to be shifted"
                exit 1
        fi
done


echo "ok - Finished"


