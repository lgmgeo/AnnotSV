#!/bin/bash -x

set -eo pipefail

###################################################################
####                                                           ####
####       Check the AnnotSV_ID for INS with same POS          ####
####                                                           ####
###################################################################

mkdir -p ./output
rm -f ./output/*.annotated.tsv


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"


# Tests with 2 INS that have the same position but different SVLENs
###################################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test_INS_1.vcf" -outputFile "./output/test_INS_1.annotated.tsv"

result=`$cut "./output/test_INS_1.annotated.tsv" "AnnotSV_ID" | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$result" != "AnnotSV_ID 18_1963003_1963003_INS_1 18_1963003_1963003_INS_2 " ]
then
        echo "Error with test_INS_1"
        exit 1
fi

echo "test_INS_1 OK"


# Tests with 2 INS that have the same position but different inserted sequences
###############################################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test_INS_2.vcf" -outputFile "./output/test_INS_2.annotated.tsv"

result=`$cut "./output/test_INS_2.annotated.tsv" "AnnotSV_ID" | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$result" != "AnnotSV_ID 18_1963003_1963003_INS_1 18_1963003_1963003_INS_2 " ]
then
        echo "Error with test_INS_2"
        exit 1
fi

echo "test_INS_2 OK"



$ANNOTSVDEV/bin/AnnotSV -SVinputFile "./input/test_INS_3.vcf" -outputFile "./output/test_INS_3.annotated.tsv"

result=`$cut "./output/test_INS_3.annotated.tsv" "AnnotSV_ID;Annotation_mode" | grep "full"  | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$result" != "2_79197378_79198779_DEL_1 2_79197378_79198779_DEL_2 " ]
then
        echo "Error with test_INS_3"
        exit 1
fi

echo "test_INS_3 OK"


################################################################
################################################################
################################################################

echo "ok - Finished"

