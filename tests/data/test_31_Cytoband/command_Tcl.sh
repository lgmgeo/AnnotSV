#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"


# Integration of cytoband annotations:
######################################

# Request sur GitHub:
# https://github.com/lgmgeo/AnnotSV/issues/52

# Input SV
# chr17	1	18960000	DEL

# GRCh37 cytoband file:
# chr17   0       3300000 p13.3   gneg
# chr17   3300000 6500000 p13.2   gpos50
# chr17   6500000 10700000        p13.1   gneg
# chr17   10700000        16000000        p12     gpos75
# chr17   16000000        22200000        p11.2   gneg

# Expected annotation:
# p13.3-p11.2


# GRCh37:
#########
rm -f "./output/test.GRCh37_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 4 -outputFile "./output/test.GRCh37_tcl.annotated.tsv" -overlap 1 -genomeBuild GRCh37

# $cut "./output/test.GRCh37_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;CytoBand" | gr full
# 17_1_18960000_DEL_1     full    p13.3-p11.2
annotations=`$cut "./output/test.GRCh37_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;CytoBand" | grep full`
if [ `echo $annotations | grep -c "p13.3-p11.2"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: cytoband annotations are not correct (GRCh37)"
        exit 1
fi


# GRCh38:
#########
rm -f "./output/test.GRCh38_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 4 -outputFile "./output/test.GRCh38_tcl.annotated.tsv" -overlap 1 -genomeBuild GRCh38

# $cut "./output/test.GRCh38_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;CytoBand" | gr full
# 17_1_18960000_DEL_1     full    p13.3-p11.2
annotations=`$cut "./output/test.GRCh38_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;CytoBand" | grep full`
if [ `echo $annotations | grep -c "p13.3-p11.2"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: cytoband annotations are not correct (GRCh38)"
        exit 1
fi


echo "ok - Finished"


