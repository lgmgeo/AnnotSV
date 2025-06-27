#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"



mkdir -p ./output
rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 4 -genomeBuild CHM13 -outputFile "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 4 -genomeBuild CHM13 -outputFile "./output/test_tcl.annotated.tsv" -Tx ENSEMBL

# $cut "./output/test_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;AnnotSV_ranking_score;AnnotSV_ranking_criteria;ACMG_class"

annotations=`$cut "./output/test_tcl.annotated.tsv" "Annotation_mode;AnnotSV_ranking_score;AnnotSV_ranking_criteria;ACMG_class" | egrep -v "full=|AnnotSV_ranking_criteria"` 

if [ `echo $annotations | grep -c "full 0.9 1A (cf Gene_count, RE_gene, +0.00);2C-1 (SLC25A24, +0.90);3A (1 gene, +0.00);5F (+0.00) 4"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Dist_nearest_SS and Nearest_SS_type values for NM_001037668 are not correct (test1)"
        exit 1
fi


echo "ok - Finished"


