#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"



###################################
# Annotation method used in AnnotSV
###################################

# In each section, only 1 criterion (from the most pathogenic to the least) is assigned.

# Example in the section 2:
###########################
# https://github.com/lgmgeo/AnnotSV/issues/41
# https://github.com/lgmgeo/AnnotSV/issues/40

# 1_108190708_108194629_DEL

# In the following example, criterion 2F is filled but not assigned because criterion 2C-1 is filled and more pathogenic.

# 2C-1
# => The SLC25A24 gene is an OMIM morbid gene. 
# => Regarding the DEL used as input, the 2C-1 criteria is verified (coding sequence of the morbid gene is involved).

# 2F: Completely contained within an established benign CNV region.

rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 4 -genomeBuild GRCh38 -outputFile "./output/test_tcl.annotated.tsv"

# $cut "./output/test_tcl.annotated.tsv" "AnnotSV_ID;Annotation_mode;AnnotSV_ranking_score;AnnotSV_ranking_criteria;ACMG_class"
# AnnotSV_ID      Annotation_mode AnnotSV_ranking_score   AnnotSV_ranking_criteria        ACMG_class
# 1_108190709_108194629_DEL_1     full    0.9     1A (cf Gene_count, RE_gene, +0.00);2C-1 (SLC25A24, +0.90);3A (1 gene, +0.00);5F (+0.00) 4
# 1_108190709_108194629_DEL_1     split                   full=4

annotations=`$cut "./output/test_tcl.annotated.tsv" "Annotation_mode;AnnotSV_ranking_score;AnnotSV_ranking_criteria;ACMG_class" | egrep -v "full=|AnnotSV_ranking_criteria"` 

if [ `echo $annotations | grep -c "full 0.9 1A (cf Gene_count, RE_gene, +0.00);2C-1 (SLC25A24, +0.90);3A (1 gene, +0.00);5F (+0.00) 4"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: Dist_nearest_SS and Nearest_SS_type values for NM_001037668 are not correct (test1)"
        exit 1
fi


echo "ok - Finished"


