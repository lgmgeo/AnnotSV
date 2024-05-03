#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"



# Special case
##############
# 2 SVs exist with the same coordinates but with a different type (DEL, DUP):
# 	1       3777548 4234731 tt      DEL
# 	1       3777548 4234731 tt      DUP

# => We do have 2 different SVs!
# We must have full + split annotations for these 2 SVs.

rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -svtBEDcol 5 -outputFile "./output/test_tcl.annotated.tsv" -genomeBuild GRCh37

# $cut "./output/test_tcl.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;;SV_type;Annotation_mode;Gene_name"
# AnnotSV_ID      SV_chrom        SV_start        SV_end  SV_length       SV_type Annotation_mode Gene_name
# 1_3777548_4234731_DEL_1 1       3777548 4234731 -457183 DEL     full    C1orf174;DFFB;LINC01134;LINC01345;LINC01346
# 1_3777548_4234731_DEL_1 1       3777548 4234731 -457183 DEL     split   C1orf174
# 1_3777548_4234731_DEL_1 1       3777548 4234731 -457183 DEL     split   DFFB
# 1_3777548_4234731_DEL_1 1       3777548 4234731 -457183 DEL     split   LINC01134
# 1_3777548_4234731_DEL_1 1       3777548 4234731 -457183 DEL     split   LINC01345
# 1_3777548_4234731_DEL_1 1       3777548 4234731 -457183 DEL     split   LINC01346
# 1_3777548_4234731_DUP_1 1       3777548 4234731 457183  DUP     full    C1orf174;DFFB;LINC01134;LINC01345;LINC01346
# 1_3777548_4234731_DUP_1 1       3777548 4234731 457183  DUP     split   C1orf174
# 1_3777548_4234731_DUP_1 1       3777548 4234731 457183  DUP     split   DFFB
# 1_3777548_4234731_DUP_1 1       3777548 4234731 457183  DUP     split   LINC01134
# 1_3777548_4234731_DUP_1 1       3777548 4234731 457183  DUP     split   LINC01345
# 1_3777548_4234731_DUP_1 1       3777548 4234731 457183  DUP     split   LINC01346

if [ `$cut "./output/test_tcl.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;Annotation_mode;Gene_name" | grep -c "1_3777549_4234731_DEL_1"` == 6 ]
then
        echo "Ok"
else
        echo "error 1: 1_3777548_4234731_DEL_1 is not correctly reported"
        exit 1
fi

if [ `$cut "./output/test_tcl.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;Annotation_mode;Gene_name" | grep -c "1_3777549_4234731_DUP_1"` == 6 ]
then
        echo "Ok"
else
        echo "error 1: 1_3777548_4234731_DUP_1 is not correctly reported"
        exit 1
fi


echo "ok - Finished"

