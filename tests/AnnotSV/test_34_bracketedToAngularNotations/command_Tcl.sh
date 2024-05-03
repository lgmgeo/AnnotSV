#!/bin/bash -x

set -eo pipefail

###################################################################
####                                                           ####
#### Compare square/angle-bracketed notations and BED notation ####
####                                                           ####
###################################################################

rm -f ./output/command_*.annotated.tsv
rm -f ./output/command_*.unannotated.tsv


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"


# command_1: del_1
##################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-del_1.vcf" -outputFile "./output/command_1_square-bracketed-del_1.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-del_1.vcf" -outputFile "./output/command_1_angle-bracketed-del_1.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-del_1.bed" -svtBEDcol 5 -outputFile "./output/command_1_BED-del_1.annotated.tsv"

# tail -2 ./input/input-square-bracketed-del_1.vcf
# 12      3000    breakend_del_1_a        T       T[12:5000[      .       PASS    SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_1_b     GT      0/0     0/1
# 12      5000    breakend_del_1_b        T       ]12:3000]T      .       PASS    SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_1_a     GT      0/0     0/1
 
# tail -1 ./input/input-angle-bracketed-del_1.vcf
# 12      3000    breakend_del_1_a        T       <DEL>   .       PASS    END=5000;EXTRA=DEL_PAIRED;MATEID=breakend_del_1_b       GT      0/0     0/1

# tail -1 ./input/input-del_1.bed
# 12      2999    5000    breakend_del_1_a        DEL

# gr breakend_del_1_ Rodrigo_tests_ignore/output.csv
# breakend_del_1_a,12,3000,12,5000,T,T[12:5000[,2000,N[,DEL

squareResult=`$cut "./output/command_1_square-bracketed-del_1.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$squareResult" != "12_3000_5000_DEL_1 breakend_del_1_a 12 3000 5000 -2000 DEL T T[12:5000[ 3 " ]
then
        echo "error 1: Error with command_1 (square bad values)"
        exit 1
fi

echo "Command_1 (square) OK"

angleResult=`$cut "./output/command_1_angle-bracketed-del_1.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$angleResult" != "12_3000_5000_DEL_1 breakend_del_1_a 12 3000 5000 -2000 DEL T <DEL> 3 " ]
then
        echo "error 1: Error with command_1 (angle bad values)"
        exit 1
fi

echo "Command_1 (angle) OK"

BEDresult=`$cut "./output/command_1_BED-del_1.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$BEDresult" != "12_3000_5000_DEL_1 12 3000 5000 -2000 DEL 3 " ]
then
        echo "error 1: Error with command_1 (BED bad values)"
        exit 1
fi

echo "Command_1 (BED) OK"



# command_2: del_partial_pass (unordered)
#########################################

$ANNOTSV/bin/AnnotSV -SVinputFile "input/input-square-bracketed-del_partial_pass.vcf" -outputFile "./output/command_2_square-bracketed-del_partial_pass.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "input/input-angle-bracketed-del_partial_pass.vcf" -outputFile "./output/command_2_angle-bracketed-del_partial_pass.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "input/input-del_partial_pass.bed" -svtBEDcol 5 -outputFile "./output/command_2_BED-del_partial_pass.annotated.tsv"

# tail -2 ./input/input-square-bracketed-del_partial_pass.vcf
# 15      5000    breakend_del_partial_pass_a     T       ]15:3000]T      .       FAIL    SVTYPE=BND;EXTRA=DEL_PAIRED_PARTIAL_PASS;MATEID=breakend_del_partial_pass_b     GT      0/0     0/1
# 15      3000    breakend_del_partial_pass_b     T       T[15:5000[      .       PASS    SVTYPE=BND;EXTRA=DEL_PAIRED_PARTIAL_PASS;MATEID=breakend_del_partial_pass_a     GT      0/0     0/1

# tail -1 ./input/input-angle-bracketed-del_partial_pass.vcf
# 15      3000    breakend_del_partial_pass_b     T       <DEL>   .       PASS    EXTRA=DEL_PAIRED_PARTIAL_PASS;MATEID=breakend_del_partial_pass_a;END=5000       GT      0/0     0/1

# tail -1 ./input/input-del_partial_pass.bed
# 15      2999    5000    breakend_del_partial_pass_b     DEL

# gr breakend_del_partial_pass_ Rodrigo_tests_ignore/output.csv
# breakend_del_partial_pass_b,15,3000,15,5000,T,T[15:5000[,2000,N[,DEL

squareResult=`$cut "./output/command_2_square-bracketed-del_partial_pass.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$squareResult" != "15_3000_5000_DEL_1 breakend_del_partial_pass_b 15 3000 5000 -2000 DEL T T[15:5000[ 3 " ]
then
        echo "error 1: Error with command_2 (square bad values)"
        exit 1
fi

echo "Command_2 (square) OK"

angleResult=`$cut "./output/command_2_angle-bracketed-del_partial_pass.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$angleResult" != "15_3000_5000_DEL_1 breakend_del_partial_pass_b 15 3000 5000 -2000 DEL T <DEL> 3 " ]
then
        echo "error 1: Error with command_2 (angle bad values)"
        exit 1
fi

echo "Command_2 (angle) OK"

BEDresult=`$cut "./output/command_2_BED-del_partial_pass.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$BEDresult" != "15_3000_5000_DEL_1 15 3000 5000 -2000 DEL 3 " ]
then
        echo "error 1: Error with command_2 (BED bad values)"
        exit 1
fi

echo "Command_2 (BED) OK"



# command_3: Square-bracketed del_unordered
###########################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-del_unordered.vcf" -outputFile "./output/command_3_square-bracketed-del_unordered.annotated.tsv"

# tail -2 ./input/input-square-bracketed-del_unordered.vcf
# 22      5000    breakend_del_unordered_b        T       ]22:3000]T      .       PASS    SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_unordered_a     GT      0/0     0/1
# 22      3000    breakend_del_unordered_a        T       T[22:5000[      .       PASS    SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_unordered_b     GT      0/0     0/1

# gr breakend_del_unordered_ Rodrigo_tests_ignore/output.csv
# breakend_del_unordered_a,22,3000,22,5000,T,T[22:5000[,2000,N[,DEL


squareResult=`$cut "./output/command_3_square-bracketed-del_unordered.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$squareResult" != "22_3000_5000_DEL_1 breakend_del_unordered_a 22 3000 5000 -2000 DEL T T[22:5000[ 3 " ]
then
        echo "error 1: Error with command_3 (square bad values)"
        exit 1
fi

echo "Command_3 (square) OK"



# command_4: dup 
################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-dup.vcf" -outputFile "./output/command_4_square-bracketed-dup.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-dup.vcf" -outputFile "./output/command_4_angle-bracketed-dup.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-dup.bed" -svtBEDcol 5 -outputFile "./output/command_4_BED-dup.annotated.tsv"

# tail -2 ./input/input-square-bracketed-dup.vcf
# 2       3000    breakend_dup_a  T       ]2:5000]T       .       PASS    SVTYPE=BND;EXTRA=DUP_PAIRED;MATEID=breakend_dup_b       GT      0/0     0/1
# 2       5000    breakend_dup_b  T       T[2:3000[       .       PASS    SVTYPE=BND;EXTRA=DUP_PAIRED;MATEID=breakend_dup_a       GT      0/0     0/1

# tail -1 ./input/input-angle-bracketed-dup.vcf
# 2       3000    breakend_dup_a  T       <DUP>   .       PASS    END=5000EXTRA=DUP_PAIRED;MATEID=breakend_dup_b  GT      0/0     0/1

# tail -1 ./input/input-dup.bed
# 2       2999    5000    breakend_dup_a  DUP

# gr breakend_dup_ Rodrigo_tests_ignore/output.csv
# breakend_dup_a,2,3000,2,5000,T,]2:5000]T,2000,]N,DUP

squareResult=`$cut "./output/command_4_square-bracketed-dup.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$squareResult" != "2_3000_5000_DUP_1 breakend_dup_a 2 3000 5000 2000 DUP T ]2:5000]T 3 " ]
then
        echo "error 1: Error with command_4 (square bad values)"
        exit 1
fi

echo "Command_4 (square) OK"

angleResult=`$cut "./output/command_4_angle-bracketed-dup.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$angleResult" != "2_3000_5000_DUP_1 breakend_dup_a 2 3000 5000 2000 DUP T <DUP> 3 " ]
then
        echo "error 1: Error with command_4 (angle bad values)"
        exit 1
fi

echo "Command_4 (angle) OK"

BEDresult=`$cut "./output/command_4_BED-dup.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$BEDresult" != "2_3000_5000_DUP_1 2 3000 5000 2000 DUP 3 " ]
then
        echo "error 1: Error with command_4 (BED bad values)"
        exit 1
fi

echo "Command_4 (BED) OK"



# command_5: empty_del
#            => the T that is in 53040041 is changed to T[15:53040042[. The T is removed, but added back later.
######################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-empty_del.vcf" -outputFile "./output/command_5_square-bracketed-empty_del.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-empty_del.vcf" -outputFile "./output/command_5_angle-bracketed-empty_del.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-empty_del.bed" -svtBEDcol 5 -outputFile "./output/command_5_BED-empty_del.annotated.tsv"

# tail -2 input-square-bracketed-empty_del.vcf
# 15      53040041        empty_del_a     T       T[15:53040042[  .       PASS    SVTYPE=BND;EXTRA=DEL_NO_SIZE;MATEID=empty_del_b GT      0/0     0/1
# 15      53040042        empty_del_b     A       ]15:53040041]A  .       PASS    SVTYPE=BND;EXTRA=DEL_NO_SIZE;MATEID=empty_del_a GT      0/0     0/1

# tail -1 input-angle-bracketed-empty_del.vcf
# 15      53040041        empty_del_a     T       <DEL>   .       PASS    END=53040042;EXTRA=DEL_NO_SIZE;MATEID=empty_del_b       GT      0/0     0/1

# tail -1 ./input/input-empty_del.bed
# 15      53040040        53040042        empty_del_a     DEL

# gr empty_del_ Rodrigo_tests_ignore/output.csv
# -> empty because SVLEN < 50
# "No SV to annotate in the SVinputFile - Exit without error."


# cat ./output/command_5_square-bracketed-empty_del.unannotated.tsv
# 5_53040041_53040042_DEL_T_T[15:53040042[: variantLength (1) < SVminSize (50) (line 2)
# 15_53040041_53040042_DEL_N_A[15:53040042[: variantLength (1) < SVminSize (50) (line 3)
if [ `grep -c "variantLength (1) < SVminSize (50)" ./output/command_5_square-bracketed-empty_del.unannotated.tsv` != 2 ]
then
        echo "error 1: Error with command_5 (square bad values)"
        exit 1
fi
echo "Command_5 (square) OK"

# cat ./output/command_5_angle-bracketed-empty_del.unannotated.tsv
# 15_53040041_53040042_DEL_T_<DEL>: variantLength (1) < SVminSize (50) (line 2)
if [ `grep -c "variantLength (1) < SVminSize (50)" ./output/command_5_angle-bracketed-empty_del.unannotated.tsv` != 1 ]
then
        echo "error 1: Error with command_5 (angle bad values)"
        exit 1
fi
echo "Command_5 (angle) OK"

# WARNING: No filtering on the SV min size (min > 50 bp) from a BED input file
BEDresult=`$cut "./output/command_5_BED-empty_del.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$BEDresult" != "15_53040041_53040042_DEL_1 15 53040041 53040042 -1 DEL full=3 " ]
then
        echo "error 1: Error with command_5 (BED bad values)"
        exit 1
fi

echo "Command_5 (BED) OK"



# command_6: ins_by_gridss
##########################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-ins_by_gridss.vcf" -SVminSize 1 -outputFile "./output/command_6_square-bracketed-ins_by_gridss.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-ins_by_gridss.vcf" -outputFile "./output/command_6_angle-bracketed-ins_by_gridss.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-ins_by_gridss.bed" -svtBEDcol 5 -outputFile "./output/command_6_BED-ins_by_gridss.annotated.tsv"

# tail -2 ./input/input-square-bracketed-ins_by_gridss.vcf
# 13      53040041        ins_by_gridss   T       TATATATATACACAC[13:53040042[    .       PASS    SVTYPE=BND;EXTRA=DEL_INS_FROM_GRIDSS    GT      0/0     0/1
# 13      53040042        ins_by_gridss   A       ]13:53040041]ATATATATACACACA    .       PASS    SVTYPE=BND;EXTRA=DEL_INS_FROM_GRIDSS    GT      0/0     0/1

# tail -1 ./input/input-angle-bracketed-ins_by_gridss.vcf
# 13      53040041        ins_by_gridss   T       <INS>   .       PASS    EXTRA=DEL_INS_FROM_GRIDSS       GT      0/0     0/1

# tail -1 ./input/input-ins_by_gridss.bed
# 13      53040040        53040041        ns_by_gridss    INS     ATATATATACACAC

# gr ins_by_gridss Rodrigo_tests_ignore/output.csv
# ins_by_gridss,13,53040041,13,53040041,T,TATATATATACACAC,14,,INS


squareResult=`$cut "./output/command_6_square-bracketed-ins_by_gridss.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$squareResult" != "13_53040041_53040041_INS_1 ins_by_gridss 13 53040041 53040041 14 INS T TATATATATACACAC[13:53040042[ NA " ]
then
        echo "error 1: Error with command_6 (square bad values)"
        exit 1
fi

echo "Command_6 (square) OK"

angleResult=`$cut "./output/command_6_angle-bracketed-ins_by_gridss.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_type;REF;ALT;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9}'`
if [ "$angleResult" != "13_53040041_53040041_INS_1 ins_by_gridss 13 53040041 53040041 INS T <INS> NA " ]
then
        echo "error 1: Error with command_6 (angle bad values)"
        exit 1
fi

echo "Command_6 (angle) OK"

BEDresult=`$cut "./output/command_6_BED-ins_by_gridss.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_type;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6}'`
if [ "$BEDresult" != "13_53040041_53040041_INS_1 13 53040041 53040041 INS NA " ]
then
        echo "error 1: Error with command_6 (BED bad values)"
        exit 1
fi

echo "Command_6 (BED) OK"



# command_7: dup
################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-dup.vcf" -outputFile "./output/command_4_square-bracketed-dup.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-dup.vcf" -outputFile "./output/command_4_angle-bracketed-dup.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-dup.bed" -svtBEDcol 5 -outputFile "./output/command_4_BED-dup.annotated.tsv"

# tail -2 ./input/input-square-bracketed-dup.vcf
# 2       3000    breakend_dup_a  T       ]2:5000]T       .       PASS    SVTYPE=BND;EXTRA=DUP_PAIRED;MATEID=breakend_dup_b       GT      0/0     0/1
# 2       5000    breakend_dup_b  T       T[2:3000[       .       PASS    SVTYPE=BND;EXTRA=DUP_PAIRED;MATEID=breakend_dup_a       GT      0/0     0/1

# tail -1 ./input/input-angle-bracketed-dup.vcf
# 2       3000    breakend_dup_a  T       <DUP>   .       PASS    END=5000EXTRA=DUP_PAIRED;MATEID=breakend_dup_b  GT      0/0     0/1

# tail -1 ./input/input-dup.bed
# 2       2999    5000    breakend_dup_a  DUP

# gr breakend_dup_ Rodrigo_tests_ignore/output.csv
# breakend_dup_a,2,3000,2,5000,T,]2:5000]T,2000,]N,DUP

squareResult=`$cut "./output/command_4_square-bracketed-dup.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$squareResult" != "2_3000_5000_DUP_1 breakend_dup_a 2 3000 5000 2000 DUP T ]2:5000]T 3 " ]
then
        echo "error 1: Error with command_7 (square bad values)"
        exit 1
fi

echo "Command_7 (square) OK"

angleResult=`$cut "./output/command_4_angle-bracketed-dup.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$angleResult" != "2_3000_5000_DUP_1 breakend_dup_a 2 3000 5000 2000 DUP T <DUP> 3 " ]
then
        echo "error 1: Error with command_7 (angle bad values)"
        exit 1
fi

echo "Command_7 (angle) OK"

BEDresult=`$cut "./output/command_4_BED-dup.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$BEDresult" != "2_3000_5000_DUP_1 2 3000 5000 2000 DUP 3 " ]
then
        echo "error 1: Error with command_7 (BED bad values)"
        exit 1
fi

echo "Command_7 (BED) OK"



# command_8: inv1 -- T]3:5000]
##############################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-inv1.vcf" -outputFile "./output/command_8_square-bracketed-inv1.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-inv1.vcf" -outputFile "./output/command_8_angle-bracketed-inv1.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-inv1.bed" -svtBEDcol 5 -outputFile "./output/command_8_BED-inv1.annotated.tsv"

# tail -2 ./input/input-square-bracketed-inv1.vcf
# 3       2999    breakend_inv_1_a        T       T]3:5000]       .       PASS    SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_1_b     GT      0/0     0/1
# 3       5000    breakend_inv_1_b        T       [3:2999[T       .       PASS    SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_1_a     GT      0/0     0/1

# tail -1 ./input/input-angle-bracketed-inv1.vcf
# 3       2999    breakend_inv_1_a        T       <INV>   .       PASS    END=5000;EXTRA=INV_PAIRED;MATEID=breakend_inv_1_b       GT      0/0     0/1

# tail -1 ./input/input-inv1.bed
# 3       2998    5000    breakend_inv_1_a        INV

# gr breakend_inv_1 Rodrigo_tests_ignore/output.csv
# breakend_inv_1_a,3,2999,3,5000,T,T]3:5000],2001,N],INV

squareResult=`$cut "./output/command_8_square-bracketed-inv1.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$squareResult" != "3_2999_5000_INV_1 breakend_inv_1_a 3 2999 5000 2001 INV T T]3:5000] NA " ]
then
        echo "error 1: Error with command_8 (square bad values)"
        exit 1
fi

echo "Command_8 (square) OK"

angleResult=`$cut "./output/command_8_angle-bracketed-inv1.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$angleResult" != "3_2999_5000_INV_1 breakend_inv_1_a 3 2999 5000 2001 INV T <INV> NA " ]
then
        echo "error 1: Error with command_8 (angle bad values)"
        exit 1
fi

echo "Command_8 (angle) OK"

BEDresult=`$cut "./output/command_8_BED-inv1.annotated.tsv" "AnnotSV_ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$BEDresult" != "3_2999_5000_INV_1 3 2999 5000 2001 INV NA " ]
then
        echo "error 1: Error with command_8 (BED bad values)"
        exit 1
fi

echo "Command_8 (BED) OK"



# command_9: inv2 -- [3:5001[T
##############################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-inv2.vcf" -outputFile "./output/command_9_square-bracketed-inv2.annotated.tsv" 
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-inv2.vcf" -outputFile "./output/command_9_angle-bracketed-inv2.annotated.tsv"

# tail -2 ./input/input-square-bracketed-inv2.vcf
# 3       3000    breakend_inv_2_a        T       [3:5001[T       .       PASS    SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_2_b     GT      0/0     0/1
# 3       5001    breakend_inv_2_b        T       T]3:3000]       .       PASS    SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_2_a     GT      0/0     0/1

# tail -1 ./input/input-angle-bracketed-inv2.vcf
# 3       3000    breakend_inv_2_a        T       <INV>   .       PASS    END=5001;EXTRA=INV_PAIRED;MATEID=breakend_inv_2_b       GT      0/0     0/1

# gr breakend_inv_2 Rodrigo_tests_ignore/output.csv
# breakend_inv_2_a,3,3000,3,5001,T,[3:5001[T,2001,[N,INV


squareResult=`$cut "./output/command_9_square-bracketed-inv2.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$squareResult" != "3_3000_5001_INV_1 breakend_inv_2_a 3 3000 5001 2001 INV T [3:5001[T NA " ]
then
        echo "error 1: Error with command_9 (square bad values)"
        exit 1
fi

echo "Command_9 (square) OK"

angleResult=`$cut "./output/command_9_angle-bracketed-inv2.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$angleResult" != "3_3000_5001_INV_1 breakend_inv_2_a 3 3000 5001 2001 INV T <INV> NA " ]
then
        echo "error 1: Error with command_9 (angle bad values)"
        exit 1
fi

echo "Command_9 (angle) OK"



# command_10: inv3 (unordered)
#             (we should have a single INV in the output (not 2))
#################################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-inv3.vcf" -outputFile "./output/command_10_square-bracketed-inv3.annotated.tsv"

# tail -2 ./input/input-square-bracketed-inv3.vcf
# 1       9135238 10_2    G       G]1:7852191]    .       PASS    SVTYPE=BND;MATEID=10_1;SVCLASS=inversion        RC:PS   0:0     10:10
# 1       7852191 10_1    A       A]1:9135238]    .       PASS    SVTYPE=BND;MATEID=10_2;SVCLASS=inversion        RC:PS   0:0     10:10

squareResult=`$cut "./output/command_10_square-bracketed-inv3.annotated.tsv" "AnnotSV_ID;ID;SV_chrom;SV_start;SV_end;SV_length;SV_type;REF;ALT;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}'`
if [ "$squareResult" != "1_7852191_9135238_INV_1 10_1 1 7852191 9135238 1283047 INV A A]1:9135238] NA " ]
then
        echo "error 1: Error with command_10 (square bad values)"
        exit 1
fi

echo "Command_10 (square) OK"



# command_11: tra1
#################################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-tra1.vcf" -outputFile "./output/command_11_square-bracketed-tra1.annotated.tsv" 
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-tra1.vcf" -outputFile "./output/command_11_angle-bracketed-tra1.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-tra1.bed" -svtBEDcol 5 -outputFile "./output/command_11_BED-tra1.annotated.tsv"

# tail -2 ./input/input-square-bracketed-tra1.vcf
# 2       321682  tra_a   T       ]13:123456]T    .       PASS    SVTYPE=BND      GT      0/0     0/1
# 13      123456  tra_b   C       C[2:321682[     .       PASS    SVTYPE=BND      GT      0/0     0/1

# tail -2 ./input/input-angle-bracketed-tra1.vcf
# 2       321682  tra_a   T       <TRA>   .       PASS    SVTYPE=BND;CHR2=13;END=123456   GT      0/0     0/1
# 13      123456  tra_b   C       <TRA>   .       PASS    SVTYPE=BND;CHR2=2;END=321682    GT      0/0     0/1

# tail -2 ./input/input-tra1.bed
# 2       321681  321682  tra_a   TRA
# 13      123455  123456  tra_b   TRA

# gr "tra_" Rodrigo_tests_ignore/output.csv
# tra_a,2,321682,13,123456,T,]13:123456]T,0,]]N,TRA
# tra_b,2,321682,13,123456,N,]13:123456]N,0,]]N,TRA



squareResult=`$cut "./output/command_11_square-bracketed-tra1.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;ACMG_class" | grep -v "full=" | tail -2 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$squareResult" != "13_123456_123456_TRA_1 tra_b 0 TRA C C[2:321682[ NA 2_321682_321682_TRA_1 tra_a 0 TRA T ]13:123456]T NA " ]
then
        echo "error 1: Error with command_11 (square bad values)"
        exit 1
fi
echo "Command_11 (square) OK"

angleResult=`$cut "./output/command_11_angle-bracketed-tra1.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;ACMG_class" | grep -v "full=" | tail -2 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$angleResult" != "13_123456_123456_TRA_1 tra_b 0 TRA C <TRA> NA 2_321682_321682_TRA_1 tra_a 0 TRA T <TRA> NA " ]
then
        echo "error 1: Error with command_11 (angle bad values)"
        exit 1
fi
echo "Command_11 (angle) OK"

BEDresult=`$cut "./output/command_11_BED-tra1.annotated.tsv" "AnnotSV_ID;SV_length;SV_type;ACMG_class" | grep -v "full=" | tail -2 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4}'`
if [ "$BEDresult" != "13_123456_123456_TRA_1 0 TRA NA 2_321682_321682_TRA_1 0 TRA NA " ]
then
        echo "error 1: Error with command_11 (BED bad values)"
        exit 1
fi
echo "Command_11 (BED) OK"



# command_12: tra2 (no_mateid)
##############################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-tra2.vcf" -outputFile "./output/command_12_square-bracketed-tra2.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-tra2.vcf" -outputFile "./output/command_12_angle-bracketed-tra2.annotated.tsv"

# tail -2 ./input/input-square-bracketed-tra2.vcf
# 17      198982  trn_no_mateid_a A       A]2:321681]     .       PASS    SVTYPE=BND;EXTRA=TRN_PAIRED_WITHOUT_MATE_ID     GT      0/0     0/1
# 2       321681  trn_no_mateid_b G       G]17:198982]    .       PASS    SVTYPE=BND;EXTRA=TRN_PAIRED_WITHOUT_MATE_ID     GT      0/0     0/1

# tail -2 ./input/input-angle-bracketed-tra2.vcf
# 17      198982  trn_no_mateid_a A       <TRA>   .       PASS    SVTYPE=BND;EXTRA=TRN_PAIRED_WITHOUT_MATE_ID     GT      0/0     0/1
# 2       321681  trn_no_mateid_b G       <TRA>   .       PASS    SVTYPE=BND;EXTRA=TRN_PAIRED_WITHOUT_MATE_ID     GT      0/0     0/1

# gr trn_no_mateid Rodrigo_tests_ignore/output.csv
# trn_no_mateid_b,2,321681,17,198982,G,G]17:198982],0,N]],TRN

squareResult=`$cut "./output/command_12_square-bracketed-tra2.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;ACMG_class" | grep -v "full=" | tail -2 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$squareResult" != "17_198982_198982_TRA_1 trn_no_mateid_a 0 TRA A A]2:321681] NA 2_321681_321681_TRA_1 trn_no_mateid_b 0 TRA G G]17:198982] NA " ]
then
        echo "error 1: Error with command_12 (square bad values)"
        exit 1
fi
echo "Command_12 (square) OK"

angleResult=`$cut "./output/command_12_angle-bracketed-tra2.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;ACMG_class" | grep -v "full=" | tail -2 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7}'`
if [ "$angleResult" != "17_198982_198982_TRA_1 trn_no_mateid_a 0 TRA A <TRA> NA 2_321681_321681_TRA_1 trn_no_mateid_b 0 TRA G <TRA> NA " ]
then
        echo "error 1: Error with command_12 (angle bad values)"
        exit 1
fi
echo "Command_12 (angle) OK"



############################################################################
##
## Square-bracketed notation with only 1 BND, and only the reciprocal BND
##
############################################################################

# The "BNDrescue" flag is added in the INFO field (for the rescue BND) 


# command_13: OnlyTheReciprocalBND_square-bracketed-del_1
#########################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-del_1.vcf" -outputFile "./output/command_13_square-bracketed-del_1.annotated.tsv" -includeci 0

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-del_1.vcf
# 12      5000    breakend_del_1_b        T       ]12:3000]T      .       PASS    SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_1_a;CIPOS=-1,1;CIEND=-5,5       GT      0/0     0/1

squareResult=`$cut "./output/command_13_square-bracketed-del_1.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;INFO;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8}'`
if [ "$squareResult" != "12_3000_5000_DEL_1 breakend_del_1_b -2000 DEL N <DEL> SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_1_a;CIPOS=-5,5;CIEND=-1,1;BNDrescue 3 " ]
then
        echo "error 1: Error with command_13 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_13 (OnlyTheReciprocalBND_square) OK"

if [ `echo $squareResult | grep -c BNDrescue` != 1 ]
then
        echo "error 1: Error with command_13 (OnlyTheReciprocalBND_square: no \"BNDrescue\" flag)"
        exit 1
fi

echo "Command_13 (OnlyTheReciprocalBND_square: \"BNDrescue\" flag present) OK"



# command_14: OnlyTheReciprocalBND_square-bracketed-del_partial_pass
####################################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-del_partial_pass.vcf" -outputFile "./output/command_14_square-bracketed-del_partial_pass.annotated.tsv"

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-del_partial_pass.vcf
# 15      5000    breakend_del_partial_pass_a     T       ]15:3000]T      .       FAIL    SVTYPE=BND;EXTRA=DEL_PAIRED_PARTIAL_PASS;MATEID=breakend_del_partial_pass_b     GT      0/0     0/1

squareResult=`$cut "./output/command_14_square-bracketed-del_partial_pass.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;INFO;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8}'`
if [ "$squareResult" != "15_3000_5000_DEL_1 breakend_del_partial_pass_a -2000 DEL N <DEL> SVTYPE=BND;EXTRA=DEL_PAIRED_PARTIAL_PASS;MATEID=breakend_del_partial_pass_b;BNDrescue 3 " ]
then
        echo "error 1: Error with command_14 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_14 (OnlyTheReciprocalBND_square) OK"



# command_15: OnlyTheReciprocalBND_square-bracketed-del_unordered
#################################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-del_unordered.vcf" -outputFile "./output/command_15_square-bracketed-del_unordered.annotated.tsv"

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-del_unordered.vcf
# 22      5000    breakend_del_unordered_b        T       ]22:3000]T      .       PASS    SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_unordered_a     GT      0/0     0/1

squareResult=`$cut "./output/command_15_square-bracketed-del_unordered.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;INFO;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8}'`
if [ "$squareResult" != "22_3000_5000_DEL_1 breakend_del_unordered_b -2000 DEL N <DEL> SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_unordered_a;BNDrescue 3 " ]
then
        echo "error 1: Error with command_15 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_15 (OnlyTheReciprocalBND_square) OK"



# command_16: OnlyTheReciprocalBND_square-bracketed-dup 
#######################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-dup.vcf" -outputFile "./output/command_16_square-bracketed-dup.annotated.tsv"

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-dup.vcf
# 2       5000    breakend_dup_b  T       T[2:3000[       .       PASS    SVTYPE=BND;EXTRA=DUP_PAIRED;MATEID=breakend_dup_a       GT      0/0     0/1

squareResult=`$cut "./output/command_16_square-bracketed-dup.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;INFO;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8}'`
if [ "$squareResult" != "2_3000_5000_DUP_1 breakend_dup_b 2000 DUP N <DUP> SVTYPE=BND;EXTRA=DUP_PAIRED;MATEID=breakend_dup_a;BNDrescue 3 " ]
then
        echo "error 1: Error with command_16 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_16 (OnlyTheReciprocalBND_square) OK"



# command_17: OnlyTheReciprocalBND_square-bracketed-empty_del
#############################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-empty_del.vcf" -outputFile "./output/command_17_square-bracketed-empty_del.annotated.tsv"

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-empty_del.vcf
# 15      53040042        empty_del_b     A       ]15:53040041]A  .       PASS    SVTYPE=BND;EXTRA=DEL_NO_SIZE;MATEID=empty_del_a GT      0/0     0/1

if [ `grep -c "variantLength (1) < SVminSize (50)" ./output/command_17_square-bracketed-empty_del.unannotated.tsv` != 1 ]
then
        echo "error 1: Error with command_17 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_17 (OnlyTheReciprocalBND_square) OK"



# command_18: OnlyTheReciprocalBND_square-bracketed-ins_by_gridss 
#################################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-ins_by_gridss.vcf" -SVminSize 1 -outputFile "./output/command_18_square-bracketed-ins_by_gridss.annotated.tsv"

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-ins_by_gridss.vcf
# 13      53040042        ins_by_gridss   A       ]13:53040041]ATATATATACACACA    .       PASS    SVTYPE=BND;EXTRA=DEL_INS_FROM_GRIDSS    GT      0/0     0/1

squareResult=`$cut "./output/command_18_square-bracketed-ins_by_gridss.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;INFO;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8}'`
if [ "$squareResult" != "13_53040041_53040041_INS_1 ins_by_gridss 14 INS N <INS> SVTYPE=BND;EXTRA=DEL_INS_FROM_GRIDSS;BNDrescue NA " ]
then
        echo "error 1: Error with command_18 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_18 (OnlyTheReciprocalBND_square) OK"



# command_19: OnlyTheReciprocalBND_square-bracketed-inv1
########################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-inv1.vcf" -outputFile "./output/command_19_square-bracketed-inv1.annotated.tsv"

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-inv1.vcf
# 3       5000    breakend_inv_1_b        T       [3:2999[T       .       PASS    SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_1_a     GT    0/0      0/1

squareResult=`$cut "./output/command_19_square-bracketed-inv1.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;INFO;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8}'`
if [ "$squareResult" != "3_2999_5000_INV_1 breakend_inv_1_b 2001 INV N <INV> SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_1_a;BNDrescue NA " ]
then
        echo "error 1: Error with command_19 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_19 (OnlyTheReciprocalBND_square) OK"



# command_20: OnlyTheReciprocalBND_square-bracketed-inv2
########################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-inv2.vcf" -outputFile "./output/command_20_square-bracketed-inv2.annotated.tsv"

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-inv2.vcf
# 3       5001    breakend_inv_2_b        T       T]3:3000]       .       PASS    SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_2_a     GT    0/0      0/1

squareResult=`$cut "./output/command_20_square-bracketed-inv2.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;INFO;ACMG_class" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8}'`
if [ "$squareResult" != "3_3000_5001_INV_1 breakend_inv_2_b 2001 INV N <INV> SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_2_a;BNDrescue NA " ]
then
        echo "error 1: Error with command_20 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_20 (OnlyTheReciprocalBND_square) OK"



# command_21: OnlyTheReciprocalBND_square-bracketed-tra1_1
##########################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-tra1_1.vcf" -outputFile "./output/input_OnlyTheReciprocalBND_square-bracketed-tra1_1.annotated.tsv"

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-tra1_1.vcf
# 2       321682  tra_a   T       ]13:123456]T    .       PASS    SVTYPE=BND      GT      0/0     0/1

squareResult=`$cut "./output/input_OnlyTheReciprocalBND_square-bracketed-tra1_1.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;INFO;ACMG_class" | grep -v "full=" | tail -2 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8}'`
if [ "$squareResult" != "13_123456_123456_TRA_1 tra_a 0 TRA N <TRA> SVTYPE=BND;BNDrescue NA 2_321682_321682_TRA_1 tra_a 0 TRA T ]13:123456]T SVTYPE=BND NA " ]
then
        echo "error 1: Error with command_21 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_21 (OnlyTheReciprocalBND_square) OK"



# command_21: OnlyTheReciprocalBND_square-bracketed-tra1_2
##########################################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input_OnlyTheReciprocalBND_square-bracketed-tra1_2.vcf" -outputFile "./output/input_OnlyTheReciprocalBND_square-bracketed-tra1_2.annotated.tsv"

# tail -1 ./input/input_OnlyTheReciprocalBND_square-bracketed-tra1_2.vcf
# 13      123456  tra_b   C       C[2:321682[     .       PASS    SVTYPE=BND      GT      0/0     0/1

squareResult=`$cut "./output/input_OnlyTheReciprocalBND_square-bracketed-tra1_2.annotated.tsv" "AnnotSV_ID;ID;SV_length;SV_type;REF;ALT;INFO;ACMG_class" | grep -v "full=" | tail -2 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1, $2, $3, $4, $5, $6, $7, $8}'`
if [ "$squareResult" != "13_123456_123456_TRA_1 tra_b 0 TRA C C[2:321682[ SVTYPE=BND NA 2_321682_321682_TRA_1 tra_b 0 TRA N <TRA> SVTYPE=BND;BNDrescue NA " ]
then
        echo "error 1: Error with command_21 (OnlyTheReciprocalBND_square bad values)"
        exit 1
fi
echo "Command_21 (OnlyTheReciprocalBND_square) OK"




################################################################
################################################################
################################################################

echo "ok - Finished"


