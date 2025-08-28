#!/bin/bash -x

set -eo pipefail

###################################################################
####                                                           ####
####                Check the SV_length in the output          ####
####                                                           ####
###################################################################

mkdir -p ./output
rm -f ./output/category_*.annotated.tsv
rm -f ./output/category_*.unannotated.tsv


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"


##############################################################
# Category 1 -
# VCF: <INS> (symbolic SV allele <=> angle-bracketed notation)
##############################################################

# ...Ajouter un exemple



##############################################################
# Category 2 -
# VCF: <DEL>, <INV>, <DUP>, <CN0>... (symbolic SV allele, other than INS <=> angle-bracketed notation)
######################################################################################################


# input-angle-bracketed-del.vcf 
###############################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-del.vcf" -outputFile "./output/category_2_angle-bracketed-del.annotated.tsv"

# tail -1 ./input/input-angle-bracketed-del.vcf
# 12      3000    breakend_del_1_a        T       <DEL>   .       PASS    END=5000;EXTRA=DEL_PAIRED;MATEID=breakend_del_1_b       GT      0/0     0/1

angleResult=`$cut "./output/category_2_angle-bracketed-del.annotated.tsv" "SV_length" | tail -1  | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$angleResult" != "-2000 " ]
then
        echo "Error with category_2 (angle del)"
        exit 1
fi

echo "Category_2 (angle del) OK"


# input-angle-bracketed-dup.vcf
###############################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-dup.vcf" -outputFile "./output/category_2_angle-bracketed-dup.annotated.tsv"

# tail -1 ./input/input-angle-bracketed-dup.vcf
# 2       3000    breakend_dup_a  T       <DUP>   .       PASS    END=5000EXTRA=DUP_PAIRED;MATEID=breakend_dup_b  GT      0/0     0/1

angleResult=`$cut "./output/category_2_angle-bracketed-dup.annotated.tsv" "SV_length" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$angleResult" != "2000 " ]
then
        echo "Error with category_2 (angle dup)"
        exit 1
fi

echo "Category_2 (angle dup) OK"


# input-angle-bracketed-inv.vcf
###############################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-inv.vcf" -outputFile "./output/category_2_angle-bracketed-inv.annotated.tsv"

# tail -1 ./input/input-angle-bracketed-inv.vcf
# 3       2999    breakend_inv_1_a        T       <INV>   .       PASS    END=5000;EXTRA=INV_PAIRED;MATEID=breakend_inv_1_b       GT      0/0     0/1

angleResult=`$cut "./output/category_2_angle-bracketed-inv.annotated.tsv" "SV_length" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$angleResult" != "2001 " ]
then
        echo "Error with category_2 (angle inv)"
        exit 1
fi

echo "Category_2 (angle inv) OK"


# input-angle-bracketed-tra.vcf
###############################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-angle-bracketed-inv.vcf" -outputFile "./output/category_2_angle-bracketed-inv.annotated.tsv"

# tail -1 ./input/input-angle-bracketed-inv.vcf
# 3       2999    breakend_inv_1_a        T       <INV>   .       PASS    END=5000;EXTRA=INV_PAIRED;MATEID=breakend_inv_1_b       GT      0/0     0/1

angleResult=`$cut "./output/category_2_angle-bracketed-inv.annotated.tsv" "SV_length" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$angleResult" != "2001 " ]
then
        echo "Error with category_2 (angle tra)"
        exit 1
fi

echo "Category_2 (angle tra) OK"




################################
# Category 3 -
# VCF: square-bracketed notation
################################


# input-square-bracketed-del.vcf
################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-del.vcf" -outputFile "./output/category_3_square-bracketed-del.annotated.tsv"

# tail -2 ./input/input-square-bracketed-del.vcf
# 12      3000    breakend_del_1_a        T       T[12:5000[      .       PASS    SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_1_b     GT      0/0     0/1
# 12      5000    breakend_del_1_b        T       ]12:3000]T      .       PASS    SVTYPE=BND;EXTRA=DEL_PAIRED;MATEID=breakend_del_1_a     GT      0/0     0/1

squareResult=`$cut "./output/category_3_square-bracketed-del.annotated.tsv" "SV_length" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$squareResult" != "-2000 " ]
then
        echo "Error with category_3 (square del)"
        exit 1
fi

echo "Category_3 (square del) OK"


# input-square-bracketed-dup.vcf
################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-dup.vcf" -outputFile "./output/category_3_square-bracketed-dup.annotated.tsv"

# tail -2 ./input/input-square-bracketed-dup.vcf
# 2       3000    breakend_dup_a  T       ]2:5000]T       .       PASS    SVTYPE=BND;EXTRA=DUP_PAIRED;MATEID=breakend_dup_b       GT      0/0     0/1
# 2       5000    breakend_dup_b  T       T[2:3000[       .       PASS    SVTYPE=BND;EXTRA=DUP_PAIRED;MATEID=breakend_dup_a       GT      0/0     0/1

squareResult=`$cut "./output/category_3_square-bracketed-dup.annotated.tsv" "SV_length" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$squareResult" != "2000 " ]
then
        echo "Error with category_3 (square dup)"
        exit 1
fi

echo "Category_3 (square dup) OK"


# input-square-bracketed-ins_by_gridss.vcf
##########################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-ins_by_gridss.vcf" -SVminSize 1 -outputFile "./output/category_3_square-bracketed-ins_by_gridss.annotated.tsv"

# tail -2 ./input/input-square-bracketed-ins_by_gridss.vcf
# 13      53040041        ins_by_gridss   T       TATATATATACACAC[13:53040042[    .       PASS    SVTYPE=BND;EXTRA=DEL_INS_FROM_GRIDSS    GT      0/0     0/1
# 13      53040042        ins_by_gridss   A       ]13:53040041]ATATATATACACACA    .       PASS    SVTYPE=BND;EXTRA=DEL_INS_FROM_GRIDSS    GT      0/0     0/1

squareResult=`$cut "./output/category_3_square-bracketed-ins_by_gridss.annotated.tsv" "SV_length" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }
'`
if [ "$squareResult" != "14 " ]
then
        echo "Error with category_3 (square ins)"
        exit 1
fi

echo "Category_3 (square ins) OK"


# input-square-bracketed-inv.vcf
################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-inv.vcf" -outputFile "./output/category_3_square-bracketed-inv.annotated.tsv"

# tail -2 ./input/input-square-bracketed-inv.vcf
# 3       2999    breakend_inv_1_a        T       T]3:5000]       .       PASS    SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_1_b     GT      0/0     0/1
# 3       5000    breakend_inv_1_b        T       [3:2999[T       .       PASS    SVTYPE=BND;EXTRA=INV_PAIRED;MATEID=breakend_inv_1_a     GT      0/0     0/1
squareResult=`$cut "./output/category_3_square-bracketed-inv.annotated.tsv" "SV_length" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$squareResult" != "2001 " ]
then
        echo "Error with category_3 (square inv)"
        exit 1
fi

echo "Category_3 (square inv) OK"


# input-square-bracketed-tra.vcf
################################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-square-bracketed-tra.vcf" -outputFile "./output/category_3_square-bracketed-tra.annotated.tsv"

# tail -2 ./input/input-square-bracketed-tra.vcf
# 2       321682  tra_a   T       ]13:123456]T    .       PASS    SVTYPE=BND      GT      0/0     0/1
# 13      123456  tra_b   C       C[2:321682[     .       PASS    SVTYPE=BND      GT      0/0     0/1

squareResult=`$cut "./output/category_3_square-bracketed-tra.annotated.tsv" "SV_length" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$squareResult" != "0 " ]
then
        echo "Error with category_3 (square tra)"
        exit 1
fi
echo "Category_3 (square tra) OK"


########################
# Category 4 -
# VCF: REF=A ALT=ATTCGTT
########################

# ...Ajouter un exemple




###########################################
# Category 5 -
# BED: CHROM POS END ... (0-based notation)
###########################################


# input-del.bed
###############

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-del.bed" -svtBEDcol 5 -outputFile "./output/category_5_BED-del.annotated.tsv"

# tail -1 ./input/input-del.bed
# 12      2999    5000    breakend_del_1_a        DEL

# gr breakend_del_ Rodrigo_tests_ignore/output.csv
# breakend_del_1_a,12,3000,12,5000,T,T[12:5000[,2000,N[,DEL


BEDresult=`$cut "./output/category_5_BED-del.annotated.tsv" "SV_length" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1}'`
if [ "$BEDresult" != "-2000 " ]
then
        echo "Error with category_5 (BED del)"
        exit 1
fi

echo "Category_5 (BED del) OK"


# input-dup.bed
###############

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-dup.bed" -svtBEDcol 5 -outputFile "./output/category_5_BED-dup.annotated.tsv"

# tail -1 ./input/input-dup.bed
# 2       2999    5000    breakend_dup_a  DUP

# gr breakend_dup_ Rodrigo_tests_ignore/output.csv
# breakend_dup_a,2,3000,2,5000,T,]2:5000]T,2000,]N,DUP

BEDresult=`$cut "./output/category_5_BED-dup.annotated.tsv" "SV_length" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$BEDresult" != "2000 " ]
then
        echo "Error with category_5 (BED dup)"
        exit 1
fi

echo "Category_5 (BED dup) OK"


# input-ins_by_gridss.bed
#########################

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-ins_by_gridss.bed" -svtBEDcol 5 -outputFile "./output/category_5_BED-ins_by_gridss.annotated.tsv"


# tail -1 ./input/input-ins_by_gridss.bed
# 13      53040040        53040041        ns_by_gridss    INS     ATATATATACACAC

# gr ins_by_gridss Rodrigo_tests_ignore/output.csv
# ins_by_gridss,13,53040041,13,53040041,T,TATATATATACACAC,14,,INS


BEDresult=`$cut "./output/category_5_BED-ins_by_gridss.annotated.tsv" "SV_length" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$BEDresult" != " " ]
then
        echo "Error with category_5 (BED ins)"
        exit 1
fi

echo "Category_5 (BED ins) OK"


# input-inv.bed
###############

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-inv.bed" -svtBEDcol 5 -outputFile "./output/category_5_BED-inv.annotated.tsv"

# tail -1 ./input/input-inv.bed
# 3       2998    5000    breakend_inv_1_a        INV

BEDresult=`$cut "./output/category_5_BED-inv.annotated.tsv" "SV_length" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$BEDresult" != "2001 " ]
then
        echo "Error with category_5 (BED inv)"
        exit 1
fi

echo "Category_5 (BED inv) OK"


# input-tra.bed
###############

$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input-tra.bed" -svtBEDcol 5 -outputFile "./output/category_5_BED-tra.annotated.tsv"

# tail -2 ./input/input-tra.bed
# 2       321681  321682  tra_a   TRA
# 13      123455  123456  tra_b   TRA

BEDresult=`$cut "./output/category_5_BED-tra.annotated.tsv" "SV_length" | grep -v "full=" | tail -1 | awk 'BEGIN{ORS=" "; OFS=" "} { print $1 }'`
if [ "$BEDresult" != "0 " ]
then
        echo "Error with category_5 (BED tra)"
        exit 1
fi
echo "Category_5 (BED tra) OK"

################################################################
################################################################
################################################################

echo "ok - Finished"


