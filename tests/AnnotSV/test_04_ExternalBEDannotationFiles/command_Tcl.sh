#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"

# Check the use of an external BED file (provided by the user)
##############################################################
#
# User_BED_file_1 = ./input/user_FtIncludedInSV.bed
# User_BED_file_2 = ./input/user_SVincludedInFt.bed


User_BED_file_1="./input/user_FtIncludedInSV.bed"
User_BED_file_2="./input/user_SVincludedInFt.bed"


cp "$User_BED_file_1" $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh37/FtIncludedInSV/
cp "$User_BED_file_2" $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh37/SVincludedInFt

rm -f "./output/test_tcl.annotated.tsv" "./output/output.log"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtBEDcol 5 -overlap 70 -genomeBuild GRCh37 &> "./output/output.log"


# Validation of the use of "User_BED_file_1"
############################################
if [ `grep -c "AAAAAAAAAAAAAAAAAAAA" "./output/output.log"` == 2 ] && [ `grep -c "FtIncludedInSV" "./output/output.log"` == 1 ]
then
	echo "Ok"
else
        echo "error 1: $User_BED_file_1 not used!"
        exit 1
fi

# Validation of the use of "User_BED_file_2"
############################################
if [ `grep -c "BBBBBBBBBBBBBBBBBBBB" "./output/output.log"` == 2 ] && [ `grep -c "SVincludedInFt" "./output/output.log"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: $User_BED_file_2 not used!"
        exit 1
fi


# The output columns should have "SVincludedInFt" and "FtIncludedInSV" annotations
##################################################################################
annotations=`$cut "./output/test_tcl.annotated.tsv" "BBBBBBBBBBBBBBBBBBBB"`
if [ `echo $annotations | tr " " "\n" | grep -c "SVincludedInFt"` == 8 ]
then
        echo "Ok"
else
        echo "error 1: SVincludedInFt annotations NOT present in the output file "
        exit 1
fi

annotations=`$cut "./output/test_tcl.annotated.tsv" "AAAAAAAAAAAAAAAAAAAA"`
if [ `echo $annotations | tr " " "\n" | grep -c "FtIncludedInSV"` == 2 ]
then
        echo "Ok"
else
        echo "error 1: FtIncludedInSV annotations NOT present in the output file "
        exit 1
fi


# Clean
#######
rm $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh37/FtIncludedInSV/user_FtIncludedInSV.*
rm $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh37/SVincludedInFt/user_SVincludedInFt.*



echo "ok - Finished"

