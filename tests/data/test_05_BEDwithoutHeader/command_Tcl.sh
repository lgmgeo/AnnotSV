#!/bin/bash -x

set -eo pipefail


# Check input BED file with header
##################################
# more input/test.bed
# ---------------------------------------------
# #chrom  start   end     tutu    SV type
# 1       2902537 3529849 ii      DEL
# ---------------------------------------------
# => We should have the "tutu" annotation column (extracted from the BED input file)

# Check input BED file without header
#####################################
# 
# more input/test.withoutHeader.bed
# ---------------------------------------------
# 1       2902537 3529849 ii      DEL
# ---------------------------------------------

# output colum names NOT given in the input file:
# => only "SV_type" should be in the output
# => we should not have "tutu" anymore, but "user#1"
# ---------------------------------------------
# user#1
# ii
# ii
# ii
# ii
# ii
# ii
# ii
# ii
# ---------------------------------------------



# Input BED file with header
############################
rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtBEDcol 5 -overlap 70 -genomeBuild GRCh37

# Check if the "tutu" annotation column from the BED file was added
if [ `grep -c "tutu" "./output/test_tcl.annotated.tsv"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: \"tutu\" annotation column from the input BED file NOT added"
        exit 1
fi


# Input BED file without header
###############################
rm -f "./output/test.withoutHeader_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.withoutHeader.bed" -SVinputInfo 1 -outputFile "./output/test.withoutHeader_tcl.annotated.tsv" -svtBEDcol 5 -overlap 70 -genomeBuild GRCh37 

# Check if the "SV_type" annotation column from the BED file was added
if [ `grep -c "SV_type" "./output/test_tcl.annotated.tsv"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: \"SV_type\" annotation column from the input BED file NOT added"
        exit 1
fi


# Check if the "user#1" annotation column from the BED file was added
if [ `grep -c "user#1" "./output/test.withoutHeader_tcl.annotated.tsv"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: \"user#1\" annotation column from the input BED file NOT added"
        exit 1
fi



echo "ok - Finished"



