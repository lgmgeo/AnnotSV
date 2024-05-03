#!/bin/bash -x

set -eo pipefail

cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"


# -reciprocal:   Use of a reciprocal overlap between SV and user features (only for annotations with features overlapping the SV)
#                Values: 0 (default) or 1
#                <=> Only for custom SVincludedInFt annotations
#
# NOTE: It is to notice that, for pathogenic SV annotations, a reciprocal overlap cannot be used. The “-reciprocal” option can only be used with custom annotations with features overlapping the SV.
# NOTE: It is to notice that, for benign SV annotations, a reciprocal overlap cannot be used. The “-reciprocal” option can only be used with custom annotations with features overlapping the SV.
#
# Custom annotations:
# - By placing the BED file in the “SVincludedInFt” directory, only the features overlapping 100% of the SV will be reported.
# In this case, a reciprocal overlap can be used (see "reciprocal" option in USAGE/OPTIONS).


cp "./input/user_SVincludedInFt.bed" $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh37/SVincludedInFt/


# Checks of the "-reciprocal 0" and "-overlap 80" options with a BED input file:
################################################################################

rm -f "./output/testBED_tcl.annotated.reciprocalNo.Overlap80.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/testBED_tcl.annotated.reciprocalNo.Overlap80.tsv"  -svtBEDcol 4 -reciprocal 0 -overlap 80 -genomeBuild GRCh37
# => "SVincludedInFt" feature reported in the "BBBBBBBBBBBBBBBBBBBB" column

annotations=`$cut "./output/testBED_tcl.annotated.reciprocalNo.Overlap80.tsv" "BBBBBBBBBBBBBBBBBBBB" | sort -u`
annotations=`echo $annotations | sed "s/ /_/"`
if [ $annotations == "BBBBBBBBBBBBBBBBBBBB_SVincludedInFt" ]
then
        echo "Ok"
else
        echo "error 1: the \"-reciprocal\" option does not work"
        exit 1
fi



# Checks of the "-reciprocal 1" and "-overlap 80" options with a BED input file:
################################################################################

rm -f "./output/testBED_tcl.annotated.reciprocalYes.Overlap80.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/testBED_tcl.annotated.reciprocalYes.Overlap80.tsv" -svtBEDcol 4 -reciprocal 1 -overlap 80 -genomeBuild GRCh37
# => "SVincludedInFt" feature not reported in the "BBBBBBBBBBBBBBBBBBBB" column

annotations=`$cut "./output/testBED_tcl.annotated.reciprocalYes.Overlap80.tsv" "BBBBBBBBBBBBBBBBBBBB" | sort -u`

if [ $annotations == "BBBBBBBBBBBBBBBBBBBB" ]
then
        echo "Ok"
else
        echo "error 1: the \"-reciprocal\" option does not work"
        exit 1
fi



# Checks of the "-reciprocal 0" and "-overlap 80" options with a VCF input file:
################################################################################

rm -f "./output/testVCF_tcl.annotated.reciprocalNo.Overlap80.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.vcf" -SVinputInfo 1 -outputFile "./output/testVCF_tcl.annotated.reciprocalNo.Overlap80.tsv" -reciprocal 0 -overlap 80 -genomeBuild GRCh37
# => "SVincludedInFt" feature reported in the "BBBBBBBBBBBBBBBBBBBB" column

annotations=`$cut "./output/testVCF_tcl.annotated.reciprocalNo.Overlap80.tsv" "BBBBBBBBBBBBBBBBBBBB" | sort -u`
annotations=`echo $annotations | sed "s/ /_/"`

if [ $annotations == "BBBBBBBBBBBBBBBBBBBB_SVincludedInFt" ]
then
        echo "Ok"
else
        echo "error 1: the \"-reciprocal\" option does not work"
        exit 1
fi



# Checks of the "-reciprocal 1" and "-overlap 80" options with a VCF input file:
################################################################################

rm -f "./output/testVCF_tcl.annotated.reciprocalYes.Overlap80.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.vcf" -SVinputInfo 1 -outputFile "./output/testVCF_tcl.annotated.reciprocalYes.Overlap80.tsv" -reciprocal 1 -overlap 80 -genomeBuild GRCh37
# => "SVincludedInFt" feature not reported in the "BBBBBBBBBBBBBBBBBBBB" column

annotations=`$cut "./output/testVCF_tcl.annotated.reciprocalYes.Overlap80.tsv" "BBBBBBBBBBBBBBBBBBBB" | sort -u`

if [ $annotations == "BBBBBBBBBBBBBBBBBBBB" ]
then
        echo "Ok"
else
        echo "error 1: the \"-reciprocal\" option does not work"
        exit 1
fi



## Nettoyage :
##############
rm $ANNOTSV/share/AnnotSV/Annotations_Human/Users/GRCh37/SVincludedInFt/user_SVincludedInFt.*




echo "ok - Finished"

