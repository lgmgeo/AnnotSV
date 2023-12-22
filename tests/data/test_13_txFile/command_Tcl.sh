#!/bin/bash -x

set -eo pipefail

cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"


################################
# Check the "-txFile" option
################################

# INFO: Transcript annotation comes from $ANNOTSV/share/AnnotSV/Annotations_Human/Genes/GRCh37/genes.RefSeq.sorted.bed


## Check 1: the CDK11B gene (with a BED input file)
###################################################

# By default, the CDK11B gene is annotated with the "NM_001787" transcript (and not the "NM_033486")
rm -f "./output/test.check1_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.check1.bed" -outputFile "./output/test.check1_tcl.annotated.tsv" -svtBEDcol 5 -genomeBuild GRCh37
annotations=`$cut "./output/test.check1_tcl.annotated.tsv" "Gene_name;Tx" | grep "CDK11B"`
if [ `echo $annotations | grep -c "NM_001787"` != 0 ]
then
        echo "Ok"
else
        echo "error 1: the CDK11B gene is not annotated with the NM_001787 transcript"
        exit 1
fi


# Changing of the transcript default:
# => "NM_033486" transcript is given in the "txFile.check1.txt" file
rm -f "./output/test.check1.WithTx_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.check1.bed" -outputFile "./output/test.check1.WithTx_tcl.annotated.tsv" -svtBEDcol 5 -txFile "./input/txFile.check1.txt" -genomeBuild GRCh37
annotations=`$cut "./output/test.check1.WithTx_tcl.annotated.tsv" "Gene_name;Tx" | grep "CDK11B"`
#=> Now, the CDK11B gene is annotated with "NM_033486"
if [ `echo $annotations | grep -c "NM_033486"` != 0 ]
then
        echo "Ok"
else
        echo "error 1: the CDK11B gene is not annotated with the NM_033486 transcript"
        exit 1
fi


## Check 2: the NAV1 gene (with a VCF input file)
#################################################

# By default, the NAV1 gene is annotated with the "NM_020443" transcript (and not the "NM_001167738")
rm -f "./output/test.check2_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.check2.vcf" -outputFile "./output/test.check2_tcl.annotated.tsv" -genomeBuild GRCh37
annotations=`$cut "./output/test.check2_tcl.annotated.tsv" "Gene_name;Tx" | grep "NAV1"`
if [ `echo $annotations | grep -c "NM_020443"` != 0 ]
then
        echo "Ok"
else
        echo "error 1: the NAV1 gene is not annotated with the NM_020443 transcript"
        exit 1
fi


# Changing of the transcript default:
# => "NM_001167738" transcript is given in the "txfile.check2.txt" file
rm -f "./output/test.check2.WithTx_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.check2.vcf" -outputFile "./output/test.check2.WithTx_tcl.annotated.tsv" -txFile "./input/txFile.check2.txt" -genomeBuild GRCh37
annotations=`$cut "./output/test.check2.WithTx_tcl.annotated.tsv" "Gene_name;Tx" | grep "NAV1"`
#=> Now, the NAV1 gene is annotated with "NM_001167738"
if [ `echo $annotations | grep -c "NM_001167738"` != 0 ]
then
        echo "Ok"
else
        echo "error 1: the NAV1 gene is not annotated with the NM_001167738 transcript"
        exit 1
fi



############################
# Check the "-tx" option
############################

## Check 3
##########

# By default, tx="RefSeq"
# Check with tx="ENSEMBL"
rm -f "./output/test.check3_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.check1.bed" -outputFile "./output/test.check3_tcl.annotated.tsv" -svtBEDcol 5 -tx "ENSEMBL" -genomeBuild GRCh38
annotations=`$cut "./output/test.check3_tcl.annotated.tsv" "Gene_name;Tx" | grep "CDK11B"`
if [ `echo $annotations | grep -c "ENST00000407249"` != 0 ]
then
        echo "Ok"
else
        echo "error 1: the CDK11B gene is not annotated with the ENST00000407249 transcript"
        exit 1
fi


# Changing of the transcript default:
# => "ENST00000341832" transcript is given in the "txfile.check3.txt" file
rm -f "./output/test.check3.WithTx_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.check1.bed" -outputFile "./output/test.check3.WithTx_tcl.annotated.tsv" -txFile "./input/txFile.check3.txt" -svtBEDcol 5 -tx "ENSEMBL" -genomeBuild GRCh38
annotations=`$cut "./output/test.check3.WithTx_tcl.annotated.tsv" "Gene_name;Tx" | grep "CDK11B"`
#=> Now, the CDK11B gene is annotated with "ENST00000341832"
if [ `echo $annotations | grep -c "ENST00000341832"` != 0 ]
then
        echo "Ok"
else
        echo "error 1: the CDK11B gene is not annotated with the ENST00000341832 transcript"
        exit 1
fi




echo "ok - Finished"







