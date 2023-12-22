#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"


# -candidateGenesFile:            Path of a file containing the candidate genes of the user (gene names can be space-separated, tabulation-separated, or line-break-separated)
# -candidateGenesFiltering:       To select only the SV annotations ("split" and "full") overlapping a gene from the "candidateGenesFile"

# candidateGenesFile.txt:
# BBS5 CDK11A SSU72


# Run without the "-candidateGenesFiltering" option
###################################################

rm -f "./output/test.without_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test.without_tcl.annotated.tsv" -svtBEDcol 5 -candidateGenesFile "./input/candidateGenesFile.txt" -genomeBuild GRCh37

nLines_Without=`wc -l "./output/test.without_tcl.annotated.tsv" | cut -d" " -f 1`


# Run with the "-candidateGenesFiltering" option
################################################

rm -f "./output/test.with_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test.with_tcl.annotated.tsv" -svtBEDcol 5 -candidateGenesFile "./input/candidateGenesFile.txt" -candidateGenesFiltering 1 -genomeBuild GRCh37 

nLines_With=`wc -l "./output/test.with_tcl.annotated.tsv" | cut -d" " -f 1`

nBadLines_With=`$cut "./output/test.with_tcl.annotated.tsv" "AnnotSV_ID;Gene_name;Annotation_mode" | egrep -v "BBS5|CDK11A|SSU72" | wc -l`
# => Only the header line
# => $nBadLines_With == 1


# Check the line filtering
##########################
if [ $nLines_Without > $nLines_With ]
then
        echo "Ok"
else
        echo "error 1: the \"-candidateGenesFiltering\" option does not work"
        exit 1
fi

if [ $nBadLines_With == 1 ]
then
        echo "Ok"
else
        echo "error 1: the \"-candidateGenesFiltering\" option does not work"
        exit 1
fi


echo "ok - Finished"

