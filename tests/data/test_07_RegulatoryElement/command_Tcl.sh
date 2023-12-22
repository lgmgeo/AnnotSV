#!/bin/bash -x

set -eo pipefail

cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"

# Case study 1:
###############

# Fœtus with retrognathia and cleft palate ( HP:0000278; HP:0000175 )
# A 3.2 Mb deletion encompassing known pathogenic ClinGen and OMIM CNV close to SOX9:
# GRCh37, chr17:66426540-69629770 DEL
# - SOX9 gene not overlapped with the SV
# - SOX9 RE overlapped with the SV


# Case study 2:
###############

# Deletion at about 600 kB of MEF2C
# MEF2C (NM_001131005), GRCh38, chr5:88914328-91802402
# - MEF2C not overlapped with the SV
# - MEF2C RE overlapped with the SV


## Check with the SOX9 gene, Exomiser and GRCh37
################################################

rm -f "./output/delSOX9_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/delSOX9.bed" -hpo "HP:0000278;HP:0000175" -outputFile "./output/delSOX9_tcl.annotated.tsv" -svtBEDcol 4 -REreport 1 -genomeBuild GRCh37

# Check that the "SOX9 RE" is annotated
val=`$cut "./output/delSOX9_tcl.annotated.tsv" "RE_gene"`
if [ `echo $val | egrep -c ";SOX9 \(EX=0.9.{3}/morbid/RE=GH_Enhancer\);"` != 1 ]
then
        echo "error1: The SOX9 RE is not annotated!"
        exit 1
fi

# Check that the "ACMG_class" is "5"
val=`$cut "./output/delSOX9_tcl.annotated.tsv" "ACMG_class"`
if [ `echo $val | grep -c "5"` != 1 ]
then
        echo "error1: The ACMG_class is not set to 5"
        exit 1
fi

# Check if the RE report exists
if [ ! -e "./output/delSOX9_tcl.SV_RE_intersect.report" ] 
then
	echo "error1: The RE report (./output/delSOX9_tcl.SV_RE_intersect.report) was not created"
        exit 1
fi


## Check with the MEF2C gene, without Exomiser, candidateGenesFile and GRCh38
#############################################################################

rm -f "./output/delMEF2C_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/delMEF2C.bed" -outputFile "./output/delMEF2C_tcl.annotated.tsv" -svtBEDcol 4 -REreport 1 -genomeBuild GRCh38 -candidateGenesFile "./input/candidateGenesFile.txt"

# Check that the "MEF2C RE" is annotated
val=`$cut "./output/delMEF2C_tcl.annotated.tsv" "RE_gene"`
if [ `echo $val | egrep -c ";MEF2C \(HI=3/morbid/candidate/RE=EA_enhancer\+GH_Enhancer\)"` != 1 ]
then
        echo "error1: The MEF2C RE is not annotated!"
        exit 1
fi

# Check that the "ACMG_class" is "5"
val=`$cut ./output/delMEF2C_tcl.annotated.tsv "ACMG_class"`
if [ `echo $val | grep -c "5"` != 1 ]
then
        echo "error1: The ACMG_class is not set to 5"
        exit 1
fi

# Check if the RE report exists
if [ ! -e "./output/delMEF2C_tcl.SV_RE_intersect.report" ]
then
        echo "error1: The RE report (./output/delMEF2C_tcl.SV_RE_intersect.report) was not created"
        exit 1
fi


echo "ok - Finished"




