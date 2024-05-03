#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"


#####################################################
## External Gene Annotation File (custom annotations)
#####################################################


# -> Copy "./input/user_annotations.tsv" dans $ANNOTSV/share/doc/AnnotSV/Annotations_Human/Users/ 
#################################################################################################
cp "./input/user_annotations.tsv" $ANNOTSV/share/AnnotSV/Annotations_Human/Users/

# ./input/user_annotations.tsv
##############################
# genes   	del.score       	dup.score       	cnv.score
# ARHGAP19      0.163145773187885       0.339478078191488       0.339117350686952
# PYROXD2 	0.422491880377013       1.23836979707209        1.10456867580225
# HPS1    	1.07868410168786        1.32773222664126        1.53322662103495



# configfile
############
# The external annotation colums will be added only if they are reported in the configfile
# => Use of ./input/configfile

rm -f "./output/input_tcl.annotated.tsv" "./output/output.log"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/input.bed" -outputFile "./output/input_tcl.annotated.tsv" -svtBEDcol 4 -genomeBuild GRCh37 &> "./output/output.log"
#                ...user_annotations.tsv
#                (3 gene identifiers and 3 annotations columns: del.score, dup.score, cnv.score)

if [ `grep -c "annotations columns: del.score, dup.score, cnv.score" "./output/output.log"` == 1 ]
then
        echo "Ok"
else
        echo "error 1: External annotation not used"
        exit 1
fi

#del.score       dup.score       cnv.score
#0.163145773187885       0.339478078191488       0.339117350686952
#0.163145773187885       0.339478078191488       0.339117350686952
if [ `$cut "./output/input_tcl.annotated.tsv" "del.score;dup.score;cnv.score;Gene_name" | grep -c "cnv.score"`  == 1 ]
then
        echo "Ok"
else
        echo "error 1: END coordinates for bracketed SV not good with -includeCI set to 1"
        exit 1
fi


# -> Nettoyage
##############
rm $ANNOTSV/share/AnnotSV/Annotations_Human/Users/user_annotations.tsv



echo "ok - Finished"

