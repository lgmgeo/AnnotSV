#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


# The "KNOTANNOTSV" environment variable should exist
#####################################################
#####################################################
if [ -z ${KNOTANNOTSV+x} ] 
then 
	echo "KNOTANNOTSV is unset"
        exit 1
	
else 
	echo "KNOTANNOTSV is set to '$KNOTANNOTSV'"
fi


# Check 1: file containing different SV types: "test.41_SV.bed"
###############################################################
###############################################################
rm -f "./output/test.41_SV_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.41_SV.bed" -svtBedCol 4 -outputFile "./output/test.41_SV_tcl.annotated.tsv" -hpo "HP:0001156,HP:0001363,HP:0011304,HP:0010055" -genomeBuild GRCh37

for v in `$cut ./output/test.41_SV_tcl.annotated.tsv "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 3 5 ACMG_class full=1 full=3 full=5 full=NA NA" " " "$v"
        then
                echo "Error with ./input/test.41_SV.bed"
                echo "error 1"
                exit 1
        fi
done

## Creation of the tmp yaml file ("./output/config_41_SV_AnnotSV.tmp.yaml")
rm -f "./output/config_41_SV_AnnotSV.tmp.yaml"
$ANNOTSV/tests/data/scripts/complete_KnotAnnotSV_YAML_test.tcl "./output/test.41_SV_tcl.annotated.tsv" "$KNOTANNOTSV/config_AnnotSV.yaml" "./output/config_41_SV_AnnotSV.tmp.yaml"

if [ ! -f "./output/config_41_SV_AnnotSV.tmp.yaml" ] 
then
	echo "File ./output/config_41_SV_AnnotSV.tmp.yaml does not exist"
        echo "error 1"
        exit 1
fi

## knotAnnotSV: creation of "./output/test.41_SV_tcl.annotated.html"
rm -f "./output/test.41_SV_tcl.annotated.html"
perl $KNOTANNOTSV/knotAnnotSV.pl --annotSVfile "./output/test.41_SV_tcl.annotated.tsv" --configFile "./output/config_41_SV_AnnotSV.tmp.yaml" --outDir "./output"

if [ ! -f "./output/test.41_SV_tcl.annotated.html" ]
then
        echo "File ./output/test.41_SV_tcl.annotated.html does not exist"
        echo "error 1"
        exit 1
fi





# Check 2: distributed knot example
###################################
###################################

# ./input/test.7_SV.bed <=> "$KNOTANNOTSV/example/example.bed"
rm -f "./output/test.7_SV_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.7_SV.bed" -svtBedCol 4 -outputFile "./output/test.7_SV_tcl.annotated.tsv" -hpo "HP:0001156,HP:0001363,HP:0011304,HP:0010055" -genomeBuild GRCh37

for v in `$cut "./output/test.7_SV_tcl.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 3 5 ACMG_class full=1 full=3 full=5 full=NA NA" " " "$v"
        then
                echo "Error with ./input/test.7_SV.bed"
                echo "error 1"
                exit 1
        fi
done

## Creation of the tmp yaml file ("./output/config_7_SV_AnnotSV.tmp.yaml")
rm -f "./output/config_7_SV_AnnotSV.tmp.yaml"
$ANNOTSV/tests/data/scripts/complete_KnotAnnotSV_YAML_test.tcl "./output/test.7_SV_tcl.annotated.tsv" "$KNOTANNOTSV/config_AnnotSV.yaml" "./output/config_7_SV_AnnotSV.tmp.yaml"

if [ ! -f "./output/config_7_SV_AnnotSV.tmp.yaml" ]
then
        echo "File ./output/config_7_SV_AnnotSV.tmp.yaml does not exist"
        echo "error 1"
        exit 1
fi

## knotAnnotSV: creation of "./output/test.7_SV_tcl.annotated.html"
rm -f "./output/test.7_SV_tcl.annotated.html"
perl $KNOTANNOTSV/knotAnnotSV.pl --annotSVfile "./output/test.7_SV_tcl.annotated.tsv" --configFile "./output/config_7_SV_AnnotSV.tmp.yaml" --outDir "./output"

if [ ! -f "./output/test.7_SV_tcl.annotated.html" ]
then
        echo "File ./output/test.7_SV_tcl.annotated.html does not exist"
        echo "error 1"
        exit 1
fi



echo "ok - Finished"




