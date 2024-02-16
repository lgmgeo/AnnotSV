#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/data/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}

##################################################################################################
#####################      HG00096 (1000g)      ##################################################
##################################################################################################


# -SVinputFile:
###############
#         HG00096.wgs.mergedSV.v8.20130502.svs.genotypes.bed
#         HG00096.wgs.mergedSV.v8.20130502.svs.genotypes.vcf
#         => 4752 SV


# Check that no column is shifted
#################################
# The "ACMG_class" column must only contain:
# 1
# 3
# 4
# 5
# ACMG_class
# full=1
# full=3
# full=4
# full=5
# full=NA
# NA

SVinputBEDfile="./input/HG00096.wgs.mergedSV.v8.20130502.svs.genotypes.bed"
SVinputVCFfile="./input/HG00096.wgs.mergedSV.v8.20130502.svs.genotypes.vcf"



rm -f ./output/HG00096_tcl_*annotated.tsv


# SVinputFile = BED
$ANNOTSV/bin/AnnotSV -SVinputFile $SVinputBEDfile -svtBEDcol 4 -SVinputInfo 1 -outputFile "./output/HG00096_tcl_1.annotated.tsv" -genomeBuild GRCh37 
for v in `$cut "./output/HG00096_tcl_1.annotated.tsv" "ACMG_class" | sort -u`
do
	if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
	then
    		echo "Error with SVinputFile = BED"
        	echo "error 1"
        	exit 1
	fi
done

# SVinputFile = VCF
$ANNOTSV/bin/AnnotSV -SVinputFile $SVinputVCFfile -SVinputInfo 1 -outputFile "./output/HG00096_tcl_2.annotated.tsv" -genomeBuild GRCh37
for v in `$cut "./output/HG00096_tcl_2.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with SVinputFile = VCF"
                echo "error 1"
                exit 1
        fi
done

# SVinputFile = zipped VCF
gzip $SVinputVCFfile
$ANNOTSV/bin/AnnotSV -SVinputFile $SVinputVCFfile.gz -SVinputInfo 1 -outputFile "./output/HG00096_tcl_3.annotated.tsv" -genomeBuild GRCh37 
for v in `$cut "./output/HG00096_tcl_3.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with SVinputFile = zipped VCF"
                echo "error 1"
                exit 1
        fi
done
gunzip $SVinputVCFfile.gz
 
# -SVinputInfo 0
$ANNOTSV/bin/AnnotSV -SVinputFile $SVinputBEDfile -svtBEDcol 4 -SVinputInfo 0 -outputFile "./output/HG00096_tcl_4.annotated.tsv" -genomeBuild GRCh37
for v in `$cut "./output/HG00096_tcl_4.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with -SVinputInfo 0"
                echo "error 1"
                exit 1
        fi
done

# -metrics fr
$ANNOTSV/bin/AnnotSV -SVinputFile $SVinputBEDfile -svtBEDcol 4 -outputFile "./output/HG00096_tcl_5.annotated.tsv" -metrics fr -genomeBuild GRCh37 
for v in `$cut "./output/HG00096_tcl_5.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with -metrics fr"
                echo "error 1"
                exit 1
        fi
done

 
# -reciprocal 1
$ANNOTSV/bin/AnnotSV -SVinputFile $SVinputBEDfile -svtBEDcol 4 -SVinputInfo 1 -reciprocal 1 -outputFile "./output/HG00096_tcl_6.annotated.tsv" -genomeBuild GRCh37
for v in `$cut "./output/HG00096_tcl_6.annotated.tsv" "ACMG_class" | sort -u`
do
        if ! exists_in_list "1 2 3 4 5 ACMG_class full=1 full=2 full=3 full=4 full=5 full=NA NA" " " "$v"
        then
                echo "Error with -reciprocal 1"
                echo "error 1"
                exit 1
        fi
done



echo "ok - Finished"


