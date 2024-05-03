#!/bin/bash -x

set -eo pipefail


cut="$ANNOTSV/tests/AnnotSV/scripts/cutWithColumnNames.tcl"

function exists_in_list() {
    LIST=$1
    DELIMITER=$2
    VALUE=$3
    [[ "$LIST" =~ ($DELIMITER|^)$VALUE($DELIMITER|$) ]]
}


## Checking the creation of the "Samples_ID" column in the AnnotSV output
## -------------------------------------------------------------------------
## -> automatic creation with a VCF as input
## -> via the "-samplesidBEDcol" option with a BED as input, else exists with "NA" values

## We check:
## - that the "Samples_ID" column exists, and that the values are correct
## - that there is no "shift" with the other columns

## Different options are to be tested (because they are involved in the creation of this feature):
## -SVinputInfo
## -svtBEDcol
## -samplesidBEDcol

## We also need to check a BED with or without header



# SV input BED file with header
###############################
rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtbedcol 4 -samplesidBEDcol 5 -genomeBuild GRCh37
# Samples_ID      ACMG_class
# sample1 	  3
# sample1,sampl2  1
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "1 3 ACMG_class" " " "$v"
        then
                echo "Error with ./input/test.bed (test 1A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID sample1 sample1,sampl2" " " "$v"
        then
                echo "Error with ./input/test.bed (test 1B)"
                echo "error 1"
                exit 1
        fi
done


rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -genomeBuild GRCh37
# => No SV ranking
# => No Samples_ID
# Samples_ID  ACMG_class
# NA          NA
# NA          NA
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 2A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 2B)"
                echo "error 1"
                exit 1
        fi
done




rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtbedcol 4 -genomeBuild GRCh37
# => No Samples_ID
# Samples_ID ACMG_class
# NA         3
# NA         1
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class 3 1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 3A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 3B)"
                echo "error 1"
                exit 1
        fi
done


rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -samplesidBEDcol 5 -genomeBuild GRCh37
# => No SV ranking
# Samples_ID      ACMG_class
# sample1,sampl2  NA
# sample1 NA
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 4A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID sample1,sampl2 sample1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 4B)"
                echo "error 1"
                exit 1
        fi
done



rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 0 -outputFile "./output/test_tcl.annotated.tsv" -svtbedcol 4 -samplesidBEDcol 5 -genomeBuild GRCh37
# Samples_ID      ACMG_class
# sample1 	  3
# sample1,sampl2  1
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class 3 1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 5A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID sample1,sampl2 sample1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 5B)"
                echo "error 1"
                exit 1
        fi
done



rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 0 -outputFile "./output/test_tcl.annotated.tsv" -genomeBuild GRCh37
# => No SV ranking
# => No Samples_ID
# Samples_ID  ACMG_class
# NA          NA
# NA          NA
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 6A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 6B)"
                echo "error 1"
                exit 1
        fi
done



rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 0 -outputFile "./output/test_tcl.annotated.tsv" -svtbedcol 4 -genomeBuild GRCh37
# => No Samples_ID
# Samples_ID    ACMG_class
# NA      	3
# NA      	1
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class 3 1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 7A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 7B)"
                echo "error 1"
                exit 1
        fi
done


rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.bed" -SVinputInfo 0 -outputFile "./output/test_tcl.annotated.tsv" -samplesidBEDcol 5 -genomeBuild GRCh37
# => No SV ranking
$cut "./output/test_tcl.annotated.tsv" "Samples_ID;ACMG_class"
# Samples_ID      ACMG_class
# sample1,sampl2  NA
# sample1         NA
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 8A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID sample1,sampl2 sample1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 8B)"
                echo "error 1"
                exit 1
        fi
done




# SV input BED file without header
##################################

rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test-withoutHeader.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtbedcol 4 -samplesidBEDcol 5 -genomeBuild GRCh37
$cut "./output/test_tcl.annotated.tsv" "Samples_ID;ACMG_class"
# Samples_ID      ACMG_class
# sample1         3
# sample1,sampl2  1
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class 3 1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 9A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID sample1,sampl2 sample1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 9B)"
                echo "error 1"
                exit 1
        fi
done


rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test-withoutHeader.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -genomeBuild GRCh37
# => No SV ranking
# => No Samples_ID
# Samples_ID   ACMG_class
# NA           NA
# NA           NA
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 10A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 10B)"
                echo "error 1"
                exit 1
        fi
done




rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test-withoutHeader.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -svtbedcol 4 -genomeBuild GRCh37
# => No Samples_ID
# Samples_ID      ACMG_class
# NA      3
# NA      1
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class 3 1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 11A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 11B)"
                echo "error 1"
                exit 1
        fi
done


rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test-withoutHeader.bed" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -samplesidBEDcol 5 -genomeBuild GRCh37
# => No SV ranking
# Samples_ID      ACMG_class
# sample1,sampl2  NA
# sample1         NA
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 12A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID sample1,sampl2 sample1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 12B)"
                echo "error 1"
                exit 1
        fi
done


rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test-withoutHeader.bed" -SVinputInfo 0 -outputFile "./output/test_tcl.annotated.tsv" -svtbedcol 4 -samplesidBEDcol 5 -genomeBuild GRCh37
# Samples_ID      ACMG_class
# sample1 	  3
# sample1,sampl2  1
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class 3 1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 13A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID sample1,sampl2 sample1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 13B)"
                echo "error 1"
                exit 1
        fi
done


rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test-withoutHeader.bed" -SVinputInfo 0 -outputFile "./output/test_tcl.annotated.tsv" -genomeBuild GRCh37
# => No SV ranking
# => No Samples_ID
# Samples_ID      ACMG_class
# NA      NA
# NA      NA
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 14A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 14B)"
                echo "error 1"
                exit 1
        fi
done



rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test-withoutHeader.bed" -SVinputInfo 0 -outputFile "./output/test_tcl.annotated.tsv" -svtbedcol 4 -genomeBuild GRCh37
# => No Samples_ID
# Samples_ID      ACMG_class
# NA              3
# NA              1
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class 3 1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 15A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 15B)"
                echo "error 1"
                exit 1
        fi
done


rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test-withoutHeader.bed" -SVinputInfo 0 -outputFile "./output/test_tcl.annotated.tsv" -samplesidBEDcol 5 -genomeBuild GRCh37
# => No SV ranking
# Samples_ID      ACMG_class
# sample1,sampl2  NA
# sample1         NA
for v in `$cut "./output/test_tcl.annotated.tsv" "ACMG_class"`
do
        if ! exists_in_list "ACMG_class NA" " " "$v"
        then
                echo "Error with ./input/test.bed (test 16A)"
                echo "error 1"
                exit 1
        fi
done
for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID"`
do
        if ! exists_in_list "Samples_ID sample1,sampl2 sample1" " " "$v"
        then
                echo "Error with ./input/test.bed (test 16B)"
                echo "error 1"
                exit 1
        fi
done



# SV input VCF file
###################

rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.vcf" -SVinputInfo 1 -outputFile "./output/test_tcl.annotated.tsv" -genomeBuild GRCh37
# $cut "./output/test_tcl.annotated.tsv" "Samples_ID;ACMG_class"| grep -v full | sort -u
# Samples_ID      ACMG_class
# sample1 3
# sample1,sample2 NA
# sample2 1
# sample2 3
# sample2 NA

annotations=`$cut "./output/test_tcl.annotated.tsv" "ACMG_class" | grep -v full | sort -u` 
if [ `echo $annotations | grep -c "1 3 ACMG_class NA"` != 1 ]
then
        echo "Error with ./input/test.bed (test 17A)"
        echo "error 1"
        exit 1
fi

for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID" | grep -v full | sort -u`
do
        if ! exists_in_list "Samples_ID sample1 sample1,sample2 sample2" " " "$v"
        then
                echo "Error with ./input/test.bed (test 17B)"
                echo "error 1"
                exit 1
        fi
done

rm -f "./output/test_tcl.annotated.tsv"
$ANNOTSV/bin/AnnotSV -SVinputFile "./input/test.vcf" -SVinputInfo 0 -outputFile "./output/test_tcl.annotated.tsv" -genomeBuild GRCh37
# $cut "./output/test_tcl.annotated.tsv" "Samples_ID;ACMG_class"| grep -v full | sort -u
# Samples_ID      ACMG_class
# sample1 3
# sample1,sample2 NA
# sample2 1
# sample2 3
# sample2 NA
annotations=`$cut "./output/test_tcl.annotated.tsv" "ACMG_class" | grep -v full | sort -u`
if [ `echo $annotations | grep -c "1 3 ACMG_class NA"` != 1 ]
then
        echo "Error with ./input/test.bed (test 18A)"
        echo "error 1"
        exit 1
fi

for v in `$cut "./output/test_tcl.annotated.tsv" "Samples_ID" | grep -v full | sort -u`
do
        if ! exists_in_list "Samples_ID sample1 sample1,sample2 sample2" " " "$v"
        then
                echo "Error with ./input/test.bed (test 18B)"
                echo "error 1"
                exit 1
        fi
done




echo "ok - Finished"



